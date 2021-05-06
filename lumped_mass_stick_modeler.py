import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from scipy import linalg
import math


class LumpedMassStickModel:
    def __init__(self, masses, stiffnesses, damp_coeffs, loadings=None, timestep=None, node_locations=None, name=''):
        self.name = name
        self.masses = masses
        self.stiffnesses = stiffnesses
        self.damp_coeffs = damp_coeffs
        self.nodes = len(masses)
        self.debug = []
        self.calc_modal_propoerties()
        self.node_locations = node_locations

        # default mass locations spaced at 1
        if node_locations == None: self.node_locations = np.arange(1, self.nodes+1, 1)

        # set_loadings if not equal to None
        if not loadings == None and not timestep==None: self.set_loadings(loadings, timestep)

        # TODO add error checking for list lengths

    def calc_modal_propoerties(self):
        
        # append zeroes to handle i+1
        m_vals = self.masses.copy()
        k_vals = self.stiffnesses.copy() + [0]
        c_vals = self.damp_coeffs.copy() + [0]

        # setup arrays
        length = len(m_vals)
        m = np.zeros((length,length))
        k = np.zeros((length,length))
        c = np.zeros((length,length))

        for i in range(0, length):
            # diagonal
            m[i,i] = m_vals[i]
            k[i,i] = k_vals[i] + k_vals[i+1]
            c[i,i] = c_vals[i] + c_vals[i+1]

            if i < length -1:
                # right of diagonal 
                k[i,i+1] = -1*k_vals[i+1]
                c[i,i+1] = -1*c_vals[i+1]

                # left of diagonal
                k[i+1,i] = -1*k_vals[i+1]
                c[i+1,i] = -1*c_vals[i+1]
        
        # get and sort eigenvalues, eigenvectors
        eig_vals, eig_vecs = linalg.eig(k,m)
        idx = eig_vals.argsort()[::1]
        phi = eig_vecs[:,idx]
        wn = np.sqrt(np.real(eig_vals[idx]))

        # normalize
        for i in range(0, length):
            phi[:,i] = phi[:,i]/phi[length-1,i]

        # calculate modal values
        mm = np.dot(phi.T,np.dot(m,phi)).flatten()[0]
        km = np.dot(phi.T,np.dot(k,phi)).flatten()[0]
        cm = np.dot(phi.T,np.dot(c,phi)).flatten()[0]
        xim = cm/(2*wn*mm)
        gam = phi[0,0]/mm*m*phi.sum()
        wnd = wn*np.sqrt(1-xim)
        
        # set matrices
        self.mass_matrix = m
        self.stiffness_matrix = k
        self.damp_coeff_matrix = c
        self.mode_shapes = phi
        self.ang_freqs = wn
        self.freqs = wn/(2*math.pi)

        # set modal values
        self.modal_mass = mm
        self.modal_stiffness = km
        self.modal_damp_coeff = cm
        self.modal_damp_ratio = xim
        self.gamma = gam
        self.ang_freqs_damped = wnd
        self.freqs_damped = wnd/(2*math.pi)
        self.freqs = wnd
        self.modes = len(wn)
    
    def set_loadings(self, loadings, timestep):
        p = loadings.copy()
        phi = self.mode_shapes
        pm = np.dot(phi.T,p)

        # calculate displacements
        mode = 0
        mm = self.modal_mass
        xim = self.modal_damp_ratio
        wnd = self.ang_freqs_damped
        nodes = self.nodes
        modes = self.modes
        dt = timestep
        steps = len(p[0])
        t = np.arange(steps) * dt
        
        # setup disp, vel, and acc matrices
        u_mat = np.zeros(shape=(modes, nodes, steps))
        v_mat = np.zeros(shape=(modes, nodes, steps))
        a_mat = np.zeros(shape=(modes, nodes, steps))

        for node in range(0, nodes):
            for mode in range(0, modes):
                h = 1/(mm*wnd[mode])*np.exp(-xim[mode]*wnd[mode]*t)*np.sin(wnd[mode]*t)
                q = (np.convolve(pm[mode].flatten(), h)*dt)[0:steps]
                self.debug.append({'q':q})
                u = phi[mode][node]*q
                v = np.gradient(u, dt)
                a = np.gradient(v, dt)

                u_mat[mode,node,:] = u
                v_mat[mode,node,:] = v
                a_mat[mode,node,:] = a
        
        # set loading properties
        self.loadings = p
        self.times = t
        self.timestep = timestep
        self.modal_loadings = pm

        # set disp, vel, acc matrix properties
        self.disps = u_mat
        self.vels = v_mat
        self.accs = a_mat

    def get_animation(self,fig, ax, mode=0):

        # get data and add leading zero
        dt = self.timestep
        xdata = [0] + np.arange(self.nodes).tolist()
        ydata = [0] + self.node_locations.tolist()

        # get maxab bounds
        xmin = self.disps[mode,:,:].min()
        xmax = self.disps[mode,:,:].max()
        xmaxabs = max(abs(xmin), abs(xmax))

        # set ax properties
        fnd = self.freqs_damped[mode]
        title = f'{self.name}\n(mode={mode}, fnd={round(fnd,3)}Hz)'
        ax.set_title(title)
        ax.set_xlabel('disp')
        ax.set_ylabel('position')
        ax.set_xlim(-1*xmaxabs, xmaxabs)

        # format line
        line, = ax.plot(xdata, ydata)
        line.set_marker('o')
        line.set_markersize(10)

        # animate line
        def init():
            line.set_xdata(xdata)
            return line,

        def animate(i):
            xdata = [0] + self.disps[mode,:,i].tolist()
            line.set_xdata(xdata)  # update the data.
            #line.set_marker('o')
            #line.set_markersize(10)
            return line,

        ani = animation.FuncAnimation(
            fig, animate, init_func=init, interval=1000*dt, blit=True, save_count=50)
        
        return ani

    def get_time_history(self, mode=0, node=0, motion='a'):
        motion = motion[0].lower()
        motions = self.accs
        if motion=='d': motions = self.disps
        if motion=='v': motions = self.vels
        if motion=='l': motions = self.loadings

        return self.times, motions[mode][node]

    def get_time_history_plt(self, ax, mode=0, node=0, motion='a'):
        title = f'Time History (mode= {mode}, node={node})'
        x, y = self.get_time_history(node, mode, motion)
        ax.set_title(title)
        ax.set_xlabel('time')
        ax.set_ylabel(motion)
        ax.plot(x, y, label=f'{self.name}')
        return ax





