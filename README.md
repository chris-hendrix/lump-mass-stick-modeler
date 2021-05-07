# Lumped Mass Stick Modeler
Simple python object that allows the modeling of MDOF lumped mass stick models (LMSM) subjected to dynamic loading.

## How to Use
See the following jupyter notebooks:
 - `toy_example.ipynb`: example of a simple SDOF toy and a tuned-mass damper
 - `building_example.ipynb`: to be completed
## Theory
### Equation of Motion SDOF
The movement of a single degree of system (Figure 2) subjected to a time varying force, $p(t)$, with a stiffness $k$, and a damping coefficient, $c$, can be  modeled with the following equation of forces:

<p align="center"><img src="https://latex.codecogs.com/gif.latex?
m\ddot{u}+c\dot{u}+ku=p(t)
" /></p>
<p align="center"><img src="images\sdof.jpg" width=50%/></p>

### Equation of Motion LMSM
A depiction of a system with n degrees of freedom is shown below. This model is displaced in lumped mass model form.
<p align="center"><img src="https://latex.codecogs.com/gif.latex?
\tilde{m}\ddot{q}+\tilde{c}\dot{q}+\tilde{k}q=\tilde{p}(t)
" /></p>


<p align="center"><img src="images\lmsm.jpg" width=25%></p>

The mode shape vector of the MDOF system can be determined by using the following eigenvalue equation, where \lambda is a dimensionless eigenvalue:
<p align="center"><img src="https://latex.codecogs.com/gif.latex?
\left(\mathbf{K}-\lambda\mathbf{M}\right)\mathbf{\Phi}=\mathbf{0}
" /></p>


Each fundamental mode shape, comprises the mode shape vector:
<p align="center"><img src="https://latex.codecogs.com/gif.latex?
\mathbf{\Phi}=[{\phi}_{1},{\phi}_{2},\ldots{\phi}_{j}\ ]
" /></p>

Modal values can be determined for the system. Below are the modal mass, modal stiffness, and modal damping values for mode j.
<p align="center"><img src="https://latex.codecogs.com/gif.latex?
{\widetilde{m}}_{j}={\phi}_{j}^T\mathbf{M}{\phi}_{j},~
{\widetilde{k}}_{j}={\phi}_{j}^T\mathbf{K}{\phi}_{j},~	
{\widetilde{c}}_{j}={\phi}_{j}^T\mathbf{C}{\phi}_{j}, ~
{\widetilde{p}}_{j}={\phi}_{j}^T\mathbf{P}
" /></p>

The mass matrix, stiffness matrix and the damping matrix can be developed as follows
<p align="center"><img src="images\matm.jpg" width=25%></p>
<p align="center"><img src="images\matk.jpg" width=100%></p>
<p align="center"><img src="images\matc.jpg" width=100%></p>
<p align="center"><img src="images\matpu.png" width=40%></p>


The displacement parameter that defines the response of the system can be calculated as follows if there are no initial values
<p align="center"><img src="https://latex.codecogs.com/gif.latex?
q\left(t\right)=\widetilde{h}\ast\widetilde{p}
" /></p>

The modal impulse response function is
<p align="center"><img src="https://latex.codecogs.com/gif.latex?
{\widetilde{h}}_j\left(t\right)=\frac{1}{{\widetilde{m}}_j\omega_d}e^{-{\widetilde{\xi}\omega}_{n,j}t}
" /></p>

Finally, the response vector can be calculated for mode j:
<p align="center"><img src="https://latex.codecogs.com/gif.latex?
\mathbf{u}\left(t\right)=q\left(t\right)\mathbf{\phi}_{j}
" /></p>