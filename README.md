# Lumped Mass Stick Modeler
Simple python tool to view the motions of MDOF lumped mass stick models (LMSM).

## Theory
### Equation of Motion SDOF
The movement of a single degree of system (Figure 2) subjected to a time varying force, $p(t)$, with a stiffness $k$, and a damping coefficient, $c$, can be  modeled with the following equation of forces:


$$
\begin{align}
m\ddot{u}+c\dot{u}+ku=p(t)
\end{align}
$$
<p align="center"><img src="images\sdof.jpg" width=50%/></p>

### Equation of Motion LMSM
A depiction of a system with n degrees of freedom is shown below. This model is displaced in lumped mass model form.

$$
\begin{align}
\tilde{m}\ddot{q}+\tilde{c}\dot{q}+\tilde{k}q=\tilde{p}(t)
\end{align}
$$

<p align="center"><img src="images\lmsm.jpg" width=25%></p>

The mode shape vector, ${\Phi}$, of the MDOF system can be determined by using the following eigenvalue equation, where \lambda is a dimensionless eigenvalue:

$$
\begin{align}
\left({K}-\lambda{M}\right){\Phi}={0}
\end{align}
$$

Each fundamental mode shape, ${\phi}_{j}=\left[\phi_{j1};\phi_{j2};\ldots;\phi_{jn}\right]$, comprises the mode shape vector:
$$
\begin{align}
{\Phi}=[{\phi}_{1},{\phi}_{2},\ldots{\phi}_{j}\ ]
\end{align}
$$

Modal values can be determined for the system. Below are the modal mass, modal stiffness, and modal damping values for mode j.
$$
\begin{align}
{\widetilde{m}}_{j}={\phi}_{j}^T{M}{\phi}_{j}\ ,	{\widetilde{k}}_{j}={\phi}_{j}^T{K}{\phi}_{j},	{\widetilde{c}}_{j}={\phi}_{j}^T{C}{\phi}_{j}, {\widetilde{p}}_{j}={\phi}_{j}^T{P}
\end{align}
$$

The mass matrix, stiffness matrix and the damping matrix can be developed as follows
<p align="center"><img src="images\matm.jpg" width=25%></p>
<p align="center"><img src="images\matk.jpg" width=100%></p>
<p align="center"><img src="images\matc.jpg" width=100%></p>
<p align="center"><img src="images\matpu.png" width=40%></p>


The displacement parameter that defines the response of the system can be calculated as follows if there are no initial values
$$
\begin{align}
q\left(t\right)=\widetilde{h}\ast\widetilde{p}
\end{align}
$$

The modal impulse response function is
$$
\begin{align}
{\widetilde{h}}_j\left(t\right)=\frac{1}{{\widetilde{m}}_j\omega_d}e^{-{\widetilde{\xi}\omega}_{n,j}t}
\end{align}
$$
Finally, the response vector can be calculated for mode j:
$$
\begin{align}
\mathbf{u}\left(t\right)=q\left(t\right)\mathbf{\phi}_\mathbf{j}
\end{align}
$$