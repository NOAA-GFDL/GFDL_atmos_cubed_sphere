The vertical Lagrangian Solver {#lagrangian}
=========================================

##Chapter 4 

GFDL will provide the additional documentation by the end of May 2021.

(This information is a reproduction of Sections 3 and 4 from HCZC20)

###4.1 Lagrangian vertical coordinates

A *Lagrangian* vertical coordinate is used in FV<sup>3</sup>. This coordinate uses the depth of each layer (in terms of mass or as geometric height) as a prognostic variable, allowing the layer interfaces to deform freely as the flow evolves. Further, the flow is constrained within the Lagrangian layers, with no flow across the layer interfaces (even for non-adiabatic flows). Instead, the flow deforms the layers themselves by advecting the layer thickness and by straining the layers by the vertical gradient of explicit vertical motion. This form is automatically consistent with the LR96 scheme, avoids the need for explicit calculation and dimensional splitting of verticaladvection, greatly reduces implicit vertical diffusion, and has no vertical Courant number restriction.

FV<sup>3</sup> uses a hybrid-pressure coordinate based on the hydrostatic surface pressure \f$p_s^*\f$:

\f[
 p_k^* = a_k + b_kp_s^*  \\  \tag {4.1}
 \f]

where *k* is the vertical index of the layer interface, counting from the top down, and \f$a_k\f$, \f$b_k\f$ are pre-defined coefficients. Typically, the top interface is at a constant pressure \f$p_T\f$, so \f$a_0 = p_T\f$ and \f$b_0 = 0\f$. The spacing of the levels depends on the particular application, and is chosen depending on how high of a model top is desired, where additional vertical resolution is applied (typically in the boundary layer, but sometimes also near the tropical tropopause), and where to transition from hybrid \f$b_k > 0\f$ to pure pressure \f$b_k = 0\f$ coordinates.

###4.2 Prognostic variables and governing equations

The mass per unit area \f$\delta m\f$ can be expressed in terms of the difference in hydrostatic pressure \f$\delta p^*\f$ between the top and bottom of the layers; and, using the hydrostatic equation, can also be written in terms of the layer depth \f$\delta z^1\f$:

\f[
 \delta m = \frac{\delta p^*}{g} = -p \delta z  \\  \tag {4.2} 
  \f]

The continuous Lagrangian equations of motion, in a layer of finite depth \f$\delta z\f$ and mass \f$\delta p^*\f$, are then given as

\f[
 D_L \delta p^* + \nabla \cdot (V \delta p^*) = 0  \\
 D_L \delta p^* \Theta_v + \nabla \cdot (V \delta p^* \Theta_v) = 0  \\
 D_L \delta p^*w + \nabla \cdot (V \delta p^*w) = -g \delta z \frac{\partial p'}{\partial z}  \\
 D_L u = \Omega u - \frac{\partial}{\partial x} K - \frac{1}{p} \frac{\partial p}{\partial x} \biggr\rvert_{z} \\
 D_L v = - \Omega u - \frac{\partial}{\partial y} K - \frac{1}{p} \frac{\partial p}{\partial y} \biggr\rvert_{z}  \\  \tag {4.3}
  \f]

Note that these equations are exact: no discretization has been made yet, and the only change from the original differential form of Euler’s equations is to integrate over an arbitrary depth \f$\delta p^*\f$. The operator \f$D_L\f$ is the “vertically-Lagrangian” derivative, formally equal to \f$\frac{\partial \psi}{\partial t} + \frac{\partial}{\partial z}(w \psi)\f$ for an arbitrary scalar \f$\psi\f$. The flow is entirely along the Lagrangian surfaces, including the vertical motion (which deforms the surfaces as appropriate, an effect included in the semi-implicit solver).

####Prognostic variables in FV<sup>3</sup>

Variable         | Description
:--------------: | :---------- 
\f$\delta p^*\f$ | Vertical difference in hydrostatic pressure, proportional to mass
\f$u\f$          | D-grid face-mean horizontal x-direction wind
\f$v\f$          | D-grid face-mean horizontal y-direction wind
\f$\Theta_v\f$   | Cell-mean virtual potential temperature
\f$w\f$          | Cell-mean vertical velocity
\f$\delta z\f$   | Geometric layer height

The vertical component of absolute vorticity is given as \f$\Omega\f$ and \f$p\f$ is the full nonhydrostatic pressure. The kinetic energy is given as \f$K = \frac{1}{2}(\tilde{u}u + \tilde{v}v)\f$: since FV<sup>3</sup> does not assume that the horizontal coordinate system is orthogonal, we use the covariant (\f$u\f$ and \f$v\f$) components of the wind vector as prognostic variables and the contravariant (\f$\tilde{u}\f$ and \f$\tilde{v}\f$) components for advection, avoiding the need to explicitly include metric terms. See PL07 and HL13 for more information about covariant and contravariant components.

The nonhydrostatic pressure gradient term in the \f$w\f$ equation is computed by the semi-implicit solver described section 5, which also computes the prognostic equation for \f$\delta z\f$. There is no projection of the vertical pressure gradient force into the horizontal; similarly, there is no projection of the horizontal winds \f$u\f$, \f$v\f$into the vertical, despite the slopes of the Lagrangian surfaces.

Finally, the ideal gas law:

\f[
 p = p^* + p' = \rho R_dT_v = - \frac{\partial p^*}{g \delta z} R_d T_v  \\  \tag {4.4}
  \f]

where

\f[
 T_v = T(1 + \in q_v)(1 - q_{cond})  \\ \tag {4.5}
  \f]

is the “condensate modified” virtual temperature, or density temperature, is used to close the system of equations. Here, \f$d_{cond}\f$ is the (moist) mixingratio of all of the liquid and solid-phase microphysical species, if present. When the gas law is used, the mass \f$p^*\f$ in this computation must be the mass of only the dry air and water vapor, and not including the mass of the condensate (non-gas) species. A rigorous derivation of the virtual and density temperatures is given in K. Emanuel, Atmospheric Convection (1994, Oxford), Sec. 4.3.

These equations are also applicable to hydrostatic flow, in which \f$w\f$ is not prognosed and \f$p = p^*\f$  is entirely hydrostatic.


______________________________
<sup>1</sup> In this document, to avoid confusion we write \f$\delta z\f$ as if it is a positive-definite quantity. In the solver itself, \f$\delta z\f$ is defined to be negative-definite, incorporating the negative sign into the definition of \f$\delta z\f$; this definition has the additional advantage of being consistent with how \f$\delta z\f$ is defined, being measured as the difference in hydrostatic pressure between the bottom and top of a layer.
