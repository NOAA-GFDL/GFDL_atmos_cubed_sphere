The horizontal discretization and Lagrangian Dynamics {#horizontal}
=========================================

##Chapter 5 

GFDL will provide the additional documentation by the end of May 2021.

(This information is a reproduction of Section 5 from HCZC20)

###5.1 Dynamics along Lagrangian surfaces

The layer-integrated equations (4.3) are discretized along the Lagrangian surfaces and integrated on the “acoustic” or “dynamical” time step \f$\delta t\f$ using forward-backward time-stepping as in LR97. The vertical velocity \f$w\f$ is a three-dimensional cell-mean value and partially advanced using the advection scheme. The geometric layer depth \f$\delta z\f$ is simply the difference of the heights of the successive layer interfaces, which with \f$\delta p^*\f$ defines the layer-mean density and the location of the Lagrangian surfaces. The air mass is the total air mass, including water vapor and condensate species.

FV<sup>3</sup> places the wind components using the D-grid (following Arakawa’s terminology), which defines the winds as face-tangential quantities. The D-grid permits us to compute the cell-mean absolute vorticity \f$\Omega\f$ exactly using Stokes’ theorem and a cell-mean value of the local Coriolis parameter \f$f\f$, without performing any averages or interpolations. The wind components themselves are face-mean values “along the cell edges” (not cell-mean values).

Following the notation from L04, PL07, and HL13, we can write the discretized forms of (4.3), excluding the vertical components, as:

\f[
 \delta p^{*(n + 1)} = \delta p^{*n} + F[\tilde{u^*}, \delta p_y] + G [\tilde{v^*}, \delta p_x]  \\ \tag {5.1}
  \f]

\f[
 \Theta^{n + 1} = \frac{1}{\delta p^{*(n+1)}} [\Theta^n \delta p^{*n} + F[X^*, \Theta_y] + G[Y^*, \Theta_x]]  \\ \tag {5.2}
  \f]

\f[
 w^* = \frac{1}{\delta p^{*(n+!)}} [w^n \delta p^{*n} + F[X^*, w_y] + G[Y^*, w_x]]  \\ \tag {5.3}
  \f]

\f[
 u^{n + 1} = u^n + \triangle \tau [Y(\tilde{v^*}, \Omega_x) - \delta_x (K^* -v \nabla^2 D) + \hat{P_x}]  \\ \tag {5.4}
  \f]

\f[
 v^{n + 1} = v^n + \triangle \tau [X(\tilde{u^*}, \Omega_y) - \delta_y (K^* -v \nabla^2 D) + \hat{P_y}]  \\ \tag {5.5}
  \f]

The quantities \f$\hat{P_x}\f$ \f$\hat{P_y}\f$ are the horizontal pressure-gradient force terms computed as in L97, the primary difference being that the forces due to hydrostatic and nonhydrostatic pressures arecomputed separately, and that the hydrostatic pressure gradient computation uses the log of the pressure to improve accuracy. The vertical nonhydrostatic pressure-gradient force is evaluated by the semi-implicit solver described in Chapter 5; only the forward advection of \f$w\f$ is performed during the Lagrangian dynamics, producing a partially-updated \f$w^*\f$.

For stability, the pressure gradient force is evaluated backwards-in-time: the flux terms in the momentum, mass, and entropy equations are evaluated forward by the advection scheme, and the resulting updated fields are used to compute the pressure gradient force. This forward-backward time-stepping is stable without needing to use predictor-corrector or Runge-Kutta methods.

In nonhydrostatic simulations it is recommended that the time off-centering for the horizontal pressure-gradient force be consistent with that used in the semi-implicit solver, which includes the vertical nonhydrostatic pressure-gradient force computation,to ensure consistency between the two. If the semi-implicit solver is run fully-implicit \f$(\alpha=1)\f$ then the pressure-gradient force should be evaluated fully backward \f$(\beta = 0)\f$; otherwise use \f$\beta = 1 - \alpha\f$
