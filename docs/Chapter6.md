The Nonhydrostatic Solver {#solver}
=========================================

##Chapter 6

GFDL will provide the additional documentation by the end of May 2021.

(This information is a reproduction of Section 6 from HCZC20)

An equation for \f$z\f$ can be derived from the definition of \f$w\f$:

\f[
 w = \frac{Dz}{Dt} = D_Lz + \vec{U} \cdot \nabla z  \\ \tag {6.1}
  \f]

The time-tendency of geopotential height is then the sum of the advective height flux along the Lagrangian interfaces and the vertical distortion of the surfaces by the gradient of \f$z\f$. Discretizing:

\f[
 z^{n + 1} = z^n + F [\tilde{u^*}, \delta z_y] + G[\tilde{v^*}, \delta z_x] + w^{n + 1} \triangle t  \\ \tag {6.2} 
  \f]

Since \f$w\f$ is solved for on the interfaces, we can then simply take the vertical difference to get \f$\delta z\f$.

Recalling that the Lagrangian dynamics in (7) only performs the forward advection of the vertical velocity, yielding \f$w^*\f$, we then need to evaluate the vertical pressure-gradient force:

\f[
 w^{n + 1} = w^* -g \delta z^{n + 1} \delta_z p'^{n + 1}  \\ \tag {6.3}
  \f]

The pressure perturbation \f$p'\f$ can be evaluated from the ideal gas law,

\f[
 p' = p  - p^* = \frac{\delta p}{g \delta z} R_d T_v - p^*  \\ \tag {6.4}
  \f]

requiring simultaneous solution of \f$w\f$, \f$p'\f$, and \f$\delta z\f$ using a tridiagonal solver.

There is an option to off-center the semi-implicit solver toreduce implicit diffusion. The parameter \f$\alpha\f$ can be varied between 0.5 and 1 to control the amount of off-centering, with \f$\alpha = 1\f$ being fully-implicit.  As discussed in Chapter 4 this off-centering parameter should be set to \f$\alpha = \beta - 1\f$, consistent with that used for the horizontal pressure-gradient force.

The boundary conditions used are \f$p' = 0\f$ at the model top, and \f$\vec{U} \widehat{\cdot} n_s = 0\f$ at the lower boundary of \f$z = z_s\f$. This is the “free-slip” boundary condition, that the lower boundary is a streamline.  The surface vertical velocity \f$w_s\f$ can be computed from (11) by advecting the surface height \f$z_s\f$:

\f[
 w_s = \frac{z^*_s - z_s}{\triangle t}  \\ \tag {6.5}
  \f]

where \f$z_s^*\f$ is the advected value and \f$z_s\f$ the height of the topography.
