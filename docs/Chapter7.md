Stabilization and filtering options {#stabilization}
=========================================

##Chapter 7 

(This information is a reproduction of Chapter 2 from LPH17)

###7.1 Divergence damping
Horizontal divergence (along a Lagrangian surface) is computed as a cell integrated quantity on the dual grid:

\f[
 D =   \frac{1}{\Delta A_c} \left[ \delta_x (u_c \Delta y_c sin\alpha) + \delta_y (v_c \Delta x_c sin \alpha) \right]    \\  \tag {7.1}
  \f]

The Laplacian of D can also be computed as a cell-integrated quantity on the dual grid:

\f[
 \nabla^2 D =   \frac{1}{\Delta A_c} \left[ \delta_x \left(\frac{\delta_x D}{\Delta x} \Delta y_c sin\alpha \right) + \delta_y \left(\frac{\delta_y D}{\Delta y} \Delta x_c sin \alpha \right) \right]    \\  \tag {7.2}
  \f]

This operator can be applied on &nabla;<sup>2</sup>D instead of D to yield &nabla;<sup>4</sup>D. The damping is then applied when the forward time step is taken for the horizontal dynamics along vertically-Lagrangian surfaces:

\f[
 u^{n+1} =  u^n + . . . +    \nu_D \frac{\delta_x  \nabla^{2N} D}{\Delta x}   \\  \tag {7.3}
  \f]

\f[
 v^{n+1} =  v^n + . . . +    \nu_D \frac{\delta_y  \nabla^{2N} D}{\Delta y}   \\  \tag {7.4}
  \f]

where N (equal to the namelist parameter `nord`) is 1 for fourth-order and 2 for sixth-order damping. The nondimensional damping coefficient is given

\f[
   \nu_D =  (d_4\Delta A_{min} )^{N+1}  \\  \tag {7.5}
  \f]

in which d<sub>4</sub> is the parameter `d4_bg` in the namelist, and &Delta;A<sub>min</sub> is the *global* minimum grid-cell area. It is recommended that this parameter be set to a value between 0.1 and 0.16, with instability likely for higher or lower values. Note that divergence damping is necessary as there is no implicit damping on the divergence in FV<sup>3</sup>. An optional second-order &nabla;<sup>2</sup> damping, in addition the higher-order divergence damping, can be applied as well; in this case the added damping is of the form &nu;<sub>2D</sub>(&delta;<sub>x</sub>D/ &Delta;x), where &nu;<sub>2D</sub> = d<sub>2</sub>&Delta;A<sub>min</sub>. Typically, the coefficient for d<sub>2</sub> should be much smaller—by at least an order of magnitude—than the higher-order coefficient, if it is used at all, since the second-order damping is only weakly scale-selective and will significantly diffuse even resolved-scale features.

The divergence damping can also be modified to add an approximate Smagorinsky-type damping, in which second-order divergence damping can be added to the flow dependent on the amount of stretching and dilation in thepflow. In this case, the d<sub>2</sub> in the expression for &nu;<sub>2D</sub> is replaced by d<sub>S</sub>&Delta;t(D<sup>2</sup> + &zeta;<sup>2</sup>)<sup>1/2</sup>, where d<sub>S</sub> is the Smagorinsky coefficient (typically set to 0.2 if used) and &zeta; is the relative vorticity interpolated to cell corners so as to be co-located with the divergence. This form of the damping coefficient is more physical than the artificial damping discussed in the rest of this chapter, and will typically be very weak except in regions of very strong flow deformation.

Divergence and flux damping (described in the next section) are both applied entirely within Lagrangian surfaces; there is no explicit diffusion applied across the surfaces. However, in regions of vertical motion the Lagrangian surfaces are distorted in the vertical as they follow the flow upward or downward. The amount of the distortion depends on the along-surface gradient of the vertical velocity; so where the distortion is largest is where there is the strongest horizontal shearing of the vertical velocity, which is also where &nabla;<sup>2n</sup> of the scalar fields should be the largest.


###7.2 Hyperdiffusion (flux, or “vorticity”) damping

Traditionally in FV<sup>3</sup> computational noise is controlled through two means: the explicit divergence damping, and the implicit, nonlinear diffusion from the monotonicity constraint used in computing the fluxes. However, the implicit diffusion may be too heavily damping of marginally-resolved flow features, and for high-resolution applications it may be beneficial to use a nonmonotonic scheme and then add a user-controllable hyperdiffusion instead. This added hyperdiffusion need not be as strong as the divergence damping, since the non-monotonic advection schemes are still weakly diffusive (while no implicit diffusion is applied to the divergence), and often the hyperdiffusion coefficient d<sub>f</sub> (`vtdm4`) should be much less than the divergence damping coefficient d<sub>4</sub>.

In FV<sup>3</sup> the hyperdiffusion is primarily meant to diffuse the kinetic energy of the rotational component of the flow, similarly to how the divergence damping dissipates the kinetic energy in the divergence component of the flow. (For this reason, the hyperdiffusion is sometimes called “vorticity” damping.) The diffusion is applied to the vorticity *flux*, allowing application of diffusion to only the rotational component fo the flow without needing to explicitly compute the rotational component.To maintain consistent advection of the other prognostic variables—w, &theta;<sub>v</sub>, and &delta;p<sup>*</sup>—the fluxes for these quantities are diffused as well, so that the potential vorticity and updraft helicity are still consistently advected as if they were scalars. (There is no way to add diffusion to the tracer mass fluxes, since higher-order diffusion cannot ensure monotonicity preservation without performing an additional flux limiting, adding even more implicit diffusion.)

The hyperdiffusion is equal to the order of the divergence damping, unless eighth-order divergence damping is used, for which the hyperdiffusion remains sixth-order. The diffusion operator itself is second-order, but higher order diffusion is easily computed by repeatedly applying the diffusion operator, as is done for the divergence damping.

Vertical vorticity is a cell-integrated quantity on the model grid. The vorticity fluxes v&Omega;/&Delta;x and -u&Omega;/&Delta;y are used to update the vector-invariant momentum equations. We can apply damping on the vorticity as well; to maintain consistent advection, the same damping is applied to the mass, heat, and vertical momentum fields, all of which are co-located with the vorticity. This additional damping is beneficial when using a non-monotonic advection scheme, which lacks the implicit diffusion of monotonic advection.

Since the diffusion added to the vorticity fluxe is known explicitly, the loss of kinetic energy due to this explicit diffusion can be computed. The lost energy optionally can be added back as heat, after applying a horizontal smoother to the energy density field (so as not to restore the grid-scale noise which the damping was intended to remove). This can greatly improve the dynamical activity on marginally-resolved scales.

###7.3 Energy-, momentum-, and mass-conserving 2&Delta;z filter

Local Richardson-number dynamic instabilities can create instabilities, especially in higher-top models, if the vertical turbulence scheme in the physical parameterizations is either disabled or insufficiently-strong to relieve these instabilities. These instabilities can grow large enough to crash the model. To relieve these instabilities, FV<sup>3</sup> has the option to use a local (2&Delta;z), vertical mixing to release the instability. This is similar to the Richardson number based subgrid-scale diffusion formulations of Lilly (1962, Tellus) and of Smagorinsky (1963), although their isotropic formulations have been simplified so as to only act on vertical gradients and perform diffusion in the vertical. This filter is completely local (2&Delta;z), diagnosing and acting only on adjacent grid cells, and is typically applied only in the stratosphere and above to avoid interference with physical dynamic instabilities in the free troposphere or boundary layer which are more accurately-simulated by a planetary boundary layer scheme or the resolved dynamics. This filter is applied at the same frequency that the physical parameterizations are updated.

We compute the local Richardson number on cell interfaces. Recall that k = 1 is the top layer of the domain and the index increases downward:

\f[
   Ri_{k-\frac{1}{2}} =  \frac{g\delta z\delta_z\theta_v}{(\theta^{k}_{v}+\theta^{k-1}_{v})((\delta_z u)^2+(\delta_z v)^2)}  \\  \tag {7.6}
  \f]


If Ri < 1, then mixing M is performed, scaled by Ri so that complete mixing occurs if R ≤ 0:

\f[
   M = max(1, (1-Ri)^2) \frac{\delta p^{*k}\delta p^{*(k-1)}}{\delta p^{*k} +\delta p^{*(k-1)}}  \\  \tag {7.7}
  \f]

The mixing is applied to each variable, including the winds interpolated to the A-grid (necessary for application of the physical paramterization tendencies) on a timescale &tau; (namelist parameter `fv_sg_adj`) which should be larger than the physics timestep to avoid suppressing resolved convective motions. The mixing is applied to the momentum (&delta;p<sup>*</sup> u<sub>a</sub>, &delta;p<sup>*</sup> v<sub>a</sub>), total energy, air mass, and all tracer masses, so that all of these quantities are conserved:

\f[
  \frac{\partial \phi^k}{\partial t} = - \frac{M}{\delta p^{*k} } \left(\phi^k - \phi^{k-1} \right) \frac{1}{\tau}  \\  \tag {7.8}
  \f]


\f[
  \frac{\partial \phi^{k-1}}{\partial t} = + \frac{M}{\delta p^{*k} } \left(\phi^k - \phi^{k-1} \right) \frac{1}{\tau}  \\  \tag {7.9}
  \f]

where  ϕ is a generic scalar. Note that since total energy and momentum are both conserved, lost kinetic energy automatically becomes heat.

This mixing is most useful for removing instabilities caused by vertically propagating waves near the top of the domain. The namelist variable `n_sponge` controls the number of levels at the top of the domain to which the filter is applied.

###7.4 Model-top sponge layer and energy-conserving Rayleigh damping

Two forms of damping are applied at the top of the domain to absorb vertically propagating waves nearing the upper boundary. The first is a diffusive sponge layer, which applies second-order damping to the divergence and to the vertical-momentum flux, and optionally also to the vorticity and mass fluxes if the hyperdiffusion is enabled. (This differs from earlier versions of FV<sup>3</sup>, which instead of adding explicit damping applied first-order upwind advection in the sponge layer, the strength of which is flow-dependent and not user-controllable.) The damping is computed in the same way as described earlier, although typically a very strong value of d<sub>2</sub> is used to ensure the vertically-propagating waves are sufficiently damped. The additional &nabla;<sup>2</sup> sponge-layer damping is applied to the top two layers of the model, with a weaker damping also applied in the third layer if d<sub>k2</sub> > 0.05. Since the model top is at a constant pressure, not constant height, it acts as a flexible lid, and therefore does not reflect as strongly as a rigid lid would.

The second form of damping is a Rayleigh damping, which damps all three components of the winds to zero with a timescale which depends on the pressure. Given a minimum timescale &tau;<sub>0</sub> and a cutoff pressure p<sub>c</sub> the damping timescale is:

\f[
    \tau(p^*) = \tau_0 sin \left( \frac{\pi}{2} \frac{log(p_c / p^*)}{log(p_c / p_T} \right)^2  \\  \tag {7.10}
  \f]

The strength of the applied damping is then determined by the magnitude of the cell-mean 3D vector wind U<sub>3D</sub>, including the vertical velocity, divided by a scaling velocity U<sub>0</sub>. The damping is only applied if the horizontal wind speed exceeds a latitude-dependent threshold (currently 25cos&theta;) or if the vertical velocity is larger than a given threshold. The damping is then applied, at the beginning of each large (physics) time step and before the Lagrangian dynamics is first called, by:

\f[
    u  \longleftarrow u(1 + \tau U_{3D} / U_0)^{-1} \\  \tag {7.11}
  \f]

The dissipated kinetic energy can then be restored as heat:

\f[
    T  \longleftarrow T + \frac{1}{2} U_{3D} * (1 - (1 + \tau U_{3D} / U_0)^{-2}) / C_v \\  \tag {7.12}
  \f]




