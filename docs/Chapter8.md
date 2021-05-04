Physics-dynamics coupling {#physics}
=========================================

##Chapter 8

(This information is a reproduction of Chapter 3 from LPH17)

###8.1 Staggered wind interpolation

The coupling to physical parameterizations and surface models is straight forward; the primary complication is interpolating between cell-centered, orthogonal winds used by most physics packages and the FV<sup>3</sup> staggered non-orthogonal D-grid. The unstaggered orthogonal wind is defined by horizontal components in longitude and latitude; the staggered non-orthogonal D-grid wind is defined by the winds tangential to the grid-cell interfaces by the horizontal covariant components in the cubed-sphere curvilinear components. A two-stage horizontal re-mapping of the wind is used when coupuling the physics packages to the dynamics. The D-grid winds are first transformed to locally orthogonal and unstaggered wind components at cell centers, as input to the physics. After the physics returns its tendencies, the wind tendencies (du/dt, dv/dt) are then remapped (using high-order algorithm) back to the native grid for updating the prognostic winds. This procedure satisfies the “no data no harm” principle — that the interpolation/remapping procedure creates no tendency on its own if there are no external sources from physics or data assimilation.

###8.2 Condensate loading and mass conservation

The mass &delta;m, and thereby also &delta;p<sup>*</sup>, in FV<sup>3</sup> is the total mass of both the dry air and of the water categories, including the vapor and condensate phases; the precise number N of water species is dependent upon the microphysics scheme used, or may be zero. This incorporates the effect of condensate loading into the air mass without a special parameterization for loading. The dry weight (per unit area) can be given as:

\f[
 g \delta m_d  =  \delta p^* \left( 1 - \sum_{m=1}^N q_m  \right)  =  \left( \delta p^* - \sum_{m=1}^N Q_m  \right) \\  \tag {8.1}
  \f]

where Q<sub>m</sub> = &delta;p<sup>*</sup> q<sub>m</sub> is the tracer mass. Dry mass should be conserved by the physical parameterizations; here we will assume this to be the case, so &delta;m<sub>d</sub> should be a constant in each grid cell during the physical tendency updates. The condition for dry mass conservation is then given by

\f[
    \delta p^{*(n+1)}  =   \delta p^{*n} + \delta \tau \sum_{m=1}^N \frac{dQ_m}{dt} = \delta p^{*n} \Delta M     \tag {8.2}        
\f] 


where &Delta;M = 1 + &delta;&tau; &Sigma;<sup>N</sup><sub>m=1</sub> dq<sub>i</sub>/dt. Physics packages usually return the rate of change in tracer mass  dQ<sub>m</sub> /dt, and so is independent of whether the solver uses total air mass or just dry air mass (as is done in many weather models). The tracer update is then done by:

\f[
    Q^{n+1}_m  =   Q^n_m + \delta \tau \frac{dQ_m}{dt}       \tag {8.3}        
\f] 


or, using (8.2)

\f[
    q^{n+1}_m  =  \left( Q^n_m + \delta \tau \frac{dq_m}{dt} \delta p^{*n}  \right) / \left( \delta p^{*(n+1)}  \right) \\  \tag {8.4}         
\f] 

\f[        
		      =	  \left( Q^n_m + \delta \tau \frac{dq_m}{dt}\delta p^{*n}  \right) / \left( \delta p^{*n} \Delta M \right) \\  \tag {8.5} 
\f] 

The full mass-conserving update algorithm is then:

\f[
    q^{*}_m  =   q^n_m + \delta \tau \frac{dq_m}{dt}       \tag {8.6}        
\f] 

\f[
    \Delta M  =   1 +  \delta \tau \sum_{m=1}^N \frac{dq_m}{dt}      \tag {8.7}        
\f] 

\f[
      \delta p^{n+1}  = \delta p^{n} \Delta M   \tag {8.8}        
\f] 

\f[
      \delta q^{n+1}_m  =  q^{*}_m / \Delta M   \tag {8.9}        
\f] 

Note that often the mass of non-water species, such as ozone or aerosols, are considered so small that they are not included in &delta;m; however, since their mixing ratio is still the quotient of the tracer mass and the total air mass, if the effects of water species are included in the total air mass their mixing ratios must still be adjusted by (8.9).



