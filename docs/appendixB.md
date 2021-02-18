Variables and notation {#variables}
===========
(This information is a reproduction of Appendix B excluding B.2 and B.3 from LPH17)

Variables  |        Notation
------------- | -------------
`u,v`              | (*) D-grid winds 
` w`               | (*) Explicit vertical velocity 
 &delta; p<sup>*</sup>  | (*) Layer hydrostatic pressure thickness, proportional to mass 
 &Theta; <sub>v</sub>   | (*) Virtual potential temperature
 &delta; <sub>z</sub>   |  (*) Layer geometric depth
&rho;  | Total (air and all water species) mass density, equal to   -&delta;p<sup>*</sup>/g&delta;z
 Q<sub>i</sub>  | (*) Density of tracer i (also written Q as a generic tracer density)
q<sub>i</sub>   | Mixing ratio of tracer i, defined with respect to total air mass; equal to  Q<sub>i</sub>/&delta;p<sup>*</sup>
q<sub>v</sub>   | Mixing ratio of water vapor    
p  | Total cell-mean pressure
p<sup>*</sup>   | Cell-mean hydrostatic pressure
p'   | Cell-mean nonhydrostatic pressure component, equal to p - p<sup>*</sup>
&Delta;A, &Delta;x, &Delta;y | D-grid cell areas and cell face lengths
&Delta;A<sub>c</sub>, &Delta;x<sub>c</sub>, &Delta;y<sub>c</sub>  |   Dual-grid cell areas and cell face lengths
 &Delta;t  | Lagrangian dynamics (or acoustic) time step
&Delta;T |  Vertical remapping interval
&Delta;&tau;  | Physics time step
c<sub>pd</sub>  | Specific heat of dry air at constant pressure 
R<sub>d</sub>  | Gas constant for dry air 
c<sub>p</sub>  | Variable specific heat of moist air 
&kappa;  | R/c<sub>p</sub>, where R is the (variable) gas constant of moist air 
i, j, k  | Spatial grid-cell indices for the local x-, y-, and z-directions, respectively, given as subscripts
n   | time index
m  | tracer index, given as a subscript
n  |  time index, given as a superscript: t<sup>n</sup>  = n&Delta;t 
 
Here, a (*) indicates a prognostic variable. All variables are cell-means (or face-means for horizontal winds) unless otherwise noted. The differencing notation used in this document follows that of LR96, LR97, and L04, in which the operator &delta;<sub>x</sub> &phi;  is defined as a centered-difference operator:

\f[
(\delta_x \phi)_{i+1/2} = \phi_{i+1} - \phi_i  \tag {B.1}
  \f]

  The indices on dependent variables are suppressed unless explicitly needed.) This definition differs from the similar operators of in the literature intended to be second-order discretizations of a derivative; to do this with our definition of &delta;<sub>x</sub> a 1/&Delta;<sub>x</sub> term would be needed to complete the discrete derivative.

========================

##B.1 Important Relations

Cell-mean pressure:
\f[
p^* =  \frac{\delta p^*}{\delta \log p^*} \; \mathrm{(hydrostatic)} \\. \tag {B.2}
  \f]
  
\f[
p =  \rho R_d T_v = - \frac{R_d}{g} \frac{\delta p^* T_v}{\delta z}  \mathrm{(nonhydrostatic)} \tag {B.3}
  \f]

Conversion from temperature *T* to virtual potential temperature &Theta;<sub>v</sub> :

\f[
T_v = T \left ( 1 + \epsilon q_v \right )  \tag {B.4}
  \f]
  
\f[
\Theta_v  = T_v / p^\kappa    \tag {B.5}
  \f]
where &epsilon; = 1 + R<sub>v</sub>/R<sub>d</sub> and &kappa; = R/c<sub>p</sub>. Note that we do *not* include the arbitrary constant factor p<sub>0</sub><sup>&kappa;</sup> in our definition of potential temperature; our form is the equivalent to setting p<sub>0</sub> = 1Pa but simpler and more computationally efficient.




 





