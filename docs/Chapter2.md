Cubed-sphere grid {#cube}
=========================================

##Chapter 2

GFDL will provide the additional documentation by the end of May 2021.

(This information is a reproduction of Appendix A from PL07)

In Section 2 of PL07 the flux-form multidimensional transport scheme is discretized in general non-orthogonal curvilinear coordinates. The covariant and contra-variant wind vector components are presented in Eqs. (4) and (5) (of Section 2 of PL07) based on the local unit vectors \f$ (\vec{e_1},\vec{e_2}) \f$ of the coordinate system. Given the angle \f$(\alpha)\f$ between the two unit vectors

\f[
 \cos \alpha = \vec{e_1} \cdot \vec{e_2}  \\ \tag {2.1}
  \f]

the covariant and contravariant components are related by the following relationships:

\f[
 u = \tilde{u} + \tilde{v} \cos \alpha  \\ \tag {2.2}
  \f]

\f[
 v = \tilde{v} + \tilde{u} \cos \alpha  \\ \tag {2.3}
  \f]

or (solving for the contravaraint components)

\f[
 \tilde{u} = \frac{1}{\sin^2 \alpha} [u - v \cos \alpha]  \\ \tag {2.4}
  \f]

\f[
 \tilde{v} = \frac{1}{\sin^2 \alpha} [v - u \cos \alpha]  \\ \tag {2.5}
  \f]

The winds on the cubed-sphere can be oriented to/from local coordinate orientation to a spherical latitude–longitude component form using the local unit vectors of the curvilinear coordinate system \f$(\vec{e_1},\vec{e_2})\f$ and the unit vector from the center of the sphere to the surface at the point of the vector location \f$(\vec{e_\lambda},\vec{e_\theta})\f$. Eqs.(2.6) and (2.7) represent the transformation from the spherical orientation \f$(u_{\lambda \theta},v_{\lambda \theta})\f$ to the local cubed-sphere form \f$(u,v)\f$ and the reverse transformation is presented in Eqs. (2.8) and (2.9).

\f[
 u = (\vec{e_1} \cdot \vec{e_\lambda}) u_{\lambda \theta} + (\vec{e_1} \cdot \vec{e_\theta}) v_{\lambda \theta}  \ \tag {2.6}
  \f]

\f[
 v = (\vec{e_2} \cdot \vec{e_\lambda}) u_{\lambda \theta} + (\vec{e_2} \cdot \vec{e_\theta}) v_{\lambda \theta}  \ \tag {2.7}
  \f]

\f[
 u_{\lambda \theta} = \frac{ (\vec{e_2} \cdot \vec{e_\theta}) u - (\vec{e_1} \cdot \vec{e_\theta}) v }{ (\vec{e_1} \cdot \vec{e_\lambda}) (\vec{e_2} \cdot \vec{e_\theta}) - (\vec{e_2} \cdot \vec{e_\lambda}) (\vec{e_1} \cdot \vec{e_\theta}) }  \\ \tag {2.8}
  \f]

\f[
 v_{\lambda \theta} = \frac{ (\vec{e_2} \cdot \vec{e_\lambda}) u - (\vec{e_1} \cdot \vec{e_\lambda}) v }{ (\vec{e_1} \cdot \vec{e_\lambda}) (\vec{e_2} \cdot \vec{e_\theta}) - (\vec{e_2} \cdot \vec{e_\lambda}) (\vec{e_1} \cdot \vec{e_\theta}) }  \\ \tag {2.9}
  \f]
