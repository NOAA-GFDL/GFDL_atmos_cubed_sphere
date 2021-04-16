Grid refinement techniques {#grid}
=========================================

##Chapter 9

(This information is a reproduction of Chapter 4 from LPH17)

There is a need for increasingly high-resolution numerical models for weather and climate simulation, but also an increasing need for coupling the newly resolved scales to the large and global-scale circulations, for which limited area models are only of limited use. However, uniformly-high resolution global models are not always practical on present-day computers. The solution to this problem is to locally refine a global grid, allowing for enhanced resolution over the area of interest while also representing the global grid. FV<sup>3</sup> has two variable-resolution methods: a simple Schmidt transformation for grid stretching, and two-way regional-to-global nesting. These methods can be combined for maximum flexibility.

FV<sup>3</sup> can also be configured as a doubly-periodic solver, in which the cubed-sphere is replaced by a Cartesian-coordinate doubly-periodic horizontal grid; otherwise the solver is unchanged. This can be useful for idealized simulations at a variety of resolutions, including very high horizontal resolutions useful for studying explicit convection.

###9.1 Grid stretching
Here we follow the development of HLT16. A relatively simple variable resolution grid can be created by taking the original cubed-sphere grid and applying the transformation of F. Schmidt (Beitr. Atmos. Phys., 1977) to “pull” grid intersections towards a “target” point, corresponding to the center point of the high-resolution region. This is done in two steps: the grid is distorted towards the south pole to get the desired degree of refinement, and then the south pole is rotated to the target point using a solid-body rotation. Distorting to the south pole means that the longitudes of the points are not changed, only the latitudes, greatly simplifying the transformation.

The transformation of the latitude &theta; to &thetasym; is given by:

\f[
 	sin\vartheta  =  \frac{D + sin\theta}{1 + Dsin\theta} \\  \tag {9.1}
  \f]

where the distortion is a function of the stretching factor c, which can be any positive number:

\f[
 	D  =  \frac{1 - c^2}{1 + c^2} \\  \tag {9.2}
  \f]

Using c = 1 causes no stretching. Note that other forms for the transformation could also be used without making any other changes to the solver.

Although the grid has been deformed, the solver still uses the assumption that the grid cells are bounded by great-circle arcs, which are not strictly identical to a Schmidt transformation of the cubed-sphere arcs of the unstretched grid.

###9.2 Grid nesting

Using grid nesting can greatly increase the flexibility of grid refinement in the model, at the cost of greater complexity in the solver. The major strength of grid nesting is its ability to use independent configurations on each grid, including different time steps and physical parameterizations, most appropriate for that particular grid. The ability to use a longer time step on the coarse grid than on the nested grid can greatly improve the efficiency of a nested-grid model; and choosing parameterizations independently allows values appropriate for each resolution without needing compromise or “scale-aware” parameterizations.

Here we follow the development of HL13, with additional updates necessary for the nonhydrostatic solver. Implementing two-way grid nesting involves two processes: interpolating the global grid variables to create boundary conditions for the nested-grid, and then updating the coarse-grid solution with nested-grid data in the region they overlap. The goal is to do so in as efficient of a manner consistent with the finite-volume methodology.

A major feature of FV<sup>3</sup>’s nesting is to use concurrent nesting, in which the nested and coarse grids run simultaneously, akin to how coupled models run their atmosphere and ocean components at the same time on different sets of processors. This can greatly reduce the amount of load imbalance between the different processors.

The entire nesting cycle is as follows, starting at the beginning of call to the solver:

- For each `p_split` step:
	- Call solver
	- Fetch boundary condition data from coarse grid
	- In Lagrangian dynamics, update boundary conditions at each ¢t by extrapolating from two earlier coarse-grid states.
	- Perform tracer transport and vertical remapping
	- Perform two-way update
- Call physics

Note that we do not do a compile cycle every coarse-grid time step, unlike many regional nested-grid models. The cycling can be carried out multiple times per physics time step, if more frequent updates of the boundary conditions and of the two-way communication are considered necessary. There is also an option to perform the last two-way update after the physics, instead of before, which changes how the physical parameterizations interact with the nested-grid solution passed to the coarse grid. Performing the update before calling the physics has been found to yield better results in real-data forecasts.

Currently, nested grids in FV<sup>3</sup> are constrained to be a proper refinement of a subset of coarse-grid cells; that is, each coarse-grid cell in the nested grid region is subdivided into N x N nested-grid cells. This greatly simplifies the nested-grid boundary condition interpolation and the two-way updating. Nested grids are also static and constrained to lie within one coarse-grid face. However, the algorithm does not require an aligned, static grid in one cube face, and any of these conditions may be relaxed in the future.

The nested-grid boundary conditions are implemented in a simple way. Coarse-grid data is interpolated from the coarse grid into the halo of the nested grid, thereby providing the nested-grid boundary conditions. Linear interpolation, although it is simple and and is not conserving, does have the advantage of not introducing new extrema in the interpolated field. The boundary conditions for staggered variables are interpolated directly from the staggered coarse grids. Boundary conditions are needed for each prognostic variable, including the tracers; also, boundary conditions are needed for the C-grid winds, available at each half-time step, and for the divergence when using fourth or higher-order divergence damping.

Finally, boundary conditions for the layer-interface nonhydrostatic pressure anomalies are also needed to evaluate the pressure-gradient force. Instead of interpolating these interface values from the coarse grid, they are instead diagnosed and interpolated from the other boundary condition variables using the same methods as the semi-implicit solver.

Most nested-grid models perform time-interpolation between two coarse grid states on each time step, but since the grids are integrated concurrently in FV<sup>3</sup>, interpolation is not possible. Instead, we can extrapolate between two earlier coarse-grid states. If interpolated coarse-grid variables are available at times t and t - &Delta; &tau;, where &Delta; &tau; = N &Delta;t, then the extrapolation for a given variable &phi; at time t + n&delta;t  (n=1,...,N) is given by:

\f[
 	\phi^{t+n\delta t}  = \left( 1 + \frac{n}{N} \right) \phi^t - \frac{n}{N} \phi^{t-\Delta \tau} \\  \tag {9.3}
  \f]

The extrapolation is constrained for positive-definite scalars so that the value of the boundary condition at t + &Delta;&tau; is non-negative, which is done by the substitution  &phi;<sup>t - &Delta;&tau; </sup>   &rarr; min(&phi;<sup>t - &Delta;&tau; </sup>, 2&phi;<sup>t</sup> ).

Two-way updates from the nested to the coarse grid are performed consistent with the finite-volume numerics. Scalars are updated to the coarse grid by performing an average of nested-grid cells, consistent with the values being cell-averages. The staggered horizontal winds are updated by averaging the winds on the faces of nested-grid cells abutting the coarse-grid cell being updated, so that the update preserves the average of the vorticity on the nested-grid cells. In FV<sup>3</sup> only the three wind components and the temperature is updated to the coarse grid; the air and tracer masses are not updated, trivially conserving their masses on the nested grid, and reducing the amount of noise created through overspecification of the solution on the coarse grid. Since the air mass determines the vertical coordinate, which will differ between the two grids, the averaged nested-grid data is remapped onto the coarse-grid’s vertical coordinate.


