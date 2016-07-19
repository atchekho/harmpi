# HARMPI exercises by Sasha Tchekhovskoy

## Hydro problems

To run problems with HARMPI and to analyze the results, please follow this [tutorial](tutorial.md).

### 1D hydro problems

* Bondi accretion

	Set `WHICHPROBLEM` to `BONDI_PROBLEM_1D` in [decs.h](decs.h). Note: a good total resolution is 256x1x1

	* Plot the profiles of density at a few times in a simulation. Determine the position of the sonic surface. 

## Magnetized problems

### 1D magnetized problems

* Monopole problem

	Set `WHICHPROBLEM` to `MONOPOLE_PROBLEM_1D` in [decs.h](decs.h). Note: a good total resolution is 768x1x1

    * Measure the ratio OmegaF/OmegaH. How does it compare with the standard value, 0.5?
    * What is the initial magnetization (near the black hole) of the plasma? Hint: look at the quantity $\sigma_0 = b^2/4\pi\rho c^2 =$ `bsq/rho` . Here $b^2/4\pi =$ `bsq` is the square of the fluid frame magnetic field and $\rho =$ `rho` is the fluid frame mass density.
    * Make a plot of Lorentz factor vs. radius, $\Gamma(r)$. Hint: $Gamma =$ `alpha * uu[0]`. What value does the Lorentz factor saturate at? How does it compare to the near-black hole magnetization, $\sigma_0$, of the magnetosphere determined above?

### 2D magnetized problems

* BZ-Michel monopole problem

	Set `WHICHPROBLEM` to `BZ_MONOPOLE_2D` in [decs.h](decs.h). Note: a good total resolution is 256x256x1

    * Plot the location of the surface in $R$--$z$ plane at which the radial contravariant component of velocity, $u^r$, vanishes. Hint: plot the contour of `uu[1] == 0`.
    * Plot the dependence of the ratio OmegaF/OmegaH at the horizon as a function of the angle, $\theta$.
    * Compare the prediction of power by the standard Blandford & Znajek (1977) power formula to the simulation results.

* 2D monopole problem

	Set `WHICHPROBLEM` to `MONOPOLE_PROBLEM_2D` in [decs.h](decs.h).  Note: a good total resolution is 1152x256x1. It will take some time for the problem to complete depending on the number of cores used.

    * Compare $\Gamma(r)$ to that in the above 1D monopole problem. In which case is the Lorentz factor higher? Why?
    * Verify that the value of the Lorentz factor at the fast surface is $\Gamma =$ `(bsq/rho)**0.5`. What is the value of $\Gamma$ at the fast surface?

* Torus problem

	Set `WHICHPROBLEM` to `TORUS_PROBLEM` in [decs.h](decs.h).  Note: a good total resolution is 256x256x1. It will take some time for the problem to complete depending on the number of cores used.

	* Check by how many cells the MRI wavelength is resolved in the initial conditions. Good resolution is $\gtrsim15$ cells per wavelength, but you could sometimes also get away with $\gtrsim 5{-}10$.
    * Make a movie of the simulation: logarithm of density shown with color contours overlaid with magnetic field lines. Feel free to ask for guidance on how to make the movie.
