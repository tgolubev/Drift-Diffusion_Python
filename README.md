# Drift-Diffusion_models

Here is a 1D model written in Python which solves the semiconductor Poisson-Drift-Diffusion equations using finite-differences. This models simulates a solar cell under illumination, but can be adapted to other semiconductor devices as well. It can be modified to solve other systems (i.e. through changing the boundary conditions, adding recombination rates, and modifying the generation rate). 

The equations are solved using the self-consistent iterative approach called the Gummel method. In order to ensure numerical stability for the continuity equations, Scharfetter Gummel discretization as well as linear mixing of old and new solutions is used. 


## C++ and Matlab implementations

You can find my C++ and Matlab implementations of this same model as well as 2D and 3D versions here: 

\url{https://github.com/tgolubev/Drift-Diffusion_models-Cpp_Matlab}
