# Drift-Diffusion_models

Here is a 1D model written in Python which solves the semiconductor Poisson-Drift-Diffusion equations using finite-differences. This models simulates a solar cell under illumination, but can be adapted to other semiconductor devices as well. It can be modified to solve other systems (i.e. through changing the boundary conditions, adding recombination rates, and modifying the generation rate). 

The equations are solved using the self-consistent iterative approach called the Gummel method. In order to ensure numerical stability for the continuity equations, Scharfetter Gummel discretization as well as linear mixing of old and new solutions is used. 

### Performance
The code has been accelerated using Numba @jit decorators. Sample CPU Times:
Without Numba: 469.7 sec
With Numba: 73.7 sec

The conclusion here is that Numba gives a significant performance boost with low effort. You may read about Numba here:
http://numba.pydata.org/

## C++ and Matlab implementations

You can find my C++ and Matlab implementations of this same model as well as 2D and 3D versions here: 

https://github.com/tgolubev/Drift-Diffusion_models-Cpp_Matlab

### Performance comparison:

For the 1D code with mesh size dx = 0.25nm and a system size of 300nm:

Python: 69.8 sec
Matlab: 40 sec
C++: 3.7 sec

Therefore, currently the C++ version is much faster, perhaps with a disadvantage of being less elegant to read.

## Citing
If you use this code in a scientific publication, I would appreciate citations to the following paper:

T. Golubev, D. Liu, R. Lunt, P. Duxbury. Understanding the impact of C60 at the interface of perovskite solar cells via drift-diffusion modeling. AIP Advances 9, 035026 (2019) https://aip.scitation.org/doi/full/10.1063/1.5068690

