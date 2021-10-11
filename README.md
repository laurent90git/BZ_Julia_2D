# BZ in Julia

Simulate the two-dimensional Belousov-Zhabotinsky oscillating reactions model, which is a very stiff reaction-diffusion system.

For the integration with Julia, a Strang splitting is implemented, using ROCK4 for the integration of the diffusive part, and ROCK4 (Radau5 or another L-stable method in the future) for the reactive part.

The post-processing is currently done in Python.

This script is heavily inspired by the work of ???? and others at CMAP, Ecole Polytechnique, France.

## Usage

Run the following commands

```bash
cd("your directory")
# Load the model and runs the simulation
include("test.jl")
# Export the solution
include("cusIO".jl")
export_sol_2dcart("save_test_singlestep.txt", sol, (ncol,nlin,dx,dy,xmin,ymin,xmax,ymax))
```
