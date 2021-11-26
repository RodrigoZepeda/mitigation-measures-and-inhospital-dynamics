# Mitigation measures and inhospital dynamics

Code to reproduce results from chapter (TBA). 

## Running code

To get the results from the paper run `run_simulations.R` which will download the data from the online repository, format it, fit the `Stan` model and evaluate the different scenarios. 

### Auxiliary functions

In order for `Stan` to work in my laptop I needed to find the starting points for the parameters. The julia function `fitmodelserver.jl` performs optimization of the parameters and generates valiable initial points. I copied them to `R` but that's where they come from. 