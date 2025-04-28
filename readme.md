# Parameter optimization and model simulations

This repository contains MATLAB source code, supporting functions, and example datasets for modeling microbial growth dynamics in a batch/chemostat setting, simulating virtual experiments, fitting parameters to data, and evaluating model performance in association with "Metabolic interplay drives population cycles in a cross-feeding microbial community."

## Getting Started

### Dependencies

* Microsoft Windows (tested on Windows 10)
* MATLAB (tested on R2024a)
* Optimization Toolbox (for ```fmincon()```)

### Installing

* Install MATLAB from the [MathWorks](https://www.mathworks.com/products/matlab.html) website.
* Select the optimization toolbox as an addon during installation.
* Download or clone this repository.
* Place all .m, .mat, and .xlsx files into the same working directory or add the folder to your MATLAB path.
* Typical installation takes Less than 10 minutes on a standard desktop computer.

### Demo
```
lambda_optimization.m
```
* Fits model parameters to the dataset across different regularization strengths.
* Expected run time: > 1 hour depending on machine specifications and script parameters (number of samples from LHS, size of regularization vector)
```
Plot_Experimental_Data_Fit.m
```
* Plots how well the fitted model matches the experimental data.
* Expected run time: < 10 seconds

## Instructions for use
Optimizing Parameters on New Data:
* Prepare a new Excel sheet formatted similarly to ```20221227_Data.xlsx```.
* Update the ```fit_data``` path in ```Lambda_Optimization.m```.
* Run ```Lambda_Optimization.m``` to generate new optimized parameter sets.

Plotting with Custom Parameter Sets:
* Load your custom theta file and run ```Plot_Experimental_Data_Fit.m``` to generate updated plots.

Simulating New Conditions:
* Modify the initial conditions and parameters inside ```Simulate_Chemostat.m``` or ```Simulate_Chemostat_minimal.m```.
* Run the script to generate new simulation results.

Reproducing Results:
* To reproduce the Fig. 1d in the manuscript: Use ```Plot_Experimental_Data_Fit.m``` to visualize fit performance with the provided optimized parameter set ```theta_full.mat```.

## Authors

Tyler Ross

## Version History

* 0.1
    * Initial Release

## License

This project is licensed under the MIT License. You are free to use, modify, and distribute the code provided proper citation to the original work is given.