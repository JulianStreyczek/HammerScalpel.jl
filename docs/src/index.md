![header](./assets/header.png)

# Replication of Chari, Kirpalani, and Phelan (2021)

> This project was conducted as a final assignment for the PhD course [Numerical Methods](https://floswald.github.io/NumericalMethods/) at Bocconi University in Fall 2021.

This package replicates the main figures of Chari, Kirpalani, and Phelan (2021): "*The hammer and the scalpel: On the economics of indiscriminate versus targeted isolation policies during pandemics*", published in the Review of Economic Dynamics. The paper is available [here](https://doi.org/10.1016/j.red.2020.11.004), and the original replication material in Matlab [here](https://ideas.repec.org/c/red/ccodes/20-237.html).

The authors develop a theoretical model that combines epidemic transmission and economic outcomes. They study the effect of different government responses, such as quarantining, testing, contact-tracing, and isolation. 

## Installation

This package is not listed as an official Julia package. If you wish to use it, you may clone [this Github repository](https://github.com/JulianStreyczek?tab=repositories) to your machine. You then have two options:

**Use as a package:** Start up Julia and go into the package editor by typing `]`. Then, type

    activate .
    instantiate 

After pressing backspace, go back into Julia's standard command mode and type
    
    using HammerScalpel

The functions should now be ready to go.

**Use as ordinary code:** Play around with the files `main.jl` or `main_parallel.jl`, which run the code in sequential and parallelized form, respectively. 

## Using the package

If you wish run the replication from beginning to end, simply type

    solveModel()

This will call all necessary functions in the correct order. 
Specifically, it first computes the no-intervention outcome using [`nopolicy`](@ref), and the respective outcome under the various policies using [`withpolicy`](@ref).

Since the computations may take a while, you have the option to run them more efficiently on multiple CPU cores at the same time. 
To do so, type

    solveModel(nprocs)

where `nprocs` specifies the number of logical cores on your CPU. 

If you just want to check my results without running the simulations yourself, you may use stored outcomes for default parameters for creating the figures by running

    createFig5()
    createFig6()
    createFig7()

## Documentation

```@docs
solveModel
nopolicy
withpolicy(::String)
createFig5(::Dict, ::Dict, ::Dict)
createFig6(::Dict, ::Dict, ::Dict)
createFig7(::Dict, ::Dict, ::Dict)
```

## Implementation details

**Optimization algorithm:** The authors originally use the algorithm _SQP_ for maximizing the Bellman equation in each iteration, which in this case requires approximating gradients and Hessian.
For some policies, the [NLopt](https://github.com/JuliaOpt/NLopt.jl) algorithm _LN\_BOBYQA_ finds the same solution with greater speed.
Therefore, although the replication results hold when using _SQP_ everywhere, I use _LN\_BOBYQA_ when applicable to improve speed.

**Dimensionality:** The authors originally use 40 grid points for each control variable. In this replication, I reduced them to 20 for computational efficiency, which is why some of the time series look a bit ragged. 

**Inconsistencies:** In my replications of Figures 6 and 7, the time series for R~0~ features a sharp increases for policies "targeted" and "isolate" around T=40 and T=50, respectively, which are not present in the published paper. I verified that these inconsistencies are also present in the original Matlab code.