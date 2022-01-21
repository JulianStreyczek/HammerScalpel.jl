
"""
    solveModel(; TT=52, ndims=20, thetai=0.38, thetas=0.0044)
    solveModel(nprocs:Int64; TT=52, ndims=20, thetai=0.38, thetas=0.0044)

Solve the model for all possible policy choices and replicate Figures 5, 6, and 7 in Chari, Kirpalani, and Phelan (2021). 
The function calls [`nopolicy`](@ref) and [`withpolicy`](@ref) with the relevant keywords, and uses the outputs to generate the figures using 
[`createFig5`](@ref), [`createFig6`](@ref), and [`createFig7`](@ref).    

# Arguments:
- `nprocs::Int64`: Optional. If specified, uses the module [Distributed](https://docs.julialang.org/en/v1/manual/distributed-computing/) for parallelizing computations. 
- `TT::Integer`: Number of periods.
- `ndims::Integer`: Number of grid points in each dimension.
- `thetai::Float`: Probability that infected person sends "infected" signal.
- `thetas::Float`: Probability that non-infected person sends "infected" signal.
"""
function solveModel(; TT=52, ndims=20, thetai=0.38, thetas=0.0044)
    println("Solving the model under no policy, and under policies 'notest', 'untarget', 'target', 'isolate'. This may take a while...")
    
    # Solve model with standard parameters
    @time noint    = nopolicy(TT);                                            # no intervention at all
    @time notest   = withpolicy("notest"; TT, ndims, thetai, thetas);         # no testing 
    @time untarget = withpolicy("untargettest"; TT, ndims, thetai, thetas);   # untargeted testing
    @time target   = withpolicy("targettest"; TT, ndims, thetai, thetas);     # targeted testing
    @time isolate  = withpolicy("isolate"; TT, ndims, thetai, thetas);        # isolation

    # Figures
    println("Generating figures...")
    createFig5(noint, notest, untarget)
    createFig6(notest, untarget, target)
    createFig7(untarget, target, isolate)
end



##########################################################################################
function solveModel(nprocs::Int64; TT=52, ndims=20, thetai=0.38, thetas=0.0044)
    println("Solving the model under no policy, and under policies 'notest', 'untarget', 'target', 'isolate'. This may take a while...")

    # Tell Julia to use several processes for distributed computing
    addprocs(nprocs)

    @everywhere begin
        # Load packages
        @eval using Pkg; Pkg.activate(".") # Activate environment so that workers "see" packages
        @eval using HammerScalpel
        
    end

    # Solve model with standard parameters
    @time noint    = nopolicy(TT);                                                   # no intervention at all
    @time notest   = withpolicy("notest", nprocs; TT, ndims, thetai, thetas);        # no testing 
    @time untarget = withpolicy("untargettest", nprocs; TT, ndims, thetai, thetas);  # untargeted testing
    @time target   = withpolicy("targettest", nprocs; TT, ndims, thetai, thetas);    # targeted testing
    @time isolate  = withpolicy("isolate", nprocs; TT, ndims, thetai, thetas);       # isolation

    # Figures
    println("Generating figures...")
    createFig5(noint, notest, untarget)
    createFig6(notest, untarget, target)
    createFig7(untarget, target, isolate)

end
