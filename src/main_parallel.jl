###############################################
#                                             #
# Main file for starting solving functions    #
# and generating figures,                     #
# parallelized with 'Distributed'             #
#                                             #
###############################################

# Load parallel computing package 
using Distributed
using SharedArrays
using NLopt
using FiniteDifferences
using GridInterpolations
using JLD2
using Plots

# Include main functions
include("nopolicy.jl")
include("withpolicy_parallel.jl")
include("figures_paper.jl")

# Include helper functions
include("helperfunctions.jl")

# Tell Julia to use several processes for distributed computing
nprocs = 8;
addprocs(nprocs)

@everywhere begin
    # Load packages
    using Pkg; Pkg.activate(".") # Activate environment so that workers "see" packages
    using Distributed
    using SharedArrays
    using NLopt
    using FiniteDifferences
    using GridInterpolations

    # Include helper functions
    include("helperfunctions.jl")
end

# Solve model with standard parameters
@time noint    = nopolicy();                          # no intervention at all
@time notest   = withpolicy("notest", nprocs;);        # no testing 
@time untarget = withpolicy("untargettest", nprocs;);  # untargeted testing
@time target   = withpolicy("targettest", nprocs;);    # targeted testing
@time isolate  = withpolicy("isolate", nprocs;);       # isolation

# Save important objects
#save("src/data/calib1_nointervention.jld2", noint);
#save("src/data/calib1_notest_small.jld2", notest);
#save("src/data/calib1_untarget_small.jld2", untarget);
#save("src/data/calib1_target_small.jld2", target);
#save("src/data/calib1_isolate_small.jld2", isolate);

# Load objects
#noint = load("src/data/calib1_nointervention.jld2");
#notest = load("src/data/calib1_notest_small.jld2");
#untarget = load("src/data/calib1_untarget_small.jld2");
#target = load("src/data/calib1_target_small.jld2");
#isolate = load("src/data/calib1_isolate_small.jld2");

# Figures (using freshly-generated data)
createFig5(noint, notest, untarget)
#savefig("src/figures/Fig5_dim20")
createFig6(notest, untarget, target)
#savefig("src/figures/Fig6_dim20")
createFig7(untarget, target, isolate)
#savefig("src/figures/Fig7_dim20")

# Figures (using stored data)
#createFig5()
#createFig6()
#createFig7()

