###############################################
#                                             #
# Main file for starting solving functions    #
# and generating figures                      #
#                                             #
###############################################

# Load packages
using NLopt
using FiniteDifferences
using GridInterpolations
using Plots
using JLD2

# Include main functions
include("nopolicy.jl")
include("withpolicy.jl")
include("figures_paper.jl")

# Include helper functions
include("helperfunctions.jl")

# Solve model with standard parameters
@time noint    = nopolicy(TT);                 # no intervention at all
@time notest   = withpolicy("notest");         # no testing 
@time untarget = withpolicy("untargettest");   # untargeted testing
@time target   = withpolicy("targettest");     # targeted testing
@time isolate  = withpolicy("isolate");        # isolation

# Save important objects
#save("src/data/calib1_nointervention.jld2", noint);
#save("src/data/calib1_notest_small.jld2", notest);
#save("src/data/calib1_untarget_small.jld2", untarget);
#save("src/data/calib1_target_small.jld2", target);
#save("src/data/calib1_isolate_small.jld2", isolate);

# Figures (using freshly-generated data)
createFig5(noint, notest, untarget)
createFig6(notest, untarget, target)
createFig7(untarget, target, isolate)

# Figures (using stored data)
createFig5()
createFig6()
createFig7()



