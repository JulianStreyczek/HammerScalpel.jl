using Test, HammerScalpel

function solveModelTest()
    
    # Solve model with standard parameters
    @time noint    = nopolicy(TT=2, ndims=3);                     # no intervention at all
    @time notest   = withpolicy("notest", TT=2, ndims=3);         # no testing 
    @time untarget = withpolicy("untargettest", TT=2, ndims=3);   # untargeted testing
    @time target   = withpolicy("targettest", TT=2, ndims=3);     # targeted testing
    @time isolate  = withpolicy("isolate", TT=2, ndims=3);        # isolation

    return("Success")

end

@test solveModelTest()=="Success"
