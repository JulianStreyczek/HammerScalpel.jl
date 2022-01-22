###############################################
#                                             #
# Four state model:                           #
# 1. Solve                                    #
# 2. Generate output for figures and tables   #
#                                             #
###############################################

"""
    withpolicy(policy; TT=52, ndims=20, thetai=0.38, thetas=0.0044)
    withpolicy(policy, nprocs; TT=52, ndims=20, thetai=0.38, thetas=0.0044)

Solve the model under a user-specified government policy. 
In all cases, the government quarantines individuals that are known to be infected.
The output is a [`Dict`](https://docs.julialang.org/en/v1/base/collections/#Dictionaries) containing the necessary variables and time-series for generating the figures. 

# Arguments:
- `policy::String`: Label for the government policy. Takes values
    - `notest`: no testing at all.
    - `untargettest`: random testing.
    - `targettest`: targeted testing under access to contact-tracing technology.
    - `isolate`: no testing, but isolation of certain fraction of population.
- `nprocs::Int64`: Optional. If specified, uses the module [Distributed](https://docs.julialang.org/en/v1/manual/distributed-computing/) for parallelizing the computations. *Requires prior setup of `nprocs` workers*.
- `TT::Integer`: Number of periods.
- `ndims::Integer`: Number of grid points in each dimension.
- `algorithm`: [NLopt](https://github.com/JuliaOpt/NLopt.jl) optimization algorithm.
- `thetai::Float`: Probability that infected person sends "infected" signal.
- `thetas::Float`: Probability that non-infected person sends "infected" signal.
"""
function withpolicy(policy::String; TT=52, ndims=20, thetai=0.38, thetas=0.0044)

    println("Solving model for policy = '$policy'. Counting iterations ($TT in total)...")
    
    # Par.ISO = 0, Par.TEST = 0: Model with no testing and isolation
    # Par.ISO = 0, Par.TEST = 1: Model with testing
    # Par.ISO = 1, Par.TEST = 0: Model with isolation and no testing

    ISO, TEST, theti, thets = [0, 0, 0, 0]; # may replace in next step
    algorithm = NLopt.LN_BOBYQA;            # may replace in next step
    
    if policy == "notest"
        ISO = 0;
        TEST = 0;
        theti = 1;
        thets = 1;
        algorithm = NLopt.LN_BOBYQA;
    elseif policy == "untargettest"
        ISO = 0;
        TEST = 1;
        theti = 1;
        thets = 1;
        algorithm = NLopt.LD_SLSQP;
    elseif policy =="targettest"
        ISO = 0;
        TEST = 1;
        theti = thetai;
        thets = thetas;
        algorithm = NLopt.LN_BOBYQA;
    elseif policy == "isolate"
        ISO = 1;
        TEST = 0;
        theti = thetai;
        thets = thetas;
        algorithm = NLopt.LD_SLSQP;
    else
        println("Invalid policy.")
    end

    # Parameters
    OPT   = 0;    # 0 for fmincon, 1 for global, 2 for interp
    beta  = .99;  # Discount factor
    gamma = 7/18; # Exit rate
    delta = .005; # Death rate
    sigma = 2;    # CES parameter on production function
    B1    = 1;    # Productivity of work sector
    B2    = .1;   # Productivity of home sector
    nu    = 1;    # Testing cost parameter
    c     = (1/.01)*50*52/63000 .* YY_noep(B1, B2, sigma, 1); # Testing cost parameter
    z     = 7500*52/63000 .* YY_noep(B1, B2, sigma, 1);       # Treatment cost
    d     = 15*3*52*(1/.96^(1/52)-1);                         # Death cost
    xi    = .8;           # Infected productivity 
    tautilde = 13/18*1/2; # Prob. of becoming symptomatic
    pi2   = .01;          # Infection rate at home
    pi1   = 1.5;          # Infection rate at work
    ut    = 1;
    # Combine into (immutable) NamedTuple
    Par = (ISO=ISO, TEST=TEST, OPT=OPT, beta=beta, gamma=gamma, delta=delta, sigma=sigma, B1=B1, B2=B2, 
    nu=nu, c=c, z=z, d=d, xi=xi, tautilde=tautilde, pi2=pi2, pi1=pi1, thets=thets, theti=theti, ut=ut)

    # Grids
    nSS = ndims; 
    nII = ndims;
    nRR = ndims;
    Icurve = -2;

    GridS = range(0, 1, length=nSS); #Susceptible
    GridR = range(0, 1, length=nRR); #Recovered
    GridI = 10 .^(range(Icurve, 0, length=nII)); #Infected but not known to be so
    GridI = GridI .- GridI[1];
    GridI[nII] = 1;

    Grid = (nS=nSS, nI=nII, nR=nRR, Icurve=Icurve, S=GridS, R=GridR, I=GridI);

    ## Solve model via backward induction
    rets = zeros(Int8, (TT,nSS,nII,nRR));
    ret_value = Dict(:SUCCESS => 1, :STOPVAL_REACHED => 2, :FTOL_REACHED => 3, :XTOL_REACHED => 4, :MAXEVAL_REACHED => 5, :MAXTIME_REACHED => 6,
    :FAILURE => -1, :INVALID_ARGS => -2, :OUT_OF_MEMORY => -3, :ROUNDOFF_LIMITED => -4, :FORCED_STOP => -5);
    ret_message = Dict("1" => :SUCCESS, "2" => :STOPVAL_REACHED, "3" => :FTOL_REACHED, "4" => :XTOL_REACHED, "5" => :MAXEVAL_REACHED, "6" => :MAXTIME_REACHED,
    "-1" => :FAILURE, "-2" => :INVALID_ARGS, "-3" => :OUT_OF_MEMORY, "-4" => :ROUNDOFF_LIMITED, "-5" => :FORCED_STOP);

    # Empty matrices for results
    V0 = zeros(TT,nSS,nII,nRR); #value function
    qq = zeros(TT,nSS,nII,nRR); #mass tested/isolated
    ll = zeros(TT,nSS,nII,nRR); #mass of unknown types assigned to work sector
    mm = zeros(TT,nSS,nII,nRR); #mass of known types assigned to work sector
    V00 = zeros(nSS,nII,nRR);
    qqg = zeros(nSS,nII,nRR);
    llg = zeros(nSS,nII,nRR);
    mmg = zeros(nSS,nII,nRR);
    V00T = zeros(nSS,nII,nRR);
    x0 = [0.0, 0.0, 0.0];

    # Value in last period
    for i in 1:nSS
        for j in 1:nII
            for k in 1:nRR
                S = Grid.S[i];
                I = Grid.I[j];
                R = Grid.R[k];
                V00T[i,j,k] = VTerm(Par,S,I,R);
            end
        end
    end

    # Backward induction
    for i in 1:TT
        println(i)
        if i>1
            V00[:,:,:] = V0[i-1,:,:,:];
            qqg[:,:,:] = qq[i-1,:,:,:];
            llg[:,:,:] = ll[i-1,:,:,:];
            mmg[:,:,:] = qq[i-1,:,:,:];
        else
            V00[:,:,:] = V00T;
        end
        for j in 1:nSS
            for k in 1:nII
                for l in 1:nRR
                    S = Grid.S[j];
                    I = Grid.I[k];
                    R = Grid.R[l];
                    if S+I>1.2 || S+I+R>1.2
                        V0[i,j,k,l] = NaN;
                        qq[i,j,k,l] = NaN;
                        ll[i,j,k,l] = NaN;
                        mm[i,j,k,l] = NaN;
                    else
                        lb = [0,0,0]; # [tau, lambda, mu]
                        if Par.ISO == 0 && Par.TEST == 0
                            ub = [0,1,1]; # tau=0 
                        else
                            ub = [1,1,1];
                        end
                        if i>1
                            x0 = [qqg[j,k,l],llg[j,k,l],mmg[j,k,l]];
                        else
                            x0 = [0,.9,.9];
                        end
                        opt = Opt(algorithm, 3);
                        # Local, derivative-based, finite differences
                        if algorithm == NLopt.LD_SLSQP
                            opt.min_objective = (x,g) -> Bellman_NLopt(Par, Grid, S, I, R, V00, x, g);
                        # Local, derivative-free
                        else
                            opt.min_objective = (x,g) -> Bellman(Par, Grid, S, I, R, V00, x);
                        end
                        opt.lower_bounds = lb;
                        opt.upper_bounds = ub;
                        opt.xtol_abs = 1.0e-8 .* ones(3);
                        opt.maxeval=2500;
                        (optf, optx, ret) = optimize(opt, x0)
                        # If local optimizer did not converge, perform global search
                        #if ret==:FORCED_STOP
                        #    # Local optimizer
                        #    local_opt = Opt(NLopt.LD_SLSQP, 3);
                        #    local_opt.xtol_abs = 1.0e-8 .* ones(3);
                        #    local_opt.maxeval=2500;
                        #    # Global optimizer
                        #    opt = Opt(NLopt.AUGLAG, 3);
                        #    opt.min_objective = (x,g) -> Bellman_NLopt(Par, Grid, S, I, R, V00, x, g);
                        #    #opt.min_objective = (x,g) -> Bellman(Par, Grid, S, I, R, V00, x);
                        #    opt.local_optimizer = local_opt;
                        #    opt.lower_bounds = lb;
                        #    opt.upper_bounds = ub;
                        #    opt.xtol_abs = 1.0e-8 .* ones(3);
                        #    opt.maxeval=2500;
                        #    #opt.population = 20;
                        #    (optf, optx, ret) = optimize(opt, x0)
                        #end
                        rets[i,j,k,l] = ret_value[ret];
                        qq[i,j,k,l] = optx[1];
                        ll[i,j,k,l] = optx[2];
                        mm[i,j,k,l] = optx[3];
                        V0[i,j,k,l] = -Bellman(Par,Grid,S,I,R,V00,optx);
                    end
                end
            end
        end
    end

    # Simulate optimal path
    I0 = 0.02;
    S0 = 1-I0/(1-Par.tautilde);
    R0  = 0;
    IT = zeros(TT);
    ST = zeros(TT);
    RT = zeros(TT);
    qT = zeros(TT);
    lamT = zeros(TT);
    muT = zeros(TT);
    Del = zeros(TT);
    VT = zeros(TT);
    YT = zeros(TT);
    L1 = zeros(TT);
    L2 = zeros(TT);
    C1T = zeros(TT);
    C2T = zeros(TT);
    Rknot = zeros(TT);
    #Test = zeros(TT);
    masstest = zeros(TT);
    mctest = zeros(TT);
    ctest = zeros(TT);
    IT[1] = I0;
    ST[1] = S0;
    RT[1] = R0;
    Del = 0;
    for i = 1:TT
        
        S = ST[i];
        I = IT[i];
        R = RT[i];
        IN = 1-S-R;
        
        if i<TT
            V00[:,:,:] = V0[TT-i,:,:,:];
        else
            V00[:,:,:] = V00T;
        end
        lb = [0,0,0];
        if Par.ISO == 0 && Par.TEST == 0
            ub = [0,1,1];
        else
            ub = [1,1,1];
        end
        if i>2
            x0 = [qT[i-1],lamT[i-1],muT[i-1]];
        else
            if Par.ISO == 0 && Par.TEST == 0
                x0 = [0,1,1];
            else
                x0 = [.3,.3,.3];
            end
        end
        opt = Opt(NLopt.LN_BOBYQA, 3);
        opt.min_objective = (x,g) -> Bellman(Par, Grid, S, I, R, V00, x);
        opt.lower_bounds = lb;
        opt.upper_bounds = ub;
        opt.xtol_abs = 1.0e-8 .* ones(3);
        opt.maxeval=2500;
        (optf, optx, ret) = optimize(opt, x0)
        qT[i] = optx[1];
        lamT[i] = optx[2];
        muT[i] = optx[3];
        VT[i] = -Bellman(Par,Grid,S,I,R,V00,optx);
        if Par.TEST == 1
            Del = qT[i].*Par.theti;
        else
            Del = 0;
        end
        if Par.ISO == 1
            LI1 = lamT[i].*((1-qT[i].*Par.thets).*S+(1-qT[i].*Par.theti).*I)+muT[i].*(1-Par.delta).*R;
            LI2 = (1-lamT[i]).*((1-qT[i].*Par.thets).*S+(1-qT[i].*Par.theti).*I)+(1-muT[i]).*(1-Par.delta).*R;
            if LI1>0
                Ipr1 = (lamT[i].^2 .*Par.pi1.*(1-qT[i].*Par.theti).*I)./LI1;
                RIpr1 = (lamT[i].^2 .*Par.pi1.*(1-qT[i].*Par.theti))./LI1;
            else
                Ipr1 = 0;
                RIpr1 = 0;
            end
            if LI2>0
                Ipr2 = ((1-lamT[i]).^2 .*Par.pi2.*(1-qT[i].*Par.theti).*I)./LI2;
                RIpr2 = ((1-lamT[i]).^2 .*Par.pi2.*(1-qT[i].*Par.theti))./LI2;
            else
                Ipr2 = 0;
            end
            Ipr = Ipr1 + Ipr2;
            Ipr = min(Ipr,1);
            if i<TT
                IT[i+1] = min(max((1-Par.tautilde).*(1-Par.gamma).*I+(1-qT[i]*Par.thets).*S.*Ipr,0),1);
                ST[i+1] = min(max((1-(1-qT[i].*Par.thets).*Ipr).*S,0),1);
                RT[i+1] = max(min(R + Par.gamma.*IN,1),0);
            end
        else
            LI1 = lamT[i].*(S+(1-Del).*I)+muT[i].*(1-Par.delta).*R;
            LI2 = (1-lamT[i]).*(S+(1-Del).*I)+(1-muT[i]).*(1-Par.delta).*R;
            if LI1>0
                Ipr1 = (lamT[i].^2 .*Par.pi1.*(1-Del).*I)./LI1;
                RIpr1 = (lamT[i].^2 .*Par.pi1)./LI1;
            else
                Ipr1 = 0;
                RIpr1 = 0;
            end
            if LI2>0
                Ipr2 = ((1-lamT[i]).^2 .*Par.pi2.*(1-Del).*I)./LI2;
                RIpr2 = ((1-lamT[i]).^2 .*Par.pi2)./LI2;
            else
                Ipr2 = 0;
                RIpr2 = 0;
            end
            Ipr = Ipr1 + Ipr2;
            Ipr = min(Ipr,1);
            if i<TT
                IT[i+1] = min(max((1-Del).*(1-Par.tautilde).*(1-Par.gamma).*I+S.*Ipr,0),1);
                ST[i+1] = min(max((1-Ipr).*S,0),1);
                RT[i+1] = max(min(R + Par.gamma.*IN,1),0);
            end
        end
        
        YT[i] = YY(Par,S,I,R,qT[i],lamT[i],muT[i]);
        C1T[i], C2T[i] = CC(Par,S,I,R,qT[i],lamT[i],muT[i]);
        if Par.ISO == 0
            L1[i] = (lamT[i].*(S+Par.xi.*(1-Del).*I)+muT[i].*(1-Par.delta).*R);
            L2[i] = ((1-lamT[i]).*(S+Par.xi.*(1-Del).*I)+(1-muT[i]).*(1-Par.delta).*R);
        else
            L1[i] = (lamT[i].*(S.*(1-qT[i].*Par.thets)+Par.xi.*(1-qT[i].*Par.theti).*I)+muT[i].*(1-Par.delta).*R);
            L2[i] = ((1-lamT[i]).*(S.*(1-qT[i].*Par.thets)+Par.xi.*(1-qT[i].*Par.theti).*I)+(1-muT[i]).*(1-Par.delta).*R);
        end
        if Par.TEST == 1
            mctest[i] = Par.c.*qT[i].*(Par.thets.*S+Par.theti.*I).*63000 ./52*1 ./YY_noep(Par,1);
            ctest[i] = Par.c.*(qT[i].*(Par.thets.*S+Par.theti.*I))^(1+Par.nu)./(1+Par.nu).*63000 ./52*1. /YY_noep(Par,1);
            masstest[i] = qT[i].*(Par.thets.*ST[i]+Par.theti.*IT[i]);
        else
            mctest[i] = 0;
            ctest[i] = 0;
            masstest[i] = qT[i].*(Par.thets.*ST[i]+Par.theti.*IT[i]);
        end
        if I>0
            Rknot[i]=(1 ./(Par.gamma)).*Ipr./(1-S-R);
        else
            Rknot[i] = 0;
        end
        
    end

    # Calculate output variables
    #cumdeaths = maximum(cumsum(Par.gamma.*Par.delta.*max.(1 .- ST .- RT, 0))); # Cum deaths
    #Ycumsum1 = YT[1]-YY_noep(Par,1);
    #for i = 2:TT
    #    Ycumsum1 = Ycumsum1 + Par.beta.^(i-1).*(YT[i]-YY_noep(Par,1));
    #end
    #Ycumsum1 = (1-Par.beta).*Ycumsum1; # Cum output loss rel to no ep. ss
    #compvar = 1/YY_noep(Par,1)*exp((1-Par.beta).*(log(YT[1]-C1T[1])-C2T[1])+Par.beta*(VT[2]));
    #compvar = 1-compvar; # Welfare loss rel to no ep. ss
    #S52 = ST[52]; # Susceptible at 52 weeks
    #I52 = 1-ST[52]-RT[52]; # Infected at 52 weeks
    #Imax = maximum(1 .-ST .-RT); #  Max number infected
    #R01 = Rknot[1]; # R0 at 1 week
    #R26 = Rknot[26]; # R0 at 26 weeks
    #R52 = Rknot[52]; # R0 at  52 weeks
    #println(ST[TT])
    #println(1-ST[TT]-RT[TT])
    #println(Rknot[TT])
    #println(Imax)
    #percentmaxeval = round(nmaxeval / (nsolve) * 100, digits=2);
    #println("MAXEVAL_REACHED $nmaxeval times ($percentmaxeval%).")
    #percentforcedstop = round(nforcedstop / (nsolve) * 100, digits=2);
    #println("FORCED_STOP $nforcedstop times ($percentforcedstop%).")

    nsolve = sum(rets .!= 0);
    for value in setdiff(unique(rets), [0])
        nvalue = sum(rets .== value);
        println("$(string(ret_message[string(value)])) $nvalue times ($(round(nvalue/nsolve*100, digits=2))%).")
    end

    return Dict("Par" => Par, "TT" => TT, "ST" => ST, "IT" => IT, "RT" => RT, 
    "YT" => YT, "C1T" => C1T, "Rknot" => Rknot, "masstest" => masstest, "mctest" => mctest)

end
    



