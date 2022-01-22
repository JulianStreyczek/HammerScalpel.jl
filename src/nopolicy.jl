####################################################
#                                                  #
# Compute the no policy / no intervention outcome  #
#                                                  #
####################################################

"""
    nopolicy(;TT=52)
Simulate the model in Chari, Kirpalani, and Phelan (2021) for `TT` periods under no government response. 
The output is a [`Dict`](https://docs.julialang.org/en/v1/base/collections/#Dictionaries) containing 
the necessary variables and time-series for generating the figures. 
"""
function nopolicy(;TT=52)
    
    # Parameters
    ISO = 0;
    TEST = 1;
    OPT = 0; #0 for fmincon, 1 for global, 2 for interp
    beta = .99; # Discount factor
    gamma = 7/18; # Exit rate
    delta = .005; # Death rate
    sigma = 2; # CES parameter on production function
    B1 = 1; # Productivity of work sector
    B2 = .1; # Productivity of home sector
    nu = 1; # Testing cost parameter
    c = (1/.01)*50*52/63000 .* YY_noep(B1, B2, sigma, 1); # Testing cost parameter
    z = 7500*52/63000 .* YY_noep(B1, B2, sigma, 1); # Treatment cost
    d = 15*3*52*(1/.96^(1/52)-1); # Death cost
    xi = .8; # Infected productivity los
    tautilde = 13/18*1/2; # Prob. of becoming symptomatic
    pi2 = .01; # Infection rate at home
    pi1 = 1.5; # Infection rate at work
    thets = .0044; #Signal probability from susceptible person
    theti = .38; #Signal probability from infected person
    ut = 1;
    
    # Combine into (immutable) NamedTuple
    Par = (ISO=ISO, TEST=TEST, OPT=OPT, beta=beta, gamma=gamma, delta=delta, sigma=sigma, B1=B1, B2=B2, 
    nu=nu, c=c, z=z, d=d, xi=xi, tautilde=tautilde, pi2=pi2, pi1=pi1, thets=thets, theti=theti, ut=ut);

    # Empty vectors for simulation results
    I0 = 0.02;
    S0 = 1-I0/(1-Par.tautilde);
    R0  = 0;
    IT = zeros(TT);
    ST = zeros(TT);
    RT = zeros(TT);
    qT = zeros(TT);
    lamT = zeros(TT);
    muT = zeros(TT);
    #Del = zeros(TT);
    VT = zeros(TT);
    YT = zeros(TT);
    L1 = zeros(TT);
    L2 = zeros(TT);
    Rknot = zeros(TT);
    Test = zeros(TT);
    IN = zeros(TT);
    C1T = zeros(TT);
    C2T = zeros(TT);
    IT[1] = I0;
    ST[1] = S0;
    RT[1] = R0;
    VV = 0.0;
    VVC = zeros(TT);
    masstest = zeros(TT);
    mctest = zeros(TT);

    # Simulate optimal path
    for i in 1:TT
        S = ST[i];
        I = IT[i];
        R = RT[i];
        IN[i] = 1-S-R;
        qT[i] = 0;
        lamT[i] = .9901;
        muT[i] = .9901;
        Del = 0;
        LI1 = lamT[i].*(S+(1-Del).*I)+muT[i].*(1-Par.delta).*R;
        LI2 = (1-lamT[i]).*(S+(1-Del).*I)+(1-muT[i]).*(1-Par.delta).*R;
        if LI1>0
            Ipr1 = (lamT[i].^2 .* Par.pi1.*(1-Del).*I)./LI1;
            RIpr1 = (lamT[i].^2 .* Par.pi1)./LI1;
        else
            Ipr1 = 0;
            RIpr1 = 0;
        end
        if LI2>0
            Ipr2 = ((1-lamT[i]).^2 .* Par.pi2.*(1-Del).*I)./LI2;
            RIpr2 = ((1-lamT[i]).^2 .* Par.pi2)./LI2;
        else
            Ipr2 = 0;
            RIpr2 = 0;
        end
        Ipr = Ipr1 + Ipr2;
        Ipr = min(Ipr,1);
        if i<TT
            IT[i+1] = min(max((1-Del).*(1-Par.tautilde).*(1-Par.gamma).*I+S.*Ipr,0),1);
            ST[i+1] = min(max((1-Ipr).*S,0),1);
            RT[i+1] = max(min(R + Par.gamma.*IN[i],1),0);
        end
        
        YT[i] = YY(Par,S,I,R,qT[i],lamT[i],muT[i]);
        C1T[i],C2T[i] = CC(Par,S,I,R,qT[i],lamT[i],muT[i]);
        L1[i] = (lamT[i].*(S+Par.xi.*(1-Del).*I)+muT[i].*(1-Par.delta).*R);
        L2[i] = ((1-lamT[i]).*(S+Par.xi.*(1-Del).*I)+(1-muT[i]).*(1-Par.delta).*R);
        
        if I>0
            Rknot[i]=(1 ./ (Par.gamma)) .* Ipr ./ (1-S-R);
        else
            Rknot[i] = 0;
        end
        VV += Par.beta.^(i-1) * (log(YT[i]-C1T[i])-C2T[i]);
        VVC[i] = VV;
    end

    # Calculate output variables
    VV = VV + Par.beta^TT/(1-Par.beta)*VTerm(Par,ST[TT],IT[TT],RT[TT]);
    cumdeaths = cumsum(Par.gamma.*Par.delta.*max.(1 .- ST .- RT,0));
    Totaldeath =cumdeaths[TT];
    Ycumsum = YT[1]-YY_noep(Par,1);
    for i = 2:TT
        Ycumsum = Ycumsum + Par.beta.^(i-1) .* (YT[i]-YY_noep(Par,1));
    end
    Ycumsum = (1-Par.beta).*Ycumsum;
    compvar2 = 1/YY_noep(Par,1)*exp((1-Par.beta)*VV);
    compvar2 = 1-compvar2;
    #suscepss = ST[TT];
    #Iss = 1-ST[TT]-RT[TT];
    #maxI = maximum(1 .- ST .- RT);
    #R01 = Rknot[1];
    #R026 = Rknot[26];
    #R052 = Rknot[52];

    return Dict("Par" => Par, "TT" => TT, "ST" => ST, "IT" => IT, "RT" => RT, 
    "YT" => YT, "C1T" => C1T, "Rknot" => Rknot, "masstest" => masstest, "mctest" => mctest)

end
    
