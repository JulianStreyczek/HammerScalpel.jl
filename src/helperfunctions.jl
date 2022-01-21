###############################################
#                                             #
# Helper functions:                           #
# Not of immediate interest, but useful       #
# for solving the model                       #
#                                             #
###############################################

# Included:
# - Bellman
# - Bellman_NLopt
# - CC
# - VTerm
# - YY_noep
# - YY



##########################################################################################
# Computes (negative of) Bellman equation
function Bellman(Par,Grid,S,I,R,V0,x)
    q = x[1];
    lam = x[2];
    mu = x[3];
    IN = max(1-S-R,0);
    if Par.TEST == 1
        Del = q.*Par.theti;
    else
        Del = 0;
    end
    if Par.ISO == 1
        LI1 = lam.*((1-q.*Par.thets).*S+(1-q.*Par.theti).*I)+mu.*(1-Par.delta).*R;
        LI2 = (1-lam).*((1-q.*Par.thets).*S+(1-q.*Par.theti).*I)+(1-mu).*(1-Par.delta).*R;
        if LI1>0
            Ipr1 = (lam.^2 .* Par.pi1 .* (1-q.*Par.theti).*I)./LI1;
        else
            Ipr1 = 0;
        end
        if LI2>0
            Ipr2 = ((1-lam).^2 .* Par.pi2 .* (1-q.*Par.theti).*I)./LI2;
        else
            Ipr2 = 0;
        end
        Ipr = Ipr1 + Ipr2;
        Ipr = min(Ipr,1);
        Ip = min(max((1-Par.tautilde).*(1-Par.gamma).*I+(1-q.*Par.thets).*S.*Ipr,0),1);
        Sp = min(max((1-(1-q.*Par.thets).*Ipr).*S,0),1);
        Rp = max(min(R + Par.gamma.*IN,1),0);
    else
        LI1 = lam.*(S+(1-Del).*I)+mu.*(1-Par.delta).*R;
        LI2 = (1-lam).*(S+(1-Del).*I)+(1-mu).*(1-Par.delta).*R;
        if LI1>0
            Ipr1 = (lam.^2 .* Par.pi1 .* (1-Del).*I)./LI1;
        else
            Ipr1 = 0;
        end
        if LI2>0
            Ipr2 = ((1-lam).^2 .* Par.pi2 .* (1-Del).*I)./LI2;
        else
            Ipr2 = 0;
        end
        Ipr = Ipr1 + Ipr2;
        Ipr = min(Ipr,1);
        Ip = min(max((1-Del).*(1-Par.tautilde).*(1-Par.gamma).*I+S.*Ipr,0),1);
        Sp = min(max((1-Ipr).*S,0),1);
        Rp = max(min(R + Par.gamma.*IN,1),0);
    end
    
    Y = YY(Par,S,I,R,q,lam,mu);
    C1,C2 = CC(Par,S,I,R,q,lam,mu);
    OO = max(Y-C1,0);
    
    if OO>0
        #grid = SimplexGrid(Grid.S, Grid.I, Grid.R);
        grid = RectangleGrid(Grid.S, Grid.I, Grid.R);
        gridData = V0;
        V_interpolated = ((1-Par.beta).*(log(OO) - C2) + Par.beta .* interpolate(grid, gridData, [Sp, Ip, Rp]));
        if isnan(V_interpolated)
           V = -VTerm(Par,0,0,0);
        else
            V = -min(VTerm(Par,0,0,0), V_interpolated);
        end
    else
        V = -(-100000000);
    end

    return V
end



##########################################################################################
# Wrapper for Bellman, as expected by NLopt optimizer
function Bellman_NLopt(Par, Grid, S, I, R, V00, x::Vector, g::Vector)
    gradient = grad(backward_fdm(2, 1), z -> Bellman(Par, Grid, S, I, R, V00, z), x)[1]
    if length(g) > 0
        g[1] = gradient[1]
        g[2] = gradient[2]
        g[3] = gradient[3]
    end
    return Bellman(Par, Grid, S, I, R, V00, x)
end
    
    

##########################################################################################
# Computes total cost (infection+testing+death)
function CC(Par,S,I,R,q,lam,mu)
    IN = max(1-S-R,0);
    if Par.TEST == 1
        x1 = Par.c.*(q.*min((Par.thets.*S+Par.theti.*I),1)).^(1+Par.nu)./(1+Par.nu);
        x2 = Par.z.*(IN-I)+Par.d.*Par.gamma.*Par.delta.*IN;
    else
        x1 = 0;
        x2 = Par.z.*(IN-I)+Par.d.*Par.gamma.*Par.delta.*IN;
    end
    return [x1, x2]  
end


##########################################################################################
# Computes terminal value function
function VTerm(Par,S,I,R)
    lam = ((Par.B1./Par.B2).^Par.sigma)./(1+(Par.B1./Par.B2).^Par.sigma);
    y = (Par.B1.*lam.^((Par.sigma-1)./(Par.sigma))+Par.B2.*(1-lam).^((Par.sigma-1)./(Par.sigma))).^(Par.sigma./(Par.sigma-1));
    return log(y.*(1-Par.delta.*R));
    end
    


##########################################################################################
# Computes output in the no epidemic steady state:
# Using Par:
function YY_noep(Par, L)
    lam = ((Par.B1./Par.B2).^Par.sigma)./(1+(Par.B1./Par.B2).^Par.sigma);
    y = L.*(Par.B1.*lam.^((Par.sigma-1)./(Par.sigma))+Par.B2.*(1-lam).^((Par.sigma-1)./(Par.sigma))).^(Par.sigma./(Par.sigma-1));
    return y
end
# Using parameters of interest directly
function YY_noep(B1, B2, sigma, L)
    lam = ((B1./B2).^sigma) ./ (1+(B1./B2).^sigma);
    y = L.*(B1.*lam.^((sigma-1)./(sigma))+B2.*(1-lam).^((sigma-1)./(sigma))).^(sigma./(sigma-1));
    return y
end



##########################################################################################
# Computes output
function YY(Par,S,I,R,q,lam,mu)
    if Par.ISO == 0 && Par.TEST == 0
        L1 = (lam.*(S+Par.xi.*I)+mu.*(1-Par.delta).*R);
        L2 = ((1-lam).*(S+Par.xi.*I)+(1-mu).*(1-Par.delta).*R);
        y = (Par.B1.*L1.^((Par.sigma-1)/Par.sigma)+Par.B2.*L2.^((Par.sigma-1)/Par.sigma)).^(Par.sigma./((Par.sigma-1)));
    elseif Par.ISO == 1
        L1 = (lam.*(S.*(1-q.*Par.thets)+Par.xi.*(1-q.*Par.theti).*I)+mu.*(1-Par.delta).*R);
        L2 = ((1-lam).*(S.*(1-q.*Par.thets)+Par.xi.*(1-q.*Par.theti).*I)+(1-mu).*(1-Par.delta).*R);
        y = (Par.B1.*L1.^((Par.sigma-1)/Par.sigma)+Par.B2.*L2.^((Par.sigma-1)/Par.sigma)).^(Par.sigma./((Par.sigma-1)));
    else
        L1 = (lam.*(S+Par.xi.*(1-q.*Par.theti).*I)+mu.*(1-Par.delta).*R);
        L2 = ((1-lam).*(S+Par.xi.*(1-q.*Par.theti).*I)+(1-mu).*(1-Par.delta).*R);
        y = (Par.B1.*L1.^((Par.sigma-1)/Par.sigma)+Par.B2.*L2.^((Par.sigma-1)/Par.sigma)).^(Par.sigma./((Par.sigma-1)));
    end
    return y
end