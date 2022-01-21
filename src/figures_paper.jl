############################
#                          #
# Generate Figures         #
#                          #
############################

### Figure 5
# Using user-provided results
"""
    createFig5(noint, notest, untarget)
    createFig5()

Replicate Figure 5 in Chari, Kirpalani, and Phelan (2021). 
If no arguments are provided, then stored results obtained with default parameters are used.

# Arguments:
- `noint::Dict`: Output from [`nopolicy`](@ref).
- `notest::Dict`: Output from [`withpolicy(::String)`](@ref) with policy `notest`.
- `untarget::Dict`: Output from [`withpolicy(::String)`](@ref) with policy `untarget`.

"""
function createFig5(noint::Dict, notest::Dict, untarget::Dict) 

    println("Replicating Figure 5...")

    # Create variables
    TT = noint["TT"];
    IN_1 = 1 .- noint["ST"] .- noint["RT"];
    IN_2 = 1 .- notest["ST"] .- notest["RT"];
    IN_3 = 1 .- untarget["ST"] .- untarget["RT"];
    death_1 = cumsum(noint["Par"].gamma .* noint["Par"].delta .* IN_1);
    death_2 = cumsum(notest["Par"].gamma .* notest["Par"].delta .* IN_2);
    death_3 = cumsum(untarget["Par"].gamma .* untarget["Par"].delta .* IN_3);
    Yss = 1.01;
    Cdev_1 = ((noint["YT"]-noint["C1T"]).-Yss)./Yss.*100;
    Cdev_2 = ((notest["YT"]-notest["C1T"]).-Yss)./Yss.*100;
    Cdev_3 = ((untarget["YT"]-untarget["C1T"]).-Yss)./Yss.*100;

    # Figure 5
    args_first = (titlefontsize=8, guidefontsize=6, linecolor=["black" "green" "red"], legendfontsize=6, label=["no policy" "no testing" "untargeted"]);
    args_rest  = (titlefontsize=8, guidefontsize=6, linecolor=["black" "green" "red"], legend=false);
    p1 = plot(1:TT, [IN_1.*100 IN_2.*100 IN_3.*100], title="Infected, I_t", ylabel="%"; args_first...)
    p2 = plot(1:TT, [noint["ST"].*100 notest["ST"].*100 untarget["ST"].*100], title="Susceptible, S_t", ylabel="%"; args_rest...)
    p3 = plot(1:TT, [death_1.*100 death_2.*100 death_3.*100], title="Cumulative Deaths", ylabel="%"; args_rest...)
    p4 = plot(1:TT, [noint["Rknot"] notest["Rknot"] untarget["Rknot"]], title="R_0"; args_rest...)
    p5 = plot(1:TT, [noint["IT"].*100 notest["IT"].*100 untarget["IT"].*100], title="Unknown Infected, ~I_t", ylabel="%"; args_rest...)
    p6 = plot(1:TT, [Cdev_1 Cdev_2 Cdev_3], title="Agg. Consumption, C_t", ylabel="% dev. from steady state"; args_rest...)
    p7 = plot(1:TT, [noint["masstest"].*100 notest["masstest"].*100 untarget["masstest"].*100], title="Mass tested", ylabel="%"; args_rest...)
    p8 = plot(1:TT, [noint["mctest"] notest["mctest"] untarget["mctest"]], title="Marginal testing cost", ylabel="Dollars"; args_rest...)
    p = plot(p1, p2, p3, p4, p5, p6, p7, p8, layout = (2, 4));
    display(p)

end


# Using stored results
function createFig5() 

    println("Replicating Figure 5...")

    # Load stored data
    noint = load("src/data/calib1_nointervention.jld2");
    notest = load("src/data/calib1_notest_small.jld2");
    untarget = load("src/data/calib1_untarget_small.jld2");
    
    # Create variables
    TT = noint["TT"];
    IN_1 = 1 .- noint["ST"] .- noint["RT"];
    IN_2 = 1 .- notest["ST"] .- notest["RT"];
    IN_3 = 1 .- untarget["ST"] .- untarget["RT"];
    death_1 = cumsum(noint["Par"].gamma .* noint["Par"].delta .* IN_1);
    death_2 = cumsum(notest["Par"].gamma .* notest["Par"].delta .* IN_2);
    death_3 = cumsum(untarget["Par"].gamma .* untarget["Par"].delta .* IN_3);
    Yss = 1.01;
    Cdev_1 = ((noint["YT"]-noint["C1T"]).-Yss)./Yss.*100;
    Cdev_2 = ((notest["YT"]-notest["C1T"]).-Yss)./Yss.*100;
    Cdev_3 = ((untarget["YT"]-untarget["C1T"]).-Yss)./Yss.*100;

    # Figure 5
    args_first = (titlefontsize=8, guidefontsize=6, linecolor=["black" "green" "red"], legendfontsize=6, label=["no policy" "no testing" "untargeted"]);
    args_rest  = (titlefontsize=8, guidefontsize=6, linecolor=["black" "green" "red"], legend=false);
    p1 = plot(1:TT, [IN_1.*100 IN_2.*100 IN_3.*100], title="Infected, I_t", ylabel="%"; args_first...)
    p2 = plot(1:TT, [noint["ST"].*100 notest["ST"].*100 untarget["ST"].*100], title="Susceptible, S_t", ylabel="%"; args_rest...)
    p3 = plot(1:TT, [death_1.*100 death_2.*100 death_3.*100], title="Cumulative Deaths", ylabel="%"; args_rest...)
    p4 = plot(1:TT, [noint["Rknot"] notest["Rknot"] untarget["Rknot"]], title="R_0"; args_rest...)
    p5 = plot(1:TT, [noint["IT"].*100 notest["IT"].*100 untarget["IT"].*100], title="Unknown Infected, ~I_t", ylabel="%"; args_rest...)
    p6 = plot(1:TT, [Cdev_1 Cdev_2 Cdev_3], title="Agg. Consumption, C_t", ylabel="% dev. from steady state"; args_rest...)
    p7 = plot(1:TT, [noint["masstest"].*100 notest["masstest"].*100 untarget["masstest"].*100], title="Mass tested", ylabel="%"; args_rest...)
    p8 = plot(1:TT, [noint["mctest"] notest["mctest"] untarget["mctest"]], title="Marginal testing cost", ylabel="Dollars"; args_rest...)
    p = plot(p1, p2, p3, p4, p5, p6, p7, p8, layout = (2, 4));
    display(p)

end



##########################################################################################
### Figure 6
# Using user-provided results
"""
    createFig6(notest, untarget, target)
    createFig6()

Analogous to [`createFig5`](@ref).
"""
function createFig6(notest::Dict, untarget::Dict, target::Dict) 

    println("Replicating Figure 6...")

    # Create variables
    TT = notest["TT"];
    IN_2 = 1 .- notest["ST"] .- notest["RT"];
    IN_3 = 1 .- untarget["ST"] .- untarget["RT"];
    IN_4 = 1 .- target["ST"] .- target["RT"];
    death_2 = cumsum(notest["Par"].gamma .* notest["Par"].delta .* IN_2);
    death_3 = cumsum(untarget["Par"].gamma .* untarget["Par"].delta .* IN_3);
    death_4 = cumsum(target["Par"].gamma .* target["Par"].delta .* IN_4);
    Yss = 1.01;
    Cdev_2 = ((notest["YT"]-notest["C1T"]).-Yss)./Yss.*100;
    Cdev_3 = ((untarget["YT"]-untarget["C1T"]).-Yss)./Yss.*100;
    Cdev_4 = ((target["YT"]-target["C1T"]).-Yss)./Yss.*100;

    # Figure 6
    args_first = (titlefontsize=8, guidefontsize=6, linecolor=["green" "red" "blue"], legendfontsize=6, label=["no testing" "untargeted" "targeted"]);
    args_rest  = (titlefontsize=8, guidefontsize=6, linecolor=["green" "red" "blue"], legend=false);
    p1 = plot(1:TT, [IN_2.*100 IN_3.*100 IN_4.*100], title="Infected, I_t", ylabel="%"; args_first...)
    p2 = plot(1:TT, [notest["ST"].*100 untarget["ST"].*100 target["ST"].*100], title="Susceptible, S_t", ylabel="%"; args_rest...)
    p3 = plot(1:TT, [death_2.*100 death_3.*100 death_4.*100], title="Cumulative Deaths", ylabel="%"; args_rest...)
    p4 = plot(1:TT, [notest["Rknot"] untarget["Rknot"] target["Rknot"]], title="R_0"; args_rest...)
    p5 = plot(1:TT, [notest["IT"].*100 untarget["IT"].*100 target["IT"].*100], title="Unknown Infected, ~I_t", ylabel="%"; args_rest...)
    p6 = plot(1:TT, [Cdev_2 Cdev_3 Cdev_4], title="Agg. Consumption, C_t", ylabel="% dev. from steady state"; args_rest...)
    p7 = plot(1:TT, [notest["masstest"].*100 untarget["masstest"].*100 target["masstest"].*100], title="Mass tested", ylabel="%"; args_rest...)
    p8 = plot(1:TT, [notest["mctest"] untarget["mctest"] target["mctest"]], title="Marginal testing cost", ylabel="Dollars"; args_rest...)
    p = plot(p1, p2, p3, p4, p5, p6, p7, p8, layout = (2, 4));
    display(p)

end



# Using stored results
function createFig6() 

    println("Replicating Figure 6...")

    # Load stored data
    notest = load("src/data/calib1_notest_small.jld2");
    untarget = load("src/data/calib1_untarget_small.jld2");
    target = load("src/data/calib1_target_small.jld2");
    
    # Create variables
    TT = notest["TT"];
    IN_2 = 1 .- notest["ST"] .- notest["RT"];
    IN_3 = 1 .- untarget["ST"] .- untarget["RT"];
    IN_4 = 1 .- target["ST"] .- target["RT"];
    death_2 = cumsum(notest["Par"].gamma .* notest["Par"].delta .* IN_2);
    death_3 = cumsum(untarget["Par"].gamma .* untarget["Par"].delta .* IN_3);
    death_4 = cumsum(target["Par"].gamma .* target["Par"].delta .* IN_4);
    Yss = 1.01;
    Cdev_2 = ((notest["YT"]-notest["C1T"]).-Yss)./Yss.*100;
    Cdev_3 = ((untarget["YT"]-untarget["C1T"]).-Yss)./Yss.*100;
    Cdev_4 = ((target["YT"]-target["C1T"]).-Yss)./Yss.*100;

    # Figure 6
    args_first = (titlefontsize=8, guidefontsize=6, linecolor=["green" "red" "blue"], legendfontsize=6, label=["no testing" "untargeted" "targeted"]);
    args_rest  = (titlefontsize=8, guidefontsize=6, linecolor=["green" "red" "blue"], legend=false);
    p1 = plot(1:TT, [IN_2.*100 IN_3.*100 IN_4.*100], title="Infected, I_t", ylabel="%"; args_first...)
    p2 = plot(1:TT, [notest["ST"].*100 untarget["ST"].*100 target["ST"].*100], title="Susceptible, S_t", ylabel="%"; args_rest...)
    p3 = plot(1:TT, [death_2.*100 death_3.*100 death_4.*100], title="Cumulative Deaths", ylabel="%"; args_rest...)
    p4 = plot(1:TT, [notest["Rknot"] untarget["Rknot"] target["Rknot"]], title="R_0"; args_rest...)
    p5 = plot(1:TT, [notest["IT"].*100 untarget["IT"].*100 target["IT"].*100], title="Unknown Infected, ~I_t", ylabel="%"; args_rest...)
    p6 = plot(1:TT, [Cdev_2 Cdev_3 Cdev_4], title="Agg. Consumption, C_t", ylabel="% dev. from steady state"; args_rest...)
    p7 = plot(1:TT, [notest["masstest"].*100 untarget["masstest"].*100 target["masstest"].*100], title="Mass tested", ylabel="%"; args_rest...)
    p8 = plot(1:TT, [notest["mctest"] untarget["mctest"] target["mctest"]], title="Marginal testing cost", ylabel="Dollars"; args_rest...)
    p = plot(p1, p2, p3, p4, p5, p6, p7, p8, layout = (2, 4));
    display(p)

end


##########################################################################################
### Figure 7
# Using user-provided results
"""
    createFig7(untarget, target, isolate)
    createFig7()

Analogous to [`createFig5`](@ref).
"""
function createFig7(untarget::Dict, target::Dict, isolate::Dict) 

    println("Replicating Figure 7...")

    # Create variables
    TT = untarget["TT"];
    IN_3 = 1 .- untarget["ST"] .- untarget["RT"];
    IN_4 = 1 .- target["ST"] .- target["RT"];
    IN_5 = 1 .- isolate["ST"] .- isolate["RT"];
    death_3 = cumsum(untarget["Par"].gamma .* untarget["Par"].delta .* IN_3);
    death_4 = cumsum(target["Par"].gamma .* target["Par"].delta .* IN_4);
    death_5 = cumsum(isolate["Par"].gamma .* isolate["Par"].delta .* IN_5);
    Yss = 1.01;
    Cdev_3 = ((untarget["YT"]-untarget["C1T"]).-Yss)./Yss.*100;
    Cdev_4 = ((target["YT"]-target["C1T"]).-Yss)./Yss.*100;
    Cdev_5 = ((isolate["YT"]-isolate["C1T"]).-Yss)./Yss.*100;

    # Figure 7
    args_first = (titlefontsize=8, guidefontsize=6, linecolor=["red" "blue" "purple"], legendfontsize=6, label=["untargeted" "targeted" "isolate"]);
    args_rest  = (titlefontsize=8, guidefontsize=6, linecolor=["red" "blue" "purple"], legend=false);
    p1 = plot(1:TT, [IN_3.*100 IN_4.*100 IN_5.*100], title="Infected, I_t", ylabel="%"; args_first...)
    p2 = plot(1:TT, [untarget["ST"].*100 target["ST"].*100 isolate["ST"].*100], title="Susceptible, S_t", ylabel="%"; args_rest...)
    p3 = plot(1:TT, [death_3.*100 death_4.*100 death_5.*100], title="Cumulative Deaths", ylabel="%"; args_rest...)
    p4 = plot(1:TT, [untarget["Rknot"] target["Rknot"] isolate["Rknot"]], title="R_0"; args_rest...)
    p5 = plot(1:TT, [untarget["IT"].*100 target["IT"].*100 isolate["IT"].*100], title="Unknown Infected, ~I_t", ylabel="%"; args_rest...)
    p6 = plot(1:TT, [Cdev_3 Cdev_4 Cdev_5], title="Agg. Consumption, C_t", ylabel="% dev. from steady state"; args_rest...)
    p7 = plot(1:TT, [untarget["masstest"].*100 target["masstest"].*100 isolate["masstest"].*100], title="Mass tested", ylabel="%"; args_rest...)
    p8 = plot(1:TT, [untarget["mctest"] target["mctest"] isolate["mctest"] ], title="Marginal testing cost", ylabel="Dollars"; args_rest...)
    p = plot(p1, p2, p3, p4, p5, p6, p7, p8, layout = (2, 4));
    display(p)

end



# Using stored results
function createFig7() 

    println("Replicating Figure 7...")

    # Load stored data
    untarget = load("src/data/calib1_untarget_small.jld2");
    target = load("src/data/calib1_target_small.jld2");
    isolate = load("src/data/calib1_isolate_small.jld2");

    # Create variables
    TT = untarget["TT"];
    IN_3 = 1 .- untarget["ST"] .- untarget["RT"];
    IN_4 = 1 .- target["ST"] .- target["RT"];
    IN_5 = 1 .- isolate["ST"] .- isolate["RT"];
    death_3 = cumsum(untarget["Par"].gamma .* untarget["Par"].delta .* IN_3);
    death_4 = cumsum(target["Par"].gamma .* target["Par"].delta .* IN_4);
    death_5 = cumsum(isolate["Par"].gamma .* isolate["Par"].delta .* IN_5);
    Yss = 1.01;
    Cdev_3 = ((untarget["YT"]-untarget["C1T"]).-Yss)./Yss.*100;
    Cdev_4 = ((target["YT"]-target["C1T"]).-Yss)./Yss.*100;
    Cdev_5 = ((isolate["YT"]-isolate["C1T"]).-Yss)./Yss.*100;

    # Figure 7
    args_first = (titlefontsize=8, guidefontsize=6, linecolor=["red" "blue" "purple"], legendfontsize=6, label=["untargeted" "targeted" "isolate"]);
    args_rest  = (titlefontsize=8, guidefontsize=6, linecolor=["red" "blue" "purple"], legend=false);
    p1 = plot(1:TT, [IN_3.*100 IN_4.*100 IN_5.*100], title="Infected, I_t", ylabel="%"; args_first...)
    p2 = plot(1:TT, [untarget["ST"].*100 target["ST"].*100 isolate["ST"].*100], title="Susceptible, S_t", ylabel="%"; args_rest...)
    p3 = plot(1:TT, [death_3.*100 death_4.*100 death_5.*100], title="Cumulative Deaths", ylabel="%"; args_rest...)
    p4 = plot(1:TT, [untarget["Rknot"] target["Rknot"] isolate["Rknot"]], title="R_0"; args_rest...)
    p5 = plot(1:TT, [untarget["IT"].*100 target["IT"].*100 isolate["IT"].*100], title="Unknown Infected, ~I_t", ylabel="%"; args_rest...)
    p6 = plot(1:TT, [Cdev_3 Cdev_4 Cdev_5], title="Agg. Consumption, C_t", ylabel="% dev. from steady state"; args_rest...)
    p7 = plot(1:TT, [untarget["masstest"].*100 target["masstest"].*100 isolate["masstest"].*100], title="Mass tested", ylabel="%"; args_rest...)
    p8 = plot(1:TT, [untarget["mctest"] target["mctest"] isolate["mctest"] ], title="Marginal testing cost", ylabel="Dollars"; args_rest...)
    p = plot(p1, p2, p3, p4, p5, p6, p7, p8, layout = (2, 4));
    display(p)

end

