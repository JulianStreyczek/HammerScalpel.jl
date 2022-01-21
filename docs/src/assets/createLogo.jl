##########################################################################################
### Logo
using JLD2, Plots
function createLogo()

    # Load stored data
    untarget = load("./src/data/calib1_untarget_small.jld2");
    target = load("./src/data/calib1_target_small.jld2");
    isolate = load("./src/data/calib1_isolate_small.jld2");

    # Create variables
    TT = untarget["TT"];
    IN_3 = 1 .- untarget["ST"] .- untarget["RT"];
    IN_4 = 1 .- target["ST"] .- target["RT"];
    IN_5 = 1 .- isolate["ST"] .- isolate["RT"];
    Yss = 1.01;
    Cdev_3 = ((untarget["YT"]-untarget["C1T"]).-Yss)./Yss.*100;
    Cdev_4 = ((target["YT"]-target["C1T"]).-Yss)./Yss.*100;
    Cdev_5 = ((isolate["YT"]-isolate["C1T"]).-Yss)./Yss.*100;

    # Figure 7
    args_first = (titlefontsize=16, guidefontsize=10, linecolor=["red" "blue" "purple"], legendfontsize=12, label=["untargeted" "targeted" "isolate"]);
    args_rest  = (titlefontsize=16, guidefontsize=10, linecolor=["red" "blue" "purple"], legend=false);
    p1 = plot(1:TT, [IN_3.*100 IN_4.*100 IN_5.*100], title="Infected, I_t", ylims=(-2,30); args_first...)
    plot!(size=(800,400))
    p6 = plot(1:TT, [Cdev_3 Cdev_4 Cdev_5], title="Agg. Consumption, C_t", ylabel="% dev. from steady state"; args_rest...)
    plot!(size=(800,400))
    plot(p1, p6, layout=(1,2))
    savefig("docs/src/assets/logo.png")

end

#createLogo()

function createHeader()
    # Load stored data
    untarget = load("./src/data/calib1_untarget_small.jld2");
    target = load("./src/data/calib1_target_small.jld2");
    isolate = load("./src/data/calib1_isolate_small.jld2");

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
    args_first = (titlefontsize=8, guidefontsize=6, linecolor=["red" "blue" "purple"], labelfontsize=1, label=["untargeted" "targeted" "isolate"]);
    args_rest  = (titlefontsize=8, guidefontsize=6, linecolor=["red" "blue" "purple"], legend=false);
    p1 = plot(1:TT, [IN_3.*100 IN_4.*100 IN_5.*100], title="Infected, I_t", ylabel=""; args_first...)
    plot!(size=(800,200))
    p2 = plot(1:TT, [untarget["ST"].*100 target["ST"].*100 isolate["ST"].*100], title="Susceptible, S_t", ylabel="%"; args_rest...)
    plot!(size=(800,200))
    p3 = plot(1:TT, [death_3.*100 death_4.*100 death_5.*100], title="Cumulative Deaths", ylabel="%"; args_rest...)
    plot!(size=(800,200))
    p6 = plot(1:TT, [Cdev_3 Cdev_4 Cdev_5], title="Agg. Consumption, C_t", ylabel="% dev. from steady state"; args_rest...)
    plot!(size=(800,200))
    p = plot(p1, p2, p3, p6, layout = (1, 4));
    savefig("docs/src/assets/header.png")

end

createHeader()