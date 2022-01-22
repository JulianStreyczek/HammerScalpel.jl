var documenterSearchIndex = {"docs":
[{"location":"","page":"Replication of Chari, Kirpalani, and Phelan (2021)","title":"Replication of Chari, Kirpalani, and Phelan (2021)","text":"(Image: header)","category":"page"},{"location":"#Replication-of-Chari,-Kirpalani,-and-Phelan-(2021)","page":"Replication of Chari, Kirpalani, and Phelan (2021)","title":"Replication of Chari, Kirpalani, and Phelan (2021)","text":"","category":"section"},{"location":"","page":"Replication of Chari, Kirpalani, and Phelan (2021)","title":"Replication of Chari, Kirpalani, and Phelan (2021)","text":"This project was conducted as a final assignment for the PhD course Numerical Methods at Bocconi University in Fall 2021.","category":"page"},{"location":"","page":"Replication of Chari, Kirpalani, and Phelan (2021)","title":"Replication of Chari, Kirpalani, and Phelan (2021)","text":"This package replicates the main figures of Chari, Kirpalani, and Phelan (2021): \"The hammer and the scalpel: On the economics of indiscriminate versus targeted isolation policies during pandemics\", published in the Review of Economic Dynamics. The paper is available here, and the original replication material in Matlab here.","category":"page"},{"location":"","page":"Replication of Chari, Kirpalani, and Phelan (2021)","title":"Replication of Chari, Kirpalani, and Phelan (2021)","text":"The authors develop a theoretical model that combines epidemic transmission and economic outcomes. They study the effect of different government responses, such as quarantining, testing, contact-tracing, and isolation. ","category":"page"},{"location":"#Installation","page":"Replication of Chari, Kirpalani, and Phelan (2021)","title":"Installation","text":"","category":"section"},{"location":"","page":"Replication of Chari, Kirpalani, and Phelan (2021)","title":"Replication of Chari, Kirpalani, and Phelan (2021)","text":"This package is not listed as an official Julia package. If you wish to use it, you may clone this Github repository to your machine. You then have two options:","category":"page"},{"location":"","page":"Replication of Chari, Kirpalani, and Phelan (2021)","title":"Replication of Chari, Kirpalani, and Phelan (2021)","text":"Use as a package: Start up Julia and go into the package editor by typing ]. Then, type","category":"page"},{"location":"","page":"Replication of Chari, Kirpalani, and Phelan (2021)","title":"Replication of Chari, Kirpalani, and Phelan (2021)","text":"activate .\ninstantiate","category":"page"},{"location":"","page":"Replication of Chari, Kirpalani, and Phelan (2021)","title":"Replication of Chari, Kirpalani, and Phelan (2021)","text":"After pressing backspace, go back into Julia's standard command mode and type","category":"page"},{"location":"","page":"Replication of Chari, Kirpalani, and Phelan (2021)","title":"Replication of Chari, Kirpalani, and Phelan (2021)","text":"using HammerScalpel","category":"page"},{"location":"","page":"Replication of Chari, Kirpalani, and Phelan (2021)","title":"Replication of Chari, Kirpalani, and Phelan (2021)","text":"The functions should now be ready to go.","category":"page"},{"location":"","page":"Replication of Chari, Kirpalani, and Phelan (2021)","title":"Replication of Chari, Kirpalani, and Phelan (2021)","text":"Use as ordinary code: Play around with the files main.jl or main_parallel.jl, which run the code in sequential and parallelized form, respectively. ","category":"page"},{"location":"#Using-the-package","page":"Replication of Chari, Kirpalani, and Phelan (2021)","title":"Using the package","text":"","category":"section"},{"location":"","page":"Replication of Chari, Kirpalani, and Phelan (2021)","title":"Replication of Chari, Kirpalani, and Phelan (2021)","text":"If you wish run the replication from beginning to end, simply type","category":"page"},{"location":"","page":"Replication of Chari, Kirpalani, and Phelan (2021)","title":"Replication of Chari, Kirpalani, and Phelan (2021)","text":"solveModel()","category":"page"},{"location":"","page":"Replication of Chari, Kirpalani, and Phelan (2021)","title":"Replication of Chari, Kirpalani, and Phelan (2021)","text":"This will call all necessary functions in the correct order.  Specifically, it first computes the no-intervention outcome using nopolicy, and the respective outcomes under the various policies using withpolicy. Then, it plots the evolution of the model under different policies over time using createFig5, createFig6, and createFig7.","category":"page"},{"location":"","page":"Replication of Chari, Kirpalani, and Phelan (2021)","title":"Replication of Chari, Kirpalani, and Phelan (2021)","text":"Since solving the model may take a while, you have the option to run them more efficiently on multiple CPU cores at the same time.  To do so, type","category":"page"},{"location":"","page":"Replication of Chari, Kirpalani, and Phelan (2021)","title":"Replication of Chari, Kirpalani, and Phelan (2021)","text":"solveModel(nprocs)","category":"page"},{"location":"","page":"Replication of Chari, Kirpalani, and Phelan (2021)","title":"Replication of Chari, Kirpalani, and Phelan (2021)","text":"where nprocs specifies the number of logical cores on your CPU. ","category":"page"},{"location":"","page":"Replication of Chari, Kirpalani, and Phelan (2021)","title":"Replication of Chari, Kirpalani, and Phelan (2021)","text":"If you just want to check my results without running the simulations yourself, you may use stored outcomes that I obtained under default parameters. You can recreate the figures by running","category":"page"},{"location":"","page":"Replication of Chari, Kirpalani, and Phelan (2021)","title":"Replication of Chari, Kirpalani, and Phelan (2021)","text":"createFig5()\ncreateFig6()\ncreateFig7()","category":"page"},{"location":"#Documentation","page":"Replication of Chari, Kirpalani, and Phelan (2021)","title":"Documentation","text":"","category":"section"},{"location":"","page":"Replication of Chari, Kirpalani, and Phelan (2021)","title":"Replication of Chari, Kirpalani, and Phelan (2021)","text":"solveModel\nnopolicy\nwithpolicy(::String)\ncreateFig5(::Dict, ::Dict, ::Dict)\ncreateFig6(::Dict, ::Dict, ::Dict)\ncreateFig7(::Dict, ::Dict, ::Dict)","category":"page"},{"location":"#HammerScalpel.solveModel","page":"Replication of Chari, Kirpalani, and Phelan (2021)","title":"HammerScalpel.solveModel","text":"solveModel(; TT=52, ndims=20, thetai=0.38, thetas=0.0044)\nsolveModel(nprocs:Int64; TT=52, ndims=20, thetai=0.38, thetas=0.0044)\n\nSolve the model for all possible policy choices and replicate Figures 5, 6, and 7 in Chari, Kirpalani, and Phelan (2021).  The function calls nopolicy and withpolicy with the relevant keywords, and uses the outputs to generate the figures using  createFig5, createFig6, and createFig7.    \n\nArguments:\n\nnprocs::Int64: Optional. If specified, uses the module Distributed for parallelizing computations. \nTT::Integer: Number of periods.\nndims::Integer: Number of grid points in each dimension.\nthetai::Float: Probability that infected person sends \"infected\" signal.\nthetas::Float: Probability that non-infected person sends \"infected\" signal.\n\n\n\n\n\n","category":"function"},{"location":"#HammerScalpel.nopolicy","page":"Replication of Chari, Kirpalani, and Phelan (2021)","title":"HammerScalpel.nopolicy","text":"nopolicy(;TT=52)\n\nSimulate the model in Chari, Kirpalani, and Phelan (2021) for TT periods under no government response.  The output is a Dict containing  the necessary variables and time-series for generating the figures. \n\n\n\n\n\n","category":"function"},{"location":"#HammerScalpel.withpolicy-Tuple{String}","page":"Replication of Chari, Kirpalani, and Phelan (2021)","title":"HammerScalpel.withpolicy","text":"withpolicy(policy; TT=52, ndims=20, thetai=0.38, thetas=0.0044)\nwithpolicy(policy, nprocs; TT=52, ndims=20, thetai=0.38, thetas=0.0044)\n\nSolve the model under a user-specified government policy.  In all cases, the government quarantines individuals that are known to be infected. The output is a Dict containing the necessary variables and time-series for generating the figures. \n\nArguments:\n\npolicy::String: Label for the government policy. Takes values\nnotest: no testing at all.\nuntargettest: random testing.\ntargettest: targeted testing under access to contact-tracing technology.\nisolate: no testing, but isolation of certain fraction of population.\nnprocs::Int64: Optional. If specified, uses the module Distributed for parallelizing the computations. Requires prior setup of nprocs workers.\nTT::Integer: Number of periods.\nndims::Integer: Number of grid points in each dimension.\nalgorithm: NLopt optimization algorithm.\nthetai::Float: Probability that infected person sends \"infected\" signal.\nthetas::Float: Probability that non-infected person sends \"infected\" signal.\n\n\n\n\n\n","category":"method"},{"location":"#HammerScalpel.createFig5-Tuple{Dict, Dict, Dict}","page":"Replication of Chari, Kirpalani, and Phelan (2021)","title":"HammerScalpel.createFig5","text":"createFig5(noint, notest, untarget)\ncreateFig5()\n\nReplicate Figure 5 in Chari, Kirpalani, and Phelan (2021).  If no arguments are provided, then stored results obtained with default parameters are used.\n\nArguments:\n\nnoint::Dict: Output from nopolicy.\nnotest::Dict: Output from withpolicy(::String) with policy notest.\nuntarget::Dict: Output from withpolicy(::String) with policy untarget.\n\n\n\n\n\n","category":"method"},{"location":"#HammerScalpel.createFig6-Tuple{Dict, Dict, Dict}","page":"Replication of Chari, Kirpalani, and Phelan (2021)","title":"HammerScalpel.createFig6","text":"createFig6(notest, untarget, target)\ncreateFig6()\n\nAnalogous to createFig5.\n\n\n\n\n\n","category":"method"},{"location":"#HammerScalpel.createFig7-Tuple{Dict, Dict, Dict}","page":"Replication of Chari, Kirpalani, and Phelan (2021)","title":"HammerScalpel.createFig7","text":"createFig7(untarget, target, isolate)\ncreateFig7()\n\nAnalogous to createFig5.\n\n\n\n\n\n","category":"method"},{"location":"#Implementation-details","page":"Replication of Chari, Kirpalani, and Phelan (2021)","title":"Implementation details","text":"","category":"section"},{"location":"","page":"Replication of Chari, Kirpalani, and Phelan (2021)","title":"Replication of Chari, Kirpalani, and Phelan (2021)","text":"Optimization algorithm: The authors originally use the algorithm SQP for maximizing the Bellman equation in each iteration, which in this case requires approximating gradients and Hessian. For some policies, the NLopt algorithm LN\\BOBYQA_ finds the same solution with greater speed. Therefore, although the replication results hold when using SQP everywhere, I use LN\\BOBYQA_ when applicable to improve speed.","category":"page"},{"location":"","page":"Replication of Chari, Kirpalani, and Phelan (2021)","title":"Replication of Chari, Kirpalani, and Phelan (2021)","text":"Dimensionality: The authors originally use 40 grid points for each control variable. In this replication, I reduced them to 20 for computational efficiency, which is why some of the time series look a bit ragged. ","category":"page"},{"location":"","page":"Replication of Chari, Kirpalani, and Phelan (2021)","title":"Replication of Chari, Kirpalani, and Phelan (2021)","text":"Inconsistencies: In my replications of Figures 6 and 7, the time series for R~0~ features a sharp increases for policies \"targeted\" and \"isolate\" around T=40 and T=50, respectively, which are not present in the published paper. I verified that these inconsistencies are also present in the original Matlab code.","category":"page"}]
}
