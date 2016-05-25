#this file is part of alk_cluster_2012
#Copyright (C) 2014 Ashok Litwin-Kumar
#see README for more information

#uncomment the line below and set doplot=true to plot a raster

include("sompsim.jl")

tic()
times,ns,Ne,Ni,Ncells,T = sompsim()
toc()

println("mean excitatory firing rate: ",mean(1000*ns[1:Ne]/T)," Hz")
println("mean inhibitory firing rate: ",mean(1000*ns[(Ne+1):(Ni+Ne)]/T)," Hz")
println("mean input firing rate: ",mean(1000*ns[(Ne+Ni+1):(Ncells)]/T)," Hz")

writedlm("times.txt", times)
writedlm("ns.txt", ns)
