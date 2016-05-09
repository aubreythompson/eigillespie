#this file is part of alk_cluster_2012
#Copyright (C) 2014 Ashok Litwin-Kumar
#see README for more information

#uncomment the line below and set doplot=true to plot a raster

using PyPlot
doplot = true

include("sompsim.jl")

tic()
times,ns,Ne,Ni,Ncells,T = sompsim()
toc()

println("mean excitatory firing rate: ",mean(1000*ns[1:Ne]/T)," Hz")
println("mean inhibitory firing rate: ",mean(1000*ns[(Ne+1):(Ni+Ne)]/T)," Hz")
println("mean input firing rate: ",mean(1000*ns[(Ne+Ni+1):(Ncells)]/T)," Hz")

if doplot
	println("creating plot")
	figure(figsize=(4,4))
	for ci = 1:(Ncells)
		vals = times[ci,1:ns[ci]]
		y = ci*ones(length(vals))
		scatter(vals,y,s=.3,c="k",marker="o",linewidths=0)
	end
	xlim(0,T)
	ylim(0,(Ncells))
	ylabel("Neuron")
	xlabel("Time")
	tight_layout()
	savefig("output.png",dpi=150)
end
