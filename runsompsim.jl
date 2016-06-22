#this file is part of alk_cluster_2012
#Copyright (C) 2014 Ashok Litwin-Kumar
#see README for more information


include("sompsim.jl")

tic()
times,ns,Ne,Ni,Ncells,T,stimrates = sompsim(0)
toc()

writedlm("times.txt", times)
writedlm("ns.txt", ns)
writedlm("stimrates.txt", stimrates)
