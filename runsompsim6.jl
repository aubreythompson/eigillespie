#this file is part of alk_cluster_2012
#Copyright (C) 2014 Ashok Litwin-Kumar
#modified by Aubrey Thompson 2016 to simulate
#Pehlevan, Sompolinsky: Selectivity and Sparseness in Randomly Connected Balanced Networks
#with excitatory, PV, and SOM populations on the left and right
#see README for more information


##Uncomment below to get the job ID
##Use this to run the same jobs with one varying parameter

# function getenv(var::ASCIIString)
#   val = ccall((:getenv, "libc"),
#               Ptr{Int8}, (Ptr{Int8},), var)
#   if val == C_NULL
#     error("getenv: undefined variable: ", var)
#   end
#   bytestring(val)
# end
#n = getenv("PBS_ARRAYID")
#stim = ccall((:getenv, "libc"), Ptr{Int64}, (Ptr{ASCIIString},), "PBS_ARRAYID")
#i=convert(Int64,parse(stim)) #this is the job array as an int
stim=0 #make a new stimulus
n=1

include("sompsim6.jl")

tic()
times,spikes,Ne,Np,Ns,Ncells,T,stimrates = sompsim(stim,n) #input 0 if you want a new stim
toc()

#make sure you don't save over current stims
println("Saving as strimrates",stim)
writedlm(string("times.txt"), times)
writedlm(string("spikes.txt"), spikes)
writedlm(string("stimrates.txt"),stimrates)
