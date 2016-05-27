using PyPlot

Nns=readdlm("Documents/Piriform/code/piriform/Nns.txt")
Ne=Nns[1]
Ni=Nns[2]
N0=Nns[3]
Ncells=sum(Nns)

times=readdlm("Documents/Piriform/code/piriform/times1.txt")
ns=readdlm("Documents/Piriform/code/piriform/spikes1.txt")
stimrates=readdlm("Documents/Piriform/code/piriform/stimrates1.txt")

T=maximum(times)

println("mean excitatory firing rate: ",mean(1000*ns[1:Ne]/T)," Hz")
println("mean inhibitory firing rate: ",mean(1000*ns[(Ne+1):(Ni+Ne)]/T)," Hz")
println("mean input firing rate: ",mean(1000*ns[(Ne+Ni+1):(Ncells)]/T)," Hz")
println("mean input firing rate (theoretically): ",mean(stimrates)," Hz")

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

figure(figsize=(4,4))
ne, binse=hist(1000*ns[1:Ne]/T)
ni, binsi=hist(1000*ns[(Ne+1):(Ni+Ne)]/T)
show()
