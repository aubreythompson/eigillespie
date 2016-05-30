using PyPlot

Nns=readdlm("Documents/Piriform/code/Data/Nns.txt")
Ne=Nns[1]
Ni=Nns[2]
N0=Nns[3]
Ncells=sum(Nns)

stimnum=1
mi=2
stimrates=readdlm(string("Documents/Piriform/code/Data/stimrates",stimnum,".txt"))


ns=readdlm(string("Documents/Piriform/code/Data/spikes",stimnum,".txt"))
times=readdlm(string("Documents/Piriform/code/Data/times",stimnum,".txt"))
T=maximum(times)
for i=2:mi
  i=2
  newtimes =readdlm(string("Documents/Piriform/code/Data/times",i,"-",stimnum,".txt"))+T
  times=[times newtimes+T] #fix this so it doesn't make them all spike at T here
  ns=ns+readdlm(string("Documents/Piriform/code/Data/spikes",i,"-",stimnum,".txt"))
  T=maximum(times)
end



println("mean excitatory firing rate: ",mean(1000*ns[1:Ne]/T)," Hz")
println("mean inhibitory firing rate: ",mean(1000*ns[(Ne+1):(Ni+Ne)]/T)," Hz")
println("mean input firing rate: ",mean(1000*ns[(Ne+Ni+1):(Ncells)]/T)," Hz")
println("mean input firing rate (theoretically): ",mean(stimrates)," Hz")

println("creating plot")
figure(figsize=(4,4))
for ci = 1:(Ncells)
  (I,J,V)=findnz(times[ci,:])
  println(V)
 # vals = V #times[ci,1:ns[ci]]
  y = ci*ones(length(V))
  scatter(V,y,s=.3,c="k",marker="o",linewidths=0)
end
xlim(0,T)
ylim(0,(Ncells))
ylabel("Neuron")
xlabel("Time")
tight_layout()
savefig("output.png",dpi=150)

figure(figsize=(4,4))
ne, binse=hist(1000*ns[1:Ne]/T)
bar(binse, ne[2:end])
figure(figsize=(4,4))
ni, binsi=hist(1000*ns[(Ne+1):(Ni+Ne)]/T)
bar(binsi, ni[2:end])

