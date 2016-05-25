using PyPlot

Ncells=45000

times=readdlm("Documents/Piriform/code/piriform/times.txt")
ns=readdlm("Documents/Piriform/code/piriform/ns.txt")

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