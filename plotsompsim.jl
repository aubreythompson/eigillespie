using PyPlot

Nns=readdlm("Documents/Piriform/code/Data/Nns.txt")
Ne=Nns[1]
Ni=Nns[2]
N0=Nns[3]
Ncells=sum(Nns)


#concatenate data from multiple iterations of one stimulus
stimnum=2
mi=3 #maximum iterations
stimrates=readdlm(string("Documents/Piriform/code/Data/stimrates",stimnum,".txt"))

ns=readdlm(string("Documents/Piriform/code/Data/spikes5-",stimnum,".txt"))
times=readdlm(string("Documents/Piriform/code/Data/times5-",stimnum,".txt"))
T=maximum(times)
# for i=2:mi
#   println(i)
#   newtimes =readdlm(string("Documents/Piriform/code/Data/times",i,"-",stimnum,".txt"))+T
#   times=[times newtimes] #fix this so it doesn't make them all spike at T here
#   ns=ns+readdlm(string("Documents/Piriform/code/Data/spikes",i,"-",stimnum,".txt"))
#   T=maximum(times)
# end
println(T)
println("mean excitatory firing rate: ",mean(1000*ns[1:Ne]/T)," Hz")
println("mean inhibitory firing rate: ",mean(1000*ns[(Ne+1):(Ni+Ne)]/T)," Hz")
println("mean input firing rate: ",mean(1000*ns[(Ne+Ni+1):(Ncells)]/T)," Hz")
println("mean input firing rate (theoretically): ",mean(stimrates)," Hz")


#plot raster

println("creating plot")
figure(figsize=(4,4))
for ci = 1:(Ncells)
  (I,J,V)=findnz(times[ci,:])
  vals = V[1:ns[ci]]
  y = ci*ones(length(vals))
  scatter(vals,y,s=.3,c="k",marker="o",linewidths=0)
end
xlim(0,T)
ylim(0,(Ncells))
ylabel("Neuron")
xlabel("Time")
tight_layout()
savefig("output.png",dpi=150)

##get CV histograms

cve=zeros(Int64(Ne))
for ci=1:Ne
  (I,J,V)=findnz(times[ci,:])
  vals = V[1:ns[ci]]
  eisis=[0]
  for vi=2:size(vals,1)
    eisis=[eisis vals[vi]-vals[vi-1]]
  end
  cve[ci]=std(eisis[2:end])/mean(eisis[2:end])
end

cvenonnan=[0]
for i=1:size(cve,1)
  if !isnan(cve[i])
    cvenonnan=[cvenonnan cve[i]]
  end
end

cvenonnan=cvenonnan[find(cvenonnan)]
figure(figsize=(4,4))
ne, binse=hist(cvenonnan,100)
bar(ne[2:end],binse,.01)
xlabel("ISI CV")
ylabel("Number of E neurons")

cvi=zeros(Int64(Ne))
for ci=(Ne+1):(Ne+Ni)
  (I,J,V)=findnz(times[ci,:])
  vals = V[1:ns[ci]]
  iisis=[0]
  for vi=2:size(vals,1)
    iisis=[iisis vals[vi]-vals[vi-1]]
  end
  cvi[ci-Ne]=std(iisis[2:end])/mean(iisis[2:end])
end

cvinonnan=[0]
for i=1:size(cvi,1)
  if !isnan(cvi[i])
    cvinonnan=[cvinonnan cvi[i]]
  end
end

cvinonnan=cvinonnan[find(cvinonnan)]

print(cvinonnan)

figure(figsize=(4,4))
ni, binsi=hist(cvinonnan,100)
bar(ni[2:end],binsi,.01)
xlabel("ISI CV")
ylabel("Number of I neurons")


figure(figsize=(4,4))
nefr, binsefr=hist(1000*ns[1:Ne]/T,50)
bar(nefr[2:end],binsefr)
xlabel("Firing rate")
ylabel("Number of E neurons")
figure(figsize=(4,4))
nifr, binsifr=hist(1000*ns[(Ne+1):(Ni+Ne)]/T,50)
bar(nifr[2:end],binsifr)
xlabel("Firing rate")
ylabel("Number of I neurons")

# figure(figsize=(4,4))
# bar(linspace(1,Ne,Ne),1000*ns[1:Ne]/T)
# figure(figsize=(4,4))
# bar(linspace(1,Ni,Ni),1000*ns[(Ne+1):(Ni+Ne)]/T)


#compute spareness index (SPI)
avgsq=(sum(1000*ns[1:Ne]/T)/Ne)^2
sqavg=(sum((1000*ns[1:Ne]/T).^2)/Ne)
espi=(1/(1-(1/Ne)))*(1-avgsq/sqavg)

avgsq=(sum(1000*ns[(Ne+1):(Ni+Ne)]/T)/Ni)^2
sqavg=(sum((1000*ns[(Ne+1):(Ni+Ne)]/T).^2)/Ni)
ispi=(1/(1-(1/Ni)))*(1-avgsq/sqavg)
