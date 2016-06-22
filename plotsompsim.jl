using PyPlot

Nns=readdlm("Documents/Piriform/code/Data/Nns.txt")
Ne=Int64(Nns[1])
Ni=Int64(Nns[2])
N0=Int64(Nns[3])
Ncells=Int64(sum(Nns))

maxrate=200
T=1000
maxTimes = Int64(maxrate*T/1000)
stimrange=3:48
s=Int64(size(stimrange,1))
stimrates=zeros(N0,s)
ns=zeros(Ncells,s)
times=zeros(Ncells,maxTimes,s)
maxts=zeros(s)

# for stimnum=1:s
#   #stimrates[:,stimnum]=readdlm(string("Documents/Piriform/code/Data/stimrates",stimrange[stimnum],".txt"))
#   #ns[:,stimnum]=readdlm(string("Documents/Piriform/code/Data/spikes",stimrange[stimnum],".txt"))
#   t=readdlm(string("Documents/Piriform/code/Data/times",stimrange[stimnum],".txt"))
#   maxts[stimnum]=size(t,2)
#   #times[:,1:size(t,2),stimnum]=t
# end

stimrates=readdlm("Documents/Piriform/code/Data/stimratesall.txt")
ns=readdlm("Documents/Piriform/code/Data/spikesall.txt")
times=readdlm("Documents/Piriform/code/Data/timesall.txt")
maxts=readdlm("Documents/Piriform/code/Data/maxtsall.txt")

# figure()
# plot(1:s,maxts)

# wns=readdlm(string("Documents/Piriform/code/Data/wspikes",stimnum,".txt"))
# wtimes=readdlm(string("Documents/Piriform/code/Data/wtimes",stimnum,".txt"))

#weights=readdlm(string("Documents/Piriform/code/Data/weights.txt"))

#wT=maximum(wtimes)
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

# println("mean excitatory firing rate: ",mean(1000*wns[1:Ne]/T)," Hz")
# println("mean inhibitory firing rate: ",mean(1000*wns[(Ne+1):(Ni+Ne)]/T)," Hz")
# println("mean input firing rate: ",mean(1000*wns[(Ne+Ni+1):(Ncells)]/T)," Hz")
#println("mean input firing rate (theoretically): ",mean(stimrates)," Hz")
# v=readdlm(string("Documents/Piriform/code/Data/wv1.txt"))
# figure()
# pcolormesh(weights)
# colorbar()

#plot raster

# wv11=readdlm(string("Documents/Piriform/code/Data/wv7",stimnum,".txt"))
# #aI11=readdlm(string("Documents/Piriform/code/Data/aI1",stimnum,".txt"))
# v1=readdlm(string("Documents/Piriform/code/Data/v1",stimnum,".txt"))
# #waI11=readdlm(string("Documents/Piriform/code/Data/waI1",stimnum,".txt"))
# wfIE1=readdlm(string("Documents/Piriform/code/Data/wfIE1",stimnum,".txt"))
# wfII1=readdlm(string("Documents/Piriform/code/Data/wfII1",stimnum,".txt"))
# fIE1=readdlm(string("Documents/Piriform/code/Data/fIE1",stimnum,".txt"))
# fII1=readdlm(string("Documents/Piriform/code/Data/fII1",stimnum,".txt"))

# wfI01=readdlm(string("Documents/Piriform/code/Data/wfI01",stimnum,".txt"))
# fI01=readdlm(string("Documents/Piriform/code/Data/fI01",stimnum,".txt"))
# fdi1=readdlm(string("Documents/Piriform/code/Data/fdi1",stimnum,".txt"))
# wfdi1=readdlm(string("Documents/Piriform/code/Data/wfdi1",stimnum,".txt"))
# fdi1=reshape(fdi1, (20000,Int64(Ncells),Int64(Ncells)))
# wfdi1=reshape(wfdi1, (20000,Int64(Ncells),Int64(Ncells)))
# maxT=size(fdi1,1)
# for i=1:size(fdi1,2)
#   figure()
#   plot(1:maxT,sum(wfdi1[:,i,1:Ne],3)[:,1,1])
#   plot(1:maxT,sum(fdi1[:,i,1:Ne],3)[:,1,1])
# end


# for i=1:size(wfIE1,2)
#   figure()
#   plot(1:maxT,fIE1[:,i])
#   plot(1:maxT,wfIE1[:,i])
# end

# figure()
# plot(1:maxT,wfII1[:,1])
# plot(1:maxT,fII1[:,1])


# # figure()
# # plot(1:2000,waI11)
# # plot(1:2000,aI11)

# for i=1:size(wfI01,2)
#   figure()
#   plot(1:maxT,fI01[:,i])
#   plot(1:maxT,wfI01[:,i])
# end


# figure()
# plot(500:2500,wfI01[500:2500,1])
# plot(500:2500,fI01[500:2500,1])

# figure()
# plot(500:2500,wfII1[500:2500,1])
# plot(500:2500,fII1[500:2500,1])

# figure()
# plot(500:2500,wfI01[500:2500,2])
# plot(500:2500,fI01[500:2500,2])


# figure()
# plot(1:2500,wfII1[1:2500,1])
# plot(1:2500,fII1[1:2500,1])
# plot(1:2500,wfI01[1:2500,1])
# plot(1:2500,fI01[1:2500,1])
# plot(1:2500,wfIE1[1:2500,2])
# plot(1:2500,fIE1[1:2500,2])
# plot(1:2500,wv11[1:2500,1])
# plot(1:2500,v1[1:2500,1])

# timerange=1:maxT

# figure()
# plot(timerange,wv11[timerange,2])
# plot(timerange,v1[timerange,2])
# plot(timerange,wv11[timerange,1])
# plot(timerange,v1[timerange,1])
# plot(timerange,wfI01[timerange,2])
# plot(timerange,fI01[timerange,2])
# plot(timerange,wfI01[timerange,1])
# plot(timerange,fI01[timerange,1])

# plot(1:maxT,wfIE1[:,2])

# plot(1:maxT,fII1[:,2])
# plot(1:maxT,fIE1[:,2])
# plot(1:maxT,v1[:,2])

# figure()
# # plot(1:maxT,fdi1[:,3,1])
# # plot(1:maxT,fdi1[:,3,2])
# plot(1:2500,fdi1[1:2500,1,2])
# plot(1:2500,wfdi1[1:2500,1,2])

# figure()
# plot(1:maxT,fI01[:,2])
# plot(1:maxT,wfI01[:,2])

# plot(1:maxT,wfII1[:,2])
# plot(1:maxT,wv11[:,2])

# figure()
# s=size(find(sum(wfdi1,1)),1)
# plot(1:s,find(sum(wfdi1,1)))

# s=size(find(sum(fdi1,1)),1)
# plot(1:s,find(sum(fdi1,1)))

println("creating plot")
figure(figsize=(4,4))
for ci = 1:(Ncells)
  (I,J,V)=findnz(times[ci,:,1])
  vals = V[1:ns[ci,1]]
  y = ci*ones(length(vals))
  scatter(vals,y,s=.3,c="k",marker="o",linewidths=0)
end
xlim(0,T)
ylim(0,(Ncells))
ylabel("Neuron")
xlabel("Time")
tight_layout()
savefig("output.png",dpi=150)

v=1
initialv=v
v=2
println(initialv)
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

espi=zeros(s)
ispi=zeros(s)
for i=1:s
  avgsq=(sum(1000*ns[1:Ne,i]/T)/Ne)^2
  sqavg=(sum((1000*ns[1:Ne,i]/T).^2)/Ne)
  espi[i]=(1/(1-(1/Ne)))*(1-avgsq/sqavg)

  avgsq=(sum(1000*ns[(Ne+1):(Ni+Ne),i]/T)/Ni)^2
  sqavg=(sum((1000*ns[(Ne+1):(Ni+Ne),i]/T).^2)/Ni)
  ispi[i]=(1/(1-(1/Ni)))*(1-avgsq/sqavg)
end

figure()
plot(1:s,espi)
plot(1:s,mean(espi)*ones(s))
plot(1:s,ispi)
plot(1:s,mean(ispi)*ones(s))

#compute selectivity index (SLI)

sli=zeros(Ne+Ni)
for i=1:(Ne+Ni)
  avgsq=(sum(1000.0*ns[i,:]/T)/s)^2
  sqavg=(sum((1000.0*ns[i,:]/T).^2)/s)
  if sqavg>0
    sli[i]=(1.0/(1-(1.0/s)))*(1.0-avgsq/sqavg)
  else
    sli[i]=0
  end
end

sorte=sortperm(sli[1:Ne])
sorti=sortperm(sli[(Ne+1):(Ne+Ni)])

figure()
plot(1:s,ns[sorte[18999],:]')

figure()
plot(1:s,ns[sorti[5000],:]')

figure()
plot(1:Ne,sort(sli[1:Ne]))
plot(1:Ne,mean(sli[1:Ne])*ones(Ne))
figure()
plot(1:Ni,sort(sli[(Ne+1):(Ne+Ni)]))
plot(1:Ni,mean(sli[(Ne+1):(Ne+Ni)])*ones(Ni))

figure(figsize=(4,4))
nesli, binsesli=hist(sli[1:Ne],20)
bar(nesli[2:end],binsesli,.01)
xlabel("SLI")
ylabel("Number of E neurons")
figure(figsize=(4,4))
nisli, binsisli=hist(sli[(Ne+1):(Ni+Ne)],20)
bar(nisli[2:end],binsisli,.01)
xlabel("SLI")
ylabel("Number of I neurons")

#SLI vs rate at preferred stim

prefstim=zeros(Ne+Ni)
ratepref=zeros(Ne+Ni)
for i=1:(Ne+Ni)
  prefstim[i]=indmax(ns[i,:])
  ratepref[i]=maximum(ns[i,:])
end

figure()
scatter(sli[1:Ne],ratepref[1:Ne])

figure()
scatter(sli[(Ne+1):(Ne+Ni)],ratepref[(Ne+1):(Ne+Ni)])

# stimtimes=times[(Ne+Ni+1):Ncells,:]

# writedlm("Documents/Piriform/code/Data/stimtimes.txt",stimtimes)

# findfirst(ones(2),1)
