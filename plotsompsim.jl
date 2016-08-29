using PyPlot

Nns=readdlm("Documents/Piriform/code/Data/Nns6doubledensesmall.txt")
Nel=Int64(Nns[1])
Ner=Int64(Nns[2])
Npl=Int64(Nns[3])
Npr=Int64(Nns[4])
Nsl=Int64(Nns[5])
Nsr=Int64(Nns[6])
N0=Int64(Nns[7])
Ncells=Int64(sum(Nns))

Nns=readdlm("Documents/Piriform/code/Data/Nns3.txt")
Ne=Int64(Nns[1])
Np=Int64(Nns[2])
Ns=Int64(Nns[3])
N0=Int64(Nns[4])
Ncells=Int64(sum(Nns))

maxrate=400
T=5000
maxTimes = Int64(maxrate*T/1000)
stimrange=[1:20]
s=Int64(size(stimrange,1))
stimrates=zeros(N0,s)
ns=zeros(Ncells,s)
# times=zeros(Ncells,maxTimes,s)
# maxts=zeros(s)

for stimnum=1:s
  stimrates[:,stimnum]=readdlm(string("Documents/Piriform/code/Data/195959stimrates",stimrange[stimnum],".txt"))
  ns[:,stimnum]=readdlm(string("Documents/Piriform/code/Data/195959spikes",stimrange[stimnum],".txt"))
#   t=readdlm(string("Documents/Piriform/code/Data/3times",stimrange[stimnum],".txt"))
#   maxts[stimnum]=size(t,2)
#   times[:,1:size(t,2),stimnum]=t
end

ind=14
stimrates=readdlm(string("Documents/Piriform/code/Data/",ind,"-full10stimrates1.txt"))
ns=readdlm(string("Documents/Piriform/code/Data/",ind,"-full10spikes1.txt"))
times=readdlm(string("Documents/Piriform/code/Data/",ind,"-full10times1.txt"))
T=maximum(times)
N=45000
Nel = convert(Int64,N*2.0/9)
Ner = convert(Int64,N*2.0/9)
Ne=Nel+Ner
Npl = convert(Int64,N*1.0/36)
Npr = convert(Int64,N*1.0/36)
Np=Npl+Npr
Nsl = convert(Int64,round(N*1.0/27))
Nsr = convert(Int64,round(N*1.0/54))
Ns=Nsl+Nsr
Ni=Np+Ns
N0 = convert(Int64,N*4.0/9)
println([Nel,Ner,Npl,Npr,Nsl,Nsr,N0])
lefr=mean(1000*ns[1:Nel,:]/T)
refr=mean(1000*ns[(Nel+1):(Nel+Ner),:]/T)
lpfr=mean(1000*ns[(Nel+Ner+1):(Nel+Ner+Npl),:]/T)
rpfr=mean(1000*ns[(Nel+Ner+Npl+1):(Nel+Ner+Npl+Npr),:]/T)
lsfr=mean(1000*ns[(Nel+Ner+Npl+Npr+1):(Nel+Ner+Npl+Npr+Nsl),:]/T)
rsfr=mean(1000*ns[(Nel+Ner+Npl+Npr+Nsl+1):(Nel+Ner+Npl+Npr+Nsl+Nsr),:]/T)
ifr=mean(1000*ns[(Nel+Ner+Npl+Npr+Nsl+Nsr):(N),:]/T)
println("mean left excitatory firing rate: ",lefr," Hz")
println("mean right excitatory firing rate: ",refr," Hz")
println("mean left PV firing rate: ",lpfr," Hz")
println("mean right PV firing rate: ",rpfr," Hz")
println("mean left SOM firing rate: ",lsfr," Hz")
println("mean right SOM firing rate: ",rsfr," Hz")
println("mean input firing rate: ",ifr," Hz")
println("mean input firing rate (theoretically): ",mean(stimrates)," Hz")
#maxts=readdlm("Documents/Piriform/code/Data/maxtsall.txt")
println(params[ind,:])

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


Nns=readdlm("Documents/Piriform/code/Data/Nns6doubledense.txt")
Nel=Int64(Nns[1])
Ner=Int64(Nns[2])
Npl=Int64(Nns[3])
Npr=Int64(Nns[4])
Nsl=Int64(Nns[5])
Nsr=Int64(Nns[6])
N0=Int64(Nns[7])
Ncells=Int64(sum(Nns))
params=readdlm("Documents/Piriform/code/Data/params2.txt")

paramnum=size(params,1)

validfulls=zeros(1,14)
#validinds=readdlm("Documents/Piriform/code/Data/validinds.txt")
#validnum=size(validinds,1)
for i=26:34
  #j=convert(Int64,validinds[i])
  println(i)
  stimrates=readdlm(string("Documents/Piriform/code/Data/",i,"-full10stimrates1.txt"))
  ns=readdlm(string("Documents/Piriform/code/Data/",i,"-full10spikes1.txt"))
  times=readdlm(string("Documents/Piriform/code/Data/",i,"-full10times1.txt"))
  T=maximum(times)

  lefr=mean(1000*ns[1:Nel,:]/T)
  refr=mean(1000*ns[(Nel+1):(Nel+Ner),:]/T)
  lpfr=mean(1000*ns[(Nel+Ner+1):(Nel+Ner+Npl),:]/T)
  rpfr=mean(1000*ns[(Nel+Ner+Npl+1):(Nel+Ner+Npl+Npr),:]/T)
  lsfr=mean(1000*ns[(Nel+Ner+Npl+Npr+1):(Nel+Ner+Npl+Npr+Nsl),:]/T)
  rsfr=mean(1000*ns[(Nel+Ner+Npl+Npr+Nsl+1):(Nel+Ner+Npl+Npr+Nsl+Nsr),:]/T)
  ifr=mean(1000*ns[(Nel+Ner+Npl+Npr+Nsl+Nsr):Ncells,:]/T)

  if lefr>refr
    validfulls=[validfulls; lefr refr lpfr rpfr lsfr rsfr ifr params[i,:] i]
  end
end

size(validfulls)
println(validfulls)
validlonginds=validlongs[:,14]
println(validlonginds)
writedlm("Documents/Piriform/code/Data/validlongs.txt",validlonginds[2:end])

allNs=[4500,9000,13500,18000,22500,27000,31500,36000,40500,45000]
n=size(allNs,1)
#frs3=zeros(4,n)
frs16asym=zeros(7,n)
#frs2=zeros(3,n)
for i=1:n
  N=allNs[i]
#   stimrates=readdlm(string("Documents/Piriform/code/Data/",i,"-2stimrates1.txt"))
#   ns=readdlm(string("Documents/Piriform/code/Data/",i,"-2spikes1.txt"))
#   times=readdlm(string("Documents/Piriform/code/Data/",i,"-2times1.txt"))
#   T=maximum(times)
#   println(size(ns,1))
#   Ne = convert(Int64,N*4.0/9)
#   Ni = convert(Int64,N*1.0/9)
#   N0 = convert(Int64,N*4.0/9)
#   efr=mean(1000*ns[1:Ne]/T)
#   inhfr=mean(1000*ns[(Ne+1):(Ni+Ne)]/T)
#   ifr=mean(1000*ns[(Ne+Ni+1):N]/T)
#   frs2[:,i]=[efr,inhfr,ifr]
#   stimrates=readdlm(string("Documents/Piriform/code/Data/",i,"-3-1stimrates1.txt"))
#   ns=readdlm(string("Documents/Piriform/code/Data/",i,"-3-1spikes1.txt"))
#   times=readdlm(string("Documents/Piriform/code/Data/",i,"-3-1times1.txt"))
#   T=maximum(times)
#   println(size(ns,1))
#   Ne = convert(Int64,N*4.0/9)
#   Np = convert(Int64,N*1.0/18)
#   Ns = convert(Int64,N*1.0/18)
#   N0 = convert(Int64,N*4.0/9)
#   efr=mean(1000*ns[1:Ne]/T)
#   pfr=mean(1000*ns[(Ne+1):(Np+Ne)]/T)
#   sfr=mean(1000*ns[(Ne+Np+1):(Np+Ns+Ne)]/T)
#   ifr=mean(1000*ns[(Ne+Np+Ns+1):N]/T)
#   frs3[:,i]=[efr,pfr,sfr,ifr]
  stimrates=readdlm(string("Documents/Piriform/code/Data/",i,"-1-asymstimrates1.txt"))
  ns=readdlm(string("Documents/Piriform/code/Data/",i,"-1-asymspikes1.txt"))
  times=readdlm(string("Documents/Piriform/code/Data/",i,"-1-asymtimes1.txt"))
  T=maximum(times)

  Nel = convert(Int64,N*2.0/9)
  Ner = convert(Int64,N*2.0/9)
  Ne=Nel+Ner
  Npl = convert(Int64,N*1.0/36)
  Npr = convert(Int64,N*1.0/36)
  Np=Npl+Npr
  Nsl = convert(Int64,round(N*1.0/27))
  Nsr = convert(Int64,round(N*1.0/54))
  Ns=Nsl+Nsr
  Ni=Np+Ns
  N0 = convert(Int64,N*4.0/9)
  println([Nel,Ner,Npl,Npr,Nsl,Nsr,N0])

  lefr=mean(1000*ns[1:Nel,:]/T)
  refr=mean(1000*ns[(Nel+1):(Nel+Ner),:]/T)
  lpfr=mean(1000*ns[(Nel+Ner+1):(Nel+Ner+Npl),:]/T)
  rpfr=mean(1000*ns[(Nel+Ner+Npl+1):(Nel+Ner+Npl+Npr),:]/T)
  lsfr=mean(1000*ns[(Nel+Ner+Npl+Npr+1):(Nel+Ner+Npl+Npr+Nsl),:]/T)
  rsfr=mean(1000*ns[(Nel+Ner+Npl+Npr+Nsl+1):(Nel+Ner+Npl+Npr+Nsl+Nsr),:]/T)
  ifr=mean(1000*ns[(Nel+Ner+Npl+Npr+Nsl+Nsr):N,:]/T)
  frs16asym[:,i]=[lefr,refr,lpfr,rpfr,lsfr,rsfr,ifr]
end


writedlm("Documents/Piriform/code/Data/frs2-1.txt",frs2)
writedlm("Documents/Piriform/code/Data/frs3-1.txt",frs3)
writedlm("Documents/Piriform/code/Data/frs6-1.txt",frs6)
writedlm("Documents/Piriform/code/Data/frs6asym-419.txt",frs6asym)


figure()
plot(allNs,frs2[1,:]')
plot(allNs,frs2[2,:]')
#plot(allNs,frs2[3,:]')

figure()
plot(allNs,frs3[1,:]')
plot(allNs,frs3[2,:]')
plot(allNs,frs3[3,:]')
#plot(allNs,frs3[4,:]')

figure()
plot(allNs,frs6[1,:]')
plot(allNs,frs6[2,:]')
plot(allNs,frs6[3,:]')
plot(allNs,frs6[4,:]')
plot(allNs,frs6[5,:]')
plot(allNs,frs6[6,:]')
#plot(allNs,frs6[7,:]')

figure()
plot(allNs,frs16asym[1,:]')
plot(allNs,frs16asym[2,:]')
plot(allNs,frs16asym[3,:]')
plot(allNs,frs16asym[4,:]')
plot(allNs,frs16asym[5,:]')
plot(allNs,frs16asym[6,:]')
#plot(allNs,frs6[7,:]')

figure()
plot(allNs,frs737[1,:]')
plot(allNs,frs737[2,:]')
plot(allNs,frs737[3,:]')
plot(allNs,frs737[4,:]')
plot(allNs,frs737[5,:]')
plot(allNs,frs737[6,:]')
#plot(allNs,frs6[7,:]')


missed=[33,45,151,152,176,206,292,333,410,412,418,419,421,437,440,492,495,520,541,549,567,603,609,617,625,626,635,639,679,701,721,751,783,818,825,835,845,862,902,912,949,965,975,1016,1018,1047]
size(valids)

print(valids)

writedlm("Documents/Piriform/code/Data/validsip2.txt",valids)

stimrates=readdlm(string("Documents/Piriform/code/Data/3-rs1stimrates1.txt"))
ns=readdlm(string("Documents/Piriform/code/Data/3-rs1spikes1.txt"))
times=readdlm(string("Documents/Piriform/code/Data/3-rs1times1.txt"))
T=maximum(times)
lefr=mean(1000*ns[1:Nel,:]/T)
refr=mean(1000*ns[(Nel+1):(Nel+Ner),:]/T)
lpfr=mean(1000*ns[(Nel+Ner+1):(Nel+Ner+Npl),:]/T)
rpfr=mean(1000*ns[(Nel+Ner+Npl+1):(Nel+Ner+Npl+Npr),:]/T)
lsfr=mean(1000*ns[(Nel+Ner+Npl+Npr+1):(Nel+Ner+Npl+Npr+Nsl),:]/T)
rsfr=mean(1000*ns[(Nel+Ner+Npl+Npr+Nsl+1):(Nel+Ner+Npl+Npr+Nsl+Nsr),:]/T)
ifr=mean(1000*ns[(Nel+Ner+Npl+Npr+Nsl+Nsr):(Ncells),:]/T)
println("mean left excitatory firing rate: ",lefr," Hz")
println("mean right excitatory firing rate: ",refr," Hz")
println("mean left PV firing rate: ",lpfr," Hz")
println("mean right PV firing rate: ",rpfr," Hz")
println("mean left SOM firing rate: ",lsfr," Hz")
println("mean right SOM firing rate: ",rsfr," Hz")
println("mean input firing rate: ",ifr," Hz")
println("mean input firing rate (theoretically): ",mean(stimrates)," Hz")

println("mean excitatory firing rate: ",mean(1000*ns[1:Ne]/T)," Hz")
println("mean PV firing rate: ",mean(1000*ns[(Ne+1):(Np+Ne)]/T)," Hz")
println("mean SOM firing rate: ",mean(1000*ns[(Ne+Np+1):(Np+Ns+Ne)]/T)," Hz")
println("mean input firing rate: ",mean(1000*ns[(Ne+Np+Ns+1):(Ncells)]/T)," Hz")
println("mean input firing rate (theoretically): ",mean(stimrates)," Hz")

allNs=[4500,9000,13500,18000,22500,27000,31500,36000,40500,45000]
writedlm("Documents/Piriform/code/Data/allNs.txt",allNs)
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
using PyPlot
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
Ni=Np+Ns
espi=zeros(s)
pspi=zeros(s)
sspi=zeros(s)
for i=1:s
  avgsq=(sum(1000*ns[1:Ne,i]/T)/Ne)^2
  sqavg=(sum((1000*ns[1:Ne,i]/T).^2)/Ne)
  espi[i]=(1.0/(1-(1.0/Ne)))*(1-avgsq/sqavg)

  avgsq=(sum(1000*ns[(Ne+1):(Np+Ne),i]/T)/Np)^2
  sqavg=(sum((1000*ns[(Ne+1):(Np+Ne),i]/T).^2)/Np)
  pspi[i]=(1.0/(1-(1.0/Np)))*(1-avgsq/sqavg)

  avgsq=(sum(1000*ns[(Ne+Np+1):(Ni+Ne),i]/T)/Ns)^2
  sqavg=(sum((1000*ns[(Ne+Np+1):(Ni+Ne),i]/T).^2)/Ns)
  sspi[i]=(1.0/(1-(1.0/Ns)))*(1-avgsq/sqavg)
end

figure()
plot(1:s,espi)
plot(1:s,mean(espi)*ones(s))
plot(1:s,pspi)
plot(1:s,mean(pspi)*ones(s))
plot(1:s,sspi)
plot(1:s,mean(sspi)*ones(s))
plot(1:s,mean([pspi,sspi])*ones(s))

s1=1
figure()
plot((Ne+1):(Ne+Np),ns[(Ne+1):(Ne+Np),s1])
pspi[s1]
espi[s1]

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
sortp=sortperm(sli[(Ne+1):(Ne+Np)])
sorts=sortperm(sli[(Ne+Np+1):(Ne+Ni)])
figure()
plot(1:s,ns[sorte[20000],:]')

figure()
plot(1:s,ns[Ne+sortp[2500],:]')

figure()
plot(1:s,ns[Ne+Np+sorts[2500],:]')

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
nislp, binsislp=hist(sli[(Ne+1):(Np+Ne)],20)
bar(nislp[2:end],binsislp,.01)
xlabel("SLI")
ylabel("Number of PV neurons")
figure(figsize=(4,4))
nisls, binsisls=hist(sli[(Ne+Np+1):(Ni+Ne)],20)
bar(nisls[2:end],binsisls,.01)
xlabel("SLI")
ylabel("Number of SOM neurons")


#SLI vs rate at preferred stim

prefstim=zeros(Ne+Ni)
ratepref=zeros(Ne+Ni)
for i=1:(Ne+Ni)
  prefstim[i]=indmax(ns[i,:])
  ratepref[i]=maximum(ns[i,:])
end

figure()
scatter(ratepref[1:Ne],sli[1:Ne])

figure()
scatter(ratepref[(Ne+1):(Ne+Np)],sli[(Ne+1):(Ne+Np)])

figure()
scatter(ratepref[(Ne+Np+1):(Ne+Ni)],sli[(Ne+Np+1):(Ne+Ni)])

# stimtimes=times[(Ne+Ni+1):Ncells,:]

# writedlm("Documents/Piriform/code/Data/stimtimes.txt",stimtimes)

# findfirst(ones(2),1)
