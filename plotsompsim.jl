using PyPlot

#import N parameters
Nns=readdlm("Documents/Piriform/code/Data/Nns6doubledensesmall.txt")
Nel=Int64(Nns[1])
Ner=Int64(Nns[2])
Npl=Int64(Nns[3])
Npr=Int64(Nns[4])
Nsl=Int64(Nns[5])
Nsr=Int64(Nns[6])
N0=Int64(Nns[7])
Ne=Nel+Ner
Np=Npl+Npr
Ns=Nsl+Nsr
Ni=Np+Ns
Ncells=Int64(sum(Nns))

#import spike times of neurons
stimrates=readdlm(string("Documents/Piriform/code/Data/11-2longstimrates1.txt"))
ns=readdlm(string("Documents/Piriform/code/Data/11-2longspikes1.txt"))
times=readdlm(string("Documents/Piriform/code/Data/11-2longtimes1.txt"))
T=maximum(times)

#plot raster
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

figure(figsize=(4,4))
ni, binsi=hist(cvinonnan,100)
bar(ni[2:end],binsi,.01)
xlabel("ISI CV")
ylabel("Number of I neurons")

#plot individual neuron firing rate histograms
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




#loop over combinations of connection probability ratios
params=readdlm("Documents/Piriform/code/Data/params2.txt")
paramnum=size(params,1)
valids=zeros(1,14) #save parameters which give higher E firing rate on the left
for i=11:paramnum
  println(i)
  stimrates=readdlm(string("Documents/Piriform/code/Data/",i,"-2longstimrates1.txt"))
  ns=readdlm(string("Documents/Piriform/code/Data/",i,"-2longspikes1.txt"))
  times=readdlm(string("Documents/Piriform/code/Data/",i,"-2longtimes1.txt"))
  T=maximum(times)

  lefr=mean(1000*ns[1:Nel,:]/T)
  refr=mean(1000*ns[(Nel+1):(Nel+Ner),:]/T)
  lpfr=mean(1000*ns[(Nel+Ner+1):(Nel+Ner+Npl),:]/T)
  rpfr=mean(1000*ns[(Nel+Ner+Npl+1):(Nel+Ner+Npl+Npr),:]/T)
  lsfr=mean(1000*ns[(Nel+Ner+Npl+Npr+1):(Nel+Ner+Npl+Npr+Nsl),:]/T)
  rsfr=mean(1000*ns[(Nel+Ner+Npl+Npr+Nsl+1):(Nel+Ner+Npl+Npr+Nsl+Nsr),:]/T)
  ifr=mean(1000*ns[(Nel+Ner+Npl+Npr+Nsl+Nsr):Ncells,:]/T)

  if lefr>refr
    valids=[valids; lefr refr lpfr rpfr lsfr rsfr ifr params[i,:] i]
  end
end

#Loop over Ns to see how firing rate converges
allNs=[4500,9000,13500,18000,22500,27000,31500,36000,40500,45000]
n=size(allNs,1)
frs=zeros(7,n) #save firing rates so we can plot them
for i=1:n
  N=allNs[i]
  stimrates=readdlm(string("Documents/Piriform/code/Data/stimrates.txt"))
  ns=readdlm(string("Documents/Piriform/code/Data/spikes.txt"))
  times=readdlm(string("Documents/Piriform/code/Data/times.txt"))
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

#plot firing rate convergence with N for all the pops
figure()
plot(allNs,frs[1,:]')
plot(allNs,frs[2,:]')
plot(allNs,frs[3,:]')
plot(allNs,frs[4,:]')
plot(allNs,frs[5,:]')
plot(allNs,frs[6,:]')


#compute sparsity of stimuli and selectivity of neurons (requires several runs with different stimuli)
stimrange=[1:20]
s=Int64(size(stimrange,1))
stimrates=zeros(N0,s)
ns=zeros(Ncells,s)

for stimnum=1:s
  stimrates[:,stimnum]=readdlm(string("Documents/Piriform/code/Data/stimrates",stimrange[stimnum],".txt"))
  ns[:,stimnum]=readdlm(string("Documents/Piriform/code/Data/spikes",stimrange[stimnum],".txt"))
end


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
