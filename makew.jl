

Ne = 20000
Ni = 5000
N0 = 20000
Nns=[Ne,Ni,N0]
Ncells = Ne+Ni+N0
K = 1000.0 #average number of E->E connections per neuron

##tau and g are inversely related UNLESS THIS IS NOT 1
CEM = 1.0 #(μF/cm^2)
CIM = 1.0


#if you change gEL, you also need to change it in sompsim!
gEL = 0.05 #(mS/cm^2)
gIL = 0.05

taue = CEM/gEL #membrane time constant for exc. neurons (ms)
taui = CIM/gIL
taua = 100
tau1 = 3
tau2 = 1

VL = -65 #(mV) leak potential
Vth = -50 #spiking threshold
VE = 0 #excitatory reversal potential
VI = -80 #inhibitory RP
Vde = VE-VL
Vdi = VI-VL

#connection probabilities
pe0 = K/N0
pi0 = pe0
pii = K/Ni
pei = pii
pee = K/Ne
pie = pee

sqrtK = sqrt(K)

#	Nepop = 80

je0 = 4.846*Vde*gEL*tau1/sqrtK
ji0 = 3.808*Vde*gIL*tau1/sqrtK
jee = 0.462*Vde*gEL*tau1/sqrtK
jie = 1.270*Vde*gIL*tau1/sqrtK
jei = 10.5*Vdi*gEL*tau1/sqrtK
jii = 9.5*Vdi*gIL*tau1/sqrtK

J=[je0 ji0 jee jie jei jii]
tau = zeros(Ncells)

tau[1:Ne] = taue
tau[(1+Ne):Ncells] = taui

#weights = zeros(Ncells,Ncells)

projections=[0]
Ks=zeros(Ncells)

#random connections
# weights[1:Ne,1:Ne] = (rand(Ne,Ne) .< pee)#*(jee/sqrtK) #excitatory --> excitatory
# weights[1:Ne,(1+Ne):(Ni+Ne)] = (rand(Ne,Ni) .< pei)#*(jei/sqrtK) #inhibitory --> excitatory
# weights[(1+Ne):(Ni+Ne),1:Ne] = (rand(Ni,Ne) .< pie)#*(jie/sqrtK) #excitatory --> inhibitory
# weights[(1+Ne):(Ni+Ne),(1+Ne):(Ni+Ne)] = (rand(Ni,Ni) .< pii)#*(jii/sqrtK) #inhibitory --> inhibitory
# weights[(1+Ne):(Ni+Ne),(1+Ni+Ne):Ncells] = (rand(Ni,N0) .< pi0)#*(ji0/sqrtK) #input layer --> inhibitory
# weights[1:Ne,(1+Ni+Ne):Ncells] = (rand(Ne,N0) .< pe0)#*(je0/sqrtK) #input layer --> excitatory

weights=readdlm("Documents/Piriform/code/Data/weights.txt")
fw=(weights.!=0)
fwn=(weights.<0)
sum(fw,1)
Ks==sum(fw,1)'
using PyPlot
figure()
plot(1:45000,Ks-sum(fw,1)')
#Ks here are out-degrees, which is why they are not on average 1000 (does that make sense?)
#a neuron RECEIVES on average K inputs from each presyn population
#which would mean I neurons would have bigger out degrees since there are fewer of them
#this is definitely what we want because these are projections.
#indices sum(K[:i]) to sum(K[:i+1]) are downstream from neuron i+1
#excitatory projections
for i=1:Ne
  weightsee = (rand(Ne,1) .< pee) #excitatory --> excitatory
  weightsee[i] = 0
  weightsie = (rand(Ni,1) .< pie) #excitatory --> inhibitory
#   weightsee = weights[1:Ne,i]
#   weightsie = weights[(Ne+1):(Ne+Ni),i]
  fee=find(weightsee)
  fie=find(weightsie)
  Ks[i]=size(fee,1)+size(fie,1)
  projections=[projections; fee; fie+Ne]
end

writedlm("Documents/Piriform/code/Data/projections-weights.txt",projections)
writedlm("Documents/Piriform/code/Data/Ks-weights.txt",Ks)

#inhibitory projections
for i=1:Ni
  weightsei = (rand(Ne,1) .< pei) #inhibitory --> excitatory
  weightsii = (rand(Ni,1) .< pii) #inhibitory --> inhibitory
  weightsii[i] = 0
#   weightsei = weights[1:Ne,Ne+i]
#   weightsii = weights[(Ne+1):(Ne+Ni),Ne+i]
  fei=find(weightsei)
  fii=find(weightsii)
  Ks[i+Ne]=size(fei,1)+size(fii,1)
  projections= [projections; fei; fii+Ne]
end



#input projections
#projections=projections2
for i=1:N0
  weightse0 = (rand(Ne,1) .< pe0) #input layer --> inhibitory
  weightsi0 = (rand(Ni,1) .< pi0) #input layer --> excitatory
#   weightse0 = weights[1:Ne,Ne+Ni+i]
#   weightsi0 = weights[(Ne+1):(Ne+Ni),Ne+Ni+i]
  fe0=find(weightse0)
  fi0=find(weightsi0)
  Ks[i+Ne+Ni]=size(fe0,1)+size(fi0,1)
  projections= [projections; fe0; fi0+Ne]
end

projections=readdlm("Documents/Piriform/code/Data/projections.txt")
Ks=readdlm("Documents/Piriform/code/Data/Ks.txt")
writedlm("Documents/Piriform/code/piriform/J.txt", J)
writedlm("Documents/Piriform/code/piriform/tau.txt", tau)
writedlm("Documents/Piriform/code/piriform/Nns.txt",Nns)

# inKs=zeros(Ncells)
# sp=size(projections,1)
# for i=1:Ncells
#   fsp=size(find(projections-i*ones(sp)),1)
#   inKs[i]=sp-fsp
# end
# fsp=size(find(projections-5*ones(sp)),1)
# sp-fsp
# ci=1
# projci=projections[(sum(Ks[1:(ci-1)])+1):(sum(Ks[1:ci]))]
# Ks[ci]
# println(projci)
# weights[44,1]
# writedlm("Documents/Piriform/code/piriform/projections2.txt", projections)
# writedlm("Documents/Piriform/code/piriform/Ks2.txt", Ks)

