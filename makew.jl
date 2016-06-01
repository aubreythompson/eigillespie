

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
# weights[1:Ne,1:Ne] = (jee/sqrtK)*(rand(Ne,Ne) .< pee) #excitatory --> excitatory
# weights[1:Ne,(1+Ne):(Ni+Ne)] = (jei/sqrtK)*(rand(Ne,Ni) .< pei) #inhibitory --> excitatory
# weights[(1+Ne):(Ni+Ne),1:Ne] = (jie/sqrtK)*(rand(Ni,Ne) .< pie) #excitatory --> inhibitory
# weights[(1+Ne):(Ni+Ne),(1+Ne):(Ni+Ne)] = (jii/sqrtK)*(rand(Ni,Ni) .< pii) #inhibitory --> inhibitory
# weights[(1+Ne):(Ni+Ne),(1+Ni+Ne):Ncells] = (ji0/sqrtK)*(rand(Ni,N0) .< pi0) #input layer --> inhibitory
# weights[1:Ne,(1+Ni+Ne):Ncells] = (je0/sqrtK)*(rand(Ne,N0) .< pe0) #input layer --> excitatory

#Ks here are out-degrees, which is why they are not on average 1000
#this is definitely what we want because these are projections.
#indices sum(K[:i]) to sum(K[:i+1]) are downstream from neuron i+1
#excitatory projections
for i=1:Ne
  weightsee = (rand(1,Ne) .< pee) #excitatory --> excitatory
  weightsee[i] = 0
  weightsie = (rand(1,Ni) .< pie) #excitatory --> inhibitory
  Ks[i]=size(find(weightsee),1)+size(find(weightsie),1)
  projections=[projections; find(weightsee); find(weightsie)+Ne]
end
#inhibitory projections
for i=1:Ni
  weightsei = (rand(1,Ne) .< pei) #inhibitory --> excitatory
  weightsii = (rand(1,Ni) .< pii) #inhibitory --> inhibitory
  weightsii[i] = 0
  Ks[i+Ne]=size(find(weightsei),1)+size(find(weightsii),1)
  projections= [projections; find(weightsei); find(weightsii)+Ne]
end
#input projections
for i=1:N0
  weightse0 = (rand(1,Ne) .< pe0) #input layer --> inhibitory
  weightsi0 = (rand(1,Ni) .< pi0) #input layer --> excitatory
  Ks[i+Ne+Ni]=size(find(weightse0),1)+size(find(weightsi0),1)
  projections= [projections; find(weightse0); find(weightsi0)+Ne]
end


projci=projections[(sum(Ks[1:0])+1):sum(Ks[1:1])]
size(projci)

writedlm("Documents/Piriform/code/piriform/projections.txt", projections)
writedlm("Documents/Piriform/code/piriform/Ks.txt", Ks)
writedlm("Documents/Piriform/code/piriform/J.txt", J)
writedlm("Documents/Piriform/code/piriform/tau.txt", tau)
writedlm("Documents/Piriform/code/piriform/Nns.txt",Nns)

