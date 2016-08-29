Ne = 20000
Np = 2500
Ns = 2500
N0 = 20000
Nns=[Ne,Np,Ns,N0]
Ncells = Ne+Np+Ns+N0
writedlm("Documents/Piriform/code/Data/Nns3.txt",Nns)
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

rpse=0.9 #ratio of jse:jpe (excitatory connection to SOM:to PV)
reps=0.2 #ratio of jes:jep (excitatory connection from SOM:from PV)
rips=0.9 #ratio of jss:jpp (interconnection within SOM:within PV)
rrps=0.2 #ratio of jsp:jps (reciprocal connections between i pops from PV:from SOM)
rps0=1.0 #ratio of js0:jp0 (input layer to SOM:to PV)
rpsp=1.0 #ratio of jps:jpp (PV input to SOM:to PV)

sqrtK = sqrt(K)

je0 = 4.846*Vde*gEL*tau1/sqrtK
jp0 = 3.808*Vde*gIL*tau1/sqrtK
js0 = jp0*rps0
jee = 0.462*Vde*gEL*tau1/sqrtK
jpe = 1.270*Vde*gIL*tau1/sqrtK
jse = jpe*rpse
jep = 10.5*Vdi*gEL*tau1/sqrtK
jes = jep*reps
jpp = 9.5*Vdi*gIL*tau1/sqrtK
jps = jpp*rpsp
jsp = jps*rrps
jss = jpp*rips

fre=(je0*(jpp*jss - jps*jsp) - jp0*(jep*jss - jes*jsp) + js0*(jep*jps - jes*jpp))/(jee*jpp*jss - jee*jps*jsp - jep*jpe*jss + jep*jps*jse + jes*jpe*jsp - jes*jpp*jse)
frs=(je0*(jpe*jsp - jpp*jse) - jp0*(jee*jsp - jep*jse) + js0*(jee*jpp - jep*jpe))/(jee*jpp*jss - jee*jps*jsp - jep*jpe*jss + jep*jps*jse + jes*jpe*jsp - jes*jpp*jse)
frp=(-je0*(jpe*jss - jps*jse) + jp0*(jee*jss - jes*jse) - js0*(jee*jps - jes*jpe))/(jee*jpp*jss - jee*jps*jsp - jep*jpe*jss + jep*jps*jse + jes*jpe*jsp - jes*jpp*jse)

#connection probabilities
pe0 = K/N0
pp0 = pe0
ps0 = pp0*rps0
ppp = K/Np
pep = ppp
pes = pep*reps
pps = K/Ns
psp = pps*rrps
pss = ppp*rips
pee = K/Ne
ppe = pee
pse = ppe*rpse

J=[je0 jp0 js0 jee jpe jse jes jep jss jps jsp jpp]
tau = zeros(Ne+Np+Ns)

tau[1:Ne] = taue
tau[(1+Ne):(Ne+Np+Ns)] = taui

#weights = zeros(Ncells,Ncells)

#random connections
# weights[1:Ne,1:Ne] = (rand(Ne,Ne) .< pee)#*(jee/sqrtK) #excitatory --> excitatory
# weights[1:Ne,(1+Ne):(Ni+Ne)] = (rand(Ne,Ni) .< pei)#*(jei/sqrtK) #inhibitory --> excitatory
# weights[(1+Ne):(Ni+Ne),1:Ne] = (rand(Ni,Ne) .< pie)#*(jie/sqrtK) #excitatory --> inhibitory
# weights[(1+Ne):(Ni+Ne),(1+Ne):(Ni+Ne)] = (rand(Ni,Ni) .< pii)#*(jii/sqrtK) #inhibitory --> inhibitory
# weights[(1+Ne):(Ni+Ne),(1+Ni+Ne):Ncells] = (rand(Ni,N0) .< pi0)#*(ji0/sqrtK) #input layer --> inhibitory
# weights[1:Ne,(1+Ni+Ne):Ncells] = (rand(Ne,N0) .< pe0)#*(je0/sqrtK) #input layer --> excitatory

#Ks here are out-degrees, which is why they are not on average 1000 (does that make sense?)
#a neuron RECEIVES on average K inputs from each presyn population
#which would mean I neurons would have bigger out degrees since there are fewer of them
#this is definitely what we want because these are projections.
#indices sum(K[:i]) to sum(K[:i+1]) are downstream from neuron i+1

projections=[0]
Ks=zeros(Ncells)

#excitatory projections

for i=1:Ne
  weightsee = (rand(Ne,1) .< pee) #excitatory --> excitatory
  weightsee[i] = 0
  weightspe = (rand(Np,1) .< ppe) #excitatory --> PV
  weightsse = (rand(Ns,1) .< pse) #excitatory --> SOM
#   weightsee = weights[1:Ne,i]
#   weightsie = weights[(Ne+1):(Ne+Ni),i]
  fee=find(weightsee)
  fpe=find(weightspe)
  fse=find(weightsse)
  Ks[i]=size(fee,1)+size(fpe,1)+size(fse,1)
  projections=[projections; fee; fpe+Ne; fse+Ne+Np]
end

projections=projections[2:end] #make sure to get rid of the leading zero

#pv projections
for i=1:Np
  weightsep = (rand(Ne,1) .< pep) #PV --> excitatory
  weightspp = (rand(Np,1) .< ppp) #PV --> PV
  weightspp[i] = 0
  weightssp = (rand(Ns,1) .< psp) #PV --> SOM
#   weightsei = weights[1:Ne,Ne+i]
#   weightsii = weights[(Ne+1):(Ne+Ni),Ne+i]
  fep=find(weightsep)
  fpp=find(weightspp)
  fsp=find(weightssp)
  Ks[i+Ne]=size(fep,1)+size(fpp,1)+size(fsp,1)
  projections= [projections; fep; fpp+Ne; fsp+Ne+Np]
end

#som projections
for i=1:Ns
  weightses = (rand(Ne,1) .< pes) #SOM --> excitatory
  weightsps = (rand(Np,1) .< pps) #SOM --> PV
  weightsps[i] = 0
  weightsss = (rand(Ns,1) .< pss) #SOM --> SOM
#   weightsei = weights[1:Ne,Ne+i]
#   weightsii = weights[(Ne+1):(Ne+Ni),Ne+i]
  fes=find(weightses)
  fps=find(weightsps)
  fss=find(weightsss)
  Ks[i+Ne+Np]=size(fes,1)+size(fps,1)+size(fss,1)
  projections= [projections; fes; fps+Ne; fss+Ne+Np]
end

#input projections
for i=1:N0
  weightse0 = (rand(Ne,1) .< pe0) #input layer --> excitatory
  weightsp0 = (rand(Np,1) .< pp0) #input layer --> PV
  weightss0 = (rand(Ns,1) .< ps0) #input layer --> PV
#   weightse0 = weights[1:Ne,Ne+Ni+i]
#   weightsi0 = weights[(Ne+1):(Ne+Ni),Ne+Ni+i]
  fe0=find(weightse0)
  fp0=find(weightsp0)
  fs0=find(weightss0)
  Ks[i+Ne+Np+Ns]=size(fe0,1)+size(fp0,1)+size(fs0,1)
  projections= [projections; fe0; fp0+Ne; fs0+Ne+Np]
end

projections=readdlm("Documents/Piriform/code/Data/projections3.txt")
Ks=readdlm("Documents/Piriform/code/Data/Ks3.txt")
writedlm("Documents/Piriform/code/Data/J3.txt", J)
writedlm("Documents/Piriform/code/Data/tau3.txt", tau)
writedlm("Documents/Piriform/code/Data/J3.txt", J)


projections=readdlm("Documents/Piriform/projections3jp.txt")

Ks3jp=zeros(Ncells)
j=1
k=0
for i=1:Ncells
  while j<size(projections,1) && projections[j+1]>projections[j]
    j=j+1
    k=k+1
  end
  Ks3jp[i]=k+1
  k=0
  j=j+1
end

println(Ks3jp)
writedlm("Documents/Piriform/code/Data/Ks3jp.txt", Ks3jp)
# Js=readdlm("Documents/Piriform/code/Data/J3.txt")

# je0=Js[1]
# jp0=Js[2]
# js0=Js[3]
# jee=Js[4]
# jpe=Js[5]
# jse=Js[6]
# jes=Js[7]
# jep=Js[8]
# jss=Js[9]
# jps=Js[10]
# jsp=Js[11]
# jpp=Js[12]

# fre=(je0*(jpp*jss - jps*jsp) - jp0*(jep*jss - jes*jsp) + js0*(jep*jps - jes*jpp))/(jee*jpp*jss - jee*jps*jsp - jep*jpe*jss + jep*jps*jse + jes*jpe*jsp - jes*jpp*jse)
# frs=(je0*(jpe*jsp - jpp*jse) - jp0*(jee*jsp - jep*jse) + js0*(jee*jpp - jep*jpe))/(jee*jpp*jss - jee*jps*jsp - jep*jpe*jss + jep*jps*jse + jes*jpe*jsp - jes*jpp*jse)
# frp=(-je0*(jpe*jss - jps*jse) + jp0*(jee*jss - jes*jse) - js0*(jee*jps - jes*jpe))/(jee*jpp*jss - jee*jps*jsp - jep*jpe*jss + jep*jps*jse + jes*jpe*jsp - jes*jpp*jse)
