##Uncomment below to get the job ID
##Use this to run the same jobs with one varying parameter

# function getenv(var::ASCIIString)
#   val = ccall((:getenv, "libc"),
#               Ptr{Int8}, (Ptr{Int8},), var)
#   if val == C_NULL
#     error("getenv: undefined variable: ", var)
#   end
#   bytestring(val)
# end
# ni = getenv("PBS_ARRAYID")
# #allNs = readdlm("allNs.txt")
# n=convert(Int64,parse(ni))
# println(n)
# valids=readdlm("params2.txt")

##Here, params.txt is an array of ratios of probability connections
##You could also loop over N:
#allNs = readdlm("Nns.txt")
#N=allNs[n]

N=45000 #default
Nel = convert(Int64,round(N*2.0/9))
Ner = convert(Int64,round(N*2.0/9))
Ne=Nel+Ner
Npl = convert(Int64,round(N*1.0/36))
Npr = convert(Int64,round(N*1.0/36))
Np=Npl+Npr
Nsl = convert(Int64,round(N*1.0/27))
Nsr = convert(Int64,round(N*1.0/54))
Ns=Nsl+Nsr
Ni=Np+Ns
N0 = convert(Int64,round(N*4.0/9))
K = convert(Int64,round(N*1.0/45)) #average number of E->E connections per neuron

Nns=[Nel,Ner,Npl,Npr,Nsl,Nsr,N0]
Ncells = sum(Nns)
Nscells = Ncells-N0
writedlm("Nns6-asym.txt",Nns)

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

params=readdlm("params2.txt")
println(size(params))
cin=convert(Int64,valids[n])
println(cin)
runparams=params[cin,:]

rpse=runparams[1] #ratio of jse:jpe (excitatory connection to SOM:to PV)
reps=runparams[2] #ratio of jes:jep (excitatory connection from SOM:from PV)
rips=runparams[3] #ratio of jss:jpp (interconnection within SOM:within PV)
rrps=runparams[4] #ratio of jsp:jps (reciprocal connections between i pops from PV:from SOM)
rps0=runparams[5] #ratio of js0:jp0 (input layer to SOM:to PV)
delta=0

sqrtK = sqrt(K)

je0 = 4.846*Vde*gEL*tau1/sqrtK
jp0 = 3.808*Vde*gIL*tau1/sqrtK
js0 = jp0
jee = 0.462*Vde*gEL*tau1/sqrtK
jpe = 1.270*Vde*gIL*tau1/sqrtK
jse = jpe
jep = 10.5*Vdi*gEL*tau1/sqrtK
jes = jep
jpp = 9.5*Vdi*gIL*tau1/sqrtK
jps = jpp
jsp = jps
jss = jpp

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
tau = zeros(Nscells)

tau[1:Ne] = taue
tau[(1+Ne):Nscells] = taui

#Ks here are out-degrees
#indices sum(K[:i]) to sum(K[:i+1]) of projections are downstream neurons from neuron i+1

projections=[0]
Ks=zeros(Ncells)

#excitatory projections
for i=1:Ne
  weightsee = (rand(Ne,1) .< pee) #excitatory (L&R) --> excitatory (L&R)
  weightsee[i] = 0
  weightspe = (rand(Np,1) .< ppe) #excitatory (L&R) --> PV (L&R)
  weightsse = (rand(Ns,1) .< pse) #excitatory (L&R)--> SOM (L&R)
  fee=find(weightsee)
  fpe=find(weightspe)
  fse=find(weightsse)
  Ks[i]=size(fee,1)+size(fpe,1)+size(fse,1)
  projections=[projections; fee; fpe+Ne; fse+Ne+Np]
end

projections=projections[2:end] #make sure to get rid of the leading zero

#pvL projections
for i=1:Npl
  weightselpl = (rand(Nel,1) .< pep+delta/2) #PVL --> EL
  weightsplpl = (rand(Npl,1) .< ppp+delta/2) #PVL --> PVL
  weightsplpl[i] = 0
  weightsslpl = (rand(Nsl,1) .< psp+delta/2) #PVL --> SOML
  fep=find(weightselpl)
  fpp=find(weightsplpl)
  fsp=find(weightsslpl)
  Ks[i+Ne]=size(fep,1)+size(fpp,1)+size(fsp,1)
  projections= [projections; fep; fpp+Ne; fsp+Ne+Np]
end

#pvR projections
for i=1:Npr
  weightserpr = (rand(Ner,1) .< pep-delta/2) #PVR --> ER
  weightsprpr = (rand(Npr,1) .< ppp-delta/2) #PVR --> PVR
  weightsprpr[i] = 0
  weightssrpr = (rand(Nsr,1) .< psp-delta/2) #PVR --> SOMR
  fep=find(weightserpr)
  fpp=find(weightsprpr)
  fsp=find(weightssrpr)
  Ks[i+Ne+Npl]=size(fep,1)+size(fpp,1)+size(fsp,1)
  projections= [projections; fep+Nel; fpp+Ne+Npl; fsp+Ne+Np+Nsl]
end

#somL projections
for i=1:Nsl
  weightselsl = (rand(Nel,1) .< pes+delta/2) #SOML --> EL
  weightsplsl = (rand(Npl,1) .< pps+delta/2) #SOML --> PVL
  weightsslsl = (rand(Nsl,1) .< pss+delta/2) #SOML --> SOML
  weightsslsl[i] = 0
  fes=find(weightselsl)
  fps=find(weightsplsl)
  fss=find(weightsslsl)
  Ks[i+Ne+Np]=size(fes,1)+size(fps,1)+size(fss,1)
  projections= [projections; fes; fps+Ne; fss+Ne+Np]
end

#somR projections
for i=1:Nsr
  weightsersr = (rand(Ner,1) .< pes-delta/2) #SOMR --> ER
  weightsprsr = (rand(Npr,1) .< pps-delta/2) #SOMR --> PVR
  weightssrsr = (rand(Nsr,1) .< pss-delta/2) #SOMR --> SOMR
  weightssrsr[i] = 0
  fes=find(weightsersr)
  fps=find(weightsprsr)
  fss=find(weightssrsr)
  Ks[i+Ne+Np+Nsl]=size(fes,1)+size(fps,1)+size(fss,1)
  projections= [projections; fes+Nel; fps+Ne+Npl; fss+Ne+Np+Nsl]
end

#input projections
for i=1:N0
  weightse0 = (rand(Ne,1) .< pe0) #input layer --> excitatory
  weightsp0 = (rand(Np,1) .< pp0) #input layer --> PV
  weightss0 = (rand(Ns,1) .< ps0) #input layer --> PV
  fe0=find(weightse0)
  fp0=find(weightsp0)
  fs0=find(weightss0)
  Ks[i+Ne+Np+Ns]=size(fe0,1)+size(fp0,1)+size(fs0,1)
  projections= [projections; fe0; fp0+Ne; fs0+Ne+Np]
end

#these get read and used in sompsim
writedlm(string("tau.txt"), tau)
writedlm(string("J.txt"), J)
writedlm(string("projections.txt"),projections)
writedlm(string("Ks.txt"), Ks)
