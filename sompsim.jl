#this file is part of alk_cluster_2012
#Copyright (C) 2014 Ashok Litwin-Kumar
#see README for more information
using Distributions

function sompsim()
	println("setting up parameters")
	Ncells = 3000
	Ne = 2000
	Ni = 500
  N0 = Ncells-Ni-Ne
  K = 1000.0 #average number of E->E connections per neuron
	T = 20 #simulation time (ms)

  CEM = 1.0 #(μF/cm^2)
  CIM = 1.0

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

  je0 = 4.846*Vde*gEL*tau1
  ji0 = 3.808*Vde*gIL*tau1
  jee = 0.462*Vde*gEL*tau1
  jie = 1.270*Vde*gIL*tau1
  jei = 10.5*Vdi*gEL*tau1
  jii = 9.5*Vdi*gIL*tau1

# 	jie = 4./(taui*sqrtK)
# 	jei = -16.*1.2/(taue*sqrtK)
# 	jii = -16./(taui*sqrtK)

# 	#set up connection probabilities within and without blocks
# 	ratioejee = 1.9
# 	jeeout = 10./(taue*sqrtK)
# 	jeein = ratioejee*10./(taue*sqrtK)

# 	ratiopee = 2.5
# 	peeout = K/(Nepop*(ratiopee-1) + Ne)
# 	peein = ratiopee*peeout

	#stimulation
# 	Nstim = 400 #number of neurons to stimulate (indices 1 to Nstim will be stimulated)
# 	stimstr = .07/taue
# 	stimstart = T-500
# 	stimend = T

	#constant bias to each neuron type
# 	muemin = 1.1
# 	muemax = 1.2
# 	muimin = 1
# 	muimax = 1.05

	vre = 0. #reset voltage

	threshe = 1 #threshold for exc. neurons
	threshi = 1

	dt = .05 #simulation timestep (ms)
	refrac = 5 #refractory period (ms)

	#synaptic time constants (ms)
# 	tauerise = 1
# 	tauedecay = 3
# 	tauirise = 1
# 	tauidecay = 2

	maxrate = 50 #(Hz) maximum average firing rate.  if the average firing rate across the simulation for any neuron exceeds this value, some of that neuron's spikes will not be saved

  #input statistics
  #input neurons fire with Poisson spike trains with rates given by an exponential distribution
  #each stimulus generates a random input firing rate pattern
  #(each stim is given to the input layer? in a way the input layer represents the stim?)
  println("getting input stats")
  r = 10.16 #mean firing rate (Hz)
  stimrates = zeros(N0)
  ud = Uniform()
  er = 200
  for si = 1:N0
    p = rand(ud)
    if p<.5
      while er>150
        sid = Exponential(1.0/r)
        er = rand(sid)
        print(er)
      end
      stimrates[si] = er
    end
  end

# 	mu = zeros(Ncells)
 	thresh = zeros(Ncells)
 	tau = zeros(Ncells)
  gL = zeros(Ncells)

# 	mu[1:Ne] = (muemax-muemin)*rand(Ne) + muemin
# 	mu[(Ne+1):Ncells] = (muimax-muimin)*rand(Ni) + muimin

	thresh[1:Ne] = Vth
	thresh[(1+Ne):Ncells] = Vth

	tau[1:Ne] = taue
	tau[(1+Ne):Ncells] = taui

  gL[1:Ne] = gEL
	gL[(1+Ne):(Ne+Ni)] = gIL

	#Npop = round(Int,Ne/Nepop)

	weights = zeros(Ncells,Ncells)

	#random connections
	weights[1:Ne,1:Ne] = (jee/sqrtK)*(rand(Ne,Ne) .< pee) #excitatory --> excitatory
	weights[1:Ne,(1+Ne):(Ni+Ne)] = (jei/sqrtK)*(rand(Ne,Ni) .< pei) #inhibitory --> excitatory
	weights[(1+Ne):(Ni+Ne),1:Ne] = (jie/sqrtK)*(rand(Ni,Ne) .< pie) #excitatory --> inhibitory
	weights[(1+Ne):(Ni+Ne),(1+Ne):(Ni+Ne)] = (jii/sqrtK)*(rand(Ni,Ni) .< pii) #inhibitory --> inhibitory
  weights[(1+Ne):(Ni+Ne),(1+Ni+Ne):Ncells] = (ji0/sqrtK)*(rand(Ni,N0) .< pi0) #input layer --> inhibitory
  weights[1:Ne,(1+Ni+Ne):Ncells] = (je0/sqrtK)*(rand(Ne,N0) .< pe0) #input layer --> excitatory


	#connections within cluster
# 	for pii = 1:Npop
# 		ipopstart = 1 + Nepop*(pii-1)
# 		ipopend = pii*Nepop

# 		weights[ipopstart:ipopend,ipopstart:ipopend] = jeein*(rand(Nepop,Nepop) .< peein)
# 	end

	for ci = 1:Ncells
		weights[ci,ci] = 0
	end

	maxTimes = round(Int,maxrate*T/1000)
	times = zeros(Ncells,maxTimes)
	ns = zeros(Int,Ncells)

  #TODO: think harder about whether the input layer needs these variables
  #or whether they need to be included in that of the other pops
  #i.e. are these variables Inputs TO 'x' or FROM 'x'

	forwardInputsE = zeros(Ne+Ni) #summed weight of incoming E spikes
	forwardInputsI = zeros(Ne+Ni)
  #forwardInputs0 = zeros(Ncells)
	forwardInputsEPrev = zeros(Ne+Ni) #as above, for previous timestep
	forwardInputsIPrev = zeros(Ne+Ni)
  #forwardInputs0Prev = zeros(Ncells)

	xerise = zeros(Ncells) #auxiliary variables for E/I currents (difference of exponentials)
	xedecay = zeros(Ncells)
	xirise = zeros(Ncells)
	xidecay = zeros(Ncells)
  x0rise = zeros(Ncells)
	x0decay = zeros(Ncells)

  adaptInput = zeros(Ne+Ni)

	v = rand(Ncells) #membrane voltage

	lastSpike = -100*ones(Ncells) #time of last spike

	Nsteps = round(Int,T/dt)

	println("starting simulation")

	#begin main simulation loop
	for ti = 1:Nsteps
		if mod(ti,Nsteps/100) == 1  #print percent complete
			@printf("\r%d%%",round(Int,100*ti/Nsteps))
		end
		t = dt*ti
		forwardInputsE[:] = 0
		forwardInputsI[:] = 0
    #forwardInputs0[:] = 0

		for ci = 1:(Ne+Ni)
			xerise[ci] += -dt*xerise[ci]/tau2 + forwardInputsEPrev[ci]
			xedecay[ci] += -dt*xedecay[ci]/tau1 + forwardInputsEPrev[ci]
			xirise[ci] += -dt*xirise[ci]/tau2 + forwardInputsIPrev[ci]
			xidecay[ci] += -dt*xidecay[ci]/tau1 + forwardInputsIPrev[ci]
      x0rise[ci] += -dt*x0rise[ci]/tau2 #+ forwardInputs0Prev[ci]
			x0decay[ci] += -dt*x0decay[ci]/tau1 #+ forwardInputs0Prev[ci]

			synInput = (xedecay[ci] - xerise[ci])/(tau1 - tau2) + (xidecay[ci] - xirise[ci])/(tau1 - tau2) + (x0decay[ci] - x0rise[ci])/(tau1 - tau2)
      if ci<Ne
        adaptInput[ci] = adaptInput[ci]-dt*adaptInput[ci]/taua
      end
# 			if (ci < Nstim) && (t > stimstart) && (t < stimend)
# 				synInput += stimstr;
# 			end
			if t > (lastSpike[ci] + refrac)  #not in refractory period
				v[ci] += dt*((1/tau[ci])*gL[ci]*(VL-v[ci]) + synInput - adaptInput[ci])
				if v[ci] > thresh[ci]  #spike occurred
          if ci<Ne+1
            adaptInput[ci] = adaptInput[ci] + 0.1*(Vth-VL)*gEL
          end
					v[ci] = vre
					lastSpike[ci] = t
					if ns[ci] < maxTimes
            ns[ci] = ns[ci]+1
						times[ci,ns[ci]] = t
					end

					for j = 1:(Ne+Ni)
						if weights[j,ci] > 0  #E synapse
							forwardInputsE[j] += weights[j,ci]
						elseif weights[j,ci] < 0  #I synapse
							forwardInputsI[j] += weights[j,ci]
						end
					end #end loop over synaptic projections
				end #end if(spike occurred)
			end #end if(not refractory)
		end #end loop over neurons
    for ci = 1:N0
      if stimrates[ci]>0
        F = 1-exp(-stimrates[ci]*dt)
        p = rand(ud)
        if F>p
          lastSpike[ci+Ni+Ne] = t
          if ns[ci+Ni+Ne] < maxTimes
            ns[ci+Ni+Ne] = ns[ci+Ni+Ne]+1
					  times[ci+Ni+Ne,ns[ci+Ni+Ne]] = t
				  end
          for j = 1:(Ne+Ni)
            forwardInputsE[j] += weights[j,ci+Ni+Ne]
					end
        end
      end
    end
		forwardInputsEPrev = copy(forwardInputsE)
		forwardInputsIPrev = copy(forwardInputsI)
    #forwardInputs0Prev = copy(forwardInputs0)
	end #end loop over time
	@printf("\r")
	times = times[:,1:maximum(ns)]

	return times,ns,Ne,Ni,Ncells,T
end
