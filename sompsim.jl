#Copyright (C) 2014 Ashok Litwin-Kumar
#modified by Aubrey Thompson 2016 to simulate
#Pehlevan, Sompolinsky: Selectivity and Sparseness in Randomly Connected Balanced Networks
#see README for more information
using Distributions

function sompsim(stim)
	println("importing shared parameters")
  projections=readdlm("Documents/Piriform/code/Data/projections.txt")
  J=readdlm("Documents/Piriform/code/Data/J.txt")
  Ks=readdlm("Documents/Piriform/code/Data/Ks.txt")
  tau=readdlm("Documents/Piriform/code/Data/tau.txt")
  Nns=readdlm("Documents/Piriform/code/Data/Nns.txt")
  Ne=convert(Int64,Nns[1])
  Ni=convert(Int64,Nns[2])
  N0=convert(Int64,Nns[3])
  Ncells=Ne+Ni+N0

  println("setting up sim parameters")
  #if you change gEL, you also need to change it in makew!
  gEL=0.05
	T = 100 #simulation time (ms)

  taua = 100
  tau1 = 3
  tau2 = 1

  VL = -65 #(mV) leak potential
  Vth = -50 #spiking threshold

	dt = .05 #simulation timestep (ms)
	refrac = 5 #refractory period (ms)

  maxrate = 400 #(Hz) maximum average firing rate.  if the average firing rate across the simulation for any neuron exceeds this value, some of that neuron's spikes will not be saved

  #input statistics
  #input neurons fire with Poisson spike trains with rates given by an exponential distribution
  #each stimulus generates a random input firing rate pattern
  #the mean population rate of all the stimuli is 10.16
  ud=Uniform()
  if stim==0
    r = 10.16 #mean firing rate (Hz)
    p = 0.5
    stimrates = zeros(N0)
    for si = 1:N0
      pr = rand(ud)
      er = 200
      if pr<p
        while er>150
          sid = Exponential(r/p)
          er = rand(sid)
        end
        stimrates[si] = er
      end
    end
  else
    stimrates=readdlm(string("stimrates",stim,".txt"))
  end

	maxTimes = round(Int,maxrate*T/1000)
	times = zeros(Ncells,maxTimes)
	ns = zeros(Int,Ncells)

  #forwardInputsX are inputs FROM X TO neuron at index
	forwardInputsE = zeros(Ne+Ni) #summed weight of incoming E spikes
	forwardInputsI = zeros(Ne+Ni)
  forwardInputs0 = zeros(Ne+Ni)
	forwardInputsEPrev = zeros(Ne+Ni) #as above, for previous timestep
	forwardInputsIPrev = zeros(Ne+Ni)
  forwardInputs0Prev = zeros(Ne+Ni)

	xerise = zeros(Ncells) #auxiliary variables for E/I currents (difference of exponentials)
	xedecay = zeros(Ncells)
	xirise = zeros(Ncells)
	xidecay = zeros(Ncells)
  x0rise = zeros(Ncells)
	x0decay = zeros(Ncells)

  adaptInput = zeros(Ne+Ni)

	v = rand(Ne+Ni)*(Vth-VL)+VL #membrane voltage

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
    forwardInputs0[:] = 0

		for ci = 1:(Ne+Ni)
			xerise[ci] += -dt*xerise[ci]/tau2 + forwardInputsEPrev[ci]
			xedecay[ci] += -dt*xedecay[ci]/tau1 + forwardInputsEPrev[ci]
			xirise[ci] += -dt*xirise[ci]/tau2 + forwardInputsIPrev[ci]
			xidecay[ci] += -dt*xidecay[ci]/tau1 + forwardInputsIPrev[ci]
      x0rise[ci] += -dt*x0rise[ci]/tau2 + forwardInputs0Prev[ci]
			x0decay[ci] += -dt*x0decay[ci]/tau1 + forwardInputs0Prev[ci]

			synInput = (xedecay[ci] - xerise[ci])/(tau1 - tau2) + (xidecay[ci] - xirise[ci])/(tau1 - tau2) + (x0decay[ci] - x0rise[ci])/(tau1 - tau2)
      if ci<Ne
        adaptInput[ci] += -dt*adaptInput[ci]/taua
      end

			if t > (lastSpike[ci] + refrac)  #not in refractory period
				v[ci] += dt*(VL-v[ci])/tau[ci] + synInput - adaptInput[ci]
				if v[ci] > Vth  #spike occurred
          if ci<Ne+1
            adaptInput[ci] = adaptInput[ci] + 0.1*(Vth-VL)*gEL
          end
					v[ci] = VL
					lastSpike[ci] = t
					if ns[ci] < maxTimes
            ns[ci] = ns[ci]+1
						times[ci,ns[ci]] = t
					end
          projci=projections[(sum(Ks[1:(ci-1)])+1):(sum(Ks[1:ci]))]
					for j =1:size(projci,1)
            if projci[j]<=Ne
              if ci <= Ne  #E->E synapse
                forwardInputsE[projci[j]] += J[3]
              elseif ci <= Ne+Ni && ci>Ne #I->E synapse
                forwardInputsI[projci[j]] += J[5]
              end
            elseif projci[j]<=Ne+Ni && projci[j]>Ne
              if ci <= Ne  #E->I synapse
							  forwardInputsE[projci[j]] += J[4]
						  elseif ci <= Ne+Ni && ci>Ne #I->I synapse
							  forwardInputsI[projci[j]] += J[6]
              end #end if(excitatory or inhibitory pre-syn cell)
            end	#end if(excitatory or inhibitory post-syn cell)
					end #end loop over synaptic projections
				end #end if(spike occurred)
			end #end if(not refractory)
		end #end loop over neurons
    for ci = 1:N0
      if stimrates[ci]>0
        F = stimrates[ci]*dt/1000 #probability that there is a spike in this interval
        pr2 = rand(ud)
        if F>pr2
          lastSpike[ci+Ni+Ne] = t
          if ns[ci+Ni+Ne] < maxTimes
            ns[ci+Ni+Ne] = ns[ci+Ni+Ne]+1
					  times[ci+Ni+Ne,ns[ci+Ni+Ne]] = t
				  end
          projci=projections[(sum(Ks[1:(ci+Ni+Ne-1)])+1):(sum(Ks[1:(ci+Ni+Ne)]))]
          for j =1:size(projci,1)
            if projci[j]<=Ne
              forwardInputs0[projci[j]] += J[1]
            elseif projci[j]<=Ne+Ni && projci[j]>Ne
              forwardInputs0[projci[j]] += J[2]
            end
					end
        end
      end
    end
		forwardInputsEPrev = copy(forwardInputsE)
		forwardInputsIPrev = copy(forwardInputsI)
    forwardInputs0Prev = copy(forwardInputs0)
	end #end loop over time
	@printf("\r")
	times = times[:,1:maximum(ns)]

	return times,ns,Ne,Ni,Ncells,T,stimrates
end

