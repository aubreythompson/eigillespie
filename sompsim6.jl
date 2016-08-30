#this file is part of alk_cluster_2012
#Copyright (C) 2014 Ashok Litwin-Kumar
#see README for more information
using Distributions

function sompsim(stim,n)
  println("importing shared parameters")
  projections=readdlm(string("projections",n,"-2.txt"))
  J=readdlm(string("J",n,"-2.txt"))
  Ks=readdlm(string("Ks",n,"-2.txt"))
  tau=readdlm(string("tau",n,"-2.txt"))
  Nns=readdlm("Nns6doubledensesmall.txt")
  println(J)
  println(Nns)
  Nel=convert(Int64,Nns[1])
  Ner=convert(Int64,Nns[2])
  Ne=Nel+Ner
  Npl=convert(Int64,Nns[3])
  Npr=convert(Int64,Nns[4])
  Np=Npl+Npr
  Nsl=convert(Int64,Nns[5])
  Nsr=convert(Int64,Nns[6])
  Ns=Nsl+Nsr
  Ni=Np+Ns
  N0=convert(Int64,Nns[7])
  Ncells=Ne+Ni+N0

  println("setting up sim parameters")
  #if you change gEL, you also need to change it in makew!
  gEL=0.05
  T = 2000 #simulation time (ms)

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
  #(each stim is given to the input layer? in a way the input layer represents the stim?)

  #the mean population rate of all the stimuli is 10.16
  #but also, population rate=avg(neuron's rates)
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
  #oh, that means forwardInputs0 should be a thing? I guess I'm crazy??
	forwardInputsE = zeros(Ne+Ni) #summed weight of incoming E spikes
	forwardInputsP = zeros(Ne+Ni)
	forwardInputsS = zeros(Ne+Ni)
 	forwardInputs0 = zeros(Ne+Ni)
	forwardInputsEPrev = zeros(Ne+Ni) #as above, for previous timestep
	forwardInputsPPrev = zeros(Ne+Ni)
	forwardInputsSPrev = zeros(Ne+Ni)
	forwardInputs0Prev = zeros(Ne+Ni)

	xerise = zeros(Ne+Ni) #auxiliary variables for E/I currents (difference of exponentials)
	xedecay = zeros(Ne+Ni)
	xprise = zeros(Ne+Ni)
	xpdecay = zeros(Ne+Ni)
	xsrise = zeros(Ne+Ni)
	xsdecay = zeros(Ne+Ni)
  x0rise = zeros(Ne+Ni)
	x0decay = zeros(Ne+Ni)

  adaptInput = zeros(Ne+Ni)

	v = rand(ud,Ne+Ni)*(Vth-VL)+VL #membrane voltage
	initialv=copy(v)
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
		forwardInputsP[:] = 0
		forwardInputsS[:] = 0
   	forwardInputs0[:] = 0

		for ci = 1:(Ne+Ni)
			xerise[ci] += -dt*xerise[ci]/tau2 + forwardInputsEPrev[ci]
			xedecay[ci] += -dt*xedecay[ci]/tau1 + forwardInputsEPrev[ci]
			xprise[ci] += -dt*xprise[ci]/tau2 + forwardInputsPPrev[ci]
			xpdecay[ci] += -dt*xpdecay[ci]/tau1 + forwardInputsPPrev[ci]
			xsrise[ci] += -dt*xsrise[ci]/tau2 + forwardInputsSPrev[ci]
			xsdecay[ci] += -dt*xsdecay[ci]/tau1 + forwardInputsSPrev[ci]
      x0rise[ci] += -dt*x0rise[ci]/tau2 + forwardInputs0Prev[ci]
			x0decay[ci] += -dt*x0decay[ci]/tau1 + forwardInputs0Prev[ci]

			synInput = (xedecay[ci] - xerise[ci])/(tau1 - tau2) + (xpdecay[ci] - xprise[ci])/(tau1 - tau2) + (xsdecay[ci] - xsrise[ci])/(tau1-tau2) + (x0decay[ci] - x0rise[ci])/(tau1 - tau2)
      if ci<Ne+1
    		adaptInput[ci] += -dt*adaptInput[ci]/taua
    	end
# 			if (ci < Nstim) && (t > stimstart) && (t < stimend)
# 				synInput += stimstr;
# 			end
			if t > (lastSpike[ci] + refrac) && ti>1 #not in refractory period
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
                forwardInputsE[projci[j]] += J[4]
    					elseif ci <=Ne+Np && ci>Ne #PV->E synapse
                forwardInputsP[projci[j]] += J[8]
              elseif ci<=(Ne+Np+Ns) && ci>(Ne+Np) #SOM->E
						    forwardInputsS[projci[j]] += J[7]
					    end
            elseif projci[j]<=Ne+Np && projci[j]>Ne
              if ci <= Ne  #E->PV synapse
                  forwardInputsE[projci[j]] += J[5]
              elseif ci <= Ne+Np && ci>Ne  #PV->PV synapse
                forwardInputsP[projci[j]] += J[12]
              elseif ci <= Ne+Np+Ns && ci>Ne+Np #S->PV
                forwardInputsS[projci[j]] += J[10]
              end
						elseif projci[j]<=Ne+Np+Ns && projci[j]>Ne+Np
							if ci <= Ne #E->SOM
								forwardInputsE[projci[j]] += J[6]
							elseif ci <=Ne+Np && ci>Ne #PV->SOM
								forwardInputsP[projci[j]] += J[11]
							elseif ci<=Ne+Np+Ns && ci>Ne+Np # SOM->SOM
								forwardInputsS[projci[j]] += J[9]
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
            elseif projci[j]<=Ne+Np && projci[j]>Ne
              forwardInputs0[projci[j]] += J[2]
            elseif projci[j]<=Ne+Np+Ns && projci[j]>Ne+Np
			          forwardInputs0[projci[j]] += J[3]
						end
					end #end projection loop
  			end #end if(input neuron spiked)
      end #end if(stim rate nonzero)
		end #end loop over neurons
		forwardInputsEPrev = copy(forwardInputsE)
   	forwardInputsPPrev = copy(forwardInputsP)
    forwardInputsSPrev = copy(forwardInputsS)
 		forwardInputs0Prev = copy(forwardInputs0)
  end #end loop over time
	@printf("\r")
  times = times[:,1:maximum(ns)]
  return times,ns,Ne,Np,Ns,Ncells,T,stimrates
end
