#Piriform cortex model
######Code written by Aubrey Thompson

Intended for internal use in the Doiron lab.

#The model
Adapted from Pehlevan and Sompolinsky (PLOS One 2014) [1] model of piriform cortex, an LIF network with E, I, and external populations and random connections.

This model splits the inhibitory population into SOM and PV, and all populations up into “left” and “right” (or rostral and caudal). The idea is that one side has more SOM neurons. We would like for that side to have a higher excitatory firing rate (due to a dis-inhibitory circuit).

The model parameters are all the same as [1] except that one can modify the relative strength of connection probabilities. There are 4 key ratios between connection probabilities:

- rpse: ratio of excitatory connection to SOM:to PV
- reps: ratio of excitatory connection from SOM:from PV
- rips: ratio of interconnection within SOM:within PV
- rrps: ratio of reciprocal connections between i pops from PV:from SOM

There is one other ratio which I allow to be varied but probably makes the most sense as 1:

- rps0=ratio of input layer to SOM:to PV

#The code
Adapted from Ashok Litwin-Kumar’s cluster model (Nature Neuro 2012) code.

I recommend using the cluster for creating parameters and running the simulation, then JUNO (Julia’s IDE) to analyze the data. You could probably also use jupyter for this if you wish.

makew.pbs: bash script which creates jobs on the cluster to make a network, run a simulation with it, and save data. you can use job arrays to vary parameters (code to handle that is included in makew6.jl and runsompsim6.jl)

makew6.jl: set network parameters (size, connection strength and probability, neural parameters) and create:
- projections.txt: a stacked vector of all the neuron’s projections. indices sum(K[:i]) to sum(K[:i+1]) of projections are downstream neurons from neuron i+1. 
- Ks.txt: Ncells x 1 vector of each neuron’s outdegree
- J.txt: 12 x 1 vector of connection strengths 
- [je0 jp0 js0 jee jpe jse jes jep jss jps jsp jpp]
- tau.txt: (Ne+Ni) x 1 vector of each neuron’s time constant

sompsim6.jl: simulate a network with imported parameters made from makew6.jl. Input:
- stim: the ID number of pre-made stimulus (each stimulus is a vector of N0 firing rates). 
- n: the ID number of Ns, for finding firing rate convergences over N

Output:
- times: spike times, one row per neuron, padded with zeros at the end
- spikes: number of spikes per neuron (Ncells x 1 vector)
- stimrates: firing rate of each neuron in the input layer (N0 x 1 vector)

runsompsim6.jl: run sompsim (with varied parameters) and save the data.

plotsompsim.jl: meant to be run on your machine in chunks. There are sections to plot rasters, find ISI CVs, and get firing rate histograms of one run; plot firing rates over N for multiple runs with varying N, find parameters which yield a higher E firing rate on the left using multiple runs with different ratio parameters; and compute selectivity and sparsity using multiple runs with differing stimuli.
