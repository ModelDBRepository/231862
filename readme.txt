This GC compartmental model accompanies the manuscript 
Beining et al (2017): T2N as a new tool for robust electrophysiological modeling demonstrated for mature and adult-born dentate granule cells. eLife

It runs in NEURON but is written in Matlab with T2N, a new Matlab-NEURON interface. The zip file contains the T2N and the TREES toolbox package and all Matlab scripts and functions that define the model and the simulations that were run to produce the figures of the paper.

Installation: Unzip the file to the desired location, then start Matlab and go to this location. Run "runthisAfterUnzip.m". All folders and files are added then to your Matlab search path.

Note: mod files are automatically compiled by T2N.

Important functions/scripts:
GC_experiments.m		in the Active GC model folder is the main script that initializes the model, runs the experiments and makes figures
t2n_Tutorial.mlx  		in the T2N/Tutorials folder is a Matlab live script that helps to get into T2N. If your Matlab version does not support live scripts(older than Matlab2016a), there is a .m script doing the same simulations with text written as comments
Documentation T2N.docx	in the T2N folder is the documentation to T2N explaining all functions/features of T2N


Scripts for understanding the model:
GC_exp_tests.m		is similar to "GC_experiments.m" but can be better used for playing with the model
GC_initModel.m		initializes the model after the parameters were defined in the first section of GC_experiments.m
GC_biophys.m 		comprises the definition of all the channels and channel densities
GC_spinedensity.m	comprises the definition of spine densitiy scaling factors


For further questions, contact me via beining@fias.uni-frankfurt.de

