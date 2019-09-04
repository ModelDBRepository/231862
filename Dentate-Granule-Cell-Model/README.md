[![GitHub tag](https://img.shields.io/github/tag/MarcelBeining/Dentate-Granule-Cell-Model.svg?style=for-the-badge)](https://github.com/MarcelBeining/Dentate-Granule-Cell-Model/releases)
![GitHub top language](https://img.shields.io/github/languages/top/MarcelBeining/Dentate-Granule-Cell-Model.svg?style=for-the-badge)
[![GitHub](https://img.shields.io/github/license/MarcelBeining/Dentate-Granule-Cell-Model.svg?style=for-the-badge)](https://github.com/MarcelBeining/Dentate-Granule-Cell-Model/blob/master/LICENSE)
[![GitHub contributors](https://img.shields.io/github/contributors/MarcelBeining/Dentate-Granule-Cell-Model.svg?style=for-the-badge)](https://github.com/MarcelBeining/Dentate-Granule-Cell-Model/graphs/contributors)
![GitHub repo size in bytes](https://img.shields.io/github/repo-size/MarcelBeining/Dentate-Granule-Cell-Model.svg?style=for-the-badge)
[![GitHub issues](https://img.shields.io/github/issues/MarcelBeining/Dentate-Granule-Cell-Model.svg?style=for-the-badge)](https://github.com/MarcelBeining/Dentate-Granule-Cell-Model/issues)

# Dentate granule cell model
This GC compartmental model accompanies the manuscript 
[Beining et al (2017): T2N as a new tool for robust electrophysiological modeling demonstrated for mature and adult-born dentate granule cells. eLife](https://elifesciences.org/articles/26517)

# Getting started
Note! Be sure to add the Modules subfolder to your Matlab path, as this folder contains the [T2N](https://github.com/MarcelBeining/T2N) and [treestoolbox](http://www.treestoolbox.org/) packages!
Matlab2009b or higher is required, some functionality in the model, requires Matlab2016b or higher!

T2N automatically compiles the .mod files in the lib_mech folder necessary for the model to be working. 

## Important functions/scripts:
GC_experiments.m	is the main script that initializes the model, runs the experiments and makes figures

## Scripts for understanding:
GC_exp_tests.m		is similar to "GC_experiments.m" but can be better used for playing with the model
GC_initModel.m		initializes the model after the parameters were defined in the first section of GC_experiments.m or GC_exp_tests.m
GC_biophys.m 		comprises the definition of all the channels and channel densities
GC_spinedensity.m	comprises the definition of spine densitiy scaling factors

*Scripts/functions that make or plot experiments start with "aGC_" and are all incorporated in "GC_experiments.m" or at least "GC_exp_tests.m"*

*Notice that you have to recreate the dll file once you change something in an mod file.*

# Licensing
This software is published under the MIT license. For further information read the LICENSE file.
