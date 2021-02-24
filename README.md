# FAUST simulation study

This repository contains scripts to reproduce the
FAUST simulation study.

To reproduce the FAUST analysis of the FlowCAP IV data, run the scripts in numeric order.

In `1_num_of_clusters_sim.R`, the script must be run twice: once with the boolean variable `gaussianSwitch` on line 40 set to FALSE, and once with it set to TRUE.

In `2_cvauc_simulation.R`, the script must be run twice: once with the boolean variable `gaussianSwitch` on line 34 set to FALSE, and once with it set to TRUE.

Libraries needed to run the simulation must be installed prior to running these scripts.

