# FAUST simulation study

This repository contains scripts to reproduce the
FAUST simulation study.

To reproduce the FAUST simulations, run the scripts in numeric order.

In `1_num_of_clusters_sim.R`, the script must be run twice: once with the boolean variable `gaussianSwitch` on line 40 set to FALSE, and once with it set to TRUE.

In `2_cvauc_simulation.R`, the script must be run twice: once with the boolean variable `gaussianSwitch` on line 34 set to FALSE, and once with it set to TRUE.

In `3_mk_figure.R`, the path on line 11 must be udpdated to point to the directory housing the simulation results.

Libraries needed to run the simulation are listed at the head of each numeric file. These must be installed prior to running the simulation scripts.

