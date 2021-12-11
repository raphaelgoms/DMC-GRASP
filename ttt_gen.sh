#!/bin/bash

./script_exp_gap.sh
perl tttplots.pl -f c_sol_file
perl tttplots.pl -f rmxdmc_sol_file
gnuplot ttt-plot-compare.gpl  