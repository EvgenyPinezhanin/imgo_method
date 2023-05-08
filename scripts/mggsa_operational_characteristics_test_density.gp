#! /usr/bin/gnuplot

datafile = "output_data/mggsa_operational_characteristics_test_density.txt"

load "output_data/mggsa_operational_characteristics_test_density_opt.txt"

set linetype 1 lc rgb "red"
set linetype 2 lc rgb "green"
set linetype 3 lc rgb "blue"
set linetype 4 lc rgb "orange"
set linetype 5 lc rgb "dark-violet"
set linetype 6 lc rgb "salmon"
set linetype 7 lc rgb "cyan"

set linetype cycle 7

set grid
set xlabel "K" font ", 15"
set ylabel "P_s(k)" font ", 15"
set tics font ", 11"
set key font ", 15"

set title "Operational characteristics for mggsa on a family of tasks ".familyNames[ARG1 + 1]." (density test)" font "Helvetica Bold, 20"
plot for [i = 1:sizeDen] datafile index sizeDen * ARG1 + i - 1 using 1:2 with lines lt i title "den = ".den[i]

bind all "alt-End" "exit gnuplot"
pause mouse close
