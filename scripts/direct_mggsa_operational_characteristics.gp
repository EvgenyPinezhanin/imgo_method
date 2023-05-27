#! /usr/bin/gnuplot

datafile_direct = "output_data/direct_mggsa_operational_characteristics_direct.txt"
datafile_mggsa = "output_data/direct_mggsa_operational_characteristics_mggsa.txt"

load "output_data/direct_operational_characteristics_opt.txt"

set linetype 1 lc rgb "red"
set linetype 2 lc rgb "green"

set linetype 3 lc rgb "blue"
set linetype 4 lc rgb "orange"
set linetype 5 lc rgb "cyan"
set linetype 6 lc rgb "dark-violet"

set linetype cycle 6

set grid

set xlabel "K" font ", 15"
set ylabel "P_s(K)" font ", 15"

set tics font ", 11"
set key font ", 15"

array AlghorithmDirect[2]
AlghorithmDirect[1] = "ORIGINAL"
AlghorithmDirect[2] = "GABLONSKY"

ind = numberKey * @ARG1

set title "Operational characteristics for DIRECT and mggsa on a family of tasks ".familyNames[ARG1 + 1] font "Helvetica Bold, 20"
plot for [i = 1:2] datafile_direct index 2 * @ARG1 + i - 1 using 1:2 with lines lt i title "DIRECT ".AlghorithmDirect[i], \
     for [i = 1:numberKey] datafile_mggsa index ind + i - 1 using 1:2 with lines lt i + 2 title "MGGSA "."r = ".r[ind + i].", key = ".key[i]

bind all "alt-End" "exit gnuplot"
pause mouse close
