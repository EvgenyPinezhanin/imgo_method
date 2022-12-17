#! /usr/bin/gnuplot

datafile="output_data/mggsa_operational_characteristics.txt"

array Name[4]
array R[3*4]
load "output_data/mggsa_operational_characteristics_opt.txt"

ind=3*ARG1

set xlabel "K" font ", 15"
set ylabel "P_s(k)" font ", 15"
set grid

set tics font ", 11"
set key font ", 15"

if (ARG1 < 4) {
     set title "Operational characteristics for imgo with peano on a family of tasks ".Name[ARG1+1] font "Helvetica Bold, 20"
     plot datafile index ind using 1:2 with lines ls 5 lc rgb "red" title "r=".R[ind+1], \
          datafile index ind+1 using 1:2 with lines ls 5 lc rgb "green" title "r=".R[ind+2], \
          datafile index ind+2 using 1:2 with lines ls 5 lc rgb "blue" title "r=".R[ind+3]
}
if (ARG1 == 4) {
     set title "Comparison of operational characteristics Grishagin and GKLS" font "Helvetica Bold, 20"
     plot datafile index 0 using 1:2 with lines ls 5 lc rgb "red" title Name[1]."-r=".R[1], \
          datafile index 1 using 1:2 with lines ls 5 lc rgb "green" title Name[1]."-r=".R[2], \
          datafile index 2 using 1:2 with lines ls 5 lc rgb "blue" title Name[1]."-r=".R[3], \
          datafile index 3 using 1:2 with lines ls 5 lc rgb "orange" title Name[2]."-r=".R[4], \
          datafile index 4 using 1:2 with lines ls 5 lc rgb "black" title Name[2]."-r=".R[5], \
          datafile index 5 using 1:2 with lines ls 5 lc rgb "violet" title Name[2]."-r=".R[6]
}
if (ARG1 == 5) {
     set title "Comparison of operational characteristics Grishagin and GKLS (constrained)" font "Helvetica Bold, 20"
     plot datafile index 6 using 1:2 with lines ls 5 lc rgb "red" title Name[3]."-r=".R[7], \
          datafile index 7 using 1:2 with lines ls 5 lc rgb "green" title Name[3]."-r=".R[8], \
          datafile index 8 using 1:2 with lines ls 5 lc rgb "blue" title Name[3]."-r=".R[9], \
          datafile index 9 using 1:2 with lines ls 5 lc rgb "orange" title Name[4]."-r=".R[10], \
          datafile index 10 using 1:2 with lines ls 5 lc rgb "black" title Name[4]."-r=".R[11], \
          datafile index 11 using 1:2 with lines ls 5 lc rgb "violet" title Name[4]."-r=".R[12]
}

bind all "alt-End" "exit gnuplot"
pause mouse close
