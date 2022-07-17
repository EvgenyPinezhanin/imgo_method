#! /usr/bin/gnuplot

load "output_data/imgo_operational_characteristics_opt.txt"
datafile="output_data/imgo_operational_characteristics.txt"

p_hill="hill-P_s(k)-r="
title_hill_1=p_hill.r1_hill
title_hill_2=p_hill.r2_hill
title_hill_3=p_hill.r3_hill

p_shekel="shekel-P_s(k)-r="
title_shekel_1=p_shekel.r1_shekel
title_shekel_2=p_shekel.r2_shekel
title_shekel_3=p_shekel.r3_shekel

set xlabel "K" font ", 13"
set ylabel "P_s(k)" font ", 13"
set grid

set tics font ", 13"
set key font ", 15"

if (ARG1 == 0) {
     set title "Operational characteristics on a family of tasks Hill" font "Helvetica Bold, 25"
     plot datafile index 0 using 1:2 with lines ls 5 lc rgb "red" title title_hill_1, \
          datafile index 1 using 1:2 with lines ls 5 lc rgb "green" title title_hill_2, \
          datafile index 2 using 1:2 with lines ls 5 lc rgb "blue" title title_hill_3
}
if (ARG1 == 1) {
     set title "Operational characteristics on a family of tasks Shekel" font "Helvetica Bold, 20"
     plot datafile index 3 using 1:2 with lines ls 5 lc rgb "red" title title_shekel_1, \
          datafile index 4 using 1:2 with lines ls 5 lc rgb "green" title title_shekel_2, \
          datafile index 5 using 1:2 with lines ls 5 lc rgb "blue" title title_shekel_3
}
if (ARG1 == 2) {
     set title "Comparison of operational characteristics" font "Helvetica Bold, 20"
     plot datafile index 0 using 1:2 with lines ls 5 lc rgb "red" title title_hill_1, \
          datafile index 1 using 1:2 with lines ls 5 lc rgb "green" title title_hill_2, \
          datafile index 2 using 1:2 with lines ls 5 lc rgb "blue" title title_hill_3, \
          datafile index 3 using 1:2 with lines ls 5 lc rgb "orange" title title_shekel_1, \
          datafile index 4 using 1:2 with lines ls 5 lc rgb "black" title title_shekel_2, \
          datafile index 5 using 1:2 with lines ls 5 lc rgb "violet" title title_shekel_3
}

bind all "alt-End" "exit gnuplot"
pause mouse close
