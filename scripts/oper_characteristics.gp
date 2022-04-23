#! /usr/bin/gnuplot

reset

load "operational_characteristics_opt.txt"
datafile="operational_characteristics.txt"

p_hill="hill-P_s(k)-r="
title_hill_1=p_hill.r1_hill
title_hill_2=p_hill.r2_hill
title_hill_3=p_hill.r3_hill

p_shekel="shekel-P_s(k)-r="
title_shekel_1=p_shekel.r1_shekel
title_shekel_2=p_shekel.r2_shekel
title_shekel_3=p_shekel.r3_shekel

set xlabel "K"
set ylabel "P_s(k)"
set grid

set title "Operational characteristics for imgo with peano" font "Helvetica Bold, 20"
plot datafile index 0 using 1:2 with lines ls 5 lc rgb "red" title title_hill_1, \
     datafile index 1 using 1:2 with lines ls 5 lc rgb "green" title title_hill_2, \
     datafile index 2 using 1:2 with lines ls 5 lc rgb "blue" title title_hill_3, \
     datafile index 3 using 1:2 with lines ls 5 lc rgb "orange" title title_shekel_1, \
     datafile index 4 using 1:2 with lines ls 5 lc rgb "black" title title_shekel_2, \
     datafile index 5 using 1:2 with lines ls 5 lc rgb "violet" title title_shekel_3, \

bind all "alt-End" "exit gnuplot"
pause mouse close
