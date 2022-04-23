#! /usr/bin/gnuplot

reset

load "peano_operational_characteristics_opt.txt"
datafile="peano_operational_characteristics.txt"

p_grish="grish-P_s(k)-r="
title_grish_1=p_grish.r1_grish
title_grish_2=p_grish.r2_grish
title_grish_3=p_grish.r3_grish

p_gkls="gkls-P_s(k)-r="
title_gkls_1=p_gkls.r1_gkls
title_gkls_2=p_gkls.r2_gkls
title_gkls_3=p_gkls.r3_gkls

set xlabel "K"
set ylabel "P_s(k)"
set grid

set title "Operational characteristics for imgo" font "Helvetica Bold, 20"
plot datafile index 0 using 1:2 with lines ls 5 lc rgb "red" title title_grish_1, \
     datafile index 1 using 1:2 with lines ls 5 lc rgb "green" title title_grish_2, \
     datafile index 2 using 1:2 with lines ls 5 lc rgb "blue" title title_grish_3, \
     datafile index 3 using 1:2 with lines ls 5 lc rgb "orange" title title_gkls_1, \
     datafile index 4 using 1:2 with lines ls 5 lc rgb "black" title title_gkls_2, \
     datafile index 5 using 1:2 with lines ls 5 lc rgb "violet" title title_gkls_3, \

bind all "alt-End" "exit gnuplot"
pause mouse close
