#! /usr/bin/gnuplot

reset

datafile = "output_data/trial_points.txt"
peanofile = "output_data/map_test.txt"

set hidden3d front
set grid
set xlabel "X"
set ylabel "Y"
set xrange [-0.5:1.0]
set yrange [-0.5:1.0]
set key box
set key opaque
set title "Map test" font "Helvetica Bold, 20"

if (ARG1 == 1) {
plot peanofile index 0 ls 5 lc rgb "blue" title "evolvent", \
     datafile index 2 ls 5 lc rgb "green" title "trial points", \
     datafile index 1 ls 5 lc rgb "red" title "X*", \
     datafile index 0 ls 5 lc rgb "orange" title "X"
}
if (ARG1 == 2) {
plot peanofile index 0 ls 5 lc rgb "blue" title "evolvent" with lines, \
     datafile index 2 ls 5 lc rgb "green" title "trial points", \
     datafile index 1 ls 5 lc rgb "red" title "X*", \
     datafile index 0 ls 5 lc rgb "orange" title "X"
}
if (ARG1 == 3) {
plot peanofile index 0 ls 5 lc rgb "blue" title "evolvent", \
     datafile index 2 ls 5 lc rgb "green" title "trial points", \
     datafile index 1 ls 5 lc rgb "red" title "X*", \
     datafile index 0 ls 5 lc rgb "orange" title "X"
}

bind all "alt-End" "exit gnuplot"
pause mouse close
