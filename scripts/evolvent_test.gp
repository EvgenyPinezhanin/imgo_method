#! /usr/bin/gnuplot

datafile = "output_data/evolvent_test_points.txt"
peanofile = "output_data/evolvent_test.txt"

graph_type(n) = (n == 2) ? "with lines" : ""

set hidden3d front
set grid

set xlabel "X"
set ylabel "Y"

set key box
set key opaque

set xrange [-0.5:1.0]
set yrange [-0.5:1.0]

type = graph_type(ARG1)

set title "Evolvent test" font "Helvetica Bold, 20"

plot peanofile index 0 ls 5 lc rgb "blue" title "evolvent" @type, \
     datafile index 2 ls 5 lc rgb "green" title "trial points", \
     datafile index 1 ls 5 lc rgb "orange" title "X", \
     datafile index 0 ls 5 lc rgb "red" title "X*"

bind all "alt-End" "exit gnuplot"
pause mouse close
