#! /usr/bin/gnuplot

nameTask = "mggsa_test_evolvent"

datafile = "output_data/".nameTask.".txt"

peanofile = "output_data/".nameTask."_peano.txt"

f(x, y) = x ** 2 + y ** 2 - cos(18.0 * x) - cos(18.0 * y)

graphType(n) = (n == 2) ? "with lines" : ""

set contour
set view map
set cntrparam bspline levels auto 14
set cntrlabel onecolor
set cntrlabel start 5 interval 150
set cntrlabel font ",10"
set contour base

set grid

set xlabel "X"
set ylabel "Y"

set isosamples 120

set key box
set key opaque

set xrange [-0.5 : 1.0]
set yrange [-0.5 : 1.0]

type = graphType(ARG1)

set title "Evolvent test, method mggsa" font "Helvetica Bold, 20"

if (ARG1 == 1) {
    set terminal pngcairo size 950, 950
    system "mkdir -p output_graph/".nameTask
    set output "output_graph/".nameTask."/".nameTask.".png"
}

splot f(x, y) title "f(x, y)" nosurface, \
      peanofile index 0 using 1:2:(0) ls 5 lc rgb "blue" ps 1  title "evolvent" @type nocontours, \
      datafile index 2  using 1:2:(0) ls 4 lc rgb "green" lw 3 title "trial points" nocontours, \
      datafile index 0  using 1:2:(0) ls 7 lc rgb "red"        title "X_{min}" nocontours, \
      datafile index 1  using 1:2:(0) ls 7 lc rgb "orange"     title "X" nocontours

if (ARG1 == 0) {
    bind all "alt-End" "exit gnuplot"
    pause mouse close
}
