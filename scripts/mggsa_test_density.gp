#! /usr/bin/gnuplot

nameTask = "mggsa_test_density"

datafile = "output_data/".nameTask.".txt"

load "output_data/".nameTask."_opt.txt"

set linetype 1 lc rgb "dark-yellow" pt 4 lw 2
set linetype 2 lc rgb "cyan"        pt 4 lw 2
set linetype 3 lc rgb "brown"       pt 4 lw 2
set linetype 4 lc rgb "orange"      pt 4 lw 2
set linetype 5 lc rgb "green"       pt 4 lw 2

set linetype 6  lc rgb "dark-yellow" pt 7 lw 3 ps 2
set linetype 7  lc rgb "cyan"        pt 7 lw 3 ps 2
set linetype 8  lc rgb "brown"       pt 7 lw 3 ps 2
set linetype 9  lc rgb "orange"      pt 7 lw 3 ps 2
set linetype 10 lc rgb "green"       pt 7 lw 3 ps 2

set linetype 11 lc rgb "dark-violet"      lw 1
set linetype 12 lc rgb "red"         pt 7 lw 5
set linetype 13 lc rgb "blue"        pt 7 lw 5 ps 1

set linetype cycle 13

f(x, y) = x ** 2 + y ** 2 - cos(18.0 * x) - cos(18.0 * y)

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

set title "Density test, method mggsa" font "Helvetica Bold, 20"

if (ARG1 == 0) {
    splot f(x, y) title "f(x, y)" lt 11 nosurface, \
          for [i = 1 : @ARG2] datafile index 2 * (ARG2 - i + 1) lt i title "trial points(m = ".den[ARG2 - i + 1].")" nocontours, \
          datafile index 0 lt 12 ps 3 title "X_{min}" nocontours, \
          for [i = 1 : @ARG2] datafile index 2 * (ARG2 - i + 1) - 1 lt i + 5 title "X(m = ".den[ARG2 - i + 1].")" nocontours

    bind all "alt-End" "exit gnuplot"
    pause mouse close
} else {
    set term gif animate optimize delay 200 size 950, 950 font "Helvetica Bold, 16"
    system "mkdir -p output_graph/".nameTask

    set output "output_graph/".nameTask."/".nameTask.".gif"

    do for [i = 1 : @ARG2] {
        splot f(x, y) title "f(x, y)" lt 11 nosurface, \
              for [j = 1 : i] datafile index 2 * (i + 1 - j) lt i + 1 - j title "trial points(m = ".den[i + 1 - j].")" nocontours, \
              datafile index 0 lt 12 ps 2 title "X_{min}" nocontours, \
              datafile index 2 * i - 1 lt 13 title "X" nocontours 
    }
}
