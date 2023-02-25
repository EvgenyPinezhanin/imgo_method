#!/usr/bin/gnuplot

trialfile = "output_data/direct_sample.txt"

load "output_data/direct_sample_opt.txt"

f1(x, y) = 1.0 - x - y
g1(x, y) = 1 / 0

f2(x, y) = (x - 1.0) ** 2 / 5.0 + (y - 1.0) ** 2 / 5.0
g2(x, y) = 1 / 0

f3(x, y) = x ** 2 / 5.0 + y ** 2 / 5.0
g3_1(x, y) = 1.0 - x - y
g3(x, y) = (g3_1(x, y) <= 0.0) ? g3_1(x, y) : 1 / 0

f4(x, y) = x ** 2 / 5.0 + y ** 2 / 5.0
g4_1(x, y) = (x - 2.0) ** 2 + (y - 2.0) ** 2 - 2.0
g4(x, y) = (g4_1(x, y) <= 0.0) ? g4_1(x, y) : 1 / 0

title_name(n) = sprintf("Graph of the sample function №%d with DIRECT", n)
function_name(s, n) = sprintf("%s%d(x, y)", s, n)

set grid
set contour
set view map
set cntrparam bspline levels auto 10
set cntrlabel onecolor
set cntrlabel start 5 interval 150
set cntrlabel font ",10"
set contour base
set isosamples 100

set xlabel "X"
set ylabel "Y"

ind = 3 * ARG1
function_f = function_name("f", ARG1 + 1)
function_g = function_name("g", ARG1 + 1)

set title title_name(ARG1 + 1) font "Helvetica Bold, 20"

set xrange [AX[ARG1 + 1]:BX[ARG1 + 1]]
set yrange [AY[ARG1 + 1]:BY[ARG1 + 1]]

splot @function_f title "φ(x, y)" nosurface, \
      @function_g lc rgb "orange" title "g(x, y)" nocontours, \
      trialfile index ind + 2 ls 5 lc rgb "green" title "trial points" nocontours, \
      trialfile index ind + 1 ls 5 lc rgb "blue" title "X" nocontours, \
      trialfile index ind ls 5 lc rgb "red" title "X*" nocontours, \
      @function_f with labels notitle nosurface

bind all "alt-End" "exit gnuplot"
pause mouse close
