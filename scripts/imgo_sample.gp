#!/usr/bin/gnuplot

trialfile = "output_data/imgo_sample.txt"

load "output_data/imgo_sample_opt.txt"

f1(x) = sin(x)

f2(x) = -2.0 * x + 3.0
g2_1(x) = sin(x)

f3(x) = 5.0 * x ** 2 + 3.0 * x - 1.0
g3_1(x) = x ** 2 - 0.05
g3_2(x) = -x + 0.1

titleName(n) = sprintf("Graph of the sample function №%d", n)

set grid
set xlabel "x"
set ylabel "y"

ind = 3 * ARG1
function = functions[ARG1 + 1]

set title titleName(ARG1 + 1) font "Helvetica Bold, 20"

plot @function, \
     trialfile index ind + 2 ls 5 lc rgb "green" title "trials", \
     trialfile index ind + 1 ls 5 lc rgb "blue" title "X", \
     trialfile index ind ls 5 lc rgb "red" title "X*"

bind all "alt-End" "exit gnuplot"
pause mouse close
