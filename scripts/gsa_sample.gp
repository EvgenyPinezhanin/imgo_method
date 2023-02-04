#!/usr/bin/gnuplot

trialfile = "output_data/gsa_sample.txt"

f1(x) = -4.0 * x + 1.0

f2(x) = 5.0 * x ** 2 + 3.0 * x - 1.0

f3(x) = x * sin(x)

f4(x) = (x != 0) ? x * sin(1 / x) : 0.0

title_name(n) = sprintf("Graph of the sample function №%d", n)
function_name(n) = sprintf("f%d(x)", n)

set grid

set xlabel "x"
set ylabel "y"

ind = 3 * ARG1
function_name = function_name(ARG1 + 1)

set title title_name(ARG1 + 1) font "Helvetica Bold, 20"

plot @function_name title "f(x)", \
     trialfile index ind + 2 ls 5 lc rgb "green" title "trials", \
     trialfile index ind + 1 ls 5 lc rgb "blue" title "X", \
     trialfile index ind ls 5 lc rgb "red" title "X*"

bind all "alt-End" "exit gnuplot"
pause mouse close
