#!/usr/bin/gnuplot

trialfile = "output_data/gsa_test.txt"

f1(x) = 1.0 / 6.0 * x ** 6 - 52.0 / 25.0 * x ** 5 + 39.0 / 80.0 * x ** 4 + \
        71.0 / 10.0 * x ** 3 - 79.0 / 20.0 * x ** 2 - x + 1.0 / 10.0

f2(x) = sin(x) + sin(10.0 / 3.0 * x)

f3(x) = -sum[i = 1:5] i * sin((i + 1.0) * x  + i)

f4(x) = (-16.0 * x ** 2 + 24.0 * x - 5.0) * exp(-x)

f5(x) = (3.0 * x - 1.4) * sin(18.0 * x)

f6(x) = -(x + sin(x)) * exp(-x ** 2)

f7(x) = sin(x) + sin(10.0 / 3.0 * x) + log(x) - 0.84 * x + 3.0

f8(x) = -sum[i = 1:5] i * cos((i + 1.0) * x  + i)

f9(x) = sin(x) + sin(2.0 / 3.0 * x)

f10(x) = -x * sin(x)

f11(x) = 2.0 * cos(x) + cos(2.0 * x)

f12(x) = sin(x) ** 3 + cos(x) ** 3

f13(x) = -x ** (2.0 / 3.0) + (x ** 2 - 1) ** (1.0 / 3.0)

f14(x) = -exp(-x) * sin(2 * M_PI * x)

f15(x) = (x ** 2 - 5.0 * x + 6.0) / (x * x + 1)

f16(x) = 2.0 * (x - 3) ** 2 + exp(- x * x / 2)

f17(x) = x ** 6 - 15.0 * x ** 4 + 27.0 * x ** 2 + 250.0

f18(x) = (x <= 3.0) ? (x - 2.0) ** 2 : 2.0 * log(x - 2.0) + 1.0

f19(x) = -x + sin(3.0 * x) - 1.0

f20(x) = (-x + sin(x)) * exp(-x ** 2)

title_name(n) = sprintf("Graph of the test function №%d", n)
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
