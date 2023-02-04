#!/usr/bin/gnuplot

trialfile = "output_data/imgo_test.txt"

load "output_data/imgo_test_opt.txt"

f1(x) = -13.0 / 6.0 * x + sin(13.0 / 4.0 * (2.0 * x + 5.0)) - 53.0 / 12.0
g1_1(x) = exp(-sin(3.0 * x)) - 1.0 / 10.0 * (x - 1.0 / 2.0) ** 2 - 1.0

f2(x) = (11.0 * x * x - 10.0 * x + 21.0) / (2.0 * (x * x + 1))
g2_1(x) = 1.0 / 20.0 - exp(-2.0 / 5.0 * (x + 5.0)) * sin(4.0 / 5.0 * pi * (x + 5.0))

f3(x) = -sum[i = 1:5] cos(i * x)
g3_1(x) = 3.0 / 2.0 * (cos(7.0 / 20.0 * (x + 10.0)) - sin(7.0 / 4.0 * (x + 10.0)) + 1.0 / 2.0)

f4_a(x) = pi / 2.0 * x + 1.0 / 10.0
f4(x) = 4.0 * sin(f4_a(x) / 2.0) * (sin(f4_a(x)) ** 3 + cos(f4_a(x)) ** 3) ** 2
g4_1(x) = 6.0 / 25.0 - sum[i = 1:5] cos(5.0 / 4.0 * (i + 1.0) * x + i)
g4_2(x) = 9.0 / 50.0 - 9.0 / 2.0 * exp(-(x - 1.0 / 10.0)) * sin(2.0 * pi * (x - 1.0 / 10.0))

f5(x) = sin(0.423531 * x + 3.13531) + sin(10.0 / 3.0 * (0.423531 * x + 3.13531)) + log(0.423531 * x + 3.13531) + 0.36634 - 0.355766 * x
g5_1(x) = 17.0 / 25.0 - 2.0 / 29763.233 * (-1.0 / 6.0 * x ** 6 + 52.0 / 25.0 * x ** 5 - 39.0 / 80.0 * x ** 4 - \
          71.0 / 10.0 * x ** 3 + 79.0 / 20.0 * x ** 2 + x - 1.0 / 10.0)
g5_2(x) = -14.0 / 125.0 * (3.0 * x - 8.0) * sin(252.0 / 125.0 * (x + 3.0 / 2.0)) - 1.0 / 2.0

f6(x) = -7.0 / 40.0 * (3.0 * x + 4.0) * sin(63.0 / 20.0 * (x + 4.0))
g6_1(x) = 40.0 * (cos(4.0 * x) * (x - sin(x)) * exp( -(x * x) / 2.0))
g6_2(x) = 2.0 / 25.0 * (x + 4.0) - sin(12.0 / 5.0 * (x + 4.0))

f7(x) = exp(-cos(4.0 * x - 3.0)) + 1.0 / 250.0 * (4.0 * x - 3.0) ** 2 - 1.0
g7_1(x) = sin(x) ** 3 * exp(-sin(3.0 * x)) + 1.0 / 2.0
g7_2(x) = cos(7.0 / 5.0 * (x + 3.0)) - sin(7.0 * (x + 3.0)) + 3.0 / 10.0

f8(x) = cos(7.0 / 4.0 * x + 241.0 / 40.0) - sin(35.0 / 4.0 * x + 241.0 / 8.0) - 5.0
g8_1(x) = exp(-sin(4.0 * x)) - 1.0 / 10.0 * pow(x - 1.0 / 2.0, 2) - 1.0
g8_2(x) = 3.0 / 10.0 - sum[i = 1:5] cos(5.0 * (i + 1.0) * (x + 1.0 / 2.0))
g8_3(x) = (-21.0 / 20.0 * x - 13.0 / 8.0) * sin(63.0 / 10.0 * x + 63.0 / 4.0) + 1.0 / 5.0

f9(x) = sum[i = 1:5] 1.0 / 5.0 * sin((i + 1.0) * x - 1.0) + 2.0
g9_1(x) = 1.0 / 40.0 * (x - 4.0) * (x - 32.0 / 5.0) * (x - 9.0) * (x - 11.0) * exp(-1.0 / 10.0 * (x - 13.0 / 2.0) ** 2)
g9_2(x) = (sin(x + 1.0) ** 3 + cos(x + 1.0) ** 3) * exp(-(x + 1.0) / 10.0)
g9_3(x) = exp(-cos(3.0 / 5.0 * (x - 5.0 / 2.0))) + 1.0 / 10.0 * (3.0 / 25.0 * x - 4.0 / 5.0) ** 2 - 1.0

f10_a(x) = 2.0 / pi * x - 1.0 / 2.0
f10_b(x) = 4.0 / pi * (x - 3.0 / 10.0) - 4.0
f10(x) = -1.0 / 500.0 * f10_b(x) ** 6 + 3.0 / 100.0 * f10_b(x) ** 4 - 27.0 / 500.0 * f10_b(x) ** 2 + 3.0 / 2.0
g10_1(x) = 2.0 * exp(-2.0 / pi * x) * sin(4.0 * x)
g10_2(x) = -f10_a(x) ** 2 * (-f10_a(x) ** 2 + 5.0 * f10_a(x) - 6.0) / (f10_a(x) ** 2 + 1.0) - 1.0 / 2.0
g10_3(x) = sin(x) ** 3 + cos(2.0 * x) ** 3 - 3.0 / 10.0

title_name(n) = sprintf("Graph of the test function №%d", n)

set grid

set samples 600

set xlabel "x"
set ylabel "y"

ind = 3 * ARG1
function = functions[ARG1 + 1]

set title title_name(ARG1 + 1) font "Helvetica Bold, 20"

plot @function, \
     trialfile index ind + 2 ls 5 lc rgb "green" title "trials", \
     trialfile index ind + 1 ls 5 lc rgb "blue" title "X", \
     trialfile index ind ls 5 lc rgb "red" title "X*"

bind all "alt-End" "exit gnuplot"
pause mouse close
