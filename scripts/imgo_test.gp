#!/usr/bin/gnuplot

taskName = "imgo_test"
trialfile = "output_data/".taskName.".txt"
load "output_data/".taskName."_opt.txt"

f1Sample(i, x) = i == 0 ? sin(x) : 1/0

f2Sample(i, x) = i == 0 ? -2.0 * x + 3.0 : \
                 i == 1 ? sin(x) : 1/0

f3Sample(i, x) = i == 0 ? 5.0 * x ** 2 + 3.0 * x - 1.0 : \
                 i == 1 ? x ** 2 - 0.05 : \
                 i == 2 ? -x + 0.1 : 1/0

f1Test(i, x) = i == 0 ? -13.0 / 6.0 * x + sin(13.0 / 4.0 * (2.0 * x + 5.0)) - 53.0 / 12.0 : \
               i == 1 ? exp(-sin(3.0 * x)) - 1.0 / 10.0 * (x - 1.0 / 2.0) ** 2 - 1.0 : 1/0

f2Test(i, x) = i == 0 ? (11.0 * x * x - 10.0 * x + 21.0) / (2.0 * (x * x + 1)) : \
               i == 1 ? 1.0 / 20.0 - exp(-2.0 / 5.0 * (x + 5.0)) * sin(4.0 / 5.0 * pi * (x + 5.0)) : 1/0

f3Test(i, x) = i == 0 ? -sum[j = 1 : 5] cos(j * x) : \
               i == 1 ? 3.0 / 2.0 * (cos(7.0 / 20.0 * (x + 10.0)) - sin(7.0 / 4.0 * (x + 10.0)) + 1.0 / 2.0) : 1/0

f4Test_a(x) = pi / 2.0 * x + 1.0 / 10.0
f4Test(i, x) = i == 0 ? 4.0 * sin(f4Test_a(x) / 2.0) * (sin(f4Test_a(x)) ** 3 + cos(f4Test_a(x)) ** 3) ** 2 : \
               i == 1 ? 6.0 / 25.0 - sum[j = 1 : 5] cos(5.0 / 4.0 * (j + 1.0) * x + j) : \
               i == 2 ? 9.0 / 50.0 - 9.0 / 2.0 * exp(-(x - 1.0 / 10.0)) * sin(2.0 * pi * (x - 1.0 / 10.0)) : 1/0

f5Test(i, x) = i == 0 ? sin(0.423531 * x + 3.13531) + sin(10.0 / 3.0 * (0.423531 * x + 3.13531)) + \
                        log(0.423531 * x + 3.13531) + 0.36634 - 0.355766 * x : \
               i == 1 ? 17.0 / 25.0 - 2.0 / 29763.233 * (-1.0 / 6.0 * x ** 6 + 52.0 / 25.0 * x ** 5 - 39.0 / 80.0 * x ** 4 - \
                        71.0 / 10.0 * x ** 3 + 79.0 / 20.0 * x ** 2 + x - 1.0 / 10.0) : \
               i == 2 ? -14.0 / 125.0 * (3.0 * x - 8.0) * sin(252.0 / 125.0 * (x + 3.0 / 2.0)) - 1.0 / 2.0 : 1/0

f6Test(i, x) = i == 0 ? -7.0 / 40.0 * (3.0 * x + 4.0) * sin(63.0 / 20.0 * (x + 4.0)) : \
               i == 1 ? 40.0 * (cos(4.0 * x) * (x - sin(x)) * exp(-(x * x) / 2.0)) : \
               i == 2 ? 2.0 / 25.0 * (x + 4.0) - sin(12.0 / 5.0 * (x + 4.0)) : 1/0

f7Test(i, x) = i == 0 ? exp(-cos(4.0 * x - 3.0)) + 1.0 / 250.0 * (4.0 * x - 3.0) ** 2 - 1.0 : \
               i == 1 ? sin(x) ** 3 * exp(-sin(3.0 * x)) + 1.0 / 2.0 : \
               i == 2 ? cos(7.0 / 5.0 * (x + 3.0)) - sin(7.0 * (x + 3.0)) + 3.0 / 10.0 : 1/0

f8Test(i, x) = i == 0 ? cos(7.0 / 4.0 * x + 241.0 / 40.0) - sin(35.0 / 4.0 * x + 241.0 / 8.0) - 5.0 : \
               i == 1 ? exp(-sin(4.0 * x)) - 1.0 / 10.0 * (x - 1.0 / 2.0) ** 2 - 1.0 : \
               i == 2 ? 3.0 / 10.0 - sum[j = 1 : 5] cos(5.0 * (j + 1.0) * (x + 1.0 / 2.0)) : \
               i == 3 ? (-21.0 / 20.0 * x - 13.0 / 8.0) * sin(63.0 / 10.0 * x + 63.0 / 4.0) + 1.0 / 5.0 : 1/0

f9Test(i, x) = i == 0 ? sum[j = 1:5] 1.0 / 5.0 * sin((j + 1.0) * x - 1.0) + 2.0 : \
               i == 1 ? 1.0 / 40.0 * (x - 4.0) * (x - 32.0 / 5.0) * (x - 9.0) * (x - 11.0) * \
                        exp(-1.0 / 10.0 * (x - 13.0 / 2.0) ** 2) : \
               i == 2 ? (sin(x + 1.0) ** 3 + cos(x + 1.0) ** 3) * exp(-(x + 1.0) / 10.0) : \
               i == 3 ? exp(-cos(3.0 / 5.0 * (x - 5.0 / 2.0))) + 1.0 / 10.0 * (3.0 / 25.0 * x - 4.0 / 5.0) ** 2 - 1.0 : 1/0

f10Test_a(x) = 2.0 / pi * x - 1.0 / 2.0
f10Test_b(x) = 4.0 / pi * (x - 3.0 / 10.0) - 4.0
f10Test(i, x) = i == 0 ? -1.0 / 500.0 * f10Test_b(x) ** 6 + 3.0 / 100.0 * f10Test_b(x) ** 4 - 27.0 / 500.0 * f10Test_b(x) ** 2 + 3.0 / 2.0 : \
            i == 1 ? 2.0 * exp(-2.0 / pi * x) * sin(4.0 * x) : \
            i == 2 ? -f10Test_a(x) ** 2 * (-f10Test_a(x) ** 2 + 5.0 * f10Test_a(x) - 6.0) / (f10Test_a(x) ** 2 + 1.0) - 1.0 / 2.0 : \
            i == 3 ? sin(x) ** 3 + cos(2.0 * x) ** 3 - 3.0 / 10.0 : 1/0

f(i, j, x) = i == 1 ? f1Sample(j, x) : \
             i == 2 ? f2Sample(j, x) : \
             i == 3 ? f3Sample(j, x) : \
             i == 4 ? f1Test(j, x) : \
             i == 5 ? f2Test(j, x) : \
             i == 6 ? f3Test(j, x) : \
             i == 7 ? f4Test(j, x) : \
             i == 8 ? f5Test(j, x) : \
             i == 9 ? f6Test(j, x) : \
             i == 10 ? f7Test(j, x) : \
             i == 11 ? f8Test(j, x) : \
             i == 12 ? f9Test(j, x) : \
             i == 13 ? f10Test(j, x) : 1/0

set samples 600

set grid

set xlabel "x" font "Helvetica, 16"
set ylabel "f(x)" font "Helvetica, 16" offset -1

set tics font "Helvetica, 16"

set key box outside right top
set key font "Helvetica, 16" spacing 1.3

title(block, n) = sprintf("Graph of the %s function №%d, method gsa", block, n)
titlePng(block, n) = ARG1 == 1 ? title(block, n) : sprintf("")
functionName(i) = (i == 0 ? sprintf("f(x)") : sprintf("g%d(x)", i))

if (ARG1 == 0) {
    set xrange [A[ARG2 + 1] : B[ARG2 + 1]]

    set title title(functionBlockName[ARG2 + 1], functionNumber[ARG2 + 1]) font "Helvetica, 16"

    ind = 3 * ARG2
    plot for [i = 0 : numberConstraints[ARG2 + 1]] f(ARG2 + 1, i, x) title functionName(i), \
         trialfile index ind + 2 ls 4 lc rgb "green" lw 2 title "trials", \
         trialfile index ind     ls 7 lc rgb "red"   ps 2 title "X_{min}", \
         trialfile index ind + 1 ls 7 lc rgb "blue"  ps 1 title "X"

    bind all "alt-End" "exit gnuplot"
    pause mouse close
} else {
    set terminal pngcairo size 1640, 950 font "Helvetica, 16"
    system "mkdir -p output_graph/".taskName

    set lmargin 12
    set rmargin 16
    set tmargin 3
    set bmargin 3

    do for [i = 1 : 13] {
        set output "output_graph/".taskName."/".functionBlockName[i]."_".functionNumber[i].".png"

        set xrange [A[i] : B[i]]

        set title titlePng(functionBlockName[i], functionNumber[i]) font "Helvetica, 16"

        ind = 3 * (i - 1)
        plot for [j = 0 : numberConstraints[i]] f(i, j, x) title functionName(j), \
             trialfile index ind + 2 ls 4 lc rgb "green" lw 2 title "trials", \
             trialfile index ind     ls 7 lc rgb "red"   lw 6 title "X_{min}", \
             trialfile index ind + 1 ls 7 lc rgb "blue"  lw 1 title "X"
    }
}
