#!/usr/bin/gnuplot

taskName = "gsa_test"
trialfile = "output_data/".taskName.".txt"
load "output_data/".taskName."_opt.txt"

f1Sample(x) = -4.0 * x + 1.0

f2Sample(x) = 5.0 * x ** 2 + 3.0 * x - 1.0

f3Sample(x) = x * sin(x)

f4Sample(x) = x != 0 ? x * sin(1 / x) : 0.0

f1Test(x) = 1.0 / 6.0 * x ** 6 - 52.0 / 25.0 * x ** 5 + 39.0 / 80.0 * x ** 4 + \
            71.0 / 10.0 * x ** 3 - 79.0 / 20.0 * x ** 2 - x + 1.0 / 10.0

f2Test(x) = sin(x) + sin(10.0 / 3.0 * x)

f3Test(x) = -sum[i = 1 : 5] i * sin((i + 1.0) * x  + i)

f4Test(x) = (-16.0 * x ** 2 + 24.0 * x - 5.0) * exp(-x)

f5Test(x) = (3.0 * x - 1.4) * sin(18.0 * x)

f6Test(x) = -(x + sin(x)) * exp(-x ** 2)

f7Test(x) = sin(x) + sin(10.0 / 3.0 * x) + log(x) - 0.84 * x + 3.0

f8Test(x) = -sum[i = 1 : 5] i * cos((i + 1.0) * x  + i)

f9Test(x) = sin(x) + sin(2.0 / 3.0 * x)

f10Test(x) = -x * sin(x)

f11Test(x) = 2.0 * cos(x) + cos(2.0 * x)

f12Test(x) = sin(x) ** 3 + cos(x) ** 3

f13Test(x) = x * x - 1 < 0 ? -x ** (2.0 / 3.0) - (-x ** 2 + 1) ** (1.0 / 3.0) : \
                             -x ** (2.0 / 3.0) + (x ** 2 - 1) ** (1.0 / 3.0)

f14Test(x) = -exp(-x) * sin(2 * pi * x)

f15Test(x) = (x ** 2 - 5.0 * x + 6.0) / (x * x + 1)

f16Test(x) = 2.0 * (x - 3) ** 2 + exp(-x * x / 2)

f17Test(x) = x ** 6 - 15.0 * x ** 4 + 27.0 * x ** 2 + 250.0

f18Test(x) = x <= 3.0 ? (x - 2.0) ** 2 : 2.0 * log(x - 2.0) + 1.0

f19Test(x) = -x + sin(3.0 * x) - 1.0

f20Test(x) = (-x + sin(x)) * exp(-x ** 2)

f(i, x) = i == 1 ? f1Sample(x) : \
          i == 2 ? f2Sample(x) : \
          i == 3 ? f3Sample(x) : \
          i == 4 ? f4Sample(x) : \
          i == 5 ? f1Test(x) : \
          i == 6 ? f2Test(x) : \
          i == 7 ? f3Test(x) : \
          i == 8 ? f4Test(x) : \
          i == 9 ? f5Test(x) : \
          i == 10 ? f6Test(x) : \
          i == 11 ? f7Test(x) : \
          i == 12 ? f8Test(x) : \
          i == 13 ? f9Test(x) : \
          i == 14 ? f10Test(x) : \
          i == 15 ? f11Test(x) : \
          i == 16 ? f12Test(x) : \
          i == 17 ? f13Test(x) : \
          i == 18 ? f14Test(x) : \
          i == 19 ? f15Test(x) : \
          i == 20 ? f16Test(x) : \
          i == 21 ? f17Test(x) : \
          i == 22 ? f18Test(x) : \
          i == 23 ? f19Test(x) : \
          i == 24 ? f20Test(x) : 1/0

set sample 400

set grid

set xlabel "x" font "Helvetica, 16"
set ylabel "f(x)" font "Helvetica, 16" offset -1

set tics font "Helvetica, 16"

set key box outside right top
set key font "Helvetica, 16" spacing 1.3

title(block, n) = sprintf("Graph of the %s function №%d, method gsa", block, n)
titlePng(block, n) = ARG1 == 1 ? title(block, n) : sprintf("")

if (ARG1 == 0) {
    set xrange [A[ARG2 + 1] : B[ARG2 + 1]]

    set title title(functionBlockName[ARG2 + 1], functionNumber[ARG2 + 1]) font "Helvetica, 16"

    ind = 3 * ARG2
    plot f(ARG2 + 1, x) title "f(x)", \
         trialfile index ind + 2 ls 4 lc rgb "green" lw 2 title "trial points", \
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

    do for [i = 1 : 24] {
        set output "output_graph/".taskName."/".functionBlockName[i]."_".functionNumber[i].".png"

        set xrange [A[i] : B[i]]

        set title titlePng(functionBlockName[i], functionNumber[i]) font "Helvetica, 16"

        ind = 3 * (i - 1)
        plot f(i, x) title "f(x)", \
             trialfile index ind + 2 ls 4 lc rgb "green" lw 2 title "trials", \
             trialfile index ind     ls 7 lc rgb "red"   lw 6 title "X_{min}", \
             trialfile index ind + 1 ls 7 lc rgb "blue"  lw 1 title "X"
    }
}
