#!/usr/bin/gnuplot

taskName = "direct_test"
trialfile = "output_data/".taskName.".txt"
load "output_data/".taskName."_opt.txt"

f1Sample(i, x, y) = (i == 0 ? 1.0 - x - y : 1/0)

f2Sample(i, x, y) = (i == 0 ? (x - 1.0) ** 2 / 5.0 + (y - 1.0) ** 2 / 5.0 : 1/0)

g3Sample_1(x, y) = 1.0 - x - y
g3Sample(x, y) = g3Sample_1(x, y) <= 0.0 ? g3Sample_1(x, y) : 1/0
f3Sample(i, x, y) = i == 0 ? x ** 2 / 5.0 + y ** 2 / 5.0 : \
                    i == 1 ? g3Sample(x, y) : 1/0

g4Sample_1(x, y) = (x - 2.0) ** 2 + (y - 2.0) ** 2 - 2.0
g4Sample(x, y) = g4Sample_1(x, y) <= 0.0 ? g4Sample_1(x, y) : 1/0
f4Sample(i, x, y) = i == 0 ? x ** 2 / 5.0 + y ** 2 / 5.0 : \
                    i == 1 ? g4Sample(x, y) : 1/0

load "scripts/functions/multidimensional/test/1.txt"
load "scripts/functions/multidimensional/test/2.txt"

f(i, j, x, y) = (i == 1 ? f1Sample(j, x, y) : \
                 i == 2 ? f2Sample(j, x, y) : \
                 i == 3 ? f3Sample(j, x, y) : \
                 i == 4 ? f4Sample(j, x, y) : \
                 i == 5 ? f1MTest(j, x, y) : \
                 i == 6 ? f2MTest(j, x, y) : 1/0)

set contour base
set cntrlabel onecolor
set cntrparam levels auto 10
set cntrlabel start 5 interval 40
set cntrlabel font "Helvetica, 10"

set view map

set isosamples 120

set grid

set xlabel "X" font "Helvetica, 16"
set ylabel "Y" font "Helvetica, 16" offset -1

set tics font "Helvetica, 16"

set key box outside right top
set key font "Helvetica, 16" spacing 1.5

title(block, n) = sprintf("Graph of the %s function №%d, method DIRECT", block, n)
titlePng(block, n) = ARG1 == 1 ? title(block, n) : sprintf("")

if (ARG1 == 0) {
    set key width 3

    set xrange [A[ARG2 + 1] : B[ARG2 + 1]]
    set yrange [A[ARG2 + 2] : B[ARG2 + 2]]

    set title title(functionBlockName[ARG2 + 1], functionNumber[ARG2 + 1]) font "Helvetica, 16"

    ind = 3 * ARG2
    splot f(ARG2 + 1, 0, x, y) title "φ(x, y)" nosurface, \
          f(ARG2 + 1, 1, x, y) lc rgb "orange" title "g(x, y)" nocontours, \
          trialfile index ind + 2 ls 4 lc rgb "green" lw 2 title "trial points" nocontours, \
          trialfile index ind     ls 7 lc rgb "red"   ps 2 title "X_{min}" nocontours, \
          trialfile index ind + 1 ls 7 lc rgb "blue"  ps 1 title "X" nocontours, \
          f(ARG2 + 1, 0, x, y) with labels notitle nosurface

    bind all "alt-End" "exit gnuplot"
    pause mouse close
} else {
    set terminal pngcairo size 1150, 950 font "Helvetica, 16"
    system "mkdir -p output_graph/".taskName

    set lmargin 0
    set rmargin 0
    set tmargin 0
    set bmargin 0

    do for [i = 1 : 6] {
        set output "output_graph/".taskName."/".functionBlockName[i]."_".functionNumber[i].".png"

        set xrange [A[2 * i - 1] : B[2 * i - 1]]
        set yrange [A[2 * i] : B[2 * i]]

        set title titlePng(functionBlockName[i], functionNumber[i]) font "Helvetica, 16"

        ind = 3 * (i - 1)
        splot f(i, 0, x, y) title "φ(x, y)" nosurface, \
              f(i, 1, x, y) lc rgb "orange" title "g(x, y)" nocontours, \
              trialfile index ind + 2 ls 4 lc rgb "green" lw 2 title "trial points" nocontours, \
              trialfile index ind     ls 7 lc rgb "red"   lw 6 title "X_{min}" nocontours, \
              trialfile index ind + 1 ls 7 lc rgb "blue"  lw 1 title "X" nocontours, \
              f(i, 0, x, y) with labels notitle nosurface
    }
}
