#!/usr/bin/gnuplot

nameTask = "mggsa_sample"

trialfile = "output_data/".nameTask.".txt"

load "output_data/".nameTask."_opt.txt"

f1(i, x, y) = (i == 0 ? 1.0 - x - y : 1/0)

f2(i, x, y) = (i == 0 ? (x - 1.0) ** 2 / 5.0 + (y - 1.0) ** 2 / 5.0 : 1/0)

g3_1(x, y) = 1.0 - x - y
g3(x, y) = (g3_1(x, y) <= 0.0) ? g3_1(x, y) : 1/0
f3(i, x, y) = (i == 0 ? x ** 2 / 5.0 + y ** 2 / 5.0 : \
               i == 1 ? g3(x, y) : 1/0)

g4_1(x, y) = (x - 2.0) ** 2 + (y - 2.0) ** 2 - 2.0
g4(x, y) = (g4_1(x, y) <= 0.0) ? g4_1(x, y) : 1/0
f4(i, x, y) = (i == 0 ? x ** 2 / 5.0 + y ** 2 / 5.0 : \
               i == 1 ? g4(x, y) : 1/0)

f(i, j, x, y) = (i == 1 ? f1(j, x, y) : \
                 i == 2 ? f2(j, x, y) : \
                 i == 3 ? f3(j, x, y) : \
                 i == 4 ? f4(j, x, y) : 1/0)

titleName(n) = sprintf("Graph of the sample function №%d, method mggsa", n)
functionName(i) = (i == 0 ? sprintf("f(x, y)") : sprintf("g%d(x, y)", i))

set contour
set view map
set cntrparam bspline levels auto 10
set cntrlabel onecolor
set cntrlabel start 5 interval 150
set cntrlabel font ",10"
set contour base

set grid

set xlabel "X"
set ylabel "Y"

set isosamples 100

if (ARG1 == 0) {
    ind = 3 * (ARG2 - 1)

    set title titleName(int(ARG2)) font "Helvetica Bold, 20"

    set xrange [AX[int(ARG2)] : BX[int(ARG2)]]
    set yrange [AY[int(ARG2)] : BY[int(ARG2)]]

    splot f(ARG2, 0, x, y) title "f(x, y)" nosurface, \
          f(ARG2, 1, x, y) lc rgb "orange" title "g(x, y)" nocontours, \
          trialfile index ind + 2 ls 4 lc rgb "green" lw 2 title "trial points" nocontours, \
          trialfile index ind     ls 7 lc rgb "red"   ps 2 title "X_{min}" nocontours, \
          trialfile index ind + 1 ls 7 lc rgb "blue"  ps 1 title "X" nocontours, \
          f(ARG2, 0, x, y) with labels notitle nosurface

    bind all "alt-End" "exit gnuplot"
    pause mouse close
} else {
    set terminal pngcairo size 950, 950
    system "mkdir -p output_graph/".nameTask

    do for [i = 1 : 4] {
        set output "output_graph/".nameTask."/".nameTask."_".i.".png"

        ind = 3 * (i - 1)

        set title titleName(i) font "Helvetica Bold, 20"

        set xrange [AX[i] : BX[i]]
        set yrange [AY[i] : BY[i]]

        splot f(i, 0, x, y) title "f(x, y)" nosurface, \
              f(i, 1, x, y) lc rgb "orange" title "g(x, y)" nocontours, \
              trialfile index ind + 2 ls 4 lc rgb "green" lw 2 title "trial points" nocontours, \
              trialfile index ind     ls 7 lc rgb "red"   lw 6 title "X_{min}" nocontours, \
              trialfile index ind + 1 ls 7 lc rgb "blue"  lw 1 title "X" nocontours, \
              f(i, 0, x, y) with labels notitle nosurface
    }
}
