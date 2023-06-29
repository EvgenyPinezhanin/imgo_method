#!/usr/bin/gnuplot

taskName = "mggsa_test"
trialfile = "output_data/".taskName.".txt"
load "output_data/".taskName."_opt.txt"

g1_1(x, y) = 0.01 * ((x - 2.2) ** 2 + (y - 1.2) ** 2 - 2.25)
g1_2(x, y) = 100.0 * (1.0 - ((x - 2.0) ** 2) / 1.44 - (0.5 * y) ** 2)
g1_3(x, y) = 10.0 * (y - 1.5 - 1.5 * sin(6.283 * (x - 1.75)))
g1(x, y) = (g1_1(x, y) <= 0.0) ? (g1_2(x, y) <= 0 ? (g1_3(x, y) <= 0 ? g1_3(x, y) : 1 / 0) : 1 / 0) : 1 / 0
f1(i, x, y) = (i == 0 ? -1.5 * x ** 2 * exp(1 - x ** 2 - 20.25 * (x - y) ** 2) - (0.5 * (x - 1) * (y - 1)) ** 4 * \
                        exp(2.0 - (0.5 * (x - 1)) ** 4 - (y - 1.0) ** 4) : \
               i == 1 ? g1(x, y) : 1/0)

array C[20]
C[1] =   75.1963666677
C[2] =   -3.8112755343
C[3] =    0.1269366345
C[4] =   -0.0020567665
C[5] =    0.0000103450
C[6] =   -6.8306567631
C[7] =    0.0302344793
C[8] =   -0.0012813448
C[9] =    0.0000352559
C[10] =  -0.0000002266
C[11] =   0.2564581253
C[12] =  -0.0034604030
C[13] =   0.0000135139
C[14] = -28.1064434908
C[15] =  -0.0000052375
C[16] =  -0.0000000063
C[17] =   0.0000000007
C[18] =   0.0003405462
C[19] =  -0.0000016638
C[20] =  -2.8673112392

g2_1(x, y) = 450.0 - x * y
g2_2(x, y) = (0.1 * x - 1.0) ** 2 - y
g2_3(x, y) = 8.0 * (x - 40.0) - (y - 30.0) * (y - 55.0)
g2_4(x, y) = y + (x - 35.0) * (x - 30.0) / 125.0 - 80.0
g2(x, y) = (g2_1(x, y) <= 0.0) ? (g2_2(x, y) <= 0 ? (g2_3(x, y) <= 0 ? (g2_4(x, y) <= 0 ? (g2_4(x, y)) : 1 / 0) : 1 / 0) : 1 / 0) : 1 / 0
f2(i, x, y) = (i == 0 ? -(C[1] + C[2] * x + C[3] * x ** 2 + C[4] * x ** 3 + C[5] * x ** 4 + C[6] * y + C[7] * x * y + C[8] * x ** 2 * y + \
                          C[9] * x ** 3 * y + C[10] * x ** 4 * y + C[11] * y ** 2 + C[12] * y ** 3 + C[13] * y ** 4 + C[14] / (y + 1) + \
                          C[15] * x ** 2 * y ** 2 + C[16] * x ** 3 * y ** 2 + C[17] * x ** 3 * y ** 3 + C[18] * x * y ** 2 + \
                          C[19] * x * y ** 3 + C[20] * exp(0.0005 * x * y)) : \
               i == 1 ? g2(x, y) : 1/0)

f(i, j, x, y) = (i == 1 ? f1(j, x, y) : \
                 i == 2 ? f2(j, x, y) : 1/0)

set contour base
set cntrlabel onecolor
set cntrlabel start 5 interval 150
set cntrlabel font ", 10"

set view map

set isosamples 120

set grid

set xlabel "X"
set ylabel "Y"

set key left bottom width -3
set key box opaque
set key spacing 1.3

title(n) = sprintf("Graph of the test function №%d, method mggsa", n)
titlePng(n) = (ARG1 == 1) ? title(n) : sprintf("")
functionName(i) = (i == 0 ? sprintf("f(x, y)") : sprintf("g%d(x, y)", i))

if (ARG1 == 0) {
    set tics font ", 13"

    set key font ", 15"

    ind = 3 * (ARG2 - 1)

    set title title(int(ARG2)) font "Helvetica, 19"

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
    set terminal pngcairo size 950, 950 font "Helvetica, 19"
    system "mkdir -p output_graph/".taskName

    set lmargin 1
    set rmargin 0
    set tmargin 0
    set bmargin 0

    do for [i = 1 : 2] {
        set output "output_graph/".taskName."/".taskName."_".i.".png"

        set cntrparam bspline levels incremental incrementalParam[3 * (i - 1) + 1], \
                                                 incrementalParam[3 * (i - 1) + 2], \
                                                 incrementalParam[3 * (i - 1) + 3]

        ind = 3 * (i - 1)

        set title titlePng(i) font "Helvetica, 19"

        set xrange [AX[i] : BX[i]]
        set yrange [AY[i] : BY[i]]

        splot f(i, 0, x, y) title "φ(x, y)" nosurface, \
              f(i, 1, x, y) lc rgb "orange" title "g(x, y)" nocontours, \
              trialfile index ind     ls 7 lc rgb "red"   lw 6 title "X_{min}" nocontours #, \
              # trialfile index ind + 2 ls 4 lc rgb "green" lw 2 title "trial points" nocontours, \
              # trialfile index ind + 1 ls 7 lc rgb "blue"  lw 1 title "X" nocontours #, \
              # f(i, 0, x, y) with labels notitle nosurface
    }
}
