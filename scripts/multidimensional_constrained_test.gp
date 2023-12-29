#!/usr/bin/gnuplot

load "scripts/functions/multidimensional_constrained.txt"

fontName = "Helvetica, 16"

set isosamples 120

set view map

set grid
set contour base
set cntrlabel onecolor
set cntrparam levels auto 10
set cntrlabel start 5 interval 40
set cntrlabel font "Helvetica, 10"

set xlabel "X" font fontName
set ylabel "Y" font fontName offset -1

set tics font fontName

set key box outside right top
set key font fontName spacing 1.5

title(block, n, method) = sprintf("Graph of the %s function №%d, %s method", block, n, method)
titlePng(block, n, method) = ARG2 == 1 ? title(block, n, method) : sprintf("")

if (ARG2 == 0) {
    trialsFile = "output_data/".ARG1."_test/".blockNames[ARG3 + 1]."_".(ARG4 + 1)
    exist = system("[ ! -f ".trialsFile." ]; echo $?")
    if (exist == 1) {
        set key width 3

        set xrange [A(ARG3 + 1, ARG4 + 1, 1) : B(ARG3 + 1, ARG4 + 1, 1)]
        set yrange [A(ARG3 + 1, ARG4 + 1, 2) : B(ARG3 + 1, ARG4 + 1, 2)]

        set title title(blockNames[ARG3 + 1], ARG4 + 1, ARG1) font fontName

        splot f(ARG3 + 1, ARG4 + 1, x, y, 0) title "φ(x, y)" nosurface, \
              f(ARG3 + 1, ARG4 + 1, x, y, 1) lc rgb "orange" title "g(x, y)" nocontours, \
              trialsFile index 2 ls 4 lc rgb "green" lw 2 title "trial points" nocontours, \
              trialsFile index 0 ls 7 lc rgb "red"   ps 2 title "X_{min}" nocontours, \
              trialsFile index 1 ls 7 lc rgb "blue"  ps 1 title "X" nocontours, \
              f(ARG3 + 1, ARG4 + 1, x, y, 0) with labels notitle nosurface

        bind all "alt-End" "exit gnuplot"
        pause mouse close
    } else {
        print("TrialsFile ".trialsFile." not exist")
    }
} else {
    set terminal pngcairo size 1150, 950 font "Helvetica, 16"
    system "mkdir -p output_graph/".ARG1."_test"

    set lmargin 1
    set rmargin 0
    set tmargin 0
    set bmargin 0

    do for [i = 1 : numberBlocks] {
        do for [j = 1 : numberFunctions[i]] {
            trialsFile = "output_data/".ARG1."_test/".blockNames[i]."_".j
            exist = system("[ ! -f ".trialsFile." ]; echo $?")
            if (exist == 1) {
                set output "output_graph/".ARG1."_test/".blockNames[i]."_".j.".png"

                set xrange [A(i, j, 1) : B(i, j, 1)]
                set yrange [A(i, j, 2) : B(i, j, 2)]

                set title titlePng(blockNames[i], j, ARG1) font fontName

                splot f(i, j, x, y, 0) title "φ(x, y)" nosurface, \
                      f(i, j, x, y, 1) lc rgb "orange" title "g(x, y)" nocontours, \
                      trialsFile index 2 ls 4 lc rgb "green" lw 2 title "trial points" nocontours, \
                      trialsFile index 0 ls 7 lc rgb "red"   lw 6 title "X_{min}" nocontours, \
                      trialsFile index 1 ls 7 lc rgb "blue"  lw 1 title "X" nocontours, \
                      f(i, j, x, y, 0) with labels notitle nosurface
            } else {
                print("TrialsFile ".trialsFile." not exist")
            }
        }
    }
}
