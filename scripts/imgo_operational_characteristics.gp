#! /usr/bin/gnuplot

taskName = "imgo_operational_characteristics"
datafile = "output_data/".taskName.".txt"
load "output_data/".taskName."_opt.txt"

set linetype 1 lc rgb "red" lw 2
set linetype 2 lc rgb "green" lw 2
set linetype 3 lc rgb "blue" lw 2

set linetype cycle 3

set grid

set tics font "Helvetica, 16"

set xlabel "K" font "Helvetica, 16"
set ylabel "P_s(K)" font "Helvetica, 16" offset -1

set key box outside right top
set key font "Helvetica, 16" spacing 1.3

title(familyName) = sprintf("Operational characteristics of the imgo method on the family %s functions", familyName)
titlePng(familyName) = ARG1 == 1 ? title(familyName) : sprintf("")

if (ARG1 == 0) {
    set title title(familyName[ARG2 + 1]) font "Helvetica, 16"
    plot for [i = 1 : 3] datafile index (ARG2 - 1) * 3 + i - 1 using 1:2 with lines lt i title "r = ".r[(ARG2 - 1) * 3 + i]

    bind all "alt-End" "exit gnuplot"
    pause mouse close
} else {
    set terminal pngcairo size 1440, 600 font "Helvetica, 16"
    system "mkdir -p output_graph/".taskName

    set lmargin 10
    set rmargin 16
    set tmargin 3
    set bmargin 3

    do for [i = 1 : 2] {
        set output "output_graph/".taskName."/".familyName[i].".png"

        set title titlePng(familyName[i]) font "Helvetica, 19"
        plot for [j = 1 : 3] datafile index (i - 1) * 3 + j - 1 using 1:2 with lines lt j title "r = ".r[(i - 1) * 3 + j]
    }
}
