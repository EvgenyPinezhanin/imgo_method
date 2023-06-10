#! /usr/bin/gnuplot

taskName = "imgo_operational_characteristics"
datafile = "output_data/".taskName.".txt"
load "output_data/".taskName."_opt.txt"

set linetype 1 lc rgb "red"
set linetype 2 lc rgb "green"
set linetype 3 lc rgb "blue"
set linetype 4 lc rgb "orange"
set linetype 5 lc rgb "brown"
set linetype 6 lc rgb "violet"

set linetype cycle 6

set grid

set xlabel "K" font ", 13"
set ylabel "P_s(k)" font ", 13"

set tics font ", 13"

set key font ", 15"

title(familyName) = sprintf("Operational characteristics of the imgo method on the family %s functions", familyName)
titlePng(familyName) = (ARG1 == 1) ? title(familyName) : sprintf("")

if (ARG1 == 0) {
    set title title(familyNames[int(ARG2)]) font "Helvetica Bold, 20"
    plot for [i = 1 : 3] datafile index (ARG2 - 1) * 3 + i - 1 using 1:2 with lines lt i title "r = ".r[(ARG2 - 1) * 3 + i]

    bind all "alt-End" "exit gnuplot"
    pause mouse close
} else {
    set terminal pngcairo size 1440, 800 font "Helvetica Bold, 15"
    system "mkdir -p output_graph/".taskName

    do for [i = 1 : 2] {
        set output "output_graph/".taskName."/".taskName."_".familyNames[i].".png"

        set title titlePng(familyNames[i]) font "Helvetica Bold, 15"
        plot for [j = 1 : 3] datafile index (i - 1) * 3 + j - 1 using 1:2 with lines lt j title "r = ".r[(i - 1) * 3 + j]
    }
}
