#! /usr/bin/gnuplot

nameTask = "imgo_operational_characteristics"

datafile = "output_data/".nameTask.".txt"

load "output_data/".nameTask."_opt.txt"

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

if (ARG1 == 0) {
    if (ARG2 == 3) {
        set title "Comparison of the operational characteristics of the imgo method on the families of Shekel and Hill functions" \
                  font "Helvetica Bold, 25"
        plot for [i = 1 : 6] datafile index i - 1 using 1:2 with lines lt i title "r = ".r[i]."(".familyNames[(i - 1) / 3 + 1].")"
    } else {
        set title "Operational characteristics of the imgo method on the family ".familyNames[int(ARG2)]." functions" \
                  font "Helvetica Bold, 25"
        plot for [i = 1 : 3] datafile index (ARG2 - 1) * 3 + i - 1 using 1:2 with lines lt i title "r = ".r[(ARG2 - 1) * 3 + i]
    }

    bind all "alt-End" "exit gnuplot"
    pause mouse close
} else {
    set terminal pngcairo size 1280, 800
    system "mkdir -p output_graph/".nameTask

    do for [i = 1 : 2] {
        set output "output_graph/".nameTask."/".nameTask."_".familyNames[i].".png"

        set title "Operational characteristics of the imgo method on the family ".familyNames[i]." functions" \
                  font "Helvetica Bold, 20"
        plot for [j = 1 : 3] datafile index (i - 1) * 3 + j - 1 using 1:2 with lines lt j title "r = ".r[(i - 1) * 3 + j]
    }

    set output "output_graph/".nameTask."/".nameTask."_comparison.png"

    set title "Comparison of the operational characteristics of the imgo method on the families of Shekel and Hill functions" \
              font "Helvetica Bold, 17"
    plot for [i = 1 : 6] datafile index i - 1 using 1:2 with lines lt i title "r = ".r[i]."(".familyNames[(i - 1) / 3 + 1].")"
}
