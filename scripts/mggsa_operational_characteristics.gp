#! /usr/bin/gnuplot

nameTask = "mggsa_operational_characteristics"

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

set xlabel "K" font ", 15"
set ylabel "P_s(k)" font ", 15"

set tics font ", 11"
set key font ", 15"
set key center right

if (ARG1 == 0) {
    if (ARG2 < 4) {
        set title "Operational characteristics of the imgo method on the family ".familyNames[int(ARG2)]." functions" \
                  font "Helvetica Bold, 20"
        plot for [i = 1 : 3] datafile index 3 * (ARG2 - 1) + i - 1 using 1:2 with lines lt i title "r = ".r[3 * (ARG2 - 1) + i]
    }
    if (ARG2 == 4) {
        set title "Comparison of the operational characteristics of the imgo method on the families of Grishagin and GKLS functions" \
                  font "Helvetica Bold, 20"
        plot for [i = 1 : 6] datafile index i - 1 using 1:2 with lines lt i title "r = ".r[i]."(".familyNames[(i - 1) / 3 + 1].")"
    }
    if (ARG2 == 5) {
        set title "Comparison of the operational characteristics of the imgo method on the families of Grishagin and GKLS functions (constrained)" \
                  font "Helvetica Bold, 20"
        plot for [i = 7 : 12] datafile index i - 1 using 1:2 with lines lt i - 6 title "r = ".r[i]."(".familyNames[(i - 1) / 3 + 1].")"
    }

    bind all "alt-End" "exit gnuplot"
    pause mouse close
} else {
    set terminal pngcairo size 1280, 800
    system "mkdir -p output_graph/".nameTask

    do for [i = 1 : 4] {
        set output "output_graph/".nameTask."/".nameTask."_".familyNames[i].".png"

        set title "Operational characteristics of the imgo method on the family ".familyNames[int(ARG2)]." functions" \
                  font "Helvetica Bold, 20"
        plot for [j = 1 : 3] datafile index 3 * (i - 1) + j - 1 using 1:2 with lines lt j title "r = ".r[3 * (i - 1) + j]
    }

    set output "output_graph/".nameTask."/".nameTask."_Grishagin_GKLS.png"

    set title "Comparison of the operational characteristics of the imgo method on the families of Grishagin and GKLS functions" \
              font "Helvetica Bold, 15"
    plot for [i = 1 : 6] datafile index i - 1 using 1:2 with lines lt i title "r = ".r[i]."(".familyNames[(i - 1) / 3 + 1].")"

    set output "output_graph/".nameTask."/".nameTask."_Grishagin_GKLS_constrained.png"

    set title "Comparison of the operational characteristics of the imgo method on the families of Grishagin and GKLS functions (constrained)" \
              font "Helvetica Bold, 14"
    plot for [i = 7 : 12] datafile index i - 1 using 1:2 with lines lt i - 6 title "r = ".r[i]."(".familyNames[(i - 1) / 3 + 1].")"
}
