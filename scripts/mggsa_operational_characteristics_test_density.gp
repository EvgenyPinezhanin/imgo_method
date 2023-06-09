#! /usr/bin/gnuplot

taskName = "mggsa_operational_characteristics_test_density"
datafile = "output_data/".taskName.".txt"
load "output_data/".taskName."_opt.txt"

set linetype 1  lc rgb "red"         lw 5 dt 1 pt 5
set linetype 2  lc rgb "blue"        lw 5 dt 1 pt 7
set linetype 3  lc rgb "orange"      lw 5 dt 1 pt 9
set linetype 4  lc rgb "green"       lw 5 dt 1 pt 11
set linetype 5  lc rgb "dark-violet" lw 5 dt 1 pt 13

set linetype 6  lc rgb "brown"       lw 2 dt 1 pt 5
set linetype 7  lc rgb "dark-yellow" lw 2 dt 1 pt 7
set linetype 8  lc rgb "dark-violet" lw 2 dt 1 pt 9
set linetype 9  lc rgb "black"       lw 2 dt 1 pt 11
set linetype 10 lc rgb "cyan"        lw 2 dt 1 pt 13

set linetype cycle 10

set grid

set xlabel "K" font ", 15"
set ylabel "P_s(K)" font ", 15"

set tics font ", 11"

set key font ", 15" right center

title(familyName) = sprintf("Operational characteristics for mggsa on a family of tasks %s (density test)", familyName)
titlePng(familyName) = (ARG1 == 1) ? title(familyName) : sprintf("")
pointSize(i) = (i <= sizeDen / 2) ? 2 : 1

if (ARG1 == 0) {
    set title title(familyNames[int(ARG2)]) font "Helvetica Bold, 20"
    plot for [i = 1 : sizeDen] datafile index sizeDen * (ARG2 - 1) + i - 1 using 1:2 with linespoints lt i ps pointSize(i) title "m = ".den[i]

    bind all "alt-End" "exit gnuplot"
    pause mouse close
} else {
    set terminal pngcairo size 1440, 800 font "Helvetica Bold, 20"
    system "mkdir -p output_graph/".taskName

    do for [i = 1 : 4] {
        set output "output_graph/".taskName."/".taskName."_".familyNames[i].".png"

        set title titlePng(familyNames[i])
        plot for [j = 1 : sizeDen] datafile index sizeDen * (i - 1) + j - 1 using 1:2 with linespoints lt j title "m = ".den[j]
    }
}
