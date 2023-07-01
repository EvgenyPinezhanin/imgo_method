#! /usr/bin/gnuplot

taskName = "mggsa_operational_characteristics_test_density"
datafile = "output_data/".taskName.".txt"
load "output_data/".taskName."_opt.txt"

set linetype 1  lc rgb "red"         lw 7 dt 1 pt 5
set linetype 2  lc rgb "blue"        lw 7 dt 1 pt 7
set linetype 3  lc rgb "orange"      lw 7 dt 1 pt 9
set linetype 4  lc rgb "green"       lw 7 dt 1 pt 11
set linetype 5  lc rgb "dark-violet" lw 7 dt 1 pt 13

set linetype 6  lc rgb "brown"       lw 2 dt 1 pt 5
set linetype 7  lc rgb "dark-yellow" lw 2 dt 1 pt 7
set linetype 8  lc rgb "dark-violet" lw 2 dt 1 pt 9
set linetype 9  lc rgb "black"       lw 2 dt 1 pt 11
set linetype 10 lc rgb "cyan"        lw 2 dt 1 pt 13

set linetype cycle 10

set grid

set tics font "Helvetica, 16"

set xlabel "K" font "Helvetica, 16"
set ylabel "P_s(K)" font "Helvetica, 16"

set key box outside right top
set key font "Helvetica, 16" spacing 1.3

title(familyName) = sprintf("Operational characteristics for mggsa on a family of tasks %s (density test)", familyName)
titlePng(familyName) = ARG1 == 1 ? title(familyName) : sprintf("")
pointSize(i) = i <= sizeDen / 2 ? 2 : 1

if (ARG1 == 0) {
    set title title(familyName[ARG2 + 1]) font "Helvetica Bold, 20"

    set lmargin 12

    ind = sizeDen * ARG2
    plot for [i = 1 : sizeDen] datafile index ind + i - 1 using 1:2 with linespoints lt i ps pointSize(i) title "m = ".den[i]

    bind all "alt-End" "exit gnuplot"
    pause mouse close
} else {
    set terminal pngcairo size 1440, 700 font "Helvetica, 20"
    system "mkdir -p output_graph/".taskName

    set key width -3

    set lmargin 8
    set rmargin 26
    set tmargin 3
    set bmargin 3

    do for [i = 1 : 4] {
        set output "output_graph/".taskName."/".familyName[i].".png"

        set title titlePng(familyName[i])

        ind = sizeDen * (i - 1)
        plot for [j = 1 : sizeDen] datafile index sizeDen * (i - 1) + j - 1 using 1:2 with linespoints lt j title "m = ".den[j]
    }
}
