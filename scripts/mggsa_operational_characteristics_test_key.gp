#! /usr/bin/gnuplot

taskName = "mggsa_operational_characteristics_test_key"
datafile = "output_data/".taskName.".txt"
load "output_data/".taskName."_opt.txt"

set linetype 1 lc rgb "red" lw 2
set linetype 2 lc rgb "green" lw 2
set linetype 3 lc rgb "blue" lw 2
set linetype 4 lc rgb "orange" lw 2

set linetype cycle 4

set grid

set xlabel "K"
set ylabel "P_s(K)"

set key spacing 1.3 center right

title(familyName) = sprintf("Operational characteristics for mggsa on a family of tasks ".familyNames[int(ARG2)]." (key test)", familyName)
titlePng(familyName) = (ARG1 == 1) ? title(familyName) : sprintf("")

if (ARG1 == 0) {
    set xlabel font ", 15"
    set ylabel font ", 15"

    set tics font ", 11"

    set key spacing 1.3 center right

    ind = numberKey * (ARG2 - 1)

    set title title(familyNames[int(ARG2)]) font "Helvetica Bold, 20"
    plot for [i = 1 : numberKey] datafile index ind + i - 1 using 1:2 with lines lt i title "r = ".r[ind + i].", key = ".key[i]

    bind all "alt-End" "exit gnuplot"
    pause mouse close
} else {
    set terminal pngcairo size 1440, 700 font "Helvetica, 20"
    system "mkdir -p output_graph/".taskName

    do for [i = 1 : 4] {
        set output "output_graph/".taskName."/".familyNames[i].".png"

        ind = numberKey * (i - 1)

        set title titlePng(familyNames[i]) font "Helvetica, 20"
        plot for [j = 1 : numberKey] datafile index ind + j - 1 using 1:2 with lines lt j title "r = ".r[ind + j].", key = ".key[j]
    }
}
