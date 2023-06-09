#! /usr/bin/gnuplot

taskName = "mggsa_operational_characteristics_test_key"
datafile = "output_data/".taskName.".txt"
load "output_data/".taskName."_opt.txt"

set linetype 1 lc rgb "red"
set linetype 2 lc rgb "green"
set linetype 3 lc rgb "blue"
set linetype 4 lc rgb "orange"

set linetype cycle 4

set grid

set xlabel "K" font ", 15"
set ylabel "P_s(K)" font ", 15"

set tics font ", 11"

set key font ", 15" center right

title(familyName) = sprintf("Operational characteristics for mggsa on a family of tasks ".familyNames[int(ARG2)]." (key test)", familyName)
titlePng(familyName) = (ARG1 == 1) ? title(familyName) : sprintf("")

if (ARG1 == 0) {
    ind = numberKey * (ARG2 - 1)

    set title title(familyNames[int(ARG2)]) font "Helvetica Bold, 20"
    plot for [i = 1 : numberKey] datafile index ind + i - 1 using 1:2 with lines lt i title "r = ".r[ind + i].", key = ".key[i]

    bind all "alt-End" "exit gnuplot"
    pause mouse close
} else {
    set terminal pngcairo size 1440, 800 font "Helvetica Bold, 20"
    system "mkdir -p output_graph/".taskName

    do for [i = 1 : 4] {
        set output "output_graph/".taskName."/".taskName."_".familyNames[i].".png"

        ind = numberKey * (i - 1)

        set title titlePng(familyNames[i])
        plot for [j = 1 : numberKey] datafile index ind + j - 1 using 1:2 with lines lt j title "r = ".r[ind + j].", key = ".key[j]
    }
}
