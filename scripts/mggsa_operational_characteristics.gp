#! /usr/bin/gnuplot

taskName = "mggsa_operational_characteristics"
datafile = "output_data/".taskName."/operational_characteristics.txt"
load "output_data/".taskName."/vars.txt"

set linetype 1 lc rgb "red" lw 2
set linetype 2 lc rgb "green" lw 2
set linetype 3 lc rgb "blue" lw 2
set linetype 4 lc rgb "orange" lw 2

set linetype cycle 4

set grid

set tics font "Helvetica, 16"

set xlabel "K" font "Helvetica, 16"
set ylabel "P_s(K)" font "Helvetica, 16" offset -1

set key box outside right top
set key font "Helvetica, 16" spacing 1.3

title(familyName) = sprintf("Operational characteristics for DIRECT and mggsa on the family %s functions", familyName)
titlePng(familyName) = ARG1 == 1 ? title(familyName) : sprintf("")

if (ARG1 == 0) {
    set title title(familyName[ARG2 + 1]) font "Helvetica, 16"

    set lmargin 12

    ind = numberKey * ARG2
    plot for [i = 1 : numberKey] datafile index ind + i - 1 using 1:2 with lines lt i title "r = ".r[ind + i].", key = ".key[i]

    bind all "alt-End" "exit gnuplot"
    pause mouse close
} else {
    // set terminal pngcairo size 1440, 600 font "Helvetica, 16"
    // system "mkdir -p output_graph/".taskName
// 
    // set key width -3
// 
    // set lmargin 10
    // set rmargin 24
    // set tmargin 3
    // set bmargin 3
// 
    // do for [i = 1 : 4] {
    //     set output "output_graph/".taskName."/".familyName[i].".png"
// 
    //     ind = numberKey * (i - 1)
    //     set title titlePng(familyName[i]) font "Helvetica, 16"
    //     plot for [j = 1 : numberKey] datafile index ind + j - 1 using 1:2 with lines lt j title "r = ".r[ind + j].", key = ".key[j]
    // }
}
