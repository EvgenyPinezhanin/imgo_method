#! /usr/bin/gnuplot

sampleName = "fitting_family_problems_operational_characteristics"

load "output_data/new_model_class/".sampleName."/vars.txt"

set linetype 1  lc rgb "red"         lw 2
set linetype 2  lc rgb "green"       lw 2
set linetype 3  lc rgb "blue"        lw 2
set linetype 4  lc rgb "orange"      lw 2
set linetype 5  lc rgb "brown"       lw 2
set linetype 6  lc rgb "dark-yellow" lw 2
set linetype 7  lc rgb "dark-violet" lw 2
set linetype 8  lc rgb "cyan"        lw 2

set linetype cycle 8

fontName = "Helvetica, 16"

set grid

set xlabel "K" font fontName
set ylabel "P_s(K)" font fontName offset -1

set tics font fontName

set key box outside right top
set key font fontName spacing 1.3

title(familyName) = sprintf("Operational characteristics on the %s", familyName)
titlePng(familyName) = ARG1 == 1 ? title(familyName) : sprintf("")
dataFile(number) = sprintf("output_data/new_model_class/%s/%s/%s_%s_%s", sampleName, methodNames[1], familyName, key[number], r[number])

if (ARG1 == 0) {
    set lmargin 10

    set title title(familyName) font fontName
    plot for [i = 1 : numberKey] dataFile(i) using 1:2 with lines lt i title "r = ".r[i].", key = ".key[i]

    bind all "alt-End" "exit gnuplot"
    pause mouse close
} else {
    set terminal pngcairo size 1440, 600 font "Helvetica, 16"
    system "mkdir -p output_graph/".sampleName

    set lmargin 10
    set rmargin 16
    set tmargin 3
    set bmargin 3

    do for [i = 1 : 2] {
        set output "output_graph/".sampleName."/".familyName[i].".png"

        set title titlePng(familyName[i]) font "Helvetica, 19"
        plot for [j = 1 : 3] datafile index (i - 1) * 3 + j - 1 using 1:2 with lines lt j title "r = ".r[(i - 1) * 3 + j]
    }
}
