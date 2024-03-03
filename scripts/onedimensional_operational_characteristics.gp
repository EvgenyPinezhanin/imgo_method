#! /usr/bin/gnuplot

sampleName = "onedimensional_operational_characteristics"

load "output_data/".sampleName."/vars.txt"

set linetype 1  lc rgb "red"         lw 2
set linetype 2  lc rgb "green"       lw 2
set linetype 3  lc rgb "blue"        lw 2
set linetype 4  lc rgb "orange"      lw 2
set linetype 5  lc rgb "brown"       lw 2
set linetype 6  lc rgb "dark-yellow" lw 2
set linetype 7  lc rgb "dark-violet" lw 2

set linetype cycle 7

fontName = "Helvetica, 16"

set grid

set xlabel "K" font fontName
set ylabel "P_s(K)" font fontName offset -1

set tics font fontName

set key box outside right top
set key font fontName spacing 1.3

title(familyName) = sprintf("Operational characteristics on the %s", familyName)
titlePng(familyName) = ARG1 == 1 ? title(familyName) : sprintf("")
dataFile(number) = number == 1 ? sprintf("output_data/%s/%s/%s", sampleName, methodNames[1], familyName[ARG2 + 1]) : \
                   number <= 4 ? sprintf("output_data/%s/%s/%s_%s", sampleName, methodNames[2], familyName[ARG2 + 1], r[number - 1]) : \
                                 sprintf("output_data/%s/%s/%s_%s", sampleName, methodNames[3], familyName[ARG2 + 1], r[number - 1])

if (ARG1 == 0) {
    set title title(familyName[ARG2 + 1]) font fontName
    plot for [i = 1 : 7] dataFile(i) using 1:2 with lines lt i title "r = "

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
