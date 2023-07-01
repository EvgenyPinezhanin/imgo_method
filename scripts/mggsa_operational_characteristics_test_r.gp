#! /usr/bin/gnuplot

taskName = "mggsa_operational_characteristics_test_r"
datafile = "output_data/".taskName.".txt"
load "output_data/".taskName."_opt.txt"

set linetype 1 lc rgb "red"    lw 3
set linetype 2 lc rgb "green"  lw 3
set linetype 3 lc rgb "orange" lw 2
set linetype 4 lc rgb "black"  lw 2
set linetype 5 lc rgb "blue"   lw 1

set linetype cycle 5

set grid

set xlabel "r" font "Helvetica, 16"
set ylabel "P_s^{max}" font "Helvetica, 16"

set tics font ", 16"

set key box outside right top
set key font "Helvetica, 16" spacing 1.3

title(familyName) = sprintf("Operational characteristics for mggsa on a family of tasks %s (r test)", familyName)
titlePng(familyName) = ARG1 == 1 ? title(familyName) : sprintf("")

array parameters[2]
parameters[1] = "key = 1"
parameters[2] = "key = 3, incr = 30"

array y2Name[3]
y2Name[1] = "trials"
y2Name[2] = "points"
y2Name[3] = "points / trials"

array fileName[4]
fileName[1] = "p"
fileName[2] = "trials"
fileName[3] = "points"
fileName[4] = "points_del_trials"

if (ARG1 == 0) {
    set title title(familyName[ARG3 + 1]) font "Helvetica, 16"

    set lmargin 12

    ind = 2 * ARG3
    if (ARG2 == 0) {
        plot for [i = 1 : 2] datafile index ind + i - 1 using 1:2 with lines lt i * 2 - 1 title parameters[i]
    } else {
        set y2tics

        set y2label y2Name[ARG2 + 1] font "Helvetica, 16"

        if (ARG2 == 1) {
            plot for [i = 1 : 2] datafile index ind + i - 1 using 1:2 with lines lt i * 2 - 1 title "P_s^{max}(".parameters[i].")" axes x1y1, \
                 for [i = 1 : 2] datafile index ind + i - 1 using 1:3 with lines lt i * 2 title "trials(".parameters[i].")" axes x1y2, \
                 datafile index ind using 1:(K[ARG3 + 1]) with lines lt 5 title "K^{max}" axes x1y2
        }
        if (ARG2 == 2) {
            plot for [i = 1 : 2] datafile index ind + i - 1 using 1:2 with lines lt i * 2 - 1 title "P_s^{max}(".parameters[i].")" axes x1y1, \
                 for [i = 1 : 2] datafile index ind + i - 1 using 1:4 with lines lt i * 2 title "points(".parameters[i].")" axes x1y2
        }
        if (ARG2 == 3) {
            plot for [i = 1 : 2] datafile index ind + i - 1 using 1:2 with lines lt i * 2 - 1 title "P_s^{max}(".parameters[i].")" axes x1y1, \
                 for [i = 1 : 2] datafile index ind + i - 1 using 1:($4 / $3) with lines lt i * 2 title "points / trials(".parameters[i].")" axes x1y2
        }
    }

    bind all "alt-End" "exit gnuplot"
    pause mouse close
} else {
    set terminal pngcairo size 1440, 800 font "Helvetica, 16"
    system "mkdir -p output_graph/".taskName

    set lmargin 10
    set rmargin 16
    set tmargin 3
    set bmargin 3

    do for [i = 1 : 4] {
        set title titlePng(familyName[i])

        ind = 2 * (i - 1)

        set output "output_graph/".taskName."/".familyName[i]."_".fileName[1].".png"
        plot for [i = 1 : 2] datafile index ind + i - 1 using 1:2 with lines lt i * 2 - 1 title parameters[i]

        set y2tics

        set y2label y2Name[ARG2 + 1] font "Helvetica, 16"

        set output "output_graph/".taskName."/".familyName[i]."_".fileName[2].".png"
        plot for [i = 1 : 2] datafile index ind + i - 1 using 1:2 with lines lt i * 2 - 1 title "P_s^{max}(".parameters[i].")" axes x1y1, \
             for [i = 1 : 2] datafile index ind + i - 1 using 1:3 with lines lt i * 2 title "trials(".parameters[i].")" axes x1y2, \
             datafile index ind using 1:(K[ARG3 + 1]) with lines lt 5 title "K^{max}" axes x1y2

        set output "output_graph/".taskName."/".familyName[i]."_".fileName[3].".png"
        plot for [i = 1 : 2] datafile index ind + i - 1 using 1:2 with lines lt i * 2 - 1 title "P_s^{max}(".parameters[i].")" axes x1y1, \
             for [i = 1 : 2] datafile index ind + i - 1 using 1:4 with lines lt i * 2 title "points(".parameters[i].")" axes x1y2

        set output "output_graph/".taskName."/".familyName[i]."_".fileName[4].".png"
        plot for [i = 1 : 2] datafile index ind + i - 1 using 1:2 with lines lt i * 2 - 1 title "P_s^{max}(".parameters[i].")" axes x1y1, \
             for [i = 1 : 2] datafile index ind + i - 1 using 1:($4 / $3) with lines lt i * 2 title "points / trials(".parameters[i].")" axes x1y2
    }
}
