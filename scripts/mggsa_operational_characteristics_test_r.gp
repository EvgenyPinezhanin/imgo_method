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

set xlabel "r" font ", 15"
set ylabel "P_s^{max}" font ", 15"

set tics font ", 11"
set key font ", 15"

title(familyName) = sprintf("Operational characteristics for mggsa on a family of tasks %s (r test)", familyName)
titlePng(familyName) = (ARG1 == 1) ? title(familyName) : sprintf("")

array Parameters[2]
Parameters[1] = "key = 1"
Parameters[2] = "key = 3, incr = 20"

array y2Name[3]
yName[1] = "trials"
yName[2] = "points"
yName[3] = "points / trials"

array fileName[3]
fileName[1] = "trials"
fileName[2] = "points"
fileName[3] = "points_trials"

if (ARG1 == 0) {
    ind = 2 * (ARG3 - 1)

    set title title(familyNames[int(ARG3)]) font "Helvetica Bold, 20"

    if (ARG2 == 0) {
        plot for [i = 1 : 2] datafile index ind + i - 1 using 1:2 with lines lt i * 2 - 1 title Parameters[i]
    } else {
        set y2tics

        set y2label y2Name[int(ARG2)] font ", 15"

        if (ARG1 == 1) {
            plot for [i = 1:2] datafile index ind + i - 1 using 1:2 with lines lt i * 2 - 1 title "P_s^{max}(".Parameters[i].")" axes x1y1, \
                 for [i = 1:2] datafile index ind + i - 1 using 1:3 with lines lt i * 2 title "trials(".Parameters[i].")" axes x1y2, \
                 datafile index ind using 1:(K[ARG2 + 1]) with lines lt 5 title "K^{max}" axes x1y2
        }
        if (ARG1 == 2) {
            plot for [i = 1:2] datafile index ind + i - 1 using 1:2 with lines lt i * 2 - 1 title "P_s^{max}(".Parameters[i].")" axes x1y1, \
                 for [i = 1:2] datafile index ind + i - 1 using 1:4 with lines lt i * 2 title "points(".Parameters[i].")" axes x1y2
        }
        if (ARG1 == 3) {
            plot for [i = 1:2] datafile index ind + i - 1 using 1:2 with lines lt i * 2 - 1 title "P_s^{max}(".Parameters[i].")" axes x1y1, \
                 for [i = 1:2] datafile index ind + i - 1 using 1:($4 / $3) with lines lt i * 2 title "points / trials(".Parameters[i].")" axes x1y2
        }
    }

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
