#! /usr/bin/gnuplot

taskName = "mggsa_operational_characteristics_test_incr"
datafile = "output_data/".taskName.".txt"
load "output_data/".taskName."_opt.txt"

set linetype 1 lc rgb "red"
set linetype 2 lc rgb "green"
set linetype 3 lc rgb "blue"
set linetype 4 lc rgb "orange"

set linetype cycle 4

set grid

set tics font "Helvetica, 16"

set xlabel "increment" font "Helvetica, 16"

set tics font ", 16"

set key box outside right top
set key font "Helvetica, 16" spacing 1.3

title(familyName) = sprintf("Operational characteristics for mggsa on a family of tasks %s, (incr test)", familyName)
titlePng(familyName) = ARG1 == 1 ? title(familyName) : sprintf("")

array yName[3]
yName[1] = "trials"
yName[2] = "point trials"
yName[3] = "point trials / trials"

array fileName[3]
fileName[1] = "trials"
fileName[2] = "point_trials"
fileName[3] = "point_trials_trials"

if (ARG1 == 0) {
    set title title(familyName[ARG3 + 1]) font "Helvetica, 16"

    set ylabel yName[ARG2 + 1] font "Helvetica, 16"

    set lmargin 12

    ind = numberKey * ARG3
    if (ARG2 == 2) {
        plot for [i = 1 : numberKey] datafile index ind + i - 1 using 1:($3 / $2) with lines lt i title "r = ".r[ind + i]
    } else {
        plot for [i = 1 : numberKey] datafile index ind + i - 1 using 1:(column(ARG2 + 2)) with lines lt i title "r = ".r[ind + i]
    }

    bind all "alt-End" "exit gnuplot"
    pause mouse close
} else {
    set terminal pngcairo size 1440, 700 font "Helvetica, 16"
    system "mkdir -p output_graph/".taskName

    set lmargin 10
    set rmargin 16
    set tmargin 3
    set bmargin 3

    do for [i = 1 : 4] {
        set title titlePng(familyName[i])
        ind = numberKey * (i - 1)

		do for [j = 1 : 3] {
			set output "output_graph/".taskName."/".familyName[i]."_".fileName[j].".png"

            set ylabel yName[j] font "Helvetica, 16"

			if (j == 3) {
    		    plot for [k = 1 : numberKey] datafile index ind + k - 1 using 1:($3 / $2) with lines lt k title "r = ".r[ind + k]
    		} else {
    		    plot for [k = 1 : numberKey] datafile index ind + k - 1 using 1:(column(j + 1)) with lines lt k title "r = ".r[ind + k]
			}
		}
    }
}
