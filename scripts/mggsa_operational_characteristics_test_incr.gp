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

set xlabel "increment" font ", 15"

set tics font ", 11"

set key font ", 15"

array yName[3]
yName[1] = "trials"
yName[2] = "point trials"
yName[3] = "point trials / trials"

array fileName[3]
fileName[1] = "trials"
fileName[2] = "point_trials"
fileName[3] = "point_trials_trials"

title(familyName) = sprintf("Operational characteristics for mggsa on a family of tasks %s, (incr test)", familyName)
titlePng(familyName) = (ARG1 == 1) ? title(familyName) : sprintf("")

if (ARG1 == 0) {
    ind = numberKey * (ARG3 - 1)

    set title title(familyNames[int(ARG3)]) font "Helvetica Bold, 20"

    set ylabel yName[ARG2 + 1] font ", 15"

    if (ARG2 == 2) {
        plot for [i = 1 : numberKey] datafile index ind + i - 1 using 1:($3 / $2) with lines lt i title "r = ".r[ind + i]
    } else {
        plot for [i = 1 : numberKey] datafile index ind + i - 1 using 1:(column(ARG2 + 2)) with lines lt i title "r = ".r[ind + i]
    }

    bind all "alt-End" "exit gnuplot"
    pause mouse close
} else {
    set terminal pngcairo size 1440, 800 font "Helvetica Bold, 20"
    system "mkdir -p output_graph/".taskName

    do for [i = 1 : 4] {
        set title titlePng(familyNames[i])
        ind = numberKey * (i - 1)

		do for [j = 1 : 3] {
			set output "output_graph/".taskName."/".taskName."_".familyNames[i]."_".fileName[j].".png"

            set ylabel yName[j] font ", 15"

			if (j == 3) {
    		    plot for [k = 1 : numberKey] datafile index ind + k - 1 using 1:($3 / $2) with lines lt k title "r = ".r[ind + k]
    		} else {
    		    plot for [k = 1 : numberKey] datafile index ind + k - 1 using 1:(column(j + 1)) with lines lt k title "r = ".r[ind + k]
			}
		}
    }
}
