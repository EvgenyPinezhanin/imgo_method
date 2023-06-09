#! /usr/bin/gnuplot

taskName = "mggsa_test_incr"
datafile = "output_data/".taskName.".txt"
load "output_data/".taskName."_opt.txt"

set linetype 1 lc rgb "#0000FF" lw 5 pt 1 dt 1
set linetype 2 lc rgb "#0020FF" lw 4 pt 1 dt 2
set linetype 3 lc rgb "#0040FF" lw 3 pt 1 dt 3
set linetype 4 lc rgb "#0060FF" lw 2 pt 1 dt 4
set linetype 5 lc rgb "#0080FF" lw 1 pt 1 dt 5

set linetype cycle 5

set xlabel "increment"

array yName[4]
yName[1] = "trials"
yName[2] = "point trials"
yName[3] = "accuracy"
yName[4] = "point trials / trials"

array fileName[4]
fileName[1] = "trials"
fileName[2] = "point_trials"
fileName[3] = "accuracy"
fileName[4] = "point_trials_trials"


if (ARG1 == 0) {
    set xrange [incrMin[ARG3 - ARG4 + 1] : incrMax[ARG3 - ARG4 + 1]]
    set title "Increment test for Rastrigin function(N = ".ARG3."), method mggsa" font "Helvetica Bold, 20"

	set ylabel yName[ARG2 + 1]

    if (ARG2 == 3) {
        plot for [i = ARG6 : ARG7] datafile index (ARG3 - ARG4) * (ARG7 - ARG6 + 1) + i - ARG6 \
             using 1:($3 / $2) with lines lt i - ARG6 + 1 title "m = ".i
    } else {
        plot for [i = ARG6 : ARG7] datafile index (ARG3 - ARG4) * (ARG7 - ARG6 + 1) + i - ARG6 \
             using 1:(column(ARG2 + 2)) with lines lt i - ARG6 + 1 title "m = ".i
	}

    bind all "alt-End" "exit gnuplot"
    pause mouse close
} else {
    set terminal pngcairo size 1280, 720
	system "mkdir -p output_graph/".taskName

    do for [i = ARG4 : ARG5] {
       	set xrange [incrMin[i - ARG4 + 1] : incrMax[i - ARG4 + 1]]
       	set title "Increment test for Rastrigin function(N = ".i."), method mggsa" font "Helvetica Bold, 20"

		do for [j = 1 : 4] {
			set output "output_graph/".taskName."/".taskName."_".i."_".fileName[j].".png"

			set ylabel yName[j]

			if (j == 4) {
    		    plot for [k = ARG6 : ARG7] datafile index (i - ARG4) * (ARG7 - ARG6 + 1) + k - ARG6 \
    		         using 1:($3 / $2) with lines lt k - ARG6 + 1 title "m = ".k
    		} else {
    		    plot for [k = ARG6 : ARG7] datafile index (i - ARG4) * (ARG7 - ARG6 + 1) + k - ARG6 \
    		         using 1:(column(j + 1)) with lines lt k - ARG6 + 1 title "m = ".k
			}
		}
    }
}
