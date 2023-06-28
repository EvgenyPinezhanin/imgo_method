#! /usr/bin/gnuplot

taskName = "mggsa_test_incr"
datafile = "output_data/".taskName.".txt"
load "output_data/".taskName."_opt.txt"

# set linetype 1 lc rgb "#0000FF" lw 5 pt 1 dt 1
# set linetype 2 lc rgb "#0020FF" lw 4 pt 1 dt 2
# set linetype 3 lc rgb "#0040FF" lw 3 pt 1 dt 3
# set linetype 4 lc rgb "#0060FF" lw 2 pt 1 dt 4
# set linetype 5 lc rgb "#0080FF" lw 1 pt 1 dt 5

set linetype 1 lc rgb "red" lw 5
set linetype 2 lc rgb "green" lw 4
set linetype 3 lc rgb "blue" lw 3

set linetype cycle 3

set xlabel "incr"

set key spacing 1.3

array yName[4]
yName[1] = "trials"
yName[2] = "point trials"
yName[3] = "accuracy"
yName[4] = "Количество точек / количество испытаний"

array fileName[4]
fileName[1] = "trials"
fileName[2] = "point_trials"
fileName[3] = "accuracy"
fileName[4] = "point_trials_div_trials"

title(n) = sprintf("Increment test for Rastrigin function(N = %d), method mggsa", n)
titlePng(n) = (ARG1 == 1) ? title(n) : sprintf("")

if (ARG1 == 0) {
	set key font ", 15"

    set xrange [incrMin[ARG3 - ARG4 + 1] : incrMax[ARG3 - ARG4 + 1]]

	set ylabel yName[ARG2 + 1]

    set title title(ARG3) font "Helvetica Bold, 20"
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
    set terminal pngcairo size 1440, 700 font "Helvetica, 18"
	system "mkdir -p output_graph/".taskName

	set lmargin 8
    set rmargin 6
    set tmargin 3
    set bmargin 3

    do for [i = ARG4 : ARG5] {
       	set xrange [incrMin[i - ARG4 + 1] : incrMax[i - ARG4 + 1]]

       	set title titlePng(i) font "Helvetica, 18"
		do for [j = 1 : 4] {
			set output "output_graph/".taskName."/".i."_".fileName[j].".png"

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
