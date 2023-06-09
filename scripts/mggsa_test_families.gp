#!/usr/bin/gnuplot

taskName = "mggsa_test_families"

trialfile = "output_data/".taskName.".txt"

load "output_data/".taskName."_opt.txt"

title(familyName, n) = sprintf("Graph of function №%d of the %s family", n, familyName)
function(taskName, familyName) = sprintf("output_data/%s_%s_1.txt", taskName, familyName)
constraints(taskName, familyName) = sprintf("output_data/%s_%s_2.txt", taskName, familyName)

set grid

set contour
set view map
set cntrparam bspline levels auto 14
set cntrlabel onecolor
set cntrlabel start 5 interval 150
set cntrlabel font ",10"
set contour base

set xlabel "X"
set ylabel "Y"

set key opaque

set isosamples 120

if (ARG1 == 0) {
    ind = 3 * (ARG2 - 1)

    set title title(familyNames[int(ARG2)], functionNumber[int(ARG2)]) font "Helvetica Bold, 20"

    set xrange [A[2 * (ARG2 - 1) + 1] : B[2 * (ARG2 - 1) + 1]]
    set yrange [A[2 * (ARG2 - 1) + 2] : B[2 * (ARG2 - 1) + 2]]
    set zrange [minValue[int(ARG2)] :]

    splot function(taskName, familyNames[int(ARG2)]) matrix nonuniform title "f(x, y)" nosurface, \
          for [i = 1 : constrained[int(ARG2)]] constraints(taskName, familyNames[int(ARG2)]) matrix nonuniform \
          with lines lc rgb "orange" title "g(x, y)" nocontours, \
          trialfile index ind + 2 ls 4 lc rgb "green" lw 2 title "trial points" nocontours, \
          trialfile index ind     ls 7 lc rgb "red"   ps 2 title "X_{min}" nocontours, \
          trialfile index ind + 1 ls 7 lc rgb "blue"  ps 1 title "X" nocontours, \
          function(taskName, familyNames[int(ARG2)]) matrix nonuniform with labels notitle nosurface

    bind all "alt-End" "exit gnuplot"
    pause mouse close
} else {
    set terminal pngcairo size 950, 950
    system "mkdir -p output_graph/".taskName

    do for [i = 1 : 4] {
        set output "output_graph/".taskName."/".taskName."_".familyNames[i].".png"

        ind = 3 * (i - 1)

        set title title(familyNames[i], functionNumber[i]) font "Helvetica Bold, 20"

        set xrange [A[2 * (i - 1) + 1] : B[2 * (i - 1) + 1]]
        set yrange [A[2 * (i - 1) + 2] : B[2 * (i - 1) + 2]]
        set zrange [minValue[i] :]

        splot function(taskName, familyNames[i]) matrix nonuniform title "f(x, y)" nosurface, \
              for [j = 1 : constrained[i]] constraints(taskName, familyNames[i]) matrix nonuniform \
              with lines lc rgb "orange" title "g(x, y)" nocontours, \
              trialfile index ind + 2 ls 4 lc rgb "green" lw 2 title "trial points" nocontours, \
              trialfile index ind     ls 7 lc rgb "red"   lw 6 title "X_{min}" nocontours, \
              trialfile index ind + 1 ls 7 lc rgb "blue"  lw 1 title "X" nocontours, \
              function(taskName, familyNames[i]) matrix nonuniform with labels notitle nosurface
    }
}
