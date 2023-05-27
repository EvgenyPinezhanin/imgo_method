#!/usr/bin/gnuplot

trialfile = "output_data/mggsa_test_families.txt"

load "output_data/mggsa_test_families_opt.txt"

titleName(s, n) = sprintf("Graph of function №%d of the %s family", n, s)
function(type, n) = sprintf("\"output_data/mggsa_test_families_%d_1.txt\" matrix nonuniform %s", n, type)
constraints(n) = (n < 3) ? sprintf("") : \
                 sprintf("\"output_data/mggsa_test_families_%d_2.txt\" matrix nonuniform with lines lc rgb \"orange\" title \"g(x, y)\" nocontours,", n)

set grid

set contour
set view map
set cntrparam bspline levels auto 14
set cntrlabel onecolor
set cntrlabel start 5 interval 150
set cntrlabel font ",10"
set contour base
set isosamples 120

set xlabel "X"
set ylabel "Y"

set datafile missing "?"

set terminal wxt size 950, 950

set title titleName(familyNames[ARG1 + 1], functionNumber[ARG1 + 1]) font "Helvetica Bold, 20"

set xrange [A[2 * ARG1 + 1]:B[2 * ARG1 + 1]]
set yrange [A[2 * ARG1 + 2]:B[2 * ARG1 + 2]]

ind = 3 * ARG1

function = function("nosurface with lines title \"φ(x, y)\"", ARG1 + 1)
constraints = constraints(ARG1 + 1)
functionLabels = function("with labels notitle nosurface", ARG1 + 1)

splot @function, @constraints \
      trialfile index ind + 2 ls 5 lc rgb "green" title "trial points" nocontours, \
      trialfile index ind + 1 ls 5 lc rgb "blue" title "X" nocontours, \
      trialfile index ind ls 5 lc rgb "red" title "X*" nocontours # , \
      # @functionLabels

bind all "alt-End" "exit gnuplot"
pause mouse close
