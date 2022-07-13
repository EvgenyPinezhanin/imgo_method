#!/usr/bin/gnuplot

reset

trialfile="peano_sample_trial_points.txt"

f_1(x, y) = 1.0-x-y

f_2(x, y) = (x-1.0)**2 / 5.0 + (y - 1.0)**2 / 5.0

f_3(x, y) = x**2 / 5.0 + y**2 / 5.0
g1_3(x, y) = 1.0-x-y
g_3(x, y) = (g1_3(x, y) <= 0.0) ? g1_3(x, y) : 1/0

f_4(x, y) = x**2 / 5.0 + y**2 / 5.0
g1_4(x, y) = (x - 2.0)**2 + (y - 2.0)**2 - 2.0
g_4(x, y) = (g1_4(x, y) <= 0.0) ? g1_4(x, y) : 1/0

set grid
set contour
set view map
set cntrparam bspline levels auto 10
set cntrlabel onecolor
set cntrlabel start 5 interval 150
set cntrlabel font ",10"
set contour base
set isosamples 60

set xlabel "X"
set ylabel "Y"

if (ARG1 == 0) {
      set title "Chart of sample function 1" font "Helvetica Bold, 20"
      set xrange [-4.0:4.0]
      set yrange [-4.0:4.0]
      splot f_1(x, y) title "φ(x, y)" nosurface, \
            trialfile index 2 ls 5 lc rgb "green" title "trial points" nocontours, \
            trialfile index 1 ls 5 lc rgb "red" title "X*" nocontours, \
            trialfile index 0 ls 5 lc rgb "blue" title "X" nocontours, \
            f_1(x, y) with labels notitle nosurface
}
if (ARG1 == 1) {
      set title "Chart of sample function 2" font "Helvetica Bold, 20"
      set xrange [-4.0:4.0]
      set yrange [-4.0:4.0]
      splot f_2(x, y) title "φ(x, y)" nosurface, \
            trialfile index 5 ls 5 lc rgb "green" title "trial points" nocontours, \
            trialfile index 4 ls 5 lc rgb "red" title "X*" nocontours, \
            trialfile index 3 ls 5 lc rgb "blue" title "X" nocontours, \
            f_2(x, y) with labels notitle nosurface
}
if (ARG1 == 2) {
      set title "Chart of sample function 3" font "Helvetica Bold, 20"
      set xrange [-1.0:1.0]
      set yrange [-1.0:1.0]
      splot f_3(x, y) title "φ(x, y)" nosurface, \
            g_3(x, y) lc rgb "orange" notitle nocontours, \
            trialfile index 8 ls 5 lc rgb "green" title "trial points" nocontours, \
            trialfile index 7 ls 5 lc rgb "red" title "X*" nocontours, \
            trialfile index 6 ls 5 lc rgb "blue" title "X" nocontours, \
            f_3(x, y) with labels notitle nosurface
}
if (ARG1 == 3) {
      set title "Chart of sample function 4" font "Helvetica Bold, 20"
      set xrange [0.0:3.0]
      set yrange [0.0:3.0]
      splot f_4(x, y) title "φ(x, y)" nosurface, \
            g_4(x, y) lc rgb "orange" notitle nocontours, \
            trialfile index 11 ls 5 lc rgb "green" title "trial points" nocontours, \
            trialfile index 10 ls 5 lc rgb "red" title "X*" nocontours, \
            trialfile index 9 ls 5 lc rgb "blue" title "X" nocontours, \
            f_4(x, y) with labels notitle nosurface
}

bind all "alt-End" "exit gnuplot"
pause mouse close
