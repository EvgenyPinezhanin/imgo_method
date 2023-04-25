#! /usr/bin/gnuplot

datafile = "output_data/mggsa_test_density.txt"

load "output_data/mggsa_test_density_opt.txt"

set linetype 1 lc rgb "salmon" pt 5 
set linetype 2 lc rgb "cyan" pt 5 
set linetype 3 lc rgb "blue" pt 5 
set linetype 4 lc rgb "orange" pt 5
set linetype 5 lc rgb "green" pt 5

set linetype 6 lc rgb "salmon" pt 2 ps 2
set linetype 7 lc rgb "cyan" pt 2 ps 2
set linetype 8 lc rgb "blue" pt 2 ps 2
set linetype 9 lc rgb "orange" pt 2 ps 2
set linetype 10 lc rgb "green" pt 2 ps 2

set linetype 11 lc rgb "dark-violet" lw 1
set linetype 12 lc rgb "red" pt 5

set linetype cycle 12

f(x, y) = x ** 2 + y ** 2 - cos(18.0 * x) - cos(18.0 * y)

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

set key box
set key opaque

set xrange [-0.5:1.0]
set yrange [-0.5:1.0]

set title "Density test" font "Helvetica Bold, 20"

splot f(x, y) title "φ(x, y)" lt 11 nosurface, \
      for [i = 1:@ARG1] datafile index 2 * (@ARG1 - i + 1) lt i title "trial points(m = ".den[@ARG1 - i + 1].")" nocontours, \
      for [i = 1:@ARG1] datafile index 2 * (@ARG1 - i + 1) - 1 lt i + 5 title "X(m = ".den[@ARG1 - i + 1].")" nocontours, \
      datafile index 0 lt 12 title "X*" nocontours

bind all "alt-End" "exit gnuplot"
pause mouse close
