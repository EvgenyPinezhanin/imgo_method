#! /usr/bin/gnuplot

reset

datafile="trial_points.txt"

set hidden3d front
set xlabel "X"
set ylabel "Y"
set zlabel "Z"
set grid

if (ARG1 eq "func_1") {
      set xrange [-4.0:4.0]
      set yrange [-4.0:4.0]
      set zrange [-5.0:10.0]
      set title "FUNC 1" font "Helvetica Bold, 20"
      splot (x-1)**2/5+(y-1)**2/5, \
            datafile index 1 ls 5 lc rgb "red" title "opt point", \
            datafile index 2 ls 5 lc rgb "green" title "trial points", \
            datafile index 0 ls 5 lc rgb "blue" title "min point"
} else {
      if (ARG1 eq "func_2") {
            set xrange [-4.0:4.0]
            set yrange [-4.0:4.0]
            set zrange [-5.0:15.0]
            set title "FUNC 2" font "Helvetica Bold, 20"
            splot x**2/5+y**2/5, \
                  @ARG2 - x - y, \
                  datafile index 1 ls 5 lc rgb "red" title "opt point", \
                  datafile index 2 ls 5 lc rgb "green" title "trial points", \
                  datafile index 0 ls 5 lc rgb "blue" title "min point"
      } else {
            if (ARG1 eq "func_3") {
                  set xrange [0.0:4.0]
                  set yrange [-1.0:3.0]
                  set zrange [-50.0:15.0]
                  set title "FUNC 3" font "Helvetica Bold, 20"
                  splot -1.5*x**2*exp(1-x**2-20.25*(x-y)**2)-(0.5*(x-1)*(y-1))**4*exp(2-(0.5*(x-1))**4-(x-1)**4), \
                        0.01*((x-2.2)**2+(y-1.2)**2.0-2.25), \
                        100*(1-((x-2)**2)/1.44-(0.5*y)**2), \
                        10*(y-1.5-1.5*sin(6.283*(x-1.75))), \
                        datafile index 1 ls 5 lc rgb "red" title "opt point", \
                        datafile index 2 ls 5 lc rgb "green" title "trial points", \
                        datafile index 0 ls 5 lc rgb "blue" title "min point"
            } else {
                  set xrange [-4.0:4.0]
                  set yrange [-4.0:4.0]
                  set zrange [-10.0:10.0]
                  set title "FUNC 4" font "Helvetica Bold, 20"
                  splot 1-x-y, \
                        datafile index 1 ls 5 lc rgb "red" title "opt point", \
                        datafile index 2 ls 5 lc rgb "green" title "trial points", \
                        datafile index 0 ls 5 lc rgb "blue" title "min point"
            }
      }
}

# set label 1 "opt" at 1.0,1.0,1.0 point pt 7

bind all "alt-End" "exit gnuplot"
pause mouse close