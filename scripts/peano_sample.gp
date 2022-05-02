#!/usr/bin/gnuplot

reset

trialfile=ARG1."_trial_points.txt"
funcfile=ARG1."_function.txt"
imagedir=ARG1."_res"

# plot_f(i) = sprintf("splot f_%d(x, y)", i)
# plot_g(i, j) = sprintf("replot g%d_%d(x, y)", j, i)  

plot_f(i, j) = (j > 0) ? gprintf("g%d_%d(x, y)", j, i) : gprintf("f_%d(x, y)", i) 

load funcfile

set terminal png size 1200,720

unset surface
set contour
set view map
set cntrparam levels auto 10

set xlabel "X"
set ylabel "Y"

do for [i=1:func_size] {
    if (Used[i] == 1) {
        name = Name[i]
        print name
        set xrange [A[2*i-1]:B[2*i-1]]
        set yrange [A[2*i]:B[2*i]]
        # eval plot_f(i)
        # do for [j=1:M[i]] {
        #     eval plot_g(i, j)
        # }
        splot for [j=0:0] plot_f(i, j)
        set output
        unset output
    }
}

name = Name[2]
print name
set output imagedir."/".name.".png"
set output
unset output

# set hidden3d front
# set xlabel "X"
# set ylabel "Y"
# set zlabel "Z"
# set grid
# 
# if (ARG1 eq "func_1") {
#       set xrange [-4.0:4.0]
#       set yrange [-4.0:4.0]
#       set zrange [-5.0:10.0]
#       set title "FUNC 1" font "Helvetica Bold, 20"
#       splot (x-1)**2/5+(y-1)**2/5, \
#             datafile index 1 ls 5 lc rgb "red" title "opt point", \
#             datafile index 2 ls 5 lc rgb "green" title "trial points", \
#             datafile index 0 ls 5 lc rgb "blue" title "min point"
# } else {
#       if (ARG1 eq "func_2") {
#             set xrange [-4.0:4.0]
#             set yrange [-4.0:4.0]
#             set zrange [-5.0:15.0]
#             set title "FUNC 2" font "Helvetica Bold, 20"
#             splot x**2/5+y**2/5, \
#                   @ARG2 - x - y, \
#                   datafile index 1 ls 5 lc rgb "red" title "opt point", \
#                   datafile index 2 ls 5 lc rgb "green" title "trial points", \
#                   datafile index 0 ls 5 lc rgb "blue" title "min point"
#       } else {
#             set xrange [-4.0:4.0]
#             set yrange [-4.0:4.0]
#             set zrange [-10.0:10.0]
#             set title "FUNC 4" font "Helvetica Bold, 20"
#             splot 1-x-y, \
#                   datafile index 1 ls 5 lc rgb "red" title "opt point", \
#                   datafile index 2 ls 5 lc rgb "green" title "trial points", \
#                   datafile index 0 ls 5 lc rgb "blue" title "min point"
#       }
# }

#                  set xrange [0.0:4.0]
#                  set yrange [-1.0:3.0]
#                  set zrange [-50.0:15.0]
#                  set title "FUNC 3" font "Helvetica Bold, 20"
#                  splot -1.5*x**2*exp(1-x**2-20.25*(x-y)**2)-(0.5*(x-1)*(y-1))**4*exp(2-(0.5*(x-1))**4-(x-1)**4), \
#                        0.01*((x-2.2)**2+(y-1.2)**2.0-2.25), \
#                        100*(1-((x-2)**2)/1.44-(0.5*y)**2), \
#                        10*(y-1.5-1.5*sin(6.283*(x-1.75))), \
#                        datafile index 1 ls 5 lc rgb "red" title "opt point", \
#                        datafile index 2 ls 5 lc rgb "green" title "trial points", \
#                        datafile index 0 ls 5 lc rgb "blue" title "min point"

# do for [FILETYPE in "png eps"] {
#  set output "plot.".FILETYPE
#  set term FILETYPE
#  plot x**2
#  set output
# }