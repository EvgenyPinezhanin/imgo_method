#!/usr/bin/gnuplot

rootDir = "output_data/sample_fitting_problem"

load rootDir."/vars.txt"

functionPointsFile = rootDir."/function_points.txt"
testPointsFile = rootDir."/test_points.txt"

fontName = "Helvetica, 16"

set sample 500

set grid

set xlabel "t" font fontName
set ylabel "u(t)" font fontName offset -1

set tics font fontName

set key box outside right top
set key font fontName spacing 1.3

title(displayType, type, method) = displayType != 2 ? \
                                          type == 1 ? sprintf("Graph of the sample test function, %s method", method)     : \
                                                      sprintf("Graph of the solution of the equation, %s method", method) : \
                                                      sprintf("")

u(method, t) = c[method * numberCoefficients + 1] + c[method * numberCoefficients + 2] * t + \
               c[method * numberCoefficients + 3] * t * t + \
               c[method * numberCoefficients + 4] * sin(xOpt[method + 1] * t) + \
               c[method * numberCoefficients + 5] * cos(xOpt[method + 1] * t)

if (ARG1 == 0) {
    set title title(ARG1, ARG2, methodNames[ARG3 + 1]) font fontName

    trialsFile = rootDir."/".methodNames[ARG3 + 1]."_trials.txt"

    if (ARG2 == 0) {
        set lmargin 10

        set key width -1

        set xrange [ARG4 : ARG5]

        plot u(ARG3, x) title "u(t)", \
             testPointsFile ls 8 lc rgb "green" lw 6 title "test points"
    } else {
        set lmargin 11

        set key width -2

        set xrange [ARG6 : ARG7]

        plot functionPointsFile lc rgb "dark-violet" with lines title "f(x)", \
             trialsFile index 2 ls 4 lc rgb "green" lw 2 title "trial points", \
             trialsFile index 0 ls 7 lc rgb "red"   ps 2 title "X_{min}", \
             trialsFile index 1 ls 7 lc rgb "blue"  ps 1 title "X"
    }

    bind all "alt-End" "exit gnuplot"
    pause mouse close
} else {
    set terminal pngcairo size 1640, 950 font fontName

    system "mkdir -p output_graph/sample_fitting_problem"

    set lmargin 12
    set rmargin 18
    set tmargin 3
    set bmargin 3

    do for [i = 0 : 1] {
        do for [j = 0 : numberMethods - 1] {
            set title title(ARG1, i, methodNames[j + 1]) font fontName

            if (i == 0) {
                set output "output_graph/sample_fitting_problem/".methodNames[j + 1]."_solution.png"

                set lmargin 10

                set key width -1

                set xrange [ARG4 : ARG5]

                plot u(j, x) title "u(t)", \
                     testPointsFile ls 8 lc rgb "green" lw 6 title "test points"
            } else {
                trialsFile = rootDir."/".methodNames[j + 1]."_trials.txt"

                set output "output_graph/sample_fitting_problem/".methodNames[j + 1]."_function.png"

                set lmargin 11

                set key width -2

                set xrange [ARG6 : ARG7]

                plot functionPointsFile lc rgb "dark-violet" with lines title "f(x)", \
                     trialsFile index 2 ls 4 lc rgb "green" lw 2 title "trial points", \
                     trialsFile index 0 ls 7 lc rgb "red"   ps 2 title "X_{min}", \
                     trialsFile index 1 ls 7 lc rgb "blue"  ps 1 title "X"
            }
        }
    }
}
