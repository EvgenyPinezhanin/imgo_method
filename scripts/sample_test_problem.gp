#!/usr/bin/gnuplot

load "output_data/sample_test_problem/vars.txt"

function_points_file = "output_data/sample_test_problem/function_points.txt"
test_points_file = "output_data/sample_test_problem/test_points.txt"
trials_file = "output_data/sample_test_problem/".method_names[ARG3 + 1]."_trials.txt"

font_name = "Helvetica, 16"

set sample 500

set grid

set xlabel "t" font font_name
set ylabel "u(t)" font font_name offset -1

set tics font font_name

set key box outside right top
set key font font_name spacing 1.3

title(display_type, type, method) = display_type != 2 ? \
                                            type == 1 ? sprintf("Graph of the sample test function, %s method", method) : \
                                                        sprintf("Graph of the solution of the equation, %s method", method) : \
                                                        sprintf("")

u(method, t) = c[method * number_coefficients + 1] + c[method * number_coefficients + 2] * t + \
               c[method * number_coefficients + 3] * t * t + \
               c[method * number_coefficients + 4] * sin(x_opt[method + 1] * t) + \
               c[method * number_coefficients + 5] * cos(x_opt[method + 1] * t)

if (ARG1 == 0) {
    set title title(ARG1, ARG2, method_names[ARG3 + 1]) font font_name

    if (ARG2 == 0) {
        set lmargin 10

        set key width -1

        set xrange [ARG4 : ARG5]

        plot u(ARG3, x) title "u(t)", \
             test_points_file ls 8 lc rgb "green" lw 6 title "test points"
    } else {
        set lmargin 11

        set key width -2

        set xrange [ARG6 : ARG7]

        plot function_points_file lc rgb "dark-violet" with lines title "f(x)", \
             trials_file index 2 ls 4 lc rgb "green" lw 2 title "trial points", \
             trials_file index 0 ls 7 lc rgb "red"   ps 2 title "X_{min}", \
             trials_file index 1 ls 7 lc rgb "blue"  ps 1 title "X"
    }

    bind all "alt-End" "exit gnuplot"
    pause mouse close
} else {
    set terminal pngcairo size 1640, 950 font font_name

    system "mkdir -p output_graph/sample_test_problem"

    set lmargin 12
    set rmargin 18
    set tmargin 3
    set bmargin 3

    do for [i = 0 : 1] {
        do for [j = 0 : number_methods - 1] {
            set title title(ARG1, i, method_names[j + 1]) font font_name

            if (i == 0) {
                set output "output_graph/sample_test_problem/".method_names[j + 1]."_solution.png"

                set lmargin 10

                set key width -1

                set xrange [ARG4 : ARG5]

                plot u(j, x) title "u(t)", \
                     test_points_file ls 8 lc rgb "green" lw 6 title "test points"
            } else {
                trials_file = "output_data/sample_test_problem/".method_names[j + 1]."_trials.txt"

                set output "output_graph/sample_test_problem/".method_names[j + 1]."_function.png"

                set lmargin 11

                set key width -2

                set xrange [ARG6 : ARG7]

                plot function_points_file lc rgb "dark-violet" with lines title "f(x)", \
                     trials_file index 2 ls 4 lc rgb "green" lw 2 title "trial points", \
                     trials_file index 0 ls 7 lc rgb "red"   ps 2 title "X_{min}", \
                     trials_file index 1 ls 7 lc rgb "blue"  ps 1 title "X"
            }
        }
    }
}
