#!/usr/bin/gnuplot

load "output_data/sample_test_problem_family/vars.txt"

# function_points_file = "output_data/sample_test_problem/function_points.txt"
test_points_file = "output_data/sample_test_problem_family/test_points.txt"
# trials_file = "output_data/sample_test_problem/".method_names[ARG3 + 1]."_trials.txt"

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

u(method, t) = c[method * number_coefficients + 1] * sin(x_opt[1] * t) + c[method * number_coefficients + 2] * cos(x_opt[1] * t) + \
               c[method * number_coefficients + 3] * sin(x_opt[2] * t) + c[method * number_coefficients + 4] * cos(x_opt[2] * t) + \
               c[method * number_coefficients + 5] * sin(x_opt[3] * t) + c[method * number_coefficients + 6] * cos(x_opt[3] * t) + \
               c[method * number_coefficients + 7] * sin(x_opt[4] * t) + c[method * number_coefficients + 8] * cos(x_opt[4] * t)

if (ARG1 == 0) {
    set title title(ARG1, ARG2, method_names[ARG3 + 1]) font font_name

    if (ARG2 == 0) {
        set lmargin 10

        set key width -1

        set xrange [ARG4 : ARG5]

        plot u(ARG3, x) title "u(t)" # , \
             # test_points_file ls 8 lc rgb "green" lw 6 notitle

        set obj 1 rect from A,-delta to B,delta front fs empty border rgb "red"
        set obj 2 rect from T[1]-delta,Q[1]-delta to T[1]+delta,Q[1]+delta front fs empty border rgb "red"
        set obj 3 rect from T[2]-delta,Q[2]-delta to T[2]+delta,Q[2]+delta front fs empty border rgb "red"
        set obj 4 rect from T[3]-delta,Q[3]-delta to T[3]+delta,Q[3]+delta front fs empty border rgb "red"

            bind all "alt-End" "exit gnuplot"
    pause mouse close
    } else {
taskName = "sample_test_problem_family"

title(f, s) = sprintf("Graph of the function slice by variables x%d and x%d", f, s)
function(taskName, f, s) = sprintf("output_data/%s/%d_%df.txt", taskName, f, s)
constraints(taskName, f, s) = sprintf("output_data/%s/%d_%dg.txt", taskName, f, s)

set grid

set contour
set view map
set cntrparam bspline levels auto 14
set cntrlabel onecolor
set cntrlabel start 5 interval 150
set cntrlabel font ",10"
set contour base

set key opaque
set key spacing 1.3

set isosamples 120

    set terminal pngcairo size 950, 950 font "Helvetica, 18"
    system "mkdir -p output_graph/".taskName

    set lmargin 1
    set rmargin 0
    set tmargin 0
    set bmargin 0

    do for [f = 1 : 3] {
        do for [s = f + 1 : 4] {
            set output "output_graph/".taskName."/".f."_".s.".png"

            set title title(f, s) font "Helvetica, 20"

            set xlabel "X".f
            set ylabel "Y".s

            set xrange [0.01 : 2.0]
            set yrange [0.01 : 2.0]
            set zrange [minValue :]

            splot function(taskName, f, s) matrix nonuniform title "f(x, y)" nosurface, \
                  constraints(taskName, f, s) matrix nonuniform with lines lc rgb "orange" title "g(x, y)" nocontours, \
                  function(taskName, f, s) matrix nonuniform with labels notitle nosurface
        }
    }
}
} else {
    // set terminal pngcairo size 1640, 950 font font_name
// 
    // system "mkdir -p output_graph/sample_test_problem"
// 
    // set lmargin 12
    // set rmargin 18
    // set tmargin 3
    // set bmargin 3
// 
    // do for [i = 0 : 1] {
    //     do for [j = 0 : number_methods - 1] {
    //         set title title(ARG1, i, method_names[j + 1]) font font_name
// 
    //         if (i == 0) {
    //             set output "output_graph/sample_test_problem/".method_names[j + 1]."_solution.png"
// 
    //             set lmargin 10
// 
    //             set key width -1
// 
    //             set xrange [ARG4 : ARG5]
// 
    //             plot u(j, x) title "u(t)", \
    //                  test_points_file ls 8 lc rgb "green" lw 6 title "test points"
    //         } else {
    //             trials_file = "output_data/sample_test_problem/".method_names[j + 1]."_trials.txt"
// 
    //             set output "output_graph/sample_test_problem/".method_names[j + 1]."_function.png"
// 
    //             set lmargin 11
// 
    //             set key width -2
// 
    //             set xrange [ARG6 : ARG7]
// 
    //             plot function_points_file lc rgb "dark-violet" with lines title "f(x)", \
    //                  trials_file index 2 ls 4 lc rgb "green" lw 2 title "trial points", \
    //                  trials_file index 0 ls 7 lc rgb "red"   ps 2 title "X_{min}", \
    //                  trials_file index 1 ls 7 lc rgb "blue"  ps 1 title "X"
    //         }
    //     }
    // }
}

        # set lmargin 11
# 
        # set key width -2
# 
        # set xrange [ARG6 : ARG7]
# 
        # plot function_points_file lc rgb "dark-violet" with lines title "f(x)", \
        #      trials_file index 2 ls 4 lc rgb "green" lw 2 title "trial points", \
        #      trials_file index 0 ls 7 lc rgb "red"   ps 2 title "X_{min}", \
        #      trials_file index 1 ls 7 lc rgb "blue"  ps 1 title "X"