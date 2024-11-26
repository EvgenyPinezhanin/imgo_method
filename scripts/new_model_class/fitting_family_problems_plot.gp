#!/usr/bin/gnuplot

sampleName = "fitting_family_problems_plot"

fontName = "Helvetica, 26"

set sample 500

set grid

set tics font fontName

set key box outside right top
set key font fontName spacing 1.3

title(displayType, graphType, problemNumber) = displayType != 2 ? \
                                                   graphType == 1 ? sprintf("Graph of the fitting problem №%d solution", problemNumber) : \
                                                   sprintf("") : \
                                               sprintf("")

if (ARG2 == 1) {
    load "output_data/new_model_class/".sampleName."/vars_solution.txt"
}

u(t, i) = coeffs[8 * i + 1] * sin(x_opt[4 * i + 1] * t) + coeffs[8 * i + 2] * cos(x_opt[4 * i + 1] * t) + \
          coeffs[8 * i + 3] * sin(x_opt[4 * i + 2] * t) + coeffs[8 * i + 4] * cos(x_opt[4 * i + 2] * t) + \
          coeffs[8 * i + 5] * sin(x_opt[4 * i + 3] * t) + coeffs[8 * i + 6] * cos(x_opt[4 * i + 3] * t) + \
          coeffs[8 * i + 7] * sin(x_opt[4 * i + 4] * t) + coeffs[8 * i + 8] * cos(x_opt[4 * i + 4] * t)

if (ARG2 == 1) {
    set xlabel "t" font fontName
    set ylabel "u(t)" font fontName

    if (ARG1 == 0) {
        set title title(ARG1, ARG2, int(ARG3)) font fontName

        set lmargin 13

        set ylabel offset -3

        set key inside left bottom
    
        set xrange [left_bound - 1.0 : right_bound + 1.0]
    
        plot u(x, ARG3) title "u(t)"
    
        set obj 1 rect from A, -delta to B, delta front fs empty border rgb "red"
        set obj 2 rect from T[3 * ARG3 + 1] - delta, Q[3 * ARG3 + 1] - delta to \
                            T[3 * ARG3 + 1] + delta, Q[3 * ARG3 + 1] + delta front fs empty border rgb "red"
        set obj 3 rect from T[3 * ARG3 + 2] - delta, Q[3 * ARG3 + 2] - delta to \
                            T[3 * ARG3 + 2] + delta, Q[3 * ARG3 + 2] + delta front fs empty border rgb "red"
        set obj 4 rect from T[3 * ARG3 + 3] - delta, Q[3 * ARG3 + 3] - delta to \
                            T[3 * ARG3 + 3] + delta, Q[3 * ARG3 + 3] + delta front fs empty border rgb "red"
    
        bind all "alt-End" "exit gnuplot"
        pause mouse close
    } else {
        set terminal pngcairo size 1800, 850 font fontName
        system "mkdir -p output_graph/new_model_class/".sampleName."/solution"

        set lmargin 8

        set ylabel offset -1

        set key inside left bottom
    
        set xrange [left_bound - 1.0 : right_bound + 1.0]

        do for [i = 0 : familySize - 1] {
            set output "output_graph/new_model_class/".sampleName."/solution/".i.".png"
            
            set title title(ARG1, ARG2, i) font fontName

            set obj 1 rect from A, -delta to B, delta front fs empty border rgb "red"
            set obj 2 rect from T[3 * i + 1] - delta, Q[3 * i + 1] - delta to \
                                T[3 * i + 1] + delta, Q[3 * i + 1] + delta front fs empty border rgb "red"
            set obj 3 rect from T[3 * i + 2] - delta, Q[3 * i + 2] - delta to \
                                T[3 * i + 2] + delta, Q[3 * i + 2] + delta front fs empty border rgb "red"
            set obj 4 rect from T[3 * i + 3] - delta, Q[3 * i + 3] - delta to \
                                T[3 * i + 3] + delta, Q[3 * i + 3] + delta front fs empty border rgb "red"

            plot u(x, i) title "u(t)"
        }
    }
} else {
    # if (ARG2 == 1) {
# 
    # }
    # if (ARG2 == 0) {
    #     set terminal pngcairo size 1800, 1000 font fontName
    #     
    #     system "mkdir -p output_graph/".sampleName
    #     
    #     set lmargin 5
    #     set rmargin 5
    #     set tmargin 2
    #     set bmargin 2
    # 
    #     set key width -11
    #     set key inside left bottom
    #     
    #     set title title(ARG1, ARG2) font fontName
    #     
    #     set xrange [A - 1 : T[3] + 1]
    #     set yrange [-14.0 : 5.0]
    #     
    #     set output "output_graph/".sampleName."/u(t).png"
    #     
    #     set obj 1 rect from A,-delta to B,delta back lw 4 dt 1 fs empty border rgb "red"
    #     set obj 2 rect from T[1]-delta,Q[1]-delta to T[1]+delta,Q[1]+delta back lw 4 dt 1 fs empty border rgb "red"
    #     set obj 3 rect from T[2]-delta,Q[2]-delta to T[2]+delta,Q[2]+delta back lw 4 dt 1 fs empty border rgb "red"
    #     set obj 4 rect from T[3]-delta,Q[3]-delta to T[3]+delta,Q[3]+delta back lw 4 dt 1 fs empty border rgb "red"
    #     
    #     plot for [i = 0 : 3] u(x, i) lw 3 title "u(t) (".num[i + 1]." испытаний)" 
# 
    #     # set terminal pngcairo size 1800, 850 font fontName
    #     # 
    #     # system "mkdir -p output_graph/".sampleName
    #     # 
    #     # set lmargin 6
    #     # set rmargin 6
    #     # set tmargin 2
    #     # set bmargin 2
    #     # 
    #     # set key width -1
    #     # set key inside left top
    #     # 
    #     # set title title(ARG1, ARG2) font fontName
    #     # 
    #     # set xrange [A - 1 : T[3] + 1]
    #     # set yrange [-9.0 : 6.0]
    #     # 
    #     # set output "output_graph/".sampleName."/u(t).png"
    #     # 
    #     # set obj 1 rect from A,-delta to B,delta back lw 4 dt 1 fs empty border rgb "red"
    #     # set obj 2 rect from T[1]-delta,Q[1]-delta to T[1]+delta,Q[1]+delta back lw 4 dt 1 fs empty border rgb "red"
    #     # set obj 3 rect from T[2]-delta,Q[2]-delta to T[2]+delta,Q[2]+delta back lw 4 dt 1 fs empty border rgb "red"
    #     # set obj 4 rect from T[3]-delta,Q[3]-delta to T[3]+delta,Q[3]+delta back lw 4 dt 1 fs empty border rgb "red"
    #     # 
    #     # plot u(x, 0) lw 3 title "u(t)"
    # } else {
    #     # title(f, s, number) = sprintf("Graph of the function slice by variables x%d and x%d, №%d", f, s, number)
    #     function(sampleName, number, f, s) = sprintf("output_data/%s/slices/%d_%d_%df.txt", sampleName, number, f, s)
    #     constraints(sampleName, number, f, s) = sprintf("output_data/%s/slices/%d_%d_%dg.txt", sampleName, number, f, s)
# 
    #     set contour
    #     set view map
    #     set cntrparam bspline levels auto 14
    #     set cntrlabel onecolor
    #     set cntrlabel start 5 interval 150
    #     set cntrlabel font ",10"
    #     set contour base
# 
    #     set key opaque
    #     set key inside top left
# 
    #     set isosamples 120
# 
    #     set terminal pngcairo size 950, 950 font fontName
    #     system "mkdir -p output_graph/".sampleName
# 
    #     set lmargin 1
    #     set rmargin 0
    #     set tmargin 0
    #     set bmargin 0
# 
    #     do for [number = 1 : 1] {
    #         system "mkdir -p output_graph/".sampleName."/".number
# 
    #         do for [f = 1 : 3] {
    #             do for [s = f + 1 : 4] {
    #                 set output "output_graph/".sampleName."/".number."/".f."_".s.".png"
# 
    #                 # set title title(f, s, number) font "Helvetica, 20"
# 
    #                 set xlabel "X".f
    #                 set ylabel "Y".s
# 
    #                 set xrange [0.01 : 2.0]
    #                 set yrange [0.01 : 2.0]
    #                 set zrange [optValues[number] :]
# 
    #                 splot function(sampleName, number, f, s) matrix nonuniform title "f(X)" nosurface, \
    #                       constraints(sampleName, number, f, s) matrix nonuniform with lines lc rgb "orange" title "g(X)" nocontours, \
    #                       function(sampleName, number, f, s) matrix nonuniform with labels notitle nosurface
    #             }
    #         }
    #     }
    # }
}
