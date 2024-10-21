#!/usr/bin/gnuplot

sampleName = "fitting_family_problems"

load "output_data/".sampleName."/vars.txt"

array optValues[100]
optValues[1] = -5.11877
optValues[2] = -5.77019
optValues[3] = -14.0882
optValues[4] = -9.18095
optValues[5] = -10.3195
optValues[6] = -11.8162
optValues[7] = -17.2099
optValues[8] = -7.6189
optValues[9] = -6.43766
optValues[10] = -3.3912
optValues[11] = -17.6577
optValues[12] = -7.35867
optValues[13] = -17.1535
optValues[14] = -8.50832
optValues[15] = -8.42136
optValues[16] = -12.4695
optValues[17] = -6.87875
optValues[18] = -7.08501
optValues[19] = -16.2055
optValues[20] = -15.6066
optValues[21] = -7.30814
optValues[22] = -15.3084
optValues[23] = -13.0607
optValues[24] = -3.12091
optValues[25] = -14.5826
optValues[26] = -5.95241
optValues[27] = -6.00309
optValues[28] = -9.77678
optValues[29] = -20.4846
optValues[30] = -7.81804
optValues[31] = -13.5663
optValues[32] = -16.0088
optValues[33] = -5.49887
optValues[34] = -5.06972
optValues[35] = -16.3472
optValues[36] = -9.28037
optValues[37] = -11.4542
optValues[38] = -9.25541
optValues[39] = -5.42001
optValues[40] = -4.18594
optValues[41] = -9.0378
optValues[42] = -7.3331
optValues[43] = -3.50313
optValues[44] = -8.93695
optValues[45] = -5.54441
optValues[46] = -10.5604
optValues[47] = -13.2762
optValues[48] = -7.41819
optValues[49] = -13.8957
optValues[50] = -3.77639
optValues[51] = -12.3072
optValues[52] = -8.63458
optValues[53] = -6.63607
optValues[54] = -4.7894
optValues[55] = -14.3287
optValues[56] = -8.67452
optValues[57] = -5.16556
optValues[58] = -12.6835
optValues[59] = -9.58347
optValues[60] = -6.61801
optValues[61] = -5.8946
optValues[62] = -4.46882
optValues[63] = -10.4205
optValues[64] = -15.6701
optValues[65] = -12.4687
optValues[66] = -4.65406
optValues[67] = -5.11613
optValues[68] = -6.59087
optValues[69] = -5.65083
optValues[70] = -13.1604
optValues[71] = -12.811
optValues[72] = -8.60641
optValues[73] = -11.1435
optValues[74] = -5.13746
optValues[75] = -4.60812
optValues[76] = -7.57139
optValues[77] = -8.33871
optValues[78] = -25.7082
optValues[79] = -5.87832
optValues[80] = -9.98561
optValues[81] = -24.6592
optValues[82] = -25.5043
optValues[83] = -11.9244
optValues[84] = -10.6042
optValues[85] = -5.93141
optValues[86] = -13.9058
optValues[87] = -7.58089
optValues[88] = -5.77671
optValues[89] = -5.50259
optValues[90] = -26.3171
optValues[91] = -11.5467
optValues[92] = -6.26045
optValues[93] = -7.03579
optValues[94] = -10.1698
optValues[95] = -13.3914
optValues[96] = -19.917
optValues[97] = -6.33359
optValues[98] = -10.9576
optValues[99] = -6.19937
optValues[100] = -8.9841

fontName = "Helvetica, 26"

set sample 500

set grid

set tics font fontName

set key box outside right top
set key font fontName spacing 1.3

title(displayType, graphType) = displayType != 2 ? \
                                  graphType == 1 ? sprintf("Graph of the sample test function, MGGSA method") : \
                                                   sprintf("Graph of the solution of the equation, MGGSA method") : \
                                                   sprintf("")

u(t, i) = coeffs[8 * i + 1] * sin(x_opt[4 * i + 1] * t) + coeffs[8 * i + 2] * cos(x_opt[4 * i + 1] * t) + \
          coeffs[8 * i + 3] * sin(x_opt[4 * i + 2] * t) + coeffs[8 * i + 4] * cos(x_opt[4 * i + 2] * t) + \
          coeffs[8 * i + 5] * sin(x_opt[4 * i + 3] * t) + coeffs[8 * i + 6] * cos(x_opt[4 * i + 3] * t) + \
          coeffs[8 * i + 7] * sin(x_opt[4 * i + 4] * t) + coeffs[8 * i + 8] * cos(x_opt[4 * i + 4] * t)

array num[4]
# num[1] = 500
# num[2] = 1000
num[1] = 50
num[2] = 500
num[3] = 5000
num[4] = 10000

if (ARG1 == 0) {
    # set xlabel "t" font fontName
    # set ylabel "u(t)" font fontName offset -1
    # set title title(ARG1, ARG2, method_names[ARG3 + 1]) font fontName

    # function_points_file = "output_data/sample_test_problem/function_points.txt"
    # test_points_file = "output_data/sample_test_problem_family/test_points.txt"
    # trials_file = "output_data/sample_test_problem/".method_names[ARG3 + 1]."_trials.txt"

    # if (ARG2 == 0) {
    #     set lmargin 10
    # 
    #     set key width -1
    # 
    #     set xrange [ARG4 : ARG5]
    # 
    #     plot u(ARG3, x) title "u(t)"
    # 
    #     set obj 1 rect from A,-delta to B,delta front fs empty border rgb "red"
    #     set obj 2 rect from T[1]-delta,Q[1]-delta to T[1]+delta,Q[1]+delta front fs empty border rgb "red"
    #     set obj 3 rect from T[2]-delta,Q[2]-delta to T[2]+delta,Q[2]+delta front fs empty border rgb "red"
    #     set obj 4 rect from T[3]-delta,Q[3]-delta to T[3]+delta,Q[3]+delta front fs empty border rgb "red"
    # 
    #     bind all "alt-End" "exit gnuplot"
    #     pause mouse close
    # } else {}
} else {
    if (ARG2 == 0) {
        set terminal pngcairo size 1800, 1000 font fontName
        
        system "mkdir -p output_graph/".sampleName
        
        set lmargin 5
        set rmargin 5
        set tmargin 2
        set bmargin 2

        set key width -11
        set key inside left bottom
        
        set title title(ARG1, ARG2) font fontName
        
        set xrange [A - 1 : T[3] + 1]
        set yrange [-14.0 : 5.0]
        
        set output "output_graph/".sampleName."/u(t).png"
        
        set obj 1 rect from A,-delta to B,delta back lw 4 dt 1 fs empty border rgb "red"
        set obj 2 rect from T[1]-delta,Q[1]-delta to T[1]+delta,Q[1]+delta back lw 4 dt 1 fs empty border rgb "red"
        set obj 3 rect from T[2]-delta,Q[2]-delta to T[2]+delta,Q[2]+delta back lw 4 dt 1 fs empty border rgb "red"
        set obj 4 rect from T[3]-delta,Q[3]-delta to T[3]+delta,Q[3]+delta back lw 4 dt 1 fs empty border rgb "red"
        
        plot for [i = 0 : 3] u(x, i) lw 3 title "u(t) (".num[i + 1]." испытаний)" 

        # set terminal pngcairo size 1800, 850 font fontName
        # 
        # system "mkdir -p output_graph/".sampleName
        # 
        # set lmargin 6
        # set rmargin 6
        # set tmargin 2
        # set bmargin 2
        # 
        # set key width -1
        # set key inside left top
        # 
        # set title title(ARG1, ARG2) font fontName
        # 
        # set xrange [A - 1 : T[3] + 1]
        # set yrange [-9.0 : 6.0]
        # 
        # set output "output_graph/".sampleName."/u(t).png"
        # 
        # set obj 1 rect from A,-delta to B,delta back lw 4 dt 1 fs empty border rgb "red"
        # set obj 2 rect from T[1]-delta,Q[1]-delta to T[1]+delta,Q[1]+delta back lw 4 dt 1 fs empty border rgb "red"
        # set obj 3 rect from T[2]-delta,Q[2]-delta to T[2]+delta,Q[2]+delta back lw 4 dt 1 fs empty border rgb "red"
        # set obj 4 rect from T[3]-delta,Q[3]-delta to T[3]+delta,Q[3]+delta back lw 4 dt 1 fs empty border rgb "red"
        # 
        # plot u(x, 0) lw 3 title "u(t)"
    } else {
        # title(f, s, number) = sprintf("Graph of the function slice by variables x%d and x%d, №%d", f, s, number)
        function(sampleName, number, f, s) = sprintf("output_data/%s/slices/%d_%d_%df.txt", sampleName, number, f, s)
        constraints(sampleName, number, f, s) = sprintf("output_data/%s/slices/%d_%d_%dg.txt", sampleName, number, f, s)

        set contour
        set view map
        set cntrparam bspline levels auto 14
        set cntrlabel onecolor
        set cntrlabel start 5 interval 150
        set cntrlabel font ",10"
        set contour base

        set key opaque
        set key inside top left

        set isosamples 120

        set terminal pngcairo size 950, 950 font fontName
        system "mkdir -p output_graph/".sampleName

        set lmargin 1
        set rmargin 0
        set tmargin 0
        set bmargin 0

        do for [number = 1 : 1] {
            system "mkdir -p output_graph/".sampleName."/".number

            do for [f = 1 : 3] {
                do for [s = f + 1 : 4] {
                    set output "output_graph/".sampleName."/".number."/".f."_".s.".png"

                    # set title title(f, s, number) font "Helvetica, 20"

                    set xlabel "X".f
                    set ylabel "Y".s

                    set xrange [0.01 : 2.0]
                    set yrange [0.01 : 2.0]
                    set zrange [optValues[number] :]

                    splot function(sampleName, number, f, s) matrix nonuniform title "f(X)" nosurface, \
                          constraints(sampleName, number, f, s) matrix nonuniform with lines lc rgb "orange" title "g(X)" nocontours, \
                          function(sampleName, number, f, s) matrix nonuniform with labels notitle nosurface
                }
            }
        }
    }
}
