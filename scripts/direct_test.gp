#!/usr/bin/gnuplot

trialfile = "output_data/direct_test.txt"

load "output_data/direct_test_opt.txt"

f1(x, y) = -1.5 * x ** 2 * exp(1 - x ** 2 - 20.25 * (x - y) ** 2) - (0.5 * (x - 1) * (y - 1)) ** 4 * \
            exp(2.0 - (0.5 * (x - 1)) ** 4 - (y - 1.0) ** 4)
g1_1(x, y) = 0.01 * ((x - 2.2) ** 2 + (y - 1.2) ** 2 - 2.25)
g1_2(x, y) = 100.0 * (1.0 - ((x - 2.0) ** 2) / 1.44 - (0.5 * y) ** 2)
g1_3(x, y) = 10.0 * (y - 1.5 - 1.5 * sin(6.283 * (x - 1.75)))
g1(x, y) = (g1_1(x, y) <= 0.0) ? (g1_2(x, y) <= 0 ? (g1_3(x, y) <= 0 ? g1_3(x, y) : 1 / 0) : 1 / 0) : 1 / 0

array C[20]
C[1]=75.1963666677
C[2]=-3.8112755343
C[3]=0.1269366345
C[4]=-0.0020567665
C[5]=0.000010345
C[6]=-6.8306567631
C[7]=0.0302344793
C[8]=-0.0012813448
C[9]=0.0000352559
C[10]=-0.0000002266
C[11]=0.2564581253
C[12]=-0.0034604030
C[13]=0.0000135139
C[14]=-28.1064434908
C[15]=-0.0000052375
C[16]=-0.0000000063
C[17]=0.0000000007
C[18]=0.0003405462
C[19]=-0.0000016638
C[20]=-2.8673112392

f2(x, y) = -(C[1] + C[2] * x + C[3] * x ** 2 + C[4] * x ** 3 + C[5] * x ** 4 + C[6] * y + C[7] * x * y + C[8] * x ** 2 * y \
           + C[9] * x ** 3 * y + C[10] * x ** 4 * y + C[11] * y ** 2 + C[12] * y ** 3 + C[13] * y ** 4 + C[14] / (y + 1) + \
           C[15] * x ** 2 * y ** 2 + C[16] * x ** 3 * y ** 2 + C[17] * x ** 3 * y ** 3 + C[18] * x * y ** 2 + \
           C[19] * x * y ** 3 + C[20] * exp(0.0005 * x * y))
g2_1(x, y) = 450.0 - x * y
g2_2(x, y) = (0.1 * x - 1.0) ** 2 - y
g2_3(x, y) = 8.0 * (x - 40.0) - (y - 30.0) * (y - 55.0)
g2_4(x, y) = y + (x - 35.0) * (x - 30.0) / 125.0 - 80.0
g2(x, y) = (g2_1(x, y) <= 0.0) ? (g2_2(x, y) <= 0 ? (g2_3(x, y) <= 0 ? (g2_4(x, y) <= 0 ? (g2_4(x, y)) : 1 / 0) : 1 / 0) : 1 / 0) : 1 / 0

title_name(n) = sprintf("Graph of the test function №%d with DIRECT", n)
function_name(s, n) = sprintf("%s%d(x, y)", s, n)

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

set terminal wxt size 950, 950

ind = 3 * ARG1
function_f = function_name("f", ARG1 + 1)
function_g = function_name("g", ARG1 + 1)

set title title_name(ARG1 + 1) font "Helvetica Bold, 20"

set xrange [AX[ARG1 + 1]:BX[ARG1 + 1]]
set yrange [AY[ARG1 + 1]:BY[ARG1 + 1]]

splot @function_f title "φ(x, y)" nosurface, \
      @function_g lc rgb "orange" title "g(x, y)" nocontours, \
      trialfile index ind + 2 ls 5 lc rgb "green" title "trial points" nocontours, \
      trialfile index ind + 1 ls 5 lc rgb "blue" title "X" nocontours, \
      trialfile index ind ls 5 lc rgb "red" title "X*" nocontours

bind all "alt-End" "exit gnuplot"
pause mouse close
