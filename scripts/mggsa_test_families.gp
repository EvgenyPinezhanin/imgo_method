#!/usr/bin/gnuplot

trialfile = "output_data/mggsa_test_families.txt"

load "output_data/mggsa_test_families_opt.txt"

array A[49]
array B[49]
array C[49]
array D[49]

A[1] = -0.391568
A[2] = -0.371071
A[3] = 0.421819
A[4] = 0.742908
A[5] = 0.851877
A[6] = -0.266324
A[7] = -0.294827
A[8] = -0.908095
A[9] = -0.00553646
A[10] = -0.507804
A[11] = -0.55012
A[12] = -0.0538802
A[13] = -0.812156
A[14] = 0.222869
A[15] = 0.816942
A[16] = 0.707656
A[17] = -0.626483
A[18] = -0.347937
A[19] = 0.784306
A[20] = 0.578099
A[21] = -0.727346
A[22] = 0.0788902
A[23] = 0.5571
A[24] = 0.924147
A[25] = 0.763278
A[26] = 0.850955
A[27] = -0.99829
A[28] = 0.426655
A[29] = 0.614509
A[30] = -0.966034
A[31] = -0.176037
A[32] = -0.715046
A[33] = 0.821229
A[34] = -0.982063
A[35] = 0.234552
A[36] = 0.511914
A[37] = -0.534883
A[38] = 0.0520894
A[39] = -0.101247
A[40] = -0.156848
A[41] = 0.894013
A[42] = 0.851274
A[43] = -0.260003
A[44] = -0.951941
A[45] = -0.609505
A[46] = 0.809819
A[47] = -0.114603
A[48] = -0.59492
A[49] = -0.39013

B[1] = -0.472228
B[2] = -0.481067
B[3] = -0.80654
B[4] = -0.130897
B[5] = -0.987827
B[6] = 0.732328
B[7] = 0.39572
B[8] = 0.515869
B[9] = -0.290875
B[10] = 0.581675
B[11] = 0.741153
B[12] = 0.449801
B[13] = 0.808699
B[14] = -0.883202
B[15] = 0.51777
B[16] = -0.0280453
B[17] = 0.745761
B[18] = -0.300433
B[19] = -0.0249191
B[20] = -0.531517
B[21] = 0.509613
B[22] = -0.364348
B[23] = -0.525241
B[24] = -0.202982
B[25] = 0.022882
B[26] = 0.250897
B[27] = 0.836626
B[28] = -0.0138833
B[29] = -0.857731
B[30] = -0.947244
B[31] = -0.421056
B[32] = 0.889776
B[33] = 0.939426
B[34] = -0.435605
B[35] = -0.542047
B[36] = -0.168104
B[37] = -0.656884
B[38] = 0.970905
B[39] = 0.98331
B[40] = -0.624119
B[41] = 0.617681
B[42] = 0.433995
B[43] = 0.120581
B[44] = 0.800537
B[45] = -0.286718
B[46] = 0.465421
B[47] = -0.937276
B[48] = 0.54388
B[49] = 0.982963

C[1] = 0.4799
C[2] = -0.868189
C[3] = -0.414678
C[4] = -0.337078
C[5] = -0.186449
C[6] = 0.323328
C[7] = -0.0327587
C[8] = 0.828282
C[9] = 0.703161
C[10] = 0.493221
C[11] = 0.900404
C[12] = 0.907243
C[13] = -0.76829
C[14] = 0.749824
C[15] = 0.252582
C[16] = 0.130364
C[17] = 0.566297
C[18] = 0.813088
C[19] = 0.675453
C[20] = -0.574807
C[21] = -0.360881
C[22] = 0.051283
C[23] = -0.243862
C[24] = 0.651978
C[25] = -0.560775
C[26] = 0.210204
C[27] = -0.77939
C[28] = 0.943773
C[29] = 0.208691
C[30] = -0.686298
C[31] = 0.645753
C[32] = -0.805916
C[33] = -0.202273
C[34] = -0.690077
C[35] = 0.225591
C[36] = -0.97876
C[37] = 0.999624
C[38] = 0.654434
C[39] = 0.922946
C[40] = 0.233188
C[41] = 0.546277
C[42] = 0.185159
C[43] = -0.979361
C[44] = -0.833939
C[45] = 0.594922
C[46] = -0.595781
C[47] = 0.711358
C[48] = 0.74251
C[49] = -0.452788

D[1] = 0.4799
D[2] = -0.868189
D[3] = -0.414678
D[4] = -0.337078
D[5] = -0.186449
D[6] = 0.323328
D[7] = -0.0327587
D[8] = 0.828282
D[9] = 0.703161
D[10] = 0.493221
D[11] = 0.900404
D[12] = 0.907243
D[13] = -0.76829
D[14] = 0.749824
D[15] = 0.252582
D[16] = 0.130364
D[17] = 0.566297
D[18] = 0.813088
D[19] = 0.675453
D[20] = -0.574807
D[21] = -0.360881
D[22] = 0.051283
D[23] = -0.243862
D[24] = 0.651978
D[25] = -0.560775
D[26] = 0.210204
D[27] = -0.77939
D[28] = 0.943773
D[29] = 0.208691
D[30] = -0.686298
D[31] = 0.645753
D[32] = -0.805916
D[33] = -0.202273
D[34] = -0.690077
D[35] = 0.225591
D[36] = -0.97876
D[37] = 0.999624
D[38] = 0.654434
D[39] = 0.922946
D[40] = 0.233188
D[41] = 0.546277
D[42] = 0.185159
D[43] = -0.979361
D[44] = -0.833939
D[45] = 0.594922
D[46] = -0.595781
D[47] = 0.711358
D[48] = 0.74251
D[49] = -0.452788

a(i, j, x, y) = sin(pi*i*x) * sin(pi*j*y)
b(i, j, x, y) = cos(pi*i*x) * cos(pi*j*y)

d1__(i, j, x, y) = A[7*(i-1) + j] * a(i, j, x, y) + B[7*(i-1) + j] * b(i, j, x, y)
d1_(i, x, y) = sum[j=1:7] d1__(i, j, x, y)
d1(x, y) = sum[i=1:7] d1_(i, x, y)

d2__(i, j, x, y) = C[7*(i-1) + j] * a(i, j, x, y) - D[7*(i-1) + j] * b(i, j, x, y)
d2_(i, x, y) = sum[j=1:7] d2__(i, j, x, y)
d2(x, y) = sum[i=1:7] d2_(i, x, y)

f(x, y) = -sqrt(d1(x, y)**2 + d2(x, y)**2)

array LOC_MIN[20]
array RHO[10]
array F[10]

LOC_MIN[1] = -0.762614
LOC_MIN[2] = 0.597254
LOC_MIN[3] = 0.0839592
LOC_MIN[4] = 0.902726
LOC_MIN[5] = 0.496543
LOC_MIN[6] = -0.939405
LOC_MIN[7] = 0.713418
LOC_MIN[8] = 0.627774
LOC_MIN[9] = -0.516797
LOC_MIN[10] = -0.605404
LOC_MIN[11] = -0.998932
LOC_MIN[12] = -0.459521
LOC_MIN[13] = 0.581651
LOC_MIN[14] = 0.54993
LOC_MIN[15] = -0.473927
LOC_MIN[16] = -0.911208
LOC_MIN[17] = 0.974159
LOC_MIN[18] = -0.021107
LOC_MIN[19] = -0.244438
LOC_MIN[20] = -0.587909

RHO[1] = 0.693
RHO[2] = 0.2
RHO[3] = 0.676827
RHO[4] = 0.0757565
RHO[5] = 0.135095
RHO[6] = 0.36359
RHO[7] = 0.0757565
RHO[8] = 0.170611
RHO[9] = 0.347901
RHO[10] = 0.135095

F[1] = 0
F[2] = -1
F[3] = 0.655211
F[4] = 1.87654
F[5] = 0.933122
F[6] = -0.0440105
F[7] = 1.52896
F[8] = 1.54059
F[9] = 1.58603
F[10] = 1.08011

dim = 2
num_minima = 11
precision = 1e-10

norma_(ind) = sqrt((LOC_MIN[1] - LOC_MIN[2 * (ind - 1) + 1])**2 + (LOC_MIN[2] - LOC_MIN[2 * (ind - 1) + 2])**2)
norma(ind, x, y) = sqrt((LOC_MIN[2 * (ind - 1) + 1] - x)**2 + (LOC_MIN[2 * (ind - 1) + 2] - y)**2)
scal(ind, x, y) = (x - LOC_MIN[2 * (ind - 1) + 1]) * (LOC_MIN[1] - LOC_MIN[2 * (ind - 1) + 1]) + (y - LOC_MIN[2 * (ind - 1) + 2]) * (LOC_MIN[2] - LOC_MIN[2 * (ind - 1) + 2])
index(x, y) = (norma(2, x, y) > RHO[2]) ? ((norma(3, x, y) > RHO[3]) ? ((norma(4, x, y) > RHO[4]) ? ((norma(5, x, y) > RHO[5]) ? ((norma(6, x, y) > RHO[6]) ? ((norma(7, x, y) > RHO[7]) ? ((norma(8, x, y) > RHO[8]) ? ((norma(9, x, y) > RHO[9]) ? ((norma(10, x, y) > RHO[10]) ? (11) : (10)) : (9)) : (8)) : (7)) : (6)) : (5)) : (4)) : (3)) : (2)
a(ind) = norma_(ind)**2 + F[1] - F[ind]

func(ind, x, y) = (2.0 / RHO[index(x, y)]**2 * scal(index(x, y), x, y) / norma(index(x, y), x, y) - 2.0 * a(index(x, y)) / RHO[index(x, y)]**3) * norma(index(x, y), x, y)**3 + (1.0 - 4.0 * scal(index(x, y), x, y) / (norma(index(x, y), x, y) * RHO[index(x, y)]) + 3.0 * a(index(x, y)) / RHO[index(x, y)]**2) * norma(index(x, y), x, y)**2 + F[index(x, y)]

f(x, y) = (index(x, y) == num_minima) ? (norma(1, x, y)**2 + F[1]) : ((norma(index(x, y), x, y) < precision) ? (F[index(x, y)]) : (func(index(x, y), x, y)))


titleName(n) = sprintf("Graph of the test function №%d", n)
functionName(s, n) = sprintf("%s%d(x, y)", s, n)

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
functionF = functionName("f", ARG1 + 1)
functionG = functionName("g", ARG1 + 1)

set title titleName(ARG1 + 1) font "Helvetica Bold, 20"

set xrange [AX[ARG1 + 1]:BX[ARG1 + 1]]
set yrange [AY[ARG1 + 1]:BY[ARG1 + 1]]

splot @functionF title "φ(x, y)" nosurface, \
      @functionG lc rgb "orange" title "g(x, y)" nocontours, \
      trialfile index ind + 2 ls 5 lc rgb "green" title "trial points" nocontours, \
      trialfile index ind + 1 ls 5 lc rgb "blue" title "X" nocontours, \
      trialfile index ind ls 5 lc rgb "red" title "X*" nocontours, \
      # @functionF with labels notitle nosurface

bind all "alt-End" "exit gnuplot"
pause mouse close
