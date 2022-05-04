#!/usr/bin/gnuplot

reset

trialfile="peano_test_trial_points.txt"

f_1(x, y) = -1.5*x**2*exp(1-x**2-20.25*(x-y)**2)-(0.5*(x-1)*(y-1))**4*exp(2-(0.5*(x-1))**4-(y-1)**4)
g1_1(x, y) = 0.01*((x-2.2)**2+(y-1.2)**2.0-2.25)
g2_1(x, y) = 100*(1-((x-2)**2)/1.44-(0.5*y)**2)
g3_1(x, y) = 10*(y-1.5-1.5*sin(6.283*(x-1.75)))

g_1(x, y) = (g1_1(x, y) <= 0.0) ? (g2_1(x, y) <= 0 ? (g3_1(x, y) <= 0 ? g3_1(x, y) : 1/0) : 1/0) : 1/0

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

f_2(x, y) = -(C[1]+C[2]*x+C[3]*x*x+C[4]*x**3+C[5]*x**4+C[6]*y+C[7]*x*y+C[8]*x*x*y+C[9]*x**3*y+C[10]*x**4*y+ \
            C[11]*y*y+C[12]*y**3+C[13]*y**4+C[14]/(y+1)+C[15]*x*x*y*y+C[16]*x**3*y*y+C[17]*x**3*y**3+C[18]*x*y*y+ \
            C[19]*x*y**3+C[20]*exp(0.0005*x*y))
g1_2(x, y) = 450.0-x*y
g2_2(x, y) = (0.1*x-1.0)**2-y
g3_2(x, y) = 8.0*(x-40.0)-(y-30.0)*(y-55.0)
g4_2(x, y) = y+(x-35.0)*(x-30.0)/125.0-80.0

g_2(x, y) = (g1_2(x, y) <= 0.0) ? (g2_2(x, y) <= 0 ? (g3_2(x, y) <= 0 ? (g4_2(x, y) <= 0 ? (g4_2(x, y)) : 1/0) : 1/0) : 1/0) : 1/0

set grid
set contour
set view map
set cntrparam bspline levels auto 10
set cntrlabel onecolor
set cntrlabel start 5 interval 150
set cntrlabel font ",10"
set contour base
set isosamples 60

set xlabel "X"
set ylabel "Y"

if (ARG1 == 0) {
      set title "Chart of test function 1" font "Helvetica Bold, 20"
      set xrange [0.0:4.0]
      set yrange [-1.0:3.0]
      splot f_1(x, y) title "φ(x, y)" nosurface, \
            g_1(x, y) lc rgb "orange" notitle nocontours, \
            trialfile index 2 ls 5 lc rgb "green" title "trial points" nocontours, \
            trialfile index 0 ls 5 lc rgb "blue" title "X" nocontours, \
            trialfile index 1 ls 5 lc rgb "red" title "X*" nocontours, \
            f_1(x, y) with labels notitle nosurface
} else {
      set title "Chart of test function 2" font "Helvetica Bold, 20"
      set xrange [0.0:80.0]
      set yrange [0.0:80.0]
      splot f_2(x, y) title "φ(x, y)" nosurface, \
            g_2(x, y) lc rgb "orange" notitle nocontours, \
            trialfile index 5 ls 5 lc rgb "green" title "trial points" nocontours, \
            trialfile index 3 ls 5 lc rgb "blue" title "X" nocontours, \
            trialfile index 4 ls 5 lc rgb "red" title "X*" nocontours, \
            f_2(x, y) with labels notitle nosurface
}

bind all "alt-End" "exit gnuplot"
pause mouse close
