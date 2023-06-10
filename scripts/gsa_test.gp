#!/usr/bin/gnuplot

taskName = "gsa_test"
trialfile = "output_data/".taskName.".txt"

f1(x) = 1.0 / 6.0 * x ** 6 - 52.0 / 25.0 * x ** 5 + 39.0 / 80.0 * x ** 4 + \
        71.0 / 10.0 * x ** 3 - 79.0 / 20.0 * x ** 2 - x + 1.0 / 10.0

f2(x) = sin(x) + sin(10.0 / 3.0 * x)

f3(x) = -sum[i = 1 : 5] i * sin((i + 1.0) * x  + i)

f4(x) = (-16.0 * x ** 2 + 24.0 * x - 5.0) * exp(-x)

f5(x) = (3.0 * x - 1.4) * sin(18.0 * x)

f6(x) = -(x + sin(x)) * exp(-x ** 2)

f7(x) = sin(x) + sin(10.0 / 3.0 * x) + log(x) - 0.84 * x + 3.0

f8(x) = -sum[i = 1 : 5] i * cos((i + 1.0) * x  + i)

f9(x) = sin(x) + sin(2.0 / 3.0 * x)

f10(x) = -x * sin(x)

f11(x) = 2.0 * cos(x) + cos(2.0 * x)

f12(x) = sin(x) ** 3 + cos(x) ** 3

f13(x) = (x * x - 1 < 0) ? -x ** (2.0 / 3.0) - (-x ** 2 + 1) ** (1.0 / 3.0) : \
                           -x ** (2.0 / 3.0) + (x ** 2 - 1) ** (1.0 / 3.0)

f14(x) = -exp(-x) * sin(2 * pi * x)

f15(x) = (x ** 2 - 5.0 * x + 6.0) / (x * x + 1)

f16(x) = 2.0 * (x - 3) ** 2 + exp(-x * x / 2)

f17(x) = x ** 6 - 15.0 * x ** 4 + 27.0 * x ** 2 + 250.0

f18(x) = (x <= 3.0) ? (x - 2.0) ** 2 : 2.0 * log(x - 2.0) + 1.0

f19(x) = -x + sin(3.0 * x) - 1.0

f20(x) = (-x + sin(x)) * exp(-x ** 2)

f(i, x) = (i == 1 ? f1(x) : \
           i == 2 ? f2(x) : \
           i == 3 ? f3(x) : \
           i == 4 ? f4(x) : \
           i == 5 ? f5(x) : \
           i == 6 ? f6(x) : \
           i == 7 ? f7(x) : \
           i == 8 ? f8(x) : \
           i == 9 ? f9(x) : \
           i == 10 ? f10(x) : \
           i == 11 ? f11(x) : \
           i == 12 ? f12(x) : \
           i == 13 ? f13(x) : \
           i == 14 ? f14(x) : \
           i == 15 ? f15(x) : \
           i == 16 ? f16(x) : \
           i == 17 ? f17(x) : \
           i == 18 ? f18(x) : \
           i == 19 ? f19(x) : \
           i == 20 ? f20(x) : 1/0)

set grid

set xlabel "x"
set ylabel "f(x)"

set sample 400

title(n) = sprintf("Graph of the test function №%d, method gsa", n)
titlePng(n) = (ARG1 == 1) ? title(n) : sprintf("")

if (ARG1 == 0) {
    ind = 3 * (ARG2 - 1)

    set title title(int(ARG2)) font "Helvetica Bold, 20"

    plot f(ARG2, x) title "f(x)", \
         trialfile index ind + 2 ls 4 lc rgb "green" lw 2 title "trials", \
         trialfile index ind     ls 7 lc rgb "red"   ps 2 title "X_{min}", \
         trialfile index ind + 1 ls 7 lc rgb "blue"  ps 1 title "X"

    bind all "alt-End" "exit gnuplot"
    pause mouse close
} else {
    set terminal pngcairo size 1440, 800 font "Helvetica Bold, 15"
    system "mkdir -p output_graph/".taskName

    do for [i = 1 : 20] {
        set output "output_graph/".taskName."/".taskName."_".i.".png"

        ind = 3 * (i - 1)

        set title titlePng(i) font "Helvetica Bold, 15"

        plot f(i, x) title "f(x)", \
             trialfile index ind + 2 ls 4 lc rgb "green" lw 2 title "trials", \
             trialfile index ind     ls 7 lc rgb "red"   lw 6 title "X_{min}", \
             trialfile index ind + 1 ls 7 lc rgb "blue"  lw 1 title "X"
    }
}
