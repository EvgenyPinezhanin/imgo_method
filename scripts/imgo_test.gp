#!/usr/bin/gnuplot

taskName = "imgo_test"
trialfile = "output_data/".taskName.".txt"
load "output_data/".taskName."_opt.txt"

f1(i, x) = (i == 0 ? -13.0 / 6.0 * x + sin(13.0 / 4.0 * (2.0 * x + 5.0)) - 53.0 / 12.0 : \
            i == 1 ? exp(-sin(3.0 * x)) - 1.0 / 10.0 * (x - 1.0 / 2.0) ** 2 - 1.0 : 1/0)

f2(i, x) = (i == 0 ? (11.0 * x * x - 10.0 * x + 21.0) / (2.0 * (x * x + 1)) : \
            i == 1 ? 1.0 / 20.0 - exp(-2.0 / 5.0 * (x + 5.0)) * sin(4.0 / 5.0 * pi * (x + 5.0)) : 1/0)

f3(i, x) = (i == 0 ? -sum[j = 1 : 5] cos(j * x) : \
            i == 1 ? 3.0 / 2.0 * (cos(7.0 / 20.0 * (x + 10.0)) - sin(7.0 / 4.0 * (x + 10.0)) + 1.0 / 2.0) : 1/0)

f4_a(x) = pi / 2.0 * x + 1.0 / 10.0
f4(i, x) = (i == 0 ? 4.0 * sin(f4_a(x) / 2.0) * (sin(f4_a(x)) ** 3 + cos(f4_a(x)) ** 3) ** 2 : \
            i == 1 ? 6.0 / 25.0 - sum[j = 1 : 5] cos(5.0 / 4.0 * (j + 1.0) * x + j) : \
            i == 2 ? 9.0 / 50.0 - 9.0 / 2.0 * exp(-(x - 1.0 / 10.0)) * sin(2.0 * pi * (x - 1.0 / 10.0)) : 1/0)

f5(i, x) = (i == 0 ? sin(0.423531 * x + 3.13531) + sin(10.0 / 3.0 * (0.423531 * x + 3.13531)) + \
                     log(0.423531 * x + 3.13531) + 0.36634 - 0.355766 * x : \
            i == 1 ? 17.0 / 25.0 - 2.0 / 29763.233 * (-1.0 / 6.0 * x ** 6 + 52.0 / 25.0 * x ** 5 - 39.0 / 80.0 * x ** 4 - \
                     71.0 / 10.0 * x ** 3 + 79.0 / 20.0 * x ** 2 + x - 1.0 / 10.0) : \
            i == 2 ? -14.0 / 125.0 * (3.0 * x - 8.0) * sin(252.0 / 125.0 * (x + 3.0 / 2.0)) - 1.0 / 2.0 : 1/0)

f6(i, x) = (i == 0 ? -7.0 / 40.0 * (3.0 * x + 4.0) * sin(63.0 / 20.0 * (x + 4.0)) : \
            i == 1 ? 40.0 * (cos(4.0 * x) * (x - sin(x)) * exp(-(x * x) / 2.0)) : \
            i == 2 ? 2.0 / 25.0 * (x + 4.0) - sin(12.0 / 5.0 * (x + 4.0)) : 1/0)

f7(i, x) = (i == 0 ? exp(-cos(4.0 * x - 3.0)) + 1.0 / 250.0 * (4.0 * x - 3.0) ** 2 - 1.0 : \
            i == 1 ? sin(x) ** 3 * exp(-sin(3.0 * x)) + 1.0 / 2.0 : \
            i == 2 ? cos(7.0 / 5.0 * (x + 3.0)) - sin(7.0 * (x + 3.0)) + 3.0 / 10.0 : 1/0)

f8(i, x) = (i == 0 ? cos(7.0 / 4.0 * x + 241.0 / 40.0) - sin(35.0 / 4.0 * x + 241.0 / 8.0) - 5.0 : \
            i == 1 ? exp(-sin(4.0 * x)) - 1.0 / 10.0 * (x - 1.0 / 2.0) ** 2 - 1.0 : \
            i == 2 ? 3.0 / 10.0 - sum[j = 1 : 5] cos(5.0 * (j + 1.0) * (x + 1.0 / 2.0)) : \
            i == 3 ? (-21.0 / 20.0 * x - 13.0 / 8.0) * sin(63.0 / 10.0 * x + 63.0 / 4.0) + 1.0 / 5.0 : 1/0)

f9(i, x) = (i == 0 ? sum[j = 1:5] 1.0 / 5.0 * sin((j + 1.0) * x - 1.0) + 2.0 : \
            i == 1 ? 1.0 / 40.0 * (x - 4.0) * (x - 32.0 / 5.0) * (x - 9.0) * (x - 11.0) * \
                     exp(-1.0 / 10.0 * (x - 13.0 / 2.0) ** 2) : \
            i == 2 ? (sin(x + 1.0) ** 3 + cos(x + 1.0) ** 3) * exp(-(x + 1.0) / 10.0) : \
            i == 3 ? exp(-cos(3.0 / 5.0 * (x - 5.0 / 2.0))) + 1.0 / 10.0 * (3.0 / 25.0 * x - 4.0 / 5.0) ** 2 - 1.0 : 1/0)

f10_a(x) = 2.0 / pi * x - 1.0 / 2.0
f10_b(x) = 4.0 / pi * (x - 3.0 / 10.0) - 4.0
f10(i, x) = (i == 0 ? -1.0 / 500.0 * f10_b(x) ** 6 + 3.0 / 100.0 * f10_b(x) ** 4 - 27.0 / 500.0 * f10_b(x) ** 2 + 3.0 / 2.0 : \
             i == 1 ? 2.0 * exp(-2.0 / pi * x) * sin(4.0 * x) : \
             i == 2 ? -f10_a(x) ** 2 * (-f10_a(x) ** 2 + 5.0 * f10_a(x) - 6.0) / (f10_a(x) ** 2 + 1.0) - 1.0 / 2.0 : \
             i == 3 ? sin(x) ** 3 + cos(2.0 * x) ** 3 - 3.0 / 10.0 : 1/0)

f(i, j, x) = (i == 1 ? f1(j, x) : \
              i == 2 ? f2(j, x) : \
              i == 3 ? f3(j, x) : \
              i == 4 ? f4(j, x) : \
              i == 5 ? f5(j, x) : \
              i == 6 ? f6(j, x) : \
              i == 7 ? f7(j, x) : \
              i == 8 ? f8(j, x) : \
              i == 9 ? f9(j, x) : \
              i == 10 ? f10(j, x) : 1/0)

set grid

set xlabel "x"
set ylabel "f(x)"

set samples 600

title(n) = sprintf("Graph of the test function №%d, method imgo", n)
titlePng(n) = (ARG1 == 1) ? title(n) : sprintf("")
functionName(i) = (i == 0 ? sprintf("f(x)") : sprintf("g%d(x)", i))

if (ARG1 == 0) {
    ind = 3 * (ARG2 - 1)

    set title title(int(ARG2)) font "Helvetica Bold, 20"

    plot for [i = 0 : numberConstraints[int(ARG2)]] f(ARG2, i, x) title functionName(i), \
         trialfile index ind + 2 ls 4 lc rgb "green" lw 2 title "trials", \
         trialfile index ind     ls 7 lc rgb "red"   ps 2 title "X_{min}", \
         trialfile index ind + 1 ls 7 lc rgb "blue"  ps 1 title "X"

    bind all "alt-End" "exit gnuplot"
    pause mouse close
} else {
    set terminal pngcairo size 1440, 800 font "Helvetica Bold, 15"
    system "mkdir -p output_graph/".taskName

    do for [i = 1 : 10] {
        set output "output_graph/".taskName."/".taskName."_".i.".png"

        ind = 3 * (i - 1)

        set title titlePng(i) font "Helvetica Bold, 15"

        plot for [j = 0 : numberConstraints[i]] f(i, j, x) title functionName(j), \
             trialfile index ind + 2 ls 4 lc rgb "green" lw 2 title "trials", \
             trialfile index ind     ls 7 lc rgb "red"   lw 6 title "X_{min}", \
             trialfile index ind + 1 ls 7 lc rgb "blue"  lw 1 title "X"
    }
}
