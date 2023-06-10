#!/usr/bin/gnuplot

taskName = "gsa_sample"
trialfile = "output_data/".taskName.".txt"

f1(x) = -4.0 * x + 1.0

f2(x) = 5.0 * x ** 2 + 3.0 * x - 1.0

f3(x) = x * sin(x)

f4(x) = (x != 0) ? x * sin(1 / x) : 0.0

f(i, x) = (i == 1 ? f1(x) : \
           i == 2 ? f2(x) : \
           i == 3 ? f3(x) : \
           i == 4 ? f4(x) : 1/0)

set grid

set xlabel "x"
set ylabel "f(x)"

set samples 400

title(n) = sprintf("Graph of the sample function №%d, method gsa", n)
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

    do for [i = 1 : 4] {
        set output "output_graph/".taskName."/".taskName."_".i.".png"

        ind = 3 * (i - 1)

        set title titlePng(i) font "Helvetica Bold, 15"

        plot f(i, x) title "f(x)", \
             trialfile index ind + 2 ls 4 lc rgb "green" lw 2 title "trials", \
             trialfile index ind     ls 7 lc rgb "red"   lw 6 title "X_{min}", \
             trialfile index ind + 1 ls 7 lc rgb "blue"  lw 1 title "X"
    }
}
