#! /usr/bin/gnuplot

reset

datafile = "lipschitz_test.txt"

set linetype 1 lc rgb "#0000FF" lw 5 pt 1 dt 1
set linetype 2 lc rgb "#0020FF" lw 5 pt 1 dt 2
set linetype 3 lc rgb "#0040FF" lw 5 pt 1 dt 3
set linetype 4 lc rgb "#0060FF" lw 4 pt 1 dt 4
set linetype 5 lc rgb "#0080FF" lw 4 pt 1 dt 5
set linetype 6 lc rgb "#00A0FF" lw 4 pt 1 dt 4
set linetype 7 lc rgb "#00C0FF" lw 3 pt 1 dt 3
set linetype 8 lc rgb "#00E0FF" lw 3 pt 1 dt 2
set linetype 9 lc rgb "#00FFFF" lw 3 pt 1 dt 1

# переделать
set linetype 10 lc rgb "#0FFF00" lw 4 pt 1 dt 1
set linetype 11 lc rgb "#2FFF00" lw 4 pt 1 dt 2
set linetype 12 lc rgb "#4FFF00" lw 4 pt 1 dt 3
set linetype 13 lc rgb "#6FFF00" lw 3 pt 1 dt 4
set linetype 14 lc rgb "#8FFF00" lw 3 pt 1 dt 5
set linetype 15 lc rgb "#AFFF00" lw 3 pt 1 dt 4
set linetype 16 lc rgb "#CFFF00" lw 2 pt 1 dt 3
set linetype 17 lc rgb "#EFFF00" lw 2 pt 1 dt 2
set linetype 18 lc rgb "#FFFF00" lw 2 pt 1 dt 2

set linetype 19 lc rgb "#FF000F" lw 3 pt 1 dt 1
set linetype 20 lc rgb "#FF002F" lw 3 pt 1 dt 2
set linetype 21 lc rgb "#FF004F" lw 3 pt 1 dt 3
set linetype 22 lc rgb "#FF006F" lw 2 pt 1 dt 4
set linetype 23 lc rgb "#FF008F" lw 2 pt 1 dt 5
set linetype 24 lc rgb "#FF00AF" lw 2 pt 1 dt 4
set linetype 25 lc rgb "#FF00CF" lw 1 pt 1 dt 3
set linetype 26 lc rgb "#FF00EF" lw 1 pt 1 dt 2
set linetype 27 lc rgb "#FF00FF" lw 1 pt 1 dt 1

set linetype cycle 27

set xlabel "incr"
set ylabel "lipschitz"
set xrange [ARG6:ARG7]

if (ARG1 == 0) {
    set title 'lipschitz test(Grishagin)'
    plot for [i = ARG2:ARG3] for [j = ARG4:ARG5] datafile index (ARG5 - ARG4 + 1) * (i - 1) + j - ARG4 using 1:2 with lines \
         lt 9 * (i - 1) + j - ARG4 + 1 title sprintf("key = %d, m = %d", i, j)
}
if (ARG1 == 1) {
    set title 'lipschitz test(GKLS)'
    plot for [i = ARG2:ARG3] for [j = ARG4:ARG5] datafile index (ARG5 - ARG4 + 1) * ARG3 + (ARG5 - ARG4 + 1) * (i - 1) + j - ARG4 using 1:2 with lines \
         lt 9 * (i - 1) + j - ARG4 + 1 title sprintf("key = %d, m = %d", i, j)
}
if (ARG1 == 2) {
    set title 'lipschitz test(GrishaginConstrained)'
    plot for [i = ARG2:ARG3] for [j = ARG4:ARG5] datafile index (ARG5 - ARG4 + 1) * ARG3 * 2 + (ARG5 - ARG4 + 1) * (i - 1) + j - ARG4 using 1:2 with lines \
         lt 9 * (i - 1) + j - ARG4 + 1 title sprintf("key = %d, m = %d", i, j)
}
if (ARG1 == 3) {
    set title 'lipschitz test(GKLSConstrained)'
    plot for [i = ARG2:ARG3] for [j = ARG4:ARG5] datafile index (ARG5 - ARG4 + 1) * ARG3 * 3 + (ARG5 - ARG4 + 1) * (i - 1) + j - ARG4 using 1:2 with lines \
         lt 9 * (i - 1) + j - ARG4 + 1 title sprintf("key = %d, m = %d", i, j)
}

bind all "alt-End" "exit gnuplot"
pause mouse close
