#set terminal png font ',14' size 1200,900
#set output 'pdf_plot.png'
set terminal x11 0 size 1200,900

reset

set style line 10 lw 2 lc rgb "black"
set style line 1 lw 1 lc rgb "blue"
set style line 2 lw 1 lc rgb "red"
set style line 3 lw 1 lc rgb "green"
set style line 4 lw 1 lc rgb "orange"
set style line 5 lw 1 lc rgb "purple"
set style line 6 lw 1 lc rgb "cyan"

set multiplot layout 2,2 rowsfirst

set title ' cloud size pdf '
set xlabel ' bin '
set ylabel ' PDF [#] '
set yrange [0:20]
set xrange [0:50]

plot "./CL_PDF2400" using 1:2 smooth csplines with boxes notitle

set title ' cloud mass flux pdf '
set xlabel ' bin '
set ylabel ' PDF [#] '
set yrange [0:8]
set xrange [0:50]

plot "./CL_PDF2400" using 1:3 with boxes ls 2 notitle, \
     "./CL_PDF2400" using 1:3 smooth acsplines with lines ls 10 notitle

set title ' updraft mass flux pdf '
set xlabel ' bin '
set ylabel ' PDF [#] '
set yrange [0:20]
set xrange [0:50]

plot "./CL_PDF2400" using 1:4 with boxes ls 3 notitle, \
     "./CL_PDF2400" using 1:4 smooth acsplines with lines ls 10 notitle

set title ' (renormalized) distance pdf '
set xlabel ' distance [km] '
set ylabel ' PDF [#] '
set yrange [0:5]
set xrange [0:50]

plot "./CL_PDF2400" using ($1*2):8 smooth acsplines with boxes ls 5 notitle, \
     "./CL_PDF2400" using ($1*2):9 smooth acsplines with boxes ls 6 notitle

unset multiplot
