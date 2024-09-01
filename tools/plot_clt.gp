#set terminal png font ',14' size 1400,900
#set output 'clt_plot.png'
set terminal x11 1 size 1400,900

reset

set style line 1 lw 1 lc rgb "blue"
set style line 2 lw 1 lc rgb "red"
set style line 3 lw 2 lc rgb "green"
set style line 4 lw 2 lc rgb "orange"
set style line 5 lw 2 lc rgb "purple"
set style line 6 lw 2 lc rgb "cyan"

set multiplot layout 3,2 rowsfirst

set title ' Number of clouds '
set xlabel ' time [h] '
set ylabel ' N clouds [#] '
set yrange [40:90]

plot "./CL_TS2400" using ($1-1728000)/3600:3 with lines ls 1 notitle

set title ' Mean Cloud size '
set xlabel ' time [h] '
set ylabel ' size [km2] '
set yrange [0:80]

plot "./CL_TS2400" using ($1-1728000)/3600:($6*2*2) with lines ls 2 notitle

set title ' Mean updraft size '
set xlabel ' time [h] '
set ylabel ' size [km2] '
set yrange [0:20]

plot "./CL_TS2400" using ($1-1728000)/3600:($11*2*2) with lines ls 6 notitle

set title ' Mean cloud mass flux '
set xlabel ' time [h] '
set ylabel ' mass flux [kg/m2/s] '
set yrange [-5E+6:2E+7]

plot "./CL_TS2400" using ($1-1728000)/3600:8 with lines ls 3 notitle

set title ' Mean updraft mass flux '
set xlabel ' time [h] '
set ylabel ' mass flux [kg/m2/s] '
set yrange [0:2E+7]

plot "./CL_TS2400" using ($1-1728000)/3600:13 with lines ls 4 notitle

set title ' Mean cloud distance '
set xlabel ' time [h] '
set ylabel ' distance [km] '
set yrange [0:64]

plot "./CL_TS2400" using ($1-1728000)/3600:($15/1000.) with lines ls 5 notitle

unset multiplot
