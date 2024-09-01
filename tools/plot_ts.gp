#set terminal png font ',14' size 2000,1200
#set output 'ts_plot.png'
set terminal x11 size 1400,900

set style line 1 lw 1 lc rgb "blue"
set style line 2 lw 1 lc rgb "red"
set style line 3 lw 2 lc rgb "green"
set style line 4 lw 2 lc rgb "orange"
set style line 5 lw 2 lc rgb "purple"
set style line 6 lw 2 lc rgb "cyan"

set multiplot layout 3,2 rowsfirst

set title ' max w '
set xlabel ' time [s] '
set ylabel ' w [m/s] '
set yrange [0:80]

plot "./OUTPUT/T_S" using 1:2 with lines ls 1 notitle

set title ' mean prec '
set xlabel ' time [s] '
set ylabel ' prec [mm/d] '
set yrange [0:120]

plot "./OUTPUT/T_S" using 1:7 with lines ls 2 notitle

set title ' qt mean '
set xlabel ' time [s] '
set ylabel ' qt [kg/m3] '
set yrange [0.0025:0.0033]

plot "./OUTPUT/T_S" using 1:14 with lines ls 3 notitle

set title ' pt mean '
set xlabel ' time [s] '
set ylabel ' pt [K] '
set yrange [338:345]

plot "./OUTPUT/T_S" using 1:16 with lines ls 4 notitle

set title ' SHF '
set xlabel ' time [s] '
set ylabel ' SHF [W/m2] '
set yrange [0:180]

plot "./OUTPUT/T_S" using 1:17 with lines ls 5 notitle

set title ' LHF '
set xlabel ' time [s] '
set ylabel ' LHF [W/m2] '
set yrange [0:500]

plot "./OUTPUT/T_S" using 1:18 with lines ls 6 notitle

unset multiplot
