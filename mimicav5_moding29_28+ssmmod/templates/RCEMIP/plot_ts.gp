#set terminal png font ',14' size 2000,1800
#set output 'ts_plot.png'
set terminal x11 size 1400,900

set style line 1 lw 1 lc rgb "blue"
set style line 2 lw 1 lc rgb "red"
set style line 3 lw 2 lc rgb "sea-green"
set style line 4 lw 2 lc rgb "gold"
set style line 5 lw 2 lc rgb "coral"
set style line 6 lw 2 lc rgb "cyan"
set style line 7 lw 2 lc rgb "dark-violet"
set style line 8 lw 2 lc rgb "grey60"

set multiplot layout 4,2 rowsfirst

set title ' max w '
set xlabel ' time [s] '
set ylabel ' w [m/s] '
set yrange [0:75]

plot "./OUTPUT/T_S" using 1:2 with lines ls 1 notitle

set title ' mean prec '
set xlabel ' time [s] '
set ylabel ' prec [mm/d] '
set yrange [0:20]

plot "./OUTPUT/T_S" using 1:7 with lines ls 2 notitle

set title ' T tropopause '
set xlabel ' time [s] '
set ylabel ' T [K] '
set yrange [190:200]

plot "./OUTPUT/T_S" using 1:16 with lines ls 3 notitle

set title ' Z tropopause '
set xlabel ' time [s] '
set ylabel ' Z [m] '
set yrange [16500:17500]

plot "./OUTPUT/T_S" using 1:17 with lines ls 4 notitle

set title ' Toposphere O3 '
set xlabel ' time [s] '
set ylabel ' O3 [kg/m2] '
set yrange [0:0.01]

plot "./OUTPUT/T_S" using 1:($18*1e-6) with lines ls 5 notitle

set title ' CT height '
set xlabel ' time [s] '
set ylabel ' H [m] '
set yrange [0:18000]

plot "./OUTPUT/T_S" using 1:4 with lines ls 6 notitle

set title ' Surface fluxes '
set xlabel ' time [s] '
set ylabel ' SHF+LHF [W/m2] '
set yrange [0:250]

plot "./OUTPUT/T_S" using 1:($19+$20) with lines ls 7 notitle

set title ' RADF '
set xlabel ' time [s] '
set ylabel ' RAD [W/m2] '
set yrange [-200:0]

plot "./OUTPUT/T_S" using 1:21 with lines ls 8 notitle

unset multiplot
