#set terminal png font ',14' size 1200,900
#set output 'jpdf_plot.png'
set terminal x11 2 size 1200,900

set pm3d map interpolate 4,4

set palette rgbformulae 7,5,15

set title ' Joint distance-size PDF '
set xlabel ' distance [km] '
set ylabel ' cloud size [km2] '
set xrange [0:64]
set yrange [0:120]
set cbrange [0:0.08]

splot "./OUTPUT/CL_JPDF2400" using ($1/1000):($2/1000000):3 with pm3d notitle
