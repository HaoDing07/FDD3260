#set terminal png font ',14' size 1200,900
#set output 'jpdf_plot.png'
set terminal x11 2 size 1200,900

set pm3d

set palette rgbformulae 7,5,15

set title ' Parcel trajectory '
set xlabel ' X [km] '
set ylabel ' Y [km] '
set zlabel ' Z [km] '
set cblabel ' Potential temperature [K] '
set xrange [0:128]
set yrange [0:128]
set zrange [0:2.5]

splot "./OUTPUT/PARCEL120" using ($2/1000):($3/1000):($4/1000):7 with lines lw 2 palette notitle, \
"./OUTPUT/PARCEL121" using ($2/1000):($3/1000):($4/1000):7 with lines lw 2 palette notitle, \
"./OUTPUT/PARCEL122" using ($2/1000):($3/1000):($4/1000):7 with lines lw 2 palette notitle, \
"./OUTPUT/PARCEL123" using ($2/1000):($3/1000):($4/1000):7 with lines lw 2 palette notitle, \
"./OUTPUT/PARCEL124" using ($2/1000):($3/1000):($4/1000):7 with lines lw 2 palette notitle, \
"./OUTPUT/PARCEL125" using ($2/1000):($3/1000):($4/1000):7 with lines lw 2 palette notitle, \
"./OUTPUT/PARCEL126" using ($2/1000):($3/1000):($4/1000):7 with lines lw 2 palette notitle, \
"./OUTPUT/PARCEL127" using ($2/1000):($3/1000):($4/1000):7 with lines lw 2 palette notitle, \
"./OUTPUT/PARCEL128" using ($2/1000):($3/1000):($4/1000):7 with lines lw 2 palette notitle, \
"./OUTPUT/PARCEL129" using ($2/1000):($3/1000):($4/1000):7 with lines lw 2 palette notitle, \
"./OUTPUT/PARCEL130" using ($2/1000):($3/1000):($4/1000):7 with lines lw 2 palette notitle, \
"./OUTPUT/PARCEL131" using ($2/1000):($3/1000):($4/1000):7 with lines lw 2 palette notitle, \
"./OUTPUT/PARCEL132" using ($2/1000):($3/1000):($4/1000):7 with lines lw 2 palette notitle, \
"./OUTPUT/PARCEL133" using ($2/1000):($3/1000):($4/1000):7 with lines lw 2 palette notitle, \
"./OUTPUT/PARCEL134" using ($2/1000):($3/1000):($4/1000):7 with lines lw 2 palette notitle, \
"./OUTPUT/PARCEL135" using ($2/1000):($3/1000):($4/1000):7 with lines lw 2 palette notitle, \
"./OUTPUT/PARCEL136" using ($2/1000):($3/1000):($4/1000):7 with lines lw 2 palette notitle, \
"./OUTPUT/PARCEL137" using ($2/1000):($3/1000):($4/1000):7 with lines lw 2 palette notitle, \
"./OUTPUT/PARCEL138" using ($2/1000):($3/1000):($4/1000):7 with lines lw 2 palette notitle, \
"./OUTPUT/PARCEL139" using ($2/1000):($3/1000):($4/1000):7 with lines lw 2 palette notitle, \
"./OUTPUT/PARCEL140" using ($2/1000):($3/1000):($4/1000):7 with lines lw 2 palette notitle, \
"./OUTPUT/PARCEL141" using ($2/1000):($3/1000):($4/1000):7 with lines lw 2 palette notitle, \
"./OUTPUT/PARCEL142" using ($2/1000):($3/1000):($4/1000):7 with lines lw 2 palette notitle, \
"./OUTPUT/PARCEL143" using ($2/1000):($3/1000):($4/1000):7 with lines lw 2 palette notitle, \
"./OUTPUT/PARCEL144" using ($2/1000):($3/1000):($4/1000):7 with lines lw 2 palette notitle, \
"./OUTPUT/PARCEL145" using ($2/1000):($3/1000):($4/1000):7 with lines lw 2 palette notitle, \
"./OUTPUT/PARCEL146" using ($2/1000):($3/1000):($4/1000):7 with lines lw 2 palette notitle, \
"./OUTPUT/PARCEL147" using ($2/1000):($3/1000):($4/1000):7 with lines lw 2 palette notitle, \
"./OUTPUT/PARCEL148" using ($2/1000):($3/1000):($4/1000):7 with lines lw 2 palette notitle, \
"./OUTPUT/PARCEL149" using ($2/1000):($3/1000):($4/1000):7 with lines lw 2 palette notitle, \
"./OUTPUT/PARCEL150" using ($2/1000):($3/1000):($4/1000):7 with lines lw 2 palette notitle, \
"./OUTPUT/PARCEL151" using ($2/1000):($3/1000):($4/1000):7 with lines lw 2 palette notitle, \
"./OUTPUT/PARCEL152" using ($2/1000):($3/1000):($4/1000):7 with lines lw 2 palette notitle, \
"./OUTPUT/PARCEL153" using ($2/1000):($3/1000):($4/1000):7 with lines lw 2 palette notitle, \
"./OUTPUT/PARCEL154" using ($2/1000):($3/1000):($4/1000):7 with lines lw 2 palette notitle, \
"./OUTPUT/PARCEL155" using ($2/1000):($3/1000):($4/1000):7 with lines lw 2 palette notitle, \
"./OUTPUT/PARCEL156" using ($2/1000):($3/1000):($4/1000):7 with lines lw 2 palette notitle, \
"./OUTPUT/PARCEL157" using ($2/1000):($3/1000):($4/1000):7 with lines lw 2 palette notitle, \
"./OUTPUT/PARCEL158" using ($2/1000):($3/1000):($4/1000):7 with lines lw 2 palette notitle, \
"./OUTPUT/PARCEL159" using ($2/1000):($3/1000):($4/1000):7 with lines lw 2 palette notitle, \
"./OUTPUT/PARCEL160" using ($2/1000):($3/1000):($4/1000):7 with lines lw 2 palette notitle, \
"./OUTPUT/PARCEL161" using ($2/1000):($3/1000):($4/1000):7 with lines lw 2 palette notitle, \
"./OUTPUT/PARCEL162" using ($2/1000):($3/1000):($4/1000):7 with lines lw 2 palette notitle, \
"./OUTPUT/PARCEL163" using ($2/1000):($3/1000):($4/1000):7 with lines lw 2 palette notitle, \
"./OUTPUT/PARCEL164" using ($2/1000):($3/1000):($4/1000):7 with lines lw 2 palette notitle, \
"./OUTPUT/PARCEL165" using ($2/1000):($3/1000):($4/1000):7 with lines lw 2 palette notitle, \
"./OUTPUT/PARCEL166" using ($2/1000):($3/1000):($4/1000):7 with lines lw 2 palette notitle, \
"./OUTPUT/PARCEL167" using ($2/1000):($3/1000):($4/1000):7 with lines lw 2 palette notitle, \
"./OUTPUT/PARCEL168" using ($2/1000):($3/1000):($4/1000):7 with lines lw 2 palette notitle, \
"./OUTPUT/PARCEL169" using ($2/1000):($3/1000):($4/1000):7 with lines lw 2 palette notitle, \
"./OUTPUT/PARCEL170" using ($2/1000):($3/1000):($4/1000):7 with lines lw 2 palette notitle, \
"./OUTPUT/PARCEL171" using ($2/1000):($3/1000):($4/1000):7 with lines lw 2 palette notitle, \
"./OUTPUT/PARCEL172" using ($2/1000):($3/1000):($4/1000):7 with lines lw 2 palette notitle, \
"./OUTPUT/PARCEL173" using ($2/1000):($3/1000):($4/1000):7 with lines lw 2 palette notitle, \
"./OUTPUT/PARCEL174" using ($2/1000):($3/1000):($4/1000):7 with lines lw 2 palette notitle, \
"./OUTPUT/PARCEL175" using ($2/1000):($3/1000):($4/1000):7 with lines lw 2 palette notitle, \
"./OUTPUT/PARCEL176" using ($2/1000):($3/1000):($4/1000):7 with lines lw 2 palette notitle, \
"./OUTPUT/PARCEL177" using ($2/1000):($3/1000):($4/1000):7 with lines lw 2 palette notitle, \
"./OUTPUT/PARCEL178" using ($2/1000):($3/1000):($4/1000):7 with lines lw 2 palette notitle, \
"./OUTPUT/PARCEL179" using ($2/1000):($3/1000):($4/1000):7 with lines lw 2 palette notitle, \
"./OUTPUT/PARCEL180" using ($2/1000):($3/1000):($4/1000):7 with lines lw 2 palette notitle, \
"./OUTPUT/PARCEL181" using ($2/1000):($3/1000):($4/1000):7 with lines lw 2 palette notitle, \
"./OUTPUT/PARCEL182" using ($2/1000):($3/1000):($4/1000):7 with lines lw 2 palette notitle, \
"./OUTPUT/PARCEL183" using ($2/1000):($3/1000):($4/1000):7 with lines lw 2 palette notitle, \
"./OUTPUT/PARCEL184" using ($2/1000):($3/1000):($4/1000):7 with lines lw 2 palette notitle, \
"./OUTPUT/PARCEL185" using ($2/1000):($3/1000):($4/1000):7 with lines lw 2 palette notitle, \
"./OUTPUT/PARCEL186" using ($2/1000):($3/1000):($4/1000):7 with lines lw 2 palette notitle, \
"./OUTPUT/PARCEL187" using ($2/1000):($3/1000):($4/1000):7 with lines lw 2 palette notitle, \
"./OUTPUT/PARCEL188" using ($2/1000):($3/1000):($4/1000):7 with lines lw 2 palette notitle, \
"./OUTPUT/PARCEL189" using ($2/1000):($3/1000):($4/1000):7 with lines lw 2 palette notitle, \
"./OUTPUT/PARCEL190" using ($2/1000):($3/1000):($4/1000):7 with lines lw 2 palette notitle, \
"./OUTPUT/PARCEL191" using ($2/1000):($3/1000):($4/1000):7 with lines lw 2 palette notitle, \
"./OUTPUT/PARCEL192" using ($2/1000):($3/1000):($4/1000):7 with lines lw 2 palette notitle, \
"./OUTPUT/PARCEL193" using ($2/1000):($3/1000):($4/1000):7 with lines lw 2 palette notitle, \
"./OUTPUT/PARCEL194" using ($2/1000):($3/1000):($4/1000):7 with lines lw 2 palette notitle, \
"./OUTPUT/PARCEL195" using ($2/1000):($3/1000):($4/1000):7 with lines lw 2 palette notitle, \
"./OUTPUT/PARCEL196" using ($2/1000):($3/1000):($4/1000):7 with lines lw 2 palette notitle, \
"./OUTPUT/PARCEL197" using ($2/1000):($3/1000):($4/1000):7 with lines lw 2 palette notitle, \
"./OUTPUT/PARCEL198" using ($2/1000):($3/1000):($4/1000):7 with lines lw 2 palette notitle, \
"./OUTPUT/PARCEL199" using ($2/1000):($3/1000):($4/1000):7 with lines lw 2 palette notitle, \
"./OUTPUT/PARCEL200" using ($2/1000):($3/1000):($4/1000):7 with lines lw 2 palette notitle
