#!/bin/python

file="e9_vpv.txt"

set size ratio -1 1

unset hidden3d
unset grid
set pm3d map
set colorbox
set pm3d corners2color mean
set palette rgbformulae 30,31,32
set title "vpv z=0"
set xlabel "x[m]"
set ylabel "y[m]"
set cblabel "Ax [Tm]"
set autoscale     
#set cbrange [2:22]

splot file using ($1):($2):($3)  notitle

set terminal wxt
replot

pause -1 "Press enter to continue:"

reset

#############################################################
a=1;
file="e9_vpv1.txt"

#set size ratio -1 1

set style data linespoints

set title "vpv z=0"
set xlabel "x[m]"
set ylabel "Ax [Tm]"
set autoscale     

#set logscale

plot file using ($1):($2*a) lt 1 title "new" , \
		 "e9_vpv1.txt" using ($1):($3*a) lt 2 title "old"


set terminal wxt
replot

pause -1 "Press enter to continue:"

reset
#############################################################

set terminal wxt
#set size ratio 1 1.1

set style data linespoints

set title "vector potential"
set xlabel "x[m]"
set ylabel "Ax [Tm]"
#set logscale y
set autoscale 

plot "e9_vpvv.txt" using ($1):($2) lt 1 title "old" , \
	   "e9_vpvv.txt" using ($1):($3) lt 2 title "new"

set terminal wxt
replot

pause -1 "press enter"
reset


