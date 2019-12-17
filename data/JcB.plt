#!/bin/python

filein="loss.txt"
filepar="par_plot.plt"

load filepar

S=x*z

########################################################################

set size ratio 1
set terminal eps
set output "JcB.eps"
set xlabel "B [T]"
set ylabel "Ic[A]"
set key on
set style data linespoints

lwidth=1
set samples 200
set style data lines             
set grid

plot "JcB1.txt" using ($1):($2) lt 1 lc 1 title " "

########################################################################

set size ratio 1
set terminal eps
set output "nB.eps"
set xlabel "B [T]"
set ylabel "Ic[A]"
set key on
set style data linespoints
set yrange[0:30]


lwidth=1
set samples 200
set style data lines             
set grid

plot "JcB1.txt" using ($1):($3) lt 1 lc 1 title " "

########################################################################

