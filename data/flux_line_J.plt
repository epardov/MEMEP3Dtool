#!/bin/bash

filepar="par_plot.plt"
filein="flux_line_J.txt"

load filepar
f= 1/T;
########################################################################

set terminal jpeg
set size ratio -1
set output "flux_line_J.jpeg"

#set size ratio -1
#set terminal postscript eps enhanced color lw 1 "Arial" 14
#set output "flux_line_J.eps"

set xlabel "x"
set ylabel "y"
set cblabel "flux lines J"
set key off
set style data linespoints

lwidth=1
set samples 200
set style data lines 
set xrange[0:0.01]
set yrange[0:0.01]            

#plot "flux_line_J.txt" using ($2):($3) lt 1 title "point"

#set term wxt
#replot
#pause -1 "Press ENTER:"

#reset


do for [i=0:step-1]{
		t= (i+1)*dt;
		Ba= Bamax * sin (2 * pi * f * t)*1e3;
		set title sprintf("Ba=%f [mT]",Ba)
		plot filein index i using ($2):($3) lt 1 title "point"
}  

reset
