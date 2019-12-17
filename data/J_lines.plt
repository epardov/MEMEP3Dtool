#!/bin/python

filepar="par_plot.plt"
filein="output3Dz.txt"

########################################################################

load filepar
w=x*1000;
t= (step+1)*dt;
Ba= Bamax * sin (2 * pi * f * t)*1e3;
a= z/(ncz);

i=0;
c=step-1;
c=18;

T0=1e8;
h=1.2;
g=Bamax*1000;
if(rel==3){Jo=Jcpe};


########################################################################

#Macros

set macros

#epsterm='set terminal postscript eps enhanced color lw 1 "Arial" 28'
epsterm='set terminal jpeg
########################################################################

do for [n=0:c]{
set size ratio -1

set terminal wxt
set tics

set xlabel "x [mm]" 
set ylabel "y [mm]" 
set cblabel "|J|/J_{c,per}"
set key off

set view map
set contour 
set nosurface
set cntrparam levels incremental -1.0,0.02, 1.0

fileJlines=sprintf("data_out/Jlines%d.txt",n)
set style data lines
lwidth=1
lt=-11
lc=1
#set dgrid3d
set dgrid3d 219,73
ts= (n+1)*dt;
b= Bamax * sin (2 * pi * f * ts);
set title sprintf("B_a=%.1f [mT]",b*1000) 

set term x11
set table
set out fileJlines
splot	filein index n using (($2*1000-w/2)/w):(($3*1000-(3*w/2))/w):($24/T0) with lines
unset table

@epsterm 
#set output sprintf("3Dz_J%.0f.eps",n)
set output sprintf("3Dz_J%.0f.jpeg",n)

set xlabel "x/w" 
set ylabel "y/w" 
set cblabel "|J|/J_{c,per}"
set key off

set palette rgbformulae 30,31,32
lwidth=1
set samples 200
set style data lines
set cbrange [-h	:h] 
set colorbox 
set xtics 0.5
set ytics 0.5	

#set xrange [0:x*1000]
#set yrange [0:y*1000]

set xrange [(-x*1000/2)/w:(x*1000/2)/w]
set yrange [(-y*1000/2)/w:(y*1000/2)/w]


ts= (n+1)*dt;
b= Bamax * sin (2 * pi * f * ts);
set title sprintf("B_a=%.1f [mT]",b*1000) 
plot \
filein index n using (($2*1000-w/2)/w):(($3*1000-(3*w/2))/w):($9/Jo) with image  notitle , \
fileJlines using 1:2 ls 1 lc 0 notitle
}
