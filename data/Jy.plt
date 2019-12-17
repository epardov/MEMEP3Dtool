#!/bin/python

filepar="par_plot.plt"
filein="Jy.txt"

########################################################################

load filepar
t= (step+1)*dt;
Ba= Bamax * sin (2 * pi * f * t)*1e3;
a= y/(ncy);

w=x*1000;
d=z*1000;

n=1;

set xrange [-x*1000/2/w:x*1000/2/w]
set yrange [-z*1000/2/d:z*1000/2/d]

g=Bamax*1000;
h=1.50;
j=0.3;

########################################################################

set size ratio -1
set terminal pdf color enhanced font 'Arial, 16'
set output sprintf("3Dy_dJy%.0f.pdf",n)

set xlabel "x/w"
set ylabel "z/d"
set cblabel "J_y/J_c"
set xtics 0.25
set ytics 0.25

set key off
unset colorbox
set style data lines 
set palette rgbformulae 30,31,32
lwidth=1
set samples 200
set style data lines     
set cbrange [-h:h]            

do for [i=0:ncy-1]{
b = a*(i+1) - (a/2);
set title sprintf("y/w=%.2f ",b*1000/w)
plot filein index i using (($3*1000-w/2)/w):(($4*1000-d/2)/d):($6/Jo) with image  notitle
} 




