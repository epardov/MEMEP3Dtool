#!/bin/python

filepar="par_plot.plt"
filein="output3Dz.txt"
#filein="outputAV.txt"
filein1="Tz.txt"

#input
#form:0-jpeg, 1-eps 
form=1;

########################################################################

load filepar
w=x*1000;

#w=w/2;##when its normalized to the radius


t= (step+1)*dt;
Ba= Bamax * sin (2 * pi * f * t)*1e3;
a= z/(ncz);
mi0=4*pi*1e-7;

n=0;

h=0.1;
g=Bamax*1000;
if(rel==3){Jo=Jcpe};


########################################################################

#Macros

set macros

if(form==0){epsterm='set terminal jpeg'}
else{epsterm='set terminal postscript eps enhanced color lw 1 "Arial" 28'}

########################################################################

set size ratio -1

set palette rgbformulae 30,31,32
lwidth=1
set samples 200
set style data lines
set cbrange [0:h] 
#unset colorbox 

set parametric
set contour surface
set view 0,0,1
unset surface
set cntrparam levels 15

set table
set out "data_out/contours.txt"
splot filein1 index n using ($1):($2):3 with line
unset table

@epsterm
if(form==0){set output sprintf("2Dz_J%.0f.jpeg",n)}
else{set output sprintf("2Dz_J%.0f.eps",n)}

set xlabel "x/w" 
set ylabel "y/w" 
#set xlabel "x/R" 
#set ylabel "y/R" 

set cblabel "|J|/J_{c,per}"
set key off

set xtics 0.5
#set ytics format " "
set ytics 0.5

set xrange [(-x*1000/2)/w:(x*1000/2)/w]
set yrange [(-y*1000/2)/w:(y*1000/2)/w]
#set xrange [0:x*1000]
#set yrange [0:y*1000]

ts= (n+1)*dt;
#b= Bamax * sin (2 * pi * f * ts);
b= Bamax * sin (2 * pi * f * ts)/(mi0*Jo*z);
set title sprintf("H_a/J_cd=%.2f",b)
#set title sprintf("B_a=%.1f [mT]",b*1000) 

#filein index n using (($2*1000-w/2)/w):(($3*1000-(1*w/2))/w):($9/Jo) with image  notitle , \
#

plot \
		filein index n using (($2*1000-w/2)/w):(($3*1000-(2*w/2))/w):(sqrt($9**2+$10**2+$11**2)/Jo) with image  notitle , \
		"data_out/contours.txt" using (($1*1000-w/2)/w):(($2*1000-(2*w/2))/w) with lines ls 1 lc "gray-50" notitle


#jpeg ls 1 lc 10

