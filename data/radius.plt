#!/bin/python

#input
#form:0-jpeg, 1-eps 
form=0;
#title:0-z distance, 1-Ba;
title=1
#smooth:0-no,1-yes
smooth=0;
#jc:0-no, 1-yes 
jc=0;
#m:0-no, 1-yes 
m=0;
#B_other:0-no, 1-yes 
Bo=0;
#T vector:0-no, 1-yes 
Tv=1;
#Dissipation factor:0-no, 1-yes 
U=0;
#plane from cube:0-no, 1-yes 
S=0;
#scale:0-h;1-autoscale
sc=0;

filepar="par_plot.plt"
if(smooth==0){
filein="output3Dz.txt"}
else{
filein="output3Dz_smooth.txt"}

########################################################################

load filepar

R=x*1000/2;
a= z/(ncz);

i=12;
n=173+i;

c=step-1;
#c=1;
n=0;

h=1.50;
g=Bamax*1000;
j=h/10;
U_scale=5e7;


set xrange [0:x*1000]
set yrange [0:y*1000]
set xtics 2
set xlabel "x [mm]" 
set ylabel "y [mm]" 
set size ratio 	-1
########################################################################

Y(x)=(R**2-(x-R)**2)**(1/2)

#set parametric

# Parametric functions for a circle
fx(t) = R*cos(t)
fy(t) = R*sin(t)

set parametric
set trange [0:2*pi]
# Parametric functions for a circle
fx(t) = R+R*cos(t)
fy(t) = R+R*sin(t)

print "R=",R

########################################################################

if(form==0){
set terminal jpeg
set output sprintf("2Dz_J%.0f.jpeg",n)}
else{
set terminal postscript eps enhanced color lw 1 "Arial" 26
set output sprintf("2Dz_J%.0f.eps",n)}

set cblabel "|J|/J_c"
set key off
set style data lines 
set palette rgbformulae 30,31,32
lwidth=1
set samples 1000
set style data lines   
if(sc==0){
set cbrange [-h	:h]}
else{
set autoscale}

set trange [0:2*pi]

if(title==0){
b = a*(i+1) - (a/2);
set title sprintf("z=%.2f [mm]",b*1000)}
else{
ts= (n+1)*dt;
b= Bamax * sin (2 * pi * f * ts);
set title sprintf("B_a=%.2f [mT]",b*1000)}

if(sc==0){
plot filein index n using ($2*1000):($3*1000):($9/Jo) with image  notitle	,	\
		 fx(t),fy(t) lt 1 lc 7 notitle

}
else{
plot filein index n using ($2*1000):($3*1000):($9) with image  notitle}


