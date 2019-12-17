#!/bin/python

filepar="par_plot.plt"
filein="output3Dy.txt"

#Jx:0-no, 1-yes
Jx=0;
#Jy:0-no, 1-yes
Jy=1;
#Jz:0-no, 1-yes
Jz=0;
#J:0-no, 1-yes
J=0;

########################################################################

load filepar
t= (step+1)*dt;
Ba= Bamax * sin (2 * pi * f * t)*1e3;
a= y/(ncy);

c=step-1;
#c=0;

time_shift=0;

set xrange [0:x*1000]
set yrange [0:z*1000]
set size ratio 1
set style data lines 
set palette rgbformulae 30,31,32
lwidth=1
set samples 200
set style data lines 

g=Bamax*1000;
h=1.3;
j=0.5;

########################################################################

if(Jx==1){
do for [n=0:c]{
set terminal pdf color enhanced font 'Arial, 16'
set output sprintf("3Dy_dJx%.0f.pdf",n)

set xlabel "x [mm]"
set ylabel "z [mm]"
set cblabel "J_x/J_c"
set key off 
set cbrange [-h:h]     

do for [i=0:ncy-1]{
b = a*(i+1) - (a/2);
set title sprintf("y=%.2f [mm]",b*1000)
plot filein index (i+time_shift) using ($2*1000):($4*1000):($9/Jc) with image  notitle
}  
time_shift=(n+1)*ncy;
}}  

########################################################################
time_shift=0;
########################################################################
#n=519;
#time_shift=519*31;

if(Jy==1){
do for [n=0:c]{
set terminal pdf color enhanced font 'Arial, 16'
set output sprintf("3Dy_dJy%.0f.pdf",n)

set xlabel "x [mm]"
set ylabel "z [mm]"
set cblabel "J_y/J_c"
set key off    
set cbrange [-h:h]            

do for [i=0:ncy-1]{
b = a*(i+1) - (a/2);
set title sprintf("y=%.2f [mm]",b*1000)
plot filein index (i+time_shift) using ($2*1000):($4*1000):($10/Jc) with image  notitle
}  
time_shift=(n+1)*ncy;
}
}  

########################################################################
time_shift=0;
########################################################################

if(Jz==1){
do for [n=0:c]{
set terminal pdf color enhanced font 'Arial, 16'
set output sprintf("3Dy_dJz%.0f.pdf",n)

set xlabel "x [mm]"
set ylabel "z [mm]"
set cblabel "J_z/J_c"
set key off
set cbrange [-j:j]        

do for [i=0:ncy-1]{
b = a*(i+1) - (a/2);
set title sprintf("y=%.2f [mm]",b*1000)
plot filein index (i+time_shift) using ($2*1000):($4*1000):($11/Jc) with image  notitle
}  
time_shift=(n+1)*ncy;
}}   

########################################################################
time_shift=0;
########################################################################

if(J==1){
do for [n=0:c]{
set terminal pdf color enhanced font 'Arial, 16'
set output sprintf("3Dy_J%.0f.pdf",n)

set xlabel "x [mm]"
set ylabel "z [mm]"
set cblabel "|J|/J_c"
set key off  
set cbrange [0:h]           

do for [i=0:ncy-1]{
b = a*(i+1) - (a/2);
set title sprintf("y=%.2f [m]",b*1000)
plot filein index (i+time_shift) using ($2*1000):($4*1000):($12/Jc) with image  notitle
}  
time_shift=(n+1)*ncy;
}} 


