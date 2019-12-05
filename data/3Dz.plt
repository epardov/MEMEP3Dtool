#!/bin/python

filepar="par_plot.plt"
filein="output3Dz.txt"

#Jx:0-no, 1-yes
Jx=1;
#Jy:0-no, 1-yes
Jy=1;
#Jz:0-no, 1-yes
Jz=1;
#J:0-no, 1-yes
J=0;
#Jcpar:0-no, 1-yes
Jcpar=0;
#Jc:0-no, 1-yes
Jc0=0;
#nc:0-no, 1-yes
nc=0;
#U:0-no, 1-yes
U=0;
#Bx:0-no, 1-yes
Bx=0;
#By:0-no, 1-yes
By=0;
#Bz:0-no, 1-yes
Bz=0;
#Ax:0-no, 1-yes
Ax=0;
#Ay:0-no, 1-yes
Ay=0;
#Az:0-no, 1-yes
Az=0;
#Tx:0-no, 1-yes
Tx=0;
#Ty:0-no, 1-yes
Ty=0;
#Tz:0-no, 1-yes
Tz=0;

########################################################################

load filepar
t= (step+1)*dt;
Ba= Bamax * sin (2 * pi * f * t)*1e3;
a= z/ncz;


c=step-1;
#c=2;

time_shift=0;

set xrange [0:x*1000]
set yrange [0:y*1000]
set style data lines 
set palette rgbformulae 30,31,32
lwidth=1
set samples 200
set style data lines 

g=Bamax*1000;
h=1.3;
hjc=1.1;
hnc=1.1;
Uc=1e5
B=0.3;
j=0.10;
A=1e-3;
T=1e5

########################################################################

if(Jx==1){
do for [n=0:c]{
set size ratio -1
set terminal pdf color enhanced font 'Arial, 16'  
set output sprintf("3Dz_dJx%.0f.pdf",n)

set xlabel "x [mm]"
set ylabel "y [mm]"
set cblabel "J_x/J_c"
set key off
set cbrange [-h:h]          

do for [i=0:ncz-1]{
b = a*(i+1) - (a/2);
set title sprintf("z=%.2f [mm]",b*1000)
plot filein index (i+time_shift) using ($2*1000):($3*1000):($9/Jc) with image  notitle
}
time_shift = (n+1)*ncz;
}}  

########################################################################
time_shift=0;
########################################################################

if(Jy==1){
do for [n=0:c]{
set size ratio -1
set terminal pdf color enhanced font 'Arial, 16'  
set output sprintf("3Dz_dJy%.0f.pdf",n)

set xlabel "x [mm]"
set ylabel "y [mm]"
set cblabel "J_y/J_c"
set key off   
set cbrange [-h:h]         

do for [i=0:ncz-1]{
b = a*(i+1) - (a/2);
set title sprintf("z=%.2f [mm]",b*1000)
plot filein index (i+time_shift) using ($2*1000):($3*1000):($10/Jc) with image  notitle
}  
time_shift = (n+1)*ncz;
}}

########################################################################
time_shift=0;
########################################################################

if(Jz==1){
do for [n=0:c]{
set size ratio -1
set terminal pdf color enhanced font 'Arial, 16'  
set output sprintf("3Dz_dJz%.0f.pdf",n)

set xlabel "x [mm]"
set ylabel "y [mm]"
set cblabel "J_z/J_c"
#unset colorbox
set key off     
set cbrange [-j:j]            

do for [i=0:ncz-1]{
b = a*(i+1) - (a/2);
set title sprintf("z=%.2f [mm]",b*1000)
plot filein index (i+time_shift) using ($2*1000):($3*1000):($11/Jc) with image  notitle
}  
time_shift = (n+1)*ncz;
}}

########################################################################
time_shift=0;
########################################################################

if(J==1){
do for [n=0:c]{
set size ratio -1
set terminal pdf color enhanced font 'Arial, 16'  
set output sprintf("3Dz_J%.0f.pdf",n)

set xlabel "x [mm]"
set ylabel "y [mm]"
set cblabel "|J|/J_c"
set key off   
set cbrange [0:h]           

do for [i=0:ncz-1]{
b = a*(i+1) - (a/2);
set title sprintf("z=%.2f [mm]",b*1000)
plot filein index (i+time_shift) using ($2*1000):($3*1000):($12/Jc) with image  notitle
}  
time_shift = (n+1)*ncz;
}}

########################################################################
time_shift=0;
########################################################################

if(Jcpar==1){
do for [n=0:c]{
set size ratio -1
set terminal pdf color enhanced font 'Arial, 16'  
set output sprintf("3Dz_Jcpar%.0f.pdf",n)

set xlabel "x [mm]"
set ylabel "y [mm]"
set cblabel "J_{c,par}/J_{c,per}"
set key off   
set cbrange [0:4]           

do for [i=0:ncz-1]{
b = a*(i+1) - (a/2);
set title sprintf("z=%.2f [mm]",b*1000)
plot filein index (i+time_shift) using ($2*1000):($3*1000):($8/Jc) with image  notitle
}  
time_shift = (n+1)*ncz;
}}

########################################################################
time_shift=0;
########################################################################

if(Jc0==1){
do for [n=0:c]{
set size ratio -1
set terminal pdf color enhanced font 'Arial, 16'  
set output sprintf("3Dz_Jc%.0f.pdf",n)

set xlabel "x [mm]"
set ylabel "y [mm]"
set cblabel "J_c/J_{c0}"
set key off   
set cbrange [0:hjc]           

do for [i=0:ncz-1]{
b = a*(i+1) - (a/2);
set title sprintf("z=%.2f [mm]",b*1000)
plot filein index (i+time_shift) using ($2*1000):($3*1000):($13/Jc) with image  notitle
}  
time_shift = (n+1)*ncz;
}}

########################################################################
time_shift=0;
########################################################################

if(nc==1){
do for [n=0:c]{
set size ratio -1
set terminal pdf color enhanced font 'Arial, 16'  
set output sprintf("3Dz_nc%.0f.pdf",n)

set xlabel "x [mm]"
set ylabel "y [mm]"
set cblabel "n/n_0"
set key off   
set cbrange [0:hnc]           

do for [i=0:ncz-1]{
b = a*(i+1) - (a/2);
set title sprintf("z=%.2f [mm]",b*1000)
plot filein index (i+time_shift) using ($2*1000):($3*1000):($15/N) with image  notitle
}  
time_shift = (n+1)*ncz;
}}

########################################################################
time_shift=0;
########################################################################

if(U==1){
do for [n=0:c]{
set size ratio -1
set terminal pdf color enhanced font 'Arial, 16'  
set output sprintf("3Dz_U%.0f.pdf",n)

set xlabel "x [mm]"
set ylabel "y [mm]"
set cblabel "U"
set key off   
set cbrange [0:Uc]           

do for [i=0:ncz-1]{
b = a*(i+1) - (a/2);
set title sprintf("z=%.2f [mm]",b*1000)
plot filein index (i+time_shift) using ($2*1000):($3*1000):($14) with image  notitle
}  
time_shift = (n+1)*ncz;
}}

########################################################################
time_shift=0;
########################################################################

if(Bx==1){
do for [n=0:c]{
set size ratio -1
set terminal pdf color enhanced font 'Arial, 16'  
set output sprintf("3Dz_Bx%.0f.pdf",n)

set xlabel "x [mm]"
set ylabel "y [mm]"
set cblabel "Bx/B_a,max"
set key off 
set cbrange [-B:B]           

do for [i=0:ncz-1]{
b = a*(i+1) - (a/2);
set title sprintf("z=%.2f [mm]",b*1000)
plot filein index (i+time_shift) using ($2*1000):($3*1000):($22/Bamax) with image  notitle
}  
time_shift = (n+1)*ncz;
}}

########################################################################
time_shift=0;
########################################################################

if(By==1){
do for [n=0:c]{
set size ratio -1
set terminal pdf color enhanced font 'Arial, 16'  
set output sprintf("3Dz_By%.0f.pdf",n)

set xlabel "x [mm]"
set ylabel "y [mm]"
set cblabel "By/B_a,max"
set key off 
set cbrange [-B:B]           

do for [i=0:ncz-1]{
b = a*(i+1) - (a/2);
set title sprintf("z=%.2f [mm]",b*1000)
plot filein index (i+time_shift) using ($2*1000):($3*1000):($23/Bamax) with image  notitle
}  
time_shift = (n+1)*ncz;
}}

########################################################################
time_shift=0;
########################################################################

if(Bz==1){
do for [n=0:c]{
set size ratio -1
set terminal pdf color enhanced font 'Arial, 16'  
set output sprintf("3Dz_Bz%.0f.pdf",n)

set xlabel "x [mm]"
set ylabel "y [mm]"
set cblabel "Bz/B_a,max"
set key off 
set cbrange [-B:B]           

do for [i=0:ncz-1]{
b = a*(i+1) - (a/2);
set title sprintf("z=%.2f [mm]",b*1000)
plot filein index (i+time_shift) using ($2*1000):($3*1000):($24/Bamax) with image  notitle
}  
time_shift = (n+1)*ncz;
}}

########################################################################
time_shift=0;
########################################################################

if(Ax==1){
do for [n=0:c]{
set size ratio -1
set terminal pdf color enhanced font 'Arial, 16'  
set output sprintf("3Dz_Ax%.0f.pdf",n)

set xlabel "x [mm]"
set ylabel "y [mm]"
set cblabel "Ax"
set key off 
set cbrange [-A:A]           

do for [i=0:ncz-1]{
b = a*(i+1) - (a/2);
set title sprintf("z=%.2f [mm]",b*1000)
plot filein index (i+time_shift) using ($2*1000):($3*1000):($16) with image  notitle
}  
time_shift = (n+1)*ncz;
}}

########################################################################
time_shift=0;
########################################################################

if(Ay==1){
do for [n=0:c]{
set size ratio -1
set terminal pdf color enhanced font 'Arial, 16'  
set output sprintf("3Dz_Ay%.0f.pdf",n)

set xlabel "x [mm]"
set ylabel "y [mm]"
set cblabel "Ay"
set key off 
set cbrange [-A:A]           

do for [i=0:ncz-1]{
b = a*(i+1) - (a/2);
set title sprintf("z=%.2f [mm]",b*1000)
plot filein index (i+time_shift) using ($2*1000):($3*1000):($17) with image  notitle
}  
time_shift = (n+1)*ncz;
}}

########################################################################
time_shift=0;
########################################################################

if(Az==1){
do for [n=0:c]{
set size ratio -1
set terminal pdf color enhanced font 'Arial, 16'  
set output sprintf("3Dz_Az%.0f.pdf",n)

set xlabel "x [mm]"
set ylabel "y [mm]"
set cblabel "Az"
set key off 
set cbrange [-A:A]           

do for [i=0:ncz-1]{
b = a*(i+1) - (a/2);
set title sprintf("z=%.2f [mm]",b*1000)
plot filein index (i+time_shift) using ($2*1000):($3*1000):($18) with image  notitle
}  
time_shift = (n+1)*ncz;
}}

########################################################################
time_shift=0;
########################################################################

if(Tx==1){
do for [n=0:c]{
set size ratio -1
set terminal pdf color enhanced font 'Arial, 16'  
set output sprintf("3Dz_Tx%.0f.pdf",n)

set xlabel "x [mm]"
set ylabel "y [mm]"
set cblabel "Tx"
set key off 
set cbrange [-T:T]           

do for [i=0:ncz-1]{
b = a*(i+1) - (a/2);
set title sprintf("z=%.2f [mm]",b*1000)
plot filein index (i+time_shift) using ($2*1000):($3*1000):($25) with image  notitle
}  
time_shift = (n+1)*ncz;
}}

########################################################################
time_shift=0;
########################################################################

if(Ty==1){
do for [n=0:c]{
set size ratio -1
set terminal pdf color enhanced font 'Arial, 16'  
set output sprintf("3Dz_Ty%.0f.pdf",n)

set xlabel "x [mm]"
set ylabel "y [mm]"
set cblabel "Ty"
set key off 
set cbrange [-T:T]           

do for [i=0:ncz-1]{
b = a*(i+1) - (a/2);
set title sprintf("z=%.2f [mm]",b*1000)
plot filein index (i+time_shift) using ($2*1000):($3*1000):($26) with image  notitle
}  
time_shift = (n+1)*ncz;
}}

########################################################################
time_shift=0;
########################################################################

if(Tz==1){
do for [n=0:c]{
set size ratio -1
set terminal pdf color enhanced font 'Arial, 16'  
set output sprintf("3Dz_Tz%.0f.pdf",n)

set xlabel "x [mm]"
set ylabel "y [mm]"
set cblabel "Tz"
set key off 
set cbrange [-T:T]           

do for [i=0:ncz-1]{
b = a*(i+1) - (a/2);
set title sprintf("z=%.2f [mm]",b*1000)
plot filein index (i+time_shift) using ($2*1000):($3*1000):($27) with image  notitle
}  
time_shift = (n+1)*ncz;
}}

########################################################################
time_shift=0;
########################################################################

