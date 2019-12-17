#!/bin/python

#input
#form:0-jpeg, 1-eps 
form=0;
#title:0-z distance, 1-Ba;
title=1;
#Jr:0-no, 1-yes 
Jr=1;
#Jf:0-no, 1-yes 
Jf=1;
#Jz:0-no, 1-yes 
Jz=1;
#U:0-no, 1-yes
U=1;

filepar="par_plot.plt"
filein="output3Dx.txt"

########################################################################

load filepar

a= z/(ncz);

i=12;
n=173+i;

c=step-1;
c=9;
#Jc=1

k  = 1000
h  = 1.3
hr = 0.05
hz = 0.1
a=10.0e-1
b1= 3.0;
b2= 1.5
g = Bamax*k;
j = h/10;
U_scale = 5e4;

if(rel==3){Jc=Jcpe}
FI1=5
set xrange [0-FI1:FI*180/pi + FI1]
set xtics 180
set yrange [0:Z*1000]
set ytics 1
set xlabel "fi [degree]" 
set ylabel "z[mm]" 
set size ratio 	1
set palette rgbformulae 30,31,32
lwidth=1
#set samples 200
#set style data lines

########################################################################
time_shift=0;
########################################################################

n=0;
if(Jr==1){
do for [n=0:c]{

set terminal pdf color enhanced font 'Arial, 16'
set output sprintf("3Dx_dJr%.0f.pdf",n)
set cblabel "J_r/J_c"

set key off          
if(Ismax!=0){set cbrange [-hr:hr]}
else{set cbrange [-hr:hr]}

#set autoscale 

ts= (n+1)*dt;
b= Bamax * sin (2 * pi * f * ts);
set title sprintf("B_a=%.2f [mT]",b*k)

do for [i=0:ncx-1]{
plot filein index (i+time_shift) using ($6*180/pi):($7*k):($9/Jc) with image  notitle
}
time_shift=(n+1)*ncx;
}}

########################################################################
time_shift=0;
########################################################################

n=0;
if(Jf==1){
do for [n=0:c]{

set terminal pdf color enhanced font 'Arial, 16'
set output sprintf("3Dx_dJf%.0f.pdf",n)
set cblabel "J_f/J_c"

set key off          
if(Ismax!=0){set cbrange [-h*0:h]}
else{set cbrange [-h:h]}

ts= (n+1)*dt;
b= Bamax * sin (2 * pi * f * ts);
set title sprintf("B_a=%.2f [mT]",b*k)

do for [i=0:ncx-1]{
plot filein index (i+time_shift) using ($6*180/pi):($7*k):($10/Jc) with image  notitle
}
time_shift=(n+1)*ncx;
}}

########################################################################
time_shift=0;
########################################################################

n=0;
if(Jz==1){
do for [n=0:c]{

set terminal pdf color enhanced font 'Arial, 16'
set output sprintf("3Dx_dJz%.0f.pdf",n)
set cblabel "J_z/J_c"

set key off          
set cbrange [-hz:hz]

ts= (n+1)*dt;
b= Bamax * sin (2 * pi * f * ts);
set title sprintf("B_a=%.2f [mT]",b*k)

do for [i=0:ncx-1]{
plot filein index (i+time_shift) using ($6*180/pi):($7*k):($11/Jc) with image  notitle
}
time_shift=(n+1)*ncx;
}}

########################################################################
time_shift=0;
########################################################################

n=0;
if(U==1){
do for [n=0:c]{

set terminal pdf color enhanced font 'Arial, 16'
set output sprintf("3Dx_U%.0f.pdf",n)
set cblabel "U[]"

set key off          
set cbrange [0:U_scale]
#set autoscale

ts= (n+1)*dt;
b= Bamax * sin (2 * pi * f * ts);
set title sprintf("B_a=%.2f [mT]",b*k)

do for [i=0:ncx-1]{
plot filein index (i+time_shift) using ($6*180/pi):($7*k):($14) with image  notitle
}
time_shift=(n+1)*ncx;
}}

########################################################################

