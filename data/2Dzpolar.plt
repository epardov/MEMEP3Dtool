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
#Ar:0-no, 1-yes 
Ar=0;
#Af:0-no, 1-yes
Af=0;
#us:0-no, 1-yes
us=0;
#at:0-no, 1-yes
at=0;

#Jx:0-no, 1-yes 
Jx=1;
#Jy:0-no, 1-yes 
Jy=1;
#J:0-no, 1-yes
J=1;
#Ax:0-no, 1-yes
Ax=0;
#Ay:0-no, 1-yes
Ay=0;
#Az:0-no, 1-yes
Az=0;
#2D:0-no, 1-yes
D=1;

filepar="par_plot.plt"
filein="output3Dz.txt"

########################################################################

load filepar

a= z/(ncz);

i=12;
n=173+i;

c=step-1;
#c=1;
#Jc=1

k = 1000;
hr= 0.1
h = 1.3;
a=10e-5
g = Bamax*k;
j = h/10;
U_scale = 5e8;

if(rel==3){Jc=Jcpe}

if(sys==0){set xrange [0:X*k]}
else{set xrange [(R1-dR)*k:R1*k]}
#set yrange [0:y*1000]
set xtics 5
set xlabel "r[mm]" 
set ylabel "fi [degree]" 
set size ratio 	1
set palette rgbformulae 30,31,32
lwidth=1
#set samples 200
#set style data lines

########################################################################

n=0;
if(Jr==1){
do for [n=0:c]{
if(form==0){
set terminal jpeg
set output sprintf("Jr%.0f.jpeg",n)}
else{
set terminal postscript eps enhanced color lw 1 "Arial" 26
set output sprintf("Jr%.0f.eps",n)}
set cblabel "J_r/J_c"

set key off          
set cbrange [-hr:hr]

ts= (n+1)*dt;
b= Bamax * sin (2 * pi * f * ts);
set title sprintf("B_a=%.2f [mT]",b*k)

plot filein index n using ($5*k):($6*180/pi):($9/Jc) with image  notitle
}}

########################################################################

n=0;
if(Jf==1){
do for [n=0:c]{
if(form==0){
set terminal jpeg
set output sprintf("Jf%.0f.jpeg",n)}
else{
set terminal postscript eps enhanced color lw 1 "Arial" 26
set output sprintf("Jf%.0f.eps",n)}
set cblabel "J_f/J_c"

set key off          
set cbrange [-h:h]


ts= (n+1)*dt;
b= Bamax * sin (2 * pi * f * ts);
set title sprintf("B_a=%.2f [mT]",b*k)

plot filein index n using ($5*k):($6*180/pi):($10/Jc) with image  notitle
}}

######################################################

n=0;
if(Ar==1){
do for [n=0:c]{
if(form==0){
set terminal jpeg
set output sprintf("Ar%.0f.jpeg",n)}
else{
set terminal postscript eps enhanced color lw 1 "Arial" 26
set output sprintf("Ar%.0f.eps",n)}
set cblabel "A_r[]"


set size ratio 	1
set key off          
#set cbrange [-h:h]
set autoscale

ts= (n+1)*dt;
b= Bamax * sin (2 * pi * f * ts);
set title sprintf("B_a=%.2f [mT]",b*k)

plot filein index n using ($5*k):($6*180/pi):($16) with image  notitle
}}

######################################################

n=0;

if(Af==1){
do for [n=0:c]{
if(form==0){
set terminal jpeg
set output sprintf("Af%.0f.jpeg",n)}
else{
set terminal postscript eps enhanced color lw 1 "Arial" 26
set output sprintf("Af%.0f.eps",n)}
set cblabel "A_f[]"

set size ratio 	1

set key off          

#set cbrange [-h:h]
set autoscale

ts= (n+1)*dt;
b= Bamax * sin (2 * pi * f * ts);
set title sprintf("B_a=%.2f [mT]",b*k)

plot filein index n using ($5*k):($6*180/pi):($17) with image  notitle
}}

######################################################

n=0;

if(at==1){
if(form==0){
set terminal jpeg
set output sprintf("at%.0f.jpeg",n)}
else{
set terminal postscript eps enhanced color lw 1 "Arial" 26
set output sprintf("at%.0f.eps",n)}
set cblabel "m-2"

set size ratio 	1

set key off          

#set cbrange [-h:h]
set autoscale

ts= (n+1)*dt;
b= Bamax * sin (2 * pi * f * ts);
set title sprintf("B_a=%.2f [mT]",b*k)

plot "output_Jy.txt" index n using ($2):($3):($11) with image  notitle
}

######################################################

if(us==1){
if(form==0){
set terminal jpeg
set output sprintf("us%.0f.jpeg",n)}
else{
set terminal postscript eps enhanced color lw 1 "Arial" 26
set output sprintf("us%.0f.eps",n)}
set cblabel "m-2"

set size ratio 	1

set key off          

#set cbrange [-h:h]
set autoscale

ts= (n+1)*dt;
b= Bamax * sin (2 * pi * f * ts);
set title sprintf("B_a=%.2f [mT]",b*k)

plot "output_Jy.txt" index n using ($2):($3):($12) with image  notitle
}

######################################################

filepar="par_plot.plt"
k=1000;
set xlabel "x [mm]" 
set ylabel "y [mm]" 
filein="data_polar.txt"

set xrange [sxl*k:sxr*k]
set yrange [syb*k:syt*k]
n=0;

set xtics 25
set xtics format "%.1f"
set ytics 25

######################################################

if(Jx==1){
do for [n=0:c]{
if(form==0){
set terminal jpeg
set output sprintf("Jx%.0f.jpeg",n)}
else{
set terminal postscript eps enhanced color lw 1 "Arial" 26
set output sprintf("Jx%.0f.eps",n)}
set cblabel "J_x/J_c"

set key off 
        
set cbrange [-h:h]
#set autoscale

ts= (n+1)*dt;
b= Bamax * sin (2 * pi * f * ts);
set title sprintf("B_a=%.2f [mT]",b*k)

plot filein index n using ($2*k):($3*k):($5/Jc) with image  notitle
}}

########################################################################

if(Jy==1){
do for [n=0:c]{
if(form==0){
set terminal jpeg
set output sprintf("Jy%.0f.jpeg",n)}
else{
set terminal postscript eps enhanced color lw 1 "Arial" 26
set output sprintf("Jy%.0f.eps",n)}
set cblabel "J_y/J_c"

set key off          
set cbrange [-h:h]
#set autoscale

ts= (n+1)*dt;
b= Bamax * sin (2 * pi * f * ts);
set title sprintf("B_a=%.2f [mT]",b*k)

plot filein index n using ($2*k):($3*k):($6/Jc) with image  notitle
}}

######################################################

if(J==1){
do for [n=0:c]{
if(form==0){
set terminal jpeg
set output sprintf("J%.0f.jpeg",n)}
else{
set terminal postscript eps enhanced color lw 1 "Arial" 26
set output sprintf("J%.0f.eps",n)}
set cblabel "J/J_c"

set key off          
set cbrange [0:h]
#set autoscale

ts= (n+1)*dt;
b= Bamax * sin (2 * pi * f * ts);
set title sprintf("B_a=%.2f [mT]",b*k)

plot filein index n using ($2*k):($3*k):(sqrt($5**2+$6**2)/Jc) with image  notitle
}}

########################################################################
time_shift=0;
########################################################################

if(ncz!=1){
if(J==1){
do for [n=0:c]{
set size ratio -1
set terminal pdf color enhanced font 'Arial, 16'  
set output sprintf("J%.0f.pdf",n)

set xlabel "x [mm]"
set ylabel "y [mm]"
set cblabel "|J|/J_c"
set key off   
set cbrange [0:h]           

do for [i=0:ncz-1]{
b = a*(i+1) - (a/2);
set title sprintf("z=%.2f [mm]",b*1000)
plot filein index (i+time_shift) using ($2*k):($3*k):(sqrt($5**2+$6**2+$7**2)/Jc) with image  notitle
}  
time_shift = (n+1)*ncz;
}}}

########################################################################

if(Ax==1){
do for [n=0:c]{
if(form==0){
set terminal jpeg
set output sprintf("Ax%.0f.jpeg",n)}
else{
set terminal postscript eps enhanced color lw 1 "Arial" 26
set output sprintf("Ax%.0f.eps",n)}
set cblabel "A_x"

set key off 
set cbrange [-a:a]         
set autoscale

ts= (n+1)*dt;
b= Bamax * sin (2 * pi * f * ts);
set title sprintf("B_a=%.2f [mT]",b*k)

plot filein index n using ($2*k):($3*k):($8) with image  notitle
}}

########################################################################

if(Ay==1){
do for [n=0:c]{
if(form==0){
set terminal jpeg
set output sprintf("Ay%.0f.jpeg",n)}
else{
set terminal postscript eps enhanced color lw 1 "Arial" 26
set output sprintf("Ay%.0f.eps",n)}
set cblabel "A_y"

set key off          
set cbrange [-a:a]         
set autoscale

ts= (n+1)*dt;
b= Bamax * sin (2 * pi * f * ts);
set title sprintf("B_a=%.2f [mT]",b*k)

plot filein index n using ($2*k):($3*k):($9) with image  notitle
}}

