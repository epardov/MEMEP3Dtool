#!/bin/python

#input
#form:0-jpeg, 1-eps 
form=0;
#title:0-z distance, 1-Ba;
title=1;
#Jf:0-no, 1-yes 
Jf=1;
#Jz:0-no, 1-yes 
Jz=1;
#Br:0-no, 1-yes
Br=0;
#Bf:0-no, 1-yes
Bf=0;
#Bz:0-no, 1-yes
Bz=0;
#Ar:0-no, 1-yes
Ar=0;
#Af:0-no, 1-yes
Af=0;
#Az:0-no, 1-yes 
Az=0;
#us:0-no, 1-yes
us=0;
#Tx:0-no, 1-yes
Tx=0;
#at:0-no, 1-yes
at=0;


filepar="par_plot.plt"
filein="output3Dx.txt"

########################################################################

load filepar

a= z/(ncz);

i=12;
n=173+i;

c=step-1;
#c=19;
#Jc=1

k  = 1000
h  = 1.5
h1 = 0.5
a=10.0e-1
b1= 3.0;
b2= 1.5
g = Bamax*k;
j = h/10;
U_scale = 5e8;

if(rel==3){Jc=Jcpe}
FI1=5
set xrange [0-FI1:FI*180/pi + FI1]
set xtics 720
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
if(Ismax!=0){set cbrange [0:h]}
else{set cbrange [-h:h]}


ts= (n+1)*dt;
b= Bamax * sin (2 * pi * f * ts);
set title sprintf("B_a=%.2f [mT]",b*k)

plot filein index n using ($6*180/pi):($7*k):($10/Jc) with image  notitle
}}

########################################################################

n=0;
if(Jz==1){
do for [n=0:c]{
if(form==0){
set terminal jpeg
set output sprintf("Jz%.0f.jpeg",n)}
else{
set terminal postscript eps enhanced color lw 1 "Arial" 26
set output sprintf("Jz%.0f.eps",n)}
set cblabel "J_z/J_c"

set key off          
set cbrange [-h1:h1]

ts= (n+1)*dt;
b= Bamax * sin (2 * pi * f * ts);
set title sprintf("B_a=%.2f [mT]",b*k)

plot filein index n using ($6*180/pi):($7*k):($11/Jc) with image  notitle
}}

########################################################################

n=0;

if(Br==1){
do for [n=0:c]{
if(form==0){
set terminal jpeg
set output sprintf("Br%.0f.jpeg",n)}
else{
set terminal postscript eps enhanced color lw 1 "Arial" 26
set output sprintf("Br%.0f.eps",n)}
set cblabel "B_r[mT]"

set key off          
set cbrange [-b1:b1]
#set autoscale

ts= (n+1)*dt;
b= Bamax * sin (2 * pi * f * ts);
set title sprintf("B_a=%.2f [mT]",b*k)

plot filein index n using ($6*180/pi):($7*k):($22*k) with image  notitle
}}

########################################################################

n=0;

if(Bf==1){
do for [n=0:c]{
if(form==0){
set terminal jpeg
set output sprintf("Bf%.0f.jpeg",n)}
else{
set terminal postscript eps enhanced color lw 1 "Arial" 26
set output sprintf("Bf%.0f.eps",n)}
set cblabel "B_f[mT]"

set key off          
set cbrange [-b1*0.10:b1*0.10]
#set autoscale

ts= (n+1)*dt;
b= Bamax * sin (2 * pi * f * ts);
set title sprintf("B_a=%.2f [mT]",b*k)

plot filein index n using ($6*180/pi):($7*k):($23*k) with image  notitle
}}

########################################################################

n=0;

if(Bz==1){
do for [n=0:c]{
if(form==0){
set terminal jpeg
set output sprintf("Bz%.0f.jpeg",n)}
else{
set terminal postscript eps enhanced color lw 1 "Arial" 26
set output sprintf("Bz%.0f.eps",n)}
set cblabel "B_z[mT]"

set key off          
set cbrange [-b2:b2]

#set autoscale

ts= (n+1)*dt;
b= Bamax * sin (2 * pi * f * ts);
set title sprintf("B_a=%.2f [mT]",b*k)

plot filein index n using ($6*180/pi):($7*k):($24*k) with image  notitle
}}

########################################################################

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

set key off          
set cbrange [-a:a]
#set autoscale

ts= (n+1)*dt;
b= Bamax * sin (2 * pi * f * ts);
set title sprintf("B_a=%.2f [mT]",b*k)

plot filein index n using ($6*180/pi):($7*k):($16) with image  notitle
}}

########################################################################

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

set key off          
set cbrange [-a:a]
#set autoscale

ts= (n+1)*dt;
b= Bamax * sin (2 * pi * f * ts);
set title sprintf("B_a=%.2f [mT]",b*k)

plot filein index n using ($6*180/pi):($7*k):($17) with image  notitle
}}

######################################################

n=0;

if(Az==1){
do for [n=0:c]{
if(form==0){
set terminal jpeg
set output sprintf("Az%.0f.jpeg",n)}
else{
set terminal postscript eps enhanced color lw 1 "Arial" 26
set output sprintf("Az%.0f.eps",n)}
set cblabel "A_z[]"

set key off          
set cbrange [-a:a]
set autoscale

ts= (n+1)*dt;
b= Bamax * sin (2 * pi * f * ts);
set title sprintf("B_a=%.2f [mT]",b*k)

plot filein index n using ($6*180/pi):($7*k):($18*k) with image  notitle
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

set key off          
#set cbrange [-h:h]
set autoscale

ts= (n+1)*dt;
b= Bamax * sin (2 * pi * f * ts);
set title sprintf("B_a=%.2f [mT]",b*k)

plot "output_Jy.txt" index n using ($2):($3):($11) with image  notitle
}

######################################################

if(Tx==1){
if(form==0){
set terminal jpeg
set output sprintf("Tx%.0f.jpeg",n)}
else{
set terminal postscript eps enhanced color lw 1 "Arial" 26
set output sprintf("Tx%.0f.eps",n)}
set cblabel "Tx"

set key off          
#set cbrange [-h:h]
set autoscale

ts= (n+1)*dt;
b= Bamax * sin (2 * pi * f * ts);
set title sprintf("B_a=%.2f [mT]",b*k)

plot filein index n using ($6*180/pi):($7*k):($25) with image  notitle
}

######################################################

if(us==1){
if(form==0){
set terminal jpeg
set output sprintf("us%.0f.jpeg",n)}
else{
set terminal postscript eps enhanced color lw 1 "Arial" 26
set output sprintf("us%.0f.eps",n)}
set cblabel "us"

set key off          
#set cbrange [-h:h]
set autoscale

ts= (n+1)*dt;
b= Bamax * sin (2 * pi * f * ts);
set title sprintf("B_a=%.2f [mT]",b*k)

plot "output_Jy.txt" index n using ($3):($4):($12) with image  notitle
}

