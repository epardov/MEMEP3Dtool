#!/bin/python

#input
#form:0-jpeg, 1-eps 
form=0;
#title:0-z distance, 1-Ba;
title=1
#Jx:0-no, 1-yes 
Jx=0;
#Jy:0-no, 1-yes 
Jy=0;
#Jz:0-no, 1-yes 
Jz=0;
#J:0-no, 1-yes 
J=1;
#T vector:0-no, 1-yes 
Tv=0;
#scale:0-h;1-autoscale
sc=0;

filepar="par_plot.plt"
filein="outputAV.txt"

########################################################################

load filepar

a= z/(ncz);

c=step-1;
#c=9;

h=1.50;
g=Bamax*1000;
j=h/10;

if(rel==3){Jo=Jcpe}


set xrange [0:x*1000]
set yrange [0:y*1000]
set xtics 2
set xlabel "x [mm]" 
set ylabel "y [mm]" 
set size ratio 	-1
########################################################################

if(Jx==1){
do for [n=0:c]{
if(form==0){
set terminal jpeg
set output sprintf("2Dz_dJx%.0f.jpeg",n)}
else{
set terminal postscript eps enhanced color lw 1 "Arial" 26
set output sprintf("2Dz_dJx%.0f.eps",n)}

set cblabel "J_x/J_c"
set key off
set style data lines 
set palette rgbformulae 30,31,32
lwidth=1
set samples 200
set style data lines          
if(sc==0){
set cbrange [-h	:h]}
else{
set autoscale}

if(title==0){
b = a*(i+1) - (a/2);
set title sprintf("z=%.2f [mm]",b*1000)}
else{
ts= (n+1)*dt;
b= Bamax * sin (2 * pi * f * ts);
set title sprintf("B_a=%.2f [mT]",b*1000)} 

if(sc==0){
plot filein index n using ($2*1000):($3*1000):($5/Jo) with image  notitle}
else{
plot filein index n using ($2*1000):($3*1000):($5) with image  notitle}
}}

########################################################################

if(Jy==1){
do for [n=0:c]{
if(form==0){
set terminal jpeg
set output sprintf("2Dz_dJy%.0f.jpeg",n)}
else{
set terminal postscript eps enhanced color lw 1 "Arial" 26
set output sprintf("2Dz_dJy%.0f.eps",n)}

set cblabel "J_y/J_c"
set key off
set style data lines 
set palette rgbformulae 30,31,32
lwidth=1
set samples 200
set style data lines        
if(sc==0){
set cbrange [-h	:h]}
else{
set autoscale}       

if(title==0){
b = a*(i+1) - (a/2);
set title sprintf("z=%.2f [mm]",b*1000)}
else{
ts= (n+1)*dt;
b= Bamax * sin (2 * pi * f * ts);
set title sprintf("B_a=%.2f [mT]",b*1000)} 

if(sc==0){
plot filein index n using ($2*1000):($3*1000):($6/Jo) with image  notitle}
else{
plot filein index n using ($2*1000):($3*1000):($6) with image  notitle}
}}

########################################################################


if(Jz==1){
do for [n=0:c]{
if(form==0){
set terminal jpeg
set output sprintf("2Dz_dJz%.0f.jpeg",n)}
else{
set terminal postscript eps enhanced color lw 1 "Arial" 26
set output sprintf("2Dz_dJz%.0f.eps",n)}

set cblabel "J_z/J_c"
set key off
set style data lines 
set palette rgbformulae 30,31,32
lwidth=1
set samples 200
set style data lines            
if(sc==0){
set cbrange [-h	:h]}
else{
set autoscale}      

if(title==0){
b = a*(i+1) - (a/2);
set title sprintf("z=%.2f [mm]",b*1000)}
else{
ts= (n+1)*dt;
b= Bamax * sin (2 * pi * f * ts);
set title sprintf("B_a=%.2f [mT]",b*1000)}

if(sc==0){
plot filein index n using ($2*1000):($3*1000):($7/Jo) with image  notitle}
else{
plot filein index n using ($2*1000):($3*1000):($7) with image  notitle}
}}

########################################################################

if(J==1){
do for [n=0:c]{
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
set samples 200
set style data lines   
if(sc==0){
set cbrange [0:h]}
else{
set autoscale}

if(title==0){
b = a*(i+1) - (a/2);
set title sprintf("z=%.2f [mm]",b*1000)}
else{
ts= (n+1)*dt;
b= Bamax * sin (2 * pi * f * ts);
set title sprintf("B_a=%.2f [mT]",b*1000)}

if(sc==0){
plot filein index n using ($2*1000):($3*1000):(sqrt($5**2+$6**2+$7**2)/Jo) with image  notitle}
else{
plot filein index n using ($2*1000):($3*1000):(sqrt($5**2+$6**2+$7**2)) with image  notitle}
}} 

########################################################################

if(Bx==1){
do for [n=0:c]{
if(form==0){
set terminal jpeg
set output sprintf("2Dz_Bx%.0f.jpeg",n)}
else{
set terminal postscript eps enhanced color lw 1 "Arial" 26
set output sprintf("2Dz_Bx%.0f.eps",n)}

set cblabel "B_x[mT]"
set key off
set style data lines 
set palette rgbformulae 30,31,32
lwidth=1
set samples 200
set style data lines   
set cbrange [-g:g]            

if(title==0){
b = a*(i+1) - (a/2);
set title sprintf("z=%.2f [mm]",b*1000)}
else{
ts= (n+1)*dt;
b= Bamax * sin (2 * pi * f * ts);
set title sprintf("B_a=%.2f [mT]",b*1000)}

plot filein index c using ($2*1000):($3*1000):($8*1000) with image  notitle
}}

########################################################################

if(By==1){
do for [n=0:c]{
if(form==0){
set terminal jpeg
set output sprintf("2Dz_By%.0f.jpeg",n)}
else{
set terminal postscript eps enhanced color lw 1 "Arial" 26
set output sprintf("2Dz_By%.0f.eps",n)}

set ylabel "y [mm]"
set cblabel "B_y[mT]"
set key off
set style data lines 
set palette rgbformulae 30,31,32
lwidth=1
set samples 200
set style data lines   
set cbrange [-g:g]             

if(title==0){
b = a*(i+1) - (a/2);
set title sprintf("z=%.2f [mm]",b*1000)}
else{
ts= (n+1)*dt;
b= Bamax * sin (2 * pi * f * ts);
set title sprintf("B_a=%.2f [mT]",b*1000)}

plot filein index c using ($2*1000):($3*1000):($9*1000) with image  notitle 
}}

########################################################################

if(Bz==1){
do for [n=0:c]{
if(form==0){
set terminal jpeg
set output sprintf("2Dz_Bz%.0f.jpeg",n)}
else{
set terminal postscript eps enhanced color lw 1 "Arial" 26
set output sprintf("2Dz_Bz%.0f.eps",n)}

set cblabel "B_z[mT]" 
set key off
set style data lines 
set palette rgbformulae 30,31,32
lwidth=1
set samples 200
set style data lines          
set cbrange [-g:g] 

if(title==0){
b = a*(i+1) - (a/2);
set title sprintf("z=%.2f [mm]",b*1000)}
else{
ts= (n+1)*dt;
b= Bamax * sin (2 * pi * f * ts);
set title sprintf("B_a=%.2f [mT]",b*1000)}

plot filein index n using ($2*1000):($3*1000):($10*1000) with image  notitle
}}

########################################################################

if(Tv==1){
do for [n=0:c]{
if(form==0){
set terminal jpeg
set output sprintf("2Dz_Tz%.0f.jpeg",n)}
else{
set terminal postscript eps enhanced color lw 1 "Arial" 26
set output sprintf("2Dz_Tz%.0f.eps",n)}

set cblabel "T [A/m2]"
set key off
set style data lines 
set palette rgbformulae 30,31,32
lwidth=1
set samples 200
set style data lines          
set cbrange [-2.5e8:2.5e8] 
#set autoscale

if(title==0){
b = a*(i+1) - (a/2);
set title sprintf("z=%.2f [mm]",b*1000)}
else{
ts= (n+1)*dt;
b= Bamax * sin (2 * pi * f * ts);
set title sprintf("B_a=%.2f [mT]",b*1000)} 

plot filein index n using ($2*1000):($3*1000):($24) with image  notitle
}}

########################################################################


