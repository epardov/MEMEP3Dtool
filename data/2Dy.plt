#!/bin/python

#input
#form:0-jpeg, 1-eps 
form=0;
#title:0-z distance, 1-Ba;
title=1
#smooth:0-no,1-yes
smooth=0;

filepar="par_plot.plt"
load filepar
filein="output3Dy.txt"

########################################################################



a= x/(ncx);

i=0;

c=step-1;
#c=1;

h=1.5;
g=Bamax*1000;

set xrange [0:x*1000]
set yrange [0:z*1000]
set xtics 1
set ytics 1
set xlabel "x [mm]" 
set ylabel "z [mm]" 
set size ratio 	-1
set palette rgbformulae 30,31,32
lwidth=1
set samples 200
set style data lines
########################################################################


do for [n=0:c]{
if(form==0){
set terminal jpeg
set output sprintf("3Dy_dJx%.0f.jpeg",n)}
else{
set terminal postscript eps enhanced color lw 1 "Arial" 26
set output sprintf("3Dy_dJx%.0f.eps",n)}

set cblabel "J_x/J_c"
set key off
set cbrange [-h	:h]

if(title==0){
b = a*(i+1) - (a/2);
set title sprintf("z=%.2f [mm]",b*1000)}
else{
ts= (n+1)*dt;
b= Bamax * sin (2 * pi * f * ts);
set title sprintf("B_a=%.1f [mT]",b*1000)} 

plot filein index n using ($2*1000):($4*1000):($9/Jo) with image  notitle
}

########################################################################
if(ncy!=1){
do for [n=0:c]{
if(form==0){
set terminal jpeg
set output sprintf("3Dy_dJy%.0f.jpeg",n)}
else{
set terminal postscript eps enhanced color lw 1 "Arial" 26
set output sprintf("3Dy_dJy%.0f.eps",n)}

set cblabel "J_y/J_c"
set key off
set cbrange [-h	:h]        

if(title==0){
b = a*(i+1) - (a/2);
set title sprintf("z=%.2f [mm]",b*1000)}
else{
ts= (n+1)*dt;
b= Bamax * sin (2 * pi * f * ts);
set title sprintf("B_a=%.1f [mT]",b*1000)} 

plot filein index n using ($2*1000):($4*1000):($10/Jo) with image  notitle
}}

########################################################################

do for [n=0:c]{
if(form==0){
set terminal jpeg
set output sprintf("3Dy_dJz%.0f.jpeg",n)}
else{
set terminal postscript eps enhanced color lw 1 "Arial" 26
set output sprintf("3Dy_dJz%.0f.eps",n)}

set cblabel "J_z/J_c"
set key off
set cbrange [-h	:h]      

if(title==0){
b = a*(i+1) - (a/2);
set title sprintf("z=%.2f [mm]",b*1000)}
else{
ts= (n+1)*dt;
b= Bamax * sin (2 * pi * f * ts);
set title sprintf("B_a=%.1f [mT]",b*1000)}

plot filein index n using ($2*1000):($4*1000):($11/Jo) with image  notitle
} 

########################################################################

do for [n=0:c]{
set size ratio 	1
if(form==0){
set terminal jpeg
set output sprintf("3Dy_J%.0f.jpeg",n)}
else{
set terminal postscript eps enhanced color lw 1 "Arial" 26
set output sprintf("3Dy_J%.0f.eps",n)}
set xtics 2

set cblabel "|J|/J_c"
set key off
set style data lines 
set palette rgbformulae 30,31,32
lwidth=1
set samples 200
set style data lines   
set cbrange [-h	:h]

if(title==0){
b = a*(i+1) - (a/2);
set title sprintf("z=%.2f [mm]",b*1000)}
else{
ts= (n+1)*dt;
b= Bamax * sin (2 * pi * f * ts);
set title sprintf("B_a=%.1f [mT]",b*1000)}

plot filein index n using ($2*1000):($3*1000):($9/Jo) with image  notitle
}

########################################################################
if(ncy!=1){
do for [n=0:c]{
set size ratio 	1
if(form==0){
set terminal jpeg
set output sprintf("3Dy_Bx%.0f.jpeg",n)}
else{
set terminal postscript eps enhanced color lw 1 "Arial" 26
set output sprintf("3Dy_Bx%.0f.eps",n)}
set xtics 2

set cblabel "B_x[mT]"
set key off
set style data lines 
set palette rgbformulae 30,31,32
lwidth=1
set samples 200
set style data lines   
set cbrange [-g	:g]            

if(title==0){
b = a*(i+1) - (a/2);
set title sprintf("z=%.2f [mm]",b*1000)}
else{
ts= (n+1)*dt;
b= Bamax * sin (2 * pi * f * ts);
set title sprintf("B_a=%.1f [mT]",b*1000)}

plot filein index n using ($2*1000):($3*1000):($19*1000) with image  notitle
}}

########################################################################


do for [n=0:c]{
set size ratio 	1
if(form==0){
set terminal jpeg
set output sprintf("3Dy_By%.0f.jpeg",n)}
else{
set terminal postscript eps enhanced color lw 1 "Arial" 26
set output sprintf("3Dy_By%.0f.eps",n)}
set xtics 2

set cblabel "B_y[mT]"
set key off
set style data lines 
set palette rgbformulae 30,31,32
lwidth=1
set samples 200
set style data lines   
set cbrange [-g	:g]             

if(title==0){
b = a*(i+1) - (a/2);
set title sprintf("z=%.2f [mm]",b*1000)}
else{
ts= (n+1)*dt;
b= Bamax * sin (2 * pi * f * ts);
set title sprintf("B_a=%.1f [mT]",b*1000)}

plot filein index n using ($2*1000):($3*1000):($20*1000) with image  notitle 
}

########################################################################

if(ncy!=1){
do for [n=0:c]{
set size ratio 	1
if(form==0){
set terminal jpeg
set output sprintf("3Dy_Bz%.0f.jpeg",n)}
else{
set terminal postscript eps enhanced color lw 1 "Arial" 26
set output sprintf("3Dy_Bz%.0f.eps",n)}
set xtics 2

set cblabel "B_z[mT]" 
set key off
set style data lines 
set palette rgbformulae 30,31,32
lwidth=1
set samples 200
set style data lines          
set cbrange [-g	:g] 

if(title==0){
b = a*(i+1) - (a/2);
set title sprintf("z=%.2f [mm]",b*1000)}
else{
ts= (n+1)*dt;
b= Bamax * sin (2 * pi * f * ts);
set title sprintf("B_a=%.1f [mT]",b*1000)}

plot filein index n using ($2*1000):($3*1000):($21*1000) with image  notitle
}}

########################################################################



