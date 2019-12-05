#!/bin/python

#input
#form:0-jpeg, 1-eps 
form=0;
#title:0-z distance, 1-Ba;
title=1
#Jx:0-no, 1-yes 
Jx=0;
#Jy:0-no, 1-yes 
Jy=1;
#Jz:0-no, 1-yes 
Jz=1;
#J:0-no, 1-yes 
J=0;
#Ax:0-no, 1-yes 
Ax=0;
#By:0-no, 1-yes 
Ay=0;
#Az:0-no, 1-yes 
Az=0;

filepar="par_plot.plt"
filein="output3Dx.txt"


########################################################################

load filepar

a= x/(ncx);

i=0;
c=step-1;
#c=1;

k=1000

h=1.3
h1=1.3
g=Bamax*k;

set xrange [0:y*k]
set yrange [0:z*k]
set size ratio 	-1
set xlabel "y [mm]" 
set ylabel "z [mm]" 
set xtics 2
set palette rgbformulae 30,31,32
lwidth=1
set samples 200
set style data lines

########################################################################

if(Jx==1){
if(ncx!=1){
do for [n=0:c]{
if(form==0){
set terminal jpeg
set output sprintf("dJx%.0f.jpeg",n)}
else{
set terminal postscript eps enhanced color lw 1 "Arial" 26
set output sprintf("dJx%.0f.eps",n)}

set cblabel "J_x/J_c"
set key off        
set cbrange [-h:h]

if(title==0){
b = a*(i+1) - (a/2);
set title sprintf("z=%.2f [mm]",b*1000)}
else{
ts= (n+1)*dt;
b= Bamax * sin (2 * pi * f * ts);
set title sprintf("B_a=%.1f [mT]",b*1000)} 

plot filein index n using ($3*k):($4*k):($9/Jc) with image  notitle
}}}

########################################################################

if(Jy==1){
do for [n=0:c]{
if(form==0){
set terminal jpeg
set output sprintf("dJy%.0f.jpeg",n)}
else{
set terminal postscript eps enhanced color lw 1 "Arial" 26
set output sprintf("dJy%.0f.eps",n)}

set cblabel "J_y/J_c"
set key off       
set cbrange [-h:h]        

if(title==0){
b = a*(i+1) - (a/2);
set title sprintf("z=%.2f [mm]",b*1000)}
else{
ts= (n+1)*dt;
b= Bamax * sin (2 * pi * f * ts);
set title sprintf("B_a=%.1f [mT]",b*1000)} 

plot filein index n using ($3*k):($4*k):($10/Jc) with image  notitle
}}

########################################################################

if(Jz==1){
do for [n=0:c]{
if(form==0){
set terminal jpeg
set output sprintf("dJz%.0f.jpeg",n)}
else{
set terminal postscript eps enhanced color lw 1 "Arial" 26
set output sprintf("dJz%.0f.eps",n)}

set cblabel "J_z/J_c"
set key off           
set cbrange [-h1:h1]      

if(title==0){
b = a*(i+1) - (a/2);
set title sprintf("z=%.2f [mm]",b*1000)}
else{
ts= (n+1)*dt;
b= Bamax * sin (2 * pi * f * ts);
set title sprintf("B_a=%.1f [mT]",b*1000)}

plot filein index n using ($3*k):($4*k):($11/Jc) with image  notitle
}}

########################################################################

if(J==1){
do for [n=0:c]{
if(form==0){
set terminal jpeg
set output sprintf("J%.0f.jpeg",n)}
else{
set terminal postscript eps enhanced color lw 1 "Arial" 26
set output sprintf("J%.0f.eps",n)}

set cblabel "|J|/J_c"
set key off
set cbrange [-h:h]

if(title==0){
b = a*(i+1) - (a/2);
set title sprintf("z=%.2f [mm]",b*1000)}
else{
ts= (n+1)*dt;
b= Bamax * sin (2 * pi * f * ts);
set title sprintf("B_a=%.1f [mT]",b*1000)}

plot filein index n using ($3*k):($4*k):($12/Jc) with image  notitle
}} 

########################################################################

if(Ax==1){
do for [n=0:c]{
if(form==0){
set terminal jpeg
set output sprintf("Ax%.0f.jpeg",n)}
else{
set terminal postscript eps enhanced color lw 1 "Arial" 26
set output sprintf("Ax%.0f.eps",n)}

set cblabel "A_x[mT/m]"
set key off 
set cbrange [-g:g]            
set autoscale

if(title==0){
b = a*(i+1) - (a/2);
set title sprintf("z=%.2f [mm]",b*1000)}
else{
ts= (n+1)*dt;
b= Bamax * sin (2 * pi * f * ts);
set title sprintf("B_a=%.1f [mT]",b*1000)}

plot filein index n using ($3*k):($4*k):($16*k) with image  notitle
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

set cblabel "A_y[mT/m]"
set key off
set cbrange [-g:g]             
set autoscale

if(title==0){
b = a*(i+1) - (a/2);
set title sprintf("z=%.2f [mm]",b*1000)}
else{
ts= (n+1)*dt;
b= Bamax * sin (2 * pi * f * ts);
set title sprintf("B_a=%.1f [mT]",b*1000)}

plot filein index n using ($3*k):($4*k):($17*k) with image  notitle 
}}

########################################################################

if(Az==1){
do for [n=0:c]{
if(form==0){
set terminal jpeg
set output sprintf("Az%.0f.jpeg",n)}
else{
set terminal postscript eps enhanced color lw 1 "Arial" 26
set output sprintf("Az%.0f.eps",n)}

set cblabel "A_z[mT/m]" 
set key off       
set cbrange [-g:g] 
set autoscale

if(title==0){
b = a*(i+1) - (a/2);
set title sprintf("z=%.2f [mm]",b*1000)}
else{
ts= (n+1)*dt;
b= Bamax * sin (2 * pi * f * ts);
set title sprintf("B_a=%.1f [mT]",b*1000)}

plot filein index n using ($3*k):($4*k):($18*k) with image  notitle
}}

########################################################################

 
