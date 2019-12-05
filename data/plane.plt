#!/bin/python

#input
#form:0-jpeg, 1-eps 
form=0;
#title:0-z distance, 1-Ba;
title=1
#BJx:0-no, 1-yes
BJx=0;
#BJy:0-no, 1-yes
BJy=0;
#BJz:0-no, 1-yes
BJz=0;
#Bx:0-no, 1-yes
Bx=0;
#By:0-no, 1-yes
By=0;
#Bz:0-no, 1-yes
Bz=1;
#Bzc:0-no, 1-yes
Bzc=0;
#scale:0-h;1-autoscale
sc=0;

filepar = "par_plot.plt"
filein  = "output_plane.txt"
filein1 = "output_plane1.txt"

########################################################################

load filepar

a= z/(ncz);

i=12;
n=173+i;

c=step-1;
#c=2;

B_scale=0.050;
Bc_scale=1.0;

set xrange [0:x*1000]
set yrange [0:y*1000]
set xtics 5
set xlabel "x [mm]" 
set ylabel "y [mm]" 
set size ratio 	-1
set key off
set style data lines 
set palette rgbformulae 30,31,32
lwidth=1
set samples 200
set style data lines 
########################################################################

if(BJx==1){
do for [n=0:c]{
if(form==0){
set terminal jpeg
set output sprintf("BJx_plane%.0f.jpeg",n)}
else{
set terminal postscript eps enhanced color lw 1 "Arial" 26
set output sprintf("BJx_plane%.0f.eps",n)}

set cblabel "BJ_x/B_a,max"
set cbrange [-B_scale:B_scale]            

if(title==0){
b = a*(i+1) - (a/2);
set title sprintf("z=%.2f [mm]",b*1000)}
else{
ts= (n+1)*dt;
b= Bamax * sin (2 * pi * f * ts);
set title sprintf("B_a=%.2f [mT]",b*1000)}

plot filein index n using ($2*1000):($3*1000):($5/Bamax) with image  notitle
}}

########################################################################
if(BJy==1){
do for [n=0:c]{
if(form==0){
set terminal jpeg
set output sprintf("BJy_plane%.0f.jpeg",n)}
else{
set terminal postscript eps enhanced color lw 1 "Arial" 26
set output sprintf("BJy_plane%.0f.eps",n)}

set cblabel "BJ_y/B_a,max"
set cbrange [-B_scale:B_scale]             

if(title==0){
b = a*(i+1) - (a/2);
set title sprintf("z=%.2f [mm]",b*1000)}
else{
ts= (n+1)*dt;
b= Bamax * sin (2 * pi * f * ts);
set title sprintf("B_a=%.2f [mT]",b*1000)}

plot filein index n using ($2*1000):($3*1000):($6/Bamax) with image  notitle 
}}

########################################################################

if(BJz==1){
do for [n=0:c]{
if(form==0){
set terminal jpeg
set output sprintf("BJz_plane%.0f.jpeg",n)}
else{
set terminal postscript eps enhanced color lw 1 "Arial" 26
set output sprintf("BJz_plane%.0f.eps",n)}

set cblabel "BJ_z/B_a,max"       
set cbrange [-B_scale:B_scale] 
#set cbrange [-1:1] 
if(title==0){
b = a*(i+1) - (a/2);
set title sprintf("z=%.2f [mm]",b*1000)}
else{
ts= (n+1)*dt;
b= Bamax * sin (2 * pi * f * ts);
set title sprintf("B_a=%.2f [mT]",b*1000)}

plot filein index n using ($2*1000):($3*1000):($7/Bamax) with image  notitle
}}

########################################################################

if(Bx==1){
do for [n=0:c]{
if(form==0){
set terminal jpeg
set output sprintf("Bx_plane%.0f.jpeg",n)}
else{
set terminal postscript eps enhanced color lw 1 "Arial" 26
set output sprintf("Bx_plane%.0f.eps",n)}

set cblabel "Bx/B_a,max" 
set cbrange [-B_scale:B_scale]            

if(title==0){
b = a*(i+1) - (a/2);
set title sprintf("z=%.2f [mm]",b*1000)}
else{
ts= (n+1)*dt;
b= Bamax * sin (2 * pi * f * ts);
#set title sprintf("B_a=%.2f [mT]",b*1000)}
set title sprintf("time step=%.0f",n)}

plot filein index n using ($2*1000):($3*1000):($8/Bamax) with image  notitle
}}

########################################################################

if(By==1){
do for [n=0:c]{
if(form==0){
set terminal jpeg
set output sprintf("By_plane%.0f.jpeg",n)}
else{
set terminal postscript eps enhanced color lw 1 "Arial" 26
set output sprintf("By_plane%.0f.eps",n)}

set cblabel "By/B_a,max"  
set cbrange [-B_scale:B_scale]             

if(title==0){
b = a*(i+1) - (a/2);
set title sprintf("z=%.2f [mm]",b*1000)}
else{
ts= (n+1)*dt;
b= Bamax * sin (2 * pi * f * ts);
#set title sprintf("B_a=%.2f [mT]",b*1000)}
set title sprintf("time step=%.0f",n)}

plot filein index n using ($2*1000):($3*1000):($9/Bamax) with image  notitle 
}}

########################################################################

if(Bz==1){
do for [n=0:c]{
if(form==0){
set terminal jpeg
set output sprintf("Bz_plane%.0f.jpeg",n)}
else{
set terminal postscript eps enhanced color lw 1 "Arial" 26
set output sprintf("Bz_plane%.0f.eps",n)}

set cblabel "Bz/B_a,max" 
set cbrange [-B_scale*0:B_scale] 

if(title==0){
b = a*(i+1) - (a/2);
set title sprintf("z=%.2f [mm]",b*1000)}
else{
ts= (n+1)*dt;
b= Bamax * sin (2 * pi * f * ts);
#set title sprintf("B_a=%.2f [mT]",b*1000)}
set title sprintf("time step=%.0f",n)}

plot filein index n using ($2*1000):($3*1000):($10/Bamax) with image  notitle
}}

########################################################################

if(Bzc==1){
do for [n=0:c]{
set size ratio 1
if(form==0){
set terminal jpeg
set output sprintf("Bz%.0f.jpeg",n)}
else{
set terminal postscript eps enhanced color lw 1 "Arial" 26
set output sprintf("Bz%.0f.eps",n)}

set ylabel "B/Bamax"
set yrange [0:Bc_scale]   
set xtics 1

if(title==0){
b = a*(i+1) - (a/2);
set title sprintf("z=%.2f [mm]",b*1000)}
else{
ts= (n+1)*dt;
b= Bamax * sin (2 * pi * f * ts);
#set title sprintf("B_a=%.2f [mT]",b*1000)} 
set title sprintf("time step=%.0f",n)}

plot filein index n every ::90::104 using ($2*1000):($10/Bamax) lt 1 lc 8 with lines title ""
}}

