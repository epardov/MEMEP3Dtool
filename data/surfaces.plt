#!/bin/python

filepar="par_plot.plt"
load filepar

#2D:0-no, 1-yes
D2=1;
#Ax:0-no, 1-yes
Ax=0;
#Ay:0-no, 1-yes
Ay=1;
#Az:0-no, 1-yes
Az=1;
#Jxs:0-no, 1-yes
Jxs=0;
#Jys:0-no, 1-yes
Jys=0;
#Jzs:0-no, 1-yes
Jzs=0;
#Jx:0-no, 1-yes
Jx=0;
#Jy:0-no, 1-yes
Jy=1;
#Jz:0-no, 1-yes
Jz=1;
#Px:0-no, 1-yes
Px=0;
#Py:0-no, 1-yes
Py=0;
#Pz:0-no, 1-yes
Pz=0;
#Sc:0-no, 1-yes
Sc=1;
#Norm:0-no, 1-yes
Norm=1;

########################################################################
########################################################################


c=step-1;
#c=0
i=0;
t= (i+1)*dt;
Ba = Bamax * sin (2 * pi * f * t)*1e3;
a= z/(ncz);

h=1.3;
g=0.1;
n=0

k=1000


set xrange[0:x*1000]
set yrange[0:y*1000]
set style data lines 
set palette rgbformulae 30,31,32
lwidth=1
set samples 200
set style data lines             
set size ratio -1

######################################################

if(Jx==1){
if(D2==1){
filein="output_Jx.txt"
do for [n=0:c]{
set terminal jpeg
set output sprintf("dJx%.0f.jpeg",n)

set xlabel "x [mm]"
set ylabel "y [mm]"
set cblabel "Jx/J0"
set key off
if(Sc==0){
set cbrange [-g:g]}
else{set autoscale}

set title sprintf("Ba=%f [mT]",Ba)

if(Norm==0){
plot filein index n using ($5*k):($6*k):($8) with image  notitle}
else{plot filein index n using ($5*k):($6*k):($8/Jc) with image  notitle}

}}
else{
filein="output_Jx.txt"
set size ratio -1
set terminal pdf color enhanced font 'Arial, 16'  
set output sprintf("3D_dJx%.0f.pdf",n)

set xlabel "x [mm]"
set ylabel "y [mm]"
set cblabel "Jx/J0"
set key off
if(Sc==0){
set cbrange [-g:g]}
else{set autoscale}         


do for [i=0:ncz-1]{
b = a*(i+1) - (a/2);
set title sprintf("z=%.2f [mm]",b*1000)
plot filein index (i) using ($5*k):($6*k):($8/Jc) with image  notitle
}}}

######################################################

if(Jy==1){
if(D2==1){
filein="output_Jy.txt"
do for [n=0:c]{
set terminal jpeg
set output sprintf("dJy%.0f.jpeg",n)

set xlabel "x [mm]"
set ylabel "y [mm]"
set cblabel "Jy/J0"
set key off
if(Sc==0){
set cbrange [-h:h]}
else{set autoscale}
#set autoscale
set title sprintf("Ba=%f [mT]",Ba)

plot filein index n using ($3):($4):($8/Jc) with image  notitle
}}
else{
filein="output_Jy.txt"
set size ratio -1
set terminal pdf color enhanced font 'Arial, 16'  
set output sprintf("3D_dJy%.0f.pdf",n)

set xlabel "x [mm]"
set ylabel "y [mm]"
set cblabel "Jy/J0"
set key off
if(Sc==1){
set cbrange [-h:h]}
else{set autoscale}         


do for [i=0:ncz-1]{
b = a*(i+1) - (a/2);
set title sprintf("z=%.2f [mm]",b*1000)
plot filein index (i) using ($2*1000):($3*1000):($5/Jc) with image  notitle
}}}

######################################################

if(Jz==1){
if(ncz!=1){
filein="output_Jz.txt"
set terminal jpeg
set output sprintf("dJz%.0f.jpeg",n)

set xlabel "x [mm]"
set ylabel "y [mm]"
set cblabel "Jz/J0"
set key off
if(Sc==0){
set cbrange [-h:h]}
else{set autoscale}

set title sprintf("Ba=%f [mT]",Ba)
plot filein index i using ($3):($4):($8/Jc) with image  notitle
}}

######################################################

#if(Jy==1){
#filein="output_Jy.txt"
#set terminal jpeg
#set output "dJy.jpeg"

#set xlabel "x [mm]"
#set ylabel "y [mm]"
#set cblabel "Jy/J0"
#set key off
#set cbrange [-h	:h]

#set title sprintf("Ba=%f [mT]",Ba)
##plot filein index i using ($2*1000):($3*1000):($5/Jo) with image  notitle
#plot filein index i using ($2):($3):($9) with image  notitle
#}

##################################################################################################################################################################
##################################################################################################################################################################

if(Ax==1){
if(D2==1){
filein="output_Jx.txt"
set terminal jpeg
set output sprintf("av_dAx%.0f.jpeg",n)

set xlabel "x [mm]"
set ylabel "y [mm]"
set cblabel "av_dAx"
set key off
#set cbrange [-k1:k1]
set autoscale

set title sprintf("Ba=%f [mT]",Ba)
plot filein index i using ($3):($4):($10) with image  notitle
}
else{
filein="output_Jx.txt"
set size ratio -1
set terminal pdf color enhanced font 'Arial, 16'  
set output sprintf("3D_dAx%.0f.pdf",n)

set xlabel "x [mm]"
set ylabel "y [mm]"
set cblabel "dAx"
set key off
#set cbrange [-k1:k1]          
set autoscale

do for [i=0:ncz-1]{
b = a*(i+1) - (a/2);
set title sprintf("z=%.2f [mm]",b*1000)
plot filein index (i) using ($2*1000):($3*1000):($7) with image  notitle
}}}

##################################################################################################################################################################

if(Ay==1){
if(D2==1){

do for [n=0:c]{
filein="output_Jy.txt"
set terminal jpeg
set output sprintf("av_dAy%.0f.jpeg",n)

set size ratio 1 

set xlabel "x [mm]"
set ylabel "y [mm]"
set cblabel "dAy"
set key off
#set cbrange [k1:k1]
set autoscale

set title sprintf("Ba=%f [mT]",Ba)
#12 snAy
#10 av_da
#9dJy
#plot filein index i using ($5*1000):($6*1000):($11) with image  notitle
plot filein index n using ($3):($4):($10) with image  notitle
}
}
else{
filein="output_Jy.txt"
set size ratio -1
set terminal pdf color enhanced font 'Arial, 16'  
set output sprintf("3D_dAy%.0f.pdf",n)

set xlabel "x [mm]"
set ylabel "y [mm]"
set cblabel "dAy"
set key off
set cbrange [-k:k]          

do for [i=0:ncz-1]{
b = a*(i+1) - (a/2);
set title sprintf("z=%.2f [mm]",b*1000)
plot filein index (i) using ($2*1000):($3*1000):($7) with image  notitle
}}}

##################################################################################################################################################################

if(Az==1){
if(D2==1){
filein="output_Jz.txt"
set terminal jpeg
set output sprintf("av_dAz%.0f.jpeg",n)

set size ratio 1

set xlabel "x [mm]"
set ylabel "y [mm]"
set cblabel "av_dAz"
set key off
#set cbrange [-h	:h]
set autoscale

set title sprintf("Ba=%f [mT]",Ba)
#11 snAz
#10 av_dAz
#9dJz
plot filein index i using ($3):($4):($10) with image  notitle
}
else{
filein="output_Jz.txt"
set size ratio -1
set terminal pdf color enhanced font 'Arial, 16'  
set output sprintf("3D_dAz%.0f.pdf",n)

set xlabel "x [mm]"
set ylabel "y [mm]"
set cblabel "dAz"
set key off
set cbrange [-k/30:k/30]          

do for [i=0:ncz-1]{
b = a*(i+1) - (a/2);
set title sprintf("z=%.2f [mm]",b*1000)
plot filein index (i) using ($2*1000):($3*1000):($7) with image  notitle
}}}


