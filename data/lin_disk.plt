#!/bin/python

#input
#form:0-jpeg, 1-eps 
form=0;
#title:0-z distance, 1-Ba;
title=1
#Jy component:0-No, 1-Yes
Jy=1 

filepar="par_plot.plt"
filein="cross_sectionX.txt"

########################################################################

load filepar
Ec=1e-4
mi0=pi*4e-7

a= z/(ncz);
#c=9;
rho=Ec/Jo

n=18
c=step-1;

########################################################################
cx=x
cy=y

ts= (n+1)*dt;

print rho
print pi*f*Bamax
print cos(2*pi*f*ts*pi/180)

Js(x)=-(pi*f*Bamax*cos(2*pi*f*ts*pi/180)*(x/1000))/rho
Js(y)=-(pi*f*Bamax*cos(2*pi*f*ts*pi/180)*(y/1000))/rho

########################################################################

set size ratio 1
if(form==0){
set terminal jpeg
set output sprintf("Jx%.0f.jpeg",n)}
else{
set terminal postscript eps enhanced color lw 1 "Arial" 26
set output sprintf("Jx%.0f.eps",n)}
set xtics 2

set xlabel "y [mm]"
set ylabel "J_x[A/m2]"
set key
set style data lines 
set palette rgbformulae 30,31,32
lwidth=1
set samples 200
set style data lines 

set xrange [-(cx/2)*1000:(cx/2)*1000]           
set autoscale

if(title==0){
b = a*(i+1) - (a/2);
set title sprintf("z=%.2f [mm]",b*1000)}
else{}
#ts= (n+1)*dt;
#b= Bamax * sin (2 * pi * f * ts);
#set title sprintf("B_a=%.2f [mT]",b*1000)} 

plot Js(x) lt 1 title "thin disk" , \
	 filein index n using ($3*1000-(cx/2)*1000):($5) lt 1 lc 7 with histeps title "Jx"	

#plot	 filein index n using ($3*1000-(cx/2)*1000):($5) lt 1 lc 7 with histeps title "Jx"	

########################################################################
if(Jy==1){
filein="cross_sectionY.txt"

set size ratio 1
if(form==0){
set terminal jpeg
set output sprintf("Jy%.0f.jpeg",n)}
else{
set terminal postscript eps enhanced color lw 1 "Arial" 26
set output sprintf("Jy%.0f.eps",n)}
set xtics 2

set xlabel "y[mm]"
set ylabel "J_y[A/m2]"
set key off
set style data lines 
set palette rgbformulae 30,31,32
lwidth=1
set samples 200
set style data lines 
set xrange [-(cx/2)*1000:(cx/2)*1000]            
set autoscale

if(title==0){
b = a*(i+1) - (a/2);
set title sprintf("z=%.2f [mm]",b*1000)}
else{}
#ts= (n+1)*dt;
#b= Bamax * sin (2 * pi * f * ts);
#set title sprintf("B_a=%.1f [mT]",b*1000)} 

plot Js(x) lt 1 title "thin disk" , \
	 filein index n using ($2*1000-(cx/2)*1000):($5) lc 7 with histeps title "Jy"	
}

