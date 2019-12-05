#!/bin/python

#input
#form:0-jpeg, 1-eps 
form=0;
#title:0-z distance, 1-Ba;
title=1;
#AC_loss:0-no, 1-yes 
AC=1;
#Mx:0-no, 1-yes 
Mx=0;
#My:0-no, 1-yes 
My=0;
#Mz:0-no, 1-yes 
Mz=1;
#cross_demag:0-no, 1-yes 
cd=1;

filein="loss.txt"
filepar="par_plot.plt"

load filepar

a=1.25;

set size ratio 1
set style data linespoints
lwidth=1
set samples 200
set style data lines   

########################################################################

if (AC==1){
if(form==0){
set terminal jpeg
set output "loss.jpeg"}
else{
set terminal postscript eps enhanced color lw 1 "Arial" 26
set output "loss.eps"}

set xlabel "time [s]"
set ylabel "P"
set y2label "Ba [mT]"
set cblabel "[W]"
set key on
          
set xrange[0:a]
set autoscale

plot "loss.txt" using ($4/T):($5) lt 1 lc "black" title "P_{total}" , \
	   "loss.txt" using ($4/T):($16) lt 1 lc "red" title "P_s" , \
	   "loss.txt" using ($4/T):($17) lt 1 lc "blue" title "P_c"
}

########################################################################

if (Mx==1){
if(form==0){
set terminal jpeg
set output "Mx.jpeg"}
else{
set terminal postscript eps enhanced color lw 1 "Arial" 26
set output "Mx.eps"}
set xlabel "Ba[mT]"
set ylabel "M_x[A/m]"
set cblabel ""
            
set autoscale

plot "loss.txt" using ($6*1000):($12) lt 1 lc 6 title " "
}

########################################################################

if (My==1){
if(form==0){
set terminal jpeg
set output "My.jpeg"}
else{
set terminal postscript eps enhanced color lw 1 "Arial" 26
set output "My.eps"}
set xlabel "Ba[mT]"
set ylabel "M_y[A/m]"
set cblabel ""

set autoscale

plot "loss.txt" using ($6*1000):($13) lt 1 lc 6 title " "
}

########################################################################

if (Mz==1){
if(form==0){
set terminal jpeg
set output "Mz.jpeg"}
else{
set terminal postscript eps enhanced color lw 1 "Arial" 26
set output "Mz.eps"}
set xlabel "Ba[T]"
set ylabel "M_z[A/m]"
set cblabel ""
          
set autoscale

plot "loss.txt" using ($8):($14) lt 1 lc 1 notitle
}
########################################################################

if(cd==1){

k=1000

set size ratio 0.5
if(form==0){
set terminal jpeg
set output "Mz.jpeg"}
else{
set terminal postscript eps enhanced color lw 1 "Arial" 26
set output "Mz.eps"}
set xlabel "t[s]"
set ylabel "B [mT]"
set cblabel ""
set key on
set key top
           
#set xrange[0:T*a]
#set x2range[0:T*a]
set xrange[0.0:500.0]
set yrange[0:1000]
set xtics 0.01
set xtics 50
set grid
#set autoscale y

plot "loss.txt" using ($2-100):($6*k) lt 1 lc 3 title "B_a_x" , \
	   "loss.txt" using ($2-100):($8*k) lt 1 lc 4 title "B_a_z" , \
	   "loss.txt" using ($2-100):($23*k) lt 1 lc 8 title "B_t"

########################################################################


set size ratio 0.5
if(form==0){
set terminal jpeg
set output "Mz_zoom.jpeg"}
else{
set terminal postscript eps enhanced color lw 1 "Arial" 26
set output "Mz_zoom.eps"}
set xlabel "t[s]"
set ylabel "B [mT]"
set cblabel ""
set key on
set key top
          
#set xrange[0:T*a]
#set x2range[0:T*a]
set xrange[400:400.02]
set yrange[-60:60]
set xtics 0.05
set ytics 10
set grid
#set autoscale y

plot "loss.txt" using ($2-100):($6*k) lt 1 lc 3 title "B_a_x" , \
	   "loss.txt" using ($2-100):($8*k) lt 1 lc 4 title "B_a_z" , \
	   "loss.txt" using ($2-100):($23*k) lt 1 lc 7 title "B_t"
}
