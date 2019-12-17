#!/bin/bash

#set key left top


filein="loss.txt"
filepar="par_plot.plt"

load filepar
k=1000;
w=x;

pi=4.0*atan(1.0);						 
mi0=4.0*pi*1e-7;						 		

Ba=Bamax*sin(2*pi*f*dt);

a=Jcpe*w;
b=Jcpe*w*mi0;

#"loss.txt" using ($1*(dt/T)):($3) lt 1 lw 2 lc 0 notitle , \
#"loss1.txt" using ($1*(dt/T)/4):($3) lt 1 lw 2 lc 7 notitle , \
########################################################################

set size ratio 1
set terminal postscript eps enhanced color lw 1 "Arial" 26
set output "loss.eps"
set xlabel "t/T"
set ylabel "P[W]"
set cblabel ""
set key on
set key at 0.80,0.16
set style data linespoints

lwidth=1
set samples 200
set style data lines             
set xrange[0:1.5]
set xtics 0.5
#set autoscale
#plot\
#"loss.txt"  using ($1*(dt/T)):($3) lt 1 lw 2 lc 0 title "fi=0" , \
#"loss3.882.txt" using ($1*(dt/T)):($3) lt 1 lw 2 lc 6 title "fi=80"
#"loss3.881.txt" using ($1*(dt/T)):($3) lt 1 lw 2 lc 7 title "fi=45" , \
#"loss3.8818.txt" using ($1*(dt/T)):($3) lt 1 lw 2 lc 7 title "fi=45" , \


########################################################################

set size ratio 1
set terminal postscript eps enhanced color lw 1 "Arial" 26
set output "magnetization.eps"
set xlabel "B_a/{/Symbol \155}_0J_cw"
set ylabel "M_z/J_cw"
set cblabel "[magnetization]"
set key on
#set key left top
set key at 0.00018,0.35
set style data linespoints
set xtics 0.0001
set ytics 0.2


lwidth=10
set samples 200
set style data lines             
set xrange[-Bamax/b:Bamax/b]
set autoscale

plot\
"loss.txt"  using ($6/b):($12/a) lt 1 lw 2 lc 0 title "fi=45" , \
#"loss3.8818.txt" using ((Bamax*sin(2*pi*f*$1*dt))/b):($10/a) lt 1 lw 2 lc 3 title "fi=45" , \
#"loss3.882.txt" using ($6/b):($12/a) lt 1 lw 2 lc 7 title "fi=80"

set terminal wxt
replot
pause -1 "Press enter"

########################################################################


set size ratio 1
set terminal eps
set output "magnetization1.eps"
set xlabel "B_a [mT]"
set ylabel "M*ea [A/m]"
set cblabel "[magnetization]"
set key off
set style data linespoints

lwidth=1
set samples 200
set style data lines             
set xrange[-Bamax*1000:Bamax*1000]
#set autoscale

#plot "loss.txt" using ($6*1000):($13) lt 1 title "point"



