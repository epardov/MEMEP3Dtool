#!/bin/python


M=1000000

k=1000

########################################################################################

set size ratio 1
set terminal postscript eps enhanced color lw 1 "Arial" 26
set output "A.eps"

set xlabel "d[um]"
set ylabel "ax"
set key on
set palette rgbformulae 30,31,32
lwidth=1
set samples 1000

#set xrange [-0.5*cx*1000/a:0.5*cx*1000/a]
#set yrange [-e*d:e*d]
set xtics 200	
set label "size of edge 20 um" at 100,0.002 
set autoscale

#set style data histeps
set style data lines             
#set title sprintf("B_a=%.2f [mT]",Ba*1000)
plot "A.txt" using ($2*M):($3) lt 1 lc 0 title "full volume-volume" , \
		 "A.txt" using ($2*M):($5) lt 1 lc 1 title "full surface-surface" , \
		 "A.txt" using ($2*M):($4) lt 1 lc 2 title "analytical film" , \
		 "A.txt" using ($2*M):($6) lt 1 lc 3 title "analytical cube" , \
#		 "A.txt" using ($2*M):($7) lt 1 lc 4 title "numerical"

#########################################################################################

#set size ratio 1
#set terminal postscript eps enhanced color lw 1 "Arial" 26
#set output "Ay.eps"

#set xlabel "d[um]"
#set ylabel "ay"
#set key on
#set palette rgbformulae 30,31,32
#lwidth=1
#set samples 1000

##set xrange [-0.5*cx*1000/a:0.5*cx*1000/a]
##set yrange [-e*d:e*d]
#set xtics 200	
#set label "size of edge 243 um" at 100,0.0003 
#set autoscale

##set style data histeps
#set style data lines             
##set title sprintf("B_a=%.2f [mT]",Ba*1000)
#plot "Ay.txt" using ($2*M):($3) lt 1 lc 0 title "full volume-volume" , \
#		 "Ay.txt" using ($2*M):($5) lt 1 lc 1 title "full surface-surface" , \
#		 "Ay.txt" using ($2*M):($4) lt 1 lc 2 title "analytical film" , \
#		 "Ay.txt" using ($2*M):($6) lt 1 lc 3 title "analytical cube" , \
#		 "Ay.txt" using ($2*M):($7) lt 1 lc 4 title "numerical"

#########################################################################################

#set size ratio 1
#set terminal postscript eps enhanced color lw 1 "Arial" 26
#set output "Az.eps"

#set xlabel "d[um]"
#set ylabel "az"
#set key on
#set palette rgbformulae 30,31,32
#lwidth=1
#set samples 1000

##set xrange [-0.5*cx*1000/a:0.5*cx*1000/a]
##set yrange [-e*d:e*d]
#set xtics 200	
#set label "size of edge 243 um" at 100,0.0003 
#set autoscale

##set style data histeps
#set style data lines             
##set title sprintf("B_a=%.2f [mT]",Ba*1000)
#plot "Az.txt" using ($2*M):($3) lt 1 lc 0 title "full volume-volume" , \
#		 "Az.txt" using ($2*M):($5) lt 1 lc 1 title "full surface-surface" , \
#		 "Az.txt" using ($2*M):($4) lt 1 lc 2 title "analytical film" , \
#		 "Az.txt" using ($2*M):($6) lt 1 lc 3 title "analytical cube" , \
#		 "Az.txt" using ($2*M):($7) lt 1 lc 4 title "numerical"

#########################################################################################


set size ratio 1
set terminal postscript eps enhanced color lw 1 "Arial" 26
set output "B.eps"

set xlabel "x[mm]"
set ylabel "a"
set key on
set palette rgbformulae 30,31,32
lwidth=1
set samples 1000

#set xrange [-0.5*cx*1000/a:0.5*cx*1000/a]
#set yrange [-e*d:e*d]
set xtics 1
set label "size of edge 20 um" at 1,0.004 
set label "d=0.1 um" at 1,0.003
#set autoscale
set xrange [0:3.0]
set yrange [0:0.0085]

#set style data histeps
set style data lines             
#set title sprintf("B_a=%.2f [mT]",Ba*1000)
plot "B.txt" using ($2*k):($3) lt 1 lc 0 title "full volume-volume" , \
		 "B.txt" using ($2*k):($4) lt 7 lc 1 with points title "full surface-surface" , \
		 "B.txt" using ($2*k):($5) lt 4 lc 4 with points title "numerical" , \
		 "B.txt" using ($2*k):($6) lt 1 lc 2 title "approx"


