#!/bin/python

filepar="par_plot.plt"
########################################################################
load filepar
########################################################################

mi0=4.0*pi*1e-7
#step
it=1;
h=1.5;

Jc=Jo;

t=dt*it;
Ba=Bamax*sin(2*pi*f*t);

cx=x;

Ha=Ba/mi0

#########################################################################################

Hd=Jc*z/2;

R=cx*1000/2.0;

a=R/cosh(Ha/Hd);

print "f=",f
print "Ba=",Ba
print "Bamax=",Bamax
print "Ha=",Ha
print "Hd*mu0=",Hd*mi0
print "Ha*mu0=",Ha*mi0
print "Jc=",Jc
print "a=",a
print "R=",R	
print "t=",t

Js_high_positive(x) = Jc/Jc
Js_high_negative(x) = -Jc/Jc

Js_low(x) = (-2*Jc/pi) * atan((x/R)*(sqrt(R**2-a**2))/(sqrt(a**2 - x**2)))/Jc

Js(x)=( abs(x)<=a ? Js_low(x) : x<0 ? Js_high_positive(x) : x>0 ? Js_high_negative(x) : 0)


#########################################################################################

#set size ratio -1
set size
set terminal postscript eps enhanced color lw 1 "Arial" 24
set output "disk.eps"
set tics font ", 24" 

set xlabel "x [mm]" font ",24"
set ylabel "J_y/J_c" font ",24"
set xtics 1
set key on
set palette rgbformulae 30,31,32
lwidth=1
set samples 1000

set xrange [-R:R]
set yrange [-h:h]

#set style data histeps
set style data lines             

plot Js(x) lt 1 title "thin disk" , \
		 "cross_sectionY.txt" index it-1 using ($2*1000-(cx/2)*1000):($5/Jo) lt 1 lc 3 title "n=100"

#########################################################################################
reset
#########################################################################################

filepar="par_plot.plt"

########################################################################
load filepar
########################################################################

mi0=4.0*pi*1e-7

Jc=Jo;

cx=x;

h0=Bamax/mi0;


Hd=Jc*z/2;

R=cx/2.0;

k0=8*R/(3*pi*z);

Ms=Jc*R/3.0

facB=1e-3

print "f=",f
print "Ba=",Ba
print "Bamax=",Bamax
print "Hd*mu0=",Hd*mi0
print "Ha*mu0=",Ha*mi0
print "h0=",h0
print "k0=",k0
print "Jc=",Jc
print "R=",R	
print "t=",t


S(x)   = (1.0/(2.0*x))*(acos(1.0/cosh(x)) + (sinh(x)/cosh(x)**2));

Mi(Ba) = (Ba>=0 ? -k0*(Ba/mi0)*S(Ba/(mi0*Hd)) : 0 )
Ha_down(Ba) = - k0*h0*S(h0/Hd) + k0*(h0 - (Ba/mi0))*S((h0 - (Ba/mi0))/(2*Hd))
Ha_up(Ba)   = 	k0*h0*S(h0/Hd) - k0*(h0 + (Ba/mi0))*S((h0 + (Ba/mi0))/(2*Hd))
line(x)			=  0


set size ratio 1
set terminal eps
set output "magnetization1.eps"
set xlabel "B_a [mT]"
set ylabel "M*ea [A/m]"
set cblabel "[magnetization]"
set key on
set style data linespoints

lwidth=1
set samples 1000
set style data lines             
set xrange[-Bamax/facB:Bamax/facB]

plot Ha_down(x*facB) lt 2 title "thin disk" , \
		 Ha_up(x*facB) lt 2 notitle , \
     Mi(x*facB) lt 2 notitle , \
		 line(Ba) lt 1 lc 7 notitle , \
		 "loss.txt" using ($6/facB):($13/(pi/4)) lt 1 title "n=100"

#     Ms lt 4 title "Ms" , \
#		 Ha_down(x) lt 1 title "thin disk" , \

set term wxt
replot

pause -1




