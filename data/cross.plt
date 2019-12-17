#!/bin/python

#input
#form:0-jpeg, 1-eps 
form=0;
#Jx:0-no, 1-yes 
Jx=0;
#Jy:0-no, 1-yes 
Jy=0;
#Bz:0-no, 1-yes 
Bz=1;
#ring:0-no, 1-yes
ring=1;

### Applied field regime ###
#Ba:0-no, 1-yes 
Ba=1;

### Transport current regime ###
#It:0-no, 1-yes 
It=0;

filepar = "par_plot.plt"
filein  = "cross_sectionX.txt"
filein1 = "cross_sectionY.txt"

########################################################################

load filepar

d= z/(ncz);

c=step-1;
#c=0;
n=0
h=1.1;
ba_r=4.5
k=1000
mi0=4.0*pi*1e-7

if(ring==1){x=dR}
w =x*k
wx=x*k
wy=y*k
wz=z*k

set xtics 0.25
set key on
set style data lines 
set palette rgbformulae 30,31,32
lwidth=1
set samples 200
set style data lines
set ytics 0.5

################################################################################################################################################

if (Jx==1){
do for [n=0:c]{
set size ratio 1
if(form==0){
set terminal jpeg
set output sprintf("Jx%.0f.jpeg",n)}
else{
set terminal postscript eps enhanced color lw 1 "Arial" 26
set output sprintf("Jx%.0f.eps",n)}

set xlabel "y [mm]"
set ylabel "J_x/J_c"
set xrange [-y*1000/2/wy:y*1000/2/wy]         
set yrange [-h:h]

ts= (n+1)*dt;
b= Bamax * sin (2 * pi * f * ts);
set title sprintf("B_a=%.2f [mT]",b*1000)

plot filein index n using (($3*1000)/wy-y*1000/2/wy):($5/Jc) lt 1 lc "black" with histeps title "model"
}}

################################################################################################################################################
################################################################################################################################################

if (Jy==1){
#w=x
w=z
d=z
d=dR

if(Ba==1){############# Applied field regime ############################################
Ha=Bamax/mi0
#Jcb=Jc*z
Jcb=Jc*x
Hc=Jcb/pi
cc(Hx)=tanh(Hx/Hc)
b(Hx)=0.5*w*1000/(cosh(Hx/Hc))

print "Bc=",Hc*mi0
print "Ba=",Ha*mi0

Js_high(x) = -(Jc*x/abs(x))/Jc
Js_low(x,Ha) = (-1)*(2*Jcb/(d*pi)) * atan(cc(Ha)*x/(sqrt(b(Hx)**2 - x**2)))/Jc
Js(x,Ha)=( abs(x)<b(Hx) ? Js_low(x,Hx) : ( abs(x)<=0.5*w*1000 ? Js_high(x) : 0 ) )}

if(It==1){############ Transport current regime #########################################
n=0
ts= (n+1)*dt;
I=Ismax*sin(2.0 * pi * f * ts);
a(I)=wx*0.5*sqrt(1-(I/Ismax)**2)

print "I ",I
print "Isamx ",Ismax

Js_high(x) = Jc/Jc
Js_low(x,I)  = (2*Jc/(pi)) * atan(sqrt(((wx*0.5)**2 - a(I)**2)/(a(I)**2 - x**2)))/Jc
Js(x,I)      = ( abs(x)<a(I) ? Js_low(x,I) : Js_high(x) )}

#########################################################################################

do for [n=0:c]{
set size ratio 1
if(form==0){
set terminal jpeg
set output sprintf("Jpy%.0f.jpeg",n)}
else{
set terminal postscript eps enhanced color lw 1 "Arial" 26
set output sprintf("Jpy%.0f.eps",n)}

#set xlabel "x/w[-]"
set xlabel "z/w[-]"
set ylabel "J_y/J_c"
set xrange [-x*k/2/wx:x*k/2/wx]          
if(Ba==1){set yrange [-h:h]}
if(It==1){set yrange [0:h]}

ts= (n+1)*dt;
b= Bamax * sin(2 * pi * f * ts);
Hx= b/mi0
I= Ismax * sin(2 * pi * f * ts);
if(Ba==1){set title sprintf("B_a=%.1f [mT]",b*1000)}
if(It==1){set title sprintf("I_t=%.1f [mT]",I)}

#if(ring==1){
#plot filein1 index n using (($2 - (dR/ncx/2.0) + (R-(dR/2)))/dR):($5/Jc) lc "black" with histeps title " model " , \
#		 Js(x*wx,I) lt 1 lc "red" title "Infinite long strip" #if(It==1)
##		 Js(x*wx,Hx) lt 1 lc "red" title "Infinite long strip" #if(Ba==1)
#}
#else{
#plot filein1 index n using (($2*k)/wx-x*k/2/wx):($5/Jc) lc "black" with histeps title " model " , \
#		 Js(x*wx,I) lt 1 lc "red" title "Infinite long strip" #if(It==1)
#		 Js(x*wx,Hx) lt 1 lc "red" title "Infinite long strip" #if(Ba==1)
#}
plot filein1 index n using (($4*k/wz)-z*k/2/wz):($5/Jc) lc "black" with histeps title " model " , \
		 Js(x*wz,Hx) lt 1 lc "red" title "Infinite long strip" #if(Ba==1)
#		 Js(x*wx,I) lt 1 lc "red" title "Infinite long strip" #if(It==1)
}}

################################################################################################################################################
################################################################################################################################################

if (Bz==1){

Jcb=Jc*x#z
Hc=Jcb/pi
Bc=Hc*mi0

print "Bc=",Hc*mi0

#########################################################################################

if(Ba==1){############# Applied field regime ############################################

#w1=0.5*x
w  =Z*k
wx =Z*k
w1 =0.5*Z
Ba =0.015
w1 =0.5*w

Bz(x,a) = (Bc*log( (abs(x)*sqrt(w1**2- a**2) + w1*sqrt(x**2-a**2))  / (a *sqrt(abs(x**2-w1**2))) ))
B(x,Ba) = (a=0.5*w/cosh(Ba/Bc), abs(x)>a  ? Bz(x,a)/Bc : 0)}

#########################################################################################

do for [n=0:c]{
set size ratio 1
if(form==0){
set terminal jpeg
set output sprintf("Bzcross%.0f.jpeg",n)}
else{
set terminal postscript eps enhanced color lw 1 "Arial" 26
set output sprintf("Bz%.0f.eps",n)}

set xlabel "x[mm]"
set ylabel "B_z"

set xrange [-Z*1000/2/w:Z*1000/2/w]
#set xrange [-x*1000/2/wx:x*1000/2/wx]     
if(Ba==1){set yrange [0:ba_r]}

#set xrange [0:Z*k]
#set xtics 1

ts= (n+1)*dt;
b= Bamax * sin (2 * pi * f * ts);
set title sprintf("B_a=%.1f [mT]",b*1000)

plot filein1 index n using (($4*k)/wx - Z*k/2/wx):($10) lc "black" title "model" , \
		 B(x*wx,b)  lt 1 lc "red" title "formula" #if(Ba==1)

#plot filein1 index n using (($4*k)/wx-x*k/2/wx):($10/Bc) lc "black" title "model" , \
#		 B(x*wx,b)  lt 1 lc "red" title "formula" #if(Ba==1)
}}

################################################################################################################################################

