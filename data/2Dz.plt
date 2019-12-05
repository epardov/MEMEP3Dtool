#!/bin/python

#input
#form:0-jpeg, 1-eps 
form=0;
#title:0-z distance, 1-Ba;
title=1;
#Jx:0-no, 1-yes 
Jx=1;
#Jy:0-no, 1-yes 
Jy=1;
#J:0-no, 1-yes 
J=1;
#Bz:0-no, 1-yes 
Bz=0;
#jc:0-no, 1-yes 
jc=0;
#m:0-no, 1-yes 
m=0;
#n-factor:0-no, 1-yes 
nf=0;
#B_other:0-no, 1-yes 
Bo=0;
#T vector:0-no, 1-yes 
Tv=0;
#Dissipation factor:0-no, 1-yes 
U=0;
#Ax:0-no, 1-yes 
Ax=0;
#Ay:0-no, 1-yes 
Ay=0;
#Az:0-no, 1-yes 
Az=0;
#plane from cube:0-no, 1-yes 
S=0;
#scale:0-h;1-autoscale
sc=0;

filepar="par_plot.plt"
filein="output3Dz.txt"

########################################################################

load filepar

a= z/(ncz);

i=12;
n=173+i;

c=step-1;
#c=3;

h=1.300;
g=Bamax*1000;g=40;
j=h/10;
U_scale=5e8;
A_scale=1.0e-4;

a1=1;													#adjustment for the range in the anisotropic case for Jx
if(rel==4){Jc=Jcpe;a1=2}

set xrange [0:x*1000]
set yrange [0:y*1000]
set xtics 4
set ytics 4
set xlabel "x [mm]" 
set ylabel "y [mm]" 
set size ratio 	1
set palette rgbformulae 30,31,32
lwidth=1
set samples 200
set style data lines

########################################################################

if (S==0){
if (Jx==1){
do for [n=0:c]{
if(form==0){
set terminal jpeg
set output sprintf("2Dz_dJx%.0f.jpeg",n)}
else{
set terminal postscript eps enhanced color lw 1 "Arial" 26
set output sprintf("2Dz_dJx%.0f.eps",n)}

#set terminal qt
set cblabel "J_x/J_c"
set key off          
if(sc==0){
set cbrange [-a1*h:a1*h]}
else{
set autoscale}

if(title==0){
b = a*(i+1) - (a/2);
set title sprintf("z=%.2f [mm]",b*1000)}
else{
ts= (n+1)*dt;
b= Bamax * sin (2 * pi * f * ts);
set title sprintf("B_a=%.2f [mT]",b*1000)} 

if(sc==0){
plot filein index n using ($2*1000):($3*1000):($9/Jc) with image  notitle}
else{
plot filein index n using ($2*1000):($3*1000):($9) with image  notitle}
}}

#pause -1

########################################################################


if (Jy==1){
do for [n=0:c]{
if(form==0){

set terminal jpeg
set output sprintf("2Dz_dJy%.0f.jpeg",n)}
else{
set terminal postscript eps enhanced color lw 1 "Arial" 26
set output sprintf("2Dz_dJy%.0f.eps",n)}

#set terminal qt
set cblabel "J_y/J_c"
set key off
set style data lines 
set palette rgbformulae 30,31,32
lwidth=1
set samples 200
set style data lines        
if(sc==0){
set cbrange [-h:h]}
else{
set autoscale}       

if(title==0){
b = a*(i+1) - (a/2);
set title sprintf("z=%.2f [mm]",b*1000)}
else{
ts= (n+1)*dt;
b= Bamax * sin (2 * pi * f * ts);
set title sprintf("B_a=%.2f [mT]",b*1000)} 

if(sc==0){
plot filein index n using ($2*1000):($3*1000):($10/Jc) with image  notitle}
else{
plot filein index n using ($2*1000):($3*1000):($10) with image  notitle}
}}

#pause -1

########################################################################

if(ncz!=1){
do for [n=0:c]{
if(form==0){
set terminal jpeg
set output sprintf("2Dz_dJz%.0f.jpeg",n)}
else{
set terminal postscript eps enhanced color lw 1 "Arial" 26
set output sprintf("2Dz_dJz%.0f.eps",n)}

set cblabel "J_z/J_c"
set key off
set style data lines 
set palette rgbformulae 30,31,32
lwidth=1
set samples 200
set style data lines            
if(sc==0){
set cbrange [-h	:h]}
else{
set autoscale}      

if(title==0){
b = a*(i+1) - (a/2);
set title sprintf("z=%.2f [mm]",b*1000)}
else{
ts= (n+1)*dt;
b= Bamax * sin (2 * pi * f * ts);
set title sprintf("B_a=%.2f [mT]",b*1000)}

if(sc==0){
plot filein index n using ($2*1000):($3*1000):($11/Jc) with image  notitle}
else{
plot filein index n using ($2*1000):($3*1000):($11) with image  notitle}
}}

########################################################################

do for [n=0:c]{
if (J==1){
if(form==0){
set terminal jpeg
set output sprintf("2Dz_J%.0f.jpeg",n)}
else{
set terminal postscript eps enhanced color lw 1 "Arial" 26
set output sprintf("2Dz_J%.0f.eps",n)}

set cblabel "|J|/J_c"
set key off
set style data lines 
set palette rgbformulae 30,31,32
lwidth=1
set samples 200
set style data lines   
if(sc==0){
set cbrange [0:h]}
else{
set autoscale}

if(title==0){
b = a*(i+1) - (a/2);
set title sprintf("z=%.2f [mm]",b*1000)}
else{
ts= (n+1)*dt;
b= Bamax * sin (2 * pi * f * ts);
set title sprintf("B_a=%.2f [mT]",b*1000)}

if(sc==0){
plot filein index n using ($2*1000):($3*1000):(sqrt($9**2+$10**2)/Jc) with image  notitle}
else{
plot filein index n using ($2*1000):($3*1000):(sqrt($9**2+$10**2)) with image  notitle}
}}

########################################################################
if(ncz!=1){
do for [n=0:c]{
if(form==0){
set terminal jpeg
set output sprintf("2Dz_Bx%.0f.jpeg",n)}
else{
set terminal postscript eps enhanced color lw 1 "Arial" 26
set output sprintf("2Dz_Bx%.0f.eps",n)}

set cblabel "B_x[mT]"
set key off
set style data lines 
set palette rgbformulae 30,31,32
lwidth=1
set samples 200
set style data lines
if(sc==0){
set cbrange[-g:g]}            
else{
set autoscale}   


if(title==0){
b = a*(i+1) - (a/2);
set title sprintf("z=%.2f [mm]",b*1000)}
else{
ts= (n+1)*dt;
b= Bamax * sin (2 * pi * f * ts);
set title sprintf("B_a=%.2f [mT]",b*1000)}

plot filein index c using ($2*1000):($3*1000):($22*1000) with image  notitle
}}

########################################################################

if(ncz!=1){
do for [n=0:c]{
if(form==0){
set terminal jpeg
set output sprintf("2Dz_By%.0f.jpeg",n)}
else{
set terminal postscript eps enhanced color lw 1 "Arial" 26
set output sprintf("2Dz_By%.0f.eps",n)}

set ylabel "y [mm]"
set cblabel "B_y[mT]"
set key off
set style data lines 
set palette rgbformulae 30,31,32
lwidth=1
set samples 200
set style data lines   
if(sc==0){
set cbrange[-g:g]}            
else{
set autoscale}            

if(title==0){
b = a*(i+1) - (a/2);
set title sprintf("z=%.2f [mm]",b*1000)}
else{
ts= (n+1)*dt;
b= Bamax * sin (2 * pi * f * ts);
set title sprintf("B_a=%.2f [mT]",b*1000)}

plot filein index c using ($2*1000):($3*1000):($23*1000) with image  notitle 
}}

########################################################################

do for [n=0:c]{
if(Bz==1){
if(form==0){
set terminal jpeg
set output sprintf("2Dz_Bz%.0f.jpeg",n)}
else{
set terminal postscript eps enhanced color lw 1 "Arial" 26
set output sprintf("2Dz_Bz%.0f.eps",n)}

set cblabel "B_z[mT]" 
set key off
set style data lines 
set palette rgbformulae 30,31,32
lwidth=1
set samples 200
set style data lines          
if(sc==0){
set cbrange[-g:g]}            
else{
set autoscale}   

if(title==0){
b = a*(i+1) - (a/2);
set title sprintf("z=%.2f [mm]",b*1000)}
else{
ts= (n+1)*dt;
b= Bamax * sin (2 * pi * f * ts);
set title sprintf("B_a=%.2f [mT]",b*1000)}

plot filein index n using ($2*1000):($3*1000):($24*1000) with image  notitle
}}

########################################################################

if(jc==1){
do for [n=0:c]{
if(form==0){
set terminal jpeg
set output sprintf("2Dz_Jc%.0f.jpeg",n)}
else{
set terminal postscript eps enhanced color lw 1 "Arial" 26
set output sprintf("2Dz_Jc%.0f.eps",n)}

set cblabel "J_c[A/m^2]" 
set key off
set style data lines 
set palette rgbformulae 30,31,32
lwidth=1
set samples 200
set style data lines          
if(sc==0){
set cbrange [-h	:h]}
else{
set autoscale} 

if(title==0){
b = a*(i+1) - (a/2);
set title sprintf("z=%.2f [mm]",b*1000)}
else{
ts= (n+1)*dt;
b= Bamax * sin (2 * pi * f * ts);
set title sprintf("B_a=%.2f [mT]",b*1000)}

plot filein index n using ($2*1000):($3*1000):($13/Jc) with image  notitle
}}

########################################################################

if(nf==1){
do for [n=0:c]{
if(form==0){
set terminal jpeg
set output sprintf("2Dz_N%.0f.jpeg",n)}
else{
set terminal postscript eps enhanced color lw 1 "Arial" 26
set output sprintf("2Dz_N%.0f.eps",n)}

set cblabel "N-factor" 
set key off
set style data lines 
set palette rgbformulae 30,31,32
lwidth=1
set samples 200
set style data lines
if(sc==0){
set cbrange [0:N]} 
else{
set autoscale}             


if(title==0){
b = a*(i+1) - (a/2);
set title sprintf("z=%.2f [mm]",b*1000)}
else{
ts= (n+1)*dt;
b= Bamax * sin (2 * pi * f * ts);
set title sprintf("B_a=%.2f [mT]",b*1000)}

plot filein index n using ($2*1000):($3*1000):($15) with image  notitle
}}

########################################################################

if(m==1){
do for [n=0:c]{
if(form==0){
set terminal jpeg
set output sprintf("2Dz_mx%.0f.jpeg",n)}
else{
set terminal postscript eps enhanced color lw 1 "Arial" 26
set output sprintf("2Dz_mx%.0f.eps",n)}

set cblabel "m_x [A.m2]"
set key off
set style data lines 
set palette rgbformulae 30,31,32
lwidth=1
set samples 200
set style data lines          
set autoscale

if(title==0){
b = a*(i+1) - (a/2);
set title sprintf("z=%.2f [mm]",b*1000)}
else{
ts= (n+1)*dt;
b= Bamax * sin (2 * pi * f * ts);
set title sprintf("B_a=%.2f [mT]",b*1000)} 

plot filein index n using ($2*1000):($3*1000):($28) with image  notitle
}

########################################################################

do for [n=0:c]{
if(form==0){
set terminal jpeg
set output sprintf("2Dz_my%.0f.jpeg",n)}
else{
set terminal postscript eps enhanced color lw 1 "Arial" 26
set output sprintf("2Dz_my%.0f.eps",n)}

set cblabel "m_y [A.m2]"
set key off
set style data lines 
set palette rgbformulae 30,31,32
lwidth=1
set samples 200
set style data lines        
set autoscale      

if(title==0){
b = a*(i+1) - (a/2);
set title sprintf("z=%.2f [mm]",b*1000)}
else{
ts= (n+1)*dt;
b= Bamax * sin (2 * pi * f * ts);
set title sprintf("B_a=%.2f [mT]",b*1000)} 

plot filein index n using ($2*1000):($3*1000):($29) with image  notitle
}

########################################################################

do for [n=0:c]{
if(form==0){
set terminal jpeg
set output sprintf("2Dz_mz%.0f.jpeg",n)}
else{
set terminal postscript eps enhanced color lw 1 "Arial" 26
set output sprintf("2Dz_mz%.0f.eps",n)}

set cblabel "m_z [A.m2]"
set key off
set style data lines 
set palette rgbformulae 30,31,32
lwidth=1
set samples 200
set style data lines        
set autoscale      

if(title==0){
b = a*(i+1) - (a/2);
set title sprintf("z=%.2f [mm]",b*1000)}
else{
ts= (n+1)*dt;
b= Bamax * sin (2 * pi * f * ts);
set title sprintf("B_a=%.2f [mT]",b*1000)} 

plot filein index n using ($2*1000):($3*1000):($30) with image  notitle
}}

########################################################################

if(Bo==1){
do for [n=0:c]{
if(form==0){
set terminal jpeg
set output sprintf("2Dz_Box%.0f.jpeg",n)}
else{
set terminal postscript eps enhanced color lw 1 "Arial" 26
set output sprintf("2Dz_Box%.0f.eps",n)}

set cblabel "Bx_other [mT]"
set key off
set style data lines 
set palette rgbformulae 30,31,32
lwidth=1
set samples 200
set style data lines          
#set cbrange [-g:g] 
set autoscale

if(title==0){
b = a*(i+1) - (a/2);
set title sprintf("z=%.2f [mm]",b*1000)}
else{
ts= (n+1)*dt;
b= Bamax * sin (2 * pi * f * ts);
set title sprintf("B_a=%.2f [mT]",b*1000)} 

plot filein index n using ($2*1000):($3*1000):($38*1000) with image  notitle
}

########################################################################

do for [n=0:c]{
if(form==0){
set terminal jpeg
set output sprintf("2Dz_Boy%.0f.jpeg",n)}
else{
set terminal postscript eps enhanced color lw 1 "Arial" 26
set output sprintf("2Dz_Boy%.0f.eps",n)}

set cblabel "By_other [mT]"
set key off
set style data lines 
set palette rgbformulae 30,31,32
lwidth=1
set samples 200
set style data lines          
#set cbrange [-g:g] 
set autoscale 

if(title==0){
b = a*(i+1) - (a/2);
set title sprintf("z=%.2f [mm]",b*1000)}
else{
ts= (n+1)*dt;
b= Bamax * sin (2 * pi * f * ts);
set title sprintf("B_a=%.2f [mT]",b*1000)} 

plot filein index n using ($2*1000):($3*1000):($39*1000) with image  notitle
}

########################################################################

do for [n=0:c]{
if(form==0){
set terminal jpeg
set output sprintf("2Dz_Boz%.0f.jpeg",n)}
else{
set terminal postscript eps enhanced color lw 1 "Arial" 26
set output sprintf("2Dz_Boz%.0f.eps",n)}

set cblabel "Bz_other [mT]"
set key off
set style data lines 
set palette rgbformulae 30,31,32
lwidth=1
set samples 200
set style data lines          
#set cbrange [-g:g] 
set autoscale

if(title==0){
b = a*(i+1) - (a/2);
set title sprintf("z=%.2f [mm]",b*1000)}
else{
ts= (n+1)*dt;
b= Bamax * sin (2 * pi * f * ts);
set title sprintf("B_a=%.2f [mT]",b*1000)} 

plot filein index n using ($2*1000):($3*1000):($40*1000) with image  notitle
}}

########################################################################

if(Tv==1){
do for [n=0:c]{
if(form==0){
set terminal jpeg
set output sprintf("2Dz_Tz%.0f.jpeg",n)}
else{
set terminal postscript eps enhanced color lw 1 "Arial" 26
set output sprintf("2Dz_Tz%.0f.eps",n)}

set cblabel "T [A/m2]"
set key off
set style data lines 
set palette rgbformulae 30,31,32
lwidth=1
set samples 200
set style data lines          
#set cbrange [-2.0e7:0] 
set autoscale

if(title==0){
b = a*(i+1) - (a/2);
set title sprintf("z=%.2f [mm]",b*1000)}
else{
ts= (n+1)*dt;
b= Bamax * sin (2 * pi * f * ts);
set title sprintf("B_a=%.2f [mT]",b*1000)} 

plot filein index n using ($2*1000):($3*1000):($27) with image  notitle
}}

########################################################################

if(U==1){
do for [n=0:c]{
if(form==0){
set terminal jpeg
set output sprintf("2Dz_U%.0f.jpeg",n)}
else{
set terminal postscript eps enhanced color lw 1 "Arial" 26
set output sprintf("2Dz_U%.0f.eps",n)}

set cblabel "U []"
set key off
set style data lines 
set palette rgbformulae 30,31,32
lwidth=1
set samples 200
set style data lines          
set cbrange [0:U_scale] 
set autoscale

if(title==0){
b = a*(i+1) - (a/2);
set title sprintf("z=%.2f [mm]",b*1000)}
else{
ts= (n+1)*dt;
b= Bamax * sin (2 * pi * f * ts);
set title sprintf("B_a=%.2f [mT]",b*1000)} 

plot filein index n using ($2*1000):($3*1000):($14) with image  notitle
}}

########################################################################

if(Ax==1){
do for [n=0:c]{
if(form==0){
set terminal jpeg
set output sprintf("2Dz_Ax%.0f.jpeg",n)}
else{
set terminal postscript eps enhanced color lw 1 "Arial" 26
set output sprintf("2Dz_Ax%.0f.eps",n)}

set cblabel "Ax"
set key off          
set style data lines 
set palette rgbformulae 30,31,32
lwidth=1
set samples 200
set style data lines          
if(sc==0){
set cbrange [-A_scale:A_scale]} 
else{
set autoscale} 

if(title==0){
b = a*(i+1) - (a/2);
set title sprintf("z=%.2f [mm]",b*1000)}
else{
ts= (n+1)*dt;
b= Bamax * sin (2 * pi * f * ts);
set title sprintf("B_a=%.2f [mT]",b*1000)} 

plot filein index n using ($2*1000):($3*1000):($16) with image  notitle
}}

########################################################################

if(Ay==1){
do for [n=0:c]{
if(form==0){
set terminal jpeg
set output sprintf("2Dz_Ay%.0f.jpeg",n)}
else{
set terminal postscript eps enhanced color lw 1 "Arial" 26
set output sprintf("2Dz_Ay%.0f.eps",n)}

set cblabel "Ay"
set key off          
set style data lines 
set palette rgbformulae 30,31,32
lwidth=1
set samples 200
set style data lines          
if(sc==0){
set cbrange [-A_scale:A_scale]} 
else{
set autoscale} 


if(title==0){
b = a*(i+1) - (a/2);
set title sprintf("z=%.2f [mm]",b*1000)}
else{
ts= (n+1)*dt;
b= Bamax * sin (2 * pi * f * ts);
set title sprintf("B_a=%.2f [mT]",b*1000)} 

plot filein index n using ($2*1000):($3*1000):($17) with image  notitle
}}

########################################################################

if(Az==1){
do for [n=0:c]{
if(form==0){
set terminal jpeg
set output sprintf("2Dz_Az%.0f.jpeg",n)}
else{
set terminal postscript eps enhanced color lw 1 "Arial" 26
set output sprintf("2Dz_Az%.0f.eps",n)}

set cblabel "Az"
set key off          
set style data lines 
set palette rgbformulae 30,31,32
lwidth=1
set samples 200
set style data lines          
if(sc==0){
set cbrange [-A_scale:A_scale]} 
else{
set autoscale} 

if(title==0){
b = a*(i+1) - (a/2);
set title sprintf("z=%.2f [mm]",b*1000)}
else{
ts= (n+1)*dt;
b= Bamax * sin (2 * pi * f * ts);
set title sprintf("B_a=%.2f [mT]",b*1000)} 

plot filein index n using ($2*1000):($3*1000):($18) with image  notitle
}}


########################################################################

}else{

########################################################################
if(form==0){
set terminal jpeg
set output sprintf("2Dz_dJx%.0f.jpeg",n)}
else{
set terminal postscript eps enhanced color lw 1 "Arial" 26
set output sprintf("2Dz_dJx%.0f.eps",n)}

set cblabel "J_x/J_c"
set key off
set style data lines 
set palette rgbformulae 30,31,32
lwidth=1
set samples 200
set style data lines          
if(sc==0){
set cbrange [-h	:h]}
else{
set autoscale}

if(title==0){
b = a*(i+1) - (a/2);
set title sprintf("z=%.2f [mm]",b*1000)}
else{
ts= (n+1)*dt;
b= Bamax * sin (2 * pi * f * ts);
set title sprintf("B_a=%.2f [mT]",b*1000)} 

plot filein index n using ($2*1000):($3*1000):($9/Jc) with image  notitle

########################################################################

if(form==0){
set terminal jpeg
set output sprintf("2Dz_dJy%.0f.jpeg",n)}
else{
set terminal postscript eps enhanced color lw 1 "Arial" 26
set output sprintf("2Dz_dJy%.0f.eps",n)}

set cblabel "J_y/J_c"
set key off
set style data lines 
set palette rgbformulae 30,31,32
lwidth=1
set samples 200
set style data lines        
if(sc==0){
set cbrange [-h	:h]}
else{
set autoscale}        

if(title==0){
b = a*(i+1) - (a/2);
set title sprintf("z=%.2f [mm]",b*1000)}
else{
ts= (n+1)*dt;
b= Bamax * sin (2 * pi * f * ts);
set title sprintf("B_a=%.2f [mT]",b*1000)} 

plot filein index n using ($2*1000):($3*1000):($10/Jc) with image  notitle

########################################################################

if(ncz!=1){
if(form==0){
set terminal jpeg
set output sprintf("2Dz_dJz%.0f.jpeg",n)}
else{
set terminal postscript eps enhanced color lw 1 "Arial" 26
set output sprintf("2Dz_dJz%.0f.eps",n)}

set cblabel "J_z/J_c"
set key off
set style data lines 
set palette rgbformulae 30,31,32
lwidth=1
set samples 200
set style data lines            
set cbrange [-j	:j]      

if(title==0){
b = a*(i+1) - (a/2);
set title sprintf("z=%.2f [mm]",b*1000)}
else{
ts= (n+1)*dt;
b= Bamax * sin (2 * pi * f * ts);
set title sprintf("B_a=%.2f [mT]",b*1000)}

plot filein index n using ($2*1000):($3*1000):($11/Jc) with image  notitle
}

########################################################################

if(form==0){
set terminal jpeg
set output sprintf("2Dz_J%.0f.jpeg",n)}
else{
set terminal postscript eps enhanced color lw 1 "Arial" 26
set output sprintf("2Dz_J%.0f.eps",n)}

set cblabel "|J|/J_c"
set key off
set style data lines 
set palette rgbformulae 30,31,32
lwidth=1
set samples 200
set style data lines   
if(sc==0){
set cbrange [-h	:h]}
else{
set autoscale}

if(title==0){
b = a*(i+1) - (a/2);
set title sprintf("z=%.2f [mm]",b*1000)}
else{
ts= (n+1)*dt;
b= Bamax * sin (2 * pi * f * ts);
set title sprintf("B_a=%.2f [mT]",b*1000)}

plot filein index n using ($2*1000):($3*1000):($12/Jc) with image  notitle

########################################################################
if(ncz!=1){
if(form==0){
set terminal jpeg
set output sprintf("2Dz_Bx%.0f.jpeg",n)}
else{
set terminal postscript eps enhanced color lw 1 "Arial" 26
set output sprintf("2Dz_Bx%.0f.eps",n)}

set cblabel "B_x[mT]"
set key off
set style data lines 
set palette rgbformulae 30,31,32
lwidth=1
set samples 200
set style data lines   
set cbrange [-g:g]            

if(title==0){
b = a*(i+1) - (a/2);
set title sprintf("z=%.2f [mm]",b*1000)}
else{
ts= (n+1)*dt;
b= Bamax * sin (2 * pi * f * ts);
set title sprintf("B_a=%.2f [mT]",b*1000)}

plot filein index n using ($2*1000):($3*1000):($22*1000) with image  notitle
}

########################################################################

if(ncz!=1){
if(form==0){
set terminal jpeg
set output sprintf("2Dz_By%.0f.jpeg",n)}
else{
set terminal postscript eps enhanced color lw 1 "Arial" 26
set output sprintf("2Dz_By%.0f.eps",n)}

set ylabel "y [mm]"
set cblabel "B_y[mT]"
set key off
set style data lines 
set palette rgbformulae 30,31,32
lwidth=1
set samples 200
set style data lines   
set cbrange [-g:g]             

if(title==0){
b = a*(i+1) - (a/2);
set title sprintf("z=%.2f [mm]",b*1000)}
else{
ts= (n+1)*dt;
b= Bamax * sin (2 * pi * f * ts);
set title sprintf("B_a=%.2f [mT]",b*1000)}

plot filein index n using ($2*1000):($3*1000):($23*1000) with image  notitle 
}

########################################################################

if(form==0){
set terminal jpeg
set output sprintf("2Dz_Bz%.0f.jpeg",n)}
else{
set terminal postscript eps enhanced color lw 1 "Arial" 26
set output sprintf("2Dz_Bz%.0f.eps",n)}

set cblabel "B_z[mT]" 
set key off
set style data lines 
set palette rgbformulae 30,31,32
lwidth=1
set samples 200
set style data lines          
set cbrange [-g:g] 

if(title==0){
b = a*(i+1) - (a/2);
set title sprintf("z=%.2f [mm]",b*1000)}
else{
ts= (n+1)*dt;
b= Bamax * sin (2 * pi * f * ts);
set title sprintf("B_a=%.2f [mT]",b*1000)}

plot filein index n using ($2*1000):($3*1000):($24*1000) with image  notitle

########################################################################

if(jc==1){
if(form==0){
set terminal jpeg
set output sprintf("2Dz_Jc%.0f.jpeg",n)}
else{
set terminal postscript eps enhanced color lw 1 "Arial" 26
set output sprintf("2Dz_Jc%.0f.eps",n)}

set cblabel "J_c[A/m^2]" 
set key off
set style data lines 
set palette rgbformulae 30,31,32
lwidth=1
set samples 200
set style data lines          
if(sc==0){
set cbrange [-h	:h]}
else{
set autoscale}

if(title==0){
b = a*(i+1) - (a/2);
set title sprintf("z=%.2f [mm]",b*1000)}
else{
ts= (n+1)*dt;
b= Bamax * sin (2 * pi * f * ts);
set title sprintf("B_a=%.2f [mT]",b*1000)}

plot filein index n using ($2*1000):($3*1000):($13) with image  notitle
}

########################################################################

if(jc==1){
if(form==0){
set terminal jpeg
set output sprintf("2Dz_N%.0f.jpeg",n)}
else{
set terminal postscript eps enhanced color lw 1 "Arial" 26
set output sprintf("2Dz_N%.0f.eps",n)}

set cblabel "N-factor" 
set key off
set style data lines 
set palette rgbformulae 30,31,32
lwidth=1
set samples 200
set style data lines          
set cbrange [0:N] 

if(title==0){
b = a*(i+1) - (a/2);
set title sprintf("z=%.2f [mm]",b*1000)}
else{
ts= (n+1)*dt;
b= Bamax * sin (2 * pi * f * ts);
set title sprintf("B_a=%.2f [mT]",b*1000)}

plot filein index n using ($2*1000):($3*1000):($15) with image  notitle
}}

