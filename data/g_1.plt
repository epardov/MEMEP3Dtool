#!/bin/bash

######################################################################################################################################

set terminal wxt

set title "Geometry cells"
set xlabel "x"
set ylabel "y"
set zlabel "z"
unset key

splot "./cells.txt" using ($5):($6):($7):(sprintf('%d', $1)) with points pt 5, \
			"./cells.txt" using ($5):($6):($7):(sprintf('%d', $11)) with labels right offset 0,1 point font ',8' , \
	 		"./cells.txt" using ($5):($6):($7):(sprintf('%d', $12)) with labels right offset 2,1 point font ',8' , \
			"./cells.txt" using ($5):($6):($7):(sprintf('%d', $13)) with labels right offset 3,1 point font ',8'

    	 	
#			"./cells.txt" using ($5):($6):($7):(sprintf('%d', $1)) with labels left offset 1,0 point font ',8' , \
#			"./cells.txt" using ($5):($6):($7):(sprintf('%d', $2)) with labels right offset 0,1 point font ',8' , \
#			"./cells.txt" using ($5):($6):($7):(sprintf('%d', $3)) with labels right offset 2,1 point font ',8' , \
#			"./cells.txt" using ($5):($6):($7):(sprintf('%d', $4)) with labels right offset 3,1 point font ',8'

#


pause -1 "Press enter to continue:"
reset

######################################################################################################################################

set terminal wxt

set title "Geometry nodes + surfaces X"
set xlabel "x"
set ylabel "y"
set zlabel "z"
unset key

splot "./nodes.txt" using ($5):($6):($7):(sprintf('%d', $1)) with points pt 5 , \
			"./surfaces_X.txt" using ($5):($6):($7):(sprintf('%d', $1)) with labels left  offset 1,0 point font ',8'

#			"./surfaces_X.txt" using ($5):($6):($7):(sprintf('%d', $2)) with labels right offset 0,1 point font ',8' , \
#			"./surfaces_X.txt" using ($5):($6):($7):(sprintf('%d', $3)) with labels right offset 2,1 point font ',8' , \
#			"./surfaces_X.txt" using ($5):($6):($7):(sprintf('%d', $4)) with labels right offset 3,1 point font ',8' 
    		      

pause -1 "Press enter to continue:"
reset

######################################################################################################################################

set terminal wxt

set title "Geometry nodes + surface Y"
set xlabel "x"
set ylabel "y"
set zlabel "z"
unset key

splot "./nodes.txt" using ($5):($6):($7):(sprintf('%d', $1)) with points pt 5 , \
			"./surfaces_Y.txt" using ($5):($6):($7):(sprintf('%d', $1)) with labels left  offset 1,0 point font ' '
#			"./surfaces_Y.txt" using ($5):($6):($7) , \
#			"./surfaces_Y.txt" using ($5):($6):($7) 


#			"./surfaces_Y.txt" using ($5):($6):($7):(sprintf('%d', $1)) with labels left  offset 1,0 point font ' ' , \
#	 		"./surfaces_Y.txt" using ($5):($6):($7):(sprintf('%d', $2)) with labels right offset 0,1 point font ' ' , \
#			"./surfaces_Y.txt" using ($5):($6):($7):(sprintf('%d', $3)) with labels right offset 2,1 point font ' ' , \
#			"./surfaces_Y.txt" using ($5):($6):($7):(sprintf('%d', $4)) with labels right offset 3,1 point font ' '
    		      

pause -1 "Press enter to continue:"
reset

######################################################################################################################################

set terminal wxt

set title "Geometry nodes + surface Z"
set xlabel "x"
set ylabel "y"
set zlabel "z"
unset key

splot "./nodes.txt" using ($5):($6):($7):(sprintf('%d', $1)) with points pt 5 , \
			"./surfaces_Z.txt" using ($5):($6):($7):(sprintf('%d', $1)) with labels left  offset 1,0 point font ',8'

#			"./surfaces_Z.txt" using ($5):($6):($7):(sprintf('%d', $2)) with labels right offset 0,1 point font ',8' , \
#			"./surfaces_Z.txt" using ($5):($6):($7):(sprintf('%d', $3)) with labels right offset 2,1 point font ',8' , \
#			"./surfaces_Z.txt" using ($5):($6):($7):(sprintf('%d', $4)) with labels right offset 3,1 point font ',8' 
    		      

pause -1 "Press enter to continue:"
reset

######################################################################################################################################

set terminal wxt

set title "Geometry edges X + surfaces X"
set xlabel "x"
set ylabel "y"
set zlabel "z"
unset key

splot "./edgesX.txt" using ($5):($6):($7):(sprintf('%d', $1)) with points pt 5 , \
			"./edgesX.txt" using ($5):($6):($7):(sprintf('%d', $1)) with labels right offset 0,1 point font ',8'
#			"./surfaces_X.txt" using ($5):($6):($7):(sprintf('%d', $1)) with labels left  offset 1,0 point font ',8'
#			"./edgesX.txt" using ($5):($6):($7):(sprintf('%d', $3)) with labels right offset 2,1 point font ',8' , \
#			"./edgesX.txt" using ($5):($6):($7):(sprintf('%d', $4)) with labels right offset 3,1 point font ',8'
#			"./edgesX.txt" using ($5):($6):($7):(sprintf('%d', $1)) with labels left offset 1,0 point font ',8' , \


    		      

pause -1 "Press enter to continue:"
reset

######################################################################################################################################

set terminal wxt

set title "Geometry edges Y + surface Y"
set xlabel "x"
set ylabel "y"
set zlabel "z"
unset key

splot "./edgesY.txt" using ($5):($6):($7):(sprintf('%d', $1)) with points pt 5 , \
			"./edgesY.txt" using ($5):($6):($7):(sprintf('%d', $1)) with labels right offset 0,1 point font ',8'
#			"./surfaces_Y.txt" using ($5):($6):($7):(sprintf('%d', $1)) with labels left  offset 1,0 point font ',8'
#			"./edgesY.txt" using ($5):($6):($7):(sprintf('%d', $3)) with labels right offset 2,1 point font ',8' , \
#			"./edgesY.txt" using ($5):($6):($7):(sprintf('%d', $4)) with labels right offset 3,1 point font ',8' , \
#			"./edgesY.txt" using ($5):($6):($7):(sprintf('%d', $1)) with labels left offset 1,0 point font ',8' , \



pause -1 "Press enter to continue:"
reset

######################################################################################################################################

set terminal wxt

set title "Geometry edges Z + surface Z"
set xlabel "x"
set ylabel "y"
set zlabel "z"
unset key

splot "./edgesZ.txt" using ($5):($6):($7):(sprintf('%d', $1)) with points pt 5 , \
			"./edgesZ.txt" using ($5):($6):($7):(sprintf('%d', $1)) with labels left offset 1,0 point font ',8'
#			"./surfaces_Z.txt" using ($5):($6):($7):(sprintf('%d', $1)) with labels left  offset 1,0 point font ',8'
#			"./edgesZ.txt" using ($5):($6):($7):(sprintf('%d', $2)) with labels right offset 0,1 point font ',8' , \
#			"./edgesZ.txt" using ($5):($6):($7):(sprintf('%d', $3)) with labels right offset 2,1 point font ',8' , \
#			"./edgesZ.txt" using ($5):($6):($7):(sprintf('%d', $4)) with labels right offset 3,1 point font ',8' , \
#			"./edgesZ.txt" using ($5):($6):($7):(sprintf('%d', $1)) with labels left offset 1,0 point font ',8' , \

    		      

pause -1 "Press enter to continue:"
reset

######################################################################################################################################

#set terminal wxt

#set title "Geometry nodes + cells"
#set xlabel "x"
#set ylabel "y"
#set zlabel "z"
#unset key

#splot "./output_plane.txt" using ($2):($3):($4):(sprintf('%d', $1)) with points pt 5
##			"./surfaces_X.txt" using ($5):($6):($7):(sprintf('%d', $1)) with labels left  offset 1,0 point font ',8' 
##			"./surfaces_X.txt" using ($5):($6):($7):(sprintf('%d', $2)) with labels right offset 0,1 point font ',8' , \
##			"./surfaces_X.txt" using ($5):($6):($7):(sprintf('%d', $3)) with labels right offset 2,1 point font ',8' , \
##			"./surfaces_X.txt" using ($5):($6):($7):(sprintf('%d', $4)) with labels right offset 3,1 point font ',8' 
#    		      

#pause -1 "Press enter to continue:"
#reset

######################################################################################################################################


