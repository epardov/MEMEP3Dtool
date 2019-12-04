
MEMEP3Dtool
===========

Three-dimensional (3D) electro-magnetic modeling tool for superconductors and normal conductors using the Minimum Electro-Magnetic Entropy Production (MEMEP) variational principle (https://doi.org/10.1016/j.jcp.2017.05.001).

Authors of the code: M. Kapolka and E. Pardo

Address: Institute of Electrical Engineering, Slovak Academy of the Sciences, Bratislava 

Email address: milan.kapolka@savba.sk and enric.pardo@savba.sk

Version of the code: 0.03

Date: 4.12.2019

WEbpage link to code: https://github.com/epardov/MEMEP3Dtool

The names and details about authors cannot be removed. Further development of the code is allowed. You may add the details about future co-authors. 

Licence
--------
GNU General Public License version 3 (GNU GPLv3).

Description of the input parameters
-----------------------------------

Currently, only thin film and bulk models in Cartesian coordinate system with hexahedral cells are supported. The modeling tool can take the following configurations into account:
- Hexahedral bulk superconductor. https://doi.org/10.1088/1361-6668/aa69ed
- Rectangular thin film (hexahedral superconductor with one element in the thickness). https://doi.org/10.1016/j.jcp.2017.05.001
- Multi-filamentary superconductor with normal conductor in between, where each filament is an hexahedron. https://arxiv.org/abs/1605.09610
- Stacks of tapes with several elements in the thickness. https://doi.org/10.1088/1361-6668/ab5aca
- Thin film disks, cylinders or spheres. https://doi.org/10.1016/j.jcp.2017.05.001

The input parameters need to be set in the input.txt file, with the following description (the input examples are in the input_example folder). 

x[m]: width of the sample

xl[m]: width of the metallic part in the sample between two filaments (shape == 3 only; see below for the 'shape' description)

y[m]: length of the sample

z[m]: thickness of the sample

full_matrix: interaction matrix: 0 - with symmetry for uniform mesh, 1 - without symmetry for nonuniform mesh (check RAM memory usage for more than 31x31x31 elements)  

nsucx[-]: number of the cells along the x axis in the superconducting material

nncx[-]: number of normal conductor joints in the striated tape along the x axis (shape == 3 only)

ncy[-]: number of the cells along the y axis

ncz[-]: number of the cells along the z axis (thin film approximation ncz=1) or total number of elements for stack (only in shape=4)

n_tapes[-]: number of superconducting layers in the stack (only for shape=4)

nc_tape[-]: number of cells per superconducting layer in the stack (only in shape=4)

nc_gap[-]: number of gaps between superconducting layers in the stack (only in shape=4)

d_tape[m]: thickness of the superconducting layer in the stack (only in shape=4)

d_gap[m]: thickness of the gap in the stack (only in shape=4)

elc[-]: 0 disable/1 enable, elongated cells in the long sample with aspect ratio greater than 2

tol_elc[-]: tolerance criterion for average vector potential of elongated cells (0.001-default)

Bamax[T]: maximum amplitude of the applied magnetic field (times the void permeability)

Bamax1[T]: maximum amplitude of the applied magnetic cross-field (times the void permeability)

Bshape[-]: waveform of the applied field : 0-sinusoidal, 1-ramp down followed by cross-field of Bamax1 (amplitude of the cross-field) and fi1 (angle of the cross-field), 2-constant ramp (triangular)

Btrape[-]: 0 disable/1 enable, calculation of the magnetic field outside of the sample in a certain plane (B-plane)

Ismax[A]: transport current

rcx_plane[m]: the center position of the B-plane, x component   (only in shape=4 and Btrape=1)

rcy_plane[m]: the center position of the B-plane, y component   (only in shape=4 and Btrape=1)

rcz_plane[m]: the center position of the B-plane, z component   (only in shape=4 and Btrape=1)

x_plane[m]: width of the B plane     (only in shape=4 and Btrape=1)

y_plane[m]: length of the B plane    (only in shape=4 and Btrape=1)

z_plane[m]: thickness of the B plane (only in shape=4 and Btrape=1)

ncx_plane[-]: number of the cells in the B plane along the x axis (only in shape=4 and Btrape=1)

ncy_plane[-]: number of the cells in the B plane along the y axis (only in shape=4 and Btrape=1)

ncz_plane[-]: number of the cells in the B plane along the z axis (only in shape=4 and Btrape=1)

theta[degree]: angle of the applied magnetic field from the x axis to the y axis

fi[degree]: angle of the applied magnetic field from z axis to the x axis

fi1[degree]: angle of the applied cross magnetic field from z axis to the x axis of amplitude Bamax1 (above) and frequency f1 (below). The cross-field is usually perpendicular to the applied magnetic field (fi).  

uni[-]: type of the mesh: 1-uniform (default)

rel[-]: type of power-law E(J) relation: 1-isotropic, 2-Jc(B) of Kim analytical model, 3-Jc(B,theta) interpolated from measured data, 4-with force-free anisotropy (https://doi.org/10.1088/1361-6668/ab016a).

nB[-]: 0 disable/1, enable power-law n(B) interpolated from measured data

measured_points[-]: Jc(B) data, number of magnetic field angles (only in rel=3)

measured_fields[-]: Jc(B) data, number of magnetic field amplitudes per one angle (only in rel=3)

sym[-]: type of minimization: 0-without sectors and symmetry, 1-sectors (default), 2-sectors with symmetry (odd input number of cells in each direction only and only rel=1!)

Ec[V/m]: critical electric field of the power-law E(J) relation

Jo[A/m2]: critical current density (ignored if rel=3,4)

Jol[A/m2]: current density for normal conducting material, defined as Jol=Ec/rho with rho being the resistivity of the normal conducting material (only in shape=1,2,3,4)

rhoR[ohm*m]: effective resistivity of the normal conducting material between filaments (only in shape=3)

dl[m]: width of the normal conducting joint (only in case=3)

Jcpa[A/m2]: parallel critical current density       (only rel=4)

Jcpe[A/m2]: perpendicular critical current density  (only rel=4)

Bo[T]: characteristic magnetic field for the Kim model

N[-]: power law exponent (ignored if nB=1)

Nl[-]: power law exponent for metallic material (1-default)

m[-]: Kim model exponent

f[Hz]: frequency of the applied field of amplitude Bamax (above)

f1[Hz]: frequency of the applied cross-field of amplitude Bamax1 (see above; only in shape=4)

ns[-]: number of time steps per cycle

step[-]: total number of time steps

tolJ[-]: tolerance of the current density (1e-4 default)

shape[-]: geometry of the sample: 0-square/rectangular, 1-disk/ball, 2-cylinder, 3-tape with filaments, 4-stack of tapes

num_threads[-]: number of parallel computing threads (recommended to be the same as the number of threads of the computer)

Building HTStool
-----------------
The code is prepared for Linux (or Unix) operating system and compiles well with the g++ 5.4.0 compiler, and probably also for higher versions. To compile the program and build the executable, run in command-line terminal:

		bash compile.bsh


Running an example
------------------

Once the program is compiled, run with the command-line:

		./main > output.txt

You may edit and save the input file input.txt to change the parameters, or copy the input file from another folder 

Results and graphs
==================

General results
---------------

To see the main results of the program check the output.txt file for AC loss and computing time at the end of file. There are more files with detailed data results and debugging information. It is practical to use the gnuplot graphs below to see the detailed results.

Geometry:
---------

Run the following command in terminal, in order to see geometry structure. You will need to previously install the GNU software gnuplot.

		gnuplot g_1.plt

Screening current density:
--------------------------

Run this command in terminal, in order to see the screening currents in the sample:

		gnuplot 2Dz.plt     

or

		gnuplot 3Dz.plt

Choose options in the 2Dz.plt file: Jx,Jy,J,Bz,...: 0 disable/1 enable to see the corresponding quantities.
Choose options in the 3Dz.plt file: Jx,Jy,J,Jc,...: 0 disable/1 enable to see the corresponding quantities. 

Ac loss:
--------

Run the command below in terminal, in order to see instantaneous AC loss and magnetization in the sample.

		gnuplot loss.plt


Citation
---------

If you use this code or data generated by this code, please cite the following:

E. Pardo and M. Kapolka "3D computation of non-linear eddy currents: variational method and superconducting cubic bulk" J. Comput. Phys., 344:339–363, 2017, https://doi.org/10.1016/j.jcp.2017.05.001



