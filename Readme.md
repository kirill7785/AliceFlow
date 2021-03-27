# AliceFlow_v0.48

The program Alice_Flow_v0.48 is intended for calculating the temperature field in three-dimensional solid-state models. In some cases, the convective transfer of the coolant is taken into account. Different nonlinearities are taken into account too. The calculation of thermal transient response is supported. To speed up the calculations, the algebraic multigrid method is used. To speed up the non-stationary calculations, Adaptive Locally Refinement Meshes (Alice) are implemented.



## Algorithms

* The curvilinear boundary of the computational domain acts as steps. The subgrid resolution method is not implemented (missing). It is recommended to use rectangular       computational areas (volumes 3D).
* 3D Temperature solver for conjugate heat transfer. Finite Volume Method. Stacionary or transient.
  Analitic or load from file or zero velocity component (vx,vy,vz) depend.
  Newton Richman or Stefan Bolcman boundary condition.
* 3D cfd Semi Implicit Method for Pressure Linked Equation (SIMPLE [1972]). Stationary or non-stationary fluid dynamics solver are available.
* Pressure monotonizer S.M. Rhee and W.L. Chow [1984] See Ibrahim Sezai's article.
* Adaptive Local Refinement Mesh (unstructured grid). Coarse MESH or Medium Mesh selector.
* High-resolution schemes for convection approximation on both structured and Adaptive Local Refinement Mesh  irregular meshes: WACEB, SMARTER, SUPER-C, etc.
* Algebraic Multigrid : in house or imported. 
imported amg: 1. CUSP NVIDIA library. 2. AMGCL ddemidov library. 3. amg1r5 (r6) Ruge Stueben[1987].
BiCGStab or FGMRes(m). smoothed aggregation or Ruge Stueben amg.  ilu0 smoother. openMP support. etc.
in house amg: РУМБАv.0.14 BiCGStab or FGMRes(m). PMIS or Ruge Stueben (RS or RS2) amg. truncation operator Prolongation,
level depend threshold, mix floating point precision, ilu0 smoother. openMP support on all operations. etc.
More effective priorite queue on Fibonacci Heap for RS coarsening.
* Turbulent models: Spalart Allmares, SST K-Omega Menter.
* Bussinesk Approach. Congruate heat transfer for Natural convection.
* User freandly GUI on Embarcadero® Delphi XE8 Version 22.0.19027.8951.
* A graph method for solving the heat equation in a nonstationary setting. The nonlinear Stefan-Boltzmann boundary condition and the model of the gap (s) with heat exchange by radiation between the gap boundaries are supported. SLAE is solved by two algebraic multigrid methods to choose from: amg1r5 or Rumba.v.0.14.
* Graphic renderer in GLFW OpenGL in C ++ built into the program code. Visualization of calculation results in 3D, animation of non-stationary CFD calculation results.


## System requirements:
1. variant a) 
* 1.1. OS Windows x64
* 1.2. compiller: Visual Studio 2015 community
* 1.3. nvidia cuda toolkit 8.0
* 1.4. nvidia cusp library 0.5.1
* 1.5. compile with option /bigobj
* 1.6. openmp off option is worked
* 1.7. GLFW OpenGL graphics library
2. variant b)
* 2.1. OS Windows x64
* 2.2. compiller: Visual Studio 2017 (or 2019) community
* 2.3. boost 1.7.0 library
* 2.4. amgcl 10.01.2021 library
* 2.5. compile with option /bigobj
* 2.6. openmp on or off option is worked.
* 2.7. GLFW OpenGL graphics library
3. variant c)
* 3.1. GNU g++ (g++ 9.1) compiller is supported. $ g++ -o AliceFlow_v0_48.exe AliceFlow_v0_48.cpp -fopenmp 2> gcc_log.txt
* 3.2. with amgcl library 4.08.2019.
* 3.3. GLFW OpenGL graphics library

4. For the exe to work, the console solver program needs to download and install the 64 bit version -
microsoft redistributable package x64 VC_redist.x64.exe if you use Microsoft Visual Studio C++ compiller.

5. To visualize the results of the calculation, you must install
https://www.tecplot.com/products/tecplot-360/
or
https://www.paraview.org/download/


## Quick compilation mini Guide

(Deprecated. The program is assembled by the Microsoft Visual Studio 2019 compiler with the written GLFW OpenGL graphics library).

### Windows

You are using a computer running Windows 10 from Microsoft. The following discusses the use of only
free (opensource) software.
2. Install the MinGW GNU g ++ compiler. To do this, you need an Internet connection on your computer.
How to install MinGW compiler (GCC / G ++) Compiller in Windows 10
https://www.youtube.com/watch?v=sXW2VLrQ3Bs
3. Download from GitHub https://github.com/kirill7785/AliceFlow
src folder and copy it to C drive. Navigate to C: \ src folder.
4. Start Windows PowerShell. Navigate to the cd C: \ src source folder.
5. Enter in terminal g++ -o AliceFlow_v0_48.exe AliceFlow_v0_48.cpp 2> gcc_log.txt This will result in the file AliceFlow_v0_48.exe.
6. Create an Alice_EXE folder on your desktop. Place the AliceMesh_v0_45.exe interface program written in Delphi into the Alice_EXE folder
(downloaded from GitHub https://github.com/kirill7785/AliceFlow).
Inside the Alice_EXE folder, create a test_pattern \ solver \ x64 folder hierarchy. Put the executable in the Alice_EXE \ test_pattern \ solver \ x64 folder
file AliceFlow_v0_48.exe.
7. On the GitHub site in the AliceEXE folder there are examples for the AliceMesh_v0_45.exe program. Place them in the AliceEXE folder on your desktop.
8. You can use. Run AliceMesh_v0_45.exe read one of the examples File-> Read. Run the Solve-> Run example.
The AliceMesh_v0_45.exe interface program will automatically call the AliceFlow_v0_48.exe solver program.
9. Install tecplot360 software or a free equivalent - ParaView software. PLT file in tecplot360 (or ParaView)
which was written to disk after the AliceFlow_v0_48.exe program terminated. Insert pictures into the report in MS Word.
10. For those who want to use the parallel version of the program, run mingw-get in the terminal (This requires an Internet connection
on your computer) and select the packages mingw32-pthreads-w32 there.
11. After installation, being in the C: \ src folder, run g ++ -o AliceFlow_v0_48.exe AliceFlow_v0_48.cpp -fopenmp 2> gcc_log.txt to create a parallel version of the program.

### Linux

1. You are working on a computer running Linux Ubuntu-20.04.1-desktop-amdx64. http://releases.ubuntu.com/20.04/
2. The AliceFlow_v0_48 solver program can be compiled and run under Linux OS. AliceMesh_v0_45.exe interface program
works only under Windows OS. written in Delphi.
3.copy to your home folder ~ / src folder from GitHub https://github.com/kirill7785/AliceFlow.
4. Change to the src folder: cd ~ / src. Enter terminal g ++ -o AliceFlow_v0_48 AliceFlow_v0_48.cpp 2> gcc_log.txt
5. If you don't have g ++, install it as prompted by Linux Ubuntu sudo apt get ...
6. Assign execute permission to the file AliceFlow_v0_48: chmod + x AliceFlow_v0_48.
7. Enter the path to the executable file in the PATH variable using the command
export PATH = "$ HOME / src: $ PATH".
8. The src folder contains an example of the premeshin.txt file - this is the input file for the AliceFlow_v0_48 solver that is generated by the AliceMesh_v0_45.exe interface.
9. In the terminal, while in the ~ / src folder, issue the command AliceFlow_v0_48
10. The solver will start and work successfully on the premeshin.txt file. To create premeshin.txt files for your tasks, you need to use
AliceMesh_v0_45.exe interface program running under Windows OS. You can use, for example, VMWare. Program interface
AliceMesh_v0_45.exe is not demanding on PC resources. used only to create a model for calculation.

## Examples

Water cooling module

![alt text](https://github.com/kirill7785/algebraic-multigrid-method/blob/master/picture/water_cooling_module.png)

![alt text](https://github.com/kirill7785/algebraic-multigrid-method/blob/master/picture/speed.jpg)

## Adaptive Locally Refinement Mesh  (Alice)

Since the summer of 2016, Adaptive Locally Refinement Meshes (ALICE) are used in this program.

When building the grid in this program, we suppose that it consists of rectangular parallelepipeds with possible local refinement, dictated by the conditions. The essence of the technology of building local adaptive grids is as follows. The initial grid is Cartesian, and all its cells are rectangular parallelepipeds. Then, in accordance with specified criteria, subregions with features of the geometry or solution are distinguished, and in these subregions a smaller grid is built. For definiteness, we consider that the distinguished feature is given by some surface. If the calculated cell lies in the zone of influence of the selected feature, (for example, it intersects the surface) then such a cell is divided into 8 equal cells. Further, if necessary, the cells are divided again, and so on until the required accuracy is achieved. The curvilinear boundary is approximated by steps. Cells of the initial grid are called as level 0 cells, cells obtained by the level 0 grinding are called as level 1 cells, etc.
When generating grids, it is necessary to impose the following additional restriction: in the neighborhood of each cell, there should not be cells that differ from it in size more than twice.


[![Watch the video](https://github.com/kirill7785/algebraic-multigrid-method/blob/master/picture/picFET3.png)](https://yadi.sk/i/Fd9L_d3bAiLD7w)

Wrap around a cube based on high resolution schemes

![alt text](https://github.com/kirill7785/algebraic-multigrid-method/blob/master/picture/Cube%20Flow.png)

Flow around wing Re = 3164, Pr = 0.7.

![alt_text](https://github.com/kirill7785/AliceFlow/blob/master/picture/speed_around_wing.png)

Modeling natural convection in laboratory condition
Ra=6.4E+7; Pr=0.7; L/H=6.

![alt text](https://github.com/kirill7785/algebraic-multigrid-method/blob/master/picture/Raley_Benar%20Natural%20Convection.png)

Temperature boundary layer on the plate Re=9736; Pr=0.7.

![alt text](https://github.com/kirill7785/algebraic-multigrid-method/blob/master/picture/Blasius%201908.png)

Calculating the thermal resistance of semiconductor structures

![alt text](https://github.com/kirill7785/algebraic-multigrid-method/blob/master/picture/picFET1.png)

![alt text](https://github.com/kirill7785/algebraic-multigrid-method/blob/master/picture/picFET2.png)

Diode thermal resistance calculation


[![Watch the video](https://github.com/kirill7785/algebraic-multigrid-method/blob/master/picture/Diod.png)](https://yadi.sk/i/3AicztJ3sbb97Q)

## Classical Algebraic Multigrid Method (CAMG)

The rate of convergence of many matrix inversion methods can be significantly increased by using a technique called "multigrid". The multigrid process involves performing early iterations on a fine grid, and later, iterations on nested more coarse “virtual” grids.

The results are then transferred back from the coarse grids to the original fine grid, more detailed. From a numerical point of view, the multigrid approach offers a significant advantage. For a given grid of a specific finite size, iterative methods are effective only at reducing errors, which have a wavelength of the order of a grid step (a tetrahedron edge in tetra meshes, a hexahedron edge in hexa dominant meshes). Thus, while shorter wavelengths errors disappear rather quickly, errors with a longer wavelength (on the order of the size of the computational domain) disappear very slowly (catastrophically slowly). The Multigrid method avoids this problem by using a number of coarse grids (they are nested in the original detailed grid as matrioshka) so that the available components of the error vector with a long wavelength are short-wave (easily suppressed by the usual iterative smoother method) on coarse grids from the nesting hierarchy. In order to avoid the need for coarse-mesh meshing of geometry using a number of different grid steps, we use the Algebraic multigrid method

nvidia CUSP 0.5.1 smoothed aggregation amg memory usage comparison:

![alt text](https://github.com/kirill7785/algebraic-multigrid-method/blob/master/picture/Comparison-of-operator-complexity-SAMG-versus-AMG1R5_ru.png)

![alt text](https://github.com/kirill7785/algebraic-multigrid-method/blob/master/picture/Comparison-of-operator-complexity-RUMBAv0-14-versus-amg1r5.png)

![alt text](https://github.com/kirill7785/algebraic-multigrid-method/blob/master/picture/RS%20coarsening.png)

![alt text](https://github.com/kirill7785/algebraic-multigrid-method/blob/master/picture/RS2.png)

![alt text](https://github.com/kirill7785/algebraic-multigrid-method/blob/master/picture/picPaper.png)

# custom visualizer (render)

![alt text](https://github.com/kirill7785/AliceFlow/blob/master/picture/IMG-20201218-WA0000.jpeg)

Kirill Andreevich Ivanov kirill7785@mail.ru MAI/2009
