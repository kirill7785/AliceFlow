# AliceFlow_v0.48

The program Alice_Flow_v0.48 is intended for calculating the temperature field in three-dimensional solid-state models. In some cases, the convective transfer of the coolant is taken into account. Different nonlinearities are taken into account too. The calculation of thermal transient response is supported. To speed up the calculations, the algebraic multigrid method is used. To speed up the non-stationary calculations, Adaptive Locally Refinement Meshes (Alice) are implemented.

System requirements:
1. variant a) 
* 1.1. OS Windows x64
* 1.2. compiller: Visual Studio 2015 community
* 1.3. nvidia cuda toolkit 8.0
* 1.4. nvidia cusp library 0.5.1
* 1.5. compile with option /bigobj
* 1.6. openmp off option is worked
2. variant b)
* 2.1. OS Windows x64
* 2.2. compiller: Visual Studio 2017 (or 2019) community
* 2.3. boost 1.7.0 library
* 2.4. amgcl 12.05.2019 library
* 2.5. compile with option /bigobj
* 2.6. openmp on or off option is worked.
3. variant c)
* 3.1. GNU gcc (g++ 9.1) compiller is supported
* 3.2. with amgcl library 4.08.2019.

* 4.1. For the exe to work, the console solver program needs to download and install the 64 bit version -
microsoft redistributable package x64 VC_redist.x64.exe

* 5.1. To visualize the results of the calculation, you must install
https://www.tecplot.com/products/tecplot-360/
or
https://www.paraview.org/download/

## Algorithms

* 3D Temperature solver on solid blocks. Finite Volume Method. Stacionary or transient.
* 3D cfd Semi Implicit Method for Pressure Linked Equation (SIMPLE [1972]). Stacionary.
* Rhie-Chow [1983].
* Adaptive Local Refinement Mesh (unstructured grid).
* High Resolution Scheme on uneven structural grid: WACEB, SMARTER, SUPER-C etc.
* Algebraic Multigrid Ruge Stueben. BiCGStab. ilu0 smoother. etc.
* Turbulent models: Spalart Allmares, SST K-Omega Menter.
* Bussinesk Approach. Congruate heat transfer.
* User freandly GUI on Embarcadero® Delphi XE8 Version 22.0.19027.8951.


## Primers

Water cooling module

![alt text](https://github.com/kirill7785/algebraic-multigrid-method/blob/master/picture/water_cooling_module.png)

![alt text](https://github.com/kirill7785/algebraic-multigrid-method/blob/master/picture/speed.jpg)

## Adaptive Locally Refinement Mesh  (Alice)

Since the summer of 2016, Adaptive Locally Refinement Meshes (ALICE) are used in this program.

When building the grid in this program, we suppose that it consists of rectangular parallelepipeds with possible local refinement, dictated by the conditions. The essence of the technology of building local adaptive grids is as follows. The initial grid is Cartesian, and all its cells are rectangular parallelepipeds. Then, in accordance with specified criteria, subregions with features of the geometry or solution are distinguished, and in these subregions a smaller grid is built. For definiteness, we consider that the distinguished feature is given by some surface. If the calculated cell lies in the zone of influence of the selected feature, (for example, it intersects the surface) then such a cell is divided into 8 equal cells. Further, if necessary, the cells are divided again, and so on until the required accuracy is achieved. The curvilinear boundary is approximated by steps. Cells of the initial grid are called as level 0 cells, cells obtained by the level 0 grinding are called as level 1 cells, etc.
When generating grids, it is necessary to impose the following additional restriction: in the neighborhood of each cell, there should not be cells that differ from it in size more than twice.

![alt text](https://github.com/kirill7785/algebraic-multigrid-method/blob/master/picture/picALICE.png)

![alt text](https://github.com/kirill7785/algebraic-multigrid-method/blob/master/picture/module.png)

[![Watch the video](https://github.com/kirill7785/algebraic-multigrid-method/blob/master/picture/picFET3.png)](https://yadi.sk/i/Fd9L_d3bAiLD7w)

Wrap around a cube based on high resolution schemes

![alt text](https://github.com/kirill7785/algebraic-multigrid-method/blob/master/picture/Cube%20Flow.png)

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

Kirill Andreevich Ivanov kirill7785@mail.ru MAI/2009
