The program Alice_Flow_v0.38 is intended for calculating the temperature field in three-dimensional solid-state models. In some cases, the convective transfer of the coolant is taken into account. Nonlinearities of different nature are taken into account. The calculation of thermal transient response is supported. To speed up the calculations, the algebraic multigrid method is used. To speed up non-stationary calculations, Adaptive Locally Refinement Mesh (Alice) have been implemented.

Water cooling module

![alt text](https://github.com/kirill7785/algebraic-multigrid-method/blob/master/picture/water_cooling_module.png)

![alt text](https://github.com/kirill7785/algebraic-multigrid-method/blob/master/picture/speed.jpg)

Adaptive Locally Refinement Mesh  (Alice)

Since the summer of 2016, Adaptive Locally Refinement Meshes (ALICE) have been programmed in this program.
When building the grid in this program, we believe that it consists of rectangular parallelepipeds with possible local refinement, dictated by the conditions. The essence of the technology of building local adaptive grids is as follows. The initial grid is Cartesian, and all its cells are rectangular parallelepipeds. Then, in accordance with the specified criteria, subregions with features of the geometry or solution are distinguished, and in these subregions a smaller grid is built, as compared with the initial one. We consider for definiteness that the distinguished feature is given by some surface. If the calculated cell lies in the zone of influence of the selected feature, for example, it intersects the surface, then such a cell is divided into 8 equal cells. Further, if necessary, the cells are divided again, and so on until the required accuracy is achieved. The curvilinear boundary is approximated by steps. Cells of the initial grid are called level 0 cells, cells obtained by level 0 grinding are called level 1 cells, etc.
When generating grids, it is necessary to impose an additional restriction: in the neighborhood of each cell there should not be cells that differ from it in size more than twice.

![alt text](https://github.com/kirill7785/algebraic-multigrid-method/blob/master/picture/picALICE.png)

![alt text](https://github.com/kirill7785/algebraic-multigrid-method/blob/master/picture/module.png)

[![Watch the video](https://github.com/kirill7785/algebraic-multigrid-method/blob/master/picture/picFET3.png)](https://yadi.sk/i/Fd9L_d3bAiLD7w)

Wrap around a cube based on high resolution schemes

![alt text](https://github.com/kirill7785/algebraic-multigrid-method/blob/master/picture/Cube%20Flow.png)

Calculating the thermal resistance of semiconductor structures

![alt text](https://github.com/kirill7785/algebraic-multigrid-method/blob/master/picture/picFET1.png)

![alt text](https://github.com/kirill7785/algebraic-multigrid-method/blob/master/picture/picFET2.png)

Diode thermal resistance calculation


[![Watch the video](https://github.com/kirill7785/algebraic-multigrid-method/blob/master/picture/Diod.png)](https://yadi.sk/i/3AicztJ3sbb97Q)


The rate of convergence of many matrix inversion methods can be significantly increased by using a technique called "multigrid". The multigrid process involves performing early iterations on a fine grid and later iterations on nested more coarse “virtual” grids.
The results are then transferred back from the coarse grids to the original more detailed fine grids. From a numerical point of view, the multigrid approach offers a significant advantage. For a given grid of a specific finite size, iterative methods are effective only at reducing errors, which have a wavelength of the order of a grid step (a tetrahedron edge in tetra meshes, a hexahedron edge in hexa dominant meshes). Thus, while shorter wavelengths errors disappear rather quickly, errors with a longer wavelength, on the order of the size of the computational domain, disappear very slowly (catastrophically slowly). The Multigrid method avoids this problem by using a number of coarse grids (they are nested in the original detailed grid as matryoshka) so that the available components of the error vector with a long wavelength are short-wave (easily suppressed by the usual iterative smoother method) on coarse grids from the nesting hierarchy. In order to avoid the need for coarse-mesh meshing of geometry using a number of different grid steps, this program uses the Algebraic multi-grid method.

Alice_Flow we know the thermal flow

Kirill Andreevich Ivanov kirill7785@mail.ru MAI/2009
