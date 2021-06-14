
//#include "cuda_runtime.h"
//#include "device_launch_parameters.h"

//#include <stdio.h>

/*
cudaError_t addWithCuda(int *c, const int *a, const int *b, unsigned int size);

__global__ void addKernel(int *c, const int *a, const int *b)
{
    int i = threadIdx.x;
    c[i] = a[i] + b[i];
}

int main()
{
    const int arraySize = 5;
    const int a[arraySize] = { 1, 2, 3, 4, 5 };
    const int b[arraySize] = { 10, 20, 30, 40, 50 };
    int c[arraySize] = { 0 };

    // Add vectors in parallel.
    cudaError_t cudaStatus = addWithCuda(c, a, b, arraySize);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "addWithCuda failed!");
        return 1;
    }

    printf("{1,2,3,4,5} + {10,20,30,40,50} = {%d,%d,%d,%d,%d}\n",
        c[0], c[1], c[2], c[3], c[4]);

    // cudaDeviceReset must be called before exiting in order for profiling and
    // tracing tools such as Nsight and Visual Profiler to show complete traces.
    cudaStatus = cudaDeviceReset();
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaDeviceReset failed!");
        return 1;
    }

    return 0;
}

// Helper function for using CUDA to add vectors in parallel.
cudaError_t addWithCuda(int *c, const int *a, const int *b, unsigned int size)
{
    int *dev_a = 0;
    int *dev_b = 0;
    int *dev_c = 0;
    cudaError_t cudaStatus;

    // Choose which GPU to run on, change this on a multi-GPU system.
    cudaStatus = cudaSetDevice(0);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?");
        goto Error;
    }

    // Allocate GPU buffers for three vectors (two input, one output)    .
    cudaStatus = cudaMalloc((void**)&dev_c, size * sizeof(int));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed!");
        goto Error;
    }

    cudaStatus = cudaMalloc((void**)&dev_a, size * sizeof(int));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed!");
        goto Error;
    }

    cudaStatus = cudaMalloc((void**)&dev_b, size * sizeof(int));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed!");
        goto Error;
    }

    // Copy input vectors from host memory to GPU buffers.
    cudaStatus = cudaMemcpy(dev_a, a, size * sizeof(int), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy failed!");
        goto Error;
    }

    cudaStatus = cudaMemcpy(dev_b, b, size * sizeof(int), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy failed!");
        goto Error;
    }

    // Launch a kernel on the GPU with one thread for each element.
    addKernel<<<1, size>>>(dev_c, dev_a, dev_b);

    // Check for any errors launching the kernel
    cudaStatus = cudaGetLastError();
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "addKernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
        goto Error;
    }
    
    // cudaDeviceSynchronize waits for the kernel to finish, and returns
    // any errors encountered during the launch.
    cudaStatus = cudaDeviceSynchronize();
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching addKernel!\n", cudaStatus);
        goto Error;
    }

    // Copy output vector from GPU buffer to host memory.
    cudaStatus = cudaMemcpy(c, dev_c, size * sizeof(int), cudaMemcpyDeviceToHost);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy failed!");
        goto Error;
    }

Error:
    cudaFree(dev_c);
    cudaFree(dev_a);
    cudaFree(dev_b);
    
    return cudaStatus;
}
*/

// AliceFlow_v0_57.cpp
// ���������� ������������ �������� 14.03.2021 ��� � ������ ���������� �������� ����
// ��� � ��� �������� ����������. 
// ��������� �������������� cfd �������� 25.12.2020.
// 22.12.2020 - � ��������������� ��������� ������ ������� ��������� ��������������
// ������������� ����� ����� v.0.14 ��� ����� �� ����� ����� ��� � amg1r5.
// 12.12.2020 - OpenGL ������������ ������� �� �++ FWGL. � Z �������.
// 17.11.2020 - ������ ��������� � ��� CAD_stl obj ��� CAD ��������������� 
// ��������� �� *.stl ������. 
// 24.06.2020 - 16 ���� 2020 �������� ������������� ������ - network_T solver.
// 2 ���� 2020 ���������� ������������� ������ �� ���� amg1r5.
// 27 ������� 2020 PMIS aggregator, FF1 amg interpolation. Acomp = 3.0;
// ��� � ������������ ����������� ������ �������� 2006. 
// ������ ���������� PMIS ��� � �� ����� 2004 ���.
// 12 ������ 2020 iluk ������������ ��� amg1r5 ���������.
// 10 ������� 2019 ���������������� ������ �������������� K-Omega SST ������� (RANS). 
// 2 ������� 2019 ���������������� ������ �������������� �������� ��������� (RANS).
// 4 ������� 2019; 03 ������ 2019 �������������� ������������ gcc (g++ 9.1). GNU project. 
// 7-8 ��� 2019 ����������� �������������� ������������� ����� amgcl ������ ��������.
// 6 ������ 2019 ��������������� � visual studio community edition 2019 (open source).
// 25.03.2019 ����� ������������ PVS-Studio 6.0
// 19 �����(03) 2019 ���������� ������������� �� ���� ������.
// 6.05.2018 LINK: fatal error LNK1102: ������������ ������ 2015 VS community.
// �����: ���������� ������� � ������ /bigobj
// ��������� ���������� ��� cuda:
// Tools->Options->Text Editor->File Extension ������ cu � ������ ������ add.
// 9 ���� 2017 ������� �� 64 ������ ����� int64_t.
// 15 ������ 2017 ��������������� � vs community edition 2017 (open source).
// 1 ������� 2016 �������������� �� nvidia cuda 8.0. 
// ���� ����� ������� � ����� ��� ������������� � ������ ����. 
// 11 ������ 2016 ���� ������� cl_agl_amg_v0_14.
// 15 ������� 2015 ����. ������ � Visual Studio 2015.
// AliceFlow_v0_21.cpp
// 15 ������� 2015. ������ � Visual Studio 2013.
// AliceFlow_v0_20.cpp
// 14 ������� 2015 ������������� ���������� ����������������� lusol �
// ilu2 decomposition �� 2 ������.
// 
// AliceFlow_v0_07.cpp �� ������  AliceFlow_v0_06.cpp, �� ������ � LES ������� ��������������.
// ��������� �� ������ ������������ � �������� ������ � ������ �������������� Germano.
// 17 ������ 2013 ����. ���������� ����������������� lusol_.
// 1 ������ 2013. ������ � Visual Studio 2012.
//
// AliceFlow_v0_06.cpp:
// 3D ��������� AliceFlow_v0_06.cpp ��������� �������� AliceFlowv0_05.cpp
// �������� �� ������, ��������� ������������ � ����������������.
// 
// ��������� AliceFlow_v0_05.cpp, 
// ���������� ������� AliceFlow_v0_03.cpp, ����������� 
// ���� �� ������� ����������� ��������������
// � ������ ����������� �� ���� �������.
// � ���� �� ��������� ������������ 
// ������������� ����������������� HEX 
// ��������� �����.
// begin one 17 ��� 2011 ����.
//
// 3D ��������� AliceFlowv0_05.cpp 
// ������������ ��������� ��������:
// 1. ��������� ��������������� �����.
// 2. �������� � ���������� ���� ��������
//    ������ ����������� ��� ������ ������� ����.
// 3. ������ ������� ����.
// 4. �������� ����.
// 5. ������� � ������������ tecplot360.
// begin two 30 ���� 2011 ����.
// begin three 14 ������� 2011 ����. ������ �� Visual Studio 2010.
// begin four 12 ����� 2012 ����. (�� ������� ����� �. ������ - ������� �� ���).
//
// �������� ������� ��������� ����� � ����������� ���� TGF2023_*
// ����� 0.2*120��� �������� 100 �������� (10-20��).

#define NO_OPENGL_GLFW

#ifndef NO_OPENGL_GLFW

#include <GLFW/glfw3.h>

#endif

//#define VISUAL_TUDIO_2008_COMPILLER

#ifdef VISUAL_TUDIO_2008_COMPILLER
#define nullptr NULL // ��� Visual Studio 2008
#endif

// ����������������� � ������ ���� ������ ���������� 
//�������������� ������������ g++ (g++ 9.1) �� GNU.
// g++ AliceFlow_v0_48.cpp -fopenmp 2> gcc_log.txt
// gcc_log.txt - ���� � ���������������� ����������� �����������.
//#define MINGW_COMPILLER 1

#ifdef MINGW_COMPILLER

// ���� float 32 ���� �� ������� - �� ������� ������� ����������� ������������ ��������������� ��������.
// ��� float (���� ���������� ��������������� �������� ����� ����� ����) ���� � ���� ������� ����� ��������
// �� ��������� � ����� double (12 ������ 7).
// ��� float128 128 ��� ������� � 32 ���� ��������� ��� ��� double 64 ����. �� ���� ���������� �������� 
// 31 ������ 38 � ���� double. ������ �� ���������� � double �� �������� � � ���� float128. ������� ������������ 
// ����������� ������ ���������� �� ���� �������.
// �� ��������� ������������� ������������ ������ ��� double.


/*
// �������� ����������
// boost.org/doc/libs/1_65_1/libs/math/example/float128_example.cpp

// ��������� ��������. ��� __float128
//#include <quadmath.h>
#include <boost/cstdfloat.hpp> // For float_64_t, float128_t. Must be first include!
//#include <boost/config.hpp>
#include <boost/multiprecision/float128.hpp>
//#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions.hpp> // For gamma function.
#include <boost/math/constants/calculate_constants.hpp> // For constants pi, e ...
#include <typeinfo> //

#include <cmath> // for pow function
*/
/*
C:\AliceFlow_v0_48_GCC>g++ -std=c++11 -fexceptions -std=gnu++17 -g -fext-numeric-literals -II:\modular-boost\libs\math\include -Ii:\modular-boost -c AliceFlow_v0_48.cpp -o obj\Debug\AliceFlow_v0_48.o
C:\AliceFlow_v0_48_GCC>g++ -o bin\Debug\AliceFlow_v0_48.exe obj\Debug\AliceFlow_v0_48.o -lquadmath
*/

//#include <iostream>

//using namespace boost::multiprecision;

//#include <stdlib.h>
//#include <stdio.h>

#endif



// /fp:except option
//#pragma float_control( except, on )

// ��� std::locale::global(std::locale("en_US.UTF-8"));
// �� ��������.
//#include <locale.h>

//using namespace std;

// ������ ������ ��� ���� cuda c.

#include "cuda_runtime.h"
#include "device_launch_parameters.h"



#include <iostream>

#include <stdio.h>
#include <string.h>

#define NOMINMAX
#include <stdlib.h> 
#include <omp.h> // OpenMP
//-Xcompiler "/openmp" 

//#define NOMINMAX // �����������������
//#include <stdio.h>
//#include <stdlib.h>
//#include <windows.h> // �����������������

// ����������������� � ������ ���� ������ NOT NUMBER �������������.
//#define MY_DEBUG_NOT_NUMBER

// ������ ����� ������ ��������� ���� ��������.
// -1 �� ������������
int BonLevelDrobim = -1;// zmeevik 7.


unsigned int calculation_main_start_time_global_Depend = 0; // ������ ����� ��.

bool CAD_GEOMETRY_OCTREE_MESHGEN = false;
bool SPARSE_MESHER = false;// true;// �������������� ���������� ��������� ����� � ������ ����������.

						   // OpenGL global variables for visualization

bool bmyconvective7248 = false;// ��������� ��� ���.

#ifndef NO_OPENGL_GLFW

const bool binverse_color_black_2_white = false;// ������ ������� ���� �� �����.

GLfloat scale_all = 1000.0f;
int SCREEN_WIDTH = 1280; // 640; // ������ ����� ��� ������ �������.
int SCREEN_HEIGHT = 960; // 480; // ������ ������ ��� ������ �������.

						 // ������� ������ ������, ��� ����������� ������������ ������.
GLfloat halfScreenWidth = SCREEN_WIDTH / 2.0f;
GLfloat halfScreenHeight = SCREEN_HEIGHT / 2.0f;

// ������� �� ��� z ��� �������� ������������ ������.
const int abbys = 500;

// n_render - ����������� ������� ����� pa_opengl ��� -1 ���� ������ ��� pa_opengl �� ��������.
int n_render = -1;

// ������ ����������� �� opengl �� ������� ������ ����������, ���������������� ����.
void keyCallback(GLFWwindow* window, int key, int scancode, int action, int mods);


GLfloat rotationX = 0.0f; // ������� ������ ��� Ox � ��������.
GLfloat rotationY = 0.0f; // ������� ������ ��� Oy � ��������.

						  // ������ � �������� ����� �������� ����������� ������� ���� ��� 
						  // � ����� ������ ������������ ��� ������� �����.
GLfloat cosAlf, cosBet, sinAlf, sinBet;

// OpenGL end global declaration.
#endif

// ���� ������ �� Oz
// ���� ����������� ������ ��������� ������ (2D) ���������� �������� (3D solver).
// ����������� �� ��� Oz ����������� � ���� ������ ��������� �������.
//2D ������ ������ ��������� � ��������� OXY.
bool b_one_cell_z = false;

// ��������� �������� ������������ ��� ������� 
// ��������� AliceFlov_v.0.48 � ������������ ��
// ��������� ����������
// AliceMesh_v.0.45.
double free_debug_parametr1 = 0.0;


// ����� ������� ����������. 
// ��������� �� ���������� ������������.
// ���� ������ ���������� ����� 4 �������.
int number_processors_global_var = 1;

int number_cores() {
	return number_processors_global_var;

	// processor xeon 2630 v4 => 4 thread.
	// ���� ������ ���������� ����� 4 �������.
}

//using namespace System;

// 0 - ������ ������ ��� ����������� ���������� ���������� windows.
#define _CUDA_IS_ACTIVE_ 1

typedef float real_mix_precision;
// ����������� ������ ����� ������ ��� ��������� � ���� unsigned int,
// ������ ��� ������� �� �������� ������ �� �������� ������������ ��� 
// int64_t.
typedef unsigned int integer_mix_precision;

// ������������ ����������.
#define doubleprecision 1
#if doubleprecision == 1
#ifdef MINGW_COMPILLER
// ��������� ��������.
//#define doublereal __float128
//#define doublereal  _Quad
//#define doublereal float128
//#define doublereal float // 32 ���
//#define doublereal double // 64 ���
typedef double doublereal;
#else
//#define doublereal  _Quad
//#define doublereal float128
//#define doublereal double // 64 ���
typedef double doublereal;
//#define doublereal long double //double // ������ ������������� ����� ������� ��������
//#define doublereal Decimal // decimal
#endif
#else
//#define doublereal float //float // ������ ������������� ����� ��������� ��������
typedef float doublereal;
#endif


// ��������� ����� ���� ��������������� �� ������ ������������ � 
// � ������������ Microsoft Visual Studio 2008.
#ifdef VISUAL_TUDIO_2008_COMPILLER

doublereal fmin(doublereal a, doublereal b) {
	doublereal r = a;
	if (b < r) r = b;
	return r;
}

doublereal fmax(doublereal a, doublereal b) {
	doublereal r = a;
	if (b > r) r = b;
	return r;
}

template <typename VType>
bool isfinite(VType x) {
	//return std::isfinite(x);

	if (x != x) {
		return false;
	}
	else {
		return true;
	}
}

#include <math.h>
doublereal expm1(doublereal x) {
	return (exp(x) - 1.0);
}

#endif


// 1 - 64 ������ �����.
// ���� ��������� �������� � ���������� �� ����� ���������
// ��������� � ������� ��������� ������ ���� int  � ��������� 
// ������������ ��� int64_t 64 ������� ������. 
#define doubleintprecision 1


#if doubleintprecision == 1
#include <cinttypes> // ��� ���� 64 ������� ������ ����� int64_t

// �������� !!! ��� ���� int64_t �� �������� ��� ������� ���������� ViennaCL.
// ���������� ������ �������� AMGCL �������� � � ����� int64_t. 
//#define integer int64_t
typedef int64_t integer;
const int64_t big_FIBO_integer_Value = 9223372036854775803; // ��� int64t
#else
//#define integer int
typedef int integer;
const integer big_FIBO_integer_Value = 4294967295; // ��� int
#endif


												   //const integer parabola_MNK = 0; // ����������� �� �������� ����������� �� ������� ���������� ���������.
												   // ��������: line_MNK ����� �� ������������, �.�. ���������������� �����������, �� ������� ������
												   // ������� tgf2023_01 ��� ���������� ��� ����� ������� ����������� �� ���������, ����������� ����� ������������� �� ����������
												   // ����� ��� 40 ���������� �������� � �����. �������� ���� ������������.
												   //const integer line_MNK = 1; // ������� ���� ������� ������, ������ ��������� �� ������� ��� �� ������ ������.
												   // ������ ����� ����� �� ������������� ������������ ! ����� ���������� ����������.
												   //const integer cubic_parabola = 2; // ���������� �������� �� 4 ������ � �� �� �����������.

												   //const unsigned int ADIABATIC_WALL_BC = 0; // �������������� ������.
												   //const unsigned int NEWTON_RICHMAN_BC = 1; // ���������� ������� �������-�������.
												   //const unsigned int STEFAN_BOLCMAN_BC = 2; // ���������� ��������� ������� �������-���������.
												   //const unsigned int MIX_CONDITION_BC = 3; // ��������� ��������� ������� ������ ������ + ������ ��������.

enum class DEFAULT_CABINET_BOUNDARY_CONDITION {
	ADIABATIC_WALL_BC = 0, // �������������, ���������� ������� �������
	NEWTON_RICHMAN_BC = 1, // ����������� ���������� � ��������� ����� �������������, ��� �������� ���������
	STEFAN_BOLCMAN_BC = 2, // ������������ ���������� �� ��������� ������ � �������� ������������
	MIX_CONDITION_BC = 3 // ��������� ������� - ������������ ����������� ���������� � ��������� ����� ������������� � ��������� �� ��������� ������ � �������� ������������
};

//const unsigned int DIRICHLET_FAMILY=1;
//const unsigned int NEIMAN_FAMILY=2;
//const unsigned int NEWTON_RICHMAN_FAMILY=3;
//const unsigned int STEFAN_BOLCMAN_FAMILY=4;

enum class WALL_BOUNDARY_CONDITION {
	DIRICHLET_FAMILY = 1, // ������� �������, �������� ����������� �� ������.
	NEIMAN_FAMILY = 2, // ������������������ ������, �� ������ ����� ���������� ������� �������.
	NEWTON_RICHMAN_FAMILY = 3, // ����������� ��� ����� ����� �� ������ � ��������� ����� �������������, ��� �������� ��������� ������ ������. ��������� �������� ������� ������� - �������.
	STEFAN_BOLCMAN_FAMILY = 4 // ����������� ���������� ��� ����� ����� �� ������  �� ��������� ������ � �������� ������������. ��������� ������� ������� ���������.
};

enum ORDER_INTERPOLATION { parabola_MNK = 0, line_MNK = 1, cubic_parabola = 2 };
enum class ORDER_DERIVATIVE { FIRST_ORDER = 1, SECOND_ORDER = 2 };

enum class INIT_SELECTOR_CASE_CAMG_RUMBAv_0_14 { ZERO_INIT = 0, RANDOM_INIT = 1 };

// ��� ������� solve_Thermal(...);
const integer bARRAYrealesation = 1; // �� ������ ���������������� �������, �� ������� � ��������������� ������ �������� �������� ��������������.
const integer bAVLrealesation = 2; // �� ������ ��� ������.

const integer iGLOBAL_RESTART_LIMIT = 2;// 28.05.2020 ���������� ���� ���(������������ ������ ���������.). // ���� 6;
bool bglobal_restart_06_10_2018 = false;
// ��� ��������� �� �� ����� ���������� �� �� ������� ��� ����� ���� �������� ��� ������������.
bool bglobal_restart_06_10_2018_stagnation[iGLOBAL_RESTART_LIMIT + 1] = { false,false,false };//{false,false,false, false,false,false, false}; // ���� 6;
integer iPnodes_count_shadow_memo = 0;

// ���������� ������ �������������� � ������������� ������, ��� �������� �� ����� ����������. 
doublereal d_GLOBAL_POWER_HEAT_GENERATION_IN_CURRENT_MODEL = 0.0; // ��
doublereal d_my_optimetric1_6_12_2019 = 0.0;// ���������� ���������� ��� �����������.
doublereal d_my_optimetric2_6_12_2019 = 0.0;// ���������� ���������� ��� �����������.
doublereal d_my_optimetric3_6_12_2019 = 0.0;// ���������� ���������� ��� �����������.

enum class LINE_DIRECTIONAL {
	X_LINE_DIRECTIONAL = 0,  // YZ plane
	Y_LINE_DIRECTIONAL = 1,  // XZ plane
	Z_LINE_DIRECTIONAL = 2   // XY plane
};

// ��������� ��������������� xyplot �������� ��� �������.
// 5.01.2018
typedef struct Tpatcher_for_print_in_report {
	// ����������� ������� �� ����������� ������� �� ��� ������� ��
	// ������������ ������� �� ��� �������
	doublereal fminimum, fmaximum;
	// � �������� ����������� ����� �� ��� ��������� ������������ ����.
	LINE_DIRECTIONAL idir; // 0 - X, 1 - Y, 2 - Z.
	Tpatcher_for_print_in_report() {
		fminimum = -1.0e+30;
		fmaximum = 1.0e+30;
		idir = LINE_DIRECTIONAL::Y_LINE_DIRECTIONAL; // 0 - X, 1 - Y, 2 - Z.
	}
} Patcher_for_print_in_report;

Patcher_for_print_in_report pfpir; // xyplot �������.

								   // 9 september 2017.
								   // ������ �� ������������ ����������� ������ � ����� ���������� �������� ������.
								   // �������� �������� �������� ��� ������������ �����, �.�. ����� ������ ��� �������� �� ������, � ���� ���������� 
								   // ������� ������. ������ ������� ������ �������� ������� �� ������������� ���� ���������� ��� ����� ����� � 
								   // ����������� ����������������� �������� ���������� �� ����� �����.
integer ireconstruction_free_construct_alloc = 1; // 0 - off, 1 - on.
												  // ���������� �� �������� � ��������� ����
												  // �� ��������� ������� ������ ���� �� �������.
integer ianimation_write_on = 0; // 0 - off, 1 - on.

								 // ��� ���������� ���������� �� ������ 300 �������� 
								 // ������� �� ��������� vel*rGradual_changes � ����� 
								 // ������ �������������������
								 // � ����������� ��� �� ��������� vel.
								 // 1.0 - �� ������������. �� ������������� 0.1!!!
doublereal rGradual_changes = 1.0;

// ��� ���������� ��������� �� ����� load.txt ��� �����������
// �� ��������� my_multiplyer_velocity_load;
doublereal my_multiplyer_velocity_load = 1.0;// 33.3333; // 14.09.2020

											 // 1 - �������� ������ amg1r6.f ������; 
											 // 0 - �������� amg1r5.f � ����.
integer AMG1R6_LABEL = 0;
integer nrd_LABEL = 1131;
integer nru_LABEL = 1131;
doublereal ecg2_LABEL = 0.25; // strong threshold amg1r5
doublereal ewt2_LABEL = 0.35; // F to F threshold amg1r5
bool b_iluk_amg1r5_LABEL_D = false; // ����
bool b_iluk_amg1r5_LABEL_U = false; // �����
									// stabilization_amg1r5_algorithm:
									// 0 - none(amg1r5); 1 - BiCGStab + amg1r5; 2 - FGMRes+amg1r5;
enum class AMG1R5_OUT_ITERATOR { NONE_only_amg1r5, BiCGStab_plus_amg1r5, FGMRes_plus_amg1r5, Non_Linear_amg1r5 };
AMG1R5_OUT_ITERATOR stabilization_amg1r5_algorithm = AMG1R5_OUT_ITERATOR::BiCGStab_plus_amg1r5; // BiCGStab + amg1r5.

																								// ������������� ��������� �������� ����������.
																								// ������ �������� ��� ���� ��������� �������.
																								// initialization value.
doublereal starting_speed_Vx = 0.0;
doublereal starting_speed_Vy = 0.0;
doublereal starting_speed_Vz = 0.0;

// ������� ����� ��� XY-Plot (variation Plot).
// �� ��������� ������� �����, ����� ������� �������� �����, � ����������� 
// ����� ����� ����� �� ���� ���������� ������������� ������� ���������.
doublereal Tochka_position_X0_for_XY_Plot = 0.0;
doublereal Tochka_position_Y0_for_XY_Plot = 0.0;
doublereal Tochka_position_Z0_for_XY_Plot = 0.0;
LINE_DIRECTIONAL idirectional_for_XY_Plot = LINE_DIRECTIONAL::X_LINE_DIRECTIONAL; // 0 - Ox axis. 

																				  // ��� iVar==TEMP && lw==1 ����� �� ������� ����� ������������ ����� ������������ ����������� ����� V ������� ���������� ����� 0.5K.
bool bPhysics_stop = false;
// ������ ��������� ���������� ����������� ������ ��� ����.
bool bPhysics_PTBSH_memory = false;
// ������ ������ ������������� � ������ ����:
bool bonly_solid_calculation = false;

// ����� ��� ������������� ������������� ������ �� ������������� �����.
// 3 ������� 2015 ����� ����� �������� ����� GUI ������������
// � ����� � ��� ���������� ��������������� �������� � ����� ������ ����.
// ������������ �����
const unsigned int UNEVEN_MUSCL = 1017;  // van Leer (1977)
const unsigned int UNEVEN_SOUCUP = 1018; // MINMOD
const unsigned int UNEVEN_HLPA = 1019;
const unsigned int UNEVEN_SMART = 1020; // Gaskell and Lau (1988)
const unsigned int UNEVEN_WACEB = 1021;
const unsigned int UNEVEN_SMARTER = 1022;
const unsigned int UNEVEN_STOIC = 1023; // Darwish (1993)
const unsigned int UNEVEN_CLAM = 1024;
const unsigned int UNEVEN_OSHER = 1025; // Chakravarthy and Osher (1983)
const unsigned int UNEVEN_VONOS = 1026;
const unsigned int UNEVEN_LPPA = 1027;
const unsigned int UNEVEN_EXPONENTIAL = 1028;
const unsigned int UNEVEN_SUPER_C = 1029;
const unsigned int UNEVEN_ISNAS = 1030;
const unsigned int UNEVEN_CUBISTA = 1031;
const unsigned int UNEVEN_GAMMA = 1032; // ����� � ���������� beta_m
const unsigned int UNEVEN_COPLA = 1033; // 1 08 2015
const unsigned int UNEVEN_SECBC = 1034; // 2 08 2015 Yu et al., (2001b) ��������, �������.
const unsigned int UNEVEN_SGSD = 1035; // 3 08 2015 Li and Tao (2002)
									   // WENO ���������� �� ������������.
const unsigned int UNEVEN_WENO5 = 1036; // 14.01.2021

bool bglobal_first_start_radiation = true;

// ���� �� ������ �������������� ������ ������������� � ������� ����.
bool bglobal_unsteady_temperature_determinant = false;

// ����� ��������� ����������:
// simplemeshgen == 0 ��� unevensimplemeshgen ==1.
// �� ��������� ������������ simplemeshgen == 0.
enum class CONFORMAL_MESH_GENERATOR_SELECTOR {
	SIMPLEMESHGEN_MESHER = 0, // ������������� ��� ����� ������� ��������� ���� ���� �����. ������������� ����� �� ��������� �������� �����.
	UNEVENSIMPLEMESHGEN_MESHER = 1, // ������������� ��� ����� ������� ��������� ���� ���� �����. ��� ����� ������������� ����� � ��� �������  ��������� �������� �����.
	COARSEMESHGEN_MESHER = 2 // ������������� ��� ������� ����������� ���������.
};
CONFORMAL_MESH_GENERATOR_SELECTOR iswitchMeshGenerator = CONFORMAL_MESH_GENERATOR_SELECTOR::SIMPLEMESHGEN_MESHER; // ������� �������� ���������.

																												  // �������������� ��� ������������ �������������.
enum class PHYSICAL_MODEL_SWITCH {
	STEADY_TEMPERATURE = 0, UNSTEADY_TEMPERATURE = 1,
	MESHER_ONLY = 2, CFD_STEADY = 3,
	STEADY_STATIC_STRUCTURAL = 5, STEADY_STATIC_STRUCTURAL_AND_TEMPERATURE = 6,
	SECOND_TEMPERATURE_SOLVER = 7,
	PREOBRAZOVATEL_FOR_REPORT = 8, CFD_UNSTEADY = 9,
	NETWORK_T = 10, NETWORK_T_UNSTEADY = 11, UNSTEADY_STATIC_STRUCTURAL = 12,
	UNSTEADY_STATIC_STRUCTURAL_AND_TEMPERATURE = 13
};
//  0 - thermal only steady state calculation, // 0 - STEADY_TEMPERATURE.
//  1 - thermal only unsteady calculation, // 1 - UNSTEADY_TEMPERATURE.
//  2 - mesh generator only.
//  3 - fluid dynamic steady state.
//  5 - Static Structural (Thermal solver #2)
//  6 - Thermal Stress
//  7 - Unsteady thermal solver #2
//  8 - Visualisation only
//  9 - cfd unsteady fluid dynamic.
// 10 - NETWORK_T steady state. �������� ����� ������� ��������� ����������������. ������������ ���������.
// 11 - NETWORK_T unsteady calculation. �������� ����� ������� ��������� ����������������. �������������� ���������.
// 12 - UNSTEADY STRUCTURAL MECHANICS. �������������� ��������,
// 13 - UNSTEADY STRUCTURAL MECHANICS AND UNSTEADY TEMPERATURE CALCULATION. �������������� �������� ��������� � �������������� ��������������.
PHYSICAL_MODEL_SWITCH steady_or_unsteady_global_determinant = PHYSICAL_MODEL_SWITCH::STEADY_TEMPERATURE;



// ������������ �� ���������� �������� ������������ ��������� �����.
bool b_on_adaptive_local_refinement_mesh = false;
enum class TYPE_ALICE_MESH { ONE_PASS_COARSE_ALICE_MESH = 0, MULTI_PASS_MEDIUM_ALICE_MESH = 1 };
TYPE_ALICE_MESH itype_ALICE_Mesh = TYPE_ALICE_MESH::MULTI_PASS_MEDIUM_ALICE_MESH;// ��� ���� �����.

																				 // ������� ��������� ��� ������
																				 // ���� ��� ��� ����������� �� ��������� ������� (���� �� ��������� ������).
typedef struct TPiecewiseConstantTimeStepLaw {
	doublereal time, timestep, m;
	TPiecewiseConstantTimeStepLaw() {
		time = 0.0; // ������ ������� � �������� ���������� ��������� ��� �� �������, �
		timestep = 0.0; // �������� ���� �� �������, �
		m = 0.0; // ��������� �� ������� ����������� �������� �������� � ������ �������� (������).
	}
} PiecewiseConstantTimeStepLawTimeStepLaw;

// 0 - Linear, 
// 1 - Square Wave,
// 2 - Square Wave 2, 
// 3 - Hot Cold (������ �� ��������� ������� �������, � ����� ��������� �������� �����. ��� �� ������� ������� ���������������.)
// 4 - Piecewise Constant 20.12.2019 - ����� ��������� �������� �������������.
enum class TIME_STEP_lAW_SELECTOR {
	LINEAR = 0, // ��� �� ������� ���������� �� ������ ������������ �������������� ��������� �� ����������� �� �������� ������������� �������� �������.
	SQUARE_WAVE = 1, // ����������� ��������� �������� �������� � � ����������. ����������� �� ��������� ������� �������. ����� � �������� ������������� ��������� tau � ����������� Q.
	SQUARE_WAVE2 = 2, // ������������������ ������ ������ ����������� ��������� �������� �������� � � ����������. ������������ ��� "�����������" �����.
	HOT_COLD = 3, // �������� �������� �������� �� ��������� ������� �������, � ����� ��������� �������� �����. ��� �� ������� ������� ���������������, ��� �������� �������� ������ ��������� �����������.
	PIECEWISE_CONSTANT = 4 // ����� ��������� ����� �� ������� � ���������� ��� �������� �������� �������� �������������.
};

// ��� ������� ������ ��������� ���� �� ������� � �������� �������� �� �������
// ��� �������������� �������������.
typedef struct TTimeStepLaw
{
	// 0 - Linear, 
	// 1 - Square Wave,
	// 2 - Square Wave 2, 
	// 3 - Hot Cold (������ �� ��������� ������� �������, � ����� ��������� �������� �����. ��� �� ������� ������� ���������������.)
	// 4 - Piecewise Constant 20.12.2019 - ����� ��������� �������� �������������.
	TIME_STEP_lAW_SELECTOR id_law;

	doublereal Factor_a_for_Linear; // ����������� ������������ �������������� ���������� = 1.0 + Factor_a_for_Linear.
	doublereal tau; // ������������ �������� ��� Square Wave
					// 06_03_2017 ���������� ����� ���� � �������.
	doublereal Q; // ���������� ��� Square Wave.
				  // ���������� ����� ��� ��������������� Square Wave 2.
				  // off_multiplyer - �������� �������� ������ ����� ���������� � ������ �����. 
				  // ��������� off_multiplyer ��������� ������� �������� �������� �� ���� � ����� � � �����.
	doublereal m1, tau1, tau2, tau_pause, T_all, off_multiplyer;
	integer n_cycle; // 20 ������.
					 // hot cold reshime (double linear)
	doublereal on_time_double_linear; // 3c ��������.

									  // 4-�� ����� ��������� ���� �� �������.
									  // ���������������� ����� ��������� ����� �� ������� � �������� ���������.
	integer n_string_PiecewiseConst;
	PiecewiseConstantTimeStepLawTimeStepLaw* table_law_piecewise_constant;

	TTimeStepLaw() {
		// 0 - Linear, 
		// 1 - Square Wave,
		// 2 - Square Wave 2, 
		// 3 - Hot Cold (���������� �.�.)
		// 4 - Piecewise Constant 20.12.2019 - ����� ��������� �������� �������������.
		id_law = TIME_STEP_lAW_SELECTOR::LINEAR;
		Factor_a_for_Linear = 0.2;
		tau = 60.0E-6; // ������������ �������� ��� Square Wave
					   // 06_03_2017 ���������� ����� ���� � �������.
		Q = 10.0; // ���������� ��� Square Wave.
				  // ���������� ����� ��� ��������������� Square Wave 2.
		m1 = 1.0; tau1 = 0.0; tau2 = 0.0; tau_pause = 0.0; T_all = 0.0; off_multiplyer = 0.0;
		n_cycle = 20; // 20 ������.
					  // hot cold reshime (double linear)
		on_time_double_linear = 3.0; // 3c ��������.
									 // 4 ����� ��������� ���� �� �������.
		n_string_PiecewiseConst = 0;
		table_law_piecewise_constant = nullptr;
	}
} TimeStepLaw;

TimeStepLaw glTSL;

// 24 ������� 2016. 
// ��� ��������� ����� ���������� ����� � ����� 0.14 ��������.
typedef struct TQuickNonlinearBoundaryCondition
{
	doublereal emissivity; // ����������� ������������� �����������.
	doublereal ViewFactor; // ������ ���������.
	doublereal Tamb, dS; // ����������� ����� � ������� �������� �� ������.
	doublereal film_coefficient; // ����������� ����������� ��� ������������ ���������� ������.
	bool bactive;
	bool bStefanBolcman_q_on;
	bool bNewtonRichman_q_on;

	TQuickNonlinearBoundaryCondition() {
		emissivity = 0.8;
		ViewFactor = 1.0; // ������ ���������.
		Tamb = 20.0, dS = 0.0;
		film_coefficient = 3.0;
		bactive = false;
		bStefanBolcman_q_on = false;
		bNewtonRichman_q_on = false;
	}
} QuickNonlinearBoundaryCondition;

QuickNonlinearBoundaryCondition* qnbc = nullptr;
integer iadd_qnbc_maxelm = 0; // ��� ����������� ������
bool b_sign_on_nonlinear_bc = false;


// ������� �� �� SIMPLE ����������.
// ��� ����� ��� ����� ������ ��������� ������� ��� ��������� �������������.
// ������ BiCGStab_internal3 ��������.
bool bSIMPLErun_now_for_temperature = false;
// ���������� �������� SIMPLE ��������� 
// ������� ����� ������������ � ����������.
unsigned int number_iteration_SIMPLE_algorithm = 0; // default - 0
													// ��� ����� ��� ����� ������ ��������� ������� ��� ��������� �������������
													// ��� ������� amg1r5 ���������� ����� � ������������ ����������.
bool bSIMPLErun_now_for_natural_convection = false;
// �������������� ������ ���������� ��� �����������.
doublereal* told_temperature_global_for_HOrelax = nullptr;

/*
��� ���������� ������� ���������� ����� ������������ ����������� �����.
��� ����� ����� ��� ���� �� ��������������� ��-����� �� ���������.
������� �������������� ������ ��������� � ���� ��� ���������������� ����� ������ ���� �����������
��� ��������� ����� ����������� ������� ����������� � ������ �����. � �������� ����������������
������ ������� ��������������.
���������� ���: ������� ���������� ������ ��������� ����� ������ ��������������.
������������ �������� �������������� ������� ����������� � ������ ����� ����. 06.01.2020
*/
bool* sourse2Dproblem = nullptr;
doublereal* conductivity2Dinsource = nullptr;

// �������������� ������ ����������.
bool bHORF = false;
bool bdontstartsolver = false;
doublereal* bPamendment_source_old = nullptr;
doublereal* bsource_term_radiation_for_relax = nullptr;
doublereal* b_buffer_correct_source = nullptr;
// �� ��������� ������������ �����_0.14 �� ������ �� ��������.
bool bfirst_start_nonlinear_process = true;

// ������� �������-������� �� ������� ��� �����������.
// 0 - adiabatic wall, 1 - Newton Richman condition, 2 - Stefan Bolcman condition, 3 - mix condition.
DEFAULT_CABINET_BOUNDARY_CONDITION adiabatic_vs_heat_transfer_coeff = DEFAULT_CABINET_BOUNDARY_CONDITION::ADIABATIC_WALL_BC;
// ��� ���������� ��������� ������� �� ������ ���� ������ ���� ��������� ����� V ������. 
bool breakRUMBAcalc_for_nonlinear_boundary_condition = false;
bool bvacuumPrism = false; // ������� ��������� �����������.
bool bdouble_vacuum_PRISM = false; // ������� ��������� ����������. ��� ����� ��� ��������� ���������� ������� - ��������� �� ������ ������.
bool bBlockStefanBolcman = false; // ���� true �� ��������� ������� ���������.
doublereal film_coefficient = 0.0; // ����������� �����������.
doublereal operating_temperature_for_film_coeff = 20.0; // Tamb for Newton-Richman condition.
														// ���� ��������� ���� ���������� ����� ���������������� ��������� operating_temperature_for_film_coeff
														// �� � ������� �������-������� �� ������� ������� �������� ����� (���������� ������� �������) ��� ���
														// ������� �������� �������������� � ���������� ������� ������� ������� � ����������� ������������ � ������������
														// ��������������� ��������. ����� ����� �������� ������������ ���������� blocker_Newton_Richman.
bool blocker_Newton_Richman = true;


// 0 - ������������ � �������� � �������� ����.
// 1 - ������������ ������ ������� ����.
enum class WHAT_VISIBLE_OPTION { FLUID_AND_SOLID_BODY_VISIBLE = 0, ONLY_SOLID_BODY_VISIBLE = 1 };
WHAT_VISIBLE_OPTION ionly_solid_visible = WHAT_VISIBLE_OPTION::FLUID_AND_SOLID_BODY_VISIBLE;

// ������������ ����� �������������� ������������� ������� � ���������� ��� ��� ������ BiCGStab+ILU2.
// 0 - �������� BiCGStab + ILU2.
// 1 - �������� ����� ���� � ������ ������� ��������������� �������������� ������ amg1r5 (r6).
// 2 - BiCGStab + ADI (Lr1sk).
// 3 - Gibrid: velocity bicgstab + ilu(lfil), Pressure - ����� v0.14.
// 4 - BiCGStab + AINV N.S.Bridson nvidia cusp 0.5.1 library.
// 5 - AMGCL bicgstab+samg ����� �������.
// 6 - Nvidia cusp 0.5.1 library BiCGStab +samg.
// 7 - Algebraic Multigrid ����� v0.14.
integer iswitchsolveramg_vs_BiCGstab_plus_ILU2 = 0; // BiCGStab + ILU2.
													// 0 - BiCGStab+ILU6, 1 - Direct, 2 - ����� 0.14, 3 - amg1r5, 4 - AMGCL_SECONT_T_SOLVER.
													// for Stress Solver
enum class SECOND_T_SOLVER_ID_SWITCH {
	BICGSTAB_PLUS_ILU6_SECOND_T_SOLVER = 0,
	DIRECT_SECOND_T_SOLVER = 1, CAMG_RUMBA_v0_14_SECOND_T_SOLVER = 2,
	AMG1R5_SECOND_T_SOLVER = 3, AMGCL_SECONT_T_SOLVER = 4
};
SECOND_T_SOLVER_ID_SWITCH iswitchsolveramg_vs_BiCGstab_plus_ILU6 = SECOND_T_SOLVER_ID_SWITCH::BICGSTAB_PLUS_ILU6_SECOND_T_SOLVER; // BiCGStab + ILU6.

bool bwait = false; // ���� false, �� �� ���������� getchar().
					// ���� ������ ����� �������� �� 1e-10 �� ��������� ������ ����� ������ ����
const doublereal admission = 1.0e-30; //1.0e-10 // ��� ����������� ���������� ���� ������������ �����.

									  // 1.0e-30 ������� ����. ��� �������� 1.0e-30 ������ �� ����������������. MCB - �� �������������. 27.11.2020.
const doublereal admission_bon_con = 1.0e-10; // ��� ������������� ������ �� ������� ������ ��������� �������.

unsigned int calculation_vorst_seach_time = 0;

// ���� ����������� ������ �������� 
// ����������� � 200 �������� ������� 
// �� ������ ������ �� ����� (������).
// � ������ ���������� ����������� ������ TEMPERATURE_FAILURE_DC
// �� ������� ���������� ��������������� ��������� � ����� ����
// ��� ������������ ������ ���������, �������������� ����� �� ���������.
const doublereal TEMPERATURE_FAILURE_DC = 5000.2;


// �������� ����� ��������� ��������� ������� �������� �������������
// ���������� �������. ��. ��� ����������� � ����� ������ �. ���������.
// BETA_PRECISION 1.0 4/3=1.333333333 6/5=1.2
const doublereal BETA_PRECISION = 1.0;

// ����� ��� ������������� ���������-��������
const unsigned int CR = 1; // ����������-����������
const unsigned int UDS = 2; // ��������������� ������� �������
const unsigned int COMB = 3; // ��������������� 
const unsigned int POLY = 4; // �������������� C. ���������
const unsigned int EXP = 5; // ���������������� �����
const unsigned int BULG = 6; // ����� �.�. ��������� ������� (23) �� ������
const unsigned int POW = 7; // �������������

							// UDS  ��. my_approx_convective.c
// ���������� ������� �����.
int iprefix_Scheme_Flow = UDS;
unsigned int iFLOWScheme = UDS; // ��������������� ������� �������
unsigned int iTEMPScheme = UDS; // ��������������� ������� �������

								// �������� ����� ������ ���������� �������� SIMPLEC
								// SIMPLEC Van Doormal and Raithby, 1984 ���.
								// SIMPLEC ���������� �� SIMPLE ������ � ���� �����:
								// 1. � SIMPLEC �� ������������ ������ ���������� ��� �������� ��� ��������� ��������, �.�. alphaP=1.0.
								// 2. � SIMPLEC ������ ����� ��������������� tau ~ alpha/(1-alpha), � � SIMPLE tau ~ alpha. 
								// � ��������� ��������� ��������� ���������.
enum class SIMPLE_CFD_ALGORITHM { SIMPLE_Carretto = 0, SIMPLEC_Van_Doormal_and_Raithby = 1 };
//const unsigned int SIMPLE_Carretto = 0; // �������� SIMPLE Carretto et al., 1973 ������������ �� ���������.
//const unsigned int SIMPLEC_Van_Doormal_and_Raithby = 1; // �������� SIMPLEC Van Doormal and Raithby, 1984 ���.
SIMPLE_CFD_ALGORITHM iSIMPLE_alg = SIMPLE_CFD_ALGORITHM::SIMPLE_Carretto;// SIMPLE_Carretto SIMPLEC_Van_Doormal_and_Raithby

																		 // �������� ������� ���� ��� ���� ������� ��� �����
																		 // (residual,residual) ��� () ��������� ������������ � R^n.
																		 // �������� ���� ��������� ������ �������� �� �������� ������������� 
																		 // ����. ��������� � ������ ����������� �� ��������� ����� �� �������
																		 // ������ ������������ ������ (FVM). ������� ��� ��������� ����� ��������������
																		 // ������������� ������� � ������ mysolverv0_03.c. ��� ������������� ���������� ��,
																		 // ��� ���� �� ����� ������ ������ ������ ��� ������������ �������� ����������� �������,
																		 // ����� ������� ��� ����������� �������� �������� ������� ��� ������� ������ �������� �������������
																		 // �� ������ ����� �������� ��������� (����������). ��� ����� ���������.
																		 // ��� ���������� ������ �������� ������� �������������������� ������� 1.0e-4
																		 // �������� ���������� ������ �� CFX �� �������.
																		 // ���� ����� �������� �� ������� ���������� ����� 1.0e-3.
																		 // �������� ����� �� ������ �� CFX.
doublereal dterminatedTResudual = 1e-5; // ��� ��� Congruate Gradients � ����� BiCGStab_internal3.

doublereal globalEndTimeUnsteadyTemperatureCalculation = 1.0; // ���������� ����� ��� �������������� ������������� ������������� � ������ ����. 

															  // ������������ ��� ��������� ������� 
															  // ������ ������������ �������������.
bool bsolid_static_only = false;

const integer inumcore = 2; // ����� ���� ����������
const bool bparallelizm_old = false;

// ��������� ������ ������� ������������ ��� ������������ �������� �
// ������������ ������� �� ��� ����������:
typedef struct TPARBOUND
{
	integer ileft_start, ileft_finish; // ������� ����� ��� ����� ������� ���������
	integer iright_start, iright_finish; // ������� ����� ��� ������ ������� ���������
	integer iseparate_start, iseparate_finish; // ������� ����� ��� ������� ����������� ������ ��������� ����� �� ��������.
	bool active; // ���������� ������������.

	TPARBOUND() {
		ileft_start = -1; ileft_finish = -2;
		iright_start = -1; iright_finish = -2;
		iseparate_start = -1; iseparate_finish = -2;
		active = false; // ���������� ������������.
	}
} PARBOUND;


// ��������� ������ ������������ ��� �����������������.
typedef struct TPARDATA {
	integer ncore; // 1, 2, 4, 8.
	integer* inumerate;
	// ��� ��� ncore==2;
	PARBOUND b0;
	// ��� ��� ncore==4;
	PARBOUND b00, b01;
	// ��� ��� ncore==8;
	PARBOUND b000, b001, b010, b011;
	TPARDATA() {
		ncore = 2; // 1, 2, 4, 8.
		inumerate = nullptr;
	}
} PARDATA;

PARDATA nd; // nd - nested desection.

			// ������������ � ����� ��� ��������� �����.
doublereal* rthdsd_no_radiosity_patch = nullptr;

// ��� ��������� (���� bdroblenie4=true)
// ������ �� ����� ������ ����� ��������� � �������� ��������� ��������.
/*
typedef struct TALICE_PARTITION {
bool bdroblenie4;
int iNODE1, iNODE2, iNODE3, iNODE4;
TALICE_PARTITION() {
bdroblenie4=false;
iNODE1=-1; iNODE2=-1; iNODE3=-1; iNODE4=-1;
}
void define_structural_mesh_neighbour(int id) {
bdroblenie4 = false;
iNODE1 = id; iNODE2 = -1; iNODE3 = -1; iNODE4 = -1;
}
} ALICE_PARTITION;
*/



// ���������� �������� ���� ����� � ���������� �����. 
// ��� ������ ��������� �������� ���������� �����������.
const bool b_thermal_source_refinement = true;

#include "adaptive_local_refinement_mesh.cpp" // ����

#ifndef NO_OPENGL_GLFW
// ��� ������������ ��������� ������ � �������������� ���������� OpenGL

// ���������� ����� pa �������������� ��� ������� ����� � ������� ��������� Z-������.
TOCHKA* pa_opengl = nullptr; // �������� ����������� �� ��������������� ����������.
TOCHKA* pa_render = nullptr; // ��������������� ���������� ������������� ������ ��� ��������� scale_all.
doublereal minimum_val_for_render_pic, maximum_val_for_render_pic; // ����������� � ������������ �������� ������� �������� �������������� ��� ���������.
int iNUMBER_FUNC_openGL = 1;
int iCURENT_FUNC_openGL = 0;
doublereal** potent_array_opengl = nullptr;
// ��������� ������� �������������� CFD ��������� ����� ������������� �� ������ ��������,
// �������� ��������.
int iCFD_animation = 0;
int iNUMBER_ANIMATION_FUNCTIONS = -1;
int iCURENT_ANIMATION_FUNCTION = 0;
int iNUMBER_ANIMATION_CADERS = -1;
int iCURENT_ANIMATION_CADER = 0;
doublereal*** animation_sequence_functions_openGL = nullptr;

// ��������� ������������ � ��������� ��������� �����. ��������� 19.12.2020.
// ������ ������ �������� ��� ������������ ������, �� �������� ���������� � ������.
int DrawZbufferColor(int maxelm, int maxnod, int**& nvtx);

// OpenGL end.
#endif

// 21.04.2021
// ���� b_adhesion_Mesh �� �� ��������� ���������� �������� ����� ���� �������� ��������� �� ���.
// ��� ��������� ���������� ��� ������� ��������� body block ��� ������������� ����� ���������� ���������.
// ��� �� �������� � ������ ��������� ����� ����� ����� ��� ����� ������������� ���� ����� ������.
// �� ������ ���� ����� � ����� ���������� b_adhesion_Mesh ���� �������� � false.
bool b_adhesion_Mesh = true;

#include "constr_struct.cpp" // ���������� �������� ������ TEMPER � FLOW
#include "uniformsimplemeshgen.cpp" // �������� ���������

#ifndef NO_OPENGL_GLFW

// ������������ ���������� ������������ ������.
int* mask_color_for_openGL = nullptr;

// ������������� �� �������������� ������� �� 0 �� 1020 
// ������� �������� ����� �� �������� �����.
void set_color_for_render_icolor(int icol)
{

	if ((0 <= icol) && (icol <= 255)) {
		// ����� �������
		glColor3f(0.0, (icol) / 255.0, 1.0);
	}
	else if ((256 <= icol) && (icol <= 510)) {
		//������� - ������
		glColor3f(0, 1.0, (255.0 - (icol - 255.0)) / 255.0);
	}
	else if ((511 <= icol) && (icol <= 765)) {

		// ������ - ������
		glColor3f((icol - 510) / 255.0, 1.0, 0.0);
	}
	else if ((766 <= icol) && (icol <= 1020)) {

		// Ƹ���� - �������
		glColor3f(1.0, (255 - (icol - 765)) / 255.0, 0.0);
	};
}

#endif

// ����� ������������ � amg1r5.
integer myi_max(integer ia, integer ib)
{
	integer ir;
	if (ia < ib) ir = ib;
	else ir = ia;
	return ir;
} // max

integer myi_min(integer ia, integer ib)
{
	integer ir;
	if (ia < ib) ir = ia;
	else ir = ib;
	return ir;
} // min

#ifndef NO_OPENGL_GLFW

  // ������������� ���� ������ ��� ���������
  // �� ������������� �������������� ������� id ��� ��������� �������� �������,
  // � ����� �������� ��� ���� �������� ��������� � �������� �������.
void set_color_for_render(int id) {

	int icol = 0;
	if (iCFD_animation == 1) {
		icol = myi_min(1020, myi_max(0, round(1020 * ((animation_sequence_functions_openGL[iCURENT_ANIMATION_FUNCTION][iCURENT_ANIMATION_CADER][id] - minimum_val_for_render_pic) / (maximum_val_for_render_pic - minimum_val_for_render_pic)))));
	}
	else {
		icol = myi_min(1020, myi_max(0, round(1020 * ((potent_array_opengl[iCURENT_FUNC_openGL][id] - minimum_val_for_render_pic) / (maximum_val_for_render_pic - minimum_val_for_render_pic)))));
	}
	//if ((0 <= icol) && (icol <= 1020)) {
	//icol = mask_color_for_openGL[icol];// ��������� ����� ��������� ������.
	//}
	set_color_for_render_icolor(icol);
}

// ������������� ���� ������ ��� ���������
// �� ������������� ������������� ��������, � ����� �������� ��� ���� �������� ��������� � ��������.
void set_color_for_render(doublereal potent_val) {

	int icol = myi_min(1020, myi_max(0, round(1020 * ((potent_val - minimum_val_for_render_pic) / (maximum_val_for_render_pic - minimum_val_for_render_pic)))));
	//if ((0 <= icol) && (icol <= 1020)) {
	//icol = mask_color_for_openGL[icol];// ��������� ����� ��������� ������.
	//}
	set_color_for_render_icolor(icol);
}

// ��������� ������������ � ��������� ��������� �����. ����� ���������� ������ ����.
int DrawZbuffer(BLOCK*& b, integer lb, int maxelm, int maxnod, int**& nvtx, int* &whot_is_block);
// ��������� ������������ � ��������� ��������� �����. ����� ���������� ������ ������������ ������� potent.
int DrawZbufferColor(BLOCK*& b, integer lb, int maxelm, int maxnod, int**& nvtx, int*& whot_is_block, doublereal*& potent);
// ��������� ������������ � ��������� ��������� ����� � � ����������. ����� ���������� ������ ������������ ������� potent.
int DrawZbufferColorLight(BLOCK*& b, integer lb, int maxelm, int maxnod, int**& nvtx, int*& whot_is_block, doublereal*& potent);

#endif

// �����, ������ � ���������� ����������� ������ (��������) � ��������� � ���������� BiCGStab.
// �� ������������� � �������������. ����� ���������� � BiCGStab.
#include "my_LR.cpp" // ������������ �����

// 8 ������ 2016.
const bool bvery_big_memory = true; // true ��� ������ � ���� �� ������ � ����������� ������. ��� ����������� ������� �� ��������.

									// ������ ����������� �������� ������ ���������� � ����� MenterSST.cpp.
const doublereal K_limiter_min = 1.0e-20; //1.0e-14; // 1.0e-14 Fluent limits
										  // 0.1 ���� �������� � ���������� ����� ����� ��� k �� ��������� �� ������� 1.0.
										  // ��� �������� 1.0E-20 ������������. �������� k �� ������ 1e-14 �� 1e-10. ��������� ���� ��� � ���������� � ��������� ����������� 1 0.5�/�.
										  // ����������� ������ ����� ���������� ���������� ddemidov AMGCL.
										  // ������� ������� (����������) omega �� ���������� ���� �������� 1.5.
const doublereal Omega_limiter_min = 0.1; // 1.0; ����� ����� ������������. // 1.0e-20 Fluent limits
const doublereal Epsilon_limiter_min = 1.0e-20; // 1.0e-14; �������� ������� ���������... 1.0e-20 Fluent limits

UNION* my_union = nullptr; // ��� �����������.

						   // ���������� ����������
TEMPER my_global_temperature_struct;
integer flow_interior = 0; // ��������� ����� FLUID ���
FLOW* f = nullptr;



#include "my_linalg.cpp" // ���������� ������� �������� �������

// ������� �������� � ��������� tecplot360
#include "my_export_tecplot3.cpp"

#include "my_material_properties.cpp" // ���������� �������� ������� ����������


// ��� �������: 
// eqsolve_simple_gauss - ������ ���� ������� ���������� ������
// eqsolv_simple_holesskii - ������ ���� ������� ���������� ����������

// ������������� ����������� ��������� ���������-��������
// �� ����������� �����
#include "pamendment3.cpp"

#include "shortest_distance.cpp" // ���������� ����������� ���������� �� ������


// ���������� � ������� ���������� ������� ������������ �� 
// icepak user guide.
typedef struct TFLUENT_RESIDUAL {
	// ������ ������� ���������� �� ������ ����� ������� ����.
	// ������� ������������� � ���������� FLUENT
	// �.�. ����������� �� ������� fluent.
	doublereal res_vx; // ������� X ��������
	doublereal res_vy; // ������� Y ��������
	doublereal res_vz; // ������� Z ��������
	doublereal res_no_balance; // ������������������ ��������� �����.
	doublereal operating_value_b; // �������� ������������������ ���������� ����� � ���������� ��������.
	doublereal res_nusha; // ������� ���������������� ������������ ������������ ��������.
	doublereal res_turb_kinetik_energy; // ������� ������������ ������� ������������ ��������� � ������ SST K-Omega.
	doublereal res_turb_omega; // ������� �������� �������� ���������� ������������ ������� ������������ ��������� � ������ SST K-Omega.
	doublereal res_turb_kinetik_energy_std_ke; // ������� ������������ ������� ������������ ��������� � ����������� ������ k-epsilon
	doublereal res_turb_epsilon; // ������� �������� ���������� ������������ ������� ������������ ��������� � ����������� ������ k-epsilon
	doublereal res_turb_gamma_Langtry_Mentor; // �������������� ��� ������ �������� �������.
	doublereal res_turb_Re_Theta_Langtry_Mentor; // ����� ���������� ���  ������ �������� �������.

	TFLUENT_RESIDUAL() {
		// ������ ������� ���������� �� ������ ����� ������� ����.
		// ������� ������������� � ���������� FLUENT
		// �.�. ����������� �� ������� fluent.
		res_vx = 1.0; // ������� X ��������
		res_vy = 1.0; // ������� Y ��������
		res_vz = 1.0; // ������� Z ��������
		res_no_balance = 1.0; // ������������������ ��������� �����.
		operating_value_b = 1.0; // �������� ������������������ ���������� ����� � ���������� ��������.
		res_nusha = 1.0; // ������� ���������������� ������������ ������������ ��������.
		res_turb_kinetik_energy = 1.0; // ������� ������������ ������� ������������ ��������� � ������ SST K-Omega.
		res_turb_omega = 1.0; // ������� �������� �������� ���������� ������������ ������� ������������ ��������� � ������ SST K-Omega.
		res_turb_kinetik_energy_std_ke = 1.0; // ������� ������������ ������� ������������ ��������� � ����������� ������ k-epsilon
		res_turb_epsilon = 1.0; // ������� �������� ���������� ������������ ������� ������������ ��������� � ����������� ������ k-epsilon
		res_turb_gamma_Langtry_Mentor = 1.0; // �������������� ��� ������ �������� �������.
		res_turb_Re_Theta_Langtry_Mentor = 1.0; // ����� ���������� ���  ������ �������� �������.
	}
} FLUENT_RESIDUAL;


// ���������� ��������
#include "mysolverv0_03.cpp"


// �������������� ������ ��� �����������
// �� ������ ������������� �������,
// � ����� �������������� ������ ��� 
// ������������� �� ������ ������������� �������.
#include "my_unsteady_temperature.cpp"

// ������������� ��� ������������ ���������.
#include "my_nested_dissection.cpp"

#include <ctime> // ��� ������ ������� ����������.

integer ltdp = 0; // ���������� �������� �������� ��������� �� ����������� � �������� �����.
TEMP_DEP_POWER* gtdps = nullptr; // Garber temperature depend power sequence. 

								 // ���� ������ ����������:
integer lmatmax = 0; // ������������ ����� ����������
TPROP* matlist = nullptr; // ��������� ���� ������ ����������

void check_data(TEMPER t) {
	if (t.potent != nullptr) {
		for (integer i = 0; i < t.maxelm + t.maxbound; i++) {
			if (t.potent[i] != t.potent[i]) {
				std::cout << "t.potent[" << i << "] is " << t.potent[i] << "\n";
				system("pause");
			}
		}
	}
} // check_data

int main_body(BLOCK* &b, integer &lb, char ch_EXPORT_ALICE_ONLY = 'y') {
	//printLOGO();

	//system("PAUSE");	

	// ���������� ������, ���������� � ������, �������.
	//integer lb = 0;
	//BLOCK* b = nullptr;// ������ ������

	integer ls = 0, lw = 0, lu = 0;
	SOURCE* s = nullptr; // ������ ����������
	WALL* w = nullptr; // ������ ������ ������

					   // ��� ��� � ������ bFULL_AUTOMATIC ������� ������������ �������� �
					   // ������� ������������ �������, �� �������� ������� ����������� ���� ���� ���, �
					   // ��� ��������� ��������� ���� ��������� � ������ ���-�������.
					   // 20mm ���� ��������� � 1��� 9� �� 53� �� ���� ������ bFULL_AUTOMATIC.
					   // ���-������� ��� automatic
					   // ��������� ����������� ������ ��� ���-�������.
	shorter_hash_X = new doublereal[isize_shorter_hash];
	shorter_hash_Y = new doublereal[isize_shorter_hash];
	shorter_hash_Z = new doublereal[isize_shorter_hash];
	bshorter_hash_X = new bool[isize_shorter_hash];
	bshorter_hash_Y = new bool[isize_shorter_hash];
	bshorter_hash_Z = new bool[isize_shorter_hash];


	// �������������, ���������� ��.
	pfpir.fmaximum = 1.0e+30;
	pfpir.fminimum = -1.0e+30;
	pfpir.idir = LINE_DIRECTIONAL::Y_LINE_DIRECTIONAL;


	// 22.01.2017 �������������.
	eqin.fluidinfo = nullptr;
	my_global_temperature_struct.rootBT = nullptr;
	my_global_temperature_struct.rootSN = nullptr;
	my_global_temperature_struct.rootWE = nullptr;

	// 29 10 2016.
	// ������������� ����� ������ � ILU ������.
	milu_gl_buffer.alu_copy = nullptr;
	milu_gl_buffer.jlu_copy = nullptr;
	milu_gl_buffer.ju_copy = nullptr;

	my_amg_manager_init();

	// ����� �������.
	calculation_main_start_time_global_Depend = 0; // ������ ����� ��.
	unsigned int calculation_main_end_time = 0; // ��������� ����� ��.
	unsigned int calculation_main_seach_time = 0; // ����� ���������� ������� ���� � ��.

	calculation_main_start_time_global_Depend = clock(); // ������ ������ �����.

	bool bextendedprint = false; // ������ �� ��������� ����� ������������ �����.


								 //std::locale::global(std::locale("en_US.UTF-8"));
								 //system("mode con cols=166 lines=12000");
								 // � ������� ��� ���� ����� �������� ���� ������ � ������� ����� ������� � ������� ���
								 //HANDLE hOCol = GetStdHandle(STD_OUTPUT_HANDLE);
								 //SetConsoleTextAttribute(hOCol, FOREGROUND_GREEN);
								 // ��������� ��������� ������� (����� ���)	
								 //SetConsoleTextAttribute(hOCol, BACKGROUND_BLUE |
								 //	BACKGROUND_GREEN |
								 //	BACKGROUND_RED |
								 //	BACKGROUND_INTENSITY);

								 //system("cls");
								 // ���������� � ���, ��� ����� ��� ���������� ���� ��� ��� ��������� �����.

								 // ��� ��������.
								 // �� ������ ������ �����������. 
	doublereal* xpos = nullptr, *ypos = nullptr, *zpos = nullptr;
	doublereal* xposadd = nullptr, *yposadd = nullptr, *zposadd = nullptr;


	std::cout << "AliceFlow 3D x64 v0.48\n";
#ifdef _OPENMP 
	omp_set_num_threads(inumcore); // ��������� ����� �������
#endif




								   //ilu0_Saadtest();
								   //printf("the end Saad ilu0 test\n");
								   //system("PAUSE");

								   // ���������� ����� �� ������ �� ����.
								   //integer inx=120, iny=64, inz=64;
	integer inx = 30, iny = 30, inz = 30;
	integer inxadd = -1, inyadd = -1, inzadd = -1;
	doublereal dgx = 0.0, dgy = 0.0, dgz = 0.0; // ���� �������
	doublereal operatingtemperature = 20.0; // Operating Temperature 20.0 ����. �.

	lu = 0;
	// lu, my_union
	loadFromFile();
	premeshin("premeshin.txt", lmatmax, lb, ls, lw, matlist, b, s, w,
		dgx, dgy, dgz, inx, iny, inz, operatingtemperature, ltdp, gtdps, lu, my_union);
	freeStringList();



	// ��������� ���� �� ����� �� ������� ��������
	// ����� ������, ������ � ���������� �����. 02.08.2019.
	// � ����� ����������. 23.07.2020
	BODY_CHECK(b, lb, w, lw, s, ls, my_union, lu);



	init_QSBid(lb, b, w, lw, s, ls); // ��� ���������� ������ ������� myisblock_id.



	if ((steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::CFD_STEADY) ||
		(steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::CFD_UNSTEADY)) {
		// ��� ������� ��������� ������������� �� ������� ������ load.txt ����.
		remove("load.txt");
	}

	if (1 && steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::PREOBRAZOVATEL_FOR_REPORT) {
		// �������������� ����� � ������������ ����������.
		// ��� ��������� �������. 05.01.2018.
		tecplot360patcher_for_print_in_report();
		exit(1);
	}

	bool bCAD = true;
	for (int i60 = 1; i60 < lb; i60++) {
		if (b[i60].g.itypegeom != CAD_STL) {
			if (b[i60].itype != PHYSICS_TYPE_IN_BODY::HOLLOW)
			{
				bCAD = false;
			}
		}
	}

	integer iCabinetMarker = 0;
	if (bCAD) {
		// ��� CAD ��������� ������ ����� ����������� ����� � ������ ��������� �������.

		int isize = inx;
		xpos = new doublereal[isize + 1];
		xpos[0] = b[0].g.xS;
		doublereal dstep = (b[0].g.xE - b[0].g.xS) / (1.0* isize);
		for (int i60 = 1; i60 <= isize; i60++) {
			xpos[i60] = b[0].g.xS + dstep * i60;
		}
		xpos[isize] = b[0].g.xE;

		isize = iny;
		ypos = new doublereal[isize + 1];
		ypos[0] = b[0].g.yS;
		dstep = (b[0].g.yE - b[0].g.yS) / (1.0 * isize);
		for (int i60 = 1; i60 <= isize; i60++) {
			ypos[i60] = b[0].g.yS + dstep * i60;
		}
		ypos[isize] = b[0].g.yE;

		isize = inz;
		zpos = new doublereal[isize + 1];
		zpos[0] = b[0].g.zS;
		dstep = (b[0].g.zE - b[0].g.zS) / (1.0 * isize);
		for (int i60 = 1; i60 <= isize; i60++) {
			zpos[i60] = b[0].g.zS + dstep * i60;
		}
		zpos[isize] = b[0].g.zE;
	}
	else {
		if (iswitchMeshGenerator == CONFORMAL_MESH_GENERATOR_SELECTOR::SIMPLEMESHGEN_MESHER) {
			simplemeshgen(xpos, ypos, zpos, inx, iny, inz, lb, ls, lw, b, s, w, lu, my_union, matlist,
				xposadd, yposadd, zposadd, inxadd, inyadd, inzadd, iCabinetMarker);
		}
		else if (iswitchMeshGenerator == CONFORMAL_MESH_GENERATOR_SELECTOR::UNEVENSIMPLEMESHGEN_MESHER) {
			unevensimplemeshgen(xpos, ypos, zpos, inx, iny, inz, lb, ls, lw, b, s, w, lu, my_union, matlist,
				dgx, dgy, dgz, xposadd, yposadd, zposadd, inxadd, inyadd, inzadd, iCabinetMarker);
			// ��������� ������������� ����� � �������������� ������������� ���������������.
		}
		else if (iswitchMeshGenerator == CONFORMAL_MESH_GENERATOR_SELECTOR::COARSEMESHGEN_MESHER) {
			// � ��������� ������� coarse Mesh ��� � Icepak.
			// �������� ������ ����������� ��������������-���������������, � 
			// ������� ������������� ���������� ����������� �����, �.�. CPU �������� � 4��� �
			// ���������� ����������������� ��������� ������� ������� ��������.
			coarsemeshgen(xpos, ypos, zpos, inx, iny, inz, lb, ls, lw, b, s, w, lu, my_union, matlist,
				xposadd, yposadd, zposadd, inxadd, inyadd, inzadd, iCabinetMarker);
		}
		else {
			switch (iswitchMeshGenerator) {
			case CONFORMAL_MESH_GENERATOR_SELECTOR::COARSEMESHGEN_MESHER:
				std::cout << "your mesh generator is undefined " << "CONFORMAL_MESH_GENERATOR_SELECTOR::COARSEMESHGEN_MESHER" << "\n";
				break;
			case CONFORMAL_MESH_GENERATOR_SELECTOR::UNEVENSIMPLEMESHGEN_MESHER:
				std::cout << "your mesh generator is undefined " << "CONFORMAL_MESH_GENERATOR_SELECTOR::UNEVENSIMPLEMESHGEN_MESHER" << "\n";
				break;
			case CONFORMAL_MESH_GENERATOR_SELECTOR::SIMPLEMESHGEN_MESHER:
				std::cout << "your mesh generator is undefined " << "CONFORMAL_MESH_GENERATOR_SELECTOR::SIMPLEMESHGEN_MESHER" << "\n";
				break;
			default:
				std::cout << "error: your mesh generator is undefined " << "\n";
				break;
			}


			system("pause");
			exit(1);
		}




		// ������ ��������� ����� � ������������.
		for (integer iu = 0; iu < lu; iu++) {
			my_union[iu].inxadd = -1;
			my_union[iu].inyadd = -1;
			my_union[iu].inzadd = -1;
			my_union[iu].xposadd = nullptr;
			my_union[iu].yposadd = nullptr;
			my_union[iu].zposadd = nullptr;
			my_union[iu].xpos = nullptr;
			my_union[iu].ypos = nullptr;
			my_union[iu].zpos = nullptr;
			integer iup1 = iu + 1;
			switch (my_union[iu].iswitchMeshGenerator) {
			case CONFORMAL_MESH_GENERATOR_SELECTOR::SIMPLEMESHGEN_MESHER: simplemeshgen(my_union[iu].xpos, my_union[iu].ypos, my_union[iu].zpos,
				my_union[iu].inx, my_union[iu].iny, my_union[iu].inz, lb, ls, lw, b, s, w, lu, my_union, matlist,
				my_union[iu].xposadd, my_union[iu].yposadd, my_union[iu].zposadd, my_union[iu].inxadd,
				my_union[iu].inyadd, my_union[iu].inzadd, iup1);
				break;
			case CONFORMAL_MESH_GENERATOR_SELECTOR::UNEVENSIMPLEMESHGEN_MESHER: unevensimplemeshgen(my_union[iu].xpos, my_union[iu].ypos, my_union[iu].zpos, my_union[iu].inx,
				my_union[iu].iny, my_union[iu].inz, lb, ls, lw, b, s, w, lu, my_union, matlist,
				dgx, dgy, dgz, my_union[iu].xposadd, my_union[iu].yposadd, my_union[iu].zposadd,
				my_union[iu].inxadd, my_union[iu].inyadd, my_union[iu].inzadd, iup1);
				// ��������� ������������� ����� � �������������� ������������� ���������������.
				break;
			case CONFORMAL_MESH_GENERATOR_SELECTOR::COARSEMESHGEN_MESHER: coarsemeshgen(my_union[iu].xpos, my_union[iu].ypos, my_union[iu].zpos,
				my_union[iu].inx, my_union[iu].iny, my_union[iu].inz, lb, ls, lw, b, s, w, lu, my_union, matlist,
				my_union[iu].xposadd, my_union[iu].yposadd, my_union[iu].zposadd, my_union[iu].inxadd,
				my_union[iu].inyadd, my_union[iu].inzadd, iup1);
				break;
			default:
				coarsemeshgen(my_union[iu].xpos, my_union[iu].ypos, my_union[iu].zpos,
					my_union[iu].inx, my_union[iu].iny, my_union[iu].inz, lb, ls, lw, b, s, w, lu, my_union, matlist,
					my_union[iu].xposadd, my_union[iu].yposadd, my_union[iu].zposadd, my_union[iu].inxadd,
					my_union[iu].inyadd, my_union[iu].inzadd, iup1);
				break;
			}
			// ��������� ������������ ������ addboundary_rudiment � �� addboundary, ��� ��� � �������� ����� 
			// �������� �� ����������� ��� ���������� ����� ��������. 14.03.2019.
			// ����������� ��������� ����� ���������� �������� ����� ����� �������. 
			// ��������� �������� ����� ����������� �������� ��� ��������� �������� �������������.
			for (integer i76 = 0; i76 <= inx; i76++) {
				// ��������� ���������� �������� ����� ��������.
				if ((xpos[i76] >= my_union[iu].xS) && (xpos[i76] <= my_union[iu].xE)) {
					addboundary_rudiment(my_union[iu].xpos, my_union[iu].inx, xpos[i76], YZ_PLANE, b, lb, w, lw, s, ls);
				}
			}
			Sort_method(my_union[iu].xpos, my_union[iu].inx);
			for (integer i76 = 0; i76 <= iny; i76++) {
				// ��������� ���������� �������� ����� ��������.
				if ((ypos[i76] >= my_union[iu].yS) && (ypos[i76] <= my_union[iu].yE)) {
					addboundary_rudiment(my_union[iu].ypos, my_union[iu].iny, ypos[i76], XZ_PLANE, b, lb, w, lw, s, ls);
				}
			}
			Sort_method(my_union[iu].ypos, my_union[iu].iny);
			for (integer i76 = 0; i76 <= inz; i76++) {
				// ��������� ���������� �������� ����� ��������.
				if ((zpos[i76] >= my_union[iu].zS) && (zpos[i76] <= my_union[iu].zE)) {
					addboundary_rudiment(my_union[iu].zpos, my_union[iu].inz, zpos[i76], XY_PLANE, b, lb, w, lw, s, ls);
				}
			}
			Sort_method(my_union[iu].zpos, my_union[iu].inz);
		}


		// � ����� add_line.txt ��� ��� ������� ���������� �������� ����� �������������� 
		// �������� ������������� ������� ��� ������.
		FILE* fp_add_line = NULL;

		// �������� ����� ��� ������.
#ifdef MINGW_COMPILLER
		fp_add_line = fopen64("add_line.txt", "r");
		int err_add_line = 0;
		if (fp_add_line != NULL) {
			err_add_line = 0;
		}
		else {
			err_add_line = 1; // ������ ��������.
		}
#else
		errno_t err_add_line;
		err_add_line = fopen_s(&fp_add_line, "add_line.txt", "r");
#endif

		if (fp_add_line != NULL) {
			if (err_add_line != 0) {
				//printf("Open File add_line Error\n");
				//getchar();
				//exit(0);
				// return bfound;
			}
			else
			{

				//printf("incomming");
				//getchar();
				// ��� ��������� �������� ������ ����������� ������������� �������.

				int ix_add0 = 0;
				int iy_add0 = 0;
				int iz_add0 = 0;
				fscanf_s(fp_add_line, "%d", &ix_add0);

				for (int i0c = 0; i0c < ix_add0; i0c++) {
					float fin0 = 0.0;
					fscanf_s(fp_add_line, "%f", &fin0);

					// ��������� ���������� �������� ����� ��������.
					if ((fin0 >= b[0].g.xS) && (fin0 <= b[0].g.xE)) {
						addboundary_rudiment(xpos, inx, fin0, YZ_PLANE, b, lb, w, lw, s, ls);
					}
					Sort_method(xpos, inx);

				}
				Sort_method(xpos, inx);


				fscanf_s(fp_add_line, "%d", &iy_add0);

				for (int i0c = 0; i0c < iy_add0; i0c++) {
					float fin0 = 0.0;
					fscanf_s(fp_add_line, "%f", &fin0);

					// ��������� ���������� �������� ����� ��������.
					if ((fin0 >= b[0].g.yS) && (fin0 <= b[0].g.yE)) {
						addboundary_rudiment(ypos, iny, fin0, XZ_PLANE, b, lb, w, lw, s, ls);
					}
					Sort_method(ypos, iny);

				}
				Sort_method(ypos, iny);

				fscanf_s(fp_add_line, "%d", &iz_add0);

				for (int i0c = 0; i0c < iz_add0; i0c++) {
					float fin0 = 0.0;
					fscanf_s(fp_add_line, "%f", &fin0);

					// ��������� ���������� �������� ����� ��������.
					if ((fin0 >= b[0].g.zS) && (fin0 <= b[0].g.zE)) {
						addboundary_rudiment(zpos, inz, fin0, XY_PLANE, b, lb, w, lw, s, ls);
					}
					Sort_method(zpos, inz);

				}
				Sort_method(zpos, inz);

				fclose(fp_add_line);
			}

		}


	}

	if (b_on_adaptive_local_refinement_mesh) {
		// ������ 1TT113 22.03.2018:
		// AliceCorse 1266911 ����, � ����������������� ���������� 93592. 
		// 40��� 33� ������ ���������� �������������� �� 10000�.
		// ����� ���������� ����� 2��� 6�. 14 ������� ��������.
		// AliceFine 1092571 ����, � ����������������� ���������� 8050.
		// AliceFine (��������� �����) 1335076 ����, � ����������������� ���������� 93448.
		// ����� ���������� ����� 4��� 25�. 14 ������� ��������.
		// ����� ���������� ����� (��������� �����) 6��� 32�. 14 ������� ��������.



		std::cout << "starting ALICE\n";
		if (0) {
			if (TYPE_ALICE_MESH::MULTI_PASS_MEDIUM_ALICE_MESH == itype_ALICE_Mesh) {
				// ��� ������ �� � ���� ������ ������ �� ������� �������� ����������� ������.
				doublereal* xpos_copy = nullptr;
				// 10 ������� ������� �������� ���������.
				const integer jiterM = my_amg_manager.nu1_Temperature;
				// ������������� ��������� ������� ��������� �������.
				for (integer jiter = 0; jiter < jiterM; jiter++) {
					xpos_copy = new doublereal[2 * (inx + 1) - 1];
#pragma omp parallel for
					for (integer i74 = 0; i74 < inx; i74++) {
						xpos_copy[2 * i74] = xpos[i74];
						xpos_copy[2 * i74 + 1] = 0.5 * (xpos[i74] + xpos[i74 + 1]);
					}
					xpos_copy[2 * (inx + 1) - 2] = xpos[inx];
					delete[] xpos;
					xpos = nullptr;
					xpos = new doublereal[2 * (inx + 1) - 1];
#pragma omp parallel for
					for (integer i74 = 0; i74 < 2 * (inx + 1) - 1; i74++) {
						xpos[i74] = xpos_copy[i74];
					}
					delete[] xpos_copy;
					xpos_copy = nullptr;
					inx = 2 * (inx + 1) - 2;
				}

				for (integer jiter = 0; jiter < jiterM; jiter++) {
					xpos_copy = new doublereal[2 * (iny + 1) - 1];
#pragma omp parallel for
					for (integer i74 = 0; i74 < iny; i74++) {
						xpos_copy[2 * i74] = ypos[i74];
						xpos_copy[2 * i74 + 1] = 0.5 * (ypos[i74] + ypos[i74 + 1]);
					}
					xpos_copy[2 * (iny + 1) - 2] = ypos[iny];
					delete[] ypos;
					ypos = nullptr;
					ypos = new doublereal[2 * (iny + 1) - 1];
#pragma omp parallel for
					for (integer i74 = 0; i74 < 2 * (iny + 1) - 1; i74++) {
						ypos[i74] = xpos_copy[i74];
					}
					delete[] xpos_copy;
					xpos_copy = nullptr;
					iny = 2 * (iny + 1) - 2;
				}

				for (integer jiter = 0; jiter < jiterM; jiter++) {
					xpos_copy = new doublereal[2 * (inz + 1) - 1];
#pragma omp parallel for
					for (integer i74 = 0; i74 < inz; i74++) {
						xpos_copy[2 * i74] = zpos[i74];
						xpos_copy[2 * i74 + 1] = 0.5 * (zpos[i74] + zpos[i74 + 1]);
					}
					xpos_copy[2 * (inz + 1) - 2] = zpos[inz];
					delete[] zpos;
					zpos = nullptr;
					zpos = new doublereal[2 * (inz + 1) - 1];
#pragma omp parallel for
					for (integer i74 = 0; i74 < 2 * (inz + 1) - 1; i74++) {
						zpos[i74] = xpos_copy[i74];
					}
					delete[] xpos_copy;
					xpos_copy = nullptr;
					inz = 2 * (inz + 1) - 2;
				}
			}
		}

		// ��� ������� ������� �������� �� ������ ������������ �������,
		// � ������ �������� �� ��������� ���� � ��������� new ��� malloc.
		integer maxelm_loc = (inx + 1) * (iny + 1) * (inz + 1);
		bool bOkal = alice_mesh(xpos, ypos, zpos, inx, iny, inz, b, lb, lw, w, s, ls, maxelm_loc, xposadd, yposadd, zposadd, inxadd, inyadd, inzadd);
		//system("PAUSE");

		

		if (0 || itype_ALICE_Mesh == TYPE_ALICE_MESH::MULTI_PASS_MEDIUM_ALICE_MESH/*1*/) {
			// �������� ��������� ���������.

			/*
			����� ���� ������� ����������� �� ��������� �������� �����
			��� ������� ������������������� �������� ����� ������������
			�� �� ����� ������������� ���� ��������� AliceMedium �����
			� ��� ��������� ���������� � ����� �������� ����� �� ����������.
			��� ���� ����� �������� ����� �� ��������� ���������� �����
			� ��� ������ ����� ����� ��������� ��������.
			17.08.2019
			*/
			doublereal dSTOP_flag1 = 1.0e+4; // ����� ���������� �����.
			doublereal dSTOP_flag2 = 1.0e+1; // ����� ���������� �����.

			while (!bOkal) {
				std::cout << "repeat call ALICE...\n";
				//system("PAUSE");

				/* 3.09.2017
				���� ����� ����� ������� ����� ��� ��������� ������ ��������������� �������������� ����� �������
				����� �������. ������ ������� ������� ������������� ������ � ���������� ���� �� ����� ����������� �������������
				����� ���� ������ ���� ����� ������ ���������. � ���������� �������� ������������� ����������� ��������. ������ �����
				��������� ������� ������ �������������� ����� ����� �� ������� ����� ����� ��� ������������ ����� (�.�. ����� ����������
				��� ���������� ����� �������� ����� ������� ����� ������ �������� ����� ���������� �� ����� ��� �� 1. ������� ������
				����� ����������� ��� ��������� ������ ������ � ������� if_disbalnce(...) � ��� ������ ����� ������ ����������� �������
				�������� ��������������� ����� �������� �����. �������� ����������� ������ � ������������ �� �������� ������� � ����������
				���� ����� ���������� ������ ������ ������ ������� ����� ��� �������� ����������� �������� �����.
				*/


				// ����� ���������� ������ �� ��� octree ������ � ����������� �����.
				std::cout << "free octree start...\n";
				//system("PAUSE");
				//system("PAUSE");
				free_octree(oc_global, maxelm_loc);
				delete[] my_ALICE_STACK;
				top_ALICE_STACK = 0;
				std::cout << "free octree end...\n";
				doublereal t_1 = NormaV(xpos, inx + 1);
				doublereal t_2 = NormaV(ypos, iny + 1);
				doublereal t_3 = NormaV(zpos, inz + 1);
				dSTOP_flag2 = sqrt(t_1 * t_1 + t_2 * t_2 + t_3 * t_3);
				std::cout << "comparison mesh " << fabs(dSTOP_flag2 - dSTOP_flag1) << "\n";
				//system("pause");

				// ����� ���������� ��������� �����.
				delete[] xpos;
				xpos = nullptr;
				inx = 0;
				delete[] ypos;
				ypos = nullptr;
				iny = 0;
				delete[] zpos;
				zpos = nullptr;
				inz = 0;

				std::cout << "free xpos, ypos, zpos\n";
				//system("PAUSE");

				integer iCabinetMarker = 0;
				if (iswitchMeshGenerator == CONFORMAL_MESH_GENERATOR_SELECTOR::SIMPLEMESHGEN_MESHER) {
					simplemeshgen(xpos, ypos, zpos, inx, iny, inz, lb, ls, lw, b, s, w, lu, my_union,
						matlist, xposadd, yposadd, zposadd, inxadd, inyadd, inzadd, iCabinetMarker);
				}
				else if (iswitchMeshGenerator == CONFORMAL_MESH_GENERATOR_SELECTOR::UNEVENSIMPLEMESHGEN_MESHER) {
					unevensimplemeshgen(xpos, ypos, zpos, inx, iny, inz, lb, ls, lw, b, s, w, lu, my_union,
						matlist, dgx, dgy, dgz, xposadd, yposadd, zposadd, inxadd, inyadd, inzadd, iCabinetMarker); // ��������� ������������� ����� � �������������� ������������� ���������������.
				}
				else if (iswitchMeshGenerator == CONFORMAL_MESH_GENERATOR_SELECTOR::COARSEMESHGEN_MESHER) {
					// � ��������� ������� coarse Mesh ��� � ANSYS Icepak.
					// �������� ������ ����������� ��������������-���������������, � 
					// ������� ������������� ���������� ����������� �����, �.�. CPU �������� � 4��� �
					// ���������� ����������������� ��������� ������� ������� ��������.
					coarsemeshgen(xpos, ypos, zpos, inx, iny, inz, lb, ls, lw, b, s, w, lu, my_union,
						matlist, xposadd, yposadd, zposadd, inxadd, inyadd, inzadd, iCabinetMarker);
				}
				else {

					printf("error: your mesh generator is undefined %d\n", iswitchMeshGenerator);

					system("pause");
					exit(1);
				}
				// ����� ������ ���� �����������.

				printf("new construct xpos, ypos, zpos\n");
				//system("PAUSE");

				bOkal = alice_mesh(xpos, ypos, zpos, inx, iny, inz, b, lb, lw, w, s, ls, maxelm_loc, xposadd, yposadd, zposadd, inxadd, inyadd, inzadd);

				if (fabs(dSTOP_flag2 - dSTOP_flag1) < 1.0e-20) {
					bOkal = true;
					break;
				}
				dSTOP_flag1 = dSTOP_flag2;

				//system("PAUSE");
			}
		}
		printf("end ALICE\n");
	}

	iCabinetMarker = 0;
	load_TEMPER_and_FLOW(my_global_temperature_struct, f, inx, iny, inz, xpos, ypos, zpos, flow_interior,
		b, lb, lw, w, s, ls, lu, my_union, operatingtemperature, matlist, bextendedprint,
		dgx, dgy, dgz, b_on_adaptive_local_refinement_mesh, false, iCabinetMarker);

	my_global_temperature_struct.operatingtemperature = operatingtemperature;


	for (integer iu = 0; iu < lu; iu++) {
		integer iup1 = iu + 1;
		load_TEMPER_and_FLOW(my_union[iu].t, my_union[iu].f,
			my_union[iu].inx, my_union[iu].iny, my_union[iu].inz,
			my_union[iu].xpos, my_union[iu].ypos, my_union[iu].zpos,
			my_union[iu].flow_interior,
			b, lb, lw, w, s, ls, lu, my_union, operatingtemperature, matlist, bextendedprint,
			dgx, dgy, dgz, b_on_adaptive_local_refinement_mesh, false, iup1);

		my_union[iu].t.operatingtemperature = operatingtemperature;
	}




	// ��� ����� ������ ����� ��� ������� �������������� ������.
	my_global_temperature_struct.inx_copy = inx;
	my_global_temperature_struct.iny_copy = iny;
	my_global_temperature_struct.inz_copy = inz;
	my_global_temperature_struct.operatingtemperature_copy = operatingtemperature;
	my_global_temperature_struct.xpos_copy = new doublereal[inx + 1];
	my_global_temperature_struct.ypos_copy = new doublereal[iny + 1];
	my_global_temperature_struct.zpos_copy = new doublereal[inz + 1];
	// ������ ���������� ����� ��� �������� ����������� ������,
	// ��������� ������ ����� ��������� �� ��� � ����� ������������� 
	// ���� ����������.
#pragma omp parallel for
	for (integer i_7 = 0; i_7 < inx + 1; i_7++) {
		my_global_temperature_struct.xpos_copy[i_7] = xpos[i_7];
	}
#pragma omp parallel for
	for (integer i_7 = 0; i_7 < iny + 1; i_7++) {
		my_global_temperature_struct.ypos_copy[i_7] = ypos[i_7];
	}
#pragma omp parallel for
	for (integer i_7 = 0; i_7 < inz + 1; i_7++) {
		my_global_temperature_struct.zpos_copy[i_7] = zpos[i_7];
	}

	for (integer iu = 0; iu < lu; iu++) {
		// ��� ����� ������ ����� ��� ������� �������������� ������.
		my_union[iu].t.inx_copy = my_union[iu].inx;
		my_union[iu].t.iny_copy = my_union[iu].iny;
		my_union[iu].t.inz_copy = my_union[iu].inz;
		my_union[iu].t.operatingtemperature_copy = operatingtemperature;
		my_union[iu].t.xpos_copy = new doublereal[my_union[iu].inx + 1];
		my_union[iu].t.ypos_copy = new doublereal[my_union[iu].iny + 1];
		my_union[iu].t.zpos_copy = new doublereal[my_union[iu].inz + 1];
		// ������ ���������� ����� ��� �������� ����������� ������,
		// ��������� ������ ����� ��������� �� ��� � ����� ������������� 
		// ���� ����������.
#pragma omp parallel for
		for (integer i_7 = 0; i_7 < my_union[iu].inx + 1; i_7++) {
			my_union[iu].t.xpos_copy[i_7] = my_union[iu].xpos[i_7];
		}
#pragma omp parallel for
		for (integer i_7 = 0; i_7 < my_union[iu].iny + 1; i_7++) {
			my_union[iu].t.ypos_copy[i_7] = my_union[iu].ypos[i_7];
		}
#pragma omp parallel for
		for (integer i_7 = 0; i_7 < my_union[iu].inz + 1; i_7++) {
			my_union[iu].t.zpos_copy[i_7] = my_union[iu].zpos[i_7];
		}
	}

	// ������������ ����������� ������ �� ��� octree ������.
	if (b_on_adaptive_local_refinement_mesh) {
		printf("free octree start...\n");
		//system("PAUSE");
		//system("PAUSE");
		integer maxelm_loc = (inx + 1) * (iny + 1) * (inz + 1);
		free_octree(oc_global, maxelm_loc);
		delete[] my_ALICE_STACK;
		top_ALICE_STACK = 0;
		printf("free octree end...\n");
		//system("PAUSE");
		//system("PAUSE");
	}

	if (0) {
		xyplot(f, flow_interior, my_global_temperature_struct);
		printf("after load temper and flow. OK.\n");
		//system("PAUSE"); // debug avtosave
		system("pause");
	}

	// �� ���� ����� ���������� ����� �� �������� ���������������.
	if (!b_on_adaptive_local_refinement_mesh) {
		// ������ ���������� � �������� ������ ��� ������������� ������ LR:
		if (2 == iswitchsolveramg_vs_BiCGstab_plus_ILU2) {
			// Lr1sk algorithm
			constr_line(f, flow_interior);  // ��� �������������
		}
		my_global_temperature_struct.rootBT = nullptr;
		my_global_temperature_struct.rootSN = nullptr;
		my_global_temperature_struct.rootWE = nullptr;
		if (2 == iswitchsolveramg_vs_BiCGstab_plus_ILU2) {
			// Lr1sk algorithm
			constr_line_temp(my_global_temperature_struct, b, lb); // ��� ����������������
			printf("LR preprocessing finish...\n");
		}
	}

	// ���������� ������ ��� ��������������� �������������� ������.

	amgGM.a = nullptr;
	amgGM.f = nullptr;
	amgGM.ia = nullptr;
	amgGM.ig = nullptr;
	amgGM.ja = nullptr;
	amgGM.u = nullptr;
	amgGM.nda = -1;
	amgGM.ndf = -1;
	amgGM.ndia = -1;
	amgGM.ndig = -1;
	amgGM.ndja = -1;
	amgGM.ndu = -1;


	//PARDATA nd;
	nd.ncore = 2; // ��� ����.
				  // �� ��������� ��� ��������� ���������.
	nd.b0.active = false;
	nd.b00.active = false;
	nd.b01.active = false;
	nd.b000.active = false;
	nd.b001.active = false;
	nd.b010.active = false;
	nd.b011.active = false;
	if (0 && (1 == flow_interior)) {
		calc_front(f, f[0], my_global_temperature_struct, flow_interior, ls, lw, w, nd, b, lb, s);
		// ���������� ��������� !
		printf("separator compleate...\n");
		//system("PAUSE");
	}



	my_global_temperature_struct.free_temper_level1 = false; // ������ ���������������� ������������ ������ ����������� ��� ������ ������� ����� �������� ������.
	my_global_temperature_struct.free_temper_level2 = false; // ������������ ������ ��� �������� ������� ��� ���������� � � SIMPLESPARSE ������.	

	printf("construction of all structures...\n");
	printf("mesh check start...\n");
	const doublereal d_zalipanie = 1.0e-23;// ������
#if doubleintprecision == 1
	doublereal minimum_gap = 1.0e60;
	for (integer i = 0; i < inx; i++) {
		if ((xpos[i + 1] - xpos[i]) < minimum_gap) minimum_gap = (xpos[i + 1] - xpos[i]);
		if ((xpos[i + 1] - xpos[i]) < d_zalipanie) {
			//printf("error: zalipanie po X: xpos[%lld]=%e xpos[%lld]=%e inx=%lld\n", i, xpos[i], i + 1, xpos[i + 1], inx);
			std::cout << "error: zalipanie po X: xpos[" << i << "]=" << xpos[i] << " xpos[" << i + 1 << "]=" << xpos[i + 1] << " inx=" << inx << std::endl;
			std::cout << "tolerance zalipanie=" << d_zalipanie << std::endl;
			system("pause");
		}
	}
	std::cout << "minimum gap X=" << minimum_gap << std::endl;
	minimum_gap = 1.0e60;
	for (integer i = 0; i < iny; i++) {
		if ((ypos[i + 1] - ypos[i]) < minimum_gap) minimum_gap = (ypos[i + 1] - ypos[i]);
		if ((ypos[i + 1] - ypos[i]) < d_zalipanie) {
			//printf("error: zalipanie po Y: ypos[%lld]=%e ypos[%lld]=%e iny=%lld\n", i, ypos[i], i + 1, ypos[i + 1], iny);
			std::cout << "error: zalipanie po Y: ypos[" << i << "]=" << ypos[i] << " ypos[" << i + 1 << "]=" << ypos[i + 1] << " iny=" << iny << std::endl;
			std::cout << "tolerance zalipanie=" << d_zalipanie << std::endl;
			system("pause");
		}
	}
	std::cout << "minimum gap Y=" << minimum_gap << std::endl;
	minimum_gap = 1.0e60;
	for (integer i = 0; i < inz; i++) {
		if ((zpos[i + 1] - zpos[i]) < minimum_gap) minimum_gap = (zpos[i + 1] - zpos[i]);
		if ((zpos[i + 1] - zpos[i]) < d_zalipanie) {
			//printf("error: zalipanie po Z: zpos[%lld]=%e zpos[%lld]=%e inz=%lld\n", i, zpos[i], i + 1, zpos[i + 1], inz);
			std::cout << "error: zalipanie po Z: zpos[" << i << "]=" << zpos[i] << " zpos[" << i + 1 << "]=" << zpos[i + 1] << " inz=" << inz << std::endl;
			std::cout << "tolerance zalipanie=" << d_zalipanie << std::endl;
			system("pause");
		}
	}
	std::cout << "minimum gap Z=" << minimum_gap << std::endl;
	for (integer iP = 0; iP < my_global_temperature_struct.maxelm; iP++) {
		if ((my_global_temperature_struct.nvtx[0][iP] == 0) || (my_global_temperature_struct.nvtx[1][iP] == 0) || (my_global_temperature_struct.nvtx[2][iP] == 0) || (my_global_temperature_struct.nvtx[3][iP] == 0) || (my_global_temperature_struct.nvtx[4][iP] == 0) || (my_global_temperature_struct.nvtx[5][iP] == 0) || (my_global_temperature_struct.nvtx[6][iP] == 0) || (my_global_temperature_struct.nvtx[7][iP] == 0)) {
			printf("nvtx[%lld]: %lld %lld %lld %lld %lld %lld %lld %lld \n", iP, my_global_temperature_struct.nvtx[0][iP] - 1, my_global_temperature_struct.nvtx[1][iP] - 1, my_global_temperature_struct.nvtx[2][iP] - 1, my_global_temperature_struct.nvtx[3][iP] - 1, my_global_temperature_struct.nvtx[4][iP] - 1, my_global_temperature_struct.nvtx[5][iP] - 1, my_global_temperature_struct.nvtx[6][iP] - 1, my_global_temperature_struct.nvtx[7][iP] - 1);
		}
	}
#else
	doublereal minimum_gap = 1.0e60;
	for (integer i = 0; i < inx; i++) {
		if ((xpos[i + 1] - xpos[i]) < minimum_gap) minimum_gap = (xpos[i + 1] - xpos[i]);
		if ((xpos[i + 1] - xpos[i]) < d_zalipanie) {
			//printf("error: zalipanie po X: xpos[%d]=%e xpos[%d]=%e inx=%d\n", i, xpos[i], i + 1, xpos[i + 1], inx);
			std::cout << "error: zalipanie po X: xpos[" << i << "]=" << xpos[i] << " xpos[" << i + 1 << "]=" << xpos[i + 1] << " inx=" << inx << std::endl;
			std::cout << "tolerance zalipanie=" << d_zalipanie << std::endl;
			system("pause");
		}
	}
	std::cout << "minimum gap X=" << minimum_gap << std::endl;
	minimum_gap = 1.0e60;
	for (integer i = 0; i < iny; i++) {
		if ((ypos[i + 1] - ypos[i]) < minimum_gap) minimum_gap = (ypos[i + 1] - ypos[i]);
		if ((ypos[i + 1] - ypos[i]) < d_zalipanie) {
			//	printf("error: zalipanie po Y: ypos[%d]=%e ypos[%d]=%e iny=%d\n", i, ypos[i], i + 1, ypos[i + 1], iny);
			std::cout << "error: zalipanie po Y: ypos[" << i << "]=" << ypos[i] << " ypos[" << i + 1 << "]=" << ypos[i + 1] << " iny=" << iny << std::endl;
			std::cout << "tolerance zalipanie=" << d_zalipanie << std::endl;
			system("pause");
		}
	}
	std::cout << "minimum gap Y=" << minimum_gap << std::endl;
	minimum_gap = 1.0e60;
	for (integer i = 0; i < inz; i++) {
		if ((zpos[i + 1] - zpos[i]) < minimum_gap) minimum_gap = (zpos[i + 1] - zpos[i]);
		if ((zpos[i + 1] - zpos[i]) < d_zalipanie) {
			//		printf("error: zalipanie po Z: zpos[%d]=%e zpos[%d]=%e inz=%d\n", i, zpos[i], i + 1, zpos[i + 1], inz);
			std::cout << "error: zalipanie po Z: zpos[" << i << "]=" << zpos[i] << " zpos[" << i + 1 << "]=" << zpos[i + 1] << " inz=" << inz << std::endl;
			std::cout << "tolerance zalipanie=" << d_zalipanie << std::endl;
			system("pause");
		}
	}
	std::cout << "minimum gap Z=" << minimum_gap << std::endl;
	for (integer iP = 0; iP < my_global_temperature_struct.maxelm; iP++) {
		if ((my_global_temperature_struct.nvtx[0][iP] == 0) || (my_global_temperature_struct.nvtx[1][iP] == 0) || (my_global_temperature_struct.nvtx[2][iP] == 0) || (my_global_temperature_struct.nvtx[3][iP] == 0) || (my_global_temperature_struct.nvtx[4][iP] == 0) || (my_global_temperature_struct.nvtx[5][iP] == 0) || (my_global_temperature_struct.nvtx[6][iP] == 0) || (my_global_temperature_struct.nvtx[7][iP] == 0)) {
			printf("nvtx[%d]: %d %d %d %d %d %d %d %d \n", iP, my_global_temperature_struct.nvtx[0][iP] - 1, my_global_temperature_struct.nvtx[1][iP] - 1, my_global_temperature_struct.nvtx[2][iP] - 1, my_global_temperature_struct.nvtx[3][iP] - 1, my_global_temperature_struct.nvtx[4][iP] - 1, my_global_temperature_struct.nvtx[5][iP] - 1, my_global_temperature_struct.nvtx[6][iP] - 1, my_global_temperature_struct.nvtx[7][iP] - 1);
		}
	}
#endif


	// �� ����� ������ ������ ���� � ��������� ����������� ������� �������������.
	// ������� ������������� O(h!2) ������� ������� ������������ ��������� ������������ ������,
	// �.�. ������� �� ����������� ��������� �����.
	for (integer i = 0; i < flow_interior; i++) {
#if doubleintprecision == 1
		printf("FLUID %lld\n", i);
#else
		printf("FLUID %d\n", i);
#endif

		// �������� � ������� ���������������� ��������� ��� �������� ��������.
		f[i].resICCG = rterminate_residual_ICCG_Oh2(f[i]); // O(h!2)
														   //printf("residual O(h!2) is equal=%e\n", f[i].resICCG);
		std::cout << "residual O(h!2) is equal=" << f[i].resICCG << std::endl;
		f[i].resLR1sk = rterminate_residual_LR1sk_Oh3(f[i]); // O(h!3)
															 //printf("residual O(h!3) is equal=%e\n", f[i].resLR1sk);
		std::cout << "residual O(h!3) is equal=" << f[i].resLR1sk << std::endl;
	}
	printf("TEMPERATURE\n");
	my_global_temperature_struct.resLR1sk = rterminate_residual_LR1sk_temp_Oh3(my_global_temperature_struct); // O(h!3)		
																											  //printf("temp residual O(h!3) is equal=%e\n", t.resLR1sk);
	std::cout << "temp residual O(h!3) is equal=" << my_global_temperature_struct.resLR1sk << std::endl;
	printf("mesh check.\n");
	if (bwait) {
		//system("PAUSE");
		system("pause");
	}

	if (steady_or_unsteady_global_determinant != PHYSICAL_MODEL_SWITCH::NETWORK_T_UNSTEADY) {
		sourse2Dproblem = new bool[my_global_temperature_struct.maxbound];
		conductivity2Dinsource = new doublereal[my_global_temperature_struct.maxbound];
		// ������ ���������� ������������� ����� ��� ������������ �������.
		bsource_term_radiation_for_relax = new doublereal[my_global_temperature_struct.maxelm];
		for (integer i_init = 0; i_init < my_global_temperature_struct.maxelm; i_init++) bsource_term_radiation_for_relax[i_init] = 0.0;
		b_buffer_correct_source = new doublereal[my_global_temperature_struct.maxelm];
	}

	// ������� continity ����� ���������� �� ��������� � ������ 1e0.
	doublereal* continity_start = nullptr;
	continity_start = new doublereal[flow_interior];
	if (continity_start == nullptr) {
		// ������������ ������ �� ������ ������������.
		printf("Problem: not enough memory on your equipment for continity start in main...\n");
		printf("Please any key to exit...\n");
		exit(1);
	}
	for (integer i = 0; i < flow_interior; i++) continity_start[i] = 1.0;

	integer* inumber_iteration_SIMPLE = nullptr;
	inumber_iteration_SIMPLE = new integer[flow_interior];
	if (nullptr == inumber_iteration_SIMPLE) {
		// ������������ ������ �� ������ ������������.
		printf("Problem: not enough memory on your equipment for inumber_iteration_SIMPLE in main...\n");
		printf("Please any key to exit...\n");
		exit(1);
	}
	for (integer i = 0; i < flow_interior; i++) inumber_iteration_SIMPLE[i] = 0; // ��������� �������� ��������� SIMPLE ��� ������ FLUID ����.

																				 // ���������� ��������� ������� �� ����� ��� ������������� �������
	bool breadOk = false;
	avtoreadvalue(f, my_global_temperature_struct, flow_interior, inumber_iteration_SIMPLE, continity_start, breadOk, b, lb, s, ls, w, lw);
	// ���� ���������� ������ ���������, �� breadOk==false � ��� ������ ��� ���� ������� ������ �� �������� �������� ��� �������������.

	if (b_on_adaptive_local_refinement_mesh) {
		// ��������� ������������ ���� �����.
		printf("the invariant correctness...\n");
		ANES_ALICE_CORRECT(my_global_temperature_struct.maxnod, my_global_temperature_struct.pa, my_global_temperature_struct.maxelm, my_global_temperature_struct.nvtx);
	}

	// ������� ���������� ���������� � ��������� tecplot360:
	// ����� ������������ ��� �������� ����������� �����.
	if (0) {
		exporttecplotxy360T_3D_part2(my_global_temperature_struct.maxelm, my_global_temperature_struct.ncell, f, my_global_temperature_struct, flow_interior, 0, bextendedprint, 0, b, lb);
		printf("read values. OK.\n");
		//system("PAUSE"); // debug avtosave
		system("pause");
	}

#ifndef NO_OPENGL_GLFW
	scale_all = fmin(SCREEN_WIDTH / (fabs(b[0].g.xE - b[0].g.xS)), SCREEN_HEIGHT / (fabs(b[0].g.yE - b[0].g.yS)));

	//halfScreenWidth -= scale_all * 0.5 * (b[0].g.xE + b[0].g.xS);
	//halfScreenHeight -= scale_all * 0.5 * (b[0].g.yE + b[0].g.yS);

	halfScreenWidth = SCREEN_WIDTH / 2.0;
	halfScreenHeight = SCREEN_HEIGHT / 2.0;

#endif

	for (integer i = 0; i < my_global_temperature_struct.maxbound; i++) {
		if (my_global_temperature_struct.border_neighbor[i].iB > my_global_temperature_struct.maxbound + my_global_temperature_struct.maxelm) {
			printf("MCB=%lld  node error %lld maxelm=%e, maxbound=%e\n", my_global_temperature_struct.border_neighbor[i].MCB, i, my_global_temperature_struct.maxelm, my_global_temperature_struct.maxbound);
			getchar();
		}
		if (my_global_temperature_struct.border_neighbor[i].iI > my_global_temperature_struct.maxelm) {
			printf("MCB=%lld node error %lld maxelm=%e, maxbound=%e\n", my_global_temperature_struct.border_neighbor[i].MCB, i, my_global_temperature_struct.maxelm, my_global_temperature_struct.maxbound);
			getchar();
		}
	}

	/*for (integer i=0; i<lw; i++) {
	printf("%e  \n",w[i].Tamb);
	}
	//exporttecplotxy360T_3D(t.maxelm, t.ncell, t.nvtx, t.nvtxcell, t.pa, t.potent);
	exporttecplotxy360T_3D_part2(t.maxelm, t.potent);
	system("PAUSE"); // debug
	*/


	// 29.01.2017
	// if (1 && steady_or_unsteady_global_determinant == MESHER_ONLY)  
	// �� �� ������ �������� ����� �� ������� �������.
	if (1 && ((steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::MESHER_ONLY) ||
		(steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::NETWORK_T))) {
		// ����� �������.
		unsigned int calculation_start_time = 0; // ������ ����� ��.
		unsigned int calculation_end_time = 0; // ��������� ����� ��.
		unsigned int calculation_seach_time = 0; // ����� ���������� ������� ���� � ��.

		calculation_start_time = clock(); // ������ ������ �����.


#ifndef NO_OPENGL_GLFW
		pa_opengl = new TOCHKA[my_global_temperature_struct.database.maxelm];
		pa_render = new TOCHKA[my_global_temperature_struct.database.maxelm];
		n_render = my_global_temperature_struct.database.maxelm;
		for (int i = 0; i < n_render; i++) {
			pa_opengl[i].x = my_global_temperature_struct.database.x[i];
			pa_opengl[i].y = my_global_temperature_struct.database.y[i];
			pa_opengl[i].z = my_global_temperature_struct.database.z[i];

			pa_render[i].x = halfScreenWidth + scale_all * my_global_temperature_struct.database.x[i];
			pa_render[i].y = halfScreenHeight + scale_all * my_global_temperature_struct.database.y[i];
			pa_render[i].z = -abbys + scale_all * my_global_temperature_struct.database.z[i];
		}

		DrawZbuffer(b, lb, my_global_temperature_struct.database.ncell, my_global_temperature_struct.database.maxelm, my_global_temperature_struct.database.nvtxcell, my_global_temperature_struct.whot_is_block);
		delete[] pa_opengl;
		delete[] pa_render;
		n_render = -1;

#endif

		// ������������ ������� ����������� ��������� ������������ � .stl �������.
		// 09.09.2019.
		// �������������� 27.05.2020
		//if (steady_or_unsteady_global_determinant == MESHER_ONLY) {
		//export_User_Geom_in_STL_format(my_global_temperature_struct);
		//}

		// ���������� ����� ������.
		massa_cabinet(my_global_temperature_struct, f, flow_interior,
			b, lb, operatingtemperature,
			matlist);

		calculation_end_time = clock(); // ������ ��������� �����.
		calculation_seach_time = calculation_end_time - calculation_start_time;
		unsigned int im = 0, is = 0, ims = 0;
		im = (unsigned int)(calculation_seach_time / 60000); // ������
		is = (unsigned int)((calculation_seach_time - 60000 * im) / 1000); // �������
		ims = (unsigned int)((calculation_seach_time - 60000 * im - 1000 * is) / 10); // ������������ ������� �� 10

																					  // 24.06.2020
																					  // �������� ����� ������� ��������� �������������.
		if ((steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::NETWORK_T)) {
			calculate_Network_T(my_global_temperature_struct, b, lb, w, lw, ls, matlist);
		}

		printf("time export to tecplot360 is:  %u minute %u second %u millisecond\n", im, is, 10 * ims);
		//system("pause");
		//exit(1);

		// 1 - solver/solid_static/
		bool bMechanical = false;
		report_temperature(flow_interior, f, my_global_temperature_struct, b, lb, s, ls, w, lw, 0, matlist, bMechanical);
		// �������� ������ �������� ����� �������� ������� ������.
		// �������� ���������� (���������) �� ������� ����� �������� ����� � ��,
		// ���������� ����� �������� ������� ������. 28,10,2019
		//report_out_boundary(f[0], my_global_temperature_struct, ls, lw, w, b, lb, matlist, f[0].OpTemp);

		// ������� ����� � tecplot 360.
		if (1) {
			if (!b_on_adaptive_local_refinement_mesh) {
				if (0 == lu) {
					// ������� ���������� ���������� � ��������� tecplot360:
					exporttecplotxy360T_3D_part2(my_global_temperature_struct.maxelm, my_global_temperature_struct.ncell, f, my_global_temperature_struct, flow_interior, 0, bextendedprint, 0, b, lb);
				}
				else {
					//exporttecplotxy360T_3D_part2_assembles(t.maxelm, t.ncell, f, t, flow_interior, 0, bextendedprint, 0, lu, my_union);
					// ����� ��� ��������.
					exporttecplot_assembles_mesh(my_global_temperature_struct, lu, my_union);
				}
			}
			else {
				// ������� � ��������� tecplot �����������.
				//� ���� �����.
				ANES_tecplot360_export_temperature(my_global_temperature_struct.maxnod, my_global_temperature_struct.pa, my_global_temperature_struct.maxelm, my_global_temperature_struct.nvtx, my_global_temperature_struct.potent, my_global_temperature_struct, f, 0, b, lb);
			}
		}


	}

	//char ch_EXPORT_ALICE_ONLY = 'y';

	// steady Temperature Finite Volume Method
	if (1 && (steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::STEADY_TEMPERATURE) &&
		(1 == eqin.itemper)) {

		// ����� �������.
		unsigned int calculation_start_time = 0; // ������ ����� ��.
		unsigned int calculation_end_time = 0; // ��������� ����� ��.
		unsigned int calculation_seach_time = 0; // ����� ���������� ������� ���� � ��.

		calculation_start_time = clock(); // ������ ������ �����.

#pragma omp parallel for
		for (integer i7 = 0; i7 < my_global_temperature_struct.maxelm + my_global_temperature_struct.maxbound; i7++)
			my_global_temperature_struct.potent[i7] = operating_temperature_for_film_coeff; // �������������.

																							// ������ ������ ������������� � ������ ����.
		bonly_solid_calculation = true;

		// �������� ����������� ���������� �� ����������� ������.
		if (1 == lw) {
			bPhysics_stop = true;
			if (lb < 11) {
				// ��� ����������� ��������:
				// MD40, AuSn, Cu, AuSn, SiC, GaN. cabinet and hollow.
				bPhysics_PTBSH_memory = true;
			}
		}

		if (adiabatic_vs_heat_transfer_coeff == DEFAULT_CABINET_BOUNDARY_CONDITION::NEWTON_RICHMAN_BC) {
			// �� �������������� ��������� ��������� ����� ��������� �� ��������������� ��� ������ ������ ����.
			//for (integer i7 = 0; i7<t.maxelm + t.maxbound; i7++) t.potent[i7] = 0.57*operating_temperature_for_film_coeff;
		}

		// ����� �������������� ��� �� ������ ������������ ������ ������ ����������������.
		bsolid_static_only = true;
		bool bcleantemp = false;
		// ������� ��� ����������� ����.
		//if (1 == eqin.itemper) 
		{
			bcleantemp = true;
			for (integer i = 0; i < flow_interior; i++) {
				if (1 == eqin.fluidinfo[i].iflow) bcleantemp = false;
			}
			// ���� bcleantemp==true �� �� ������ ������ ������ ������������� ��� ����� ���������.
		}

		if (1 || bcleantemp) {
			// ������� ������������ ���������� (��� ��������) ������ ������ ���������������� � ��������� �������. 
			printf("solution of pure heat...\n");
			printf("please, press any key to continue...\n");
			if (bwait) {
				//system("PAUSE");
				system("pause");
			}

			// ��� ������������ ������������� ����������� ��������.
			bool bprintmessage = true; // �������� �� ��������� �� �������.

			doublereal dbeta = 1.0; // ������ ������� ������������� �� �������.
			bool bmyconvective = false;
			if (starting_speed_Vx * starting_speed_Vx + starting_speed_Vy * starting_speed_Vy + starting_speed_Vz * starting_speed_Vz > 1.0e-30) {
				if (f[0].maxelm > 0) {
					bmyconvective = true;
				}
			}
			else {
				// �������� ������������� ��������� ��������.

				FILE* fp_inicialization_data = NULL;
#ifdef MINGW_COMPILLER
				int err_inicialization_data = 0;
				fp_inicialization_data = fopen64("load.txt", "r");
				if (fp_inicialization_data == NULL) err_inicialization_data = 1;
#else
				errno_t err_inicialization_data = 0;
				err_inicialization_data = fopen_s(&fp_inicialization_data, "load.txt", "r");
#endif
				if (0 == err_inicialization_data) {
					// �������� ������ � ���� ������������.
					if (f[0].maxelm > 0) {
						// ������ � ������ ������� ������ ����� �� ������ ������ � ����������.
						bmyconvective = true;
					}
					fclose(fp_inicialization_data);
				}
				//printf("maxelm flow = %lld\n",f[0].maxelm);
				//system("PAUSE");
			}

			if (bmyconvective) {
				bmyconvective7248 = true;
			}

			// if (flow_interior>0) bmyconvective=true;
			// ������ ���������� ����������,
			// ��������� ��� �������� ������� ���-��� 1983
			doublereal** rhie_chow = nullptr;
			QuickMemVorst m;
			m.ballocCRSt = false; // �������� ������
			m.bsignalfreeCRSt = true; // � ����� �����������.
									  // ������������� ����������.
			m.tval = nullptr;
			m.tcol_ind = nullptr;
			m.trow_ptr = nullptr;
			m.tri = nullptr;
			m.troc = nullptr;
			m.ts = nullptr;
			m.tt = nullptr;
			m.tvi = nullptr;
			m.tpi = nullptr;
			m.tdx = nullptr;
			m.tdax = nullptr;
			m.ty = nullptr;
			m.tz = nullptr;
			m.ta = nullptr;
			m.tja = nullptr;
			m.tia = nullptr;
			m.talu = nullptr;
			m.tjlu = nullptr;
			m.tju = nullptr;
			m.tiw = nullptr;
			m.tlevs = nullptr;
			m.tw = nullptr;
			m.tjw = nullptr;
			m.icount_vel = 100000; // ����� ������� �����.			

			{
				integer* color = nullptr;
				integer dist_max = 3;
				calculate_color_for_temperature(color, my_global_temperature_struct, inx, xpos);

				// ���� flow_interior == 0 �� f[0] ������ ���������� ��������  
				solve_nonlinear_temp(f[0], f, my_global_temperature_struct,
					rhie_chow,
					b, lb, s, ls, w, lw,
					dbeta, flow_interior,
					bmyconvective, nullptr, 0.001, 0.001,
					false,
					matlist, 0,
					bprintmessage,
					gtdps, ltdp, 1.0, 1.0, // ��������� �������� ������ 1.0 �������� ��� �������� �������.
					m, nullptr, // �������� � ����������� ���������� ����. 
					nullptr,
					lu, my_union, color, dist_max); // �������� ����� ����� ������� � ����������� ���������� ����.
				delete[] color;
			}

			while ((bglobal_restart_06_10_2018)) {

				// ��������� ����� ���� ����������� ���������. ��������� ������������ ������.
				if (my_global_temperature_struct.potent != nullptr) {
					delete[] my_global_temperature_struct.potent;
					my_global_temperature_struct.potent = nullptr;
					my_global_temperature_struct.potent = new doublereal[my_global_temperature_struct.maxelm + my_global_temperature_struct.maxbound];
					for (integer i7 = 0; i7 < my_global_temperature_struct.maxelm + my_global_temperature_struct.maxbound; i7++) my_global_temperature_struct.potent[i7] = operating_temperature_for_film_coeff; // �������������.
					if (bsource_term_radiation_for_relax != nullptr) {
						delete[] bsource_term_radiation_for_relax;
						bsource_term_radiation_for_relax = nullptr;
						bsource_term_radiation_for_relax = new doublereal[my_global_temperature_struct.maxelm];
					}
					if (sourse2Dproblem != nullptr) {
						delete[] sourse2Dproblem;
						sourse2Dproblem = nullptr;
						sourse2Dproblem = new bool[my_global_temperature_struct.maxbound];
					}
					if (conductivity2Dinsource != nullptr) {
						delete[] conductivity2Dinsource;
						conductivity2Dinsource = nullptr;
						conductivity2Dinsource = new doublereal[my_global_temperature_struct.maxbound];
					}

					if (rthdsd_no_radiosity_patch != nullptr) {
						delete[]	rthdsd_no_radiosity_patch;
						rthdsd_no_radiosity_patch = nullptr;
					}

					if (b_buffer_correct_source != nullptr) {
						delete[] b_buffer_correct_source;
						b_buffer_correct_source = nullptr;
						b_buffer_correct_source = new doublereal[my_global_temperature_struct.maxelm];
					}

					if (my_global_temperature_struct.slau != nullptr) {
						delete[] my_global_temperature_struct.slau;
						my_global_temperature_struct.slau = nullptr;
						my_global_temperature_struct.slau = new equation3D[my_global_temperature_struct.maxelm]; // ������������ ������� ���� ��� ���������� ��.
						if (my_global_temperature_struct.slau == nullptr) {
							// ������������ ������ �� ������ ������������.
							printf("Problem: not enough memory on your equipment for slau temperature constr struct...\n");
							printf("Please any key to exit...\n");
							//system("PAUSE");
							system("pause");
							exit(1);
						}
					}

					if (my_global_temperature_struct.slau_bon != nullptr) {
						delete[] my_global_temperature_struct.slau_bon;
						my_global_temperature_struct.slau_bon = nullptr;
						my_global_temperature_struct.slau_bon = new equation3D_bon[my_global_temperature_struct.maxbound]; // ������������ ������� ���� ��� ��������� ��
						if (my_global_temperature_struct.slau_bon == nullptr) {
							// ������������ ������ �� ������ ������������.
							printf("Problem: not enough memory on your equipment for slau boundary temperature constr struct...\n");
							printf("Please any key to exit...\n");
							//system("PAUSE");
							system("pause");
							exit(1);
						}
					}
				}
				//bglobal_restart_06_10_2018 = false;

				integer* color = nullptr;
				integer dist_max = 3;
				calculate_color_for_temperature(color, my_global_temperature_struct, inx, xpos);


				// ���� flow_interior == 0 �� f[0] ������ ���������� ��������  
				solve_nonlinear_temp(f[0], f, my_global_temperature_struct,
					rhie_chow,
					b, lb, s, ls, w, lw,
					dbeta, flow_interior,
					bmyconvective, nullptr, 0.001, 0.001,
					false,
					matlist, 0,
					bprintmessage,
					gtdps, ltdp, 1.0, 1.0, m,
					nullptr, // �������� � ����������� ���������� ����. 
					nullptr,
					lu, my_union, color, dist_max); // �������� ����� ����� ������� � ����������� ���������� ����.
													// ��������� �������� ������ 1.0 �������� ��� �������� �������.
				delete[] color;

			}

			// ��������� ������� ����������� � ������ �������.
			doublereal Tavg_fluid = 0.0;
			doublereal ic_avg_Temp_fluid = 0.0;

#pragma omp parallel for reduction(+ : Tavg_fluid, ic_avg_Temp_fluid)
			for (integer i_7 = 0; i_7 < f[0].maxelm; i_7++) {
				if ((f[0].ptr[i_7] >= 0) && (f[0].ptr[i_7] < my_global_temperature_struct.maxelm)) {
					Tavg_fluid += my_global_temperature_struct.potent[f[0].ptr[i_7]];
					ic_avg_Temp_fluid += 1.0;
				}
			}
			if (ic_avg_Temp_fluid > 1.0e-30) {
				Tavg_fluid /= ic_avg_Temp_fluid;
				printf("average fluid temperature is %e\n", Tavg_fluid);
			}
			else {
				//printf("no fluid cell\n");
			}



			// ���������� ����� ������.
			massa_cabinet(my_global_temperature_struct, f, flow_interior,
				b, lb, operatingtemperature,
				matlist);

			// 10.10.2017
			// ���������� ���������� ������� ����� ���������.
			xyplot_temp(my_global_temperature_struct, my_global_temperature_struct.potent);
			//printf("graphics writing sucseful\n");
			//system("PAUSE");





			if (1) {
				if (!b_on_adaptive_local_refinement_mesh) {

#ifndef NO_OPENGL_GLFW
					pa_opengl = new TOCHKA[my_global_temperature_struct.database.maxelm];
					pa_render = new TOCHKA[my_global_temperature_struct.database.maxelm];
					n_render = my_global_temperature_struct.database.maxelm;
					for (int i = 0; i < n_render; i++) {
						pa_opengl[i].x = my_global_temperature_struct.database.x[i];
						pa_opengl[i].y = my_global_temperature_struct.database.y[i];
						pa_opengl[i].z = my_global_temperature_struct.database.z[i];

						pa_render[i].x = halfScreenWidth + scale_all * my_global_temperature_struct.database.x[i];
						pa_render[i].y = halfScreenHeight + scale_all * my_global_temperature_struct.database.y[i];
						pa_render[i].z = -abbys + scale_all * my_global_temperature_struct.database.z[i];
					}


					doublereal dmin = 1.0e30;
					doublereal dmax = -1.0e30;

					for (int i37 = 0; i37 < my_global_temperature_struct.maxelm; i37++) {
						if (my_global_temperature_struct.potent[i37] > dmax) {
							dmax = my_global_temperature_struct.potent[i37];
						}
						if (my_global_temperature_struct.potent[i37] < dmin) {
							dmin = my_global_temperature_struct.potent[i37];
						}
					}

					minimum_val_for_render_pic = dmin;
					maximum_val_for_render_pic = dmax;

					//DrawZbufferColor(b, lb, my_global_temperature_struct.database.ncell, my_global_temperature_struct.database.maxelm, my_global_temperature_struct.database.nvtxcell, my_global_temperature_struct.whot_is_block, my_global_temperature_struct.potent);
					DrawZbufferColorLight(b, lb, my_global_temperature_struct.database.ncell, my_global_temperature_struct.database.maxelm, my_global_temperature_struct.database.nvtxcell, my_global_temperature_struct.whot_is_block, my_global_temperature_struct.potent);
					delete[] pa_opengl;
					delete[] pa_render;
					n_render = -1;

#endif

					calculation_end_time = clock(); // ������ ��������� �����.
					calculation_seach_time = calculation_end_time - calculation_start_time;
					unsigned int im = 0, is = 0, ims = 0;
					im = (unsigned int)(calculation_seach_time / 60000); // ������
					is = (unsigned int)((calculation_seach_time - 60000 * im) / 1000); // �������
					ims = (unsigned int)((calculation_seach_time - 60000 * im - 1000 * is) / 10); // ������������ ������� �� 10

					printf("time calculation is:  %u minute %u second %u millisecond\n", im, is, 10 * ims);
					printf("export to tecplot start...\n");

					// ������� ���������� ���������� � ��������� tecplot360:
					exporttecplotxy360T_3D_part2(my_global_temperature_struct.maxelm, my_global_temperature_struct.ncell, f, my_global_temperature_struct, flow_interior, 0, bextendedprint, 0, b, lb);
				}
				else {
					if (b_on_adaptive_local_refinement_mesh) {

						calculation_main_end_time = clock(); // ������ ��������� �����.
						unsigned int calculation_main_seach_time = calculation_main_end_time - calculation_main_start_time_global_Depend;

						// ����� ����� ���������� �� ������������� ������������ ���������.
						int im = 0, is = 0, ims = 0;
						im = (int)(calculation_main_seach_time / 60000); // ������
						is = (int)((calculation_main_seach_time - 60000 * im) / 1000); // �������
						ims = (int)((calculation_main_seach_time - 60000 * im - 1000 * is) / 10); // ������������ ������� �� 10

						printf("time calculation is:  %d minute %d second %d millisecond\n", im, is, 10 * ims);
						printf("export to tecplot start...\n");

						//printf("Would you like to save the result on the ALICE grid ? y/n\n");
						//ch_EXPORT_ALICE_ONLY = getchar(); // ����� ������ getchar();
						//ch_EXPORT_ALICE_ONLY = 'y';
					}

					if (ch_EXPORT_ALICE_ONLY == 'y') {
						// ������� � ��������� tecplot �����������.
						//� ���� �����.
						ANES_tecplot360_export_temperature(my_global_temperature_struct.maxnod, my_global_temperature_struct.pa, my_global_temperature_struct.maxelm, my_global_temperature_struct.nvtx, my_global_temperature_struct.potent, my_global_temperature_struct, f, 0, b, lb);
					}

				}
			}

		}

		doublereal tmaxfinish = -273.15; // ���������� ����.
										 // ���������� �������� ������������ ����������� ������ ��������� ������� � �� � ��������:
										 //for (integer i = 0; i < t.maxelm + t.maxbound; i++) tmaxfinish = fmax(tmaxfinish, fabs(t.potent[i]));
										 // 23 ������� 2015
										 // �� ��������� ������ ���������� ����� �� ����� ��������� ������� �����������, �������
										 // �������� �� ������� ����� � ��������� ����������� ������ �� ���������� ��. 
		for (integer i = 0; i < my_global_temperature_struct.maxelm; i++) tmaxfinish = fmax(tmaxfinish, my_global_temperature_struct.potent[i]);

		FILE* fp = NULL;

#ifdef MINGW_COMPILLER
		int err1 = 0;
		fp = fopen64("report.txt", "w");
		if (fp == NULL) err1 = 1;
#else
		errno_t err1 = 0;
		err1 = fopen_s(&fp, "report.txt", "w");
#endif
		// �������� ����� ��� ������.
		if ((err1) != 0) {
			printf("Create File report.txt Error\n");
			//system("PAUSE");
			system("pause");
		}
		else {
			// ������ ���������
			fprintf(fp, "Maximum Temperature %.2f\n", tmaxfinish);
			fclose(fp);
		}
		// 1 - solver/solid_static/
		bool bMechanical = false;
		report_temperature(flow_interior, f, my_global_temperature_struct, b, lb, s, ls, w, lw, 0, matlist, bMechanical);
		// �������� ������ �������� ����� �������� ������� ������.
		// �������� ���������� (���������) �� ������� ����� �������� ����� � ��,
		// ���������� ����� �������� ������� ������. 28.10.2019
		report_out_boundary(f[0], my_global_temperature_struct, ls, lw, w, b, lb, matlist, f[0].OpTemp);


		if (ch_EXPORT_ALICE_ONLY != 'y') {

			calculation_end_time = clock(); // ������ ��������� �����.
			calculation_seach_time = calculation_end_time - calculation_start_time;
			unsigned int im = 0, is = 0, ims = 0;
			im = (unsigned int)(calculation_seach_time / 60000); // ������
			is = (unsigned int)((calculation_seach_time - 60000 * im) / 1000); // �������
			ims = (unsigned int)((calculation_seach_time - 60000 * im - 1000 * is) / 10); // ������������ ������� �� 10

			printf("time calculation is:  %u minute %u second %u millisecond\n", im, is, 10 * ims);

			// �� �������� ������ ������� ����������� (11.4� ����� �� ���� Coarse )
			// ����������� ����������� �����. 

			// 25.11.2017
			// 1. �������� ������ � ����������� � ����� �� ����.
			// 2. ��������� �� � ����������� ������.
			// 3. ���������� ������.
			// 4. ��������� ������� ���������� ������������� �����.
			// 5. ��������� ������ � ����������� � ���� �� ������� ���������� ������������� �����.
			// 6. ��������������� ����������� � ����������� �� ����������������� ����� �������� ������.
			// 7.1 Tecplot �������� ������������� �� ��� ���� �� ����������������� ����� � � ������� � � ������. 
			// 7.2 ������������ �� ���� ����� � ������ ���� �������� ����� ���������� ���� �� ����������������� (���� tecplot�).


			if (b_on_adaptive_local_refinement_mesh) {
				// 1. ��������� x,y,z,T,nvtx, m_sizeT, m_size_nvtx.
				doublereal* x_buf = nullptr;
				doublereal* y_buf = nullptr;
				doublereal* z_buf = nullptr;
				doublereal* t_buf = nullptr;
				integer** nvtx_buf = nullptr;
				integer m_sizeT = 0, m_size_nvtx = 0;

				ANES_tecplot360_export_temperature_preobrazovatel(my_global_temperature_struct.maxnod, my_global_temperature_struct.pa, my_global_temperature_struct.maxelm, my_global_temperature_struct.nvtx, my_global_temperature_struct.potent, my_global_temperature_struct, x_buf, y_buf, z_buf, t_buf, nvtx_buf, m_sizeT, m_size_nvtx, operatingtemperature);

				// 2. ������������ ������.
				// ������������ ����������� ������.
				if (my_global_temperature_struct.xpos_copy != nullptr) {
					delete[] my_global_temperature_struct.xpos_copy;
					my_global_temperature_struct.xpos_copy = nullptr;
				}
				if (my_global_temperature_struct.ypos_copy != nullptr) {
					delete[] my_global_temperature_struct.ypos_copy;
					my_global_temperature_struct.ypos_copy = nullptr;
				}
				if (my_global_temperature_struct.zpos_copy != nullptr) {
					delete[] my_global_temperature_struct.zpos_copy;
					my_global_temperature_struct.zpos_copy = nullptr;
				}


				if (bsource_term_radiation_for_relax != nullptr) {
					delete[] bsource_term_radiation_for_relax; // ���������� ������������ ������ ������������ �������.
					bsource_term_radiation_for_relax = nullptr;
				}
				if (b_buffer_correct_source != nullptr) {
					delete[] b_buffer_correct_source;
					b_buffer_correct_source = nullptr;
				}

				if (rthdsd_no_radiosity_patch != nullptr) {
					delete[] rthdsd_no_radiosity_patch;
					rthdsd_no_radiosity_patch = nullptr;
				}


				// ������� ��������� ���������� ��������� ������� � ����� 0.14 ��������.
				if (qnbc != nullptr) {
					delete[] qnbc;
					qnbc = nullptr;
					iadd_qnbc_maxelm = 0;
				}

				// ����� ���������� ����������� ������ �� ��� ���� �������� ������:
				free_level1_temp(my_global_temperature_struct);
				free_level2_temp(my_global_temperature_struct); // ������������ ������ �� ��� ������.
																// ����������� ������ ��� LR ������.
				if (my_global_temperature_struct.rootWE != nullptr) {
					free_root(my_global_temperature_struct.rootWE, my_global_temperature_struct.iWE);
				}
				if (my_global_temperature_struct.rootSN != nullptr) {
					free_root(my_global_temperature_struct.rootSN, my_global_temperature_struct.iSN);
				}
				if (my_global_temperature_struct.rootBT != nullptr) {
					free_root(my_global_temperature_struct.rootBT, my_global_temperature_struct.iBT);
				}
				if (my_global_temperature_struct.rootWE != nullptr) {
					delete[] my_global_temperature_struct.rootWE;
					my_global_temperature_struct.rootWE = nullptr;
				}
				if (my_global_temperature_struct.rootSN != nullptr) {
					delete[] my_global_temperature_struct.rootSN;
					my_global_temperature_struct.rootSN = nullptr;
				}
				if (my_global_temperature_struct.rootBT != nullptr) {
					delete[] my_global_temperature_struct.rootBT;
					my_global_temperature_struct.rootBT = nullptr;
				}

				if (bvery_big_memory) {
					if (my_global_temperature_struct.database.x != nullptr) {
						free(my_global_temperature_struct.database.x);
					}
					if (my_global_temperature_struct.database.y != nullptr) {
						free(my_global_temperature_struct.database.y);
					}
					if (my_global_temperature_struct.database.z != nullptr) {
						free(my_global_temperature_struct.database.z);
					}
					if (my_global_temperature_struct.database.nvtxcell != nullptr) {
						for (integer i = 0; i <= 7; i++) {
							delete[] my_global_temperature_struct.database.nvtxcell[i];
						}
						delete[] my_global_temperature_struct.database.nvtxcell;
					}
					if (my_global_temperature_struct.database.ptr != nullptr) {
						if (my_global_temperature_struct.database.ptr[0] != nullptr) {
							delete[] my_global_temperature_struct.database.ptr[0];
						}
						if (my_global_temperature_struct.database.ptr[1] != nullptr) {
							delete[] my_global_temperature_struct.database.ptr[1];
						}
						delete[] my_global_temperature_struct.database.ptr;
					}
				}

				// ������������ ������ ��� LR �����.
				free_level1_flow(f, flow_interior);
				free_level2_flow(f, flow_interior); // ������������ ������ �� ��� ������.

				if (sourse2Dproblem != nullptr) {
					delete[] sourse2Dproblem;
					sourse2Dproblem = nullptr;
				}
				if (conductivity2Dinsource != nullptr) {
					delete[] conductivity2Dinsource;
					conductivity2Dinsource = nullptr;
				}

				if (x_jacoby_buffer != nullptr) {
					// 30 ������� 2016. 
					// � seidelsor2 ������ ������������� �� ����� ������ ���������� �.�. �����.
					// ������������ ������ �� ��� Jacobi buffer.
					delete[] x_jacoby_buffer;
				}



				// ������������ ����� ������ � ILU �������.
				//if (milu_gl_buffer.alu_copy != nullptr) delete[] milu_gl_buffer.alu_copy;
				//if (milu_gl_buffer.jlu_copy != nullptr) delete[] milu_gl_buffer.jlu_copy;
				//if (milu_gl_buffer.ju_copy != nullptr) delete[] milu_gl_buffer.ju_copy;
				//milu_gl_buffer.alu_copy = nullptr;
				//milu_gl_buffer.jlu_copy = nullptr;
				//milu_gl_buffer.ju_copy = nullptr;

				flow_interior = 0;

				// 3. ���������� ������� �����.

				b_on_adaptive_local_refinement_mesh = false;
				iCabinetMarker = 0;
				load_TEMPER_and_FLOW(my_global_temperature_struct, f, inx, iny, inz, xpos, ypos, zpos, flow_interior,
					b, lb, lw, w, s, ls, lu, my_union, operatingtemperature, matlist, bextendedprint,
					dgx, dgy, dgz, b_on_adaptive_local_refinement_mesh, false, iCabinetMarker);

				my_global_temperature_struct.operatingtemperature = operatingtemperature;


				// ��� ����� ������ ����� ��� ������� �������������� ������.
				my_global_temperature_struct.inx_copy = inx;
				my_global_temperature_struct.iny_copy = iny;
				my_global_temperature_struct.inz_copy = inz;
				my_global_temperature_struct.operatingtemperature_copy = operatingtemperature;
				my_global_temperature_struct.xpos_copy = new doublereal[inx + 1];
				my_global_temperature_struct.ypos_copy = new doublereal[iny + 1];
				my_global_temperature_struct.zpos_copy = new doublereal[inz + 1];
				// ������ ���������� ����� ��� �������� ����������� ������,
				// ��������� ������ ����� ��������� �� ��� � ����� ������������� 
				// ���� ����c�����.
				for (integer i_7 = 0; i_7 < inx + 1; i_7++) {
					my_global_temperature_struct.xpos_copy[i_7] = xpos[i_7];
				}
				for (integer i_7 = 0; i_7 < iny + 1; i_7++) {
					my_global_temperature_struct.ypos_copy[i_7] = ypos[i_7];
				}
				for (integer i_7 = 0; i_7 < inz + 1; i_7++) {
					my_global_temperature_struct.zpos_copy[i_7] = zpos[i_7];
				}

				my_global_temperature_struct.free_temper_level1 = false; // ������ ���������������� ������������ ������ ����������� ��� ������ ������� ����� �������� ������.
				my_global_temperature_struct.free_temper_level2 = false; // ������������ ������ ��� �������� ������� ��� ���������� � � SIMPLESPARSE ������.	


																		 // 4. ������������ ��� �����������.
				ALICE_2_Structural(my_global_temperature_struct.maxnod, my_global_temperature_struct.pa, my_global_temperature_struct.maxelm, my_global_temperature_struct.nvtx, my_global_temperature_struct.potent, x_buf, y_buf, z_buf, t_buf, nvtx_buf, m_sizeT, m_size_nvtx, my_global_temperature_struct.operatingtemperature_copy);


				if (x_buf != nullptr) {
					delete[] x_buf;
					x_buf = nullptr;
				}
				if (y_buf != nullptr) {
					delete[] y_buf;
					y_buf = nullptr;
				}
				if (z_buf != nullptr) {
					delete[] z_buf;
					z_buf = nullptr;
				}
				if (t_buf != nullptr) {
					delete[] t_buf;
					t_buf = nullptr;
				}
				if (nvtx_buf != nullptr) {
					for (integer i_1 = 0; i_1 < 8; i_1++) {
						if (nvtx_buf[i_1] != nullptr) {
							delete[] nvtx_buf[i_1];
							nvtx_buf[i_1] = nullptr;
						}
					}
					delete[] nvtx_buf;
					nvtx_buf = nullptr;
				}
				m_sizeT = 0, m_size_nvtx = 0;
				// 5. ������� ������� � tecplot.
				exporttecplotxy360T_3D_part2(my_global_temperature_struct.maxelm, my_global_temperature_struct.ncell, f, my_global_temperature_struct, flow_interior, 0, bextendedprint, 0, b, lb);
			}
		}
		else {
			// �������� ������������ �� ���� ����� ������ ���� � ������ �����, �������
			// ��� ������������ �������� ������������ ����� ������� �� ����������������� �����.
			// � ������� �� ���� ����� ��� ������ �� tecplot�.
			// �������� �� ���� ����� �����. ����� �������� ������� �������� ��� ���� ����������, ��
			// ��������� �������� ������� ������� �� ����� ������� (������������ �����) ������ � ������
			// ��������������� ������ ���� ������ ����� ����� ���� ����� ��������� (������ ������� �����������).

			// ������� � ��������� tecplot �����������.
			// � ���� �����.
			//ANES_tecplot360_export_temperature(t.maxnod, t.pa, t.maxelm, t.nvtx, t.potent,t,0,b,lb);
		}



	}

	// steady Temperature Finite Element Method
	if (1 && (steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::STEADY_TEMPERATURE) && (eqin.itemper == 2)) {

		// ����� �������.
		unsigned int calculation_start_time = 0; // ������ ����� ��.
		unsigned int calculation_end_time = 0; // ��������� ����� ��.
		unsigned int calculation_seach_time = 0; // ����� ���������� ������� ���� � ��.

		calculation_start_time = clock(); // ������ ������ �����.

#pragma omp parallel for
		for (integer i7 = 0; i7 < my_global_temperature_struct.maxelm + my_global_temperature_struct.maxbound; i7++) {
			my_global_temperature_struct.potent[i7] = operating_temperature_for_film_coeff; // �������������.
		}

		// ������ ������ Static Structural.
		bonly_solid_calculation = true;

		// �������� ����������� ���������� �� ����������� ������.
		if (lw == 1) {
			bPhysics_stop = true;
			if (lb < 11) {
				// ��� ����������� ��������:
				// MD40, AuSn, Cu, AuSn, SiC, GaN. cabinet and hollow.
				bPhysics_PTBSH_memory = true;
			}
		}

		if (adiabatic_vs_heat_transfer_coeff == DEFAULT_CABINET_BOUNDARY_CONDITION::NEWTON_RICHMAN_BC) {
			// �� �������������� ��������� ��������� ����� ��������� �� ��������������� ��� ������ ������ ����.
			//for (integer i7 = 0; i7<t.maxelm + t.maxbound; i7++) t.potent[i7] = 0.57*operating_temperature_for_film_coeff;
		}

		// ����� �������������� ��� �� ������ ������������ ������ ������ ����������������.
		bsolid_static_only = true;
		bool bcleantemp = false;
		// �� ���� ������ ������ ������ ���������������� ������
		// ������ ������������� ��������.
		//if (2 == eqin.itemper)
		{
			bcleantemp = true;
			for (integer i = 0; i < flow_interior; i++) {
				if (eqin.fluidinfo[i].iflow == 1) bcleantemp = false;
			}
			// ���� bcleantemp==true �� �� ������ ������ ������ ������������� ��� ����� ���������.
		}
		if (1 || bcleantemp) {
			// ������� ������������ ���������� (��� ��������) ������ ������ ���������������� � ��������� �������. 
			printf("solution of pure Static Structural...\n");
			printf("please, press any key to continue...\n");
			if (bwait) {
				//system("PAUSE");
				system("pause");
			}

			// ��� ������������ ������������� ����������� ��������.
			bool bprintmessage = true; // �������� �� ��������� �� �������.

			doublereal dbeta = 1.0; // ������ ������� ������������� �� �������.
			bool bmyconvective = false;
			if (starting_speed_Vx * starting_speed_Vx + starting_speed_Vy * starting_speed_Vy + starting_speed_Vz * starting_speed_Vz > 1.0e-30) {
				if (f[0].maxelm > 0) {
					bmyconvective = true;
				}
			}
			else {
				// �������� ������������� ��������� ��������.

				FILE* fp_inicialization_data = NULL;
#ifdef MINGW_COMPILLER
				int err_inicialization_data = 0;
				fp_inicialization_data = fopen64("load.txt", "r");
				if (fp_inicialization_data == NULL) err_inicialization_data = 1;
#else
				errno_t err_inicialization_data = 0;
				err_inicialization_data = fopen_s(&fp_inicialization_data, "load.txt", "r");
#endif

				if (err_inicialization_data == 0) {
					// �������� ������ � ���� ������������.
					if (f[0].maxelm > 0) {
						bmyconvective = true;
					}
					fclose(fp_inicialization_data);
				}
			}

			// if (flow_interior>0) bmyconvective=true;
			// ������ ���������� ����������,
			// ��������� ��� �������� ������� ���-��� 1983
			doublereal** rhie_chow = nullptr;
			QuickMemVorst m;
			m.ballocCRSt = false; // �������� ������
			m.bsignalfreeCRSt = true; // � ����� �����������.

									  // ������������� ����������.
			m.tval = nullptr;
			m.tcol_ind = nullptr;
			m.trow_ptr = nullptr;
			m.tri = nullptr;
			m.troc = nullptr;
			m.ts = nullptr;
			m.tt = nullptr;
			m.tvi = nullptr;
			m.tpi = nullptr;
			m.tdx = nullptr;
			m.tdax = nullptr;
			m.ty = nullptr;
			m.tz = nullptr;
			m.ta = nullptr;
			m.tja = nullptr;
			m.tia = nullptr;
			m.talu = nullptr;
			m.tjlu = nullptr;
			m.tju = nullptr;
			m.tiw = nullptr;
			m.tlevs = nullptr;
			m.tw = nullptr;
			m.tjw = nullptr;
			m.icount_vel = 100000; // ����� ������� �����.

			bPhysics_stop = false;

			// ����������� 19.05.2018
		    

			doublereal* lstub = nullptr;
			integer maxelm_global_ret = 0;
			doublereal* t_for_Mechanical = nullptr;
			solve_Thermal(my_global_temperature_struct, f, matlist, w, lw, lu, b, lb, ls, m, false,
				operatingtemperature, false, 0.0, lstub, lstub,
				maxelm_global_ret, 1.0, 1.0, bAVLrealesation, t_for_Mechanical);


			/*
			// ���� flow_interior == 0 �� f[0] ������ ���������� ��������
			solve_nonlinear_temp(f[0], f, t,
			rhie_chow,
			b, lb, s, ls, w, lw,
			dbeta, flow_interior,
			bmyconvective, nullptr, 0.001, 0.001,
			false,
			matlist, 0,
			bprintmessage,
			gtdps, ltdp, 1.0, m,
			nullptr, // �������� � ����������� ���������� ����.
			nullptr); // �������� ����� ����� ������� � ����������� ���������� ����.
			// ��������� �������� ������ 1.0 �������� ��� �������� �������.
			*/
			// ���������� ����� ������.
			massa_cabinet(my_global_temperature_struct, f, flow_interior,
				b, lb, operatingtemperature,
				matlist);


			calculation_end_time = clock(); // ������ ��������� �����.
			calculation_seach_time = calculation_end_time - calculation_start_time;
			unsigned int im = 0, is = 0, ims = 0;
			im = (unsigned int)(calculation_seach_time / 60000); // ������
			is = (unsigned int)((calculation_seach_time - 60000 * im) / 1000); // �������
			ims = (unsigned int)((calculation_seach_time - 60000 * im - 1000 * is) / 10); // ������������ ������� �� 10

			printf("time calculation is:  %u minute %u second %u millisecond\n", im, is, 10 * ims);

			if (1) {
				if (!b_on_adaptive_local_refinement_mesh) {
					// ������� ���������� ���������� � ��������� tecplot360:
					exporttecplotxy360T_3D_part2(my_global_temperature_struct.maxelm, my_global_temperature_struct.ncell, f, my_global_temperature_struct, flow_interior, 0, bextendedprint, 0, b, lb);
				}
				else {
					// ������� � ��������� tecplot �����������.
					//� ���� �����.
					ANES_tecplot360_export_temperature(my_global_temperature_struct.maxnod, my_global_temperature_struct.pa, my_global_temperature_struct.maxelm, my_global_temperature_struct.nvtx, my_global_temperature_struct.potent, my_global_temperature_struct, f, 0, b, lb);
				}
			}

		}

		doublereal tmaxfinish = -273.15; // ���������� ����.
										 // ���������� �������� ������������ ����������� ������ ��������� ������� � �� � ��������:
										 //for (integer i = 0; i < t.maxelm + t.maxbound; i++) tmaxfinish = fmax(tmaxfinish, fabs(t.potent[i]));
										 // 23 ������� 2015
										 // �� ��������� ������ ���������� ����� �� ����� ��������� ������� �����������, �������
										 // �������� �� ������� ����� � ��������� ����������� ������ �� ���������� ��. 
		for (integer i = 0; i < my_global_temperature_struct.maxelm; i++) tmaxfinish = fmax(tmaxfinish, my_global_temperature_struct.potent[i]);

		doublereal totaldeform_max = -1.0e+30;
		if ((steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::STEADY_STATIC_STRUCTURAL) ||
			(steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::STEADY_STATIC_STRUCTURAL_AND_TEMPERATURE) ||
			(steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::UNSTEADY_STATIC_STRUCTURAL) ||
			(steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::UNSTEADY_STATIC_STRUCTURAL_AND_TEMPERATURE))
		{
			for (integer i = 0; i < my_global_temperature_struct.maxelm; i++) totaldeform_max = fmax(totaldeform_max, my_global_temperature_struct.total_deformation[TOTALDEFORMATION][i]);
		}

		FILE* fp = NULL;

#ifdef MINGW_COMPILLER
		int err1 = 0;
		fp = fopen64("report.txt", "w");
		if (fp == NULL) err1 = 1;
#else
		errno_t err1 = 0;
		err1 = fopen_s(&fp, "report.txt", "w");
#endif
		// �������� ����� ��� ������.
		if ((err1) != 0) {
			printf("Create File report.txt Error\n");
			//system("PAUSE");
			system("pause");
		}
		else {
			// ������ ���������
			fprintf(fp, "Maximum Temperature %.2f\n", tmaxfinish);
			fclose(fp);
		}
		// 1 - solver/solid_static/
		bool bMechanical = false;
		report_temperature(flow_interior, f, my_global_temperature_struct, b, lb, s, ls, w, lw, 0, matlist, bMechanical);


	}

	// steady Static Structural
	if (1 && steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::STEADY_STATIC_STRUCTURAL) {


		// ����� �������.
		unsigned int calculation_start_time = 0; // ������ ����� ��.
		unsigned int calculation_end_time = 0; // ��������� ����� ��.
		unsigned int calculation_seach_time = 0; // ����� ���������� ������� ���� � ��.

		calculation_start_time = clock(); // ������ ������ �����.

#pragma omp parallel for
		for (integer i7 = 0; i7 < my_global_temperature_struct.maxelm + my_global_temperature_struct.maxbound; i7++) {
			my_global_temperature_struct.potent[i7] = operating_temperature_for_film_coeff; // �������������.
		}

		// ������ ������ Static Structural.
		bonly_solid_calculation = true;


		// ����� �������������� ��� �� ������ ������������ ������ ������ �������� �� ���������� ���������������� ���������.
		bsolid_static_only = true;
		bool bcleantemp = false;

		if (1 || bcleantemp) {
			// ������� ������������ ���������� (��� ��������) ������ ������ 
			// �������� �� ���������� ���������������� ��������� � ��������� �������. 
			printf("solution of pure Static Structural...\n");
			printf("please, press any key to continue...\n");
			if (bwait) {
				//system("PAUSE");
				system("pause");
			}

			// ��� ������������ ������������� ����������� ��������.
			bool bprintmessage = true; // �������� �� ��������� �� �������.

			doublereal dbeta = 1.0; // ������ ������� ������������� �� �������.
			bool bmyconvective = false;
			if (starting_speed_Vx * starting_speed_Vx +
				starting_speed_Vy * starting_speed_Vy +
				starting_speed_Vz * starting_speed_Vz > 1.0e-30) {
				if (f[0].maxelm > 0) {
					bmyconvective = true;
				}
			}
			else {
				// �������� ������������� ��������� ��������.

				FILE* fp_inicialization_data = NULL;
#ifdef MINGW_COMPILLER
				int err_inicialization_data = 0;
				fp_inicialization_data = fopen64("load.txt", "r");
				if (fp_inicialization_data == NULL) err_inicialization_data = 1;
#else
				errno_t err_inicialization_data = 0;
				err_inicialization_data = fopen_s(&fp_inicialization_data, "load.txt", "r");
#endif

				if (err_inicialization_data == 0) {
					// �������� ������ � ���� ������������.
					if (f[0].maxelm > 0) {
						bmyconvective = true;
					}
					fclose(fp_inicialization_data);
				}
			}


			bPhysics_stop = false;
			// ����� ������� Static Structural.
			// �������� 19.05.2018
			doublereal* stub_Mechanical = nullptr;
			solve_Structural(my_global_temperature_struct, w, lw, false, operatingtemperature, b, lb, lu,
				false, 1.0e9, 1.0e9, 1.0e9, 1.0e9, stub_Mechanical, stub_Mechanical, stub_Mechanical, stub_Mechanical, 1.0, matlist, stub_Mechanical);
			bPhysics_stop = true;


			// ���������� ����� ������.
			massa_cabinet(my_global_temperature_struct, f, flow_interior,
				b, lb, operatingtemperature,
				matlist);


			calculation_end_time = clock(); // ������ ��������� �����.
			calculation_seach_time = calculation_end_time - calculation_start_time;
			unsigned int im = 0, is = 0, ims = 0;
			im = (unsigned int)(calculation_seach_time / 60000); // ������
			is = (unsigned int)((calculation_seach_time - 60000 * im) / 1000); // �������
			ims = (unsigned int)((calculation_seach_time - 60000 * im - 1000 * is) / 10); // ������������ ������� �� 10

			printf("time calculation is:  %u minute %u second %u millisecond\n", im, is, 10 * ims);

			if (1) {
				if (!b_on_adaptive_local_refinement_mesh) {
					// ������� ���������� ���������� � ��������� tecplot360:
					exporttecplotxy360T_3D_part2(my_global_temperature_struct.maxelm, my_global_temperature_struct.ncell, f, my_global_temperature_struct, flow_interior, 0, bextendedprint, 0, b, lb);
				}
				else {
					// ������� � ��������� tecplot �����������.
					//� ���� �����.
					ANES_tecplot360_export_temperature(my_global_temperature_struct.maxnod, my_global_temperature_struct.pa, my_global_temperature_struct.maxelm, my_global_temperature_struct.nvtx, my_global_temperature_struct.potent, my_global_temperature_struct, f, 0, b, lb);
				}
			}

		}

		doublereal tmaxfinish = -273.15; // ���������� ����.
										 // ���������� �������� ������������ ����������� ������ ��������� ������� � �� � ��������:
										 //for (integer i = 0; i < t.maxelm + t.maxbound; i++) tmaxfinish = fmax(tmaxfinish, fabs(t.potent[i]));
										 // 23 ������� 2015
										 // �� ��������� ������ ���������� ����� �� ����� ��������� ������� �����������, �������
										 // �������� �� ������� ����� � ��������� ����������� ������ �� ���������� ��. 
		for (integer i = 0; i < my_global_temperature_struct.maxelm; i++) tmaxfinish = fmax(tmaxfinish, my_global_temperature_struct.potent[i]);

		doublereal totaldeform_max = -1.0e+30;
		if ((steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::STEADY_STATIC_STRUCTURAL) ||
			(steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::STEADY_STATIC_STRUCTURAL_AND_TEMPERATURE) ||
			(steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::UNSTEADY_STATIC_STRUCTURAL) ||
			(steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::UNSTEADY_STATIC_STRUCTURAL_AND_TEMPERATURE))
		{
			for (integer i = 0; i < my_global_temperature_struct.maxelm; i++) totaldeform_max = fmax(totaldeform_max, my_global_temperature_struct.total_deformation[TOTALDEFORMATION][i]);
		}

		FILE* fp = NULL;

#ifdef MINGW_COMPILLER
		int err1 = 0;
		fp = fopen64("report.txt", "w");
		if (fp == NULL) err1 = 1;
#else
		errno_t err1 = 0;
		err1 = fopen_s(&fp, "report.txt", "w");
#endif
		// �������� ����� ��� ������.
		if ((err1) != 0) {
			printf("Create File report.txt Error\n");
			//system("PAUSE");
			system("pause");
		}
		else {
			// ������ ���������
			fprintf(fp, "Maximum Temperature %.2f\n", tmaxfinish);
			fclose(fp);
		}
		// 1 - solver/solid_static/
		bool bMechanical = true;
		report_temperature(flow_interior, f, my_global_temperature_struct, b, lb, s, ls, w, lw, 0, matlist, bMechanical);

	}

	// steady Static Structural and Temperature (Thermal Stress).
	if (1 && steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::STEADY_STATIC_STRUCTURAL_AND_TEMPERATURE) {

		// ����� �������.
		unsigned int calculation_start_time = 0; // ������ ����� ��.
		unsigned int calculation_end_time = 0; // ��������� ����� ��.
		unsigned int calculation_seach_time = 0; // ����� ���������� ������� ���� � ��.

		calculation_start_time = clock(); // ������ ������ �����.

		for (integer i7 = 0; i7 < my_global_temperature_struct.maxelm + my_global_temperature_struct.maxbound; i7++)
			my_global_temperature_struct.potent[i7] = operating_temperature_for_film_coeff; // �������������.

																							// ������ �������������, � ����� Static Structural.
		bonly_solid_calculation = true;

		// �������� ����������� ���������� �� ����������� ������.
		if (lw == 1) {
			bPhysics_stop = true;
			if (lb < 11) {
				// ��� ����������� ��������:
				// MD40, AuSn, Cu, AuSn, SiC, GaN. cabinet and hollow.
				bPhysics_PTBSH_memory = true;
			}
		}

		if (adiabatic_vs_heat_transfer_coeff == DEFAULT_CABINET_BOUNDARY_CONDITION::NEWTON_RICHMAN_BC) {
			// �� �������������� ��������� ��������� ����� ��������� �� ��������������� ��� ������ ������ ����.
			//for (integer i7 = 0; i7<t.maxelm + t.maxbound; i7++) t.potent[i7] = 0.57*operating_temperature_for_film_coeff;
		}

		// ����� �������������� ��� �� ������ ������������ ������ ������ ����������������.
		bsolid_static_only = true;
		bool bcleantemp = false;
		if (eqin.itemper == 1) {
			bcleantemp = true;
			for (integer i = 0; i < flow_interior; i++) {
				if (eqin.fluidinfo[i].iflow == 1) bcleantemp = false;
			}
			// ���� bcleantemp==true �� �� ������ ������ ������ ������������� ��� ����� ���������.
		}
		if (1 || bcleantemp) {
			// ������� ������������ ���������� (��� ��������) ������ ������ ���������������� � ��������� �������. 
			printf("solution of pure Static Structural...\n");
			printf("please, press any key to continue...\n");
			if (bwait) {
				//system("PAUSE");
				system("pause");
			}

			// ��� ������������ ������������� ����������� ��������.
			bool bprintmessage = true; // �������� �� ��������� �� �������.

			doublereal dbeta = 1.0; // ������ ������� ������������� �� �������.
			bool bmyconvective = false;
			if (starting_speed_Vx * starting_speed_Vx + starting_speed_Vy * starting_speed_Vy + starting_speed_Vz * starting_speed_Vz > 1.0e-30) {
				if (f[0].maxelm > 0) {
					bmyconvective = true;
				}
			}
			else {
				// �������� ������������� ��������� ��������.

				FILE* fp_inicialization_data = NULL;
#ifdef MINGW_COMPILLER
				int err_inicialization_data = 0;
				fp_inicialization_data = fopen64("load.txt", "r");
				if (fp_inicialization_data == NULL) err_inicialization_data = 1;
#else
				errno_t err_inicialization_data = 0;
				err_inicialization_data = fopen_s(&fp_inicialization_data, "load.txt", "r");
#endif

				if (err_inicialization_data == 0) {
					// �������� ������ � ���� ������������.
					if (f[0].maxelm > 0) {
						bmyconvective = true;
					}
					fclose(fp_inicialization_data);
				}
			}

			// if (flow_interior>0) bmyconvective=true;
			// ������ ���������� ����������,
			// ��������� ��� �������� ������� ���-��� 1983
			doublereal** rhie_chow = nullptr;
			QuickMemVorst m;
			m.ballocCRSt = false; // �������� ������
			m.bsignalfreeCRSt = true; // � ����� �����������.

									  // ������������� ����������.
			m.tval = nullptr;
			m.tcol_ind = nullptr;
			m.trow_ptr = nullptr;
			m.tri = nullptr;
			m.troc = nullptr;
			m.ts = nullptr;
			m.tt = nullptr;
			m.tvi = nullptr;
			m.tpi = nullptr;
			m.tdx = nullptr;
			m.tdax = nullptr;
			m.ty = nullptr;
			m.tz = nullptr;
			m.ta = nullptr;
			m.tja = nullptr;
			m.tia = nullptr;
			m.talu = nullptr;
			m.tjlu = nullptr;
			m.tju = nullptr;
			m.tiw = nullptr;
			m.tlevs = nullptr;
			m.tw = nullptr;
			m.tjw = nullptr;
			m.icount_vel = 100000; // ����� ������� �����.



			integer* color = nullptr;
			integer dist_max = 3;
			calculate_color_for_temperature(color, my_global_temperature_struct, inx, xpos);


			// ���� flow_interior == 0 �� f[0] ������ ���������� ��������
			solve_nonlinear_temp(f[0], f, my_global_temperature_struct,
				rhie_chow,
				b, lb, s, ls, w, lw,
				dbeta, flow_interior,
				bmyconvective, nullptr, 0.001, 0.001,
				false,
				matlist, 0,
				bprintmessage,
				gtdps, ltdp, 1.0, 1.0, // ��������� �������� ������ 1.0 �������� ��� �������� �������.
				m, nullptr, // �������� � ����������� ���������� ����.
				nullptr,
				lu, my_union, color, dist_max); // �������� ����� ����� ������� � ����������� ���������� ����.

			delete[] color;



			// ����������� 19.05.2018

			doublereal* t_for_Mechanical = nullptr;
			//doublereal* lstub = nullptr;
			//integer maxelm_global_ret = 0;

			//solve_Thermal(my_global_temperature_struct, f, matlist, w, lw, lu, b, lb, ls, m, false,
			//operatingtemperature, false, 0.0, lstub, lstub,
			//maxelm_global_ret, 1.0, bAVLrealesation, t_for_Mechanical);

			t_for_Mechanical = new doublereal[my_global_temperature_struct.maxnod + 2];
			{
				// ����� ��������� �������.
				doublereal min_x = 1e60;
				doublereal min_y = 1e60;
				doublereal min_z = 1e60;
				doublereal max_x = -1e60;
				doublereal max_y = -1e60;
				doublereal max_z = -1e60;

				for (integer i = 0; i < my_global_temperature_struct.maxnod; i++) {
					if (my_global_temperature_struct.pa[i].x < min_x) {
						min_x = my_global_temperature_struct.pa[i].x;
					}
					if (my_global_temperature_struct.pa[i].y < min_y) {
						min_y = my_global_temperature_struct.pa[i].y;
					}
					if (my_global_temperature_struct.pa[i].z < min_z) {
						min_z = my_global_temperature_struct.pa[i].z;
					}
					if (my_global_temperature_struct.pa[i].x > max_x) {
						max_x = my_global_temperature_struct.pa[i].x;
					}
					if (my_global_temperature_struct.pa[i].y > max_y) {
						max_y = my_global_temperature_struct.pa[i].y;
					}
					if (my_global_temperature_struct.pa[i].z > max_z) {
						max_z = my_global_temperature_struct.pa[i].z;
					}
				}

				min_x = 1.05 * fabs(max_x - min_x);
				if (min_x < 1.0e-30) {
					min_x = 1.05 * fabs(max_x);
				}
				min_y = 1.05 * fabs(max_y - min_y);
				if (min_y < 1.0e-30) {
					min_y = 1.05 * fabs(max_y);
				}
				min_z = 1.05 * fabs(max_z - min_z);
				if (min_z < 1.0e-30) {
					min_z = 1.05 * fabs(max_z);
				}
				doublereal eps_mashine = 1.0e-308; // double

				doublereal* vol = new doublereal[my_global_temperature_struct.maxnod];

				for (integer i = 0; i < my_global_temperature_struct.maxnod; i++) {
					vol[i] = 0.0;
				}

				// �������������� ����������� � ����� ��� �� ����� ���.
				SECOND_ORDER_QUADRATIC_RECONSTRUCTA(my_global_temperature_struct.maxnod,
					my_global_temperature_struct.maxelm, my_global_temperature_struct.pa,
					my_global_temperature_struct.nvtx, vol, t_for_Mechanical, min_x, min_y, min_z, my_global_temperature_struct.potent,
					my_global_temperature_struct, eps_mashine, false, 1.0e-2);

				delete[] vol;

			}
			//bPhysics_stop = false;
			// ����� ������� Static Structural.
			doublereal* stub_Mechanical = nullptr;
			solve_Structural(my_global_temperature_struct, w, lw, true, operatingtemperature, b, lb, lu,
				false, 1.0e9, 1.0e9, 1.0e9, 1.0e9, stub_Mechanical, stub_Mechanical, stub_Mechanical, stub_Mechanical, 1.0,  matlist, t_for_Mechanical);
			//bPhysics_stop = true;

			delete[] t_for_Mechanical;

			// ���������� ����� ������.
			massa_cabinet(my_global_temperature_struct, f, flow_interior,
				b, lb, operatingtemperature,
				matlist);

			calculation_end_time = clock(); // ������ ��������� �����.
			calculation_seach_time = calculation_end_time - calculation_start_time;
			unsigned int im = 0, is = 0, ims = 0;
			im = (unsigned int)(calculation_seach_time / 60000); // ������
			is = (unsigned int)((calculation_seach_time - 60000 * im) / 1000); // �������
			ims = (unsigned int)((calculation_seach_time - 60000 * im - 1000 * is) / 10); // ������������ ������� �� 10

			printf("time calculation is:  %u minute %u second %u millisecond\n", im, is, 10 * ims);

			if (1) {
				if (!b_on_adaptive_local_refinement_mesh) {
					// ������� ���������� ���������� � ��������� tecplot360:
					exporttecplotxy360T_3D_part2(my_global_temperature_struct.maxelm, my_global_temperature_struct.ncell, f, my_global_temperature_struct, flow_interior, 0, bextendedprint, 0, b, lb);
				}
				else {
					// ������� � ��������� tecplot �����������.
					//� ���� �����.
					ANES_tecplot360_export_temperature(my_global_temperature_struct.maxnod, my_global_temperature_struct.pa, my_global_temperature_struct.maxelm, my_global_temperature_struct.nvtx, my_global_temperature_struct.potent, my_global_temperature_struct, f, 0, b, lb);
				}
			}

		}

		doublereal tmaxfinish = -273.15; // ���������� ����.
										 // ���������� �������� ������������ ����������� ������ ��������� ������� � �� � ��������:
										 //for (integer i = 0; i < t.maxelm + t.maxbound; i++) tmaxfinish = fmax(tmaxfinish, fabs(t.potent[i]));
										 // 23 ������� 2015
										 // �� ��������� ������ ���������� ����� �� ����� ��������� ������� �����������, �������
										 // �������� �� ������� ����� � ��������� ����������� ������ �� ���������� ��. 
		for (integer i = 0; i < my_global_temperature_struct.maxelm; i++) tmaxfinish = fmax(tmaxfinish, my_global_temperature_struct.potent[i]);

		doublereal totaldeform_max = -1.0e+30;
		if ((steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::STEADY_STATIC_STRUCTURAL) ||
			(steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::STEADY_STATIC_STRUCTURAL_AND_TEMPERATURE) ||
			(steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::UNSTEADY_STATIC_STRUCTURAL) ||
			(steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::UNSTEADY_STATIC_STRUCTURAL_AND_TEMPERATURE))
		{
			for (integer i = 0; i < my_global_temperature_struct.maxelm; i++) {
				totaldeform_max = fmax(totaldeform_max, my_global_temperature_struct.total_deformation[TOTALDEFORMATION][i]);
			}
		}

		FILE* fp = NULL;

#ifdef MINGW_COMPILLER
		int err1 = 0;
		fp = fopen64("report.txt", "w");
		if (fp == NULL) err1 = 1;
#else
		errno_t err1 = 0;
		err1 = fopen_s(&fp, "report.txt", "w");
#endif

		// �������� ����� ��� ������.
		if ((err1) != 0) {
			printf("Create File report.txt Error\n");
			//system("PAUSE");
			system("pause");
		}
		else {
			// ������ ���������
			fprintf(fp, "Maximum Temperature %.2f\n", tmaxfinish);
			fclose(fp);
		}
		// 1 - solver/solid_static/
		bool bMechanical = true;
		report_temperature(flow_interior, f, my_global_temperature_struct, b, lb, s, ls, w, lw, 0, matlist, bMechanical);
		// �������� ������ �������� ����� �������� ������� ������.
		// �������� ���������� (���������) �� ������� ����� �������� ����� � ��,
		// ���������� ����� �������� ������� ������. 28,10,2019
		report_out_boundary(f[0], my_global_temperature_struct, ls, lw, w, b, lb, matlist, f[0].OpTemp);


	}



	/*
	if (b_on_adaptive_local_refinement_mesh) {
	printf("t.maxbound=%d\n", t.maxbound);
	printf("v dvuch shagah ot ALICE sborki. \n");
	system("PAUSE");
	exit(1);  // ����� ����������� ������������ ��������� ���������� ���������� �������� ������������ �����.
	}

	if (b_on_adaptive_local_refinement_mesh) {
	printf("Solve temperature is compleate. \n");
	system("PAUSE");
	exit(1);  // ����� ����������� ������������ ��������� ���������� ���������� �������� ������������ �����.
	}
	*/

	//system("pause");

	// �������������� ����������������.
	if (1 && ((steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::UNSTEADY_TEMPERATURE) ||
		(steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::SECOND_TEMPERATURE_SOLVER) ||
		(steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::NETWORK_T_UNSTEADY) ||
		(steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::UNSTEADY_STATIC_STRUCTURAL) ||
		(steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::UNSTEADY_STATIC_STRUCTURAL_AND_TEMPERATURE))) {


		// ������ ������ ������������� � ������ ����.
		bonly_solid_calculation = true;

		// �������� ����������� ���������� �� ����������� ������.
		// ��� ���������� �������� ������ �� ���� �������� �������� ����������.
		if (lw == 1) {
			bPhysics_stop = true;
			if (lb < 11) {
				// ��� ����������� ��������:
				// MD40, AuSn, Cu, AuSn, SiC, GaN. cabinet and hollow.
				bPhysics_PTBSH_memory = true;
			}
		}

		bglobal_unsteady_temperature_determinant = true;
		// ����� �������������� ��� �� ������ ������ ������ ����������������.
		// ��������� 13 �������� bconvective
		// �� ����������� �������� ��������, �������� ������� ����������� ����������,
		// ����������� �� ����������� ��������� �����. ������� ������� ������������� �� 
		// ������� ������� ����������� ������ �����.
		// ��� ������������� ����� � ������������� ��������� 
		// ���������� ��������� ������� ������������� � ������� �� ������
		// ��� ��������� ����������� � ������ ������� �� 5 �������� �� ���� ��������� � 167 �������� � �������.
		// ��� ���������, �������� � icepak ����� 120 ��������, ��� ������ � ��������� ����������� Rt 
		// (�������� 6.875��, RT=16K/W) � �����������������
		// ������.
		doublereal dbeta = 1.3333333;//1.0; // ���� 1.0 �� ������ ������� ������������� �� �������.
		dbeta = 1.0; // ����� ���������� ��������.
					 // ������ ���������� ����������,
					 // ��������� ��� �������� ������� ���-��� 1983
		doublereal** rhie_chow = nullptr;

		// 29.06.2020
		if (steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::NETWORK_T_UNSTEADY) {
			calculate_Network_T_unsteady(my_global_temperature_struct, f,
				b, lb, w, lw, s, ls, matlist);
			//printf("unsteady temperature calculation is finished...\n");
			//system("PAUSE");
			//exit(1);
		}

		//std::cout << "steady_or_unsteady_global_determinant=" << steady_or_unsteady_global_determinant << std::endl;
		//getchar();

		//solve_nonlinear_temp(f[0], f, t, rhie_chow, b, lb, s, ls, w, lw, dbeta, flow_interior, false, nullptr, 0.001, false);
		bool bsecond_T_solver = false;
		if (steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::SECOND_TEMPERATURE_SOLVER) {
			// ������������� ������ �� ������ ���������� ������ ������� 10.11.2018.
			bsecond_T_solver = true;
		}

		if ((steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::UNSTEADY_TEMPERATURE) ||
			(steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::SECOND_TEMPERATURE_SOLVER) ||
			(steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::UNSTEADY_STATIC_STRUCTURAL) ||
			(steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::UNSTEADY_STATIC_STRUCTURAL_AND_TEMPERATURE)) {

			bool bMechanical = false; // true - ������ ������������ ���������� �������������� ���.
			bool bTemperature = true; // true - ������ �������������� �������������.
									  // bMechanical   && bTemperature   - ������ �������������� �������� � ������������� ���������.

			if (steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::UNSTEADY_STATIC_STRUCTURAL) {
				bMechanical = true;
				bTemperature = false;
			}
			if (steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::UNSTEADY_STATIC_STRUCTURAL_AND_TEMPERATURE) {
				bMechanical = true;
				bTemperature = true;
			}


			unsteady_temperature_calculation(f[0], f, my_global_temperature_struct,
				rhie_chow,
				b, lb, s, ls, w, lw,
				dbeta, flow_interior,
				matlist,
				operatingtemperature,
				gtdps, ltdp, lu, my_union, bsecond_T_solver, inx, xpos, bTemperature, bMechanical); // �������������� ������������� ������

		}

		// ��������� ������� ����������� � ������ �������.
		doublereal Tavg_fluid = 0.0;
		doublereal ic_avg_Temp_fluid = 0.0;
		for (integer i_7 = 0; i_7 < f[0].maxelm; i_7++) {
			if ((f[0].ptr[i_7] >= 0) && (f[0].ptr[i_7] < my_global_temperature_struct.maxelm)) {
				Tavg_fluid += my_global_temperature_struct.potent[f[0].ptr[i_7]];
				ic_avg_Temp_fluid += 1.0;
			}
		}
		if (ic_avg_Temp_fluid > 1.0e-30) {
			Tavg_fluid /= ic_avg_Temp_fluid;
			printf("average fluid temperature is %e\n", Tavg_fluid);
		}
		else {
			//printf("no fluid cell\n");
		}

		// ���������� ����� ������.
		massa_cabinet(my_global_temperature_struct, f, flow_interior,
			b, lb, operatingtemperature,
			matlist);

		// 10.10.2017
		// ���������� ���������� ������� ����������� ����� ���������
		// �� �������� ������ �������.
		xyplot_temp(my_global_temperature_struct, my_global_temperature_struct.potent);

		if (!bsecond_T_solver) {
			if (!b_on_adaptive_local_refinement_mesh) {
				// ������� ���������� ���������� � ��������� tecplot360:
				exporttecplotxy360T_3D_part2(my_global_temperature_struct.maxelm, my_global_temperature_struct.ncell, f, my_global_temperature_struct, flow_interior, 0, bextendedprint, 0, b, lb);
			}
			else {
				// ������� � ��������� tecplot �����������.
				//� ���� �����.
				ANES_tecplot360_export_temperature(my_global_temperature_struct.maxnod, my_global_temperature_struct.pa, my_global_temperature_struct.maxelm, my_global_temperature_struct.nvtx, my_global_temperature_struct.potent, my_global_temperature_struct, f, 0, b, lb);
			}


			doublereal tmaxfinish = -273.15;
			// ���������� �������� ������������ ����������� ������ ��������� ������� � �� � ��������:
			for (integer i = 0; i < my_global_temperature_struct.maxelm + my_global_temperature_struct.maxbound; i++) {
				tmaxfinish = fmax(tmaxfinish, my_global_temperature_struct.potent[i]);
			}
			FILE* fp = NULL;

#ifdef MINGW_COMPILLER
			int err1 = 0;
			fp = fopen64("report.txt", "w");
			if (fp == NULL) err1 = 1;
#else
			errno_t err1 = 0;
			err1 = fopen_s(&fp, "report.txt", "w");
#endif
			// �������� ����� ��� ������.
			if ((err1) != 0) {
				printf("Create File report.txt Error\n");
				//system("PAUSE");
				system("pause");
			}
			else {
				// ������ ���������
				fprintf(fp, "Maximum Temperature %.2f\n", tmaxfinish);
				fclose(fp);
			}
			// 1 - solver/solid_static/
			bool bMechanical = true;
			report_temperature(flow_interior, f, my_global_temperature_struct, b, lb, s, ls, w, lw, 0, matlist, bMechanical);
			// �������� ������ �������� ����� �������� ������� ������.
			// �������� ���������� (���������) �� ������� ����� �������� ����� � ��,
			// ���������� ����� �������� ������� ������. 28,10,2019
			report_out_boundary(f[0], my_global_temperature_struct, ls, lw, w, b, lb, matlist, f[0].OpTemp);
		}
		else {
			printf("THIS IS SECOND UNSTEADY TEMPERATURE SOLVER ON ALL MESHES.\n");
			printf("NO EXPOPRT TECPLOT.\n");
			printf("NO PRINT REPORT.\n");
		}
		printf("calculation complete...\n");
		// system("PAUSE");
	}



	// ������� ���������� ���������� � ��������� tecplot360:
	// ����� ������������ ��� �������� ����������� �����.
	if (false) {
		exporttecplotxy360T_3D_part2(my_global_temperature_struct.maxelm,
			my_global_temperature_struct.ncell,
			f, my_global_temperature_struct,
			flow_interior, 0, bextendedprint,
			0, b, lb);
		printf("read values. OK.\n");
		if (bwait) {
			//system("PAUSE"); // debug avtosave
			system("pause");
		}
	}

	// Fluid dynamic ������������ �������������� �������.
	if ((1 && steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::CFD_STEADY)) {

#ifdef _OPENMP
		// �������������� � ������������� ������� cfd �� openMP 08.05.2019.
		printf("CFD not work in OPENMP ON and bparallelismold is true.\n");
		printf("uskorenie ot OPENMP otsutstvuet. Rabotaen odnopotochnaq versiq.\n");
		printf("variable bparallelismold must be equal false.\n");
		//system("PAUSE");
#endif

		told_temperature_global_for_HOrelax = new doublereal[my_global_temperature_struct.maxelm + my_global_temperature_struct.maxbound];
		bSIMPLErun_now_for_temperature = true;

		if (dgx * dgx + dgy * dgy + dgz * dgz > 1.0e-20) {
			// ���� ����� ��������� �������� �� ��� fluid ���������� ����������� ����������.
			bool bbussinesk_7 = false;
#pragma omp parallel for 
			for (integer i_8 = 0; i_8 < f[0].maxelm; i_8++) {
				integer ib = my_global_temperature_struct.whot_is_block[f[0].ptr[i_8]];
				if (ib > -1) {
					if (b[ib].itype == PHYSICS_TYPE_IN_BODY::FLUID) {
						integer i_7 = b[ib].imatid;
						if (matlist[i_7].bBussineskApproach) {
#pragma omp critical
							{
								if (bbussinesk_7 == false) {
									bbussinesk_7 = true;
								}
							}
						}
					}
				}
			}
			if (bbussinesk_7) {
				bSIMPLErun_now_for_natural_convection = true;
			}
		}

		const integer ISIZEf = f[0].maxelm + f[0].maxbound;

		bHORF = true;
		bPamendment_source_old = new doublereal[ISIZEf];
#pragma omp parallel for
		for (integer i5 = 0; i5 < ISIZEf; i5++) {
			bPamendment_source_old[i5] = 0.0;
		}
		// exporttecplotxy360T_3D_part2(t.maxelm, t.ncell, f, t, flow_interior, 0, bextendedprint);
		//system("PAUSE");
		if (dgx * dgx + dgy * dgy + dgz * dgz > 1.0e-20) {
			// ���� ����� ��������� �������� �� ��� fluid ���������� ����������� ����������.
			bool bbussinesk_7 = false;
#pragma omp parallel for 
			for (integer i_8 = 0; i_8 < f[0].maxelm; i_8++) {
				integer ib = my_global_temperature_struct.whot_is_block[f[0].ptr[i_8]];
				if (ib > -1) {
					if (b[ib].itype == PHYSICS_TYPE_IN_BODY::FLUID) {
						integer i_7 = b[ib].imatid;
						if (matlist[i_7].bBussineskApproach) {
#pragma omp critical 
							{
								if (bbussinesk_7 == false) {
									bbussinesk_7 = true;
								}
							}
						}
					}
				}
			}
			if (bbussinesk_7) {
				printf("Bussinesk approach Operating Temperature=%e\n", f[0].OpTemp); // Operating Temperature);
			}
		}

		//exporttecplotxy360T_3D_part2(t.maxelm, t.ncell, f, t, flow_interior, 0, bextendedprint, 0);
		//system("PAUSE");

		// ������������ ����������������� ��������.
		steady_cfd_calculation(breadOk,
			eqin, dgx, dgy, dgz,
			continity_start,
			inumber_iteration_SIMPLE,
			flow_interior, f, my_global_temperature_struct, b, lb,
			s, ls, w, lw, matlist,
			gtdps, ltdp, bextendedprint, lu, my_union, inx, xpos);
		// xyplot( f, 0, t);
		// boundarylayer_info(f, t, flow_interior, w, lw);
		// 2 - solver/conjugate_heat_transfer_static/
		bool bMechanical = false;

		report_temperature(flow_interior, f, my_global_temperature_struct, b, lb, s, ls, w, lw, 0/*2*/, matlist, bMechanical);
		// �������� ������ �������� ����� �������� ������� ������.
		// �������� ���������� (���������) �� ������� ����� �������� ����� � ��,
		// ���������� ����� �������� ������� ������. 28,10,2019
		report_out_boundary(f[0], my_global_temperature_struct, ls, lw, w, b, lb, matlist, f[0].OpTemp);

		// ���������� ����� ������.
		massa_cabinet(my_global_temperature_struct, f, flow_interior,
			b, lb, operatingtemperature,
			matlist);

		// ���������� ���������� ������� ����������������� �������.
		// ������� ������� ����� �������� �������� ������ ��������� 
		// �� ������������ ����������. AliceMesh v.0.45.
		xyplot(f, flow_interior, my_global_temperature_struct);

		// 10.10.2017
		// ���������� ���������� ������� ����������� ����� ���������
		// �� �������� ������ �������.
		xyplot_temp(my_global_temperature_struct, my_global_temperature_struct.potent);

		// ������� ���������� ���������� � ��������� tecplot360:
		if (!b_on_adaptive_local_refinement_mesh) {
			exporttecplotxy360T_3D_part2(my_global_temperature_struct.maxelm, my_global_temperature_struct.ncell, f, my_global_temperature_struct, flow_interior, 0, bextendedprint, 0, b, lb);
		}
		else {
			ANES_tecplot360_export_temperature(my_global_temperature_struct.maxnod, my_global_temperature_struct.pa, my_global_temperature_struct.maxelm, my_global_temperature_struct.nvtx, my_global_temperature_struct.potent, my_global_temperature_struct, f, 0, b, lb);
		}

		save_velocity_for_init(my_global_temperature_struct.maxelm, my_global_temperature_struct.ncell, f, my_global_temperature_struct, flow_interior);

		// exporttecplotxy360T_3D_part2_rev(t.maxelm, t.ncell, f, t, flow_interior, 0, bextendedprint,b,lb);
		delete[] bPamendment_source_old;
		bPamendment_source_old = nullptr;
		delete[] told_temperature_global_for_HOrelax;
		told_temperature_global_for_HOrelax = nullptr;
	}


	// �������������� ����������������� ��������: 
	if ((1 && (steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::CFD_UNSTEADY))) {


#ifdef _OPENMP
		// �������������� � ������������� ������� cfd �� openMP 08.05.2019.
		//printf("Unsteady CFD not work in OPENMP ON and bparallelismold is true.\n");
		//printf("uskorenie ot OPENMP otsutstvuet. Rabotaen odnopotochnaq versiq.\n");
		//printf("variable bparallelismold must be equal false.\n");
		//system("PAUSE");
#endif

		told_temperature_global_for_HOrelax = new doublereal[my_global_temperature_struct.maxelm + my_global_temperature_struct.maxbound];
		bSIMPLErun_now_for_temperature = true;


		if (dgx * dgx + dgy * dgy + dgz * dgz > 1.0e-20) {
			// ���� ����� ��������� �������� �� ��� fluid ���������� ����������� ����������.
			bool bbussinesk_7 = false;
#pragma omp parallel for
			for (integer i_8 = 0; i_8 < f[0].maxelm; i_8++) {
				integer ib = my_global_temperature_struct.whot_is_block[f[0].ptr[i_8]];
				if (ib > -1) {
					if (b[ib].itype == PHYSICS_TYPE_IN_BODY::FLUID) {
						integer i_7 = b[ib].imatid;
						if (matlist[i_7].bBussineskApproach) {
#pragma omp critical
							{
								if (bbussinesk_7 == false) {
									bbussinesk_7 = true;
								}
							}
						}
					}
				}
			}
			if (bbussinesk_7) {
				bSIMPLErun_now_for_natural_convection = true;
			}
		}

		const integer ISIZEf = f[0].maxelm + f[0].maxbound;
		bHORF = true;
		bPamendment_source_old = new doublereal[ISIZEf];
#pragma omp parallel for
		for (integer i5 = 0; i5 < ISIZEf; i5++) {
			bPamendment_source_old[i5] = 0.0;
		}

		// exporttecplotxy360T_3D_part2(t.maxelm, t.ncell, f, t, flow_interior, 0, bextendedprint);
		//system("PAUSE");
		if (dgx * dgx + dgy * dgy + dgz * dgz > 1.0e-20) {
			// ���� ����� ��������� �������� �� ��� fluid ���������� ����������� ����������.
			bool bbussinesk_7 = false;
#pragma omp parallel for 
			for (integer i_8 = 0; i_8 < f[0].maxelm; i_8++) {
				integer ib = my_global_temperature_struct.whot_is_block[f[0].ptr[i_8]];
				if (ib > -1) {
					if (b[ib].itype == PHYSICS_TYPE_IN_BODY::FLUID) {
						integer i_7 = b[ib].imatid;
						if (matlist[i_7].bBussineskApproach) {
							if (bbussinesk_7 == false)
							{
#pragma omp critical
								{
									bbussinesk_7 = true;
								}
							}
						}
					}
				}
			}
			if (bbussinesk_7) {
				printf("Bussinesk approach Operating Temperature=%e\n", f[0].OpTemp); // Operating Temperature);
			}
		}

		//exporttecplotxy360T_3D_part2(t.maxelm, t.ncell, f, t, flow_interior, 0, bextendedprint, 0);
		//system("PAUSE");

		if (0) {
			printf("Sorry unsteady cfd calcuation dont support... 21.07.2019\n");
			printf("Your may send your message to kirill7785@mail.ru.\n");
			system("pause");
		}
		else {
			// �������������� ����������������� ��������:
			usteady_cfd_calculation(breadOk, eqin,
				dgx, dgy, dgz,
				continity_start,
				inumber_iteration_SIMPLE,
				flow_interior,
				f, my_global_temperature_struct,
				b, lb, s, ls,
				w, lw, matlist, gtdps, ltdp, bextendedprint, lu, my_union, inx, xpos);
		}

		//xyplot( f, 0, t);
		// boundarylayer_info(f, t, flow_interior, w, lw);
		// 2 - solver/conjugate_heat_transfer_static/
		bool bMechanical = false;
		report_temperature(flow_interior, f, my_global_temperature_struct, b, lb, s, ls, w, lw, 0/*2*/, matlist, bMechanical);

		// ���������� ����� ������.
		massa_cabinet(my_global_temperature_struct, f, flow_interior,
			b, lb, operatingtemperature,
			matlist);

		// ���������� ���������� ������� ����������������� �������.
		// ������� ������� ����� �������� �������� ������ ��������� 
		// �� ������������ ����������. AliceMesh v.0.45.
		xyplot(f, flow_interior, my_global_temperature_struct);

		// 10.10.2017
		// ���������� ���������� ������� ����������� ����� ���������
		// �� �������� ������ �������.
		xyplot_temp(my_global_temperature_struct, my_global_temperature_struct.potent);

		// ������� ���������� ���������� � ��������� tecplot360:
		if (!b_on_adaptive_local_refinement_mesh) {
			exporttecplotxy360T_3D_part2(my_global_temperature_struct.maxelm, my_global_temperature_struct.ncell, f, my_global_temperature_struct, flow_interior, 0, bextendedprint, 0, b, lb);
		}
		else {
			ANES_tecplot360_export_temperature(my_global_temperature_struct.maxnod, my_global_temperature_struct.pa, my_global_temperature_struct.maxelm, my_global_temperature_struct.nvtx, my_global_temperature_struct.potent, my_global_temperature_struct, f, 0, b, lb);
		}

		save_velocity_for_init(my_global_temperature_struct.maxelm, my_global_temperature_struct.ncell, f, my_global_temperature_struct, flow_interior);
		// exporttecplotxy360T_3D_part2_rev(t.maxelm, t.ncell, f, t, flow_interior, 0, bextendedprint,b,lb);

		delete[] bPamendment_source_old;
		bPamendment_source_old = nullptr;
		delete[] told_temperature_global_for_HOrelax;
		told_temperature_global_for_HOrelax = nullptr;

	}


	if (continity_start != nullptr) {
		delete[] continity_start;
		continity_start = nullptr;
	}

	if (inumber_iteration_SIMPLE != nullptr) {
		delete[] inumber_iteration_SIMPLE;
		inumber_iteration_SIMPLE = nullptr;
	}

	
	// ������������ ����������� ������.
	if (xpos != nullptr) {
		delete[] xpos;
		xpos = nullptr;
	}
	if (ypos != nullptr) {
		delete[] ypos;
		ypos = nullptr;
	}
	if (zpos != nullptr) {
		delete[] zpos;
		zpos = nullptr;
	}

	if (bsource_term_radiation_for_relax != nullptr) {
		delete[] bsource_term_radiation_for_relax; // ���������� ������������ ������ ������������ �������.
		bsource_term_radiation_for_relax = nullptr;
	}
	if (b_buffer_correct_source != nullptr) {
		delete[] b_buffer_correct_source;
		b_buffer_correct_source = nullptr;
	}

	printf("free memory begin...\n");
	if (bwait) {
		//system("PAUSE");
		system("pause");
	}

	if (rthdsd_no_radiosity_patch != nullptr) {
		delete[] rthdsd_no_radiosity_patch;
		rthdsd_no_radiosity_patch = nullptr;
	}


	// ������� ��������� ���������� ��������� ������� � ����� v. 0.14 ��������.
	if (qnbc != nullptr) {
		delete[] qnbc;
		qnbc = nullptr;
		iadd_qnbc_maxelm = 0;
	}

	

	/*
	for (integer i_7 = 0; i_7 < lb; i_7++) {
	if (b[i_7].temp_Sc != nullptr) {
	delete[] b[i_7].temp_Sc;
	b[i_7].temp_Sc = nullptr;
	}
	if (b[i_7].arr_Sc != nullptr) {
	delete[] b[i_7].arr_Sc;
	b[i_7].arr_Sc = nullptr;
	}
	if (b[i_7].g.hi != nullptr) {
	delete[] b[i_7].g.hi;
	b[i_7].g.hi = nullptr;
	}
	if (b[i_7].g.xi != nullptr) {
	delete[] b[i_7].g.xi;
	b[i_7].g.xi = nullptr;
	}
	if (b[i_7].g.yi != nullptr) {
	delete[] b[i_7].g.yi;
	b[i_7].g.yi = nullptr;
	}
	if (b[i_7].g.zi != nullptr) {
	delete[] b[i_7].g.zi;
	b[i_7].g.zi = nullptr;
	}
	}
	delete[] b;
	b = nullptr;
	*/

	delete[] s; delete[] w; // ������������ ������
	s = nullptr;
	w = nullptr;
	
	delete[] gtdps;
	gtdps = nullptr;
	if (eqin.fluidinfo != nullptr) {
		delete[] eqin.fluidinfo;
		eqin.fluidinfo = nullptr;
	}
	// ����� ���������� ����������� ������ �� ��� ���� �������� ������:
	free_level1_temp(my_global_temperature_struct);
	free_level2_temp(my_global_temperature_struct); // ������������ ������ �� ��� ������.
													// ����������� ������ ��� LR ������.
	free_root(my_global_temperature_struct.rootWE, my_global_temperature_struct.iWE);
	free_root(my_global_temperature_struct.rootSN, my_global_temperature_struct.iSN);
	free_root(my_global_temperature_struct.rootBT, my_global_temperature_struct.iBT);
	if (my_global_temperature_struct.rootWE != nullptr) {
		delete[] my_global_temperature_struct.rootWE;
		my_global_temperature_struct.rootWE = nullptr;
	}
	if (my_global_temperature_struct.rootSN != nullptr) {
		delete[] my_global_temperature_struct.rootSN;
		my_global_temperature_struct.rootSN = nullptr;
	}
	if (my_global_temperature_struct.rootBT != nullptr) {
		delete[] my_global_temperature_struct.rootBT;
		my_global_temperature_struct.rootBT = nullptr;
	}
	// ������������ ������ ��� LR �����.
	free_level1_flow(f, flow_interior);
	free_level2_flow(f, flow_interior); // ������������ ������ �� ��� ������.

	delete[] f;
	f = nullptr;

	if (sourse2Dproblem != nullptr) {
		delete[] sourse2Dproblem;
		sourse2Dproblem = nullptr;
	}
	if (conductivity2Dinsource != nullptr) {
		delete[] conductivity2Dinsource;
		conductivity2Dinsource = nullptr;
	}

	free_global_1(matlist, lmatmax, x_jacoby_buffer, bvery_big_memory, my_global_temperature_struct, amgGM, xposadd, yposadd, zposadd);

	

	for (integer i63 = 0; i63 < lu; i63++) {
		// ����� ���������� ����������� ������ �� ��� ���� �������� ������:
		free_level1_temp(my_union[i63].t);
		free_level2_temp(my_union[i63].t); // ������������ ������ �� ��� ������.

										   // ����������� ������ ��� LR ������.
		free_root(my_union[i63].t.rootWE, my_union[i63].t.iWE);
		free_root(my_union[i63].t.rootSN, my_union[i63].t.iSN);
		free_root(my_union[i63].t.rootBT, my_union[i63].t.iBT);
		if (my_union[i63].t.rootWE != nullptr) {
			delete[] my_union[i63].t.rootWE;
			my_union[i63].t.rootWE = nullptr;
		}
		if (my_union[i63].t.rootSN != nullptr) {
			delete[] my_union[i63].t.rootSN;
			my_union[i63].t.rootSN = nullptr;
		}
		if (my_union[i63].t.rootBT != nullptr) {
			delete[] my_union[i63].t.rootBT;
			my_union[i63].t.rootBT = nullptr;
		}

		// ������������ ������ ��� LR �����.
		free_level1_flow(my_union[i63].f, my_union[i63].flow_interior);
		free_level2_flow(my_union[i63].f, my_union[i63].flow_interior); // ������������ ������ �� ��� ������.

		delete[] my_union[i63].f;
		my_union[i63].f = nullptr;

		if (bvery_big_memory) {
			if (my_union[i63].t.database.x != nullptr) {
				free(my_union[i63].t.database.x);
				my_union[i63].t.database.x = nullptr;
			}
			if (my_union[i63].t.database.y != nullptr) {
				free(my_union[i63].t.database.y);
				my_union[i63].t.database.y = nullptr;
			}
			if (my_union[i63].t.database.z != nullptr) {
				free(my_union[i63].t.database.z);
				my_union[i63].t.database.z = nullptr;
			}
			if (my_union[i63].t.database.nvtxcell != nullptr) {
				for (integer i = 0; i <= 7; i++) {
					delete[] my_union[i63].t.database.nvtxcell[i];
				}
				delete[] my_union[i63].t.database.nvtxcell;
				my_union[i63].t.database.nvtxcell = nullptr;
			}
			if (my_union[i63].t.database.ptr != nullptr) {
				if (my_union[i63].t.database.ptr[0] != nullptr) {
					delete[] my_union[i63].t.database.ptr[0];
				}
				if (my_union[i63].t.database.ptr[1] != nullptr) {
					delete[] my_union[i63].t.database.ptr[1];
				}
				delete[] my_union[i63].t.database.ptr;
				my_union[i63].t.database.ptr = nullptr;
			}
		}
	}

	// ������������ ����� ������ � ILU ������.
	if (milu_gl_buffer.alu_copy != nullptr) delete[] milu_gl_buffer.alu_copy;
	if (milu_gl_buffer.jlu_copy != nullptr) delete[] milu_gl_buffer.jlu_copy;
	if (milu_gl_buffer.ju_copy != nullptr) delete[] milu_gl_buffer.ju_copy;
	milu_gl_buffer.alu_copy = nullptr;
	milu_gl_buffer.jlu_copy = nullptr;
	milu_gl_buffer.ju_copy = nullptr;

	free_QSBid(); // ��� ���������� ������ ������� myisblock_id.

	flow_interior = 0;
	printf("free memory finish...\n");

	if (1 && steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::MESHER_ONLY) {
		// ��� ������ ����� ��������� ����������.
		printf("Mesh generation procedure is finish.\n");
	}
	else {
		printf("Calculation procedure is finish.\n");
	}

	if (bwait) {
		printf("Please, press any key to exit...\n");
		//system("PAUSE");
		system("pause");
	}

	// ������������ ����������� ������ �� ��� ������ ���������
	// ���� �� �������.
	if (glTSL.table_law_piecewise_constant != nullptr) {
		delete[] glTSL.table_law_piecewise_constant;
		glTSL.table_law_piecewise_constant = nullptr;
	}

	// ��� ��� � ������ bFULL_AUTOMATIC ������� ������������ �������� �
	// ������� ������������ �������, �� �������� ������� ����������� ���� ���� ���, �
	// ��� ��������� ��������� ���� ��������� � ������ ���-�������.
	// 20mm ���� ��������� � 1��� 9� �� 53� �� ���� ������ bFULL_AUTOMATIC.
	// ���-������� ��� automatic
	// ������������ ����������� ������ �� ��� ���-�������.
	if (shorter_hash_X != nullptr) {
		delete[] shorter_hash_X;
		shorter_hash_X = nullptr;
	}
	if (shorter_hash_Y != nullptr) {
		delete[] shorter_hash_Y;
		shorter_hash_Y = nullptr;
	}
	if (shorter_hash_Z != nullptr) {
		delete[] shorter_hash_Z;
		shorter_hash_Z = nullptr;
	}
	if (bshorter_hash_X != nullptr) {
		delete[] bshorter_hash_X;
		bshorter_hash_X = nullptr;
	}
	if (bshorter_hash_Y != nullptr) {
		delete[] bshorter_hash_Y;
		bshorter_hash_Y = nullptr;
	}
	if (bshorter_hash_Z != nullptr) {
		delete[] bshorter_hash_Z;
		bshorter_hash_Z = nullptr;
	}

	// ����������� ������ �� ��� ������� �����.
	//delete[] StringList;

	calculation_main_end_time = clock();
	calculation_main_seach_time = calculation_main_end_time - calculation_main_start_time_global_Depend;


	/*printf("time=%d statistic vorst=%3.2f %% \n",calculation_main_seach_time,(float)(100.0*calculation_vorst_seach_time/calculation_main_seach_time));
	system("PAUSE");
	*/



	// ����� ����� ����������.
	int im = 0, is = 0, ims = 0;
	im = (int)(calculation_main_seach_time / 60000); // ������
	is = (int)((calculation_main_seach_time - 60000 * im) / 1000); // �������
	ims = (int)((calculation_main_seach_time - 60000 * im - 1000 * is) / 10); // ������������ ������� �� 10

	printf("time calculation is:  %d minute %d second %d millisecond\n", im, is, 10 * ims);

	if (1 && (steady_or_unsteady_global_determinant != PHYSICAL_MODEL_SWITCH::PREOBRAZOVATEL_FOR_REPORT)) {
		//system("pause");
		//46.321
		//31.655
	}

	return calculation_main_seach_time;
} // main_body

int main_solverv0_48(BLOCK*& b, integer& lb)
{

	// ��������� ��������� � ������� Windows
	// ������ ����������� ����� � ������� �����.
	//setlocale(LC_ALL, "");
	system("mode con cols=166 lines=12000");

	char ch_EXPORT_ALICE_ONLY = 'y';

	if (0) {
		// �����������. ������� �� ������ ��� ����������:
		// d_my_optimetric1_6_12_2019;
		// d_my_optimetric2_6_12_2019;
		// d_my_optimetric3_6_12_2019.
		// � ���������� �������� ���������.
		// ��� ����������� ��� ���������� ��������� �������� � ��� � 
		// � ������ ��������� �� ���������.
		// ��� ����������
		//int in_optimetric = 21;
		int* time_optimetric = new int[1000];
		doublereal* op_comp = new doublereal[1000];
		int iscan_optimetric = 0;
		for (d_my_optimetric1_6_12_2019 = 0.33; d_my_optimetric1_6_12_2019 < 0.7; d_my_optimetric1_6_12_2019 += 0.01) {
			//for (d_my_optimetric2_6_12_2019 = 0.09; d_my_optimetric2_6_12_2019 < 0.95; d_my_optimetric2_6_12_2019 += 0.005) {
			//for (d_my_optimetric3_6_12_2019 = 0.11; d_my_optimetric3_6_12_2019 < 0.12; d_my_optimetric3_6_12_2019 += 0.01) {
			//printf("%e %e %e %d\n", d_my_optimetric1_6_12_2019, d_my_optimetric2_6_12_2019, d_my_optimetric3_6_12_2019, time_optimetric[iscan_optimetric]);
			//printf("%e %d\n", d_my_optimetric1_6_12_2019, time_optimetric[iscan_optimetric]);
			time_optimetric[iscan_optimetric] = main_body(b, lb, ch_EXPORT_ALICE_ONLY);
			op_comp[iscan_optimetric] = d_my_optimetric2_6_12_2019;
			iscan_optimetric++;

			//}
			//}
		}

		printf("\n\n\n");
		iscan_optimetric = 0;
		for (d_my_optimetric1_6_12_2019 = 0.33; d_my_optimetric1_6_12_2019 < 0.7; d_my_optimetric1_6_12_2019 += 0.01) {
			//for (d_my_optimetric2_6_12_2019 = 0.09; d_my_optimetric2_6_12_2019 < 0.95; d_my_optimetric2_6_12_2019 += 0.005) {
			//for (d_my_optimetric3_6_12_2019 = 0.11; d_my_optimetric3_6_12_2019 < 0.12; d_my_optimetric3_6_12_2019 += 0.01) {
			//printf("%e %e %e %d\n", d_my_optimetric1_6_12_2019, d_my_optimetric2_6_12_2019, d_my_optimetric3_6_12_2019, time_optimetric[iscan_optimetric]);
			printf("%e CA=%e %d\n", d_my_optimetric1_6_12_2019, op_comp[iscan_optimetric], time_optimetric[iscan_optimetric]);
			iscan_optimetric++;
			//}
			//}
		}
		delete[] time_optimetric;
		delete[] op_comp;
	}
	if (0) {
		// �����������. ������� �� ������ ����� ����������.
		// d_my_optimetric1_6_12_2019.
		// � ���������� �������� ���������.
		// ��� ����������� ��� ���������� ��������� �������� � ��� � 
		// � ������ �������� � ���������.
		// ���� ����������

		int* time_optimetric = new int[300];
		int iscan_optimetric = 0;
		for (d_my_optimetric1_6_12_2019 = 0.65; d_my_optimetric1_6_12_2019 < 0.8; d_my_optimetric1_6_12_2019 += 0.01) {
			printf("%e %d\n", d_my_optimetric1_6_12_2019, time_optimetric[iscan_optimetric]);
			time_optimetric[iscan_optimetric] = main_body(b, lb, ch_EXPORT_ALICE_ONLY);
			iscan_optimetric++;
		}

		printf("\n\n\n");
		iscan_optimetric = 0;
		for (d_my_optimetric1_6_12_2019 = 0.65; d_my_optimetric1_6_12_2019 < 0.8; d_my_optimetric1_6_12_2019 += 0.01) {
			printf("%e %d\n", d_my_optimetric1_6_12_2019, time_optimetric[iscan_optimetric]);
			iscan_optimetric++;
		}
		delete[] time_optimetric;
	}
	if (1) {
		main_body(b, lb, ch_EXPORT_ALICE_ONLY);
	}

	// ����������� ������ �� ��� ������� �����.
	delete[] StringList;


	system("pause");
	return 0;
}

#ifndef NO_OPENGL_GLFW

void DrawPrism(BLOCK*& b, integer lb, integer id, GLfloat centerPosX, GLfloat centerPosY, GLfloat  centerPosZ, GLfloat scale);
void DrawCylinder(BLOCK*& b, integer lb, integer id, GLfloat centerPosX, GLfloat centerPosY, GLfloat  centerPosZ, GLfloat scale);
void DrawPolygon(BLOCK*& b, integer lb, integer id, GLfloat centerPosX, GLfloat centerPosY, GLfloat  centerPosZ, GLfloat scale);
void DrawCADobj(BLOCK*& b, integer lb, integer id, GLfloat centerPosX, GLfloat centerPosY, GLfloat  centerPosZ, GLfloat scale);
void DrawKarkas(BLOCK*& b, integer lb, GLfloat centerPosX, GLfloat centerPosY, GLfloat  centerPosZ, GLfloat scale) {

	bool b_flag = true;
	for (integer i = 0; i < lb; i++) {
		if (b_flag) {
			if (binverse_color_black_2_white) {
				glColor3f(0.0, 0.0, 0.0); // ������
			}
			else {
				glColor3f(1.0, 1.0, 1.0); // �����
			}
			b_flag = false;
		}
		if (b[i].g.itypegeom == PRISM) {
			DrawPrism(b, lb, i, centerPosX, centerPosY, centerPosZ, scale);
		}
		if (b[i].g.itypegeom == CYLINDER) {
			DrawCylinder(b, lb, i, centerPosX, centerPosY, centerPosZ, scale);
		}
		if (b[i].g.itypegeom == POLYGON) {
			DrawPolygon(b, lb, i, centerPosX, centerPosY, centerPosZ, scale);
		}
		if (b[i].g.itypegeom == CAD_STL)
		{
			glColor3f(0.502, 0.502, 0.502); // ����� ��� CAD_STL obj
			DrawCADobj(b, lb, i, centerPosX, centerPosY, centerPosZ, scale);
			b_flag = true;
		}
	}

}

#endif

int main(void)
{

	

	int icol_now = 0;
#ifndef NO_OPENGL_GLFW
	mask_color_for_openGL = new int[1021];
	for (int i = 0; i < 1021; i++) {
		mask_color_for_openGL[i] = icol_now;
		if (i % 72 == 0) icol_now = i;// ����� 14 ��������� ������.
	}
#endif

	// ���������� ������, ���������� � ������, �������.
	integer lb = 0;
	BLOCK* b = nullptr;// ������ ������

					   // �������� ������.
	main_solverv0_48(b, lb);

#ifndef NO_OPENGL_GLFW

	GLFWwindow* window;

	/* Initialize the library */
	if (!glfwInit())
		return -1;

	/* Create a windowed mode window and its OpenGL context */
	window = glfwCreateWindow(SCREEN_WIDTH, SCREEN_HEIGHT, "AliceFlow_v0_49", NULL, NULL);


	glfwSetKeyCallback(window, keyCallback);
	glfwSetInputMode(window, GLFW_STICKY_KEYS, 1);

	int screenWidth, screenHeight;
	glfwGetFramebufferSize(window, &screenWidth, &screenHeight);

	if (!window)
	{
		glfwTerminate();
		return -1;
	}

	/* Make the window's context current */
	glfwMakeContextCurrent(window);


	glViewport(0.0f, 0.0f, screenWidth, screenHeight); // specifies the part of the window to which OpenGL will draw (in pixels), convert from
													   // normalised to pixels
	glMatrixMode(GL_PROJECTION); // projection matrix defines the properties of the camera that views the objects in the world coordinate frame.
								 // Here you thalfScreenHeight+scale_all*pa.yally set the zoom factor, aspect ratio and the near and far clipping planes
	glLoadIdentity(); // replace the current matrix with the identity matrix and starts us a fresh because matrix transforms such as glOrtho and
					  // glRotate cumulate, besically puts us at (0, 0, 0)
	glOrtho(0, SCREEN_WIDTH, 0, SCREEN_HEIGHT, 0, 1000); // essentially set coordinate system
	glMatrixMode(GL_MODELVIEW); // (default matrix mode) modelview matrix defines how objects are transformed (meaning translation, rotation
								// and scaling) in your world

	glLoadIdentity(); // same as above comment 

					  //GLfloat halfScreenWidth = SCREEN_WIDTH / 2;
					  //GLfloat halfScreenHeight = SCREEN_HEIGHT / 2;

	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LESS);
	glEnable(GL_LIGHTING); // ��������� ����.
	glEnable(GL_LIGHT0); // �������� ����� ������ 0.
	glEnable(GL_COLOR_MATERIAL); // ����������� ����� �������� ��������.

								 /* Loop until the user closes the window */
	while (!glfwWindowShouldClose(window))
	{

		if (binverse_color_black_2_white) {
			glClearColor(1.0f, 1.0f, 1.0f, 0.0f);
		}
		//glClearColor(0.7f, 1.0f, 0.7f, 0.0f);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		/* Render OpenGL here */

		glPushMatrix();

		float position[] = { 0,0,1,0 };
		glLightfv(GL_LIGHT0, GL_POSITION, position);

		glTranslatef(halfScreenWidth, halfScreenHeight, -500);
		glRotatef(rotationX, 1, 0, 0);
		glRotatef(rotationY, 0, 1, 0);
		glTranslatef(-halfScreenWidth, -halfScreenHeight, 500);


		DrawKarkas(b, lb, halfScreenWidth, halfScreenHeight, -500.0, scale_all);


		glPopMatrix();

		/* Swap front and back buffers */
		glfwSwapBuffers(window);

		/* Poll for and process events */
		glfwPollEvents();
	}

	glfwTerminate();

#endif

	for (integer i_7 = 0; i_7 < lb; i_7++) {
		if (b[i_7].temp_Sc != nullptr) {
			delete[] b[i_7].temp_Sc;
			b[i_7].temp_Sc = nullptr;
		}
		if (b[i_7].arr_Sc != nullptr) {
			delete[] b[i_7].arr_Sc;
			b[i_7].arr_Sc = nullptr;
		}
		if (b[i_7].g.hi != nullptr) {
			delete[] b[i_7].g.hi;
			b[i_7].g.hi = nullptr;
		}
		if (b[i_7].g.xi != nullptr) {
			delete[] b[i_7].g.xi;
			b[i_7].g.xi = nullptr;
		}
		if (b[i_7].g.yi != nullptr) {
			delete[] b[i_7].g.yi;
			b[i_7].g.yi = nullptr;
		}
		if (b[i_7].g.zi != nullptr) {
			delete[] b[i_7].g.zi;
			b[i_7].g.zi = nullptr;
		}
	}
	delete[] b;
	b = nullptr;

#ifndef NO_OPENGL_GLFW
	delete[] mask_color_for_openGL;
#endif

	return 0;
}



#ifndef NO_OPENGL_GLFW

// 16.12.2020
// ��������� �� ��� ����� ����� ������ ������ ��� ����������.
// ��������� ������������ ��������� � ������� ��������� Z ������.
bool CheckLine1(int il1, int il2)
{
	if ((il1 != 4) &&
		(il1 != 8) &&
		(il2 != 4) &&
		(il2 != 8) &&
		(il1 == il2))
	{
		return true;
	}
	else
	{
		if (((il1 == 1) && ((il2 == 2) || (il2 == 3) || (il2 == 5)))
			|| ((il2 == 1) && ((il1 == 2) || (il1 == 3) || (il1 == 5))))
		{
			return true;
		}
		else if ((il1 == 2) && (il2 == 2)) {
			return true;
		}
		else if (((il1 == 3) && (il2 == 2)) ||
			((il2 == 3) && (il1 == 2))) {

			return true;
		}
		else if (((il1 == 7) && (il2 == 5)) ||
			((il2 == 7) && (il1 == 5)))
		{
			return true;
		}
		else if (((il1 == 7) && (il2 == 3)) ||
			((il2 == 7) && (il1 == 3)))
		{
			return  true;
		}
		//if (((il1==6)&&(il2==3))||
		//((il2==6)&&(il1==3))) 
		//{
		// return true;
		//}
		//else
		else if (((il1 == 6) && (il2 == 5)) ||
			((il2 == 6) && (il1 == 5)))
		{
			return  true;
		}
		else
		{
			if (((il1 == 6) && (il2 == 7)) ||
				((il2 == 6) && (il1 == 7)))
			{
				return true;
			}
			else
			{
				return  false;
			}
		}

	}
	return false;
} // CheckLine

  // ���������� �������� ���������� ����� � ������� ��� �������.
const int CHECK_LINE_SIZE = 9;
bool* checkline1 = nullptr, *checkline2 = nullptr;
bool CheckLine(int il1, int il2)
{
	int key = il1 * CHECK_LINE_SIZE + il2;
	// ��� ��������� ���������� ��� ������� ������������� ������� ���� ���.
	//if (checkline1[key]) {
	return checkline2[key];
	//}
	/*else {
	checkline1[key] = true;
	checkline2[key] = CheckLine1(il1, il2);
	return checkline2[key];
	}
	*/
}

//  �������������� �������������� � ����� �������� � ���������� 18.12.2020.
// ������ ������� ����������� ���������� ��� �������, � ������ �������� � ��������������.
bool binvisible_face_detect(doublereal nx, doublereal ny, doublereal nz)
{
	//if (0) {
	//return true;
	//}
	//else {


	//bool br;
	//GLfloat mx[3][3];
	//GLfloat my[3][3];
	//doublereal mz0[3][3];
	//doublereal mz[3][3];
	//GLfloat mr[3][3];
	//doublereal mr1[3][3];

	//int i88, j88, k88;
	GLfloat /*nx1, ny1,*/ nz1;

	// (nx,ny,nz) - �������.
	//
	//br = true;
	//

	//sinAlf, cosAlf, sinBet, cosBet - ����������� ���� ���.

	// glRotatef(180.0*Gam0/3.141,0.0,0.0,1.0); // z apriory
	//glRotatef(180.0*Alf/3.141,1.0,0.0,0.0); // x
	//glRotatef(180.0*Bet/3.141,0.0,1.0,0.0); // y
	//glRotatef(180.0*Gam/3.141,0.0,0.0,1.0); // z
	// ������� �������� ������ ��� Oz
	//mz0[0][0]: = cos(Gam0); mz0[0][1]: = -sin(Gam0); mz0[0][2]: = 0.0;
	//mz0[1][0]: = sin(Gam0); mz0[1][1]: = cos(Gam0); mz0[1][2]: = 0.0;
	//mz0[2][0]: = 0.0; mz0[2][1]: = 0.0; mz0[2][2]: = 1.0;
	/*
	mz0[0][0] = 1.0; mz0[0][1] = 0.0; mz0[0][2] = 0.0;
	mz0[1][0] = 0.0; mz0[1][1] = 1.0; mz0[1][2] = 0.0;
	mz0[2][0] = 0.0; mz0[2][1] = 0.0; mz0[2][2] = 1.0;
	*/
	// ������� �������� ������ ��� Ox
	//mx[0][0] = 1.0; mx[0][1] = 0.0; mx[0][2] = 0.0;
	//mx[1][0] = 0.0; mx[1][1] = cosAlf; mx[1][2] = -sinAlf;
	//mx[2][0] = 0.0; mx[2][1] = sinAlf; mx[2][2] = cosAlf;
	// ������� �������� ������ ��� Oy
	//my[0][0] = cosBet; my[0][1] = 0.0; my[0][2] = sinBet;
	//my[1][0] = 0.0; my[1][1] = 1.0; my[1][2] = 0.0;
	//my[2][0] = -sinBet; my[2][1] = 0.0; my[2][2] = cosBet;
	// ������� �������� ������ ��� Oz
	//mz[0][0] = cos(Gam); mz[0][1] = -sin(Gam); mz[0][2] = 0.0;
	//mz[1][0] = sin(Gam); mz[1][1] = cos(Gam); mz[1][2] = 0.0;
	//mz[2][0] = 0.0; mz[2][1] = 0.0; mz[2][2] = 1.0;

	/*
	mz[0][0] = 1.0; mz[0][1] = 0.0; mz[0][2] = 0.0;
	mz[1][0] = 0.0; mz[1][1] = 1.0; mz[1][2] = 0.0;
	mz[2][0] = 0.0; mz[2][1] = 0.0; mz[2][2] = 1.0;
	*/

	//mr[0][0] = 0.0; mr[0][1] = 0.0; mr[0][2] = 0.0;
	//mr[1][0] = 0.0; mr[1][1] = 0.0; mr[1][2] = 0.0;
	/*mr[2][0] = 0.0; mr[2][1] = 0.0; mr[2][2] = 0.0;


	//for (i88 = 0; i88 < 3; i88++)
	i88 = 2; // ����� ������ Z ����������.
	{
	for (j88 = 0; j88 < 3; j88++) {

	for (k88 = 0; k88 < 3; k88++) {

	mr[i88][j88] = mr[i88][j88] + mx[i88][k88] * my[k88][j88];
	}
	}
	}*/
	/*
	for (i88 = 0; i88 < 3; i88++) {

	for (j88 = 0; j88 < 3; j88++) {

	for (k88 = 0; k88 < 3; k88++) {

	mr[i88][j88] = mr[i88][j88] + mz0[i88][k88] * mx[k88][j88];
	}
	}
	}

	mr1[0][0] = 0.0; mr1[0][1] = 0.0; mr1[0][2] = 0.0;
	mr1[1][0] = 0.0; mr1[1][1] = 0.0; mr1[1][2] = 0.0;
	mr1[2][0] = 0.0; mr1[2][1] = 0.0; mr1[2][2] = 0.0;
	for (i88 = 0; i88 < 3; i88++) {

	for (j88 = 0; j88 < 3; j88++) {

	for (k88 = 0; k88 < 3; k88++) {

	mr1[i88][j88] = mr1[i88][j88] + mr[i88][k88] * my[k88][j88];
	}
	}
	}
	for (i88 = 0; i88 < 3; i88++) {

	for (j88 = 0; j88 < 3; j88++) {

	mr[i88][j88] = mr1[i88][j88];
	mr1[i88][j88] = 0.0;
	}
	}
	for (i88 = 0; i88 < 3; i88++) {

	for (j88 = 0; j88 < 3; j88++) {

	for (k88 = 0; k88 < 3; k88++) {

	mr1[i88][j88] = mr1[i88][j88] + mr[i88][k88] * mz[k88][j88];
	}
	}
	}
	nx1 = 0.0;
	ny1 = 0.0;
	nz1 = 0.0;
	nx1 = mr1[0][0] * nx + mr1[0][1] * ny + mr1[0][2] * nz;
	ny1 = mr1[1][0] * nx + mr1[1][1] * ny + mr1[1][2] * nz;
	nz1 = mr1[2][0] * nx + mr1[2][1] * ny + mr1[2][2] * nz;
	*/
	//nx1 = 0.0;
	//ny1 = 0.0;
	//nz1 = 0.0;
	//nx1 = mr[0][0] * nx + mr[0][1] * ny + mr[0][2] * nz;
	//ny1 = mr[1][0] * nx + mr[1][1] * ny + mr[1][2] * nz;
	//nz1 = mr[2][0] * nx + mr[2][1] * ny + mr[2][2] * nz;

	// �� ����������� �������, �� ������ ������� ��������.
	// ��� ���� ������ � �������� ��������� ������� � ���������.
	// ������ ������� ���������� ����������� ����� ��� � �������������� ����� �����.
	// 18.12.2020.
	nz1 = -sinBet * cosAlf * nx + sinAlf * ny + cosAlf * cosBet * nz;

	//if (nx*nx1+ny*ny1+nz*nz1>-1.0e-3)
	if (nz1 > -1.0e-3)
	{
		return true;
	}
	else
	{
		return false;
	}


	//}
}

//  �������� ��������.
bool binvisible_face_detectold(doublereal nx, doublereal ny, doublereal nz)
{
	if (0) {
		return true;
	}
	else {


		bool br;
		doublereal mx[3][3];
		doublereal my[3][3];
		doublereal mz0[3][3];
		doublereal mz[3][3];
		doublereal mr[3][3];
		doublereal mr1[3][3];

		int i88, j88, k88;
		doublereal nx1, ny1, nz1;

		// (nx,ny,nz) - �������.
		//
		br = true;
		//

		//sinAlf, cosAlf, sinBet, cosBet - ����������� ���� ���.

		// glRotatef(180.0*Gam0/3.141,0.0,0.0,1.0); // z apriory
		//glRotatef(180.0*Alf/3.141,1.0,0.0,0.0); // x
		//glRotatef(180.0*Bet/3.141,0.0,1.0,0.0); // y
		//glRotatef(180.0*Gam/3.141,0.0,0.0,1.0); // z
		// ������� �������� ������ ��� Oz
		//mz0[0][0]: = cos(Gam0); mz0[0][1]: = -sin(Gam0); mz0[0][2]: = 0.0;
		//mz0[1][0]: = sin(Gam0); mz0[1][1]: = cos(Gam0); mz0[1][2]: = 0.0;
		//mz0[2][0]: = 0.0; mz0[2][1]: = 0.0; mz0[2][2]: = 1.0;

		mz0[0][0] = 1.0; mz0[0][1] = 0.0; mz0[0][2] = 0.0;
		mz0[1][0] = 0.0; mz0[1][1] = 1.0; mz0[1][2] = 0.0;
		mz0[2][0] = 0.0; mz0[2][1] = 0.0; mz0[2][2] = 1.0;
		// ������� �������� ������ ��� Ox
		mx[0][0] = 1.0; mx[0][1] = 0.0; mx[0][2] = 0.0;
		mx[1][0] = 0.0; mx[1][1] = cosAlf; mx[1][2] = -sinAlf;
		mx[2][0] = 0.0; mx[2][1] = sinAlf; mx[2][2] = cosAlf;
		// ������� �������� ������ ��� Oy
		my[0][0] = cosBet; my[0][1] = 0.0; my[0][2] = sinBet;
		my[1][0] = 0.0; my[1][1] = 1.0; my[1][2] = 0.0;
		my[2][0] = -sinBet; my[2][1] = 0.0; my[2][2] = cosBet;
		// ������� �������� ������ ��� Oz
		//mz[0][0] = cos(Gam); mz[0][1] = -sin(Gam); mz[0][2] = 0.0;
		//mz[1][0] = sin(Gam); mz[1][1] = cos(Gam); mz[1][2] = 0.0;
		//mz[2][0] = 0.0; mz[2][1] = 0.0; mz[2][2] = 1.0;

		mz[0][0] = 1.0; mz[0][1] = 0.0; mz[0][2] = 0.0;
		mz[1][0] = 0.0; mz[1][1] = 1.0; mz[1][2] = 0.0;
		mz[2][0] = 0.0; mz[2][1] = 0.0; mz[2][2] = 1.0;

		mr[0][0] = 0.0; mr[0][1] = 0.0; mr[0][2] = 0.0;
		mr[1][0] = 0.0; mr[1][1] = 0.0; mr[1][2] = 0.0;
		mr[2][0] = 0.0; mr[2][1] = 0.0; mr[2][2] = 0.0;
		for (i88 = 0; i88 < 3; i88++) {

			for (j88 = 0; j88 < 3; j88++) {

				for (k88 = 0; k88 < 3; k88++) {

					mr[i88][j88] = mr[i88][j88] + mz0[i88][k88] * mx[k88][j88];
				}
			}
		}

		mr1[0][0] = 0.0; mr1[0][1] = 0.0; mr1[0][2] = 0.0;
		mr1[1][0] = 0.0; mr1[1][1] = 0.0; mr1[1][2] = 0.0;
		mr1[2][0] = 0.0; mr1[2][1] = 0.0; mr1[2][2] = 0.0;
		for (i88 = 0; i88 < 3; i88++) {

			for (j88 = 0; j88 < 3; j88++) {

				for (k88 = 0; k88 < 3; k88++) {

					mr1[i88][j88] = mr1[i88][j88] + mr[i88][k88] * my[k88][j88];
				}
			}
		}
		for (i88 = 0; i88 < 3; i88++) {

			for (j88 = 0; j88 < 3; j88++) {

				mr[i88][j88] = mr1[i88][j88];
				mr1[i88][j88] = 0.0;
			}
		}
		for (i88 = 0; i88 < 3; i88++) {

			for (j88 = 0; j88 < 3; j88++) {

				for (k88 = 0; k88 < 3; k88++) {

					mr1[i88][j88] = mr1[i88][j88] + mr[i88][k88] * mz[k88][j88];
				}
			}
		}
		nx1 = 0.0;
		ny1 = 0.0;
		nz1 = 0.0;
		nx1 = mr1[0][0] * nx + mr1[0][1] * ny + mr1[0][2] * nz;
		ny1 = mr1[1][0] * nx + mr1[1][1] * ny + mr1[1][2] * nz;
		nz1 = mr1[2][0] * nx + mr1[2][1] * ny + mr1[2][2] * nz;

		//if (nx*nx1+ny*ny1+nz*nz1>-1.0e-3)
		if (nz1 > -1.0e-3)
		{
			br = true;
		}
		else
		{
			br = false;
		}

		if (br)
		{
			return true;
		}
		else
		{
			return false;
		}
	}
}

void Drawnvtx(int**& nvtx, TOCHKA*& pa, int id, GLfloat centerPosX, GLfloat centerPosY, GLfloat  centerPosZ, GLfloat scale);


// ��������� ������������ � ��������� ��������� �����. ��������� 19.12.2020.
// ������ ������ �������� ��� ������������ ������, �� �������� ���������� � ������.
int DrawZbufferColor(int maxelm, int maxnod, int**& nvtx) {

	int CHECK_LINE_SIZE = 9;
	checkline1 = new bool[CHECK_LINE_SIZE * CHECK_LINE_SIZE];
	checkline2 = new bool[CHECK_LINE_SIZE * CHECK_LINE_SIZE];
	for (int i = 0; i < CHECK_LINE_SIZE * CHECK_LINE_SIZE; i++) {
		checkline1[i] = false;
		checkline2[i] = false;
	}

	for (int i = 0; i < CHECK_LINE_SIZE; i++) {
		for (int j = 0; j < CHECK_LINE_SIZE; j++) {
			int key = i * CHECK_LINE_SIZE + j;

			checkline1[key] = true;
			checkline2[key] = CheckLine1(i, j);
		}
	}

	bool** bvisible_gran = new bool*[7];
	for (int i = 0; i < 7; i++) {
		bvisible_gran[i] = new bool[maxelm];
	}

	int* ipa_count = new int[maxnod];
	for (int i = 0; i < maxnod; i++) {
		ipa_count[i] = 0;
	}
	for (int i = 0; i < maxelm; i++) {

		bvisible_gran[E_SIDE][i] = false;
		bvisible_gran[W_SIDE][i] = false;
		bvisible_gran[N_SIDE][i] = false;
		bvisible_gran[S_SIDE][i] = false;
		bvisible_gran[T_SIDE][i] = false;
		bvisible_gran[B_SIDE][i] = false;
		bvisible_gran[6][i] = false; // ��� ������ �������!!!


		integer inode1, inode2, inode3, inode4, inode5, inode6, inode7, inode8;
		inode1 = nvtx[0][i] - 1;
		inode2 = nvtx[1][i] - 1;
		inode3 = nvtx[2][i] - 1;
		inode4 = nvtx[3][i] - 1;
		inode5 = nvtx[4][i] - 1;
		inode6 = nvtx[5][i] - 1;
		inode7 = nvtx[6][i] - 1;
		inode8 = nvtx[7][i] - 1;

		TOCHKA ptest;
		//center_cord3D(inode1, nvtx, pa, ptest, 100);
		ptest.x = pa_opengl[inode1].x;
		ptest.y = pa_opengl[inode1].y;
		ptest.z = pa_opengl[inode1].z;
		bool bvisible = false;

		if ((fabs(pfpir.fminimum) < 1.0e-30) && (fabs(pfpir.fmaximum) < 1.0e-30)) {
			pfpir.fminimum = -1.0e30;
			pfpir.fmaximum = 1.0e30;
		}

		switch (pfpir.idir) {
		case LINE_DIRECTIONAL::X_LINE_DIRECTIONAL:
			if ((ptest.x >= pfpir.fminimum) && (ptest.x <= pfpir.fmaximum)) {
				bvisible = true;
			}
			break;
		case LINE_DIRECTIONAL::Y_LINE_DIRECTIONAL:
			if ((ptest.y >= pfpir.fminimum) && (ptest.y <= pfpir.fmaximum)) {
				bvisible = true;
			}
			break;
		case LINE_DIRECTIONAL::Z_LINE_DIRECTIONAL:
			if ((ptest.z >= pfpir.fminimum) && (ptest.z <= pfpir.fmaximum)) {
				bvisible = true;
			}
			break;
		}

		if (bvisible) {

			for (int j = 0; j < NUMBER_OF_VERTEX_FINITE_ELEMENT(); j++) {
				//std::cout << nvtx[j][i] - 1 << " ";
				ipa_count[nvtx[j][i] - 1]++;
			}

			//getchar();
		}
	}

	const int bVisibleCount = 7;

	for (int j = 0; j < maxelm; j++) {

		integer inode1, inode2, inode3, inode4, inode5, inode6, inode7, inode8;
		inode1 = nvtx[0][j] - 1;
		inode2 = nvtx[1][j] - 1;
		inode3 = nvtx[2][j] - 1;
		inode4 = nvtx[3][j] - 1;
		inode5 = nvtx[4][j] - 1;
		inode6 = nvtx[5][j] - 1;
		inode7 = nvtx[6][j] - 1;
		inode8 = nvtx[7][j] - 1;


		TOCHKA ptest;
		//center_cord3D(inode1, nvtx, pa, ptest, 100);
		ptest.x = pa_opengl[inode1].x;
		ptest.y = pa_opengl[inode1].y;
		ptest.z = pa_opengl[inode1].z;
		bool bvisible = false;

		switch (pfpir.idir) {
		case LINE_DIRECTIONAL::X_LINE_DIRECTIONAL:
			if ((ptest.x >= pfpir.fminimum) && (ptest.x <= pfpir.fmaximum)) {
				bvisible = true;
			}
			break;
		case LINE_DIRECTIONAL::Y_LINE_DIRECTIONAL:
			if ((ptest.y >= pfpir.fminimum) && (ptest.y <= pfpir.fmaximum)) {
				bvisible = true;
			}
			break;
		case LINE_DIRECTIONAL::Z_LINE_DIRECTIONAL:
			if ((ptest.z >= pfpir.fminimum) && (ptest.z <= pfpir.fmaximum)) {
				bvisible = true;
			}
			break;
		}

		if (bvisible) {


			bvisible_gran[6][j] = true;

			// XY bottom
			if ((ipa_count[inode1] <= bVisibleCount) &&
				(ipa_count[inode2] <= bVisibleCount) &&
				(ipa_count[inode3] <= bVisibleCount) &&
				(ipa_count[inode4] <= bVisibleCount)) {

				bvisible_gran[B_SIDE][j] = true;
			}

			// XY Top
			if ((ipa_count[inode5] <= bVisibleCount) &&
				(ipa_count[inode6] <= bVisibleCount) &&
				(ipa_count[inode7] <= bVisibleCount) &&
				(ipa_count[inode8] <= bVisibleCount)) {

				bvisible_gran[T_SIDE][j] = true;
			}

			// XZ SSIDE min Y
			if ((ipa_count[inode1] <= bVisibleCount) &&
				(ipa_count[inode2] <= bVisibleCount) &&
				(ipa_count[inode6] <= bVisibleCount) &&
				(ipa_count[inode5] <= bVisibleCount)) {

				bvisible_gran[S_SIDE][j] = true;
			}

			// XZ SSIDE max Y
			if ((ipa_count[inode4] <= bVisibleCount) &&
				(ipa_count[inode8] <= bVisibleCount) &&
				(ipa_count[inode7] <= bVisibleCount) &&
				(ipa_count[inode3] <= bVisibleCount)) {

				bvisible_gran[N_SIDE][j] = true;
			}

			// YZ SSIDE max X
			if ((ipa_count[inode2] <= bVisibleCount) &&
				(ipa_count[inode3] <= bVisibleCount) &&
				(ipa_count[inode7] <= bVisibleCount) &&
				(ipa_count[inode6] <= bVisibleCount)) {

				bvisible_gran[E_SIDE][j] = true;
			}

			// YZ SSIDE min X
			if ((ipa_count[inode1] <= bVisibleCount) &&
				(ipa_count[inode5] <= bVisibleCount) &&
				(ipa_count[inode8] <= bVisibleCount) &&
				(ipa_count[inode4] <= bVisibleCount))
			{
				bvisible_gran[W_SIDE][j] = true;
			}


		}
	}

	GLFWwindow* window;

	/* Initialize the library */
	if (!glfwInit())
		return -1;

	/* Create a windowed mode window and its OpenGL context */
	window = glfwCreateWindow(SCREEN_WIDTH, SCREEN_HEIGHT, "AliceFlow_v0_49", NULL, NULL);


	glfwSetKeyCallback(window, keyCallback);
	glfwSetInputMode(window, GLFW_STICKY_KEYS, 1);

	int screenWidth, screenHeight;
	glfwGetFramebufferSize(window, &screenWidth, &screenHeight);

	if (!window)
	{
		glfwTerminate();
		return -1;
	}

	/* Make the window's context current */
	glfwMakeContextCurrent(window);


	glViewport(0.0f, 0.0f, screenWidth, screenHeight); // specifies the part of the window to which OpenGL will draw (in pixels), convert from
													   // normalised to pixels
	glMatrixMode(GL_PROJECTION); // projection matrix defines the properties of the camera that views the objects in the world coordinate frame.
								 // Here you thalfScreenHeight+scale_all*pa.yally set the zoom factor, aspect ratio and the near and far clipping planes
	glLoadIdentity(); // replace the current matrix with the identity matrix and starts us a fresh because matrix transforms such as glOrtho and
					  // glRotate cumulate, besically puts us at (0, 0, 0)
	glOrtho(0, SCREEN_WIDTH, 0, SCREEN_HEIGHT, 0.1, 10000); // essentially set coordinate system
	glMatrixMode(GL_MODELVIEW); // (default matrix mode) modelview matrix defines how objects are transformed (meaning translation, rotation
								// and scaling) in your world

	glLoadIdentity(); // same as above comment 



	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LESS);

	/* Loop until the user closes the window */
	while (!glfwWindowShouldClose(window))
	{
		if (binverse_color_black_2_white) {
			glClearColor(1.0f, 1.0f, 1.0f, 0.0f);
		}
		//glColor3f(0.0, 0.0, 0.0); // ׸����
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		/* Render OpenGL here */

		glPushMatrix();
		glTranslatef(halfScreenWidth, halfScreenHeight, -abbys);
		glRotatef(rotationX, 1, 0, 0);
		glRotatef(rotationY, 0, 1, 0);
		glTranslatef(-halfScreenWidth, -halfScreenHeight, abbys);

		glLineWidth(1);

		// �������� ���.
		for (int j = 0; j < maxelm; j++) {

			if (bvisible_gran[6][j]) {

				int inode1, inode2, inode3, inode4, inode5, inode6, inode7, inode8;
				inode1 = nvtx[0][j] - 1;
				inode2 = nvtx[1][j] - 1;
				inode3 = nvtx[2][j] - 1;
				inode4 = nvtx[3][j] - 1;
				inode5 = nvtx[4][j] - 1;
				inode6 = nvtx[5][j] - 1;
				inode7 = nvtx[6][j] - 1;
				inode8 = nvtx[7][j] - 1;


				// XY bottom
				if (bvisible_gran[B_SIDE][j])
				{

					if (binvisible_face_detect(0.0, 0.0, -1.0)) {


						//glColor3f(0.84,0.84,0.84);

						//glColor3f(0.0, 0.0, 0.0); // ׸����


						//Drawnvtx(nvtx, pa, j, halfScreenWidth, halfScreenHeight, -abbys, scale_all);


						glBegin(GL_QUADS);

						set_color_for_render(inode4);
						//glColor3f(1.0, 1.0, 1.0);
						//glNormal3f(0.0, 0.0, -1.0);
						glVertex3f(pa_render[inode4].x, pa_render[inode4].y, pa_render[inode4].z);
						set_color_for_render(inode3);
						//glColor3f(1.0, 1.0, 1.0);
						//glNormal3f(0.0, 0.0, -1.0);
						glVertex3f(pa_render[inode3].x, pa_render[inode3].y, pa_render[inode3].z);
						set_color_for_render(inode2);
						//glColor3f(1.0, 1.0, 1.0);
						//glNormal3f(0.0, 0.0, -1.0);
						glVertex3f(pa_render[inode2].x, pa_render[inode2].y, pa_render[inode2].z);
						set_color_for_render(inode1);
						//glColor3f(1.0, 1.0, 1.0);
						//glNormal3f(0.0, 0.0, -1.0);
						glVertex3f(pa_render[inode1].x, pa_render[inode1].y, pa_render[inode1].z);
						//glColor3f(0.0, 0.0, 0.0);
						glEnd();


						if (binverse_color_black_2_white) {
							glColor3f(0.0, 0.0, 0.0); // ׸����
						}
						else {
							glColor3f(1.0, 1.0, 1.0);// �����
						}

						if (CheckLine(ipa_count[inode1],
							ipa_count[inode2])
							)
						{

							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode1].x, pa_render[inode1].y, pa_render[inode1].z);
							glVertex3f(pa_render[inode2].x, pa_render[inode2].y, pa_render[inode2].z);
							glEnd();
						}

						if (CheckLine(ipa_count[inode3],
							ipa_count[inode2])
							)
						{

							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode2].x, pa_render[inode2].y, pa_render[inode2].z);
							glVertex3f(pa_render[inode3].x, pa_render[inode3].y, pa_render[inode3].z);
							glEnd();
						}

						if (CheckLine(ipa_count[inode3],
							ipa_count[inode4])
							)
						{
							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode3].x, pa_render[inode3].y, pa_render[inode3].z);
							glVertex3f(pa_render[inode4].x, pa_render[inode4].y, pa_render[inode4].z);
							glEnd();
						}

						if (CheckLine(ipa_count[inode1],
							ipa_count[inode4])
							)
						{
							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode4].x, pa_render[inode4].y, pa_render[inode4].z);
							glVertex3f(pa_render[inode1].x, pa_render[inode1].y, pa_render[inode1].z);
							glEnd();
						}

						//glLineWidth(1);
						//glColor3f(1.0, 1.0, 1.0);
						//glLineWidth(1);

					}
				}


				// XY Top
				if (bvisible_gran[T_SIDE][j])
				{


					if (binvisible_face_detect(0.0, 0.0, 1.0)) {


						//glColor3f(0.84,0.84,0.84);
						//glColor3f(0.0, 0.0, 0.0); // ׸����

						//Drawnvtx(nvtx, pa, j, halfScreenWidth, halfScreenHeight, -abbys, scale_all);


						glBegin(GL_QUADS);
						set_color_for_render(inode8);
						//glColor3f(1.0, 1.0, 1.0);
						//glNormal3f(0.0, 0.0, 1.0);
						glVertex3f(pa_render[inode8].x, pa_render[inode8].y, pa_render[inode8].z);
						set_color_for_render(inode7);
						//glColor3f(1.0, 1.0, 1.0);
						//glNormal3f(0.0, 0.0, 1.0);						
						glVertex3f(pa_render[inode7].x, pa_render[inode7].y, pa_render[inode7].z);
						set_color_for_render(inode6);
						//glColor3f(1.0, 1.0, 1.0);
						//glNormal3f(0.0, 0.0, 1.0);						
						glVertex3f(pa_render[inode6].x, pa_render[inode6].y, pa_render[inode6].z);
						set_color_for_render(inode5);
						//glColor3f(1.0, 1.0, 1.0);
						//glNormal3f(0.0, 0.0, 1.0);
						glVertex3f(pa_render[inode5].x, pa_render[inode5].y, pa_render[inode5].z);
						//glColor3f(0.0, 0.0, 0.0);
						glEnd();


						//glLineWidth(ComboBoxlineWidth.ItemIndex + 1);
						if (binverse_color_black_2_white) {
							glColor3f(0.0, 0.0, 0.0); // ׸����
						}
						else {
							glColor3f(1.0, 1.0, 1.0);// �����
						}

						if (CheckLine(ipa_count[inode5],
							ipa_count[inode6])
							)
						{

							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode5].x, pa_render[inode5].y, pa_render[inode5].z);
							glVertex3f(pa_render[inode6].x, pa_render[inode6].y, pa_render[inode6].z);
							glEnd();
						}

						if (CheckLine(ipa_count[inode6],
							ipa_count[inode7])
							) {

							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode6].x, pa_render[inode6].y, pa_render[inode6].z);
							glVertex3f(pa_render[inode7].x, pa_render[inode7].y, pa_render[inode7].z);
							glEnd();
						}

						if (CheckLine(ipa_count[inode7],
							ipa_count[inode8])
							)
						{
							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode7].x, pa_render[inode7].y, pa_render[inode7].z);
							glVertex3f(pa_render[inode8].x, pa_render[inode8].y, pa_render[inode8].z);
							glEnd();
						}

						if (CheckLine(ipa_count[inode5],
							ipa_count[inode8])
							)
						{
							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode8].x, pa_render[inode8].y, pa_render[inode8].z);
							glVertex3f(pa_render[inode5].x, pa_render[inode5].y, pa_render[inode5].z);
							glEnd();
						}

						// glLineWidth(1);
						//glColor3f(1.0, 1.0, 1.0);
						//glLineWidth(1);

					}
				}



				// XZ SSIDE min Y
				if (bvisible_gran[S_SIDE][j])
				{


					if (binvisible_face_detect(0.0, -1.0, 0.0)) {

						//glColor3f(0.84,0.84,0.84);
						//glColor3f(0.0, 0.0, 0.0); // ׸����

						//Drawnvtx(nvtx, pa, j, halfScreenWidth, halfScreenHeight, -abbys, scale_all);


						glBegin(GL_QUADS);
						set_color_for_render(inode5);
						//glColor3f(1.0, 1.0, 1.0);
						//glNormal3f(0.0, -1.0, 0.0);
						glVertex3f(pa_render[inode5].x, pa_render[inode5].y, pa_render[inode5].z);
						set_color_for_render(inode1);
						//glColor3f(1.0, 1.0, 1.0);
						//glNormal3f(0.0, -1.0, 0.0);
						glVertex3f(pa_render[inode1].x, pa_render[inode1].y, pa_render[inode1].z);
						set_color_for_render(inode2);
						//glColor3f(1.0, 1.0, 1.0);
						//glNormal3f(0.0, -1.0, 0.0);
						glVertex3f(pa_render[inode2].x, pa_render[inode2].y, pa_render[inode2].z);
						set_color_for_render(inode6);
						//glColor3f(1.0, 1.0, 1.0);
						///glNormal3f(0.0, -1.0, 0.0);
						glVertex3f(pa_render[inode6].x, pa_render[inode6].y, pa_render[inode6].z);
						//glColor3f(0.0, 0.0, 0.0);
						glEnd();


						//glLineWidth(ComboBoxlineWidth.ItemIndex + 1);
						if (binverse_color_black_2_white) {
							glColor3f(0.0, 0.0, 0.0); // ׸����
						}
						else {
							glColor3f(1.0, 1.0, 1.0);// �����
						}

						if (CheckLine(ipa_count[inode1],
							ipa_count[inode2])
							) {

							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode1].x, pa_render[inode1].y, pa_render[inode1].z);
							glVertex3f(pa_render[inode2].x, pa_render[inode2].y, pa_render[inode2].z);
							glEnd();
						}

						if (CheckLine(ipa_count[inode2],
							ipa_count[inode6])
							) {

							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode2].x, pa_render[inode2].y, pa_render[inode2].z);
							glVertex3f(pa_render[inode6].x, pa_render[inode6].y, pa_render[inode6].z);
							glEnd();
						}

						if (CheckLine(ipa_count[inode6],
							ipa_count[inode5])
							) {

							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode6].x, pa_render[inode6].y, pa_render[inode6].z);
							glVertex3f(pa_render[inode5].x, pa_render[inode5].y, pa_render[inode5].z);
							glEnd();
						}

						if (CheckLine(ipa_count[inode1],
							ipa_count[inode5])
							) {

							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode5].x, pa_render[inode5].y, pa_render[inode5].z);
							glVertex3f(pa_render[inode1].x, pa_render[inode1].y, pa_render[inode1].z);
							glEnd();
						}

						//glLineWidth(1);
						//glColor3f(1.0, 1.0, 1.0);
						//glLineWidth(1);
					}
				}


				// XZ SSIDE max Y
				if (bvisible_gran[N_SIDE][j])
				{

					if (binvisible_face_detect(0.0, 1.0, 0.0)) {

						//glColor3f(0.84,0.84,0.84);
						//glColor3f(0.0, 0.0, 0.0); // ׸����

						//Drawnvtx(nvtx, pa, j, halfScreenWidth, halfScreenHeight, -abbys, scale_all);


						glBegin(GL_QUADS);
						set_color_for_render(inode8);
						//glColor3f(1.0, 1.0, 1.0);
						//glNormal3f(0.0, 1.0, 0.0);
						glVertex3f(pa_render[inode8].x, pa_render[inode8].y, pa_render[inode8].z);
						set_color_for_render(inode4);
						//glColor3f(1.0, 1.0, 1.0);
						//glNormal3f(0.0, 1.0, 0.0);
						glVertex3f(pa_render[inode4].x, pa_render[inode4].y, pa_render[inode4].z);
						set_color_for_render(inode3);
						//glColor3f(1.0, 1.0, 1.0);
						//glNormal3f(0.0, 1.0, 0.0);
						glVertex3f(pa_render[inode3].x, pa_render[inode3].y, pa_render[inode3].z);
						set_color_for_render(inode7);
						//glColor3f(1.0, 1.0, 1.0);
						//glNormal3f(0.0, 1.0, 0.0);
						glVertex3f(pa_render[inode7].x, pa_render[inode7].y, pa_render[inode7].z);

						//glColor3f(0.0, 0.0, 0.0);
						glEnd();


						//glLineWidth(ComboBoxlineWidth.ItemIndex + 1);
						if (binverse_color_black_2_white) {
							glColor3f(0.0, 0.0, 0.0); // ׸����
						}
						else {
							glColor3f(1.0, 1.0, 1.0);// �����
						}

						if (CheckLine(ipa_count[inode3],
							ipa_count[inode4])
							) {

							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode4].x, pa_render[inode4].y, pa_render[inode4].z);
							glVertex3f(pa_render[inode3].x, pa_render[inode3].y, pa_render[inode3].z);
							glEnd();
						}

						if (CheckLine(ipa_count[inode3],
							ipa_count[inode7])
							) {

							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode3].x, pa_render[inode3].y, pa_render[inode3].z);
							glVertex3f(pa_render[inode7].x, pa_render[inode7].y, pa_render[inode7].z);
							glEnd();
						}

						if (CheckLine(ipa_count[inode7],
							ipa_count[inode8])
							) {
							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode7].x, pa_render[inode7].y, pa_render[inode7].z);
							glVertex3f(pa_render[inode8].x, pa_render[inode8].y, pa_render[inode8].z);
							glEnd();
						}

						if (CheckLine(ipa_count[inode8],
							ipa_count[inode4])
							) {
							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode8].x, pa_render[inode8].y, pa_render[inode8].z);
							glVertex3f(pa_render[inode4].x, pa_render[inode4].y, pa_render[inode4].z);
							glEnd();
						}

						// glLineWidth(1);
						//glColor3f(1.0, 1.0, 1.0);
						//glLineWidth(1);
					}
				}



				// YZ SSIDE max X
				if (bvisible_gran[E_SIDE][j])
				{

					if (binvisible_face_detect(1.0, 0.0, 0.0)) {


						//glColor3f(0.84,0.84,0.84);
						//glColor3f(0.0, 0.0, 0.0); // ׸����

						//Drawnvtx(nvtx, pa, j, halfScreenWidth, halfScreenHeight, -abbys, scale_all);


						glBegin(GL_QUADS);
						set_color_for_render(inode6);
						//glColor3f(1.0, 1.0, 1.0);
						//glNormal3f(1.0, 0.0, 0.0);
						glVertex3f(pa_render[inode6].x, pa_render[inode6].y, pa_render[inode6].z);
						set_color_for_render(inode2);
						//glColor3f(1.0, 1.0, 1.0);
						//glNormal3f(1.0, 0.0, 0.0);
						glVertex3f(pa_render[inode2].x, pa_render[inode2].y, pa_render[inode2].z);
						set_color_for_render(inode3);
						//glColor3f(1.0, 1.0, 1.0);
						//glNormal3f(1.0, 0.0, 0.0);
						glVertex3f(pa_render[inode3].x, pa_render[inode3].y, pa_render[inode3].z);
						set_color_for_render(inode7);
						//glColor3f(1.0, 1.0, 1.0);
						//glNormal3f(1.0, 0.0, 0.0);
						glVertex3f(pa_render[inode7].x, pa_render[inode7].y, pa_render[inode7].z);
						//glColor3f(0.0, 0.0, 0.0);
						glEnd();


						//glLineWidth(ComboBoxlineWidth.ItemIndex + 1);
						if (binverse_color_black_2_white) {
							glColor3f(0.0, 0.0, 0.0); // ׸����
						}
						else {
							glColor3f(1.0, 1.0, 1.0);// �����
						}


						if (CheckLine(ipa_count[inode3],
							ipa_count[inode2])
							)
						{
							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode2].x, pa_render[inode2].y, pa_render[inode2].z);
							glVertex3f(pa_render[inode3].x, pa_render[inode3].y, pa_render[inode3].z);
							glEnd();
						}

						if (CheckLine(ipa_count[inode3],
							ipa_count[inode7])
							)
						{
							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode3].x, pa_render[inode3].y, pa_render[inode3].z);
							glVertex3f(pa_render[inode7].x, pa_render[inode7].y, pa_render[inode7].z);
							glEnd();
						}

						if (CheckLine(ipa_count[inode7],
							ipa_count[inode6])
							)
						{
							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode7].x, pa_render[inode7].y, pa_render[inode7].z);
							glVertex3f(pa_render[inode6].x, pa_render[inode6].y, pa_render[inode6].z);
							glEnd();
						}

						if (CheckLine(ipa_count[inode6],
							ipa_count[inode2])
							)
						{
							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode6].x, pa_render[inode6].y, pa_render[inode6].z);
							glVertex3f(pa_render[inode2].x, pa_render[inode2].y, pa_render[inode2].z);
							glEnd();
						}

						//glLineWidth(1);
						//glColor3f(1.0, 1.0, 1.0);
						//glLineWidth(1);

					}
				}


				// YZ SSIDE min X
				if (bvisible_gran[W_SIDE][j])
				{

					if (binvisible_face_detect(-1.0, 0.0, 0.0))
					{

						//glColor3f(0.84,0.84,0.84);

						//glColor3f(0.0, 0.0, 0.0); // ׸����

						//Drawnvtx(nvtx, pa, j, halfScreenWidth, halfScreenHeight, -abbys, scale_all);


						glBegin(GL_QUADS);
						set_color_for_render(inode8);
						//glColor3f(1.0, 1.0, 1.0);
						//glNormal3f(-1.0, 0.0, 0.0);
						glVertex3f(pa_render[inode8].x, pa_render[inode8].y, pa_render[inode8].z);
						set_color_for_render(inode4);
						//glColor3f(1.0, 1.0, 1.0);
						//glNormal3f(-1.0, 0.0, 0.0);
						glVertex3f(pa_render[inode4].x, pa_render[inode4].y, pa_render[inode4].z);
						set_color_for_render(inode1);
						//glColor3f(1.0, 1.0, 1.0);
						//glNormal3f(-1.0, 0.0, 0.0);
						glVertex3f(pa_render[inode1].x, pa_render[inode1].y, pa_render[inode1].z);
						set_color_for_render(inode5);
						//glColor3f(1.0, 1.0, 1.0);
						//glNormal3f(-1.0, 0.0, 0.0);
						glVertex3f(pa_render[inode5].x, pa_render[inode5].y, pa_render[inode5].z);
						//glColor3f(0.0, 0.0, 0.0);
						glEnd();



						//glLineWidth(ComboBoxlineWidth.ItemIndex + 1);

						if (binverse_color_black_2_white) {
							glColor3f(0.0, 0.0, 0.0); // ׸����
						}
						else {
							glColor3f(1.0, 1.0, 1.0);// �����
						}

						if (CheckLine(ipa_count[inode1],
							ipa_count[inode4])
							) {

							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode1].x, pa_render[inode1].y, pa_render[inode1].z);
							glVertex3f(pa_render[inode4].x, pa_render[inode4].y, pa_render[inode4].z);
							glEnd();
						}

						if (CheckLine(ipa_count[inode4],
							ipa_count[inode8])
							)
						{
							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode4].x, pa_render[inode4].y, pa_render[inode4].z);
							glVertex3f(pa_render[inode8].x, pa_render[inode8].y, pa_render[inode8].z);
							glEnd();
						}

						if (CheckLine(ipa_count[inode5],
							ipa_count[inode8])
							)
						{
							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode8].x, pa_render[inode8].y, pa_render[inode8].z);
							glVertex3f(pa_render[inode5].x, pa_render[inode5].y, pa_render[inode5].z);
							glEnd();
						}

						if (CheckLine(ipa_count[inode1],
							ipa_count[inode5])
							)
						{
							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode5].x, pa_render[inode5].y, pa_render[inode5].z);
							glVertex3f(pa_render[inode1].x, pa_render[inode1].y, pa_render[inode1].z);
							glEnd();
						}

						//glLineWidth(1);
						//glColor3f(1.0, 1.0, 1.0);
						//glLineWidth(1);
					}
				}

			}
		}


		glPopMatrix();

		/* Swap front and back buffers */
		glfwSwapBuffers(window);

		/* Poll for and process events */
		glfwPollEvents();
	}

	delete[] ipa_count;

	for (int i = 0; i < 7; i++) {
		delete[] bvisible_gran[i];
	}
	delete[] bvisible_gran;

	delete[] checkline1;
	delete[] checkline2;

	glfwTerminate();

	return 0;

} // DrawZbufferColor


  // ��������� ������������ � ��������� ��������� ����� � � ����������. 
int DrawZbufferColorLight(BLOCK*& b, integer lb, int maxelm, int maxnod, int**& nvtx, int*& whot_is_block, doublereal*& potent) {

	int CHECK_LINE_SIZE = 9;
	checkline1 = new bool[CHECK_LINE_SIZE * CHECK_LINE_SIZE];
	checkline2 = new bool[CHECK_LINE_SIZE * CHECK_LINE_SIZE];
	for (int i = 0; i < CHECK_LINE_SIZE * CHECK_LINE_SIZE; i++) {
		checkline1[i] = false;
		checkline2[i] = false;
	}

	for (int i = 0; i < CHECK_LINE_SIZE; i++) {
		for (int j = 0; j < CHECK_LINE_SIZE; j++) {
			int key = i * CHECK_LINE_SIZE + j;

			checkline1[key] = true;
			checkline2[key] = CheckLine1(i, j);
		}
	}

	bool** bvisible_gran = new bool*[7];
	for (int i = 0; i < 7; i++) {
		bvisible_gran[i] = new bool[maxelm];
	}

	int* ipa_count = new int[maxnod];
	for (int i = 0; i < maxnod; i++) {
		ipa_count[i] = 0;
	}
	for (int i = 0; i < maxelm; i++) {

		bvisible_gran[E_SIDE][i] = false;
		bvisible_gran[W_SIDE][i] = false;
		bvisible_gran[N_SIDE][i] = false;
		bvisible_gran[S_SIDE][i] = false;
		bvisible_gran[T_SIDE][i] = false;
		bvisible_gran[B_SIDE][i] = false;
		bvisible_gran[6][i] = false; // ��� ������ �������!!!


		integer inode1, inode2, inode3, inode4, inode5, inode6, inode7, inode8;
		inode1 = nvtx[0][i] - 1;
		inode2 = nvtx[1][i] - 1;
		inode3 = nvtx[2][i] - 1;
		inode4 = nvtx[3][i] - 1;
		inode5 = nvtx[4][i] - 1;
		inode6 = nvtx[5][i] - 1;
		inode7 = nvtx[6][i] - 1;
		inode8 = nvtx[7][i] - 1;

		TOCHKA ptest;
		//center_cord3D(inode1, nvtx, pa, ptest, 100);
		ptest.x = pa_opengl[inode1].x;
		ptest.y = pa_opengl[inode1].y;
		ptest.z = pa_opengl[inode1].z;
		bool bvisible = false;

		if ((fabs(pfpir.fminimum) < 1.0e-30) && (fabs(pfpir.fmaximum) < 1.0e-30)) {
			pfpir.fminimum = -1.0e30;
			pfpir.fmaximum = 1.0e30;
		}

		switch (pfpir.idir) {
		case LINE_DIRECTIONAL::X_LINE_DIRECTIONAL:
			if ((ptest.x >= pfpir.fminimum) && (ptest.x <= pfpir.fmaximum)) {
				bvisible = true;
			}
			break;
		case LINE_DIRECTIONAL::Y_LINE_DIRECTIONAL:
			if ((ptest.y >= pfpir.fminimum) && (ptest.y <= pfpir.fmaximum)) {
				bvisible = true;
			}
			break;
		case LINE_DIRECTIONAL::Z_LINE_DIRECTIONAL:
			if ((ptest.z >= pfpir.fminimum) && (ptest.z <= pfpir.fmaximum)) {
				bvisible = true;
			}
			break;
		}

		if (bvisible) {

			//TOCHKA p1, p2, p3, p4, p5, p6, p7, p8;
			//TOCHKA pall;
			//center_cord3D(inode1, nvtx, pa, p1, 100);
			//center_cord3D(inode2, nvtx, pa, p2, 100);
			//center_cord3D(inode3, nvtx, pa, p3, 100);
			//center_cord3D(inode4, nvtx, pa, p4, 100);
			//center_cord3D(inode5, nvtx, pa, p5, 100);
			//center_cord3D(inode6, nvtx, pa, p6, 100);
			//center_cord3D(inode7, nvtx, pa, p7, 100);
			//center_cord3D(inode8, nvtx, pa, p8, 100);

			integer ib1, ib2, ib3, ib4, ib5, ib6, ib7, ib8;
			// ������������� ����������������� ��� ������� whot_is_block �����������
			// ������� � �� �������� � ����� �� �� �� ���� �������� �� ������� �����.
			ib1 = whot_is_block[inode1];
			//in_model_temp(p1, ib1, b, lb);	
			//if (ib1 != t.whot_is_block[inode1]) {
			//printf("ib1=%d whot_is_block=%d\n", ib1, t.whot_is_block[inode1]);
			//getchar();
			//}
			ib2 = whot_is_block[inode2];
			//in_model_temp(p2, ib2, b, lb);
			//if (ib2 != t.whot_is_block[inode2]) {
			//printf("ib2=%d whot_is_block=%d\n", ib2, t.whot_is_block[inode2]);
			//getchar();
			//}
			ib3 = whot_is_block[inode3];
			//in_model_temp(p3, ib3, b, lb);
			//if (ib3 != t.whot_is_block[inode3]) {
			//printf("ib3=%d whot_is_block=%d\n", ib3, t.whot_is_block[inode3]);
			//getchar();
			//}
			ib4 = whot_is_block[inode4];
			//in_model_temp(p4, ib4, b, lb);
			//if (ib4 != t.whot_is_block[inode4]) {
			//printf("ib4=%d whot_is_block=%d\n", ib4, t.whot_is_block[inode4]);
			//getchar();
			//}
			ib5 = whot_is_block[inode5];
			//in_model_temp(p5, ib5, b, lb);
			//if (ib5 != t.whot_is_block[inode5]) {
			//printf("ib5=%d whot_is_block=%d\n", ib5, t.whot_is_block[inode5]);
			//getchar();
			//}
			ib6 = whot_is_block[inode6];
			//in_model_temp(p6, ib6, b, lb);
			//if (ib6 != t.whot_is_block[inode6]) {
			//printf("ib6=%d whot_is_block=%d\n", ib6, t.whot_is_block[inode6]);
			//getchar();
			//}
			ib7 = whot_is_block[inode7];
			//in_model_temp(p7, ib7, b, lb);
			//if (ib7 != t.whot_is_block[inode7]) {
			//printf("ib7=%d whot_is_block=%d\n", ib7, t.whot_is_block[inode7]);
			//getchar();
			//}
			ib8 = whot_is_block[inode8];
			//in_model_temp(p8, ib8, b, lb);
			//if (ib8 != t.whot_is_block[inode8]) {
			//printf("ib8=%d whot_is_block=%d\n", ib8, t.whot_is_block[inode8]);
			//getchar();
			//}


			if ((b[ib1].bvisible) && (b[ib2].bvisible) && (b[ib3].bvisible) && (b[ib4].bvisible) && (b[ib5].bvisible) && (b[ib6].bvisible) && (b[ib7].bvisible) && (b[ib8].bvisible)) {


				for (int j = 0; j < NUMBER_OF_VERTEX_FINITE_ELEMENT(); j++) {
					//std::cout << nvtx[j][i] - 1 << " ";
					ipa_count[nvtx[j][i] - 1]++;
				}
			}
			//getchar();
		}
	}

	const int bVisibleCount = 7;

	for (int j = 0; j < maxelm; j++) {

		integer inode1, inode2, inode3, inode4, inode5, inode6, inode7, inode8;
		inode1 = nvtx[0][j] - 1;
		inode2 = nvtx[1][j] - 1;
		inode3 = nvtx[2][j] - 1;
		inode4 = nvtx[3][j] - 1;
		inode5 = nvtx[4][j] - 1;
		inode6 = nvtx[5][j] - 1;
		inode7 = nvtx[6][j] - 1;
		inode8 = nvtx[7][j] - 1;


		TOCHKA ptest;
		//center_cord3D(inode1, nvtx, pa, ptest, 100);
		ptest.x = pa_opengl[inode1].x;
		ptest.y = pa_opengl[inode1].y;
		ptest.z = pa_opengl[inode1].z;
		bool bvisible = false;

		switch (pfpir.idir) {
		case LINE_DIRECTIONAL::X_LINE_DIRECTIONAL:
			if ((ptest.x >= pfpir.fminimum) && (ptest.x <= pfpir.fmaximum)) {
				bvisible = true;
			}
			break;
		case LINE_DIRECTIONAL::Y_LINE_DIRECTIONAL:
			if ((ptest.y >= pfpir.fminimum) && (ptest.y <= pfpir.fmaximum)) {
				bvisible = true;
			}
			break;
		case LINE_DIRECTIONAL::Z_LINE_DIRECTIONAL:
			if ((ptest.z >= pfpir.fminimum) && (ptest.z <= pfpir.fmaximum)) {
				bvisible = true;
			}
			break;
		}

		if (bvisible) {

			//TOCHKA p1, p2, p3, p4, p5, p6, p7, p8;
			//TOCHKA pall;
			//center_cord3D(inode1, nvtx, pa, p1, 100);
			//center_cord3D(inode2, nvtx, pa, p2, 100);
			//center_cord3D(inode3, nvtx, pa, p3, 100);
			//center_cord3D(inode4, nvtx, pa, p4, 100);
			//center_cord3D(inode5, nvtx, pa, p5, 100);
			//center_cord3D(inode6, nvtx, pa, p6, 100);
			//center_cord3D(inode7, nvtx, pa, p7, 100);
			//center_cord3D(inode8, nvtx, pa, p8, 100);

			integer ib1, ib2, ib3, ib4, ib5, ib6, ib7, ib8;
			// ������������� ����������������� ��� ������� whot_is_block �����������
			// ������� � �� �������� � ����� �� �� �� ���� �������� �� ������� �����.
			ib1 = whot_is_block[inode1];
			//in_model_temp(p1, ib1, b, lb);	
			//if (ib1 != t.whot_is_block[inode1]) {
			//printf("ib1=%d whot_is_block=%d\n", ib1, t.whot_is_block[inode1]);
			//getchar();
			//}
			ib2 = whot_is_block[inode2];
			//in_model_temp(p2, ib2, b, lb);
			//if (ib2 != t.whot_is_block[inode2]) {
			//printf("ib2=%d whot_is_block=%d\n", ib2, t.whot_is_block[inode2]);
			//getchar();
			//}
			ib3 = whot_is_block[inode3];
			//in_model_temp(p3, ib3, b, lb);
			//if (ib3 != t.whot_is_block[inode3]) {
			//printf("ib3=%d whot_is_block=%d\n", ib3, t.whot_is_block[inode3]);
			//getchar();
			//}
			ib4 = whot_is_block[inode4];
			//in_model_temp(p4, ib4, b, lb);
			//if (ib4 != t.whot_is_block[inode4]) {
			//printf("ib4=%d whot_is_block=%d\n", ib4, t.whot_is_block[inode4]);
			//getchar();
			//}
			ib5 = whot_is_block[inode5];
			//in_model_temp(p5, ib5, b, lb);
			//if (ib5 != t.whot_is_block[inode5]) {
			//printf("ib5=%d whot_is_block=%d\n", ib5, t.whot_is_block[inode5]);
			//getchar();
			//}
			ib6 = whot_is_block[inode6];
			//in_model_temp(p6, ib6, b, lb);
			//if (ib6 != t.whot_is_block[inode6]) {
			//printf("ib6=%d whot_is_block=%d\n", ib6, t.whot_is_block[inode6]);
			//getchar();
			//}
			ib7 = whot_is_block[inode7];
			//in_model_temp(p7, ib7, b, lb);
			//if (ib7 != t.whot_is_block[inode7]) {
			//printf("ib7=%d whot_is_block=%d\n", ib7, t.whot_is_block[inode7]);
			//getchar();
			//}
			ib8 = whot_is_block[inode8];
			//in_model_temp(p8, ib8, b, lb);
			//if (ib8 != t.whot_is_block[inode8]) {
			//printf("ib8=%d whot_is_block=%d\n", ib8, t.whot_is_block[inode8]);
			//getchar();
			//}


			if ((b[ib1].bvisible) && (b[ib2].bvisible) && (b[ib3].bvisible) && (b[ib4].bvisible) && (b[ib5].bvisible) && (b[ib6].bvisible) && (b[ib7].bvisible) && (b[ib8].bvisible)) {

				bvisible_gran[6][j] = true;

				// XY bottom
				if ((ipa_count[inode1] <= bVisibleCount) &&
					(ipa_count[inode2] <= bVisibleCount) &&
					(ipa_count[inode3] <= bVisibleCount) &&
					(ipa_count[inode4] <= bVisibleCount)) {

					bvisible_gran[B_SIDE][j] = true;
				}

				// XY Top
				if ((ipa_count[inode5] <= bVisibleCount) &&
					(ipa_count[inode6] <= bVisibleCount) &&
					(ipa_count[inode7] <= bVisibleCount) &&
					(ipa_count[inode8] <= bVisibleCount)) {

					bvisible_gran[T_SIDE][j] = true;
				}

				// XZ SSIDE min Y
				if ((ipa_count[inode1] <= bVisibleCount) &&
					(ipa_count[inode2] <= bVisibleCount) &&
					(ipa_count[inode6] <= bVisibleCount) &&
					(ipa_count[inode5] <= bVisibleCount)) {

					bvisible_gran[S_SIDE][j] = true;
				}

				// XZ SSIDE max Y
				if ((ipa_count[inode4] <= bVisibleCount) &&
					(ipa_count[inode8] <= bVisibleCount) &&
					(ipa_count[inode7] <= bVisibleCount) &&
					(ipa_count[inode3] <= bVisibleCount)) {

					bvisible_gran[N_SIDE][j] = true;
				}

				// YZ SSIDE max X
				if ((ipa_count[inode2] <= bVisibleCount) &&
					(ipa_count[inode3] <= bVisibleCount) &&
					(ipa_count[inode7] <= bVisibleCount) &&
					(ipa_count[inode6] <= bVisibleCount)) {

					bvisible_gran[E_SIDE][j] = true;
				}

				// YZ SSIDE min X
				if ((ipa_count[inode1] <= bVisibleCount) &&
					(ipa_count[inode5] <= bVisibleCount) &&
					(ipa_count[inode8] <= bVisibleCount) &&
					(ipa_count[inode4] <= bVisibleCount))
				{
					bvisible_gran[W_SIDE][j] = true;
				}

			}
		}
	}

	GLFWwindow* window;

	/* Initialize the library */
	if (!glfwInit())
		return -1;

	/* Create a windowed mode window and its OpenGL context */
	window = glfwCreateWindow(SCREEN_WIDTH, SCREEN_HEIGHT, "AliceFlow_v0_49", NULL, NULL);


	glfwSetKeyCallback(window, keyCallback);
	glfwSetInputMode(window, GLFW_STICKY_KEYS, 1);

	int screenWidth, screenHeight;
	glfwGetFramebufferSize(window, &screenWidth, &screenHeight);

	if (!window)
	{
		glfwTerminate();
		return -1;
	}

	/* Make the window's context current */
	glfwMakeContextCurrent(window);


	glViewport(0.0f, 0.0f, screenWidth, screenHeight); // specifies the part of the window to which OpenGL will draw (in pixels), convert from
													   // normalised to pixels
	glMatrixMode(GL_PROJECTION); // projection matrix defines the properties of the camera that views the objects in the world coordinate frame.
								 // Here you thalfScreenHeight+scale_all*pa.yally set the zoom factor, aspect ratio and the near and far clipping planes
	glLoadIdentity(); // replace the current matrix with the identity matrix and starts us a fresh because matrix transforms such as glOrtho and
					  // glRotate cumulate, besically puts us at (0, 0, 0)
	glOrtho(0, SCREEN_WIDTH, 0, SCREEN_HEIGHT, 0.1, 10000); // essentially set coordinate system
	glMatrixMode(GL_MODELVIEW); // (default matrix mode) modelview matrix defines how objects are transformed (meaning translation, rotation
								// and scaling) in your world

	glLoadIdentity(); // same as above comment 



	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LESS);
	glEnable(GL_LIGHTING); // ��������� ����.
	glEnable(GL_LIGHT0); // �������� ����� ������ 0.
	glEnable(GL_COLOR_MATERIAL); // ����������� ����� �������� ��������.

								 /* Loop until the user closes the window */
	while (!glfwWindowShouldClose(window))
	{


		if (binverse_color_black_2_white) {
			glClearColor(1.0f, 1.0f, 1.0f, 0.0f);// ����� ���
		}
		//glColor3f(0.0, 0.0, 0.0); // ׸����
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		/* Render OpenGL here */

		glPushMatrix();

		float position[] = { 0,0,1,0 };
		glLightfv(GL_LIGHT0, GL_POSITION, position);

		glTranslatef(halfScreenWidth, halfScreenHeight, -abbys);
		glRotatef(rotationX, 1, 0, 0);
		glRotatef(rotationY, 0, 1, 0);
		glTranslatef(-halfScreenWidth, -halfScreenHeight, abbys);



		glLineWidth(1);

		// �������� ���.
		for (int j = 0; j < maxelm; j++) {

			if (bvisible_gran[6][j]) {

				int inode1, inode2, inode3, inode4, inode5, inode6, inode7, inode8;
				inode1 = nvtx[0][j] - 1;
				inode2 = nvtx[1][j] - 1;
				inode3 = nvtx[2][j] - 1;
				inode4 = nvtx[3][j] - 1;
				inode5 = nvtx[4][j] - 1;
				inode6 = nvtx[5][j] - 1;
				inode7 = nvtx[6][j] - 1;
				inode8 = nvtx[7][j] - 1;


				// XY bottom
				if (bvisible_gran[B_SIDE][j])
				{

					if (binvisible_face_detect(0.0, 0.0, -1.0)) {


						//glColor3f(0.84,0.84,0.84);

						//glColor3f(0.0, 0.0, 0.0); // ׸����


						//Drawnvtx(nvtx, pa, j, halfScreenWidth, halfScreenHeight, -abbys, scale_all);


						glBegin(GL_QUADS);

						set_color_for_render(potent[inode4]);
						//glColor3f(1.0, 1.0, 1.0);
						glNormal3f(0.0, 0.0, -1.0);
						glVertex3f(pa_render[inode4].x, pa_render[inode4].y, pa_render[inode4].z);
						set_color_for_render(potent[inode3]);
						//glColor3f(1.0, 1.0, 1.0);
						glNormal3f(0.0, 0.0, -1.0);
						glVertex3f(pa_render[inode3].x, pa_render[inode3].y, pa_render[inode3].z);
						set_color_for_render(potent[inode2]);
						//glColor3f(1.0, 1.0, 1.0);
						glNormal3f(0.0, 0.0, -1.0);
						glVertex3f(pa_render[inode2].x, pa_render[inode2].y, pa_render[inode2].z);
						set_color_for_render(potent[inode1]);
						//glColor3f(1.0, 1.0, 1.0);
						glNormal3f(0.0, 0.0, -1.0);
						glVertex3f(pa_render[inode1].x, pa_render[inode1].y, pa_render[inode1].z);
						//glColor3f(0.0, 0.0, 0.0);
						glEnd();


						/*
						if (binverse_color_black_2_white) {
						glColor3f(0.0, 0.0, 0.0); // ׸����
						}
						else {
						glColor3f(1.0, 1.0, 1.0);// �����
						}

						if (CheckLine(ipa_count[inode1],
						ipa_count[inode2])
						)
						{

						glBegin(GL_LINE_LOOP);
						glVertex3f(pa_render[inode1].x, pa_render[inode1].y, pa_render[inode1].z);
						glVertex3f(pa_render[inode2].x, pa_render[inode2].y, pa_render[inode2].z);
						glEnd();
						}

						if (CheckLine(ipa_count[inode3],
						ipa_count[inode2])
						)
						{

						glBegin(GL_LINE_LOOP);
						glVertex3f(pa_render[inode2].x, pa_render[inode2].y, pa_render[inode2].z);
						glVertex3f(pa_render[inode3].x, pa_render[inode3].y, pa_render[inode3].z);
						glEnd();
						}

						if (CheckLine(ipa_count[inode3],
						ipa_count[inode4])
						)
						{
						glBegin(GL_LINE_LOOP);
						glVertex3f(pa_render[inode3].x, pa_render[inode3].y, pa_render[inode3].z);
						glVertex3f(pa_render[inode4].x, pa_render[inode4].y, pa_render[inode4].z);
						glEnd();
						}

						if (CheckLine(ipa_count[inode1],
						ipa_count[inode4])
						)
						{
						glBegin(GL_LINE_LOOP);
						glVertex3f(pa_render[inode4].x, pa_render[inode4].y, pa_render[inode4].z);
						glVertex3f(pa_render[inode1].x, pa_render[inode1].y, pa_render[inode1].z);
						glEnd();
						}

						//glLineWidth(1);
						//glColor3f(1.0, 1.0, 1.0);
						//glLineWidth(1);
						*/

					}
				}


				// XY Top
				if (bvisible_gran[T_SIDE][j])
				{


					if (binvisible_face_detect(0.0, 0.0, 1.0)) {


						//glColor3f(0.84,0.84,0.84);
						//glColor3f(0.0, 0.0, 0.0); // ׸����

						//Drawnvtx(nvtx, pa, j, halfScreenWidth, halfScreenHeight, -abbys, scale_all);


						glBegin(GL_QUADS);
						set_color_for_render(potent[inode8]);
						//glColor3f(1.0, 1.0, 1.0);
						glNormal3f(0.0, 0.0, 1.0);
						glVertex3f(pa_render[inode8].x, pa_render[inode8].y, pa_render[inode8].z);
						set_color_for_render(potent[inode7]);
						//glColor3f(1.0, 1.0, 1.0);
						glNormal3f(0.0, 0.0, 1.0);
						glVertex3f(pa_render[inode7].x, pa_render[inode7].y, pa_render[inode7].z);
						set_color_for_render(potent[inode6]);
						//glColor3f(1.0, 1.0, 1.0);
						glNormal3f(0.0, 0.0, 1.0);
						glVertex3f(pa_render[inode6].x, pa_render[inode6].y, pa_render[inode6].z);
						set_color_for_render(potent[inode5]);
						//glColor3f(1.0, 1.0, 1.0);
						glNormal3f(0.0, 0.0, 1.0);
						glVertex3f(pa_render[inode5].x, pa_render[inode5].y, pa_render[inode5].z);
						//glColor3f(0.0, 0.0, 0.0);
						glEnd();

						/*
						//glLineWidth(ComboBoxlineWidth.ItemIndex + 1);
						if (binverse_color_black_2_white) {
						glColor3f(0.0, 0.0, 0.0); // ׸����
						}
						else {
						glColor3f(1.0, 1.0, 1.0);// �����
						}

						if (CheckLine(ipa_count[inode5],
						ipa_count[inode6])
						)
						{

						glBegin(GL_LINE_LOOP);
						glVertex3f(pa_render[inode5].x, pa_render[inode5].y, pa_render[inode5].z);
						glVertex3f(pa_render[inode6].x, pa_render[inode6].y, pa_render[inode6].z);
						glEnd();
						}

						if (CheckLine(ipa_count[inode6],
						ipa_count[inode7])
						) {

						glBegin(GL_LINE_LOOP);
						glVertex3f(pa_render[inode6].x, pa_render[inode6].y, pa_render[inode6].z);
						glVertex3f(pa_render[inode7].x, pa_render[inode7].y, pa_render[inode7].z);
						glEnd();
						}

						if (CheckLine(ipa_count[inode7],
						ipa_count[inode8])
						)
						{
						glBegin(GL_LINE_LOOP);
						glVertex3f(pa_render[inode7].x, pa_render[inode7].y, pa_render[inode7].z);
						glVertex3f(pa_render[inode8].x, pa_render[inode8].y, pa_render[inode8].z);
						glEnd();
						}

						if (CheckLine(ipa_count[inode5],
						ipa_count[inode8])
						)
						{
						glBegin(GL_LINE_LOOP);
						glVertex3f(pa_render[inode8].x, pa_render[inode8].y, pa_render[inode8].z);
						glVertex3f(pa_render[inode5].x, pa_render[inode5].y, pa_render[inode5].z);
						glEnd();
						}

						// glLineWidth(1);
						//glColor3f(1.0, 1.0, 1.0);
						//glLineWidth(1);
						*/
					}
				}



				// XZ SSIDE min Y
				if (bvisible_gran[S_SIDE][j])
				{


					if (binvisible_face_detect(0.0, -1.0, 0.0)) {

						//glColor3f(0.84,0.84,0.84);


						//Drawnvtx(nvtx, pa, j, halfScreenWidth, halfScreenHeight, -abbys, scale_all);


						glBegin(GL_QUADS);
						set_color_for_render(potent[inode5]);
						//glColor3f(1.0, 1.0, 1.0);
						glNormal3f(0.0, -1.0, 0.0);
						glVertex3f(pa_render[inode5].x, pa_render[inode5].y, pa_render[inode5].z);
						set_color_for_render(potent[inode1]);
						//glColor3f(1.0, 1.0, 1.0);
						glNormal3f(0.0, -1.0, 0.0);
						glVertex3f(pa_render[inode1].x, pa_render[inode1].y, pa_render[inode1].z);
						set_color_for_render(potent[inode2]);
						//glColor3f(1.0, 1.0, 1.0);
						glNormal3f(0.0, -1.0, 0.0);
						glVertex3f(pa_render[inode2].x, pa_render[inode2].y, pa_render[inode2].z);
						set_color_for_render(potent[inode6]);
						//glColor3f(1.0, 1.0, 1.0);
						glNormal3f(0.0, -1.0, 0.0);
						glVertex3f(pa_render[inode6].x, pa_render[inode6].y, pa_render[inode6].z);
						//glColor3f(0.0, 0.0, 0.0);
						glEnd();

						/*
						//glLineWidth(ComboBoxlineWidth.ItemIndex + 1);
						if (binverse_color_black_2_white) {
						glColor3f(0.0, 0.0, 0.0); // ׸����
						}
						else {
						glColor3f(1.0, 1.0, 1.0);// �����
						}

						if (CheckLine(ipa_count[inode1],
						ipa_count[inode2])
						) {

						glBegin(GL_LINE_LOOP);
						glVertex3f(pa_render[inode1].x, pa_render[inode1].y, pa_render[inode1].z);
						glVertex3f(pa_render[inode2].x, pa_render[inode2].y, pa_render[inode2].z);
						glEnd();
						}

						if (CheckLine(ipa_count[inode2],
						ipa_count[inode6])
						) {

						glBegin(GL_LINE_LOOP);
						glVertex3f(pa_render[inode2].x, pa_render[inode2].y, pa_render[inode2].z);
						glVertex3f(pa_render[inode6].x, pa_render[inode6].y, pa_render[inode6].z);
						glEnd();
						}

						if (CheckLine(ipa_count[inode6],
						ipa_count[inode5])
						) {

						glBegin(GL_LINE_LOOP);
						glVertex3f(pa_render[inode6].x, pa_render[inode6].y, pa_render[inode6].z);
						glVertex3f(pa_render[inode5].x, pa_render[inode5].y, pa_render[inode5].z);
						glEnd();
						}

						if (CheckLine(ipa_count[inode1],
						ipa_count[inode5])
						) {

						glBegin(GL_LINE_LOOP);
						glVertex3f(pa_render[inode5].x, pa_render[inode5].y, pa_render[inode5].z);
						glVertex3f(pa_render[inode1].x, pa_render[inode1].y, pa_render[inode1].z);
						glEnd();
						}

						//glLineWidth(1);
						//glColor3f(1.0, 1.0, 1.0);
						//glLineWidth(1);
						*/
					}
				}


				// XZ SSIDE max Y
				if (bvisible_gran[N_SIDE][j])
				{

					if (binvisible_face_detect(0.0, 1.0, 0.0)) {

						//glColor3f(0.84,0.84,0.84);
						//glColor3f(0.0, 0.0, 0.0); // ׸����

						//Drawnvtx(nvtx, pa, j, halfScreenWidth, halfScreenHeight, -abbys, scale_all);


						glBegin(GL_QUADS);
						set_color_for_render(potent[inode8]);
						//glColor3f(1.0, 1.0, 1.0);
						glNormal3f(0.0, 1.0, 0.0);
						glVertex3f(pa_render[inode8].x, pa_render[inode8].y, pa_render[inode8].z);
						set_color_for_render(potent[inode4]);
						//glColor3f(1.0, 1.0, 1.0);
						glNormal3f(0.0, 1.0, 0.0);
						glVertex3f(pa_render[inode4].x, pa_render[inode4].y, pa_render[inode4].z);
						set_color_for_render(potent[inode3]);
						//glColor3f(1.0, 1.0, 1.0);
						glNormal3f(0.0, 1.0, 0.0);
						glVertex3f(pa_render[inode3].x, pa_render[inode3].y, pa_render[inode3].z);
						set_color_for_render(potent[inode7]);
						//glColor3f(1.0, 1.0, 1.0);
						glNormal3f(0.0, 1.0, 0.0);
						glVertex3f(pa_render[inode7].x, pa_render[inode7].y, pa_render[inode7].z);

						//glColor3f(0.0, 0.0, 0.0);
						glEnd();

						/*
						//glLineWidth(ComboBoxlineWidth.ItemIndex + 1);
						if (binverse_color_black_2_white) {
						glColor3f(0.0, 0.0, 0.0); // ׸����
						}
						else {
						glColor3f(1.0, 1.0, 1.0);// �����
						}

						if (CheckLine(ipa_count[inode3],
						ipa_count[inode4])
						) {

						glBegin(GL_LINE_LOOP);
						glVertex3f(pa_render[inode4].x, pa_render[inode4].y, pa_render[inode4].z);
						glVertex3f(pa_render[inode3].x, pa_render[inode3].y, pa_render[inode3].z);
						glEnd();
						}

						if (CheckLine(ipa_count[inode3],
						ipa_count[inode7])
						) {

						glBegin(GL_LINE_LOOP);
						glVertex3f(pa_render[inode3].x, pa_render[inode3].y, pa_render[inode3].z);
						glVertex3f(pa_render[inode7].x, pa_render[inode7].y, pa_render[inode7].z);
						glEnd();
						}

						if (CheckLine(ipa_count[inode7],
						ipa_count[inode8])
						) {
						glBegin(GL_LINE_LOOP);
						glVertex3f(pa_render[inode7].x, pa_render[inode7].y, pa_render[inode7].z);
						glVertex3f(pa_render[inode8].x, pa_render[inode8].y, pa_render[inode8].z);
						glEnd();
						}

						if (CheckLine(ipa_count[inode8],
						ipa_count[inode4])
						) {
						glBegin(GL_LINE_LOOP);
						glVertex3f(pa_render[inode8].x, pa_render[inode8].y, pa_render[inode8].z);
						glVertex3f(pa_render[inode4].x, pa_render[inode4].y, pa_render[inode4].z);
						glEnd();
						}

						// glLineWidth(1);
						//glColor3f(1.0, 1.0, 1.0);
						//glLineWidth(1);
						*/
					}
				}



				// YZ SSIDE max X
				if (bvisible_gran[E_SIDE][j])
				{

					if (binvisible_face_detect(1.0, 0.0, 0.0)) {


						//glColor3f(0.84,0.84,0.84);
						//glColor3f(0.0, 0.0, 0.0); // ׸����

						//Drawnvtx(nvtx, pa, j, halfScreenWidth, halfScreenHeight, -abbys, scale_all);


						glBegin(GL_QUADS);
						set_color_for_render(potent[inode6]);
						//glColor3f(1.0, 1.0, 1.0);
						glNormal3f(1.0, 0.0, 0.0);
						glVertex3f(pa_render[inode6].x, pa_render[inode6].y, pa_render[inode6].z);
						set_color_for_render(potent[inode2]);
						//glColor3f(1.0, 1.0, 1.0);
						glNormal3f(1.0, 0.0, 0.0);
						glVertex3f(pa_render[inode2].x, pa_render[inode2].y, pa_render[inode2].z);
						set_color_for_render(potent[inode3]);
						//glColor3f(1.0, 1.0, 1.0);
						glNormal3f(1.0, 0.0, 0.0);
						glVertex3f(pa_render[inode3].x, pa_render[inode3].y, pa_render[inode3].z);
						set_color_for_render(potent[inode7]);
						//glColor3f(1.0, 1.0, 1.0);
						glNormal3f(1.0, 0.0, 0.0);
						glVertex3f(pa_render[inode7].x, pa_render[inode7].y, pa_render[inode7].z);
						//glColor3f(0.0, 0.0, 0.0);
						glEnd();

						/*
						//glLineWidth(ComboBoxlineWidth.ItemIndex + 1);
						if (binverse_color_black_2_white) {
						glColor3f(0.0, 0.0, 0.0); // ׸����
						}
						else {
						glColor3f(1.0, 1.0, 1.0);// �����
						}


						if (CheckLine(ipa_count[inode3],
						ipa_count[inode2])
						)
						{
						glBegin(GL_LINE_LOOP);
						glVertex3f(pa_render[inode2].x, pa_render[inode2].y, pa_render[inode2].z);
						glVertex3f(pa_render[inode3].x, pa_render[inode3].y, pa_render[inode3].z);
						glEnd();
						}

						if (CheckLine(ipa_count[inode3],
						ipa_count[inode7])
						)
						{
						glBegin(GL_LINE_LOOP);
						glVertex3f(pa_render[inode3].x, pa_render[inode3].y, pa_render[inode3].z);
						glVertex3f(pa_render[inode7].x, pa_render[inode7].y, pa_render[inode7].z);
						glEnd();
						}

						if (CheckLine(ipa_count[inode7],
						ipa_count[inode6])
						)
						{
						glBegin(GL_LINE_LOOP);
						glVertex3f(pa_render[inode7].x, pa_render[inode7].y, pa_render[inode7].z);
						glVertex3f(pa_render[inode6].x, pa_render[inode6].y, pa_render[inode6].z);
						glEnd();
						}

						if (CheckLine(ipa_count[inode6],
						ipa_count[inode2])
						)
						{
						glBegin(GL_LINE_LOOP);
						glVertex3f(pa_render[inode6].x, pa_render[inode6].y, pa_render[inode6].z);
						glVertex3f(pa_render[inode2].x, pa_render[inode2].y, pa_render[inode2].z);
						glEnd();
						}

						//glLineWidth(1);
						//glColor3f(1.0, 1.0, 1.0);
						//glLineWidth(1);
						*/
					}
				}


				// YZ SSIDE min X
				if (bvisible_gran[W_SIDE][j])
				{

					if (binvisible_face_detect(-1.0, 0.0, 0.0))
					{

						//glColor3f(0.84,0.84,0.84);

						//glColor3f(0.0, 0.0, 0.0); // ׸����

						//Drawnvtx(nvtx, pa, j, halfScreenWidth, halfScreenHeight, -abbys, scale_all);


						glBegin(GL_QUADS);
						set_color_for_render(potent[inode8]);
						//glColor3f(1.0, 1.0, 1.0);
						glNormal3f(-1.0, 0.0, 0.0);
						glVertex3f(pa_render[inode8].x, pa_render[inode8].y, pa_render[inode8].z);
						set_color_for_render(potent[inode4]);
						//glColor3f(1.0, 1.0, 1.0);
						glNormal3f(-1.0, 0.0, 0.0);
						glVertex3f(pa_render[inode4].x, pa_render[inode4].y, pa_render[inode4].z);
						set_color_for_render(potent[inode1]);
						//glColor3f(1.0, 1.0, 1.0);
						glNormal3f(-1.0, 0.0, 0.0);
						glVertex3f(pa_render[inode1].x, pa_render[inode1].y, pa_render[inode1].z);
						set_color_for_render(potent[inode5]);
						//glColor3f(1.0, 1.0, 1.0);
						glNormal3f(-1.0, 0.0, 0.0);
						glVertex3f(pa_render[inode5].x, pa_render[inode5].y, pa_render[inode5].z);
						//glColor3f(0.0, 0.0, 0.0);
						glEnd();


						/*
						//glLineWidth(ComboBoxlineWidth.ItemIndex + 1);

						if (binverse_color_black_2_white) {
						glColor3f(0.0, 0.0, 0.0); // ׸����
						}
						else {
						glColor3f(1.0, 1.0, 1.0);// �����
						}

						if (CheckLine(ipa_count[inode1],
						ipa_count[inode4])
						) {

						glBegin(GL_LINE_LOOP);
						glVertex3f(pa_render[inode1].x, pa_render[inode1].y, pa_render[inode1].z);
						glVertex3f(pa_render[inode4].x, pa_render[inode4].y, pa_render[inode4].z);
						glEnd();
						}

						if (CheckLine(ipa_count[inode4],
						ipa_count[inode8])
						)
						{
						glBegin(GL_LINE_LOOP);
						glVertex3f(pa_render[inode4].x, pa_render[inode4].y, pa_render[inode4].z);
						glVertex3f(pa_render[inode8].x, pa_render[inode8].y, pa_render[inode8].z);
						glEnd();
						}

						if (CheckLine(ipa_count[inode5],
						ipa_count[inode8])
						)
						{
						glBegin(GL_LINE_LOOP);
						glVertex3f(pa_render[inode8].x, pa_render[inode8].y, pa_render[inode8].z);
						glVertex3f(pa_render[inode5].x, pa_render[inode5].y, pa_render[inode5].z);
						glEnd();
						}

						if (CheckLine(ipa_count[inode1],
						ipa_count[inode5])
						)
						{
						glBegin(GL_LINE_LOOP);
						glVertex3f(pa_render[inode5].x, pa_render[inode5].y, pa_render[inode5].z);
						glVertex3f(pa_render[inode1].x, pa_render[inode1].y, pa_render[inode1].z);
						glEnd();
						}

						//glLineWidth(1);
						//glColor3f(1.0, 1.0, 1.0);
						//glLineWidth(1);
						*/
					}
				}

			}
		}


		glPopMatrix();

		/* Swap front and back buffers */
		glfwSwapBuffers(window);

		/* Poll for and process events */
		glfwPollEvents();
	}

	delete[] ipa_count;

	for (int i = 0; i < 7; i++) {
		delete[] bvisible_gran[i];
	}
	delete[] bvisible_gran;

	delete[] checkline1;
	delete[] checkline2;

	glfwTerminate();

	return 0;

} // DrawZbufferColorLight


  // ��������� ������������ � ��������� ��������� �����. ��������� 19.12.2020.
int DrawZbufferColor(BLOCK*& b, integer lb, int maxelm, int maxnod, int**& nvtx, int*& whot_is_block, doublereal* &potent) {

	int CHECK_LINE_SIZE = 9;
	checkline1 = new bool[CHECK_LINE_SIZE * CHECK_LINE_SIZE];
	checkline2 = new bool[CHECK_LINE_SIZE * CHECK_LINE_SIZE];
	for (int i = 0; i < CHECK_LINE_SIZE * CHECK_LINE_SIZE; i++) {
		checkline1[i] = false;
		checkline2[i] = false;
	}

	for (int i = 0; i < CHECK_LINE_SIZE; i++) {
		for (int j = 0; j < CHECK_LINE_SIZE; j++) {
			int key = i * CHECK_LINE_SIZE + j;

			checkline1[key] = true;
			checkline2[key] = CheckLine1(i, j);
		}
	}

	bool** bvisible_gran = new bool*[7];
	for (int i = 0; i < 7; i++) {
		bvisible_gran[i] = new bool[maxelm];
	}

	int* ipa_count = new int[maxnod];
	for (int i = 0; i < maxnod; i++) {
		ipa_count[i] = 0;
	}
	for (int i = 0; i < maxelm; i++) {

		bvisible_gran[E_SIDE][i] = false;
		bvisible_gran[W_SIDE][i] = false;
		bvisible_gran[N_SIDE][i] = false;
		bvisible_gran[S_SIDE][i] = false;
		bvisible_gran[T_SIDE][i] = false;
		bvisible_gran[B_SIDE][i] = false;
		bvisible_gran[6][i] = false; // ��� ������ �������!!!


		integer inode1, inode2, inode3, inode4, inode5, inode6, inode7, inode8;
		inode1 = nvtx[0][i] - 1;
		inode2 = nvtx[1][i] - 1;
		inode3 = nvtx[2][i] - 1;
		inode4 = nvtx[3][i] - 1;
		inode5 = nvtx[4][i] - 1;
		inode6 = nvtx[5][i] - 1;
		inode7 = nvtx[6][i] - 1;
		inode8 = nvtx[7][i] - 1;

		TOCHKA ptest;
		//center_cord3D(inode1, nvtx, pa, ptest, 100);
		ptest.x = pa_opengl[inode1].x;
		ptest.y = pa_opengl[inode1].y;
		ptest.z = pa_opengl[inode1].z;
		bool bvisible = false;

		if ((fabs(pfpir.fminimum) < 1.0e-30) && (fabs(pfpir.fmaximum) < 1.0e-30)) {
			pfpir.fminimum = -1.0e30;
			pfpir.fmaximum = 1.0e30;
		}

		switch (pfpir.idir) {
		case LINE_DIRECTIONAL::X_LINE_DIRECTIONAL:
			if ((ptest.x >= pfpir.fminimum) && (ptest.x <= pfpir.fmaximum)) {
				bvisible = true;
			}
			break;
		case LINE_DIRECTIONAL::Y_LINE_DIRECTIONAL:
			if ((ptest.y >= pfpir.fminimum) && (ptest.y <= pfpir.fmaximum)) {
				bvisible = true;
			}
			break;
		case LINE_DIRECTIONAL::Z_LINE_DIRECTIONAL:
			if ((ptest.z >= pfpir.fminimum) && (ptest.z <= pfpir.fmaximum)) {
				bvisible = true;
			}
			break;
		}

		if (bvisible) {

			//TOCHKA p1, p2, p3, p4, p5, p6, p7, p8;
			//TOCHKA pall;
			//center_cord3D(inode1, nvtx, pa, p1, 100);
			//center_cord3D(inode2, nvtx, pa, p2, 100);
			//center_cord3D(inode3, nvtx, pa, p3, 100);
			//center_cord3D(inode4, nvtx, pa, p4, 100);
			//center_cord3D(inode5, nvtx, pa, p5, 100);
			//center_cord3D(inode6, nvtx, pa, p6, 100);
			//center_cord3D(inode7, nvtx, pa, p7, 100);
			//center_cord3D(inode8, nvtx, pa, p8, 100);

			integer ib1, ib2, ib3, ib4, ib5, ib6, ib7, ib8;
			// ������������� ����������������� ��� ������� whot_is_block �����������
			// ������� � �� �������� � ����� �� �� �� ���� �������� �� ������� �����.
			ib1 = whot_is_block[inode1];
			//in_model_temp(p1, ib1, b, lb);	
			//if (ib1 != t.whot_is_block[inode1]) {
			//printf("ib1=%d whot_is_block=%d\n", ib1, t.whot_is_block[inode1]);
			//getchar();
			//}
			ib2 = whot_is_block[inode2];
			//in_model_temp(p2, ib2, b, lb);
			//if (ib2 != t.whot_is_block[inode2]) {
			//printf("ib2=%d whot_is_block=%d\n", ib2, t.whot_is_block[inode2]);
			//getchar();
			//}
			ib3 = whot_is_block[inode3];
			//in_model_temp(p3, ib3, b, lb);
			//if (ib3 != t.whot_is_block[inode3]) {
			//printf("ib3=%d whot_is_block=%d\n", ib3, t.whot_is_block[inode3]);
			//getchar();
			//}
			ib4 = whot_is_block[inode4];
			//in_model_temp(p4, ib4, b, lb);
			//if (ib4 != t.whot_is_block[inode4]) {
			//printf("ib4=%d whot_is_block=%d\n", ib4, t.whot_is_block[inode4]);
			//getchar();
			//}
			ib5 = whot_is_block[inode5];
			//in_model_temp(p5, ib5, b, lb);
			//if (ib5 != t.whot_is_block[inode5]) {
			//printf("ib5=%d whot_is_block=%d\n", ib5, t.whot_is_block[inode5]);
			//getchar();
			//}
			ib6 = whot_is_block[inode6];
			//in_model_temp(p6, ib6, b, lb);
			//if (ib6 != t.whot_is_block[inode6]) {
			//printf("ib6=%d whot_is_block=%d\n", ib6, t.whot_is_block[inode6]);
			//getchar();
			//}
			ib7 = whot_is_block[inode7];
			//in_model_temp(p7, ib7, b, lb);
			//if (ib7 != t.whot_is_block[inode7]) {
			//printf("ib7=%d whot_is_block=%d\n", ib7, t.whot_is_block[inode7]);
			//getchar();
			//}
			ib8 = whot_is_block[inode8];
			//in_model_temp(p8, ib8, b, lb);
			//if (ib8 != t.whot_is_block[inode8]) {
			//printf("ib8=%d whot_is_block=%d\n", ib8, t.whot_is_block[inode8]);
			//getchar();
			//}


			if ((b[ib1].bvisible) && (b[ib2].bvisible) && (b[ib3].bvisible) && (b[ib4].bvisible) && (b[ib5].bvisible) && (b[ib6].bvisible) && (b[ib7].bvisible) && (b[ib8].bvisible)) {


				for (int j = 0; j < NUMBER_OF_VERTEX_FINITE_ELEMENT(); j++) {
					//std::cout << nvtx[j][i] - 1 << " ";
					ipa_count[nvtx[j][i] - 1]++;
				}
			}
			//getchar();
		}
	}

	const int bVisibleCount = 7;

	for (int j = 0; j < maxelm; j++) {

		integer inode1, inode2, inode3, inode4, inode5, inode6, inode7, inode8;
		inode1 = nvtx[0][j] - 1;
		inode2 = nvtx[1][j] - 1;
		inode3 = nvtx[2][j] - 1;
		inode4 = nvtx[3][j] - 1;
		inode5 = nvtx[4][j] - 1;
		inode6 = nvtx[5][j] - 1;
		inode7 = nvtx[6][j] - 1;
		inode8 = nvtx[7][j] - 1;


		TOCHKA ptest;
		//center_cord3D(inode1, nvtx, pa, ptest, 100);
		ptest.x = pa_opengl[inode1].x;
		ptest.y = pa_opengl[inode1].y;
		ptest.z = pa_opengl[inode1].z;
		bool bvisible = false;

		switch (pfpir.idir) {
		case LINE_DIRECTIONAL::X_LINE_DIRECTIONAL:
			if ((ptest.x >= pfpir.fminimum) && (ptest.x <= pfpir.fmaximum)) {
				bvisible = true;
			}
			break;
		case LINE_DIRECTIONAL::Y_LINE_DIRECTIONAL:
			if ((ptest.y >= pfpir.fminimum) && (ptest.y <= pfpir.fmaximum)) {
				bvisible = true;
			}
			break;
		case LINE_DIRECTIONAL::Z_LINE_DIRECTIONAL:
			if ((ptest.z >= pfpir.fminimum) && (ptest.z <= pfpir.fmaximum)) {
				bvisible = true;
			}
			break;
		}

		if (bvisible) {

			//TOCHKA p1, p2, p3, p4, p5, p6, p7, p8;
			//TOCHKA pall;
			//center_cord3D(inode1, nvtx, pa, p1, 100);
			//center_cord3D(inode2, nvtx, pa, p2, 100);
			//center_cord3D(inode3, nvtx, pa, p3, 100);
			//center_cord3D(inode4, nvtx, pa, p4, 100);
			//center_cord3D(inode5, nvtx, pa, p5, 100);
			//center_cord3D(inode6, nvtx, pa, p6, 100);
			//center_cord3D(inode7, nvtx, pa, p7, 100);
			//center_cord3D(inode8, nvtx, pa, p8, 100);

			integer ib1, ib2, ib3, ib4, ib5, ib6, ib7, ib8;
			// ������������� ����������������� ��� ������� whot_is_block �����������
			// ������� � �� �������� � ����� �� �� �� ���� �������� �� ������� �����.
			ib1 = whot_is_block[inode1];
			//in_model_temp(p1, ib1, b, lb);	
			//if (ib1 != t.whot_is_block[inode1]) {
			//printf("ib1=%d whot_is_block=%d\n", ib1, t.whot_is_block[inode1]);
			//getchar();
			//}
			ib2 = whot_is_block[inode2];
			//in_model_temp(p2, ib2, b, lb);
			//if (ib2 != t.whot_is_block[inode2]) {
			//printf("ib2=%d whot_is_block=%d\n", ib2, t.whot_is_block[inode2]);
			//getchar();
			//}
			ib3 = whot_is_block[inode3];
			//in_model_temp(p3, ib3, b, lb);
			//if (ib3 != t.whot_is_block[inode3]) {
			//printf("ib3=%d whot_is_block=%d\n", ib3, t.whot_is_block[inode3]);
			//getchar();
			//}
			ib4 = whot_is_block[inode4];
			//in_model_temp(p4, ib4, b, lb);
			//if (ib4 != t.whot_is_block[inode4]) {
			//printf("ib4=%d whot_is_block=%d\n", ib4, t.whot_is_block[inode4]);
			//getchar();
			//}
			ib5 = whot_is_block[inode5];
			//in_model_temp(p5, ib5, b, lb);
			//if (ib5 != t.whot_is_block[inode5]) {
			//printf("ib5=%d whot_is_block=%d\n", ib5, t.whot_is_block[inode5]);
			//getchar();
			//}
			ib6 = whot_is_block[inode6];
			//in_model_temp(p6, ib6, b, lb);
			//if (ib6 != t.whot_is_block[inode6]) {
			//printf("ib6=%d whot_is_block=%d\n", ib6, t.whot_is_block[inode6]);
			//getchar();
			//}
			ib7 = whot_is_block[inode7];
			//in_model_temp(p7, ib7, b, lb);
			//if (ib7 != t.whot_is_block[inode7]) {
			//printf("ib7=%d whot_is_block=%d\n", ib7, t.whot_is_block[inode7]);
			//getchar();
			//}
			ib8 = whot_is_block[inode8];
			//in_model_temp(p8, ib8, b, lb);
			//if (ib8 != t.whot_is_block[inode8]) {
			//printf("ib8=%d whot_is_block=%d\n", ib8, t.whot_is_block[inode8]);
			//getchar();
			//}


			if ((b[ib1].bvisible) && (b[ib2].bvisible) && (b[ib3].bvisible) && (b[ib4].bvisible) && (b[ib5].bvisible) && (b[ib6].bvisible) && (b[ib7].bvisible) && (b[ib8].bvisible)) {

				bvisible_gran[6][j] = true;

				// XY bottom
				if ((ipa_count[inode1] <= bVisibleCount) &&
					(ipa_count[inode2] <= bVisibleCount) &&
					(ipa_count[inode3] <= bVisibleCount) &&
					(ipa_count[inode4] <= bVisibleCount)) {

					bvisible_gran[B_SIDE][j] = true;
				}

				// XY Top
				if ((ipa_count[inode5] <= bVisibleCount) &&
					(ipa_count[inode6] <= bVisibleCount) &&
					(ipa_count[inode7] <= bVisibleCount) &&
					(ipa_count[inode8] <= bVisibleCount)) {

					bvisible_gran[T_SIDE][j] = true;
				}

				// XZ SSIDE min Y
				if ((ipa_count[inode1] <= bVisibleCount) &&
					(ipa_count[inode2] <= bVisibleCount) &&
					(ipa_count[inode6] <= bVisibleCount) &&
					(ipa_count[inode5] <= bVisibleCount)) {

					bvisible_gran[S_SIDE][j] = true;
				}

				// XZ SSIDE max Y
				if ((ipa_count[inode4] <= bVisibleCount) &&
					(ipa_count[inode8] <= bVisibleCount) &&
					(ipa_count[inode7] <= bVisibleCount) &&
					(ipa_count[inode3] <= bVisibleCount)) {

					bvisible_gran[N_SIDE][j] = true;
				}

				// YZ SSIDE max X
				if ((ipa_count[inode2] <= bVisibleCount) &&
					(ipa_count[inode3] <= bVisibleCount) &&
					(ipa_count[inode7] <= bVisibleCount) &&
					(ipa_count[inode6] <= bVisibleCount)) {

					bvisible_gran[E_SIDE][j] = true;
				}

				// YZ SSIDE min X
				if ((ipa_count[inode1] <= bVisibleCount) &&
					(ipa_count[inode5] <= bVisibleCount) &&
					(ipa_count[inode8] <= bVisibleCount) &&
					(ipa_count[inode4] <= bVisibleCount))
				{
					bvisible_gran[W_SIDE][j] = true;
				}

			}
		}
	}

	GLFWwindow* window;

	/* Initialize the library */
	if (!glfwInit())
		return -1;

	/* Create a windowed mode window and its OpenGL context */
	window = glfwCreateWindow(SCREEN_WIDTH, SCREEN_HEIGHT, "AliceFlow_v0_49", NULL, NULL);


	glfwSetKeyCallback(window, keyCallback);
	glfwSetInputMode(window, GLFW_STICKY_KEYS, 1);

	int screenWidth, screenHeight;
	glfwGetFramebufferSize(window, &screenWidth, &screenHeight);

	if (!window)
	{
		glfwTerminate();
		return -1;
	}

	/* Make the window's context current */
	glfwMakeContextCurrent(window);


	glViewport(0.0f, 0.0f, screenWidth, screenHeight); // specifies the part of the window to which OpenGL will draw (in pixels), convert from
													   // normalised to pixels
	glMatrixMode(GL_PROJECTION); // projection matrix defines the properties of the camera that views the objects in the world coordinate frame.
								 // Here you thalfScreenHeight+scale_all*pa.yally set the zoom factor, aspect ratio and the near and far clipping planes
	glLoadIdentity(); // replace the current matrix with the identity matrix and starts us a fresh because matrix transforms such as glOrtho and
					  // glRotate cumulate, besically puts us at (0, 0, 0)
	glOrtho(0, SCREEN_WIDTH, 0, SCREEN_HEIGHT, 0.1, 10000); // essentially set coordinate system
	glMatrixMode(GL_MODELVIEW); // (default matrix mode) modelview matrix defines how objects are transformed (meaning translation, rotation
								// and scaling) in your world

	glLoadIdentity(); // same as above comment 



	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LESS);

	/* Loop until the user closes the window */
	while (!glfwWindowShouldClose(window))
	{


		if (binverse_color_black_2_white) {
			glClearColor(1.0f, 1.0f, 1.0f, 0.0f);// ����� ���
		}
		//glColor3f(0.0, 0.0, 0.0); // ׸����
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		/* Render OpenGL here */

		glPushMatrix();
		glTranslatef(halfScreenWidth, halfScreenHeight, -abbys);
		glRotatef(rotationX, 1, 0, 0);
		glRotatef(rotationY, 0, 1, 0);
		glTranslatef(-halfScreenWidth, -halfScreenHeight, abbys);

		glLineWidth(1);

		// �������� ���.
		for (int j = 0; j < maxelm; j++) {

			if (bvisible_gran[6][j]) {

				int inode1, inode2, inode3, inode4, inode5, inode6, inode7, inode8;
				inode1 = nvtx[0][j] - 1;
				inode2 = nvtx[1][j] - 1;
				inode3 = nvtx[2][j] - 1;
				inode4 = nvtx[3][j] - 1;
				inode5 = nvtx[4][j] - 1;
				inode6 = nvtx[5][j] - 1;
				inode7 = nvtx[6][j] - 1;
				inode8 = nvtx[7][j] - 1;


				// XY bottom
				if (bvisible_gran[B_SIDE][j])
				{

					if (binvisible_face_detect(0.0, 0.0, -1.0)) {


						//glColor3f(0.84,0.84,0.84);

						//glColor3f(0.0, 0.0, 0.0); // ׸����


						//Drawnvtx(nvtx, pa, j, halfScreenWidth, halfScreenHeight, -abbys, scale_all);


						glBegin(GL_QUADS);

						set_color_for_render(potent[inode4]);
						//glColor3f(1.0, 1.0, 1.0);
						//glNormal3f(0.0, 0.0, -1.0);
						glVertex3f(pa_render[inode4].x, pa_render[inode4].y, pa_render[inode4].z);
						set_color_for_render(potent[inode3]);
						//glColor3f(1.0, 1.0, 1.0);
						//glNormal3f(0.0, 0.0, -1.0);
						glVertex3f(pa_render[inode3].x, pa_render[inode3].y, pa_render[inode3].z);
						set_color_for_render(potent[inode2]);
						//glColor3f(1.0, 1.0, 1.0);
						//glNormal3f(0.0, 0.0, -1.0);
						glVertex3f(pa_render[inode2].x, pa_render[inode2].y, pa_render[inode2].z);
						set_color_for_render(potent[inode1]);
						//glColor3f(1.0, 1.0, 1.0);
						//glNormal3f(0.0, 0.0, -1.0);
						glVertex3f(pa_render[inode1].x, pa_render[inode1].y, pa_render[inode1].z);
						//glColor3f(0.0, 0.0, 0.0);
						glEnd();



						if (binverse_color_black_2_white) {
							glColor3f(0.0, 0.0, 0.0); // ׸����
						}
						else {
							glColor3f(1.0, 1.0, 1.0);// �����
						}

						if (CheckLine(ipa_count[inode1],
							ipa_count[inode2])
							)
						{

							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode1].x, pa_render[inode1].y, pa_render[inode1].z);
							glVertex3f(pa_render[inode2].x, pa_render[inode2].y, pa_render[inode2].z);
							glEnd();
						}

						if (CheckLine(ipa_count[inode3],
							ipa_count[inode2])
							)
						{

							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode2].x, pa_render[inode2].y, pa_render[inode2].z);
							glVertex3f(pa_render[inode3].x, pa_render[inode3].y, pa_render[inode3].z);
							glEnd();
						}

						if (CheckLine(ipa_count[inode3],
							ipa_count[inode4])
							)
						{
							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode3].x, pa_render[inode3].y, pa_render[inode3].z);
							glVertex3f(pa_render[inode4].x, pa_render[inode4].y, pa_render[inode4].z);
							glEnd();
						}

						if (CheckLine(ipa_count[inode1],
							ipa_count[inode4])
							)
						{
							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode4].x, pa_render[inode4].y, pa_render[inode4].z);
							glVertex3f(pa_render[inode1].x, pa_render[inode1].y, pa_render[inode1].z);
							glEnd();
						}

						//glLineWidth(1);
						//glColor3f(1.0, 1.0, 1.0);
						//glLineWidth(1);

					}
				}


				// XY Top
				if (bvisible_gran[T_SIDE][j])
				{


					if (binvisible_face_detect(0.0, 0.0, 1.0)) {


						//glColor3f(0.84,0.84,0.84);
						//glColor3f(0.0, 0.0, 0.0); // ׸����

						//Drawnvtx(nvtx, pa, j, halfScreenWidth, halfScreenHeight, -abbys, scale_all);


						glBegin(GL_QUADS);
						set_color_for_render(potent[inode8]);
						//glColor3f(1.0, 1.0, 1.0);
						//glNormal3f(0.0, 0.0, 1.0);
						glVertex3f(pa_render[inode8].x, pa_render[inode8].y, pa_render[inode8].z);
						set_color_for_render(potent[inode7]);
						//glColor3f(1.0, 1.0, 1.0);
						//glNormal3f(0.0, 0.0, 1.0);						
						glVertex3f(pa_render[inode7].x, pa_render[inode7].y, pa_render[inode7].z);
						set_color_for_render(potent[inode6]);
						//glColor3f(1.0, 1.0, 1.0);
						//glNormal3f(0.0, 0.0, 1.0);						
						glVertex3f(pa_render[inode6].x, pa_render[inode6].y, pa_render[inode6].z);
						set_color_for_render(potent[inode5]);
						//glColor3f(1.0, 1.0, 1.0);
						//glNormal3f(0.0, 0.0, 1.0);
						glVertex3f(pa_render[inode5].x, pa_render[inode5].y, pa_render[inode5].z);
						//glColor3f(0.0, 0.0, 0.0);
						glEnd();


						//glLineWidth(ComboBoxlineWidth.ItemIndex + 1);
						if (binverse_color_black_2_white) {
							glColor3f(0.0, 0.0, 0.0); // ׸����
						}
						else {
							glColor3f(1.0, 1.0, 1.0);// �����
						}

						if (CheckLine(ipa_count[inode5],
							ipa_count[inode6])
							)
						{

							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode5].x, pa_render[inode5].y, pa_render[inode5].z);
							glVertex3f(pa_render[inode6].x, pa_render[inode6].y, pa_render[inode6].z);
							glEnd();
						}

						if (CheckLine(ipa_count[inode6],
							ipa_count[inode7])
							) {

							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode6].x, pa_render[inode6].y, pa_render[inode6].z);
							glVertex3f(pa_render[inode7].x, pa_render[inode7].y, pa_render[inode7].z);
							glEnd();
						}

						if (CheckLine(ipa_count[inode7],
							ipa_count[inode8])
							)
						{
							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode7].x, pa_render[inode7].y, pa_render[inode7].z);
							glVertex3f(pa_render[inode8].x, pa_render[inode8].y, pa_render[inode8].z);
							glEnd();
						}

						if (CheckLine(ipa_count[inode5],
							ipa_count[inode8])
							)
						{
							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode8].x, pa_render[inode8].y, pa_render[inode8].z);
							glVertex3f(pa_render[inode5].x, pa_render[inode5].y, pa_render[inode5].z);
							glEnd();
						}

						// glLineWidth(1);
						//glColor3f(1.0, 1.0, 1.0);
						//glLineWidth(1);

					}
				}



				// XZ SSIDE min Y
				if (bvisible_gran[S_SIDE][j])
				{


					if (binvisible_face_detect(0.0, -1.0, 0.0)) {

						//glColor3f(0.84,0.84,0.84);


						//Drawnvtx(nvtx, pa, j, halfScreenWidth, halfScreenHeight, -abbys, scale_all);


						glBegin(GL_QUADS);
						set_color_for_render(potent[inode5]);
						//glColor3f(1.0, 1.0, 1.0);
						//glNormal3f(0.0, -1.0, 0.0);
						glVertex3f(pa_render[inode5].x, pa_render[inode5].y, pa_render[inode5].z);
						set_color_for_render(potent[inode1]);
						//glColor3f(1.0, 1.0, 1.0);
						//glNormal3f(0.0, -1.0, 0.0);
						glVertex3f(pa_render[inode1].x, pa_render[inode1].y, pa_render[inode1].z);
						set_color_for_render(potent[inode2]);
						//glColor3f(1.0, 1.0, 1.0);
						//glNormal3f(0.0, -1.0, 0.0);
						glVertex3f(pa_render[inode2].x, pa_render[inode2].y, pa_render[inode2].z);
						set_color_for_render(potent[inode6]);
						//glColor3f(1.0, 1.0, 1.0);
						///glNormal3f(0.0, -1.0, 0.0);
						glVertex3f(pa_render[inode6].x, pa_render[inode6].y, pa_render[inode6].z);
						//glColor3f(0.0, 0.0, 0.0);
						glEnd();


						//glLineWidth(ComboBoxlineWidth.ItemIndex + 1);
						if (binverse_color_black_2_white) {
							glColor3f(0.0, 0.0, 0.0); // ׸����
						}
						else {
							glColor3f(1.0, 1.0, 1.0);// �����
						}

						if (CheckLine(ipa_count[inode1],
							ipa_count[inode2])
							) {

							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode1].x, pa_render[inode1].y, pa_render[inode1].z);
							glVertex3f(pa_render[inode2].x, pa_render[inode2].y, pa_render[inode2].z);
							glEnd();
						}

						if (CheckLine(ipa_count[inode2],
							ipa_count[inode6])
							) {

							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode2].x, pa_render[inode2].y, pa_render[inode2].z);
							glVertex3f(pa_render[inode6].x, pa_render[inode6].y, pa_render[inode6].z);
							glEnd();
						}

						if (CheckLine(ipa_count[inode6],
							ipa_count[inode5])
							) {

							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode6].x, pa_render[inode6].y, pa_render[inode6].z);
							glVertex3f(pa_render[inode5].x, pa_render[inode5].y, pa_render[inode5].z);
							glEnd();
						}

						if (CheckLine(ipa_count[inode1],
							ipa_count[inode5])
							) {

							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode5].x, pa_render[inode5].y, pa_render[inode5].z);
							glVertex3f(pa_render[inode1].x, pa_render[inode1].y, pa_render[inode1].z);
							glEnd();
						}

						//glLineWidth(1);
						//glColor3f(1.0, 1.0, 1.0);
						//glLineWidth(1);
					}
				}


				// XZ SSIDE max Y
				if (bvisible_gran[N_SIDE][j])
				{

					if (binvisible_face_detect(0.0, 1.0, 0.0)) {

						//glColor3f(0.84,0.84,0.84);
						//glColor3f(0.0, 0.0, 0.0); // ׸����

						//Drawnvtx(nvtx, pa, j, halfScreenWidth, halfScreenHeight, -abbys, scale_all);


						glBegin(GL_QUADS);
						set_color_for_render(potent[inode8]);
						//glColor3f(1.0, 1.0, 1.0);
						//glNormal3f(0.0, 1.0, 0.0);
						glVertex3f(pa_render[inode8].x, pa_render[inode8].y, pa_render[inode8].z);
						set_color_for_render(potent[inode4]);
						//glColor3f(1.0, 1.0, 1.0);
						//glNormal3f(0.0, 1.0, 0.0);
						glVertex3f(pa_render[inode4].x, pa_render[inode4].y, pa_render[inode4].z);
						set_color_for_render(potent[inode3]);
						//glColor3f(1.0, 1.0, 1.0);
						//glNormal3f(0.0, 1.0, 0.0);
						glVertex3f(pa_render[inode3].x, pa_render[inode3].y, pa_render[inode3].z);
						set_color_for_render(potent[inode7]);
						//glColor3f(1.0, 1.0, 1.0);
						//glNormal3f(0.0, 1.0, 0.0);
						glVertex3f(pa_render[inode7].x, pa_render[inode7].y, pa_render[inode7].z);

						//glColor3f(0.0, 0.0, 0.0);
						glEnd();


						//glLineWidth(ComboBoxlineWidth.ItemIndex + 1);
						if (binverse_color_black_2_white) {
							glColor3f(0.0, 0.0, 0.0); // ׸����
						}
						else {
							glColor3f(1.0, 1.0, 1.0);// �����
						}

						if (CheckLine(ipa_count[inode3],
							ipa_count[inode4])
							) {

							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode4].x, pa_render[inode4].y, pa_render[inode4].z);
							glVertex3f(pa_render[inode3].x, pa_render[inode3].y, pa_render[inode3].z);
							glEnd();
						}

						if (CheckLine(ipa_count[inode3],
							ipa_count[inode7])
							) {

							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode3].x, pa_render[inode3].y, pa_render[inode3].z);
							glVertex3f(pa_render[inode7].x, pa_render[inode7].y, pa_render[inode7].z);
							glEnd();
						}

						if (CheckLine(ipa_count[inode7],
							ipa_count[inode8])
							) {
							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode7].x, pa_render[inode7].y, pa_render[inode7].z);
							glVertex3f(pa_render[inode8].x, pa_render[inode8].y, pa_render[inode8].z);
							glEnd();
						}

						if (CheckLine(ipa_count[inode8],
							ipa_count[inode4])
							) {
							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode8].x, pa_render[inode8].y, pa_render[inode8].z);
							glVertex3f(pa_render[inode4].x, pa_render[inode4].y, pa_render[inode4].z);
							glEnd();
						}

						// glLineWidth(1);
						//glColor3f(1.0, 1.0, 1.0);
						//glLineWidth(1);
					}
				}



				// YZ SSIDE max X
				if (bvisible_gran[E_SIDE][j])
				{

					if (binvisible_face_detect(1.0, 0.0, 0.0)) {


						//glColor3f(0.84,0.84,0.84);
						//glColor3f(0.0, 0.0, 0.0); // ׸����

						//Drawnvtx(nvtx, pa, j, halfScreenWidth, halfScreenHeight, -abbys, scale_all);


						glBegin(GL_QUADS);
						set_color_for_render(potent[inode6]);
						//glColor3f(1.0, 1.0, 1.0);
						//glNormal3f(1.0, 0.0, 0.0);
						glVertex3f(pa_render[inode6].x, pa_render[inode6].y, pa_render[inode6].z);
						set_color_for_render(potent[inode2]);
						//glColor3f(1.0, 1.0, 1.0);
						//glNormal3f(1.0, 0.0, 0.0);
						glVertex3f(pa_render[inode2].x, pa_render[inode2].y, pa_render[inode2].z);
						set_color_for_render(potent[inode3]);
						//glColor3f(1.0, 1.0, 1.0);
						//glNormal3f(1.0, 0.0, 0.0);
						glVertex3f(pa_render[inode3].x, pa_render[inode3].y, pa_render[inode3].z);
						set_color_for_render(potent[inode7]);
						//glColor3f(1.0, 1.0, 1.0);
						//glNormal3f(1.0, 0.0, 0.0);
						glVertex3f(pa_render[inode7].x, pa_render[inode7].y, pa_render[inode7].z);
						//glColor3f(0.0, 0.0, 0.0);
						glEnd();


						//glLineWidth(ComboBoxlineWidth.ItemIndex + 1);
						if (binverse_color_black_2_white) {
							glColor3f(0.0, 0.0, 0.0); // ׸����
						}
						else {
							glColor3f(1.0, 1.0, 1.0);// �����
						}


						if (CheckLine(ipa_count[inode3],
							ipa_count[inode2])
							)
						{
							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode2].x, pa_render[inode2].y, pa_render[inode2].z);
							glVertex3f(pa_render[inode3].x, pa_render[inode3].y, pa_render[inode3].z);
							glEnd();
						}

						if (CheckLine(ipa_count[inode3],
							ipa_count[inode7])
							)
						{
							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode3].x, pa_render[inode3].y, pa_render[inode3].z);
							glVertex3f(pa_render[inode7].x, pa_render[inode7].y, pa_render[inode7].z);
							glEnd();
						}

						if (CheckLine(ipa_count[inode7],
							ipa_count[inode6])
							)
						{
							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode7].x, pa_render[inode7].y, pa_render[inode7].z);
							glVertex3f(pa_render[inode6].x, pa_render[inode6].y, pa_render[inode6].z);
							glEnd();
						}

						if (CheckLine(ipa_count[inode6],
							ipa_count[inode2])
							)
						{
							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode6].x, pa_render[inode6].y, pa_render[inode6].z);
							glVertex3f(pa_render[inode2].x, pa_render[inode2].y, pa_render[inode2].z);
							glEnd();
						}

						//glLineWidth(1);
						//glColor3f(1.0, 1.0, 1.0);
						//glLineWidth(1);

					}
				}


				// YZ SSIDE min X
				if (bvisible_gran[W_SIDE][j])
				{

					if (binvisible_face_detect(-1.0, 0.0, 0.0))
					{

						//glColor3f(0.84,0.84,0.84);

						//glColor3f(0.0, 0.0, 0.0); // ׸����

						//Drawnvtx(nvtx, pa, j, halfScreenWidth, halfScreenHeight, -abbys, scale_all);


						glBegin(GL_QUADS);
						set_color_for_render(potent[inode8]);
						//glColor3f(1.0, 1.0, 1.0);
						//glNormal3f(-1.0, 0.0, 0.0);
						glVertex3f(pa_render[inode8].x, pa_render[inode8].y, pa_render[inode8].z);
						set_color_for_render(potent[inode4]);
						//glColor3f(1.0, 1.0, 1.0);
						//glNormal3f(-1.0, 0.0, 0.0);
						glVertex3f(pa_render[inode4].x, pa_render[inode4].y, pa_render[inode4].z);
						set_color_for_render(potent[inode1]);
						//glColor3f(1.0, 1.0, 1.0);
						//glNormal3f(-1.0, 0.0, 0.0);
						glVertex3f(pa_render[inode1].x, pa_render[inode1].y, pa_render[inode1].z);
						set_color_for_render(potent[inode5]);
						//glColor3f(1.0, 1.0, 1.0);
						//glNormal3f(-1.0, 0.0, 0.0);
						glVertex3f(pa_render[inode5].x, pa_render[inode5].y, pa_render[inode5].z);
						//glColor3f(0.0, 0.0, 0.0);
						glEnd();



						//glLineWidth(ComboBoxlineWidth.ItemIndex + 1);

						if (binverse_color_black_2_white) {
							glColor3f(0.0, 0.0, 0.0); // ׸����
						}
						else {
							glColor3f(1.0, 1.0, 1.0);// �����
						}

						if (CheckLine(ipa_count[inode1],
							ipa_count[inode4])
							) {

							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode1].x, pa_render[inode1].y, pa_render[inode1].z);
							glVertex3f(pa_render[inode4].x, pa_render[inode4].y, pa_render[inode4].z);
							glEnd();
						}

						if (CheckLine(ipa_count[inode4],
							ipa_count[inode8])
							)
						{
							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode4].x, pa_render[inode4].y, pa_render[inode4].z);
							glVertex3f(pa_render[inode8].x, pa_render[inode8].y, pa_render[inode8].z);
							glEnd();
						}

						if (CheckLine(ipa_count[inode5],
							ipa_count[inode8])
							)
						{
							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode8].x, pa_render[inode8].y, pa_render[inode8].z);
							glVertex3f(pa_render[inode5].x, pa_render[inode5].y, pa_render[inode5].z);
							glEnd();
						}

						if (CheckLine(ipa_count[inode1],
							ipa_count[inode5])
							)
						{
							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode5].x, pa_render[inode5].y, pa_render[inode5].z);
							glVertex3f(pa_render[inode1].x, pa_render[inode1].y, pa_render[inode1].z);
							glEnd();
						}

						//glLineWidth(1);
						//glColor3f(1.0, 1.0, 1.0);
						//glLineWidth(1);
					}
				}

			}
		}


		glPopMatrix();

		/* Swap front and back buffers */
		glfwSwapBuffers(window);

		/* Poll for and process events */
		glfwPollEvents();
	}

	delete[] ipa_count;

	for (int i = 0; i < 7; i++) {
		delete[] bvisible_gran[i];
	}
	delete[] bvisible_gran;

	delete[] checkline1;
	delete[] checkline2;

	glfwTerminate();

	return 0;

} // DrawZbufferColor

  // ��������� ������������ � ��������� ��������� �����. ��������� 17.12.2020.
int DrawZbuffer(BLOCK*& b, integer lb, int maxelm, int maxnod, int** &nvtx, int* &whot_is_block) {

	int CHECK_LINE_SIZE = 9;
	checkline1 = new bool[CHECK_LINE_SIZE* CHECK_LINE_SIZE];
	checkline2 = new bool[CHECK_LINE_SIZE * CHECK_LINE_SIZE];
	for (int i = 0; i < CHECK_LINE_SIZE * CHECK_LINE_SIZE; i++) {
		checkline1[i] = false;
		checkline2[i] = false;
	}

	for (int i = 0; i < CHECK_LINE_SIZE; i++) {
		for (int j = 0; j < CHECK_LINE_SIZE; j++) {
			int key = i * CHECK_LINE_SIZE + j;

			checkline1[key] = true;
			checkline2[key] = CheckLine1(i, j);
		}
	}

	bool** bvisible_gran = new bool*[7];
	for (int i = 0; i < 7; i++) {
		bvisible_gran[i] = new bool[maxelm];
	}

	int* ipa_count = new int[maxnod];
	for (int i = 0; i < maxnod; i++) {
		ipa_count[i] = 0;
	}
	for (int i = 0; i < maxelm; i++) {

		bvisible_gran[E_SIDE][i] = false;
		bvisible_gran[W_SIDE][i] = false;
		bvisible_gran[N_SIDE][i] = false;
		bvisible_gran[S_SIDE][i] = false;
		bvisible_gran[T_SIDE][i] = false;
		bvisible_gran[B_SIDE][i] = false;
		bvisible_gran[6][i] = false; // ��� ������ �������!!!


		integer inode1, inode2, inode3, inode4, inode5, inode6, inode7, inode8;
		inode1 = nvtx[0][i] - 1;
		inode2 = nvtx[1][i] - 1;
		inode3 = nvtx[2][i] - 1;
		inode4 = nvtx[3][i] - 1;
		inode5 = nvtx[4][i] - 1;
		inode6 = nvtx[5][i] - 1;
		inode7 = nvtx[6][i] - 1;
		inode8 = nvtx[7][i] - 1;

		TOCHKA ptest;
		//center_cord3D(inode1, nvtx, pa, ptest, 100);
		ptest.x = pa_opengl[inode1].x;
		ptest.y = pa_opengl[inode1].y;
		ptest.z = pa_opengl[inode1].z;
		bool bvisible = false;

		if ((fabs(pfpir.fminimum) < 1.0e-30) && (fabs(pfpir.fmaximum) < 1.0e-30)) {
			pfpir.fminimum = -1.0e30;
			pfpir.fmaximum = 1.0e30;
		}

		switch (pfpir.idir) {
		case LINE_DIRECTIONAL::X_LINE_DIRECTIONAL:
			if ((ptest.x >= pfpir.fminimum) && (ptest.x <= pfpir.fmaximum)) {
				bvisible = true;
			}
			break;
		case LINE_DIRECTIONAL::Y_LINE_DIRECTIONAL:
			if ((ptest.y >= pfpir.fminimum) && (ptest.y <= pfpir.fmaximum)) {
				bvisible = true;
			}
			break;
		case LINE_DIRECTIONAL::Z_LINE_DIRECTIONAL:
			if ((ptest.z >= pfpir.fminimum) && (ptest.z <= pfpir.fmaximum)) {
				bvisible = true;
			}
			break;
		}

		if (bvisible) {

			//TOCHKA p1, p2, p3, p4, p5, p6, p7, p8;
			//TOCHKA pall;
			//center_cord3D(inode1, nvtx, pa, p1, 100);
			//center_cord3D(inode2, nvtx, pa, p2, 100);
			//center_cord3D(inode3, nvtx, pa, p3, 100);
			//center_cord3D(inode4, nvtx, pa, p4, 100);
			//center_cord3D(inode5, nvtx, pa, p5, 100);
			//center_cord3D(inode6, nvtx, pa, p6, 100);
			//center_cord3D(inode7, nvtx, pa, p7, 100);
			//center_cord3D(inode8, nvtx, pa, p8, 100);

			integer ib1, ib2, ib3, ib4, ib5, ib6, ib7, ib8;
			// ������������� ����������������� ��� ������� whot_is_block �����������
			// ������� � �� �������� � ����� �� �� �� ���� �������� �� ������� �����.
			ib1 = whot_is_block[inode1];
			//in_model_temp(p1, ib1, b, lb);	
			//if (ib1 != t.whot_is_block[inode1]) {
			//printf("ib1=%d whot_is_block=%d\n", ib1, t.whot_is_block[inode1]);
			//getchar();
			//}
			ib2 = whot_is_block[inode2];
			//in_model_temp(p2, ib2, b, lb);
			//if (ib2 != t.whot_is_block[inode2]) {
			//printf("ib2=%d whot_is_block=%d\n", ib2, t.whot_is_block[inode2]);
			//getchar();
			//}
			ib3 = whot_is_block[inode3];
			//in_model_temp(p3, ib3, b, lb);
			//if (ib3 != t.whot_is_block[inode3]) {
			//printf("ib3=%d whot_is_block=%d\n", ib3, t.whot_is_block[inode3]);
			//getchar();
			//}
			ib4 = whot_is_block[inode4];
			//in_model_temp(p4, ib4, b, lb);
			//if (ib4 != t.whot_is_block[inode4]) {
			//printf("ib4=%d whot_is_block=%d\n", ib4, t.whot_is_block[inode4]);
			//getchar();
			//}
			ib5 = whot_is_block[inode5];
			//in_model_temp(p5, ib5, b, lb);
			//if (ib5 != t.whot_is_block[inode5]) {
			//printf("ib5=%d whot_is_block=%d\n", ib5, t.whot_is_block[inode5]);
			//getchar();
			//}
			ib6 = whot_is_block[inode6];
			//in_model_temp(p6, ib6, b, lb);
			//if (ib6 != t.whot_is_block[inode6]) {
			//printf("ib6=%d whot_is_block=%d\n", ib6, t.whot_is_block[inode6]);
			//getchar();
			//}
			ib7 = whot_is_block[inode7];
			//in_model_temp(p7, ib7, b, lb);
			//if (ib7 != t.whot_is_block[inode7]) {
			//printf("ib7=%d whot_is_block=%d\n", ib7, t.whot_is_block[inode7]);
			//getchar();
			//}
			ib8 = whot_is_block[inode8];
			//in_model_temp(p8, ib8, b, lb);
			//if (ib8 != t.whot_is_block[inode8]) {
			//printf("ib8=%d whot_is_block=%d\n", ib8, t.whot_is_block[inode8]);
			//getchar();
			//}


			if ((b[ib1].bvisible) && (b[ib2].bvisible) && (b[ib3].bvisible) && (b[ib4].bvisible) && (b[ib5].bvisible) && (b[ib6].bvisible) && (b[ib7].bvisible) && (b[ib8].bvisible)) {


				for (int j = 0; j < NUMBER_OF_VERTEX_FINITE_ELEMENT(); j++) {
					//std::cout << nvtx[j][i] - 1 << " ";
					ipa_count[nvtx[j][i] - 1]++;
				}
			}
			//getchar();
		}
	}

	const int bVisibleCount = 7;

	for (int j = 0; j < maxelm; j++) {

		integer inode1, inode2, inode3, inode4, inode5, inode6, inode7, inode8;
		inode1 = nvtx[0][j] - 1;
		inode2 = nvtx[1][j] - 1;
		inode3 = nvtx[2][j] - 1;
		inode4 = nvtx[3][j] - 1;
		inode5 = nvtx[4][j] - 1;
		inode6 = nvtx[5][j] - 1;
		inode7 = nvtx[6][j] - 1;
		inode8 = nvtx[7][j] - 1;


		TOCHKA ptest;
		//center_cord3D(inode1, nvtx, pa, ptest, 100);
		ptest.x = pa_opengl[inode1].x;
		ptest.y = pa_opengl[inode1].y;
		ptest.z = pa_opengl[inode1].z;
		bool bvisible = false;

		switch (pfpir.idir) {
		case LINE_DIRECTIONAL::X_LINE_DIRECTIONAL:
			if ((ptest.x >= pfpir.fminimum) && (ptest.x <= pfpir.fmaximum)) {
				bvisible = true;
			}
			break;
		case LINE_DIRECTIONAL::Y_LINE_DIRECTIONAL:
			if ((ptest.y >= pfpir.fminimum) && (ptest.y <= pfpir.fmaximum)) {
				bvisible = true;
			}
			break;
		case LINE_DIRECTIONAL::Z_LINE_DIRECTIONAL:
			if ((ptest.z >= pfpir.fminimum) && (ptest.z <= pfpir.fmaximum)) {
				bvisible = true;
			}
			break;
		}

		if (bvisible) {

			//TOCHKA p1, p2, p3, p4, p5, p6, p7, p8;
			//TOCHKA pall;
			//center_cord3D(inode1, nvtx, pa, p1, 100);
			//center_cord3D(inode2, nvtx, pa, p2, 100);
			//center_cord3D(inode3, nvtx, pa, p3, 100);
			//center_cord3D(inode4, nvtx, pa, p4, 100);
			//center_cord3D(inode5, nvtx, pa, p5, 100);
			//center_cord3D(inode6, nvtx, pa, p6, 100);
			//center_cord3D(inode7, nvtx, pa, p7, 100);
			//center_cord3D(inode8, nvtx, pa, p8, 100);

			integer ib1, ib2, ib3, ib4, ib5, ib6, ib7, ib8;
			// ������������� ����������������� ��� ������� whot_is_block �����������
			// ������� � �� �������� � ����� �� �� �� ���� �������� �� ������� �����.
			ib1 = whot_is_block[inode1];
			//in_model_temp(p1, ib1, b, lb);	
			//if (ib1 != t.whot_is_block[inode1]) {
			//printf("ib1=%d whot_is_block=%d\n", ib1, t.whot_is_block[inode1]);
			//getchar();
			//}
			ib2 = whot_is_block[inode2];
			//in_model_temp(p2, ib2, b, lb);
			//if (ib2 != t.whot_is_block[inode2]) {
			//printf("ib2=%d whot_is_block=%d\n", ib2, t.whot_is_block[inode2]);
			//getchar();
			//}
			ib3 = whot_is_block[inode3];
			//in_model_temp(p3, ib3, b, lb);
			//if (ib3 != t.whot_is_block[inode3]) {
			//printf("ib3=%d whot_is_block=%d\n", ib3, t.whot_is_block[inode3]);
			//getchar();
			//}
			ib4 = whot_is_block[inode4];
			//in_model_temp(p4, ib4, b, lb);
			//if (ib4 != t.whot_is_block[inode4]) {
			//printf("ib4=%d whot_is_block=%d\n", ib4, t.whot_is_block[inode4]);
			//getchar();
			//}
			ib5 = whot_is_block[inode5];
			//in_model_temp(p5, ib5, b, lb);
			//if (ib5 != t.whot_is_block[inode5]) {
			//printf("ib5=%d whot_is_block=%d\n", ib5, t.whot_is_block[inode5]);
			//getchar();
			//}
			ib6 = whot_is_block[inode6];
			//in_model_temp(p6, ib6, b, lb);
			//if (ib6 != t.whot_is_block[inode6]) {
			//printf("ib6=%d whot_is_block=%d\n", ib6, t.whot_is_block[inode6]);
			//getchar();
			//}
			ib7 = whot_is_block[inode7];
			//in_model_temp(p7, ib7, b, lb);
			//if (ib7 != t.whot_is_block[inode7]) {
			//printf("ib7=%d whot_is_block=%d\n", ib7, t.whot_is_block[inode7]);
			//getchar();
			//}
			ib8 = whot_is_block[inode8];
			//in_model_temp(p8, ib8, b, lb);
			//if (ib8 != t.whot_is_block[inode8]) {
			//printf("ib8=%d whot_is_block=%d\n", ib8, t.whot_is_block[inode8]);
			//getchar();
			//}


			if ((b[ib1].bvisible) && (b[ib2].bvisible) && (b[ib3].bvisible) && (b[ib4].bvisible) && (b[ib5].bvisible) && (b[ib6].bvisible) && (b[ib7].bvisible) && (b[ib8].bvisible)) {

				bvisible_gran[6][j] = true;

				// XY bottom
				if ((ipa_count[inode1] <= bVisibleCount) &&
					(ipa_count[inode2] <= bVisibleCount) &&
					(ipa_count[inode3] <= bVisibleCount) &&
					(ipa_count[inode4] <= bVisibleCount)) {

					bvisible_gran[B_SIDE][j] = true;
				}

				// XY Top
				if ((ipa_count[inode5] <= bVisibleCount) &&
					(ipa_count[inode6] <= bVisibleCount) &&
					(ipa_count[inode7] <= bVisibleCount) &&
					(ipa_count[inode8] <= bVisibleCount)) {

					bvisible_gran[T_SIDE][j] = true;
				}

				// XZ SSIDE min Y
				if ((ipa_count[inode1] <= bVisibleCount) &&
					(ipa_count[inode2] <= bVisibleCount) &&
					(ipa_count[inode6] <= bVisibleCount) &&
					(ipa_count[inode5] <= bVisibleCount)) {

					bvisible_gran[S_SIDE][j] = true;
				}

				// XZ SSIDE max Y
				if ((ipa_count[inode4] <= bVisibleCount) &&
					(ipa_count[inode8] <= bVisibleCount) &&
					(ipa_count[inode7] <= bVisibleCount) &&
					(ipa_count[inode3] <= bVisibleCount)) {

					bvisible_gran[N_SIDE][j] = true;
				}

				// YZ SSIDE max X
				if ((ipa_count[inode2] <= bVisibleCount) &&
					(ipa_count[inode3] <= bVisibleCount) &&
					(ipa_count[inode7] <= bVisibleCount) &&
					(ipa_count[inode6] <= bVisibleCount)) {

					bvisible_gran[E_SIDE][j] = true;
				}

				// YZ SSIDE min X
				if ((ipa_count[inode1] <= bVisibleCount) &&
					(ipa_count[inode5] <= bVisibleCount) &&
					(ipa_count[inode8] <= bVisibleCount) &&
					(ipa_count[inode4] <= bVisibleCount))
				{
					bvisible_gran[W_SIDE][j] = true;
				}

			}
		}
	}

	GLFWwindow* window;

	/* Initialize the library */
	if (!glfwInit())
		return -1;

	/* Create a windowed mode window and its OpenGL context */
	window = glfwCreateWindow(SCREEN_WIDTH, SCREEN_HEIGHT, "AliceFlow_v0_49", NULL, NULL);


	glfwSetKeyCallback(window, keyCallback);
	glfwSetInputMode(window, GLFW_STICKY_KEYS, 1);

	int screenWidth, screenHeight;
	glfwGetFramebufferSize(window, &screenWidth, &screenHeight);

	if (!window)
	{
		glfwTerminate();
		return -1;
	}

	/* Make the window's context current */
	glfwMakeContextCurrent(window);


	glViewport(0.0f, 0.0f, screenWidth, screenHeight); // specifies the part of the window to which OpenGL will draw (in pixels), convert from
													   // normalised to pixels
	glMatrixMode(GL_PROJECTION); // projection matrix defines the properties of the camera that views the objects in the world coordinate frame.
								 // Here you thalfScreenHeight+scale_all*pa.yally set the zoom factor, aspect ratio and the near and far clipping planes
	glLoadIdentity(); // replace the current matrix with the identity matrix and starts us a fresh because matrix transforms such as glOrtho and
					  // glRotate cumulate, besically puts us at (0, 0, 0)
	glOrtho(0, SCREEN_WIDTH, 0, SCREEN_HEIGHT, 0.1, 10000); // essentially set coordinate system
	glMatrixMode(GL_MODELVIEW); // (default matrix mode) modelview matrix defines how objects are transformed (meaning translation, rotation
								// and scaling) in your world

	glLoadIdentity(); // same as above comment 



	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LESS);

	/* Loop until the user closes the window */
	while (!glfwWindowShouldClose(window))
	{
		if (binverse_color_black_2_white) {
			glClearColor(1.0f, 1.0f, 1.0f, 0.0f);// ����� ���.
		}
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		/* Render OpenGL here */

		glPushMatrix();
		glTranslatef(halfScreenWidth, halfScreenHeight, -abbys);
		glRotatef(rotationX, 1, 0, 0);
		glRotatef(rotationY, 0, 1, 0);
		glTranslatef(-halfScreenWidth, -halfScreenHeight, abbys);

		glLineWidth(3);

		// �������� ���.
		for (int j = 0; j < maxelm; j++) {

			if (bvisible_gran[6][j]) {

				int inode1, inode2, inode3, inode4, inode5, inode6, inode7, inode8;
				inode1 = nvtx[0][j] - 1;
				inode2 = nvtx[1][j] - 1;
				inode3 = nvtx[2][j] - 1;
				inode4 = nvtx[3][j] - 1;
				inode5 = nvtx[4][j] - 1;
				inode6 = nvtx[5][j] - 1;
				inode7 = nvtx[6][j] - 1;
				inode8 = nvtx[7][j] - 1;


				// XY bottom
				if (bvisible_gran[B_SIDE][j])
				{

					if (binvisible_face_detect(0.0, 0.0, -1.0)) {


						//glColor3f(0.84,0.84,0.84);

						//glColor3f(0.0, 0.0, 0.0); // ׸����
						if (binverse_color_black_2_white) {
							glColor3f(1.0, 1.0, 1.0);// �����
						}
						else {
							glColor3f(0.0, 0.0, 0.0); // ׸����
						}


						//Drawnvtx(nvtx, pa, j, halfScreenWidth, halfScreenHeight, -abbys, scale_all);


						glBegin(GL_QUADS);

						//glColor3f(1.0, 1.0, 1.0);
						//glNormal3f(0.0, 0.0, -1.0);
						glVertex3f(pa_render[inode4].x, pa_render[inode4].y, pa_render[inode4].z);
						//glColor3f(1.0, 1.0, 1.0);
						//glNormal3f(0.0, 0.0, -1.0);
						glVertex3f(pa_render[inode3].x, pa_render[inode3].y, pa_render[inode3].z);
						//glColor3f(1.0, 1.0, 1.0);
						//glNormal3f(0.0, 0.0, -1.0);
						glVertex3f(pa_render[inode2].x, pa_render[inode2].y, pa_render[inode2].z);
						//glColor3f(1.0, 1.0, 1.0);
						//glNormal3f(0.0, 0.0, -1.0);
						glVertex3f(pa_render[inode1].x, pa_render[inode1].y, pa_render[inode1].z);
						//glColor3f(0.0, 0.0, 0.0);
						glEnd();



						if (binverse_color_black_2_white) {
							glColor3f(0.0, 0.0, 0.0); // ׸����
						}
						else {
							glColor3f(1.0, 1.0, 1.0);// �����
						}

						if (CheckLine(ipa_count[inode1],
							ipa_count[inode2])
							)
						{

							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode1].x, pa_render[inode1].y, pa_render[inode1].z);
							glVertex3f(pa_render[inode2].x, pa_render[inode2].y, pa_render[inode2].z);
							glEnd();
						}

						if (CheckLine(ipa_count[inode3],
							ipa_count[inode2])
							)
						{

							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode2].x, pa_render[inode2].y, pa_render[inode2].z);
							glVertex3f(pa_render[inode3].x, pa_render[inode3].y, pa_render[inode3].z);
							glEnd();
						}

						if (CheckLine(ipa_count[inode3],
							ipa_count[inode4])
							)
						{
							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode3].x, pa_render[inode3].y, pa_render[inode3].z);
							glVertex3f(pa_render[inode4].x, pa_render[inode4].y, pa_render[inode4].z);
							glEnd();
						}

						if (CheckLine(ipa_count[inode1],
							ipa_count[inode4])
							)
						{
							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode4].x, pa_render[inode4].y, pa_render[inode4].z);
							glVertex3f(pa_render[inode1].x, pa_render[inode1].y, pa_render[inode1].z);
							glEnd();
						}

						//glLineWidth(1);
						//glColor3f(1.0, 1.0, 1.0);
						//glLineWidth(1);

					}
				}


				// XY Top
				if (bvisible_gran[T_SIDE][j])
				{


					if (binvisible_face_detect(0.0, 0.0, 1.0)) {


						//glColor3f(0.84,0.84,0.84);
						if (binverse_color_black_2_white) {
							glColor3f(1.0, 1.0, 1.0);// �����
						}
						else {
							glColor3f(0.0, 0.0, 0.0); // ׸����
						}

						//Drawnvtx(nvtx, pa, j, halfScreenWidth, halfScreenHeight, -abbys, scale_all);


						glBegin(GL_QUADS);
						//glColor3f(1.0, 1.0, 1.0);
						//glNormal3f(0.0, 0.0, 1.0);
						glVertex3f(pa_render[inode8].x, pa_render[inode8].y, pa_render[inode8].z);
						//glColor3f(1.0, 1.0, 1.0);
						//glNormal3f(0.0, 0.0, 1.0);
						glVertex3f(pa_render[inode7].x, pa_render[inode7].y, pa_render[inode7].z);
						//glColor3f(1.0, 1.0, 1.0);
						//glNormal3f(0.0, 0.0, 1.0);
						glVertex3f(pa_render[inode6].x, pa_render[inode6].y, pa_render[inode6].z);
						//glColor3f(1.0, 1.0, 1.0);
						//glNormal3f(0.0, 0.0, 1.0);
						glVertex3f(pa_render[inode5].x, pa_render[inode5].y, pa_render[inode5].z);
						//glColor3f(0.0, 0.0, 0.0);
						glEnd();


						//glLineWidth(ComboBoxlineWidth.ItemIndex + 1);
						if (binverse_color_black_2_white) {
							glColor3f(0.0, 0.0, 0.0); // ׸����
						}
						else {
							glColor3f(1.0, 1.0, 1.0);// �����
						}

						if (CheckLine(ipa_count[inode5],
							ipa_count[inode6])
							)
						{

							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode5].x, pa_render[inode5].y, pa_render[inode5].z);
							glVertex3f(pa_render[inode6].x, pa_render[inode6].y, pa_render[inode6].z);
							glEnd();
						}

						if (CheckLine(ipa_count[inode6],
							ipa_count[inode7])
							) {

							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode6].x, pa_render[inode6].y, pa_render[inode6].z);
							glVertex3f(pa_render[inode7].x, pa_render[inode7].y, pa_render[inode7].z);
							glEnd();
						}

						if (CheckLine(ipa_count[inode7],
							ipa_count[inode8])
							)
						{
							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode7].x, pa_render[inode7].y, pa_render[inode7].z);
							glVertex3f(pa_render[inode8].x, pa_render[inode8].y, pa_render[inode8].z);
							glEnd();
						}

						if (CheckLine(ipa_count[inode5],
							ipa_count[inode8])
							)
						{
							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode8].x, pa_render[inode8].y, pa_render[inode8].z);
							glVertex3f(pa_render[inode5].x, pa_render[inode5].y, pa_render[inode5].z);
							glEnd();
						}

						// glLineWidth(1);
						//glColor3f(1.0, 1.0, 1.0);
						//glLineWidth(1);

					}
				}



				// XZ SSIDE min Y
				if (bvisible_gran[S_SIDE][j])
				{


					if (binvisible_face_detect(0.0, -1.0, 0.0)) {

						//glColor3f(0.84,0.84,0.84);
						if (binverse_color_black_2_white) {
							glColor3f(1.0, 1.0, 1.0);// �����
						}
						else {
							glColor3f(0.0, 0.0, 0.0); // ׸����
						}

						//Drawnvtx(nvtx, pa, j, halfScreenWidth, halfScreenHeight, -abbys, scale_all);


						glBegin(GL_QUADS);
						//glColor3f(1.0, 1.0, 1.0);
						//glNormal3f(0.0, -1.0, 0.0);
						glVertex3f(pa_render[inode5].x, pa_render[inode5].y, pa_render[inode5].z);
						//glColor3f(1.0, 1.0, 1.0);
						//glNormal3f(0.0, -1.0, 0.0);
						glVertex3f(pa_render[inode1].x, pa_render[inode1].y, pa_render[inode1].z);
						//glColor3f(1.0, 1.0, 1.0);
						//glNormal3f(0.0, -1.0, 0.0);
						glVertex3f(pa_render[inode2].x, pa_render[inode2].y, pa_render[inode2].z);
						//glColor3f(1.0, 1.0, 1.0);
						///glNormal3f(0.0, -1.0, 0.0);
						glVertex3f(pa_render[inode6].x, pa_render[inode6].y, pa_render[inode6].z);
						//glColor3f(0.0, 0.0, 0.0);
						glEnd();


						//glLineWidth(ComboBoxlineWidth.ItemIndex + 1);
						if (binverse_color_black_2_white) {
							glColor3f(0.0, 0.0, 0.0); // ׸����
						}
						else {
							glColor3f(1.0, 1.0, 1.0);// �����
						}

						if (CheckLine(ipa_count[inode1],
							ipa_count[inode2])
							) {

							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode1].x, pa_render[inode1].y, pa_render[inode1].z);
							glVertex3f(pa_render[inode2].x, pa_render[inode2].y, pa_render[inode2].z);
							glEnd();
						}

						if (CheckLine(ipa_count[inode2],
							ipa_count[inode6])
							) {

							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode2].x, pa_render[inode2].y, pa_render[inode2].z);
							glVertex3f(pa_render[inode6].x, pa_render[inode6].y, pa_render[inode6].z);
							glEnd();
						}

						if (CheckLine(ipa_count[inode6],
							ipa_count[inode5])
							) {

							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode6].x, pa_render[inode6].y, pa_render[inode6].z);
							glVertex3f(pa_render[inode5].x, pa_render[inode5].y, pa_render[inode5].z);
							glEnd();
						}

						if (CheckLine(ipa_count[inode1],
							ipa_count[inode5])
							) {

							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode5].x, pa_render[inode5].y, pa_render[inode5].z);
							glVertex3f(pa_render[inode1].x, pa_render[inode1].y, pa_render[inode1].z);
							glEnd();
						}

						//glLineWidth(1);
						//glColor3f(1.0, 1.0, 1.0);
						//glLineWidth(1);
					}
				}


				// XZ SSIDE max Y
				if (bvisible_gran[N_SIDE][j])
				{

					if (binvisible_face_detect(0.0, 1.0, 0.0)) {

						//glColor3f(0.84,0.84,0.84);
						if (binverse_color_black_2_white) {
							glColor3f(1.0, 1.0, 1.0);// �����
						}
						else {
							glColor3f(0.0, 0.0, 0.0); // ׸����
						}

						//Drawnvtx(nvtx, pa, j, halfScreenWidth, halfScreenHeight, -abbys, scale_all);


						glBegin(GL_QUADS);

						//glColor3f(1.0, 1.0, 1.0);
						//glNormal3f(0.0, 1.0, 0.0);
						glVertex3f(pa_render[inode8].x, pa_render[inode8].y, pa_render[inode8].z);
						//glColor3f(1.0, 1.0, 1.0);
						//glNormal3f(0.0, 1.0, 0.0);
						glVertex3f(pa_render[inode4].x, pa_render[inode4].y, pa_render[inode4].z);
						//glColor3f(1.0, 1.0, 1.0);
						//glNormal3f(0.0, 1.0, 0.0);
						glVertex3f(pa_render[inode3].x, pa_render[inode3].y, pa_render[inode3].z);
						//glColor3f(1.0, 1.0, 1.0);
						//glNormal3f(0.0, 1.0, 0.0);
						glVertex3f(pa_render[inode7].x, pa_render[inode7].y, pa_render[inode7].z);

						//glColor3f(0.0, 0.0, 0.0);
						glEnd();


						//glLineWidth(ComboBoxlineWidth.ItemIndex + 1);
						if (binverse_color_black_2_white) {
							glColor3f(0.0, 0.0, 0.0); // ׸����
						}
						else {
							glColor3f(1.0, 1.0, 1.0);// �����
						}

						if (CheckLine(ipa_count[inode3],
							ipa_count[inode4])
							) {

							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode4].x, pa_render[inode4].y, pa_render[inode4].z);
							glVertex3f(pa_render[inode3].x, pa_render[inode3].y, pa_render[inode3].z);
							glEnd();
						}

						if (CheckLine(ipa_count[inode3],
							ipa_count[inode7])
							) {

							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode3].x, pa_render[inode3].y, pa_render[inode3].z);
							glVertex3f(pa_render[inode7].x, pa_render[inode7].y, pa_render[inode7].z);
							glEnd();
						}

						if (CheckLine(ipa_count[inode7],
							ipa_count[inode8])
							) {
							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode7].x, pa_render[inode7].y, pa_render[inode7].z);
							glVertex3f(pa_render[inode8].x, pa_render[inode8].y, pa_render[inode8].z);
							glEnd();
						}

						if (CheckLine(ipa_count[inode8],
							ipa_count[inode4])
							) {
							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode8].x, pa_render[inode8].y, pa_render[inode8].z);
							glVertex3f(pa_render[inode4].x, pa_render[inode4].y, pa_render[inode4].z);
							glEnd();
						}

						// glLineWidth(1);
						//glColor3f(1.0, 1.0, 1.0);
						//glLineWidth(1);
					}
				}



				// YZ SSIDE max X
				if (bvisible_gran[E_SIDE][j])
				{

					if (binvisible_face_detect(1.0, 0.0, 0.0)) {


						//glColor3f(0.84,0.84,0.84);
						if (binverse_color_black_2_white) {
							glColor3f(1.0, 1.0, 1.0);// �����
						}
						else {
							glColor3f(0.0, 0.0, 0.0); // ׸����
						}

						//Drawnvtx(nvtx, pa, j, halfScreenWidth, halfScreenHeight, -abbys, scale_all);


						glBegin(GL_QUADS);
						//glColor3f(1.0, 1.0, 1.0);
						//glNormal3f(1.0, 0.0, 0.0);
						glVertex3f(pa_render[inode6].x, pa_render[inode6].y, pa_render[inode6].z);
						//glColor3f(1.0, 1.0, 1.0);
						//glNormal3f(1.0, 0.0, 0.0);
						glVertex3f(pa_render[inode2].x, pa_render[inode2].y, pa_render[inode2].z);
						//glColor3f(1.0, 1.0, 1.0);
						//glNormal3f(1.0, 0.0, 0.0);
						glVertex3f(pa_render[inode3].x, pa_render[inode3].y, pa_render[inode3].z);
						//glColor3f(1.0, 1.0, 1.0);
						//glNormal3f(1.0, 0.0, 0.0);
						glVertex3f(pa_render[inode7].x, pa_render[inode7].y, pa_render[inode7].z);
						//glColor3f(0.0, 0.0, 0.0);
						glEnd();


						//glLineWidth(ComboBoxlineWidth.ItemIndex + 1);
						if (binverse_color_black_2_white) {
							glColor3f(0.0, 0.0, 0.0); // ׸����
						}
						else {
							glColor3f(1.0, 1.0, 1.0);// �����
						}



						if (CheckLine(ipa_count[inode3],
							ipa_count[inode2])
							)
						{
							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode2].x, pa_render[inode2].y, pa_render[inode2].z);
							glVertex3f(pa_render[inode3].x, pa_render[inode3].y, pa_render[inode3].z);
							glEnd();
						}

						if (CheckLine(ipa_count[inode3],
							ipa_count[inode7])
							)
						{
							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode3].x, pa_render[inode3].y, pa_render[inode3].z);
							glVertex3f(pa_render[inode7].x, pa_render[inode7].y, pa_render[inode7].z);
							glEnd();
						}

						if (CheckLine(ipa_count[inode7],
							ipa_count[inode6])
							)
						{
							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode7].x, pa_render[inode7].y, pa_render[inode7].z);
							glVertex3f(pa_render[inode6].x, pa_render[inode6].y, pa_render[inode6].z);
							glEnd();
						}

						if (CheckLine(ipa_count[inode6],
							ipa_count[inode2])
							)
						{
							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode6].x, pa_render[inode6].y, pa_render[inode6].z);
							glVertex3f(pa_render[inode2].x, pa_render[inode2].y, pa_render[inode2].z);
							glEnd();
						}

						//glLineWidth(1);
						//glColor3f(1.0, 1.0, 1.0);
						//glLineWidth(1);

					}
				}


				// YZ SSIDE min X
				if (bvisible_gran[W_SIDE][j])
				{

					if (binvisible_face_detect(-1.0, 0.0, 0.0))
					{

						//glColor3f(0.84,0.84,0.84);

						if (binverse_color_black_2_white) {
							glColor3f(1.0, 1.0, 1.0);// �����
						}
						else {
							glColor3f(0.0, 0.0, 0.0); // ׸����
						}

						//Drawnvtx(nvtx, pa, j, halfScreenWidth, halfScreenHeight, -abbys, scale_all);


						glBegin(GL_QUADS);
						//glColor3f(1.0, 1.0, 1.0);
						//glNormal3f(-1.0, 0.0, 0.0);
						glVertex3f(pa_render[inode8].x, pa_render[inode8].y, pa_render[inode8].z);
						//glColor3f(1.0, 1.0, 1.0);
						//glNormal3f(-1.0, 0.0, 0.0);
						glVertex3f(pa_render[inode4].x, pa_render[inode4].y, pa_render[inode4].z);
						//glColor3f(1.0, 1.0, 1.0);
						//glNormal3f(-1.0, 0.0, 0.0);
						glVertex3f(pa_render[inode1].x, pa_render[inode1].y, pa_render[inode1].z);
						//glColor3f(1.0, 1.0, 1.0);
						//glNormal3f(-1.0, 0.0, 0.0);
						glVertex3f(pa_render[inode5].x, pa_render[inode5].y, pa_render[inode5].z);
						//glColor3f(0.0, 0.0, 0.0);
						glEnd();



						//glLineWidth(ComboBoxlineWidth.ItemIndex + 1);

						if (binverse_color_black_2_white) {
							glColor3f(0.0, 0.0, 0.0); // ׸����
						}
						else {
							glColor3f(1.0, 1.0, 1.0);// �����
						}

						if (CheckLine(ipa_count[inode1],
							ipa_count[inode4])
							) {

							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode1].x, pa_render[inode1].y, pa_render[inode1].z);
							glVertex3f(pa_render[inode4].x, pa_render[inode4].y, pa_render[inode4].z);
							glEnd();
						}

						if (CheckLine(ipa_count[inode4],
							ipa_count[inode8])
							)
						{
							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode4].x, pa_render[inode4].y, pa_render[inode4].z);
							glVertex3f(pa_render[inode8].x, pa_render[inode8].y, pa_render[inode8].z);
							glEnd();
						}

						if (CheckLine(ipa_count[inode5],
							ipa_count[inode8])
							)
						{
							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode8].x, pa_render[inode8].y, pa_render[inode8].z);
							glVertex3f(pa_render[inode5].x, pa_render[inode5].y, pa_render[inode5].z);
							glEnd();
						}

						if (CheckLine(ipa_count[inode1],
							ipa_count[inode5])
							)
						{
							glBegin(GL_LINE_LOOP);
							glVertex3f(pa_render[inode5].x, pa_render[inode5].y, pa_render[inode5].z);
							glVertex3f(pa_render[inode1].x, pa_render[inode1].y, pa_render[inode1].z);
							glEnd();
						}

						//glLineWidth(1);
						//glColor3f(1.0, 1.0, 1.0);
						//glLineWidth(1);
					}
				}

			}
		}


		glPopMatrix();

		/* Swap front and back buffers */
		glfwSwapBuffers(window);

		/* Poll for and process events */
		glfwPollEvents();
	}

	delete[] ipa_count;

	for (int i = 0; i < 7; i++) {
		delete[] bvisible_gran[i];
	}
	delete[] bvisible_gran;

	delete[] checkline1;
	delete[] checkline2;

	glfwTerminate();

	return 0;

} // DrawZbuffer

#endif

  // ������� ����������� ����� ���������� �������.
int main_copy1(void)
{

	// ���������� ������, ���������� � ������, �������.
	integer lb = 0;
	BLOCK* b = nullptr;// ������ ������

					   // �������� ������.
	main_solverv0_48(b, lb);

#ifndef NO_OPENGL_GLFW

	GLFWwindow* window;

	/* Initialize the library */
	if (!glfwInit())
		return -1;

	/* Create a windowed mode window and its OpenGL context */
	window = glfwCreateWindow(SCREEN_WIDTH, SCREEN_HEIGHT, "AliceFlow_v0_49", NULL, NULL);


	glfwSetKeyCallback(window, keyCallback);
	glfwSetInputMode(window, GLFW_STICKY_KEYS, 1);

	int screenWidth, screenHeight;
	glfwGetFramebufferSize(window, &screenWidth, &screenHeight);

	if (!window)
	{
		glfwTerminate();
		return -1;
	}

	/* Make the window's context current */
	glfwMakeContextCurrent(window);


	glViewport(0.0f, 0.0f, screenWidth, screenHeight); // specifies the part of the window to which OpenGL will draw (in pixels), convert from
													   // normalised to pixels
	glMatrixMode(GL_PROJECTION); // projection matrix defines the properties of the camera that views the objects in the world coordinate frame.
								 // Here you thalfScreenHeight+scale_all*pa.yally set the zoom factor, aspect ratio and the near and far clipping planes
	glLoadIdentity(); // replace the current matrix with the identity matrix and starts us a fresh because matrix transforms such as glOrtho and
					  // glRotate cumulate, besically puts us at (0, 0, 0)
	glOrtho(0, SCREEN_WIDTH, 0, SCREEN_HEIGHT, 0, 1000); // essentially set coordinate system
	glMatrixMode(GL_MODELVIEW); // (default matrix mode) modelview matrix defines how objects are transformed (meaning translation, rotation
								// and scaling) in your world

	glLoadIdentity(); // same as above comment 

					  //GLfloat halfScreenWidth = SCREEN_WIDTH / 2;
					  //GLfloat halfScreenHeight = SCREEN_HEIGHT / 2;



					  /* Loop until the user closes the window */
	while (!glfwWindowShouldClose(window))
	{

		//glClearColor(0.7f, 1.0f, 0.7f, 0.0f);
		glClear(GL_COLOR_BUFFER_BIT);

		/* Render OpenGL here */

		glPushMatrix();
		glTranslatef(halfScreenWidth, halfScreenHeight, -500);
		glRotatef(rotationX, 1, 0, 0);
		glRotatef(rotationY, 0, 1, 0);
		glTranslatef(-halfScreenWidth, -halfScreenHeight, 500);

		DrawKarkas(b, lb, halfScreenWidth, halfScreenHeight, -500.0, scale_all);

		glPopMatrix();

		/* Swap front and back buffers */
		glfwSwapBuffers(window);

		/* Poll for and process events */
		glfwPollEvents();
	}

	glfwTerminate();

#endif


	for (integer i_7 = 0; i_7 < lb; i_7++) {
		if (b[i_7].temp_Sc != nullptr) {
			delete[] b[i_7].temp_Sc;
			b[i_7].temp_Sc = nullptr;
		}
		if (b[i_7].arr_Sc != nullptr) {
			delete[] b[i_7].arr_Sc;
			b[i_7].arr_Sc = nullptr;
		}
		if (b[i_7].g.hi != nullptr) {
			delete[] b[i_7].g.hi;
			b[i_7].g.hi = nullptr;
		}
		if (b[i_7].g.xi != nullptr) {
			delete[] b[i_7].g.xi;
			b[i_7].g.xi = nullptr;
		}
		if (b[i_7].g.yi != nullptr) {
			delete[] b[i_7].g.yi;
			b[i_7].g.yi = nullptr;
		}
		if (b[i_7].g.zi != nullptr) {
			delete[] b[i_7].g.zi;
			b[i_7].g.zi = nullptr;
		}
	}
	delete[] b;
	b = nullptr;

	return 0;
} // ���������� ����� ������� �1

#ifndef NO_OPENGL_GLFW

void keyCallback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
	// std::cout << key << std::endl;

	const GLfloat rotationSpeed = 10;
	doublereal dmin = 1.0e30;
	doublereal dmax = -1.0e30;

	// actions are GLFW_PRESS, GLFW_RELEASE or GLFW_REPEAT
	if (action == GLFW_PRESS || action == GLFW_REPEAT)
	{
		switch (key)
		{
		case GLFW_KEY_UP:
			rotationX -= rotationSpeed;
			rotationX = ((int)rotationX) % 360;
			break;
		case GLFW_KEY_DOWN:
			rotationX += rotationSpeed;
			rotationX = ((int)rotationX) % 360;
			break;
		case GLFW_KEY_RIGHT:
			rotationY += rotationSpeed;
			rotationY = ((int)rotationY) % 360;
			break;
		case GLFW_KEY_LEFT:
			rotationY -= rotationSpeed;
			rotationY = ((int)rotationY) % 360;
			break;
		case GLFW_KEY_R:
			// ��������� ����.
			if (animation_sequence_functions_openGL != nullptr)
			{
				if (iCURENT_ANIMATION_CADER < iNUMBER_ANIMATION_CADERS - 1) {
					iCURENT_ANIMATION_CADER++;
				}
				else {
					iCURENT_ANIMATION_CADER = 0;
				}

				dmin = 1.0e30;
				dmax = -1.0e30;


				for (int i37 = 0; i37 < n_render; i37++) {
					if (animation_sequence_functions_openGL[iCURENT_ANIMATION_FUNCTION][iCURENT_ANIMATION_CADER][i37] > dmax) {
						dmax = animation_sequence_functions_openGL[iCURENT_ANIMATION_FUNCTION][iCURENT_ANIMATION_CADER][i37];
					}
					if (animation_sequence_functions_openGL[iCURENT_ANIMATION_FUNCTION][iCURENT_ANIMATION_CADER][i37] < dmin) {
						dmin = animation_sequence_functions_openGL[iCURENT_ANIMATION_FUNCTION][iCURENT_ANIMATION_CADER][i37];
					}
				}

				minimum_val_for_render_pic = dmin;
				maximum_val_for_render_pic = dmax;
			}
			break;
		case GLFW_KEY_F:
			// ���������� ����.
			if (animation_sequence_functions_openGL != nullptr)
			{
				if (iCURENT_ANIMATION_CADER > 0) {
					iCURENT_ANIMATION_CADER--;
				}
				else {
					iCURENT_ANIMATION_CADER = iNUMBER_ANIMATION_CADERS - 1;
				}

				dmin = 1.0e30;
				dmax = -1.0e30;


				for (int i37 = 0; i37 < n_render; i37++) {
					if (animation_sequence_functions_openGL[iCURENT_ANIMATION_FUNCTION][iCURENT_ANIMATION_CADER][i37] > dmax) {
						dmax = animation_sequence_functions_openGL[iCURENT_ANIMATION_FUNCTION][iCURENT_ANIMATION_CADER][i37];
					}
					if (animation_sequence_functions_openGL[iCURENT_ANIMATION_FUNCTION][iCURENT_ANIMATION_CADER][i37] < dmin) {
						dmin = animation_sequence_functions_openGL[iCURENT_ANIMATION_FUNCTION][iCURENT_ANIMATION_CADER][i37];
					}
				}

				minimum_val_for_render_pic = dmin;
				maximum_val_for_render_pic = dmax;
			}
			break;
			/*
			// ��� ��� �� ��������. ��� ��������� ��� ������.
			case GLFW_KEY_O:
			halfScreenHeight += 0.05 * SCREEN_HEIGHT;
			break;
			case GLFW_KEY_L:
			halfScreenHeight -= 0.05 * SCREEN_HEIGHT;
			break;
			case GLFW_KEY_U:
			halfScreenWidth += 0.05 * SCREEN_WIDTH;
			break;
			case GLFW_KEY_Y:
			halfScreenWidth -= 0.05 * SCREEN_WIDTH;
			break;
			*/
		case GLFW_KEY_S:
			scale_all = scale_all / 1.2f;
			for (int i = 0; i < n_render; i++) {
				pa_render[i].x = halfScreenWidth + scale_all * pa_opengl[i].x;
				pa_render[i].y = halfScreenHeight + scale_all * pa_opengl[i].y;
				pa_render[i].z = -abbys + scale_all * pa_opengl[i].z;
			}
			break;
		case GLFW_KEY_W:
			scale_all = scale_all * 1.2f;
			for (int i = 0; i < n_render; i++) {
				pa_render[i].x = halfScreenWidth + scale_all * pa_opengl[i].x;
				pa_render[i].y = halfScreenHeight + scale_all * pa_opengl[i].y;
				pa_render[i].z = -abbys + scale_all * pa_opengl[i].z;
			}
			break;
		case GLFW_KEY_E:

			if (iCFD_animation == 1) {

				if (animation_sequence_functions_openGL != nullptr)
				{

					if (iCURENT_ANIMATION_FUNCTION < iNUMBER_ANIMATION_FUNCTIONS - 1) {
						iCURENT_ANIMATION_FUNCTION++;
					}
					else {
						iCURENT_ANIMATION_FUNCTION = 0;
					}

					dmin = 1.0e30;
					dmax = -1.0e30;


					for (int i37 = 0; i37 < n_render; i37++) {
						if (animation_sequence_functions_openGL[iCURENT_ANIMATION_FUNCTION][iCURENT_ANIMATION_CADER][i37] > dmax) {
							dmax = animation_sequence_functions_openGL[iCURENT_ANIMATION_FUNCTION][iCURENT_ANIMATION_CADER][i37];
						}
						if (animation_sequence_functions_openGL[iCURENT_ANIMATION_FUNCTION][iCURENT_ANIMATION_CADER][i37] < dmin) {
							dmin = animation_sequence_functions_openGL[iCURENT_ANIMATION_FUNCTION][iCURENT_ANIMATION_CADER][i37];
						}
					}

					minimum_val_for_render_pic = dmin;
					maximum_val_for_render_pic = dmax;
				}

			}
			else {
				if (potent_array_opengl != nullptr) {

					if (iCURENT_FUNC_openGL < iNUMBER_FUNC_openGL - 1) {
						iCURENT_FUNC_openGL++;
					}
					else {
						iCURENT_FUNC_openGL = 0;
					}
					dmin = 1.0e30;
					dmax = -1.0e30;


					for (int i37 = 0; i37 < n_render; i37++) {
						if (potent_array_opengl[iCURENT_FUNC_openGL][i37] > dmax) {
							dmax = potent_array_opengl[iCURENT_FUNC_openGL][i37];
						}
						if (potent_array_opengl[iCURENT_FUNC_openGL][i37] < dmin) {
							dmin = potent_array_opengl[iCURENT_FUNC_openGL][i37];
						}
					}

					minimum_val_for_render_pic = dmin;
					maximum_val_for_render_pic = dmax;
				}
			}
			break;

		}
		GLfloat Alf = 3.14159265 * rotationX / 180.0;// ���� � ��������.
		GLfloat Bet = 3.14159265 * rotationY / 180.0;// ���� � ��������.
		cosAlf = cos(Alf);
		sinAlf = sin(Alf);
		cosBet = cos(Bet);
		sinBet = sin(Bet);

	}

}


void DrawCADobj(BLOCK*& b, integer lb, integer id, GLfloat centerPosX, GLfloat centerPosY, GLfloat  centerPosZ, GLfloat scale) {
	//if (!b[id].g.bbigCADmodel) 
	{

		bool bvisible = true;

		if (id == 0) {
			bvisible = false;
		}
		if (b[id].itype == PHYSICS_TYPE_IN_BODY::HOLLOW) {
			bvisible = false;
		}


		if (bvisible) {

			glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

			gCAD_STL* tmp = b[id].g.root_CAD_STL;

			while (tmp != nullptr) {

				glBegin(GL_TRIANGLES); // ��������� ��� ��������
				glNormal3f(tmp->n.x, tmp->n.y, tmp->n.z);
				glVertex3f(scale * tmp->pa.x + centerPosX, scale * tmp->pa.y + centerPosY, scale * tmp->pa.z + centerPosZ);
				glVertex3f(scale * tmp->pb.x + centerPosX, scale * tmp->pb.y + centerPosY, scale * tmp->pb.z + centerPosZ);
				glVertex3f(scale * tmp->pc.x + centerPosX, scale * tmp->pc.y + centerPosY, scale * tmp->pc.z + centerPosZ);
				glEnd(); // ��������� ��������

				tmp = tmp->next;
			}
		}
	}
}

void DrawPrism(BLOCK*& b, integer lb, integer id, GLfloat centerPosX, GLfloat centerPosY, GLfloat  centerPosZ, GLfloat scale) {

	bool bCAD = false;
	for (int i = 0; i < lb; i++) {
		if (b[i].g.itypegeom == CAD_STL) {
			bCAD = true;
		}
	}

	bool bvisible = true;
	if (bCAD) {
		if (id == 0) {
			bvisible = false;
		}
		if (b[id].itype == PHYSICS_TYPE_IN_BODY::HOLLOW) {
			bvisible = false;
		}
	}

	if (bvisible) {

		GLfloat vertices[] =
		{
			// front face
			scale * b[id].g.xS + centerPosX, scale * b[id].g.yE + centerPosY, scale * b[id].g.zE + centerPosZ, // top left
			scale * b[id].g.xE + centerPosX, scale * b[id].g.yE + centerPosY, scale * b[id].g.zE + centerPosZ, // top right
			scale * b[id].g.xE + centerPosX, scale * b[id].g.yS + centerPosY, scale * b[id].g.zE + centerPosZ, // bottom right
			scale * b[id].g.xS + centerPosX, scale * b[id].g.yS + centerPosY, scale * b[id].g.zE + centerPosZ, // bottom left

																											   // back face
																											   scale * b[id].g.xS + centerPosX, scale * b[id].g.yE + centerPosY, scale * b[id].g.zS + centerPosZ, // top left
																											   scale * b[id].g.xE + centerPosX, scale * b[id].g.yE + centerPosY,  scale * b[id].g.zS + centerPosZ, // top right
																											   scale * b[id].g.xE + centerPosX, scale * b[id].g.yS + centerPosY, scale * b[id].g.zS + centerPosZ, // bottom right
																											   scale * b[id].g.xS + centerPosX, scale * b[id].g.yS + centerPosY, scale * b[id].g.zS + centerPosZ, // bottom left

																																																				  // left face
																																																				  scale * b[id].g.xS + centerPosX,  scale * b[id].g.yE + centerPosY,  scale * b[id].g.zE + centerPosZ, // top left
																																																				  scale * b[id].g.xS + centerPosX, scale * b[id].g.yE + centerPosY,  scale * b[id].g.zS + centerPosZ, // top right
																																																				  scale * b[id].g.xS + centerPosX,  scale * b[id].g.yS + centerPosY, scale * b[id].g.zS + centerPosZ, // bottom right
																																																				  scale * b[id].g.xS + centerPosX, scale * b[id].g.yS + centerPosY, scale * b[id].g.zE + centerPosZ, // bottom left

																																																																													 // right face
																																																																													 scale * b[id].g.xE + centerPosX,  scale * b[id].g.yE + centerPosY, scale * b[id].g.zE + centerPosZ, // top left
																																																																													 scale * b[id].g.xE + centerPosX,  scale * b[id].g.yE + centerPosY,  scale * b[id].g.zS + centerPosZ, // top right
																																																																													 scale * b[id].g.xE + centerPosX, scale * b[id].g.yS + centerPosY,  scale * b[id].g.zS + centerPosZ, // bottom right
																																																																													 scale * b[id].g.xE + centerPosX,  scale * b[id].g.yS + centerPosY,  scale * b[id].g.zE + centerPosZ, // bottom left

																																																																																																						  // top face
																																																																																																						  scale * b[id].g.xS + centerPosX,  scale * b[id].g.yE + centerPosY,  scale * b[id].g.zE + centerPosZ, // top left
																																																																																																						  scale * b[id].g.xS + centerPosX,  scale * b[id].g.yE + centerPosY,  scale * b[id].g.zS + centerPosZ, // top right
																																																																																																						  scale * b[id].g.xE + centerPosX,  scale * b[id].g.yE + centerPosY,  scale * b[id].g.zS + centerPosZ, // bottom right
																																																																																																						  scale * b[id].g.xE + centerPosX,  scale * b[id].g.yE + centerPosY,  scale * b[id].g.zE + centerPosZ, // bottom left

																																																																																																																															   // bottom face
																																																																																																																															   scale * b[id].g.xS + centerPosX, scale * b[id].g.yS + centerPosY,  scale * b[id].g.zE + centerPosZ, // top left
																																																																																																																															   scale * b[id].g.xS + centerPosX, scale * b[id].g.yS + centerPosY,  scale * b[id].g.zS + centerPosZ, // top right
																																																																																																																															   scale * b[id].g.xE + centerPosX,  scale * b[id].g.yS + centerPosY,  scale * b[id].g.zS + centerPosZ, // bottom right
																																																																																																																															   scale * b[id].g.xE + centerPosX,  scale * b[id].g.yS + centerPosY,  scale * b[id].g.zE + centerPosZ, // bottom left

		};

		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		glEnableClientState(GL_VERTEX_ARRAY);
		glVertexPointer(3, GL_FLOAT, 0, vertices);
		glDrawArrays(GL_QUADS, 0, 24);
		glDisableClientState(GL_VERTEX_ARRAY);

	}

}

void Drawnvtx(int** &nvtx, TOCHKA* &pa, int id, GLfloat centerPosX, GLfloat centerPosY, GLfloat  centerPosZ, GLfloat scale) {



	GLfloat vertices[] =
	{
		// front face
		scale * pa[nvtx[0][id] - 1].x + centerPosX, scale * pa[nvtx[2][id] - 1].y + centerPosY, scale * pa[nvtx[4][id] - 1].z + centerPosZ, // top left
		scale * pa[nvtx[1][id] - 1].x + centerPosX, scale * pa[nvtx[2][id] - 1].y + centerPosY, scale * pa[nvtx[4][id] - 1].z + centerPosZ, // top right
		scale * pa[nvtx[1][id] - 1].x + centerPosX, scale * pa[nvtx[0][id] - 1].y + centerPosY, scale * pa[nvtx[4][id] - 1].z + centerPosZ, // bottom right
		scale * pa[nvtx[0][id] - 1].x + centerPosX, scale * pa[nvtx[0][id] - 1].y + centerPosY, scale * pa[nvtx[4][id] - 1].z + centerPosZ, // bottom left

																																			// back face
																																			scale * pa[nvtx[0][id] - 1].x + centerPosX, scale * pa[nvtx[2][id] - 1].y + centerPosY, scale * pa[nvtx[0][id] - 1].z + centerPosZ, // top left
																																			scale * pa[nvtx[1][id] - 1].x + centerPosX, scale * pa[nvtx[2][id] - 1].y + centerPosY,  scale * pa[nvtx[0][id] - 1].z + centerPosZ, // top right
																																			scale * pa[nvtx[1][id] - 1].x + centerPosX, scale * pa[nvtx[0][id] - 1].y + centerPosY, scale * pa[nvtx[0][id] - 1].z + centerPosZ, // bottom right
																																			scale * pa[nvtx[0][id] - 1].x + centerPosX, scale * pa[nvtx[0][id] - 1].y + centerPosY, scale * pa[nvtx[0][id] - 1].z + centerPosZ, // bottom left

																																																																				// left face
																																																																				scale * pa[nvtx[0][id] - 1].x + centerPosX,  scale * pa[nvtx[2][id] - 1].y + centerPosY,  scale * pa[nvtx[4][id] - 1].z + centerPosZ, // top left
																																																																				scale * pa[nvtx[0][id] - 1].x + centerPosX, scale * pa[nvtx[2][id] - 1].y + centerPosY,  scale * pa[nvtx[0][id] - 1].z + centerPosZ, // top right
																																																																				scale * pa[nvtx[0][id] - 1].x + centerPosX,  scale * pa[nvtx[0][id] - 1].y + centerPosY, scale * pa[nvtx[0][id] - 1].z + centerPosZ, // bottom right
																																																																				scale * pa[nvtx[0][id] - 1].x + centerPosX, scale * pa[nvtx[0][id] - 1].y + centerPosY, scale * pa[nvtx[4][id] - 1].z + centerPosZ, // bottom left

																																																																																																					// right face
																																																																																																					scale * pa[nvtx[1][id] - 1].x + centerPosX,  scale * pa[nvtx[2][id] - 1].y + centerPosY, scale * pa[nvtx[4][id] - 1].z + centerPosZ, // top left
																																																																																																					scale * pa[nvtx[1][id] - 1].x + centerPosX,  scale * pa[nvtx[2][id] - 1].y + centerPosY,  scale * pa[nvtx[0][id] - 1].z + centerPosZ, // top right
																																																																																																					scale * pa[nvtx[1][id] - 1].x + centerPosX, scale * pa[nvtx[0][id] - 1].y + centerPosY,  scale * pa[nvtx[0][id] - 1].z + centerPosZ, // bottom right
																																																																																																					scale * pa[nvtx[1][id] - 1].x + centerPosX,  scale * pa[nvtx[0][id] - 1].y + centerPosY,  scale * pa[nvtx[4][id] - 1].z + centerPosZ, // bottom left

																																																																																																																																						  // top face
																																																																																																																																						  scale * pa[nvtx[0][id] - 1].x + centerPosX,  scale * pa[nvtx[2][id] - 1].y + centerPosY,  scale * pa[nvtx[4][id] - 1].z + centerPosZ, // top left
																																																																																																																																						  scale * pa[nvtx[0][id] - 1].x + centerPosX,  scale * pa[nvtx[2][id] - 1].y + centerPosY,  scale * pa[nvtx[0][id] - 1].z + centerPosZ, // top right
																																																																																																																																						  scale * pa[nvtx[1][id] - 1].x + centerPosX,  scale * pa[nvtx[2][id] - 1].y + centerPosY,  scale * pa[nvtx[0][id] - 1].z + centerPosZ, // bottom right
																																																																																																																																						  scale * pa[nvtx[1][id] - 1].x + centerPosX,  scale * pa[nvtx[2][id] - 1].y + centerPosY,  scale * pa[nvtx[4][id] - 1].z + centerPosZ, // bottom left

																																																																																																																																																																								// bottom face
																																																																																																																																																																								scale * pa[nvtx[0][id] - 1].x + centerPosX, scale * pa[nvtx[0][id] - 1].y + centerPosY,  scale * pa[nvtx[4][id] - 1].z + centerPosZ, // top left
																																																																																																																																																																								scale * pa[nvtx[0][id] - 1].x + centerPosX, scale * pa[nvtx[0][id] - 1].y + centerPosY,  scale * pa[nvtx[0][id] - 1].z + centerPosZ, // top right
																																																																																																																																																																								scale * pa[nvtx[1][id] - 1].x + centerPosX,  scale * pa[nvtx[0][id] - 1].y + centerPosY,  scale * pa[nvtx[0][id] - 1].z + centerPosZ, // bottom right
																																																																																																																																																																								scale * pa[nvtx[1][id] - 1].x + centerPosX,  scale * pa[nvtx[0][id] - 1].y + centerPosY,  scale * pa[nvtx[4][id] - 1].z + centerPosZ, // bottom left

	};

	//glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	glEnableClientState(GL_VERTEX_ARRAY);
	glVertexPointer(3, GL_FLOAT, 0, vertices);
	glDrawArrays(GL_QUADS, 0, 24);
	glDisableClientState(GL_VERTEX_ARRAY);

}

void DrawCylinder(BLOCK*& b, integer lb, integer id, GLfloat centerPosX, GLfloat centerPosY, GLfloat  centerPosZ, GLfloat scale)
{

	GLfloat angle33, dx33, dy33;
	GLfloat Hcyl = scale * b[id].g.Hcyl;
	GLfloat xC = scale * b[id].g.xC + centerPosX;
	GLfloat yC = scale * b[id].g.yC + centerPosY;
	GLfloat zC = scale * b[id].g.zC + centerPosZ;
	GLfloat R_out_cyl = scale * b[id].g.R_out_cyl;
	GLfloat R_out_cyl2 = 0.0;// scale* b[id].g.R_out_cyl;
	GLfloat R_in_cyl = scale * b[id].g.R_in_cyl;
	GLfloat R_in_cyl2 = 0.0; // scale* b[id].g.R_in_cyl;

							 // Cylinder
	switch (b[id].g.iPlane) {
	case XY_PLANE:
		// XY
		glBegin(GL_LINE_LOOP);
		for (int i33 = 0; i33 <= 29; i33++) {

			angle33 = 2.0 * 3.1415926 * i33 / 29.0;
			dx33 = R_out_cyl * cos(angle33);
			dy33 = R_out_cyl * sin(angle33);
			glVertex3f(xC + dx33, yC + dy33, zC);
		}
		glEnd();
		if (R_out_cyl2 <= 0.0) {

			glBegin(GL_LINE_LOOP);
			for (int i33 = 0; i33 <= 29; i33++) {

				angle33 = 2.0 * 3.1415926 * i33 / 29.0;
				dx33 = R_out_cyl * cos(angle33);
				dy33 = R_out_cyl * sin(angle33);
				glVertex3f(xC + dx33, yC + dy33, zC + Hcyl);
			}
			glEnd();
		}
		else
		{
			glBegin(GL_LINE_LOOP);
			for (int i33 = 0; i33 <= 29; i33++) {

				angle33 = 2.0 * 3.1415926 * i33 / 29.0;
				dx33 = R_out_cyl2 * cos(angle33);
				dy33 = R_out_cyl2 * sin(angle33);
				glVertex3f(xC + dx33, yC + dy33, zC + Hcyl);
			}
			glEnd();
		}

		if (R_in_cyl > 1.0e-20) {

			glBegin(GL_LINE_LOOP);
			for (int i33 = 0; i33 <= 29; i33++) {

				angle33 = 2.0 * 3.1415926 * i33 / 29.0;
				dx33 = R_in_cyl * cos(angle33);
				dy33 = R_in_cyl * sin(angle33);
				glVertex3f(xC + dx33, yC + dy33, zC);
			}
			glEnd();

			if (R_in_cyl2 <= 0.0) {

				glBegin(GL_LINE_LOOP);
				for (int i33 = 0; i33 <= 29; i33++) {

					angle33 = 2.0 * 3.1415926 * i33 / 29.0;
					dx33 = R_in_cyl * cos(angle33);
					dy33 = R_in_cyl * sin(angle33);
					glVertex3f(xC + dx33, yC + dy33, zC + Hcyl);
				}
				glEnd();
			}
			else
			{
				glBegin(GL_LINE_LOOP);
				for (int i33 = 0; i33 <= 29; i33++) {

					angle33 = 2.0 * 3.1415926 * i33 / 29.0;
					dx33 = R_in_cyl2 * cos(angle33);
					dy33 = R_in_cyl2 * sin(angle33);
					glVertex3f(xC + dx33, yC + dy33, zC + Hcyl);
				}
				glEnd();
			}
		}

		glBegin(GL_LINE_LOOP);
		angle33 = 2.0 * 3.1415926 * 0.125;
		dx33 = R_out_cyl * cos(angle33);
		dy33 = R_out_cyl * sin(angle33);
		glVertex3f(xC + dx33, yC + dy33, zC);
		if (R_out_cyl2 <= 0.0) {

			glVertex3f(xC + dx33, yC + dy33, zC + Hcyl);
		}
		else
		{
			dx33 = R_out_cyl2 * cos(angle33);
			dy33 = R_out_cyl2 * sin(angle33);
			glVertex3f(xC + dx33, yC + dy33, zC + Hcyl);
		}
		glEnd();
		glBegin(GL_LINE_LOOP);
		angle33 = 2.0 * 3.1415926 * 0.375;
		dx33 = R_out_cyl * cos(angle33);
		dy33 = R_out_cyl * sin(angle33);
		glVertex3f(xC + dx33, yC + dy33, zC);
		if (R_out_cyl2 <= 0.0) {

			glVertex3f(xC + dx33, yC + dy33, zC + Hcyl);
		}
		else
		{
			dx33 = R_out_cyl2 * cos(angle33);
			dy33 = R_out_cyl2 * sin(angle33);
			glVertex3f(xC + dx33, yC + dy33, zC + Hcyl);
		}
		glEnd();
		glBegin(GL_LINE_LOOP);
		angle33 = 2.0 * 3.1415926 * 0.625;
		dx33 = R_out_cyl * cos(angle33);
		dy33 = R_out_cyl * sin(angle33);
		glVertex3f(xC + dx33, yC + dy33, zC);
		if (R_out_cyl2 <= 0.0) {

			glVertex3f(xC + dx33, yC + dy33, zC + Hcyl);
		}
		else
		{
			dx33 = R_out_cyl2 * cos(angle33);
			dy33 = R_out_cyl2 * sin(angle33);
			glVertex3f(xC + dx33, yC + dy33, zC + Hcyl);
		}
		glEnd();
		glBegin(GL_LINE_LOOP);
		angle33 = 2.0 * 3.1415926 * 0.875;
		dx33 = R_out_cyl * cos(angle33);
		dy33 = R_out_cyl * sin(angle33);
		glVertex3f(xC + dx33, yC + dy33, zC);
		if (R_out_cyl2 <= 0.0) {

			glVertex3f(xC + dx33, yC + dy33, zC + Hcyl);
		}
		else
		{
			dx33 = R_out_cyl2 * cos(angle33);
			dy33 = R_out_cyl2 * sin(angle33);
			glVertex3f(xC + dx33, yC + dy33, zC + Hcyl);
		}
		glEnd();

		if (R_in_cyl > 1.0e-20) {

			glBegin(GL_LINE_LOOP);
			angle33 = 2.0 * 3.1415926 * 0.125;
			dx33 = R_in_cyl * cos(angle33);
			dy33 = R_in_cyl * sin(angle33);
			glVertex3f(xC + dx33, yC + dy33, zC);
			if (R_in_cyl2 <= 0.0) {

				glVertex3f(xC + dx33, yC + dy33, zC + Hcyl);
			}
			else
			{
				dx33 = R_in_cyl2 * cos(angle33);
				dy33 = R_in_cyl2 * sin(angle33);
				glVertex3f(xC + dx33, yC + dy33, zC + Hcyl);
			}
			glEnd();
			glBegin(GL_LINE_LOOP);
			angle33 = 2.0 * 3.1415926 * 0.375;
			dx33 = R_in_cyl * cos(angle33);
			dy33 = R_in_cyl * sin(angle33);
			glVertex3f(xC + dx33, yC + dy33, zC);
			if (R_in_cyl2 <= 0.0) {

				glVertex3f(xC + dx33, yC + dy33, zC + Hcyl);
			}
			else
			{
				dx33 = R_in_cyl2 * cos(angle33);
				dy33 = R_in_cyl2 * sin(angle33);
				glVertex3f(xC + dx33, yC + dy33, zC + Hcyl);
			}
			glEnd();
			glBegin(GL_LINE_LOOP);
			angle33 = 2.0 * 3.1415926 * 0.625;
			dx33 = R_in_cyl * cos(angle33);
			dy33 = R_in_cyl * sin(angle33);
			glVertex3f(xC + dx33, yC + dy33, zC);
			if (R_in_cyl2 <= 0.0) {

				glVertex3f(xC + dx33, yC + dy33, zC + Hcyl);
			}
			else
			{
				dx33 = R_in_cyl2 * cos(angle33);
				dy33 = R_in_cyl2 * sin(angle33);
				glVertex3f(xC + dx33, yC + dy33, zC + Hcyl);
			}
			glEnd();
			glBegin(GL_LINE_LOOP);
			angle33 = 2.0 * 3.1415926 * 0.875;
			dx33 = R_in_cyl * cos(angle33);
			dy33 = R_in_cyl * sin(angle33);
			glVertex3f(xC + dx33, yC + dy33, zC);
			if (R_in_cyl2 <= 0.0) {

				glVertex3f(xC + dx33, yC + dy33, zC + Hcyl);
			}
			else
			{
				dx33 = R_in_cyl2 * cos(angle33);
				dy33 = R_in_cyl2 * sin(angle33);
				glVertex3f(xC + dx33, yC + dy33, zC + Hcyl);
			}
			glEnd();
		}

		break;
	case XZ_PLANE:
		// XZ
		glBegin(GL_LINE_LOOP);
		for (int i33 = 0; i33 <= 29; i33++) {

			angle33 = 2.0 * 3.1415926 * i33 / 29.0;
			dx33 = R_out_cyl * cos(angle33);
			dy33 = R_out_cyl * sin(angle33);
			glVertex3f(xC + dx33, yC, zC + dy33);
		}
		glEnd();
		if (R_out_cyl2 <= 0.0) {

			glBegin(GL_LINE_LOOP);
			for (int i33 = 0; i33 <= 29; i33++) {

				angle33 = 2.0 * 3.1415926 * i33 / 29.0;
				dx33 = R_out_cyl * cos(angle33);
				dy33 = R_out_cyl * sin(angle33);
				glVertex3f(xC + dx33, yC + Hcyl, zC + dy33);
			}
			glEnd();
		}
		else
		{
			glBegin(GL_LINE_LOOP);
			for (int i33 = 0; i33 <= 29; i33++) {

				angle33 = 2.0 * 3.1415926 * i33 / 29.0;
				dx33 = R_out_cyl2 * cos(angle33);
				dy33 = R_out_cyl2 * sin(angle33);
				glVertex3f(xC + dx33, yC + Hcyl, zC + dy33);
			}
			glEnd();
		}

		if (R_in_cyl > 1.0e-20) {

			glBegin(GL_LINE_LOOP);
			for (int i33 = 0; i33 <= 29; i33++) {

				angle33 = 2.0 * 3.1415926 * i33 / 29.0;
				dx33 = R_in_cyl * cos(angle33);
				dy33 = R_in_cyl * sin(angle33);
				glVertex3f(xC + dx33, yC, zC + dy33);
			}
			glEnd();
			if (R_in_cyl2 <= 0.0) {

				glBegin(GL_LINE_LOOP);
				for (int i33 = 0; i33 <= 29; i33++) {

					angle33 = 2.0 * 3.1415926 * i33 / 29.0;
					dx33 = R_in_cyl * cos(angle33);
					dy33 = R_in_cyl * sin(angle33);
					glVertex3f(xC + dx33, yC + Hcyl, zC + dy33);
				}
				glEnd();
			}
			else
			{
				glBegin(GL_LINE_LOOP);
				for (int i33 = 0; i33 <= 29; i33++) {

					angle33 = 2.0 * 3.1415926 * i33 / 29.0;
					dx33 = R_in_cyl2 * cos(angle33);
					dy33 = R_in_cyl2 * sin(angle33);
					glVertex3f(xC + dx33, yC + Hcyl, zC + dy33);
				}
				glEnd();
			}
		}

		glBegin(GL_LINE_LOOP);
		angle33 = 2.0 * 3.1415926 * 0.125;
		dx33 = R_out_cyl * cos(angle33);
		dy33 = R_out_cyl * sin(angle33);
		glVertex3f(xC + dx33, yC, zC + dy33);
		if (R_out_cyl2 <= 0.0) {

			glVertex3f(xC + dx33, yC + Hcyl, zC + dy33);
		}
		else
		{
			dx33 = R_out_cyl2 * cos(angle33);
			dy33 = R_out_cyl2 * sin(angle33);
			glVertex3f(xC + dx33, yC + Hcyl, zC + dy33);
		}
		glEnd();
		glBegin(GL_LINE_LOOP);
		angle33 = 2.0 * 3.1415926 * 0.375;
		dx33 = R_out_cyl * cos(angle33);
		dy33 = R_out_cyl * sin(angle33);
		glVertex3f(xC + dx33, yC, zC + dy33);
		if (R_out_cyl2 <= 0.0) {

			glVertex3f(xC + dx33, yC + Hcyl, zC + dy33);
		}
		else
		{
			dx33 = R_out_cyl2 * cos(angle33);
			dy33 = R_out_cyl2 * sin(angle33);
			glVertex3f(xC + dx33, yC + Hcyl, zC + dy33);
		}
		glEnd();
		glBegin(GL_LINE_LOOP);
		angle33 = 2.0 * 3.1415926 * 0.625;
		dx33 = R_out_cyl * cos(angle33);
		dy33 = R_out_cyl * sin(angle33);
		glVertex3f(xC + dx33, yC, zC + dy33);
		if (R_out_cyl2 <= 0.0) {

			glVertex3f(xC + dx33, yC + Hcyl, zC + dy33);
		}
		else
		{
			dx33 = R_out_cyl2 * cos(angle33);
			dy33 = R_out_cyl2 * sin(angle33);
			glVertex3f(xC + dx33, yC + Hcyl, zC + dy33);
		}
		glEnd();
		glBegin(GL_LINE_LOOP);
		angle33 = 2.0 * 3.1415926 * 0.875;
		dx33 = R_out_cyl * cos(angle33);
		dy33 = R_out_cyl * sin(angle33);
		glVertex3f(xC + dx33, yC, zC + dy33);
		if (R_out_cyl2 <= 0.0) {

			glVertex3f(xC + dx33, yC + Hcyl, zC + dy33);
		}
		else
		{
			dx33 = R_out_cyl2 * cos(angle33);
			dy33 = R_out_cyl2 * sin(angle33);
			glVertex3f(xC + dx33, yC + Hcyl, zC + dy33);
		}
		glEnd();

		if (R_in_cyl > 1.0e-20) {

			glBegin(GL_LINE_LOOP);
			angle33 = 2.0 * 3.1415926 * 0.125;
			dx33 = R_in_cyl * cos(angle33);
			dy33 = R_in_cyl * sin(angle33);
			glVertex3f(xC + dx33, yC, zC + dy33);
			if (R_in_cyl2 <= 0.0) {

				glVertex3f(xC + dx33, yC + Hcyl, zC + dy33);
			}
			else
			{
				dx33 = R_in_cyl2 * cos(angle33);
				dy33 = R_in_cyl2 * sin(angle33);
				glVertex3f(xC + dx33, yC + Hcyl, zC + dy33);
			};
			glEnd();
			glBegin(GL_LINE_LOOP);
			angle33 = 2.0 * 3.1415926 * 0.375;
			dx33 = R_in_cyl * cos(angle33);
			dy33 = R_in_cyl * sin(angle33);
			glVertex3f(xC + dx33, yC, zC + dy33);
			if (R_in_cyl2 <= 0.0) {

				glVertex3f(xC + dx33, yC + Hcyl, zC + dy33);
			}
			else
			{
				dx33 = R_in_cyl2 * cos(angle33);
				dy33 = R_in_cyl2 * sin(angle33);
				glVertex3f(xC + dx33, yC + Hcyl, zC + dy33);
			}
			glEnd();
			glBegin(GL_LINE_LOOP);
			angle33 = 2.0 * 3.1415926 * 0.625;
			dx33 = R_in_cyl * cos(angle33);
			dy33 = R_in_cyl * sin(angle33);
			glVertex3f(xC + dx33, yC, zC + dy33);
			if (R_in_cyl2 <= 0.0) {

				glVertex3f(xC + dx33, yC + Hcyl, zC + dy33);
			}
			else
			{
				dx33 = R_in_cyl2 * cos(angle33);
				dy33 = R_in_cyl2 * sin(angle33);
				glVertex3f(xC + dx33, yC + Hcyl, zC + dy33);
			}
			glEnd();
			glBegin(GL_LINE_LOOP);
			angle33 = 2.0 * 3.1415926 * 0.875;
			dx33 = R_in_cyl * cos(angle33);
			dy33 = R_in_cyl * sin(angle33);
			glVertex3f(xC + dx33, yC, zC + dy33);
			if (R_in_cyl2 <= 0.0) {

				glVertex3f(xC + dx33, yC + Hcyl, zC + dy33);
			}
			else
			{
				dx33 = R_in_cyl2 * cos(angle33);
				dy33 = R_in_cyl2 * sin(angle33);
				glVertex3f(xC + dx33, yC + Hcyl, zC + dy33);
			}
			glEnd();
		}


		break;
	case YZ_PLANE:
		// YZ
		glBegin(GL_LINE_LOOP);
		for (int i33 = 0; i33 <= 29; i33++) {

			angle33 = 2.0 * 3.1415926 * i33 / 29.0;
			dx33 = R_out_cyl * cos(angle33);
			dy33 = R_out_cyl * sin(angle33);
			glVertex3f(xC, yC + dx33, zC + dy33);
		}
		glEnd();
		if (R_out_cyl2 <= 0.0) {

			glBegin(GL_LINE_LOOP);
			for (int i33 = 0; i33 <= 29; i33++) {

				angle33 = 2.0 * 3.1415926 * i33 / 29.0;
				dx33 = R_out_cyl * cos(angle33);
				dy33 = R_out_cyl * sin(angle33);
				glVertex3f(xC + Hcyl, yC + dx33, zC + dy33);
			}
			glEnd();
		}
		else
		{
			glBegin(GL_LINE_LOOP);
			for (int i33 = 0; i33 <= 29; i33++) {

				angle33 = 2.0 * 3.1415926 * i33 / 29.0;
				dx33 = R_out_cyl2 * cos(angle33);
				dy33 = R_out_cyl2 * sin(angle33);
				glVertex3f(xC + Hcyl, yC + dx33, zC + dy33);
			}
			glEnd();
		}
		if (R_in_cyl > 1.0e-20) {

			glBegin(GL_LINE_LOOP);
			for (int i33 = 0; i33 <= 29; i33++) {

				angle33 = 2.0 * 3.1415926 * i33 / 29.0;
				dx33 = R_in_cyl * cos(angle33);
				dy33 = R_in_cyl * sin(angle33);
				glVertex3f(xC, yC + dx33, zC + dy33);
			}
			glEnd();
			if (R_in_cyl2 <= 0.0) {

				glBegin(GL_LINE_LOOP);
				for (int i33 = 0; i33 <= 29; i33++) {

					angle33 = 2.0 * 3.1415926 * i33 / 29.0;
					dx33 = R_in_cyl * cos(angle33);
					dy33 = R_in_cyl * sin(angle33);
					glVertex3f(xC + Hcyl, yC + dx33, zC + dy33);
				}
				glEnd();
			}
			else
			{
				glBegin(GL_LINE_LOOP);
				for (int i33 = 0; i33 <= 29; i33++) {

					angle33 = 2.0 * 3.1415926 * i33 / 29.0;
					dx33 = R_in_cyl2 * cos(angle33);
					dy33 = R_in_cyl2 * sin(angle33);
					glVertex3f(xC + Hcyl, yC + dx33, zC + dy33);
				}
				glEnd();
			}
		}

		glBegin(GL_LINE_LOOP);
		angle33 = 2.0 * 3.1415926 * 0.125;
		dx33 = R_out_cyl * cos(angle33);
		dy33 = R_out_cyl * sin(angle33);
		glVertex3f(xC, yC + dx33, zC + dy33);
		if (R_out_cyl2 <= 0.0) {

			glVertex3f(xC + Hcyl, yC + dx33, zC + dy33);
		}
		else
		{
			dx33 = R_out_cyl * cos(angle33);
			dy33 = R_out_cyl * sin(angle33);
			glVertex3f(xC + Hcyl, yC + dx33, zC + dy33);
		}
		glEnd();
		glBegin(GL_LINE_LOOP);
		angle33 = 2.0 * 3.1415926 * 0.375;
		dx33 = R_out_cyl * cos(angle33);
		dy33 = R_out_cyl * sin(angle33);
		glVertex3f(xC, yC + dx33, zC + dy33);
		if (R_out_cyl2 <= 0.0) {

			glVertex3f(xC + Hcyl, yC + dx33, zC + dy33);
		}
		else
		{
			dx33 = R_out_cyl * cos(angle33);
			dy33 = R_out_cyl * sin(angle33);
			glVertex3f(xC + Hcyl, yC + dx33, zC + dy33);
		}
		glEnd();
		glBegin(GL_LINE_LOOP);
		angle33 = 2.0 * 3.1415926 * 0.625;
		dx33 = R_out_cyl * cos(angle33);
		dy33 = R_out_cyl * sin(angle33);
		glVertex3f(xC, yC + dx33, zC + dy33);
		if (R_out_cyl2 <= 0.0) {

			glVertex3f(xC + Hcyl, yC + dx33, zC + dy33);
		}
		else
		{
			dx33 = R_out_cyl * cos(angle33);
			dy33 = R_out_cyl * sin(angle33);
			glVertex3f(xC + Hcyl, yC + dx33, zC + dy33);
		}
		glEnd();
		glBegin(GL_LINE_LOOP);
		angle33 = 2.0 * 3.1415926 * 0.875;
		dx33 = R_out_cyl * cos(angle33);
		dy33 = R_out_cyl * sin(angle33);
		glVertex3f(xC, yC + dx33, zC + dy33);
		if (R_out_cyl2 <= 0.0) {

			glVertex3f(xC + Hcyl, yC + dx33, zC + dy33);
		}
		else
		{
			dx33 = R_out_cyl2 * cos(angle33);
			dy33 = R_out_cyl2 * sin(angle33);
			glVertex3f(xC + Hcyl, yC + dx33, zC + dy33);
		}
		glEnd();

		if (R_in_cyl > 1.0e-20) {

			glBegin(GL_LINE_LOOP);
			angle33 = 2.0 * 3.1415926 * 0.125;
			dx33 = R_in_cyl * cos(angle33);
			dy33 = R_in_cyl * sin(angle33);
			glVertex3f(xC, yC + dx33, zC + dy33);
			if (R_in_cyl2 <= 0.0) {

				glVertex3f(xC + Hcyl, yC + dx33, zC + dy33);
			}
			else
			{
				dx33 = R_in_cyl2 * cos(angle33);
				dy33 = R_in_cyl2 * sin(angle33);
				glVertex3f(xC + Hcyl, yC + dx33, zC + dy33);
			}
			glEnd();
			glBegin(GL_LINE_LOOP);
			angle33 = 2.0 * 3.1415926 * 0.375;
			dx33 = R_in_cyl * cos(angle33);
			dy33 = R_in_cyl * sin(angle33);
			glVertex3f(xC, yC + dx33, zC + dy33);
			if (R_in_cyl2 <= 0.0) {

				glVertex3f(xC + Hcyl, yC + dx33, zC + dy33);
			}
			else
			{
				dx33 = R_in_cyl2 * cos(angle33);
				dy33 = R_in_cyl2 * sin(angle33);
				glVertex3f(xC + Hcyl, yC + dx33, zC + dy33);
			}
			glEnd();
			glBegin(GL_LINE_LOOP);
			angle33 = 2.0 * 3.1415926 * 0.625;
			dx33 = R_in_cyl * cos(angle33);
			dy33 = R_in_cyl * sin(angle33);
			glVertex3f(xC, yC + dx33, zC + dy33);
			if (R_in_cyl2 <= 0.0) {

				glVertex3f(xC + Hcyl, yC + dx33, zC + dy33);
			}
			else
			{
				dx33 = R_in_cyl2 * cos(angle33);
				dy33 = R_in_cyl2 * sin(angle33);
				glVertex3f(xC + Hcyl, yC + dx33, zC + dy33);
			}
			glEnd();
			glBegin(GL_LINE_LOOP);
			angle33 = 2.0 * 3.1415926 * 0.875;
			dx33 = R_in_cyl * cos(angle33);
			dy33 = R_in_cyl * sin(angle33);
			glVertex3f(xC, yC + dx33, zC + dy33);
			if (R_in_cyl2 <= 0.0) {

				glVertex3f(xC + Hcyl, yC + dx33, zC + dy33);
			}
			else
			{
				dx33 = R_in_cyl2 * cos(angle33);
				dy33 = R_in_cyl2 * sin(angle33);
				glVertex3f(xC + Hcyl, yC + dx33, zC + dy33);
			}
			glEnd();
		}
		break;
	}

}

void DrawPolygon(BLOCK*& b, integer lb, integer id, GLfloat centerPosX, GLfloat centerPosY, GLfloat  centerPosZ, GLfloat scale) {



	// Polygon
	switch (b[id].g.iPlane_obj2) {
	case XY_PLANE:
		// XY
		for (int j3 = 0; j3 <= b[id].g.nsizei - 2; j3++)
		{
			// GL_LINE_LOOP
			glBegin(GL_LINE_LOOP);
			glVertex3f(centerPosX + scale * b[id].g.xi[j3], centerPosY + scale * b[id].g.yi[j3], centerPosZ + scale * b[id].g.zi[j3]);
			glVertex3f(centerPosX + scale * b[id].g.xi[j3 + 1], centerPosY + scale * b[id].g.yi[j3 + 1], centerPosZ + scale * b[id].g.zi[j3 + 1]);
			glVertex3f(centerPosX + scale * b[id].g.xi[j3 + 1], centerPosY + scale * b[id].g.yi[j3 + 1], centerPosZ + scale * b[id].g.zi[j3 + 1] + scale * b[id].g.hi[j3 + 1]);
			glVertex3f(centerPosX + scale * b[id].g.xi[j3], centerPosY + scale * b[id].g.yi[j3], centerPosZ + scale * b[id].g.zi[j3] + scale * b[id].g.hi[j3]);
			glEnd();
		}
		glBegin(GL_LINE_LOOP);
		glVertex3f(centerPosX + scale * b[id].g.xi[b[id].g.nsizei - 1], centerPosY + scale * b[id].g.yi[b[id].g.nsizei - 1], centerPosZ + scale * b[id].g.zi[b[id].g.nsizei - 1]);
		glVertex3f(centerPosX + scale * b[id].g.xi[0], centerPosY + scale * b[id].g.yi[0], centerPosZ + scale * b[id].g.zi[0]);
		glVertex3f(centerPosX + scale * b[id].g.xi[0], centerPosY + scale * b[id].g.yi[0], centerPosZ + scale * b[id].g.zi[0] + scale * b[id].g.hi[0]);
		glVertex3f(centerPosX + scale * b[id].g.xi[b[id].g.nsizei - 1], centerPosY + scale * b[id].g.yi[b[id].g.nsizei - 1], centerPosZ + scale * b[id].g.zi[b[id].g.nsizei - 1] + scale * b[id].g.hi[b[id].g.nsizei - 1]);
		glEnd();
		break;
	case XZ_PLANE:
		// XZ
		for (int j3 = 0; j3 <= b[id].g.nsizei - 2; j3++)
		{

			glBegin(GL_LINE_LOOP);
			glVertex3f(centerPosX + scale * b[id].g.xi[j3], centerPosY + scale * b[id].g.yi[j3], centerPosZ + scale * b[id].g.zi[j3]);
			glVertex3f(centerPosX + scale * b[id].g.xi[j3 + 1], centerPosY + scale * b[id].g.yi[j3 + 1], centerPosZ + scale * b[id].g.zi[j3 + 1]);
			glVertex3f(centerPosX + scale * b[id].g.xi[j3 + 1], centerPosY + scale * b[id].g.yi[j3 + 1] + scale * b[id].g.hi[j3 + 1], centerPosZ + scale * b[id].g.zi[j3 + 1]);
			glVertex3f(centerPosX + scale * b[id].g.xi[j3], centerPosY + scale * b[id].g.yi[j3] + scale * b[id].g.hi[j3], centerPosZ + scale * b[id].g.zi[j3]);
			glEnd();
		}
		glBegin(GL_LINE_LOOP);
		glVertex3f(centerPosX + scale * b[id].g.xi[b[id].g.nsizei - 1], centerPosY + scale * b[id].g.yi[b[id].g.nsizei - 1], centerPosZ + scale * b[id].g.zi[b[id].g.nsizei - 1]);
		glVertex3f(centerPosX + scale * b[id].g.xi[0], centerPosY + scale * b[id].g.yi[0], centerPosZ + scale * b[id].g.zi[0]);
		glVertex3f(centerPosX + scale * b[id].g.xi[0], centerPosY + scale * b[id].g.yi[0] + scale * b[id].g.hi[0], centerPosZ + scale * b[id].g.zi[0]);
		glVertex3f(centerPosX + scale * b[id].g.xi[b[id].g.nsizei - 1], centerPosY + scale * b[id].g.yi[b[id].g.nsizei - 1] + scale * b[id].g.hi[b[id].g.nsizei - 1], centerPosZ + scale * b[id].g.zi[b[id].g.nsizei - 1]);
		glEnd();

		break;
	case YZ_PLANE:
		// YZ
		for (int j3 = 0; j3 <= b[id].g.nsizei - 2; j3++)
		{
			glBegin(GL_LINE_LOOP);
			glVertex3f(centerPosX + scale * b[id].g.xi[j3], centerPosY + scale * b[id].g.yi[j3], centerPosZ + scale * b[id].g.zi[j3]);
			glVertex3f(centerPosX + scale * b[id].g.xi[j3 + 1], centerPosY + scale * b[id].g.yi[j3 + 1], centerPosZ + scale * b[id].g.zi[j3 + 1]);
			glVertex3f(centerPosX + scale * b[id].g.xi[j3 + 1] + scale * b[id].g.hi[j3 + 1], centerPosY + scale * b[id].g.yi[j3 + 1], centerPosZ + scale * b[id].g.zi[j3 + 1]);
			glVertex3f(centerPosX + scale * b[id].g.xi[j3] + scale * b[id].g.hi[j3], centerPosY + scale * b[id].g.yi[j3], centerPosZ + scale * b[id].g.zi[j3]);
			glEnd();
		}
		glBegin(GL_LINE_LOOP);
		glVertex3f(centerPosX + scale * b[id].g.xi[b[id].g.nsizei - 1], centerPosY + scale * b[id].g.yi[b[id].g.nsizei - 1], centerPosZ + scale * b[id].g.zi[b[id].g.nsizei - 1]);
		glVertex3f(centerPosX + scale * b[id].g.xi[0], centerPosY + scale * b[id].g.yi[0], centerPosZ + scale * b[id].g.zi[0]);
		glVertex3f(centerPosX + scale * b[id].g.xi[0] + scale * b[id].g.hi[0], centerPosY + scale * b[id].g.yi[0], centerPosZ + scale * b[id].g.zi[0]);
		glVertex3f(centerPosX + scale * b[id].g.xi[b[id].g.nsizei - 1] + scale * b[id].g.hi[b[id].g.nsizei - 1], centerPosY + scale * b[id].g.yi[b[id].g.nsizei - 1], centerPosZ + scale * b[id].g.zi[b[id].g.nsizei - 1]);
		glEnd();
		break;

	}
}

#endif