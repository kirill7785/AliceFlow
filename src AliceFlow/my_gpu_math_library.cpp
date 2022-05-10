// my_gpu_math_library.cu
// Реализация операций с длинными векторами на видеокарте NVIDIA.
// 30.07.2021


#pragma once
#ifndef MY_GPU_MATH_LIBRARY_CU
#define MY_GPU_MATH_LIBRARY_CU 1



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


// false - используется центральный процессор.
 // Вычисление произведения матрицы на вектор происходит на видеокарте.
// Каждый раз данные передаются на видеокарту и освобождаются.
// Работает при включенном  MY_GPU_MATH_LIBRARY_CU_ON
#define MY_GPU_MATH_LIBRARY_CU_REALLOC_ON false
// Данные на видеокарту передаются минимальное число раз фактически единожды. 
// Во время потребности в вычислении произведения матрицы на вектор память видеокарты постоянно занята.
#define MY_GPU_MATH_LIBRARY_CU_ON  false 
#define EllFormat  true

// true если девайс ещё не инициализирован.
bool init_b_first_device = true;

#if MY_GPU_MATH_LIBRARY_CU_ON

#include <device_functions.h> // для __syncthreads();

#include <thrust/reduce.h>
#include <thrust/functional.h>

#include <thrust/inner_product.h>
#include <thrust/execution_policy.h>

cudaError_t cudaStatus;

__global__ void initKernel(doublereal* x, integer n, doublereal x0)
{
	int i = threadIdx.x + blockIdx.x * blockDim.x;
	while (i < n) {
		x[i] = x0;
		i += blockDim.x * gridDim.x;
	}

}
#endif

// bdevice == true запуск на GPU.
// инициализирует вектор x значением x0.
void init_v(doublereal*& x, integer n, doublereal x0, bool bdevice)
{

	if ((MY_GPU_MATH_LIBRARY_CU_ON) && (bdevice)) {

#if MY_GPU_MATH_LIBRARY_CU_ON	

		if (init_b_first_device) {

			

			int device;
			device = idevice_Tesla;

			// Choose which GPU to run on, change this on a multi-GPU system.
			cudaStatus = cudaSetDevice(device);
			if (cudaStatus != cudaSuccess) {
				fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?");
				system("PAUSE");
				exit(1);
			}

			cudaDeviceProp prop;
			cudaGetDeviceProperties(&prop, device);
			printf("%s\n", prop.name);
			// GeForce 840M 384потока 1ГГц каждый март 2014 года. 28нм.

			init_b_first_device = false;
		}

		// на видеокарте GPU.

		doublereal* dev_x = nullptr;

		// Allocate GPU buffers for three vectors (two input, one output)    .
		cudaStatus = cudaMalloc((void**)&dev_x, (n) * sizeof(doublereal));
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMalloc dev_x failed!");
			system("PAUSE");
			exit(1);
		}

		// Copy input vectors from host memory to GPU buffers.
		cudaStatus = cudaMemcpy(dev_x, x, n * sizeof(doublereal), cudaMemcpyHostToDevice);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMemcpy dev_x, x HostToDevice failed!");
			system("PAUSE");
			exit(1);
		}

		initKernel<<<128, 128 >>>(dev_x,n,x0);
		// Check for any errors launching the kernel
		cudaStatus = cudaGetLastError();
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "initKernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
			system("PAUSE");
			exit(1);
		}

		cudaStatus = cudaMemcpy(x, dev_x, n * sizeof(doublereal), cudaMemcpyDeviceToHost);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMemcpy x, dev_x DeviceToHost failed!");
			system("PAUSE");
			exit(1);
		}

		cudaFree(dev_x);

#endif

	}
	else {
		// на процессоре CPU.
#pragma omp parallel for
		for (integer i = 0; i < n; ++i) {
			x[i] = x0;
		}
	}
}

#if MY_GPU_MATH_LIBRARY_CU_ON

// реализация неэффективна как пишут в интернете
// Скалярная реализация.
__global__ void MatrixCRSByVectorKernelScalar(const doublereal* val, const integer* col_ind, const integer* row_ptr, const doublereal* V, doublereal* tmp, const integer n)
{
	int i = threadIdx.x + blockIdx.x * blockDim.x;
	while (i < n) {
		//tmp[i] = 0.0;

		doublereal sum=0.0;
		const integer rowend = row_ptr[i + 1];
		const integer rowbeg = row_ptr[i];


		for (integer j = rowbeg; j < rowend; ++j)
		{
			sum += val[j] * V[col_ind[j]];
		}
		tmp[i] = sum;


		i += blockDim.x * gridDim.x;
	}

}

#define FULL_WARP_MASK 0xFFFFFFFF

template <class T>
__device__ T warp_reduce(T val)
{

	for (integer offset = warpSize / 2; offset > 0; offset /= 2)
	{
		val += __shfl_down_sync(FULL_WARP_MASK, val, offset);
	}

	return val;
}


// Векторная реализация.
__global__ void MatrixCRSByVectorKernel(const doublereal* val, const integer* col_ind, const integer* row_ptr, const doublereal* V, doublereal* tmp, const integer n)
{
	int i = threadIdx.x + blockIdx.x * blockDim.x;
    int warp_id = i / 32;
	int lane = i % 32;

	int row = warp_id; //// One warp per raw

	

	while (row < n) {
		//tmp[i] = 0.0;

		doublereal sum = 0.0;

		const integer rowend = row_ptr[row + 1]; 
		const integer rowbeg = row_ptr[row];

		
		for (integer j = rowbeg + lane; j < rowend; j += 32)
		{
			sum += val[j] * V[col_ind[j]];
		}

		sum = warp_reduce(sum);

		if (lane == 0 && row < n) {
			tmp[row] = sum;
		}	


		i += blockDim.x * gridDim.x;
		warp_id = i / 32;
		lane = i % 32;

		row = warp_id; //// One warp per raw

	}

}


// Ellpack Itpack
__global__ void MatrixEllPackItpackByVectorKernel(const doublereal* data, const integer* indices, const integer string_size, const doublereal* V, doublereal* tmp, const integer n, const integer iadd)
{
	int i = threadIdx.x + blockIdx.x * blockDim.x;
	while (i < n) {
		//tmp[i] = 0.0;

		doublereal sum = 0.0;
		
		for (integer j = 0; j < string_size; ++j)
		{
			const integer element_offset = i + j * n;
			const doublereal val = data[element_offset];
			//if (val != 0)
			sum +=  val * V[indices[element_offset]];
		}
		tmp[i+ iadd] = sum;


		i += blockDim.x * gridDim.x;
	}

}
#endif

void MatrixCRSByVector(doublereal*& val, integer*& col_ind, integer*& row_ptr, doublereal*& V, doublereal*& tmp, const integer n)
{

	if (MY_GPU_MATH_LIBRARY_CU_ON) {

#if MY_GPU_MATH_LIBRARY_CU_ON

		if (init_b_first_device) {

			int device;
			device = idevice_Tesla;

			// Choose which GPU to run on, change this on a multi-GPU system.
			cudaStatus = cudaSetDevice(device);
			if (cudaStatus != cudaSuccess) {
				fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?");
				system("PAUSE");
				exit(1);
			}

			cudaDeviceProp prop;
			cudaGetDeviceProperties(&prop, device);
			printf("%s\n", prop.name);
			// GeForce 840M 384потока 1ГГц каждый март 2014 года. 28нм.

			init_b_first_device = false;
		}


#if EllFormat
		doublereal* data = nullptr;
		integer* indices = nullptr;
		integer string_size = -1;

		

		for (integer i77 = 0; i77 < n; ++i77) {
			if (row_ptr[i77 + 1] - row_ptr[i77] > string_size) {
				string_size = row_ptr[i77 + 1] - row_ptr[i77];
			}
			//if (i77>= maxelm) {
				//std::cout << "\nstart hvost=" << i77 << "  maxelm=" << maxelm << "start "<< row_ptr75[i77]  << "end=" << row_ptr75[i77 + 1] << std::endl;
				//system("pause");
			//}
		}

		data = new doublereal[string_size * n];
		indices = new integer[string_size * n];

		


		// Перепаковка в Ell формат.
#pragma omp parallel for
		for (integer i77 = 0; i77 < string_size * n; ++i77) {
			data[i77] = 0.0;
		}



#pragma omp parallel for
		for (integer i77 = 0; i77 < n; ++i77) {

			const integer rowend = row_ptr[i77 + 1];
			const integer rowbeg = row_ptr[i77];

			integer jstart = 0;
			for (integer j77 = rowbeg; j77 < rowend; ++j77)
			{
				const integer element_offset = i77 + jstart * n;
				data[element_offset] = val[j77];
				indices[element_offset] = col_ind[j77];

				++jstart;
			}
			for (integer j77 = jstart; j77 < string_size; ++j77)
			{
				const integer element_offset = i77 + j77 * n;

				if (indices[element_offset - 1] < n - 2) {
					indices[element_offset] = indices[element_offset - 1] + 1;
					data[element_offset] = 0.0;
				}
				else {
					indices[element_offset] = indices[element_offset - 1];
					data[element_offset] = 0.0;
				}

			}
		}

		doublereal* dev_data = nullptr;
		integer* dev_indices = nullptr;	

		doublereal* dev_V = nullptr;
		doublereal* dev_tmp = nullptr;

		// Allocate GPU buffers for three vectors (two input, one output)    .
		cudaStatus = cudaMalloc((void**)&dev_data, (n * string_size) * sizeof(doublereal));
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMalloc dev_val failed!");
			system("PAUSE");
			exit(1);
		}

		// Allocate GPU buffers for three vectors (two input, one output)    .
		cudaStatus = cudaMalloc((void**)&dev_indices, (n * string_size) * sizeof(integer));
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMalloc dev_col_ind failed!");
			system("PAUSE");
			exit(1);
		}

		// Allocate GPU buffers for three vectors (two input, one output)    .
		cudaStatus = cudaMalloc((void**)&dev_V, (n) * sizeof(doublereal));
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMalloc dev_V failed!");
			system("PAUSE");
			exit(1);
		}

		// Allocate GPU buffers for three vectors (two input, one output)    .
		cudaStatus = cudaMalloc((void**)&dev_tmp, (n) * sizeof(doublereal));
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMalloc dev_tmp failed!");
			system("PAUSE");
			exit(1);
		}

		// Copy input vectors from host memory to GPU buffers.
		cudaStatus = cudaMemcpy(dev_data, data, n * string_size * sizeof(doublereal), cudaMemcpyHostToDevice);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMemcpy dev_val, val HostToDevice failed!");
			system("PAUSE");
			exit(1);
		}

		// Copy input vectors from host memory to GPU buffers.
		cudaStatus = cudaMemcpy(dev_indices, indices, n * string_size * sizeof(integer), cudaMemcpyHostToDevice);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMemcpy dev_col_ind, val HostToDevice failed!");
			system("PAUSE");
			exit(1);
		}

		delete[] data;
		data = nullptr;
		delete[] indices;
		indices = nullptr;

		// Copy input vectors from host memory to GPU buffers.
		cudaStatus = cudaMemcpy(dev_V, V, n * sizeof(doublereal), cudaMemcpyHostToDevice);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMemcpy dev_V, val HostToDevice failed!");
			system("PAUSE");
			exit(1);
		}

		MatrixEllPackItpackByVectorKernel<<<128, 128>>>(dev_data, dev_indices, string_size, dev_V, dev_tmp, n, 0);
		// Check for any errors launching the kernel
		cudaStatus = cudaGetLastError();
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, " MatrixEllPackItpackByVectorKernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
			system("PAUSE");
			exit(1);
		}

		cudaStatus = cudaMemcpy(tmp, dev_tmp, n * sizeof(doublereal), cudaMemcpyDeviceToHost);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMemcpy tmp, dev_tmp DeviceToHost failed!");
			system("PAUSE");
			exit(1);
		}


		cudaFree(dev_data);
		cudaFree(dev_indices);

		cudaFree(dev_V);
		cudaFree(dev_tmp);

#else

		doublereal* dev_val = nullptr;
		integer* dev_col_ind = nullptr;
		integer* dev_row_ptr = nullptr;
		doublereal* dev_V = nullptr;
		doublereal* dev_tmp = nullptr;

		integer nnz = row_ptr[n];

		// Allocate GPU buffers for three vectors (two input, one output)    .
		cudaStatus = cudaMalloc((void**)&dev_val, (nnz) * sizeof(doublereal));
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMalloc dev_val failed!");
			system("PAUSE");
			exit(1);
		}

		// Allocate GPU buffers for three vectors (two input, one output)    .
		cudaStatus = cudaMalloc((void**)&dev_col_ind, (nnz) * sizeof(integer));
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMalloc dev_col_ind failed!");
			system("PAUSE");
			exit(1);
		}

		// Allocate GPU buffers for three vectors (two input, one output)    .
		cudaStatus = cudaMalloc((void**)&dev_row_ptr, (n+1) * sizeof(integer));
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMalloc dev_row_ptr failed!");
			system("PAUSE");
			exit(1);
		}

		// Allocate GPU buffers for three vectors (two input, one output)    .
		cudaStatus = cudaMalloc((void**)&dev_V, (n) * sizeof(doublereal));
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMalloc dev_V failed!");
			system("PAUSE");
			exit(1);
		}

		// Allocate GPU buffers for three vectors (two input, one output)    .
		cudaStatus = cudaMalloc((void**)&dev_tmp, (n) * sizeof(doublereal));
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMalloc dev_tmp failed!");
			system("PAUSE");
			exit(1);
		}

		// Copy input vectors from host memory to GPU buffers.
		cudaStatus = cudaMemcpy(dev_val, val, nnz * sizeof(doublereal), cudaMemcpyHostToDevice);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMemcpy dev_val, val HostToDevice failed!");
			system("PAUSE");
			exit(1);
		}

		// Copy input vectors from host memory to GPU buffers.
		cudaStatus = cudaMemcpy(dev_col_ind, col_ind, nnz * sizeof(integer), cudaMemcpyHostToDevice);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMemcpy dev_col_ind, val HostToDevice failed!");
			system("PAUSE");
			exit(1);
		}


		// Copy input vectors from host memory to GPU buffers.
		cudaStatus = cudaMemcpy(dev_row_ptr, row_ptr, (n+1) * sizeof(integer), cudaMemcpyHostToDevice);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMemcpy dev_rho_ptr, val HostToDevice failed!");
			system("PAUSE");
			exit(1);
		}

		// Copy input vectors from host memory to GPU buffers.
		cudaStatus = cudaMemcpy(dev_V, V, n * sizeof(doublereal), cudaMemcpyHostToDevice);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMemcpy dev_V, val HostToDevice failed!");
			system("PAUSE");
			exit(1);
		}

		// Copy input vectors from host memory to GPU buffers.
		cudaStatus = cudaMemcpy(dev_tmp, tmp, n * sizeof(doublereal), cudaMemcpyHostToDevice);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMemcpy dev_tmp, val HostToDevice failed!");
			system("PAUSE");
			exit(1);
		}

		MatrixCRSByVectorKernel<<<128,128>>>(dev_val, dev_col_ind, dev_row_ptr,dev_V, dev_tmp,n);
		// Check for any errors launching the kernel
		cudaStatus = cudaGetLastError();
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "MatrixCRSByVectorKernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
			system("PAUSE");
			exit(1);
		}

		cudaStatus = cudaMemcpy(tmp, dev_tmp, n * sizeof(doublereal), cudaMemcpyDeviceToHost);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMemcpy tmp, dev_tmp DeviceToHost failed!");
			system("PAUSE");
			exit(1);
		}

		cudaFree(dev_tmp);
		cudaFree(dev_val);
		cudaFree(dev_col_ind);
		cudaFree(dev_V);
		cudaFree(dev_row_ptr);

#endif

#endif

	}
	else {
		// на процессоре CPU.

		if (number_cores() == 1) {

			// Если у нас только один поток то лучше вообще не писать #pragma omp parallel for
			// т.к. с #pragma omp parallel for будет медленней

			for (integer i = 0; i < n; ++i) {
				doublereal sum = 0.0;
				const integer rowend = row_ptr[i + 1];
				const integer rowbeg = row_ptr[i];
				
				for (integer j = rowbeg; j < rowend; ++j)
				{
					sum += val[j] * V[col_ind[j]];
				}
				tmp[i] = sum;
			}
		}
		else {

			// вектор tmp индексируется начиная с нуля так же как и вектор V
	//#pragma omp parallel for
		//	for (integer i = 0; i < n; ++i) tmp[i] = 0.0;

			// В целях увеличения быстродействия
			// вся необходимая память выделяется заранее.
			//if (tmp == nullptr)
			//{
			//printf("malloc: out of memory for vector tmp in MatrixCRSByVector\n"); // нехватка памяти
			//system("pause");
			//exit(0);  // завершение программы
			//}



			//omp_set_num_threads(inumcore);

#pragma omp parallel for  schedule (guided)
			for (integer i = 0; i < n; ++i) {
				doublereal sum = 0.0;
				const integer rowend = row_ptr[i + 1];
				const integer rowbeg = row_ptr[i];

				
				for (integer j = rowbeg; j < rowend; ++j)
				{
					sum += val[j] * V[col_ind[j]];
				}
				tmp[i] = sum;
			}

			//return tmp;
		}
	}
} // MatrixCRSByVector


#if MY_GPU_MATH_LIBRARY_CU_ON

__global__ void Jacoby_kernel_internal(const doublereal* data, const int* col_ind, const doublereal* b,
	doublereal* x_new, const doublereal* x, doublereal* residual, integer n)
{

	int i = threadIdx.x + blockIdx.x * blockDim.x;
	while (i < n) {

		const int is = i * 7;
		const int is1 = is + 1;
		const int is2 = is + 2;
		const int is3 = is + 3;
		const int is4 = is + 4;
		const int is5 = is + 5;
		const int is6 = is + 6;
		const int it6 = col_ind[is6];

		const double ptildaQ = (data[is] * x[col_ind[is]] +
			data[is1] * x[col_ind[is1]] +
			data[is2] * x[col_ind[is2]] +
			data[is3] * x[col_ind[is3]] +
			data[is4] * x[col_ind[is4]] +
			data[is5] * x[col_ind[is5]] +
			b[it6]) / data[is6];

		const doublereal ptilda = static_cast<doublereal>(ptildaQ);
		residual[it6] = fabs((ptilda - x[it6])); // сумма модулей
		x_new[it6] = x[it6] + (ptilda - x[it6]);

		i += blockDim.x * gridDim.x;
	}
}

__global__ void Jacoby_kernel_bound(const doublereal* data_b, const int* col_ind_b, const doublereal* b,
	doublereal* x_new, const doublereal* x, doublereal* residual, integer n)
{

	int i = threadIdx.x + blockIdx.x * blockDim.x;
	while (i < n) {

		const int is = i * 2;
		const int it1 = col_ind_b[is];//iI
		const int it2 = col_ind_b[is + 1];//iW 

		const doublereal ptilda = (data_b[is] * x[it1] + b[it2]) / data_b[is + 1];
		residual[it2] = fabs((ptilda - x[it2])); // сумма модулей
		if (it1 == it2) x_new[it2] = (ptilda);
		else x_new[it2] = x[it2] + (ptilda - x[it2]);

		i += blockDim.x * gridDim.x;
	}
}

__global__ void Chebyshev_kernel_bound(const doublereal* data_b, const int* col_ind_b, const doublereal* b,
	doublereal* x_new, const doublereal* x,  doublereal rURF, integer n)
{
	int i = threadIdx.x + blockIdx.x * blockDim.x;
	while (i < n) {

		const int is = i * 2;
		const int it1 = col_ind_b[is];//iI
		const int it2 = col_ind_b[is + 1];//iW                   


		const doublereal residual1 = (data_b[is] * x[it1] + b[it2]) - data_b[is + 1] * x[it2];
		// Чебышева
		
	

		if (it1 == it2) x_new[it2] = (b[it2] / data_b[is + 1]);
		else x_new[it2] = x[it2] + rURF * (residual1);

		i += blockDim.x * gridDim.x;
	}

}

__global__ void Chebyshev_kernel_bound2(const doublereal* data_b, const int* col_ind_b, const doublereal* b,
	doublereal* x_new, const doublereal* x, doublereal* residual, doublereal rURF, integer n)
{
	int i = threadIdx.x + blockIdx.x * blockDim.x;
	while (i < n) {

		const int is = i * 2;
		const int it1 = col_ind_b[is];//iI
		const int it2 = col_ind_b[is + 1];//iW                   


		const doublereal residual1 = (data_b[is] * x[it1] + b[it2]) - data_b[is + 1] * x[it2];
		// Чебышева
		residual[it2] = fabs(residual1); // сумма модулей.


		if (it1 == it2) x_new[it2] = (b[it2] / data_b[is + 1]);
		else x_new[it2] = x[it2] + rURF * (residual1);

		i += blockDim.x * gridDim.x;
	}

}


__global__ void Chebyshev_kernel_bound_residual2(const doublereal* data_b, const int* col_ind_b, const doublereal* b,
	 const doublereal* x, doublereal* residual,  integer n)
{
	int i = threadIdx.x + blockIdx.x * blockDim.x;
	while (i < n) {

		const int is = i * 2;
		const int it1 = col_ind_b[is];//iI
		const int it2 = col_ind_b[is + 1];//iW                   


		const doublereal residual1 = (data_b[is] * x[it1] + b[it2]) - data_b[is + 1] * x[it2];
		// Чебышева
		residual[it2] = (residual1); // сумма модулей.


		if (it1 == it2) residual[it2] = 0.0;
		

		i += blockDim.x * gridDim.x;
	}

}


__global__ void  Chebyshev_kernel_internal(const doublereal* data, const int* col_ind, const doublereal* b,
	doublereal* x_new, const doublereal* x,  const doublereal  rURF, const integer n)
{
	int i = threadIdx.x + blockIdx.x * blockDim.x;
	while (i < n) {

		const int is = i * 7;
		const int is1 = is + 1;
		const int is2 = is + 2;
		const int is3 = is + 3;
		const int is4 = is + 4;
		const int is5 = is + 5;
		const int is6 = is + 6;
		const int it6 = col_ind[is6];

		const double residualQ = (data[is] * x[col_ind[is]] +
			data[is1] * x[col_ind[is1]] +
			data[is2] * x[col_ind[is2]] +
			data[is3] * x[col_ind[is3]] +
			data[is4] * x[col_ind[is4]] +
			data[is5] * x[col_ind[is5]] +
			b[it6]) - data[is6] * x[it6];

		const doublereal residual = static_cast<doublereal>(residualQ);


		x_new[it6] = x[it6] + rURF * (residual);

		i += blockDim.x * gridDim.x;
	}
}


__global__ void  Chebyshev_kernel_internal2(const doublereal* data, const int* col_ind, const doublereal* b,
	doublereal* x_new, const doublereal* x,  doublereal* residual, const doublereal  rURF, const integer n)
{
	int i = threadIdx.x + blockIdx.x * blockDim.x;
	while (i < n) {

		const int is = i * 7;
		const int is1 = is + 1;
		const int is2 = is + 2;
		const int is3 = is + 3;
		const int is4 = is + 4;
		const int is5 = is + 5;
		const int is6 = is + 6;
		const int it6 = col_ind[is6];

		const double residualQ = (data[is] * x[col_ind[is]] +
			data[is1] * x[col_ind[is1]] +
			data[is2] * x[col_ind[is2]] +
			data[is3] * x[col_ind[is3]] +
			data[is4] * x[col_ind[is4]] +
			data[is5] * x[col_ind[is5]] +
			b[it6]) - data[is6] * x[it6];

		const doublereal residual1 = static_cast<doublereal>(residualQ);

		residual[it6] = fabs(residual1);

		x_new[it6] = x[it6] + rURF * (residual1);

		i += blockDim.x * gridDim.x;
	}
}

__global__ void  Chebyshev_kernel_internal_residual2(const doublereal* data, const int* col_ind, const doublereal* b,
	 const doublereal* x, doublereal* residual,  const integer n)
{
	int i = threadIdx.x + blockIdx.x * blockDim.x;
	while (i < n) {

		const int is = i * 7;
		const int is1 = is + 1;
		const int is2 = is + 2;
		const int is3 = is + 3;
		const int is4 = is + 4;
		const int is5 = is + 5;
		const int is6 = is + 6;
		const int it6 = col_ind[is6];

		const double residualQ = (data[is] * x[col_ind[is]] +
			data[is1] * x[col_ind[is1]] +
			data[is2] * x[col_ind[is2]] +
			data[is3] * x[col_ind[is3]] +
			data[is4] * x[col_ind[is4]] +
			data[is5] * x[col_ind[is5]] +
			b[it6]) - data[is6] * x[it6];

		const doublereal residual1 = static_cast<doublereal>(residualQ);

		residual[it6] = (residual1);

		

		i += blockDim.x * gridDim.x;
	}
}

__global__ void copy_Jacoby_kernel(const doublereal* x_new, doublereal* x, integer n)
{
	int i = threadIdx.x + blockIdx.x * blockDim.x;
	while (i < n) {

		x[i] = x_new[i];

		i += blockDim.x * gridDim.x;
	}
}

__global__ void sqr_vec_kernel(doublereal* x, integer n)
{
	int i = threadIdx.x + blockIdx.x * blockDim.x;
	while (i < n) {

		x[i] *= x[i];

		i += blockDim.x * gridDim.x;
	}
}

const int  threadsPerBlock = 128;
#define blocksPerGrid(size0) (static_cast<unsigned int>(32 < ((size0)+threadsPerBlock-1)/threadsPerBlock ? 32 : ((size0)+threadsPerBlock-1)/threadsPerBlock))

__global__ void reduction_sum_kernel(doublereal* x, integer n, doublereal* sum)
{
	__shared__ doublereal cache[threadsPerBlock];
	int i = threadIdx.x + blockIdx.x * blockDim.x;
	int cacheIndex = threadIdx.x;
	doublereal temp = 0.0;
	while (i < n) {

		temp += x[i];

		i += blockDim.x * gridDim.x;
	}

	// сохранить значения в кэше
	cache[cacheIndex] = temp;
	//синхронизировать нити в этом блоке
	__syncthreads();

	// для успешной редукции threadsPerBlock должно быть степенью 2
	// из за следующего кода
	int j = blockDim.x / 2;
	while (j != 0) {
		if (cacheIndex < j)
			cache[cacheIndex] += cache[cacheIndex + j];
		__syncthreads();
		j /= 2;
	}

	if (cacheIndex == 0)
		sum[blockIdx.x] = cache[0];
}

#endif

// применяется для компонент скорости
// Достаточно относительной точности 0.1.
// Полностью параллельный код.
// Метод Якоби без использования параметра нижней релаксации.
// Судя по статьям программного комплекса ЛОГОС для компонент 
// скорости может использоваться метод Якоби т.к. число 
// обусловленности матриц на компоненты скорости примерно 1.0e+1;
// Для сравнения число обусловленности для давления 1e+7 - 1e+8.
// 
// Работает только для структурированной расчётной сетки.
// 
// Работает на видеокарте. 19.02.2022
// 23.02.2022 scal(.,.),  reduction  на cuda
//
void specialSpeed_cuda(equation3D*& sl, equation3D_bon*& slb, doublereal*& b,
	doublereal*& x, integer maxelm, integer maxbound, integer iVar,
	bool bprintMessage, doublereal alpharelax) {
	// К величине xcor - будет осуществляться нижняя релаксация, т.е. x_cor - это
	// скорректированная компонента скорости удовлетворяющая уравнению неразрывности.
	//printf("SOR3D incomming...\n"); // debug.
	//system("pause");

	// Параметр релаксации лежит в интервале от 0.0 до 2.0.
	// Верхняя релаксация способна существенно ускорить вычислительный процесс. 
	// Патанкар рекомендует брать коэффициент верхней релаксации равный 1.5;
	// rURF=1.5 можно использовать в уравнении на поправку давления.
	// В уравнениях на скорость мы имеем двоякую релаксацию. С одной стороны это нижняя релаксация
	// к скорректированной скорости с коэффициентом alpha. Эта нижняя релаксация учитывается в матрице СЛАУ,
	// в результате чего мы получаем новую матрицу СЛАУ.  К этой новой преоблразованной матрице СЛАУ казалось бы в целях ускорения сходимости
	// можно применить верхнюю релаксацию, пусть она будет опять с коэффициентом rURF=1.5. Вычисления показывают что при введении в уравнения на скорость
	// коэффициента верхней релаксации 1.5 точность решения уравнений на скорость за 4000 итераций падает (или вообще имеем расходимость),
	// поэтому наверное лучше не вводить коээффициент 
	// верхней релаксации равный 1.5 в уравнениях на скорость. 

#ifdef _OPENMP 
	omp_set_num_threads(number_cores()); // установка числа потоков
#endif

	// пороговое значение невязки
	doublereal eps = 1e-40;
	//doublereal sRF = 1.0;// alpharelax;


	integer j = 0, kend = 10000000;//100; // Для целей пост сглаживания должно хватить 40 итераций.

	const integer size0 = maxelm + maxbound;

	doublereal* x_new = new doublereal[size0];
	doublereal* residual = new doublereal[size0];
#pragma omp parallel for
	for (int i = 0; i < size0; ++i)
	{
		x_new[i] = 0.0;
		residual[i] = 0.0;
	}

	// Улучшает попадание в кеш ускорение на 12.5%.
	doublereal* data = new doublereal[maxelm * 7];
	int* col_ind = new int[maxelm * 7];

#pragma omp parallel for 
	for (int i = 0; i < maxelm; ++i) {

		const int is = i * 7;
		data[is] = sl[i].ae;
		data[is + 1] = sl[i].aw;
		data[is + 2] = sl[i].an;
		data[is + 3] = sl[i].as;
		data[is + 4] = sl[i].at;
		data[is + 5] = sl[i].ab;
		data[is + 6] = sl[i].ap / alpharelax;

		col_ind[is] = sl[i].iE;
		col_ind[is + 1] = sl[i].iW;
		col_ind[is + 2] = sl[i].iN;
		col_ind[is + 3] = sl[i].iS;
		col_ind[is + 4] = sl[i].iT;
		col_ind[is + 5] = sl[i].iB;
		col_ind[is + 6] = sl[i].iP;
	}


	doublereal* data_b = new doublereal[maxbound * 2];
	int* col_ind_b = new int[maxbound * 2];

#pragma omp parallel for 
	for (int i = 0; i < maxbound; ++i) {

		const int is = i * 2;
		data_b[is] = slb[i].ai;
		data_b[is + 1] = slb[i].aw;

		if (slb[i].iI > -1) {
			col_ind_b[is] = slb[i].iI;
		}
		else {
			data_b[is] = 0.0;
			col_ind_b[is] = slb[i].iW;
		}
		col_ind_b[is + 1] = slb[i].iW;
	}

#if MY_GPU_MATH_LIBRARY_CU_ON

	if (init_b_first_device) {

		int device;
		device = idevice_Tesla;

		// Choose which GPU to run on, change this on a multi-GPU system.
		cudaStatus = cudaSetDevice(device);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?");
			system("PAUSE");
			exit(1);
		}

		cudaDeviceProp prop;
		cudaGetDeviceProperties(&prop, device);
		printf("%s\n", prop.name);
		// GeForce 840M 384потока 1ГГц каждый март 2014 года. 28нм.

		init_b_first_device = false;
	}

	doublereal* dev_data = nullptr;
	int* dev_col_ind = nullptr;

	doublereal* dev_data_b = nullptr;
	int* dev_col_ind_b = nullptr;


	doublereal* dev_x_new = nullptr;
	doublereal* dev_x = nullptr;
	doublereal* dev_residual = nullptr;
	doublereal* dev_b = nullptr;

	doublereal* sum = (doublereal*)malloc(blocksPerGrid(size0) * sizeof(doublereal));
	doublereal* dev_sum;

	// Allocate GPU buffers for three vectors (two input, one output)    .
	cudaStatus = cudaMalloc((void**)&dev_sum, blocksPerGrid(size0) * sizeof(doublereal));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc dev_sum failed!");
		system("PAUSE");
		exit(1);
	}

	// Allocate GPU buffers for three vectors (two input, one output)    .
	cudaStatus = cudaMalloc((void**)&dev_data, (maxelm * 7) * sizeof(doublereal));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc dev_data failed!");
		system("PAUSE");
		exit(1);
	}

	// Allocate GPU buffers for three vectors (two input, one output)    .
	cudaStatus = cudaMalloc((void**)&dev_col_ind, (maxelm * 7) * sizeof(int));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc dev_col_ind failed!");
		system("PAUSE");
		exit(1);
	}

	// Allocate GPU buffers for three vectors (two input, one output)    .
	cudaStatus = cudaMalloc((void**)&dev_data_b, (maxbound * 2) * sizeof(doublereal));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc dev_data_b failed!");
		system("PAUSE");
		exit(1);
	}

	// Allocate GPU buffers for three vectors (two input, one output)    .
	cudaStatus = cudaMalloc((void**)&dev_col_ind_b, (maxbound * 2) * sizeof(int));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc dev_col_ind_b failed!");
		system("PAUSE");
		exit(1);
	}

	// Allocate GPU buffers for three vectors (two input, one output)    .
	cudaStatus = cudaMalloc((void**)&dev_x_new, (size0) * sizeof(doublereal));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc dev_xnew failed!");
		system("PAUSE");
		exit(1);
	}

	// Allocate GPU buffers for three vectors (two input, one output)    .
	cudaStatus = cudaMalloc((void**)&dev_b, (size0) * sizeof(doublereal));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc dev_b failed!");
		system("PAUSE");
		exit(1);
	}

	// Allocate GPU buffers for three vectors (two input, one output)    .
	cudaStatus = cudaMalloc((void**)&dev_x, (size0) * sizeof(doublereal));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc dev_x failed!");
		system("PAUSE");
		exit(1);
	}

	// Allocate GPU buffers for three vectors (two input, one output)    .
	cudaStatus = cudaMalloc((void**)&dev_residual, (size0) * sizeof(doublereal));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc dev_residual failed!");
		system("PAUSE");
		exit(1);
	}

	// Copy input vectors from host memory to GPU buffers.
	cudaStatus = cudaMemcpy(dev_data, data, maxelm * 7 * sizeof(doublereal), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy dev_data, val HostToDevice failed!");
		system("PAUSE");
		exit(1);
	}

	// Copy input vectors from host memory to GPU buffers.
	cudaStatus = cudaMemcpy(dev_col_ind, col_ind, maxelm * 7 * sizeof(int), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy dev_col_ind, val HostToDevice failed!");
		system("PAUSE");
		exit(1);
	}

	// Copy input vectors from host memory to GPU buffers.
	cudaStatus = cudaMemcpy(dev_data_b, data_b, maxbound * 2 * sizeof(doublereal), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy dev_data_b, val HostToDevice failed!");
		system("PAUSE");
		exit(1);
	}

	// Copy input vectors from host memory to GPU buffers.
	cudaStatus = cudaMemcpy(dev_col_ind_b, col_ind_b, maxbound * 2 * sizeof(int), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy dev_col_ind_b, val HostToDevice failed!");
		system("PAUSE");
		exit(1);
	}

	// Copy input vectors from host memory to GPU buffers.
	cudaStatus = cudaMemcpy(dev_x_new, x_new, size0 * sizeof(doublereal), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy dev_xnew, val HostToDevice failed!");
		system("PAUSE");
		exit(1);
	}


	// Copy input vectors from host memory to GPU buffers.
	cudaStatus = cudaMemcpy(dev_b, b, size0 * sizeof(doublereal), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy dev_b val HostToDevice failed!");
		system("PAUSE");
		exit(1);
	}


	// Copy input vectors from host memory to GPU buffers.
	cudaStatus = cudaMemcpy(dev_x, x, size0 * sizeof(doublereal), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy dev_x, val HostToDevice failed!");
		system("PAUSE");
		exit(1);
	}

	// Copy input vectors from host memory to GPU buffers.
	cudaStatus = cudaMemcpy(dev_residual, residual, size0 * sizeof(doublereal), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy dev_residual, val HostToDevice failed!");
		system("PAUSE");
		exit(1);
	}

	doublereal dmax = 1.0;
	while ((dmax > eps) && (j < kend)) {
		dmax = 0.0;

		doublereal dmax_loc = 0.0;


		Jacoby_kernel_internal << <128, 128 >> > (dev_data, dev_col_ind, dev_b,
			dev_x_new, dev_x, dev_residual, maxelm);

		// Check for any errors launching the kernel
		cudaStatus = cudaGetLastError();
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "Jacoby_kernel_internal launch failed: %s\n", cudaGetErrorString(cudaStatus));
			system("PAUSE");
			exit(1);
		}

		//cudaStatus = cudaMemcpy(residual, dev_residual, size0 * sizeof(doublereal), cudaMemcpyDeviceToHost);
		//if (cudaStatus != cudaSuccess) {
			//fprintf(stderr, "cudaMemcpy residual, dev_residual DeviceToHost failed!");
			//system("PAUSE");
			//exit(1);
		//}

//#pragma omp parallel for  reduction(+: dmax_loc)
	//	for (int i = 0; i < maxelm; ++i)
		//{
			//const int is = i * 7;
			//const int is6 = is + 6;
			//const int it6 = col_ind[is6];

			//dmax_loc += residual[it6];
		//}

		//dmax += dmax_loc;
		//dmax_loc = 0.0;


		Jacoby_kernel_bound << <128, 128 >> > (dev_data_b, dev_col_ind_b, dev_b,
			dev_x_new, dev_x, dev_residual, maxbound);

		// Check for any errors launching the kernel
		cudaStatus = cudaGetLastError();
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "Jacoby_kernel_bound launch failed: %s\n", cudaGetErrorString(cudaStatus));
			system("PAUSE");
			exit(1);
		}

		//cudaStatus = cudaMemcpy(residual, dev_residual, size0 * sizeof(doublereal), cudaMemcpyDeviceToHost);
		//if (cudaStatus != cudaSuccess) {
			//fprintf(stderr, "cudaMemcpy residual, dev_residual DeviceToHost failed!");
			//system("PAUSE");
			//exit(1);
		//}

//#pragma omp parallel for reduction(+: dmax_loc)
	//	for (int i = 0; i < maxbound; ++i) {

		//	const int is = i * 2;
			//const int it2 = col_ind_b[is + 1];//iW 

//			dmax_loc += residual[it2];
	//	}

		//dmax_loc = thrust::reduce(residual, residual + size0);

		reduction_sum_kernel << <blocksPerGrid(size0), threadsPerBlock >> > (dev_residual, size0, dev_sum);

		// Check for any errors launching the kernel
		cudaStatus = cudaGetLastError();
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "reduction_sum_kernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
			system("PAUSE");
			exit(1);
		}

		cudaStatus = cudaMemcpy(sum, dev_sum, blocksPerGrid(size0) * sizeof(doublereal), cudaMemcpyDeviceToHost);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMemcpy residual, dev_residual DeviceToHost failed!");
			system("PAUSE");
			exit(1);
		}

		dmax_loc = 0.0;
		const int isize_cycle = static_cast<int>(blocksPerGrid(size0));

#pragma omp parallel for reduction(+: dmax_loc)
		for (int i89 = 0; i89 < isize_cycle; ++i89) {
			dmax_loc += sum[i89];
		}

		dmax += dmax_loc;

		if (bprintMessage)
		{
			if (iVar == PAM)
			{
				//dmax/=maxelm;
				if (j % 1 == 0) {
#if doubleintprecision == 1
					printf("%lld %e \n", j + 1, dmax);
#else
					printf("%d %e \n", j + 1, dmax);
#endif

				}
			}
		}

		if (j == 0) eps = 0.1 * dmax; // Важнейшее условие выхода из итерационного процесса.

		copy_Jacoby_kernel << <128, 128 >> > (dev_x_new, dev_x, size0);

		// Check for any errors launching the kernel
		cudaStatus = cudaGetLastError();
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "copy_Jacoby_kernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
			system("PAUSE");
			exit(1);
		}

		++j;

	}

	cudaStatus = cudaMemcpy(x, dev_x, size0 * sizeof(doublereal), cudaMemcpyDeviceToHost);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy x, dev_x DeviceToHost failed!");
		system("PAUSE");
		exit(1);
	}

	cudaFree(dev_data);
	cudaFree(dev_col_ind);

	cudaFree(dev_data_b);
	cudaFree(dev_col_ind_b);

	cudaFree(dev_b);
	cudaFree(dev_x_new);
	cudaFree(dev_x);
	cudaFree(dev_residual);

	cudaFree(dev_sum);
	free(sum);

#else

	doublereal dmax = 1.0;
	while ((dmax > eps) && (j < kend)) {
		dmax = 0.0;

		doublereal dmax_loc = 0.0;

#pragma omp parallel for  reduction(+: dmax_loc)
		for (int i = 0; i < maxelm; ++i)
		{

			const int is = i * 7;
			const int is1 = is + 1;
			const int is2 = is + 2;
			const int is3 = is + 3;
			const int is4 = is + 4;
			const int is5 = is + 5;
			const int is6 = is + 6;
			const int it6 = col_ind[is6];

			const long double ptildaQ = (data[is] * x[col_ind[is]] +
				data[is1] * x[col_ind[is1]] +
				data[is2] * x[col_ind[is2]] +
				data[is3] * x[col_ind[is3]] +
				data[is4] * x[col_ind[is4]] +
				data[is5] * x[col_ind[is5]] +
				b[it6]) / data[is6];

			const doublereal ptilda = static_cast<doublereal>(ptildaQ);
			dmax_loc += fabs((ptilda - x[it6])); // сумма модулей
			x_new[it6] = x[it6] + (ptilda - x[it6]);


			/*
			if (0&&j==0) {// debug
				printf("ae=%e, aw=%e, an=%e, as=%e, at=%e, ab=%e, ap=%e, b=%e\n",sl[i].ae,sl[i].aw,sl[i].an,sl[i].as,sl[i].at,sl[i].ab,sl[i].ap,sl[i].b);
				system("pause");
			}
			*/
		}

		dmax += dmax_loc;
		dmax_loc = 0.0;

#pragma omp parallel for reduction(+: dmax_loc)
		for (int i = 0; i < maxbound; ++i) {

			const int is = i * 2;
			const int it1 = col_ind_b[is];//iI
			const int it2 = col_ind_b[is + 1];//iW 

			const doublereal ptilda = (data_b[is] * x[it1] + b[it2]) / data_b[is + 1];
			dmax_loc += fabs((ptilda - x[it2])); // сумма модулей
			if (it1 == it2) x_new[it2] = (ptilda);
			else x_new[it2] = x[it2] +  (ptilda - x[it2]);

		}

		dmax += dmax_loc;

		if (bprintMessage)
		{
			if (iVar == PAM)
			{
				//dmax/=maxelm;
				//if (j % 1 == 0) Печатаем всегда. Раскоментировать и внести правки если нужно печатать реже.
				{
#if doubleintprecision == 1
					printf("%lld %e \n", j + 1, dmax);
#else
					printf("%d %e \n", j + 1, dmax);
#endif

				}
			}
		}

		if (j == 0) eps = 0.1 * dmax; // Важнейшее условие выхода из итерационного процесса.

#pragma omp parallel for
		for (integer i = 0; i < size0; ++i)
		{
			x[i] = x_new[i];
		}


		j++;
	}

#endif

	delete[] data;
	delete[] col_ind;

	delete[] data_b;
	delete[] col_ind_b;

	delete[] x_new;
	delete[] residual;

	if (bprintMessage) {
		if (iVar == PAM) {
			printf("calc complete...\n");
			system("pause");
		}
	}
	//printf("4000 %e \n", dmax);
	//system("pause");

} // specialSpeed_cuda


// 08.01.2022
typedef struct TLEVEL_CHEBYSHEV_INFO {

	// 11.01.2022
	// Для метода Чебышева нужен положительно определенный SPD оператор.
	// Т.е. для задач с конвекцией метод Чебышева неприменим. 
	// Метод может быть использован для теплопередачи в твёрдом теле, для механики,
	// для поправки давления при решении уравнений гидродинамики.
	// Не стоит прменять метод Чебышева для компонент скорости, уравнений переноса
	// турбулентных характеристик, конвективной теплопередачи.



	// Сначала должно быть инициализировано значением false чтобы применить теорему Гершгорина.
	bool bGershgorin_calc_Temp, bGershgorin_calc_Speed, bGershgorin_calc_Pressure, bGershgorin_calc_Stress;
	bool bGershgorin_calc_Vx, bGershgorin_calc_Vy, bGershgorin_calc_Vz;
	bool bGershgorin_calc_nu, bGershgorin_calc_k, bGershgorin_calc_omega;
	bool bGershgorin_calc_gamma, bGershgorin_calc_ReTheta, bGershgorin_calc_epsilon;
	bool bGershgorin_calc_k_for_ke;

	doublereal hi_Temp, hi_Speed, hi_Pressure, hi_Stress; // храним спектральный радиус.
	doublereal hi_Vx, hi_Vy, hi_Vz;
	doublereal hi_nu, hi_k, hi_omega;
	doublereal hi_gamma, hi_ReTheta, hi_epsilon;
	doublereal hi_k_for_ke;

	// Сначала должно быть инициализировано значением false чтобы применить теорему Релея - Ритца.
	bool bstart_lo_Temp, bstart_lo_Speed, bstart_lo_Pressure, bstart_lo_Stress;
	bool bstart_lo_Vx, bstart_lo_Vy, bstart_lo_Vz;
	bool bstart_lo_nu, bstart_lo_k, bstart_lo_omega;
	bool bstart_lo_gamma, bstart_lo_ReTheta, bstart_lo_epsilon;
	bool bstart_lo_k_for_ke;

	doublereal lo_Temp, lo_Speed, lo_Pressure, lo_Stress; // храним текущую нижнию границу спектра. Она будет уточнятся в ходе вычислений.
	doublereal lo_Vx, lo_Vy, lo_Vz;
	doublereal lo_nu, lo_k, lo_omega;
	doublereal lo_gamma, lo_ReTheta, lo_epsilon;
	doublereal lo_k_for_ke;

} LEVEL_CHEBYSHEV_INFO;

// В многосеточном методе сперктральный радиус и нижняя граница спектра свои 
// для каждого уровня многосеточной иерархии.
LEVEL_CHEBYSHEV_INFO* level_Chebyshev_info = nullptr;



// Estimate spectral radius of the matrix.
// Use Gershgorin disk theorem when power_iters == 0,
// Use Power method when power_iters > 0.
// When scale = true, scale the matrix by its inverse diagonal.

doublereal spectral_radius(equation3D*& sl, equation3D_bon*& slb,
	integer maxelm, integer maxbound)
{
	// Теорема Гершгорина о кругах.

	doublereal radius;

	//if (power_iters <= 0)
	{
		// Use Gershgorin disk theorem.
		radius = 0;

#pragma omp parallel
		{
			doublereal emax = 0;
			// doublereal  dia = 1.0;



#pragma omp for nowait
			for (integer i = 0; i < maxelm; ++i) {

				doublereal s = 0;

				s += sl[i].ae + sl[i].aw + sl[i].an + sl[i].as + sl[i].at + sl[i].ab + sl[i].ap;
				//s += fabs(sl[i].ae) + fabs(sl[i].aw) + fabs(sl[i].an) + fabs(sl[i].as) + fabs(sl[i].at) + fabs(sl[i].ab) + fabs(sl[i].ap);
				//dia = sl[i].ap;

				emax = std::max(emax, s);
			}

#pragma omp for nowait
			for (integer i = 0; i < maxbound; ++i) {

				doublereal s = 0;

				if (slb[i].iI == -1) {
					// s += fabs(slb[i].aw);
					s += slb[i].aw;
				}
				else {
					//s += fabs(slb[i].aw) + fabs(slb[i].ai);
					s += slb[i].aw + slb[i].ai;
				}
				//dia = slb[i].aw;

				emax = std::max(emax, s);
			}

#pragma omp critical
			radius = std::max(radius, emax);
		}
	}
	//std::cout << "spectral radius" << std::endl;

	return  (radius < 0 ? static_cast<doublereal>(2.0) : radius);
}

// Функция переупорядочивания нулей многочлена Чебышева.
void Lebedev_Samarskii_Nikolaev(integer degree, doublereal*& mu)
{

	// Без переупорядочивания.
	//int theta4[4] = { 1,3,5,7 };
	//int theta5[5] = { 1,3,5,7,9 };
	//int theta6[6] = { 1,3,5,7,9,11 };
	//int theta6[7] = { 1,3,5,7,9,11,13 };
	//int theta8[8] = { 1,3,5,7,9,11,13,15 };
	//int theta16[16] = { 1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31 };

	// С переупорядочиванием корней многочлена Чебышева.
	// 29.01.2022
	// int theta4[4] = { 1,7,3,5 };
	int theta5[5] = { 1, 9, 3, 7, 5 };
	// int theta6[6] = { 1,11,5,7,3,9 };
	// int theta7[7] = { 1,13,5,9,3,11,7 };
	int theta8[8] = { 1,15,7,9,3,13,5,11 };
	// int theta12[12] = {1, 23, 11, 13, 5, 19, 7, 17, 3, 21, 9, 15};
	// int theta14[14] = { 1,29,13,17,5,25,9,21,3,27,11,19,7,23 };
	// int theta15[15] = { 1,29,13,17,5,25,9,21,3,27,11,19,7,23, 15 };
	int theta16[16] = { 1,31, 15,17,7,25,9,23,3,29,13,19,5,27,11,21 };
	// int theta24[24] = { 1, 49, 23, 27, 11, 39, 13, 37, 5, 45, 19, 31, 7, 43, 17, 33, 3, 47, 21, 29, 9, 41, 15, 35};
	int theta25[25] = { 1, 49, 23, 27, 11, 39, 13, 37, 5, 45, 19, 31, 7, 43, 17, 33, 3, 47, 21, 29, 9, 41, 15, 35, 25 };
	int theta32[32] = { 1, 63, 31, 33, 15, 49, 17, 47, 7, 57, 25, 39, 9, 55, 23, 41, 3, 61, 29, 35, 13, 51, 19, 45, 5, 59, 27, 37, 11, 53, 21, 43 };
	int theta64[64] = { 1, 127, 63, 65, 31, 97, 33, 95, 15, 113, 49, 79, 17, 111, 47, 81, 7, 121, 57, 71, 25, 103, 39, 89, 9, 119, 55, 73, 23, 105, 41, 87, 3, 125, 61, 67, 29, 99, 35, 93, 13, 115, 51, 77, 19, 109, 45, 83, 5, 123, 59, 69, 27, 101, 37, 91, 11, 117, 53, 75, 21, 107, 43, 85 };
	int theta128[128] = { 1, 255, 127,129, 63,193, 65,191, 31,225, 97,159, 33,223, 95,161, 15,241, 113,143, 49,207, 79,177, 17,239, 111,145, 47,209, 81,175, 7,249, 121,135, 57,199, 71,185, 25,231, 103,153, 39,217, 89,167, 9,247, 119,137, 55,201, 73,183, 23,233, 105,151, 41,215, 87,169, 3,253, 125,131, 61,195, 67,189, 29,227, 99,157, 35,221, 93,163, 13,243, 115,141, 51,205, 77,179, 19,237, 109,147, 45,211, 83,173, 5,251, 123,133, 59,197, 69,187, 27,229, 101,155, 37,219, 91,165, 11,245, 117,139, 53,203, 75,181, 21,235, 107,149, 43,213, 85, 171 };

	// 06.02.2022
	int theta256[256];
	for (int k = 0; k < 128; ++k) {
		const int is = 2 * k;
		theta256[is] = theta128[k];
		theta256[is + 1] = 512 - theta256[is];
	}

	int theta512[512];
	for (int k = 0; k < 256; ++k) {
		const int is = 2 * k;
		theta512[is] = theta256[k];
		theta512[is + 1] = 1024 - theta512[is];
	}

	switch (degree) {
	case 5:
		for (int k = 0; k < degree; ++k) {
			mu[k] = -cos(M_PI * theta5[k] / (2.0 * degree));
		}
		break;
	case 8:
		for (int k = 0; k < degree; ++k) {
			mu[k] = -cos(M_PI * theta8[k] / (2.0 * degree));
		}
		break;
	case 16:
		for (int k = 0; k < degree; ++k) {
			mu[k] = -cos(M_PI * theta16[k] / (2.0 * degree));
		}
		break;
	case 25:
		for (int k = 0; k < degree; ++k) {
			mu[k] = -cos(M_PI * theta25[k] / (2.0 * degree));
		}
		break;
	case 32:
#pragma omp parallel for
		for (int k = 0; k < degree; ++k) {
			mu[k] = -cos(M_PI * theta32[k] / (2.0 * degree));
		}
		break;
	case 64:
#pragma omp parallel for
		for (int k = 0; k < degree; ++k) {
			mu[k] = -cos(M_PI * theta64[k] / (2.0 * degree));
		}
		break;
	case 128:
#pragma omp parallel for
		for (int k = 0; k < degree; ++k) {
			mu[k] = -cos(M_PI * theta128[k] / (2.0 * degree));
		}
		break;
	case 256:
#pragma omp parallel for
		for (int k = 0; k < degree; ++k) {
			mu[k] = -cos(M_PI * theta256[k] / (2.0 * degree));
		}
		break;
	case 512:
#pragma omp parallel for
		for (int k = 0; k < degree; ++k) {
			mu[k] = -cos(M_PI * theta512[k] / (2.0 * degree));
		}
		break;
	default:

		std::cout << "fatal error !!! unknown Chebyshev polynom degree." << std::endl;
		std::cout << "Chebyshev degree may be equal : 5, 8, 16, 25, 32, 64, 128, 256, 512." << std::endl;
		system("pause");
		exit(1);
		break;
	}
}



// Корректировка нижней границы спектра линейного оператора А.
// См. Жуков, Феодоритова.
// 13.02.2022
void correct_lo(doublereal& lo, bool& bmem, const integer degree,
	const doublereal hi, const doublereal delta, const bool bprintMessage) {
	// Итерационное уточнение нижней границы спектра.
		   // 1.
	doublereal etta = lo / hi;
	doublereal rho1 = (1.0 - sqrt(etta)) / (1.0 + sqrt(etta));
	long double qp = 2.0 * pow(rho1, 1.0 * degree) / (1.0 + pow(rho1, 2.0 * degree));

	if (qp > 1.0e-40) {

		// 2.
		long double y1 = delta / qp;

		if ((y1 > 1.00005)/*&&(y1<1.0e71)*/) {

			bmem = true;

			long double y2 = log(y1 + sqrt(y1 * y1 - 1.0));
			// 3.
			long double x_zv = cosh(y2 / (1.0 * degree));
			// 4.
			long double etta_new = 0.5 * (1.0 + etta) - 0.5 * (1.0 - etta) * x_zv;
			// 5. 
			if (bprintMessage) {
				std::cout << "\n x_zv=" << x_zv << " y2=" << y2 << " y1=" << y1 << " qp=" << qp << " rho1=" << rho1 << " etta=" << etta << " delta=" << delta << std::endl;
				std::cout << "apostoriory lo = " << lo << " hi = " << hi << "  etta_new=" << etta_new << std::endl;
				// getchar();
			}

			// Уточнение нижней границы спектра один раз в degree.
			if (etta_new != etta_new) {
				std::cout << "\n x_zv=" << x_zv << " y2=" << y2 << " y1=" << y1 << " qp=" << qp << " rho1=" << rho1 << " etta=" << etta << " delta=" << delta << std::endl;
				std::cout << "apostoriory lo = " << lo << " hi = " << hi << "  etta_new=" << etta_new << std::endl;
				system("pause");

				lo = 1.0e-10;
			}
			else {
				lo = (1.0e-20 > etta_new * hi ? 1.0e-20 : etta_new * hi);
			}
		}
		else {
			if (bprintMessage) {
				std::cout << "\n etta=" << etta << " rho1=" << rho1 << std::endl;
				std::cout << "\n delta= " << delta << "  qp=" << qp << " y1=" << y1 << std::endl;
			}
			/*if (y1 > 1.0e71) {
				bmem = true;
			}
			else*/ {
				bmem = false;
			}
		}
	}
	else {

		if (bprintMessage) {
			long double y1 = delta / qp;

			std::cout << "\n etta=" << etta << " rho1=" << rho1 << std::endl;
			std::cout << "\n delta= " << delta << "  qp=" << qp << " y1=" << y1 << std::endl;
		}

		lo = lo / 5.0;
	}
}


// применяется для уравнения поправки давления. 
// Работает только на структурированной сетке данная версия.
// 05.02.2022.
// 
// Основываясь на идеях алгоритма Федоренко этот метод можно применять как сглаживатель.
//
//  Работает на видеокарте nvidia. 19.02.2022
//
bool chebyshev3Dnow_cuda(equation3D*& sl, equation3D_bon*& slb, doublereal*& b,
	doublereal*& x, integer maxelm, integer maxbound, integer iVar,
	bool bprintMessage, bool bsolidT) {
	// К величине xcor - будет осуществляться нижняя релаксация, т.е. x_cor - это
	// скорректированная компонента скорости удовлетворяющая уравнению неразрывности.
	//printf("chebyshev3Dnow incomming...\n"); // debug.
	//system("pause");

	// Параметр релаксации лежит в интервале от 0.0 до 2.0.
	// Верхняя релаксация способна существенно ускорить вычислительный процесс. 
	// Патанкар рекомендует брать коэффициент верхней релаксации равный 1.5;
	// rURF=1.5 можно использовать в уравнении на поправку давления.
	// В уравнениях на скорость мы имеем двоякую релаксацию. С одной стороны это нижняя релаксация
	// к скорректированной скорости с коэффициентом alpha. Эта нижняя релаксация учитывается в матрице СЛАУ,
	// в результате чего мы получаем новую матрицу СЛАУ.  К этой новой преоблразованной матрице СЛАУ казалось бы в целях ускорения сходимости
	// можно применить верхнюю релаксацию, пусть она будет опять с коэффициентом rURF=1.5. Вычисления показывают что при введении в уравнения на скорость
	// коэффициента верхней релаксации 1.5 точность решения уравнений на скорость за 4000 итераций падает (или вообще имеем расходимость),
	// поэтому наверное лучше не вводить коээффициент 
	// верхней релаксации равный 1.5 в уравнениях на скорость. 

	 // А.А. Самарский, Е.С. Николаев Методы решения сеточных уравнений. М. Наука, 1978. страница 270-271.

#ifdef _OPENMP 
	omp_set_num_threads(number_cores()); // установка числа потоков
#endif

	if (bsolidT) {
		// диагональное предобуславливание
		// для температуры в твёрдом теле.
		for (int i = 0; i < maxelm; ++i) {
			sl[i].ae /= sl[i].ap;
			sl[i].aw /= sl[i].ap;
			sl[i].an /= sl[i].ap;
			sl[i].as /= sl[i].ap;
			sl[i].at /= sl[i].ap;
			sl[i].ab /= sl[i].ap;

			b[sl[i].iP] /= sl[i].ap;
			sl[i].ap = 1.0;
		}

		for (int i = 0; i < maxbound; ++i) {
			slb[i].ai /= slb[i].aw;
			b[slb[i].iW] /= slb[i].aw;

			slb[i].aw = 1.0;
		}
	}

	/// Chebyshev polynomial degree.
	const integer degree = my_amg_manager.Chebyshev_degree;// 5;  64; 128;

	/// highest eigen value safety upscaling.
	// use boosting factor for a more conservative upper bound estimate
	// See: Adams, Brezina, Hu, Tuminaro,
	//      PARALLEL MULTIGRID SMOOTHING: POLYNOMIAL VERSUS
	//      GAUSS-SEIDEL, J. Comp. Phys. 188 (2003) 593-610.
	//
	float higher;

	/// Lowest-to-highest eigen value ratio.
	float lower;

	// Number of power iterations to apply for the spectral radius
	// estimation. When 0, use Gershgorin disk theorem to estimate
	// spectral radius.
	//integer power_iters;

	// Scale the system matrix
   // bool scale;


	higher = 1.0f;
	lower = 1.0f / 3.0f;// 1.0f / 30;
	//if (inumiterSIMPLE371 > 30) {
	  //  lower = free_debug_parametr1;
   // }
	//lower = 1.0f / 60; // достаточно одной тридцатой.
	//power_iters = 0;
   // scale = false;


	doublereal hi, lo;

	hi = spectral_radius(sl, slb, maxelm, maxbound);

	lo = hi * lower;


	doublereal rURF = 1.0; // параметр верхней релаксации

	hi *= higher;

	// Centre of ellipse containing the eigenvalues of A:
	//doublereal d = 0.5 * (hi + lo);

	// Semi-major axis of ellipse containing the eigenvalues of A:
	//doublereal c = 0.5 * (hi - lo);


	//static const scalar_type one = math::identity<scalar_type>();
	//static const scalar_type zero = math::zero<scalar_type>();







	// Переупорядоченное множество корней многочлена П.Л. Чебышева
	doublereal* mu = new doublereal[degree];
	Lebedev_Samarskii_Nikolaev(degree, mu);

	// пороговое значение невязки
	doublereal eps = 1e-40;


	integer  j = 0, kend = 20;//100; // Для целей пост сглаживания должно хватить 40 итераций.


	const integer size0 = maxelm + maxbound;

	doublereal* x_new = new doublereal[size0];
	doublereal* residual = new doublereal[size0];
#pragma omp parallel for
	for (integer i = 0; i < size0; ++i)
	{
		x_new[i] = 0.0;
		residual[i] = 0.0;
	}

	{
		// Диапазон изменения коэффициентов диффузии в уравнении на поправку давления.
		// опираемся на опыт теплопередачи.
		doublereal c1 = 0.026; // воздух
		doublereal c2 = 2000; // алмаз

		// см. Николаев
		if (bprintMessage) {
			std::cout << " number equations " << maxelm + maxbound << std::endl;
		}
		kend = static_cast<integer>(0.32 * sqrt(c2 / c1) * sqrt(size0) * log(2.0 / 1.0e-4) / degree);
		if (bprintMessage) {
			std::cout << "number out iteration = " << kend << "  " << kend * degree << std::endl;
		}
	}



	// Улучшает попадание в кеш ускорение на 12.5%.
	doublereal* data = new doublereal[maxelm * 7];
	int* col_ind = new int[maxelm * 7];

#pragma omp parallel for 
	for (int i = 0; i < maxelm; ++i) {

		const int is = i * 7;
		data[is] = sl[i].ae;
		data[is + 1] = sl[i].aw;
		data[is + 2] = sl[i].an;
		data[is + 3] = sl[i].as;
		data[is + 4] = sl[i].at;
		data[is + 5] = sl[i].ab;
		data[is + 6] = sl[i].ap;

		col_ind[is] = sl[i].iE;
		col_ind[is + 1] = sl[i].iW;
		col_ind[is + 2] = sl[i].iN;
		col_ind[is + 3] = sl[i].iS;
		col_ind[is + 4] = sl[i].iT;
		col_ind[is + 5] = sl[i].iB;
		col_ind[is + 6] = sl[i].iP;
	}


	doublereal* data_b = new doublereal[maxbound * 2];
	int* col_ind_b = new int[maxbound * 2];

#pragma omp parallel for 
	for (int i = 0; i < maxbound; ++i) {

		const int is = i * 2;
		data_b[is] = slb[i].ai;
		data_b[is + 1] = slb[i].aw;

		if (slb[i].iI > -1) {
			col_ind_b[is] = slb[i].iI;
		}
		else {
			data_b[is] = 0.0;
			col_ind_b[is] = slb[i].iW;
		}
		col_ind_b[is + 1] = slb[i].iW;
	}



#if MY_GPU_MATH_LIBRARY_CU_ON

	if (init_b_first_device) {

		int device;
		device = idevice_Tesla;

		// Choose which GPU to run on, change this on a multi-GPU system.
		cudaStatus = cudaSetDevice(device);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?");
			system("PAUSE");
			exit(1);
		}

		cudaDeviceProp prop;
		cudaGetDeviceProperties(&prop, device);
		printf("%s\n", prop.name);
		// GeForce 840M 384потока 1ГГц каждый март 2014 года. 28нм.

		init_b_first_device = false;
	}

	doublereal* dev_data = nullptr;
	int* dev_col_ind = nullptr;

	doublereal* dev_data_b = nullptr;
	int* dev_col_ind_b = nullptr;


	doublereal* dev_x_new = nullptr;
	doublereal* dev_x = nullptr;
	doublereal* dev_residual = nullptr;
	doublereal* dev_b = nullptr;
	doublereal* sum = (doublereal*)malloc(blocksPerGrid(size0) * sizeof(doublereal));
	doublereal* dev_sum;

	// Allocate GPU buffers for three vectors (two input, one output)    .
	cudaStatus = cudaMalloc((void**)&dev_sum, blocksPerGrid(size0) * sizeof(doublereal));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc dev_sum failed!");
		system("PAUSE");
		exit(1);
	}

	// Allocate GPU buffers for three vectors (two input, one output)    .
	cudaStatus = cudaMalloc((void**)&dev_data, (maxelm * 7) * sizeof(doublereal));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc dev_data failed!");
		system("PAUSE");
		exit(1);
	}

	// Allocate GPU buffers for three vectors (two input, one output)    .
	cudaStatus = cudaMalloc((void**)&dev_col_ind, (maxelm * 7) * sizeof(int));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc dev_col_ind failed!");
		system("PAUSE");
		exit(1);
	}

	// Allocate GPU buffers for three vectors (two input, one output)    .
	cudaStatus = cudaMalloc((void**)&dev_data_b, (maxbound * 2) * sizeof(doublereal));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc dev_data_b failed!");
		system("PAUSE");
		exit(1);
	}

	// Allocate GPU buffers for three vectors (two input, one output)    .
	cudaStatus = cudaMalloc((void**)&dev_col_ind_b, (maxbound * 2) * sizeof(int));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc dev_col_ind_b failed!");
		system("PAUSE");
		exit(1);
	}

	// Allocate GPU buffers for three vectors (two input, one output)    .
	cudaStatus = cudaMalloc((void**)&dev_x_new, (size0) * sizeof(doublereal));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc dev_xnew failed!");
		system("PAUSE");
		exit(1);
	}

	// Allocate GPU buffers for three vectors (two input, one output)    .
	cudaStatus = cudaMalloc((void**)&dev_b, (size0) * sizeof(doublereal));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc dev_b failed!");
		system("PAUSE");
		exit(1);
	}

	// Allocate GPU buffers for three vectors (two input, one output)    .
	cudaStatus = cudaMalloc((void**)&dev_x, (size0) * sizeof(doublereal));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc dev_x failed!");
		system("PAUSE");
		exit(1);
	}

	// Allocate GPU buffers for three vectors (two input, one output)    .
	cudaStatus = cudaMalloc((void**)&dev_residual, (size0) * sizeof(doublereal));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc dev_residual failed!");
		system("PAUSE");
		exit(1);
	}

	// Copy input vectors from host memory to GPU buffers.
	cudaStatus = cudaMemcpy(dev_data, data, maxelm * 7 * sizeof(doublereal), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy dev_data, val HostToDevice failed!");
		system("PAUSE");
		exit(1);
	}

	// Copy input vectors from host memory to GPU buffers.
	cudaStatus = cudaMemcpy(dev_col_ind, col_ind, maxelm * 7 * sizeof(int), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy dev_col_ind, val HostToDevice failed!");
		system("PAUSE");
		exit(1);
	}

	// Copy input vectors from host memory to GPU buffers.
	cudaStatus = cudaMemcpy(dev_data_b, data_b, maxbound * 2 * sizeof(doublereal), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy dev_data_b, val HostToDevice failed!");
		system("PAUSE");
		exit(1);
	}

	// Copy input vectors from host memory to GPU buffers.
	cudaStatus = cudaMemcpy(dev_col_ind_b, col_ind_b, maxbound * 2 * sizeof(int), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy dev_col_ind_b, val HostToDevice failed!");
		system("PAUSE");
		exit(1);
	}

	// Copy input vectors from host memory to GPU buffers.
	cudaStatus = cudaMemcpy(dev_x_new, x_new, size0 * sizeof(doublereal), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy dev_xnew, val HostToDevice failed!");
		system("PAUSE");
		exit(1);
	}


	// Copy input vectors from host memory to GPU buffers.
	cudaStatus = cudaMemcpy(dev_b, b, size0 * sizeof(doublereal), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy dev_b val HostToDevice failed!");
		system("PAUSE");
		exit(1);
	}


	// Copy input vectors from host memory to GPU buffers.
	cudaStatus = cudaMemcpy(dev_x, x, size0 * sizeof(doublereal), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy dev_x, val HostToDevice failed!");
		system("PAUSE");
		exit(1);
	}

	// Copy input vectors from host memory to GPU buffers.
	cudaStatus = cudaMemcpy(dev_residual, residual, size0 * sizeof(doublereal), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy dev_residual, val HostToDevice failed!");
		system("PAUSE");
		exit(1);
	}

	doublereal dmax = 1.0, dmax2 = 0.0, delta = 0.0;
	integer iout = 0;
	while ((dmax > eps) && (j < kend * degree)) {

		if (bsolidT) {
			if (dmax < 1.0e-4) break; // Вычисление сошлось.
		}

		bool bmem = true;

		doublereal r_start;

		const doublereal tau0 = 2.0 / (hi + lo);
		const doublereal ksi = lo / hi;
		const doublereal rho0 = (1.0 - ksi) / (1.0 + ksi);

		for (int k = 0; k < degree; ++k) {



			rURF = tau0 / (1.0 + rho0 * mu[k]);

			dmax = 0.0;
			dmax2 = 0.0;
			delta = 0.0;

			doublereal dmax_loc = 0.0, dmax2_loc = 0.0;

			if ((k == 0) || (k == degree - 1)) {

				Chebyshev_kernel_internal2 << <128, 128 >> > (dev_data, dev_col_ind, dev_b,
					dev_x_new, dev_x, dev_residual, rURF, maxelm);

				// Check for any errors launching the kernel
				cudaStatus = cudaGetLastError();
				if (cudaStatus != cudaSuccess) {
					fprintf(stderr, "Jacoby_kernel_internal launch failed: %s\n", cudaGetErrorString(cudaStatus));
					system("PAUSE");
					exit(1);
				}

				//cudaStatus = cudaMemcpy(residual, dev_residual, size0 * sizeof(doublereal), cudaMemcpyDeviceToHost);
				//if (cudaStatus != cudaSuccess) {
					//fprintf(stderr, "cudaMemcpy residual, dev_residual DeviceToHost failed!");
					//system("PAUSE");
					//exit(1);
				//}

//#pragma omp parallel for  reduction(+: dmax_loc, dmax2_loc)
	///			for (int i = 0; i < maxelm; ++i)
		//		{
			//		const int is = i * 7;
				//	const int is6 = is + 6;
					//const int it6 = col_ind[is6];

					//dmax_loc += residual[it6];
					//dmax2_loc += residual[it6] * residual[it6];

				//}

				//dmax += dmax_loc;
				//dmax2 += dmax2_loc;

				//dmax_loc = 0.0;
				//dmax2_loc = 0.0;

			}
			else {
				Chebyshev_kernel_internal << <128, 128 >> > (dev_data, dev_col_ind, dev_b,
					dev_x_new, dev_x, rURF, maxelm);

				// Check for any errors launching the kernel
				cudaStatus = cudaGetLastError();
				if (cudaStatus != cudaSuccess) {
					fprintf(stderr, "Jacoby_kernel_internal launch failed: %s\n", cudaGetErrorString(cudaStatus));
					system("PAUSE");
					exit(1);
				}
			}


			if ((k == 0) || (k == degree - 1)) {

				Chebyshev_kernel_bound2 << <128, 128 >> > (dev_data_b, dev_col_ind_b, dev_b, dev_x_new, dev_x, dev_residual, rURF, maxbound);

				// Check for any errors launching the kernel
				cudaStatus = cudaGetLastError();
				if (cudaStatus != cudaSuccess) {
					fprintf(stderr, "Jacoby_kernel_internal launch failed: %s\n", cudaGetErrorString(cudaStatus));
					system("PAUSE");
					exit(1);
				}

				//cudaStatus = cudaMemcpy(residual, dev_residual, size0 * sizeof(doublereal), cudaMemcpyDeviceToHost);
				//if (cudaStatus != cudaSuccess) {
					//fprintf(stderr, "cudaMemcpy residual, dev_residual DeviceToHost failed!");
					//system("PAUSE");
					//exit(1);
				//}

//#pragma omp parallel for reduction(+: dmax_loc, dmax2_loc)
	//			for (int i = 0; i < maxbound; ++i) {


		//			const int is = i * 2;
			//		const int it1 = col_ind_b[is];//iI
				//	const int it2 = col_ind_b[is + 1];//iW                  


					// Чебышева
					//dmax_loc += (residual[it2]); // сумма модулей.
					//dmax2_loc += (residual[it2]) * (residual[it2]);

				//}

				//dmax_loc = thrust::reduce(residual, residual + size0);
				//dmax_loc = thrust::reduce(dev_residual, dev_residual + size0);
				//dmax_loc = thrust::reduce(&dev_residual[0], &dev_residual[size0 - 1], (doublereal) 0.0, thrust::plus<doublereal>());

				reduction_sum_kernel << <blocksPerGrid(size0), threadsPerBlock >> > (dev_residual, size0, dev_sum);

				// Check for any errors launching the kernel
				cudaStatus = cudaGetLastError();
				if (cudaStatus != cudaSuccess) {
					fprintf(stderr, "reduction_sum_kernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
					system("PAUSE");
					exit(1);
				}

				cudaStatus = cudaMemcpy(sum, dev_sum, blocksPerGrid(size0) * sizeof(doublereal), cudaMemcpyDeviceToHost);
				if (cudaStatus != cudaSuccess) {
					fprintf(stderr, "cudaMemcpy residual, dev_residual DeviceToHost failed!");
					system("PAUSE");
					exit(1);
				}

				dmax_loc = 0.0;
				const int isize_cycle = static_cast<int>(blocksPerGrid(size0));

#pragma omp parallel for reduction(+: dmax_loc)
				for (int i89 = 0; i89 < isize_cycle; ++i89) {
					dmax_loc += sum[i89];
				}

				


				sqr_vec_kernel << <128, 128 >> > (dev_residual, size0);

				// Check for any errors launching the kernel
				cudaStatus = cudaGetLastError();
				if (cudaStatus != cudaSuccess) {
					fprintf(stderr, "sqr_vec_kernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
					system("PAUSE");
					exit(1);
				}

				//cudaStatus = cudaMemcpy(residual, dev_residual, size0 * sizeof(doublereal), cudaMemcpyDeviceToHost);
				//if (cudaStatus != cudaSuccess) {
					//fprintf(stderr, "cudaMemcpy residual, dev_residual DeviceToHost failed!");
					//system("PAUSE");
					//exit(1);
				//}

				//dmax2_loc = thrust::reduce(residual, residual + size0);
				//dmax2_loc = thrust::reduce(dev_residual, dev_residual + size0);
				//dmax2_loc = thrust::reduce(&dev_residual[0], &dev_residual[size0-1], (doublereal)0.0, thrust::plus<doublereal>());
			
							

				reduction_sum_kernel << <blocksPerGrid(size0), threadsPerBlock >> > (dev_residual, size0, dev_sum);

				// Check for any errors launching the kernel
				cudaStatus = cudaGetLastError();
				if (cudaStatus != cudaSuccess) {
					fprintf(stderr, "reduction_sum_kernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
					system("PAUSE");
					exit(1);
				}

				cudaStatus = cudaMemcpy(sum, dev_sum, blocksPerGrid(size0) * sizeof(doublereal), cudaMemcpyDeviceToHost);
				if (cudaStatus != cudaSuccess) {
					fprintf(stderr, "cudaMemcpy residual, dev_residual DeviceToHost failed!");
					system("PAUSE");
					exit(1);
				}

				dmax2_loc = 0.0;
				

#pragma omp parallel for reduction(+: dmax2_loc)
				for (int i89 = 0; i89 < isize_cycle; ++i89) {
					dmax2_loc += sum[i89];
				}
			}
			else {

				Chebyshev_kernel_bound << <128, 128 >> > (dev_data_b, dev_col_ind_b, dev_b, dev_x_new, dev_x, rURF, maxbound);

				// Check for any errors launching the kernel
				cudaStatus = cudaGetLastError();
				if (cudaStatus != cudaSuccess) {
					fprintf(stderr, "Jacoby_kernel_internal launch failed: %s\n", cudaGetErrorString(cudaStatus));
					system("PAUSE");
					exit(1);
				}
			}

			dmax += dmax_loc;
			dmax2 += dmax2_loc;

			//*std::min(1.0, 17000.0 / size0) 
			if (bsolidT) {
				if (iout == 0) {
					eps = 1.0e-7 * dmax; // Важнейшее условие выхода из итерационного процесса.
					if (j == 0) {
						std::cout << "start residual =" << dmax << std::endl;
					}
				}
			}
			else {
				if (iout == 0) eps = 0.1 * dmax; // Важнейшее условие выхода из итерационного процесса.
			}

			if (bprintMessage) {
				if ((iVar == PAM) || (iVar == TEMP)) {
					if ((k == 0) || (k == degree - 1)) {
						//dmax/=maxelm;
						if (j % degree == 0) {
#if doubleintprecision == 1
							printf("%lld %lld %e ", j + 1, iout + 1, dmax);
#else
							printf("%d %d %e ", j + 1, iout + 1, dmax);
#endif

						}
					}
				}
			}

			copy_Jacoby_kernel << <128, 128 >> > (dev_x_new, dev_x, size0);

			// Check for any errors launching the kernel
			cudaStatus = cudaGetLastError();
			if (cudaStatus != cudaSuccess) {
				fprintf(stderr, "copy_Jacoby_kernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
				system("PAUSE");
				exit(1);
			}



			if ((bmem) && (k == 0)) r_start = sqrt(dmax2 / (size0));
			if (k == degree - 1) delta = sqrt(dmax2 / (size0));

			if (bprintMessage) {
				if ((k == 0) || (k == degree - 1)) {
					if (j % degree == 0) {
						std::cout << std::endl;
					}
				}
			}

			j++;
		}// degree chebyshev


		delta /= r_start;


		// Корректировка нижней границы спектра lo линейного оператора А.
			// Жуков, Феодоритова. 13.02.2022
		correct_lo(lo, bmem, degree, hi, delta, bprintMessage);


		


		//getchar();

		++iout;
	}


	


	cudaStatus = cudaMemcpy(x, dev_x, size0 * sizeof(doublereal), cudaMemcpyDeviceToHost);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy x, dev_x DeviceToHost failed!");
		system("PAUSE");
		exit(1);
	}

	cudaFree(dev_data);
	cudaFree(dev_col_ind);

	cudaFree(dev_data_b);
	cudaFree(dev_col_ind_b);

	cudaFree(dev_b);
	cudaFree(dev_x_new);
	cudaFree(dev_x);
	cudaFree(dev_residual);

	cudaFree(dev_sum);
	free(sum);


#else



	doublereal dmax = 1.0, dmax2 = 0.0, delta = 0.0;
	integer iout = 0;
	while ((dmax > eps) && (j < kend * degree)) {

		if (bsolidT) {
			if (dmax < 1.0e-4) break; // Вычисление сошлось.
		}

		bool bmem = true;

		doublereal r_start;

		const doublereal tau0 = 2.0 / (hi + lo);
		const doublereal ksi = lo / hi;
		const doublereal rho0 = (1.0 - ksi) / (1.0 + ksi);

		for (int k = 0; k < degree; ++k) {



			rURF = tau0 / (1.0 + rho0 * mu[k]);

			dmax = 0.0;
			dmax2 = 0.0;
			delta = 0.0;

			doublereal dmax_loc = 0.0, dmax2_loc = 0.0;

			if ((k == 0) || (k == degree - 1)) {

#pragma omp parallel for reduction(+: dmax_loc, dmax2_loc)
				for (int i = 0; i < maxelm; ++i) {



					//long double residualQ = (sl[i].ae * x[sl[i].iE] +
					  //              sl[i].aw * x[sl[i].iW] +
						//            sl[i].an * x[sl[i].iN] +
						  //          sl[i].as * x[sl[i].iS] +
							//        sl[i].at * x[sl[i].iT] +
							  //      sl[i].ab * x[sl[i].iB] +
								//    b[sl[i].iP]) - sl[i].ap * x[sl[i].iP];

					const int is = i * 7;
					const int is1 = is + 1;
					const int is2 = is + 2;
					const int is3 = is + 3;
					const int is4 = is + 4;
					const int is5 = is + 5;
					const int is6 = is + 6;
					const int it6 = col_ind[is6];

					const long double residualQ = (data[is] * x[col_ind[is]] +
						data[is1] * x[col_ind[is1]] +
						data[is2] * x[col_ind[is2]] +
						data[is3] * x[col_ind[is3]] +
						data[is4] * x[col_ind[is4]] +
						data[is5] * x[col_ind[is5]] +
						b[it6]) - data[is6] * x[it6];

					const doublereal residual = static_cast<doublereal>(residualQ);

					// Чебышева
					dmax_loc += fabs((residual)); // сумма модулей
					dmax2_loc += residual * residual;
					x_new[it6] = x[it6] + rURF * (residual);


				}
			}
			else {

				// Редукция медленная и здесь мы от неё отказываемся в пользу скорости
				// вычисления т.к. внутри она ненужна.

#pragma omp parallel for 
				for (int i = 0; i < maxelm; ++i) {



					//long double residualQ = (sl[i].ae * x[sl[i].iE] +
					  //              sl[i].aw * x[sl[i].iW] +
						//            sl[i].an * x[sl[i].iN] +
						  //          sl[i].as * x[sl[i].iS] +
							//        sl[i].at * x[sl[i].iT] +
							  //      sl[i].ab * x[sl[i].iB] +
								//    b[sl[i].iP]) - sl[i].ap * x[sl[i].iP];

					const int is = i * 7;
					const int is1 = is + 1;
					const int is2 = is + 2;
					const int is3 = is + 3;
					const int is4 = is + 4;
					const int is5 = is + 5;
					const int is6 = is + 6;
					const int it6 = col_ind[is6];

					const long double residualQ = (data[is] * x[col_ind[is]] +
						data[is1] * x[col_ind[is1]] +
						data[is2] * x[col_ind[is2]] +
						data[is3] * x[col_ind[is3]] +
						data[is4] * x[col_ind[is4]] +
						data[is5] * x[col_ind[is5]] +
						b[it6]) - data[is6] * x[it6];

					const doublereal residual = static_cast<doublereal>(residualQ);


					x_new[it6] = x[it6] + rURF * (residual);

				}

			}



			dmax += dmax_loc;
			dmax2 += dmax2_loc;

			dmax_loc = 0.0;
			dmax2_loc = 0.0;


			if ((k == 0) || (k == degree - 1)) {

#pragma omp parallel for reduction(+: dmax_loc, dmax2_loc)
				for (int i = 0; i < maxbound; ++i) {


					const int is = i * 2;
					const int it1 = col_ind_b[is];//iI
					const int it2 = col_ind_b[is + 1];//iW                   


					const doublereal residual = (data_b[is] * x[it1] + b[it2]) - data_b[is + 1] * x[it2];
					// Чебышева
					dmax_loc += fabs(residual); // сумма модулей.
					dmax2_loc += residual * residual;

					if (it1 == it2) x_new[it2] = (b[it2] / data_b[is + 1]);
					else x_new[it2] = x[it2] + rURF * (residual);


				}

			}
			else {
				// Редукция медленная и здесь мы от неё отказываемся в пользу скорости
				// вычисления т.к. внутри она ненужна.

#pragma omp parallel for
				for (int i = 0; i < maxbound; ++i) {


					const int is = i * 2;
					const int it1 = col_ind_b[is];//iI
					const int it2 = col_ind_b[is + 1];//iW                   


					const doublereal residual = (data_b[is] * x[it1] + b[it2]) - data_b[is + 1] * x[it2];

					if (it1 == it2) x_new[it2] = (b[it2] / data_b[is + 1]);
					else x_new[it2] = x[it2] + rURF * (residual);


				}
			}



			dmax += dmax_loc;
			dmax2 += dmax2_loc;

			//*std::min(1.0, 17000.0 / size0) 
			if (bsolidT) {
				if (iout == 0) {
					eps = 1.0e-7 * dmax; // Важнейшее условие выхода из итерационного процесса.
					if (j == 0) {
						std::cout << "start residual =" << dmax << std::endl;
					}
				}
			}
			else {
				if (iout == 0) eps = 0.1 * dmax; // Важнейшее условие выхода из итерационного процесса.
			}

			if (bprintMessage) {
				if ((iVar == PAM) || (iVar == TEMP)) {
					if ((k == 0) || (k == degree - 1)) {
						//dmax/=maxelm;
						if (j % degree == 0) {
#if doubleintprecision == 1
							printf("%lld %lld %e ", j + 1, iout + 1, dmax);
#else
							printf("%d %d %e ", j + 1, iout + 1, dmax);
#endif

						}
					}
				}
			}

#pragma omp parallel for
			for (integer i = 0; i < size0; ++i)
			{
				x[i] = x_new[i];
			}

			if ((bmem) && (k == 0)) r_start = sqrt(dmax2 / (size0));
			if (k == degree - 1) delta = sqrt(dmax2 / (size0));

			if (bprintMessage) {
				if ((k == 0) || (k == degree - 1)) {
					if (j % degree == 0) {
						std::cout << std::endl;
					}
				}
			}

			j++;
		}// degree chebyshev


		delta /= r_start;


		// Корректировка нижней границы спектра lo линейного оператора А.
			// Жуков, Феодоритова. 13.02.2022
		correct_lo(lo, bmem, degree, hi, delta, bprintMessage);		

		//getchar();

		++iout;
	}


#endif

	delete[] data;
	delete[] col_ind;

	delete[] data_b;
	delete[] col_ind_b;

	std::cout << " " << lo << " " << j << " " << dmax << " ";

	if (j > kend * degree - 100) {
		std::cout << "\n chebyshev divergence " << dmax << " ";
		system("pause");
	}

	delete[] x_new;
	delete[] residual;

	if (bprintMessage) {
		if ((iVar == PAM) || (iVar == TEMP)) {
			printf("calc complete...\n");
			system("pause");
		}
	}
	//printf("4000 %e \n", dmax);
	//system("pause");

	if ((dmax < eps)&&(dmax<99.0)) {
		return false; // Сошлось
	}
	else {
		return true; // Разошлось
	}

} // chebyshev3Dnow_cuda

/*
void precontidioner_cuda(const integer degree,
	const doublereal tau0,
	const doublereal rho0,
	doublereal*& mu,
	integer maxelm, integer maxbound,
	doublereal& dmax, doublereal& dmax2, doublereal& delta,
	doublereal& r_start,
	bool bprintMessage, bool bsolidT, integer& j,
	integer iout, integer iVar, integer size0,
	bool bmem, doublereal eps,
	doublereal* dev_data, int* dev_col_ind,
	doublereal* dev_data_b, int* dev_col_ind_b,
	doublereal*& x, doublereal*& b, doublereal*& x_new,
	doublereal*& residual,
	doublereal* dev_x, doublereal* dev_x_new,
	doublereal* dev_residual, doublereal* dev_b)
{
	
	// Copy input vectors from host memory to GPU buffers.
	cudaStatus = cudaMemcpy(dev_b, b, size0 * sizeof(doublereal), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy dev_b val HostToDevice failed!");
		system("PAUSE");
		exit(1);
	}

	// Copy input vectors from host memory to GPU buffers.
	cudaStatus = cudaMemcpy(dev_x, x, size0 * sizeof(doublereal), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy dev_b val HostToDevice failed!");
		system("PAUSE");
		exit(1);
	}

	// Copy input vectors from host memory to GPU buffers.
	cudaStatus = cudaMemcpy(dev_x_new, x_new, size0 * sizeof(doublereal), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy dev_b val HostToDevice failed!");
		system("PAUSE");
		exit(1);
	}

	// std::cout << "degree="<< degree << std::endl;

	for (int k = 0; k < degree; ++k) {



		const doublereal rURF = tau0 / (1.0 + rho0 * mu[k]);

		dmax = 0.0;
		dmax2 = 0.0;
		delta = 0.0;

		doublereal dmax_loc = 0.0, dmax2_loc = 0.0;

		if ((k == 0) || (k == degree - 1)) {

			Chebyshev_kernel_internal2 << <128, 128 >> > (dev_data, dev_col_ind, dev_b,
				dev_x_new, dev_x, dev_residual, rURF, maxelm);

			// Check for any errors launching the kernel
			cudaStatus = cudaGetLastError();
			if (cudaStatus != cudaSuccess) {
				fprintf(stderr, "Jacoby_kernel_internal launch failed: %s\n", cudaGetErrorString(cudaStatus));
				system("PAUSE");
				exit(1);
			}

			//cudaStatus = cudaMemcpy(residual, dev_residual, size0 * sizeof(doublereal), cudaMemcpyDeviceToHost);
			//if (cudaStatus != cudaSuccess) {
				//fprintf(stderr, "cudaMemcpy residual, dev_residual DeviceToHost failed!");
				//system("PAUSE");
				//exit(1);
			//}

//#pragma omp parallel for  reduction(+: dmax_loc, dmax2_loc)
	///			for (int i = 0; i < maxelm; ++i)
		//		{
			//		const int is = i * 7;
				//	const int is6 = is + 6;
					//const int it6 = col_ind[is6];

					//dmax_loc += residual[it6];
					//dmax2_loc += residual[it6] * residual[it6];

				//}

				//dmax += dmax_loc;
				//dmax2 += dmax2_loc;

				//dmax_loc = 0.0;
				//dmax2_loc = 0.0;

		}
		else {
			Chebyshev_kernel_internal << <128, 128 >> > (dev_data, dev_col_ind, dev_b,
				dev_x_new, dev_x, rURF, maxelm);

			// Check for any errors launching the kernel
			cudaStatus = cudaGetLastError();
			if (cudaStatus != cudaSuccess) {
				fprintf(stderr, "Jacoby_kernel_internal launch failed: %s\n", cudaGetErrorString(cudaStatus));
				system("PAUSE");
				exit(1);
			}
		}


		if ((k == 0) || (k == degree - 1)) {

			Chebyshev_kernel_bound2 << <128, 128 >> > (dev_data_b, dev_col_ind_b, dev_b, dev_x_new, dev_x, dev_residual, rURF, maxbound);

			// Check for any errors launching the kernel
			cudaStatus = cudaGetLastError();
			if (cudaStatus != cudaSuccess) {
				fprintf(stderr, "Jacoby_kernel_internal launch failed: %s\n", cudaGetErrorString(cudaStatus));
				system("PAUSE");
				exit(1);
			}

			cudaStatus = cudaMemcpy(residual, dev_residual, size0 * sizeof(doublereal), cudaMemcpyDeviceToHost);
			if (cudaStatus != cudaSuccess) {
				fprintf(stderr, "cudaMemcpy residual, dev_residual DeviceToHost failed!");
				system("PAUSE");
				exit(1);
			}

			//#pragma omp parallel for reduction(+: dmax_loc, dmax2_loc)
				//			for (int i = 0; i < maxbound; ++i) {


					//			const int is = i * 2;
						//		const int it1 = col_ind_b[is];//iI
							//	const int it2 = col_ind_b[is + 1];//iW                  


								// Чебышева
								//dmax_loc += (residual[it2]); // сумма модулей.
								//dmax2_loc += (residual[it2]) * (residual[it2]);

							//}

			dmax_loc = thrust::reduce(residual, residual + size0);

			sqr_vec_kernel << <128, 128 >> > (dev_residual, size0);

			// Check for any errors launching the kernel
			cudaStatus = cudaGetLastError();
			if (cudaStatus != cudaSuccess) {
				fprintf(stderr, "sqr_vec_kernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
				system("PAUSE");
				exit(1);
			}

			cudaStatus = cudaMemcpy(residual, dev_residual, size0 * sizeof(doublereal), cudaMemcpyDeviceToHost);
			if (cudaStatus != cudaSuccess) {
				fprintf(stderr, "cudaMemcpy residual, dev_residual DeviceToHost failed!");
				system("PAUSE");
				exit(1);
			}

			dmax2_loc = thrust::reduce(residual, residual + size0);
		}
		else {

			Chebyshev_kernel_bound << <128, 128 >> > (dev_data_b, dev_col_ind_b, dev_b, dev_x_new, dev_x, rURF, maxbound);

			// Check for any errors launching the kernel
			cudaStatus = cudaGetLastError();
			if (cudaStatus != cudaSuccess) {
				fprintf(stderr, "Jacoby_kernel_internal launch failed: %s\n", cudaGetErrorString(cudaStatus));
				system("PAUSE");
				exit(1);
			}
		}

		dmax += dmax_loc;
		dmax2 += dmax2_loc;



		//std::min(1.0, 17000.0 / size0) 
		if (bsolidT) {
			if (iout == 0) eps = 1.0e-7 * dmax; // Важнейшее условие выхода из итерационного процесса.
		}
		else {
			if (iout == 0) eps = 0.1 * dmax; // Важнейшее условие выхода из итерационного процесса.
		}

		if (bprintMessage)
		{
			if ((iVar == PAM) || (iVar == TEMP))
			{
				if ((k == 0) || (k == degree - 1))
				{
					//dmax/=maxelm;
					if (j % degree == 0)
					{
#if doubleintprecision == 1
						printf("%lld %lld %e ", j + 1, iout + 1, dmax);
#else
						printf("%d %d %e ", j + 1, iout + 1, dmax);
#endif

					}
				}
			}
		}

		copy_Jacoby_kernel << <128, 128 >> > (dev_x_new, dev_x, size0);

		// Check for any errors launching the kernel
		cudaStatus = cudaGetLastError();
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "copy_Jacoby_kernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
			system("PAUSE");
			exit(1);
		}

		if ((bmem) && (k == 0)) r_start = sqrt(dmax2 / (size0));
		if (k == degree - 1) delta = sqrt(dmax2 / (size0));

		if (bprintMessage) {
			if ((k == 0) || (k == degree - 1)) {
				if (j % degree == 0) {
					std::cout << std::endl;
				}
			}
		}

		j++;
	}// degree chebyshev


	cudaStatus = cudaMemcpy(x, dev_x, size0 * sizeof(doublereal), cudaMemcpyDeviceToHost);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy x, dev_x DeviceToHost failed!");
		system("PAUSE");
		exit(1);
	}

	//std::wcout << "finish degree:"<< size0 << r_start << "  " << delta << "\n";
	//getchar();
}
*/




// Норма вектора
// как корень квадратный из суммы квадратов.
double NormaV_for_gmres(double const* const dV, const integer isize)
{

	double dnorma, dsum;

	// инициализация переменных
	dsum = 0.0;

#if MY_GPU_MATH_LIBRARY_CU_ON

	if (init_b_first_device) {

		int device;
		device = idevice_Tesla;

		// Choose which GPU to run on, change this on a multi-GPU system.
		cudaStatus = cudaSetDevice(device);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?");
			system("PAUSE");
			exit(1);
		}

		cudaDeviceProp prop;
		cudaGetDeviceProperties(&prop, device);
		printf("%s\n", prop.name);
		// GeForce 840M 384потока 1ГГц каждый март 2014 года. 28нм.

		init_b_first_device = false;
	}

	dsum = thrust::inner_product(thrust::host, dV, dV + isize, dV, 0.0f);

	dnorma = sqrt(dsum); // норма вектора
#else


	const integer ISIZE = (isize - 1);

#pragma omp parallel for reduction(+ : dsum) 
	for (integer i = 0; i <= ISIZE; ++i)
	{
		dsum += dV[i] * dV[i];
	}
	dnorma = sqrt(dsum); // норма вектора

#endif

	return dnorma;
}// NormaV_for_gmres


double Scal(double const* const v1, double const* const v2, const integer n) {
	double sum_squares = 0.0;


#if MY_GPU_MATH_LIBRARY_CU_ON

	if (init_b_first_device) {

		int device;
		device = idevice_Tesla;

		// Choose which GPU to run on, change this on a multi-GPU system.
		cudaStatus = cudaSetDevice(device);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?");
			system("PAUSE");
			exit(1);
		}

		cudaDeviceProp prop;
		cudaGetDeviceProperties(&prop, device);
		printf("%s\n", prop.name);
		// GeForce 840M 384потока 1ГГц каждый март 2014 года. 28нм.

		init_b_first_device = false;
	}

	sum_squares = thrust::inner_product(thrust::host, v1, v1 + n, v2, 0.0f);

#else


	// Возможно не имеет смысла так часто постоянно задавать количество потоков.
	//#ifdef _OPENMP
		//omp_set_num_threads(inumcore);
	//#endif


#pragma omp parallel for  reduction (+: sum_squares) 
	for (integer i = 0; i < n; ++i)
	{
		sum_squares += v1[i] * v2[i];
	}

#endif

	return sum_squares;
} // Scal


void residual_cg(
	integer maxelm, integer maxbound,
	doublereal*& data, int*& col_ind,
	doublereal*& data_b, int*& col_ind_b,
	doublereal*& x, doublereal*& b, doublereal*& r)
{


#pragma omp parallel for 
	for (int i = 0; i < maxelm; ++i) {


		//long double residualQ = (sl[i].ae * x[sl[i].iE] +
		  //              sl[i].aw * x[sl[i].iW] +
			//            sl[i].an * x[sl[i].iN] +
			  //          sl[i].as * x[sl[i].iS] +
				//        sl[i].at * x[sl[i].iT] +
				  //      sl[i].ab * x[sl[i].iB] +
					//    b[sl[i].iP]) - sl[i].ap * x[sl[i].iP];

		const int is = i * 7;
		const int is1 = is + 1;
		const int is2 = is + 2;
		const int is3 = is + 3;
		const int is4 = is + 4;
		const int is5 = is + 5;
		const int is6 = is + 6;
		const int it6 = col_ind[is6];

		const long double residualQ = (data[is] * x[col_ind[is]] +
			data[is1] * x[col_ind[is1]] +
			data[is2] * x[col_ind[is2]] +
			data[is3] * x[col_ind[is3]] +
			data[is4] * x[col_ind[is4]] +
			data[is5] * x[col_ind[is5]] +
			b[it6]) - data[is6] * x[it6];

		const doublereal residual = static_cast<doublereal>(residualQ);

		// Чебышева

		r[it6] = (residual);

	}


#pragma omp parallel for 
	for (int i = 0; i < maxbound; ++i) {

		const int is = i * 2;
		const int it1 = col_ind_b[is];//iI
		const int it2 = col_ind_b[is + 1];//iW                   


		const doublereal residual = (data_b[is] * x[it1] + b[it2]) - data_b[is + 1] * x[it2];
		// Чебышева


		if (it1 == it2) r[it2] = 0.0;// (b[it2]) - data_b[is + 1] * x[it2];
		else r[it2] = (residual);

	}

} // residual_cg


/*
void residual_cg_cuda(
	integer maxelm, integer maxbound,
	doublereal* data, int* col_ind,
	doublereal* data_b, int* col_ind_b,
	doublereal*& x, doublereal*& b, doublereal*& r,
	doublereal* dev_x, doublereal* dev_b, doublereal* dev_residual)
{

	integer size0 = maxelm + maxbound;

	// Copy input vectors from host memory to GPU buffers.
	cudaStatus = cudaMemcpy(dev_x, x, size0 * sizeof(doublereal), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy dev_x, x HostToDevice failed!");
		system("PAUSE");
		exit(1);
	}

	// Copy input vectors from host memory to GPU buffers.
	cudaStatus = cudaMemcpy(dev_b, b, size0 * sizeof(doublereal), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy dev_b, b HostToDevice failed!");
		system("PAUSE");
		exit(1);
	}

	Chebyshev_kernel_internal_residual2 << <128, 128 >> > (data, col_ind, dev_b, dev_x, dev_residual, maxelm);

	// Check for any errors launching the kernel
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "Chebyshev_kernel_internal_residual2 launch failed: %s\n", cudaGetErrorString(cudaStatus));
		system("PAUSE");
		exit(1);
	}

	Chebyshev_kernel_bound_residual2 << <128, 128 >> > (data_b, col_ind_b, dev_b, dev_x, dev_residual, maxbound);

	// Check for any errors launching the kernel
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "Chebyshev_kernel_bound_residual2 launch failed: %s\n", cudaGetErrorString(cudaStatus));
		system("PAUSE");
		exit(1);
	}

	cudaStatus = cudaMemcpy(r, dev_residual, size0 * sizeof(doublereal), cudaMemcpyDeviceToHost);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy residual, dev_residual DeviceToHost failed!");
		system("PAUSE");
		exit(1);
	}

} // residual_cg_cuda
*/

void precontidioner(const integer degree,
	const doublereal tau0,
	const doublereal rho0,
	doublereal*& mu,
	integer maxelm, integer maxbound,
	doublereal& dmax, doublereal& dmax2, doublereal& delta,
	doublereal& r_start,
	bool bprintMessage, bool bsolidT, integer& j,
	integer iout, integer iVar, integer size0,
	bool bmem, doublereal eps,
	doublereal*& data, int*& col_ind,
	doublereal*& data_b, int*& col_ind_b,
	doublereal*& x, doublereal*& b, doublereal*& x_new)
{

	// std::cout << "degree="<< degree << std::endl;

	for (int k = 0; k < degree; ++k) {



		const doublereal rURF = tau0 / (1.0 + rho0 * mu[k]);

		dmax = 0.0;
		dmax2 = 0.0;
		delta = 0.0;

		doublereal dmax_loc = 0.0, dmax2_loc = 0.0;

		if ((k == 0) || (k == degree - 1)) {

#pragma omp parallel for reduction(+: dmax_loc, dmax2_loc)  //schedule(static, 4)//schedule(dynamic)
			for (int i = 0; i < maxelm; ++i) {



				//long double residualQ = (sl[i].ae * x[sl[i].iE] +
				  //              sl[i].aw * x[sl[i].iW] +
					//            sl[i].an * x[sl[i].iN] +
					  //          sl[i].as * x[sl[i].iS] +
						//        sl[i].at * x[sl[i].iT] +
						  //      sl[i].ab * x[sl[i].iB] +
							//    b[sl[i].iP]) - sl[i].ap * x[sl[i].iP];

				const int is = i * 7;
				const int is1 = is + 1;
				const int is2 = is + 2;
				const int is3 = is + 3;
				const int is4 = is + 4;
				const int is5 = is + 5;
				const int is6 = is + 6;
				const int it6 = col_ind[is6];

				const long double residualQ = (data[is] * x[col_ind[is]] +
					data[is1] * x[col_ind[is1]] +
					data[is2] * x[col_ind[is2]] +
					data[is3] * x[col_ind[is3]] +
					data[is4] * x[col_ind[is4]] +
					data[is5] * x[col_ind[is5]] +
					b[it6]) - data[is6] * x[it6];

				const doublereal residual = static_cast<doublereal>(residualQ);

				// Чебышева
				dmax_loc += fabs((residual)); // сумма модулей
				dmax2_loc += residual * residual;
				x_new[it6] = x[it6] + rURF * (residual);


			}
		}
		else {

			// Редукция медленная и здесь мы от неё отказываемся в пользу скорости
			// вычисления т.к. внутри она ненужна.

#pragma omp parallel for //schedule(static, 4) //schedule(dynamic)
			for (int i = 0; i < maxelm; ++i) {



				//long double residualQ = (sl[i].ae * x[sl[i].iE] +
				  //              sl[i].aw * x[sl[i].iW] +
					//            sl[i].an * x[sl[i].iN] +
					  //          sl[i].as * x[sl[i].iS] +
						//        sl[i].at * x[sl[i].iT] +
						  //      sl[i].ab * x[sl[i].iB] +
							//    b[sl[i].iP]) - sl[i].ap * x[sl[i].iP];

				const int is = i * 7;
				const int is1 = is + 1;
				const int is2 = is + 2;
				const int is3 = is + 3;
				const int is4 = is + 4;
				const int is5 = is + 5;
				const int is6 = is + 6;
				const int it6 = col_ind[is6];

				const long double residualQ = (data[is] * x[col_ind[is]] +
					data[is1] * x[col_ind[is1]] +
					data[is2] * x[col_ind[is2]] +
					data[is3] * x[col_ind[is3]] +
					data[is4] * x[col_ind[is4]] +
					data[is5] * x[col_ind[is5]] +
					b[it6]) - data[is6] * x[it6];

				const doublereal residual = static_cast<doublereal>(residualQ);


				x_new[it6] = x[it6] + rURF * (residual);

			}

		}



		dmax += dmax_loc;
		dmax2 += dmax2_loc;

		dmax_loc = 0.0;
		dmax2_loc = 0.0;


		if ((k == 0) || (k == degree - 1)) {

#pragma omp parallel for reduction(+: dmax_loc, dmax2_loc)// schedule(static, 4)//schedule(dynamic)
			for (int i = 0; i < maxbound; ++i) {


				const int is = i * 2;
				const int it1 = col_ind_b[is];//iI
				const int it2 = col_ind_b[is + 1];//iW                   


				const doublereal residual = (data_b[is] * x[it1] + b[it2]) - data_b[is + 1] * x[it2];
				// Чебышева
				dmax_loc += fabs(residual); // сумма модулей.
				dmax2_loc += residual * residual;

				if (it1 == it2) x_new[it2] = (b[it2] / data_b[is + 1]);
				else x_new[it2] = x[it2] + rURF * (residual);


			}

		}
		else {
			// Редукция медленная и здесь мы от неё отказываемся в пользу скорости
			// вычисления т.к. внутри она ненужна.

#pragma omp parallel for //schedule(static, 4) //schedule(dynamic)
			for (int i = 0; i < maxbound; ++i) {


				const int is = i * 2;
				const int it1 = col_ind_b[is];//iI
				const int it2 = col_ind_b[is + 1];//iW                   


				const doublereal residual = (data_b[is] * x[it1] + b[it2]) - data_b[is + 1] * x[it2];

				if (it1 == it2) x_new[it2] = (b[it2] / data_b[is + 1]);
				else x_new[it2] = x[it2] + rURF * (residual);


			}
		}



		dmax += dmax_loc;
		dmax2 += dmax2_loc;

		//*std::min(1.0, 17000.0 / size0) 
		if (bsolidT) {
			if (iout == 0) eps = 1.0e-7 * dmax; // Важнейшее условие выхода из итерационного процесса.
		}
		else {
			if (iout == 0) eps = 0.1 * dmax; // Важнейшее условие выхода из итерационного процесса.
		}

		if (bprintMessage)
		{
			if ((iVar == PAM) || (iVar == TEMP))
			{
				if ((k == 0) || (k == degree - 1))
				{
					//dmax/=maxelm;
					if (j % degree == 0)
					{
#if doubleintprecision == 1
						printf("%lld %lld %e ", j + 1, iout + 1, dmax);
#else
						printf("%d %d %e ", j + 1, iout + 1, dmax);
#endif

					}
				}
			}
		}

#pragma omp parallel for
		for (integer i = 0; i < size0; ++i)
		{

			x[i] = x_new[i];
		}

		if ((bmem) && (k == 0)) r_start = sqrt(dmax2 / (size0));
		if (k == degree - 1) delta = sqrt(dmax2 / (size0));

		if (bprintMessage) {
			if ((k == 0) || (k == degree - 1)) {
				if (j % degree == 0) {
					std::cout << std::endl;
				}
			}
		}

		j++;
	}// degree chebyshev

	//std::wcout << "finish degree:"<< size0 << r_start << "  " << delta << "\n";
	//getchar();
}


doublereal dres_init_old_Simple_iteration = 1.0e30;

// . 
// Сделан на структурированной сетке. 
// 05.02.2022.
// 
// Метод BiCGStab предобусловленный методом Чебышева. 12.02.2022.
//
// Работает на видеокарте NVIDIA. 19.02.2022
//
bool bicgstab_chebyshev_3Dnow_cuda(equation3D*& sl, equation3D_bon*& slb, doublereal*& b,
	doublereal*& x, integer maxelm, integer maxbound, integer iVar,
	bool bprintMessage, bool bsolidT, int lb)
{
	// К величине xcor - будет осуществляться нижняя релаксация, т.е. x_cor - это
	// скорректированная компонента скорости удовлетворяющая уравнению неразрывности.
	//printf("chebyshev3Dnow incomming...\n"); // debug.
	//system("pause");

	// Параметр релаксации лежит в интервале от 0.0 до 2.0.
	// Верхняя релаксация способна существенно ускорить вычислительный процесс. 
	// Патанкар рекомендует брать коэффициент верхней релаксации равный 1.5;
	// rURF=1.5 можно использовать в уравнении на поправку давления.
	// В уравнениях на скорость мы имеем двоякую релаксацию. С одной стороны это нижняя релаксация
	// к скорректированной скорости с коэффициентом alpha. Эта нижняя релаксация учитывается в матрице СЛАУ,
	// в результате чего мы получаем новую матрицу СЛАУ.  К этой новой преоблразованной матрице СЛАУ казалось бы в целях ускорения сходимости
	// можно применить верхнюю релаксацию, пусть она будет опять с коэффициентом rURF=1.5. Вычисления показывают что при введении в уравнения на скорость
	// коэффициента верхней релаксации 1.5 точность решения уравнений на скорость за 4000 итераций падает (или вообще имеем расходимость),
	// поэтому наверное лучше не вводить коээффициент 
	// верхней релаксации равный 1.5 в уравнениях на скорость. 

	 // А.А. Самарский, Е.С. Николаев Методы решения сеточных уравнений. М. Наука, 1978. страница 270-271.

#ifdef _OPENMP 
	omp_set_num_threads(number_cores()); // установка числа потоков
#endif


	// Чистый cg лучше сходится без деления на диагональ.
	if (1 && bsolidT) {
		// диагональное предобуславливание
		// для температуры в твёрдом теле.
		for (int i = 0; i < maxelm; ++i) {
			sl[i].ae /= sl[i].ap;
			sl[i].aw /= sl[i].ap;
			sl[i].an /= sl[i].ap;
			sl[i].as /= sl[i].ap;
			sl[i].at /= sl[i].ap;
			sl[i].ab /= sl[i].ap;

			b[sl[i].iP] /= sl[i].ap;
			sl[i].ap = 1.0;
		}

		for (int i = 0; i < maxbound; ++i) {
			slb[i].ai /= slb[i].aw;
			b[slb[i].iW] /= slb[i].aw;

			slb[i].aw = 1.0;
		}
	}

	std::vector<double> D(maxelm + maxbound, 1.0);

	if (0 && bsolidT) {

		// Scale the matrix so that it has the unit diagonal.
   // First, find the diagonal values:
#pragma omp parallel for 
		for (int i = 0; i < maxelm; ++i) {
			D[sl[i].iP] = 1 / sqrt(fabs(sl[i].ap));
		}

#pragma omp parallel for 
		for (int i = 0; i < maxbound; ++i) {
			D[slb[i].iW] = 1 / sqrt(fabs(slb[i].aw));
		}

		// Then, apply the scaling in-place:
#pragma omp parallel for 
		for (int i = 0; i < maxelm; ++i) {
			sl[i].ae *= D[sl[i].iP] * D[sl[i].iE];
			sl[i].aw *= D[sl[i].iP] * D[sl[i].iW];
			sl[i].an *= D[sl[i].iP] * D[sl[i].iN];
			sl[i].as *= D[sl[i].iP] * D[sl[i].iS];
			sl[i].at *= D[sl[i].iP] * D[sl[i].iT];
			sl[i].ab *= D[sl[i].iP] * D[sl[i].iB];
			sl[i].ap *= D[sl[i].iP] * D[sl[i].iP];

			b[sl[i].iP] *= D[sl[i].iP];
		}

#pragma omp parallel for 
		for (int i = 0; i < maxbound; ++i) {
			if (slb[i].iI > -1) {
				slb[i].ai *= D[slb[i].iW] * D[slb[i].iI];
			}
			slb[i].aw *= D[slb[i].iW] * D[slb[i].iW];

			b[slb[i].iW] *= D[slb[i].iW];
		}



	}

	/// Chebyshev polynomial degree.
	const integer degree = my_amg_manager.Chebyshev_degree;// 5;  64; 128;

	/// highest eigen value safety upscaling.
	// use boosting factor for a more conservative upper bound estimate
	// See: Adams, Brezina, Hu, Tuminaro,
	//      PARALLEL MULTIGRID SMOOTHING: POLYNOMIAL VERSUS
	//      GAUSS-SEIDEL, J. Comp. Phys. 188 (2003) 593-610.
	//
	float higher;

	/// Lowest-to-highest eigen value ratio.
	float lower;

	// Number of power iterations to apply for the spectral radius
	// estimation. When 0, use Gershgorin disk theorem to estimate
	// spectral radius.
	//integer power_iters;

	// Scale the system matrix
   // bool scale;


	higher = 1.0f;
	lower = 1.0f / 3.0f;// 1.0f / 30;
	//if (inumiterSIMPLE371 > 30) {
	  //  lower = free_debug_parametr1;
   // }
	//lower = 1.0f / 60; // достаточно одной тридцатой.
	//power_iters = 0;
   // scale = false;


	doublereal hi, lo; // верхняя hi и нижняя lo границы спектра.

	// Вычисление верхней границы спектра по теореме Гершгорина о кругах.
	hi = spectral_radius(sl, slb, maxelm, maxbound);

	lo = hi * lower;


	

	hi *= higher;

	// Centre of ellipse containing the eigenvalues of A:
	//doublereal d = 0.5 * (hi + lo);

	// Semi-major axis of ellipse containing the eigenvalues of A:
	//doublereal c = 0.5 * (hi - lo);


	//static const scalar_type one = math::identity<scalar_type>();
	//static const scalar_type zero = math::zero<scalar_type>();



	// Переупорядоченное множество корней многочлена П.Л. Чебышева
	doublereal* mu = new doublereal[degree];
	Lebedev_Samarskii_Nikolaev(degree, mu);

	// пороговое значение невязки
	doublereal eps = 1e-40;


	integer  j = 0, kend = 20;//100; // Для целей пост сглаживания должно хватить 40 итераций.

	if (iVar == TEMP) {
		if (maxbound < maxelm) {
			// Т.к. данный метод предназначен только для структурированной сетки то мы 
			// переливаем данные из slb на свободные места из-за АЛИС в sl не тратя ни грамма лишненй памяти.
#pragma omp parallel for 
			for (int i = 0; i < maxbound; ++i) {

				sl[i].ab2 = slb[i].ai;
				sl[i].ab3 = slb[i].aw;
				sl[i].iB2 = slb[i].iI;
				sl[i].iB3 = slb[i].iW;
			}

			// Сразу освобождаем оперативную память.
			delete[] slb;		

		}
	}

	doublereal* data_b = new doublereal[maxbound * 2];
	int* col_ind_b = new int[maxbound * 2];

	if ((iVar == PAM) || (maxbound >= maxelm)) {
#pragma omp parallel for
		for (int i = 0; i < maxbound; ++i) {

			const int is = i * 2;
			data_b[is] = slb[i].ai;
			data_b[is + 1] = slb[i].aw;

			if (slb[i].iI > -1) {
				col_ind_b[is] = slb[i].iI;
			}
			else {
				data_b[is] = 0.0;
				col_ind_b[is] = slb[i].iW;
			}
			col_ind_b[is + 1] = slb[i].iW;
		}

		// Сразу освобождаем оперативную память.
		if (iVar == TEMP) {
			delete[] slb;
		}
	}
	else {

#pragma omp parallel for
		for (int i = 0; i < maxbound; ++i) {

			const int is = i * 2;
			data_b[is] = sl[i].ab2;
			data_b[is + 1] = sl[i].ab3;

			if (sl[i].iB2 > -1) {
				col_ind_b[is] = sl[i].iB2;
			}
			else {
				data_b[is] = 0.0;
				col_ind_b[is] = sl[i].iB3;
			}
			col_ind_b[is + 1] = sl[i].iB3;
		}
	}


	// Улучшает попадание в кеш ускорение на 12.5%.
	doublereal* data = new doublereal[maxelm * 7];
	int* col_ind = new int[maxelm * 7];

#pragma omp parallel for 
	for (int i = 0; i < maxelm; ++i) {

		const int is = i * 7;
		data[is] = sl[i].ae;
		data[is + 1] = sl[i].aw;
		data[is + 2] = sl[i].an;
		data[is + 3] = sl[i].as;
		data[is + 4] = sl[i].at;
		data[is + 5] = sl[i].ab;
		data[is + 6] = sl[i].ap;

		col_ind[is] = sl[i].iE;
		col_ind[is + 1] = sl[i].iW;
		col_ind[is + 2] = sl[i].iN;
		col_ind[is + 3] = sl[i].iS;
		col_ind[is + 4] = sl[i].iT;
		col_ind[is + 5] = sl[i].iB;
		col_ind[is + 6] = sl[i].iP;
	}

	if (iVar == TEMP) {
		delete[] sl;
	}



	const integer size0 = maxelm + maxbound;

	// 16 векторов это 2.3 матрицы СЛАУ.
	// 7 векторов дополнительная матрица для кеша.
	// Всего 3.3 матрицы. 23 вектора.


	// Вектора необходимые для работы предобусловленного метода сопряжённых градиентов.
   // Вектора необходимые для работы BiCGStab.

	doublereal* residual = nullptr;
	doublereal* x_new = nullptr;

	doublereal* ri75 = nullptr;
	doublereal* roc75 = nullptr;
	doublereal* s75 = nullptr;
	doublereal* t75 = nullptr;
	doublereal* vec75 = nullptr;
	doublereal* vi75 = nullptr;
	doublereal* pi75 = nullptr;
	doublereal* dx75 = nullptr;

	doublereal* y75 = nullptr;
	doublereal* z75 = nullptr;
	// Первое предобуславливание:
	doublereal* pi76 = nullptr;

	// Второе предобуславливание:
	doublereal* z76 = nullptr;
	doublereal* s76 = nullptr;
	doublereal* z77 = nullptr;

	residual = new doublereal[size0];
	x_new = new doublereal[size0];
	
	pi76 = new doublereal[size0];
	ri75 = new doublereal[size0];
	roc75 = new doublereal[size0];
	s75 = new doublereal[size0];
	t75 = new doublereal[size0];
	vec75 = new doublereal[size0];
	vi75 = new doublereal[size0];
	pi75 = new doublereal[size0];
	dx75 = new doublereal[size0];
	y75 = new doublereal[size0];
	z75 = new doublereal[size0];
	z77 = new doublereal[size0];
	z76 = new doublereal[size0];
	s76 = new doublereal[size0];
	z77 = new doublereal[size0];

#pragma omp parallel for
	for (integer i = 0; i < size0; ++i)
	{

		residual[i] = 0.0;
		x_new[i] = 0.0;

		s75[i] = 0.0;
		t75[i] = 0.0;

		vi75[i] = 0.0;
		pi75[i] = 0.0;

		y75[i] = 0.0;
		z75[i] = 0.0;

		// Начальное приближение.
		dx75[i] = x[i];
	}

	
	doublereal delta075 = 1.0e30;// deltai75 = 1.0e30;
	doublereal bet75 = 0.0, roi75 = 0.0;
	doublereal roim175 = 1.0, al75 = 1.0, wi75 = 1.0;

	{
		// Диапазон изменения коэффициентов диффузии в уравнении на поправку давления.
		// опираемся на опыт теплопередачи.
		doublereal c1 = 0.026; // воздух
		doublereal c2 = 2000; // алмаз

		// см. Николаев
		if (bprintMessage) {
			std::cout << " number equations " << maxelm + maxbound << std::endl;
		}
		kend = static_cast<integer>(0.32 * sqrt(c2 / c1) * sqrt(size0) * log(2.0 / 1.0e-4) / degree);
		if (bprintMessage) {
			std::cout << "number out iteration = " << kend << "  " << kend * degree << std::endl;
		}
	}


	doublereal dmax = 1.0, dmax2 = 0.0, delta = 0.0;
	integer iout = 0;


	bool bdivergence_chebyshev = false;

#if MY_GPU_MATH_LIBRARY_CU_ON

	if (init_b_first_device) {

		int device;
		device = idevice_Tesla;

		// Choose which GPU to run on, change this on a multi-GPU system.
		cudaStatus = cudaSetDevice(device);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?");
			system("PAUSE");
			exit(1);
		}

		cudaDeviceProp prop;
		cudaGetDeviceProperties(&prop, device);
		printf("%s\n", prop.name);
		// GeForce 840M 384потока 1ГГц каждый март 2014 года. 28нм.

		init_b_first_device = false;
	}

	doublereal* dev_data = nullptr;
	int* dev_col_ind = nullptr;

	doublereal* dev_data_b = nullptr;
	int* dev_col_ind_b = nullptr;


	doublereal* dev_x_new = nullptr;
	doublereal* dev_x = nullptr;
	doublereal* dev_residual = nullptr;
	doublereal* dev_b = nullptr;
	doublereal* sum = (doublereal*)malloc(blocksPerGrid(size0) * sizeof(doublereal));
	doublereal* dev_sum;

	// Allocate GPU buffers for three vectors (two input, one output)    .
	cudaStatus = cudaMalloc((void**)&dev_sum, blocksPerGrid(size0) * sizeof(doublereal));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc dev_sum failed!");
		system("PAUSE");
		exit(1);
	}

	// Allocate GPU buffers for three vectors (two input, one output)    .
	cudaStatus = cudaMalloc((void**)&dev_data, (maxelm * 7) * sizeof(doublereal));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc dev_data failed!");
		system("PAUSE");
		exit(1);
	}

	// Allocate GPU buffers for three vectors (two input, one output)    .
	cudaStatus = cudaMalloc((void**)&dev_col_ind, (maxelm * 7) * sizeof(int));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc dev_col_ind failed!");
		system("PAUSE");
		exit(1);
	}

	// Allocate GPU buffers for three vectors (two input, one output)    .
	cudaStatus = cudaMalloc((void**)&dev_data_b, (maxbound * 2) * sizeof(doublereal));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc dev_data_b failed!");
		system("PAUSE");
		exit(1);
	}

	// Allocate GPU buffers for three vectors (two input, one output)    .
	cudaStatus = cudaMalloc((void**)&dev_col_ind_b, (maxbound * 2) * sizeof(int));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc dev_col_ind_b failed!");
		system("PAUSE");
		exit(1);
	}

	// Allocate GPU buffers for three vectors (two input, one output)    .
	cudaStatus = cudaMalloc((void**)&dev_x_new, (size0) * sizeof(doublereal));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc dev_xnew failed!");
		system("PAUSE");
		exit(1);
	}

	// Allocate GPU buffers for three vectors (two input, one output)    .
	cudaStatus = cudaMalloc((void**)&dev_b, (size0) * sizeof(doublereal));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc dev_b failed!");
		system("PAUSE");
		exit(1);
	}

	// Allocate GPU buffers for three vectors (two input, one output)    .
	cudaStatus = cudaMalloc((void**)&dev_x, (size0) * sizeof(doublereal));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc dev_x failed!");
		system("PAUSE");
		exit(1);
	}

	// Allocate GPU buffers for three vectors (two input, one output)    .
	cudaStatus = cudaMalloc((void**)&dev_residual, (size0) * sizeof(doublereal));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc dev_residual failed!");
		system("PAUSE");
		exit(1);
	}

	// Copy input vectors from host memory to GPU buffers.
	cudaStatus = cudaMemcpy(dev_data, data, maxelm * 7 * sizeof(doublereal), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy dev_data, val HostToDevice failed!");
		system("PAUSE");
		exit(1);
	}

	// Copy input vectors from host memory to GPU buffers.
	cudaStatus = cudaMemcpy(dev_col_ind, col_ind, maxelm * 7 * sizeof(int), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy dev_col_ind, val HostToDevice failed!");
		system("PAUSE");
		exit(1);
	}

	// Copy input vectors from host memory to GPU buffers.
	cudaStatus = cudaMemcpy(dev_data_b, data_b, maxbound * 2 * sizeof(doublereal), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy dev_data_b, val HostToDevice failed!");
		system("PAUSE");
		exit(1);
	}

	// Copy input vectors from host memory to GPU buffers.
	cudaStatus = cudaMemcpy(dev_col_ind_b, col_ind_b, maxbound * 2 * sizeof(int), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy dev_col_ind_b, val HostToDevice failed!");
		system("PAUSE");
		exit(1);
	}

	// Copy input vectors from host memory to GPU buffers.
	cudaStatus = cudaMemcpy(dev_x_new, x_new, size0 * sizeof(doublereal), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy dev_xnew, val HostToDevice failed!");
		system("PAUSE");
		exit(1);
	}


	// Copy input vectors from host memory to GPU buffers.
	//cudaStatus = cudaMemcpy(dev_b, b, size0 * sizeof(doublereal), cudaMemcpyHostToDevice);
	//if (cudaStatus != cudaSuccess) {
		//fprintf(stderr, "cudaMemcpy dev_b val HostToDevice failed!");
		//system("PAUSE");
		//exit(1);
	//}


	// Copy input vectors from host memory to GPU buffers.
	//cudaStatus = cudaMemcpy(dev_x, x, size0 * sizeof(doublereal), cudaMemcpyHostToDevice);
	//if (cudaStatus != cudaSuccess) {
		//fprintf(stderr, "cudaMemcpy dev_x, val HostToDevice failed!");
		//system("PAUSE");
		//exit(1);
	//}

	// Copy input vectors from host memory to GPU buffers.
	cudaStatus = cudaMemcpy(dev_residual, residual, size0 * sizeof(doublereal), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy dev_residual, val HostToDevice failed!");
		system("PAUSE");
		exit(1);
	}


	// результат записан в r75.
	// r75 - вектор невязки.
	//residual_cg_cuda(
		//maxelm, maxbound,
		//dev_data, dev_col_ind,
		//dev_data_b, dev_col_ind_b,
		//x, b, ri75,
		//dev_x, dev_b, dev_residual);


	// Copy input vectors from host memory to GPU buffers.
	cudaStatus = cudaMemcpy(dev_x, dx75, size0 * sizeof(doublereal), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy dev_x, x HostToDevice failed!");
		system("PAUSE");
		exit(1);
	}

	// Copy input vectors from host memory to GPU buffers.
	cudaStatus = cudaMemcpy(dev_b, b, size0 * sizeof(doublereal), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy dev_b, b HostToDevice failed!");
		system("PAUSE");
		exit(1);
	}

	Chebyshev_kernel_internal_residual2 << <128, 128 >> > (dev_data, dev_col_ind, dev_b, dev_x, dev_residual, maxelm);

	// Check for any errors launching the kernel
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "Chebyshev_kernel_internal_residual2 launch failed: %s\n", cudaGetErrorString(cudaStatus));
		system("PAUSE");
		exit(1);
	}

	Chebyshev_kernel_bound_residual2 << <128, 128 >> > (dev_data_b, dev_col_ind_b, dev_b, dev_x, dev_residual, maxbound);

	// Check for any errors launching the kernel
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "Chebyshev_kernel_bound_residual2 launch failed: %s\n", cudaGetErrorString(cudaStatus));
		system("PAUSE");
		exit(1);
	}

	cudaStatus = cudaMemcpy(residual, dev_residual, size0 * sizeof(doublereal), cudaMemcpyDeviceToHost);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "1. cudaMemcpy residual, dev_residual DeviceToHost failed!");
		system("PAUSE");
		exit(1);
	}

#pragma omp parallel for 
	for (integer i75 = 0; i75 < size0; ++i75) {

		ri75[i75] = residual[i75];

		roc75[i75] = 1.0;

		s76[i75] = ri75[i75];
		// r75[i75] = s76[i75 + 1];
	}

	delta075 = NormaV(ri75, size0);
	if (iVar == TEMP) {
		std::wcout << delta075 << std::endl;
	}

	// Используется для досрочного прерывания вычислительного процесса
		// как в алгоритме FGMRES Юсефа Саада и Мартина Г. Шульца.
	doublereal norma_b = NormaV_for_gmres(b, size0);

	doublereal dres00 = (NormaV_for_gmres(ri75, size0) / norma_b);

	integer iprecond = 0;


	

	iout = 0;

	bool bmem = true;
	doublereal dres0;



	do {

		

		roi75 = Scal(roc75, ri75, size0);
		bet75 = (roi75 / roim175) * (al75 / wi75);

#pragma omp parallel for 
		for (integer i75 = 0; i75 < size0; ++i75) {
			doublereal pibuf75 = ri75[i75] + (pi75[i75] - vi75[i75] * wi75) * bet75;
			pi75[i75] = pibuf75;
		}


		// Первое предобуславливание.
		// Ky=pi
#pragma omp parallel for 
		for (integer i75 = 0; i75 < size0; ++i75) {
			y75[i75] = 0.0; // Если начинать не с нуля то небудет сходимости для PAM !.

		   //y75[i75] = (rand() % 2 == 0 ? -1 : 1) * (1.0e-7*(rand()%1000));


			z77[i75] = 0.0;
			pi76[i75] = pi75[i75];
		}

		if (iprecond == 0)
		{

			doublereal r_start;

			// do {



			const doublereal tau0 = 2.0 / (hi + lo);
			const doublereal ksi = lo / hi;
			const doublereal rho0 = (1.0 - ksi) / (1.0 + ksi);


			//precontidioner
			// A, b, x -> z77
			// A*z77=pi76;
			//precontidioner_cuda(degree,
				//tau0,
				//rho0,
				//mu,
				//maxelm, maxbound,
				//dmax, dmax2, delta, r_start,
				//bprintMessage, bsolidT, j,
				//iout, iVar, size0,
				//bmem, eps,
				//dev_data, dev_col_ind,
				//dev_data_b, dev_col_ind_b, y75, pi76, z77,
				//residual, dev_x, dev_x_new, dev_residual, dev_b);

			{

				// Copy input vectors from host memory to GPU buffers.
				cudaStatus = cudaMemcpy(dev_b, pi76, size0 * sizeof(doublereal), cudaMemcpyHostToDevice);
				if (cudaStatus != cudaSuccess) {
					fprintf(stderr, "cudaMemcpy dev_b val HostToDevice failed!");
					system("PAUSE");
					exit(1);
				}

				// Copy input vectors from host memory to GPU buffers.
				cudaStatus = cudaMemcpy(dev_x, y75, size0 * sizeof(doublereal), cudaMemcpyHostToDevice);
				if (cudaStatus != cudaSuccess) {
					fprintf(stderr, "cudaMemcpy dev_b val HostToDevice failed!");
					system("PAUSE");
					exit(1);
				}

				// Copy input vectors from host memory to GPU buffers.
				cudaStatus = cudaMemcpy(dev_x_new, x_new, size0 * sizeof(doublereal), cudaMemcpyHostToDevice);
				if (cudaStatus != cudaSuccess) {
					fprintf(stderr, "cudaMemcpy dev_b val HostToDevice failed!");
					system("PAUSE");
					exit(1);
				}

				// std::cout << "degree="<< degree << std::endl;

				for (int k = 0; k < degree; ++k) {



					const doublereal rURF = tau0 / (1.0 + rho0 * mu[k]);

					dmax = 0.0;
					dmax2 = 0.0;
					delta = 0.0;

					doublereal dmax_loc = 0.0, dmax2_loc = 0.0;

					if ((k == 0) || (k == degree - 1)) {

						Chebyshev_kernel_internal2 << <128, 128 >> > (dev_data, dev_col_ind, dev_b,
							dev_x_new, dev_x, dev_residual, rURF, maxelm);

						// Check for any errors launching the kernel
						cudaStatus = cudaGetLastError();
						if (cudaStatus != cudaSuccess) {
							fprintf(stderr, "Jacoby_kernel_internal launch failed: %s\n", cudaGetErrorString(cudaStatus));
							system("PAUSE");
							exit(1);
						}

						//cudaStatus = cudaMemcpy(residual, dev_residual, size0 * sizeof(doublereal), cudaMemcpyDeviceToHost);
						//if (cudaStatus != cudaSuccess) {
							//fprintf(stderr, "cudaMemcpy residual, dev_residual DeviceToHost failed!");
							//system("PAUSE");
							//exit(1);
						//}

			//#pragma omp parallel for  reduction(+: dmax_loc, dmax2_loc)
				///			for (int i = 0; i < maxelm; ++i)
					//		{
						//		const int is = i * 7;
							//	const int is6 = is + 6;
								//const int it6 = col_ind[is6];

								//dmax_loc += residual[it6];
								//dmax2_loc += residual[it6] * residual[it6];

							//}

							//dmax += dmax_loc;
							//dmax2 += dmax2_loc;

							//dmax_loc = 0.0;
							//dmax2_loc = 0.0;

					}
					else {
						Chebyshev_kernel_internal << <128, 128 >> > (dev_data, dev_col_ind, dev_b,
							dev_x_new, dev_x, rURF, maxelm);

						// Check for any errors launching the kernel
						cudaStatus = cudaGetLastError();
						if (cudaStatus != cudaSuccess) {
							fprintf(stderr, "Jacoby_kernel_internal launch failed: %s\n", cudaGetErrorString(cudaStatus));
							system("PAUSE");
							exit(1);
						}
					}


					if ((k == 0) || (k == degree - 1)) {

						Chebyshev_kernel_bound2 << <128, 128 >> > (dev_data_b, dev_col_ind_b, dev_b, dev_x_new, dev_x, dev_residual, rURF, maxbound);

						// Check for any errors launching the kernel
						cudaStatus = cudaGetLastError();
						if (cudaStatus != cudaSuccess) {
							fprintf(stderr, "Jacoby_kernel_internal launch failed: %s\n", cudaGetErrorString(cudaStatus));
							system("PAUSE");
							exit(1);
						}

						//cudaStatus = cudaMemcpy(residual, dev_residual, size0 * sizeof(doublereal), cudaMemcpyDeviceToHost);
						//if (cudaStatus != cudaSuccess) {
							//fprintf(stderr, "2. cudaMemcpy residual, dev_residual DeviceToHost failed! %d",k);
							//system("PAUSE");
							//exit(1);
						//}

						//#pragma omp parallel for reduction(+: dmax_loc, dmax2_loc)
							//			for (int i = 0; i < maxbound; ++i) {


								//			const int is = i * 2;
									//		const int it1 = col_ind_b[is];//iI
										//	const int it2 = col_ind_b[is + 1];//iW                  


											// Чебышева
											//dmax_loc += (residual[it2]); // сумма модулей.
											//dmax2_loc += (residual[it2]) * (residual[it2]);

										//}

						//dmax_loc = thrust::reduce(residual, residual + size0);


						reduction_sum_kernel << <blocksPerGrid(size0), threadsPerBlock >> > (dev_residual, size0, dev_sum);

						// Check for any errors launching the kernel
						cudaStatus = cudaGetLastError();
						if (cudaStatus != cudaSuccess) {
							fprintf(stderr, "reduction_sum_kernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
							system("PAUSE");
							exit(1);
						}

						cudaStatus = cudaMemcpy(sum, dev_sum, blocksPerGrid(size0) * sizeof(doublereal), cudaMemcpyDeviceToHost);
						if (cudaStatus != cudaSuccess) {
							fprintf(stderr, "cudaMemcpy residual, dev_residual DeviceToHost failed!");
							system("PAUSE");
							exit(1);
						}

						dmax_loc = 0.0;
						const int isize_cycle = static_cast<int>(blocksPerGrid(size0));

#pragma omp parallel for reduction(+: dmax_loc)
						for (int i89 = 0; i89 < isize_cycle; ++i89) {
							dmax_loc += sum[i89];
						}


						sqr_vec_kernel << <128, 128 >> > (dev_residual, size0);

						// Check for any errors launching the kernel
						cudaStatus = cudaGetLastError();
						if (cudaStatus != cudaSuccess) {
							fprintf(stderr, "sqr_vec_kernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
							system("PAUSE");
							exit(1);
						}

						//cudaStatus = cudaMemcpy(residual, dev_residual, size0 * sizeof(doublereal), cudaMemcpyDeviceToHost);
						//if (cudaStatus != cudaSuccess) {
							//fprintf(stderr, "3. cudaMemcpy residual, dev_residual DeviceToHost failed!");
							//system("PAUSE");
							//exit(1);
						//}

						//dmax2_loc = thrust::reduce(residual, residual + size0);

						reduction_sum_kernel << <blocksPerGrid(size0), threadsPerBlock >> > (dev_residual, size0, dev_sum);

						// Check for any errors launching the kernel
						cudaStatus = cudaGetLastError();
						if (cudaStatus != cudaSuccess) {
							fprintf(stderr, "reduction_sum_kernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
							system("PAUSE");
							exit(1);
						}

						cudaStatus = cudaMemcpy(sum, dev_sum, blocksPerGrid(size0) * sizeof(doublereal), cudaMemcpyDeviceToHost);
						if (cudaStatus != cudaSuccess) {
							fprintf(stderr, "cudaMemcpy residual, dev_residual DeviceToHost failed!");
							system("PAUSE");
							exit(1);
						}

						dmax2_loc = 0.0;
						

#pragma omp parallel for reduction(+: dmax2_loc)
						for (int i89 = 0; i89 < isize_cycle; ++i89) {
							dmax2_loc += sum[i89];
						}


					}
					else {

						Chebyshev_kernel_bound << <128, 128 >> > (dev_data_b, dev_col_ind_b, dev_b, dev_x_new, dev_x, rURF, maxbound);

						// Check for any errors launching the kernel
						cudaStatus = cudaGetLastError();
						if (cudaStatus != cudaSuccess) {
							fprintf(stderr, "Jacoby_kernel_internal launch failed: %s\n", cudaGetErrorString(cudaStatus));
							system("PAUSE");
							exit(1);
						}
					}

					dmax += dmax_loc;
					dmax2 += dmax2_loc;



					//*std::min(1.0, 17000.0 / size0) 
					if (bsolidT) {
						if (iout == 0) eps = 1.0e-7 * dmax; // Важнейшее условие выхода из итерационного процесса.
					}
					else {
						if (iout == 0) eps = 0.1 * dmax; // Важнейшее условие выхода из итерационного процесса.
					}

					if (bprintMessage)
					{
						if ((iVar == PAM) || (iVar == TEMP))
						{
							if ((k == 0) || (k == degree - 1))
							{
								//dmax/=maxelm;
								if (j % degree == 0)
								{
#if doubleintprecision == 1
									printf("%lld %lld %e ", j + 1, iout + 1, dmax);
#else
									printf("%d %d %e ", j + 1, iout + 1, dmax);
#endif

								}
							}
						}
					}

					copy_Jacoby_kernel << <128, 128 >> > (dev_x_new, dev_x, size0);

					// Check for any errors launching the kernel
					cudaStatus = cudaGetLastError();
					if (cudaStatus != cudaSuccess) {
						fprintf(stderr, "copy_Jacoby_kernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
						system("PAUSE");
						exit(1);
					}

					if ((bmem) && (k == 0)) r_start = sqrt(dmax2 / (size0));
					if (k == degree - 1) delta = sqrt(dmax2 / (size0));

					if (bprintMessage) {
						if ((k == 0) || (k == degree - 1)) {
							if (j % degree == 0) {
								std::cout << std::endl;
							}
						}
					}

					j++;
				}// degree chebyshev


				cudaStatus = cudaMemcpy(z77, dev_x, size0 * sizeof(doublereal), cudaMemcpyDeviceToHost);
				if (cudaStatus != cudaSuccess) {
					fprintf(stderr, "cudaMemcpy x, dev_x DeviceToHost failed!");
					system("PAUSE");
					exit(1);
				}

				//std::wcout << "finish degree:"<< size0 << r_start << "  " << delta << "\n";
				//getchar();

			}



			delta /= r_start;

			// Корректировка нижней границы спектра lo линейного оператора А.
			// Жуков, Феодоритова. 13.02.2022
			correct_lo(lo, bmem, degree, hi, delta, bprintMessage);

			// std::wcout << dmax << std::endl;

		 //} while (dmax >= 0.0001);
		 //std::cout << "lo = " << lo << std::endl;

		 //getchar();
		}
		else {

			//residual_cg_cuda(
				//maxelm, maxbound,
				//dev_data, dev_col_ind,
				//dev_data_b, dev_col_ind_b,
				//y75, pi76, z77,
				//dev_x, dev_b, dev_residual);


			// Copy input vectors from host memory to GPU buffers.
			cudaStatus = cudaMemcpy(dev_x, y75, size0 * sizeof(doublereal), cudaMemcpyHostToDevice);
			if (cudaStatus != cudaSuccess) {
				fprintf(stderr, "cudaMemcpy dev_x, x HostToDevice failed!");
				system("PAUSE");
				exit(1);
			}

			// Copy input vectors from host memory to GPU buffers.
			cudaStatus = cudaMemcpy(dev_b, pi76, size0 * sizeof(doublereal), cudaMemcpyHostToDevice);
			if (cudaStatus != cudaSuccess) {
				fprintf(stderr, "cudaMemcpy dev_b, b HostToDevice failed!");
				system("PAUSE");
				exit(1);
			}

			Chebyshev_kernel_internal_residual2 << <128, 128 >> > (dev_data, dev_col_ind, dev_b, dev_x, dev_residual, maxelm);

			// Check for any errors launching the kernel
			cudaStatus = cudaGetLastError();
			if (cudaStatus != cudaSuccess) {
				fprintf(stderr, "Chebyshev_kernel_internal_residual2 launch failed: %s\n", cudaGetErrorString(cudaStatus));
				system("PAUSE");
				exit(1);
			}

			Chebyshev_kernel_bound_residual2 << <128, 128 >> > (dev_data_b, dev_col_ind_b, dev_b, dev_x, dev_residual, maxbound);

			// Check for any errors launching the kernel
			cudaStatus = cudaGetLastError();
			if (cudaStatus != cudaSuccess) {
				fprintf(stderr, "Chebyshev_kernel_bound_residual2 launch failed: %s\n", cudaGetErrorString(cudaStatus));
				system("PAUSE");
				exit(1);
			}

			cudaStatus = cudaMemcpy(z77, dev_residual, size0 * sizeof(doublereal), cudaMemcpyDeviceToHost);
			if (cudaStatus != cudaSuccess) {
				fprintf(stderr, "4. cudaMemcpy residual, dev_residual DeviceToHost failed!");
				system("PAUSE");
				exit(1);
			}

		}

#pragma omp parallel for 
		for (integer i75 = 0; i75 < size0; ++i75) {

			// Возвращаем результат.            
			y75[i75] = z77[i75];
			z77[i75] = 0.0;
		}



		// vi==A*y;
		// нумерация  y75 и vi75 начинается с нуля.
		//Matrix_by_vector_for_Cheb  
		//residual_cg_cuda(
			//maxelm, maxbound,
			//dev_data, dev_col_ind,
			//dev_data_b, dev_col_ind_b,
			//y75, z77, vi75,
			//dev_x, dev_b, dev_residual);


		// Copy input vectors from host memory to GPU buffers.
		cudaStatus = cudaMemcpy(dev_x, y75, size0 * sizeof(doublereal), cudaMemcpyHostToDevice);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMemcpy dev_x, x HostToDevice failed!");
			system("PAUSE");
			exit(1);
		}

		// Copy input vectors from host memory to GPU buffers.
		cudaStatus = cudaMemcpy(dev_b, z77, size0 * sizeof(doublereal), cudaMemcpyHostToDevice);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMemcpy dev_b, b HostToDevice failed!");
			system("PAUSE");
			exit(1);
		}

		Chebyshev_kernel_internal_residual2 << <128, 128 >> > (dev_data, dev_col_ind, dev_b, dev_x, dev_residual, maxelm);

		// Check for any errors launching the kernel
		cudaStatus = cudaGetLastError();
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "Chebyshev_kernel_internal_residual2 launch failed: %s\n", cudaGetErrorString(cudaStatus));
			system("PAUSE");
			exit(1);
		}

		Chebyshev_kernel_bound_residual2 << <128, 128 >> > (dev_data_b, dev_col_ind_b, dev_b, dev_x, dev_residual, maxbound);

		// Check for any errors launching the kernel
		cudaStatus = cudaGetLastError();
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "Chebyshev_kernel_bound_residual2 launch failed: %s\n", cudaGetErrorString(cudaStatus));
			system("PAUSE");
			exit(1);
		}

		cudaStatus = cudaMemcpy(vi75, dev_residual, size0 * sizeof(doublereal), cudaMemcpyDeviceToHost);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "5. cudaMemcpy residual, dev_residual DeviceToHost failed!");
			system("PAUSE");
			exit(1);
		}



#pragma omp parallel for 
		for (integer i75 = 0; i75 < size0; ++i75) {
			vi75[i75] *= -1;// результат занесён в  dax75
		}


		doublereal scal_val1 = Scal(roc75, vi75, size0);

		if ((fabs(roi75) < 1e-30) && (fabs(scal_val1) < 1e-30)) {
			al75 = 1.0;
		}
		else if (fabs(roi75) < 1e-30) {
			al75 = 0.0;
		}
		else {
			al75 = roi75 / scal_val1;
		}

#pragma omp parallel for
		for (integer i75 = 0; i75 < size0; ++i75)
		{
			s75[i75] = ri75[i75] - al75 * vi75[i75];
			z75[i75] = 0.0;
		}

		// Предобуславливание.
#pragma omp parallel for 
		for (integer i75 = 0; i75 < size0; ++i75)
		{
			z77[i75] = 0.0;// только ноль
			//z77[i75] = (rand() % 2 == 0 ? -1 : 1) * (1.0e-7 * (rand() % 1000));


			vec75[i75] = s75[i75];
			z76[i75] = 0.0;
			s76[i75] = s75[i75];
		}

		// A*z76=s76;

		if (iprecond == 0)
		{

			doublereal r_start;

			// do {



			const doublereal tau0 = 2.0 / (hi + lo);
			const doublereal ksi = lo / hi;
			const doublereal rho0 = (1.0 - ksi) / (1.0 + ksi);


			//std::wcout << Scal(s76, s76, size0) << std::endl;
		   // getchar();

			//precontidioner
			// A, b, x -> z76
			// A*z76=s76;
			//precontidioner_cuda(degree,
				//tau0,
				//rho0,
				//mu,
				//maxelm, maxbound,
				//dmax, dmax2, delta, r_start,
				//bprintMessage, bsolidT, j,
				//iout, iVar, size0,
				//bmem, eps,
				//dev_data, dev_col_ind,
				//dev_data_b, dev_col_ind_b, z77, s76, z76,
				//residual, dev_x, dev_x_new, dev_residual, dev_b);

			{


				// Copy input vectors from host memory to GPU buffers.
				cudaStatus = cudaMemcpy(dev_b, s76, size0 * sizeof(doublereal), cudaMemcpyHostToDevice);
				if (cudaStatus != cudaSuccess) {
					fprintf(stderr, "cudaMemcpy dev_b val HostToDevice failed!");
					system("PAUSE");
					exit(1);
				}

				// Copy input vectors from host memory to GPU buffers.
				cudaStatus = cudaMemcpy(dev_x, z77, size0 * sizeof(doublereal), cudaMemcpyHostToDevice);
				if (cudaStatus != cudaSuccess) {
					fprintf(stderr, "cudaMemcpy dev_b val HostToDevice failed!");
					system("PAUSE");
					exit(1);
				}

				// Copy input vectors from host memory to GPU buffers.
				cudaStatus = cudaMemcpy(dev_x_new, x_new, size0 * sizeof(doublereal), cudaMemcpyHostToDevice);
				if (cudaStatus != cudaSuccess) {
					fprintf(stderr, "cudaMemcpy dev_b val HostToDevice failed!");
					system("PAUSE");
					exit(1);
				}

				// std::cout << "degree="<< degree << std::endl;

				for (int k = 0; k < degree; ++k) {



					const doublereal rURF = tau0 / (1.0 + rho0 * mu[k]);

					dmax = 0.0;
					dmax2 = 0.0;
					delta = 0.0;

					doublereal dmax_loc = 0.0, dmax2_loc = 0.0;

					if ((k == 0) || (k == degree - 1)) {

						Chebyshev_kernel_internal2 << <128, 128 >> > (dev_data, dev_col_ind, dev_b,
							dev_x_new, dev_x, dev_residual, rURF, maxelm);

						// Check for any errors launching the kernel
						cudaStatus = cudaGetLastError();
						if (cudaStatus != cudaSuccess) {
							fprintf(stderr, "Jacoby_kernel_internal launch failed: %s\n", cudaGetErrorString(cudaStatus));
							system("PAUSE");
							exit(1);
						}

						//cudaStatus = cudaMemcpy(residual, dev_residual, size0 * sizeof(doublereal), cudaMemcpyDeviceToHost);
						//if (cudaStatus != cudaSuccess) {
							//fprintf(stderr, "cudaMemcpy residual, dev_residual DeviceToHost failed!");
							//system("PAUSE");
							//exit(1);
						//}

			//#pragma omp parallel for  reduction(+: dmax_loc, dmax2_loc)
				///			for (int i = 0; i < maxelm; ++i)
					//		{
						//		const int is = i * 7;
							//	const int is6 = is + 6;
								//const int it6 = col_ind[is6];

								//dmax_loc += residual[it6];
								//dmax2_loc += residual[it6] * residual[it6];

							//}

							//dmax += dmax_loc;
							//dmax2 += dmax2_loc;

							//dmax_loc = 0.0;
							//dmax2_loc = 0.0;

					}
					else {
						Chebyshev_kernel_internal << <128, 128 >> > (dev_data, dev_col_ind, dev_b,
							dev_x_new, dev_x, rURF, maxelm);

						// Check for any errors launching the kernel
						cudaStatus = cudaGetLastError();
						if (cudaStatus != cudaSuccess) {
							fprintf(stderr, "Jacoby_kernel_internal launch failed: %s\n", cudaGetErrorString(cudaStatus));
							system("PAUSE");
							exit(1);
						}
					}


					if ((k == 0) || (k == degree - 1)) {

						Chebyshev_kernel_bound2 << <128, 128 >> > (dev_data_b, dev_col_ind_b, dev_b, dev_x_new, dev_x, dev_residual, rURF, maxbound);

						// Check for any errors launching the kernel
						cudaStatus = cudaGetLastError();
						if (cudaStatus != cudaSuccess) {
							fprintf(stderr, "Jacoby_kernel_internal launch failed: %s\n", cudaGetErrorString(cudaStatus));
							system("PAUSE");
							exit(1);
						}

						//cudaStatus = cudaMemcpy(residual, dev_residual, size0 * sizeof(doublereal), cudaMemcpyDeviceToHost);
						//if (cudaStatus != cudaSuccess) {
							//fprintf(stderr, "6. cudaMemcpy residual, dev_residual DeviceToHost failed! %d",k);
							//system("PAUSE");
							//exit(1);
						//}

						//#pragma omp parallel for reduction(+: dmax_loc, dmax2_loc)
							//			for (int i = 0; i < maxbound; ++i) {


								//			const int is = i * 2;
									//		const int it1 = col_ind_b[is];//iI
										//	const int it2 = col_ind_b[is + 1];//iW                  


											// Чебышева
											//dmax_loc += (residual[it2]); // сумма модулей.
											//dmax2_loc += (residual[it2]) * (residual[it2]);

										//}

						//dmax_loc = thrust::reduce(residual, residual + size0);


						reduction_sum_kernel << <blocksPerGrid(size0), threadsPerBlock >> > (dev_residual, size0, dev_sum);

						// Check for any errors launching the kernel
						cudaStatus = cudaGetLastError();
						if (cudaStatus != cudaSuccess) {
							fprintf(stderr, "reduction_sum_kernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
							system("PAUSE");
							exit(1);
						}

						cudaStatus = cudaMemcpy(sum, dev_sum, blocksPerGrid(size0) * sizeof(doublereal), cudaMemcpyDeviceToHost);
						if (cudaStatus != cudaSuccess) {
							fprintf(stderr, "cudaMemcpy residual, dev_residual DeviceToHost failed!");
							system("PAUSE");
							exit(1);
						}

						dmax_loc = 0.0;
						const int isize_cycle = static_cast<int>(blocksPerGrid(size0));

#pragma omp parallel for reduction(+: dmax_loc)
						for (int i89 = 0; i89 < isize_cycle; ++i89) {
							dmax_loc += sum[i89];
						}


						sqr_vec_kernel << <128, 128 >> > (dev_residual, size0);

						// Check for any errors launching the kernel
						cudaStatus = cudaGetLastError();
						if (cudaStatus != cudaSuccess) {
							fprintf(stderr, "sqr_vec_kernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
							system("PAUSE");
							exit(1);
						}

						//cudaStatus = cudaMemcpy(residual, dev_residual, size0 * sizeof(doublereal), cudaMemcpyDeviceToHost);
						//if (cudaStatus != cudaSuccess) {
							//fprintf(stderr, "7. cudaMemcpy residual, dev_residual DeviceToHost failed! %d",k);
							//system("PAUSE");
							//exit(1);
						//}

						//dmax2_loc = thrust::reduce(residual, residual + size0);

						reduction_sum_kernel << <blocksPerGrid(size0), threadsPerBlock >> > (dev_residual, size0, dev_sum);

						// Check for any errors launching the kernel
						cudaStatus = cudaGetLastError();
						if (cudaStatus != cudaSuccess) {
							fprintf(stderr, "reduction_sum_kernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
							system("PAUSE");
							exit(1);
						}

						cudaStatus = cudaMemcpy(sum, dev_sum, blocksPerGrid(size0) * sizeof(doublereal), cudaMemcpyDeviceToHost);
						if (cudaStatus != cudaSuccess) {
							fprintf(stderr, "cudaMemcpy residual, dev_residual DeviceToHost failed!");
							system("PAUSE");
							exit(1);
						}

						dmax2_loc = 0.0;
						

#pragma omp parallel for reduction(+: dmax2_loc)
						for (int i89 = 0; i89 < isize_cycle; ++i89) {
							dmax2_loc += sum[i89];
						}

					}
					else {

						Chebyshev_kernel_bound << <128, 128 >> > (dev_data_b, dev_col_ind_b, dev_b, dev_x_new, dev_x, rURF, maxbound);

						// Check for any errors launching the kernel
						cudaStatus = cudaGetLastError();
						if (cudaStatus != cudaSuccess) {
							fprintf(stderr, "Jacoby_kernel_internal launch failed: %s\n", cudaGetErrorString(cudaStatus));
							system("PAUSE");
							exit(1);
						}
					}

					dmax += dmax_loc;
					dmax2 += dmax2_loc;



					//*std::min(1.0, 17000.0 / size0) 
					if (bsolidT) {
						if (iout == 0) eps = 1.0e-7 * dmax; // Важнейшее условие выхода из итерационного процесса.
					}
					else {
						if (iout == 0) eps = 0.1 * dmax; // Важнейшее условие выхода из итерационного процесса.
					}

					if (bprintMessage)
					{
						if ((iVar == PAM) || (iVar == TEMP))
						{
							if ((k == 0) || (k == degree - 1))
							{
								//dmax/=maxelm;
								if (j % degree == 0)
								{
#if doubleintprecision == 1
									printf("%lld %lld %e ", j + 1, iout + 1, dmax);
#else
									printf("%d %d %e ", j + 1, iout + 1, dmax);
#endif

								}
							}
						}
					}

					copy_Jacoby_kernel << <128, 128 >> > (dev_x_new, dev_x, size0);

					// Check for any errors launching the kernel
					cudaStatus = cudaGetLastError();
					if (cudaStatus != cudaSuccess) {
						fprintf(stderr, "copy_Jacoby_kernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
						system("PAUSE");
						exit(1);
					}

					if ((bmem) && (k == 0)) r_start = sqrt(dmax2 / (size0));
					if (k == degree - 1) delta = sqrt(dmax2 / (size0));

					if (bprintMessage) {
						if ((k == 0) || (k == degree - 1)) {
							if (j % degree == 0) {
								std::cout << std::endl;
							}
						}
					}

					j++;
				}// degree chebyshev


				cudaStatus = cudaMemcpy(z76, dev_x, size0 * sizeof(doublereal), cudaMemcpyDeviceToHost);
				if (cudaStatus != cudaSuccess) {
					fprintf(stderr, "cudaMemcpy x, dev_x DeviceToHost failed!");
					system("PAUSE");
					exit(1);
				}

				//std::wcout << "finish degree:"<< size0 << r_start << "  " << delta << "\n";
				//getchar();
			}

			delta /= r_start;

			// Корректировка нижней границы спектра lo линейного оператора А.
			// Жуков, Феодоритова. 13.02.2022
			correct_lo(lo, bmem, degree, hi, delta, bprintMessage);


			// std::wcout << dmax << std::endl;

			// } while (dmax >= 0.0001);

			 //getchar();

			 //std::cout << "lo = " << lo << std::endl;
		}
		else {

			//residual_cg_cuda(
				//maxelm, maxbound,
				//dev_data, dev_col_ind,
				//dev_data_b, dev_col_ind_b,
				//z77, s76, z76,
				//dev_x, dev_b, dev_residual);

			// Copy input vectors from host memory to GPU buffers.
			cudaStatus = cudaMemcpy(dev_x, z77, size0 * sizeof(doublereal), cudaMemcpyHostToDevice);
			if (cudaStatus != cudaSuccess) {
				fprintf(stderr, "cudaMemcpy dev_x, x HostToDevice failed!");
				system("PAUSE");
				exit(1);
			}

			// Copy input vectors from host memory to GPU buffers.
			cudaStatus = cudaMemcpy(dev_b, s76, size0 * sizeof(doublereal), cudaMemcpyHostToDevice);
			if (cudaStatus != cudaSuccess) {
				fprintf(stderr, "cudaMemcpy dev_b, b HostToDevice failed!");
				system("PAUSE");
				exit(1);
			}

			Chebyshev_kernel_internal_residual2 << <128, 128 >> > (dev_data, dev_col_ind, dev_b, dev_x, dev_residual, maxelm);

			// Check for any errors launching the kernel
			cudaStatus = cudaGetLastError();
			if (cudaStatus != cudaSuccess) {
				fprintf(stderr, "Chebyshev_kernel_internal_residual2 launch failed: %s\n", cudaGetErrorString(cudaStatus));
				system("PAUSE");
				exit(1);
			}

			Chebyshev_kernel_bound_residual2 << <128, 128 >> > (dev_data_b, dev_col_ind_b, dev_b, dev_x, dev_residual, maxbound);

			// Check for any errors launching the kernel
			cudaStatus = cudaGetLastError();
			if (cudaStatus != cudaSuccess) {
				fprintf(stderr, "Chebyshev_kernel_bound_residual2 launch failed: %s\n", cudaGetErrorString(cudaStatus));
				system("PAUSE");
				exit(1);
			}

			cudaStatus = cudaMemcpy(z76, dev_residual, size0 * sizeof(doublereal), cudaMemcpyDeviceToHost);
			if (cudaStatus != cudaSuccess) {
				fprintf(stderr, "8. cudaMemcpy residual, dev_residual DeviceToHost failed!");
				system("PAUSE");
				exit(1);
			}


		}

#pragma omp parallel for 
		for (integer i75 = 0; i75 < size0; ++i75) {

			z77[i75] = 0.0;
			// Возвращаем результат.
			s75[i75] = vec75[i75];
			// Возвращаем результат.
			z75[i75] = z76[i75];
		}

		// t==A*z;
		 // нумерация  y75 и vi75 начинается с нуля.
		//Matrix_by_vector_for_Cheb  
		//residual_cg_cuda(
			//maxelm, maxbound,
			//dev_data, dev_col_ind,
			//dev_data_b, dev_col_ind_b,
			//z75, z77, t75,
			//dev_x, dev_b, dev_residual);

		// Copy input vectors from host memory to GPU buffers.
		cudaStatus = cudaMemcpy(dev_x, z75, size0 * sizeof(doublereal), cudaMemcpyHostToDevice);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMemcpy dev_x, x HostToDevice failed!");
			system("PAUSE");
			exit(1);
		}

		// Copy input vectors from host memory to GPU buffers.
		cudaStatus = cudaMemcpy(dev_b, z77, size0 * sizeof(doublereal), cudaMemcpyHostToDevice);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMemcpy dev_b, b HostToDevice failed!");
			system("PAUSE");
			exit(1);
		}

		Chebyshev_kernel_internal_residual2 << <128, 128 >> > (dev_data, dev_col_ind, dev_b, dev_x, dev_residual, maxelm);

		// Check for any errors launching the kernel
		cudaStatus = cudaGetLastError();
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "Chebyshev_kernel_internal_residual2 launch failed: %s\n", cudaGetErrorString(cudaStatus));
			system("PAUSE");
			exit(1);
		}

		Chebyshev_kernel_bound_residual2 << <128, 128 >> > (dev_data_b, dev_col_ind_b, dev_b, dev_x, dev_residual, maxbound);

		// Check for any errors launching the kernel
		cudaStatus = cudaGetLastError();
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "Chebyshev_kernel_bound_residual2 launch failed: %s\n", cudaGetErrorString(cudaStatus));
			system("PAUSE");
			exit(1);
		}

		cudaStatus = cudaMemcpy(t75, dev_residual, size0 * sizeof(doublereal), cudaMemcpyDeviceToHost);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "9. cudaMemcpy residual, dev_residual DeviceToHost failed!");
			system("PAUSE");
			exit(1);
		}



#pragma omp parallel for 
		for (integer i75 = 0; i75 < size0; ++i75) {
			t75[i75] *= -1;// результат занесён в  dax75
		}


		wi75 = Scal(t75, s75, size0) / Scal(t75, t75, size0);
		// std::cout << "Scal(t75,s75,n75)==" << Scal(t75,s75,n75) << ", Scal(t75,t75,n75)=="<< Scal(t75,t75,n75) << std::endl;

#pragma omp parallel for 
		for (integer i75 = 0; i75 < size0; ++i75) {
			//dx75[i75]+=al75*pi75[i75]+wi75*s75[i75]; // так было без предобуславливателя
			dx75[i75] += al75 * y75[i75] + wi75 * z75[i75]; // так стало с предобуславливателем
			ri75[i75] = s75[i75] - wi75 * t75[i75];
		}
		//dres0 = deltai75 = NormaV(ri75, size0);

		//dres0 = (NormaV_for_gmres(ri75, size0) / norma_b)/ size0;
		dres0 = (NormaV_for_gmres(ri75, size0) / norma_b);
		//dres0 = (NormaV_for_gmres(ri75, size0));

		if (iVar == TEMP) {
			std::cout << iout << "  " << dres0 << std::endl;
		}

		if (iVar == PAM) {
			// CFD break for PAM.

			doublereal dres000 = std::min(dres_init_old_Simple_iteration, dres00);

			if (((dres0 / dres000 < 0.01) && ((size0 < 6000000) || (iout > 3)))) {

				if (lb < 23) {
					break;
				}
				else {

					if (iout > 6) break;
				}
			}

			if (iout > (size0 > 6000000 ? 15 : 5)) {
				if (dres0 / dres000 < (size0 > 6000000 ? 0.01 : 0.01)) {
					break;
				}
				else {
					if (iout > (size0 > 6000000 ? 21 : 7)) {
						std::cout << " PAM chebyshev degree = " << degree << " divergence!#it = " << (size0 > 6000000 ? 21 : 7) << "\n";
						//system("PAUSE");
						bdivergence_chebyshev = true;
						break;
					}
				}
			}
		}

		//getchar();

		++iout;
	} while (dres0 >= 1.0e-7);



	cudaFree(dev_data);
	cudaFree(dev_col_ind);

	cudaFree(dev_data_b);
	cudaFree(dev_col_ind_b);

	cudaFree(dev_b);
	cudaFree(dev_x_new);
	cudaFree(dev_x);
	cudaFree(dev_residual);

	cudaFree(dev_sum);
	free(sum);

#else 



	// результат записан в r75.
	// r75 - вектор невязки.
	residual_cg(
		maxelm, maxbound,
		data, col_ind,
		data_b, col_ind_b,
		x, b, ri75);

#pragma omp parallel for 
	for (integer i75 = 0; i75 < size0; ++i75) {

		roc75[i75] = 1.0;

		s76[i75] = ri75[i75];
		// r75[i75] = s76[i75 + 1];
	}

	delta075 = NormaV(ri75, size0);
	if (iVar == TEMP) {
		std::wcout << delta075 << std::endl;
	}

	// Используется для досрочного прерывания вычислительного процесса
		// как в алгоритме FGMRES Юсефа Саада и Мартина Г. Шульца.
	doublereal norma_b = NormaV_for_gmres(b, size0);

	doublereal dres00 = (NormaV_for_gmres(ri75, size0) / norma_b);

	integer iprecond = 0;


	

	iout = 0;

	bool bmem = true;
	doublereal dres0;

	

	do {


		roi75 = Scal(roc75, ri75, size0);
		bet75 = (roi75 / roim175) * (al75 / wi75);

#pragma omp parallel for 
		for (integer i75 = 0; i75 < size0; ++i75) {
			doublereal pibuf75 = ri75[i75] + (pi75[i75] - vi75[i75] * wi75) * bet75;
			pi75[i75] = pibuf75;
		}


		// Первое предобуславливание.
		// Ky=pi
#pragma omp parallel for 
		for (integer i75 = 0; i75 < size0; ++i75) {
			y75[i75] = 0.0; // Если начинать не с нуля то небудет сходимости для PAM !.

		   //y75[i75] = (rand() % 2 == 0 ? -1 : 1) * (1.0e-7*(rand()%1000));


			z77[i75] = 0.0;
			pi76[i75] = pi75[i75];
		}

		if (iprecond == 0)
		{

			doublereal r_start;

			// do {



			const doublereal tau0 = 2.0 / (hi + lo);
			const doublereal ksi = lo / hi;
			const doublereal rho0 = (1.0 - ksi) / (1.0 + ksi);


			//precontidioner
			// A, b, x -> z77
			// A*z77=pi76;
			precontidioner(degree,
				tau0,
				rho0,
				mu,
				maxelm, maxbound,
				dmax, dmax2, delta, r_start,
				bprintMessage, bsolidT, j,
				iout, iVar, size0,
				bmem, eps,
				data, col_ind,
				data_b, col_ind_b, y75, pi76, z77);



			delta /= r_start;

			// Корректировка нижней границы спектра lo линейного оператора А.
			// Жуков, Феодоритова. 13.02.2022
			correct_lo(lo, bmem, degree, hi, delta, bprintMessage);

			// std::wcout << dmax << std::endl;

		 //} while (dmax >= 0.0001);
		 //std::cout << "lo = " << lo << std::endl;

		 //getchar();
		}
		else {

			residual_cg(
				maxelm, maxbound,
				data, col_ind,
				data_b, col_ind_b,
				y75, pi76, z77);
		}

#pragma omp parallel for 
		for (integer i75 = 0; i75 < size0; ++i75) {

			// Возвращаем результат.            
			y75[i75] = z77[i75];
			z77[i75] = 0.0;
		}



		// vi==A*y;
		// нумерация  y75 и vi75 начинается с нуля.
		//Matrix_by_vector_for_Cheb  
		residual_cg(
			maxelm, maxbound,
			data, col_ind,
			data_b, col_ind_b,
			y75, z77, vi75);

#pragma omp parallel for 
		for (integer i75 = 0; i75 < size0; ++i75) {
			vi75[i75] *= -1;// результат занесён в  dax75
		}


		doublereal scal_val1 = Scal(roc75, vi75, size0);

		if ((fabs(roi75) < 1e-30) && (fabs(scal_val1) < 1e-30)) {
			al75 = 1.0;
		}
		else if (fabs(roi75) < 1e-30) {
			al75 = 0.0;
		}
		else {
			al75 = roi75 / scal_val1;
		}

#pragma omp parallel for
		for (integer i75 = 0; i75 < size0; ++i75)
		{
			s75[i75] = ri75[i75] - al75 * vi75[i75];
			z75[i75] = 0.0;
		}

		// Предобуславливание.
#pragma omp parallel for 
		for (integer i75 = 0; i75 < size0; ++i75)
		{
			z77[i75] = 0.0;// только ноль
			//z77[i75] = (rand() % 2 == 0 ? -1 : 1) * (1.0e-7 * (rand() % 1000));


			vec75[i75] = s75[i75];
			z76[i75] = 0.0;
			s76[i75] = s75[i75];
		}

		// A*z76=s76;

		if (iprecond == 0)
		{

			doublereal r_start;

			// do {



			const doublereal tau0 = 2.0 / (hi + lo);
			const doublereal ksi = lo / hi;
			const doublereal rho0 = (1.0 - ksi) / (1.0 + ksi);


			//std::wcout << Scal(s76, s76, size0) << std::endl;
		   // getchar();

			//precontidioner
			// A, b, x -> z76
			// A*z76=s76;
			precontidioner(degree,
				tau0,
				rho0,
				mu,
				maxelm, maxbound,
				dmax, dmax2, delta, r_start,
				bprintMessage, bsolidT, j,
				iout, iVar, size0,
				bmem, eps,
				data, col_ind,
				data_b, col_ind_b, z77, s76, z76);



			delta /= r_start;

			// Корректировка нижней границы спектра lo линейного оператора А.
			// Жуков, Феодоритова. 13.02.2022
			correct_lo(lo, bmem, degree, hi, delta, bprintMessage);


			// std::wcout << dmax << std::endl;

			// } while (dmax >= 0.0001);

			 //getchar();

			 //std::cout << "lo = " << lo << std::endl;
		}
		else {

			residual_cg(
				maxelm, maxbound,
				data, col_ind,
				data_b, col_ind_b,
				z77, s76, z76);
		}

#pragma omp parallel for 
		for (integer i75 = 0; i75 < size0; ++i75) {

			z77[i75] = 0.0;
			// Возвращаем результат.
			s75[i75] = vec75[i75];
			// Возвращаем результат.
			z75[i75] = z76[i75];
		}

		// t==A*z;
		 // нумерация  y75 и vi75 начинается с нуля.
		//Matrix_by_vector_for_Cheb  
		residual_cg(
			maxelm, maxbound,
			data, col_ind,
			data_b, col_ind_b,
			z75, z77, t75);

#pragma omp parallel for 
		for (integer i75 = 0; i75 < size0; ++i75) {
			t75[i75] *= -1;// результат занесён в  dax75
		}


		wi75 = Scal(t75, s75, size0) / Scal(t75, t75, size0);
		// std::cout << "Scal(t75,s75,n75)==" << Scal(t75,s75,n75) << ", Scal(t75,t75,n75)=="<< Scal(t75,t75,n75) << std::endl;

#pragma omp parallel for 
		for (integer i75 = 0; i75 < size0; ++i75) {
			//dx75[i75]+=al75*pi75[i75]+wi75*s75[i75]; // так было без предобуславливателя
			dx75[i75] += al75 * y75[i75] + wi75 * z75[i75]; // так стало с предобуславливателем
			ri75[i75] = s75[i75] - wi75 * t75[i75];
		}
		//dres0 = deltai75 = NormaV(ri75, size0);

		//dres0 = (NormaV_for_gmres(ri75, size0) / norma_b)/ size0;
		dres0 = (NormaV_for_gmres(ri75, size0) / norma_b);
		//dres0 = (NormaV_for_gmres(ri75, size0));

		if (iVar == TEMP) {
			std::cout << iout << "  " << dres0 << std::endl;
		}

		if (iVar == PAM) {
			// CFD break for PAM.

			doublereal dres000 = std::min(dres_init_old_Simple_iteration, dres00);

			
				if (((dres0 / dres000 < 0.01) && ((size0 < 6000000) || (iout > 3)))) {

					if (lb < 23) {
						break;
					}
					else {

						if (iout > 6) break;
					}
				}

				if (iout > ((size0 > 6000000) || (lb >= 23) ? 15 : 5)) {
					if (dres0 / dres000 < (size0 > 6000000 ? 0.01 : 0.01)) {
						break;
					}
					else {
						if (iout > (size0 > 6000000 ? 21 : 7)) {
							std::cout << " PAM chebyshev degree = " << degree << " divergence!#it = " << (size0 > 6000000 ? 21 : 7) << "\n";
							//system("PAUSE");
							bdivergence_chebyshev = true;
							break;
						}
					}
				}
			
		}

		//getchar();

		++iout;
	} while (dres0 >= 1.0e-7);


	

#endif


	if (dres00 <= dres_init_old_Simple_iteration)
	{
		dres_init_old_Simple_iteration = dres00;
	}

	delete[] x_new;
	delete[] residual;

	if (iVar == PAM) {
		std::cout << " " << iout << " ";
	}

	if (!bdivergence_chebyshev) {

#pragma omp parallel for
		for (integer i75 = 0; i75 < size0; ++i75) {

			x[i75] = dx75[i75];
		}



		for (int i = 0; i < maxelm; ++i) {
			const int is = i * 7;
			const int iP = col_ind[is + 6];
			x[iP] *= D[iP];
		}

#pragma omp parallel for 
		for (int i = 0; i < maxbound; ++i) {
			const int is = i * 2;
			const int iW = col_ind_b[is + 1];
			x[iW] *= D[iW];
		}
	}

	delete[] data;
	delete[] col_ind;

	delete[] data_b;
	delete[] col_ind_b;

	if (iVar == TEMP) {
		std::cout << " " << lo << " " << j << " ";
	}

	if (j > kend * degree - 100) {
		std::cout << "\n chebyshev divergence " << dmax << " ";
		system("pause");
	}


	// Освобождение оперативной памяти.
	// Первое предобуславливание

	delete[] pi76;
	pi76 = nullptr;


	// Второе предобуславливание

	delete[] z76;
	z76 = nullptr;
	delete[] s76;
	s76 = nullptr;

	delete[] ri75;
	ri75 = nullptr;
	delete[] roc75;
	roc75 = nullptr;

	delete[] s75;
	s75 = nullptr;
	delete[] t75;
	t75 = nullptr;

	delete[] vec75;
	vec75 = nullptr;
	delete[] vi75;
	vi75 = nullptr;


	delete[] pi75;
	pi75 = nullptr;
	delete[] dx75;
	dx75 = nullptr;


	delete[] y75;
	y75 = nullptr;

	delete[] z75;
	z75 = nullptr;

	delete[] z77;
	z77 = nullptr;

	if (iVar == TEMP) {
		// Реанимируем то что было выделено до решения СЛАУ,
		// а в момент решения СЛАУ мы экономили оперативную память.
		slb = new equation3D_bon[maxbound];
		sl = new equation3D[maxelm];
	}

	if (bprintMessage) {
		if ((iVar == PAM) || (iVar == TEMP)) {
			printf("calc complete...\n");
			system("pause");
		}
	}
	//printf("4000 %e \n", dmax);
	//system("pause");

	return bdivergence_chebyshev;

} // bicgstab_chebyshev_3Dnow_cuda



#endif