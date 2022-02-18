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

#endif