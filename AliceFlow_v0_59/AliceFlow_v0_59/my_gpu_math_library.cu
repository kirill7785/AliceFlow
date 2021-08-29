// my_gpu_math_library.cu
// Реализация операций с длинными векторами на видеокарте.
// 30.07.2021


#pragma once
#ifndef MY_GPU_MATH_LIBRARY_CU
#define MY_GPU_MATH_LIBRARY_CU 1

// false - используется центральный процессор.
#define MY_GPU_MATH_LIBRARY_CU_ON  false 
#define EllFormat  false

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

	if ((MY_GPU_MATH_LIBRARY_CU_ON)&&(bdevice)) {

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
		for (integer i = 0; i < n; i++) {
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
			//getchar();
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