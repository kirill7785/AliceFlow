
#pragma once
#ifndef BASIC_FUNCTIONS_MY_AGREGAT_AMG_RELAXATION_SPAI0_CPP
#define BASIC_FUNCTIONS_MY_AGREGAT_AMG_RELAXATION_SPAI0_CPP 1

// Seidel.
// n Ц число неизвестных,
// x Ц вектор решени€ (текущее приближение к решению)
// val, col_ind, row_ptr матрица в CRS формате.
// rthdsd Ц права€ часть —Ћј”.
// ћетод практически не сходитс€, или очень очень долго сходитс€.
// 07.01.2022
template <typename doublerealT>
void Seidel(integer n, doublerealT*& val,
	integer*& col_ind, integer*& row_ptr,
	doublerealT*& x, doublerealT*& rthdsd)
{
	const int nu = 1;

	// ¬нешние итерации.
	for (integer k = 0; k < nu; ++k) {

#pragma omp parallel for 
		for (integer i = 0; i < n; ++i) {
			doublereal diag = 0.0;
			doublereal sum = rthdsd[i];

			for (integer j = row_ptr[i]; j <= row_ptr[i + 1] - 1; ++j) {

				if (col_ind[j] == i) {
					diag = val[j]; // диагональ
				}
				else {
					sum -= val[j] * x[col_ind[j]];
				}
			}
			x[i] += 0.8 * ((sum) / diag - x[i]);
		}
	}
} // Seidel



// Sparse approximate inverse relaxation scheme.
// n Ц число неизвестных,
// x Ц вектор решени€ (текущее приближение к решению)
// val, col_ind, row_ptr матрица в CRS формате.
// rthdsd Ц права€ часть —Ћј”.
// ћетод практически не сходитс€, или очень очень долго сходитс€.
template <typename doublerealT>
void spai0_smoother(integer n, doublerealT*& val,
	integer*& col_ind, integer*& row_ptr,
	doublerealT*& x, doublerealT*& rthdsd)
{

#pragma omp parallel for 
	for (integer i = 0; i < n; ++i) {
		doublereal num = 0.0;
		doublereal den = 0.0;
		doublereal sum = 0.0;

		for (integer j = row_ptr[i]; j <= row_ptr[i + 1] - 1; ++j) {

			den += val[j] * val[j];

			if (col_ind[j] == i) {
				num += val[j]; // диагональ
			}
			else {
				sum += val[j] * x[col_ind[j]];
			}
		}
		doublereal M_precond = num / den;
		x[i] += M_precond * ((rthdsd[i] - sum) / num - x[i]);
	}
} // spai0_smoother

// 21.12.2019
// Sparse approximate inverse relaxation scheme.
// n Ц число неизвестных,
// x Ц вектор решени€ (текущее приближение к решению)
// val, col_ind, row_ptr матрица в CRS формате.
// rthdsd Ц права€ часть —Ћј”.
// ћетод практически не сходитс€, или очень очень долго сходитс€.
template <typename doublerealT>
void spai0_smoother(Ak2& Amat, integer istartq, integer iendq,
	doublerealT*& x, doublerealT*& b, integer*& row_ptr_start, integer*& row_ptr_end, integer iadd)
{

	std::cout << "incomming spai0_smoother no diagonal\n";
	system("pause");

	integer startpos = istartq + iadd;
	integer endpos = iendq + iadd;

	//#pragma omp parallel for 
	for (integer ii = startpos; ii <= endpos; ++ii) {
		integer istr = ii - iadd;

		// Ќачальное приближение всегда нулевое:
		x[istr] = 0.0;

		//doublereal num = 0.0;
		doublereal den = 0.0;
		doublereal sum = 0.0;
		doublereal E1 = 1.0 / Amat.aij[row_ptr_start[ii]];
		//num = E1;
		doublereal E2 = E1 * E1;// ƒл€ нормировки иначе происходит неправильна€ обработка.

		integer is1 = row_ptr_start[ii] + 1;
		integer is2 = row_ptr_end[ii];

		//num =  1.0 / Amat.aij[row_ptr_start[ii]];
		// 1.0 т.к. нормировка на диагональ.
		den += 1.0;// num* num;
		//den = num * num;

		// –аспараллеливание почемуто тормозит очень сильно.
		//#pragma omp parallel for reduction(+:rsum)
		for (integer ii1 = is1; ii1 <= is2; ++ii1)
		{

			den += (Amat.aij[ii1] * Amat.aij[ii1]) / E2;
			//den += (Amat.aij[ii1] * Amat.aij[ii1]);

			sum += Amat.aij[ii1] * x[Amat.j[ii1]];

		}
		doublereal M_precond = 1.0 / den; //num / den;
		x[istr] += M_precond * (Amat.aij[row_ptr_start[ii]] * (b[istr] - sum) - x[istr]);
		/*// debug
		if (x[istr] != x[istr]) {
			printf("%e 1\n", x[istr]);
			system("pause");
		}
		*/
		if (x[istr] != x[istr]) {
			std::cout << "no diag !!! M=" << M_precond << " E1=" << E1 << " den=" << den << "b= " << b[istr] << "diagminusone otsutstvuet" << std::endl;
			for (integer ii1 = is1 - 1; ii1 <= is2; ++ii1)
			{
				std::cout << "is1-1=" << is1 - 1 << " is2=" << is2 << std::endl;
				std::cout << "ii= " << ii;
				std::cout << " indx= " << row_ptr_start[ii] << std::endl;
				std::cout << " i=" << Amat.j[row_ptr_start[ii]];
				std::cout << " j=" << Amat.j[ii1] << " " << (ii1 == is1 - 1 ? E1 : Amat.aij[ii1]) << std::endl;
			}
			system("pause");
		}
	}
} // spai0_smoother

// 21.12.2019
// Sparse approximate inverse relaxation scheme.
// n Ц число неизвестных,
// x Ц вектор решени€ (текущее приближение к решению)
// val, col_ind, row_ptr матрица в CRS формате.
// rthdsd Ц права€ часть —Ћј”.
// ћетод практически не сходитс€, или очень очень долго сходитс€.
template <typename doublerealT>
void spai0_smoother(Ak2& Amat, integer istartq, integer iendq,
	doublerealT*& x, doublerealT*& b, integer*& row_ptr_start,
	integer*& row_ptr_end, integer iadd, doublerealT*& diag_minus_one)
{

	integer startpos = istartq + iadd;
	integer endpos = iendq + iadd;

	//#pragma omp parallel for 
	for (integer ii = startpos; ii <= endpos; ++ii) {
		integer istr = ii - iadd;

		// Ќачальное приближение всегда нулевое:
		x[istr] = 0.0;

		if (x[istr] != x[istr]) {
			std::cout << "apriory\n";
		}
		if (b[istr] != b[istr]) {
			std::cout << "apriory istr=" << istr << "\n";
		}

		//doublereal num = 0.0;
		doublereal den = 0.0;
		doublereal sum = 0.0;
		doublereal E1 = 1.0 / Amat.aij[row_ptr_start[ii]];
		//num = E1;// ƒиагональный элемент.
		doublereal E2 = Amat.aij[row_ptr_start[ii]] * Amat.aij[row_ptr_start[ii]];

		//E1 *= E1;// ƒл€ нормировки иначе происходит неправильна€ обработка.

		integer is1 = row_ptr_start[ii] + 1;
		integer is2 = row_ptr_end[ii];


		//num = 1.0;
		// 1.0 т.к. нормировка на диагональ. 
		den += 1.0;// num* num;
		//den += num * num;
		//den += fabs(num);

		//if (num != num) {
			//for (integer ii1 = is1 - 1; ii1 <= is2; ++ii1)
			//{
				//std::cout << "ii=" << ii<< Amat.j[ii1] << " " << Amat.aij[ii1] << std::endl;
			//}
			//system("pause");
		//}

		// –аспараллеливание почемуто тормозит очень сильно.
		//#pragma omp parallel for reduction(+:rsum)
		for (integer ii1 = is1; ii1 <= is2; ++ii1)
		{

			integer ipos = Amat.j[ii1];

			//den += (Amat.aij[ii1] * Amat.aij[ii1])/ E1;
			den += (Amat.aij[ii1] * Amat.aij[ii1]) * E2;
			//den += (Amat.aij[ii1] * Amat.aij[ii1]);
			//den += fabs(Amat.aij[ii1]);

			sum += Amat.aij[ii1] * x[ipos];

		}
		//if (is2<is1) den = 1.0;

		//doublereal M_precond = 1.0 / den;//num/den;
		//x[istr] = M_precond * (diag_minus_one[istr] * b[istr] - sum) / num;
		//x[istr] = M_precond * (diag_minus_one[istr] * b[istr] - sum)*Amat.aij[row_ptr_start[ii]];

		doublereal M_precond = 1.0 / den;
		//doublereal M_precond = E1 / den;
		//x[istr] += M_precond * (Amat.aij[row_ptr_start[ii]]*(b[istr] - sum)-x[istr]);
		x[istr] += M_precond * (Amat.aij[row_ptr_start[ii]] * (b[istr] - sum) - x[istr]);
		//x[istr] += M_precond * (diag_minus_one[istr] * b[istr] - sum);
		/* // debug
		if (x[istr] != x[istr]) {
			printf("%e M=%e num=%e den=%e b=%e sum=%e 2\n", x[istr], M_precond,num,den, b[istr],sum);
			system("pause");
		}
		*/
		if (x[istr] != x[istr]) {
			std::cout << "M=" << M_precond << " E1=" << E1 << " den=" << den << "b= " << b[istr] << "diagminusone=" << diag_minus_one[istr] << std::endl;
			for (integer ii1 = is1 - 1; ii1 <= is2; ++ii1)
			{
				std::cout << "is1-1=" << is1 - 1 << " is2=" << is2 << std::endl;
				std::cout << "ii= " << ii;
				std::cout << " indx= " << row_ptr_start[ii] << std::endl;
				std::cout << " i=" << Amat.j[row_ptr_start[ii]];
				std::cout << " j=" << Amat.j[ii1] << " " << (ii1 == is1 - 1 ? E1 : Amat.aij[ii1]) << std::endl;
			}
			system("pause");
		}
	}
} // spai0_smoother

#endif 