
#pragma once
#ifndef BASIC_FUNCTIONS_MY_AGREGAT_AMG_CPP
#define BASIC_FUNCTIONS_MY_AGREGAT_AMG_CPP 1

// 01.07.2021 CF порядок обхода теперь можно отключить или включить в меню пользователя.
// Вид нового условия во всех циклах.
// if ((!my_amg_manager.bcf_reorder) || (F_false_C_true[ii] == false)) { // F nodes
// (!my_amg_manager.bcf_reorder) служит для отключения режима CF избирательности.
// для thresholda theta теперь тожно есть возможэность использовать константу пользователя из меню на всех уровнях
// или поставить галочку для full automatic режима при котором theta по формуле вычисляется на каждом уровне индивидуально.

#include "my_sort_algorithm2.cpp"

// 15 11 2016 метод полностью на malloc с контролем через handle_error.
// Мы усилим РУМБА0.14 алгоритм ILU2 предобуславливанием на каждом уровне вложенности.
// Преобразует equation3D  формат хранения в CRS формат.
// Цель написания этого преобразователя: экономия оперативной памяти компьютера.
// Т.к. формат SIMPLESPARSE требует слишком много памяти.
void equation3DtoCRSRUMBA1(LEVEL_ADDITIONAL_DATA &milu2,
	bool ballocmemory, Ak2 &Amat, integer istartq, integer iendq,
	integer* &row_ptr_start, integer* &row_ptr_end, integer iadd, integer ilevel, 
	integer iVar) {

	// Если ballocmemory равен true то происходит выделение памяти.
	const bool bonly_negative_connections = false;
	bool flag = true;
	integer n = 0; // число ненулевых элементов
	integer maxelm_plus_maxbound = 0; // Число уравнений в СЛАУ.

	const doublereal nonzeroEPS = 1e-37; // для отделения вещественного нуля

	integer startpos = istartq + iadd;
	integer endpos = iendq + iadd;


	// подсчёт количества ненулевых элементов
	// во всех точках расчётной области как внутренних так и граничных.
	for (integer ii = startpos; ii <= endpos; ++ii) {
		maxelm_plus_maxbound++;
		integer istr = ii - iadd;
		if (fabs(Amat.aij[row_ptr_start[ii]]) < 1.0e-20)
		{
#if doubleintprecision == 1
			printf("zero diagonal coefficient in level = %lld, istr=%lld\n", ilevel, istr);
#else
			printf("zero diagonal coefficient in level = %d, istr=%d\n", ilevel, istr);
#endif

			//system("pause");
			system("PAUSE");
			exit(1);
		}
		doublereal ap_now = 1.0 / Amat.aij[row_ptr_start[ii]];

		//if (fabs(ap_now)> 1e10*nonzeroEPS) n++; // диагональный элемент
		if (ap_now > nonzeroEPS) n++; // Диагональный элемент.
		else {
			// 5 августа 2016. 
			flag = false;
			printf("internal zero diagonal element.\n");
			//printf("ae=%e aw=%e an=%e as=%e at=%e ab=%e sum_nb=%e", sl[k].ae, sl[k].aw, sl[k].an, sl[k].as, sl[k].at, sl[k].ab, sl[k].ae + sl[k].aw + sl[k].an + sl[k].as + sl[k].at + sl[k].ab);
			if (ap_now < 0.0) {
				printf("found negativ diagonal coefficient=%e...\n", ap_now);
			}
			printf("fatal error equation3DtoCRS in RUMBA...\n");
			//system("pause");
			system("PAUSE");
			exit(1);
			//n++;
			//sl[k].ap = fabs(sl[k].ae) + fabs(sl[k].aw) + fabs(sl[k].an) + fabs(sl[k].as) + fabs(sl[k].at) + fabs(sl[k].ab);
		}

		integer is1 = row_ptr_start[ii] + 1;
		integer is2 = row_ptr_end[ii];

		for (integer ii1 = is1; ii1 <= is2; ++ii1)
		{
			if (bonly_negative_connections) {
				if (Amat.aij[ii1] < -nonzeroEPS) {
					n++;
				}
			}
			else {
				if (fabs(Amat.aij[ii1]) > nonzeroEPS) n++;
			}
		}
	}


	if (flag) {
		// memory +15N
		// Теперь выделение памяти будет происходить централизованно, вне данного кода.
		// Это сделано для кода BICGSTAB_internal3. дата изменения 12 апреля 2013.
		// Другой код, использующий equation3dtoCRS может оказаться неработоспособным после этого изменения.
		if (ballocmemory) {
			// Важно выделить память с запасом, т.к. одна и таже память используется и для компонент скорости и для поправки давления.
			//val = new doublereal[7 * (maxelm + maxbound) + 2 * maxbound + 2];
			//col_ind = new integer[7 * (maxelm + maxbound) + 2 * maxbound + 2];

			//milu2.val = new doublereal[n + 2];
			milu2.val = (doublereal*)malloc((n + 2) * sizeof(doublereal));
			char c1[11] = "milu2.val";
			char c2[23] = "equation3DtoCRSRUMBA1";
			handle_error(milu2.val, c1, c2, (n + 2));
			//milu2.col_ind = new integer[n + 2];
			milu2.col_ind = (integer*)malloc((n + 2) * sizeof(integer));
			char c3[14] = "milu2.col_ind";
			handle_error(milu2.col_ind, c3, c2, (n + 2));
			//row_ptr = new integer[(maxelm + maxbound) + 1];
			//milu2.row_ptr = new integer[maxelm_plus_maxbound + 1];
			milu2.row_ptr = (integer*)malloc((maxelm_plus_maxbound + 1) * sizeof(integer));
			char c4[14] = "milu2.row_ptr";
			handle_error(milu2.row_ptr, c4, c2, (maxelm_plus_maxbound + 1));

			if ((milu2.val == nullptr) || (milu2.col_ind == nullptr) || (milu2.row_ptr == nullptr)) {
				// недостаточно памяти на данном оборудовании.
				printf("Problem: not enough memory on your equipment...\n");
				printf("Please any key to exit...\n");
				exit(1);
			}
		}


		// инициализация
#pragma omp parallel for
		for (integer k = 0; k<(n); ++k) {
			milu2.val[k] = 0.0;
			milu2.col_ind[k] = -1;
		}
#pragma omp parallel for
		for (integer k = 0; k <= (maxelm_plus_maxbound); ++k) {
			milu2.row_ptr[k] = n; // присваиваем количество ненулевых элементов плюс 1 с учётом того что нумерация массива начинается с 0
		}

		// Быстрая Сортировка Хоара.
		// упорядочивание по строкам
		//QuickSort(...); не требуется,
		// т.к. сама структура хранения 
		// подразумевает упорядочивание по строкам.

		/*
		// заполнение разреженной матрицы
		for (integer k=0; k<M.n; ++k) {
		val[k]=M.a[k].aij;
		col_ind[k]=M.a[k].j;
		row_ptr[M.a[k].i]=min(k,row_ptr[M.a[k].i]);
		}
		*/
		integer ik = 0; // счётчик ненулевых элементов СЛАУ

						// для всех узлов расчётной области как внутренних так и внешних:
		for (integer ii = startpos; ii <= endpos; ++ii) {
			integer istr = ii - iadd;
			integer k = istr - 1; // нумерация с нуля.
			doublereal ap_now = 1.0 / Amat.aij[row_ptr_start[ii]];


			// Диагональный коэффициент.
			if (fabs(ap_now) > nonzeroEPS) {
				milu2.val[ik] = ap_now;
				milu2.col_ind[ik] = istr - 1; // Нуменрация начинается с нуля.
				milu2.row_ptr[k] = (ik < milu2.row_ptr[k] ? ik : milu2.row_ptr[k]);
				ik++;
			}

			integer is1 = row_ptr_start[ii] + 1;
			integer is2 = row_ptr_end[ii];

			for (integer ii1 = is1; ii1 <= is2; ++ii1)
			{
				if (bonly_negative_connections) {
					if (Amat.aij[ii1] < -nonzeroEPS) {
						milu2.val[ik] = Amat.aij[ii1]; // уже с учётом знака.
						milu2.col_ind[ik] = Amat.j[ii1] - 1; // Нумерация начинается с нуля.
						milu2.row_ptr[k] = (ik < milu2.row_ptr[k] ? ik : milu2.row_ptr[k]);
						ik++;
					}
				}
				else {
					if (fabs(Amat.aij[ii1]) > nonzeroEPS) {
						milu2.val[ik] = Amat.aij[ii1]; // уже с учётом знака.
						milu2.col_ind[ik] = Amat.j[ii1] - 1; // Нумерация начинается с нуля.
						milu2.row_ptr[k] = (ik < milu2.row_ptr[k] ? ik : milu2.row_ptr[k]);
						ik++;
					}
				}
			}
		}

		// в каждой строке элементы отсортированы по номерам столбцов:
#pragma omp parallel for
		for (integer k = 0; k < (maxelm_plus_maxbound); ++k) {
			QuickSortCSIR(milu2.col_ind, milu2.val, milu2.row_ptr[k] + 1, milu2.row_ptr[k + 1] - 1);
		}

#if doubleintprecision == 1
		//printf("n==%lld row_ptr=%lld\n", n, milu2.row_ptr[maxelm_plus_maxbound]);
#else
		//printf("n==%d row_ptr=%d\n", n, milu2.row_ptr[maxelm_plus_maxbound]);
#endif

		//system("pause");

		milu2.maxelm_plus_maxbound = maxelm_plus_maxbound;
#pragma omp parallel for
		for (integer i = 0; i<milu2.row_ptr[maxelm_plus_maxbound]; ++i) {
			milu2.col_ind[i] = milu2.col_ind[i] + 1;
		}
#pragma omp parallel for
		for (integer i = 0; i<maxelm_plus_maxbound + 1; ++i) {
			milu2.row_ptr[i] = milu2.row_ptr[i] + 1;
		}

		/*
		FILE *fp;
		errno_t err;
		// создание файла для записи.
		if ((err = fopen_s( &fp, "matr.txt", "w")) != 0) {
		printf("Create File Error\n");
		}
		else {

		// debug
		for (k=0; k<=maxelm+maxbound; ++k) {
		#if doubleintprecision == 1
		fprintf(fp,"%lld ",milu2.row_ptr[k]);
		#else
		fprintf(fp,"%d ",milu2.row_ptr[k]);
		#endif

		}
		fprintf(fp,"\n");
		for (k=0; k<milu2.row_ptr[maxelm+maxbound]; ++k) {
		#if doubleintprecision == 1
		fprintf(fp, "%e %lld\n",milu2.val[k],milu2.col_ind[k]);
		#else
		fprintf(fp, "%e %d\n",milu2.val[k],milu2.col_ind[k]);
		#endif

		}

		fclose(fp);
		}
		printf("ready");
		system("pause");
		*/
		//doublereal radd = 0.0;
		/*
		if (ilevel == 0) {
		radd = 1.8;
		}
		else if (ilevel == 1) {
		radd = 4.2;
		}
		else if ((ilevel >=2)&&(ilevel<=5)) {
		radd = 7.0;
		}
		else if (ilevel==6) {
		radd = 7.0;
		}
		else {
		radd = 8.0;
		}
		*/
		if (1) {

			// Параметры подбирал на АЛИС сетке.
			integer ilevel_loc = ilevel;
			if (ilevel == 0) ilevel_loc = 1;
			//radd = 3.3*(n / (5.0*maxelm_plus_maxbound));
			//milu2.iwk = static_cast<integer>((milu2.lfil + 1 + radd) * n + 4 * maxelm_plus_maxbound);
			//milu2.iwk = n + 4; // ilu0 не требует лишней памяти.
			integer lfil = my_amg_manager.lfil;

			// Усиливаем на первом самом разреженном уровне так как имеют место быть проблемы со сходимостью.
			if ((ilevel == 0) && (iVar == TOTALDEFORMATIONVAR)) lfil = 4;
			if (lfil == 0) {
				milu2.iwk = static_cast<integer>((lfil + 1) * (n + 2 * maxelm_plus_maxbound + 2) + 4 * maxelm_plus_maxbound);
			}
			else if (lfil==1) {
				// 2.4 +; 2.2 -; 2.3 + opt; 
				doublereal m_1 = 4.3;
				if (ilevel_loc == 10) m_1 = 4.0;
				if (ilevel_loc >= 11) m_1 = 5.0;
				if (ilevel_loc >= 15) m_1 = 8.0;
				milu2.iwk = static_cast<integer>((m_1 * ilevel_loc * lfil + 1) * (n + 2 * maxelm_plus_maxbound + 2) + 4 * maxelm_plus_maxbound);
				//std::cout << "m_1=" <<  m_1 << std::endl;
			}
			else if (lfil == 2) {
				// milu2.iwk = (40 * lfil + 1) * (n + 2 * maxelm_plus_maxbound + 2) + 4 * maxelm_plus_maxbound;
				// ilevel
				// lfil==1 3 is work Ok
				// milu2.iwk = (3* ilevel * lfil + 1) * (n + 2 * maxelm_plus_maxbound + 2) + 4 * maxelm_plus_maxbound;
				// 2.5 -; 2.6 + opt;
				doublereal m_1 = 2.6;
				if (ilevel_loc == 10) m_1 = 4.3;
				if (ilevel_loc >= 11) m_1 = 5.3;
				if (ilevel_loc >= 15) m_1 = 8.0;
				milu2.iwk = static_cast<integer>((m_1 * ilevel_loc * lfil + 1) * (n + 2 * maxelm_plus_maxbound + 2) + 4 * maxelm_plus_maxbound);
			}
			else if (lfil == 3) {
				// 3 -; 3.2 + opt;
				milu2.iwk = static_cast<integer>((3.2 * ilevel_loc*lfil + 1) * (n + 2 * maxelm_plus_maxbound + 2) + 4 * maxelm_plus_maxbound);
			}
			else if (lfil == 4) {
				// 3.8 -; 3.9 + opt;
				milu2.iwk = static_cast<integer>((3.9 * ilevel_loc * lfil + 1) * (n + 2 * maxelm_plus_maxbound + 2) + 4 * maxelm_plus_maxbound);
			}
			else if (lfil == 5) {
				// 4.4 -; 
				milu2.iwk = static_cast<integer>((4.6 * ilevel_loc * lfil + 1) * (n + 2 * maxelm_plus_maxbound + 2) + 4 * maxelm_plus_maxbound);
			}
			else if (lfil == 6) {
				milu2.iwk = static_cast<integer>((5.2 * ilevel_loc * lfil + 1) * (n + 2 * maxelm_plus_maxbound + 2) + 4 * maxelm_plus_maxbound);
			}
			else if (lfil >= 7) {
				milu2.iwk = static_cast<integer>((5.8 * ilevel_loc * lfil + 1) * (n + 2 * maxelm_plus_maxbound + 2) + 4 * maxelm_plus_maxbound);
			}
		}

		// размерность памяти под матрицу предобуславливания.
		// Выделение оперативной памяти может быть 
		// повторным, поэтому необходимо освободить оперативную память,
		// перед повторным выделением.
		if (milu2.alu != nullptr) {
			//delete[] milu2.alu;
			free(milu2.alu);
			milu2.alu = nullptr;
		}
		//milu2.alu = new doublereal[milu2.iwk + 2]; // +2 запас по памяти.
		milu2.alu = (doublereal*)malloc((milu2.iwk + 2) * sizeof(doublereal));
		char c1[23] = "equation3DtoCRSRUMBA1";
		char c2[11] = "milu2.alu";
		handle_error(milu2.alu, c2, c1, (milu2.iwk + 2));
		for (integer i_63 = 0; i_63 < (milu2.iwk + 2); ++i_63) {
			milu2.alu[i_63] = 0.0; // Инициализация матрицы.
		}
		// Выделение оперативной памяти может быть 
		// повторным, поэтому необходимо освободить оперативную память,
		// перед повторным выделением.
		if (milu2.jlu != nullptr) {
			//delete[] milu2.jlu;
			free(milu2.jlu);
			milu2.jlu = nullptr;
		}
		//milu2.jlu = new integer[milu2.iwk + 2];
		milu2.jlu = (integer*)malloc((milu2.iwk + 2) * sizeof(integer));
		char c3[11] = "milu2.jlu";
		handle_error(milu2.jlu, c3, c1, (milu2.iwk + 2));
		for (integer i_63 = 0; i_63 < (milu2.iwk + 2); ++i_63) {
			milu2.jlu[i_63] = -1; // Инициализация матрицы. Не существующее значение столбца.
		}
		// Выделение оперативной памяти может быть 
		// повторным, поэтому необходимо освободить оперативную память,
		// перед повторным выделением.
		if (milu2.ju != nullptr) {
			//delete[] milu2.ju;
			free(milu2.ju);
			milu2.ju = nullptr;
		}
		//milu2.ju = new integer[maxelm_plus_maxbound + 2];
		milu2.ju = (integer*)malloc((maxelm_plus_maxbound + 2) * sizeof(integer));
		char c4[10] = "milu2.ju";
		handle_error(milu2.ju, c4, c1, (maxelm_plus_maxbound + 2));
		for (integer i_63 = 0; i_63 < (maxelm_plus_maxbound + 2); ++i_63) {
			milu2.ju[i_63] = -1; // Инициализация матрицы. Не существующее значение.
		}
		//milu2.levs = new integer[milu2.iwk + 2]; // уровень.
		milu2.levs = (integer*)malloc((milu2.iwk + 2) * sizeof(integer));
		char c5[12] = "milu2.levs";
		handle_error(milu2.levs, c5, c1, (milu2.iwk + 2));
		//milu2.w = new doublereal[maxelm_plus_maxbound + 2]; // +2 запас по памяти.
		milu2.w = (doublereal*)malloc((maxelm_plus_maxbound + 2) * sizeof(doublereal));
		char c6[9] = "milu2.w";
		handle_error(milu2.w, c6, c1, (maxelm_plus_maxbound + 2));
		char c7[10] = "milu2.jw";
		if (1) {
			//milu2.jw = new integer[3 * maxelm_plus_maxbound + 2]; // +2 запас по памяти.
			milu2.jw = (integer*)malloc((4 * maxelm_plus_maxbound + 2) * sizeof(integer));
			handle_error(milu2.jw, c7, c1, (4 * maxelm_plus_maxbound + 2));
		}
		else {
			if (ilevel < 6) {
				///milu2.jw = new integer[5 * maxelm_plus_maxbound + 2]; // +2 запас по памяти.
				milu2.jw = (integer*)malloc((5 * maxelm_plus_maxbound + 2) * sizeof(integer));
				handle_error(milu2.jw, c7, c1, (5 * maxelm_plus_maxbound + 2));
			}
			else {
				//milu2.jw = new integer[12 * maxelm_plus_maxbound + 2]; // +2 запас по памяти.
				milu2.jw = (integer*)malloc((12 * maxelm_plus_maxbound + 2) * sizeof(integer));
				handle_error(milu2.jw, c7, c1, (12 * maxelm_plus_maxbound + 2));
			}
		}
		if ((milu2.alu == nullptr) || (milu2.jlu == nullptr) || (milu2.levs == nullptr) || (milu2.ju == nullptr) || (milu2.w == nullptr) || (milu2.jw == nullptr)) {
			// недостаточно памяти на данном оборудовании.
			printf("Problem: not enough memory on your equipment...\n");
			printf("Please any key to exit...\n");
			exit(1);
		}

		// копии объектов
		// 9 ноября 2016 Копии объектов теперь хранятся централизованно:
		// см. LEVEL_ADDITIONAL_DATA_BUFER milu_gl_buffer; in ilut.c module.
		//milu2.alu_copy = new doublereal[milu2.iwk + 2]; // +2 запас по памяти.
		//milu2.jlu_copy = new integer[milu2.iwk + 2];
		//milu2.ju_copy = new integer[maxelm_plus_maxbound + 2];
		//milu2.b_copy = new doublereal[maxelm_plus_maxbound + 2];
		char c8[13] = "milu2.b_copy";
		milu2.b_copy = (doublereal*)malloc((maxelm_plus_maxbound + 2) * sizeof(doublereal));
		handle_error(milu2.b_copy, c8, c1, (maxelm_plus_maxbound + 2));
		//milu2.x_copy = new doublereal[maxelm_plus_maxbound + 2];
		milu2.x_copy = (doublereal*)malloc((maxelm_plus_maxbound + 2) * sizeof(doublereal));
		char c9[13] = "milu2.x_copy";
		handle_error(milu2.x_copy, c9, c1, (maxelm_plus_maxbound + 2));
		//milu2.zbuf = new doublereal[maxelm_plus_maxbound + 2];
		milu2.zbuf = (doublereal*)malloc((maxelm_plus_maxbound + 2) * sizeof(doublereal));
		char c10[11] = "milu2.zbuf";
		handle_error(milu2.zbuf, c10, c1, (maxelm_plus_maxbound + 2));
		//milu2.zbuf2 = new doublereal[maxelm_plus_maxbound + 2];
		milu2.zbuf2 = (doublereal*)malloc((maxelm_plus_maxbound + 2) * sizeof(doublereal));
		char c11[12] = "milu2.zbuf2";
		handle_error(milu2.zbuf2, c11, c1, (maxelm_plus_maxbound + 2));
		for (integer i_63 = 0; i_63 < (maxelm_plus_maxbound + 2); ++i_63) {
			milu2.zbuf2[i_63] = 0.0;
		}

		//if ((milu2.alu_copy == nullptr) || (milu2.jlu_copy == nullptr) || (milu2.ju_copy == nullptr) || (milu2.b_copy == nullptr) || (milu2.x_copy == nullptr) || (milu2.zbuf==nullptr)||(milu2.zbuf2==nullptr)) {
		if ((milu2.b_copy == nullptr) || (milu2.x_copy == nullptr) || (milu2.zbuf == nullptr) || (milu2.zbuf2 == nullptr)) {
			// недостаточно памяти на данном оборудовании.
			printf("Problem: not enough memory on your equipment...\n");
			printf("Please any key to exit...\n");
			exit(1);
		}

		//milu2.lfil = 2;
		integer ierr = 0;
		iluk_(maxelm_plus_maxbound, milu2.val, milu2.col_ind, milu2.row_ptr, my_amg_manager.lfil/*milu2.lfil*/, milu2.alu, milu2.jlu, milu2.ju, milu2.levs, milu2.iwk, milu2.w, milu2.jw, ierr);
		if (ierr != 0) {
#if doubleintprecision == 1
			printf("ierr=%lld\n", ierr);
#else
			printf("ierr=%d\n", ierr);
#endif

			//system("pause");
			system("PAUSE");
			exit(1);
		}

		if ((ierr == -2) || (ierr == -3)) {

			integer ipassage = 1;
			do {
				printf("\nPlease WAIT... ... ...\n");

				// задаче не хватило памяти, значит нужно перевыделить !
				//if (milu2.alu != nullptr) delete[] milu2.alu;
				//if (milu2.jlu != nullptr) delete[] milu2.jlu;
				//if (milu2.levs != nullptr) delete[] milu2.levs;
				if (milu2.alu != nullptr) free(milu2.alu);
				if (milu2.jlu != nullptr) free(milu2.jlu);
				if (milu2.levs != nullptr) free(milu2.levs);
				//if (milu2.alu_copy != nullptr) delete[] milu2.alu_copy;
				//if (milu2.jlu_copy != nullptr) delete[] milu2.jlu_copy;

				// инициализация !
				milu2.alu = nullptr;
				milu2.jlu = nullptr;
				milu2.levs = nullptr;
				//milu2.alu_copy = nullptr;
				//milu2.jlu_copy = nullptr;
				integer ilevel_loc = ilevel;
				if (ilevel == 0) ilevel_loc = 1;
				integer lfil = my_amg_manager.lfil;
				if (lfil == 0) {
					milu2.iwk = static_cast<integer>((lfil + 1) * (n + 2 * maxelm_plus_maxbound + 2) + (4 + 3 * ipassage) * maxelm_plus_maxbound);
				}
				else if (lfil == 1) {
					milu2.iwk = static_cast<integer>((2.3* ilevel_loc * lfil + 1) * (n + 2 * maxelm_plus_maxbound + 2) + (4 + 3 * ipassage) * maxelm_plus_maxbound);
				}
				else if (lfil == 2) {
					milu2.iwk = static_cast<integer>((2.6* ilevel_loc * lfil + 1) * (n + 2 * maxelm_plus_maxbound + 2) + (4 + 3 * ipassage) * maxelm_plus_maxbound);
				}
				else if (lfil == 3) {
					milu2.iwk = static_cast<integer>((3.2* ilevel_loc * lfil + 1) * (n + 2 * maxelm_plus_maxbound + 2) + (4 + 3 * ipassage) * maxelm_plus_maxbound);
				}
				else if (lfil == 4) {
					milu2.iwk = static_cast<integer>((3.9 * ilevel_loc * lfil + 1) * (n + 2 * maxelm_plus_maxbound + 2) + (4 + 3 * ipassage) * maxelm_plus_maxbound);
				}
				else if (lfil == 5) {
					milu2.iwk = static_cast<integer>((4.4 * ilevel_loc * lfil + 1) * (n + 2 * maxelm_plus_maxbound + 2) + (4 + 3 * ipassage) * maxelm_plus_maxbound);
				}
				else if (lfil == 6) {
					milu2.iwk = static_cast<integer>((5.0 * ilevel_loc * lfil + 1) * (n + 2 * maxelm_plus_maxbound + 2) + (4 + 3 * ipassage) * maxelm_plus_maxbound);
				}
				else if (lfil >= 7) {
					milu2.iwk = static_cast<integer>((5.6 * ilevel_loc * lfil + 1) * (n + 2 * maxelm_plus_maxbound + 2) + (4 + 3 * ipassage) * maxelm_plus_maxbound);
				}
				
				// Выделение оперативной памяти может быть 
		        // повторным, поэтому необходимо освободить оперативную память,
		        // перед повторным выделением.
				if (milu2.alu != nullptr) {
					//delete[] milu2.alu;
					free(milu2.alu);
					milu2.alu = nullptr;
				}
				//milu2.alu = new doublereal[milu2.iwk + 2]; // +2 запас по памяти.
				milu2.alu = (doublereal*)malloc((milu2.iwk + 2) * sizeof(doublereal));
				char c12[11] = "milu2.alu";
				handle_error(milu2.alu, c12, c1, (milu2.iwk + 2));
				// Выделение оперативной памяти может быть 
		        // повторным, поэтому необходимо освободить оперативную память,
	        	// перед повторным выделением.
				if (milu2.jlu != nullptr) {
					//delete[] milu2.jlu;
					free(milu2.jlu);
					milu2.jlu = nullptr;
				}
				//milu2.jlu = new integer[milu2.iwk + 2];
				milu2.jlu = (integer*)malloc((milu2.iwk + 2) * sizeof(integer));
				char c13[11] = "milu2.jlu";
				handle_error(milu2.jlu, c13, c1, (milu2.iwk + 2));
				//milu2.levs = new integer[milu2.iwk + 2]; // уровень.
				milu2.levs = (integer*)malloc((milu2.iwk + 2) * sizeof(integer));
				char c14[11] = "milu2.levs";
				handle_error(milu2.levs, c14, c1, (milu2.iwk + 2));
				//milu2.alu_copy = new doublereal[milu2.iwk + 2]; // +2 запас по памяти.
				//milu2.jlu_copy = new integer[milu2.iwk + 2];

				//if ((milu2.alu_copy != nullptr) && (milu2.jlu_copy != nullptr) && (milu2.alu != nullptr) && (milu2.jlu != nullptr) && (milu2.levs != nullptr)) {
				if ((milu2.alu != nullptr) && (milu2.jlu != nullptr) && (milu2.levs != nullptr)) {
					iluk_(maxelm_plus_maxbound, milu2.val, milu2.col_ind, milu2.row_ptr, lfil, milu2.alu, milu2.jlu, milu2.ju, milu2.levs, milu2.iwk, milu2.w, milu2.jw, ierr);
				}
				else {
					// недостаточно памяти на данном оборудовании.
					ipassage = 4;
					printf("Problem: not enough memory on your equipment...\n");
					printf("Please any key to exit...\n");
					exit(1);
				}

				ipassage++;
				// Пока просто контролируем, в добавок хорошо бы выводить информацию на каком уровне это произошло.
#if doubleintprecision == 1
				printf("control memory allocation for ilu2 in RUMBA0.14 in level=%lld\n", ilevel);
#else
				printf("control memory allocation for ilu2 in RUMBA0.14 in level=%d\n", ilevel);
#endif

				system("PAUSE");
			} while ((ierr != 0) && (ipassage<4));

			if (ipassage == 4) {
				printf("Error memory alloc !!!\n");
				printf("failed to obtain an expansion for the 4 approaches...\n");
				printf("Please any key to exit...\n");
				//system("pause");
				system("PAUSE");
				exit(1);
			}
		}

		for (integer k = 0; k<n; ++k) if (milu2.col_ind[k] == (-1)) {
			printf("Error equation3D to CRS.\n");
			//system("pause");
			system("pause");
		}

		// Нумерация начинается с единицы (это важно).
		for (integer k = 0; k<maxelm_plus_maxbound; ++k) {
			if (milu2.val[milu2.row_ptr[k] - 1]<nonzeroEPS) {
#if doubleintprecision == 1
				printf("negativ diagonal element equation3DtoCRS %lld\n", k);
#else
				printf("negativ diagonal element equation3DtoCRS %d\n", k);
#endif

				//system("pause");
				system("pause");
			}
		}

		if (1 || milu2.lfil == 0) {
			// экономим оперативную память.
			//if (milu2.val != nullptr) delete[] milu2.val;
			//if (milu2.col_ind != nullptr) delete[] milu2.col_ind;
			//if (milu2.row_ptr != nullptr) delete[] milu2.row_ptr;
			if (milu2.val != nullptr) free(milu2.val);
			if (milu2.col_ind != nullptr) free(milu2.col_ind);
			if (milu2.row_ptr != nullptr) free(milu2.row_ptr);
			milu2.val = nullptr;
			milu2.col_ind = nullptr;
			milu2.row_ptr = nullptr;
		}

	}

	if (!flag) {
		printf("Error equation 3D to CRS: zero diagonal element...\n");
		//system("pause");
		system("pause");
	}


} // equation3DtoCRSRUMBA1


// 31.12.2019 iluk сглаживатель для amg1r5 алгоритма. Окончание 1.01.2020.
// 15.11.2016 метод полностью на malloc с контролем через handle_error.
// Мы усилим РУМБА0.14 алгоритм ILU2 предобуславливанием на каждом уровне вложенности.
// Преобразует equation3D  формат хранения в CRS формат.
// Цель написания этого преобразователя: экономия оперативной памяти компьютера.
// Т.к. формат SIMPLESPARSE требует слишком много памяти.
void equation3DtoCRSRUMBA1_amg1r5(LEVEL_ADDITIONAL_DATA& milu2,
	doublereal* a, integer* ia, integer* ja,
	integer* imin, integer* imax, integer ilevel) {

	// Если ballocmemory равен true то происходит выделение памяти.
	//const bool bonly_negative_connections = false;
	
	integer n = 0; // число ненулевых элементов
	integer maxelm_plus_maxbound = 0; // Число уравнений в СЛАУ.

	const doublereal nonzeroEPS = 1e-37; // для отделения вещественного нуля

	//integer istartq = imin[ilevel + 1];
	//integer iendq = imax[ilevel + 1];
	maxelm_plus_maxbound = 1 + imax[ilevel + 1] - imin[ilevel + 1];
	n = ia[imax[ilevel + 1] + 1] - 1 - ia[imin[ilevel + 1]] + 1;
	

	//integer startpos = istartq;
	//integer endpos = iendq;

	bool ballocmemory = false;
	//if (milu2.jw == nullptr) ballocmemory = true;
	if (bflag_visit_amg1r5[ilevel] == false) {
		ballocmemory = true;
		bflag_visit_amg1r5[ilevel] = true;
	}
	
	// memory +15N
	// Теперь выделение памяти будет происходить централизованно, вне данного кода.
	// Это сделано для кода BICGSTAB_internal3. дата изменения 12 апреля 2013.
	// Другой код, использующий equation3dtoCRS может оказаться неработоспособным после этого изменения.
	if (ballocmemory) {

		// Замер времени.
		unsigned int calculation_main_start_time = 0; // начало счёта мс.
		unsigned int calculation_main_end_time = 0; // окончание счёта мс.


		calculation_main_start_time = clock(); // момент начала счёта.

		if (!(steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::NETWORK_T_UNSTEADY)) {
			printf("apply ilu%lld smoother for number ", my_amg_manager.lfil);
			printf("%2lld level %3d ", ilevel, (int)(n / (maxelm_plus_maxbound)));
		}

		// Важно выделить память с запасом, т.к. одна и таже память используется 
		// и для компонент скорости и для поправки давления.
		//val = new doublereal[7 * (maxelm + maxbound) + 2 * maxbound + 2];
		//col_ind = new integer[7 * (maxelm + maxbound) + 2 * maxbound + 2];
		
		//milu2.val = new doublereal[n + 2];
		milu2.val = (doublereal*)malloc((n + 2) * sizeof(doublereal));
		char c1[11] = "milu2.val";
		char c2[23] = "equation3DtoCRSRUMBA1";
		handle_error(milu2.val, c1, c2, (n + 2));
		//milu2.col_ind = new integer[n + 2];
		milu2.col_ind = (integer*)malloc((n + 2) * sizeof(integer));
		char c3[14] = "milu2.col_ind";
		handle_error(milu2.col_ind, c3, c2, (n + 2));
		//row_ptr = new integer[(maxelm + maxbound) + 1];
		//milu2.row_ptr = new integer[maxelm_plus_maxbound + 1];
		milu2.row_ptr = (integer*)malloc((maxelm_plus_maxbound + 1) * sizeof(integer));
		char c4[14] = "milu2.row_ptr";
		handle_error(milu2.row_ptr, c4, c2, (maxelm_plus_maxbound + 1));

		if ((milu2.val == nullptr) || (milu2.col_ind == nullptr) || (milu2.row_ptr == nullptr)) {
			// недостаточно памяти на данном оборудовании.
			printf("Problem: not enough memory on your equipment...\n");
			printf("Please any key to exit...\n");
			exit(1);
		}

		integer iminus = 0;
		if (ilevel > 0) iminus = ja[ia[imin[ilevel + 1]]]-1;

		integer k = 0, k_1 = 0;
		integer i__1 = imax[ilevel + 1];
		for (integer i__ = imin[ilevel + 1]; i__ <= i__1; ++i__) {
			milu2.row_ptr[k_1] = k;
			k_1++;

			integer i__2 = ia[i__ + 1] - 1;
			for (integer j = ia[i__]; j <= i__2; ++j) {
				milu2.val[k] = a[j];
				milu2.col_ind[k] = ja[j] -1 -iminus;
				k++;
			}
		}
		milu2.row_ptr[k_1] = k; // n==k гарантированно.
	

		// в каждой строке элементы отсортированы по номерам столбцов:
#pragma omp parallel for
		for (integer k = 0; k < (maxelm_plus_maxbound); ++k) {
			QuickSortCSIR(milu2.col_ind, milu2.val, milu2.row_ptr[k] + 1, milu2.row_ptr[k + 1] - 1);
		}

		milu2.maxelm_plus_maxbound = maxelm_plus_maxbound;
		
#pragma omp parallel for
		for (integer i = 0; i < milu2.row_ptr[maxelm_plus_maxbound]; ++i) {
			milu2.col_ind[i] = milu2.col_ind[i] + 1;
		}
#pragma omp parallel for
		for (integer i = 0; i < maxelm_plus_maxbound + 1; ++i) {
			milu2.row_ptr[i] = milu2.row_ptr[i] + 1;
		}
		
			
		//doublereal radd = 0.0;
			
		if (1) {

			// Параметры подбирал на АЛИС сетке.
			integer ilevel_loc = ilevel;
			if (ilevel == 0) ilevel_loc = 1;
		   	//radd = 3.3*(n / (5.0*maxelm_plus_maxbound));
		    //milu2.iwk = static_cast<integer>((milu2.lfil + 1 + radd) * n + 4 * maxelm_plus_maxbound);
			//milu2.iwk = n + 4; // ilu0 не требует лишней памяти.
			integer lfil = my_amg_manager.lfil;
			if (lfil == 0) {
				milu2.iwk = static_cast<integer>((lfil + 1) * (n + 2 * maxelm_plus_maxbound + 2) + 4 * maxelm_plus_maxbound);
			}
			else if (lfil == 1) {
				// 2.4 +; 2.2 -; 2.3 + opt; 
				doublereal m_1 = 4.3;
				if (ilevel_loc == 10) m_1 = 4.0;
				if (ilevel_loc >= 11) m_1 = 5.0;
				if (ilevel_loc >= 15) m_1 = 8.0;
				milu2.iwk = static_cast<integer>((m_1 * ilevel_loc * lfil + 1) * (n + 2 * maxelm_plus_maxbound + 2) + 4 * maxelm_plus_maxbound);
				//std::cout << "m_1=" << m_1 << std::endl; // Механика на основе amg1r5.
			}
			else if (lfil == 2) {
				// milu2.iwk = (40 * lfil + 1) * (n + 2 * maxelm_plus_maxbound + 2) + 4 * maxelm_plus_maxbound;
				// ilevel
				// lfil==1 3 is work Ok
				// milu2.iwk = (3* ilevel * lfil + 1) * (n + 2 * maxelm_plus_maxbound + 2) + 4 * maxelm_plus_maxbound;
				// 2.5 -; 2.6 + opt;
				doublereal m_1 = 2.6;
				if (ilevel_loc == 10) m_1 = 4.3;
				if (ilevel_loc >= 11) m_1 = 5.3;
				if (ilevel_loc >= 15) m_1 = 8.0;
				milu2.iwk = static_cast<integer>((m_1 * ilevel_loc * lfil + 1) * (n + 2 * maxelm_plus_maxbound + 2) + 4 * maxelm_plus_maxbound);
			}
			else if (lfil == 3) {
				// 3 -; 3.2 + opt;
				milu2.iwk = static_cast<integer>((3.2 * ilevel_loc * lfil + 1) * (n + 2 * maxelm_plus_maxbound + 2) + 4 * maxelm_plus_maxbound);
			}
			else if (lfil == 4) {
				// 3.8 -; 3.9 + opt;
				milu2.iwk = static_cast<integer>((3.9 * ilevel_loc * lfil + 1) * (n + 2 * maxelm_plus_maxbound + 2) + 4 * maxelm_plus_maxbound);
			}
			else if (lfil == 5) {
				// 4.4 -; 
				milu2.iwk = static_cast<integer>((4.6 * ilevel_loc * lfil + 1) * (n + 2 * maxelm_plus_maxbound + 2) + 4 * maxelm_plus_maxbound);
			}
			else if (lfil == 6) {
				milu2.iwk = static_cast<integer>((5.2 * ilevel_loc * lfil + 1) * (n + 2 * maxelm_plus_maxbound + 2) + 4 * maxelm_plus_maxbound);
			}
			else if (lfil >= 7) {
				milu2.iwk = static_cast<integer>((5.8 * ilevel_loc * lfil + 1) * (n + 2 * maxelm_plus_maxbound + 2) + 4 * maxelm_plus_maxbound);
			}
		}

		// Выделение оперативной памяти может быть 
		// повторным, поэтому необходимо освободить оперативную память,
		// перед повторным выделением.
		if (milu2.alu != nullptr) {
			//delete[] milu2.alu;
			free(milu2.alu);
			milu2.alu = nullptr;
		}
		// размерность памяти под матрицу предобуславливания.
		//milu2.alu = new doublereal[milu2.iwk + 2]; // +2 запас по памяти.
		milu2.alu = (doublereal*)malloc((milu2.iwk + 2) * sizeof(doublereal));
		char c11[23] = "equation3DtoCRSRUMBA1";
		char c21[11] = "milu2.alu";
		handle_error(milu2.alu, c21, c11, (milu2.iwk + 2));
		// Инициализация матрицы нулевыми значениями.
		for (integer i_63 = 0; i_63 < (milu2.iwk + 2); ++i_63) {
			milu2.alu[i_63] = 0.0;
		}
		// Выделение оперативной памяти может быть 
		// повторным, поэтому необходимо освободить оперативную память,
		// перед повторным выделением.
		if (milu2.jlu != nullptr) {
			//delete[] milu2.jlu;
			free(milu2.jlu);
			milu2.jlu = nullptr;
		}
		//milu2.jlu = new integer[milu2.iwk + 2];
		milu2.jlu = (integer*)malloc((milu2.iwk + 2) * sizeof(integer));
		char c31[11] = "milu2.jlu";
		handle_error(milu2.jlu, c31, c11, (milu2.iwk + 2));
		// Инициализация матрицы несуществующими номерами столбцов.
		for (integer i_63 = 0; i_63 < (milu2.iwk + 2); ++i_63) {
			milu2.jlu[i_63] = -1; // Столбец с номером -1 не существует в матрице.
		}
		// Выделение оперативной памяти может быть 
		// повторным, поэтому необходимо освободить оперативную память,
		// перед повторным выделением.
		if (milu2.ju != nullptr) {
			//delete[] milu2.ju;
			free(milu2.ju);
			milu2.ju = nullptr;
		}
		//milu2.ju = new integer[maxelm_plus_maxbound + 2];
		milu2.ju = (integer*)malloc((maxelm_plus_maxbound + 2) * sizeof(integer));
		char c41[10] = "milu2.ju";
		handle_error(milu2.ju, c41, c11, (maxelm_plus_maxbound + 2));
		// Инициализация матрицы несуществующими номерами столбцов.
		for (integer i_63 = 0; i_63 < (maxelm_plus_maxbound + 2); ++i_63) {
			milu2.ju[i_63] = -1; // Значение -1 не существует в матрице.
		}
		//milu2.levs = new integer[milu2.iwk + 2]; // уровень.
		milu2.levs = (integer*)malloc((milu2.iwk + 2) * sizeof(integer));
		char c5[12] = "milu2.levs";
		handle_error(milu2.levs, c5, c1, (milu2.iwk + 2));
		//milu2.w = new doublereal[maxelm_plus_maxbound + 2]; // +2 запас по памяти.
		milu2.w = (doublereal*)malloc((maxelm_plus_maxbound + 2) * sizeof(doublereal));
		char c6[9] = "milu2.w";
		handle_error(milu2.w, c6, c1, (maxelm_plus_maxbound + 2));
		char c7[10] = "milu2.jw";
		if (1) {
			//milu2.jw = new integer[3 * maxelm_plus_maxbound + 2]; // +2 запас по памяти.
			milu2.jw = (integer*)malloc((4 * maxelm_plus_maxbound + 2) * sizeof(integer));
			handle_error(milu2.jw, c7, c1, (4 * maxelm_plus_maxbound + 2));
		}
		else {
			if (ilevel < 6) {
				///milu2.jw = new integer[5 * maxelm_plus_maxbound + 2]; // +2 запас по памяти.
				milu2.jw = (integer*)malloc((5 * maxelm_plus_maxbound + 2) * sizeof(integer));
				handle_error(milu2.jw, c7, c1, (5 * maxelm_plus_maxbound + 2));
			}
			else {
				//milu2.jw = new integer[12 * maxelm_plus_maxbound + 2]; // +2 запас по памяти.
				milu2.jw = (integer*)malloc((12 * maxelm_plus_maxbound + 2) * sizeof(integer));
				handle_error(milu2.jw, c7, c1, (12 * maxelm_plus_maxbound + 2));
			}
		}
		if ((milu2.alu == nullptr) || (milu2.jlu == nullptr) ||
			(milu2.levs == nullptr) || (milu2.ju == nullptr) ||
			(milu2.w == nullptr) || (milu2.jw == nullptr)) {
			// недостаточно памяти на данном оборудовании.
			printf("Problem: not enough memory on your equipment...\n");
			printf("Please any key to exit...\n");
			exit(1);
		}

		// копии объектов
		// 9 ноября 2016 Копии объектов теперь хранятся централизованно:
		// см. LEVEL_ADDITIONAL_DATA_BUFER milu_gl_buffer; in ilut.c module.
		//milu2.alu_copy = new doublereal[milu2.iwk + 2]; // +2 запас по памяти.
		//milu2.jlu_copy = new integer[milu2.iwk + 2];
		//milu2.ju_copy = new integer[maxelm_plus_maxbound + 2];
		//milu2.b_copy = new doublereal[maxelm_plus_maxbound + 2];
		char c8[13] = "milu2.b_copy";
		milu2.b_copy = (doublereal*)malloc((maxelm_plus_maxbound + 2) * sizeof(doublereal));
		handle_error(milu2.b_copy, c8, c1, (maxelm_plus_maxbound + 2));
		//milu2.x_copy = new doublereal[maxelm_plus_maxbound + 2];
		milu2.x_copy = (doublereal*)malloc((maxelm_plus_maxbound + 2) * sizeof(doublereal));
		char c9[13] = "milu2.x_copy";
		handle_error(milu2.x_copy, c9, c1, (maxelm_plus_maxbound + 2));
		// Выделение оперативной памяти может быть 
		// повторным, поэтому необходимо освободить оперативную память,
		// перед повторным выделением.
		if (milu2.zbuf != nullptr) {
			//delete[] milu2.zbuf;
			free(milu2.zbuf);
			milu2.zbuf = nullptr;
		}
		//milu2.zbuf = new doublereal[maxelm_plus_maxbound + 2];
		milu2.zbuf = (doublereal*)malloc((maxelm_plus_maxbound + 2) * sizeof(doublereal));
		char c10[11] = "milu2.zbuf";
		handle_error(milu2.zbuf, c10, c1, (maxelm_plus_maxbound + 2));
		// Выделение оперативной памяти может быть 
		// повторным, поэтому необходимо освободить оперативную память,
		// перед повторным выделением.
		if (milu2.zbuf2 != nullptr) {
			//delete[] milu2.zbuf2;
			free(milu2.zbuf2);
			milu2.zbuf2 = nullptr;
		}
		//milu2.zbuf2 = new doublereal[maxelm_plus_maxbound + 2];
		milu2.zbuf2 = (doublereal*)malloc((maxelm_plus_maxbound + 2) * sizeof(doublereal));
		char c111[12] = "milu2.zbuf2";
		handle_error(milu2.zbuf2, c111, c1, (maxelm_plus_maxbound + 2));
		for (integer i_63 = 0; i_63 < (maxelm_plus_maxbound + 2); ++i_63) {
			milu2.zbuf2[i_63] = 0.0;
		}

		//if ((milu2.alu_copy == nullptr) || (milu2.jlu_copy == nullptr) ||
		//  (milu2.ju_copy == nullptr) || (milu2.b_copy == nullptr) || 
		//  (milu2.x_copy == nullptr) || (milu2.zbuf==nullptr)||(milu2.zbuf2==nullptr)) {
		if ((milu2.b_copy == nullptr) || (milu2.x_copy == nullptr) || 
			(milu2.zbuf == nullptr) || (milu2.zbuf2 == nullptr)) {
			// недостаточно памяти на данном оборудовании.
			printf("Problem: not enough memory on your equipment...\n");
			printf("Please any key to exit...\n");
			exit(1);
		}

		//milu2.lfil = 2;
		integer ierr = 0;
		iluk_(maxelm_plus_maxbound, milu2.val, milu2.col_ind, milu2.row_ptr, my_amg_manager.lfil/*milu2.lfil*/, milu2.alu, milu2.jlu, milu2.ju, milu2.levs, milu2.iwk, milu2.w, milu2.jw, ierr);
		if (ierr != 0) {
#if doubleintprecision == 1
		    printf("ierr=%lld\n", ierr);
#else
		    printf("ierr=%d\n", ierr);
#endif
			//system("pause");
			system("PAUSE");
			exit(1);
		}

		if ((ierr == -2) || (ierr == -3)) {

			integer ipassage = 1;
			do {
				printf("\nPlease WAIT... ... ...\n");

				// задаче не хватило памяти, значит нужно перевыделить !
				//if (milu2.alu != nullptr) delete[] milu2.alu;
				//if (milu2.jlu != nullptr) delete[] milu2.jlu;
		     	//if (milu2.levs != nullptr) delete[] milu2.levs;
				if (milu2.alu != nullptr) free(milu2.alu);
				if (milu2.jlu != nullptr) free(milu2.jlu);
				if (milu2.levs != nullptr) free(milu2.levs);
				//if (milu2.alu_copy != nullptr) delete[] milu2.alu_copy;
				//if (milu2.jlu_copy != nullptr) delete[] milu2.jlu_copy;

				// инициализация !
				milu2.alu = nullptr;
				milu2.jlu = nullptr;
				milu2.levs = nullptr;
				//milu2.alu_copy = nullptr;
				//milu2.jlu_copy = nullptr;
				integer ilevel_loc = ilevel;
				if (ilevel == 0) ilevel_loc = 1;
				integer lfil = my_amg_manager.lfil;
				if (lfil == 0) {
					milu2.iwk = static_cast<integer>((lfil + 1) * (n + 2 * maxelm_plus_maxbound + 2) + (4 + 3 * ipassage) * maxelm_plus_maxbound);
				}
				else if (lfil == 1) {
					milu2.iwk = static_cast<integer>((2.3 * ilevel_loc * lfil + 1) * (n + 2 * maxelm_plus_maxbound + 2) + (4 + 3 * ipassage) * maxelm_plus_maxbound);
				}
				else if (lfil == 2) {
					milu2.iwk = static_cast<integer>((2.6 * ilevel_loc * lfil + 1) * (n + 2 * maxelm_plus_maxbound + 2) + (4 + 3 * ipassage) * maxelm_plus_maxbound);
				}
				else if (lfil == 3) {
					milu2.iwk = static_cast<integer>((3.2 * ilevel_loc * lfil + 1) * (n + 2 * maxelm_plus_maxbound + 2) + (4 + 3 * ipassage) * maxelm_plus_maxbound);
				}
				else if (lfil == 4) {
					milu2.iwk = static_cast<integer>((3.9 * ilevel_loc * lfil + 1) * (n + 2 * maxelm_plus_maxbound + 2) + (4 + 3 * ipassage) * maxelm_plus_maxbound);
				}
				else if (lfil == 5) {
					milu2.iwk = static_cast<integer>((4.4 * ilevel_loc * lfil + 1) * (n + 2 * maxelm_plus_maxbound + 2) + (4 + 3 * ipassage) * maxelm_plus_maxbound);
				}
				else if (lfil == 6) {
					milu2.iwk = static_cast<integer>((5.0 * ilevel_loc * lfil + 1) * (n + 2 * maxelm_plus_maxbound + 2) + (4 + 3 * ipassage) * maxelm_plus_maxbound);
				}
				else if (lfil >= 7) {
					milu2.iwk = static_cast<integer>((5.6 * ilevel_loc * lfil + 1) * (n + 2 * maxelm_plus_maxbound + 2) + (4 + 3 * ipassage) * maxelm_plus_maxbound);
				}

				//milu2.alu = new doublereal[milu2.iwk + 2]; // +2 запас по памяти.
				milu2.alu = (doublereal*)malloc((milu2.iwk + 2) * sizeof(doublereal));
				char c12[11] = "milu2.alu";
				handle_error(milu2.alu, c12, c1, (milu2.iwk + 2));
				//milu2.jlu = new integer[milu2.iwk + 2];
				milu2.jlu = (integer*)malloc((milu2.iwk + 2) * sizeof(integer));
				char c13[11] = "milu2.jlu";
				handle_error(milu2.jlu, c13, c1, (milu2.iwk + 2));
				//milu2.levs = new integer[milu2.iwk + 2]; // уровень.
				milu2.levs = (integer*)malloc((milu2.iwk + 2) * sizeof(integer));
				char c14[11] = "milu2.levs";
				handle_error(milu2.levs, c14, c1, (milu2.iwk + 2));
				//milu2.alu_copy = new doublereal[milu2.iwk + 2]; // +2 запас по памяти.
				//milu2.jlu_copy = new integer[milu2.iwk + 2];

				//if ((milu2.alu_copy != nullptr) && (milu2.jlu_copy != nullptr) && (milu2.alu != nullptr) && (milu2.jlu != nullptr) && (milu2.levs != nullptr)) {
				if ((milu2.alu != nullptr) && (milu2.jlu != nullptr) && (milu2.levs != nullptr)) {
					iluk_(maxelm_plus_maxbound, milu2.val, milu2.col_ind, milu2.row_ptr, lfil, milu2.alu, milu2.jlu, milu2.ju, milu2.levs, milu2.iwk, milu2.w, milu2.jw, ierr);
				}
				else {
					// недостаточно памяти на данном оборудовании.
					ipassage = 4;
					printf("Problem: not enough memory on your equipment...\n");
					printf("Please any key to exit...\n");
					exit(1);
				}

				ipassage++;
				// Пока просто контролируем, в добавок хорошо бы выводить информацию на каком уровне это произошло.
#if doubleintprecision == 1
				printf("control memory allocation for ilu2 in RUMBA0.14 in level=%lld\n", ilevel);
#else
				printf("control memory allocation for ilu2 in RUMBA0.14 in level=%d\n", ilevel);
#endif

				system("PAUSE");
				
			} while ((ierr != 0) && (ipassage < 4));

			if (ipassage == 4) {
				printf("Error memory alloc !!!\n");
				printf("failed to obtain an expansion for the 4 approaches...\n");
			    printf("Please any key to exit...\n");
				//system("pause");
				system("PAUSE");
				exit(1);
			}
		}

		for (integer k = 0; k < n; ++k) if (milu2.col_ind[k] == (-1)) {
			printf("Error equation3D to CRS.\n");
			//system("pause");
			system("pause");
		}

		// Нумерация начинается с единицы (это важно).
		for (integer k = 0; k < maxelm_plus_maxbound; ++k) {
			if (milu2.val[milu2.row_ptr[k] - 1] < nonzeroEPS) {
#if doubleintprecision == 1
		    	printf("negativ diagonal element equation3DtoCRS %lld\n", k);
#else
				printf("negativ diagonal element equation3DtoCRS %d\n", k);
#endif

				//system("pause");
				system("pause");
		    }
	    }

		if (1 || milu2.lfil == 0) {
			// экономим оперативную память.
			//if (milu2.val != nullptr) delete[] milu2.val;
			//if (milu2.col_ind != nullptr) delete[] milu2.col_ind;
			//if (milu2.row_ptr != nullptr) delete[] milu2.row_ptr;
			if (milu2.val != nullptr) free(milu2.val);
			if (milu2.col_ind != nullptr) free(milu2.col_ind);
			if (milu2.row_ptr != nullptr) free(milu2.row_ptr);
			milu2.val = nullptr;
			milu2.col_ind = nullptr;
			milu2.row_ptr = nullptr;
		}

		calculation_main_end_time = clock(); // момент окончания счёта.
		if (!(steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::NETWORK_T_UNSTEADY)) {
			printf("%6u \n", calculation_main_end_time - calculation_main_start_time);
		}
    }

} // equation3DtoCRSRUMBA1_amg1r5

  // Мы усилим РУМБА0.14 алгоритм ILU0 предобуславливанием (сглаживателем) на каждом уровне вложенности.
  // Преобразует equation3D  формат хранения в CRS формат.
  // Цель написания этого преобразователя: экономия оперативной памяти компьютера.
  // Т.к. формат SIMPLESPARSE требует слишком много памяти.
  // Для многих задач сходимость мульгрида на базе Зейдель сглаживателя недостаточна. 
  // Судя по документации ILU0 smoother more effective по сравнению с Зейдель сглаживателем.
  // ILU0 - универсален и не накладывает ограничений на структуру матрицы СЛАУ. ILU0 потребляет мало пользователей.
void equation3DtoCRSRUMBA0(LEVEL_ADDITIONAL_DATA0 &milu0,
	bool ballocmemory, Ak2 &Amat, integer istartq, integer iendq,
	integer* &row_ptr_start, integer* &row_ptr_end, integer iadd, integer ilevel) {


	// iadd_now=n_a[0]+...+n_a[ilevel_detector-1];
	//equation3DtoCRSRUMBA0(milu0[ilevel_detector], true,
	//Amat, 1, n_a[ilevel_detector], row_ptr_start, row_ptr_end, iadd_now, ilevel_detector);

	// Если ballocmemory равен true то происходит выделение памяти.
	const bool bonly_negative_connections = false;
	bool flag = true;
	integer n = 0; // число ненулевых элементов
	integer maxelm_plus_maxbound = 0; // Число уравнений в СЛАУ.

	const doublereal nonzeroEPS = 1e-37; // для отделения вещественного нуля

	integer startpos = istartq + iadd;
	integer endpos = iendq + iadd;


	// подсчёт количества ненулевых элементов
	// во всех точках расчётной области как внутренних так и граничных.
	for (integer ii = startpos; ii <= endpos; ++ii) {
		maxelm_plus_maxbound++;
		integer istr = ii - iadd;
		if (fabs(Amat.aij[row_ptr_start[ii]]) < 1.0e-20)
		{
#if doubleintprecision == 1
			printf("zero diagonal coefficient in level = %lld, istr=%lld\n", ilevel, istr);
#else
			printf("zero diagonal coefficient in level = %d, istr=%d\n", ilevel, istr);
#endif

			//system("pause");
			system("PAUSE");
			exit(1);
		}
		doublereal ap_now = 1.0 / Amat.aij[row_ptr_start[ii]];

		//if (fabs(ap_now)> 1e10*nonzeroEPS) n++; // диагональный элемент
		if (ap_now > nonzeroEPS) n++; // Диагональный элемент.
		else {
			// 5 августа 2016. 
			flag = false;
			printf("internal zero diagonal element.\n");
			//printf("ae=%e aw=%e an=%e as=%e at=%e ab=%e sum_nb=%e", sl[k].ae, sl[k].aw, sl[k].an, sl[k].as, sl[k].at, sl[k].ab, sl[k].ae + sl[k].aw + sl[k].an + sl[k].as + sl[k].at + sl[k].ab);
			if (ap_now < 0.0) {
				printf("found negativ diagonal coefficient=%e...\n", ap_now);
			}
			printf("fatal error equation3DtoCRS in RUMBA...\n");
			//system("pause");
			system("PAUSE");
			exit(1);
			//n++;
			//sl[k].ap = fabs(sl[k].ae) + fabs(sl[k].aw) + fabs(sl[k].an) + fabs(sl[k].as) + fabs(sl[k].at) + fabs(sl[k].ab);
		}

		integer is1 = row_ptr_start[ii] + 1;
		integer is2 = row_ptr_end[ii];

		for (integer ii1 = is1; ii1 <= is2; ++ii1)
		{
			if (bonly_negative_connections) {
				if (Amat.aij[ii1] < -nonzeroEPS) {
					n++;
				}
			}
			else {
				if (fabs(Amat.aij[ii1]) > nonzeroEPS) {
					n++;
				}
				else {
#if doubleintprecision == 1
					printf("igmoring a[%d][%d]=%e istr=%lld\n", Amat.i[ii1], Amat.j[ii1], fabs(Amat.aij[ii1]), istr);
#else
					printf("igmoring a[%d][%d]=%e istr=%d\n", Amat.i[ii1], Amat.j[ii1], fabs(Amat.aij[ii1]), istr);
#endif

				}
			}
		}
	}

#if doubleintprecision == 1
	printf("nnz=%lld, %lld startpos=%lld endpos=%lld\n", n, row_ptr_end[endpos] - row_ptr_start[startpos] + 1, startpos, endpos);
#else
	printf("nnz=%d, %d startpos=%d endpos=%d\n", n, row_ptr_end[endpos] - row_ptr_start[startpos] + 1, startpos, endpos);
#endif

	system("PAUSE");

	if (flag) {
		// memory +15N
		// Теперь выделение памяти будет происходить централизованно, вне данного кода.
		// Это сделано для кода BICGSTAB_internal3. дата изменения 12 апреля 2013.
		// Другой код, использующий equation3dtoCRS может оказаться неработоспособным после этого изменения.
		if (ballocmemory) {
			// Важно выделить память с запасом, т.к. одна и таже память используется и для компонент скорости и для поправки давления.
			//val = new doublereal[7 * (maxelm + maxbound) + 2 * maxbound + 2];
			//col_ind = new integer[7 * (maxelm + maxbound) + 2 * maxbound + 2];
			milu0.val = new doublereal[n + 2];
			milu0.col_ind = new integer[n + 2];
			//row_ptr = new integer[(maxelm + maxbound) + 1];
			milu0.row_ptr = new integer[maxelm_plus_maxbound + 1];
			if ((milu0.val == nullptr) || (milu0.col_ind == nullptr) || (milu0.row_ptr == nullptr)) {
				// недостаточно памяти на данном оборудовании.
				printf("Problem: not enough memory on your equipment for ILU0 decomposition...\n");
				printf("Please any key to exit...\n");
				exit(1);
			}
		}


		// инициализация
#pragma omp parallel for
		for (integer k = 0; k<(n); ++k) {
			milu0.val[k] = 0.0;
			milu0.col_ind[k] = -1;
		}
#pragma omp parallel for
		for (integer k = 0; k <= (maxelm_plus_maxbound); ++k) {
			milu0.row_ptr[k] = n; // присваиваем количество ненулевых элементов плюс 1 с учётом того что нумерация массива начинается с 0
		}

		// Быстрая Сортировка Ч. Хоара.
		// упорядочивание по строкам
		//QuickSort(...); не требуется,
		// т.к. сама структура хранения 
		// подразумевает упорядочивание по строкам.

		/*
		// заполнение разреженной матрицы
		for (integer k=0; k<M.n; ++k) {
		val[k]=M.a[k].aij;
		col_ind[k]=M.a[k].j;
		row_ptr[M.a[k].i]=min(k,row_ptr[M.a[k].i]);
		}
		*/
		integer ik = 0; // счётчик ненулевых элементов СЛАУ

						// для всех узлов расчётной области как внутренних так и внешних:
		for (integer ii = startpos; ii <= endpos; ++ii) {
			integer istr = ii - iadd;
			integer k = istr - 1; // нумерация с нуля.
			doublereal ap_now = 1.0 / Amat.aij[row_ptr_start[ii]];


			// Диагональный коэффициент.
			if (fabs(ap_now) > nonzeroEPS) {
				milu0.val[ik] = ap_now;
				milu0.col_ind[ik] = istr - 1; // Нуменрация начинается с нуля.
				milu0.row_ptr[k] = (ik < milu0.row_ptr[k] ? ik : milu0.row_ptr[k]);
				ik++;
			}

			integer is1 = row_ptr_start[ii] + 1;
			integer is2 = row_ptr_end[ii];

			for (integer ii1 = is1; ii1 <= is2; ++ii1)
			{
				if (bonly_negative_connections) {
					if (Amat.aij[ii1] < -nonzeroEPS) {
						milu0.val[ik] = Amat.aij[ii1]; // уже с учётом знака.
						milu0.col_ind[ik] = Amat.j[ii1] - 1; // Нумерация начинается с нуля.
						milu0.row_ptr[k] = (ik < milu0.row_ptr[k] ? ik : milu0.row_ptr[k]);
						ik++;
					}
				}
				else {
					if (fabs(Amat.aij[ii1]) > nonzeroEPS) {
						milu0.val[ik] = Amat.aij[ii1]; // уже с учётом знака.
						milu0.col_ind[ik] = Amat.j[ii1] - 1; // Нумерация начинается с нуля.
						milu0.row_ptr[k] = (ik < milu0.row_ptr[k] ? ik : milu0.row_ptr[k]);
						ik++;
					}
				}
			}
		}

		// в каждой строке элементы отсортированы по номерам столбцов:
#pragma omp parallel for
		for (integer k = 0; k<(maxelm_plus_maxbound); ++k) QuickSortCSIR(milu0.col_ind, milu0.val, milu0.row_ptr[k] + 1, milu0.row_ptr[k + 1] - 1);

#if doubleintprecision == 1
		//printf("n==%lld row_ptr=%lld\n", n, milu0.row_ptr[maxelm_plus_maxbound]);
#else
		//printf("n==%d row_ptr=%d\n", n, milu0.row_ptr[maxelm_plus_maxbound]);
#endif

		//system("pause");
#pragma omp parallel for
		for (integer i = 0; i<milu0.row_ptr[maxelm_plus_maxbound]; ++i) {
			milu0.col_ind[i] = milu0.col_ind[i] + 1;
		}
#pragma omp parallel for
		for (integer i = 0; i<maxelm_plus_maxbound + 1; ++i) {
			milu0.row_ptr[i] = milu0.row_ptr[i] + 1;
		}

		/*
		FILE *fp;
		errno_t err;
		// создание файла для записи.
		if ((err = fopen_s( &fp, "matr.txt", "w")) != 0) {
		printf("Create File Error\n");
		}
		else {
		#if doubleintprecision == 1
		// debug
		for (k=0; k<=maxelm+maxbound; ++k) {
		fprintf(fp,"%lld ",milu0.row_ptr[k]);
		}
		fprintf(fp,"\n");
		for (k=0; k<milu0.row_ptr[maxelm+maxbound]; ++k) {
		fprintf(fp, "%e %lld\n",milu0.val[k],milu0.col_ind[k]);
		}
		#else
		// debug
		for (k=0; k<=maxelm+maxbound; ++k) {
		fprintf(fp,"%d ",milu0.row_ptr[k]);
		}
		fprintf(fp,"\n");
		for (k=0; k<milu0.row_ptr[maxelm+maxbound]; ++k) {
		fprintf(fp, "%e %d\n",milu0.val[k],milu0.col_ind[k]);
		}
		#endif



		fclose(fp);
		}
		printf("ready");
		system("pause");
		*/
		//doublereal radd = 0.0;
		/*
		if (ilevel == 0) {
		radd = 1.8;
		}
		else if (ilevel == 1) {
		radd = 4.2;
		}
		else if ((ilevel >= 2) && (ilevel <= 5)) {
		radd = 6.5;
		}
		else if (ilevel == 6) {
		radd = 7.0;
		}
		else {
		radd = 8.0;
		}
		*/

		// ILU0 не требует памяти больше чем размер исходных матриц. 
		// Новых ненулевых элементов для него не появляется.
		milu0.iwk = milu0.row_ptr[maxelm_plus_maxbound];
		milu0.alu = new doublereal[milu0.row_ptr[maxelm_plus_maxbound] + 2];
		milu0.jlu = new integer[milu0.row_ptr[maxelm_plus_maxbound] + 2];

		milu0.ju = new integer[maxelm_plus_maxbound + 2];


		milu0.iw = new integer[maxelm_plus_maxbound + 2]; // +2 запас по памяти.

		if ((milu0.alu == nullptr) || (milu0.jlu == nullptr) || (milu0.ju == nullptr) || (milu0.iw == nullptr)) {
			// недостаточно памяти на данном оборудовании.
			printf("Problem: not enough memory on your equipment for ILU0 decomposition: alu, jlu, ju, iw...\n");
			printf("Please any key to exit...\n");
			exit(1);
		}

		// копии объектов

		// ILU0 не требует памяти больше чем размер исходных матриц. 
		// Новых ненулевых элементов для него не появляется.
		milu0.alu_copy = new doublereal[milu0.row_ptr[maxelm_plus_maxbound] + 2];
		milu0.jlu_copy = new integer[milu0.row_ptr[maxelm_plus_maxbound] + 2];

		milu0.ju_copy = new integer[maxelm_plus_maxbound + 2];

		milu0.b_copy = new doublereal[maxelm_plus_maxbound + 2];
		milu0.x_copy = new doublereal[maxelm_plus_maxbound + 2];
		milu0.zbuf = new doublereal[maxelm_plus_maxbound + 2];
		milu0.zbuf2 = new doublereal[maxelm_plus_maxbound + 2];

		if ((milu0.alu_copy == nullptr) || (milu0.jlu_copy == nullptr) || (milu0.ju_copy == nullptr) || (milu0.b_copy == nullptr) || (milu0.x_copy == nullptr) || (milu0.zbuf == nullptr) || (milu0.zbuf2 == nullptr)) {
			// недостаточно памяти на данном оборудовании.
			printf("Problem: not enough memory on your equipment for ILU0 decomposition: alu_copy, jlu_copy, ju_copy etc...\n");
			printf("Please any key to exit...\n");
			exit(1);
		}

		// Построение ILU0 декомпозиции.
		integer ierr = 0;
		ilu0_(maxelm_plus_maxbound, milu0.val, milu0.col_ind, milu0.row_ptr, milu0.alu, milu0.jlu, milu0.ju, milu0.iw, ierr);
		if (ierr != 0) {
#if doubleintprecision == 1
			printf("ierr=%lld\n", ierr);
#else
			printf("ierr=%d\n", ierr);
#endif

			//system("pause");
			system("PAUSE");
			exit(1);
		}

		if ((ierr == -2) || (ierr == -3)) {

			integer ipassage = 1;
			do {
				printf("\nPlease WAIT... ... ...\n");

				// задаче не хватило памяти, значит нужно перевыделить !
				if (milu0.alu != nullptr) delete milu0.alu;
				if (milu0.jlu != nullptr) delete milu0.jlu;
				if (milu0.alu_copy != nullptr) delete milu0.alu_copy;
				if (milu0.jlu_copy != nullptr) delete milu0.jlu_copy;

				// инициализация !
				milu0.alu = nullptr;
				milu0.jlu = nullptr;
				milu0.alu_copy = nullptr;
				milu0.jlu_copy = nullptr;

				// При использовании ILU0 новой памяти не требуется должно хватать размера исходной матрицы.
				integer iwk = milu0.row_ptr[maxelm_plus_maxbound] + 2; // ИЗМЕНИТЬ ЭТУ СТРОКУ для увеличения памяти.
				milu0.iwk = iwk - 2;
				milu0.alu = new doublereal[iwk + 2]; // +2 запас по памяти.
				milu0.jlu = new integer[iwk + 2];
				milu0.alu_copy = new doublereal[iwk + 2]; // +2 запас по памяти.
				milu0.jlu_copy = new integer[iwk + 2];

				if ((milu0.alu_copy != nullptr) && (milu0.jlu_copy != nullptr) && (milu0.alu != nullptr) && (milu0.jlu != nullptr)) {
					ilu0_(maxelm_plus_maxbound, milu0.val, milu0.col_ind, milu0.row_ptr, milu0.alu, milu0.jlu, milu0.ju, milu0.iw, ierr);
				}
				else {
					// недостаточно памяти на данном оборудовании.
					ipassage = 4;
					printf("Problem: not enough memory on your equipment...\n");
					printf("Please any key to exit...\n");
					exit(1);
				}

				ipassage++;
				// Пока просто контролируем, в добавок хорошо бы выводить информацию на каком уровне это произошло.
#if doubleintprecision == 1
				printf("control memory allocation for ilu0 in RUMBA0.14 in level=%lld\n", ilevel);
#else
				printf("control memory allocation for ilu0 in RUMBA0.14 in level=%d\n", ilevel);
#endif

				system("PAUSE");
			} while ((ierr != 0) && (ipassage<4));

			if (ipassage == 4) {
				printf("Error memory alloc !!!\n");
				printf("failed to obtain an expansion for the 4 approaches...\n");
				printf("Please any key to exit...\n");
				//system("pause");
				system("PAUSE");
				exit(1);
			}
		}

		// По идее здесь уже можно освободить ОЗУ от milu0.val, milu0.col_ind, milu0.row_ptr.
		// 7 ноября 2016.

	}

	if (!flag) {
		printf("Error equation 3D to CRS: zero diagonal element...\n");
		//system("pause");
		system("pause");
	}

	for (integer k = 0; k<n; ++k) if (milu0.col_ind[k] == (-1)) {
		printf("Error equation3D to CRS.\n");
		//system("pause");
		system("pause");
	}

	// Нумерация начинается с единицы (это важно).
	for (integer k = 0; k<maxelm_plus_maxbound; ++k) {
		if (milu0.val[milu0.row_ptr[k] - 1]<nonzeroEPS) {
#if doubleintprecision == 1
			printf("negativ diagonal element equation3DtoCRS %lld\n", k);
#else
			printf("negativ diagonal element equation3DtoCRS %d\n", k);
#endif

			//system("pause");
			system("pause");
		}
	}


} // equation3DtoCRSRUMBA0

void nested_desection_patch(Ak2 &Amat, integer isize_na, bool* &nested_desection,
	integer* &row_ptr_start, integer* &row_ptr_end, integer iadd_now) {

	// Данный метод самописный и очень хреновый. Его не рекомендуется использовать в частности из-за граничных узлов.
	// Лучше производить разбиения по геометрическому признаку. Т.е. если есть координаты узлов, то слева выделить всё с координатами меньше
	// чем центральная координата, а справа всё что больше. По моему так делают в Open Foam X.

	// Описание способа работы.
	// Метод требует подачи на вход матрицы в Ak1 формате.
	// Данная функция выдаёт булевый вектор nested desection [1..isize_na]:
	// Ячейки равные false относятся к границе сепаратора (разделителя двух подобластей).
	// Введем число integer middle = static_cast<integer>(0.5*(1 + isize_na));
	// Ячейки помеченные true индекс которых [1..middle] Могут выполнятся первым потоком.
	// Ячейки помеченные true индекс которых [middle+1..isize_na] Могут выполнятся вторым потоком независимо от первого.

	// Готовим nested desection
	// для двух потоков.
	// Самая подробная матрица 0.
	// nested_desection0
	bool *nested_desection1 = nullptr;
	nested_desection1 = (bool*)malloc((isize_na + 1) * sizeof(bool));
	char c1[19] = "nested_desection1";
	char c2[24] = "nested_desection_patch";
	char c3[19] = "nested_desection2";
	handle_error(nested_desection1, c1, c2, (isize_na + 1));
	bool *nested_desection2 = nullptr;
	nested_desection2 = (bool*)malloc((isize_na + 1) * sizeof(bool));
	handle_error(nested_desection2, c3, c2, (isize_na + 1));

	// инициализация.
	for (integer i_73 = 1; i_73 <= isize_na; ++i_73) {
		nested_desection1[i_73] = false;
		nested_desection2[i_73] = false;
		nested_desection[i_73] = true;
	}
	integer middle = (1 + isize_na)/2;
	for (integer i_72 = 1; i_72 <= middle; ++i_72) {
		for (integer i_82 = row_ptr_start[i_72 + iadd_now]; i_82 <= row_ptr_end[i_72 + iadd_now]; ++i_82) {
			nested_desection1[Amat.j[i_82]] = true;
		}
	}
	for (integer i_72 = middle + 1; i_72 <= isize_na; ++i_72) {
		for (integer i_82 = row_ptr_start[i_72 + iadd_now]; i_82 <= row_ptr_end[i_72 + iadd_now]; ++i_82) {
			nested_desection2[Amat.j[i_82]] = true;
		}
	}
	// На пересечении двух множеств находятся узлы сепаратора.
	for (integer i_73 = 1; i_73 <= isize_na; ++i_73) {
		if ((nested_desection1[i_73]) && (nested_desection2[i_73])) {
			nested_desection[i_73] = false;
		}
	}
	// запоминаем узлов из сепаратора.
	for (integer i_73 = 1; i_73 <= isize_na; ++i_73) {
		nested_desection1[i_73] = nested_desection[i_73];
	}
	// Если какой-либо сосед принадлежит сепаратору то и сам узел принадлежит сепаратору.
	for (integer i_72 = 1; i_72 <= isize_na; ++i_72) {
		for (integer i_82 = row_ptr_start[i_72 + iadd_now]; i_82 <= row_ptr_end[i_72 + iadd_now]; ++i_82) {
			if (!nested_desection1[Amat.j[i_82]]) nested_desection[Amat.i[i_82]] = false;
		}
	}

	free(nested_desection1);
	free(nested_desection2);
} // nested_desection_patch

  // residual.
  // 3 jan 2016 ещё быстрее.
  // 9 september 2015.
  // q - quick.
template <typename doublerealT1, typename doublerealT2>
void residualq2(Ak2 &Amat, integer istartq, integer iendq, doublerealT1* &x, doublerealT1* &b, integer* &row_ptr_start, integer* &row_ptr_end, integer iadd, doublerealT2* &residual, doublerealT2* &my_diag)
{
	// istart - начальная позиция ненулевых элементов в матрице А.
	// iend - конечная позиция ненулевых элементов в матрице А.
	integer startpos = istartq + iadd;
	integer endpos = iendq + iadd;
	doublerealT2 dsum = 0.0;

#pragma omp parallel for 
	for (integer ii = startpos; ii <= endpos; ++ii) {
		integer istr = ii - iadd;
		residual[istr] = b[istr];
		for (integer ii1 = row_ptr_start[ii] + 1; ii1 <= row_ptr_end[ii]; ++ii1) {
			if (1 || (Amat.aij[ii1] < 0.0)) {
				residual[istr] += -Amat.aij[ii1] *x[Amat.j[ii1]];
			}
			else {
				// 16 августа 2016.
				// не работает.
				//dsum += Amat.aij[ii1];
			}
		}
		residual[istr] -= (1.0 / Amat.aij[row_ptr_start[ii]])*x[istr];  // 1 april 2017
		//residual[istr] -= (my_diag[Amat.i[row_ptr_start[ii]]] + dsum)* x[istr]; // 3 jan 2016
	}


} // residualq2

 // residual.
  // 3 jan 2016 ещё быстрее.
  // 9 september 2015.
  // q - quick.
template <typename doublerealT1, typename doublerealT2>
void residualq2/*_old_version*/(Ak2 &Amat, integer istartq, integer iendq, doublerealT1* &x, doublerealT1* &b,
	 integer* &row_ptr_start, integer* &row_ptr_end, integer iadd, doublerealT2* &residual,
	 doublerealT2* &my_diag, doublerealT2* &diag_minus_one)
{

#ifdef _OPENMP
	int i_my_num_core_parallelesation = omp_get_max_threads();
	omp_set_num_threads(number_cores());
#endif

	// istart - начальная позиция ненулевых элементов в матрице А.
	// iend - конечная позиция ненулевых элементов в матрице А.
	const integer startpos = istartq + iadd;
	const integer endpos = iendq + iadd;
	//doublerealT2 dsum = 0.0;

#pragma omp parallel for 
	for (integer ii = startpos; ii <= endpos; ++ii) {
		integer istr = ii - iadd;
		
		// домножаем на минус 1.0 если требуется.
		residual[istr] = diag_minus_one[istr]*b[istr];
		const integer istart = row_ptr_start[ii] + 1;
		const integer iend = row_ptr_end[ii];
		for (integer ii1 = istart; ii1 <= iend; ++ii1) {
			residual[istr] += -Amat.aij[ii1] * x[Amat.j[ii1]];
			//if (1 || (Amat.aij[ii1] < 0.0)) {
				//residual[istr] += -Amat.aij[ii1] * x[Amat.j[ii1]];
			//}
			//else {
				// 16 августа 2016.
				// не работает.
				//dsum += Amat.aij[ii1];
			//}
		}
		residual[istr] -= (1.0 / Amat.aij[row_ptr_start[ii]]) * x[istr];  // 1 april 2017
		//residual[istr] -= (my_diag[Amat.i[row_ptr_start[ii]]] + dsum)* x[istr]; // 3 jan 2016

		/*if (residual[istr] != residual[istr]) {
			std::cout << "residualq2\n";

			//for (integer ii1 = startpos; ii1 <= endpos; ++ii1) {
				//integer istr1 = ii1 - iadd;
				//std::cout << " istr1=" << istr1<<"/n";
				//system("pause");
			//}


			std::cout << "istr=" << istr<< " ap= " << (1.0 / Amat.aij[row_ptr_start[ii]]);
			std::cout << " 1/Ap= " << Amat.aij[row_ptr_start[ii]];
			std::cout << " b= " << b[istr]<< "  x[istr]="<<x[istr]<<"\n";
			for (integer ii1 = istart; ii1 <= iend; ++ii1) {
				std::cout << " aij= " << Amat.aij[ii1] << " j=" << Amat.j[ii1] << " x[j]=" << x[Amat.j[ii1]] <<" \n";
			}
			system("pause");
		}*/

	}


#ifdef _OPENMP
	omp_set_num_threads(i_my_num_core_parallelesation);
#endif

} // residualq2

// residual.
// 20.06.2021 Быстрее см. 
// https://habr.com/ru/post/309796/
// Оптимизация кода: процессор
  // 3 jan 2016 ещё быстрее.
  // 9 september 2015.
  // q - quick.
template <typename doublerealT1, typename doublerealT2>
void residualq2_optimizer(Ak2& Amat, integer istartq, integer iendq, doublerealT1*& x, doublerealT1*& b,
	integer*& row_ptr_start, integer*& row_ptr_end, integer iadd, doublerealT2*& residual,
	doublerealT2*& my_diag, doublerealT2*& diag_minus_one)
{
	// istart - начальная позиция ненулевых элементов в матрице А.
	// iend - конечная позиция ненулевых элементов в матрице А.
	const integer startpos = istartq + iadd;
	const integer endpos = iendq + iadd;
	//doublerealT2 dsum = 0.0;

#pragma omp parallel for 
	for (integer ii = startpos; ii <= endpos; ++ii) {
		integer istr = ii - iadd;
		// домножаем на минус 1.0 если требуется.
		residual[istr] = diag_minus_one[istr] * b[istr];
		const integer istart = row_ptr_start[ii] + 1;
		const integer iend = row_ptr_end[ii];

		doublerealT2 sum = 0.0;
		if (iend - istart + 1 >= 8) {
			integer ii1;
			doublerealT2 sum1 = 0.0;
			doublerealT2 sum2 = 0.0;
			doublerealT2 sum3 = 0.0;
			for (ii1 = istart; ii1+3 <= iend; ii1+=4) {
				sum += -Amat.aij[ii1] * x[Amat.j[ii1]];
				sum1 += -Amat.aij[ii1+1] * x[Amat.j[ii1+1]];
				sum2 += -Amat.aij[ii1+2] * x[Amat.j[ii1+2]];
				sum3 += -Amat.aij[ii1+3] * x[Amat.j[ii1+3]];
			}
			while (ii1 <= iend) {
				sum += -Amat.aij[ii1] * x[Amat.j[ii1]];
				++ii1;
			}
			sum += sum1 + sum2 + sum3;
		}
		else {

			for (integer ii1 = istart; ii1 <= iend; ++ii1) {
				sum += -Amat.aij[ii1] * x[Amat.j[ii1]];
			}
		}
		residual[istr] += sum - (1.0 / Amat.aij[row_ptr_start[ii]]) * x[istr];  // 1 april 2017
		//residual[istr] -= (my_diag[Amat.i[row_ptr_start[ii]]] + dsum)* x[istr]; // 3 jan 2016
	}


} // residualq2

  // residual.
  // Анализ структуры невязки в случае проблемы с солвером.
  // 3 jan 2016 ещё быстрее.
  // 9 september 2015.
  // q - quick.
template <typename doublerealT1, typename doublerealT2>
void residualq2_analysys(Ak2 &Amat, integer istartq, integer iendq,
	doublerealT1* &x, doublerealT1* &b, 
	integer* &row_ptr_start, integer* &row_ptr_end,
	integer iadd, doublerealT2* &residual, doublerealT2* &my_diag)
{
	// Анализ структуры невязки.
	// Мы делим вектор невязки на 20 равных частей и смотрим где невязка наиболее велика.



	printf("residual2 analysys start:\n");
	// istart - начальная позиция ненулевых элементов в матрице А.
	// iend - конечная позиция ненулевых элементов в матрице А.
	integer startpos = istartq + iadd;
	integer endpos = iendq + iadd;
	doublerealT2 dsum = 0.0;
	doublerealT2 statistics_sum = 0.0;
	integer icount_log = 1;
	integer i_srez = ((endpos - startpos) / 20);
	for (integer ii = startpos; ii <= endpos; ++ii) {
		integer istr = ii - iadd;
		residual[istr] = b[istr];
		for (integer ii1 = row_ptr_start[ii] + 1; ii1 <= row_ptr_end[ii]; ++ii1) {
			if (1 || (Amat.aij[ii1] < 0.0)) {
				residual[istr] += -Amat.aij[ii1] *x[Amat.j[ii1]];
			}
			else {
				// 16 августа 2016.
				// не работает.
				//dsum += Amat.aij[ii1];
			}
		}
		//residual[istr] += (-1.0 / Amat.aij[row_ptr_start[ii]])*x[istr];
		//residual[istr] -= (my_diag[Amat.i[row_ptr_start[ii]]] + dsum)* x[istr]; // 3 jan 2016
		residual[istr] -= (my_diag[istr] + dsum)* x[istr]; // 8 feb 2019
		statistics_sum += fabs(residual[istr]);
		if ((ii - startpos + 1) % i_srez == 0) {
#if doubleintprecision == 1
			printf("%lld %e ", icount_log++, statistics_sum);
#else
			printf("%d %e ", icount_log++, statistics_sum);
#endif

			statistics_sum = 0.0;
		}
	}
	printf("\n");

} // residualq2_analysys


// residual.
  // Анализ структуры невязки в случае проблемы с солвером.
  // 3 jan 2016 ещё быстрее.
  // 9 september 2015.
  // q - quick.
template <typename doublerealT1, typename doublerealT2>
void residualq2_analysys(Ak2& Amat, integer istartq, integer iendq, 
	doublerealT1*& x, doublerealT1*& b, 
	integer*& row_ptr_start, integer*& row_ptr_end, 
	integer iadd, doublerealT2*& residual,
	doublerealT2*& my_diag, doublerealT2*& diag_minus_one)
{
	// Анализ структуры невязки.
	// Мы делим вектор невязки на 20 равных частей и смотрим где невязка наиболее велика.



	printf("residual2 analysys start:\n");
	// istart - начальная позиция ненулевых элементов в матрице А.
	// iend - конечная позиция ненулевых элементов в матрице А.
	integer startpos = istartq + iadd;
	integer endpos = iendq + iadd;
	doublerealT2 dsum = 0.0;
	doublerealT2 statistics_sum = 0.0;
	integer icount_log = 1;
	integer i_srez = ((endpos - startpos) / 20);
	for (integer ii = startpos; ii <= endpos; ++ii) {
		integer istr = ii - iadd;
		residual[istr] = diag_minus_one[istr]*b[istr];
		for (integer ii1 = row_ptr_start[ii] + 1; ii1 <= row_ptr_end[ii]; ++ii1) {
			if (1 || (Amat.aij[ii1] < 0.0)) {
				residual[istr] += -Amat.aij[ii1] * x[Amat.j[ii1]];
			}
			else {
				// 16 августа 2016.
				// не работает.
				//dsum += Amat.aij[ii1];
			}
		}
		//residual[istr] += (-1.0 / Amat.aij[row_ptr_start[ii]])*x[istr];
		//residual[istr] -= (my_diag[Amat.i[row_ptr_start[ii]]] + dsum)* x[istr]; // 3 jan 2016
		residual[istr] -= (my_diag[istr] + dsum) * x[istr]; // 8 feb 2019
		statistics_sum += fabs(residual[istr]);
		if ((ii - startpos + 1) % i_srez == 0) {
#if doubleintprecision == 1
			printf("%lld %e ", icount_log++, statistics_sum);
#else
			printf("%d %e ", icount_log++, statistics_sum);
#endif

			statistics_sum = 0.0;
		}
	}
	printf("\n");

} // residualq2_analysys


#include "relaxation/spai0.cpp"
#include "relaxation/Runge_Kutt_3or5_order.cpp"
#include "relaxation/seidelqSor.cpp"

  

  // smoother.
  // 9 september 2015.
  // q - quick.
template <typename doublerealT>
void seidelq(Ak2 &Amat, integer istartq, integer iendq,
	doublerealT* &x, doublerealT* &b, 
	integer* &row_ptr_start, integer* &row_ptr_end, integer iadd)
{

	//seidelqstable<doublerealT>(Amat, istartq, iendq, x, b, row_ptr_start, row_ptr_end, iadd);

	// SOR!!!
	// 3 jan 2016
	//seidelqsor(Amat, istartq, iendq, x, b, row_ptr_start, row_ptr_end, iadd);
	// лучший выбор
	seidelqsor2(Amat, istartq, iendq, x, b, row_ptr_start, row_ptr_end, iadd);
	// 14 января 2016 смена направлений сканирования.
	// оказалось значительно хуже.
	//seidelqsor3(Amat, istartq, iendq, x, b, row_ptr_start, row_ptr_end, iadd);
	//seidelqsor2Pcpu(Amat, istartq, iendq, x, b, row_ptr_start, row_ptr_end, iadd);
	// early_naive_relaxation_method2.
    // Этот код (данная функция) больше не используется. Это одна из ранних
    // наивных попыток написать метод релаксации.
	//early_naive_relaxation_method2(Amat, istartq, iendq, x, b, row_ptr_start, row_ptr_end, iadd);
	// Экспериментальный вариант слишком медленный:
	//bool* flag = new bool[iendq - istartq + 2];
	//early_naive_relaxation_method(Amat, istartq, iendq, x, b, flag, iendq - istartq + 1);
	//delete[] flag;

} // seidelq

// smoother.
  // 9 september 2015.
  // q - quick.
template <typename doublerealT>
void seidelq(Ak2& Amat, integer istartq, integer iendq, doublerealT*& x, doublerealT*& b, integer*& row_ptr_start, integer*& row_ptr_end, integer iadd, doublerealT* &diag_minus_one)
{

	//seidelqstable<doublerealT>(Amat, istartq, iendq, x, b, row_ptr_start, row_ptr_end, iadd);

	// SOR!!!
	// 3 jan 2016
	//seidelqsor(Amat, istartq, iendq, x, b, row_ptr_start, row_ptr_end, iadd);
	// лучший выбор
	seidelqsor2(Amat, istartq, iendq, x, b, row_ptr_start, row_ptr_end, iadd, diag_minus_one);
	// 14 января 2016 смена направлений сканирования.
	// оказалось значительно хуже.
	//seidelqsor3(Amat, istartq, iendq, x, b, row_ptr_start, row_ptr_end, iadd);
	//seidelqsor2Pcpu(Amat, istartq, iendq, x, b, row_ptr_start, row_ptr_end, iadd);
	// early_naive_relaxation_method2.
	//early_naive_relaxation_method2(Amat, istartq, iendq, x, b, row_ptr_start, row_ptr_end, iadd);
	// Экспериментальный вариант слишком медленный:
	//bool* flag = new bool[iendq - istartq + 2];
	//early_naive_relaxation_method(Amat, istartq, iendq, x, b, flag, iendq - istartq + 1);
	//delete[] flag;

} // seidelq

  // smoother.
  // 5 января 2016 с использованием формулы из книги Патрика Роуча.
  // 9 september 2015.
  // q - quick.
  // 9 июня 2016. 
// early_naive_relaxation_method2.
// Этот код (данная функция) больше не используется. Это одна из ранних
// наивных попыток написать метод релаксации.
template <typename doublerealT>
void early_naive_relaxation_method2(Ak2 &Amat, integer istartq, integer iendq, doublerealT* &x, doublerealT* &b, integer* &row_ptr_start, integer* &row_ptr_end, integer iadd)
{
	// istart - начальная позиция ненулевых элементов в матрице А.
	// iend - конечная позиция ненулевых элементов в матрице А.
	// sor 1.855 использовать нельзя мы имеем расходимость вычислительного процесса.
	// Попробуем нижнюю релаксацию.
	// Расходимость имеем и при нижней релаксации. 0.8

	// BSKDmitrii
	// omega   iter  time,s
	// 1.0 106 43
	// 1.1 98 42
	// 1.15 94 40 best
	// 1.2 90 40
	// 1.225 413 1min 37s
	// 1.25 divergence detected
	// 1.3 divergence detected

	doublerealT omega = 1.0; // initialize.

							 // За подробностями смотри книгу Патрика Роуча стр. 183.
	doublerealT rn = static_cast<doublerealT>(iendq - istartq + 1);
	optimal_omega(rn, omega);
	omega = 1.0;

	// это нужне сделать единожды на фазе препроцессинга.
	doublerealT* diagonal = nullptr;
	diagonal = new doublerealT[iendq - istartq + 2];
	if (diagonal == nullptr) {
		// недостаточно памяти на данном оборудовании.
		printf("Problem: not enough memory on your equipment for diagonal in classical_ilu2...\n");
		printf("Please any key to exit...\n");
		//system("pause");
		system("pause");
		exit(1);
	}
	// initialization:
	for (integer i_6 = 0; i_6 < iendq - istartq + 2; ++i_6) diagonal[i_6] = 1.0;
	integer startpos = istartq + iadd;
	integer endpos = iendq + iadd;
	for (integer ii = startpos; ii <= endpos; ++ii) {
		{
			integer istr = ii - iadd;
			integer ic = ii;
			integer is1 = row_ptr_start[ii] + 1;
			integer is2 = row_ptr_end[ii];
			doublerealT ap = 0.0;
			doublerealT mult = 0.0;
			for (ic = is1; ic <= is2; ++ic) {
				if (Amat.j[ic] != istr) {
					if (Amat.j[ic] < istr) {
						// aij*aji/djj;
						// найти aji;
						doublerealT aji = 0.0;
						// линейный поиск очень медленный особенно на каждой итерации.
						//for (integer ii1 = startpos; ii1 <= endpos; ++ii1) {
						integer ii1 = iadd + Amat.j[ic];
						if (Amat.i[ii1] == Amat.j[ic]) {
							integer ic1 = ii1;
							integer is11 = row_ptr_start[ii1] + 1;
							integer is21 = row_ptr_end[ii1];
							for (ic1 = is11; ic1 <= is21; ++ic1) {
								if (Amat.j[ic1] != Amat.j[ic]) {
									if (Amat.j[ic1] == istr) {
										aji = Amat.aij[ic1];
									}
								}
							}
							break;
						}
						//}
						mult += (Amat.aij[ic] *aji) / diagonal[Amat.j[ic]];
					}
				}
			}
			ap = 1.0 / Amat.aij[is1 - 1];
			diagonal[istr] = ap - mult;
		}
	}

	//if (isorintmemo == iadd) {
	// Это точно не первый раз
	//bfirst = false;
	//}
	startpos = istartq + iadd;
	endpos = iendq + iadd;
	for (integer ii = startpos; ii <= endpos; ++ii) {
		integer istr = ii - iadd;
		doublerealT rold = x[istr];

		x[istr] = b[istr];

		doublerealT rsum = 0.0;
		integer is1 = row_ptr_start[ii] + 1;
		integer is2 = row_ptr_end[ii];
		// Распараллеливание почемуто тормозит очень сильно.
		//#pragma omp parallel for reduction(+:rsum)
		for (integer ii1 = is1; ii1 <= is2; ++ii1)
		{
			//x[istr] += -Amat.aij[ii1]*x[Amat.j[ii1]];
			integer ipos = Amat.j[ii1];
			rsum += -Amat.aij[ii1] *x[ipos];
		}
		x[istr] += rsum;
		//x[istr] *= Amat.aij[row_ptr_start[ii]];
		// incomplete lower upwer decomposition.
		x[istr] /= diagonal[istr];

		// Запускается только после первого раза сразу запускать нельзя
		// т.к. начальное значение не является адекватным.
		//x[istr] = omega*x[istr] + (1.0 - omega)*rold; // this is SOR
	}

	delete[] diagonal;

} // early_naive_relaxation_method2

  
  
// smoother.
// 5 jan 2016.
// 1 september 2015.
// 9 июня 2016.
// early_naive_relaxation_method.
// Этот код (данная функция) больше не используется. Это одна из ранних
// наивных попыток написать метод релаксации.
template <typename doublerealT>
void early_naive_relaxation_method(Ak2 &Amat, integer istart, integer iend,
	doublerealT* &x, doublerealT* &b, bool* &flag, integer n)
{

	for (integer i = 1; i <= n; ++i) {
		flag[i] = false;
	}
	doublerealT rn = 1.0;
	for (integer ii = istart; ii <= iend; ++ii) {
		if (flag[Amat.i[ii]] == false) {
			flag[Amat.i[ii]] = true;
			rn += 1.0;
		}
	}

	doublerealT omega = 1.0; // initialize.
	for (integer i = 1; i <= n; ++i) {
		flag[i] = false;
	}

	doublerealT* diagonal = new doublerealT[static_cast<integer>(rn)];
	for (integer ii = istart; ii <= iend; ++ii) {
		if (flag[Amat.i[ii]] == false) {
			integer istr = Amat.i[ii];
			integer ic = ii;
			doublerealT ap = 0.0;
			doublerealT mult = 0.0;
			x[istr] = b[istr];
			while ((ic <= iend) && (Amat.i[ic] == istr)) {
				if (Amat.j[ic] != istr) {
					if (Amat.j[ic] < istr) {
						// aij*aji/djj;
						// найти aji;
						doublerealT aji = 0.0;
						// Здесь нужно сделать хеш таблицу которая позволяет избавиться от линейного поиска.
						// линейный поиск очень медленный особенно на каждой итерации.
						for (integer ii1 = istart; ii1 <= iend; ++ii1) {
							if (Amat.i[ii1] == Amat.j[ic]) {
								integer ic1 = ii1;
								while ((ic1 <= iend) && (Amat.i[ic1] == Amat.j[ic])) {
									if (Amat.j[ic1] != Amat.j[ic]) {
										if (Amat.j[ic1] == istr) {
											aji = Amat.aij[ic1];
										}
									}
									ic1++;
								}
								break;
							}
						}
						mult += (Amat.aij[ic] *aji) / diagonal[Amat.j[ic]];
					}
				}
				else ap = Amat.aij[ic];
				ic++;
			}
			diagonal[istr] = ap - mult;
			flag[Amat.i[ii]] = true;
		}
	}


	// За подробностями смотри книгу Патрика Роуча стр. 183.
	optimal_omega(rn, omega);
	omega = 1.0;

	// istart - начальная позиция ненулевых элементов в матрице А.
	// iend - конечная позиция ненулевых элементов в матрице А.
	for (integer i = 1; i <= n; ++i) {
		flag[i] = false;
	}
	for (integer ii = istart; ii <= iend; ++ii) {
		if (flag[Amat.i[ii]] == false) {
			integer istr = Amat.i[ii];
			integer ic = ii;
			doublerealT ap = 0.0;
			doublerealT rold = x[istr];
			x[istr] = b[istr];
			while ((ic <= iend) && (Amat.i[ic] == istr)) {
				if (Amat.j[ic] != istr) {
					x[istr] += -Amat.aij[ic] *x[Amat.j[ic]];
				}
				else ap = Amat.aij[ic];
				ic++;
			}
			/*
			if (fabs(ap) < 1.0e-30) {
			#if doubleintprecision == 1
				printf("zero diagonal elements in string %lld",istr);
			#else
				printf("zero diagonal elements in string %d",istr);
			#endif

			system("pause");
			exit(1);
			}
			else */ {
			//x[istr] /= ap;
				x[istr] /= diagonal[istr];
				//x[istr] = omega*x[istr] + (1.0 - omega)*rold; // this is SOR
			}
			flag[Amat.i[ii]] = true;
		}
	}

	delete[] diagonal;

} // early_naive_relaxation_method




typedef struct TDATA_BASE_GMRES {
	doublereal* val;
	integer* col_ind;
	integer* row_ptr;

	TDATA_BASE_GMRES() {
		val = nullptr;
		col_ind = nullptr;
		row_ptr = nullptr;
	}
} DATA_BASE_GMRES;

DATA_BASE_GMRES* dbgmres_smoother = nullptr;

// 13.09.2017.
/*
* 30.12.2021
Да, можно использовать сглаживатели из подпространств Крылова в алгебраическом многосеточном методе,
но методы Крылова, как правило, не обладают хорошими сглаживающими свойствами. Это потому, что они
нацелены на весь спектр адаптивным способом, который минимизирует остаточную или подходящую норму
ошибки. Как правило, это будет включать в себя некоторые низкочастотные (длинноволновые) моды, с
которыми грубые сетки могли бы нормально работать. Сглаживатели Крылова также делают многосеточный
цикл нелинейным, поэтому, если мультисетка используется в качестве предварительного условия для 
внешнего метода Крылова, внешний метод должен быть «гибким» (например, GCR или FGMRES).

Использование сглаживателей Крылова также значительно увеличивает число точечных произведений,
которые должны быть вычислены, что становится существенным узким местом параллельно. Однако даже
с этими непривлекательными свойствами сглаживатели Крылова иногда полезны, особенно для сложных
задач, в которых хорошие операторы интерполяции недоступны.
*/
template <typename doublerealT>
void gmres_smoother(Ak2 &Amat, integer istartq, integer iendq, 
	doublerealT* &x, doublerealT* &b, 
	integer* &row_ptr_start, integer* &row_ptr_end, integer iadd, integer iVar, integer ilevel)
{

	if (my_amg_manager.b_ChebyshevSmoother) {
		// Выяснено что заходит сюда.

		// 26.01.2022
		// Сглаживатель Пафнутия Львовича Чебышева.
		// Никаких лишних аллокаций и освобождений оперативной памяти.
		// На примере механического расчета фермы в 137К ячеек обнаружено, 
		// что с аллокациями памяти и освобождениями памяти время расчёта составляет 10минут 21с 320мс.
		// А при работе на месте без дополнительных аллокаций и освобождения памяти время расчёта
		// фермы сократилось до 8 минут 27с 470мс. При степени полинома Чебышева 64 число итераций составляло 400.
		// При степени полинома Чебышева 25 время расчёта было 15мин, было сделано 1260 итераций V циклов. AMG Румба PMIS2, интерполяция amg1r5.

		// Без использования полиномов Чебышева время расчёта фермы составляло более часа на ilu сглаживателе, потребление оперативной памяти еле влезало в 32Гб машину.

		chebyshev(Amat, istartq, iendq, x, b, row_ptr_start, row_ptr_end, iadd, iVar, ilevel);
	}
	else {

		integer n = iendq - istartq + 1;

		// Разреженная матрица СЛАУ
		// в CRS формате.

		//doublerealT* val = nullptr;
		//integer* col_ind = nullptr, * row_ptr = nullptr;
		doublerealT* dX0 = new doublerealT[n];
		doublerealT* dV = new doublerealT[n];

		integer startpos = istartq + iadd;
		integer endpos = iendq + iadd;
#pragma omp parallel for
		for (integer ii = startpos; ii <= endpos; ++ii) {
			integer istr = ii - iadd;
			dX0[ii - startpos] = x[istr];
			dV[ii - startpos] = b[istr];
		}

		integer nnz = row_ptr_end[endpos] - row_ptr_start[startpos] + 1;
		if (dbgmres_smoother[ilevel].val == nullptr) {
			dbgmres_smoother[ilevel].val = new doublerealT[nnz];
			dbgmres_smoother[ilevel].col_ind = new integer[nnz];
			dbgmres_smoother[ilevel].row_ptr = new integer[n + 1];
#pragma omp parallel for
			for (integer i_1 = 0; i_1 <= n; ++i_1) dbgmres_smoother[ilevel].row_ptr[i_1] = nnz;

			for (integer ii = startpos; ii <= endpos; ++ii) {
				integer istr = ii - iadd;
				integer is1 = row_ptr_start[ii];// +1;
				integer is2 = row_ptr_end[ii];


				// Распараллеливание почемуто тормозит очень сильно.
				//#pragma omp parallel for reduction(+:rsum)
				for (integer ii1 = is1; ii1 <= is2; ++ii1)
				{
					///printf("is1=%d is2=%d\n",is1,is2);
					if (iadd > 20) {
						//printf("i==%d j==%d\n", Amat.i[ii1], Amat.j[ii1]); system("pause");

						//printf("%d \n", ii1 - row_ptr_start[startpos]);
						//system("pause");
					}
					dbgmres_smoother[ilevel].col_ind[ii1 - row_ptr_start[startpos]] = Amat.j[ii1] - 1;
					dbgmres_smoother[ilevel].val[ii1 - row_ptr_start[startpos]] = Amat.aij[ii1];
					if (Amat.j[ii1] != istr/*Amat.i[ii1]*/) {
						dbgmres_smoother[ilevel].val[ii1 - row_ptr_start[startpos]] *= 1.0;
					}
					else {
						dbgmres_smoother[ilevel].val[ii1 - row_ptr_start[startpos]] = 1.0 / dbgmres_smoother[ilevel].val[ii1 - row_ptr_start[startpos]];
					}
					dbgmres_smoother[ilevel].row_ptr[istr/*Amat.i[ii1]*/ - 1] = (dbgmres_smoother[ilevel].row_ptr[istr/*Amat.i[ii1]*/ - 1] < ii1 - row_ptr_start[startpos] ? dbgmres_smoother[ilevel].row_ptr[istr/*Amat.i[ii1]*/ - 1] : ii1 - row_ptr_start[startpos]);
				}
			}
		}

		/*
		if (0) {
			// Матрица col_ind, row_ptr, val собрана верно!!! 14 сентября 2017.

			doublerealT omega = 1.0; // initialize.

			// За подробностями смотри книгу Патрика Роуча стр. 183.
			doublerealT rn = static_cast<doublerealT>(iendq - istartq + 1);
			optimal_omega(rn, omega);

			for (integer i_1 = 0; i_1 < n; ++i_1) {
				doublerealT r = dV[i_1];
				doublerealT ap = 0.0;
				for (integer j_1 = row_ptr[i_1]; j_1 <= row_ptr[i_1 + 1] - 1; ++j_1) {
					if (i_1 != col_ind[j_1]) {
						r += -val[j_1] * dX0[col_ind[j_1]];
						//printf("%e ",);
					}
					else {
						ap = val[j_1];
					}
				}
				dX0[i_1] = (1.0 - omega)*dX0[i_1] + ((omega)*(r)) / ap;
			}
		}
		*/

		// Методы GMRES и BiCGStab расходятся в качестве сглаживателей и 
		// не обеспечивают не то что сглаживания а даже дают расходимость.
		// GMRES Саад и Шульц. [1986]
		integer maxit = 20;
		integer m_restart = my_amg_manager.m_restart;// 20; // Так рекомендовано в статье.

		if (my_amg_manager.b_ChebyshevSmoother) {
			// Выяснено что заходит сюда.
			//std::cout << "basic_func\n";
			//getchar();

			//Seidel(n, dbgmres_smoother[ilevel].val, dbgmres_smoother[ilevel].col_ind, dbgmres_smoother[ilevel].row_ptr, dX0, dV); // работает.
			// BiCGStab просто не работает даже при 40 итерациях.
			//maxit = 40;// Только для BiCGStab
			//Bi_CGStabCRS_smoother(n, dbgmres_smoother[ilevel].val, dbgmres_smoother[ilevel].col_ind, dbgmres_smoother[ilevel].row_ptr, dV, dX0, maxit);
			//chebyshev(dbgmres_smoother[ilevel].val, dbgmres_smoother[ilevel].col_ind, dbgmres_smoother[ilevel].row_ptr, dV, dX0, n);
			chebyshev(dbgmres_smoother[ilevel].val, dbgmres_smoother[ilevel].col_ind, dbgmres_smoother[ilevel].row_ptr, dV, dX0, n, iVar, ilevel);

			//P_Rouch(dbgmres_smoother[ilevel].val, dbgmres_smoother[ilevel].col_ind, dbgmres_smoother[ilevel].row_ptr, dV, dX0, n);
		}
		else {

			gmres(n, dbgmres_smoother[ilevel].val, dbgmres_smoother[ilevel].col_ind, dbgmres_smoother[ilevel].row_ptr, dV, dX0, maxit, m_restart);
		}


		// Возвращение результата расчёта.
#pragma omp parallel for
		for (integer ii = startpos; ii <= endpos; ++ii) {
			integer istr = ii - iadd;
			x[istr] = dX0[ii - startpos];
		}

		delete[] dX0;
		delete[] dV;
		//delete[] val;
		//delete[] col_ind;
		//delete[] row_ptr;
	}
}


// smoother.
  // 9 september 2015 and 4 june 2017.
  // q - quick.
template <typename doublerealT>
void seidelq(Ak2& Amat, integer istartq, integer iendq,
	doublerealT*& x, doublerealT*& b,
	integer*& row_ptr_start, integer*& row_ptr_end,
	integer iadd, bool*& F_false_C_true, integer idirect, integer iVar, integer ilevel)
{

	if (my_amg_manager.b_spai0) {
		for (integer i = istartq + iadd; i <= iendq + iadd; ++i) {
			integer istr = i - iadd;
			if (b[istr] != b[istr]) {
				std::cout << "ilevel= " << ilevel << " istr= " << istr;
				system("pause");
			}
		}

		spai0_smoother<doublerealT>(Amat, istartq, iendq, x, b, row_ptr_start, row_ptr_end, iadd);
	}
	else if (my_amg_manager.b_gmres) {
		if (my_amg_manager.bCFJacoby) {
			seidelqsor2<doublerealT>(Amat, istartq, iendq, x, b, row_ptr_start, row_ptr_end, iadd, F_false_C_true, idirect);
		}
		else {
			seidelqsor2<doublerealT>(Amat, istartq, iendq, x, b, row_ptr_start, row_ptr_end, iadd);
		}
		// gmres smoother или bicgstab smmother если my_amg_manager.b_gmres&&my_amg_manager.b_ChebyshevSmoother
		gmres_smoother<doublerealT>(Amat, istartq, iendq, x, b, row_ptr_start, row_ptr_end, iadd, iVar, ilevel);
	}
	else {
		// лучший выбор
		if (my_amg_manager.bCFJacoby) {
			if ((my_amg_manager.iRunge_Kutta_smoother == 3) || (my_amg_manager.iRunge_Kutta_smoother == 5)) {
				// Трёхшаговый метод Рунге-Кутты.
				integer iorder = my_amg_manager.iRunge_Kutta_smoother;
				Runge_Kutt_3or5<doublerealT>(Amat, istartq, iendq, x, b, row_ptr_start, row_ptr_end, iadd, iorder, F_false_C_true, idirect);
			}
			else {

				seidelqsor2<doublerealT>(Amat, istartq, iendq, x, b, row_ptr_start, row_ptr_end, iadd, F_false_C_true, idirect);
			}
		}
		else {
			if ((my_amg_manager.iRunge_Kutta_smoother == 3) || (my_amg_manager.iRunge_Kutta_smoother == 5)) {
				// Трёхшаговый метод Рунге-Кутты.
				integer iorder = my_amg_manager.iRunge_Kutta_smoother;
				Runge_Kutt_3or5<doublerealT>(Amat, istartq, iendq, x, b, row_ptr_start, row_ptr_end, iadd, iorder);
			}
			else {

				seidelqsor2<doublerealT>(Amat, istartq, iendq, x, b, row_ptr_start, row_ptr_end, iadd);
			}

		}
	}

} // seidelq


// smoother.
  // 9 september 2015 and 4 june 2017.
  // q - quick.
template <typename doublerealT>
void seidelq(Ak2& Amat, integer istartq, integer iendq,
	doublerealT*& x, doublerealT*& b,
	integer*& row_ptr_start, integer*& row_ptr_end,
	integer iadd, bool*& F_false_C_true, integer idirect,
	doublerealT*& diag_minus_one, integer iVar, integer ilevel, integer ibsp_length)
{
	if (my_amg_manager.b_spai0) {

		for (integer i = istartq + iadd; i <= iendq + iadd; ++i) {
			integer istr = i - iadd;
			if (b[istr] != b[istr]) {
				std::cout << "ibsp_length ilevel= " << ilevel << " istr= " << istr;
				system("pause");
			}
		}

		spai0_smoother<doublerealT>(Amat, istartq, iendq, x, b, row_ptr_start, row_ptr_end, iadd, diag_minus_one);
	}
	else if (my_amg_manager.b_gmres) {
		if (my_amg_manager.bCFJacoby) {
			seidelqsor2<doublerealT>(Amat, istartq, iendq, x, b, row_ptr_start, row_ptr_end, iadd, F_false_C_true, idirect, diag_minus_one, ibsp_length);
		}
		else {
			seidelqsor2<doublerealT>(Amat, istartq, iendq, x, b, row_ptr_start, row_ptr_end, iadd, diag_minus_one);
		}
		// gmres smoother или bicgstab smmother если my_amg_manager.b_gmres&&my_amg_manager.b_ChebyshevSmoother
		gmres_smoother<doublerealT>(Amat, istartq, iendq, x, b, row_ptr_start, row_ptr_end, iadd, iVar, ilevel);
	}
	else {
		// лучший выбор
		if (my_amg_manager.bCFJacoby) {
			if ((my_amg_manager.iRunge_Kutta_smoother == 3) || (my_amg_manager.iRunge_Kutta_smoother == 5)) {
				// Трёхшаговый метод Рунге-Кутты.
				integer iorder = my_amg_manager.iRunge_Kutta_smoother;
				Runge_Kutt_3or5<doublerealT>(Amat, istartq, iendq, x, b, row_ptr_start, row_ptr_end, iadd, iorder, F_false_C_true, idirect, diag_minus_one, ibsp_length);
			}
			else {

				seidelqsor2<doublerealT>(Amat, istartq, iendq, x, b, row_ptr_start, row_ptr_end, iadd, F_false_C_true, idirect, diag_minus_one, ibsp_length);
			}
		}
		else {
			if ((my_amg_manager.iRunge_Kutta_smoother == 3) || (my_amg_manager.iRunge_Kutta_smoother == 5)) {
				// Трёхшаговый метод Рунге-Кутты.
				integer iorder = my_amg_manager.iRunge_Kutta_smoother;
				Runge_Kutt_3or5<doublerealT>(Amat, istartq, iendq, x, b, row_ptr_start, row_ptr_end, iadd, iorder, diag_minus_one, ibsp_length);
			}
			else {

				seidelqsor2<doublerealT>(Amat, istartq, iendq, x, b, row_ptr_start, row_ptr_end, iadd, diag_minus_one);
			}

		}
	}

} // seidelq

 


// 08.01.2018 Перенесено в отдельную функцию, т.к. используется неоднократно. Передаётся большое число параметров.
//  A*z76=s76; 891 строка.
template <typename doublerealT>
void V_cycle_solve(Ak2 &Amat, doublereal* &z76, doublereal* &s76, bool process_flow_logic, integer* &row_ptr_start,
	integer* &row_ptr_end, doublerealT** &residual_fine, doublerealT** &diag, integer* n_a, bool bonly_serial,
	doublerealT process_flow_beta, bool* &F_false_C_true, integer &nu1, integer &nu2, integer bILU2smoother,
	integer ilevel, integer inumberVcyclelocbicgstab, INIT_SELECTOR_CASE_CAMG_RUMBAv_0_14 imyinit, const integer idim_diag,
	LEVEL_ADDITIONAL_DATA* &milu2, LEVEL_ADDITIONAL_DATA0* milu0, bool** &nested_desection, 
	Ak1* &P, // prolongation он же restriction (метод не чувствителен к сортировке).
	integer* nnz_aRP, doublerealT** &residual_coarse, integer igam, integer* nnz_a,
	doublerealT** &error_approx_coarse, doublerealT dapply_ilu_max_pattern_size,
	doublerealT process_flow_alpha, doublerealT** &error_approx_fine,
	integer nFinestSweeps, integer ibsp_length, integer iVar) {

	// Один V - цикл алгебраического многосеточного метода.
	// A*z76=s76;

	// надо увеличивать nu1, nu2 с 1,2 до 5 наверно.

	for (integer i_1 = 0; i_1 < inumberVcyclelocbicgstab; ++i_1)
	{


		doublerealT R0_0 = 0.0;
		doublerealT Rprev_0 = 0.0, Rnext_0 = 0.0;
		if (process_flow_logic) {
			// calculate initial residual.
			//residualq(Amat, 1, n_a[0], z76, s76, row_ptr_start, row_ptr_end, 0, residual_fine[0]);
			residualq2(Amat, static_cast<integer>(1), n_a[0], z76, s76, row_ptr_start, row_ptr_end, static_cast<integer>(0), residual_fine[0], diag[0]);
			R0_0 = norma(residual_fine[0], n_a[0]);
			Rprev_0 = R0_0;

			// smother
			integer iter = 0;
			for (iter = 0; iter < nu1; ++iter) {
				//quick seidel
				if (bonly_serial) {
					seidelq<doublereal>(Amat, static_cast<integer>(1), n_a[0], z76, s76, row_ptr_start, row_ptr_end, 0, F_false_C_true, 0, iVar, 0);
				}
				else {
					seidelq<doublereal>(Amat, static_cast<integer>(1), n_a[0], z76, s76, nested_desection[0], row_ptr_start, row_ptr_end, 0, F_false_C_true, 0);
				}
				//residualq(Amat, 1, n_a[0], z76, s76, row_ptr_start, row_ptr_end, 0, residual_fine[0]);
				residualq2(Amat, static_cast<integer>(1), n_a[0], z76, s76, row_ptr_start, row_ptr_end, 0, residual_fine[0], diag[0]);
				Rnext_0 = norma(residual_fine[0], n_a[0]);
				// this is process flow logic
				if (Rnext_0 > process_flow_beta*Rprev_0) {
					// Смысл модификации в том что мы экономим итерации на пресмутере.
					break; // досрочно опускаемся на следующий уровень если он есть конечно.
				}
				else {
					Rprev_0 = Rnext_0;
				}
			}
			if (iter == nu1) {
				printf("level 0 limit presmother iteration is reached\n");
			}

		}
		else {
			// smoother
			for (integer iter = 0; iter < nu1; ++iter) {
				//seidel(Amat, 1, nnz_a[0], z76, s76, flag, n_a[0]);
				//quick seidel
				if (bonly_serial) {
					if (bILU2smoother == 1) {
						// ILU0
						seidelq<doublereal>(Amat, static_cast<integer>(1), n_a[0], z76, s76, row_ptr_start, row_ptr_end, 0, F_false_C_true, 0, iVar, 0);
						residualq2(Amat, static_cast<integer>(1), n_a[0], z76, s76, row_ptr_start, row_ptr_end, 0, residual_fine[0], diag[0]);
#pragma omp parallel for
						for (integer i = 1; i <= n_a[0]; ++i) {
							milu0[0].zbuf[i] = residual_fine[0][i];
						}
						lusol_1patchforRUMBA(n_a[0], milu0[0].zbuf, milu0[0].zbuf2, milu0[0]);
#pragma omp parallel for
						for (integer i = 1; i <= n_a[0]; ++i) {
							z76[i] += milu0[0].zbuf2[i];
						}
					}
					else if ((bILU2smoother == 2) || (my_amg_manager.iFinnest_ilu == 1)) {
						seidelq<doublereal>(Amat, static_cast<integer>(1), n_a[0], z76, s76, row_ptr_start, row_ptr_end, 0, F_false_C_true, 0, iVar, 0);
						residualq2(Amat, static_cast<integer>(1), n_a[0], z76, s76, row_ptr_start, row_ptr_end, 0, residual_fine[0], diag[0]);
#pragma omp parallel for
						for (integer i = 1; i <= n_a[0]; ++i) {
							milu2[0].zbuf[i] = residual_fine[0][i];
						}
						lusol_1patchforRUMBA(n_a[0], milu2[0].zbuf, milu2[0].zbuf2, milu2[0]);
#pragma omp parallel for
						for (integer i = 1; i <= n_a[0]; ++i) {
							z76[i] += milu2[0].zbuf2[i];
						}
					}
					else {
						seidelq<doublereal>(Amat, static_cast<integer>(1), n_a[0], z76, s76, row_ptr_start, row_ptr_end, 0, F_false_C_true, 0, iVar, 0);
					}
				}
				else {
					seidelq<doublereal>(Amat, static_cast<integer>(1), n_a[0], z76, s76, nested_desection[0], row_ptr_start, row_ptr_end, 0, F_false_C_true, 0);
				}
			}
		}

		//exporttecplot(x, n);

		move_down(nu1, nu2);

		if (!process_flow_logic) {
			// residual_r
			//doublerealT *residual_fine[0] = new doublerealT[n_a[0] + 1];
			//residual<doublereal>(Amat, static_cast<integer>(1), nnz_a[0], z76, s76, flag, n_a[0], residual_fine[0]);
			//residualq<doublereal>(Amat, static_cast<integer>(1), n_a[0], z76, s76, row_ptr_start, row_ptr_end, 0, residual_fine[0]);
			residualq2(Amat, static_cast<integer>(1), n_a[0], z76, s76, row_ptr_start, row_ptr_end, 0, residual_fine[0], diag[0]);
		}



		//if (iprohod%5==0) system("pause");
		if (ilevel > 1) {

			//doublerealT *residual_coarse = new doublerealT[n_a[1] + 1];
#pragma omp parallel for
			for (integer ii = 1; ii <= n_a[1]; ++ii) {
				residual_coarse[0][ii] = 0.0;
			}

			// restriction
			restriction(P, static_cast<integer>(1), nnz_aRP[0],  residual_fine[0], residual_coarse[0],  n_a[1]);


			// Amat*e=r;
			//doublerealT* error_approx_coarse = new doublerealT[n_a[1] + 1];
			if ((imyinit == INIT_SELECTOR_CASE_CAMG_RUMBAv_0_14::ZERO_INIT)) {
#pragma omp parallel for
				for (integer ii = 1; ii <= n_a[1]; ++ii) {
					error_approx_coarse[0][ii] = 0.0;
				}
			}

			if ((imyinit == INIT_SELECTOR_CASE_CAMG_RUMBAv_0_14::RANDOM_INIT)) {
				// (1,110); (0.8, 37); (0.7, 29); (0.6, 25); (0.5, 20); (0.4, 17); (0.3, 18); (0.0, 19);
#pragma omp parallel for
				for (integer ii = 1; ii <= n_a[1]; ++ii) {
					error_approx_coarse[0][ii] = 0.4*fabs(residual_coarse[0][ii])*(rand() / (static_cast<doublerealT>(RAND_MAX)));
				}
			}


			for (integer i_2 = 1; i_2 <= igam; ++i_2) {
				doublerealT R0_1 = 0.0;
				doublerealT Rprev_1 = 0.0, Rnext_1 = 0.0;
				if (process_flow_logic) {
					// calculate initial residual.
					//residualq(Amat, static_cast<integer>(1), n_a[1], error_approx_coarse[0], residual_coarse[0], row_ptr_start, row_ptr_end, n_a[0], residual_fine[1]);
					residualq2(Amat, static_cast<integer>(1), n_a[1], error_approx_coarse[0], residual_coarse[0], row_ptr_start, row_ptr_end, n_a[0], residual_fine[1], diag[1]);

					R0_1 = norma(residual_fine[1], n_a[1]);
					Rprev_1 = R0_1;

					// smother
					integer iter = 0;
					for (iter = 0; iter < nu1; ++iter) {
						//quick seidel
						if (bonly_serial) {
							seidelq(Amat, static_cast<integer>(1), n_a[1], error_approx_coarse[0], residual_coarse[0], row_ptr_start, row_ptr_end, n_a[0], F_false_C_true, 0,1);
						}
						else {
							seidelq(Amat, static_cast<integer>(1), n_a[1], error_approx_coarse[0], residual_coarse[0], nested_desection[1], row_ptr_start, row_ptr_end, n_a[0], F_false_C_true, 0);
						}

						//residualq(Amat, static_cast<integer>(1), n_a[1], error_approx_coarse[0], residual_coarse[0], row_ptr_start, row_ptr_end, n_a[0], residual_fine[1]);
						residualq2(Amat, static_cast<integer>(1), n_a[1], error_approx_coarse[0], residual_coarse[0], row_ptr_start, row_ptr_end, n_a[0], residual_fine[1], diag[1]);

						Rnext_1 = norma(residual_fine[1], n_a[1]);
						// this is process flow logic
						if (Rnext_1 > process_flow_beta*Rprev_1) {
							// Смысл модификации в том что мы экономим итерации на пресмутере.
							break; // досрочно опускаемся на следующий уровень если он есть конечно.
						}
						else {
							Rprev_1 = Rnext_1;
						}
					}

					if (iter == nu1) {
						printf("level 1 limit presmother iteration is reached\n");
					}

				}
				else {

					// pre smothing
					for (integer iter = 0; iter < nu1; ++iter) {
						//seidel(Amat, 1 + 2 * nnz_a[0], 2 * nnz_a[0] + nnz_a[1], error_approx_coarse, residual_coarse, flag, n_a[1]);
						if (bonly_serial) {
							if (bILU2smoother == 1) {
								seidelq(Amat, static_cast<integer>(1), n_a[1], error_approx_coarse[0], residual_coarse[0], row_ptr_start, row_ptr_end, n_a[0], F_false_C_true, 0, static_cast<integer>(1));
								residualq2(Amat, static_cast<integer>(1), n_a[1], error_approx_coarse[0], residual_coarse[0], row_ptr_start, row_ptr_end, n_a[0], residual_fine[1], diag[1]);
#pragma omp parallel for
								for (integer i = 1; i <= n_a[1]; ++i) {
									milu0[1].zbuf[i] = residual_fine[1][i];
								}
								lusol_1patchforRUMBA(n_a[1], milu0[1].zbuf, milu0[1].zbuf2, milu0[1]);
#pragma omp parallel for
								for (integer i = 1; i <= n_a[1]; ++i) {
									error_approx_coarse[0][i + 1] += milu0[1].zbuf2[i + 1];
								}


							}
							else if (1 && bILU2smoother == 2) {
								seidelq(Amat, static_cast<integer>(1), n_a[1], error_approx_coarse[0], residual_coarse[0], row_ptr_start, row_ptr_end, n_a[0], F_false_C_true, 0, static_cast<integer>(1));
								residualq2(Amat, static_cast<integer>(1), n_a[1], error_approx_coarse[0], residual_coarse[0], row_ptr_start, row_ptr_end, n_a[0], residual_fine[1], diag[1]);
#pragma omp parallel for
								for (integer i = 1; i <= n_a[1]; ++i) {
									milu2[1].zbuf[i] = residual_fine[1][i];
								}
								lusol_1patchforRUMBA(n_a[1], milu2[1].zbuf, milu2[1].zbuf2, milu2[1]);
#pragma omp parallel for
								for (integer i = 1; i <= n_a[1]; ++i) {
									error_approx_coarse[0][i] += milu2[1].zbuf2[i];
								}


							}
							else {
								seidelq(Amat, static_cast<integer>(1), n_a[1], error_approx_coarse[0], residual_coarse[0], row_ptr_start, row_ptr_end, n_a[0], F_false_C_true, 0, static_cast<integer>(1));
							}
						}
						else {
							seidelq(Amat, static_cast<integer>(1), n_a[1], error_approx_coarse[0], residual_coarse[0], nested_desection[1], row_ptr_start, row_ptr_end, n_a[0], F_false_C_true, 0);
						}
					}
				}





				doublerealT* R0_21 = new doublerealT[idim_diag];
				doublerealT* Rprev_21 = new doublerealT[idim_diag];
				doublerealT* Rnext_21 = new doublerealT[idim_diag];

				//3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21
				for (integer i_id_level_local = 2; i_id_level_local < idim_diag; ++i_id_level_local) {

					move_down(nu1, nu2);


					if (ilevel > i_id_level_local) {
						// Впервые создан универсальный програмный код на произвольное количество уровней вложенности 
						// в solution phase.
						// Нужно только написать универсальный обработчик для R0_21.
						// 4 декабря 2016. 

						// residual
						if (!process_flow_logic) {
							integer in_a_loc = 0;
							for (integer i_72 = 0; i_72 < i_id_level_local - 1; ++i_72) {
								// sum(n_a[0]+n_a[19];
								in_a_loc += n_a[i_72];
							}
							//doublerealT *residual_fine[i_id_level_local-1] = new doublerealT[n_a[i_id_level_local-1] + 1];
							//residual(Amat, 1 + nnz_a[0] + nnz_a[1] + nnz_a[2] + nnz_a[3] + nnz_a[4] + nnz_a[5] + nnz_a[6]+ nnz_a[7]+ nnz_a[8]+ nnz_a[9]+ nnz_a[10]+ nnz_a[11]+ nnz_a[12]+ nnz_a[13]+ nnz_a[14]+ nnz_a[15]+ nnz_a[16]+nnz_a[17]+nnz_a[18]+nnz_a[19], nnz_a[0] + nnz_a[1] + nnz_a[2] + nnz_a[3] + nnz_a[4] + nnz_a[5] + nnz_a[6] + nnz_a[7]+ nnz_a[8]+ nnz_a[9]+ nnz_a[10]+ nnz_a[11]+ nnz_a[12]+ nnz_a[13]+ nnz_a[14]+ nnz_a[15]+ nnz_a[16]+ nnz_a[17]+nnz_a[18]+nnz_a[19]+nnz_a[20], error_approx_coarse[i_id_level_local - 2], residual_coarse[i_id_level_local - 2], flag, n_a[i_id_level_local - 1], residual_fine[i_id_level_local - 1]);
							//residualq(Amat, static_cast<integer>(1), n_a[i_id_level_local-1], error_approx_coarse[i_id_level_local - 2], residual_coarse[i_id_level_local - 2], row_ptr_start, row_ptr_end,  in_a_loc , residual_fine[i_id_level_local - 1]);
							//residualq2(Amat, static_cast<integer>(1), n_a[i_id_level_local-1], error_approx_coarse[i_id_level_local - 2], residual_coarse[i_id_level_local - 2], row_ptr_start, row_ptr_end, in_a_loc, residual_fine[i_id_level_local - 1], diag[i_id_level_local - 1]);
							residualq2(Amat, static_cast<integer>(1), n_a[i_id_level_local - 1], error_approx_coarse[i_id_level_local - 2], residual_coarse[i_id_level_local - 2], row_ptr_start, row_ptr_end, in_a_loc, residual_fine[i_id_level_local - 1], diag[i_id_level_local - 1]);

						}

						//doublerealT *residual_coarse[i_id_level_local-1] = new doublerealT[n_a[i_id_level_local] + 1];

						// restriction

						integer innz_aRP_loc = 1;
						for (integer i_72 = 0; i_72 < i_id_level_local - 1; ++i_72) {
							innz_aRP_loc += nnz_aRP[i_72];
						}
						//restriction(P, 1 + nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2] + nnz_aRP[3] + nnz_aRP[4] + nnz_aRP[5] + nnz_aRP[6] + nnz_aRP[7] + nnz_aRP[8] + nnz_aRP[9] + nnz_aRP[10] + nnz_aRP[11] + nnz_aRP[12] + nnz_aRP[13] + nnz_aRP[14] + nnz_aRP[15] + nnz_aRP[16] + nnz_aRP[17] + nnz_aRP[18] + nnz_aRP[19], nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2] + nnz_aRP[3] + nnz_aRP[4] + nnz_aRP[5] + nnz_aRP[6] + nnz_aRP[7] + nnz_aRP[8] + nnz_aRP[9] + nnz_aRP[10] + nnz_aRP[11] + nnz_aRP[12] + nnz_aRP[13] + nnz_aRP[14] + nnz_aRP[15] + nnz_aRP[16] + nnz_aRP[17] + nnz_aRP[18] + nnz_aRP[19] + nnz_aRP[20], flag, residual_fine[i_id_level_local - 1], residual_coarse[i_id_level_local - 1], n_a[i_id_level_local - 1], n_a[i_id_level_local]);
						restriction(P, innz_aRP_loc, innz_aRP_loc - 1 + nnz_aRP[i_id_level_local - 1],  residual_fine[i_id_level_local - 1], residual_coarse[i_id_level_local - 1], n_a[i_id_level_local]);



						// Amat*e=r;
						//doublerealT* error_approx_coarse[i_id_level_local-1] = new doublerealT[n_a[i_id_level_local] + 1];

						if (0) {
							//if (icount_V_cycle == 1) {
							//	for (integer ii = 1; ii <= n_a[i_id_level_local]; ++ii) {
							// обнуление только на первом проходе.
							//		error_approx_coarse[i_id_level_local - 1][ii] = 0.0;
							//error_approx_coarse[i_id_level_local - 1][ii] = (rand() / (static_cast<doublerealT>(RAND_MAX)));
							//	}
							//}
							// На задаче BSK_Dmitrii random оптимизация оказалась гораздо хуже, но сходимость тоже присутствовала.
#pragma omp parallel for
							for (integer ii = 1; ii <= n_a[i_id_level_local]; ++ii) {
								// 0.4*fabs(residual_coarse[i_id_level_local - 1][ii]) - амплитуда.
								// 0.4 - демпфирующий множитель.
								// diag[i_id_level_local]
								// (1,110); (0.8, 37); (0.7, 29); (0.6, 25); (0.5, 20); (0.4, 17); (0.3, 18); (0.0, 19);
								error_approx_coarse[i_id_level_local - 1][ii] = 0.4*fabs(residual_coarse[i_id_level_local - 1][ii])*(rand() / (static_cast<doublerealT>(RAND_MAX)));
							}
						}
						else {
							if (imyinit == INIT_SELECTOR_CASE_CAMG_RUMBAv_0_14::ZERO_INIT) {
#pragma omp parallel for
								for (integer ii = 1; ii <= n_a[i_id_level_local]; ++ii) {
									error_approx_coarse[i_id_level_local - 1][ii] = 0.0;
								}
							}
							if (imyinit == INIT_SELECTOR_CASE_CAMG_RUMBAv_0_14::RANDOM_INIT) {
#pragma omp parallel for
								for (integer ii = 1; ii <= n_a[i_id_level_local]; ++ii) {
									// (1,110); (0.8, 37); (0.7, 29); (0.6, 25); (0.5, 20); (0.4, 17); (0.3, 18); (0.0, 19);
									error_approx_coarse[i_id_level_local - 1][ii] = 0.4*fabs(residual_coarse[i_id_level_local - 1][ii])*(rand() / (static_cast<doublerealT>(RAND_MAX)));
								}
							}
						}


						//for (integer i_2 = 1; i_2 <= igam; ++i_2)
						{
							// pre smothing
							//doublerealT R0_21 = 0.0;
							//doublerealT Rprev_21 = 0.0, Rnext_21 = 0.0;
							R0_21[i_id_level_local] = 0.0;
							Rprev_21[i_id_level_local] = 0.0;
							Rnext_21[i_id_level_local] = 0.0;
							if (process_flow_logic) {

								integer in_a_loc = 0;
								for (integer i_72 = 0; i_72 < i_id_level_local - 1; ++i_72) {
									// sum(n_a[0]+n_a[19];
									in_a_loc += n_a[i_72];
								}

								// calculate initial residual.
								//residualq(Amat, static_cast<integer>(1), n_a[i_id_level_local], error_approx_coarse[i_id_level_local-1], residual_coarse[i_id_level_local-1], row_ptr_start, row_ptr_end,in_a_loc + n_a[i_id_level_local-1] , residual_fine[i_id_level_local]);
								residualq2(Amat, static_cast<integer>(1), n_a[i_id_level_local], error_approx_coarse[i_id_level_local - 1], residual_coarse[i_id_level_local - 1], row_ptr_start, row_ptr_end, in_a_loc + n_a[i_id_level_local - 1], residual_fine[i_id_level_local], diag[i_id_level_local]);



								R0_21[i_id_level_local] = norma(residual_fine[i_id_level_local], n_a[i_id_level_local]);
								Rprev_21[i_id_level_local] = R0_21[i_id_level_local];



								// smother
								integer iter = 0;
								integer nu1_count = nu1;
								if (i_id_level_local == ilevel - 1) {
									// на самом грубом уровне нам требуется точное решение.
									//nu1_count = 100;

									printf("Direct method is start.../n");

									IMatrix sparseS; // разреженная матрица в формате IMatrix
									initIMatrix(&sparseS, n_a[i_id_level_local]);

									const integer jstart = 1 + in_a_loc + n_a[i_id_level_local - 1];
									const integer jend = n_a[i_id_level_local] + in_a_loc + n_a[i_id_level_local - 1];

									for (integer j = jstart; j <= jend; ++j) {

										//setValueIMatrix(&sparseS, Amat.i[row_ptr_start[j]] - 1, Amat.i[row_ptr_start[j]] - 1, 1.0 / Amat.aij[row_ptr_start[j]]);
										setValueIMatrix(&sparseS, j - jstart, j - jstart, 1.0 / Amat.aij[row_ptr_start[j]]);//3.02.2019
										//if (Amat.i[row_ptr_start[j]] != j - jstart + 1) {
											//printf("Amat.i[row_ptr_start[j]]=%lld j - jstart + 1=%lld\n", Amat.i[row_ptr_start[j]],j - jstart + 1);
											//system("pause");
										//}
										const doublerealT nonzeroEPS = 1e-37; // для отделения вещественного нуля

										const integer istart = row_ptr_start[j] + 1;
										const integer iend = row_ptr_end[j];

										for (integer i = istart; i <= iend; ++i)
										{

											if ((fabs(Amat.aij[i]) > nonzeroEPS)) {
												//setValueIMatrix(&sparseS, Amat.i[i] - 1, Amat.j[i] - 1, Amat.aij[i]);
												setValueIMatrix(&sparseS, i - istart, Amat.j[i] - 1, Amat.aij[i]);
											}

										}
									}

									doublereal* dX025 = new doublereal[n_a[i_id_level_local]];
									doublereal* dV25 = new doublereal[n_a[i_id_level_local]];
#pragma omp parallel for
									for (integer i = 0; i < n_a[i_id_level_local]; ++i) {
										dX025[i] = error_approx_coarse[i_id_level_local - 1][i + 1];
										dV25[i] = residual_coarse[i_id_level_local - 1][i + 1];
									}

									// главный метод, возвращающий решение x,
									// принимает вектор свободных членов b и 
									// квадратную матрицу xO в специальном разреженном формате.
									// реализация без барьера и итерационного уточнения.
									calculateSPARSEgaussArray(&sparseS, dX025, dV25);
#pragma omp parallel for
									for (integer i = 0; i < n_a[i_id_level_local]; ++i) {
										error_approx_coarse[i_id_level_local - 1][i + 1] = dX025[i];
									}

									delete[] dX025;
									delete[] dV25;


									freeIMatrix(&sparseS);

									//residualq(Amat, static_cast<integer>(1), n_a[i_id_level_local], error_approx_coarse[i_id_level_local-1], residual_coarse[i_id_level_local-1], row_ptr_start, row_ptr_end, in_a_loc + n_a[i_id_level_local - 1], residual_fine[i_id_level_local]);
									residualq2(Amat, static_cast<integer>(1), n_a[i_id_level_local], error_approx_coarse[i_id_level_local - 1], residual_coarse[i_id_level_local - 1], row_ptr_start, row_ptr_end, in_a_loc + n_a[i_id_level_local - 1], residual_fine[i_id_level_local], diag[i_id_level_local]);


									Rnext_21[i_id_level_local] = norma(residual_fine[i_id_level_local], n_a[i_id_level_local]);
									// this is process flow logic
									if (Rnext_21[i_id_level_local] > process_flow_beta*Rprev_21[i_id_level_local]) {
										// Смысл модификации в том что мы экономим итерации на пресмутере.
										break; // досрочно опускаемся на следующий уровень если он есть конечно.
									}
									else {
										Rprev_21[i_id_level_local] = Rnext_21[i_id_level_local];
									}
								}
								else {
									for (iter = 0; iter < nu1_count; ++iter) {
										//quick seidel
										if (bonly_serial) {
											seidelq(Amat, static_cast<integer>(1), n_a[i_id_level_local], error_approx_coarse[i_id_level_local - 1], residual_coarse[i_id_level_local - 1], row_ptr_start, row_ptr_end, in_a_loc + n_a[i_id_level_local - 1], F_false_C_true, 0, i_id_level_local);
										}
										else {
											seidelq(Amat, static_cast<integer>(1), n_a[i_id_level_local], error_approx_coarse[i_id_level_local - 1], residual_coarse[i_id_level_local - 1], nested_desection[i_id_level_local], row_ptr_start, row_ptr_end, in_a_loc + n_a[i_id_level_local - 1], F_false_C_true, 0);
										}

										//residualq(Amat, static_cast<integer>(1), n_a[i_id_level_local], error_approx_coarse[i_id_level_local-1], residual_coarse[i_id_level_local-1], row_ptr_start, row_ptr_end, in_a_loc + n_a[i_id_level_local - 1], residual_fine[i_id_level_local]);
										residualq2(Amat, static_cast<integer>(1), n_a[i_id_level_local], error_approx_coarse[i_id_level_local - 1], residual_coarse[i_id_level_local - 1], row_ptr_start, row_ptr_end, in_a_loc + n_a[i_id_level_local - 1], residual_fine[i_id_level_local], diag[i_id_level_local]);


										Rnext_21[i_id_level_local] = norma(residual_fine[i_id_level_local], n_a[i_id_level_local]);
										// this is process flow logic
										if (Rnext_21[i_id_level_local] > process_flow_beta*Rprev_21[i_id_level_local]) {
											// Смысл модификации в том что мы экономим итерации на пресмутере.
											break; // досрочно опускаемся на следующий уровень если он есть конечно.
										}
										else {
											Rprev_21[i_id_level_local] = Rnext_21[i_id_level_local];
										}
									}
								}

								if (iter == nu1) {
#if doubleintprecision == 1
									printf("level %lld limit presmother iteration is reached\n", i_id_level_local);
#else
									printf("level %d limit presmother iteration is reached\n", i_id_level_local);
#endif

								}

							}
							else {
								integer nu1_count = nu1;
								if (i_id_level_local == ilevel - 1) {
									// на самом грубом уровне нам требуется точное решение.
									//nu1_count = 100;
								}
								for (integer iter = 0; iter < nu1_count; ++iter) {

									integer in_a_loc = 0;
									for (integer i_72 = 0; i_72 < i_id_level_local - 1; ++i_72) {
										// sum(n_a[0]+n_a[19];
										in_a_loc += n_a[i_72];
									}

									if (bonly_serial) {

										bool bflag56 = false;
										if (my_amg_manager.iFinnest_ilu == 1) {
											if (my_amg_manager.b_ilu_smoothers_in_nnz_n_LE_6) {
												doublerealT dn = 1.0*n_a[i_id_level_local];
												doublerealT dnnz = 1.0*nnz_a[i_id_level_local];
												if (dnnz / dn <= dapply_ilu_max_pattern_size) {
													bflag56 = true;
												}
											}
										}

										if (bILU2smoother == 1) {
											seidelq(Amat, static_cast<integer>(1), n_a[i_id_level_local], error_approx_coarse[i_id_level_local - 1], residual_coarse[i_id_level_local - 1], row_ptr_start, row_ptr_end, in_a_loc + n_a[i_id_level_local - 1], F_false_C_true, 0, i_id_level_local);
											residualq2(Amat, static_cast<integer>(1), n_a[i_id_level_local], error_approx_coarse[i_id_level_local - 1], residual_coarse[i_id_level_local - 1], row_ptr_start, row_ptr_end, in_a_loc + n_a[i_id_level_local - 1], residual_fine[i_id_level_local], diag[i_id_level_local]);
#pragma omp parallel for
											for (integer i = 1; i <= n_a[i_id_level_local]; ++i) {
												milu0[i_id_level_local].zbuf[i] = residual_fine[i_id_level_local][i];
											}
											lusol_1patchforRUMBA(n_a[i_id_level_local], milu0[i_id_level_local].zbuf, milu0[i_id_level_local].zbuf2, milu0[i_id_level_local]);
#pragma omp parallel for
											for (integer i = 1; i <= n_a[i_id_level_local]; ++i) {
												error_approx_coarse[i_id_level_local - 1][i] += milu0[i_id_level_local].zbuf2[i];
											}

										}
										else if (1 && ((bILU2smoother == 2) || bflag56)) {
											seidelq(Amat, static_cast<integer>(1), n_a[i_id_level_local], error_approx_coarse[i_id_level_local - 1], residual_coarse[i_id_level_local - 1], row_ptr_start, row_ptr_end, in_a_loc + n_a[i_id_level_local - 1], F_false_C_true, 0, i_id_level_local);
											residualq2(Amat, static_cast<integer>(1), n_a[i_id_level_local], error_approx_coarse[i_id_level_local - 1], residual_coarse[i_id_level_local - 1], row_ptr_start, row_ptr_end, in_a_loc + n_a[i_id_level_local - 1], residual_fine[i_id_level_local], diag[i_id_level_local]);
#pragma omp parallel for
											for (integer i = 1; i <= n_a[i_id_level_local]; ++i) {
												milu2[i_id_level_local].zbuf[i] = residual_fine[i_id_level_local][i];
											}
											lusol_1patchforRUMBA(n_a[i_id_level_local], milu2[i_id_level_local].zbuf, milu2[i_id_level_local].zbuf2, milu2[i_id_level_local]);
#pragma omp parallel for
											for (integer i = 1; i <= n_a[i_id_level_local]; ++i) {
												error_approx_coarse[i_id_level_local - 1][i] += milu2[i_id_level_local].zbuf2[i];
											}

										}
										else {
											seidelq(Amat, static_cast<integer>(1), n_a[i_id_level_local], error_approx_coarse[i_id_level_local - 1], residual_coarse[i_id_level_local - 1], row_ptr_start, row_ptr_end, in_a_loc + n_a[i_id_level_local - 1], F_false_C_true, 0, i_id_level_local);
										}
									}
									else {
										seidelq(Amat, static_cast<integer>(1), n_a[i_id_level_local], error_approx_coarse[i_id_level_local - 1], residual_coarse[i_id_level_local - 1], nested_desection[i_id_level_local], row_ptr_start, row_ptr_end, in_a_loc + n_a[i_id_level_local - 1], F_false_C_true, 0);
									}
								}
							}
						}

					}

				}







				//center
				// ЭТО сердцевина SOLUTION PHASE.

				//  нижний 3.12.2016 (осталось один нижний и один верхний).


				// 21
				for (integer i_id_level_local = idim_diag - 1; i_id_level_local >= 2; i_id_level_local--) {

					if (ilevel > i_id_level_local) {

						{

							integer in_a_loc = 0;
							for (integer i_72 = 0; i_72 < i_id_level_local - 1; ++i_72) {
								// sum(n_a[0]+n_a[19];
								in_a_loc += n_a[i_72];
							}

							// post smoothing
							// doublerealT R0_20 = 0.0;
							///doublerealT Rprev_20 = 0.0, Rnext_20 = 0.0;
							if (process_flow_logic) {
								// calculate initial residual.
								//residualq(Amat, 1, n_a[i_id_level_local], error_approx_coarse[i_id_level_local - 1], residual_coarse[i_id_level_local - 1], row_ptr_start, row_ptr_end, in_a_loc  + n_a[i_id_level_local - 1], residual_fine[i_id_level_local]);



								//Rprev_21[i_id_level_local] = norma(residual_fine[i_id_level_local], n_a[i_id_level_local]);


								// smother
								integer iter = 0;
								integer nu2_count = nu2;
								if (i_id_level_local == ilevel - 1) {
									// на самом грубом уровне нам требуется точное решение.
									//nu2_count = 100;																												
								}
								else {
									for (iter = 0; iter < nu2_count; ++iter) {
										//quick seidel
										if (bonly_serial) {
											seidelq(Amat, static_cast<integer>(1), n_a[i_id_level_local], error_approx_coarse[i_id_level_local - 1], residual_coarse[i_id_level_local - 1], row_ptr_start, row_ptr_end, in_a_loc + n_a[i_id_level_local - 1], F_false_C_true, static_cast<integer>(1), i_id_level_local);
										}
										else {
											seidelq(Amat, static_cast<integer>(1), n_a[i_id_level_local], error_approx_coarse[i_id_level_local - 1], residual_coarse[i_id_level_local - 1], nested_desection[i_id_level_local], row_ptr_start, row_ptr_end, in_a_loc + n_a[i_id_level_local - 1], F_false_C_true, static_cast<integer>(1));
										}

										//residualq(Amat, 1, n_a[i_id_level_local], error_approx_coarse[i_id_level_local - 1], residual_coarse[i_id_level_local - 1], row_ptr_start, row_ptr_end, in_a_loc + n_a[i_id_level_local-1], residual_fine[i_id_level_local]);
										residualq2(Amat, static_cast<integer>(1), n_a[i_id_level_local], error_approx_coarse[i_id_level_local - 1], residual_coarse[i_id_level_local - 1], row_ptr_start, row_ptr_end, in_a_loc + n_a[i_id_level_local - 1], residual_fine[i_id_level_local], diag[i_id_level_local]);


										Rnext_21[i_id_level_local] = norma(residual_fine[i_id_level_local], n_a[i_id_level_local]);
										// this is process flow logic
										if (Rnext_21[i_id_level_local] < process_flow_alpha*R0_21[i_id_level_local]) {
											// Смысл модификации в том что мы экономим итерации на пресмутере.
											break; // досрочно опускаемся на следующий уровень если он есть конечно.
										}
										else {
											Rprev_21[i_id_level_local] = Rnext_21[i_id_level_local];
										}
									}
								}

								if (iter == nu2) {
#if doubleintprecision == 1
									printf("level %lld limit postsmother iteration is reached\n", i_id_level_local);
#else
									printf("level %d limit postsmother iteration is reached\n", i_id_level_local);
#endif

								}

							}
							else {
								integer nu2_count = nu2;
								if (i_id_level_local == ilevel - 1) {
									// на самом грубом уровне нам требуется точное решение.
									//nu2_count = 100;																																					
								}
								for (integer iter = 0; iter < nu2_count; ++iter) {
									if (bonly_serial) {

										bool bflag56 = false;
										if (my_amg_manager.iFinnest_ilu == 1) {
											if (my_amg_manager.b_ilu_smoothers_in_nnz_n_LE_6) {
												doublerealT dn = 1.0*n_a[i_id_level_local];
												doublerealT dnnz = 1.0*nnz_a[i_id_level_local];
												if (dnnz / dn <= dapply_ilu_max_pattern_size) {
													bflag56 = true;
												}
											}
										}

										if (bILU2smoother == 1) {
											seidelq(Amat, static_cast<integer>(1), n_a[i_id_level_local], error_approx_coarse[i_id_level_local - 1], residual_coarse[i_id_level_local - 1], row_ptr_start, row_ptr_end, in_a_loc + n_a[i_id_level_local - 1], F_false_C_true, static_cast<integer>(1), i_id_level_local);
											residualq2(Amat, static_cast<integer>(1), n_a[i_id_level_local], error_approx_coarse[i_id_level_local - 1], residual_coarse[i_id_level_local - 1], row_ptr_start, row_ptr_end, in_a_loc + n_a[i_id_level_local - 1], residual_fine[i_id_level_local], diag[i_id_level_local]);
#pragma omp parallel for
											for (integer i = 1; i <= n_a[i_id_level_local]; ++i) {
												milu0[i_id_level_local].zbuf[i] = residual_fine[i_id_level_local][i];
											}
											lusol_1patchforRUMBA(n_a[i_id_level_local], milu0[i_id_level_local].zbuf, milu0[i_id_level_local].zbuf2, milu0[i_id_level_local]);
#pragma omp parallel for
											for (integer i = 1; i <= n_a[i_id_level_local]; ++i) {
												error_approx_coarse[i_id_level_local - 1][i] += milu0[i_id_level_local].zbuf2[i];
											}
										}
										else if (1 && ((bILU2smoother == 2) || bflag56)) {
											seidelq(Amat, static_cast<integer>(1), n_a[i_id_level_local], error_approx_coarse[i_id_level_local - 1], residual_coarse[i_id_level_local - 1], row_ptr_start, row_ptr_end, in_a_loc + n_a[i_id_level_local - 1], F_false_C_true, static_cast<integer>(1), i_id_level_local);
											residualq2(Amat, static_cast<integer>(1), n_a[i_id_level_local], error_approx_coarse[i_id_level_local - 1], residual_coarse[i_id_level_local - 1], row_ptr_start, row_ptr_end, in_a_loc + n_a[i_id_level_local - 1], residual_fine[i_id_level_local], diag[i_id_level_local]);
#pragma omp parallel for
											for (integer i = 1; i <= n_a[i_id_level_local]; ++i) {
												milu2[i_id_level_local].zbuf[i] = residual_fine[i_id_level_local][i];
											}
											lusol_1patchforRUMBA(n_a[i_id_level_local], milu2[i_id_level_local].zbuf, milu2[i_id_level_local].zbuf2, milu2[i_id_level_local]);
#pragma omp parallel for
											for (integer i = 1; i <= n_a[i_id_level_local]; ++i) {
												error_approx_coarse[i_id_level_local - 1][i] += milu2[i_id_level_local].zbuf2[i];
											}
										}
										else {
											seidelq(Amat, static_cast<integer>(1), n_a[i_id_level_local], error_approx_coarse[i_id_level_local - 1], residual_coarse[i_id_level_local - 1], row_ptr_start, row_ptr_end, in_a_loc + n_a[i_id_level_local - 1], F_false_C_true, static_cast<integer>(1), i_id_level_local);
										}
									}
									else {
										seidelq(Amat, static_cast<integer>(1), n_a[i_id_level_local], error_approx_coarse[i_id_level_local - 1], residual_coarse[i_id_level_local - 1], nested_desection[i_id_level_local], row_ptr_start, row_ptr_end, in_a_loc + n_a[i_id_level_local - 1], F_false_C_true, static_cast<integer>(1));
									}
								}
							}


						}

						// prolongation
						// residual_r
						//doublerealT *error_approx_fine[i_id_level_local - 1] = new doublerealT[n_a[i_id_level_local - 1] + 1];
#pragma omp parallel for
						for (integer ii = 1; ii <= n_a[i_id_level_local - 1]; ++ii) {
							error_approx_fine[i_id_level_local - 1][ii] = 0.0;
						}

						integer innz_aRP_loc = 1;
						for (integer i_72 = 0; i_72 < i_id_level_local - 1; ++i_72) {
							innz_aRP_loc += nnz_aRP[i_72];
						}
						prolongation(P, innz_aRP_loc, innz_aRP_loc - 1 + nnz_aRP[i_id_level_local - 1],  error_approx_fine[i_id_level_local - 1], error_approx_coarse[i_id_level_local - 1], n_a[i_id_level_local - 1]);



						// correction
#pragma omp parallel for
						for (integer ii = 1; ii <= n_a[i_id_level_local - 1]; ++ii) {
							error_approx_coarse[i_id_level_local - 2][ii] += error_approx_fine[i_id_level_local - 1][ii];
						}

						// free
						//delete[] error_approx_fine[i_id_level_local - 1];
						//delete[] error_approx_coarse[i_id_level_local - 1];
						//delete[] residual_coarse[i_id_level_local - 1];
						//delete[] residual_fine[i_id_level_local - 1];

					} // 21
				}






				// post smothing
				if (process_flow_logic) {


					// smother
					integer iter = 0;
					for (iter = 0; iter < nu2; ++iter) {
						//quick seidel
						if (bonly_serial) {
							seidelq(Amat, static_cast<integer>(1), n_a[1], error_approx_coarse[0], residual_coarse[0], row_ptr_start, row_ptr_end, n_a[0], F_false_C_true, static_cast<integer>(1),1);
						}
						else {
							seidelq(Amat, static_cast<integer>(1), n_a[1], error_approx_coarse[0], residual_coarse[0], nested_desection[1], row_ptr_start, row_ptr_end, n_a[0], F_false_C_true, static_cast<integer>(1));
						}

						//residualq(Amat, static_cast<integer>(1), n_a[1], error_approx_coarse[0], residual_coarse[0], row_ptr_start, row_ptr_end, n_a[0] , residual_fine[1]);
						residualq2(Amat, static_cast<integer>(1), n_a[1], error_approx_coarse[0], residual_coarse[0], row_ptr_start, row_ptr_end, n_a[0], residual_fine[1], diag[1]);


						Rnext_1 = norma(residual_fine[1], n_a[1]);
						// this is process flow logic
						if (Rnext_1 < process_flow_alpha*R0_1) {
							// Смысл модификации в том что мы экономим итерации на пресмутере.
							break; // досрочно опускаемся на следующий уровень если он есть конечно.
						}
						else {
							Rprev_1 = Rnext_1;
						}
					}

					if (iter == nu2) {
						printf("level 1 limit postsmother iteration is reached\n");
					}

				}
				else {
					for (integer iter = 0; iter < nu2; ++iter) {
						if (bonly_serial) {
							if (bILU2smoother == 1) {
								seidelq(Amat, static_cast<integer>(1), n_a[1], error_approx_coarse[0], residual_coarse[0], row_ptr_start, row_ptr_end, n_a[0], F_false_C_true, static_cast<integer>(1),1);
								residualq2(Amat, static_cast<integer>(1), n_a[1], error_approx_coarse[0], residual_coarse[0], row_ptr_start, row_ptr_end, n_a[0], residual_fine[1], diag[1]);
#pragma omp parallel for
								for (integer i = 1; i <= n_a[1]; ++i) {
									milu0[1].zbuf[i] = residual_fine[1][i];
								}
								lusol_1patchforRUMBA(n_a[1], milu0[1].zbuf, milu0[1].zbuf2, milu0[1]);
#pragma omp parallel for
								for (integer i = 1; i <= n_a[1]; ++i) {
									error_approx_coarse[0][i] += milu0[1].zbuf2[i];
								}
							}
							else if (1 && bILU2smoother == 2) {
								seidelq(Amat, static_cast<integer>(1), n_a[1], error_approx_coarse[0], residual_coarse[0], row_ptr_start, row_ptr_end, n_a[0], F_false_C_true, static_cast<integer>(1),1);
								residualq2(Amat, static_cast<integer>(1), n_a[1], error_approx_coarse[0], residual_coarse[0], row_ptr_start, row_ptr_end, n_a[0], residual_fine[1], diag[1]);
#pragma omp parallel for
								for (integer i = 1; i <= n_a[1]; ++i) {
									milu2[1].zbuf[i] = residual_fine[1][i];
								}
								lusol_1patchforRUMBA(n_a[1], milu2[1].zbuf, milu2[1].zbuf2, milu2[1]);
#pragma omp parallel for
								for (integer i = 1; i <= n_a[1]; ++i) {
									error_approx_coarse[0][i] += milu2[1].zbuf2[i];
								}
							}
							else {
								seidelq(Amat, static_cast<integer>(1), n_a[1], error_approx_coarse[0], residual_coarse[0], row_ptr_start, row_ptr_end, n_a[0], F_false_C_true, static_cast<integer>(1),1);
							}
						}
						else {
							seidelq(Amat, static_cast<integer>(1), n_a[1], error_approx_coarse[0], residual_coarse[0], nested_desection[1], row_ptr_start, row_ptr_end, n_a[0], F_false_C_true, static_cast<integer>(1));
						}
					}
				}

				delete[] R0_21;
				delete[] Rprev_21;
				delete[] Rnext_21;

				R0_21 = nullptr;
				Rprev_21 = nullptr;
				Rnext_21 = nullptr;

			}

			move_up(nu1, nu2);

			// prolongation
			// residual_r
			//doublerealT *error_approx_fine[0] = new doublerealT[n_a[0] + 1];
#pragma omp parallel for
			for (integer ii = 1; ii <= n_a[0]; ++ii) {
				error_approx_fine[0][ii] = 0.0;
			}

			prolongation(P, static_cast<integer>(1), nnz_aRP[0],  error_approx_fine[0], error_approx_coarse[0], n_a[0]);

			// correction
#pragma omp parallel for
			for (integer ii = 1; ii <= n_a[0]; ++ii) {
				//if (row_ptr_start[ii] != row_ptr_end[ii]) {
					// Не условие Дирихле.
				z76[ii] += error_approx_fine[0][ii];
				//	}
			}

			// free
			//delete[] error_approx_fine[0];
			//delete[] error_approx_coarse[0];
			//delete[] residual_coarse[0];
			//delete[] residual_fine[0];
		}



		//doublerealT R0_0 = 0.0;
		//doublerealT Rprev_0 = 0.0, Rnext_0 = 0.0;
		// post smothing
		if (process_flow_logic) {
			// calculate initial residual.
			//residualq<doublereal>(Amat, static_cast<integer>(1), n_a[0], z76, s76, row_ptr_start, row_ptr_end, 0, residual_fine[0]);
			residualq2(Amat, static_cast<integer>(1), n_a[0], z76, s76, row_ptr_start, row_ptr_end, 0, residual_fine[0], diag[0]);
			Rprev_0 = norma(residual_fine[0], n_a[0]);

			// smother
			integer iter = 0;
			for (iter = 0; iter < nFinestSweeps; ++iter) {
				//quick seidel
				integer i_end_n0 = n_a[0];
				if (bonly_serial) {
					seidelq<doublereal>(Amat, static_cast<integer>(1), i_end_n0, z76, s76, row_ptr_start, row_ptr_end, 0, F_false_C_true, static_cast<integer>(1), iVar, 0);
				}
				else {
					seidelq<doublereal>(Amat, static_cast<integer>(1), i_end_n0, z76, s76, nested_desection[0], row_ptr_start, row_ptr_end, 0, F_false_C_true, static_cast<integer>(1));
				}
				//residualq<doublereal>(Amat, static_cast<integer>(1), i_end_n0, z76, s76, row_ptr_start, row_ptr_end, 0, residual_fine[0]);
				residualq2(Amat, static_cast<integer>(1), n_a[0], z76, s76, row_ptr_start, row_ptr_end, 0, residual_fine[0], diag[0]);
				Rnext_0 = norma(residual_fine[0], n_a[0]);
				// this is process flow logic
				if (Rnext_0 < process_flow_alpha*Rprev_0) {
					// Смысл модификации в том что мы экономим итерации на пресмутере.
					break; // досрочно опускаемся на следующий уровень если он есть конечно.
				}
				else {
					Rprev_0 = Rnext_0;
				}
			}
			if (iter == nFinestSweeps) {
				printf("level 0 limit postsmother iteration is reached\n");
			}

		}
		else {
			// nFinnestSweeps new logic 14 jan 2016.
			// smother
			for (integer iter = 0; iter < nFinestSweeps; ++iter) {
				//seidel<doublereal>(Amat, 1, nnz_a[0], z76, s76, flag, n_a[0]);
				//quick seidel
				if (bonly_serial) {
					if (bILU2smoother == 1) {
						seidelq<doublereal>(Amat, static_cast<integer>(1), n_a[0], z76, s76, row_ptr_start, row_ptr_end, 0, F_false_C_true, static_cast<integer>(1), iVar, 0);
						residualq2(Amat, static_cast<integer>(1), n_a[0], z76, s76, row_ptr_start, row_ptr_end, 0, residual_fine[0], diag[0]);
#pragma omp parallel for
						for (integer i = 1; i <= n_a[0]; ++i) {
							milu0[0].zbuf[i] = residual_fine[0][i];
						}
						lusol_1patchforRUMBA(n_a[0], milu0[0].zbuf, milu0[0].zbuf2, milu0[0]);
#pragma omp parallel for
						for (integer i = 1; i <= n_a[0]; ++i) {
							z76[i] += milu0[0].zbuf2[i];
						}
					}
					else if ((bILU2smoother == 2) || (my_amg_manager.iFinnest_ilu == 1)) {
						seidelq<doublereal>(Amat, static_cast<integer>(1), n_a[0], z76, s76, row_ptr_start, row_ptr_end, 0, F_false_C_true, static_cast<integer>(1), iVar, 0);
						residualq2(Amat, 1, n_a[0], z76, s76, row_ptr_start, row_ptr_end, 0, residual_fine[0], diag[0]);
#pragma omp parallel for
						for (integer i = 1; i <= n_a[0]; ++i) {
							milu2[0].zbuf[i] = residual_fine[0][i];
						}
						lusol_1patchforRUMBA(n_a[0], milu2[0].zbuf, milu2[0].zbuf2, milu2[0]);
#pragma omp parallel for
						for (integer i = 1; i <= n_a[0]; ++i) {
							z76[i] += milu2[0].zbuf2[i];
						}
					}
					else {
						seidelq<doublereal>(Amat, static_cast<integer>(1), n_a[0], z76, s76, row_ptr_start, row_ptr_end, 0, F_false_C_true, static_cast<integer>(1), iVar, 0);
					}
				}
				else {
					seidelq<doublereal>(Amat, static_cast<integer>(1), n_a[0], z76, s76, nested_desection[0], row_ptr_start, row_ptr_end, static_cast<integer>(0), F_false_C_true, static_cast<integer>(1));
				}
			}
		}

	}
} // V_cycle_solve

// 08.01.2018 Перенесено в отдельную функцию, т.к. используется неоднократно. Передаётся большое число параметров.
//  A*z76=s76; 891 строка.
// Для плохих строк при их наличии мы домножаем правую часть на минус 1.0
template <typename doublerealT>
void V_cycle_solve(Ak2& Amat, doublereal*& z76, doublereal*& s76, bool process_flow_logic, integer*& row_ptr_start,
	integer*& row_ptr_end, doublereal**& residual_fine, doublereal**& diag, 
	doublereal**& diag_minus_one, integer* n_a, bool bonly_serial,
	doublerealT process_flow_beta, bool*& F_false_C_true, integer& nu1, integer& nu2, integer bILU2smoother,
	integer ilevel, integer inumberVcyclelocbicgstab, INIT_SELECTOR_CASE_CAMG_RUMBAv_0_14 imyinit, const integer idim_diag,
	LEVEL_ADDITIONAL_DATA*& milu2, LEVEL_ADDITIONAL_DATA0* milu0, bool**& nested_desection,
	Ak1*& P, // prolongation он же restriction (метод не чувствителен к сортировке).
	integer* nnz_aRP, doublereal**& residual_coarse, integer igam, integer* nnz_a,
	doublereal**& error_approx_coarse, doublerealT dapply_ilu_max_pattern_size,
	doublerealT process_flow_alpha, doublereal**& error_approx_fine,
	integer nFinestSweeps, integer ibsp_length, integer iVar) {
	

	// Один V - цикл алгебраического многосеточного метода.
	// A*z76=s76;

	// надо увеличивать nu1, nu2 с 1,2 до 5 наверно.

	for (integer i_1 = 0; i_1 < inumberVcyclelocbicgstab; ++i_1)
	{


		doublerealT R0_0 = 0.0;
		doublerealT Rprev_0 = 0.0, Rnext_0 = 0.0;
		if (process_flow_logic) {
			// calculate initial residual.
			//residualq(Amat, 1, n_a[0], z76, s76, row_ptr_start, row_ptr_end, 0, residual_fine[0]);
			residualq2(Amat, 1, n_a[0], z76, s76, row_ptr_start, row_ptr_end, 0, residual_fine[0], diag[0], diag_minus_one[0]);
			R0_0 = static_cast<doublerealT>(norma(residual_fine[0], n_a[0]));
			Rprev_0 = R0_0;
			

			// smother
			integer iter = 0;
			for (iter = 0; iter < nu1; ++iter) {
				//quick seidel
				if (bonly_serial) {
					seidelq<doublereal>(Amat, 1, n_a[0], z76, s76, row_ptr_start, row_ptr_end, 0, F_false_C_true, 0, diag_minus_one[0], iVar, 0, ibsp_length);
				}
				else {
					seidelq<doublereal>(Amat, 1, n_a[0], z76, s76, nested_desection[0], row_ptr_start, row_ptr_end, 0, F_false_C_true, 0, diag_minus_one[0]);
				}
				//residualq(Amat, 1, n_a[0], z76, s76, row_ptr_start, row_ptr_end, 0, residual_fine[0]);
				residualq2(Amat, 1, n_a[0], z76, s76, row_ptr_start, row_ptr_end, 0, residual_fine[0], diag[0], diag_minus_one[0]);
				Rnext_0 = static_cast<doublerealT>(norma(residual_fine[0], n_a[0]));
				// this is process flow logic
				if (Rnext_0 > process_flow_beta * Rprev_0) {
					// Смысл модификации в том что мы экономим итерации на пресмутере.
					break; // досрочно опускаемся на следующий уровень если он есть конечно.
				}
				else {
					Rprev_0 = Rnext_0;
				}
			}
			if (iter == nu1) {
				printf("level 0 limit presmother iteration is reached\n");
			}

		}
		else {
			// smoother
			for (integer iter = 0; iter < nu1; ++iter) {
				//seidel(Amat, 1, nnz_a[0], z76, s76, flag, n_a[0]);
				//quick seidel
				if (bonly_serial) {
					if (bILU2smoother == 1) {
						// ILU0
						seidelq<doublereal>(Amat, 1, n_a[0], z76, s76, row_ptr_start, row_ptr_end, 0, F_false_C_true, 0, diag_minus_one[0], iVar, 0, ibsp_length);
						residualq2(Amat, 1, n_a[0], z76, s76, row_ptr_start, row_ptr_end, 0, residual_fine[0], diag[0], diag_minus_one[0]);
#pragma omp parallel for
						for (integer i = 1; i <= n_a[0]; ++i) {
							milu0[0].zbuf[i] = residual_fine[0][i];
						}
						lusol_1patchforRUMBA(n_a[0], milu0[0].zbuf, milu0[0].zbuf2, milu0[0]);
#pragma omp parallel for
						for (integer i = 1; i <= n_a[0]; ++i) {
							z76[i] += milu0[0].zbuf2[i];
						}
					}
					else if ((bILU2smoother == 2) || (my_amg_manager.iFinnest_ilu == 1)) {
						seidelq<doublereal>(Amat, 1, n_a[0], z76, s76, row_ptr_start, row_ptr_end, 0, F_false_C_true, 0, diag_minus_one[0], iVar, 0, ibsp_length);
						residualq2(Amat, 1, n_a[0], z76, s76, row_ptr_start, row_ptr_end, 0, residual_fine[0], diag[0], diag_minus_one[0]);
#pragma omp parallel for
						for (integer i = 1; i <= n_a[0]; ++i) {
							milu2[0].zbuf[i] = residual_fine[0][i];
						}
						lusol_1patchforRUMBA(n_a[0], milu2[0].zbuf, milu2[0].zbuf2, milu2[0]);
#pragma omp parallel for
						for (integer i = 1; i <= n_a[0]; ++i) {
							z76[i] += milu2[0].zbuf2[i];
						}
					}
					else {
						seidelq<doublereal>(Amat, 1, n_a[0], z76, s76, row_ptr_start, row_ptr_end, 0, F_false_C_true, 0, diag_minus_one[0], iVar, 0, ibsp_length);
					}
				}
				else {
					seidelq<doublereal>(Amat, 1, n_a[0], z76, s76, nested_desection[0], row_ptr_start, row_ptr_end, 0, F_false_C_true, 0, diag_minus_one[0]);
				}
			}
		}

		//exporttecplot(x, n);

		move_down(nu1, nu2);

		if (!process_flow_logic) {
			// residual_r
			//doublerealT *residual_fine[0] = new doublerealT[n_a[0] + 1];
			//residual<doublereal>(Amat, 1, nnz_a[0], z76, s76, flag, n_a[0], residual_fine[0]);
			//residualq<doublereal>(Amat, 1, n_a[0], z76, s76, row_ptr_start, row_ptr_end, 0, residual_fine[0]);
			residualq2(Amat, 1, n_a[0], z76, s76, row_ptr_start, row_ptr_end, 0, residual_fine[0], diag[0], diag_minus_one[0]);
		}



		//if (iprohod%5==0) system("pause");
		if (ilevel > 1) {

			//doublerealT *residual_coarse = new doublerealT[n_a[1] + 1];
			

#pragma omp parallel for
				for (integer ii = 1; ii <= n_a[1]; ++ii) {
					residual_coarse[0][ii] = 0.0;
				}
			

			// restriction
			restriction(P, 1, nnz_aRP[0], residual_fine[0], residual_coarse[0], n_a[1]);

			//std::cout <<"ilevel= "<< 1 << std::endl;
			for (integer ii = 1; ii <= n_a[1]; ++ii) {
				if (residual_coarse[0][ii]!= residual_coarse[0][ii]) {
					std::cout << "residual_coarse = " << residual_coarse[0][ii] << std::endl;
					system("pause");
				}
			}

			// Amat*e=r;
			//doublerealT* error_approx_coarse = new doublerealT[n_a[1] + 1];			

#pragma omp parallel for
			for (integer ii = 1; ii <= n_a[1]; ++ii) {
				error_approx_coarse[0][ii] = 0.0;
			}


			for (integer i_2 = 1; i_2 <= igam; ++i_2) {
				doublerealT R0_1 = 0.0;
				doublerealT Rprev_1 = 0.0, Rnext_1 = 0.0;
				if (process_flow_logic) {
					// calculate initial residual.
					//residualq(Amat, 1, n_a[1], error_approx_coarse[0], residual_coarse[0], row_ptr_start, row_ptr_end, n_a[0], residual_fine[1]);
					residualq2(Amat, 1, n_a[1], error_approx_coarse[0], residual_coarse[0], row_ptr_start, row_ptr_end, n_a[0], residual_fine[1], diag[1], diag_minus_one[1]);

					

					R0_1 = static_cast<doublerealT>(norma(residual_fine[1], n_a[1]));
					Rprev_1 = R0_1;

					// smother
					integer iter = 0;
					for (iter = 0; iter < nu1; ++iter) {
						//quick seidel
						if (bonly_serial) {
							seidelq(Amat, 1, n_a[1], error_approx_coarse[0], residual_coarse[0], row_ptr_start, row_ptr_end, n_a[0], F_false_C_true, 0, diag_minus_one[1], iVar, 1, ibsp_length);
						}
						else {
							seidelq(Amat, 1, n_a[1], error_approx_coarse[0], residual_coarse[0], nested_desection[1], row_ptr_start, row_ptr_end, n_a[0], F_false_C_true, 0, diag_minus_one[1]);
						}

						//residualq(Amat, 1, n_a[1], error_approx_coarse[0], residual_coarse[0], row_ptr_start, row_ptr_end, n_a[0], residual_fine[1]);
						residualq2(Amat, 1, n_a[1], error_approx_coarse[0], residual_coarse[0], row_ptr_start, row_ptr_end, n_a[0], residual_fine[1], diag[1], diag_minus_one[1]);

						Rnext_1 = static_cast<doublerealT>(norma(residual_fine[1], n_a[1]));
						// this is process flow logic
						if (Rnext_1 > process_flow_beta * Rprev_1) {
							// Смысл модификации в том что мы экономим итерации на пресмутере.
							break; // досрочно опускаемся на следующий уровень если он есть конечно.
						}
						else {
							Rprev_1 = Rnext_1;
						}
					}

					if (iter == nu1) {
						printf("level 1 limit presmother iteration is reached\n");
					}

				}
				else {

					// pre smothing
					for (integer iter = 0; iter < nu1; ++iter) {
						//seidel(Amat, 1 + 2 * nnz_a[0], 2 * nnz_a[0] + nnz_a[1], error_approx_coarse, residual_coarse, flag, n_a[1]);
						if (bonly_serial) {
							if (bILU2smoother == 1) {
								seidelq(Amat, 1, n_a[1], error_approx_coarse[0], residual_coarse[0], row_ptr_start, row_ptr_end, n_a[0], F_false_C_true, 0, diag_minus_one[1],iVar, 1, ibsp_length);
								residualq2(Amat, 1, n_a[1], error_approx_coarse[0], residual_coarse[0], row_ptr_start, row_ptr_end, n_a[0], residual_fine[1], diag[1], diag_minus_one[1]);
#pragma omp parallel for
								for (integer i = 1; i <= n_a[1]; ++i) {
									milu0[1].zbuf[i] = residual_fine[1][i];
								}
								lusol_1patchforRUMBA(n_a[1], milu0[1].zbuf, milu0[1].zbuf2, milu0[1]);
#pragma omp parallel for
								for (integer i = 1; i <= n_a[1]; ++i) {
									error_approx_coarse[0][i] += milu0[1].zbuf2[i];
								}


							}
							else if (1 && bILU2smoother == 2) {
								seidelq(Amat, 1, n_a[1], error_approx_coarse[0], residual_coarse[0], row_ptr_start, row_ptr_end, n_a[0], F_false_C_true, 0, diag_minus_one[1], iVar, 1, ibsp_length);
								residualq2(Amat, 1, n_a[1], error_approx_coarse[0], residual_coarse[0], row_ptr_start, row_ptr_end, n_a[0], residual_fine[1], diag[1], diag_minus_one[1]);
#pragma omp parallel for
								for (integer i = 1; i <= n_a[1]; ++i) {
									milu2[1].zbuf[i] = residual_fine[1][i];
								}
								lusol_1patchforRUMBA(n_a[1], milu2[1].zbuf, milu2[1].zbuf2, milu2[1]);
#pragma omp parallel for
								for (integer i = 1; i <= n_a[1]; ++i) {
									error_approx_coarse[0][i] += milu2[1].zbuf2[i];
								}


							}
							else {
								seidelq(Amat, 1, n_a[1], error_approx_coarse[0], residual_coarse[0], row_ptr_start, row_ptr_end, n_a[0], F_false_C_true, 0, diag_minus_one[1], iVar, 1, ibsp_length);
							}
						}
						else {
							seidelq(Amat, 1, n_a[1], error_approx_coarse[0], residual_coarse[0], nested_desection[1], row_ptr_start, row_ptr_end, n_a[0], F_false_C_true, 0, diag_minus_one[1]);
						}
					}
				}


				for (integer ii = 1; ii <= n_a[1]; ++ii) {
					if (error_approx_coarse[0][ii] != error_approx_coarse[0][ii]) {
						std::cout << "posle seidel q: error_approx_coarse = " << error_approx_coarse[0][ii] << std::endl;
						system("pause");
					}
				}

				//std::cout << "idim_diag=" << idim_diag << std::endl;

				doublerealT* R0_21 = new doublerealT[idim_diag];
				doublerealT* Rprev_21 = new doublerealT[idim_diag];
				doublerealT* Rnext_21 = new doublerealT[idim_diag];

				//3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21
				for (integer i_id_level_local = 2; i_id_level_local < idim_diag; ++i_id_level_local) {

					move_down(nu1, nu2);

					//std::cout << "ilevel" << ilevel << " i_id_level_local=" << i_id_level_local << std::endl;

					if (ilevel > i_id_level_local) {
						// Впервые создан универсальный програмный код на произвольное количество уровней вложенности 
						// в solution phase.
						// Нужно только написать универсальный обработчик для R0_21.
						// 4 декабря 2016. 

						for (integer ii = 1; ii <= n_a[i_id_level_local-1]; ++ii) {
							if (residual_coarse[i_id_level_local - 2][ii] != residual_coarse[i_id_level_local - 2][ii]) {
								std::cout << "residual_coarse ilevel-2 apriory = " << residual_coarse[i_id_level_local - 2][ii] << " ii=" << ii << std::endl;
								system("pause");
							}
							if (error_approx_coarse[i_id_level_local - 2][ii] != error_approx_coarse[i_id_level_local - 2][ii]) {
								std::cout << "error_approx_coarse ilevel -2 apriory = " << error_approx_coarse[i_id_level_local - 2][ii] << " ii=" << ii << std::endl;
								system("pause");
							}
						}

						// residual
						if (!process_flow_logic) {
							integer in_a_loc = 0;
							for (integer i_72 = 0; i_72 < i_id_level_local - 1; ++i_72) {
								// sum(n_a[0]+n_a[19];
								in_a_loc += n_a[i_72];
							}

							//for (int i_94 = 78831; i_94 < 78831 + 10; ++i_94) {
								//std::cout << "x="<<error_approx_coarse[i_id_level_local - 2][i_94] << " \n";
							//}

							//doublerealT *residual_fine[i_id_level_local-1] = new doublerealT[n_a[i_id_level_local-1] + 1];
							//residual(Amat, 1 + nnz_a[0] + nnz_a[1] + nnz_a[2] + nnz_a[3] + nnz_a[4] + nnz_a[5] + nnz_a[6]+ nnz_a[7]+ nnz_a[8]+ nnz_a[9]+ nnz_a[10]+ nnz_a[11]+ nnz_a[12]+ nnz_a[13]+ nnz_a[14]+ nnz_a[15]+ nnz_a[16]+nnz_a[17]+nnz_a[18]+nnz_a[19], nnz_a[0] + nnz_a[1] + nnz_a[2] + nnz_a[3] + nnz_a[4] + nnz_a[5] + nnz_a[6] + nnz_a[7]+ nnz_a[8]+ nnz_a[9]+ nnz_a[10]+ nnz_a[11]+ nnz_a[12]+ nnz_a[13]+ nnz_a[14]+ nnz_a[15]+ nnz_a[16]+ nnz_a[17]+nnz_a[18]+nnz_a[19]+nnz_a[20], error_approx_coarse[i_id_level_local - 2], residual_coarse[i_id_level_local - 2], flag, n_a[i_id_level_local - 1], residual_fine[i_id_level_local - 1]);
							//residualq(Amat, 1, n_a[i_id_level_local-1], error_approx_coarse[i_id_level_local - 2], residual_coarse[i_id_level_local - 2], row_ptr_start, row_ptr_end,  in_a_loc , residual_fine[i_id_level_local - 1]);
							//residualq2(Amat, 1, n_a[i_id_level_local-1], error_approx_coarse[i_id_level_local - 2], residual_coarse[i_id_level_local - 2], row_ptr_start, row_ptr_end, in_a_loc, residual_fine[i_id_level_local - 1], diag[i_id_level_local - 1]);
							residualq2(Amat, 1, n_a[i_id_level_local - 1], error_approx_coarse[i_id_level_local - 2],
								residual_coarse[i_id_level_local - 2], row_ptr_start, row_ptr_end, in_a_loc,
								residual_fine[i_id_level_local - 1], diag[i_id_level_local - 1],
								diag_minus_one[i_id_level_local - 1]);

						}

						//doublerealT *residual_coarse[i_id_level_local-1] = new doublerealT[n_a[i_id_level_local] + 1];

						// restriction

						integer innz_aRP_loc = 1;
						for (integer i_72 = 0; i_72 < i_id_level_local - 1; ++i_72) {
							innz_aRP_loc += nnz_aRP[i_72];
						}

						for (integer ii = 1; ii <= n_a[i_id_level_local]; ++ii) {
							if (residual_coarse[i_id_level_local - 1][ii] != residual_coarse[i_id_level_local - 1][ii]) {
								std::cout << "residual_coarse apriory = " << residual_coarse[i_id_level_local - 1][ii] << " ii=" << ii << "i_id_level_local = " << i_id_level_local << std::endl;
								system("pause");
							}
							if (residual_fine[i_id_level_local - 1][ii] != residual_fine[i_id_level_local - 1][ii]) {
								std::cout << "residual_fine apriory = " << residual_fine[i_id_level_local - 1][ii] << " ii=" << ii << std::endl;
								system("pause");
							}
						}

						//restriction(P, 1 + nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2] + nnz_aRP[3] + nnz_aRP[4] + nnz_aRP[5] + nnz_aRP[6] + nnz_aRP[7] + nnz_aRP[8] + nnz_aRP[9] + nnz_aRP[10] + nnz_aRP[11] + nnz_aRP[12] + nnz_aRP[13] + nnz_aRP[14] + nnz_aRP[15] + nnz_aRP[16] + nnz_aRP[17] + nnz_aRP[18] + nnz_aRP[19], nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2] + nnz_aRP[3] + nnz_aRP[4] + nnz_aRP[5] + nnz_aRP[6] + nnz_aRP[7] + nnz_aRP[8] + nnz_aRP[9] + nnz_aRP[10] + nnz_aRP[11] + nnz_aRP[12] + nnz_aRP[13] + nnz_aRP[14] + nnz_aRP[15] + nnz_aRP[16] + nnz_aRP[17] + nnz_aRP[18] + nnz_aRP[19] + nnz_aRP[20], flag, residual_fine[i_id_level_local - 1], residual_coarse[i_id_level_local - 1], n_a[i_id_level_local - 1], n_a[i_id_level_local]);
						restriction(P, innz_aRP_loc, innz_aRP_loc - 1 + nnz_aRP[i_id_level_local - 1], residual_fine[i_id_level_local - 1], residual_coarse[i_id_level_local - 1], n_a[i_id_level_local]);

						//std::cout << "ilevel= " << i_id_level_local << std::endl;
						for (integer ii = 1; ii <= n_a[i_id_level_local]; ++ii) {
							if (residual_coarse[i_id_level_local - 1][ii] != residual_coarse[i_id_level_local - 1][ii]) {
								std::cout << "residual_coarse = " << residual_coarse[i_id_level_local - 1][ii] << " ii="<<ii<<std::endl;
								system("pause");
							}
						}

						// Amat*e=r;
						//doublerealT* error_approx_coarse[i_id_level_local-1] = new doublerealT[n_a[i_id_level_local] + 1];

						
							

#pragma omp parallel for
					                for (integer ii = 1; ii <= n_a[i_id_level_local]; ++ii) {
										error_approx_coarse[i_id_level_local - 1][ii] = 0.0;
									}
						
						


						//for (integer i_2 = 1; i_2 <= igam; ++i_2)
						{
							// pre smothing
							//doublerealT R0_21 = 0.0;
							//doublerealT Rprev_21 = 0.0, Rnext_21 = 0.0;
							R0_21[i_id_level_local] = 0.0;
							Rprev_21[i_id_level_local] = 0.0;
							Rnext_21[i_id_level_local] = 0.0;
							if (process_flow_logic) {

								integer in_a_loc = 0;
								for (integer i_72 = 0; i_72 < i_id_level_local - 1; ++i_72) {
									// sum(n_a[0]+n_a[19];
									in_a_loc += n_a[i_72];
								}

								// calculate initial residual.
								//residualq(Amat, 1, n_a[i_id_level_local], error_approx_coarse[i_id_level_local-1], residual_coarse[i_id_level_local-1], row_ptr_start, row_ptr_end,in_a_loc + n_a[i_id_level_local-1] , residual_fine[i_id_level_local]);
								residualq2(Amat, 1, n_a[i_id_level_local], error_approx_coarse[i_id_level_local - 1], residual_coarse[i_id_level_local - 1], row_ptr_start, row_ptr_end, in_a_loc + n_a[i_id_level_local - 1], residual_fine[i_id_level_local], diag[i_id_level_local], diag_minus_one[i_id_level_local - 1]);



								R0_21[i_id_level_local] = static_cast<doublerealT>(norma(residual_fine[i_id_level_local], n_a[i_id_level_local]));
								Rprev_21[i_id_level_local] = R0_21[i_id_level_local];



								// smother
								integer iter = 0;
								integer nu1_count = nu1;
								if (i_id_level_local == ilevel - 1) {
									// на самом грубом уровне нам требуется точное решение.
									//nu1_count = 100;

									printf("Direct method is start.../n");

									IMatrix sparseS; // разреженная матрица в формате IMatrix
									initIMatrix(&sparseS, n_a[i_id_level_local]);

									const integer jstart = 1 + in_a_loc + n_a[i_id_level_local - 1];
									const integer jend = n_a[i_id_level_local] + in_a_loc + n_a[i_id_level_local - 1];

									for (integer j = jstart; j <= jend; ++j) {

										//setValueIMatrix(&sparseS, Amat.i[row_ptr_start[j]] - 1, Amat.i[row_ptr_start[j]] - 1, 1.0 / Amat.aij[row_ptr_start[j]]);
										setValueIMatrix(&sparseS, j - jstart, j - jstart, 1.0 / Amat.aij[row_ptr_start[j]]);//3.02.2019
										//if (Amat.i[row_ptr_start[j]] != j - jstart + 1) {
											//printf("Amat.i[row_ptr_start[j]]=%lld j - jstart + 1=%lld\n", Amat.i[row_ptr_start[j]],j - jstart + 1);
											//system("pause");
										//}
										const doublerealT nonzeroEPS = static_cast <doublerealT>(1e-37); // для отделения вещественного нуля

										const integer istart = row_ptr_start[j] + 1;
										const integer iend = row_ptr_end[j];

										for (integer i = istart; i <= iend; ++i)
										{

											if ((fabs(Amat.aij[i]) > nonzeroEPS)) {
												//setValueIMatrix(&sparseS, Amat.i[i] - 1, Amat.j[i] - 1, Amat.aij[i]);
												setValueIMatrix(&sparseS, i - istart, Amat.j[i] - 1, Amat.aij[i]);
											}

										}
									}

									doublereal* dX025 = new doublereal[n_a[i_id_level_local]];
									doublereal* dV25 = new doublereal[n_a[i_id_level_local]];
#pragma omp parallel for
									for (integer i = 0; i < n_a[i_id_level_local]; ++i) {
										dX025[i] = error_approx_coarse[i_id_level_local - 1][i + 1];
										dV25[i] = diag_minus_one[i_id_level_local - 1][i + 1] *residual_coarse[i_id_level_local - 1][i + 1];
									}

									// главный метод, возвращающий решение x,
									// принимает вектор свободных членов b и 
									// квадратную матрицу xO в специальном разреженном формате.
									// реализация без барьера и итерационного уточнения.
									calculateSPARSEgaussArray(&sparseS, dX025, dV25);
#pragma omp parallel for
									for (integer i = 0; i < n_a[i_id_level_local]; ++i) {
										error_approx_coarse[i_id_level_local - 1][i + 1] = dX025[i];
									}

									delete[] dX025;
									delete[] dV25;


									freeIMatrix(&sparseS);

									//residualq(Amat, 1, n_a[i_id_level_local], error_approx_coarse[i_id_level_local-1], residual_coarse[i_id_level_local-1], row_ptr_start, row_ptr_end, in_a_loc + n_a[i_id_level_local - 1], residual_fine[i_id_level_local]);
									residualq2(Amat, 1, n_a[i_id_level_local], error_approx_coarse[i_id_level_local - 1], residual_coarse[i_id_level_local - 1], row_ptr_start, row_ptr_end, in_a_loc + n_a[i_id_level_local - 1], residual_fine[i_id_level_local], diag[i_id_level_local], diag_minus_one[i_id_level_local]);


									Rnext_21[i_id_level_local] = static_cast<doublerealT>(norma(residual_fine[i_id_level_local], n_a[i_id_level_local]));
									// this is process flow logic
									if (Rnext_21[i_id_level_local] > process_flow_beta * Rprev_21[i_id_level_local]) {
										// Смысл модификации в том что мы экономим итерации на пресмутере.
										break; // досрочно опускаемся на следующий уровень если он есть конечно.
									}
									else {
										Rprev_21[i_id_level_local] = Rnext_21[i_id_level_local];
									}
								}
								else {
									for (iter = 0; iter < nu1_count; ++iter) {
										//quick seidel
										if (bonly_serial) {
											seidelq(Amat, 1, n_a[i_id_level_local], error_approx_coarse[i_id_level_local - 1], residual_coarse[i_id_level_local - 1], row_ptr_start, row_ptr_end, in_a_loc + n_a[i_id_level_local - 1], F_false_C_true, 0, diag_minus_one[i_id_level_local], iVar,  i_id_level_local, ibsp_length);
										}
										else {
											seidelq(Amat, 1, n_a[i_id_level_local], error_approx_coarse[i_id_level_local - 1], residual_coarse[i_id_level_local - 1], nested_desection[i_id_level_local], row_ptr_start, row_ptr_end, in_a_loc + n_a[i_id_level_local - 1], F_false_C_true, 0, diag_minus_one[i_id_level_local]);
										}

										//residualq(Amat, 1, n_a[i_id_level_local], error_approx_coarse[i_id_level_local-1], residual_coarse[i_id_level_local-1], row_ptr_start, row_ptr_end, in_a_loc + n_a[i_id_level_local - 1], residual_fine[i_id_level_local]);
										residualq2(Amat, 1, n_a[i_id_level_local], error_approx_coarse[i_id_level_local - 1], residual_coarse[i_id_level_local - 1], row_ptr_start, row_ptr_end, in_a_loc + n_a[i_id_level_local - 1], residual_fine[i_id_level_local], diag[i_id_level_local], diag_minus_one[i_id_level_local]);


										Rnext_21[i_id_level_local] = static_cast<doublerealT>(norma(residual_fine[i_id_level_local], n_a[i_id_level_local]));
										// this is process flow logic
										if (Rnext_21[i_id_level_local] > process_flow_beta * Rprev_21[i_id_level_local]) {
											// Смысл модификации в том что мы экономим итерации на пресмутере.
											break; // досрочно опускаемся на следующий уровень если он есть конечно.
										}
										else {
											Rprev_21[i_id_level_local] = Rnext_21[i_id_level_local];
										}
									}
								}

								if (iter == nu1) {
#if doubleintprecision == 1
									printf("level %lld limit presmother iteration is reached\n", i_id_level_local);
#else
									printf("level %d limit presmother iteration is reached\n", i_id_level_local);
#endif

								}

							}
							else {
								integer nu1_count = nu1;
								if (i_id_level_local == ilevel - 1) {
									// на самом грубом уровне нам требуется точное решение.
									//nu1_count = 100;
								}
								for (integer iter = 0; iter < nu1_count; ++iter) {

									integer in_a_loc = 0;
									for (integer i_72 = 0; i_72 < i_id_level_local - 1; ++i_72) {
										// sum(n_a[0]+n_a[19];
										in_a_loc += n_a[i_72];
									}

									if (bonly_serial) {

										bool bflag56 = false;
										if (my_amg_manager.iFinnest_ilu == 1) {
											if (my_amg_manager.b_ilu_smoothers_in_nnz_n_LE_6) {
												doublerealT dn = static_cast<doublerealT>(1.0 * n_a[i_id_level_local]);
												doublerealT dnnz = static_cast<doublerealT>(1.0 * nnz_a[i_id_level_local]);
												if (dnnz / dn <= dapply_ilu_max_pattern_size) {
													bflag56 = true;
												}
											}
										}

										if (bILU2smoother == 1) {
											seidelq(Amat, 1, n_a[i_id_level_local], error_approx_coarse[i_id_level_local - 1], residual_coarse[i_id_level_local - 1], row_ptr_start, row_ptr_end, in_a_loc + n_a[i_id_level_local - 1], F_false_C_true, 0, diag_minus_one[i_id_level_local], iVar, i_id_level_local, ibsp_length);
											residualq2(Amat, 1, n_a[i_id_level_local], error_approx_coarse[i_id_level_local - 1], residual_coarse[i_id_level_local - 1], row_ptr_start, row_ptr_end, in_a_loc + n_a[i_id_level_local - 1], residual_fine[i_id_level_local], diag[i_id_level_local], diag_minus_one[i_id_level_local]);
#pragma omp parallel for
											for (integer i = 1; i <= n_a[i_id_level_local]; ++i) {
												milu0[i_id_level_local].zbuf[i] = residual_fine[i_id_level_local][i];
											}
											lusol_1patchforRUMBA(n_a[i_id_level_local], milu0[i_id_level_local].zbuf, milu0[i_id_level_local].zbuf2, milu0[i_id_level_local]);
#pragma omp parallel for
											for (integer i = 1; i <= n_a[i_id_level_local]; ++i) {
												error_approx_coarse[i_id_level_local - 1][i] += milu0[i_id_level_local].zbuf2[i];
											}

										}
										else if (1 && ((bILU2smoother == 2) || bflag56)) {
											seidelq(Amat, 1, n_a[i_id_level_local], error_approx_coarse[i_id_level_local - 1], residual_coarse[i_id_level_local - 1], row_ptr_start, row_ptr_end, in_a_loc + n_a[i_id_level_local - 1], F_false_C_true, 0, diag_minus_one[i_id_level_local], iVar, i_id_level_local, ibsp_length);
											residualq2(Amat, 1, n_a[i_id_level_local], error_approx_coarse[i_id_level_local - 1], residual_coarse[i_id_level_local - 1], row_ptr_start, row_ptr_end, in_a_loc + n_a[i_id_level_local - 1], residual_fine[i_id_level_local], diag[i_id_level_local], diag_minus_one[i_id_level_local]);
#pragma omp parallel for
											for (integer i = 1; i <= n_a[i_id_level_local]; ++i) {
												milu2[i_id_level_local].zbuf[i] = residual_fine[i_id_level_local][i];
											}
											lusol_1patchforRUMBA(n_a[i_id_level_local], milu2[i_id_level_local].zbuf, milu2[i_id_level_local].zbuf2, milu2[i_id_level_local]);
#pragma omp parallel for
											for (integer i = 1; i <= n_a[i_id_level_local]; ++i) {
												error_approx_coarse[i_id_level_local - 1][i] += milu2[i_id_level_local].zbuf2[i];
											}

										}
										else {
											seidelq(Amat, 1, n_a[i_id_level_local], error_approx_coarse[i_id_level_local - 1], residual_coarse[i_id_level_local - 1], row_ptr_start, row_ptr_end, in_a_loc + n_a[i_id_level_local - 1], F_false_C_true, 0, diag_minus_one[i_id_level_local], iVar,  i_id_level_local, ibsp_length);
										}
									}
									else {
										seidelq(Amat, 1, n_a[i_id_level_local], error_approx_coarse[i_id_level_local - 1], residual_coarse[i_id_level_local - 1], nested_desection[i_id_level_local], row_ptr_start, row_ptr_end, in_a_loc + n_a[i_id_level_local - 1], F_false_C_true, 0, diag_minus_one[i_id_level_local]);
									}
								}
							}
						}

					}

				}







				//center
				// ЭТО сердцевина SOLUTION PHASE.

				// 3.12.2016 (осталось один нижний и один верхний).


				// 21
				for (integer i_id_level_local = idim_diag - 1; i_id_level_local >= 2; i_id_level_local--) {

					if (ilevel > i_id_level_local) {

						{

							integer in_a_loc = 0;
							for (integer i_72 = 0; i_72 < i_id_level_local - 1; ++i_72) {
								// sum(n_a[0]+n_a[19];
								in_a_loc += n_a[i_72];
							}

							// post smoothing
							// doublerealT R0_20 = 0.0;
							///doublerealT Rprev_20 = 0.0, Rnext_20 = 0.0;
							if (process_flow_logic) {
								// calculate initial residual.
								//residualq(Amat, 1, n_a[i_id_level_local], error_approx_coarse[i_id_level_local - 1], residual_coarse[i_id_level_local - 1], row_ptr_start, row_ptr_end, in_a_loc  + n_a[i_id_level_local - 1], residual_fine[i_id_level_local]);



								//Rprev_21[i_id_level_local] = norma(residual_fine[i_id_level_local], n_a[i_id_level_local]);


								// smother
								integer iter = 0;
								integer nu2_count = nu2;
								if (i_id_level_local == ilevel - 1) {
									// на самом грубом уровне нам требуется точное решение.
									//nu2_count = 100;																												
								}
								else {
									for (iter = 0; iter < nu2_count; ++iter) {
										//quick seidel
										if (bonly_serial) {
											seidelq(Amat, 1, n_a[i_id_level_local], error_approx_coarse[i_id_level_local - 1], residual_coarse[i_id_level_local - 1], row_ptr_start, row_ptr_end, in_a_loc + n_a[i_id_level_local - 1], F_false_C_true, 1, diag_minus_one[i_id_level_local], iVar, i_id_level_local, ibsp_length);
										}
										else {
											seidelq(Amat, 1, n_a[i_id_level_local], error_approx_coarse[i_id_level_local - 1], residual_coarse[i_id_level_local - 1], nested_desection[i_id_level_local], row_ptr_start, row_ptr_end, in_a_loc + n_a[i_id_level_local - 1], F_false_C_true, 1, diag_minus_one[i_id_level_local]);
										}

										//residualq(Amat, 1, n_a[i_id_level_local], error_approx_coarse[i_id_level_local - 1], residual_coarse[i_id_level_local - 1], row_ptr_start, row_ptr_end, in_a_loc + n_a[i_id_level_local-1], residual_fine[i_id_level_local]);
										residualq2(Amat, 1, n_a[i_id_level_local], error_approx_coarse[i_id_level_local - 1], residual_coarse[i_id_level_local - 1], row_ptr_start, row_ptr_end, in_a_loc + n_a[i_id_level_local - 1], residual_fine[i_id_level_local], diag[i_id_level_local], diag_minus_one[i_id_level_local]);


										Rnext_21[i_id_level_local] = static_cast<doublerealT>(norma(residual_fine[i_id_level_local], n_a[i_id_level_local]));
										// this is process flow logic
										if (Rnext_21[i_id_level_local] < process_flow_alpha * R0_21[i_id_level_local]) {
											// Смысл модификации в том что мы экономим итерации на пресмутере.
											break; // досрочно опускаемся на следующий уровень если он есть конечно.
										}
										else {
											Rprev_21[i_id_level_local] = Rnext_21[i_id_level_local];
										}
									}
								}

								if (iter == nu2) {
#if doubleintprecision == 1
									printf("level %lld limit postsmother iteration is reached\n", i_id_level_local);
#else
									printf("level %d limit postsmother iteration is reached\n", i_id_level_local);
#endif

								}

							}
							else {
								integer nu2_count = nu2;
								if (i_id_level_local == ilevel - 1) {
									// на самом грубом уровне нам требуется точное решение.
									//nu2_count = 100;																																					
								}
								for (integer iter = 0; iter < nu2_count; ++iter) {
									if (bonly_serial) {

										bool bflag56 = false;
										if (my_amg_manager.iFinnest_ilu == 1) {
											if (my_amg_manager.b_ilu_smoothers_in_nnz_n_LE_6) {
												doublerealT dn = static_cast<doublerealT>(1.0 * n_a[i_id_level_local]);
												doublerealT dnnz = static_cast<doublerealT>(1.0 * nnz_a[i_id_level_local]);
												if (dnnz / dn <= dapply_ilu_max_pattern_size) {
													bflag56 = true;
												}
											}
										}

										if (bILU2smoother == 1) {
											seidelq(Amat, 1, n_a[i_id_level_local], error_approx_coarse[i_id_level_local - 1], residual_coarse[i_id_level_local - 1], row_ptr_start, row_ptr_end, in_a_loc + n_a[i_id_level_local - 1], F_false_C_true, 1, diag_minus_one[i_id_level_local], iVar, i_id_level_local, ibsp_length);
											residualq2(Amat, 1, n_a[i_id_level_local], error_approx_coarse[i_id_level_local - 1], residual_coarse[i_id_level_local - 1], row_ptr_start, row_ptr_end, in_a_loc + n_a[i_id_level_local - 1], residual_fine[i_id_level_local], diag[i_id_level_local], diag_minus_one[i_id_level_local]);
#pragma omp parallel for
											for (integer i = 1; i <= n_a[i_id_level_local]; ++i) {
												milu0[i_id_level_local].zbuf[i] = residual_fine[i_id_level_local][i];
											}
											lusol_1patchforRUMBA(n_a[i_id_level_local], milu0[i_id_level_local].zbuf, milu0[i_id_level_local].zbuf2, milu0[i_id_level_local]);
#pragma omp parallel for
											for (integer i = 1; i <= n_a[i_id_level_local]; ++i) {
												error_approx_coarse[i_id_level_local - 1][i] += milu0[i_id_level_local].zbuf2[i];
											}
										}
										else if (1 && ((bILU2smoother == 2) || bflag56)) {
											seidelq(Amat, 1, n_a[i_id_level_local], error_approx_coarse[i_id_level_local - 1], residual_coarse[i_id_level_local - 1], row_ptr_start, row_ptr_end, in_a_loc + n_a[i_id_level_local - 1], F_false_C_true, 1, diag_minus_one[i_id_level_local], iVar,  i_id_level_local, ibsp_length);
											residualq2(Amat, 1, n_a[i_id_level_local], error_approx_coarse[i_id_level_local - 1], residual_coarse[i_id_level_local - 1], row_ptr_start, row_ptr_end, in_a_loc + n_a[i_id_level_local - 1], residual_fine[i_id_level_local], diag[i_id_level_local], diag_minus_one[i_id_level_local]);
#pragma omp parallel for
											for (integer i = 1; i <= n_a[i_id_level_local]; ++i) {
												milu2[i_id_level_local].zbuf[i] = residual_fine[i_id_level_local][i];
											}
											lusol_1patchforRUMBA(n_a[i_id_level_local], milu2[i_id_level_local].zbuf, milu2[i_id_level_local].zbuf2, milu2[i_id_level_local]);
#pragma omp parallel for
											for (integer i = 1; i <= n_a[i_id_level_local]; ++i) {
												error_approx_coarse[i_id_level_local - 1][i] += milu2[i_id_level_local].zbuf2[i];
											}
										}
										else {
											seidelq(Amat, 1, n_a[i_id_level_local], error_approx_coarse[i_id_level_local - 1], residual_coarse[i_id_level_local - 1], row_ptr_start, row_ptr_end, in_a_loc + n_a[i_id_level_local - 1], F_false_C_true, 1, diag_minus_one[i_id_level_local], iVar,  i_id_level_local, ibsp_length);
										}
									}
									else {
										seidelq(Amat, 1, n_a[i_id_level_local], error_approx_coarse[i_id_level_local - 1], residual_coarse[i_id_level_local - 1], nested_desection[i_id_level_local], row_ptr_start, row_ptr_end, in_a_loc + n_a[i_id_level_local - 1], F_false_C_true, 1, diag_minus_one[i_id_level_local]);
									}
								}
							}


						}

						// prolongation
						// residual_r
						//doublerealT *error_approx_fine[i_id_level_local - 1] = new doublerealT[n_a[i_id_level_local - 1] + 1];
#pragma omp parallel for
						for (integer ii = 1; ii <= n_a[i_id_level_local - 1]; ++ii) {
							error_approx_fine[i_id_level_local - 1][ii] = 0.0;
						}

						integer innz_aRP_loc = 1;
						for (integer i_72 = 0; i_72 < i_id_level_local - 1; ++i_72) {
							innz_aRP_loc += nnz_aRP[i_72];
						}
						prolongation(P, innz_aRP_loc, innz_aRP_loc - 1 + nnz_aRP[i_id_level_local - 1], error_approx_fine[i_id_level_local - 1], error_approx_coarse[i_id_level_local - 1], n_a[i_id_level_local - 1]);



						// correction
#pragma omp parallel for
						for (integer ii = 1; ii <= n_a[i_id_level_local - 1]; ++ii) {
							error_approx_coarse[i_id_level_local - 2][ii] += error_approx_fine[i_id_level_local - 1][ii];
						}

						// free
						//delete[] error_approx_fine[i_id_level_local - 1];
						//delete[] error_approx_coarse[i_id_level_local - 1];
						//delete[] residual_coarse[i_id_level_local - 1];
						//delete[] residual_fine[i_id_level_local - 1];

					} // 21
				}






				// post smothing
				if (process_flow_logic) {


					// smother
					integer iter = 0;
					for (iter = 0; iter < nu2; ++iter) {
						//quick seidel
						if (bonly_serial) {
							seidelq(Amat, 1, n_a[1], error_approx_coarse[0], residual_coarse[0], row_ptr_start, row_ptr_end, n_a[0], F_false_C_true, 1, diag_minus_one[1], iVar,  1, ibsp_length);
						}
						else {
							seidelq(Amat, 1, n_a[1], error_approx_coarse[0], residual_coarse[0], nested_desection[1], row_ptr_start, row_ptr_end, n_a[0], F_false_C_true, 1, diag_minus_one[1]);
						}

						//residualq(Amat, 1, n_a[1], error_approx_coarse[0], residual_coarse[0], row_ptr_start, row_ptr_end, n_a[0] , residual_fine[1]);
						residualq2(Amat, 1, n_a[1], error_approx_coarse[0], residual_coarse[0], row_ptr_start, row_ptr_end, n_a[0], residual_fine[1], diag[1], diag_minus_one[1]);


						Rnext_1 = static_cast<doublerealT>(norma(residual_fine[1], n_a[1]));
						// this is process flow logic
						if (Rnext_1 < process_flow_alpha * R0_1) {
							// Смысл модификации в том что мы экономим итерации на пресмутере.
							break; // досрочно опускаемся на следующий уровень если он есть конечно.
						}
						else {
							Rprev_1 = Rnext_1;
						}
					}

					if (iter == nu2) {
						printf("level 1 limit postsmother iteration is reached\n");
					}

				}
				else {
					for (integer iter = 0; iter < nu2; ++iter) {
						if (bonly_serial) {
							if (bILU2smoother == 1) {
								seidelq(Amat, 1, n_a[1], error_approx_coarse[0], residual_coarse[0], row_ptr_start, row_ptr_end, n_a[0], F_false_C_true, 1, diag_minus_one[1], iVar, 1, ibsp_length);
								residualq2(Amat, 1, n_a[1], error_approx_coarse[0], residual_coarse[0], row_ptr_start, row_ptr_end, n_a[0], residual_fine[1], diag[1], diag_minus_one[1]);
#pragma omp parallel for
								for (integer i = 1; i <= n_a[1]; ++i) {
									milu0[1].zbuf[i] = residual_fine[1][i];
								}
								lusol_1patchforRUMBA(n_a[1], milu0[1].zbuf, milu0[1].zbuf2, milu0[1]);
#pragma omp parallel for
								for (integer i = 1; i <= n_a[1]; ++i) {
									error_approx_coarse[0][i] += milu0[1].zbuf2[i];
								}
							}
							else if (1 && bILU2smoother == 2) {
								seidelq(Amat, 1, n_a[1], error_approx_coarse[0], residual_coarse[0], row_ptr_start, row_ptr_end, n_a[0], F_false_C_true, 1, diag_minus_one[1], iVar, 1, ibsp_length);
								residualq2(Amat, 1, n_a[1], error_approx_coarse[0], residual_coarse[0], row_ptr_start, row_ptr_end, n_a[0], residual_fine[1], diag[1], diag_minus_one[1]);
#pragma omp parallel for
								for (integer i = 1; i <= n_a[1]; ++i) {
									milu2[1].zbuf[i] = residual_fine[1][i];
								}
								lusol_1patchforRUMBA(n_a[1], milu2[1].zbuf, milu2[1].zbuf2, milu2[1]);
#pragma omp parallel for
								for (integer i = 1; i <= n_a[1]; ++i) {
									error_approx_coarse[0][i] += milu2[1].zbuf2[i];
								}
							}
							else {
								seidelq(Amat, 1, n_a[1], error_approx_coarse[0], residual_coarse[0], row_ptr_start, row_ptr_end, n_a[0], F_false_C_true, 1, diag_minus_one[1], iVar, 1, ibsp_length);
							}
						}
						else {
							seidelq(Amat, 1, n_a[1], error_approx_coarse[0], residual_coarse[0], nested_desection[1], row_ptr_start, row_ptr_end, n_a[0], F_false_C_true, 1, diag_minus_one[1]);
						}
					}
				}

				delete[] R0_21;
				delete[] Rprev_21;
				delete[] Rnext_21;

				R0_21 = nullptr;
				Rprev_21 = nullptr;
				Rnext_21 = nullptr;

			}

			move_up(nu1, nu2);

			// prolongation
			// residual_r
			//doublerealT *error_approx_fine[0] = new doublerealT[n_a[0] + 1];
#pragma omp parallel for
			for (integer ii = 1; ii <= n_a[0]; ++ii) {
				error_approx_fine[0][ii] = 0.0;
			}

			prolongation(P, 1, nnz_aRP[0], error_approx_fine[0], error_approx_coarse[0], n_a[0]);

			// correction
#pragma omp parallel for
			for (integer ii = 1; ii <= n_a[0]; ++ii) {
				//if (row_ptr_start[ii] != row_ptr_end[ii]) {
					// Не условие Дирихле.
				z76[ii] += error_approx_fine[0][ii];
				//	}
			}

			// free
			//delete[] error_approx_fine[0];
			//delete[] error_approx_coarse[0];
			//delete[] residual_coarse[0];
			//delete[] residual_fine[0];
		}



		//doublerealT R0_0 = 0.0;
		//doublerealT Rprev_0 = 0.0, Rnext_0 = 0.0;
		// post smothing
		if (process_flow_logic) {
			// calculate initial residual.
			//residualq<doublereal>(Amat, 1, n_a[0], z76, s76, row_ptr_start, row_ptr_end, 0, residual_fine[0]);
			residualq2(Amat, 1, n_a[0], z76, s76, row_ptr_start, row_ptr_end, 0, residual_fine[0], diag[0], diag_minus_one[0]);
			Rprev_0 = static_cast<doublerealT>(norma(residual_fine[0], n_a[0]));

			// smother
			integer iter = 0;
			for (iter = 0; iter < nFinestSweeps; ++iter) {
				//quick seidel
				if (bonly_serial) {
					seidelq<doublereal>(Amat, 1, n_a[0], z76, s76, row_ptr_start, row_ptr_end, 0, F_false_C_true, 1, diag_minus_one[0], iVar,  0, ibsp_length);
				}
				else {
					seidelq<doublereal>(Amat, 1, n_a[0], z76, s76, nested_desection[0], row_ptr_start, row_ptr_end, 0, F_false_C_true, 1, diag_minus_one[0]);
				}
				//residualq<doublereal>(Amat, 1, n_a[0], z76, s76, row_ptr_start, row_ptr_end, 0, residual_fine[0]);
				residualq2(Amat, 1, n_a[0], z76, s76, row_ptr_start, row_ptr_end, 0, residual_fine[0], diag[0], diag_minus_one[0]);
				Rnext_0 = static_cast<doublerealT>(norma(residual_fine[0], n_a[0]));
				// this is process flow logic
				if (Rnext_0 < process_flow_alpha * Rprev_0) {
					// Смысл модификации в том что мы экономим итерации на пресмутере.
					break; // досрочно опускаемся на следующий уровень если он есть конечно.
				}
				else {
					Rprev_0 = Rnext_0;
				}
			}
			if (iter == nFinestSweeps) {
				printf("level 0 limit postsmother iteration is reached\n");
			}

		}
		else {
			// nFinnestSweeps new logic 14 jan 2016.
			// smother
			for (integer iter = 0; iter < nFinestSweeps; ++iter) {
				//seidel<doublereal>(Amat, 1, nnz_a[0], z76, s76, flag, n_a[0]);
				//quick seidel
				if (bonly_serial) {
					if (bILU2smoother == 1) {
						seidelq<doublereal>(Amat, 1, n_a[0], z76, s76, row_ptr_start, row_ptr_end, 0, F_false_C_true, 1, diag_minus_one[0], iVar, 0, ibsp_length);
						residualq2(Amat, 1, n_a[0], z76, s76, row_ptr_start, row_ptr_end, 0, residual_fine[0], diag[0], diag_minus_one[0]);
#pragma omp parallel for
						for (integer i = 1; i <= n_a[0]; ++i) {
							milu0[0].zbuf[i] = residual_fine[0][i];
						}
						lusol_1patchforRUMBA(n_a[0], milu0[0].zbuf, milu0[0].zbuf2, milu0[0]);
#pragma omp parallel for
						for (integer i = 1; i <= n_a[0]; ++i) {
							z76[i] += milu0[0].zbuf2[i];
						}
					}
					else if ((bILU2smoother == 2) || (my_amg_manager.iFinnest_ilu == 1)) {
						seidelq<doublereal>(Amat, 1, n_a[0], z76, s76, row_ptr_start, row_ptr_end, 0, F_false_C_true, 1, diag_minus_one[0], iVar, 0, ibsp_length);
						residualq2(Amat, 1, n_a[0], z76, s76, row_ptr_start, row_ptr_end, 0, residual_fine[0], diag[0], diag_minus_one[0]);
#pragma omp parallel for
						for (integer i = 1; i <= n_a[0]; ++i) {
							milu2[0].zbuf[i] = residual_fine[0][i];
						}
						lusol_1patchforRUMBA(n_a[0], milu2[0].zbuf, milu2[0].zbuf2, milu2[0]);
#pragma omp parallel for
						for (integer i = 1; i <= n_a[0]; ++i) {
							z76[i] += milu2[0].zbuf2[i];
						}
					}
					else {
						seidelq<doublereal>(Amat, 1, n_a[0], z76, s76, row_ptr_start, row_ptr_end, 0, F_false_C_true, 1, diag_minus_one[0], iVar,  0, ibsp_length);
					}
				}
				else {
					seidelq<doublereal>(Amat, 1, n_a[0], z76, s76, nested_desection[0], row_ptr_start, row_ptr_end, 0, F_false_C_true, 1, diag_minus_one[0]);
				}
			}
		}

	}
} // V_cycle_solve

#endif /*BASIC_FUNCTIONS_MY_AGREGAT_AMG_CPP*/

