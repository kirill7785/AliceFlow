// Сходимость гораздо хуже, считает гораздо медленее. Надежды не оправдались.
// seidel 62
// P_Rouch1 107
// P_Rouch2 71
// P_Rouch3  замедляется катострофически сильно.




// Саусвел. 7.01.2022
void P_Rouch1(doublereal* val, integer* col_ind, integer* row_ptr, doublereal* b, doublereal* x, int n)
{
	doublereal* r = new doublereal[n];

	// Целочисленная очередь по приоритетам.
	PQ<doublereal>  pq(n, n);

	for (integer i = 0; i < n; ++i) {

		r[i] = b[i];
		for (integer j = row_ptr[i]; j < row_ptr[i + 1]; ++j) {
			r[i] -= x[col_ind[j]] * val[j];
		}

		// push r[i]
		pq.insert(fabs(r[i]), i);

	}

	for (integer it = 0; it < 1; ++it) {
		for (integer i = 0; i < n; ++i) {

			integer k = pq.readkeymaxelm();

			pq.remove(k);

			doublereal s = b[k];
			doublereal diag=0.0;
			for (integer j = row_ptr[k]; j < row_ptr[k + 1]; ++j) {
				if (col_ind[j] == k) {
					diag = val[j];
				}
				else {
					s -= x[col_ind[j]] * val[j];
				}
			}
			if (diag > 1.0e-30) {
				x[k] = s / diag;
			}
			else {
				std::cout << "zero diagonal " << k << std::endl;
				for (integer j = row_ptr[k]; j < row_ptr[k + 1]; ++j) {
					std::cout <<"val=" << val[j] << " col_ind=" << col_ind[j] << std::endl;
				}
				system("pause");
			}
			//x[k] += 0.8 * ((s) / diag - x[k]);
			doublereal omega = 1.0; // initialize.
			if (0)
			{
				// За подробностями смотри книгу Патрика Роуча стр. 183.
				doublereal rn = static_cast<doublereal>(1.0 * n);
				//optimal_omega(rn, omega); //28.07.2016

				doublereal rarg = powf(static_cast<float>(rn), 0.333333f);
				doublereal ksi = cosf(static_cast<float>(M_PI) / (static_cast<float>(rarg)));
				ksi *= ksi; // pow(ksi,2.0);

				omega = static_cast<doublereal>(2.0 * (1.0 - sqrtf(static_cast<float>(1.0 - ksi))) / ksi);// лучший выбор.

				/*if ((rn > 3.5e6) ) {
					omega = static_cast<doublereal>(1.0 + 0.08 * (omega - 1.0)); // optimum zero init 0.08
				}
				else if ((rn > 1.8e6) && (rn < 3.5e6)) {
					omega = static_cast<doublereal>(1.0 + 0.21 * (omega - 1.0)); // optimum zero init 0.21
				}
				else {
					omega = static_cast<doublereal>(1.0 + my_amg_manager.gold_const * (omega - 1.0)); // optimum zero init 0.75
				}*/

				omega = static_cast<doublereal>(1.0 + 0.14 * (omega - 1.0));
			}
			//x[k] += omega * ((s) / diag - x[k]);

			doublereal r0 = b[k];
			for (integer j = row_ptr[k]; j < row_ptr[k + 1]; ++j) {
				r0 -= x[col_ind[j]] * val[j];
			}

			// push r
			pq.insert(fabs(r0), k);

		}
	}

	delete[] r;
	r = nullptr;

}

// Саусвел. 7.01.2022
void P_Rouch2(doublereal* val, integer* col_ind, integer* row_ptr, doublereal* b, doublereal* x, int n)
{
	doublereal* r = new doublereal[n];
	doublereal* v = new doublereal[n];

	// Целочисленная очередь по приоритетам.
	PQ<doublereal>  pq(n, n);

	for (integer i = 0; i < n; ++i) {

		r[i] = b[i];
		for (integer j = row_ptr[i]; j < row_ptr[i + 1]; ++j) {
			r[i] -= x[col_ind[j]] * val[j];
		}

		// push r[i]
		//pq.insert(fabs(r[i]), i);

	}

	for (integer i = 0; i < n; ++i) {
		v[i] = 0.0;
		for (integer j = row_ptr[i]; j < row_ptr[i + 1]; ++j) {
			if (col_ind[j] != i) {
				v[i] += fabs(r[i]- r[col_ind[j]]);
			}
		}

		// push v[i]
		pq.insert(fabs(v[i]), i);
	}


	for (integer it = 0; it < 1; ++it) {
		for (integer i = 0; i < n; ++i) {

			integer k = pq.readkeymaxelm();

			pq.remove(k);

			doublereal s = b[k];
			doublereal diag=0.0;
			for (integer j = row_ptr[k]; j < row_ptr[k + 1]; ++j) {
				if (col_ind[j] == k) {
					diag = val[j];
				}
				else {
					s -= x[col_ind[j]] * val[j];
				}
			}
			if (diag > 1.0e-30) {
				x[k] = s / diag;
			}
			else {
				std::cout << "zero diagonal " << k << std::endl;
				for (integer j = row_ptr[k]; j < row_ptr[k + 1]; ++j) {
					std::cout << "val=" << val[j] << " col_ind=" << col_ind[j] << std::endl;
				}
				system("pause");
			}
			//x[k] += 0.8 * ((s) / diag - x[k]);
			doublereal omega = 1.0; // initialize.
			if (0)
			{
				// За подробностями смотри книгу Патрика Роуча стр. 183.
				doublereal rn = static_cast<doublereal>(1.0 * n);
				//optimal_omega(rn, omega); //28.07.2016

				doublereal rarg = powf(static_cast<float>(rn), 0.333333f);
				doublereal ksi = cosf(static_cast<float>(M_PI) / (static_cast<float>(rarg)));
				ksi *= ksi; // pow(ksi,2.0);

				omega = static_cast<doublereal>(2.0 * (1.0 - sqrtf(static_cast<float>(1.0 - ksi))) / ksi);// лучший выбор.

				/*if ((rn > 3.5e6) ) {
					omega = static_cast<doublereal>(1.0 + 0.08 * (omega - 1.0)); // optimum zero init 0.08
				}
				else if ((rn > 1.8e6) && (rn < 3.5e6)) {
					omega = static_cast<doublereal>(1.0 + 0.21 * (omega - 1.0)); // optimum zero init 0.21
				}
				else {
					omega = static_cast<doublereal>(1.0 + my_amg_manager.gold_const * (omega - 1.0)); // optimum zero init 0.75
				}*/

				omega = static_cast<doublereal>(1.0 + 0.14 * (omega - 1.0));
			}
			//x[k] += omega * ((s) / diag - x[k]);

			doublereal r0 = b[k];
			for (integer j = row_ptr[k]; j < row_ptr[k + 1]; ++j) {
				r0 -= x[col_ind[j]] * val[j];
			}

			r[k] = r0;


			v[k] = 0.0;
			for (integer j = row_ptr[k]; j < row_ptr[k + 1]; ++j) {
				if (col_ind[j] != k) {
					v[k] += fabs(r[k] - r[col_ind[j]]);

					// Модифицировать вариации всех соседей.


				}

				
			}

			// push v[k]
			pq.insert(fabs(v[k]), k);

			// push r
			//pq.insert(fabs(r), k);

		}
	}

	delete[] r;
	r = nullptr;

	delete[] v;
	v = nullptr;

}

// Саусвел. 7.01.2022
void P_Rouch3(doublereal* val, integer* col_ind, integer* row_ptr, doublereal* b, doublereal* x, int n)
{
	doublereal* r = new doublereal[n];
	doublereal* v = new doublereal[n];

	// Целочисленная очередь по приоритетам.
	PQ<doublereal>  pq(n, n);

	for (integer i = 0; i < n; ++i) {

		r[i] = b[i];
		for (integer j = row_ptr[i]; j < row_ptr[i + 1]; ++j) {
			r[i] -= x[col_ind[j]] * val[j];
		}

		// push r[i]
		//pq.insert(fabs(r[i]), i);

	}

	for (integer i = 0; i < n; ++i) {
		v[i] = 0.0;
		for (integer j = row_ptr[i]; j < row_ptr[i + 1]; ++j) {
			if (col_ind[j] != i) {
				v[i] += fabs(r[i] - r[col_ind[j]]);
			}
		}

		// push v[i]
		pq.insert(fabs(v[i]), i);
	}


	for (integer it = 0; it < 1; ++it) {
		for (integer i = 0; i < n; ++i) {

			integer k = pq.readkeymaxelm();

			pq.remove(k);

			doublereal s = b[k];
			doublereal diag=0.0;
			for (integer j = row_ptr[k]; j < row_ptr[k + 1]; ++j) {
				if (col_ind[j] == k) {
					diag = val[j];
				}
				else {
					s -= x[col_ind[j]] * val[j];
				}
			}
			if (diag > 1.0e-30) {
				x[k] = s / diag;
			}
			else {
				std::cout << "zero diagonal " << k << std::endl;
				for (integer j = row_ptr[k]; j < row_ptr[k + 1]; ++j) {
					std::cout << "val=" << val[j] << " col_ind=" << col_ind[j] << std::endl;
				}
				system("pause");
			}
			//x[k] += 0.8 * ((s) / diag - x[k]);
			doublereal omega = 1.0; // initialize.
			if (0)
			{
				// За подробностями смотри книгу Патрика Роуча стр. 183.
				doublereal rn = static_cast<doublereal>(1.0 * n);
				//optimal_omega(rn, omega); //28.07.2016

				doublereal rarg = powf(static_cast<float>(rn), 0.333333f);
				doublereal ksi = cosf(static_cast<float>(M_PI) / (static_cast<float>(rarg)));
				ksi *= ksi; // pow(ksi,2.0);

				omega = static_cast<doublereal>(2.0 * (1.0 - sqrtf(static_cast<float>(1.0 - ksi))) / ksi);// лучший выбор.

				/*if ((rn > 3.5e6) ) {
					omega = static_cast<doublereal>(1.0 + 0.08 * (omega - 1.0)); // optimum zero init 0.08
				}
				else if ((rn > 1.8e6) && (rn < 3.5e6)) {
					omega = static_cast<doublereal>(1.0 + 0.21 * (omega - 1.0)); // optimum zero init 0.21
				}
				else {
					omega = static_cast<doublereal>(1.0 + my_amg_manager.gold_const * (omega - 1.0)); // optimum zero init 0.75
				}*/

				omega = static_cast<doublereal>(1.0 + 0.14 * (omega - 1.0));
			}
			//x[k] += omega * ((s) / diag - x[k]);

			doublereal r0 = b[k];
			for (integer j = row_ptr[k]; j < row_ptr[k + 1]; ++j) {
				r0 -= x[col_ind[j]] * val[j];
			}

			r[k] = r0;


			v[k] = 0.0;
			for (integer j = row_ptr[k]; j < row_ptr[k + 1]; ++j) {
				if (col_ind[j] != k) {
					v[k] += fabs(r[k] - r[col_ind[j]]);

					// Модифицировать вариации всех соседей.


				}
			}

			// push v[k]
			pq.insert(fabs(v[k]), k);

			for (integer j = row_ptr[k]; j < row_ptr[k + 1]; ++j) {
				if (col_ind[j] != k) {

					integer z = col_ind[j];

					v[z] = 0.0;
					for (integer p = row_ptr[z]; p < row_ptr[z + 1]; ++p) {
						if (col_ind[p] != z) {
							v[z] += fabs(r[z] - r[col_ind[p]]);
						}
					}

					pq.remove(z);
					// push v[z]
					pq.insert(fabs(v[z]), z);

					// Модифицировать вариации всех соседей.


				}


			}

			// push r
			//pq.insert(fabs(r), k);

		}
	}

	delete[] r;
	r = nullptr;

	delete[] v;
	v = nullptr;

}