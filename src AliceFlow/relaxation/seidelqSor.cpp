

#pragma once
#ifndef BASIC_FUNCTIONS_MY_AGREGAT_AMG_RELAXATION_SEIDELQSOR_CPP
#define BASIC_FUNCTIONS_MY_AGREGAT_AMG_RELAXATION_SEIDELQSOR_CPP 1


// smoother.
  // 5 ������ 2016 � �������������� ������� �� ����� ������� �����.
  // 9 september 2015.
  // q - quick.
template <typename doublerealT>
void seidelqsor2(Ak2& Amat, integer istartq, integer iendq,
	doublerealT*& x, doublerealT*& b,
	integer*& row_ptr_start, integer*& row_ptr_end, integer iadd)
{
	// istart - ��������� ������� ��������� ��������� � ������� �.
	// iend - �������� ������� ��������� ��������� � ������� �.
	// sor 1.855 ������������ ������ �� ����� ������������ ��������������� ��������.
	// ��������� ������ ����������.
	// ������������ ����� � ��� ������ ����������. 0.8

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

							 // �� ������������� ������ ����� ������� ����� ���. 183.
	doublerealT rn = (doublerealT)(iendq - istartq + 1);
	optimal_omega(rn, omega); //28.07.2016
							  //omega = 0.7;

							  //if (isorintmemo == iadd) {
							  // ��� ����� �� ������ ���
							  //bfirst = false;
							  //}
	integer startpos = istartq + iadd;
	integer endpos = iendq + iadd;

	if (omega < 1.0) {
		if (bfirst_jacoby_start) {
			x_jacoby_buffer = new doublereal[3 * (endpos - startpos + 1)];
			i_x_jacoby_buffer_pool_size = 3 * (endpos - startpos + 1);
			bfirst_jacoby_start = false;
		}
		else {
			// ������������� ����������� ������ � ������ nu1==0.
			if (i_x_jacoby_buffer_pool_size < 3 * (endpos - startpos + 1)) {
				if (x_jacoby_buffer != nullptr) {
					delete[] x_jacoby_buffer;
					x_jacoby_buffer = nullptr;
					x_jacoby_buffer = new doublereal[3 * (endpos - startpos + 1)];
					i_x_jacoby_buffer_pool_size = 3 * (endpos - startpos + 1);
					bfirst_jacoby_start = false;
				}
			}
		}
		// copy

		if (x_jacoby_buffer == nullptr) {
			printf("ERROR: x_jacoby_buffer == nullptr.\n");
			system("PAUSE");
			exit(1);
		}

		//#pragma loop(hint_parallel(8))
#pragma omp parallel for
		for (integer ii = startpos; ii <= endpos; ++ii) {
			integer istr = ii - iadd;
			x_jacoby_buffer[istr] = x[istr];
		}


		//#pragma loop(hint_parallel(8))
#pragma omp parallel for
		for (integer ii = startpos; ii <= endpos; ++ii) {
			integer istr = ii - iadd;
			doublerealT rold = x_jacoby_buffer[istr];

			// 13.07.2016
			doublerealT ap_now = 1.0 / Amat.aij[row_ptr_start[ii]];

			x[istr] = b[istr];

			doublerealT rsum = 0.0;
			integer is1 = row_ptr_start[ii] + 1;
			integer is2 = row_ptr_end[ii];
			// ����������������� �������� �������� ����� ������.
			//#pragma omp parallel for reduction(+:rsum)
			for (integer ii1 = is1; ii1 <= is2; ++ii1)
			{
				//x[istr] += -Amat.aij[ii1]*x_jacoby_buffer[Amat.j[ii1]];
				integer ipos = Amat.j[ii1];
				// 13.07.2016
				// ������������� positive connections.
				//if ((Amat.aij[ii1] < 0.0)) {
				rsum += -Amat.aij[ii1] * x_jacoby_buffer[ipos];
				//}
				//else {
				// �� �������.
				//	ap_now += Amat.aij[ii1];
				//}
			}
			x[istr] += rsum;
			//x[istr] *= Amat.aij[row_ptr_start[ii]];
			// 13.07.2016
			x[istr] /= ap_now;

			// ����������� ������ ����� ������� ���� ����� ��������� ������
			// �.�. ��������� �������� �� �������� ����������.
			x[istr] = omega * x[istr] + (1.0 - omega) * rold; // this is SOR
		}
	}
	else {

		// 3 ������ 2016. ������������ ����� ������-�������.
		if (isimmetricGS_switch == 0) {
			// 3 ������ 2016 ���� ���������������� �������� �� BSKDmitrii ��� ������������ ����� ������ - ������� ����������.
			// ������ ������������ ������������ ����� ������ -�������. ����������� ������� ����� �������.


#pragma omp parallel for
			for (integer ii = startpos; ii <= endpos; ++ii) {
				integer istr = ii - iadd;
				doublerealT rold = x[istr];

				// 13.07.2016
				doublerealT ap_now = 1.0 / Amat.aij[row_ptr_start[ii]];
				// 28.01.2017
				//ap_now = 0.0;

				x[istr] = b[istr];

				doublerealT rsum = 0.0;
				integer is1 = row_ptr_start[ii] + 1;
				integer is2 = row_ptr_end[ii];
				// ����������������� �������� �������� ����� ������.
				//#pragma omp parallel for reduction(+:rsum)
				for (integer ii1 = is1; ii1 <= is2; ++ii1)
				{
					//x[istr] += -Amat.aij[ii1]*x[Amat.j[ii1]];
					integer ipos = Amat.j[ii1];
					// 13.07.2016
					if (1) {
						// ������������� positive connections.
						//if ((Amat.aij[ii1] < 0.0)) {
						// ����� �� � ���� � ������� � ������� ������������.
						rsum += -Amat.aij[ii1] * x[ipos];
						//}
						//else {
						// �� ��������.
						//	ap_now += Amat.aij[ii1];
						//}
					}
					else {
						// ����� � ����������. �� �������� ����� �������: TKM, TKM1, TKM2.

						// ������������� positive connections.
						if ((Amat.aij[ii1] < 0.0)) {
							rsum += -Amat.aij[ii1] * x[ipos];
							//ap_now += fabs(Amat.aij[ii1]);
						}
						else {
							// �� ��������.
							// �������� ��-�� ���� ��� ��� ������� ������.
							//ap_now += fabs(Amat.aij[ii1]);
							if (fabs(x[ipos]) > fabs(x[istr])) {
								ap_now += fabs(Amat.aij[ii1]);
							}
							else rsum += -Amat.aij[ii1] * x[ipos];
						}
					}
				}
				x[istr] += rsum;
				//x[istr] *= Amat.aij[row_ptr_start[ii]];
				// 13.07.2016
				x[istr] /= ap_now;

				// ����������� ������ ����� ������� ���� ����� ��������� ������
				// �.�. ��������� �������� �� �������� ����������.
				x[istr] = omega * x[istr] + (1.0 - omega) * rold; // this is SOR
			}
			// �� � ���� ������ �� �����������. 3 ������ 2016.
			//isimmetricGS_switch = 1;

		}
		else {

#pragma omp parallel for
			for (integer ii = endpos; ii >= startpos; ii--) {
				integer istr = ii - iadd;
				doublerealT rold = x[istr];

				// 13.07.2016
				doublerealT ap_now = 1.0 / Amat.aij[row_ptr_start[ii]];

				x[istr] = b[istr];

				doublerealT rsum = 0.0;
				integer is1 = row_ptr_start[ii] + 1;
				integer is2 = row_ptr_end[ii];
				// ����������������� �������� �������� ����� ������.
				//#pragma omp parallel for reduction(+:rsum)
				for (integer ii1 = is1; ii1 <= is2; ++ii1)
				{
					//x[istr] += -Amat.aij[ii1]*x[Amat.j[ii1]];
					integer ipos = Amat.j[ii1];
					// 13.07.2016
					// ������������� positive connections.
					//if ((Amat.aij[ii1] < 0.0)) {
					rsum += -Amat.aij[ii1] * x[ipos];
					//}
					//else {
					// �� �������.
					//	ap_now += Amat.aij[ii1];
					//}
				}
				x[istr] += rsum;
				//x[istr] *= Amat.aij[row_ptr_start[ii]];
				// 13.07.2016
				x[istr] /= ap_now;

				// ����������� ������ ����� ������� ���� ����� ��������� ������
				// �.�. ��������� �������� �� �������� ����������.
				x[istr] = omega * x[istr] + (1.0 - omega) * rold; // this is SOR
			}
			//isimmetricGS_switch = 0;
		}
	}


} // seidelqsor2

// smoother.
  // 5 ������ 2016 � �������������� ������� �� ����� ������� �����.
  // 9 september 2015.
  // q - quick.
template <typename doublerealT>
void seidelqsor2(Ak2& Amat, integer istartq, integer iendq,
	doublerealT*& x, doublerealT*& b,
	integer*& row_ptr_start, integer*& row_ptr_end,
	integer iadd, doublerealT*& diag_minus_one)
{
	// istart - ��������� ������� ��������� ��������� � ������� �.
	// iend - �������� ������� ��������� ��������� � ������� �.
	// sor 1.855 ������������ ������ �� ����� ������������ ��������������� ��������.
	// ��������� ������ ����������.
	// ������������ ����� � ��� ������ ����������. 0.8

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

							 // �� ������������� ������ ����� ������� ����� ���. 183.
	doublerealT rn = (doublerealT)(iendq - istartq + 1);
	optimal_omega(rn, omega); //28.07.2016
							  //omega = 0.7;

							  //if (isorintmemo == iadd) {
							  // ��� ����� �� ������ ���
							  //bfirst = false;
							  //}
	integer startpos = istartq + iadd;
	integer endpos = iendq + iadd;

	if (omega < 1.0) {
		if (bfirst_jacoby_start) {
			x_jacoby_buffer = new doublereal[3 * (endpos - startpos + 1)];
			i_x_jacoby_buffer_pool_size = 3 * (endpos - startpos + 1);
			bfirst_jacoby_start = false;
		}
		else {
			// ������������� ����������� ������ � ������ nu1==0.
			if (i_x_jacoby_buffer_pool_size < 3 * (endpos - startpos + 1)) {
				if (x_jacoby_buffer != nullptr) {
					delete[] x_jacoby_buffer;
					x_jacoby_buffer = nullptr;
					x_jacoby_buffer = new doublereal[3 * (endpos - startpos + 1)];
					i_x_jacoby_buffer_pool_size = 3 * (endpos - startpos + 1);
					bfirst_jacoby_start = false;
				}
			}
		}
		// copy

		if (x_jacoby_buffer == nullptr) {
			printf("ERROR: x_jacoby_buffer == nullptr.\n");
			system("PAUSE");
			exit(1);
		}

		//#pragma loop(hint_parallel(8))
#pragma omp parallel for
		for (integer ii = startpos; ii <= endpos; ++ii) {
			integer istr = ii - iadd;
			x_jacoby_buffer[istr] = x[istr];
		}


		//#pragma loop(hint_parallel(8))
#pragma omp parallel for
		for (integer ii = startpos; ii <= endpos; ++ii) {
			integer istr = ii - iadd;
			doublerealT rold = x_jacoby_buffer[istr];

			// 13.07.2016
			doublerealT ap_now = 1.0 / Amat.aij[row_ptr_start[ii]];

			x[istr] = diag_minus_one[istr] * b[istr];

			doublerealT rsum = 0.0;
			integer is1 = row_ptr_start[ii] + 1;
			integer is2 = row_ptr_end[ii];
			// ����������������� �������� �������� ����� ������.
			//#pragma omp parallel for reduction(+:rsum)
			for (integer ii1 = is1; ii1 <= is2; ++ii1)
			{
				//x[istr] += -Amat.aij[ii1]*x_jacoby_buffer[Amat.j[ii1]];
				integer ipos = Amat.j[ii1];
				// 13.07.2016
				// ������������� positive connections.
				//if ((Amat.aij[ii1] < 0.0)) {
				rsum += -Amat.aij[ii1] * x_jacoby_buffer[ipos];
				//}
				//else {
				// �� �������.
				//	ap_now += Amat.aij[ii1];
				//}
			}
			x[istr] += rsum;
			//x[istr] *= Amat.aij[row_ptr_start[ii]];
			// 13.07.2016
			x[istr] /= ap_now;

			// ����������� ������ ����� ������� ���� ����� ��������� ������
			// �.�. ��������� �������� �� �������� ����������.
			x[istr] = omega * x[istr] + (1.0 - omega) * rold; // this is SOR
		}
	}
	else {

		// 3 ������ 2016. ������������ ����� ������-�������.
		if (isimmetricGS_switch == 0) {
			// 3 ������ 2016 ���� ���������������� �������� �� BSKDmitrii ��� ������������ ����� ������ - ������� ����������.
			// ������ ������������ ������������ ����� ������ -�������. ����������� ������� ����� �������.


#pragma omp parallel for
			for (integer ii = startpos; ii <= endpos; ++ii) {
				integer istr = ii - iadd;
				doublerealT rold = x[istr];

				// 13.07.2016
				doublerealT ap_now = 1.0 / Amat.aij[row_ptr_start[ii]];
				// 28.01.2017
				//ap_now = 0.0;

				x[istr] = diag_minus_one[istr] * b[istr];

				doublerealT rsum = 0.0;
				integer is1 = row_ptr_start[ii] + 1;
				integer is2 = row_ptr_end[ii];
				// ����������������� �������� �������� ����� ������.
				//#pragma omp parallel for reduction(+:rsum)
				for (integer ii1 = is1; ii1 <= is2; ++ii1)
				{
					//x[istr] += -Amat.aij[ii1]*x[Amat.j[ii1]];
					integer ipos = Amat.j[ii1];
					// 13.07.2016
					if (1) {
						// ������������� positive connections.
						//if ((Amat.aij[ii1] < 0.0)) {
						// ����� �� � ���� � ������� � ������� ������������.
						rsum += -Amat.aij[ii1] * x[ipos];
						//}
						//else {
						// �� ��������.
						//	ap_now += Amat.aij[ii1];
						//}
					}
					else {
						// ����� � ����������. �� �������� ����� �������: TKM, TKM1, TKM2.

						// ������������� positive connections.
						if ((Amat.aij[ii1] < 0.0)) {
							rsum += -Amat.aij[ii1] * x[ipos];
							//ap_now += fabs(Amat.aij[ii1]);
						}
						else {
							// �� ��������.
							// �������� ��-�� ���� ��� ��� ������� ������.
							//ap_now += fabs(Amat.aij[ii1]);
							if (fabs(x[ipos]) > fabs(x[istr])) {
								ap_now += fabs(Amat.aij[ii1]);
							}
							else rsum += -Amat.aij[ii1] * x[ipos];
						}
					}
				}
				x[istr] += rsum;
				//x[istr] *= Amat.aij[row_ptr_start[ii]];
				// 13.07.2016
				x[istr] /= ap_now;

				// ����������� ������ ����� ������� ���� ����� ��������� ������
				// �.�. ��������� �������� �� �������� ����������.
				x[istr] = omega * x[istr] + (1.0 - omega) * rold; // this is SOR
			}
			// �� � ���� ������ �� �����������. 3 ������ 2016.
			//isimmetricGS_switch = 1;

		}
		else {

#pragma omp parallel for
			for (integer ii = endpos; ii >= startpos; ii--) {
				integer istr = ii - iadd;
				doublerealT rold = x[istr];

				// 13.07.2016
				doublerealT ap_now = 1.0 / Amat.aij[row_ptr_start[ii]];

				x[istr] = diag_minus_one[istr] * b[istr];

				doublerealT rsum = 0.0;
				integer is1 = row_ptr_start[ii] + 1;
				integer is2 = row_ptr_end[ii];
				// ����������������� �������� �������� ����� ������.
				//#pragma omp parallel for reduction(+:rsum)
				for (integer ii1 = is1; ii1 <= is2; ++ii1)
				{
					//x[istr] += -Amat.aij[ii1]*x[Amat.j[ii1]];
					integer ipos = Amat.j[ii1];
					// 13.07.2016
					// ������������� positive connections.
					//if ((Amat.aij[ii1] < 0.0)) {
					rsum += -Amat.aij[ii1] * x[ipos];
					//}
					//else {
					// �� �������.
					//	ap_now += Amat.aij[ii1];
					//}
				}
				x[istr] += rsum;
				//x[istr] *= Amat.aij[row_ptr_start[ii]];
				// 13.07.2016
				x[istr] /= ap_now;

				// ����������� ������ ����� ������� ���� ����� ��������� ������
				// �.�. ��������� �������� �� �������� ����������.
				x[istr] = omega * x[istr] + (1.0 - omega) * rold; // this is SOR
			}
			//isimmetricGS_switch = 0;
		}
	}


} // seidelqsor2


// smoother.
  // 14 ������ 2015 ������ ��� �������������� ����� ����������� ������������.
  // �������� ������ � ���� �������: nFinestSweeps=2, nPreSweeps=0, nPostSweeps=2.
  // ����� ����������� ��������� ����������� ���� �� ��������������.
  // 5 ������ 2016 � �������������� ������� �� ����� ������� �����.
  // 9 september 2015.
  // q - quick.
template <typename doublerealT>
void seidelqsor3(Ak2& Amat, integer istartq, integer iendq, doublerealT*& x, doublerealT*& b, integer*& row_ptr_start, integer*& row_ptr_end, integer iadd)
{
	// istart - ��������� ������� ��������� ��������� � ������� �.
	// iend - �������� ������� ��������� ��������� � ������� �.
	// sor 1.855 ������������ ������ �� ����� ������������ ��������������� ��������.
	// ��������� ������ ����������.
	// ������������ ����� � ��� ������ ����������. 0.8

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

	// �� ������������� ������ ����� ������� ����� ���. 183.
	doublerealT rn = (doublerealT)(iendq - istartq + 1);
	optimal_omega(rn, omega);

	//if (isorintmemo == iadd) {
	// ��� ����� �� ������ ���
	//bfirst = false;
	//}
	integer startpos = istartq + iadd;
	integer endpos = iendq + iadd;
	if (bswitch_direct_seidelqsor3) {
		for (integer ii = startpos; ii <= endpos; ++ii) {
			integer istr = ii - iadd;
			doublerealT rold = x[istr];

			x[istr] = b[istr];

			doublerealT rsum = 0.0;
			integer is1 = row_ptr_start[ii] + 1;
			integer is2 = row_ptr_end[ii];
			// ����������������� �������� �������� ����� ������.
			//#pragma omp parallel for reduction(+:rsum)
			for (integer ii1 = is1; ii1 <= is2; ++ii1)
			{
				//x[istr] += -Amat.aij[ii1]*x[Amat.j[ii1]];
				integer ipos = Amat.j[ii1];
				rsum += -Amat.aij[ii1] * x[ipos];
			}
			x[istr] += rsum;
			x[istr] *= Amat.aij[row_ptr_start[ii]];

			// ����������� ������ ����� ������� ���� ����� ��������� ������
			// �.�. ��������� �������� �� �������� ����������.
			x[istr] = omega * x[istr] + (1.0 - omega) * rold; // this is SOR
		}
	}
	else {
		// ������ ����������� ������������.

		for (integer ii = endpos; ii >= startpos; ii--) {
			integer istr = ii - iadd;
			doublerealT rold = x[istr];

			x[istr] = b[istr];

			doublerealT rsum = 0.0;
			integer is1 = row_ptr_start[ii] + 1;
			integer is2 = row_ptr_end[ii];
			// ����������������� �������� �������� ����� ������.
			//#pragma omp parallel for reduction(+:rsum)
			for (integer ii1 = is1; ii1 <= is2; ++ii1)
			{
				//x[istr] += -Amat.aij[ii1]*x[Amat.j[ii1]];
				integer ipos = Amat.j[ii1];
				rsum += -Amat.aij[ii1] * x[ipos];
			}
			x[istr] += rsum;
			x[istr] *= Amat.aij[row_ptr_start[ii]];

			// ����������� ������ ����� ������� ���� ����� ��������� ������
			// �.�. ��������� �������� �� �������� ����������.
			x[istr] = omega * x[istr] + (1.0 - omega) * rold; // this is SOR
		}
	}

	// ����� ����������� ������������ �� ��������.
	bswitch_direct_seidelqsor3 = !bswitch_direct_seidelqsor3;

} // seidelqsor3




// smoother.
  // 16 ������ 2016 ����������������� �� ����������� ����������.
  // 5 ������ 2016 � �������������� ������� �� ����� ������� �����.
  // 9 september 2015.
  // q - quick.
template <typename doublerealT>
void seidelqsor2Pcpu(Ak2& Amat, integer istartq, integer iendq, doublerealT*& x, doublerealT*& b, integer*& row_ptr_start, integer*& row_ptr_end, integer iadd)
{
	// istart - ��������� ������� ��������� ��������� � ������� �.
	// iend - �������� ������� ��������� ��������� � ������� �.
	// sor 1.855 ������������ ������ �� ����� ������������ ��������������� ��������.
	// ��������� ������ ����������.
	// ������������ ����� � ��� ������ ����������. 0.8

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

							 // �� ������������� ������ ����� ������� ����� ���. 183.
	doublerealT rn = (doublerealT)(iendq - istartq + 1);
	optimal_omega(rn, omega);

	//if (isorintmemo == iadd) {
	// ��� ����� �� ������ ���
	//bfirst = false;
	//}
	const integer startpos = istartq + iadd;
	const integer endpos = iendq + iadd;

	const integer inumcore_loc = 2;

if (inumcore_loc == 1)
	{

		for (integer ii = startpos; ii <= endpos; ++ii) {
			integer istr = ii - iadd;
			doublerealT rold = x[istr];

			x[istr] = b[istr];

			doublerealT rsum = 0.0;
			integer is1 = row_ptr_start[ii] + 1;
			integer is2 = row_ptr_end[ii];
			// ����������������� �������� �������� ����� ������.
			//#pragma omp parallel for reduction(+:rsum)
			for (integer ii1 = is1; ii1 <= is2; ++ii1)
			{
				//x[istr] += -Amat.aij[ii1]*x[Amat.j[ii1]];
				integer ipos = Amat.j[ii1];
				rsum += -Amat.aij[ii1] * x[ipos];
			}
			x[istr] += rsum;
			x[istr] *= Amat.aij[row_ptr_start[ii]];

			// ����������� ������ ����� ������� ���� ����� ��������� ������
			// �.�. ��������� �������� �� �������� ����������.
			x[istr] = omega * x[istr] + (1.0 - omega) * rold; // this is SOR
		}

	}


if (inumcore_loc == 2)
	{

		// ����� ����������� ���������� � nested_desection
		// ������� ���� ������� �������� �� ��������� �������������.

		integer middle = (startpos + endpos) / 2;
		integer iadd1 = iadd;
		integer iadd2 = iadd;
		integer middle1 = middle;

#pragma omp parallel sections
		{
#pragma omp section
			{
				// ������ ������� ������.
				for (integer ii = startpos; ii <= middle; ++ii) {
					integer istr = ii - iadd1;
					doublerealT rold = x[istr];

					doublerealT x1buf = 0.0;
					x1buf = b[istr];

					doublerealT rsum1 = 0.0;
					integer is1 = row_ptr_start[ii] + 1;
					integer is2 = row_ptr_end[ii];
					// ����������������� �������� �������� ����� ������.
					//#pragma omp parallel for reduction(+:rsum)
					for (integer ii1 = is1; ii1 <= is2; ++ii1)
					{
						//x1buf += -Amat.aij[ii1]*x[Amat.j[ii1]];
						integer ipos = Amat.j[ii1];
						rsum1 = rsum1 - Amat.aij[ii1] * x[ipos];
					}
					x1buf = x1buf + rsum1;
					x1buf = x1buf * Amat.aij[row_ptr_start[ii]];

					// ����������� ������ ����� ������� ���� ����� ��������� ������
					// �.�. ��������� �������� �� �������� ����������.
					x1buf = omega * x1buf + (1.0 - omega) * rold; // this is SOR
					x[istr] = x1buf;
				}
			}

#pragma omp section
			{
				// ������ ������� ������. 
				for (integer ii_1 = middle1 + 1; ii_1 <= endpos; ii_1++) {
					integer istr1 = ii_1 - iadd2;
					doublerealT rold1 = x[istr1];

					doublerealT x2buf = 0.0;
					x2buf = b[istr1];

					doublerealT rsum2 = 0.0;
					integer is3 = row_ptr_start[ii_1] + 1;
					integer is4 = row_ptr_end[ii_1];
					// ����������������� �������� �������� ����� ������.
					//#pragma omp parallel for reduction(+:rsum)
					for (integer ii2 = is3; ii2 <= is4; ii2++)
					{
						//x[istr1] += -Amat.aij[ii2]*x[Amat.j[ii2]];
						integer ipos = Amat.j[ii2];
						rsum2 = rsum2 - Amat.aij[ii2] * x[ipos];
					}
					x2buf = x2buf + rsum2;
					x2buf = x2buf * Amat.aij[row_ptr_start[ii_1]];

					// ����������� ������ ����� ������� ���� ����� ��������� ������
					// �.�. ��������� �������� �� �������� ����������.
					x2buf = omega * x2buf + (1.0 - omega) * rold1; // this is SOR
					x[istr1] = x2buf;
				}
			}
		}

	}


} // seidelqsor2Pcpu



  // smoother.
  // 16 ������ 2016 ����������������� �� ����������� ����������.
  // 5 ������ 2016 � �������������� ������� �� ����� ������� �����.
  // 9 september 2015.
  // q - quick.
template <typename doublerealT>
void seidelqsor2Pcpu(Ak2& Amat, integer istartq, integer iendq, doublerealT*& x, doublerealT*& b, bool*& bnested_desection, integer*& row_ptr_start, integer*& row_ptr_end, integer iadd)
{
	// istart - ��������� ������� ��������� ��������� � ������� �.
	// iend - �������� ������� ��������� ��������� � ������� �.
	// sor 1.855 ������������ ������ �� ����� ������������ ��������������� ��������.
	// ��������� ������ ����������.
	// ������������ ����� � ��� ������ ����������. 0.8

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

	// �� ������������� ������ ����� ������� ����� ���. 183.
	doublerealT rn = (doublerealT)(iendq - istartq + 1);
	optimal_omega(rn, omega);

	//if (isorintmemo == iadd) {
	// ��� ����� �� ������ ���
	//bfirst = false;
	//}
	const integer startpos = istartq + iadd;
	const integer endpos = iendq + iadd;

	const integer inumcore_loc = 1;

if (inumcore_loc == 1)
	{

		// ������������ ������� ���������.
		for (integer ii = startpos; ii <= endpos; ++ii) {
			integer istr = ii - iadd;
			doublerealT rold = x[istr];

			x[istr] = b[istr];

			doublerealT rsum = 0.0;
			integer is1 = row_ptr_start[ii] + 1;
			integer is2 = row_ptr_end[ii];
			// ����������������� �������� �������� ����� ������.
			//#pragma omp parallel for reduction(+:rsum)
			for (integer ii1 = is1; ii1 <= is2; ++ii1)
			{
				//x[istr] += -Amat.aij[ii1]*x[Amat.j[ii1]];
				integer ipos = Amat.j[ii1];
				rsum += -Amat.aij[ii1] * x[ipos];
			}
			x[istr] += rsum;
			x[istr] *= Amat.aij[row_ptr_start[ii]];

			// ����������� ������ ����� ������� ���� ����� ��������� ������
			// �.�. ��������� �������� �� �������� ����������.
			x[istr] = omega * x[istr] + (1.0 - omega) * rold; // this is SOR
		}

	}


if (inumcore_loc == 2) 
	{

		integer middle = static_cast<integer>(0.5 * (startpos + endpos));
		doublerealT omega1 = omega;
		doublerealT omega2 = omega;
		integer iadd1 = iadd;
		integer iadd2 = iadd;
		integer middle1 = middle;
		doublerealT ic_all = 0.0;
		doublerealT ic_separator = 0.0;
		for (integer i_51 = startpos - iadd1; i_51 <= endpos - iadd1; i_51++) {
			bnested_desection_global_amg[i_51] = bnested_desection[i_51];
			ic_all += 1.0;
			if (!bnested_desection_global_amg[i_51]) ic_separator += 1.0;
		}
		//printf("all=%e separator=%e\n",ic_all,ic_separator);
		//system("pause");

		if (2.0 * ic_separator < ic_all) {
			// ����������� ������ ��� ������� �������.



			//#pragma omp parallel sections num_threads(4)
			//{

			//printf_s("Hello from thread %d\n", omp_get_thread_num());
			//#pragma omp section
			//	printf_s("Hello from thread %d\n", omp_get_thread_num());
			///}

			// ������ � nesteddesection ��������� �� ��������� ������������� � ��
			// �������� �������� ��������� ������������������ ��� ���������� ����� ����.

			//default(shared)

#pragma omp parallel  shared(bnested_desection_global_amg, bnested_desection,x,row_ptr_start,row_ptr_end,b,Amat)
			{
#pragma omp	sections
				{
#pragma omp section
					{
						//#pragma omp parallel sections num_threads(2)
						//	{

						// ������ ������� ������.
						for (integer ii = startpos; ii <= middle; ++ii) {
							integer istr = ii - iadd1;
							if (bnested_desection_global_amg[istr]) {
								doublerealT rold = x[istr];

								doublerealT x1buf = 0.0;
								x1buf = b[istr];

								doublerealT rsum1 = 0.0;
								integer is1 = row_ptr_start[ii] + 1;
								integer is2 = row_ptr_end[ii];
								// ����������������� �������� �������� ����� ������.
								//#pragma omp parallel for reduction(+:rsum)
								for (integer ii1 = is1; ii1 <= is2; ++ii1)
								{
									//x1buf += -Amat.aij[ii1]*x[Amat.j[ii1]];
									integer ipos = Amat.j[ii1];
									rsum1 = rsum1 - Amat.aij[ii1] * x[ipos];
								}
								x1buf = x1buf + rsum1;
								x1buf = x1buf * Amat.aij[row_ptr_start[ii]];

								// ����������� ������ ����� ������� ���� ����� ��������� ������
								// �.�. ��������� �������� �� �������� ����������.
								x1buf = omega1 * x1buf + (1.0 - omega1) * rold; // this is SOR
								x[istr] = x1buf;
							}
						}
					}

#pragma omp section
					{
						// ������ ������� ������. 
						for (integer ii_1 = middle1 + 1; ii_1 <= endpos; ii_1++) {
							integer istr1 = ii_1 - iadd2;
							if (bnested_desection[istr1]) {
								doublerealT rold1 = x[istr1];

								doublerealT x2buf = 0.0;
								x2buf = b[istr1];

								doublerealT rsum2 = 0.0;
								integer is3 = row_ptr_start[ii_1] + 1;
								integer is4 = row_ptr_end[ii_1];
								// ����������������� �������� �������� ����� ������.
								//#pragma omp parallel for reduction(+:rsum)
								for (integer ii2 = is3; ii2 <= is4; ii2++)
								{
									//x[istr1] += -Amat.aij[ii2]*x[Amat.j[ii2]];
									integer ipos = Amat.j[ii2];
									rsum2 = rsum2 - Amat.aij[ii2] * x[ipos];
								}
								x2buf = x2buf + rsum2;
								x2buf = x2buf * Amat.aij[row_ptr_start[ii_1]];

								// ����������� ������ ����� ������� ���� ����� ��������� ������
								// �.�. ��������� �������� �� �������� ����������.
								x2buf = omega2 * x2buf + (1.0 - omega2) * rold1; // this is SOR
								x[istr1] = x2buf;
							}
						}
					}
				}
			}

			// ������������ ��������� �����.
			for (integer ii = startpos; ii <= endpos; ++ii) {
				integer istr = ii - iadd;
				if (!bnested_desection[istr]) {
					doublerealT rold = x[istr];

					x[istr] = b[istr];

					doublerealT rsum = 0.0;
					integer is1 = row_ptr_start[ii] + 1;
					integer is2 = row_ptr_end[ii];
					// ����������������� �������� �������� ����� ������.
					//#pragma omp parallel for reduction(+:rsum)
					for (integer ii1 = is1; ii1 <= is2; ++ii1)
					{
						//x[istr] += -Amat.aij[ii1]*x[Amat.j[ii1]];
						integer ipos = Amat.j[ii1];
						rsum += -Amat.aij[ii1] * x[ipos];
					}
					x[istr] += rsum;
					x[istr] *= Amat.aij[row_ptr_start[ii]];

					// ����������� ������ ����� ������� ���� ����� ��������� ������
					// �.�. ��������� �������� �� �������� ����������.
					x[istr] = omega * x[istr] + (1.0 - omega) * rold; // this is SOR
				}
			}
		}
		else {
			// �������
			// ������������ ������� ���������.
			for (integer ii = startpos; ii <= endpos; ++ii) {
				integer istr = ii - iadd;
				doublerealT rold = x[istr];

				x[istr] = b[istr];

				doublerealT rsum = 0.0;
				integer is1 = row_ptr_start[ii] + 1;
				integer is2 = row_ptr_end[ii];
				// ����������������� �������� �������� ����� ������.
				//#pragma omp parallel for reduction(+:rsum)
				for (integer ii1 = is1; ii1 <= is2; ++ii1)
				{
					//x[istr] += -Amat.aij[ii1]*x[Amat.j[ii1]];
					integer ipos = Amat.j[ii1];
					rsum += -Amat.aij[ii1] * x[ipos];
				}
				x[istr] += rsum;
				x[istr] *= Amat.aij[row_ptr_start[ii]];

				// ����������� ������ ����� ������� ���� ����� ��������� ������
				// �.�. ��������� �������� �� �������� ����������.
				x[istr] = omega * x[istr] + (1.0 - omega) * rold; // this is SOR
			}
		}

	}


} // seidelqsor2Pcpu+nested desection


// smoother.
  // 16 ������ 2016 ����������������� �� ����������� ����������.
  // 5 ������ 2016 � �������������� ������� �� ����� ������� �����.
  // 9 september 2015.
  // q - quick.
template <typename doublerealT>
void seidelqsor2Pcpu(Ak2& Amat, integer istartq, integer iendq, doublerealT*& x, doublerealT*& b, integer*& row_ptr_start, integer*& row_ptr_end, integer iadd, doublerealT*& diag_minus_one)
{
	// istart - ��������� ������� ��������� ��������� � ������� �.
	// iend - �������� ������� ��������� ��������� � ������� �.
	// sor 1.855 ������������ ������ �� ����� ������������ ��������������� ��������.
	// ��������� ������ ����������.
	// ������������ ����� � ��� ������ ����������. 0.8

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

							 // �� ������������� ������ ����� ������� ����� ���. 183.
	doublerealT rn = (doublerealT)(iendq - istartq + 1);
	optimal_omega(rn, omega);

	//if (isorintmemo == iadd) {
	// ��� ����� �� ������ ���
	//bfirst = false;
	//}
	const integer startpos = istartq + iadd;
	const integer endpos = iendq + iadd;

	const integer inumcore_loc = 2;

if (inumcore_loc == 1)
	{

		for (integer ii = startpos; ii <= endpos; ++ii) {
			integer istr = ii - iadd;
			doublerealT rold = x[istr];

			x[istr] = diag_minus_one[istr] * b[istr];

			doublerealT rsum = 0.0;
			integer is1 = row_ptr_start[ii] + 1;
			integer is2 = row_ptr_end[ii];
			// ����������������� �������� �������� ����� ������.
			//#pragma omp parallel for reduction(+:rsum)
			for (integer ii1 = is1; ii1 <= is2; ++ii1)
			{
				//x[istr] += -Amat.aij[ii1]*x[Amat.j[ii1]];
				integer ipos = Amat.j[ii1];
				rsum += -Amat.aij[ii1] * x[ipos];
			}
			x[istr] += rsum;
			x[istr] *= Amat.aij[row_ptr_start[ii]];

			// ����������� ������ ����� ������� ���� ����� ��������� ������
			// �.�. ��������� �������� �� �������� ����������.
			x[istr] = omega * x[istr] + (1.0 - omega) * rold; // this is SOR
		}

	}


if (inumcore_loc == 2) 
	{

		// ����� ����������� ���������� � nested_desection
		// ������� ���� ������� �������� �� ��������� �������������.

		integer middle = (startpos + endpos) / 2;
		integer iadd1 = iadd;
		integer iadd2 = iadd;
		integer middle1 = middle;

#pragma omp parallel sections
		{
#pragma omp section
			{
				// ������ ������� ������.
				for (integer ii = startpos; ii <= middle; ++ii) {
					integer istr = ii - iadd1;
					doublerealT rold = x[istr];

					doublerealT x1buf = 0.0;
					x1buf = diag_minus_one[istr] * b[istr];

					doublerealT rsum1 = 0.0;
					integer is1 = row_ptr_start[ii] + 1;
					integer is2 = row_ptr_end[ii];
					// ����������������� �������� �������� ����� ������.
					//#pragma omp parallel for reduction(+:rsum)
					for (integer ii1 = is1; ii1 <= is2; ++ii1)
					{
						//x1buf += -Amat.aij[ii1]*x[Amat.j[ii1]];
						integer ipos = Amat.j[ii1];
						rsum1 = rsum1 - Amat.aij[ii1] * x[ipos];
					}
					x1buf = x1buf + rsum1;
					x1buf = x1buf * Amat.aij[row_ptr_start[ii]];

					// ����������� ������ ����� ������� ���� ����� ��������� ������
					// �.�. ��������� �������� �� �������� ����������.
					x1buf = omega * x1buf + (1.0 - omega) * rold; // this is SOR
					x[istr] = x1buf;
				}
			}

#pragma omp section
			{
				// ������ ������� ������. 
				for (integer ii_1 = middle1 + 1; ii_1 <= endpos; ii_1++) {
					integer istr1 = ii_1 - iadd2;
					doublerealT rold1 = x[istr1];

					doublerealT x2buf = 0.0;
					x2buf = diag_minus_one[istr1] * b[istr1];

					doublerealT rsum2 = 0.0;
					integer is3 = row_ptr_start[ii_1] + 1;
					integer is4 = row_ptr_end[ii_1];
					// ����������������� �������� �������� ����� ������.
					//#pragma omp parallel for reduction(+:rsum)
					for (integer ii2 = is3; ii2 <= is4; ii2++)
					{
						//x[istr1] += -Amat.aij[ii2]*x[Amat.j[ii2]];
						integer ipos = Amat.j[ii2];
						rsum2 = rsum2 - Amat.aij[ii2] * x[ipos];
					}
					x2buf = x2buf + rsum2;
					x2buf = x2buf * Amat.aij[row_ptr_start[ii_1]];

					// ����������� ������ ����� ������� ���� ����� ��������� ������
					// �.�. ��������� �������� �� �������� ����������.
					x2buf = omega * x2buf + (1.0 - omega) * rold1; // this is SOR
					x[istr1] = x2buf;
				}
			}
		}

	}


} // seidelqsor2Pcpu





// smoother.
 // 16 ������ 2016 ����������������� �� ����������� ����������.
 // 5 ������ 2016 � �������������� ������� �� ����� ������� �����.
 // 9 september 2015.
 // q - quick.
template <typename doublerealT>
void seidelqsor2Pcpu(Ak2& Amat, integer istartq, integer iendq, doublerealT*& x, doublerealT*& b, bool*& bnested_desection, integer*& row_ptr_start, integer*& row_ptr_end, integer iadd, doublerealT*& diag_minus_one)
{
	// istart - ��������� ������� ��������� ��������� � ������� �.
	// iend - �������� ������� ��������� ��������� � ������� �.
	// sor 1.855 ������������ ������ �� ����� ������������ ��������������� ��������.
	// ��������� ������ ����������.
	// ������������ ����� � ��� ������ ����������. 0.8

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

	// �� ������������� ������ ����� ������� ����� ���. 183.
	doublerealT rn = (doublerealT)(iendq - istartq + 1);
	optimal_omega(rn, omega);

	//if (isorintmemo == iadd) {
	// ��� ����� �� ������ ���
	//bfirst = false;
	//}
	const integer startpos = istartq + iadd;
	const integer endpos = iendq + iadd;

	const integer inumcore_loc = 1;

	if (inumcore_loc == 1) 
	{

		// ������������ ������� ���������.
		for (integer ii = startpos; ii <= endpos; ++ii) {
			integer istr = ii - iadd;
			doublerealT rold = x[istr];

			x[istr] = diag_minus_one[istr] * b[istr];

			doublerealT rsum = 0.0;
			integer is1 = row_ptr_start[ii] + 1;
			integer is2 = row_ptr_end[ii];
			// ����������������� �������� �������� ����� ������.
			//#pragma omp parallel for reduction(+:rsum)
			for (integer ii1 = is1; ii1 <= is2; ++ii1)
			{
				//x[istr] += -Amat.aij[ii1]*x[Amat.j[ii1]];
				integer ipos = Amat.j[ii1];
				rsum += -Amat.aij[ii1] * x[ipos];
			}
			x[istr] += rsum;
			x[istr] *= Amat.aij[row_ptr_start[ii]];

			// ����������� ������ ����� ������� ���� ����� ��������� ������
			// �.�. ��������� �������� �� �������� ����������.
			x[istr] = omega * x[istr] + (1.0 - omega) * rold; // this is SOR
		}

	}


	if (inumcore_loc == 2)
	{

		integer middle = static_cast<integer>(0.5 * (startpos + endpos));
		doublerealT omega1 = omega;
		doublerealT omega2 = omega;
		integer iadd1 = iadd;
		integer iadd2 = iadd;
		integer middle1 = middle;
		doublerealT ic_all = 0.0;
		doublerealT ic_separator = 0.0;
		for (integer i_51 = startpos - iadd1; i_51 <= endpos - iadd1; i_51++) {
			bnested_desection_global_amg[i_51] = bnested_desection[i_51];
			ic_all += 1.0;
			if (!bnested_desection_global_amg[i_51]) ic_separator += 1.0;
		}
		//printf("all=%e separator=%e\n",ic_all,ic_separator);
		//system("pause");

		if (2.0 * ic_separator < ic_all) {
			// ����������� ������ ��� ������� �������.



			//#pragma omp parallel sections num_threads(4)
			//{

			//printf_s("Hello from thread %d\n", omp_get_thread_num());
			//#pragma omp section
			//	printf_s("Hello from thread %d\n", omp_get_thread_num());
			///}

			// ������ � nesteddesection ��������� �� ��������� ������������� � ��
			// �������� �������� ��������� ������������������ ��� ���������� ����� ����.

			//default(shared)

#pragma omp parallel  shared(bnested_desection_global_amg, bnested_desection,x,row_ptr_start,row_ptr_end,b,Amat)
			{
#pragma omp	sections
				{
#pragma omp section
					{
						//#pragma omp parallel sections num_threads(2)
						//	{

						// ������ ������� ������.
						for (integer ii = startpos; ii <= middle; ++ii) {
							integer istr = ii - iadd1;
							if (bnested_desection_global_amg[istr]) {
								doublerealT rold = x[istr];

								doublerealT x1buf = 0.0;
								x1buf = diag_minus_one[istr] * b[istr];

								doublerealT rsum1 = 0.0;
								integer is1 = row_ptr_start[ii] + 1;
								integer is2 = row_ptr_end[ii];
								// ����������������� �������� �������� ����� ������.
								//#pragma omp parallel for reduction(+:rsum)
								for (integer ii1 = is1; ii1 <= is2; ++ii1)
								{
									//x1buf += -Amat.aij[ii1]*x[Amat.j[ii1]];
									integer ipos = Amat.j[ii1];
									rsum1 = rsum1 - Amat.aij[ii1] * x[ipos];
								}
								x1buf = x1buf + rsum1;
								x1buf = x1buf * Amat.aij[row_ptr_start[ii]];

								// ����������� ������ ����� ������� ���� ����� ��������� ������
								// �.�. ��������� �������� �� �������� ����������.
								x1buf = omega1 * x1buf + (1.0 - omega1) * rold; // this is SOR
								x[istr] = x1buf;
							}
						}
					}

#pragma omp section
					{
						// ������ ������� ������. 
						for (integer ii_1 = middle1 + 1; ii_1 <= endpos; ii_1++) {
							integer istr1 = ii_1 - iadd2;
							if (bnested_desection[istr1]) {
								doublerealT rold1 = x[istr1];

								doublerealT x2buf = 0.0;
								x2buf = diag_minus_one[istr1] * b[istr1];

								doublerealT rsum2 = 0.0;
								integer is3 = row_ptr_start[ii_1] + 1;
								integer is4 = row_ptr_end[ii_1];
								// ����������������� �������� �������� ����� ������.
								//#pragma omp parallel for reduction(+:rsum)
								for (integer ii2 = is3; ii2 <= is4; ii2++)
								{
									//x[istr1] += -Amat.aij[ii2]*x[Amat.j[ii2]];
									integer ipos = Amat.j[ii2];
									rsum2 = rsum2 - Amat.aij[ii2] * x[ipos];
								}
								x2buf = x2buf + rsum2;
								x2buf = x2buf * Amat.aij[row_ptr_start[ii_1]];

								// ����������� ������ ����� ������� ���� ����� ��������� ������
								// �.�. ��������� �������� �� �������� ����������.
								x2buf = omega2 * x2buf + (1.0 - omega2) * rold1; // this is SOR
								x[istr1] = x2buf;
							}
						}
					}
				}
			}

			// ������������ ��������� �����.
			for (integer ii = startpos; ii <= endpos; ++ii) {
				integer istr = ii - iadd;
				if (!bnested_desection[istr]) {
					doublerealT rold = x[istr];

					x[istr] = diag_minus_one[istr] * b[istr];

					doublerealT rsum = 0.0;
					integer is1 = row_ptr_start[ii] + 1;
					integer is2 = row_ptr_end[ii];
					// ����������������� �������� �������� ����� ������.
					//#pragma omp parallel for reduction(+:rsum)
					for (integer ii1 = is1; ii1 <= is2; ++ii1)
					{
						//x[istr] += -Amat.aij[ii1]*x[Amat.j[ii1]];
						integer ipos = Amat.j[ii1];
						rsum += -Amat.aij[ii1] * x[ipos];
					}
					x[istr] += rsum;
					x[istr] *= Amat.aij[row_ptr_start[ii]];

					// ����������� ������ ����� ������� ���� ����� ��������� ������
					// �.�. ��������� �������� �� �������� ����������.
					x[istr] = omega * x[istr] + (1.0 - omega) * rold; // this is SOR
				}
			}
		}
		else {
			// �������
			// ������������ ������� ���������.
			for (integer ii = startpos; ii <= endpos; ++ii) {
				integer istr = ii - iadd;
				doublerealT rold = x[istr];

				x[istr] = diag_minus_one[istr] * b[istr];

				doublerealT rsum = 0.0;
				integer is1 = row_ptr_start[ii] + 1;
				integer is2 = row_ptr_end[ii];
				// ����������������� �������� �������� ����� ������.
				//#pragma omp parallel for reduction(+:rsum)
				for (integer ii1 = is1; ii1 <= is2; ++ii1)
				{
					//x[istr] += -Amat.aij[ii1]*x[Amat.j[ii1]];
					integer ipos = Amat.j[ii1];
					rsum += -Amat.aij[ii1] * x[ipos];
				}
				x[istr] += rsum;
				x[istr] *= Amat.aij[row_ptr_start[ii]];

				// ����������� ������ ����� ������� ���� ����� ��������� ������
				// �.�. ��������� �������� �� �������� ����������.
				x[istr] = omega * x[istr] + (1.0 - omega) * rold; // this is SOR
			}
		}

	}


} // seidelqsor2Pcpu+nested desection

// smoother.
// 1 september 2015.
template <typename doublerealT>
void seidel1(Ak2& Amat, integer istart, integer iend, doublerealT*& x, doublerealT*& b, bool*& flag, integer n)
{
	// istart - ��������� ������� ��������� ��������� � ������� �.
	// iend - �������� ������� ��������� ��������� � ������� �.
	for (integer i = 1; i <= n; ++i) {
		flag[i] = false;
	}
	for (integer ii = istart; ii <= iend; ++ii) {
		if (flag[Amat.i[ii]] == false) {
			integer istr = Amat.i[ii];
			integer ic = ii;
			doublerealT ap = 0.0;
			x[istr] = b[istr];
			while ((ic <= iend) && (Amat.i[ic] == istr)) {
				if (Amat.j[ic] != istr) {
					x[istr] += -Amat.aij[ic] * x[Amat.j[ic]];
				}
				else ap = Amat.aij[ic];
				ic++;
			}
			/*
			if (fabs(ap) < 1.0e-30) {
			
				std::cout << "zero diagonal elements in string " << istr << std::endl;

				system("pause");
				exit(1);
			}
			else */ {
				x[istr] /= ap;
			}
			flag[Amat.i[ii]] = true;
		}
	}


} // seidel1


//seidel_for_cg(Amat, row_ptr_start, row_ptr_end, z76, s76, flag_seidel, n_a[0]);
// smoother.
// 1 september 2015. 02.02.2022.
template <typename doublerealT>
void seidel_for_cg(Ak2& Amat, integer* &row_ptr_start, integer* &row_ptr_end, doublerealT*& x, doublerealT*& b, integer n)
{
	// istart - ��������� ������� ��������� ��������� � ������� �.
	// iend - �������� ������� ��������� ��������� � ������� �.
	
	//for (integer ii = istart; ii <= iend; ++ii)
	for (integer i = 1; i <= n; ++i) {
		
		
				integer istr = i;
				integer ic = row_ptr_start[i];
				doublerealT ap = 0.0;
				x[istr] = b[istr];
				while ((ic <= row_ptr_end[i])) {
					if (Amat.j[ic] != istr) {
						x[istr] += -Amat.aij[ic] * x[Amat.j[ic]];
					}
					else ap = Amat.aij[ic]; // �������� �������� ��������� ��������.
					ic++;
				}
				/*
				if (fabs(ap) < 1.0e-30) {

					std::cout << "zero diagonal elements in string " << istr << std::endl;

					system("pause");
					exit(1);
				}
				else */ {
					x[istr] *= ap;
				}
	
	}


} // seidel_for_cg

// smoother.
// 5 jan 2016.
// 1 september 2015.
template <typename doublerealT>
void seidelsor(Ak2& Amat, integer istart, integer iend, doublerealT*& x, doublerealT*& b, bool*& flag, integer n)
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

	// �� ������������� ������ ����� ������� ����� ���. 183.
	optimal_omega(rn, omega);//28.07.2016
	//omega = 0.7;

	// istart - ��������� ������� ��������� ��������� � ������� �.
	// iend - �������� ������� ��������� ��������� � ������� �.
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
					if (Amat.aij[ic] < 0.0) {
						x[istr] += -Amat.aij[ic] * x[Amat.j[ic]];
					}
					else {
						ap += Amat.aij[ic];
					}
				}
				else ap += Amat.aij[ic];
				ic++;
			}
			/*
			if (fabs(ap) < 1.0e-30) {
			
					std::cout << "zero diagonal elements in string " << istr << std::endl;

					system("pause");
					exit(1);
			}
			else */ {
				x[istr] /= ap;
				x[istr] = omega * x[istr] + (1.0 - omega) * rold; // this is SOR
			}
			flag[Amat.i[ii]] = true;
		}
	}


} // seidelsor




// smoother.
// 5 jan 2016.
// 1 september 2015.
template <typename doublerealT>
void seidel(Ak2& Amat, integer istart, integer iend, doublerealT*& x, doublerealT*& b, bool*& flag, integer n)
{
	//seidel1<doublerealT>(Amat,  istart,  iend, x, b, flag, n);
	seidelsor<doublerealT>(Amat, istart, iend, x, b, flag, n);
	//early_naive_relaxation_method<doublerealT>(Amat, istart, iend, x, b, flag, n);
}

// smoother.
// 9 september 2015.
// q - quick.
template <typename doublerealT>
void seidelqstable(Ak2& Amat, integer istartq, integer iendq, doublerealT*& x, doublerealT*& b, integer*& row_ptr_start, integer*& row_ptr_end, integer iadd)
{
	// istart - ��������� ������� ��������� ��������� � ������� �.
	// iend - �������� ������� ��������� ��������� � ������� �.
	integer startpos = istartq + iadd;
	integer endpos = iendq + iadd;
	for (integer ii = startpos; ii <= endpos; ++ii) {
		integer istr = ii - iadd;
		x[istr] = b[istr];

		for (integer ii1 = row_ptr_start[ii] + 1; ii1 <= row_ptr_end[ii]; ++ii1)
		{
			x[istr] += -Amat.aij[ii1] * x[Amat.j[ii1]];
		}
		x[istr] *= Amat.aij[row_ptr_start[ii]];
	}


} // seidelq


template <typename doublerealT>
void seidelqsor(Ak2& Amat, integer istartq, integer iendq, doublerealT*& x, doublerealT*& b, integer*& row_ptr_start, integer*& row_ptr_end, integer iadd)
{
	// istart - ��������� ������� ��������� ��������� � ������� �.
	// iend - �������� ������� ��������� ��������� � ������� �.
	// sor 1.855 ������������ ������ �� ����� ������������ ��������������� ��������.
	// ��������� ������ ����������.
	// ������������ ����� � ��� ������ ����������. 0.8

	// BSKDmitrii
	// omega   iter  time,s
	// 1.0 106 43
	// 1.1 98 42
	// 1.15 94 40 best
	// 1.2 90 40
	// 1.225 413 1min 37s
	// 1.25 divergence detected
	// 1.3 divergence detected

	doublerealT omega = 1.15; // ������ �����.
	bool bfirst = false;
	//if (isorintmemo == iadd) {
	// ��� ����� �� ������ ���
	//bfirst = false;
	//}
	integer startpos = istartq + iadd;
	integer endpos = iendq + iadd;
	for (integer ii = startpos; ii <= endpos; ++ii) {
		integer istr = ii - iadd;
		doublerealT rold = x[istr];

		x[istr] = b[istr];

		doublerealT rsum = 0.0;
		// �������������� ������-�� �������� ����� ������.
		//#pragma omp parallel for reduction(+:rsum)
		for (integer ii1 = row_ptr_start[ii] + 1; ii1 <= row_ptr_end[ii]; ++ii1)
		{
			rsum += -Amat.aij[ii1] * x[Amat.j[ii1]];
			//x[istr] += -Amat.aij[ii1]*x[Amat.j[ii1]];
		}
		x[istr] += rsum;

		x[istr] *= Amat.aij[row_ptr_start[ii]];
		if (bfirst) {
			bfirst = false;
		}
		else {
			// ����������� ������ ����� ������� ���� ����� ��������� ������
			// �.�. ��������� �������� �� �������� ����������.
			x[istr] = omega * x[istr] + (1.0 - omega) * rold; // this is SOR
		}
	}


} // seidelqsor

// smoother.
// 5 ���� 2017 ��������� CF-Jacobi smoothing (F - smoothing).
// 5 ������ 2016 � �������������� ������� �� ����� ������� �����.
// 9 september 2015.
// q - quick.
template <typename doublerealT>
void seidelqsor2(Ak2& Amat, integer istartq, integer iendq, doublerealT*& x, doublerealT*& b, integer*& row_ptr_start, integer*& row_ptr_end, integer iadd, bool*& F_false_C_true, integer idirect)
{
	// F_false_C_true - ��������� ���������� � 1.
	// idirect==0 douwn
	// idirect==1 up

	// istart - ��������� ������� ��������� ��������� � ������� �.
	// iend - �������� ������� ��������� ��������� � ������� �.
	// sor 1.855 ������������ ������ �� ����� ������������ ��������������� ��������.
	// ��������� ������ ����������.
	// ������������ ����� � ��� ������ ����������. 0.8

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

							// �� ������������� ������ ����� ������� ����� ���. 183.
	doublerealT rn = (doublerealT)(iendq - istartq + 1);
	optimal_omega(rn, omega); //28.07.2016
							  //omega = 0.7;

							  //if (isorintmemo == iadd) {
							  // ��� ����� �� ������ ���
							  //bfirst = false;
							  //}
	integer startpos = istartq + iadd;
	integer endpos = iendq + iadd;

	if (omega < 1.0) {
		if (bfirst_jacoby_start) {
			x_jacoby_buffer = new doublereal[3 * (endpos - startpos + 1)];
			i_x_jacoby_buffer_pool_size = 3 * (endpos - startpos + 1);
			bfirst_jacoby_start = false;
		}
		else {
			// ������������� ����������� ������ � ������ nu1==0.
			if (i_x_jacoby_buffer_pool_size < 3 * (endpos - startpos + 1)) {
				if (x_jacoby_buffer != nullptr) {
					delete[] x_jacoby_buffer;
					x_jacoby_buffer = nullptr;
					x_jacoby_buffer = new doublereal[3 * (endpos - startpos + 1)];
					i_x_jacoby_buffer_pool_size = 3 * (endpos - startpos + 1);
					bfirst_jacoby_start = false;
				}
			}
		}
		// copy



		if (idirect == 1) {
			// ���������� �����.

			// ������� F ����� C.

			//#pragma loop(hint_parallel(8))
#pragma omp parallel for
			for (integer ii = startpos; ii <= endpos; ++ii) {
				integer istr = ii - iadd;
				x_jacoby_buffer[istr] = x[istr];
			}

			//#pragma loop(hint_parallel(8))
#pragma omp parallel for
			for (integer ii = startpos; ii <= endpos; ++ii) {
				if ((!my_amg_manager.bcf_reorder) || (F_false_C_true[ii] == false)) { // F nodes

					integer istr = ii - iadd;
					doublerealT rold = x_jacoby_buffer[istr];

					// 13.07.2016
					doublerealT ap_now = 1.0 / Amat.aij[row_ptr_start[ii]];

					x[istr] = b[istr];

					doublerealT rsum = 0.0;
					integer is1 = row_ptr_start[ii] + 1;
					integer is2 = row_ptr_end[ii];
					// ����������������� �������� �������� ����� ������.
					//#pragma omp parallel for reduction(+:rsum)
					for (integer ii1 = is1; ii1 <= is2; ++ii1)
					{
						//x[istr] += -Amat.aij[ii1]*x_jacoby_buffer[Amat.j[ii1]];
						integer ipos = Amat.j[ii1];
						// 13.07.2016
						// ������������� positive connections.
						//if ((Amat.aij[ii1] < 0.0)) {
						rsum += -Amat.aij[ii1] * x_jacoby_buffer[ipos];
						//}
						//else {
						// �� �������.
						//	ap_now += Amat.aij[ii1];
						//}
					}
					x[istr] += rsum;
					//x[istr] *= Amat.aij[row_ptr_start[ii]];
					// 13.07.2016
					x[istr] /= ap_now;

					// ����������� ������ ����� ������� ���� ����� ��������� ������
					// �.�. ��������� �������� �� �������� ����������.
					//if (is1 <= is2) {
						// ������ ���� ��� �� ������� ������� ��������� ����������.
					x[istr] = omega * x[istr] + (1.0 - omega) * rold; // this is SOR
				//}
				}
			}

			// update.
			//#pragma loop(hint_parallel(8))
#pragma omp parallel for
			for (integer ii = startpos; ii <= endpos; ++ii) {
				integer istr = ii - iadd;
				x_jacoby_buffer[istr] = x[istr];
			}

			//#pragma loop(hint_parallel(8))
#pragma omp parallel for
			for (integer ii = startpos; ii <= endpos; ++ii) {
				if ((!my_amg_manager.bcf_reorder) || (F_false_C_true[ii])) { // C nodes

					integer istr = ii - iadd;
					doublerealT rold = x_jacoby_buffer[istr];

					// 13.07.2016
					doublerealT ap_now = 1.0 / Amat.aij[row_ptr_start[ii]];

					x[istr] = b[istr];

					doublerealT rsum = 0.0;
					integer is1 = row_ptr_start[ii] + 1;
					integer is2 = row_ptr_end[ii];
					// ����������������� �������� �������� ����� ������.
					//#pragma omp parallel for reduction(+:rsum)
					for (integer ii1 = is1; ii1 <= is2; ++ii1)
					{
						//x[istr] += -Amat.aij[ii1]*x_jacoby_buffer[Amat.j[ii1]];
						integer ipos = Amat.j[ii1];
						// 13.07.2016
						// ������������� positive connections.
						//if ((Amat.aij[ii1] < 0.0)) {
						rsum += -Amat.aij[ii1] * x_jacoby_buffer[ipos];
						//}
						//else {
						// �� �������.
						//	ap_now += Amat.aij[ii1];
						//}
					}
					x[istr] += rsum;
					//x[istr] *= Amat.aij[row_ptr_start[ii]];
					// 13.07.2016
					x[istr] /= ap_now;

					// ����������� ������ ����� ������� ���� ����� ��������� ������
					// �.�. ��������� �������� �� �������� ����������.
					//if (is1 <= is2) {
						// ������ ���� ��� �� ������� ������� ��������� ����������.
					x[istr] = omega * x[istr] + (1.0 - omega) * rold; // this is SOR
				//}
				}
			}



		}
		else {
			// idirect==0
			// ������� � ����� F.

			//#pragma loop(hint_parallel(8))
#pragma omp parallel for
			for (integer ii = startpos; ii <= endpos; ++ii) {
				integer istr = ii - iadd;
				x_jacoby_buffer[istr] = x[istr];
			}

			//#pragma loop(hint_parallel(8))
#pragma omp parallel for
			for (integer ii = startpos; ii <= endpos; ++ii) {
				if ((!my_amg_manager.bcf_reorder) || (F_false_C_true[ii])) { // C nodes

					integer istr = ii - iadd;
					doublerealT rold = x_jacoby_buffer[istr];

					// 13.07.2016
					doublerealT ap_now = 1.0 / Amat.aij[row_ptr_start[ii]];

					x[istr] = b[istr];

					doublerealT rsum = 0.0;
					integer is1 = row_ptr_start[ii] + 1;
					integer is2 = row_ptr_end[ii];
					// ����������������� �������� �������� ����� ������.
					//#pragma omp parallel for reduction(+:rsum)
					for (integer ii1 = is1; ii1 <= is2; ++ii1)
					{
						//x[istr] += -Amat.aij[ii1]*x_jacoby_buffer[Amat.j[ii1]];
						integer ipos = Amat.j[ii1];
						// 13.07.2016
						// ������������� positive connections.
						//if ((Amat.aij[ii1] < 0.0)) {
						rsum += -Amat.aij[ii1] * x_jacoby_buffer[ipos];
						//}
						//else {
						// �� �������.
						//	ap_now += Amat.aij[ii1];
						//}
					}
					x[istr] += rsum;
					//x[istr] *= Amat.aij[row_ptr_start[ii]];
					// 13.07.2016
					x[istr] /= ap_now;

					// ����������� ������ ����� ������� ���� ����� ��������� ������
					// �.�. ��������� �������� �� �������� ����������.
					//if (is1 <= is2) {
						// ������ ���� ��� �� ������� ������� ��������� ����������.
					x[istr] = omega * x[istr] + (1.0 - omega) * rold; // this is SOR
				//}
				}
			}


			//#pragma loop(hint_parallel(8))
#pragma omp parallel for
			for (integer ii = startpos; ii <= endpos; ++ii) {
				integer istr = ii - iadd;
				x_jacoby_buffer[istr] = x[istr];
			}

			//#pragma loop(hint_parallel(8))
#pragma omp parallel for
			for (integer ii = startpos; ii <= endpos; ++ii) {
				if ((!my_amg_manager.bcf_reorder) || (F_false_C_true[ii] == false)) { // F nodes

					integer istr = ii - iadd;
					doublerealT rold = x_jacoby_buffer[istr];

					// 13.07.2016
					doublerealT ap_now = 1.0 / Amat.aij[row_ptr_start[ii]];

					x[istr] = b[istr];

					doublerealT rsum = 0.0;
					integer is1 = row_ptr_start[ii] + 1;
					integer is2 = row_ptr_end[ii];
					// ����������������� �������� �������� ����� ������.
					//#pragma omp parallel for reduction(+:rsum)
					for (integer ii1 = is1; ii1 <= is2; ++ii1)
					{
						//x[istr] += -Amat.aij[ii1]*x_jacoby_buffer[Amat.j[ii1]];
						integer ipos = Amat.j[ii1];
						// 13.07.2016
						// ������������� positive connections.
						//if ((Amat.aij[ii1] < 0.0)) {
						rsum += -Amat.aij[ii1] * x_jacoby_buffer[ipos];
						//}
						//else {
						// �� �������.
						//	ap_now += Amat.aij[ii1];
						//}
					}
					x[istr] += rsum;
					//x[istr] *= Amat.aij[row_ptr_start[ii]];
					// 13.07.2016
					x[istr] /= ap_now;

					// ����������� ������ ����� ������� ���� ����� ��������� ������
					// �.�. ��������� �������� �� �������� ����������.
					//if (is1 <= is2) {
						// ������ ���� ��� �� ������� ������� ��������� ����������.
					x[istr] = omega * x[istr] + (1.0 - omega) * rold; // this is SOR
				//}
				}
			}

		}


	}
	else {

		// 3 ������ 2016. ������������ ����� ������-�������.
		if (isimmetricGS_switch == 0) {
			// 3 ������ 2016 ���� ���������������� �������� �� BSKDmitrii ��� ������������ ����� ������ - ������� ����������.
			// ������ ������������ ������������ ����� ������ -�������. ����������� ������� ����� �������.

			if (idirect == 1) {
				// ���������� �����.



				// ������� F ����� C.

//----->#pragma omp parallel for
				for (integer ii = startpos; ii <= endpos; ++ii) {
					if ((!my_amg_manager.bcf_reorder) || (F_false_C_true[ii] == false)) { // F nodes

						integer istr = ii - iadd;
						doublerealT rold = x[istr];

						// 13.07.2016
						doublerealT ap_now = 1.0 / Amat.aij[row_ptr_start[ii]];
						// 28.01.2017
						//ap_now = 0.0;

						x[istr] = b[istr];

						doublerealT rsum = 0.0;
						integer is1 = row_ptr_start[ii] + 1;
						integer is2 = row_ptr_end[ii];
						// ����������������� �������� �������� ����� ������.
						//#pragma omp parallel for reduction(+:rsum)
						for (integer ii1 = is1; ii1 <= is2; ++ii1)
						{
							//x[istr] += -Amat.aij[ii1]*x[Amat.j[ii1]];
							integer ipos = Amat.j[ii1];
							// 13.07.2016
							if (1) {
								// ������������� positive connections.
								//if ((Amat.aij[ii1] < 0.0)) {
								// ����� �� � ���� � ������� � ������� ������������.
								rsum += -Amat.aij[ii1] * x[ipos];
								//}
								//else {
								// �� ��������.
								//	ap_now += Amat.aij[ii1];
								//}
							}
							else {
								// ����� � ����������. �� �������� ����� �������: TKM, TKM1, TKM2.

								// ������������� positive connections.
								if ((Amat.aij[ii1] < 0.0)) {
									rsum += -Amat.aij[ii1] * x[ipos];
									//ap_now += fabs(Amat.aij[ii1]);
								}
								else {
									// �� ��������.
									// �������� ��-�� ���� ��� ��� ������� ������.
									//ap_now += fabs(Amat.aij[ii1]);
									if (fabs(x[ipos]) > fabs(x[istr])) {
										ap_now += fabs(Amat.aij[ii1]);
									}
									else rsum += -Amat.aij[ii1] * x[ipos];
								}
							}
						}
						x[istr] += rsum;
						//x[istr] *= Amat.aij[row_ptr_start[ii]];
						// 13.07.2016
						x[istr] /= ap_now;

						// ����������� ������ ����� ������� ���� ����� ��������� ������
						// �.�. ��������� �������� �� �������� ����������.
						//if (is1 <= is2) {
							// ������ ���� ��� �� ������� ������� ��������� ����������.
						x[istr] = omega * x[istr] + (1.0 - omega) * rold; // this is SOR
					//}
					}
				}

				//---->#pragma omp parallel for
				for (integer ii = startpos; ii <= endpos; ++ii) {
					if ((!my_amg_manager.bcf_reorder) || (F_false_C_true[ii])) { // C nodes

						integer istr = ii - iadd;
						doublerealT rold = x[istr];

						// 13.07.2016
						doublerealT ap_now = 1.0 / Amat.aij[row_ptr_start[ii]];
						// 28.01.2017
						//ap_now = 0.0;

						x[istr] = b[istr];

						doublerealT rsum = 0.0;
						integer is1 = row_ptr_start[ii] + 1;
						integer is2 = row_ptr_end[ii];
						// ����������������� �������� �������� ����� ������.
						//#pragma omp parallel for reduction(+:rsum)
						for (integer ii1 = is1; ii1 <= is2; ++ii1)
						{
							//x[istr] += -Amat.aij[ii1]*x[Amat.j[ii1]];
							integer ipos = Amat.j[ii1];
							// 13.07.2016
							if (1) {
								// ������������� positive connections.
								//if ((Amat.aij[ii1] < 0.0)) {
								// ����� �� � ���� � ������� � ������� ������������.
								rsum += -Amat.aij[ii1] * x[ipos];
								//}
								//else {
								// �� ��������.
								//	ap_now += Amat.aij[ii1];
								//}
							}
							else {
								// ����� � ����������. �� �������� ����� �������: TKM, TKM1, TKM2.

								// ������������� positive connections.
								if ((Amat.aij[ii1] < 0.0)) {
									rsum += -Amat.aij[ii1] * x[ipos];
									//ap_now += fabs(Amat.aij[ii1]);
								}
								else {
									// �� ��������.
									// �������� ��-�� ���� ��� ��� ������� ������.
									//ap_now += fabs(Amat.aij[ii1]);
									if (fabs(x[ipos]) > fabs(x[istr])) {
										ap_now += fabs(Amat.aij[ii1]);
									}
									else rsum += -Amat.aij[ii1] * x[ipos];
								}
							}
						}
						x[istr] += rsum;
						//x[istr] *= Amat.aij[row_ptr_start[ii]];
						// 13.07.2016
						x[istr] /= ap_now;

						// ����������� ������ ����� ������� ���� ����� ��������� ������
						// �.�. ��������� �������� �� �������� ����������.
						//if (is1 <= is2) {
							// ������ ���� ��� �� ������� ������� ��������� ����������.
						x[istr] = omega * x[istr] + (1.0 - omega) * rold; // this is SOR
					//}
					}
				}

			}
			else {

				// idirect==0
				// ������� � ����� F.

//---->#pragma omp parallel for
				for (integer ii = startpos; ii <= endpos; ++ii) {
					if ((!my_amg_manager.bcf_reorder) || (F_false_C_true[ii])) { // C nodes

						integer istr = ii - iadd;
						doublerealT rold = x[istr];

						// 13.07.2016
						doublerealT ap_now = 1.0 / Amat.aij[row_ptr_start[ii]];
						// 28.01.2017
						//ap_now = 0.0;

						x[istr] = b[istr];

						doublerealT rsum = 0.0;
						integer is1 = row_ptr_start[ii] + 1;
						integer is2 = row_ptr_end[ii];
						// ����������������� �������� �������� ����� ������.
						//#pragma omp parallel for reduction(+:rsum)
						for (integer ii1 = is1; ii1 <= is2; ++ii1)
						{
							//x[istr] += -Amat.aij[ii1]*x[Amat.j[ii1]];
							integer ipos = Amat.j[ii1];
							// 13.07.2016
							if (1) {
								// ������������� positive connections.
								//if ((Amat.aij[ii1] < 0.0)) {
								// ����� �� � ���� � ������� � ������� ������������.
								rsum += -Amat.aij[ii1] * x[ipos];
								//}
								//else {
								// �� ��������.
								//	ap_now += Amat.aij[ii1];
								//}
							}
							else {
								// ����� � ����������. �� �������� ����� �������: TKM, TKM1, TKM2.

								// ������������� positive connections.
								if ((Amat.aij[ii1] < 0.0)) {
									rsum += -Amat.aij[ii1] * x[ipos];
									//ap_now += fabs(Amat.aij[ii1]);
								}
								else {
									// �� ��������.
									// �������� ��-�� ���� ��� ��� ������� ������.
									//ap_now += fabs(Amat.aij[ii1]);
									if (fabs(x[ipos]) > fabs(x[istr])) {
										ap_now += fabs(Amat.aij[ii1]);
									}
									else rsum += -Amat.aij[ii1] * x[ipos];
								}
							}
						}
						x[istr] += rsum;
						//x[istr] *= Amat.aij[row_ptr_start[ii]];
						// 13.07.2016
						x[istr] /= ap_now;

						// ����������� ������ ����� ������� ���� ����� ��������� ������
						// �.�. ��������� �������� �� �������� ����������.
						//if (is1 <= is2) {
							// ������ ���� ��� �� ������� ������� ��������� ����������.
						x[istr] = omega * x[istr] + (1.0 - omega) * rold; // this is SOR
					//}
					}
				}

				//---->#pragma omp parallel for
				for (integer ii = startpos; ii <= endpos; ++ii) {
					if ((!my_amg_manager.bcf_reorder) || (F_false_C_true[ii] == false)) { // F nodes

						integer istr = ii - iadd;
						doublerealT rold = x[istr];

						// 13.07.2016
						doublerealT ap_now = 1.0 / Amat.aij[row_ptr_start[ii]];
						// 28.01.2017
						//ap_now = 0.0;

						x[istr] = b[istr];

						doublerealT rsum = 0.0;
						integer is1 = row_ptr_start[ii] + 1;
						integer is2 = row_ptr_end[ii];
						// ����������������� �������� �������� ����� ������.
						//#pragma omp parallel for reduction(+:rsum)
						for (integer ii1 = is1; ii1 <= is2; ++ii1)
						{
							//x[istr] += -Amat.aij[ii1]*x[Amat.j[ii1]];
							integer ipos = Amat.j[ii1];
							// 13.07.2016
							if (1) {
								// ������������� positive connections.
								//if ((Amat.aij[ii1] < 0.0)) {
								// ����� �� � ���� � ������� � ������� ������������.
								rsum += -Amat.aij[ii1] * x[ipos];
								//}
								//else {
								// �� ��������.
								//	ap_now += Amat.aij[ii1];
								//}
							}
							else {
								// ����� � ����������. �� �������� ����� �������: TKM, TKM1, TKM2.

								// ������������� positive connections.
								if ((Amat.aij[ii1] < 0.0)) {
									rsum += -Amat.aij[ii1] * x[ipos];
									//ap_now += fabs(Amat.aij[ii1]);
								}
								else {
									// �� ��������.
									// �������� ��-�� ���� ��� ��� ������� ������.
									//ap_now += fabs(Amat.aij[ii1]);
									if (fabs(x[ipos]) > fabs(x[istr])) {
										ap_now += fabs(Amat.aij[ii1]);
									}
									else rsum += -Amat.aij[ii1] * x[ipos];
								}
							}
						}
						x[istr] += rsum;
						//x[istr] *= Amat.aij[row_ptr_start[ii]];
						// 13.07.2016
						x[istr] /= ap_now;

						// ����������� ������ ����� ������� ���� ����� ��������� ������
						// �.�. ��������� �������� �� �������� ����������.
						//if (is1 <= is2) {
							// ������ ���� ��� �� ������� ������� ��������� ����������.
						x[istr] = omega * x[istr] + (1.0 - omega) * rold; // this is SOR
					//}
					}
				}

			}

			// �� � ���� ������ �� �����������. 3 ������ 2016.
			//isimmetricGS_switch = 1;

		}
		else {

			if (idirect == 1) {
				// ���������� �����.

				// ������� F ����� C.

//---->#pragma omp parallel for
				for (integer ii = endpos; ii >= startpos; ii--) {
					if ((!my_amg_manager.bcf_reorder) || (F_false_C_true[ii] == false)) { // F nodes

						integer istr = ii - iadd;
						doublerealT rold = x[istr];

						// 13.07.2016
						doublerealT ap_now = 1.0 / Amat.aij[row_ptr_start[ii]];

						x[istr] = b[istr];

						doublerealT rsum = 0.0;
						integer is1 = row_ptr_start[ii] + 1;
						integer is2 = row_ptr_end[ii];
						// ����������������� �������� �������� ����� ������.
						//#pragma omp parallel for reduction(+:rsum)
						for (integer ii1 = is1; ii1 <= is2; ++ii1)
						{
							//x[istr] += -Amat.aij[ii1]*x[Amat.j[ii1]];
							integer ipos = Amat.j[ii1];
							// 13.07.2016
							// ������������� positive connections.
							//if ((Amat.aij[ii1] < 0.0)) {
							rsum += -Amat.aij[ii1] * x[ipos];
							//}
							//else {
							// �� �������.
							//	ap_now += Amat.aij[ii1];
							//}
						}
						x[istr] += rsum;
						//x[istr] *= Amat.aij[row_ptr_start[ii]];
						// 13.07.2016
						x[istr] /= ap_now;

						// ����������� ������ ����� ������� ���� ����� ��������� ������
						// �.�. ��������� �������� �� �������� ����������.
						//if (is1 <= is2) {
							// ������ ���� ��� �� ������� ������� ��������� ����������.
						x[istr] = omega * x[istr] + (1.0 - omega) * rold; // this is SOR
					//}
					}
				}


				//----->#pragma omp parallel for
				for (integer ii = endpos; ii >= startpos; ii--) {
					if ((!my_amg_manager.bcf_reorder) || (F_false_C_true[ii])) { // C nodes

						integer istr = ii - iadd;
						doublerealT rold = x[istr];

						// 13.07.2016
						doublerealT ap_now = 1.0 / Amat.aij[row_ptr_start[ii]];

						x[istr] = b[istr];

						doublerealT rsum = 0.0;
						integer is1 = row_ptr_start[ii] + 1;
						integer is2 = row_ptr_end[ii];
						// ����������������� �������� �������� ����� ������.
						//#pragma omp parallel for reduction(+:rsum)
						for (integer ii1 = is1; ii1 <= is2; ++ii1)
						{
							//x[istr] += -Amat.aij[ii1]*x[Amat.j[ii1]];
							integer ipos = Amat.j[ii1];
							// 13.07.2016
							// ������������� positive connections.
							//if ((Amat.aij[ii1] < 0.0)) {
							rsum += -Amat.aij[ii1] * x[ipos];
							//}
							//else {
							// �� �������.
							//	ap_now += Amat.aij[ii1];
							//}
						}
						x[istr] += rsum;
						//x[istr] *= Amat.aij[row_ptr_start[ii]];
						// 13.07.2016
						x[istr] /= ap_now;

						// ����������� ������ ����� ������� ���� ����� ��������� ������
						// �.�. ��������� �������� �� �������� ����������.
						//if (is1 <= is2) {
							// ������ ���� ��� �� ������� ������� ��������� ����������.
						x[istr] = omega * x[istr] + (1.0 - omega) * rold; // this is SOR
					//}
					}
				}



			}
			else {
				// idirect==0
				// ������� � ����� F.

//--->#pragma omp parallel for
				for (integer ii = endpos; ii >= startpos; ii--) {
					if ((!my_amg_manager.bcf_reorder) || (F_false_C_true[ii])) { // C nodes

						integer istr = ii - iadd;
						doublerealT rold = x[istr];

						// 13.07.2016
						doublerealT ap_now = 1.0 / Amat.aij[row_ptr_start[ii]];

						x[istr] = b[istr];

						doublerealT rsum = 0.0;
						integer is1 = row_ptr_start[ii] + 1;
						integer is2 = row_ptr_end[ii];
						// ����������������� �������� �������� ����� ������.
						//#pragma omp parallel for reduction(+:rsum)
						for (integer ii1 = is1; ii1 <= is2; ++ii1)
						{
							//x[istr] += -Amat.aij[ii1]*x[Amat.j[ii1]];
							integer ipos = Amat.j[ii1];
							// 13.07.2016
							// ������������� positive connections.
							//if ((Amat.aij[ii1] < 0.0)) {
							rsum += -Amat.aij[ii1] * x[ipos];
							//}
							//else {
							// �� �������.
							//	ap_now += Amat.aij[ii1];
							//}
						}
						x[istr] += rsum;
						//x[istr] *= Amat.aij[row_ptr_start[ii]];
						// 13.07.2016
						x[istr] /= ap_now;

						// ����������� ������ ����� ������� ���� ����� ��������� ������
						// �.�. ��������� �������� �� �������� ����������.
						//if (is1 <= is2) {
							// ������ ���� ��� �� ������� ������� ��������� ����������.
						x[istr] = omega * x[istr] + (1.0 - omega) * rold; // this is SOR
					//}
					}
				}


				//--->#pragma omp parallel for
				for (integer ii = endpos; ii >= startpos; ii--) {
					if ((!my_amg_manager.bcf_reorder) || (F_false_C_true[ii] == false)) { // F nodes

						integer istr = ii - iadd;
						doublerealT rold = x[istr];

						// 13.07.2016
						doublerealT ap_now = 1.0 / Amat.aij[row_ptr_start[ii]];

						x[istr] = b[istr];

						doublerealT rsum = 0.0;
						integer is1 = row_ptr_start[ii] + 1;
						integer is2 = row_ptr_end[ii];
						// ����������������� �������� �������� ����� ������.
						//#pragma omp parallel for reduction(+:rsum)
						for (integer ii1 = is1; ii1 <= is2; ++ii1)
						{
							//x[istr] += -Amat.aij[ii1]*x[Amat.j[ii1]];
							integer ipos = Amat.j[ii1];
							// 13.07.2016
							// ������������� positive connections.
							//if ((Amat.aij[ii1] < 0.0)) {
							rsum += -Amat.aij[ii1] * x[ipos];
							//}
							//else {
							// �� �������.
							//	ap_now += Amat.aij[ii1];
							//}
						}
						x[istr] += rsum;
						//x[istr] *= Amat.aij[row_ptr_start[ii]];
						// 13.07.2016
						x[istr] /= ap_now;

						// ����������� ������ ����� ������� ���� ����� ��������� ������
						// �.�. ��������� �������� �� �������� ����������.
						//if (is1 <= is2) {
							// ������ ���� ��� �� ������� ������� ��������� ����������.
						x[istr] = omega * x[istr] + (1.0 - omega) * rold; // this is SOR
					//}
					}
				}

			}

			//isimmetricGS_switch = 0;  // ����� �����������.
		}
	}


} // seidelqsor2

// smoother.
// 5 ���� 2017 ��������� CF-Jacobi smoothing (F - smoothing).
// 5 ������ 2016 � �������������� ������� �� ����� ������� �����.
// 9 september 2015.
// q - quick.
template <typename doublerealT>
void seidelqsor2(Ak2& Amat, integer istartq, integer iendq, doublerealT*& x, doublerealT*& b,
	integer*& row_ptr_start, integer*& row_ptr_end, integer iadd, bool*& F_false_C_true,
	integer idirect, doublerealT*& diag_minus_one, integer ibsp_length)
{
	// F_false_C_true - ��������� ���������� � 1.
	// idirect==0 douwn
	// idirect==1 up

	// istart - ��������� ������� ��������� ��������� � ������� �.
	// iend - �������� ������� ��������� ��������� � ������� �.
	// sor 1.855 ������������ ������ �� ����� ������������ ��������������� ��������.
	// ��������� ������ ����������.
	// ������������ ����� � ��� ������ ����������. 0.8

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

							// �� ������������� ������ ����� ������� ����� ���. 183.
	doublerealT rn = (doublerealT)(iendq - istartq + 1);
	optimal_omega(rn, omega); //28.07.2016
							  //omega = 0.7;

							  //if (isorintmemo == iadd) {
							  // ��� ����� �� ������ ���
							  //bfirst = false;
							  //}
	integer startpos = istartq + iadd;
	integer endpos = iendq + iadd;

	if (omega < 1.0) {
		if (bfirst_jacoby_start) {
			x_jacoby_buffer = new doublereal[3 * (endpos - startpos + 1)];
			i_x_jacoby_buffer_pool_size = 3 * (endpos - startpos + 1);
			bfirst_jacoby_start = false;
		}
		else {
			// ������������� ����������� ������ � ������ nu1==0.
			if (i_x_jacoby_buffer_pool_size < 3 * (endpos - startpos + 1)) {
				if (x_jacoby_buffer != nullptr) {
					delete[] x_jacoby_buffer;
					x_jacoby_buffer = nullptr;
					x_jacoby_buffer = new doublereal[3 * (endpos - startpos + 1)];
					i_x_jacoby_buffer_pool_size = 3 * (endpos - startpos + 1);
					bfirst_jacoby_start = false;
				}
			}
		}
		// copy

		if (ibsp_length == 0)
		{
			if (idirect == 1) {
				// ���������� �����.

				// ������� F ����� C.

				//#pragma loop(hint_parallel(8))
#pragma omp parallel for
				for (integer ii = startpos; ii <= endpos; ++ii) {
					integer istr = ii - iadd;
					x_jacoby_buffer[istr] = x[istr];
				}

				//#pragma loop(hint_parallel(8))
#pragma omp parallel for
				for (integer ii = startpos; ii <= endpos; ++ii) {
					if ((!my_amg_manager.bcf_reorder) || (F_false_C_true[ii])) { // F nodes

						integer istr = ii - iadd;
						doublerealT rold = x_jacoby_buffer[istr];

						// 13.07.2016
						doublerealT ap_now = 1.0 / Amat.aij[row_ptr_start[ii]];

						// ibsp_length==0
						x[istr] = b[istr];

						doublerealT rsum = 0.0;
						integer is1 = row_ptr_start[ii] + 1;
						integer is2 = row_ptr_end[ii];
						// ����������������� �������� �������� ����� ������.
						//#pragma omp parallel for reduction(+:rsum)
						for (integer ii1 = is1; ii1 <= is2; ++ii1)
						{
							//x[istr] += -Amat.aij[ii1]*x_jacoby_buffer[Amat.j[ii1]];
							integer ipos = Amat.j[ii1];
							// 13.07.2016
							// ������������� positive connections.
							//if ((Amat.aij[ii1] < 0.0)) {
							rsum += -Amat.aij[ii1] * x_jacoby_buffer[ipos];
							//}
							//else {
							// �� �������.
							//	ap_now += Amat.aij[ii1];
							//}
						}
						x[istr] += rsum;
						//x[istr] *= Amat.aij[row_ptr_start[ii]];
						// 13.07.2016
						x[istr] /= ap_now;

						// ����������� ������ ����� ������� ���� ����� ��������� ������
						// �.�. ��������� �������� �� �������� ����������.
						//if (is1 <= is2) {
							// ������ ���� ��� �� ������� ������� ��������� ����������.
						x[istr] = omega * x[istr] + (1.0 - omega) * rold; // this is SOR
					//}
					}
				}

				// update.
				//#pragma loop(hint_parallel(8))
#pragma omp parallel for
				for (integer ii = startpos; ii <= endpos; ++ii) {
					integer istr = ii - iadd;
					x_jacoby_buffer[istr] = x[istr];
				}

				//#pragma loop(hint_parallel(8))
#pragma omp parallel for
				for (integer ii = startpos; ii <= endpos; ++ii) {
					if ((!my_amg_manager.bcf_reorder) || (F_false_C_true[ii])) { // C nodes

						integer istr = ii - iadd;
						doublerealT rold = x_jacoby_buffer[istr];

						// 13.07.2016
						doublerealT ap_now = 1.0 / Amat.aij[row_ptr_start[ii]];

						// ibsp_length==0
						x[istr] = b[istr];

						doublerealT rsum = 0.0;
						integer is1 = row_ptr_start[ii] + 1;
						integer is2 = row_ptr_end[ii];
						// ����������������� �������� �������� ����� ������.
						//#pragma omp parallel for reduction(+:rsum)
						for (integer ii1 = is1; ii1 <= is2; ++ii1)
						{
							//x[istr] += -Amat.aij[ii1]*x_jacoby_buffer[Amat.j[ii1]];
							integer ipos = Amat.j[ii1];
							// 13.07.2016
							// ������������� positive connections.
							//if ((Amat.aij[ii1] < 0.0)) {
							rsum += -Amat.aij[ii1] * x_jacoby_buffer[ipos];
							//}
							//else {
							// �� �������.
							//	ap_now += Amat.aij[ii1];
							//}
						}
						x[istr] += rsum;
						//x[istr] *= Amat.aij[row_ptr_start[ii]];
						// 13.07.2016
						x[istr] /= ap_now;

						// ����������� ������ ����� ������� ���� ����� ��������� ������
						// �.�. ��������� �������� �� �������� ����������.
						//if (is1 <= is2) {
							// ������ ���� ��� �� ������� ������� ��������� ����������.
						x[istr] = omega * x[istr] + (1.0 - omega) * rold; // this is SOR
					//}
					}
				}



			}
			else {
				// idirect==0
				// ������� � ����� F.

				//#pragma loop(hint_parallel(8))
#pragma omp parallel for
				for (integer ii = startpos; ii <= endpos; ++ii) {
					integer istr = ii - iadd;
					x_jacoby_buffer[istr] = x[istr];
				}

				//#pragma loop(hint_parallel(8))
#pragma omp parallel for
				for (integer ii = startpos; ii <= endpos; ++ii) {
					if ((!my_amg_manager.bcf_reorder) || (F_false_C_true[ii])) { // C nodes

						integer istr = ii - iadd;
						doublerealT rold = x_jacoby_buffer[istr];

						// 13.07.2016
						doublerealT ap_now = 1.0 / Amat.aij[row_ptr_start[ii]];

						// ibsp_length==0
						x[istr] = b[istr];

						doublerealT rsum = 0.0;
						integer is1 = row_ptr_start[ii] + 1;
						integer is2 = row_ptr_end[ii];
						// ����������������� �������� �������� ����� ������.
						//#pragma omp parallel for reduction(+:rsum)
						for (integer ii1 = is1; ii1 <= is2; ++ii1)
						{
							//x[istr] += -Amat.aij[ii1]*x_jacoby_buffer[Amat.j[ii1]];
							integer ipos = Amat.j[ii1];
							// 13.07.2016
							// ������������� positive connections.
							//if ((Amat.aij[ii1] < 0.0)) {
							rsum += -Amat.aij[ii1] * x_jacoby_buffer[ipos];
							//}
							//else {
							// �� �������.
							//	ap_now += Amat.aij[ii1];
							//}
						}
						x[istr] += rsum;
						//x[istr] *= Amat.aij[row_ptr_start[ii]];
						// 13.07.2016
						x[istr] /= ap_now;

						// ����������� ������ ����� ������� ���� ����� ��������� ������
						// �.�. ��������� �������� �� �������� ����������.
						//if (is1 <= is2) {
							// ������ ���� ��� �� ������� ������� ��������� ����������.
						x[istr] = omega * x[istr] + (1.0 - omega) * rold; // this is SOR
					//}
					}
				}


				//#pragma loop(hint_parallel(8))
#pragma omp parallel for
				for (integer ii = startpos; ii <= endpos; ++ii) {
					integer istr = ii - iadd;
					x_jacoby_buffer[istr] = x[istr];
				}

				//#pragma loop(hint_parallel(8))
#pragma omp parallel for
				for (integer ii = startpos; ii <= endpos; ++ii) {
					if ((!my_amg_manager.bcf_reorder) || (F_false_C_true[ii] == false)) { // F nodes

						integer istr = ii - iadd;
						doublerealT rold = x_jacoby_buffer[istr];

						// 13.07.2016
						doublerealT ap_now = 1.0 / Amat.aij[row_ptr_start[ii]];

						// ibsp_length==0
						x[istr] = b[istr];

						doublerealT rsum = 0.0;
						integer is1 = row_ptr_start[ii] + 1;
						integer is2 = row_ptr_end[ii];
						// ����������������� �������� �������� ����� ������.
						//#pragma omp parallel for reduction(+:rsum)
						for (integer ii1 = is1; ii1 <= is2; ++ii1)
						{
							//x[istr] += -Amat.aij[ii1]*x_jacoby_buffer[Amat.j[ii1]];
							integer ipos = Amat.j[ii1];
							// 13.07.2016
							// ������������� positive connections.
							//if ((Amat.aij[ii1] < 0.0)) {
							rsum += -Amat.aij[ii1] * x_jacoby_buffer[ipos];
							//}
							//else {
							// �� �������.
							//	ap_now += Amat.aij[ii1];
							//}
						}
						x[istr] += rsum;
						//x[istr] *= Amat.aij[row_ptr_start[ii]];
						// 13.07.2016
						x[istr] /= ap_now;

						// ����������� ������ ����� ������� ���� ����� ��������� ������
						// �.�. ��������� �������� �� �������� ����������.
						//if (is1 <= is2) {
							// ������ ���� ��� �� ������� ������� ��������� ����������.
						x[istr] = omega * x[istr] + (1.0 - omega) * rold; // this is SOR
					//}
					}
				}

			}
		}
		else {

			if (idirect == 1) {
				// ���������� �����.

				// ������� F ����� C.

				//#pragma loop(hint_parallel(8))
#pragma omp parallel for
				for (integer ii = startpos; ii <= endpos; ++ii) {
					integer istr = ii - iadd;
					x_jacoby_buffer[istr] = x[istr];
				}

				//#pragma loop(hint_parallel(8))
#pragma omp parallel for
				for (integer ii = startpos; ii <= endpos; ++ii) {
					if ((!my_amg_manager.bcf_reorder) || (F_false_C_true[ii])) { // F nodes

						integer istr = ii - iadd;
						doublerealT rold = x_jacoby_buffer[istr];

						// 13.07.2016
						doublerealT ap_now = 1.0 / Amat.aij[row_ptr_start[ii]];

						x[istr] = diag_minus_one[istr] * b[istr];

						doublerealT rsum = 0.0;
						integer is1 = row_ptr_start[ii] + 1;
						integer is2 = row_ptr_end[ii];
						// ����������������� �������� �������� ����� ������.
						//#pragma omp parallel for reduction(+:rsum)
						for (integer ii1 = is1; ii1 <= is2; ++ii1)
						{
							//x[istr] += -Amat.aij[ii1]*x_jacoby_buffer[Amat.j[ii1]];
							integer ipos = Amat.j[ii1];
							// 13.07.2016
							// ������������� positive connections.
							//if ((Amat.aij[ii1] < 0.0)) {
							rsum += -Amat.aij[ii1] * x_jacoby_buffer[ipos];
							//}
							//else {
							// �� �������.
							//	ap_now += Amat.aij[ii1];
							//}
						}
						x[istr] += rsum;
						//x[istr] *= Amat.aij[row_ptr_start[ii]];
						// 13.07.2016
						x[istr] /= ap_now;

						// ����������� ������ ����� ������� ���� ����� ��������� ������
						// �.�. ��������� �������� �� �������� ����������.
						//if (is1 <= is2) {
							// ������ ���� ��� �� ������� ������� ��������� ����������.
						x[istr] = omega * x[istr] + (1.0 - omega) * rold; // this is SOR
					//}
					}
				}

				// update.
				//#pragma loop(hint_parallel(8))
#pragma omp parallel for
				for (integer ii = startpos; ii <= endpos; ++ii) {
					integer istr = ii - iadd;
					x_jacoby_buffer[istr] = x[istr];
				}

				//#pragma loop(hint_parallel(8))
#pragma omp parallel for
				for (integer ii = startpos; ii <= endpos; ++ii) {
					if ((!my_amg_manager.bcf_reorder) || (F_false_C_true[ii])) { // C nodes

						integer istr = ii - iadd;
						doublerealT rold = x_jacoby_buffer[istr];

						// 13.07.2016
						doublerealT ap_now = 1.0 / Amat.aij[row_ptr_start[ii]];

						x[istr] = diag_minus_one[istr] * b[istr];

						doublerealT rsum = 0.0;
						integer is1 = row_ptr_start[ii] + 1;
						integer is2 = row_ptr_end[ii];
						// ����������������� �������� �������� ����� ������.
						//#pragma omp parallel for reduction(+:rsum)
						for (integer ii1 = is1; ii1 <= is2; ++ii1)
						{
							//x[istr] += -Amat.aij[ii1]*x_jacoby_buffer[Amat.j[ii1]];
							integer ipos = Amat.j[ii1];
							// 13.07.2016
							// ������������� positive connections.
							//if ((Amat.aij[ii1] < 0.0)) {
							rsum += -Amat.aij[ii1] * x_jacoby_buffer[ipos];
							//}
							//else {
							// �� �������.
							//	ap_now += Amat.aij[ii1];
							//}
						}
						x[istr] += rsum;
						//x[istr] *= Amat.aij[row_ptr_start[ii]];
						// 13.07.2016
						x[istr] /= ap_now;

						// ����������� ������ ����� ������� ���� ����� ��������� ������
						// �.�. ��������� �������� �� �������� ����������.
						//if (is1 <= is2) {
							// ������ ���� ��� �� ������� ������� ��������� ����������.
						x[istr] = omega * x[istr] + (1.0 - omega) * rold; // this is SOR
					//}
					}
				}



			}
			else {
				// idirect==0
				// ������� � ����� F.

				//#pragma loop(hint_parallel(8))
#pragma omp parallel for
				for (integer ii = startpos; ii <= endpos; ++ii) {
					integer istr = ii - iadd;
					x_jacoby_buffer[istr] = x[istr];
				}

				//#pragma loop(hint_parallel(8))
#pragma omp parallel for
				for (integer ii = startpos; ii <= endpos; ++ii) {
					if ((!my_amg_manager.bcf_reorder) || (F_false_C_true[ii])) { // C nodes

						integer istr = ii - iadd;
						doublerealT rold = x_jacoby_buffer[istr];

						// 13.07.2016
						doublerealT ap_now = 1.0 / Amat.aij[row_ptr_start[ii]];

						x[istr] = diag_minus_one[istr] * b[istr];

						doublerealT rsum = 0.0;
						integer is1 = row_ptr_start[ii] + 1;
						integer is2 = row_ptr_end[ii];
						// ����������������� �������� �������� ����� ������.
						//#pragma omp parallel for reduction(+:rsum)
						for (integer ii1 = is1; ii1 <= is2; ++ii1)
						{
							//x[istr] += -Amat.aij[ii1]*x_jacoby_buffer[Amat.j[ii1]];
							integer ipos = Amat.j[ii1];
							// 13.07.2016
							// ������������� positive connections.
							//if ((Amat.aij[ii1] < 0.0)) {
							rsum += -Amat.aij[ii1] * x_jacoby_buffer[ipos];
							//}
							//else {
							// �� �������.
							//	ap_now += Amat.aij[ii1];
							//}
						}
						x[istr] += rsum;
						//x[istr] *= Amat.aij[row_ptr_start[ii]];
						// 13.07.2016
						x[istr] /= ap_now;

						// ����������� ������ ����� ������� ���� ����� ��������� ������
						// �.�. ��������� �������� �� �������� ����������.
						//if (is1 <= is2) {
							// ������ ���� ��� �� ������� ������� ��������� ����������.
						x[istr] = omega * x[istr] + (1.0 - omega) * rold; // this is SOR
					//}
					}
				}


				//#pragma loop(hint_parallel(8))
#pragma omp parallel for
				for (integer ii = startpos; ii <= endpos; ++ii) {
					integer istr = ii - iadd;
					x_jacoby_buffer[istr] = x[istr];
				}

				//#pragma loop(hint_parallel(8))
#pragma omp parallel for
				for (integer ii = startpos; ii <= endpos; ++ii) {
					if ((!my_amg_manager.bcf_reorder) || (F_false_C_true[ii] == false)) { // F nodes

						integer istr = ii - iadd;
						doublerealT rold = x_jacoby_buffer[istr];

						// 13.07.2016
						doublerealT ap_now = 1.0 / Amat.aij[row_ptr_start[ii]];

						x[istr] = diag_minus_one[istr] * b[istr];

						doublerealT rsum = 0.0;
						integer is1 = row_ptr_start[ii] + 1;
						integer is2 = row_ptr_end[ii];
						// ����������������� �������� �������� ����� ������.
						//#pragma omp parallel for reduction(+:rsum)
						for (integer ii1 = is1; ii1 <= is2; ++ii1)
						{
							//x[istr] += -Amat.aij[ii1]*x_jacoby_buffer[Amat.j[ii1]];
							integer ipos = Amat.j[ii1];
							// 13.07.2016
							// ������������� positive connections.
							//if ((Amat.aij[ii1] < 0.0)) {
							rsum += -Amat.aij[ii1] * x_jacoby_buffer[ipos];
							//}
							//else {
							// �� �������.
							//	ap_now += Amat.aij[ii1];
							//}
						}
						x[istr] += rsum;
						//x[istr] *= Amat.aij[row_ptr_start[ii]];
						// 13.07.2016
						x[istr] /= ap_now;

						// ����������� ������ ����� ������� ���� ����� ��������� ������
						// �.�. ��������� �������� �� �������� ����������.
						//if (is1 <= is2) {
							// ������ ���� ��� �� ������� ������� ��������� ����������.
						x[istr] = omega * x[istr] + (1.0 - omega) * rold; // this is SOR
					//}
					}
				}

			}

		}
	}
	else {

		if (ibsp_length == 0) {
			// 3 ������ 2016. ������������ ����� ������-�������.
			if (isimmetricGS_switch == 0) {
				// 3 ������ 2016 ���� ���������������� �������� �� BSKDmitrii ��� ������������ ����� ������ - ������� ����������.
				// ������ ������������ ������������ ����� ������ -�������. ����������� ������� ����� �������.

				if (idirect == 1) {
					// ���������� �����.



					// ������� F ����� C.

	//----->#pragma omp parallel for
					for (integer ii = startpos; ii <= endpos; ++ii) {
						if ((!my_amg_manager.bcf_reorder) || (F_false_C_true[ii] == false)) { // F nodes

							integer istr = ii - iadd;
							doublerealT rold = x[istr];

							// 13.07.2016
							doublerealT ap_now = 1.0 / Amat.aij[row_ptr_start[ii]];
							// 28.01.2017
							//ap_now = 0.0;

							// ibsp_length==0
							x[istr] = b[istr];

							doublerealT rsum = 0.0;
							integer is1 = row_ptr_start[ii] + 1;
							integer is2 = row_ptr_end[ii];
							// ����������������� �������� �������� ����� ������.
							//#pragma omp parallel for reduction(+:rsum)
							for (integer ii1 = is1; ii1 <= is2; ++ii1)
							{
								//x[istr] += -Amat.aij[ii1]*x[Amat.j[ii1]];
								integer ipos = Amat.j[ii1];
								// 13.07.2016
								if (1) {
									// ������������� positive connections.
									//if ((Amat.aij[ii1] < 0.0)) {
									// ����� �� � ���� � ������� � ������� ������������.
									rsum += -Amat.aij[ii1] * x[ipos];
									//}
									//else {
									// �� ��������.
									//	ap_now += Amat.aij[ii1];
									//}
								}
								else {
									// ����� � ����������. �� �������� ����� �������: TKM, TKM1, TKM2.

									// ������������� positive connections.
									if ((Amat.aij[ii1] < 0.0)) {
										rsum += -Amat.aij[ii1] * x[ipos];
										//ap_now += fabs(Amat.aij[ii1]);
									}
									else {
										// �� ��������.
										// �������� ��-�� ���� ��� ��� ������� ������.
										//ap_now += fabs(Amat.aij[ii1]);
										if (fabs(x[ipos]) > fabs(x[istr])) {
											ap_now += fabs(Amat.aij[ii1]);
										}
										else rsum += -Amat.aij[ii1] * x[ipos];
									}
								}
							}
							x[istr] += rsum;
							//x[istr] *= Amat.aij[row_ptr_start[ii]];
							// 13.07.2016
							x[istr] /= ap_now;

							// ����������� ������ ����� ������� ���� ����� ��������� ������
							// �.�. ��������� �������� �� �������� ����������.
							//if (is1 <= is2) {
								// ������ ���� ��� �� ������� ������� ��������� ����������.
							x[istr] = omega * x[istr] + (1.0 - omega) * rold; // this is SOR
						//}
						}
					}

					//---->#pragma omp parallel for
					for (integer ii = startpos; ii <= endpos; ++ii) {
						if ((!my_amg_manager.bcf_reorder) || (F_false_C_true[ii])) { // C nodes

							integer istr = ii - iadd;
							doublerealT rold = x[istr];

							// 13.07.2016
							doublerealT ap_now = 1.0 / Amat.aij[row_ptr_start[ii]];
							// 28.01.2017
							//ap_now = 0.0;

							// ibsp_length==0
							x[istr] = b[istr];

							doublerealT rsum = 0.0;
							integer is1 = row_ptr_start[ii] + 1;
							integer is2 = row_ptr_end[ii];
							// ����������������� �������� �������� ����� ������.
							//#pragma omp parallel for reduction(+:rsum)
							for (integer ii1 = is1; ii1 <= is2; ++ii1)
							{
								//x[istr] += -Amat.aij[ii1]*x[Amat.j[ii1]];
								integer ipos = Amat.j[ii1];
								// 13.07.2016
								if (1) {
									// ������������� positive connections.
									//if ((Amat.aij[ii1] < 0.0)) {
									// ����� �� � ���� � ������� � ������� ������������.
									rsum += -Amat.aij[ii1] * x[ipos];
									//}
									//else {
									// �� ��������.
									//	ap_now += Amat.aij[ii1];
									//}
								}
								else {
									// ����� � ����������. �� �������� ����� �������: TKM, TKM1, TKM2.

									// ������������� positive connections.
									if ((Amat.aij[ii1] < 0.0)) {
										rsum += -Amat.aij[ii1] * x[ipos];
										//ap_now += fabs(Amat.aij[ii1]);
									}
									else {
										// �� ��������.
										// �������� ��-�� ���� ��� ��� ������� ������.
										//ap_now += fabs(Amat.aij[ii1]);
										if (fabs(x[ipos]) > fabs(x[istr])) {
											ap_now += fabs(Amat.aij[ii1]);
										}
										else rsum += -Amat.aij[ii1] * x[ipos];
									}
								}
							}
							x[istr] += rsum;
							//x[istr] *= Amat.aij[row_ptr_start[ii]];
							// 13.07.2016
							x[istr] /= ap_now;

							// ����������� ������ ����� ������� ���� ����� ��������� ������
							// �.�. ��������� �������� �� �������� ����������.
							//if (is1 <= is2) {
								// ������ ���� ��� �� ������� ������� ��������� ����������.
							x[istr] = omega * x[istr] + (1.0 - omega) * rold; // this is SOR
						//}
						}
					}

				}
				else {

					// idirect==0
					// ������� � ����� F.

	//---->#pragma omp parallel for
					for (integer ii = startpos; ii <= endpos; ++ii) {
						if ((!my_amg_manager.bcf_reorder) || (F_false_C_true[ii])) { // C nodes

							integer istr = ii - iadd;
							doublerealT rold = x[istr];

							// 13.07.2016
							doublerealT ap_now = 1.0 / Amat.aij[row_ptr_start[ii]];
							// 28.01.2017
							//ap_now = 0.0;

							// ibsp_length==0
							x[istr] = b[istr];

							doublerealT rsum = 0.0;
							integer is1 = row_ptr_start[ii] + 1;
							integer is2 = row_ptr_end[ii];
							// ����������������� �������� �������� ����� ������.
							//#pragma omp parallel for reduction(+:rsum)
							for (integer ii1 = is1; ii1 <= is2; ++ii1)
							{
								//x[istr] += -Amat.aij[ii1]*x[Amat.j[ii1]];
								integer ipos = Amat.j[ii1];
								// 13.07.2016
								if (1) {
									// ������������� positive connections.
									//if ((Amat.aij[ii1] < 0.0)) {
									// ����� �� � ���� � ������� � ������� ������������.
									rsum += -Amat.aij[ii1] * x[ipos];
									//}
									//else {
									// �� ��������.
									//	ap_now += Amat.aij[ii1];
									//}
								}
								else {
									// ����� � ����������. �� �������� ����� �������: TKM, TKM1, TKM2.

									// ������������� positive connections.
									if ((Amat.aij[ii1] < 0.0)) {
										rsum += -Amat.aij[ii1] * x[ipos];
										//ap_now += fabs(Amat.aij[ii1]);
									}
									else {
										// �� ��������.
										// �������� ��-�� ���� ��� ��� ������� ������.
										//ap_now += fabs(Amat.aij[ii1]);
										if (fabs(x[ipos]) > fabs(x[istr])) {
											ap_now += fabs(Amat.aij[ii1]);
										}
										else rsum += -Amat.aij[ii1] * x[ipos];
									}
								}
							}
							x[istr] += rsum;
							//x[istr] *= Amat.aij[row_ptr_start[ii]];
							// 13.07.2016
							x[istr] /= ap_now;

							// ����������� ������ ����� ������� ���� ����� ��������� ������
							// �.�. ��������� �������� �� �������� ����������.
							//if (is1 <= is2) {
								// ������ ���� ��� �� ������� ������� ��������� ����������.
							x[istr] = omega * x[istr] + (1.0 - omega) * rold; // this is SOR
						//}
						}
					}

					//---->#pragma omp parallel for
					for (integer ii = startpos; ii <= endpos; ++ii) {
						if ((!my_amg_manager.bcf_reorder) || (F_false_C_true[ii] == false)) { // F nodes

							integer istr = ii - iadd;
							doublerealT rold = x[istr];

							// 13.07.2016
							doublerealT ap_now = 1.0 / Amat.aij[row_ptr_start[ii]];
							// 28.01.2017
							//ap_now = 0.0;

							// ibsp_length==0
							x[istr] = b[istr];

							doublerealT rsum = 0.0;
							integer is1 = row_ptr_start[ii] + 1;
							integer is2 = row_ptr_end[ii];
							// ����������������� �������� �������� ����� ������.
							//#pragma omp parallel for reduction(+:rsum)
							for (integer ii1 = is1; ii1 <= is2; ++ii1)
							{
								//x[istr] += -Amat.aij[ii1]*x[Amat.j[ii1]];
								integer ipos = Amat.j[ii1];
								// 13.07.2016
								if (1) {
									// ������������� positive connections.
									//if ((Amat.aij[ii1] < 0.0)) {
									// ����� �� � ���� � ������� � ������� ������������.
									rsum += -Amat.aij[ii1] * x[ipos];
									//}
									//else {
									// �� ��������.
									//	ap_now += Amat.aij[ii1];
									//}
								}
								else {
									// ����� � ����������. �� �������� ����� �������: TKM, TKM1, TKM2.

									// ������������� positive connections.
									if ((Amat.aij[ii1] < 0.0)) {
										rsum += -Amat.aij[ii1] * x[ipos];
										//ap_now += fabs(Amat.aij[ii1]);
									}
									else {
										// �� ��������.
										// �������� ��-�� ���� ��� ��� ������� ������.
										//ap_now += fabs(Amat.aij[ii1]);
										if (fabs(x[ipos]) > fabs(x[istr])) {
											ap_now += fabs(Amat.aij[ii1]);
										}
										else rsum += -Amat.aij[ii1] * x[ipos];
									}
								}
							}
							x[istr] += rsum;
							//x[istr] *= Amat.aij[row_ptr_start[ii]];
							// 13.07.2016
							x[istr] /= ap_now;

							// ����������� ������ ����� ������� ���� ����� ��������� ������
							// �.�. ��������� �������� �� �������� ����������.
							//if (is1 <= is2) {
								// ������ ���� ��� �� ������� ������� ��������� ����������.
							x[istr] = omega * x[istr] + (1.0 - omega) * rold; // this is SOR
						//}
						}
					}

				}

				// �� � ���� ������ �� �����������. 3 ������ 2016.
				//isimmetricGS_switch = 1;

			}
			else {

				if (idirect == 1) {
					// ���������� �����.

					// ������� F ����� C.

	//---->#pragma omp parallel for
					for (integer ii = endpos; ii >= startpos; --ii) {
						if ((!my_amg_manager.bcf_reorder) || (F_false_C_true[ii] == false)) { // F nodes

							integer istr = ii - iadd;
							doublerealT rold = x[istr];

							// 13.07.2016
							doublerealT ap_now = 1.0 / Amat.aij[row_ptr_start[ii]];

							// ibsp_length==0
							x[istr] = b[istr];

							doublerealT rsum = 0.0;
							integer is1 = row_ptr_start[ii] + 1;
							integer is2 = row_ptr_end[ii];
							// ����������������� �������� �������� ����� ������.
							//#pragma omp parallel for reduction(+:rsum)
							for (integer ii1 = is1; ii1 <= is2; ++ii1)
							{
								//x[istr] += -Amat.aij[ii1]*x[Amat.j[ii1]];
								integer ipos = Amat.j[ii1];
								// 13.07.2016
								// ������������� positive connections.
								//if ((Amat.aij[ii1] < 0.0)) {
								rsum += -Amat.aij[ii1] * x[ipos];
								//}
								//else {
								// �� �������.
								//	ap_now += Amat.aij[ii1];
								//}
							}
							x[istr] += rsum;
							//x[istr] *= Amat.aij[row_ptr_start[ii]];
							// 13.07.2016
							x[istr] /= ap_now;

							// ����������� ������ ����� ������� ���� ����� ��������� ������
							// �.�. ��������� �������� �� �������� ����������.
							//if (is1 <= is2) {
								// ������ ���� ��� �� ������� ������� ��������� ����������.
							x[istr] = omega * x[istr] + (1.0 - omega) * rold; // this is SOR
						//}
						}
					}


					//----->#pragma omp parallel for
					for (integer ii = endpos; ii >= startpos; --ii) {
						if ((!my_amg_manager.bcf_reorder) || (F_false_C_true[ii])) { // C nodes

							integer istr = ii - iadd;
							doublerealT rold = x[istr];

							// 13.07.2016
							doublerealT ap_now = 1.0 / Amat.aij[row_ptr_start[ii]];

							// ibsp_length==0
							x[istr] = b[istr];

							doublerealT rsum = 0.0;
							integer is1 = row_ptr_start[ii] + 1;
							integer is2 = row_ptr_end[ii];
							// ����������������� �������� �������� ����� ������.
							//#pragma omp parallel for reduction(+:rsum)
							for (integer ii1 = is1; ii1 <= is2; ++ii1)
							{
								//x[istr] += -Amat.aij[ii1]*x[Amat.j[ii1]];
								integer ipos = Amat.j[ii1];
								// 13.07.2016
								// ������������� positive connections.
								//if ((Amat.aij[ii1] < 0.0)) {
								rsum += -Amat.aij[ii1] * x[ipos];
								//}
								//else {
								// �� �������.
								//	ap_now += Amat.aij[ii1];
								//}
							}
							x[istr] += rsum;
							//x[istr] *= Amat.aij[row_ptr_start[ii]];
							// 13.07.2016
							x[istr] /= ap_now;

							// ����������� ������ ����� ������� ���� ����� ��������� ������
							// �.�. ��������� �������� �� �������� ����������.
							//if (is1 <= is2) {
								// ������ ���� ��� �� ������� ������� ��������� ����������.
							x[istr] = omega * x[istr] + (1.0 - omega) * rold; // this is SOR
						//}
						}
					}



				}
				else {
					// idirect==0
					// ������� � ����� F.

	//--->#pragma omp parallel for
					for (integer ii = endpos; ii >= startpos; --ii) {
						if ((!my_amg_manager.bcf_reorder) || (F_false_C_true[ii])) { // C nodes

							integer istr = ii - iadd;
							doublerealT rold = x[istr];

							// 13.07.2016
							doublerealT ap_now = 1.0 / Amat.aij[row_ptr_start[ii]];

							// ibsp_length==0
							x[istr] = b[istr];

							doublerealT rsum = 0.0;
							integer is1 = row_ptr_start[ii] + 1;
							integer is2 = row_ptr_end[ii];
							// ����������������� �������� �������� ����� ������.
							//#pragma omp parallel for reduction(+:rsum)
							for (integer ii1 = is1; ii1 <= is2; ++ii1)
							{
								//x[istr] += -Amat.aij[ii1]*x[Amat.j[ii1]];
								integer ipos = Amat.j[ii1];
								// 13.07.2016
								// ������������� positive connections.
								//if ((Amat.aij[ii1] < 0.0)) {
								rsum += -Amat.aij[ii1] * x[ipos];
								//}
								//else {
								// �� �������.
								//	ap_now += Amat.aij[ii1];
								//}
							}
							x[istr] += rsum;
							//x[istr] *= Amat.aij[row_ptr_start[ii]];
							// 13.07.2016
							x[istr] /= ap_now;

							// ����������� ������ ����� ������� ���� ����� ��������� ������
							// �.�. ��������� �������� �� �������� ����������.
							//if (is1 <= is2) {
								// ������ ���� ��� �� ������� ������� ��������� ����������.
							x[istr] = omega * x[istr] + (1.0 - omega) * rold; // this is SOR
						//}
						}
					}


					//--->#pragma omp parallel for
					for (integer ii = endpos; ii >= startpos; --ii) {
						if ((!my_amg_manager.bcf_reorder) || (F_false_C_true[ii] == false)) { // F nodes

							integer istr = ii - iadd;
							doublerealT rold = x[istr];

							// 13.07.2016
							doublerealT ap_now = 1.0 / Amat.aij[row_ptr_start[ii]];

							// ibsp_length==0
							x[istr] = b[istr];

							doublerealT rsum = 0.0;
							integer is1 = row_ptr_start[ii] + 1;
							integer is2 = row_ptr_end[ii];
							// ����������������� �������� �������� ����� ������.
							//#pragma omp parallel for reduction(+:rsum)
							for (integer ii1 = is1; ii1 <= is2; ++ii1)
							{
								//x[istr] += -Amat.aij[ii1]*x[Amat.j[ii1]];
								integer ipos = Amat.j[ii1];
								// 13.07.2016
								// ������������� positive connections.
								//if ((Amat.aij[ii1] < 0.0)) {
								rsum += -Amat.aij[ii1] * x[ipos];
								//}
								//else {
								// �� �������.
								//	ap_now += Amat.aij[ii1];
								//}
							}
							x[istr] += rsum;
							//x[istr] *= Amat.aij[row_ptr_start[ii]];
							// 13.07.2016
							x[istr] /= ap_now;

							// ����������� ������ ����� ������� ���� ����� ��������� ������
							// �.�. ��������� �������� �� �������� ����������.
							//if (is1 <= is2) {
								// ������ ���� ��� �� ������� ������� ��������� ����������.
							x[istr] = omega * x[istr] + (1.0 - omega) * rold; // this is SOR
						//}
						}
					}

				}

				//isimmetricGS_switch = 0;  // ����� �����������.
			}
		}
		else {
			// 3 ������ 2016. ������������ ����� ������-�������.
			if (isimmetricGS_switch == 0) {
				// 3 ������ 2016 ���� ���������������� �������� �� BSKDmitrii ��� ������������ ����� ������ - ������� ����������.
				// ������ ������������ ������������ ����� ������ -�������. ����������� ������� ����� �������.

				if (idirect == 1) {
					// ���������� �����.



					// ������� F ����� C.

	//----->#pragma omp parallel for
					for (integer ii = startpos; ii <= endpos; ++ii) {
						if ((!my_amg_manager.bcf_reorder) || (F_false_C_true[ii] == false)) { // F nodes

							integer istr = ii - iadd;
							doublerealT rold = x[istr];

							// 13.07.2016
							doublerealT ap_now = 1.0 / Amat.aij[row_ptr_start[ii]];
							// 28.01.2017
							//ap_now = 0.0;

							x[istr] = diag_minus_one[istr] * b[istr];

							doublerealT rsum = 0.0;
							integer is1 = row_ptr_start[ii] + 1;
							integer is2 = row_ptr_end[ii];
							// ����������������� �������� �������� ����� ������.
							//#pragma omp parallel for reduction(+:rsum)
							for (integer ii1 = is1; ii1 <= is2; ++ii1)
							{
								//x[istr] += -Amat.aij[ii1]*x[Amat.j[ii1]];
								integer ipos = Amat.j[ii1];
								// 13.07.2016
								if (1) {
									// ������������� positive connections.
									//if ((Amat.aij[ii1] < 0.0)) {
									// ����� �� � ���� � ������� � ������� ������������.
									rsum += -Amat.aij[ii1] * x[ipos];
									//}
									//else {
									// �� ��������.
									//	ap_now += Amat.aij[ii1];
									//}
								}
								else {
									// ����� � ����������. �� �������� ����� �������: TKM, TKM1, TKM2.

									// ������������� positive connections.
									if ((Amat.aij[ii1] < 0.0)) {
										rsum += -Amat.aij[ii1] * x[ipos];
										//ap_now += fabs(Amat.aij[ii1]);
									}
									else {
										// �� ��������.
										// �������� ��-�� ���� ��� ��� ������� ������.
										//ap_now += fabs(Amat.aij[ii1]);
										if (fabs(x[ipos]) > fabs(x[istr])) {
											ap_now += fabs(Amat.aij[ii1]);
										}
										else rsum += -Amat.aij[ii1] * x[ipos];
									}
								}
							}
							x[istr] += rsum;
							//x[istr] *= Amat.aij[row_ptr_start[ii]];
							// 13.07.2016
							x[istr] /= ap_now;

							// ����������� ������ ����� ������� ���� ����� ��������� ������
							// �.�. ��������� �������� �� �������� ����������.
							//if (is1 <= is2) {
								// ������ ���� ��� �� ������� ������� ��������� ����������.
							x[istr] = omega * x[istr] + (1.0 - omega) * rold; // this is SOR
						//}
						}
					}

					//---->#pragma omp parallel for
					for (integer ii = startpos; ii <= endpos; ++ii) {
						if ((!my_amg_manager.bcf_reorder) || (F_false_C_true[ii])) { // C nodes

							integer istr = ii - iadd;
							doublerealT rold = x[istr];

							// 13.07.2016
							doublerealT ap_now = 1.0 / Amat.aij[row_ptr_start[ii]];
							// 28.01.2017
							//ap_now = 0.0;

							x[istr] = diag_minus_one[istr] * b[istr];

							doublerealT rsum = 0.0;
							integer is1 = row_ptr_start[ii] + 1;
							integer is2 = row_ptr_end[ii];
							// ����������������� �������� �������� ����� ������.
							//#pragma omp parallel for reduction(+:rsum)
							for (integer ii1 = is1; ii1 <= is2; ++ii1)
							{
								//x[istr] += -Amat.aij[ii1]*x[Amat.j[ii1]];
								integer ipos = Amat.j[ii1];
								// 13.07.2016
								if (1) {
									// ������������� positive connections.
									//if ((Amat.aij[ii1] < 0.0)) {
									// ����� �� � ���� � ������� � ������� ������������.
									rsum += -Amat.aij[ii1] * x[ipos];
									//}
									//else {
									// �� ��������.
									//	ap_now += Amat.aij[ii1];
									//}
								}
								else {
									// ����� � ����������. �� �������� ����� �������: TKM, TKM1, TKM2.

									// ������������� positive connections.
									if ((Amat.aij[ii1] < 0.0)) {
										rsum += -Amat.aij[ii1] * x[ipos];
										//ap_now += fabs(Amat.aij[ii1]);
									}
									else {
										// �� ��������.
										// �������� ��-�� ���� ��� ��� ������� ������.
										//ap_now += fabs(Amat.aij[ii1]);
										if (fabs(x[ipos]) > fabs(x[istr])) {
											ap_now += fabs(Amat.aij[ii1]);
										}
										else rsum += -Amat.aij[ii1] * x[ipos];
									}
								}
							}
							x[istr] += rsum;
							//x[istr] *= Amat.aij[row_ptr_start[ii]];
							// 13.07.2016
							x[istr] /= ap_now;

							// ����������� ������ ����� ������� ���� ����� ��������� ������
							// �.�. ��������� �������� �� �������� ����������.
							//if (is1 <= is2) {
								// ������ ���� ��� �� ������� ������� ��������� ����������.
							x[istr] = omega * x[istr] + (1.0 - omega) * rold; // this is SOR
						//}
						}
					}

				}
				else {

					// idirect==0
					// ������� � ����� F.

	//---->#pragma omp parallel for
					for (integer ii = startpos; ii <= endpos; ++ii) {
						if ((!my_amg_manager.bcf_reorder) || (F_false_C_true[ii])) { // C nodes

							integer istr = ii - iadd;
							doublerealT rold = x[istr];

							// 13.07.2016
							doublerealT ap_now = 1.0 / Amat.aij[row_ptr_start[ii]];
							// 28.01.2017
							//ap_now = 0.0;

							x[istr] = diag_minus_one[istr] * b[istr];

							doublerealT rsum = 0.0;
							integer is1 = row_ptr_start[ii] + 1;
							integer is2 = row_ptr_end[ii];
							// ����������������� �������� �������� ����� ������.
							//#pragma omp parallel for reduction(+:rsum)
							for (integer ii1 = is1; ii1 <= is2; ++ii1)
							{
								//x[istr] += -Amat.aij[ii1]*x[Amat.j[ii1]];
								integer ipos = Amat.j[ii1];
								// 13.07.2016
								if (1) {
									// ������������� positive connections.
									//if ((Amat.aij[ii1] < 0.0)) {
									// ����� �� � ���� � ������� � ������� ������������.
									rsum += -Amat.aij[ii1] * x[ipos];
									//}
									//else {
									// �� ��������.
									//	ap_now += Amat.aij[ii1];
									//}
								}
								else {
									// ����� � ����������. �� �������� ����� �������: TKM, TKM1, TKM2.

									// ������������� positive connections.
									if ((Amat.aij[ii1] < 0.0)) {
										rsum += -Amat.aij[ii1] * x[ipos];
										//ap_now += fabs(Amat.aij[ii1]);
									}
									else {
										// �� ��������.
										// �������� ��-�� ���� ��� ��� ������� ������.
										//ap_now += fabs(Amat.aij[ii1]);
										if (fabs(x[ipos]) > fabs(x[istr])) {
											ap_now += fabs(Amat.aij[ii1]);
										}
										else rsum += -Amat.aij[ii1] * x[ipos];
									}
								}
							}
							x[istr] += rsum;
							//x[istr] *= Amat.aij[row_ptr_start[ii]];
							// 13.07.2016
							x[istr] /= ap_now;

							// ����������� ������ ����� ������� ���� ����� ��������� ������
							// �.�. ��������� �������� �� �������� ����������.
							//if (is1 <= is2) {
								// ������ ���� ��� �� ������� ������� ��������� ����������.
							x[istr] = omega * x[istr] + (1.0 - omega) * rold; // this is SOR
						//}
						}
					}

					//---->#pragma omp parallel for
					for (integer ii = startpos; ii <= endpos; ++ii) {
						if ((!my_amg_manager.bcf_reorder) || (F_false_C_true[ii] == false)) { // F nodes

							integer istr = ii - iadd;
							doublerealT rold = x[istr];

							// 13.07.2016
							doublerealT ap_now = 1.0 / Amat.aij[row_ptr_start[ii]];
							// 28.01.2017
							//ap_now = 0.0;

							x[istr] = diag_minus_one[istr] * b[istr];

							doublerealT rsum = 0.0;
							integer is1 = row_ptr_start[ii] + 1;
							integer is2 = row_ptr_end[ii];
							// ����������������� �������� �������� ����� ������.
							//#pragma omp parallel for reduction(+:rsum)
							for (integer ii1 = is1; ii1 <= is2; ++ii1)
							{
								//x[istr] += -Amat.aij[ii1]*x[Amat.j[ii1]];
								integer ipos = Amat.j[ii1];
								// 13.07.2016
								if (1) {
									// ������������� positive connections.
									//if ((Amat.aij[ii1] < 0.0)) {
									// ����� �� � ���� � ������� � ������� ������������.
									rsum += -Amat.aij[ii1] * x[ipos];
									//}
									//else {
									// �� ��������.
									//	ap_now += Amat.aij[ii1];
									//}
								}
								else {
									// ����� � ����������. �� �������� ����� �������: TKM, TKM1, TKM2.

									// ������������� positive connections.
									if ((Amat.aij[ii1] < 0.0)) {
										rsum += -Amat.aij[ii1] * x[ipos];
										//ap_now += fabs(Amat.aij[ii1]);
									}
									else {
										// �� ��������.
										// �������� ��-�� ���� ��� ��� ������� ������.
										//ap_now += fabs(Amat.aij[ii1]);
										if (fabs(x[ipos]) > fabs(x[istr])) {
											ap_now += fabs(Amat.aij[ii1]);
										}
										else rsum += -Amat.aij[ii1] * x[ipos];
									}
								}
							}
							x[istr] += rsum;
							//x[istr] *= Amat.aij[row_ptr_start[ii]];
							// 13.07.2016
							x[istr] /= ap_now;

							// ����������� ������ ����� ������� ���� ����� ��������� ������
							// �.�. ��������� �������� �� �������� ����������.
							//if (is1 <= is2) {
								// ������ ���� ��� �� ������� ������� ��������� ����������.
							x[istr] = omega * x[istr] + (1.0 - omega) * rold; // this is SOR
						//}
						}
					}

				}

				// �� � ���� ������ �� �����������. 3 ������ 2016.
				//isimmetricGS_switch = 1;

			}
			else {

				if (idirect == 1) {
					// ���������� �����.

					// ������� F ����� C.

	//---->#pragma omp parallel for
					for (integer ii = endpos; ii >= startpos; --ii) {
						if ((!my_amg_manager.bcf_reorder) || (F_false_C_true[ii] == false)) { // F nodes

							integer istr = ii - iadd;
							doublerealT rold = x[istr];

							// 13.07.2016
							doublerealT ap_now = 1.0 / Amat.aij[row_ptr_start[ii]];

							x[istr] = diag_minus_one[istr] * b[istr];

							doublerealT rsum = 0.0;
							integer is1 = row_ptr_start[ii] + 1;
							integer is2 = row_ptr_end[ii];
							// ����������������� �������� �������� ����� ������.
							//#pragma omp parallel for reduction(+:rsum)
							for (integer ii1 = is1; ii1 <= is2; ++ii1)
							{
								//x[istr] += -Amat.aij[ii1]*x[Amat.j[ii1]];
								integer ipos = Amat.j[ii1];
								// 13.07.2016
								// ������������� positive connections.
								//if ((Amat.aij[ii1] < 0.0)) {
								rsum += -Amat.aij[ii1] * x[ipos];
								//}
								//else {
								// �� �������.
								//	ap_now += Amat.aij[ii1];
								//}
							}
							x[istr] += rsum;
							//x[istr] *= Amat.aij[row_ptr_start[ii]];
							// 13.07.2016
							x[istr] /= ap_now;

							// ����������� ������ ����� ������� ���� ����� ��������� ������
							// �.�. ��������� �������� �� �������� ����������.
							//if (is1 <= is2) {
								// ������ ���� ��� �� ������� ������� ��������� ����������.
							x[istr] = omega * x[istr] + (1.0 - omega) * rold; // this is SOR
						//}
						}
					}


					//----->#pragma omp parallel for
					for (integer ii = endpos; ii >= startpos; --ii) {
						if ((!my_amg_manager.bcf_reorder) || (F_false_C_true[ii])) { // C nodes

							integer istr = ii - iadd;
							doublerealT rold = x[istr];

							// 13.07.2016
							doublerealT ap_now = 1.0 / Amat.aij[row_ptr_start[ii]];

							x[istr] = diag_minus_one[istr] * b[istr];

							doublerealT rsum = 0.0;
							integer is1 = row_ptr_start[ii] + 1;
							integer is2 = row_ptr_end[ii];
							// ����������������� �������� �������� ����� ������.
							//#pragma omp parallel for reduction(+:rsum)
							for (integer ii1 = is1; ii1 <= is2; ++ii1)
							{
								//x[istr] += -Amat.aij[ii1]*x[Amat.j[ii1]];
								integer ipos = Amat.j[ii1];
								// 13.07.2016
								// ������������� positive connections.
								//if ((Amat.aij[ii1] < 0.0)) {
								rsum += -Amat.aij[ii1] * x[ipos];
								//}
								//else {
								// �� �������.
								//	ap_now += Amat.aij[ii1];
								//}
							}
							x[istr] += rsum;
							//x[istr] *= Amat.aij[row_ptr_start[ii]];
							// 13.07.2016
							x[istr] /= ap_now;

							// ����������� ������ ����� ������� ���� ����� ��������� ������
							// �.�. ��������� �������� �� �������� ����������.
							//if (is1 <= is2) {
								// ������ ���� ��� �� ������� ������� ��������� ����������.
							x[istr] = omega * x[istr] + (1.0 - omega) * rold; // this is SOR
						//}
						}
					}



				}
				else {
					// idirect==0
					// ������� � ����� F.

	//--->#pragma omp parallel for
					for (integer ii = endpos; ii >= startpos; --ii) {
						if ((!my_amg_manager.bcf_reorder) || (F_false_C_true[ii])) { // C nodes

							integer istr = ii - iadd;
							doublerealT rold = x[istr];

							// 13.07.2016
							doublerealT ap_now = 1.0 / Amat.aij[row_ptr_start[ii]];

							x[istr] = diag_minus_one[istr] * b[istr];

							doublerealT rsum = 0.0;
							integer is1 = row_ptr_start[ii] + 1;
							integer is2 = row_ptr_end[ii];
							// ����������������� �������� �������� ����� ������.
							//#pragma omp parallel for reduction(+:rsum)
							for (integer ii1 = is1; ii1 <= is2; ++ii1)
							{
								//x[istr] += -Amat.aij[ii1]*x[Amat.j[ii1]];
								integer ipos = Amat.j[ii1];
								// 13.07.2016
								// ������������� positive connections.
								//if ((Amat.aij[ii1] < 0.0)) {
								rsum += -Amat.aij[ii1] * x[ipos];
								//}
								//else {
								// �� �������.
								//	ap_now += Amat.aij[ii1];
								//}
							}
							x[istr] += rsum;
							//x[istr] *= Amat.aij[row_ptr_start[ii]];
							// 13.07.2016
							x[istr] /= ap_now;

							// ����������� ������ ����� ������� ���� ����� ��������� ������
							// �.�. ��������� �������� �� �������� ����������.
							//if (is1 <= is2) {
								// ������ ���� ��� �� ������� ������� ��������� ����������.
							x[istr] = omega * x[istr] + (1.0 - omega) * rold; // this is SOR
						//}
						}
					}


					//--->#pragma omp parallel for
					for (integer ii = endpos; ii >= startpos; --ii) {
						if ((!my_amg_manager.bcf_reorder) || (F_false_C_true[ii] == false)) { // F nodes

							integer istr = ii - iadd;
							doublerealT rold = x[istr];

							// 13.07.2016
							doublerealT ap_now = 1.0 / Amat.aij[row_ptr_start[ii]];

							x[istr] = diag_minus_one[istr] * b[istr];

							doublerealT rsum = 0.0;
							integer is1 = row_ptr_start[ii] + 1;
							integer is2 = row_ptr_end[ii];
							// ����������������� �������� �������� ����� ������.
							//#pragma omp parallel for reduction(+:rsum)
							for (integer ii1 = is1; ii1 <= is2; ++ii1)
							{
								//x[istr] += -Amat.aij[ii1]*x[Amat.j[ii1]];
								integer ipos = Amat.j[ii1];
								// 13.07.2016
								// ������������� positive connections.
								//if ((Amat.aij[ii1] < 0.0)) {
								rsum += -Amat.aij[ii1] * x[ipos];
								//}
								//else {
								// �� �������.
								//	ap_now += Amat.aij[ii1];
								//}
							}
							x[istr] += rsum;
							//x[istr] *= Amat.aij[row_ptr_start[ii]];
							// 13.07.2016
							x[istr] /= ap_now;

							// ����������� ������ ����� ������� ���� ����� ��������� ������
							// �.�. ��������� �������� �� �������� ����������.
							//if (is1 <= is2) {
								// ������ ���� ��� �� ������� ������� ��������� ����������.
							x[istr] = omega * x[istr] + (1.0 - omega) * rold; // this is SOR
						//}
						}
					}

				}

				//isimmetricGS_switch = 0;  // ����� �����������.
			}
		}
	}


} // seidelqsor2


// smoother.
// 16 jan 2016.  Seidel q -quick SOR + parallel
// 9 september 2015.
// q - quick.
template <typename doublerealT>
void seidelq(Ak2& Amat, integer istartq, integer iendq, doublerealT*& x, doublerealT*& b, bool*& bnested_desection, integer*& row_ptr_start, integer*& row_ptr_end, integer iadd)
{
	seidelqsor2Pcpu<doublerealT>(Amat, istartq, iendq, x, b, bnested_desection, row_ptr_start, row_ptr_end, iadd);
	//seidelqsor2Pcpu(Amat, istartq, iendq, x, b, row_ptr_start, row_ptr_end, iadd);

} // seidelq

// smoother.
// 16 jan 2016.  Seidel q -quick SOR + parallel
// 9 september 2015.
// q - quick.
template <typename doublerealT>
void seidelq(Ak2& Amat, integer istartq, integer iendq, doublerealT*& x, doublerealT*& b, bool*& bnested_desection, integer*& row_ptr_start, integer*& row_ptr_end, integer iadd, doublerealT*& diag_minus_one)
{
	seidelqsor2Pcpu<doublerealT>(Amat, istartq, iendq, x, b, bnested_desection, row_ptr_start, row_ptr_end, iadd, diag_minus_one);
	//seidelqsor2Pcpu(Amat, istartq, iendq, x, b, row_ptr_start, row_ptr_end, iadd);

} // seidelq

  // smoother.
  // 16 jan 2016.  Seidel q -quick SOR + parallel
  // 9 september 2015.
  // q - quick.
template <typename doublerealT>
void seidelq(Ak2& Amat, integer istartq, integer iendq, doublerealT*& x, doublerealT*& b, bool*& bnested_desection, integer*& row_ptr_start, integer*& row_ptr_end, integer iadd, bool*& F_false_C_true, integer idirect)
{
	// , bool* &F_false_C_true, integer idirect ��������, ��������� �� ������������.
	// �������� �������� �������������.

	seidelqsor2Pcpu<doublerealT>(Amat, istartq, iendq, x, b, bnested_desection, row_ptr_start, row_ptr_end, iadd);
	//seidelqsor2Pcpu<doublerealT>(Amat, istartq, iendq, x, b, row_ptr_start, row_ptr_end, iadd);

} // seidelq


// smoother.
  // 16 jan 2016.  Seidel q -quick SOR + parallel
  // 9 september 2015.
  // q - quick.
template <typename doublerealT>
void seidelq(Ak2& Amat, integer istartq, integer iendq, doublerealT*& x, doublerealT*& b,
	bool*& bnested_desection, integer*& row_ptr_start, integer*& row_ptr_end, integer iadd,
	bool*& F_false_C_true, integer idirect, doublerealT*& diag_minus_one)
{
	// , bool* &F_false_C_true, integer idirect ��������, ��������� �� ������������.
	// �������� �������� �������������.

	seidelqsor2Pcpu<doublerealT>(Amat, istartq, iendq, x, b, bnested_desection, row_ptr_start, row_ptr_end, iadd, diag_minus_one);
	//seidelqsor2Pcpu<doublerealT>(Amat, istartq, iendq, x, b, row_ptr_start, row_ptr_end, iadd);

} // seidelq


#endif 