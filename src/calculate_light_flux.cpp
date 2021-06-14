// Файл calculate_light_flux.cpp содержит функции необходимые для расчёта освещенности,
// с помощью алгоритма трассировки лучей.

#pragma once
#ifndef CALCULATE_LIGHT_FLUX_CPP
#define CALCULATE_LIGHT_FLUX_CPP 1


// Расстояние от точки до отрезка.
// p - точка. Отрезок [sp0,sp1]. s - segment.
doublereal distance_Point_to_Segment(TOCHKA p, TOCHKA sp0, TOCHKA sp1) {
	doublereal t = ((p.x - sp0.x) * (sp1.x - sp0.x) + (p.y - sp0.y) * (sp1.y - sp0.y) + (p.z - sp0.z) * (sp1.z - sp0.z)) /
		((sp1.x - sp0.x)* (sp1.x - sp0.x) + (sp1.y - sp0.y) * (sp1.y - sp0.y) + (sp1.z - sp0.z) * (sp1.z - sp0.z));
	if (t < 0) t = 0.0;
	if (t > 1.0) t = 1.0;
	return sqrt((sp0.x - p.x + (sp1.x - sp0.x)*t) * (sp0.x - p.x + (sp1.x - sp0.x) * t) +
		(sp0.y - p.y + (sp1.y - sp0.y) * t) * (sp0.y - p.y + (sp1.y - sp0.y) * t) +
		(sp0.z - p.z + (sp1.z - sp0.z) * t) * (sp0.z - p.z + (sp1.z - sp0.z) * t));
} // distance_Point_to_Segment

  // Квадрат расстояния от точки до отрезка.
  // p - точка. Отрезок [sp0,sp1]. s - segment.
doublereal distance_Point_to_Segment2(TOCHKA p, TOCHKA sp0, TOCHKA sp1) {

	doublereal lx = (sp1.x - sp0.x);
	doublereal ly = (sp1.y - sp0.y);
	doublereal lz = (sp1.z - sp0.z);

	doublereal t = ((p.x - sp0.x) * lx + (p.y - sp0.y) * ly + (p.z - sp0.z) * lz) /
		(lx * lx + ly * ly + lz * lz);
	if (t < 0) t = 0.0;
	if (t > 1.0) t = 1.0;
	doublereal kx = (sp0.x - p.x + lx * t);
	doublereal ky = (sp0.y - p.y + ly * t);
	doublereal kz = (sp0.z - p.z + lz * t);
	return (kx * kx + ky * ky + kz * kz);
} // distance_Point_to_Segment

  // Квадрат расстояния от точки до отрезка.
  // p - точка. Отрезок [sp0,sp1]. s - segment.
doublereal distance_Point_to_Segment3(TOCHKA p, TOCHKA sp0, TOCHKA sp1,
	doublereal lx, doublereal ly, doublereal lz, doublereal ld2) {

	doublereal t = ((p.x - sp0.x) * lx + (p.y - sp0.y) * ly + (p.z - sp0.z) * lz) / (ld2);
	if (t < 0) t = 0.0;
	if (t > 1.0) t = 1.0;
	doublereal kx = (sp0.x - p.x + lx * t);
	doublereal ky = (sp0.y - p.y + ly * t);
	doublereal kz = (sp0.z - p.z + lz * t);
	return (kx * kx + ky * ky + kz * kz);
} // distance_Point_to_Segment


void calculate_light_flux2(doublereal* & myF, TEMPER& t,
	BLOCK*& b, integer &lb) 
{
	doublereal* distance = new doublereal[t.maxelm + t.maxbound];
	for (integer i = 0; i < t.maxelm; i++) {
		myF[i] = 0.0;
	}
	for (integer ibid = 0; ibid < lb; ibid++) {
		for (integer i = 0; i < t.maxelm; i++) {
			distance[i] = 1.0e30;
			if (b[t.whot_is_block[i]].itype == PHYSICS_TYPE_IN_BODY::SOLID) {
				distance[i] = -0.1e-3;// одна десятая мм
			}
		}

		if (b[ibid].arr_Sc[0] > 1.0e-30) {
			printf("power %e [W], ", b[ibid].arr_Sc[0]);
			doublereal x0 = 0.5 * (b[ibid].g.xS + b[ibid].g.xE);
			doublereal y0 = 0.5 * (b[ibid].g.yS + b[ibid].g.yE);
			doublereal z0 = 0.5 * (b[ibid].g.zS + b[ibid].g.zE);
			printf("x0=%e [m], y0=%e [m], z0=%e [m], maxelm=%lld\n",x0,y0,z0, t.maxelm);
			for (integer i = 0; i < t.maxelm; i++) {
				TOCHKA p1;
				center_cord3D(i, t.nvtx, t.pa, p1, 100);
				doublereal dist0 = sqrt((p1.x - x0) * (p1.x - x0) + (p1.y - y0) * (p1.y - y0) + (p1.z - z0) * (p1.z - z0));
				if (dist0 < distance[i]) {
					distance[i] = dist0;
				}
			}
			for (integer i = 0; i < t.maxelm; i++) {
				if (distance[i] < 0.0) {
					distance[i] = 1.0e10;
				}
			}
			for (integer i = 0; i < t.maxelm; i++) {
				myF[i] += 1.0 / (M_PI * distance[i] * distance[i]);
			}
		}
	}
	delete[] distance;

} // calculate_light_flux2



bool permissible758 = true;
integer* stack_ray_trayser = nullptr;
integer itop_stack_ray_trayser = -1;
bool* visit_ray_traser = nullptr;
doublereal* size_cell = nullptr;
const doublereal kR_size_cell = 0.7; //0,49; 1.44;
bool bray_tracing_fluid_only = false;

// p0, i0 - diod, p1, i1 - point.
// p3 - центр отрезка, radius2 - квадрат половины длины отрезка.
void ray_tracing(TOCHKA p3, TOCHKA p0, TOCHKA p1, integer i0, integer i1, TEMPER& t, BLOCK*& b, doublereal radius2) {


	doublereal lx = (p1.x - p0.x);
	doublereal ly = (p1.y - p0.y);
	doublereal lz = (p1.z - p0.z);
	doublereal ld2 = (lx * lx + ly * ly + lz * lz);


	visit_ray_traser = new bool[t.maxelm];
	for (integer i = 0; i < t.maxelm; i++) visit_ray_traser[i] = false;

	while ((itop_stack_ray_trayser > -1) && (permissible758)) {
		integer iP = stack_ray_trayser[itop_stack_ray_trayser];
		itop_stack_ray_trayser--;

		if (permissible758) {


			// вызов.

			if ((t.neighbors_for_the_internal_node[E_SIDE][0] != nullptr)&&
				(t.neighbors_for_the_internal_node[E_SIDE][0][iP] < t.maxelm)) {
				integer i = t.neighbors_for_the_internal_node[E_SIDE][0][iP];

				if (visit_ray_traser[i] == false) {

					TOCHKA p2;
					center_cord3D_ray_tracing(i, t.nvtx, t.pa, p2, 100);


					if (distance_Point_to_Segment3(p2, p0, p1, lx, ly, lz, ld2) < size_cell[i])
					{

						// точка p2 принадлежит отрезку [p0,p1].
						if (!(((t.whot_is_block[i] == t.whot_is_block[i0]) ||
							(b[t.whot_is_block[i]].imatid == b[t.whot_is_block[i1]].imatid)))) {
							permissible758 = false; // не достижима.
						}
						else {
							// Продолжаем сканирование.
							itop_stack_ray_trayser++;
							stack_ray_trayser[itop_stack_ray_trayser] = i;
							visit_ray_traser[i] = true;
						}
					}

				}


			}

		}
		if (permissible758) {


			if ((t.neighbors_for_the_internal_node[W_SIDE][0] != nullptr)&&
				(t.neighbors_for_the_internal_node[W_SIDE][0][iP] < t.maxelm)) {
				integer i = t.neighbors_for_the_internal_node[W_SIDE][0][iP];

				if (visit_ray_traser[i] == false) {

					TOCHKA p2;
					center_cord3D_ray_tracing(i, t.nvtx, t.pa, p2, 100);

					if (distance_Point_to_Segment3(p2, p0, p1, lx, ly, lz, ld2) < size_cell[i])
					{

						// точка p2 принадлежит отрезку [p0,p1].
						if (!(((t.whot_is_block[i] == t.whot_is_block[i0]) ||
							(b[t.whot_is_block[i]].imatid == b[t.whot_is_block[i1]].imatid)))) {
							permissible758 = false; // не достижима.
						}
						else {
							// Продолжаем сканирование.
							itop_stack_ray_trayser++;
							stack_ray_trayser[itop_stack_ray_trayser] = i;
							visit_ray_traser[i] = true;
						}
					}

				}
			}
		}
		if (permissible758) {


			if ((t.neighbors_for_the_internal_node[N_SIDE][0] != nullptr) &&
				(t.neighbors_for_the_internal_node[N_SIDE][0][iP] < t.maxelm)) {
				integer i = t.neighbors_for_the_internal_node[N_SIDE][0][iP];

				if (visit_ray_traser[i] == false) {

					TOCHKA p2;
					center_cord3D_ray_tracing(i, t.nvtx, t.pa, p2, 100);


					if (distance_Point_to_Segment3(p2, p0, p1, lx, ly, lz, ld2) < size_cell[i])
					{

						// точка p2 принадлежит отрезку [p0,p1].
						if (!(((t.whot_is_block[i] == t.whot_is_block[i0]) ||
							(b[t.whot_is_block[i]].imatid == b[t.whot_is_block[i1]].imatid)))) {
							permissible758 = false; // не достижима.
						}
						else {
							// Продолжаем сканирование.
							itop_stack_ray_trayser++;
							stack_ray_trayser[itop_stack_ray_trayser] = i;
							visit_ray_traser[i] = true;
						}
					}

				}
			}
		}
		if (permissible758) {


			if ((t.neighbors_for_the_internal_node[S_SIDE][0] != nullptr)&&
				(t.neighbors_for_the_internal_node[S_SIDE][0][iP] < t.maxelm)) {
				integer i = t.neighbors_for_the_internal_node[S_SIDE][0][iP];

				if (visit_ray_traser[i] == false) {

					TOCHKA p2;
					center_cord3D_ray_tracing(i, t.nvtx, t.pa, p2, 100);

					if (distance_Point_to_Segment3(p2, p0, p1, lx, ly, lz, ld2) < size_cell[i])
					{

						// точка p2 принадлежит отрезку [p0,p1].
						if (!(((t.whot_is_block[i] == t.whot_is_block[i0]) ||
							(b[t.whot_is_block[i]].imatid == b[t.whot_is_block[i1]].imatid)))) {
							permissible758 = false; // не достижима.
						}
						else {
							// Продолжаем сканирование.
							itop_stack_ray_trayser++;
							stack_ray_trayser[itop_stack_ray_trayser] = i;
							visit_ray_traser[i] = true;
						}
					}

				}
			}
		}

		if (permissible758) {

			if ((t.neighbors_for_the_internal_node[T_SIDE][0] != nullptr)&&
				(t.neighbors_for_the_internal_node[T_SIDE][0][iP] < t.maxelm)) {
				integer i = t.neighbors_for_the_internal_node[T_SIDE][0][iP];

				if (visit_ray_traser[i] == false) {

					TOCHKA p2;
					center_cord3D_ray_tracing(i, t.nvtx, t.pa, p2, 100);

					if (distance_Point_to_Segment3(p2, p0, p1, lx, ly, lz, ld2) < size_cell[i])
					{
						// точка p2 принадлежит отрезку [p0,p1].
						if (!(((t.whot_is_block[i] == t.whot_is_block[i0]) ||
							(b[t.whot_is_block[i]].imatid == b[t.whot_is_block[i1]].imatid)))) {
							permissible758 = false; // не достижима.
						}
						else {
							// Продолжаем сканирование.
							itop_stack_ray_trayser++;
							stack_ray_trayser[itop_stack_ray_trayser] = i;
							visit_ray_traser[i] = true;

						}
					}

				}
			}
		}
		if (permissible758) {


			if ((t.neighbors_for_the_internal_node[B_SIDE][0] != nullptr)&&
				(t.neighbors_for_the_internal_node[B_SIDE][0][iP] < t.maxelm)) {
				integer i = t.neighbors_for_the_internal_node[B_SIDE][0][iP];

				if (visit_ray_traser[i] == false) {

					TOCHKA p2;
					center_cord3D_ray_tracing(i, t.nvtx, t.pa, p2, 100);

					if (distance_Point_to_Segment3(p2, p0, p1, lx, ly, lz, ld2) < size_cell[i])
					{
						// точка p2 принадлежит отрезку [p0,p1].
						if (!(((t.whot_is_block[i] == t.whot_is_block[i0]) ||
							(b[t.whot_is_block[i]].imatid == b[t.whot_is_block[i1]].imatid)))) {
							permissible758 = false; // не достижима.
						}
						else {
							// Продолжаем сканирование.
							itop_stack_ray_trayser++;
							stack_ray_trayser[itop_stack_ray_trayser] = i;
							visit_ray_traser[i] = true;
						}
					}
				}
			}
		}
	}


	if (visit_ray_traser[i1] == false) {
		// Точка i1 не была достигнута, путь преградил hollow block.
		permissible758 = false; // не достижимо.
	}

	delete[] visit_ray_traser;
} // ray_tracing

  // p0, i0 - diod, p1, i1 - point.
  // p3 - центр отрезка, radius2 - квадрат половины длины отрезка.
  // работает на АЛИС сетке.
void ray_tracing_Alice(TOCHKA p3, TOCHKA p0, TOCHKA p1, integer i0, integer i1,
	TEMPER& t, BLOCK*& b, doublereal radius2) {


	doublereal lx = (p1.x - p0.x);
	doublereal ly = (p1.y - p0.y);
	doublereal lz = (p1.z - p0.z);
	doublereal ld2 = (lx * lx + ly * ly + lz * lz);


	visit_ray_traser = new bool[t.maxelm];
	for (integer i = 0; i < t.maxelm; i++) visit_ray_traser[i] = false;

	while ((itop_stack_ray_trayser > -1) && (permissible758)) {
		integer iP = stack_ray_trayser[itop_stack_ray_trayser];
		itop_stack_ray_trayser--;

		if (permissible758) {

			// вызов.

			if ((t.neighbors_for_the_internal_node[E_SIDE][0] != nullptr) &&
				(t.neighbors_for_the_internal_node[E_SIDE][0][iP] > -1) &&
				(t.neighbors_for_the_internal_node[E_SIDE][0][iP] < t.maxelm)) {
				integer i = t.neighbors_for_the_internal_node[E_SIDE][0][iP];

				if (visit_ray_traser[i] == false) {

					TOCHKA p2;
					center_cord3D_ray_tracing(i, t.nvtx, t.pa, p2, 100);


					if (distance_Point_to_Segment3(p2, p0, p1, lx, ly, lz, ld2) < size_cell[i])
					{

						// точка p2 принадлежит отрезку [p0,p1].
						if (!(((t.whot_is_block[i] == t.whot_is_block[i0]) ||
							(b[t.whot_is_block[i]].imatid == b[t.whot_is_block[i1]].imatid)))) {
							permissible758 = false; // не достижима.
						}
						else {
							// Продолжаем сканирование.
							itop_stack_ray_trayser++;
							stack_ray_trayser[itop_stack_ray_trayser] = i;
							visit_ray_traser[i] = true;
						}
					}

				}

			}

		}

		if (permissible758) {

			// вызов.

			if ((t.neighbors_for_the_internal_node[E_SIDE][1] != nullptr) &&
				(t.neighbors_for_the_internal_node[E_SIDE][1][iP] > -1) &&
				(t.neighbors_for_the_internal_node[E_SIDE][1][iP] < t.maxelm)) {
				integer i = t.neighbors_for_the_internal_node[E_SIDE][1][iP];

				if (visit_ray_traser[i] == false) {

					TOCHKA p2;
					center_cord3D_ray_tracing(i, t.nvtx, t.pa, p2, 100);


					if (distance_Point_to_Segment3(p2, p0, p1, lx, ly, lz, ld2) < size_cell[i])
					{

						// точка p2 принадлежит отрезку [p0,p1].
						if (!(((t.whot_is_block[i] == t.whot_is_block[i0]) ||
							(b[t.whot_is_block[i]].imatid == b[t.whot_is_block[i1]].imatid)))) {
							permissible758 = false; // не достижима.
						}
						else {
							// Продолжаем сканирование.
							itop_stack_ray_trayser++;
							stack_ray_trayser[itop_stack_ray_trayser] = i;
							visit_ray_traser[i] = true;
						}
					}

				}

			}

		}

		if (permissible758) {

			// вызов.

			if ((t.neighbors_for_the_internal_node[E_SIDE][2] != nullptr) &&
				(t.neighbors_for_the_internal_node[E_SIDE][2][iP] > -1) &&
				(t.neighbors_for_the_internal_node[E_SIDE][2][iP] < t.maxelm)) {
				integer i = t.neighbors_for_the_internal_node[E_SIDE][2][iP];

				if (visit_ray_traser[i] == false) {

					TOCHKA p2;
					center_cord3D_ray_tracing(i, t.nvtx, t.pa, p2, 100);


					if (distance_Point_to_Segment3(p2, p0, p1, lx, ly, lz, ld2) < size_cell[i])
					{

						// точка p2 принадлежит отрезку [p0,p1].
						if (!(((t.whot_is_block[i] == t.whot_is_block[i0]) ||
							(b[t.whot_is_block[i]].imatid == b[t.whot_is_block[i1]].imatid)))) {
							permissible758 = false; // не достижима.
						}
						else {
							// Продолжаем сканирование.
							itop_stack_ray_trayser++;
							stack_ray_trayser[itop_stack_ray_trayser] = i;
							visit_ray_traser[i] = true;
						}
					}

				}

			}

		}

		if (permissible758) {

			// вызов.

			if ((t.neighbors_for_the_internal_node[E_SIDE][3] != nullptr) &&
				(t.neighbors_for_the_internal_node[E_SIDE][3][iP] > -1) &&
				(t.neighbors_for_the_internal_node[E_SIDE][3][iP] < t.maxelm)) {
				integer i = t.neighbors_for_the_internal_node[E_SIDE][3][iP];

				if (visit_ray_traser[i] == false) {

					TOCHKA p2;
					center_cord3D_ray_tracing(i, t.nvtx, t.pa, p2, 100);


					if (distance_Point_to_Segment3(p2, p0, p1, lx, ly, lz, ld2) < size_cell[i])
					{

						// точка p2 принадлежит отрезку [p0,p1].
						if (!(((t.whot_is_block[i] == t.whot_is_block[i0]) ||
							(b[t.whot_is_block[i]].imatid == b[t.whot_is_block[i1]].imatid)))) {
							permissible758 = false; // не достижима.
						}
						else {
							// Продолжаем сканирование.
							itop_stack_ray_trayser++;
							stack_ray_trayser[itop_stack_ray_trayser] = i;
							visit_ray_traser[i] = true;
						}
					}

				}

			}

		}

		if (permissible758) {


			if ((t.neighbors_for_the_internal_node[W_SIDE][0] != nullptr) &&
				(t.neighbors_for_the_internal_node[W_SIDE][0][iP] > -1) &&
				(t.neighbors_for_the_internal_node[W_SIDE][0][iP] < t.maxelm)) {
				integer i = t.neighbors_for_the_internal_node[W_SIDE][0][iP];

				if (visit_ray_traser[i] == false) {

					TOCHKA p2;
					center_cord3D_ray_tracing(i, t.nvtx, t.pa, p2, 100);

					if (distance_Point_to_Segment3(p2, p0, p1, lx, ly, lz, ld2) < size_cell[i])
					{

						// точка p2 принадлежит отрезку [p0,p1].
						if (!(((t.whot_is_block[i] == t.whot_is_block[i0]) ||
							(b[t.whot_is_block[i]].imatid == b[t.whot_is_block[i1]].imatid)))) {
							permissible758 = false; // не достижима.
						}
						else {
							// Продолжаем сканирование.
							itop_stack_ray_trayser++;
							stack_ray_trayser[itop_stack_ray_trayser] = i;
							visit_ray_traser[i] = true;
						}
					}

				}
			}
		}

		if (permissible758) {


			if ((t.neighbors_for_the_internal_node[W_SIDE][1] != nullptr) &&
				(t.neighbors_for_the_internal_node[W_SIDE][1][iP] > -1) &&
				(t.neighbors_for_the_internal_node[W_SIDE][1][iP] < t.maxelm)) {
				integer i = t.neighbors_for_the_internal_node[W_SIDE][1][iP];

				if (visit_ray_traser[i] == false) {

					TOCHKA p2;
					center_cord3D_ray_tracing(i, t.nvtx, t.pa, p2, 100);

					if (distance_Point_to_Segment3(p2, p0, p1, lx, ly, lz, ld2) < size_cell[i])
					{

						// точка p2 принадлежит отрезку [p0,p1].
						if (!(((t.whot_is_block[i] == t.whot_is_block[i0]) ||
							(b[t.whot_is_block[i]].imatid == b[t.whot_is_block[i1]].imatid)))) {
							permissible758 = false; // не достижима.
						}
						else {
							// Продолжаем сканирование.
							itop_stack_ray_trayser++;
							stack_ray_trayser[itop_stack_ray_trayser] = i;
							visit_ray_traser[i] = true;
						}
					}

				}
			}
		}


		if (permissible758) {


			if ((t.neighbors_for_the_internal_node[W_SIDE][2] != nullptr) && 
				(t.neighbors_for_the_internal_node[W_SIDE][2][iP] > -1) &&
				(t.neighbors_for_the_internal_node[W_SIDE][2][iP] < t.maxelm)) {
				integer i = t.neighbors_for_the_internal_node[W_SIDE][2][iP];

				if (visit_ray_traser[i] == false) {

					TOCHKA p2;
					center_cord3D_ray_tracing(i, t.nvtx, t.pa, p2, 100);

					if (distance_Point_to_Segment3(p2, p0, p1, lx, ly, lz, ld2) < size_cell[i])
					{

						// точка p2 принадлежит отрезку [p0,p1].
						if (!(((t.whot_is_block[i] == t.whot_is_block[i0]) ||
							(b[t.whot_is_block[i]].imatid == b[t.whot_is_block[i1]].imatid)))) {
							permissible758 = false; // не достижима.
						}
						else {
							// Продолжаем сканирование.
							itop_stack_ray_trayser++;
							stack_ray_trayser[itop_stack_ray_trayser] = i;
							visit_ray_traser[i] = true;
						}
					}
				}
			}
		}

		if (permissible758) {


			if ((t.neighbors_for_the_internal_node[W_SIDE][3] != nullptr) &&
				(t.neighbors_for_the_internal_node[W_SIDE][3][iP] > -1) &&
				(t.neighbors_for_the_internal_node[W_SIDE][3][iP] < t.maxelm)) {
				integer i = t.neighbors_for_the_internal_node[W_SIDE][3][iP];

				if (visit_ray_traser[i] == false) {

					TOCHKA p2;
					center_cord3D_ray_tracing(i, t.nvtx, t.pa, p2, 100);

					if (distance_Point_to_Segment3(p2, p0, p1, lx, ly, lz, ld2) < size_cell[i])
					{

						// точка p2 принадлежит отрезку [p0,p1].
						if (!(((t.whot_is_block[i] == t.whot_is_block[i0]) ||
							(b[t.whot_is_block[i]].imatid == b[t.whot_is_block[i1]].imatid)))) {
							permissible758 = false; // не достижима.
						}
						else {
							// Продолжаем сканирование.
							itop_stack_ray_trayser++;
							stack_ray_trayser[itop_stack_ray_trayser] = i;
							visit_ray_traser[i] = true;
						}
					}

				}
			}
		}

		if (permissible758) {


			if ((t.neighbors_for_the_internal_node[N_SIDE][0] != nullptr) &&
				(t.neighbors_for_the_internal_node[N_SIDE][0][iP] > -1)&&
				(t.neighbors_for_the_internal_node[N_SIDE][0][iP] < t.maxelm)) {
				integer i = t.neighbors_for_the_internal_node[N_SIDE][0][iP];

				if (visit_ray_traser[i] == false) {

					TOCHKA p2;
					center_cord3D_ray_tracing(i, t.nvtx, t.pa, p2, 100);


					if (distance_Point_to_Segment3(p2, p0, p1, lx, ly, lz, ld2) < size_cell[i])
					{

						// точка p2 принадлежит отрезку [p0,p1].
						if (!(((t.whot_is_block[i] == t.whot_is_block[i0]) ||
							(b[t.whot_is_block[i]].imatid == b[t.whot_is_block[i1]].imatid)))) {
							permissible758 = false; // не достижима.
						}
						else {
							// Продолжаем сканирование.
							itop_stack_ray_trayser++;
							stack_ray_trayser[itop_stack_ray_trayser] = i;
							visit_ray_traser[i] = true;
						}
					}
				}
			}
		}

		if (permissible758) {


			if ((t.neighbors_for_the_internal_node[N_SIDE][1] != nullptr) && 
				(t.neighbors_for_the_internal_node[N_SIDE][1][iP] > -1) &&
				(t.neighbors_for_the_internal_node[N_SIDE][1][iP] < t.maxelm)) {
				integer i = t.neighbors_for_the_internal_node[N_SIDE][1][iP];

				if (visit_ray_traser[i] == false) {

					TOCHKA p2;
					center_cord3D_ray_tracing(i, t.nvtx, t.pa, p2, 100);


					if (distance_Point_to_Segment3(p2, p0, p1, lx, ly, lz, ld2) < size_cell[i])
					{

						// точка p2 принадлежит отрезку [p0,p1].
						if (!(((t.whot_is_block[i] == t.whot_is_block[i0]) ||
							(b[t.whot_is_block[i]].imatid == b[t.whot_is_block[i1]].imatid)))) {
							permissible758 = false; // не достижима.
						}
						else {
							// Продолжаем сканирование.
							itop_stack_ray_trayser++;
							stack_ray_trayser[itop_stack_ray_trayser] = i;
							visit_ray_traser[i] = true;
						}
					}
				}
			}
		}

		if (permissible758) {


			if ((t.neighbors_for_the_internal_node[N_SIDE][2] != nullptr) && 
				(t.neighbors_for_the_internal_node[N_SIDE][2][iP] > -1) &&
				(t.neighbors_for_the_internal_node[N_SIDE][2][iP] < t.maxelm)) {
				integer i = t.neighbors_for_the_internal_node[N_SIDE][2][iP];

				if (visit_ray_traser[i] == false) {

					TOCHKA p2;
					center_cord3D_ray_tracing(i, t.nvtx, t.pa, p2, 100);


					if (distance_Point_to_Segment3(p2, p0, p1, lx, ly, lz, ld2) < size_cell[i])
					{

						// точка p2 принадлежит отрезку [p0,p1].
						if (!(((t.whot_is_block[i] == t.whot_is_block[i0]) ||
							(b[t.whot_is_block[i]].imatid == b[t.whot_is_block[i1]].imatid)))) {
							permissible758 = false; // не достижима.
						}
						else {
							// Продолжаем сканирование.
							itop_stack_ray_trayser++;
							stack_ray_trayser[itop_stack_ray_trayser] = i;
							visit_ray_traser[i] = true;
						}
					}
				}
			}
		}

		if (permissible758) {


			if ((t.neighbors_for_the_internal_node[N_SIDE][3] != nullptr) &&
				(t.neighbors_for_the_internal_node[N_SIDE][3][iP] > -1) &&
				(t.neighbors_for_the_internal_node[N_SIDE][3][iP] < t.maxelm)) {
				integer i = t.neighbors_for_the_internal_node[N_SIDE][3][iP];

				if (visit_ray_traser[i] == false) {

					TOCHKA p2;
					center_cord3D_ray_tracing(i, t.nvtx, t.pa, p2, 100);


					if (distance_Point_to_Segment3(p2, p0, p1, lx, ly, lz, ld2) < size_cell[i])
					{

						// точка p2 принадлежит отрезку [p0,p1].
						if (!(((t.whot_is_block[i] == t.whot_is_block[i0]) ||
							(b[t.whot_is_block[i]].imatid == b[t.whot_is_block[i1]].imatid)))) {
							permissible758 = false; // не достижима.
						}
						else {
							// Продолжаем сканирование.
							itop_stack_ray_trayser++;
							stack_ray_trayser[itop_stack_ray_trayser] = i;
							visit_ray_traser[i] = true;
						}
					}
				}
			}
		}

		if (permissible758) {


			if ((t.neighbors_for_the_internal_node[S_SIDE][0] != nullptr) && 
				(t.neighbors_for_the_internal_node[S_SIDE][0][iP] > -1) &&
				(t.neighbors_for_the_internal_node[S_SIDE][0][iP] < t.maxelm)) {
				integer i = t.neighbors_for_the_internal_node[S_SIDE][0][iP];

				if (visit_ray_traser[i] == false) {

					TOCHKA p2;
					center_cord3D_ray_tracing(i, t.nvtx, t.pa, p2, 100);

					if (distance_Point_to_Segment3(p2, p0, p1, lx, ly, lz, ld2) < size_cell[i])
					{

						// точка p2 принадлежит отрезку [p0,p1].
						if (!(((t.whot_is_block[i] == t.whot_is_block[i0]) ||
							(b[t.whot_is_block[i]].imatid == b[t.whot_is_block[i1]].imatid)))) {
							permissible758 = false; // не достижима.
						}
						else {
							// Продолжаем сканирование.
							itop_stack_ray_trayser++;
							stack_ray_trayser[itop_stack_ray_trayser] = i;
							visit_ray_traser[i] = true;
						}
					}
				}
			}
		}

		if (permissible758) {


			if ((t.neighbors_for_the_internal_node[S_SIDE][1] != nullptr) &&
				(t.neighbors_for_the_internal_node[S_SIDE][1][iP] > -1) &&
				(t.neighbors_for_the_internal_node[S_SIDE][1][iP] < t.maxelm)) {
				integer i = t.neighbors_for_the_internal_node[S_SIDE][1][iP];

				if (visit_ray_traser[i] == false) {

					TOCHKA p2;
					center_cord3D_ray_tracing(i, t.nvtx, t.pa, p2, 100);

					if (distance_Point_to_Segment3(p2, p0, p1, lx, ly, lz, ld2) < size_cell[i])
					{

						// точка p2 принадлежит отрезку [p0,p1].
						if (!(((t.whot_is_block[i] == t.whot_is_block[i0]) ||
							(b[t.whot_is_block[i]].imatid == b[t.whot_is_block[i1]].imatid)))) {
							permissible758 = false; // не достижима.
						}
						else {
							// Продолжаем сканирование.
							itop_stack_ray_trayser++;
							stack_ray_trayser[itop_stack_ray_trayser] = i;
							visit_ray_traser[i] = true;
						}
					}
				}
			}
		}

		if (permissible758) {


			if ((t.neighbors_for_the_internal_node[S_SIDE][2] != nullptr) &&
				(t.neighbors_for_the_internal_node[S_SIDE][2][iP] > -1) &&
				(t.neighbors_for_the_internal_node[S_SIDE][2][iP] < t.maxelm)) {
				integer i = t.neighbors_for_the_internal_node[S_SIDE][2][iP];

				if (visit_ray_traser[i] == false) {

					TOCHKA p2;
					center_cord3D_ray_tracing(i, t.nvtx, t.pa, p2, 100);

					if (distance_Point_to_Segment3(p2, p0, p1, lx, ly, lz, ld2) < size_cell[i])
					{

						// точка p2 принадлежит отрезку [p0,p1].
						if (!(((t.whot_is_block[i] == t.whot_is_block[i0]) ||
							(b[t.whot_is_block[i]].imatid == b[t.whot_is_block[i1]].imatid)))) {
							permissible758 = false; // не достижима.
						}
						else {
							// Продолжаем сканирование.
							itop_stack_ray_trayser++;
							stack_ray_trayser[itop_stack_ray_trayser] = i;
							visit_ray_traser[i] = true;
						}
					}
				}
			}
		}

		if (permissible758) {


			if ((t.neighbors_for_the_internal_node[S_SIDE][3] != nullptr) && 
				(t.neighbors_for_the_internal_node[S_SIDE][3][iP] > -1) &&
				(t.neighbors_for_the_internal_node[S_SIDE][3][iP] < t.maxelm)) {
				integer i = t.neighbors_for_the_internal_node[S_SIDE][3][iP];

				if (visit_ray_traser[i] == false) {

					TOCHKA p2;
					center_cord3D_ray_tracing(i, t.nvtx, t.pa, p2, 100);

					if (distance_Point_to_Segment3(p2, p0, p1, lx, ly, lz, ld2) < size_cell[i])
					{

						// точка p2 принадлежит отрезку [p0,p1].
						if (!(((t.whot_is_block[i] == t.whot_is_block[i0]) ||
							(b[t.whot_is_block[i]].imatid == b[t.whot_is_block[i1]].imatid)))) {
							permissible758 = false; // не достижима.
						}
						else {
							// Продолжаем сканирование.
							itop_stack_ray_trayser++;
							stack_ray_trayser[itop_stack_ray_trayser] = i;
							visit_ray_traser[i] = true;
						}
					}
				}
			}
		}

		if (permissible758) {

			if ((t.neighbors_for_the_internal_node[T_SIDE][0] != nullptr)&&
				(t.neighbors_for_the_internal_node[T_SIDE][0][iP] > -1) &&
				(t.neighbors_for_the_internal_node[T_SIDE][0][iP] < t.maxelm)) {
				integer i = t.neighbors_for_the_internal_node[T_SIDE][0][iP];

				if (visit_ray_traser[i] == false) {

					TOCHKA p2;
					center_cord3D_ray_tracing(i, t.nvtx, t.pa, p2, 100);

					if (distance_Point_to_Segment3(p2, p0, p1, lx, ly, lz, ld2) < size_cell[i])
					{
						// точка p2 принадлежит отрезку [p0,p1].
						if (!(((t.whot_is_block[i] == t.whot_is_block[i0]) ||
							(b[t.whot_is_block[i]].imatid == b[t.whot_is_block[i1]].imatid)))) {
							permissible758 = false; // не достижима.
						}
						else {
							// Продолжаем сканирование.
							itop_stack_ray_trayser++;
							stack_ray_trayser[itop_stack_ray_trayser] = i;
							visit_ray_traser[i] = true;

						}
					}

				}
			}
		}

		if (permissible758) {

			if ((t.neighbors_for_the_internal_node[T_SIDE][1] != nullptr) &&
				(t.neighbors_for_the_internal_node[T_SIDE][1][iP] > -1) &&
				(t.neighbors_for_the_internal_node[T_SIDE][1][iP] < t.maxelm)) {
				integer i = t.neighbors_for_the_internal_node[T_SIDE][1][iP];

				if (visit_ray_traser[i] == false) {

					TOCHKA p2;
					center_cord3D_ray_tracing(i, t.nvtx, t.pa, p2, 100);

					if (distance_Point_to_Segment3(p2, p0, p1, lx, ly, lz, ld2) < size_cell[i])
					{
						// точка p2 принадлежит отрезку [p0,p1].
						if (!(((t.whot_is_block[i] == t.whot_is_block[i0]) ||
							(b[t.whot_is_block[i]].imatid == b[t.whot_is_block[i1]].imatid)))) {
							permissible758 = false; // не достижима.
						}
						else {
							// Продолжаем сканирование.
							itop_stack_ray_trayser++;
							stack_ray_trayser[itop_stack_ray_trayser] = i;
							visit_ray_traser[i] = true;

						}
					}

				}
			}
		}

		if (permissible758) {

			if ((t.neighbors_for_the_internal_node[T_SIDE][2] != nullptr) &&
				(t.neighbors_for_the_internal_node[T_SIDE][2][iP] > -1) &&
				(t.neighbors_for_the_internal_node[T_SIDE][2][iP] < t.maxelm)) {
				integer i = t.neighbors_for_the_internal_node[T_SIDE][2][iP];

				if (visit_ray_traser[i] == false) {

					TOCHKA p2;
					center_cord3D_ray_tracing(i, t.nvtx, t.pa, p2, 100);

					if (distance_Point_to_Segment3(p2, p0, p1, lx, ly, lz, ld2) < size_cell[i])
					{
						// точка p2 принадлежит отрезку [p0,p1].
						if (!(((t.whot_is_block[i] == t.whot_is_block[i0]) ||
							(b[t.whot_is_block[i]].imatid == b[t.whot_is_block[i1]].imatid)))) {
							permissible758 = false; // не достижима.
						}
						else {
							// Продолжаем сканирование.
							itop_stack_ray_trayser++;
							stack_ray_trayser[itop_stack_ray_trayser] = i;
							visit_ray_traser[i] = true;

						}
					}

				}
			}
		}

		if (permissible758) {

			if ((t.neighbors_for_the_internal_node[T_SIDE][3] != nullptr) &&
				(t.neighbors_for_the_internal_node[T_SIDE][3][iP] > -1) &&
				(t.neighbors_for_the_internal_node[T_SIDE][3][iP] < t.maxelm)) {
				integer i = t.neighbors_for_the_internal_node[T_SIDE][3][iP];

				if (visit_ray_traser[i] == false) {

					TOCHKA p2;
					center_cord3D_ray_tracing(i, t.nvtx, t.pa, p2, 100);

					if (distance_Point_to_Segment3(p2, p0, p1, lx, ly, lz, ld2) < size_cell[i])
					{
						// точка p2 принадлежит отрезку [p0,p1].
						if (!(((t.whot_is_block[i] == t.whot_is_block[i0]) ||
							(b[t.whot_is_block[i]].imatid == b[t.whot_is_block[i1]].imatid)))) {
							permissible758 = false; // не достижима.
						}
						else {
							// Продолжаем сканирование.
							itop_stack_ray_trayser++;
							stack_ray_trayser[itop_stack_ray_trayser] = i;
							visit_ray_traser[i] = true;

						}
					}

				}
			}
		}


		if (permissible758) {


			if ((t.neighbors_for_the_internal_node[B_SIDE][0] != nullptr) &&
				(t.neighbors_for_the_internal_node[B_SIDE][0][iP] > -1)&&
				(t.neighbors_for_the_internal_node[B_SIDE][0][iP] < t.maxelm)) 
			{
				integer i = t.neighbors_for_the_internal_node[B_SIDE][0][iP];

				if (visit_ray_traser[i] == false) {

					TOCHKA p2;
					center_cord3D_ray_tracing(i, t.nvtx, t.pa, p2, 100);

					if (distance_Point_to_Segment3(p2, p0, p1, lx, ly, lz, ld2) < size_cell[i])
					{
						// точка p2 принадлежит отрезку [p0,p1].
						if (!(((t.whot_is_block[i] == t.whot_is_block[i0]) ||
							(b[t.whot_is_block[i]].imatid == b[t.whot_is_block[i1]].imatid)))) {
							permissible758 = false; // не достижима.
						}
						else {
							// Продолжаем сканирование.
							itop_stack_ray_trayser++;
							stack_ray_trayser[itop_stack_ray_trayser] = i;
							visit_ray_traser[i] = true;
						}
					}
				}
			}
		}

		if (permissible758) {


			if ((t.neighbors_for_the_internal_node[B_SIDE][1] != nullptr) &&
				(t.neighbors_for_the_internal_node[B_SIDE][1][iP] > -1) &&
				(t.neighbors_for_the_internal_node[B_SIDE][1][iP] < t.maxelm))
			{
				integer i = t.neighbors_for_the_internal_node[B_SIDE][1][iP];

				if (visit_ray_traser[i] == false) {

					TOCHKA p2;
					center_cord3D_ray_tracing(i, t.nvtx, t.pa, p2, 100);

					if (distance_Point_to_Segment3(p2, p0, p1, lx, ly, lz, ld2) < size_cell[i])
					{
						// точка p2 принадлежит отрезку [p0,p1].
						if (!(((t.whot_is_block[i] == t.whot_is_block[i0]) ||
							(b[t.whot_is_block[i]].imatid == b[t.whot_is_block[i1]].imatid)))) {
							permissible758 = false; // не достижима.
						}
						else {
							// Продолжаем сканирование.
							itop_stack_ray_trayser++;
							stack_ray_trayser[itop_stack_ray_trayser] = i;
							visit_ray_traser[i] = true;
						}
					}
				}
			}
		}

		if (permissible758) {


			if ((t.neighbors_for_the_internal_node[B_SIDE][2] != nullptr) &&
				(t.neighbors_for_the_internal_node[B_SIDE][2][iP] > -1) &&
				(t.neighbors_for_the_internal_node[B_SIDE][2][iP] < t.maxelm))
			{
				integer i = t.neighbors_for_the_internal_node[B_SIDE][2][iP];

				if (visit_ray_traser[i] == false) {

					TOCHKA p2;
					center_cord3D_ray_tracing(i, t.nvtx, t.pa, p2, 100);

					if (distance_Point_to_Segment3(p2, p0, p1, lx, ly, lz, ld2) < size_cell[i])
					{
						// точка p2 принадлежит отрезку [p0,p1].
						if (!(((t.whot_is_block[i] == t.whot_is_block[i0]) ||
							(b[t.whot_is_block[i]].imatid == b[t.whot_is_block[i1]].imatid)))) {
							permissible758 = false; // не достижима.
						}
						else {
							// Продолжаем сканирование.
							itop_stack_ray_trayser++;
							stack_ray_trayser[itop_stack_ray_trayser] = i;
							visit_ray_traser[i] = true;
						}
					}
				}
			}
		}

		if (permissible758) {


			if ((t.neighbors_for_the_internal_node[B_SIDE][3] != nullptr) &&
				(t.neighbors_for_the_internal_node[B_SIDE][3][iP] > -1) &&
				(t.neighbors_for_the_internal_node[B_SIDE][3][iP] < t.maxelm))
			{
				integer i = t.neighbors_for_the_internal_node[B_SIDE][3][iP];

				if (visit_ray_traser[i] == false) {

					TOCHKA p2;
					center_cord3D_ray_tracing(i, t.nvtx, t.pa, p2, 100);

					if (distance_Point_to_Segment3(p2, p0, p1, lx, ly, lz, ld2) < size_cell[i])
					{
						// точка p2 принадлежит отрезку [p0,p1].
						if (!(((t.whot_is_block[i] == t.whot_is_block[i0]) ||
							(b[t.whot_is_block[i]].imatid == b[t.whot_is_block[i1]].imatid)))) {
							permissible758 = false; // не достижима.
						}
						else {
							// Продолжаем сканирование.
							itop_stack_ray_trayser++;
							stack_ray_trayser[itop_stack_ray_trayser] = i;
							visit_ray_traser[i] = true;
						}
					}
				}
			}
		}

	}


	if (visit_ray_traser[i1] == false) {
		// Точка i1 не была достигнута, путь преградил hollow block.
		permissible758 = false; // не достижимо.
	}

	delete[] visit_ray_traser;
} // ray_tracing_Alice


  // p0, i0 - diod, p1, i1 - point.
  // p3 - центр отрезка, radius2 - квадрат половины длины отрезка.
  // _fluid_only - более быстродействующий вариант за счёт избавления от ряда проверок.
  // Данный код корректно работает только если все ячейки расчётной области жидкие, есть hollow блоки,
  // но solid блоков вообще нет.
  // Работает только на структурированной сетке.
void ray_tracing_fluid_only(TOCHKA p3, TOCHKA p0, TOCHKA p1, integer i0, integer i1, TEMPER& t, BLOCK*& b, doublereal radius2) {


	doublereal lx = (p1.x - p0.x);
	doublereal ly = (p1.y - p0.y);
	doublereal lz = (p1.z - p0.z);
	doublereal ld2 = (lx * lx + ly * ly + lz * lz);


	visit_ray_traser = new bool[t.maxelm];
	for (integer i = 0; i < t.maxelm; i++) visit_ray_traser[i] = false;

	while (itop_stack_ray_trayser > -1) {
		integer iP = stack_ray_trayser[itop_stack_ray_trayser];
		itop_stack_ray_trayser--;


		// вызов.

		if ((t.neighbors_for_the_internal_node[E_SIDE][0] != nullptr)&&
			(t.neighbors_for_the_internal_node[E_SIDE][0][iP] < t.maxelm)) {
			integer i = t.neighbors_for_the_internal_node[E_SIDE][0][iP];

			if (visit_ray_traser[i] == false) {

				TOCHKA p2;
				center_cord3D_ray_tracing(i, t.nvtx, t.pa, p2, 100);


				if (distance_Point_to_Segment3(p2, p0, p1, lx, ly, lz, ld2) < size_cell[i])
				{

					// точка p2 принадлежит отрезку [p0,p1].

					// Продолжаем сканирование.
					itop_stack_ray_trayser++;
					stack_ray_trayser[itop_stack_ray_trayser] = i;
					visit_ray_traser[i] = true;

				}
			}
		}

		if ((t.neighbors_for_the_internal_node[W_SIDE][0] != nullptr) &&
			(t.neighbors_for_the_internal_node[W_SIDE][0][iP] < t.maxelm)) {
			integer i = t.neighbors_for_the_internal_node[W_SIDE][0][iP];

			if (visit_ray_traser[i] == false) {

				TOCHKA p2;
				center_cord3D_ray_tracing(i, t.nvtx, t.pa, p2, 100);

				if (distance_Point_to_Segment3(p2, p0, p1, lx, ly, lz, ld2) < size_cell[i])
				{

					// точка p2 принадлежит отрезку [p0,p1].

					// Продолжаем сканирование.
					itop_stack_ray_trayser++;
					stack_ray_trayser[itop_stack_ray_trayser] = i;
					visit_ray_traser[i] = true;

				}
			}
		}


		if ((t.neighbors_for_the_internal_node[N_SIDE][0] != nullptr)&&
			(t.neighbors_for_the_internal_node[N_SIDE][0][iP] < t.maxelm)) {
			integer i = t.neighbors_for_the_internal_node[N_SIDE][0][iP];

			if (visit_ray_traser[i] == false) {

				TOCHKA p2;
				center_cord3D_ray_tracing(i, t.nvtx, t.pa, p2, 100);


				if (distance_Point_to_Segment3(p2, p0, p1, lx, ly, lz, ld2) < size_cell[i])
				{

					// точка p2 принадлежит отрезку [p0,p1].

					// Продолжаем сканирование.
					itop_stack_ray_trayser++;
					stack_ray_trayser[itop_stack_ray_trayser] = i;
					visit_ray_traser[i] = true;

				}
			}
		}



		if ((t.neighbors_for_the_internal_node[S_SIDE][0] != nullptr)&&
			(t.neighbors_for_the_internal_node[S_SIDE][0][iP] < t.maxelm)) {
			integer i = t.neighbors_for_the_internal_node[S_SIDE][0][iP];

			if (visit_ray_traser[i] == false) {

				TOCHKA p2;
				center_cord3D_ray_tracing(i, t.nvtx, t.pa, p2, 100);

				if (distance_Point_to_Segment3(p2, p0, p1, lx, ly, lz, ld2) < size_cell[i])
				{

					// точка p2 принадлежит отрезку [p0,p1].

					// Продолжаем сканирование.
					itop_stack_ray_trayser++;
					stack_ray_trayser[itop_stack_ray_trayser] = i;
					visit_ray_traser[i] = true;

				}
			}
		}


		if ((t.neighbors_for_the_internal_node[T_SIDE][0] != nullptr)&&
			(t.neighbors_for_the_internal_node[T_SIDE][0][iP] < t.maxelm)) {
			integer i = t.neighbors_for_the_internal_node[T_SIDE][0][iP];

			if (visit_ray_traser[i] == false) {

				TOCHKA p2;
				center_cord3D_ray_tracing(i, t.nvtx, t.pa, p2, 100);

				if (distance_Point_to_Segment3(p2, p0, p1, lx, ly, lz, ld2) < size_cell[i])
				{
					// точка p2 принадлежит отрезку [p0,p1].

					// Продолжаем сканирование.
					itop_stack_ray_trayser++;
					stack_ray_trayser[itop_stack_ray_trayser] = i;
					visit_ray_traser[i] = true;


				}
			}
		}



		if ((t.neighbors_for_the_internal_node[B_SIDE][0] != nullptr)&&
			(t.neighbors_for_the_internal_node[B_SIDE][0][iP] < t.maxelm)) {
			integer i = t.neighbors_for_the_internal_node[B_SIDE][0][iP];

			if (visit_ray_traser[i] == false) {

				TOCHKA p2;
				center_cord3D_ray_tracing(i, t.nvtx, t.pa, p2, 100);

				if (distance_Point_to_Segment3(p2, p0, p1, lx, ly, lz, ld2) < size_cell[i])
				{
					// точка p2 принадлежит отрезку [p0,p1].

					// Продолжаем сканирование.
					itop_stack_ray_trayser++;
					stack_ray_trayser[itop_stack_ray_trayser] = i;
					visit_ray_traser[i] = true;

				}
			}
		}
	}



	if (visit_ray_traser[i1] == false) {
		// Точка i1 не была достигнута, путь преградил hollow block.
		permissible758 = false; // не достижимо.
	}

	delete[] visit_ray_traser;
} // ray_tracing_fluid_only

  // p0, i0 - diod, p1, i1 - point.
  // p3 - центр отрезка, radius2 - квадрат половины длины отрезка.
  // _fluid_only - более быстродействующий вариант за счёт избавления от ряда проверок.
  // Данный код корректно работает только если все ячейки расчётной области жидкие, есть hollow блоки,
  // но solid блоков вообще нет.
  // Работает только на структурированной сетке.
void ray_tracing_fluid_only2(TOCHKA p3, TOCHKA p0, TOCHKA p1, integer i0, integer i1, TEMPER& t, BLOCK*& b, doublereal radius2) {


	doublereal lx = (p1.x - p0.x);
	doublereal ly = (p1.y - p0.y);
	doublereal lz = (p1.z - p0.z);
	doublereal ld2 = (lx * lx + ly * ly + lz * lz);


	visit_ray_traser = new bool[t.maxelm];
	for (integer i = 0; i < t.maxelm; i++) visit_ray_traser[i] = false;

	doublereal* distance_arr = nullptr;
	distance_arr = new doublereal[6];
	integer* marker_arr = nullptr;
	marker_arr = new integer[6];

	while (itop_stack_ray_trayser > -1) {
		integer iP = stack_ray_trayser[itop_stack_ray_trayser];
		itop_stack_ray_trayser--;


		// вызов.
		for (integer j = 0; j < 6; j++) {
			distance_arr[j] = 1.0e30;
			marker_arr[j] = -1;
		}

		if ((t.neighbors_for_the_internal_node[E_SIDE][0] != nullptr)&&
			(t.neighbors_for_the_internal_node[E_SIDE][0][iP] < t.maxelm)) {
			integer i = t.neighbors_for_the_internal_node[E_SIDE][0][iP];

			if (visit_ray_traser[i] == false) {

				TOCHKA p2;
				center_cord3D_ray_tracing(i, t.nvtx, t.pa, p2, 100);

				distance_arr[0] = distance_Point_to_Segment3(p2, p0, p1, lx, ly, lz, ld2);
				marker_arr[0] = i;
			}
		}

		if ((t.neighbors_for_the_internal_node[W_SIDE][0] != nullptr) &&
			(t.neighbors_for_the_internal_node[W_SIDE][0][iP] < t.maxelm)) {
			integer i = t.neighbors_for_the_internal_node[W_SIDE][0][iP];

			if (visit_ray_traser[i] == false) {

				TOCHKA p2;
				center_cord3D_ray_tracing(i, t.nvtx, t.pa, p2, 100);

				distance_arr[1] = distance_Point_to_Segment3(p2, p0, p1, lx, ly, lz, ld2);
				marker_arr[1] = i;

			}
		}


		if ((t.neighbors_for_the_internal_node[N_SIDE][0] != nullptr) &&
			(t.neighbors_for_the_internal_node[N_SIDE][0][iP] < t.maxelm)) {
			integer i = t.neighbors_for_the_internal_node[N_SIDE][0][iP];

			if (visit_ray_traser[i] == false) {

				TOCHKA p2;
				center_cord3D_ray_tracing(i, t.nvtx, t.pa, p2, 100);

				distance_arr[2] = distance_Point_to_Segment3(p2, p0, p1, lx, ly, lz, ld2);
				marker_arr[2] = i;

			}
		}



		if ((t.neighbors_for_the_internal_node[S_SIDE][0] != nullptr)&&
			(t.neighbors_for_the_internal_node[S_SIDE][0][iP] < t.maxelm)) {
			integer i = t.neighbors_for_the_internal_node[S_SIDE][0][iP];

			if (visit_ray_traser[i] == false) {

				TOCHKA p2;
				center_cord3D_ray_tracing(i, t.nvtx, t.pa, p2, 100);

				distance_arr[3] = distance_Point_to_Segment3(p2, p0, p1, lx, ly, lz, ld2);
				marker_arr[3] = i;
			}
		}


		if ((t.neighbors_for_the_internal_node[T_SIDE][0] != nullptr)&&
			(t.neighbors_for_the_internal_node[T_SIDE][0][iP] < t.maxelm)) {
			integer i = t.neighbors_for_the_internal_node[T_SIDE][0][iP];

			if (visit_ray_traser[i] == false) {

				TOCHKA p2;
				center_cord3D_ray_tracing(i, t.nvtx, t.pa, p2, 100);

				distance_arr[4] = distance_Point_to_Segment3(p2, p0, p1, lx, ly, lz, ld2);
				marker_arr[4] = i;
			}
		}



		if ((t.neighbors_for_the_internal_node[B_SIDE][0] != nullptr)&&
			(t.neighbors_for_the_internal_node[B_SIDE][0][iP] < t.maxelm)) {
			integer i = t.neighbors_for_the_internal_node[B_SIDE][0][iP];

			if (visit_ray_traser[i] == false) {

				TOCHKA p2;
				center_cord3D_ray_tracing(i, t.nvtx, t.pa, p2, 100);

				distance_arr[5] = distance_Point_to_Segment3(p2, p0, p1, lx, ly, lz, ld2);
				marker_arr[5] = i;
			}
		}

		// Ищем точку с минимальным расстоянием, среди шести соседей.

		doublereal dmin = 1.0e29;
		integer jmin = -1;
		for (integer j = 0; j < 6; j++) {
			if (distance_arr[j] < dmin) {
				dmin = distance_arr[j];
				jmin = j;
			}
		}

		// Запускаемся снова из точки с минимальным расстонием,
		// в остальные точки больше не заходим.

		for (integer j = 0; j < 6; j++) {

			integer i = marker_arr[j];

			if (j == jmin) {				

				// точка p2 принадлежит отрезку [p0,p1].

				// Продолжаем сканирование.
				itop_stack_ray_trayser++;
				stack_ray_trayser[itop_stack_ray_trayser] = i;
				visit_ray_traser[i] = true;
			}
			else {
				visit_ray_traser[i] = true;
			}

		}
	}


	delete[] distance_arr;
	delete[] marker_arr;

	if (visit_ray_traser[i1] == false) {
		// Точка i1 не была достигнута, путь преградил hollow block.
		permissible758 = false; // не достижимо.
	}

	delete[] visit_ray_traser;
} // ray_tracing_fluid_only2



  // p0, i0 - diod, p1, i1 - point.
  // p3 - центр отрезка, radius2 - квадрат половины длины отрезка.
  // _fluid_only - более быстродействующий вариант за счёт избавления от ряда проверок.
  // Данный код корректно работает только если все ячейки расчётной области жидкие, есть hollow блоки,
  // но solid блоков вообще нет.
  // Работает на АЛИС сетке.
void ray_tracing_fluid_only_Alice(TOCHKA p3, TOCHKA p0, TOCHKA p1, integer i0, integer i1,
	TEMPER& t, BLOCK*& b, doublereal radius2)
{


	doublereal lx = (p1.x - p0.x);
	doublereal ly = (p1.y - p0.y);
	doublereal lz = (p1.z - p0.z);
	doublereal ld2 = (lx * lx + ly * ly + lz * lz);


	visit_ray_traser = new bool[t.maxelm];
	for (integer i = 0; i < t.maxelm; i++) visit_ray_traser[i] = false;

	while (itop_stack_ray_trayser > -1) {
		integer iP = stack_ray_trayser[itop_stack_ray_trayser];
		itop_stack_ray_trayser--;


		// вызов.

		if ((t.neighbors_for_the_internal_node[E_SIDE][0] != nullptr)&&
			(t.neighbors_for_the_internal_node[E_SIDE][0][iP] > -1) &&
			(t.neighbors_for_the_internal_node[E_SIDE][0][iP] < t.maxelm))
		{
			integer i = t.neighbors_for_the_internal_node[E_SIDE][0][iP];

			if (visit_ray_traser[i] == false) {

				TOCHKA p2;
				center_cord3D_ray_tracing(i, t.nvtx, t.pa, p2, 100);


				if (distance_Point_to_Segment3(p2, p0, p1, lx, ly, lz, ld2) < size_cell[i])
				{

					// точка p2 принадлежит отрезку [p0,p1].

					// Продолжаем сканирование.
					itop_stack_ray_trayser++;
					stack_ray_trayser[itop_stack_ray_trayser] = i;
					visit_ray_traser[i] = true;

				}
			}
		}

		if ((t.neighbors_for_the_internal_node[E_SIDE][1] != nullptr) &&
			(t.neighbors_for_the_internal_node[E_SIDE][1][iP] > -1) &&
			(t.neighbors_for_the_internal_node[E_SIDE][1][iP] < t.maxelm))
		{
			integer i = t.neighbors_for_the_internal_node[E_SIDE][1][iP];

			if (visit_ray_traser[i] == false) {

				TOCHKA p2;
				center_cord3D_ray_tracing(i, t.nvtx, t.pa, p2, 100);


				if (distance_Point_to_Segment3(p2, p0, p1, lx, ly, lz, ld2) < size_cell[i])
				{

					// точка p2 принадлежит отрезку [p0,p1].

					// Продолжаем сканирование.
					itop_stack_ray_trayser++;
					stack_ray_trayser[itop_stack_ray_trayser] = i;
					visit_ray_traser[i] = true;

				}
			}
		}

		if ((t.neighbors_for_the_internal_node[E_SIDE][2] != nullptr) && 
			(t.neighbors_for_the_internal_node[E_SIDE][2][iP] > -1) &&
			(t.neighbors_for_the_internal_node[E_SIDE][2][iP] < t.maxelm))
		{
			integer i = t.neighbors_for_the_internal_node[E_SIDE][2][iP];

			if (visit_ray_traser[i] == false) {

				TOCHKA p2;
				center_cord3D_ray_tracing(i, t.nvtx, t.pa, p2, 100);


				if (distance_Point_to_Segment3(p2, p0, p1, lx, ly, lz, ld2) < size_cell[i])
				{

					// точка p2 принадлежит отрезку [p0,p1].

					// Продолжаем сканирование.
					itop_stack_ray_trayser++;
					stack_ray_trayser[itop_stack_ray_trayser] = i;
					visit_ray_traser[i] = true;

				}
			}
		}

		if ((t.neighbors_for_the_internal_node[E_SIDE][3] != nullptr) && 
			(t.neighbors_for_the_internal_node[E_SIDE][3][iP] > -1) &&
			(t.neighbors_for_the_internal_node[E_SIDE][3][iP] < t.maxelm))
		{
			integer i = t.neighbors_for_the_internal_node[E_SIDE][3][iP];

			if (visit_ray_traser[i] == false) {

				TOCHKA p2;
				center_cord3D_ray_tracing(i, t.nvtx, t.pa, p2, 100);


				if (distance_Point_to_Segment3(p2, p0, p1, lx, ly, lz, ld2) < size_cell[i])
				{

					// точка p2 принадлежит отрезку [p0,p1].

					// Продолжаем сканирование.
					itop_stack_ray_trayser++;
					stack_ray_trayser[itop_stack_ray_trayser] = i;
					visit_ray_traser[i] = true;

				}
			}
		}


		if ((t.neighbors_for_the_internal_node[W_SIDE][0] != nullptr) &&
			(t.neighbors_for_the_internal_node[W_SIDE][0][iP] > -1) &&
			(t.neighbors_for_the_internal_node[W_SIDE][0][iP] < t.maxelm))
		{
			integer i = t.neighbors_for_the_internal_node[W_SIDE][0][iP];

			if (visit_ray_traser[i] == false) {

				TOCHKA p2;
				center_cord3D_ray_tracing(i, t.nvtx, t.pa, p2, 100);

				if (distance_Point_to_Segment3(p2, p0, p1, lx, ly, lz, ld2) < size_cell[i])
				{

					// точка p2 принадлежит отрезку [p0,p1].

					// Продолжаем сканирование.
					itop_stack_ray_trayser++;
					stack_ray_trayser[itop_stack_ray_trayser] = i;
					visit_ray_traser[i] = true;

				}
			}
		}

		if ((t.neighbors_for_the_internal_node[W_SIDE][1] != nullptr) && 
			(t.neighbors_for_the_internal_node[W_SIDE][1][iP] > -1) &&
			(t.neighbors_for_the_internal_node[W_SIDE][1][iP] < t.maxelm))
		{
			integer i = t.neighbors_for_the_internal_node[W_SIDE][1][iP];

			if (visit_ray_traser[i] == false) {

				TOCHKA p2;
				center_cord3D_ray_tracing(i, t.nvtx, t.pa, p2, 100);

				if (distance_Point_to_Segment3(p2, p0, p1, lx, ly, lz, ld2) < size_cell[i])
				{

					// точка p2 принадлежит отрезку [p0,p1].

					// Продолжаем сканирование.
					itop_stack_ray_trayser++;
					stack_ray_trayser[itop_stack_ray_trayser] = i;
					visit_ray_traser[i] = true;

				}
			}
		}

		if ((t.neighbors_for_the_internal_node[W_SIDE][2] != nullptr) && 
			(t.neighbors_for_the_internal_node[W_SIDE][2][iP] > -1) &&
			(t.neighbors_for_the_internal_node[W_SIDE][2][iP] < t.maxelm))
		{
			integer i = t.neighbors_for_the_internal_node[W_SIDE][2][iP];

			if (visit_ray_traser[i] == false) {

				TOCHKA p2;
				center_cord3D_ray_tracing(i, t.nvtx, t.pa, p2, 100);

				if (distance_Point_to_Segment3(p2, p0, p1, lx, ly, lz, ld2) < size_cell[i])
				{

					// точка p2 принадлежит отрезку [p0,p1].

					// Продолжаем сканирование.
					itop_stack_ray_trayser++;
					stack_ray_trayser[itop_stack_ray_trayser] = i;
					visit_ray_traser[i] = true;

				}
			}
		}

		if ((t.neighbors_for_the_internal_node[W_SIDE][3] != nullptr) && 
			(t.neighbors_for_the_internal_node[W_SIDE][3][iP] > -1) &&
			(t.neighbors_for_the_internal_node[W_SIDE][3][iP] < t.maxelm))
		{
			integer i = t.neighbors_for_the_internal_node[W_SIDE][3][iP];

			if (visit_ray_traser[i] == false) {

				TOCHKA p2;
				center_cord3D_ray_tracing(i, t.nvtx, t.pa, p2, 100);

				if (distance_Point_to_Segment3(p2, p0, p1, lx, ly, lz, ld2) < size_cell[i])
				{

					// точка p2 принадлежит отрезку [p0,p1].

					// Продолжаем сканирование.
					itop_stack_ray_trayser++;
					stack_ray_trayser[itop_stack_ray_trayser] = i;
					visit_ray_traser[i] = true;

				}
			}
		}


		if ((t.neighbors_for_the_internal_node[N_SIDE][0] != nullptr) && 
			(t.neighbors_for_the_internal_node[N_SIDE][0][iP] > -1) &&
			(t.neighbors_for_the_internal_node[N_SIDE][0][iP] < t.maxelm)) {
			integer i = t.neighbors_for_the_internal_node[N_SIDE][0][iP];

			if (visit_ray_traser[i] == false) {

				TOCHKA p2;
				center_cord3D_ray_tracing(i, t.nvtx, t.pa, p2, 100);


				if (distance_Point_to_Segment3(p2, p0, p1, lx, ly, lz, ld2) < size_cell[i])
				{

					// точка p2 принадлежит отрезку [p0,p1].

					// Продолжаем сканирование.
					itop_stack_ray_trayser++;
					stack_ray_trayser[itop_stack_ray_trayser] = i;
					visit_ray_traser[i] = true;

				}
			}
		}

		if ((t.neighbors_for_the_internal_node[N_SIDE][1] != nullptr) && 
			(t.neighbors_for_the_internal_node[N_SIDE][1][iP] > -1) &&
			(t.neighbors_for_the_internal_node[N_SIDE][1][iP] < t.maxelm)) {
			integer i = t.neighbors_for_the_internal_node[N_SIDE][1][iP];

			if (visit_ray_traser[i] == false) {

				TOCHKA p2;
				center_cord3D_ray_tracing(i, t.nvtx, t.pa, p2, 100);


				if (distance_Point_to_Segment3(p2, p0, p1, lx, ly, lz, ld2) < size_cell[i])
				{

					// точка p2 принадлежит отрезку [p0,p1].

					// Продолжаем сканирование.
					itop_stack_ray_trayser++;
					stack_ray_trayser[itop_stack_ray_trayser] = i;
					visit_ray_traser[i] = true;

				}
			}
		}

		if ((t.neighbors_for_the_internal_node[N_SIDE][2] != nullptr) && 
			(t.neighbors_for_the_internal_node[N_SIDE][2][iP] > -1) &&
			(t.neighbors_for_the_internal_node[N_SIDE][2][iP] < t.maxelm)) {
			integer i = t.neighbors_for_the_internal_node[N_SIDE][2][iP];

			if (visit_ray_traser[i] == false) {

				TOCHKA p2;
				center_cord3D_ray_tracing(i, t.nvtx, t.pa, p2, 100);


				if (distance_Point_to_Segment3(p2, p0, p1, lx, ly, lz, ld2) < size_cell[i])
				{

					// точка p2 принадлежит отрезку [p0,p1].

					// Продолжаем сканирование.
					itop_stack_ray_trayser++;
					stack_ray_trayser[itop_stack_ray_trayser] = i;
					visit_ray_traser[i] = true;

				}
			}
		}

		if ((t.neighbors_for_the_internal_node[N_SIDE][3] != nullptr) &&
			(t.neighbors_for_the_internal_node[N_SIDE][3][iP] > -1) &&
			(t.neighbors_for_the_internal_node[N_SIDE][3][iP] < t.maxelm)) {
			integer i = t.neighbors_for_the_internal_node[N_SIDE][3][iP];

			if (visit_ray_traser[i] == false) {

				TOCHKA p2;
				center_cord3D_ray_tracing(i, t.nvtx, t.pa, p2, 100);


				if (distance_Point_to_Segment3(p2, p0, p1, lx, ly, lz, ld2) < size_cell[i])
				{

					// точка p2 принадлежит отрезку [p0,p1].

					// Продолжаем сканирование.
					itop_stack_ray_trayser++;
					stack_ray_trayser[itop_stack_ray_trayser] = i;
					visit_ray_traser[i] = true;

				}
			}
		}


		if ((t.neighbors_for_the_internal_node[S_SIDE][0] != nullptr )&&
			(t.neighbors_for_the_internal_node[S_SIDE][0][iP] > -1 )&&
			(t.neighbors_for_the_internal_node[S_SIDE][0][iP] < t.maxelm)) {
			integer i = t.neighbors_for_the_internal_node[S_SIDE][0][iP];

			if (visit_ray_traser[i] == false) {

				TOCHKA p2;
				center_cord3D_ray_tracing(i, t.nvtx, t.pa, p2, 100);

				if (distance_Point_to_Segment3(p2, p0, p1, lx, ly, lz, ld2) < size_cell[i])
				{

					// точка p2 принадлежит отрезку [p0,p1].

					// Продолжаем сканирование.
					itop_stack_ray_trayser++;
					stack_ray_trayser[itop_stack_ray_trayser] = i;
					visit_ray_traser[i] = true;

				}
			}
		}

		if ((t.neighbors_for_the_internal_node[S_SIDE][1] != nullptr) && 
			(t.neighbors_for_the_internal_node[S_SIDE][1][iP] > -1) &&
			(t.neighbors_for_the_internal_node[S_SIDE][1][iP] < t.maxelm)) {
			integer i = t.neighbors_for_the_internal_node[S_SIDE][1][iP];

			if (visit_ray_traser[i] == false) {

				TOCHKA p2;
				center_cord3D_ray_tracing(i, t.nvtx, t.pa, p2, 100);

				if (distance_Point_to_Segment3(p2, p0, p1, lx, ly, lz, ld2) < size_cell[i])
				{

					// точка p2 принадлежит отрезку [p0,p1].

					// Продолжаем сканирование.
					itop_stack_ray_trayser++;
					stack_ray_trayser[itop_stack_ray_trayser] = i;
					visit_ray_traser[i] = true;

				}
			}
		}

		if ((t.neighbors_for_the_internal_node[S_SIDE][2] != nullptr) &&
			(t.neighbors_for_the_internal_node[S_SIDE][2][iP] > -1) &&
			(t.neighbors_for_the_internal_node[S_SIDE][2][iP] < t.maxelm)) {
			integer i = t.neighbors_for_the_internal_node[S_SIDE][2][iP];

			if (visit_ray_traser[i] == false) {

				TOCHKA p2;
				center_cord3D_ray_tracing(i, t.nvtx, t.pa, p2, 100);

				if (distance_Point_to_Segment3(p2, p0, p1, lx, ly, lz, ld2) < size_cell[i])
				{

					// точка p2 принадлежит отрезку [p0,p1].

					// Продолжаем сканирование.
					itop_stack_ray_trayser++;
					stack_ray_trayser[itop_stack_ray_trayser] = i;
					visit_ray_traser[i] = true;

				}
			}
		}


		if ((t.neighbors_for_the_internal_node[S_SIDE][3] != nullptr) && 
			(t.neighbors_for_the_internal_node[S_SIDE][3][iP] > -1) &&
			(t.neighbors_for_the_internal_node[S_SIDE][3][iP] < t.maxelm)) {
			integer i = t.neighbors_for_the_internal_node[S_SIDE][3][iP];

			if (visit_ray_traser[i] == false) {

				TOCHKA p2;
				center_cord3D_ray_tracing(i, t.nvtx, t.pa, p2, 100);

				if (distance_Point_to_Segment3(p2, p0, p1, lx, ly, lz, ld2) < size_cell[i])
				{

					// точка p2 принадлежит отрезку [p0,p1].

					// Продолжаем сканирование.
					itop_stack_ray_trayser++;
					stack_ray_trayser[itop_stack_ray_trayser] = i;
					visit_ray_traser[i] = true;

				}
			}
		}

		if ((t.neighbors_for_the_internal_node[T_SIDE][0] != nullptr) &&
			(t.neighbors_for_the_internal_node[T_SIDE][0][iP] > -1)&&
			(t.neighbors_for_the_internal_node[T_SIDE][0][iP] < t.maxelm)) {
			integer i = t.neighbors_for_the_internal_node[T_SIDE][0][iP];

			if (visit_ray_traser[i] == false) {

				TOCHKA p2;
				center_cord3D_ray_tracing(i, t.nvtx, t.pa, p2, 100);

				if (distance_Point_to_Segment3(p2, p0, p1, lx, ly, lz, ld2) < size_cell[i])
				{
					// точка p2 принадлежит отрезку [p0,p1].

					// Продолжаем сканирование.
					itop_stack_ray_trayser++;
					stack_ray_trayser[itop_stack_ray_trayser] = i;
					visit_ray_traser[i] = true;


				}
			}
		}

		if ((t.neighbors_for_the_internal_node[T_SIDE][1] != nullptr) &&
			(t.neighbors_for_the_internal_node[T_SIDE][1][iP] > -1) &&
			(t.neighbors_for_the_internal_node[T_SIDE][1][iP] < t.maxelm)) {
			integer i = t.neighbors_for_the_internal_node[T_SIDE][1][iP];

			if (visit_ray_traser[i] == false) {

				TOCHKA p2;
				center_cord3D_ray_tracing(i, t.nvtx, t.pa, p2, 100);

				if (distance_Point_to_Segment3(p2, p0, p1, lx, ly, lz, ld2) < size_cell[i])
				{
					// точка p2 принадлежит отрезку [p0,p1].

					// Продолжаем сканирование.
					itop_stack_ray_trayser++;
					stack_ray_trayser[itop_stack_ray_trayser] = i;
					visit_ray_traser[i] = true;


				}
			}
		}

		if ((t.neighbors_for_the_internal_node[T_SIDE][2] != nullptr) && 
			(t.neighbors_for_the_internal_node[T_SIDE][2][iP] > -1) &&
			(t.neighbors_for_the_internal_node[T_SIDE][2][iP] < t.maxelm)) {
			integer i = t.neighbors_for_the_internal_node[T_SIDE][2][iP];

			if (visit_ray_traser[i] == false) {

				TOCHKA p2;
				center_cord3D_ray_tracing(i, t.nvtx, t.pa, p2, 100);

				if (distance_Point_to_Segment3(p2, p0, p1, lx, ly, lz, ld2) < size_cell[i])
				{
					// точка p2 принадлежит отрезку [p0,p1].

					// Продолжаем сканирование.
					itop_stack_ray_trayser++;
					stack_ray_trayser[itop_stack_ray_trayser] = i;
					visit_ray_traser[i] = true;

				}
			}
		}

		if ((t.neighbors_for_the_internal_node[T_SIDE][3] != nullptr) && 
			(t.neighbors_for_the_internal_node[T_SIDE][3][iP] > -1) &&
			(t.neighbors_for_the_internal_node[T_SIDE][3][iP] < t.maxelm)) {
			integer i = t.neighbors_for_the_internal_node[T_SIDE][3][iP];

			if (visit_ray_traser[i] == false) {

				TOCHKA p2;
				center_cord3D_ray_tracing(i, t.nvtx, t.pa, p2, 100);

				if (distance_Point_to_Segment3(p2, p0, p1, lx, ly, lz, ld2) < size_cell[i])
				{
					// точка p2 принадлежит отрезку [p0,p1].

					// Продолжаем сканирование.
					itop_stack_ray_trayser++;
					stack_ray_trayser[itop_stack_ray_trayser] = i;
					visit_ray_traser[i] = true;

				}
			}
		}

		if ((t.neighbors_for_the_internal_node[B_SIDE][0] != nullptr) && 
			(t.neighbors_for_the_internal_node[B_SIDE][0][iP] > -1)&&
			(t.neighbors_for_the_internal_node[B_SIDE][0][iP] < t.maxelm)) {
			integer i = t.neighbors_for_the_internal_node[B_SIDE][0][iP];

			if (visit_ray_traser[i] == false) {

				TOCHKA p2;
				center_cord3D_ray_tracing(i, t.nvtx, t.pa, p2, 100);

				if (distance_Point_to_Segment3(p2, p0, p1, lx, ly, lz, ld2) < size_cell[i])
				{
					// точка p2 принадлежит отрезку [p0,p1].

					// Продолжаем сканирование.
					itop_stack_ray_trayser++;
					stack_ray_trayser[itop_stack_ray_trayser] = i;
					visit_ray_traser[i] = true;

				}
			}
		}


		if ((t.neighbors_for_the_internal_node[B_SIDE][1] != nullptr) && 
			(t.neighbors_for_the_internal_node[B_SIDE][1][iP] > -1) &&
			(t.neighbors_for_the_internal_node[B_SIDE][1][iP] < t.maxelm)) {
			integer i = t.neighbors_for_the_internal_node[B_SIDE][1][iP];

			if (visit_ray_traser[i] == false) {

				TOCHKA p2;
				center_cord3D_ray_tracing(i, t.nvtx, t.pa, p2, 100);

				if (distance_Point_to_Segment3(p2, p0, p1, lx, ly, lz, ld2) < size_cell[i])
				{
					// точка p2 принадлежит отрезку [p0,p1].

					// Продолжаем сканирование.
					itop_stack_ray_trayser++;
					stack_ray_trayser[itop_stack_ray_trayser] = i;
					visit_ray_traser[i] = true;

				}
			}
		}


		if ((t.neighbors_for_the_internal_node[B_SIDE][2] != nullptr) && 
			(t.neighbors_for_the_internal_node[B_SIDE][2][iP] > -1) &&
			(t.neighbors_for_the_internal_node[B_SIDE][2][iP] < t.maxelm)) {
			integer i = t.neighbors_for_the_internal_node[B_SIDE][2][iP];

			if (visit_ray_traser[i] == false) {

				TOCHKA p2;
				center_cord3D_ray_tracing(i, t.nvtx, t.pa, p2, 100);

				if (distance_Point_to_Segment3(p2, p0, p1, lx, ly, lz, ld2) < size_cell[i])
				{
					// точка p2 принадлежит отрезку [p0,p1].

					// Продолжаем сканирование.
					itop_stack_ray_trayser++;
					stack_ray_trayser[itop_stack_ray_trayser] = i;
					visit_ray_traser[i] = true;

				}
			}
		}

		if ((t.neighbors_for_the_internal_node[B_SIDE][3] != nullptr) && 
			(t.neighbors_for_the_internal_node[B_SIDE][3][iP] > -1) &&
			(t.neighbors_for_the_internal_node[B_SIDE][3][iP] < t.maxelm)) {
			integer i = t.neighbors_for_the_internal_node[B_SIDE][3][iP];

			if (visit_ray_traser[i] == false) {

				TOCHKA p2;
				center_cord3D_ray_tracing(i, t.nvtx, t.pa, p2, 100);

				if (distance_Point_to_Segment3(p2, p0, p1, lx, ly, lz, ld2) < size_cell[i])
				{
					// точка p2 принадлежит отрезку [p0,p1].

					// Продолжаем сканирование.
					itop_stack_ray_trayser++;
					stack_ray_trayser[itop_stack_ray_trayser] = i;
					visit_ray_traser[i] = true;

				}
			}
		}

	} // до тех пор пока стек не пуст.



	if (visit_ray_traser[i1] == false) {
		// Точка i1 не была достигнута, путь преградил hollow block.
		permissible758 = false; // не достижимо.
	}

	delete[] visit_ray_traser;
} // ray_tracing_fluid_only_Alice

  // Достижима ли точка p1 из центра светодиода p0.
bool is_permisseble(TOCHKA p0, TOCHKA p1, integer i0, integer i1, TEMPER& t, BLOCK*& b)
{
	permissible758 = true; // По умолчанию точка достижима.
	if (b[t.whot_is_block[i1]].itype == PHYSICS_TYPE_IN_BODY::FLUID) {

		TOCHKA p3;
		p3.x = 0.5 * (p0.x + p1.x);
		p3.y = 0.5 * (p0.y + p1.y);
		p3.z = 0.5 * (p0.z + p1.z);
		doublereal radius2 = (p3.x - p0.x) * (p3.x - p0.x) + (p3.y - p0.y) * (p3.y - p0.y) + (p3.z - p0.z) * (p3.z - p0.z);

		// Ищем точки принадлежащие отрезку [p0,p1].

		itop_stack_ray_trayser = 0;
		stack_ray_trayser[itop_stack_ray_trayser] = i0;

		if (bray_tracing_fluid_only) {
			if (b_on_adaptive_local_refinement_mesh) {
				ray_tracing_fluid_only_Alice(p3, p0, p1, i0, i1, t, b, radius2);
			}
			else {
				ray_tracing_fluid_only(p3, p0, p1, i0, i1, t, b, radius2);
				// ray_tracing_fluid_only2() ни в коем случае использовать нельзя 18,05,2020.
				//ray_tracing_fluid_only2(p3, p0, p1, i0, i1, t, b, radius2); // не работает.
			}
		}
		else {
			if (b_on_adaptive_local_refinement_mesh) {
				ray_tracing_Alice(p3, p0, p1, i0, i1, t, b, radius2);
			}
			else {
				ray_tracing(p3, p0, p1, i0, i1, t, b, radius2);
			}
		}

	}
	else {
		permissible758 = false;
	}

	return permissible758;
} // is_permisseble

void calculate_light_flux(doublereal*& myF, TEMPER& t,
	BLOCK*& b, integer& lb) {

	size_cell = new doublereal[t.maxelm];

	stack_ray_trayser = new integer[t.maxelm + 2];
	doublereal* distance = new doublereal[t.maxelm + t.maxbound];
	bray_tracing_fluid_only = true;
	for (integer i = 0; i < t.maxelm; i++) {
		myF[i] = 0.0;
		size_cell[i] = kR_size_cell * volume3D_ray_tracing(i, t.nvtx, t.pa);
		if (b[t.whot_is_block[i]].itype != PHYSICS_TYPE_IN_BODY::FLUID) bray_tracing_fluid_only = false;
	}


#ifdef _OPENMP 
	omp_set_num_threads(inumcore); // установка числа потоков
#endif

	for (integer ibid = 0; ibid < lb; ibid++) {
		for (integer i = 0; i < t.maxelm; i++) {
			distance[i] = 1.0e30;
			if (b[t.whot_is_block[i]].itype == PHYSICS_TYPE_IN_BODY::SOLID) {
				distance[i] = -0.1e-3;// одна десятая мм
			}
		}

		if (b[ibid].arr_Sc[0] > 1.0e-30) {
			TOCHKA p0;
			integer i0=-1;// Если не будет инициализировано далее, то со значением -1  намеренно вызовет ошибку памяти.

			printf("power %e\n", b[ibid].arr_Sc[0]);
			doublereal x0 = 0.5 * (b[ibid].g.xS + b[ibid].g.xE);
			doublereal y0 = 0.5 * (b[ibid].g.yS + b[ibid].g.yE);
			doublereal z0 = 0.5 * (b[ibid].g.zS + b[ibid].g.zE);

			doublereal dist32 = 1.0e30;
			for (integer i = 0; i < t.maxelm; i++) {
				TOCHKA p2;
				center_cord3D(i, t.nvtx, t.pa, p2, 100);
				doublereal distance = sqrt((p2.x - x0) * (p2.x - x0) + (p2.y - y0) * (p2.y - y0) + (p2.z - z0) * (p2.z - z0));

				if (distance < dist32) {
					dist32 = distance;
					p0 = p2;
					i0 = i;
				}
			}
			// В p0 хранится центр текущего диода.

#pragma omp parallel for
			for (integer i = 0; i < t.maxelm; i++) {
				TOCHKA p1;
				center_cord3D(i, t.nvtx, t.pa, p1, 100);
				doublereal dist0 = sqrt((p1.x - x0) * (p1.x - x0) + (p1.y - y0) * (p1.y - y0) + (p1.z - z0) * (p1.z - z0));
				if (dist0 < distance[i]) {
					bool bpermisseble = is_permisseble(p0, p1, i0, i, t, b);
					if (bpermisseble) {
						distance[i] = dist0;
					}
				}
			}
			for (integer i = 0; i < t.maxelm; i++) {
				if (distance[i] < 0.0) {
					distance[i] = 1.0e10;
				}
			}
			for (integer i = 0; i < t.maxelm; i++) {
				myF[i] += 1.0 / (M_PI * distance[i] * distance[i]);
			}
		}
	}

#ifdef _OPENMP 
	omp_set_num_threads(1); // установка одного потока.
#endif

	delete[] distance;
	delete[] stack_ray_trayser;
	delete[] size_cell;

} // calculate_light_flux


#endif
