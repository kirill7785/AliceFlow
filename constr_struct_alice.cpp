#pragma once
#ifndef _CONSTR_STRUCT_CPP_ALICE_
#define _CONSTR_STRUCT_CPP_ALICE_ 1

// Не существующий узел. 
const integer NON_EXISTENT_NODE = - 1;

// 13 апреля 2015 года.
// вычисление градиентов Температуры с помощью теоремы Грина-Гаусса. 
void green_gaussTemperature(integer iP, doublereal* &potent, integer** &nvtx, TOCHKA* &pa,
	ALICE_PARTITION** &sosedi, integer maxelm, bool bbond,
	BOUND* &sosedb, doublereal* &Tx, doublereal* &Ty, doublereal* &Tz,
	integer *ilevel_alice);

// Модуль constr_struct_alice.cpp повторяет функционал модуля constr_struct, НО с учётом использования АЛИС сетки.
// начало разработки 9 сентября 2016. На данный момент уже 21 сентября 2016. На данный момент разработка не завершена 
// но сделан хороший задел и в основном идет уточнение деталей алгоритма которые приводят к существенным правкам кода.
// 21 сентября 2016 10911 строк кода.
// 22 сентября 2016 21529 строк кода. 
// 25 сентября 2016. 
// Учтено что есть flow.ptr и temp.ptr. На данный момент заполнено и flow.ptr и temp.ptr.
// 22237 строк кода.
// 26 сентября 2016. Не тестировалось. 
// 01.10.2016 Учтено iI1 и iI2 для внутреннего плоского источника тепла.
// 22484 строк кода. Первый успешный запуск данных алгоритмов 28 сентября 2016 года.

// Вычисляет количество элементов модели.
void calculate_max_elm(octTree* &oc, integer &maxelm, integer iflag, BLOCK* b, integer lb, bool binit_TD) {

	// Если binit_TD true то нужно обнулять octree1->inum_TD = 0; Иначе мы его не трогаем.

	maxelm = 0; // инициализация.
	top_ALICE_STACK = 0;
	if (oc->link0 != NULL) {
		my_ALICE_STACK[top_ALICE_STACK].link = (oc->link0);
		my_ALICE_STACK[top_ALICE_STACK].maxx = oc->link0->maxx;
		my_ALICE_STACK[top_ALICE_STACK].minx = oc->link0->minx;
		my_ALICE_STACK[top_ALICE_STACK].maxy = oc->link0->maxy;
		my_ALICE_STACK[top_ALICE_STACK].miny = oc->link0->miny;
		my_ALICE_STACK[top_ALICE_STACK].maxz = oc->link0->maxz;
		my_ALICE_STACK[top_ALICE_STACK].minz = oc->link0->minz;
		top_ALICE_STACK++;
	}
	if (oc->link1 != NULL) {
		my_ALICE_STACK[top_ALICE_STACK].link = (oc->link1);
		my_ALICE_STACK[top_ALICE_STACK].maxx = oc->link1->maxx;
		my_ALICE_STACK[top_ALICE_STACK].minx = oc->link1->minx;
		my_ALICE_STACK[top_ALICE_STACK].maxy = oc->link1->maxy;
		my_ALICE_STACK[top_ALICE_STACK].miny = oc->link1->miny;
		my_ALICE_STACK[top_ALICE_STACK].maxz = oc->link1->maxz;
		my_ALICE_STACK[top_ALICE_STACK].minz = oc->link1->minz;
		top_ALICE_STACK++;
	}
	if (oc->link2 != NULL) {
		my_ALICE_STACK[top_ALICE_STACK].link = (oc->link2);
		my_ALICE_STACK[top_ALICE_STACK].maxx = oc->link2->maxx;
		my_ALICE_STACK[top_ALICE_STACK].minx = oc->link2->minx;
		my_ALICE_STACK[top_ALICE_STACK].maxy = oc->link2->maxy;
		my_ALICE_STACK[top_ALICE_STACK].miny = oc->link2->miny;
		my_ALICE_STACK[top_ALICE_STACK].maxz = oc->link2->maxz;
		my_ALICE_STACK[top_ALICE_STACK].minz = oc->link2->minz;
		top_ALICE_STACK++;
	}
	if (oc->link3 != NULL) {
		my_ALICE_STACK[top_ALICE_STACK].link = (oc->link3);
		my_ALICE_STACK[top_ALICE_STACK].maxx = oc->link3->maxx;
		my_ALICE_STACK[top_ALICE_STACK].minx = oc->link3->minx;
		my_ALICE_STACK[top_ALICE_STACK].maxy = oc->link3->maxy;
		my_ALICE_STACK[top_ALICE_STACK].miny = oc->link3->miny;
		my_ALICE_STACK[top_ALICE_STACK].maxz = oc->link3->maxz;
		my_ALICE_STACK[top_ALICE_STACK].minz = oc->link3->minz;
		top_ALICE_STACK++;
	}
	if (oc->link4 != NULL) {
		my_ALICE_STACK[top_ALICE_STACK].link = (oc->link4);
		my_ALICE_STACK[top_ALICE_STACK].maxx = oc->link4->maxx;
		my_ALICE_STACK[top_ALICE_STACK].minx = oc->link4->minx;
		my_ALICE_STACK[top_ALICE_STACK].maxy = oc->link4->maxy;
		my_ALICE_STACK[top_ALICE_STACK].miny = oc->link4->miny;
		my_ALICE_STACK[top_ALICE_STACK].maxz = oc->link4->maxz;
		my_ALICE_STACK[top_ALICE_STACK].minz = oc->link4->minz;
		top_ALICE_STACK++;
	}
	if (oc->link5 != NULL) {
		my_ALICE_STACK[top_ALICE_STACK].link = (oc->link5);
		my_ALICE_STACK[top_ALICE_STACK].maxx = oc->link5->maxx;
		my_ALICE_STACK[top_ALICE_STACK].minx = oc->link5->minx;
		my_ALICE_STACK[top_ALICE_STACK].maxy = oc->link5->maxy;
		my_ALICE_STACK[top_ALICE_STACK].miny = oc->link5->miny;
		my_ALICE_STACK[top_ALICE_STACK].maxz = oc->link5->maxz;
		my_ALICE_STACK[top_ALICE_STACK].minz = oc->link5->minz;
		top_ALICE_STACK++;
	}
	if (oc->link6 != NULL) {
		my_ALICE_STACK[top_ALICE_STACK].link = (oc->link6);
		my_ALICE_STACK[top_ALICE_STACK].maxx = oc->link6->maxx;
		my_ALICE_STACK[top_ALICE_STACK].minx = oc->link6->minx;
		my_ALICE_STACK[top_ALICE_STACK].maxy = oc->link6->maxy;
		my_ALICE_STACK[top_ALICE_STACK].miny = oc->link6->miny;
		my_ALICE_STACK[top_ALICE_STACK].maxz = oc->link6->maxz;
		my_ALICE_STACK[top_ALICE_STACK].minz = oc->link6->minz;
		top_ALICE_STACK++;
	}
	if (oc->link7 != NULL) {
		my_ALICE_STACK[top_ALICE_STACK].link = (oc->link7);
		my_ALICE_STACK[top_ALICE_STACK].maxx = oc->link7->maxx;
		my_ALICE_STACK[top_ALICE_STACK].minx = oc->link7->minx;
		my_ALICE_STACK[top_ALICE_STACK].maxy = oc->link7->maxy;
		my_ALICE_STACK[top_ALICE_STACK].miny = oc->link7->miny;
		my_ALICE_STACK[top_ALICE_STACK].maxz = oc->link7->maxz;
		my_ALICE_STACK[top_ALICE_STACK].minz = oc->link7->minz;
		top_ALICE_STACK++;
	}
	while (top_ALICE_STACK > 0) {
		if (my_ALICE_STACK[top_ALICE_STACK - 1].link != NULL) {
			if (my_ALICE_STACK[top_ALICE_STACK - 1].link->dlist == true) {
				// Гасим информацию о посещениях.
				octTree* octree1 = my_ALICE_STACK[top_ALICE_STACK - 1].link;
				TOCHKA p;
				p.x = 0.125*(octree1->p0.x + octree1->p1.x + octree1->p2.x + octree1->p3.x + octree1->p4.x + octree1->p5.x + octree1->p6.x + octree1->p7.x);
				p.y = 0.125*(octree1->p0.y + octree1->p1.y + octree1->p2.y + octree1->p3.y + octree1->p4.y + octree1->p5.y + octree1->p6.y + octree1->p7.y);
				p.z = 0.125*(octree1->p0.z + octree1->p1.z + octree1->p2.z + octree1->p3.z + octree1->p4.z + octree1->p5.z + octree1->p6.z + octree1->p7.z);
				integer ib;
				bool inDomain = false;
				switch (iflag) {
				case TEMPERATURE: inDomain = in_model_temp(p, ib, b, lb); break;
				case HYDRODINAMIC: inDomain = in_model_flow(p, ib, b, lb);break;
				}
				if (inDomain) {
					// Точка вышла за границы кабинета.
					if ((p.x < b[0].g.xS) || (p.x > b[0].g.xE) || (p.y < b[0].g.yS) || (p.y > b[0].g.yE) || (p.z < b[0].g.zS) || (p.z > b[0].g.zE)) {
						printf("ERROR!!! ib=%lld,%e %e %e\n", ib, p.x, p.y, p.z);
						printf("%e %e %e %e %e %e %e %e \n", octree1->p0.x, octree1->p1.x, octree1->p2.x, octree1->p3.x, octree1->p4.x, octree1->p5.x, octree1->p6.x, octree1->p7.x);
						//system("PAUSE");
						// Причина этого бага в дефектности массивов xpos, ypos или zpos. Помоему это проявляется на этапе coarcemeshgen.
						// Здесь эти сообщения об ошибках погашены.
					}
				}
				if (inDomain) {
					if ((p.x >= b[0].g.xS) && (p.x <= b[0].g.xE) && (p.y >= b[0].g.yS) && (p.y <= b[0].g.yE) && (p.z >= b[0].g.zS) && (p.z <= b[0].g.zE)) {
						// Точка точно внутри кабинета.
						maxelm++;
					}
				}
				if (binit_TD) {
					octree1->inum_TD = 0; // По умолчанию не принадлежит расчётной области.
				}
				octree1->inum_FD = 0; // По умолчанию не принадлежит расчётной области.
				octree1 = NULL;
				my_ALICE_STACK[top_ALICE_STACK - 1].link = NULL;
				top_ALICE_STACK--;
			}
			else {
				// продолжаем добираться до листьев.
				STACK_ALICE buf1 = my_ALICE_STACK[top_ALICE_STACK - 1];
				STACK_ALICE* buf = &buf1;
				top_ALICE_STACK--;
				if (buf->link->link0 != NULL) {
					my_ALICE_STACK[top_ALICE_STACK].link = (buf->link->link0);
					my_ALICE_STACK[top_ALICE_STACK].maxx = buf->link->link0->maxx;
					my_ALICE_STACK[top_ALICE_STACK].minx = buf->link->link0->minx;
					my_ALICE_STACK[top_ALICE_STACK].maxy = buf->link->link0->maxy;
					my_ALICE_STACK[top_ALICE_STACK].miny = buf->link->link0->miny;
					my_ALICE_STACK[top_ALICE_STACK].maxz = buf->link->link0->maxz;
					my_ALICE_STACK[top_ALICE_STACK].minz = buf->link->link0->minz;
					top_ALICE_STACK++;
				}
				if (buf->link->link1 != NULL) {
					my_ALICE_STACK[top_ALICE_STACK].link = (buf->link->link1);
					my_ALICE_STACK[top_ALICE_STACK].maxx = buf->link->link1->maxx;
					my_ALICE_STACK[top_ALICE_STACK].minx = buf->link->link1->minx;
					my_ALICE_STACK[top_ALICE_STACK].maxy = buf->link->link1->maxy;
					my_ALICE_STACK[top_ALICE_STACK].miny = buf->link->link1->miny;
					my_ALICE_STACK[top_ALICE_STACK].maxz = buf->link->link1->maxz;
					my_ALICE_STACK[top_ALICE_STACK].minz = buf->link->link1->minz;
					top_ALICE_STACK++;
				}
				if (buf->link->link2 != NULL) {
					my_ALICE_STACK[top_ALICE_STACK].link = (buf->link->link2);
					my_ALICE_STACK[top_ALICE_STACK].maxx = buf->link->link2->maxx;
					my_ALICE_STACK[top_ALICE_STACK].minx = buf->link->link2->minx;
					my_ALICE_STACK[top_ALICE_STACK].maxy = buf->link->link2->maxy;
					my_ALICE_STACK[top_ALICE_STACK].miny = buf->link->link2->miny;
					my_ALICE_STACK[top_ALICE_STACK].maxz = buf->link->link2->maxz;
					my_ALICE_STACK[top_ALICE_STACK].minz = buf->link->link2->minz;
					top_ALICE_STACK++;
				}
				if (buf->link->link3 != NULL) {
					my_ALICE_STACK[top_ALICE_STACK].link = (buf->link->link3);
					my_ALICE_STACK[top_ALICE_STACK].maxx = buf->link->link3->maxx;
					my_ALICE_STACK[top_ALICE_STACK].minx = buf->link->link3->minx;
					my_ALICE_STACK[top_ALICE_STACK].maxy = buf->link->link3->maxy;
					my_ALICE_STACK[top_ALICE_STACK].miny = buf->link->link3->miny;
					my_ALICE_STACK[top_ALICE_STACK].maxz = buf->link->link3->maxz;
					my_ALICE_STACK[top_ALICE_STACK].minz = buf->link->link3->minz;
					top_ALICE_STACK++;
				}
				if (buf->link->link4 != NULL) {
					my_ALICE_STACK[top_ALICE_STACK].link = (buf->link->link4);
					my_ALICE_STACK[top_ALICE_STACK].maxx = buf->link->link4->maxx;
					my_ALICE_STACK[top_ALICE_STACK].minx = buf->link->link4->minx;
					my_ALICE_STACK[top_ALICE_STACK].maxy = buf->link->link4->maxy;
					my_ALICE_STACK[top_ALICE_STACK].miny = buf->link->link4->miny;
					my_ALICE_STACK[top_ALICE_STACK].maxz = buf->link->link4->maxz;
					my_ALICE_STACK[top_ALICE_STACK].minz = buf->link->link4->minz;
					top_ALICE_STACK++;
				}
				if (buf->link->link5 != NULL) {
					my_ALICE_STACK[top_ALICE_STACK].link = (buf->link->link5);
					my_ALICE_STACK[top_ALICE_STACK].maxx = buf->link->link5->maxx;
					my_ALICE_STACK[top_ALICE_STACK].minx = buf->link->link5->minx;
					my_ALICE_STACK[top_ALICE_STACK].maxy = buf->link->link5->maxy;
					my_ALICE_STACK[top_ALICE_STACK].miny = buf->link->link5->miny;
					my_ALICE_STACK[top_ALICE_STACK].maxz = buf->link->link5->maxz;
					my_ALICE_STACK[top_ALICE_STACK].minz = buf->link->link5->minz;
					top_ALICE_STACK++;
				}
				if (buf->link->link6 != NULL) {
					my_ALICE_STACK[top_ALICE_STACK].link = (buf->link->link6);
					my_ALICE_STACK[top_ALICE_STACK].maxx = buf->link->link6->maxx;
					my_ALICE_STACK[top_ALICE_STACK].minx = buf->link->link6->minx;
					my_ALICE_STACK[top_ALICE_STACK].maxy = buf->link->link6->maxy;
					my_ALICE_STACK[top_ALICE_STACK].miny = buf->link->link6->miny;
					my_ALICE_STACK[top_ALICE_STACK].maxz = buf->link->link6->maxz;
					my_ALICE_STACK[top_ALICE_STACK].minz = buf->link->link6->minz;
					top_ALICE_STACK++;
				}
				if (buf->link->link7 != NULL) {
					my_ALICE_STACK[top_ALICE_STACK].link = (buf->link->link7);
					my_ALICE_STACK[top_ALICE_STACK].maxx = buf->link->link7->maxx;
					my_ALICE_STACK[top_ALICE_STACK].minx = buf->link->link7->minx;
					my_ALICE_STACK[top_ALICE_STACK].maxy = buf->link->link7->maxy;
					my_ALICE_STACK[top_ALICE_STACK].miny = buf->link->link7->miny;
					my_ALICE_STACK[top_ALICE_STACK].maxz = buf->link->link7->maxz;
					my_ALICE_STACK[top_ALICE_STACK].minz = buf->link->link7->minz;
					top_ALICE_STACK++;
				}
			}
		}
		//}
		//getchar();
	}
} // calculate_max_elm


// визуализация в tecplot 360 с учётом hollow блоков в программной модели.
// Построение nodes, nvtx, prop. Частей 1..6 в программной модели.
void constr_nodes_nvtx_prop_alice(octTree* &oc, integer inx, integer iny, integer inz, integer &maxelm, doublereal* &xpos, doublereal* &ypos, doublereal* &zpos,
	integer iflag, BLOCK* b, integer lb, integer* &whot_is_block, TOCHKA* &pa, integer &maxnode, integer** &nvtx, 
	doublereal** &prop, doublereal* &Sc, integer* &ipower_time_depend, TPROP* matlist, integer* &ilevel_alice) {

	integer maxelm_loc = (inx + 1)*(iny + 1)*(inz + 1);
	// Вычисление maxelm.
	calculate_max_elm(oc, maxelm, iflag, b, lb, true);
#if doubleintprecision == 1
	//printf("maxelm=%lld inx=%lld iny=%lld inz=%lld all_in=%lld\n",maxelm,inx,iny,inz,inx*iny*inz);
#else
	//printf("maxelm=%d inx=%d iny=%d inz=%d all_in=%d\n",maxelm,inx,iny,inz,inx*iny*inz);
#endif
	
	//system("PAUSE");
	
	ilevel_alice = NULL;
	ilevel_alice = new integer[maxelm];
	if (ilevel_alice == NULL) {
		// недостаточно памяти на данном оборудовании.
		printf("Problem : not enough memory on your equipment for ilevel_alice constr struct_alice...\n");
		printf("Please any key to exit...\n");
		exit(1);
	}

	whot_is_block = NULL;
	whot_is_block = new integer[maxelm];
	if (whot_is_block == NULL) {
		// недостаточно памяти на данном оборудовании.
		printf("Problem : not enough memory on your equipment for whot_is_block constr struct_alice...\n");
		printf("Please any key to exit...\n");
		exit(1);
	}

	prop = NULL;
	prop = new doublereal*[9];
	if (prop == NULL) {
		// недостаточно памяти на данном оборудовании.
		printf("Problem : not enough memory on your equipment for prop constr struct_alice...\n");
		printf("Please any key to exit...\n");
		//getchar();
		system("pause");
		exit(1);
	}
	for (integer i = 0; i<9; i++) prop[i] = NULL;
	for (integer i = 0; i<9; i++) {
		prop[i] = new doublereal[maxelm];
		if (prop[i] == NULL) {
			// недостаточно памяти на данном оборудовании.
#if doubleintprecision == 1
			printf("Problem : not enough memory on your equipment for prop[%lld] constr struct_alice...\n", i);
#else
			printf("Problem : not enough memory on your equipment for prop[%d] constr struct_alice...\n", i);
#endif

			
			printf("Please any key to exit...\n");
			//getchar();
			system("pause");
			exit(1);
		}
	}
	Sc = NULL;
	Sc = new doublereal[maxelm];
	if (Sc == NULL) {
		// недостаточно памяти на данном оборудовании.
		printf("Problem : not enough memory on your equipment for Sc constr struct_alice...\n");
		printf("Please any key to exit...\n");
		//getchar();
		system("pause");
		exit(1);
	}
	ipower_time_depend = NULL;
	ipower_time_depend = new integer[maxelm];
	if (ipower_time_depend == NULL) {
		// недостаточно памяти на данном оборудовании.
		printf("Problem : not enough memory on your equipment for ipower_time_depend constr struct_alice...\n");
		printf("Please any key to exit...\n");
		//getchar();
		system("pause");
		exit(1);
	}


	// Вычисление допуска.
	doublereal epsTolx = 1.0e40;
	doublereal epsToly = 1.0e40;
	doublereal epsTolz = 1.0e40;
	const doublereal mdop = 0.75;
	for (integer i = 0; i < inx; i++) {
		if (fabs(xpos[i + 1] - xpos[i]) < epsTolx) {
			epsTolx = mdop*fabs(xpos[i + 1] - xpos[i]);
		}
	}
	for (integer i = 0; i < iny; i++) {
		if (fabs(ypos[i + 1] - ypos[i]) < epsToly) {
			epsToly = mdop*fabs(ypos[i + 1] - ypos[i]);
		}
	}
	for (integer i = 0; i < inz; i++) {
		if (fabs(zpos[i + 1] - zpos[i]) < epsTolz) {
			epsTolz = mdop*fabs(zpos[i + 1] - zpos[i]);
		}
	}

	printf("geometric precision tolerance : epsTolx=%e epsToly=%e epsTolz=%e\n", epsTolx, epsToly, epsTolz);
	//system("PAUSE");

	
	
	HASH_POLE* hash_table_export = new HASH_POLE[(inx + 1)*(iny + 1)*(inz + 1)];
	for (integer i_1 = 0; i_1 < (inx + 1)*(iny + 1)*(inz + 1); i_1++) {
		hash_table_export[i_1].flag = false;
		hash_table_export[i_1].inum = NON_EXISTENT_NODE;
	}


	// Вычисление необходимого объёма оперативной памяти.
	// 25 марта 2017

	integer marker_pa_local = 0;

	top_ALICE_STACK = 0;
	if (oc->link0 != NULL) {
		my_ALICE_STACK[top_ALICE_STACK].link = (oc->link0);
		my_ALICE_STACK[top_ALICE_STACK].maxx = oc->link0->maxx;
		my_ALICE_STACK[top_ALICE_STACK].minx = oc->link0->minx;
		my_ALICE_STACK[top_ALICE_STACK].maxy = oc->link0->maxy;
		my_ALICE_STACK[top_ALICE_STACK].miny = oc->link0->miny;
		my_ALICE_STACK[top_ALICE_STACK].maxz = oc->link0->maxz;
		my_ALICE_STACK[top_ALICE_STACK].minz = oc->link0->minz;
		top_ALICE_STACK++;
	}
	if (oc->link1 != NULL) {
		my_ALICE_STACK[top_ALICE_STACK].link = (oc->link1);
		my_ALICE_STACK[top_ALICE_STACK].maxx = oc->link1->maxx;
		my_ALICE_STACK[top_ALICE_STACK].minx = oc->link1->minx;
		my_ALICE_STACK[top_ALICE_STACK].maxy = oc->link1->maxy;
		my_ALICE_STACK[top_ALICE_STACK].miny = oc->link1->miny;
		my_ALICE_STACK[top_ALICE_STACK].maxz = oc->link1->maxz;
		my_ALICE_STACK[top_ALICE_STACK].minz = oc->link1->minz;
		top_ALICE_STACK++;
	}
	if (oc->link2 != NULL) {
		my_ALICE_STACK[top_ALICE_STACK].link = (oc->link2);
		my_ALICE_STACK[top_ALICE_STACK].maxx = oc->link2->maxx;
		my_ALICE_STACK[top_ALICE_STACK].minx = oc->link2->minx;
		my_ALICE_STACK[top_ALICE_STACK].maxy = oc->link2->maxy;
		my_ALICE_STACK[top_ALICE_STACK].miny = oc->link2->miny;
		my_ALICE_STACK[top_ALICE_STACK].maxz = oc->link2->maxz;
		my_ALICE_STACK[top_ALICE_STACK].minz = oc->link2->minz;
		top_ALICE_STACK++;
	}
	if (oc->link3 != NULL) {
		my_ALICE_STACK[top_ALICE_STACK].link = (oc->link3);
		my_ALICE_STACK[top_ALICE_STACK].maxx = oc->link3->maxx;
		my_ALICE_STACK[top_ALICE_STACK].minx = oc->link3->minx;
		my_ALICE_STACK[top_ALICE_STACK].maxy = oc->link3->maxy;
		my_ALICE_STACK[top_ALICE_STACK].miny = oc->link3->miny;
		my_ALICE_STACK[top_ALICE_STACK].maxz = oc->link3->maxz;
		my_ALICE_STACK[top_ALICE_STACK].minz = oc->link3->minz;
		top_ALICE_STACK++;
	}
	if (oc->link4 != NULL) {
		my_ALICE_STACK[top_ALICE_STACK].link = (oc->link4);
		my_ALICE_STACK[top_ALICE_STACK].maxx = oc->link4->maxx;
		my_ALICE_STACK[top_ALICE_STACK].minx = oc->link4->minx;
		my_ALICE_STACK[top_ALICE_STACK].maxy = oc->link4->maxy;
		my_ALICE_STACK[top_ALICE_STACK].miny = oc->link4->miny;
		my_ALICE_STACK[top_ALICE_STACK].maxz = oc->link4->maxz;
		my_ALICE_STACK[top_ALICE_STACK].minz = oc->link4->minz;
		top_ALICE_STACK++;
	}
	if (oc->link5 != NULL) {
		my_ALICE_STACK[top_ALICE_STACK].link = (oc->link5);
		my_ALICE_STACK[top_ALICE_STACK].maxx = oc->link5->maxx;
		my_ALICE_STACK[top_ALICE_STACK].minx = oc->link5->minx;
		my_ALICE_STACK[top_ALICE_STACK].maxy = oc->link5->maxy;
		my_ALICE_STACK[top_ALICE_STACK].miny = oc->link5->miny;
		my_ALICE_STACK[top_ALICE_STACK].maxz = oc->link5->maxz;
		my_ALICE_STACK[top_ALICE_STACK].minz = oc->link5->minz;
		top_ALICE_STACK++;
	}
	if (oc->link6 != NULL) {
		my_ALICE_STACK[top_ALICE_STACK].link = (oc->link6);
		my_ALICE_STACK[top_ALICE_STACK].maxx = oc->link6->maxx;
		my_ALICE_STACK[top_ALICE_STACK].minx = oc->link6->minx;
		my_ALICE_STACK[top_ALICE_STACK].maxy = oc->link6->maxy;
		my_ALICE_STACK[top_ALICE_STACK].miny = oc->link6->miny;
		my_ALICE_STACK[top_ALICE_STACK].maxz = oc->link6->maxz;
		my_ALICE_STACK[top_ALICE_STACK].minz = oc->link6->minz;
		top_ALICE_STACK++;
	}
	if (oc->link7 != NULL) {
		my_ALICE_STACK[top_ALICE_STACK].link = (oc->link7);
		my_ALICE_STACK[top_ALICE_STACK].maxx = oc->link7->maxx;
		my_ALICE_STACK[top_ALICE_STACK].minx = oc->link7->minx;
		my_ALICE_STACK[top_ALICE_STACK].maxy = oc->link7->maxy;
		my_ALICE_STACK[top_ALICE_STACK].miny = oc->link7->miny;
		my_ALICE_STACK[top_ALICE_STACK].maxz = oc->link7->maxz;
		my_ALICE_STACK[top_ALICE_STACK].minz = oc->link7->minz;
		top_ALICE_STACK++;
	}


	while (top_ALICE_STACK > 0) {

		if (my_ALICE_STACK[top_ALICE_STACK - 1].link != NULL) {
			if (my_ALICE_STACK[top_ALICE_STACK - 1].link->dlist == true) {
				// Гасим информацию о посещениях.
				octTree* octree1 = my_ALICE_STACK[top_ALICE_STACK - 1].link;

				// это лист update pa.
				integer i0, i1, i2, i3, i4, i5, i6, i7;

				bool bfound = false;

				TOCHKA p;
				p.x = 0.125*(octree1->p0.x + octree1->p1.x + octree1->p2.x + octree1->p3.x + octree1->p4.x + octree1->p5.x + octree1->p6.x + octree1->p7.x);
				p.y = 0.125*(octree1->p0.y + octree1->p1.y + octree1->p2.y + octree1->p3.y + octree1->p4.y + octree1->p5.y + octree1->p6.y + octree1->p7.y);
				p.z = 0.125*(octree1->p0.z + octree1->p1.z + octree1->p2.z + octree1->p3.z + octree1->p4.z + octree1->p5.z + octree1->p6.z + octree1->p7.z);
				integer ib;
				bool inDomain = false;
				switch (iflag) {
				case TEMPERATURE: inDomain = in_model_temp(p, ib, b, lb);
					break;
				case HYDRODINAMIC: inDomain = in_model_flow(p, ib, b, lb);
					break;
				}
				if ((p.x < b[0].g.xS) || (p.x > b[0].g.xE) || (p.y < b[0].g.yS) || (p.y > b[0].g.yE) || (p.z < b[0].g.zS) || (p.z > b[0].g.zE)) {
					// Лежит за пределами кабинета.
					inDomain = false;
				}
				if (inDomain) {

					// Мы будем оперировать только над точувми принадлежащими РО.
					// Т.е. точки в HOLLOW блоках мы игнорируем.
					// Этот кусок кода занесён внутрь 08 ноября 2016.

					bfound = false;
					integer key_now = hash_key_alice33(inx, iny, inz, xpos, ypos, zpos, octree1->p0, epsTolx, epsToly, epsTolz);
					bfound = hash_table_export[key_now].flag;

					if (!bfound) {
						i0 = marker_pa_local;
						
						hash_table_export[key_now].flag = true;
						hash_table_export[key_now].inum = marker_pa_local;
						marker_pa_local++;
					}
					else {
						i0 = hash_table_export[key_now].inum;
					}
					bfound = false;
					key_now = hash_key_alice33(inx, iny, inz, xpos, ypos, zpos, octree1->p1, epsTolx, epsToly, epsTolz);
					bfound = hash_table_export[key_now].flag;
					if (!bfound) {
						i1 = marker_pa_local;
						
						hash_table_export[key_now].flag = true;
						hash_table_export[key_now].inum = marker_pa_local;
						marker_pa_local++;
					}
					else {
						i1 = hash_table_export[key_now].inum;
					}
					bfound = false;
					key_now = hash_key_alice33(inx, iny, inz, xpos, ypos, zpos, octree1->p2, epsTolx, epsToly, epsTolz);
					bfound = hash_table_export[key_now].flag;
					if (!bfound) {
						i2 = marker_pa_local;
						
						hash_table_export[key_now].flag = true;
						hash_table_export[key_now].inum = marker_pa_local;
						marker_pa_local++;
					}
					else {
						i2 = hash_table_export[key_now].inum;
					}
					bfound = false;
					key_now = hash_key_alice33(inx, iny, inz, xpos, ypos, zpos, octree1->p3, epsTolx, epsToly, epsTolz);
					bfound = hash_table_export[key_now].flag;
					if (!bfound) {
						i3 = marker_pa_local;
						
						hash_table_export[key_now].flag = true;
						hash_table_export[key_now].inum = marker_pa_local;
						marker_pa_local++;
					}
					else {
						i3 = hash_table_export[key_now].inum;
					}
					bfound = false;
					key_now = hash_key_alice33(inx, iny, inz, xpos, ypos, zpos, octree1->p4, epsTolx, epsToly, epsTolz);
					bfound = hash_table_export[key_now].flag;
					if (!bfound) {
						i4 = marker_pa_local;
						
						hash_table_export[key_now].flag = true;
						hash_table_export[key_now].inum = marker_pa_local;
						marker_pa_local++;
					}
					else {
						i4 = hash_table_export[key_now].inum;
					}
					bfound = false;
					key_now = hash_key_alice33(inx, iny, inz, xpos, ypos, zpos, octree1->p5, epsTolx, epsToly, epsTolz);
					bfound = hash_table_export[key_now].flag;
					if (!bfound) {
						i5 = marker_pa_local;
						hash_table_export[key_now].flag = true;
						hash_table_export[key_now].inum = marker_pa_local;
						marker_pa_local++;
					}
					else {
						i5 = hash_table_export[key_now].inum;
					}
					bfound = false;
					key_now = hash_key_alice33(inx, iny, inz, xpos, ypos, zpos, octree1->p6, epsTolx, epsToly, epsTolz);
					bfound = hash_table_export[key_now].flag;
					if (!bfound) {
						i6 = marker_pa_local;
						hash_table_export[key_now].flag = true;
						hash_table_export[key_now].inum = marker_pa_local;
						marker_pa_local++;
					}
					else {
						i6 = hash_table_export[key_now].inum;
					}
					bfound = false;
					key_now = hash_key_alice33(inx, iny, inz, xpos, ypos, zpos, octree1->p7, epsTolx, epsToly, epsTolz);
					bfound = hash_table_export[key_now].flag;
					if (!bfound) {
						i7 = marker_pa_local;
						hash_table_export[key_now].flag = true;
						hash_table_export[key_now].inum = marker_pa_local;
						marker_pa_local++;
					}
					else {
						i7 = hash_table_export[key_now].inum;
					}
									
				}

				octree1 = NULL;
				my_ALICE_STACK[top_ALICE_STACK - 1].link = NULL;
				top_ALICE_STACK--;


			}
			else {

				// продолжаем добираться до листьев.
				STACK_ALICE buf1 = my_ALICE_STACK[top_ALICE_STACK - 1];
				STACK_ALICE* buf = &buf1;
				top_ALICE_STACK--;
				if (buf->link->link0 != NULL) {
					my_ALICE_STACK[top_ALICE_STACK].link = (buf->link->link0);
					my_ALICE_STACK[top_ALICE_STACK].maxx = buf->link->link0->maxx;
					my_ALICE_STACK[top_ALICE_STACK].minx = buf->link->link0->minx;
					my_ALICE_STACK[top_ALICE_STACK].maxy = buf->link->link0->maxy;
					my_ALICE_STACK[top_ALICE_STACK].miny = buf->link->link0->miny;
					my_ALICE_STACK[top_ALICE_STACK].maxz = buf->link->link0->maxz;
					my_ALICE_STACK[top_ALICE_STACK].minz = buf->link->link0->minz;
					top_ALICE_STACK++;
				}
				if (buf->link->link1 != NULL) {
					my_ALICE_STACK[top_ALICE_STACK].link = (buf->link->link1);
					my_ALICE_STACK[top_ALICE_STACK].maxx = buf->link->link1->maxx;
					my_ALICE_STACK[top_ALICE_STACK].minx = buf->link->link1->minx;
					my_ALICE_STACK[top_ALICE_STACK].maxy = buf->link->link1->maxy;
					my_ALICE_STACK[top_ALICE_STACK].miny = buf->link->link1->miny;
					my_ALICE_STACK[top_ALICE_STACK].maxz = buf->link->link1->maxz;
					my_ALICE_STACK[top_ALICE_STACK].minz = buf->link->link1->minz;
					top_ALICE_STACK++;
				}
				if (buf->link->link2 != NULL) {
					my_ALICE_STACK[top_ALICE_STACK].link = (buf->link->link2);
					my_ALICE_STACK[top_ALICE_STACK].maxx = buf->link->link2->maxx;
					my_ALICE_STACK[top_ALICE_STACK].minx = buf->link->link2->minx;
					my_ALICE_STACK[top_ALICE_STACK].maxy = buf->link->link2->maxy;
					my_ALICE_STACK[top_ALICE_STACK].miny = buf->link->link2->miny;
					my_ALICE_STACK[top_ALICE_STACK].maxz = buf->link->link2->maxz;
					my_ALICE_STACK[top_ALICE_STACK].minz = buf->link->link2->minz;
					top_ALICE_STACK++;
				}
				if (buf->link->link3 != NULL) {
					my_ALICE_STACK[top_ALICE_STACK].link = (buf->link->link3);
					my_ALICE_STACK[top_ALICE_STACK].maxx = buf->link->link3->maxx;
					my_ALICE_STACK[top_ALICE_STACK].minx = buf->link->link3->minx;
					my_ALICE_STACK[top_ALICE_STACK].maxy = buf->link->link3->maxy;
					my_ALICE_STACK[top_ALICE_STACK].miny = buf->link->link3->miny;
					my_ALICE_STACK[top_ALICE_STACK].maxz = buf->link->link3->maxz;
					my_ALICE_STACK[top_ALICE_STACK].minz = buf->link->link3->minz;
					top_ALICE_STACK++;
				}
				if (buf->link->link4 != NULL) {
					my_ALICE_STACK[top_ALICE_STACK].link = (buf->link->link4);
					my_ALICE_STACK[top_ALICE_STACK].maxx = buf->link->link4->maxx;
					my_ALICE_STACK[top_ALICE_STACK].minx = buf->link->link4->minx;
					my_ALICE_STACK[top_ALICE_STACK].maxy = buf->link->link4->maxy;
					my_ALICE_STACK[top_ALICE_STACK].miny = buf->link->link4->miny;
					my_ALICE_STACK[top_ALICE_STACK].maxz = buf->link->link4->maxz;
					my_ALICE_STACK[top_ALICE_STACK].minz = buf->link->link4->minz;
					top_ALICE_STACK++;
				}
				if (buf->link->link5 != NULL) {
					my_ALICE_STACK[top_ALICE_STACK].link = (buf->link->link5);
					my_ALICE_STACK[top_ALICE_STACK].maxx = buf->link->link5->maxx;
					my_ALICE_STACK[top_ALICE_STACK].minx = buf->link->link5->minx;
					my_ALICE_STACK[top_ALICE_STACK].maxy = buf->link->link5->maxy;
					my_ALICE_STACK[top_ALICE_STACK].miny = buf->link->link5->miny;
					my_ALICE_STACK[top_ALICE_STACK].maxz = buf->link->link5->maxz;
					my_ALICE_STACK[top_ALICE_STACK].minz = buf->link->link5->minz;
					top_ALICE_STACK++;
				}
				if (buf->link->link6 != NULL) {
					my_ALICE_STACK[top_ALICE_STACK].link = (buf->link->link6);
					my_ALICE_STACK[top_ALICE_STACK].maxx = buf->link->link6->maxx;
					my_ALICE_STACK[top_ALICE_STACK].minx = buf->link->link6->minx;
					my_ALICE_STACK[top_ALICE_STACK].maxy = buf->link->link6->maxy;
					my_ALICE_STACK[top_ALICE_STACK].miny = buf->link->link6->miny;
					my_ALICE_STACK[top_ALICE_STACK].maxz = buf->link->link6->maxz;
					my_ALICE_STACK[top_ALICE_STACK].minz = buf->link->link6->minz;
					top_ALICE_STACK++;
				}
				if (buf->link->link7 != NULL) {
					my_ALICE_STACK[top_ALICE_STACK].link = (buf->link->link7);
					my_ALICE_STACK[top_ALICE_STACK].maxx = buf->link->link7->maxx;
					my_ALICE_STACK[top_ALICE_STACK].minx = buf->link->link7->minx;
					my_ALICE_STACK[top_ALICE_STACK].maxy = buf->link->link7->maxy;
					my_ALICE_STACK[top_ALICE_STACK].miny = buf->link->link7->miny;
					my_ALICE_STACK[top_ALICE_STACK].maxz = buf->link->link7->maxz;
					my_ALICE_STACK[top_ALICE_STACK].minz = buf->link->link7->minz;
					top_ALICE_STACK++;
				}

			}
		}

	}
	delete[] hash_table_export;


	// сформировать pa.
	// заодно посчитать количество узловых точек.
	// сформировать nvtx которые ссылаются на pa.
	// визуализировать сетку.
	TOCHKA* pa_alice = NULL;
	//pa_alice = new TOCHKA[(inx + 1)*(iny + 1)*(inz + 1)];
	pa_alice = new TOCHKA[marker_pa_local+2];
	// Оператор new не требует проверки.
	//if (pa_alice == NULL) {
		// недостаточно памяти на данном оборудовании.
		//printf("Problem : not enough memory on your equipment for pa_alice in adaptive_local_refinement_mesh generator...\n");
		//printf("Please any key to exit...\n");
		//exit(1);
	//}
	integer marker_pa = 0;
	// И тут же сразу формируем nvtx:
	nvtx = NULL;
	nvtx = new integer*[8];
	// Оператор new не требует проверки.
	//if (nvtx == NULL) {
		// недостаточно памяти на данном оборудовании.
		//printf("Problem : not enough memory on your equipment for nvtx in adaptive_local_refinement_mesh generator...\n");
		//printf("Please any key to exit...\n");
		//exit(1);
	//}
	for (integer k_1 = 0; k_1 < 8; k_1++) {
		nvtx[k_1] = NULL;
		nvtx[k_1] = new integer[maxelm];
		//nvtx[k_1] = new integer[marker_pa_local + 2];
		// Оператор new не требует проверки.
		//if (nvtx[k_1] == NULL) {
			// недостаточно памяти на данном оборудовании.
//#if doubleintprecision == 1
	//		printf("Problem : not enough memory on your equipment for nvtx[%lld] in adaptive_local_refinement_mesh generator...\n", k_1);
//#else
	//		printf("Problem : not enough memory on your equipment for nvtx[%d] in adaptive_local_refinement_mesh generator...\n", k_1);
//#endif
			
	//		printf("Please any key to exit...\n");
		//	exit(1);
		//}
	}
	integer imarker_nvtx = 0;

	// Конец вычисления необходимого объема оперативной памяти.
	const integer size_HASH_POLE = (inx + 1)*(iny + 1)*(inz + 1);
	hash_table_export = new HASH_POLE[size_HASH_POLE];
	for (integer i_1 = 0; i_1 < size_HASH_POLE; i_1++) {
		hash_table_export[i_1].flag = false;
		hash_table_export[i_1].inum = NON_EXISTENT_NODE;
	}


	top_ALICE_STACK = 0;
	if (oc->link0 != NULL) {
		my_ALICE_STACK[top_ALICE_STACK].link = (oc->link0);
		my_ALICE_STACK[top_ALICE_STACK].maxx = oc->link0->maxx;
		my_ALICE_STACK[top_ALICE_STACK].minx = oc->link0->minx;
		my_ALICE_STACK[top_ALICE_STACK].maxy = oc->link0->maxy;
		my_ALICE_STACK[top_ALICE_STACK].miny = oc->link0->miny;
		my_ALICE_STACK[top_ALICE_STACK].maxz = oc->link0->maxz;
		my_ALICE_STACK[top_ALICE_STACK].minz = oc->link0->minz;
		top_ALICE_STACK++;
	}
	if (oc->link1 != NULL) {
		my_ALICE_STACK[top_ALICE_STACK].link = (oc->link1);
		my_ALICE_STACK[top_ALICE_STACK].maxx = oc->link1->maxx;
		my_ALICE_STACK[top_ALICE_STACK].minx = oc->link1->minx;
		my_ALICE_STACK[top_ALICE_STACK].maxy = oc->link1->maxy;
		my_ALICE_STACK[top_ALICE_STACK].miny = oc->link1->miny;
		my_ALICE_STACK[top_ALICE_STACK].maxz = oc->link1->maxz;
		my_ALICE_STACK[top_ALICE_STACK].minz = oc->link1->minz;
		top_ALICE_STACK++;
	}
	if (oc->link2 != NULL) {
		my_ALICE_STACK[top_ALICE_STACK].link = (oc->link2);
		my_ALICE_STACK[top_ALICE_STACK].maxx = oc->link2->maxx;
		my_ALICE_STACK[top_ALICE_STACK].minx = oc->link2->minx;
		my_ALICE_STACK[top_ALICE_STACK].maxy = oc->link2->maxy;
		my_ALICE_STACK[top_ALICE_STACK].miny = oc->link2->miny;
		my_ALICE_STACK[top_ALICE_STACK].maxz = oc->link2->maxz;
		my_ALICE_STACK[top_ALICE_STACK].minz = oc->link2->minz;
		top_ALICE_STACK++;
	}
	if (oc->link3 != NULL) {
		my_ALICE_STACK[top_ALICE_STACK].link = (oc->link3);
		my_ALICE_STACK[top_ALICE_STACK].maxx = oc->link3->maxx;
		my_ALICE_STACK[top_ALICE_STACK].minx = oc->link3->minx;
		my_ALICE_STACK[top_ALICE_STACK].maxy = oc->link3->maxy;
		my_ALICE_STACK[top_ALICE_STACK].miny = oc->link3->miny;
		my_ALICE_STACK[top_ALICE_STACK].maxz = oc->link3->maxz;
		my_ALICE_STACK[top_ALICE_STACK].minz = oc->link3->minz;
		top_ALICE_STACK++;
	}
	if (oc->link4 != NULL) {
		my_ALICE_STACK[top_ALICE_STACK].link = (oc->link4);
		my_ALICE_STACK[top_ALICE_STACK].maxx = oc->link4->maxx;
		my_ALICE_STACK[top_ALICE_STACK].minx = oc->link4->minx;
		my_ALICE_STACK[top_ALICE_STACK].maxy = oc->link4->maxy;
		my_ALICE_STACK[top_ALICE_STACK].miny = oc->link4->miny;
		my_ALICE_STACK[top_ALICE_STACK].maxz = oc->link4->maxz;
		my_ALICE_STACK[top_ALICE_STACK].minz = oc->link4->minz;
		top_ALICE_STACK++;
	}
	if (oc->link5 != NULL) {
		my_ALICE_STACK[top_ALICE_STACK].link = (oc->link5);
		my_ALICE_STACK[top_ALICE_STACK].maxx = oc->link5->maxx;
		my_ALICE_STACK[top_ALICE_STACK].minx = oc->link5->minx;
		my_ALICE_STACK[top_ALICE_STACK].maxy = oc->link5->maxy;
		my_ALICE_STACK[top_ALICE_STACK].miny = oc->link5->miny;
		my_ALICE_STACK[top_ALICE_STACK].maxz = oc->link5->maxz;
		my_ALICE_STACK[top_ALICE_STACK].minz = oc->link5->minz;
		top_ALICE_STACK++;
	}
	if (oc->link6 != NULL) {
		my_ALICE_STACK[top_ALICE_STACK].link = (oc->link6);
		my_ALICE_STACK[top_ALICE_STACK].maxx = oc->link6->maxx;
		my_ALICE_STACK[top_ALICE_STACK].minx = oc->link6->minx;
		my_ALICE_STACK[top_ALICE_STACK].maxy = oc->link6->maxy;
		my_ALICE_STACK[top_ALICE_STACK].miny = oc->link6->miny;
		my_ALICE_STACK[top_ALICE_STACK].maxz = oc->link6->maxz;
		my_ALICE_STACK[top_ALICE_STACK].minz = oc->link6->minz;
		top_ALICE_STACK++;
	}
	if (oc->link7 != NULL) {
		my_ALICE_STACK[top_ALICE_STACK].link = (oc->link7);
		my_ALICE_STACK[top_ALICE_STACK].maxx = oc->link7->maxx;
		my_ALICE_STACK[top_ALICE_STACK].minx = oc->link7->minx;
		my_ALICE_STACK[top_ALICE_STACK].maxy = oc->link7->maxy;
		my_ALICE_STACK[top_ALICE_STACK].miny = oc->link7->miny;
		my_ALICE_STACK[top_ALICE_STACK].maxz = oc->link7->maxz;
		my_ALICE_STACK[top_ALICE_STACK].minz = oc->link7->minz;
		top_ALICE_STACK++;
	}

	
	while (top_ALICE_STACK > 0) {
		
		if (my_ALICE_STACK[top_ALICE_STACK - 1].link != NULL) {
			if (my_ALICE_STACK[top_ALICE_STACK - 1].link->dlist == true) {
				// Гасим информацию о посещениях.
				octTree* octree1 = my_ALICE_STACK[top_ALICE_STACK - 1].link;

				// это лист update pa.
				integer i0, i1, i2, i3, i4, i5, i6, i7;
				
				bool bfound = false;
				
				TOCHKA p;
				p.x = 0.125*(octree1->p0.x + octree1->p1.x + octree1->p2.x + octree1->p3.x + octree1->p4.x + octree1->p5.x + octree1->p6.x + octree1->p7.x);
				p.y = 0.125*(octree1->p0.y + octree1->p1.y + octree1->p2.y + octree1->p3.y + octree1->p4.y + octree1->p5.y + octree1->p6.y + octree1->p7.y);
				p.z = 0.125*(octree1->p0.z + octree1->p1.z + octree1->p2.z + octree1->p3.z + octree1->p4.z + octree1->p5.z + octree1->p6.z + octree1->p7.z);
				integer ib;
				bool inDomain = false;
				switch (iflag) {
				case TEMPERATURE: inDomain = in_model_temp(p, ib, b, lb);
					break;
				case HYDRODINAMIC: inDomain = in_model_flow(p, ib, b, lb);
					break;
				}
				if ((p.x < b[0].g.xS) || (p.x > b[0].g.xE) || (p.y < b[0].g.yS) || (p.y > b[0].g.yE) || (p.z < b[0].g.zS) || (p.z > b[0].g.zE)) {
					// Лежит за пределами кабинета.
					inDomain = false;
				}
				if (inDomain) {

					// Мы будем оперировать только над точками принадлежащими РО.
					// Т.е. точки в HOLLOW блоках мы игнорируем.
					// Этот кусок кода занесён внутрь 08 ноября 2016.

					bfound = false;
					integer key_now = hash_key_alice33(inx, iny, inz, xpos, ypos, zpos, octree1->p0, epsTolx, epsToly, epsTolz);
					bfound = hash_table_export[key_now].flag;

					if (!bfound) {
						i0 = marker_pa;
						pa_alice[marker_pa] = octree1->p0;
						hash_table_export[key_now].flag = true;
						hash_table_export[key_now].inum = marker_pa;
						marker_pa++;
					}
					else {
						i0 = hash_table_export[key_now].inum;
					}
					bfound = false;
					key_now = hash_key_alice33(inx, iny, inz, xpos, ypos, zpos, octree1->p1, epsTolx, epsToly, epsTolz);
					bfound = hash_table_export[key_now].flag;
					if (!bfound) {
						i1 = marker_pa;
						pa_alice[marker_pa] = octree1->p1;
						hash_table_export[key_now].flag = true;
						hash_table_export[key_now].inum = marker_pa;
						marker_pa++;
					}
					else {
						i1 = hash_table_export[key_now].inum;
					}
					bfound = false;
					key_now = hash_key_alice33(inx, iny, inz, xpos, ypos, zpos, octree1->p2, epsTolx, epsToly, epsTolz);
					bfound = hash_table_export[key_now].flag;
					if (!bfound) {
						i2 = marker_pa;
						pa_alice[marker_pa] = octree1->p2;
						hash_table_export[key_now].flag = true;
						hash_table_export[key_now].inum = marker_pa;
						marker_pa++;
					}
					else {
						i2 = hash_table_export[key_now].inum;
					}
					bfound = false;
					key_now = hash_key_alice33(inx, iny, inz, xpos, ypos, zpos, octree1->p3, epsTolx, epsToly, epsTolz);
					bfound = hash_table_export[key_now].flag;
					if (!bfound) {
						i3 = marker_pa;
						pa_alice[marker_pa] = octree1->p3;
						hash_table_export[key_now].flag = true;
						hash_table_export[key_now].inum = marker_pa;
						marker_pa++;
					}
					else {
						i3 = hash_table_export[key_now].inum;
					}
					bfound = false;
					key_now = hash_key_alice33(inx, iny, inz, xpos, ypos, zpos, octree1->p4, epsTolx, epsToly, epsTolz);
					bfound = hash_table_export[key_now].flag;
					if (!bfound) {
						i4 = marker_pa;
						pa_alice[marker_pa] = octree1->p4;
						hash_table_export[key_now].flag = true;
						hash_table_export[key_now].inum = marker_pa;
						marker_pa++;
					}
					else {
						i4 = hash_table_export[key_now].inum;
					}
					bfound = false;
					key_now = hash_key_alice33(inx, iny, inz, xpos, ypos, zpos, octree1->p5, epsTolx, epsToly, epsTolz);
					bfound = hash_table_export[key_now].flag;
					if (!bfound) {
						i5 = marker_pa;
						pa_alice[marker_pa] = octree1->p5;
						hash_table_export[key_now].flag = true;
						hash_table_export[key_now].inum = marker_pa;
						marker_pa++;
					}
					else {
						i5 = hash_table_export[key_now].inum;
					}
					bfound = false;
					key_now = hash_key_alice33(inx, iny, inz, xpos, ypos, zpos, octree1->p6, epsTolx, epsToly, epsTolz);
					bfound = hash_table_export[key_now].flag;
					if (!bfound) {
						i6 = marker_pa;
						pa_alice[marker_pa] = octree1->p6;
						hash_table_export[key_now].flag = true;
						hash_table_export[key_now].inum = marker_pa;
						marker_pa++;
					}
					else {
						i6 = hash_table_export[key_now].inum;
					}
					bfound = false;
					key_now = hash_key_alice33(inx, iny, inz, xpos, ypos, zpos, octree1->p7, epsTolx, epsToly, epsTolz);
					bfound = hash_table_export[key_now].flag;
					if (!bfound) {
						i7 = marker_pa;
						pa_alice[marker_pa] = octree1->p7;
						hash_table_export[key_now].flag = true;
						hash_table_export[key_now].inum = marker_pa;
						marker_pa++;
					}
					else {
						i7 = hash_table_export[key_now].inum;
					}





					integer l = imarker_nvtx;
					switch (iflag) {
					case TEMPERATURE: prop[RHO][l] = matlist[b[ib].imatid].rho;
						//prop[HEAT_CAPACITY][l] = matlist[b[ib].imatid].cp;
						//prop[LAM][l] = matlist[b[ib].imatid].lam;
						prop[HEAT_CAPACITY][l] = get_lam(matlist[b[ib].imatid].n_cp, matlist[b[ib].imatid].temp_cp, matlist[b[ib].imatid].arr_cp, 25.0);
						prop[LAM][l] = get_lam(matlist[b[ib].imatid].n_lam, matlist[b[ib].imatid].temp_lam, matlist[b[ib].imatid].arr_lam, 25.0);

						prop[MULT_LAM_X][l] = matlist[b[ib].imatid].orthotropy_multiplyer_x;
						prop[MULT_LAM_Y][l] = matlist[b[ib].imatid].orthotropy_multiplyer_y;
						prop[MULT_LAM_Z][l] = matlist[b[ib].imatid].orthotropy_multiplyer_z;
						// Коэффициенты Лямэ.
						prop[MU_LAME][l] = matlist[b[ib].imatid].mu_Lame;
						prop[LAMBDA_LAME][l] = matlist[b[ib].imatid].lambda_Lame;
						prop[BETA_T_MECHANICAL][l] = matlist[b[ib].imatid].beta_t_solid;

						//Sc[l] = b[ib].Sc;
						Sc[l] = get_power(b[ib].n_Sc, b[ib].temp_Sc, b[ib].arr_Sc, 20.0);
						ipower_time_depend[l] = b[ib].ipower_time_depend;
						break;
					case HYDRODINAMIC:prop[RHO][l] = matlist[b[ib].imatid].rho;
						prop[MU][l] = matlist[b[ib].imatid].mu;
						prop[BETA_T][l] = matlist[b[ib].imatid].beta_t;
						break;
					}
					octree1->inum_TD = imarker_nvtx + 1; // Нумерация начинается с 1.

					// Нумерация начинается с единицы .
					nvtx[0][imarker_nvtx] = i0 + 1;
					nvtx[1][imarker_nvtx] = i1 + 1;
					nvtx[2][imarker_nvtx] = i2 + 1;
					nvtx[3][imarker_nvtx] = i3 + 1;
					nvtx[4][imarker_nvtx] = i4 + 1;
					nvtx[5][imarker_nvtx] = i5 + 1;
					nvtx[6][imarker_nvtx] = i6 + 1;
					nvtx[7][imarker_nvtx] = i7 + 1;
					// конструируем whot_is_block.
					whot_is_block[imarker_nvtx] = ib;
					ilevel_alice[imarker_nvtx] = octree1->ilevel;
					imarker_nvtx++;
				}

				octree1 = NULL;
				my_ALICE_STACK[top_ALICE_STACK - 1].link = NULL;
				top_ALICE_STACK--;


			}
			else {

				// продолжаем добираться до листьев.
				STACK_ALICE buf1 = my_ALICE_STACK[top_ALICE_STACK - 1];
				STACK_ALICE* buf = &buf1;
				top_ALICE_STACK--;
				if (buf->link->link0 != NULL) {
					my_ALICE_STACK[top_ALICE_STACK].link = (buf->link->link0);
					my_ALICE_STACK[top_ALICE_STACK].maxx = buf->link->link0->maxx;
					my_ALICE_STACK[top_ALICE_STACK].minx = buf->link->link0->minx;
					my_ALICE_STACK[top_ALICE_STACK].maxy = buf->link->link0->maxy;
					my_ALICE_STACK[top_ALICE_STACK].miny = buf->link->link0->miny;
					my_ALICE_STACK[top_ALICE_STACK].maxz = buf->link->link0->maxz;
					my_ALICE_STACK[top_ALICE_STACK].minz = buf->link->link0->minz;
					top_ALICE_STACK++;
				}
				if (buf->link->link1 != NULL) {
					my_ALICE_STACK[top_ALICE_STACK].link = (buf->link->link1);
					my_ALICE_STACK[top_ALICE_STACK].maxx = buf->link->link1->maxx;
					my_ALICE_STACK[top_ALICE_STACK].minx = buf->link->link1->minx;
					my_ALICE_STACK[top_ALICE_STACK].maxy = buf->link->link1->maxy;
					my_ALICE_STACK[top_ALICE_STACK].miny = buf->link->link1->miny;
					my_ALICE_STACK[top_ALICE_STACK].maxz = buf->link->link1->maxz;
					my_ALICE_STACK[top_ALICE_STACK].minz = buf->link->link1->minz;
					top_ALICE_STACK++;
				}
				if (buf->link->link2 != NULL) {
					my_ALICE_STACK[top_ALICE_STACK].link = (buf->link->link2);
					my_ALICE_STACK[top_ALICE_STACK].maxx = buf->link->link2->maxx;
					my_ALICE_STACK[top_ALICE_STACK].minx = buf->link->link2->minx;
					my_ALICE_STACK[top_ALICE_STACK].maxy = buf->link->link2->maxy;
					my_ALICE_STACK[top_ALICE_STACK].miny = buf->link->link2->miny;
					my_ALICE_STACK[top_ALICE_STACK].maxz = buf->link->link2->maxz;
					my_ALICE_STACK[top_ALICE_STACK].minz = buf->link->link2->minz;
					top_ALICE_STACK++;
				}
				if (buf->link->link3 != NULL) {
					my_ALICE_STACK[top_ALICE_STACK].link = (buf->link->link3);
					my_ALICE_STACK[top_ALICE_STACK].maxx = buf->link->link3->maxx;
					my_ALICE_STACK[top_ALICE_STACK].minx = buf->link->link3->minx;
					my_ALICE_STACK[top_ALICE_STACK].maxy = buf->link->link3->maxy;
					my_ALICE_STACK[top_ALICE_STACK].miny = buf->link->link3->miny;
					my_ALICE_STACK[top_ALICE_STACK].maxz = buf->link->link3->maxz;
					my_ALICE_STACK[top_ALICE_STACK].minz = buf->link->link3->minz;
					top_ALICE_STACK++;
				}
				if (buf->link->link4 != NULL) {
					my_ALICE_STACK[top_ALICE_STACK].link = (buf->link->link4);
					my_ALICE_STACK[top_ALICE_STACK].maxx = buf->link->link4->maxx;
					my_ALICE_STACK[top_ALICE_STACK].minx = buf->link->link4->minx;
					my_ALICE_STACK[top_ALICE_STACK].maxy = buf->link->link4->maxy;
					my_ALICE_STACK[top_ALICE_STACK].miny = buf->link->link4->miny;
					my_ALICE_STACK[top_ALICE_STACK].maxz = buf->link->link4->maxz;
					my_ALICE_STACK[top_ALICE_STACK].minz = buf->link->link4->minz;
					top_ALICE_STACK++;
				}
				if (buf->link->link5 != NULL) {
					my_ALICE_STACK[top_ALICE_STACK].link = (buf->link->link5);
					my_ALICE_STACK[top_ALICE_STACK].maxx = buf->link->link5->maxx;
					my_ALICE_STACK[top_ALICE_STACK].minx = buf->link->link5->minx;
					my_ALICE_STACK[top_ALICE_STACK].maxy = buf->link->link5->maxy;
					my_ALICE_STACK[top_ALICE_STACK].miny = buf->link->link5->miny;
					my_ALICE_STACK[top_ALICE_STACK].maxz = buf->link->link5->maxz;
					my_ALICE_STACK[top_ALICE_STACK].minz = buf->link->link5->minz;
					top_ALICE_STACK++;
				}
				if (buf->link->link6 != NULL) {
					my_ALICE_STACK[top_ALICE_STACK].link = (buf->link->link6);
					my_ALICE_STACK[top_ALICE_STACK].maxx = buf->link->link6->maxx;
					my_ALICE_STACK[top_ALICE_STACK].minx = buf->link->link6->minx;
					my_ALICE_STACK[top_ALICE_STACK].maxy = buf->link->link6->maxy;
					my_ALICE_STACK[top_ALICE_STACK].miny = buf->link->link6->miny;
					my_ALICE_STACK[top_ALICE_STACK].maxz = buf->link->link6->maxz;
					my_ALICE_STACK[top_ALICE_STACK].minz = buf->link->link6->minz;
					top_ALICE_STACK++;
				}
				if (buf->link->link7 != NULL) {
					my_ALICE_STACK[top_ALICE_STACK].link = (buf->link->link7);
					my_ALICE_STACK[top_ALICE_STACK].maxx = buf->link->link7->maxx;
					my_ALICE_STACK[top_ALICE_STACK].minx = buf->link->link7->minx;
					my_ALICE_STACK[top_ALICE_STACK].maxy = buf->link->link7->maxy;
					my_ALICE_STACK[top_ALICE_STACK].miny = buf->link->link7->miny;
					my_ALICE_STACK[top_ALICE_STACK].maxz = buf->link->link7->maxz;
					my_ALICE_STACK[top_ALICE_STACK].minz = buf->link->link7->minz;
					top_ALICE_STACK++;
				}

			}
		}

	}
	delete[] hash_table_export;

	maxnode = marker_pa;
	pa = NULL;
	pa = new TOCHKA[maxnode];
	if (pa == NULL) {
		// недостаточно памяти на данном оборудовании.
		printf("Problem : not enough memory on your equipment for pa constr struct alice...\n");
		printf("Please any key to exit...\n");
		//getchar();
		system("pause");
		exit(1);
	}
	for (integer i = 0; i < maxnode; i++) {
		pa[i] = pa_alice[i];
	}
	delete[] pa_alice;
	pa_alice = NULL;

	// nvtx && pa сформированы, можно экспортировать в tecplot360
	FILE *fp_4 = NULL;
	errno_t err_4=0;
#ifdef MINGW_COMPILLER
	fp_4=fopen64("ALICEFLOW0_24ALICEMESH.PLT", "w");
#else
	err_4 = fopen_s(&fp_4, "ALICEFLOW0_24ALICEMESH.PLT", "w");
#endif

	if ((err_4) != 0) {
		printf("Create File temp Error\n");
		//getchar();
		system("pause");

	}
	else {

		if (fp_4 != NULL) {
			fprintf(fp_4, "TITLE = \"ALICEFLOW0_24\"\n");
			fprintf(fp_4, "VARIABLES = x, y, z\n");
#if doubleintprecision == 1
			fprintf(fp_4, "ZONE T=\"Rampant\", N=%lld, E=%lld, ET=BRICK, F=FEBLOCK\n\n", marker_pa, imarker_nvtx);
#else
			fprintf(fp_4, "ZONE T=\"Rampant\", N=%d, E=%d, ET=BRICK, F=FEBLOCK\n\n", marker_pa, imarker_nvtx);
#endif
			
			// запись x
			for (integer i = 0; i < marker_pa; i++) {
				fprintf(fp_4, "%+.16f ", pa[i].x);
				if (i % 10 == 0) fprintf(fp_4, "\n");
			}
			fprintf(fp_4, "\n");
			// запись y
			for (integer i = 0; i < marker_pa; i++) {
				fprintf(fp_4, "%+.16f ", pa[i].y);
				if (i % 10 == 0) fprintf(fp_4, "\n");
			}
			fprintf(fp_4, "\n");
			// запись z
			for (integer i = 0; i < marker_pa; i++) {
				fprintf(fp_4, "%+.16f ", pa[i].z);
				if (i % 10 == 0) fprintf(fp_4, "\n");
			}
			fprintf(fp_4, "\n");
			for (integer i = 0; i <= imarker_nvtx - 1; i++) {
#if doubleintprecision == 1
				fprintf(fp_4, "%lld %lld %lld %lld %lld %lld %lld %lld \n", nvtx[0][i], nvtx[1][i], nvtx[2][i], nvtx[3][i], nvtx[4][i], nvtx[5][i], nvtx[6][i], nvtx[7][i]);
#else
				fprintf(fp_4, "%d %d %d %d %d %d %d %d \n", nvtx[0][i], nvtx[1][i], nvtx[2][i], nvtx[3][i], nvtx[4][i], nvtx[5][i], nvtx[6][i], nvtx[7][i]);
#endif
				}
			fclose(fp_4);
			//WinExec("C:\\Program Files\\Tecplot\\Tecplot 360 EX 2014 R1\\bin\\tec360.exe ALICEFLOW0_24ALICEMESH.PLT", SW_NORMAL);
		}
	}
} // constr_nodes_nvtx_prop_alice


  // данная параллельная версия метода гаусса не даёт отличий от серийной версии
  // Это было проверено массовым тестированием (см. код eqsolve_simple_gauss).
  // Данная версия правильная.
  // Существует еще одна версия параллельной программы. В книге Анализ алгоритмов вводный курс
  // есть параллельный алгоритмом SIMD. Он пригоден для очень большого числа потока и возможно
  // более хорош чем данный код. 
bool my_version_gauss1(doublereal** &Amatr, integer nodesmatr, doublereal* &bmatr, doublereal* &xmatr, bool bparallel, integer id) {

	// если bparallel==true то параллельный режим работы иначе серийный.
	// можно использовать как то так и другой код они оба правильные.

	integer i = 0, j = 0, k = 0; // счётчики цикла for
	const doublereal epsilon = 1e-30;
	doublereal M, sum, akk;

	//omp_set_num_threads(inumcore); // установка числа потоков

	// приведение к треугольному виду:
	for (k = 0; k<nodesmatr; k++) {
		akk = Amatr[k][k];
		if (fabs(akk)<epsilon) {
			// решение не может быть получено, т.к.
			// на диагонали находится ноль.
			// Это может возникнуть в случае линейно зависимой системы.
#if doubleintprecision == 1
			printf("\nSolution is not exist! Gauss divizion by zero...k=%lld akk=%e id=%lld\n", k, akk, id);
#else
			printf("\nSolution is not exist! Gauss divizion by zero...k=%d akk=%e id=%d\n", k, akk, id);
#endif
			
			for (integer j1 = 0; j1 < nodesmatr; j1++) {
				for (integer j2 = 0; j2 < nodesmatr; j2++) {
					printf("%e ",Amatr[j1][j2]);
				}
				printf("%e ", bmatr[j1]);
				printf("\n");
			}
			system("PAUSE");
			//system("pause");
			//exit(0);
			return false;
		}

		if (bparallel) {
#pragma omp parallel for shared(k, nodesmatr, Amatr, bmatr) private(i,j,M) firstprivate(akk)
			for (i = k + 1; i<nodesmatr; i++) {

				M = Amatr[i][k] / akk;
				for (j = k; j<nodesmatr; j++) {
					Amatr[i][j] -= M * Amatr[k][j];
				}
				bmatr[i] -= M*bmatr[k];
			}
		}
		else {
			for (i = k + 1; i<nodesmatr; i++) {

				M = Amatr[i][k] / akk;
				for (j = k; j<nodesmatr; j++) {
					Amatr[i][j] -= M * Amatr[k][j];
				}
				bmatr[i] -= M*bmatr[k];
			}
		}
	}
	// процесс обратного исключения
	xmatr[nodesmatr - 1] = bmatr[nodesmatr - 1] / Amatr[nodesmatr - 1][nodesmatr - 1];
	for (i = nodesmatr - 2; i >= 0; i--) {
		sum = 0.0;
		if (bparallel) {
#pragma omp parallel for shared(Amatr,xmatr,i,nodesmatr) private(j) reduction (+: sum)
			for (j = i + 1; j<nodesmatr; j++) {
				sum += Amatr[i][j] * xmatr[j];
			}
		}
		else {
			for (j = i + 1; j<nodesmatr; j++) {
				sum += Amatr[i][j] * xmatr[j];
			}
		}
		xmatr[i] = (bmatr[i] - sum) / Amatr[i][i];
		if (xmatr[i] != xmatr[i]) return false;
	}
	return true;
} // my_version_gauss1



  // данная параллельная версия метода гаусса не даёт отличий от серийной версии
  // Это было проверено массовым тестированием (см. код eqsolve_simple_gauss).
  // Данная версия правильная.
  // Существует еще одна версия параллельной программы. В книге Анализ алгоритмов вводный курс
  // есть параллельный алгоритмом SIMD. Он пригоден для очень большого числа потока и возможно
  // более хорош чем данный код. 
bool my_version_gauss2(doublereal **A, integer nodes, doublereal *b, doublereal* &x, bool bparallel) {

	// если bparallel==true то параллельный режим работы иначе серийный.
	// можно использовать как то так и другой код они оба правильные.

	integer i = 0, j = 0, k = 0; // счётчики цикла for
	const doublereal epsilon = 1e-30;
	doublereal M, sum, akk;

	//omp_set_num_threads(inumcore); // установка числа потоков

	// приведение к треугольному виду:
	for (k = 0; k<nodes; k++) {
		akk = A[k][k];
		if (fabs(akk)<epsilon) {
			// решение не может быть получено, т.к.
			// на диагонали находится ноль.
#if doubleintprecision == 1
			printf("\nSolution is not exist! Gauss divizion by zero...k=%lld akk=%e\n", k, akk);
#else
			printf("\nSolution is not exist! Gauss divizion by zero...k=%d akk=%e\n", k, akk);
#endif
			
			system("PAUSE");
			//system("pause");
			//exit(0);
			return false;
		}

		if (bparallel) {
#pragma omp parallel for shared(k, nodes, A, b) private(i,j,M) firstprivate(akk)
			for (i = k + 1; i<nodes; i++) {

				M = A[i][k] / akk;
				for (j = k; j<nodes; j++) {
					A[i][j] -= M * A[k][j];
				}
				b[i] -= M*b[k];
			}
		}
		else {
			for (i = k + 1; i<nodes; i++) {

				M = A[i][k] / akk;
				for (j = k; j<nodes; j++) {
					A[i][j] -= M * A[k][j];
				}
				b[i] -= M*b[k];
			}
		}
	}
	// процесс обратного исключения
	x[nodes - 1] = b[nodes - 1] / A[nodes - 1][nodes - 1];
	for (i = nodes - 2; i >= 0; i--) {
		sum = 0.0;
		if (bparallel) {
#pragma omp parallel for shared(A,x,i,nodes) private(j) reduction (+: sum)
			for (j = i + 1; j<nodes; j++) {
				sum += A[i][j] * x[j];
			}
		}
		else {
			for (j = i + 1; j<nodes; j++) {
				sum += A[i][j] * x[j];
			}
		}
		x[i] = (b[i] - sum) / A[i][i];
		if (x[i] != x[i]) return false;
	}
	return true;
} // my_version_gauss2

// Рефакторинг функции ANES_tecplot360_export_temperature. 02.04.2019.
// Сигнатура вызова:
// ZERO_ORDER_RECONSTRUCT(maxnod, maxelm, pa, nvtx, vol, lam, eps_mashine, f1, f2,f3, false);
void ZERO_ORDER_RECONSTRUCT(integer maxnod, integer maxelm, 
	TOCHKA* &pa, integer** &nvtx, doublereal* &vol, doublereal* &lam, 
	const doublereal eps_mashine, doublereal* f1, doublereal* f2,
	doublereal* f3, integer bfirst_or_sqrt) {
	// Метод нулевого порядка.

	doublereal maximum = -1.0e60;
	doublereal maximum_loc = -1.0e60;

	// TODO Метод повидимому неточный и чтобы увеличить точность надо организовать 
	// преобразование сохраняющее максимальное значение. 02.04.2019

	for (integer i = 0; i < maxelm; i++) {
		if (bfirst_or_sqrt) {
			if (f1[i] > maximum) maximum = f1[i];
		}
		else {
			doublereal d1 = sqrt(f1[i] * f1[i] + f2[i] * f2[i] + f3[i] * f3[i]);
			if (d1 > maximum) maximum = d1;
		}
	}

	for (integer i = 0; i < maxnod; i++) {
		vol[i] = 0.0;
		lam[i] = 0.0;
	}

	for (integer i = 0; i <= maxelm - 1; i++) {
		// вычисление размеров текущего контрольного объёма:
		doublereal dx = 0.0, dy = 0.0, dz = 0.0;// объём текущего контроольного объёма
		volume3D(i, nvtx, pa, dx, dy, dz);

		for (integer j = 0; j <= 7; j++) {
			vol[nvtx[j][i] - 1] += dx * dy*dz;
			if (bfirst_or_sqrt) {
				lam[nvtx[j][i] - 1] += dx * dy*dz*f1[i];
			}
			else {
				lam[nvtx[j][i] - 1] += dx * dy*dz*sqrt(f1[i] * f1[i] + f2[i] * f2[i] + f3[i] * f3[i]);
			}
		}
	}
	for (integer i = 0; i < maxnod; i++) {
		if (fabs(vol[i]) > eps_mashine) {
			lam[i] = lam[i] / vol[i];
			if (lam[i] > maximum_loc)  maximum_loc = lam[i];
		}
		else {
#if doubleintprecision == 1
			printf("fatal error! ANES_tecplot_export in module constr_struct_alice.cpp. i=%lld maxnod=%lld\n", i, maxnod);
			printf("vol[%lld]==%e\n", i, vol[i]);
#else
			printf("fatal error! ANES_tecplot_export in module constr_struct_alice.cpp. i=%d maxnod=%d\n", i, maxnod);
			printf("vol[%d]==%e\n", i, vol[i]);
#endif

			//getchar();
			system("PAUSE");
			lam[i] = 0.0;
		}
	}
	if (fabs(maximum_loc) > 1.0e-20) {
		for (integer i = 0; i < maxnod; i++) {
			// Обеспечиваем сохранение модуля, 
			// т.к. интерполчция сглаживает и 
			// занижает пиковые значения.
			lam[i] *= (maximum) / (maximum_loc);
		}
	}
}// ZERO_ORDER_RECONSTRUCT


 // Рефакторинг функции ANES_tecplot360_export_temperature. 02.04.2019.
 // Сигнатура вызова:
 // ZERO_ORDER_RECONSTRUCT_ORTHOTROPY(maxnod, maxelm, pa, nvtx, vol, lamx, lamy, lamz, eps_mashine, f1, f2,f3, false,fx,fy,fz);
void ZERO_ORDER_RECONSTRUCT_ORTHOTROPY(integer maxnod, integer maxelm,
	TOCHKA* &pa, integer** &nvtx, doublereal* &vol, doublereal* &lamx, doublereal* &lamy, doublereal* &lamz,
	const doublereal eps_mashine, doublereal* f1, doublereal* f2,
	doublereal* f3, integer bfirst_or_sqrt, doublereal* fx, doublereal* fy,
	doublereal* fz ) {
	// Метод нулевого порядка.

	doublereal maximum1 = -1.0e60;
	doublereal maximum_loc1 = -1.0e60;
	doublereal maximum2 = -1.0e60;
	doublereal maximum_loc2 = -1.0e60;
	doublereal maximum3 = -1.0e60;
	doublereal maximum_loc3 = -1.0e60;


	// TODO Метод повидимому неточный и чтобы увеличить точность надо организовать 
	// преобразование сохраняющее максимальное значение. 02.04.2019

	for (integer i = 0; i < maxelm; i++) {
		if (bfirst_or_sqrt) {
			if (f1[i]* fx[i] > maximum1) maximum1 = f1[i] * fx[i];
			if (f1[i] * fy[i] > maximum2) maximum2 = f1[i] * fy[i];
			if (f1[i] * fz[i] > maximum3) maximum3 = f1[i] * fz[i];
		}
		else {
			doublereal d1 = sqrt(f1[i] * f1[i] + f2[i] * f2[i] + f3[i] * f3[i]);
			if (d1 > maximum1) maximum1 = d1;
		}
	}

	for (integer i = 0; i < maxnod; i++) {
		vol[i] = 0.0;
		lamx[i] = 0.0;
		lamy[i] = 0.0;
		lamz[i] = 0.0;
	}

	for (integer i = 0; i <= maxelm - 1; i++) {
		// вычисление размеров текущего контрольного объёма:
		doublereal dx = 0.0, dy = 0.0, dz = 0.0;// объём текущего контроольного объёма
		volume3D(i, nvtx, pa, dx, dy, dz);

		for (integer j = 0; j <= 7; j++) {
			vol[nvtx[j][i] - 1] += dx * dy*dz;
			if (bfirst_or_sqrt) {
				lamx[nvtx[j][i] - 1] += dx * dy * dz*f1[i] * fx[i];
				lamy[nvtx[j][i] - 1] += dx * dy * dz*f1[i] * fy[i];
				lamz[nvtx[j][i] - 1] += dx * dy * dz*f1[i] * fz[i];
			}
			else {
				lamx[nvtx[j][i] - 1] += dx * dy*dz*sqrt(f1[i] * f1[i] + f2[i] * f2[i] + f3[i] * f3[i]);
			}
		}
	}
	for (integer i = 0; i < maxnod; i++) {
		if (fabs(vol[i]) > eps_mashine) {
			lamx[i] = lamx[i] / vol[i];
			lamy[i] = lamy[i] / vol[i];
			lamz[i] = lamz[i] / vol[i];
			if (lamx[i] > maximum_loc1)  maximum_loc1 = lamx[i];
			if (lamy[i] > maximum_loc2)  maximum_loc2 = lamy[i];
			if (lamz[i] > maximum_loc3)  maximum_loc3 = lamz[i];
		}
		else {
#if doubleintprecision == 1
			printf("fatal error! ANES_tecplot_export in module constr_struct_alice.cpp. i=%lld maxnod=%lld\n", i, maxnod);
			printf("vol[%lld]==%e\n", i, vol[i]);
#else
			printf("fatal error! ANES_tecplot_export in module constr_struct_alice.cpp. i=%d maxnod=%d\n", i, maxnod);
			printf("vol[%d]==%e\n", i, vol[i]);
#endif

			//getchar();
			system("PAUSE");
			lamx[i] = 0.0;
			lamy[i] = 0.0;
			lamz[i] = 0.0;
		}
	}
	if (bfirst_or_sqrt) {
		if ((fabs(maximum_loc1) > 1.0e-20)&& (fabs(maximum_loc2) > 1.0e-20) && (fabs(maximum_loc3) > 1.0e-20)) {
			for (integer i = 0; i < maxnod; i++) {
				// Обеспечиваем сохранение модуля, 
				// т.к. интерполчция сглаживает и 
				// занижает пиковые значения.
				lamx[i] *= (maximum1) / (maximum_loc1);
				lamy[i] *= (maximum2) / (maximum_loc2);
				lamz[i] *= (maximum3) / (maximum_loc3);
			}
		}
	}
	else {
		if (fabs(maximum_loc1) > 1.0e-20) {
			for (integer i = 0; i < maxnod; i++) {
				// Обеспечиваем сохранение модуля, 
				// т.к. интерполчция сглаживает и 
				// занижает пиковые значения.
				lamx[i] *= (maximum1) / (maximum_loc1);
			}
		}
	}
}// ZERO_ORDER_RECONSTRUCT_ORTHOTROPY

//fglobal[0].potent[PAM]
// Сигнатура вызова.
// FIRST_ORDER_LINEAR_RECONSTRUCT(fp_4, maxnod, maxelm, pa, nvtx, vol, temp, min_x, min_y, min_z, fglobal[0].potent[PAM], t, eps_mashine);
//
void FIRST_ORDER_LINEAR_RECONSTRUCT(FILE* &fp_4, 
	integer &maxnod, integer &maxelm, TOCHKA* &pa, integer** &nvtx,
	doublereal* &vol, doublereal* &temp, 
	doublereal& min_x, doublereal& min_y, doublereal& min_z, 
	doublereal* &potent, TEMPER &t, const doublereal eps_mashine,
	bool bptr_rule, doublereal* &gradX, doublereal* &gradY, doublereal* &gradZ) {

	

	if ((bptr_rule && bSIMPLErun_now_for_temperature)||(!bptr_rule)) {
		// Для ускорения сканирования в методе наименьших квадратов 8.07.2017.


		doublereal maximum = -1.0e60;
		doublereal maximum1 = -1.0e60;

		integer** q_hash = NULL;
		q_hash = new integer*[maxnod + 1];
		integer* q_ic = NULL;
		q_ic = new integer[maxnod + 1];
		for (integer j = 0; j <= maxnod; j++) {
			q_hash[j] = new integer[9];
			q_ic[j] = 0;
		}
		for (integer j = 0; j <= maxnod; j++) {
			for (integer i_1 = 0; i_1 < 9; i_1++) {
				q_hash[j][i_1] = NON_EXISTENT_NODE;
			}
		}
		for (integer i_1 = 0; i_1 <= maxelm - 1; i_1++) {
			for (integer j_1 = 0; j_1 <= 7; j_1++) {
				q_hash[nvtx[j_1][i_1]][q_ic[nvtx[j_1][i_1]]] = i_1;
				q_ic[nvtx[j_1][i_1]]++;
			}
		}

		integer* inum_now = new integer[maxnod];

		for (integer i = 0; i < maxnod; i++) {
			temp[i] = 0.0;
			vol[i] = 0.0;
			inum_now[i] = 0;
		}

		// Идея расширения шаблона: мы рассматриваем не только текущий nvtx, но и всех
		// его соседей имеющих с ним хоть одну общую вершину. Этим самым шаблон используемый
		// для реконструкции будет расширен новыми точками.
		// Модификация 30.05.2017.
		for (integer i = 0; i <= maxelm - 1; i++) {
			//if (((10 * i) % maxelm) == 0) printf("complete %lld\n", (100 * i / maxelm));
			for (integer j = 0; j <= 7; j++) {
				inum_now[nvtx[j][i] - 1] += 1;
				//for (integer i_1 = 0; i_1 <= maxelm - 1; i_1++) {
				integer i_1 = NON_EXISTENT_NODE;
				for (integer i_2 = 0; i_2 < 9; i_2++) {
					i_1 = q_hash[nvtx[j][i]][i_2];
					if (i_1 >= 0) {
						for (integer j_1 = 0; j_1 <= 7; j_1++) {
							if (i_1 != i) {
								if (nvtx[j][i] == nvtx[j_1][i_1]) {
									inum_now[nvtx[j][i] - 1] += 1;
								}
							}
						}
					}
				}
			}
		}

		for (integer i = 0; i < maxnod; i++) {
			if (inum_now[i] < 1) {
#if doubleintprecision == 1
				printf("i=%lld maxnod=%lld inum_now[%lld]=%lld\n", i, maxnod, i, inum_now[i]);
#else
				printf("i=%d maxnod=%d inum_now[%d]=%d\n", i, maxnod, i, inum_now[i]);
#endif

				system("PAUSE");
			}
		}

		TOCHKA** pointerlist = new TOCHKA*[maxnod];
		doublereal** rthdsd_Gauss = new doublereal*[maxnod];
		for (integer i = 0; i < maxnod; i++) {
			pointerlist[i] = new TOCHKA[(inum_now[i])];
			rthdsd_Gauss[i] = new doublereal[(inum_now[i])];
		}




		for (integer i = 0; i < maxnod; i++) {
			inum_now[i] = 0;
		}

		

		// Непосредственно само вычисление.
		for (integer i = 0; i <= maxelm - 1; i++) {
			//if (((10 * i) % maxelm) == 0) printf("complete %lld\n", (100 * i / maxelm));

			if (bptr_rule) {
				if (potent[t.ptr[ENUMERATECONTVOL][i]] > maximum) {
					maximum = potent[t.ptr[ENUMERATECONTVOL][i]];
				}
			}
			else {
				if (potent[i] > maximum) {
					maximum = potent[i];
				}
			}

			// вычисление размеров текущего контрольного объёма:
			doublereal dx = 0.0, dy = 0.0, dz = 0.0;// объём текущего контрольного объёма
			volume3D(i, nvtx, pa, dx, dy, dz);

			for (integer j = 0; j <= 7; j++) {
				vol[nvtx[j][i] - 1] += dx * dy*dz;
				//temp[nvtx[j][i] - 1] += dx * dy*dz*fglobal[0].potent[PAM][t.ptr[ENUMERATECONTVOL][i]];
			}

			TOCHKA p;
			center_cord3D(i, nvtx, pa, p, 100);
			p.x = p.x + min_x;
			p.y = p.y + min_y;
			p.z = p.z + min_z;


			for (integer j = 0; j <= 7; j++) {
				pointerlist[nvtx[j][i] - 1][inum_now[nvtx[j][i] - 1]] = p;
				if (fabs(p.x) < 1.0e-30) {
					printf("problem x=%e\n", p.x);
					system("PAUSE");
				}
				if (fabs(p.y) < 1.0e-30) {
					printf("problem y=%e\n", p.y);
					system("PAUSE");
				}
				if (fabs(p.z) < 1.0e-30) {
					printf("problem z=%e\n", p.z);
					system("PAUSE");
				}
				//fglobal[0].potent[PAM][t.ptr[ENUMERATECONTVOL][i]]
				if (bptr_rule) {
					rthdsd_Gauss[nvtx[j][i] - 1][inum_now[nvtx[j][i] - 1]] = potent[t.ptr[ENUMERATECONTVOL][i]];// potent[i];
				}
				else {
					rthdsd_Gauss[nvtx[j][i] - 1][inum_now[nvtx[j][i] - 1]] = potent[i];
				}
				inum_now[nvtx[j][i] - 1]++;
			}





			for (integer j = 0; j <= 7; j++) {
				//for (integer i_1 = 0; i_1 <= maxelm - 1; i_1++) {
				integer i_1 = NON_EXISTENT_NODE;
				for (integer i_2 = 0; i_2 < 9; i_2++) {
					i_1 = q_hash[nvtx[j][i]][i_2];
					if (i_1 >= 0) {
						for (integer j_1 = 0; j_1 <= 7; j_1++) {
							if (i_1 != i) {
								if (nvtx[j][i] == nvtx[j_1][i_1]) {
									TOCHKA p_1;
									center_cord3D(i_1, nvtx, pa, p_1, 100);
									p_1.x = p_1.x + min_x;
									p_1.y = p_1.y + min_y;
									p_1.z = p_1.z + min_z;

									pointerlist[nvtx[j][i] - 1][inum_now[nvtx[j][i] - 1]] = p_1;
									if (fabs(p_1.x) < 1.0e-30) {
										printf("problem x=%e\n", p_1.x);
										system("PAUSE");
									}
									if (fabs(p_1.y) < 1.0e-30) {
										printf("problem y=%e\n", p_1.y);
										system("PAUSE");
									}
									if (fabs(p_1.z) < 1.0e-30) {
										printf("problem z=%e\n", p_1.z);
										system("PAUSE");
									}
									//fglobal[0].potent[PAM][t.ptr[ENUMERATECONTVOL][i_1]]
									if (bptr_rule) {
										rthdsd_Gauss[nvtx[j][i] - 1][inum_now[nvtx[j][i] - 1]] = potent[t.ptr[ENUMERATECONTVOL][i_1]];// potent[i_1];
									}
									else {
										rthdsd_Gauss[nvtx[j][i] - 1][inum_now[nvtx[j][i] - 1]] = potent[i_1];
									}
									inum_now[nvtx[j][i] - 1]++;
								}
							}
						}
					}
				}
			}


			/*
			for (integer j = 0; j <= 7; j++) {
			vol[nvtx[j][i] - 1] += dx*dy*dz;
			temp[nvtx[j][i] - 1] += dx*dy*dz*potent[i];
			}
			*/
		}


		for (integer j = 0; j <= maxnod; j++) {
			delete[] q_hash[j];
		}
		delete[] q_ic;
		delete[] q_hash;
		q_ic = NULL;
		q_hash = NULL;

		//integer jcontrol = 0;
		for (integer i = 0; i < maxnod; i++) {
			//if (((10 * i) % maxnod) == 0) printf("complete %lld\n", (100 * i / maxnod));
			if (fabs(vol[i]) > eps_mashine) {


				doublereal** Xmatr = new doublereal*[4];
				for (integer j = 0; j <= 3; j++) {
					Xmatr[j] = new doublereal[4];
				}


				doublereal* bmatr = new doublereal[4];
				doublereal* koefmatr = new doublereal[4];

				for (integer j1 = 0; j1 <= 3; j1++) {
					for (integer j2 = 0; j2 <= 3; j2++) {
						Xmatr[j1][j2] = 0.0;
					}
					bmatr[j1] = 0.0;
					koefmatr[j1] = 0.0;
				}



				for (integer j = 0; j < inum_now[i]; j++) {

					Xmatr[0][0] += 1.0;
					Xmatr[0][1] += pointerlist[i][j].x;
					Xmatr[0][2] += pointerlist[i][j].y;
					Xmatr[0][3] += pointerlist[i][j].z;

					Xmatr[1][0] += pointerlist[i][j].x;
					Xmatr[1][1] += pointerlist[i][j].x*pointerlist[i][j].x;
					Xmatr[1][2] += pointerlist[i][j].x*pointerlist[i][j].y;
					Xmatr[1][3] += pointerlist[i][j].x*pointerlist[i][j].z;

					Xmatr[2][0] += pointerlist[i][j].y;
					Xmatr[2][1] += pointerlist[i][j].y*pointerlist[i][j].x;
					Xmatr[2][2] += pointerlist[i][j].y*pointerlist[i][j].y;
					Xmatr[2][3] += pointerlist[i][j].y*pointerlist[i][j].z;

					Xmatr[3][0] += pointerlist[i][j].z;
					Xmatr[3][1] += pointerlist[i][j].z*pointerlist[i][j].x;
					Xmatr[3][2] += pointerlist[i][j].z*pointerlist[i][j].y;
					Xmatr[3][3] += pointerlist[i][j].z*pointerlist[i][j].z;

					bmatr[0] += rthdsd_Gauss[i][j];
					bmatr[1] += pointerlist[i][j].x*rthdsd_Gauss[i][j];
					bmatr[2] += pointerlist[i][j].y*rthdsd_Gauss[i][j];
					bmatr[3] += pointerlist[i][j].z*rthdsd_Gauss[i][j];

				}

				if ((fabs(Xmatr[0][0]) < 1.e-30) || (fabs(Xmatr[1][1]) < 1.e-30) || (fabs(Xmatr[2][2]) < 1.e-30) || (fabs(Xmatr[3][3]) < 1.e-30)) {
#if doubleintprecision == 1
					printf("inum_now[%lld]=%lld\n", i, inum_now[i]);
#else
					printf("inum_now[%d]=%d\n", i, inum_now[i]);
#endif

					system("PAUSE");
				}

				//Xmatr*koefmatr = bmatr;
				/*
				// Метод Гаусса не работает т.к. система линейно зависима.
				if (!my_version_gauss1(Xmatr, 4, bmatr, koefmatr, false, i)) {
				temp[i] = temp[i] / vol[i];
				}
				else {
				temp[i] = koefmatr[0] + koefmatr[1] * (pa[i].x+min_x) + koefmatr[2] * (pa[i].y+min_y) + koefmatr[3] * (pa[i].z+min_z);
				}
				*/
				for (integer j1 = 0; j1 <= 3; j1++) {
					koefmatr[j1] = 0.0;
				}
				for (integer j1 = 0; j1 <= 250; j1++) {
					doublereal alpha = 0.2;
					doublereal d_0 = koefmatr[0];
					doublereal d_1 = koefmatr[1];
					doublereal d_2 = koefmatr[2];
					doublereal d_3 = koefmatr[3];
					koefmatr[0] = (1.0 - alpha)*d_0 + alpha * ((bmatr[0] - Xmatr[0][1] * koefmatr[1] - Xmatr[0][2] * koefmatr[2] - Xmatr[0][3] * koefmatr[3]) / Xmatr[0][0]);
					koefmatr[1] = (1.0 - alpha)*d_1 + alpha * ((bmatr[1] - Xmatr[1][0] * koefmatr[0] - Xmatr[1][2] * koefmatr[2] - Xmatr[1][3] * koefmatr[3]) / Xmatr[1][1]);
					koefmatr[2] = (1.0 - alpha)*d_2 + alpha * ((bmatr[2] - Xmatr[2][0] * koefmatr[0] - Xmatr[2][1] * koefmatr[1] - Xmatr[2][3] * koefmatr[3]) / Xmatr[2][2]);
					koefmatr[3] = (1.0 - alpha)*d_3 + alpha * ((bmatr[3] - Xmatr[3][0] * koefmatr[0] - Xmatr[3][1] * koefmatr[1] - Xmatr[3][2] * koefmatr[2]) / Xmatr[3][3]);
				}
				temp[i] = koefmatr[0] + koefmatr[1] * (pa[i].x + min_x) + koefmatr[2] * (pa[i].y + min_y) + koefmatr[3] * (pa[i].z + min_z);
				if (temp[i] > maximum1) {
					maximum1 = temp[i];
				}
				//temp[i] = koefmatr[0] + koefmatr[1] * (pa[i].x ) + koefmatr[2] * (pa[i].y ) + koefmatr[3] * (pa[i].z );
				//heat_flux_X[i] = koefmatr[1];
				//heat_flux_Y[i] = koefmatr[2];
				//heat_flux_Z[i] = koefmatr[3];
				gradX[i] = koefmatr[1];
				gradY[i] = koefmatr[2];
				gradZ[i] = koefmatr[3];
				//heat_flux_X[i] = 0.0;
				//heat_flux_Y[i] = 0.0;
				//heat_flux_Z[i] = 0.0;
				// вычисление размеров текущего контрольного объёма:
				/*
				doublereal h_1= 1.0e-4*(max_x-min_x1);
				h_1 = pow(fabs(vol[i]), 0.333);
				//h_1 = 0.5*dx1;
				heat_flux_X[i] = ((koefmatr[0] + koefmatr[1] * (pa[i].x + min_x+h_1) + koefmatr[2] * (pa[i].y + min_y) + koefmatr[3] * (pa[i].z + min_z)
				)-(koefmatr[0] + koefmatr[1] * (pa[i].x + min_x-h_1) + koefmatr[2] * (pa[i].y + min_y) + koefmatr[3] * (pa[i].z + min_z)
				)) / (2 * h_1);
				//h_1 = 1.0e-4*(max_y - min_y1);
				//h_1 = 0.5*dy1;
				heat_flux_Y[i] = ((koefmatr[0] + koefmatr[1] * (pa[i].x + min_x) + koefmatr[2] * (pa[i].y + min_y+h_1) + koefmatr[3] * (pa[i].z + min_z)
				) - (koefmatr[0] + koefmatr[1] * (pa[i].x + min_x) + koefmatr[2] * (pa[i].y + min_y-h_1) + koefmatr[3] * (pa[i].z + min_z)
				)) / (2 * h_1);
				//h_1 = 1.0e-4*(max_z - min_z1);
				//h_1 = 0.5*dz1;
				heat_flux_Z[i] = ((koefmatr[0] + koefmatr[1] * (pa[i].x + min_x) + koefmatr[2] * (pa[i].y + min_y) + koefmatr[3] * (pa[i].z + min_z+h_1)
				) - (koefmatr[0] + koefmatr[1] * (pa[i].x + min_x) + koefmatr[2] * (pa[i].y + min_y) + koefmatr[3] * (pa[i].z + min_z-h_1)
				)) / (2 * h_1);
				*/

				for (integer j = 0; j <= 3; j++) {
					delete[] Xmatr[j];
				}
				delete[] Xmatr;
				delete[] bmatr;
				delete[] koefmatr;

			}
			else {
#if doubleintprecision == 1
				printf("fatal error! ANES_tecplot_export in module constr_struct_alice.cpp. i=%lld maxnod=%lld\n", i, maxnod);
				printf("vol[%lld]==%e\n", i, vol[i]);
#else
				printf("fatal error! ANES_tecplot_export in module constr_struct_alice.cpp. i=%d maxnod=%d\n", i, maxnod);
				printf("vol[%d]==%e\n", i, vol[i]);
#endif

				//getchar();
				system("PAUSE");
				temp[i] = 0.0;
			}



		}

		// При преобразовании сохраняем модуль величины неизменным.
		if (fabs(maximum1) > 1.0e-20) {
			for (integer i = 0; i < maxnod; i++) {
				temp[i] *= (maximum / maximum1);
			}
		}

		delete[] inum_now;
		for (integer i = 0; i < maxnod; i++) {
			delete[] pointerlist[i];
			delete[] rthdsd_Gauss[i];
		}
		delete[] pointerlist;
		delete[] rthdsd_Gauss;



		// запись PAM
		for (integer i = 0; i < maxnod; i++) {
			fprintf(fp_4, "%+.16f ", temp[i]);
			if (i % 10 == 0) fprintf(fp_4, "\n");
		}
		fprintf(fp_4, "\n");
	}
} // FIRST_ORDER_LINEAR_RECONSTRUCT

// 05.04.2019
  //fglobal[0].potent[PAM]
  // Сигнатура вызова.
  // SECOND_ORDER_QUADRATIC_RECONSTRUCT(fp_4, maxnod, maxelm, pa, nvtx, vol, temp, min_x, min_y, min_z, fglobal[0].potent[PAM], t, eps_mashine);
  //
void SECOND_ORDER_QUADRATIC_RECONSTRUCT(FILE* &fp_4,
	integer &maxnod, integer &maxelm, TOCHKA* &pa, integer** &nvtx,
	doublereal* &vol, doublereal* &temp,
	doublereal& min_x, doublereal& min_y, doublereal& min_z,
	doublereal* &potent, TEMPER &t, const doublereal eps_mashine,
	bool bptr_rule, doublereal* &gradX, doublereal* &gradY, doublereal* &gradZ) {



	if ((bptr_rule && bSIMPLErun_now_for_temperature) || (!bptr_rule)) {
		// Для ускорения сканирования в методе наименьших квадратов 8.07.2017.


		doublereal maximum = -1.0e60;
		doublereal maximum1 = -1.0e60;

		integer** q_hash = NULL;
		q_hash = new integer*[maxnod + 1];
		integer* q_ic = NULL;
		q_ic = new integer[maxnod + 1];
		for (integer j = 0; j <= maxnod; j++) {
			q_hash[j] = new integer[9];
			q_ic[j] = 0;
		}
		for (integer j = 0; j <= maxnod; j++) {
			for (integer i_1 = 0; i_1 < 9; i_1++) {
				q_hash[j][i_1] = NON_EXISTENT_NODE;
			}
		}
		for (integer i_1 = 0; i_1 <= maxelm - 1; i_1++) {
			for (integer j_1 = 0; j_1 <= 7; j_1++) {
				q_hash[nvtx[j_1][i_1]][q_ic[nvtx[j_1][i_1]]] = i_1;
				q_ic[nvtx[j_1][i_1]]++;
			}
		}

		integer* inum_now = new integer[maxnod];

		for (integer i = 0; i < maxnod; i++) {
			temp[i] = 0.0;
			vol[i] = 0.0;
			inum_now[i] = 0;
		}

		// Идея расширения шаблона: мы рассматриваем не только текущий nvtx, но и всех
		// его соседей имеющих с ним хоть одну общую вершину. Этим самым шаблон используемый
		// для реконструкции будет расширен новыми точками.
		// Модификация 30.05.2017.
		for (integer i = 0; i <= maxelm - 1; i++) {
			//if (((10 * i) % maxelm) == 0) printf("complete %lld\n", (100 * i / maxelm));
			for (integer j = 0; j <= 7; j++) {
				inum_now[nvtx[j][i] - 1] += 1;
				//for (integer i_1 = 0; i_1 <= maxelm - 1; i_1++) {
				integer i_1 = NON_EXISTENT_NODE;
				for (integer i_2 = 0; i_2 < 9; i_2++) {
					i_1 = q_hash[nvtx[j][i]][i_2];
					if (i_1 >= 0) {
						for (integer j_1 = 0; j_1 <= 7; j_1++) {
							if (i_1 != i) {
								if (nvtx[j][i] == nvtx[j_1][i_1]) {
									inum_now[nvtx[j][i] - 1] += 1;
								}
							}
						}
					}
				}
			}
		}

		for (integer i = 0; i < maxnod; i++) {
			if (inum_now[i] < 1) {
#if doubleintprecision == 1
				printf("i=%lld maxnod=%lld inum_now[%lld]=%lld\n", i, maxnod, i, inum_now[i]);
#else
				printf("i=%d maxnod=%d inum_now[%d]=%d\n", i, maxnod, i, inum_now[i]);
#endif

				system("PAUSE");
			}
		}

		TOCHKA** pointerlist = new TOCHKA*[maxnod];
		doublereal** rthdsd_Gauss = new doublereal*[maxnod];
		for (integer i = 0; i < maxnod; i++) {
			pointerlist[i] = new TOCHKA[(inum_now[i])];
			rthdsd_Gauss[i] = new doublereal[(inum_now[i])];
		}




		for (integer i = 0; i < maxnod; i++) {
			inum_now[i] = 0;
		}



		// Непосредственно само вычисление.
		for (integer i = 0; i <= maxelm - 1; i++) {
			//if (((10 * i) % maxelm) == 0) printf("complete %lld\n", (100 * i / maxelm));

			if (bptr_rule) {
				if (potent[t.ptr[ENUMERATECONTVOL][i]] > maximum) {
					maximum = potent[t.ptr[ENUMERATECONTVOL][i]];
				}
			}
			else {
				if (potent[i] > maximum) {
					maximum = potent[i];
				}
			}

			// вычисление размеров текущего контрольного объёма:
			doublereal dx = 0.0, dy = 0.0, dz = 0.0;// объём текущего контрольного объёма
			volume3D(i, nvtx, pa, dx, dy, dz);

			for (integer j = 0; j <= 7; j++) {
				vol[nvtx[j][i] - 1] += dx * dy*dz;
				//temp[nvtx[j][i] - 1] += dx * dy*dz*fglobal[0].potent[PAM][t.ptr[ENUMERATECONTVOL][i]];
			}

			TOCHKA p;
			center_cord3D(i, nvtx, pa, p, 100);
			p.x = p.x + min_x;
			p.y = p.y + min_y;
			p.z = p.z + min_z;


			for (integer j = 0; j <= 7; j++) {
				pointerlist[nvtx[j][i] - 1][inum_now[nvtx[j][i] - 1]] = p;
				if (fabs(p.x) < 1.0e-30) {
					printf("problem x=%e\n", p.x);
					system("PAUSE");
				}
				if (fabs(p.y) < 1.0e-30) {
					printf("problem y=%e\n", p.y);
					system("PAUSE");
				}
				if (fabs(p.z) < 1.0e-30) {
					printf("problem z=%e\n", p.z);
					system("PAUSE");
				}
				//fglobal[0].potent[PAM][t.ptr[ENUMERATECONTVOL][i]]
				if (bptr_rule) {
					rthdsd_Gauss[nvtx[j][i] - 1][inum_now[nvtx[j][i] - 1]] = potent[t.ptr[ENUMERATECONTVOL][i]];// potent[i];
				}
				else {
					rthdsd_Gauss[nvtx[j][i] - 1][inum_now[nvtx[j][i] - 1]] = potent[i];
				}
				inum_now[nvtx[j][i] - 1]++;
			}





			for (integer j = 0; j <= 7; j++) {
				//for (integer i_1 = 0; i_1 <= maxelm - 1; i_1++) {
				integer i_1 = NON_EXISTENT_NODE;
				for (integer i_2 = 0; i_2 < 9; i_2++) {
					i_1 = q_hash[nvtx[j][i]][i_2];
					if (i_1 >= 0) {
						for (integer j_1 = 0; j_1 <= 7; j_1++) {
							if (i_1 != i) {
								if (nvtx[j][i] == nvtx[j_1][i_1]) {
									TOCHKA p_1;
									center_cord3D(i_1, nvtx, pa, p_1, 100);
									p_1.x = p_1.x + min_x;
									p_1.y = p_1.y + min_y;
									p_1.z = p_1.z + min_z;

									pointerlist[nvtx[j][i] - 1][inum_now[nvtx[j][i] - 1]] = p_1;
									if (fabs(p_1.x) < 1.0e-30) {
										printf("problem x=%e\n", p_1.x);
										system("PAUSE");
									}
									if (fabs(p_1.y) < 1.0e-30) {
										printf("problem y=%e\n", p_1.y);
										system("PAUSE");
									}
									if (fabs(p_1.z) < 1.0e-30) {
										printf("problem z=%e\n", p_1.z);
										system("PAUSE");
									}
									//fglobal[0].potent[PAM][t.ptr[ENUMERATECONTVOL][i_1]]
									if (bptr_rule) {
										rthdsd_Gauss[nvtx[j][i] - 1][inum_now[nvtx[j][i] - 1]] = potent[t.ptr[ENUMERATECONTVOL][i_1]];// potent[i_1];
									}
									else {
										rthdsd_Gauss[nvtx[j][i] - 1][inum_now[nvtx[j][i] - 1]] = potent[i_1];
									}
									inum_now[nvtx[j][i] - 1]++;
								}
							}
						}
					}
				}
			}


			/*
			for (integer j = 0; j <= 7; j++) {
			vol[nvtx[j][i] - 1] += dx*dy*dz;
			temp[nvtx[j][i] - 1] += dx*dy*dz*potent[i];
			}
			*/
		}


		for (integer j = 0; j <= maxnod; j++) {
			delete[] q_hash[j];
		}
		delete[] q_ic;
		delete[] q_hash;
		q_ic = NULL;
		q_hash = NULL;

		//integer jcontrol = 0;
		for (integer i = 0; i < maxnod; i++) {
			//if (((10 * i) % maxnod) == 0) printf("complete %lld\n", (100 * i / maxnod));
			if (fabs(vol[i]) > eps_mashine) {


				doublereal** Xmatr = new doublereal*[10];
				for (integer j = 0; j <= 9; j++) {
					Xmatr[j] = new doublereal[10];
				}


				doublereal* bmatr = new doublereal[10];
				doublereal* koefmatr = new doublereal[10];

				for (integer j1 = 0; j1 <= 9; j1++) {
					for (integer j2 = 0; j2 <= 9; j2++) {
						Xmatr[j1][j2] = 0.0;
					}
					bmatr[j1] = 0.0;
					koefmatr[j1] = 0.0;
				}
				
				for (integer j = 0; j < inum_now[i]; j++) {

					Xmatr[0][0] += 1.0;
					Xmatr[0][1] += pointerlist[i][j].x;
					Xmatr[0][2] += pointerlist[i][j].y;
					Xmatr[0][3] += pointerlist[i][j].z;
					Xmatr[0][4] += pointerlist[i][j].x*pointerlist[i][j].x;
					Xmatr[0][5] += pointerlist[i][j].y*pointerlist[i][j].y;
					Xmatr[0][6] += pointerlist[i][j].z*pointerlist[i][j].z;
					Xmatr[0][7] += pointerlist[i][j].x*pointerlist[i][j].y;
					Xmatr[0][8] += pointerlist[i][j].x*pointerlist[i][j].z;
					Xmatr[0][9] += pointerlist[i][j].y*pointerlist[i][j].z;

					Xmatr[1][0] += pointerlist[i][j].x;
					Xmatr[1][1] += pointerlist[i][j].x*pointerlist[i][j].x;
					Xmatr[1][2] += pointerlist[i][j].x*pointerlist[i][j].y;
					Xmatr[1][3] += pointerlist[i][j].x*pointerlist[i][j].z;
					Xmatr[1][4] += pointerlist[i][j].x*pointerlist[i][j].x*pointerlist[i][j].x;
					Xmatr[1][5] += pointerlist[i][j].x*pointerlist[i][j].y*pointerlist[i][j].y;
					Xmatr[1][6] += pointerlist[i][j].x*pointerlist[i][j].z*pointerlist[i][j].z;
					Xmatr[1][7] += pointerlist[i][j].x*pointerlist[i][j].x*pointerlist[i][j].y;
					Xmatr[1][8] += pointerlist[i][j].x*pointerlist[i][j].x*pointerlist[i][j].z;
					Xmatr[1][9] += pointerlist[i][j].x*pointerlist[i][j].y*pointerlist[i][j].z;

					Xmatr[2][0] += pointerlist[i][j].y;
					Xmatr[2][1] += pointerlist[i][j].y*pointerlist[i][j].x;
					Xmatr[2][2] += pointerlist[i][j].y*pointerlist[i][j].y;
					Xmatr[2][3] += pointerlist[i][j].y*pointerlist[i][j].z;
					Xmatr[2][4] += pointerlist[i][j].y*pointerlist[i][j].x*pointerlist[i][j].x;
					Xmatr[2][5] += pointerlist[i][j].y*pointerlist[i][j].y*pointerlist[i][j].y;
					Xmatr[2][6] += pointerlist[i][j].y*pointerlist[i][j].z*pointerlist[i][j].z;
					Xmatr[2][7] += pointerlist[i][j].y*pointerlist[i][j].x*pointerlist[i][j].y;
					Xmatr[2][8] += pointerlist[i][j].y*pointerlist[i][j].x*pointerlist[i][j].z;
					Xmatr[2][9] += pointerlist[i][j].y*pointerlist[i][j].y*pointerlist[i][j].z;

					Xmatr[3][0] += pointerlist[i][j].z;
					Xmatr[3][1] += pointerlist[i][j].z*pointerlist[i][j].x;
					Xmatr[3][2] += pointerlist[i][j].z*pointerlist[i][j].y;
					Xmatr[3][3] += pointerlist[i][j].z*pointerlist[i][j].z;
					Xmatr[3][4] += pointerlist[i][j].z*pointerlist[i][j].x*pointerlist[i][j].x;
					Xmatr[3][5] += pointerlist[i][j].z*pointerlist[i][j].y*pointerlist[i][j].y;
					Xmatr[3][6] += pointerlist[i][j].z*pointerlist[i][j].z*pointerlist[i][j].z;
					Xmatr[3][7] += pointerlist[i][j].z*pointerlist[i][j].x*pointerlist[i][j].y;
					Xmatr[3][8] += pointerlist[i][j].z*pointerlist[i][j].x*pointerlist[i][j].z;
					Xmatr[3][9] += pointerlist[i][j].z*pointerlist[i][j].y*pointerlist[i][j].z;

					Xmatr[4][0] += pointerlist[i][j].x*pointerlist[i][j].x;
					Xmatr[4][1] += pointerlist[i][j].x*pointerlist[i][j].x*pointerlist[i][j].x;
					Xmatr[4][2] += pointerlist[i][j].x*pointerlist[i][j].x*pointerlist[i][j].y;
					Xmatr[4][3] += pointerlist[i][j].x*pointerlist[i][j].x*pointerlist[i][j].z;
					Xmatr[4][4] += pointerlist[i][j].x*pointerlist[i][j].x*pointerlist[i][j].x*pointerlist[i][j].x;
					Xmatr[4][5] += pointerlist[i][j].x*pointerlist[i][j].x*pointerlist[i][j].y*pointerlist[i][j].y;
					Xmatr[4][6] += pointerlist[i][j].x*pointerlist[i][j].x*pointerlist[i][j].z*pointerlist[i][j].z;
					Xmatr[4][7] += pointerlist[i][j].x*pointerlist[i][j].x*pointerlist[i][j].x*pointerlist[i][j].y;
					Xmatr[4][8] += pointerlist[i][j].x*pointerlist[i][j].x*pointerlist[i][j].x*pointerlist[i][j].z;
					Xmatr[4][9] += pointerlist[i][j].x*pointerlist[i][j].x*pointerlist[i][j].y*pointerlist[i][j].z;

					Xmatr[5][0] += pointerlist[i][j].y*pointerlist[i][j].y;
					Xmatr[5][1] += pointerlist[i][j].y*pointerlist[i][j].y*pointerlist[i][j].x;
					Xmatr[5][2] += pointerlist[i][j].y*pointerlist[i][j].y*pointerlist[i][j].y;
					Xmatr[5][3] += pointerlist[i][j].y*pointerlist[i][j].y*pointerlist[i][j].z;
					Xmatr[5][4] += pointerlist[i][j].y*pointerlist[i][j].y*pointerlist[i][j].x*pointerlist[i][j].x;
					Xmatr[5][5] += pointerlist[i][j].y*pointerlist[i][j].y*pointerlist[i][j].y*pointerlist[i][j].y;
					Xmatr[5][6] += pointerlist[i][j].y*pointerlist[i][j].y*pointerlist[i][j].z*pointerlist[i][j].z;
					Xmatr[5][7] += pointerlist[i][j].y*pointerlist[i][j].y*pointerlist[i][j].x*pointerlist[i][j].y;
					Xmatr[5][8] += pointerlist[i][j].y*pointerlist[i][j].y*pointerlist[i][j].x*pointerlist[i][j].z;
					Xmatr[5][9] += pointerlist[i][j].y*pointerlist[i][j].y*pointerlist[i][j].y*pointerlist[i][j].z;

					Xmatr[6][0] += pointerlist[i][j].z*pointerlist[i][j].z;
					Xmatr[6][1] += pointerlist[i][j].z*pointerlist[i][j].z*pointerlist[i][j].x;
					Xmatr[6][2] += pointerlist[i][j].z*pointerlist[i][j].z*pointerlist[i][j].y;
					Xmatr[6][3] += pointerlist[i][j].z*pointerlist[i][j].z*pointerlist[i][j].z;
					Xmatr[6][4] += pointerlist[i][j].z*pointerlist[i][j].z*pointerlist[i][j].x*pointerlist[i][j].x;
					Xmatr[6][5] += pointerlist[i][j].z*pointerlist[i][j].z*pointerlist[i][j].y*pointerlist[i][j].y;
					Xmatr[6][6] += pointerlist[i][j].z*pointerlist[i][j].z*pointerlist[i][j].z*pointerlist[i][j].z;
					Xmatr[6][7] += pointerlist[i][j].z*pointerlist[i][j].z*pointerlist[i][j].x*pointerlist[i][j].y;
					Xmatr[6][8] += pointerlist[i][j].z*pointerlist[i][j].z*pointerlist[i][j].x*pointerlist[i][j].z;
					Xmatr[6][9] += pointerlist[i][j].z*pointerlist[i][j].z*pointerlist[i][j].y*pointerlist[i][j].z;

					Xmatr[7][0] += pointerlist[i][j].x*pointerlist[i][j].y;
					Xmatr[7][1] += pointerlist[i][j].x*pointerlist[i][j].y*pointerlist[i][j].x;
					Xmatr[7][2] += pointerlist[i][j].x*pointerlist[i][j].y*pointerlist[i][j].y;
					Xmatr[7][3] += pointerlist[i][j].x*pointerlist[i][j].y*pointerlist[i][j].z;
					Xmatr[7][4] += pointerlist[i][j].x*pointerlist[i][j].y*pointerlist[i][j].x*pointerlist[i][j].x;
					Xmatr[7][5] += pointerlist[i][j].x*pointerlist[i][j].y*pointerlist[i][j].y*pointerlist[i][j].y;
					Xmatr[7][6] += pointerlist[i][j].x*pointerlist[i][j].y*pointerlist[i][j].z*pointerlist[i][j].z;
					Xmatr[7][7] += pointerlist[i][j].x*pointerlist[i][j].y*pointerlist[i][j].x*pointerlist[i][j].y;
					Xmatr[7][8] += pointerlist[i][j].x*pointerlist[i][j].y*pointerlist[i][j].x*pointerlist[i][j].z;
					Xmatr[7][9] += pointerlist[i][j].x*pointerlist[i][j].y*pointerlist[i][j].y*pointerlist[i][j].z;

					Xmatr[8][0] += pointerlist[i][j].x*pointerlist[i][j].z;
					Xmatr[8][1] += pointerlist[i][j].x*pointerlist[i][j].z*pointerlist[i][j].x;
					Xmatr[8][2] += pointerlist[i][j].x*pointerlist[i][j].z*pointerlist[i][j].y;
					Xmatr[8][3] += pointerlist[i][j].x*pointerlist[i][j].z*pointerlist[i][j].z;
					Xmatr[8][4] += pointerlist[i][j].x*pointerlist[i][j].z*pointerlist[i][j].x*pointerlist[i][j].x;
					Xmatr[8][5] += pointerlist[i][j].x*pointerlist[i][j].z*pointerlist[i][j].y*pointerlist[i][j].y;
					Xmatr[8][6] += pointerlist[i][j].x*pointerlist[i][j].z*pointerlist[i][j].z*pointerlist[i][j].z;
					Xmatr[8][7] += pointerlist[i][j].x*pointerlist[i][j].z*pointerlist[i][j].x*pointerlist[i][j].y;
					Xmatr[8][8] += pointerlist[i][j].x*pointerlist[i][j].z*pointerlist[i][j].x*pointerlist[i][j].z;
					Xmatr[8][9] += pointerlist[i][j].x*pointerlist[i][j].z*pointerlist[i][j].y*pointerlist[i][j].z;

					Xmatr[9][0] += pointerlist[i][j].y*pointerlist[i][j].z;
					Xmatr[9][1] += pointerlist[i][j].y*pointerlist[i][j].z*pointerlist[i][j].x;
					Xmatr[9][2] += pointerlist[i][j].y*pointerlist[i][j].z*pointerlist[i][j].y;
					Xmatr[9][3] += pointerlist[i][j].y*pointerlist[i][j].z*pointerlist[i][j].z;
					Xmatr[9][4] += pointerlist[i][j].y*pointerlist[i][j].z*pointerlist[i][j].x*pointerlist[i][j].x;
					Xmatr[9][5] += pointerlist[i][j].y*pointerlist[i][j].z*pointerlist[i][j].y*pointerlist[i][j].y;
					Xmatr[9][6] += pointerlist[i][j].y*pointerlist[i][j].z*pointerlist[i][j].z*pointerlist[i][j].z;
					Xmatr[9][7] += pointerlist[i][j].y*pointerlist[i][j].z*pointerlist[i][j].x*pointerlist[i][j].y;
					Xmatr[9][8] += pointerlist[i][j].y*pointerlist[i][j].z*pointerlist[i][j].x*pointerlist[i][j].z;
					Xmatr[9][9] += pointerlist[i][j].y*pointerlist[i][j].z*pointerlist[i][j].y*pointerlist[i][j].z;

					bmatr[0] += rthdsd_Gauss[i][j];
					bmatr[1] += pointerlist[i][j].x*rthdsd_Gauss[i][j];
					bmatr[2] += pointerlist[i][j].y*rthdsd_Gauss[i][j];
					bmatr[3] += pointerlist[i][j].z*rthdsd_Gauss[i][j];
					bmatr[4] += pointerlist[i][j].x*pointerlist[i][j].x*rthdsd_Gauss[i][j];
					bmatr[5] += pointerlist[i][j].y*pointerlist[i][j].y*rthdsd_Gauss[i][j];
					bmatr[6] += pointerlist[i][j].z*pointerlist[i][j].z*rthdsd_Gauss[i][j];
					bmatr[7] += pointerlist[i][j].x*pointerlist[i][j].y*rthdsd_Gauss[i][j];
					bmatr[8] += pointerlist[i][j].x*pointerlist[i][j].z*rthdsd_Gauss[i][j];
					bmatr[9] += pointerlist[i][j].y*pointerlist[i][j].z*rthdsd_Gauss[i][j];

				}

				if ((fabs(Xmatr[0][0]) < 1.e-30) || (fabs(Xmatr[1][1]) < 1.e-30) || (fabs(Xmatr[2][2]) < 1.e-30) || (fabs(Xmatr[3][3]) < 1.e-30)||
					(fabs(Xmatr[4][4]) < 1.e-30) || (fabs(Xmatr[5][5]) < 1.e-30) || (fabs(Xmatr[6][6]) < 1.e-30) || (fabs(Xmatr[7][7]) < 1.e-30)||
					(fabs(Xmatr[8][8]) < 1.e-30) || (fabs(Xmatr[9][9]) < 1.e-30)) {
#if doubleintprecision == 1
					printf("inum_now[%lld]=%lld\n", i, inum_now[i]);
#else
					printf("inum_now[%d]=%d\n", i, inum_now[i]);
#endif

					system("PAUSE");
				}

				//Xmatr*koefmatr = bmatr;
				/*
				// Метод Гаусса не работает т.к. система линейно зависима.
				if (!my_version_gauss1(Xmatr, 10, bmatr, koefmatr, false, i)) {
				temp[i] = temp[i] / vol[i];
				}
				else {
				temp[i] = koefmatr[0] + koefmatr[1] * (pa[i].x+min_x) + koefmatr[2] * (pa[i].y+min_y) + koefmatr[3] * (pa[i].z+min_z)+
				koefmatr[4] * (pa[i].x+min_x)* (pa[i].x+min_x) + koefmatr[5] * (pa[i].y+min_y) * (pa[i].y+min_y) + koefmatr[6] * (pa[i].z+min_z)* (pa[i].z+min_z)+
				koefmatr[7] * (pa[i].x+min_x)* (pa[i].y+min_y) + koefmatr[8] * (pa[i].x+min_x) * (pa[i].z+min_z) + koefmatr[9] * (pa[i].y+min_y)* (pa[i].z+min_z);
				}
				*/
				for (integer j1 = 0; j1 <= 3; j1++) {
					koefmatr[j1] = 0.0;
				}
				for (integer j1 = 0; j1 <= 250; j1++) {
					doublereal alpha = 0.2;
					doublereal d_0 = koefmatr[0];
					doublereal d_1 = koefmatr[1];
					doublereal d_2 = koefmatr[2];
					doublereal d_3 = koefmatr[3];
					doublereal d_4 = koefmatr[4];
					doublereal d_5 = koefmatr[5];
					doublereal d_6 = koefmatr[6];
					doublereal d_7 = koefmatr[7];
					doublereal d_8 = koefmatr[8];
					doublereal d_9 = koefmatr[9];
					koefmatr[0] = (1.0 - alpha)*d_0 + alpha * ((bmatr[0] - Xmatr[0][1] * koefmatr[1] - Xmatr[0][2] * koefmatr[2] - Xmatr[0][3] * koefmatr[3] - Xmatr[0][4] * koefmatr[4] - Xmatr[0][5] * koefmatr[5] - Xmatr[0][6] * koefmatr[6] - Xmatr[0][7] * koefmatr[7] - Xmatr[0][8] * koefmatr[8] - Xmatr[0][9] * koefmatr[9]) / Xmatr[0][0]);
					koefmatr[1] = (1.0 - alpha)*d_1 + alpha * ((bmatr[1] - Xmatr[1][0] * koefmatr[0] - Xmatr[1][2] * koefmatr[2] - Xmatr[1][3] * koefmatr[3] - Xmatr[1][4] * koefmatr[4] - Xmatr[1][5] * koefmatr[5] - Xmatr[1][6] * koefmatr[6] - Xmatr[1][7] * koefmatr[7] - Xmatr[1][8] * koefmatr[8] - Xmatr[1][9] * koefmatr[9]) / Xmatr[1][1]);
					koefmatr[2] = (1.0 - alpha)*d_2 + alpha * ((bmatr[2] - Xmatr[2][0] * koefmatr[0] - Xmatr[2][1] * koefmatr[1] - Xmatr[2][3] * koefmatr[3] - Xmatr[2][4] * koefmatr[4] - Xmatr[2][5] * koefmatr[5] - Xmatr[2][6] * koefmatr[6] - Xmatr[2][7] * koefmatr[7] - Xmatr[2][8] * koefmatr[8] - Xmatr[2][9] * koefmatr[9]) / Xmatr[2][2]);
					koefmatr[3] = (1.0 - alpha)*d_3 + alpha * ((bmatr[3] - Xmatr[3][0] * koefmatr[0] - Xmatr[3][1] * koefmatr[1] - Xmatr[3][2] * koefmatr[2] - Xmatr[3][4] * koefmatr[4] - Xmatr[3][5] * koefmatr[5] - Xmatr[3][6] * koefmatr[6] - Xmatr[3][7] * koefmatr[7] - Xmatr[3][8] * koefmatr[8] - Xmatr[3][9] * koefmatr[9]) / Xmatr[3][3]);
					koefmatr[4] = (1.0 - alpha)*d_4 + alpha * ((bmatr[4] - Xmatr[4][0] * koefmatr[0] - Xmatr[4][1] * koefmatr[1] - Xmatr[4][2] * koefmatr[2] - Xmatr[4][3] * koefmatr[3] - Xmatr[4][5] * koefmatr[5] - Xmatr[4][6] * koefmatr[6] - Xmatr[4][7] * koefmatr[7] - Xmatr[4][8] * koefmatr[8] - Xmatr[4][9] * koefmatr[9]) / Xmatr[4][4]);
					koefmatr[5] = (1.0 - alpha)*d_5 + alpha * ((bmatr[5] - Xmatr[5][0] * koefmatr[0] - Xmatr[5][1] * koefmatr[1] - Xmatr[5][2] * koefmatr[2] - Xmatr[5][3] * koefmatr[3] - Xmatr[5][4] * koefmatr[4] - Xmatr[5][6] * koefmatr[6] - Xmatr[5][7] * koefmatr[7] - Xmatr[5][8] * koefmatr[8] - Xmatr[5][9] * koefmatr[9]) / Xmatr[5][5]);
					koefmatr[6] = (1.0 - alpha)*d_6 + alpha * ((bmatr[6] - Xmatr[6][0] * koefmatr[0] - Xmatr[6][1] * koefmatr[1] - Xmatr[6][2] * koefmatr[2] - Xmatr[6][3] * koefmatr[3] - Xmatr[6][4] * koefmatr[4] - Xmatr[6][5] * koefmatr[5] - Xmatr[6][7] * koefmatr[7] - Xmatr[6][8] * koefmatr[8] - Xmatr[6][9] * koefmatr[9]) / Xmatr[6][6]);
					koefmatr[7] = (1.0 - alpha)*d_7 + alpha * ((bmatr[7] - Xmatr[7][0] * koefmatr[0] - Xmatr[7][1] * koefmatr[1] - Xmatr[7][2] * koefmatr[2] - Xmatr[7][3] * koefmatr[3] - Xmatr[7][4] * koefmatr[4] - Xmatr[7][5] * koefmatr[5] - Xmatr[7][6] * koefmatr[6] - Xmatr[7][8] * koefmatr[8] - Xmatr[7][9] * koefmatr[9]) / Xmatr[7][7]);
					koefmatr[8] = (1.0 - alpha)*d_8 + alpha * ((bmatr[8] - Xmatr[8][0] * koefmatr[0] - Xmatr[8][1] * koefmatr[1] - Xmatr[8][2] * koefmatr[2] - Xmatr[8][3] * koefmatr[3] - Xmatr[8][4] * koefmatr[4] - Xmatr[8][5] * koefmatr[5] - Xmatr[8][6] * koefmatr[6] - Xmatr[8][7] * koefmatr[7] - Xmatr[8][9] * koefmatr[9]) / Xmatr[8][8]);
					koefmatr[9] = (1.0 - alpha)*d_9 + alpha * ((bmatr[9] - Xmatr[9][0] * koefmatr[0] - Xmatr[9][1] * koefmatr[1] - Xmatr[9][2] * koefmatr[2] - Xmatr[9][3] * koefmatr[3] - Xmatr[9][4] * koefmatr[4] - Xmatr[9][5] * koefmatr[5] - Xmatr[9][6] * koefmatr[6] - Xmatr[9][7] * koefmatr[7] - Xmatr[9][8] * koefmatr[8]) / Xmatr[9][9]);
									
				}
				temp[i] = koefmatr[0] + koefmatr[1] * (pa[i].x + min_x) + koefmatr[2] * (pa[i].y + min_y) + koefmatr[3] * (pa[i].z + min_z)+
					koefmatr[4] * (pa[i].x + min_x)* (pa[i].x + min_x) + koefmatr[5] * (pa[i].y + min_y)* (pa[i].y + min_y) + koefmatr[6] * (pa[i].z + min_z)* (pa[i].z + min_z)+
					koefmatr[7] * (pa[i].x + min_x)* (pa[i].y + min_y) + koefmatr[8] * (pa[i].x + min_x)* (pa[i].z + min_z) + koefmatr[9] * (pa[i].y + min_y)* (pa[i].z + min_z);
				if (temp[i] > maximum1) {
					maximum1 = temp[i];
				}
				//temp[i] = koefmatr[0] + koefmatr[1] * (pa[i].x ) + koefmatr[2] * (pa[i].y ) + koefmatr[3] * (pa[i].z );
				//heat_flux_X[i] = koefmatr[1]+ koefmatr[4] * (pa[i].x + min_x)*2.0+ koefmatr[7] * (pa[i].y + min_y)+ koefmatr[8] * (pa[i].z + min_z);
				//heat_flux_Y[i] = koefmatr[2]+ koefmatr[5] * (pa[i].y + min_y)*2.0+ koefmatr[7] * (pa[i].x + min_x)+ koefmatr[9] * (pa[i].z + min_z);
				//heat_flux_Z[i] = koefmatr[3]+ koefmatr[6] * (pa[i].z + min_z)*2.0+ koefmatr[8] * (pa[i].x + min_x)+ koefmatr[9] * (pa[i].y + min_y);
				gradX[i] = koefmatr[1] + koefmatr[4] * (pa[i].x + min_x)*2.0+ koefmatr[7] * (pa[i].y + min_y)+ koefmatr[8] * (pa[i].z + min_z);
				gradY[i] = koefmatr[2] + koefmatr[5] * (pa[i].y + min_y)*2.0+ koefmatr[7] * (pa[i].x + min_x)+ koefmatr[9] * (pa[i].z + min_z);
				gradZ[i] = koefmatr[3] + koefmatr[6] * (pa[i].z + min_z)*2.0+ koefmatr[8] * (pa[i].x + min_x)+ koefmatr[9] * (pa[i].y + min_y);
				//heat_flux_X[i] = 0.0;
				//heat_flux_Y[i] = 0.0;
				//heat_flux_Z[i] = 0.0;
				// вычисление размеров текущего контрольного объёма:
				/*
				doublereal h_1= 1.0e-4*(max_x-min_x1);
				h_1 = pow(fabs(vol[i]), 0.333);
				//h_1 = 0.5*dx1;
				heat_flux_X[i] = ((koefmatr[0] + koefmatr[1] * (pa[i].x + min_x+h_1) + koefmatr[2] * (pa[i].y + min_y) + koefmatr[3] * (pa[i].z + min_z)
				)-(koefmatr[0] + koefmatr[1] * (pa[i].x + min_x-h_1) + koefmatr[2] * (pa[i].y + min_y) + koefmatr[3] * (pa[i].z + min_z)
				)) / (2 * h_1);
				//h_1 = 1.0e-4*(max_y - min_y1);
				//h_1 = 0.5*dy1;
				heat_flux_Y[i] = ((koefmatr[0] + koefmatr[1] * (pa[i].x + min_x) + koefmatr[2] * (pa[i].y + min_y+h_1) + koefmatr[3] * (pa[i].z + min_z)
				) - (koefmatr[0] + koefmatr[1] * (pa[i].x + min_x) + koefmatr[2] * (pa[i].y + min_y-h_1) + koefmatr[3] * (pa[i].z + min_z)
				)) / (2 * h_1);
				//h_1 = 1.0e-4*(max_z - min_z1);
				//h_1 = 0.5*dz1;
				heat_flux_Z[i] = ((koefmatr[0] + koefmatr[1] * (pa[i].x + min_x) + koefmatr[2] * (pa[i].y + min_y) + koefmatr[3] * (pa[i].z + min_z+h_1)
				) - (koefmatr[0] + koefmatr[1] * (pa[i].x + min_x) + koefmatr[2] * (pa[i].y + min_y) + koefmatr[3] * (pa[i].z + min_z-h_1)
				)) / (2 * h_1);
				*/

				for (integer j = 0; j <= 9; j++) {
					delete[] Xmatr[j];
				}
				delete[] Xmatr;
				delete[] bmatr;
				delete[] koefmatr;

			}
			else {
#if doubleintprecision == 1
				printf("fatal error! ANES_tecplot_export in module constr_struct_alice.cpp. i=%lld maxnod=%lld\n", i, maxnod);
				printf("vol[%lld]==%e\n", i, vol[i]);
#else
				printf("fatal error! ANES_tecplot_export in module constr_struct_alice.cpp. i=%d maxnod=%d\n", i, maxnod);
				printf("vol[%d]==%e\n", i, vol[i]);
#endif

				//getchar();
				system("PAUSE");
				temp[i] = 0.0;
			}



		}

		// При преобразовании сохраняем модуль величины неизменным.
		if (fabs(maximum1) > 1.0e-20) {
			for (integer i = 0; i < maxnod; i++) {
				temp[i] *= (maximum / maximum1);
			}
		}

		delete[] inum_now;
		for (integer i = 0; i < maxnod; i++) {
			delete[] pointerlist[i];
			delete[] rthdsd_Gauss[i];
		}
		delete[] pointerlist;
		delete[] rthdsd_Gauss;



		// запись PAM
		for (integer i = 0; i < maxnod; i++) {
			fprintf(fp_4, "%+.16f ", temp[i]);
			if (i % 10 == 0) fprintf(fp_4, "\n");
		}
		fprintf(fp_4, "\n");
	}
} // SECOND_ORDER_QUADRATIC_RECONSTRUCT

// 05.04.2019
  //fglobal[0].potent[PAM]
  // Сигнатура вызова.
  // SECOND_ORDER_QUADRATIC_RECONSTRUCTA( maxnod, maxelm, pa, nvtx, vol, temp, min_x, min_y, min_z, fglobal[0].potent[PAM], t, eps_mashine,false);
  //
void SECOND_ORDER_QUADRATIC_RECONSTRUCTA(
	integer& maxnod, integer& maxelm, TOCHKA*& pa, integer**& nvtx,
	doublereal*& vol, doublereal*& temp,
	doublereal& min_x, doublereal& min_y, doublereal& min_z,
	doublereal*& potent, TEMPER& t, const doublereal eps_mashine,
	bool bptr_rule) {



	if ((bptr_rule && bSIMPLErun_now_for_temperature) || (!bptr_rule)) {
		// Для ускорения сканирования в методе наименьших квадратов 8.07.2017.


		doublereal maximum = -1.0e60;
		doublereal maximum1 = -1.0e60;

		integer** q_hash = NULL;
		q_hash = new integer * [maxnod + 1];
		integer* q_ic = NULL;
		q_ic = new integer[maxnod + 1];
		for (integer j = 0; j <= maxnod; j++) {
			q_hash[j] = new integer[9];
			q_ic[j] = 0;
		}
		for (integer j = 0; j <= maxnod; j++) {
			for (integer i_1 = 0; i_1 < 9; i_1++) {
				q_hash[j][i_1] = NON_EXISTENT_NODE;
			}
		}
		for (integer i_1 = 0; i_1 <= maxelm - 1; i_1++) {
			for (integer j_1 = 0; j_1 <= 7; j_1++) {
				q_hash[nvtx[j_1][i_1]][q_ic[nvtx[j_1][i_1]]] = i_1;
				q_ic[nvtx[j_1][i_1]]++;
			}
		}

		integer* inum_now = new integer[maxnod];

		for (integer i = 0; i < maxnod; i++) {
			temp[i] = 0.0;
			vol[i] = 0.0;
			inum_now[i] = 0;
		}

		// Идея расширения шаблона: мы рассматриваем не только текущий nvtx, но и всех
		// его соседей имеющих с ним хоть одну общую вершину. Этим самым шаблон используемый
		// для реконструкции будет расширен новыми точками.
		// Модификация 30.05.2017.
		for (integer i = 0; i <= maxelm - 1; i++) {
			//if (((10 * i) % maxelm) == 0) printf("complete %lld\n", (100 * i / maxelm));
			for (integer j = 0; j <= 7; j++) {
				inum_now[nvtx[j][i] - 1] += 1;
				//for (integer i_1 = 0; i_1 <= maxelm - 1; i_1++) {
				integer i_1 = NON_EXISTENT_NODE;
				for (integer i_2 = 0; i_2 < 9; i_2++) {
					i_1 = q_hash[nvtx[j][i]][i_2];
					if (i_1 >= 0) {
						for (integer j_1 = 0; j_1 <= 7; j_1++) {
							if (i_1 != i) {
								if (nvtx[j][i] == nvtx[j_1][i_1]) {
									inum_now[nvtx[j][i] - 1] += 1;
								}
							}
						}
					}
				}
			}
		}

		for (integer i = 0; i < maxnod; i++) {
			if (inum_now[i] < 1) {
#if doubleintprecision == 1
				printf("i=%lld maxnod=%lld inum_now[%lld]=%lld\n", i, maxnod, i, inum_now[i]);
#else
				printf("i=%d maxnod=%d inum_now[%d]=%d\n", i, maxnod, i, inum_now[i]);
#endif

				system("pause");
			}
		}

		TOCHKA** pointerlist = new TOCHKA * [maxnod];
		doublereal** rthdsd_Gauss = new doublereal * [maxnod];
		for (integer i = 0; i < maxnod; i++) {
			pointerlist[i] = new TOCHKA[(inum_now[i])];
			rthdsd_Gauss[i] = new doublereal[(inum_now[i])];
		}




		for (integer i = 0; i < maxnod; i++) {
			inum_now[i] = 0;
		}



		// Непосредственно само вычисление.
		for (integer i = 0; i <= maxelm - 1; i++) {
			//if (((10 * i) % maxelm) == 0) printf("complete %lld\n", (100 * i / maxelm));

			if (bptr_rule) {
				if (potent[t.ptr[ENUMERATECONTVOL][i]] > maximum) {
					maximum = potent[t.ptr[ENUMERATECONTVOL][i]];
				}
			}
			else {
				if (potent[i] > maximum) {
					maximum = potent[i];
				}
			}

			// вычисление размеров текущего контрольного объёма:
			doublereal dx = 0.0, dy = 0.0, dz = 0.0;// объём текущего контрольного объёма
			volume3D(i, nvtx, pa, dx, dy, dz);

			for (integer j = 0; j <= 7; j++) {
				vol[nvtx[j][i] - 1] += dx * dy * dz;
				//temp[nvtx[j][i] - 1] += dx * dy*dz*fglobal[0].potent[PAM][t.ptr[ENUMERATECONTVOL][i]];
			}

			TOCHKA p;
			center_cord3D(i, nvtx, pa, p, 100);
			p.x = p.x + min_x;
			p.y = p.y + min_y;
			p.z = p.z + min_z;


			for (integer j = 0; j <= 7; j++) {
				pointerlist[nvtx[j][i] - 1][inum_now[nvtx[j][i] - 1]] = p;
				if (fabs(p.x) < 1.0e-30) {
					printf("problem x=%e\n", p.x);
					system("PAUSE");
				}
				if (fabs(p.y) < 1.0e-30) {
					printf("problem y=%e\n", p.y);
					system("PAUSE");
				}
				if (fabs(p.z) < 1.0e-30) {
					printf("problem z=%e\n", p.z);
					system("PAUSE");
				}
				//fglobal[0].potent[PAM][t.ptr[ENUMERATECONTVOL][i]]
				if (bptr_rule) {
					rthdsd_Gauss[nvtx[j][i] - 1][inum_now[nvtx[j][i] - 1]] = potent[t.ptr[ENUMERATECONTVOL][i]];// potent[i];
				}
				else {
					rthdsd_Gauss[nvtx[j][i] - 1][inum_now[nvtx[j][i] - 1]] = potent[i];
				}
				inum_now[nvtx[j][i] - 1]++;
			}





			for (integer j = 0; j <= 7; j++) {
				//for (integer i_1 = 0; i_1 <= maxelm - 1; i_1++) {
				integer i_1 = NON_EXISTENT_NODE;
				for (integer i_2 = 0; i_2 < 9; i_2++) {
					i_1 = q_hash[nvtx[j][i]][i_2];
					if (i_1 >= 0) {
						for (integer j_1 = 0; j_1 <= 7; j_1++) {
							if (i_1 != i) {
								if (nvtx[j][i] == nvtx[j_1][i_1]) {
									TOCHKA p_1;
									center_cord3D(i_1, nvtx, pa, p_1, 100);
									p_1.x = p_1.x + min_x;
									p_1.y = p_1.y + min_y;
									p_1.z = p_1.z + min_z;

									pointerlist[nvtx[j][i] - 1][inum_now[nvtx[j][i] - 1]] = p_1;
									if (fabs(p_1.x) < 1.0e-30) {
										printf("problem x=%e\n", p_1.x);
										system("PAUSE");
									}
									if (fabs(p_1.y) < 1.0e-30) {
										printf("problem y=%e\n", p_1.y);
										system("PAUSE");
									}
									if (fabs(p_1.z) < 1.0e-30) {
										printf("problem z=%e\n", p_1.z);
										system("PAUSE");
									}
									//fglobal[0].potent[PAM][t.ptr[ENUMERATECONTVOL][i_1]]
									if (bptr_rule) {
										rthdsd_Gauss[nvtx[j][i] - 1][inum_now[nvtx[j][i] - 1]] = potent[t.ptr[ENUMERATECONTVOL][i_1]];// potent[i_1];
									}
									else {
										rthdsd_Gauss[nvtx[j][i] - 1][inum_now[nvtx[j][i] - 1]] = potent[i_1];
									}
									inum_now[nvtx[j][i] - 1]++;
								}
							}
						}
					}
				}
			}


			/*
			for (integer j = 0; j <= 7; j++) {
			vol[nvtx[j][i] - 1] += dx*dy*dz;
			temp[nvtx[j][i] - 1] += dx*dy*dz*potent[i];
			}
			*/
		}


		for (integer j = 0; j <= maxnod; j++) {
			delete[] q_hash[j];
		}
		delete[] q_ic;
		delete[] q_hash;
		q_ic = NULL;
		q_hash = NULL;

		//integer jcontrol = 0;
		for (integer i = 0; i < maxnod; i++) {
			//if (((10 * i) % maxnod) == 0) printf("complete %lld\n", (100 * i / maxnod));
			if (fabs(vol[i]) > eps_mashine) {


				doublereal** Xmatr = new doublereal * [10];
				for (integer j = 0; j <= 9; j++) {
					Xmatr[j] = new doublereal[10];
				}


				doublereal* bmatr = new doublereal[10];
				doublereal* koefmatr = new doublereal[10];

				for (integer j1 = 0; j1 <= 9; j1++) {
					for (integer j2 = 0; j2 <= 9; j2++) {
						Xmatr[j1][j2] = 0.0;
					}
					bmatr[j1] = 0.0;
					koefmatr[j1] = 0.0;
				}

				for (integer j = 0; j < inum_now[i]; j++) {

					Xmatr[0][0] += 1.0;
					Xmatr[0][1] += pointerlist[i][j].x;
					Xmatr[0][2] += pointerlist[i][j].y;
					Xmatr[0][3] += pointerlist[i][j].z;
					Xmatr[0][4] += pointerlist[i][j].x * pointerlist[i][j].x;
					Xmatr[0][5] += pointerlist[i][j].y * pointerlist[i][j].y;
					Xmatr[0][6] += pointerlist[i][j].z * pointerlist[i][j].z;
					Xmatr[0][7] += pointerlist[i][j].x * pointerlist[i][j].y;
					Xmatr[0][8] += pointerlist[i][j].x * pointerlist[i][j].z;
					Xmatr[0][9] += pointerlist[i][j].y * pointerlist[i][j].z;

					Xmatr[1][0] += pointerlist[i][j].x;
					Xmatr[1][1] += pointerlist[i][j].x * pointerlist[i][j].x;
					Xmatr[1][2] += pointerlist[i][j].x * pointerlist[i][j].y;
					Xmatr[1][3] += pointerlist[i][j].x * pointerlist[i][j].z;
					Xmatr[1][4] += pointerlist[i][j].x * pointerlist[i][j].x * pointerlist[i][j].x;
					Xmatr[1][5] += pointerlist[i][j].x * pointerlist[i][j].y * pointerlist[i][j].y;
					Xmatr[1][6] += pointerlist[i][j].x * pointerlist[i][j].z * pointerlist[i][j].z;
					Xmatr[1][7] += pointerlist[i][j].x * pointerlist[i][j].x * pointerlist[i][j].y;
					Xmatr[1][8] += pointerlist[i][j].x * pointerlist[i][j].x * pointerlist[i][j].z;
					Xmatr[1][9] += pointerlist[i][j].x * pointerlist[i][j].y * pointerlist[i][j].z;

					Xmatr[2][0] += pointerlist[i][j].y;
					Xmatr[2][1] += pointerlist[i][j].y * pointerlist[i][j].x;
					Xmatr[2][2] += pointerlist[i][j].y * pointerlist[i][j].y;
					Xmatr[2][3] += pointerlist[i][j].y * pointerlist[i][j].z;
					Xmatr[2][4] += pointerlist[i][j].y * pointerlist[i][j].x * pointerlist[i][j].x;
					Xmatr[2][5] += pointerlist[i][j].y * pointerlist[i][j].y * pointerlist[i][j].y;
					Xmatr[2][6] += pointerlist[i][j].y * pointerlist[i][j].z * pointerlist[i][j].z;
					Xmatr[2][7] += pointerlist[i][j].y * pointerlist[i][j].x * pointerlist[i][j].y;
					Xmatr[2][8] += pointerlist[i][j].y * pointerlist[i][j].x * pointerlist[i][j].z;
					Xmatr[2][9] += pointerlist[i][j].y * pointerlist[i][j].y * pointerlist[i][j].z;

					Xmatr[3][0] += pointerlist[i][j].z;
					Xmatr[3][1] += pointerlist[i][j].z * pointerlist[i][j].x;
					Xmatr[3][2] += pointerlist[i][j].z * pointerlist[i][j].y;
					Xmatr[3][3] += pointerlist[i][j].z * pointerlist[i][j].z;
					Xmatr[3][4] += pointerlist[i][j].z * pointerlist[i][j].x * pointerlist[i][j].x;
					Xmatr[3][5] += pointerlist[i][j].z * pointerlist[i][j].y * pointerlist[i][j].y;
					Xmatr[3][6] += pointerlist[i][j].z * pointerlist[i][j].z * pointerlist[i][j].z;
					Xmatr[3][7] += pointerlist[i][j].z * pointerlist[i][j].x * pointerlist[i][j].y;
					Xmatr[3][8] += pointerlist[i][j].z * pointerlist[i][j].x * pointerlist[i][j].z;
					Xmatr[3][9] += pointerlist[i][j].z * pointerlist[i][j].y * pointerlist[i][j].z;

					Xmatr[4][0] += pointerlist[i][j].x * pointerlist[i][j].x;
					Xmatr[4][1] += pointerlist[i][j].x * pointerlist[i][j].x * pointerlist[i][j].x;
					Xmatr[4][2] += pointerlist[i][j].x * pointerlist[i][j].x * pointerlist[i][j].y;
					Xmatr[4][3] += pointerlist[i][j].x * pointerlist[i][j].x * pointerlist[i][j].z;
					Xmatr[4][4] += pointerlist[i][j].x * pointerlist[i][j].x * pointerlist[i][j].x * pointerlist[i][j].x;
					Xmatr[4][5] += pointerlist[i][j].x * pointerlist[i][j].x * pointerlist[i][j].y * pointerlist[i][j].y;
					Xmatr[4][6] += pointerlist[i][j].x * pointerlist[i][j].x * pointerlist[i][j].z * pointerlist[i][j].z;
					Xmatr[4][7] += pointerlist[i][j].x * pointerlist[i][j].x * pointerlist[i][j].x * pointerlist[i][j].y;
					Xmatr[4][8] += pointerlist[i][j].x * pointerlist[i][j].x * pointerlist[i][j].x * pointerlist[i][j].z;
					Xmatr[4][9] += pointerlist[i][j].x * pointerlist[i][j].x * pointerlist[i][j].y * pointerlist[i][j].z;

					Xmatr[5][0] += pointerlist[i][j].y * pointerlist[i][j].y;
					Xmatr[5][1] += pointerlist[i][j].y * pointerlist[i][j].y * pointerlist[i][j].x;
					Xmatr[5][2] += pointerlist[i][j].y * pointerlist[i][j].y * pointerlist[i][j].y;
					Xmatr[5][3] += pointerlist[i][j].y * pointerlist[i][j].y * pointerlist[i][j].z;
					Xmatr[5][4] += pointerlist[i][j].y * pointerlist[i][j].y * pointerlist[i][j].x * pointerlist[i][j].x;
					Xmatr[5][5] += pointerlist[i][j].y * pointerlist[i][j].y * pointerlist[i][j].y * pointerlist[i][j].y;
					Xmatr[5][6] += pointerlist[i][j].y * pointerlist[i][j].y * pointerlist[i][j].z * pointerlist[i][j].z;
					Xmatr[5][7] += pointerlist[i][j].y * pointerlist[i][j].y * pointerlist[i][j].x * pointerlist[i][j].y;
					Xmatr[5][8] += pointerlist[i][j].y * pointerlist[i][j].y * pointerlist[i][j].x * pointerlist[i][j].z;
					Xmatr[5][9] += pointerlist[i][j].y * pointerlist[i][j].y * pointerlist[i][j].y * pointerlist[i][j].z;

					Xmatr[6][0] += pointerlist[i][j].z * pointerlist[i][j].z;
					Xmatr[6][1] += pointerlist[i][j].z * pointerlist[i][j].z * pointerlist[i][j].x;
					Xmatr[6][2] += pointerlist[i][j].z * pointerlist[i][j].z * pointerlist[i][j].y;
					Xmatr[6][3] += pointerlist[i][j].z * pointerlist[i][j].z * pointerlist[i][j].z;
					Xmatr[6][4] += pointerlist[i][j].z * pointerlist[i][j].z * pointerlist[i][j].x * pointerlist[i][j].x;
					Xmatr[6][5] += pointerlist[i][j].z * pointerlist[i][j].z * pointerlist[i][j].y * pointerlist[i][j].y;
					Xmatr[6][6] += pointerlist[i][j].z * pointerlist[i][j].z * pointerlist[i][j].z * pointerlist[i][j].z;
					Xmatr[6][7] += pointerlist[i][j].z * pointerlist[i][j].z * pointerlist[i][j].x * pointerlist[i][j].y;
					Xmatr[6][8] += pointerlist[i][j].z * pointerlist[i][j].z * pointerlist[i][j].x * pointerlist[i][j].z;
					Xmatr[6][9] += pointerlist[i][j].z * pointerlist[i][j].z * pointerlist[i][j].y * pointerlist[i][j].z;

					Xmatr[7][0] += pointerlist[i][j].x * pointerlist[i][j].y;
					Xmatr[7][1] += pointerlist[i][j].x * pointerlist[i][j].y * pointerlist[i][j].x;
					Xmatr[7][2] += pointerlist[i][j].x * pointerlist[i][j].y * pointerlist[i][j].y;
					Xmatr[7][3] += pointerlist[i][j].x * pointerlist[i][j].y * pointerlist[i][j].z;
					Xmatr[7][4] += pointerlist[i][j].x * pointerlist[i][j].y * pointerlist[i][j].x * pointerlist[i][j].x;
					Xmatr[7][5] += pointerlist[i][j].x * pointerlist[i][j].y * pointerlist[i][j].y * pointerlist[i][j].y;
					Xmatr[7][6] += pointerlist[i][j].x * pointerlist[i][j].y * pointerlist[i][j].z * pointerlist[i][j].z;
					Xmatr[7][7] += pointerlist[i][j].x * pointerlist[i][j].y * pointerlist[i][j].x * pointerlist[i][j].y;
					Xmatr[7][8] += pointerlist[i][j].x * pointerlist[i][j].y * pointerlist[i][j].x * pointerlist[i][j].z;
					Xmatr[7][9] += pointerlist[i][j].x * pointerlist[i][j].y * pointerlist[i][j].y * pointerlist[i][j].z;

					Xmatr[8][0] += pointerlist[i][j].x * pointerlist[i][j].z;
					Xmatr[8][1] += pointerlist[i][j].x * pointerlist[i][j].z * pointerlist[i][j].x;
					Xmatr[8][2] += pointerlist[i][j].x * pointerlist[i][j].z * pointerlist[i][j].y;
					Xmatr[8][3] += pointerlist[i][j].x * pointerlist[i][j].z * pointerlist[i][j].z;
					Xmatr[8][4] += pointerlist[i][j].x * pointerlist[i][j].z * pointerlist[i][j].x * pointerlist[i][j].x;
					Xmatr[8][5] += pointerlist[i][j].x * pointerlist[i][j].z * pointerlist[i][j].y * pointerlist[i][j].y;
					Xmatr[8][6] += pointerlist[i][j].x * pointerlist[i][j].z * pointerlist[i][j].z * pointerlist[i][j].z;
					Xmatr[8][7] += pointerlist[i][j].x * pointerlist[i][j].z * pointerlist[i][j].x * pointerlist[i][j].y;
					Xmatr[8][8] += pointerlist[i][j].x * pointerlist[i][j].z * pointerlist[i][j].x * pointerlist[i][j].z;
					Xmatr[8][9] += pointerlist[i][j].x * pointerlist[i][j].z * pointerlist[i][j].y * pointerlist[i][j].z;

					Xmatr[9][0] += pointerlist[i][j].y * pointerlist[i][j].z;
					Xmatr[9][1] += pointerlist[i][j].y * pointerlist[i][j].z * pointerlist[i][j].x;
					Xmatr[9][2] += pointerlist[i][j].y * pointerlist[i][j].z * pointerlist[i][j].y;
					Xmatr[9][3] += pointerlist[i][j].y * pointerlist[i][j].z * pointerlist[i][j].z;
					Xmatr[9][4] += pointerlist[i][j].y * pointerlist[i][j].z * pointerlist[i][j].x * pointerlist[i][j].x;
					Xmatr[9][5] += pointerlist[i][j].y * pointerlist[i][j].z * pointerlist[i][j].y * pointerlist[i][j].y;
					Xmatr[9][6] += pointerlist[i][j].y * pointerlist[i][j].z * pointerlist[i][j].z * pointerlist[i][j].z;
					Xmatr[9][7] += pointerlist[i][j].y * pointerlist[i][j].z * pointerlist[i][j].x * pointerlist[i][j].y;
					Xmatr[9][8] += pointerlist[i][j].y * pointerlist[i][j].z * pointerlist[i][j].x * pointerlist[i][j].z;
					Xmatr[9][9] += pointerlist[i][j].y * pointerlist[i][j].z * pointerlist[i][j].y * pointerlist[i][j].z;

					bmatr[0] += rthdsd_Gauss[i][j];
					bmatr[1] += pointerlist[i][j].x * rthdsd_Gauss[i][j];
					bmatr[2] += pointerlist[i][j].y * rthdsd_Gauss[i][j];
					bmatr[3] += pointerlist[i][j].z * rthdsd_Gauss[i][j];
					bmatr[4] += pointerlist[i][j].x * pointerlist[i][j].x * rthdsd_Gauss[i][j];
					bmatr[5] += pointerlist[i][j].y * pointerlist[i][j].y * rthdsd_Gauss[i][j];
					bmatr[6] += pointerlist[i][j].z * pointerlist[i][j].z * rthdsd_Gauss[i][j];
					bmatr[7] += pointerlist[i][j].x * pointerlist[i][j].y * rthdsd_Gauss[i][j];
					bmatr[8] += pointerlist[i][j].x * pointerlist[i][j].z * rthdsd_Gauss[i][j];
					bmatr[9] += pointerlist[i][j].y * pointerlist[i][j].z * rthdsd_Gauss[i][j];

				}

				if ((fabs(Xmatr[0][0]) < 1.e-30) || (fabs(Xmatr[1][1]) < 1.e-30) || (fabs(Xmatr[2][2]) < 1.e-30) || (fabs(Xmatr[3][3]) < 1.e-30) ||
					(fabs(Xmatr[4][4]) < 1.e-30) || (fabs(Xmatr[5][5]) < 1.e-30) || (fabs(Xmatr[6][6]) < 1.e-30) || (fabs(Xmatr[7][7]) < 1.e-30) ||
					(fabs(Xmatr[8][8]) < 1.e-30) || (fabs(Xmatr[9][9]) < 1.e-30)) {
#if doubleintprecision == 1
					printf("inum_now[%lld]=%lld\n", i, inum_now[i]);
#else
					printf("inum_now[%d]=%d\n", i, inum_now[i]);
#endif

					system("PAUSE");
				}

				//Xmatr*koefmatr = bmatr;
				/*
				// Метод Гаусса не работает т.к. система линейно зависима.
				if (!my_version_gauss1(Xmatr, 10, bmatr, koefmatr, false, i)) {
				temp[i] = temp[i] / vol[i];
				}
				else {
				temp[i] = koefmatr[0] + koefmatr[1] * (pa[i].x+min_x) + koefmatr[2] * (pa[i].y+min_y) + koefmatr[3] * (pa[i].z+min_z)+
				koefmatr[4] * (pa[i].x+min_x)* (pa[i].x+min_x) + koefmatr[5] * (pa[i].y+min_y) * (pa[i].y+min_y) + koefmatr[6] * (pa[i].z+min_z)* (pa[i].z+min_z)+
				koefmatr[7] * (pa[i].x+min_x)* (pa[i].y+min_y) + koefmatr[8] * (pa[i].x+min_x) * (pa[i].z+min_z) + koefmatr[9] * (pa[i].y+min_y)* (pa[i].z+min_z);
				}
				*/
				for (integer j1 = 0; j1 <= 3; j1++) {
					koefmatr[j1] = 0.0;
				}
				for (integer j1 = 0; j1 <= 250; j1++) {
					doublereal alpha = 0.2;
					doublereal d_0 = koefmatr[0];
					doublereal d_1 = koefmatr[1];
					doublereal d_2 = koefmatr[2];
					doublereal d_3 = koefmatr[3];
					doublereal d_4 = koefmatr[4];
					doublereal d_5 = koefmatr[5];
					doublereal d_6 = koefmatr[6];
					doublereal d_7 = koefmatr[7];
					doublereal d_8 = koefmatr[8];
					doublereal d_9 = koefmatr[9];
					koefmatr[0] = (1.0 - alpha) * d_0 + alpha * ((bmatr[0] - Xmatr[0][1] * koefmatr[1] - Xmatr[0][2] * koefmatr[2] - Xmatr[0][3] * koefmatr[3] - Xmatr[0][4] * koefmatr[4] - Xmatr[0][5] * koefmatr[5] - Xmatr[0][6] * koefmatr[6] - Xmatr[0][7] * koefmatr[7] - Xmatr[0][8] * koefmatr[8] - Xmatr[0][9] * koefmatr[9]) / Xmatr[0][0]);
					koefmatr[1] = (1.0 - alpha) * d_1 + alpha * ((bmatr[1] - Xmatr[1][0] * koefmatr[0] - Xmatr[1][2] * koefmatr[2] - Xmatr[1][3] * koefmatr[3] - Xmatr[1][4] * koefmatr[4] - Xmatr[1][5] * koefmatr[5] - Xmatr[1][6] * koefmatr[6] - Xmatr[1][7] * koefmatr[7] - Xmatr[1][8] * koefmatr[8] - Xmatr[1][9] * koefmatr[9]) / Xmatr[1][1]);
					koefmatr[2] = (1.0 - alpha) * d_2 + alpha * ((bmatr[2] - Xmatr[2][0] * koefmatr[0] - Xmatr[2][1] * koefmatr[1] - Xmatr[2][3] * koefmatr[3] - Xmatr[2][4] * koefmatr[4] - Xmatr[2][5] * koefmatr[5] - Xmatr[2][6] * koefmatr[6] - Xmatr[2][7] * koefmatr[7] - Xmatr[2][8] * koefmatr[8] - Xmatr[2][9] * koefmatr[9]) / Xmatr[2][2]);
					koefmatr[3] = (1.0 - alpha) * d_3 + alpha * ((bmatr[3] - Xmatr[3][0] * koefmatr[0] - Xmatr[3][1] * koefmatr[1] - Xmatr[3][2] * koefmatr[2] - Xmatr[3][4] * koefmatr[4] - Xmatr[3][5] * koefmatr[5] - Xmatr[3][6] * koefmatr[6] - Xmatr[3][7] * koefmatr[7] - Xmatr[3][8] * koefmatr[8] - Xmatr[3][9] * koefmatr[9]) / Xmatr[3][3]);
					koefmatr[4] = (1.0 - alpha) * d_4 + alpha * ((bmatr[4] - Xmatr[4][0] * koefmatr[0] - Xmatr[4][1] * koefmatr[1] - Xmatr[4][2] * koefmatr[2] - Xmatr[4][3] * koefmatr[3] - Xmatr[4][5] * koefmatr[5] - Xmatr[4][6] * koefmatr[6] - Xmatr[4][7] * koefmatr[7] - Xmatr[4][8] * koefmatr[8] - Xmatr[4][9] * koefmatr[9]) / Xmatr[4][4]);
					koefmatr[5] = (1.0 - alpha) * d_5 + alpha * ((bmatr[5] - Xmatr[5][0] * koefmatr[0] - Xmatr[5][1] * koefmatr[1] - Xmatr[5][2] * koefmatr[2] - Xmatr[5][3] * koefmatr[3] - Xmatr[5][4] * koefmatr[4] - Xmatr[5][6] * koefmatr[6] - Xmatr[5][7] * koefmatr[7] - Xmatr[5][8] * koefmatr[8] - Xmatr[5][9] * koefmatr[9]) / Xmatr[5][5]);
					koefmatr[6] = (1.0 - alpha) * d_6 + alpha * ((bmatr[6] - Xmatr[6][0] * koefmatr[0] - Xmatr[6][1] * koefmatr[1] - Xmatr[6][2] * koefmatr[2] - Xmatr[6][3] * koefmatr[3] - Xmatr[6][4] * koefmatr[4] - Xmatr[6][5] * koefmatr[5] - Xmatr[6][7] * koefmatr[7] - Xmatr[6][8] * koefmatr[8] - Xmatr[6][9] * koefmatr[9]) / Xmatr[6][6]);
					koefmatr[7] = (1.0 - alpha) * d_7 + alpha * ((bmatr[7] - Xmatr[7][0] * koefmatr[0] - Xmatr[7][1] * koefmatr[1] - Xmatr[7][2] * koefmatr[2] - Xmatr[7][3] * koefmatr[3] - Xmatr[7][4] * koefmatr[4] - Xmatr[7][5] * koefmatr[5] - Xmatr[7][6] * koefmatr[6] - Xmatr[7][8] * koefmatr[8] - Xmatr[7][9] * koefmatr[9]) / Xmatr[7][7]);
					koefmatr[8] = (1.0 - alpha) * d_8 + alpha * ((bmatr[8] - Xmatr[8][0] * koefmatr[0] - Xmatr[8][1] * koefmatr[1] - Xmatr[8][2] * koefmatr[2] - Xmatr[8][3] * koefmatr[3] - Xmatr[8][4] * koefmatr[4] - Xmatr[8][5] * koefmatr[5] - Xmatr[8][6] * koefmatr[6] - Xmatr[8][7] * koefmatr[7] - Xmatr[8][9] * koefmatr[9]) / Xmatr[8][8]);
					koefmatr[9] = (1.0 - alpha) * d_9 + alpha * ((bmatr[9] - Xmatr[9][0] * koefmatr[0] - Xmatr[9][1] * koefmatr[1] - Xmatr[9][2] * koefmatr[2] - Xmatr[9][3] * koefmatr[3] - Xmatr[9][4] * koefmatr[4] - Xmatr[9][5] * koefmatr[5] - Xmatr[9][6] * koefmatr[6] - Xmatr[9][7] * koefmatr[7] - Xmatr[9][8] * koefmatr[8]) / Xmatr[9][9]);

				}
				temp[i] = koefmatr[0] + koefmatr[1] * (pa[i].x + min_x) + koefmatr[2] * (pa[i].y + min_y) + koefmatr[3] * (pa[i].z + min_z) +
					koefmatr[4] * (pa[i].x + min_x) * (pa[i].x + min_x) + koefmatr[5] * (pa[i].y + min_y) * (pa[i].y + min_y) + koefmatr[6] * (pa[i].z + min_z) * (pa[i].z + min_z) +
					koefmatr[7] * (pa[i].x + min_x) * (pa[i].y + min_y) + koefmatr[8] * (pa[i].x + min_x) * (pa[i].z + min_z) + koefmatr[9] * (pa[i].y + min_y) * (pa[i].z + min_z);
				if (temp[i] > maximum1) {
					maximum1 = temp[i];
				}
				//temp[i] = koefmatr[0] + koefmatr[1] * (pa[i].x ) + koefmatr[2] * (pa[i].y ) + koefmatr[3] * (pa[i].z );
				//heat_flux_X[i] = koefmatr[1]+ koefmatr[4] * (pa[i].x + min_x)*2.0+ koefmatr[7] * (pa[i].y + min_y)+ koefmatr[8] * (pa[i].z + min_z);
				//heat_flux_Y[i] = koefmatr[2]+ koefmatr[5] * (pa[i].y + min_y)*2.0+ koefmatr[7] * (pa[i].x + min_x)+ koefmatr[9] * (pa[i].z + min_z);
				//heat_flux_Z[i] = koefmatr[3]+ koefmatr[6] * (pa[i].z + min_z)*2.0+ koefmatr[8] * (pa[i].x + min_x)+ koefmatr[9] * (pa[i].y + min_y);
				//gradX[i] = koefmatr[1] + koefmatr[4] * (pa[i].x + min_x) * 2.0 + koefmatr[7] * (pa[i].y + min_y) + koefmatr[8] * (pa[i].z + min_z);
				//gradY[i] = koefmatr[2] + koefmatr[5] * (pa[i].y + min_y) * 2.0 + koefmatr[7] * (pa[i].x + min_x) + koefmatr[9] * (pa[i].z + min_z);
				//gradZ[i] = koefmatr[3] + koefmatr[6] * (pa[i].z + min_z) * 2.0 + koefmatr[8] * (pa[i].x + min_x) + koefmatr[9] * (pa[i].y + min_y);
				//heat_flux_X[i] = 0.0;
				//heat_flux_Y[i] = 0.0;
				//heat_flux_Z[i] = 0.0;
				// вычисление размеров текущего контрольного объёма:
				/*
				doublereal h_1= 1.0e-4*(max_x-min_x1);
				h_1 = pow(fabs(vol[i]), 0.333);
				//h_1 = 0.5*dx1;
				heat_flux_X[i] = ((koefmatr[0] + koefmatr[1] * (pa[i].x + min_x+h_1) + koefmatr[2] * (pa[i].y + min_y) + koefmatr[3] * (pa[i].z + min_z)
				)-(koefmatr[0] + koefmatr[1] * (pa[i].x + min_x-h_1) + koefmatr[2] * (pa[i].y + min_y) + koefmatr[3] * (pa[i].z + min_z)
				)) / (2 * h_1);
				//h_1 = 1.0e-4*(max_y - min_y1);
				//h_1 = 0.5*dy1;
				heat_flux_Y[i] = ((koefmatr[0] + koefmatr[1] * (pa[i].x + min_x) + koefmatr[2] * (pa[i].y + min_y+h_1) + koefmatr[3] * (pa[i].z + min_z)
				) - (koefmatr[0] + koefmatr[1] * (pa[i].x + min_x) + koefmatr[2] * (pa[i].y + min_y-h_1) + koefmatr[3] * (pa[i].z + min_z)
				)) / (2 * h_1);
				//h_1 = 1.0e-4*(max_z - min_z1);
				//h_1 = 0.5*dz1;
				heat_flux_Z[i] = ((koefmatr[0] + koefmatr[1] * (pa[i].x + min_x) + koefmatr[2] * (pa[i].y + min_y) + koefmatr[3] * (pa[i].z + min_z+h_1)
				) - (koefmatr[0] + koefmatr[1] * (pa[i].x + min_x) + koefmatr[2] * (pa[i].y + min_y) + koefmatr[3] * (pa[i].z + min_z-h_1)
				)) / (2 * h_1);
				*/

				for (integer j = 0; j <= 9; j++) {
					delete[] Xmatr[j];
				}
				delete[] Xmatr;
				delete[] bmatr;
				delete[] koefmatr;

			}
			else {
#if doubleintprecision == 1
				printf("fatal error! ANES_tecplot_export in module constr_struct_alice.cpp. i=%lld maxnod=%lld\n", i, maxnod);
				printf("vol[%lld]==%e\n", i, vol[i]);
#else
				printf("fatal error! ANES_tecplot_export in module constr_struct_alice.cpp. i=%d maxnod=%d\n", i, maxnod);
				printf("vol[%d]==%e\n", i, vol[i]);
#endif

				//getchar();
				system("PAUSE");
				temp[i] = 0.0;
			}



		}

		// При преобразовании сохраняем модуль величины неизменным.
		if (fabs(maximum1) > 1.0e-20) {
			for (integer i = 0; i < maxnod; i++) {
				temp[i] *= (maximum / maximum1);
			}
		}

		delete[] inum_now;
		for (integer i = 0; i < maxnod; i++) {
			delete[] pointerlist[i];
			delete[] rthdsd_Gauss[i];
		}
		delete[] pointerlist;
		delete[] rthdsd_Gauss;



		
	}
} // SECOND_ORDER_QUADRATIC_RECONSTRUCTA




// Экспорт в программу техплот температуры.
//С АЛИС сетки.
void ANES_tecplot360_export_temperature(integer maxnod, TOCHKA* pa,
	integer maxelm, integer** nvtx, doublereal* potent, TEMPER &t, FLOW* &fglobal, integer i_754,
	BLOCK* b, integer lb) {

	// 2 ноября 2016. машинное эпсилон.
	//doublereal eps_mashine = 1.0e-44; // float
	doublereal eps_mashine = 1.0e-308; // double

	
	// nvtx && pa сформированы, можно экспортировать в tecplot360
	FILE *fp_4 = NULL;
	errno_t err_4=0;
#ifdef MINGW_COMPILLER
	if (i_754 == 1) {
		fp_4 = fopen64("ALICEFLOW0_07_temp_apparat_hot.PLT", "w");
	}
	else {
		fp_4 = fopen64("ALICEFLOW0_07_temp.PLT", "w");
	}
#else
	if (i_754 == 1) {
		err_4 = fopen_s(&fp_4, "ALICEFLOW0_07_temp_apparat_hot.PLT", "w");
	}
	else {
		err_4 = fopen_s(&fp_4, "ALICEFLOW0_07_temp.PLT", "w");
	}
#endif

	if ((err_4) != 0) {
		printf("Create File temp Error\n");
		//getchar();
		system("pause");

	}
	else {

		integer maxelm_loc = 0;
		for (integer i = 0; i < maxelm; i++) {
			if ((t.whot_is_block[i] < 0) || (t.whot_is_block[i] >= lb)) {
				printf("ERROR in function ANES_tecplot360_export_temperature.\n");
				printf("(t.whot_is_block[i] < 0) || (t.whot_is_block[i] >= lb)\n");
				system("PAUSE");
				exit(1);
			}
			if (b[t.whot_is_block[i]].bvisible) {
				// Учитываем только те ячейки которые видимы, 
				// а именно user указал для них в интерфейсе 
				// bvisible == true.
				maxelm_loc++;
			}
		}


		fprintf(fp_4, "TITLE = \"ALICEFLOW0_24\"\n");
		if (bSIMPLErun_now_for_temperature) {
			// CFD
			fprintf(fp_4, "VARIABLES = x, y, z, Temp, Lam, log10_heat_flux_X, log10_heat_flux_Y, log10_heat_flux_Z, log10_heat_flux_mag, PAM, PRESS, VX, VY, VZ, SPEED, Viscosity_ratio, Distance_Wall, Curl, dVx_dx, dVx_dy, dVx_dz, dVy_dx, dVy_dy, dVy_dz, dVz_dx, dVz_dy, dVz_dz, total_deformation, x_deformation, y_deformation, z_deformation\n");
		}
		else {
			fprintf(fp_4, "VARIABLES = x, y, z, Temp, Lam, log10_heat_flux_X, log10_heat_flux_Y, log10_heat_flux_Z, log10_heat_flux_mag, total_deformation, x_deformation, y_deformation, z_deformation\n");
		}
#if doubleintprecision == 1
		fprintf(fp_4, "ZONE T=\"Rampant\", N=%lld, E=%lld, ET=BRICK, F=FEBLOCK\n\n", maxnod, maxelm_loc);
#else
		fprintf(fp_4, "ZONE T=\"Rampant\", N=%d, E=%d, ET=BRICK, F=FEBLOCK\n\n", maxnod, maxelm_loc);
#endif
		
		// запись x
		for (integer i = 0; i < maxnod; i++) {
			fprintf(fp_4, "%+.16f ", pa[i].x);
			if (i % 10 == 0) fprintf(fp_4, "\n");
		}
		fprintf(fp_4, "\n");
		// запись y
		for (integer i = 0; i < maxnod; i++) {
			fprintf(fp_4, "%+.16f ", pa[i].y);
			if (i % 10 == 0) fprintf(fp_4, "\n");
		}
		fprintf(fp_4, "\n");
		// запись z
		for (integer i = 0; i < maxnod; i++) {
			fprintf(fp_4, "%+.16f ", pa[i].z);
			if (i % 10 == 0) fprintf(fp_4, "\n");
		}
		fprintf(fp_4, "\n");
		doublereal* temp = new doublereal[maxnod];
		doublereal* lam = new doublereal[maxnod];
		doublereal* vol = new doublereal[maxnod];
		doublereal* heat_flux_X = new doublereal[maxnod];
		doublereal* heat_flux_Y = new doublereal[maxnod];
		doublereal* heat_flux_Z = new doublereal[maxnod];
		doublereal* heat_flux_mag = new doublereal[maxnod];

		for (integer i = 0; i < maxnod; i++) {
			temp[i] = 0.0;
			vol[i] = 0.0;
			lam[i] = 0.0;
			heat_flux_X[i] = 0.0;
			heat_flux_Y[i] = 0.0;
			heat_flux_Z[i] = 0.0;
			heat_flux_mag[i] = 0.0;
		}
		/*
		// Здесь мы ловим конкретную исключительную ситуацию.
		bool bfound = false;
		for (integer i = 0; i <= maxelm - 1; i++) {
			for (integer j = 0; j <= 7; j++) {
				if (nvtx[j][i] - 1 == 69462) bfound = true;
				if (nvtx[j][i] - 1 == 69463) bfound = true;
				if (nvtx[j][i] - 1 == 69464) bfound = true;
				if (nvtx[j][i] - 1 == 69477) bfound = true;
				if (nvtx[j][i] - 1 == 69478) bfound = true;
			}
		}
		if (bfound) {
			printf("bfound\n");
		}
		else
		{
			printf("notfound\n");
		}
		*/

		const integer iRECONSTRUCTION_METHOD = 3;
		if (0==iRECONSTRUCTION_METHOD) {
			// Метод нулевого порядка.
			ZERO_ORDER_RECONSTRUCT(maxnod, maxelm, pa, nvtx, vol, temp, eps_mashine,
				potent, NULL, NULL, true);			
		}
		else if (1== iRECONSTRUCTION_METHOD) {

			// Метод линейного порядка.
			doublereal min_x = 1e60;
			doublereal min_y = 1e60;
			doublereal min_z = 1e60;
			doublereal max_x = -1e60;
			doublereal max_y = -1e60;
			doublereal max_z = -1e60;

			for (integer i = 0; i < maxnod; i++) {
				if (pa[i].x < min_x) {
					min_x=pa[i].x;
				}
				if (pa[i].y < min_y) {
					min_y = pa[i].y;
				}
				if (pa[i].z < min_z) {
					min_z = pa[i].z;
				}
				if (pa[i].x > max_x) {
					max_x = pa[i].x;
				}
				if (pa[i].y > max_y) {
					max_y = pa[i].y;
				}
				if (pa[i].z > max_z) {
					max_z = pa[i].z;
				}
			}

			//min_x *= 1.2;
			//min_y *= 1.2;
			//min_z *= 1.2;

			

			min_x = 1.05*fabs(max_x - min_x);
			if (min_x < 1.0e-30) {
				min_x = 1.05*fabs(max_x);
			}
			min_y = 1.05*fabs(max_y - min_y);
			if (min_y < 1.0e-30) {
				min_y = 1.05*fabs(max_y);
			}
			min_z = 1.05*fabs(max_z - min_z);
			if (min_z < 1.0e-30) {
				min_z = 1.05*fabs(max_z);
			}


			/*
			if (min_x < 1.0e-30) {
				printf("error!!! negative min_x MNK!\n");
				printf("min_x=%e max_x=%e\n",min_x,max_x);
			}
			if (min_y < 1.0e-30) {
				printf("error!!! negative min_y MNK!\n");
				printf("min_y=%e max_y=%e\n", min_y, max_y);
			}
			if (min_z < 1.0e-30) {
				printf("error!!! negative min_z MNK!\n");
				printf("min_z=%e max_z=%e\n", min_z, max_z);
			}
			*/

			integer* inum_now = new integer[maxnod];

			for (integer i = 0; i < maxnod; i++) {
				temp[i] = 0.0;
				vol[i] = 0.0;
				inum_now[i] = 0;
			}

			for (integer i = 0; i <= maxelm - 1; i++) {
				for (integer j = 0; j <= 7; j++) {
					inum_now[nvtx[j][i] - 1] += 1;
				}
			}

			for (integer i = 0; i < maxnod; i++) {
				if (inum_now[i] < 1) {
#if doubleintprecision == 1
					printf("i=%lld maxnod=%lld inum_now[%lld]=%lld\n", i, maxnod, i, inum_now[i]);
#else
					printf("i=%d maxnod=%d inum_now[%d]=%d\n", i, maxnod, i, inum_now[i]);
#endif
					
					system("pause");
				}
			}

			TOCHKA** pointerlist = new TOCHKA*[maxnod];
			doublereal** rthdsd_Gauss = new doublereal*[maxnod];
			for (integer i = 0; i < maxnod; i++) {
				pointerlist[i] = new TOCHKA[(inum_now[i])];
				rthdsd_Gauss[i] = new doublereal[(inum_now[i])];
			}

		
			for (integer i = 0; i < maxnod; i++) {
				inum_now[i] = 0;
			}

			for (integer i = 0; i <= maxelm - 1; i++) {
				// вычисление размеров текущего контрольного объёма:
				doublereal dx = 0.0, dy = 0.0, dz = 0.0;// объём текущего контроольного объёма
				volume3D(i, nvtx, pa, dx, dy, dz);

				for (integer j = 0; j <= 7; j++) {
					vol[nvtx[j][i] - 1] += dx*dy*dz;
					temp[nvtx[j][i] - 1] += dx*dy*dz*potent[i];
				}

				TOCHKA p;
				center_cord3D(i, nvtx, pa, p, 100);
				p.x = p.x + min_x;
				p.y = p.y + min_y;
				p.z = p.z + min_z;

				
				for (integer j = 0; j <= 7; j++) {
					pointerlist[nvtx[j][i] - 1][inum_now[nvtx[j][i] - 1]] = p;
					if (fabs(p.x) < 1.0e-30) {
						printf("problem x=%e\n",p.x);
						system("PAUSE");
					}
					if (fabs(p.y) < 1.0e-30) {
						printf("problem y=%e\n",p.y);
						system("PAUSE");
					}
					if (fabs(p.z) < 1.0e-30) {
						printf("problem z=%e\n",p.z);
						system("PAUSE");
					}
					rthdsd_Gauss[nvtx[j][i] - 1][inum_now[nvtx[j][i] - 1]] = potent[i];
					inum_now[nvtx[j][i] - 1]++;
				}


				/*
				for (integer j = 0; j <= 7; j++) {
					vol[nvtx[j][i] - 1] += dx*dy*dz;
					temp[nvtx[j][i] - 1] += dx*dy*dz*potent[i];
				}
				*/
			}
			

			//integer jcontrol = 0;
			for (integer i = 0; i < maxnod; i++) {
				if (fabs(vol[i]) > eps_mashine) {
					

					doublereal** Xmatr = new doublereal*[4];
					for (integer j = 0; j <= 3; j++) {
						Xmatr[j] = new doublereal[4];
					}


					doublereal* bmatr = new doublereal[4];
					doublereal* koefmatr = new doublereal[4];

					for (integer j1 = 0; j1 <= 3; j1++) {
						for (integer j2 = 0; j2 <= 3; j2++) {
							Xmatr[j1][j2] = 0.0;
						}
						bmatr[j1] = 0.0;
						koefmatr[j1] = 0.0;
					}

					

					for (integer j = 0; j < inum_now[i]; j++) {
						
						Xmatr[0][0] += 1.0; 
						Xmatr[0][1] += pointerlist[i][j].x;
						Xmatr[0][2] += pointerlist[i][j].y;
						Xmatr[0][3] += pointerlist[i][j].z;

						Xmatr[1][0] += pointerlist[i][j].x;
						Xmatr[1][1] += pointerlist[i][j].x*pointerlist[i][j].x;
						Xmatr[1][2] += pointerlist[i][j].x*pointerlist[i][j].y;
						Xmatr[1][3] += pointerlist[i][j].x*pointerlist[i][j].z;
						
						Xmatr[2][0] += pointerlist[i][j].y;
						Xmatr[2][1] += pointerlist[i][j].y*pointerlist[i][j].x;
						Xmatr[2][2] += pointerlist[i][j].y*pointerlist[i][j].y;
						Xmatr[2][3] += pointerlist[i][j].y*pointerlist[i][j].z;

						Xmatr[3][0] += pointerlist[i][j].z;
						Xmatr[3][1] += pointerlist[i][j].z*pointerlist[i][j].x;
						Xmatr[3][2] += pointerlist[i][j].z*pointerlist[i][j].y;
						Xmatr[3][3] += pointerlist[i][j].z*pointerlist[i][j].z;

						bmatr[0] += rthdsd_Gauss[i][j];
						bmatr[1] += pointerlist[i][j].x*rthdsd_Gauss[i][j];
						bmatr[2] += pointerlist[i][j].y*rthdsd_Gauss[i][j];
						bmatr[3] += pointerlist[i][j].z*rthdsd_Gauss[i][j];
						
						/*
						Xmatr[0][0] += 1.0;
						Xmatr[0][1] += pointerlist[i][j].x-min_x;
						Xmatr[0][2] += pointerlist[i][j].y-min_y;
						Xmatr[0][3] += pointerlist[i][j].z-min_z;

						Xmatr[1][0] += pointerlist[i][j].x-min_x;
						Xmatr[1][1] += (pointerlist[i][j].x-min_x)*(pointerlist[i][j].x-min_x);
						Xmatr[1][2] += (pointerlist[i][j].x-min_x)*(pointerlist[i][j].y-min_y);
						Xmatr[1][3] += (pointerlist[i][j].x-min_x)*(pointerlist[i][j].z-min_z);

						Xmatr[2][0] += pointerlist[i][j].y-min_y;
						Xmatr[2][1] += (pointerlist[i][j].y-min_y)*(pointerlist[i][j].x-min_x);
						Xmatr[2][2] += (pointerlist[i][j].y-min_y)*(pointerlist[i][j].y-min_y);
						Xmatr[2][3] += (pointerlist[i][j].y-min_y)*(pointerlist[i][j].z-min_z);

						Xmatr[3][0] += (pointerlist[i][j].z-min_z);
						Xmatr[3][1] += (pointerlist[i][j].z-min_z)*(pointerlist[i][j].x-min_x);
						Xmatr[3][2] += (pointerlist[i][j].z-min_z)*(pointerlist[i][j].y-min_y);
						Xmatr[3][3] += (pointerlist[i][j].z-min_z)*(pointerlist[i][j].z-min_z);

						bmatr[0] += rthdsd_Gauss[i][j];
						bmatr[1] += (pointerlist[i][j].x-min_x)*rthdsd_Gauss[i][j];
						bmatr[2] += (pointerlist[i][j].y-min_y)*rthdsd_Gauss[i][j];
						bmatr[3] += (pointerlist[i][j].z-min_z)*rthdsd_Gauss[i][j];
						*/
					}

					if ((fabs(Xmatr[0][0]) < 1.e-30) || (fabs(Xmatr[1][1]) < 1.e-30) || (fabs(Xmatr[2][2]) < 1.e-30) || (fabs(Xmatr[3][3]) < 1.e-30)) {
#if doubleintprecision == 1
						printf("inum_now[%lld]=%lld\n", i, inum_now[i]);
#else
						printf("inum_now[%d]=%d\n", i, inum_now[i]);
#endif
						
						system("pause");
					}
				
					//Xmatr*koefmatr = bmatr;
					/*
					if (!my_version_gauss1(Xmatr, 4, bmatr, koefmatr, false, i)) {
						temp[i] = temp[i] / vol[i];
					}
					else {
						temp[i] = koefmatr[0] + koefmatr[1] * (pa[i].x+min_x) + koefmatr[2] * (pa[i].y+min_y) + koefmatr[3] * (pa[i].z+min_z);
					}
					*/
					for (integer j1 = 0; j1 <= 100; j1++) {
						koefmatr[0] = (bmatr[0] - Xmatr[0][1]*koefmatr[1] - Xmatr[0][2]*koefmatr[2] - Xmatr[0][3]*koefmatr[3])/ Xmatr[0][0];
						koefmatr[1] = (bmatr[1] - Xmatr[1][0] * koefmatr[0] - Xmatr[1][2] * koefmatr[2] - Xmatr[1][3] * koefmatr[3]) / Xmatr[1][1];
						koefmatr[2] = (bmatr[2] - Xmatr[2][0] * koefmatr[0] - Xmatr[2][1] * koefmatr[1] - Xmatr[2][3] * koefmatr[3]) / Xmatr[2][2];
						koefmatr[3] = (bmatr[3] - Xmatr[3][0] * koefmatr[0] - Xmatr[3][1] * koefmatr[1] - Xmatr[3][2] * koefmatr[2]) / Xmatr[3][3];
					}
					temp[i] = koefmatr[0] + koefmatr[1] * (pa[i].x + min_x) + koefmatr[2] * (pa[i].y + min_y) + koefmatr[3] * (pa[i].z + min_z);
					heat_flux_X[i] = koefmatr[1];
					heat_flux_Y[i] = koefmatr[2];
					heat_flux_Z[i] = koefmatr[3];


					for (integer j = 0; j <= 3; j++) {
						delete[] Xmatr[j];
					}
					delete[] Xmatr;
					delete[] bmatr;
					delete[] koefmatr;

				}
				else {
#if doubleintprecision == 1
					printf("fatal error! ANES_tecplot_export in module constr_struct_alice.cpp. i=%lld maxnod=%lld\n", i, maxnod);
					printf("vol[%lld]==%e\n", i, vol[i]);
#else
					printf("fatal error! ANES_tecplot_export in module constr_struct_alice.cpp. i=%d maxnod=%d\n", i, maxnod);
					printf("vol[%d]==%e\n", i, vol[i]);
#endif
					
					//getchar();
					system("PAUSE");
					temp[i] = 0.0;
				}
			}

			

			delete[] inum_now;
			for (integer i = 0; i < maxnod; i++) {
				delete[] pointerlist[i];
				delete[] rthdsd_Gauss[i];
			}
			delete[] pointerlist;
			delete[] rthdsd_Gauss;
		}
		else if (2== iRECONSTRUCTION_METHOD) {

			// Метод второго порядка.
			doublereal min_x = 1e60;
			doublereal min_y = 1e60;
			doublereal min_z = 1e60;
			doublereal max_x = -1e60;
			doublereal max_y = -1e60;
			doublereal max_z = -1e60;

			for (integer i = 0; i < maxnod; i++) {
				if (pa[i].x < min_x) {
					min_x = pa[i].x;
				}
				if (pa[i].y < min_y) {
					min_y = pa[i].y;
				}
				if (pa[i].z < min_z) {
					min_z = pa[i].z;
				}
				if (pa[i].x > max_x) {
					max_x = pa[i].x;
				}
				if (pa[i].y > max_y) {
					max_y = pa[i].y;
				}
				if (pa[i].z > max_z) {
					max_z = pa[i].z;
				}
			}

			//min_x *= 1.2;
			//min_y *= 1.2;
			//min_z *= 1.2;

			

			min_x = 1.05*fabs(max_x - min_x);
			if (min_x < 1.0e-30) {
				min_x = 1.05*fabs(max_x);
			}
			min_y = 1.05*fabs(max_y - min_y);
			if (min_y < 1.0e-30) {
				min_y = 1.05*fabs(max_y);
			}
			min_z = 1.05*fabs(max_z - min_z);
			if (min_z < 1.0e-30) {
				min_z = 1.05*fabs(max_z);
			}


			/*
			if (min_x < 1.0e-30) {
				printf("error!!! negative min_x MNK!\n");
				printf("min_x=%e max_x=%e\n", min_x, max_x);
			}
			if (min_y < 1.0e-30) {
				printf("error!!! negative min_y MNK!\n");
				printf("min_y=%e max_y=%e\n", min_y, max_y);
			}
			if (min_z < 1.0e-30) {
				printf("error!!! negative min_z MNK!\n");
				printf("min_z=%e max_z=%e\n", min_z, max_z);
			}
			*/

			integer* inum_now = new integer[maxnod];

			for (integer i = 0; i < maxnod; i++) {
				temp[i] = 0.0;
				vol[i] = 0.0;
				inum_now[i] = 0;
			}

			for (integer i = 0; i <= maxelm - 1; i++) {
				for (integer j = 0; j <= 7; j++) {
					inum_now[nvtx[j][i] - 1] += 1;
				}
			}

			for (integer i = 0; i < maxnod; i++) {
				if (inum_now[i] < 1) {
#if doubleintprecision == 1
					printf("i=%lld maxnod=%lld inum_now[%lld]=%lld\n", i, maxnod, i, inum_now[i]);
#else
					printf("i=%d maxnod=%d inum_now[%d]=%d\n", i, maxnod, i, inum_now[i]);
#endif
					
					system("pause");
				}
			}

			TOCHKA** pointerlist = new TOCHKA*[maxnod];
			doublereal** rthdsd_Gauss = new doublereal*[maxnod];
			for (integer i = 0; i < maxnod; i++) {
				pointerlist[i] = new TOCHKA[(inum_now[i])];
				rthdsd_Gauss[i] = new doublereal[(inum_now[i])];
			}


			for (integer i = 0; i < maxnod; i++) {
				inum_now[i] = 0;
			}

			for (integer i = 0; i <= maxelm - 1; i++) {
				// вычисление размеров текущего контрольного объёма:
				doublereal dx = 0.0, dy = 0.0, dz = 0.0;// объём текущего контроольного объёма
				volume3D(i, nvtx, pa, dx, dy, dz);

				for (integer j = 0; j <= 7; j++) {
					vol[nvtx[j][i] - 1] += dx*dy*dz;
					temp[nvtx[j][i] - 1] += dx*dy*dz*potent[i];
				}

				TOCHKA p;
				center_cord3D(i, nvtx, pa, p, 100);
				p.x = p.x + min_x;
				p.y = p.y + min_y;
				p.z = p.z + min_z;


				for (integer j = 0; j <= 7; j++) {
					pointerlist[nvtx[j][i] - 1][inum_now[nvtx[j][i] - 1]] = p;
					if (fabs(p.x) < 1.0e-30) {
						printf("problem x=%e\n", p.x);
						system("PAUSE");
					}
					if (fabs(p.y) < 1.0e-30) {
						printf("problem y=%e\n", p.y);
						system("PAUSE");
					}
					if (fabs(p.z) < 1.0e-30) {
						printf("problem z=%e\n", p.z);
						system("PAUSE");
					}
					rthdsd_Gauss[nvtx[j][i] - 1][inum_now[nvtx[j][i] - 1]] = potent[i];
					inum_now[nvtx[j][i] - 1]++;
				}


				/*
				for (integer j = 0; j <= 7; j++) {
				vol[nvtx[j][i] - 1] += dx*dy*dz;
				temp[nvtx[j][i] - 1] += dx*dy*dz*potent[i];
				}
				*/
			}


			//integer jcontrol = 0;
			for (integer i = 0; i < maxnod; i++) {
				if (fabs(vol[i]) > eps_mashine) {


					doublereal** Xmatr = new doublereal*[7];
					for (integer j = 0; j <= 6; j++) {
						Xmatr[j] = new doublereal[7];
					}


					doublereal* bmatr = new doublereal[7];
					doublereal* koefmatr = new doublereal[7];

					for (integer j1 = 0; j1 <= 6; j1++) {
						for (integer j2 = 0; j2 <= 6; j2++) {
							Xmatr[j1][j2] = 0.0;
						}
						bmatr[j1] = 0.0;
						koefmatr[j1] = 0.0;
					}



					for (integer j = 0; j < inum_now[i]; j++) {

						Xmatr[0][0] += 1.0;
						Xmatr[0][1] += pointerlist[i][j].x;
						Xmatr[0][2] += pointerlist[i][j].y;
						Xmatr[0][3] += pointerlist[i][j].z;
						Xmatr[0][4] += pointerlist[i][j].x*pointerlist[i][j].x;
						Xmatr[0][5] += pointerlist[i][j].y*pointerlist[i][j].y;
						Xmatr[0][6] += pointerlist[i][j].z*pointerlist[i][j].z;

						Xmatr[1][0] += pointerlist[i][j].x;
						Xmatr[1][1] += pointerlist[i][j].x*pointerlist[i][j].x;
						Xmatr[1][2] += pointerlist[i][j].x*pointerlist[i][j].y;
						Xmatr[1][3] += pointerlist[i][j].x*pointerlist[i][j].z;
						Xmatr[1][4] += pointerlist[i][j].x*pointerlist[i][j].x*pointerlist[i][j].x;
						Xmatr[1][5] += pointerlist[i][j].x*pointerlist[i][j].y*pointerlist[i][j].y;
						Xmatr[1][6] += pointerlist[i][j].x*pointerlist[i][j].z*pointerlist[i][j].z;

						Xmatr[2][0] += pointerlist[i][j].y;
						Xmatr[2][1] += pointerlist[i][j].y*pointerlist[i][j].x;
						Xmatr[2][2] += pointerlist[i][j].y*pointerlist[i][j].y;
						Xmatr[2][3] += pointerlist[i][j].y*pointerlist[i][j].z;
						Xmatr[2][4] += pointerlist[i][j].y*pointerlist[i][j].x*pointerlist[i][j].x;
						Xmatr[2][5] += pointerlist[i][j].y*pointerlist[i][j].y*pointerlist[i][j].y;
						Xmatr[2][6] += pointerlist[i][j].y*pointerlist[i][j].z*pointerlist[i][j].z;

						Xmatr[3][0] += pointerlist[i][j].z;
						Xmatr[3][1] += pointerlist[i][j].z*pointerlist[i][j].x;
						Xmatr[3][2] += pointerlist[i][j].z*pointerlist[i][j].y;
						Xmatr[3][3] += pointerlist[i][j].z*pointerlist[i][j].z;
						Xmatr[3][4] += pointerlist[i][j].z*pointerlist[i][j].x*pointerlist[i][j].x;
						Xmatr[3][5] += pointerlist[i][j].z*pointerlist[i][j].y*pointerlist[i][j].y;
						Xmatr[3][6] += pointerlist[i][j].z*pointerlist[i][j].z*pointerlist[i][j].z;

						Xmatr[4][0] += pointerlist[i][j].x*pointerlist[i][j].x;
						Xmatr[4][1] += pointerlist[i][j].x*pointerlist[i][j].x*pointerlist[i][j].x;
						Xmatr[4][2] += pointerlist[i][j].x*pointerlist[i][j].x*pointerlist[i][j].y;
						Xmatr[4][3] += pointerlist[i][j].x*pointerlist[i][j].x*pointerlist[i][j].z;
						Xmatr[4][4] += pointerlist[i][j].x*pointerlist[i][j].x*pointerlist[i][j].x*pointerlist[i][j].x;
						Xmatr[4][5] += pointerlist[i][j].x*pointerlist[i][j].x*pointerlist[i][j].y*pointerlist[i][j].y;
						Xmatr[4][6] += pointerlist[i][j].x*pointerlist[i][j].x*pointerlist[i][j].z*pointerlist[i][j].z;

						Xmatr[5][0] += pointerlist[i][j].y*pointerlist[i][j].y;
						Xmatr[5][1] += pointerlist[i][j].y*pointerlist[i][j].y*pointerlist[i][j].x;
						Xmatr[5][2] += pointerlist[i][j].y*pointerlist[i][j].y*pointerlist[i][j].y;
						Xmatr[5][3] += pointerlist[i][j].y*pointerlist[i][j].y*pointerlist[i][j].z;
						Xmatr[5][4] += pointerlist[i][j].y*pointerlist[i][j].y*pointerlist[i][j].x*pointerlist[i][j].x;
						Xmatr[5][5] += pointerlist[i][j].y*pointerlist[i][j].y*pointerlist[i][j].y*pointerlist[i][j].y;
						Xmatr[5][6] += pointerlist[i][j].y*pointerlist[i][j].y*pointerlist[i][j].z*pointerlist[i][j].z;

						Xmatr[6][0] += pointerlist[i][j].z*pointerlist[i][j].z;
						Xmatr[6][1] += pointerlist[i][j].z*pointerlist[i][j].z*pointerlist[i][j].x;
						Xmatr[6][2] += pointerlist[i][j].z*pointerlist[i][j].z*pointerlist[i][j].y;
						Xmatr[6][3] += pointerlist[i][j].z*pointerlist[i][j].z*pointerlist[i][j].z;
						Xmatr[6][4] += pointerlist[i][j].z*pointerlist[i][j].z*pointerlist[i][j].x*pointerlist[i][j].x;
						Xmatr[6][5] += pointerlist[i][j].z*pointerlist[i][j].z*pointerlist[i][j].y*pointerlist[i][j].y;
						Xmatr[6][6] += pointerlist[i][j].z*pointerlist[i][j].z*pointerlist[i][j].z*pointerlist[i][j].z;

						bmatr[0] += rthdsd_Gauss[i][j];
						bmatr[1] += pointerlist[i][j].x*rthdsd_Gauss[i][j];
						bmatr[2] += pointerlist[i][j].y*rthdsd_Gauss[i][j];
						bmatr[3] += pointerlist[i][j].z*rthdsd_Gauss[i][j];
						bmatr[4] += pointerlist[i][j].x*pointerlist[i][j].x*rthdsd_Gauss[i][j];
						bmatr[5] += pointerlist[i][j].y*pointerlist[i][j].y*rthdsd_Gauss[i][j];
						bmatr[6] += pointerlist[i][j].z*pointerlist[i][j].z*rthdsd_Gauss[i][j];

					
					}

					if ((fabs(Xmatr[0][0]) < 1.e-30) || (fabs(Xmatr[1][1]) < 1.e-30) || (fabs(Xmatr[2][2]) < 1.e-30) || (fabs(Xmatr[3][3]) < 1.e-30)
						|| (fabs(Xmatr[4][4]) < 1.e-30) || (fabs(Xmatr[5][5]) < 1.e-30) || (fabs(Xmatr[6][6]) < 1.e-30)) {
#if doubleintprecision == 1
						printf("inum_now[%lld]=%lld\n", i, inum_now[i]);
#else
						printf("inum_now[%d]=%d\n", i, inum_now[i]);
#endif
						
						system("pause");
					}

					//Xmatr*koefmatr = bmatr;
					/*
					if (!my_version_gauss1(Xmatr, 4, bmatr, koefmatr, false, i)) {
					temp[i] = temp[i] / vol[i];
					}
					else {
					temp[i] = koefmatr[0] + koefmatr[1] * (pa[i].x+min_x) + koefmatr[2] * (pa[i].y+min_y) + koefmatr[3] * (pa[i].z+min_z);
					}
					*/
					for (integer j1 = 0; j1 <= 2000; j1++) {
						koefmatr[0] = (bmatr[0] - Xmatr[0][1] * koefmatr[1] - Xmatr[0][2] * koefmatr[2] - Xmatr[0][3] * koefmatr[3] - Xmatr[0][4] * koefmatr[4] - Xmatr[0][5] * koefmatr[5] - Xmatr[0][6] * koefmatr[6]) / Xmatr[0][0];
						koefmatr[1] = (bmatr[1] - Xmatr[1][0] * koefmatr[0] - Xmatr[1][2] * koefmatr[2] - Xmatr[1][3] * koefmatr[3] - Xmatr[1][4] * koefmatr[4] - Xmatr[1][5] * koefmatr[5] - Xmatr[1][6] * koefmatr[6]) / Xmatr[1][1];
						koefmatr[2] = (bmatr[2] - Xmatr[2][0] * koefmatr[0] - Xmatr[2][1] * koefmatr[1] - Xmatr[2][3] * koefmatr[3] - Xmatr[2][4] * koefmatr[4] - Xmatr[2][5] * koefmatr[5] - Xmatr[2][6] * koefmatr[6]) / Xmatr[2][2];
						koefmatr[3] = (bmatr[3] - Xmatr[3][0] * koefmatr[0] - Xmatr[3][1] * koefmatr[1] - Xmatr[3][2] * koefmatr[2] - Xmatr[3][4] * koefmatr[4] - Xmatr[3][5] * koefmatr[5] - Xmatr[3][6] * koefmatr[6]) / Xmatr[3][3];
						koefmatr[4] = (bmatr[4] - Xmatr[4][0] * koefmatr[0] - Xmatr[4][1] * koefmatr[1] - Xmatr[4][2] * koefmatr[2] - Xmatr[4][3] * koefmatr[3] - Xmatr[4][5] * koefmatr[5] - Xmatr[4][6] * koefmatr[6]) / Xmatr[4][4];
						koefmatr[5] = (bmatr[5] - Xmatr[5][0] * koefmatr[0] - Xmatr[5][1] * koefmatr[1] - Xmatr[5][2] * koefmatr[2] - Xmatr[5][4] * koefmatr[4] - Xmatr[5][3] * koefmatr[3] - Xmatr[5][6] * koefmatr[6]) / Xmatr[5][5];
						koefmatr[6] = (bmatr[6] - Xmatr[6][0] * koefmatr[0] - Xmatr[6][1] * koefmatr[1] - Xmatr[6][2] * koefmatr[2] - Xmatr[6][4] * koefmatr[4] - Xmatr[6][5] * koefmatr[5] - Xmatr[6][3] * koefmatr[3]) / Xmatr[6][6];
					}
					temp[i] = koefmatr[0] + koefmatr[1] * (pa[i].x + min_x) + koefmatr[2] * (pa[i].y + min_y) + koefmatr[3] * (pa[i].z + min_z) + koefmatr[4] * (pa[i].x + min_x)*(pa[i].x + min_x) + koefmatr[5] * (pa[i].y + min_y) * (pa[i].y + min_y) + koefmatr[6] * (pa[i].z + min_z) * (pa[i].z + min_z);
					heat_flux_X[i] = koefmatr[1] + 2.0*koefmatr[4] * (pa[i].x + min_x);
					heat_flux_Y[i] = koefmatr[2] + 2.0*koefmatr[5] * (pa[i].y + min_y);
					heat_flux_Z[i] = koefmatr[3] + 2.0*koefmatr[6] * (pa[i].z + min_z);


					for (integer j = 0; j <= 6; j++) {
						delete[] Xmatr[j];
					}
					delete[] Xmatr;
					delete[] bmatr;
					delete[] koefmatr;

				}
				else {
#if doubleintprecision == 1
					printf("fatal error! ANES_tecplot_export in module constr_struct_alice.cpp. i=%lld maxnod=%lld\n", i, maxnod);
					printf("vol[%lld]==%e\n", i, vol[i]);
#else
					printf("fatal error! ANES_tecplot_export in module constr_struct_alice.cpp. i=%d maxnod=%d\n", i, maxnod);
					printf("vol[%d]==%e\n", i, vol[i]);
#endif
					
					//getchar();
					system("PAUSE");
					temp[i] = 0.0;
				}
			}



			delete[] inum_now;
			for (integer i = 0; i < maxnod; i++) {
				delete[] pointerlist[i];
				delete[] rthdsd_Gauss[i];
			}
			delete[] pointerlist;
			delete[] rthdsd_Gauss;
		}
		else if (3 == iRECONSTRUCTION_METHOD) {
			// линейная реконструкция, НО по расширенному шаблону.
			// Расширенный шаблон содержит больше опорных точек интерполляции, 
			// поэтому в теории реконструкция должна быть более точная.

			// Метод линейного порядка.
			doublereal min_x = 1e60;
			doublereal min_y = 1e60;
			doublereal min_z = 1e60;
			doublereal max_x = -1e60;
			doublereal max_y = -1e60;
			doublereal max_z = -1e60;

			for (integer i = 0; i < maxnod; i++) {
				if (pa[i].x < min_x) {
					min_x = pa[i].x;
				}
				if (pa[i].y < min_y) {
					min_y = pa[i].y;
				}
				if (pa[i].z < min_z) {
					min_z = pa[i].z;
				}
				if (pa[i].x > max_x) {
					max_x = pa[i].x;
				}
				if (pa[i].y > max_y) {
					max_y = pa[i].y;
				}
				if (pa[i].z > max_z) {
					max_z = pa[i].z;
				}
			}

			doublereal min_x1 = min_x;
			doublereal min_y1 = min_y;
			doublereal min_z1 = min_z;
			//min_x *= 1.2;
			//min_y *= 1.2;
			//min_z *= 1.2;

			// 05.07.2017

			min_x = 1.05*fabs(max_x - min_x);
			if (min_x < 1.0e-30) {
				min_x = 1.05*fabs(max_x);
			}
			min_y = 1.05*fabs(max_y - min_y);
			if (min_y < 1.0e-30) {
				min_y = 1.05*fabs(max_y);
			}
			min_z = 1.05*fabs(max_z - min_z);
			if (min_z < 1.0e-30) {
				min_z = 1.05*fabs(max_z);
			}

			/*
			if (min_x < 1.0e-30) {
				printf("error!!! negative min_x MNK!\n");
				printf("min_x=%e max_x=%e\n", min_x, max_x);
			}
			if (min_y < 1.0e-30) {
				printf("error!!! negative min_y MNK!\n");
				printf("min_y=%e max_y=%e\n", min_y, max_y);
			}
			if (min_z < 1.0e-30) {
				printf("error!!! negative min_z MNK!\n");
				printf("min_z=%e max_z=%e\n", min_z, max_z);
			}
			*/

			//*******************START***************************
					
			// Refactoring 04.04.2019
			//FIRST_ORDER_LINEAR_RECONSTRUCT(fp_4, maxnod, maxelm, pa, nvtx, vol, temp, min_x, min_y, min_z, potent, t, eps_mashine, false, heat_flux_X, heat_flux_Y, heat_flux_Z);
			SECOND_ORDER_QUADRATIC_RECONSTRUCT(fp_4, maxnod, maxelm, pa, nvtx, vol, temp, min_x, min_y, min_z, potent, t, eps_mashine, false, heat_flux_X, heat_flux_Y, heat_flux_Z);

			doublereal *Tx = NULL;
			doublereal *Ty = NULL;
			doublereal *Tz = NULL;
			Tx = new doublereal[t.maxelm + t.maxbound];
			Ty = new doublereal[t.maxelm + t.maxbound];
			Tz = new doublereal[t.maxelm + t.maxbound];

			// инициализация нулём.
			for (integer i_9 = 0; i_9 < t.maxelm + t.maxbound; i_9++) {
				Tx[i_9] = 0.0;
				Ty[i_9] = 0.0;
				Tz[i_9] = 0.0;
			}

			// нахождение градиентов.
			for (integer i_9 = 0; i_9 < t.maxelm; i_9++) {
				// Только внутренние узлы.
				green_gaussTemperature(i_9, t.potent, t.nvtx, t.pa,
					t.sosedi, t.maxelm, false,
					t.sosedb, Tx, Ty, Tz, t.ilevel_alice);
			}

			for (integer i_9 = 0; i_9 < t.maxelm; i_9++) {
				// Только граничные узлы.
				green_gaussTemperature(i_9, t.potent, t.nvtx, t.pa,
					t.sosedi, t.maxelm, true,
					t.sosedb, Tx, Ty, Tz, t.ilevel_alice);
			}


			doublereal *Txq = NULL;
			doublereal *Tyq = NULL;
			doublereal *Tzq = NULL;
			Txq = new doublereal[t.maxnod];
			Tyq = new doublereal[t.maxnod];
			Tzq = new doublereal[t.maxnod];

			// инициализация нулём.
			for (integer i_9 = 0; i_9 < t.maxnod; i_9++) {
				Txq[i_9] = 0.0;
				Tyq[i_9] = 0.0;
				Tzq[i_9] = 0.0;
			}

			SECOND_ORDER_QUADRATIC_RECONSTRUCTA(maxnod, maxelm, pa, nvtx, vol, Txq, min_x, min_y, min_z, Tx, t, eps_mashine, false);
			SECOND_ORDER_QUADRATIC_RECONSTRUCTA(maxnod, maxelm, pa, nvtx, vol, Tyq, min_x, min_y, min_z, Ty, t, eps_mashine, false);
			SECOND_ORDER_QUADRATIC_RECONSTRUCTA(maxnod, maxelm, pa, nvtx, vol, Tzq, min_x, min_y, min_z, Tz, t, eps_mashine, false);


			//*****************************END***********************

			if (1) {
				// Метод нулевого порядка.
				ZERO_ORDER_RECONSTRUCT(maxnod, maxelm, pa, nvtx, vol, lam, eps_mashine,
					t.prop[LAM], NULL, NULL, true);				
			}
			// запись теплопроводности
			for (integer i = 0; i < maxnod; i++) {
				fprintf(fp_4, "%+.16f ", lam[i]);
				if (i % 10 == 0) fprintf(fp_4, "\n");
			}
			fprintf(fp_4, "\n");

			// Учитываем ортотропность.
			doublereal* lamx = new doublereal[maxnod];
			doublereal* lamy = new doublereal[maxnod];
			doublereal* lamz = new doublereal[maxnod];
			ZERO_ORDER_RECONSTRUCT_ORTHOTROPY(maxnod, maxelm, pa, nvtx, vol, lamx, lamy, lamz, eps_mashine, t.prop[LAM], NULL, NULL, true, t.prop[MULT_LAM_X], t.prop[MULT_LAM_Y], t.prop[MULT_LAM_Z]);

			// запись теплового потока по Х
			for (integer i = 0; i < maxnod; i++) {
				heat_flux_X[i] = -lamx[i]*Txq[i];
				doublereal d_1 = heat_flux_X[i];
				if (d_1 > 2.0) {
					d_1 = log10(d_1);
				}
				else if (d_1 < -2.0) {
					d_1 = -log10(fabs(d_1));
				}
				else {
					d_1 = 0.0;
				}
				fprintf(fp_4, "%+.16f ", d_1);
				if (i % 10 == 0) fprintf(fp_4, "\n");
			}
			fprintf(fp_4, "\n");

			// запись теплового потока по Y
			for (integer i = 0; i < maxnod; i++) {
				heat_flux_Y[i] = -lamy[i]*Tyq[i];
				doublereal d_1 = heat_flux_Y[i];
				if (d_1 > 2.0) {
					d_1 = log10(d_1);
				}
				else if (d_1 < -2.0) {
					d_1 = -log10(fabs(d_1));
				}
				else {
					d_1 = 0.0;
				}
				fprintf(fp_4, "%+.16f ", d_1);
				if (i % 10 == 0) fprintf(fp_4, "\n");
			}
			fprintf(fp_4, "\n");

			// запись теплового потока по Z
			for (integer i = 0; i < maxnod; i++) {
				heat_flux_Z[i] = -lamz[i]*Tzq[i];
				doublereal d_1 = heat_flux_Z[i];
				if (d_1 > 2.0) {
					d_1 = log10(d_1);
				}
				else if (d_1 < -2.0) {
					d_1 = -log10(fabs(d_1));
				}
				else {
					d_1 = 0.0;
				}
				fprintf(fp_4, "%+.16f ", d_1);
				if (i % 10 == 0) fprintf(fp_4, "\n");
			}
			fprintf(fp_4, "\n");

			delete[] lamx;
			delete[] lamy;
			delete[] lamz;

			// запись модуля теплового потока
			for (integer i = 0; i < maxnod; i++) {
				doublereal d_1 = sqrt(heat_flux_X[i] * heat_flux_X[i] + heat_flux_Y[i] * heat_flux_Y[i] + heat_flux_Z[i] * heat_flux_Z[i]);
				if (d_1 > 2.0) {
					d_1 = log10(d_1);
				}
				else {
					d_1 = 0.0;
				}
				fprintf(fp_4, "%+.16f ", d_1);
				if (i % 10 == 0) fprintf(fp_4, "\n");
			}
			fprintf(fp_4, "\n");

			// Освобождение оперативной памяти.
			if (Tx != NULL) {
				delete[] Tx;
			}
			if (Ty != NULL) {
				delete[] Ty;
			}
			if (Tz != NULL) {
				delete[] Tz;
			}

			// Освобождение оперативной памяти.
			if (Txq != NULL) {
				delete[] Txq;
			}
			if (Tyq != NULL) {
				delete[] Tyq;
			}
			if (Tzq != NULL) {
				delete[] Tzq;
			}

			if (bSIMPLErun_now_for_temperature) {
				//***************** WRITE CFD DATA******************
				// Refactoring 04.04.2019
				//FIRST_ORDER_LINEAR_RECONSTRUCT(fp_4, maxnod, maxelm, pa, nvtx, vol, temp, min_x, min_y, min_z, fglobal[0].potent[PAM], t, eps_mashine,true, heat_flux_X, heat_flux_Y, heat_flux_Z);
				SECOND_ORDER_QUADRATIC_RECONSTRUCT(fp_4, maxnod, maxelm, pa, nvtx, vol, temp, min_x, min_y, min_z, fglobal[0].potent[PAM], t, eps_mashine, true, heat_flux_X, heat_flux_Y, heat_flux_Z);
				//FIRST_ORDER_LINEAR_RECONSTRUCT(fp_4, maxnod, maxelm, pa, nvtx, vol, temp, min_x, min_y, min_z, fglobal[0].potent[PRESS], t, eps_mashine, true, heat_flux_X, heat_flux_Y, heat_flux_Z);
				SECOND_ORDER_QUADRATIC_RECONSTRUCT(fp_4, maxnod, maxelm, pa, nvtx, vol, temp, min_x, min_y, min_z, fglobal[0].potent[PRESS], t, eps_mashine, true, heat_flux_X, heat_flux_Y, heat_flux_Z);
				//FIRST_ORDER_LINEAR_RECONSTRUCT(fp_4, maxnod, maxelm, pa, nvtx, vol, temp, min_x, min_y, min_z, fglobal[0].potent[VX], t, eps_mashine, true, heat_flux_X, heat_flux_Y, heat_flux_Z);
				SECOND_ORDER_QUADRATIC_RECONSTRUCT(fp_4, maxnod, maxelm, pa, nvtx, vol, temp, min_x, min_y, min_z, fglobal[0].potent[VX], t, eps_mashine, true, heat_flux_X, heat_flux_Y, heat_flux_Z);
				//FIRST_ORDER_LINEAR_RECONSTRUCT(fp_4, maxnod, maxelm, pa, nvtx, vol, temp, min_x, min_y, min_z, fglobal[0].potent[VY], t, eps_mashine, true, heat_flux_X, heat_flux_Y, heat_flux_Z);
				SECOND_ORDER_QUADRATIC_RECONSTRUCT(fp_4, maxnod, maxelm, pa, nvtx, vol, temp, min_x, min_y, min_z, fglobal[0].potent[VY], t, eps_mashine, true, heat_flux_X, heat_flux_Y, heat_flux_Z);
				//FIRST_ORDER_LINEAR_RECONSTRUCT(fp_4, maxnod, maxelm, pa, nvtx, vol, temp, min_x, min_y, min_z, fglobal[0].potent[VZ], t, eps_mashine, true, heat_flux_X, heat_flux_Y, heat_flux_Z);
				SECOND_ORDER_QUADRATIC_RECONSTRUCT(fp_4, maxnod, maxelm, pa, nvtx, vol, temp, min_x, min_y, min_z, fglobal[0].potent[VZ], t, eps_mashine, true, heat_flux_X, heat_flux_Y, heat_flux_Z);
				doublereal* tmp_speed = new doublereal[fglobal[0].maxelm + fglobal[0].maxbound];
				if (bSIMPLErun_now_for_temperature) {
					for (integer i_1 = 0; i_1 < fglobal[0].maxelm + fglobal[0].maxbound; i_1++) {
						tmp_speed[i_1] = sqrt(fglobal[0].potent[VX][i_1] * fglobal[0].potent[VX][i_1] + fglobal[0].potent[VY][i_1] * fglobal[0].potent[VY][i_1] + fglobal[0].potent[VZ][i_1] * fglobal[0].potent[VZ][i_1]);
					}
				}
				//FIRST_ORDER_LINEAR_RECONSTRUCT(fp_4, maxnod, maxelm, pa, nvtx, vol, temp, min_x, min_y, min_z, tmp_speed, t, eps_mashine, true, heat_flux_X, heat_flux_Y, heat_flux_Z);
				SECOND_ORDER_QUADRATIC_RECONSTRUCT(fp_4, maxnod, maxelm, pa, nvtx, vol, temp, min_x, min_y, min_z, tmp_speed, t, eps_mashine, true, heat_flux_X, heat_flux_Y, heat_flux_Z);
				for (integer i_1 = 0; i_1 < fglobal[0].maxelm + fglobal[0].maxbound; i_1++) {
					tmp_speed[i_1] = 0.0;
				}
				for (integer i_1 = 0; i_1 < fglobal[0].maxelm; i_1++) {
					tmp_speed[i_1] = fglobal[0].potent[MUT][i_1] / fglobal[0].prop[MU][i_1];
				}
				for (integer i_1 = 0; i_1 < fglobal[0].maxbound; i_1++) {
					tmp_speed[fglobal[0].maxelm + i_1] = fglobal[0].potent[MUT][fglobal[0].maxelm + i_1] / fglobal[0].prop_b[MU][i_1];
				}
				//FIRST_ORDER_LINEAR_RECONSTRUCT(fp_4, maxnod, maxelm, pa, nvtx, vol, temp, min_x, min_y, min_z, tmp_speed, t, eps_mashine, true, heat_flux_X, heat_flux_Y, heat_flux_Z);
				SECOND_ORDER_QUADRATIC_RECONSTRUCT(fp_4, maxnod, maxelm, pa, nvtx, vol, temp, min_x, min_y, min_z, tmp_speed, t, eps_mashine, true, heat_flux_X, heat_flux_Y, heat_flux_Z);
				delete[] tmp_speed;
				if (fglobal[0].rdistWall != NULL) {
					// rdistWall
					//FIRST_ORDER_LINEAR_RECONSTRUCT(fp_4, maxnod, maxelm, pa, nvtx, vol, temp, min_x, min_y, min_z, fglobal[0].rdistWall, t, eps_mashine, true, heat_flux_X, heat_flux_Y, heat_flux_Z);
					SECOND_ORDER_QUADRATIC_RECONSTRUCT(fp_4, maxnod, maxelm, pa, nvtx, vol, temp, min_x, min_y, min_z, fglobal[0].rdistWall, t, eps_mashine, true, heat_flux_X, heat_flux_Y, heat_flux_Z);
				}
				else {
					// В ламинарном режиме расстояние до стенки не рассчитывается.
					for (integer i = 0; i < maxnod; i++) {
						fprintf(fp_4, "%+.16f ", 0.0);
						if (i % 10 == 0) fprintf(fp_4, "\n");
					}
					fprintf(fp_4, "\n");
				}
				//FIRST_ORDER_LINEAR_RECONSTRUCT(fp_4, maxnod, maxelm, pa, nvtx, vol, temp, min_x, min_y, min_z, fglobal[0].potent[CURL], t, eps_mashine, true, heat_flux_X, heat_flux_Y, heat_flux_Z);
				SECOND_ORDER_QUADRATIC_RECONSTRUCT(fp_4, maxnod, maxelm, pa, nvtx, vol, temp, min_x, min_y, min_z, fglobal[0].potent[CURL], t, eps_mashine, true, heat_flux_X, heat_flux_Y, heat_flux_Z);
				//FIRST_ORDER_LINEAR_RECONSTRUCT(fp_4, maxnod, maxelm, pa, nvtx, vol, temp, min_x, min_y, min_z, fglobal[0].potent[GRADXVX], t, eps_mashine, true, heat_flux_X, heat_flux_Y, heat_flux_Z);
				SECOND_ORDER_QUADRATIC_RECONSTRUCT(fp_4, maxnod, maxelm, pa, nvtx, vol, temp, min_x, min_y, min_z, fglobal[0].potent[GRADXVX], t, eps_mashine, true, heat_flux_X, heat_flux_Y, heat_flux_Z);
				//FIRST_ORDER_LINEAR_RECONSTRUCT(fp_4, maxnod, maxelm, pa, nvtx, vol, temp, min_x, min_y, min_z, fglobal[0].potent[GRADYVX], t, eps_mashine, true, heat_flux_X, heat_flux_Y, heat_flux_Z);
				SECOND_ORDER_QUADRATIC_RECONSTRUCT(fp_4, maxnod, maxelm, pa, nvtx, vol, temp, min_x, min_y, min_z, fglobal[0].potent[GRADYVX], t, eps_mashine, true, heat_flux_X, heat_flux_Y, heat_flux_Z);
				//FIRST_ORDER_LINEAR_RECONSTRUCT(fp_4, maxnod, maxelm, pa, nvtx, vol, temp, min_x, min_y, min_z, fglobal[0].potent[GRADZVX], t, eps_mashine, true, heat_flux_X, heat_flux_Y, heat_flux_Z);
				SECOND_ORDER_QUADRATIC_RECONSTRUCT(fp_4, maxnod, maxelm, pa, nvtx, vol, temp, min_x, min_y, min_z, fglobal[0].potent[GRADZVX], t, eps_mashine, true, heat_flux_X, heat_flux_Y, heat_flux_Z);
				//FIRST_ORDER_LINEAR_RECONSTRUCT(fp_4, maxnod, maxelm, pa, nvtx, vol, temp, min_x, min_y, min_z, fglobal[0].potent[GRADXVY], t, eps_mashine, true, heat_flux_X, heat_flux_Y, heat_flux_Z);
				SECOND_ORDER_QUADRATIC_RECONSTRUCT(fp_4, maxnod, maxelm, pa, nvtx, vol, temp, min_x, min_y, min_z, fglobal[0].potent[GRADXVY], t, eps_mashine, true, heat_flux_X, heat_flux_Y, heat_flux_Z);
				//FIRST_ORDER_LINEAR_RECONSTRUCT(fp_4, maxnod, maxelm, pa, nvtx, vol, temp, min_x, min_y, min_z, fglobal[0].potent[GRADYVY], t, eps_mashine, true, heat_flux_X, heat_flux_Y, heat_flux_Z);
				SECOND_ORDER_QUADRATIC_RECONSTRUCT(fp_4, maxnod, maxelm, pa, nvtx, vol, temp, min_x, min_y, min_z, fglobal[0].potent[GRADYVY], t, eps_mashine, true, heat_flux_X, heat_flux_Y, heat_flux_Z);
				//FIRST_ORDER_LINEAR_RECONSTRUCT(fp_4, maxnod, maxelm, pa, nvtx, vol, temp, min_x, min_y, min_z, fglobal[0].potent[GRADZVY], t, eps_mashine, true, heat_flux_X, heat_flux_Y, heat_flux_Z);
				SECOND_ORDER_QUADRATIC_RECONSTRUCT(fp_4, maxnod, maxelm, pa, nvtx, vol, temp, min_x, min_y, min_z, fglobal[0].potent[GRADZVY], t, eps_mashine, true, heat_flux_X, heat_flux_Y, heat_flux_Z);
				//FIRST_ORDER_LINEAR_RECONSTRUCT(fp_4, maxnod, maxelm, pa, nvtx, vol, temp, min_x, min_y, min_z, fglobal[0].potent[GRADXVZ], t, eps_mashine, true, heat_flux_X, heat_flux_Y, heat_flux_Z);
				SECOND_ORDER_QUADRATIC_RECONSTRUCT(fp_4, maxnod, maxelm, pa, nvtx, vol, temp, min_x, min_y, min_z, fglobal[0].potent[GRADXVZ], t, eps_mashine, true, heat_flux_X, heat_flux_Y, heat_flux_Z);
				//FIRST_ORDER_LINEAR_RECONSTRUCT(fp_4, maxnod, maxelm, pa, nvtx, vol, temp, min_x, min_y, min_z, fglobal[0].potent[GRADYVZ], t, eps_mashine, true, heat_flux_X, heat_flux_Y, heat_flux_Z);
				SECOND_ORDER_QUADRATIC_RECONSTRUCT(fp_4, maxnod, maxelm, pa, nvtx, vol, temp, min_x, min_y, min_z, fglobal[0].potent[GRADYVZ], t, eps_mashine, true, heat_flux_X, heat_flux_Y, heat_flux_Z);
				//FIRST_ORDER_LINEAR_RECONSTRUCT(fp_4, maxnod, maxelm, pa, nvtx, vol, temp, min_x, min_y, min_z, fglobal[0].potent[GRADZVZ], t, eps_mashine, true, heat_flux_X, heat_flux_Y, heat_flux_Z);
				SECOND_ORDER_QUADRATIC_RECONSTRUCT(fp_4, maxnod, maxelm, pa, nvtx, vol, temp, min_x, min_y, min_z, fglobal[0].potent[GRADZVZ], t, eps_mashine, true, heat_flux_X, heat_flux_Y, heat_flux_Z);
				//**************END WRITE CFD DATA************************************
			}
		}

		// TODO 5.04.2019 повысить точность аппроксимации деформации в узлах.

			// TOTAL DEFORMATION
			if (1) {
				// Метод нулевого порядка.
				ZERO_ORDER_RECONSTRUCT(maxnod, maxelm, pa, nvtx, vol, lam, eps_mashine, 
					t.total_deformation[XDEFORMATION], t.total_deformation[YDEFORMATION], t.total_deformation[ZDEFORMATION], false);				
			}
			// запись TOTAL DEFORMATION
			for (integer i = 0; i < maxnod; i++) {
				fprintf(fp_4, "%+.16f ", lam[i]);
				if (i % 10 == 0) fprintf(fp_4, "\n");
			}
			fprintf(fp_4, "\n");

			// Deformation X
			if (1) {
				// Метод нулевого порядка.
				ZERO_ORDER_RECONSTRUCT(maxnod, maxelm, pa, nvtx, vol, lam, eps_mashine,
					t.total_deformation[XDEFORMATION], NULL, NULL, true);				
			}
			// запись XDEFORMATION
			for (integer i = 0; i < maxnod; i++) {
				fprintf(fp_4, "%+.16f ", lam[i]);
				if (i % 10 == 0) fprintf(fp_4, "\n");
			}
			fprintf(fp_4, "\n");

			// Deformation Y
			if (1) {
				// Метод нулевого порядка.
				ZERO_ORDER_RECONSTRUCT(maxnod, maxelm, pa, nvtx, vol, lam, eps_mashine,
					t.total_deformation[YDEFORMATION], NULL, NULL, true);				
			}
			// запись YDEFORMATION
			for (integer i = 0; i < maxnod; i++) {
				fprintf(fp_4, "%+.16f ", lam[i]);
				if (i % 10 == 0) fprintf(fp_4, "\n");
			}
			fprintf(fp_4, "\n");

			// Deformation Z
			if (1) {
				// Метод нулевого порядка.
				ZERO_ORDER_RECONSTRUCT(maxnod, maxelm, pa, nvtx, vol, lam, eps_mashine,
					t.total_deformation[ZDEFORMATION], NULL, NULL, true);				
			}
			// запись ZDeformation.
			for (integer i = 0; i < maxnod; i++) {
				fprintf(fp_4, "%+.16f ", lam[i]);
				if (i % 10 == 0) fprintf(fp_4, "\n");
			}
			fprintf(fp_4, "\n");

			delete[] temp;
			delete[] vol;
			delete[] lam;
			delete[] heat_flux_X;
			delete[] heat_flux_Y;
			delete[] heat_flux_Z;
			delete[] heat_flux_mag;
		

		for (integer i = 0; i <= maxelm - 1; i++) {


			//integer ib = t.whot_is_block[i];
			//if (b[ib].bvisible) 
			{


#if doubleintprecision == 1			

				if (1) {
					// Контрольные объёмы состоят из кубиков. Вершины 
					// каждого такого кубика должны быть перечислены в строго порядке сортировки.
					// Сортировка : упорядочивание восьмерки по возракстанию координаты z.
					// Сортировка: упорядочивание каждой четверки по возрастанию координаты y. Всего 2 четверки.
					// Сортировка четырехз двоек по возрастанию х.
					if ((t.whot_is_block[i] < 0) || (t.whot_is_block[i] >= lb)) {
						printf("ERROR in function ANES_tecplot360_export_temperature.\n");
						printf("(t.whot_is_block[i] < 0) || (t.whot_is_block[i] >= lb)\n");
						system("PAUSE");
						exit(1);
					}
					if (b[t.whot_is_block[i]].bvisible) {
						// Учитываем только те ячейки которые видимы, 
						// а именно user указал для них в интерфейсе 
						// bvisible == true.


						integer invtx[8];
						doublereal znvtx[8];
						doublereal ynvtx[8];
						doublereal xnvtx[8];
						for (integer j28 = 0; j28 < 8; j28++) {
							invtx[j28] = nvtx[j28][i];
							znvtx[j28] = pa[nvtx[j28][i] - 1].z;
							ynvtx[j28] = pa[nvtx[j28][i] - 1].y;
							xnvtx[j28] = pa[nvtx[j28][i] - 1].x;
						}

						// Сортировка по возрастанию z.
						for (integer i28 = 1; i28 < 8; i28++) {
							integer inewelm = invtx[i28];
							doublereal znewelm = znvtx[i28];
							doublereal ynewelm = ynvtx[i28];
							doublereal xnewelm = xnvtx[i28];
							integer location = i28 - 1;
							while ((location >= 0) && (znvtx[location] > znewelm)) {
								// сдвигаем все элементы большие очередного.
								invtx[location + 1] = invtx[location];
								znvtx[location + 1] = znvtx[location];
								ynvtx[location + 1] = ynvtx[location];
								xnvtx[location + 1] = xnvtx[location];
								location--;
							}
							invtx[location + 1] = inewelm;
							znvtx[location + 1] = znewelm;
							ynvtx[location + 1] = ynewelm;
							xnvtx[location + 1] = xnewelm;
						}
						// Сортировка по возрастанию y. Часть 1.
						//0,1,2,3<4
						for (integer i28 = 1; i28 < 4; i28++) {
							integer inewelm = invtx[i28];
							doublereal znewelm = znvtx[i28];
							doublereal ynewelm = ynvtx[i28];
							doublereal xnewelm = xnvtx[i28];
							integer location = i28 - 1;
							while ((location >= 0) && (ynvtx[location] > ynewelm)) {
								// сдвигаем все элементы большие очередного.
								invtx[location + 1] = invtx[location];
								znvtx[location + 1] = znvtx[location];
								ynvtx[location + 1] = ynvtx[location];
								xnvtx[location + 1] = xnvtx[location];
								location--;
							}
							invtx[location + 1] = inewelm;
							znvtx[location + 1] = znewelm;
							ynvtx[location + 1] = ynewelm;
							xnvtx[location + 1] = xnewelm;
						}
						// 4,5,6,7<8
						// Сортировка по возрастанию y. Часть 2.
						for (integer i28 = 5; i28 < 8; i28++) {
							integer inewelm = invtx[i28];
							doublereal znewelm = znvtx[i28];
							doublereal ynewelm = ynvtx[i28];
							doublereal xnewelm = xnvtx[i28];
							integer location = i28 - 1;
							while ((location >= 4) && (ynvtx[location] > ynewelm)) {
								// сдвигаем все элементы большие очередного.
								invtx[location + 1] = invtx[location];
								znvtx[location + 1] = znvtx[location];
								ynvtx[location + 1] = ynvtx[location];
								xnvtx[location + 1] = xnvtx[location];
								location--;
							}
							invtx[location + 1] = inewelm;
							znvtx[location + 1] = znewelm;
							ynvtx[location + 1] = ynewelm;
							xnvtx[location + 1] = xnewelm;
						}
						//0,1<2
						//2,3<4
						//4,5<6
						//6,7<8
						// Сортировка по возрастанию X. Часть 1.
						for (integer i28 = 1; i28 < 2; i28++) {
							integer inewelm = invtx[i28];
							doublereal znewelm = znvtx[i28];
							doublereal ynewelm = ynvtx[i28];
							doublereal xnewelm = xnvtx[i28];
							integer location = i28 - 1;
							while ((location >= 0) && (xnvtx[location] > xnewelm)) {
								// сдвигаем все элементы большие очередного.
								invtx[location + 1] = invtx[location];
								znvtx[location + 1] = znvtx[location];
								ynvtx[location + 1] = ynvtx[location];
								xnvtx[location + 1] = xnvtx[location];
								location--;
							}
							invtx[location + 1] = inewelm;
							znvtx[location + 1] = znewelm;
							ynvtx[location + 1] = ynewelm;
							xnvtx[location + 1] = xnewelm;
						}
						// Сортировка по возрастанию X. Часть 2.
						for (integer i28 = 3; i28 < 4; i28++) {
							integer inewelm = invtx[i28];
							doublereal znewelm = znvtx[i28];
							doublereal ynewelm = ynvtx[i28];
							doublereal xnewelm = xnvtx[i28];
							integer location = i28 - 1;
							while ((location >= 2) && (xnvtx[location] < xnewelm)) {
								// сдвигаем все элементы большие очередного.
								invtx[location + 1] = invtx[location];
								znvtx[location + 1] = znvtx[location];
								ynvtx[location + 1] = ynvtx[location];
								xnvtx[location + 1] = xnvtx[location];
								location--;
							}
							invtx[location + 1] = inewelm;
							znvtx[location + 1] = znewelm;
							ynvtx[location + 1] = ynewelm;
							xnvtx[location + 1] = xnewelm;
						}
						// Сортировка по возрастанию X. Часть 3.
						for (integer i28 = 5; i28 < 6; i28++) {
							integer inewelm = invtx[i28];
							doublereal znewelm = znvtx[i28];
							doublereal ynewelm = ynvtx[i28];
							doublereal xnewelm = xnvtx[i28];
							integer location = i28 - 1;
							while ((location >= 4) && (xnvtx[location] > xnewelm)) {
								// сдвигаем все элементы большие очередного.
								invtx[location + 1] = invtx[location];
								znvtx[location + 1] = znvtx[location];
								ynvtx[location + 1] = ynvtx[location];
								xnvtx[location + 1] = xnvtx[location];
								location--;
							}
							invtx[location + 1] = inewelm;
							znvtx[location + 1] = znewelm;
							ynvtx[location + 1] = ynewelm;
							xnvtx[location + 1] = xnewelm;
						}
						// Сортировка по возрастанию X. Часть 4.
						for (integer i28 = 7; i28 < 8; i28++) {
							integer inewelm = invtx[i28];
							doublereal znewelm = znvtx[i28];
							doublereal ynewelm = ynvtx[i28];
							doublereal xnewelm = xnvtx[i28];
							integer location = i28 - 1;
							while ((location >= 6) && (xnvtx[location] < xnewelm)) {
								// сдвигаем все элементы большие очередного.
								invtx[location + 1] = invtx[location];
								znvtx[location + 1] = znvtx[location];
								ynvtx[location + 1] = ynvtx[location];
								xnvtx[location + 1] = xnvtx[location];
								location--;
							}
							invtx[location + 1] = inewelm;
							znvtx[location + 1] = znewelm;
							ynvtx[location + 1] = ynewelm;
							xnvtx[location + 1] = xnewelm;
						}

						fprintf(fp_4, "%lld %lld %lld %lld %lld %lld %lld %lld \n", invtx[0], invtx[1], invtx[2], invtx[3], invtx[4], invtx[5], invtx[6], invtx[7]);
						if (nvtx[0][i] < 1) printf("bad nvtx[0][%lld]=%lld", i, nvtx[0][i]);
						if (nvtx[1][i] < 1) printf("bad nvtx[1][%lld]=%lld", i, nvtx[1][i]);
						if (nvtx[2][i] < 1) printf("bad nvtx[2][%lld]=%lld", i, nvtx[2][i]);
						if (nvtx[3][i] < 1) printf("bad nvtx[3][%lld]=%lld", i, nvtx[3][i]);
						if (nvtx[4][i] < 1) printf("bad nvtx[4][%lld]=%lld", i, nvtx[4][i]);
						if (nvtx[5][i] < 1) printf("bad nvtx[5][%lld]=%lld", i, nvtx[5][i]);
						if (nvtx[6][i] < 1) printf("bad nvtx[6][%lld]=%lld", i, nvtx[6][i]);
						if (nvtx[7][i] < 1) printf("bad nvtx[7][%lld]=%lld", i, nvtx[7][i]);
					}
				}
				else {
					if ((t.whot_is_block[i] < 0) || (t.whot_is_block[i] >= lb)) {
						printf("ERROR in function ANES_tecplot360_export_temperature.\n");
						printf("(t.whot_is_block[i] < 0) || (t.whot_is_block[i] >= lb)\n");
						system("PAUSE");
						exit(1);
					}
					if (b[t.whot_is_block[i]].bvisible) {
						// Учитываем только те ячейки которые видимы, 
						// а именно user указал для них в интерфейсе 
						// bvisible == true.
						fprintf(fp_4, "%lld %lld %lld %lld %lld %lld %lld %lld \n", nvtx[0][i], nvtx[1][i], nvtx[2][i], nvtx[3][i], nvtx[4][i], nvtx[5][i], nvtx[6][i], nvtx[7][i]);
					}
				}


#else

				if (1) {
					integer invtx[8];
					doublereal znvtx[8];
					doublereal ynvtx[8];
					doublereal xnvtx[8];
					for (integer j28 = 0; j28 < 8; j28++) {
						invtx[j28] = nvtx[j28][i];
						znvtx[j28] = pa[nvtx[j28][i] - 1].z;
						ynvtx[j28] = pa[nvtx[j28][i] - 1].y;
						xnvtx[j28] = pa[nvtx[j28][i] - 1].x;
					}

					// Сортировка по возрастанию z.
					for (integer i28 = 1; i28 < 8; i28++) {
						integer inewelm = invtx[i28];
						doublereal znewelm = znvtx[i28];
						doublereal ynewelm = ynvtx[i28];
						doublereal xnewelm = xnvtx[i28];
						integer location = i28 - 1;
						while ((location >= 0) && (znvtx[location] > znewelm)) {
							// сдвигаем все элементы большие очередного.
							invtx[location + 1] = invtx[location];
							znvtx[location + 1] = znvtx[location];
							ynvtx[location + 1] = ynvtx[location];
							xnvtx[location + 1] = xnvtx[location];
							location--;
						}
						invtx[location + 1] = inewelm;
						znvtx[location + 1] = znewelm;
						ynvtx[location + 1] = ynewelm;
						xnvtx[location + 1] = xnewelm;
					}
					// Сортировка по возрастанию y. Часть 1.
					for (integer i28 = 1; i28 < 4; i28++) {
						integer inewelm = invtx[i28];
						doublereal znewelm = znvtx[i28];
						doublereal ynewelm = ynvtx[i28];
						doublereal xnewelm = xnvtx[i28];
						integer location = i28 - 1;
						while ((location >= 0) && (ynvtx[location] > ynewelm)) {
							// сдвигаем все элементы большие очередного.
							invtx[location + 1] = invtx[location];
							znvtx[location + 1] = znvtx[location];
							ynvtx[location + 1] = ynvtx[location];
							xnvtx[location + 1] = xnvtx[location];
							location--;
						}
						invtx[location + 1] = inewelm;
						znvtx[location + 1] = znewelm;
						ynvtx[location + 1] = ynewelm;
						xnvtx[location + 1] = xnewelm;
					}
					// Сортировка по возрастанию y. Часть 2.
					for (integer i28 = 5; i28 < 8; i28++) {
						integer inewelm = invtx[i28];
						doublereal znewelm = znvtx[i28];
						doublereal ynewelm = ynvtx[i28];
						doublereal xnewelm = xnvtx[i28];
						integer location = i28 - 1;
						while ((location >= 4) && (ynvtx[location] > ynewelm)) {
							// сдвигаем все элементы большие очередного.
							invtx[location + 1] = invtx[location];
							znvtx[location + 1] = znvtx[location];
							ynvtx[location + 1] = ynvtx[location];
							xnvtx[location + 1] = xnvtx[location];
							location--;
						}
						invtx[location + 1] = inewelm;
						znvtx[location + 1] = znewelm;
						ynvtx[location + 1] = ynewelm;
						xnvtx[location + 1] = xnewelm;
					}
					// Сортировка по возрастанию X. Часть 1.
					for (integer i28 = 1; i28 < 2; i28++) {
						integer inewelm = invtx[i28];
						doublereal znewelm = znvtx[i28];
						doublereal ynewelm = ynvtx[i28];
						doublereal xnewelm = xnvtx[i28];
						integer location = i28 - 1;
						while ((location >= 0) && (xnvtx[location] > xnewelm)) {
							// сдвигаем все элементы большие очередного.
							invtx[location + 1] = invtx[location];
							znvtx[location + 1] = znvtx[location];
							ynvtx[location + 1] = ynvtx[location];
							xnvtx[location + 1] = xnvtx[location];
							location--;
						}
						invtx[location + 1] = inewelm;
						znvtx[location + 1] = znewelm;
						ynvtx[location + 1] = ynewelm;
						xnvtx[location + 1] = xnewelm;
					}
					// Сортировка по возрастанию X. Часть 2.
					for (integer i28 = 3; i28 < 4; i28++) {
						integer inewelm = invtx[i28];
						doublereal znewelm = znvtx[i28];
						doublereal ynewelm = ynvtx[i28];
						doublereal xnewelm = xnvtx[i28];
						integer location = i28 - 1;
						while ((location >= 2) && (xnvtx[location] < xnewelm)) {
							// сдвигаем все элементы большие очередного.
							invtx[location + 1] = invtx[location];
							znvtx[location + 1] = znvtx[location];
							ynvtx[location + 1] = ynvtx[location];
							xnvtx[location + 1] = xnvtx[location];
							location--;
						}
						invtx[location + 1] = inewelm;
						znvtx[location + 1] = znewelm;
						ynvtx[location + 1] = ynewelm;
						xnvtx[location + 1] = xnewelm;
					}
					// Сортировка по возрастанию X. Часть 3.
					for (integer i28 = 5; i28 < 6; i28++) {
						integer inewelm = invtx[i28];
						doublereal znewelm = znvtx[i28];
						doublereal ynewelm = ynvtx[i28];
						doublereal xnewelm = xnvtx[i28];
						integer location = i28 - 1;
						while ((location >= 4) && (xnvtx[location] > xnewelm)) {
							// сдвигаем все элементы большие очередного.
							invtx[location + 1] = invtx[location];
							znvtx[location + 1] = znvtx[location];
							ynvtx[location + 1] = ynvtx[location];
							xnvtx[location + 1] = xnvtx[location];
							location--;
						}
						invtx[location + 1] = inewelm;
						znvtx[location + 1] = znewelm;
						ynvtx[location + 1] = ynewelm;
						xnvtx[location + 1] = xnewelm;
					}
					// Сортировка по возрастанию X. Часть 4.
					for (integer i28 = 7; i28 < 8; i28++) {
						integer inewelm = invtx[i28];
						doublereal znewelm = znvtx[i28];
						doublereal ynewelm = ynvtx[i28];
						doublereal xnewelm = xnvtx[i28];
						integer location = i28 - 1;
						while ((location >= 6) && (xnvtx[location] < xnewelm)) {
							// сдвигаем все элементы большие очередного.
							invtx[location + 1] = invtx[location];
							znvtx[location + 1] = znvtx[location];
							ynvtx[location + 1] = ynvtx[location];
							xnvtx[location + 1] = xnvtx[location];
							location--;
						}
						invtx[location + 1] = inewelm;
						znvtx[location + 1] = znewelm;
						ynvtx[location + 1] = ynewelm;
						xnvtx[location + 1] = xnewelm;
					}

					fprintf_s(fp_4, "%d %d %d %d %d %d %d %d \n", invtx[0], invtx[1], invtx[2], invtx[3], invtx[4], invtx[5], invtx[6], invtx[7]);

					//debug
					//printf("%d %d %d %d %d %d %d %d \n", invtx[0], invtx[1], invtx[2], invtx[3], invtx[4], invtx[5], invtx[6], invtx[7]);
					//printf( "%d %d %d %d %d %d %d %d \n", nvtx[0][i], nvtx[1][i], nvtx[2][i], nvtx[3][i], nvtx[4][i], nvtx[5][i], nvtx[6][i], nvtx[7][i]);
					//getchar();

					if (nvtx[0][i] < 1) printf("bad nvtx[0][%d]=%d", i, nvtx[0][i]);
					if (nvtx[1][i] < 1) printf("bad nvtx[1][%d]=%d", i, nvtx[1][i]);
					if (nvtx[2][i] < 1) printf("bad nvtx[2][%d]=%d", i, nvtx[2][i]);
					if (nvtx[3][i] < 1) printf("bad nvtx[3][%d]=%d", i, nvtx[3][i]);
					if (nvtx[4][i] < 1) printf("bad nvtx[4][%d]=%d", i, nvtx[4][i]);
					if (nvtx[5][i] < 1) printf("bad nvtx[5][%d]=%d", i, nvtx[5][i]);
					if (nvtx[6][i] < 1) printf("bad nvtx[6][%d]=%d", i, nvtx[6][i]);
					if (nvtx[7][i] < 1) printf("bad nvtx[7][%d]=%d", i, nvtx[7][i]);
			}
				else {
					fprintf(fp_4, "%d %d %d %d %d %d %d %d \n", nvtx[0][i], nvtx[1][i], nvtx[2][i], nvtx[3][i], nvtx[4][i], nvtx[5][i], nvtx[6][i], nvtx[7][i]);
				}
#endif
		}
			
		}
		fclose(fp_4);
		//WinExec("C:\\Program Files\\Tecplot\\Tecplot 360 EX 2014 R1\\bin\\tec360.exe ALICEFLOW0_24ALICEMESH.PLT", SW_NORMAL);
	}

} // ANES_tecplot360_export_temperature

// Переинтерполляция температуры с сетки на сетку. Конкретно с АЛИС сетки на обычную сетку.
// t_buf->potent
void ALICE_2_Structural(integer maxnod, TOCHKA* pa,
	integer maxelm, integer** nvtx, doublereal* &potent,
	doublereal* &x47, doublereal* &y47, doublereal* &z47,
	doublereal* &temp47, integer** &nvtx47,
	integer &m_sizeT, integer &maxelm47, 
	doublereal operatingtemperature_loc) {

	// Инициализация компонент скорости во внутренности расчётной области.
	// 26 марта 2017.

	

	// Метод линейного порядка.
	doublereal min_x = 1e60;
	doublereal min_y = 1e60;
	doublereal min_z = 1e60;
	doublereal max_x = -1e60;
	doublereal max_y = -1e60;
	doublereal max_z = -1e60;

	for (integer i_48 = 0; i_48 < maxelm47; i_48++) {
		for (integer i_47 = 0; i_47 < 8; i_47++) {
			if (x47[nvtx47[i_47][i_48]] < min_x) {
				min_x = x47[nvtx47[i_47][i_48]];
			}
			if (y47[nvtx47[i_47][i_48]] < min_y) {
				min_y = y47[nvtx47[i_47][i_48]];
			}
			if (z47[nvtx47[i_47][i_48]] < min_z) {
				min_z = z47[nvtx47[i_47][i_48]];
			}
			if (x47[nvtx47[i_47][i_48]] > max_x) {
				max_x = x47[nvtx47[i_47][i_48]];
			}
			if (y47[nvtx47[i_47][i_48]] > max_y) {
				max_y = y47[nvtx47[i_47][i_48]];
			}
			if (z47[nvtx47[i_47][i_48]] > max_z) {
				max_z = z47[nvtx47[i_47][i_48]];
			}
		}
	}

	min_x = 1.05*fabs(max_x - min_x);
	if (min_x < 1.0e-40) {
		min_x = 1.05*fabs(max_x);
	}
	min_y = 1.05*fabs(max_y - min_y);
	if (min_y < 1.0e-40) {
		min_y = 1.05*fabs(max_y);
	}
	min_z = 1.05*fabs(max_z - min_z);
	if (min_z < 1.0e-40) {
		min_z = 1.05*fabs(max_z);
	}


	TOCHKA** pointerlist = new TOCHKA*[maxelm];
	doublereal** rthdsd_Gauss = new doublereal*[maxelm];
	for (integer i_47 = 0; i_47 < maxelm; i_47++) {
		pointerlist[i_47] = new TOCHKA[8];
		rthdsd_Gauss[i_47] = new doublereal[8];
	}

	integer istart_i47 = maxelm47 - 1, ic76 = 0, ih64=0;

	for (integer i = 0; i < maxelm; i++) {
		doublereal xc47, yc47, zc47;

		TOCHKA p;
		center_cord3D(i, nvtx, pa, p, 100);
		xc47 = p.x;
		yc47 = p.y;
		zc47 = p.z;

		// Найдем номер конечного изопараметрического элемента 
		// в котором содержится точка p.
		bool bfound = false;
		integer ifound = NON_EXISTENT_NODE;

		// Эвристика 02.01.2018 для существенного ускорения поиска при интерполляции.
		// Наиболее часто, следующий найденный элемент оказывается меньше либо равен предыдущему.
		// Первый вызов основан на этой эвристике и хочет быстро отсечь большинство сканирований.
		// Только в том случае когда он ничего не находит запускается повторный вызов.
		// Это работает потому что элеменеты -nvtx ячейки были занумерованы до интерполляции 
		// согласно определенному порядку обхода. Этот же порядок обхода в общем и целом справелив и для ячеек 
		// полученных после на обычной декартовой прямоугольной структурированной сетке. 
		// Здесь мы переводим информацию с АЛИС сетки на обычную декартовую
		// прямоугольную структурированную сетку с главной целью правильно и точно вычислить градиенты. Что проблематично на АЛИС.
		// 2.01.2018 Необходимо увеличить точность (доработать) сборку матрицы на АЛИС сетке.
		// АЛИС сетка необходима для увеличения скорости проводимого моделирования без существенной потери в точности.
		const doublereal tolerance_eps = 0.05; // 5%

		for (integer i_47 = istart_i47; i_47 >= 0; i_47--) {
			doublereal Xmin, Xmax, Ymin, Ymax, Zmin, Zmax; 	
			if (0) {
				Xmin = 1.0e60;
				Xmax = -1.0e60;
				for (integer i_7 = 0; i_7 < 8; i_7++) {
					if (x47[nvtx47[i_7][i_47]] < Xmin) Xmin = x47[nvtx47[i_7][i_47]];
					if (x47[nvtx47[i_7][i_47]] > Xmax) Xmax = x47[nvtx47[i_7][i_47]];
				}
				Ymin = 1.0e60;
				Ymax = -1.0e60;
				for (integer i_7 = 0; i_7 < 8; i_7++) {
					if (y47[nvtx47[i_7][i_47]] < Ymin) Ymin = y47[nvtx47[i_7][i_47]];
					if (y47[nvtx47[i_7][i_47]] > Ymax) Ymax = y47[nvtx47[i_7][i_47]];
				}
				Zmin = 1.0e60;
				Zmax = -1.0e60;
				for (integer i_7 = 0; i_7 < 8; i_7++) {
					if (z47[nvtx47[i_7][i_47]] < Zmin) Zmin = z47[nvtx47[i_7][i_47]];
					if (z47[nvtx47[i_7][i_47]] > Zmax) Zmax = z47[nvtx47[i_7][i_47]];
				}
			}
			else {
				Xmin = x47[nvtx47[0][i_47]];
				Xmax = x47[nvtx47[1][i_47]];
				Ymin = y47[nvtx47[0][i_47]];
				Ymax = y47[nvtx47[3][i_47]];
				Zmin = z47[nvtx47[0][i_47]];
				Zmax = z47[nvtx47[4][i_47]];
			}

			

			doublereal distance=0.0;
			distance = fabs(Xmax - Xmin);
			Xmax += tolerance_eps*distance;
			Xmin -= tolerance_eps*distance;
			distance = fabs(Ymax - Ymin);
			Ymax += tolerance_eps*distance;
			Ymin -= tolerance_eps*distance;
			distance = fabs(Zmax - Zmin);
			Zmax += tolerance_eps*distance;
			Zmin -= tolerance_eps*distance;

			
			if ((xc47 >= Xmin) && (xc47 <= Xmax) && (yc47 >= Ymin) && (yc47 <= Ymax) &&
				(zc47 >= Zmin) && (zc47 <= Zmax))
			{
				ifound = i_47;
				bfound = true;
				break;
			}
		}
		// Повторный вызов в случае если ничего найти не удалось. 
		// Ресурсоёмко но запускается редко.
		if (!(bfound)) {
			for (integer i_47 = 0; i_47 < maxelm47; i_47++) {
				doublereal Xmin, Xmax, Ymin, Ymax, Zmin, Zmax;
				if (0) {
					Xmin = 1.0e60;
					Xmax = -1.0e60;
					for (integer i_7 = 0; i_7 < 8; i_7++) {
						if (x47[nvtx47[i_7][i_47]] < Xmin) Xmin = x47[nvtx47[i_7][i_47]];
						if (x47[nvtx47[i_7][i_47]] > Xmax) Xmax = x47[nvtx47[i_7][i_47]];
					}
					Ymin = 1.0e60;
					Ymax = -1.0e60;
					for (integer i_7 = 0; i_7 < 8; i_7++) {
						if (y47[nvtx47[i_7][i_47]] < Ymin) Ymin = y47[nvtx47[i_7][i_47]];
						if (y47[nvtx47[i_7][i_47]] > Ymax) Ymax = y47[nvtx47[i_7][i_47]];
					}
					Zmin = 1.0e60;
					Zmax = -1.0e60;
					for (integer i_7 = 0; i_7 < 8; i_7++) {
						if (z47[nvtx47[i_7][i_47]] < Zmin) Zmin = z47[nvtx47[i_7][i_47]];
						if (z47[nvtx47[i_7][i_47]] > Zmax) Zmax = z47[nvtx47[i_7][i_47]];
					}
				}
				else {
					Xmin = x47[nvtx47[0][i_47]];
					Xmax = x47[nvtx47[1][i_47]];
					Ymin = y47[nvtx47[0][i_47]];
					Ymax = y47[nvtx47[3][i_47]];
					Zmin = z47[nvtx47[0][i_47]];
					Zmax = z47[nvtx47[4][i_47]];
				}

				doublereal distance = 0.0;
				distance = fabs(Xmax - Xmin);
				Xmax += tolerance_eps*distance;
				Xmin -= tolerance_eps*distance;
				distance = fabs(Ymax - Ymin);
				Ymax += tolerance_eps*distance;
				Ymin -= tolerance_eps*distance;
				distance = fabs(Zmax - Zmin);
				Zmax += tolerance_eps*distance;
				Zmin -= tolerance_eps*distance;


				if ((xc47 >= Xmin) && (xc47 <= Xmax) && (yc47 >= Ymin) && (yc47 <= Ymax) &&
					(zc47 >= Zmin) && (zc47 <= Zmax))
				{
					ifound = i_47;
					bfound = true;
					break;
				}
			}
		}
		/*
		// Экспериментально установлено что ih64 почти в точности равно maxelm47 и с этой точки зрения ничего сэкономить нельзя.
		if (ifound > istart_i47) {
			ic76++;
			//printf("monotonnost narushena\n");
			//getchar();
			if (fabs(ifound - istart_i47) > ih64) {
				ih64 = ifound - istart_i47;
			}
		}
		*/
		// Запоминаем номер последнего найденного чтобы существенно ускорить поиск.
		istart_i47 = ifound;


		// МНК
		// Три раза, для каждой компоненты скорости.
		doublereal temp_47 = 0.0;

		if (bfound) {

			//TOCHKA p;
			//center_cord3D(i, nvtx, pa, p, 100);
			p.x = p.x + min_x;
			p.y = p.y + min_y;
			p.z = p.z + min_z;

			// VX

			for (integer j = 0; j <= 7; j++) {
				TOCHKA p1;
				p1.x = x47[nvtx47[j][ifound]];
				p1.y = y47[nvtx47[j][ifound]];
				p1.z = z47[nvtx47[j][ifound]];
				p1.x = p1.x + min_x;
				p1.y = p1.y + min_y;
				p1.z = p1.z + min_z;

				pointerlist[i][j] = p1;
				if (fabs(p1.x) < 1.0e-40) {
					printf("problem x=%e\n", p1.x);
					system("PAUSE");
				}
				if (fabs(p1.y) < 1.0e-40) {
					printf("problem y=%e\n", p1.y);
					system("PAUSE");
				}
				if (fabs(p1.z) < 1.0e-40) {
					printf("problem z=%e\n", p1.z);
					system("PAUSE");
				}
				rthdsd_Gauss[i][j] = temp47[nvtx47[j][ifound]];
			}

			doublereal** Xmatr = new doublereal*[4];
			for (integer j = 0; j <= 3; j++) {
				Xmatr[j] = new doublereal[4];
			}


			doublereal* bmatr = new doublereal[4];
			doublereal* koefmatr = new doublereal[4];

			for (integer j1 = 0; j1 <= 3; j1++) {
				for (integer j2 = 0; j2 <= 3; j2++) {
					Xmatr[j1][j2] = 0.0;
				}
				bmatr[j1] = 0.0;
				koefmatr[j1] = 0.0;
			}




			for (integer j = 0; j < 8; j++) {

				Xmatr[0][0] += 1.0;
				Xmatr[0][1] += pointerlist[i][j].x;
				Xmatr[0][2] += pointerlist[i][j].y;
				Xmatr[0][3] += pointerlist[i][j].z;

				Xmatr[1][0] += pointerlist[i][j].x;
				Xmatr[1][1] += pointerlist[i][j].x*pointerlist[i][j].x;
				Xmatr[1][2] += pointerlist[i][j].x*pointerlist[i][j].y;
				Xmatr[1][3] += pointerlist[i][j].x*pointerlist[i][j].z;

				Xmatr[2][0] += pointerlist[i][j].y;
				Xmatr[2][1] += pointerlist[i][j].y*pointerlist[i][j].x;
				Xmatr[2][2] += pointerlist[i][j].y*pointerlist[i][j].y;
				Xmatr[2][3] += pointerlist[i][j].y*pointerlist[i][j].z;

				Xmatr[3][0] += pointerlist[i][j].z;
				Xmatr[3][1] += pointerlist[i][j].z*pointerlist[i][j].x;
				Xmatr[3][2] += pointerlist[i][j].z*pointerlist[i][j].y;
				Xmatr[3][3] += pointerlist[i][j].z*pointerlist[i][j].z;

				bmatr[0] += rthdsd_Gauss[i][j];
				bmatr[1] += pointerlist[i][j].x*rthdsd_Gauss[i][j];
				bmatr[2] += pointerlist[i][j].y*rthdsd_Gauss[i][j];
				bmatr[3] += pointerlist[i][j].z*rthdsd_Gauss[i][j];
			}


			for (integer j1 = 0; j1 <= 100; j1++) {
				koefmatr[0] = (bmatr[0] - Xmatr[0][1] * koefmatr[1] - Xmatr[0][2] * koefmatr[2] - Xmatr[0][3] * koefmatr[3]) / Xmatr[0][0];
				koefmatr[1] = (bmatr[1] - Xmatr[1][0] * koefmatr[0] - Xmatr[1][2] * koefmatr[2] - Xmatr[1][3] * koefmatr[3]) / Xmatr[1][1];
				koefmatr[2] = (bmatr[2] - Xmatr[2][0] * koefmatr[0] - Xmatr[2][1] * koefmatr[1] - Xmatr[2][3] * koefmatr[3]) / Xmatr[2][2];
				koefmatr[3] = (bmatr[3] - Xmatr[3][0] * koefmatr[0] - Xmatr[3][1] * koefmatr[1] - Xmatr[3][2] * koefmatr[2]) / Xmatr[3][3];
			}
			temp_47 = koefmatr[0] + koefmatr[1] * (p.x) + koefmatr[2] * (p.y) + koefmatr[3] * (p.z);

			for (integer j = 0; j <= 3; j++) {
				delete[] Xmatr[j];
			}
			delete[] Xmatr;
			delete[] bmatr;
			delete[] koefmatr;

			potent[i] = temp_47;
		}
		else {
			printf("WARNING : ALICE_2_Structural node not found!!! bfound==false...\n");
			//system("PAUSE");
			potent[i] = operatingtemperature_loc;
		}

		
		
		}

		// Диагностическая печать для исследования закономерностей в процессе интерполляции.
		// Закономерности используются чтобы существенно ускорить процесс интерполляции.
		//if (ic76 > 0) {
			//printf("monotonnost narushena=%d %d h=%d\n",ic76, maxelm47,ih64);
			//getchar();
		//}

		for (integer i = 0; i < maxelm; i++) {
			delete[] pointerlist[i];
			delete[] rthdsd_Gauss[i];
		}
		delete[] pointerlist;
		delete[] rthdsd_Gauss;

		// Освобождение оперативной памяти.
		if (x47 != NULL) {
			delete[] x47;
			x47 = NULL;
		}
		if (y47 != NULL) {
			delete[] y47;
			y47 = NULL;
		}
		if (z47 != NULL) {
			delete[] z47;
			z47 = NULL;
		}
		if (temp47 != NULL) {
			delete[] temp47;
			temp47 = NULL;
		}
		
		if (nvtx47 != NULL) {
			for (integer i_47 = 0; i_47 < 8; i_47++) {
				if (nvtx47[i_47] != NULL) {
					delete[] nvtx47[i_47];
					nvtx47[i_47] = NULL;
				}
			}
			delete[] nvtx47;
			nvtx47 = NULL;
		}

		m_sizeT = 0;
		maxelm47 = 0;

		printf("done.\n");

}//ALICE_2_Structural

  // Экспорт в программу техплот температуры.
  //С АЛИС сетки.
void ANES_tecplot360_export_temperature_preobrazovatel(integer maxnod, TOCHKA* pa,
	integer maxelm, integer** nvtx, doublereal* potent, TEMPER &t, 
	doublereal* &x_buf, doublereal* &y_buf, doublereal* &z_buf, doublereal* &t_buf, integer** &nvtx_buf, integer &m_sizeT, integer &m_size_nvtx,
	doublereal operating_temperature) {

	// 2 ноября 2016. машинное эпсилон.
	//doublereal eps_mashine = 1.0e-44; // float
	doublereal eps_mashine = 1.0e-308; // double


									   // nvtx && pa сформированы, можно экспортировать в tecplot360
	//FILE *fp_4 = NULL;
	//errno_t err_4;
	//err_4 = fopen_s(&fp_4, "ALICEFLOW0_07_temp.PLT", "w");

	//if ((err_4) != 0) {
		//printf("Create File temp Error\n");
		//getchar();
		//system("pause");

	//}
	//else {

		//fprintf(fp_4, "TITLE = \"ALICEFLOW0_24\"\n");
		//fprintf(fp_4, "VARIABLES = x, y, z, Temp, Lam, log10_heat_flux_X, log10_heat_flux_Y, log10_heat_flux_Z, log10_heat_flux_mag\n");
#if doubleintprecision == 1
		//fprintf(fp_4, "ZONE T=\"Rampant\", N=%lld, E=%lld, ET=BRICK, F=FEBLOCK\n\n", maxnod, maxelm);
#else
		//fprintf(fp_4, "ZONE T=\"Rampant\", N=%d, E=%d, ET=BRICK, F=FEBLOCK\n\n", maxnod, maxelm);
#endif

	m_sizeT = maxnod;
	m_size_nvtx = maxelm;

	x_buf = new doublereal[maxnod];
	y_buf = new doublereal[maxnod];
	z_buf = new doublereal[maxnod];
	t_buf = new doublereal[maxnod];

	nvtx_buf = new integer*[8];
	for (integer i = 0; i < 8; i++) {
		nvtx_buf[i] = new integer[maxelm];
	}

		// запись x
		for (integer i = 0; i < maxnod; i++) {
			//fprintf(fp_4, "%+.16f ", pa[i].x);
			x_buf[i] = pa[i].x;
			//if (i % 10 == 0) fprintf(fp_4, "\n");
		}
		//fprintf(fp_4, "\n");
		// запись y
		for (integer i = 0; i < maxnod; i++) {
			//fprintf(fp_4, "%+.16f ", pa[i].y);
			//if (i % 10 == 0) fprintf(fp_4, "\n");
			y_buf[i] = pa[i].y;
		}
		//fprintf(fp_4, "\n");
		// запись z
		for (integer i = 0; i < maxnod; i++) {
			//fprintf(fp_4, "%+.16f ", pa[i].z);
			//if (i % 10 == 0) fprintf(fp_4, "\n");
			z_buf[i] = pa[i].z;
		}
		//fprintf(fp_4, "\n");
		doublereal* temp = new doublereal[maxnod];
		doublereal* lam = new doublereal[maxnod];
		doublereal* vol = new doublereal[maxnod];
		doublereal* heat_flux_X = new doublereal[maxnod];
		doublereal* heat_flux_Y = new doublereal[maxnod];
		doublereal* heat_flux_Z = new doublereal[maxnod];
		doublereal* heat_flux_mag = new doublereal[maxnod];

		for (integer i = 0; i < maxnod; i++) {
			temp[i] = 0.0;
			vol[i] = 0.0;
			lam[i] = 0.0;
			heat_flux_X[i] = 0.0;
			heat_flux_Y[i] = 0.0;
			heat_flux_Z[i] = 0.0;
			heat_flux_mag[i] = 0.0;
		}
		/*
		// Здесь мы ловим конкретную исключительную ситуацию.
		bool bfound = false;
		for (integer i = 0; i <= maxelm - 1; i++) {
		for (integer j = 0; j <= 7; j++) {
		if (nvtx[j][i] - 1 == 69462) bfound = true;
		if (nvtx[j][i] - 1 == 69463) bfound = true;
		if (nvtx[j][i] - 1 == 69464) bfound = true;
		if (nvtx[j][i] - 1 == 69477) bfound = true;
		if (nvtx[j][i] - 1 == 69478) bfound = true;
		}
		}
		if (bfound) {
		printf("bfound\n");
		}
		else
		{
		printf("notfound\n");
		}
		*/

		const integer IORDER_METHOD = 3;

		if (0== IORDER_METHOD) {
			// Метод нулевого порядка.

			for (integer i = 0; i <= maxelm - 1; i++) {
				// вычисление размеров текущего контрольного объёма:
				doublereal dx = 0.0, dy = 0.0, dz = 0.0;// объём текущего контроольного объёма
				volume3D(i, nvtx, pa, dx, dy, dz);

				for (integer j = 0; j <= 7; j++) {
					vol[nvtx[j][i] - 1] += dx*dy*dz;
					temp[nvtx[j][i] - 1] += dx*dy*dz*potent[i];
				}
			}
			for (integer i = 0; i < maxnod; i++) {
				if (fabs(vol[i]) > eps_mashine) {
					temp[i] = temp[i] / vol[i];
				}
				else {
#if doubleintprecision == 1
					printf("fatal error! ANES_tecplot_export in module constr_struct_alice.cpp. i=%lld maxnod=%lld\n", i, maxnod);
					printf("vol[%lld]==%e\n", i, vol[i]);
#else
					printf("fatal error! ANES_tecplot_export in module constr_struct_alice.cpp. i=%d maxnod=%d\n", i, maxnod);
					printf("vol[%d]==%e\n", i, vol[i]);
#endif

					//getchar();
					system("PAUSE");
					temp[i] = 0.0;
				}
			}
		}
		else if (1 == IORDER_METHOD) {

			// Метод линейного порядка.
			doublereal min_x = 1e60;
			doublereal min_y = 1e60;
			doublereal min_z = 1e60;
			doublereal max_x = -1e60;
			doublereal max_y = -1e60;
			doublereal max_z = -1e60;

			for (integer i = 0; i < maxnod; i++) {
				if (pa[i].x < min_x) {
					min_x = pa[i].x;
				}
				if (pa[i].y < min_y) {
					min_y = pa[i].y;
				}
				if (pa[i].z < min_z) {
					min_z = pa[i].z;
				}
				if (pa[i].x > max_x) {
					max_x = pa[i].x;
				}
				if (pa[i].y > max_y) {
					max_y = pa[i].y;
				}
				if (pa[i].z > max_z) {
					max_z = pa[i].z;
				}
			}

			//min_x *= 1.2;
			//min_y *= 1.2;
			//min_z *= 1.2;



			min_x = 1.05*fabs(max_x - min_x);
			if (min_x < 1.0e-30) {
				min_x = 1.05*fabs(max_x);
			}
			min_y = 1.05*fabs(max_y - min_y);
			if (min_y < 1.0e-30) {
				min_y = 1.05*fabs(max_y);
			}
			min_z = 1.05*fabs(max_z - min_z);
			if (min_z < 1.0e-30) {
				min_z = 1.05*fabs(max_z);
			}


			/*
			if (min_x < 1.0e-30) {
			printf("error!!! negative min_x MNK!\n");
			printf("min_x=%e max_x=%e\n",min_x,max_x);
			}
			if (min_y < 1.0e-30) {
			printf("error!!! negative min_y MNK!\n");
			printf("min_y=%e max_y=%e\n", min_y, max_y);
			}
			if (min_z < 1.0e-30) {
			printf("error!!! negative min_z MNK!\n");
			printf("min_z=%e max_z=%e\n", min_z, max_z);
			}
			*/

			integer* inum_now = new integer[maxnod];

			for (integer i = 0; i < maxnod; i++) {
				temp[i] = 0.0;
				vol[i] = 0.0;
				inum_now[i] = 0;
			}

			for (integer i = 0; i <= maxelm - 1; i++) {
				for (integer j = 0; j <= 7; j++) {
					inum_now[nvtx[j][i] - 1] += 1;
				}
			}

			for (integer i = 0; i < maxnod; i++) {
				if (inum_now[i] < 1) {
#if doubleintprecision == 1
					printf("i=%lld maxnod=%lld inum_now[%lld]=%lld\n", i, maxnod, i, inum_now[i]);
#else
					printf("i=%d maxnod=%d inum_now[%d]=%d\n", i, maxnod, i, inum_now[i]);
#endif

					system("pause");
				}
			}

			TOCHKA** pointerlist = new TOCHKA*[maxnod];
			doublereal** rthdsd_Gauss = new doublereal*[maxnod];
			for (integer i = 0; i < maxnod; i++) {
				pointerlist[i] = new TOCHKA[(inum_now[i])];
				rthdsd_Gauss[i] = new doublereal[(inum_now[i])];
			}


			for (integer i = 0; i < maxnod; i++) {
				inum_now[i] = 0;
			}

			for (integer i = 0; i <= maxelm - 1; i++) {
				// вычисление размеров текущего контрольного объёма:
				doublereal dx = 0.0, dy = 0.0, dz = 0.0;// объём текущего контроольного объёма
				volume3D(i, nvtx, pa, dx, dy, dz);

				for (integer j = 0; j <= 7; j++) {
					vol[nvtx[j][i] - 1] += dx*dy*dz;
					temp[nvtx[j][i] - 1] += dx*dy*dz*potent[i];
				}

				TOCHKA p;
				center_cord3D(i, nvtx, pa, p, 100);
				p.x = p.x + min_x;
				p.y = p.y + min_y;
				p.z = p.z + min_z;


				for (integer j = 0; j <= 7; j++) {
					pointerlist[nvtx[j][i] - 1][inum_now[nvtx[j][i] - 1]] = p;
					if (fabs(p.x) < 1.0e-30) {
						printf("problem x=%e\n", p.x);
						system("pause");
					}
					if (fabs(p.y) < 1.0e-30) {
						printf("problem y=%e\n", p.y);
						system("pause");
					}
					if (fabs(p.z) < 1.0e-30) {
						printf("problem z=%e\n", p.z);
						system("pause");
					}
					rthdsd_Gauss[nvtx[j][i] - 1][inum_now[nvtx[j][i] - 1]] = potent[i];
					inum_now[nvtx[j][i] - 1]++;
				}


				/*
				for (integer j = 0; j <= 7; j++) {
				vol[nvtx[j][i] - 1] += dx*dy*dz;
				temp[nvtx[j][i] - 1] += dx*dy*dz*potent[i];
				}
				*/
			}


			//integer jcontrol = 0;
			for (integer i = 0; i < maxnod; i++) {
				if (fabs(vol[i]) > eps_mashine) {


					doublereal** Xmatr = new doublereal*[4];
					for (integer j = 0; j <= 3; j++) {
						Xmatr[j] = new doublereal[4];
					}


					doublereal* bmatr = new doublereal[4];
					doublereal* koefmatr = new doublereal[4];

					for (integer j1 = 0; j1 <= 3; j1++) {
						for (integer j2 = 0; j2 <= 3; j2++) {
							Xmatr[j1][j2] = 0.0;
						}
						bmatr[j1] = 0.0;
						koefmatr[j1] = 0.0;
					}



					for (integer j = 0; j < inum_now[i]; j++) {

						Xmatr[0][0] += 1.0;
						Xmatr[0][1] += pointerlist[i][j].x;
						Xmatr[0][2] += pointerlist[i][j].y;
						Xmatr[0][3] += pointerlist[i][j].z;

						Xmatr[1][0] += pointerlist[i][j].x;
						Xmatr[1][1] += pointerlist[i][j].x*pointerlist[i][j].x;
						Xmatr[1][2] += pointerlist[i][j].x*pointerlist[i][j].y;
						Xmatr[1][3] += pointerlist[i][j].x*pointerlist[i][j].z;

						Xmatr[2][0] += pointerlist[i][j].y;
						Xmatr[2][1] += pointerlist[i][j].y*pointerlist[i][j].x;
						Xmatr[2][2] += pointerlist[i][j].y*pointerlist[i][j].y;
						Xmatr[2][3] += pointerlist[i][j].y*pointerlist[i][j].z;

						Xmatr[3][0] += pointerlist[i][j].z;
						Xmatr[3][1] += pointerlist[i][j].z*pointerlist[i][j].x;
						Xmatr[3][2] += pointerlist[i][j].z*pointerlist[i][j].y;
						Xmatr[3][3] += pointerlist[i][j].z*pointerlist[i][j].z;

						bmatr[0] += rthdsd_Gauss[i][j];
						bmatr[1] += pointerlist[i][j].x*rthdsd_Gauss[i][j];
						bmatr[2] += pointerlist[i][j].y*rthdsd_Gauss[i][j];
						bmatr[3] += pointerlist[i][j].z*rthdsd_Gauss[i][j];

						/*
						Xmatr[0][0] += 1.0;
						Xmatr[0][1] += pointerlist[i][j].x-min_x;
						Xmatr[0][2] += pointerlist[i][j].y-min_y;
						Xmatr[0][3] += pointerlist[i][j].z-min_z;

						Xmatr[1][0] += pointerlist[i][j].x-min_x;
						Xmatr[1][1] += (pointerlist[i][j].x-min_x)*(pointerlist[i][j].x-min_x);
						Xmatr[1][2] += (pointerlist[i][j].x-min_x)*(pointerlist[i][j].y-min_y);
						Xmatr[1][3] += (pointerlist[i][j].x-min_x)*(pointerlist[i][j].z-min_z);

						Xmatr[2][0] += pointerlist[i][j].y-min_y;
						Xmatr[2][1] += (pointerlist[i][j].y-min_y)*(pointerlist[i][j].x-min_x);
						Xmatr[2][2] += (pointerlist[i][j].y-min_y)*(pointerlist[i][j].y-min_y);
						Xmatr[2][3] += (pointerlist[i][j].y-min_y)*(pointerlist[i][j].z-min_z);

						Xmatr[3][0] += (pointerlist[i][j].z-min_z);
						Xmatr[3][1] += (pointerlist[i][j].z-min_z)*(pointerlist[i][j].x-min_x);
						Xmatr[3][2] += (pointerlist[i][j].z-min_z)*(pointerlist[i][j].y-min_y);
						Xmatr[3][3] += (pointerlist[i][j].z-min_z)*(pointerlist[i][j].z-min_z);

						bmatr[0] += rthdsd_Gauss[i][j];
						bmatr[1] += (pointerlist[i][j].x-min_x)*rthdsd_Gauss[i][j];
						bmatr[2] += (pointerlist[i][j].y-min_y)*rthdsd_Gauss[i][j];
						bmatr[3] += (pointerlist[i][j].z-min_z)*rthdsd_Gauss[i][j];
						*/
					}

					if ((fabs(Xmatr[0][0]) < 1.e-30) || (fabs(Xmatr[1][1]) < 1.e-30) || (fabs(Xmatr[2][2]) < 1.e-30) || (fabs(Xmatr[3][3]) < 1.e-30)) {
#if doubleintprecision == 1
						printf("inum_now[%lld]=%lld\n", i, inum_now[i]);
#else
						printf("inum_now[%d]=%d\n", i, inum_now[i]);
#endif

						system("pause");
					}

					//Xmatr*koefmatr = bmatr;
					/*
					if (!my_version_gauss1(Xmatr, 4, bmatr, koefmatr, false, i)) {
					temp[i] = temp[i] / vol[i];
					}
					else {
					temp[i] = koefmatr[0] + koefmatr[1] * (pa[i].x+min_x) + koefmatr[2] * (pa[i].y+min_y) + koefmatr[3] * (pa[i].z+min_z);
					}
					*/
					for (integer j1 = 0; j1 <= 100; j1++) {
						koefmatr[0] = (bmatr[0] - Xmatr[0][1] * koefmatr[1] - Xmatr[0][2] * koefmatr[2] - Xmatr[0][3] * koefmatr[3]) / Xmatr[0][0];
						koefmatr[1] = (bmatr[1] - Xmatr[1][0] * koefmatr[0] - Xmatr[1][2] * koefmatr[2] - Xmatr[1][3] * koefmatr[3]) / Xmatr[1][1];
						koefmatr[2] = (bmatr[2] - Xmatr[2][0] * koefmatr[0] - Xmatr[2][1] * koefmatr[1] - Xmatr[2][3] * koefmatr[3]) / Xmatr[2][2];
						koefmatr[3] = (bmatr[3] - Xmatr[3][0] * koefmatr[0] - Xmatr[3][1] * koefmatr[1] - Xmatr[3][2] * koefmatr[2]) / Xmatr[3][3];
					}
					temp[i] = koefmatr[0] + koefmatr[1] * (pa[i].x + min_x) + koefmatr[2] * (pa[i].y + min_y) + koefmatr[3] * (pa[i].z + min_z);
					heat_flux_X[i] = koefmatr[1];
					heat_flux_Y[i] = koefmatr[2];
					heat_flux_Z[i] = koefmatr[3];


					for (integer j = 0; j <= 3; j++) {
						delete[] Xmatr[j];
					}
					delete[] Xmatr;
					delete[] bmatr;
					delete[] koefmatr;

				}
				else {
#if doubleintprecision == 1
					printf("fatal error! ANES_tecplot_export in module constr_struct_alice.cpp. i=%lld maxnod=%lld\n", i, maxnod);
					printf("vol[%lld]==%e\n", i, vol[i]);
#else
					printf("fatal error! ANES_tecplot_export in module constr_struct_alice.cpp. i=%d maxnod=%d\n", i, maxnod);
					printf("vol[%d]==%e\n", i, vol[i]);
#endif

					//getchar();
					system("PAUSE");
					temp[i] = 0.0;
				}
			}



			delete[] inum_now;
			for (integer i = 0; i < maxnod; i++) {
				delete[] pointerlist[i];
				delete[] rthdsd_Gauss[i];
			}
			delete[] pointerlist;
			delete[] rthdsd_Gauss;
		}
		else if (2 == IORDER_METHOD) {

			// Метод второго порядка.
			doublereal min_x = 1e60;
			doublereal min_y = 1e60;
			doublereal min_z = 1e60;
			doublereal max_x = -1e60;
			doublereal max_y = -1e60;
			doublereal max_z = -1e60;

			for (integer i = 0; i < maxnod; i++) {
				if (pa[i].x < min_x) {
					min_x = pa[i].x;
				}
				if (pa[i].y < min_y) {
					min_y = pa[i].y;
				}
				if (pa[i].z < min_z) {
					min_z = pa[i].z;
				}
				if (pa[i].x > max_x) {
					max_x = pa[i].x;
				}
				if (pa[i].y > max_y) {
					max_y = pa[i].y;
				}
				if (pa[i].z > max_z) {
					max_z = pa[i].z;
				}
			}

			//min_x *= 1.2;
			//min_y *= 1.2;
			//min_z *= 1.2;



			min_x = 1.05*fabs(max_x - min_x);
			if (min_x < 1.0e-30) {
				min_x = 1.05*fabs(max_x);
			}
			min_y = 1.05*fabs(max_y - min_y);
			if (min_y < 1.0e-30) {
				min_y = 1.05*fabs(max_y);
			}
			min_z = 1.05*fabs(max_z - min_z);
			if (min_z < 1.0e-30) {
				min_z = 1.05*fabs(max_z);
			}


			/*
			if (min_x < 1.0e-30) {
			printf("error!!! negative min_x MNK!\n");
			printf("min_x=%e max_x=%e\n", min_x, max_x);
			}
			if (min_y < 1.0e-30) {
			printf("error!!! negative min_y MNK!\n");
			printf("min_y=%e max_y=%e\n", min_y, max_y);
			}
			if (min_z < 1.0e-30) {
			printf("error!!! negative min_z MNK!\n");
			printf("min_z=%e max_z=%e\n", min_z, max_z);
			}
			*/

			integer* inum_now = new integer[maxnod];

			for (integer i = 0; i < maxnod; i++) {
				temp[i] = 0.0;
				vol[i] = 0.0;
				inum_now[i] = 0;
			}

			for (integer i = 0; i <= maxelm - 1; i++) {
				for (integer j = 0; j <= 7; j++) {
					inum_now[nvtx[j][i] - 1] += 1;
				}
			}

			for (integer i = 0; i < maxnod; i++) {
				if (inum_now[i] < 1) {
#if doubleintprecision == 1
					printf("i=%lld maxnod=%lld inum_now[%lld]=%lld\n", i, maxnod, i, inum_now[i]);
#else
					printf("i=%d maxnod=%d inum_now[%d]=%d\n", i, maxnod, i, inum_now[i]);
#endif

					system("pause");
				}
			}

			TOCHKA** pointerlist = new TOCHKA*[maxnod];
			doublereal** rthdsd_Gauss = new doublereal*[maxnod];
			for (integer i = 0; i < maxnod; i++) {
				pointerlist[i] = new TOCHKA[(inum_now[i])];
				rthdsd_Gauss[i] = new doublereal[(inum_now[i])];
			}


			for (integer i = 0; i < maxnod; i++) {
				inum_now[i] = 0;
			}

			for (integer i = 0; i <= maxelm - 1; i++) {
				// вычисление размеров текущего контрольного объёма:
				doublereal dx = 0.0, dy = 0.0, dz = 0.0;// объём текущего контроольного объёма
				volume3D(i, nvtx, pa, dx, dy, dz);

				for (integer j = 0; j <= 7; j++) {
					vol[nvtx[j][i] - 1] += dx*dy*dz;
					temp[nvtx[j][i] - 1] += dx*dy*dz*potent[i];
				}

				TOCHKA p;
				center_cord3D(i, nvtx, pa, p, 100);
				p.x = p.x + min_x;
				p.y = p.y + min_y;
				p.z = p.z + min_z;


				for (integer j = 0; j <= 7; j++) {
					pointerlist[nvtx[j][i] - 1][inum_now[nvtx[j][i] - 1]] = p;
					if (fabs(p.x) < 1.0e-30) {
						printf("problem x=%e\n", p.x);
						system("pause");
					}
					if (fabs(p.y) < 1.0e-30) {
						printf("problem y=%e\n", p.y);
						system("pause");
					}
					if (fabs(p.z) < 1.0e-30) {
						printf("problem z=%e\n", p.z);
						system("pause");
					}
					rthdsd_Gauss[nvtx[j][i] - 1][inum_now[nvtx[j][i] - 1]] = potent[i];
					inum_now[nvtx[j][i] - 1]++;
				}


				/*
				for (integer j = 0; j <= 7; j++) {
				vol[nvtx[j][i] - 1] += dx*dy*dz;
				temp[nvtx[j][i] - 1] += dx*dy*dz*potent[i];
				}
				*/
			}


			//integer jcontrol = 0;
			for (integer i = 0; i < maxnod; i++) {
				if (fabs(vol[i]) > eps_mashine) {


					doublereal** Xmatr = new doublereal*[7];
					for (integer j = 0; j <= 6; j++) {
						Xmatr[j] = new doublereal[7];
					}


					doublereal* bmatr = new doublereal[7];
					doublereal* koefmatr = new doublereal[7];

					for (integer j1 = 0; j1 <= 6; j1++) {
						for (integer j2 = 0; j2 <= 6; j2++) {
							Xmatr[j1][j2] = 0.0;
						}
						bmatr[j1] = 0.0;
						koefmatr[j1] = 0.0;
					}



					for (integer j = 0; j < inum_now[i]; j++) {

						Xmatr[0][0] += 1.0;
						Xmatr[0][1] += pointerlist[i][j].x;
						Xmatr[0][2] += pointerlist[i][j].y;
						Xmatr[0][3] += pointerlist[i][j].z;
						Xmatr[0][4] += pointerlist[i][j].x*pointerlist[i][j].x;
						Xmatr[0][5] += pointerlist[i][j].y*pointerlist[i][j].y;
						Xmatr[0][6] += pointerlist[i][j].z*pointerlist[i][j].z;

						Xmatr[1][0] += pointerlist[i][j].x;
						Xmatr[1][1] += pointerlist[i][j].x*pointerlist[i][j].x;
						Xmatr[1][2] += pointerlist[i][j].x*pointerlist[i][j].y;
						Xmatr[1][3] += pointerlist[i][j].x*pointerlist[i][j].z;
						Xmatr[1][4] += pointerlist[i][j].x*pointerlist[i][j].x*pointerlist[i][j].x;
						Xmatr[1][5] += pointerlist[i][j].x*pointerlist[i][j].y*pointerlist[i][j].y;
						Xmatr[1][6] += pointerlist[i][j].x*pointerlist[i][j].z*pointerlist[i][j].z;

						Xmatr[2][0] += pointerlist[i][j].y;
						Xmatr[2][1] += pointerlist[i][j].y*pointerlist[i][j].x;
						Xmatr[2][2] += pointerlist[i][j].y*pointerlist[i][j].y;
						Xmatr[2][3] += pointerlist[i][j].y*pointerlist[i][j].z;
						Xmatr[2][4] += pointerlist[i][j].y*pointerlist[i][j].x*pointerlist[i][j].x;
						Xmatr[2][5] += pointerlist[i][j].y*pointerlist[i][j].y*pointerlist[i][j].y;
						Xmatr[2][6] += pointerlist[i][j].y*pointerlist[i][j].z*pointerlist[i][j].z;

						Xmatr[3][0] += pointerlist[i][j].z;
						Xmatr[3][1] += pointerlist[i][j].z*pointerlist[i][j].x;
						Xmatr[3][2] += pointerlist[i][j].z*pointerlist[i][j].y;
						Xmatr[3][3] += pointerlist[i][j].z*pointerlist[i][j].z;
						Xmatr[3][4] += pointerlist[i][j].z*pointerlist[i][j].x*pointerlist[i][j].x;
						Xmatr[3][5] += pointerlist[i][j].z*pointerlist[i][j].y*pointerlist[i][j].y;
						Xmatr[3][6] += pointerlist[i][j].z*pointerlist[i][j].z*pointerlist[i][j].z;

						Xmatr[4][0] += pointerlist[i][j].x*pointerlist[i][j].x;
						Xmatr[4][1] += pointerlist[i][j].x*pointerlist[i][j].x*pointerlist[i][j].x;
						Xmatr[4][2] += pointerlist[i][j].x*pointerlist[i][j].x*pointerlist[i][j].y;
						Xmatr[4][3] += pointerlist[i][j].x*pointerlist[i][j].x*pointerlist[i][j].z;
						Xmatr[4][4] += pointerlist[i][j].x*pointerlist[i][j].x*pointerlist[i][j].x*pointerlist[i][j].x;
						Xmatr[4][5] += pointerlist[i][j].x*pointerlist[i][j].x*pointerlist[i][j].y*pointerlist[i][j].y;
						Xmatr[4][6] += pointerlist[i][j].x*pointerlist[i][j].x*pointerlist[i][j].z*pointerlist[i][j].z;

						Xmatr[5][0] += pointerlist[i][j].y*pointerlist[i][j].y;
						Xmatr[5][1] += pointerlist[i][j].y*pointerlist[i][j].y*pointerlist[i][j].x;
						Xmatr[5][2] += pointerlist[i][j].y*pointerlist[i][j].y*pointerlist[i][j].y;
						Xmatr[5][3] += pointerlist[i][j].y*pointerlist[i][j].y*pointerlist[i][j].z;
						Xmatr[5][4] += pointerlist[i][j].y*pointerlist[i][j].y*pointerlist[i][j].x*pointerlist[i][j].x;
						Xmatr[5][5] += pointerlist[i][j].y*pointerlist[i][j].y*pointerlist[i][j].y*pointerlist[i][j].y;
						Xmatr[5][6] += pointerlist[i][j].y*pointerlist[i][j].y*pointerlist[i][j].z*pointerlist[i][j].z;

						Xmatr[6][0] += pointerlist[i][j].z*pointerlist[i][j].z;
						Xmatr[6][1] += pointerlist[i][j].z*pointerlist[i][j].z*pointerlist[i][j].x;
						Xmatr[6][2] += pointerlist[i][j].z*pointerlist[i][j].z*pointerlist[i][j].y;
						Xmatr[6][3] += pointerlist[i][j].z*pointerlist[i][j].z*pointerlist[i][j].z;
						Xmatr[6][4] += pointerlist[i][j].z*pointerlist[i][j].z*pointerlist[i][j].x*pointerlist[i][j].x;
						Xmatr[6][5] += pointerlist[i][j].z*pointerlist[i][j].z*pointerlist[i][j].y*pointerlist[i][j].y;
						Xmatr[6][6] += pointerlist[i][j].z*pointerlist[i][j].z*pointerlist[i][j].z*pointerlist[i][j].z;

						bmatr[0] += rthdsd_Gauss[i][j];
						bmatr[1] += pointerlist[i][j].x*rthdsd_Gauss[i][j];
						bmatr[2] += pointerlist[i][j].y*rthdsd_Gauss[i][j];
						bmatr[3] += pointerlist[i][j].z*rthdsd_Gauss[i][j];
						bmatr[4] += pointerlist[i][j].x*pointerlist[i][j].x*rthdsd_Gauss[i][j];
						bmatr[5] += pointerlist[i][j].y*pointerlist[i][j].y*rthdsd_Gauss[i][j];
						bmatr[6] += pointerlist[i][j].z*pointerlist[i][j].z*rthdsd_Gauss[i][j];


					}

					if ((fabs(Xmatr[0][0]) < 1.e-30) || (fabs(Xmatr[1][1]) < 1.e-30) || (fabs(Xmatr[2][2]) < 1.e-30) || (fabs(Xmatr[3][3]) < 1.e-30)
						|| (fabs(Xmatr[4][4]) < 1.e-30) || (fabs(Xmatr[5][5]) < 1.e-30) || (fabs(Xmatr[6][6]) < 1.e-30)) {
#if doubleintprecision == 1
						printf("inum_now[%lld]=%lld\n", i, inum_now[i]);
#else
						printf("inum_now[%d]=%d\n", i, inum_now[i]);
#endif

						system("pause");
					}

					//Xmatr*koefmatr = bmatr;
					/*
					if (!my_version_gauss1(Xmatr, 4, bmatr, koefmatr, false, i)) {
					temp[i] = temp[i] / vol[i];
					}
					else {
					temp[i] = koefmatr[0] + koefmatr[1] * (pa[i].x+min_x) + koefmatr[2] * (pa[i].y+min_y) + koefmatr[3] * (pa[i].z+min_z);
					}
					*/
					for (integer j1 = 0; j1 <= 2000; j1++) {
						koefmatr[0] = (bmatr[0] - Xmatr[0][1] * koefmatr[1] - Xmatr[0][2] * koefmatr[2] - Xmatr[0][3] * koefmatr[3] - Xmatr[0][4] * koefmatr[4] - Xmatr[0][5] * koefmatr[5] - Xmatr[0][6] * koefmatr[6]) / Xmatr[0][0];
						koefmatr[1] = (bmatr[1] - Xmatr[1][0] * koefmatr[0] - Xmatr[1][2] * koefmatr[2] - Xmatr[1][3] * koefmatr[3] - Xmatr[1][4] * koefmatr[4] - Xmatr[1][5] * koefmatr[5] - Xmatr[1][6] * koefmatr[6]) / Xmatr[1][1];
						koefmatr[2] = (bmatr[2] - Xmatr[2][0] * koefmatr[0] - Xmatr[2][1] * koefmatr[1] - Xmatr[2][3] * koefmatr[3] - Xmatr[2][4] * koefmatr[4] - Xmatr[2][5] * koefmatr[5] - Xmatr[2][6] * koefmatr[6]) / Xmatr[2][2];
						koefmatr[3] = (bmatr[3] - Xmatr[3][0] * koefmatr[0] - Xmatr[3][1] * koefmatr[1] - Xmatr[3][2] * koefmatr[2] - Xmatr[3][4] * koefmatr[4] - Xmatr[3][5] * koefmatr[5] - Xmatr[3][6] * koefmatr[6]) / Xmatr[3][3];
						koefmatr[4] = (bmatr[4] - Xmatr[4][0] * koefmatr[0] - Xmatr[4][1] * koefmatr[1] - Xmatr[4][2] * koefmatr[2] - Xmatr[4][3] * koefmatr[3] - Xmatr[4][5] * koefmatr[5] - Xmatr[4][6] * koefmatr[6]) / Xmatr[4][4];
						koefmatr[5] = (bmatr[5] - Xmatr[5][0] * koefmatr[0] - Xmatr[5][1] * koefmatr[1] - Xmatr[5][2] * koefmatr[2] - Xmatr[5][4] * koefmatr[4] - Xmatr[5][3] * koefmatr[3] - Xmatr[5][6] * koefmatr[6]) / Xmatr[5][5];
						koefmatr[6] = (bmatr[6] - Xmatr[6][0] * koefmatr[0] - Xmatr[6][1] * koefmatr[1] - Xmatr[6][2] * koefmatr[2] - Xmatr[6][4] * koefmatr[4] - Xmatr[6][5] * koefmatr[5] - Xmatr[6][3] * koefmatr[3]) / Xmatr[6][6];
					}
					temp[i] = koefmatr[0] + koefmatr[1] * (pa[i].x + min_x) + koefmatr[2] * (pa[i].y + min_y) + koefmatr[3] * (pa[i].z + min_z) + koefmatr[4] * (pa[i].x + min_x)*(pa[i].x + min_x) + koefmatr[5] * (pa[i].y + min_y) * (pa[i].y + min_y) + koefmatr[6] * (pa[i].z + min_z) * (pa[i].z + min_z);
					heat_flux_X[i] = koefmatr[1] + 2.0*koefmatr[4] * (pa[i].x + min_x);
					heat_flux_Y[i] = koefmatr[2] + 2.0*koefmatr[5] * (pa[i].y + min_y);
					heat_flux_Z[i] = koefmatr[3] + 2.0*koefmatr[6] * (pa[i].z + min_z);


					for (integer j = 0; j <= 6; j++) {
						delete[] Xmatr[j];
					}
					delete[] Xmatr;
					delete[] bmatr;
					delete[] koefmatr;

				}
				else {
#if doubleintprecision == 1
					printf("fatal error! ANES_tecplot_export in module constr_struct_alice.cpp. i=%lld maxnod=%lld\n", i, maxnod);
					printf("vol[%lld]==%e\n", i, vol[i]);
#else
					printf("fatal error! ANES_tecplot_export in module constr_struct_alice.cpp. i=%d maxnod=%d\n", i, maxnod);
					printf("vol[%d]==%e\n", i, vol[i]);
#endif

					//getchar();
					system("PAUSE");
					temp[i] = 0.0;
				}
			}



			delete[] inum_now;
			for (integer i = 0; i < maxnod; i++) {
				delete[] pointerlist[i];
				delete[] rthdsd_Gauss[i];
			}
			delete[] pointerlist;
			delete[] rthdsd_Gauss;
		}
		else if (3 == IORDER_METHOD) {
			// линейная реконструкция, НО по расширенному шаблону.
			// Расширенный шаблон содержит больше опорных точек интерполляции, 
			// поэтому в теории реконструкция должна быть более точная.

			// Метод линейного порядка.
			doublereal min_x = 1e60;
			doublereal min_y = 1e60;
			doublereal min_z = 1e60;
			doublereal max_x = -1e60;
			doublereal max_y = -1e60;
			doublereal max_z = -1e60;

			for (integer i = 0; i < maxnod; i++) {
				if (pa[i].x < min_x) {
					min_x = pa[i].x;
				}
				if (pa[i].y < min_y) {
					min_y = pa[i].y;
				}
				if (pa[i].z < min_z) {
					min_z = pa[i].z;
				}
				if (pa[i].x > max_x) {
					max_x = pa[i].x;
				}
				if (pa[i].y > max_y) {
					max_y = pa[i].y;
				}
				if (pa[i].z > max_z) {
					max_z = pa[i].z;
				}
			}

			doublereal min_x1 = min_x;
			doublereal min_y1 = min_y;
			doublereal min_z1 = min_z;
			//min_x *= 1.2;
			//min_y *= 1.2;
			//min_z *= 1.2;

			// 05.07.2017

			min_x = 1.05*fabs(max_x - min_x);
			if (min_x < 1.0e-30) {
				min_x = 1.05*fabs(max_x);
			}
			min_y = 1.05*fabs(max_y - min_y);
			if (min_y < 1.0e-30) {
				min_y = 1.05*fabs(max_y);
			}
			min_z = 1.05*fabs(max_z - min_z);
			if (min_z < 1.0e-30) {
				min_z = 1.05*fabs(max_z);
			}

			/*
			if (min_x < 1.0e-30) {
			printf("error!!! negative min_x MNK!\n");
			printf("min_x=%e max_x=%e\n", min_x, max_x);
			}
			if (min_y < 1.0e-30) {
			printf("error!!! negative min_y MNK!\n");
			printf("min_y=%e max_y=%e\n", min_y, max_y);
			}
			if (min_z < 1.0e-30) {
			printf("error!!! negative min_z MNK!\n");
			printf("min_z=%e max_z=%e\n", min_z, max_z);
			}
			*/

			// Для ускорения сканирования в методе наименьших квадратов 8.07.2017.
			integer** q_hash = NULL;
			q_hash = new integer*[maxnod + 1];
			integer* q_ic = NULL;
			q_ic = new integer[maxnod + 1];
			for (integer j = 0; j <= maxnod; j++) {
				q_hash[j] = new integer[9];
				q_ic[j] = 0;
			}
			for (integer j = 0; j <= maxnod; j++) {
				for (integer i_1 = 0; i_1 < 9; i_1++) {
					q_hash[j][i_1] = NON_EXISTENT_NODE;
				}
			}
			for (integer i_1 = 0; i_1 <= maxelm - 1; i_1++) {
				for (integer j_1 = 0; j_1 <= 7; j_1++) {
					q_hash[nvtx[j_1][i_1]][q_ic[nvtx[j_1][i_1]]] = i_1;
					q_ic[nvtx[j_1][i_1]]++;
				}
			}

			integer* inum_now = new integer[maxnod];

			for (integer i = 0; i < maxnod; i++) {
				temp[i] = 0.0;
				vol[i] = 0.0;
				inum_now[i] = 0;
			}

			// Идея расширения шаблона: мы рассматриваем не только текущий nvtx, но и всех
			// его соседей имеющих с ним хоть одну общую вершину. Этим самым шаблон используемый
			// для реконструкции будет расширен новыми точками.
			// Модификация 30.05.2017.
			for (integer i = 0; i <= maxelm - 1; i++) {
				for (integer j = 0; j <= 7; j++) {
					inum_now[nvtx[j][i] - 1] += 1;
					//for (integer i_1 = 0; i_1 <= maxelm - 1; i_1++) {
					integer i_1 = NON_EXISTENT_NODE;
					for (integer i_2 = 0; i_2<9; i_2++) {
						i_1 = q_hash[nvtx[j][i]][i_2];
						if (i_1 >= 0) {
							for (integer j_1 = 0; j_1 <= 7; j_1++) {
								if (i_1 != i) {
									if (nvtx[j][i] == nvtx[j_1][i_1]) {
										inum_now[nvtx[j][i] - 1] += 1;
									}
								}
							}
						}
					}
				}
			}

			for (integer i = 0; i < maxnod; i++) {
				if (inum_now[i] < 1) {
#if doubleintprecision == 1
					printf("i=%lld maxnod=%lld inum_now[%lld]=%lld\n", i, maxnod, i, inum_now[i]);
#else
					printf("i=%d maxnod=%d inum_now[%d]=%d\n", i, maxnod, i, inum_now[i]);
#endif

					system("pause");
				}
			}

			TOCHKA** pointerlist = new TOCHKA*[maxnod];
			doublereal** rthdsd_Gauss = new doublereal*[maxnod];
			for (integer i = 0; i < maxnod; i++) {
				pointerlist[i] = new TOCHKA[(inum_now[i])];
				rthdsd_Gauss[i] = new doublereal[(inum_now[i])];
			}


			for (integer i = 0; i < maxnod; i++) {
				inum_now[i] = 0;
			}



			// Непосредственно само вычисление.
			for (integer i = 0; i <= maxelm - 1; i++) {
				// вычисление размеров текущего контрольного объёма:
				doublereal dx = 0.0, dy = 0.0, dz = 0.0;// объём текущего контрольного объёма
				volume3D(i, nvtx, pa, dx, dy, dz);

				for (integer j = 0; j <= 7; j++) {
					vol[nvtx[j][i] - 1] += dx*dy*dz;
					temp[nvtx[j][i] - 1] += dx*dy*dz*potent[i];
				}

				TOCHKA p;
				center_cord3D(i, nvtx, pa, p, 100);
				p.x = p.x + min_x;
				p.y = p.y + min_y;
				p.z = p.z + min_z;


				for (integer j = 0; j <= 7; j++) {
					pointerlist[nvtx[j][i] - 1][inum_now[nvtx[j][i] - 1]] = p;
					if (fabs(p.x) < 1.0e-30) {
						printf("problem x=%e\n", p.x);
						system("pause");
					}
					if (fabs(p.y) < 1.0e-30) {
						printf("problem y=%e\n", p.y);
						system("pause");
					}
					if (fabs(p.z) < 1.0e-30) {
						printf("problem z=%e\n", p.z);
						system("pause");
					}
					rthdsd_Gauss[nvtx[j][i] - 1][inum_now[nvtx[j][i] - 1]] = potent[i];
					inum_now[nvtx[j][i] - 1]++;
				}





				for (integer j = 0; j <= 7; j++) {
					//for (integer i_1 = 0; i_1 <= maxelm - 1; i_1++) {
					integer i_1 = NON_EXISTENT_NODE;
					for (integer i_2 = 0; i_2<9; i_2++) {
						i_1 = q_hash[nvtx[j][i]][i_2];
						if (i_1 >= 0) {
							for (integer j_1 = 0; j_1 <= 7; j_1++) {
								if (i_1 != i) {
									if (nvtx[j][i] == nvtx[j_1][i_1]) {
										TOCHKA p_1;
										center_cord3D(i_1, nvtx, pa, p_1, 100);
										p_1.x = p_1.x + min_x;
										p_1.y = p_1.y + min_y;
										p_1.z = p_1.z + min_z;

										pointerlist[nvtx[j][i] - 1][inum_now[nvtx[j][i] - 1]] = p_1;
										if (fabs(p_1.x) < 1.0e-30) {
											printf("problem x=%e\n", p_1.x);
											system("pause");
										}
										if (fabs(p_1.y) < 1.0e-30) {
											printf("problem y=%e\n", p_1.y);
											system("pause");
										}
										if (fabs(p_1.z) < 1.0e-30) {
											printf("problem z=%e\n", p_1.z);
											system("pause");
										}

										rthdsd_Gauss[nvtx[j][i] - 1][inum_now[nvtx[j][i] - 1]] = potent[i_1];
										inum_now[nvtx[j][i] - 1]++;
									}
								}
							}
						}
					}
				}


				/*
				for (integer j = 0; j <= 7; j++) {
				vol[nvtx[j][i] - 1] += dx*dy*dz;
				temp[nvtx[j][i] - 1] += dx*dy*dz*potent[i];
				}
				*/
			}


			for (integer j = 0; j <= maxnod; j++) {
				delete[] q_hash[j];
			}
			delete[] q_ic;
			delete[] q_hash;
			q_ic = NULL;
			q_hash = NULL;

			//integer jcontrol = 0;
			for (integer i = 0; i < maxnod; i++) {
				if (fabs(vol[i]) > eps_mashine) {


					doublereal** Xmatr = new doublereal*[4];
					for (integer j = 0; j <= 3; j++) {
						Xmatr[j] = new doublereal[4];
					}


					doublereal* bmatr = new doublereal[4];
					doublereal* koefmatr = new doublereal[4];

					for (integer j1 = 0; j1 <= 3; j1++) {
						for (integer j2 = 0; j2 <= 3; j2++) {
							Xmatr[j1][j2] = 0.0;
						}
						bmatr[j1] = 0.0;
						koefmatr[j1] = 0.0;
					}



					for (integer j = 0; j < inum_now[i]; j++) {

						Xmatr[0][0] += 1.0;
						Xmatr[0][1] += pointerlist[i][j].x;
						Xmatr[0][2] += pointerlist[i][j].y;
						Xmatr[0][3] += pointerlist[i][j].z;

						Xmatr[1][0] += pointerlist[i][j].x;
						Xmatr[1][1] += pointerlist[i][j].x*pointerlist[i][j].x;
						Xmatr[1][2] += pointerlist[i][j].x*pointerlist[i][j].y;
						Xmatr[1][3] += pointerlist[i][j].x*pointerlist[i][j].z;

						Xmatr[2][0] += pointerlist[i][j].y;
						Xmatr[2][1] += pointerlist[i][j].y*pointerlist[i][j].x;
						Xmatr[2][2] += pointerlist[i][j].y*pointerlist[i][j].y;
						Xmatr[2][3] += pointerlist[i][j].y*pointerlist[i][j].z;

						Xmatr[3][0] += pointerlist[i][j].z;
						Xmatr[3][1] += pointerlist[i][j].z*pointerlist[i][j].x;
						Xmatr[3][2] += pointerlist[i][j].z*pointerlist[i][j].y;
						Xmatr[3][3] += pointerlist[i][j].z*pointerlist[i][j].z;

						bmatr[0] += rthdsd_Gauss[i][j];
						bmatr[1] += pointerlist[i][j].x*rthdsd_Gauss[i][j];
						bmatr[2] += pointerlist[i][j].y*rthdsd_Gauss[i][j];
						bmatr[3] += pointerlist[i][j].z*rthdsd_Gauss[i][j];

					}

					if ((fabs(Xmatr[0][0]) < 1.e-30) || (fabs(Xmatr[1][1]) < 1.e-30) || (fabs(Xmatr[2][2]) < 1.e-30) || (fabs(Xmatr[3][3]) < 1.e-30)) {
#if doubleintprecision == 1
						printf("inum_now[%lld]=%lld\n", i, inum_now[i]);
#else
						printf("inum_now[%d]=%d\n", i, inum_now[i]);
#endif

						system("pause");
					}

					//Xmatr*koefmatr = bmatr;
					/*
					// Метод Гаусса не работает т.к. система линейно зависима.
					if (!my_version_gauss1(Xmatr, 4, bmatr, koefmatr, false, i)) {
					temp[i] = temp[i] / vol[i];
					}
					else {
					temp[i] = koefmatr[0] + koefmatr[1] * (pa[i].x+min_x) + koefmatr[2] * (pa[i].y+min_y) + koefmatr[3] * (pa[i].z+min_z);
					}
					*/
					for (integer j1 = 0; j1 <= 3; j1++) {
						koefmatr[j1] = 0.0;
					}
					for (integer j1 = 0; j1 <= 250; j1++) {
						doublereal alpha = 0.2;
						doublereal d_0 = koefmatr[0];
						doublereal d_1 = koefmatr[1];
						doublereal d_2 = koefmatr[2];
						doublereal d_3 = koefmatr[3];
						koefmatr[0] = (1.0 - alpha)*d_0 + alpha*((bmatr[0] - Xmatr[0][1] * koefmatr[1] - Xmatr[0][2] * koefmatr[2] - Xmatr[0][3] * koefmatr[3]) / Xmatr[0][0]);
						koefmatr[1] = (1.0 - alpha)*d_1 + alpha*((bmatr[1] - Xmatr[1][0] * koefmatr[0] - Xmatr[1][2] * koefmatr[2] - Xmatr[1][3] * koefmatr[3]) / Xmatr[1][1]);
						koefmatr[2] = (1.0 - alpha)*d_2 + alpha*((bmatr[2] - Xmatr[2][0] * koefmatr[0] - Xmatr[2][1] * koefmatr[1] - Xmatr[2][3] * koefmatr[3]) / Xmatr[2][2]);
						koefmatr[3] = (1.0 - alpha)*d_3 + alpha*((bmatr[3] - Xmatr[3][0] * koefmatr[0] - Xmatr[3][1] * koefmatr[1] - Xmatr[3][2] * koefmatr[2]) / Xmatr[3][3]);
					}
					temp[i] = koefmatr[0] + koefmatr[1] * (pa[i].x + min_x) + koefmatr[2] * (pa[i].y + min_y) + koefmatr[3] * (pa[i].z + min_z);
					//temp[i] = koefmatr[0] + koefmatr[1] * (pa[i].x ) + koefmatr[2] * (pa[i].y ) + koefmatr[3] * (pa[i].z );
					//heat_flux_X[i] = koefmatr[1];
					//heat_flux_Y[i] = koefmatr[2];
					//heat_flux_Z[i] = koefmatr[3];
					heat_flux_X[i] = 0.0;
					heat_flux_Y[i] = 0.0;
					heat_flux_Z[i] = 0.0;
					// вычисление размеров текущего контрольного объёма:
					/*
					doublereal h_1= 1.0e-4*(max_x-min_x1);
					h_1 = pow(fabs(vol[i]), 0.333);
					//h_1 = 0.5*dx1;
					heat_flux_X[i] = ((koefmatr[0] + koefmatr[1] * (pa[i].x + min_x+h_1) + koefmatr[2] * (pa[i].y + min_y) + koefmatr[3] * (pa[i].z + min_z)
					)-(koefmatr[0] + koefmatr[1] * (pa[i].x + min_x-h_1) + koefmatr[2] * (pa[i].y + min_y) + koefmatr[3] * (pa[i].z + min_z)
					)) / (2 * h_1);
					//h_1 = 1.0e-4*(max_y - min_y1);
					//h_1 = 0.5*dy1;
					heat_flux_Y[i] = ((koefmatr[0] + koefmatr[1] * (pa[i].x + min_x) + koefmatr[2] * (pa[i].y + min_y+h_1) + koefmatr[3] * (pa[i].z + min_z)
					) - (koefmatr[0] + koefmatr[1] * (pa[i].x + min_x) + koefmatr[2] * (pa[i].y + min_y-h_1) + koefmatr[3] * (pa[i].z + min_z)
					)) / (2 * h_1);
					//h_1 = 1.0e-4*(max_z - min_z1);
					//h_1 = 0.5*dz1;
					heat_flux_Z[i] = ((koefmatr[0] + koefmatr[1] * (pa[i].x + min_x) + koefmatr[2] * (pa[i].y + min_y) + koefmatr[3] * (pa[i].z + min_z+h_1)
					) - (koefmatr[0] + koefmatr[1] * (pa[i].x + min_x) + koefmatr[2] * (pa[i].y + min_y) + koefmatr[3] * (pa[i].z + min_z-h_1)
					)) / (2 * h_1);
					*/

					for (integer j = 0; j <= 3; j++) {
						delete[] Xmatr[j];
					}
					delete[] Xmatr;
					delete[] bmatr;
					delete[] koefmatr;

				}
				else {
#if doubleintprecision == 1
					printf("fatal error! ANES_tecplot_export in module constr_struct_alice.cpp. i=%lld maxnod=%lld\n", i, maxnod);
					printf("vol[%lld]==%e\n", i, vol[i]);
#else
					printf("fatal error! ANES_tecplot_export in module constr_struct_alice.cpp. i=%d maxnod=%d\n", i, maxnod);
					printf("vol[%d]==%e\n", i, vol[i]);
#endif

					//getchar();
					system("PAUSE");
					temp[i] = 0.0;
				}
			}



			delete[] inum_now;
			for (integer i = 0; i < maxnod; i++) {
				delete[] pointerlist[i];
				delete[] rthdsd_Gauss[i];
			}
			delete[] pointerlist;
			delete[] rthdsd_Gauss;


		}
		
		// запись temp
		for (integer i = 0; i < maxnod; i++) {
			//fprintf(fp_4, "%+.16f ", temp[i]);
			//if (i % 10 == 0) fprintf(fp_4, "\n");
			t_buf[i] = temp[i];
		}
		//fprintf(fp_4, "\n");

/*
		if (1) {
			// Метод нулевого порядка.

			for (integer i = 0; i < maxnod; i++) {
				vol[i] = 0.0;
			}

			for (integer i = 0; i <= maxelm - 1; i++) {
				// вычисление размеров текущего контрольного объёма:
				doublereal dx = 0.0, dy = 0.0, dz = 0.0;// объём текущего контроольного объёма
				volume3D(i, nvtx, pa, dx, dy, dz);

				for (integer j = 0; j <= 7; j++) {
					vol[nvtx[j][i] - 1] += dx*dy*dz;
					lam[nvtx[j][i] - 1] += dx*dy*dz*t.prop[LAM][i];
				}
			}
			for (integer i = 0; i < maxnod; i++) {
				if (fabs(vol[i]) > eps_mashine) {
					lam[i] = lam[i] / vol[i];
				}
				else {
#if doubleintprecision == 1
					printf("fatal error! ANES_tecplot_export in module constr_struct_alice.cpp. i=%lld maxnod=%lld\n", i, maxnod);
					printf("vol[%lld]==%e\n", i, vol[i]);
#else
					printf("fatal error! ANES_tecplot_export in module constr_struct_alice.cpp. i=%d maxnod=%d\n", i, maxnod);
					printf("vol[%d]==%e\n", i, vol[i]);
#endif

					//getchar();
					system("PAUSE");
					temp[i] = 0.0;
				}
			}
		}
		// запись теплопроводности
		for (integer i = 0; i < maxnod; i++) {
			fprintf(fp_4, "%+.16f ", lam[i]);
			if (i % 10 == 0) fprintf(fp_4, "\n");
		}
		fprintf(fp_4, "\n");

		// запись теплового потока по Х
		for (integer i = 0; i < maxnod; i++) {
			heat_flux_X[i] *= -lam[i];
			doublereal d_1 = heat_flux_X[i];
			if (d_1 > 2.0) {
				d_1 = log10(d_1);
			}
			else if (d_1<-2.0) {
				d_1 = -log10(fabs(d_1));
			}
			else {
				d_1 = 0.0;
			}
			fprintf(fp_4, "%+.16f ", d_1);
			if (i % 10 == 0) fprintf(fp_4, "\n");
		}
		fprintf(fp_4, "\n");

		// запись теплового потока по Y
		for (integer i = 0; i < maxnod; i++) {
			heat_flux_Y[i] *= -lam[i];
			doublereal d_1 = heat_flux_Y[i];
			if (d_1 > 2.0) {
				d_1 = log10(d_1);
			}
			else if (d_1<-2.0) {
				d_1 = -log10(fabs(d_1));
			}
			else {
				d_1 = 0.0;
			}
			fprintf(fp_4, "%+.16f ", d_1);
			if (i % 10 == 0) fprintf(fp_4, "\n");
		}
		fprintf(fp_4, "\n");

		// запись теплового потока по Z
		for (integer i = 0; i < maxnod; i++) {
			heat_flux_Z[i] *= -lam[i];
			doublereal d_1 = heat_flux_Z[i];
			if (d_1 > 2.0) {
				d_1 = log10(d_1);
			}
			else if (d_1<-2.0) {
				d_1 = -log10(fabs(d_1));
			}
			else {
				d_1 = 0.0;
			}
			fprintf(fp_4, "%+.16f ", d_1);
			if (i % 10 == 0) fprintf(fp_4, "\n");
		}
		fprintf(fp_4, "\n");


		// запись модуля теплового потока
		for (integer i = 0; i < maxnod; i++) {
			doublereal d_1 = sqrt(heat_flux_X[i] * heat_flux_X[i] + heat_flux_Y[i] * heat_flux_Y[i] + heat_flux_Z[i] * heat_flux_Z[i]);
			if (d_1 > 2.0) {
				d_1 = log10(d_1);
			}
			else {
				d_1 = 0.0;
			}
			fprintf(fp_4, "%+.16f ", d_1);
			if (i % 10 == 0) fprintf(fp_4, "\n");
		}
		fprintf(fp_4, "\n");
		*/

		delete[] temp;
		delete[] vol;
		delete[] lam;
		delete[] heat_flux_X;
		delete[] heat_flux_Y;
		delete[] heat_flux_Z;
		delete[] heat_flux_mag;


		for (integer i = 0; i <= maxelm - 1; i++) {
#if doubleintprecision == 1

			if (0) {
				// Контрольные объёмы состоят из кубиков. Вершины 
				// каждого такого кубика должны быть перечислены в строго порядке сортировки.
				// Сортировка : упорядочивание восьмерки по возракстанию координаты z.
				// Сортировка: упорядочивание каждой четверки по возрастанию координаты y. Всего 2 четверки.
				// Сортировка четырехз двоек по возрастанию х.


				integer invtx[8];
				doublereal znvtx[8];
				doublereal ynvtx[8];
				doublereal xnvtx[8];
				for (integer j28 = 0; j28 < 8; j28++) {
					invtx[j28] = nvtx[j28][i];
					znvtx[j28] = pa[nvtx[j28][i] - 1].z;
					ynvtx[j28] = pa[nvtx[j28][i] - 1].y;
					xnvtx[j28] = pa[nvtx[j28][i] - 1].x;
				}

				// Сортировка по возрастанию z.
				for (integer i28 = 1; i28 < 8; i28++) {
					integer inewelm = invtx[i28];
					doublereal znewelm = znvtx[i28];
					doublereal ynewelm = ynvtx[i28];
					doublereal xnewelm = xnvtx[i28];
					integer location = i28 - 1;
					while ((location >= 0) && (znvtx[location] > znewelm)) {
						// сдвигаем все элементы большие очередного.
						invtx[location + 1] = invtx[location];
						znvtx[location + 1] = znvtx[location];
						ynvtx[location + 1] = ynvtx[location];
						xnvtx[location + 1] = xnvtx[location];
						location--;
					}
					invtx[location + 1] = inewelm;
					znvtx[location + 1] = znewelm;
					ynvtx[location + 1] = ynewelm;
					xnvtx[location + 1] = xnewelm;
				}
				// Сортировка по возрастанию y. Часть 1.
				//0,1,2,3<4
				for (integer i28 = 1; i28 < 4; i28++) {
					integer inewelm = invtx[i28];
					doublereal znewelm = znvtx[i28];
					doublereal ynewelm = ynvtx[i28];
					doublereal xnewelm = xnvtx[i28];
					integer location = i28 - 1;
					while ((location >= 0) && (ynvtx[location] > ynewelm)) {
						// сдвигаем все элементы большие очередного.
						invtx[location + 1] = invtx[location];
						znvtx[location + 1] = znvtx[location];
						ynvtx[location + 1] = ynvtx[location];
						xnvtx[location + 1] = xnvtx[location];
						location--;
					}
					invtx[location + 1] = inewelm;
					znvtx[location + 1] = znewelm;
					ynvtx[location + 1] = ynewelm;
					xnvtx[location + 1] = xnewelm;
				}
				// 4,5,6,7<8
				// Сортировка по возрастанию y. Часть 2.
				for (integer i28 = 5; i28 < 8; i28++) {
					integer inewelm = invtx[i28];
					doublereal znewelm = znvtx[i28];
					doublereal ynewelm = ynvtx[i28];
					doublereal xnewelm = xnvtx[i28];
					integer location = i28 - 1;
					while ((location >= 4) && (ynvtx[location] > ynewelm)) {
						// сдвигаем все элементы большие очередного.
						invtx[location + 1] = invtx[location];
						znvtx[location + 1] = znvtx[location];
						ynvtx[location + 1] = ynvtx[location];
						xnvtx[location + 1] = xnvtx[location];
						location--;
					}
					invtx[location + 1] = inewelm;
					znvtx[location + 1] = znewelm;
					ynvtx[location + 1] = ynewelm;
					xnvtx[location + 1] = xnewelm;
				}
				//0,1<2
				//2,3<4
				//4,5<6
				//6,7<8
				// Сортировка по возрастанию X. Часть 1.
				for (integer i28 = 1; i28 < 2; i28++) {
					integer inewelm = invtx[i28];
					doublereal znewelm = znvtx[i28];
					doublereal ynewelm = ynvtx[i28];
					doublereal xnewelm = xnvtx[i28];
					integer location = i28 - 1;
					while ((location >= 0) && (xnvtx[location] > xnewelm)) {
						// сдвигаем все элементы большие очередного.
						invtx[location + 1] = invtx[location];
						znvtx[location + 1] = znvtx[location];
						ynvtx[location + 1] = ynvtx[location];
						xnvtx[location + 1] = xnvtx[location];
						location--;
					}
					invtx[location + 1] = inewelm;
					znvtx[location + 1] = znewelm;
					ynvtx[location + 1] = ynewelm;
					xnvtx[location + 1] = xnewelm;
				}
				// Сортировка по возрастанию X. Часть 2.
				for (integer i28 = 3; i28 < 4; i28++) {
					integer inewelm = invtx[i28];
					doublereal znewelm = znvtx[i28];
					doublereal ynewelm = ynvtx[i28];
					doublereal xnewelm = xnvtx[i28];
					integer location = i28 - 1;
					while ((location >= 2) && (xnvtx[location] < xnewelm)) {
						// сдвигаем все элементы большие очередного.
						invtx[location + 1] = invtx[location];
						znvtx[location + 1] = znvtx[location];
						ynvtx[location + 1] = ynvtx[location];
						xnvtx[location + 1] = xnvtx[location];
						location--;
					}
					invtx[location + 1] = inewelm;
					znvtx[location + 1] = znewelm;
					ynvtx[location + 1] = ynewelm;
					xnvtx[location + 1] = xnewelm;
				}
				// Сортировка по возрастанию X. Часть 3.
				for (integer i28 = 5; i28 < 6; i28++) {
					integer inewelm = invtx[i28];
					doublereal znewelm = znvtx[i28];
					doublereal ynewelm = ynvtx[i28];
					doublereal xnewelm = xnvtx[i28];
					integer location = i28 - 1;
					while ((location >= 4) && (xnvtx[location] > xnewelm)) {
						// сдвигаем все элементы большие очередного.
						invtx[location + 1] = invtx[location];
						znvtx[location + 1] = znvtx[location];
						ynvtx[location + 1] = ynvtx[location];
						xnvtx[location + 1] = xnvtx[location];
						location--;
					}
					invtx[location + 1] = inewelm;
					znvtx[location + 1] = znewelm;
					ynvtx[location + 1] = ynewelm;
					xnvtx[location + 1] = xnewelm;
				}
				// Сортировка по возрастанию X. Часть 4.
				for (integer i28 = 7; i28 < 8; i28++) {
					integer inewelm = invtx[i28];
					doublereal znewelm = znvtx[i28];
					doublereal ynewelm = ynvtx[i28];
					doublereal xnewelm = xnvtx[i28];
					integer location = i28 - 1;
					while ((location >= 6) && (xnvtx[location] < xnewelm)) {
						// сдвигаем все элементы большие очередного.
						invtx[location + 1] = invtx[location];
						znvtx[location + 1] = znvtx[location];
						ynvtx[location + 1] = ynvtx[location];
						xnvtx[location + 1] = xnvtx[location];
						location--;
					}
					invtx[location + 1] = inewelm;
					znvtx[location + 1] = znewelm;
					ynvtx[location + 1] = ynewelm;
					xnvtx[location + 1] = xnewelm;
				}

				//fprintf(fp_4, "%lld %lld %lld %lld %lld %lld %lld %lld \n", invtx[0], invtx[1], invtx[2], invtx[3], invtx[4], invtx[5], invtx[6], invtx[7]);
				
				nvtx_buf[0][i] = invtx[0]-1;
				nvtx_buf[1][i] = invtx[1]-1;
				nvtx_buf[2][i] = invtx[2]-1;
				nvtx_buf[3][i] = invtx[3]-1;
				nvtx_buf[4][i] = invtx[4]-1;
				nvtx_buf[5][i] = invtx[5]-1;
				nvtx_buf[6][i] = invtx[6]-1;
				nvtx_buf[7][i] = invtx[7]-1;
				
				/*
				nvtx_buf[0][i] = invtx[0];
				nvtx_buf[1][i] = invtx[1];
				nvtx_buf[2][i] = invtx[2];
				nvtx_buf[3][i] = invtx[3];
				nvtx_buf[4][i] = invtx[4];
				nvtx_buf[5][i] = invtx[5];
				nvtx_buf[6][i] = invtx[6];
				nvtx_buf[7][i] = invtx[7];
				*/

				if (nvtx[0][i] < 1) printf("bad nvtx[0][%lld]=%lld", i, nvtx[0][i]);
				if (nvtx[1][i] < 1) printf("bad nvtx[1][%lld]=%lld", i, nvtx[1][i]);
				if (nvtx[2][i] < 1) printf("bad nvtx[2][%lld]=%lld", i, nvtx[2][i]);
				if (nvtx[3][i] < 1) printf("bad nvtx[3][%lld]=%lld", i, nvtx[3][i]);
				if (nvtx[4][i] < 1) printf("bad nvtx[4][%lld]=%lld", i, nvtx[4][i]);
				if (nvtx[5][i] < 1) printf("bad nvtx[5][%lld]=%lld", i, nvtx[5][i]);
				if (nvtx[6][i] < 1) printf("bad nvtx[6][%lld]=%lld", i, nvtx[6][i]);
				if (nvtx[7][i] < 1) printf("bad nvtx[7][%lld]=%lld", i, nvtx[7][i]);
			}
			else {
				//fprintf(fp_4, "%lld %lld %lld %lld %lld %lld %lld %lld \n", nvtx[0][i], nvtx[1][i], nvtx[2][i], nvtx[3][i], nvtx[4][i], nvtx[5][i], nvtx[6][i], nvtx[7][i]);
				
				nvtx_buf[0][i] = nvtx[0][i]-1;
				nvtx_buf[1][i] = nvtx[1][i]-1;
				nvtx_buf[2][i] = nvtx[2][i]-1;
				nvtx_buf[3][i] = nvtx[3][i]-1;
				nvtx_buf[4][i] = nvtx[4][i]-1;
				nvtx_buf[5][i] = nvtx[5][i]-1;
				nvtx_buf[6][i] = nvtx[6][i]-1;
				nvtx_buf[7][i] = nvtx[7][i]-1;
			
				/*
				nvtx_buf[0][i] = nvtx[0][i];
				nvtx_buf[1][i] = nvtx[1][i];
				nvtx_buf[2][i] = nvtx[2][i];
				nvtx_buf[3][i] = nvtx[3][i];
				nvtx_buf[4][i] = nvtx[4][i];
				nvtx_buf[5][i] = nvtx[5][i];
				nvtx_buf[6][i] = nvtx[6][i];
				nvtx_buf[7][i] = nvtx[7][i];
				*/
			}
#else

			if (1) {
				integer invtx[8];
				doublereal znvtx[8];
				doublereal ynvtx[8];
				doublereal xnvtx[8];
				for (integer j28 = 0; j28 < 8; j28++) {
					invtx[j28] = nvtx[j28][i];
					znvtx[j28] = pa[nvtx[j28][i] - 1].z;
					ynvtx[j28] = pa[nvtx[j28][i] - 1].y;
					xnvtx[j28] = pa[nvtx[j28][i] - 1].x;
				}

				// Сортировка по возрастанию z.
				for (integer i28 = 1; i28 < 8; i28++) {
					integer inewelm = invtx[i28];
					doublereal znewelm = znvtx[i28];
					doublereal ynewelm = ynvtx[i28];
					doublereal xnewelm = xnvtx[i28];
					integer location = i28 - 1;
					while ((location >= 0) && (znvtx[location] > znewelm)) {
						// сдвигаем все элементы большие очередного.
						invtx[location + 1] = invtx[location];
						znvtx[location + 1] = znvtx[location];
						ynvtx[location + 1] = ynvtx[location];
						xnvtx[location + 1] = xnvtx[location];
						location--;
					}
					invtx[location + 1] = inewelm;
					znvtx[location + 1] = znewelm;
					ynvtx[location + 1] = ynewelm;
					xnvtx[location + 1] = xnewelm;
				}
				// Сортировка по возрастанию y. Часть 1.
				for (integer i28 = 1; i28 < 4; i28++) {
					integer inewelm = invtx[i28];
					doublereal znewelm = znvtx[i28];
					doublereal ynewelm = ynvtx[i28];
					doublereal xnewelm = xnvtx[i28];
					integer location = i28 - 1;
					while ((location >= 0) && (ynvtx[location] > ynewelm)) {
						// сдвигаем все элементы большие очередного.
						invtx[location + 1] = invtx[location];
						znvtx[location + 1] = znvtx[location];
						ynvtx[location + 1] = ynvtx[location];
						xnvtx[location + 1] = xnvtx[location];
						location--;
					}
					invtx[location + 1] = inewelm;
					znvtx[location + 1] = znewelm;
					ynvtx[location + 1] = ynewelm;
					xnvtx[location + 1] = xnewelm;
				}
				// Сортировка по возрастанию y. Часть 2.
				for (integer i28 = 5; i28 < 8; i28++) {
					integer inewelm = invtx[i28];
					doublereal znewelm = znvtx[i28];
					doublereal ynewelm = ynvtx[i28];
					doublereal xnewelm = xnvtx[i28];
					integer location = i28 - 1;
					while ((location >= 4) && (ynvtx[location] > ynewelm)) {
						// сдвигаем все элементы большие очередного.
						invtx[location + 1] = invtx[location];
						znvtx[location + 1] = znvtx[location];
						ynvtx[location + 1] = ynvtx[location];
						xnvtx[location + 1] = xnvtx[location];
						location--;
					}
					invtx[location + 1] = inewelm;
					znvtx[location + 1] = znewelm;
					ynvtx[location + 1] = ynewelm;
					xnvtx[location + 1] = xnewelm;
				}
				// Сортировка по возрастанию X. Часть 1.
				for (integer i28 = 1; i28 < 2; i28++) {
					integer inewelm = invtx[i28];
					doublereal znewelm = znvtx[i28];
					doublereal ynewelm = ynvtx[i28];
					doublereal xnewelm = xnvtx[i28];
					integer location = i28 - 1;
					while ((location >= 0) && (xnvtx[location] > xnewelm)) {
						// сдвигаем все элементы большие очередного.
						invtx[location + 1] = invtx[location];
						znvtx[location + 1] = znvtx[location];
						ynvtx[location + 1] = ynvtx[location];
						xnvtx[location + 1] = xnvtx[location];
						location--;
					}
					invtx[location + 1] = inewelm;
					znvtx[location + 1] = znewelm;
					ynvtx[location + 1] = ynewelm;
					xnvtx[location + 1] = xnewelm;
				}
				// Сортировка по возрастанию X. Часть 2.
				for (integer i28 = 3; i28 < 4; i28++) {
					integer inewelm = invtx[i28];
					doublereal znewelm = znvtx[i28];
					doublereal ynewelm = ynvtx[i28];
					doublereal xnewelm = xnvtx[i28];
					integer location = i28 - 1;
					while ((location >= 2) && (xnvtx[location] < xnewelm)) {
						// сдвигаем все элементы большие очередного.
						invtx[location + 1] = invtx[location];
						znvtx[location + 1] = znvtx[location];
						ynvtx[location + 1] = ynvtx[location];
						xnvtx[location + 1] = xnvtx[location];
						location--;
					}
					invtx[location + 1] = inewelm;
					znvtx[location + 1] = znewelm;
					ynvtx[location + 1] = ynewelm;
					xnvtx[location + 1] = xnewelm;
				}
				// Сортировка по возрастанию X. Часть 3.
				for (integer i28 = 5; i28 < 6; i28++) {
					integer inewelm = invtx[i28];
					doublereal znewelm = znvtx[i28];
					doublereal ynewelm = ynvtx[i28];
					doublereal xnewelm = xnvtx[i28];
					integer location = i28 - 1;
					while ((location >= 4) && (xnvtx[location] > xnewelm)) {
						// сдвигаем все элементы большие очередного.
						invtx[location + 1] = invtx[location];
						znvtx[location + 1] = znvtx[location];
						ynvtx[location + 1] = ynvtx[location];
						xnvtx[location + 1] = xnvtx[location];
						location--;
					}
					invtx[location + 1] = inewelm;
					znvtx[location + 1] = znewelm;
					ynvtx[location + 1] = ynewelm;
					xnvtx[location + 1] = xnewelm;
				}
				// Сортировка по возрастанию X. Часть 4.
				for (integer i28 = 7; i28 < 8; i28++) {
					integer inewelm = invtx[i28];
					doublereal znewelm = znvtx[i28];
					doublereal ynewelm = ynvtx[i28];
					doublereal xnewelm = xnvtx[i28];
					integer location = i28 - 1;
					while ((location >= 6) && (xnvtx[location] < xnewelm)) {
						// сдвигаем все элементы большие очередного.
						invtx[location + 1] = invtx[location];
						znvtx[location + 1] = znvtx[location];
						ynvtx[location + 1] = ynvtx[location];
						xnvtx[location + 1] = xnvtx[location];
						location--;
					}
					invtx[location + 1] = inewelm;
					znvtx[location + 1] = znewelm;
					ynvtx[location + 1] = ynewelm;
					xnvtx[location + 1] = xnewelm;
				}

				//fprintf_s(fp_4, "%d %d %d %d %d %d %d %d \n", invtx[0], invtx[1], invtx[2], invtx[3], invtx[4], invtx[5], invtx[6], invtx[7]);

				nvtx_buf[0][i] = invtx[0]-1;
				nvtx_buf[1][i] = invtx[1]-1;
				nvtx_buf[2][i] = invtx[2]-1;
				nvtx_buf[3][i] = invtx[3]-1;
				nvtx_buf[4][i] = invtx[4]-1;
				nvtx_buf[5][i] = invtx[5]-1;
				nvtx_buf[6][i] = invtx[6]-1;
				nvtx_buf[7][i] = invtx[7]-1;

				//debug
				//printf("%d %d %d %d %d %d %d %d \n", invtx[0], invtx[1], invtx[2], invtx[3], invtx[4], invtx[5], invtx[6], invtx[7]);
				//printf( "%d %d %d %d %d %d %d %d \n", nvtx[0][i], nvtx[1][i], nvtx[2][i], nvtx[3][i], nvtx[4][i], nvtx[5][i], nvtx[6][i], nvtx[7][i]);
				//getchar();

				if (nvtx[0][i] < 1) printf("bad nvtx[0][%d]=%d", i, nvtx[0][i]);
				if (nvtx[1][i] < 1) printf("bad nvtx[1][%d]=%d", i, nvtx[1][i]);
				if (nvtx[2][i] < 1) printf("bad nvtx[2][%d]=%d", i, nvtx[2][i]);
				if (nvtx[3][i] < 1) printf("bad nvtx[3][%d]=%d", i, nvtx[3][i]);
				if (nvtx[4][i] < 1) printf("bad nvtx[4][%d]=%d", i, nvtx[4][i]);
				if (nvtx[5][i] < 1) printf("bad nvtx[5][%d]=%d", i, nvtx[5][i]);
				if (nvtx[6][i] < 1) printf("bad nvtx[6][%d]=%d", i, nvtx[6][i]);
				if (nvtx[7][i] < 1) printf("bad nvtx[7][%d]=%d", i, nvtx[7][i]);
			}
			else {
				//fprintf(fp_4, "%d %d %d %d %d %d %d %d \n", nvtx[0][i], nvtx[1][i], nvtx[2][i], nvtx[3][i], nvtx[4][i], nvtx[5][i], nvtx[6][i], nvtx[7][i]);
				
				nvtx_buf[0][i] = nvtx[0][i]-1;
				nvtx_buf[1][i] = nvtx[1][i]-1;
				nvtx_buf[2][i] = nvtx[2][i]-1;
				nvtx_buf[3][i] = nvtx[3][i]-1;
				nvtx_buf[4][i] = nvtx[4][i]-1;
				nvtx_buf[5][i] = nvtx[5][i]-1;
				nvtx_buf[6][i] = nvtx[6][i]-1;
				nvtx_buf[7][i] = nvtx[7][i]-1;
			}
#endif

		}
		//fclose(fp_4);
		//WinExec("C:\\Program Files\\Tecplot\\Tecplot 360 EX 2014 R1\\bin\\tec360.exe ALICEFLOW0_24ALICEMESH.PLT", SW_NORMAL);
	//}

} // ANES_tecplot360_export_temperature_preobrazovatel


  // Специальный проверяющий корректность код.
// Он должен проходить без предупреждающих сообщений.
void ANES_ALICE_CORRECT(integer maxnod, TOCHKA* pa,
	integer maxelm, integer** nvtx) {

	// Семантика вызова:
	//ANES_ALICE_CORRECT(maxnod, pa, maxelm, nvtx);

	// 2 ноября 2016. машинное эпсилон.
	//doublereal eps_mashine = 1.0e-44; // float
	doublereal eps_mashine = 1.0e-308; // double


	// nvtx && pa сформированы, можно экспортировать в tecplot360
	

		
		doublereal* vol = new doublereal[maxnod];
		for (integer i = 0; i < maxnod; i++) {
			vol[i] = 0.0;
		}
		/*
		// Здесь мы ловим конкретную исключительную ситуацию.
		bool bfound = false;
		for (integer i = 0; i <= maxelm - 1; i++) {
		for (integer j = 0; j <= 7; j++) {
		if (nvtx[j][i] - 1 == 69462) bfound = true;
		if (nvtx[j][i] - 1 == 69463) bfound = true;
		if (nvtx[j][i] - 1 == 69464) bfound = true;
		if (nvtx[j][i] - 1 == 69477) bfound = true;
		if (nvtx[j][i] - 1 == 69478) bfound = true;
		}
		}
		if (bfound) {
		printf("bfound\n");
		}
		else
		{
		printf("notfound\n");
		}
		*/

		for (integer i = 0; i <= maxelm - 1; i++) {
			// вычисление размеров текущего контрольного объёма:
			doublereal dx = 0.0, dy = 0.0, dz = 0.0;// объём текущего контроольного объёма
			volume3D(i, nvtx, pa, dx, dy, dz);

			// Проверка корректности нахождения объёма ячейки.
#if doubleintprecision == 1
			if (dx < 1.0e-37) {
				printf("error !!! i=%lld  dx=%e\n", i, dx);
			}
			if (dy < 1.0e-37) {
				printf("error !!! i=%lld dy=%e\n", i, dy);
			}
			if (dz < 1.0e-37) {
				printf("error !!! i=%lld dz=%e\n", i, dz);
			}
#else
			if (dx < 1.0e-37) {
				printf("error !!! i=%d  dx=%e\n", i, dx);
			}
			if (dy < 1.0e-37) {
				printf("error !!! i=%d dy=%e\n", i, dy);
			}
			if (dz < 1.0e-37) {
				printf("error !!! i=%d dz=%e\n", i, dz);
			}
#endif
			

			for (integer j = 0; j <= 7; j++) {
				// nvtx[j][i]-1 - это по определению номер в pa указывающий на координаты узла.
				vol[nvtx[j][i] - 1] += dx*dy*dz;
			}
		}
		for (integer i = 0; i < maxnod; i++) {
			if (fabs(vol[i])>eps_mashine) {
				// Здесь содержится код корректной обработки.
			}
			else {
#if doubleintprecision == 1
				printf("fatal error! ANES_tecplot_export in module constr_struct_alice.cpp. i=%lld maxnod=%lld\n", i, maxnod);
				printf("vol[%lld]==%e x=%e y=%e z=%e\n", i, vol[i], pa[i].x, pa[i].y, pa[i].z);
#else
				printf("fatal error! ANES_tecplot_export in module constr_struct_alice.cpp. i=%d maxnod=%d\n", i, maxnod);
				printf("vol[%d]==%e x=%e y=%e z=%e\n", i, vol[i], pa[i].x, pa[i].y, pa[i].z);
#endif
				
				//getchar();
				system("PAUSE");
			}
		}
		
		delete[] vol;

		for (integer i = 0; i <= maxelm - 1; i++) {
#if doubleintprecision == 1
			if (nvtx[0][i] < 1) printf("bad nvtx[0][%lld]=%lld", i, nvtx[0][i]);
			if (nvtx[1][i] < 1) printf("bad nvtx[1][%lld]=%lld", i, nvtx[1][i]);
			if (nvtx[2][i] < 1) printf("bad nvtx[2][%lld]=%lld", i, nvtx[2][i]);
			if (nvtx[3][i] < 1) printf("bad nvtx[3][%lld]=%lld", i, nvtx[3][i]);
			if (nvtx[4][i] < 1) printf("bad nvtx[4][%lld]=%lld", i, nvtx[4][i]);
			if (nvtx[5][i] < 1) printf("bad nvtx[5][%lld]=%lld", i, nvtx[5][i]);
			if (nvtx[6][i] < 1) printf("bad nvtx[6][%lld]=%lld", i, nvtx[6][i]);
			if (nvtx[7][i] < 1) printf("bad nvtx[7][%lld]=%lld", i, nvtx[7][i]);
#else
			if (nvtx[0][i] < 1) printf("bad nvtx[0][%d]=%d", i, nvtx[0][i]);
			if (nvtx[1][i] < 1) printf("bad nvtx[1][%d]=%d", i, nvtx[1][i]);
			if (nvtx[2][i] < 1) printf("bad nvtx[2][%d]=%d", i, nvtx[2][i]);
			if (nvtx[3][i] < 1) printf("bad nvtx[3][%d]=%d", i, nvtx[3][i]);
			if (nvtx[4][i] < 1) printf("bad nvtx[4][%d]=%d", i, nvtx[4][i]);
			if (nvtx[5][i] < 1) printf("bad nvtx[5][%d]=%d", i, nvtx[5][i]);
			if (nvtx[6][i] < 1) printf("bad nvtx[6][%d]=%d", i, nvtx[6][i]);
			if (nvtx[7][i] < 1) printf("bad nvtx[7][%d]=%d", i, nvtx[7][i]);
#endif
			
		}
		

} // ANES_ALICE_CORRECT

// визуализация в tecplot 360 с учётом hollow блоков в программной модели.
// Построение nodes, nvtx, prop для гидродинамической подобласти. Частей 19..22 в программной модели. 
void constr_nodes_nvtx_prop_flow_alice(octTree* &oc, integer inx, integer iny, integer inz, integer &maxelm, doublereal* &xpos, doublereal* &ypos, doublereal* &zpos,
	integer iflag, BLOCK* b, integer lb,  TOCHKA* &pa, integer &maxnode, integer** &nvtx,
	doublereal** &prop, TPROP* matlist, integer* &ptr, integer* &whot_is_block, integer** &ptr_temp, integer maxelm_temp) {

	integer maxelm_loc = (inx + 1)*(iny + 1)*(inz + 1);
	// Вычисление maxelm. (maxelm_flow).
	calculate_max_elm(oc, maxelm, iflag, b, lb, false);
#if doubleintprecision == 1
	//printf("%lld \n",maxelm);
#else
	//printf("%d \n",maxelm);
#endif
	
	//getchar();


	// Выделение памяти :
	ptr_temp = NULL;
	ptr_temp = new integer*[2];
	if (ptr_temp == NULL) {
		// недостаточно памяти на данном оборудовании.
		printf("Problem : not enough memory on your equipment for ptr_temp constr struct_alice...\n");
		printf("Please any key to exit...\n");
		//getchar();
		system("pause");
		exit(1);
	}
	for (integer i = 0; i<2; i++) ptr_temp[i] = NULL;

	for (integer i = 0; i<2; i++) {
		ptr_temp[i] = new integer[maxelm_temp];
		if (ptr_temp[i] == NULL) {
			// недостаточно памяти на данном оборудовании.
#if doubleintprecision == 1
			printf("Problem : not enough memory on your equipment for ptr[%lld] constr struct_alice...\n", i);
#else
			printf("Problem : not enough memory on your equipment for ptr[%d] constr struct_alice...\n", i);
#endif
			
			printf("Please any key to exit...\n");
			//getchar();
			system("pause");
			exit(1);
		}
	}

	for (integer i = 0; i < 2; i++) {
		for (integer i87 = 0; i87 < maxelm_temp; i87++) {
			ptr_temp[i][i87] = NON_EXISTENT_NODE; // инициализация твердым телом.
		}
	}

	// ptr[1][temp_elm_id]==fluid_domain_id or -1.


	whot_is_block = NULL;
	whot_is_block = new integer[maxelm];
	if (whot_is_block == NULL) {
		// недостаточно памяти на данном оборудовании.
		printf("Problem : not enough memory on your equipment for whot_is_block constr struct_alice...\n");
		printf("Please any key to exit...\n");
		exit(1);
	}

	ptr = NULL; // maxelm_flow
	ptr = new integer[maxelm];
	if (ptr == NULL) { // -2N
		// недостаточно памяти на данном оборудовании.
		printf("Problem : not enough memory on your equipment for ptr flow0 constr struct_alice...\n");
		printf("Please any key to exit...\n");
		//getchar();
		system("pause");
		exit(1);
	}

	prop = NULL;
	switch (iflag) {
	case TEMPERATURE: 
		prop = new doublereal*[9];
		if (prop == NULL) {
			// недостаточно памяти на данном оборудовании.
			printf("Problem : not enough memory on your equipment for prop constr struct_alice...\n");
			printf("Please any key to exit...\n");
			//getchar();
			system("pause");
			exit(1);
		}
		for (integer i = 0; i<9; i++) prop[i] = NULL;
		for (integer i = 0; i<9; i++) {
			prop[i] = new doublereal[maxelm];
			if (prop[i] == NULL) {
				// недостаточно памяти на данном оборудовании.
#if doubleintprecision == 1
				printf("Problem : not enough memory on your equipment for prop[%lld] constr struct_alice...\n", i);
#else
				printf("Problem : not enough memory on your equipment for prop[%d] constr struct_alice...\n", i);
#endif
				
				printf("Please any key to exit...\n");
				//getchar();
				system("pause");
				exit(1);
			}
		}
		break;
		// Работает универсально и для гидродинамической подобласти.
	case HYDRODINAMIC: 
		prop = new doublereal*[3];
		if (prop == NULL) {
			// недостаточно памяти на данном оборудовании.
			printf("Problem : not enough memory on your equipment for prop flow constr struct_alice...\n");
			printf("Please any key to exit...\n");
			//getchar();
			system("pause");
			exit(1);
		}
		for (integer i = 0; i<3; i++) prop[i] = NULL;
		for (integer i = 0; i<3; i++) {
			prop[i] = new doublereal[maxelm];
			if (prop[i] == NULL) {
				// недостаточно памяти на данном оборудовании.
#if doubleintprecision == 1
				printf("Problem : not enough memory on your equipment for prop[%lld] flow constr struct_alice...\n", i);
#else
				printf("Problem : not enough memory on your equipment for prop[%d] flow constr struct_alice...\n", i);
#endif
				
				printf("Please any key to exit...\n");
				//getchar();
				system("pause");
				exit(1);
			}
		}
		break;
	}


	
	
	


	// Вычисление допуска.
	doublereal epsTolx = 1.0e40;
	doublereal epsToly = 1.0e40;
	doublereal epsTolz = 1.0e40;
	const doublereal mdop = 0.75;
	for (integer i = 0; i < inx; i++) {
		if (fabs(xpos[i + 1] - xpos[i]) < epsTolx) {
			epsTolx = mdop*fabs(xpos[i + 1] - xpos[i]);
		}
	}
	for (integer i = 0; i < iny; i++) {
		if (fabs(ypos[i + 1] - ypos[i]) < epsToly) {
			epsToly = mdop*fabs(ypos[i + 1] - ypos[i]);
		}
	}
	for (integer i = 0; i < inz; i++) {
		if (fabs(zpos[i + 1] - zpos[i]) < epsTolz) {
			epsTolz = mdop*fabs(zpos[i + 1] - zpos[i]);
		}
	}

	printf("geometric precision tolerance : epsTolx=%e epsToly=%e epsTolz=%e\n", epsTolx, epsToly, epsTolz);
	//system("PAUSE");


	// сформировать pa.
	// заодно посчитать количество узловых точек.
	// сформировать nvtx которые ссылаются на pa.
	// визуализировать сетку.
	TOCHKA* pa_alice = NULL;
	pa_alice = new TOCHKA[(inx + 1)*(iny + 1)*(inz + 1)];
	// Оператор new не требует проверки.
	//if (pa_alice == NULL) {
		// недостаточно памяти на данном оборудовании.
		//printf("Problem : not enough memory on your equipment for pa_alice in adaptive_local_refinement_mesh generator...\n");
		//printf("Please any key to exit...\n");
		//exit(1);
	//}
	integer marker_pa = 0;
	// И тут же сразу формируем nvtx:
	nvtx = NULL;
	nvtx = new integer*[8];
	// Оператор new не требует проверки.
	//if (nvtx == NULL) {
		// недостаточно памяти на данном оборудовании.
		//printf("Problem : not enough memory on your equipment for nvtx in adaptive_local_refinement_mesh generator...\n");
		//printf("Please any key to exit...\n");
		//exit(1);
	//}
	for (integer k_1 = 0; k_1 < 8; k_1++) {
		nvtx[k_1] = NULL;
		nvtx[k_1] = new integer[maxelm];
		// Оператор new не требует проверки.
		//if (nvtx[k_1] == NULL) {
			// недостаточно памяти на данном оборудовании.
//#if doubleintprecision == 1
	//		printf("Problem : not enough memory on your equipment for nvtx[%lld] in adaptive_local_refinement_mesh generator...\n", k_1);
//#else
	//		printf("Problem : not enough memory on your equipment for nvtx[%d] in adaptive_local_refinement_mesh generator...\n", k_1);
//#endif
			
	//		printf("Please any key to exit...\n");
		//	exit(1);
	//	}
	}
	integer imarker_nvtx = 0;

	const integer size_HASH_POLE = (inx + 1)*(iny + 1)*(inz + 1);
	HASH_POLE* hash_table_export = new HASH_POLE[size_HASH_POLE];
	for (integer i_1 = 0; i_1 < size_HASH_POLE; i_1++) {
		hash_table_export[i_1].flag = false;
		hash_table_export[i_1].inum = NON_EXISTENT_NODE;
	}

	top_ALICE_STACK = 0;
	if (oc->link0 != NULL) {
		my_ALICE_STACK[top_ALICE_STACK].link = (oc->link0);
		my_ALICE_STACK[top_ALICE_STACK].maxx = oc->link0->maxx;
		my_ALICE_STACK[top_ALICE_STACK].minx = oc->link0->minx;
		my_ALICE_STACK[top_ALICE_STACK].maxy = oc->link0->maxy;
		my_ALICE_STACK[top_ALICE_STACK].miny = oc->link0->miny;
		my_ALICE_STACK[top_ALICE_STACK].maxz = oc->link0->maxz;
		my_ALICE_STACK[top_ALICE_STACK].minz = oc->link0->minz;
		top_ALICE_STACK++;
	}
	if (oc->link1 != NULL) {
		my_ALICE_STACK[top_ALICE_STACK].link = (oc->link1);
		my_ALICE_STACK[top_ALICE_STACK].maxx = oc->link1->maxx;
		my_ALICE_STACK[top_ALICE_STACK].minx = oc->link1->minx;
		my_ALICE_STACK[top_ALICE_STACK].maxy = oc->link1->maxy;
		my_ALICE_STACK[top_ALICE_STACK].miny = oc->link1->miny;
		my_ALICE_STACK[top_ALICE_STACK].maxz = oc->link1->maxz;
		my_ALICE_STACK[top_ALICE_STACK].minz = oc->link1->minz;
		top_ALICE_STACK++;
	}
	if (oc->link2 != NULL) {
		my_ALICE_STACK[top_ALICE_STACK].link = (oc->link2);
		my_ALICE_STACK[top_ALICE_STACK].maxx = oc->link2->maxx;
		my_ALICE_STACK[top_ALICE_STACK].minx = oc->link2->minx;
		my_ALICE_STACK[top_ALICE_STACK].maxy = oc->link2->maxy;
		my_ALICE_STACK[top_ALICE_STACK].miny = oc->link2->miny;
		my_ALICE_STACK[top_ALICE_STACK].maxz = oc->link2->maxz;
		my_ALICE_STACK[top_ALICE_STACK].minz = oc->link2->minz;
		top_ALICE_STACK++;
	}
	if (oc->link3 != NULL) {
		my_ALICE_STACK[top_ALICE_STACK].link = (oc->link3);
		my_ALICE_STACK[top_ALICE_STACK].maxx = oc->link3->maxx;
		my_ALICE_STACK[top_ALICE_STACK].minx = oc->link3->minx;
		my_ALICE_STACK[top_ALICE_STACK].maxy = oc->link3->maxy;
		my_ALICE_STACK[top_ALICE_STACK].miny = oc->link3->miny;
		my_ALICE_STACK[top_ALICE_STACK].maxz = oc->link3->maxz;
		my_ALICE_STACK[top_ALICE_STACK].minz = oc->link3->minz;
		top_ALICE_STACK++;
	}
	if (oc->link4 != NULL) {
		my_ALICE_STACK[top_ALICE_STACK].link = (oc->link4);
		my_ALICE_STACK[top_ALICE_STACK].maxx = oc->link4->maxx;
		my_ALICE_STACK[top_ALICE_STACK].minx = oc->link4->minx;
		my_ALICE_STACK[top_ALICE_STACK].maxy = oc->link4->maxy;
		my_ALICE_STACK[top_ALICE_STACK].miny = oc->link4->miny;
		my_ALICE_STACK[top_ALICE_STACK].maxz = oc->link4->maxz;
		my_ALICE_STACK[top_ALICE_STACK].minz = oc->link4->minz;
		top_ALICE_STACK++;
	}
	if (oc->link5 != NULL) {
		my_ALICE_STACK[top_ALICE_STACK].link = (oc->link5);
		my_ALICE_STACK[top_ALICE_STACK].maxx = oc->link5->maxx;
		my_ALICE_STACK[top_ALICE_STACK].minx = oc->link5->minx;
		my_ALICE_STACK[top_ALICE_STACK].maxy = oc->link5->maxy;
		my_ALICE_STACK[top_ALICE_STACK].miny = oc->link5->miny;
		my_ALICE_STACK[top_ALICE_STACK].maxz = oc->link5->maxz;
		my_ALICE_STACK[top_ALICE_STACK].minz = oc->link5->minz;
		top_ALICE_STACK++;
	}
	if (oc->link6 != NULL) {
		my_ALICE_STACK[top_ALICE_STACK].link = (oc->link6);
		my_ALICE_STACK[top_ALICE_STACK].maxx = oc->link6->maxx;
		my_ALICE_STACK[top_ALICE_STACK].minx = oc->link6->minx;
		my_ALICE_STACK[top_ALICE_STACK].maxy = oc->link6->maxy;
		my_ALICE_STACK[top_ALICE_STACK].miny = oc->link6->miny;
		my_ALICE_STACK[top_ALICE_STACK].maxz = oc->link6->maxz;
		my_ALICE_STACK[top_ALICE_STACK].minz = oc->link6->minz;
		top_ALICE_STACK++;
	}
	if (oc->link7 != NULL) {
		my_ALICE_STACK[top_ALICE_STACK].link = (oc->link7);
		my_ALICE_STACK[top_ALICE_STACK].maxx = oc->link7->maxx;
		my_ALICE_STACK[top_ALICE_STACK].minx = oc->link7->minx;
		my_ALICE_STACK[top_ALICE_STACK].maxy = oc->link7->maxy;
		my_ALICE_STACK[top_ALICE_STACK].miny = oc->link7->miny;
		my_ALICE_STACK[top_ALICE_STACK].maxz = oc->link7->maxz;
		my_ALICE_STACK[top_ALICE_STACK].minz = oc->link7->minz;
		top_ALICE_STACK++;
	}


	while (top_ALICE_STACK > 0) {

		if (my_ALICE_STACK[top_ALICE_STACK - 1].link != NULL) {
			if (my_ALICE_STACK[top_ALICE_STACK - 1].link->dlist == true) {
				// Гасим информацию о посещениях.
				octTree* octree1 = my_ALICE_STACK[top_ALICE_STACK - 1].link;

				// это лист update pa.
				integer i0, i1, i2, i3, i4, i5, i6, i7;

				bool bfound = false;
				integer key_now = hash_key_alice33(inx, iny, inz, xpos, ypos, zpos, octree1->p0, epsTolx, epsToly, epsTolz);
				bfound = hash_table_export[key_now].flag;

				if (!bfound) {
					i0 = marker_pa;
					pa_alice[marker_pa] = octree1->p0;
					hash_table_export[key_now].flag = true;
					hash_table_export[key_now].inum = marker_pa;
					marker_pa++;
				}
				else {
					i0 = hash_table_export[key_now].inum;
				}
				bfound = false;
				key_now = hash_key_alice33(inx, iny, inz, xpos, ypos, zpos, octree1->p1, epsTolx, epsToly, epsTolz);
				bfound = hash_table_export[key_now].flag;
				if (!bfound) {
					i1 = marker_pa;
					pa_alice[marker_pa] = octree1->p1;
					hash_table_export[key_now].flag = true;
					hash_table_export[key_now].inum = marker_pa;
					marker_pa++;
				}
				else {
					i1 = hash_table_export[key_now].inum;
				}
				bfound = false;
				key_now = hash_key_alice33(inx, iny, inz, xpos, ypos, zpos, octree1->p2, epsTolx, epsToly, epsTolz);
				bfound = hash_table_export[key_now].flag;
				if (!bfound) {
					i2 = marker_pa;
					pa_alice[marker_pa] = octree1->p2;
					hash_table_export[key_now].flag = true;
					hash_table_export[key_now].inum = marker_pa;
					marker_pa++;
				}
				else {
					i2 = hash_table_export[key_now].inum;
				}
				bfound = false;
				key_now = hash_key_alice33(inx, iny, inz, xpos, ypos, zpos, octree1->p3, epsTolx, epsToly, epsTolz);
				bfound = hash_table_export[key_now].flag;
				if (!bfound) {
					i3 = marker_pa;
					pa_alice[marker_pa] = octree1->p3;
					hash_table_export[key_now].flag = true;
					hash_table_export[key_now].inum = marker_pa;
					marker_pa++;
				}
				else {
					i3 = hash_table_export[key_now].inum;
				}
				bfound = false;
				key_now = hash_key_alice33(inx, iny, inz, xpos, ypos, zpos, octree1->p4, epsTolx, epsToly, epsTolz);
				bfound = hash_table_export[key_now].flag;
				if (!bfound) {
					i4 = marker_pa;
					pa_alice[marker_pa] = octree1->p4;
					hash_table_export[key_now].flag = true;
					hash_table_export[key_now].inum = marker_pa;
					marker_pa++;
				}
				else {
					i4 = hash_table_export[key_now].inum;
				}
				bfound = false;
				key_now = hash_key_alice33(inx, iny, inz, xpos, ypos, zpos, octree1->p5, epsTolx, epsToly, epsTolz);
				bfound = hash_table_export[key_now].flag;
				if (!bfound) {
					i5 = marker_pa;
					pa_alice[marker_pa] = octree1->p5;
					hash_table_export[key_now].flag = true;
					hash_table_export[key_now].inum = marker_pa;
					marker_pa++;
				}
				else {
					i5 = hash_table_export[key_now].inum;
				}
				bfound = false;
				key_now = hash_key_alice33(inx, iny, inz, xpos, ypos, zpos, octree1->p6, epsTolx, epsToly, epsTolz);
				bfound = hash_table_export[key_now].flag;
				if (!bfound) {
					i6 = marker_pa;
					pa_alice[marker_pa] = octree1->p6;
					hash_table_export[key_now].flag = true;
					hash_table_export[key_now].inum = marker_pa;
					marker_pa++;
				}
				else {
					i6 = hash_table_export[key_now].inum;
				}
				bfound = false;
				key_now = hash_key_alice33(inx, iny, inz, xpos, ypos, zpos, octree1->p7, epsTolx, epsToly, epsTolz);
				bfound = hash_table_export[key_now].flag;
				if (!bfound) {
					i7 = marker_pa;
					pa_alice[marker_pa] = octree1->p7;
					hash_table_export[key_now].flag = true;
					hash_table_export[key_now].inum = marker_pa;
					marker_pa++;
				}
				else {
					i7 = hash_table_export[key_now].inum;
				}
				TOCHKA p;
				p.x = 0.125*(octree1->p0.x + octree1->p1.x + octree1->p2.x + octree1->p3.x + octree1->p4.x + octree1->p5.x + octree1->p6.x + octree1->p7.x);
				p.y = 0.125*(octree1->p0.y + octree1->p1.y + octree1->p2.y + octree1->p3.y + octree1->p4.y + octree1->p5.y + octree1->p6.y + octree1->p7.y);
				p.z = 0.125*(octree1->p0.z + octree1->p1.z + octree1->p2.z + octree1->p3.z + octree1->p4.z + octree1->p5.z + octree1->p6.z + octree1->p7.z);
				integer ib;
				bool inDomain = false;
				switch (iflag) {
				case TEMPERATURE: inDomain = in_model_temp(p, ib, b, lb);
					break;
					// Работает универсально и для гидродинамической подобласти.
				case HYDRODINAMIC: inDomain = in_model_flow(p, ib, b, lb);
					break;
				}
				if ((p.x < b[0].g.xS) || (p.x > b[0].g.xE) || (p.y < b[0].g.yS) || (p.y > b[0].g.yE) || (p.z < b[0].g.zS) || (p.z > b[0].g.zE)) {
					// Лежит за пределами кабинета.
					inDomain = false;
				}
				if (inDomain) {
					integer l = imarker_nvtx;
					switch (iflag) {
					case TEMPERATURE: prop[RHO][l] = matlist[b[ib].imatid].rho;
						//prop[HEAT_CAPACITY][l] = matlist[b[ib].imatid].cp;
						//prop[LAM][l] = matlist[b[ib].imatid].lam;
						prop[HEAT_CAPACITY][l] = get_lam(matlist[b[ib].imatid].n_cp, matlist[b[ib].imatid].temp_cp, matlist[b[ib].imatid].arr_cp, 25.0);
						prop[LAM][l] = get_lam(matlist[b[ib].imatid].n_lam, matlist[b[ib].imatid].temp_lam, matlist[b[ib].imatid].arr_lam, 25.0);

						prop[MULT_LAM_X][l] = matlist[b[ib].imatid].orthotropy_multiplyer_x;
						prop[MULT_LAM_Y][l] = matlist[b[ib].imatid].orthotropy_multiplyer_y;
						prop[MULT_LAM_Z][l] = matlist[b[ib].imatid].orthotropy_multiplyer_z;
						// Эти объекты при гидродинамической обработке не передаются к заполнению.
						//Sc[l] = b[ib].Sc;
						//ipower_time_depend[l] = b[ib].ipower_time_depend;
						break;
					case HYDRODINAMIC:prop[RHO][l] = matlist[b[ib].imatid].rho;
						prop[MU][l] = matlist[b[ib].imatid].mu;
						prop[BETA_T][l] = matlist[b[ib].imatid].beta_t;
						break;
					}

					// Связь гидродинамики с теплопроводностью
					// для задач сопряжённого теплообмена.
					// ptr[fluid_elm_id]=соответствующий temper_elm_id.
					ptr[imarker_nvtx] = octree1->inum_TD-1;
					// Связь теплового внутреннего элемента с гидродинамическим внутреним элементом или -1.
					ptr_temp[0][octree1->inum_TD - 1] = imarker_nvtx; // 25.09.2016.
					// Теперь у нас все гидродинамические подобласти слиты в одну и это позволяет считать 
					// в одной модели сразу несколько несвязанных гидродинамически fluid областей которые 
					// связаны лишь уравнением теплопередачи. Поэтому заполняем следующую структуру : 
					ptr_temp[1][octree1->inum_TD - 1] = 0; // Одна единственная fluid зона с идентификатором 0.
					octree1->inum_FD = imarker_nvtx + 1; // Нумерация начинается с 1.

					// Нумерация начинается с единицы .
					nvtx[0][imarker_nvtx] = i0 + 1;
					nvtx[1][imarker_nvtx] = i1 + 1;
					nvtx[2][imarker_nvtx] = i2 + 1;
					nvtx[3][imarker_nvtx] = i3 + 1;
					nvtx[4][imarker_nvtx] = i4 + 1;
					nvtx[5][imarker_nvtx] = i5 + 1;
					nvtx[6][imarker_nvtx] = i6 + 1;
					nvtx[7][imarker_nvtx] = i7 + 1;

					whot_is_block[imarker_nvtx] = ib;

					imarker_nvtx++;
				}

				octree1 = NULL;
				my_ALICE_STACK[top_ALICE_STACK - 1].link = NULL;
				top_ALICE_STACK--;


			}
			else {

				// продолжаем добираться до листьев.
				STACK_ALICE buf1 = my_ALICE_STACK[top_ALICE_STACK - 1];
				STACK_ALICE* buf = &buf1;
				top_ALICE_STACK--;
				if (buf->link->link0 != NULL) {
					my_ALICE_STACK[top_ALICE_STACK].link = (buf->link->link0);
					my_ALICE_STACK[top_ALICE_STACK].maxx = buf->link->link0->maxx;
					my_ALICE_STACK[top_ALICE_STACK].minx = buf->link->link0->minx;
					my_ALICE_STACK[top_ALICE_STACK].maxy = buf->link->link0->maxy;
					my_ALICE_STACK[top_ALICE_STACK].miny = buf->link->link0->miny;
					my_ALICE_STACK[top_ALICE_STACK].maxz = buf->link->link0->maxz;
					my_ALICE_STACK[top_ALICE_STACK].minz = buf->link->link0->minz;
					top_ALICE_STACK++;
				}
				if (buf->link->link1 != NULL) {
					my_ALICE_STACK[top_ALICE_STACK].link = (buf->link->link1);
					my_ALICE_STACK[top_ALICE_STACK].maxx = buf->link->link1->maxx;
					my_ALICE_STACK[top_ALICE_STACK].minx = buf->link->link1->minx;
					my_ALICE_STACK[top_ALICE_STACK].maxy = buf->link->link1->maxy;
					my_ALICE_STACK[top_ALICE_STACK].miny = buf->link->link1->miny;
					my_ALICE_STACK[top_ALICE_STACK].maxz = buf->link->link1->maxz;
					my_ALICE_STACK[top_ALICE_STACK].minz = buf->link->link1->minz;
					top_ALICE_STACK++;
				}
				if (buf->link->link2 != NULL) {
					my_ALICE_STACK[top_ALICE_STACK].link = (buf->link->link2);
					my_ALICE_STACK[top_ALICE_STACK].maxx = buf->link->link2->maxx;
					my_ALICE_STACK[top_ALICE_STACK].minx = buf->link->link2->minx;
					my_ALICE_STACK[top_ALICE_STACK].maxy = buf->link->link2->maxy;
					my_ALICE_STACK[top_ALICE_STACK].miny = buf->link->link2->miny;
					my_ALICE_STACK[top_ALICE_STACK].maxz = buf->link->link2->maxz;
					my_ALICE_STACK[top_ALICE_STACK].minz = buf->link->link2->minz;
					top_ALICE_STACK++;
				}
				if (buf->link->link3 != NULL) {
					my_ALICE_STACK[top_ALICE_STACK].link = (buf->link->link3);
					my_ALICE_STACK[top_ALICE_STACK].maxx = buf->link->link3->maxx;
					my_ALICE_STACK[top_ALICE_STACK].minx = buf->link->link3->minx;
					my_ALICE_STACK[top_ALICE_STACK].maxy = buf->link->link3->maxy;
					my_ALICE_STACK[top_ALICE_STACK].miny = buf->link->link3->miny;
					my_ALICE_STACK[top_ALICE_STACK].maxz = buf->link->link3->maxz;
					my_ALICE_STACK[top_ALICE_STACK].minz = buf->link->link3->minz;
					top_ALICE_STACK++;
				}
				if (buf->link->link4 != NULL) {
					my_ALICE_STACK[top_ALICE_STACK].link = (buf->link->link4);
					my_ALICE_STACK[top_ALICE_STACK].maxx = buf->link->link4->maxx;
					my_ALICE_STACK[top_ALICE_STACK].minx = buf->link->link4->minx;
					my_ALICE_STACK[top_ALICE_STACK].maxy = buf->link->link4->maxy;
					my_ALICE_STACK[top_ALICE_STACK].miny = buf->link->link4->miny;
					my_ALICE_STACK[top_ALICE_STACK].maxz = buf->link->link4->maxz;
					my_ALICE_STACK[top_ALICE_STACK].minz = buf->link->link4->minz;
					top_ALICE_STACK++;
				}
				if (buf->link->link5 != NULL) {
					my_ALICE_STACK[top_ALICE_STACK].link = (buf->link->link5);
					my_ALICE_STACK[top_ALICE_STACK].maxx = buf->link->link5->maxx;
					my_ALICE_STACK[top_ALICE_STACK].minx = buf->link->link5->minx;
					my_ALICE_STACK[top_ALICE_STACK].maxy = buf->link->link5->maxy;
					my_ALICE_STACK[top_ALICE_STACK].miny = buf->link->link5->miny;
					my_ALICE_STACK[top_ALICE_STACK].maxz = buf->link->link5->maxz;
					my_ALICE_STACK[top_ALICE_STACK].minz = buf->link->link5->minz;
					top_ALICE_STACK++;
				}
				if (buf->link->link6 != NULL) {
					my_ALICE_STACK[top_ALICE_STACK].link = (buf->link->link6);
					my_ALICE_STACK[top_ALICE_STACK].maxx = buf->link->link6->maxx;
					my_ALICE_STACK[top_ALICE_STACK].minx = buf->link->link6->minx;
					my_ALICE_STACK[top_ALICE_STACK].maxy = buf->link->link6->maxy;
					my_ALICE_STACK[top_ALICE_STACK].miny = buf->link->link6->miny;
					my_ALICE_STACK[top_ALICE_STACK].maxz = buf->link->link6->maxz;
					my_ALICE_STACK[top_ALICE_STACK].minz = buf->link->link6->minz;
					top_ALICE_STACK++;
				}
				if (buf->link->link7 != NULL) {
					my_ALICE_STACK[top_ALICE_STACK].link = (buf->link->link7);
					my_ALICE_STACK[top_ALICE_STACK].maxx = buf->link->link7->maxx;
					my_ALICE_STACK[top_ALICE_STACK].minx = buf->link->link7->minx;
					my_ALICE_STACK[top_ALICE_STACK].maxy = buf->link->link7->maxy;
					my_ALICE_STACK[top_ALICE_STACK].miny = buf->link->link7->miny;
					my_ALICE_STACK[top_ALICE_STACK].maxz = buf->link->link7->maxz;
					my_ALICE_STACK[top_ALICE_STACK].minz = buf->link->link7->minz;
					top_ALICE_STACK++;
				}

			}
		}

	}
	delete[] hash_table_export;

	maxnode = marker_pa;
	pa = NULL;
	pa = new TOCHKA[maxnode];
	if (pa == NULL) {
		// недостаточно памяти на данном оборудовании.
		printf("Problem : not enough memory on your equipment for pa constr struct alice...\n");
		printf("Please any key to exit...\n");
		//getchar();
		system("pause");
		exit(1);
	}
	for (integer i = 0; i < maxnode; i++) {
		pa[i] = pa_alice[i];
	}
	delete[] pa_alice;
	pa_alice = NULL;

	// nvtx && pa сформированы, можно экспортировать в tecplot360
	FILE *fp_4 = NULL;
	errno_t err_4=0;
#ifdef MINGW_COMPILLER
	fp_4= fopen64("ALICEFLOW0_24ALICEMESH_FLOW.PLT", "w");
#else
	err_4 = fopen_s(&fp_4, "ALICEFLOW0_24ALICEMESH_FLOW.PLT", "w");
#endif
	if (err_4 != 0) {
		printf("Create File temp Error\n");
		//getchar();
		system("pause");

	}
	else {
		if (fp_4 != NULL) {
			fprintf(fp_4, "TITLE = \"ALICEFLOW0_24\"\n");
			fprintf(fp_4, "VARIABLES = x, y, z\n");
#if doubleintprecision == 1
			fprintf(fp_4, "ZONE T=\"Rampant\", N=%lld, E=%lld, ET=BRICK, F=FEBLOCK\n\n", marker_pa, imarker_nvtx);
#else
			fprintf(fp_4, "ZONE T=\"Rampant\", N=%d, E=%d, ET=BRICK, F=FEBLOCK\n\n", marker_pa, imarker_nvtx);
#endif
			
			// запись x
			for (integer i = 0; i < marker_pa; i++) {
				fprintf(fp_4, "%+.16f ", pa[i].x);
				if (i % 10 == 0) fprintf(fp_4, "\n");
			}
			fprintf(fp_4, "\n");
			// запись y
			for (integer i = 0; i < marker_pa; i++) {
				fprintf(fp_4, "%+.16f ", pa[i].y);
				if (i % 10 == 0) fprintf(fp_4, "\n");
			}
			fprintf(fp_4, "\n");
			// запись z
			for (integer i = 0; i < marker_pa; i++) {
				fprintf(fp_4, "%+.16f ", pa[i].z);
				if (i % 10 == 0) fprintf(fp_4, "\n");
			}
			fprintf(fp_4, "\n");
			for (integer i = 0; i <= imarker_nvtx - 1; i++) {
#if doubleintprecision == 1
				fprintf(fp_4, "%lld %lld %lld %lld %lld %lld %lld %lld \n", nvtx[0][i], nvtx[1][i], nvtx[2][i], nvtx[3][i], nvtx[4][i], nvtx[5][i], nvtx[6][i], nvtx[7][i]);
#else
				fprintf(fp_4, "%d %d %d %d %d %d %d %d \n", nvtx[0][i], nvtx[1][i], nvtx[2][i], nvtx[3][i], nvtx[4][i], nvtx[5][i], nvtx[6][i], nvtx[7][i]);
#endif
			}
			fclose(fp_4);
			//WinExec("C:\\Program Files\\Tecplot\\Tecplot 360 EX 2014 R1\\bin\\tec360.exe ALICEFLOW0_24ALICEMESH_FLOW.PLT", SW_NORMAL);
		}
	}
} // constr_nodes_nvtx_prop_flow_alice

// Увеличивает количество граничных граней при нахождении внутреннего источника тепла.
void patch_maxbound(integer iplane, SOURCE* s, integer ls, doublereal x_c, doublereal y_c, doublereal z_c, integer &maxbound, bool &binc, integer &lsid) {
	bool bfind = false;
	binc = false;
	for (integer j = 0; j<ls; j++) {
		if (s[j].iPlane == iplane) {
			switch (iplane) {
			case XY: s[j].g.zE = s[j].g.zS;
				if ((x_c>s[j].g.xS) && (x_c<s[j].g.xE) && (y_c>s[j].g.yS) && (y_c<s[j].g.yE) && (fabs(z_c - s[j].g.zE)<admission)) bfind = true;
				break;
			case XZ:  s[j].g.yE = s[j].g.yS;
				if ((x_c>s[j].g.xS) && (x_c<s[j].g.xE) && (z_c>s[j].g.zS) && (z_c<s[j].g.zE) && (fabs(y_c - s[j].g.yE)<admission)) bfind = true;
				break;
			case YZ: s[j].g.xE = s[j].g.xS;
				if ((z_c>s[j].g.zS) && (z_c<s[j].g.zE) && (y_c>s[j].g.yS) && (y_c<s[j].g.yE) && (fabs(x_c - s[j].g.xE)<admission)) bfind = true;
				break;
			}
		}
		if (bfind) {
			lsid = j;
			break; // досрочный выход из цикла for.
		}
	}
	if (bfind) {
		// нужно присвоить грани соответствующий номер maxbound.
		// gran[G][i]=maxbound++;
		maxbound++;
		binc = true;
	}
} // patch_maxbound

// Вычисляет количество граничных элементов модели для температуры.
// Внутренние источники тепла также пронумерованы.
void calculate_max_bound_temp(octTree* &oc, integer &maxbound, integer maxelm_memo, BLOCK* b, integer lb, 
	SOURCE* s, integer ls) {
	maxbound = 0; // инициализация.
	integer maxelm = 0;
	bool *bvisit = NULL;
	bvisit = new bool[maxelm_memo];
	// оператор new не требует проверки на null.
	//if (bvisit == NULL) {
		// недостаточно памяти на данном оборудовании.
		//printf("Problem : not enough memory on your equipment for bvisit constr struct_alice...\n");
		//printf("Please any key to exit...\n");
		//getchar();
		//system("pause");
		//exit(1);
	//}
	for (integer j = 0; j<maxelm_memo; j++) bvisit[j] = false; // признак посещения узла.
	doublereal x_c, y_c, z_c;
	integer iplane;

	top_ALICE_STACK = 0;	
	if (oc->link0 != NULL) {
		my_ALICE_STACK[top_ALICE_STACK].link = (oc->link0);
		my_ALICE_STACK[top_ALICE_STACK].maxx = oc->link0->maxx;
		my_ALICE_STACK[top_ALICE_STACK].minx = oc->link0->minx;
		my_ALICE_STACK[top_ALICE_STACK].maxy = oc->link0->maxy;
		my_ALICE_STACK[top_ALICE_STACK].miny = oc->link0->miny;
		my_ALICE_STACK[top_ALICE_STACK].maxz = oc->link0->maxz;
		my_ALICE_STACK[top_ALICE_STACK].minz = oc->link0->minz;
		top_ALICE_STACK++;
	}
	if (oc->link1 != NULL) {
		my_ALICE_STACK[top_ALICE_STACK].link = (oc->link1);
		my_ALICE_STACK[top_ALICE_STACK].maxx = oc->link1->maxx;
		my_ALICE_STACK[top_ALICE_STACK].minx = oc->link1->minx;
		my_ALICE_STACK[top_ALICE_STACK].maxy = oc->link1->maxy;
		my_ALICE_STACK[top_ALICE_STACK].miny = oc->link1->miny;
		my_ALICE_STACK[top_ALICE_STACK].maxz = oc->link1->maxz;
		my_ALICE_STACK[top_ALICE_STACK].minz = oc->link1->minz;
		top_ALICE_STACK++;
	}
	if (oc->link2 != NULL) {
		my_ALICE_STACK[top_ALICE_STACK].link = (oc->link2);
		my_ALICE_STACK[top_ALICE_STACK].maxx = oc->link2->maxx;
		my_ALICE_STACK[top_ALICE_STACK].minx = oc->link2->minx;
		my_ALICE_STACK[top_ALICE_STACK].maxy = oc->link2->maxy;
		my_ALICE_STACK[top_ALICE_STACK].miny = oc->link2->miny;
		my_ALICE_STACK[top_ALICE_STACK].maxz = oc->link2->maxz;
		my_ALICE_STACK[top_ALICE_STACK].minz = oc->link2->minz;
		top_ALICE_STACK++;
	}
	if (oc->link3 != NULL) {
		my_ALICE_STACK[top_ALICE_STACK].link = (oc->link3);
		my_ALICE_STACK[top_ALICE_STACK].maxx = oc->link3->maxx;
		my_ALICE_STACK[top_ALICE_STACK].minx = oc->link3->minx;
		my_ALICE_STACK[top_ALICE_STACK].maxy = oc->link3->maxy;
		my_ALICE_STACK[top_ALICE_STACK].miny = oc->link3->miny;
		my_ALICE_STACK[top_ALICE_STACK].maxz = oc->link3->maxz;
		my_ALICE_STACK[top_ALICE_STACK].minz = oc->link3->minz;
		top_ALICE_STACK++;
	}
	if (oc->link4 != NULL) {
		my_ALICE_STACK[top_ALICE_STACK].link = (oc->link4);
		my_ALICE_STACK[top_ALICE_STACK].maxx = oc->link4->maxx;
		my_ALICE_STACK[top_ALICE_STACK].minx = oc->link4->minx;
		my_ALICE_STACK[top_ALICE_STACK].maxy = oc->link4->maxy;
		my_ALICE_STACK[top_ALICE_STACK].miny = oc->link4->miny;
		my_ALICE_STACK[top_ALICE_STACK].maxz = oc->link4->maxz;
		my_ALICE_STACK[top_ALICE_STACK].minz = oc->link4->minz;
		top_ALICE_STACK++;
	}
	if (oc->link5 != NULL) {
		my_ALICE_STACK[top_ALICE_STACK].link = (oc->link5);
		my_ALICE_STACK[top_ALICE_STACK].maxx = oc->link5->maxx;
		my_ALICE_STACK[top_ALICE_STACK].minx = oc->link5->minx;
		my_ALICE_STACK[top_ALICE_STACK].maxy = oc->link5->maxy;
		my_ALICE_STACK[top_ALICE_STACK].miny = oc->link5->miny;
		my_ALICE_STACK[top_ALICE_STACK].maxz = oc->link5->maxz;
		my_ALICE_STACK[top_ALICE_STACK].minz = oc->link5->minz;
		top_ALICE_STACK++;
	}
	if (oc->link6 != NULL) {
		my_ALICE_STACK[top_ALICE_STACK].link = (oc->link6);
		my_ALICE_STACK[top_ALICE_STACK].maxx = oc->link6->maxx;
		my_ALICE_STACK[top_ALICE_STACK].minx = oc->link6->minx;
		my_ALICE_STACK[top_ALICE_STACK].maxy = oc->link6->maxy;
		my_ALICE_STACK[top_ALICE_STACK].miny = oc->link6->miny;
		my_ALICE_STACK[top_ALICE_STACK].maxz = oc->link6->maxz;
		my_ALICE_STACK[top_ALICE_STACK].minz = oc->link6->minz;
		top_ALICE_STACK++;
	}
	if (oc->link7 != NULL) {
		my_ALICE_STACK[top_ALICE_STACK].link = (oc->link7);
		my_ALICE_STACK[top_ALICE_STACK].maxx = oc->link7->maxx;
		my_ALICE_STACK[top_ALICE_STACK].minx = oc->link7->minx;
		my_ALICE_STACK[top_ALICE_STACK].maxy = oc->link7->maxy;
		my_ALICE_STACK[top_ALICE_STACK].miny = oc->link7->miny;
		my_ALICE_STACK[top_ALICE_STACK].maxz = oc->link7->maxz;
		my_ALICE_STACK[top_ALICE_STACK].minz = oc->link7->minz;
		top_ALICE_STACK++;
	}
	while (top_ALICE_STACK > 0) {
		if (my_ALICE_STACK[top_ALICE_STACK - 1].link != NULL) {
			if (my_ALICE_STACK[top_ALICE_STACK - 1].link->dlist == true) {
				// Гасим информацию о посещениях.
				octTree* octree1 = my_ALICE_STACK[top_ALICE_STACK - 1].link;

				if (!((octree1->link0 == NULL) && (octree1->link1 == NULL) && (octree1->link2 == NULL) && (octree1->link3 == NULL) && (octree1->link4 == NULL) && (octree1->link5 == NULL) && (octree1->link6 == NULL) && (octree1->link7 == NULL))) {
					printf("eto ne list error!!!\n");
					//getchar();
					system("PAUSE");
				}

				TOCHKA p;
				p.x = 0.125*(octree1->p0.x + octree1->p1.x + octree1->p2.x + octree1->p3.x + octree1->p4.x + octree1->p5.x + octree1->p6.x + octree1->p7.x);
				p.y = 0.125*(octree1->p0.y + octree1->p1.y + octree1->p2.y + octree1->p3.y + octree1->p4.y + octree1->p5.y + octree1->p6.y + octree1->p7.y);
				p.z = 0.125*(octree1->p0.z + octree1->p1.z + octree1->p2.z + octree1->p3.z + octree1->p4.z + octree1->p5.z + octree1->p6.z + octree1->p7.z);
				integer ib;
				bool inDomain = false;
				
				inDomain = in_model_temp(p, ib, b, lb); // TEMPERATURE
				
				if ((p.x < b[0].g.xS) || (p.x > b[0].g.xE) || (p.y < b[0].g.yS) || (p.y > b[0].g.yE) || (p.z < b[0].g.zS) || (p.z > b[0].g.zE)) {
					// Лежит за пределами кабинета.
					inDomain = false;
				}
				if (inDomain) {
					// Узел принадлекжит расчётной области и его номер maxelm==octree1->inum_TD.

					// сосед не существует
					// записываем идентификатор грани и
					// увеличиваем счётчик граней.
					if (!octree1->b4E) {
						if (octree1->linkE == NULL) {
							maxbound++;
						}
						else {
							if (octree1->linkE->inum_TD == 0) {
								maxbound++;
							}
							else {
								// сосед существует.
								if (bvisit[octree1->linkE->inum_TD - 1]) {
									// узел уже был посещён
									// Здесь неким образом модифицировалось gran_t
								}
								else {
									bool binc;
									integer lsid;
									// узнать координаты центра грани и ориентацию в пространстве
									x_c = octree1->p1.x;  
									y_c = 0.5*(octree1->p1.y+octree1->p2.y);
									z_c = 0.5*(octree1->p1.z+octree1->p5.z);
									iplane = YZ;
									patch_maxbound(iplane, s, ls, x_c, y_c, z_c, maxbound,binc,lsid);
								}
							}
						}
					}
					else {

						if ((octree1->linkE1 == NULL) && (octree1->linkE2 == NULL) && (octree1->linkE5 == NULL) && (octree1->linkE6 == NULL)) {
							// Зачем дробить если граница заканчивается NULL.
							maxbound++; // Вместо четырёх один.
						}
						else if ((octree1->linkE1 != NULL) && (octree1->linkE2 != NULL) && (octree1->linkE5 != NULL) && (octree1->linkE6 != NULL)&&
							(octree1->linkE1->inum_TD == 0) && (octree1->linkE2->inum_TD == 0) && (octree1->linkE5->inum_TD == 0) && (octree1->linkE6->inum_TD == 0)) {
							// Зачем дробить если граница заканчивается NULL.
							maxbound++; // Вместо четырёх один.
						}
						else if ((octree1->linkE1 != NULL) && (octree1->linkE2 == NULL) && (octree1->linkE5 != NULL) && (octree1->linkE6 == NULL) &&
							(octree1->linkE1->inum_TD == 0) && (octree1->linkE5->inum_TD == 0) ) {
							// Зачем дробить если граница заканчивается NULL.
							maxbound++; // Вместо двух один.
						}
						else if ((octree1->linkE1 != NULL) && (octree1->linkE2 != NULL) && (octree1->linkE5 == NULL) && (octree1->linkE6 == NULL) &&
							(octree1->linkE1->inum_TD == 0) && (octree1->linkE2->inum_TD == 0)) {
							// Зачем дробить если граница заканчивается NULL.
							maxbound++; // Вместо двух один.
						}
						else if ((octree1->linkE1 != NULL) && (octree1->linkE2 == NULL) && (octree1->linkE5 == NULL) && (octree1->linkE6 == NULL) &&
							(octree1->linkE1->inum_TD == 0)) {
							// Зачем дробить если граница заканчивается NULL.
							maxbound++; // Вместо двух один.
						}
						else {

							bool binc1, binc2, binc3, binc4;
							integer lsid1, lsid2, lsid3, lsid4;

							if (octree1->linkE1 == NULL) {
								//maxbound++;
								printf("Error!!! E1 is NULL\n");
								//getchar();
								system("PAUSE");
							}
							else {
								if (octree1->linkE1->inum_TD == 0) {
									maxbound++;
									
								}
								else {
									// сосед существует.
									if (bvisit[octree1->linkE1->inum_TD - 1]) {
										// узел уже был посещён
										// Здесь неким образом модифицировалось gran_t
									}
									else {
										// узнать координаты центра грани и ориентацию в пространстве
										x_c = octree1->p1.x;
										y_c = 0.5*(octree1->linkE1->p0.y + octree1->linkE1->p3.y);
										z_c = 0.5*(octree1->linkE1->p0.z + octree1->linkE1->p4.z);
										iplane = YZ;
										patch_maxbound(iplane, s, ls, x_c, y_c, z_c, maxbound, binc1, lsid1);
										
									}
								}
							}

							if (octree1->linkE2 == NULL) {
								//maxbound++;
								// 21 сентября 2016 (вырождение ячейки).
							}
							else {
								if (octree1->linkE2->inum_TD == 0) {
									maxbound++;
									
								}
								else {
									// сосед существует.
									if (bvisit[octree1->linkE2->inum_TD - 1]) {
										// узел уже был посещён
										// Здесь неким образом модифицировалось gran_t
									}
									else {
										// узнать координаты центра грани и ориентацию в пространстве
										x_c = octree1->p1.x;
										y_c = 0.5*(octree1->linkE2->p0.y + octree1->linkE2->p3.y);
										z_c = 0.5*(octree1->linkE2->p0.z + octree1->linkE2->p4.z);
										iplane = YZ;
										patch_maxbound(iplane, s, ls, x_c, y_c, z_c, maxbound, binc2, lsid2);
										
									}
								}
							}

							if (octree1->linkE5 == NULL) {
								//maxbound++;
								// 21 сентября 2016 (вырождение ячейки).
							}
							else {
								if (octree1->linkE5->inum_TD == 0) {
									maxbound++;
								}
								else {
									// сосед существует.
									if (bvisit[octree1->linkE5->inum_TD - 1]) {
										// узел уже был посещён
										// Здесь неким образом модифицировалось gran_t
									}
									else {
										// узнать координаты центра грани и ориентацию в пространстве
										x_c = octree1->p1.x;
										y_c = 0.5*(octree1->linkE5->p0.y + octree1->linkE5->p3.y);
										z_c = 0.5*(octree1->linkE5->p0.z + octree1->linkE5->p4.z);
										iplane = YZ;
										patch_maxbound(iplane, s, ls, x_c, y_c, z_c, maxbound, binc3, lsid3);
									}
								}
							}

							if (octree1->linkE6 == NULL) {
								//maxbound++;
								// 21 сентября 2016 (вырождение ячейки).
							}
							else {
								if (octree1->linkE6->inum_TD == 0) {
									maxbound++;
								}
								else {
									// сосед существует.
									if (bvisit[octree1->linkE6->inum_TD - 1]) {
										// узел уже был посещён
										// Здесь неким образом модифицировалось gran_t
									}
									else {
										// узнать координаты центра грани и ориентацию в пространстве
										x_c = octree1->p1.x;
										y_c = 0.5*(octree1->linkE6->p0.y + octree1->linkE6->p3.y);
										z_c = 0.5*(octree1->linkE6->p0.z + octree1->linkE6->p4.z);
										iplane = YZ;
										patch_maxbound(iplane, s, ls, x_c, y_c, z_c, maxbound, binc4, lsid4);
									}
								}
							}
						}
					}
					
					// сосед не существует
					// записываем идентификатор грани и
					// увеличиваем счётчик граней.
					if (!octree1->b4W) {
						if (octree1->linkW == NULL) {
							maxbound++;
						}
						else {
							if (octree1->linkW->inum_TD == 0) {
								maxbound++;
							}
							else {
								// сосед существует.
								if (bvisit[octree1->linkW->inum_TD - 1]) {
									// узел уже был посещён
									// Здесь неким образом модифицировалось gran_t
								}
								else {
									bool binc;
									integer lsid;
									// узнать координаты центра грани и ориентацию в пространстве
									x_c = octree1->p0.x;
									y_c = 0.5*(octree1->p0.y + octree1->p3.y);
									z_c = 0.5*(octree1->p0.z + octree1->p4.z);
									iplane = YZ;
									patch_maxbound(iplane, s, ls, x_c, y_c, z_c, maxbound,binc,lsid);
								}
							}
						}
					}
					else {
						if ((octree1->linkW0 == NULL) && (octree1->linkW3 == NULL) && (octree1->linkW4 == NULL) && (octree1->linkW7 == NULL)) {
							// Зачем дробить если граница заканчивается NULL.
							maxbound++; // Вместо четырёх один.
						}
						else if ((octree1->linkW0 != NULL) && (octree1->linkW3 != NULL) && (octree1->linkW4 != NULL) && (octree1->linkW7 != NULL) &&
							(octree1->linkW0->inum_TD == 0) && (octree1->linkW3->inum_TD == 0) && (octree1->linkW4->inum_TD == 0) && (octree1->linkW7->inum_TD == 0)) {
							// Зачем дробить если граница заканчивается NULL.
							maxbound++; // Вместо четырёх один.
						}
						else if ((octree1->linkW0 != NULL) && (octree1->linkW3 == NULL) && (octree1->linkW4 != NULL) && (octree1->linkW7 == NULL) &&
							(octree1->linkW0->inum_TD == 0)  && (octree1->linkW4->inum_TD == 0) ) {
							// Зачем дробить если граница заканчивается NULL.
							maxbound++; // Вместо двух один.
						}						
						else if ((octree1->linkW0 != NULL) && (octree1->linkW3 != NULL) && (octree1->linkW4 == NULL) && (octree1->linkW7 == NULL) &&
							(octree1->linkW0->inum_TD == 0) && (octree1->linkW3->inum_TD == 0)) {
							// Зачем дробить если граница заканчивается NULL.
							maxbound++; // Вместо двух один.
						}
						else if ((octree1->linkW0 != NULL) && (octree1->linkW3 == NULL) && (octree1->linkW4 == NULL) && (octree1->linkW7 == NULL) &&
							(octree1->linkW0->inum_TD == 0)) {
							// Зачем дробить если граница заканчивается NULL.
							maxbound++; // Вместо двух один.
						}
						else {

							bool binc1, binc2, binc3, binc4;
							integer lsid1, lsid2, lsid3, lsid4;

							if (octree1->linkW0 == NULL) {
								//maxbound++;
								printf("Error!!! W0 is NULL\n");
								//getchar();
								system("PAUSE");
							}
							else {
								if (octree1->linkW0->inum_TD == 0) {
									maxbound++;
								}
								else {
									// сосед существует.
									if (bvisit[octree1->linkW0->inum_TD - 1]) {
										// узел уже был посещён
										// Здесь неким образом модифицировалось gran_t
									}
									else {
										// узнать координаты центра грани и ориентацию в пространстве
										x_c = octree1->p0.x;
										y_c = 0.5*(octree1->linkW0->p1.y + octree1->linkW0->p2.y);
										z_c = 0.5*(octree1->linkW0->p1.z + octree1->linkW0->p5.z);
										iplane = YZ;
										patch_maxbound(iplane, s, ls, x_c, y_c, z_c, maxbound, binc1, lsid1);
									}
								}
							}

							if (octree1->linkW3 == NULL) {
								//maxbound++;
								// 21 сентября 2016 (вырождение ячейки).
							}
							else {
								if (octree1->linkW3->inum_TD == 0) {
									maxbound++;
								}
								else {
									// сосед существует.
									if (bvisit[octree1->linkW3->inum_TD - 1]) {
										// узел уже был посещён
										// Здесь неким образом модифицировалось gran_t
									}
									else {
										// узнать координаты центра грани и ориентацию в пространстве
										x_c = octree1->p0.x;
										y_c = 0.5*(octree1->linkW3->p1.y + octree1->linkW3->p2.y);
										z_c = 0.5*(octree1->linkW3->p1.z + octree1->linkW3->p5.z);
										iplane = YZ;
										patch_maxbound(iplane, s, ls, x_c, y_c, z_c, maxbound, binc2, lsid2);
									}
								}
							}

							if (octree1->linkW4 == NULL) {
								//maxbound++;
								// 21 сентября 2016 (вырождение ячейки).
							}
							else {
								if (octree1->linkW4->inum_TD == 0) {
									maxbound++;
								}
								else {
									// сосед существует.
									if (bvisit[octree1->linkW4->inum_TD - 1]) {
										// узел уже был посещён
										// Здесь неким образом модифицировалось gran_t
									}
									else {
										// узнать координаты центра грани и ориентацию в пространстве
										x_c = octree1->p0.x;
										y_c = 0.5*(octree1->linkW4->p1.y + octree1->linkW4->p2.y);
										z_c = 0.5*(octree1->linkW4->p1.z + octree1->linkW4->p5.z);
										iplane = YZ;
										patch_maxbound(iplane, s, ls, x_c, y_c, z_c, maxbound, binc3, lsid3);
									}
								}
							}

							if (octree1->linkW7 == NULL) {
								//maxbound++;
								// 21 сентября 2016 (вырождение ячейки).
							}
							else {
								if (octree1->linkW7->inum_TD == 0) {
									maxbound++;
								}
								else {
									// сосед существует.
									if (bvisit[octree1->linkW7->inum_TD - 1]) {
										// узел уже был посещён
										// Здесь неким образом модифицировалось gran_t
									}
									else {
										// узнать координаты центра грани и ориентацию в пространстве
										x_c = octree1->p0.x;
										y_c = 0.5*(octree1->linkW7->p1.y + octree1->linkW7->p2.y);
										z_c = 0.5*(octree1->linkW7->p1.z + octree1->linkW7->p5.z);
										iplane = YZ;
										patch_maxbound(iplane, s, ls, x_c, y_c, z_c, maxbound, binc4, lsid4);
									}
								}
							}
						}
					}

					
					// сосед не существует
					// записываем идентификатор грани и
					// увеличиваем счётчик граней.
					if (!octree1->b4N) {
						if (octree1->linkN == NULL) {
							maxbound++;
						}
						else {
							if (octree1->linkN->inum_TD == 0) {
								maxbound++;
							}
							else {
								// сосед существует.
								if (bvisit[octree1->linkN->inum_TD - 1]) {
									// узел уже был посещён
									// Здесь неким образом модифицировалось gran_t
								}
								else {
									bool binc;
									integer lsid;
									// узнать координаты центра грани и ориентацию в пространстве
									x_c =  0.5*(octree1->p2.x + octree1->p3.x);
									y_c = octree1->p3.y;
									z_c = 0.5*(octree1->p3.z + octree1->p7.z);
									iplane = XZ;
									patch_maxbound(iplane, s, ls, x_c, y_c, z_c, maxbound,binc,lsid);
								}
							}
						}
					}
					else {

						if ((octree1->linkN2 == NULL) && (octree1->linkN3 == NULL) && (octree1->linkN6 == NULL) && (octree1->linkN7 == NULL)) {
							// Зачем дробить если граница заканчивается NULL.
							maxbound++; // Вместо четырёх один.
						}
						else if ((octree1->linkN2 != NULL) && (octree1->linkN3 != NULL) && (octree1->linkN6 != NULL) && (octree1->linkN7 != NULL) &&
							(octree1->linkN2->inum_TD == 0) && (octree1->linkN3->inum_TD == 0) && (octree1->linkN6->inum_TD == 0) && (octree1->linkN7->inum_TD == 0)) {
							// Зачем дробить если граница заканчивается NULL.
							maxbound++; // Вместо четырёх один.
						}
						else if ((octree1->linkN2 == NULL) && (octree1->linkN3 != NULL) && (octree1->linkN6 == NULL) && (octree1->linkN7 != NULL) &&
							 (octree1->linkN3->inum_TD == 0)  && (octree1->linkN7->inum_TD == 0)) {
							// Зачем дробить если граница заканчивается NULL.
							maxbound++; // Вместо двух один.
						}
						else if ((octree1->linkN2 != NULL) && (octree1->linkN3 != NULL) && (octree1->linkN6 == NULL) && (octree1->linkN7 == NULL) &&
							(octree1->linkN2->inum_TD == 0) && (octree1->linkN3->inum_TD == 0) ) {
							// Зачем дробить если граница заканчивается NULL.
							maxbound++; // Вместо двух один.
						}
						else if ((octree1->linkN2 == NULL) && (octree1->linkN3 != NULL) && (octree1->linkN6 == NULL) && (octree1->linkN7 == NULL) &&
							 (octree1->linkN3->inum_TD == 0)) {
							// Зачем дробить если граница заканчивается NULL.
							maxbound++; // Вместо двух один.
						}
						else {

							bool binc1, binc2, binc3, binc4;
							integer lsid1, lsid2, lsid3, lsid4;

							if (octree1->linkN2 == NULL) {
								//maxbound++;
								// 21 сентября 2016 (вырождение ячейки).
							}
							else {
								if (octree1->linkN2->inum_TD == 0) {
									maxbound++;
								}
								else {
									// сосед существует.
									if (bvisit[octree1->linkN2->inum_TD - 1]) {
										// узел уже был посещён
										// Здесь неким образом модифицировалось gran_t
									}
									else {
										// узнать координаты центра грани и ориентацию в пространстве
										x_c = 0.5*(octree1->linkN2->p0.x + octree1->linkN2->p1.x);
										y_c = octree1->p3.y;
										z_c = 0.5*(octree1->linkN2->p0.z + octree1->linkN2->p4.z);
										iplane = XZ;
										patch_maxbound(iplane, s, ls, x_c, y_c, z_c, maxbound, binc1, lsid1);
									}
								}
							}

							if (octree1->linkN3 == NULL) {
								//maxbound++;
								printf("Error!!! N3 is NULL\n");
								//getchar();
								system("PAUSE");
							}
							else {
								if (octree1->linkN3->inum_TD == 0) {
									maxbound++;
								}
								else {
									// сосед существует.
									if (bvisit[octree1->linkN3->inum_TD - 1]) {
										// узел уже был посещён
										// Здесь неким образом модифицировалось gran_t
									}
									else {
										// узнать координаты центра грани и ориентацию в пространстве
										x_c = 0.5*(octree1->linkN3->p0.x + octree1->linkN3->p1.x);
										y_c = octree1->p3.y;
										z_c = 0.5*(octree1->linkN3->p0.z + octree1->linkN3->p4.z);
										iplane = XZ;
										patch_maxbound(iplane, s, ls, x_c, y_c, z_c, maxbound, binc2, lsid2);
									}
								}
							}

							if (octree1->linkN6 == NULL) {
								//maxbound++;
								// 21 сентября 2016 (вырождение ячейки).
							}
							else {
								if (octree1->linkN6->inum_TD == 0) {
									maxbound++;
								}
								else {
									// сосед существует.
									if (bvisit[octree1->linkN6->inum_TD - 1]) {
										// узел уже был посещён
										// Здесь неким образом модифицировалось gran_t
									}
									else {
										// узнать координаты центра грани и ориентацию в пространстве
										x_c = 0.5*(octree1->linkN6->p0.x + octree1->linkN6->p1.x);
										y_c = octree1->p3.y;
										z_c = 0.5*(octree1->linkN6->p0.z + octree1->linkN6->p4.z);
										iplane = XZ;
										patch_maxbound(iplane, s, ls, x_c, y_c, z_c, maxbound, binc3, lsid3);
									}
								}
							}

							if (octree1->linkN7 == NULL) {
								//maxbound++;
								// 21 сентября 2016 (вырождение ячейки).
							}
							else {
								if (octree1->linkN7->inum_TD == 0) {
									maxbound++;
								}
								else {
									// сосед существует.
									if (bvisit[octree1->linkN7->inum_TD - 1]) {
										// узел уже был посещён
										// Здесь неким образом модифицировалось gran_t
									}
									else {
										// узнать координаты центра грани и ориентацию в пространстве
										x_c = 0.5*(octree1->linkN7->p0.x + octree1->linkN7->p1.x);
										y_c = octree1->p3.y;
										z_c = 0.5*(octree1->linkN7->p0.z + octree1->linkN7->p4.z);
										iplane = XZ;
										patch_maxbound(iplane, s, ls, x_c, y_c, z_c, maxbound, binc4, lsid4);
									}
								}
							}
						}
					}

					// сосед не существует
					// записываем идентификатор грани и
					// увеличиваем счётчик граней.
					if (!octree1->b4S) {
						if (octree1->linkS == NULL) {
							maxbound++;
						}
						else {
							if (octree1->linkS->inum_TD == 0) {
								maxbound++;
							}
							else {
								// сосед существует.
								if (bvisit[octree1->linkS->inum_TD - 1]) {
									// узел уже был посещён
									// Здесь неким образом модифицировалось gran_t
								}
								else {
									bool binc;
									integer lsid;
									// узнать координаты центра грани и ориентацию в пространстве
									x_c =  0.5*(octree1->p0.x + octree1->p1.x);
									y_c = octree1->p0.y;
									z_c = 0.5*(octree1->p0.z + octree1->p4.z);
									iplane = XZ;
									patch_maxbound(iplane, s, ls, x_c, y_c, z_c, maxbound,binc, lsid);
								}
							}
						}
					}
					else {
						if ((octree1->linkS0 == NULL) && (octree1->linkS1 == NULL) && (octree1->linkS4 == NULL) && (octree1->linkS5 == NULL)) {
							// Зачем дробить если граница заканчивается NULL.
							maxbound++; // Вместо четырёх один.
						}
						else if ((octree1->linkS0 != NULL) && (octree1->linkS1 != NULL) && (octree1->linkS4 != NULL) && (octree1->linkS5 != NULL) &&
							(octree1->linkS0->inum_TD == 0) && (octree1->linkS1->inum_TD == 0) && (octree1->linkS4->inum_TD == 0) && (octree1->linkS5->inum_TD == 0)) {
							// Зачем дробить если граница заканчивается NULL.
							maxbound++; // Вместо четырёх один.
						}
						else if ((octree1->linkS0 != NULL) && (octree1->linkS1 == NULL) && (octree1->linkS4 != NULL) && (octree1->linkS5 == NULL) &&
							(octree1->linkS0->inum_TD == 0)  && (octree1->linkS4->inum_TD == 0) ) {
							// Зачем дробить если граница заканчивается NULL.
							maxbound++; // Вместо двух один.
						}
						else if ((octree1->linkS0 != NULL) && (octree1->linkS1 != NULL) && (octree1->linkS4 == NULL) && (octree1->linkS5 == NULL) &&
							(octree1->linkS0->inum_TD == 0) && (octree1->linkS1->inum_TD == 0) ) {
							// Зачем дробить если граница заканчивается NULL.
							maxbound++; // Вместо двух один.
						}
						else if ((octree1->linkS0 != NULL) && (octree1->linkS1 == NULL) && (octree1->linkS4 == NULL) && (octree1->linkS5 == NULL) &&
							(octree1->linkS0->inum_TD == 0)) {
							// Зачем дробить если граница заканчивается NULL.
							maxbound++; // Вместо двух один.
						}
						else {

							bool binc1, binc2, binc3, binc4;
							integer lsid1, lsid2, lsid3, lsid4;

							if (octree1->linkS0 == NULL) {
								//maxbound++;
								printf("Error!!! S0 is NULL\n");
								//getchar();
								system("PAUSE");
							}
							else {
								if (octree1->linkS0->inum_TD == 0) {
									maxbound++;
								}
								else {
									// сосед существует.
									if (bvisit[octree1->linkS0->inum_TD - 1]) {
										// узел уже был посещён
										// Здесь неким образом модифицировалось gran_t
									}
									else {
										// узнать координаты центра грани и ориентацию в пространстве
										x_c = 0.5*(octree1->linkS0->p3.x + octree1->linkS0->p2.x);
										y_c = octree1->p0.y;
										z_c = 0.5*(octree1->linkS0->p3.z + octree1->linkS0->p7.z);
										iplane = XZ;
										patch_maxbound(iplane, s, ls, x_c, y_c, z_c, maxbound, binc1, lsid1);
									}
								}
							}

							if (octree1->linkS1 == NULL) {
								//maxbound++;
								// 21 сентября 2016 (вырождение ячейки).
							}
							else {
								if (octree1->linkS1->inum_TD == 0) {
									maxbound++;
								}
								else {
									// сосед существует.
									if (bvisit[octree1->linkS1->inum_TD - 1]) {
										// узел уже был посещён
										// Здесь неким образом модифицировалось gran_t
									}
									else {
										// узнать координаты центра грани и ориентацию в пространстве
										x_c = 0.5*(octree1->linkS1->p3.x + octree1->linkS1->p2.x);
										y_c = octree1->p0.y;
										z_c = 0.5*(octree1->linkS1->p3.z + octree1->linkS1->p7.z);
										iplane = XZ;
										patch_maxbound(iplane, s, ls, x_c, y_c, z_c, maxbound, binc2, lsid2);
									}
								}
							}

							if (octree1->linkS4 == NULL) {
								//maxbound++;
								// 21 сентября 2016 (вырождение ячейки).
							}
							else {
								if (octree1->linkS4->inum_TD == 0) {
									maxbound++;
								}
								else {
									// сосед существует.
									if (bvisit[octree1->linkS4->inum_TD - 1]) {
										// узел уже был посещён
										// Здесь неким образом модифицировалось gran_t
									}
									else {
										// узнать координаты центра грани и ориентацию в пространстве
										x_c = 0.5*(octree1->linkS4->p3.x + octree1->linkS4->p2.x);
										y_c = octree1->p0.y;
										z_c = 0.5*(octree1->linkS4->p3.z + octree1->linkS4->p7.z);
										iplane = XZ;
										patch_maxbound(iplane, s, ls, x_c, y_c, z_c, maxbound, binc3, lsid3);
									}
								}
							}

							if (octree1->linkS5 == NULL) {
								//maxbound++;
								// 21 сентября 2016 (вырождение ячейки).
							}
							else {
								if (octree1->linkS5->inum_TD == 0) {
									maxbound++;
								}
								else {
									// сосед существует.
									if (bvisit[octree1->linkS5->inum_TD - 1]) {
										// узел уже был посещён
										// Здесь неким образом модифицировалось gran_t
									}
									else {
										// узнать координаты центра грани и ориентацию в пространстве
										x_c = 0.5*(octree1->linkS5->p3.x + octree1->linkS5->p2.x);
										y_c = octree1->p0.y;
										z_c = 0.5*(octree1->linkS5->p3.z + octree1->linkS5->p7.z);
										iplane = XZ;
										patch_maxbound(iplane, s, ls, x_c, y_c, z_c, maxbound, binc4, lsid4);
									}
								}
							}
						}
					}

					// сосед не существует
					// записываем идентификатор грани и
					// увеличиваем счётчик граней.
					if (!octree1->b4T) {
						if (octree1->linkT == NULL) {
							maxbound++;
						}
						else {
							if (octree1->linkT->inum_TD == 0) {
								maxbound++;
							}
							else {
								// сосед существует.
								if (bvisit[octree1->linkT->inum_TD - 1]) {
									// узел уже был посещён
									// Здесь неким образом модифицировалось gran_t
								}
								else {
									bool binc;
									integer lsid;
									// узнать координаты центра грани и ориентацию в пространстве
									x_c = 0.5*(octree1->p4.x + octree1->p5.x);
									y_c = 0.5*(octree1->p4.y + octree1->p7.y);
									z_c = octree1->p4.z;
									iplane = XY;
									patch_maxbound(iplane, s, ls, x_c, y_c, z_c, maxbound, binc, lsid);
								}
							}
						}
					}
					else {

						if ((octree1->linkT4 == NULL) && (octree1->linkT5 == NULL) && (octree1->linkT6 == NULL) && (octree1->linkT7 == NULL)) {
							// Зачем дробить если граница заканчивается NULL.
							maxbound++; // Вместо четырёх один.
						}
						else if ((octree1->linkT4 != NULL) && (octree1->linkT5 != NULL) && (octree1->linkT6 != NULL) && (octree1->linkT7 != NULL) &&
							(octree1->linkT4->inum_TD == 0) && (octree1->linkT5->inum_TD == 0) && (octree1->linkT6->inum_TD == 0) && (octree1->linkT7->inum_TD == 0)) {
							// Зачем дробить если граница заканчивается NULL.
							maxbound++; // Вместо четырёх один.
						}
						else if ((octree1->linkT4 != NULL) && (octree1->linkT5 == NULL) && (octree1->linkT6 == NULL) && (octree1->linkT7 != NULL) &&
							(octree1->linkT4->inum_TD == 0)  && (octree1->linkT7->inum_TD == 0)) {
							// Зачем дробить если граница заканчивается NULL.
							maxbound++; // Вместо двух один.
						}
						else if ((octree1->linkT4 != NULL) && (octree1->linkT5 != NULL) && (octree1->linkT6 == NULL) && (octree1->linkT7 == NULL) &&
							(octree1->linkT4->inum_TD == 0) && (octree1->linkT5->inum_TD == 0)) {
							// Зачем дробить если граница заканчивается NULL.
							maxbound++; // Вместо двух один.
						}
						else if ((octree1->linkT4 != NULL) && (octree1->linkT5 == NULL) && (octree1->linkT6 == NULL) && (octree1->linkT7 == NULL) &&
							(octree1->linkT4->inum_TD == 0)) {
							// Зачем дробить если граница заканчивается NULL.
							maxbound++; // Вместо двух один.
						}
						else {

							bool binc1, binc2, binc3, binc4;
							integer lsid1, lsid2, lsid3, lsid4;

							if (octree1->linkT4 == NULL) {
								//maxbound++;
								printf("Error!!! T4 is NULL\n");
								//getchar();
								system("PAUSE");
							}
							else {
								if (octree1->linkT4->inum_TD == 0) {
									maxbound++;
								}
								else {
									// сосед существует.
									if (bvisit[octree1->linkT4->inum_TD - 1]) {
										// узел уже был посещён
										// Здесь неким образом модифицировалось gran_t
									}
									else {
										// узнать координаты центра грани и ориентацию в пространстве
										x_c = 0.5*(octree1->linkT4->p0.x + octree1->linkT4->p1.x);
										y_c = 0.5*(octree1->linkT4->p0.y + octree1->linkT4->p3.y);
										z_c = octree1->p4.z;
										iplane = XY;
										patch_maxbound(iplane, s, ls, x_c, y_c, z_c, maxbound, binc1, lsid1);
									}
								}
							}

							if (octree1->linkT5 == NULL) {
								//maxbound++;
								// 21 сентября 2016 (вырождение ячейки).
							}
							else {
								if (octree1->linkT5->inum_TD == 0) {
									maxbound++;
								}
								else {
									// сосед существует.
									if (bvisit[octree1->linkT5->inum_TD - 1]) {
										// узел уже был посещён
										// Здесь неким образом модифицировалось gran_t
									}
									else {
										// узнать координаты центра грани и ориентацию в пространстве
										x_c = 0.5*(octree1->linkT5->p0.x + octree1->linkT5->p1.x);
										y_c = 0.5*(octree1->linkT5->p0.y + octree1->linkT5->p3.y);
										z_c = octree1->p4.z;
										iplane = XY;
										patch_maxbound(iplane, s, ls, x_c, y_c, z_c, maxbound, binc2, lsid2);
									}
								}
							}

							if (octree1->linkT6 == NULL) {
								//maxbound++;
								// 21 сентября 2016 (вырождение ячейки).
							}
							else {
								if (octree1->linkT6->inum_TD == 0) {
									maxbound++;
								}
								else {
									// сосед существует.
									if (bvisit[octree1->linkT6->inum_TD - 1]) {
										// узел уже был посещён
										// Здесь неким образом модифицировалось gran_t
									}
									else {
										// узнать координаты центра грани и ориентацию в пространстве
										x_c = 0.5*(octree1->linkT6->p0.x + octree1->linkT6->p1.x);
										y_c = 0.5*(octree1->linkT6->p0.y + octree1->linkT6->p3.y);
										z_c = octree1->p4.z;
										iplane = XY;
										patch_maxbound(iplane, s, ls, x_c, y_c, z_c, maxbound, binc3, lsid3);
									}
								}
							}

							if (octree1->linkT7 == NULL) {
								//maxbound++;
								// 21 сентября 2016 (вырождение ячейки).
							}
							else {
								if (octree1->linkT7->inum_TD == 0) {
									maxbound++;
								}
								else {
									// сосед существует.
									if (bvisit[octree1->linkT7->inum_TD - 1]) {
										// узел уже был посещён
										// Здесь неким образом модифицировалось gran_t
									}
									else {
										// узнать координаты центра грани и ориентацию в пространстве
										x_c = 0.5*(octree1->linkT7->p0.x + octree1->linkT7->p1.x);
										y_c = 0.5*(octree1->linkT7->p0.y + octree1->linkT7->p3.y);
										z_c = octree1->p4.z;
										iplane = XY;
										patch_maxbound(iplane, s, ls, x_c, y_c, z_c, maxbound, binc4, lsid4);
									}
								}
							}
						}
					}

					// сосед не существует
					// записываем идентификатор грани и
					// увеличиваем счётчик граней.
					if (!octree1->b4B) {
						if (octree1->linkB == NULL) {
							maxbound++;
						}
						else {
							if (octree1->linkB->inum_TD == 0) {
								maxbound++;
							}
							else {
								// сосед существует.
								if (bvisit[octree1->linkB->inum_TD - 1]) {
									// узел уже был посещён
									// Здесь неким образом модифицировалось gran_t
								}
								else {
									bool binc;
									integer lsid;
									// узнать координаты центра грани и ориентацию в пространстве
									x_c = 0.5*(octree1->p0.x + octree1->p1.x);
									y_c = 0.5*(octree1->p0.y + octree1->p3.y);
									z_c = octree1->p0.z;
									iplane = XY;
									patch_maxbound(iplane, s, ls, x_c, y_c, z_c, maxbound, binc, lsid);
								}
							}
						}
					}
					else {

						if ((octree1->linkB0 == NULL) && (octree1->linkB1 == NULL) && (octree1->linkB2 == NULL) && (octree1->linkB3 == NULL)) {
							// Зачем дробить если граница заканчивается NULL.
							maxbound++; // Вместо четырёх один.
						}
						else if ((octree1->linkB0 != NULL) && (octree1->linkB1 != NULL) && (octree1->linkB2 != NULL) && (octree1->linkB3 != NULL) &&
							(octree1->linkB0->inum_TD == 0) && (octree1->linkB1->inum_TD == 0) && (octree1->linkB2->inum_TD == 0) && (octree1->linkB3->inum_TD == 0)) {
							// Зачем дробить если граница заканчивается NULL.
							maxbound++; // Вместо четырёх один.
						}
						else if ((octree1->linkB0 != NULL) && (octree1->linkB1 == NULL) && (octree1->linkB2 == NULL) && (octree1->linkB3 != NULL) &&
							(octree1->linkB0->inum_TD == 0)  && (octree1->linkB3->inum_TD == 0)) {
							// Зачем дробить если граница заканчивается NULL.
							maxbound++; // Вместо двух один.
						}
						else if ((octree1->linkB0 != NULL) && (octree1->linkB1 != NULL) && (octree1->linkB2 == NULL) && (octree1->linkB3 == NULL) &&
							(octree1->linkB0->inum_TD == 0) && (octree1->linkB1->inum_TD == 0) ) {
							// Зачем дробить если граница заканчивается NULL.
							maxbound++; // Вместо двух один.
						}
						else if ((octree1->linkB0 != NULL) && (octree1->linkB1 == NULL) && (octree1->linkB2 == NULL) && (octree1->linkB3 == NULL) &&
							(octree1->linkB0->inum_TD == 0)) {
							// Зачем дробить если граница заканчивается NULL.
							maxbound++; // Вместо двух один.
						}
						else {

							bool binc1, binc2, binc3, binc4;
							integer lsid1, lsid2, lsid3, lsid4;

							if (octree1->linkB0 == NULL) {
								//maxbound++;
								printf("Error!!! B0 is NULL\n");
								//getchar();
								system("PAUSE");
							}
							else {
								if (octree1->linkB0->inum_TD == 0) {
									maxbound++;
								}
								else {
									// сосед существует.
									if (bvisit[octree1->linkB0->inum_TD - 1]) {
										// узел уже был посещён
										// Здесь неким образом модифицировалось gran_t
									}
									else {
										// узнать координаты центра грани и ориентацию в пространстве
										x_c = 0.5*(octree1->linkB0->p4.x + octree1->linkB0->p5.x);
										y_c = 0.5*(octree1->linkB0->p4.y + octree1->linkB0->p7.y);
										z_c = octree1->p0.z;
										iplane = XY;
										patch_maxbound(iplane, s, ls, x_c, y_c, z_c, maxbound, binc1, lsid1);
									}
								}
							}

							if (octree1->linkB1 == NULL) {
								//maxbound++;
								// 21 сентября 2016 (вырождение ячейки).
							}
							else {
								if (octree1->linkB1->inum_TD == 0) {
									maxbound++;
								}
								else {
									// сосед существует.
									if (bvisit[octree1->linkB1->inum_TD - 1]) {
										// узел уже был посещён
										// Здесь неким образом модифицировалось gran_t
									}
									else {
										// узнать координаты центра грани и ориентацию в пространстве
										x_c = 0.5*(octree1->linkB1->p4.x + octree1->linkB1->p5.x);
										y_c = 0.5*(octree1->linkB1->p4.y + octree1->linkB1->p7.y);
										z_c = octree1->p0.z;
										iplane = XY;
										patch_maxbound(iplane, s, ls, x_c, y_c, z_c, maxbound, binc2, lsid2);
									}
								}
							}

							if (octree1->linkB2 == NULL) {
								//maxbound++;
								// 21 сентября 2016 (вырождение ячейки).
							}
							else {
								if (octree1->linkB2->inum_TD == 0) {
									maxbound++;
								}
								else {
									// сосед существует.
									if (bvisit[octree1->linkB2->inum_TD - 1]) {
										// узел уже был посещён
										// Здесь неким образом модифицировалось gran_t
									}
									else {
										// узнать координаты центра грани и ориентацию в пространстве
										x_c = 0.5*(octree1->linkB2->p4.x + octree1->linkB2->p5.x);
										y_c = 0.5*(octree1->linkB2->p4.y + octree1->linkB2->p7.y);
										z_c = octree1->p0.z;
										iplane = XY;
										patch_maxbound(iplane, s, ls, x_c, y_c, z_c, maxbound, binc3, lsid3);
									}
								}
							}

							if (octree1->linkB3 == NULL) {
								//maxbound++;
								// 21 сентября 2016 (вырождение ячейки).
							}
							else {
								if (octree1->linkB3->inum_TD == 0) {
									maxbound++;
								}
								else {
									// сосед существует.
									if (bvisit[octree1->linkB3->inum_TD - 1]) {
										// узел уже был посещён
										// Здесь неким образом модифицировалось gran_t
									}
									else {
										// узнать координаты центра грани и ориентацию в пространстве
										x_c = 0.5*(octree1->linkB3->p4.x + octree1->linkB3->p5.x);
										y_c = 0.5*(octree1->linkB3->p4.y + octree1->linkB3->p7.y);
										z_c = octree1->p0.z;
										iplane = XY;
										patch_maxbound(iplane, s, ls, x_c, y_c, z_c, maxbound, binc4, lsid4);
									}
								}
							}
						}
					}


					bvisit[maxelm] = true; // Узел был посещён.
					maxelm++;
				}
				//octree1->inum_TD = 0; // По умолчанию не принадлежит расчётной области.
				octree1 = NULL;
				my_ALICE_STACK[top_ALICE_STACK - 1].link = NULL;
				top_ALICE_STACK--;
			}
			else {
				// продолжаем добираться до листьев.
				STACK_ALICE buf1 = my_ALICE_STACK[top_ALICE_STACK - 1];
				STACK_ALICE* buf = &buf1;
				top_ALICE_STACK--;
				if (buf->link->link0 != NULL) {
					my_ALICE_STACK[top_ALICE_STACK].link = (buf->link->link0);
					my_ALICE_STACK[top_ALICE_STACK].maxx = buf->link->link0->maxx;
					my_ALICE_STACK[top_ALICE_STACK].minx = buf->link->link0->minx;
					my_ALICE_STACK[top_ALICE_STACK].maxy = buf->link->link0->maxy;
					my_ALICE_STACK[top_ALICE_STACK].miny = buf->link->link0->miny;
					my_ALICE_STACK[top_ALICE_STACK].maxz = buf->link->link0->maxz;
					my_ALICE_STACK[top_ALICE_STACK].minz = buf->link->link0->minz;
					top_ALICE_STACK++;
				}
				if (buf->link->link1 != NULL) {
					my_ALICE_STACK[top_ALICE_STACK].link = (buf->link->link1);
					my_ALICE_STACK[top_ALICE_STACK].maxx = buf->link->link1->maxx;
					my_ALICE_STACK[top_ALICE_STACK].minx = buf->link->link1->minx;
					my_ALICE_STACK[top_ALICE_STACK].maxy = buf->link->link1->maxy;
					my_ALICE_STACK[top_ALICE_STACK].miny = buf->link->link1->miny;
					my_ALICE_STACK[top_ALICE_STACK].maxz = buf->link->link1->maxz;
					my_ALICE_STACK[top_ALICE_STACK].minz = buf->link->link1->minz;
					top_ALICE_STACK++;
				}
				if (buf->link->link2 != NULL) {
					my_ALICE_STACK[top_ALICE_STACK].link = (buf->link->link2);
					my_ALICE_STACK[top_ALICE_STACK].maxx = buf->link->link2->maxx;
					my_ALICE_STACK[top_ALICE_STACK].minx = buf->link->link2->minx;
					my_ALICE_STACK[top_ALICE_STACK].maxy = buf->link->link2->maxy;
					my_ALICE_STACK[top_ALICE_STACK].miny = buf->link->link2->miny;
					my_ALICE_STACK[top_ALICE_STACK].maxz = buf->link->link2->maxz;
					my_ALICE_STACK[top_ALICE_STACK].minz = buf->link->link2->minz;
					top_ALICE_STACK++;
				}
				if (buf->link->link3 != NULL) {
					my_ALICE_STACK[top_ALICE_STACK].link = (buf->link->link3);
					my_ALICE_STACK[top_ALICE_STACK].maxx = buf->link->link3->maxx;
					my_ALICE_STACK[top_ALICE_STACK].minx = buf->link->link3->minx;
					my_ALICE_STACK[top_ALICE_STACK].maxy = buf->link->link3->maxy;
					my_ALICE_STACK[top_ALICE_STACK].miny = buf->link->link3->miny;
					my_ALICE_STACK[top_ALICE_STACK].maxz = buf->link->link3->maxz;
					my_ALICE_STACK[top_ALICE_STACK].minz = buf->link->link3->minz;
					top_ALICE_STACK++;
				}
				if (buf->link->link4 != NULL) {
					my_ALICE_STACK[top_ALICE_STACK].link = (buf->link->link4);
					my_ALICE_STACK[top_ALICE_STACK].maxx = buf->link->link4->maxx;
					my_ALICE_STACK[top_ALICE_STACK].minx = buf->link->link4->minx;
					my_ALICE_STACK[top_ALICE_STACK].maxy = buf->link->link4->maxy;
					my_ALICE_STACK[top_ALICE_STACK].miny = buf->link->link4->miny;
					my_ALICE_STACK[top_ALICE_STACK].maxz = buf->link->link4->maxz;
					my_ALICE_STACK[top_ALICE_STACK].minz = buf->link->link4->minz;
					top_ALICE_STACK++;
				}
				if (buf->link->link5 != NULL) {
					my_ALICE_STACK[top_ALICE_STACK].link = (buf->link->link5);
					my_ALICE_STACK[top_ALICE_STACK].maxx = buf->link->link5->maxx;
					my_ALICE_STACK[top_ALICE_STACK].minx = buf->link->link5->minx;
					my_ALICE_STACK[top_ALICE_STACK].maxy = buf->link->link5->maxy;
					my_ALICE_STACK[top_ALICE_STACK].miny = buf->link->link5->miny;
					my_ALICE_STACK[top_ALICE_STACK].maxz = buf->link->link5->maxz;
					my_ALICE_STACK[top_ALICE_STACK].minz = buf->link->link5->minz;
					top_ALICE_STACK++;
				}
				if (buf->link->link6 != NULL) {
					my_ALICE_STACK[top_ALICE_STACK].link = (buf->link->link6);
					my_ALICE_STACK[top_ALICE_STACK].maxx = buf->link->link6->maxx;
					my_ALICE_STACK[top_ALICE_STACK].minx = buf->link->link6->minx;
					my_ALICE_STACK[top_ALICE_STACK].maxy = buf->link->link6->maxy;
					my_ALICE_STACK[top_ALICE_STACK].miny = buf->link->link6->miny;
					my_ALICE_STACK[top_ALICE_STACK].maxz = buf->link->link6->maxz;
					my_ALICE_STACK[top_ALICE_STACK].minz = buf->link->link6->minz;
					top_ALICE_STACK++;
				}
				if (buf->link->link7 != NULL) {
					my_ALICE_STACK[top_ALICE_STACK].link = (buf->link->link7);
					my_ALICE_STACK[top_ALICE_STACK].maxx = buf->link->link7->maxx;
					my_ALICE_STACK[top_ALICE_STACK].minx = buf->link->link7->minx;
					my_ALICE_STACK[top_ALICE_STACK].maxy = buf->link->link7->maxy;
					my_ALICE_STACK[top_ALICE_STACK].miny = buf->link->link7->miny;
					my_ALICE_STACK[top_ALICE_STACK].maxz = buf->link->link7->maxz;
					my_ALICE_STACK[top_ALICE_STACK].minz = buf->link->link7->minz;
					top_ALICE_STACK++;
				}
			}
		}
		//}
		//getchar();
	}

	delete[] bvisit;
	bvisit = NULL;
} // calculate_max_bound_temp


// Вычисляет количество граничных элементов модели для температуры.
// Внутренние источники тепла также пронумерованы.
void calculate_max_bound_flow(octTree* &oc, integer &maxbound, integer maxelm_memo, BLOCK* b, integer lb,
	SOURCE* s, integer ls) {
	maxbound = 0; // инициализация.
	integer maxelm = 0;
	bool *bvisit = NULL;
	bvisit = new bool[maxelm_memo];
	// оператор new не требует проверки на null.
	//if (bvisit == NULL) {
		// недостаточно памяти на данном оборудовании.
		//printf("Problem : not enough memory on your equipment for bvisit constr struct_alice...\n");
		//printf("Please any key to exit...\n");
		//getchar();
		//system("pause");
		//exit(1);
	//}
	for (integer j = 0; j<maxelm_memo; j++) bvisit[j] = false; // признак посещения узла.
	doublereal x_c, y_c, z_c;
	integer iplane;

	top_ALICE_STACK = 0;
	if (oc->link0 != NULL) {
		my_ALICE_STACK[top_ALICE_STACK].link = (oc->link0);
		my_ALICE_STACK[top_ALICE_STACK].maxx = oc->link0->maxx;
		my_ALICE_STACK[top_ALICE_STACK].minx = oc->link0->minx;
		my_ALICE_STACK[top_ALICE_STACK].maxy = oc->link0->maxy;
		my_ALICE_STACK[top_ALICE_STACK].miny = oc->link0->miny;
		my_ALICE_STACK[top_ALICE_STACK].maxz = oc->link0->maxz;
		my_ALICE_STACK[top_ALICE_STACK].minz = oc->link0->minz;
		top_ALICE_STACK++;
	}
	if (oc->link1 != NULL) {
		my_ALICE_STACK[top_ALICE_STACK].link = (oc->link1);
		my_ALICE_STACK[top_ALICE_STACK].maxx = oc->link1->maxx;
		my_ALICE_STACK[top_ALICE_STACK].minx = oc->link1->minx;
		my_ALICE_STACK[top_ALICE_STACK].maxy = oc->link1->maxy;
		my_ALICE_STACK[top_ALICE_STACK].miny = oc->link1->miny;
		my_ALICE_STACK[top_ALICE_STACK].maxz = oc->link1->maxz;
		my_ALICE_STACK[top_ALICE_STACK].minz = oc->link1->minz;
		top_ALICE_STACK++;
	}
	if (oc->link2 != NULL) {
		my_ALICE_STACK[top_ALICE_STACK].link = (oc->link2);
		my_ALICE_STACK[top_ALICE_STACK].maxx = oc->link2->maxx;
		my_ALICE_STACK[top_ALICE_STACK].minx = oc->link2->minx;
		my_ALICE_STACK[top_ALICE_STACK].maxy = oc->link2->maxy;
		my_ALICE_STACK[top_ALICE_STACK].miny = oc->link2->miny;
		my_ALICE_STACK[top_ALICE_STACK].maxz = oc->link2->maxz;
		my_ALICE_STACK[top_ALICE_STACK].minz = oc->link2->minz;
		top_ALICE_STACK++;
	}
	if (oc->link3 != NULL) {
		my_ALICE_STACK[top_ALICE_STACK].link = (oc->link3);
		my_ALICE_STACK[top_ALICE_STACK].maxx = oc->link3->maxx;
		my_ALICE_STACK[top_ALICE_STACK].minx = oc->link3->minx;
		my_ALICE_STACK[top_ALICE_STACK].maxy = oc->link3->maxy;
		my_ALICE_STACK[top_ALICE_STACK].miny = oc->link3->miny;
		my_ALICE_STACK[top_ALICE_STACK].maxz = oc->link3->maxz;
		my_ALICE_STACK[top_ALICE_STACK].minz = oc->link3->minz;
		top_ALICE_STACK++;
	}
	if (oc->link4 != NULL) {
		my_ALICE_STACK[top_ALICE_STACK].link = (oc->link4);
		my_ALICE_STACK[top_ALICE_STACK].maxx = oc->link4->maxx;
		my_ALICE_STACK[top_ALICE_STACK].minx = oc->link4->minx;
		my_ALICE_STACK[top_ALICE_STACK].maxy = oc->link4->maxy;
		my_ALICE_STACK[top_ALICE_STACK].miny = oc->link4->miny;
		my_ALICE_STACK[top_ALICE_STACK].maxz = oc->link4->maxz;
		my_ALICE_STACK[top_ALICE_STACK].minz = oc->link4->minz;
		top_ALICE_STACK++;
	}
	if (oc->link5 != NULL) {
		my_ALICE_STACK[top_ALICE_STACK].link = (oc->link5);
		my_ALICE_STACK[top_ALICE_STACK].maxx = oc->link5->maxx;
		my_ALICE_STACK[top_ALICE_STACK].minx = oc->link5->minx;
		my_ALICE_STACK[top_ALICE_STACK].maxy = oc->link5->maxy;
		my_ALICE_STACK[top_ALICE_STACK].miny = oc->link5->miny;
		my_ALICE_STACK[top_ALICE_STACK].maxz = oc->link5->maxz;
		my_ALICE_STACK[top_ALICE_STACK].minz = oc->link5->minz;
		top_ALICE_STACK++;
	}
	if (oc->link6 != NULL) {
		my_ALICE_STACK[top_ALICE_STACK].link = (oc->link6);
		my_ALICE_STACK[top_ALICE_STACK].maxx = oc->link6->maxx;
		my_ALICE_STACK[top_ALICE_STACK].minx = oc->link6->minx;
		my_ALICE_STACK[top_ALICE_STACK].maxy = oc->link6->maxy;
		my_ALICE_STACK[top_ALICE_STACK].miny = oc->link6->miny;
		my_ALICE_STACK[top_ALICE_STACK].maxz = oc->link6->maxz;
		my_ALICE_STACK[top_ALICE_STACK].minz = oc->link6->minz;
		top_ALICE_STACK++;
	}
	if (oc->link7 != NULL) {
		my_ALICE_STACK[top_ALICE_STACK].link = (oc->link7);
		my_ALICE_STACK[top_ALICE_STACK].maxx = oc->link7->maxx;
		my_ALICE_STACK[top_ALICE_STACK].minx = oc->link7->minx;
		my_ALICE_STACK[top_ALICE_STACK].maxy = oc->link7->maxy;
		my_ALICE_STACK[top_ALICE_STACK].miny = oc->link7->miny;
		my_ALICE_STACK[top_ALICE_STACK].maxz = oc->link7->maxz;
		my_ALICE_STACK[top_ALICE_STACK].minz = oc->link7->minz;
		top_ALICE_STACK++;
	}
	while (top_ALICE_STACK > 0) {
		if (my_ALICE_STACK[top_ALICE_STACK - 1].link != NULL) {
			if (my_ALICE_STACK[top_ALICE_STACK - 1].link->dlist == true) {
				// Гасим информацию о посещениях.
				octTree* octree1 = my_ALICE_STACK[top_ALICE_STACK - 1].link;

				if (!((octree1->link0 == NULL) && (octree1->link1 == NULL) && (octree1->link2 == NULL) && (octree1->link3 == NULL) && (octree1->link4 == NULL) && (octree1->link5 == NULL) && (octree1->link6 == NULL) && (octree1->link7 == NULL))) {
					printf("eto ne list error!!!\n");
					//getchar();
					system("PAUSE");
				}

				TOCHKA p;
				p.x = 0.125*(octree1->p0.x + octree1->p1.x + octree1->p2.x + octree1->p3.x + octree1->p4.x + octree1->p5.x + octree1->p6.x + octree1->p7.x);
				p.y = 0.125*(octree1->p0.y + octree1->p1.y + octree1->p2.y + octree1->p3.y + octree1->p4.y + octree1->p5.y + octree1->p6.y + octree1->p7.y);
				p.z = 0.125*(octree1->p0.z + octree1->p1.z + octree1->p2.z + octree1->p3.z + octree1->p4.z + octree1->p5.z + octree1->p6.z + octree1->p7.z);
				integer ib;
				bool inDomain = false;

				inDomain = in_model_flow(p, ib, b, lb); // HYDRODINAMIC  CFD!!!

				if ((p.x < b[0].g.xS) || (p.x > b[0].g.xE) || (p.y < b[0].g.yS) || (p.y > b[0].g.yE) || (p.z < b[0].g.zS) || (p.z > b[0].g.zE)) {
					// Лежит за пределами кабинета.
					inDomain = false;
				}
				if (inDomain) {
					// Узел принадлекжит расчётной области и его номер maxelm==octree1->inum_FD.

					// сосед не существует
					// записываем идентификатор грани и
					// увеличиваем счётчик граней.
					if (!octree1->b4E) {
						if (octree1->linkE == NULL) {
							maxbound++;
						}
						else {
							if (octree1->linkE->inum_FD == 0) {
								maxbound++;
							}
							else {
								// сосед существует.
								if (bvisit[octree1->linkE->inum_FD - 1]) {
									// узел уже был посещён
									// Здесь неким образом модифицировалось gran_t
								}
								else {
									bool binc;
									integer lsid;
									// узнать координаты центра грани и ориентацию в пространстве
									x_c = octree1->p1.x;
									y_c = 0.5*(octree1->p1.y + octree1->p2.y);
									z_c = 0.5*(octree1->p1.z + octree1->p5.z);
									iplane = YZ;
									patch_maxbound(iplane, s, ls, x_c, y_c, z_c, maxbound, binc, lsid);
									if (binc) {
										printf("error in calculate_max_bound_flow...\n");
										printf("model is incorrect. FATALL error!!!\n");
										printf("internal flat heat source can not be located\n");
										printf("inside the liquid.\n");
										printf("Please, press any key to exit...\n");
										//getchar();
										system("PAUSE");
										exit(1);
									}

								}
							}
						}
					}
					else {

						if ((octree1->linkE1 == NULL) && (octree1->linkE2 == NULL) && (octree1->linkE5 == NULL) && (octree1->linkE6 == NULL)) {
							// Зачем дробить если граница заканчивается NULL.
							maxbound++; // Вместо четырёх один.
						}
						else if ((octree1->linkE1 != NULL) && (octree1->linkE2 != NULL) && (octree1->linkE5 != NULL) && (octree1->linkE6 != NULL) &&
							(octree1->linkE1->inum_FD == 0) && (octree1->linkE2->inum_FD == 0) && (octree1->linkE5->inum_FD == 0) && (octree1->linkE6->inum_FD == 0)) {
							// Зачем дробить если граница заканчивается NULL.
							maxbound++; // Вместо четырёх один.
						}
						else if ((octree1->linkE1 != NULL) && (octree1->linkE2 == NULL) && (octree1->linkE5 != NULL) && (octree1->linkE6 == NULL) &&
							(octree1->linkE1->inum_FD == 0) && (octree1->linkE5->inum_FD == 0)) {
							// Зачем дробить если граница заканчивается NULL.
							maxbound++; // Вместо двух один.
						}
						else if ((octree1->linkE1 != NULL) && (octree1->linkE2 != NULL) && (octree1->linkE5 == NULL) && (octree1->linkE6 == NULL) &&
							(octree1->linkE1->inum_FD == 0) && (octree1->linkE2->inum_FD == 0)) {
							// Зачем дробить если граница заканчивается NULL.
							maxbound++; // Вместо двух один.
						}
						else if ((octree1->linkE1 != NULL) && (octree1->linkE2 == NULL) && (octree1->linkE5 == NULL) && (octree1->linkE6 == NULL) &&
							(octree1->linkE1->inum_FD == 0)) {
							// Зачем дробить если граница заканчивается NULL.
							maxbound++; // Вместо двух один.
						}
						else {

							bool binc1, binc2, binc3, binc4;
							integer lsid1, lsid2, lsid3, lsid4;

							if (octree1->linkE1 == NULL) {
								//maxbound++;
								printf("Error!!! E1 is NULL\n");
								//getchar();
								system("PAUSE");
							}
							else {
								if (octree1->linkE1->inum_FD == 0) {
									maxbound++;

								}
								else {
									// сосед существует.
									if (bvisit[octree1->linkE1->inum_FD - 1]) {
										// узел уже был посещён
										// Здесь неким образом модифицировалось gran_t
									}
									else {
										// узнать координаты центра грани и ориентацию в пространстве
										x_c = octree1->p1.x;
										y_c = 0.5*(octree1->linkE1->p0.y + octree1->linkE1->p3.y);
										z_c = 0.5*(octree1->linkE1->p0.z + octree1->linkE1->p4.z);
										iplane = YZ;
										patch_maxbound(iplane, s, ls, x_c, y_c, z_c, maxbound, binc1, lsid1);

									}
								}
							}

							if (octree1->linkE2 == NULL) {
								//maxbound++;
								// 21 сентября 2016 (вырождение ячейки).
							}
							else {
								if (octree1->linkE2->inum_FD == 0) {
									maxbound++;

								}
								else {
									// сосед существует.
									if (bvisit[octree1->linkE2->inum_FD - 1]) {
										// узел уже был посещён
										// Здесь неким образом модифицировалось gran_t
									}
									else {
										// узнать координаты центра грани и ориентацию в пространстве
										x_c = octree1->p1.x;
										y_c = 0.5*(octree1->linkE2->p0.y + octree1->linkE2->p3.y);
										z_c = 0.5*(octree1->linkE2->p0.z + octree1->linkE2->p4.z);
										iplane = YZ;
										patch_maxbound(iplane, s, ls, x_c, y_c, z_c, maxbound, binc2, lsid2);

									}
								}
							}

							if (octree1->linkE5 == NULL) {
								//maxbound++;
								// 21 сентября 2016 (вырождение ячейки).
							}
							else {
								if (octree1->linkE5->inum_FD == 0) {
									maxbound++;
								}
								else {
									// сосед существует.
									if (bvisit[octree1->linkE5->inum_FD - 1]) {
										// узел уже был посещён
										// Здесь неким образом модифицировалось gran_t
									}
									else {
										// узнать координаты центра грани и ориентацию в пространстве
										x_c = octree1->p1.x;
										y_c = 0.5*(octree1->linkE5->p0.y + octree1->linkE5->p3.y);
										z_c = 0.5*(octree1->linkE5->p0.z + octree1->linkE5->p4.z);
										iplane = YZ;
										patch_maxbound(iplane, s, ls, x_c, y_c, z_c, maxbound, binc3, lsid3);
									}
								}
							}

							if (octree1->linkE6 == NULL) {
								//maxbound++;
								// 21 сентября 2016 (вырождение ячейки).
							}
							else {
								if (octree1->linkE6->inum_FD == 0) {
									maxbound++;
								}
								else {
									// сосед существует.
									if (bvisit[octree1->linkE6->inum_FD - 1]) {
										// узел уже был посещён
										// Здесь неким образом модифицировалось gran_t
									}
									else {
										// узнать координаты центра грани и ориентацию в пространстве
										x_c = octree1->p1.x;
										y_c = 0.5*(octree1->linkE6->p0.y + octree1->linkE6->p3.y);
										z_c = 0.5*(octree1->linkE6->p0.z + octree1->linkE6->p4.z);
										iplane = YZ;
										patch_maxbound(iplane, s, ls, x_c, y_c, z_c, maxbound, binc4, lsid4);
									}
								}
							}
						}
					}

					// сосед не существует
					// записываем идентификатор грани и
					// увеличиваем счётчик граней.
					if (!octree1->b4W) {
						if (octree1->linkW == NULL) {
							maxbound++;
						}
						else {
							if (octree1->linkW->inum_FD == 0) {
								maxbound++;
							}
							else {
								// сосед существует.
								if (bvisit[octree1->linkW->inum_FD - 1]) {
									// узел уже был посещён
									// Здесь неким образом модифицировалось gran_t
								}
								else {
									bool binc;
									integer lsid;
									// узнать координаты центра грани и ориентацию в пространстве
									x_c = octree1->p0.x;
									y_c = 0.5*(octree1->p0.y + octree1->p3.y);
									z_c = 0.5*(octree1->p0.z + octree1->p4.z);
									iplane = YZ;
									patch_maxbound(iplane, s, ls, x_c, y_c, z_c, maxbound, binc, lsid);
									if (binc) {
										printf("error in calculate_max_bound_flow...\n");
										printf("model is incorrect. FATALL error!!!\n");
										printf("internal flat heat source can not be located\n");
										printf("inside the liquid.\n");
										printf("Please, press any key to exit...\n");
										//getchar();
										system("PAUSE");
										exit(1);
									}

								}
							}
						}
					}
					else {
						if ((octree1->linkW0 == NULL) && (octree1->linkW3 == NULL) && (octree1->linkW4 == NULL) && (octree1->linkW7 == NULL)) {
							// Зачем дробить если граница заканчивается NULL.
							maxbound++; // Вместо четырёх один.
						}
						else if ((octree1->linkW0 != NULL) && (octree1->linkW3 != NULL) && (octree1->linkW4 != NULL) && (octree1->linkW7 != NULL) &&
							(octree1->linkW0->inum_FD == 0) && (octree1->linkW3->inum_FD == 0) && (octree1->linkW4->inum_FD == 0) && (octree1->linkW7->inum_FD == 0)) {
							// Зачем дробить если граница заканчивается NULL.
							maxbound++; // Вместо четырёх один.
						}
						else if ((octree1->linkW0 != NULL) && (octree1->linkW3 == NULL) && (octree1->linkW4 != NULL) && (octree1->linkW7 == NULL) &&
							(octree1->linkW0->inum_FD == 0) && (octree1->linkW4->inum_FD == 0)) {
							// Зачем дробить если граница заканчивается NULL.
							maxbound++; // Вместо двух один.
						}
						else if ((octree1->linkW0 != NULL) && (octree1->linkW3 != NULL) && (octree1->linkW4 == NULL) && (octree1->linkW7 == NULL) &&
							(octree1->linkW0->inum_FD == 0) && (octree1->linkW3->inum_FD == 0)) {
							// Зачем дробить если граница заканчивается NULL.
							maxbound++; // Вместо двух один.
						}
						else if ((octree1->linkW0 != NULL) && (octree1->linkW3 == NULL) && (octree1->linkW4 == NULL) && (octree1->linkW7 == NULL) &&
							(octree1->linkW0->inum_FD == 0)) {
							// Зачем дробить если граница заканчивается NULL.
							maxbound++; // Вместо двух один.
						}
						else {

							bool binc1, binc2, binc3, binc4;
							integer lsid1, lsid2, lsid3, lsid4;

							if (octree1->linkW0 == NULL) {
								//maxbound++;
								printf("Error!!! W0 is NULL\n");
								//getchar();
								system("PAUSE");
							}
							else {
								if (octree1->linkW0->inum_FD == 0) {
									maxbound++;
								}
								else {
									// сосед существует.
									if (bvisit[octree1->linkW0->inum_FD - 1]) {
										// узел уже был посещён
										// Здесь неким образом модифицировалось gran_t
									}
									else {
										// узнать координаты центра грани и ориентацию в пространстве
										x_c = octree1->p0.x;
										y_c = 0.5*(octree1->linkW0->p1.y + octree1->linkW0->p2.y);
										z_c = 0.5*(octree1->linkW0->p1.z + octree1->linkW0->p5.z);
										iplane = YZ;
										patch_maxbound(iplane, s, ls, x_c, y_c, z_c, maxbound, binc1, lsid1);
									}
								}
							}

							if (octree1->linkW3 == NULL) {
								//maxbound++;
								// 21 сентября 2016 (вырождение ячейки).
							}
							else {
								if (octree1->linkW3->inum_FD == 0) {
									maxbound++;
								}
								else {
									// сосед существует.
									if (bvisit[octree1->linkW3->inum_FD - 1]) {
										// узел уже был посещён
										// Здесь неким образом модифицировалось gran_t
									}
									else {
										// узнать координаты центра грани и ориентацию в пространстве
										x_c = octree1->p0.x;
										y_c = 0.5*(octree1->linkW3->p1.y + octree1->linkW3->p2.y);
										z_c = 0.5*(octree1->linkW3->p1.z + octree1->linkW3->p5.z);
										iplane = YZ;
										patch_maxbound(iplane, s, ls, x_c, y_c, z_c, maxbound, binc2, lsid2);
									}
								}
							}

							if (octree1->linkW4 == NULL) {
								//maxbound++;
								// 21 сентября 2016 (вырождение ячейки).
							}
							else {
								if (octree1->linkW4->inum_FD == 0) {
									maxbound++;
								}
								else {
									// сосед существует.
									if (bvisit[octree1->linkW4->inum_FD - 1]) {
										// узел уже был посещён
										// Здесь неким образом модифицировалось gran_t
									}
									else {
										// узнать координаты центра грани и ориентацию в пространстве
										x_c = octree1->p0.x;
										y_c = 0.5*(octree1->linkW4->p1.y + octree1->linkW4->p2.y);
										z_c = 0.5*(octree1->linkW4->p1.z + octree1->linkW4->p5.z);
										iplane = YZ;
										patch_maxbound(iplane, s, ls, x_c, y_c, z_c, maxbound, binc3, lsid3);
									}
								}
							}

							if (octree1->linkW7 == NULL) {
								//maxbound++;
								// 21 сентября 2016 (вырождение ячейки).
							}
							else {
								if (octree1->linkW7->inum_FD == 0) {
									maxbound++;
								}
								else {
									// сосед существует.
									if (bvisit[octree1->linkW7->inum_FD - 1]) {
										// узел уже был посещён
										// Здесь неким образом модифицировалось gran_t
									}
									else {
										// узнать координаты центра грани и ориентацию в пространстве
										x_c = octree1->p0.x;
										y_c = 0.5*(octree1->linkW7->p1.y + octree1->linkW7->p2.y);
										z_c = 0.5*(octree1->linkW7->p1.z + octree1->linkW7->p5.z);
										iplane = YZ;
										patch_maxbound(iplane, s, ls, x_c, y_c, z_c, maxbound, binc4, lsid4);
									}
								}
							}
						}
					}


					// сосед не существует
					// записываем идентификатор грани и
					// увеличиваем счётчик граней.
					if (!octree1->b4N) {
						if (octree1->linkN == NULL) {
							maxbound++;
						}
						else {
							if (octree1->linkN->inum_FD == 0) {
								maxbound++;
							}
							else {
								// сосед существует.
								if (bvisit[octree1->linkN->inum_FD - 1]) {
									// узел уже был посещён
									// Здесь неким образом модифицировалось gran_t
								}
								else {
									bool binc;
									integer lsid;
									// узнать координаты центра грани и ориентацию в пространстве
									x_c = 0.5*(octree1->p2.x + octree1->p3.x);
									y_c = octree1->p3.y;
									z_c = 0.5*(octree1->p3.z + octree1->p7.z);
									iplane = XZ;
									patch_maxbound(iplane, s, ls, x_c, y_c, z_c, maxbound, binc, lsid);
									if (binc) {
										printf("error in calculate_max_bound_flow...\n");
										printf("model is incorrect. FATALL error!!!\n");
										printf("internal flat heat source can not be located\n");
										printf("inside the liquid.\n");
										printf("Please, press any key to exit...\n");
										//getchar();
										system("PAUSE");
										exit(1);
									}

								}
							}
						}
					}
					else {

						if ((octree1->linkN2 == NULL) && (octree1->linkN3 == NULL) && (octree1->linkN6 == NULL) && (octree1->linkN7 == NULL)) {
							// Зачем дробить если граница заканчивается NULL.
							maxbound++; // Вместо четырёх один.
						}
						else if ((octree1->linkN2 != NULL) && (octree1->linkN3 != NULL) && (octree1->linkN6 != NULL) && (octree1->linkN7 != NULL) &&
							(octree1->linkN2->inum_FD == 0) && (octree1->linkN3->inum_FD == 0) && (octree1->linkN6->inum_FD == 0) && (octree1->linkN7->inum_FD == 0)) {
							// Зачем дробить если граница заканчивается NULL.
							maxbound++; // Вместо четырёх один.
						}
						else if ((octree1->linkN2 == NULL) && (octree1->linkN3 != NULL) && (octree1->linkN6 == NULL) && (octree1->linkN7 != NULL) &&
							(octree1->linkN3->inum_FD == 0) && (octree1->linkN7->inum_FD == 0)) {
							// Зачем дробить если граница заканчивается NULL.
							maxbound++; // Вместо двух один.
						}
						else if ((octree1->linkN2 != NULL) && (octree1->linkN3 != NULL) && (octree1->linkN6 == NULL) && (octree1->linkN7 == NULL) &&
							(octree1->linkN2->inum_FD == 0) && (octree1->linkN3->inum_FD == 0)) {
							// Зачем дробить если граница заканчивается NULL.
							maxbound++; // Вместо двух один.
						}
						else if ((octree1->linkN2 == NULL) && (octree1->linkN3 != NULL) && (octree1->linkN6 == NULL) && (octree1->linkN7 == NULL) &&
							(octree1->linkN3->inum_FD == 0)) {
							// Зачем дробить если граница заканчивается NULL.
							maxbound++; // Вместо двух один.
						}
						else {

							bool binc1, binc2, binc3, binc4;
							integer lsid1, lsid2, lsid3, lsid4;

							if (octree1->linkN2 == NULL) {
								//maxbound++;
								// 21 сентября 2016 (вырождение ячейки).
							}
							else {
								if (octree1->linkN2->inum_FD == 0) {
									maxbound++;
								}
								else {
									// сосед существует.
									if (bvisit[octree1->linkN2->inum_FD - 1]) {
										// узел уже был посещён
										// Здесь неким образом модифицировалось gran_t
									}
									else {
										// узнать координаты центра грани и ориентацию в пространстве
										x_c = 0.5*(octree1->linkN2->p0.x + octree1->linkN2->p1.x);
										y_c = octree1->p3.y;
										z_c = 0.5*(octree1->linkN2->p0.z + octree1->linkN2->p4.z);
										iplane = XZ;
										patch_maxbound(iplane, s, ls, x_c, y_c, z_c, maxbound, binc1, lsid1);
									}
								}
							}

							if (octree1->linkN3 == NULL) {
								//maxbound++;
								printf("Error!!! N3 is NULL\n");
								//getchar();
								system("PAUSE");
							}
							else {
								if (octree1->linkN3->inum_FD == 0) {
									maxbound++;
								}
								else {
									// сосед существует.
									if (bvisit[octree1->linkN3->inum_FD - 1]) {
										// узел уже был посещён
										// Здесь неким образом модифицировалось gran_t
									}
									else {
										// узнать координаты центра грани и ориентацию в пространстве
										x_c = 0.5*(octree1->linkN3->p0.x + octree1->linkN3->p1.x);
										y_c = octree1->p3.y;
										z_c = 0.5*(octree1->linkN3->p0.z + octree1->linkN3->p4.z);
										iplane = XZ;
										patch_maxbound(iplane, s, ls, x_c, y_c, z_c, maxbound, binc2, lsid2);
									}
								}
							}

							if (octree1->linkN6 == NULL) {
								//maxbound++;
								// 21 сентября 2016 (вырождение ячейки).
							}
							else {
								if (octree1->linkN6->inum_FD == 0) {
									maxbound++;
								}
								else {
									// сосед существует.
									if (bvisit[octree1->linkN6->inum_FD - 1]) {
										// узел уже был посещён
										// Здесь неким образом модифицировалось gran_t
									}
									else {
										// узнать координаты центра грани и ориентацию в пространстве
										x_c = 0.5*(octree1->linkN6->p0.x + octree1->linkN6->p1.x);
										y_c = octree1->p3.y;
										z_c = 0.5*(octree1->linkN6->p0.z + octree1->linkN6->p4.z);
										iplane = XZ;
										patch_maxbound(iplane, s, ls, x_c, y_c, z_c, maxbound, binc3, lsid3);
									}
								}
							}

							if (octree1->linkN7 == NULL) {
								//maxbound++;
								// 21 сентября 2016 (вырождение ячейки).
							}
							else {
								if (octree1->linkN7->inum_FD == 0) {
									maxbound++;
								}
								else {
									// сосед существует.
									if (bvisit[octree1->linkN7->inum_FD - 1]) {
										// узел уже был посещён
										// Здесь неким образом модифицировалось gran_t
									}
									else {
										// узнать координаты центра грани и ориентацию в пространстве
										x_c = 0.5*(octree1->linkN7->p0.x + octree1->linkN7->p1.x);
										y_c = octree1->p3.y;
										z_c = 0.5*(octree1->linkN7->p0.z + octree1->linkN7->p4.z);
										iplane = XZ;
										patch_maxbound(iplane, s, ls, x_c, y_c, z_c, maxbound, binc4, lsid4);
									}
								}
							}
						}
					}

					// сосед не существует
					// записываем идентификатор грани и
					// увеличиваем счётчик граней.
					if (!octree1->b4S) {
						if (octree1->linkS == NULL) {
							maxbound++;
						}
						else {
							if (octree1->linkS->inum_FD == 0) {
								maxbound++;
							}
							else {
								// сосед существует.
								if (bvisit[octree1->linkS->inum_FD - 1]) {
									// узел уже был посещён
									// Здесь неким образом модифицировалось gran_t
								}
								else {
									bool binc;
									integer lsid;
									// узнать координаты центра грани и ориентацию в пространстве
									x_c = 0.5*(octree1->p0.x + octree1->p1.x);
									y_c = octree1->p0.y;
									z_c = 0.5*(octree1->p0.z + octree1->p4.z);
									iplane = XZ;
									patch_maxbound(iplane, s, ls, x_c, y_c, z_c, maxbound, binc, lsid);
									if (binc) {
										printf("error in calculate_max_bound_flow...\n");
										printf("model is incorrect. FATALL error!!!\n");
										printf("internal flat heat source can not be located\n");
										printf("inside the liquid.\n");
										printf("Please, press any key to exit...\n");
										//getchar();
										system("PAUSE");
										exit(1);
									}

								}
							}
						}
					}
					else {
						if ((octree1->linkS0 == NULL) && (octree1->linkS1 == NULL) && (octree1->linkS4 == NULL) && (octree1->linkS5 == NULL)) {
							// Зачем дробить если граница заканчивается NULL.
							maxbound++; // Вместо четырёх один.
						}
						else if ((octree1->linkS0 != NULL) && (octree1->linkS1 != NULL) && (octree1->linkS4 != NULL) && (octree1->linkS5 != NULL) &&
							(octree1->linkS0->inum_FD == 0) && (octree1->linkS1->inum_FD == 0) && (octree1->linkS4->inum_FD == 0) && (octree1->linkS5->inum_FD == 0)) {
							// Зачем дробить если граница заканчивается NULL.
							maxbound++; // Вместо четырёх один.
						}
						else if ((octree1->linkS0 != NULL) && (octree1->linkS1 == NULL) && (octree1->linkS4 != NULL) && (octree1->linkS5 == NULL) &&
							(octree1->linkS0->inum_FD == 0) && (octree1->linkS4->inum_FD == 0)) {
							// Зачем дробить если граница заканчивается NULL.
							maxbound++; // Вместо двух один.
						}
						else if ((octree1->linkS0 != NULL) && (octree1->linkS1 != NULL) && (octree1->linkS4 == NULL) && (octree1->linkS5 == NULL) &&
							(octree1->linkS0->inum_FD == 0) && (octree1->linkS1->inum_FD == 0)) {
							// Зачем дробить если граница заканчивается NULL.
							maxbound++; // Вместо двух один.
						}
						else if ((octree1->linkS0 != NULL) && (octree1->linkS1 == NULL) && (octree1->linkS4 == NULL) && (octree1->linkS5 == NULL) &&
							(octree1->linkS0->inum_FD == 0)) {
							// Зачем дробить если граница заканчивается NULL.
							maxbound++; // Вместо двух один.
						}
						else {

							bool binc1, binc2, binc3, binc4;
							integer lsid1, lsid2, lsid3, lsid4;

							if (octree1->linkS0 == NULL) {
								//maxbound++;
								printf("Error!!! S0 is NULL\n");
								//getchar();
								system("PAUSE");
							}
							else {
								if (octree1->linkS0->inum_FD == 0) {
									maxbound++;
								}
								else {
									// сосед существует.
									if (bvisit[octree1->linkS0->inum_FD - 1]) {
										// узел уже был посещён
										// Здесь неким образом модифицировалось gran_t
									}
									else {
										// узнать координаты центра грани и ориентацию в пространстве
										x_c = 0.5*(octree1->linkS0->p3.x + octree1->linkS0->p2.x);
										y_c = octree1->p0.y;
										z_c = 0.5*(octree1->linkS0->p3.z + octree1->linkS0->p7.z);
										iplane = XZ;
										patch_maxbound(iplane, s, ls, x_c, y_c, z_c, maxbound, binc1, lsid1);
									}
								}
							}

							if (octree1->linkS1 == NULL) {
								//maxbound++;
								// 21 сентября 2016 (вырождение ячейки).
							}
							else {
								if (octree1->linkS1->inum_FD == 0) {
									maxbound++;
								}
								else {
									// сосед существует.
									if (bvisit[octree1->linkS1->inum_FD - 1]) {
										// узел уже был посещён
										// Здесь неким образом модифицировалось gran_t
									}
									else {
										// узнать координаты центра грани и ориентацию в пространстве
										x_c = 0.5*(octree1->linkS1->p3.x + octree1->linkS1->p2.x);
										y_c = octree1->p0.y;
										z_c = 0.5*(octree1->linkS1->p3.z + octree1->linkS1->p7.z);
										iplane = XZ;
										patch_maxbound(iplane, s, ls, x_c, y_c, z_c, maxbound, binc2, lsid2);
									}
								}
							}

							if (octree1->linkS4 == NULL) {
								//maxbound++;
								// 21 сентября 2016 (вырождение ячейки).
							}
							else {
								if (octree1->linkS4->inum_FD == 0) {
									maxbound++;
								}
								else {
									// сосед существует.
									if (bvisit[octree1->linkS4->inum_FD - 1]) {
										// узел уже был посещён
										// Здесь неким образом модифицировалось gran_t
									}
									else {
										// узнать координаты центра грани и ориентацию в пространстве
										x_c = 0.5*(octree1->linkS4->p3.x + octree1->linkS4->p2.x);
										y_c = octree1->p0.y;
										z_c = 0.5*(octree1->linkS4->p3.z + octree1->linkS4->p7.z);
										iplane = XZ;
										patch_maxbound(iplane, s, ls, x_c, y_c, z_c, maxbound, binc3, lsid3);
									}
								}
							}

							if (octree1->linkS5 == NULL) {
								//maxbound++;
								// 21 сентября 2016 (вырождение ячейки).
							}
							else {
								if (octree1->linkS5->inum_FD == 0) {
									maxbound++;
								}
								else {
									// сосед существует.
									if (bvisit[octree1->linkS5->inum_FD - 1]) {
										// узел уже был посещён
										// Здесь неким образом модифицировалось gran_t
									}
									else {
										// узнать координаты центра грани и ориентацию в пространстве
										x_c = 0.5*(octree1->linkS5->p3.x + octree1->linkS5->p2.x);
										y_c = octree1->p0.y;
										z_c = 0.5*(octree1->linkS5->p3.z + octree1->linkS5->p7.z);
										iplane = XZ;
										patch_maxbound(iplane, s, ls, x_c, y_c, z_c, maxbound, binc4, lsid4);
									}
								}
							}
						}
					}

					// сосед не существует
					// записываем идентификатор грани и
					// увеличиваем счётчик граней.
					if (!octree1->b4T) {
						if (octree1->linkT == NULL) {
							maxbound++;
						}
						else {
							if (octree1->linkT->inum_FD == 0) {
								maxbound++;
							}
							else {
								// сосед существует.
								if (bvisit[octree1->linkT->inum_FD - 1]) {
									// узел уже был посещён
									// Здесь неким образом модифицировалось gran_t
								}
								else {
									bool binc;
									integer lsid;
									// узнать координаты центра грани и ориентацию в пространстве
									x_c = 0.5*(octree1->p4.x + octree1->p5.x);
									y_c = 0.5*(octree1->p4.y + octree1->p7.y);
									z_c = octree1->p4.z;
									iplane = XY;
									patch_maxbound(iplane, s, ls, x_c, y_c, z_c, maxbound, binc, lsid);
									if (binc) {
										printf("error in calculate_max_bound_flow...\n");
										printf("model is incorrect. FATALL error!!!\n");
										printf("internal flat heat source can not be located\n");
										printf("inside the liquid.\n");
										printf("Please, press any key to exit...\n");
										//getchar();
										system("PAUSE");
										exit(1);
									}

								}
							}
						}
					}
					else {

						if ((octree1->linkT4 == NULL) && (octree1->linkT5 == NULL) && (octree1->linkT6 == NULL) && (octree1->linkT7 == NULL)) {
							// Зачем дробить если граница заканчивается NULL.
							maxbound++; // Вместо четырёх один.
						}
						else if ((octree1->linkT4 != NULL) && (octree1->linkT5 != NULL) && (octree1->linkT6 != NULL) && (octree1->linkT7 != NULL) &&
							(octree1->linkT4->inum_FD == 0) && (octree1->linkT5->inum_FD == 0) && (octree1->linkT6->inum_FD == 0) && (octree1->linkT7->inum_FD == 0)) {
							// Зачем дробить если граница заканчивается NULL.
							maxbound++; // Вместо четырёх один.
						}
						else if ((octree1->linkT4 != NULL) && (octree1->linkT5 == NULL) && (octree1->linkT6 == NULL) && (octree1->linkT7 != NULL) &&
							(octree1->linkT4->inum_FD == 0) && (octree1->linkT7->inum_FD == 0)) {
							// Зачем дробить если граница заканчивается NULL.
							maxbound++; // Вместо двух один.
						}
						else if ((octree1->linkT4 != NULL) && (octree1->linkT5 != NULL) && (octree1->linkT6 == NULL) && (octree1->linkT7 == NULL) &&
							(octree1->linkT4->inum_FD == 0) && (octree1->linkT5->inum_FD == 0)) {
							// Зачем дробить если граница заканчивается NULL.
							maxbound++; // Вместо двух один.
						}
						else if ((octree1->linkT4 != NULL) && (octree1->linkT5 == NULL) && (octree1->linkT6 == NULL) && (octree1->linkT7 == NULL) &&
							(octree1->linkT4->inum_FD == 0)) {
							// Зачем дробить если граница заканчивается NULL.
							maxbound++; // Вместо двух один.
						}
						else {

							bool binc1, binc2, binc3, binc4;
							integer lsid1, lsid2, lsid3, lsid4;

							if (octree1->linkT4 == NULL) {
								//maxbound++;
								printf("Error!!! T4 is NULL\n");
								//getchar();
								system("PAUSE");
							}
							else {
								if (octree1->linkT4->inum_FD == 0) {
									maxbound++;
								}
								else {
									// сосед существует.
									if (bvisit[octree1->linkT4->inum_FD - 1]) {
										// узел уже был посещён
										// Здесь неким образом модифицировалось gran_t
									}
									else {
										// узнать координаты центра грани и ориентацию в пространстве
										x_c = 0.5*(octree1->linkT4->p0.x + octree1->linkT4->p1.x);
										y_c = 0.5*(octree1->linkT4->p0.y + octree1->linkT4->p3.y);
										z_c = octree1->p4.z;
										iplane = XY;
										patch_maxbound(iplane, s, ls, x_c, y_c, z_c, maxbound, binc1, lsid1);
									}
								}
							}

							if (octree1->linkT5 == NULL) {
								//maxbound++;
								// 21 сентября 2016 (вырождение ячейки).
							}
							else {
								if (octree1->linkT5->inum_FD == 0) {
									maxbound++;
								}
								else {
									// сосед существует.
									if (bvisit[octree1->linkT5->inum_FD - 1]) {
										// узел уже был посещён
										// Здесь неким образом модифицировалось gran_t
									}
									else {
										// узнать координаты центра грани и ориентацию в пространстве
										x_c = 0.5*(octree1->linkT5->p0.x + octree1->linkT5->p1.x);
										y_c = 0.5*(octree1->linkT5->p0.y + octree1->linkT5->p3.y);
										z_c = octree1->p4.z;
										iplane = XY;
										patch_maxbound(iplane, s, ls, x_c, y_c, z_c, maxbound, binc2, lsid2);
									}
								}
							}

							if (octree1->linkT6 == NULL) {
								//maxbound++;
								// 21 сентября 2016 (вырождение ячейки).
							}
							else {
								if (octree1->linkT6->inum_FD == 0) {
									maxbound++;
								}
								else {
									// сосед существует.
									if (bvisit[octree1->linkT6->inum_FD - 1]) {
										// узел уже был посещён
										// Здесь неким образом модифицировалось gran_t
									}
									else {
										// узнать координаты центра грани и ориентацию в пространстве
										x_c = 0.5*(octree1->linkT6->p0.x + octree1->linkT6->p1.x);
										y_c = 0.5*(octree1->linkT6->p0.y + octree1->linkT6->p3.y);
										z_c = octree1->p4.z;
										iplane = XY;
										patch_maxbound(iplane, s, ls, x_c, y_c, z_c, maxbound, binc3, lsid3);
									}
								}
							}

							if (octree1->linkT7 == NULL) {
								//maxbound++;
								// 21 сентября 2016 (вырождение ячейки).
							}
							else {
								if (octree1->linkT7->inum_FD == 0) {
									maxbound++;
								}
								else {
									// сосед существует.
									if (bvisit[octree1->linkT7->inum_FD - 1]) {
										// узел уже был посещён
										// Здесь неким образом модифицировалось gran_t
									}
									else {
										// узнать координаты центра грани и ориентацию в пространстве
										x_c = 0.5*(octree1->linkT7->p0.x + octree1->linkT7->p1.x);
										y_c = 0.5*(octree1->linkT7->p0.y + octree1->linkT7->p3.y);
										z_c = octree1->p4.z;
										iplane = XY;
										patch_maxbound(iplane, s, ls, x_c, y_c, z_c, maxbound, binc4, lsid4);
									}
								}
							}
						}
					}

					// сосед не существует
					// записываем идентификатор грани и
					// увеличиваем счётчик граней.
					if (!octree1->b4B) {
						if (octree1->linkB == NULL) {
							maxbound++;
						}
						else {
							if (octree1->linkB->inum_FD == 0) {
								maxbound++;
							}
							else {
								// сосед существует.
								if (bvisit[octree1->linkB->inum_FD - 1]) {
									// узел уже был посещён
									// Здесь неким образом модифицировалось gran_t
								}
								else {
									bool binc;
									integer lsid;
									// узнать координаты центра грани и ориентацию в пространстве
									x_c = 0.5*(octree1->p0.x + octree1->p1.x);
									y_c = 0.5*(octree1->p0.y + octree1->p3.y);
									z_c = octree1->p0.z;
									iplane = XY;
									patch_maxbound(iplane, s, ls, x_c, y_c, z_c, maxbound, binc, lsid);
									if (binc) {
										printf("error in calculate_max_bound_flow...\n");
										printf("model is incorrect. FATALL error!!!\n");
										printf("internal flat heat source can not be located\n");
										printf("inside the liquid.\n");
										printf("Please, press any key to exit...\n");
										//getchar();
										system("PAUSE");
										exit(1);
									}

								}
							}
						}
					}
					else {

						if ((octree1->linkB0 == NULL) && (octree1->linkB1 == NULL) && (octree1->linkB2 == NULL) && (octree1->linkB3 == NULL)) {
							// Зачем дробить если граница заканчивается NULL.
							maxbound++; // Вместо четырёх один.
						}
						else if ((octree1->linkB0 != NULL) && (octree1->linkB1 != NULL) && (octree1->linkB2 != NULL) && (octree1->linkB3 != NULL) &&
							(octree1->linkB0->inum_FD == 0) && (octree1->linkB1->inum_FD == 0) && (octree1->linkB2->inum_FD == 0) && (octree1->linkB3->inum_FD == 0)) {
							// Зачем дробить если граница заканчивается NULL.
							maxbound++; // Вместо четырёх один.
						}
						else if ((octree1->linkB0 != NULL) && (octree1->linkB1 == NULL) && (octree1->linkB2 == NULL) && (octree1->linkB3 != NULL) &&
							(octree1->linkB0->inum_FD == 0) && (octree1->linkB3->inum_FD == 0)) {
							// Зачем дробить если граница заканчивается NULL.
							maxbound++; // Вместо двух один.
						}
						else if ((octree1->linkB0 != NULL) && (octree1->linkB1 != NULL) && (octree1->linkB2 == NULL) && (octree1->linkB3 == NULL) &&
							(octree1->linkB0->inum_FD == 0) && (octree1->linkB1->inum_FD == 0)) {
							// Зачем дробить если граница заканчивается NULL.
							maxbound++; // Вместо двух один.
						}
						else if ((octree1->linkB0 != NULL) && (octree1->linkB1 == NULL) && (octree1->linkB2 == NULL) && (octree1->linkB3 == NULL) &&
							(octree1->linkB0->inum_FD == 0)) {
							// Зачем дробить если граница заканчивается NULL.
							maxbound++; // Вместо двух один.
						}
						else {

							bool binc1, binc2, binc3, binc4;
							integer lsid1, lsid2, lsid3, lsid4;

							if (octree1->linkB0 == NULL) {
								//maxbound++;
								printf("Error!!! B0 is NULL\n");
								//getchar();
								system("PAUSE");
							}
							else {
								if (octree1->linkB0->inum_FD == 0) {
									maxbound++;
								}
								else {
									// сосед существует.
									if (bvisit[octree1->linkB0->inum_FD - 1]) {
										// узел уже был посещён
										// Здесь неким образом модифицировалось gran_t
									}
									else {
										// узнать координаты центра грани и ориентацию в пространстве
										x_c = 0.5*(octree1->linkB0->p4.x + octree1->linkB0->p5.x);
										y_c = 0.5*(octree1->linkB0->p4.y + octree1->linkB0->p7.y);
										z_c = octree1->p0.z;
										iplane = XY;
										patch_maxbound(iplane, s, ls, x_c, y_c, z_c, maxbound, binc1, lsid1);
									}
								}
							}

							if (octree1->linkB1 == NULL) {
								//maxbound++;
								// 21 сентября 2016 (вырождение ячейки).
							}
							else {
								if (octree1->linkB1->inum_FD == 0) {
									maxbound++;
								}
								else {
									// сосед существует.
									if (bvisit[octree1->linkB1->inum_FD - 1]) {
										// узел уже был посещён
										// Здесь неким образом модифицировалось gran_t
									}
									else {
										// узнать координаты центра грани и ориентацию в пространстве
										x_c = 0.5*(octree1->linkB1->p4.x + octree1->linkB1->p5.x);
										y_c = 0.5*(octree1->linkB1->p4.y + octree1->linkB1->p7.y);
										z_c = octree1->p0.z;
										iplane = XY;
										patch_maxbound(iplane, s, ls, x_c, y_c, z_c, maxbound, binc2, lsid2);
									}
								}
							}

							if (octree1->linkB2 == NULL) {
								//maxbound++;
								// 21 сентября 2016 (вырождение ячейки).
							}
							else {
								if (octree1->linkB2->inum_FD == 0) {
									maxbound++;
								}
								else {
									// сосед существует.
									if (bvisit[octree1->linkB2->inum_FD - 1]) {
										// узел уже был посещён
										// Здесь неким образом модифицировалось gran_t
									}
									else {
										// узнать координаты центра грани и ориентацию в пространстве
										x_c = 0.5*(octree1->linkB2->p4.x + octree1->linkB2->p5.x);
										y_c = 0.5*(octree1->linkB2->p4.y + octree1->linkB2->p7.y);
										z_c = octree1->p0.z;
										iplane = XY;
										patch_maxbound(iplane, s, ls, x_c, y_c, z_c, maxbound, binc3, lsid3);
									}
								}
							}

							if (octree1->linkB3 == NULL) {
								//maxbound++;
								// 21 сентября 2016 (вырождение ячейки).
							}
							else {
								if (octree1->linkB3->inum_FD == 0) {
									maxbound++;
								}
								else {
									// сосед существует.
									if (bvisit[octree1->linkB3->inum_FD - 1]) {
										// узел уже был посещён
										// Здесь неким образом модифицировалось gran_t
									}
									else {
										// узнать координаты центра грани и ориентацию в пространстве
										x_c = 0.5*(octree1->linkB3->p4.x + octree1->linkB3->p5.x);
										y_c = 0.5*(octree1->linkB3->p4.y + octree1->linkB3->p7.y);
										z_c = octree1->p0.z;
										iplane = XY;
										patch_maxbound(iplane, s, ls, x_c, y_c, z_c, maxbound, binc4, lsid4);
									}
								}
							}
						}
					}


					bvisit[maxelm] = true; // Узел был посещён.
					maxelm++;
				}
				//octree1->inum_FD = 0; // По умолчанию не принадлежит расчётной области.
				octree1 = NULL;
				my_ALICE_STACK[top_ALICE_STACK - 1].link = NULL;
				top_ALICE_STACK--;
			}
			else {
				// продолжаем добираться до листьев.
				STACK_ALICE buf1 = my_ALICE_STACK[top_ALICE_STACK - 1];
				STACK_ALICE* buf = &buf1;
				top_ALICE_STACK--;
				if (buf->link->link0 != NULL) {
					my_ALICE_STACK[top_ALICE_STACK].link = (buf->link->link0);
					my_ALICE_STACK[top_ALICE_STACK].maxx = buf->link->link0->maxx;
					my_ALICE_STACK[top_ALICE_STACK].minx = buf->link->link0->minx;
					my_ALICE_STACK[top_ALICE_STACK].maxy = buf->link->link0->maxy;
					my_ALICE_STACK[top_ALICE_STACK].miny = buf->link->link0->miny;
					my_ALICE_STACK[top_ALICE_STACK].maxz = buf->link->link0->maxz;
					my_ALICE_STACK[top_ALICE_STACK].minz = buf->link->link0->minz;
					top_ALICE_STACK++;
				}
				if (buf->link->link1 != NULL) {
					my_ALICE_STACK[top_ALICE_STACK].link = (buf->link->link1);
					my_ALICE_STACK[top_ALICE_STACK].maxx = buf->link->link1->maxx;
					my_ALICE_STACK[top_ALICE_STACK].minx = buf->link->link1->minx;
					my_ALICE_STACK[top_ALICE_STACK].maxy = buf->link->link1->maxy;
					my_ALICE_STACK[top_ALICE_STACK].miny = buf->link->link1->miny;
					my_ALICE_STACK[top_ALICE_STACK].maxz = buf->link->link1->maxz;
					my_ALICE_STACK[top_ALICE_STACK].minz = buf->link->link1->minz;
					top_ALICE_STACK++;
				}
				if (buf->link->link2 != NULL) {
					my_ALICE_STACK[top_ALICE_STACK].link = (buf->link->link2);
					my_ALICE_STACK[top_ALICE_STACK].maxx = buf->link->link2->maxx;
					my_ALICE_STACK[top_ALICE_STACK].minx = buf->link->link2->minx;
					my_ALICE_STACK[top_ALICE_STACK].maxy = buf->link->link2->maxy;
					my_ALICE_STACK[top_ALICE_STACK].miny = buf->link->link2->miny;
					my_ALICE_STACK[top_ALICE_STACK].maxz = buf->link->link2->maxz;
					my_ALICE_STACK[top_ALICE_STACK].minz = buf->link->link2->minz;
					top_ALICE_STACK++;
				}
				if (buf->link->link3 != NULL) {
					my_ALICE_STACK[top_ALICE_STACK].link = (buf->link->link3);
					my_ALICE_STACK[top_ALICE_STACK].maxx = buf->link->link3->maxx;
					my_ALICE_STACK[top_ALICE_STACK].minx = buf->link->link3->minx;
					my_ALICE_STACK[top_ALICE_STACK].maxy = buf->link->link3->maxy;
					my_ALICE_STACK[top_ALICE_STACK].miny = buf->link->link3->miny;
					my_ALICE_STACK[top_ALICE_STACK].maxz = buf->link->link3->maxz;
					my_ALICE_STACK[top_ALICE_STACK].minz = buf->link->link3->minz;
					top_ALICE_STACK++;
				}
				if (buf->link->link4 != NULL) {
					my_ALICE_STACK[top_ALICE_STACK].link = (buf->link->link4);
					my_ALICE_STACK[top_ALICE_STACK].maxx = buf->link->link4->maxx;
					my_ALICE_STACK[top_ALICE_STACK].minx = buf->link->link4->minx;
					my_ALICE_STACK[top_ALICE_STACK].maxy = buf->link->link4->maxy;
					my_ALICE_STACK[top_ALICE_STACK].miny = buf->link->link4->miny;
					my_ALICE_STACK[top_ALICE_STACK].maxz = buf->link->link4->maxz;
					my_ALICE_STACK[top_ALICE_STACK].minz = buf->link->link4->minz;
					top_ALICE_STACK++;
				}
				if (buf->link->link5 != NULL) {
					my_ALICE_STACK[top_ALICE_STACK].link = (buf->link->link5);
					my_ALICE_STACK[top_ALICE_STACK].maxx = buf->link->link5->maxx;
					my_ALICE_STACK[top_ALICE_STACK].minx = buf->link->link5->minx;
					my_ALICE_STACK[top_ALICE_STACK].maxy = buf->link->link5->maxy;
					my_ALICE_STACK[top_ALICE_STACK].miny = buf->link->link5->miny;
					my_ALICE_STACK[top_ALICE_STACK].maxz = buf->link->link5->maxz;
					my_ALICE_STACK[top_ALICE_STACK].minz = buf->link->link5->minz;
					top_ALICE_STACK++;
				}
				if (buf->link->link6 != NULL) {
					my_ALICE_STACK[top_ALICE_STACK].link = (buf->link->link6);
					my_ALICE_STACK[top_ALICE_STACK].maxx = buf->link->link6->maxx;
					my_ALICE_STACK[top_ALICE_STACK].minx = buf->link->link6->minx;
					my_ALICE_STACK[top_ALICE_STACK].maxy = buf->link->link6->maxy;
					my_ALICE_STACK[top_ALICE_STACK].miny = buf->link->link6->miny;
					my_ALICE_STACK[top_ALICE_STACK].maxz = buf->link->link6->maxz;
					my_ALICE_STACK[top_ALICE_STACK].minz = buf->link->link6->minz;
					top_ALICE_STACK++;
				}
				if (buf->link->link7 != NULL) {
					my_ALICE_STACK[top_ALICE_STACK].link = (buf->link->link7);
					my_ALICE_STACK[top_ALICE_STACK].maxx = buf->link->link7->maxx;
					my_ALICE_STACK[top_ALICE_STACK].minx = buf->link->link7->minx;
					my_ALICE_STACK[top_ALICE_STACK].maxy = buf->link->link7->maxy;
					my_ALICE_STACK[top_ALICE_STACK].miny = buf->link->link7->miny;
					my_ALICE_STACK[top_ALICE_STACK].maxz = buf->link->link7->maxz;
					my_ALICE_STACK[top_ALICE_STACK].minz = buf->link->link7->minz;
					top_ALICE_STACK++;
				}
			}
		}
		//}
		//getchar();
	}

	delete[] bvisit;
	bvisit = NULL;
} // calculate_max_bound_flow


// G - грань по которой ставится граничное условие.
// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
// bound_id - номер граничного узла, начиная с нуля.
// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
void obrabotka_granichnoi_grani(integer G, BOUND* &sosedb, integer **nvtx, TOCHKA* pa, integer elm_id, integer bound_id,
	integer elm_id_inverse, integer maxelm_out, integer* &whot_is_block, integer ls, integer lw, WALL* w, SOURCE* s, BLOCK* b,
	doublereal dS, bool binternalg, bool bvisit_g)
{
	// G - грань по которой ставится граничное условие.
	// elm_id - номер конечного элемента, начиная с нуля для которого рассматривается граничный узел.
	// bound_id - номер граничного узла, начиная с нуля.
	// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.

	doublereal x_c = 0.0, y_c = 0.0, z_c = 0.0; // координаты центра грани
	integer iplane; // плоскость в которой лежит грань.

	// Эта информация совсем неактивна при данной сборке.
	for (integer j = 0; j<6; j++) sosedb[bound_id].iW[j] = NON_EXISTENT_NODE; // инициализация

	// узнать координаты центра грани и ориентацию в пространстве
	switch (G) {
	case ESIDE: x_c = pa[nvtx[1][elm_id] - 1].x;
		y_c = 0.5*(pa[nvtx[1][elm_id] - 1].y + pa[nvtx[3][elm_id] - 1].y);
		z_c = 0.5*(pa[nvtx[0][elm_id] - 1].z + pa[nvtx[4][elm_id] - 1].z);
		iplane = YZ;
		break;
	case WSIDE: x_c = pa[nvtx[0][elm_id] - 1].x;
		y_c = 0.5*(pa[nvtx[1][elm_id] - 1].y + pa[nvtx[3][elm_id] - 1].y);
		z_c = 0.5*(pa[nvtx[0][elm_id] - 1].z + pa[nvtx[4][elm_id] - 1].z);
		iplane = YZ;
		break;
	case NSIDE: x_c = 0.5*(pa[nvtx[0][elm_id] - 1].x + pa[nvtx[1][elm_id] - 1].x);
		y_c = pa[nvtx[2][elm_id] - 1].y;
		z_c = 0.5*(pa[nvtx[0][elm_id] - 1].z + pa[nvtx[4][elm_id] - 1].z);
		iplane = XZ;
		break;
	case SSIDE: x_c = 0.5*(pa[nvtx[0][elm_id] - 1].x + pa[nvtx[1][elm_id] - 1].x);
		y_c = pa[nvtx[0][elm_id] - 1].y;
		z_c = 0.5*(pa[nvtx[0][elm_id] - 1].z + pa[nvtx[4][elm_id] - 1].z);
		iplane = XZ;
		break;
	case TSIDE: x_c = 0.5*(pa[nvtx[0][elm_id] - 1].x + pa[nvtx[1][elm_id] - 1].x);
		y_c = 0.5*(pa[nvtx[1][elm_id] - 1].y + pa[nvtx[3][elm_id] - 1].y);
		z_c = pa[nvtx[4][elm_id] - 1].z;
		iplane = XY;
		break;
	case BSIDE: x_c = 0.5*(pa[nvtx[0][elm_id] - 1].x + pa[nvtx[1][elm_id] - 1].x);
		y_c = 0.5*(pa[nvtx[1][elm_id] - 1].y + pa[nvtx[3][elm_id] - 1].y);
		z_c = pa[nvtx[0][elm_id] - 1].z;
		iplane = XY;
		break;
	} // end case

	if ((fabs(x_c) > 1.0e20) || (fabs(y_c) > 1.0e20) || (fabs(z_c) > 1.0e20)) 
	{
		printf("ERROR in function obrabotka_granichnoi_grani in module constr_struct_alice.cpp\n");
		printf("if ((fabs(x_c) > 1.0e20) || (fabs(y_c) > 1.0e20) || (fabs(z_c) > 1.0e20)) \n");
#if doubleintprecision == 1
		printf("elm_id=%lld maxelm=%lld bound_id=%lld\n", elm_id, maxelm_out, bound_id);
#else
		printf("elm_id=%d maxelm=%d bound_id=%d\n", elm_id, maxelm_out, bound_id);
#endif		
		//getchar();
		system("PAUSE");
	}

	// Площадь грани передаётся извне.
	sosedb[bound_id].dS = dS; // Присваиваем площадь граничной грани.

	// грань лежит на границе расчётной области.
	switch (G) {
	case ESIDE: sosedb[bound_id].Norm = WSIDE;
		break;
	case WSIDE: sosedb[bound_id].Norm = ESIDE;
		break;
	case NSIDE: sosedb[bound_id].Norm = SSIDE;
		break;
	case SSIDE: sosedb[bound_id].Norm = NSIDE;
		break;
	case TSIDE: sosedb[bound_id].Norm = BSIDE;
		break;
	case BSIDE: sosedb[bound_id].Norm = TSIDE;	
		break;
	} // end определеиие внутренней нормали
	sosedb[bound_id].iII = elm_id_inverse;
	sosedb[bound_id].iB = maxelm_out + bound_id;
	sosedb[bound_id].iI = elm_id;
	
	if (binternalg) {
		// 1 октября 2016 (суббота).

		// Управляющая логика здесь не нужна т.к. она уже была выполнена при вызове obrabotka_granichnoi_grani.
		//TOCHKA p;
		//p.x = 0.125*(pa[nvtx[0][elm_id] - 1].x + pa[nvtx[1][elm_id] - 1].x + pa[nvtx[2][elm_id] - 1].x + pa[nvtx[3][elm_id] - 1].x + pa[nvtx[4][elm_id] - 1].x + pa[nvtx[5][elm_id] - 1].x + pa[nvtx[6][elm_id] - 1].x + pa[nvtx[7][elm_id] - 1].x);
		//p.y = 0.125*(pa[nvtx[0][elm_id] - 1].y + pa[nvtx[1][elm_id] - 1].y + pa[nvtx[2][elm_id] - 1].y + pa[nvtx[3][elm_id] - 1].y + pa[nvtx[4][elm_id] - 1].y + pa[nvtx[5][elm_id] - 1].y + pa[nvtx[6][elm_id] - 1].y + pa[nvtx[7][elm_id] - 1].y);
		//p.z = 0.125*(pa[nvtx[0][elm_id] - 1].z + pa[nvtx[1][elm_id] - 1].z + pa[nvtx[2][elm_id] - 1].z + pa[nvtx[3][elm_id] - 1].z + pa[nvtx[4][elm_id] - 1].z + pa[nvtx[5][elm_id] - 1].z + pa[nvtx[6][elm_id] - 1].z + pa[nvtx[7][elm_id] - 1].z);
		//integer ib;
		//if (in_model_flow(p, ib, b, lb)) {
			//sosedb[bound_id].iI1 = elm_id;
		//} else {sosedb[bound_id].iI1 = elm_id; }
		// В любом случае мы присваиваем iI позицию узла.
		if (!bvisit_g) {
			// грань еще не была посещена.
			sosedb[bound_id].iI1 = elm_id;
		}
		else {
			// грань уже была посещена.
			sosedb[bound_id].iI2 = elm_id;
		}
	}
	else {
		// Не внутренний источник тепла.
		sosedb[bound_id].iI1 = NON_EXISTENT_NODE;
		sosedb[bound_id].iI2 = NON_EXISTENT_NODE;
	}
	//bvisit[bound_id] = true; // грань была посещена
	// Вычисляем emissivity:
	//integer ibfound = NON_EXISTENT_NODE;
	integer ibfound = whot_is_block[elm_id];
	if (ibfound >= 0) {
		// блок найден.
		// Определяем внутреннюю нормаль:
		switch (sosedb[bound_id].Norm) {
		case WSIDE: sosedb[bound_id].emissivity = b[ibfound].radiation.emissE;
			break;
		case ESIDE: sosedb[bound_id].emissivity = b[ibfound].radiation.emissW;
			break;
		case SSIDE: sosedb[bound_id].emissivity = b[ibfound].radiation.emissN;
			break;
		case NSIDE: sosedb[bound_id].emissivity = b[ibfound].radiation.emissS;
			break;
		case BSIDE: sosedb[bound_id].emissivity = b[ibfound].radiation.emissT;
			break;
		case TSIDE: sosedb[bound_id].emissivity = b[ibfound].radiation.emissB;
			break;
		}
	}
	else {
		printf("error: emissivity calculation block unfound. in constr_sosedb_temp.\n");
		system("PAUSE");
		exit(1);
	}
	integer jpos;

	bool bfind = false;
	sosedb[bound_id].MCB = ls + lw; // Инициализация.
	for (integer j = 0; j<ls; j++) {
		if (s[j].iPlane == iplane) {
			// Важно не только попадание в фокус объекта но и нахождение с объектом на одном уровне:
			switch (iplane) {
			case XY: if ((x_c>s[j].g.xS) && (x_c<s[j].g.xE) && (y_c>s[j].g.yS) && (y_c<s[j].g.yE) && (fabs(z_c - s[j].g.zE)<admission)) { bfind = true; jpos = j; } break;
			case YZ: if ((z_c>s[j].g.zS) && (z_c<s[j].g.zE) && (y_c>s[j].g.yS) && (y_c<s[j].g.yE) && (fabs(x_c - s[j].g.xE)<admission)) { bfind = true; jpos = j; } break;
			case XZ: if ((x_c>s[j].g.xS) && (x_c<s[j].g.xE) && (z_c>s[j].g.zS) && (z_c<s[j].g.zE) && (fabs(y_c - s[j].g.yE)<admission)) { bfind = true; jpos = j; } break;
			}
		}
	}
	if (bfind) {
		// внешний источник тепла
		//printf("source out found...\n"); // debug
		//getchar();
		sosedb[bound_id].MCB = jpos;
	}
	else {
		for (integer j = 0; j<lw; j++) {
			if (w[j].iPlane == iplane) {
				// Важно не только попадание в фокус объекта но и нахождение с объектом на одном уровне:
				switch (iplane) {
				case XY: if ((x_c>w[j].g.xS) && (x_c<w[j].g.xE) && (y_c>w[j].g.yS) && (y_c<w[j].g.yE) && (fabs(z_c - w[j].g.zE)<admission)) { bfind = true; jpos = j; } break;
				case YZ: if ((z_c>w[j].g.zS) && (z_c<w[j].g.zE) && (y_c>w[j].g.yS) && (y_c<w[j].g.yE) && (fabs(x_c - w[j].g.xE)<admission)) { bfind = true; jpos = j; } break;
				case XZ: if ((x_c>w[j].g.xS) && (x_c<w[j].g.xE) && (z_c>w[j].g.zS) && (z_c<w[j].g.zE) && (fabs(y_c - w[j].g.yE)<admission)) { bfind = true; jpos = j; } break;
				}
			}
		}
		if (bfind) {
			// найдена стека
			//printf("ambient wall found...\n"); // debug
			//getchar();
			sosedb[bound_id].MCB = ls + jpos;
		}
		else {
			// граничное условие по умолчанию
			sosedb[bound_id].MCB = ls + lw;
		}

	}


} // obrabotka_granichnoi_grani

  // G - грань по которой ставится граничное условие.
  // elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
  // bound_id - номер граничного узла, начиная с нуля.
  // elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
void obrabotka_granichnoi_grani(integer G, BOUND* &sosedb, integer **nvtx, TOCHKA* pa, integer elm_id, integer bound_id,
	integer elm_id_inverse, integer maxelm_out, integer* &whot_is_block, integer ls, integer lw, WALL* w, SOURCE* s, BLOCK* b,
	doublereal dS, bool binternalg, bool bvisit_g, TOCHKA &p_centerG)
{
	// G - грань по которой ставится граничное условие.
	// elm_id - номер конечного элемента, начиная с нуля для которого рассматривается граничный узел.
	// bound_id - номер граничного узла, начиная с нуля.
	// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
	// p_centerG - декартовы координаты центра грани ячейки.

	
	integer iplane; // плоскость в которой лежит грань.

					// Эта информация совсем неактивна при данной сборке.
	for (integer j = 0; j<6; j++) sosedb[bound_id].iW[j] = NON_EXISTENT_NODE; // инициализация

																			  // узнать координаты центра грани и ориентацию в пространстве
	switch (G) {
	case ESIDE: 
		iplane = YZ;
		break;
	case WSIDE: 
		iplane = YZ;
		break;
	case NSIDE: 
		iplane = XZ;
		break;
	case SSIDE: 
		iplane = XZ;
		break;
	case TSIDE: 
		iplane = XY;
		break;
	case BSIDE: 
		iplane = XY;
		break;
	} // end case

	if ((fabs(p_centerG.x) > 1.0e20) || (fabs(p_centerG.y) > 1.0e20) || (fabs(p_centerG.z) > 1.0e20))
	{
		printf("ERROR in function obrabotka_granichnoi_grani in module constr_struct_alice.cpp\n");
		printf("if ((fabs(p_centerG.x) > 1.0e20) || (fabs(p_centerG.y) > 1.0e20) || (fabs(p_centerG.z) > 1.0e20)) \n");
#if doubleintprecision == 1
		printf("elm_id=%lld maxelm=%lld bound_id=%lld\n", elm_id, maxelm_out, bound_id);
#else
		printf("elm_id=%d maxelm=%d bound_id=%d\n", elm_id, maxelm_out, bound_id);
#endif		
		//getchar();
		system("PAUSE");
	}

	// Площадь грани передаётся извне.
	sosedb[bound_id].dS = dS; // Присваиваем площадь граничной грани.
	sosedb[bound_id].p_c.x = p_centerG.x;
	sosedb[bound_id].p_c.y = p_centerG.y;
	sosedb[bound_id].p_c.z = p_centerG.z;

							  // грань лежит на границе расчётной области.
	switch (G) {
	case ESIDE: sosedb[bound_id].Norm = WSIDE;
		break;
	case WSIDE: sosedb[bound_id].Norm = ESIDE;
		break;
	case NSIDE: sosedb[bound_id].Norm = SSIDE;
		break;
	case SSIDE: sosedb[bound_id].Norm = NSIDE;
		break;
	case TSIDE: sosedb[bound_id].Norm = BSIDE;
		break;
	case BSIDE: sosedb[bound_id].Norm = TSIDE;
		break;
	} // end определеиие внутренней нормали
	sosedb[bound_id].iII = elm_id_inverse;
	sosedb[bound_id].iB = maxelm_out + bound_id;
	sosedb[bound_id].iI = elm_id;

	if (binternalg) {
		// 1 октября 2016 (суббота).

		// Управляющая логика здесь не нужна т.к. она уже была выполнена при вызове obrabotka_granichnoi_grani.
		//TOCHKA p;
		//p.x = 0.125*(pa[nvtx[0][elm_id] - 1].x + pa[nvtx[1][elm_id] - 1].x + pa[nvtx[2][elm_id] - 1].x + pa[nvtx[3][elm_id] - 1].x + pa[nvtx[4][elm_id] - 1].x + pa[nvtx[5][elm_id] - 1].x + pa[nvtx[6][elm_id] - 1].x + pa[nvtx[7][elm_id] - 1].x);
		//p.y = 0.125*(pa[nvtx[0][elm_id] - 1].y + pa[nvtx[1][elm_id] - 1].y + pa[nvtx[2][elm_id] - 1].y + pa[nvtx[3][elm_id] - 1].y + pa[nvtx[4][elm_id] - 1].y + pa[nvtx[5][elm_id] - 1].y + pa[nvtx[6][elm_id] - 1].y + pa[nvtx[7][elm_id] - 1].y);
		//p.z = 0.125*(pa[nvtx[0][elm_id] - 1].z + pa[nvtx[1][elm_id] - 1].z + pa[nvtx[2][elm_id] - 1].z + pa[nvtx[3][elm_id] - 1].z + pa[nvtx[4][elm_id] - 1].z + pa[nvtx[5][elm_id] - 1].z + pa[nvtx[6][elm_id] - 1].z + pa[nvtx[7][elm_id] - 1].z);
		//integer ib;
		//if (in_model_flow(p, ib, b, lb)) {
		//sosedb[bound_id].iI1 = elm_id;
		//} else {sosedb[bound_id].iI1 = elm_id; }
		// В любом случае мы присваиваем iI позицию узла.
		if (!bvisit_g) {
			// грань еще не была посещена.
			sosedb[bound_id].iI1 = elm_id;
		}
		else {
			// грань уже была посещена.
			sosedb[bound_id].iI2 = elm_id;
		}
	}
	else {
		// Не внутренний источник тепла.
		sosedb[bound_id].iI1 = NON_EXISTENT_NODE;
		sosedb[bound_id].iI2 = NON_EXISTENT_NODE;
	}
	//bvisit[bound_id] = true; // грань была посещена
	// Вычисляем emissivity:
	//integer ibfound = NON_EXISTENT_NODE;
	integer ibfound = whot_is_block[elm_id];
	if (ibfound >= 0) {
		// блок найден.
		// Определяем внутреннюю нормаль:
		switch (sosedb[bound_id].Norm) {
		case WSIDE: sosedb[bound_id].emissivity = b[ibfound].radiation.emissE;
			break;
		case ESIDE: sosedb[bound_id].emissivity = b[ibfound].radiation.emissW;
			break;
		case SSIDE: sosedb[bound_id].emissivity = b[ibfound].radiation.emissN;
			break;
		case NSIDE: sosedb[bound_id].emissivity = b[ibfound].radiation.emissS;
			break;
		case BSIDE: sosedb[bound_id].emissivity = b[ibfound].radiation.emissT;
			break;
		case TSIDE: sosedb[bound_id].emissivity = b[ibfound].radiation.emissB;
			break;
		}
	}
	else {
		printf("error: emissivity calculation block unfound. in constr_sosedb_temp.\n");
		system("PAUSE");
		exit(1);
	}
	integer jpos;

	bool bfind = false;
	sosedb[bound_id].MCB = ls + lw; // Инициализация.
	for (integer j = 0; j<ls; j++) {
		if (s[j].iPlane == iplane) {
			// Важно не только попадание в фокус объекта но и нахождение с объектом на одном уровне:
			switch (iplane) {
			case XY: if ((p_centerG.x>s[j].g.xS) && (p_centerG.x<s[j].g.xE) && (p_centerG.y>s[j].g.yS) && (p_centerG.y<s[j].g.yE) && (fabs(p_centerG.z - s[j].g.zE)<admission)) { bfind = true; jpos = j; } break;
			case YZ: if ((p_centerG.z>s[j].g.zS) && (p_centerG.z<s[j].g.zE) && (p_centerG.y>s[j].g.yS) && (p_centerG.y<s[j].g.yE) && (fabs(p_centerG.x - s[j].g.xE)<admission)) { bfind = true; jpos = j; } break;
			case XZ: if ((p_centerG.x>s[j].g.xS) && (p_centerG.x<s[j].g.xE) && (p_centerG.z>s[j].g.zS) && (p_centerG.z<s[j].g.zE) && (fabs(p_centerG.y - s[j].g.yE)<admission)) { bfind = true; jpos = j; } break;
			}
		}
	}
	if (bfind) {
		// внешний источник тепла
		//printf("source out found...\n"); // debug
		//getchar();
		sosedb[bound_id].MCB = jpos;
	}
	else {
		for (integer j = 0; j<lw; j++) {
			if (w[j].iPlane == iplane) {
				// Важно не только попадание в фокус объекта но и нахождение с объектом на одном уровне:
				switch (iplane) {
				case XY: if ((p_centerG.x>w[j].g.xS) && (p_centerG.x<w[j].g.xE) && (p_centerG.y>w[j].g.yS) && (p_centerG.y<w[j].g.yE) && (fabs(p_centerG.z - w[j].g.zE)<admission)) { bfind = true; jpos = j; } break;
				case YZ: if ((p_centerG.z>w[j].g.zS) && (p_centerG.z<w[j].g.zE) && (p_centerG.y>w[j].g.yS) && (p_centerG.y<w[j].g.yE) && (fabs(p_centerG.x - w[j].g.xE)<admission)) { bfind = true; jpos = j; } break;
				case XZ: if ((p_centerG.x>w[j].g.xS) && (p_centerG.x<w[j].g.xE) && (p_centerG.z>w[j].g.zS) && (p_centerG.z<w[j].g.zE) && (fabs(p_centerG.y - w[j].g.yE)<admission)) { bfind = true; jpos = j; } break;
				}
			}
		}
		if (bfind) {
			// найдена стека
			//printf("ambient wall found...\n"); // debug
			//getchar();
			sosedb[bound_id].MCB = ls + jpos;
		}
		else {
			// граничное условие по умолчанию
			sosedb[bound_id].MCB = ls + lw;
		}

	}


} // obrabotka_granichnoi_grani

// 9.03.2019
// CALC_GG(G,GG); // Вычисляет грань GG по грани G.
void CALC_GG(integer G, integer &GG)
{
	// после maxelm_memo внутренних КО.
	// Вычисление дальнего соседа
	switch (G) {
	case ESIDE: GG = EE; break; // ESIDE
	case NSIDE: GG = NN; break; // NSIDE
	case TSIDE: GG = TTSIDE; break; // TSIDE
	case WSIDE: GG = WW; break; // WSIDE
	case SSIDE: GG = SS; break; // SSIDE
	case BSIDE: GG = BB; break; // BSIDE
	}
} // CALC_GG

// 9.03.2019
// INIT_SOSEDI(sosedi, maxelm, ESIDE, GG, false, maxelm_memo + maxbound,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE, false,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE); // заполнение sosedi[G][maxelm].
void INIT_SOSEDI(ALICE_PARTITION** &sosedi, integer maxelm, integer G,integer &GG, bool bdroblenie4G,
	integer iNODE1G, integer iNODE2G, integer iNODE3G, integer iNODE4G,
	bool bdroblenie4GG, integer iNODE1GG, integer iNODE2GG, integer iNODE3GG, integer iNODE4GG)
{
	// граничная грань:
	sosedi[G][maxelm].bdroblenie4 = bdroblenie4G;
	sosedi[G][maxelm].iNODE1 = iNODE1G; // граничные КО нумеруются в последнюю очередь,	
	sosedi[G][maxelm].iNODE2 = iNODE2G;
	sosedi[G][maxelm].iNODE3 = iNODE3G;
	sosedi[G][maxelm].iNODE4 = iNODE4G;
	// после maxelm_memo внутренних КО.
	// Вычисление дальнего соседа
	CALC_GG(G, GG); // Вычисляет грань GG по грани G.
	sosedi[GG][maxelm].bdroblenie4 = bdroblenie4GG;
	sosedi[GG][maxelm].iNODE1 = iNODE1GG; // соседа нет.
	sosedi[GG][maxelm].iNODE2 = iNODE2GG;
	sosedi[GG][maxelm].iNODE3 = iNODE3GG;
	sosedi[GG][maxelm].iNODE4 = iNODE4GG;
} // CALC_SOSEDI

// 9.03.2019
// DEFINE_BOUNDARY_PROPERTIES_TEMP(prop, prop_b, maxelm, maxbound); // В граничный узел сносятся свойства прилегающего КО.
void DEFINE_BOUNDARY_PROPERTIES_TEMP(doublereal **prop, doublereal** &prop_b, integer maxelm, integer maxbound) 
{
	// В граничный узел сносятся свойства прилегающего КО:
	prop_b[RHO][maxbound] = prop[RHO][maxelm];
	prop_b[HEAT_CAPACITY][maxbound] = prop[HEAT_CAPACITY][maxelm];
	prop_b[LAM][maxbound] = prop[LAM][maxelm];
	prop_b[MULT_LAM_X][maxbound] = prop[MULT_LAM_X][maxelm];
	prop_b[MULT_LAM_Y][maxbound] = prop[MULT_LAM_Y][maxelm];
	prop_b[MULT_LAM_Z][maxbound] = prop[MULT_LAM_Z][maxelm];
} // DEFINE_BOUNDARY_PROPERTIES_TEMP

// 9.03.2019
// DEFINE_BOUNDARY_PROPERTIES_FLOW(prop, prop_b, maxelm, maxbound); // В граничный узел сносятся свойства прилегающего КО.
void DEFINE_BOUNDARY_PROPERTIES_FLOW(doublereal **prop, doublereal** &prop_b, integer maxelm, integer maxbound) 
{
	// В граничный узел сносятся свойства прилегающего КО:
	prop_b[RHO][maxbound] = prop[RHO][maxelm];
	prop_b[MU][maxbound] = prop[MU][maxelm];
	prop_b[BETA_T][maxbound] = prop[BETA_T][maxelm];
}// DEFINE_BOUNDARY_PROPERTIES_FLOW

// 9.03.2019
// CALC_SQUARE(G, maxelm, nvtx, pa); // Вычисление площади грани ячейки.
doublereal CALC_SQUARE(integer G, integer maxelm, integer** nvtx, TOCHKA* &pa) 
{
	// вычисление размеров текущего контрольного объёма:
	doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контроольного объёма
	volume3D(maxelm, nvtx, pa, dx_loc, dy_loc, dz_loc);
	doublereal dSloc = 0.0;

	switch (G) {
	case ESIDE: dSloc = dy_loc * dz_loc; break; // ESIDE
	case NSIDE: dSloc = dx_loc * dz_loc; break; // NSIDE
	case TSIDE: dSloc = dx_loc * dy_loc; break; // TSIDE
	case WSIDE: dSloc = dy_loc * dz_loc; break; // WSIDE
	case SSIDE: dSloc = dx_loc * dz_loc; break; // SSIDE
	case BSIDE: dSloc = dx_loc * dy_loc; break; // BSIDE
	}

	dSloc = fabs(dSloc); // Площадь всегда >=0..

	return (dSloc);
} // CALC_SQUARE

  // 13.04.2019
  // CALC_CENTERG(G, maxelm, nvtx, pa); // Вычисление координат центра грани ячейки.
TOCHKA CALC_CENTERG(integer G, integer iP, integer** nvtx, TOCHKA* &pa)
{
	// вычисление размеров текущего контрольного объёма:
	doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контроольного объёма
	volume3D(iP, nvtx, pa, dx_loc, dy_loc, dz_loc);
	TOCHKA p;

	p.x = 0.5*(pa[nvtx[1][iP] - 1].x + pa[nvtx[0][iP] - 1].x);
	p.y = 0.5*(pa[nvtx[2][iP] - 1].y + pa[nvtx[0][iP] - 1].y);
	p.z = 0.5*(pa[nvtx[4][iP] - 1].z + pa[nvtx[0][iP] - 1].z);

	switch (G) {
		case ESIDE: p.x += 0.5*dx_loc; break; // ESIDE
		case NSIDE: p.y += 0.5*dy_loc; break; // NSIDE
		case TSIDE: p.z += 0.5*dz_loc; break; // TSIDE
		case WSIDE: p.x -= 0.5*dx_loc; break; // WSIDE
		case SSIDE: p.y -= 0.5*dy_loc; break; // SSIDE
		case BSIDE: p.z -= 0.5*dz_loc; break; // BSIDE
	}
	
	return (p);
} // CALC_CENTERG

// Вычисление соседей для каждого внутреннего узла. 
// Причём множество соседей выбирается как среди 
// внутренних КО так и среди граничных КО.
// Создана на базе calculate_max_bound_temp.
// Самое начало реализации. 10.сентября.2016 14_20. 
// Середина разработки 21 сентября 2016.
// Окончил писать и жду запуска для тестирования 22 сетября 2016.
// Написано очень плохо, много дублирующегося кода, за большим
// объёмом кода не прослеживается логика исполнения данной функции.
// Начало рефакторинга 9.марта.2019. (сократил в этом модуле в двух этих 
// функциях 6945 строк кода с сохранением работоспособности.
//  Был объём модуля 29601 строк стало 22656 строк.
// Метод работает верно. Проверено 10.03.2019.
void constr_sosedi_prop_b_alice(octTree* &oc, BOUND* &sosedb, bool* &binternalsource, 
	ALICE_PARTITION** &sosedi, doublereal **prop, doublereal** &prop_b,
	integer maxbound_out, integer maxelm_memo, integer maxelm_memo1,
	BLOCK* b, integer lb, SOURCE* s, integer ls, WALL* w, integer lw, 
	integer* &whot_is_block,
	integer** nvtx, TOCHKA* &pa)
{

	integer G, GG; // грань ячейки дискретизации.

	// Выделение оперативной памяти.
	sosedi = NULL;
	sosedi = new ALICE_PARTITION*[12];
	if (sosedi == NULL) {
		// недостаточно памяти на данном оборудовании.
		printf("Problem : not enough memory on your equipment for sosedi constr struct_alice...\n");
		printf("Please any key to exit...\n");
		//getchar();
		system("pause");
		exit(1);
	}
	for (integer i = 0; i<12; i++) sosedi[i] = NULL;
	for (integer i = 0; i<12; i++) {
		sosedi[i] = new ALICE_PARTITION[maxelm_memo+2];
		if (sosedi[i] == NULL) {
			// недостаточно памяти на данном оборудовании.
#if doubleintprecision == 1
			printf("Problem : not enough memory on your equipment for sosedi[%lld] constr struct_alice...\n", i);
#else
			printf("Problem : not enough memory on your equipment for sosedi[%d] constr struct_alice...\n", i);
#endif
			
			printf("Please any key to exit...\n");
			//getchar();
			system("pause");
			exit(1);
		}
	}
	// Инициализация. 09.03.2019
	for (integer i = 0; i < maxelm_memo + 2; i++) {
		integer GGloc;
		INIT_SOSEDI(sosedi, i, ESIDE, GGloc, false,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE, false,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE); // заполнение sosedi[G][maxelm].
		INIT_SOSEDI(sosedi, i, WSIDE, GGloc, false,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE, false,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE); // заполнение sosedi[G][maxelm].
		INIT_SOSEDI(sosedi, i, NSIDE, GGloc, false,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE, false,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE); // заполнение sosedi[G][maxelm].
		INIT_SOSEDI(sosedi, i, WSIDE, GGloc, false,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE, false,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE); // заполнение sosedi[G][maxelm].
		INIT_SOSEDI(sosedi, i, TSIDE, GGloc, false,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE, false,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE); // заполнение sosedi[G][maxelm].
		INIT_SOSEDI(sosedi, i, BSIDE, GGloc, false,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE, false,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE); // заполнение sosedi[G][maxelm].
	}

	// Выделение оперативной памяти:
	prop_b = NULL;
	prop_b = new doublereal*[6];
	if (prop_b == NULL) {
		// недостаточно памяти на данном оборудовании.
		printf("Problem : not enough memory on your equipment for prop_b constr struct_alice...\n");
		printf("Please any key to exit...\n");
		//getchar();
		system("pause");
		exit(1);
	}
	
	for (integer i = 0; i<6; i++) prop_b[i] = NULL;
	for (integer i = 0; i<6; i++) {
		prop_b[i] = new doublereal[maxbound_out];
		if (prop_b[i] == NULL) {
			// недостаточно памяти на данном оборудовании.
#if doubleintprecision == 1
			printf("Problem : not enough memory on your equipment for prop_b[%lld] constr struct_alice...\n", i);
#else
			printf("Problem : not enough memory on your equipment for prop_b[%d] constr struct_alice...\n", i);
#endif
			
			printf("Please any key to exit...\n");
			//getchar();
			system("pause");
			exit(1);
		}
	}
	// инициализация.
	for (integer i = 0; i < maxbound_out; i++) {
		// Дюраль Д16Т.
		prop_b[RHO][i] = 2800.0;
		prop_b[HEAT_CAPACITY][i] = 921.0;
		prop_b[LAM][i] = 164.0;
		prop_b[MULT_LAM_X][i] = 1.0;
		prop_b[MULT_LAM_Y][i] = 1.0;
		prop_b[MULT_LAM_Z][i] = 1.0;
	}


	// Алгоритм:
	// Выделение оперативной памяти:
	// gran[0..5][0..maxbound-1]
	sosedb = NULL;
	sosedb = new BOUND[maxbound_out];
	if (sosedb == NULL) {
		// недостаточно памяти на данном оборудовании.
		printf("Problem : not enough memory on your equipment for sosedb constr struct_alice...\n");
		printf("Please any key to exit...\n");
		//getchar();
		system("pause");
		exit(1);
	}
	binternalsource = NULL;
	binternalsource = new bool[maxbound_out];
	// оператор new не требует проверки на null.
	//if (binternalsource == NULL) {
		// недостаточно памяти на данном оборудовании.
		//printf("Problem : not enough memory on your equipment for binternalsource constr struct_alice...\n");
		//printf("Please any key to exit...\n");
		//getchar();
		//system("pause");
		//exit(1);
	//}
	for (integer i = 0; i<maxbound_out; i++) binternalsource[i] = false; // инициализация.

	integer maxbound = 0; // инициализация.
	integer maxelm = 0;
	bool *bvisit = NULL;
	bvisit = new bool[maxelm_memo+2];
	// оператор new не требует проверки на null.
	//if (bvisit == NULL) {
		// недостаточно памяти на данном оборудовании.
		//printf("Problem : not enough memory on your equipment for bvisit constr struct_alice...\n");
		//printf("Please any key to exit...\n");
		//getchar();
		//system("pause");
		//exit(1);
	//}
	for (integer j = 0; j<maxelm_memo+2; j++) bvisit[j] = false; // признак посещения узла.
	doublereal x_c, y_c, z_c;
	integer iplane;

	top_ALICE_STACK = 0;
	if (oc->link0 != NULL) {
		my_ALICE_STACK[top_ALICE_STACK].link = (oc->link0);
		my_ALICE_STACK[top_ALICE_STACK].maxx = oc->link0->maxx;
		my_ALICE_STACK[top_ALICE_STACK].minx = oc->link0->minx;
		my_ALICE_STACK[top_ALICE_STACK].maxy = oc->link0->maxy;
		my_ALICE_STACK[top_ALICE_STACK].miny = oc->link0->miny;
		my_ALICE_STACK[top_ALICE_STACK].maxz = oc->link0->maxz;
		my_ALICE_STACK[top_ALICE_STACK].minz = oc->link0->minz;
		top_ALICE_STACK++;
	}
	if (oc->link1 != NULL) {
		my_ALICE_STACK[top_ALICE_STACK].link = (oc->link1);
		my_ALICE_STACK[top_ALICE_STACK].maxx = oc->link1->maxx;
		my_ALICE_STACK[top_ALICE_STACK].minx = oc->link1->minx;
		my_ALICE_STACK[top_ALICE_STACK].maxy = oc->link1->maxy;
		my_ALICE_STACK[top_ALICE_STACK].miny = oc->link1->miny;
		my_ALICE_STACK[top_ALICE_STACK].maxz = oc->link1->maxz;
		my_ALICE_STACK[top_ALICE_STACK].minz = oc->link1->minz;
		top_ALICE_STACK++;
	}
	if (oc->link2 != NULL) {
		my_ALICE_STACK[top_ALICE_STACK].link = (oc->link2);
		my_ALICE_STACK[top_ALICE_STACK].maxx = oc->link2->maxx;
		my_ALICE_STACK[top_ALICE_STACK].minx = oc->link2->minx;
		my_ALICE_STACK[top_ALICE_STACK].maxy = oc->link2->maxy;
		my_ALICE_STACK[top_ALICE_STACK].miny = oc->link2->miny;
		my_ALICE_STACK[top_ALICE_STACK].maxz = oc->link2->maxz;
		my_ALICE_STACK[top_ALICE_STACK].minz = oc->link2->minz;
		top_ALICE_STACK++;
	}
	if (oc->link3 != NULL) {
		my_ALICE_STACK[top_ALICE_STACK].link = (oc->link3);
		my_ALICE_STACK[top_ALICE_STACK].maxx = oc->link3->maxx;
		my_ALICE_STACK[top_ALICE_STACK].minx = oc->link3->minx;
		my_ALICE_STACK[top_ALICE_STACK].maxy = oc->link3->maxy;
		my_ALICE_STACK[top_ALICE_STACK].miny = oc->link3->miny;
		my_ALICE_STACK[top_ALICE_STACK].maxz = oc->link3->maxz;
		my_ALICE_STACK[top_ALICE_STACK].minz = oc->link3->minz;
		top_ALICE_STACK++;
	}
	if (oc->link4 != NULL) {
		my_ALICE_STACK[top_ALICE_STACK].link = (oc->link4);
		my_ALICE_STACK[top_ALICE_STACK].maxx = oc->link4->maxx;
		my_ALICE_STACK[top_ALICE_STACK].minx = oc->link4->minx;
		my_ALICE_STACK[top_ALICE_STACK].maxy = oc->link4->maxy;
		my_ALICE_STACK[top_ALICE_STACK].miny = oc->link4->miny;
		my_ALICE_STACK[top_ALICE_STACK].maxz = oc->link4->maxz;
		my_ALICE_STACK[top_ALICE_STACK].minz = oc->link4->minz;
		top_ALICE_STACK++;
	}
	if (oc->link5 != NULL) {
		my_ALICE_STACK[top_ALICE_STACK].link = (oc->link5);
		my_ALICE_STACK[top_ALICE_STACK].maxx = oc->link5->maxx;
		my_ALICE_STACK[top_ALICE_STACK].minx = oc->link5->minx;
		my_ALICE_STACK[top_ALICE_STACK].maxy = oc->link5->maxy;
		my_ALICE_STACK[top_ALICE_STACK].miny = oc->link5->miny;
		my_ALICE_STACK[top_ALICE_STACK].maxz = oc->link5->maxz;
		my_ALICE_STACK[top_ALICE_STACK].minz = oc->link5->minz;
		top_ALICE_STACK++;
	}
	if (oc->link6 != NULL) {
		my_ALICE_STACK[top_ALICE_STACK].link = (oc->link6);
		my_ALICE_STACK[top_ALICE_STACK].maxx = oc->link6->maxx;
		my_ALICE_STACK[top_ALICE_STACK].minx = oc->link6->minx;
		my_ALICE_STACK[top_ALICE_STACK].maxy = oc->link6->maxy;
		my_ALICE_STACK[top_ALICE_STACK].miny = oc->link6->miny;
		my_ALICE_STACK[top_ALICE_STACK].maxz = oc->link6->maxz;
		my_ALICE_STACK[top_ALICE_STACK].minz = oc->link6->minz;
		top_ALICE_STACK++;
	}
	if (oc->link7 != NULL) {
		my_ALICE_STACK[top_ALICE_STACK].link = (oc->link7);
		my_ALICE_STACK[top_ALICE_STACK].maxx = oc->link7->maxx;
		my_ALICE_STACK[top_ALICE_STACK].minx = oc->link7->minx;
		my_ALICE_STACK[top_ALICE_STACK].maxy = oc->link7->maxy;
		my_ALICE_STACK[top_ALICE_STACK].miny = oc->link7->miny;
		my_ALICE_STACK[top_ALICE_STACK].maxz = oc->link7->maxz;
		my_ALICE_STACK[top_ALICE_STACK].minz = oc->link7->minz;
		top_ALICE_STACK++;
	}
	while (top_ALICE_STACK > 0) {
		if (my_ALICE_STACK[top_ALICE_STACK - 1].link != NULL) {
			if (my_ALICE_STACK[top_ALICE_STACK - 1].link->dlist == true) {
				
				octTree* octree1 = my_ALICE_STACK[top_ALICE_STACK - 1].link;
				TOCHKA p;
				p.x = 0.125*(octree1->p0.x + octree1->p1.x + octree1->p2.x + octree1->p3.x + octree1->p4.x + octree1->p5.x + octree1->p6.x + octree1->p7.x);
				p.y = 0.125*(octree1->p0.y + octree1->p1.y + octree1->p2.y + octree1->p3.y + octree1->p4.y + octree1->p5.y + octree1->p6.y + octree1->p7.y);
				p.z = 0.125*(octree1->p0.z + octree1->p1.z + octree1->p2.z + octree1->p3.z + octree1->p4.z + octree1->p5.z + octree1->p6.z + octree1->p7.z);
				integer ib;
				bool inDomain = false;
				inDomain = in_model_temp(p, ib, b, lb); // TEMPERATURE
				//integer ib1;
				//bool bi_fluid = in_model_flow(p, ib1, b, lb);
				bool bi_fluid = (b[ib].itype == FLUID);

				if ((p.x < b[0].g.xS) || (p.x > b[0].g.xE) || (p.y < b[0].g.yS) || (p.y > b[0].g.yE) || (p.z < b[0].g.zS) || (p.z > b[0].g.zE)) {
					// Лежит за пределами кабинета.
					inDomain = false;
				}
				if (inDomain) {
					// Узел принадлекжит расчётной области и его номер maxelm==octree1->inum_TD.

					// сосед не существует
					// записываем идентификатор грани и
					// увеличиваем счётчик граней.
					if (!octree1->b4E) {
						if ((octree1->linkE == NULL)||((octree1->linkE != NULL)&&(octree1->linkE->inum_TD == 0))) {
							// Если мы находимся на границе расчётной области или если мы граничим с HOLLOW блоком.
							// И у нас нет дробления на 4 а только один сосед.

							G = ESIDE;
							// граничная грань:
							INIT_SOSEDI(sosedi, maxelm, ESIDE, GG, false, maxelm_memo + maxbound,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE, false,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE, NON_EXISTENT_NODE); // заполнение sosedi[G][maxelm].
							DEFINE_BOUNDARY_PROPERTIES_TEMP(prop, prop_b, maxelm, maxbound); // В граничный узел сносятся свойства прилегающего КО.
							
							doublereal dSloc = CALC_SQUARE(G, maxelm, nvtx, pa); // Вычисление площади грани ячейки.
							TOCHKA p_centerG = CALC_CENTERG(G, maxelm, nvtx, pa); // Вычисление геометрического центра грани.

							// G - грань по которой ставится граничное условие.
							// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
							// bound_id - номер граничного узла, начиная с нуля.
							// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
							if (octree1->b4W) {
								obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
									octree1->linkW0->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, false, false, p_centerG);
							}
							else {
								if (octree1->linkW == NULL) {
									obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
										 NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, false, false, p_centerG);
								}
								else {
									obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
										octree1->linkW->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, false, false, p_centerG);
								}
							}
							maxbound++;
						}
						else if (octree1->linkE != NULL) {// сосед существует.
							
								
								if (bvisit[octree1->linkE->inum_TD - 1]) {// Сосед был посещен ранее
									
									// return_current_node_number - Заглянули в соседа и вернулись из соседа в текущий узел.
									// Номер должен быть тем же самым и поэтому этот код бессмысленен.
									// Строго внутренние плоские бесконечно тонкие источники тепла внутри расчётной области больше не поддерживаются.
									// Поэтому от этого кода можно безболезненно избавиться. Поддерживаются только объёмные источники тепла.
									// 10.03.2019
									integer return_current_node_number = sosedi[WSIDE][octree1->linkE->inum_TD - 1].iNODE1;
									if (octree1->linkE->b4W) {
										integer current_node_number = octree1->inum_TD - 1;// Текущий внутренний номер.
										bool bcontinue = true;
										if (bcontinue) {
											if (octree1->linkE->linkW0 != NULL) {
												if (octree1->linkE->linkW0->inum_TD - 1 == current_node_number) {
													bcontinue = false; // iNODE1
													if (sosedi[WSIDE][octree1->linkE->inum_TD - 1].iNODE1 == current_node_number) {
														return_current_node_number = sosedi[WSIDE][octree1->linkE->inum_TD - 1].iNODE1;
													}
													else {
														printf("this can not be E W0\n");
														//getchar();
														system("PAUSE");
													}
												}
											}
										}
										if (bcontinue) {
											if (octree1->linkE->linkW3!=NULL) {
												if (octree1->linkE->linkW3->inum_TD - 1 == current_node_number) {
													bcontinue = false; // iNODE2
													if (sosedi[WSIDE][octree1->linkE->inum_TD - 1].iNODE2 == current_node_number) {
														return_current_node_number = sosedi[WSIDE][octree1->linkE->inum_TD - 1].iNODE2;
													}
													else {
														printf("this can not be E W3\n");
														//getchar();
														system("PAUSE");
													}
												}
											}
										}
										if (bcontinue) {
											if (octree1->linkE->linkW4!=NULL) {
												if (octree1->linkE->linkW4->inum_TD - 1 == current_node_number) {
													bcontinue = false; // iNODE3
													if (sosedi[WSIDE][octree1->linkE->inum_TD - 1].iNODE3 == current_node_number) {
														return_current_node_number = sosedi[WSIDE][octree1->linkE->inum_TD - 1].iNODE3;
													}
													else {
														printf("this can not be E W4\n");
														//getchar();
														system("PAUSE");
													}
												}
											}
										}
										if (bcontinue) {
											if (octree1->linkE->linkW7!=NULL) {
												if (octree1->linkE->linkW7->inum_TD - 1 == current_node_number) {
													bcontinue = false; // iNODE4
													if (sosedi[WSIDE][octree1->linkE->inum_TD - 1].iNODE4 == current_node_number) {
														return_current_node_number = sosedi[WSIDE][octree1->linkE->inum_TD - 1].iNODE4;
													}
													else {
														printf("this can not be E W7\n");
														//getchar();
														system("PAUSE");
													}
												}
											}
										}

									}


									// узел уже был посещён
									// Здесь неким образом модифицировалось gran_t
									if (return_current_node_number >= maxelm_memo) {

										printf("WARNING!!! internal source gran ESIDE 10.03.2019\n");
										system("PAUSE");

										// Это граничная грань внутреннего источника тепла.
										G = ESIDE;
										// граничная грань:
										INIT_SOSEDI(sosedi, maxelm, ESIDE, GG, false, sosedi[WSIDE][octree1->linkE->inum_TD - 1].iNODE1,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE, false,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE); // заполнение sosedi[G][maxelm].

										if (bi_fluid) {
											// FLUID
											// В граничный узел сносятся свойства прилегающего КО.
											DEFINE_BOUNDARY_PROPERTIES_TEMP(prop, prop_b, octree1->linkE->inum_TD - 1, sosedi[G][maxelm].iNODE1 - maxelm_memo); 
										}
										else {
											// SOLID
											// В граничный узел сносятся свойства прилегающего КО.
											DEFINE_BOUNDARY_PROPERTIES_TEMP(prop, prop_b, maxelm, sosedi[G][maxelm].iNODE1 - maxelm_memo); 
										}

										binternalsource[sosedi[WSIDE][octree1->linkE->inum_TD - 1].iNODE1 - maxelm_memo] = true; // внутренний источник тепла на границе жидкости и твёрдого тела

										doublereal dSloc = CALC_SQUARE(G, maxelm, nvtx, pa); // Вычисление площади грани ячейки.
										TOCHKA p_centerG = CALC_CENTERG(G, maxelm, nvtx, pa); // Вычисление геометрического центра грани.

										// G - грань по которой ставится граничное условие.
										// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
										// bound_id - номер граничного узла, начиная с нуля.
										// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
										if (octree1->b4W) {
											obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, sosedi[WSIDE][octree1->linkE->inum_TD - 1].iNODE1 - maxelm_memo,
												octree1->linkW0->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, true, true, p_centerG);
										}
										else {
											if (octree1->linkW == NULL) {
												obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, sosedi[WSIDE][octree1->linkE->inum_TD - 1].iNODE1 - maxelm_memo,
													 NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, true, true, p_centerG);
											}
											else {
												obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, sosedi[WSIDE][octree1->linkE->inum_TD - 1].iNODE1 - maxelm_memo,
													octree1->linkW->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, true, true, p_centerG);
											}
										}

									}
									else {
										// Это строго внутренний узел.
										// Внутренний узел.
										G = ESIDE;
										// внутренняя грань (нумерация начинается с нуля):
										// Ссылка на строго внутреннего соседа который точно существует.
										sosedi[G][maxelm].iNODE1 = octree1->linkE->inum_TD - 1; // граничные КО нумеруются в последнюю очередь,
										sosedi[G][maxelm].bdroblenie4 = false;
										sosedi[G][maxelm].iNODE2 = NON_EXISTENT_NODE;
										sosedi[G][maxelm].iNODE3 = NON_EXISTENT_NODE;
										sosedi[G][maxelm].iNODE4 = NON_EXISTENT_NODE;
										// после maxelm_memo внутренних КО.
										//CALC_GG(G, GG); // Вычисляет грань GG по грани G.
										// Определение дальнего соседа :
										// TODO 13 сентября 2016.
										//if ((octree1->linkE->linkE == NULL) || ((octree1->linkE->linkE != NULL)&&(octree1->linkE->linkE->inum_TD==0))) {
											//sosedi[GG][maxelm].iNODE1 = NON_EXISTENT_NODE;
										//}
										//else {
											//sosedi[GG][maxelm].iNODE1 = уникальный номер граничного узла;
										// уникальный номер граничного узла неизвестен, нужен второй проход для дальнего соседа.
										//}

									}

								}
								else { // Узел еще не был посещен ранее.
								// Это строго внутренний узел.

									bool binc;
									integer lsid;
									// узнать координаты центра грани и ориентацию в пространстве
									x_c = octree1->p1.x;
									y_c = 0.5*(octree1->p1.y + octree1->p2.y);
									z_c = 0.5*(octree1->p1.z + octree1->p5.z);
									iplane = YZ;
									patch_maxbound(iplane, s, ls, x_c, y_c, z_c, maxbound, binc, lsid);
									if (binc) {
										// Внутренний источник тепла.
										// Внутренние плоские бесконечно тонкие источники тепла скорее всего больше не поддерживаются. 10,03,2019
										G = ESIDE;
										// граничная грань:
										INIT_SOSEDI(sosedi, maxelm, ESIDE, GG, false, maxelm_memo + maxbound - 1,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE, false,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE); // заполнение sosedi[G][maxelm].

										if (bi_fluid) {
											// FLUID
											DEFINE_BOUNDARY_PROPERTIES_TEMP(prop, prop_b, octree1->linkE->inum_TD - 1, maxbound-1); // В граничный узел сносятся свойства прилегающего КО.
										}
										else {
											// SOLID
											DEFINE_BOUNDARY_PROPERTIES_TEMP(prop, prop_b, maxelm, maxbound-1); // В граничный узел сносятся свойства прилегающего КО.
										}

										binternalsource[maxbound - 1] = true; // внутренний источник тепла на границе жидкости и твёрдого тела

										doublereal dSloc = CALC_SQUARE(G, maxelm, nvtx, pa); // Вычисление площади грани ячейки.
										TOCHKA p_centerG = CALC_CENTERG(G, maxelm, nvtx, pa); // Вычисление геометрического центра грани.

										// G - грань по которой ставится граничное условие.
										// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
										// bound_id - номер граничного узла, начиная с нуля.
										// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
										if (octree1->b4W) {
											obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound - 1,
												octree1->linkW0->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, true, false, p_centerG);
										}
										else {
											if (octree1->linkW == NULL) {
												obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound - 1,
													 NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, true, false, p_centerG);
											}
											else {
												obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound - 1,
													octree1->linkW->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, true, false, p_centerG);
											}
										}

									}
									else {
										// Внутренний узел который не был посещен ранее.

										G = ESIDE;
										// внутренняя грань (нумерация начинается с нуля):
										sosedi[G][maxelm].iNODE1 = octree1->linkE->inum_TD - 1; // граничные КО нумеруются в последнюю очередь,
										sosedi[G][maxelm].bdroblenie4 = false;
										sosedi[G][maxelm].iNODE2 = NON_EXISTENT_NODE;
										sosedi[G][maxelm].iNODE3 = NON_EXISTENT_NODE;
										sosedi[G][maxelm].iNODE4 = NON_EXISTENT_NODE;
										// после maxelm_memo внутренних КО.
										// Вычисление дальнего соседа
										// CALC_GG(G, GG); // Вычисляет грань GG по грани G.
										// Определение дальнего соседа :
										// TODO 13 сентября 2016.

									}

								}
							
						}
					}
					else {
					// Присутствуют все 4 соседа: E: 1,2,5,6.

						if (((octree1->linkE1 == NULL) && (octree1->linkE2 == NULL) && (octree1->linkE5 == NULL) && (octree1->linkE6 == NULL))||
						(((octree1->linkE1 != NULL) && (octree1->linkE2 != NULL) && (octree1->linkE5 != NULL) && (octree1->linkE6 != NULL) &&
							(octree1->linkE1->inum_TD == 0) && (octree1->linkE2->inum_TD == 0) && (octree1->linkE5->inum_TD == 0) && (octree1->linkE6->inum_TD == 0)))||
							(((octree1->linkE1 != NULL) && (octree1->linkE2 == NULL) && (octree1->linkE5 != NULL) && (octree1->linkE6 == NULL) &&
							(octree1->linkE1->inum_TD == 0) && (octree1->linkE5->inum_TD == 0)))||
							(((octree1->linkE1 != NULL) && (octree1->linkE2 != NULL) && (octree1->linkE5 == NULL) && (octree1->linkE6 == NULL) &&
							(octree1->linkE1->inum_TD == 0) && (octree1->linkE2->inum_TD == 0)))||
							(((octree1->linkE1 != NULL) && (octree1->linkE2 == NULL) && (octree1->linkE5 == NULL) && (octree1->linkE6 == NULL) &&
							(octree1->linkE1->inum_TD == 0))))
						{
							// 1 - Все четыре соседа на границе расчётной области.
							// 2 - Все четыре соседа лежат на границе HOLLOW блока.
							// 3 - E2 && E6 граница РО, E1 && E5 граница HOLLOW блока.
							// 4 - E1 && E2 граница HOLLOW блока, E5 && E6 граница РО.
							// 5 - E1 граница HOLLOW блока, E2 && E5 && E6 граница РО. 

							G = ESIDE;
							// граничная грань:
							INIT_SOSEDI(sosedi, maxelm, ESIDE, GG, false, maxelm_memo + maxbound,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE, false,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE); // заполнение sosedi[G][maxelm].							
							DEFINE_BOUNDARY_PROPERTIES_TEMP(prop, prop_b, maxelm, maxbound); // В граничный узел сносятся свойства прилегающего КО.
							
							doublereal dSloc = CALC_SQUARE(G, maxelm, nvtx, pa); // Вычисление площади грани ячейки.
							TOCHKA p_centerG = CALC_CENTERG(G, maxelm, nvtx, pa); // Вычисление геометрического центра грани.
							
							// G - грань по которой ставится граничное условие.
							// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
							// bound_id - номер граничного узла, начиная с нуля.
							// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
							if (octree1->b4W) {
								obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
									octree1->linkW0->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, false, false, p_centerG);
							}
							else {
								if (octree1->linkW == NULL) {
									obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
										 NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, false, false, p_centerG);
								}
								else {

									obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
										octree1->linkW->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, false, false, p_centerG);
								}
							}
							// Зачем дробить если граница заканчивается NULL.
							maxbound++; // Вместо четырёх один.
						}						
						else {

							bool binc1, binc2, binc3, binc4;
							integer lsid1, lsid2, lsid3, lsid4;

							if (octree1->linkE1 == NULL) {
#if doubleintprecision == 1
								printf("error NULL sosedi E1 elm=%lld\n", maxelm);
#else
								printf("error NULL sosedi E1 elm=%d\n", maxelm);
#endif
								
								//getchar();
								system("PAUSE");
								/*
								G = ESIDE;
								// граничная грань:
								sosedi[G][maxelm].iNODE1 = maxelm_memo + maxbound; // граничные КО нумеруются в последнюю очередь,
								sosedi[G][maxelm].bdroblenie4 = true;

								// после maxelm_memo внутренних КО.
								// Вычисление дальнего соседа
								CALC_GG(G,GG); // Вычисляет грань GG по грани G.
								sosedi[GG][maxelm].iNODE1 = NON_EXISTENT_NODE; // соседа нет.
								sosedi[GG][maxelm].bdroblenie4 = true;
								#if doubleintprecision == 1
									printf("error NULL sosedi E1 elm=%lld\n", maxelm);
								#else
									printf("error NULL sosedi E1 elm=%d\n", maxelm);
								#endif
								
								getchar();

								// Эта ситуация возможна и здесь нужно предусмотреть обработку.
								DEFINE_BOUNDARY_PROPERTIES_TEMP(prop, prop_b, maxelm, maxbound); // В граничный узел сносятся свойства прилегающего КО.
								
								doublereal dS_loc = 0.0;

								// G - грань по которой ставится граничное условие.
								// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
								// bound_id - номер граничного узла, начиная с нуля.
								// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
								if (octree1->b4W) {
									obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
										octree1->linkW0->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
								}
								else {
									if (octree1->linkW == NULL) {
										obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
											 NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
									}
									else {
										obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
											octree1->linkW->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
									}
								}

								maxbound++;
								*/
							}
							else {
								if (octree1->linkE1->inum_TD == 0) {
									G = ESIDE;
									// граничная грань:
									sosedi[G][maxelm].iNODE1 = maxelm_memo + maxbound; // граничные КО нумеруются в последнюю очередь,
									sosedi[G][maxelm].bdroblenie4 = true;

									// после maxelm_memo внутренних КО.
									// Вычисление дальнего соседа
									CALC_GG(G, GG); // Вычисляет грань GG по грани G.
									sosedi[GG][maxelm].iNODE1 = NON_EXISTENT_NODE; // соседа нет.
									sosedi[GG][maxelm].bdroblenie4 = true;

#if doubleintprecision == 1
									//printf("error inum_TD == 0 sosedi E1 elm=%lld\n", maxelm);
#else
									//printf("error inum_TD == 0 sosedi E1 elm=%d\n", maxelm);
#endif
									
									//getchar();

									// Эта ситуация возможна и здесь нужно предусмотреть обработку.
									DEFINE_BOUNDARY_PROPERTIES_TEMP(prop, prop_b, maxelm, maxbound); // В граничный узел сносятся свойства прилегающего КО.
									

									// Sy*Sz
									doublereal dS_loc = fabs(octree1->linkE1->p4.z - octree1->linkE1->p0.z)*fabs(octree1->linkE1->p3.y - octree1->linkE1->p0.y);
									TOCHKA p_centerG;
									p_centerG.z = 0.5*(octree1->linkE1->p4.z + octree1->linkE1->p0.z);
									p_centerG.y = 0.5*(octree1->linkE1->p3.y + octree1->linkE1->p0.y);
									p_centerG.x = octree1->linkE1->p0.x;


									// G - грань по которой ставится граничное условие.
									// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
									// bound_id - номер граничного узла, начиная с нуля.
									// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
									if (octree1->b4W) {
										obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
											octree1->linkW0->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
									}
									else {
										if (octree1->linkW == NULL) {
											obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
												 NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
										}
										else {
											obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
												octree1->linkW->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
										}
									}

									maxbound++;
								}
								else {
									// сосед существует.
									if (bvisit[octree1->linkE1->inum_TD - 1]) {
										// узел уже был посещён
										// Здесь неким образом модифицировалось gran_t
										// Здесь не может быть внутреннего источника тепла иначе 
										// выдавалось бы предупреждение 
#if doubleintprecision == 1
										//printf("error internal source E1 elm=%lld\n", maxelm);
#else
										//printf("error internal source E1 elm=%d\n", maxelm);
#endif
										
										//getchar();
										// В узле octree1->linkE1->inum_TD - 1 Который был посещен ранее.
										// Здесь возможен только внутренний узел.
										// Внутренний узел.
										G = ESIDE;
										// Внутренняя грань:
										sosedi[G][maxelm].iNODE1 = octree1->linkE1->inum_TD - 1; // граничные КО нумеруются в последнюю очередь,
										sosedi[G][maxelm].bdroblenie4 = true;

										// после maxelm_memo внутренних КО.
										// Вычисление дальнего соседа
										//CALC_GG(G, GG); // Вычисляет грань GG по грани G.

										// Определение дальнего соседа :
										// TODO 13 сентября 2016.

									}
									else {
										// узнать координаты центра грани и ориентацию в пространстве
										x_c = octree1->p1.x;
										y_c = 0.5*(octree1->linkE1->p0.y + octree1->linkE1->p3.y);
										z_c = 0.5*(octree1->linkE1->p0.z + octree1->linkE1->p4.z);
										iplane = YZ;
										patch_maxbound(iplane, s, ls, x_c, y_c, z_c, maxbound, binc1, lsid1);
										if (binc1) {
#if doubleintprecision == 1
											printf("error internal source E1 elm=%lld\n", maxelm);
#else
											printf("error internal source E1 elm=%d\n", maxelm);
#endif
											
											//getchar();
											system("PAUSE");
										}
										else {
											// Внутренний узел.
											G = ESIDE;
											// Внутренняя грань:
											sosedi[G][maxelm].iNODE1 = octree1->linkE1->inum_TD - 1; // граничные КО нумеруются в последнюю очередь,
											sosedi[G][maxelm].bdroblenie4 = true;

											// после maxelm_memo внутренних КО.
											// Вычисление дальнего соседа
											//CALC_GG(G, GG); // Вычисляет грань GG по грани G.

											// Определение дальнего соседа :
											// TODO 13 сентября 2016.

										}
									}
								}
							}

							if (octree1->linkE2 == NULL) {
								// 21 сентября 2016 (вырождение ячейки).

								// Заглушка.
								G = ESIDE;
								// граничная грань:
								sosedi[G][maxelm].iNODE2 = NON_EXISTENT_NODE; // граничные КО нумеруются в последнюю очередь,
								sosedi[G][maxelm].bdroblenie4 = true;

								// после maxelm_memo внутренних КО.
								// Вычисление дальнего соседа
								CALC_GG(G, GG); // Вычисляет грань GG по грани G.
								sosedi[GG][maxelm].iNODE2 = NON_EXISTENT_NODE; // соседа нет.
								sosedi[GG][maxelm].bdroblenie4 = true;


								/*
								G = ESIDE;
								// граничная грань:
								sosedi[G][maxelm].iNODE2 = maxelm_memo + maxbound; // граничные КО нумеруются в последнюю очередь,
								sosedi[G][maxelm].bdroblenie4 = true;

								// после maxelm_memo внутренних КО.
								// Вычисление дальнего соседа
								CALC_GG(G,GG); // Вычисляет грань GG по грани G.
								sosedi[GG][maxelm].iNODE2 = NON_EXISTENT_NODE; // соседа нет.
								sosedi[GG][maxelm].bdroblenie4 = true;

								#if doubleintprecision == 1
									printf("error NULL sosedi E2 elm=%lld\n", maxelm);
								#else
									printf("error NULL sosedi E2 elm=%d\n", maxelm);
								#endif
								
								getchar();

								// Эта ситуация возможна и здесь нужно предусмотреть обработку.
								DEFINE_BOUNDARY_PROPERTIES_TEMP(prop, prop_b, maxelm, maxbound); // В граничный узел сносятся свойства прилегающего КО.
								
								doublereal dS_loc = 0.0;

								// G - грань по которой ставится граничное условие.
								// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
								// bound_id - номер граничного узла, начиная с нуля.
								// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
								if (octree1->b4W) {
									obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
										octree1->linkW0->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
								}
								else {
									if (octree1->linkW == NULL) {
										obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
											 NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
									}
									else {
										obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
											octree1->linkW->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
									}
								}



								maxbound++;
								*/
							}
							else {
								if (octree1->linkE2->inum_TD == 0) {
									G = ESIDE;
									// граничная грань:
									sosedi[G][maxelm].iNODE2 = maxelm_memo + maxbound; // граничные КО нумеруются в последнюю очередь,
									sosedi[G][maxelm].bdroblenie4 = true;

									// после maxelm_memo внутренних КО.
									// Вычисление дальнего соседа
									CALC_GG(G, GG); // Вычисляет грань GG по грани G.
									sosedi[GG][maxelm].iNODE2 = NON_EXISTENT_NODE; // соседа нет.
									sosedi[GG][maxelm].bdroblenie4 = true;

#if doubleintprecision == 1
									//printf("error inum_TD == 0 sosedi E2 elm=%lld\n", maxelm);
#else
									//printf("error inum_TD == 0 sosedi E2 elm=%d\n", maxelm);
#endif
									
									//getchar();

									// Эта ситуация возможна и здесь нужно предусмотреть обработку.
									DEFINE_BOUNDARY_PROPERTIES_TEMP(prop, prop_b, maxelm, maxbound); // В граничный узел сносятся свойства прилегающего КО.

									// Sy*Sz
									doublereal dS_loc = fabs(octree1->linkE2->p4.z - octree1->linkE2->p0.z)*fabs(octree1->linkE2->p3.y - octree1->linkE2->p0.y);
									TOCHKA p_centerG;
									p_centerG.z = 0.5*(octree1->linkE2->p4.z + octree1->linkE2->p0.z);
									p_centerG.y = 0.5*(octree1->linkE2->p3.y + octree1->linkE2->p0.y);
									p_centerG.x = octree1->linkE2->p0.x;

									// G - грань по которой ставится граничное условие.
									// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
									// bound_id - номер граничного узла, начиная с нуля.
									// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
									if (octree1->b4W) {
										obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
											octree1->linkW0->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
									}
									else {
										if (octree1->linkW == NULL) {
											obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
												 NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
										}
										else {
											obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
												octree1->linkW->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
										}
									}

									maxbound++;
								}
								else {
									// сосед существует.
									if (bvisit[octree1->linkE2->inum_TD - 1]) {
										// узел уже был посещён
										// Здесь неким образом модифицировалось gran_t
										// Здесь тоже возможен лишь строго внутренний узел.
										// Внутренний узел.
										G = ESIDE;
										// Внутренняя грань:
										sosedi[G][maxelm].iNODE2 = octree1->linkE2->inum_TD - 1; // граничные КО нумеруются в последнюю очередь,
										sosedi[G][maxelm].bdroblenie4 = true;

										// после maxelm_memo внутренних КО.
										// Вычисление дальнего соседа
										CALC_GG(G, GG); // Вычисляет грань GG по грани G.

										// Определение дальнего соседа :
										// TODO 13 сентября 2016.

									}
									else {
										// узнать координаты центра грани и ориентацию в пространстве
										x_c = octree1->p1.x;
										y_c = 0.5*(octree1->linkE2->p0.y + octree1->linkE2->p3.y);
										z_c = 0.5*(octree1->linkE2->p0.z + octree1->linkE2->p4.z);
										iplane = YZ;
										patch_maxbound(iplane, s, ls, x_c, y_c, z_c, maxbound, binc2, lsid2);
										if (binc2) {
#if doubleintprecision == 1
											printf("error internal source E2 elm=%lld\n", maxelm);
#else
											printf("error internal source E2 elm=%d\n", maxelm);
#endif
											
											//getchar();
											system("PAUSE");
										}
										else {
											// Внутренний узел.
											G = ESIDE;
											// Внутренняя грань:
											sosedi[G][maxelm].iNODE2 = octree1->linkE2->inum_TD - 1; // граничные КО нумеруются в последнюю очередь,
											sosedi[G][maxelm].bdroblenie4 = true;

											// после maxelm_memo внутренних КО.
											// Вычисление дальнего соседа
											CALC_GG(G, GG); // Вычисляет грань GG по грани G.

											// Определение дальнего соседа :
											// TODO 13 сентября 2016.

										}
									}
								}
							}

							if (octree1->linkE5 == NULL) {

								// 21 сентября 2016 (вырождение ячейки).

								// Заглушка.
								G = ESIDE;
								// граничная грань:
								sosedi[G][maxelm].iNODE3 = NON_EXISTENT_NODE; // граничные КО нумеруются в последнюю очередь,
								sosedi[G][maxelm].bdroblenie4 = true;

								// после maxelm_memo внутренних КО.
								// Вычисление дальнего соседа
								CALC_GG(G, GG); // Вычисляет грань GG по грани G.
								sosedi[GG][maxelm].iNODE3 = NON_EXISTENT_NODE; // соседа нет.
								sosedi[GG][maxelm].bdroblenie4 = true;

								/*
								G = ESIDE;
								// граничная грань:
								sosedi[G][maxelm].iNODE3 = maxelm_memo + maxbound; // граничные КО нумеруются в последнюю очередь,
								sosedi[G][maxelm].bdroblenie4 = true;

								// после maxelm_memo внутренних КО.
								// Вычисление дальнего соседа
								CALC_GG(G,GG); // Вычисляет грань GG по грани G.
								sosedi[GG][maxelm].iNODE3 = NON_EXISTENT_NODE; // соседа нет.
								sosedi[GG][maxelm].bdroblenie4 = true;
								#if doubleintprecision == 1
									printf("error NULL sosedi E5 elm=%lld\n", maxelm);
								#else
									printf("error NULL sosedi E5 elm=%d\n", maxelm);
								#endif
								
								getchar();

								// Эта ситуация возможна и здесь нужно предусмотреть обработку.
								DEFINE_BOUNDARY_PROPERTIES_TEMP(prop, prop_b, maxelm, maxbound); // В граничный узел сносятся свойства прилегающего КО.

								doublereal dS_loc = 0.0;

								// G - грань по которой ставится граничное условие.
								// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
								// bound_id - номер граничного узла, начиная с нуля.
								// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
								if (octree1->b4W) {
									obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
										octree1->linkW0->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
								}
								else {
									if (octree1->linkW == NULL) {
										obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
											 NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
									}
									else {
										obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
											octree1->linkW->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
									}
								}

								maxbound++;
								*/
							}
							else {
								if (octree1->linkE5->inum_TD == 0) {
									G = ESIDE;
									// граничная грань:
									sosedi[G][maxelm].iNODE3 = maxelm_memo + maxbound; // граничные КО нумеруются в последнюю очередь,
									sosedi[G][maxelm].bdroblenie4 = true;

									// после maxelm_memo внутренних КО.
									// Вычисление дальнего соседа
									CALC_GG(G, GG); // Вычисляет грань GG по грани G.
									sosedi[GG][maxelm].iNODE3 = NON_EXISTENT_NODE; // соседа нет.
									sosedi[GG][maxelm].bdroblenie4 = true;
#if doubleintprecision == 1
									//printf("error inum_TD == 0 sosedi E5 elm=%lld\n", maxelm);
#else
									//printf("error inum_TD == 0 sosedi E5 elm=%d\n", maxelm);
#endif
									
									//getchar();

									// Эта ситуация возможна и здесь нужно предусмотреть обработку.
									DEFINE_BOUNDARY_PROPERTIES_TEMP(prop, prop_b, maxelm, maxbound); // В граничный узел сносятся свойства прилегающего КО.
									
									// Sy*Sz
									doublereal dS_loc = fabs(octree1->linkE5->p4.z - octree1->linkE5->p0.z)*fabs(octree1->linkE5->p3.y - octree1->linkE5->p0.y);
									TOCHKA p_centerG;
									p_centerG.z = 0.5*(octree1->linkE5->p4.z + octree1->linkE5->p0.z);
									p_centerG.y = 0.5*(octree1->linkE5->p3.y + octree1->linkE5->p0.y);
									p_centerG.x = octree1->linkE5->p0.x;

									// G - грань по которой ставится граничное условие.
									// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
									// bound_id - номер граничного узла, начиная с нуля.
									// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
									if (octree1->b4W) {
										obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
											octree1->linkW0->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
									}
									else {
										if (octree1->linkW == NULL) {
											obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
												 NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
										}
										else {
											obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
												octree1->linkW->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
										}
									}

									maxbound++;
								}
								else {
									// сосед существует.
									if (bvisit[octree1->linkE5->inum_TD - 1]) {
										// узел уже был посещён
										// Здесь неким образом модифицировалось gran_t
										// Здесь тоже возможен лишь строго внутренний узел.
										// Внутренний узел.
										G = ESIDE;
										// Внутренняя грань:
										sosedi[G][maxelm].iNODE3 = octree1->linkE5->inum_TD - 1; // граничные КО нумеруются в последнюю очередь,
										sosedi[G][maxelm].bdroblenie4 = true;

										// после maxelm_memo внутренних КО.
										// Вычисление дальнего соседа
										CALC_GG(G, GG); // Вычисляет грань GG по грани G.

										// Определение дальнего соседа :
										// TODO 13 сентября 2016.
									}
									else {
										// узнать координаты центра грани и ориентацию в пространстве
										x_c = octree1->p1.x;
										y_c = 0.5*(octree1->linkE5->p0.y + octree1->linkE5->p3.y);
										z_c = 0.5*(octree1->linkE5->p0.z + octree1->linkE5->p4.z);
										iplane = YZ;
										patch_maxbound(iplane, s, ls, x_c, y_c, z_c, maxbound, binc3, lsid3);
										if (binc3) {
#if doubleintprecision == 1
											printf("error internal source E5 elm=%lld\n", maxelm);
#else
											printf("error internal source E5 elm=%d\n", maxelm);
#endif
											
											//getchar();
											system("PAUSE");
										}
										else {
											// Внутренний узел.
											G = ESIDE;
											// Внутренняя грань:
											sosedi[G][maxelm].iNODE3 = octree1->linkE5->inum_TD - 1; // граничные КО нумеруются в последнюю очередь,
											sosedi[G][maxelm].bdroblenie4 = true;

											// после maxelm_memo внутренних КО.
											// Вычисление дальнего соседа
											CALC_GG(G, GG); // Вычисляет грань GG по грани G.

											// Определение дальнего соседа :
											// TODO 13 сентября 2016.
										}
									}
								}
							}

							if (octree1->linkE6 == NULL) {

								// 21 сентября 2016 (вырождение ячейки).

								// Заглушка.
								G = ESIDE;
								// граничная грань:
								sosedi[G][maxelm].iNODE4 = NON_EXISTENT_NODE; // граничные КО нумеруются в последнюю очередь,
								sosedi[G][maxelm].bdroblenie4 = true;

								// после maxelm_memo внутренних КО.
								// Вычисление дальнего соседа
								CALC_GG(G, GG); // Вычисляет грань GG по грани G.
								sosedi[GG][maxelm].iNODE4 = NON_EXISTENT_NODE; // соседа нет.
								sosedi[GG][maxelm].bdroblenie4 = true;


								/*
								G = ESIDE;
								// граничная грань:
								sosedi[G][maxelm].iNODE4 = maxelm_memo + maxbound; // граничные КО нумеруются в последнюю очередь,
								sosedi[G][maxelm].bdroblenie4 = true;

								// после maxelm_memo внутренних КО.
								// Вычисление дальнего соседа
								CALC_GG(G,GG); // Вычисляет грань GG по грани G.
								sosedi[GG][maxelm].iNODE4 = NON_EXISTENT_NODE; // соседа нет.
								sosedi[GG][maxelm].bdroblenie4 = true;
								#if doubleintprecision == 1
										printf("error NULL sosedi E6 elm=%lld\n", maxelm);
								#else
										printf("error NULL sosedi E6 elm=%d\n", maxelm);
								#endif
								
								getchar();

								// Эта ситуация возможна и здесь нужно предусмотреть обработку.
								DEFINE_BOUNDARY_PROPERTIES_TEMP(prop, prop_b, maxelm, maxbound); // В граничный узел сносятся свойства прилегающего КО.
								

								doublereal dS_loc = 0.0;

								// G - грань по которой ставится граничное условие.
								// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
								// bound_id - номер граничного узла, начиная с нуля.
								// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
								if (octree1->b4W) {
									obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
										octree1->linkW0->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
								}
								else {
									if (octree1->linkW == NULL) {
										obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
											 NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
									}
									else {
										obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
											octree1->linkW->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
									}
								}


								maxbound++;
								*/
							}
							else {
								if (octree1->linkE6->inum_TD == 0) {
									G = ESIDE;
									// граничная грань:
									sosedi[G][maxelm].iNODE4 = maxelm_memo + maxbound; // граничные КО нумеруются в последнюю очередь,
									sosedi[G][maxelm].bdroblenie4 = true;

									// после maxelm_memo внутренних КО.
									// Вычисление дальнего соседа
									CALC_GG(G, GG); // Вычисляет грань GG по грани G.
									sosedi[GG][maxelm].iNODE4 = NON_EXISTENT_NODE; // соседа нет.
									sosedi[GG][maxelm].bdroblenie4 = true;
#if doubleintprecision == 1
									//printf("error inum_TD == 0 sosedi E6 elm=%lld\n", maxelm);
#else
									//printf("error inum_TD == 0 sosedi E6 elm=%d\n", maxelm);
#endif
									
									//getchar();

									// Эта ситуация возможна и здесь нужно предусмотреть обработку.
									DEFINE_BOUNDARY_PROPERTIES_TEMP(prop, prop_b, maxelm, maxbound); // В граничный узел сносятся свойства прилегающего КО.
		
									// Sy*Sz
									doublereal dS_loc = fabs(octree1->linkE6->p4.z - octree1->linkE6->p0.z)*fabs(octree1->linkE6->p3.y - octree1->linkE6->p0.y);
									TOCHKA p_centerG;
									p_centerG.z = 0.5*(octree1->linkE6->p4.z + octree1->linkE6->p0.z);
									p_centerG.y = 0.5*(octree1->linkE6->p3.y + octree1->linkE6->p0.y);
									p_centerG.x = octree1->linkE6->p0.x;



									// G - грань по которой ставится граничное условие.
									// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
									// bound_id - номер граничного узла, начиная с нуля.
									// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
									if (octree1->b4W) {
										obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
											octree1->linkW0->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
									}
									else {
										if (octree1->linkW == NULL) {
											obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
												 NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
										}
										else {
											obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
												octree1->linkW->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
										}
									}

									maxbound++;
								}
								else {
									// сосед существует.
									if (bvisit[octree1->linkE6->inum_TD - 1]) {
										// узел уже был посещён
										// Здесь неким образом модифицировалось gran_t
										// Здесь тоже возможен лишь строго внутренний узел.
										// Внутренний узел.
										G = ESIDE;
										// Внутренняя грань:
										sosedi[G][maxelm].iNODE4 = octree1->linkE6->inum_TD - 1; // граничные КО нумеруются в последнюю очередь,
										sosedi[G][maxelm].bdroblenie4 = true;

										// после maxelm_memo внутренних КО.
										// Вычисление дальнего соседа
										CALC_GG(G, GG); // Вычисляет грань GG по грани G.

										// Определение дальнего соседа :
										// TODO 13 сентября 2016.
									}
									else {
										// узнать координаты центра грани и ориентацию в пространстве
										x_c = octree1->p1.x;
										y_c = 0.5*(octree1->linkE6->p0.y + octree1->linkE6->p3.y);
										z_c = 0.5*(octree1->linkE6->p0.z + octree1->linkE6->p4.z);
										iplane = YZ;
										patch_maxbound(iplane, s, ls, x_c, y_c, z_c, maxbound, binc4, lsid4);
										if (binc4) {
#if doubleintprecision == 1
											printf("error internal source E6 elm=%lld\n", maxelm);
#else
											printf("error internal source E6 elm=%d\n", maxelm);
#endif
											
											//getchar();
											system("PAUSE");
										}
										else {
											// Внутренний узел.
											G = ESIDE;
											// Внутренняя грань:
											sosedi[G][maxelm].iNODE4 = octree1->linkE6->inum_TD - 1; // граничные КО нумеруются в последнюю очередь,
											sosedi[G][maxelm].bdroblenie4 = true;

											// после maxelm_memo внутренних КО.
											// Вычисление дальнего соседа
											CALC_GG(G, GG); // Вычисляет грань GG по грани G.

											// Определение дальнего соседа :
											// TODO 13 сентября 2016.
										}
									}
								}
							}
						}
					}

					// сосед не существует
					// записываем идентификатор грани и
					// увеличиваем счётчик граней.
					if (!octree1->b4W) {
						if ((octree1->linkW == NULL)||((octree1->linkW != NULL)&&(octree1->linkW->inum_TD == 0))) {
							// Если мы находимся на границе расчётной области или если мы граничим с HOLLOW блоком.

							G = WSIDE;
							// граничная грань:
							INIT_SOSEDI(sosedi, maxelm, WSIDE, GG, false, maxelm_memo + maxbound,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE, false,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE); // заполнение sosedi[G][maxelm].	
							DEFINE_BOUNDARY_PROPERTIES_TEMP(prop, prop_b, maxelm, maxbound); // В граничный узел сносятся свойства прилегающего КО.
							
							doublereal dSloc = CALC_SQUARE(G, maxelm, nvtx, pa); // Вычисление площади грани ячейки.
							TOCHKA p_centerG = CALC_CENTERG(G, maxelm, nvtx, pa); // Вычисление геометрического центра грани.
							
							// G - грань по которой ставится граничное условие.
							// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
							// bound_id - номер граничного узла, начиная с нуля.
							// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
							if (octree1->b4E) {
								
								obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
									octree1->linkE1->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, false, false, p_centerG);
							}
							else {
								if (octree1->linkE == NULL) {
									obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
										 NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, false, false, p_centerG);
								}
								else {
									obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
										octree1->linkE->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, false, false, p_centerG);
								}
							}

							maxbound++;
						}
						else {
							 if (octree1->linkW != NULL) {
								// сосед существует.
								if (bvisit[octree1->linkW->inum_TD - 1]) {

									// return_current_node_number - Заглянули в соседа и вернулись из соседа в текущий узел.
									// Номер должен быть тем же самым и поэтому этот код бессмысленен.
									// Строго внутренние плоские бесконечно тонкие источники тепла внутри расчётной области больше не поддерживаются.
									// Поэтому от этого кода можно безболезненно избавиться. Поддерживаются только объёмные источники тепла.
									// 10.03.2019
									integer return_current_node_number = sosedi[ESIDE][octree1->linkW->inum_TD - 1].iNODE1;
									if (octree1->linkW->b4E) {
										integer current_node_number = octree1->inum_TD - 1;
										bool bcontinue = true;
										if (bcontinue) {
											if (octree1->linkW->linkE1 != NULL) {
												if (octree1->linkW->linkE1->inum_TD - 1 == current_node_number) {
													bcontinue = false; // iNODE1
													if (sosedi[ESIDE][octree1->linkW->inum_TD - 1].iNODE1 == current_node_number) {
														return_current_node_number = sosedi[ESIDE][octree1->linkW->inum_TD - 1].iNODE1;
													}
													else {
														printf("this can not be W E1\n");
														//getchar();
														system("PAUSE");
													}
												}
											}
										}
										if (bcontinue) {
											if (octree1->linkW->linkE2 != NULL) {
												if (octree1->linkW->linkE2->inum_TD - 1 == current_node_number) {
													bcontinue = false; // iNODE2
													if (sosedi[ESIDE][octree1->linkW->inum_TD - 1].iNODE2 == current_node_number) {
														return_current_node_number = sosedi[ESIDE][octree1->linkW->inum_TD - 1].iNODE2;
													}
													else {
														printf("this can not be W E2\n");
														//getchar();
														system("PAUSE");
													}
												}
											}
										}
										if (bcontinue) {
											if (octree1->linkW->linkE5 != NULL) {
												if (octree1->linkW->linkE5->inum_TD - 1 == current_node_number) {
													bcontinue = false; // iNODE3
													if (sosedi[ESIDE][octree1->linkW->inum_TD - 1].iNODE3 == current_node_number) {
														return_current_node_number = sosedi[ESIDE][octree1->linkW->inum_TD - 1].iNODE3;
													}
													else {
														printf("this can not be W E5\n");
														//getchar();
														system("PAUSE");
													}
												}
											}
										}
										if (bcontinue) {
											if (octree1->linkW->linkE6 != NULL) {
												if (octree1->linkW->linkE6->inum_TD - 1 == current_node_number) {
													bcontinue = false; // iNODE4
													if (sosedi[ESIDE][octree1->linkW->inum_TD - 1].iNODE4 == current_node_number) {
														return_current_node_number = sosedi[ESIDE][octree1->linkW->inum_TD - 1].iNODE4;
													}
													else {
														printf("this can not be W E6\n");
														//getchar();
														system("PAUSE");
													}
												}
											}
										}

									}

									// узел уже был посещён
									// Здесь неким образом модифицировалось gran_t
									if (return_current_node_number >= maxelm_memo) {

										printf("WARNING!!! internal source gran WSIDE 10.03.2019\n");
										system("PAUSE");

										// Это граничная грань внутреннего источника тепла.
										// Внутренний источник тепла.
										G = WSIDE;
										// граничная грань:
										INIT_SOSEDI(sosedi, maxelm, WSIDE, GG, false, sosedi[ESIDE][octree1->linkW->inum_TD - 1].iNODE1,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE, false,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE); // заполнение sosedi[G][maxelm].	
										
										if (bi_fluid) {
											// FLUID
											// В граничный узел сносятся свойства прилегающего КО.
											DEFINE_BOUNDARY_PROPERTIES_TEMP(prop, prop_b, octree1->linkW->inum_TD - 1, sosedi[G][maxelm].iNODE1 - maxelm_memo); 
										}
										else {
											// SOLID
											// В граничный узел сносятся свойства прилегающего КО.
											DEFINE_BOUNDARY_PROPERTIES_TEMP(prop, prop_b, maxelm, sosedi[G][maxelm].iNODE1 - maxelm_memo);
										}

										doublereal dSloc = CALC_SQUARE(G, maxelm, nvtx, pa); // Вычисление площади грани ячейки.
										TOCHKA p_centerG = CALC_CENTERG(G, maxelm, nvtx, pa); // Вычисление геометрического центра грани.
										
										// G - грань по которой ставится граничное условие.
										// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
										// bound_id - номер граничного узла, начиная с нуля.
										// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
										if (octree1->b4E) {
											obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, sosedi[ESIDE][octree1->linkW->inum_TD - 1].iNODE1 - maxelm_memo,
												octree1->linkE1->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, true, true, p_centerG);
										}
										else {
											if (octree1->linkE == NULL) {
												obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, sosedi[ESIDE][octree1->linkW->inum_TD - 1].iNODE1 - maxelm_memo,
													 NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc,true, true);
											}
											else {
												obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, sosedi[ESIDE][octree1->linkW->inum_TD - 1].iNODE1 - maxelm_memo,
													octree1->linkE->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc,true, true);
											}
										}
										binternalsource[sosedi[ESIDE][octree1->linkW->inum_TD - 1].iNODE1 - maxelm_memo] = true; // внутренний источник тепла на границе жидкости и твёрдого тела

									}
									else {
										// Это строго внутренний узел.
										G = WSIDE;
										// внутренняя грань:
										sosedi[G][maxelm].iNODE1 = octree1->linkW->inum_TD - 1; // граничные КО нумеруются в последнюю очередь,
										sosedi[G][maxelm].bdroblenie4 = false;
										sosedi[G][maxelm].iNODE2 = NON_EXISTENT_NODE;
										sosedi[G][maxelm].iNODE3 = NON_EXISTENT_NODE;
										sosedi[G][maxelm].iNODE4 = NON_EXISTENT_NODE;
										// после maxelm_memo внутренних КО.
										// Вычисление дальнего соседа
										//CALC_GG(G, GG); // Вычисляет грань GG по грани G.
										// Определение дальнего соседа :
										// TODO 13 сентября 2016.
									}

								}
								else {
									bool binc;
									integer lsid;
									// узнать координаты центра грани и ориентацию в пространстве
									x_c = octree1->p0.x;
									y_c = 0.5*(octree1->p0.y + octree1->p3.y);
									z_c = 0.5*(octree1->p0.z + octree1->p4.z);
									iplane = YZ;
									patch_maxbound(iplane, s, ls, x_c, y_c, z_c, maxbound,binc, lsid);
									if (binc) {
										// Внутренний источник тепла.
										// Внутренние плоские бесконечно тонкие источники тепла скорее всего больше не поддерживаются. 10,03,2019

										G = WSIDE;
										// граничная грань:
										INIT_SOSEDI(sosedi, maxelm, WSIDE, GG, false, maxelm_memo + maxbound-1,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE, false,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE); // заполнение sosedi[G][maxelm].
										
										if (bi_fluid) {
											// FLUID
											// В граничный узел сносятся свойства прилегающего КО.
											DEFINE_BOUNDARY_PROPERTIES_TEMP(prop, prop_b, octree1->linkW->inum_TD - 1, maxbound-1); 
										}
										else {
											// SOLID
											// В граничный узел сносятся свойства прилегающего КО.
											DEFINE_BOUNDARY_PROPERTIES_TEMP(prop, prop_b, maxelm, maxbound-1); 
										}

										doublereal dSloc = CALC_SQUARE(G, maxelm, nvtx, pa); // Вычисление площади грани ячейки.
										TOCHKA p_centerG = CALC_CENTERG(G, maxelm, nvtx, pa); // Вычисление геометрического центра грани.
										
										// G - грань по которой ставится граничное условие.
										// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
										// bound_id - номер граничного узла, начиная с нуля.
										// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
										if (octree1->b4E) {
											obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound-1,
												octree1->linkE1->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, true, false, p_centerG);
										}
										else {
											if (octree1->linkE == NULL) {
												obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound-1,
													 NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, true, false, p_centerG);
											}
											else {
												obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound-1,
													octree1->linkE->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, true, false, p_centerG);
											}
										}

										binternalsource[maxbound - 1] = true; // внутренний источник тепла на границе жидкости и твёрдого тела

									}
									else {
										// Внутренний узел.
										G = WSIDE;
										// внутренняя грань:
										sosedi[G][maxelm].iNODE1 = octree1->linkW->inum_TD - 1; // граничные КО нумеруются в последнюю очередь,
										sosedi[G][maxelm].bdroblenie4 = false;
										sosedi[G][maxelm].iNODE2 = NON_EXISTENT_NODE;
										sosedi[G][maxelm].iNODE3 = NON_EXISTENT_NODE;
										sosedi[G][maxelm].iNODE4 = NON_EXISTENT_NODE;
										// после maxelm_memo внутренних КО.
										// Вычисление дальнего соседа
										//CALC_GG(G, GG); // Вычисляет грань GG по грани G.
										// Определение дальнего соседа :
										// TODO 13 сентября 2016.

									}
								}
							}
						}
					}
					else {
					// Присутствуют все 4 соседа: W: 0,3,4,7.

						if (((octree1->linkW0 == NULL) && (octree1->linkW3 == NULL) && (octree1->linkW4 == NULL) && (octree1->linkW7 == NULL))||
						(((octree1->linkW0 != NULL) && (octree1->linkW3 != NULL) && (octree1->linkW4 != NULL) && (octree1->linkW7 != NULL) &&
							(octree1->linkW0->inum_TD == 0) && (octree1->linkW3->inum_TD == 0) && (octree1->linkW4->inum_TD == 0) && (octree1->linkW7->inum_TD == 0)))||
							(((octree1->linkW0 != NULL) && (octree1->linkW3 == NULL) && (octree1->linkW4 != NULL) && (octree1->linkW7 == NULL) &&
							(octree1->linkW0->inum_TD == 0) && (octree1->linkW4->inum_TD == 0)))||
							(((octree1->linkW0 != NULL) && (octree1->linkW3 != NULL) && (octree1->linkW4 == NULL) && (octree1->linkW7 == NULL) &&
							(octree1->linkW0->inum_TD == 0) && (octree1->linkW3->inum_TD == 0)))||
							(((octree1->linkW0 != NULL) && (octree1->linkW3 == NULL) && (octree1->linkW4 == NULL) && (octree1->linkW7 == NULL) &&
							(octree1->linkW0->inum_TD == 0))))
						{
							// 1 - W0&&W3&&W4&&W7 - на границе расчётной области (РО).
							// 2 - W0&&W3&&W4&&W7 - на гранеице HOLLOW блока.
							// 3 - W0&&W4 на границе HOLLOW блока, W3&&W7 на границе РО.
							// 4 - W0&&W3 на границе HOLLOW блока, W4&&W7 на границе РО.
							// 5 - W0 на границе HOLLOW блока, W3&&W4&&W7  на границе РО.

							G = WSIDE;
							// граничная грань:
							INIT_SOSEDI(sosedi, maxelm, WSIDE, GG, false, maxelm_memo + maxbound,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE, false,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE); // заполнение sosedi[G][maxelm].
							DEFINE_BOUNDARY_PROPERTIES_TEMP(prop, prop_b, maxelm, maxbound); // В граничный узел сносятся свойства прилегающего КО.

							doublereal dSloc = CALC_SQUARE(G, maxelm, nvtx, pa); // Вычисление площади грани ячейки.
							TOCHKA p_centerG = CALC_CENTERG(G, maxelm, nvtx, pa); // Вычисление геометрического центра грани.
							
							// G - грань по которой ставится граничное условие.
							// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
							// bound_id - номер граничного узла, начиная с нуля.
							// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
							if (octree1->b4E) {
								obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
									octree1->linkE1->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, false, false, p_centerG);
							}
							else {
								if (octree1->linkE == NULL) {
									obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
										 NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, false, false, p_centerG);
								}
								else {
									obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
										octree1->linkE->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, false, false, p_centerG);
								}
							}

							// Зачем дробить если граница заканчивается NULL.
							maxbound++; // Вместо четырёх один.
						}						
						else {

							bool binc1, binc2, binc3, binc4;
							integer lsid1, lsid2, lsid3, lsid4;

							if (octree1->linkW0 == NULL) {

#if doubleintprecision == 1
								printf("error NULL sosedi W0 elm=%lld\n", maxelm);
#else
								printf("error NULL sosedi W0 elm=%d\n", maxelm);
#endif
								
								//getchar();
								system("PAUSE");

								/*
								G = WSIDE;
								// граничная грань:
								sosedi[G][maxelm].iNODE1 = maxelm_memo + maxbound; // граничные КО нумеруются в последнюю очередь,
								sosedi[G][maxelm].bdroblenie4 = true;
								// после maxelm_memo внутренних КО.
								// Вычисление дальнего соседа
								CALC_GG(G,GG); // Вычисляет грань GG по грани G.
								sosedi[GG][maxelm].iNODE1 = NON_EXISTENT_NODE; // соседа нет.
								sosedi[GG][maxelm].bdroblenie4 = true;
								#if doubleintprecision == 1
									printf("error NULL sosedi W0 elm=%lld\n", maxelm);
								#else
									printf("error NULL sosedi W0 elm=%d\n", maxelm);
								#endif
								
								getchar();

								// Эта ситуация возможна и здесь надо предусмотреть её обработку.
								DEFINE_BOUNDARY_PROPERTIES_TEMP(prop, prop_b, maxelm, maxbound); // В граничный узел сносятся свойства прилегающего КО.
								
								doublereal dS_loc = 0.0;

								// G - грань по которой ставится граничное условие.
								// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
								// bound_id - номер граничного узла, начиная с нуля.
								// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
								if (octree1->b4E) {
									obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
										octree1->linkE1->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
								}
								else {
									if (octree1->linkE == NULL) {
										obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
											 NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
									}
									else {
										obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
											octree1->linkE->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
									}
								}

								maxbound++;
								*/
							}
							else {
								if (octree1->linkW0->inum_TD == 0) {
									G = WSIDE;
									// граничная грань:
									sosedi[G][maxelm].iNODE1 = maxelm_memo + maxbound; // граничные КО нумеруются в последнюю очередь,
									sosedi[G][maxelm].bdroblenie4 = true;
									// после maxelm_memo внутренних КО.
									// Вычисление дальнего соседа
									CALC_GG(G, GG); // Вычисляет грань GG по грани G.
									sosedi[GG][maxelm].iNODE1 = NON_EXISTENT_NODE; // соседа нет.
									sosedi[GG][maxelm].bdroblenie4 = true;
#if doubleintprecision == 1
									//printf("error inum_TD == 0 sosedi W0 elm=%lld\n", maxelm);
#else
									//printf("error inum_TD == 0 sosedi W0 elm=%d\n", maxelm);
#endif
									
									//getchar();

									// Эта ситуация возможна и здесь надо предусмотреть её обработку.
									DEFINE_BOUNDARY_PROPERTIES_TEMP(prop, prop_b, maxelm, maxbound); // В граничный узел сносятся свойства прилегающего КО.
									
									// Sy*Sz
									doublereal dS_loc = fabs(octree1->linkW0->p5.z - octree1->linkW0->p1.z)*fabs(octree1->linkW0->p2.y - octree1->linkW0->p1.y);
									TOCHKA p_centerG;
									p_centerG.z = 0.5*(octree1->linkW0->p5.z + octree1->linkW0->p1.z);
									p_centerG.y = 0.5*(octree1->linkW0->p2.y + octree1->linkW0->p1.y);
									p_centerG.x = octree1->linkW0->p1.x;

									// G - грань по которой ставится граничное условие.
									// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
									// bound_id - номер граничного узла, начиная с нуля.
									// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
									if (octree1->b4E) {
										obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
											octree1->linkE1->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
									}
									else {
										if (octree1->linkE == NULL) {
											obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
												 NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
										}
										else {
											obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
												octree1->linkE->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
										}
									}

									maxbound++;
								}
								else {
									// сосед существует.
									if (bvisit[octree1->linkW0->inum_TD - 1]) {
										// узел уже был посещён
										// Здесь неким образом модифицировалось gran_t
										// Только внутренняя грань.
										G = WSIDE;
										// граничная грань:
										sosedi[G][maxelm].iNODE1 = octree1->linkW0->inum_TD - 1; // граничные КО нумеруются в последнюю очередь,
										sosedi[G][maxelm].bdroblenie4 = true;
										// после maxelm_memo внутренних КО.
										// Вычисление дальнего соседа
										CALC_GG(G, GG); // Вычисляет грань GG по грани G.

										// Определение дальнего соседа :
										// TODO 13 сентября 2016.

									}
									else {
										// узнать координаты центра грани и ориентацию в пространстве
										x_c = octree1->p0.x;
										y_c = 0.5*(octree1->linkW0->p1.y + octree1->linkW0->p2.y);
										z_c = 0.5*(octree1->linkW0->p1.z + octree1->linkW0->p5.z);
										iplane = YZ;
										patch_maxbound(iplane, s, ls, x_c, y_c, z_c, maxbound, binc1, lsid1);
										if (binc1) {
#if doubleintprecision == 1
											printf("error internal source W0 elm=%lld\n", maxelm);
#else
											printf("error internal source W0 elm=%d\n", maxelm);
#endif
											
											//getchar();
											system("PAUSE");
										}
										else {
											// Только внутренняя грань.
											G = WSIDE;
											// граничная грань:
											sosedi[G][maxelm].iNODE1 = octree1->linkW0->inum_TD - 1; // граничные КО нумеруются в последнюю очередь,
											sosedi[G][maxelm].bdroblenie4 = true;
											// после maxelm_memo внутренних КО.
											// Вычисление дальнего соседа
											CALC_GG(G, GG); // Вычисляет грань GG по грани G.

											// Определение дальнего соседа :
											// TODO 13 сентября 2016.
										}
									}
								}
							}

							if (octree1->linkW3 == NULL) {

								// 21 сентября 2016 (вырождение ячейки).

								// Заглушка.
								G = WSIDE;
								// граничная грань:
								sosedi[G][maxelm].iNODE2 = NON_EXISTENT_NODE; // граничные КО нумеруются в последнюю очередь,
								sosedi[G][maxelm].bdroblenie4 = true;
								// после maxelm_memo внутренних КО.
								// Вычисление дальнего соседа
								CALC_GG(G, GG); // Вычисляет грань GG по грани G.
								sosedi[GG][maxelm].iNODE2 = NON_EXISTENT_NODE; // соседа нет.
								sosedi[GG][maxelm].bdroblenie4 = true;

								/*
								G = WSIDE;
								// граничная грань:
								sosedi[G][maxelm].iNODE2 = maxelm_memo + maxbound; // граничные КО нумеруются в последнюю очередь,
								sosedi[G][maxelm].bdroblenie4 = true;
								// после maxelm_memo внутренних КО.
								// Вычисление дальнего соседа
								CALC_GG(G,GG); // Вычисляет грань GG по грани G.
								sosedi[GG][maxelm].iNODE2 = NON_EXISTENT_NODE; // соседа нет.
								sosedi[GG][maxelm].bdroblenie4 = true;

								#if doubleintprecision == 1
									printf("error NULL sosedi W3 elm=%lld\n", maxelm);
								#else
									printf("error NULL sosedi W3 elm=%d\n", maxelm);
								#endif
								
								getchar();

								// Эта ситуация возможна и здесь надо предусмотреть её обработку.
								DEFINE_BOUNDARY_PROPERTIES_TEMP(prop, prop_b, maxelm, maxbound); // В граничный узел сносятся свойства прилегающего КО.
								
								doublereal dS_loc = 0.0;

								// G - грань по которой ставится граничное условие.
								// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
								// bound_id - номер граничного узла, начиная с нуля.
								// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
								if (octree1->b4E) {
									obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
										octree1->linkE1->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
								}
								else {
									if (octree1->linkE == NULL) {
										obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
											 NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
									}
									else {
										obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
											octree1->linkE->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
									}
								}

								maxbound++;
								*/

							}
							else {
								if (octree1->linkW3->inum_TD == 0) {
									G = WSIDE;
									// граничная грань:
									sosedi[G][maxelm].iNODE2 = maxelm_memo + maxbound; // граничные КО нумеруются в последнюю очередь,
									sosedi[G][maxelm].bdroblenie4 = true;
									// после maxelm_memo внутренних КО.
									// Вычисление дальнего соседа
									CALC_GG(G, GG); // Вычисляет грань GG по грани G.
									sosedi[GG][maxelm].iNODE2 = NON_EXISTENT_NODE; // соседа нет.
									sosedi[GG][maxelm].bdroblenie4 = true;
#if doubleintprecision == 1
									//printf("error inum_TD == 0 sosedi W3 elm=%lld\n", maxelm);
#else
									//printf("error inum_TD == 0 sosedi W3 elm=%d\n", maxelm);
#endif
									
									//getchar();

									// Эта ситуация возможна и здесь надо предусмотреть её обработку.
									DEFINE_BOUNDARY_PROPERTIES_TEMP(prop, prop_b, maxelm, maxbound); // В граничный узел сносятся свойства прилегающего КО.
									
									// Sy*Sz
									doublereal dS_loc = fabs(octree1->linkW3->p5.z - octree1->linkW3->p1.z)*fabs(octree1->linkW3->p2.y - octree1->linkW3->p1.y);
									TOCHKA p_centerG;
									p_centerG.z = 0.5*(octree1->linkW3->p5.z + octree1->linkW3->p1.z);
									p_centerG.y = 0.5*(octree1->linkW3->p2.y + octree1->linkW3->p1.y);
									p_centerG.x = octree1->linkW3->p1.x;


									// G - грань по которой ставится граничное условие.
									// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
									// bound_id - номер граничного узла, начиная с нуля.
									// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
									if (octree1->b4E) {
										obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
											octree1->linkE1->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
									}
									else {
										if (octree1->linkE == NULL) {
											obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
												 NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
										}
										else {
											obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
												octree1->linkE->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
										}
									}


									maxbound++;
								}
								else {
									// сосед существует.
									if (bvisit[octree1->linkW3->inum_TD - 1]) {
										// узел уже был посещён
										// Здесь неким образом модифицировалось gran_t
										// Только внутренняя грань.
										G = WSIDE;
										// граничная грань:
										sosedi[G][maxelm].iNODE2 = octree1->linkW3->inum_TD - 1; // граничные КО нумеруются в последнюю очередь,
										sosedi[G][maxelm].bdroblenie4 = true;
										// после maxelm_memo внутренних КО.
										// Вычисление дальнего соседа
										CALC_GG(G, GG); // Вычисляет грань GG по грани G.

										// Определение дальнего соседа :
										// TODO 13 сентября 2016.
									}
									else {
										// узнать координаты центра грани и ориентацию в пространстве
										x_c = octree1->p0.x;
										y_c = 0.5*(octree1->linkW3->p1.y + octree1->linkW3->p2.y);
										z_c = 0.5*(octree1->linkW3->p1.z + octree1->linkW3->p5.z);
										iplane = YZ;
										patch_maxbound(iplane, s, ls, x_c, y_c, z_c, maxbound, binc2, lsid2);
										if (binc2) {
#if doubleintprecision == 1
											printf("error internal source W3 elm=%lld\n", maxelm);
#else
											printf("error internal source W3 elm=%d\n", maxelm);
#endif
											
											//getchar();
											system("PAUSE");
										}
										else {
											// Только внутренняя грань.
											G = WSIDE;
											// граничная грань:
											sosedi[G][maxelm].iNODE2 = octree1->linkW3->inum_TD - 1; // граничные КО нумеруются в последнюю очередь,
											sosedi[G][maxelm].bdroblenie4 = true;
											// после maxelm_memo внутренних КО.
											// Вычисление дальнего соседа
											CALC_GG(G, GG); // Вычисляет грань GG по грани G.

											// Определение дальнего соседа :
											// TODO 13 сентября 2016.
										}
									}
								}
							}

							if (octree1->linkW4 == NULL) {

								// 21 сентября 2016 (вырождение ячейки).

								// Заглушка.
								G = WSIDE;
								// граничная грань:
								sosedi[G][maxelm].iNODE3 = NON_EXISTENT_NODE; // граничные КО нумеруются в последнюю очередь,
								sosedi[G][maxelm].bdroblenie4 = true;
								// после maxelm_memo внутренних КО.
								// Вычисление дальнего соседа
								CALC_GG(G, GG); // Вычисляет грань GG по грани G.
								sosedi[GG][maxelm].iNODE3 = NON_EXISTENT_NODE; // соседа нет.
								sosedi[GG][maxelm].bdroblenie4 = true;

								/*
								G = WSIDE;
								// граничная грань:
								sosedi[G][maxelm].iNODE3 = maxelm_memo + maxbound; // граничные КО нумеруются в последнюю очередь,
								sosedi[G][maxelm].bdroblenie4 = true;
								// после maxelm_memo внутренних КО.
								// Вычисление дальнего соседа
								CALC_GG(G,GG); // Вычисляет грань GG по грани G.
								sosedi[GG][maxelm].iNODE3 = NON_EXISTENT_NODE; // соседа нет.
								sosedi[GG][maxelm].bdroblenie4 = true;

								#if doubleintprecision == 1
										printf("error NULL sosedi W4 elm=%lld\n", maxelm);
								#else
										printf("error NULL sosedi W4 elm=%d\n", maxelm);
								#endif
								
								getchar();

								// Эта ситуация возможна и здесь надо предусмотреть её обработку.
								DEFINE_BOUNDARY_PROPERTIES_TEMP(prop, prop_b, maxelm, maxbound); // В граничный узел сносятся свойства прилегающего КО.
								
								doublereal dS_loc = 0.0;

								// G - грань по которой ставится граничное условие.
								// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
								// bound_id - номер граничного узла, начиная с нуля.
								// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
								if (octree1->b4E) {
									obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
										octree1->linkE1->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
								}
								else {
									if (octree1->linkE == NULL) {
										obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
											 NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
									}
									else {
										obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
											octree1->linkE->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
									}
								}

								maxbound++;
								*/

							}
							else {
								if (octree1->linkW4->inum_TD == 0) {
									G = WSIDE;
									// граничная грань:
									sosedi[G][maxelm].iNODE3 = maxelm_memo + maxbound; // граничные КО нумеруются в последнюю очередь,
									sosedi[G][maxelm].bdroblenie4 = true;
									// после maxelm_memo внутренних КО.
									// Вычисление дальнего соседа
									CALC_GG(G, GG); // Вычисляет грань GG по грани G.
									sosedi[GG][maxelm].iNODE3 = NON_EXISTENT_NODE; // соседа нет.
									sosedi[GG][maxelm].bdroblenie4 = true;

#if doubleintprecision == 1
									//printf("error inum_TD == 0 sosedi W4 elm=%lld\n", maxelm);
#else
									//printf("error inum_TD == 0 sosedi W4 elm=%d\n", maxelm);
#endif
									
									//getchar();

									// Эта ситуация возможна и здесь надо предусмотреть её обработку.
									DEFINE_BOUNDARY_PROPERTIES_TEMP(prop, prop_b, maxelm, maxbound); // В граничный узел сносятся свойства прилегающего КО.
									

									// Sy*Sz
									doublereal dS_loc = fabs(octree1->linkW4->p5.z - octree1->linkW4->p1.z)*fabs(octree1->linkW4->p2.y - octree1->linkW4->p1.y);
									TOCHKA p_centerG;
									p_centerG.z = 0.5*(octree1->linkW4->p5.z + octree1->linkW4->p1.z);
									p_centerG.y = 0.5*(octree1->linkW4->p2.y + octree1->linkW4->p1.y);
									p_centerG.x = octree1->linkW4->p1.x;


									// G - грань по которой ставится граничное условие.
									// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
									// bound_id - номер граничного узла, начиная с нуля.
									// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
									if (octree1->b4E) {
										obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
											octree1->linkE1->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
									}
									else {
										if (octree1->linkE == NULL) {
											obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
												 NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
										}
										else {
											obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
												octree1->linkE->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
										}
									}

									maxbound++;
								}
								else {
									// сосед существует.
									if (bvisit[octree1->linkW4->inum_TD - 1]) {
										// узел уже был посещён
										// Здесь неким образом модифицировалось gran_t
										// Только внутренняя грань.
										G = WSIDE;
										// граничная грань:
										sosedi[G][maxelm].iNODE3 = octree1->linkW4->inum_TD - 1; // граничные КО нумеруются в последнюю очередь,
										sosedi[G][maxelm].bdroblenie4 = true;
										// после maxelm_memo внутренних КО.
										// Вычисление дальнего соседа
										CALC_GG(G, GG); // Вычисляет грань GG по грани G.

										// Определение дальнего соседа :
										// TODO 13 сентября 2016.
									}
									else {
										// узнать координаты центра грани и ориентацию в пространстве
										x_c = octree1->p0.x;
										y_c = 0.5*(octree1->linkW4->p1.y + octree1->linkW4->p2.y);
										z_c = 0.5*(octree1->linkW4->p1.z + octree1->linkW4->p5.z);
										iplane = YZ;
										patch_maxbound(iplane, s, ls, x_c, y_c, z_c, maxbound, binc3, lsid3);
										if (binc3) {
#if doubleintprecision == 1
											printf("error internal source W4 elm=%lld\n", maxelm);
#else
											printf("error internal source W4 elm=%d\n", maxelm);
#endif
											
											//getchar();
											system("PAUSE");
										}
										else {
											// Только внутренняя грань.
											G = WSIDE;
											// граничная грань:
											sosedi[G][maxelm].iNODE3 = octree1->linkW4->inum_TD - 1; // граничные КО нумеруются в последнюю очередь,
											sosedi[G][maxelm].bdroblenie4 = true;
											// после maxelm_memo внутренних КО.
											// Вычисление дальнего соседа
											CALC_GG(G, GG); // Вычисляет грань GG по грани G.

											// Определение дальнего соседа :
											// TODO 13 сентября 2016.
										}
									}
								}
							}

							if (octree1->linkW7 == NULL) {

								// 21 сентября 2016 (вырождение ячейки).

								// Заглушка.
								G = WSIDE;
								// граничная грань:
								sosedi[G][maxelm].iNODE4 = NON_EXISTENT_NODE; // граничные КО нумеруются в последнюю очередь,
								sosedi[G][maxelm].bdroblenie4 = true;
								// после maxelm_memo внутренних КО.
								// Вычисление дальнего соседа
								CALC_GG(G, GG); // Вычисляет грань GG по грани G.
								sosedi[GG][maxelm].iNODE4 = NON_EXISTENT_NODE; // соседа нет.
								sosedi[GG][maxelm].bdroblenie4 = true;

								/*
								G = WSIDE;
								// граничная грань:
								sosedi[G][maxelm].iNODE4 = maxelm_memo + maxbound; // граничные КО нумеруются в последнюю очередь,
								sosedi[G][maxelm].bdroblenie4 = true;
								// после maxelm_memo внутренних КО.
								// Вычисление дальнего соседа
								CALC_GG(G,GG); // Вычисляет грань GG по грани G.
								sosedi[GG][maxelm].iNODE4 = NON_EXISTENT_NODE; // соседа нет.
								sosedi[GG][maxelm].bdroblenie4 = true;

								#if doubleintprecision == 1
										printf("error NULL sosedi W7 elm=%lld\n", maxelm);
								#else
										printf("error NULL sosedi W7 elm=%d\n", maxelm);
								#endif
								
								getchar();

								// Эта ситуация возможна и здесь надо предусмотреть её обработку.
								DEFINE_BOUNDARY_PROPERTIES_TEMP(prop, prop_b, maxelm, maxbound); // В граничный узел сносятся свойства прилегающего КО.
								
								doublereal dS_loc = 0.0;

								// G - грань по которой ставится граничное условие.
								// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
								// bound_id - номер граничного узла, начиная с нуля.
								// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
								if (octree1->b4E) {
									obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
										octree1->linkE1->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
								}
								else {
									if (octree1->linkE == NULL) {
										obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
											 NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
									}
									else {
										obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
											octree1->linkE->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
									}
								}

								maxbound++;
								*/

							}
							else {
								if (octree1->linkW7->inum_TD == 0) {
									G = WSIDE;
									// граничная грань:
									sosedi[G][maxelm].iNODE4 = maxelm_memo + maxbound; // граничные КО нумеруются в последнюю очередь,
									sosedi[G][maxelm].bdroblenie4 = true;
									// после maxelm_memo внутренних КО.
									// Вычисление дальнего соседа
									CALC_GG(G, GG); // Вычисляет грань GG по грани G.
									sosedi[GG][maxelm].iNODE4 = NON_EXISTENT_NODE; // соседа нет.
									sosedi[GG][maxelm].bdroblenie4 = true;

#if doubleintprecision == 1
									//printf("error inum_TD == 0 sosedi W7 elm=%lld\n", maxelm);
#else
									//printf("error inum_TD == 0 sosedi W7 elm=%d\n", maxelm);
#endif
									
									//getchar();

									// Эта ситуация возможна и здесь надо предусмотреть её обработку.
									DEFINE_BOUNDARY_PROPERTIES_TEMP(prop, prop_b, maxelm, maxbound); // В граничный узел сносятся свойства прилегающего КО.
									
									// Sy*Sz
									doublereal dS_loc = fabs(octree1->linkW7->p5.z - octree1->linkW7->p1.z)*fabs(octree1->linkW7->p2.y - octree1->linkW7->p1.y);
									TOCHKA p_centerG;
									p_centerG.z = 0.5*(octree1->linkW7->p5.z + octree1->linkW7->p1.z);
									p_centerG.y = 0.5*(octree1->linkW7->p2.y + octree1->linkW7->p1.y);
									p_centerG.x = octree1->linkW7->p1.x;

									// G - грань по которой ставится граничное условие.
									// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
									// bound_id - номер граничного узла, начиная с нуля.
									// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
									if (octree1->b4E) {
										obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
											octree1->linkE1->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
									}
									else {
										if (octree1->linkE == NULL) {
											obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
												 NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
										}
										else {
											obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
												octree1->linkE->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
										}
									}

									maxbound++;
								}
								else {
									// сосед существует.
									if (bvisit[octree1->linkW7->inum_TD - 1]) {
										// узел уже был посещён
										// Здесь неким образом модифицировалось gran_t
										// Только внутренняя грань.
										G = WSIDE;
										// граничная грань:
										sosedi[G][maxelm].iNODE4 = octree1->linkW7->inum_TD - 1; // граничные КО нумеруются в последнюю очередь,
										sosedi[G][maxelm].bdroblenie4 = true;
										// после maxelm_memo внутренних КО.
										// Вычисление дальнего соседа
										CALC_GG(G, GG); // Вычисляет грань GG по грани G.

										// Определение дальнего соседа :
										// TODO 13 сентября 2016.
									}
									else {
										// узнать координаты центра грани и ориентацию в пространстве
										x_c = octree1->p0.x;
										y_c = 0.5*(octree1->linkW7->p1.y + octree1->linkW7->p2.y);
										z_c = 0.5*(octree1->linkW7->p1.z + octree1->linkW7->p5.z);
										iplane = YZ;
										patch_maxbound(iplane, s, ls, x_c, y_c, z_c, maxbound, binc4, lsid4);
										if (binc4) {
#if doubleintprecision == 1
											printf("error internal source W7 elm=%lld\n", maxelm);
#else
											printf("error internal source W7 elm=%d\n", maxelm);
#endif
											
											//getchar();
											system("PAUSE");
										}
										else {
											// Только внутренняя грань.
											G = WSIDE;
											// граничная грань:
											sosedi[G][maxelm].iNODE4 = octree1->linkW7->inum_TD - 1; // граничные КО нумеруются в последнюю очередь,
											sosedi[G][maxelm].bdroblenie4 = true;
											// после maxelm_memo внутренних КО.
											// Вычисление дальнего соседа
											CALC_GG(G, GG); // Вычисляет грань GG по грани G.

											// Определение дальнего соседа :
											// TODO 13 сентября 2016.
										}
									}
								}
							}
						}
					}


					// сосед не существует
					// записываем идентификатор грани и
					// увеличиваем счётчик граней.
					if (!octree1->b4N) {
						if ((octree1->linkN == NULL)||((octree1->linkN != NULL)&&(octree1->linkN->inum_TD == 0))) {
							// Если мы находимся на границе расчётной области или если мы граничим с HOLLOW блоком.

							G = NSIDE;
							// граничная грань:
							INIT_SOSEDI(sosedi, maxelm, NSIDE, GG, false, maxelm_memo + maxbound,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE, false,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE); // заполнение sosedi[G][maxelm].
							DEFINE_BOUNDARY_PROPERTIES_TEMP(prop, prop_b, maxelm, maxbound); // В граничный узел сносятся свойства прилегающего КО.

							doublereal dSloc = CALC_SQUARE(G, maxelm, nvtx, pa); // Вычисление площади грани ячейки.
							TOCHKA p_centerG = CALC_CENTERG(G, maxelm, nvtx, pa); // Вычисление геометрического центра грани.
							
							// G - грань по которой ставится граничное условие.
							// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
							// bound_id - номер граничного узла, начиная с нуля.
							// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
							if (octree1->b4S) {
								obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
									octree1->linkS0->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, false, false, p_centerG);
							}
							else {
								if (octree1->linkS == NULL) {
									obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
										 NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, false, false, p_centerG);
								}
								else {
									obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
										octree1->linkS->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, false, false, p_centerG);
								}
							}

							maxbound++;
						}
						else {
							if (octree1->linkN != NULL) {
								// сосед существует.
								if (bvisit[octree1->linkN->inum_TD - 1]) {

									// return_current_node_number - Заглянули в соседа и вернулись из соседа в текущий узел.
									// Номер должен быть тем же самым и поэтому этот код бессмысленен.
									// Строго внутренние плоские бесконечно тонкие источники тепла внутри расчётной области больше не поддерживаются.
									// Поэтому от этого кода можно безболезненно избавиться. Поддерживаются только объёмные источники тепла.
									// 10.03.2019
									integer return_current_node_number = sosedi[SSIDE][octree1->linkN->inum_TD - 1].iNODE1;
									//S0 NODE1,
									//S1 NODE2,
									//S4 NODE3,
									//S5 NODE4.
									// TODO 23 сентября 2016.
									if (octree1->linkN->b4S) {
										integer current_node_number = octree1->inum_TD - 1;
										bool bcontinue = true;
										if (bcontinue) {
											if (octree1->linkN->linkS0 != NULL) {
												if (octree1->linkN->linkS0->inum_TD - 1 == current_node_number) {
													bcontinue = false; // iNODE1
													if (sosedi[SSIDE][octree1->linkN->inum_TD - 1].iNODE1 == current_node_number) {
														return_current_node_number = sosedi[SSIDE][octree1->linkN->inum_TD - 1].iNODE1;
													}
													else {
														printf("this can not be N S0\n");
														//getchar();
														system("PAUSE");
													}
												}
											}
										}
										if (bcontinue) {
											if (octree1->linkN->linkS1 != NULL) {
												if (octree1->linkN->linkS1->inum_TD - 1 == current_node_number) {
													bcontinue = false; // iNODE2
													if (sosedi[SSIDE][octree1->linkN->inum_TD - 1].iNODE2 == current_node_number) {
														return_current_node_number = sosedi[SSIDE][octree1->linkN->inum_TD - 1].iNODE2;
													}
													else {
														printf("this can not be N S1\n");
														//getchar();
														system("PAUSE");
													}
												}
											}
										}
										if (bcontinue) {
											if (octree1->linkN->linkS4 != NULL) {
												if (octree1->linkN->linkS4->inum_TD - 1 == current_node_number) {
													bcontinue = false; // iNODE3
													if (sosedi[SSIDE][octree1->linkN->inum_TD - 1].iNODE3 == current_node_number) {
														return_current_node_number = sosedi[SSIDE][octree1->linkN->inum_TD - 1].iNODE3;
													}
													else {
														printf("this can not be N S4\n");
														//getchar();
														system("PAUSE");
													}
												}
											}
										}
										if (bcontinue) {
											if (octree1->linkN->linkS5 != NULL) {
												if (octree1->linkN->linkS5->inum_TD - 1 == current_node_number) {
													bcontinue = false; // iNODE4
													if (sosedi[SSIDE][octree1->linkN->inum_TD - 1].iNODE4 == current_node_number) {
														return_current_node_number = sosedi[SSIDE][octree1->linkN->inum_TD - 1].iNODE4;
													}
													else {
														printf("this can not be N S5\n");
														//getchar();
														system("PAUSE");
													}
												}
											}
										}

									}

									// узел уже был посещён
									// Здесь неким образом модифицировалось gran_t
									if (return_current_node_number >= maxelm_memo) {

										printf("WARNING!!! internal source gran NSIDE 10.03.2019\n");
										system("PAUSE");

										// Это граничная грань внутреннего источника тепла.
										// Внутренний источник тепла.
										G = NSIDE;
										// граничная грань:
										INIT_SOSEDI(sosedi, maxelm, NSIDE, GG, false, sosedi[SSIDE][octree1->linkN->inum_TD - 1].iNODE1,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE, false,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE); // заполнение sosedi[G][maxelm].
										
										if (bi_fluid) {
											// FLUID
											// В граничный узел сносятся свойства прилегающего КО.
											DEFINE_BOUNDARY_PROPERTIES_TEMP(prop, prop_b, octree1->linkN->inum_TD - 1, sosedi[G][maxelm].iNODE1 - maxelm_memo);
										}
										else {
											// SOLID
											// В граничный узел сносятся свойства прилегающего КО.
											DEFINE_BOUNDARY_PROPERTIES_TEMP(prop, prop_b, maxelm, sosedi[G][maxelm].iNODE1 - maxelm_memo);
										}

										doublereal dSloc = CALC_SQUARE(G, maxelm, nvtx, pa); // Вычисление площади грани ячейки.
										TOCHKA p_centerG = CALC_CENTERG(G, maxelm, nvtx, pa); // Вычисление геометрического центра грани.
										
										// G - грань по которой ставится граничное условие.
										// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
										// bound_id - номер граничного узла, начиная с нуля.
										// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
										if (octree1->b4S) {
											obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, sosedi[SSIDE][octree1->linkN->inum_TD - 1].iNODE1 - maxelm_memo,
												octree1->linkS0->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, true, true, p_centerG);
										}
										else {
											if (octree1->linkS == NULL) {
												obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, sosedi[SSIDE][octree1->linkN->inum_TD - 1].iNODE1 - maxelm_memo,
													 NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, true, true, p_centerG);
											}
											else {
												obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, sosedi[SSIDE][octree1->linkN->inum_TD - 1].iNODE1 - maxelm_memo,
													octree1->linkS->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, true, true, p_centerG);
											}
										}
										binternalsource[sosedi[SSIDE][octree1->linkN->inum_TD - 1].iNODE1 - maxelm_memo] = true; // внутренний источник тепла на границе жидкости и твёрдого тела

									}
									else {
										// Это строго внутренний узел.
										// Внутренний узел.
										G = NSIDE;
										// внутренняя грань:
										sosedi[G][maxelm].iNODE1 = octree1->linkN->inum_TD - 1; // граничные КО нумеруются в последнюю очередь,
										sosedi[G][maxelm].bdroblenie4 = false;
										sosedi[G][maxelm].iNODE2 = NON_EXISTENT_NODE;
										sosedi[G][maxelm].iNODE3 = NON_EXISTENT_NODE;
										sosedi[G][maxelm].iNODE4 = NON_EXISTENT_NODE;
										// после maxelm_memo внутренних КО.
										// Вычисление дальнего соседа.
										//CALC_GG(G, GG); // Вычисляет грань GG по грани G.
										// Определение дальнего соседа :
										// TODO 13 сентября 2016.
									}

								}
								else {
									bool binc;
									integer lsid;
									// узнать координаты центра грани и ориентацию в пространстве
									x_c = 0.5*(octree1->p2.x + octree1->p3.x);
									y_c = octree1->p3.y;
									z_c = 0.5*(octree1->p3.z + octree1->p7.z);
									iplane = XZ;
									patch_maxbound(iplane, s, ls, x_c, y_c, z_c, maxbound,binc,lsid);
									if (binc) {
										// Внутренний источник тепла.
										// Внутренние плоские бесконечно тонкие источники тепла скорее всего больше не поддерживаются. 10,03,2019

										G = NSIDE;
										// граничная грань:
										INIT_SOSEDI(sosedi, maxelm, NSIDE, GG, false, maxelm_memo + maxbound -1,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE, false,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE); // заполнение sosedi[G][maxelm].
										
										if (bi_fluid) {
											// FLUID
											// В граничный узел сносятся свойства прилегающего КО.
											DEFINE_BOUNDARY_PROPERTIES_TEMP(prop, prop_b, octree1->linkN->inum_TD - 1, maxbound-1);
										}
										else {
											// SOLID
											// В граничный узел сносятся свойства прилегающего КО.
											DEFINE_BOUNDARY_PROPERTIES_TEMP(prop, prop_b, maxelm, maxbound-1);
										}

										doublereal dSloc = CALC_SQUARE(G, maxelm, nvtx, pa); // Вычисление площади грани ячейки.
										TOCHKA p_centerG = CALC_CENTERG(G, maxelm, nvtx, pa); // Вычисление геометрического центра грани.
										
										// G - грань по которой ставится граничное условие.
										// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
										// bound_id - номер граничного узла, начиная с нуля.
										// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
										if (octree1->b4S) {
											obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound - 1,
												octree1->linkS0->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, true, false, p_centerG);
										}
										else {
											if (octree1->linkS == NULL) {
												obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound - 1,
													 NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, true, false, p_centerG);
											}
											else {
												obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound - 1,
													octree1->linkS->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, true, false, p_centerG);
											}
										}
										binternalsource[maxbound - 1] = true; // внутренний источник тепла на границе жидкости и твёрдого тела

									}
									else {
										// Внутренний узел.
										G = NSIDE;
										// внутренняя грань:
										sosedi[G][maxelm].iNODE1 = octree1->linkN->inum_TD - 1; // граничные КО нумеруются в последнюю очередь,
										sosedi[G][maxelm].bdroblenie4 = false;
										sosedi[G][maxelm].iNODE2 = NON_EXISTENT_NODE;
										sosedi[G][maxelm].iNODE3 = NON_EXISTENT_NODE;
										sosedi[G][maxelm].iNODE4 = NON_EXISTENT_NODE;
										// после maxelm_memo внутренних КО.
										// Вычисление дальнего соседа.
										// CALC_GG(G, GG); // Вычисляет грань GG по грани G.
										// Определение дальнего соседа :
										// TODO 13 сентября 2016.

									}
									
								}
							}
						}
					}
					else {
					// Присутствуют все 4 соседа: N: 2,3,6,7.


						if (((octree1->linkN2 == NULL) && (octree1->linkN3 == NULL) && (octree1->linkN6 == NULL) && (octree1->linkN7 == NULL))||
							(((octree1->linkN2 != NULL) && (octree1->linkN3 != NULL) && (octree1->linkN6 != NULL) && (octree1->linkN7 != NULL) &&
							(octree1->linkN2->inum_TD == 0) && (octree1->linkN3->inum_TD == 0) && (octree1->linkN6->inum_TD == 0) && (octree1->linkN7->inum_TD == 0)))||
							(((octree1->linkN2 == NULL) && (octree1->linkN3 != NULL) && (octree1->linkN6 == NULL) && (octree1->linkN7 != NULL) &&
							(octree1->linkN3->inum_TD == 0) && (octree1->linkN7->inum_TD == 0)))||
							(((octree1->linkN2 != NULL) && (octree1->linkN3 != NULL) && (octree1->linkN6 == NULL) && (octree1->linkN7 == NULL) &&
							(octree1->linkN2->inum_TD == 0) && (octree1->linkN3->inum_TD == 0)))||
								(((octree1->linkN2 == NULL) && (octree1->linkN3 != NULL) && (octree1->linkN6 == NULL) && (octree1->linkN7 == NULL) &&
							(octree1->linkN3->inum_TD == 0))))						
						{
							// 1 - N2&&N3&&N6&&N7 - на границе расчётной области (РО).
							// 2 -  N2&&N3&&N6&&N7 - на гранеице HOLLOW блока.
							// 3 - N3&&N7 на границе HOLLOW блока, N2&&N6 на границе РО.
							// 4 - N2&&N3 на границе HOLLOW блока, N6&N7 на границе РО.
							// 5 - N3 на границе HOLLOW блока, N2&&N6&&N7  на границе РО.

							G = NSIDE;
							// граничная грань:
							INIT_SOSEDI(sosedi, maxelm, NSIDE, GG, false, maxelm_memo + maxbound,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE, false,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE); // заполнение sosedi[G][maxelm].
							DEFINE_BOUNDARY_PROPERTIES_TEMP(prop, prop_b, maxelm, maxbound); // В граничный узел сносятся свойства прилегающего КО.
							
							doublereal dSloc = CALC_SQUARE(G, maxelm, nvtx, pa); // Вычисление площади грани ячейки.
							TOCHKA p_centerG = CALC_CENTERG(G, maxelm, nvtx, pa); // Вычисление геометрического центра грани.
							
							// G - грань по которой ставится граничное условие.
							// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
							// bound_id - номер граничного узла, начиная с нуля.
							// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
							if (octree1->b4S) {
								obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
									octree1->linkS0->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, false, false, p_centerG);
							}
							else {
								if (octree1->linkS == NULL) {
									obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
										 NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, false, false, p_centerG);
								}
								else {
									obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
										octree1->linkS->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, false, false, p_centerG);
								}
							}

							// Зачем дробить если граница заканчивается NULL.
							maxbound++; // Вместо четырёх один.
						}
						else {

							bool binc1, binc2, binc3, binc4;
							integer lsid1, lsid2, lsid3, lsid4;

							if (octree1->linkN2 == NULL) {

								// 21 сентября 2016 (вырождение ячейки).
								
								// Заглушка.
								G = NSIDE;
								// граничная грань:
								sosedi[G][maxelm].iNODE1 = NON_EXISTENT_NODE; // граничные КО нумеруются в последнюю очередь,
								sosedi[G][maxelm].bdroblenie4 = true;
								// после maxelm_memo внутренних КО.
								// Вычисление дальнего соседа
								CALC_GG(G, GG); // Вычисляет грань GG по грани G.

								sosedi[GG][maxelm].iNODE1 = NON_EXISTENT_NODE; // нет соседа,
								sosedi[GG][maxelm].bdroblenie4 = true;


								/*
								G = NSIDE;
								// граничная грань:
								sosedi[G][maxelm].iNODE1 = maxelm_memo + maxbound; // граничные КО нумеруются в последнюю очередь,
								sosedi[G][maxelm].bdroblenie4 = true;
								// после maxelm_memo внутренних КО.
								// Вычисление дальнего соседа
								CALC_GG(G,GG); // Вычисляет грань GG по грани G.

								sosedi[GG][maxelm].iNODE1 = NON_EXISTENT_NODE; // нет соседа,
								sosedi[GG][maxelm].bdroblenie4 = true;
								#if doubleintprecision == 1
									printf("error NULL sosedi N2 elm=%lld\n", maxelm);
								#else
									printf("error NULL sosedi N2 elm=%d\n", maxelm);
								#endif
								
								getchar();

								// Эта ситуация возможна и здесь надо предусмотреть её обработку.
								DEFINE_BOUNDARY_PROPERTIES_TEMP(prop, prop_b, maxelm, maxbound); // В граничный узел сносятся свойства прилегающего КО.
								
								doublereal dS_loc = 0.0;

								// G - грань по которой ставится граничное условие.
								// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
								// bound_id - номер граничного узла, начиная с нуля.
								// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
								if (octree1->b4S) {
									obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
										octree1->linkS0->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
								}
								else {
									if (octree1->linkS == NULL) {
										obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
											 NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
									}
									else {
										obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
											octree1->linkS->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
									}
								}


								maxbound++;
								*/
							}
							else {
								if (octree1->linkN2->inum_TD == 0) {

									G = NSIDE;
									// граничная грань:
									sosedi[G][maxelm].iNODE1 = maxelm_memo + maxbound; // граничные КО нумеруются в последнюю очередь,
									sosedi[G][maxelm].bdroblenie4 = true;
									// после maxelm_memo внутренних КО.
									// Вычисление дальнего соседа
									CALC_GG(G, GG); // Вычисляет грань GG по грани G.

									sosedi[GG][maxelm].iNODE1 = NON_EXISTENT_NODE; // нет соседа,
									sosedi[GG][maxelm].bdroblenie4 = true;
#if doubleintprecision == 1
									//printf("error inum_TD == 0 sosedi N2 elm=%lld\n", maxelm);
#else
									//printf("error inum_TD == 0 sosedi N2 elm=%d\n", maxelm);
#endif
									
									//getchar();

									// Эта ситуация возможна и здесь надо предусмотреть её обработку.
									DEFINE_BOUNDARY_PROPERTIES_TEMP(prop, prop_b, maxelm, maxbound); // В граничный узел сносятся свойства прилегающего КО.
									
									// Sx*Sz
									doublereal dS_loc = fabs(octree1->linkN2->p4.z - octree1->linkN2->p0.z)*fabs(octree1->linkN2->p1.x - octree1->linkN2->p0.x);
									TOCHKA p_centerG;
									p_centerG.z = 0.5*(octree1->linkN2->p4.z + octree1->linkN2->p0.z);
									p_centerG.y = octree1->linkN2->p1.y;
									p_centerG.x = 0.5*(octree1->linkN2->p1.x + octree1->linkN2->p0.x);


									// G - грань по которой ставится граничное условие.
									// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
									// bound_id - номер граничного узла, начиная с нуля.
									// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
									if (octree1->b4S) {
										obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
											octree1->linkS0->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
									}
									else {
										if (octree1->linkS == NULL) {
											obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
												 NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
										}
										else {
											obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
												octree1->linkS->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
										}
									}

									maxbound++;
								}
								else {
									// сосед существует.
									if (bvisit[octree1->linkN2->inum_TD - 1]) {
										// узел уже был посещён
										// Здесь неким образом модифицировалось gran_t
										// Только внутренняя грань.
										G = NSIDE;
										// граничная грань:
										sosedi[G][maxelm].iNODE1 = octree1->linkN2->inum_TD - 1; // граничные КО нумеруются в последнюю очередь,
										sosedi[G][maxelm].bdroblenie4 = true;
										// после maxelm_memo внутренних КО.
										// Вычисление дальнего соседа
										CALC_GG(G, GG); // Вычисляет грань GG по грани G.

										// Определение дальнего соседа :
										// TODO 13 сентября 2016.
									}
									else {
										// узнать координаты центра грани и ориентацию в пространстве
										x_c = 0.5*(octree1->linkN2->p0.x + octree1->linkN2->p1.x);
										y_c = octree1->p3.y;
										z_c = 0.5*(octree1->linkN2->p0.z + octree1->linkN2->p4.z);
										iplane = XZ;
										patch_maxbound(iplane, s, ls, x_c, y_c, z_c, maxbound, binc1, lsid1);
										if (binc1) {
#if doubleintprecision == 1
											printf("error internal source N2 elm=%lld\n", maxelm);
#else
											printf("error internal source N2 elm=%d\n", maxelm);
#endif
											
											///getchar();
											system("PAUSE");
										}
										else {
											// Только внутренняя грань.
											G = NSIDE;
											// граничная грань:
											sosedi[G][maxelm].iNODE1 = octree1->linkN2->inum_TD - 1; // граничные КО нумеруются в последнюю очередь,
											sosedi[G][maxelm].bdroblenie4 = true;
											// после maxelm_memo внутренних КО.
											// Вычисление дальнего соседа
											CALC_GG(G, GG); // Вычисляет грань GG по грани G.

											// Определение дальнего соседа :
											// TODO 13 сентября 2016.
										}
									}
								}
							}

							if (octree1->linkN3 == NULL) {
#if doubleintprecision == 1
								printf("error NULL sosedi N3 elm=%lld\n", maxelm);
#else
								printf("error NULL sosedi N3 elm=%d\n", maxelm);
#endif
								
								//getchar();
								system("PAUSE");

								/*
								G = NSIDE;
								// граничная грань:
								sosedi[G][maxelm].iNODE2 = maxelm_memo + maxbound; // граничные КО нумеруются в последнюю очередь,
								sosedi[G][maxelm].bdroblenie4 = true;
								// после maxelm_memo внутренних КО.
								// Вычисление дальнего соседа
								CALC_GG(G,GG); // Вычисляет грань GG по грани G.

								sosedi[GG][maxelm].iNODE2 = NON_EXISTENT_NODE; // нет соседей,
								sosedi[GG][maxelm].bdroblenie4 = true;

								#if doubleintprecision == 1
										printf("error NULL sosedi N3 elm=%lld\n", maxelm);
								#else
										printf("error NULL sosedi N3 elm=%d\n", maxelm);
								#endif
								
								getchar();

								// Эта ситуация возможна и здесь надо предусмотреть её обработку.
								DEFINE_BOUNDARY_PROPERTIES_TEMP(prop, prop_b, maxelm, maxbound); // В граничный узел сносятся свойства прилегающего КО.
								


								doublereal dS_loc = 0.0;

								// G - грань по которой ставится граничное условие.
								// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
								// bound_id - номер граничного узла, начиная с нуля.
								// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
								if (octree1->b4S) {
									obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
										octree1->linkS0->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
								}
								else {
									if (octree1->linkS == NULL) {
										obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
											 NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
									}
									else {
										obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
											octree1->linkS->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
									}
								}


								maxbound++;
								*/

							}
							else {
								if (octree1->linkN3->inum_TD == 0) {

									G = NSIDE;
									// граничная грань:
									sosedi[G][maxelm].iNODE2 = maxelm_memo + maxbound; // граничные КО нумеруются в последнюю очередь,
									sosedi[G][maxelm].bdroblenie4 = true;
									// после maxelm_memo внутренних КО.
									// Вычисление дальнего соседа
									CALC_GG(G, GG); // Вычисляет грань GG по грани G.

									sosedi[GG][maxelm].iNODE2 = NON_EXISTENT_NODE; // нет соседей,
									sosedi[GG][maxelm].bdroblenie4 = true;

#if doubleintprecision == 1
									//printf("error inum_TD == 0 sosedi N3 elm=%lld\n", maxelm);
#else
									//printf("error inum_TD == 0 sosedi N3 elm=%d\n", maxelm);
#endif
									
									//getchar();

									// Эта ситуация возможна и здесь надо предусмотреть её обработку.
									DEFINE_BOUNDARY_PROPERTIES_TEMP(prop, prop_b, maxelm, maxbound); // В граничный узел сносятся свойства прилегающего КО.
									
									// Sx*Sz
									doublereal dS_loc = fabs(octree1->linkN3->p4.z - octree1->linkN3->p0.z)*fabs(octree1->linkN3->p1.x - octree1->linkN3->p0.x);
									TOCHKA p_centerG;
									p_centerG.z = 0.5*(octree1->linkN3->p4.z + octree1->linkN3->p0.z);
									p_centerG.y = octree1->linkN3->p1.y;
									p_centerG.x = 0.5*(octree1->linkN3->p1.x + octree1->linkN3->p0.x);

									// G - грань по которой ставится граничное условие.
									// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
									// bound_id - номер граничного узла, начиная с нуля.
									// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
									if (octree1->b4S) {
										obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
											octree1->linkS0->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
									}
									else {
										if (octree1->linkS == NULL) {
											obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
												 NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
										}
										else {
											obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
												octree1->linkS->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
										}
									}


									maxbound++;
								}
								else {
									// сосед существует.
									if (bvisit[octree1->linkN3->inum_TD - 1]) {
										// узел уже был посещён
										// Здесь неким образом модифицировалось gran_t
										// Только внутренняя грань.
										G = NSIDE;
										// граничная грань:
										sosedi[G][maxelm].iNODE2 = octree1->linkN3->inum_TD - 1; // граничные КО нумеруются в последнюю очередь,
										sosedi[G][maxelm].bdroblenie4 = true;
										// после maxelm_memo внутренних КО.
										// Вычисление дальнего соседа
										CALC_GG(G, GG); // Вычисляет грань GG по грани G.

										// Определение дальнего соседа :
										// TODO 13 сентября 2016.
									}
									else {
										// узнать координаты центра грани и ориентацию в пространстве
										x_c = 0.5*(octree1->linkN3->p0.x + octree1->linkN3->p1.x);
										y_c = octree1->p3.y;
										z_c = 0.5*(octree1->linkN3->p0.z + octree1->linkN3->p4.z);
										iplane = XZ;
										patch_maxbound(iplane, s, ls, x_c, y_c, z_c, maxbound, binc2, lsid2);
										if (binc2) {
#if doubleintprecision == 1
											printf("error internal source N3 elm=%lld\n", maxelm);
#else
											printf("error internal source N3 elm=%d\n", maxelm);
#endif
											
											//getchar();
											system("PAUSE");
										}
										else {
											// Только внутренняя грань.
											G = NSIDE;
											// граничная грань:
											sosedi[G][maxelm].iNODE2 = octree1->linkN3->inum_TD - 1; // граничные КО нумеруются в последнюю очередь,
											sosedi[G][maxelm].bdroblenie4 = true;
											// после maxelm_memo внутренних КО.
											// Вычисление дальнего соседа
											CALC_GG(G, GG); // Вычисляет грань GG по грани G.

											// Определение дальнего соседа :
											// TODO 13 сентября 2016.
										}
									}
								}
							}

							if (octree1->linkN6 == NULL) {
								
								// 21 сентября 2016 (вырождение ячейки).
								
								// Заглушка.
								G = NSIDE;
								// граничная грань:
								sosedi[G][maxelm].iNODE3 = NON_EXISTENT_NODE; // граничные КО нумеруются в последнюю очередь,
								sosedi[G][maxelm].bdroblenie4 = true;
								// после maxelm_memo внутренних КО.
								// Вычисление дальнего соседа
								CALC_GG(G, GG); // Вычисляет грань GG по грани G.

								sosedi[GG][maxelm].iNODE3 = NON_EXISTENT_NODE; // сосед отсутствует,
								sosedi[GG][maxelm].bdroblenie4 = true;

								/*
								G = NSIDE;
								// граничная грань:
								sosedi[G][maxelm].iNODE3 = maxelm_memo + maxbound; // граничные КО нумеруются в последнюю очередь,
								sosedi[G][maxelm].bdroblenie4 = true;
								// после maxelm_memo внутренних КО.
								// Вычисление дальнего соседа
								CALC_GG(G,GG); // Вычисляет грань GG по грани G.
								
								sosedi[GG][maxelm].iNODE3 = NON_EXISTENT_NODE; // сосед отсутствует,
								sosedi[GG][maxelm].bdroblenie4 = true;

								#if doubleintprecision == 1
										printf("error NULL sosedi N6 elm=%lld\n", maxelm);
								#else
										printf("error NULL sosedi N6 elm=%d\n", maxelm);
								#endif
								
								getchar();

								// Эта ситуация возможна и здесь надо предусмотреть её обработку.
								DEFINE_BOUNDARY_PROPERTIES_TEMP(prop, prop_b, maxelm, maxbound); // В граничный узел сносятся свойства прилегающего КО.
								
								doublereal dS_loc = 0.0;

								// G - грань по которой ставится граничное условие.
								// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
								// bound_id - номер граничного узла, начиная с нуля.
								// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
								if (octree1->b4S) {
									obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
										octree1->linkS0->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
								}
								else {
									if (octree1->linkS == NULL) {
										obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
											 NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
									}
									else {
										obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
											octree1->linkS->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
									}
								}

								maxbound++;
								*/

							}
							else {
								if (octree1->linkN6->inum_TD == 0) {

									G = NSIDE;
									// граничная грань:
									sosedi[G][maxelm].iNODE3 = maxelm_memo + maxbound; // граничные КО нумеруются в последнюю очередь,
									sosedi[G][maxelm].bdroblenie4 = true;
									// после maxelm_memo внутренних КО.
									// Вычисление дальнего соседа
									CALC_GG(G, GG); // Вычисляет грань GG по грани G.

									sosedi[GG][maxelm].iNODE3 = NON_EXISTENT_NODE; // сосед отсутствует,
									sosedi[GG][maxelm].bdroblenie4 = true;

#if doubleintprecision == 1
									//printf("error inum_TD == 0 sosedi N6 elm=%lld\n", maxelm);
#else
									//printf("error inum_TD == 0 sosedi N6 elm=%d\n", maxelm);
#endif
									
									//getchar();

									// Эта ситуация возможна и здесь надо предусмотреть её обработку.
									DEFINE_BOUNDARY_PROPERTIES_TEMP(prop, prop_b, maxelm, maxbound); // В граничный узел сносятся свойства прилегающего КО.
									
									// Sx*Sz
									doublereal dS_loc = fabs(octree1->linkN6->p4.z - octree1->linkN6->p0.z)*fabs(octree1->linkN6->p1.x - octree1->linkN6->p0.x);
									TOCHKA p_centerG;
									p_centerG.z = 0.5*(octree1->linkN6->p4.z + octree1->linkN6->p0.z);
									p_centerG.y = octree1->linkN6->p1.y;
									p_centerG.x = 0.5*(octree1->linkN6->p1.x + octree1->linkN6->p0.x);


									// G - грань по которой ставится граничное условие.
									// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
									// bound_id - номер граничного узла, начиная с нуля.
									// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
									if (octree1->b4S) {
										obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
											octree1->linkS0->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
									}
									else {
										if (octree1->linkS == NULL) {
											obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
												 NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
										}
										else {
											obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
												octree1->linkS->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
										}
									}

									maxbound++;
								}
								else {
									// сосед существует.
									if (bvisit[octree1->linkN6->inum_TD - 1]) {
										// узел уже был посещён
										// Здесь неким образом модифицировалось gran_t
										// Только внутренняя грань.
										G = NSIDE;
										// граничная грань:
										sosedi[G][maxelm].iNODE3 = octree1->linkN6->inum_TD - 1; // граничные КО нумеруются в последнюю очередь,
										sosedi[G][maxelm].bdroblenie4 = true;
										// после maxelm_memo внутренних КО.
										// Вычисление дальнего соседа
										CALC_GG(G, GG); // Вычисляет грань GG по грани G.

										// Определение дальнего соседа :
										// TODO 13 сентября 2016.
									}
									else {
										// узнать координаты центра грани и ориентацию в пространстве
										x_c = 0.5*(octree1->linkN6->p0.x + octree1->linkN6->p1.x);
										y_c = octree1->p3.y;
										z_c = 0.5*(octree1->linkN6->p0.z + octree1->linkN6->p4.z);
										iplane = XZ;
										patch_maxbound(iplane, s, ls, x_c, y_c, z_c, maxbound, binc3, lsid3);
										if (binc3) {
#if doubleintprecision == 1
											printf("error internal source N6 elm=%lld\n", maxelm);
#else
											printf("error internal source N6 elm=%d\n", maxelm);
#endif
											
											//getchar();
											system("PAUSE");
										}
										else {
											// Только внутренняя грань.
											G = NSIDE;
											// граничная грань:
											sosedi[G][maxelm].iNODE3 = octree1->linkN6->inum_TD - 1; // граничные КО нумеруются в последнюю очередь,
											sosedi[G][maxelm].bdroblenie4 = true;
											// после maxelm_memo внутренних КО.
											// Вычисление дальнего соседа
											CALC_GG(G, GG); // Вычисляет грань GG по грани G.

											// Определение дальнего соседа :
											// TODO 13 сентября 2016.
										}
									}
								}
							}

							if (octree1->linkN7 == NULL) {

								// 21 сентября 2016 (вырождение ячейки).

								// Заглушка.
								G = NSIDE;
								// граничная грань:
								sosedi[G][maxelm].iNODE4 = NON_EXISTENT_NODE; // граничные КО нумеруются в последнюю очередь,
								sosedi[G][maxelm].bdroblenie4 = true;
								// после maxelm_memo внутренних КО.
								// Вычисление дальнего соседа
								CALC_GG(G, GG); // Вычисляет грань GG по грани G.

								sosedi[GG][maxelm].iNODE4 = NON_EXISTENT_NODE; // граничные КО нумеруются в последнюю очередь,
								sosedi[GG][maxelm].bdroblenie4 = true;


								/*
								G = NSIDE;
								// граничная грань:
								sosedi[G][maxelm].iNODE4 = maxelm_memo + maxbound; // граничные КО нумеруются в последнюю очередь,
								sosedi[G][maxelm].bdroblenie4 = true;
								// после maxelm_memo внутренних КО.
								// Вычисление дальнего соседа
								CALC_GG(G,GG); // Вычисляет грань GG по грани G.

								sosedi[GG][maxelm].iNODE4 = - 1; // граничные КО нумеруются в последнюю очередь,
								sosedi[GG][maxelm].bdroblenie4 = true;
								#if doubleintprecision == 1
										printf("error NULL sosedi N7 elm=%lld\n", maxelm);
#else
		printf("error NULL sosedi N7 elm=%d\n", maxelm);
#endif
								
								getchar();

								// Эта ситуация возможна и здесь надо предусмотреть её обработку.
								DEFINE_BOUNDARY_PROPERTIES_TEMP(prop, prop_b, maxelm, maxbound); // В граничный узел сносятся свойства прилегающего КО.
								
								doublereal dS_loc = 0.0;

								// G - грань по которой ставится граничное условие.
								// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
								// bound_id - номер граничного узла, начиная с нуля.
								// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
								if (octree1->b4S) {
									obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
										octree1->linkS0->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
								}
								else {
									if (octree1->linkS == NULL) {
										obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
											 NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
									}
									else {
										obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
											octree1->linkS->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
									}
								}


								maxbound++;
								*/
							}
							else {
								if (octree1->linkN7->inum_TD == 0) {

									G = NSIDE;
									// граничная грань:
									sosedi[G][maxelm].iNODE4 = maxelm_memo + maxbound; // граничные КО нумеруются в последнюю очередь,
									sosedi[G][maxelm].bdroblenie4 = true;
									// после maxelm_memo внутренних КО.
									// Вычисление дальнего соседа
									CALC_GG(G, GG); // Вычисляет грань GG по грани G.

									sosedi[GG][maxelm].iNODE4 = NON_EXISTENT_NODE; // граничные КО нумеруются в последнюю очередь,
									sosedi[GG][maxelm].bdroblenie4 = true;

#if doubleintprecision == 1
									//printf("error inum_TD == 0 sosedi N7 elm=%lld\n", maxelm);
#else
									//printf("error inum_TD == 0 sosedi N7 elm=%d\n", maxelm);
#endif
									
									//getchar();

									// Эта ситуация возможна и здесь надо предусмотреть её обработку.
									DEFINE_BOUNDARY_PROPERTIES_TEMP(prop, prop_b, maxelm, maxbound); // В граничный узел сносятся свойства прилегающего КО.
									
									// Sx*Sz
									doublereal dS_loc = fabs(octree1->linkN7->p4.z - octree1->linkN7->p0.z)*fabs(octree1->linkN7->p1.x - octree1->linkN7->p0.x);
									TOCHKA p_centerG;
									p_centerG.z = 0.5*(octree1->linkN7->p4.z + octree1->linkN7->p0.z);
									p_centerG.y = octree1->linkN7->p1.y;
									p_centerG.x = 0.5*(octree1->linkN7->p1.x + octree1->linkN7->p0.x);


									// G - грань по которой ставится граничное условие.
									// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
									// bound_id - номер граничного узла, начиная с нуля.
									// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
									if (octree1->b4S) {
										obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
											octree1->linkS0->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
									}
									else {
										if (octree1->linkS == NULL) {
											obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
												 NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
										}
										else {
											obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
												octree1->linkS->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
										}
									}

									maxbound++;
								}
								else {
									// сосед существует.
									if (bvisit[octree1->linkN7->inum_TD - 1]) {
										// узел уже был посещён
										// Здесь неким образом модифицировалось gran_t
										// Только внутренняя грань.
										G = NSIDE;
										// граничная грань:
										sosedi[G][maxelm].iNODE4 = octree1->linkN7->inum_TD - 1; // граничные КО нумеруются в последнюю очередь,
										sosedi[G][maxelm].bdroblenie4 = true;
										// после maxelm_memo внутренних КО.
										// Вычисление дальнего соседа
										CALC_GG(G, GG); // Вычисляет грань GG по грани G.

										// Определение дальнего соседа :
										// TODO 13 сентября 2016.
									}
									else {
										// узнать координаты центра грани и ориентацию в пространстве
										x_c = 0.5*(octree1->linkN7->p0.x + octree1->linkN7->p1.x);
										y_c = octree1->p3.y;
										z_c = 0.5*(octree1->linkN7->p0.z + octree1->linkN7->p4.z);
										iplane = XZ;
										patch_maxbound(iplane, s, ls, x_c, y_c, z_c, maxbound, binc4, lsid4);
										if (binc4) {
#if doubleintprecision == 1
											printf("error internal source N7 elm=%lld\n", maxelm);
#else
											printf("error internal source N7 elm=%d\n", maxelm);
#endif
											
											//getchar();
											system("PAUSE");
										}
										else {
											// Только внутренняя грань.
											G = NSIDE;
											// граничная грань:
											sosedi[G][maxelm].iNODE4 = octree1->linkN7->inum_TD - 1; // граничные КО нумеруются в последнюю очередь,
											sosedi[G][maxelm].bdroblenie4 = true;
											// после maxelm_memo внутренних КО.
											// Вычисление дальнего соседа
											CALC_GG(G, GG); // Вычисляет грань GG по грани G.

											// Определение дальнего соседа :
											// TODO 13 сентября 2016.
										}
									}
								}
							}
						}
					}

					// сосед не существует
					// записываем идентификатор грани и
					// увеличиваем счётчик граней.
					if (!octree1->b4S) {
						if ((octree1->linkS == NULL)||((octree1->linkS != NULL)&&(octree1->linkS->inum_TD == 0))) {
							// Если мы находимся на границе расчётной области или если мы граничим с HOLLOW блоком.

							G = SSIDE;
							// граничная грань:
							INIT_SOSEDI(sosedi, maxelm, SSIDE, GG, false, maxelm_memo + maxbound,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE, false,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE); // заполнение sosedi[G][maxelm].
							DEFINE_BOUNDARY_PROPERTIES_TEMP(prop, prop_b, maxelm, maxbound); // В граничный узел сносятся свойства прилегающего КО.
							
							doublereal dSloc = CALC_SQUARE(G, maxelm, nvtx, pa); // Вычисление площади грани ячейки.
							TOCHKA p_centerG = CALC_CENTERG(G, maxelm, nvtx, pa); // Вычисление геометрического центра грани.
							
							// G - грань по которой ставится граничное условие.
							// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
							// bound_id - номер граничного узла, начиная с нуля.
							// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
							if (octree1->b4N) {
								obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
									octree1->linkN3->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, false, false, p_centerG);
							}
							else {
								if (octree1->linkN == NULL) {
									obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
										 NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, false, false, p_centerG);
								}
								else {
									obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
										octree1->linkN->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, false, false, p_centerG);
								}
							}

							maxbound++;
						}
						else {
							if (octree1->linkS != NULL) {
								//bool binc;
								//integer lsid;
								// сосед существует.
								if (bvisit[octree1->linkS->inum_TD - 1]) {


									// return_current_node_number - Заглянули в соседа и вернулись из соседа в текущий узел.
									// Номер должен быть тем же самым и поэтому этот код бессмысленен.
									// Строго внутренние плоские бесконечно тонкие источники тепла внутри расчётной области больше не поддерживаются.
									// Поэтому от этого кода можно безболезненно избавиться. Поддерживаются только объёмные источники тепла.
									// 10.03.2019
									integer return_current_node_number = sosedi[NSIDE][octree1->linkS->inum_TD - 1].iNODE1;
									//N2 NODE1,
									//N3 NODE2,
									//N6 NODE3,
									//N7 NODE4.
									// TODO 23 сентября 2016.
									if (octree1->linkS->b4N) {
										integer current_node_number = octree1->inum_TD - 1;
										bool bcontinue = true;
										if (bcontinue) {
											if (octree1->linkS->linkN2!=NULL) {
												if (octree1->linkS->linkN2->inum_TD - 1 == current_node_number) {
													bcontinue = false; // iNODE1
													if (sosedi[NSIDE][octree1->linkS->inum_TD - 1].iNODE1 == current_node_number) {
														return_current_node_number = sosedi[NSIDE][octree1->linkS->inum_TD - 1].iNODE1;
													}
													else {
														printf("this can not be S N2\n");
														//getchar();
														system("PAUSE");
													}
												}
											}
										}
										if (bcontinue) {
											if (octree1->linkS->linkN3!=NULL) {
												if (octree1->linkS->linkN3->inum_TD - 1 == current_node_number) {
													bcontinue = false; // iNODE2
													if (sosedi[NSIDE][octree1->linkS->inum_TD - 1].iNODE2 == current_node_number) {
														return_current_node_number = sosedi[NSIDE][octree1->linkS->inum_TD - 1].iNODE2;
													}
													else {
														printf("this can not be S N3\n");
														//getchar();
														system("PAUSE");
													}
												}
											}
										}
										if (bcontinue) {
											if (octree1->linkS->linkN6!=NULL) {
												if (octree1->linkS->linkN6->inum_TD - 1 == current_node_number) {
													bcontinue = false; // iNODE3
													if (sosedi[NSIDE][octree1->linkS->inum_TD - 1].iNODE3 == current_node_number) {
														return_current_node_number = sosedi[NSIDE][octree1->linkS->inum_TD - 1].iNODE3;
													}
													else {
														printf("this can not be S N6\n");
														//getchar();
														system("PAUSE");
													}
												}
											}
										}
										if (bcontinue) {
											if (octree1->linkS->linkN7 != NULL) {
												if (octree1->linkS->linkN7->inum_TD - 1 == current_node_number) {
													bcontinue = false; // iNODE4
													if (sosedi[NSIDE][octree1->linkS->inum_TD - 1].iNODE4 == current_node_number) {
														return_current_node_number = sosedi[NSIDE][octree1->linkS->inum_TD - 1].iNODE4;
													}
													else {
														printf("this can not be S N7\n");
														//getchar();
														system("PAUSE");
													}
												}
											}
										}

									}


									// узел уже был посещён
									// Здесь неким образом модифицировалось gran_t
									if (return_current_node_number >= maxelm_memo) {

										printf("WARNING!!! internal source gran SSIDE 10.03.2019\n");
										system("PAUSE");

										// Это граничная грань внутреннего источника тепла.
										// Внутренний источник тепла.
										G = SSIDE;
										// граничная грань:
										INIT_SOSEDI(sosedi, maxelm, SSIDE, GG, false, sosedi[NSIDE][octree1->linkS->inum_TD - 1].iNODE1,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE, false,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE); // заполнение sosedi[G][maxelm].

										if (bi_fluid) {
											// FLUID
											 // В граничный узел сносятся свойства прилегающего КО.
											DEFINE_BOUNDARY_PROPERTIES_TEMP(prop, prop_b, octree1->linkS->inum_TD - 1, sosedi[G][maxelm].iNODE1 - maxelm_memo);
										}
										else {
											// SOLID
											 // В граничный узел сносятся свойства прилегающего КО.
											DEFINE_BOUNDARY_PROPERTIES_TEMP(prop, prop_b, maxelm, sosedi[G][maxelm].iNODE1 - maxelm_memo);
										}

										doublereal dSloc = CALC_SQUARE(G, maxelm, nvtx, pa); // Вычисление площади грани ячейки.
										TOCHKA p_centerG = CALC_CENTERG(G, maxelm, nvtx, pa); // Вычисление геометрического центра грани.
										
										// G - грань по которой ставится граничное условие.
										// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
										// bound_id - номер граничного узла, начиная с нуля.
										// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
										if (octree1->b4N) {
											obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, sosedi[NSIDE][octree1->linkS->inum_TD - 1].iNODE1 - maxelm_memo,
												octree1->linkN3->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, true, true, p_centerG);
										}
										else {
											if (octree1->linkN == NULL) {
												obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, sosedi[NSIDE][octree1->linkS->inum_TD - 1].iNODE1 - maxelm_memo,
													 NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, true, true, p_centerG);
											}
											else {
												obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, sosedi[NSIDE][octree1->linkS->inum_TD - 1].iNODE1 - maxelm_memo,
													octree1->linkN->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, true, true, p_centerG);
											}
										}
										binternalsource[sosedi[NSIDE][octree1->linkS->inum_TD - 1].iNODE1 - maxelm_memo] = true; // внутренний источник тепла на границе жидкости и твёрдого тела

									}
									else {
										// Это строго внутренний узел.
										// Внутренний узел.
										G = SSIDE;
										// внутренняя грань:
										sosedi[G][maxelm].iNODE1 = octree1->linkS->inum_TD - 1; // граничные КО нумеруются в последнюю очередь,
										sosedi[G][maxelm].bdroblenie4 = false;
										sosedi[G][maxelm].iNODE2 = NON_EXISTENT_NODE;
										sosedi[G][maxelm].iNODE3 = NON_EXISTENT_NODE;
										sosedi[G][maxelm].iNODE4 = NON_EXISTENT_NODE;
										// после maxelm_memo внутренних КО.
										// Вычисление дальнего соседа
										//CALC_GG(G, GG); // Вычисляет грань GG по грани G.

										// Определение дальнего соседа :
										// TODO 13 сентября 2016.
									}

								}
								else {
									bool binc;
									integer lsid;
									// узнать координаты центра грани и ориентацию в пространстве
									x_c = 0.5*(octree1->p0.x + octree1->p1.x);
									y_c = octree1->p0.y;
									z_c = 0.5*(octree1->p0.z + octree1->p4.z);
									iplane = XZ;
									patch_maxbound(iplane, s, ls, x_c, y_c, z_c, maxbound,binc,lsid);
									if (binc) {
										// Внутренний источник тепла.
										// Внутренние плоские бесконечно тонкие источники тепла скорее всего больше не поддерживаются. 10,03,2019

										G = SSIDE;
										// граничная грань:
										INIT_SOSEDI(sosedi, maxelm, SSIDE, GG, false, maxelm_memo + maxbound-1,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE, false,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE); // заполнение sosedi[G][maxelm].
										
										if (bi_fluid) {
											// FLUID
											 // В граничный узел сносятся свойства прилегающего КО.
											DEFINE_BOUNDARY_PROPERTIES_TEMP(prop, prop_b, octree1->linkS->inum_TD - 1, maxbound-1);
										}
										else {
											// SOLID
											 // В граничный узел сносятся свойства прилегающего КО.
											DEFINE_BOUNDARY_PROPERTIES_TEMP(prop, prop_b, maxelm, maxbound-1);
										}

										doublereal dSloc = CALC_SQUARE(G, maxelm, nvtx, pa); // Вычисление площади грани ячейки.
										TOCHKA p_centerG = CALC_CENTERG(G, maxelm, nvtx, pa); // Вычисление геометрического центра грани.
										
										// G - грань по которой ставится граничное условие.
										// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
										// bound_id - номер граничного узла, начиная с нуля.
										// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
										if (octree1->b4N) {
											obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound-1,
												octree1->linkN3->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, true, false, p_centerG);
										}
										else {
											if (octree1->linkN == NULL) {
												obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound-1,
													 NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, true, false, p_centerG);
											}
											else {
												obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound-1,
													octree1->linkN->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, true, false, p_centerG);
											}
										}
										binternalsource[maxbound - 1] = true; // внутренний источник тепла на границе жидкости и твёрдого тела

									}
									else {
										// Внутренний узел.
										G = SSIDE;
										// внутренняя грань:
										sosedi[G][maxelm].iNODE1 = octree1->linkS->inum_TD - 1; // граничные КО нумеруются в последнюю очередь,
										sosedi[G][maxelm].bdroblenie4 = false;
										sosedi[G][maxelm].iNODE2 = NON_EXISTENT_NODE;
										sosedi[G][maxelm].iNODE3 = NON_EXISTENT_NODE;
										sosedi[G][maxelm].iNODE4 = NON_EXISTENT_NODE;
										// после maxelm_memo внутренних КО.
										// Вычисление дальнего соседа
										//CALC_GG(G, GG); // Вычисляет грань GG по грани G.

										// Определение дальнего соседа :
										// TODO 13 сентября 2016.
									}
								}
							}
						}
					}
					else {
					// Присутствуют все 4 соседа: S: 0,1,4,5.


						if (((octree1->linkS0 == NULL) && (octree1->linkS1 == NULL) && (octree1->linkS4 == NULL) && (octree1->linkS5 == NULL)) ||
							(((octree1->linkS0 != NULL) && (octree1->linkS1 != NULL) && (octree1->linkS4 != NULL) && (octree1->linkS5 != NULL) &&
							(octree1->linkS0->inum_TD == 0) && (octree1->linkS1->inum_TD == 0) && (octree1->linkS4->inum_TD == 0) && (octree1->linkS5->inum_TD == 0)))||
							(((octree1->linkS0 != NULL) && (octree1->linkS1 == NULL) && (octree1->linkS4 != NULL) && (octree1->linkS5 == NULL) &&
							(octree1->linkS0->inum_TD == 0) && (octree1->linkS4->inum_TD == 0)))||
							(((octree1->linkS0 != NULL) && (octree1->linkS1 != NULL) && (octree1->linkS4 == NULL) && (octree1->linkS5 == NULL) &&
							(octree1->linkS0->inum_TD == 0) && (octree1->linkS1->inum_TD == 0)))||
							(((octree1->linkS0 != NULL) && (octree1->linkS1 == NULL) && (octree1->linkS4 == NULL) && (octree1->linkS5 == NULL) &&
							(octree1->linkS0->inum_TD == 0))))
						{
							// 1 - S0&&S1&&S4&&S5 - на границе расчётной области (РО).
							// 2 - S0&&S1&&S4&&S5 - на гранеице HOLLOW блока.
							// 3 - S0&&S4 на границе HOLLOW блока, S1&&S5 на границе РО.
							// 4 - S0&&S1 на границе HOLLOW блока, S4&S5 на границе РО.
							// 5 - S0 на границе HOLLOW блока, S1&&S4&&S5  на границе РО.

							G = SSIDE;
							// граничная грань:
							INIT_SOSEDI(sosedi, maxelm, SSIDE, GG, false, maxelm_memo + maxbound,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE, false,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE); // заполнение sosedi[G][maxelm].
							DEFINE_BOUNDARY_PROPERTIES_TEMP(prop, prop_b, maxelm, maxbound); // В граничный узел сносятся свойства прилегающего КО.
								
							doublereal dSloc = CALC_SQUARE(G, maxelm, nvtx, pa); // Вычисление площади грани ячейки.
							TOCHKA p_centerG = CALC_CENTERG(G, maxelm, nvtx, pa); // Вычисление геометрического центра грани.
							
							// G - грань по которой ставится граничное условие.
							// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
							// bound_id - номер граничного узла, начиная с нуля.
							// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
							if (octree1->b4N) {
								obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
									octree1->linkN3->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, false, false, p_centerG);
							}
							else {
								if (octree1->linkN == NULL) {
									obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
										 NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, false, false, p_centerG);
								}
								else {
									obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
										octree1->linkN->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, false, false, p_centerG);
								}
							}

							// Зачем дробить если граница заканчивается NULL.
							maxbound++; // Вместо четырёх один.
						}					
						else {

							bool binc1, binc2, binc3, binc4;
							integer lsid1, lsid2, lsid3, lsid4;

							if (octree1->linkS0 == NULL) {
#if doubleintprecision == 1
								printf("error NULL sosedi S0 elm=%lld\n", maxelm);
#else
								printf("error NULL sosedi S0 elm=%d\n", maxelm);
#endif
								
								//getchar();
								system("PAUSE");

								/*

								G = SSIDE;
								// граничная грань:
								sosedi[G][maxelm].iNODE1 = maxelm_memo + maxbound; // граничные КО нумеруются в последнюю очередь,
								sosedi[G][maxelm].bdroblenie4 = true;
								// после maxelm_memo внутренних КО.
								// Вычисление дальнего соседа
								CALC_GG(G,GG); // Вычисляет грань GG по грани G.

								sosedi[GG][maxelm].iNODE1 = NON_EXISTENT_NODE; // соседа нет,
								sosedi[GG][maxelm].bdroblenie4 = true;
								#if doubleintprecision == 1
									printf("error NULL sosedi S0 elm=%lld\n", maxelm);
								#else
									printf("error NULL sosedi S0 elm=%d\n", maxelm);
								#endif
								
								getchar();

								// Эта ситуация возможна и здесь надо предусмотреть её обработку.
								DEFINE_BOUNDARY_PROPERTIES_TEMP(prop, prop_b, maxelm, maxbound); // В граничный узел сносятся свойства прилегающего КО.
								
								doublereal dS_loc = 0.0;

								// G - грань по которой ставится граничное условие.
								// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
								// bound_id - номер граничного узла, начиная с нуля.
								// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
								if (octree1->b4N) {
									obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
										octree1->linkN3->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
								}
								else {
									if (octree1->linkN == NULL) {
										obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
											 NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
									}
									else {
										obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
											octree1->linkN->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
									}
								}


								maxbound++;
								*/
							}
							else {
								if (octree1->linkS0->inum_TD == 0) {
									
									G = SSIDE;
									// граничная грань:
									sosedi[G][maxelm].iNODE1 = maxelm_memo + maxbound; // граничные КО нумеруются в последнюю очередь,
									sosedi[G][maxelm].bdroblenie4 = true;
									// после maxelm_memo внутренних КО.
									// Вычисление дальнего соседа
									CALC_GG(G, GG); // Вычисляет грань GG по грани G.

									sosedi[GG][maxelm].iNODE1 = NON_EXISTENT_NODE; // соседа нет,
									sosedi[GG][maxelm].bdroblenie4 = true;
									
#if doubleintprecision == 1
									//printf("error inum_TD == 0 sosedi S0 elm=%lld\n", maxelm);
#else
									//printf("error inum_TD == 0 sosedi S0 elm=%d\n", maxelm);
#endif
									
									//getchar();

									// Эта ситуация возможна и здесь надо предусмотреть её обработку.
									DEFINE_BOUNDARY_PROPERTIES_TEMP(prop, prop_b, maxelm, maxbound); // В граничный узел сносятся свойства прилегающего КО.
									
									// Sx*Sz
									doublereal dS_loc = fabs(octree1->linkS0->p7.z - octree1->linkS0->p3.z)*fabs(octree1->linkS0->p2.x - octree1->linkS0->p3.x);
									TOCHKA p_centerG;
									p_centerG.z = 0.5*(octree1->linkS0->p7.z + octree1->linkS0->p3.z);
									p_centerG.y = octree1->linkS0->p2.y;
									p_centerG.x = 0.5*(octree1->linkS0->p2.x + octree1->linkS0->p3.x);


									// G - грань по которой ставится граничное условие.
									// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
									// bound_id - номер граничного узла, начиная с нуля.
									// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
									if (octree1->b4N) {
										obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
											octree1->linkN3->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
									}
									else {
										if (octree1->linkN == NULL) {
											obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
												 NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
										}
										else {
											obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
												octree1->linkN->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
										}
									}

									maxbound++;
								}
								else {
									// сосед существует.
									if (bvisit[octree1->linkS0->inum_TD - 1]) {
										// узел уже был посещён
										// Здесь неким образом модифицировалось gran_t
										// Только внутренняя грань.
										G = SSIDE;
										// граничная грань:
										sosedi[G][maxelm].iNODE1 = octree1->linkS0->inum_TD - 1; // граничные КО нумеруются в последнюю очередь,
										sosedi[G][maxelm].bdroblenie4 = true;
										// после maxelm_memo внутренних КО.
										// Вычисление дальнего соседа
										CALC_GG(G, GG); // Вычисляет грань GG по грани G.

										// Определение дальнего соседа :
										// TODO 13 сентября 2016.
									}
									else {
										// узнать координаты центра грани и ориентацию в пространстве
										x_c = 0.5*(octree1->linkS0->p3.x + octree1->linkS0->p2.x);
										y_c = octree1->p0.y;
										z_c = 0.5*(octree1->linkS0->p3.z + octree1->linkS0->p7.z);
										iplane = XZ;
										patch_maxbound(iplane, s, ls, x_c, y_c, z_c, maxbound, binc1, lsid1);
										if (binc1) {
#if doubleintprecision == 1
											printf("error internal source S0 elm=%lld\n", maxelm);
#else
											printf("error internal source S0 elm=%d\n", maxelm);
#endif
											
											//getchar();
											system("PAUSE");
										}
										else {
											// Только внутренняя грань.
											G = SSIDE;
											// граничная грань:
											sosedi[G][maxelm].iNODE1 = octree1->linkS0->inum_TD - 1; // граничные КО нумеруются в последнюю очередь,
											sosedi[G][maxelm].bdroblenie4 = true;
											// после maxelm_memo внутренних КО.
											// Вычисление дальнего соседа
											CALC_GG(G, GG); // Вычисляет грань GG по грани G.

											// Определение дальнего соседа :
											// TODO 13 сентября 2016.
										}
									}
								}
							}

							if (octree1->linkS1 == NULL) {

								// 21 сентября 2016 (вырождение ячейки).

								// Заглушка.
								G = SSIDE;
								// граничная грань:
								sosedi[G][maxelm].iNODE2 = NON_EXISTENT_NODE; // граничные КО нумеруются в последнюю очередь,
								sosedi[G][maxelm].bdroblenie4 = true;
								// после maxelm_memo внутренних КО.
								// Вычисление дальнего соседа
								CALC_GG(G, GG); // Вычисляет грань GG по грани G.

								sosedi[GG][maxelm].iNODE2 = NON_EXISTENT_NODE; // сосед отсутствует,
								sosedi[GG][maxelm].bdroblenie4 = true;


								/*
								G = SSIDE;
								// граничная грань:
								sosedi[G][maxelm].iNODE2 = maxelm_memo + maxbound; // граничные КО нумеруются в последнюю очередь,
								sosedi[G][maxelm].bdroblenie4 = true;
								// после maxelm_memo внутренних КО.
								// Вычисление дальнего соседа
								CALC_GG(G,GG); // Вычисляет грань GG по грани G.

								sosedi[GG][maxelm].iNODE2 = NON_EXISTENT_NODE; // сосед отсутствует,
								sosedi[GG][maxelm].bdroblenie4 = true;
								#if doubleintprecision == 1
									printf("error NULL sosedi S1 elm=%lld\n", maxelm);
								#else
									printf("error NULL sosedi S1 elm=%d\n", maxelm);
								#endif
								
								getchar();

								// Эта ситуация возможна и здесь надо предусмотреть её обработку.
								DEFINE_BOUNDARY_PROPERTIES_TEMP(prop, prop_b, maxelm, maxbound); // В граничный узел сносятся свойства прилегающего КО.
								
								doublereal dS_loc = 0.0;

								// G - грань по которой ставится граничное условие.
								// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
								// bound_id - номер граничного узла, начиная с нуля.
								// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
								if (octree1->b4N) {
									obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
										octree1->linkN3->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
								}
								else {
									if (octree1->linkN == NULL) {
										obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
											 NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
									}
									else {
										obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
											octree1->linkN->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
									}
								}

								maxbound++;
								*/
							}
							else {
								if (octree1->linkS1->inum_TD == 0) {
									
									G = SSIDE;
									// граничная грань:
									sosedi[G][maxelm].iNODE2 = maxelm_memo + maxbound; // граничные КО нумеруются в последнюю очередь,
									sosedi[G][maxelm].bdroblenie4 = true;
									// после maxelm_memo внутренних КО.
									// Вычисление дальнего соседа
									CALC_GG(G, GG); // Вычисляет грань GG по грани G.

									sosedi[GG][maxelm].iNODE2 = NON_EXISTENT_NODE; // сосед отсутствует,
									sosedi[GG][maxelm].bdroblenie4 = true;
#if doubleintprecision == 1
									//printf("error inum_TD == 0 sosedi S1 elm=%lld\n", maxelm);
#else
									//printf("error inum_TD == 0 sosedi S1 elm=%d\n", maxelm);
#endif
									
									//getchar();

									// Эта ситуация возможна и здесь надо предусмотреть её обработку.
									DEFINE_BOUNDARY_PROPERTIES_TEMP(prop, prop_b, maxelm, maxbound); // В граничный узел сносятся свойства прилегающего КО.
									
									// Sx*Sz
									doublereal dS_loc = fabs(octree1->linkS1->p7.z - octree1->linkS1->p3.z)*fabs(octree1->linkS1->p2.x - octree1->linkS1->p3.x);
									TOCHKA p_centerG;
									p_centerG.z = 0.5*(octree1->linkS1->p7.z + octree1->linkS1->p3.z);
									p_centerG.y = octree1->linkS1->p2.y;
									p_centerG.x = 0.5*(octree1->linkS1->p2.x + octree1->linkS1->p3.x);

									// G - грань по которой ставится граничное условие.
									// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
									// bound_id - номер граничного узла, начиная с нуля.
									// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
									if (octree1->b4N) {
										obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
											octree1->linkN3->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
									}
									else {
										if (octree1->linkN == NULL) {
											obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
												 NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
										}
										else {
											obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
												octree1->linkN->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
										}
									}


									maxbound++;
								}
								else {
									// сосед существует.
									if (bvisit[octree1->linkS1->inum_TD - 1]) {
										// узел уже был посещён
										// Здесь неким образом модифицировалось gran_t
										// Только внутренняя грань.
										G = SSIDE;
										// граничная грань:
										sosedi[G][maxelm].iNODE2 = octree1->linkS1->inum_TD - 1; // граничные КО нумеруются в последнюю очередь,
										sosedi[G][maxelm].bdroblenie4 = true;
										// после maxelm_memo внутренних КО.
										// Вычисление дальнего соседа
										CALC_GG(G, GG); // Вычисляет грань GG по грани G.

										// Определение дальнего соседа :
										// TODO 13 сентября 2016.
									}
									else {
										// узнать координаты центра грани и ориентацию в пространстве
										x_c = 0.5*(octree1->linkS1->p3.x + octree1->linkS1->p2.x);
										y_c = octree1->p0.y;
										z_c = 0.5*(octree1->linkS1->p3.z + octree1->linkS1->p7.z);
										iplane = XZ;
										patch_maxbound(iplane, s, ls, x_c, y_c, z_c, maxbound, binc2, lsid2);
										if (binc2) {
#if doubleintprecision == 1
											printf("error internal source S1 elm=%lld\n", maxelm);
#else
											printf("error internal source S1 elm=%d\n", maxelm);
#endif
											
											//getchar();
											system("PAUSE");
										}
										else {
											// Только внутренняя грань.
											G = SSIDE;
											// граничная грань:
											sosedi[G][maxelm].iNODE2 = octree1->linkS1->inum_TD - 1; // граничные КО нумеруются в последнюю очередь,
											sosedi[G][maxelm].bdroblenie4 = true;
											// после maxelm_memo внутренних КО.
											// Вычисление дальнего соседа
											CALC_GG(G, GG); // Вычисляет грань GG по грани G.

											// Определение дальнего соседа :
											// TODO 13 сентября 2016.
										}
									}
								}
							}

							if (octree1->linkS4 == NULL) {

								// 21 сентября 2016 (вырождение ячейки).

								// Заглушка.
								G = SSIDE;
								// граничная грань:
								sosedi[G][maxelm].iNODE3 = NON_EXISTENT_NODE; // граничные КО нумеруются в последнюю очередь,
								sosedi[G][maxelm].bdroblenie4 = true;
								// после maxelm_memo внутренних КО.
								// Вычисление дальнего соседа
								CALC_GG(G, GG); // Вычисляет грань GG по грани G.
								sosedi[GG][maxelm].iNODE3 = NON_EXISTENT_NODE; // сосед отсутствует,
								sosedi[GG][maxelm].bdroblenie4 = true;


								/*
								G = SSIDE;
								// граничная грань:
								sosedi[G][maxelm].iNODE3 = maxelm_memo + maxbound; // граничные КО нумеруются в последнюю очередь,
								sosedi[G][maxelm].bdroblenie4 = true;
								// после maxelm_memo внутренних КО.
								// Вычисление дальнего соседа
								CALC_GG(G,GG); // Вычисляет грань GG по грани G.
								sosedi[GG][maxelm].iNODE3 = NON_EXISTENT_NODE; // сосед отсутствует,
								sosedi[GG][maxelm].bdroblenie4 = true;

								#if doubleintprecision == 1
										printf("error NULL sosedi S4 elm=%lld\n", maxelm);
								#else
										printf("error NULL sosedi S4 elm=%d\n", maxelm);
								#endif
								
								getchar();

								// Эта ситуация возможна и здесь надо предусмотреть её обработку.
								DEFINE_BOUNDARY_PROPERTIES_TEMP(prop, prop_b, maxelm, maxbound); // В граничный узел сносятся свойства прилегающего КО.
								
								doublereal dS_loc = 0.0;

								// G - грань по которой ставится граничное условие.
								// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
								// bound_id - номер граничного узла, начиная с нуля.
								// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
								if (octree1->b4N) {
									obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
										octree1->linkN3->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
								}
								else {
									if (octree1->linkN == NULL) {
										obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
											 NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
									}
									else {
										obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
											octree1->linkN->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
									}
								}

								maxbound++;
								*/
							}
							else {
								if (octree1->linkS4->inum_TD == 0) {
									
									G = SSIDE;
									// граничная грань:
									sosedi[G][maxelm].iNODE3 = maxelm_memo + maxbound; // граничные КО нумеруются в последнюю очередь,
									sosedi[G][maxelm].bdroblenie4 = true;
									// после maxelm_memo внутренних КО.
									// Вычисление дальнего соседа
									CALC_GG(G, GG); // Вычисляет грань GG по грани G.
									sosedi[GG][maxelm].iNODE3 = NON_EXISTENT_NODE; // сосед отсутствует,
									sosedi[GG][maxelm].bdroblenie4 = true;
									
#if doubleintprecision == 1
									//printf("error inum_TD == 0 sosedi S4 elm=%lld\n", maxelm);
#else
									//printf("error inum_TD == 0 sosedi S4 elm=%d\n", maxelm);
#endif
									
									//getchar();

									// Эта ситуация возможна и здесь надо предусмотреть её обработку.
									DEFINE_BOUNDARY_PROPERTIES_TEMP(prop, prop_b, maxelm, maxbound); // В граничный узел сносятся свойства прилегающего КО.
									
									// Sx*Sz
									doublereal dS_loc = fabs(octree1->linkS4->p7.z - octree1->linkS4->p3.z)*fabs(octree1->linkS4->p2.x - octree1->linkS4->p3.x);
									TOCHKA p_centerG;
									p_centerG.z = 0.5*(octree1->linkS4->p7.z + octree1->linkS4->p3.z);
									p_centerG.y = octree1->linkS4->p2.y;
									p_centerG.x = 0.5*(octree1->linkS4->p2.x + octree1->linkS4->p3.x);

									// G - грань по которой ставится граничное условие.
									// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
									// bound_id - номер граничного узла, начиная с нуля.
									// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
									if (octree1->b4N) {
										obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
											octree1->linkN3->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
									}
									else {
										if (octree1->linkN == NULL) {
											obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
												 NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
										
										}
										else {
											obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
												octree1->linkN->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
										}
									}

									maxbound++;
								}
								else {
									// сосед существует.
									if (bvisit[octree1->linkS4->inum_TD - 1]) {
										// узел уже был посещён
										// Здесь неким образом модифицировалось gran_t
										// Только внутренняя грань.
										G = SSIDE;
										// граничная грань:
										sosedi[G][maxelm].iNODE3 = octree1->linkS4->inum_TD - 1; // граничные КО нумеруются в последнюю очередь,
										sosedi[G][maxelm].bdroblenie4 = true;
										// после maxelm_memo внутренних КО.
										// Вычисление дальнего соседа
										CALC_GG(G, GG); // Вычисляет грань GG по грани G.

										// Определение дальнего соседа :
										// TODO 13 сентября 2016.
									}
									else {
										// узнать координаты центра грани и ориентацию в пространстве
										x_c = 0.5*(octree1->linkS4->p3.x + octree1->linkS4->p2.x);
										y_c = octree1->p0.y;
										z_c = 0.5*(octree1->linkS4->p3.z + octree1->linkS4->p7.z);
										iplane = XZ;
										patch_maxbound(iplane, s, ls, x_c, y_c, z_c, maxbound, binc3, lsid3);
										if (binc3) {
#if doubleintprecision == 1
											printf("error internal source S4 elm=%lld\n", maxelm);
#else
											printf("error internal source S4 elm=%d\n", maxelm);
#endif
											
											//getchar();
											system("PAUSE");
										}
										else {
											// Только внутренняя грань.
											G = SSIDE;
											// граничная грань:
											sosedi[G][maxelm].iNODE3 = octree1->linkS4->inum_TD - 1; // граничные КО нумеруются в последнюю очередь,
											sosedi[G][maxelm].bdroblenie4 = true;
											// после maxelm_memo внутренних КО.
											// Вычисление дальнего соседа
											CALC_GG(G, GG); // Вычисляет грань GG по грани G.

											// Определение дальнего соседа :
											// TODO 13 сентября 2016.
										}
									}
								}
							}

							if (octree1->linkS5 == NULL) {
								
								// 21 сентября 2016 (вырождение ячейки).

								// Заглушка.
								G = SSIDE;
								// граничная грань:
								sosedi[G][maxelm].iNODE4 = NON_EXISTENT_NODE; // граничные КО нумеруются в последнюю очередь,
								sosedi[G][maxelm].bdroblenie4 = true;
								// после maxelm_memo внутренних КО.
								// Вычисление дальнего соседа
								CALC_GG(G, GG); // Вычисляет грань GG по грани G.

								sosedi[GG][maxelm].iNODE4 = NON_EXISTENT_NODE; // сосед отсутствует,
								sosedi[GG][maxelm].bdroblenie4 = true;

								/*
								G = SSIDE;
								// граничная грань:
								sosedi[G][maxelm].iNODE4 = maxelm_memo + maxbound; // граничные КО нумеруются в последнюю очередь,
								sosedi[G][maxelm].bdroblenie4 = true;
								// после maxelm_memo внутренних КО.
								// Вычисление дальнего соседа
								CALC_GG(G,GG); // Вычисляет грань GG по грани G.
								
								sosedi[GG][maxelm].iNODE4 = NON_EXISTENT_NODE; // сосед отсутствует,
								sosedi[GG][maxelm].bdroblenie4 = true;

								#if doubleintprecision == 1
										printf("error NULL sosedi S5 elm=%lld\n", maxelm);
								#else
										printf("error NULL sosedi S5 elm=%d\n", maxelm);
								#endif
								
								getchar();

								// Эта ситуация возможна и здесь надо предусмотреть её обработку.
								DEFINE_BOUNDARY_PROPERTIES_TEMP(prop, prop_b, maxelm, maxbound); // В граничный узел сносятся свойства прилегающего КО.
								
								doublereal dS_loc = 0.0;

								// G - грань по которой ставится граничное условие.
								// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
								// bound_id - номер граничного узла, начиная с нуля.
								// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
								if (octree1->b4N) {
									obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
										octree1->linkN3->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
								}
								else {
									if (octree1->linkN == NULL) {
										obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
											 NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
									}
									else {
										obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
											octree1->linkN->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
									}
								}

								maxbound++;
								*/
							}
							else {
								if (octree1->linkS5->inum_TD == 0) {
									
									G = SSIDE;
									// граничная грань:
									sosedi[G][maxelm].iNODE4 = maxelm_memo + maxbound; // граничные КО нумеруются в последнюю очередь,
									sosedi[G][maxelm].bdroblenie4 = true;
									// после maxelm_memo внутренних КО.
									// Вычисление дальнего соседа
									CALC_GG(G, GG); // Вычисляет грань GG по грани G.

									sosedi[GG][maxelm].iNODE4 = NON_EXISTENT_NODE; // сосед отсутствует,
									sosedi[GG][maxelm].bdroblenie4 = true;

#if doubleintprecision == 1
									//printf("error inum_TD == 0 sosedi S5 elm=%lld\n", maxelm);
#else
									//printf("error inum_TD == 0 sosedi S5 elm=%d\n", maxelm);
#endif
									
									//getchar();

									// Эта ситуация возможна и здесь надо предусмотреть её обработку.
									DEFINE_BOUNDARY_PROPERTIES_TEMP(prop, prop_b, maxelm, maxbound); // В граничный узел сносятся свойства прилегающего КО.
									
									// Sx*Sz
									doublereal dS_loc = fabs(octree1->linkS5->p7.z - octree1->linkS5->p3.z)*fabs(octree1->linkS5->p2.x - octree1->linkS5->p3.x);
									TOCHKA p_centerG;
									p_centerG.z = 0.5*(octree1->linkS5->p7.z + octree1->linkS5->p3.z);
									p_centerG.y = octree1->linkS5->p2.y;
									p_centerG.x = 0.5*(octree1->linkS5->p2.x + octree1->linkS5->p3.x);

									// G - грань по которой ставится граничное условие.
									// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
									// bound_id - номер граничного узла, начиная с нуля.
									// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
									if (octree1->b4N) {
										obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
											octree1->linkN3->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
									}
									else {
										if (octree1->linkN == NULL) {
											obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
												 NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
										}
										else {
											obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
												octree1->linkN->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
										}
									}

									maxbound++;
								}
								else {
									// сосед существует.
									if (bvisit[octree1->linkS5->inum_TD - 1]) {
										// узел уже был посещён
										// Здесь неким образом модифицировалось gran_t
										// Только внутренняя грань.
										G = SSIDE;
										// граничная грань:
										sosedi[G][maxelm].iNODE4 = octree1->linkS5->inum_TD - 1; // граничные КО нумеруются в последнюю очередь,
										sosedi[G][maxelm].bdroblenie4 = true;
										// после maxelm_memo внутренних КО.
										// Вычисление дальнего соседа
										CALC_GG(G, GG); // Вычисляет грань GG по грани G.

										// Определение дальнего соседа :
										// TODO 13 сентября 2016.
									}
									else {
										// узнать координаты центра грани и ориентацию в пространстве
										x_c = 0.5*(octree1->linkS5->p3.x + octree1->linkS5->p2.x);
										y_c = octree1->p0.y;
										z_c = 0.5*(octree1->linkS5->p3.z + octree1->linkS5->p7.z);
										iplane = XZ;
										patch_maxbound(iplane, s, ls, x_c, y_c, z_c, maxbound, binc4, lsid4);
										if (binc4) {
#if doubleintprecision == 1
											printf("error internal source S5 elm=%lld\n", maxelm);
#else
											printf("error internal source S5 elm=%d\n", maxelm);
#endif
											
											//getchar();
											system("PAUSE");
										}
										else {
											// Только внутренняя грань.
											G = SSIDE;
											// граничная грань:
											sosedi[G][maxelm].iNODE4 = octree1->linkS5->inum_TD - 1; // граничные КО нумеруются в последнюю очередь,
											sosedi[G][maxelm].bdroblenie4 = true;
											// после maxelm_memo внутренних КО.
											// Вычисление дальнего соседа
											CALC_GG(G, GG); // Вычисляет грань GG по грани G.

											// Определение дальнего соседа :
											// TODO 13 сентября 2016.
										}
									}
								}
							}
						}
					}

					// сосед не существует
					// записываем идентификатор грани и
					// увеличиваем счётчик граней.
					if (!octree1->b4T) {
						if ((octree1->linkT == NULL)||((octree1->linkT != NULL)&&(octree1->linkT->inum_TD == 0))) {
							// Если мы находимся на границе расчётной области или если мы граничим с HOLLOW блоком.

							G = TSIDE;
							// граничная грань:
							INIT_SOSEDI(sosedi, maxelm, TSIDE, GG, false, maxelm_memo + maxbound,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE, false,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE); // заполнение sosedi[G][maxelm].
							DEFINE_BOUNDARY_PROPERTIES_TEMP(prop, prop_b, maxelm, maxbound); // В граничный узел сносятся свойства прилегающего КО.
							
							doublereal dSloc = CALC_SQUARE(G, maxelm, nvtx, pa); // Вычисление площади грани ячейки.
							TOCHKA p_centerG = CALC_CENTERG(G, maxelm, nvtx, pa); // Вычисление геометрического центра грани.
							
							// G - грань по которой ставится граничное условие.
							// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
							// bound_id - номер граничного узла, начиная с нуля.
							// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
							if (octree1->b4B) {
								obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
									octree1->linkB0->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, false, false, p_centerG);
							}
							else {
								if (octree1->linkB == NULL) {
									obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
										 NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, false, false, p_centerG);
								}
								else {
									obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
										octree1->linkB->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, false, false, p_centerG);
								}
								
							}

							maxbound++;
						}
						else {
							if (octree1->linkT != NULL) {
								// сосед существует.
								if (bvisit[octree1->linkT->inum_TD - 1]) {


									// return_current_node_number - Заглянули в соседа и вернулись из соседа в текущий узел.
									// Номер должен быть тем же самым и поэтому этот код бессмысленен.
									// Строго внутренние плоские бесконечно тонкие источники тепла внутри расчётной области больше не поддерживаются.
									// Поэтому от этого кода можно безболезненно избавиться. Поддерживаются только объёмные источники тепла.
									// 10.03.2019
									integer return_current_node_number = sosedi[BSIDE][octree1->linkT->inum_TD - 1].iNODE1;
									if (octree1->linkT->b4B) {
										integer current_node_number = octree1->inum_TD - 1;
										bool bcontinue = true;
										if (bcontinue) {
											if (octree1->linkT->linkB0 != NULL) {
												if (octree1->linkT->linkB0->inum_TD - 1 == current_node_number) {
													bcontinue = false; // iNODE1
													if (sosedi[BSIDE][octree1->linkT->inum_TD - 1].iNODE1 == current_node_number) {
														return_current_node_number = sosedi[BSIDE][octree1->linkT->inum_TD - 1].iNODE1;
													}
													else {
														printf("this can not be T B0\n");
														//getchar();
														system("PAUSE");
													}
												}
											}
										}
										if (bcontinue) {
											if (octree1->linkT->linkB1 != NULL) {
												if (octree1->linkT->linkB1->inum_TD - 1 == current_node_number) {
													bcontinue = false; // iNODE2
													if (sosedi[BSIDE][octree1->linkT->inum_TD - 1].iNODE2 == current_node_number) {
														return_current_node_number = sosedi[BSIDE][octree1->linkT->inum_TD - 1].iNODE2;
													}
													else {
														printf("this can not be T B1\n");
														//getchar();
														system("PAUSE");
													}
												}
											}
										}
										if (bcontinue) {
											if (octree1->linkT->linkB2 != NULL) {
												if (octree1->linkT->linkB2->inum_TD - 1 == current_node_number) {
													bcontinue = false; // iNODE3
													if (sosedi[BSIDE][octree1->linkT->inum_TD - 1].iNODE3 == current_node_number) {
														return_current_node_number = sosedi[BSIDE][octree1->linkT->inum_TD - 1].iNODE3;
													}
													else {
														printf("this can not be T B2\n");
														//getchar();
														system("PAUSE");
													}
												}
											}
										}
										if (bcontinue) {
											if (octree1->linkT->linkB3 != NULL) {
												if (octree1->linkT->linkB3->inum_TD - 1 == current_node_number) {
													bcontinue = false; // iNODE4
													if (sosedi[BSIDE][octree1->linkT->inum_TD - 1].iNODE4 == current_node_number) {
														return_current_node_number = sosedi[BSIDE][octree1->linkT->inum_TD - 1].iNODE4;
													}
													else {
														printf("this can not be T B3\n");
														//getchar();
														system("PAUSE");
													}
												}
											}
										}

									}

									// узел уже был посещён
									// Здесь неким образом модифицировалось gran_t
									if (return_current_node_number >= maxelm_memo) {

										printf("WARNING!!! internal source gran TSIDE 10.03.2019\n");
										system("PAUSE");

										// Это граничная грань внутреннего источника тепла.
										// Внутренний источник тепла.
										G = TSIDE;
										// граничная грань:
										INIT_SOSEDI(sosedi, maxelm, TSIDE, GG, false, sosedi[BSIDE][octree1->linkT->inum_TD - 1].iNODE1,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE, false,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE); // заполнение sosedi[G][maxelm].

										if (bi_fluid) {
											// FLUID
											// В граничный узел сносятся свойства прилегающего КО.
											DEFINE_BOUNDARY_PROPERTIES_TEMP(prop, prop_b, octree1->linkT->inum_TD - 1, sosedi[G][maxelm].iNODE1 - maxelm_memo);
										}
										else {
											// SOLID
											// В граничный узел сносятся свойства прилегающего КО.
											DEFINE_BOUNDARY_PROPERTIES_TEMP(prop, prop_b, maxelm, sosedi[G][maxelm].iNODE1 - maxelm_memo);
										}

										doublereal dSloc = CALC_SQUARE(G, maxelm, nvtx, pa); // Вычисление площади грани ячейки.
										TOCHKA p_centerG = CALC_CENTERG(G, maxelm, nvtx, pa); // Вычисление геометрического центра грани.
										
										// G - грань по которой ставится граничное условие.
										// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
										// bound_id - номер граничного узла, начиная с нуля.
										// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
										if (octree1->b4B) {
											obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, sosedi[BSIDE][octree1->linkT->inum_TD - 1].iNODE1 - maxelm_memo,
												octree1->linkB0->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, true, true, p_centerG);
										}
										else {
											if (octree1->linkB == NULL) {
												obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, sosedi[BSIDE][octree1->linkT->inum_TD - 1].iNODE1 - maxelm_memo,
													 NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, true, true, p_centerG);
											}
											else {
												obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, sosedi[BSIDE][octree1->linkT->inum_TD - 1].iNODE1 - maxelm_memo,
													octree1->linkB->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, true, true, p_centerG);
											}
										}
										binternalsource[sosedi[BSIDE][octree1->linkT->inum_TD - 1].iNODE1 - maxelm_memo] = true; // внутренний источник тепла на границе жидкости и твёрдого тела

									}
									else {
										// Это строго внутренний узел.
										// Внутренняя грань.
										G = TSIDE;
										// Внутренняя грань:
										sosedi[G][maxelm].iNODE1 = octree1->linkT->inum_TD - 1; // граничные КО нумеруются в последнюю очередь,
										sosedi[G][maxelm].bdroblenie4 = false;
										sosedi[G][maxelm].iNODE2 = NON_EXISTENT_NODE;
										sosedi[G][maxelm].iNODE3 = NON_EXISTENT_NODE;
										sosedi[G][maxelm].iNODE4 = NON_EXISTENT_NODE;
										// после maxelm_memo внутренних КО.
										// Вычисление дальнего соседа
										//CALC_GG(G, GG); // Вычисляет грань GG по грани G.

										// Определение дальнего соседа :
										// TODO 13 сентября 2016.
									}

								}
								else {
									bool binc;
									integer lsid;
									// узнать координаты центра грани и ориентацию в пространстве
									x_c = 0.5*(octree1->p4.x + octree1->p5.x);
									y_c = 0.5*(octree1->p4.y + octree1->p7.y);
									z_c = octree1->p4.z;
									iplane = XY;
									patch_maxbound(iplane, s, ls, x_c, y_c, z_c, maxbound,binc,lsid);
									if (binc) {
										// Внутренний источник тепла.
										// Внутренние плоские бесконечно тонкие источники тепла скорее всего больше не поддерживаются. 10,03,2019

										G = TSIDE;
										// граничная грань:
										INIT_SOSEDI(sosedi, maxelm, TSIDE, GG, false, maxelm_memo + maxbound-1,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE, false,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE); // заполнение sosedi[G][maxelm].
										
										if (bi_fluid) {
											// FLUID
											// В граничный узел сносятся свойства прилегающего КО.
											DEFINE_BOUNDARY_PROPERTIES_TEMP(prop, prop_b, octree1->linkT->inum_TD - 1, maxbound-1);
										}
										else {
											// SOLID
											// В граничный узел сносятся свойства прилегающего КО.
											DEFINE_BOUNDARY_PROPERTIES_TEMP(prop, prop_b, maxelm, maxbound-1);
										}

										doublereal dSloc = CALC_SQUARE(G, maxelm, nvtx, pa); // Вычисление площади грани ячейки.
										TOCHKA p_centerG = CALC_CENTERG(G, maxelm, nvtx, pa); // Вычисление геометрического центра грани.
										
										// G - грань по которой ставится граничное условие.
										// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
										// bound_id - номер граничного узла, начиная с нуля.
										// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
										if (octree1->b4B) {
											obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound - 1,
												octree1->linkB0->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, true, false, p_centerG);
										}
										else {
											if (octree1->linkB == NULL) {
												obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound - 1,
													 NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, true, false, p_centerG);
											}
											else {
												obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound - 1,
													octree1->linkB->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, true, false, p_centerG);
											}
										}
										binternalsource[maxbound - 1] = true; // внутренний источник тепла на границе жидкости и твёрдого тела
									}
									else {
										// Внутренняя грань.
										G = TSIDE;
										// Внутренняя грань:
										sosedi[G][maxelm].iNODE1 = octree1->linkT->inum_TD - 1; // граничные КО нумеруются в последнюю очередь,
										sosedi[G][maxelm].bdroblenie4 = false;
										sosedi[G][maxelm].iNODE2 = NON_EXISTENT_NODE;
										sosedi[G][maxelm].iNODE3 = NON_EXISTENT_NODE;
										sosedi[G][maxelm].iNODE4 = NON_EXISTENT_NODE;
										// после maxelm_memo внутренних КО.
										// Вычисление дальнего соседа
										//CALC_GG(G, GG); // Вычисляет грань GG по грани G.

										// Определение дальнего соседа :
										// TODO 13 сентября 2016.
									}
								}
							}
						}
					}
					else {
						// Присутствуют все 4 соседа: T: 4,5,6,7.

						if (((octree1->linkT4 == NULL) && (octree1->linkT5 == NULL) && (octree1->linkT6 == NULL) && (octree1->linkT7 == NULL))||
							(((octree1->linkT4 != NULL) && (octree1->linkT5 != NULL) && (octree1->linkT6 != NULL) && (octree1->linkT7 != NULL) &&
							(octree1->linkT4->inum_TD == 0) && (octree1->linkT5->inum_TD == 0) && (octree1->linkT6->inum_TD == 0) && (octree1->linkT7->inum_TD == 0)))||
								(((octree1->linkT4 != NULL) && (octree1->linkT5 == NULL) && (octree1->linkT6 == NULL) && (octree1->linkT7 != NULL) &&
							(octree1->linkT4->inum_TD == 0) && (octree1->linkT7->inum_TD == 0)))||
									(((octree1->linkT4 != NULL) && (octree1->linkT5 != NULL) && (octree1->linkT6 == NULL) && (octree1->linkT7 == NULL) &&
							(octree1->linkT4->inum_TD == 0) && (octree1->linkT5->inum_TD == 0)))||
							(((octree1->linkT4 != NULL) && (octree1->linkT5 == NULL) && (octree1->linkT6 == NULL) && (octree1->linkT7 == NULL) &&
							(octree1->linkT4->inum_TD == 0))))
						{

							// 1 - T4&&T5&&T6&&T7 - на границе расчётной области (РО).
							// 2 - T4&&T5&&T6&&T7 - на гранеице HOLLOW блока.
							// 3 - T4&&T7 на границе HOLLOW блока, T5&&T6 на границе РО.
							// 4 - T4&&T5 на границе HOLLOW блока, T6&T7 на границе РО.
							// 5 - T4 на границе HOLLOW блока, T5&&T6&&T7  на границе РО.

							G = TSIDE;
							// граничная грань:
							INIT_SOSEDI(sosedi, maxelm, TSIDE, GG, false, maxelm_memo + maxbound,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE, false,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE); // заполнение sosedi[G][maxelm].
							DEFINE_BOUNDARY_PROPERTIES_TEMP(prop, prop_b, maxelm, maxbound); // В граничный узел сносятся свойства прилегающего КО.
							
							doublereal dSloc = CALC_SQUARE(G, maxelm, nvtx, pa); // Вычисление площади грани ячейки.
							TOCHKA p_centerG = CALC_CENTERG(G, maxelm, nvtx, pa); // Вычисление геометрического центра грани.
							
							// G - грань по которой ставится граничное условие.
							// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
							// bound_id - номер граничного узла, начиная с нуля.
							// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
							if (octree1->b4B) {
								obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
									octree1->linkB0->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, false, false, p_centerG);
							}
							else {
								if (octree1->linkB == NULL) {
									obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
										 NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, false, false, p_centerG);
								}
								else {
									obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
										octree1->linkB->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, false, false, p_centerG);
								}
							}

							// Зачем дробить если граница заканчивается NULL.
							maxbound++; // Вместо четырёх один.
						}						
						else {

							bool binc1, binc2, binc3, binc4;
							integer lsid1, lsid2, lsid3, lsid4;

							if (octree1->linkT4 == NULL) {
#if doubleintprecision == 1
								printf("error NULL sosedi T4 elm=%lld\n", maxelm);
#else
								printf("error NULL sosedi T4 elm=%d\n", maxelm);
#endif
								
								//getchar();
								system("PAUSE");

								/*

								G = TSIDE;
								// граничная грань:
								sosedi[G][maxelm].iNODE1 = maxelm_memo + maxbound; // граничные КО нумеруются в последнюю очередь,
								sosedi[G][maxelm].bdroblenie4 = true;
								// после maxelm_memo внутренних КО.
								// Вычисление дальнего соседа
								CALC_GG(G,GG); // Вычисляет грань GG по грани G.

								sosedi[GG][maxelm].iNODE1 = NON_EXISTENT_NODE; // сосед отсутствует,
								sosedi[GG][maxelm].bdroblenie4 = true;

								#if doubleintprecision == 1
										printf("error NULL sosedi T4 elm=%lld\n", maxelm);
								#else
										printf("error NULL sosedi T4 elm=%d\n", maxelm);
								#endif
								
								getchar();

								// Эта ситуация возможна и здесь надо предусмотреть её обработку.
								DEFINE_BOUNDARY_PROPERTIES_TEMP(prop, prop_b, maxelm, maxbound); // В граничный узел сносятся свойства прилегающего КО.
								
								doublereal dS_loc = 0.0;

								// G - грань по которой ставится граничное условие.
								// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
								// bound_id - номер граничного узла, начиная с нуля.
								// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
								if (octree1->b4B) {
									obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
										octree1->linkB0->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
								}
								else {
									if (octree1->linkB == NULL) {
										obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
											 NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
									}
									else {
										obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
											octree1->linkB->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
									}

								}

								maxbound++;
								*/
							}
							else {
								if (octree1->linkT4->inum_TD == 0) {

									G = TSIDE;
									// граничная грань:
									sosedi[G][maxelm].iNODE1 = maxelm_memo + maxbound; // граничные КО нумеруются в последнюю очередь,
									sosedi[G][maxelm].bdroblenie4 = true;
									// после maxelm_memo внутренних КО.
									// Вычисление дальнего соседа
									CALC_GG(G, GG); // Вычисляет грань GG по грани G.

									sosedi[GG][maxelm].iNODE1 = NON_EXISTENT_NODE; // сосед отсутствует,
									sosedi[GG][maxelm].bdroblenie4 = true;
#if doubleintprecision == 1
									//printf("error inum_TD == 0 sosedi T4 elm=%lld\n", maxelm);
#else
									//printf("error inum_TD == 0 sosedi T4 elm=%d\n", maxelm);
#endif
									
									//getchar();

									// Эта ситуация возможна и здесь надо предусмотреть её обработку.
									DEFINE_BOUNDARY_PROPERTIES_TEMP(prop, prop_b, maxelm, maxbound); // В граничный узел сносятся свойства прилегающего КО.
									
									// Sx*Sy
									doublereal dS_loc = fabs(octree1->linkT4->p3.y - octree1->linkT4->p0.y)*fabs(octree1->linkT4->p1.x - octree1->linkT4->p0.x);
									TOCHKA p_centerG;
									p_centerG.z = octree1->linkT4->p0.z;
									p_centerG.y = 0.5*(octree1->linkT4->p3.y + octree1->linkT4->p0.y);
									p_centerG.x = 0.5*(octree1->linkT4->p1.x + octree1->linkT4->p0.x);

									// G - грань по которой ставится граничное условие.
									// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
									// bound_id - номер граничного узла, начиная с нуля.
									// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
									if (octree1->b4B) {
										obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
											octree1->linkB0->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
									}
									else {
										if (octree1->linkB == NULL) {
											obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
												 NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
										}
										else {
											obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
												octree1->linkB->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
										}

									}

									maxbound++;
								}
								else {
									// сосед существует.
									if (bvisit[octree1->linkT4->inum_TD - 1]) {
										// узел уже был посещён
										// Здесь неким образом модифицировалось gran_t
										// Только внутренняя грань.
										G = TSIDE;
										// граничная грань:
										sosedi[G][maxelm].iNODE1 = octree1->linkT4->inum_TD - 1; // граничные КО нумеруются в последнюю очередь,
										sosedi[G][maxelm].bdroblenie4 = true;
										// после maxelm_memo внутренних КО.
										// Вычисление дальнего соседа
										CALC_GG(G, GG); // Вычисляет грань GG по грани G.

										// Определение дальнего соседа :
										// TODO 13 сентября 2016.
									}
									else {
										// узнать координаты центра грани и ориентацию в пространстве
										x_c = 0.5*(octree1->linkT4->p0.x + octree1->linkT4->p1.x);
										y_c = 0.5*(octree1->linkT4->p0.y + octree1->linkT4->p3.y);
										z_c = octree1->p4.z;
										iplane = XY;
										patch_maxbound(iplane, s, ls, x_c, y_c, z_c, maxbound, binc1, lsid1);
										if (binc1) {
#if doubleintprecision == 1
											printf("error internal source T4 elm=%lld\n", maxelm);
#else
											printf("error internal source T4 elm=%d\n", maxelm);
#endif
											
											//getchar();
											system("PAUSE");
										}
										else {
											// Только внутренняя грань.
											G = TSIDE;
											// граничная грань:
											sosedi[G][maxelm].iNODE1 = octree1->linkT4->inum_TD - 1; // граничные КО нумеруются в последнюю очередь,
											sosedi[G][maxelm].bdroblenie4 = true;
											// после maxelm_memo внутренних КО.
											// Вычисление дальнего соседа
											CALC_GG(G, GG); // Вычисляет грань GG по грани G.

											// Определение дальнего соседа :
											// TODO 13 сентября 2016.
										}
									}
								}
							}

							if (octree1->linkT5 == NULL) {

								// 21 сентября 2016 (вырождение ячейки).

								// Заглушка.
								G = TSIDE;
								// граничная грань:
								sosedi[G][maxelm].iNODE2 = NON_EXISTENT_NODE; // граничные КО нумеруются в последнюю очередь,
								sosedi[G][maxelm].bdroblenie4 = true;
								// после maxelm_memo внутренних КО.
								// Вычисление дальнего соседа
								CALC_GG(G, GG); // Вычисляет грань GG по грани G.

								sosedi[GG][maxelm].iNODE2 = NON_EXISTENT_NODE; // соседа нет,
								sosedi[GG][maxelm].bdroblenie4 = true;

								/*
								G = TSIDE;
								// граничная грань:
								sosedi[G][maxelm].iNODE2 = maxelm_memo + maxbound; // граничные КО нумеруются в последнюю очередь,
								sosedi[G][maxelm].bdroblenie4 = true;
								// после maxelm_memo внутренних КО.
								// Вычисление дальнего соседа
								CALC_GG(G,GG); // Вычисляет грань GG по грани G.

								sosedi[GG][maxelm].iNODE2 = NON_EXISTENT_NODE; // соседа нет,
								sosedi[GG][maxelm].bdroblenie4 = true;

								#if doubleintprecision == 1
										printf("error NULL sosedi T5 elm=%lld\n", maxelm);
								#else
										printf("error NULL sosedi T5 elm=%d\n", maxelm);
								#endif
								
								getchar();
								
								// Эта ситуация возможна и здесь надо предусмотреть её обработку.
								DEFINE_BOUNDARY_PROPERTIES_TEMP(prop, prop_b, maxelm, maxbound); // В граничный узел сносятся свойства прилегающего КО.
								
								doublereal dS_loc = 0.0;

								// G - грань по которой ставится граничное условие.
								// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
								// bound_id - номер граничного узла, начиная с нуля.
								// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
								if (octree1->b4B) {
									obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
										octree1->linkB0->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
								}
								else {
									if (octree1->linkB == NULL) {
										obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
											 NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
									}
									else {
										obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
											octree1->linkB->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
									}

								}

								maxbound++;
								*/
							}
							else {
								if (octree1->linkT5->inum_TD == 0) {
									
									G = TSIDE;
									// граничная грань:
									sosedi[G][maxelm].iNODE2 = maxelm_memo + maxbound; // граничные КО нумеруются в последнюю очередь,
									sosedi[G][maxelm].bdroblenie4 = true;
									// после maxelm_memo внутренних КО.
									// Вычисление дальнего соседа
									CALC_GG(G, GG); // Вычисляет грань GG по грани G.

									sosedi[GG][maxelm].iNODE2 = NON_EXISTENT_NODE; // соседа нет,
									sosedi[GG][maxelm].bdroblenie4 = true;

#if doubleintprecision == 1
									//printf("error inum_TD == 0 sosedi T5 elm=%lld\n", maxelm);
#else
									//printf("error inum_TD == 0 sosedi T5 elm=%d\n", maxelm);
#endif
									
									//getchar();

									// Эта ситуация возможна и здесь надо предусмотреть её обработку.
									DEFINE_BOUNDARY_PROPERTIES_TEMP(prop, prop_b, maxelm, maxbound); // В граничный узел сносятся свойства прилегающего КО.
									
									// Sx*Sy
									doublereal dS_loc = fabs(octree1->linkT5->p3.y - octree1->linkT5->p0.y)*fabs(octree1->linkT5->p1.x - octree1->linkT5->p0.x);
									TOCHKA p_centerG;
									p_centerG.z = octree1->linkT5->p0.z;
									p_centerG.y = 0.5*(octree1->linkT5->p3.y + octree1->linkT5->p0.y);
									p_centerG.x = 0.5*(octree1->linkT5->p1.x + octree1->linkT5->p0.x);

									// G - грань по которой ставится граничное условие.
									// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
									// bound_id - номер граничного узла, начиная с нуля.
									// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
									if (octree1->b4B) {
										obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
											octree1->linkB0->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
									}
									else {
										if (octree1->linkB == NULL) {
											obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
												 NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
										}
										else {
											obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
												octree1->linkB->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
										}

									}

									maxbound++;
								}
								else {
									// сосед существует.
									if (bvisit[octree1->linkT5->inum_TD - 1]) {
										// узел уже был посещён
										// Здесь неким образом модифицировалось gran_t
										// Только внутренняя грань.
										G = TSIDE;
										// граничная грань:
										sosedi[G][maxelm].iNODE2 = octree1->linkT5->inum_TD - 1; // граничные КО нумеруются в последнюю очередь,
										sosedi[G][maxelm].bdroblenie4 = true;
										// после maxelm_memo внутренних КО.
										// Вычисление дальнего соседа
										CALC_GG(G, GG); // Вычисляет грань GG по грани G.

										// Определение дальнего соседа :
										// TODO 13 сентября 2016.
									}
									else {
										// узнать координаты центра грани и ориентацию в пространстве
										x_c = 0.5*(octree1->linkT5->p0.x + octree1->linkT5->p1.x);
										y_c = 0.5*(octree1->linkT5->p0.y + octree1->linkT5->p3.y);
										z_c = octree1->p4.z;
										iplane = XY;
										patch_maxbound(iplane, s, ls, x_c, y_c, z_c, maxbound, binc2, lsid2);
										if (binc2) {
#if doubleintprecision == 1
											printf("error internal source T5 elm=%lld\n", maxelm);
#else
											printf("error internal source T5 elm=%d\n", maxelm);
#endif
											
											//getchar();
											system("PAUSE");
										}
										else {
											// Только внутренняя грань.
											G = TSIDE;
											// граничная грань:
											sosedi[G][maxelm].iNODE2 = octree1->linkT5->inum_TD - 1; // граничные КО нумеруются в последнюю очередь,
											sosedi[G][maxelm].bdroblenie4 = true;
											// после maxelm_memo внутренних КО.
											// Вычисление дальнего соседа
											CALC_GG(G, GG); // Вычисляет грань GG по грани G.

											// Определение дальнего соседа :
											// TODO 13 сентября 2016.
										}
									}
								}
							}

							if (octree1->linkT6 == NULL) {

								// 21 сентября 2016 (вырождение ячейки).

								// Заглушка.
								G = TSIDE;
								// граничная грань:
								sosedi[G][maxelm].iNODE3 = NON_EXISTENT_NODE; // граничные КО нумеруются в последнюю очередь,
								sosedi[G][maxelm].bdroblenie4 = true;
								// после maxelm_memo внутренних КО.
								// Вычисление дальнего соседа
								CALC_GG(G, GG); // Вычисляет грань GG по грани G.

								sosedi[GG][maxelm].iNODE3 = NON_EXISTENT_NODE; // сосед отсутствует,
								sosedi[GG][maxelm].bdroblenie4 = true;

								/*
								G = TSIDE;
								// граничная грань:
								sosedi[G][maxelm].iNODE3 = maxelm_memo + maxbound; // граничные КО нумеруются в последнюю очередь,
								sosedi[G][maxelm].bdroblenie4 = true;
								// после maxelm_memo внутренних КО.
								// Вычисление дальнего соседа
								CALC_GG(G,GG); // Вычисляет грань GG по грани G.

								sosedi[GG][maxelm].iNODE3 = NON_EXISTENT_NODE; // сосед отсутствует,
								sosedi[GG][maxelm].bdroblenie4 = true;

								#if doubleintprecision == 1
										printf("error NULL sosedi T6 elm=%lld\n", maxelm);
								#else
										printf("error NULL sosedi T6 elm=%d\n", maxelm);
								#endif
								
								getchar();

								// Эта ситуация возможна и здесь надо предусмотреть её обработку.
								DEFINE_BOUNDARY_PROPERTIES_TEMP(prop, prop_b, maxelm, maxbound); // В граничный узел сносятся свойства прилегающего КО.
								
								doublereal dS_loc = 0.0;

								// G - грань по которой ставится граничное условие.
								// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
								// bound_id - номер граничного узла, начиная с нуля.
								// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
								if (octree1->b4B) {
									obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
										octree1->linkB0->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
								}
								else {
									if (octree1->linkB == NULL) {
										obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
											 NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
									}
									else {
										obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
											octree1->linkB->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
									}

								}

								maxbound++;
								*/
							}
							else {
								if (octree1->linkT6->inum_TD == 0) {
									
									G = TSIDE;
									// граничная грань:
									sosedi[G][maxelm].iNODE3 = maxelm_memo + maxbound; // граничные КО нумеруются в последнюю очередь,
									sosedi[G][maxelm].bdroblenie4 = true;
									// после maxelm_memo внутренних КО.
									// Вычисление дальнего соседа
									CALC_GG(G, GG); // Вычисляет грань GG по грани G.

									sosedi[GG][maxelm].iNODE3 = NON_EXISTENT_NODE; // сосед отсутствует,
									sosedi[GG][maxelm].bdroblenie4 = true;

#if doubleintprecision == 1
									//printf("error inum_TD == 0 sosedi T6 elm=%lld\n", maxelm);
#else
									//printf("error inum_TD == 0 sosedi T6 elm=%d\n", maxelm);
#endif
									
									//getchar();

									// Эта ситуация возможна и здесь надо предусмотреть её обработку.
									DEFINE_BOUNDARY_PROPERTIES_TEMP(prop, prop_b, maxelm, maxbound); // В граничный узел сносятся свойства прилегающего КО.
									
									// Sx*Sy
									doublereal dS_loc = fabs(octree1->linkT6->p3.y - octree1->linkT6->p0.y)*fabs(octree1->linkT6->p1.x - octree1->linkT6->p0.x);
									TOCHKA p_centerG;
									p_centerG.z = octree1->linkT6->p0.z;
									p_centerG.y = 0.5*(octree1->linkT6->p3.y + octree1->linkT6->p0.y);
									p_centerG.x = 0.5*(octree1->linkT6->p1.x + octree1->linkT6->p0.x);

									// G - грань по которой ставится граничное условие.
									// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
									// bound_id - номер граничного узла, начиная с нуля.
									// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
									if (octree1->b4B) {
										obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
											octree1->linkB0->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
									}
									else {
										if (octree1->linkB == NULL) {
											obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
												 NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
										}
										else {
											obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
												octree1->linkB->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
										}

									}


									maxbound++;
								}
								else {
									// сосед существует.
									if (bvisit[octree1->linkT6->inum_TD - 1]) {
										// узел уже был посещён
										// Здесь неким образом модифицировалось gran_t
										// Только внутренняя грань.
										G = TSIDE;
										// граничная грань:
										sosedi[G][maxelm].iNODE3 = octree1->linkT6->inum_TD - 1; // граничные КО нумеруются в последнюю очередь,
										sosedi[G][maxelm].bdroblenie4 = true;
										// после maxelm_memo внутренних КО.
										// Вычисление дальнего соседа
										CALC_GG(G, GG); // Вычисляет грань GG по грани G.

										// Определение дальнего соседа :
										// TODO 13 сентября 2016.
									}
									else {
										// узнать координаты центра грани и ориентацию в пространстве
										x_c = 0.5*(octree1->linkT6->p0.x + octree1->linkT6->p1.x);
										y_c = 0.5*(octree1->linkT6->p0.y + octree1->linkT6->p3.y);
										z_c = octree1->p4.z;
										iplane = XY;
										patch_maxbound(iplane, s, ls, x_c, y_c, z_c, maxbound, binc3, lsid3);
										if (binc3) {
#if doubleintprecision == 1
											printf("error internal source T6 elm=%lld\n", maxelm);
#else
											printf("error internal source T6 elm=%d\n", maxelm);
#endif
											
											//getchar();
											system("PAUSE");
										}
										else {
											// Только внутренняя грань.
											G = TSIDE;
											// граничная грань:
											sosedi[G][maxelm].iNODE3 = octree1->linkT6->inum_TD - 1; // граничные КО нумеруются в последнюю очередь,
											sosedi[G][maxelm].bdroblenie4 = true;
											// после maxelm_memo внутренних КО.
											// Вычисление дальнего соседа
											CALC_GG(G, GG); // Вычисляет грань GG по грани G.

											// Определение дальнего соседа :
											// TODO 13 сентября 2016.
										}
									}
								}
							}

							if (octree1->linkT7 == NULL) {

								// 21 сентября 2016 (вырождение ячейки).

								// Заглушка.
								G = TSIDE;
								// граничная грань:
								sosedi[G][maxelm].iNODE4 = NON_EXISTENT_NODE; // граничные КО нумеруются в последнюю очередь,
								sosedi[G][maxelm].bdroblenie4 = true;
								// после maxelm_memo внутренних КО.
								// Вычисление дальнего соседа
								CALC_GG(G, GG); // Вычисляет грань GG по грани G.

								sosedi[GG][maxelm].iNODE4 = NON_EXISTENT_NODE; // соседа нет,
								sosedi[GG][maxelm].bdroblenie4 = true;

								/*
								G = TSIDE;
								// граничная грань:
								sosedi[G][maxelm].iNODE4 = maxelm_memo + maxbound; // граничные КО нумеруются в последнюю очередь,
								sosedi[G][maxelm].bdroblenie4 = true;
								// после maxelm_memo внутренних КО.
								// Вычисление дальнего соседа
								CALC_GG(G,GG); // Вычисляет грань GG по грани G.

								sosedi[GG][maxelm].iNODE4 = NON_EXISTENT_NODE; // соседа нет,
								sosedi[GG][maxelm].bdroblenie4 = true;

								#if doubleintprecision == 1
										printf("error NULL sosedi T7 elm=%lld\n", maxelm);
								#else
										printf("error NULL sosedi T7 elm=%d\n", maxelm);
								#endif
								
								getchar();

								// Эта ситуация возможна и здесь надо предусмотреть её обработку.
								DEFINE_BOUNDARY_PROPERTIES_TEMP(prop, prop_b, maxelm, maxbound); // В граничный узел сносятся свойства прилегающего КО.
								
								doublereal dS_loc = 0.0;

								// G - грань по которой ставится граничное условие.
								// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
								// bound_id - номер граничного узла, начиная с нуля.
								// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
								if (octree1->b4B) {
									obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
										octree1->linkB0->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
								}
								else {
									if (octree1->linkB == NULL) {
										obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
											 NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
									}
									else {
										obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
											octree1->linkB->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
									}

								}

								maxbound++;
								*/
							}
							else {
								if (octree1->linkT7->inum_TD == 0) {
									
									G = TSIDE;
									// граничная грань:
									sosedi[G][maxelm].iNODE4 = maxelm_memo + maxbound; // граничные КО нумеруются в последнюю очередь,
									sosedi[G][maxelm].bdroblenie4 = true;
									// после maxelm_memo внутренних КО.
									// Вычисление дальнего соседа
									CALC_GG(G, GG); // Вычисляет грань GG по грани G.

									sosedi[GG][maxelm].iNODE4 = NON_EXISTENT_NODE; // соседа нет,
									sosedi[GG][maxelm].bdroblenie4 = true;

#if doubleintprecision == 1
									//printf("error inum_TD == 0 sosedi T7 elm=%lld\n", maxelm);
#else
									//printf("error inum_TD == 0 sosedi T7 elm=%d\n", maxelm);
#endif
									
									//getchar();
									
									// Эта ситуация возможна и здесь надо предусмотреть её обработку.
									DEFINE_BOUNDARY_PROPERTIES_TEMP(prop, prop_b, maxelm, maxbound); // В граничный узел сносятся свойства прилегающего КО.
									
									// Sx*Sy
									doublereal dS_loc = fabs(octree1->linkT7->p3.y - octree1->linkT7->p0.y)*fabs(octree1->linkT7->p1.x - octree1->linkT7->p0.x);
									TOCHKA p_centerG;
									p_centerG.z = octree1->linkT7->p0.z;
									p_centerG.y = 0.5*(octree1->linkT7->p3.y + octree1->linkT7->p0.y);
									p_centerG.x = 0.5*(octree1->linkT7->p1.x + octree1->linkT7->p0.x);

									// G - грань по которой ставится граничное условие.
									// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
									// bound_id - номер граничного узла, начиная с нуля.
									// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
									if (octree1->b4B) {
										obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
											octree1->linkB0->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
									}
									else {
										if (octree1->linkB == NULL) {
											obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
												 NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
										}
										else {
											obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
												octree1->linkB->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
										}

									}

									maxbound++;
								}
								else {
									// сосед существует.
									if (bvisit[octree1->linkT7->inum_TD - 1]) {
										// узел уже был посещён
										// Здесь неким образом модифицировалось gran_t
										// Только внутренняя грань.
										G = TSIDE;
										// граничная грань:
										sosedi[G][maxelm].iNODE4 = octree1->linkT7->inum_TD - 1; // граничные КО нумеруются в последнюю очередь,
										sosedi[G][maxelm].bdroblenie4 = true;
										// после maxelm_memo внутренних КО.
										// Вычисление дальнего соседа
										CALC_GG(G, GG); // Вычисляет грань GG по грани G.

										// Определение дальнего соседа :
										// TODO 13 сентября 2016.
									}
									else {
										// узнать координаты центра грани и ориентацию в пространстве
										x_c = 0.5*(octree1->linkT7->p0.x + octree1->linkT7->p1.x);
										y_c = 0.5*(octree1->linkT7->p0.y + octree1->linkT7->p3.y);
										z_c = octree1->p4.z;
										iplane = XY;
										patch_maxbound(iplane, s, ls, x_c, y_c, z_c, maxbound, binc4, lsid4);
										if (binc4) {
#if doubleintprecision == 1
											printf("error internal source T7 elm=%lld\n", maxelm);
#else
											printf("error internal source T7 elm=%d\n", maxelm);
#endif
											
											//getchar();
											system("PAUSE");
										}
										else {
											// Только внутренняя грань.
											G = TSIDE;
											// граничная грань:
											sosedi[G][maxelm].iNODE4 = octree1->linkT7->inum_TD - 1; // граничные КО нумеруются в последнюю очередь,
											sosedi[G][maxelm].bdroblenie4 = true;
											// после maxelm_memo внутренних КО.
											// Вычисление дальнего соседа
											CALC_GG(G, GG); // Вычисляет грань GG по грани G.

											// Определение дальнего соседа :
											// TODO 13 сентября 2016.
										}
									}
								}
							}
						}
					}

					// сосед не существует
					// записываем идентификатор грани и
					// увеличиваем счётчик граней.
					if (!octree1->b4B) {
						if ((octree1->linkB == NULL)||((octree1->linkB != NULL)&&(octree1->linkB->inum_TD == 0))) {
							// Если мы находимся на границе расчётной области или если мы граничим с HOLLOW блоком.

							G = BSIDE;
							// граничная грань:
							INIT_SOSEDI(sosedi, maxelm, BSIDE, GG, false, maxelm_memo + maxbound,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE, false,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE); // заполнение sosedi[G][maxelm].
							DEFINE_BOUNDARY_PROPERTIES_TEMP(prop, prop_b, maxelm, maxbound); // В граничный узел сносятся свойства прилегающего КО.
						
							doublereal dSloc = CALC_SQUARE(G, maxelm, nvtx, pa); // Вычисление площади грани ячейки.
							TOCHKA p_centerG = CALC_CENTERG(G, maxelm, nvtx, pa); // Вычисление геометрического центра грани.
							
							// G - грань по которой ставится граничное условие.
							// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
							// bound_id - номер граничного узла, начиная с нуля.
							// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
							if (octree1->b4T) {
								obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
									octree1->linkT4->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, false, false, p_centerG);
							}
							else {
								if (octree1->linkT == NULL) {
									obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
										 NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, false, false, p_centerG);
								}
								else {
									obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
										octree1->linkT->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, false, false, p_centerG);
								}
							}

							maxbound++;
						}
						else {
							if (octree1->linkB != NULL) {
								// сосед существует.
								if (bvisit[octree1->linkB->inum_TD - 1]) {

									// return_current_node_number - Заглянули в соседа и вернулись из соседа в текущий узел.
									// Номер должен быть тем же самым и поэтому этот код бессмысленен.
									// Строго внутренние плоские бесконечно тонкие источники тепла внутри расчётной области больше не поддерживаются.
									// Поэтому от этого кода можно безболезненно избавиться. Поддерживаются только объёмные источники тепла.
									// 10.03.2019
									integer return_current_node_number = sosedi[TSIDE][octree1->linkB->inum_TD - 1].iNODE1;
									if (octree1->linkB->b4B) {
										integer current_node_number = octree1->inum_TD - 1;
										bool bcontinue = true;
										if (bcontinue) {
											if (octree1->linkB->linkT4 != NULL) {
												if (octree1->linkB->linkT4->inum_TD - 1 == current_node_number) {
													bcontinue = false; // iNODE1
													if (sosedi[TSIDE][octree1->linkB->inum_TD - 1].iNODE1 == current_node_number) {
														return_current_node_number = sosedi[TSIDE][octree1->linkB->inum_TD - 1].iNODE1;
													}
													else {
														printf("this can not be B T4\n");
														//getchar();
														system("PAUSE");
													}
												}
											}
										}
										if (bcontinue) {
											if (octree1->linkB->linkT5 != NULL) {
												if (octree1->linkB->linkT5->inum_TD - 1 == current_node_number) {
													bcontinue = false; // iNODE2
													if (sosedi[TSIDE][octree1->linkB->inum_TD - 1].iNODE2 == current_node_number) {
														return_current_node_number = sosedi[TSIDE][octree1->linkB->inum_TD - 1].iNODE2;
													}
													else {
														printf("this can not be B T5\n");
														//getchar();
														system("PAUSE");
													}
												}
											}
										}
										if (bcontinue) {
											if (octree1->linkB->linkT6 != NULL) {
												if (octree1->linkB->linkT6->inum_TD - 1 == current_node_number) {
													bcontinue = false; // iNODE3
													if (sosedi[TSIDE][octree1->linkB->inum_TD - 1].iNODE3 == current_node_number) {
														return_current_node_number = sosedi[TSIDE][octree1->linkB->inum_TD - 1].iNODE3;
													}
													else {
														printf("this can not be B T6\n");
														//getchar();
														system("PAUSE");
													}
												}
											}
										}
										if (bcontinue) {
											if (octree1->linkB->linkT7 != NULL) {
												if (octree1->linkB->linkT7->inum_TD - 1 == current_node_number) {
													bcontinue = false; // iNODE4
													if (sosedi[TSIDE][octree1->linkB->inum_TD - 1].iNODE4 == current_node_number) {
														return_current_node_number = sosedi[TSIDE][octree1->linkB->inum_TD - 1].iNODE4;
													}
													else {
														printf("this can not be B T7\n");
														//getchar();
														system("PAUSE");
													}
												}
											}
										}

									}

									// узел уже был посещён
									// Здесь неким образом модифицировалось gran_t
									if (return_current_node_number >= maxelm_memo) {

										printf("WARNING!!! internal source gran BSIDE 10.03.2019\n");
										system("PAUSE");

										// Это граничная грань внутреннего источника тепла.
										// Внутренний источник тепла.
										G = BSIDE;
										// граничная грань:
										INIT_SOSEDI(sosedi, maxelm, BSIDE, GG, false, sosedi[TSIDE][octree1->linkB->inum_TD - 1].iNODE1,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE, false,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE); // заполнение sosedi[G][maxelm].
										
										if (bi_fluid) {
											// FLUID
											// В граничный узел сносятся свойства прилегающего КО.
											DEFINE_BOUNDARY_PROPERTIES_TEMP(prop, prop_b, octree1->linkB->inum_TD - 1, sosedi[G][maxelm].iNODE1 - maxelm_memo);
										}
										else {
											// SOLID
											// В граничный узел сносятся свойства прилегающего КО.
											DEFINE_BOUNDARY_PROPERTIES_TEMP(prop, prop_b, maxelm, sosedi[G][maxelm].iNODE1 - maxelm_memo);
										}

										doublereal dSloc = CALC_SQUARE(G, maxelm, nvtx, pa); // Вычисление площади грани ячейки.
										TOCHKA p_centerG = CALC_CENTERG(G, maxelm, nvtx, pa); // Вычисление геометрического центра грани.
										
										// G - грань по которой ставится граничное условие.
										// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
										// bound_id - номер граничного узла, начиная с нуля.
										// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
										if (octree1->b4T) {
											obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, sosedi[TSIDE][octree1->linkB->inum_TD - 1].iNODE1 - maxelm_memo,
												octree1->linkT4->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, true, true, p_centerG);
										}
										else {
											if (octree1->linkT == NULL) {
												obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, sosedi[TSIDE][octree1->linkB->inum_TD - 1].iNODE1 - maxelm_memo,
													 NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, true, true, p_centerG);
											}
											else {
												obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, sosedi[TSIDE][octree1->linkB->inum_TD - 1].iNODE1 - maxelm_memo,
													octree1->linkT->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, true, true, p_centerG);
											}
										}
										binternalsource[sosedi[TSIDE][octree1->linkB->inum_TD - 1].iNODE1 - maxelm_memo] = true; // внутренний источник тепла на границе жидкости и твёрдого тела

									}
									else {
										// Это строго внутренний узел.
										// Внутренняя грань.
										G = BSIDE;
										// Внутренняя грань:
										sosedi[G][maxelm].iNODE1 = octree1->linkB->inum_TD - 1; // граничные КО нумеруются в последнюю очередь,
										sosedi[G][maxelm].bdroblenie4 = false;
										sosedi[G][maxelm].iNODE2 = NON_EXISTENT_NODE;
										sosedi[G][maxelm].iNODE3 = NON_EXISTENT_NODE;
										sosedi[G][maxelm].iNODE4 = NON_EXISTENT_NODE;
										// после maxelm_memo внутренних КО.
										// Вычисление дальнего соседа
										//CALC_GG(G, GG); // Вычисляет грань GG по грани G.

										// Определение дальнего соседа :
										// TODO 13 сентября 2016.
									}

								}
								else {
									bool binc;
									integer lsid;
									// узнать координаты центра грани и ориентацию в пространстве
									x_c = 0.5*(octree1->p0.x + octree1->p1.x);
									y_c = 0.5*(octree1->p0.y + octree1->p3.y);
									z_c = octree1->p0.z;
									iplane = XY;
									patch_maxbound(iplane, s, ls, x_c, y_c, z_c, maxbound,binc,lsid);
									if (binc) {
										// Внутренний источник тепла.
										// Внутренние плоские бесконечно тонкие источники тепла скорее всего больше не поддерживаются. 10,03,2019

										G = BSIDE;
										// граничная грань:
										INIT_SOSEDI(sosedi, maxelm, BSIDE, GG, false, maxelm_memo + maxbound -1,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE, false,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE); // заполнение sosedi[G][maxelm].
										
										if (bi_fluid) {
											// FLUID
											// В граничный узел сносятся свойства прилегающего КО.
											DEFINE_BOUNDARY_PROPERTIES_TEMP(prop, prop_b, octree1->linkB->inum_TD - 1, maxbound-1);
										}
										else {
											// SOLID
											// В граничный узел сносятся свойства прилегающего КО.
											DEFINE_BOUNDARY_PROPERTIES_TEMP(prop, prop_b, maxelm, maxbound-1);
										}

										doublereal dSloc = CALC_SQUARE(G, maxelm, nvtx, pa); // Вычисление площади грани ячейки.
										TOCHKA p_centerG = CALC_CENTERG(G, maxelm, nvtx, pa); // Вычисление геометрического центра грани.
										
										// G - грань по которой ставится граничное условие.
										// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
										// bound_id - номер граничного узла, начиная с нуля.
										// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
										if (octree1->b4T) {
											obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound - 1,
												octree1->linkT4->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, true, false, p_centerG);
										}
										else {
											if (octree1->linkT == NULL) {
												obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound - 1,
													 NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, true, false, p_centerG);
											}
											else {
												obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound - 1,
													octree1->linkT->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, true, false, p_centerG);
											}
										}
										binternalsource[maxbound - 1] = true; // внутренний источник тепла на границе жидкости и твёрдого тела

									}
									else {
										// Внутренняя грань.
										G = BSIDE;
										// Внутренняя грань:
										sosedi[G][maxelm].iNODE1 = octree1->linkB->inum_TD - 1; // граничные КО нумеруются в последнюю очередь,
										sosedi[G][maxelm].bdroblenie4 = false;
										sosedi[G][maxelm].iNODE2 = NON_EXISTENT_NODE;
										sosedi[G][maxelm].iNODE3 = NON_EXISTENT_NODE;
										sosedi[G][maxelm].iNODE4 = NON_EXISTENT_NODE;
										// после maxelm_memo внутренних КО.
										// Вычисление дальнего соседа
										//CALC_GG(G, GG); // Вычисляет грань GG по грани G.

										// Определение дальнего соседа :
										// TODO 13 сентября 2016.
									}
								}
							}
						}
					}
					else {
					     // Присутствуют все 4 соседа: B: 0,1,2,3.

						if (((octree1->linkB0 == NULL) && (octree1->linkB1 == NULL) && (octree1->linkB2 == NULL) && (octree1->linkB3 == NULL))||
							(((octree1->linkB0 != NULL) && (octree1->linkB1 != NULL) && (octree1->linkB2 != NULL) && (octree1->linkB3 != NULL) &&
							(octree1->linkB0->inum_TD == 0) && (octree1->linkB1->inum_TD == 0) && (octree1->linkB2->inum_TD == 0) && (octree1->linkB3->inum_TD == 0)))||
							(((octree1->linkB0 != NULL) && (octree1->linkB1 == NULL) && (octree1->linkB2 == NULL) && (octree1->linkB3 != NULL) &&
							(octree1->linkB0->inum_TD == 0) && (octree1->linkB3->inum_TD == 0)))||
								(((octree1->linkB0 != NULL) && (octree1->linkB1 != NULL) && (octree1->linkB2 == NULL) && (octree1->linkB3 == NULL) &&
								(octree1->linkB0->inum_TD == 0) && (octree1->linkB1->inum_TD == 0)))||
								(((octree1->linkB0 != NULL) && (octree1->linkB1 == NULL) && (octree1->linkB2 == NULL) && (octree1->linkB3 == NULL) &&
							(octree1->linkB0->inum_TD == 0))))
						
						{
							// 1 - B0&&B1&&B2&&B3 - на границе расчётной области (РО).
							// 2 - B0&&B1&&B2&&B3 - на гранеице HOLLOW блока.
							// 3 - B0&&B3 на границе HOLLOW блока, B1&&B2 на границе РО.
							// 4 - B0&&B1 на границе HOLLOW блока, B2&B3 на границе РО.
							// 5 - B0 на границе HOLLOW блока, B1&&B2&&B3  на границе РО.

							G = BSIDE;
							// граничная грань:
							INIT_SOSEDI(sosedi, maxelm, BSIDE, GG, false, maxelm_memo + maxbound,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE, false,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE,  NON_EXISTENT_NODE); // заполнение sosedi[G][maxelm].
							DEFINE_BOUNDARY_PROPERTIES_TEMP(prop, prop_b, maxelm, maxbound); // В граничный узел сносятся свойства прилегающего КО.

							doublereal dSloc = CALC_SQUARE(G, maxelm, nvtx, pa); // Вычисление площади грани ячейки.
							TOCHKA p_centerG = CALC_CENTERG(G, maxelm, nvtx, pa); // Вычисление геометрического центра грани.
							
							// G - грань по которой ставится граничное условие.
							// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
							// bound_id - номер граничного узла, начиная с нуля.
							// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
							if (octree1->b4T) {
								obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
									octree1->linkT4->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, false, false, p_centerG);
							}
							else {
								if (octree1->linkT == NULL) {
									obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
										 NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, false, false, p_centerG);
								}
								else {
									obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
										octree1->linkT->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, false, false, p_centerG);
								}
							}

							// Зачем дробить если граница заканчивается NULL.
							maxbound++; // Вместо четырёх один.
						}						
						else {

							bool binc1, binc2, binc3, binc4;
							integer lsid1, lsid2, lsid3, lsid4;

							if (octree1->linkB0 == NULL) {
#if doubleintprecision == 1
								printf("error NULL sosedi B0 elm=%lld\n", maxelm);
#else
								printf("error NULL sosedi B0 elm=%d\n", maxelm);
#endif
								
								//getchar();
								system("PAUSE");

								/*
								G = BSIDE;
								// граничная грань:
								sosedi[G][maxelm].iNODE1 = maxelm_memo + maxbound; // граничные КО нумеруются в последнюю очередь,
								sosedi[G][maxelm].bdroblenie4 = true;
								// после maxelm_memo внутренних КО.
								// Вычисление дальнего соседа
								CALC_GG(G,GG); // Вычисляет грань GG по грани G.

								sosedi[GG][maxelm].iNODE1 = NON_EXISTENT_NODE; // граничные КО нумеруются в последнюю очередь,
								sosedi[GG][maxelm].bdroblenie4 = true;

								#if doubleintprecision == 1
										printf("error NULL sosedi B0 elm=%lld\n", maxelm);
								#else
										printf("error NULL sosedi B0 elm=%d\n", maxelm);
								#endif
								
								getchar();

								// Эта ситуация возможна и здесь надо предусмотреть её обработку.
								DEFINE_BOUNDARY_PROPERTIES_TEMP(prop, prop_b, maxelm, maxbound); // В граничный узел сносятся свойства прилегающего КО.
								
								doublereal dS_loc = 0.0;

								// G - грань по которой ставится граничное условие.
								// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
								// bound_id - номер граничного узла, начиная с нуля.
								// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
								if (octree1->b4T) {
									obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
										octree1->linkT4->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
								}
								else {
									if (octree1->linkT == NULL) {
										obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
											 NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
									}
									else {
										obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
											octree1->linkT->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
									}
								}

								maxbound++;
								*/
							}
							else {
								if (octree1->linkB0->inum_TD == 0) {

									G = BSIDE;
									// граничная грань:
									sosedi[G][maxelm].iNODE1 = maxelm_memo + maxbound; // граничные КО нумеруются в последнюю очередь,
									sosedi[G][maxelm].bdroblenie4 = true;
									// после maxelm_memo внутренних КО.
									// Вычисление дальнего соседа
									CALC_GG(G, GG); // Вычисляет грань GG по грани G.

									sosedi[GG][maxelm].iNODE1 = NON_EXISTENT_NODE; // граничные КО нумеруются в последнюю очередь,
									sosedi[GG][maxelm].bdroblenie4 = true;

#if doubleintprecision == 1
									//printf("error inum_TD == 0 sosedi B0 elm=%lld\n", maxelm);
#else
									//printf("error inum_TD == 0 sosedi B0 elm=%d\n", maxelm);
#endif
									
									//getchar();

									// Эта ситуация возможна и здесь надо предусмотреть её обработку.
									DEFINE_BOUNDARY_PROPERTIES_TEMP(prop, prop_b, maxelm, maxbound); // В граничный узел сносятся свойства прилегающего КО.
									
									// Sx*Sy
									doublereal dS_loc = fabs(octree1->linkB0->p7.y - octree1->linkB0->p4.y)*fabs(octree1->linkB0->p6.x - octree1->linkB0->p7.x);
									TOCHKA p_centerG;
									p_centerG.z = octree1->linkB0->p4.z;
									p_centerG.y = 0.5*(octree1->linkB0->p7.y + octree1->linkB0->p4.y);
									p_centerG.x = 0.5*(octree1->linkB0->p6.x + octree1->linkB0->p7.x);


									// G - грань по которой ставится граничное условие.
									// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
									// bound_id - номер граничного узла, начиная с нуля.
									// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
									if (octree1->b4T) {
										obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
											octree1->linkT4->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
									}
									else {
										if (octree1->linkT == NULL) {
											obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
												 NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
										}
										else {
											obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
												octree1->linkT->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
										}
									}

									maxbound++;
								}
								else {
									// сосед существует.
									if (bvisit[octree1->linkB0->inum_TD - 1]) {
										// узел уже был посещён
										// Здесь неким образом модифицировалось gran_t
										// Только внутренняя грань.
										G = BSIDE;
										// граничная грань:
										sosedi[G][maxelm].iNODE1 = octree1->linkB0->inum_TD - 1; // граничные КО нумеруются в последнюю очередь,
										sosedi[G][maxelm].bdroblenie4 = true;
										// после maxelm_memo внутренних КО.
										// Вычисление дальнего соседа
										CALC_GG(G, GG); // Вычисляет грань GG по грани G.

										// Определение дальнего соседа :
										// TODO 13 сентября 2016.
									}
									else {
										// узнать координаты центра грани и ориентацию в пространстве
										x_c = 0.5*(octree1->linkB0->p4.x + octree1->linkB0->p5.x);
										y_c = 0.5*(octree1->linkB0->p4.y + octree1->linkB0->p7.y);
										z_c = octree1->p0.z;
										iplane = XY;
										patch_maxbound(iplane, s, ls, x_c, y_c, z_c, maxbound, binc1, lsid1);
										if (binc1) {
#if doubleintprecision == 1
											printf("error internal source B0 elm=%lld\n", maxelm);
#else
											printf("error internal source B0 elm=%d\n", maxelm);
#endif
											
											//getchar();
											system("PAUSE");
										}
										else {
											// Только внутренняя грань.
											G = BSIDE;
											// граничная грань:
											sosedi[G][maxelm].iNODE1 = octree1->linkB0->inum_TD - 1; // граничные КО нумеруются в последнюю очередь,
											sosedi[G][maxelm].bdroblenie4 = true;
											// после maxelm_memo внутренних КО.
											// Вычисление дальнего соседа
											CALC_GG(G, GG); // Вычисляет грань GG по грани G.

											// Определение дальнего соседа :
											// TODO 13 сентября 2016.
										}
									}
								}
							}

							if (octree1->linkB1 == NULL) {

								// 21 сентября 2016 (вырождение ячейки).

								// Заглушка.
								G = BSIDE;
								// граничная грань:
								sosedi[G][maxelm].iNODE2 = NON_EXISTENT_NODE; // граничные КО нумеруются в последнюю очередь,
								sosedi[G][maxelm].bdroblenie4 = true;
								// после maxelm_memo внутренних КО.
								// Вычисление дальнего соседа
								CALC_GG(G, GG); // Вычисляет грань GG по грани G.

								sosedi[GG][maxelm].iNODE2 = NON_EXISTENT_NODE; // сосед отсутствует,
								sosedi[GG][maxelm].bdroblenie4 = true;

								/*
								G = BSIDE;
								// граничная грань:
								sosedi[G][maxelm].iNODE2 = maxelm_memo + maxbound; // граничные КО нумеруются в последнюю очередь,
								sosedi[G][maxelm].bdroblenie4 = true;
								// после maxelm_memo внутренних КО.
								// Вычисление дальнего соседа
								CALC_GG(G,GG); // Вычисляет грань GG по грани G.

								sosedi[GG][maxelm].iNODE2 = NON_EXISTENT_NODE; // сосед отсутствует,
								sosedi[GG][maxelm].bdroblenie4 = true;

								#if doubleintprecision == 1
										printf("error NULL sosedi B1 elm=%lld\n", maxelm);
								#else
										printf("error NULL sosedi B1 elm=%d\n", maxelm);
								#endif
								
								getchar();

								// Эта ситуация возможна и здесь надо предусмотреть её обработку.
								DEFINE_BOUNDARY_PROPERTIES_TEMP(prop, prop_b, maxelm, maxbound); // В граничный узел сносятся свойства прилегающего КО.
								
								doublereal dS_loc = 0.0;

								// G - грань по которой ставится граничное условие.
								// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
								// bound_id - номер граничного узла, начиная с нуля.
								// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
								if (octree1->b4T) {
									obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
										octree1->linkT4->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
								}
								else {
									if (octree1->linkT == NULL) {
										obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
											 NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
									}
									else {
										obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
											octree1->linkT->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
									}
								}

								maxbound++;
								*/
							}
							else {
								if (octree1->linkB1->inum_TD == 0) {
									
									G = BSIDE;
									// граничная грань:
									sosedi[G][maxelm].iNODE2 = maxelm_memo + maxbound; // граничные КО нумеруются в последнюю очередь,
									sosedi[G][maxelm].bdroblenie4 = true;
									// после maxelm_memo внутренних КО.
									// Вычисление дальнего соседа
									CALC_GG(G, GG); // Вычисляет грань GG по грани G.

									sosedi[GG][maxelm].iNODE2 = NON_EXISTENT_NODE; // сосед отсутствует,
									sosedi[GG][maxelm].bdroblenie4 = true;

#if doubleintprecision == 1
									//printf("error inum_TD == 0 sosedi B1 elm=%lld\n", maxelm);
#else
									//printf("error inum_TD == 0 sosedi B1 elm=%d\n", maxelm);
#endif
									
									//getchar();

									// Эта ситуация возможна и здесь надо предусмотреть её обработку.
									DEFINE_BOUNDARY_PROPERTIES_TEMP(prop, prop_b, maxelm, maxbound); // В граничный узел сносятся свойства прилегающего КО.
									
									// Sx*Sy
									doublereal dS_loc = fabs(octree1->linkB1->p7.y - octree1->linkB1->p4.y)*fabs(octree1->linkB1->p6.x - octree1->linkB1->p7.x);
									TOCHKA p_centerG;
									p_centerG.z = octree1->linkB1->p4.z;
									p_centerG.y = 0.5*(octree1->linkB1->p7.y + octree1->linkB1->p4.y);
									p_centerG.x = 0.5*(octree1->linkB1->p6.x + octree1->linkB1->p7.x);

									// G - грань по которой ставится граничное условие.
									// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
									// bound_id - номер граничного узла, начиная с нуля.
									// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
									if (octree1->b4T) {
										obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
											octree1->linkT4->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
									}
									else {
										if (octree1->linkT == NULL) {
											obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
												 NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
										}
										else {
											obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
												octree1->linkT->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
										}
									}

									maxbound++;
								}
								else {
									// сосед существует.
									if (bvisit[octree1->linkB1->inum_TD - 1]) {
										// узел уже был посещён
										// Здесь неким образом модифицировалось gran_t
										// Только внутренняя грань.
										G = BSIDE;
										// граничная грань:
										sosedi[G][maxelm].iNODE2 = octree1->linkB1->inum_TD - 1; // граничные КО нумеруются в последнюю очередь,
										sosedi[G][maxelm].bdroblenie4 = true;
										// после maxelm_memo внутренних КО.
										// Вычисление дальнего соседа
										CALC_GG(G, GG); // Вычисляет грань GG по грани G.

										// Определение дальнего соседа :
										// TODO 13 сентября 2016.
									}
									else {
										// узнать координаты центра грани и ориентацию в пространстве
										x_c = 0.5*(octree1->linkB1->p4.x + octree1->linkB1->p5.x);
										y_c = 0.5*(octree1->linkB1->p4.y + octree1->linkB1->p7.y);
										z_c = octree1->p0.z;
										iplane = XY;
										patch_maxbound(iplane, s, ls, x_c, y_c, z_c, maxbound, binc2, lsid2);
										if (binc2) {
#if doubleintprecision == 1
											printf("error internal source B1 elm=%lld\n", maxelm);
#else
											printf("error internal source B1 elm=%d\n", maxelm);
#endif
											
											//getchar();
											system("PAUSE");
										}
										else {
											// Только внутренняя грань.
											G = BSIDE;
											// граничная грань:
											sosedi[G][maxelm].iNODE2 = octree1->linkB1->inum_TD - 1; // граничные КО нумеруются в последнюю очередь,
											sosedi[G][maxelm].bdroblenie4 = true;
											// после maxelm_memo внутренних КО.
											// Вычисление дальнего соседа
											CALC_GG(G, GG); // Вычисляет грань GG по грани G.

											// Определение дальнего соседа :
											// TODO 13 сентября 2016.
										}
									}
								}
							}

							if (octree1->linkB2 == NULL) {

								// 21 сентября 2016 (вырождение ячейки).
								
								// Заглушка.
								G = BSIDE;
								// граничная грань:
								sosedi[G][maxelm].iNODE3 = NON_EXISTENT_NODE; // граничные КО нумеруются в последнюю очередь,
								sosedi[G][maxelm].bdroblenie4 = true;
								// после maxelm_memo внутренних КО.
								// Вычисление дальнего соседа
								CALC_GG(G, GG); // Вычисляет грань GG по грани G.

								sosedi[GG][maxelm].iNODE3 = NON_EXISTENT_NODE; // граничные КО нумеруются в последнюю очередь,
								sosedi[GG][maxelm].bdroblenie4 = true;

								/*
								G = BSIDE;
								// граничная грань:
								sosedi[G][maxelm].iNODE3 = maxelm_memo + maxbound; // граничные КО нумеруются в последнюю очередь,
								sosedi[G][maxelm].bdroblenie4 = true;
								// после maxelm_memo внутренних КО.
								// Вычисление дальнего соседа
								CALC_GG(G,GG); // Вычисляет грань GG по грани G.
								
								sosedi[GG][maxelm].iNODE3 = NON_EXISTENT_NODE; // граничные КО нумеруются в последнюю очередь,
								sosedi[GG][maxelm].bdroblenie4 = true;

								#if doubleintprecision == 1
										printf("error NULL sosedi B2 elm=%lld\n", maxelm);
								#else
										printf("error NULL sosedi B2 elm=%d\n", maxelm);
								#endif
								
								getchar();

								// Эта ситуация возможна и здесь надо предусмотреть её обработку.
								DEFINE_BOUNDARY_PROPERTIES_TEMP(prop, prop_b, maxelm, maxbound); // В граничный узел сносятся свойства прилегающего КО.
								
								doublereal dS_loc = 0.0;

								// G - грань по которой ставится граничное условие.
								// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
								// bound_id - номер граничного узла, начиная с нуля.
								// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
								if (octree1->b4T) {
									obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
										octree1->linkT4->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
								}
								else {
									if (octree1->linkT == NULL) {
										obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
											 NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
									}
									else {
										obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
											octree1->linkT->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
									}
								}


								maxbound++;
								*/
							}
							else {
								if (octree1->linkB2->inum_TD == 0) {
									
									G = BSIDE;
									// граничная грань:
									sosedi[G][maxelm].iNODE3 = maxelm_memo + maxbound; // граничные КО нумеруются в последнюю очередь,
									sosedi[G][maxelm].bdroblenie4 = true;
									// после maxelm_memo внутренних КО.
									// Вычисление дальнего соседа
									CALC_GG(G, GG); // Вычисляет грань GG по грани G.

									sosedi[GG][maxelm].iNODE3 = NON_EXISTENT_NODE; // граничные КО нумеруются в последнюю очередь,
									sosedi[GG][maxelm].bdroblenie4 = true;
									
#if doubleintprecision == 1
									//printf("error inum_TD == 0 sosedi B2 elm=%lld\n", maxelm);
#else
									//printf("error inum_TD == 0 sosedi B2 elm=%d\n", maxelm);
#endif
									
									//getchar();


									// Эта ситуация возможна и здесь надо предусмотреть её обработку.
									DEFINE_BOUNDARY_PROPERTIES_TEMP(prop, prop_b, maxelm, maxbound); // В граничный узел сносятся свойства прилегающего КО.
									
									// Sx*Sy
									doublereal dS_loc = fabs(octree1->linkB2->p7.y - octree1->linkB2->p4.y)*fabs(octree1->linkB2->p6.x - octree1->linkB2->p7.x);
									TOCHKA p_centerG;
									p_centerG.z = octree1->linkB2->p4.z;
									p_centerG.y = 0.5*(octree1->linkB2->p7.y + octree1->linkB2->p4.y);
									p_centerG.x = 0.5*(octree1->linkB2->p6.x + octree1->linkB2->p7.x);

									// G - грань по которой ставится граничное условие.
									// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
									// bound_id - номер граничного узла, начиная с нуля.
									// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
									if (octree1->b4T) {
										obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
											octree1->linkT4->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
									}
									else {
										if (octree1->linkT == NULL) {
											obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
												 NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
										}
										else {
											obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
												octree1->linkT->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
										}
									}

									maxbound++;
								}
								else {
									// сосед существует.
									if (bvisit[octree1->linkB2->inum_TD - 1]) {
										// узел уже был посещён
										// Здесь неким образом модифицировалось gran_t
										// Только внутренняя грань.
										G = BSIDE;
										// граничная грань:
										sosedi[G][maxelm].iNODE3 = octree1->linkB2->inum_TD - 1; // граничные КО нумеруются в последнюю очередь,
										sosedi[G][maxelm].bdroblenie4 = true;
										// после maxelm_memo внутренних КО.
										// Вычисление дальнего соседа
										CALC_GG(G, GG); // Вычисляет грань GG по грани G.

										// Определение дальнего соседа :
										// TODO 13 сентября 2016.
									}
									else {
										// узнать координаты центра грани и ориентацию в пространстве
										x_c = 0.5*(octree1->linkB2->p4.x + octree1->linkB2->p5.x);
										y_c = 0.5*(octree1->linkB2->p4.y + octree1->linkB2->p7.y);
										z_c = octree1->p0.z;
										iplane = XY;
										patch_maxbound(iplane, s, ls, x_c, y_c, z_c, maxbound, binc3, lsid3);
										if (binc3) {
#if doubleintprecision == 1
											printf("error internal source B2 elm=%lld\n", maxelm);
#else
											printf("error internal source B2 elm=%d\n", maxelm);
#endif
											
											//getchar();
											system("PAUSE");
										}
										else {
											// Только внутренняя грань.
											G = BSIDE;
											// граничная грань:
											sosedi[G][maxelm].iNODE3 = octree1->linkB2->inum_TD - 1; // граничные КО нумеруются в последнюю очередь,
											sosedi[G][maxelm].bdroblenie4 = true;
											// после maxelm_memo внутренних КО.
											// Вычисление дальнего соседа
											CALC_GG(G, GG); // Вычисляет грань GG по грани G.

											// Определение дальнего соседа :
											// TODO 13 сентября 2016.
										}
									}
								}
							}

							if (octree1->linkB3 == NULL) {

								// 21 сентября 2016 (вырождение ячейки).

								// Заглушка.
								G = BSIDE;
								// граничная грань:
								sosedi[G][maxelm].iNODE4 = NON_EXISTENT_NODE; // граничные КО нумеруются в последнюю очередь,
								sosedi[G][maxelm].bdroblenie4 = true;
								// после maxelm_memo внутренних КО.
								// Вычисление дальнего соседа
								CALC_GG(G, GG); // Вычисляет грань GG по грани G.

								sosedi[GG][maxelm].iNODE4 = NON_EXISTENT_NODE; // сосед отсутствует,
								sosedi[GG][maxelm].bdroblenie4 = true;

								/*
								G = BSIDE;
								// граничная грань:
								sosedi[G][maxelm].iNODE4 = maxelm_memo + maxbound; // граничные КО нумеруются в последнюю очередь,
								sosedi[G][maxelm].bdroblenie4 = true;
								// после maxelm_memo внутренних КО.
								// Вычисление дальнего соседа
								CALC_GG(G,GG); // Вычисляет грань GG по грани G.

								sosedi[GG][maxelm].iNODE4 = NON_EXISTENT_NODE; // сосед отсутствует,
								sosedi[GG][maxelm].bdroblenie4 = true;
								#if doubleintprecision == 1
										printf("error NULL sosedi B3 elm=%lld\n", maxelm);
								#else
										printf("error NULL sosedi B3 elm=%d\n", maxelm);
								#endif
								
								getchar();

								// Эта ситуация возможна и здесь надо предусмотреть её обработку.
								DEFINE_BOUNDARY_PROPERTIES_TEMP(prop, prop_b, maxelm, maxbound); // В граничный узел сносятся свойства прилегающего КО.
								
								doublereal dS_loc = 0.0;

								// G - грань по которой ставится граничное условие.
								// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
								// bound_id - номер граничного узла, начиная с нуля.
								// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
								if (octree1->b4T) {
									obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
										octree1->linkT4->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
								}
								else {
									if (octree1->linkT == NULL) {
										obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
											 NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
									}
									else {
										obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
											octree1->linkT->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
									}
								}

								maxbound++;
								*/
							}
							else {
								if (octree1->linkB3->inum_TD == 0) {
									

									G = BSIDE;
									// граничная грань:
									sosedi[G][maxelm].iNODE4 = maxelm_memo + maxbound; // граничные КО нумеруются в последнюю очередь,
									sosedi[G][maxelm].bdroblenie4 = true;
									// после maxelm_memo внутренних КО.
									// Вычисление дальнего соседа
									CALC_GG(G, GG); // Вычисляет грань GG по грани G.

									sosedi[GG][maxelm].iNODE4 = NON_EXISTENT_NODE; // сосед отсутствует,
									sosedi[GG][maxelm].bdroblenie4 = true;

#if doubleintprecision == 1
									//printf("error inum_TD == 0 sosedi B3 elm=%lld\n", maxelm);
#else
									//printf("error inum_TD == 0 sosedi B3 elm=%d\n", maxelm);
#endif
									
									//getchar();

									// Эта ситуация возможна и здесь надо предусмотреть её обработку.
									DEFINE_BOUNDARY_PROPERTIES_TEMP(prop, prop_b, maxelm, maxbound); // В граничный узел сносятся свойства прилегающего КО.
									
									// Sx*Sy
									doublereal dS_loc = fabs(octree1->linkB3->p7.y - octree1->linkB3->p4.y)*fabs(octree1->linkB3->p6.x - octree1->linkB3->p7.x);
									TOCHKA p_centerG;
									p_centerG.z = octree1->linkB3->p4.z;
									p_centerG.y = 0.5*(octree1->linkB3->p7.y + octree1->linkB3->p4.y);
									p_centerG.x = 0.5*(octree1->linkB3->p6.x + octree1->linkB3->p7.x);

									// G - грань по которой ставится граничное условие.
									// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
									// bound_id - номер граничного узла, начиная с нуля.
									// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
									if (octree1->b4T) {
										obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
											octree1->linkT4->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
									}
									else {
										if (octree1->linkT == NULL) {
											obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
												 NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
										}
										else {
											obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
												octree1->linkT->inum_TD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
										}
									}

									maxbound++;
								}
								else {
									// сосед существует.
									if (bvisit[octree1->linkB3->inum_TD - 1]) {
										// узел уже был посещён
										// Здесь неким образом модифицировалось gran_t
										// Только внутренняя грань.
										G = BSIDE;
										// граничная грань:
										sosedi[G][maxelm].iNODE4 = octree1->linkB3->inum_TD - 1; // граничные КО нумеруются в последнюю очередь,
										sosedi[G][maxelm].bdroblenie4 = true;
										// после maxelm_memo внутренних КО.
										// Вычисление дальнего соседа
										CALC_GG(G, GG); // Вычисляет грань GG по грани G.

										// Определение дальнего соседа :
										// TODO 13 сентября 2016.
									}
									else {
										// узнать координаты центра грани и ориентацию в пространстве
										x_c = 0.5*(octree1->linkB3->p4.x + octree1->linkB3->p5.x);
										y_c = 0.5*(octree1->linkB3->p4.y + octree1->linkB3->p7.y);
										z_c = octree1->p0.z;
										iplane = XY;
										patch_maxbound(iplane, s, ls, x_c, y_c, z_c, maxbound, binc4, lsid4);
										if (binc4) {
#if doubleintprecision == 1
											printf("error internal source B3 elm=%lld\n", maxelm);
#else
											printf("error internal source B3 elm=%d\n", maxelm);
#endif
											
											//getchar();
											system("PAUSE");
										}
										else {
											// Только внутренняя грань.
											G = BSIDE;
											// граничная грань:
											sosedi[G][maxelm].iNODE4 = octree1->linkB3->inum_TD - 1; // граничные КО нумеруются в последнюю очередь,
											sosedi[G][maxelm].bdroblenie4 = true;
											// после maxelm_memo внутренних КО.
											// Вычисление дальнего соседа
											CALC_GG(G, GG); // Вычисляет грань GG по грани G.

											// Определение дальнего соседа :
											// TODO 13 сентября 2016.
										}
									}
								}
							}
						}
					}


					bvisit[maxelm] = true; // Узел был посещён.
					maxelm++;
				}
				//octree1->inum_TD = 0; // По умолчанию не принадлежит расчётной области.
				octree1 = NULL;
				my_ALICE_STACK[top_ALICE_STACK - 1].link = NULL;
				top_ALICE_STACK--;
			}
			else {
				// продолжаем добираться до листьев.
				STACK_ALICE buf1 = my_ALICE_STACK[top_ALICE_STACK - 1];
				STACK_ALICE* buf = &buf1;
				top_ALICE_STACK--;
				if (buf->link->link0 != NULL) {
					my_ALICE_STACK[top_ALICE_STACK].link = (buf->link->link0);
					my_ALICE_STACK[top_ALICE_STACK].maxx = buf->link->link0->maxx;
					my_ALICE_STACK[top_ALICE_STACK].minx = buf->link->link0->minx;
					my_ALICE_STACK[top_ALICE_STACK].maxy = buf->link->link0->maxy;
					my_ALICE_STACK[top_ALICE_STACK].miny = buf->link->link0->miny;
					my_ALICE_STACK[top_ALICE_STACK].maxz = buf->link->link0->maxz;
					my_ALICE_STACK[top_ALICE_STACK].minz = buf->link->link0->minz;
					top_ALICE_STACK++;
				}
				if (buf->link->link1 != NULL) {
					my_ALICE_STACK[top_ALICE_STACK].link = (buf->link->link1);
					my_ALICE_STACK[top_ALICE_STACK].maxx = buf->link->link1->maxx;
					my_ALICE_STACK[top_ALICE_STACK].minx = buf->link->link1->minx;
					my_ALICE_STACK[top_ALICE_STACK].maxy = buf->link->link1->maxy;
					my_ALICE_STACK[top_ALICE_STACK].miny = buf->link->link1->miny;
					my_ALICE_STACK[top_ALICE_STACK].maxz = buf->link->link1->maxz;
					my_ALICE_STACK[top_ALICE_STACK].minz = buf->link->link1->minz;
					top_ALICE_STACK++;
				}
				if (buf->link->link2 != NULL) {
					my_ALICE_STACK[top_ALICE_STACK].link = (buf->link->link2);
					my_ALICE_STACK[top_ALICE_STACK].maxx = buf->link->link2->maxx;
					my_ALICE_STACK[top_ALICE_STACK].minx = buf->link->link2->minx;
					my_ALICE_STACK[top_ALICE_STACK].maxy = buf->link->link2->maxy;
					my_ALICE_STACK[top_ALICE_STACK].miny = buf->link->link2->miny;
					my_ALICE_STACK[top_ALICE_STACK].maxz = buf->link->link2->maxz;
					my_ALICE_STACK[top_ALICE_STACK].minz = buf->link->link2->minz;
					top_ALICE_STACK++;
				}
				if (buf->link->link3 != NULL) {
					my_ALICE_STACK[top_ALICE_STACK].link = (buf->link->link3);
					my_ALICE_STACK[top_ALICE_STACK].maxx = buf->link->link3->maxx;
					my_ALICE_STACK[top_ALICE_STACK].minx = buf->link->link3->minx;
					my_ALICE_STACK[top_ALICE_STACK].maxy = buf->link->link3->maxy;
					my_ALICE_STACK[top_ALICE_STACK].miny = buf->link->link3->miny;
					my_ALICE_STACK[top_ALICE_STACK].maxz = buf->link->link3->maxz;
					my_ALICE_STACK[top_ALICE_STACK].minz = buf->link->link3->minz;
					top_ALICE_STACK++;
				}
				if (buf->link->link4 != NULL) {
					my_ALICE_STACK[top_ALICE_STACK].link = (buf->link->link4);
					my_ALICE_STACK[top_ALICE_STACK].maxx = buf->link->link4->maxx;
					my_ALICE_STACK[top_ALICE_STACK].minx = buf->link->link4->minx;
					my_ALICE_STACK[top_ALICE_STACK].maxy = buf->link->link4->maxy;
					my_ALICE_STACK[top_ALICE_STACK].miny = buf->link->link4->miny;
					my_ALICE_STACK[top_ALICE_STACK].maxz = buf->link->link4->maxz;
					my_ALICE_STACK[top_ALICE_STACK].minz = buf->link->link4->minz;
					top_ALICE_STACK++;
				}
				if (buf->link->link5 != NULL) {
					my_ALICE_STACK[top_ALICE_STACK].link = (buf->link->link5);
					my_ALICE_STACK[top_ALICE_STACK].maxx = buf->link->link5->maxx;
					my_ALICE_STACK[top_ALICE_STACK].minx = buf->link->link5->minx;
					my_ALICE_STACK[top_ALICE_STACK].maxy = buf->link->link5->maxy;
					my_ALICE_STACK[top_ALICE_STACK].miny = buf->link->link5->miny;
					my_ALICE_STACK[top_ALICE_STACK].maxz = buf->link->link5->maxz;
					my_ALICE_STACK[top_ALICE_STACK].minz = buf->link->link5->minz;
					top_ALICE_STACK++;
				}
				if (buf->link->link6 != NULL) {
					my_ALICE_STACK[top_ALICE_STACK].link = (buf->link->link6);
					my_ALICE_STACK[top_ALICE_STACK].maxx = buf->link->link6->maxx;
					my_ALICE_STACK[top_ALICE_STACK].minx = buf->link->link6->minx;
					my_ALICE_STACK[top_ALICE_STACK].maxy = buf->link->link6->maxy;
					my_ALICE_STACK[top_ALICE_STACK].miny = buf->link->link6->miny;
					my_ALICE_STACK[top_ALICE_STACK].maxz = buf->link->link6->maxz;
					my_ALICE_STACK[top_ALICE_STACK].minz = buf->link->link6->minz;
					top_ALICE_STACK++;
				}
				if (buf->link->link7 != NULL) {
					my_ALICE_STACK[top_ALICE_STACK].link = (buf->link->link7);
					my_ALICE_STACK[top_ALICE_STACK].maxx = buf->link->link7->maxx;
					my_ALICE_STACK[top_ALICE_STACK].minx = buf->link->link7->minx;
					my_ALICE_STACK[top_ALICE_STACK].maxy = buf->link->link7->maxy;
					my_ALICE_STACK[top_ALICE_STACK].miny = buf->link->link7->miny;
					my_ALICE_STACK[top_ALICE_STACK].maxz = buf->link->link7->maxz;
					my_ALICE_STACK[top_ALICE_STACK].minz = buf->link->link7->minz;
					top_ALICE_STACK++;
				}
			}
		}
		//}
		//getchar();
	}

	//if (bvisit != NULL) {
		// оператор delete может быть вызван повторно в том числе и к нулевому указателю.
		delete[] bvisit;
		bvisit = NULL;
	//}
} // constr_sosedi_prop_b_alice


// Вычисление соседей для каждого внутреннего узла. 
// Причём множество соседей выбирается и среди внутренних КО и 
// среди граничных КО.
// Создана на базе calculate_max_bound_temp.
// Самое начало реализации. 10.сентября.2016 14_20. 
//  Середина разработки 21 сентября 2016.
// Окончил писать и жду запуска для тестирования 22 сетября 2016.
// Написано очень плохо, много дублирующегося кода, за большим
// объёмом кода не прослеживается логика исполнения данной функции.
// Начало рефакторинга 9.марта.2019. (сократил в этом модуле в двух этих 
// функциях 6945 строк кода с сохранением работоспособности.
//  Был объём модуля 29601 строк стало 22656 строк.
void constr_sosedi_prop_b_flow_alice(octTree* &oc, BOUND* &sosedb,
	ALICE_PARTITION** &sosedi, doublereal **prop, doublereal** &prop_b,
	integer maxbound_out, integer maxelm_memo, integer maxelm_memo1,
	BLOCK* b, integer lb, SOURCE* s, integer ls, WALL* w, integer lw,
	integer* &whot_is_block,
	integer** nvtx, TOCHKA* &pa, TPROP* matlist) {

	integer G, GG; // грань ячейки дискретизации.

	// Выделение оперативной памяти.
	sosedi = NULL;
	sosedi = new ALICE_PARTITION*[12];
	if (sosedi == NULL) {
		// недостаточно памяти на данном оборудовании.
		printf("Problem : not enough memory on your equipment for sosedi constr struct_alice...\n");
		printf("Please any key to exit...\n");
		//getchar();
		system("pause");
		exit(1);
	}
	for (integer i = 0; i<12; i++) sosedi[i] = NULL;
	for (integer i = 0; i<12; i++) {
		sosedi[i] = new ALICE_PARTITION[maxelm_memo + 2];
		if (sosedi[i] == NULL) {
			// недостаточно памяти на данном оборудовании.
#if doubleintprecision == 1
			printf("Problem : not enough memory on your equipment for sosedi[%lld] constr struct_alice...\n", i);
#else
			printf("Problem : not enough memory on your equipment for sosedi[%d] constr struct_alice...\n", i);
#endif
			
			printf("Please any key to exit...\n");
			//getchar();
			system("pause");
			exit(1);
		}
	}
	// Инициализация. 09.03.2019
	for (integer i = 0; i < maxelm_memo + 2; i++) {
		integer GGloc;
		INIT_SOSEDI(sosedi, i, ESIDE, GGloc, false, NON_EXISTENT_NODE, NON_EXISTENT_NODE, NON_EXISTENT_NODE, NON_EXISTENT_NODE, false, NON_EXISTENT_NODE, NON_EXISTENT_NODE, NON_EXISTENT_NODE, NON_EXISTENT_NODE); // заполнение sosedi[G][maxelm].
		INIT_SOSEDI(sosedi, i, WSIDE, GGloc, false, NON_EXISTENT_NODE, NON_EXISTENT_NODE, NON_EXISTENT_NODE, NON_EXISTENT_NODE, false, NON_EXISTENT_NODE, NON_EXISTENT_NODE, NON_EXISTENT_NODE, NON_EXISTENT_NODE); // заполнение sosedi[G][maxelm].
		INIT_SOSEDI(sosedi, i, NSIDE, GGloc, false, NON_EXISTENT_NODE, NON_EXISTENT_NODE, NON_EXISTENT_NODE, NON_EXISTENT_NODE, false, NON_EXISTENT_NODE, NON_EXISTENT_NODE, NON_EXISTENT_NODE, NON_EXISTENT_NODE); // заполнение sosedi[G][maxelm].
		INIT_SOSEDI(sosedi, i, WSIDE, GGloc, false, NON_EXISTENT_NODE, NON_EXISTENT_NODE, NON_EXISTENT_NODE, NON_EXISTENT_NODE, false, NON_EXISTENT_NODE, NON_EXISTENT_NODE, NON_EXISTENT_NODE, NON_EXISTENT_NODE); // заполнение sosedi[G][maxelm].
		INIT_SOSEDI(sosedi, i, TSIDE, GGloc, false, NON_EXISTENT_NODE, NON_EXISTENT_NODE, NON_EXISTENT_NODE, NON_EXISTENT_NODE, false, NON_EXISTENT_NODE, NON_EXISTENT_NODE, NON_EXISTENT_NODE, NON_EXISTENT_NODE); // заполнение sosedi[G][maxelm].
		INIT_SOSEDI(sosedi, i, BSIDE, GGloc, false, NON_EXISTENT_NODE, NON_EXISTENT_NODE, NON_EXISTENT_NODE, NON_EXISTENT_NODE, false, NON_EXISTENT_NODE, NON_EXISTENT_NODE, NON_EXISTENT_NODE, NON_EXISTENT_NODE); // заполнение sosedi[G][maxelm].
	}


	// Выделение оперативной памяти:
	prop_b = NULL;
	prop_b = new doublereal*[3];
	if (prop_b == NULL) {
		// недостаточно памяти на данном оборудовании.
		printf("Problem : not enough memory on your equipment for prop_b constr struct_alice...\n");
		printf("Please any key to exit...\n");
		//getchar();
		system("pause");
		exit(1);
	}

	for (integer i = 0; i<3; i++) prop_b[i] = NULL;
	for (integer i = 0; i<3; i++) {
		prop_b[i] = new doublereal[maxbound_out];
		if (prop_b[i] == NULL) {
			// недостаточно памяти на данном оборудовании.
#if doubleintprecision == 1
			printf("Problem : not enough memory on your equipment for prop_b[%lld] constr struct_alice...\n", i);
#else
			printf("Problem : not enough memory on your equipment for prop_b[%d] constr struct_alice...\n", i);
#endif
			
			printf("Please any key to exit...\n");
			//getchar();
			system("pause");
			exit(1);
		}
	}
	integer iblock_FLUID_propb = NON_EXISTENT_NODE;// По умолчанию у нас только SOLID и HOLLOW блоки в проекте.
	// Ищем идентификатор FLUID блока в проекте.
	for (integer ib_scan = 0; ib_scan < lb; ib_scan++) {
		if (b[ib_scan].itype == FLUID) {
			// FLUID блок проекта найден.
			iblock_FLUID_propb = ib_scan;
			break;
		}
	}
	// инициализация.
	for (integer i = 0; i < maxbound_out; i++) {
		// Присваиваем свойства воздуха.
		if (iblock_FLUID_propb == -1) {
			prop_b[RHO][i] = 1.1614;
			prop_b[MU][i] = 1.7894e-5;
			prop_b[BETA_T][i] = 0.003331;
		}
		else {
			// iblock_FLUID_propb - хранит уникальный номер найденного в проекте FLUID блока.
			prop_b[RHO][i] = matlist[b[iblock_FLUID_propb].imatid].rho; 
			prop_b[MU][i] = matlist[b[iblock_FLUID_propb].imatid].mu;
			prop_b[BETA_T][i] = matlist[b[iblock_FLUID_propb].imatid].beta_t;
		}
	}

	// Алгоритм:
	// Выделение оперативной памяти:
	// gran[0..5][0..maxbound-1]
	sosedb = NULL;
	sosedb = new BOUND[maxbound_out];
	// оператор new не требует проверки на null.
	//if (sosedb == NULL) {
		// недостаточно памяти на данном оборудовании.
		//printf("Problem : not enough memory on your equipment for sosedb constr struct_alice...\n");
		//printf("Please any key to exit...\n");
		//getchar();
		//system("pause");
		//exit(1);
	//}


	integer maxbound = 0; // инициализация.
	integer maxelm = 0;
	bool *bvisit = NULL;
	bvisit = new bool[maxelm_memo + 2];
	// оператор new не требует проверки на null.
	//if (bvisit == NULL) {
		// недостаточно памяти на данном оборудовании.
		//printf("Problem : not enough memory on your equipment for bvisit constr struct_alice...\n");
		//printf("Please any key to exit...\n");
		//getchar();
		//system("pause");
		//exit(1);
	//}
	for (integer j = 0; j<maxelm_memo + 2; j++) bvisit[j] = false; // признак посещения узла.
	doublereal x_c=0.0, y_c=0.0, z_c=0.0;
	integer iplane=XY;

	top_ALICE_STACK = 0;
	if (oc->link0 != NULL) {
		my_ALICE_STACK[top_ALICE_STACK].link = (oc->link0);
		my_ALICE_STACK[top_ALICE_STACK].maxx = oc->link0->maxx;
		my_ALICE_STACK[top_ALICE_STACK].minx = oc->link0->minx;
		my_ALICE_STACK[top_ALICE_STACK].maxy = oc->link0->maxy;
		my_ALICE_STACK[top_ALICE_STACK].miny = oc->link0->miny;
		my_ALICE_STACK[top_ALICE_STACK].maxz = oc->link0->maxz;
		my_ALICE_STACK[top_ALICE_STACK].minz = oc->link0->minz;
		top_ALICE_STACK++;
	}
	if (oc->link1 != NULL) {
		my_ALICE_STACK[top_ALICE_STACK].link = (oc->link1);
		my_ALICE_STACK[top_ALICE_STACK].maxx = oc->link1->maxx;
		my_ALICE_STACK[top_ALICE_STACK].minx = oc->link1->minx;
		my_ALICE_STACK[top_ALICE_STACK].maxy = oc->link1->maxy;
		my_ALICE_STACK[top_ALICE_STACK].miny = oc->link1->miny;
		my_ALICE_STACK[top_ALICE_STACK].maxz = oc->link1->maxz;
		my_ALICE_STACK[top_ALICE_STACK].minz = oc->link1->minz;
		top_ALICE_STACK++;
	}
	if (oc->link2 != NULL) {
		my_ALICE_STACK[top_ALICE_STACK].link = (oc->link2);
		my_ALICE_STACK[top_ALICE_STACK].maxx = oc->link2->maxx;
		my_ALICE_STACK[top_ALICE_STACK].minx = oc->link2->minx;
		my_ALICE_STACK[top_ALICE_STACK].maxy = oc->link2->maxy;
		my_ALICE_STACK[top_ALICE_STACK].miny = oc->link2->miny;
		my_ALICE_STACK[top_ALICE_STACK].maxz = oc->link2->maxz;
		my_ALICE_STACK[top_ALICE_STACK].minz = oc->link2->minz;
		top_ALICE_STACK++;
	}
	if (oc->link3 != NULL) {
		my_ALICE_STACK[top_ALICE_STACK].link = (oc->link3);
		my_ALICE_STACK[top_ALICE_STACK].maxx = oc->link3->maxx;
		my_ALICE_STACK[top_ALICE_STACK].minx = oc->link3->minx;
		my_ALICE_STACK[top_ALICE_STACK].maxy = oc->link3->maxy;
		my_ALICE_STACK[top_ALICE_STACK].miny = oc->link3->miny;
		my_ALICE_STACK[top_ALICE_STACK].maxz = oc->link3->maxz;
		my_ALICE_STACK[top_ALICE_STACK].minz = oc->link3->minz;
		top_ALICE_STACK++;
	}
	if (oc->link4 != NULL) {
		my_ALICE_STACK[top_ALICE_STACK].link = (oc->link4);
		my_ALICE_STACK[top_ALICE_STACK].maxx = oc->link4->maxx;
		my_ALICE_STACK[top_ALICE_STACK].minx = oc->link4->minx;
		my_ALICE_STACK[top_ALICE_STACK].maxy = oc->link4->maxy;
		my_ALICE_STACK[top_ALICE_STACK].miny = oc->link4->miny;
		my_ALICE_STACK[top_ALICE_STACK].maxz = oc->link4->maxz;
		my_ALICE_STACK[top_ALICE_STACK].minz = oc->link4->minz;
		top_ALICE_STACK++;
	}
	if (oc->link5 != NULL) {
		my_ALICE_STACK[top_ALICE_STACK].link = (oc->link5);
		my_ALICE_STACK[top_ALICE_STACK].maxx = oc->link5->maxx;
		my_ALICE_STACK[top_ALICE_STACK].minx = oc->link5->minx;
		my_ALICE_STACK[top_ALICE_STACK].maxy = oc->link5->maxy;
		my_ALICE_STACK[top_ALICE_STACK].miny = oc->link5->miny;
		my_ALICE_STACK[top_ALICE_STACK].maxz = oc->link5->maxz;
		my_ALICE_STACK[top_ALICE_STACK].minz = oc->link5->minz;
		top_ALICE_STACK++;
	}
	if (oc->link6 != NULL) {
		my_ALICE_STACK[top_ALICE_STACK].link = (oc->link6);
		my_ALICE_STACK[top_ALICE_STACK].maxx = oc->link6->maxx;
		my_ALICE_STACK[top_ALICE_STACK].minx = oc->link6->minx;
		my_ALICE_STACK[top_ALICE_STACK].maxy = oc->link6->maxy;
		my_ALICE_STACK[top_ALICE_STACK].miny = oc->link6->miny;
		my_ALICE_STACK[top_ALICE_STACK].maxz = oc->link6->maxz;
		my_ALICE_STACK[top_ALICE_STACK].minz = oc->link6->minz;
		top_ALICE_STACK++;
	}
	if (oc->link7 != NULL) {
		my_ALICE_STACK[top_ALICE_STACK].link = (oc->link7);
		my_ALICE_STACK[top_ALICE_STACK].maxx = oc->link7->maxx;
		my_ALICE_STACK[top_ALICE_STACK].minx = oc->link7->minx;
		my_ALICE_STACK[top_ALICE_STACK].maxy = oc->link7->maxy;
		my_ALICE_STACK[top_ALICE_STACK].miny = oc->link7->miny;
		my_ALICE_STACK[top_ALICE_STACK].maxz = oc->link7->maxz;
		my_ALICE_STACK[top_ALICE_STACK].minz = oc->link7->minz;
		top_ALICE_STACK++;
	}
	while (top_ALICE_STACK > 0) {
		if (my_ALICE_STACK[top_ALICE_STACK - 1].link != NULL) {
			if (my_ALICE_STACK[top_ALICE_STACK - 1].link->dlist == true) {

				octTree* octree1 = my_ALICE_STACK[top_ALICE_STACK - 1].link;
				TOCHKA p;
				p.x = 0.125*(octree1->p0.x + octree1->p1.x + octree1->p2.x + octree1->p3.x + octree1->p4.x + octree1->p5.x + octree1->p6.x + octree1->p7.x);
				p.y = 0.125*(octree1->p0.y + octree1->p1.y + octree1->p2.y + octree1->p3.y + octree1->p4.y + octree1->p5.y + octree1->p6.y + octree1->p7.y);
				p.z = 0.125*(octree1->p0.z + octree1->p1.z + octree1->p2.z + octree1->p3.z + octree1->p4.z + octree1->p5.z + octree1->p6.z + octree1->p7.z);
				integer ib;
				bool inDomain = false;

				inDomain = in_model_flow(p, ib, b, lb); // HYDRODINAMIC  CFD!!!

				if ((p.x < b[0].g.xS) || (p.x > b[0].g.xE) || (p.y < b[0].g.yS) || (p.y > b[0].g.yE) || (p.z < b[0].g.zS) || (p.z > b[0].g.zE)) {
					// Лежит за пределами кабинета.
					inDomain = false;
				}
				if (inDomain) {
					// Узел принадлекжит расчётной области и его номер maxelm==octree1->inum_FD.				

					// сосед не существует
					// записываем идентификатор грани и
					// увеличиваем счётчик граней.
					if (!octree1->b4E) {
						if ((octree1->linkE == NULL)||((octree1->linkE != NULL)&&(octree1->linkE->inum_FD == 0))) {
							G = ESIDE;
							// граничная грань:
							INIT_SOSEDI(sosedi, maxelm, ESIDE, GG, false, maxelm_memo + maxbound, NON_EXISTENT_NODE, NON_EXISTENT_NODE, NON_EXISTENT_NODE, false, NON_EXISTENT_NODE, NON_EXISTENT_NODE, NON_EXISTENT_NODE, NON_EXISTENT_NODE); // заполнение sosedi[G][maxelm].
							DEFINE_BOUNDARY_PROPERTIES_FLOW(prop, prop_b, maxelm, maxbound); // В граничный узел сносятся свойства прилегающего КО.
													
							doublereal dSloc = CALC_SQUARE(G, maxelm, nvtx, pa); // Вычисление площади грани ячейки.
							TOCHKA p_centerG = CALC_CENTERG(G, maxelm, nvtx, pa); // Вычисление геометрического центра грани.
							
							// G - грань по которой ставится граничное условие.
							// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
							// bound_id - номер граничного узла, начиная с нуля.
							// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
							if (octree1->b4W) {
								obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
									octree1->linkW0->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, false, false, p_centerG);
							}
							else {
								if (octree1->linkW == NULL) {
									obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
										NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, false, false, p_centerG);
								}
								else {
									obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
										octree1->linkW->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, false, false, p_centerG);
								}
							}
							maxbound++;
						}
						else {
							if (octree1->linkE != NULL) {
								// сосед существует.
								if (bvisit[octree1->linkE->inum_FD - 1]) {

									integer return_current_node_number = sosedi[WSIDE][octree1->linkE->inum_FD - 1].iNODE1;
									if (octree1->linkE->b4W) {
										integer current_node_number = octree1->inum_FD - 1;
										bool bcontinue = true;
										if (bcontinue) {
											if (octree1->linkE->linkW0!=NULL) {
												if (octree1->linkE->linkW0->inum_FD - 1 == current_node_number) {
													bcontinue = false; // iNODE1
													if (sosedi[WSIDE][octree1->linkE->inum_FD - 1].iNODE1 == current_node_number) {
														return_current_node_number = sosedi[WSIDE][octree1->linkE->inum_FD - 1].iNODE1;
													}
													else {
														printf("this can not be E W0\n");
														//getchar();
														system("PAUSE");
													}
												}
											}
										}
										if (bcontinue) {
											if (octree1->linkE->linkW3 != NULL) {
												if (octree1->linkE->linkW3->inum_FD - 1 == current_node_number) {
													bcontinue = false; // iNODE2
													if (sosedi[WSIDE][octree1->linkE->inum_FD - 1].iNODE2 == current_node_number) {
														return_current_node_number = sosedi[WSIDE][octree1->linkE->inum_FD - 1].iNODE2;
													}
													else {
														printf("this can not be E W3\n");
														//getchar();
														system("PAUSE");
													}
												}
											}
										}
										if (bcontinue) {
											if (octree1->linkE->linkW4 != NULL) {
												if (octree1->linkE->linkW4->inum_FD - 1 == current_node_number) {
													bcontinue = false; // iNODE3
													if (sosedi[WSIDE][octree1->linkE->inum_FD - 1].iNODE3 == current_node_number) {
														return_current_node_number = sosedi[WSIDE][octree1->linkE->inum_FD - 1].iNODE3;
													}
													else {
														printf("this can not be E W4\n");
														//getchar();
														system("PAUSE");
													}
												}
											}
										}
										if (bcontinue) {
											if (octree1->linkE->linkW7 != NULL) {
												if (octree1->linkE->linkW7->inum_FD - 1 == current_node_number) {
													bcontinue = false; // iNODE4
													if (sosedi[WSIDE][octree1->linkE->inum_FD - 1].iNODE4 == current_node_number) {
														return_current_node_number = sosedi[WSIDE][octree1->linkE->inum_FD - 1].iNODE4;
													}
													else {
														printf("this can not be E W7\n");
														//getchar();
														system("PAUSE");
													}
												}
											}
										}

									}


									// узел уже был посещён
									// Здесь неким образом модифицировалось gran_t
									if (return_current_node_number >= maxelm_memo) {

										printf("bvisit error  E ...\n");
										printf("model is incorrect. FATALL error!!!\n");
										printf("internal flat heat source can not be located\n");
										printf("inside the liquid.\n");
										printf("Please, press any key to exit...\n");
										//getchar();
										system("PAUSE");
										//exit(1);

										/*
										// Это граничная грань внутреннего источника тепла.
										G = ESIDE;
										// граничная грань:
										INIT_SOSEDI(sosedi, maxelm, ESIDE, GG, false, sosedi[WSIDE][octree1->linkE->inum_FD - 1].iNODE1, NON_EXISTENT_NODE, NON_EXISTENT_NODE, NON_EXISTENT_NODE, false, NON_EXISTENT_NODE, NON_EXISTENT_NODE, NON_EXISTENT_NODE, NON_EXISTENT_NODE); // заполнение sosedi[G][maxelm].
										
										integer ib1;

										prop_b[RHO][sosedi[G][maxelm].iNODE1 - maxelm_memo] = prop[RHO][maxelm];
										prop_b[MU][sosedi[G][maxelm].iNODE1 - maxelm_memo] = prop[MU][maxelm];
										prop_b[BETA_T][sosedi[G][maxelm].iNODE1 - maxelm_memo] = prop[BETA_T][maxelm];

										// TODO 22.09.2016
										// Проблема. Надо подумать над ситуацией.
										printf("nado podumat nad situaciei.\n");
										getchar();

										
											// FLUID
											prop_b[RHO][sosedi[G][maxelm].iNODE1 - maxelm_memo] = prop[RHO][octree1->linkE->inum_FD - 1];
											prop_b[CP][sosedi[G][maxelm].iNODE1 - maxelm_memo] = prop[CP][octree1->linkE->inum_FD - 1];
											prop_b[LAM][sosedi[G][maxelm].iNODE1 - maxelm_memo] = prop[LAM][octree1->linkE->inum_FD - 1];
											prop_b[MULT_LAM_X][sosedi[G][maxelm].iNODE1 - maxelm_memo] = prop[MULT_LAM_X][octree1->linkE->inum_FD - 1];
											prop_b[MULT_LAM_Y][sosedi[G][maxelm].iNODE1 - maxelm_memo] = prop[MULT_LAM_Y][octree1->linkE->inum_FD - 1];
											prop_b[MULT_LAM_Z][sosedi[G][maxelm].iNODE1 - maxelm_memo] = prop[MULT_LAM_Z][octree1->linkE->inum_FD - 1];
										


										doublereal dSloc = CALC_SQUARE(G, maxelm, nvtx, pa); // Вычисление площади грани ячейки.
										
										// G - грань по которой ставится граничное условие.
										// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
										// bound_id - номер граничного узла, начиная с нуля.
										// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
										if (octree1->b4W) {
											obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, sosedi[WSIDE][octree1->linkE->inum_FD - 1].iNODE1 - maxelm_memo,
												octree1->linkW0->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, true, true, p_centerG);
										}
										else {
											if (octree1->linkW == NULL) {
												obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, sosedi[WSIDE][octree1->linkE->inum_FD - 1].iNODE1 - maxelm_memo,
													NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, true, true, p_centerG);
											}
											else {
												obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, sosedi[WSIDE][octree1->linkE->inum_FD - 1].iNODE1 - maxelm_memo,
													octree1->linkW->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, true, true, p_centerG);
											}
										}

										*/

									}
									else {
										// Это строго внутренний узел.
										// Внутренний узел.
										G = ESIDE;
										// внутренняя грань (нумерация начинается с нуля):
										sosedi[G][maxelm].iNODE1 = octree1->linkE->inum_FD - 1; // граничные КО нумеруются в последнюю очередь,
										sosedi[G][maxelm].bdroblenie4 = false;
										sosedi[G][maxelm].iNODE2 = NON_EXISTENT_NODE;
										sosedi[G][maxelm].iNODE3 = NON_EXISTENT_NODE;
										sosedi[G][maxelm].iNODE4 = NON_EXISTENT_NODE;
										// после maxelm_memo внутренних КО.
										// Вычисление дальнего соседа
										CALC_GG(G, GG); // Вычисляет грань GG по грани G.
										// Определение дальнего соседа :
										// TODO 13 сентября 2016.
									}

								}
								else {
									bool binc;
									integer lsid;
									// узнать координаты центра грани и ориентацию в пространстве
									x_c = octree1->p1.x;
									y_c = 0.5*(octree1->p1.y + octree1->p2.y);
									z_c = 0.5*(octree1->p1.z + octree1->p5.z);
									iplane = YZ;
									patch_maxbound(iplane, s, ls, x_c, y_c, z_c, maxbound, binc, lsid);
									if (binc) {

										printf("model is incorrect. FATALL error!!!\n");
										printf("internal flat heat source can not be located\n");
										printf("inside the liquid.\n");
										printf("Please, press any key to exit...\n");
										//getchar();
										system("PAUSE");
										exit(1);

										/*
										// Внутренний источник тепла.
										G = ESIDE;
										// граничная грань:
										INIT_SOSEDI(sosedi, maxelm, ESIDE, GG, false, maxelm_memo + maxbound-1, NON_EXISTENT_NODE, NON_EXISTENT_NODE, NON_EXISTENT_NODE, false, NON_EXISTENT_NODE, NON_EXISTENT_NODE, NON_EXISTENT_NODE, NON_EXISTENT_NODE); // заполнение sosedi[G][maxelm].

										// Сосед E есть внутренний узел и он еще не был посещен.
										// На границе этих двух внутренних узлов расположен источник тепла.

										
											// FLUID
											prop_b[RHO][maxbound - 1] = prop[RHO][octree1->linkE->inum_FD - 1];
											prop_b[CP][maxbound - 1] = prop[CP][octree1->linkE->inum_FD - 1];
											prop_b[LAM][maxbound - 1] = prop[LAM][octree1->linkE->inum_FD - 1];
											prop_b[MULT_LAM_X][maxbound - 1] = prop[MULT_LAM_X][octree1->linkE->inum_FD - 1];
											prop_b[MULT_LAM_Y][maxbound - 1] = prop[MULT_LAM_Y][octree1->linkE->inum_FD - 1];
											prop_b[MULT_LAM_Z][maxbound - 1] = prop[MULT_LAM_Z][octree1->linkE->inum_FD - 1];
										


										doublereal dSloc = CALC_SQUARE(G, maxelm, nvtx, pa); // Вычисление площади грани ячейки.
										
										// G - грань по которой ставится граничное условие.
										// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
										// bound_id - номер граничного узла, начиная с нуля.
										// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
										if (octree1->b4W) {
											obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound - 1,
												octree1->linkW0->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, true, false, p_centerG);
										}
										else {
											if (octree1->linkW == NULL) {
												obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound - 1,
													NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, true, false, p_centerG);
											}
											else {
												obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound - 1,
													octree1->linkW->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, true, false, p_centerG);
											}
										}

										*/

									}
									else {
										// Внутренний узел.
										G = ESIDE;
										// внутренняя грань (нумерация начинается с нуля):
										sosedi[G][maxelm].iNODE1 = octree1->linkE->inum_FD - 1; // граничные КО нумеруются в последнюю очередь,
										sosedi[G][maxelm].bdroblenie4 = false;
										sosedi[G][maxelm].iNODE2 = NON_EXISTENT_NODE;
										sosedi[G][maxelm].iNODE3 = NON_EXISTENT_NODE;
										sosedi[G][maxelm].iNODE4 = NON_EXISTENT_NODE;
										// после maxelm_memo внутренних КО.
										// Вычисление дальнего соседа
										CALC_GG(G, GG); // Вычисляет грань GG по грани G.
										// Определение дальнего соседа :
										// TODO 13 сентября 2016.

									}

								}
							}
						}
					}
					else {
					    // Присутствуют все 4 соседа: E: 1,2,5,6.
					
						if (((octree1->linkE1 == NULL) && (octree1->linkE2 == NULL) && (octree1->linkE5 == NULL) && (octree1->linkE6 == NULL))||
							(((octree1->linkE1 != NULL) && (octree1->linkE2 != NULL) && (octree1->linkE5 != NULL) && (octree1->linkE6 != NULL) &&
							(octree1->linkE1->inum_FD == 0) && (octree1->linkE2->inum_FD == 0) && (octree1->linkE5->inum_FD == 0) && (octree1->linkE6->inum_FD == 0)))||
							(((octree1->linkE1 != NULL) && (octree1->linkE2 == NULL) && (octree1->linkE5 != NULL) && (octree1->linkE6 == NULL) &&
							(octree1->linkE1->inum_FD == 0) && (octree1->linkE5->inum_FD == 0))) ||
							(((octree1->linkE1 != NULL) && (octree1->linkE2 != NULL) && (octree1->linkE5 == NULL) && (octree1->linkE6 == NULL) &&
							(octree1->linkE1->inum_FD == 0) && (octree1->linkE2->inum_FD == 0)))||
							(((octree1->linkE1 != NULL) && (octree1->linkE2 == NULL) && (octree1->linkE5 == NULL) && (octree1->linkE6 == NULL) &&
							(octree1->linkE1->inum_FD == 0))))						
						{
							// 1 - Все четыре соседа на границе расчётной области.
							// 2 - Все четыре соседа лежат на границе HOLLOW блока.
							// 3 - E2 && E6 граница РО, E1 && E5 граница HOLLOW блока.
							// 4 - E1 && E2 граница HOLLOW блока, E5 && E6 граница РО.
							// 5 - E1 граница HOLLOW блока, E2 && E5 && E6 граница РО. 

							G = ESIDE;
							// граничная грань:
							INIT_SOSEDI(sosedi, maxelm, ESIDE, GG, false, maxelm_memo + maxbound, NON_EXISTENT_NODE, NON_EXISTENT_NODE, NON_EXISTENT_NODE, false, NON_EXISTENT_NODE, NON_EXISTENT_NODE, NON_EXISTENT_NODE, NON_EXISTENT_NODE); // заполнение sosedi[G][maxelm].
							DEFINE_BOUNDARY_PROPERTIES_FLOW(prop, prop_b, maxelm, maxbound); // В граничный узел сносятся свойства прилегающего КО.

							doublereal dSloc = CALC_SQUARE(G, maxelm, nvtx, pa); // Вычисление площади грани ячейки.
							TOCHKA p_centerG = CALC_CENTERG(G, maxelm, nvtx, pa); // Вычисление геометрического центра грани.

							// G - грань по которой ставится граничное условие.
							// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
							// bound_id - номер граничного узла, начиная с нуля.
							// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
							if (octree1->b4W) {
								obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
									octree1->linkW0->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, false, false, p_centerG);
							}
							else {
								if (octree1->linkW == NULL) {
									obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
										NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, false, false, p_centerG);
								}
								else {

									obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
										octree1->linkW->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, false, false, p_centerG);
								}
							}
							// Зачем дробить если граница заканчивается NULL.
							maxbound++; // Вместо четырёх один.
						}						
						else {

							bool binc1, binc2, binc3, binc4;
							integer lsid1, lsid2, lsid3, lsid4;

							if (octree1->linkE1 == NULL) {
#if doubleintprecision == 1
								printf("error NULL sosedi E1 elm=%lld\n", maxelm);
#else
								printf("error NULL sosedi E1 elm=%d\n", maxelm);
#endif
								
								//getchar();
								system("PAUSE");
								/*
								G = ESIDE;
								// граничная грань:
								sosedi[G][maxelm].iNODE1 = maxelm_memo + maxbound; // граничные КО нумеруются в последнюю очередь,
								sosedi[G][maxelm].bdroblenie4 = true;

								// после maxelm_memo внутренних КО.
								// Вычисление дальнего соседа
								CALC_GG(G,GG); // Вычисляет грань GG по грани G.
								sosedi[GG][maxelm].iNODE1 = NON_EXISTENT_NODE; // соседа нет.
								sosedi[GG][maxelm].bdroblenie4 = true;

								#if doubleintprecision == 1
										printf("error NULL sosedi E1 elm=%lld\n", maxelm);
								#else
										printf("error NULL sosedi E1 elm=%d\n", maxelm);
								#endif
								
								getchar();

								// Эта ситуация возможна и здесь нужно предусмотреть обработку.
								 DEFINE_BOUNDARY_PROPERTIES_FLOW(prop, prop_b, maxelm, maxbound); // В граничный узел сносятся свойства прилегающего КО.
								
								doublereal dS_loc = 0.0;

								// G - грань по которой ставится граничное условие.
								// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
								// bound_id - номер граничного узла, начиная с нуля.
								// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
								if (octree1->b4W) {
								obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
								octree1->linkW0->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
								}
								else {
								if (octree1->linkW == NULL) {
								obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
								NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
								}
								else {
								obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
								octree1->linkW->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
								}
								}

								maxbound++;
								*/
							}
							else {
								if (octree1->linkE1->inum_FD == 0) {
									G = ESIDE;
									// граничная грань:
									sosedi[G][maxelm].iNODE1 = maxelm_memo + maxbound; // граничные КО нумеруются в последнюю очередь,
									sosedi[G][maxelm].bdroblenie4 = true;

									// после maxelm_memo внутренних КО.
									// Вычисление дальнего соседа
									CALC_GG(G, GG); // Вычисляет грань GG по грани G.
									sosedi[GG][maxelm].iNODE1 = NON_EXISTENT_NODE; // соседа нет.
									sosedi[GG][maxelm].bdroblenie4 = true;

#if doubleintprecision == 1
									//printf("error inum_FD == 0 sosedi E1 elm=%lld\n", maxelm);
#else
									//printf("error inum_FD == 0 sosedi E1 elm=%d\n", maxelm);
#endif
									
									//getchar();

									// Эта ситуация возможна и здесь нужно предусмотреть обработку.
									DEFINE_BOUNDARY_PROPERTIES_FLOW(prop, prop_b, maxelm, maxbound); // В граничный узел сносятся свойства прилегающего КО.
									
									// Sy*Sz
									doublereal dS_loc = fabs(octree1->linkE1->p4.z - octree1->linkE1->p0.z)*fabs(octree1->linkE1->p3.y - octree1->linkE1->p0.y);
									TOCHKA p_centerG;
									p_centerG.z = 0.5*(octree1->linkE1->p4.z + octree1->linkE1->p0.z);
									p_centerG.y = 0.5*(octree1->linkE1->p3.y + octree1->linkE1->p0.y);
									p_centerG.x = octree1->linkE1->p0.x;

									// G - грань по которой ставится граничное условие.
									// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
									// bound_id - номер граничного узла, начиная с нуля.
									// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
									if (octree1->b4W) {
										obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
											octree1->linkW0->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
									}
									else {
										if (octree1->linkW == NULL) {
											obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
												NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
										}
										else {
											obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
												octree1->linkW->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
										}
									}

									maxbound++;
								}
								else {
									// сосед существует.
									if (bvisit[octree1->linkE1->inum_FD - 1]) {
										// узел уже был посещён
										// Здесь неким образом модифицировалось gran_t
										// Здесь не может быть внутреннего источника тепла иначе 
										// выдавалось бы предупреждение 
#if doubleintprecision == 1
										//printf("error internal source E1 elm=%lld\n", maxelm);
#else
										//printf("error internal source E1 elm=%d\n", maxelm);
#endif
										
										//getchar();
										// В узле octree1->linkE1->inum_FD - 1 Который был посещен ранее.
										// Здесь возможен только внутренний узел.
										// Внутренний узел.
										G = ESIDE;
										// Внутренняя грань:
										sosedi[G][maxelm].iNODE1 = octree1->linkE1->inum_FD - 1; // граничные КО нумеруются в последнюю очередь,
										sosedi[G][maxelm].bdroblenie4 = true;

										// после maxelm_memo внутренних КО.
										// Вычисление дальнего соседа
										CALC_GG(G, GG); // Вычисляет грань GG по грани G.

										// Определение дальнего соседа :
										// TODO 13 сентября 2016.

									}
									else {
										// узнать координаты центра грани и ориентацию в пространстве
										x_c = octree1->p1.x;
										y_c = 0.5*(octree1->linkE1->p0.y + octree1->linkE1->p3.y);
										z_c = 0.5*(octree1->linkE1->p0.z + octree1->linkE1->p4.z);
										iplane = YZ;
										patch_maxbound(iplane, s, ls, x_c, y_c, z_c, maxbound, binc1, lsid1);
										if (binc1) {
#if doubleintprecision == 1
											printf("error internal source E1 elm=%lld\n", maxelm);
#else
											printf("error internal source E1 elm=%d\n", maxelm);
#endif
											
											//getchar();
											system("PAUSE");
										}
										else {
											// Внутренний узел.
											G = ESIDE;
											// Внутренняя грань:
											sosedi[G][maxelm].iNODE1 = octree1->linkE1->inum_FD - 1; // граничные КО нумеруются в последнюю очередь,
											sosedi[G][maxelm].bdroblenie4 = true;

											// после maxelm_memo внутренних КО.
											// Вычисление дальнего соседа
											CALC_GG(G, GG); // Вычисляет грань GG по грани G.

											// Определение дальнего соседа :
											// TODO 13 сентября 2016.

										}
									}
								}
							}

							if (octree1->linkE2 == NULL) {
								// 21 сентября 2016 (вырождение ячейки).

								// Заглушка.
								G = ESIDE;
								// граничная грань:
								sosedi[G][maxelm].iNODE2 = NON_EXISTENT_NODE; // граничные КО нумеруются в последнюю очередь,
								sosedi[G][maxelm].bdroblenie4 = true;

								// после maxelm_memo внутренних КО.
								// Вычисление дальнего соседа
								CALC_GG(G, GG); // Вычисляет грань GG по грани G.
								sosedi[GG][maxelm].iNODE2 = NON_EXISTENT_NODE; // соседа нет.
								sosedi[GG][maxelm].bdroblenie4 = true;


								/*
								G = ESIDE;
								// граничная грань:
								sosedi[G][maxelm].iNODE2 = maxelm_memo + maxbound; // граничные КО нумеруются в последнюю очередь,
								sosedi[G][maxelm].bdroblenie4 = true;

								// после maxelm_memo внутренних КО.
								// Вычисление дальнего соседа
								CALC_GG(G,GG); // Вычисляет грань GG по грани G.
								sosedi[GG][maxelm].iNODE2 = NON_EXISTENT_NODE; // соседа нет.
								sosedi[GG][maxelm].bdroblenie4 = true;

								#if doubleintprecision == 1
									printf("error NULL sosedi E2 elm=%lld\n", maxelm);
								#else
									printf("error NULL sosedi E2 elm=%d\n", maxelm);
								#endif
								
								getchar();

								// Эта ситуация возможна и здесь нужно предусмотреть обработку.
								 DEFINE_BOUNDARY_PROPERTIES_FLOW(prop, prop_b, maxelm, maxbound); // В граничный узел сносятся свойства прилегающего КО.
								
								doublereal dS_loc = 0.0;

								// G - грань по которой ставится граничное условие.
								// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
								// bound_id - номер граничного узла, начиная с нуля.
								// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
								if (octree1->b4W) {
								obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
								octree1->linkW0->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
								}
								else {
								if (octree1->linkW == NULL) {
								obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
								NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
								}
								else {
								obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
								octree1->linkW->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
								}
								}



								maxbound++;
								*/
							}
							else {
								if (octree1->linkE2->inum_FD == 0) {
									G = ESIDE;
									// граничная грань:
									sosedi[G][maxelm].iNODE2 = maxelm_memo + maxbound; // граничные КО нумеруются в последнюю очередь,
									sosedi[G][maxelm].bdroblenie4 = true;

									// после maxelm_memo внутренних КО.
									// Вычисление дальнего соседа
									CALC_GG(G, GG); // Вычисляет грань GG по грани G.
									sosedi[GG][maxelm].iNODE2 = NON_EXISTENT_NODE; // соседа нет.
									sosedi[GG][maxelm].bdroblenie4 = true;

#if doubleintprecision == 1
									//printf("error inum_FD == 0 sosedi E2 elm=%lld\n", maxelm);
#else
									//printf("error inum_FD == 0 sosedi E2 elm=%d\n", maxelm);
#endif
									
									//getchar();

									// Эта ситуация возможна и здесь нужно предусмотреть обработку.
									DEFINE_BOUNDARY_PROPERTIES_FLOW(prop, prop_b, maxelm, maxbound); // В граничный узел сносятся свойства прилегающего КО.
									
									// Sy*Sz
									doublereal dS_loc = fabs(octree1->linkE2->p4.z - octree1->linkE2->p0.z)*fabs(octree1->linkE2->p3.y - octree1->linkE2->p0.y);
									TOCHKA p_centerG;
									p_centerG.z = 0.5*(octree1->linkE2->p4.z + octree1->linkE2->p0.z);
									p_centerG.y = 0.5*(octree1->linkE2->p3.y + octree1->linkE2->p0.y);
									p_centerG.x = octree1->linkE2->p0.x;


									// G - грань по которой ставится граничное условие.
									// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
									// bound_id - номер граничного узла, начиная с нуля.
									// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
									if (octree1->b4W) {
										obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
											octree1->linkW0->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
									}
									else {
										if (octree1->linkW == NULL) {
											obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
												NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
										}
										else {
											obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
												octree1->linkW->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
										}
									}

									maxbound++;
								}
								else {
									// сосед существует.
									if (bvisit[octree1->linkE2->inum_FD - 1]) {
										// узел уже был посещён
										// Здесь неким образом модифицировалось gran_t
										// Здесь тоже возможен лишь строго внутренний узел.
										// Внутренний узел.
										G = ESIDE;
										// Внутренняя грань:
										sosedi[G][maxelm].iNODE2 = octree1->linkE2->inum_FD - 1; // граничные КО нумеруются в последнюю очередь,
										sosedi[G][maxelm].bdroblenie4 = true;

										// после maxelm_memo внутренних КО.
										// Вычисление дальнего соседа
										CALC_GG(G, GG); // Вычисляет грань GG по грани G.

										// Определение дальнего соседа :
										// TODO 13 сентября 2016.

									}
									else {
										// узнать координаты центра грани и ориентацию в пространстве
										x_c = octree1->p1.x;
										y_c = 0.5*(octree1->linkE2->p0.y + octree1->linkE2->p3.y);
										z_c = 0.5*(octree1->linkE2->p0.z + octree1->linkE2->p4.z);
										iplane = YZ;
										patch_maxbound(iplane, s, ls, x_c, y_c, z_c, maxbound, binc2, lsid2);
										if (binc2) {
#if doubleintprecision == 1
											printf("error internal source E2 elm=%lld\n", maxelm);
#else
											printf("error internal source E2 elm=%d\n", maxelm);
#endif
											
											//getchar();
											system("PAUSE");
										}
										else {
											// Внутренний узел.
											G = ESIDE;
											// Внутренняя грань:
											sosedi[G][maxelm].iNODE2 = octree1->linkE2->inum_FD - 1; // граничные КО нумеруются в последнюю очередь,
											sosedi[G][maxelm].bdroblenie4 = true;

											// после maxelm_memo внутренних КО.
											// Вычисление дальнего соседа
											CALC_GG(G, GG); // Вычисляет грань GG по грани G.

											// Определение дальнего соседа :
											// TODO 13 сентября 2016.

										}
									}
								}
							}

							if (octree1->linkE5 == NULL) {

								// 21 сентября 2016 (вырождение ячейки).

								// Заглушка.
								G = ESIDE;
								// граничная грань:
								sosedi[G][maxelm].iNODE3 = NON_EXISTENT_NODE; // граничные КО нумеруются в последнюю очередь,
								sosedi[G][maxelm].bdroblenie4 = true;

								// после maxelm_memo внутренних КО.
								// Вычисление дальнего соседа
								CALC_GG(G, GG); // Вычисляет грань GG по грани G.
								sosedi[GG][maxelm].iNODE3 = NON_EXISTENT_NODE; // соседа нет.
								sosedi[GG][maxelm].bdroblenie4 = true;

								/*
								G = ESIDE;
								// граничная грань:
								sosedi[G][maxelm].iNODE3 = maxelm_memo + maxbound; // граничные КО нумеруются в последнюю очередь,
								sosedi[G][maxelm].bdroblenie4 = true;

								// после maxelm_memo внутренних КО.
								// Вычисление дальнего соседа
								CALC_GG(G,GG); // Вычисляет грань GG по грани G.
								sosedi[GG][maxelm].iNODE3 = NON_EXISTENT_NODE; // соседа нет.
								sosedi[GG][maxelm].bdroblenie4 = true;

								#if doubleintprecision == 1
										printf("error NULL sosedi E5 elm=%lld\n", maxelm);
								#else
										printf("error NULL sosedi E5 elm=%d\n", maxelm);
								#endif
								
								getchar();

								// Эта ситуация возможна и здесь нужно предусмотреть обработку.
								 DEFINE_BOUNDARY_PROPERTIES_FLOW(prop, prop_b, maxelm, maxbound); // В граничный узел сносятся свойства прилегающего КО.
								
								doublereal dS_loc = 0.0;

								// G - грань по которой ставится граничное условие.
								// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
								// bound_id - номер граничного узла, начиная с нуля.
								// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
								if (octree1->b4W) {
								obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
								octree1->linkW0->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
								}
								else {
								if (octree1->linkW == NULL) {
								obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
								NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
								}
								else {
								obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
								octree1->linkW->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
								}
								}

								maxbound++;
								*/
							}
							else {
								if (octree1->linkE5->inum_FD == 0) {
									G = ESIDE;
									// граничная грань:
									sosedi[G][maxelm].iNODE3 = maxelm_memo + maxbound; // граничные КО нумеруются в последнюю очередь,
									sosedi[G][maxelm].bdroblenie4 = true;

									// после maxelm_memo внутренних КО.
									// Вычисление дальнего соседа
									CALC_GG(G, GG); // Вычисляет грань GG по грани G.
									sosedi[GG][maxelm].iNODE3 = NON_EXISTENT_NODE; // соседа нет.
									sosedi[GG][maxelm].bdroblenie4 = true;

#if doubleintprecision == 1
									//printf("error inum_FD == 0 sosedi E5 elm=%lld\n", maxelm);
#else
									//printf("error inum_FD == 0 sosedi E5 elm=%d\n", maxelm);
#endif
									
									//getchar();

									// Эта ситуация возможна и здесь нужно предусмотреть обработку.
									DEFINE_BOUNDARY_PROPERTIES_FLOW(prop, prop_b, maxelm, maxbound); // В граничный узел сносятся свойства прилегающего КО.
									
									// Sy*Sz
									doublereal dS_loc = fabs(octree1->linkE5->p4.z - octree1->linkE5->p0.z)*fabs(octree1->linkE5->p3.y - octree1->linkE5->p0.y);
									TOCHKA p_centerG;
									p_centerG.z = 0.5*(octree1->linkE5->p4.z + octree1->linkE5->p0.z);
									p_centerG.y = 0.5*(octree1->linkE5->p3.y + octree1->linkE5->p0.y);
									p_centerG.x = octree1->linkE5->p0.x;

									// G - грань по которой ставится граничное условие.
									// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
									// bound_id - номер граничного узла, начиная с нуля.
									// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
									if (octree1->b4W) {
										obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
											octree1->linkW0->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
									}
									else {
										if (octree1->linkW == NULL) {
											obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
												NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
										}
										else {
											obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
												octree1->linkW->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
										}
									}

									maxbound++;
								}
								else {
									// сосед существует.
									if (bvisit[octree1->linkE5->inum_FD - 1]) {
										// узел уже был посещён
										// Здесь неким образом модифицировалось gran_t
										// Здесь тоже возможен лишь строго внутренний узел.
										// Внутренний узел.
										G = ESIDE;
										// Внутренняя грань:
										sosedi[G][maxelm].iNODE3 = octree1->linkE5->inum_FD - 1; // граничные КО нумеруются в последнюю очередь,
										sosedi[G][maxelm].bdroblenie4 = true;

										// после maxelm_memo внутренних КО.
										// Вычисление дальнего соседа
										CALC_GG(G, GG); // Вычисляет грань GG по грани G.

										// Определение дальнего соседа :
										// TODO 13 сентября 2016.
									}
									else {
										// узнать координаты центра грани и ориентацию в пространстве
										x_c = octree1->p1.x;
										y_c = 0.5*(octree1->linkE5->p0.y + octree1->linkE5->p3.y);
										z_c = 0.5*(octree1->linkE5->p0.z + octree1->linkE5->p4.z);
										iplane = YZ;
										patch_maxbound(iplane, s, ls, x_c, y_c, z_c, maxbound, binc3, lsid3);
										if (binc3) {
#if doubleintprecision == 1
											printf("error internal source E5 elm=%lld\n", maxelm);
#else
											printf("error internal source E5 elm=%d\n", maxelm);
#endif
											
											//getchar();
											system("PAUSE");
										}
										else {
											// Внутренний узел.
											G = ESIDE;
											// Внутренняя грань:
											sosedi[G][maxelm].iNODE3 = octree1->linkE5->inum_FD - 1; // граничные КО нумеруются в последнюю очередь,
											sosedi[G][maxelm].bdroblenie4 = true;

											// после maxelm_memo внутренних КО.
											// Вычисление дальнего соседа
											CALC_GG(G, GG); // Вычисляет грань GG по грани G.

											// Определение дальнего соседа :
											// TODO 13 сентября 2016.
										}
									}
								}
							}

							if (octree1->linkE6 == NULL) {

								// 21 сентября 2016 (вырождение ячейки).

								// Заглушка.
								G = ESIDE;
								// граничная грань:
								sosedi[G][maxelm].iNODE4 = NON_EXISTENT_NODE; // граничные КО нумеруются в последнюю очередь,
								sosedi[G][maxelm].bdroblenie4 = true;

								// после maxelm_memo внутренних КО.
								// Вычисление дальнего соседа
								CALC_GG(G, GG); // Вычисляет грань GG по грани G.
								sosedi[GG][maxelm].iNODE4 = NON_EXISTENT_NODE; // соседа нет.
								sosedi[GG][maxelm].bdroblenie4 = true;


								/*
								G = ESIDE;
								// граничная грань:
								sosedi[G][maxelm].iNODE4 = maxelm_memo + maxbound; // граничные КО нумеруются в последнюю очередь,
								sosedi[G][maxelm].bdroblenie4 = true;

								// после maxelm_memo внутренних КО.
								// Вычисление дальнего соседа
								CALC_GG(G,GG); // Вычисляет грань GG по грани G.
								sosedi[GG][maxelm].iNODE4 = NON_EXISTENT_NODE; // соседа нет.
								sosedi[GG][maxelm].bdroblenie4 = true;

								#if doubleintprecision == 1
									printf("error NULL sosedi E6 elm=%lld\n", maxelm);
								#else
									printf("error NULL sosedi E6 elm=%d\n", maxelm);
								#endif
								
								getchar();

								// Эта ситуация возможна и здесь нужно предусмотреть обработку.
								 DEFINE_BOUNDARY_PROPERTIES_FLOW(prop, prop_b, maxelm, maxbound); // В граничный узел сносятся свойства прилегающего КО.
								
								doublereal dS_loc = 0.0;

								// G - грань по которой ставится граничное условие.
								// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
								// bound_id - номер граничного узла, начиная с нуля.
								// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
								if (octree1->b4W) {
								obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
								octree1->linkW0->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
								}
								else {
								if (octree1->linkW == NULL) {
								obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
								NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
								}
								else {
								obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
								octree1->linkW->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
								}
								}


								maxbound++;
								*/
							}
							else {
								if (octree1->linkE6->inum_FD == 0) {
									G = ESIDE;
									// граничная грань:
									sosedi[G][maxelm].iNODE4 = maxelm_memo + maxbound; // граничные КО нумеруются в последнюю очередь,
									sosedi[G][maxelm].bdroblenie4 = true;

									// после maxelm_memo внутренних КО.
									// Вычисление дальнего соседа
									CALC_GG(G, GG); // Вычисляет грань GG по грани G.
									sosedi[GG][maxelm].iNODE4 = NON_EXISTENT_NODE; // соседа нет.
									sosedi[GG][maxelm].bdroblenie4 = true;

#if doubleintprecision == 1
									//printf("error inum_FD == 0 sosedi E6 elm=%lld\n", maxelm);
#else
									//printf("error inum_FD == 0 sosedi E6 elm=%d\n", maxelm);
#endif
									
									//getchar();

									// Эта ситуация возможна и здесь нужно предусмотреть обработку.
									DEFINE_BOUNDARY_PROPERTIES_FLOW(prop, prop_b, maxelm, maxbound); // В граничный узел сносятся свойства прилегающего КО.
									
									// Sy*Sz
									doublereal dS_loc = fabs(octree1->linkE6->p4.z - octree1->linkE6->p0.z)*fabs(octree1->linkE6->p3.y - octree1->linkE6->p0.y);
									TOCHKA p_centerG;
									p_centerG.z = 0.5*(octree1->linkE6->p4.z + octree1->linkE6->p0.z);
									p_centerG.y = 0.5*(octree1->linkE6->p3.y + octree1->linkE6->p0.y);
									p_centerG.x = octree1->linkE6->p0.x;

									// G - грань по которой ставится граничное условие.
									// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
									// bound_id - номер граничного узла, начиная с нуля.
									// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
									if (octree1->b4W) {
										obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
											octree1->linkW0->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
									}
									else {
										if (octree1->linkW == NULL) {
											obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
												NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
										}
										else {
											obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
												octree1->linkW->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
										}
									}

									maxbound++;
								}
								else {
									// сосед существует.
									if (bvisit[octree1->linkE6->inum_FD - 1]) {
										// узел уже был посещён
										// Здесь неким образом модифицировалось gran_t
										// Здесь тоже возможен лишь строго внутренний узел.
										// Внутренний узел.
										G = ESIDE;
										// Внутренняя грань:
										sosedi[G][maxelm].iNODE4 = octree1->linkE6->inum_FD - 1; // граничные КО нумеруются в последнюю очередь,
										sosedi[G][maxelm].bdroblenie4 = true;

										// после maxelm_memo внутренних КО.
										// Вычисление дальнего соседа
										CALC_GG(G, GG); // Вычисляет грань GG по грани G.

										// Определение дальнего соседа :
										// TODO 13 сентября 2016.
									}
									else {
										// узнать координаты центра грани и ориентацию в пространстве
										x_c = octree1->p1.x;
										y_c = 0.5*(octree1->linkE6->p0.y + octree1->linkE6->p3.y);
										z_c = 0.5*(octree1->linkE6->p0.z + octree1->linkE6->p4.z);
										iplane = YZ;
										patch_maxbound(iplane, s, ls, x_c, y_c, z_c, maxbound, binc4, lsid4);
										if (binc4) {
#if doubleintprecision == 1
											printf("error internal source E6 elm=%lld\n", maxelm);
#else
											printf("error internal source E6 elm=%d\n", maxelm);
#endif
											
											//getchar();
											system("PAUSE");
										}
										else {
											// Внутренний узел.
											G = ESIDE;
											// Внутренняя грань:
											sosedi[G][maxelm].iNODE4 = octree1->linkE6->inum_FD - 1; // граничные КО нумеруются в последнюю очередь,
											sosedi[G][maxelm].bdroblenie4 = true;

											// после maxelm_memo внутренних КО.
											// Вычисление дальнего соседа
											CALC_GG(G, GG); // Вычисляет грань GG по грани G.

											// Определение дальнего соседа :
											// TODO 13 сентября 2016.
										}
									}
								}
							}
						}
					}

					// сосед не существует
					// записываем идентификатор грани и
					// увеличиваем счётчик граней.
					if (!octree1->b4W) {
						if ((octree1->linkW == NULL)||((octree1->linkW != NULL)&&(octree1->linkW->inum_FD == 0))) {

							G = WSIDE;
							// граничная грань:
							INIT_SOSEDI(sosedi, maxelm, WSIDE, GG, false, maxelm_memo + maxbound, NON_EXISTENT_NODE, NON_EXISTENT_NODE, NON_EXISTENT_NODE, false, NON_EXISTENT_NODE, NON_EXISTENT_NODE, NON_EXISTENT_NODE, NON_EXISTENT_NODE); // заполнение sosedi[G][maxelm].							
							DEFINE_BOUNDARY_PROPERTIES_FLOW(prop, prop_b, maxelm, maxbound); // В граничный узел сносятся свойства прилегающего КО.
							
							doublereal dSloc = CALC_SQUARE(G, maxelm, nvtx, pa); // Вычисление площади грани ячейки.
							TOCHKA p_centerG = CALC_CENTERG(G, maxelm, nvtx, pa); // Вычисление геометрического центра грани.

							// G - грань по которой ставится граничное условие.
							// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
							// bound_id - номер граничного узла, начиная с нуля.
							// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
							if (octree1->b4E) {
								obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
									octree1->linkE1->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, false, false, p_centerG);
							}
							else {
								if (octree1->linkE == NULL) {
									obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
										NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, false, false, p_centerG);
								}
								else {
									obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
										octree1->linkE->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, false, false, p_centerG);
								}
							}

							maxbound++;
						}
						else {
							if (octree1->linkW != NULL) {

								// сосед существует.
								if (bvisit[octree1->linkW->inum_FD - 1]) {


									integer return_current_node_number = sosedi[ESIDE][octree1->linkW->inum_FD - 1].iNODE1;
									if (octree1->linkW->b4E) {
										integer current_node_number = octree1->inum_FD - 1;
										bool bcontinue = true;
										if (bcontinue) {
											if (octree1->linkW->linkE1 != NULL) {
												if (octree1->linkW->linkE1->inum_FD - 1 == current_node_number) {
													bcontinue = false; // iNODE1
													if (sosedi[ESIDE][octree1->linkW->inum_FD - 1].iNODE1 == current_node_number) {
														return_current_node_number = sosedi[ESIDE][octree1->linkW->inum_FD - 1].iNODE1;
													}
													else {
														printf("this can not be W E1\n");
														//getchar();
														system("PAUSE");
													}
												}
											}
										}
										if (bcontinue) {
											if (octree1->linkW->linkE2 != NULL) {
												if (octree1->linkW->linkE2->inum_FD - 1 == current_node_number) {
													bcontinue = false; // iNODE2
													if (sosedi[ESIDE][octree1->linkW->inum_FD - 1].iNODE2 == current_node_number) {
														return_current_node_number = sosedi[ESIDE][octree1->linkW->inum_FD - 1].iNODE2;
													}
													else {
														printf("this can not be W E2\n");
														//getchar();
														system("PAUSE");
													}
												}
											}
										}
										if (bcontinue) {
											if (octree1->linkW->linkE5 != NULL) {
												if (octree1->linkW->linkE5->inum_FD - 1 == current_node_number) {
													bcontinue = false; // iNODE3
													if (sosedi[ESIDE][octree1->linkW->inum_FD - 1].iNODE3 == current_node_number) {
														return_current_node_number = sosedi[ESIDE][octree1->linkW->inum_FD - 1].iNODE3;
													}
													else {
														printf("this can not be W E5\n");
														//getchar();
														system("PAUSE");
													}
												}
											}
										}
										if (bcontinue) {
											if (octree1->linkW->linkE6 != NULL) {
												if (octree1->linkW->linkE6->inum_FD - 1 == current_node_number) {
													bcontinue = false; // iNODE4
													if (sosedi[ESIDE][octree1->linkW->inum_FD - 1].iNODE4 == current_node_number) {
														return_current_node_number = sosedi[ESIDE][octree1->linkW->inum_FD - 1].iNODE4;
													}
													else {
														printf("this can not be W E6\n");
														//getchar();
														system("PAUSE");
													}
												}
											}
										}								

									}

									// узел уже был посещён
									// Здесь неким образом модифицировалось gran_t
									if (return_current_node_number >= maxelm_memo) {

										printf("bvisit error  W ...\n");
										printf("model is incorrect. FATALL error!!!\n");
										printf("internal flat heat source can not be located\n");
										printf("inside the liquid.\n");
										printf("Please, press any key to exit...\n");
										//getchar();
										system("PAUSE");
										//exit(1);

										/*

										// Это граничная грань внутреннего источника тепла.
										// Внутренний источник тепла.
										G = WSIDE;
										// граничная грань:
										INIT_SOSEDI(sosedi, maxelm, WSIDE, GG, false, sosedi[ESIDE][octree1->linkW->inum_FD - 1].iNODE1, NON_EXISTENT_NODE, NON_EXISTENT_NODE, NON_EXISTENT_NODE, false, NON_EXISTENT_NODE, NON_EXISTENT_NODE, NON_EXISTENT_NODE, NON_EXISTENT_NODE); // заполнение sosedi[G][maxelm].
																				
										
											// FLUID
											prop_b[RHO][sosedi[G][maxelm].iNODE1 - maxelm_memo] = prop[RHO][octree1->linkW->inum_FD - 1];
											prop_b[CP][sosedi[G][maxelm].iNODE1 - maxelm_memo] = prop[CP][octree1->linkW->inum_FD - 1];
											prop_b[LAM][sosedi[G][maxelm].iNODE1 - maxelm_memo] = prop[LAM][octree1->linkW->inum_FD - 1];
											prop_b[MULT_LAM_X][sosedi[G][maxelm].iNODE1 - maxelm_memo] = prop[MULT_LAM_X][octree1->linkW->inum_FD - 1];
											prop_b[MULT_LAM_Y][sosedi[G][maxelm].iNODE1 - maxelm_memo] = prop[MULT_LAM_Y][octree1->linkW->inum_FD - 1];
											prop_b[MULT_LAM_Z][sosedi[G][maxelm].iNODE1 - maxelm_memo] = prop[MULT_LAM_Z][octree1->linkW->inum_FD - 1];
										

										doublereal dSloc = CALC_SQUARE(G, maxelm, nvtx, pa); // Вычисление площади грани ячейки.
										TOCHKA p_centerG = CALC_CENTERG(G, maxelm, nvtx, pa); // Вычисление геометрического центра грани.

										// G - грань по которой ставится граничное условие.
										// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
										// bound_id - номер граничного узла, начиная с нуля.
										// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
										if (octree1->b4E) {
											obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, sosedi[ESIDE][octree1->linkW->inum_FD - 1].iNODE1 - maxelm_memo,
												octree1->linkE1->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, true, true, p_centerG);
										}
										else {
											if (octree1->linkE == NULL) {
												obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, sosedi[ESIDE][octree1->linkW->inum_FD - 1].iNODE1 - maxelm_memo,
													NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, true, true, p_centerG);
											}
											else {
												obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, sosedi[ESIDE][octree1->linkW->inum_FD - 1].iNODE1 - maxelm_memo,
													octree1->linkE->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, true, true, p_centerG);
											}
										}

										*/

									}
									else {
										// Это строго внутренний узел.
										G = WSIDE;
										// внутренняя грань:
										sosedi[G][maxelm].iNODE1 = octree1->linkW->inum_FD - 1; // граничные КО нумеруются в последнюю очередь,
										sosedi[G][maxelm].bdroblenie4 = false;
										sosedi[G][maxelm].iNODE2 = NON_EXISTENT_NODE;
										sosedi[G][maxelm].iNODE3 = NON_EXISTENT_NODE;
										sosedi[G][maxelm].iNODE4 = NON_EXISTENT_NODE;
										// после maxelm_memo внутренних КО.
										// Вычисление дальнего соседа
										CALC_GG(G, GG); // Вычисляет грань GG по грани G.
										// Определение дальнего соседа :
										// TODO 13 сентября 2016.
									}

								}
								else {
									bool binc;
									integer lsid;
									// узнать координаты центра грани и ориентацию в пространстве
									x_c = octree1->p0.x;
									y_c = 0.5*(octree1->p0.y + octree1->p3.y);
									z_c = 0.5*(octree1->p0.z + octree1->p4.z);
									iplane = YZ;
									patch_maxbound(iplane, s, ls, x_c, y_c, z_c, maxbound, binc, lsid);
									if (binc) {

										printf("model is incorrect. FATALL error!!!\n");
										printf("internal flat heat source can not be located\n");
										printf("inside the liquid.\n");
										printf("Please, press any key to exit...\n");
										//getchar();
										system("PAUSE");
										exit(1);

										/*
										// Внутренний источник тепла.
										G = WSIDE;
										// граничная грань:
										INIT_SOSEDI(sosedi, maxelm, WSIDE, GG, false, maxelm_memo + maxbound-1, NON_EXISTENT_NODE, NON_EXISTENT_NODE, NON_EXISTENT_NODE, false, NON_EXISTENT_NODE, NON_EXISTENT_NODE, NON_EXISTENT_NODE, NON_EXISTENT_NODE); // заполнение sosedi[G][maxelm].

										
											// FLUID
											prop_b[RHO][maxbound - 1] = prop[RHO][octree1->linkW->inum_FD - 1];
											prop_b[CP][maxbound - 1] = prop[CP][octree1->linkW->inum_FD - 1];
											prop_b[LAM][maxbound - 1] = prop[LAM][octree1->linkW->inum_FD - 1];
											prop_b[MULT_LAM_X][maxbound - 1] = prop[MULT_LAM_X][octree1->linkW->inum_FD - 1];
											prop_b[MULT_LAM_Y][maxbound - 1] = prop[MULT_LAM_Y][octree1->linkW->inum_FD - 1];
											prop_b[MULT_LAM_Z][maxbound - 1] = prop[MULT_LAM_Z][octree1->linkW->inum_FD - 1];
										

										doublereal dSloc = CALC_SQUARE(G, maxelm, nvtx, pa); // Вычисление площади грани ячейки.
										TOCHKA p_centerG = CALC_CENTERG(G, maxelm, nvtx, pa); // Вычисление геометрического центра грани.
										
										// G - грань по которой ставится граничное условие.
										// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
										// bound_id - номер граничного узла, начиная с нуля.
										// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
										if (octree1->b4E) {
											obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound - 1,
												octree1->linkE1->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, true, false, p_centerG);
										}
										else {
											if (octree1->linkE == NULL) {
												obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound - 1,
													NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, true, false, p_centerG);
											}
											else {
												obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound - 1,
													octree1->linkE->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, true, false, p_centerG);
											}
										}

										*/

									}
									else {
										// Внутренний узел.
										G = WSIDE;
										// внутренняя грань:
										sosedi[G][maxelm].iNODE1 = octree1->linkW->inum_FD - 1; // граничные КО нумеруются в последнюю очередь,
										sosedi[G][maxelm].bdroblenie4 = false;
										sosedi[G][maxelm].iNODE2 = NON_EXISTENT_NODE;
										sosedi[G][maxelm].iNODE3 = NON_EXISTENT_NODE;
										sosedi[G][maxelm].iNODE4 = NON_EXISTENT_NODE;
										// после maxelm_memo внутренних КО.
										// Вычисление дальнего соседа
										CALC_GG(G, GG); // Вычисляет грань GG по грани G.
										// Определение дальнего соседа :
										// TODO 13 сентября 2016.

									}
								}
							}
						}
					}
					else {
					// Присутствуют все 4 соседа: W: 0,3,4,7.

						if (((octree1->linkW0 == NULL) && (octree1->linkW3 == NULL) && (octree1->linkW4 == NULL) && (octree1->linkW7 == NULL))||
							(((octree1->linkW0 != NULL) && (octree1->linkW3 != NULL) && (octree1->linkW4 != NULL) && (octree1->linkW7 != NULL) &&
							(octree1->linkW0->inum_FD == 0) && (octree1->linkW3->inum_FD == 0) && (octree1->linkW4->inum_FD == 0) && (octree1->linkW7->inum_FD == 0)))||
							(((octree1->linkW0 != NULL) && (octree1->linkW3 == NULL) && (octree1->linkW4 != NULL) && (octree1->linkW7 == NULL) &&
							(octree1->linkW0->inum_FD == 0) && (octree1->linkW4->inum_FD == 0))) || 
							(((octree1->linkW0 != NULL) && (octree1->linkW3 != NULL) && (octree1->linkW4 == NULL) && (octree1->linkW7 == NULL) &&
							(octree1->linkW0->inum_FD == 0) && (octree1->linkW3->inum_FD == 0)))||
							(((octree1->linkW0 != NULL) && (octree1->linkW3 == NULL) && (octree1->linkW4 == NULL) && (octree1->linkW7 == NULL) &&
							(octree1->linkW0->inum_FD == 0))))
						{
							// 1 - W0&&W3&&W4&&W7 - на границе расчётной области (РО).
							// 2 -  W0&&W3&&W4&&W7 - на гранеице HOLLOW блока.
							// 3 - W0&&W4 на границе HOLLOW блока, W3&&W7 на границе РО.
							// 4 - W0&&W3 на границе HOLLOW блока, W4&&W7 на границе РО.
							// 5 - W0 на границе HOLLOW блока, W3&&W4&&W7  на границе РО.

							G = WSIDE;
							// граничная грань:
							INIT_SOSEDI(sosedi, maxelm, WSIDE, GG, false, maxelm_memo + maxbound, NON_EXISTENT_NODE, NON_EXISTENT_NODE, NON_EXISTENT_NODE, false, NON_EXISTENT_NODE, NON_EXISTENT_NODE, NON_EXISTENT_NODE, NON_EXISTENT_NODE); // заполнение sosedi[G][maxelm].
							DEFINE_BOUNDARY_PROPERTIES_FLOW(prop, prop_b, maxelm, maxbound); // В граничный узел сносятся свойства прилегающего КО.
													
							doublereal dSloc = CALC_SQUARE(G, maxelm, nvtx, pa); // Вычисление площади грани ячейки.
							TOCHKA p_centerG = CALC_CENTERG(G, maxelm, nvtx, pa); // Вычисление геометрического центра грани.
							
							// G - грань по которой ставится граничное условие.
							// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
							// bound_id - номер граничного узла, начиная с нуля.
							// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
							if (octree1->b4E) {
								obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
									octree1->linkE1->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, false, false, p_centerG);
							}
							else {
								if (octree1->linkE == NULL) {
									obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
										NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, false, false, p_centerG);
								}
								else {
									obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
										octree1->linkE->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, false, false, p_centerG);
								}
							}

							// Зачем дробить если граница заканчивается NULL.
							maxbound++; // Вместо четырёх один.
						}						
						else {

							bool binc1, binc2, binc3, binc4;
							integer lsid1, lsid2, lsid3, lsid4;

							if (octree1->linkW0 == NULL) {

#if doubleintprecision == 1
								printf("error NULL sosedi W0 elm=%lld\n", maxelm);
#else
								printf("error NULL sosedi W0 elm=%d\n", maxelm);
#endif
								
								//getchar();
								system("PAUSE");

								/*
								G = WSIDE;
								// граничная грань:
								sosedi[G][maxelm].iNODE1 = maxelm_memo + maxbound; // граничные КО нумеруются в последнюю очередь,
								sosedi[G][maxelm].bdroblenie4 = true;
								// после maxelm_memo внутренних КО.
								// Вычисление дальнего соседа
								CALC_GG(G,GG); // Вычисляет грань GG по грани G.
								sosedi[GG][maxelm].iNODE1 = NON_EXISTENT_NODE; // соседа нет.
								sosedi[GG][maxelm].bdroblenie4 = true;
								#if doubleintprecision == 1
										printf("error NULL sosedi W0 elm=%lld\n", maxelm);
								#else
										printf("error NULL sosedi W0 elm=%d\n", maxelm);
								#endif
								
								getchar();

								// Эта ситуация возможна и здесь надо предусмотреть её обработку.
								 DEFINE_BOUNDARY_PROPERTIES_FLOW(prop, prop_b, maxelm, maxbound); // В граничный узел сносятся свойства прилегающего КО.
								
								doublereal dS_loc = 0.0;

								// G - грань по которой ставится граничное условие.
								// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
								// bound_id - номер граничного узла, начиная с нуля.
								// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
								if (octree1->b4E) {
								obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
								octree1->linkE1->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
								}
								else {
								if (octree1->linkE == NULL) {
								obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
								NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
								}
								else {
								obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
								octree1->linkE->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
								}
								}

								maxbound++;
								*/
							}
							else {
								if (octree1->linkW0->inum_FD == 0) {
									G = WSIDE;
									// граничная грань:
									sosedi[G][maxelm].iNODE1 = maxelm_memo + maxbound; // граничные КО нумеруются в последнюю очередь,
									sosedi[G][maxelm].bdroblenie4 = true;
									// после maxelm_memo внутренних КО.
									// Вычисление дальнего соседа
									CALC_GG(G, GG); // Вычисляет грань GG по грани G.
									sosedi[GG][maxelm].iNODE1 = NON_EXISTENT_NODE; // соседа нет.
									sosedi[GG][maxelm].bdroblenie4 = true;
#if doubleintprecision == 1
									//printf("error inum_FD == 0 sosedi W0 elm=%lld\n", maxelm);
#else
									//printf("error inum_FD == 0 sosedi W0 elm=%d\n", maxelm);
#endif
									
									//getchar();

									// Эта ситуация возможна и здесь надо предусмотреть её обработку.
									DEFINE_BOUNDARY_PROPERTIES_FLOW(prop, prop_b, maxelm, maxbound); // В граничный узел сносятся свойства прилегающего КО.
									
									// Sy*Sz
									doublereal dS_loc = fabs(octree1->linkW0->p5.z - octree1->linkW0->p1.z)*fabs(octree1->linkW0->p2.y - octree1->linkW0->p1.y);
									TOCHKA p_centerG;
									p_centerG.z = 0.5*(octree1->linkW0->p5.z + octree1->linkW0->p1.z);
									p_centerG.y = 0.5*(octree1->linkW0->p2.y + octree1->linkW0->p1.y);
									p_centerG.x = octree1->linkW0->p1.x;

									// G - грань по которой ставится граничное условие.
									// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
									// bound_id - номер граничного узла, начиная с нуля.
									// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
									if (octree1->b4E) {
										obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
											octree1->linkE1->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
									}
									else {
										if (octree1->linkE == NULL) {
											obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
												NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
										}
										else {
											obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
												octree1->linkE->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
										}
									}

									maxbound++;
								}
								else {
									// сосед существует.
									if (bvisit[octree1->linkW0->inum_FD - 1]) {
										// узел уже был посещён
										// Здесь неким образом модифицировалось gran_t
										// Только внутренняя грань.
										G = WSIDE;
										// граничная грань:
										sosedi[G][maxelm].iNODE1 = octree1->linkW0->inum_FD - 1; // граничные КО нумеруются в последнюю очередь,
										sosedi[G][maxelm].bdroblenie4 = true;
										// после maxelm_memo внутренних КО.
										// Вычисление дальнего соседа
										CALC_GG(G, GG); // Вычисляет грань GG по грани G.

										// Определение дальнего соседа :
										// TODO 13 сентября 2016.

									}
									else {
										// узнать координаты центра грани и ориентацию в пространстве
										x_c = octree1->p0.x;
										y_c = 0.5*(octree1->linkW0->p1.y + octree1->linkW0->p2.y);
										z_c = 0.5*(octree1->linkW0->p1.z + octree1->linkW0->p5.z);
										iplane = YZ;
										patch_maxbound(iplane, s, ls, x_c, y_c, z_c, maxbound, binc1, lsid1);
										if (binc1) {
#if doubleintprecision == 1
											printf("error internal source W0 elm=%lld\n", maxelm);
#else
											printf("error internal source W0 elm=%d\n", maxelm);
#endif
											
											//getchar();
											system("PAUSE");
										}
										else {
											// Только внутренняя грань.
											G = WSIDE;
											// граничная грань:
											sosedi[G][maxelm].iNODE1 = octree1->linkW0->inum_FD - 1; // граничные КО нумеруются в последнюю очередь,
											sosedi[G][maxelm].bdroblenie4 = true;
											// после maxelm_memo внутренних КО.
											// Вычисление дальнего соседа
											CALC_GG(G, GG); // Вычисляет грань GG по грани G.

											// Определение дальнего соседа :
											// TODO 13 сентября 2016.
										}
									}
								}
							}

							if (octree1->linkW3 == NULL) {

								// 21 сентября 2016 (вырождение ячейки).

								// Заглушка.
								G = WSIDE;
								// граничная грань:
								sosedi[G][maxelm].iNODE2 = NON_EXISTENT_NODE; // граничные КО нумеруются в последнюю очередь,
								sosedi[G][maxelm].bdroblenie4 = true;
								// после maxelm_memo внутренних КО.
								// Вычисление дальнего соседа
								CALC_GG(G, GG); // Вычисляет грань GG по грани G.
								sosedi[GG][maxelm].iNODE2 = NON_EXISTENT_NODE; // соседа нет.
								sosedi[GG][maxelm].bdroblenie4 = true;

								/*
								G = WSIDE;
								// граничная грань:
								sosedi[G][maxelm].iNODE2 = maxelm_memo + maxbound; // граничные КО нумеруются в последнюю очередь,
								sosedi[G][maxelm].bdroblenie4 = true;
								// после maxelm_memo внутренних КО.
								// Вычисление дальнего соседа
								CALC_GG(G,GG); // Вычисляет грань GG по грани G.
								sosedi[GG][maxelm].iNODE2 = NON_EXISTENT_NODE; // соседа нет.
								sosedi[GG][maxelm].bdroblenie4 = true;

								#if doubleintprecision == 1
									printf("error NULL sosedi W3 elm=%lld\n", maxelm);
								#else
									printf("error NULL sosedi W3 elm=%d\n", maxelm);
								#endif
								
								getchar();

								// Эта ситуация возможна и здесь надо предусмотреть её обработку.
								 DEFINE_BOUNDARY_PROPERTIES_FLOW(prop, prop_b, maxelm, maxbound); // В граничный узел сносятся свойства прилегающего КО.
								
								doublereal dS_loc = 0.0;

								// G - грань по которой ставится граничное условие.
								// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
								// bound_id - номер граничного узла, начиная с нуля.
								// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
								if (octree1->b4E) {
								obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
								octree1->linkE1->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
								}
								else {
								if (octree1->linkE == NULL) {
								obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
								NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
								}
								else {
								obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
								octree1->linkE->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
								}
								}

								maxbound++;
								*/

							}
							else {
								if (octree1->linkW3->inum_FD == 0) {
									G = WSIDE;
									// граничная грань:
									sosedi[G][maxelm].iNODE2 = maxelm_memo + maxbound; // граничные КО нумеруются в последнюю очередь,
									sosedi[G][maxelm].bdroblenie4 = true;
									// после maxelm_memo внутренних КО.
									// Вычисление дальнего соседа
									CALC_GG(G, GG); // Вычисляет грань GG по грани G.
									sosedi[GG][maxelm].iNODE2 = NON_EXISTENT_NODE; // соседа нет.
									sosedi[GG][maxelm].bdroblenie4 = true;

#if doubleintprecision == 1
									//printf("error inum_FD == 0 sosedi W3 elm=%lld\n", maxelm);
#else
									//printf("error inum_FD == 0 sosedi W3 elm=%d\n", maxelm);
#endif
									
									//getchar();

									// Эта ситуация возможна и здесь надо предусмотреть её обработку.
									DEFINE_BOUNDARY_PROPERTIES_FLOW(prop, prop_b, maxelm, maxbound); // В граничный узел сносятся свойства прилегающего КО.
									
									// Sy*Sz
									doublereal dS_loc = fabs(octree1->linkW3->p5.z - octree1->linkW3->p1.z)*fabs(octree1->linkW3->p2.y - octree1->linkW3->p1.y);
									TOCHKA p_centerG;
									p_centerG.z = 0.5*(octree1->linkW3->p5.z + octree1->linkW3->p1.z);
									p_centerG.y = 0.5*(octree1->linkW3->p2.y + octree1->linkW3->p1.y);
									p_centerG.x = octree1->linkW3->p1.x;


									// G - грань по которой ставится граничное условие.
									// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
									// bound_id - номер граничного узла, начиная с нуля.
									// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
									if (octree1->b4E) {
										obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
											octree1->linkE1->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
									}
									else {
										if (octree1->linkE == NULL) {
											obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
												NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
										}
										else {
											obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
												octree1->linkE->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
										}
									}


									maxbound++;
								}
								else {
									// сосед существует.
									if (bvisit[octree1->linkW3->inum_FD - 1]) {
										// узел уже был посещён
										// Здесь неким образом модифицировалось gran_t
										// Только внутренняя грань.
										G = WSIDE;
										// граничная грань:
										sosedi[G][maxelm].iNODE2 = octree1->linkW3->inum_FD - 1; // граничные КО нумеруются в последнюю очередь,
										sosedi[G][maxelm].bdroblenie4 = true;
										// после maxelm_memo внутренних КО.
										// Вычисление дальнего соседа
										CALC_GG(G, GG); // Вычисляет грань GG по грани G.

										// Определение дальнего соседа :
										// TODO 13 сентября 2016.
									}
									else {
										// узнать координаты центра грани и ориентацию в пространстве
										x_c = octree1->p0.x;
										y_c = 0.5*(octree1->linkW3->p1.y + octree1->linkW3->p2.y);
										z_c = 0.5*(octree1->linkW3->p1.z + octree1->linkW3->p5.z);
										iplane = YZ;
										patch_maxbound(iplane, s, ls, x_c, y_c, z_c, maxbound, binc2, lsid2);
										if (binc2) {
#if doubleintprecision == 1
											printf("error internal source W3 elm=%lld\n", maxelm);
#else
											printf("error internal source W3 elm=%d\n", maxelm);
#endif
											
											//getchar();
											system("PAUSE");
										}
										else {
											// Только внутренняя грань.
											G = WSIDE;
											// граничная грань:
											sosedi[G][maxelm].iNODE2 = octree1->linkW3->inum_FD - 1; // граничные КО нумеруются в последнюю очередь,
											sosedi[G][maxelm].bdroblenie4 = true;
											// после maxelm_memo внутренних КО.
											// Вычисление дальнего соседа
											CALC_GG(G, GG); // Вычисляет грань GG по грани G.

											// Определение дальнего соседа :
											// TODO 13 сентября 2016.
										}
									}
								}
							}

							if (octree1->linkW4 == NULL) {

								// 21 сентября 2016 (вырождение ячейки).

								// Заглушка.
								G = WSIDE;
								// граничная грань:
								sosedi[G][maxelm].iNODE3 = NON_EXISTENT_NODE; // граничные КО нумеруются в последнюю очередь,
								sosedi[G][maxelm].bdroblenie4 = true;
								// после maxelm_memo внутренних КО.
								// Вычисление дальнего соседа
								CALC_GG(G, GG); // Вычисляет грань GG по грани G.
								sosedi[GG][maxelm].iNODE3 = NON_EXISTENT_NODE; // соседа нет.
								sosedi[GG][maxelm].bdroblenie4 = true;

								/*
								G = WSIDE;
								// граничная грань:
								sosedi[G][maxelm].iNODE3 = maxelm_memo + maxbound; // граничные КО нумеруются в последнюю очередь,
								sosedi[G][maxelm].bdroblenie4 = true;
								// после maxelm_memo внутренних КО.
								// Вычисление дальнего соседа
								CALC_GG(G,GG); // Вычисляет грань GG по грани G.
								sosedi[GG][maxelm].iNODE3 = NON_EXISTENT_NODE; // соседа нет.
								sosedi[GG][maxelm].bdroblenie4 = true;

								#if doubleintprecision == 1
										printf("error NULL sosedi W4 elm=%lld\n", maxelm);
								#else
										printf("error NULL sosedi W4 elm=%d\n", maxelm);
								#endif
								
								getchar();

								// Эта ситуация возможна и здесь надо предусмотреть её обработку.
								 DEFINE_BOUNDARY_PROPERTIES_FLOW(prop, prop_b, maxelm, maxbound); // В граничный узел сносятся свойства прилегающего КО.
								
								doublereal dS_loc = 0.0;

								// G - грань по которой ставится граничное условие.
								// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
								// bound_id - номер граничного узла, начиная с нуля.
								// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
								if (octree1->b4E) {
								obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
								octree1->linkE1->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
								}
								else {
								if (octree1->linkE == NULL) {
								obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
								NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
								}
								else {
								obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
								octree1->linkE->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
								}
								}

								maxbound++;
								*/

							}
							else {
								if (octree1->linkW4->inum_FD == 0) {
									G = WSIDE;
									// граничная грань:
									sosedi[G][maxelm].iNODE3 = maxelm_memo + maxbound; // граничные КО нумеруются в последнюю очередь,
									sosedi[G][maxelm].bdroblenie4 = true;
									// после maxelm_memo внутренних КО.
									// Вычисление дальнего соседа
									CALC_GG(G, GG); // Вычисляет грань GG по грани G.
									sosedi[GG][maxelm].iNODE3 = NON_EXISTENT_NODE; // соседа нет.
									sosedi[GG][maxelm].bdroblenie4 = true;

#if doubleintprecision == 1
									//printf("error inum_FD == 0 sosedi W4 elm=%lld\n", maxelm);
#else
									//printf("error inum_FD == 0 sosedi W4 elm=%d\n", maxelm);
#endif
									
									//getchar();

									// Эта ситуация возможна и здесь надо предусмотреть её обработку.
									DEFINE_BOUNDARY_PROPERTIES_FLOW(prop, prop_b, maxelm, maxbound); // В граничный узел сносятся свойства прилегающего КО.
									
									// Sy*Sz
									doublereal dS_loc = fabs(octree1->linkW4->p5.z - octree1->linkW4->p1.z)*fabs(octree1->linkW4->p2.y - octree1->linkW4->p1.y);
									TOCHKA p_centerG;
									p_centerG.z = 0.5*(octree1->linkW4->p5.z + octree1->linkW4->p1.z);
									p_centerG.y = 0.5*(octree1->linkW4->p2.y + octree1->linkW4->p1.y);
									p_centerG.x = octree1->linkW4->p1.x;


									// G - грань по которой ставится граничное условие.
									// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
									// bound_id - номер граничного узла, начиная с нуля.
									// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
									if (octree1->b4E) {
										obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
											octree1->linkE1->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
									}
									else {
										if (octree1->linkE == NULL) {
											obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
												NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
										}
										else {
											obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
												octree1->linkE->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
										}
									}

									maxbound++;
								}
								else {
									// сосед существует.
									if (bvisit[octree1->linkW4->inum_FD - 1]) {
										// узел уже был посещён
										// Здесь неким образом модифицировалось gran_t
										// Только внутренняя грань.
										G = WSIDE;
										// граничная грань:
										sosedi[G][maxelm].iNODE3 = octree1->linkW4->inum_FD - 1; // граничные КО нумеруются в последнюю очередь,
										sosedi[G][maxelm].bdroblenie4 = true;
										// после maxelm_memo внутренних КО.
										// Вычисление дальнего соседа
										CALC_GG(G, GG); // Вычисляет грань GG по грани G.

										// Определение дальнего соседа :
										// TODO 13 сентября 2016.
									}
									else {
										// узнать координаты центра грани и ориентацию в пространстве
										x_c = octree1->p0.x;
										y_c = 0.5*(octree1->linkW4->p1.y + octree1->linkW4->p2.y);
										z_c = 0.5*(octree1->linkW4->p1.z + octree1->linkW4->p5.z);
										iplane = YZ;
										patch_maxbound(iplane, s, ls, x_c, y_c, z_c, maxbound, binc3, lsid3);
										if (binc3) {
#if doubleintprecision == 1
											printf("error internal source W4 elm=%lld\n", maxelm);
#else
											printf("error internal source W4 elm=%d\n", maxelm);
#endif
											
											//getchar();
											system("PAUSE");
										}
										else {
											// Только внутренняя грань.
											G = WSIDE;
											// граничная грань:
											sosedi[G][maxelm].iNODE3 = octree1->linkW4->inum_FD - 1; // граничные КО нумеруются в последнюю очередь,
											sosedi[G][maxelm].bdroblenie4 = true;
											// после maxelm_memo внутренних КО.
											// Вычисление дальнего соседа
											CALC_GG(G, GG); // Вычисляет грань GG по грани G.

											// Определение дальнего соседа :
											// TODO 13 сентября 2016.
										}
									}
								}
							}

							if (octree1->linkW7 == NULL) {

								// 21 сентября 2016 (вырождение ячейки).

								// Заглушка.
								G = WSIDE;
								// граничная грань:
								sosedi[G][maxelm].iNODE4 = NON_EXISTENT_NODE; // граничные КО нумеруются в последнюю очередь,
								sosedi[G][maxelm].bdroblenie4 = true;
								// после maxelm_memo внутренних КО.
								// Вычисление дальнего соседа
								CALC_GG(G, GG); // Вычисляет грань GG по грани G.
								sosedi[GG][maxelm].iNODE4 = NON_EXISTENT_NODE; // соседа нет.
								sosedi[GG][maxelm].bdroblenie4 = true;

								/*
								G = WSIDE;
								// граничная грань:
								sosedi[G][maxelm].iNODE4 = maxelm_memo + maxbound; // граничные КО нумеруются в последнюю очередь,
								sosedi[G][maxelm].bdroblenie4 = true;
								// после maxelm_memo внутренних КО.
								// Вычисление дальнего соседа
								CALC_GG(G,GG); // Вычисляет грань GG по грани G.
								sosedi[GG][maxelm].iNODE4 = NON_EXISTENT_NODE; // соседа нет.
								sosedi[GG][maxelm].bdroblenie4 = true;
								#if doubleintprecision == 1
										printf("error NULL sosedi W7 elm=%lld\n", maxelm);
								#else
										printf("error NULL sosedi W7 elm=%d\n", maxelm);
								#endif
								
								getchar();

								// Эта ситуация возможна и здесь надо предусмотреть её обработку.
								 DEFINE_BOUNDARY_PROPERTIES_FLOW(prop, prop_b, maxelm, maxbound); // В граничный узел сносятся свойства прилегающего КО.
								
								doublereal dS_loc = 0.0;

								// G - грань по которой ставится граничное условие.
								// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
								// bound_id - номер граничного узла, начиная с нуля.
								// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
								if (octree1->b4E) {
								obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
								octree1->linkE1->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
								}
								else {
								if (octree1->linkE == NULL) {
								obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
								NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
								}
								else {
								obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
								octree1->linkE->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
								}
								}

								maxbound++;
								*/

							}
							else {
								if (octree1->linkW7->inum_FD == 0) {
									G = WSIDE;
									// граничная грань:
									sosedi[G][maxelm].iNODE4 = maxelm_memo + maxbound; // граничные КО нумеруются в последнюю очередь,
									sosedi[G][maxelm].bdroblenie4 = true;
									// после maxelm_memo внутренних КО.
									// Вычисление дальнего соседа
									CALC_GG(G, GG); // Вычисляет грань GG по грани G.
									sosedi[GG][maxelm].iNODE4 = NON_EXISTENT_NODE; // соседа нет.
									sosedi[GG][maxelm].bdroblenie4 = true;
#if doubleintprecision == 1
									//printf("error inum_FD == 0 sosedi W7 elm=%lld\n", maxelm);
#else
									//printf("error inum_FD == 0 sosedi W7 elm=%d\n", maxelm);
#endif
									
									//getchar();

									// Эта ситуация возможна и здесь надо предусмотреть её обработку.
									DEFINE_BOUNDARY_PROPERTIES_FLOW(prop, prop_b, maxelm, maxbound); // В граничный узел сносятся свойства прилегающего КО.
									
									// Sy*Sz
									doublereal dS_loc = fabs(octree1->linkW7->p5.z - octree1->linkW7->p1.z)*fabs(octree1->linkW7->p2.y - octree1->linkW7->p1.y);
									TOCHKA p_centerG;
									p_centerG.z = 0.5*(octree1->linkW7->p5.z + octree1->linkW7->p1.z);
									p_centerG.y = 0.5*(octree1->linkW7->p2.y + octree1->linkW7->p1.y);
									p_centerG.x = octree1->linkW7->p1.x;



									// G - грань по которой ставится граничное условие.
									// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
									// bound_id - номер граничного узла, начиная с нуля.
									// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
									if (octree1->b4E) {
										obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
											octree1->linkE1->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
									}
									else {
										if (octree1->linkE == NULL) {
											obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
												NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
										}
										else {
											obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
												octree1->linkE->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
										}
									}

									maxbound++;
								}
								else {
									// сосед существует.
									if (bvisit[octree1->linkW7->inum_FD - 1]) {
										// узел уже был посещён
										// Здесь неким образом модифицировалось gran_t
										// Только внутренняя грань.
										G = WSIDE;
										// граничная грань:
										sosedi[G][maxelm].iNODE4 = octree1->linkW7->inum_FD - 1; // граничные КО нумеруются в последнюю очередь,
										sosedi[G][maxelm].bdroblenie4 = true;
										// после maxelm_memo внутренних КО.
										// Вычисление дальнего соседа
										CALC_GG(G, GG); // Вычисляет грань GG по грани G.

										// Определение дальнего соседа :
										// TODO 13 сентября 2016.
									}
									else {
										// узнать координаты центра грани и ориентацию в пространстве
										x_c = octree1->p0.x;
										y_c = 0.5*(octree1->linkW7->p1.y + octree1->linkW7->p2.y);
										z_c = 0.5*(octree1->linkW7->p1.z + octree1->linkW7->p5.z);
										iplane = YZ;
										patch_maxbound(iplane, s, ls, x_c, y_c, z_c, maxbound, binc4, lsid4);
										if (binc4) {
#if doubleintprecision == 1
											printf("error internal source W7 elm=%lld\n", maxelm);
#else
											printf("error internal source W7 elm=%d\n", maxelm);
#endif
											
											//getchar();
											system("PAUSE");
										}
										else {
											// Только внутренняя грань.
											G = WSIDE;
											// граничная грань:
											sosedi[G][maxelm].iNODE4 = octree1->linkW7->inum_FD - 1; // граничные КО нумеруются в последнюю очередь,
											sosedi[G][maxelm].bdroblenie4 = true;
											// после maxelm_memo внутренних КО.
											// Вычисление дальнего соседа
											CALC_GG(G, GG); // Вычисляет грань GG по грани G.

											// Определение дальнего соседа :
											// TODO 13 сентября 2016.
										}
									}
								}
							}
						}
					}


					// сосед не существует
					// записываем идентификатор грани и
					// увеличиваем счётчик граней.
					if (!octree1->b4N) {

						if ((octree1->linkN == NULL)||((octree1->linkN != NULL)&&(octree1->linkN->inum_FD == 0))) {
							G = NSIDE;
							// граничная грань:
							INIT_SOSEDI(sosedi, maxelm, NSIDE, GG, false, maxelm_memo + maxbound, NON_EXISTENT_NODE, NON_EXISTENT_NODE, NON_EXISTENT_NODE, false, NON_EXISTENT_NODE, NON_EXISTENT_NODE, NON_EXISTENT_NODE, NON_EXISTENT_NODE); // заполнение sosedi[G][maxelm].
							DEFINE_BOUNDARY_PROPERTIES_FLOW(prop, prop_b, maxelm, maxbound); // В граничный узел сносятся свойства прилегающего КО.
							
							doublereal dSloc = CALC_SQUARE(G, maxelm, nvtx, pa); // Вычисление площади грани ячейки.
							TOCHKA p_centerG = CALC_CENTERG(G, maxelm, nvtx, pa); // Вычисление геометрического центра грани.
							
							// G - грань по которой ставится граничное условие.
							// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
							// bound_id - номер граничного узла, начиная с нуля.
							// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
							if (octree1->b4S) {
								obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
									octree1->linkS0->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, false, false, p_centerG);
							}
							else {
								if (octree1->linkS == NULL) {
									obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
										NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, false, false, p_centerG);
								}
								else {
									obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
										octree1->linkS->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, false, false, p_centerG);
								}
							}

							maxbound++;
						}
						else {
							if (octree1->linkN != NULL) {
								// сосед существует.
								if (bvisit[octree1->linkN->inum_FD - 1]) {

									integer return_current_node_number = sosedi[SSIDE][octree1->linkN->inum_FD - 1].iNODE1;
									//S0 NODE1,
									//S1 NODE2,
									//S4 NODE3,
									//S5 NODE4.
									 // TODO 23 сентября 2016.
									if (octree1->linkN->b4S) {
										integer current_node_number = octree1->inum_FD - 1;
										bool bcontinue = true;
										if (bcontinue) {
											if (octree1->linkN->linkS0 != NULL) {
												if (octree1->linkN->linkS0->inum_FD - 1 == current_node_number) {
													bcontinue = false; // iNODE1
													if (sosedi[SSIDE][octree1->linkN->inum_FD - 1].iNODE1 == current_node_number) {
														return_current_node_number = sosedi[SSIDE][octree1->linkN->inum_FD - 1].iNODE1;
													}
													else {
														printf("this can not be N S0\n");
														//getchar();
														system("PAUSE");
													}
												}
											}
										}
										if (bcontinue) {
											if (octree1->linkN->linkS1 != NULL) {
												if (octree1->linkN->linkS1->inum_FD - 1 == current_node_number) {
													bcontinue = false; // iNODE2
													if (sosedi[SSIDE][octree1->linkN->inum_FD - 1].iNODE2 == current_node_number) {
														return_current_node_number = sosedi[SSIDE][octree1->linkN->inum_FD - 1].iNODE2;
													}
													else {
														printf("this can not be N S1\n");
														//getchar();
														system("PAUSE");
													}
												}
											}
										}
										if (bcontinue) {
											if (octree1->linkN->linkS4 != NULL) {
												if (octree1->linkN->linkS4->inum_FD - 1 == current_node_number) {
													bcontinue = false; // iNODE3
													if (sosedi[SSIDE][octree1->linkN->inum_FD - 1].iNODE3 == current_node_number) {
														return_current_node_number = sosedi[SSIDE][octree1->linkN->inum_FD - 1].iNODE3;
													}
													else {
														printf("this can not be N S4\n");
														//getchar();
														system("PAUSE");
													}
												}
											}
										}
										if (bcontinue) {
											if (octree1->linkN->linkS5 != NULL) {
												if (octree1->linkN->linkS5->inum_FD - 1 == current_node_number) {
													bcontinue = false; // iNODE4
													if (sosedi[SSIDE][octree1->linkN->inum_FD - 1].iNODE4 == current_node_number) {
														return_current_node_number = sosedi[SSIDE][octree1->linkN->inum_FD - 1].iNODE4;
													}
													else {
														printf("this can not be N S5\n");
														//getchar();
														system("PAUSE");
													}
												}
											}
										}

									}
									

									// узел уже был посещён
									// Здесь неким образом модифицировалось gran_t
									if (return_current_node_number >= maxelm_memo) {

										printf("bvisit error  N ...\n");
										printf("model is incorrect. FATALL error!!!\n");
										printf("internal flat heat source can not be located\n");
										printf("inside the liquid.\n");
										printf("Please, press any key to exit...\n");
										//getchar();
										system("PAUSE");
										//exit(1);

										/*

										// Это граничная грань внутреннего источника тепла.
										// Внутренний источник тепла.
										G = NSIDE;
										// граничная грань:
										INIT_SOSEDI(sosedi, maxelm, NSIDE, GG, false, sosedi[SSIDE][octree1->linkN->inum_FD - 1].iNODE1, NON_EXISTENT_NODE, NON_EXISTENT_NODE, NON_EXISTENT_NODE, false, NON_EXISTENT_NODE, NON_EXISTENT_NODE, NON_EXISTENT_NODE, NON_EXISTENT_NODE); // заполнение sosedi[G][maxelm].
										
										
											// FLUID
											prop_b[RHO][sosedi[G][maxelm].iNODE1 - maxelm_memo] = prop[RHO][octree1->linkN->inum_FD - 1];
											prop_b[CP][sosedi[G][maxelm].iNODE1 - maxelm_memo] = prop[CP][octree1->linkN->inum_FD - 1];
											prop_b[LAM][sosedi[G][maxelm].iNODE1 - maxelm_memo] = prop[LAM][octree1->linkN->inum_FD - 1];
											prop_b[MULT_LAM_X][sosedi[G][maxelm].iNODE1 - maxelm_memo] = prop[MULT_LAM_X][octree1->linkN->inum_FD - 1];
											prop_b[MULT_LAM_Y][sosedi[G][maxelm].iNODE1 - maxelm_memo] = prop[MULT_LAM_Y][octree1->linkN->inum_FD - 1];
											prop_b[MULT_LAM_Z][sosedi[G][maxelm].iNODE1 - maxelm_memo] = prop[MULT_LAM_Z][octree1->linkN->inum_FD - 1];
										

										doublereal dSloc = CALC_SQUARE(G, maxelm, nvtx, pa); // Вычисление площади грани ячейки.
										TOCHKA p_centerG = CALC_CENTERG(G, maxelm, nvtx, pa); // Вычисление геометрического центра грани.
										
										// G - грань по которой ставится граничное условие.
										// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
										// bound_id - номер граничного узла, начиная с нуля.
										// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
										if (octree1->b4S) {
											obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, sosedi[SSIDE][octree1->linkN->inum_FD - 1].iNODE1 - maxelm_memo,
												octree1->linkS0->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, true, true, p_centerG);
										}
										else {
											if (octree1->linkS == NULL) {
												obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, sosedi[SSIDE][octree1->linkN->inum_FD - 1].iNODE1 - maxelm_memo,
													NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, true, true, p_centerG);
											}
											else {
												obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, sosedi[SSIDE][octree1->linkN->inum_FD - 1].iNODE1 - maxelm_memo,
													octree1->linkS->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, true, true, p_centerG);
											}
										}

										*/

									}
									else {
										// Это строго внутренний узел.
										// Внутренний узел.
										G = NSIDE;
										// внутренняя грань:
										sosedi[G][maxelm].iNODE1 = octree1->linkN->inum_FD - 1; // граничные КО нумеруются в последнюю очередь,
										sosedi[G][maxelm].bdroblenie4 = false;
										sosedi[G][maxelm].iNODE2 = NON_EXISTENT_NODE;
										sosedi[G][maxelm].iNODE3 = NON_EXISTENT_NODE;
										sosedi[G][maxelm].iNODE4 = NON_EXISTENT_NODE;
										// после maxelm_memo внутренних КО.
										// Вычисление дальнего соседа.
										CALC_GG(G, GG); // Вычисляет грань GG по грани G.
										// Определение дальнего соседа :
										// TODO 13 сентября 2016.
									}

								}
								else {
									bool binc;
									integer lsid;
									// узнать координаты центра грани и ориентацию в пространстве
									x_c = 0.5*(octree1->p2.x + octree1->p3.x);
									y_c = octree1->p3.y;
									z_c = 0.5*(octree1->p3.z + octree1->p7.z);
									iplane = XZ;
									patch_maxbound(iplane, s, ls, x_c, y_c, z_c, maxbound, binc, lsid);
									if (binc) {

										printf("model is incorrect. FATALL error!!!\n");
										printf("internal flat heat source can not be located\n");
										printf("inside the liquid.\n");
										printf("Please, press any key to exit...\n");
										//getchar();
										system("PAUSE");
										exit(1);

										/*

										// Внутренний источник тепла.
										G = NSIDE;
										// граничная грань:
										INIT_SOSEDI(sosedi, maxelm, NSIDE, GG, false, maxelm_memo + maxbound-1, NON_EXISTENT_NODE, NON_EXISTENT_NODE, NON_EXISTENT_NODE, false, NON_EXISTENT_NODE, NON_EXISTENT_NODE, NON_EXISTENT_NODE, NON_EXISTENT_NODE); // заполнение sosedi[G][maxelm].
										
										
											// FLUID
											prop_b[RHO][maxbound - 1] = prop[RHO][octree1->linkN->inum_FD - 1];
											prop_b[CP][maxbound - 1] = prop[CP][octree1->linkN->inum_FD - 1];
											prop_b[LAM][maxbound - 1] = prop[LAM][octree1->linkN->inum_FD - 1];
											prop_b[MULT_LAM_X][maxbound - 1] = prop[MULT_LAM_X][octree1->linkN->inum_FD - 1];
											prop_b[MULT_LAM_Y][maxbound - 1] = prop[MULT_LAM_Y][octree1->linkN->inum_FD - 1];
											prop_b[MULT_LAM_Z][maxbound - 1] = prop[MULT_LAM_Z][octree1->linkN->inum_FD - 1];
										

										doublereal dSloc = CALC_SQUARE(G, maxelm, nvtx, pa); // Вычисление площади грани ячейки.
										TOCHKA p_centerG = CALC_CENTERG(G, maxelm, nvtx, pa); // Вычисление геометрического центра грани.
										
										// G - грань по которой ставится граничное условие.
										// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
										// bound_id - номер граничного узла, начиная с нуля.
										// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
										if (octree1->b4S) {
											obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound - 1,
												octree1->linkS0->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, true, false, p_centerG);
										}
										else {
											if (octree1->linkS == NULL) {
												obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound - 1,
													NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, true, false, p_centerG);
											}
											else {
												obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound - 1,
													octree1->linkS->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, true, false, p_centerG);
											}
										}
                                        
										*/

									}
									else {
										// Внутренний узел.
										G = NSIDE;
										// внутренняя грань:
										sosedi[G][maxelm].iNODE1 = octree1->linkN->inum_FD - 1; // граничные КО нумеруются в последнюю очередь,
										sosedi[G][maxelm].bdroblenie4 = false;
										sosedi[G][maxelm].iNODE2 = NON_EXISTENT_NODE;
										sosedi[G][maxelm].iNODE3 = NON_EXISTENT_NODE;
										sosedi[G][maxelm].iNODE4 = NON_EXISTENT_NODE;
										// после maxelm_memo внутренних КО.
										// Вычисление дальнего соседа.
										CALC_GG(G, GG); // Вычисляет грань GG по грани G.
										// Определение дальнего соседа :
										// TODO 13 сентября 2016.

									}

								}
							}
						}
					}
					else {
					// Присутствуют все 4 соседа: N: 2,3,6,7.

						if (((octree1->linkN2 == NULL) && (octree1->linkN3 == NULL) && (octree1->linkN6 == NULL) && (octree1->linkN7 == NULL))||
							(((octree1->linkN2 != NULL) && (octree1->linkN3 != NULL) && (octree1->linkN6 != NULL) && (octree1->linkN7 != NULL) &&
							(octree1->linkN2->inum_FD == 0) && (octree1->linkN3->inum_FD == 0) && (octree1->linkN6->inum_FD == 0) && (octree1->linkN7->inum_FD == 0)))||
							(((octree1->linkN2 == NULL) && (octree1->linkN3 != NULL) && (octree1->linkN6 == NULL) && (octree1->linkN7 != NULL) &&
							(octree1->linkN3->inum_FD == 0) && (octree1->linkN7->inum_FD == 0)))||
							(((octree1->linkN2 != NULL) && (octree1->linkN3 != NULL) && (octree1->linkN6 == NULL) && (octree1->linkN7 == NULL) &&
							(octree1->linkN2->inum_FD == 0) && (octree1->linkN3->inum_FD == 0)))||
							(((octree1->linkN2 == NULL) && (octree1->linkN3 != NULL) && (octree1->linkN6 == NULL) && (octree1->linkN7 == NULL) &&
							(octree1->linkN3->inum_FD == 0))))
						{
							// 1 - N2&&N3&&N6&&N7 - на границе расчётной области (РО).
							// 2 -  N2&&N3&&N6&&N7 - на гранеице HOLLOW блока.
							// 3 - N3&&N7 на границе HOLLOW блока, N2&&N6 на границе РО.
							// 4 - N2&&N3 на границе HOLLOW блока, N6&N7 на границе РО.
							// 5 - N3 на границе HOLLOW блока, N2&&N6&&N7  на границе РО.
							
							G = NSIDE;
							// граничная грань:
							INIT_SOSEDI(sosedi, maxelm, NSIDE, GG, false, maxelm_memo + maxbound, NON_EXISTENT_NODE, NON_EXISTENT_NODE, NON_EXISTENT_NODE, false, NON_EXISTENT_NODE, NON_EXISTENT_NODE, NON_EXISTENT_NODE, NON_EXISTENT_NODE); // заполнение sosedi[G][maxelm].
							DEFINE_BOUNDARY_PROPERTIES_FLOW(prop, prop_b, maxelm, maxbound); // В граничный узел сносятся свойства прилегающего КО.
							
							doublereal dSloc = CALC_SQUARE(G, maxelm, nvtx, pa); // Вычисление площади грани ячейки.
							TOCHKA p_centerG = CALC_CENTERG(G, maxelm, nvtx, pa); // Вычисление геометрического центра грани.
							
							// G - грань по которой ставится граничное условие.
							// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
							// bound_id - номер граничного узла, начиная с нуля.
							// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
							if (octree1->b4S) {
								obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
									octree1->linkS0->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, false, false, p_centerG);
							}
							else {
								if (octree1->linkS == NULL) {
									obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
										NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, false, false, p_centerG);
								}
								else {
									obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
										octree1->linkS->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, false, false, p_centerG);
								}
							}

							// Зачем дробить если граница заканчивается NULL.
							maxbound++; // Вместо четырёх один.
						}						
						else {

							bool binc1, binc2, binc3, binc4;
							integer lsid1, lsid2, lsid3, lsid4;

							if (octree1->linkN2 == NULL) {

								// 21 сентября 2016 (вырождение ячейки).

								// Заглушка.
								G = NSIDE;
								// граничная грань:
								sosedi[G][maxelm].iNODE1 = NON_EXISTENT_NODE; // граничные КО нумеруются в последнюю очередь,
								sosedi[G][maxelm].bdroblenie4 = true;
								// после maxelm_memo внутренних КО.
								// Вычисление дальнего соседа
								CALC_GG(G, GG); // Вычисляет грань GG по грани G.

								sosedi[GG][maxelm].iNODE1 = NON_EXISTENT_NODE; // нет соседа,
								sosedi[GG][maxelm].bdroblenie4 = true;


								/*
								G = NSIDE;
								// граничная грань:
								sosedi[G][maxelm].iNODE1 = maxelm_memo + maxbound; // граничные КО нумеруются в последнюю очередь,
								sosedi[G][maxelm].bdroblenie4 = true;
								// после maxelm_memo внутренних КО.
								// Вычисление дальнего соседа
								CALC_GG(G,GG); // Вычисляет грань GG по грани G.

								sosedi[GG][maxelm].iNODE1 = NON_EXISTENT_NODE; // нет соседа,
								sosedi[GG][maxelm].bdroblenie4 = true;
								#if doubleintprecision == 1
										printf("error NULL sosedi N2 elm=%lld\n", maxelm);
								#else
										printf("error NULL sosedi N2 elm=%d\n", maxelm);
								#endif
								
								getchar();

								// Эта ситуация возможна и здесь надо предусмотреть её обработку.
								// В граничный узел сносятся свойства прилегающего КО:
								prop_b[RHO][maxbound] = prop[RHO][maxelm];
								prop_b[MU][maxbound] = prop[MU][maxelm];
								prop_b[BETA_T][maxbound] = prop[BETA_T][maxelm];


								doublereal dS_loc = 0.0;

								// G - грань по которой ставится граничное условие.
								// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
								// bound_id - номер граничного узла, начиная с нуля.
								// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
								if (octree1->b4S) {
								obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
								octree1->linkS0->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
								}
								else {
								if (octree1->linkS == NULL) {
								obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
								NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
								}
								else {
								obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
								octree1->linkS->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
								}
								}


								maxbound++;
								*/
							}
							else {
								if (octree1->linkN2->inum_FD == 0) {

									G = NSIDE;
									// граничная грань:
									sosedi[G][maxelm].iNODE1 = maxelm_memo + maxbound; // граничные КО нумеруются в последнюю очередь,
									sosedi[G][maxelm].bdroblenie4 = true;
									// после maxelm_memo внутренних КО.
									// Вычисление дальнего соседа
									CALC_GG(G, GG); // Вычисляет грань GG по грани G.

									sosedi[GG][maxelm].iNODE1 = NON_EXISTENT_NODE; // нет соседа,
									sosedi[GG][maxelm].bdroblenie4 = true;

#if doubleintprecision == 1
									//printf("error inum_FD == 0 sosedi N2 elm=%lld\n", maxelm);
#else
									//printf("error inum_FD == 0 sosedi N2 elm=%d\n", maxelm);
#endif
									
									//getchar();

									// Эта ситуация возможна и здесь надо предусмотреть её обработку.
									DEFINE_BOUNDARY_PROPERTIES_FLOW(prop, prop_b, maxelm, maxbound); // В граничный узел сносятся свойства прилегающего КО.
									
									// Sx*Sz
									doublereal dS_loc = fabs(octree1->linkN2->p4.z - octree1->linkN2->p0.z)*fabs(octree1->linkN2->p1.x - octree1->linkN2->p0.x);
									TOCHKA p_centerG;
									p_centerG.z = 0.5*(octree1->linkN2->p4.z + octree1->linkN2->p0.z);
									p_centerG.y = octree1->linkN2->p1.y;
									p_centerG.x = 0.5*(octree1->linkN2->p1.x + octree1->linkN2->p0.x);



									// G - грань по которой ставится граничное условие.
									// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
									// bound_id - номер граничного узла, начиная с нуля.
									// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
									if (octree1->b4S) {
										obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
											octree1->linkS0->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
									}
									else {
										if (octree1->linkS == NULL) {
											obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
												NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
										}
										else {
											obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
												octree1->linkS->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
										}
									}

									maxbound++;
								}
								else {
									// сосед существует.
									if (bvisit[octree1->linkN2->inum_FD - 1]) {
										// узел уже был посещён
										// Здесь неким образом модифицировалось gran_t
										// Только внутренняя грань.
										G = NSIDE;
										// граничная грань:
										sosedi[G][maxelm].iNODE1 = octree1->linkN2->inum_FD - 1; // граничные КО нумеруются в последнюю очередь,
										sosedi[G][maxelm].bdroblenie4 = true;
										// после maxelm_memo внутренних КО.
										// Вычисление дальнего соседа
										CALC_GG(G, GG); // Вычисляет грань GG по грани G.

										// Определение дальнего соседа :
										// TODO 13 сентября 2016.
									}
									else {
										// узнать координаты центра грани и ориентацию в пространстве
										x_c = 0.5*(octree1->linkN2->p0.x + octree1->linkN2->p1.x);
										y_c = octree1->p3.y;
										z_c = 0.5*(octree1->linkN2->p0.z + octree1->linkN2->p4.z);
										iplane = XZ;
										patch_maxbound(iplane, s, ls, x_c, y_c, z_c, maxbound, binc1, lsid1);
										if (binc1) {
#if doubleintprecision == 1
											printf("error internal source N2 elm=%lld\n", maxelm);
#else
											printf("error internal source N2 elm=%d\n", maxelm);
#endif
											
											//getchar();
											system("PAUSE");
										}
										else {
											// Только внутренняя грань.
											G = NSIDE;
											// граничная грань:
											sosedi[G][maxelm].iNODE1 = octree1->linkN2->inum_FD - 1; // граничные КО нумеруются в последнюю очередь,
											sosedi[G][maxelm].bdroblenie4 = true;
											// после maxelm_memo внутренних КО.
											// Вычисление дальнего соседа
											CALC_GG(G, GG); // Вычисляет грань GG по грани G.

											// Определение дальнего соседа :
											// TODO 13 сентября 2016.
										}
									}
								}
							}

							if (octree1->linkN3 == NULL) {

#if doubleintprecision == 1
								printf("error NULL sosedi N3 elm=%lld\n", maxelm);
#else
								printf("error NULL sosedi N3 elm=%d\n", maxelm);
#endif
								
								//getchar();
								system("PAUSE");

								/*
								G = NSIDE;
								// граничная грань:
								sosedi[G][maxelm].iNODE2 = maxelm_memo + maxbound; // граничные КО нумеруются в последнюю очередь,
								sosedi[G][maxelm].bdroblenie4 = true;
								// после maxelm_memo внутренних КО.
								// Вычисление дальнего соседа
								CALC_GG(G,GG); // Вычисляет грань GG по грани G.

								sosedi[GG][maxelm].iNODE2 = NON_EXISTENT_NODE; // нет соседей,
								sosedi[GG][maxelm].bdroblenie4 = true;

								#if doubleintprecision == 1
										printf("error NULL sosedi N3 elm=%lld\n", maxelm);
								#else
										printf("error NULL sosedi N3 elm=%d\n", maxelm);
								#endif
								
								getchar();

								// Эта ситуация возможна и здесь надо предусмотреть её обработку.
								DEFINE_BOUNDARY_PROPERTIES_FLOW(prop, prop_b, maxelm, maxbound); // В граничный узел сносятся свойства прилегающего КО.
								
								doublereal dS_loc = 0.0;

								// G - грань по которой ставится граничное условие.
								// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
								// bound_id - номер граничного узла, начиная с нуля.
								// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
								if (octree1->b4S) {
								obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
								octree1->linkS0->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
								}
								else {
								if (octree1->linkS == NULL) {
								obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
								NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
								}
								else {
								obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
								octree1->linkS->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
								}
								}


								maxbound++;
								*/

							}
							else {
								if (octree1->linkN3->inum_FD == 0) {

									G = NSIDE;
									// граничная грань:
									sosedi[G][maxelm].iNODE2 = maxelm_memo + maxbound; // граничные КО нумеруются в последнюю очередь,
									sosedi[G][maxelm].bdroblenie4 = true;
									// после maxelm_memo внутренних КО.
									// Вычисление дальнего соседа
									CALC_GG(G, GG); // Вычисляет грань GG по грани G.

									sosedi[GG][maxelm].iNODE2 = NON_EXISTENT_NODE; // нет соседей,
									sosedi[GG][maxelm].bdroblenie4 = true;

#if doubleintprecision == 1
									//printf("error inum_FD == 0 sosedi N3 elm=%lld\n", maxelm);
#else
									//printf("error inum_FD == 0 sosedi N3 elm=%d\n", maxelm);
#endif
									
									//getchar();

									// Эта ситуация возможна и здесь надо предусмотреть её обработку.
									DEFINE_BOUNDARY_PROPERTIES_FLOW(prop, prop_b, maxelm, maxbound); // В граничный узел сносятся свойства прилегающего КО.
									
									// Sx*Sz
									doublereal dS_loc = fabs(octree1->linkN3->p4.z - octree1->linkN3->p0.z)*fabs(octree1->linkN3->p1.x - octree1->linkN3->p0.x);
									TOCHKA p_centerG;
									p_centerG.z = 0.5*(octree1->linkN3->p4.z + octree1->linkN3->p0.z);
									p_centerG.y = octree1->linkN3->p1.y;
									p_centerG.x = 0.5*(octree1->linkN3->p1.x + octree1->linkN3->p0.x);

									// G - грань по которой ставится граничное условие.
									// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
									// bound_id - номер граничного узла, начиная с нуля.
									// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
									if (octree1->b4S) {
										obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
											octree1->linkS0->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
									}
									else {
										if (octree1->linkS == NULL) {
											obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
												NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
										}
										else {
											obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
												octree1->linkS->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
										}
									}


									maxbound++;
								}
								else {
									// сосед существует.
									if (bvisit[octree1->linkN3->inum_FD - 1]) {
										// узел уже был посещён
										// Здесь неким образом модифицировалось gran_t
										// Только внутренняя грань.
										G = NSIDE;
										// граничная грань:
										sosedi[G][maxelm].iNODE2 = octree1->linkN3->inum_FD - 1; // граничные КО нумеруются в последнюю очередь,
										sosedi[G][maxelm].bdroblenie4 = true;
										// после maxelm_memo внутренних КО.
										// Вычисление дальнего соседа
										CALC_GG(G, GG); // Вычисляет грань GG по грани G.

										// Определение дальнего соседа :
										// TODO 13 сентября 2016.
									}
									else {
										// узнать координаты центра грани и ориентацию в пространстве
										x_c = 0.5*(octree1->linkN3->p0.x + octree1->linkN3->p1.x);
										y_c = octree1->p3.y;
										z_c = 0.5*(octree1->linkN3->p0.z + octree1->linkN3->p4.z);
										iplane = XZ;
										patch_maxbound(iplane, s, ls, x_c, y_c, z_c, maxbound, binc2, lsid2);
										if (binc2) {
#if doubleintprecision == 1
											printf("error internal source N3 elm=%lld\n", maxelm);
#else
											printf("error internal source N3 elm=%d\n", maxelm);
#endif
											
											//getchar();
											system("PAUSE");
										}
										else {
											// Только внутренняя грань.
											G = NSIDE;
											// граничная грань:
											sosedi[G][maxelm].iNODE2 = octree1->linkN3->inum_FD - 1; // граничные КО нумеруются в последнюю очередь,
											sosedi[G][maxelm].bdroblenie4 = true;
											// после maxelm_memo внутренних КО.
											// Вычисление дальнего соседа
											CALC_GG(G, GG); // Вычисляет грань GG по грани G.

											// Определение дальнего соседа :
											// TODO 13 сентября 2016.
										}
									}
								}
							}

							if (octree1->linkN6 == NULL) {

								// 21 сентября 2016 (вырождение ячейки).

								// Заглушка.
								G = NSIDE;
								// граничная грань:
								sosedi[G][maxelm].iNODE3 = NON_EXISTENT_NODE; // граничные КО нумеруются в последнюю очередь,
								sosedi[G][maxelm].bdroblenie4 = true;
								// после maxelm_memo внутренних КО.
								// Вычисление дальнего соседа
								CALC_GG(G, GG); // Вычисляет грань GG по грани G.

								sosedi[GG][maxelm].iNODE3 = NON_EXISTENT_NODE; // сосед отсутствует,
								sosedi[GG][maxelm].bdroblenie4 = true;

								/*
								G = NSIDE;
								// граничная грань:
								sosedi[G][maxelm].iNODE3 = maxelm_memo + maxbound; // граничные КО нумеруются в последнюю очередь,
								sosedi[G][maxelm].bdroblenie4 = true;
								// после maxelm_memo внутренних КО.
								// Вычисление дальнего соседа
								CALC_GG(G,GG); // Вычисляет грань GG по грани G.

								sosedi[GG][maxelm].iNODE3 = NON_EXISTENT_NODE; // сосед отсутствует,
								sosedi[GG][maxelm].bdroblenie4 = true;

								#if doubleintprecision == 1
										printf("error NULL sosedi N6 elm=%lld\n", maxelm);
								#else
										printf("error NULL sosedi N6 elm=%d\n", maxelm);
								#endif
								
								getchar();

								// Эта ситуация возможна и здесь надо предусмотреть её обработку.
								DEFINE_BOUNDARY_PROPERTIES_FLOW(prop, prop_b, maxelm, maxbound); // В граничный узел сносятся свойства прилегающего КО.
								
								doublereal dS_loc = 0.0;

								// G - грань по которой ставится граничное условие.
								// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
								// bound_id - номер граничного узла, начиная с нуля.
								// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
								if (octree1->b4S) {
								obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
								octree1->linkS0->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
								}
								else {
								if (octree1->linkS == NULL) {
								obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
								NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
								}
								else {
								obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
								octree1->linkS->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
								}
								}

								maxbound++;
								*/

							}
							else {
								if (octree1->linkN6->inum_FD == 0) {

									G = NSIDE;
									// граничная грань:
									sosedi[G][maxelm].iNODE3 = maxelm_memo + maxbound; // граничные КО нумеруются в последнюю очередь,
									sosedi[G][maxelm].bdroblenie4 = true;
									// после maxelm_memo внутренних КО.
									// Вычисление дальнего соседа
									CALC_GG(G, GG); // Вычисляет грань GG по грани G.

									sosedi[GG][maxelm].iNODE3 = NON_EXISTENT_NODE; // сосед отсутствует,
									sosedi[GG][maxelm].bdroblenie4 = true;

#if doubleintprecision == 1
									//printf("error inum_FD == 0 sosedi N6 elm=%lld\n", maxelm);
#else
									//printf("error inum_FD == 0 sosedi N6 elm=%d\n", maxelm);
#endif
									
									//getchar();

									// Эта ситуация возможна и здесь надо предусмотреть её обработку.
									DEFINE_BOUNDARY_PROPERTIES_FLOW(prop, prop_b, maxelm, maxbound); // В граничный узел сносятся свойства прилегающего КО.
									
									// Sx*Sz
									doublereal dS_loc = fabs(octree1->linkN6->p4.z - octree1->linkN6->p0.z)*fabs(octree1->linkN6->p1.x - octree1->linkN6->p0.x);
									TOCHKA p_centerG;
									p_centerG.z = 0.5*(octree1->linkN6->p4.z + octree1->linkN6->p0.z);
									p_centerG.y = octree1->linkN6->p1.y;
									p_centerG.x = 0.5*(octree1->linkN6->p1.x + octree1->linkN6->p0.x);

									// G - грань по которой ставится граничное условие.
									// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
									// bound_id - номер граничного узла, начиная с нуля.
									// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
									if (octree1->b4S) {
										obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
											octree1->linkS0->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
									}
									else {
										if (octree1->linkS == NULL) {
											obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
												NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
										}
										else {
											obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
												octree1->linkS->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
										}
									}

									maxbound++;
								}
								else {
									// сосед существует.
									if (bvisit[octree1->linkN6->inum_FD - 1]) {
										// узел уже был посещён
										// Здесь неким образом модифицировалось gran_t
										// Только внутренняя грань.
										G = NSIDE;
										// граничная грань:
										sosedi[G][maxelm].iNODE3 = octree1->linkN6->inum_FD - 1; // граничные КО нумеруются в последнюю очередь,
										sosedi[G][maxelm].bdroblenie4 = true;
										// после maxelm_memo внутренних КО.
										// Вычисление дальнего соседа
										CALC_GG(G, GG); // Вычисляет грань GG по грани G.

										// Определение дальнего соседа :
										// TODO 13 сентября 2016.
									}
									else {
										// узнать координаты центра грани и ориентацию в пространстве
										x_c = 0.5*(octree1->linkN6->p0.x + octree1->linkN6->p1.x);
										y_c = octree1->p3.y;
										z_c = 0.5*(octree1->linkN6->p0.z + octree1->linkN6->p4.z);
										iplane = XZ;
										patch_maxbound(iplane, s, ls, x_c, y_c, z_c, maxbound, binc3, lsid3);
										if (binc3) {
#if doubleintprecision == 1
											printf("error internal source N6 elm=%lld\n", maxelm);
#else
											printf("error internal source N6 elm=%d\n", maxelm);
#endif
											
											//getchar();
											system("PAUSE");
										}
										else {
											// Только внутренняя грань.
											G = NSIDE;
											// граничная грань:
											sosedi[G][maxelm].iNODE3 = octree1->linkN6->inum_FD - 1; // граничные КО нумеруются в последнюю очередь,
											sosedi[G][maxelm].bdroblenie4 = true;
											// после maxelm_memo внутренних КО.
											// Вычисление дальнего соседа
											CALC_GG(G, GG); // Вычисляет грань GG по грани G.

											// Определение дальнего соседа :
											// TODO 13 сентября 2016.
										}
									}
								}
							}

							if (octree1->linkN7 == NULL) {

								// 21 сентября 2016 (вырождение ячейки).

								// Заглушка.
								G = NSIDE;
								// граничная грань:
								sosedi[G][maxelm].iNODE4 = NON_EXISTENT_NODE; // граничные КО нумеруются в последнюю очередь,
								sosedi[G][maxelm].bdroblenie4 = true;
								// после maxelm_memo внутренних КО.
								// Вычисление дальнего соседа
								CALC_GG(G, GG); // Вычисляет грань GG по грани G.

								sosedi[GG][maxelm].iNODE4 = NON_EXISTENT_NODE; // граничные КО нумеруются в последнюю очередь,
								sosedi[GG][maxelm].bdroblenie4 = true;


								/*
								G = NSIDE;
								// граничная грань:
								sosedi[G][maxelm].iNODE4 = maxelm_memo + maxbound; // граничные КО нумеруются в последнюю очередь,
								sosedi[G][maxelm].bdroblenie4 = true;
								// после maxelm_memo внутренних КО.
								// Вычисление дальнего соседа
								CALC_GG(G,GG); // Вычисляет грань GG по грани G.

								sosedi[GG][maxelm].iNODE4 = - 1; // граничные КО нумеруются в последнюю очередь,
								sosedi[GG][maxelm].bdroblenie4 = true;

								#if doubleintprecision == 1
										printf("error NULL sosedi N7 elm=%lld\n", maxelm);
								#else
										printf("error NULL sosedi N7 elm=%d\n", maxelm);
								#endif
								
								getchar();

								// Эта ситуация возможна и здесь надо предусмотреть её обработку.
								DEFINE_BOUNDARY_PROPERTIES_FLOW(prop, prop_b, maxelm, maxbound); // В граничный узел сносятся свойства прилегающего КО.
								
								doublereal dS_loc = 0.0;

								// G - грань по которой ставится граничное условие.
								// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
								// bound_id - номер граничного узла, начиная с нуля.
								// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
								if (octree1->b4S) {
								obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
								octree1->linkS0->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
								}
								else {
								if (octree1->linkS == NULL) {
								obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
								NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
								}
								else {
								obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
								octree1->linkS->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
								}
								}


								maxbound++;
								*/
							}
							else {
								if (octree1->linkN7->inum_FD == 0) {

									G = NSIDE;
									// граничная грань:
									sosedi[G][maxelm].iNODE4 = maxelm_memo + maxbound; // граничные КО нумеруются в последнюю очередь,
									sosedi[G][maxelm].bdroblenie4 = true;
									// после maxelm_memo внутренних КО.
									// Вычисление дальнего соседа
									CALC_GG(G, GG); // Вычисляет грань GG по грани G.

									sosedi[GG][maxelm].iNODE4 = NON_EXISTENT_NODE; // граничные КО нумеруются в последнюю очередь,
									sosedi[GG][maxelm].bdroblenie4 = true;

#if doubleintprecision == 1
									//printf("error inum_FD == 0 sosedi N7 elm=%lld\n", maxelm);
#else
									//printf("error inum_FD == 0 sosedi N7 elm=%d\n", maxelm);
#endif
									
									//getchar();

									// Эта ситуация возможна и здесь надо предусмотреть её обработку.
									DEFINE_BOUNDARY_PROPERTIES_FLOW(prop, prop_b, maxelm, maxbound); // В граничный узел сносятся свойства прилегающего КО.
									
									// Sx*Sz
									doublereal dS_loc = fabs(octree1->linkN7->p4.z - octree1->linkN7->p0.z)*fabs(octree1->linkN7->p1.x - octree1->linkN7->p0.x);
									TOCHKA p_centerG;
									p_centerG.z = 0.5*(octree1->linkN7->p4.z + octree1->linkN7->p0.z);
									p_centerG.y = octree1->linkN7->p1.y;
									p_centerG.x = 0.5*(octree1->linkN7->p1.x + octree1->linkN7->p0.x);

									// G - грань по которой ставится граничное условие.
									// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
									// bound_id - номер граничного узла, начиная с нуля.
									// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
									if (octree1->b4S) {
										obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
											octree1->linkS0->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
									}
									else {
										if (octree1->linkS == NULL) {
											obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
												NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
										}
										else {
											obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
												octree1->linkS->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
										}
									}

									maxbound++;
								}
								else {
									// сосед существует.
									if (bvisit[octree1->linkN7->inum_FD - 1]) {
										// узел уже был посещён
										// Здесь неким образом модифицировалось gran_t
										// Только внутренняя грань.
										G = NSIDE;
										// граничная грань:
										sosedi[G][maxelm].iNODE4 = octree1->linkN7->inum_FD - 1; // граничные КО нумеруются в последнюю очередь,
										sosedi[G][maxelm].bdroblenie4 = true;
										// после maxelm_memo внутренних КО.
										// Вычисление дальнего соседа
										CALC_GG(G, GG); // Вычисляет грань GG по грани G.

										// Определение дальнего соседа :
										// TODO 13 сентября 2016.
									}
									else {
										// узнать координаты центра грани и ориентацию в пространстве
										x_c = 0.5*(octree1->linkN7->p0.x + octree1->linkN7->p1.x);
										y_c = octree1->p3.y;
										z_c = 0.5*(octree1->linkN7->p0.z + octree1->linkN7->p4.z);
										iplane = XZ;
										patch_maxbound(iplane, s, ls, x_c, y_c, z_c, maxbound, binc4, lsid4);
										if (binc4) {
#if doubleintprecision == 1
											printf("error internal source N7 elm=%lld\n", maxelm);
#else
											printf("error internal source N7 elm=%d\n", maxelm);
#endif
											
											//getchar();
											system("PAUSE");
										}
										else {
											// Только внутренняя грань.
											G = NSIDE;
											// граничная грань:
											sosedi[G][maxelm].iNODE4 = octree1->linkN7->inum_FD - 1; // граничные КО нумеруются в последнюю очередь,
											sosedi[G][maxelm].bdroblenie4 = true;
											// после maxelm_memo внутренних КО.
											// Вычисление дальнего соседа
											CALC_GG(G, GG); // Вычисляет грань GG по грани G.

											// Определение дальнего соседа :
											// TODO 13 сентября 2016.
										}
									}
								}
							}
						}
					}

					// сосед не существует
					// записываем идентификатор грани и
					// увеличиваем счётчик граней.
					if (!octree1->b4S) {
						if ((octree1->linkS == NULL)||((octree1->linkS != NULL)&&(octree1->linkS->inum_FD == 0))) {
							G = SSIDE;
							// граничная грань:
							INIT_SOSEDI(sosedi, maxelm, SSIDE, GG, false, maxelm_memo + maxbound, NON_EXISTENT_NODE, NON_EXISTENT_NODE, NON_EXISTENT_NODE, false, NON_EXISTENT_NODE, NON_EXISTENT_NODE, NON_EXISTENT_NODE, NON_EXISTENT_NODE); // заполнение sosedi[G][maxelm].
							DEFINE_BOUNDARY_PROPERTIES_FLOW(prop, prop_b, maxelm, maxbound); // В граничный узел сносятся свойства прилегающего КО.

							doublereal dSloc = CALC_SQUARE(G, maxelm, nvtx, pa); // Вычисление площади грани ячейки.
							TOCHKA p_centerG = CALC_CENTERG(G, maxelm, nvtx, pa); // Вычисление геометрического центра грани.
							
							// G - грань по которой ставится граничное условие.
							// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
							// bound_id - номер граничного узла, начиная с нуля.
							// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
							if (octree1->b4N) {
								obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
									octree1->linkN3->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, false, false, p_centerG);
							}
							else {
								if (octree1->linkN == NULL) {
									obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
										NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, false, false, p_centerG);
								}
								else {
									obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
										octree1->linkN->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, false, false, p_centerG);
								}
							}

							maxbound++;
						}
						else {
							if (octree1->linkS != NULL) {
								//bool binc;
								//integer lsid;
								// сосед существует.
								if (bvisit[octree1->linkS->inum_FD - 1]) {

									integer return_current_node_number = sosedi[NSIDE][octree1->linkS->inum_FD - 1].iNODE1;
									//N2 NODE1,
									//N3 NODE2,
									//N6 NODE3,
									//N7 NODE4.
									// TODO 23 сентября 2016.
									if (octree1->linkS->b4N) {
										integer current_node_number = octree1->inum_FD - 1;
										bool bcontinue = true;
										if (bcontinue) {
											if (octree1->linkS->linkN2 != NULL) {
												if (octree1->linkS->linkN2->inum_FD - 1 == current_node_number) {
													bcontinue = false; // iNODE1
													if (sosedi[NSIDE][octree1->linkS->inum_FD - 1].iNODE1 == current_node_number) {
														return_current_node_number = sosedi[NSIDE][octree1->linkS->inum_FD - 1].iNODE1;
													}
													else {
														printf("this can not be S N2\n");
														//getchar();
														system("PAUSE");
													}
												}
											}
										}
										if (bcontinue) {
											if (octree1->linkS->linkN3 != NULL) {
												if (octree1->linkS->linkN3->inum_FD - 1 == current_node_number) {
													bcontinue = false; // iNODE2
													if (sosedi[NSIDE][octree1->linkS->inum_FD - 1].iNODE2 == current_node_number) {
														return_current_node_number = sosedi[NSIDE][octree1->linkS->inum_FD - 1].iNODE2;
													}
													else {
														printf("this can not be S N3\n");
														//getchar();
														system("PAUSE");
													}
												}
											}
										}
										if (bcontinue) {
											if (octree1->linkS->linkN6 != NULL) {
												if (octree1->linkS->linkN6->inum_FD - 1 == current_node_number) {
													bcontinue = false; // iNODE3
													if (sosedi[NSIDE][octree1->linkS->inum_FD - 1].iNODE3 == current_node_number) {
														return_current_node_number = sosedi[NSIDE][octree1->linkS->inum_FD - 1].iNODE3;
													}
													else {
														printf("this can not be S N6\n");
														//getchar();
														system("PAUSE");
													}
												}
											}
										}
										if (bcontinue) {
											if (octree1->linkS->linkN7 != NULL) {
												if (octree1->linkS->linkN7->inum_FD - 1 == current_node_number) {
													bcontinue = false; // iNODE4
													if (sosedi[NSIDE][octree1->linkS->inum_FD - 1].iNODE4 == current_node_number) {
														return_current_node_number = sosedi[NSIDE][octree1->linkS->inum_FD - 1].iNODE4;
													}
													else {
														printf("this can not be S N7\n");
														//getchar();
														system("PAUSE");
													}
												}
											}
										}

									}


									// узел уже был посещён
									// Здесь неким образом модифицировалось gran_t
									if (return_current_node_number >= maxelm_memo) {

										printf("bvisit error S ...\n");
										printf("model is incorrect. FATALL error!!!\n");
										printf("internal flat heat source can not be located\n");
										printf("inside the liquid.\n");
										printf("Please, press any key to exit...\n");
										//getchar();
										system("PAUSE");
										//exit(1);


									    /*

										// Это граничная грань внутреннего источника тепла.
										// Внутренний источник тепла.
										G = SSIDE;
										// граничная грань:
										INIT_SOSEDI(sosedi, maxelm, SSIDE, GG, false, sosedi[NSIDE][octree1->linkS->inum_FD - 1].iNODE1, NON_EXISTENT_NODE, NON_EXISTENT_NODE, NON_EXISTENT_NODE, false, NON_EXISTENT_NODE, NON_EXISTENT_NODE, NON_EXISTENT_NODE, NON_EXISTENT_NODE); // заполнение sosedi[G][maxelm].
										
									
											// FLUID
											prop_b[RHO][sosedi[G][maxelm].iNODE1 - maxelm_memo] = prop[RHO][octree1->linkS->inum_FD - 1];
											prop_b[CP][sosedi[G][maxelm].iNODE1 - maxelm_memo] = prop[CP][octree1->linkS->inum_FD - 1];
											prop_b[LAM][sosedi[G][maxelm].iNODE1 - maxelm_memo] = prop[LAM][octree1->linkS->inum_FD - 1];
											prop_b[MULT_LAM_X][sosedi[G][maxelm].iNODE1 - maxelm_memo] = prop[MULT_LAM_X][octree1->linkS->inum_FD - 1];
											prop_b[MULT_LAM_Y][sosedi[G][maxelm].iNODE1 - maxelm_memo] = prop[MULT_LAM_Y][octree1->linkS->inum_FD - 1];
											prop_b[MULT_LAM_Z][sosedi[G][maxelm].iNODE1 - maxelm_memo] = prop[MULT_LAM_Z][octree1->linkS->inum_FD - 1];
										

										doublereal dSloc = CALC_SQUARE(G, maxelm, nvtx, pa); // Вычисление площади грани ячейки.
										TOCHKA p_centerG = CALC_CENTERG(G, maxelm, nvtx, pa); // Вычисление геометрического центра грани.
										
										// G - грань по которой ставится граничное условие.
										// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
										// bound_id - номер граничного узла, начиная с нуля.
										// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
										if (octree1->b4N) {
											obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, sosedi[NSIDE][octree1->linkS->inum_FD - 1].iNODE1 - maxelm_memo,
												octree1->linkN3->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc,true, true);
										}
										else {
											if (octree1->linkN == NULL) {
												obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, sosedi[NSIDE][octree1->linkS->inum_FD - 1].iNODE1 - maxelm_memo,
													NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc,true, true);
											}
											else {
												obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, sosedi[NSIDE][octree1->linkS->inum_FD - 1].iNODE1 - maxelm_memo,
													octree1->linkN->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, true, true, p_centerG);
											}
										}

										*/

									}
									else {
										// Это строго внутренний узел.
										// Внутренний узел.
										G = SSIDE;
										// внутренняя грань:
										sosedi[G][maxelm].iNODE1 = octree1->linkS->inum_FD - 1; // граничные КО нумеруются в последнюю очередь,
										sosedi[G][maxelm].bdroblenie4 = false;
										sosedi[G][maxelm].iNODE2 = NON_EXISTENT_NODE;
										sosedi[G][maxelm].iNODE3 = NON_EXISTENT_NODE;
										sosedi[G][maxelm].iNODE4 = NON_EXISTENT_NODE;
										// после maxelm_memo внутренних КО.
										// Вычисление дальнего соседа
										CALC_GG(G, GG); // Вычисляет грань GG по грани G.

										// Определение дальнего соседа :
										// TODO 13 сентября 2016.
									}

								}
								else {
									bool binc;
									integer lsid;
									// узнать координаты центра грани и ориентацию в пространстве
									x_c = 0.5*(octree1->p0.x + octree1->p1.x);
									y_c = octree1->p0.y;
									z_c = 0.5*(octree1->p0.z + octree1->p4.z);
									iplane = XZ;
									patch_maxbound(iplane, s, ls, x_c, y_c, z_c, maxbound, binc, lsid);
									if (binc) {

										printf("model is incorrect. FATALL error!!!\n");
										printf("internal flat heat source can not be located\n");
										printf("inside the liquid.\n");
										printf("Please, press any key to exit...\n");
										//getchar();
										system("PAUSE");
										exit(1);

										/*

										// Внутренний источник тепла.
										G = SSIDE;
										// граничная грань:
										INIT_SOSEDI(sosedi, maxelm, SSIDE, GG, false, maxelm_memo + maxbound-1, NON_EXISTENT_NODE, NON_EXISTENT_NODE, NON_EXISTENT_NODE, false, NON_EXISTENT_NODE, NON_EXISTENT_NODE, NON_EXISTENT_NODE, NON_EXISTENT_NODE); // заполнение sosedi[G][maxelm].
										
										
											// FLUID
											prop_b[RHO][maxbound - 1] = prop[RHO][octree1->linkS->inum_FD - 1];
											prop_b[CP][maxbound - 1] = prop[CP][octree1->linkS->inum_FD - 1];
											prop_b[LAM][maxbound - 1] = prop[LAM][octree1->linkS->inum_FD - 1];
											prop_b[MULT_LAM_X][maxbound - 1] = prop[MULT_LAM_X][octree1->linkS->inum_FD - 1];
											prop_b[MULT_LAM_Y][maxbound - 1] = prop[MULT_LAM_Y][octree1->linkS->inum_FD - 1];
											prop_b[MULT_LAM_Z][maxbound - 1] = prop[MULT_LAM_Z][octree1->linkS->inum_FD - 1];
										

										doublereal dSloc = CALC_SQUARE(G, maxelm, nvtx, pa); // Вычисление площади грани ячейки.
										TOCHKA p_centerG = CALC_CENTERG(G, maxelm, nvtx, pa); // Вычисление геометрического центра грани.
										
										// G - грань по которой ставится граничное условие.
										// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
										// bound_id - номер граничного узла, начиная с нуля.
										// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
										if (octree1->b4N) {
											obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound - 1,
												octree1->linkN3->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, true, false, p_centerG);
										}
										else {
											if (octree1->linkN == NULL) {
												obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound - 1,
													NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, true, false, p_centerG);
											}
											else {
												obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound - 1,
													octree1->linkN->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, true, false, p_centerG);
											}
										}

										*/

									}
									else {
										// Внутренний узел.
										G = SSIDE;
										// внутренняя грань:
										sosedi[G][maxelm].iNODE1 = octree1->linkS->inum_FD - 1; // граничные КО нумеруются в последнюю очередь,
										sosedi[G][maxelm].bdroblenie4 = false;
										sosedi[G][maxelm].iNODE2 = NON_EXISTENT_NODE;
										sosedi[G][maxelm].iNODE3 = NON_EXISTENT_NODE;
										sosedi[G][maxelm].iNODE4 = NON_EXISTENT_NODE;
										// после maxelm_memo внутренних КО.
										// Вычисление дальнего соседа
										CALC_GG(G, GG); // Вычисляет грань GG по грани G.

										// Определение дальнего соседа :
										// TODO 13 сентября 2016.
									}
								}
							}
						}
					}
					else {
					// Присутствуют все 4 соседа: S: 0,1,4,5.

						if (((octree1->linkS0 == NULL) && (octree1->linkS1 == NULL) && (octree1->linkS4 == NULL) && (octree1->linkS5 == NULL))||
							(((octree1->linkS0 != NULL) && (octree1->linkS1 != NULL) && (octree1->linkS4 != NULL) && (octree1->linkS5 != NULL) &&
							(octree1->linkS0->inum_FD == 0) && (octree1->linkS1->inum_FD == 0) && (octree1->linkS4->inum_FD == 0) && (octree1->linkS5->inum_FD == 0)))||
							(((octree1->linkS0 != NULL) && (octree1->linkS1 == NULL) && (octree1->linkS4 != NULL) && (octree1->linkS5 == NULL) &&
							(octree1->linkS0->inum_FD == 0) && (octree1->linkS4->inum_FD == 0)))||
							(((octree1->linkS0 != NULL) && (octree1->linkS1 != NULL) && (octree1->linkS4 == NULL) && (octree1->linkS5 == NULL) &&
							(octree1->linkS0->inum_FD == 0) && (octree1->linkS1->inum_FD == 0)))||
								(((octree1->linkS0 != NULL) && (octree1->linkS1 == NULL) && (octree1->linkS4 == NULL) && (octree1->linkS5 == NULL) &&
							(octree1->linkS0->inum_FD == 0))))
						{
							// 1 - S0&&S1&&S4&&S5 - на границе расчётной области (РО).
							// 2 - S0&&S1&&S4&&S5 - на гранеице HOLLOW блока.
							// 3 - S0&&S4 на границе HOLLOW блока, S1&&S5 на границе РО.
							// 4 - S0&&S1 на границе HOLLOW блока, S4&S5 на границе РО.
							// 5 - S0 на границе HOLLOW блока, S1&&S4&&S5  на границе РО.

							G = SSIDE;
							// граничная грань:
							INIT_SOSEDI(sosedi, maxelm, SSIDE, GG, false, maxelm_memo + maxbound, NON_EXISTENT_NODE, NON_EXISTENT_NODE, NON_EXISTENT_NODE, false, NON_EXISTENT_NODE, NON_EXISTENT_NODE, NON_EXISTENT_NODE, NON_EXISTENT_NODE); // заполнение sosedi[G][maxelm].
							DEFINE_BOUNDARY_PROPERTIES_FLOW(prop, prop_b, maxelm, maxbound); // В граничный узел сносятся свойства прилегающего КО.

							doublereal dSloc = CALC_SQUARE(G, maxelm, nvtx, pa); // Вычисление площади грани ячейки.
							TOCHKA p_centerG = CALC_CENTERG(G, maxelm, nvtx, pa); // Вычисление геометрического центра грани.
							
							// G - грань по которой ставится граничное условие.
							// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
							// bound_id - номер граничного узла, начиная с нуля.
							// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
							if (octree1->b4N) {
								obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
									octree1->linkN3->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, false, false, p_centerG);
							}
							else {
								if (octree1->linkN == NULL) {
									obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
										NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, false, false, p_centerG);
								}
								else {
									obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
										octree1->linkN->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, false, false, p_centerG);
								}
							}

							// Зачем дробить если граница заканчивается NULL.
							maxbound++; // Вместо четырёх один.
						}						
						else {


							bool binc1, binc2, binc3, binc4;
							integer lsid1, lsid2, lsid3, lsid4;

							if (octree1->linkS0 == NULL) {
#if doubleintprecision == 1
								printf("error NULL sosedi S0 elm=%lld\n", maxelm);
#else
								printf("error NULL sosedi S0 elm=%d\n", maxelm);
#endif
								
								//getchar();
								system("PAUSE");

								/*

								G = SSIDE;
								// граничная грань:
								sosedi[G][maxelm].iNODE1 = maxelm_memo + maxbound; // граничные КО нумеруются в последнюю очередь,
								sosedi[G][maxelm].bdroblenie4 = true;
								// после maxelm_memo внутренних КО.
								// Вычисление дальнего соседа
								CALC_GG(G,GG); // Вычисляет грань GG по грани G.

								sosedi[GG][maxelm].iNODE1 = NON_EXISTENT_NODE; // соседа нет,
								sosedi[GG][maxelm].bdroblenie4 = true;

								#if doubleintprecision == 1
										printf("error NULL sosedi S0 elm=%lld\n", maxelm);
								#else
										printf("error NULL sosedi S0 elm=%d\n", maxelm);
								#endif
								
								getchar();

								// Эта ситуация возможна и здесь надо предусмотреть её обработку.
								DEFINE_BOUNDARY_PROPERTIES_FLOW(prop, prop_b, maxelm, maxbound); // В граничный узел сносятся свойства прилегающего КО.
								
								doublereal dS_loc = 0.0;

								// G - грань по которой ставится граничное условие.
								// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
								// bound_id - номер граничного узла, начиная с нуля.
								// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
								if (octree1->b4N) {
								obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
								octree1->linkN3->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
								}
								else {
								if (octree1->linkN == NULL) {
								obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
								NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
								}
								else {
								obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
								octree1->linkN->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
								}
								}


								maxbound++;
								*/
							}
							else {
								if (octree1->linkS0->inum_FD == 0) {

									G = SSIDE;
									// граничная грань:
									sosedi[G][maxelm].iNODE1 = maxelm_memo + maxbound; // граничные КО нумеруются в последнюю очередь,
									sosedi[G][maxelm].bdroblenie4 = true;
									// после maxelm_memo внутренних КО.
									// Вычисление дальнего соседа
									CALC_GG(G, GG); // Вычисляет грань GG по грани G.

									sosedi[GG][maxelm].iNODE1 = NON_EXISTENT_NODE; // соседа нет,
									sosedi[GG][maxelm].bdroblenie4 = true;

#if doubleintprecision == 1
									//printf("error inum_FD == 0 sosedi S0 elm=%lld\n", maxelm);
#else
									//printf("error inum_FD == 0 sosedi S0 elm=%d\n", maxelm);
#endif
									
									//getchar();

									// Эта ситуация возможна и здесь надо предусмотреть её обработку.
									DEFINE_BOUNDARY_PROPERTIES_FLOW(prop, prop_b, maxelm, maxbound); // В граничный узел сносятся свойства прилегающего КО.
									
									// Sx*Sz
									doublereal dS_loc = fabs(octree1->linkS0->p7.z - octree1->linkS0->p3.z)*fabs(octree1->linkS0->p2.x - octree1->linkS0->p3.x);
									TOCHKA p_centerG;
									p_centerG.z = 0.5*(octree1->linkS0->p7.z + octree1->linkS0->p3.z);
									p_centerG.y = octree1->linkS0->p2.y;
									p_centerG.x = 0.5*(octree1->linkS0->p2.x + octree1->linkS0->p3.x);


									// G - грань по которой ставится граничное условие.
									// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
									// bound_id - номер граничного узла, начиная с нуля.
									// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
									if (octree1->b4N) {
										obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
											octree1->linkN3->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
									}
									else {
										if (octree1->linkN == NULL) {
											obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
												NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
										}
										else {
											obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
												octree1->linkN->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
										}
									}

									maxbound++;
								}
								else {
									// сосед существует.
									if (bvisit[octree1->linkS0->inum_FD - 1]) {
										// узел уже был посещён
										// Здесь неким образом модифицировалось gran_t
										// Только внутренняя грань.
										G = SSIDE;
										// граничная грань:
										sosedi[G][maxelm].iNODE1 = octree1->linkS0->inum_FD - 1; // граничные КО нумеруются в последнюю очередь,
										sosedi[G][maxelm].bdroblenie4 = true;
										// после maxelm_memo внутренних КО.
										// Вычисление дальнего соседа
										CALC_GG(G, GG); // Вычисляет грань GG по грани G.

										// Определение дальнего соседа :
										// TODO 13 сентября 2016.
									}
									else {
										// узнать координаты центра грани и ориентацию в пространстве
										x_c = 0.5*(octree1->linkS0->p3.x + octree1->linkS0->p2.x);
										y_c = octree1->p0.y;
										z_c = 0.5*(octree1->linkS0->p3.z + octree1->linkS0->p7.z);
										iplane = XZ;
										patch_maxbound(iplane, s, ls, x_c, y_c, z_c, maxbound, binc1, lsid1);
										if (binc1) {
#if doubleintprecision == 1
											printf("error internal source S0 elm=%lld\n", maxelm);
#else
											printf("error internal source S0 elm=%d\n", maxelm);
#endif
											
											//getchar();
											system("PAUSE");
										}
										else {
											// Только внутренняя грань.
											G = SSIDE;
											// граничная грань:
											sosedi[G][maxelm].iNODE1 = octree1->linkS0->inum_FD - 1; // граничные КО нумеруются в последнюю очередь,
											sosedi[G][maxelm].bdroblenie4 = true;
											// после maxelm_memo внутренних КО.
											// Вычисление дальнего соседа
											CALC_GG(G, GG); // Вычисляет грань GG по грани G.

											// Определение дальнего соседа :
											// TODO 13 сентября 2016.
										}
									}
								}
							}

							if (octree1->linkS1 == NULL) {

								// 21 сентября 2016 (вырождение ячейки).

								// Заглушка.
								G = SSIDE;
								// граничная грань:
								sosedi[G][maxelm].iNODE2 = NON_EXISTENT_NODE; // граничные КО нумеруются в последнюю очередь,
								sosedi[G][maxelm].bdroblenie4 = true;
								// после maxelm_memo внутренних КО.
								// Вычисление дальнего соседа
								CALC_GG(G, GG); // Вычисляет грань GG по грани G.

								sosedi[GG][maxelm].iNODE2 = NON_EXISTENT_NODE; // сосед отсутствует,
								sosedi[GG][maxelm].bdroblenie4 = true;


								/*
								G = SSIDE;
								// граничная грань:
								sosedi[G][maxelm].iNODE2 = maxelm_memo + maxbound; // граничные КО нумеруются в последнюю очередь,
								sosedi[G][maxelm].bdroblenie4 = true;
								// после maxelm_memo внутренних КО.
								// Вычисление дальнего соседа
								CALC_GG(G,GG); // Вычисляет грань GG по грани G.

								sosedi[GG][maxelm].iNODE2 = NON_EXISTENT_NODE; // сосед отсутствует,
								sosedi[GG][maxelm].bdroblenie4 = true;

								#if doubleintprecision == 1
										printf("error NULL sosedi S1 elm=%lld\n", maxelm);
								#else
										printf("error NULL sosedi S1 elm=%d\n", maxelm);
								#endif
								
								getchar();

								// Эта ситуация возможна и здесь надо предусмотреть её обработку.
								DEFINE_BOUNDARY_PROPERTIES_FLOW(prop, prop_b, maxelm, maxbound); // В граничный узел сносятся свойства прилегающего КО.
								
								doublereal dS_loc = 0.0;

								// G - грань по которой ставится граничное условие.
								// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
								// bound_id - номер граничного узла, начиная с нуля.
								// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
								if (octree1->b4N) {
								obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
								octree1->linkN3->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
								}
								else {
								if (octree1->linkN == NULL) {
								obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
								NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
								}
								else {
								obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
								octree1->linkN->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
								}
								}

								maxbound++;
								*/
							}
							else {
								if (octree1->linkS1->inum_FD == 0) {

									G = SSIDE;
									// граничная грань:
									sosedi[G][maxelm].iNODE2 = maxelm_memo + maxbound; // граничные КО нумеруются в последнюю очередь,
									sosedi[G][maxelm].bdroblenie4 = true;
									// после maxelm_memo внутренних КО.
									// Вычисление дальнего соседа
									CALC_GG(G, GG); // Вычисляет грань GG по грани G.

									sosedi[GG][maxelm].iNODE2 = NON_EXISTENT_NODE; // сосед отсутствует,
									sosedi[GG][maxelm].bdroblenie4 = true;

#if doubleintprecision == 1
									//printf("error inum_FD == 0 sosedi S1 elm=%lld\n", maxelm);
#else
									//printf("error inum_FD == 0 sosedi S1 elm=%d\n", maxelm);
#endif
									
									//getchar();

									// Эта ситуация возможна и здесь надо предусмотреть её обработку.
									DEFINE_BOUNDARY_PROPERTIES_FLOW(prop, prop_b, maxelm, maxbound); // В граничный узел сносятся свойства прилегающего КО.
									
									// Sx*Sz
									doublereal dS_loc = fabs(octree1->linkS1->p7.z - octree1->linkS1->p3.z)*fabs(octree1->linkS1->p2.x - octree1->linkS1->p3.x);
									TOCHKA p_centerG;
									p_centerG.z = 0.5*(octree1->linkS1->p7.z + octree1->linkS1->p3.z);
									p_centerG.y = octree1->linkS1->p2.y;
									p_centerG.x = 0.5*(octree1->linkS1->p2.x + octree1->linkS1->p3.x);


									// G - грань по которой ставится граничное условие.
									// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
									// bound_id - номер граничного узла, начиная с нуля.
									// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
									if (octree1->b4N) {
										obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
											octree1->linkN3->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
									}
									else {
										if (octree1->linkN == NULL) {
											obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
												NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
										}
										else {
											obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
												octree1->linkN->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
										}
									}


									maxbound++;
								}
								else {
									// сосед существует.
									if (bvisit[octree1->linkS1->inum_FD - 1]) {
										// узел уже был посещён
										// Здесь неким образом модифицировалось gran_t
										// Только внутренняя грань.
										G = SSIDE;
										// граничная грань:
										sosedi[G][maxelm].iNODE2 = octree1->linkS1->inum_FD - 1; // граничные КО нумеруются в последнюю очередь,
										sosedi[G][maxelm].bdroblenie4 = true;
										// после maxelm_memo внутренних КО.
										// Вычисление дальнего соседа
										CALC_GG(G, GG); // Вычисляет грань GG по грани G.

										// Определение дальнего соседа :
										// TODO 13 сентября 2016.
									}
									else {
										// узнать координаты центра грани и ориентацию в пространстве
										x_c = 0.5*(octree1->linkS1->p3.x + octree1->linkS1->p2.x);
										y_c = octree1->p0.y;
										z_c = 0.5*(octree1->linkS1->p3.z + octree1->linkS1->p7.z);
										iplane = XZ;
										patch_maxbound(iplane, s, ls, x_c, y_c, z_c, maxbound, binc2, lsid2);
										if (binc2) {
#if doubleintprecision == 1
											printf("error internal source S1 elm=%lld\n", maxelm);
#else
											printf("error internal source S1 elm=%d\n", maxelm);
#endif
											
											//getchar();
											system("PAUSE");
										}
										else {
											// Только внутренняя грань.
											G = SSIDE;
											// граничная грань:
											sosedi[G][maxelm].iNODE2 = octree1->linkS1->inum_FD - 1; // граничные КО нумеруются в последнюю очередь,
											sosedi[G][maxelm].bdroblenie4 = true;
											// после maxelm_memo внутренних КО.
											// Вычисление дальнего соседа
											CALC_GG(G, GG); // Вычисляет грань GG по грани G.

											// Определение дальнего соседа :
											// TODO 13 сентября 2016.
										}
									}
								}
							}

							if (octree1->linkS4 == NULL) {

								// 21 сентября 2016 (вырождение ячейки).

								// Заглушка.
								G = SSIDE;
								// граничная грань:
								sosedi[G][maxelm].iNODE3 = NON_EXISTENT_NODE; // граничные КО нумеруются в последнюю очередь,
								sosedi[G][maxelm].bdroblenie4 = true;
								// после maxelm_memo внутренних КО.
								// Вычисление дальнего соседа
								CALC_GG(G, GG); // Вычисляет грань GG по грани G.
								sosedi[GG][maxelm].iNODE3 = NON_EXISTENT_NODE; // сосед отсутствует,
								sosedi[GG][maxelm].bdroblenie4 = true;


								/*
								G = SSIDE;
								// граничная грань:
								sosedi[G][maxelm].iNODE3 = maxelm_memo + maxbound; // граничные КО нумеруются в последнюю очередь,
								sosedi[G][maxelm].bdroblenie4 = true;
								// после maxelm_memo внутренних КО.
								// Вычисление дальнего соседа
								CALC_GG(G,GG); // Вычисляет грань GG по грани G.
								sosedi[GG][maxelm].iNODE3 = NON_EXISTENT_NODE; // сосед отсутствует,
								sosedi[GG][maxelm].bdroblenie4 = true;

								#if doubleintprecision == 1
										printf("error NULL sosedi S4 elm=%lld\n", maxelm);
								#else
										printf("error NULL sosedi S4 elm=%d\n", maxelm);
								#endif
								
								getchar();

								// Эта ситуация возможна и здесь надо предусмотреть её обработку.
								DEFINE_BOUNDARY_PROPERTIES_FLOW(prop, prop_b, maxelm, maxbound); // В граничный узел сносятся свойства прилегающего КО.
								
								doublereal dS_loc = 0.0;

								// G - грань по которой ставится граничное условие.
								// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
								// bound_id - номер граничного узла, начиная с нуля.
								// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
								if (octree1->b4N) {
								obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
								octree1->linkN3->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
								}
								else {
								if (octree1->linkN == NULL) {
								obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
								NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
								}
								else {
								obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
								octree1->linkN->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
								}
								}

								maxbound++;
								*/
							}
							else {
								if (octree1->linkS4->inum_FD == 0) {

									G = SSIDE;
									// граничная грань:
									sosedi[G][maxelm].iNODE3 = maxelm_memo + maxbound; // граничные КО нумеруются в последнюю очередь,
									sosedi[G][maxelm].bdroblenie4 = true;
									// после maxelm_memo внутренних КО.
									// Вычисление дальнего соседа
									CALC_GG(G, GG); // Вычисляет грань GG по грани G.
									sosedi[GG][maxelm].iNODE3 = NON_EXISTENT_NODE; // сосед отсутствует,
									sosedi[GG][maxelm].bdroblenie4 = true;

#if doubleintprecision == 1
									//printf("error inum_FD == 0 sosedi S4 elm=%lld\n", maxelm);
#else
									//printf("error inum_FD == 0 sosedi S4 elm=%d\n", maxelm);
#endif
									
									//getchar();

									// Эта ситуация возможна и здесь надо предусмотреть её обработку.
									DEFINE_BOUNDARY_PROPERTIES_FLOW(prop, prop_b, maxelm, maxbound); // В граничный узел сносятся свойства прилегающего КО.
									
									// Sx*Sz
									doublereal dS_loc = fabs(octree1->linkS4->p7.z - octree1->linkS4->p3.z)*fabs(octree1->linkS4->p2.x - octree1->linkS4->p3.x);
									TOCHKA p_centerG;
									p_centerG.z = 0.5*(octree1->linkS4->p7.z + octree1->linkS4->p3.z);
									p_centerG.y = octree1->linkS4->p2.y;
									p_centerG.x = 0.5*(octree1->linkS4->p2.x + octree1->linkS4->p3.x);


									// G - грань по которой ставится граничное условие.
									// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
									// bound_id - номер граничного узла, начиная с нуля.
									// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
									if (octree1->b4N) {
										obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
											octree1->linkN3->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
									}
									else {
										if (octree1->linkN == NULL) {
											obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
												NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
										}
										else {
											obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
												octree1->linkN->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
										}
									}

									maxbound++;
								}
								else {
									// сосед существует.
									if (bvisit[octree1->linkS4->inum_FD - 1]) {
										// узел уже был посещён
										// Здесь неким образом модифицировалось gran_t
										// Только внутренняя грань.
										G = SSIDE;
										// граничная грань:
										sosedi[G][maxelm].iNODE3 = octree1->linkS4->inum_FD - 1; // граничные КО нумеруются в последнюю очередь,
										sosedi[G][maxelm].bdroblenie4 = true;
										// после maxelm_memo внутренних КО.
										// Вычисление дальнего соседа
										CALC_GG(G, GG); // Вычисляет грань GG по грани G.

										// Определение дальнего соседа :
										// TODO 13 сентября 2016.
									}
									else {
										// узнать координаты центра грани и ориентацию в пространстве
										x_c = 0.5*(octree1->linkS4->p3.x + octree1->linkS4->p2.x);
										y_c = octree1->p0.y;
										z_c = 0.5*(octree1->linkS4->p3.z + octree1->linkS4->p7.z);
										iplane = XZ;
										patch_maxbound(iplane, s, ls, x_c, y_c, z_c, maxbound, binc3, lsid3);
										if (binc3) {
#if doubleintprecision == 1
											printf("error internal source S4 elm=%lld\n", maxelm);
#else
											printf("error internal source S4 elm=%d\n", maxelm);
#endif
											
											//getchar();
											system("PAUSE");
										}
										else {
											// Только внутренняя грань.
											G = SSIDE;
											// граничная грань:
											sosedi[G][maxelm].iNODE3 = octree1->linkS4->inum_FD - 1; // граничные КО нумеруются в последнюю очередь,
											sosedi[G][maxelm].bdroblenie4 = true;
											// после maxelm_memo внутренних КО.
											// Вычисление дальнего соседа
											CALC_GG(G, GG); // Вычисляет грань GG по грани G.

											// Определение дальнего соседа :
											// TODO 13 сентября 2016.
										}
									}
								}
							}

							if (octree1->linkS5 == NULL) {

								// 21 сентября 2016 (вырождение ячейки).

								// Заглушка.
								G = SSIDE;
								// граничная грань:
								sosedi[G][maxelm].iNODE4 = NON_EXISTENT_NODE; // граничные КО нумеруются в последнюю очередь,
								sosedi[G][maxelm].bdroblenie4 = true;
								// после maxelm_memo внутренних КО.
								// Вычисление дальнего соседа
								CALC_GG(G, GG); // Вычисляет грань GG по грани G.

								sosedi[GG][maxelm].iNODE4 = NON_EXISTENT_NODE; // сосед отсутствует,
								sosedi[GG][maxelm].bdroblenie4 = true;

								/*
								G = SSIDE;
								// граничная грань:
								sosedi[G][maxelm].iNODE4 = maxelm_memo + maxbound; // граничные КО нумеруются в последнюю очередь,
								sosedi[G][maxelm].bdroblenie4 = true;
								// после maxelm_memo внутренних КО.
								// Вычисление дальнего соседа
								CALC_GG(G,GG); // Вычисляет грань GG по грани G.

								sosedi[GG][maxelm].iNODE4 = NON_EXISTENT_NODE; // сосед отсутствует,
								sosedi[GG][maxelm].bdroblenie4 = true;

								#if doubleintprecision == 1
										printf("error NULL sosedi S5 elm=%lld\n", maxelm);
								#else
										printf("error NULL sosedi S5 elm=%d\n", maxelm);
								#endif
								
								getchar();

								// Эта ситуация возможна и здесь надо предусмотреть её обработку.
								DEFINE_BOUNDARY_PROPERTIES_FLOW(prop, prop_b, maxelm, maxbound); // В граничный узел сносятся свойства прилегающего КО.
								
								doublereal dS_loc = 0.0;

								// G - грань по которой ставится граничное условие.
								// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
								// bound_id - номер граничного узла, начиная с нуля.
								// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
								if (octree1->b4N) {
								obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
								octree1->linkN3->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
								}
								else {
								if (octree1->linkN == NULL) {
								obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
								NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
								}
								else {
								obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
								octree1->linkN->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
								}
								}

								maxbound++;
								*/
							}
							else {
								if (octree1->linkS5->inum_FD == 0) {

									G = SSIDE;
									// граничная грань:
									sosedi[G][maxelm].iNODE4 = maxelm_memo + maxbound; // граничные КО нумеруются в последнюю очередь,
									sosedi[G][maxelm].bdroblenie4 = true;
									// после maxelm_memo внутренних КО.
									// Вычисление дальнего соседа
									CALC_GG(G, GG); // Вычисляет грань GG по грани G.

									sosedi[GG][maxelm].iNODE4 = NON_EXISTENT_NODE; // сосед отсутствует,
									sosedi[GG][maxelm].bdroblenie4 = true;

#if doubleintprecision == 1
									//printf("error inum_FD == 0 sosedi S5 elm=%lld\n", maxelm);
#else
									//printf("error inum_FD == 0 sosedi S5 elm=%d\n", maxelm);
#endif
									
									//getchar();

									// Эта ситуация возможна и здесь надо предусмотреть её обработку.
									DEFINE_BOUNDARY_PROPERTIES_FLOW(prop, prop_b, maxelm, maxbound); // В граничный узел сносятся свойства прилегающего КО.
									
									// Sx*Sz
									doublereal dS_loc = fabs(octree1->linkS5->p7.z - octree1->linkS5->p3.z)*fabs(octree1->linkS5->p2.x - octree1->linkS5->p3.x);
									TOCHKA p_centerG;
									p_centerG.z = 0.5*(octree1->linkS5->p7.z + octree1->linkS5->p3.z);
									p_centerG.y = octree1->linkS5->p2.y;
									p_centerG.x = 0.5*(octree1->linkS5->p2.x + octree1->linkS5->p3.x);

									// G - грань по которой ставится граничное условие.
									// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
									// bound_id - номер граничного узла, начиная с нуля.
									// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
									if (octree1->b4N) {
										obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
											octree1->linkN3->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
									}
									else {
										if (octree1->linkN == NULL) {
											obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
												NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
										}
										else {
											obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
												octree1->linkN->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
										}
									}

									maxbound++;
								}
								else {
									// сосед существует.
									if (bvisit[octree1->linkS5->inum_FD - 1]) {
										// узел уже был посещён
										// Здесь неким образом модифицировалось gran_t
										// Только внутренняя грань.
										G = SSIDE;
										// граничная грань:
										sosedi[G][maxelm].iNODE4 = octree1->linkS5->inum_FD - 1; // граничные КО нумеруются в последнюю очередь,
										sosedi[G][maxelm].bdroblenie4 = true;
										// после maxelm_memo внутренних КО.
										// Вычисление дальнего соседа
										CALC_GG(G, GG); // Вычисляет грань GG по грани G.

										// Определение дальнего соседа :
										// TODO 13 сентября 2016.
									}
									else {
										// узнать координаты центра грани и ориентацию в пространстве
										x_c = 0.5*(octree1->linkS5->p3.x + octree1->linkS5->p2.x);
										y_c = octree1->p0.y;
										z_c = 0.5*(octree1->linkS5->p3.z + octree1->linkS5->p7.z);
										iplane = XZ;
										patch_maxbound(iplane, s, ls, x_c, y_c, z_c, maxbound, binc4, lsid4);
										if (binc4) {
#if doubleintprecision == 1
											printf("error internal source S5 elm=%lld\n", maxelm);
#else
											printf("error internal source S5 elm=%d\n", maxelm);
#endif
											
											//getchar();
											system("PAUSE");
										}
										else {
											// Только внутренняя грань.
											G = SSIDE;
											// граничная грань:
											sosedi[G][maxelm].iNODE4 = octree1->linkS5->inum_FD - 1; // граничные КО нумеруются в последнюю очередь,
											sosedi[G][maxelm].bdroblenie4 = true;
											// после maxelm_memo внутренних КО.
											// Вычисление дальнего соседа
											CALC_GG(G, GG); // Вычисляет грань GG по грани G.

											// Определение дальнего соседа :
											// TODO 13 сентября 2016.
										}
									}
								}
							}
						}
					}

					// сосед не существует
					// записываем идентификатор грани и
					// увеличиваем счётчик граней.
					if (!octree1->b4T) {
						if ((octree1->linkT == NULL)||((octree1->linkT != NULL)&&(octree1->linkT->inum_FD == 0))) {
							G = TSIDE;
							// граничная грань:
							INIT_SOSEDI(sosedi, maxelm, TSIDE, GG, false, maxelm_memo + maxbound, NON_EXISTENT_NODE, NON_EXISTENT_NODE, NON_EXISTENT_NODE, false, NON_EXISTENT_NODE, NON_EXISTENT_NODE, NON_EXISTENT_NODE, NON_EXISTENT_NODE); // заполнение sosedi[G][maxelm].
							DEFINE_BOUNDARY_PROPERTIES_FLOW(prop, prop_b, maxelm, maxbound); // В граничный узел сносятся свойства прилегающего КО.
							
							doublereal dSloc = CALC_SQUARE(G, maxelm, nvtx, pa); // Вычисление площади грани ячейки.
							TOCHKA p_centerG = CALC_CENTERG(G, maxelm, nvtx, pa); // Вычисление геометрического центра грани.
							
							// G - грань по которой ставится граничное условие.
							// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
							// bound_id - номер граничного узла, начиная с нуля.
							// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
							if (octree1->b4B) {
								obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
									octree1->linkB0->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, false, false, p_centerG);
							}
							else {
								if (octree1->linkB == NULL) {
									obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
										NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, false, false, p_centerG);
								}
								else {
									obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
										octree1->linkB->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, false, false, p_centerG);
								}

							}

							maxbound++;
						}
						else {
							if (octree1->linkT != NULL) {
								// сосед существует.
								if (bvisit[octree1->linkT->inum_FD - 1]) {
									
									integer return_current_node_number = sosedi[BSIDE][octree1->linkT->inum_FD - 1].iNODE1;
									if (octree1->linkT->b4B) {
										integer current_node_number = octree1->inum_FD - 1;
										bool bcontinue = true;
										if (bcontinue) {
											if (octree1->linkT->linkB0 != NULL) {
												if (octree1->linkT->linkB0->inum_FD - 1 == current_node_number) {
													bcontinue = false; // iNODE1
													if (sosedi[BSIDE][octree1->linkT->inum_FD - 1].iNODE1 == current_node_number) {
														return_current_node_number = sosedi[BSIDE][octree1->linkT->inum_FD - 1].iNODE1;
													}
													else {
														printf("this can not be T B0\n");
														//getchar();
														system("PAUSE");
													}
												}
											}
										}
										if (bcontinue) {
											if (octree1->linkT->linkB1 != NULL) {
												if (octree1->linkT->linkB1->inum_FD - 1 == current_node_number) {
													bcontinue = false; // iNODE2
													if (sosedi[BSIDE][octree1->linkT->inum_FD - 1].iNODE2 == current_node_number) {
														return_current_node_number = sosedi[BSIDE][octree1->linkT->inum_FD - 1].iNODE2;
													}
													else {
														printf("this can not be T B1\n");
														//getchar();
														system("PAUSE");
													}
												}
											}
										}
										if (bcontinue) {
											if (octree1->linkT->linkB2 != NULL) {
												if (octree1->linkT->linkB2->inum_FD - 1 == current_node_number) {
													bcontinue = false; // iNODE3
													if (sosedi[BSIDE][octree1->linkT->inum_FD - 1].iNODE3 == current_node_number) {
														return_current_node_number = sosedi[BSIDE][octree1->linkT->inum_FD - 1].iNODE3;
													}
													else {
														printf("this can not be T B2\n");
														//getchar();
														system("PAUSE");
													}
												}
											}
										}
										if (bcontinue) {
											if (octree1->linkT->linkB3 != NULL) {
												if (octree1->linkT->linkB3->inum_FD - 1 == current_node_number) {
													bcontinue = false; // iNODE4
													if (sosedi[BSIDE][octree1->linkT->inum_FD - 1].iNODE4 == current_node_number) {
														return_current_node_number = sosedi[BSIDE][octree1->linkT->inum_FD - 1].iNODE4;
													}
													else {
														printf("this can not be T B3\n");
														//getchar();
														system("PAUSE");
													}
												}
											}
										}

									}


									// узел уже был посещён
									// Здесь неким образом модифицировалось gran_t
									if (return_current_node_number >= maxelm_memo) {

#if doubleintprecision == 1
										printf("maxelm=%lld sosedi[BSIDE][octree1->linkT->inum_FD - 1].iNODE1=%lld\n", maxelm_memo, sosedi[BSIDE][octree1->linkT->inum_FD - 1].iNODE1);
										printf("inum_FD-1=%lld\n", octree1->linkT->inum_FD - 1);

										octTree* oc2 = octree1->linkT;
										if (oc2->b4B) {
											printf("b4B situation.\n");
											if (oc2->linkB0 == NULL) {
												printf("B0 is null\n");
											}
											else {
												printf("inum_FD-1=%lld\n", oc2->linkB0->inum_FD - 1);
											}
											if (oc2->linkB1 == NULL) {
												printf("B1 is null\n");
											}
											else {
												printf("inum_FD-1=%lld\n", oc2->linkB1->inum_FD - 1);
											}
											if (oc2->linkB2 == NULL) {
												printf("B2 is null\n");
											}
											else {
												printf("inum_FD-1=%lld\n", oc2->linkB2->inum_FD - 1);
											}
											if (oc2->linkB3 == NULL) {
												printf("B3 is null\n");
											}
											else {
												printf("inum_FD-1=%lld\n", oc2->linkB3->inum_FD - 1);
											}
											printf("octree1->inum_FD-1=%lld\n", octree1->inum_FD - 1);

										}
										else {
											if (oc2->linkB == NULL) {
												printf("linkB is null. \n");
											}
										}
#else
										printf("maxelm=%d sosedi[BSIDE][octree1->linkT->inum_FD - 1].iNODE1=%d\n", maxelm_memo, sosedi[BSIDE][octree1->linkT->inum_FD - 1].iNODE1);
										printf("inum_FD-1=%d\n", octree1->linkT->inum_FD - 1);

										octTree* oc2 = octree1->linkT;
										if (oc2->b4B) {
											printf("b4B situation.\n");
											if (oc2->linkB0 == NULL) {
												printf("B0 is null\n");
											}
											else {
												printf("inum_FD-1=%d\n", oc2->linkB0->inum_FD - 1);
											}
											if (oc2->linkB1 == NULL) {
												printf("B1 is null\n");
											}
											else {
												printf("inum_FD-1=%d\n", oc2->linkB1->inum_FD - 1);
											}
											if (oc2->linkB2 == NULL) {
												printf("B2 is null\n");
											}
											else {
												printf("inum_FD-1=%d\n", oc2->linkB2->inum_FD - 1);
											}
											if (oc2->linkB3 == NULL) {
												printf("B3 is null\n");
											}
											else {
												printf("inum_FD-1=%d\n", oc2->linkB3->inum_FD - 1);
											}
											printf("octree1->inum_FD-1=%d\n", octree1->inum_FD - 1);

										}
										else {
											if (oc2->linkB == NULL) {
												printf("linkB is null. \n");
											}
										}
#endif
										

										printf("bvisit error T ...\n");
										printf("model is incorrect. FATALL error!!!\n");
										printf("internal flat heat source can not be located\n");
										printf("inside the liquid.\n");
										printf("Please, press any key to exit...\n");
										//getchar();
										system("PAUSE");
										//exit(1);

										/*

										// Это граничная грань внутреннего источника тепла.
										// Внутренний источник тепла.
										G = TSIDE;
										// граничная грань:
										INIT_SOSEDI(sosedi, maxelm, TSIDE, GG, false, sosedi[BSIDE][octree1->linkT->inum_FD - 1].iNODE1, NON_EXISTENT_NODE, NON_EXISTENT_NODE, NON_EXISTENT_NODE, false, NON_EXISTENT_NODE, NON_EXISTENT_NODE, NON_EXISTENT_NODE, NON_EXISTENT_NODE); // заполнение sosedi[G][maxelm].
										
										
											// FLUID
											prop_b[RHO][sosedi[G][maxelm].iNODE1 - maxelm_memo] = prop[RHO][octree1->linkT->inum_FD - 1];
											prop_b[CP][sosedi[G][maxelm].iNODE1 - maxelm_memo] = prop[CP][octree1->linkT->inum_FD - 1];
											prop_b[LAM][sosedi[G][maxelm].iNODE1 - maxelm_memo] = prop[LAM][octree1->linkT->inum_FD - 1];
											prop_b[MULT_LAM_X][sosedi[G][maxelm].iNODE1 - maxelm_memo] = prop[MULT_LAM_X][octree1->linkT->inum_FD - 1];
											prop_b[MULT_LAM_Y][sosedi[G][maxelm].iNODE1 - maxelm_memo] = prop[MULT_LAM_Y][octree1->linkT->inum_FD - 1];
											prop_b[MULT_LAM_Z][sosedi[G][maxelm].iNODE1 - maxelm_memo] = prop[MULT_LAM_Z][octree1->linkT->inum_FD - 1];
									

										doublereal dSloc = CALC_SQUARE(G, maxelm, nvtx, pa); // Вычисление площади грани ячейки.
										TOCHKA p_centerG = CALC_CENTERG(G, maxelm, nvtx, pa); // Вычисление геометрического центра грани.
										
										// G - грань по которой ставится граничное условие.
										// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
										// bound_id - номер граничного узла, начиная с нуля.
										// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
										if (octree1->b4B) {
											obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, sosedi[BSIDE][octree1->linkT->inum_FD - 1].iNODE1 - maxelm_memo,
												octree1->linkB0->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, true, true, p_centerG);
										}
										else {
											if (octree1->linkB == NULL) {
												obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, sosedi[BSIDE][octree1->linkT->inum_FD - 1].iNODE1 - maxelm_memo,
													NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, true, true, p_centerG);
											}
											else {
												obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, sosedi[BSIDE][octree1->linkT->inum_FD - 1].iNODE1 - maxelm_memo,
													octree1->linkB->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, true, true, p_centerG);
											}
										}

										*/

									}
									else {
										// Это строго внутренний узел.
										// Внутренняя грань.
										G = TSIDE;
										// Внутренняя грань:
										sosedi[G][maxelm].iNODE1 = octree1->linkT->inum_FD - 1; // граничные КО нумеруются в последнюю очередь,
										sosedi[G][maxelm].bdroblenie4 = false;
										sosedi[G][maxelm].iNODE2 = NON_EXISTENT_NODE;
										sosedi[G][maxelm].iNODE3 = NON_EXISTENT_NODE;
										sosedi[G][maxelm].iNODE4 = NON_EXISTENT_NODE;
										// после maxelm_memo внутренних КО.
										// Вычисление дальнего соседа
										CALC_GG(G, GG); // Вычисляет грань GG по грани G.

										// Определение дальнего соседа :
										// TODO 13 сентября 2016.
									}

								}
								else {
									bool binc;
									integer lsid;
									// узнать координаты центра грани и ориентацию в пространстве
									x_c = 0.5*(octree1->p4.x + octree1->p5.x);
									y_c = 0.5*(octree1->p4.y + octree1->p7.y);
									z_c = octree1->p4.z;
									iplane = XY;
									patch_maxbound(iplane, s, ls, x_c, y_c, z_c, maxbound, binc, lsid);
									if (binc) {

										printf("model is incorrect. FATALL error!!!\n");
										printf("internal flat heat source can not be located\n");
										printf("inside the liquid.\n");
										printf("Please, press any key to exit...\n");
										//getchar();
										system("PAUSE");
										exit(1);

										/*

										// Внутренний источник тепла.
										G = TSIDE;
										// граничная грань:
										INIT_SOSEDI(sosedi, maxelm, TSIDE, GG, false, maxelm_memo + maxbound-1, NON_EXISTENT_NODE, NON_EXISTENT_NODE, NON_EXISTENT_NODE, false, NON_EXISTENT_NODE, NON_EXISTENT_NODE, NON_EXISTENT_NODE, NON_EXISTENT_NODE); // заполнение sosedi[G][maxelm].

										
											// FLUID
											prop_b[RHO][maxbound - 1] = prop[RHO][octree1->linkT->inum_FD - 1];
											prop_b[CP][maxbound - 1] = prop[CP][octree1->linkT->inum_FD - 1];
											prop_b[LAM][maxbound - 1] = prop[LAM][octree1->linkT->inum_FD - 1];
											prop_b[MULT_LAM_X][maxbound - 1] = prop[MULT_LAM_X][octree1->linkT->inum_FD - 1];
											prop_b[MULT_LAM_Y][maxbound - 1] = prop[MULT_LAM_Y][octree1->linkT->inum_FD - 1];
											prop_b[MULT_LAM_Z][maxbound - 1] = prop[MULT_LAM_Z][octree1->linkT->inum_FD - 1];
										
										doublereal dSloc = CALC_SQUARE(G, maxelm, nvtx, pa); // Вычисление площади грани ячейки.
										TOCHKA p_centerG = CALC_CENTERG(G, maxelm, nvtx, pa); // Вычисление геометрического центра грани.
										
										// G - грань по которой ставится граничное условие.
										// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
										// bound_id - номер граничного узла, начиная с нуля.
										// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
										if (octree1->b4B) {
											obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound - 1,
												octree1->linkB0->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, true, false, p_centerG);
										}
										else {
											if (octree1->linkB == NULL) {
												obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound - 1,
													NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, true, false, p_centerG);
											}
											else {
												obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound - 1,
													octree1->linkB->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, true, false, p_centerG);
											}
										}

										*/

									}
									else {
										// Внутренняя грань.
										G = TSIDE;
										// Внутренняя грань:
										sosedi[G][maxelm].iNODE1 = octree1->linkT->inum_FD - 1; // граничные КО нумеруются в последнюю очередь,
										sosedi[G][maxelm].bdroblenie4 = false;
										sosedi[G][maxelm].iNODE2 = NON_EXISTENT_NODE;
										sosedi[G][maxelm].iNODE3 = NON_EXISTENT_NODE;
										sosedi[G][maxelm].iNODE4 = NON_EXISTENT_NODE;
										// после maxelm_memo внутренних КО.
										// Вычисление дальнего соседа
										CALC_GG(G, GG); // Вычисляет грань GG по грани G.

										// Определение дальнего соседа :
										// TODO 13 сентября 2016.
									}
								}
							}
						}
					}
					else {
					// Присутствуют все 4 соседа: T: 4,5,6,7.

						if (((octree1->linkT4 == NULL) && (octree1->linkT5 == NULL) && (octree1->linkT6 == NULL) && (octree1->linkT7 == NULL))||
						(((octree1->linkT4 != NULL) && (octree1->linkT5 != NULL) && (octree1->linkT6 != NULL) && (octree1->linkT7 != NULL) &&
							(octree1->linkT4->inum_FD == 0) && (octree1->linkT5->inum_FD == 0) && (octree1->linkT6->inum_FD == 0) && (octree1->linkT7->inum_FD == 0)))||
							(((octree1->linkT4 != NULL) && (octree1->linkT5 == NULL) && (octree1->linkT6 == NULL) && (octree1->linkT7 != NULL) &&
							(octree1->linkT4->inum_FD == 0) && (octree1->linkT7->inum_FD == 0)))||
							(((octree1->linkT4 != NULL) && (octree1->linkT5 != NULL) && (octree1->linkT6 == NULL) && (octree1->linkT7 == NULL) &&
							(octree1->linkT4->inum_FD == 0) && (octree1->linkT5->inum_FD == 0)))||
							(((octree1->linkT4 != NULL) && (octree1->linkT5 == NULL) && (octree1->linkT6 == NULL) && (octree1->linkT7 == NULL) &&
							(octree1->linkT4->inum_FD == 0))))
						{
							// 1 - T4&&T5&&T6&&T7 - на границе расчётной области (РО).
							// 2 - T4&&T5&&T6&&T7 - на гранеице HOLLOW блока.
							// 3 - T4&&T7 на границе HOLLOW блока, T5&&T6 на границе РО.
							// 4 - T4&&T5 на границе HOLLOW блока, T6&T7 на границе РО.
							// 5 - T4 на границе HOLLOW блока, T5&&T6&&T7  на границе РО.

							G = TSIDE;
							// граничная грань:
							INIT_SOSEDI(sosedi, maxelm, TSIDE, GG, false, maxelm_memo + maxbound, NON_EXISTENT_NODE, NON_EXISTENT_NODE, NON_EXISTENT_NODE, false, NON_EXISTENT_NODE, NON_EXISTENT_NODE, NON_EXISTENT_NODE, NON_EXISTENT_NODE); // заполнение sosedi[G][maxelm].
							DEFINE_BOUNDARY_PROPERTIES_FLOW(prop, prop_b, maxelm, maxbound); // В граничный узел сносятся свойства прилегающего КО.

							doublereal dSloc = CALC_SQUARE(G, maxelm, nvtx, pa); // Вычисление площади грани ячейки.
							TOCHKA p_centerG = CALC_CENTERG(G, maxelm, nvtx, pa); // Вычисление геометрического центра грани.
							
							// G - грань по которой ставится граничное условие.
							// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
							// bound_id - номер граничного узла, начиная с нуля.
							// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
							if (octree1->b4B) {
								obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
									octree1->linkB0->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, false, false, p_centerG);
							}
							else {
								if (octree1->linkB == NULL) {
									obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
										NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, false, false, p_centerG);
								}
								else {
									obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
										octree1->linkB->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, false, false, p_centerG);
								}
							}

							// Зачем дробить если граница заканчивается NULL.
							maxbound++; // Вместо четырёх один.
						}						
						else {

							bool binc1, binc2, binc3, binc4;
							integer lsid1, lsid2, lsid3, lsid4;

							if (octree1->linkT4 == NULL) {

#if doubleintprecision == 1
								printf("error NULL sosedi T4 elm=%lld\n", maxelm);
#else
								printf("error NULL sosedi T4 elm=%d\n", maxelm);
#endif
								
								//getchar();
								system("PAUSE");

								/*

								G = TSIDE;
								// граничная грань:
								sosedi[G][maxelm].iNODE1 = maxelm_memo + maxbound; // граничные КО нумеруются в последнюю очередь,
								sosedi[G][maxelm].bdroblenie4 = true;
								// после maxelm_memo внутренних КО.
								// Вычисление дальнего соседа
								CALC_GG(G,GG); // Вычисляет грань GG по грани G.

								sosedi[GG][maxelm].iNODE1 = NON_EXISTENT_NODE; // сосед отсутствует,
								sosedi[GG][maxelm].bdroblenie4 = true;

								#if doubleintprecision == 1
										printf("error NULL sosedi T4 elm=%lld\n", maxelm);
								#else
										printf("error NULL sosedi T4 elm=%d\n", maxelm);
								#endif
								
								getchar();

								// Эта ситуация возможна и здесь надо предусмотреть её обработку.
								DEFINE_BOUNDARY_PROPERTIES_FLOW(prop, prop_b, maxelm, maxbound); // В граничный узел сносятся свойства прилегающего КО.
								
								doublereal dS_loc = 0.0;

								// G - грань по которой ставится граничное условие.
								// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
								// bound_id - номер граничного узла, начиная с нуля.
								// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
								if (octree1->b4B) {
								obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
								octree1->linkB0->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
								}
								else {
								if (octree1->linkB == NULL) {
								obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
								NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
								}
								else {
								obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
								octree1->linkB->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
								}

								}

								maxbound++;
								*/
							}
							else {
								if (octree1->linkT4->inum_FD == 0) {

									G = TSIDE;
									// граничная грань:
									sosedi[G][maxelm].iNODE1 = maxelm_memo + maxbound; // граничные КО нумеруются в последнюю очередь,
									sosedi[G][maxelm].bdroblenie4 = true;
									// после maxelm_memo внутренних КО.
									// Вычисление дальнего соседа
									CALC_GG(G, GG); // Вычисляет грань GG по грани G.

									sosedi[GG][maxelm].iNODE1 = NON_EXISTENT_NODE; // сосед отсутствует,
									sosedi[GG][maxelm].bdroblenie4 = true;

#if doubleintprecision == 1
									//printf("error inum_FD == 0 sosedi T4 elm=%lld\n", maxelm);
#else
									//printf("error inum_FD == 0 sosedi T4 elm=%d\n", maxelm);
#endif
									
									//getchar();

									// Эта ситуация возможна и здесь надо предусмотреть её обработку.
									DEFINE_BOUNDARY_PROPERTIES_FLOW(prop, prop_b, maxelm, maxbound); // В граничный узел сносятся свойства прилегающего КО.
									
									// Sx*Sy
									doublereal dS_loc = fabs(octree1->linkT4->p3.y - octree1->linkT4->p0.y)*fabs(octree1->linkT4->p1.x - octree1->linkT4->p0.x);
									TOCHKA p_centerG;
									p_centerG.z =  octree1->linkT4->p0.z;
									p_centerG.y = 0.5*(octree1->linkT4->p3.y + octree1->linkT4->p0.y);
									p_centerG.x = 0.5*(octree1->linkT4->p1.x + octree1->linkT4->p0.x);

									// G - грань по которой ставится граничное условие.
									// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
									// bound_id - номер граничного узла, начиная с нуля.
									// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
									if (octree1->b4B) {
										obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
											octree1->linkB0->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
									}
									else {
										if (octree1->linkB == NULL) {
											obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
												NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
										}
										else {
											obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
												octree1->linkB->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
										}

									}

									maxbound++;
								}
								else {
									// сосед существует.
									if (bvisit[octree1->linkT4->inum_FD - 1]) {
										// узел уже был посещён
										// Здесь неким образом модифицировалось gran_t
										// Только внутренняя грань.
										G = TSIDE;
										// граничная грань:
										sosedi[G][maxelm].iNODE1 = octree1->linkT4->inum_FD - 1; // граничные КО нумеруются в последнюю очередь,
										sosedi[G][maxelm].bdroblenie4 = true;
										// после maxelm_memo внутренних КО.
										// Вычисление дальнего соседа
										CALC_GG(G, GG); // Вычисляет грань GG по грани G.

										// Определение дальнего соседа :
										// TODO 13 сентября 2016.
									}
									else {
										// узнать координаты центра грани и ориентацию в пространстве
										x_c = 0.5*(octree1->linkT4->p0.x + octree1->linkT4->p1.x);
										y_c = 0.5*(octree1->linkT4->p0.y + octree1->linkT4->p3.y);
										z_c = octree1->p4.z;
										iplane = XY;
										patch_maxbound(iplane, s, ls, x_c, y_c, z_c, maxbound, binc1, lsid1);
										if (binc1) {
#if doubleintprecision == 1
											printf("error internal source T4 elm=%lld\n", maxelm);
#else
											printf("error internal source T4 elm=%d\n", maxelm);
#endif
											
											//getchar();
											system("PAUSE");
										}
										else {
											// Только внутренняя грань.
											G = TSIDE;
											// граничная грань:
											sosedi[G][maxelm].iNODE1 = octree1->linkT4->inum_FD - 1; // граничные КО нумеруются в последнюю очередь,
											sosedi[G][maxelm].bdroblenie4 = true;
											// после maxelm_memo внутренних КО.
											// Вычисление дальнего соседа
											CALC_GG(G, GG); // Вычисляет грань GG по грани G.

											// Определение дальнего соседа :
											// TODO 13 сентября 2016.
										}
									}
								}
							}

							if (octree1->linkT5 == NULL) {

								// 21 сентября 2016 (вырождение ячейки).

								// Заглушка.
								G = TSIDE;
								// граничная грань:
								sosedi[G][maxelm].iNODE2 = NON_EXISTENT_NODE; // граничные КО нумеруются в последнюю очередь,
								sosedi[G][maxelm].bdroblenie4 = true;
								// после maxelm_memo внутренних КО.
								// Вычисление дальнего соседа
								CALC_GG(G, GG); // Вычисляет грань GG по грани G.

								sosedi[GG][maxelm].iNODE2 = NON_EXISTENT_NODE; // соседа нет,
								sosedi[GG][maxelm].bdroblenie4 = true;

								/*
								G = TSIDE;
								// граничная грань:
								sosedi[G][maxelm].iNODE2 = maxelm_memo + maxbound; // граничные КО нумеруются в последнюю очередь,
								sosedi[G][maxelm].bdroblenie4 = true;
								// после maxelm_memo внутренних КО.
								// Вычисление дальнего соседа
								CALC_GG(G,GG); // Вычисляет грань GG по грани G.

								sosedi[GG][maxelm].iNODE2 = NON_EXISTENT_NODE; // соседа нет,
								sosedi[GG][maxelm].bdroblenie4 = true;

								#if doubleintprecision == 1
										printf("error NULL sosedi T5 elm=%lld\n", maxelm);
								#else
										printf("error NULL sosedi T5 elm=%d\n", maxelm);
								#endif
								
								getchar();

								// Эта ситуация возможна и здесь надо предусмотреть её обработку.
								DEFINE_BOUNDARY_PROPERTIES_FLOW(prop, prop_b, maxelm, maxbound); // В граничный узел сносятся свойства прилегающего КО.

								doublereal dS_loc = 0.0;

								// G - грань по которой ставится граничное условие.
								// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
								// bound_id - номер граничного узла, начиная с нуля.
								// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
								if (octree1->b4B) {
								obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
								octree1->linkB0->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
								}
								else {
								if (octree1->linkB == NULL) {
								obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
								NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
								}
								else {
								obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
								octree1->linkB->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
								}

								}

								maxbound++;
								*/
							}
							else {
								if (octree1->linkT5->inum_FD == 0) {

									G = TSIDE;
									// граничная грань:
									sosedi[G][maxelm].iNODE2 = maxelm_memo + maxbound; // граничные КО нумеруются в последнюю очередь,
									sosedi[G][maxelm].bdroblenie4 = true;
									// после maxelm_memo внутренних КО.
									// Вычисление дальнего соседа
									CALC_GG(G, GG); // Вычисляет грань GG по грани G.

									sosedi[GG][maxelm].iNODE2 = NON_EXISTENT_NODE; // соседа нет,
									sosedi[GG][maxelm].bdroblenie4 = true;

#if doubleintprecision == 1
									//printf("error inum_FD == 0 sosedi T5 elm=%lld\n", maxelm);
#else
									//printf("error inum_FD == 0 sosedi T5 elm=%d\n", maxelm);
#endif
									
									//getchar();

									// Эта ситуация возможна и здесь надо предусмотреть её обработку.
									DEFINE_BOUNDARY_PROPERTIES_FLOW(prop, prop_b, maxelm, maxbound); // В граничный узел сносятся свойства прилегающего КО.

									// Sx*Sy
									doublereal dS_loc = fabs(octree1->linkT5->p3.y - octree1->linkT5->p0.y)*fabs(octree1->linkT5->p1.x - octree1->linkT5->p0.x);
									TOCHKA p_centerG;
									p_centerG.z = octree1->linkT5->p0.z;
									p_centerG.y = 0.5*(octree1->linkT5->p3.y + octree1->linkT5->p0.y);
									p_centerG.x = 0.5*(octree1->linkT5->p1.x + octree1->linkT5->p0.x);


									// G - грань по которой ставится граничное условие.
									// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
									// bound_id - номер граничного узла, начиная с нуля.
									// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
									if (octree1->b4B) {
										obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
											octree1->linkB0->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
									}
									else {
										if (octree1->linkB == NULL) {
											obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
												NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
										}
										else {
											obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
												octree1->linkB->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
										}

									}

									maxbound++;
								}
								else {
									// сосед существует.
									if (bvisit[octree1->linkT5->inum_FD - 1]) {
										// узел уже был посещён
										// Здесь неким образом модифицировалось gran_t
										// Только внутренняя грань.
										G = TSIDE;
										// граничная грань:
										sosedi[G][maxelm].iNODE2 = octree1->linkT5->inum_FD - 1; // граничные КО нумеруются в последнюю очередь,
										sosedi[G][maxelm].bdroblenie4 = true;
										// после maxelm_memo внутренних КО.
										// Вычисление дальнего соседа
										CALC_GG(G, GG); // Вычисляет грань GG по грани G.

										// Определение дальнего соседа :
										// TODO 13 сентября 2016.
									}
									else {
										// узнать координаты центра грани и ориентацию в пространстве
										x_c = 0.5*(octree1->linkT5->p0.x + octree1->linkT5->p1.x);
										y_c = 0.5*(octree1->linkT5->p0.y + octree1->linkT5->p3.y);
										z_c = octree1->p4.z;
										iplane = XY;
										patch_maxbound(iplane, s, ls, x_c, y_c, z_c, maxbound, binc2, lsid2);
										if (binc2) {
#if doubleintprecision == 1
											printf("error internal source T5 elm=%lld\n", maxelm);
#else
											printf("error internal source T5 elm=%d\n", maxelm);
#endif
											
											//getchar();
											system("PAUSE");
										}
										else {
											// Только внутренняя грань.
											G = TSIDE;
											// граничная грань:
											sosedi[G][maxelm].iNODE2 = octree1->linkT5->inum_FD - 1; // граничные КО нумеруются в последнюю очередь,
											sosedi[G][maxelm].bdroblenie4 = true;
											// после maxelm_memo внутренних КО.
											// Вычисление дальнего соседа
											CALC_GG(G, GG); // Вычисляет грань GG по грани G.

											// Определение дальнего соседа :
											// TODO 13 сентября 2016.
										}
									}
								}
							}

							if (octree1->linkT6 == NULL) {

								// 21 сентября 2016 (вырождение ячейки).

								// Заглушка.
								G = TSIDE;
								// граничная грань:
								sosedi[G][maxelm].iNODE3 = NON_EXISTENT_NODE; // граничные КО нумеруются в последнюю очередь,
								sosedi[G][maxelm].bdroblenie4 = true;
								// после maxelm_memo внутренних КО.
								// Вычисление дальнего соседа
								CALC_GG(G, GG); // Вычисляет грань GG по грани G.

								sosedi[GG][maxelm].iNODE3 = NON_EXISTENT_NODE; // сосед отсутствует,
								sosedi[GG][maxelm].bdroblenie4 = true;

								/*
								G = TSIDE;
								// граничная грань:
								sosedi[G][maxelm].iNODE3 = maxelm_memo + maxbound; // граничные КО нумеруются в последнюю очередь,
								sosedi[G][maxelm].bdroblenie4 = true;
								// после maxelm_memo внутренних КО.
								// Вычисление дальнего соседа
								CALC_GG(G,GG); // Вычисляет грань GG по грани G.

								sosedi[GG][maxelm].iNODE3 = NON_EXISTENT_NODE; // сосед отсутствует,
								sosedi[GG][maxelm].bdroblenie4 = true;

								#if doubleintprecision == 1
										printf("error NULL sosedi T6 elm=%lld\n", maxelm);
								#else
										printf("error NULL sosedi T6 elm=%d\n", maxelm);
								#endif
								
								getchar();

								// Эта ситуация возможна и здесь надо предусмотреть её обработку.
								DEFINE_BOUNDARY_PROPERTIES_FLOW(prop, prop_b, maxelm, maxbound); // В граничный узел сносятся свойства прилегающего КО.

								doublereal dS_loc = 0.0;

								// G - грань по которой ставится граничное условие.
								// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
								// bound_id - номер граничного узла, начиная с нуля.
								// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
								if (octree1->b4B) {
								obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
								octree1->linkB0->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
								}
								else {
								if (octree1->linkB == NULL) {
								obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
								NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
								}
								else {
								obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
								octree1->linkB->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
								}

								}

								maxbound++;
								*/
							}
							else {
								if (octree1->linkT6->inum_FD == 0) {

									G = TSIDE;
									// граничная грань:
									sosedi[G][maxelm].iNODE3 = maxelm_memo + maxbound; // граничные КО нумеруются в последнюю очередь,
									sosedi[G][maxelm].bdroblenie4 = true;
									// после maxelm_memo внутренних КО.
									// Вычисление дальнего соседа
									CALC_GG(G, GG); // Вычисляет грань GG по грани G.

									sosedi[GG][maxelm].iNODE3 = NON_EXISTENT_NODE; // сосед отсутствует,
									sosedi[GG][maxelm].bdroblenie4 = true;

#if doubleintprecision == 1
									//printf("error inum_FD == 0 sosedi T6 elm=%lld\n", maxelm);
#else
									//printf("error inum_FD == 0 sosedi T6 elm=%d\n", maxelm);
#endif
									
									//getchar();

									// Эта ситуация возможна и здесь надо предусмотреть её обработку.
									DEFINE_BOUNDARY_PROPERTIES_FLOW(prop, prop_b, maxelm, maxbound); // В граничный узел сносятся свойства прилегающего КО.

									// Sx*Sy
									doublereal dS_loc = fabs(octree1->linkT6->p3.y - octree1->linkT6->p0.y)*fabs(octree1->linkT6->p1.x - octree1->linkT6->p0.x);
									TOCHKA p_centerG;
									p_centerG.z = octree1->linkT6->p0.z;
									p_centerG.y = 0.5*(octree1->linkT6->p3.y + octree1->linkT6->p0.y);
									p_centerG.x = 0.5*(octree1->linkT6->p1.x + octree1->linkT6->p0.x);

									// G - грань по которой ставится граничное условие.
									// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
									// bound_id - номер граничного узла, начиная с нуля.
									// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
									if (octree1->b4B) {
										obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
											octree1->linkB0->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
									}
									else {
										if (octree1->linkB == NULL) {
											obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
												NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
										}
										else {
											obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
												octree1->linkB->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
										}

									}


									maxbound++;
								}
								else {
									// сосед существует.
									if (bvisit[octree1->linkT6->inum_FD - 1]) {
										// узел уже был посещён
										// Здесь неким образом модифицировалось gran_t
										// Только внутренняя грань.
										G = TSIDE;
										// граничная грань:
										sosedi[G][maxelm].iNODE3 = octree1->linkT6->inum_FD - 1; // граничные КО нумеруются в последнюю очередь,
										sosedi[G][maxelm].bdroblenie4 = true;
										// после maxelm_memo внутренних КО.
										// Вычисление дальнего соседа
										CALC_GG(G, GG); // Вычисляет грань GG по грани G.

										// Определение дальнего соседа :
										// TODO 13 сентября 2016.
									}
									else {
										// узнать координаты центра грани и ориентацию в пространстве
										x_c = 0.5*(octree1->linkT6->p0.x + octree1->linkT6->p1.x);
										y_c = 0.5*(octree1->linkT6->p0.y + octree1->linkT6->p3.y);
										z_c = octree1->p4.z;
										iplane = XY;
										patch_maxbound(iplane, s, ls, x_c, y_c, z_c, maxbound, binc3, lsid3);
										if (binc3) {
#if doubleintprecision == 1
											printf("error internal source T6 elm=%lld\n", maxelm);
#else
											printf("error internal source T6 elm=%d\n", maxelm);
#endif
											
											//getchar();
											system("PAUSE");
										}
										else {
											// Только внутренняя грань.
											G = TSIDE;
											// граничная грань:
											sosedi[G][maxelm].iNODE3 = octree1->linkT6->inum_FD - 1; // граничные КО нумеруются в последнюю очередь,
											sosedi[G][maxelm].bdroblenie4 = true;
											// после maxelm_memo внутренних КО.
											// Вычисление дальнего соседа
											CALC_GG(G, GG); // Вычисляет грань GG по грани G.

											// Определение дальнего соседа :
											// TODO 13 сентября 2016.
										}
									}
								}
							}

							if (octree1->linkT7 == NULL) {

								// 21 сентября 2016 (вырождение ячейки).

								// Заглушка.
								G = TSIDE;
								// граничная грань:
								sosedi[G][maxelm].iNODE4 = NON_EXISTENT_NODE; // граничные КО нумеруются в последнюю очередь,
								sosedi[G][maxelm].bdroblenie4 = true;
								// после maxelm_memo внутренних КО.
								// Вычисление дальнего соседа
								CALC_GG(G, GG); // Вычисляет грань GG по грани G.

								sosedi[GG][maxelm].iNODE4 = NON_EXISTENT_NODE; // соседа нет,
								sosedi[GG][maxelm].bdroblenie4 = true;

								/*
								G = TSIDE;
								// граничная грань:
								sosedi[G][maxelm].iNODE4 = maxelm_memo + maxbound; // граничные КО нумеруются в последнюю очередь,
								sosedi[G][maxelm].bdroblenie4 = true;
								// после maxelm_memo внутренних КО.
								// Вычисление дальнего соседа
								CALC_GG(G,GG); // Вычисляет грань GG по грани G.

								sosedi[GG][maxelm].iNODE4 = NON_EXISTENT_NODE; // соседа нет,
								sosedi[GG][maxelm].bdroblenie4 = true;

								#if doubleintprecision == 1
										printf("error NULL sosedi T7 elm=%lld\n", maxelm);
								#else
										printf("error NULL sosedi T7 elm=%d\n", maxelm);
								#endif
								
								getchar();

								// Эта ситуация возможна и здесь надо предусмотреть её обработку.
								DEFINE_BOUNDARY_PROPERTIES_FLOW(prop, prop_b, maxelm, maxbound); // В граничный узел сносятся свойства прилегающего КО.

								doublereal dS_loc = 0.0;

								// G - грань по которой ставится граничное условие.
								// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
								// bound_id - номер граничного узла, начиная с нуля.
								// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
								if (octree1->b4B) {
								obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
								octree1->linkB0->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
								}
								else {
								if (octree1->linkB == NULL) {
								obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
								NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
								}
								else {
								obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
								octree1->linkB->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
								}

								}

								maxbound++;
								*/
							}
							else {
								if (octree1->linkT7->inum_FD == 0) {

									G = TSIDE;
									// граничная грань:
									sosedi[G][maxelm].iNODE4 = maxelm_memo + maxbound; // граничные КО нумеруются в последнюю очередь,
									sosedi[G][maxelm].bdroblenie4 = true;
									// после maxelm_memo внутренних КО.
									// Вычисление дальнего соседа
									CALC_GG(G, GG); // Вычисляет грань GG по грани G.

									sosedi[GG][maxelm].iNODE4 = NON_EXISTENT_NODE; // соседа нет,
									sosedi[GG][maxelm].bdroblenie4 = true;

#if doubleintprecision == 1
									//printf("error inum_FD == 0 sosedi T7 elm=%lld\n", maxelm);
#else
									//printf("error inum_FD == 0 sosedi T7 elm=%d\n", maxelm);
#endif
									
									//getchar();

									// Эта ситуация возможна и здесь надо предусмотреть её обработку.
									DEFINE_BOUNDARY_PROPERTIES_FLOW(prop, prop_b, maxelm, maxbound); // В граничный узел сносятся свойства прилегающего КО.

									// Sx*Sy
									doublereal dS_loc = fabs(octree1->linkT7->p3.y - octree1->linkT7->p0.y)*fabs(octree1->linkT7->p1.x - octree1->linkT7->p0.x);
									TOCHKA p_centerG;
									p_centerG.z = octree1->linkT7->p0.z;
									p_centerG.y = 0.5*(octree1->linkT7->p3.y + octree1->linkT7->p0.y);
									p_centerG.x = 0.5*(octree1->linkT7->p1.x + octree1->linkT7->p0.x);

									// G - грань по которой ставится граничное условие.
									// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
									// bound_id - номер граничного узла, начиная с нуля.
									// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
									if (octree1->b4B) {
										obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
											octree1->linkB0->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
									}
									else {
										if (octree1->linkB == NULL) {
											obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
												NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
										}
										else {
											obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
												octree1->linkB->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
										}

									}

									maxbound++;
								}
								else {
									// сосед существует.
									if (bvisit[octree1->linkT7->inum_FD - 1]) {
										// узел уже был посещён
										// Здесь неким образом модифицировалось gran_t
										// Только внутренняя грань.
										G = TSIDE;
										// граничная грань:
										sosedi[G][maxelm].iNODE4 = octree1->linkT7->inum_FD - 1; // граничные КО нумеруются в последнюю очередь,
										sosedi[G][maxelm].bdroblenie4 = true;
										// после maxelm_memo внутренних КО.
										// Вычисление дальнего соседа
										CALC_GG(G, GG); // Вычисляет грань GG по грани G.

										// Определение дальнего соседа :
										// TODO 13 сентября 2016.
									}
									else {
										// узнать координаты центра грани и ориентацию в пространстве
										x_c = 0.5*(octree1->linkT7->p0.x + octree1->linkT7->p1.x);
										y_c = 0.5*(octree1->linkT7->p0.y + octree1->linkT7->p3.y);
										z_c = octree1->p4.z;
										iplane = XY;
										patch_maxbound(iplane, s, ls, x_c, y_c, z_c, maxbound, binc4, lsid4);
										if (binc4) {
#if doubleintprecision == 1
											printf("error internal source T7 elm=%lld\n", maxelm);
#else
											printf("error internal source T7 elm=%d\n", maxelm);
#endif
											
											//getchar();
											system("PAUSE");
										}
										else {
											// Только внутренняя грань.
											G = TSIDE;
											// граничная грань:
											sosedi[G][maxelm].iNODE4 = octree1->linkT7->inum_FD - 1; // граничные КО нумеруются в последнюю очередь,
											sosedi[G][maxelm].bdroblenie4 = true;
											// после maxelm_memo внутренних КО.
											// Вычисление дальнего соседа
											CALC_GG(G, GG); // Вычисляет грань GG по грани G.

											// Определение дальнего соседа :
											// TODO 13 сентября 2016.
										}
									}
								}
							}
						}
					}

					// сосед не существует
					// записываем идентификатор грани и
					// увеличиваем счётчик граней.
					if (!octree1->b4B) {
						if ((octree1->linkB == NULL)||((octree1->linkB != NULL)&&(octree1->linkB->inum_FD == 0))) {
							G = BSIDE;
							// граничная грань:
							INIT_SOSEDI(sosedi, maxelm, BSIDE, GG, false, maxelm_memo + maxbound, NON_EXISTENT_NODE, NON_EXISTENT_NODE, NON_EXISTENT_NODE, false, NON_EXISTENT_NODE, NON_EXISTENT_NODE, NON_EXISTENT_NODE, NON_EXISTENT_NODE); // заполнение sosedi[G][maxelm].
							DEFINE_BOUNDARY_PROPERTIES_FLOW(prop, prop_b, maxelm, maxbound); // В граничный узел сносятся свойства прилегающего КО.

							doublereal dSloc = CALC_SQUARE(G, maxelm, nvtx, pa); // Вычисление площади грани ячейки.
							TOCHKA p_centerG = CALC_CENTERG(G, maxelm, nvtx, pa); // Вычисление геометрического центра грани.
							
							// G - грань по которой ставится граничное условие.
							// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
							// bound_id - номер граничного узла, начиная с нуля.
							// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
							if (octree1->b4T) {
								obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
									octree1->linkT4->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, false, false, p_centerG);
							}
							else {
								if (octree1->linkT == NULL) {
									obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
										NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, false, false, p_centerG);
								}
								else {
									obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
										octree1->linkT->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, false, false, p_centerG);
								}
							}

							maxbound++;
						}
						else {
							if (octree1->linkB != NULL) {
								// сосед существует.
								if (bvisit[octree1->linkB->inum_FD - 1]) {


									integer return_current_node_number = sosedi[TSIDE][octree1->linkB->inum_FD - 1].iNODE1;
									if (octree1->linkB->b4B) {
										integer current_node_number = octree1->inum_FD - 1;
										bool bcontinue = true;
										if (bcontinue) {
											if (octree1->linkB->linkT4 != NULL) {
												if (octree1->linkB->linkT4->inum_FD - 1 == current_node_number) {
													bcontinue = false; // iNODE1
													if (sosedi[TSIDE][octree1->linkB->inum_FD - 1].iNODE1 == current_node_number) {
														return_current_node_number = sosedi[TSIDE][octree1->linkB->inum_FD - 1].iNODE1;
													}
													else {
														printf("this can not be B T4\n");
														//getchar();
														system("PAUSE");
													}
												}
											}
										}
										if (bcontinue) {
											if (octree1->linkB->linkT5 != NULL) {
												if (octree1->linkB->linkT5->inum_FD - 1 == current_node_number) {
													bcontinue = false; // iNODE2
													if (sosedi[TSIDE][octree1->linkB->inum_FD - 1].iNODE2 == current_node_number) {
														return_current_node_number = sosedi[TSIDE][octree1->linkB->inum_FD - 1].iNODE2;
													}
													else {
														printf("this can not be B T5\n");
														//getchar();
														system("PAUSE");
													}
												}
											}
										}
										if (bcontinue) {
											if (octree1->linkB->linkT6 != NULL) {
												if (octree1->linkB->linkT6->inum_FD - 1 == current_node_number) {
													bcontinue = false; // iNODE3
													if (sosedi[TSIDE][octree1->linkB->inum_FD - 1].iNODE3 == current_node_number) {
														return_current_node_number = sosedi[TSIDE][octree1->linkB->inum_FD - 1].iNODE3;
													}
													else {
														printf("this can not be B T6\n");
														//getchar();
														system("PAUSE");
													}
												}
											}
										}
										if (bcontinue) {
											if (octree1->linkB->linkT7 != NULL) {
												if (octree1->linkB->linkT7->inum_FD - 1 == current_node_number) {
													bcontinue = false; // iNODE4
													if (sosedi[TSIDE][octree1->linkB->inum_FD - 1].iNODE4 == current_node_number) {
														return_current_node_number = sosedi[TSIDE][octree1->linkB->inum_FD - 1].iNODE4;
													}
													else {
														printf("this can not be B T7\n");
														//getchar();
														system("PAUSE");
													}
												}
											}
										}

									}

									// узел уже был посещён
									// Здесь неким образом модифицировалось gran_t
									if (return_current_node_number >= maxelm_memo) {

										printf("bvisit error B ...\n");
										printf("model is incorrect. FATALL error!!!\n");
										printf("internal flat heat source can not be located\n");
										printf("inside the liquid.\n");
										printf("Please, press any key to exit...\n");
										//getchar();
										system("PAUSE");
										//exit(1);

										/*
										// Это граничная грань внутреннего источника тепла.
										// Внутренний источник тепла.
										G = BSIDE;
										// граничная грань:
										INIT_SOSEDI(sosedi, maxelm, BSIDE, GG, false, sosedi[TSIDE][octree1->linkB->inum_FD - 1].iNODE1, NON_EXISTENT_NODE, NON_EXISTENT_NODE, NON_EXISTENT_NODE, false, NON_EXISTENT_NODE, NON_EXISTENT_NODE, NON_EXISTENT_NODE, NON_EXISTENT_NODE); // заполнение sosedi[G][maxelm].
										
										
											// FLUID
											prop_b[RHO][sosedi[G][maxelm].iNODE1 - maxelm_memo] = prop[RHO][octree1->linkB->inum_FD - 1];
											prop_b[CP][sosedi[G][maxelm].iNODE1 - maxelm_memo] = prop[CP][octree1->linkB->inum_FD - 1];
											prop_b[LAM][sosedi[G][maxelm].iNODE1 - maxelm_memo] = prop[LAM][octree1->linkB->inum_FD - 1];
											prop_b[MULT_LAM_X][sosedi[G][maxelm].iNODE1 - maxelm_memo] = prop[MULT_LAM_X][octree1->linkB->inum_FD - 1];
											prop_b[MULT_LAM_Y][sosedi[G][maxelm].iNODE1 - maxelm_memo] = prop[MULT_LAM_Y][octree1->linkB->inum_FD - 1];
											prop_b[MULT_LAM_Z][sosedi[G][maxelm].iNODE1 - maxelm_memo] = prop[MULT_LAM_Z][octree1->linkB->inum_FD - 1];
										

										doublereal dSloc = CALC_SQUARE(G, maxelm, nvtx, pa); // Вычисление площади грани ячейки.
										TOCHKA p_centerG = CALC_CENTERG(G, maxelm, nvtx, pa); // Вычисление геометрического центра грани.
										
										// G - грань по которой ставится граничное условие.
										// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
										// bound_id - номер граничного узла, начиная с нуля.
										// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
										if (octree1->b4T) {
											obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, sosedi[TSIDE][octree1->linkB->inum_FD - 1].iNODE1 - maxelm_memo,
												octree1->linkT4->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, true, true, p_centerG);
										}
										else {
											if (octree1->linkT == NULL) {
												obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, sosedi[TSIDE][octree1->linkB->inum_FD - 1].iNODE1 - maxelm_memo,
													NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, true, true, p_centerG);
											}
											else {
												obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, sosedi[TSIDE][octree1->linkB->inum_FD - 1].iNODE1 - maxelm_memo,
													octree1->linkT->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, true, true, p_centerG);
											}
										}

										*/

									}
									else {
										// Это строго внутренний узел.
										// Внутренняя грань.
										G = BSIDE;
										// Внутренняя грань:
										sosedi[G][maxelm].iNODE1 = octree1->linkB->inum_FD - 1; // граничные КО нумеруются в последнюю очередь,
										sosedi[G][maxelm].bdroblenie4 = false;
										sosedi[G][maxelm].iNODE2 = NON_EXISTENT_NODE;
										sosedi[G][maxelm].iNODE3 = NON_EXISTENT_NODE;
										sosedi[G][maxelm].iNODE4 = NON_EXISTENT_NODE;
										// после maxelm_memo внутренних КО.
										// Вычисление дальнего соседа
										CALC_GG(G, GG); // Вычисляет грань GG по грани G.

										// Определение дальнего соседа :
										// TODO 13 сентября 2016.
									}

								}
								else {
									bool binc;
									integer lsid;
									// узнать координаты центра грани и ориентацию в пространстве
									x_c = 0.5*(octree1->p0.x + octree1->p1.x);
									y_c = 0.5*(octree1->p0.y + octree1->p3.y);
									z_c = octree1->p0.z;
									iplane = XY;
									patch_maxbound(iplane, s, ls, x_c, y_c, z_c, maxbound, binc, lsid);
									if (binc) {

										printf("model is incorrect. FATALL error!!!\n");
										printf("internal flat heat source can not be located\n");
										printf("inside the liquid.\n");
										printf("Please, press any key to exit...\n");
										//getchar();
										system("PAUSE");
										exit(1);

										/*

										// Внутренний источник тепла.
										G = BSIDE;
										// граничная грань:
										INIT_SOSEDI(sosedi, maxelm, BSIDE, GG, false, maxelm_memo + maxbound-1, NON_EXISTENT_NODE, NON_EXISTENT_NODE, NON_EXISTENT_NODE, false, NON_EXISTENT_NODE, NON_EXISTENT_NODE, NON_EXISTENT_NODE, NON_EXISTENT_NODE); // заполнение sosedi[G][maxelm].

										
											// FLUID
											prop_b[RHO][maxbound - 1] = prop[RHO][octree1->linkB->inum_FD - 1];
											prop_b[CP][maxbound - 1] = prop[CP][octree1->linkB->inum_FD - 1];
											prop_b[LAM][maxbound - 1] = prop[LAM][octree1->linkB->inum_FD - 1];
											prop_b[MULT_LAM_X][maxbound - 1] = prop[MULT_LAM_X][octree1->linkB->inum_FD - 1];
											prop_b[MULT_LAM_Y][maxbound - 1] = prop[MULT_LAM_Y][octree1->linkB->inum_FD - 1];
											prop_b[MULT_LAM_Z][maxbound - 1] = prop[MULT_LAM_Z][octree1->linkB->inum_FD - 1];
										

										doublereal dSloc = CALC_SQUARE(G, maxelm, nvtx, pa); // Вычисление площади грани ячейки.
										TOCHKA p_centerG = CALC_CENTERG(G, maxelm, nvtx, pa); // Вычисление геометрического центра грани.
										
										// G - грань по которой ставится граничное условие.
										// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
										// bound_id - номер граничного узла, начиная с нуля.
										// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
										if (octree1->b4T) {
											obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound - 1,
												octree1->linkT4->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, true, false, p_centerG);
										}
										else {
											if (octree1->linkT == NULL) {
												obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound - 1,
													NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, true, false, p_centerG);
											}
											else {
												obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound - 1,
													octree1->linkT->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, true, false, p_centerG);
											}
										}

										*/

									}
									else {
										// Внутренняя грань.
										G = BSIDE;
										// Внутренняя грань:
										sosedi[G][maxelm].iNODE1 = octree1->linkB->inum_FD - 1; // граничные КО нумеруются в последнюю очередь,
										sosedi[G][maxelm].bdroblenie4 = false;
										sosedi[G][maxelm].iNODE2 = NON_EXISTENT_NODE;
										sosedi[G][maxelm].iNODE3 = NON_EXISTENT_NODE;
										sosedi[G][maxelm].iNODE4 = NON_EXISTENT_NODE;
										// после maxelm_memo внутренних КО.
										// Вычисление дальнего соседа
										CALC_GG(G, GG); // Вычисляет грань GG по грани G.

										// Определение дальнего соседа :
										// TODO 13 сентября 2016.
									}
								}
							}
						}
					}
					else {
					// Присутствуют все 4 соседа: B: 0,1,2,3.

						if (((octree1->linkB0 == NULL) && (octree1->linkB1 == NULL) && (octree1->linkB2 == NULL) && (octree1->linkB3 == NULL))||
						(((octree1->linkB0 != NULL) && (octree1->linkB1 != NULL) && (octree1->linkB2 != NULL) && (octree1->linkB3 != NULL) &&
							(octree1->linkB0->inum_FD == 0) && (octree1->linkB1->inum_FD == 0) && (octree1->linkB2->inum_FD == 0) && (octree1->linkB3->inum_FD == 0)))||
							(((octree1->linkB0 != NULL) && (octree1->linkB1 == NULL) && (octree1->linkB2 == NULL) && (octree1->linkB3 != NULL) &&
							(octree1->linkB0->inum_FD == 0) && (octree1->linkB3->inum_FD == 0)))||
								(((octree1->linkB0 != NULL) && (octree1->linkB1 != NULL) && (octree1->linkB2 == NULL) && (octree1->linkB3 == NULL) &&
							(octree1->linkB0->inum_FD == 0) && (octree1->linkB1->inum_FD == 0)))||
									(((octree1->linkB0 != NULL) && (octree1->linkB1 == NULL) && (octree1->linkB2 == NULL) && (octree1->linkB3 == NULL) &&
							(octree1->linkB0->inum_FD == 0))))
						{
							// 1 - B0&&B1&&B2&&B3 - на границе расчётной области (РО).
							// 2 - B0&&B1&&B2&&B3 - на гранеице HOLLOW блока.
							// 3 - B0&&B3 на границе HOLLOW блока, B1&&B2 на границе РО.
							// 4 - B0&&B1 на границе HOLLOW блока, B2&B3 на границе РО.
							// 5 - B0 на границе HOLLOW блока, B1&&B2&&B3  на границе РО.

							G = BSIDE;
							// граничная грань:
							INIT_SOSEDI(sosedi, maxelm, BSIDE, GG, false, maxelm_memo + maxbound, NON_EXISTENT_NODE, NON_EXISTENT_NODE, NON_EXISTENT_NODE, false, NON_EXISTENT_NODE, NON_EXISTENT_NODE, NON_EXISTENT_NODE, NON_EXISTENT_NODE); // заполнение sosedi[G][maxelm].
							DEFINE_BOUNDARY_PROPERTIES_FLOW(prop, prop_b, maxelm, maxbound); // В граничный узел сносятся свойства прилегающего КО.

							doublereal dSloc = CALC_SQUARE(G, maxelm, nvtx, pa); // Вычисление площади грани ячейки.
							TOCHKA p_centerG = CALC_CENTERG(G, maxelm, nvtx, pa); // Вычисление геометрического центра грани.
							
							// G - грань по которой ставится граничное условие.
							// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
							// bound_id - номер граничного узла, начиная с нуля.
							// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
							if (octree1->b4T) {
								obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
									octree1->linkT4->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, false, false, p_centerG);
							}
							else {
								if (octree1->linkT == NULL) {
									obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
										NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, false, false, p_centerG);
								}
								else {
									obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
										octree1->linkT->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dSloc, false, false, p_centerG);
								}
							}

							// Зачем дробить если граница заканчивается NULL.
							maxbound++; // Вместо четырёх один.
						}						
						else {

							bool binc1, binc2, binc3, binc4;
							integer lsid1, lsid2, lsid3, lsid4;

							if (octree1->linkB0 == NULL) {
#if doubleintprecision == 1
								printf("error NULL sosedi B0 elm=%lld\n", maxelm);
#else
								printf("error NULL sosedi B0 elm=%d\n", maxelm);
#endif
								
								//getchar();
								system("PAUSE");

								/*
								G = BSIDE;
								// граничная грань:
								sosedi[G][maxelm].iNODE1 = maxelm_memo + maxbound; // граничные КО нумеруются в последнюю очередь,
								sosedi[G][maxelm].bdroblenie4 = true;
								// после maxelm_memo внутренних КО.
								// Вычисление дальнего соседа
								CALC_GG(G,GG); // Вычисляет грань GG по грани G.

								sosedi[GG][maxelm].iNODE1 = NON_EXISTENT_NODE; // граничные КО нумеруются в последнюю очередь,
								sosedi[GG][maxelm].bdroblenie4 = true;

								#if doubleintprecision == 1
										printf("error NULL sosedi B0 elm=%lld\n", maxelm);
								#else
										printf("error NULL sosedi B0 elm=%d\n", maxelm);
								#endif
								
								getchar();

								// Эта ситуация возможна и здесь надо предусмотреть её обработку.
								DEFINE_BOUNDARY_PROPERTIES_FLOW(prop, prop_b, maxelm, maxbound); // В граничный узел сносятся свойства прилегающего КО.

								doublereal dS_loc = 0.0;

								// G - грань по которой ставится граничное условие.
								// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
								// bound_id - номер граничного узла, начиная с нуля.
								// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
								if (octree1->b4T) {
								obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
								octree1->linkT4->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
								}
								else {
								if (octree1->linkT == NULL) {
								obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
								NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
								}
								else {
								obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
								octree1->linkT->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
								}
								}

								maxbound++;
								*/
							}
							else {
								if (octree1->linkB0->inum_FD == 0) {

									G = BSIDE;
									// граничная грань:
									sosedi[G][maxelm].iNODE1 = maxelm_memo + maxbound; // граничные КО нумеруются в последнюю очередь,
									sosedi[G][maxelm].bdroblenie4 = true;
									// после maxelm_memo внутренних КО.
									// Вычисление дальнего соседа
									CALC_GG(G, GG); // Вычисляет грань GG по грани G.

									sosedi[GG][maxelm].iNODE1 = NON_EXISTENT_NODE; // граничные КО нумеруются в последнюю очередь,
									sosedi[GG][maxelm].bdroblenie4 = true;

#if doubleintprecision == 1
									//printf("error inum_FD == 0 sosedi B0 elm=%lld\n", maxelm);
#else
									//printf("error inum_FD == 0 sosedi B0 elm=%d\n", maxelm);
#endif
									
									//getchar();

									// Эта ситуация возможна и здесь надо предусмотреть её обработку.
									DEFINE_BOUNDARY_PROPERTIES_FLOW(prop, prop_b, maxelm, maxbound); // В граничный узел сносятся свойства прилегающего КО.

									// Sx*Sy
									doublereal dS_loc = fabs(octree1->linkB0->p7.y - octree1->linkB0->p4.y)*fabs(octree1->linkB0->p6.x - octree1->linkB0->p7.x);
									TOCHKA p_centerG;
									p_centerG.z = octree1->linkB0->p4.z;
									p_centerG.y = 0.5*(octree1->linkB0->p7.y + octree1->linkB0->p4.y);
									p_centerG.x = 0.5*(octree1->linkB0->p6.x + octree1->linkB0->p7.x);


									// G - грань по которой ставится граничное условие.
									// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
									// bound_id - номер граничного узла, начиная с нуля.
									// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
									if (octree1->b4T) {
										obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
											octree1->linkT4->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
									}
									else {
										if (octree1->linkT == NULL) {
											obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
												NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
										}
										else {
											obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
												octree1->linkT->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
										}
									}

									maxbound++;
								}
								else {
									// сосед существует.
									if (bvisit[octree1->linkB0->inum_FD - 1]) {
										// узел уже был посещён
										// Здесь неким образом модифицировалось gran_t
										// Только внутренняя грань.
										G = BSIDE;
										// граничная грань:
										sosedi[G][maxelm].iNODE1 = octree1->linkB0->inum_FD - 1; // граничные КО нумеруются в последнюю очередь,
										sosedi[G][maxelm].bdroblenie4 = true;
										// после maxelm_memo внутренних КО.
										// Вычисление дальнего соседа
										CALC_GG(G, GG); // Вычисляет грань GG по грани G.

										// Определение дальнего соседа :
										// TODO 13 сентября 2016.
									}
									else {
										// узнать координаты центра грани и ориентацию в пространстве
										x_c = 0.5*(octree1->linkB0->p4.x + octree1->linkB0->p5.x);
										y_c = 0.5*(octree1->linkB0->p4.y + octree1->linkB0->p7.y);
										z_c = octree1->p0.z;
										iplane = XY;
										patch_maxbound(iplane, s, ls, x_c, y_c, z_c, maxbound, binc1, lsid1);
										if (binc1) {
#if doubleintprecision == 1
											printf("error internal source B0 elm=%lld\n", maxelm);
#else
											printf("error internal source B0 elm=%d\n", maxelm);
#endif
											
											//getchar();
											system("PAUSE");
										}
										else {
											// Только внутренняя грань.
											G = BSIDE;
											// граничная грань:
											sosedi[G][maxelm].iNODE1 = octree1->linkB0->inum_FD - 1; // граничные КО нумеруются в последнюю очередь,
											sosedi[G][maxelm].bdroblenie4 = true;
											// после maxelm_memo внутренних КО.
											// Вычисление дальнего соседа
											CALC_GG(G, GG); // Вычисляет грань GG по грани G.

											// Определение дальнего соседа :
											// TODO 13 сентября 2016.
										}
									}
								}
							}

							if (octree1->linkB1 == NULL) {

								// 21 сентября 2016 (вырождение ячейки).

								// Заглушка.
								G = BSIDE;
								// граничная грань:
								sosedi[G][maxelm].iNODE2 = NON_EXISTENT_NODE; // граничные КО нумеруются в последнюю очередь,
								sosedi[G][maxelm].bdroblenie4 = true;
								// после maxelm_memo внутренних КО.
								// Вычисление дальнего соседа
								CALC_GG(G, GG); // Вычисляет грань GG по грани G.

								sosedi[GG][maxelm].iNODE2 = NON_EXISTENT_NODE; // сосед отсутствует,
								sosedi[GG][maxelm].bdroblenie4 = true;

								/*
								G = BSIDE;
								// граничная грань:
								sosedi[G][maxelm].iNODE2 = maxelm_memo + maxbound; // граничные КО нумеруются в последнюю очередь,
								sosedi[G][maxelm].bdroblenie4 = true;
								// после maxelm_memo внутренних КО.
								// Вычисление дальнего соседа
								CALC_GG(G,GG); // Вычисляет грань GG по грани G.

								sosedi[GG][maxelm].iNODE2 = NON_EXISTENT_NODE; // сосед отсутствует,
								sosedi[GG][maxelm].bdroblenie4 = true;

								#if doubleintprecision == 1
										printf("error NULL sosedi B1 elm=%lld\n", maxelm);
								#else
										printf("error NULL sosedi B1 elm=%d\n", maxelm);
								#endif
								
								getchar();

								// Эта ситуация возможна и здесь надо предусмотреть её обработку.
								DEFINE_BOUNDARY_PROPERTIES_FLOW(prop, prop_b, maxelm, maxbound); // В граничный узел сносятся свойства прилегающего КО.

								doublereal dS_loc = 0.0;

								// G - грань по которой ставится граничное условие.
								// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
								// bound_id - номер граничного узла, начиная с нуля.
								// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
								if (octree1->b4T) {
								obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
								octree1->linkT4->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
								}
								else {
								if (octree1->linkT == NULL) {
								obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
								NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
								}
								else {
								obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
								octree1->linkT->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
								}
								}

								maxbound++;
								*/
							}
							else {
								if (octree1->linkB1->inum_FD == 0) {

									G = BSIDE;
									// граничная грань:
									sosedi[G][maxelm].iNODE2 = maxelm_memo + maxbound; // граничные КО нумеруются в последнюю очередь,
									sosedi[G][maxelm].bdroblenie4 = true;
									// после maxelm_memo внутренних КО.
									// Вычисление дальнего соседа
									CALC_GG(G, GG); // Вычисляет грань GG по грани G.

									sosedi[GG][maxelm].iNODE2 = NON_EXISTENT_NODE; // сосед отсутствует,
									sosedi[GG][maxelm].bdroblenie4 = true;

#if doubleintprecision == 1
									//printf("error inum_FD == 0 sosedi B1 elm=%lld\n", maxelm);
#else
									//printf("error inum_FD == 0 sosedi B1 elm=%d\n", maxelm);
#endif
									
									//getchar();

									// Эта ситуация возможна и здесь надо предусмотреть её обработку.
									DEFINE_BOUNDARY_PROPERTIES_FLOW(prop, prop_b, maxelm, maxbound); // В граничный узел сносятся свойства прилегающего КО.

									// Sx*Sy
									doublereal dS_loc = fabs(octree1->linkB1->p7.y - octree1->linkB1->p4.y)*fabs(octree1->linkB1->p6.x - octree1->linkB1->p7.x);
									TOCHKA p_centerG;
									p_centerG.z = octree1->linkB1->p4.z;
									p_centerG.y = 0.5*(octree1->linkB1->p7.y + octree1->linkB1->p4.y);
									p_centerG.x = 0.5*(octree1->linkB1->p6.x + octree1->linkB1->p7.x);

									// G - грань по которой ставится граничное условие.
									// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
									// bound_id - номер граничного узла, начиная с нуля.
									// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
									if (octree1->b4T) {
										obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
											octree1->linkT4->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
									}
									else {
										if (octree1->linkT == NULL) {
											obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
												NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
										}
										else {
											obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
												octree1->linkT->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
										}
									}

									maxbound++;
								}
								else {
									// сосед существует.
									if (bvisit[octree1->linkB1->inum_FD - 1]) {
										// узел уже был посещён
										// Здесь неким образом модифицировалось gran_t
										// Только внутренняя грань.
										G = BSIDE;
										// граничная грань:
										sosedi[G][maxelm].iNODE2 = octree1->linkB1->inum_FD - 1; // граничные КО нумеруются в последнюю очередь,
										sosedi[G][maxelm].bdroblenie4 = true;
										// после maxelm_memo внутренних КО.
										// Вычисление дальнего соседа
										CALC_GG(G, GG); // Вычисляет грань GG по грани G.

										// Определение дальнего соседа :
										// TODO 13 сентября 2016.
									}
									else {
										// узнать координаты центра грани и ориентацию в пространстве
										x_c = 0.5*(octree1->linkB1->p4.x + octree1->linkB1->p5.x);
										y_c = 0.5*(octree1->linkB1->p4.y + octree1->linkB1->p7.y);
										z_c = octree1->p0.z;
										iplane = XY;
										patch_maxbound(iplane, s, ls, x_c, y_c, z_c, maxbound, binc2, lsid2);
										if (binc2) {
#if doubleintprecision == 1
											printf("error internal source B1 elm=%lld\n", maxelm);
#else
											printf("error internal source B1 elm=%d\n", maxelm);
#endif
											
											//getchar();
											system("PAUSE");
										}
										else {
											// Только внутренняя грань.
											G = BSIDE;
											// граничная грань:
											sosedi[G][maxelm].iNODE2 = octree1->linkB1->inum_FD - 1; // граничные КО нумеруются в последнюю очередь,
											sosedi[G][maxelm].bdroblenie4 = true;
											// после maxelm_memo внутренних КО.
											// Вычисление дальнего соседа
											CALC_GG(G, GG); // Вычисляет грань GG по грани G.

											// Определение дальнего соседа :
											// TODO 13 сентября 2016.
										}
									}
								}
							}

							if (octree1->linkB2 == NULL) {

								// 21 сентября 2016 (вырождение ячейки).

								// Заглушка.
								G = BSIDE;
								// граничная грань:
								sosedi[G][maxelm].iNODE3 = NON_EXISTENT_NODE; // граничные КО нумеруются в последнюю очередь,
								sosedi[G][maxelm].bdroblenie4 = true;
								// после maxelm_memo внутренних КО.
								// Вычисление дальнего соседа
								CALC_GG(G, GG); // Вычисляет грань GG по грани G.

								sosedi[GG][maxelm].iNODE3 = NON_EXISTENT_NODE; // граничные КО нумеруются в последнюю очередь,
								sosedi[GG][maxelm].bdroblenie4 = true;

								/*
								G = BSIDE;
								// граничная грань:
								sosedi[G][maxelm].iNODE3 = maxelm_memo + maxbound; // граничные КО нумеруются в последнюю очередь,
								sosedi[G][maxelm].bdroblenie4 = true;
								// после maxelm_memo внутренних КО.
								// Вычисление дальнего соседа
								CALC_GG(G,GG); // Вычисляет грань GG по грани G.

								sosedi[GG][maxelm].iNODE3 = NON_EXISTENT_NODE; // граничные КО нумеруются в последнюю очередь,
								sosedi[GG][maxelm].bdroblenie4 = true;

								#if doubleintprecision == 1
										printf("error NULL sosedi B2 elm=%lld\n", maxelm);
								#else
										printf("error NULL sosedi B2 elm=%d\n", maxelm);
								#endif
								
								getchar();

								// Эта ситуация возможна и здесь надо предусмотреть её обработку.
								DEFINE_BOUNDARY_PROPERTIES_FLOW(prop, prop_b, maxelm, maxbound); // В граничный узел сносятся свойства прилегающего КО.

								doublereal dS_loc = 0.0;

								// G - грань по которой ставится граничное условие.
								// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
								// bound_id - номер граничного узла, начиная с нуля.
								// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
								if (octree1->b4T) {
								obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
								octree1->linkT4->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
								}
								else {
								if (octree1->linkT == NULL) {
								obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
								NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
								}
								else {
								obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
								octree1->linkT->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
								}
								}


								maxbound++;
								*/
							}
							else {
								if (octree1->linkB2->inum_FD == 0) {

									G = BSIDE;
									// граничная грань:
									sosedi[G][maxelm].iNODE3 = maxelm_memo + maxbound; // граничные КО нумеруются в последнюю очередь,
									sosedi[G][maxelm].bdroblenie4 = true;
									// после maxelm_memo внутренних КО.
									// Вычисление дальнего соседа
									CALC_GG(G, GG); // Вычисляет грань GG по грани G.

									sosedi[GG][maxelm].iNODE3 = NON_EXISTENT_NODE; // граничные КО нумеруются в последнюю очередь,
									sosedi[GG][maxelm].bdroblenie4 = true;

#if doubleintprecision == 1
									//printf("error inum_FD == 0 sosedi B2 elm=%lld\n", maxelm);
#else
									//printf("error inum_FD == 0 sosedi B2 elm=%d\n", maxelm);
#endif
									
									//getchar();


									// Эта ситуация возможна и здесь надо предусмотреть её обработку.
									DEFINE_BOUNDARY_PROPERTIES_FLOW(prop, prop_b, maxelm, maxbound); // В граничный узел сносятся свойства прилегающего КО.

									// Sx*Sy
									doublereal dS_loc = fabs(octree1->linkB2->p7.y - octree1->linkB2->p4.y)*fabs(octree1->linkB2->p6.x - octree1->linkB2->p7.x);
									TOCHKA p_centerG;
									p_centerG.z = octree1->linkB2->p4.z;
									p_centerG.y = 0.5*(octree1->linkB2->p7.y + octree1->linkB2->p4.y);
									p_centerG.x = 0.5*(octree1->linkB2->p6.x + octree1->linkB2->p7.x);

									// G - грань по которой ставится граничное условие.
									// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
									// bound_id - номер граничного узла, начиная с нуля.
									// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
									if (octree1->b4T) {
										obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
											octree1->linkT4->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
									}
									else {
										if (octree1->linkT == NULL) {
											obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
												NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
										}
										else {
											obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
												octree1->linkT->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
										}
									}

									maxbound++;
								}
								else {
									// сосед существует.
									if (bvisit[octree1->linkB2->inum_FD - 1]) {
										// узел уже был посещён
										// Здесь неким образом модифицировалось gran_t
										// Только внутренняя грань.
										G = BSIDE;
										// граничная грань:
										sosedi[G][maxelm].iNODE3 = octree1->linkB2->inum_FD - 1; // граничные КО нумеруются в последнюю очередь,
										sosedi[G][maxelm].bdroblenie4 = true;
										// после maxelm_memo внутренних КО.
										// Вычисление дальнего соседа
										CALC_GG(G, GG); // Вычисляет грань GG по грани G.

										// Определение дальнего соседа :
										// TODO 13 сентября 2016.
									}
									else {
										// узнать координаты центра грани и ориентацию в пространстве
										x_c = 0.5*(octree1->linkB2->p4.x + octree1->linkB2->p5.x);
										y_c = 0.5*(octree1->linkB2->p4.y + octree1->linkB2->p7.y);
										z_c = octree1->p0.z;
										iplane = XY;
										patch_maxbound(iplane, s, ls, x_c, y_c, z_c, maxbound, binc3, lsid3);
										if (binc3) {
#if doubleintprecision == 1
											printf("error internal source B2 elm=%lld\n", maxelm);
#else
											printf("error internal source B2 elm=%d\n", maxelm);
#endif
											
											//getchar();
											system("PAUSE");
										}
										else {
											// Только внутренняя грань.
											G = BSIDE;
											// граничная грань:
											sosedi[G][maxelm].iNODE3 = octree1->linkB2->inum_FD - 1; // граничные КО нумеруются в последнюю очередь,
											sosedi[G][maxelm].bdroblenie4 = true;
											// после maxelm_memo внутренних КО.
											// Вычисление дальнего соседа
											CALC_GG(G, GG); // Вычисляет грань GG по грани G.

											// Определение дальнего соседа :
											// TODO 13 сентября 2016.
										}
									}
								}
							}

							if (octree1->linkB3 == NULL) {

								// 21 сентября 2016 (вырождение ячейки).

								// Заглушка.
								G = BSIDE;
								// граничная грань:
								sosedi[G][maxelm].iNODE4 = NON_EXISTENT_NODE; // граничные КО нумеруются в последнюю очередь,
								sosedi[G][maxelm].bdroblenie4 = true;
								// после maxelm_memo внутренних КО.
								// Вычисление дальнего соседа
								CALC_GG(G, GG); // Вычисляет грань GG по грани G.

								sosedi[GG][maxelm].iNODE4 = NON_EXISTENT_NODE; // сосед отсутствует,
								sosedi[GG][maxelm].bdroblenie4 = true;

								/*
								G = BSIDE;
								// граничная грань:
								sosedi[G][maxelm].iNODE4 = maxelm_memo + maxbound; // граничные КО нумеруются в последнюю очередь,
								sosedi[G][maxelm].bdroblenie4 = true;
								// после maxelm_memo внутренних КО.
								// Вычисление дальнего соседа
								CALC_GG(G,GG); // Вычисляет грань GG по грани G.

								sosedi[GG][maxelm].iNODE4 = NON_EXISTENT_NODE; // сосед отсутствует,
								sosedi[GG][maxelm].bdroblenie4 = true;

								#if doubleintprecision == 1
										printf("error NULL sosedi B3 elm=%lld\n", maxelm);
								#else
										printf("error NULL sosedi B3 elm=%d\n", maxelm);
								#endif
								
								getchar();

								// Эта ситуация возможна и здесь надо предусмотреть её обработку.
								DEFINE_BOUNDARY_PROPERTIES_FLOW(prop, prop_b, maxelm, maxbound); // В граничный узел сносятся свойства прилегающего КО.

								doublereal dS_loc = 0.0;

								// G - грань по которой ставится граничное условие.
								// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
								// bound_id - номер граничного узла, начиная с нуля.
								// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
								if (octree1->b4T) {
								obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
								octree1->linkT4->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
								}
								else {
								if (octree1->linkT == NULL) {
								obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
								NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
								}
								else {
								obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
								octree1->linkT->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
								}
								}

								maxbound++;
								*/
							}
							else {
								if (octree1->linkB3->inum_FD == 0) {


									G = BSIDE;
									// граничная грань:
									sosedi[G][maxelm].iNODE4 = maxelm_memo + maxbound; // граничные КО нумеруются в последнюю очередь,
									sosedi[G][maxelm].bdroblenie4 = true;
									// после maxelm_memo внутренних КО.
									// Вычисление дальнего соседа
									CALC_GG(G, GG); // Вычисляет грань GG по грани G.

									sosedi[GG][maxelm].iNODE4 = NON_EXISTENT_NODE; // сосед отсутствует,
									sosedi[GG][maxelm].bdroblenie4 = true;

#if doubleintprecision == 1
									//printf("error inum_FD == 0 sosedi B3 elm=%lld\n", maxelm);
#else
									//printf("error inum_FD == 0 sosedi B3 elm=%d\n", maxelm);
#endif
									
									//getchar();

									// Эта ситуация возможна и здесь надо предусмотреть её обработку.
									DEFINE_BOUNDARY_PROPERTIES_FLOW(prop, prop_b, maxelm, maxbound); // В граничный узел сносятся свойства прилегающего КО.

									// Sx*Sy
									doublereal dS_loc = fabs(octree1->linkB3->p7.y - octree1->linkB3->p4.y)*fabs(octree1->linkB3->p6.x - octree1->linkB3->p7.x);
									TOCHKA p_centerG;
									p_centerG.z = octree1->linkB3->p4.z;
									p_centerG.y = 0.5*(octree1->linkB3->p7.y + octree1->linkB3->p4.y);
									p_centerG.x = 0.5*(octree1->linkB3->p6.x + octree1->linkB3->p7.x);

									// G - грань по которой ставится граничное условие.
									// elm_id - номер конечного элемента, начиная с нуля для которго рассматривается граничный узел.
									// bound_id - номер граничного узла, начиная с нуля.
									// elm_id_inverse - это номер конечного элемента по другую сторону от границы для текущего элемента.
									if (octree1->b4T) {
										obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
											octree1->linkT4->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
									}
									else {
										if (octree1->linkT == NULL) {
											obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
												NON_EXISTENT_NODE, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
										}
										else {
											obrabotka_granichnoi_grani(G, sosedb, nvtx, pa, maxelm, maxbound,
												octree1->linkT->inum_FD - 1, maxelm_memo, whot_is_block, ls, lw, w, s, b, dS_loc, false, false, p_centerG);
										}
									}

									maxbound++;
								}
								else {
									// сосед существует.
									if (bvisit[octree1->linkB3->inum_FD - 1]) {
										// узел уже был посещён
										// Здесь неким образом модифицировалось gran_t
										// Только внутренняя грань.
										G = BSIDE;
										// граничная грань:
										sosedi[G][maxelm].iNODE4 = octree1->linkB3->inum_FD - 1; // граничные КО нумеруются в последнюю очередь,
										sosedi[G][maxelm].bdroblenie4 = true;
										// после maxelm_memo внутренних КО.
										// Вычисление дальнего соседа
										CALC_GG(G, GG); // Вычисляет грань GG по грани G.

										// Определение дальнего соседа :
										// TODO 13 сентября 2016.
									}
									else {
										// узнать координаты центра грани и ориентацию в пространстве
										x_c = 0.5*(octree1->linkB3->p4.x + octree1->linkB3->p5.x);
										y_c = 0.5*(octree1->linkB3->p4.y + octree1->linkB3->p7.y);
										z_c = octree1->p0.z;
										iplane = XY;
										patch_maxbound(iplane, s, ls, x_c, y_c, z_c, maxbound, binc4, lsid4);
										if (binc4) {
#if doubleintprecision == 1
											printf("error internal source B3 elm=%lld\n", maxelm);
#else
											printf("error internal source B3 elm=%d\n", maxelm);
#endif
											
											//getchar();
											system("PAUSE");
										}
										else {
											// Только внутренняя грань.
											G = BSIDE;
											// граничная грань:
											sosedi[G][maxelm].iNODE4 = octree1->linkB3->inum_FD - 1; // граничные КО нумеруются в последнюю очередь,
											sosedi[G][maxelm].bdroblenie4 = true;
											// после maxelm_memo внутренних КО.
											// Вычисление дальнего соседа
											CALC_GG(G, GG); // Вычисляет грань GG по грани G.

											// Определение дальнего соседа :
											// TODO 13 сентября 2016.
										}
									}
								}
							}
						}
					}


					bvisit[maxelm] = true; // Узел был посещён.
					maxelm++;
				}
				//octree1->inum_FD = 0; // По умолчанию не принадлежит расчётной области.
				octree1 = NULL;
				my_ALICE_STACK[top_ALICE_STACK - 1].link = NULL;
				top_ALICE_STACK--;
			}
			else {
				// продолжаем добираться до листьев.
				STACK_ALICE buf1 = my_ALICE_STACK[top_ALICE_STACK - 1];
				STACK_ALICE* buf = &buf1;
				top_ALICE_STACK--;
				if (buf->link->link0 != NULL) {
					my_ALICE_STACK[top_ALICE_STACK].link = (buf->link->link0);
					my_ALICE_STACK[top_ALICE_STACK].maxx = buf->link->link0->maxx;
					my_ALICE_STACK[top_ALICE_STACK].minx = buf->link->link0->minx;
					my_ALICE_STACK[top_ALICE_STACK].maxy = buf->link->link0->maxy;
					my_ALICE_STACK[top_ALICE_STACK].miny = buf->link->link0->miny;
					my_ALICE_STACK[top_ALICE_STACK].maxz = buf->link->link0->maxz;
					my_ALICE_STACK[top_ALICE_STACK].minz = buf->link->link0->minz;
					top_ALICE_STACK++;
				}
				if (buf->link->link1 != NULL) {
					my_ALICE_STACK[top_ALICE_STACK].link = (buf->link->link1);
					my_ALICE_STACK[top_ALICE_STACK].maxx = buf->link->link1->maxx;
					my_ALICE_STACK[top_ALICE_STACK].minx = buf->link->link1->minx;
					my_ALICE_STACK[top_ALICE_STACK].maxy = buf->link->link1->maxy;
					my_ALICE_STACK[top_ALICE_STACK].miny = buf->link->link1->miny;
					my_ALICE_STACK[top_ALICE_STACK].maxz = buf->link->link1->maxz;
					my_ALICE_STACK[top_ALICE_STACK].minz = buf->link->link1->minz;
					top_ALICE_STACK++;
				}
				if (buf->link->link2 != NULL) {
					my_ALICE_STACK[top_ALICE_STACK].link = (buf->link->link2);
					my_ALICE_STACK[top_ALICE_STACK].maxx = buf->link->link2->maxx;
					my_ALICE_STACK[top_ALICE_STACK].minx = buf->link->link2->minx;
					my_ALICE_STACK[top_ALICE_STACK].maxy = buf->link->link2->maxy;
					my_ALICE_STACK[top_ALICE_STACK].miny = buf->link->link2->miny;
					my_ALICE_STACK[top_ALICE_STACK].maxz = buf->link->link2->maxz;
					my_ALICE_STACK[top_ALICE_STACK].minz = buf->link->link2->minz;
					top_ALICE_STACK++;
				}
				if (buf->link->link3 != NULL) {
					my_ALICE_STACK[top_ALICE_STACK].link = (buf->link->link3);
					my_ALICE_STACK[top_ALICE_STACK].maxx = buf->link->link3->maxx;
					my_ALICE_STACK[top_ALICE_STACK].minx = buf->link->link3->minx;
					my_ALICE_STACK[top_ALICE_STACK].maxy = buf->link->link3->maxy;
					my_ALICE_STACK[top_ALICE_STACK].miny = buf->link->link3->miny;
					my_ALICE_STACK[top_ALICE_STACK].maxz = buf->link->link3->maxz;
					my_ALICE_STACK[top_ALICE_STACK].minz = buf->link->link3->minz;
					top_ALICE_STACK++;
				}
				if (buf->link->link4 != NULL) {
					my_ALICE_STACK[top_ALICE_STACK].link = (buf->link->link4);
					my_ALICE_STACK[top_ALICE_STACK].maxx = buf->link->link4->maxx;
					my_ALICE_STACK[top_ALICE_STACK].minx = buf->link->link4->minx;
					my_ALICE_STACK[top_ALICE_STACK].maxy = buf->link->link4->maxy;
					my_ALICE_STACK[top_ALICE_STACK].miny = buf->link->link4->miny;
					my_ALICE_STACK[top_ALICE_STACK].maxz = buf->link->link4->maxz;
					my_ALICE_STACK[top_ALICE_STACK].minz = buf->link->link4->minz;
					top_ALICE_STACK++;
				}
				if (buf->link->link5 != NULL) {
					my_ALICE_STACK[top_ALICE_STACK].link = (buf->link->link5);
					my_ALICE_STACK[top_ALICE_STACK].maxx = buf->link->link5->maxx;
					my_ALICE_STACK[top_ALICE_STACK].minx = buf->link->link5->minx;
					my_ALICE_STACK[top_ALICE_STACK].maxy = buf->link->link5->maxy;
					my_ALICE_STACK[top_ALICE_STACK].miny = buf->link->link5->miny;
					my_ALICE_STACK[top_ALICE_STACK].maxz = buf->link->link5->maxz;
					my_ALICE_STACK[top_ALICE_STACK].minz = buf->link->link5->minz;
					top_ALICE_STACK++;
				}
				if (buf->link->link6 != NULL) {
					my_ALICE_STACK[top_ALICE_STACK].link = (buf->link->link6);
					my_ALICE_STACK[top_ALICE_STACK].maxx = buf->link->link6->maxx;
					my_ALICE_STACK[top_ALICE_STACK].minx = buf->link->link6->minx;
					my_ALICE_STACK[top_ALICE_STACK].maxy = buf->link->link6->maxy;
					my_ALICE_STACK[top_ALICE_STACK].miny = buf->link->link6->miny;
					my_ALICE_STACK[top_ALICE_STACK].maxz = buf->link->link6->maxz;
					my_ALICE_STACK[top_ALICE_STACK].minz = buf->link->link6->minz;
					top_ALICE_STACK++;
				}
				if (buf->link->link7 != NULL) {
					my_ALICE_STACK[top_ALICE_STACK].link = (buf->link->link7);
					my_ALICE_STACK[top_ALICE_STACK].maxx = buf->link->link7->maxx;
					my_ALICE_STACK[top_ALICE_STACK].minx = buf->link->link7->minx;
					my_ALICE_STACK[top_ALICE_STACK].maxy = buf->link->link7->maxy;
					my_ALICE_STACK[top_ALICE_STACK].miny = buf->link->link7->miny;
					my_ALICE_STACK[top_ALICE_STACK].maxz = buf->link->link7->maxz;
					my_ALICE_STACK[top_ALICE_STACK].minz = buf->link->link7->minz;
					top_ALICE_STACK++;
				}
			}
		}
		//}
		//getchar();
	}

	delete[] bvisit;
	bvisit = NULL;
} // constr_sosedi_prop_b_flow_alice

#endif