// Файл avtosave.c автосохранение и возобновление расчёта
// с прерванного места.
// начало 16 октября 2011 г.

#ifndef AVTO_SAVE_C
#define AVTO_SAVE_C 1

// Автосохранение рассчитанных полевых величин в файл avtosave.txt
// для возобновления счёта после прерывания.
void avtosave(FLOW* &f, TEMPER &t, int flow_interior, int* &inumber_iteration_SIMPLE, doublereal* &continity_start) {

	// flow_interior - общее число независимых FLUID зон.
	// inumber_iteration_SIMPLE - номер итерации SIMPLE алгоритма на момент которой происходит запись в файл (для каждой жидкой зоны).
	// continity_start - начальное значение невязки continity (для каждой жидкой зоны).


	FILE *fp_avtosave=NULL; // файл в который будут записываться вычисленные поля
	
#ifdef MINGW_COMPILLER
	int err = 0;
	fp_avtosave=fopen64("avtosave.txt", "w");
	if (fp_avtosave == NULL) {
		err = 1;
	}
#else
	errno_t err = 0;
	err = fopen_s(&fp_avtosave, "avtosave.txt", "w");
#endif
	
	if (((err) !=0)||(fp_avtosave == NULL)) {
		// Файл открывается для полной перезаписи,
		// если такого файла раньше не было то он будет создан.
	   printf("Create File avtosave.txt Error\n");
       //system("pause");
	   system("pause");
       exit(0);
	}
	else {

#if doubleintprecision == 1
		fprintf(fp_avtosave, "%d %d ", t.maxelm + t.maxbound, flow_interior);
		integer i = 0; // счётчик жидких FLUID зон
		for (i = 0; i<flow_interior; ++i) fprintf(fp_avtosave, "%d ", f[i].maxelm + f[i].maxbound);
		fprintf(fp_avtosave, "\n");

		for (i = 0; i<flow_interior; ++i) fprintf(fp_avtosave, "%d ", inumber_iteration_SIMPLE[i]);
		for (i = 0; i<flow_interior; ++i) fprintf(fp_avtosave, "%+.16f ", continity_start[i]);
		fprintf(fp_avtosave, "\n"); // переход к новой строке
#else
		fprintf(fp_avtosave, "%d %d ", t.maxelm + t.maxbound, flow_interior);
		integer i = 0; // счётчик жидких FLUID зон
		for (i = 0; i<flow_interior; ++i) fprintf(fp_avtosave, "%d ", f[i].maxelm + f[i].maxbound);
		fprintf(fp_avtosave, "\n");

		for (i = 0; i<flow_interior; ++i) fprintf(fp_avtosave, "%d ", inumber_iteration_SIMPLE[i]);
		for (i = 0; i<flow_interior; ++i) fprintf(fp_avtosave, "%+.16f ", continity_start[i]);
		fprintf(fp_avtosave, "\n"); // переход к новой строке
#endif

		

		// TEMP
		for (i=0; i<t.maxelm+t.maxbound; ++i) {
           fprintf(fp_avtosave,"%+.16f ",t.potent[i]);
		   if ((i+1)%20==0) fprintf(fp_avtosave,"\n");
		}
		if ((t.maxelm+t.maxbound)%20 != 0) fprintf(fp_avtosave,"\n");

		integer j=0;
		for (j=0; j<flow_interior; ++j) {
			// VX
			for (i=0; i<(f[j].maxelm+f[j].maxbound); ++i) {
                fprintf(fp_avtosave,"%+.16f ",f[j].potent[VELOCITY_X_COMPONENT][i]);
		        if ((i+1)%20==0) fprintf(fp_avtosave,"\n");
		    }
		    if ((f[j].maxelm+f[j].maxbound)%20 != 0) fprintf(fp_avtosave,"\n");
			// VY
            for (i=0; i<(f[j].maxelm+f[j].maxbound); ++i) {
                fprintf(fp_avtosave,"%+.16f ",f[j].potent[VELOCITY_Y_COMPONENT][i]);
		        if ((i+1)%20==0) fprintf(fp_avtosave,"\n");
		    }
		    if ((f[j].maxelm+f[j].maxbound)%20 != 0) fprintf(fp_avtosave,"\n");
			// VZ
			for (i=0; i<(f[j].maxelm+f[j].maxbound); ++i) {
                fprintf(fp_avtosave,"%+.16f ",f[j].potent[VELOCITY_Z_COMPONENT][i]);
		        if ((i+1)%20==0) fprintf(fp_avtosave,"\n");
		    }
		    if ((f[j].maxelm+f[j].maxbound)%20 != 0) fprintf(fp_avtosave,"\n");
			// PRESSURE
			for (i=0; i<(f[j].maxelm+f[j].maxbound); ++i) {
                fprintf(fp_avtosave,"%+.16f ",f[j].potent[PRESS][i]);
		        if ((i+1)%20==0) fprintf(fp_avtosave,"\n");
		    }
		    if ((f[j].maxelm+f[j].maxbound)%20 != 0) fprintf(fp_avtosave,"\n");

		} // конец цикла по всем жидким зонам

		fclose(fp_avtosave); // закрытие файла для записи рассчитанных полевых величин.
	}


} // avtosave

// возобновление счёта после прерывания:
// считывание значений полевых величин из файла avtosave.txt 
// при его наличии и корректности.
void avtoreadvalue(FLOW* &f, TEMPER &t, int flow_interior,
	int* &inumber_iteration_SIMPLE, doublereal* &continity_start, bool &breadOk,
	BLOCK* &b, int &lb, SOURCE* &s_loc, int &ls, WALL* &w, int &lw) {

	// flow_interior - общее число независимых FLUID зон.
	// inumber_iteration_SIMPLE - номер итерации SIMPLE алгоритма на момент которой происходит запись в файл (для каждой жидкой зоны).
	// continity_start - начальное значение невязки continity (для каждой жидкой зоны).

	FILE *fp_avtosave=NULL; // файл в котором записаны сохранённые полевые величины.
	
#ifdef MINGW_COMPILLER
	int err1 = 0;
	fp_avtosave=fopen64("avtosave.txt", "r");
	if (fp_avtosave == NULL) err1 = 1;
#else
	errno_t err1 = 0;
	err1 = fopen_s(&fp_avtosave, "avtosave.txt", "r");
#endif

	if ((err1)!=0) {
		// Файла автосохранения изначально просто нет,
		// значит счёт начнётся от печки.
		printf("No input File avtosave.txt \n");
		breadOk=false;
		//system("pause");
	}
	else if (fp_avtosave!=NULL)
	{
		bool bcontinue=true; // произвести ли считывание полей

		//float fin=0.0;
		doublereal fin=0.0;
		integer din=0;

		// Сначала небольшая проверка соответствия файлов:
		// premeshin.txt и avtosave.txt
		// bcontinue==true если проверка прошла успешно.
#ifdef MINGW_COMPILLER
#if doubleintprecision == 1
		fscanf(fp_avtosave, "%lld", &din);
#else
		fscanf(fp_avtosave, "%d", &din);
#endif
#else
#if doubleintprecision == 1
		fscanf_s(fp_avtosave, "%lld", &din);
#else
		fscanf_s(fp_avtosave, "%d", &din);
#endif
#endif
        
		if (din != (t.maxelm+t.maxbound)) bcontinue=false; // проверка на совпадение числа КО в задаче теплопроводности

#ifdef MINGW_COMPILLER
#if doubleintprecision == 1
		fscanf(fp_avtosave, "%lld", &din);
#else
		fscanf(fp_avtosave, "%d", &din);
#endif
#else
#if doubleintprecision == 1
		fscanf_s(fp_avtosave, "%lld", &din);
#else
		fscanf_s(fp_avtosave, "%d", &din);
#endif
#endif

		
		if (din != flow_interior) bcontinue=false; // проверка на количество жидких FLUID зон

		integer i=0; // счётчик FLUID зон
		for (i=0; i<flow_interior; ++i) {

#ifdef MINGW_COMPILLER
#if doubleintprecision == 1
			fscanf(fp_avtosave, "%lld", &din);
#else
			fscanf(fp_avtosave, "%d", &din);
#endif
#else
#if doubleintprecision == 1
			fscanf_s(fp_avtosave, "%lld", &din);
#else
			fscanf_s(fp_avtosave, "%d", &din);
#endif
#endif

			
			if (din != (f[i].maxelm+f[i].maxbound)) bcontinue=false; // проверка на количество КО в FLUID зоне номер i.
		}

		if (!bcontinue) {
			// файл premeshin.txt не соответствует файлу avtosave.txt
			printf("premeshin.txt does not correspond to avtosave.txt...\n");
			// вычисление начнётся заново
			printf("computation begins a new! Please press any key to continue...\n");
			breadOk=false;
			//system("pause");
			system("pause");
		}
		else {

			for (i=0; i<flow_interior; ++i) {

#ifdef MINGW_COMPILLER
#if doubleintprecision == 1
				fscanf(fp_avtosave, "%lld", &din);
#else
				fscanf(fp_avtosave, "%d", &din);
#endif
#else
#if doubleintprecision == 1
				fscanf_s(fp_avtosave, "%lld", &din);
#else
				fscanf_s(fp_avtosave, "%d", &din);
#endif
#endif

				
				inumber_iteration_SIMPLE[i]=(int)(din);
			}
			for (i=0; i<flow_interior; ++i) {

#ifdef MINGW_COMPILLER
				fscanf(fp_avtosave, "%lf", &fin);
#else
				fscanf_s(fp_avtosave, "%lf", &fin);
#endif

				
				continity_start[i]=fin;
			}

			integer j=0; // счётчик КО.
			// единое поле температур
			printf("read TEMP...\n");
			for (j=0; j<(t.maxelm+t.maxbound); ++j) {

#ifdef MINGW_COMPILLER
				fscanf(fp_avtosave, "%lf", &fin);
#else
				fscanf_s(fp_avtosave, "%lf", &fin);
#endif

               
				t.potent[j]=fin;
			}
			// Цикл по всем жидким зонам:
			for (i=0; i<flow_interior; ++i) {
#if doubleintprecision == 1
				printf("read FLUID%lld...\n", i);
#else
				printf("read FLUID%d...\n", i);
#endif
				
				// VX
				printf("read VX...\n");
				for (j=0; j<(f[i].maxelm+f[i].maxbound); ++j) {

#ifdef MINGW_COMPILLER
					fscanf(fp_avtosave, "%lf", &fin);
#else
					fscanf_s(fp_avtosave, "%lf", &fin);
#endif

					
					f[i].potent[VELOCITY_X_COMPONENT][j]=fin;
					f[i].potent[VXCOR][j]=fin;
				}
				// VY
				printf("read VY...\n");
				for (j=0; j<(f[i].maxelm+f[i].maxbound); ++j) {

#ifdef MINGW_COMPILLER
					fscanf(fp_avtosave, "%lf", &fin);
#else
					fscanf_s(fp_avtosave, "%lf", &fin);
#endif

					
					f[i].potent[VELOCITY_Y_COMPONENT][j]=fin;
					f[i].potent[VYCOR][j]=fin;
				}
				// VZ
				printf("read VZ...\n");
				for (j=0; j<(f[i].maxelm+f[i].maxbound); ++j) {

#ifdef MINGW_COMPILLER
					fscanf(fp_avtosave, "%lf", &fin);
#else
					fscanf_s(fp_avtosave, "%lf", &fin);
#endif

					
					f[i].potent[VELOCITY_Z_COMPONENT][j]=fin;
					f[i].potent[VZCOR][j]=fin;
				}
				// PRESSURE
				printf("read PRESSURE...\n");
				for (j=0; j<(f[i].maxelm+f[i].maxbound); ++j) {

#ifdef MINGW_COMPILLER
					fscanf(fp_avtosave, "%lf", &fin);
#else
					fscanf_s(fp_avtosave, "%lf", &fin);
#endif

					
					f[i].potent[PRESS][j]=fin;
				}

				// Вычисление градиентов скорости:
				for (j=0; j<f[i].maxelm; ++j) {
					// градиенты скоростей для внутренних КО.
				    green_gauss(j, f[i].potent, f[i].nvtx, f[i].pa, f[i].neighbors_for_the_internal_node, f[i].maxelm, false,f[i], f[i].border_neighbor, t.ilevel_alice);
				}
				for (j=0; j<f[i].maxelm; ++j) {
					// градиенты скоростей для граничных КО.
				    green_gauss(j, f[i].potent, f[i].nvtx, f[i].pa, f[i].neighbors_for_the_internal_node, f[i].maxelm, true, f[i], f[i].border_neighbor, t.ilevel_alice);
				}

				// Вычисление градиентов давления:
	            // на основе скорректированного поля давления.
            	// Градиенты давления понадобятся при вычислении поправки Рхи-Чоу 1983г.
				for (integer j_1 = 0; j_1 < f[i].maxelm; ++j_1) {
					// градиенты давления для внутренних КО.
					green_gaussPRESS(j_1, f[i].potent, f[i].nvtx, f[i].pa, f[i].neighbors_for_the_internal_node, f[i].maxelm, false, f[i].border_neighbor, ls, lw, w, f[i].bLR1free, t.ilevel_alice, f[i].ptr, f[i].volume);
				}
				for (integer j_1 = 0; j_1 < f[i].maxelm; ++j_1) {
					// градиенты давления для граничных КО.
					green_gaussPRESS(j_1, f[i].potent, f[i].nvtx, f[i].pa, f[i].neighbors_for_the_internal_node, f[i].maxelm, true, f[i].border_neighbor, ls, lw, w, f[i].bLR1free, t.ilevel_alice, f[i].ptr, f[i].volume);
				}



				// На твёрдой стенке турбулентная динамическая вязкость равна нулю.
	            // Вычисление S инварианта тензора скоростей-деформаций для всех
	            // внутренних и граничных контрольных объёмов расчётной области.
	            #pragma omp parallel for shared (f) schedule (guided)
		        for (integer j_1=0; j_1<(f[i].maxelm+f[i].maxbound); ++j_1) {
			        // по поводу правильности формулы см. user_manual.
			        doublereal sum=0.0;
			        sum+=2.0*f[i].potent[GRADXVX][j_1]*f[i].potent[GRADXVX][j_1];
			        sum+=2.0*f[i].potent[GRADYVY][j_1]*f[i].potent[GRADYVY][j_1];
			        sum+=2.0*f[i].potent[GRADZVZ][j_1]*f[i].potent[GRADZVZ][j_1];
			        sum+=(f[i].potent[GRADYVX][j_1]+f[i].potent[GRADXVY][j_1])*(f[i].potent[GRADYVX][j_1]+f[i].potent[GRADXVY][j_1]);
			        sum+=(f[i].potent[GRADZVX][j_1]+f[i].potent[GRADXVZ][j_1])*(f[i].potent[GRADZVX][j_1]+f[i].potent[GRADXVZ][j_1]);
			        sum+=(f[i].potent[GRADYVZ][j_1]+f[i].potent[GRADZVY][j_1])*(f[i].potent[GRADYVZ][j_1]+f[i].potent[GRADZVY][j_1]);
					// следующее слагаемое может сделать подкоренное выражение отрицательным.
		            // вычитаем две трети квадрата дивергенции.
					sum-=(2.0/3.0)*(f[i].potent[GRADXVX][j_1]+f[i].potent[GRADYVY][j_1]+f[i].potent[GRADZVZ][j_1])*(f[i].potent[GRADXVX][j_1]+f[i].potent[GRADYVY][j_1]+f[i].potent[GRADZVZ][j_1]); // добавок связанный с несжимаемостью/сжимаемостью
			        f[i].SInvariantStrainRateTensor[j_1]=sqrt(fmax(0.0,sum));

                    // Вихрь (модуль ротора скорости).
					// по моему это тоже самое что и модуль тензора вращения.
		            sum=0.0;
		            sum+=(f[i].potent[GRADYVZ][j_1]-f[i].potent[GRADZVY][j_1])*(f[i].potent[GRADYVZ][j_1]-f[i].potent[GRADZVY][j_1]); // проверено !
		            sum+=(f[i].potent[GRADZVX][j_1]-f[i].potent[GRADXVZ][j_1])*(f[i].potent[GRADZVX][j_1]-f[i].potent[GRADXVZ][j_1]);
		            sum+=(f[i].potent[GRADXVY][j_1]-f[i].potent[GRADYVX][j_1])*(f[i].potent[GRADXVY][j_1]-f[i].potent[GRADYVX][j_1]);
		            f[i].potent[CURL][j_1]=sqrt(sum);
		        }

				breadOk=true; // считывание прошло успешно.

			} // конец цикла по жидким зонам
			
		    

		}

		

		fclose(fp_avtosave); // закрытие файла из которого считывались сохранённые полевые величины
	}

} // avtoreadvalue


#endif