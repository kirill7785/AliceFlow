// Файл avtosave.c автосохранение и возобновление расчёта
// с прерванного места.
// начало 16 октября 2011 г.

#ifndef AVTO_SAVE_C
#define AVTO_SAVE_C 1

// Автосохранение расчитанных полевых величин в файл avtosave.txt
// для возобновления счёта после прерывания.
void avtosave(FLOW* &f, TEMPER &t, integer flow_interior, integer* &inumber_iteration_SIMPLE, doublereal* &continity_start) {

	// flow_interior - общее число независимых FLUID зон.
	// inumber_iteration_SIMPLE - номер итерации SIMPLE алгоритма на момент которой происходит запись в файл (для каждой жидкой зоны).
	// continity_start - начальное значение невязки continity (для каждой жидкой зоны).


	FILE *fp_avtosave; // файл в который будут записываться вычисленные поля
	errno_t err;
	err = fopen_s(&fp_avtosave, "avtosave.txt", "w");
	
	if ((err) !=0) {
		// Файл открывается для полной перезаписи,
		// если такого файла раньше не было то он будет создан.
	   printf("Create File avtosave.txt Error\n");
       //getchar();
	   system("pause");
       exit(0);
	}
	else {

#if doubleintprecision == 1
		fprintf(fp_avtosave, "%lld %lld ", t.maxelm + t.maxbound, flow_interior);
		integer i = 0; // счётчик жидких FLUID зон
		for (i = 0; i<flow_interior; i++) fprintf(fp_avtosave, "%lld ", f[i].maxelm + f[i].maxbound);
		fprintf(fp_avtosave, "\n");

		for (i = 0; i<flow_interior; i++) fprintf(fp_avtosave, "%lld ", inumber_iteration_SIMPLE[i]);
		for (i = 0; i<flow_interior; i++) fprintf(fp_avtosave, "%+.16f ", continity_start[i]);
		fprintf(fp_avtosave, "\n"); // переход к новой строке
#else
		fprintf(fp_avtosave, "%d %d ", t.maxelm + t.maxbound, flow_interior);
		integer i = 0; // счётчик жидких FLUID зон
		for (i = 0; i<flow_interior; i++) fprintf(fp_avtosave, "%d ", f[i].maxelm + f[i].maxbound);
		fprintf(fp_avtosave, "\n");

		for (i = 0; i<flow_interior; i++) fprintf(fp_avtosave, "%d ", inumber_iteration_SIMPLE[i]);
		for (i = 0; i<flow_interior; i++) fprintf(fp_avtosave, "%+.16f ", continity_start[i]);
		fprintf(fp_avtosave, "\n"); // переход к новой строке
#endif

		

		// TEMP
		for (i=0; i<t.maxelm+t.maxbound; i++) {
           fprintf(fp_avtosave,"%+.16f ",t.potent[i]);
		   if ((i+1)%20==0) fprintf(fp_avtosave,"\n");
		}
		if ((t.maxelm+t.maxbound)%20 != 0) fprintf(fp_avtosave,"\n");

		integer j=0;
		for (j=0; j<flow_interior; j++) {
			// VX
			for (i=0; i<(f[j].maxelm+f[j].maxbound); i++) {
                fprintf(fp_avtosave,"%+.16f ",f[j].potent[VX][i]);
		        if ((i+1)%20==0) fprintf(fp_avtosave,"\n");
		    }
		    if ((f[j].maxelm+f[j].maxbound)%20 != 0) fprintf(fp_avtosave,"\n");
			// VY
            for (i=0; i<(f[j].maxelm+f[j].maxbound); i++) {
                fprintf(fp_avtosave,"%+.16f ",f[j].potent[VY][i]);
		        if ((i+1)%20==0) fprintf(fp_avtosave,"\n");
		    }
		    if ((f[j].maxelm+f[j].maxbound)%20 != 0) fprintf(fp_avtosave,"\n");
			// VZ
			for (i=0; i<(f[j].maxelm+f[j].maxbound); i++) {
                fprintf(fp_avtosave,"%+.16f ",f[j].potent[VZ][i]);
		        if ((i+1)%20==0) fprintf(fp_avtosave,"\n");
		    }
		    if ((f[j].maxelm+f[j].maxbound)%20 != 0) fprintf(fp_avtosave,"\n");
			// PRESSURE
			for (i=0; i<(f[j].maxelm+f[j].maxbound); i++) {
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
void avtoreadvalue(FLOW* &f, TEMPER &t, integer flow_interior, integer* &inumber_iteration_SIMPLE, doublereal* &continity_start, bool &breadOk) {

	// flow_interior - общее число независимых FLUID зон.
	// inumber_iteration_SIMPLE - номер итерации SIMPLE алгоритма на момент которой происходит запись в файл (для каждой жидкой зоны).
	// continity_start - начальное значение невязки continity (для каждой жидкой зоны).

	FILE *fp_avtosave=NULL; // файл в котором записаны сохранённые полевые величины.
	errno_t err1;
	err1 = fopen_s(&fp_avtosave, "avtosave.txt", "r");

	if ((err1)!=0) {
		// Файла автосохранения изначально просто нет,
		// значит счёт начнётся от печки.
		printf("No input File avtosave.txt \n");
		breadOk=false;
		//getchar();
	}
	else if (fp_avtosave!=NULL)
	{
		bool bcontinue=true; // произвести ли считывание полей

		//float fin=0.0;
		doublereal fin=0.0;
		integer din=0;

		// Сначала небольшая проверка соответствия файлов :
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
		for (i=0; i<flow_interior; i++) {

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
			//getchar();
			system("pause");
		}
		else {

			for (i=0; i<flow_interior; i++) {

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

				
				inumber_iteration_SIMPLE[i]=din;
			}
			for (i=0; i<flow_interior; i++) {

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
			for (j=0; j<(t.maxelm+t.maxbound); j++) {

#ifdef MINGW_COMPILLER
				fscanf(fp_avtosave, "%lf", &fin);
#else
				fscanf_s(fp_avtosave, "%lf", &fin);
#endif

               
				t.potent[j]=fin;
			}
			// Цикл по всем жидким зонам:
			for (i=0; i<flow_interior; i++) {
#if doubleintprecision == 1
				printf("read FLUID%lld...\n", i);
#else
				printf("read FLUID%d...\n", i);
#endif
				
				// VX
				printf("read VX...\n");
				for (j=0; j<(f[i].maxelm+f[i].maxbound); j++) {

#ifdef MINGW_COMPILLER
					fscanf(fp_avtosave, "%lf", &fin);
#else
					fscanf_s(fp_avtosave, "%lf", &fin);
#endif

					
					f[i].potent[VX][j]=fin;
					f[i].potent[VXCOR][j]=fin;
				}
				// VY
				printf("read VY...\n");
				for (j=0; j<(f[i].maxelm+f[i].maxbound); j++) {

#ifdef MINGW_COMPILLER
					fscanf(fp_avtosave, "%lf", &fin);
#else
					fscanf_s(fp_avtosave, "%lf", &fin);
#endif

					
					f[i].potent[VY][j]=fin;
					f[i].potent[VYCOR][j]=fin;
				}
				// VZ
				printf("read VZ...\n");
				for (j=0; j<(f[i].maxelm+f[i].maxbound); j++) {

#ifdef MINGW_COMPILLER
					fscanf(fp_avtosave, "%lf", &fin);
#else
					fscanf_s(fp_avtosave, "%lf", &fin);
#endif

					
					f[i].potent[VZ][j]=fin;
					f[i].potent[VZCOR][j]=fin;
				}
				// PRESSURE
				printf("read PRESSURE...\n");
				for (j=0; j<(f[i].maxelm+f[i].maxbound); j++) {

#ifdef MINGW_COMPILLER
					fscanf(fp_avtosave, "%lf", &fin);
#else
					fscanf_s(fp_avtosave, "%lf", &fin);
#endif

					
					f[i].potent[PRESS][j]=fin;
				}

				// Вычисление градиентов скорости :
				for (j=0; j<f[i].maxelm; j++) {
					// градиенты скоростей для внутренних КО.
				    green_gauss(j, f[i].potent, f[i].nvtx, f[i].pa, f[i].sosedi, f[i].maxelm, false,f[i]);
				}
				for (j=0; j<f[i].maxelm; j++) {
					// градиенты скоростей для граничных КО.
				    green_gauss(j, f[i].potent, f[i].nvtx, f[i].pa, f[i].sosedi, f[i].maxelm, true, f[i]);
				}

				// На твёрдой стенке турбулентная динамическая вязкость равна нулю.
	            // Вычисление S инварианта тензора скоростей-деформаций для всех
	            // внутренних и граничных контрольных объёмов расчётной области.
	            #pragma omp parallel for shared (f) private (j) schedule (guided)
		        for (j=0; j<(f[i].maxelm+f[i].maxbound); j++) {
			        // по поводу правильности формулы см. user_manual.
			        doublereal sum=0.0;
			        sum+=2.0*f[i].potent[GRADXVX][j]*f[i].potent[GRADXVX][j];
			        sum+=2.0*f[i].potent[GRADYVY][j]*f[i].potent[GRADYVY][j];
			        sum+=2.0*f[i].potent[GRADZVZ][j]*f[i].potent[GRADZVZ][j];
			        sum+=(f[i].potent[GRADYVX][j]+f[i].potent[GRADXVY][j])*(f[i].potent[GRADYVX][j]+f[i].potent[GRADXVY][j]);
			        sum+=(f[i].potent[GRADZVX][j]+f[i].potent[GRADXVZ][j])*(f[i].potent[GRADZVX][j]+f[i].potent[GRADXVZ][j]);
			        sum+=(f[i].potent[GRADYVZ][j]+f[i].potent[GRADZVY][j])*(f[i].potent[GRADYVZ][j]+f[i].potent[GRADZVY][j]);
					// следующее слагаемое может сделать подкоренное выражение отрицательным.
		            // вычитаем две трети квадрата дивергенции.
					sum-=(2.0/3.0)*(f[i].potent[GRADXVX][j]+f[i].potent[GRADYVY][j]+f[i].potent[GRADZVZ][j])*(f[i].potent[GRADXVX][j]+f[i].potent[GRADYVY][j]+f[i].potent[GRADZVZ][j]); // добавок связанный с несжимаемостью/сжимаемостью
			        f[i].SInvariantStrainRateTensor[j]=sqrt(fmax(0.0,sum));

                    // Вихрь (модуль ротора скорости).
					// помоему это тоже самое что и модуль тензора вращения.
		            sum=0.0;
		            sum+=(f[i].potent[GRADYVZ][j]-f[i].potent[GRADZVY][j])*(f[i].potent[GRADYVZ][j]-f[i].potent[GRADZVY][j]); // проверено !
		            sum+=(f[i].potent[GRADZVX][j]-f[i].potent[GRADXVZ][j])*(f[i].potent[GRADZVX][j]-f[i].potent[GRADXVZ][j]);
		            sum+=(f[i].potent[GRADXVY][j]-f[i].potent[GRADYVX][j])*(f[i].potent[GRADXVY][j]-f[i].potent[GRADYVX][j]);
		            f[i].potent[CURL][j]=sqrt(sum);
		        }

				breadOk=true; // считывание прошло успешно.

			} // конец цикла по жидким зонам
			
		    

		}

		

		fclose(fp_avtosave); // закрытие файла из которого считывались сохранённые полевые величины
	}

} // avtoreadvalue


#endif