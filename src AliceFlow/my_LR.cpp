// begin 20 сентября 2011. (для задач гидродинамики)
// begin 22 октября 2011 года для задачи теплопроводности.
// begin2 21 ноября 2011 года для задачи теплопроводности
// в связи с проблемой плоского бесконечно-тонкого источника тепла.
// begin 3 end 29, 30 декабря распаралеливание метода прогонки по стандарту OpenMP.
// my_LR.c - реализация полилинейного метода.
// 12 апреля 2013 года. В связи с переходом на более подробные сетки до 10млн узлов было
// обнаружено что ILU(2) предобуславливание на много лучше чем полинейный метод Зейделя.
// от Lr1sk алгоритма можно отказываться в пользу ilu(2) (FlowVision рекомендует).

// Прямоугольная HEX геометрия позволяет использовать полилинейный метод решения СЛАУ.
// Об эффективности полилинейного метода пишет С. Патанкар. Именно его он рекомендует использовать.
// В статье Л.Н. Фоминой из Кемеровского государственного университета, Кемерово
// сначала показывается что по отдельности методы неэффективны. Т.е. отдельно использование LR1 даёт
// снижение невязки порядка на 4-5 а потом метод буксует. Отдельное применение BiCGStab приводит к колеблющейся
// по невязкам сходимости на уровне хуже чем LR1. Явный метод Н.И. Булеева сходится к тому же уровню что и LR1
// но показывает в разы более худшую скорость сходимости.
// Предлагается гибридное решение основанное на применении прямого сочетания LR1 и BiCGStab.
// В данном файле реализуется LR1.

// 5 мая 2012 года.
// Конфликт доступа к данным:
// Ошибка возникает при одновременном выполнении следующих условий:
// A) Две или более нитей обращаются к одной и той же ячейке памяти.
// B) По крайней мере, один из этих доступов к памяти является записью.
// C) Нити не синхронизируют свой доступ к данной ячейке памяти.
// При одновременном выполнении всех трёх условий порядок доступа становится неопределённым.
// Вывод несколько нитей могут одновременно читать из одной и той-же ячейки памяти если о записи речи нет.

#ifndef  MY_LR_CPP
#define  MY_LR_CPP 1

//#include <math.h>
//#include <omp.h>

// Для сошедшейся задачи подходит уровень среднеквадратических невязок 1.0e-4
// Источник информации мануал по CFX на русском.
// Если норма Чебышёва то невязка сошедшаяся равна 1.0e-3.
// Источник опять же мануал по CFX.

// Евклидова норма вектора
doublereal NormaV(doublereal *V, integer n){
	doublereal norma;
	doublereal s=0.0;
    doublereal dsize= static_cast<doublereal>(1.0*n);

#pragma omp parallel for shared (V,dsize) /*schedule (guided)*/ reduction (+:s)
	for (integer i = 0; i < n; ++i) {
		s += V[i] * V[i];
	}

	norma=sqrt(s / dsize);
	return norma;
} // NormaV

// создаёт информацию о сеточных линиях для гидродинамической части:
void constr_line(FLOW* &f, integer flow_interior) {

	printf("LR preprocessing FLUID start...\n");

	bool *bmask; // был ли посещён ранее внутренний КО.

	// цикл по всем flow_interior
	integer i,j; // счётчики цикла
	for (i=0; i<flow_interior; ++i) {

#if doubleintprecision == 1
		printf("LR preprocessing FLUID%lld...\n", i);
#else
		printf("LR preprocessing FLUID%d...\n", i);
#endif

		
		NODELR_BASE *rootWE=nullptr;
	    NODELR_BASE *rootSN=nullptr;
	    NODELR_BASE *rootBT=nullptr;

		// Для ускорения поиска.
	    // данная модификация введена 12 апреля 2013 года для ускорения приложения,
	    // так как была обнаружена проблема производительности на сетке из 10млн узлов.
	    NODELR_BASE *rootWE_end=nullptr;
	    NODELR_BASE *rootSN_end=nullptr;
	    NODELR_BASE *rootBT_end=nullptr;

		bmask=new bool[f[i].maxelm];
		for (j=0; j<f[i].maxelm; ++j) bmask[j]=false; // инициализация
		integer iline=0; // номер сеточной линии начиная с нуля.
		// цикл по всем внутренним контрольным объёмам:
		for (j=0; j<f[i].maxelm; ++j) {
			integer iP=j, idubl=j;
			if (!bmask[iP]) {
				// если узел ещё не был посещён
				while (iP < f[i].maxelm) {
					iP=f[i].neighbors_for_the_internal_node[W_SIDE][0][iP];
					if (iP < f[i].maxelm) idubl=iP;
				}
				iP=idubl; // самый W-ый внутренний сеточный узел в сеточной линии WE
				NODELR_BASE *r1=rootWE;
				integer iN=2; // два граничных узла

				if (r1 != nullptr) {
					//while (r1->next !=nullptr) {
						//r1=r1->next;
					//}
					r1=rootWE_end;
					r1->next=new NODELR_BASE;
					r1=r1->next;
					r1->ilineid=iline++;
					r1->next=nullptr;
					r1->root=nullptr;
					rootWE_end=r1;
				}
				else {
					rootWE=new NODELR_BASE;
					rootWE->ilineid=iline++;
					rootWE->next=nullptr; // следующая сеточная линия
					rootWE->root=nullptr; // текущая сеточная линия
					r1=rootWE;
					rootWE_end=r1;
				}
				// r1->root - указывает на необходимый корень.
				r1->root=new NODELR;
				NODELR* r2=r1->root;
				while (iP < f[i].maxelm) {
					bmask[iP]=true; // узел был посещён.
					iN++;
					r2->id=iP;
					iP = f[i].neighbors_for_the_internal_node[E_SIDE][0][iP];
					if (iP < f[i].maxelm) {
						idubl=iP;
						r2->next=new NODELR;
						r2=r2->next;
					}
					else r2->next=nullptr;
				} // while

				r1->iN=iN;

				r1=nullptr; r2=nullptr; // линия построена

			}
		}
		
		rootWE_end=nullptr;
		if (bmask != nullptr) {
			delete[] bmask;
			bmask = nullptr;
		}

		bmask=new bool[f[i].maxelm];
		for (j=0; j<f[i].maxelm; ++j) bmask[j]=false; // инициализация
	    iline=0; // номер сеточной линии начиная с нуля.
		// цикл по всем внутренним контрольным объёмам:
		for (j=0; j<f[i].maxelm; ++j) {
			integer iP=j, idubl=j;
			if (!bmask[iP]) {
				// если узел ещё не был посещён
				while (iP < f[i].maxelm) {
					iP = f[i].neighbors_for_the_internal_node[S_SIDE][0][iP];
					if (iP < f[i].maxelm) idubl=iP;
				}
				iP=idubl; // самый S-ый внутренний сеточный узел в сеточной линии SN
				NODELR_BASE *r1=rootSN;
				integer iN=2; // два граничных узла

				if (r1 != nullptr) {
					//while (r1->next !=nullptr) {
						//r1=r1->next;
					//}
					r1=rootSN_end;
					r1->next=new NODELR_BASE;
					r1=r1->next;
					r1->ilineid=iline++;
					r1->next=nullptr;
					r1->root=nullptr;
					rootSN_end=r1;
				}
				else {
					rootSN=new NODELR_BASE;
					rootSN->ilineid=iline++;
					rootSN->next=nullptr; // следующая сеточная линия
					rootSN->root=nullptr; // текущая сеточная линия
					r1=rootSN;
					rootSN_end=r1;
				}
				// r1->root - указывает на необходимый корень.
				r1->root=new NODELR;
				NODELR* r2=r1->root;
				while (iP < f[i].maxelm) {
					bmask[iP]=true; // узел был посещён.
					iN++;
					r2->id=iP;
					iP = f[i].neighbors_for_the_internal_node[N_SIDE][0][iP];
					if (iP < f[i].maxelm) {
						idubl=iP;
						r2->next=new NODELR;
						r2=r2->next;
					}
					else r2->next=nullptr;
				} // while

                r1->iN=iN;   

				r1=nullptr; r2=nullptr; // линия построена

			}
		}
		rootSN_end=nullptr;
		delete[] bmask;

		bmask=new bool[f[i].maxelm];
		for (j=0; j<f[i].maxelm; ++j) bmask[j]=false; // инициализация
		iline=0; // номер сеточной линии начиная с нуля.
		// цикл по всем внутренним контрольным объёмам:
		for (j=0; j<f[i].maxelm; ++j) {
			integer iP=j, idubl=j;
			if (!bmask[iP]) {
				// если узел ещё не был посещён
				while (iP < f[i].maxelm) {
					iP = f[i].neighbors_for_the_internal_node[B_SIDE][0][iP];
					if (iP < f[i].maxelm) idubl=iP;
				}
				iP=idubl; // самый S-ый внутренний сеточный узел в сеточной линии SN
				NODELR_BASE *r1=rootBT;
				integer iN=2; // два граничных узла

				if (r1 != nullptr) {
					//while (r1->next !=nullptr) {
						//r1=r1->next;
					//}
					r1=rootBT_end;
					r1->next=new NODELR_BASE;
					r1=r1->next;
					r1->ilineid=iline++;
					r1->next=nullptr;
					r1->root=nullptr;
					rootBT_end=r1;
				}
				else {
					rootBT=new NODELR_BASE;
					rootBT->ilineid=iline++;
					rootBT->next=nullptr; // следующая сеточная линия
					rootBT->root=nullptr; // текущая сеточная линия
					r1=rootBT;
					rootBT_end=r1;
				}
				// r1->root - указывает на необходимый корень.
				r1->root=new NODELR;
				NODELR* r2=r1->root;
				while (iP < f[i].maxelm) {
					bmask[iP]=true; // узел был посещён.
					iN++;
					r2->id=iP;
					iP = f[i].neighbors_for_the_internal_node[T_SIDE][0][iP];
					if (iP < f[i].maxelm) {
						idubl=iP;
						r2->next=new NODELR;
						r2=r2->next;
					}
					else r2->next=nullptr;
				} // while

				r1->iN=iN;

				r1=nullptr; r2=nullptr; // линия построена

			}
		}
		rootBT_end=nullptr;
		delete[] bmask;


		integer imax=-1;
		NODELR_BASE *r1;
		r1=rootWE;
		while (r1!=nullptr) {
			if (r1->iN>imax) imax=r1->iN;
			r1=r1->next;
		}
		r1=rootSN;
		while (r1!=nullptr) {
			if (r1->iN>imax) imax=r1->iN;
			r1=r1->next;
		}
		r1=rootBT;
		while (r1!=nullptr) {
			if (r1->iN>imax) imax=r1->iN;
			r1=r1->next;
		}

		NODELR_BASE *rootscan=rootWE;
		integer il=0;
		while (rootscan!=nullptr) {
			il++;
			rootscan=rootscan->next;
		}
		f[i].iWE=il;

		
		rootscan=rootSN;
	    il=0;
		while (rootscan!=nullptr) {
			il++;
			rootscan=rootscan->next;
		}
		f[i].iSN=il;

		
		rootscan=rootBT;
		il=0;
		while (rootscan!=nullptr) {
			il++;
			rootscan=rootscan->next;
		}
		f[i].iBT=il;

		f[i].iN=new integer*[3];
		f[i].id=new integer**[3];
		for (integer i1=0; i1<3; i1++) {
			switch (i1) {
			  case 0: f[i].iN[0]=new integer[f[i].iWE];
				       f[i].id[0]=new integer*[f[i].iWE]; 
					   break;
			  case 1: f[i].iN[1]=new integer[f[i].iSN]; 
				       f[i].id[1]=new integer*[f[i].iSN];
					   break;
			  case 2: f[i].iN[2]=new integer[f[i].iBT];
				       f[i].id[2]=new integer*[f[i].iBT]; 
					   break;
			}
		}

		rootscan=rootWE;
		for (integer ic=0; ic<f[i].iWE; ic++) {
			f[i].iN[0][ic]=rootscan->iN;
			f[i].id[0][ic]=new integer[rootscan->iN-2]; // только внутренние узлы 
			rootscan=rootscan->next;
		}
		rootscan=rootSN;
		for (integer ic=0; ic<f[i].iSN; ic++) {
			f[i].iN[1][ic]=rootscan->iN;
			f[i].id[1][ic]=new integer[rootscan->iN-2]; // только внутренние узлы
			rootscan=rootscan->next;
		}
		rootscan=rootBT;
		for (integer ic=0; ic<f[i].iBT; ic++) {
			f[i].iN[2][ic]=rootscan->iN;
			f[i].id[2][ic]=new integer[rootscan->iN-2]; // только внутренние узлы
			rootscan=rootscan->next;
		}

		
		
		for (integer ic=0; ic<f[i].iWE; ic++) {
			rootscan=rootWE;
			rootWE=rootWE->next;
	
			NODELR* rootlocscan;
			NODELR* rootclean;
			rootlocscan=rootscan->root;
			rootscan->root=nullptr; // обрыв связи
			for (integer i1=0; i1<rootscan->iN-2; i1++) { // т.к. граничные узлы тоже посчитаны а здесь интересуют тока внутренние
				rootclean=rootlocscan;
#if doubleintprecision == 1
				//printf("%lld %lld \n",rootclean->id, rootclean->next->id); // debug
#else
				//printf("%d %d \n",rootclean->id, rootclean->next->id); // debug
#endif

				
				//system("pause");
				f[i].id[0][ic][i1]=rootlocscan->id;
				rootlocscan=rootlocscan->next;
				rootclean->next=nullptr;
				delete rootclean;
			}

			rootscan->next=nullptr;
			delete rootscan; // освобождение оперативной памяти.
		}

		for (integer ic=0; ic<f[i].iSN; ic++) {
			rootscan=rootSN;
			rootSN=rootSN->next;
			
			NODELR* rootlocscan;
			NODELR* rootclean;
			rootlocscan=rootscan->root;
			rootscan->root=nullptr; // обрыв связи
			for (integer i1=0; i1<rootscan->iN-2; i1++) {
				rootclean=rootlocscan;
				f[i].id[1][ic][i1]=rootlocscan->id;
				rootlocscan=rootlocscan->next;
				rootclean->next=nullptr;
				delete rootclean;
			}

			rootscan->next=nullptr;
			delete rootscan; // освобождение оперативной памяти.
		}

		for (integer ic=0; ic<f[i].iBT; ic++) {
			rootscan=rootBT;
			rootBT=rootBT->next;
			
			NODELR* rootlocscan;
			NODELR* rootclean;
			rootlocscan=rootscan->root;
			rootscan->root=nullptr; // обрыв связи
			for (integer i1=0; i1<rootscan->iN-2; i1++) {
				rootclean=rootlocscan;
				f[i].id[2][ic][i1]=rootlocscan->id;
				rootlocscan=rootlocscan->next;
				rootclean->next=nullptr;
				delete rootclean;
			}

			rootscan->next=nullptr;
			delete rootscan; // освобождение оперативной памяти.
		}


	} // конец цикла по всем жидким зонам

} // constr_line

// создаёт информацию о сеточных линиях для тепловой части:
void constr_line_temp(TEMPER &t, BLOCK* b, int lb) {

	printf("LR preprocessing TEMPER start...\n");

	bool *bmask; // был ли посещён ранее внутренний КО.

	
	integer j; // счётчик цикла
	

	bmask=new bool[t.maxelm];
	for (j=0; j<t.maxelm; ++j) bmask[j]=false; // инициализация
	integer iline=0; // номер сеточной линии начиная с нуля.

	NODELR_BASE *rootWE=nullptr;
	NODELR_BASE *rootSN=nullptr;
	NODELR_BASE *rootBT=nullptr;

	// Для ускорения поиска.
	// данная модификация введена 12 апреля 2013 года для ускорения приложения,
	// так как была обнаружена проблема производительности на сетке из 10млн узлов.
	NODELR_BASE *rootWE_end=nullptr;
	NODELR_BASE *rootSN_end=nullptr;
	NODELR_BASE *rootBT_end=nullptr;

	// В каждой сеточной линии перечислены лишь идентификаторы внутренних
	// контрольных объёмов.
	// В алгоритме учитывается что источник тепла является одной из границ сеточной линии.
	// Т.е. каждый источник тепла разделяет линию на две подлинии естественным образом.

	// цикл по всем внутренним контрольным объёмам:
	for (j=0; j<t.maxelm; ++j) {
		integer iP=j, idubl=j;
		if (!bmask[iP]) {
			// если узел ещё не был посещён
			while (iP < t.maxelm) {
				iP = t.neighbors_for_the_internal_node[W_SIDE][0][iP];
				if (iP < t.maxelm) idubl=iP;
				// случай источника тепла также предусмотрен т.к. его идентификатор контрольного объёма меньше t.maxelm.
			}
			iP=idubl; // самый W-ый внутренний сеточный узел в сеточной линии WE
			NODELR_BASE *r1=rootWE;
			integer iN=2; // два граничных узла

			if (r1 != nullptr) {
				// поиск последней сеточной линии.
				//while (r1->next !=nullptr) {
					//r1=r1->next;
				//}
				r1=rootWE_end;
				r1->next=new NODELR_BASE;
				r1=r1->next;
				r1->ilineid=iline++;
				r1->next=nullptr;
				r1->root=nullptr;
				rootWE_end=r1; // запоминаем последнюю сеточную линию.
			}
			else {
				rootWE=new NODELR_BASE;
				rootWE->ilineid=iline++;
				rootWE->next=nullptr; // следующая сеточная линия
				rootWE->root=nullptr; // текущая сеточная линия
				r1=rootWE;
				rootWE_end=rootWE; // указатель на последнюю сеточную линию.
			}
			// r1->root - указывает на необходимый корень.
			r1->root=new NODELR;
			NODELR* r2=r1->root;
			r2->next = nullptr;
			while (iP < t.maxelm) {
				bmask[iP]=true; // узел был посещён.
				iN++;
				r2->id=iP;
				iP = t.neighbors_for_the_internal_node[E_SIDE][0][iP];
				if (iP < t.maxelm) {
					//idubl=iP; // это лишнее TODO
					if (r2->next == nullptr) {
						r2->next = new NODELR;
						r2->next->next = nullptr; // инициализация.
					}
					r2=r2->next;
				}
				else r2->next=nullptr;
			} // while

			r1->iN=iN;

			r1=nullptr; r2=nullptr; // линия построена

		}
	}
	rootWE_end=nullptr; // вспомогательный указатель для ускорения поиска.
		
	if (bmask != nullptr) {
		delete[] bmask;
		bmask = nullptr;
	}

	bmask=new bool[t.maxelm];
	for (j=0; j<t.maxelm; ++j) bmask[j]=false; // инициализация
	iline=0; // номер сеточной линии начиная с нуля.
	// цикл по всем внутренним контрольным объёмам:
	for (j=0; j<t.maxelm; ++j) {
		integer iP=j, idubl=j;
		if (!bmask[iP]) {
			// если узел ещё не был посещён
			while (iP < t.maxelm) {
				iP = t.neighbors_for_the_internal_node[S_SIDE][0][iP];
				if (iP < t.maxelm) idubl=iP;
			}
			iP=idubl; // самый S-ый внутренний сеточный узел в сеточной линии SN
			NODELR_BASE *r1=rootSN;
			integer iN=2; // два граничных узла

			if (r1 != nullptr) {
				//while (r1->next !=nullptr) {
					//r1=r1->next;
				//}
				r1=rootSN_end;
				r1->next=new NODELR_BASE;
				r1=r1->next;
				r1->ilineid=iline++;
				r1->next=nullptr;
				r1->root=nullptr;
				rootSN_end=r1;
			}
			else {
				rootSN=new NODELR_BASE;
				rootSN->ilineid=iline++;
				rootSN->next=nullptr; // следующая сеточная линия
				rootSN->root=nullptr; // текущая сеточная линия
				r1=rootSN;
				rootSN_end=rootSN; // указатель на последнюю сеточную линию.
			}
			// r1->root - указывает на необходимый корень.
			r1->root=new NODELR;
			NODELR* r2=r1->root;
			r2->next = nullptr;
			while (iP < t.maxelm) {
				bmask[iP]=true; // узел был посещён.
				iN++;
				r2->id=iP;
				iP = t.neighbors_for_the_internal_node[N_SIDE][0][iP];
				if (iP < t.maxelm) {
					idubl=iP;
					if (r2->next == nullptr) {
						r2->next = new NODELR;
						r2->next->next = nullptr;
					}
					r2=r2->next;
				}
				else r2->next=nullptr;
			} // while

            r1->iN=iN;   

			r1=nullptr; r2=nullptr; // линия построена

		}
	}

	rootSN_end=nullptr;
	if (bmask != nullptr) {
		delete[] bmask;
		bmask = nullptr;
	}
	

	bmask=new bool[t.maxelm];
	for (j=0; j<t.maxelm; ++j) bmask[j]=false; // инициализация
	iline=0; // номер сеточной линии начиная с нуля.
	// цикл по всем внутренним контрольным объёмам:
	for (j=0; j<t.maxelm; ++j) {
		integer iP=j, idubl=j;
		if (!bmask[iP]) {
			// если узел ещё не был посещён
			while (iP < t.maxelm) {
				iP = t.neighbors_for_the_internal_node[B_SIDE][0][iP];
				if (iP < t.maxelm) idubl=iP;
			}
			iP=idubl; // самый S-ый внутренний сеточный узел в сеточной линии SN
			NODELR_BASE *r1=rootBT;
			integer iN=2; // два граничных узла

			if (r1 != nullptr) {
				//while (r1->next !=nullptr) {
					//r1=r1->next;
				//}
				r1=rootBT_end; // ускоряющее присваивание.
				r1->next=new NODELR_BASE;
				r1=r1->next;
				r1->ilineid=iline++;
				r1->next=nullptr;
				r1->root=nullptr;
				rootBT_end=r1;
			}
			else {
				rootBT=new NODELR_BASE;
				rootBT->ilineid=iline++;
				rootBT->next=nullptr; // следующая сеточная линия
				rootBT->root=nullptr; // текущая сеточная линия
				r1=rootBT;
				rootBT_end=rootBT;
			}
			// r1->root - указывает на необходимый корень.
			r1->root=new NODELR;
			NODELR* r2=r1->root;
			r2->next = nullptr;
			while (iP < t.maxelm) {
				bmask[iP]=true; // узел был посещён.
				iN++;
				r2->id=iP;
				iP = t.neighbors_for_the_internal_node[T_SIDE][0][iP];
				if (iP < t.maxelm) {
					idubl=iP;
					if (r2->next == nullptr) {
						r2->next = new NODELR;
						r2->next->next = nullptr; // инициализация.
					}					
					r2=r2->next;
				}
				else r2->next=nullptr;
			} // while

			r1->iN=iN;

			r1=nullptr; r2=nullptr; // линия построена

		}
	}
	rootBT_end=nullptr;

	delete[] bmask;


	integer imax=-1;
	NODELR_BASE *r1;
	r1=rootWE;
	while (r1!=nullptr) {
		if (r1->iN>imax) imax=r1->iN;
		r1=r1->next;
	}
	r1=rootSN;
	while (r1!=nullptr) {
		if (r1->iN>imax) imax=r1->iN;
		r1=r1->next;
	}
	r1=rootBT;
	while (r1!=nullptr) {
		if (r1->iN>imax) imax=r1->iN;
		r1=r1->next;
	}

	//t.imaxsl=imax;

	/*
	// трехдиагональная матрица
	t.a=new doublereal[imax];
	t.b=new doublereal[imax];
	t.c=new doublereal[imax];
	t.d=new doublereal[imax];
	// прогоночные коэффициенты
	t.P=new doublereal[imax];
	t.Q=new doublereal[imax];
	// глобальный номер перменной (связь)
	t.ind= new integer[imax];
	*/

	// дополнение в связи с особой обработкой
	// плоского источника тепла.

	
    TOCHKA p_c;
	bool bi_fluid;
	int ib; // номер блока

	//ось Оx
	r1=rootWE;
	while (r1 != nullptr) {
		
		// инициализация
		r1->bNeimanStart=false;
		r1->bNeimanEnd=false; 
		// Если будет стоять true то связь с источником с этой стороны будет разорвана.
		// Это нужно для внутренних источников тепла, когда с одной стороны источника жидкость
		// а с другой твёрдое тело.

	    NODELR* r2=r1->root;

		p_c.x=0.5*(t.pa[t.nvtx[0][r2->id]-1].x+t.pa[t.nvtx[1][r2->id]-1].x);
	    p_c.y=0.5*(t.pa[t.nvtx[1][r2->id]-1].y+t.pa[t.nvtx[3][r2->id]-1].y);
	    p_c.z=0.5*(t.pa[t.nvtx[0][r2->id]-1].z+t.pa[t.nvtx[4][r2->id]-1].z);

		// В ib возвращается номер блока которому принадлежит контрольный объём КО.
	    bi_fluid=in_model_flow(p_c,ib,b,lb); // принадлежит ли КО жидкой зоне.

		if ((bi_fluid) && (t.binternalsource[t.neighbors_for_the_internal_node[W_SIDE][0][r2->id] - t.maxelm])) r1->bNeimanStart = true;
		
		while (r2->next!=nullptr) r2=r2->next; // перемотка в конец списка

		p_c.x=0.5*(t.pa[t.nvtx[0][r2->id]-1].x+t.pa[t.nvtx[1][r2->id]-1].x);
	    p_c.y=0.5*(t.pa[t.nvtx[1][r2->id]-1].y+t.pa[t.nvtx[3][r2->id]-1].y);
	    p_c.z=0.5*(t.pa[t.nvtx[0][r2->id]-1].z+t.pa[t.nvtx[4][r2->id]-1].z);

	    bi_fluid=in_model_flow(p_c,ib,b,lb); // принадлежит ли КО жидкой зоне.

		if ((bi_fluid) && (t.binternalsource[t.neighbors_for_the_internal_node[E_SIDE][0][r2->id] - t.maxelm])) r1->bNeimanEnd = true;

		r1=r1->next; // переход к следующей линиии
	}

	//ось Оy
	r1=rootSN;
	while (r1 != nullptr) {
		
		// инициализация
		r1->bNeimanStart=false;
		r1->bNeimanEnd=false; 

	    NODELR* r2=r1->root;

		p_c.x=0.5*(t.pa[t.nvtx[0][r2->id]-1].x+t.pa[t.nvtx[1][r2->id]-1].x);
	    p_c.y=0.5*(t.pa[t.nvtx[1][r2->id]-1].y+t.pa[t.nvtx[3][r2->id]-1].y);
	    p_c.z=0.5*(t.pa[t.nvtx[0][r2->id]-1].z+t.pa[t.nvtx[4][r2->id]-1].z);

		// В ib возвращается номер блока которому принадлежит контрольный объём КО.
	    bi_fluid=in_model_flow(p_c,ib,b,lb); // принадлежит ли КО жидкой зоне.

		if ((bi_fluid) && (t.binternalsource[t.neighbors_for_the_internal_node[S_SIDE][0][r2->id] - t.maxelm])) r1->bNeimanStart = true;
		
		while (r2->next!=nullptr) r2=r2->next; // перемотка в конец списка

		p_c.x=0.5*(t.pa[t.nvtx[0][r2->id]-1].x+t.pa[t.nvtx[1][r2->id]-1].x);
	    p_c.y=0.5*(t.pa[t.nvtx[1][r2->id]-1].y+t.pa[t.nvtx[3][r2->id]-1].y);
	    p_c.z=0.5*(t.pa[t.nvtx[0][r2->id]-1].z+t.pa[t.nvtx[4][r2->id]-1].z);

	    bi_fluid=in_model_flow(p_c,ib,b,lb); // принадлежит ли КО жидкой зоне.

		if ((bi_fluid) && (t.binternalsource[t.neighbors_for_the_internal_node[N_SIDE][0][r2->id] - t.maxelm])) r1->bNeimanEnd = true;

		r1=r1->next; // переход к следующей линиии
	}

	// ось Оz
	r1=rootBT;
	while (r1 != nullptr) {
		
		// инициализация
		r1->bNeimanStart=false;
		r1->bNeimanEnd=false; 

	    NODELR* r2=r1->root;

		p_c.x=0.5*(t.pa[t.nvtx[0][r2->id]-1].x+t.pa[t.nvtx[1][r2->id]-1].x);
	    p_c.y=0.5*(t.pa[t.nvtx[1][r2->id]-1].y+t.pa[t.nvtx[3][r2->id]-1].y);
	    p_c.z=0.5*(t.pa[t.nvtx[0][r2->id]-1].z+t.pa[t.nvtx[4][r2->id]-1].z);

		// В ib возвращается номер блока которому принадлежит контрольный объём КО.
	    bi_fluid=in_model_flow(p_c,ib,b,lb); // принадлежит ли КО жидкой зоне.

		if ((bi_fluid) && (t.binternalsource[t.neighbors_for_the_internal_node[B_SIDE][0][r2->id] - t.maxelm])) r1->bNeimanStart = true;
		
		while (r2->next!=nullptr) r2=r2->next; // перемотка в конец списка

		p_c.x=0.5*(t.pa[t.nvtx[0][r2->id]-1].x+t.pa[t.nvtx[1][r2->id]-1].x);
	    p_c.y=0.5*(t.pa[t.nvtx[1][r2->id]-1].y+t.pa[t.nvtx[3][r2->id]-1].y);
	    p_c.z=0.5*(t.pa[t.nvtx[0][r2->id]-1].z+t.pa[t.nvtx[4][r2->id]-1].z);

	    bi_fluid=in_model_flow(p_c,ib,b,lb); // принадлежит ли КО жидкой зоне.

		if ((bi_fluid) && (t.binternalsource[t.neighbors_for_the_internal_node[T_SIDE][0][r2->id] - t.maxelm])) r1->bNeimanEnd = true;

		r1=r1->next; // переход к следующей линиии
	}

	NODELR_BASE *rootscan=rootWE;

	integer il=0;
	while (rootscan!=nullptr) {
		il++;
		rootscan=rootscan->next;
	}
	t.iWE=il;

	il=0;
	rootscan=rootSN;
	while (rootscan!=nullptr) {
		il++;
		rootscan=rootscan->next;
	}
	t.iSN=il;

	il=0;
	rootscan=rootBT;
	while (rootscan!=nullptr) {
		il++;
		rootscan=rootscan->next;
	}
	t.iBT=il;

	// Копирование построенных структур данных в рабочую структуру.
	// Список по связи next преобразовался в одномерный массив.
	// Рабочая структура это одномерный массив линейных списков.
	if (t.rootWE != nullptr) {
		delete t.rootWE;
		t.rootWE = nullptr;
	}
	t.rootWE=new NODELR_BASE[t.iWE];
	for (integer i1=0; i1<t.iWE; i1++) {
		rootscan=rootWE;
		t.rootWE[i1].bNeimanStart=rootscan->bNeimanStart;
		t.rootWE[i1].bNeimanEnd=rootscan->bNeimanEnd;
		t.rootWE[i1].ilineid=rootscan->ilineid;
		t.rootWE[i1].iN=rootscan->iN;
		t.rootWE[i1].root=rootscan->root;
		t.rootWE[i1].next=nullptr;
		rootWE=rootWE->next;
		rootscan->next=nullptr;
		rootscan->root=nullptr;
		delete rootscan;
	}

	if (t.rootSN != nullptr) {
		delete t.rootSN;
		t.rootSN = nullptr;
	}
	t.rootSN=new NODELR_BASE[t.iSN];
	for (integer i1=0; i1<t.iSN; i1++) {
		rootscan=rootSN;
		t.rootSN[i1].bNeimanStart=rootscan->bNeimanStart;
		t.rootSN[i1].bNeimanEnd=rootscan->bNeimanEnd;
		t.rootSN[i1].ilineid=rootscan->ilineid;
		t.rootSN[i1].iN=rootscan->iN;
		t.rootSN[i1].root=rootscan->root;
		t.rootSN[i1].next=nullptr;
		rootSN=rootSN->next;
		rootscan->next=nullptr;
		rootscan->root=nullptr;
		delete rootscan;
	}

	if (t.rootBT != nullptr) {
		delete t.rootBT;
		t.rootBT = nullptr;
	}
	t.rootBT=new NODELR_BASE[t.iBT];
	for (integer i1=0; i1<t.iBT; i1++) {
		rootscan=rootBT;
		t.rootBT[i1].bNeimanStart=rootscan->bNeimanStart;
		t.rootBT[i1].bNeimanEnd=rootscan->bNeimanEnd;
		t.rootBT[i1].ilineid=rootscan->ilineid;
		t.rootBT[i1].iN=rootscan->iN;
		t.rootBT[i1].root=rootscan->root;
		t.rootBT[i1].next=nullptr;
		rootBT=rootBT->next;
		rootscan->next=nullptr;
		rootscan->root=nullptr;
		delete rootscan;
	}

} // constr_line_temp

// Освобождает память
void free_root(NODELR_BASE* &root, integer &isize) {
	if (root!=nullptr) {
		NODELR  *r1, *r2;
		for (integer i=0; i<isize; ++i) {
			r1=root[i].root;
			root[i].root=nullptr;
			// уничтожение нити
			while (r1->next!=nullptr) {
				r2=r1;
				r1=r1->next;
				r2->next=nullptr;
				delete r2;
				r2=nullptr;
			}
			delete r1;
			r1=nullptr;
		}
		delete root;
		root=nullptr;
		isize=0;
	}
} 

// полилинейный солвер LR1:
// Л.Н. Фомина Кемеровский государственный университет, Кемерово
// О повышении эффективности полилинейного рекурретного метода решения 
// разностных эллиптических уравнений.
void solveLR1(doublereal* &x, doublereal* &rthdsd, integer ns, integer iVar, bool bnorelax,
	int*** neighbors_for_the_internal_node, integer maxelm, equation3D** slau, equation3D_bon** slau_bon,
			   integer** iN, integer***id, integer iWE, integer iSN, integer iBT, doublereal *alpha, integer maxbound) {

	// Если bnorelax==true то нижняя релаксация не применяется.

    doublereal alphaP=1.0;
	if (!bnorelax) {
		switch (iVar) {
		    case VELOCITY_X_COMPONENT: alphaP=alpha[iVar]; break;
		    case VELOCITY_Y_COMPONENT: alphaP=alpha[iVar]; break;
		    case VELOCITY_Z_COMPONENT: alphaP=alpha[iVar]; break;
	    }
	}
	else alphaP=1.0;
	doublereal alphaPb=alphaP; // 1.0; norelax boundary! // релаксация граничных узлов

	// ns - размерность вектора x и правой части rthdsd.
	// f.slau - хранитель матрцы СЛАУ.
	// x - вектор решения.
	// rthdsd - вектор правой части.

	//f.a, f.b, f.c, f.d - трёхдиагональная матрица СЛАУ и правая часть
	// f.P, f.Q -  прогоночные коэффициенты
	// f.ind - глобальный номер переменной

    // Если стоит условие Дирихле то нижней релаксации нет, а если стоит условие Неймана то нижняя релаксация есть.

	// WE direct

	// проход по всем линиям:

	// гонка данных: одна и таже переменная читается и записывается одновременно.
	// Во избежании гонки данных введём новый массив:
	doublereal* xc=new doublereal[maxelm+maxbound];
	#pragma omp parallel for shared (x, xc) schedule (guided) // 
	for (integer j=0; j<maxelm+maxbound; ++j) {
		xc[j]=x[j]; // xc - x copy
	}

    #pragma omp parallel for shared (neighbors_for_the_internal_node, maxelm, slau_bon, iVar, slau, rthdsd, alphaP, alphaPb, x, iN, id) schedule (guided) // 
	for (integer j=0; j<iWE; ++j) {
		/* глобальные переменные:
		*  f.neighbors_for_the_internal_node, f.maxelm, f.slau_bon, iVar,
		*  f.slau, rthdsd, alphaP, alphaPb, x, 
		*  f.rootWE
		*/
		
		integer n=iN[0][j];

		// Вставка данных для паралельной обработки:
		doublereal *a, *b, *c, *d; // трёхдиагональная матрица и правая часть
	    doublereal *P, *Q; // прогоночные коэффициенты
	    integer *ind; // связь с глобальной нумерацией узлов.
		// трехдиагональная матрица
	    a=new doublereal[n];
	    b=new doublereal[n];
	    c=new doublereal[n];
	    d=new doublereal[n];
	    // прогоночные коэффициенты
	    P=new doublereal[n];
	    Q=new doublereal[n];
	    // глобальный номер перменной (связь)
	    ind= new integer[n];

		
		// т.к. граничные условия:
		c[0]=0.0;
		b[n-1]=0.0;

		integer iP=neighbors_for_the_internal_node[W_SIDE][0][id[0][j][0]]; // номер граничного узла
		a[0]=slau_bon[iVar][iP-maxelm].aw;///alphaP; // граничные узлы не релаксируются если на них стоит условие Дирихле
		if (slau_bon[iVar][iP-maxelm].iI>-1) {
			b[0]=slau_bon[iVar][iP-maxelm].ai;
			a[0]=a[0]/alphaPb;
		} else b[0]=0.0;
		d[0]=rthdsd[iP]; //slau_bon[iVar][iP-maxelm].b;
		ind[0]=iP;

		for (integer i=1; i<n-1; ++i) {
			iP=id[0][j][i-1];
			ind[i]=iP;
			a[i]=slau[iVar][iP].ap/alphaP;
			b[i]=slau[iVar][iP].ae;
			c[i]=slau[iVar][iP].aw;
			d[i] = slau[iVar][iP].an*xc[neighbors_for_the_internal_node[N_SIDE][0][iP]];
			d[i] += slau[iVar][iP].as*xc[neighbors_for_the_internal_node[S_SIDE][0][iP]];
			d[i] += slau[iVar][iP].at*xc[neighbors_for_the_internal_node[T_SIDE][0][iP]];
			d[i] += slau[iVar][iP].ab*xc[neighbors_for_the_internal_node[B_SIDE][0][iP]];
			d[i]+=rthdsd[iP]; //f.slau[iVar][iP].b;
		}

		iP = neighbors_for_the_internal_node[E_SIDE][0][id[0][j][n - 3]]; // номер граничного узла
		a[n-1]=slau_bon[iVar][iP-maxelm].aw;///alphaP; // граничные узлы не релаксируются
        if (slau_bon[iVar][iP-maxelm].iI>-1) {
		    c[n-1]=slau_bon[iVar][iP-maxelm].ai;
			a[n-1]=a[n-1]/alphaPb;
		} else c[n-1]=0.0;
		d[n-1]=rthdsd[iP]; //f.slau_bon[iVar][iP-f.maxelm].b;
        ind[n-1]=iP;

		// собственно сама прогонка или алгоритм Томаса:
        P[0]=b[0]/a[0];
		Q[0]=d[0]/a[0];
		doublereal rdiv;
		for (integer i=1; i<n; ++i) {
            rdiv=1.0/(a[i]-c[i]*P[i-1]);
			P[i]=b[i]*rdiv;
			Q[i]=(d[i]+c[i]*Q[i-1])*rdiv;
		}
		x[ind[n-1]]=Q[n-1];
		for (integer i=n-2; i>=0; i--) {
			x[ind[i]]=P[i]*x[ind[i+1]]+Q[i];
		}       

		// Освобождение оперативной памяти
		delete[] a;
	    delete[] b;
	    delete[] c;
	    delete[] d;
	    // прогоночные коэффициенты
	    delete[] P;
	    delete[] Q;
	    // глобальный номер перменной (связь)
	    delete[] ind;

	} // for

//#pragma omp barrier


	// SN direct
	
	// проход по всем линиям:
	#pragma omp parallel for shared (x, xc) schedule (guided) // 
	for (integer j=0; j<maxelm+maxbound; ++j) {
		xc[j]=x[j]; // xc - x copy
	}

    #pragma omp parallel for shared (neighbors_for_the_internal_node, maxelm, slau_bon, iVar, slau, rthdsd, alphaP, alphaPb, x, iN, id) schedule (guided) // 
	for (integer j=0; j<iSN; ++j) {

		// printf("Ok\n"); // debug

		integer n=iN[1][j];

		// Вставка данных для паралельной обработки:
		doublereal *a, *b, *c, *d; // трёхдиагональная матрица и правая часть
	    doublereal *P, *Q; // прогоночные коэффициенты
	    integer *ind; // связь с глобальной нумерацией узлов.
		// трехдиагональная матрица
	    a=new doublereal[n];
	    b=new doublereal[n];
	    c=new doublereal[n];
	    d=new doublereal[n];
	    // прогоночные коэффициенты
	    P=new doublereal[n];
	    Q=new doublereal[n];
	    // глобальный номер перменной (связь)
	    ind= new integer[n];
		
		// т.к. граничные условия:
		c[0]=0.0;
		b[n-1]=0.0;

		integer iP = neighbors_for_the_internal_node[S_SIDE][0][id[1][j][0]]; // номер граничного узла
		a[0]=slau_bon[iVar][iP-maxelm].aw;///alphaP;
		if (slau_bon[iVar][iP-maxelm].iI>-1) {
			b[0]=slau_bon[iVar][iP-maxelm].ai;
			a[0]=a[0]/alphaPb;
		} else b[0]=0.0;
		d[0]=rthdsd[iP]; //slau_bon[iVar][iP-maxelm].b;
		ind[0]=iP;

		for (integer i=1; i<n-1; ++i) {
			iP=id[1][j][i-1];
			ind[i]=iP;
			a[i]=slau[iVar][iP].ap/alphaP;
			b[i]=slau[iVar][iP].an;
			c[i]=slau[iVar][iP].as;
			d[i] = slau[iVar][iP].ae*xc[neighbors_for_the_internal_node[E_SIDE][0][iP]];
			d[i] += slau[iVar][iP].aw*xc[neighbors_for_the_internal_node[W_SIDE][0][iP]];
			d[i] += slau[iVar][iP].at*xc[neighbors_for_the_internal_node[T_SIDE][0][iP]];
			d[i] += slau[iVar][iP].ab*xc[neighbors_for_the_internal_node[B_SIDE][0][iP]];
			d[i]+=rthdsd[iP]; //f.slau[iVar][iP].b
		}


		iP = neighbors_for_the_internal_node[N_SIDE][0][id[1][j][n - 3]]; // номер граничного узла
		a[n-1]=slau_bon[iVar][iP-maxelm].aw;///alphaP;
        if (slau_bon[iVar][iP-maxelm].iI>-1) {
		    c[n-1]=slau_bon[iVar][iP-maxelm].ai;
			a[n-1]=a[n-1]/alphaPb;
		} else c[n-1]=0.0;
		d[n-1]=rthdsd[iP]; //f.slau_bon[iVar][iP-f.maxelm].b;
        ind[n-1]=iP;
		

		// собственно сама прогонка или алгоритм Томаса:
        P[0]=b[0]/a[0];
		Q[0]=d[0]/a[0];
		doublereal rdiv;
		for (integer i=1; i<n; ++i) {
            rdiv=1.0/(a[i]-c[i]*P[i-1]);
			P[i]=b[i]*rdiv;
			Q[i]=(d[i]+c[i]*Q[i-1])*rdiv;
		}
		// printf("Ok2\n"); // debug
		x[ind[n-1]]=Q[n-1];
		for (integer i=n-2; i>=0; i--) {
			x[ind[i]]=P[i]*x[ind[i+1]]+Q[i];
		}       

		// Освобождение оперативной памяти
		delete[] a;
	    delete[] b;
	    delete[] c;
	    delete[] d;
	    // прогоночные коэффициенты
	    delete[] P;
	    delete[] Q;
	    // глобальный номер перменной (связь)
	    delete[] ind;

	} // for
	


	// BT direct

	// проход по всем линиям:
	#pragma omp parallel for shared (x, xc) schedule (guided) // 
	for (integer j=0; j<maxelm+maxbound; ++j) {
		xc[j]=x[j]; // xc - x copy
	}

    #pragma omp parallel for shared (neighbors_for_the_internal_node, maxelm, slau_bon, iVar, slau, rthdsd, alphaP, alphaPb, x, iN, id) schedule (guided) // 
	for (integer j=0; j<iBT; ++j) {

		integer n=iN[2][j];

		// Вставка данных для паралельной обработки:
		doublereal *a=nullptr, *b=nullptr, *c=nullptr, *d=nullptr; // трёхдиагональная матрица и правая часть
	    doublereal *P=nullptr, *Q=nullptr; // прогоночные коэффициенты
	    integer *ind=nullptr; // связь с глобальной нумерацией узлов.
		// трехдиагональная матрица
	    a=new doublereal[n];
	    b=new doublereal[n];
	    c=new doublereal[n];
	    d=new doublereal[n];
	    // прогоночные коэффициенты
	    P=new doublereal[n];
	    Q=new doublereal[n];
	    // глобальный номер перменной (связь)
	    ind= new integer[n];
		
		// т.к. граничные условия:
		c[0]=0.0;
		b[n-1]=0.0;
 
		integer iP = neighbors_for_the_internal_node[B_SIDE][0][id[2][j][0]]; // номер граничного узла
		a[0]=slau_bon[iVar][iP-maxelm].aw;///alphaP;
		if (slau_bon[iVar][iP-maxelm].iI>-1) {
			b[0]=slau_bon[iVar][iP-maxelm].ai;
			a[0]=a[0]/alphaPb;
		} else b[0]=0.0;
		d[0]=rthdsd[iP]; //slau_bon[iVar][iP-maxelm].b;
		ind[0]=iP;

		for (integer i=1; i<n-1; ++i) {
			iP=id[2][j][i-1];
			ind[i]=iP;
			a[i]=slau[iVar][iP].ap/alphaP;
			b[i]=slau[iVar][iP].at;
			c[i]=slau[iVar][iP].ab;
			d[i] = slau[iVar][iP].ae*xc[neighbors_for_the_internal_node[E_SIDE][0][iP]];
			d[i] += slau[iVar][iP].aw*xc[neighbors_for_the_internal_node[W_SIDE][0][iP]];
			d[i] += slau[iVar][iP].an*xc[neighbors_for_the_internal_node[N_SIDE][0][iP]];
			d[i] += slau[iVar][iP].as*xc[neighbors_for_the_internal_node[S_SIDE][0][iP]];
			d[i]+=rthdsd[iP]; //f.slau[iVar][iP].b;
		}
		
		iP = neighbors_for_the_internal_node[T_SIDE][0][id[2][j][n - 3]]; // номер граничного узла
		a[n-1]=slau_bon[iVar][iP-maxelm].aw;//alphaP;
        if (slau_bon[iVar][iP-maxelm].iI>-1) {
		    c[n-1]=slau_bon[iVar][iP-maxelm].ai;
			a[n-1]=a[n-1]/alphaPb;
		} else c[n-1]=0.0;
		d[n-1]=rthdsd[iP]; //slau_bon[iVar][iP-maxelm].b;
        ind[n-1]=iP;

		// собственно сама прогонка или алгоритм Томаса:
        P[0]=b[0]/a[0];
		Q[0]=d[0]/a[0];
		doublereal rdiv;
		for (integer i=1; i<n; ++i) {
            rdiv=1.0/(a[i]-c[i]*P[i-1]);
			P[i]=b[i]*rdiv;
			Q[i]=(d[i]+c[i]*Q[i-1])*rdiv;
		}
		x[ind[n-1]]=Q[n-1];
		for (integer i=n-2; i>=0; i--) {
			x[ind[i]]=P[i]*x[ind[i+1]]+Q[i];
		}       

		// Освобождение оперативной памяти
		//if (a != nullptr) {
			// оператор delete может быть применен повторно в том числе и к нулевому указателю. Проверки на null излишни.
			delete[] a;
			a = nullptr;
		//}
		//if (b != nullptr) {
		    // оператор delete может быть применен повторно в том числе и к нулевому указателю. Проверки на null излишни.
			delete[] b;
			b = nullptr;
		//}
		//if (c != nullptr) {
			// оператор delete может быть применен повторно в том числе и к нулевому указателю. Проверки на null излишни.
			delete[] c;
			c = nullptr;
		//}
		//if (d != nullptr) {
			// оператор delete может быть применен повторно в том числе и к нулевому указателю. Проверки на null излишни.
			delete[] d;
			d = nullptr;
		//}
	    // прогоночные коэффициенты
		//if (P != nullptr) {
		// оператор delete может быть применен повторно в том числе и к нулевому указателю. Проверки на null излишни.
			delete[] P;
			P = nullptr;
		//}
		//if (Q != nullptr) {
		// оператор delete может быть применен повторно в том числе и к нулевому указателю. Проверки на null излишни.
			delete[] Q;
			Q = nullptr;
		//}
	    // глобальный номер перменной (связь)
		//if (ind != nullptr) {
		// оператор delete может быть применен повторно в том числе и к нулевому указателю. Проверки на null излишни.
			delete[] ind;
			ind = nullptr;
		//}

	} // for
 
	if (xc != nullptr) {
		delete[] xc; // Освобождение оперативной памяти
	}
    
} // solveLR1

// Однопоточная версия.
// полилинейный солвер LR1:
// Л.Н. Фомина Кемеровский государственный университет, Кемерово
// О повышении эффективности полилинейного рекурретного метода решения 
// разностных эллиптических уравнений.
void solveLR1_serial_sor(doublereal* &x, doublereal* &rthdsd, integer ns, integer iVar, bool bnorelax,
	int*** neighbors_for_the_internal_node, integer maxelm, equation3D** slau, equation3D_bon** slau_bon,
	integer** iN, integer***id, integer iWE, integer iSN, integer iBT, doublereal *alpha, integer maxbound, doublereal omega) {

	// Если bnorelax==true то нижняя релаксация не применяется.
	//doublereal omega = 1.855;//1.855 Нив коем случае не применять для компонент скорости.

	doublereal alphaP = 1.0;
	if (!bnorelax) {
		switch (iVar) {
		case VELOCITY_X_COMPONENT: alphaP = alpha[iVar]; break;
		case VELOCITY_Y_COMPONENT: alphaP = alpha[iVar]; break;
		case VELOCITY_Z_COMPONENT: alphaP = alpha[iVar]; break;
		}
	}
	else alphaP = 1.0;
	doublereal alphaPb = alphaP; // 1.0; norelax boundary! // релаксация граничных узлов

	// ns - размерность вектора x и правой части rthdsd.
	// f.slau - хранитель матрцы СЛАУ.
	// x - вектор решения.
	// rthdsd - вектор правой части.

	//f.a, f.b, f.c, f.d - трёхдиагональная матрица СЛАУ и правая часть
	// f.P, f.Q -  прогоночные коэффициенты
	// f.ind - глобальный номер переменной

	// Если стоит условие Дирихле то нижней релаксации нет, а если стоит условие Неймана то нижняя релаксация есть.

	// WE direct

	// проход по всем линиям:

	

	// гонка данных: одна и таже переменная читается и записывается одновременно.
	// Во избежании гонки данных введём новый массив:
	doublereal* xc = new doublereal[maxelm + maxbound];
//#pragma omp parallel for shared (x, xc) schedule (guided) // 
	for (integer j = 0; j<maxelm + maxbound; ++j) {
		xc[j] = x[j]; // xc - x copy
	}

//#pragma omp parallel for shared (neighbors_for_the_internal_node, maxelm, slau_bon, iVar, slau, rthdsd, alphaP, alphaPb, x, iN, id) schedule (guided) // 
	for (integer j = 0; j<iWE; ++j) {
		

		/* глобальные переменные:
		*  f.neighbors_for_the_internal_node, f.maxelm, f.slau_bon, iVar,
		*  f.slau, rthdsd, alphaP, alphaPb, x,
		*  f.rootWE
		*/

		integer n = iN[0][j];

		// Вставка данных для паралельной обработки:
		doublereal *a, *b, *c, *d; // трёхдиагональная матрица и правая часть
		doublereal *P, *Q; // прогоночные коэффициенты
		integer *ind; // связь с глобальной нумерацией узлов.
		// трехдиагональная матрица
		a = new doublereal[n];
		b = new doublereal[n];
		c = new doublereal[n];
		d = new doublereal[n];
		// прогоночные коэффициенты
		P = new doublereal[n];
		Q = new doublereal[n];
		// глобальный номер перменной (связь)
		ind = new integer[n];


		// т.к. граничные условия:
		c[0] = 0.0;
		b[n - 1] = 0.0;

		

		integer iP = neighbors_for_the_internal_node[W_SIDE][0][id[0][j][0]]; // номер граничного узла
		a[0] = slau_bon[iVar][iP - maxelm].aw;///alphaP; // граничные узлы не релаксируются если на них стоит условие Дирихле
		if (slau_bon[iVar][iP - maxelm].iI>-1) {
			b[0] = slau_bon[iVar][iP - maxelm].ai;
			a[0] = a[0] / alphaPb;
		}
		else b[0] = 0.0;
		d[0] = rthdsd[iP]; //slau_bon[iVar][iP-maxelm].b;
		ind[0] = iP;

		

		for (integer i = 1; i<n - 1; ++i) {
			// n-2-1==n-3
			iP = id[0][j][i - 1];
			ind[i] = iP;
			a[i] = slau[iVar][iP].ap / alphaP;
			b[i] = slau[iVar][iP].ae;
			c[i] = slau[iVar][iP].aw;
			//d[i] = slau[iVar][iP].an*xc[neighbors_for_the_internal_node[NSIDE][0][iP]];
			//d[i] += slau[iVar][iP].as*xc[neighbors_for_the_internal_node[SSIDE][0][iP]];
			//d[i] += slau[iVar][iP].at*xc[neighbors_for_the_internal_node[TSIDE][0][iP]];
			//d[i] += slau[iVar][iP].ab*xc[neighbors_for_the_internal_node[BSIDE][0][iP]];
			// serial
			
			
			d[i] = slau[iVar][iP].an*x[neighbors_for_the_internal_node[N_SIDE][0][iP]];
			d[i] += slau[iVar][iP].as*x[neighbors_for_the_internal_node[S_SIDE][0][iP]];
			d[i] += slau[iVar][iP].at*x[neighbors_for_the_internal_node[T_SIDE][0][iP]];
			d[i] += slau[iVar][iP].ab*x[neighbors_for_the_internal_node[B_SIDE][0][iP]];
			d[i] += rthdsd[iP]; //f.slau[iVar][iP].b;

			
		}

		

		iP = neighbors_for_the_internal_node[E_SIDE][0][id[0][j][n - 3]]; // номер граничного узла
		
		
		a[n - 1] = slau_bon[iVar][iP - maxelm].aw;///alphaP; // граничные узлы не релаксируются
		if (slau_bon[iVar][iP - maxelm].iI>-1) {
			c[n - 1] = slau_bon[iVar][iP - maxelm].ai;
			a[n - 1] = a[n - 1] / alphaPb;
		}
		else c[n - 1] = 0.0;
			
		d[n - 1] = rthdsd[iP]; //f.slau_bon[iVar][iP-f.maxelm].b;
		ind[n - 1] = iP;

		

		// собственно сама прогонка или алгоритм Томаса:
		P[0] = b[0] / a[0];
		Q[0] = d[0] / a[0];
		doublereal rdiv;
		for (integer i = 1; i<n; ++i) {
			rdiv = 1.0 / (a[i] - c[i] * P[i - 1]);
			P[i] = b[i] * rdiv;
			Q[i] = (d[i] + c[i] * Q[i - 1])*rdiv;
		}

		
		//x[ind[n - 1]] = Q[n - 1];
		//for (integer i = n - 2; i >= 0; i--) {
			//x[ind[i]] = P[i] * x[ind[i + 1]] + Q[i];
		//}
		xc[ind[n - 1]] = Q[n - 1];
		for (integer i = n - 2; i >= 0; i--) {
			xc[ind[i]] = P[i] * xc[ind[i + 1]] + Q[i];
		}
		for (integer i = n - 2; i >= 1; i--) {
			// завышае новые значения.
			x[ind[i]] = omega*xc[ind[i]] + (1.0 - omega)*x[ind[i]];
		}
		// Граничные узлы не подвержены процессу релаксации.
		x[ind[0]] = xc[ind[0]];
		x[ind[n - 1]] = xc[ind[n - 1]];

		// Освобождение оперативной памяти
		delete[] a;
		delete[] b;
		delete[] c;
		delete[] d;
		// прогоночные коэффициенты
		delete[] P;
		delete[] Q;
		// глобальный номер перменной (связь)
		delete[] ind;

	} // for

	

	//#pragma omp barrier

	
	// SN direct

	// проход по всем линиям:
//#pragma omp parallel for shared (x, xc) schedule (guided) // 
	for (integer j = 0; j<maxelm + maxbound; ++j) {
		xc[j] = x[j]; // xc - x copy
	}

//#pragma omp parallel for shared (neighbors_for_the_internal_node, maxelm, slau_bon, iVar, slau, rthdsd, alphaP, alphaPb, x, iN, id) schedule (guided) // 
	for (integer j = 0; j<iSN; ++j) {

		// printf("Ok\n"); // debug

		integer n = iN[1][j];

		// Вставка данных для паралельной обработки:
		doublereal *a, *b, *c, *d; // трёхдиагональная матрица и правая часть
		doublereal *P, *Q; // прогоночные коэффициенты
		integer *ind; // связь с глобальной нумерацией узлов.
		// трехдиагональная матрица
		a = new doublereal[n];
		b = new doublereal[n];
		c = new doublereal[n];
		d = new doublereal[n];
		// прогоночные коэффициенты
		P = new doublereal[n];
		Q = new doublereal[n];
		// глобальный номер перменной (связь)
		ind = new integer[n];

		// т.к. граничные условия:
		c[0] = 0.0;
		b[n - 1] = 0.0;

		integer iP = neighbors_for_the_internal_node[S_SIDE][0][id[1][j][0]]; // номер граничного узла
		a[0] = slau_bon[iVar][iP - maxelm].aw;///alphaP;
		if (slau_bon[iVar][iP - maxelm].iI>-1) {
			b[0] = slau_bon[iVar][iP - maxelm].ai;
			a[0] = a[0] / alphaPb;
		}
		else b[0] = 0.0;
		d[0] = rthdsd[iP]; //slau_bon[iVar][iP-maxelm].b;
		ind[0] = iP;

		for (integer i = 1; i<n - 1; ++i) {
			iP = id[1][j][i - 1];
			ind[i] = iP;
			a[i] = slau[iVar][iP].ap / alphaP;
			b[i] = slau[iVar][iP].an;
			c[i] = slau[iVar][iP].as;
			//d[i] = slau[iVar][iP].ae*xc[neighbors_for_the_internal_node[ESIDE][0][iP]];
			//d[i] += slau[iVar][iP].aw*xc[neighbors_for_the_internal_node[WSIDE][0][iP]];
			//d[i] += slau[iVar][iP].at*xc[neighbors_for_the_internal_node[TSIDE][0][iP]];
			//d[i] += slau[iVar][iP].ab*xc[neighbors_for_the_internal_node[BSIDE][0][iP]];
			d[i] = slau[iVar][iP].ae*x[neighbors_for_the_internal_node[E_SIDE][0][iP]];
			d[i] += slau[iVar][iP].aw*x[neighbors_for_the_internal_node[W_SIDE][0][iP]];
			d[i] += slau[iVar][iP].at*x[neighbors_for_the_internal_node[T_SIDE][0][iP]];
			d[i] += slau[iVar][iP].ab*x[neighbors_for_the_internal_node[B_SIDE][0][iP]];
			d[i] += rthdsd[iP]; //f.slau[iVar][iP].b
		}


		iP = neighbors_for_the_internal_node[N_SIDE][0][id[1][j][n - 3]]; // номер граничного узла
		a[n - 1] = slau_bon[iVar][iP - maxelm].aw;///alphaP;
		if (slau_bon[iVar][iP - maxelm].iI>-1) {
			c[n - 1] = slau_bon[iVar][iP - maxelm].ai;
			a[n - 1] = a[n - 1] / alphaPb;
		}
		else c[n - 1] = 0.0;
		d[n - 1] = rthdsd[iP]; //f.slau_bon[iVar][iP-f.maxelm].b;
		ind[n - 1] = iP;


		// собственно сама прогонка или алгоритм Томаса:
		P[0] = b[0] / a[0];
		Q[0] = d[0] / a[0];
		doublereal rdiv;
		for (integer i = 1; i<n; ++i) {
			rdiv = 1.0 / (a[i] - c[i] * P[i - 1]);
			P[i] = b[i] * rdiv;
			Q[i] = (d[i] + c[i] * Q[i - 1])*rdiv;
		}
		// printf("Ok2\n"); // debug
		//x[ind[n - 1]] = Q[n - 1];
		//for (integer i = n - 2; i >= 0; i--) {
			//x[ind[i]] = P[i] * x[ind[i + 1]] + Q[i];
		//}
		xc[ind[n - 1]] = Q[n - 1];
		for (integer i = n - 2; i >= 0; i--) {
			xc[ind[i]] = P[i] * xc[ind[i + 1]] + Q[i];
		}
		for (integer i = n - 2; i >= 1; i--) {
			// завышае новые значения.
			x[ind[i]] = omega*xc[ind[i]] + (1.0 - omega)*x[ind[i]];
		}
		// Граничные узлы не подвержены процессу релаксации.
		x[ind[0]] = xc[ind[0]];
		x[ind[n - 1]] = xc[ind[n - 1]];

		// Освобождение оперативной памяти
		delete[] a;
		delete[] b;
		delete[] c;
		delete[] d;
		// прогоночные коэффициенты
		delete[] P;
		delete[] Q;
		// глобальный номер перменной (связь)
		delete[] ind;

	} // for



	// BT direct

	// проход по всем линиям:
//#pragma omp parallel for shared (x, xc) schedule (guided) // 
	for (integer j = 0; j<maxelm + maxbound; ++j) {
		xc[j] = x[j]; // xc - x copy
	}

//#pragma omp parallel for shared (neighbors_for_the_internal_node, maxelm, slau_bon, iVar, slau, rthdsd, alphaP, alphaPb, x, iN, id) schedule (guided) // 
	for (integer j = 0; j<iBT; ++j) {

		integer n = iN[2][j];

		// Вставка данных для паралельной обработки:
		doublereal *a = nullptr, *b = nullptr, *c = nullptr, *d = nullptr; // трёхдиагональная матрица и правая часть
		doublereal *P = nullptr, *Q = nullptr; // прогоночные коэффициенты
		integer *ind = nullptr; // связь с глобальной нумерацией узлов.
		// трехдиагональная матрица
		a = new doublereal[n];
		b = new doublereal[n];
		c = new doublereal[n];
		d = new doublereal[n];
		// прогоночные коэффициенты
		P = new doublereal[n];
		Q = new doublereal[n];
		// глобальный номер перменной (связь)
		ind = new integer[n];

		// т.к. граничные условия:
		c[0] = 0.0;
		b[n - 1] = 0.0;

		integer iP = neighbors_for_the_internal_node[B_SIDE][0][id[2][j][0]]; // номер граничного узла
		a[0] = slau_bon[iVar][iP - maxelm].aw;///alphaP;
		if (slau_bon[iVar][iP - maxelm].iI>-1) {
			b[0] = slau_bon[iVar][iP - maxelm].ai;
			a[0] = a[0] / alphaPb;
		}
		else b[0] = 0.0;
		d[0] = rthdsd[iP]; //slau_bon[iVar][iP-maxelm].b;
		ind[0] = iP;

		for (integer i = 1; i<n - 1; ++i) {
			iP = id[2][j][i - 1];
			ind[i] = iP;
			a[i] = slau[iVar][iP].ap / alphaP;
			b[i] = slau[iVar][iP].at;
			c[i] = slau[iVar][iP].ab;
			//d[i] = slau[iVar][iP].ae*xc[neighbors_for_the_internal_node[ESIDE][0][iP]];
			//d[i] += slau[iVar][iP].aw*xc[neighbors_for_the_internal_node[WSIDE][0][iP]];
			//d[i] += slau[iVar][iP].an*xc[neighbors_for_the_internal_node[NSIDE][0][iP]];
			//d[i] += slau[iVar][iP].as*xc[neighbors_for_the_internal_node[SSIDE][0][iP]];

			d[i] = slau[iVar][iP].ae*x[neighbors_for_the_internal_node[E_SIDE][0][iP]];
			d[i] += slau[iVar][iP].aw*x[neighbors_for_the_internal_node[W_SIDE][0][iP]];
			d[i] += slau[iVar][iP].an*x[neighbors_for_the_internal_node[N_SIDE][0][iP]];
			d[i] += slau[iVar][iP].as*x[neighbors_for_the_internal_node[S_SIDE][0][iP]];
			d[i] += rthdsd[iP]; //f.slau[iVar][iP].b;
		}

		iP = neighbors_for_the_internal_node[T_SIDE][0][id[2][j][n - 3]]; // номер граничного узла
		a[n - 1] = slau_bon[iVar][iP - maxelm].aw;//alphaP;
		if (slau_bon[iVar][iP - maxelm].iI>-1) {
			c[n - 1] = slau_bon[iVar][iP - maxelm].ai;
			a[n - 1] = a[n - 1] / alphaPb;
		}
		else c[n - 1] = 0.0;
		d[n - 1] = rthdsd[iP]; //slau_bon[iVar][iP-maxelm].b;
		ind[n - 1] = iP;

		// собственно сама прогонка или алгоритм Томаса:
		P[0] = b[0] / a[0];
		Q[0] = d[0] / a[0];
		doublereal rdiv;
		for (integer i = 1; i<n; ++i) {
			rdiv = 1.0 / (a[i] - c[i] * P[i - 1]);
			P[i] = b[i] * rdiv;
			Q[i] = (d[i] + c[i] * Q[i - 1])*rdiv;
		}
		//x[ind[n - 1]] = Q[n - 1];
		//for (integer i = n - 2; i >= 0; i--) {
			//x[ind[i]] = P[i] * x[ind[i + 1]] + Q[i];
		//}
		xc[ind[n - 1]] = Q[n - 1];
		for (integer i = n - 2; i >= 0; i--) {
			xc[ind[i]] = P[i] * xc[ind[i + 1]] + Q[i];
		}
		for (integer i = n - 2; i >= 1; i--) {
			// завышае новые значения.
			x[ind[i]] = omega*xc[ind[i]] + (1.0 - omega)*x[ind[i]];
		}
		// Граничные узлы не подвержены процессу релаксации.
		x[ind[0]] = xc[ind[0]];
		x[ind[n - 1]] = xc[ind[n - 1]];

		// Освобождение оперативной памяти
		//if (a != nullptr) {
			// оператор delete может быть применен повторно в том числе и к нулевому указателю. Проверки на null излишни.
			delete[] a;
			a = nullptr;
		//}
		//if (b != nullptr) {
			// оператор delete может быть применен повторно в том числе и к нулевому указателю. Проверки на null излишни.
			delete[] b;
			b = nullptr;
		//}
		//if (c != nullptr) {
			// оператор delete может быть применен повторно в том числе и к нулевому указателю. Проверки на null излишни.
			delete[] c;
			c = nullptr;
		//}
		//if (d != nullptr) {
			// оператор delete может быть применен повторно в том числе и к нулевому указателю. Проверки на null излишни.
			delete[] d;
			d = nullptr;
		//}
		// прогоночные коэффициенты
		//if (P != nullptr) {
			// оператор delete может быть применен повторно в том числе и к нулевому указателю. Проверки на null излишни.
			delete[] P;
			P = nullptr;
		//}
		//if (Q != nullptr) {
			// оператор delete может быть применен повторно в том числе и к нулевому указателю. Проверки на null излишни.
			delete[] Q;
			Q = nullptr;
		//}
		// глобальный номер перменной (связь)
		//if (ind != nullptr) {
			// оператор delete может быть применен повторно в том числе и к нулевому указателю. Проверки на null излишни.
			delete[] ind;
			ind = nullptr;
		//}

	} // for

	//if (xc != nullptr) {
		// оператор delete может быть применен повторно в том числе и к нулевому указателю. Проверки на null излишни.
		delete[] xc; // Освобождение оперативной памяти
		xc = nullptr;
	//}

} // solveLR1_serial


// полилинейный солвер LR1:
// Л.Н. Фомина Кемеровский государственный университет, Кемерово
// О повышении эффективности полилинейного рекурретного метода решения 
// разностных эллиптических уравнений.
void solveLR1_temp(TEMPER &t, doublereal* &x, doublereal* &rthdsd, integer ns) {

	// Каждый источник тепла пересекающий сеточную линию перпендикулярную источнику
	// разбивает сеточную линию на две подлинии. Это нужно учесть в программе. TODO.

    doublereal alphaP=1.0;
	
	// ns - размерность вектора x и правой части rthdsd.
	// f.slau - хранитель матрцы СЛАУ.
	// x - вектор решения.
	// rthdsd - вектор правой части.

	//f.a, f.b, f.c, f.d - трёхдиагональная матрица СЛАУ и правая часть
	// f.P, f.Q -  прогоночные коэффициенты
	// f.ind - глобальный номер переменной

    // Если стоит условие Дирихле то нижней релаксации нет, а если стоит условие Неймана то нижняя релаксация есть.

	// WE direct
	
	doublereal* xc=new doublereal[static_cast<integer>(t.maxelm)+ static_cast<integer>(t.maxbound)];

	#pragma omp parallel for shared (x, xc) schedule (guided) // 
	for (integer i1=0; i1<t.maxelm+t.maxbound; i1++) {
		xc[i1]=x[i1]; // copy x
	}

	// проход по всем линиям:
	#pragma omp parallel for shared (t, rthdsd, alphaP, x, xc) schedule (guided) // 
	for (integer j=0; j<t.iWE; ++j) {		

		integer n=t.rootWE[j].iN; // инициализация, в случае разрыва связей это значение переопределяется.

		// Вставка данных для паралельной обработки:
		doublereal *a=nullptr, *b=nullptr, *c=nullptr, *d=nullptr; // трёхдиагональная матрица и правая часть
	    doublereal *P=nullptr, *Q=nullptr; // прогоночные коэффициенты
	    integer *ind=nullptr; // связь с глобальной нумерацией узлов.
		// трехдиагональная матрица
	    a=new doublereal[n];
	    b=new doublereal[n];
	    c=new doublereal[n];
	    d=new doublereal[n];
	    // прогоночные коэффициенты
	    P=new doublereal[n];
	    Q=new doublereal[n];
	    // глобальный номер перменной (связь)
	    ind= new integer[n];

		integer iP=-1, i=-1;

		if ((a != nullptr) && (b != nullptr) && (c != nullptr) && (d != nullptr) && (P != nullptr) && (Q != nullptr) && (ind != nullptr)) {

			

			NODELR* r2 = t.rootWE[j].root;


			if (!t.rootWE[j].bNeimanStart) {
				iP = t.neighbors_for_the_internal_node[W_SIDE][0][r2->id]; // номер граничного узла
				a[0] = t.slau_bon[iP - t.maxelm].aw;///alphaP; // граничные узлы не релаксируются если на них стоит условие Дирихле
				if (t.slau_bon[iP - t.maxelm].iI > -1) {
					b[0] = t.slau_bon[iP - t.maxelm].ai;
					a[0] = a[0] / alphaP;
				}
				else b[0] = 0.0;
				d[0] = rthdsd[iP]; //f.slau_bon[iVar][iP-f.maxelm].b;
				ind[0] = iP;
				c[0] = 0.0;
			}
			else {
				//printf("bug\n"); system("pause"); // debug
				// разорвана связь с источником
				iP = r2->id;
				ind[0] = iP;
				a[0] = t.slau[iP].ap / alphaP;
				b[0] = t.slau[iP].ae;
				c[0] = 0.0; //t.slau[iP].aw; // связь была разорвана
				d[0] = t.slau[iP].an*xc[t.neighbors_for_the_internal_node[N_SIDE][0][iP]];
				d[0] += t.slau[iP].as*xc[t.neighbors_for_the_internal_node[S_SIDE][0][iP]];
				d[0] += t.slau[iP].at*xc[t.neighbors_for_the_internal_node[T_SIDE][0][iP]];
				d[0] += t.slau[iP].ab*xc[t.neighbors_for_the_internal_node[B_SIDE][0][iP]];
				d[0] += rthdsd[iP]; //f.slau[iVar][iP].b;

				r2 = r2->next;
			}

			i = 1;

			while (r2->next != nullptr) {
				iP = r2->id;

				if (i < t.rootWE[j].iN) {
					ind[i] = iP;
					a[i] = t.slau[iP].ap / alphaP;
					b[i] = t.slau[iP].ae;
					c[i] = t.slau[iP].aw;
					d[i] = t.slau[iP].an*xc[t.neighbors_for_the_internal_node[N_SIDE][0][iP]];
					d[i] += t.slau[iP].as*xc[t.neighbors_for_the_internal_node[S_SIDE][0][iP]];
					d[i] += t.slau[iP].at*xc[t.neighbors_for_the_internal_node[T_SIDE][0][iP]];
					d[i] += t.slau[iP].ab*xc[t.neighbors_for_the_internal_node[B_SIDE][0][iP]];
					d[i] += rthdsd[iP]; //f.slau[iVar][iP].b;
				}
				else {
					printf("error 1: i >= t.rootWE[j].iN in my_LR.c file in solveLR1_temp.\n");
					system("pause");
					exit(1);
				}

				r2 = r2->next;


				i++;
			}



			if (!t.rootWE[j].bNeimanEnd) {
				if (r2 != nullptr) {
					iP = r2->id;
					if (i < t.rootWE[j].iN) {
						ind[i] = iP;
						a[i] = t.slau[iP].ap / alphaP;
						b[i] = t.slau[iP].ae;
						c[i] = t.slau[iP].aw;
						d[i] = t.slau[iP].an*xc[t.neighbors_for_the_internal_node[N_SIDE][0][iP]];
						d[i] += t.slau[iP].as*xc[t.neighbors_for_the_internal_node[S_SIDE][0][iP]];
						d[i] += t.slau[iP].at*xc[t.neighbors_for_the_internal_node[T_SIDE][0][iP]];
						d[i] += t.slau[iP].ab*xc[t.neighbors_for_the_internal_node[B_SIDE][0][iP]];
						d[i] += rthdsd[iP]; //f.slau[iVar][iP].b;
					}
					else {
						printf("error 2: i >= t.rootWE[j].iN in my_LR.c file in solveLR1_temp.\n");
						system("pause");
						exit(1);
					}


					i++;
					iP = t.neighbors_for_the_internal_node[E_SIDE][0][r2->id]; // номер граничного узла
					a[i] = t.slau_bon[iP - t.maxelm].aw;///alphaP; // граничные узлы не релаксируются
					if (t.slau_bon[iP - t.maxelm].iI > -1) {
						c[i] = t.slau_bon[iP - t.maxelm].ai;
						a[i] = a[i] / alphaP;
					}
					else c[i] = 0.0;
					d[i] = rthdsd[iP]; //t.slau_bon[iP-t.maxelm].b;
					ind[i] = iP;
					b[i] = 0.0;
				}
				else {
					printf("error 3: i >= t.rootWE[j].iN in my_LR.c file in solveLR1_temp.\n");
					system("pause");
					exit(1);
				}
			}
			else
			{
				//printf("bug\n"); system("pause"); // debug
				// Этот случай соответствует 
				// разорванной с источником связи из жидкости.
				if (r2 != nullptr) {
					iP = r2->id;
					if (i < t.rootWE[j].iN) {
						ind[i] = iP;
						a[i] = t.slau[iP].ap / alphaP;
						b[i] = 0.0; // t.slau[iP].ae; связь разорвана
						c[i] = t.slau[iP].aw;
						d[i] = t.slau[iP].an*xc[t.neighbors_for_the_internal_node[N_SIDE][0][iP]];
						d[i] += t.slau[iP].as*xc[t.neighbors_for_the_internal_node[S_SIDE][0][iP]];
						d[i] += t.slau[iP].at*xc[t.neighbors_for_the_internal_node[T_SIDE][0][iP]];
						d[i] += t.slau[iP].ab*xc[t.neighbors_for_the_internal_node[B_SIDE][0][iP]];
						d[i] += rthdsd[iP]; //f.slau[iVar][iP].b;
					}
					else {
						printf("error: i >= t.rootWE[j].iN in my_LR.c file in solveLR1_temp.\n");
						system("pause");
						exit(1);
					}
				}
				else {
					printf("error 2: i >= t.rootWE[j].iN in my_LR.c file in solveLR1_temp.\n");
					system("pause");
					exit(1);
				}
			}
			n = i + 1;

			r2 = nullptr; // этот указатель больше не нужен

			// собственно сама прогонка или алгоритм Томаса:
			P[0] = b[0] / a[0];
			Q[0] = d[0] / a[0];
			doublereal rdiv;

			if (n <= t.rootWE[j].iN) {

			    for (i = 1; i < n; ++i) {
				    rdiv = 1.0 / (a[i] - c[i] * P[i - 1]);
				    P[i] = b[i] * rdiv;
				    Q[i] = (d[i] + c[i] * Q[i - 1])*rdiv;
			    }			

				x[ind[n - 1]] = Q[n - 1];
				for (i = n - 2; i >= 0; i--) {
					x[ind[i]] = P[i] * x[ind[i + 1]] + Q[i];
				}
			}
			else {
				printf("error: n>t.rootWE[j].iN\n");
				printf("see function solveLR1_temp in my_LR.c file.\n");
				system("pause");
				exit(1);
			}

			// Освобождение оперативной памяти
			//if (a != nullptr) {
				// оператор delete может быть применен повторно в том числе и к нулевому указателю. Проверки на null излишни.
				delete[] a;
				a = nullptr;
			//}
			//if (b != nullptr) {
			// оператор delete может быть применен повторно в том числе и к нулевому указателю. Проверки на null излишни.
				delete[] b;
				b = nullptr;
			//}
			//if (c != nullptr) {
				// оператор delete может быть применен повторно в том числе и к нулевому указателю. Проверки на null излишни.
				delete[] c;
				c = nullptr;
			//}
			//if (d != nullptr) {
				// оператор delete может быть применен повторно в том числе и к нулевому указателю. Проверки на null излишни.
				delete[] d;
				d = nullptr;
			//}
			// прогоночные коэффициенты
			//if (P != nullptr) {
				// оператор delete может быть применен повторно в том числе и к нулевому указателю. Проверки на null излишни.
				delete[] P;
				P = nullptr;
			//}
			//if (Q != nullptr) {
			// оператор delete может быть применен повторно в том числе и к нулевому указателю. Проверки на null излишни.
				delete[] Q;
				Q = nullptr;
			//}
			// глобальный номер перменной (связь)
			//if (ind != nullptr) {
				// оператор delete может быть применен повторно в том числе и к нулевому указателю. Проверки на null излишни.
				delete[] ind;
				ind = nullptr;
			//}
		}

	} // while


	// SN direct

	#pragma omp parallel for shared (x, xc) schedule (guided) // 
	for (integer i1=0; i1<t.maxelm+t.maxbound; i1++) {
		xc[i1]=x[i1]; // copy x
	}
	
	// проход по всем линиям:
	#pragma omp parallel for shared (t, rthdsd, alphaP, x, xc) schedule (guided) // 
	for (integer j=0; j<t.iSN; ++j) {
		// printf("Ok\n"); // debug

		integer n=t.rootSN[j].iN; // инициализация в случае разрыва связей это значение переопределяется.

		// Вставка данных для паралельной обработки:
		doublereal *a=nullptr, *b=nullptr, *c=nullptr, *d=nullptr; // трёхдиагональная матрица и правая часть
	    doublereal *P=nullptr, *Q=nullptr; // прогоночные коэффициенты
	    integer *ind=nullptr; // связь с глобальной нумерацией узлов.
		// трехдиагональная матрица
	    a=new doublereal[n];
	    b=new doublereal[n];
	    c=new doublereal[n];
	    d=new doublereal[n];
	    // прогоночные коэффициенты
	    P=new doublereal[n];
	    Q=new doublereal[n];
	    // глобальный номер перменной (связь)
	    ind= new integer[n];

		integer iP=-1, i=-1;

		if ((a != nullptr) && (b != nullptr) && (c != nullptr) && (d != nullptr) && (P != nullptr) && (Q != nullptr) && (ind != nullptr)) {

			

			NODELR* r2 = t.rootSN[j].root;

			if (!t.rootSN[j].bNeimanStart) {
				iP = t.neighbors_for_the_internal_node[S_SIDE][0][r2->id]; // номер граничного узла
				a[0] = t.slau_bon[iP - t.maxelm].aw;///alphaP;
				if (t.slau_bon[iP - t.maxelm].iI > -1) {
					b[0] = t.slau_bon[iP - t.maxelm].ai;
					a[0] = a[0] / alphaP;
				}
				else b[0] = 0.0;
				d[0] = rthdsd[iP]; //f.slau_bon[iVar][iP-f.maxelm].b;
				ind[0] = iP;
				c[0] = 0.0;// т.к. граничные условия:
			}
			else {
				//printf("bug\n"); system("pause"); // debug
				// разорвана связь с источником
				iP = r2->id;
				ind[0] = iP;
				a[0] = t.slau[iP].ap / alphaP;
				b[0] = t.slau[iP].an;
				c[0] = 0.0; //t.slau[iP].as; // разорвана связь с источником
				d[0] = t.slau[iP].ae*xc[t.neighbors_for_the_internal_node[E_SIDE][0][iP]];
				d[0] += t.slau[iP].aw*xc[t.neighbors_for_the_internal_node[W_SIDE][0][iP]];
				d[0] += t.slau[iP].at*xc[t.neighbors_for_the_internal_node[T_SIDE][0][iP]];
				d[0] += t.slau[iP].ab*xc[t.neighbors_for_the_internal_node[B_SIDE][0][iP]];
				d[0] += rthdsd[iP]; //f.slau[iVar][iP].b;

				r2 = r2->next;
			}

			i = 1;

			while (r2->next != nullptr) {
				iP = r2->id;
				if (i < t.rootSN[j].iN) {
					ind[i] = iP;
					a[i] = t.slau[iP].ap / alphaP;
					b[i] = t.slau[iP].an;
					c[i] = t.slau[iP].as;
					d[i] = t.slau[iP].ae*xc[t.neighbors_for_the_internal_node[E_SIDE][0][iP]];
					d[i] += t.slau[iP].aw*xc[t.neighbors_for_the_internal_node[W_SIDE][0][iP]];
					d[i] += t.slau[iP].at*xc[t.neighbors_for_the_internal_node[T_SIDE][0][iP]];
					d[i] += t.slau[iP].ab*xc[t.neighbors_for_the_internal_node[B_SIDE][0][iP]];
					d[i] += rthdsd[iP]; //f.slau[iVar][iP].b;
				}
				else {
					printf("i>=t.rootSN[j].iN in solveLR1_temp in my_LR.c\n");
					printf("model incorrect\n");
					system("pause");
					exit(1);
				}

				r2 = r2->next;


				i++;
			}


			if (!t.rootSN[j].bNeimanEnd) {
				if (r2 != nullptr) {
					iP = r2->id;
					if (i<t.rootSN[j].iN) {
						ind[i] = iP;
						a[i] = t.slau[iP].ap / alphaP;
						b[i] = t.slau[iP].an;
						c[i] = t.slau[iP].as;
						d[i] = t.slau[iP].ae*xc[t.neighbors_for_the_internal_node[E_SIDE][0][iP]];
						d[i] += t.slau[iP].aw*xc[t.neighbors_for_the_internal_node[W_SIDE][0][iP]];
						d[i] += t.slau[iP].at*xc[t.neighbors_for_the_internal_node[T_SIDE][0][iP]];
						d[i] += t.slau[iP].ab*xc[t.neighbors_for_the_internal_node[B_SIDE][0][iP]];
						d[i] += rthdsd[iP]; //f.slau[iVar][iP].b;
						i++;
						iP = t.neighbors_for_the_internal_node[N_SIDE][0][r2->id]; // номер граничного узла
						a[i] = t.slau_bon[iP - t.maxelm].aw;///alphaP;
						if (t.slau_bon[iP - t.maxelm].iI > -1) {
							c[i] = t.slau_bon[iP - t.maxelm].ai;
							a[i] = a[i] / alphaP;
						}
						else c[i] = 0.0;
						d[i] = rthdsd[iP]; //f.slau_bon[iVar][iP-f.maxelm].b;
						//printf("source=%e\n",t.d[i]); system("pause"); // debug
						ind[i] = iP;
						b[i] = 0.0; // т.к. граничное условие
					}
				}
				else {
					printf("error 2: r2==nullptr: in solveLR1_temp in my_LR.c file\n");
					system("pause");
					exit(1);
				}
			}
			else {
				//printf("bug\n"); system("pause"); // debug
				if (r2 != nullptr) {
					iP = r2->id;
					if (i < t.rootSN[j].iN) {
						ind[i] = iP;
						a[i] = t.slau[iP].ap / alphaP;
						b[i] = 0.0; //t.slau[iP].an; разрыв связи с источником
						c[i] = t.slau[iP].as;
						d[i] = t.slau[iP].ae*xc[t.neighbors_for_the_internal_node[E_SIDE][0][iP]];
						d[i] += t.slau[iP].aw*xc[t.neighbors_for_the_internal_node[W_SIDE][0][iP]];
						d[i] += t.slau[iP].at*xc[t.neighbors_for_the_internal_node[T_SIDE][0][iP]];
						d[i] += t.slau[iP].ab*xc[t.neighbors_for_the_internal_node[B_SIDE][0][iP]];
						d[i] += rthdsd[iP]; //f.slau[iVar][iP].b;
					}
				}
				else {
					printf("error r2==nullptr in solveLR1_temp in my_LR.c file\n");
					system("pause");
					exit(1);
				}
			}
			n = i + 1;

			r2 = nullptr; // этот указатель больше не нужен


			// собственно сама прогонка или алгоритм Томаса:
			P[0] = b[0] / a[0];
			Q[0] = d[0] / a[0];
			doublereal rdiv=0.0;
			if (n <= t.rootSN[j].iN) {

			    for (i = 1; i < n; ++i) {
				    rdiv = 1.0 / (a[i] - c[i] * P[i - 1]);
				    P[i] = b[i] * rdiv;
				    Q[i] = (d[i] + c[i] * Q[i - 1])*rdiv;
			    }
			    // printf("Ok2\n"); // debug
			
				if (ind[n - 1] < t.maxelm + t.maxbound) {
					x[ind[n - 1]] = Q[n - 1];
				}
				else {
					printf("my_LR.c solveLR1_temp ind[n-1]>=t.maxelm+t.maxbound\n");
					system("pause");
					exit(1);
				}
			
			    //x[t.ind[n-1]]=(t.d[n-1]+t.c[n-1]*t.Q[n-2])/(t.a[n-1]-t.c[n-1]*t.P[n-2]); // равносильно значению t.Q[n-1].
			    //printf("T==%e, %e, d==%e\n",x[t.ind[n-1]],t.Q[n-1],t.d[n-1]); system("pause"); // debug
			    for (i = n - 2; i >= 0; i--) {
			    	x[ind[i]] = P[i] * x[ind[i + 1]] + Q[i];
			    }
			}
			else {
				printf("error: n>t.rootSN[j].iN  in solveLR1_temp in my_LR.c file\n");
				system("pause");
				exit(1);
			}

			// Освобождение оперативной памяти
			//if (a != nullptr) {
				// оператор delete может быть применен повторно в том числе и к нулевому указателю. Проверки на null излишни.
				delete[] a;
				a = nullptr;
			//}
			//if (b != nullptr) {
			// оператор delete может быть применен повторно в том числе и к нулевому указателю. Проверки на null излишни.
				delete[] b;
				b = nullptr;
			//}
			//if (c != nullptr) {
			// оператор delete может быть применен повторно в том числе и к нулевому указателю. Проверки на null излишни.
				delete[] c;
				c = nullptr;
			//}
			//if (d != nullptr) {
				// оператор delete может быть применен повторно в том числе и к нулевому указателю. Проверки на null излишни.
				delete[] d;
				d = nullptr;
			//}
			// прогоночные коэффициенты
			//if (P != nullptr) {
				// оператор delete может быть применен повторно в том числе и к нулевому указателю. Проверки на null излишни.
				delete[] P;
				P = nullptr;
			//}
			//if (Q != nullptr) {
				// оператор delete может быть применен повторно в том числе и к нулевому указателю. Проверки на null излишни.
				delete[] Q;
				Q = nullptr;
			//}
			// глобальный номер перменной (связь)
			//if (ind != nullptr) {
				// оператор delete может быть применен повторно в том числе и к нулевому указателю. Проверки на null излишни.
				delete[] ind;
				ind = nullptr;
			//}
		}

	} // for


	// BT direct

	#pragma omp parallel for shared (x, xc) schedule (guided) // 
	for (integer i1=0; i1<t.maxelm+t.maxbound; i1++) {
		xc[i1]=x[i1]; // copy x
	}

	// проход по всем линиям:
	#pragma omp parallel for shared (t, rthdsd, alphaP, x, xc) schedule (guided) // 
	for (integer j=0; j<t.iBT; ++j) {

		integer n=t.rootBT[j].iN; // инициализация. В случае разрыва связей это значеение будет переопределено.

		// Вставка данных для паралельной обработки:
		doublereal *a=nullptr, *b=nullptr, *c=nullptr, *d=nullptr; // трёхдиагональная матрица и правая часть
	    doublereal *P=nullptr, *Q=nullptr; // прогоночные коэффициенты
	    integer *ind=nullptr; // связь с глобальной нумерацией узлов.
		// трехдиагональная матрица
	    a=new doublereal[n];
	    b=new doublereal[n];
	    c=new doublereal[n];
	    d=new doublereal[n];
	    // прогоночные коэффициенты
	    P=new doublereal[n];
	    Q=new doublereal[n];
	    // глобальный номер перменной (связь)
	    ind= new integer[n];

		integer iP = -1, i = -1;

		// После применения оператора new не требуется делать проверку на null.
		//if ((a != nullptr) && (b != nullptr) && (c != nullptr) && (d != nullptr) && (P != nullptr) && (Q != nullptr) && (ind != nullptr))
		{

			

			NODELR* r2 = t.rootBT[j].root;

			if (!t.rootBT[j].bNeimanStart) {
				iP = t.neighbors_for_the_internal_node[B_SIDE][0][r2->id]; // номер граничного узла
				a[0] = t.slau_bon[iP - t.maxelm].aw;///alphaP;
				if (t.slau_bon[iP - t.maxelm].iI > -1) {
					b[0] = t.slau_bon[iP - t.maxelm].ai;
					a[0] = a[0] / alphaP;
				}
				else b[0] = 0.0;
				d[0] = rthdsd[iP]; //f.slau_bon[iVar][iP-f.maxelm].b;
				ind[0] = iP;
				c[0] = 0.0; // т.к. граничные условия:

			}
			else {
				//printf("bug\n"); system("pause"); // debug
				// разрыв связи с источником тепла.
				iP = r2->id;
				ind[0] = iP;
				a[0] = t.slau[iP].ap / alphaP;
				b[0] = t.slau[iP].at;
				c[0] = 0.0; // t.slau[iP].ab; // разрыв связи с источником тепла
				d[0] = t.slau[iP].ae*xc[t.neighbors_for_the_internal_node[E_SIDE][0][iP]];
				d[0] += t.slau[iP].aw*xc[t.neighbors_for_the_internal_node[W_SIDE][0][iP]];
				d[0] += t.slau[iP].an*xc[t.neighbors_for_the_internal_node[N_SIDE][0][iP]];
				d[0] += t.slau[iP].as*xc[t.neighbors_for_the_internal_node[S_SIDE][0][iP]];
				d[0] += rthdsd[iP]; //f.slau[iVar][iP].b;

				r2 = r2->next;
			}

			i = 1;

			while (r2->next != nullptr) {
				iP = r2->id;
				if (i < n) {
					ind[i] = iP;
					a[i] = t.slau[iP].ap / alphaP;
					b[i] = t.slau[iP].at;
					c[i] = t.slau[iP].ab;
					d[i] = t.slau[iP].ae*xc[t.neighbors_for_the_internal_node[E_SIDE][0][iP]];
					d[i] += t.slau[iP].aw*xc[t.neighbors_for_the_internal_node[W_SIDE][0][iP]];
					d[i] += t.slau[iP].an*xc[t.neighbors_for_the_internal_node[N_SIDE][0][iP]];
					d[i] += t.slau[iP].as*xc[t.neighbors_for_the_internal_node[S_SIDE][0][iP]];
					d[i] += rthdsd[iP]; //f.slau[iVar][iP].b;
				}
				else {
					printf("error i>=n in my_LR.c in solveLR1_temp\n");
					system("pause");
					exit(1);
				}

				r2 = r2->next;


				i++;
			}

			if (i<n) {

				if (!t.rootBT[j].bNeimanEnd) {
					iP = r2->id;
					ind[i] = iP;
					a[i] = t.slau[iP].ap / alphaP;
					b[i] = t.slau[iP].at;
					c[i] = t.slau[iP].ab;
					d[i] = t.slau[iP].ae*xc[t.neighbors_for_the_internal_node[E_SIDE][0][iP]];
					d[i] += t.slau[iP].aw*xc[t.neighbors_for_the_internal_node[W_SIDE][0][iP]];
					d[i] += t.slau[iP].an*xc[t.neighbors_for_the_internal_node[N_SIDE][0][iP]];
					d[i] += t.slau[iP].as*xc[t.neighbors_for_the_internal_node[S_SIDE][0][iP]];
					d[i] += rthdsd[iP]; //f.slau[iVar][iP].b;
					i++;
					iP = t.neighbors_for_the_internal_node[T_SIDE][0][r2->id]; // номер граничного узла
					a[i] = t.slau_bon[iP - t.maxelm].aw;//alphaP;
					if (t.slau_bon[iP - t.maxelm].iI > -1) {
						c[i] = t.slau_bon[iP - t.maxelm].ai;
						a[i] = a[i] / alphaP;
					}
					else c[n - 1] = 0.0;
					d[i] = rthdsd[iP]; //f.slau_bon[iVar][iP-f.maxelm].b;
					ind[i] = iP;
					b[i] = 0.0; // т.к. граничные условия:
				}
				else
				{
					//printf("bug\n"); system("pause"); // debug
					iP = r2->id;
					ind[i] = iP;
					a[i] = t.slau[iP].ap / alphaP;
					b[i] = 0.0; //t.slau[iP].at; // разрыв связи с источником
					c[i] = t.slau[iP].ab;
					d[i] = t.slau[iP].ae*xc[t.neighbors_for_the_internal_node[E_SIDE][0][iP]];
					d[i] += t.slau[iP].aw*xc[t.neighbors_for_the_internal_node[W_SIDE][0][iP]];
					d[i] += t.slau[iP].an*xc[t.neighbors_for_the_internal_node[N_SIDE][0][iP]];
					d[i] += t.slau[iP].as*xc[t.neighbors_for_the_internal_node[S_SIDE][0][iP]];
					d[i] += rthdsd[iP]; //f.slau[iVar][iP].b;
				}
				n = i + 1;
			}
			else {
				printf("i index >=n error\n");
				printf("my_LR.c function solveLR1_temp\n");
				system("pause");
				exit(1);
			}

			r2 = nullptr; // этот указатель больше не нужен

			// собственно сама прогонка или алгоритм Томаса:
			P[0] = b[0] / a[0];
			Q[0] = d[0] / a[0];
			doublereal rdiv;

			if (n <= t.rootBT[j].iN) {

			    for (i = 1; i < n; ++i) {
				    rdiv = 1.0 / (a[i] - c[i] * P[i - 1]);
			    	P[i] = b[i] * rdiv;
			    	Q[i] = (d[i] + c[i] * Q[i - 1])*rdiv;
			    }


			
				if ((ind[n - 1] < t.maxelm + t.maxbound) && (ind[n - 1] >= 0)) {
					x[ind[n - 1]] = Q[n - 1];
				}
				else {
					printf("ind[n-1] ne leshit v dopustimjh predelax\n");
					printf("see function solveLR1_temp in my_LR.c file in source code\n");
					system("pause");
					exit(1);
				}
				for (i = n - 2; i >= 0; i--) {
					x[ind[i]] = P[i] * x[ind[i + 1]] + Q[i];
				}
			}
			else {
				printf("error: n>t.rootBT[j].iN in my_LR.c file in solveLR1_temp function\n");
				system("pause");
				exit(1);
			}


			// Освобождение оперативной памяти
			if (a != nullptr) {
				delete[] a;
			}
			if (b != nullptr) {
				delete[] b;
			}
			if (c != nullptr) {
				delete[] c;
			}
			if (d != nullptr) {
				delete[] d;
			}
			// прогоночные коэффициенты
			if (P != nullptr) {
				delete[] P;
			}
			if (Q != nullptr) {
				delete[] Q;
			}
			// глобальный номер перменной (связь)
			if (ind != nullptr) {
				delete[] ind;
			}

		}
		//else {
			//printf("memory no allocate for: a,b,c,d,P,Q,ind in my_LR.c file in function solveLR1_temp\n");
			//system("pause");
			//exit(1);
		//}

	} // for

	if (xc != nullptr) {
		delete[] xc; // освобождение оперативной памяти
	}
    
} // solveLR1_temp

doublereal fmaxloc(doublereal dA, doublereal dB) {
	doublereal dr=dA;
	if (dB>dA) dr=dB;
	return(dr);
}

// делаем несколько итераций полилинейного метода:
void solveLRn(doublereal* &x, doublereal* &rthdsd, integer ns, integer iVar, integer maxit, bool bprintf, bool bnorelax,
	int*** neighbors_for_the_internal_node, integer maxelm, equation3D** slau, equation3D_bon** slau_bon,
			   integer** iN, integer***id, integer iWE, integer iSN, integer iBT, doublereal *alpha, integer maxbound) {
    // bprintf - печатать диагнотику в файл или нет.
	// если полилинейный метод используется самостоятельно то печатать, а
	// если он используется как предобуславливатель то не печатать.

	

	bool only_serial = true;

	// переменная bdopfinsh включает зарекомендовавше себя дополнительное
	// условие выхода из цикла предобуславливания.
	// Если bdopfinish==false то дополнительное условие выключено.
	// Если bdopfinish==true то дополнительное условие работает.
	bool bdopfinish=true;

	integer i=0; // счётчик количества итераций
	if (bprintf) {
		printf("LR1 start solution...\n");
	}
	doublereal e=dterminatedTResudual;
	doublereal dresgl=1.0;
	doublereal divres=1.0; //  разница в невязках между двумя ближайшими итерациями
	// делает не менее одной итераций
	while (((i<maxit) && (!bdopfinish || ((fabs(dresgl)>e) && (divres>(e*1e-25))))) || (i<1)) {
		// При решении CFD задач, например радиатор водяного охлаждения 3л/мин было обнаружено что 
		// sor и даже Зейдель хуже для сходимости чем метод Якоби.
		if (only_serial) {
			doublereal omega = 1.0;
			switch (iVar) {
			case VELOCITY_X_COMPONENT:  omega = 1.0; break;
			case VELOCITY_Y_COMPONENT:  omega = 1.0;  break;
			case VELOCITY_Z_COMPONENT:  omega = 1.0; break;
			case PAM: omega = 1.7; break; // 1.855
			default:  omega = 1.0;  break;
			}
			omega = 0.7; // это очень хорошо (более быстрая сходимость в глобальном смысле).
			// нижняя релаксация 0.7 избавляет практически от всех всплесков на графиках сходимости cfd алгоритма.
			// если всплески всё ещё присутствуют то можно немного уменьшить параметр нижней релаксации omega.
			solveLR1_serial_sor(x, rthdsd, ns, iVar, bnorelax, neighbors_for_the_internal_node, maxelm, slau, slau_bon, iN, id, iWE, iSN, iBT, alpha, maxbound, omega);
		}
		else {
			// Метод Якоби да еще и распараллеленный оказывается куда лучше для сходимости.
			solveLR1(x, rthdsd, ns, iVar, bnorelax, neighbors_for_the_internal_node, maxelm, slau, slau_bon, iN, id, iWE, iSN, iBT, alpha, maxbound);
		}

		
		//doublereal rmax=-1.0; // значение невязки
		doublereal dsum=0.0;

		// Невязка согласованная с LR1sk солвером.
		doublereal* residual = nullptr;
		residual = new doublereal[maxelm + maxbound];
#pragma omp parallel for
		for (integer j=0; j<maxelm+maxbound; ++j) residual[j]=0.0; // инициализация.

		// Вычисляем невязку:
		// 1.
		// Внутренние контрольные объёмы:
        #pragma omp parallel for shared (maxelm, slau, neighbors_for_the_internal_node, rthdsd, x, residual) schedule (guided) //reduction (+:dsum)  
		for (integer j=0; j<maxelm; ++j) {
			doublereal dbuf=0.0;
			//rmax=fmaxloc(rmax,fabs(f.slau[iVar][j].ap*x[j]-f.slau[iVar][j].ab*x[f.neighbors_for_the_internal_node[BSIDE][0][j]]-f.slau[iVar][j].at*x[f.neighbors_for_the_internal_node[TSIDE][0][j]]-f.slau[iVar][j].an*x[f.neighbors_for_the_internal_node[NSIDE][0][j]]-f.slau[iVar][j].as*x[f.neighbors_for_the_internal_node[SSIDE][0][j]]-f.slau[iVar][j].ae*x[f.neighbors_for_the_internal_node[ESIDE][0][j]]-f.slau[iVar][j].aw*x[f.neighbors_for_the_internal_node[WSIDE][0][j]]-f.slau[iVar][j].b));
			dbuf = fabs(slau[iVar][j].ap*x[j] - slau[iVar][j].ab*x[neighbors_for_the_internal_node[B_SIDE][0][j]] - slau[iVar][j].at*x[neighbors_for_the_internal_node[T_SIDE][0][j]] - slau[iVar][j].an*x[neighbors_for_the_internal_node[N_SIDE][0][j]] - slau[iVar][j].as*x[neighbors_for_the_internal_node[S_SIDE][0][j]] - slau[iVar][j].ae*x[neighbors_for_the_internal_node[E_SIDE][0][j]] - slau[iVar][j].aw*x[neighbors_for_the_internal_node[W_SIDE][0][j]] - rthdsd[j]); //-f.slau[iVar][j].b
			//dbuf=dbuf*dbuf;
			//dsum+=dbuf;
			residual[j]=dbuf;
		}
		// 2.
		// Граничные контрольные объёмы:
		//doublereal dsumbuf=dsum;
		//dsum=0.0; // см. присваивание далее.
		#pragma omp parallel for shared (maxelm, maxbound, slau_bon, rthdsd, x, residual)  schedule (guided) //reduction (+:dsum) 
		for (integer j=maxelm; j<maxelm+maxbound; ++j) {
			doublereal dbuf=0.0;
			if (slau_bon[iVar][j-maxelm].iI>-1) {
				dbuf=fabs(slau_bon[iVar][j-maxelm].aw*x[j]-slau_bon[iVar][j-maxelm].ai*x[slau_bon[iVar][j-maxelm].iI]-rthdsd[j]); // -f.slau_bon[iVar][j-f.maxelm].b
			}
			else {
				dbuf=fabs(slau_bon[iVar][j-maxelm].aw*x[j]-rthdsd[j]); // -f.slau_bon[iVar][j-f.maxelm].b
			}
			//dbuf=dbuf*dbuf;
			/*// Если потребуется нижняя релаксация к граничным узлам то  для её проверки можно использовать этот код
			  if (dbuf >1e+6) { // debug

			  #if doubleintprecision == 1
					printf("Error node =%lld, maxelm=%lld\n",j,f.maxelm);
					printf("iI=%lld, aw=%e, ai=%e, b=%e\n",f.slau_bon[iVar][j-f.maxelm].iI,f.slau_bon[iVar][j-f.maxelm].aw,f.slau_bon[iVar][j-f.maxelm].ai,f.slau_bon[iVar][j-f.maxelm].b);
			  #else
					printf("Error node =%d, maxelm=%d\n",j,f.maxelm);
					printf("iI=%d, aw=%e, ai=%e, b=%e\n",f.slau_bon[iVar][j-f.maxelm].iI,f.slau_bon[iVar][j-f.maxelm].aw,f.slau_bon[iVar][j-f.maxelm].ai,f.slau_bon[iVar][j-f.maxelm].b);
			  #endif

				printf("rthdsd=%e, x=%e\n",rthdsd[j],x[j]);
				system("pause");
			}*/
			//dsum+=dbuf;
            residual[j]=dbuf;
		}
		//dsum+=dsumbuf;

		// нормы невязки по видимому должны быть согласованы с основным алгоритмом.
		dsum=NormaV(residual, maxelm+maxbound); // вычисление нормы невязки.
		
		delete[] residual; // освобождение памяти из под вектора невязки.
		residual = nullptr;

		if (bprintf) {
		    if (i%10==0) {
				std::cout<<"iter residual"<<std::endl;
			}
			std::cout << " " << i+1 << "  " << dsum << std::endl;   
		}
		divres=fabs(dresgl-dsum); // разница в невязках между двумя ближайшими итерациями
		dresgl=dsum;

		i++;
	}

} // solveLRn

// делаем несколько итераций полилинейного метода:
void solveLRn_temp(TEMPER &t, doublereal* &x, doublereal* &rthdsd, integer ns, integer maxit, bool bprintf) {
    // bprintf - печатать диагнотику в файл или нет.
	// если полилинейный метод используется самостоятельно то печатать, а
	// если он используется как предобуславливатель то не печатать.

	// переменная bdopfinsh включает зарекомендовавше себя дополнительное
	// условие выхода из цикла предобуславливания.
	// Если bdopfinish==false то дополнительное условие выключено.
	// Если bdopfinish==true то дополнительное условие работает.
	bool bdopfinish=true;

	// По видимому требуется ещё одно условие. смысл этого
	// условия в том чтобы мы не вышли за пределы вещественной арифметики, т.е.
	// чтобы не было переполнения. По ряду наблюдений невязка предобуславливателя 
	// равна квадрату невязки в алгоритме LR1sk. Нужно следить чтобы это значение 
	// не было лишком малым.


	integer i=0; // счётчик количества итераций
	if (bprintf) {
		printf("LR1 start solution...\n");
	}
	doublereal e=dterminatedTResudual; // Норма должна быть согласована с LR1sk солвером.
	doublereal dresgl=1.0;
	doublereal divres=1.0; //  разница в невязках между двумя ближайшими итерациями
	// делает не менее одной итераций
	while (((i<maxit) && (!bdopfinish || ((fabs(dresgl)>e) && (divres>(e*1e-25))))) || (i<1)) {
		solveLR1_temp(t, x, rthdsd, ns);
		
		//doublereal rmax=-1.0; // значение невязки
		doublereal dsum=0.0;
		
		// Невязка согласованная с LR1sk солвером.
		doublereal* residual=new doublereal[t.maxelm+t.maxbound];
#pragma omp parallel for
		for (integer j=0; j<t.maxelm+t.maxbound; ++j) residual[j]=0.0; // инициализация.

		// Вычисляем невязку:
		// 1.
		// Внутренние контрольные объёмы:
		#pragma omp parallel for shared (t,rthdsd,x,residual)  schedule (guided) //reduction (+:dsum)
		for (integer j=0; j<t.maxelm; ++j) {
			doublereal dbuf=0.0;

			//rmax=fmaxloc(rmax,fabs(f.slau[iVar][j].ap*x[j]-
			// f.slau[iVar][j].ab*x[f.neighbors_for_the_internal_node[BSIDE][0][j]]-
			// f.slau[iVar][j].at*x[f.neighbors_for_the_internal_node[TSIDE][0][j]]-
			// f.slau[iVar][j].an*x[f.neighbors_for_the_internal_node[NSIDE][0][j]]-
			// f.slau[iVar][j].as*x[f.neighbors_for_the_internal_node[SSIDE][0][j]]-
			// f.slau[iVar][j].ae*x[f.neighbors_for_the_internal_node[ESIDE][0][j]]-
			// f.slau[iVar][j].aw*x[f.neighbors_for_the_internal_node[WSIDE][0][j]]-
			// f.slau[iVar][j].b)); // или rthdsd[j]);
			dbuf=fabs(t.slau[j].ap*x[j]-
				t.slau[j].ab*x[t.neighbors_for_the_internal_node[B_SIDE][0][j]] -
				t.slau[j].at*x[t.neighbors_for_the_internal_node[T_SIDE][0][j]] -
				t.slau[j].an*x[t.neighbors_for_the_internal_node[N_SIDE][0][j]] -
				t.slau[j].as*x[t.neighbors_for_the_internal_node[S_SIDE][0][j]] -
				t.slau[j].ae*x[t.neighbors_for_the_internal_node[E_SIDE][0][j]] -
				t.slau[j].aw*x[t.neighbors_for_the_internal_node[W_SIDE][0][j]] -
					  rthdsd[j]);

			//dbuf=dbuf*dbuf;
			//dsum+=dbuf;
			residual[j]=dbuf;
		}
		//printf("sum=%e\n",dsum); system("pause");
		// 2.
		//doublereal dbufsum=dsum;
		//dsum=0.0;
		// Граничные контрольные объёмы:
		#pragma omp parallel for shared (t,rthdsd,x,residual)  schedule (guided) //reduction (+:dsum)
		for (integer j=t.maxelm; j<t.maxelm+t.maxbound; ++j) {
			
			doublereal dbuf=0.0;
			if (t.slau_bon[j-t.maxelm].iI>-1) {
				dbuf=fabs(t.slau_bon[j-t.maxelm].aw*x[j]-
					      t.slau_bon[j-t.maxelm].ai*x[t.slau_bon[j-t.maxelm].iI]-
						  rthdsd[j]);
			}
			else {
				dbuf=fabs(t.slau_bon[j-t.maxelm].aw*x[j]-rthdsd[j]);
			}
			//dbuf=dbuf*dbuf;
			/*// Если потребуется нижняя релаксация к граничным узлам то  для её проверки можно использовать этот код
			  if (dbuf >1e+6) { // debug
			  #if doubleintprecision == 1
					 printf("Error node =%lld, maxelm=%lld\n",j,f.maxelm);
					 printf("iI=%lld, aw=%e, ai=%e, b=%e\n",f.slau_bon[iVar][j-f.maxelm].iI,f.slau_bon[iVar][j-f.maxelm].aw,f.slau_bon[iVar][j-f.maxelm].ai,f.slau_bon[iVar][j-f.maxelm].b);
			  #else
					 printf("Error node =%d, maxelm=%d\n",j,f.maxelm);
				     printf("iI=%d, aw=%e, ai=%e, b=%e\n",f.slau_bon[iVar][j-f.maxelm].iI,f.slau_bon[iVar][j-f.maxelm].aw,f.slau_bon[iVar][j-f.maxelm].ai,f.slau_bon[iVar][j-f.maxelm].b);
			  #endif

				printf("rthdsd=%e, x=%e\n",rthdsd[j],x[j]);
				system("pause");
			}*/
			//dsum+=dbuf;
			residual[j]=dbuf;
		}
		//dsum+=dbufsum;

		// нормы невязки по видимому должны быть согласованы с основным алгоритмом.
		dsum=NormaV(residual, t.maxelm+t.maxbound); // вычисление нормы невязки.
		
		delete[] residual; // освобождение памяти из под вектора невязки.
		residual = nullptr;
		


		if (bprintf) {
			if (i % 10 == 0) {
				std::cout << "iter residual"<< std::endl;
			}
			std::cout << " " << i+1 << "  " << dsum << std::endl; //system("pause");
		}
		divres=fabs(dresgl-dsum); // разница в невязках между двумя ближайшими итерациями
		dresgl=dsum;

		i++;
	}

} // solveLRn_temp

#endif