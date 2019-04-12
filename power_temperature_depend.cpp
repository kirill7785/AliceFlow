// В файле power_temperature_depend.cpp содержится
// реализация зависимости мощности источника тепла от
// текущего значения температуры и смещения стока.
// предполагается задание в виде таблицы для 
// набора дискретных значений в стиле:
//        0В      5В     10В    15В    20В    25В    30В
//  20    0.0    10Вт   20Вт   30Вт   35Вт   40Вт   50Вт 
//  30    0.0    11Вт   22Вт   31Вт   37Вт   ...
//  40    0.0    12Вт   24Вт   32Вт   ...
//  70    ...    ...    ...    ...
//  100
//  120
//  150
// Данная таблица была получена путём расчётов в другой программе.
// Где в самом левом столбце перечислены температуры всего intemp штук,
// в первой строке перечислены возможные смещения стоков в количестве inoffsetdrain.
// Таблица должна быть полностью заполнена пропуски не обрабатываются и вызывают сбой работы приложения.
// Остальные числа (Х) на пересечении соответствующего  строки и столбца это соответствующее
// данной температуре и смещению стока значение рассеиваемой в тепло мощности.
// Т.к. значения температуры и возможные значения смещения стока принимают свои значения непрерывно,
// то для любой заданной пары значений (температура, смещение стока) соответствующее значение рассеиваемой
// в тепло мощности предлагается  находить с помощью сплайновой интерполляции.
// А именно для каждой строки таблицы (каждому дискретному  значению температуры) предлагается сопоставить свой сплайн 
// аппроксимирующий значения в данной строке. Далее имея intemp таких сплайнов можно найти столбец значений соответствующих
// произвольному значению смещения стока. Тогда построив сплайн для этого нового столбца мощно определить искомую мощность
// в этом столбце соответствующую произвольному значению температуры.
// Задача кода приводимого в этом файле состоит в реализации функций которые смогут :
// 1. считывать таблицу из текстового файла.
// 2. Находить моменты сплайна для выделенной сеточной линии.
// 3. Находить значения сплайна в произвольной точке его определения зная значения эго моментов.
// 4. зная таблицу чисел power_table[i][j] находить значение рассеиваемой мощности, соответствующее
// заданной температуре и смещению стока.
// 
// Импульсный режим работы транзистора управляется смещением стока, например, от 0В до 30В.
// При нулевом смещении стока имеем нулевую рассеиваемую мощность. При положительном 
// смещении стока имеем положительную рассеиваемую в тепло мощность. Смещение стока,
// а значит и рассеиваемая в тепло мощность меняются с течением времени - имеем серию
// прямоугольных импульсов (соответствует работе в импульсном режиме с некоторой скважностью Q).
// На самом деле импульсы не прямоугольные а имеют форму трапеций : передний фронт 70нс, задний 100нс.
// Длительности импульса могут быть : 2мкс, 20мкс, 200мкс, 2000мкс. Смещение стока меняется только на фронтах,
// а в режиме молчания оно нулевое, а в течении длительности разогревающего импульса tau оно постоянно (равно константе)
// и равно например 30В > 0В. Я думаю фронтами можно принебречь и рассматривать прямоугольный импульс.
// begin 23 апреля 2012 года.

// Функции обработки сплайнов заимствованы из книги Г.З. Гарбера 
// Основы программирования на VBA Exel и численных методов.

#ifndef MY_POWER_TEMPERATYRE_DEPEND_CPP
#define MY_POWER_TEMPERATYRE_DEPEND_CPP 1

#include <stdlib.h>

// считывание таблицы из текстового файла :
// power_table.txt
void my_read_power_table(char* sname, integer &intemp, integer &inoffset_drain, 
	                     doublereal* &rtemp, doublereal* &roffset_drain,
						 doublereal** &rpower_table) 
{
	// sname - уникальное имя файла с табличной зависимостью.
	intemp=0; inoffset_drain=0; // эти значения останутся нулевыми если файл не сможет быть открыт.
	
	FILE *fpt=NULL; // файл из которого будет считываться таблица мощностей.
	errno_t err;
	err = fopen_s(&fpt, sname, "r");

	if ((err ) != 0) {
		// Файл открывается только для чтения,
		// если такого файла раньше не было то будет сообщено об ошибке приложения.
	    printf("Open File %s Error\n",sname);
       // system("PAUSE");
		system("pause");
        exit(0);
	}
	else {
		if (fpt != NULL) {
			// Структура считываемого файла следующая:
			// сначала идёт intemp пробел inoffset_drain переход на новую строку
			// с новой строки печатается таблица мощностей как обычно.

			//printf("comming...\n"); // debug

			integer din = 0;
			float fin = 0.0;


#ifdef MINGW_COMPILLER
			fscanf(fpt, "%d", &din);
			intemp = din; // количество различных дискретных значений температуры.
			fscanf(fpt, "%d", &din);
			inoffset_drain = din; // количество различных дискретных значений температуры.

#elseif sizeof(integer) == 4
			fscanf_s(fpt, "%d", &din);
			intemp = din; // количество различных дискретных значений температуры.
			fscanf_s(fpt, "%d", &din);
			inoffset_drain = din; // количество различных дискретных значений температуры.
#else
			fscanf_s(fpt, "%lld", &din);
			intemp = din; // количество различных дискретных значений температуры.
			fscanf_s(fpt, "%lld", &din);
			inoffset_drain = din; // количество различных дискретных значений температуры.
#endif


		//printf("%d %d \n",intemp,din); // debug

		// Выделение оперативной памяти под таблицу:
			rtemp = new doublereal[intemp];
			roffset_drain = new doublereal[inoffset_drain];
			rpower_table = new doublereal*[intemp];
			for (integer i = 0; i < intemp; i++) rpower_table[i] = new doublereal[inoffset_drain];

			// Считывание значений из файла :

			// Значения смещения стока.
			for (integer i = 0; i < inoffset_drain; i++) {

#ifdef MINGW_COMPILLER
				fscanf(fpt, "%f", &fin);
#else

				fscanf_s(fpt, "%f", &fin);
#endif

				roffset_drain[i] = fin;
			}

			// соответствующее значение температуры и 
			// список мощностей для разных значений смещения стока.
			for (integer j = 0; j < intemp; j++) {
				// Считывание таблицы построчно.
				for (integer i = 0; i <= inoffset_drain; i++) {

#ifdef MINGW_COMPILLER
					fscanf(fpt, "%f", &fin);
#else
					fscanf_s(fpt, "%f", &fin);
#endif
					if (i == 0) {
						rtemp[j] = fin;
					}
					else {
						rpower_table[j][i - 1] = fin;
					}
				}
			}

			fclose(fpt); // освобождение файлового дескриптора.
		}
	}
} // my_read_power_table

// печатает таблицу на консоль.
void my_print_table(integer intemp, integer inoffset_drain, 
	                     doublereal* rtemp, doublereal* roffset_drain,
						 doublereal** rpower_table) {
	printf("\t   ");
	for (integer i=0; i<inoffset_drain; i++) {
		printf("%.2e ",roffset_drain[i]);
	}
	printf("\n");

	for (integer i=0; i<intemp; i++) {
		printf("%.3e ",rtemp[i]);
		for (integer j=0; j<inoffset_drain; j++) printf("%.3e ",rpower_table[i][j]);
		printf("\n");
	}
}

// вычисление моментов сплайна по книге Г.З.Гарбера.
// стр. 238-246. "Основы программирования на VBA Exel и численных методов."
void mos(integer n, doublereal* &x, doublereal* &f, doublereal gamma_left, doublereal delta_left,
	     doublereal alpha_right, doublereal delta_right, doublereal* &mom) 
{

	if (n >= 4) {

		// имеем одномерную таблично заданную функцию одного вещественного
		// аргумента.
		// n - количество дискретных значений функции,
		// x - значения аргументов функции,
		// f - соответствующие значения функции.
		// mom - соответствующие моменты сплайна.
		// Для нахождения моментов имеем СЛАУ с трёх диагональной матрицей :
		// alpha[i]*mom[i-1]+2.0*mom[i]+gamma[i]*mom[i+1]=delta[i]; // i - ое уравнение 3-х диагональной СЛАУ.
		// с граничными условиями:
		// 2*mom[0]+gamma_left*mom[1]=delta_left; // левое граничное условие.
		// alpha_right*mom[n-2]+2.0*mom[n-1]=delta_right; // правое граничное условие.

		doublereal ralpha=0.0, rgamma=0.0, rdelta=0.0, rw=0.0; // r - doublereal;
		integer i=0;
		doublereal* rH = NULL;
		rH = new doublereal[n];
		for (i = 1; i < n; i++) rH[i] = x[i] - x[i - 1]; // шаги сетки
		doublereal* rP = NULL;
		doublereal* rQ = NULL;

		rP = new doublereal[n];
		rQ = new doublereal[n];

		// Прямая прогонка:
		rP[1] = -0.5*gamma_left;
		rQ[1] = 0.5*delta_left;
		for (i = 1; i < n - 1; i++) {
			rw = rH[i] + rH[i + 1];
			ralpha = rH[i] / rw;
			rgamma = 1.0 - ralpha;
			rdelta = 6.0*(((f[i + 1] - f[i]) / rH[i + 1]) - ((f[i] - f[i - 1]) / rH[i])) / rw;
			rw = ralpha*rP[i] + 2.0;
			rP[i + 1] = -rgamma / rw;
			rQ[i + 1] = (rdelta - ralpha*rQ[i]) / rw;
		} // end for i

		// Обратная прогонка:
		mom[n - 1] = (delta_right - alpha_right*rQ[n - 1]) / (alpha_right*rP[n - 1] + 2.0);
		for (i = n - 1; i > 0; i--) {
			mom[i - 1] = rP[i] * mom[i] + rQ[i];
		} // end for i

		// Освобождение оперативной памяти:
		// Не имеет смысла тестировать указатель на NULL так как память
		// была точно выделена оператором new.
		//if (rH != NULL) {
			delete[] rH;
		//}
		//if (rP != NULL) {
			delete[] rP;
		//}
		//if (rQ != NULL) {
			delete[] rQ;
		//}

	}
	else {
		printf("Error dimensional for mos function. \n");
		printf("parametr n<=3. error. \n");
		system("pause");
		exit(1);
	}
} // mos

// подпрограмма сплайновой интерполяции.
// Вычисляет значение сплайна в заданной точке.
// реализовано по книге Г.З.Гарбера.
// стр. 238-246. "Основы программирования на VBA Exel и численных методов."
void si(integer n, doublereal* x, doublereal* f, doublereal* mom, doublereal chi,
	    doublereal &s, doublereal &s1, doublereal &s2) {
	 // n - количество дискретных точек по которым построен сплайн,
     // x - массив значений аргументов функции,
     // f - массив значений функции,
	 // mom - массив моментов сплайна,
	 // chi - координата точки из области определения сплайна в которой нужно вычислить значение сплайна,
	 // s - рассчитываемое значение сплайна,
	 // s1, s2 - рассчитываемые значения первой и второй производных сплайна.

	integer i;
	doublereal rh, rhh, rh1, rh1h1, rh2, rh2h2; // r - doublereal
	// Нахождение элементарного отрезка, содержащего chi:
	for (i=1; i<n; i++) {
		if (x[i]>chi) break; // принудительный выход из цикла
	}
	if (i>n-1) i=n-1;
	// Расчёт значения сплайна в точке chi :
	rh=x[i]-x[i-1];
	rhh=rh*rh;
	rh1=chi-x[i-1];
	rh1h1=rh1*rh1;
	rh2=x[i]-chi;
	rh2h2=rh2*rh2;
	s=(mom[i-1]*rh2h2*rh2+mom[i]*rh1h1*rh1)/(6.0*rh)+((f[i-1]-mom[i-1]*rhh/6.0)*rh2+(f[i]-mom[i]*rhh/6.0)*rh1)/rh;
	// Расчёт первой производной сплайна в точке chi :
	s1=(-mom[i-1]*rh2h2+mom[i]*rh1h1)/(2.0*rh)+(f[i]-f[i-1])/rh-rh*(mom[i]-mom[i-1])/6.0;
	// Расчёт второй производной сплайна в точке chi :
	s2=(mom[i-1]*rh2+mom[i]*rh1)/rh;
} // si


// по заданной по точкам таблице мощностей, вычисление значения мощности
// в произвольной точке(rtemp_current, roffset_drain_current) путём 
// сплайновой интерполляции.
doublereal my_splain_interpol_power_table(integer intemp, integer inoffset_drain, 
	                     doublereal* rtemp, doublereal* roffset_drain, doublereal** rpower_table, 
						 doublereal rtemp_current, doublereal roffset_drain_current)  
{

	// Алгоритм:
	// сначала строятся сплайны для каждой строки таблицы.
	// Зная сплайн для строки таблицы производится вычисление значения сплайна
	// для заданного значения смещения стока в рамках данной строки, т.е. при заданной температуре.
	// проделая это для всех строк мы получем столбец соответствующий заданному смещению стока.
	// Проделав сплайновую интерполяцию для полученного столбца можно без труда вычислить значение
	// сплайна построенного по этому столбцу для заданной температуры - это и будет искомая мощность.

	doublereal rpower_curent;

	doublereal* m2=new doublereal[inoffset_drain];
	doublereal** mm=new doublereal*[intemp];
	for (integer i=0; i<intemp; i++) mm[i]=new doublereal[inoffset_drain];
	doublereal* f2=new doublereal[inoffset_drain];

	// Расчёт двумерного массива моментов по значениям смещения стока для всех строк таблицы
	// последовательно сверху вниз.
	for (integer i=0; i<intemp; i++) {
		for (integer j=0; j<inoffset_drain; j++) {
			// запоминаем строку мощностей в векторе f2:
			f2[j]=rpower_table[i][j];
		}  // end for j
		mos(inoffset_drain,roffset_drain,f2,0.0,0.0,0.0,0.0,m2);
		for (integer j=0; j<inoffset_drain; j++) {
			mm[i][j]=m2[j];
		} // end for j
	}


	doublereal* f1=new doublereal[intemp];
	for (integer i=0; i<intemp; i++) {
		// для всех значений температуры построчно
		// соответствующих заданному значению смещения стока:
		for (integer j=0; j<inoffset_drain; j++) {
			f2[j]=rpower_table[i][j];
			m2[j]=mm[i][j];
		}
		// вычисление стобца соответствующего заданному значению смещения стока.
		doublereal s1, s2;
		si(inoffset_drain, roffset_drain, f2, m2, roffset_drain_current,  f1[i], s1, s2); 
	}
	doublereal* m1=new doublereal[intemp];
	mos(intemp,rtemp,f1,0.0,0.0,0.0,0.0,m1);
	doublereal s1, s2;
	// Вычисление искомого значения мощности путём повторной сплайновой интерполяции по 
	// столбцу, полученному в результате интерполяции первого прохода.
    si(intemp, rtemp, f1, m1, rtemp_current,  rpower_curent, s1, s2);

	// Освобождение оперативной памяти:
	delete[] m2;
	delete[] f2;
	for (integer i=0; i<intemp; i++) delete[] mm[i];
	delete[] mm;
	delete[] f1; 
	delete[] m1;

	//printf("power is %e\n",rpower_curent); // debug

	return rpower_curent;
}

#endif
