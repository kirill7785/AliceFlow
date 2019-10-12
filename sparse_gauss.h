// файл sparse_gauss.h 
// Реализация метода Гаусса для разреженной матрицы.

#pragma once
#ifndef SPARSE_GAUSS_H
#define SPARSE_GAUSS_H 1

// Объявление модели вещественной арифметики содержится в самом начале программы в главно модуле 
// AliceFlow_v0_27
//#define doublereal double

// тип элемент
typedef struct tagTERM {
	integer key;
	doublereal val;
	// специальные поля для
	// АВЛ дерева, при реализации
	// интерфейса строк и стобцов
	// на массиве они не используются.
	// left и right убраны 10,10,2019
	//struct tagTERM *left=NULL;
	//struct tagTERM *right=NULL;
	integer bal;
} TERM;

// строка или столбец в матрице СЛАУ
typedef struct tagIRow{
    doublereal eps0; // для определения вещественного нуля

	// следующие два поля используются
	// только при реализации интерфейса 
	// строк и столбцов на динамическом массиве.
	integer POOL_SIZE; // размер массива включая нулевые элементы
	integer n; // число ненулевых элементов

	TERM *elm=NULL;

} IRow;

// разреженная матрица СЛАУ:
// Квадратная nxn с диагональным преобладанием,
// возможно несимметричная. Нумерация начинается 
// с нуля.
typedef struct tagIMatrix{
    doublereal eps0; // для определения вещественного нуля

	integer n; // размерность матрицы nxn.
    // jp - строки верхней полуматрицы
    // jm - столбцы нижней полуматрицы
    // dd - главная диагональ
	IRow *jp=NULL, *jm=NULL;
	doublereal *dd;
} IMatrix;

// поиск ячейки с индексом key  
integer search_i(TERM* list, integer n, integer key);

// удаляет элемент с ключём равным key
void deleteAt(TERM* &list, integer key, integer &n, integer &pool);

// добаляет элемент со значениями : num, val.
void add(TERM* &list, integer &n, integer &pool, integer num, doublereal val);

// возвращает значение ячейки ключ
// которой равен key
doublereal get_val(TERM* list, integer n, integer key);

// добавляет число value в ячейку
// с ключём key
void modify_add(TERM* &list, integer n, integer key, doublereal value);

// устанавливает число value в ячейку
// с ключём key
void modify_set(TERM* &list, integer n, integer key, doublereal value);

// зависит от внутреннего представления
// Возвращает все ненулевые элементы
void get_values(TERM *list, integer n, integer* &indexes, doublereal* &values);

// возвращает 0 при i<=0, i при i>0
integer sigma(integer i);

// преобразование row,col в координаты полуматриц (d,j)
// для верхней полуматрицы: d=0..(n-1), j=-(n-d-1)..-1
// для нижней полуматрицы:  d=0..(n-1), j=1..(n-d-1)
integer getD(integer row, integer col);
integer getJ(integer row, integer col);

// обратное преобразование координат полуматриц (d,j) в row,col
integer getRow(integer d, integer j);
integer getCol(integer d, integer j);


// записывает значение value в ячейку с индексом num
void setValueIRow(IRow *xO, integer num, doublereal value);

// добавляет value к существующему значению в ячейке num
void addValueIRow(IRow *xO, integer num, double value);

// возвращает значение ячейки num
doublereal getValueIRow(IRow *xO, integer num);

// возвращает все ненулевые ячейки строки/столбца: 
// индексы ячеек - в indexes, значения в values
integer getValuesIRow(IRow *xO, integer* &indexes, doublereal* &values);

// выделение памяти под разреженную матрицу
void initIMatrix(IMatrix *xO, integer n);

// освобождение памяти из под объекта
void freeIMatrix(IMatrix* xO);

// устанавливает значение value в ячейку с координатами [row,col];
// row - номер строки матрицы
// col - номер столбца матрицы
void setValueIMatrix(IMatrix *xO, integer row, integer col, doublereal value);

// добавляет значение value к ячейке [row,col]
void addValueIMatrix(IMatrix *xO, integer row, integer col, double value);

// возвращает значение ячейки [row,col]
doublereal  getValueIMatrix(IMatrix *xO, integer  row, integer  col);

// возвращает ненулевые значения и индексы ячеек строки d,
// которые находятся правее главной диагонали
integer  getJRowIMatrix(IMatrix *xO, integer  d, integer* &indexes, doublereal* &values);

// возвращает ненулевые значения и индексы ячеек столбца d, 
// которые находятся ниже главной диагонали
integer  getJColIMatrix(IMatrix *xO, integer  d, integer* &indexes, doublereal* &values);

// главный метод, возвращающий решение x,
// принимает вектор свободных членов b и 
// квадратную матрицу xO в специальном разреженном формате.
void calculateSPARSEgaussArray(IMatrix *xO, doublereal *x, doublereal *b);

// Преобразует матрицу в формате IMatrix в CSIR формат 
// совместимый с библиотекой ITL для реализации ILU разложения.
void convertIMatrixtoCSIR_ILU_ITL(IMatrix *xO, doublereal* &U_val, integer* &U_ind, integer* &U_ptr, doublereal* &L_val, integer* &L_ind, integer* &L_ptr);

#endif // !SPARSE_GAUSS_H