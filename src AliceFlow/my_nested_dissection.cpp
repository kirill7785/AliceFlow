// Файл my_nested_dissection.cpp содержит код для 
// реализации алгоритма вложенных сечений. 17 апреля 2013 года.
// Алгоритм построения вложенных сечений должен использовать обход графа в ширину,
// если мы имеем дело с теорией графов. В терминах графов мы будем смотреть как движется фронт.
// И используя этот фронт строить вложенные сечения постепенно добавляя ячейки лежащие на одном фронте, 
// фронт за фронтом. Здесь же мы решили не связываться с теорией графов а использовать алгоритм вычисления расстояния до стенки.
// Данный модуль базируется на модуле shortest_distance.cpp.

// Суть вычисления кратчайшего расстояния до стенки: в приложении 2 книги А.Ю.Снегирёва
// Высокопроизводительные вычисления в технической физике.
// Численное моделирование турбулентных течений. 
// Санкт-Петербург, Издательство Политехнического университета 2009.
// говорится, что в ANSYS CFX для вычисления кратчайшего расстояния до
// ближайшей стенки решается эллиптическое уравнение для функции FI в правой
// части уравнения стоит -1. На твердой стенке расстояние до которой требуется
// вычислить стоит нулевое условие Дирихле, а на остальных границах которые не
// являются стенкой стоят однородные условия Неймана. Само кратчайшее расстояние 
// вычисляется как некая функция зависящая от градиентов найденной функции FI.
// Данное уравнение является эллиптическим с положительно определённой симметричной
// матрицей для обращения которой может быть применён ICCG solver.
// В данном модуле реализуется программный код находящий кратчайшее расстояние до стенки
// путём решения эллиптического уравнения в частных производных.
// begin 10 апреля 2012 года.

// Кратчайшее расстояние до стенки используется в частности в Zero Equation Turbulence Model.
// Кратчайшее расстояние используется также в алгоритме вложенных сечений.
// Цель.  Занумеровать расчётную область(получить массив inumerate), 
// чтобы параллельно выполнять обращение треугольных матриц и многое другое. 
// В общем данный модуль служит целям распараллеливания.

#ifndef MY_NESTED_DISSECTION_CPP
#define MY_NESTED_DISSECTION_CPP 1

// Нам важна не форма контрольного объёма а их количество
// и взаимное расположение.
const bool bunion_8_08_2015=true;

// Отключение экспорта результата разбиения в программу tecplot резко увеличивает
// быстродействие кода.
const bool debug_export_tecplot=false;

// Если мы делим первое деление (ilevel_info_nested_desection==1) то его надо сделать хорошо,
// далее всё менее и менее хорошо.
integer ilevel_info_nested_desection=1;

//#define FIDISTW 0 // потенциал FI
//#define GRADXFI 1 // и его градиенты
//#define GRADYFI 2
//#define GRADZFI 3

/* нам нужна глобальная область видимости для этой структуры данных.
// Структура границ деления:
typedef struct TPARBOUND {
	integer ileft_start, ileft_finish, iright_start, iright_finish, iseparate_start, iseparate_finish;
	bool active; // активность декомпозиции.
} PARBOUND;


// Структура данных используемая для распараллеливания.
typedef struct TPARDATA {
	integer ncore; // 1, 2, 4, 8.
	integer *inumerate;
	// это для ncore==2;
	PARBOUND b0;
	// это для ncore==4;
    PARBOUND b00, b01;
	// это для ncore==8;
	PARBOUND b000, b001, b010, b011;
} PARDATA;
*/
// Вычисление градиентов потенциала FI в центрах внутренних КО
// и на границах с помощью линейной интерполяции.
void green_gauss_FI_(integer iP, doublereal** &potent, int** nvtx, TOCHKA* pa,
	int*** neighbors_for_the_internal_node, integer maxelm, bool bbond) {
	// maxelm - число внутренних КО.
	// Вычисляет градиенты скоростей для внутренних КО.
	// если bbond   то будут вычислены значения в граничных КО, иначе только во внутренних.
    // Замечание во внутренних КО значения градиентов должны быть вычислены в первую очередь. Т.е.
	// необходимо совершить два последовательных запуска данной функции.

	// iP - номер внутреннего контрольного объёма
	// iP изменяется от 0 до maxelm-1.
	integer iE, iN, iT, iW, iS, iB; // номера соседних контрольных объёмов
	iE=neighbors_for_the_internal_node[E_SIDE][0][iP]; iN=neighbors_for_the_internal_node[N_SIDE][0][iP]; iT=neighbors_for_the_internal_node[T_SIDE][0][iP]; iW=neighbors_for_the_internal_node[W_SIDE][0][iP]; iS=neighbors_for_the_internal_node[S_SIDE][0][iP]; iB=neighbors_for_the_internal_node[B_SIDE][0][iP];

	// Если с одной из сторон стоит граница расчётной области
	// то соответствующая переменная равна true
	bool bE=false, bN=false, bT=false, bW=false, bS=false, bB=false;
    
	if (iE>=maxelm) bE=true;
	if (iN>=maxelm) bN=true;
	if (iT>=maxelm) bT=true;
    if (iW>=maxelm) bW=true;
	if (iS>=maxelm) bS=true;
	if (iB>=maxelm) bB=true;

	// вычисление размеров текущего контрольного объёма:
	doublereal dx=0.0, dy=0.0, dz=0.0;// объём текущего контрольного объёма
	volume3D(iP, nvtx, pa, dx, dy, dz);

	doublereal dxe=0.5*dx, dxw=0.5*dx, dyn=0.5*dy, dys=0.5*dy, dzt=0.5*dz, dzb=0.5*dz;
    // т.к. известна нумерация вершин куба, то здесь она используется
	// x - direction
    if (!bE) dxe=0.5*(pa[nvtx[1][iE]-1].x+pa[nvtx[0][iE]-1].x);
	if (!bE) dxe-=0.5*(pa[nvtx[1][iP]-1].x+pa[nvtx[0][iP]-1].x);
	if (!bW) dxw=0.5*(pa[nvtx[1][iP]-1].x+pa[nvtx[0][iP]-1].x);
	if (!bW) dxw-=0.5*(pa[nvtx[1][iW]-1].x+pa[nvtx[0][iW]-1].x);
    // y - direction
	if (!bN) dyn=0.5*(pa[nvtx[2][iN]-1].y+pa[nvtx[0][iN]-1].y);
	if (!bN) dyn-=0.5*(pa[nvtx[2][iP]-1].y+pa[nvtx[0][iP]-1].y);
	if (!bS) dys=0.5*(pa[nvtx[2][iP]-1].y+pa[nvtx[0][iP]-1].y);
	if (!bS) dys-=0.5*(pa[nvtx[2][iS]-1].y+pa[nvtx[0][iS]-1].y);
    // z - direction
	if (!bT) dzt=0.5*(pa[nvtx[4][iT]-1].z+pa[nvtx[0][iT]-1].z);
	if (!bT) dzt-=0.5*(pa[nvtx[4][iP]-1].z+pa[nvtx[0][iP]-1].z);
	if (!bB) dzb=0.5*(pa[nvtx[4][iP]-1].z+pa[nvtx[0][iP]-1].z);
	if (!bB) dzb-=0.5*(pa[nvtx[4][iB]-1].z+pa[nvtx[0][iB]-1].z);

	// Учёт неравномерности расчётной сетки:
	doublereal feplus, fwplus, fnplus, fsplus, ftplus, fbplus;
	// x-direction
	feplus=0.5*dx/dxe;
	fwplus=0.5*dx/dxw;
	// y-direction
	fnplus=0.5*dy/dyn;
	fsplus=0.5*dy/dys;
	// z-direction
	ftplus=0.5*dz/dzt;
	fbplus=0.5*dz/dzb;

	// линейная интерполяция скорости VX на грань КО.
    doublereal fe, fw, fn, fs, ft, fb;
	if (!bbond) {
		// внутренние КО.

	    if (!bE) fe=feplus*potent[FIDISTW][iE]+(1.0-feplus)*potent[FIDISTW][iP]; else fe=potent[FIDISTW][iE];
        if (!bW) fw=fwplus*potent[FIDISTW][iW]+(1.0-fwplus)*potent[FIDISTW][iP]; else fw=potent[FIDISTW][iW];
	    if (!bN) fn=fnplus*potent[FIDISTW][iN]+(1.0-fnplus)*potent[FIDISTW][iP]; else fn=potent[FIDISTW][iN];
        if (!bS) fs=fsplus*potent[FIDISTW][iS]+(1.0-fsplus)*potent[FIDISTW][iP]; else fs=potent[FIDISTW][iS];
        if (!bT) ft=ftplus*potent[FIDISTW][iT]+(1.0-ftplus)*potent[FIDISTW][iP]; else ft=potent[FIDISTW][iT];
        if (!bB) fb=fbplus*potent[FIDISTW][iB]+(1.0-fbplus)*potent[FIDISTW][iP]; else fb=potent[FIDISTW][iB];

	    potent[GRADXFI][iP]=(fe-fw)/dx;
	    potent[GRADYFI][iP]=(fn-fs)/dy;
	    potent[GRADZFI][iP]=(ft-fb)/dz;
	}
	else {
		// граничные узлы.
		// градиенты в граничных узлах восстанавливаются с помощью линейной интерполяции.

		if (bE) {
			potent[GRADXFI][iE]=potent[GRADXFI][iP]+(dxe/dxw)*(potent[GRADXFI][iP]-potent[GRADXFI][iW]);
			potent[GRADYFI][iE]=potent[GRADYFI][iP]+(dxe/dxw)*(potent[GRADYFI][iP]-potent[GRADYFI][iW]);
			potent[GRADZFI][iE]=potent[GRADZFI][iP]+(dxe/dxw)*(potent[GRADZFI][iP]-potent[GRADZFI][iW]);
		}

		if (bW) {
			potent[GRADXFI][iW]=potent[GRADXFI][iP]+(dxw/dxe)*(potent[GRADXFI][iP]-potent[GRADXFI][iE]);
			potent[GRADYFI][iW]=potent[GRADYFI][iP]+(dxw/dxe)*(potent[GRADYFI][iP]-potent[GRADYFI][iE]);
			potent[GRADZFI][iW]=potent[GRADZFI][iP]+(dxw/dxe)*(potent[GRADZFI][iP]-potent[GRADZFI][iE]);
		}

		if (bN) {
			potent[GRADXFI][iN]=potent[GRADXFI][iP]+(dyn/dys)*(potent[GRADXFI][iP]-potent[GRADXFI][iS]);
			potent[GRADYFI][iN]=potent[GRADYFI][iP]+(dyn/dys)*(potent[GRADYFI][iP]-potent[GRADYFI][iS]);
			potent[GRADZFI][iN]=potent[GRADZFI][iP]+(dyn/dys)*(potent[GRADZFI][iP]-potent[GRADZFI][iS]);
		}

		if (bS) {
			potent[GRADXFI][iS]=potent[GRADXFI][iP]+(dys/dyn)*(potent[GRADXFI][iP]-potent[GRADXFI][iN]);
			potent[GRADYFI][iS]=potent[GRADYFI][iP]+(dys/dyn)*(potent[GRADYFI][iP]-potent[GRADYFI][iN]);
			potent[GRADZFI][iS]=potent[GRADZFI][iP]+(dys/dyn)*(potent[GRADZFI][iP]-potent[GRADZFI][iN]);
		}

		if (bT) {
			potent[GRADXFI][iT]=potent[GRADXFI][iP]+(dzt/dzb)*(potent[GRADXFI][iP]-potent[GRADXFI][iB]);
			potent[GRADYFI][iT]=potent[GRADYFI][iP]+(dzt/dzb)*(potent[GRADYFI][iP]-potent[GRADYFI][iB]);
			potent[GRADZFI][iT]=potent[GRADZFI][iP]+(dzt/dzb)*(potent[GRADZFI][iP]-potent[GRADZFI][iB]);
		}

		if (bB) {
			potent[GRADXFI][iB]=potent[GRADXFI][iP]+(dzb/dzt)*(potent[GRADXFI][iP]-potent[GRADXFI][iT]);
			potent[GRADYFI][iB]=potent[GRADYFI][iP]+(dzb/dzt)*(potent[GRADYFI][iP]-potent[GRADYFI][iT]);
			potent[GRADZFI][iB]=potent[GRADZFI][iP]+(dzb/dzt)*(potent[GRADZFI][iP]-potent[GRADZFI][iT]);
		}
	}

} // green_gauss_FI

// Вычисление градиентов потенциала FI в центрах внутренних КО
// и на границах с помощью линейной интерполяции.
// Для процедуры nested desection нам важно количество контрольных объёмов а не их 
// геометрическая форма, поэтому мы считает что сетка единичная равномерная, т.е. 
// состоит из единичных кубов с единичными рёбрами.
// 8 августа 2015.
void green_gauss_FI_union(integer iP, doublereal** &potent, int** nvtx, TOCHKA* pa,
	int*** neighbors_for_the_internal_node, integer maxelm, bool bbond) {
	// maxelm - число внутренних КО.
	// Вычисляет градиенты скоростей для внутренних КО.
	// если bbond   то будут вычислены значения в граничных КО, иначе только во внутренних.
    // Замечание во внутренних КО значения градиентов должны быть вычислены в первую очередь. Т.е.
	// необходимо совершить два последовательных запуска данной функции.

	// iP - номер внутреннего контрольного объёма
	// iP изменяется от 0 до maxelm-1.
	integer iE, iN, iT, iW, iS, iB; // номера соседних контрольных объёмов
	iE=neighbors_for_the_internal_node[E_SIDE][0][iP]; iN=neighbors_for_the_internal_node[N_SIDE][0][iP]; iT=neighbors_for_the_internal_node[T_SIDE][0][iP]; iW=neighbors_for_the_internal_node[W_SIDE][0][iP]; iS=neighbors_for_the_internal_node[S_SIDE][0][iP]; iB=neighbors_for_the_internal_node[B_SIDE][0][iP];

	// Если с одной из сторон стоит граница расчётной области
	// то соответствующая переменная равна true
	bool bE=false, bN=false, bT=false, bW=false, bS=false, bB=false;
    
	if (iE>=maxelm) bE=true;
	if (iN>=maxelm) bN=true;
	if (iT>=maxelm) bT=true;
    if (iW>=maxelm) bW=true;
	if (iS>=maxelm) bS=true;
	if (iB>=maxelm) bB=true;

	// вычисление размеров текущего контрольного объёма:
	//doublereal dx=0.0, dy=0.0, dz=0.0;// объём текущего контрольного объёма
	//volume3D(iP, nvtx, pa, dx, dy, dz);
	/*
	doublereal dxe=0.5*dx, dxw=0.5*dx, dyn=0.5*dy, dys=0.5*dy, dzt=0.5*dz, dzb=0.5*dz;
    // т.к. известна нумерация вершин куба, то здесь она используется
	// x - direction
    if (!bE) dxe=0.5*(pa[nvtx[1][iE]-1].x+pa[nvtx[0][iE]-1].x);
	if (!bE) dxe-=0.5*(pa[nvtx[1][iP]-1].x+pa[nvtx[0][iP]-1].x);
	if (!bW) dxw=0.5*(pa[nvtx[1][iP]-1].x+pa[nvtx[0][iP]-1].x);
	if (!bW) dxw-=0.5*(pa[nvtx[1][iW]-1].x+pa[nvtx[0][iW]-1].x);
    // y - direction
	if (!bN) dyn=0.5*(pa[nvtx[2][iN]-1].y+pa[nvtx[0][iN]-1].y);
	if (!bN) dyn-=0.5*(pa[nvtx[2][iP]-1].y+pa[nvtx[0][iP]-1].y);
	if (!bS) dys=0.5*(pa[nvtx[2][iP]-1].y+pa[nvtx[0][iP]-1].y);
	if (!bS) dys-=0.5*(pa[nvtx[2][iS]-1].y+pa[nvtx[0][iS]-1].y);
    // z - direction
	if (!bT) dzt=0.5*(pa[nvtx[4][iT]-1].z+pa[nvtx[0][iT]-1].z);
	if (!bT) dzt-=0.5*(pa[nvtx[4][iP]-1].z+pa[nvtx[0][iP]-1].z);
	if (!bB) dzb=0.5*(pa[nvtx[4][iP]-1].z+pa[nvtx[0][iP]-1].z);
	if (!bB) dzb-=0.5*(pa[nvtx[4][iB]-1].z+pa[nvtx[0][iB]-1].z);
	*/

	doublereal dx=1.0, dy=1.0, dz=1.0;
	doublereal dxe=1.0, dxw=1.0, dyn=1.0, dys=1.0, dzt=1.0, dzb=1.0;
	// Учёт неравномерности расчётной сетки:
	doublereal feplus, fwplus, fnplus, fsplus, ftplus, fbplus;
	/*
	// x-direction
	feplus=0.5*dx/dxe;
	fwplus=0.5*dx/dxw;
	// y-direction
	fnplus=0.5*dy/dyn;
	fsplus=0.5*dy/dys;
	// z-direction
	ftplus=0.5*dz/dzt;
	fbplus=0.5*dz/dzb;
	*/
	feplus=0.5;
	fwplus=0.5;
	fnplus=0.5;
	fsplus=0.5;
	ftplus=0.5;
	fbplus=0.5;

	// линейная интерполяция скорости VX на грань КО.
    doublereal fe, fw, fn, fs, ft, fb;
	if (!bbond) {
		// внутренние КО.

	    if (!bE) fe=feplus*potent[FIDISTW][iE]+(1.0-feplus)*potent[FIDISTW][iP]; else fe=potent[FIDISTW][iE];
        if (!bW) fw=fwplus*potent[FIDISTW][iW]+(1.0-fwplus)*potent[FIDISTW][iP]; else fw=potent[FIDISTW][iW];
	    if (!bN) fn=fnplus*potent[FIDISTW][iN]+(1.0-fnplus)*potent[FIDISTW][iP]; else fn=potent[FIDISTW][iN];
        if (!bS) fs=fsplus*potent[FIDISTW][iS]+(1.0-fsplus)*potent[FIDISTW][iP]; else fs=potent[FIDISTW][iS];
        if (!bT) ft=ftplus*potent[FIDISTW][iT]+(1.0-ftplus)*potent[FIDISTW][iP]; else ft=potent[FIDISTW][iT];
        if (!bB) fb=fbplus*potent[FIDISTW][iB]+(1.0-fbplus)*potent[FIDISTW][iP]; else fb=potent[FIDISTW][iB];

	    potent[GRADXFI][iP]=(fe-fw)/dx;
	    potent[GRADYFI][iP]=(fn-fs)/dy;
	    potent[GRADZFI][iP]=(ft-fb)/dz;
	}
	else {
		// граничные узлы.
		// градиенты в граничных узлах восстанавливаются с помощью линейной интерполяции.

		if (bE) {
			potent[GRADXFI][iE]=potent[GRADXFI][iP]+(dxe/dxw)*(potent[GRADXFI][iP]-potent[GRADXFI][iW]);
			potent[GRADYFI][iE]=potent[GRADYFI][iP]+(dxe/dxw)*(potent[GRADYFI][iP]-potent[GRADYFI][iW]);
			potent[GRADZFI][iE]=potent[GRADZFI][iP]+(dxe/dxw)*(potent[GRADZFI][iP]-potent[GRADZFI][iW]);
		}

		if (bW) {
			potent[GRADXFI][iW]=potent[GRADXFI][iP]+(dxw/dxe)*(potent[GRADXFI][iP]-potent[GRADXFI][iE]);
			potent[GRADYFI][iW]=potent[GRADYFI][iP]+(dxw/dxe)*(potent[GRADYFI][iP]-potent[GRADYFI][iE]);
			potent[GRADZFI][iW]=potent[GRADZFI][iP]+(dxw/dxe)*(potent[GRADZFI][iP]-potent[GRADZFI][iE]);
		}

		if (bN) {
			potent[GRADXFI][iN]=potent[GRADXFI][iP]+(dyn/dys)*(potent[GRADXFI][iP]-potent[GRADXFI][iS]);
			potent[GRADYFI][iN]=potent[GRADYFI][iP]+(dyn/dys)*(potent[GRADYFI][iP]-potent[GRADYFI][iS]);
			potent[GRADZFI][iN]=potent[GRADZFI][iP]+(dyn/dys)*(potent[GRADZFI][iP]-potent[GRADZFI][iS]);
		}

		if (bS) {
			potent[GRADXFI][iS]=potent[GRADXFI][iP]+(dys/dyn)*(potent[GRADXFI][iP]-potent[GRADXFI][iN]);
			potent[GRADYFI][iS]=potent[GRADYFI][iP]+(dys/dyn)*(potent[GRADYFI][iP]-potent[GRADYFI][iN]);
			potent[GRADZFI][iS]=potent[GRADZFI][iP]+(dys/dyn)*(potent[GRADZFI][iP]-potent[GRADZFI][iN]);
		}

		if (bT) {
			potent[GRADXFI][iT]=potent[GRADXFI][iP]+(dzt/dzb)*(potent[GRADXFI][iP]-potent[GRADXFI][iB]);
			potent[GRADYFI][iT]=potent[GRADYFI][iP]+(dzt/dzb)*(potent[GRADYFI][iP]-potent[GRADYFI][iB]);
			potent[GRADZFI][iT]=potent[GRADZFI][iP]+(dzt/dzb)*(potent[GRADZFI][iP]-potent[GRADZFI][iB]);
		}

		if (bB) {
			potent[GRADXFI][iB]=potent[GRADXFI][iP]+(dzb/dzt)*(potent[GRADXFI][iP]-potent[GRADXFI][iT]);
			potent[GRADYFI][iB]=potent[GRADYFI][iP]+(dzb/dzt)*(potent[GRADYFI][iP]-potent[GRADYFI][iT]);
			potent[GRADZFI][iB]=potent[GRADZFI][iP]+(dzb/dzt)*(potent[GRADZFI][iP]-potent[GRADZFI][iT]);
		}
	}

} // green_gauss_FI_union

// учёт граничных условий для вычисления расстояния до ближайшей
// твёрдой неподвижной поверхности.
void my_nested_dissection_bound(integer inumber, integer maxelm, 
							  bool bDirichlet, BOUND* border_neighbor, integer ls, integer lw,
							  WALL* w, equation3D_bon* &slb, doublereal dbeta,
							  TOCHKA* pa, int** nvtx, doublereal* potent, integer &iMCB
							  ) 
{

	// Только на одной входной границе стоит условие Дирихле на остальных нет.

	 // potent - потенциал FI на основе градиентов которого 
	 // вычисляется кратчайшее расстояние до одной из входных границ.
	 // В данном методе значения potent могут использоваться для граничного условия
	 // Неймана повышенного порядка точности при соответствующем значении dbeta.

     // bDirichlet   осуществляется сборка только граничных условий Дирихле.
     // bDirichlet == false осуществляется сборка только однородных условий Неймана.

     // inumber - номер граничного КО.
	 // inumber изменяется от 0..maxbound-1

     bool bDirichleti=false; // bDirichleti,  i - internal


     // Сначала запишем граничные условия Дирихле
	 if ((border_neighbor[inumber].MCB<(ls+lw)) && (border_neighbor[inumber].MCB>=ls) && (!w[border_neighbor[inumber].MCB-ls].bsymmetry) && (!w[border_neighbor[inumber].MCB-ls].bpressure)) {
		// граничное условие Дирихле
		// Задана скорость на границе
        // Это не граница симметрии и не выходная граница.

		 doublereal vel_mag=sqrt(w[border_neighbor[inumber].MCB-ls].Vx*w[border_neighbor[inumber].MCB-ls].Vx+
			          w[border_neighbor[inumber].MCB-ls].Vy*w[border_neighbor[inumber].MCB-ls].Vy+
					  w[border_neighbor[inumber].MCB-ls].Vz*w[border_neighbor[inumber].MCB-ls].Vz);

		 doublereal epsilon0=1e-32;

		 

		 if (fabs(vel_mag)<epsilon0) {
			 // скорость на стенке равна нулю.
			 bDirichleti=false;
		 } else {
			 if ((iMCB==-1) || (iMCB==border_neighbor[inumber].MCB))
			 {
				 if (iMCB==-1) iMCB=border_neighbor[inumber].MCB;
			     bDirichleti=true;
			 }
		 }
        
	}
	

	 if (bDirichlet&&bDirichleti) {
			 // условие Дирихле функция равна нулю.
			 slb[inumber].aw=1.0;
		     slb[inumber].ai=0.0;
			 slb[inumber].b=0.0; // на твёрдой неподвижной стенке функция принимает значение 0.
			 slb[inumber].iI=-1; // не присутствует в матрице
		     slb[inumber].iW=border_neighbor[inumber].iB;
	}
	if ((!bDirichlet)&&(!bDirichleti)) {
		// однородное условие Неймана.
		doublereal dl, dS, deltal;// , fiplus;

		switch (border_neighbor[inumber].Norm) {
		case E_SIDE: dl=pa[nvtx[1][border_neighbor[inumber].iI]-1].x-pa[nvtx[0][border_neighbor[inumber].iI]-1].x;
                 dS=pa[nvtx[2][border_neighbor[inumber].iI]-1].y-pa[nvtx[1][border_neighbor[inumber].iI]-1].y; 
				 dS*=(pa[nvtx[4][border_neighbor[inumber].iI]-1].z-pa[nvtx[0][border_neighbor[inumber].iI]-1].z); // площадь грани
                 slb[inumber].ai=2.0*dbeta*dS/dl;
				 slb[inumber].iI=border_neighbor[inumber].iI;
				 slb[inumber].aw=slb[inumber].ai;
				 slb[inumber].iW=border_neighbor[inumber].iB;
                 deltal=0.5*(pa[nvtx[1][border_neighbor[inumber].iII]-1].x+pa[nvtx[0][border_neighbor[inumber].iII]-1].x);
				 deltal-=0.5*(pa[nvtx[1][border_neighbor[inumber].iI]-1].x+pa[nvtx[0][border_neighbor[inumber].iI]-1].x);
                 //fiplus=0.5*dl/deltal;
				 // правая часть
				 slb[inumber].b=(dbeta-1.0)*dS*(potent[border_neighbor[inumber].iI]-potent[border_neighbor[inumber].iII])/deltal;
			    break;
		case N_SIDE: 
			     dl = pa[nvtx[2][border_neighbor[inumber].iI]-1].y-pa[nvtx[0][border_neighbor[inumber].iI]-1].y;
                 dS=pa[nvtx[1][border_neighbor[inumber].iI]-1].x-pa[nvtx[0][border_neighbor[inumber].iI]-1].x; 
				 dS*=(pa[nvtx[4][border_neighbor[inumber].iI]-1].z-pa[nvtx[0][border_neighbor[inumber].iI]-1].z); // площадь грани
				 slb[inumber].ai=2.0*dbeta*dS/dl;
				 slb[inumber].iI=border_neighbor[inumber].iI;
				 slb[inumber].aw=slb[inumber].ai;
				 slb[inumber].iW=border_neighbor[inumber].iB;
				 deltal=0.5*(pa[nvtx[2][border_neighbor[inumber].iII]-1].y+pa[nvtx[0][border_neighbor[inumber].iII]-1].y);
				 deltal-=0.5*(pa[nvtx[2][border_neighbor[inumber].iI]-1].y+pa[nvtx[0][border_neighbor[inumber].iI]-1].y);
                 //fiplus=0.5*dl/deltal;
				 // правая часть
				 slb[inumber].b=(dbeta-1.0)*dS*(potent[border_neighbor[inumber].iI]-potent[border_neighbor[inumber].iII])/deltal;
				 break;
      case T_SIDE:  
			     dl = pa[nvtx[4][border_neighbor[inumber].iI]-1].z-pa[nvtx[0][border_neighbor[inumber].iI]-1].z;
                 dS=pa[nvtx[1][border_neighbor[inumber].iI]-1].x-pa[nvtx[0][border_neighbor[inumber].iI]-1].x; 
				 dS*=(pa[nvtx[2][border_neighbor[inumber].iI]-1].y-pa[nvtx[0][border_neighbor[inumber].iI]-1].y); // площадь грани
				 slb[inumber].ai=2.0*dbeta*dS/dl;
				 slb[inumber].iI=border_neighbor[inumber].iI;
				 slb[inumber].aw=slb[inumber].ai;
				 slb[inumber].iW=border_neighbor[inumber].iB;
				 deltal=0.5*(pa[nvtx[4][border_neighbor[inumber].iII]-1].z+pa[nvtx[0][border_neighbor[inumber].iII]-1].z);
				 deltal-=0.5*(pa[nvtx[4][border_neighbor[inumber].iI]-1].z+pa[nvtx[0][border_neighbor[inumber].iI]-1].z);
                 //fiplus=0.5*dl/deltal;
				 // правая часть
				 slb[inumber].b=(dbeta-1.0)*dS*(potent[border_neighbor[inumber].iI]-potent[border_neighbor[inumber].iII])/deltal;
                 break;
		case W_SIDE: 
			     dl = pa[nvtx[1][border_neighbor[inumber].iI]-1].x-pa[nvtx[0][border_neighbor[inumber].iI]-1].x;
                 dS=pa[nvtx[2][border_neighbor[inumber].iI]-1].y-pa[nvtx[1][border_neighbor[inumber].iI]-1].y; 
				 dS*=(pa[nvtx[4][border_neighbor[inumber].iI]-1].z-pa[nvtx[0][border_neighbor[inumber].iI]-1].z); // площадь грани
    			 slb[inumber].ai=2.0*dbeta*dS/dl;
				 slb[inumber].iI=border_neighbor[inumber].iI;
				 slb[inumber].aw=slb[inumber].ai;
				 slb[inumber].iW=border_neighbor[inumber].iB;
				 deltal=-0.5*(pa[nvtx[1][border_neighbor[inumber].iII]-1].x+pa[nvtx[0][border_neighbor[inumber].iII]-1].x);
				 deltal+=0.5*(pa[nvtx[1][border_neighbor[inumber].iI]-1].x+pa[nvtx[0][border_neighbor[inumber].iI]-1].x);
                 //fiplus=0.5*dl/deltal;
     			 // правая часть
				 slb[inumber].b=(dbeta-1.0)*dS*(potent[border_neighbor[inumber].iI]-potent[border_neighbor[inumber].iII])/deltal;
    			 break;
         case S_SIDE:
			     dl = pa[nvtx[2][border_neighbor[inumber].iI]-1].y-pa[nvtx[0][border_neighbor[inumber].iI]-1].y;
                 dS=pa[nvtx[1][border_neighbor[inumber].iI]-1].x-pa[nvtx[0][border_neighbor[inumber].iI]-1].x; 
				 dS*=(pa[nvtx[4][border_neighbor[inumber].iI]-1].z-pa[nvtx[0][border_neighbor[inumber].iI]-1].z); // площадь грани
				 slb[inumber].ai=2.0*dbeta*dS/dl;
				 slb[inumber].iI=border_neighbor[inumber].iI;
				 slb[inumber].aw=slb[inumber].ai;
				 slb[inumber].iW=border_neighbor[inumber].iB;
				 deltal=-0.5*(pa[nvtx[2][border_neighbor[inumber].iII]-1].y+pa[nvtx[0][border_neighbor[inumber].iII]-1].y);
				 deltal+=0.5*(pa[nvtx[2][border_neighbor[inumber].iI]-1].y+pa[nvtx[0][border_neighbor[inumber].iI]-1].y);
                 //fiplus=0.5*dl/deltal;
				 // правая часть
				 slb[inumber].b=(dbeta-1.0)*dS*(potent[border_neighbor[inumber].iI]-potent[border_neighbor[inumber].iII])/deltal;
				 break;
		 case B_SIDE: 
			     dl = pa[nvtx[4][border_neighbor[inumber].iI]-1].z-pa[nvtx[0][border_neighbor[inumber].iI]-1].z;
                 dS=pa[nvtx[1][border_neighbor[inumber].iI]-1].x-pa[nvtx[0][border_neighbor[inumber].iI]-1].x; 
				 dS*=(pa[nvtx[2][border_neighbor[inumber].iI]-1].y-pa[nvtx[0][border_neighbor[inumber].iI]-1].y); // площадь грани
				 slb[inumber].ai=2.0*dbeta*dS/dl;
				 slb[inumber].iI=border_neighbor[inumber].iI;
				 slb[inumber].aw=slb[inumber].ai;
				 slb[inumber].iW=border_neighbor[inumber].iB;
                 deltal=-0.5*(pa[nvtx[4][border_neighbor[inumber].iII]-1].z+pa[nvtx[0][border_neighbor[inumber].iII]-1].z);
    			 deltal+=0.5*(pa[nvtx[4][border_neighbor[inumber].iI]-1].z+pa[nvtx[0][border_neighbor[inumber].iI]-1].z);
                 //fiplus=0.5*dl/deltal;
				 // правая часть
				 slb[inumber].b=(dbeta-1.0)*dS*(potent[border_neighbor[inumber].iI]-potent[border_neighbor[inumber].iII])/deltal;
				break;
		} // end switch
	}

}

// учёт граничных условий для вычисления расстояния до ближайшей
// твёрдой неподвижной поверхности.
// Собирается в предположении всюду единичной сетки.
// 8 августа 2015.
void my_nested_dissection_bound_union(integer inumber, integer maxelm, 
							  bool bDirichlet, BOUND* border_neighbor, integer ls, integer lw,
							  WALL* w, equation3D_bon* &slb, doublereal dbeta,
							  TOCHKA* pa, int** nvtx, doublereal* potent, integer &iMCB
							  ) 
{

	// Только на одной входной границе стоит условие Дирихле на остальных нет.

	 // potent - потенциал FI на основе градиентов которого 
	 // вычисляется кратчайшее расстояние до одной из входных границ.
	 // В данном методе значения potent могут использоваться для граничного условия
	 // Неймана повышенного порядка точности при соответствующем значении dbeta.

     // bDirichlet   осуществляется сборка только граничных условий Дирихле.
     // bDirichlet == false осуществляется сборка только однородных условий Неймана.

     // inumber - номер граничного КО.
	 // inumber изменяется от 0..maxbound-1

     bool bDirichleti=false; // bDirichleti,  i - internal


     // Сначала запишем граничные условия Дирихле
	 if ((border_neighbor[inumber].MCB<(ls+lw)) && (border_neighbor[inumber].MCB>=ls) && (!w[border_neighbor[inumber].MCB-ls].bsymmetry) && (!w[border_neighbor[inumber].MCB-ls].bpressure)) {
		// граничное условие Дирихле
		// Задана скорость на границе
        // Это не граница симметрии и не выходная граница.

		 doublereal vel_mag=sqrt(w[border_neighbor[inumber].MCB-ls].Vx*w[border_neighbor[inumber].MCB-ls].Vx+
			          w[border_neighbor[inumber].MCB-ls].Vy*w[border_neighbor[inumber].MCB-ls].Vy+
					  w[border_neighbor[inumber].MCB-ls].Vz*w[border_neighbor[inumber].MCB-ls].Vz);

		 doublereal epsilon0=1e-32;

		 

		 if (fabs(vel_mag)<epsilon0) {
			 // скорость на стенке равна нулю.
			 bDirichleti=false;
		 } else {
			 if ((iMCB==-1) || (iMCB==border_neighbor[inumber].MCB))
			 {
				 if (iMCB==-1) iMCB=border_neighbor[inumber].MCB;
			     bDirichleti=true;
			 }
		 }
        
	}
	

	 if (bDirichlet&&bDirichleti) {
			 // условие Дирихле функция равна нулю.
			 slb[inumber].aw=1.0;
		     slb[inumber].ai=0.0;
			 slb[inumber].b=0.0; // на твёрдой неподвижной стенке функция принимает значение 0.
			 slb[inumber].iI=-1; // не присутствует в матрице
		     slb[inumber].iW=border_neighbor[inumber].iB;
	}
	if ((!bDirichlet)&&(!bDirichleti)) {
		// однородное условие Неймана.
		//doublereal dl, dS, deltal, fiplus;

		switch (border_neighbor[inumber].Norm) {
		case E_SIDE: //dl=pa[nvtx[1][border_neighbor[inumber].iI]-1].x-pa[nvtx[0][border_neighbor[inumber].iI]-1].x;
                 //dS=pa[nvtx[2][border_neighbor[inumber].iI]-1].y-pa[nvtx[1][border_neighbor[inumber].iI]-1].y; 
				 //dS*=(pa[nvtx[4][border_neighbor[inumber].iI]-1].z-pa[nvtx[0][border_neighbor[inumber].iI]-1].z); // площадь грани
                 slb[inumber].ai=1.0;//2.0*dbeta*dS/dl;
				 slb[inumber].iI=border_neighbor[inumber].iI;
				 slb[inumber].aw=1.0;//slb[inumber].ai;
				 slb[inumber].iW=border_neighbor[inumber].iB;
                // deltal=0.5*(pa[nvtx[1][border_neighbor[inumber].iII]-1].x+pa[nvtx[0][border_neighbor[inumber].iII]-1].x);
				// deltal-=0.5*(pa[nvtx[1][border_neighbor[inumber].iI]-1].x+pa[nvtx[0][border_neighbor[inumber].iI]-1].x);
                 //fiplus=0.5*dl/deltal;
				 // правая часть
				 slb[inumber].b=0.0;//(dbeta-1.0)*dS*(potent[border_neighbor[inumber].iI]-potent[border_neighbor[inumber].iII])/deltal;
			    break;
		case N_SIDE: 
			     //dl = pa[nvtx[2][border_neighbor[inumber].iI]-1].y-pa[nvtx[0][border_neighbor[inumber].iI]-1].y;
                 //dS=pa[nvtx[1][border_neighbor[inumber].iI]-1].x-pa[nvtx[0][border_neighbor[inumber].iI]-1].x; 
				 //dS*=(pa[nvtx[4][border_neighbor[inumber].iI]-1].z-pa[nvtx[0][border_neighbor[inumber].iI]-1].z); // площадь грани
				 slb[inumber].ai=1.0;//2.0*dbeta*dS/dl;
				 slb[inumber].iI=border_neighbor[inumber].iI;
				 slb[inumber].aw=1.0;//slb[inumber].ai;
				 slb[inumber].iW=border_neighbor[inumber].iB;
				 //deltal=0.5*(pa[nvtx[2][border_neighbor[inumber].iII]-1].y+pa[nvtx[0][border_neighbor[inumber].iII]-1].y);
				 //deltal-=0.5*(pa[nvtx[2][border_neighbor[inumber].iI]-1].y+pa[nvtx[0][border_neighbor[inumber].iI]-1].y);
                 //fiplus=0.5*dl/deltal;
				 // правая часть
				 slb[inumber].b=0.0;//(dbeta-1.0)*dS*(potent[border_neighbor[inumber].iI]-potent[border_neighbor[inumber].iII])/deltal;
				 break;
      case T_SIDE:  
			     //dl = pa[nvtx[4][border_neighbor[inumber].iI]-1].z-pa[nvtx[0][border_neighbor[inumber].iI]-1].z;
                 //dS=pa[nvtx[1][border_neighbor[inumber].iI]-1].x-pa[nvtx[0][border_neighbor[inumber].iI]-1].x; 
				 //dS*=(pa[nvtx[2][border_neighbor[inumber].iI]-1].y-pa[nvtx[0][border_neighbor[inumber].iI]-1].y); // площадь грани
				 slb[inumber].ai=1.0;//2.0*dbeta*dS/dl;
				 slb[inumber].iI=border_neighbor[inumber].iI;
				 slb[inumber].aw=1.0;//slb[inumber].ai;
				 slb[inumber].iW=border_neighbor[inumber].iB;
				 //deltal=0.5*(pa[nvtx[4][border_neighbor[inumber].iII]-1].z+pa[nvtx[0][border_neighbor[inumber].iII]-1].z);
				 //deltal-=0.5*(pa[nvtx[4][border_neighbor[inumber].iI]-1].z+pa[nvtx[0][border_neighbor[inumber].iI]-1].z);
                 //fiplus=0.5*dl/deltal;
				 // правая часть
				 slb[inumber].b=0.0;//(dbeta-1.0)*dS*(potent[border_neighbor[inumber].iI]-potent[border_neighbor[inumber].iII])/deltal;
                 break;
		case W_SIDE: 
			     //dl = pa[nvtx[1][border_neighbor[inumber].iI]-1].x-pa[nvtx[0][border_neighbor[inumber].iI]-1].x;
                 //dS=pa[nvtx[2][border_neighbor[inumber].iI]-1].y-pa[nvtx[1][border_neighbor[inumber].iI]-1].y; 
				 //dS*=(pa[nvtx[4][border_neighbor[inumber].iI]-1].z-pa[nvtx[0][border_neighbor[inumber].iI]-1].z); // площадь грани
    			 slb[inumber].ai=1.0;//2.0*dbeta*dS/dl;
				 slb[inumber].iI=border_neighbor[inumber].iI;
				 slb[inumber].aw=1.0;//slb[inumber].ai;
				 slb[inumber].iW=border_neighbor[inumber].iB;
				 //deltal=-0.5*(pa[nvtx[1][border_neighbor[inumber].iII]-1].x+pa[nvtx[0][border_neighbor[inumber].iII]-1].x);
				 //deltal+=0.5*(pa[nvtx[1][border_neighbor[inumber].iI]-1].x+pa[nvtx[0][border_neighbor[inumber].iI]-1].x);
                 //fiplus=0.5*dl/deltal;
     			 // правая часть
				 slb[inumber].b=0.0;//(dbeta-1.0)*dS*(potent[border_neighbor[inumber].iI]-potent[border_neighbor[inumber].iII])/deltal;
    			 break;
         case S_SIDE:
			     //dl = pa[nvtx[2][border_neighbor[inumber].iI]-1].y-pa[nvtx[0][border_neighbor[inumber].iI]-1].y;
                 //dS=pa[nvtx[1][border_neighbor[inumber].iI]-1].x-pa[nvtx[0][border_neighbor[inumber].iI]-1].x; 
				 //dS*=(pa[nvtx[4][border_neighbor[inumber].iI]-1].z-pa[nvtx[0][border_neighbor[inumber].iI]-1].z); // площадь грани
				 slb[inumber].ai=1.0;//2.0*dbeta*dS/dl;
				 slb[inumber].iI=border_neighbor[inumber].iI;
				 slb[inumber].aw=1.0;//slb[inumber].ai;
				 slb[inumber].iW=border_neighbor[inumber].iB;
				 //deltal=-0.5*(pa[nvtx[2][border_neighbor[inumber].iII]-1].y+pa[nvtx[0][border_neighbor[inumber].iII]-1].y);
				 //deltal+=0.5*(pa[nvtx[2][border_neighbor[inumber].iI]-1].y+pa[nvtx[0][border_neighbor[inumber].iI]-1].y);
                 //fiplus=0.5*dl/deltal;
				 // правая часть
				 slb[inumber].b=0.0;//(dbeta-1.0)*dS*(potent[border_neighbor[inumber].iI]-potent[border_neighbor[inumber].iII])/deltal;
				 break;
		 case B_SIDE: 
			     //dl = pa[nvtx[4][border_neighbor[inumber].iI]-1].z-pa[nvtx[0][border_neighbor[inumber].iI]-1].z;
                 //dS=pa[nvtx[1][border_neighbor[inumber].iI]-1].x-pa[nvtx[0][border_neighbor[inumber].iI]-1].x; 
				 //dS*=(pa[nvtx[2][border_neighbor[inumber].iI]-1].y-pa[nvtx[0][border_neighbor[inumber].iI]-1].y); // площадь грани
				 slb[inumber].ai=1.0;//2.0*dbeta*dS/dl;
				 slb[inumber].iI=border_neighbor[inumber].iI;
				 slb[inumber].aw=1.0;//slb[inumber].ai;
				 slb[inumber].iW=border_neighbor[inumber].iB;
                 //deltal=-0.5*(pa[nvtx[4][border_neighbor[inumber].iII]-1].z+pa[nvtx[0][border_neighbor[inumber].iII]-1].z);
    			 //deltal+=0.5*(pa[nvtx[4][border_neighbor[inumber].iI]-1].z+pa[nvtx[0][border_neighbor[inumber].iI]-1].z);
                 //fiplus=0.5*dl/deltal;
				 // правая часть
				 slb[inumber].b=0.0;//(dbeta-1.0)*dS*(potent[border_neighbor[inumber].iI]-potent[border_neighbor[inumber].iII])/deltal;
				break;
		} // end switch
	}

} // my_nested_dissection_bound_union


// Составляет матрицу для уравнения 
// поправки давления
void my_elmatr_quad_nested_dissection(int iP, equation3D* &sl, equation3D_bon* &slb,  
	TOCHKA* pa, int** nvtx, int*** neighbors_for_the_internal_node, integer maxelm, doublereal dbeta) {

   
	doublereal eps=1e-37; // для отделения вещественного нуля.
	

	// Внутренний узел и его соседи:

    // iP - номер центрального контрольного объёма
	int iE, iN, iT, iW, iS, iB; // номера соседних контрольных объёмов
	iE=neighbors_for_the_internal_node[E_SIDE][0][iP]; iN=neighbors_for_the_internal_node[N_SIDE][0][iP]; iT=neighbors_for_the_internal_node[T_SIDE][0][iP];
	iW=neighbors_for_the_internal_node[W_SIDE][0][iP]; iS=neighbors_for_the_internal_node[S_SIDE][0][iP]; iB=neighbors_for_the_internal_node[B_SIDE][0][iP];
	sl[iP].iP=iP;
	sl[iP].iE=iE; sl[iP].iN=iN; 
	sl[iP].iS=iS; sl[iP].iW=iW;
    sl[iP].iT=iT; sl[iP].iB=iB;


	// если с одной из сторон граница расчётной области 
	// то переменная равна true
	bool bE=false, bN=false, bT=false, bW=false, bS=false, bB=false;

    if (iE>=maxelm) bE=true;
	if (iN>=maxelm) bN=true;
	if (iT>=maxelm) bT=true;
    if (iW>=maxelm) bW=true;
	if (iS>=maxelm) bS=true;
	if (iB>=maxelm) bB=true;

	// вычисление размеров текущего контрольного объёма:
	doublereal dx=0.0, dy=0.0, dz=0.0; // размеры контрольного объёма
    volume3D(iP, nvtx, pa, dx, dy, dz);
	    

    doublereal dxe=0.5*dx, dxw=0.5*dx, dyn=0.5*dy, dys=0.5*dy, dzt=0.5*dz, dzb=0.5*dz;
    // т.к. известна нумерация вершин куба, то здесь она используется
	// x - direction
    if (!bE) dxe=0.5*(pa[nvtx[1][iE]-1].x+pa[nvtx[0][iE]-1].x);
	if (!bE) dxe-=0.5*(pa[nvtx[1][iP]-1].x+pa[nvtx[0][iP]-1].x);
	if (!bW) dxw=0.5*(pa[nvtx[1][iP]-1].x+pa[nvtx[0][iP]-1].x);
	if (!bW) dxw-=0.5*(pa[nvtx[1][iW]-1].x+pa[nvtx[0][iW]-1].x);
    // y - direction
	if (!bN) dyn=0.5*(pa[nvtx[2][iN]-1].y+pa[nvtx[0][iN]-1].y);
	if (!bN) dyn-=0.5*(pa[nvtx[2][iP]-1].y+pa[nvtx[0][iP]-1].y);
	if (!bS) dys=0.5*(pa[nvtx[2][iP]-1].y+pa[nvtx[0][iP]-1].y);
	if (!bS) dys-=0.5*(pa[nvtx[2][iS]-1].y+pa[nvtx[0][iS]-1].y);
    // z - direction
	if (!bT) dzt=0.5*(pa[nvtx[4][iT]-1].z+pa[nvtx[0][iT]-1].z);
	if (!bT) dzt-=0.5*(pa[nvtx[4][iP]-1].z+pa[nvtx[0][iP]-1].z);
	if (!bB) dzb=0.5*(pa[nvtx[4][iP]-1].z+pa[nvtx[0][iP]-1].z);
	if (!bB) dzb-=0.5*(pa[nvtx[4][iB]-1].z+pa[nvtx[0][iB]-1].z);


    doublereal De, Dw, Ds, Dn, Dt, Db; // диффузионный поток через грани КО.
	

	if (!bE) {
		if (bW) De=dbeta*dy*dz/dxe;
		else De=dy*dz/dxe;
	}
	else De=dbeta*dy*dz/dxe;

	if (!bW) {
		if (bE) Dw=dbeta*dy*dz/dxw;
		else Dw=dy*dz/dxw; 
	}
	else Dw=dbeta*dy*dz/dxw;

	if (!bN) {
		if (bS) Dn=dbeta*dx*dz/dyn;
		else Dn=dx*dz/dyn; 
	}
	else Dn=dbeta*dx*dz/dyn;


	if (!bS) {
		if (bN) Ds=dbeta*dx*dz/dys;
		else Ds=dx*dz/dys;
	}
	else Ds=dbeta*dx*dz/dys;


	if (!bT) {
		if (bB) Dt=dbeta*dx*dy/dzt;
		else Dt=dx*dy/dzt;
	} 
	else Dt=dbeta*dx*dy/dzt;


	if (!bB) {
		if (bT) Db=dbeta*dx*dy/dzb;
		else Db=dx*dy/dzb;
	}
	else Db=dbeta*dx*dy/dzb;

	// Число Пекле равно нулю т.к. это чисто эллиптическое уравнение.
	sl[iP].ae=De*1.0; // при Pe==0.0 величина fD(0.0, EXP2, true, feplus); равна строго 1.0;
	sl[iP].aw=Dw; // *fD(0.0, EXP2, true, fwplus); равна строго 1.0;
	sl[iP].an=Dn; // *fD(0.0, EXP2, true, fnplus); равна строго 1.0;
	sl[iP].as=Ds; // *fD(0.0, EXP2, true, fsplus); равна строго 1.0;
	sl[iP].at=Dt; // *fD(0.0, EXP2, true, ftplus); равна строго 1.0;
	sl[iP].ab=Db; // *fD(0.0, EXP2, true, fbplus); равна строго 1.0;

	sl[iP].ap=sl[iP].ae+sl[iP].aw+sl[iP].an+sl[iP].as+sl[iP].at+sl[iP].ab;

	doublereal dSc=-1.0;

	sl[iP].b=dSc*dx*dy*dz;

	// Симметризация СЛАУ:

	// для потенциала FI получается эллиптическое уравнение с SPD матрицей.
	
	// Строка матрицы выглядит примерно следующим образом:
	// -ab ... -as ... -aw ... +ap ... -ae ... -an ... -at == b

	    // 1. Учёт условия Дирихле:
		if ((iE>=maxelm) && (fabs(slb[iE-maxelm].ai)<eps)) {
			sl[iP].b+=sl[iP].ae*slb[iE-maxelm].b/slb[iE-maxelm].aw;
			sl[iP].ae=0.0;
			sl[iP].iE=-1; // не входит в матрицу СЛАУ.
		}
		if ((iW>=maxelm) && (fabs(slb[iW-maxelm].ai)<eps)) {
			sl[iP].b+=sl[iP].aw*slb[iW-maxelm].b/slb[iW-maxelm].aw;
			sl[iP].aw=0.0;
			sl[iP].iW=-1; // не входит в матрицу СЛАУ.
		}
		if ((iN>=maxelm) && (fabs(slb[iN-maxelm].ai)<eps)) {
			sl[iP].b+=sl[iP].an*slb[iN-maxelm].b/slb[iN-maxelm].aw;
			sl[iP].an=0.0;
			sl[iP].iN=-1; // не входит в матрицу СЛАУ.
		}
		if ((iS>=maxelm) && (fabs(slb[iS-maxelm].ai)<eps)) {
			sl[iP].b+=sl[iP].as*slb[iS-maxelm].b/slb[iS-maxelm].aw;
			sl[iP].as=0.0;
			sl[iP].iS=-1; // не входит в матрицу СЛАУ.
		}
		if ((iT>=maxelm) && (fabs(slb[iT-maxelm].ai)<eps)) {
			sl[iP].b+=sl[iP].at*slb[iT-maxelm].b/slb[iT-maxelm].aw;
			sl[iP].at=0.0;
			sl[iP].iT=-1; // не входит в матрицу СЛАУ.
		}
		if ((iB>=maxelm) && (fabs(slb[iB-maxelm].ai)<eps)) {
			sl[iP].b+=sl[iP].ab*slb[iB-maxelm].b/slb[iB-maxelm].aw;
			sl[iP].ab=0.0;
			sl[iP].iB=-1; // не входит в матрицу СЛАУ.
		}

} // внутренность матрицы для вычисления кратчайшего расстояния до стенки.


// Составляет матрицу для уравнения 
// nested desection: Т.к. для nested desection важно лишь число контрольных объёмов но не их 
// размеры то матрица собирается как будто сетка состоит из единичных кубов все стороны которой имеют 
// длину единица. Сборка матрицы для единичной такой равномерной сетке должна очень плодотворно повлиять 
// на качество разбиения nested desection.
// 8 августа 2015. 
void my_elmatr_quad_nested_dissection_union(int iP, equation3D* &sl, equation3D_bon* &slb,  
	TOCHKA* pa, int** nvtx, int*** neighbors_for_the_internal_node, integer maxelm, doublereal dbeta) {

   
	doublereal eps=1e-37; // для отделения вещественного нуля.
	

	// Внутренний узел и его соседи:

    // iP - номер центрального контрольного объёма
	int iE, iN, iT, iW, iS, iB; // номера соседних контрольных объёмов
	iE=neighbors_for_the_internal_node[E_SIDE][0][iP]; iN=neighbors_for_the_internal_node[N_SIDE][0][iP]; iT=neighbors_for_the_internal_node[T_SIDE][0][iP];
	iW=neighbors_for_the_internal_node[W_SIDE][0][iP]; iS=neighbors_for_the_internal_node[S_SIDE][0][iP]; iB=neighbors_for_the_internal_node[B_SIDE][0][iP];
	sl[iP].iP=iP;
	sl[iP].iE=iE; sl[iP].iN=iN; 
	sl[iP].iS=iS; sl[iP].iW=iW;
    sl[iP].iT=iT; sl[iP].iB=iB;


	// если с одной из сторон граница расчётной области 
	// то переменная равна true
	//bool bE=false, bN=false, bT=false, bW=false, bS=false, bB=false;

    //if (iE>=maxelm) bE=true;
	//if (iN>=maxelm) bN=true;
	//if (iT>=maxelm) bT=true;
    //if (iW>=maxelm) bW=true;
	//if (iS>=maxelm) bS=true;
	//if (iB>=maxelm) bB=true;

	/*
	// вычисление размеров текущего контрольного объёма:
	doublereal dx=0.0, dy=0.0, dz=0.0; // размеры контрольного объёма
    volume3D(iP, nvtx, pa, dx, dy, dz);
	    

    doublereal dxe=0.5*dx, dxw=0.5*dx, dyn=0.5*dy, dys=0.5*dy, dzt=0.5*dz, dzb=0.5*dz;
    // т.к. известна нумерация вершин куба, то здесь она используется
	// x - direction
    if (!bE) dxe=0.5*(pa[nvtx[1][iE]-1].x+pa[nvtx[0][iE]-1].x);
	if (!bE) dxe-=0.5*(pa[nvtx[1][iP]-1].x+pa[nvtx[0][iP]-1].x);
	if (!bW) dxw=0.5*(pa[nvtx[1][iP]-1].x+pa[nvtx[0][iP]-1].x);
	if (!bW) dxw-=0.5*(pa[nvtx[1][iW]-1].x+pa[nvtx[0][iW]-1].x);
    // y - direction
	if (!bN) dyn=0.5*(pa[nvtx[2][iN]-1].y+pa[nvtx[0][iN]-1].y);
	if (!bN) dyn-=0.5*(pa[nvtx[2][iP]-1].y+pa[nvtx[0][iP]-1].y);
	if (!bS) dys=0.5*(pa[nvtx[2][iP]-1].y+pa[nvtx[0][iP]-1].y);
	if (!bS) dys-=0.5*(pa[nvtx[2][iS]-1].y+pa[nvtx[0][iS]-1].y);
    // z - direction
	if (!bT) dzt=0.5*(pa[nvtx[4][iT]-1].z+pa[nvtx[0][iT]-1].z);
	if (!bT) dzt-=0.5*(pa[nvtx[4][iP]-1].z+pa[nvtx[0][iP]-1].z);
	if (!bB) dzb=0.5*(pa[nvtx[4][iP]-1].z+pa[nvtx[0][iP]-1].z);
	if (!bB) dzb-=0.5*(pa[nvtx[4][iB]-1].z+pa[nvtx[0][iB]-1].z);
	*/

    doublereal De, Dw, Ds, Dn, Dt, Db; // диффузионный поток через грани КО.
	
/*
	if (!bE) {
		if (bW) De=dbeta*dy*dz/dxe;
		else De=dy*dz/dxe;
	}
	else De=dbeta*dy*dz/dxe;

	if (!bW) {
		if (bE) Dw=dbeta*dy*dz/dxw;
		else Dw=dy*dz/dxw; 
	}
	else Dw=dbeta*dy*dz/dxw;

	if (!bN) {
		if (bS) Dn=dbeta*dx*dz/dyn;
		else Dn=dx*dz/dyn; 
	}
	else Dn=dbeta*dx*dz/dyn;


	if (!bS) {
		if (bN) Ds=dbeta*dx*dz/dys;
		else Ds=dx*dz/dys;
	}
	else Ds=dbeta*dx*dz/dys;


	if (!bT) {
		if (bB) Dt=dbeta*dx*dy/dzt;
		else Dt=dx*dy/dzt;
	} 
	else Dt=dbeta*dx*dy/dzt;


	if (!bB) {
		if (bT) Db=dbeta*dx*dy/dzb;
		else Db=dx*dy/dzb;
	}
	else Db=dbeta*dx*dy/dzb;

	*/
	De=1.0;
	Dw=1.0;
	Dn=1.0;
	Ds=1.0;
	Dt=1.0;
	Db=1.0;

	// Число Пекле равно нулю т.к. это чисто эллиптическое уравнение.
	sl[iP].ae=De*1.0; // при Pe==0.0 величина fD(0.0, EXP2, true, feplus); равна строго 1.0;
	sl[iP].aw=Dw; // *fD(0.0, EXP2, true, fwplus); равна строго 1.0;
	sl[iP].an=Dn; // *fD(0.0, EXP2, true, fnplus); равна строго 1.0;
	sl[iP].as=Ds; // *fD(0.0, EXP2, true, fsplus); равна строго 1.0;
	sl[iP].at=Dt; // *fD(0.0, EXP2, true, ftplus); равна строго 1.0;
	sl[iP].ab=Db; // *fD(0.0, EXP2, true, fbplus); равна строго 1.0;

	sl[iP].ap=sl[iP].ae+sl[iP].aw+sl[iP].an+sl[iP].as+sl[iP].at+sl[iP].ab;

	doublereal dSc=-1.0;

	//sl[iP].b=dSc*dx*dy*dz;
	sl[iP].b=dSc*1.0;

	// Симметризация СЛАУ:

	// для потенциала FI получается эллиптическое уравнение с SPD матрицей.
	
	// Строка матрицы выглядит примерно следующим образом:
	// -ab ... -as ... -aw ... +ap ... -ae ... -an ... -at == b

	    // 1. Учёт условия Дирихле:
		if ((iE>=maxelm) && (fabs(slb[iE-maxelm].ai)<eps)) {
			sl[iP].b+=sl[iP].ae*slb[iE-maxelm].b/slb[iE-maxelm].aw;
			sl[iP].ae=0.0;
			sl[iP].iE=-1; // не входит в матрицу СЛАУ.
		}
		if ((iW>=maxelm) && (fabs(slb[iW-maxelm].ai)<eps)) {
			sl[iP].b+=sl[iP].aw*slb[iW-maxelm].b/slb[iW-maxelm].aw;
			sl[iP].aw=0.0;
			sl[iP].iW=-1; // не входит в матрицу СЛАУ.
		}
		if ((iN>=maxelm) && (fabs(slb[iN-maxelm].ai)<eps)) {
			sl[iP].b+=sl[iP].an*slb[iN-maxelm].b/slb[iN-maxelm].aw;
			sl[iP].an=0.0;
			sl[iP].iN=-1; // не входит в матрицу СЛАУ.
		}
		if ((iS>=maxelm) && (fabs(slb[iS-maxelm].ai)<eps)) {
			sl[iP].b+=sl[iP].as*slb[iS-maxelm].b/slb[iS-maxelm].aw;
			sl[iP].as=0.0;
			sl[iP].iS=-1; // не входит в матрицу СЛАУ.
		}
		if ((iT>=maxelm) && (fabs(slb[iT-maxelm].ai)<eps)) {
			sl[iP].b+=sl[iP].at*slb[iT-maxelm].b/slb[iT-maxelm].aw;
			sl[iP].at=0.0;
			sl[iP].iT=-1; // не входит в матрицу СЛАУ.
		}
		if ((iB>=maxelm) && (fabs(slb[iB-maxelm].ai)<eps)) {
			sl[iP].b+=sl[iP].ab*slb[iB-maxelm].b/slb[iB-maxelm].aw;
			sl[iP].ab=0.0;
			sl[iP].iB=-1; // не входит в матрицу СЛАУ.
		}

} // внутренность матрицы для вычисления кратчайшего расстояния до стенки.



// проверка построенной сетки
// экспорт результата расчёта в программу tecplot360
// часть 2.
void exporttecplotxy360T_3D_part2nd(integer maxelm, integer ncell, FLOW* &f, TEMPER &t, integer flow_interior_count, integer ianimate,
									doublereal *fvisible)
{
    // ianimate - номер добавляемый к имени файла для анимации.
	bool bprintmessage=false;

	FILE *fp=nullptr;
    FILE *fp1=nullptr; // часть 1 или 3
	
#ifdef MINGW_COMPILLER
	int err = 0;
	fp=fopen64("ALICEFLOW0_07_temp.PLT", "w");
	if (fp == nullptr) err = 1;
#else
	errno_t err = 0;
	err = fopen_s(&fp, "ALICEFLOW0_07_temp.PLT", "w");
#endif
	// создание файла для записи:
	// файл состоит из трёх частей: 
	// 1 и 3 часть записываются сразу
	// вторая часть с результатами расчёта записывается
	// после расчёта. Такая трёхэтапная запись файла выбрана в целях
	// сокращения объёма используемой оперативной памяти.
	// Экономия памяти 19N.

	
	// чтение частей 1 и 3 и запись всех трёх частей в итоговый файл.
	// 
	if ((err) != 0) {
		printf("Create File temp Error\n");
		//system("pause");
		system("pause");

	}
	else {

		if (fp != nullptr) {

			// ИМЕННО INT.!!!
			int c; // читаемый символ
			integer ivarexport = 1; // по умолчанию только поле температур:
			integer i = 0; // счётчик цикла

#ifdef MINGW_COMPILLER
			err = 0;
			fp1=fopen64("ALICEFLOW0_06_temp_part1.txt", "r");
			if (fp1 == nullptr) err = 1;
#else
			err = fopen_s(&fp1, "ALICEFLOW0_06_temp_part1.txt", "r");
#endif

			if ((err) != 0) {
				printf("Open File temp part1 Error\n");
				//system("pause");
				system("pause");

			}
			else {

				if (fp1!=nullptr) { 

				// копирование первой части в итоговый файл
				// Особенность: иногда необходимо изменить вторую строку в файле:
				if (flow_interior_count > 0) {
					// есть жидкие зоны. Теперь нужно проверить активность жидких зон.
					for (i = 0; i < flow_interior_count; ++i) if (f[i].bactive) {
						ivarexport = 3; // считаем что температура всегда активна
					}
				}

				if (ivarexport == 1) {
					// запись заголовка
					fprintf(fp, "TITLE = \"ALICEFLOW0_07\"\n");

					// запись имён переменных
					fprintf(fp, "VARIABLES = x, y, z, Temp\n");  // печатается только поле температур

					// запись информации о зонах
#if doubleintprecision == 1
					fprintf(fp, "ZONE T=\"Rampant\", N=%lld, E=%lld, ET=BRICK, F=FEBLOCK\n\n", maxelm, ncell);
#else
					fprintf(fp, "ZONE T=\"Rampant\", N=%d, E=%d, ET=BRICK, F=FEBLOCK\n\n", maxelm, ncell);
#endif
					


					while ((c = fgetc(fp1)) != EOF) fputc(c, fp);
				}
				else if (ivarexport == 3) {
					// запись заголовка
					fprintf(fp, "TITLE = \"ALICEFLOW0_07\"\n");

					// Полный набор искомых величин и теплопередача и гидродинамика:
					fprintf(fp, "\nVARIABLES = x, y, z, Temp, data\n");

					// запись информации о зонах
#if doubleintprecision == 1
					fprintf(fp, "ZONE T=\"Rampant\", N=%lld, E=%lld, ET=BRICK, F=FEBLOCK\n\n", maxelm, ncell);
#else
					fprintf(fp, "ZONE T=\"Rampant\", N=%d, E=%d, ET=BRICK, F=FEBLOCK\n\n", maxelm, ncell);
#endif
					

					while ((c = fgetc(fp1)) != EOF) fputc(c, fp);
				}
				fclose(fp1); // закрытие файла
				if (bprintmessage) {
					printf("export tecplot part1 is successfully reading and written...OK.\n");
				}

			}
			}

			// запись второй части

			// Запись поля температур производится всегда.

			// запись температуры
			for (i = 0; i < maxelm; ++i) {
				fprintf(fp, "%+.16f ", t.potent[i]);
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			fprintf(fp, "\n");



			// Запись гидродинамических величин если необходимо:
			if (ivarexport == 3) {
				// Печатаем переданные данные.
				for (i = 0; i<maxelm; ++i) {
					if (t.ptr[1][i]>-1) {
						fprintf(fp, "%+.16f ", fvisible[t.ptr[0][i]]);
					}
					else fprintf(fp, "%+.16f ", 0.0);
					if (i % 10 == 0) fprintf(fp, "\n");
				}

				fprintf(fp, "\n");


			}
#ifdef MINGW_COMPILLER
			err = 0;
			fp1=fopen64("ALICEFLOW0_06_temp_part3.txt", "r");
			if (fp1 == nullptr) err = 1;
#else
			err = fopen_s(&fp1, "ALICEFLOW0_06_temp_part3.txt", "r");
#endif
			if ((err) != 0) {
				printf("Open File temp part3 Error\n");
				//system("pause");
				system("pause");

			}
			else {
				if (fp1 != nullptr) {

					// копирование третьей части в итоговый файл
					while ((c = fgetc(fp1)) != EOF) fputc(c, fp);
					fclose(fp1); // закрытие файла
					if (bprintmessage) {
						printf("export tecplot part1 is successfully reading and written...OK.\n");
					}
				}
			}

			fclose(fp); // закрытие файла
			if (bprintmessage) {
				printf("export tecplot is successfully written...OK.\n");
			}
			else printf("export tecplot 360... "); // короткое сообщение без перехода на новую строку.
		}
	}

	// WinExec("C:\\Program Files (x86)\\Tecplot\\Tec360 2008\\bin\\tec360.exe ALICEFLOW0_03.PLT",SW_NORMAL);
	  //WinExec("C:\\Program Files (x86)\\Tecplot\\Tec360 2009\\bin\\tec360.exe ALICEFLOW0_03.PLT", SW_NORMAL);

}

// используется в nested desection.
const integer xnorm=0;
const integer ynorm=1;
const integer znorm=2;

// Рекурсивный алгоритм заливки данной поверхности определяемой точкой iP и нормалью iNorm.
// помечает всю связную окрестность связную с заданной нормалью в окрестности точки iP.
// 0 - не помечено, 1 - помечено, 2 - запрет на пометку.
void fill_marker_surface(integer iP, integer* &color, integer iNorm, FLOW &f) {
	// 0 - не помечено, 1 - помечено, 2 - запрет на пометку.
	integer iE, iN, iT, iW, iS, iB; // номера соседних контрольных объёмов
	iE=f.neighbors_for_the_internal_node[E_SIDE][0][iP]; iN=f.neighbors_for_the_internal_node[N_SIDE][0][iP]; iT=f.neighbors_for_the_internal_node[T_SIDE][0][iP]; iW=f.neighbors_for_the_internal_node[W_SIDE][0][iP]; iS=f.neighbors_for_the_internal_node[S_SIDE][0][iP]; iB=f.neighbors_for_the_internal_node[B_SIDE][0][iP];

	switch (iNorm) {
	case xnorm: 
		if ((color[iN]!=2)&&(color[iN]!=1)) {
		    color[iN]=1; // помечено.
		    if (iN<f.maxelm) fill_marker_surface(iN, color, iNorm, f);
	    }
		if ((color[iS]!=2)&&(color[iS]!=1)) {
		    color[iS]=1; // помечено.
		    if (iS<f.maxelm) fill_marker_surface(iS, color, iNorm, f);
	    }
		if ((color[iT]!=2)&&(color[iT]!=1)) {
		    color[iT]=1; // помечено.
		    if (iT<f.maxelm) fill_marker_surface(iT, color, iNorm, f);
	    }
		if ((color[iB]!=2)&&(color[iB]!=1)) {
		    color[iB]=1; // помечено.
		    if (iB<f.maxelm) fill_marker_surface(iB, color, iNorm, f);
	    }
		// Также нужно пометить граничные узлы в направлениях перпендикулярных плоскости разделителя.
		if (iE>=f.maxelm) {
			if ((color[iE]!=2)&&(color[iE]!=1)) {
				color[iE]=1;
			}
		}
		if (iW>=f.maxelm) {
			if ((color[iW]!=2)&&(color[iW]!=1)) {
				color[iW]=1;
			}
		}
		break;
	case ynorm: 
		if ((color[iE]!=2)&&(color[iE]!=1)) {
		    color[iE]=1; // помечено.
		    if (iE<f.maxelm) fill_marker_surface(iE, color, iNorm, f);
	    }
		if ((color[iW]!=2)&&(color[iW]!=1)) {
		    color[iW]=1; // помечено.
		    if (iW<f.maxelm) fill_marker_surface(iW, color, iNorm, f);
	    }
		if ((color[iT]!=2)&&(color[iT]!=1)) {
		    color[iT]=1; // помечено.
		    if (iT<f.maxelm) fill_marker_surface(iT, color, iNorm, f);
	    }
		if ((color[iB]!=2)&&(color[iB]!=1)) {
		    color[iB]=1; // помечено.
		    if (iB<f.maxelm) fill_marker_surface(iB, color, iNorm, f);
	    }
		// Также нужно пометить граничные узлы в направлениях перпендикулярных плоскости разделителя.
		if (iN>=f.maxelm) {
			if ((color[iN]!=2)&&(color[iN]!=1)) {
				color[iN]=1;
			}
		}
		if (iS>=f.maxelm) {
			if ((color[iS]!=2)&&(color[iS]!=1)) {
				color[iS]=1;
			}
		}
		break;
	case znorm: 
		if ((color[iN]!=2)&&(color[iN]!=1)) {
		    color[iN]=1; // помечено.
		    if (iN<f.maxelm) fill_marker_surface(iN, color, iNorm, f);
	    }
		if ((color[iS]!=2)&&(color[iS]!=1)) {
		    color[iS]=1; // помечено.
		    if (iS<f.maxelm) fill_marker_surface(iS, color, iNorm, f);
	    }
		if ((color[iE]!=2)&&(color[iE]!=1)) {
		    color[iE]=1; // помечено.
		    if (iE<f.maxelm) fill_marker_surface(iE, color, iNorm, f);
	    }
		if ((color[iW]!=2)&&(color[iW]!=1)) {
		    color[iW]=1; // помечено.
		    if (iW<f.maxelm) fill_marker_surface(iW, color, iNorm, f);
	    }
		// Также нужно пометить граничные узлы в направлениях перпендикулярных плоскости разделителя.
		if (iT>=f.maxelm) {
			if ((color[iT]!=2)&&(color[iT]!=1)) {
				color[iT]=1;
			}
		}
		if (iB>=f.maxelm) {
			if ((color[iB]!=2)&&(color[iB]!=1)) {
				color[iB]=1;
			}
		}
		break;
	}
} // fill_marker_surface

// Раскрашивает связный объём.
// помечает всю связную окрестность связную с заданной нормалью в окрестности точки iP.
// 0 - не помечено, 1 - помечено, 2 - запрет на пометку.
// Внимание от рекурсивного алгоритма пришлось отказаться так он не подходит из-за большого количества
// рекурсивных вызовов.
void fill_marker_volume(integer iP, integer* &color, FLOW &f) {
	// 0 - не помечено, 1 - помечено, 2 - запрет на пометку.
	integer iE, iN, iT, iW, iS, iB; // номера соседних контрольных объёмов
	iE=f.neighbors_for_the_internal_node[E_SIDE][0][iP]; iN=f.neighbors_for_the_internal_node[N_SIDE][0][iP]; iT=f.neighbors_for_the_internal_node[T_SIDE][0][iP];
	iW=f.neighbors_for_the_internal_node[W_SIDE][0][iP]; iS=f.neighbors_for_the_internal_node[S_SIDE][0][iP]; iB=f.neighbors_for_the_internal_node[B_SIDE][0][iP];

	if ((color[iE]!=2)&&(color[iE]!=1)) {
	    color[iE]=1; // помечено.
	    if (iE<f.maxelm) fill_marker_volume(iE, color, f);
	}
	if ((color[iW]!=2)&&(color[iW]!=1)) {
	    color[iW]=1; // помечено.
	    if (iW<f.maxelm) fill_marker_volume(iW, color, f);
	}
	if ((color[iN]!=2)&&(color[iN]!=1)) {
		color[iN]=1; // помечено.
		if (iN<f.maxelm) fill_marker_volume(iN, color, f);
	}
	if ((color[iS]!=2)&&(color[iS]!=1)) {
	    color[iS]=1; // помечено.
	    if (iS<f.maxelm) fill_marker_volume(iS, color, f);
	}
	if ((color[iT]!=2)&&(color[iT]!=1)) {
	    color[iT]=1; // помечено.
	    if (iT<f.maxelm) fill_marker_volume(iT, color, f);
	}
	if ((color[iB]!=2)&&(color[iB]!=1)) {
	    color[iB]=1; // помечено.
	    if (iB<f.maxelm) fill_marker_volume(iB, color, f);
	}
} // fill_marker_volume

// Раскрашивает связный объём.
// помечает всю связную окрестность связную с заданной нормалью в окрестности точки iP.
// 0 - не помечено, 1 - помечено, 2 - запрет на пометку.
void fill_marker_volume_iterat(integer iP, integer* &color, FLOW &f) {
	// 0 - не помечено, 1 - помечено, 2 - запрет на пометку.
	integer iE, iN, iT, iW, iS, iB; // номера соседних контрольных объёмов

	// Итеративный алгоритм.
	// Пометим стартовую точку.
	// Будем перебирать все точки до тех пор пока все соседи каждой помеченой точки не станут также помечены.

	color[iP]=1;
	bool bweshouldbecontinue=true;
	while (bweshouldbecontinue) {
		bweshouldbecontinue=false;
		for (iP=0; iP<f.maxelm; iP++) {
			if (color[iP]==1) {
			   iE=f.neighbors_for_the_internal_node[E_SIDE][0][iP]; iN=f.neighbors_for_the_internal_node[N_SIDE][0][iP]; iT=f.neighbors_for_the_internal_node[T_SIDE][0][iP];
	           iW=f.neighbors_for_the_internal_node[W_SIDE][0][iP]; iS=f.neighbors_for_the_internal_node[S_SIDE][0][iP]; iB=f.neighbors_for_the_internal_node[B_SIDE][0][iP];

			   if ((color[iE]!=2)&&(color[iE]!=1)) {
	              color[iE]=1; // помечено.
	              bweshouldbecontinue=true;
	           }
	           if ((color[iW]!=2)&&(color[iW]!=1)) {
	               color[iW]=1; // помечено.
	               bweshouldbecontinue=true;
	           }
	           if ((color[iN]!=2)&&(color[iN]!=1)) {
		          color[iN]=1; // помечено.
		          bweshouldbecontinue=true;
	           }
	           if ((color[iS]!=2)&&(color[iS]!=1)) {
	              color[iS]=1; // помечено.
	              bweshouldbecontinue=true;
	           }
	           if ((color[iT]!=2)&&(color[iT]!=1)) {
	              color[iT]=1; // помечено.
	              bweshouldbecontinue=true;
	           }
	           if ((color[iB]!=2)&&(color[iB]!=1)) {
	              color[iB]=1; // помечено.
	              bweshouldbecontinue=true;
	           }

			}
		}
	}	
} // fill_marker_volume_iterat

// убираем лишние разделяющие элементы.
// Очистка лишних приповерхностных разделителей.
// 30 апреля 2013 года.
void correct_marker_surface(integer* &color, integer iNorm, FLOW &f) 
{
	bool brepeat=true;
	while (brepeat) {
		brepeat=false;
		for (integer i=0; i<f.maxelm+f.maxbound; ++i) {
			if ((i<f.maxelm) && (color[i]==1)) {
				integer i1, i2, i3, i4;
				switch (iNorm) {
		          case xnorm:
					  i1=f.neighbors_for_the_internal_node[W_SIDE][0][i];
					  if (color[i1]==1) {
						  if (i1>=f.maxelm) {
							  // граничный узел.
							  i2=i;
							  i3=f.neighbors_for_the_internal_node[E_SIDE][0][i];
							  if (color[i3]==1) {
								  if (i3<f.maxelm) {
									  // стираем выделение.
                                      color[i1]=0;
									  color[i2]=0;
									  color[i3]=0;
									  brepeat=true;
									  i4=f.neighbors_for_the_internal_node[E_SIDE][0][i3];
									  if (color[i4]==1) {
										  color[i4]=0;
									  }
								  }
							  }
						  }
						  else {
							  // узел i и сосед W узла i внутренние и отмеченные.
							  i2=i1;
							  i1=f.neighbors_for_the_internal_node[W_SIDE][0][i2];
							  if (color[i1]==1) {
								  if (i1>=f.maxelm) {
									  // i1 - граничный узел.
								     i3=i;
									 // стираем выделение.
                                     color[i1]=0;
									 color[i2]=0;
									 color[i3]=0;
									 brepeat=true;
									 if (i3<f.maxelm) {
									    i4=f.neighbors_for_the_internal_node[E_SIDE][0][i3];
									    if (color[i4]==1) {
										    color[i4]=0;
									    }
									 }
								  }
							  }
							  else {
								  i3=i2;
								  i2=i;
								  i1=f.neighbors_for_the_internal_node[E_SIDE][0][i];
								  if (color[i1]==1) {
								     if (i1>=f.maxelm) {
									     // i1 - граничный узел.
								     	 // стираем выделение.
                                         color[i1]=0;
									     color[i2]=0;
									     color[i3]=0;
									     brepeat=true;
										 if (i3<f.maxelm) {
									        i4=f.neighbors_for_the_internal_node[W_SIDE][0][i3];
									        if (color[i4]==1) {
										        color[i4]=0;
									        }
									     }
								     }
							      }
							  }
						  }
					  }
					  else {
						  i1=i;
						  i2=f.neighbors_for_the_internal_node[E_SIDE][0][i];
						  if (color[i2]==1) {
							  i3=f.neighbors_for_the_internal_node[E_SIDE][0][i2];
							  if (color[i3]==1) {
								  if (i3>=f.maxelm) {
									  // граничный узел.
									  // стираем выделение.
									  color[i1]=0;
									  color[i2]=0;
									  color[i3]=0;
									  brepeat=true;
									   if (i1<f.maxelm) {
									      i4=f.neighbors_for_the_internal_node[W_SIDE][0][i1];
									      if (color[i4]==1) {
										       color[i4]=0;
									      }
									   }
								  }
							  }
						  }
					  }
			          break;
		          case ynorm:
					  i1=f.neighbors_for_the_internal_node[S_SIDE][0][i];
					  if (color[i1]==1) {
						  if (i1>=f.maxelm) {
							  // граничный узел.
							  i2=i;
							  i3=f.neighbors_for_the_internal_node[N_SIDE][0][i];
							  if (color[i3]==1) {
								  if (i3<f.maxelm) {
									  // стираем выделение.
                                      color[i1]=0;
									  color[i2]=0;
									  color[i3]=0;
									  brepeat=true;
									  i4=f.neighbors_for_the_internal_node[N_SIDE][0][i3];
									  if (color[i4]==1) {
										  color[i4]=0;
									  }
								  }
							  }
						  }
						  else {
							  // узел i и сосед W узла i внутренние и отмеченные.
							  i2=i1;
							  i1=f.neighbors_for_the_internal_node[S_SIDE][0][i2];
							  if (color[i1]==1) {
								  if (i1>=f.maxelm) {
									  // i1 - граничный узел.
								     i3=i;
									 // стираем выделение.
                                     color[i1]=0;
									 color[i2]=0;
									 color[i3]=0;
									 brepeat=true;
									 if (i3<f.maxelm) {
									    i4=f.neighbors_for_the_internal_node[N_SIDE][0][i3];
									    if (color[i4]==1) {
										    color[i4]=0;
									    }
									 }
								  }
							  }
							  else {
								  i3=i2;
								  i2=i;
								  i1=f.neighbors_for_the_internal_node[N_SIDE][0][i];
								  if (color[i1]==1) {
								     if (i1>=f.maxelm) {
									     // i1 - граничный узел.
								     	 // стираем выделение.
                                         color[i1]=0;
									     color[i2]=0;
									     color[i3]=0;
									     brepeat=true;
										 if (i3<f.maxelm) {
									        i4=f.neighbors_for_the_internal_node[S_SIDE][0][i3];
									        if (color[i4]==1) {
										        color[i4]=0;
									        }
									     }
								     }
							      }
							  }
						  }
					  }
					  else {
						  i1=i;
						  i2=f.neighbors_for_the_internal_node[N_SIDE][0][i];
						  if (color[i2]==1) {
							  i3=f.neighbors_for_the_internal_node[N_SIDE][0][i2];
							  if (color[i3]==1) {
								  if (i3>=f.maxelm) {
									  // граничный узел.
									  // стираем выделение.
									  color[i1]=0;
									  color[i2]=0;
									  color[i3]=0;
									  brepeat=true;
									   if (i1<f.maxelm) {
									      i4=f.neighbors_for_the_internal_node[S_SIDE][0][i1];
									      if (color[i4]==1) {
										       color[i4]=0;
									      }
									   }
								  }
							  }
						  }
					  }
			          break;
		          case znorm:
					  i1=f.neighbors_for_the_internal_node[B_SIDE][0][i];
					  if (color[i1]==1) {
						  if (i1>=f.maxelm) {
							  // граничный узел.
							  i2=i;
							  i3=f.neighbors_for_the_internal_node[T_SIDE][0][i];
							  if (color[i3]==1) {
								  if (i3<f.maxelm) {
									  // стираем выделение.
                                      color[i1]=0;
									  color[i2]=0;
									  color[i3]=0;
									  brepeat=true;
									  i4=f.neighbors_for_the_internal_node[T_SIDE][0][i3];
									  if (color[i4]==1) {
										  color[i4]=0;
									  }
								  }
							  }
						  }
						  else {
							  // узел i и сосед W узла i внутренние и отмеченные.
							  i2=i1;
							  i1=f.neighbors_for_the_internal_node[B_SIDE][0][i2];
							  if (color[i1]==1) {
								  if (i1>=f.maxelm) {
									  // i1 - граничный узел.
								     i3=i;
									 // стираем выделение.
                                     color[i1]=0;
									 color[i2]=0;
									 color[i3]=0;
									 brepeat=true;
									 if (i3<f.maxelm) {
									    i4=f.neighbors_for_the_internal_node[T_SIDE][0][i3];
									    if (color[i4]==1) {
										    color[i4]=0;
									    }
									 }
								  }
							  }
							  else {
								  i3=i2;
								  i2=i;
								  i1=f.neighbors_for_the_internal_node[T_SIDE][0][i];
								  if (color[i1]==1) {
								     if (i1>=f.maxelm) {
									     // i1 - граничный узел.
								     	 // стираем выделение.
                                         color[i1]=0;
									     color[i2]=0;
									     color[i3]=0;
									     brepeat=true;
										 if (i3<f.maxelm) {
									        i4=f.neighbors_for_the_internal_node[B_SIDE][0][i3];
									        if (color[i4]==1) {
										        color[i4]=0;
									        }
									     }
								     }
							      }
							  }
						  }
					  }
					  else {
						  i1=i;
						  i2=f.neighbors_for_the_internal_node[T_SIDE][0][i];
						  if (color[i2]==1) {
							  i3=f.neighbors_for_the_internal_node[T_SIDE][0][i2];
							  if (color[i3]==1) {
								  if (i3>=f.maxelm) {
									  // граничный узел.
									  // стираем выделение.
									  color[i1]=0;
									  color[i2]=0;
									  color[i3]=0;
									  brepeat=true;
									   if (i1<f.maxelm) {
									      i4=f.neighbors_for_the_internal_node[B_SIDE][0][i1];
									      if (color[i4]==1) {
										       color[i4]=0;
									      }
									   }
								  }
							  }
						  }
					  }
			          break;
				}

			}
		}
	}


	// Финишная зачистка.
	// Если есть внутренний непомеченный граничный КО, а у него есть 
	// помеченный граничный то снять выделение с граничного !!!.
	for (integer iP=0; iP<f.maxelm; iP++) {
			if (color[iP]==0) {
                integer iE, iN, iT, iW, iS, iB; // номера соседних контрольных объёмов
	            iE=f.neighbors_for_the_internal_node[E_SIDE][0][iP]; iN=f.neighbors_for_the_internal_node[N_SIDE][0][iP]; iT=f.neighbors_for_the_internal_node[T_SIDE][0][iP];
				iW=f.neighbors_for_the_internal_node[W_SIDE][0][iP]; iS=f.neighbors_for_the_internal_node[S_SIDE][0][iP]; iB=f.neighbors_for_the_internal_node[B_SIDE][0][iP];

			    if ((iE>=f.maxelm)&&(color[iE]==1)) {
					color[iE]=0;
				}
				if ((iW>=f.maxelm)&&(color[iW]==1)) {
					color[iW]=0;
				}
				if ((iN>=f.maxelm)&&(color[iN]==1)) {
					color[iN]=0;
				}
				if ((iS>=f.maxelm)&&(color[iS]==1)) {
					color[iS]=0;
				}
				if ((iT>=f.maxelm)&&(color[iT]==1)) {
					color[iT]=0;
				}
				if ((iB>=f.maxelm)&&(color[iB]==1)) {
					color[iB]=0;
				}
			}
		}
	
} // correct_marker_surface

// нужно дополнить разделитель до тройного слоя
// для ILU(2) (также это подойдёт с запасом и для ILU(0)).
// Дополнение будет симметричным: если закрашена лишь одна линия,
// то дополнительно будут закрашены ещё две линии - сверху и снизу.
// 6 апреля 2013 года.
//void correct_marker_surface_fill(integer* &color, integer iNorm, FLOW &f) 
//{
	// функция не реализована 27.09.2020
//} // correct_marker_surface_fill

integer uluchenie_katchestva_razbieniq(FLOW &f, int* &ifrontregulation, integer igold, 
								   const integer istart, const integer ifinish, 
								   doublereal* &front, integer iPseparator_init) {

        integer iPseparator=iPseparator_init;

        if (bunion_8_08_2015) {
				doublereal delta=1.1; // импирическая постоянная.
		        // добавка для улучшения стабильности.
				// более корректное определение iPseparatorl.
				// 8 августа 2015.
				doublereal avgxpos=0.0;
				doublereal avgypos=0.0;
				doublereal avgzpos=0.0;
				doublereal inumbervol=0;
                for (integer i=istart; i<=ifinish; ++i) {
					if (fabs(front[i]-front[igold])<5.0*delta) {
						integer iPl=ifrontregulation[i];
						if (iPl<f.maxelm) {
							 TOCHKA pointP;
			                 center_cord3D(iPl,f.nvtx,f.pa,pointP,100);
							 avgxpos+=pointP.x;
							 avgypos+=pointP.y;
							 avgzpos+=pointP.z;
							 inumbervol++;
						}
					}
                }
				doublereal mindistance2=1.0e35;
				avgxpos/=inumbervol;
				avgypos/=inumbervol;
				avgzpos/=inumbervol;
				for (integer i=istart; i<=ifinish; ++i) {
					if (fabs(front[i]-front[igold])<5.0*delta) {
						integer iPl=ifrontregulation[i];
						if (iPl<f.maxelm) {
							 TOCHKA pointP;
			                 center_cord3D(iPl,f.nvtx,f.pa,pointP,100);
							 
							 if (((pointP.x-avgxpos)*(pointP.x-avgxpos)+(pointP.y-avgypos)*(pointP.y-avgypos)+(pointP.z-avgzpos)*(pointP.z-avgzpos))<mindistance2) {
                                 mindistance2=((pointP.x-avgxpos)*(pointP.x-avgxpos)+(pointP.y-avgypos)*(pointP.y-avgypos)+(pointP.z-avgzpos)*(pointP.z-avgzpos));
								 iPseparator=iPl;
							 }
						}
					}
				}

         }

        // В качестве iPseparator нам нужен строго внутренний контрольный объём.
		if (iPseparator>=f.maxelm) {
			 for (integer i=istart; i<=ifinish; ++i) {
			   integer iP=ifrontregulation[i];
			   if (iP<f.maxelm) {
				   integer iE, iN, iT, iW, iS, iB; // номера соседних контрольных объёмов
	               iE=f.neighbors_for_the_internal_node[E_SIDE][0][iP]; iN=f.neighbors_for_the_internal_node[N_SIDE][0][iP]; iT=f.neighbors_for_the_internal_node[T_SIDE][0][iP]; iW=f.neighbors_for_the_internal_node[W_SIDE][0][iP]; iS=f.neighbors_for_the_internal_node[S_SIDE][0][iP]; iB=f.neighbors_for_the_internal_node[B_SIDE][0][iP];

				   if ((iPseparator==iE)||(iPseparator==iW)||(iPseparator==iN)||(iPseparator==iS)||(iPseparator==iT)||(iPseparator==iB)) {
					   iPseparator=iP; // именно внутренний контрольный объём.
					   break;
				   }
			   }
			 }
		}

		return iPseparator;
}

// три клетки в разделителе.
integer my_separator3(FLOW &f, int* &ifrontregulation, const integer istart, const integer ifinish, integer iNorm,
					  doublereal* &front, doublereal* &front_orig,  PARBOUND &b0, bool bseparate, doublereal &acor,
					  bool &brestart, TEMPER &t, FLOW* &fglobal, integer flow_interior, doublereal mul, 
					  integer &irules_restart, doublereal *rules_step, integer isize_rules_step,
					  integer &ibestvar, bool &bbestvar, integer &minimumdisbalance, doublereal &acorbestvar) {

		// Если brestart==true, то мы поднимаемся вверх и осуществляем перезапуск с новым параметром acor.
		// Если false то всё впорядке разбиение прошло успешно.

		// Если bseparate равен истине то мы финишно делим, а если ложь то просто вычисляем новые границы деления.
	    // acor - корректирующая поправка на деление если она равна нулю то коррекции нету, если она меньше нуля то правая часть меньше левой
		// и левую надо уменьшить, если она больше нуля то левая часть меньше правой и левую часть нужно увеличить.
        //acor=0.05; 
		integer igold=istart+static_cast<integer>((0.5+acor)*(ifinish-istart));
		integer iPseparator=ifrontregulation[igold];
		iPseparator=uluchenie_katchestva_razbieniq(f, ifrontregulation, igold, istart, ifinish, front, iPseparator);

		
        integer iPsepNext, iPsepBack;
		//integer iPsepNext2, iPsepBack2;
		integer iPleft, iPright; // левый и правый кусок расчетной области.
		// В качестве iPleft нам нужен строго внутренний контрольный объём.
		iPleft=ifrontregulation[istart+static_cast<integer>(0.25*(ifinish-istart))];
		iPleft=uluchenie_katchestva_razbieniq(f, ifrontregulation, istart+static_cast<integer>(0.25*(ifinish-istart)), istart, ifinish, front, iPleft);

		
		// В качестве iPright нам нужен строго внутренний контрольный объём.
        iPright=ifrontregulation[istart+static_cast<integer>(0.75*(ifinish-istart))];
		iPright=uluchenie_katchestva_razbieniq(f, ifrontregulation, istart+static_cast<integer>(0.75*(ifinish-istart)), istart, ifinish, front, iPright);

		// Если мы упёрлись в стенку то отодвигаемся от неё.
		switch (iNorm) {
		  case xnorm: iPsepNext=f.neighbors_for_the_internal_node[E_SIDE][0][iPseparator];
			           iPsepBack=f.neighbors_for_the_internal_node[W_SIDE][0][iPseparator];
					   if (iPsepNext>=f.maxelm) {
						   // сдвиг влево в сторону West.
						   iPsepNext=iPseparator;
						   iPseparator=iPsepBack;
						   iPsepBack=f.neighbors_for_the_internal_node[W_SIDE][0][iPseparator];
					   }
					   if (iPsepBack>=f.maxelm) {
						   // сдвиг вправо в сторону East.
						   iPsepBack=iPseparator;
						   iPseparator=iPsepNext;
						   iPsepNext=f.neighbors_for_the_internal_node[E_SIDE][0][iPseparator];
					   }
					   /*
					   iPsepNext2=f.neighbors_for_the_internal_node[ESIDE][0][iPsepNext];
			           iPsepBack2=f.neighbors_for_the_internal_node[WSIDE][ iPsepBack].iNODE1;
					   if (iPsepNext2>=f.maxelm) {
						   // сдвиг влево в сторону West.
						   iPsepNext2=iPsepNext;
						   iPsepNext=iPseparator;
						   iPseparator=iPsepBack;
						   iPsepBack=f.neighbors_for_the_internal_node[WSIDE][0][iPseparator];
						   iPsepBack2=f.neighbors_for_the_internal_node[WSIDE][0][iPsepBack];
					   }
					   if (iPsepBack2>=f.maxelm) {
						   // сдвиг вправо в сторону East.
						   iPsepBack2=iPsepBack;
						   iPsepBack=iPseparator;
						   iPseparator=iPsepNext;
						   iPsepNext=f.neighbors_for_the_internal_node[ESIDE][0][iPseparator];
						   iPsepNext2=f.neighbors_for_the_internal_node[ESIDE][0][iPsepNext];
					   }
					   */
			           break;
		  case ynorm: iPsepNext=f.neighbors_for_the_internal_node[N_SIDE][0][iPseparator];
			           iPsepBack=f.neighbors_for_the_internal_node[S_SIDE][0][iPseparator];
					   if (iPsepNext>=f.maxelm) {
						   // сдвиг влево в сторону S.
						   iPsepNext=iPseparator;
						   iPseparator=iPsepBack;
						   iPsepBack=f.neighbors_for_the_internal_node[S_SIDE][0][iPseparator];
					   }
					   if (iPsepBack>=f.maxelm) {
						   // сдвиг вправо в сторону N.
						   iPsepBack=iPseparator;
						   iPseparator=iPsepNext;
						   iPsepNext = f.neighbors_for_the_internal_node[N_SIDE][0][iPseparator];
					   }
					   /*
					   iPsepNext2=f.neighbors_for_the_internal_node[NSIDE][0][iPsepNext];
			           iPsepBack2=f.neighbors_for_the_internal_node[SSIDE][ iPsepBack].iNODE1;
					   if (iPsepNext2>=f.maxelm) {
						   // сдвиг влево в сторону S.
						   iPsepNext2=iPsepNext;
						   iPsepNext=iPseparator;
						   iPseparator=iPsepBack;
						   iPsepBack=f.neighbors_for_the_internal_node[SSIDE][0][iPseparator];
						   iPsepBack2=f.neighbors_for_the_internal_node[SSIDE][0][iPsepBack];
					   }
					   if (iPsepBack2>=f.maxelm) {
						   // сдвиг вправо в сторону N.
						   iPsepBack2=iPsepBack;
						   iPsepBack=iPseparator;
						   iPseparator=iPsepNext;
						   iPsepNext=f.neighbors_for_the_internal_node[NSIDE][0][iPseparator];
						   iPsepNext2=f.neighbors_for_the_internal_node[NSIDE][0][iPsepNext];
					   }
					   */
			           break;
		  case znorm: iPsepNext = f.neighbors_for_the_internal_node[T_SIDE][0][iPseparator];
			  iPsepBack = f.neighbors_for_the_internal_node[B_SIDE][0][iPseparator];
					   if (iPsepNext>=f.maxelm) {
						   // сдвиг влево в сторону Bottom.
						   iPsepNext=iPseparator;
						   iPseparator=iPsepBack;
						   iPsepBack = f.neighbors_for_the_internal_node[B_SIDE][0][iPseparator];
					   }
					   if (iPsepBack>=f.maxelm) {
						   // сдвиг вправо в сторону Top.
						   iPsepBack=iPseparator;
						   iPseparator=iPsepNext;
						   iPsepNext = f.neighbors_for_the_internal_node[T_SIDE][0][iPseparator];
					   }
					   /*
					   iPsepNext2=f.neighbors_for_the_internal_node[TSIDE][0][iPsepNext];
			           iPsepBack2=f.neighbors_for_the_internal_node[BSIDE][ iPsepBack].iNODE1;
					    if (iPsepNext2>=f.maxelm) {
						   // сдвиг влево в сторону Bottom.
                           iPsepNext2=iPsepNext;
						   iPsepNext=iPseparator;
						   iPseparator=iPsepBack;
						   iPsepBack=f.neighbors_for_the_internal_node[BSIDE][0][iPseparator];
						   iPsepBack2=f.neighbors_for_the_internal_node[BSIDE][0][iPsepBack];
					   }
					   if (iPsepBack2>=f.maxelm) {
						   // сдвиг вправо в сторону Top.
						   iPsepBack2=iPsepBack;
						   iPsepBack=iPseparator;
						   iPseparator=iPsepNext;
						   iPsepNext=f.neighbors_for_the_internal_node[TSIDE][0][iPseparator];
						   iPsepNext2=f.neighbors_for_the_internal_node[TSIDE][0][iPsepNext];
					   }
					   */
			           break;
		}
		// Коррекция igold:
		for (integer i=0; i<f.maxelm+f.maxbound; ++i) {
			integer iP=ifrontregulation[i];
			if (iP==iPseparator) {
				igold=i;
				break;
			}
		}

	    integer *color=nullptr;
		color=new integer[f.maxelm+f.maxbound];
		integer *color_left=nullptr;
		color_left=new integer[f.maxelm+f.maxbound];
		integer *color_right=nullptr;
		color_right=new integer[f.maxelm+f.maxbound];
		if ((color_right==nullptr)||(color_left==nullptr)||(color==nullptr)) 
		{
             // недостаточно памяти на данном оборудовании.
		     printf("Problem: not enough memory on your equipment...\n");
			 printf("Please any key to exit...\n");
			 exit(1);
		}
		bool bOk=false;
		bool bfirst=true; // первый заход.
		while (!bOk) { // 1
		     for (integer i=0; i<f.maxelm+f.maxbound; ++i) {
			      if (bfirst) color[i]=2; // нетрогать  точки помеченные двойкой.
			      color_left[i]=2;
			      color_right[i]=2;
		     }

		     for (integer i=istart; i<=ifinish; ++i) {
			      integer iP=ifrontregulation[i];
			      if (bfirst) color[iP]=0; // можно помечать !
			      color_left[iP]=0;
			      color_right[iP]=0;
		     }


			 integer iPs, iPn, iPb;
			// integer iPnn, iPbb;
			 if (bfirst) {
				 iPs=iPseparator;
				 iPn=iPsepNext;
				 iPb=iPsepBack;
				 //iPnn=iPsepNext2;
				 //iPbb=iPsepBack2;

			 } else { //2
				 // Нужно найти новую точку iPs и её соседей так чтобы граница была разделяющей расчётную область на две половины.
				 // Это случай многосвязной границы разделителя, такое тоже может быть!.
				 // Задача хорошего разделителя поместить всю разделяющую границу внутрь разделителя color.

				 TOCHKA pointP;
			     center_cord3D(iPseparator,f.nvtx,f.pa,pointP,100);


				 doublereal riso_value=0.0;
				 switch (iNorm) {
				    case xnorm: riso_value=pointP.x;
			                     break;
		            case ynorm: riso_value=pointP.y;
			                     break;
					case znorm: riso_value=pointP.z;
			                     break;
		         }

				 integer iPcandidate=iPseparator;

				 bool bfound=false;
				 for (integer i1=1; ((igold+i1)<=ifinish)&&((igold-i1)>=istart); i1++) {
					  const doublereal epsilon=1.0e-35; // для отделения вещественного нуля.
					 integer iP=ifrontregulation[igold+i1];
					 if (iP<f.maxelm) {
					      TOCHKA pointP;
			              center_cord3D(iP,f.nvtx,f.pa,pointP,100);
					 
					      doublereal rval=0.0;
					      switch (iNorm) {
				            case xnorm: rval=pointP.x;
			                             break;
		                    case ynorm: rval=pointP.y;
			                             break;
					        case znorm: rval=pointP.z;
			                             break;
		                  }

				          if (fabs(rval-riso_value)<epsilon) {
					          if ((iP<f.maxelm)&&(color[iP]==0)) {
						          iPcandidate=iP;
						          bfound=true;
						          break;
					          }
					      }
					 }
					 iP=ifrontregulation[igold-i1];
					 if (iP<f.maxelm) {
					     TOCHKA pointP;
			             center_cord3D(iP,f.nvtx,f.pa,pointP,100);
					     doublereal rval=0.0;
					     switch (iNorm) {
				            case xnorm: rval=pointP.x;
			                           break;
		                    case ynorm: rval=pointP.y;
			                           break;
					        case znorm: rval=pointP.z;
			                           break;
		                 }

				         if (fabs(rval-riso_value)<epsilon) {

					           if ((iP<f.maxelm)&&(color[iP]==0)) {
						          iPcandidate=iP;
						          bfound=true;
						          break;
					           }
					      }
					 }
					 
				 }

				 if (bfound) {
				     iPs=iPcandidate;
				     switch (iNorm) {
					 case xnorm: iPn = f.neighbors_for_the_internal_node[E_SIDE][0][iPs];
						 iPb = f.neighbors_for_the_internal_node[W_SIDE][0][iPs];
									 if (iPn>=f.maxelm) {
						                // сдвиг влево в сторону West.
						                iPn=iPs;
						                iPs=iPb;
										iPb = f.neighbors_for_the_internal_node[W_SIDE][0][iPs];
					                 }
					                 if (iPb>=f.maxelm) {
						                // сдвиг вправо в сторону East.
						                iPb=iPs;
						                iPs=iPn;
										iPn = f.neighbors_for_the_internal_node[E_SIDE][0][iPs];
					                 }
									 /*
									 iPnn=f.neighbors_for_the_internal_node[ESIDE][0][iPn];
			                         iPbb=f.neighbors_for_the_internal_node[WSIDE][0][iPb];
									 if (iPnn>=f.maxelm) {
						                // сдвиг влево в сторону West.
										 iPnn=iPn;
						                iPn=iPs;
						                iPs=iPb;
						                iPb=f.neighbors_for_the_internal_node[WSIDE][0][iPs];
										iPbb=f.neighbors_for_the_internal_node[WSIDE][0][iPb];
					                 }
					                 if (iPbb>=f.maxelm) {
						                // сдвиг вправо в сторону East.
										iPbb=iPb;
						                iPb=iPs;
						                iPs=iPn;
						                iPn=f.neighbors_for_the_internal_node[ESIDE][0][iPs];
										iPnn=f.neighbors_for_the_internal_node[ESIDE][0][iPn];
					                 }
									 */
			                         break;
					 case ynorm: iPn = f.neighbors_for_the_internal_node[N_SIDE][0][iPs];
						 iPb = f.neighbors_for_the_internal_node[S_SIDE][0][iPs];
									 if (iPn>=f.maxelm) {
						                // сдвиг влево в сторону West.
						                iPn=iPs;
						                iPs=iPb;
										iPb = f.neighbors_for_the_internal_node[S_SIDE][0][iPs];
					                 }
					                 if (iPb>=f.maxelm) {
						                // сдвиг вправо в сторону East.
						                iPb=iPs;
						                iPs=iPn;
										iPn = f.neighbors_for_the_internal_node[N_SIDE][0][iPs];
					                 }
									 /*
									 iPnn=f.neighbors_for_the_internal_node[NSIDE][0][iPn];
			                         iPbb=f.neighbors_for_the_internal_node[SSIDE][0][iPb];
									 if (iPnn>=f.maxelm) {
						                // сдвиг влево в сторону West.
										iPnn=iPn;
						                iPn=iPs;
						                iPs=iPb;
						                iPb=f.neighbors_for_the_internal_node[SSIDE][0][iPs];
										iPbb=f.neighbors_for_the_internal_node[SSIDE][0][iPb];
					                 }
					                 if (iPbb>=f.maxelm) {
						                // сдвиг вправо в сторону East.
										iPbb=iPb;
						                iPb=iPs;
						                iPs=iPn;
						                iPn=f.neighbors_for_the_internal_node[NSIDE][0][iPs];
										iPnn=f.neighbors_for_the_internal_node[NSIDE][0][iPn];
					                 }
									 */
			                         break;
					 case znorm: iPn = f.neighbors_for_the_internal_node[T_SIDE][0][iPs];
						 iPb = f.neighbors_for_the_internal_node[B_SIDE][0][iPs];
									 if (iPn>=f.maxelm) {
						                // сдвиг влево в сторону West.
						                iPn=iPs;
						                iPs=iPb;
										iPb = f.neighbors_for_the_internal_node[B_SIDE][0][iPs];
					                 }
					                 if (iPb>=f.maxelm) {
						                // сдвиг вправо в сторону East.
						                iPb=iPs;
						                iPs=iPn;
										iPn = f.neighbors_for_the_internal_node[T_SIDE][0][iPs];
					                 }
									 /*
									 iPnn=f.neighbors_for_the_internal_node[TSIDE][0][iPn];
			                         iPbb=f.neighbors_for_the_internal_node[BSIDE][0][iPb];
									 if (iPnn>=f.maxelm) {
						                // сдвиг влево в сторону West.
										iPnn=iPn;
						                iPn=iPs;
						                iPs=iPb;
						                iPb=f.neighbors_for_the_internal_node[BSIDE][0][iPs];
										iPbb=f.neighbors_for_the_internal_node[BSIDE][0][iPb];
					                 }
					                 if (iPbb>=f.maxelm) {
						                // сдвиг вправо в сторону East.
										iPbb=iPb;
						                iPb=iPs;
						                iPs=iPn;
						                iPn=f.neighbors_for_the_internal_node[TSIDE][0][iPs];
										iPnn=f.neighbors_for_the_internal_node[TSIDE][0][iPn];
					                 }
									 */
			                         break;
		             }
				 }
				 else { // 3
					 // Мы не смогли найти кандидатов на продление сепаратора с таким же riso_value.

					 /* Случай когда границу сепаратора невозможно определить 
					 может свидетельствовать и тому, что сепаратор прилегает к
					 стенке расчётной области и фактически не делит расчётную область на две части.
					 Лечение этого может быть следующим, пошевелить границу сепаратора на 1% выполнив рестарт.
					 */

					 printf("can not determine separator surface...\n");
					 printf("WARNING ! may be cfd model is incorrect...\n");
					 printf("Please, press any key to continue calculation.\n");
					 //system("pause");

					 acor+=mul*0.005; // сдвигаем разделитель на 1%.
					 brestart=true;

					 delete[] color;
					 delete[] color_left;
					 delete[] color_right;

					 return 1; // досрочное прекращение работы, выход из функции.


				 } // 3 if

			 } // 2 if

			
		     // Рекурсивный алгоритм заливки данной поверхности определяемой точкой iP и нормалью iNorm.
             // помечает всю связную окрестность связную с заданной нормалью в окрестности точки iP.
             // 0 - не помечено, 1 - помечено, 2 - запрет на пометку.
		     // Помечаем разделитель для ILU(2) разложения (подойдёт и для ILU(0)
		     // маркером единицей.
             fill_marker_surface(iPs, color, iNorm, f); // центральный слой.
			 //TOCHKA pointP; // debug solution
			 //center_cord3D(iPs,f.nvtx,f.pa,pointP);
			 //printf("eto on=%e\n",pointP.x);
#if doubleintprecision == 1
			 //printf("iPs=%lld, maxelm=%lld, maxbound=%lld\n",iPs,f.maxelm,f.maxbound);
#else
			 //printf("iPs=%d, maxelm=%d, maxbound=%d\n",iPs,f.maxelm,f.maxbound);
#endif
			 
		     printf("1");
		     //system("pause");
		     fill_marker_surface(iPn, color, iNorm, f); // правый слой.
			 //fill_marker_surface(iPnn, color, iNorm, f);
		     printf("2");
		     //system("pause");
		     fill_marker_surface(iPb, color, iNorm, f); // левый слой.
			 //fill_marker_surface(iPbb, color,iNorm,f);
		     printf("3");
		     //system("pause");

			 // очистка лишних приповерхностных разделителей и уменьшение числа разделяющих элементов
			 // без нарушения разделяющих свойств сепаратора.
			 correct_marker_surface(color, iNorm, f); // убираем лишние разделяющие элементы.

			 bool breturn=true;
			 for (integer i=0; i<f.maxelm+f.maxbound; ++i) {
			     if (color[i]==1) {
					 breturn=false;
				 }
			 }
			 if (breturn) {
				 // В качестве разделителя нет ни одного внутреннего контрольного объёма !!!
				 // и ни одного граничного. Эта ситуация возникает в случае когда разделитель полностью
				 // прилегает к границе расчётной области и по сути не является разделителем.
				 // Внимание перед использованием данного кода обязательно необходимо вызвать correct_marker_surface; !!!
				 acor+=mul*0.005; // сдвигаем разделитель на 1%.
				 printf("Situation empty delimiter. \n");
			     printf("separator is moving plus length 1%% \n");
			     printf("Please, press any key to continue calculation...\n");
			     //system("pause");

				 // Перезапуск.

				 
				delete[] color;
				color = nullptr;				 
				 
				delete[] color_left;
				color_left = nullptr;				 
				
				delete[] color_right;
				color_right = nullptr;
				

				 return 1;

			 }

			 ///*
			 // Для графической визуализации разделителя:
			 doublereal *viewsep = nullptr;
			 viewsep=new doublereal[f.maxelm + f.maxbound];
			 for (integer i=0; i<f.maxelm+f.maxbound; ++i) {
				 if (color[i]==1) {
					 viewsep[i]=1.0;
				 }
				 else {
					  viewsep[i]=0.0;
				 }
			 }
			 if (debug_export_tecplot) {
				 // экспортируем в tecplot только если это действительно необъходимо.
			    exporttecplotxy360T_3D_part2nd(t.maxelm, t.ncell, fglobal, t, flow_interior, 0, viewsep);
			 }
			 printf("separator visible compleate...\n");
		     //system("pause");
			 if (viewsep != nullptr) {
				 delete[] viewsep;
			 }
			 // для графической визуализации разделителя . finish.
			//*/
			 


			 //integer ik7=0;
		     for (integer i=0; i<f.maxelm+f.maxbound; ++i) {
			     if (color[i]==1) {
				    color_left[i]=2;
			        color_right[i]=2;
					//ik7++;
#if doubleintprecision == 1
					//printf("i==%lld, number=%lld\n",i,ik7);
#else
					//printf("i==%d, number=%d\n",i,ik7);
#endif
					
					//system("pause");
			     }
		     }

		     // Выделяем левую подобласть:
		     // Раскрашивает связный объём.
             // помечает всю связную окрестность связную с заданной нормалью в окрестности точки iP.
             // 0 - не помечено, 1 - помечено, 2 - запрет на пометку.
		     // Мы полностью отказались от рекурсии и создали полностью итеративный алгоритм. 24 апреля 2012.
             fill_marker_volume_iterat(iPleft, color_left, f);
			 //fill_marker_volume(iPleft, color_left, f);
		
			 integer icount_non_markerVolume=0;
		     for (integer i=0; i<f.maxelm+f.maxbound; ++i) {
			     if (color_left[i]==0) {
				    bOk=true;
					//printf("incomming\n");
					//system("pause");
					icount_non_markerVolume++;
#if doubleintprecision == 1
					/*if (i>f.maxelm) {
						printf("bound %lld maxelm+maxbound-1=%lld\n",i,f.maxelm+f.maxbound-1);
					}
					else {
						printf("internal %lld maxelm+maxbound-1=%lld\n",i,f.maxelm+f.maxbound-1);
						system("pause");
					}*/
#else
					/*if (i>f.maxelm) {
						printf("bound %d maxelm+maxbound-1=%d\n",i,f.maxelm+f.maxbound-1);
					}
					else {
						printf("internal %d maxelm+maxbound-1=%d\n",i,f.maxelm+f.maxbound-1);
						system("pause");
					}*/
#endif
					
				    //break; // Обязательно раскоментировать !!!!
			     }
		     }
			 bfirst=false;
#if doubleintprecision == 1
			 printf("icount_non_markerVolume=%lld\n", icount_non_markerVolume);
#else
			 printf("icount_non_markerVolume=%d\n", icount_non_markerVolume);
#endif
			 

			 /* bOk - остаётся false если граница сепаратора определена неполностью, а именно
			 граница разделителя имеет дырку через которую осуществлена закраска всей расчётной области.
			 В этом случае(осуществлена закраска всей расчётной области) требуется продолжить нахождение 
			 границы разделителя с тем же iso_value в окрестности данного расстояния от входной границы.
			 */

		} // 1 while 


		printf("4");
		//system("pause");
		// Выделяем правую подобласть:
        // Раскрашивает связный объём.
        // помечает всю связную окрестность связную с заданной нормалью в окрестности точки iP.
        // 0 - не помечено, 1 - помечено, 2 - запрет на пометку.
		// Мы полностью отказались от рекурсии и создали полностью итеративный алгоритм. 24 апреля 2012.
        fill_marker_volume_iterat(iPright, color_right, f);
		printf("5\n");
		//system("pause");

		/* Здесь ещё остаётся возможность неверного разделения:
		Случай неверного разделения определяется следующим образом:
		Все узлы минус узлы помеченные разделителем, минус узлы помеченные левым сепаратором минус узлы помеченные правым сепаратором
		не равно пусто значит разбиение произведено ошбочно. Это может быть ели сепаратор отрезал изолированный кусок расчётной области,
		алевый и правый маркеры находятся на одном куске разбиения отличного от отрезанного обособленного куска.
		*/
		integer *ifillcorect = nullptr;
		ifillcorect = new integer[f.maxelm + f.maxbound];
		for (integer i=0; i<f.maxelm+f.maxbound; ++i) {
			  ifillcorect[i]=2;
		}
		for (integer i=istart; i<=ifinish; ++i) {
		    integer iP=ifrontregulation[i];
		    ifillcorect[iP]=0;
		}
		for (integer i=0; i<f.maxelm+f.maxbound; ++i) {
			if (color[i]==1) {
			   ifillcorect[i]=1;
			}
			if (color_left[i]==1) {
			   ifillcorect[i]=1;
			}
			if (color_right[i]==1) {
			   ifillcorect[i]=1;
			}
		}
		// Наконец поверка.
		bool bfound0=false;
		for (integer i=0; i<f.maxelm+f.maxbound; ++i) {
			if (ifillcorect[i]==0) {
				// разделение некорректно !!!
				bfound0=true;
				break;
			}
		}
		if (!bfound0) {
			// разделение полностью корректно и рестарт ненужен.
			brestart=false;
		}
		if (ifillcorect != nullptr) {
			delete[] ifillcorect; // Освобождение памяти.
		}

		if (!brestart)  { //1
		    // Финишная стадия: переупорядочивание.
		    int *iforig=nullptr;
		    iforig=new int[f.maxelm+f.maxbound];
		    doublereal *rfront_orig=nullptr;
		    rfront_orig=new doublereal[f.maxelm+f.maxbound];
		    for (integer i=0; i<f.maxelm+f.maxbound; ++i) {
			    iforig[i]=ifrontregulation[i]; // создаём копию ifrontregularization.
			    rfront_orig[i]=front[i];
		    }

		    // Препроцессинг для упорядочивания, чтобы выяснить  нужно он или нет.
		    // левая часть.
		    integer icl=0;
		    for (integer i=istart; i<=ifinish; ++i) {
			    integer iP=iforig[i];
			    if (color_left[iP]==1) {
			    	//ifrontregulation[istart+icl]=iP;
			    	// также нужно переупорядочить фронт.
				    //front[istart+icl]=rfront_orig[i];
			    	icl++;
			    }
		    }
		    
		    // правая часть.
            integer icr=0;
		    for (integer i=istart; i<=ifinish; ++i) {
			    integer iP=iforig[i];
			    if (color_right[iP]==1) {
			    	//ifrontregulation[b0.ileft_finish+1+icr]=iP;
				    // также нужно переупорядочить фронт.
			    	//front[b0.ileft_finish+1+icr]=rfront_orig[i];
			    	icr++;
			    }
		    }
			

		    // разделитель.
		    integer ics=0;
		    for (integer i=istart; i<=ifinish; ++i) {
		    	integer iP=iforig[i];
			    if (color[iP]==1) {
			    	//ifrontregulation[b0.iright_finish+1+ics]=iP;
			    	// также нужно переупорядочить фронт.
			     	//front[b0.iright_finish+1+ics]=rfront_orig[i];
			    	ics++;
			    }
		    }

			// 
			// 8 августа 2015 года.
			// этот код предназначен для улучшения качества разбиения.
			// проверка баланса: если нарушен баланс (левая и правая части отличаются более чем
			// на 20% ) то будут осуществлены подвижки границы сепаратора немного пошевелены 
			// до тех пора пока левая  правая части будут равны с точностью выше чем 20%.

			if ((100.0*(icl)/(ifinish-istart+1)<=40.0)||(100.0*(icr)/(ifinish-istart+1)<=40.0)) {
			    irules_restart++;
				if (irules_restart<isize_rules_step) {
					if (irules_restart==isize_rules_step-1) {
						// делаем лучший вариант и заканчиваем.
#if doubleintprecision == 1
						printf("rules move best iteration %lld getready...\n", ibestvar);
#else
						printf("rules move best iteration %d getready...\n", ibestvar);
#endif
						
						acor=acorbestvar;
						brestart=true;
					}
					else {
		                  acor+=3.0*rules_step[irules_restart];
#if doubleintprecision == 1
						  printf("rules move iteration %lld getready...\n", irules_restart);
#else
						  printf("rules move iteration %d getready...\n", irules_restart);
#endif
				          
				          brestart=true;
				          // накопление информации о лучшем варианте.
				          if ((icl-icr)*(icl-icr)<minimumdisbalance) {
					         minimumdisbalance=(icl-icr)*(icl-icr);
					         ibestvar=irules_restart;
							 acorbestvar=acor;
				          }
					}
				}
				else {
					irules_restart=0;
					brestart=false;
				}
			}
			else {
                irules_restart=0; // всё успешно и улучшений не требуется.
			}

			if (!brestart) { //2

			b0.ileft_start=istart;
		    b0.ileft_finish=istart+icl-1;
		    b0.iright_start=b0.ileft_finish+1;
		    b0.iright_finish=b0.ileft_finish+1+icr-1;
		    b0.iseparate_start=b0.iright_finish+1;
		    b0.iseparate_finish=b0.iright_finish+1+ics-1;
#if doubleintprecision == 1
			printf("icl=%lld, icr=%lld, ics=%lld\n", icl, icr, ics);
#else
			printf("icl=%d, icr=%d, ics=%d\n", icl, icr, ics);
#endif
		   
			

		    //system("pause");
		    doublereal rsep=100.0*(b0.iseparate_finish-b0.iseparate_start+1)/(ifinish-istart+1);
		    doublereal rleft=100.0*(b0.ileft_finish-b0.ileft_start+1)/(ifinish-istart+1);
		    doublereal rright=100.0*(b0.iright_finish-b0.iright_start+1)/(ifinish-istart+1);

		

		    if (bseparate) {
		       if ((rleft>rsep)&&(rright>rsep)) {
			       b0.active=true;
		           // Наконец переупорядочивание.
		           // левая часть.
		           icl=0;
		           for (integer i=istart; i<=ifinish; ++i) {
			           int iP=iforig[i];
			           if (color_left[iP]==1) {
				            ifrontregulation[istart+icl]=iP;
				            // также нужно переупорядочить фронт.
				            //front[istart+icl]=rfront_orig[i];
				            front[istart+icl]=front_orig[iP];
				            icl++;
			           }
		           }
		           b0.ileft_start=istart;
		           b0.ileft_finish=istart+icl-1;
		           // правая часть.
                   icr=0;
		           for (integer i=istart; i<=ifinish; ++i) {
			           int iP=iforig[i];
			           if (color_right[iP]==1) {
				           ifrontregulation[b0.ileft_finish+1+icr]=iP;
				           // также нужно переупорядочить фронт.
				           //front[b0.ileft_finish+1+icr]=rfront_orig[i];
				           front[b0.ileft_finish+1+icr]=front_orig[iP];
				           icr++;
			           }
		           }
		           b0.iright_start=b0.ileft_finish+1;
		           b0.iright_finish=b0.ileft_finish+1+icr-1;

		           // разделитель.
		           ics=0;
		           for (integer i=istart; i<=ifinish; ++i) {
			           int iP=iforig[i];
			           if (color[iP]==1) {
				          ifrontregulation[b0.iright_finish+1+ics]=iP;
				          // также нужно переупорядочить фронт.
				          //front[b0.iright_finish+1+ics]=rfront_orig[i];
				          front[b0.iright_finish+1+ics]=front_orig[iP];
				          ics++;
			           }
		           }
		           b0.iseparate_start=b0.iright_finish+1;
		           b0.iseparate_finish=b0.iright_finish+1+ics-1;
		      }
		   }

			} // 2

			if (rfront_orig != nullptr) {
				delete[] rfront_orig;
			}
			if (iforig != nullptr) {
				delete[] iforig;
			}
		}  // 1
		else {
			// Подготовка к рестарту !!!
			acor+=mul*0.005; // сдвигаем разделитель на 1%.
			printf("separator is moving plus length 1%% \n");
			printf("Please, press any key to continue calculation...\n");
			//system("pause");
		}

		// метка для рестарта с другим значением acor.

		delete[] color;
		delete[] color_left;
		delete[] color_right;

		return 0; // успешное деление. выход из функции.
		
		/*

		     // Внешний цикл нужен для улучшения характеристик разбиения.
		     // Всего три улучшающих подхода.
		     for (integer istep=0; istep<3; istep++) {
		        // Вычисляем корректирующую поправку разбиения.
		        if (rleft<rright) {
			        acor=fmin(50.0-rleft,rright-50.0);
			        acor=fabs(acor)/100.0;
		        }
		        else {
			        acor=fmin(50.0-rright,rleft-50.0);
			        acor=-fabs(acor)/100.0;
		        }
             }
		*/
} // my_separator3


// пять клеток в разделителе.
integer my_separator5(FLOW &f, int* &ifrontregulation, const integer istart, const integer ifinish, integer iNorm,
					  doublereal* &front, doublereal* &front_orig,  PARBOUND &b0, bool bseparate, doublereal &acor,
					  bool &brestart, TEMPER &t, FLOW* &fglobal, integer flow_interior, doublereal mul, 
					  integer &irules_restart, doublereal *rules_step, integer isize_rules_step,
					  integer &ibestvar, bool &bbestvar, integer &minimumdisbalance, doublereal &acorbestvar) {

		// Если brestart==true, то мы поднимаемся вверх и осуществляем перезапуск с новым параметром acor.
		// Если false то всё впорядке разбиение прошло успешно.

		// Если bseparate равен истине то мы финишно делим, а если ложь то просто вычисляем новые границы деления.
	    // acor - корректирующая поправка на деление если она равна нулю то коррекции нету, если она меньше нуля то правая часть меньше левой
		// и левую надо уменьшить, если она больше нуля то левая часть меньше правой и левую часть нужно увеличить.
        //acor=0.05; 
		integer igold=istart+static_cast<integer>((0.5+acor)*(ifinish-istart));
		integer iPseparator=ifrontregulation[igold];
		iPseparator=uluchenie_katchestva_razbieniq(f, ifrontregulation, igold, istart, ifinish, front, iPseparator);

		
        integer iPsepNext, iPsepBack;
		integer iPsepNext2, iPsepBack2;
		integer iPleft, iPright; // левый и правый кусок расчетной области.
		// В качестве iPleft нам нужен строго внутренний контрольный объём.
		iPleft=ifrontregulation[istart+static_cast<integer>(0.25*(ifinish-istart))];
		iPleft=uluchenie_katchestva_razbieniq(f, ifrontregulation, istart+static_cast<integer>(0.25*(ifinish-istart)), istart, ifinish, front, iPleft);

		
		// В качестве iPright нам нужен строго внутренний контрольный объём.
        iPright=ifrontregulation[istart+static_cast<integer>(0.75*(ifinish-istart))];
		iPright=uluchenie_katchestva_razbieniq(f, ifrontregulation, istart+static_cast<integer>(0.75*(ifinish-istart)), istart, ifinish, front, iPright);

		// Если мы упёрлись в стенку то отодвигаемся от неё.
		switch (iNorm) {
		case xnorm: iPsepNext = f.neighbors_for_the_internal_node[E_SIDE][0][iPseparator];
			iPsepBack = f.neighbors_for_the_internal_node[W_SIDE][0][iPseparator];
					   if (iPsepNext>=f.maxelm) {
						   // сдвиг влево в сторону West.
						   iPsepNext=iPseparator;
						   iPseparator=iPsepBack;
						   iPsepBack = f.neighbors_for_the_internal_node[W_SIDE][0][iPseparator];
					   }
					   if (iPsepBack>=f.maxelm) {
						   // сдвиг вправо в сторону East.
						   iPsepBack=iPseparator;
						   iPseparator=iPsepNext;
						   iPsepNext = f.neighbors_for_the_internal_node[E_SIDE][0][iPseparator];
					   }
					   iPsepNext2 = f.neighbors_for_the_internal_node[E_SIDE][0][iPsepNext];
					   iPsepBack2 = f.neighbors_for_the_internal_node[W_SIDE][0][iPsepBack];
					   if (iPsepNext2>=f.maxelm) {
						   // сдвиг влево в сторону West.
						   iPsepNext2=iPsepNext;
						   iPsepNext=iPseparator;
						   iPseparator=iPsepBack;
						   iPsepBack = f.neighbors_for_the_internal_node[W_SIDE][0][iPseparator];
						   iPsepBack2 = f.neighbors_for_the_internal_node[W_SIDE][0][iPsepBack];
					   }
					   if (iPsepBack2>=f.maxelm) {
						   // сдвиг вправо в сторону East.
						   iPsepBack2=iPsepBack;
						   iPsepBack=iPseparator;
						   iPseparator=iPsepNext;
						   iPsepNext = f.neighbors_for_the_internal_node[E_SIDE][0][iPseparator];
						   iPsepNext2 = f.neighbors_for_the_internal_node[E_SIDE][0][iPsepNext];
					   }
			           break;
		case ynorm: iPsepNext = f.neighbors_for_the_internal_node[N_SIDE][0][iPseparator];
			iPsepBack = f.neighbors_for_the_internal_node[S_SIDE][0][iPseparator];
					   if (iPsepNext>=f.maxelm) {
						   // сдвиг влево в сторону S.
						   iPsepNext=iPseparator;
						   iPseparator=iPsepBack;
						   iPsepBack = f.neighbors_for_the_internal_node[S_SIDE][0][iPseparator];
					   }
					   if (iPsepBack>=f.maxelm) {
						   // сдвиг вправо в сторону N.
						   iPsepBack=iPseparator;
						   iPseparator=iPsepNext;
						   iPsepNext = f.neighbors_for_the_internal_node[N_SIDE][0][iPseparator];
					   }
					   iPsepNext2 = f.neighbors_for_the_internal_node[N_SIDE][0][iPsepNext];
					   iPsepBack2 = f.neighbors_for_the_internal_node[S_SIDE][0][iPsepBack];
					   if (iPsepNext2>=f.maxelm) {
						   // сдвиг влево в сторону S.
						   iPsepNext2=iPsepNext;
						   iPsepNext=iPseparator;
						   iPseparator=iPsepBack;
						   iPsepBack = f.neighbors_for_the_internal_node[S_SIDE][0][iPseparator];
						   iPsepBack2 = f.neighbors_for_the_internal_node[S_SIDE][0][iPsepBack];
					   }
					   if (iPsepBack2>=f.maxelm) {
						   // сдвиг вправо в сторону N.
						   iPsepBack2=iPsepBack;
						   iPsepBack=iPseparator;
						   iPseparator=iPsepNext;
						   iPsepNext = f.neighbors_for_the_internal_node[N_SIDE][0][iPseparator];
						   iPsepNext2 = f.neighbors_for_the_internal_node[N_SIDE][0][iPsepNext];
					   }
			           break;
		case znorm: iPsepNext = f.neighbors_for_the_internal_node[T_SIDE][0][iPseparator];
			iPsepBack = f.neighbors_for_the_internal_node[B_SIDE][0][iPseparator];
					   if (iPsepNext>=f.maxelm) {
						   // сдвиг влево в сторону Bottom.
						   iPsepNext=iPseparator;
						   iPseparator=iPsepBack;
						   iPsepBack = f.neighbors_for_the_internal_node[B_SIDE][0][iPseparator];
					   }
					   if (iPsepBack>=f.maxelm) {
						   // сдвиг вправо в сторону Top.
						   iPsepBack=iPseparator;
						   iPseparator=iPsepNext;
						   iPsepNext = f.neighbors_for_the_internal_node[T_SIDE][0][iPseparator];
					   }
					   iPsepNext2 = f.neighbors_for_the_internal_node[T_SIDE][0][iPsepNext];
					   iPsepBack2 = f.neighbors_for_the_internal_node[B_SIDE][0][iPsepBack];
					    if (iPsepNext2>=f.maxelm) {
						   // сдвиг влево в сторону Bottom.
                           iPsepNext2=iPsepNext;
						   iPsepNext=iPseparator;
						   iPseparator=iPsepBack;
						   iPsepBack = f.neighbors_for_the_internal_node[B_SIDE][0][iPseparator];
						   iPsepBack2 = f.neighbors_for_the_internal_node[B_SIDE][0][iPsepBack];
					   }
					   if (iPsepBack2>=f.maxelm) {
						   // сдвиг вправо в сторону Top.
						   iPsepBack2=iPsepBack;
						   iPsepBack=iPseparator;
						   iPseparator=iPsepNext;
						   iPsepNext = f.neighbors_for_the_internal_node[T_SIDE][0][iPseparator];
						   iPsepNext2 = f.neighbors_for_the_internal_node[T_SIDE][0][iPsepNext];
					   }

			           break;
		}
		// Коррекция igold:
		for (integer i=0; i<f.maxelm+f.maxbound; ++i) {
			integer iP=ifrontregulation[i];
			if (iP==iPseparator) {
				igold=i;
				break;
			}
		}

	    integer *color=nullptr;
		color=new integer[static_cast<integer>(f.maxelm)+ static_cast<integer>(f.maxbound)];
		integer *color_left=nullptr;
		color_left=new integer[f.maxelm+f.maxbound];
		integer *color_right=nullptr;
		color_right=new integer[f.maxelm+f.maxbound];
		if ((color_right==nullptr)||(color_left==nullptr)||(color==nullptr)) 
		{
             // недостаточно памяти на данном оборудовании.
		     printf("Problem: not enough memory on your equipment...\n");
			 printf("Please any key to exit...\n");
			 exit(1);
		}
		bool bOk=false;
		bool bfirst=true; // первый заход.
		while (!bOk) {
		     for (integer i=0; i<f.maxelm+f.maxbound; ++i) {
			      if (bfirst) color[i]=2; // нетрогать  точки помеченные двойкой.
			      color_left[i]=2;
			      color_right[i]=2;
		     }

		     for (integer i=istart; i<=ifinish; ++i) {
			      integer iP=ifrontregulation[i];
			      if (bfirst) color[iP]=0; // можно помечать !
			      color_left[iP]=0;
			      color_right[iP]=0;
		     }


			 integer iPs, iPn, iPb;
			 integer iPnn, iPbb;
			 if (bfirst) {
				 iPs=iPseparator;
				 iPn=iPsepNext;
				 iPb=iPsepBack;
				 iPnn=iPsepNext2;
				 iPbb=iPsepBack2;

			 } else { 
				 // Нужно найти новую точку iPs и её соседей так чтобы граница была разделяющей расчётную область на две половины.
				 // Это случай многосвязной границы разделителя, такое тоже может быть!.
				 // Задача хорошего разделителя поместить всю разделяющую границу внутрь разделителя color.

				 TOCHKA pointP;
			     center_cord3D(iPseparator,f.nvtx,f.pa,pointP,100);


				 doublereal riso_value=0.0;
				 switch (iNorm) {
				    case xnorm: riso_value=pointP.x;
			                     break;
		            case ynorm: riso_value=pointP.y;
			                     break;
					case znorm: riso_value=pointP.z;
			                     break;
		         }

				 integer iPcandidate=iPseparator;

				 bool bfound=false;
				 for (integer i1=1; ((igold+i1)<=ifinish)&&((igold-i1)>=istart); i1++) {
					  const doublereal epsilon=1.0e-35; // для отделения вещественного нуля.
					 integer iP=ifrontregulation[igold+i1];
					 if (iP<f.maxelm) {
					      TOCHKA pointP;
			              center_cord3D(iP,f.nvtx,f.pa,pointP,100);
					 
					      doublereal rval=0.0;
					      switch (iNorm) {
				            case xnorm: rval=pointP.x;
			                             break;
		                    case ynorm: rval=pointP.y;
			                             break;
					        case znorm: rval=pointP.z;
			                             break;
		                  }

				          if (fabs(rval-riso_value)<epsilon) {
					          if ((iP<f.maxelm)&&(color[iP]==0)) {
						          iPcandidate=iP;
						          bfound=true;
						          break;
					          }
					      }
					 }
					 iP=ifrontregulation[igold-i1];
					 if (iP<f.maxelm) {
					     TOCHKA pointP;
			             center_cord3D(iP,f.nvtx,f.pa,pointP,100);
					     doublereal rval=0.0;
					     switch (iNorm) {
				            case xnorm: rval=pointP.x;
			                           break;
		                    case ynorm: rval=pointP.y;
			                           break;
					        case znorm: rval=pointP.z;
			                           break;
		                 }

				         if (fabs(rval-riso_value)<epsilon) {

					           if ((iP<f.maxelm)&&(color[iP]==0)) {
						          iPcandidate=iP;
						          bfound=true;
						          break;
					           }
					      }
					 }
					 
				 }

				 if (bfound) {
				     iPs=iPcandidate;
				     switch (iNorm) {
					 case xnorm: iPn = f.neighbors_for_the_internal_node[E_SIDE][0][iPs];
						 iPb = f.neighbors_for_the_internal_node[W_SIDE][0][iPs];
									 if (iPn>=f.maxelm) {
						                // сдвиг влево в сторону West.
						                iPn=iPs;
						                iPs=iPb;
										iPb = f.neighbors_for_the_internal_node[W_SIDE][0][iPs];
					                 }
					                 if (iPb>=f.maxelm) {
						                // сдвиг вправо в сторону East.
						                iPb=iPs;
						                iPs=iPn;
										iPn = f.neighbors_for_the_internal_node[E_SIDE][0][iPs];
					                 }
									 iPnn = f.neighbors_for_the_internal_node[E_SIDE][0][iPn];
									 iPbb = f.neighbors_for_the_internal_node[W_SIDE][0][iPb];
									 if (iPnn>=f.maxelm) {
						                // сдвиг влево в сторону West.
										 iPnn=iPn;
						                iPn=iPs;
						                iPs=iPb;
										iPb = f.neighbors_for_the_internal_node[W_SIDE][0][iPs];
										iPbb = f.neighbors_for_the_internal_node[W_SIDE][0][iPb];
					                 }
					                 if (iPbb>=f.maxelm) {
						                // сдвиг вправо в сторону East.
										iPbb=iPb;
						                iPb=iPs;
						                iPs=iPn;
										iPn = f.neighbors_for_the_internal_node[E_SIDE][0][iPs];
										iPnn = f.neighbors_for_the_internal_node[E_SIDE][0][iPn];
					                 }
			                         break;
					 case ynorm: iPn = f.neighbors_for_the_internal_node[N_SIDE][0][iPs];
						 iPb = f.neighbors_for_the_internal_node[S_SIDE][0][iPs];
									 if (iPn>=f.maxelm) {
						                // сдвиг влево в сторону West.
						                iPn=iPs;
						                iPs=iPb;
										iPb = f.neighbors_for_the_internal_node[S_SIDE][0][iPs];
					                 }
					                 if (iPb>=f.maxelm) {
						                // сдвиг вправо в сторону East.
						                iPb=iPs;
						                iPs=iPn;
										iPn = f.neighbors_for_the_internal_node[N_SIDE][0][iPs];
					                 }
									 iPnn = f.neighbors_for_the_internal_node[N_SIDE][0][iPn];
									 iPbb = f.neighbors_for_the_internal_node[S_SIDE][0][iPb];
									 if (iPnn>=f.maxelm) {
						                // сдвиг влево в сторону West.
										iPnn=iPn;
						                iPn=iPs;
						                iPs=iPb;
										iPb = f.neighbors_for_the_internal_node[S_SIDE][0][iPs];
										iPbb = f.neighbors_for_the_internal_node[S_SIDE][0][iPb];
					                 }
					                 if (iPbb>=f.maxelm) {
						                // сдвиг вправо в сторону East.
										iPbb=iPb;
						                iPb=iPs;
						                iPs=iPn;
										iPn = f.neighbors_for_the_internal_node[N_SIDE][0][iPs];
										iPnn = f.neighbors_for_the_internal_node[N_SIDE][0][iPn];
					                 }
			                         break;
					 case znorm: iPn = f.neighbors_for_the_internal_node[T_SIDE][0][iPs];
						 iPb = f.neighbors_for_the_internal_node[B_SIDE][0][iPs];
									 if (iPn>=f.maxelm) {
						                // сдвиг влево в сторону West.
						                iPn=iPs;
						                iPs=iPb;
										iPb = f.neighbors_for_the_internal_node[B_SIDE][0][iPs];
					                 }
					                 if (iPb>=f.maxelm) {
						                // сдвиг вправо в сторону East.
						                iPb=iPs;
						                iPs=iPn;
										iPn = f.neighbors_for_the_internal_node[T_SIDE][0][iPs];
					                 }
									 iPnn = f.neighbors_for_the_internal_node[T_SIDE][0][iPn];
									 iPbb = f.neighbors_for_the_internal_node[B_SIDE][0][iPb];
									 if (iPnn>=f.maxelm) {
						                // сдвиг влево в сторону West.
										iPnn=iPn;
						                iPn=iPs;
						                iPs=iPb;
										iPb = f.neighbors_for_the_internal_node[B_SIDE][0][iPs];
										iPbb = f.neighbors_for_the_internal_node[B_SIDE][0][iPb];
					                 }
					                 if (iPbb>=f.maxelm) {
						                // сдвиг вправо в сторону East.
										iPbb=iPb;
						                iPb=iPs;
						                iPs=iPn;
										iPn = f.neighbors_for_the_internal_node[T_SIDE][0][iPs];
										iPnn = f.neighbors_for_the_internal_node[T_SIDE][0][iPn];
					                 }
			                         break;
		             }
				 }
				 else {
					 // Мы не смогли найти кандидатов на продление сепаратора с таким же riso_value.

					 /* Случай когда границу сепаратора невозможно определить 
					 может свидетельствовать и тому, что сепаратор прилегает к
					 стенке расчётной области и фактически не делит расчётную область на две части.
					 Лечение этого может быть следующим, пошевелить границу сепаратора на 1% выполнив рестарт.
					 */

					 printf("can not determine separator surface...\n");
					 printf("WARNING ! may be cfd model is incorrect...\n");
					 printf("Please, press any key to continue calculation.\n");
					// system("pause");

					 acor+=mul*0.005; // сдвигаем разделитель на 1%.
					 brestart=true;

					 if (color != nullptr) {
						 delete[] color;
					 }
					 if (color_left != nullptr) {
						 delete[] color_left;
					 }
					 if (color_right != nullptr) {
						 delete[] color_right;
					 }

					 return 1;


				 }

			 }

			
		     // Рекурсивный алгоритм заливки данной поверхности определяемой точкой iP и нормалью iNorm.
             // помечает всю связную окрестность связную с заданной нормалью в окрестности точки iP.
             // 0 - не помечено, 1 - помечено, 2 - запрет на пометку.
		     // Помечаем разделитель для ILU(2) разложения (подойдёт и для ILU(0)
		     // маркером единицей.
             fill_marker_surface(iPs, color, iNorm, f); // центральный слой.
			 //TOCHKA pointP; // debug solution
			 //center_cord3D(iPs,f.nvtx,f.pa,pointP);
			 //printf("eto on=%e\n",pointP.x);
#if doubleintprecision == 1
			 //printf("iPs=%lld, maxelm=%lld, maxbound=%lld\n",iPs,f.maxelm,f.maxbound);
#else
			 //printf("iPs=%d, maxelm=%d, maxbound=%d\n",iPs,f.maxelm,f.maxbound);
#endif
			
		     printf("1");
		     //system("pause");
		     fill_marker_surface(iPn, color, iNorm, f); // правый слой.
			 fill_marker_surface(iPnn, color, iNorm, f);
		     printf("2");
		     //system("pause");
		     fill_marker_surface(iPb, color, iNorm, f); // левый слой.
			 fill_marker_surface(iPbb, color,iNorm,f);
		     printf("3");
		     //system("pause");

			 // очистка лишних приповерхностных разделителей и уменьшение числа разделяющих элементов
			 // без нарушения разделяющих свойств сепаратора.
			 correct_marker_surface(color, iNorm, f); // убираем лишние разделяющие элементы.

			 bool breturn=true;
			 if (color != nullptr) {
				 for (integer i = 0; i < f.maxelm + f.maxbound; ++i) {
					 if (color[i] == 1) {
						 breturn = false;
					 }
				 }
			 }
			 if (breturn) {
				 // В качестве разделителя нет ни одного внутреннего контрольного объёма !!!
				 // и ни одного граничного. Эта ситуация возникает в случае когда разделитель полностью
				 // прилегает к границе расчётной области и по сути не является разделителем.
				 // Внимание перед использованием данного кода обязательно необходимо вызвать correct_marker_surface; !!!
				 acor+=mul*0.005; // сдвигаем разделитель на 1%.
				 printf("Situation empty delimiter. \n");
			     printf("separator is moving plus length 1%% \n");
			     printf("Please, press any key to continue calculation...\n");
			     //system("pause");
				
				 if (color != nullptr) {
					 delete[] color;
				 }
				 if (color_left != nullptr) {
					 delete[] color_left;
				 }
				 if (color_right != nullptr) {
					 delete[] color_right;
				 }

				 return 0;

			 }

			 ///*
			 // Для графической визуализации разделителя:
			 doublereal *viewsep = nullptr;
			 viewsep=new doublereal[f.maxelm + f.maxbound];
			 for (integer i=0; i<f.maxelm+f.maxbound; ++i) {
				 if (color[i]==1) {
					 viewsep[i]=1.0;
				 }
				 else {
					  viewsep[i]=0.0;
				 }
			 }
			  if (debug_export_tecplot) {
				 // экспортируем в tecplot только если это действительно необъходимо.
			     exporttecplotxy360T_3D_part2nd(t.maxelm, t.ncell, fglobal, t, flow_interior, 0, viewsep);
			  }
			 printf("separator visible compleate...\n");
		     //system("pause");
			 if (viewsep != nullptr) {
				 delete[] viewsep;
			 }
			 // для графической визуализации разделителя . finish.
			//*/
			 


			 //integer ik7=0;
		     for (integer i=0; i<f.maxelm+f.maxbound; ++i) {
			     if (color[i]==1) {
				    color_left[i]=2;
			        color_right[i]=2;
					//ik7++;
#if doubleintprecision == 1
					//printf("i==%lld, number=%lld\n",i,ik7);
#else
					//printf("i==%d, number=%d\n",i,ik7);
#endif
					
					//system("pause");
			     }
		     }

		     // Выделяем левую подобласть:
		     // Раскрашивает связный объём.
             // помечает всю связную окрестность связную с заданной нормалью в окрестности точки iP.
             // 0 - не помечено, 1 - помечено, 2 - запрет на пометку.
		     // Мы полностью отказались от рекурсии и создали полностью итеративный алгоритм. 24 апреля 2012.
             fill_marker_volume_iterat(iPleft, color_left, f);
			 //fill_marker_volume(iPleft, color_left, f);
		
			 integer icount_non_markerVolume=0;
		     for (integer i=0; i<f.maxelm+f.maxbound; ++i) {
			     if (color_left[i]==0) {
				    bOk=true;
					//printf("incomming\n");
					//system("pause");
					icount_non_markerVolume++;
#if doubleintprecision == 1
					/*if (i>f.maxelm) {
						printf("bound %lld maxelm+maxbound-1=%lld\n",i,f.maxelm+f.maxbound-1);
					}
					else {
						printf("internal %lld maxelm+maxbound-1=%lld\n",i,f.maxelm+f.maxbound-1);
						system("pause");
					}*/
#else
					/*if (i>f.maxelm) {
						printf("bound %d maxelm+maxbound-1=%d\n",i,f.maxelm+f.maxbound-1);
					}
					else {
						printf("internal %d maxelm+maxbound-1=%d\n",i,f.maxelm+f.maxbound-1);
						system("pause");
					}*/
#endif
					
				    //break; // Обязательно раскоментировать !!!!
			     }
		     }
			 bfirst=false;
#if doubleintprecision == 1
			 printf("icount_non_markerVolume=%lld\n", icount_non_markerVolume);
#else
			 printf("icount_non_markerVolume=%d\n", icount_non_markerVolume);
#endif
			 

			 /* bOk - остаётся false если граница сепаратора определена неполностью, а именно
			 граница разделителя имеет дырку через которую осуществлена закраска всей расчётной области.
			 В этом случае(осуществлена закраска всей расчётной области) требуется продолжить нахождение 
			 границы разделителя с тем же iso_value в окрестности данного расстояния от входной границы.
			 */

		}
		printf("4");
		//system("pause");
		// Выделяем правую подобласть:
        // Раскрашивает связный объём.
        // помечает всю связную окрестность связную с заданной нормалью в окрестности точки iP.
        // 0 - не помечено, 1 - помечено, 2 - запрет на пометку.
		// Мы полностью отказались от рекурсии и создали полностью итеративный алгоритм. 24 апреля 2012.
        fill_marker_volume_iterat(iPright, color_right, f);
		printf("5\n");
		//system("pause");

		/* Здесь ещё остаётся возможность неверного разделения:
		Случай неверного разделения определяется следующим образом:
		Все узлы минус узлы помеченные разделителем, минус узлы помеченные левым сепаратором минус узлы помеченные правым сепаратором
		не равно пусто значит разбиение произведено ошбочно. Это может быть ели сепаратор отрезал изолированный кусок расчётной области,
		алевый и правый маркеры находятся на одном куске разбиения отличного от отрезанного обособленного куска.
		*/
		integer *ifillcorect=new integer[f.maxelm+f.maxbound];
		for (integer i=0; i<f.maxelm+f.maxbound; ++i) {
			  ifillcorect[i]=2;
		}
		for (integer i=istart; i<=ifinish; ++i) {
		    integer iP=ifrontregulation[i];
		    ifillcorect[iP]=0;
		}
		for (integer i=0; i<f.maxelm+f.maxbound; ++i) {
			if (color[i]==1) {
			   ifillcorect[i]=1;
			}
			if (color_left[i]==1) {
			   ifillcorect[i]=1;
			}
			if (color_right[i]==1) {
			   ifillcorect[i]=1;
			}
		}
		// Наконец поверка.
		bool bfound0=false;
		for (integer i=0; i<f.maxelm+f.maxbound; ++i) {
			if (ifillcorect[i]==0) {
				// разделение некорректно !!!
				bfound0=true;
				break;
			}
		}
		if (!bfound0) {
			// разделение полностью корректно и рестарт ненужен.
			brestart=false;
		}
		if (ifillcorect != nullptr) {
			delete[] ifillcorect; // Освобождение памяти.
		}

		if (!brestart) {
		    // Финишная стадия: переупорядочивание.
		    int *iforig=nullptr;
		    iforig=new int[f.maxelm+f.maxbound];
		    doublereal *rfront_orig=nullptr;
		    rfront_orig=new doublereal[f.maxelm+f.maxbound];
		    for (integer i=0; i<f.maxelm+f.maxbound; ++i) {
			    iforig[i]=ifrontregulation[i]; // создаём копию ifrontregularization.
			    rfront_orig[i]=front[i];
		    }

		    // Препроцессинг для упорядочивания, чтобы выяснить  нужно он или нет.
		    // левая часть.
		    integer icl=0;
		    for (integer i=istart; i<=ifinish; ++i) {
			    integer iP=iforig[i];
			    if (color_left[iP]==1) {
			    	//ifrontregulation[istart+icl]=iP;
			    	// также нужно переупорядочить фронт.
				    //front[istart+icl]=rfront_orig[i];
			    	icl++;
			    }
		    }
		    
		    // правая часть.
            integer icr=0;
		    for (integer i=istart; i<=ifinish; ++i) {
			    integer iP=iforig[i];
			    if (color_right[iP]==1) {
			    	//ifrontregulation[b0.ileft_finish+1+icr]=iP;
				    // также нужно переупорядочить фронт.
			    	//front[b0.ileft_finish+1+icr]=rfront_orig[i];
			    	icr++;
			    }
		    }
			

		    // разделитель.
		    integer ics=0;
		    for (integer i=istart; i<=ifinish; ++i) {
		    	integer iP=iforig[i];
			    if (color[iP]==1) {
			    	//ifrontregulation[b0.iright_finish+1+ics]=iP;
			    	// также нужно переупорядочить фронт.
			     	//front[b0.iright_finish+1+ics]=rfront_orig[i];
			    	ics++;
			    }
		    }

			// 
			// 8 августа 2015 года.
			// этот код предназначен для улучшения качества разбиения.
			// проверка баланса: если нарушен баланс (левая и правая части отличаются более чем
			// на 20% ) то будут осуществлены подвижки границы сепаратора немного пошевелены 
			// до тех пора пока левая  правая части будут равны с точностью выше чем 20%.
			if (ilevel_info_nested_desection==1) {
				// Первое деление очень качественное т.к. от него многое зависит.
				if ((100.0*(icl)/(ifinish-istart+1)<=33.0)||(100.0*(icr)/(ifinish-istart+1)<=33.0)) {
			    irules_restart++;
				if (irules_restart<isize_rules_step) {
					if (irules_restart==isize_rules_step-1) {
						// делаем лучший вариант и заканчиваем.
#if doubleintprecision == 1
						printf("rules move best iteration %lld getready...\n", ibestvar);
#else
						printf("rules move best iteration %d getready...\n", ibestvar);
#endif
						
						acor=acorbestvar;
						brestart=true;
					}
					else {
		                  acor+=3.0*rules_step[irules_restart];
#if doubleintprecision == 1
						  printf("rules move iteration %lld getready...\n", irules_restart);
#else
						  printf("rules move iteration %d getready...\n", irules_restart);
#endif
				          
				          brestart=true;
				          // накопление информации о лучшем варианте.
				          if ((icl-icr)*(icl-icr)<minimumdisbalance) {
					         minimumdisbalance=((icl-icr)*(icl-icr));
					         ibestvar=irules_restart;
							 acorbestvar=acor;
				          }
					}
				}
				else {
					irules_restart=0;
					brestart=false;
				}
			}
			else {
                irules_restart=0; // всё успешно и улучшений не требуется.
			}
			}
			else {
			if ((100.0*(icl)/(ifinish-istart+1)<=33.0)||(100.0*(icr)/(ifinish-istart+1)<=33.0)) {
			    irules_restart++;
				if (irules_restart<isize_rules_step) {
					if (irules_restart==isize_rules_step-1) {
						// делаем лучший вариант и заканчиваем.
#if doubleintprecision == 1
						printf("rules move best iteration %lld getready...\n", ibestvar);
#else
						printf("rules move best iteration %d getready...\n", ibestvar);
#endif
						
						acor=acorbestvar;
						brestart=true;
					}
					else {
		                  acor+=3.0*rules_step[irules_restart];
#if doubleintprecision == 1
						  printf("rules move iteration %lld getready...\n", irules_restart);
#else
						  printf("rules move iteration %d getready...\n", irules_restart);
#endif
				          
				          brestart=true;
				          // накопление информации о лучшем варианте.
				          if ((icl-icr)*(icl-icr)<minimumdisbalance) {
					         minimumdisbalance=((icl-icr)*(icl-icr));
					         ibestvar=irules_restart;
							 acorbestvar=acor;
				          }
					}
				}
				else {
					irules_restart=0;
					brestart=false;
				}
			}
			else {
                irules_restart=0; // всё успешно и улучшений не требуется.
			}
			}

			if (!brestart) {

			b0.ileft_start=istart;
		    b0.ileft_finish=istart+icl-1;
		    b0.iright_start=b0.ileft_finish+1;
		    b0.iright_finish=b0.ileft_finish+1+icr-1;
		    b0.iseparate_start=b0.iright_finish+1;
		    b0.iseparate_finish=b0.iright_finish+1+ics-1;
#if doubleintprecision == 1
			printf("icl=%lld, icr=%lld, ics=%lld\n", icl, icr, ics);
#else
			printf("icl=%d, icr=%d, ics=%d\n", icl, icr, ics);
#endif
		    
			

		    //system("pause");
		    doublereal rsep=100.0*(b0.iseparate_finish-b0.iseparate_start+1)/(ifinish-istart+1);
		    doublereal rleft=100.0*(b0.ileft_finish-b0.ileft_start+1)/(ifinish-istart+1);
		    doublereal rright=100.0*(b0.iright_finish-b0.iright_start+1)/(ifinish-istart+1);

		

		    if (bseparate) {
		       if ((rleft>rsep)&&(rright>rsep)) {
			       b0.active=true;
		           // Наконец переупорядочивание.
		           // левая часть.
		           icl=0;
		           for (integer i=istart; i<=ifinish; ++i) {
			           int iP=iforig[i];
			           if (color_left[iP]==1) {
				            ifrontregulation[istart+icl]=iP;
				            // также нужно переупорядочить фронт.
				            //front[istart+icl]=rfront_orig[i];
				            front[istart+icl]=front_orig[iP];
				            icl++;
			           }
		           }
		           b0.ileft_start=istart;
		           b0.ileft_finish=istart+icl-1;
		           // правая часть.
                   icr=0;
		           for (integer i=istart; i<=ifinish; ++i) {
			           int iP=iforig[i];
			           if (color_right[iP]==1) {
				           ifrontregulation[b0.ileft_finish+1+icr]=iP;
				           // также нужно переупорядочить фронт.
				           //front[b0.ileft_finish+1+icr]=rfront_orig[i];
				           front[b0.ileft_finish+1+icr]=front_orig[iP];
				           icr++;
			           }
		           }
		           b0.iright_start=b0.ileft_finish+1;
		           b0.iright_finish=b0.ileft_finish+1+icr-1;

		           // разделитель.
		           ics=0;
		           for (integer i=istart; i<=ifinish; ++i) {
			           int iP=iforig[i];
			           if (color[iP]==1) {
				          ifrontregulation[b0.iright_finish+1+ics]=iP;
				          // также нужно переупорядочить фронт.
				          //front[b0.iright_finish+1+ics]=rfront_orig[i];
				          front[b0.iright_finish+1+ics]=front_orig[iP];
				          ics++;
			           }
		           }
		           b0.iseparate_start=b0.iright_finish+1;
		           b0.iseparate_finish=b0.iright_finish+1+ics-1;
		      }
		   }

			}

		   if (rfront_orig != nullptr) {
				delete[] rfront_orig;
		   }
		   if (iforig != nullptr) {
			   delete[] iforig;
		   }
		}
		else {
			// Подготовка к рестарту !!!
			acor+=mul*0.005; // сдвигаем разделитель на 1%.
			printf("separator is moving plus length 1%% \n");
			printf("Please, press any key to continue calculation...\n");
			//system("pause");
		}

		// метка для рестарта с другим значением acor.

		if (color != nullptr) {
			delete[] color;
		}
		if (color_left != nullptr) {
			delete[] color_left;
		}
		if (color_right != nullptr) {
			delete[] color_right;
		}
		
		return 0;

		/*

		     // Внешний цикл нужен для улучшения характеристик разбиения.
		     // Всего три улучшающих подхода.
		     for (integer istep=0; istep<3; istep++) {
		        // Вычисляем корректирующую поправку разбиения.
		        if (rleft<rright) {
			        acor=fmin(50.0-rleft,rright-50.0);
			        acor=fabs(acor)/100.0;
		        }
		        else {
			        acor=fmin(50.0-rright,rleft-50.0);
			        acor=-fabs(acor)/100.0;
		        }
             }
		*/
} // my_separator5

void my_separator(FLOW &f, int* &ifrontregulation, const integer istart, const integer ifinish, integer iNorm,
					  doublereal* &front, doublereal* &front_orig,  PARBOUND &b0, bool bseparate, doublereal &acor,
					  bool &brestart, TEMPER &t, FLOW* &fglobal, integer flow_interior, doublereal mul, 
					  integer &irules_restart, doublereal *rules_step, integer isize_rules_step,
					  integer &ibestvar, bool &bbestvar, integer &minimumdisbalance, doublereal &acorbestvar)  {

						  if (1) {
							  // Этот факт проверен 15 августа 2015.
							  // Если сделать три клетки то потоки перекрываются
							  // и растут накладные расходы на организацию распараллеливания.
							  // Вообще непонятно, может и три клетки нормально. 



							  // Ширина разделителя 5 клеток,
							  // возможно это нужно для корректной работы всего что основано
							  // на ilu2 включая получение самого ilu2.
							  my_separator5(f, ifrontregulation, istart, ifinish, iNorm,
					                        front, front_orig,  b0, bseparate, acor,
					                        brestart, t, fglobal, flow_interior, mul, 
					                        irules_restart, rules_step, isize_rules_step,
					                        ibestvar, bbestvar, minimumdisbalance, acorbestvar) ;
						  }
						  else {
							  // Ширина разделителя 3 клетки.
							  my_separator3(f, ifrontregulation, istart, ifinish, iNorm,
					                        front, front_orig,  b0, bseparate, acor,
					                        brestart, t, fglobal, flow_interior, mul, 
					                        irules_restart, rules_step, isize_rules_step,
					                        ibestvar, bbestvar, minimumdisbalance, acorbestvar) ;
						  }

}

// разбивает на два потока.
void separate2(PARBOUND &b0, doublereal* front, doublereal* front_orig, FLOW &f,
			   int *ifrontregulation,const doublereal dmax, integer istart, integer ifinish,
			   TEMPER &t, FLOW* &fglobal, integer flow_interior) {

        // на начальном этапе istart==0, ifinish==f.maxelm+f.maxbound-1;
		bool brestart=true;
		doublereal acor=0.0; // без коррекции.
		doublereal mul=1.0; 
		integer icountadapt=0;
		
        integer irules_restart=0;
		integer isize_rules_step=201; //23
		//doublereal rules_step[7]={0.0,0.005,-0.01,0.02,-0.03,0.045,-0.06};
		doublereal *rules_step=nullptr;
		rules_step=new doublereal[isize_rules_step];
		rules_step[0]=0.0;
		doublereal rpattern=0.001;
		for (integer irs=1; irs<isize_rules_step; irs++) {
			if (irs%2==0) {
				rules_step[irs]=-rpattern;
			}
			else {
			    rules_step[irs]=rpattern;
			}
			rpattern+=0.001;
		}
		//{0.0,0.001,-0.002,0.003,-0.004,0.005,-0.006,0.007,-0.008,0.009,-0.01,0.011,-0.012,0.013,-0.014,0.015,-0.016,0.017,-0.018,0.019,-0.02,0.021,-0.022};
		// 9 августа 2015 
		// Мы запомним лучшее разбиение из isize_rules_step-1 попыток.
		integer ibestvar=0;
		bool bbestvar=false;
		integer minimumdisbalance=0;
#if doubleintprecision == 1
			minimumdisbalance = 140737488355326;//48 бит
#else
			minimumdisbalance = 2147483646;//32 бита
#endif
		doublereal acorbestvar=0.0;

		while ((brestart)&&(icountadapt<1000)) {

			if (acor>0.25) {
				acor=0.0;
				mul=-1.0;
				// теперь будем двигаться вниз.
			}
			if (acor<-0.25) {
				acor=0.0;
				mul=1.0;
				// теперь будем двигаться вверх.
			}

	        // пополам делим только в том случае если n==2 или n==4.
            integer igold=istart+static_cast<integer>((0.5+acor)*(ifinish-istart));
		    // Следующее деление более предпочтительно !:
		    doublereal l1=front[igold]; // делим пополам по количеству контрольных объёмов. // именно front всё верно Ок.
		    printf("separator is: %e\n",l1);
		    //system("pause");
		    // Вычисление среднеобъёмного линейного размера контрольного объёма.
		    doublereal avgdx=0.0;
		    integer ic=0;
			doublereal delta=2.0;
			if (!bunion_8_08_2015) {
		        for (integer i=istart; i<=ifinish; ++i) { // Только внутренние КО !
			        integer iP=ifrontregulation[i];
			        if (iP<f.maxelm) {
			            // вычисление размеров текущего контрольного объёма:
	                    doublereal dx=0.0, dy=0.0, dz=0.0;// объём текущего контрольного объёма
	                    volume3D(iP, f.nvtx, f.pa, dx, dy, dz);
			            avgdx+=dx*dy*dz;
			            ic++; // счётчик числа элементов.
			        }
		        }
		        avgdx/= static_cast<doublereal>(ic);
		        avgdx=exp(1.0/3.0*log(avgdx)); // средний линейный размер контрольного объёма по всей расчётной области.
		        printf(" avgVol=%e, \n",avgdx);
		        delta=2.0*avgdx;
		        // Вычисление среднеобъёмного линейного размера контрольного объёма в окрестности сепаратора.
		        avgdx=0.0; // уточнение avgdx
		        integer icount2=0;
		        for (integer i=istart; i<=ifinish; ++i) { // Только внутренние КО !
			        integer iP=ifrontregulation[i];
			        if (iP<f.maxelm) {
			            if (fabs(front_orig[iP]-l1)<delta) {
				            // вычисление размеров текущего контрольного объёма:
	                        doublereal dx=0.0, dy=0.0, dz=0.0;// объём текущего контрольного объёма
	                        volume3D(iP, f.nvtx, f.pa, dx, dy, dz);
			                avgdx+=dx*dy*dz;
				            icount2++;
			            }
			        }
		        }
		        avgdx/= static_cast<doublereal>(icount2);
		        avgdx=exp(1.0/3.0*log(avgdx)); // Средний линейный размер контрольного объёма в окрестности сепаратора.
		        printf(" avgVol separator=%e\n",avgdx);
	            icount2=0;
		        delta=2.0*avgdx;
		
			}
			else {
				// для фронта все объёмы одинаковые единичные !!!
				delta=1.1;
			}

		    // Более точное определение элементов входящих в разделитель.
		    // Вычисление Нормали к поверхности разделителя.
		    integer iX=0, iY=0, iZ=0; // количество граничных элементов по каждой из оси.
		    // Нормаль будет присвоена тому элементу у которого минимальное значение iG.
		    for (integer i=istart; i<=ifinish; ++i) { // Только внутренние КО !
			    integer iP=ifrontregulation[i];
			    if (iP<f.maxelm) {
		
			        if (fabs(front_orig[iP]-l1)<delta) {
				         integer iE, iN, iT, iW, iS, iB; // номера соседних контрольных объёмов
						 iE = f.neighbors_for_the_internal_node[E_SIDE][0][iP]; iN = f.neighbors_for_the_internal_node[N_SIDE][0][iP]; iT = f.neighbors_for_the_internal_node[T_SIDE][0][iP]; iW = f.neighbors_for_the_internal_node[W_SIDE][0][iP]; iS = f.neighbors_for_the_internal_node[S_SIDE][0][iP]; iB = f.neighbors_for_the_internal_node[B_SIDE][0][iP];
	
	                     // Внутренний КО.	

	                     // Если с одной из сторон стоит граница расчётной области
	                     // то соответствующая переменная равна true
	   
	                     if (iE>=f.maxelm) iX++;
	                     if (iN>=f.maxelm) iY++;
	                     if (iT>=f.maxelm) iZ++;
                         if (iW>=f.maxelm) iX++;
	                     if (iS>=f.maxelm) iY++;
	                     if (iB>=f.maxelm) iZ++;
			        }
			    }
		    }
		    //const integer xnorm=0;
		    //const integer ynorm=1;
		    //const integer znorm=2;
		    bool undefined=false;
		    integer iNorm=-1; // Направление перпендикулярное сечению разделителя.
			bool boldstyle=false;
			if (boldstyle) {
				//  этот метод датируется 2013 годом.

				// по какому направлению среди x, y или z меньше всего граничных узлов то направление
				// и является направлением нормали.
				if ((iX<=iY)&&(iX<=iZ)) {
					iNorm=xnorm;
				}
				if ((iY<=iX)&&(iY<=iZ)) {
					iNorm=ynorm;
				}
				if ((iZ<=iY)&&(iZ<=iX)) iNorm=znorm;
			}
			else {
				// этот метод 2015 года.

				// Мы берём точку А принадлежащую сепаратору (области разделения двух подобластей).
				// Наша цель на данном этапе найти ориентацию нормали к сепаратору.
				// Нормаль будет коллинеарна наиболее длинной линии построенной из точки А до границ 
				// расчётной области в пределах делимой в данный раз геометрической подобласти.

				/*{
				// 4 августа 2015.
				integer igoldl=istart+static_cast<integer>((0.5+acor)*(ifinish-istart));
		        integer iPseparatorl=ifrontregulation[igoldl];
		        // В качестве iPseparatorl нам нужен строго внутренний контрольный объём.
		        if (iPseparatorl>=f.maxelm) {
			        for (integer i=istart; i<=ifinish; ++i) {
			        integer iPl=ifrontregulation[i];
			        if (iPl<f.maxelm) {
				         integer iE, iN, iT, iW, iS, iB; // номера соседних контрольных объёмов
	                     iE=f.neighbors_for_the_internal_node[ESIDE][iPl].iNODE1; iN=f.neighbors_for_the_internal_node[NSIDE][iPl].iNODE1; iT=f.neighbors_for_the_internal_node[TSIDE][iPl].iNODE1; iW=f.neighbors_for_the_internal_node[WSIDE][iPl].iNODE1; iS=f.neighbors_for_the_internal_node[SSIDE][iPl].iNODE1; iB=f.neighbors_for_the_internal_node[BSIDE][iPl].iNODE1;

				         if ((iPseparatorl==iE)||(iPseparatorl==iW)||(iPseparatorl==iN)||(iPseparatorl==iS)||(iPseparatorl==iT)||(iPseparatorl==iB)) {
					         iPseparatorl=iPl; // именно внутренний контрольный объём.
					         break;
				         }
			         }
			     }
		        }
				integer ilx=1;
				integer iPsearch=f.neighbors_for_the_internal_node[ESIDE][iPseparatorl].iNODE1;
				while ((f.neighbors_for_the_internal_node[ESIDE][0][iPsearch]<f.maxelm)&&(front[istart]<front_orig[iPsearch])&&(front_orig[iPsearch]<front[ifinish])) { 
					iPsearch=f.neighbors_for_the_internal_node[ESIDE][0][iPsearch];
					ilx++;
				}
				 iPsearch=f.neighbors_for_the_internal_node[WSIDE][iPseparatorl].iNODE1;
				while ((f.neighbors_for_the_internal_node[WSIDE][0][iPsearch]<f.maxelm)&&(front[istart]<front_orig[iPsearch])&&(front_orig[iPsearch]<front[ifinish])) { 
					iPsearch=f.neighbors_for_the_internal_node[WSIDE][0][iPsearch];
					ilx++;
				}
				integer ily=1;
				
                    iPsearch=f.neighbors_for_the_internal_node[NSIDE][iPseparatorl].iNODE1;
				    while ((f.neighbors_for_the_internal_node[NSIDE][0][iPsearch]<f.maxelm)&&(front[istart]<front_orig[iPsearch])&&(front_orig[iPsearch]<front[ifinish])) { 
					   iPsearch=f.neighbors_for_the_internal_node[NSIDE][0][iPsearch];
					   ily++;
				    }
				     iPsearch=f.neighbors_for_the_internal_node[SSIDE][iPseparatorl].iNODE1;
				    while ((f.neighbors_for_the_internal_node[SSIDE][0][iPsearch]<f.maxelm)&&(front[istart]<front_orig[iPsearch])&&(front_orig[iPsearch]<front[ifinish])) { 
					   iPsearch=f.neighbors_for_the_internal_node[SSIDE][0][iPsearch];
					   ily++;
				    }
					
					
				
				integer ilz=1;
				
                     iPsearch=f.neighbors_for_the_internal_node[TSIDE][iPseparatorl].iNODE1;
				    while ((f.neighbors_for_the_internal_node[TSIDE][0][iPsearch]<f.maxelm)&&(front[istart]<front_orig[iPsearch])&&(front_orig[iPsearch]<front[ifinish])) { 
					   iPsearch=f.neighbors_for_the_internal_node[TSIDE][0][iPsearch];
					   ilz++;
				    }
				     iPsearch=f.neighbors_for_the_internal_node[BSIDE][iPseparatorl].iNODE1;
				    while ((f.neighbors_for_the_internal_node[BSIDE][0][iPsearch]<f.maxelm)&&(front[istart]<front_orig[iPsearch])&&(front_orig[iPsearch]<front[ifinish])) { 
					   iPsearch=f.neighbors_for_the_internal_node[BSIDE][0][iPsearch];
					   ilz++;
				    }
					if ((ilz>=ilx)&&(ilz>=ily)) iNorm=znorm;
					if ((ily>=ilx)&&(ily>=ilz)) iNorm=ynorm;
					if ((ilx>=ilz)&&(ilx>=ily)) iNorm=xnorm;
					
			
			}*/

				{ // 8 августа 2015
					bool* bactivcolor = nullptr;
					bactivcolor = new bool[f.maxelm + f.maxbound];
					if (bactivcolor == nullptr) {
						// недостаточно памяти на данном оборудовании.
						printf("Problem: not enough memory on your equipment for ent constr struct...\n");
						printf("Please any key to exit...\n");
						//system("pause");
						system("pause");
						exit(1);
					}
					else {

						if (bactivcolor != nullptr) {

							for (integer i = 0; i < f.maxelm + f.maxbound; ++i) bactivcolor[i] = false;
							for (integer i = istart; i <= ifinish; ++i) {
								integer iP = ifrontregulation[i];
								bactivcolor[iP] = true;
							}


							// 4 августа 2015.
							integer igoldl = istart + static_cast<integer>((0.5 + acor)*(ifinish - istart));
							integer iPseparatorl = ifrontregulation[igoldl];
							// 8 августа 2015.
							iPseparatorl = uluchenie_katchestva_razbieniq(f, ifrontregulation, igoldl, istart, ifinish, front, iPseparatorl);


							integer ilx = 1;
							integer iPsearch = iPseparatorl;
							while ((f.neighbors_for_the_internal_node[E_SIDE][0][iPsearch] < f.maxelm) && (bactivcolor[f.neighbors_for_the_internal_node[E_SIDE][0][iPsearch]])) {
								iPsearch = f.neighbors_for_the_internal_node[E_SIDE][0][iPsearch];
								ilx++;
							}
							if (f.neighbors_for_the_internal_node[E_SIDE][0][iPsearch] < f.maxelm + f.maxbound) {
								if ((f.neighbors_for_the_internal_node[E_SIDE][0][iPsearch] >= f.maxelm) && (bactivcolor[f.neighbors_for_the_internal_node[E_SIDE][0][iPsearch]])) {
									ilx++; // учёт граничного узла.
								}
							}

							iPsearch = iPseparatorl;
							while ((f.neighbors_for_the_internal_node[W_SIDE][0][iPsearch] < f.maxelm) && (bactivcolor[f.neighbors_for_the_internal_node[W_SIDE][0][iPsearch]])) {
								iPsearch = f.neighbors_for_the_internal_node[W_SIDE][0][iPsearch];
								ilx++;
							}
							if (f.neighbors_for_the_internal_node[W_SIDE][0][iPsearch] < f.maxelm + f.maxbound) {
								if ((f.neighbors_for_the_internal_node[W_SIDE][0][iPsearch] >= f.maxelm) && (bactivcolor[f.neighbors_for_the_internal_node[W_SIDE][0][iPsearch]])) {
									ilx++; // учёт граничного узла.
								}
							}

							integer ily = 1;

							iPsearch = iPseparatorl;
							while ((f.neighbors_for_the_internal_node[N_SIDE][0][iPsearch] < f.maxelm) && (bactivcolor[f.neighbors_for_the_internal_node[N_SIDE][0][iPsearch]])) {
								iPsearch = f.neighbors_for_the_internal_node[N_SIDE][0][iPsearch];
								ily++;
							}
							if (f.neighbors_for_the_internal_node[N_SIDE][0][iPsearch] < f.maxelm + f.maxbound) {
								if ((f.neighbors_for_the_internal_node[N_SIDE][0][iPsearch] >= f.maxelm) && (bactivcolor[f.neighbors_for_the_internal_node[N_SIDE][0][iPsearch]])) {
									ily++; // учёт граничного узла.
								}
							}
							iPsearch = iPseparatorl;
							while ((f.neighbors_for_the_internal_node[S_SIDE][0][iPsearch] < f.maxelm) && (bactivcolor[f.neighbors_for_the_internal_node[S_SIDE][0][iPsearch]])) {
								iPsearch = f.neighbors_for_the_internal_node[S_SIDE][0][iPsearch];
								ily++;
							}
							if (f.neighbors_for_the_internal_node[S_SIDE][0][iPsearch] < f.maxelm + f.maxbound) {
								if ((f.neighbors_for_the_internal_node[S_SIDE][0][iPsearch] >= f.maxelm) && (bactivcolor[f.neighbors_for_the_internal_node[S_SIDE][0][iPsearch]])) {
									ily++; // учёт граничного узла.
								}
							}


							integer ilz = 1;

							iPsearch = iPseparatorl;
							while ((f.neighbors_for_the_internal_node[T_SIDE][0][iPsearch] < f.maxelm) && (bactivcolor[f.neighbors_for_the_internal_node[T_SIDE][0][iPsearch]])) {
								iPsearch = f.neighbors_for_the_internal_node[T_SIDE][0][iPsearch];
								ilz++;
							}
							if (f.neighbors_for_the_internal_node[T_SIDE][0][iPsearch] < f.maxelm + f.maxbound) {
								if ((f.neighbors_for_the_internal_node[T_SIDE][0][iPsearch] >= f.maxelm) && (bactivcolor[f.neighbors_for_the_internal_node[T_SIDE][0][iPsearch]])) {
									ilz++; // учёт граничного узла.
								}
							}
							iPsearch = iPseparatorl;
							while ((f.neighbors_for_the_internal_node[B_SIDE][0][iPsearch] < f.maxelm) && (bactivcolor[f.neighbors_for_the_internal_node[B_SIDE][0][iPsearch]])) {
								iPsearch = f.neighbors_for_the_internal_node[B_SIDE][0][iPsearch];
								ilz++;
							}
							if (f.neighbors_for_the_internal_node[B_SIDE][0][iPsearch] < f.maxelm + f.maxbound) {
								if ((f.neighbors_for_the_internal_node[B_SIDE][0][iPsearch] >= f.maxelm) && (bactivcolor[f.neighbors_for_the_internal_node[B_SIDE][0][iPsearch]])) {
									ilz++; // учёт граничного узла.
								}
							}
							if ((ilz >= ilx) && (ilz >= ily)) iNorm = znorm;
							if ((ily >= ilx) && (ily >= ilz)) iNorm = ynorm;
							if ((ilx >= ilz) && (ilx >= ily)) iNorm = xnorm;
							printf("new normal detection diagnostic mesage:\n");
#if doubleintprecision == 1
							printf("ilx=%lld, ily=%lld, ilz=%lld\n", ilx, ily, ilz);
#else
							printf("ilx=%d, ily=%d, ilz=%d\n", ilx, ily, ilz);
#endif
							


							if (bactivcolor != nullptr) {
								delete[] bactivcolor;
							}
						}
					}

				
				}

			}
		    
				
		    
#if doubleintprecision == 1
			printf("iX=%lld, iY=%lld, iZ=%lld, normal is:", iX, iY, iZ);
#else
			printf("iX=%d, iY=%d, iZ=%d, normal is:", iX, iY, iZ);
#endif
		    
		    switch(iNorm) {
		       case xnorm: printf(" X\n");
			           if (iX!=0) {
						   undefined=true;
					   }
			           break;
		       case ynorm: printf(" Y\n");
			           if (iY!=0) {
						   undefined=true;
					   }
			           break;
		       case znorm: printf(" Z\n");
			           if (iZ!=0) {
						   undefined=true;
					   } 
			           break;
		       default: printf("Error identification norm \n");
				         printf("Please, press any key to halt calculation process...\n");
			             //system("pause"); 
						 system("pause");
						 exit(1);
			             break;
		    }
		    // Случай undefined==true: 
		    /* сложная граница разделителя не лежащая в одной плоскости. 
		       Этот случай может често встречатьс в реальной геометрии и требует особой обработки.
		    */
		    if (undefined) {
			     printf("Warning! Undifined situation found...\n");
		    }

		    // вычисление линейного размера контрольного объёма в окрестности 
		    // сепаратора в направлении перпендикулярном плоскости сепаратора.
		    avgdx=0.0; // уточнение avgdx теперь с учётом выделенного направления.
		    integer icount2=0;
			if (!bunion_8_08_2015) {
		    for (integer i=istart; i<=ifinish; ++i) { // Только внутренние КО !
			    integer iP=ifrontregulation[i];
			    if (iP<f.maxelm) {
		
			        if (fabs(front_orig[iP]-l1)<delta) {
				        // вычисление размеров текущего контрольного объёма:
	                    doublereal dx=0.0, dy=0.0, dz=0.0;// объём текущего контрольного объёма
	                    volume3D(iP, f.nvtx, f.pa, dx, dy, dz);
				        switch (iNorm) {
				           case xnorm: avgdx=fmax(avgdx,dx); break; 
				           case ynorm: avgdx=fmax(avgdx,dy); break;
				           case znorm: avgdx=fmax(avgdx,dz); break;
				           default: avgdx+=dx*dy*dz;
					                break;
				        }			    
				        icount2++;
			        }
			    }
		    }
		    if ((iNorm!=xnorm)&&(iNorm!=ynorm)&&(iNorm!=znorm)) {
		        avgdx/= static_cast<doublereal>(icount2);
		        avgdx=exp(1.0/3.0*log(avgdx)); // Средний линейный размер контрольного объёма в окрестности сепаратора.
		    }
		    printf(" avg direct separator=%e\n",avgdx);
	        icount2=0;
		    delta=2.0*avgdx;
			}
			else {
				// Для nested desection важна не неравномерность сетки и геометрические размеры контрольного
				// объёма, а количество контрольных объёмов и их связи друг с другом. 
				// В версии после 8 августа 2015 года все объёмы есть единичные кубы все рёбра которых имеют длину единица.
				avgdx=1.0;
				delta=1.1;
			}

		    //system("pause");

		    // Внимание ! знание нормали позволит более точно определить элементы сепаратора !.

		
		    integer imin=ifinish+1;
		
		    // Здесь будет код my_separator самое последнее достижение в сепараторах.
		    //  20 апреля 2013 года.
		    bool bseparate=true; // делим.
		
		    // Параметры t, fglobal, flow_interior нужны для графической визуализации.
		    my_separator(f, ifrontregulation, istart, ifinish, iNorm, front, front_orig, b0, bseparate, 
				acor, brestart, t,fglobal,flow_interior,mul,irules_restart,rules_step,isize_rules_step,
				ibestvar, bbestvar, minimumdisbalance, acorbestvar);
			
			
			
			if (!brestart) {
				// если мы не рестартим, то это финишное разбиение.
		       printf("disbalance: left: %3.0f %%, right: %3.0f %%, separation: %3.0f %%\n",100.0*(b0.ileft_finish-b0.ileft_start+1)/(ifinish-istart+1),100.0*(b0.iright_finish-b0.iright_start+1)/(ifinish-istart+1),100.0*(b0.iseparate_finish-b0.iseparate_start+1)/(ifinish-istart+1));
			}
            // упорядочим элементы в соответствии с возрастанием элементов в массиве front.
		    //if (b0.active) {
		       // Сортировка ненужна так как массивы и так уже упорядочены.
		       //printf("Quick sort start...\n");
		       // QuickSortnd(ifrontregulation,  front,  b0.ileft_start, b0.ileft_finish); // Неэкономично, т.к. требует много флопов.
		       // QuickSortnd(ifrontregulation,  front,  b0.iright_start, b0.iright_finish); 
		       // QuickSortnd(ifrontregulation,  front,  b0.iseparate_start, b0.iseparate_finish); 
		       // printf("Quick sort finish...\n");
		    //}
		     


		    // Приближённый алгоритм.
		    // Этот алгоритм скорее всего никогда не будет применён !
		    // т.к. эта ситуация по-видимому не может встретиться.
		    if (iNorm==-1) {
			    // Криволинейный разделитель !!!.
		        // Это приближённое определение элементов разделителя,
		        // нужно более точно определить какие элементы попадут в разделитель.
		        for (integer i=istart; i<=ifinish; ++i) {
			        if (fabs(front[i]-l1)<delta) { // здесь именно переупорядоченный фронт.
				        front[i]=dmax+2.0; //  возможно потребуется увеличить это значение, но смысл в том что эти значения записываются в конец массива inumerate. 
				        icount2++;
				        imin=(imin < i ? imin : i);
			        }
		        }
		    }

		    // Если потребуется осуществим перезапуск.
			icountadapt++;
		} // brestart

		

		// Перезапуски завершены.
		// Работа по разбиению на два закончена.

		if (icountadapt==1000) {
			printf("partition found for 1000 passes...\n");
			printf("Please, press any key to exit...\n");
			//system("pause");
			system("pause");
			exit(1);
		}

		delete[] rules_step;

} // separate2

// Визуализация разбиения !
void visible_partition(FLOW* &fglobal, FLOW &f, integer flow_interior, TEMPER &t, int* &ifrontregulation, PARDATA &nd) {
	// Раскраска для контроля путём визуализации разбиения.
	doublereal* parcolor = nullptr;
	parcolor=new doublereal[f.maxelm + f.maxbound];

        if (nd.ncore==2) {
			if (nd.b0.active) {
			   for (integer i=nd.b0.ileft_start; i<=nd.b0.ileft_finish; ++i) {
				   integer iP=ifrontregulation[i];
                   parcolor[iP]=1.0;
			   }
			   for (integer i=nd.b0.iright_start; i<=nd.b0.iright_finish; ++i) {
				   integer iP=ifrontregulation[i];
                   parcolor[iP]=2.0;
			   }
			   for (integer i=nd.b0.iseparate_start; i<=nd.b0.iseparate_finish; ++i) {
				   integer iP=ifrontregulation[i];
                   parcolor[iP]=3.0;
			   }
			}
			else {
				for (integer i=0; i<f.maxelm+f.maxbound; ++i) {
				   integer iP=ifrontregulation[i];
                   parcolor[iP]=1.0;
			     }
			}
        }
		if (nd.ncore==4) {
			doublereal rc=1.0;
			if (nd.b0.active) {
			if (nd.b00.active) {
			    for (integer i=nd.b00.ileft_start; i<=nd.b00.ileft_finish; ++i) {
				    integer iP=ifrontregulation[i];
                    parcolor[iP]=rc;
			    }
				rc+=1.0;
			    for (integer i=nd.b00.iright_start; i<=nd.b00.iright_finish; ++i) {
				    integer iP=ifrontregulation[i];
                    parcolor[iP]=rc;
			    }
				rc+=1.0;
			}
			else {
				for (integer i=nd.b0.ileft_start; i<=nd.b0.ileft_finish; ++i) {
				    integer iP=ifrontregulation[i];
                    parcolor[iP]=rc;
			    }
				rc+=1.0;
			}
			if (nd.b01.active) {
			   for (integer i=nd.b01.ileft_start; i<=nd.b01.ileft_finish; ++i) {
				   integer iP=ifrontregulation[i];
                   parcolor[iP]=rc;
			   }
			   rc+=1.0;
			   for (integer i=nd.b01.iright_start; i<=nd.b01.iright_finish; ++i) {
				   integer iP=ifrontregulation[i];
                   parcolor[iP]=rc;
			   }
			   rc+=1.0;
			}
			else {
				 for (integer i=nd.b0.iright_start; i<=nd.b0.iright_finish; ++i) {
				   integer iP=ifrontregulation[i];
                   parcolor[iP]=rc;
			   }
			   rc+=1.0;
			}
			}
			else {
				for (integer i=0; i<f.maxelm+f.maxbound; ++i) {
					integer iP=ifrontregulation[i];
                    parcolor[iP]=rc;
				}
				rc+=1.0;
			}
			if (nd.b00.active) {
			   for (integer i=nd.b00.iseparate_start; i<=nd.b00.iseparate_finish; ++i) {
				   integer iP=ifrontregulation[i];
                   parcolor[iP]=rc;
			   }
			   rc+=1.0;
			}
			if (nd.b01.active) {
			   for (integer i=nd.b01.iseparate_start; i<=nd.b01.iseparate_finish; ++i) {
				   integer iP=ifrontregulation[i];
                   parcolor[iP]=rc;
			   }
			   rc+=1.0;
			}
			if (nd.b0.active) {
		    	for (integer i=nd.b0.iseparate_start; i<=nd.b0.iseparate_finish; ++i) {
			    	integer iP=ifrontregulation[i];
                    parcolor[iP]=rc;
			    }
			}
        }
		if (nd.ncore==8) {
			doublereal rc=1.0;
			if (nd.b0.active) {
			if (nd.b00.active) {
			   if (nd.b000.active) {
			      for (integer i=nd.b000.ileft_start; i<=nd.b000.ileft_finish; ++i) {
				      integer iP=ifrontregulation[i];
                      parcolor[iP]=rc;
			      }
			      rc+=1.0;
			      for (integer i=nd.b000.iright_start; i<=nd.b000.iright_finish; ++i) {
				      integer iP=ifrontregulation[i];
                      parcolor[iP]=rc;
			      }
			      rc+=1.0;
			   } 
			   else {
				   for (integer i=nd.b00.ileft_start; i<=nd.b00.ileft_finish; ++i) {
				       integer iP=ifrontregulation[i];
                       parcolor[iP]=rc;
			       }
				   rc+=1.0;
			   }
			   if (nd.b001.active) {
			      for (integer i=nd.b001.ileft_start; i<=nd.b001.ileft_finish; ++i) {
				      integer iP=ifrontregulation[i];
                      parcolor[iP]=rc;
			      }
			      rc+=1.0;
			      for (integer i=nd.b001.iright_start; i<=nd.b001.iright_finish; ++i) {
				      integer iP=ifrontregulation[i];
                      parcolor[iP]=rc;
			      }
			      rc+=1.0;
			   }
			   else {
				   for (integer i=nd.b00.iright_start; i<=nd.b00.iright_finish; ++i) {
				       integer iP=ifrontregulation[i];
                       parcolor[iP]=rc;
			       }
				   rc+=1.0;
			   }
			}
			else {
				for (integer i=nd.b0.ileft_start; i<=nd.b0.ileft_finish; ++i) {
				    integer iP=ifrontregulation[i];
                    parcolor[iP]=rc;
			    }
				rc+=1.0;
			}
			if (nd.b01.active) {
			   if (nd.b010.active) {
			      for (integer i=nd.b010.ileft_start; i<=nd.b010.ileft_finish; ++i) {
				      integer iP=ifrontregulation[i];
                      parcolor[iP]=rc;
			      }
			      rc+=1.0;
			      for (integer i=nd.b010.iright_start; i<=nd.b010.iright_finish; ++i) {
				      integer iP=ifrontregulation[i];
                      parcolor[iP]=rc;
			      }
			      rc+=1.0;
			   }
			   else {
				   for (integer i=nd.b01.ileft_start; i<=nd.b01.ileft_finish; ++i) {
				       integer iP=ifrontregulation[i];
                       parcolor[iP]=rc;
			       }
				   rc+=1.0;
			   }
			   if (nd.b011.active) {
			      for (integer i=nd.b011.ileft_start; i<=nd.b011.ileft_finish; ++i) {
				     integer iP=ifrontregulation[i];
                     parcolor[iP]=rc;
			      }
			      rc+=1.0;
			      for (integer i=nd.b011.iright_start; i<=nd.b011.iright_finish; ++i) {
				     integer iP=ifrontregulation[i];
                     parcolor[iP]=rc;
			      }
			      rc+=1.0;
			   }
			   else {
				   for (integer i=nd.b01.iright_start; i<=nd.b01.iright_finish; ++i) {
				       integer iP=ifrontregulation[i];
                       parcolor[iP]=rc;
			       }
				   rc+=1.0;
			   }
			}
			else {
				for (integer i=nd.b0.iright_start; i<=nd.b0.iright_finish; ++i) {
				   integer iP=ifrontregulation[i];
                   parcolor[iP]=rc;
			    }
			    rc+=1.0;
			}
			}
			else {
				for (integer i=0; i<f.maxelm+f.maxbound; ++i) {
					integer iP=ifrontregulation[i];
                    parcolor[iP]=rc;
				}
				rc+=1.0;
			}
			if (nd.b000.active) {
			   for (integer i=nd.b000.iseparate_start; i<=nd.b000.iseparate_finish; ++i) {
				   integer iP=ifrontregulation[i];
                   parcolor[iP]=rc;
			   }
			   rc+=1.0;
			}
			if (nd.b001.active) {
			   for (integer i=nd.b001.iseparate_start; i<=nd.b001.iseparate_finish; ++i) {
				   integer iP=ifrontregulation[i];
                   parcolor[iP]=rc;
			   }
			   rc+=1.0;
			}
			if (nd.b010.active) {
			   for (integer i=nd.b010.iseparate_start; i<=nd.b010.iseparate_finish; ++i) {
				   integer iP=ifrontregulation[i];
                   parcolor[iP]=rc;
			   }
			   rc+=1.0;
			}
			if (nd.b011.active) {
			   for (integer i=nd.b011.iseparate_start; i<=nd.b011.iseparate_finish; ++i) {
				   integer iP=ifrontregulation[i];
                   parcolor[iP]=rc;
			   }
			   rc+=1.0;
			}
			if (nd.b00.active) {
			   for (integer i=nd.b00.iseparate_start; i<=nd.b00.iseparate_finish; ++i) {
				   integer iP=ifrontregulation[i];
                   parcolor[iP]=rc;
		 	   }  
			   rc+=1.0;
			}
			if (nd.b01.active) {
			   for (integer i=nd.b01.iseparate_start; i<=nd.b01.iseparate_finish; ++i) {
				   integer iP=ifrontregulation[i];
                   parcolor[iP]=rc;
			   }
			   rc+=1.0;
			}
			if (nd.b0.active) {
			   for (integer i=nd.b0.iseparate_start; i<=nd.b0.iseparate_finish; ++i) {
				   integer iP=ifrontregulation[i];
                   parcolor[iP]=rc;
			   }
			   rc+=1.0;
			}
        }
		/*
		// Шестнадцать ядер. (когда процессоры от Intel станут 8 ядерными. 20 апреля 2013.)
		if (nd.ncore==16) {
			doublereal rc=1.0;
			if (nd.b000.active) {
			   for (integer i=nd.b000.ileft_start; i<=nd.b000.ileft_finish; ++i) {
				   integer iP=ifrontregulation[i];
                   parcolor[iP]=rc;
			   }
			   rc+=1.0;
			   for (integer i=nd.b000.iright_start; i<=nd.b000.iright_finish; ++i) {
				   integer iP=ifrontregulation[i];
                   parcolor[iP]=rc;
			   }
			   rc+=1.0;
			}
			for (integer i=nd.b001.ileft_start; i<=nd.b001.ileft_finish; ++i) {
				integer iP=ifrontregulation[i];
                parcolor[iP]=3.0;
			}
			for (integer i=nd.b001.iright_start; i<=nd.b001.iright_finish; ++i) {
				integer iP=ifrontregulation[i];
                parcolor[iP]=4.0;
			}
			for (integer i=nd.b010.ileft_start; i<=nd.b010.ileft_finish; ++i) {
				integer iP=ifrontregulation[i];
                parcolor[iP]=5.0;
			}
			for (integer i=nd.b010.iright_start; i<=nd.b010.iright_finish; ++i) {
				integer iP=ifrontregulation[i];
                parcolor[iP]=6.0;
			}
			for (integer i=nd.b011.ileft_start; i<=nd.b011.ileft_finish; ++i) {
				integer iP=ifrontregulation[i];
                parcolor[iP]=7.0;
			}
			for (integer i=nd.b011.iright_start; i<=nd.b011.iright_finish; ++i) {
				integer iP=ifrontregulation[i];
                parcolor[iP]=8.0;
			}
			for (integer i=nd.b000.iseparate_start; i<=nd.b000.iseparate_finish; ++i) {
				integer iP=ifrontregulation[i];
                parcolor[iP]=9.0;
			}
			for (integer i=nd.b001.iseparate_start; i<=nd.b001.iseparate_finish; ++i) {
				integer iP=ifrontregulation[i];
                parcolor[iP]=10.0;
			}
			for (integer i=nd.b010.iseparate_start; i<=nd.b010.iseparate_finish; ++i) {
				integer iP=ifrontregulation[i];
                parcolor[iP]=11.0;
			}
			for (integer i=nd.b011.iseparate_start; i<=nd.b011.iseparate_finish; ++i) {
				integer iP=ifrontregulation[i];
                parcolor[iP]=12.0;
			}
			for (integer i=nd.b00.iseparate_start; i<=nd.b00.iseparate_finish; ++i) {
				integer iP=ifrontregulation[i];
                parcolor[iP]=13.0;
			}
			for (integer i=nd.b01.iseparate_start; i<=nd.b01.iseparate_finish; ++i) {
				integer iP=ifrontregulation[i];
                parcolor[iP]=14.0;
			}
			for (integer i=nd.b0.iseparate_start; i<=nd.b0.iseparate_finish; ++i) {
				integer iP=ifrontregulation[i];
                parcolor[iP]=15.0;
			}
        }
		*/

		// Визуализация критериальной функции.
		// Это очень важно для понимания работы алгоритма.
		 if (debug_export_tecplot) {
			// экспортируем в tecplot только если это действительно необъходимо.
		    exporttecplotxy360T_3D_part2nd(t.maxelm, t.ncell, fglobal, t, flow_interior, 0, parcolor);
		 }
		//system("pause");

		 if (parcolor != nullptr) {
			 delete[] parcolor;
		 }
}

// Быстрая сортировка Ч. Хоара. Время n*log2(n).
void quickSort(int* &ifrQ, int* &indx, integer left, integer right)
{
	integer i=left, j=right;
	int tmp;
	integer pivot=ifrQ[(left+right)/2];

	/* partition */
	while (i<=j) {
		while (ifrQ[i]<pivot)
			i++;
		while (ifrQ[j]>pivot)
			j--;
		if (i<=j) {
			tmp=ifrQ[i];
			ifrQ[i]=ifrQ[j];
			ifrQ[j]=tmp;
			tmp=indx[i];
			indx[i]=indx[j];
			indx[j]=tmp;
			i++;
			j--;
		}
	};

	/* recursion */
	if (left<j)
		quickSort(ifrQ,indx,left,j);
	if (i<right) 
        quickSort(ifrQ,indx,i,right);

	
}

// Двоичный поиск.
integer BinarySearchi(int* &A, integer key, integer n)
{
	integer left=0, right=n, mid;
	while (left<=right)
	{
		mid=left+(right-left)/2;
		if (key<A[mid]) right=mid-1;
		else if (key>A[mid]) left=mid+1;
		else return mid;
	}
	return -1;
}

// Алгоритм построения вложенных сечений.
void calc_front(FLOW* &fglobal, FLOW &f, TEMPER &t, integer flow_interior,
	int ls, int lw, WALL* w, PARDATA &nd, BLOCK* &b, int &lb, SOURCE* &s_loc) {

	bool bdiagnostic_pause=false;

	// nd - nested dissection.
	if ((nd.ncore==2)||(nd.ncore==4)||(nd.ncore==8)) {

	    bool bprintmessage=true; // печатать ли сообщения на консоль.

	    // Выделение памяти
	    doublereal** potent;
	    potent=new doublereal*[4];
	    for (integer i=0; i<4; ++i) potent[i]=new doublereal[f.maxelm+f.maxbound];
	    // матрица
	    equation3D_bon* slb=new equation3D_bon[f.maxbound];
	    equation3D* sl=new equation3D[f.maxelm];
	    doublereal* rthdsd = new doublereal[f.maxelm+f.maxbound];
	    // инициализация
	    for (integer i=0; i<4; ++i) for (integer j=0; j<f.maxelm+f.maxbound; ++j) potent[i][j]=0.0;
	    // сборка матрицы

		

	    integer iMCB=-1;
	    // Граничные условия Дирихле обязательно 
	    // должны собираться в первую очередь
	    doublereal dbeta=1.0; // первый порядок точности.
        for (integer i=0; i<f.maxbound; ++i) {
		    if (!bunion_8_08_2015) {
			my_nested_dissection_bound(i, f.maxelm, 
			              true, f.border_neighbor, ls, lw, w,
						  slb, dbeta, f.pa, f.nvtx,
						  potent[FIDISTW], iMCB);
			}
			else {
			my_nested_dissection_bound_union(i, f.maxelm, 
			              true, f.border_neighbor, ls, lw, w,
						  slb, dbeta, f.pa, f.nvtx,
						  potent[FIDISTW], iMCB);
			}
	    }

	    // Собираем однородные условия Неймана.
        // последний параметр bool bDirichlet равен false.
        for (integer i=0; i<f.maxbound; ++i) {
		    if (!bunion_8_08_2015) {
			     my_nested_dissection_bound(i, f.maxelm, 
			              false, f.border_neighbor, ls, lw, w,
						  slb, dbeta, f.pa, f.nvtx,
						  potent[FIDISTW], iMCB);
			}
			else {
			   my_nested_dissection_bound_union(i, f.maxelm, 
			              false, f.border_neighbor, ls, lw, w,
						  slb, dbeta, f.pa, f.nvtx,
						  potent[FIDISTW], iMCB);
			}
	    }

	    // Собираем матрицу для внутренних контрольных объёмов.
        for (int iP=0; iP<f.maxelm; ++iP) {
		    if (!bunion_8_08_2015 ) {
			    my_elmatr_quad_nested_dissection(iP, sl, slb,  
						f.pa, f.nvtx, f.neighbors_for_the_internal_node,
						f.maxelm, dbeta);
			}
			else {
			    my_elmatr_quad_nested_dissection_union(iP, sl, slb,  
						f.pa, f.nvtx, f.neighbors_for_the_internal_node,
						f.maxelm, dbeta);
			}
	    }

	    SIMPLESPARSE sparseM; // разреженная матрица
	    IMatrix sparseS; // разреженная матрица в формате IMatrix

	    // выделение памяти и инициализация для 
	    // простейшей разреженной матрицы.
	    initsimplesparse(sparseM, f.maxelm + f.maxbound );
	    initIMatrix(&sparseS, f.maxelm + f.maxbound);

	    // Для внутренних узлов расчётной сетки:
	    for (integer i=0; i<f.maxelm; ++i) {
		     addelmsimplesparse(sparseM, sl[i].ap, sl[i].iP, sl[i].iP, true);
		     setValueIMatrix(&sparseS, sl[i].iP, sl[i].iP, sl[i].ap);
		     rthdsd[sl[i].iP]=sl[i].b;

		     const doublereal nonzeroEPS=1e-37; // для отделения вещественного нуля

		          if ((sl[i].iE>-1) && (fabs(sl[i].ae) > nonzeroEPS)) {
				       addelmsimplesparse(sparseM, -sl[i].ae, sl[i].iP, sl[i].iE, true);
				       setValueIMatrix(&sparseS, sl[i].iP, sl[i].iE, -sl[i].ae);
			      }
			      if ((sl[i].iN>-1) && (fabs(sl[i].an) > nonzeroEPS)) {
				     addelmsimplesparse(sparseM, -sl[i].an, sl[i].iP, sl[i].iN, true);
				     setValueIMatrix(&sparseS, sl[i].iP, sl[i].iN, -sl[i].an);
			      }
			      if ((sl[i].iT>-1) && (fabs(sl[i].at) > nonzeroEPS)) {
				      addelmsimplesparse(sparseM, -sl[i].at, sl[i].iP, sl[i].iT, true);
                      setValueIMatrix(&sparseS, sl[i].iP, sl[i].iT, -sl[i].at);
			      }
			      if ((sl[i].iS>-1) && (fabs(sl[i].as) > nonzeroEPS)) {
				      addelmsimplesparse(sparseM, -sl[i].as, sl[i].iP, sl[i].iS, true);
                      setValueIMatrix(&sparseS, sl[i].iP, sl[i].iS, -sl[i].as);
			      }
			      if ((sl[i].iW>-1) && (fabs(sl[i].aw) > nonzeroEPS)) {
				     addelmsimplesparse(sparseM, -sl[i].aw, sl[i].iP, sl[i].iW, true);
                     setValueIMatrix(&sparseS, sl[i].iP, sl[i].iW, -sl[i].aw);
			      }
			      if ((sl[i].iB>-1) && (fabs(sl[i].ab) > nonzeroEPS)) {
				     addelmsimplesparse(sparseM, -sl[i].ab, sl[i].iP, sl[i].iB, true);
				     setValueIMatrix(&sparseS, sl[i].iP, sl[i].iB, -sl[i].ab);
			      }
	    }

	    // Для граничных узлов расчётной сетки:
	    for (integer i=0; i<f.maxbound; ++i) {
		     addelmsimplesparse(sparseM, slb[i].aw, slb[i].iW, slb[i].iW, true);
		     setValueIMatrix(&sparseS, slb[i].iW, slb[i].iW, slb[i].aw);
		     rthdsd[slb[i].iW]=slb[i].b;

		     const doublereal nonzeroEPS=1e-37; // для отделения вещественного нуля

		     if ((slb[i].iI>-1) && (fabs(slb[i].ai) > nonzeroEPS)) {
			     addelmsimplesparse(sparseM, -slb[i].ai, slb[i].iW, slb[i].iI, true);
                 setValueIMatrix(&sparseS, slb[i].iW, slb[i].iI, -slb[i].ai);
		     }
	    }


	    // решение СЛАУ
	    for (integer i=0; i<(f.maxelm+f.maxbound); ++i) potent[FIDISTW][i]=0.0;
        integer maxiter=10000; //120
	    //doublereal *val;
        //integer *col_ind, *row_ptr;

	    freeIMatrix(&sparseS); // освобождение оперативной памяти.
	    //simplesparsetoCRS(sparseM, val, col_ind, row_ptr, (f.maxelm+f.maxbound)); // преобразование матрицы из одного формата хранения в другой.

	    //delete val; delete col_ind; delete row_ptr;
	    //SOR3D(sl, slb, potent[FIDISTW], f.maxelm, f.maxbound, PAM); // проверка кода
	    // system("pause"); отладка
	    // MICCG решатель.
		if (f.maxelm + f.maxbound < 60000) {
	       ICCG(FIDISTW,sparseM, rthdsd, potent[FIDISTW], f.maxelm + f.maxbound ,bprintmessage,true,maxiter); //->//
		   // при использовании ICCG память из под sparseM освобождается автоматом
		}
		else {

			   // Для больших размерностей нам нужен очень робастый решатель, т.к. ICCG не справляется.
	           QuickMemVorst my_memory_bicgstab;
						  my_memory_bicgstab.ballocCRScfd=false; // выделяем память.
						  my_memory_bicgstab.bsignalfreeCRScfd=false; // не уничтожаем память. (еще рано).
						  // Инициализация указателей !
						  my_memory_bicgstab.val=nullptr;
			              my_memory_bicgstab.col_ind=nullptr;
			              my_memory_bicgstab.row_ptr=nullptr;
			              my_memory_bicgstab.ri=nullptr;
			              my_memory_bicgstab.roc=nullptr;
			              my_memory_bicgstab.s=nullptr;
			              my_memory_bicgstab.t=nullptr;
			              my_memory_bicgstab.vi=nullptr;
			              my_memory_bicgstab.pi=nullptr;
			              my_memory_bicgstab.dx=nullptr;
			              my_memory_bicgstab.dax=nullptr;
			              my_memory_bicgstab.y=nullptr;
			              my_memory_bicgstab.z=nullptr;
			              my_memory_bicgstab.a=nullptr;
			              my_memory_bicgstab.ja=nullptr;
			              my_memory_bicgstab.ia=nullptr;
			              my_memory_bicgstab.alu=nullptr;
			              my_memory_bicgstab.jlu=nullptr;
			              my_memory_bicgstab.ju=nullptr;
						   my_memory_bicgstab.alu1=nullptr;
			              my_memory_bicgstab.jlu1=nullptr;
			              my_memory_bicgstab.ju1=nullptr;
						  my_memory_bicgstab.x1=nullptr;
			              my_memory_bicgstab.iw=nullptr;
			              my_memory_bicgstab.levs=nullptr;
			              my_memory_bicgstab.w=nullptr;
			              my_memory_bicgstab.jw=nullptr;
						  // Иногда совместно с уравнениями гидродинамики решается и уравнение теплопередачи.
			              my_memory_bicgstab.ballocCRSt=false; // Выделять память
			              my_memory_bicgstab.bsignalfreeCRSt=false; // и сразу не освобождать.
			              // инициализация указателей.
                          my_memory_bicgstab.tval=nullptr;
			              my_memory_bicgstab.tcol_ind=nullptr;
			              my_memory_bicgstab.trow_ptr=nullptr;
			              my_memory_bicgstab.tri=nullptr;
			              my_memory_bicgstab.troc=nullptr;
			              my_memory_bicgstab.ts=nullptr;
			              my_memory_bicgstab.tt=nullptr;
			              my_memory_bicgstab.tvi=nullptr;
			              my_memory_bicgstab.tpi=nullptr;
			              my_memory_bicgstab.tdx=nullptr;
			              my_memory_bicgstab.tdax=nullptr;
			              my_memory_bicgstab.ty=nullptr;
			              my_memory_bicgstab.tz=nullptr;
			              my_memory_bicgstab.ta=nullptr;
			              my_memory_bicgstab.tja=nullptr;
			              my_memory_bicgstab.tia=nullptr;
			              my_memory_bicgstab.talu=nullptr;
			              my_memory_bicgstab.tjlu=nullptr;
			              my_memory_bicgstab.tju=nullptr;
			              my_memory_bicgstab.tiw=nullptr;
			              my_memory_bicgstab.tlevs=nullptr;
			              my_memory_bicgstab.tw=nullptr;
			              my_memory_bicgstab.tjw=nullptr;
						  my_memory_bicgstab.icount_vel=100000; // очень большое число.

						  bool worked_successfully = false;
						  const integer iHAVorstModification_id = 0;// 0 не используется.
						  const bool bprint_message = true;
		   amg(sl, slb, f.maxelm, f.maxbound, rthdsd, potent[FIDISTW], 1.0, FIDISTW,false,my_memory_bicgstab, 
			   f.ifrontregulationgl, f.ibackregulationgl, iHAVorstModification_id, worked_successfully,b,lb,
			   s_loc,ls, w, lw, t.whot_is_block, bprint_message);
		}

	    // Вычисление градиентов.
	    // на основе поля скорости удовлетворяющего уравнению неразрывности.
	    for (integer i=0; i<f.maxelm; ++i) {
		    // градиенты потенциала FI для внутренних КО.
			if (!bunion_8_08_2015) {
	            green_gauss_FI_(i, potent, f.nvtx, f.pa,
	                        f.neighbors_for_the_internal_node, f.maxelm, false);
			}
			else {
			    green_gauss_FI_union(i, potent, f.nvtx, f.pa,
	                        f.neighbors_for_the_internal_node, f.maxelm, false);
			}
	    }
	    for (integer i=0; i<f.maxelm; ++i) {
		     // градиенты потенциала FI для граничных КО.
			if (!bunion_8_08_2015) {
	             green_gauss_FI_(i, potent, f.nvtx, f.pa,
	                         f.neighbors_for_the_internal_node, f.maxelm, true);
			}
			else {
			    green_gauss_FI_union(i, potent, f.nvtx, f.pa,
	                         f.neighbors_for_the_internal_node, f.maxelm, true);
			}
        }

	    // Вычисление расстояния до ближайшей стенки.
	    doublereal* front=new doublereal[f.maxelm+f.maxbound];
		doublereal* front_orig=new doublereal[f.maxelm+f.maxbound];

	    for (integer iP=0; iP<f.maxelm+f.maxbound; iP++) {
		      // формула (223) из приложения 2 в книге Снегирёва.
		      // А также смотри: Transport equation based wall distance calculation
	 	      // www.cfd-online.com/Wiki/
		      // В обоих источниках этой формулы содержится неточность.
		      // Так как мы берём отрицательный источниковый член -1.0*dx*dy*dz
		      // То это эквивалентно "холодильнику" и поэтому внутри расчётной области
		      // значения потенциала FI будут отрицателдьными. Но подкоренное выражение
		      // не может быть отрицательным,  поэтому нужно перед 2.0 поставить минус,
	          // тем самым -2.0*FI будет положительной величиной и с подкоренным выражением 
		      // будет всё впорядке.
		      front[iP]=-sqrt(potent[GRADXFI][iP]*potent[GRADXFI][iP]+
			                  potent[GRADYFI][iP]*potent[GRADYFI][iP]+
							  potent[GRADZFI][iP]*potent[GRADZFI][iP])
					    +sqrt(potent[GRADXFI][iP]*potent[GRADXFI][iP]+
							  potent[GRADYFI][iP]*potent[GRADYFI][iP]+
							  potent[GRADZFI][iP]*potent[GRADZFI][iP]-
							  2.0*potent[FIDISTW][iP]);
			  front_orig[iP]=front[iP]; // оригинальный фронт связанный с изначальнной нумерацией.
	    }

	    doublereal dmax=-1.0;
	    for (integer i=0; i<f.maxelm+f.maxbound; ++i) {
	    	dmax=fmax(dmax,front[i]);
	    }

		// Визуализация критериальной функции.
		// Это очень важно для понимания работы алгоритма.
		printf("distance found\n");
		if (debug_export_tecplot) {
		   // экспортируем в tecplot только если это действительно необъходимо.
		   exporttecplotxy360T_3D_part2nd(t.maxelm, t.ncell, fglobal, t, flow_interior, 0, front);
		}
		if (bdiagnostic_pause) {
			//system("pause");
			system("pause"); 
		}
		 /*debug
		for (integer i=0; i<f.maxelm+f.maxbound; ++i) {
			if ((i>=8820)&&(i<=9959)) {
				front[i]=1.25*dmax; // dmax + 25% * dmax;
			}
		}
		exporttecplotxy360T_3D_part2nd(t.maxelm, t.ncell, fglobal, t, flow_interior, 0, front);
		system("pause");
		//*/

	    // максимальное расстояние.
	    // В случае нескольких входных границ с заданной на них нетривиальной скоростью может потребоваться
	    // выбрать границу для которой величина dmax максимальна. Массив front даёт гораздо более лучшее упорядочивание в случае наибольшего dmax. // TODO.
	
	    // Освобождение памяти.
	    for (integer i=0; i<4; ++i) delete[] potent[i];
	    delete[] potent;
	    // освобождение памяти из под матрицы.
	    delete[] sl; delete[] slb;
	    delete[] rthdsd;

	    // Построение вложенных сечений на основе массива front.
	    printf("construct nested dissection start...\n");
	    int *ifrontregulation=new int[f.maxelm+f.maxbound];
	    for (int i=0; i<f.maxelm+f.maxbound; ++i) {
		      ifrontregulation[i]=i; // естественная нумерация в программе.
	    }
	    // упорядочим элементы в соответствии с возрастанием элементов в массиве front.
		printf("Quick sort start...\n");
	    QuickSortnd(ifrontregulation,  front, 0, f.maxelm+f.maxbound-1);
		/*
		bool swapped=true;
		while (swapped) {
			swapped=false;
			for (integer i=1; i<f.maxelm+f.maxbound; ++i) {
				if (front[i-1]>front[i]) {
					doublereal r=front[i-1];
					front[i-1]=front[i];
					front[i]=r;
					integer ir=ifrontregulation[i-1];
					ifrontregulation[i-1]=ifrontregulation[i];
					ifrontregulation[i]=ir;
					swapped=true;
				}
			}
		}*/
		printf("Quick sort finish...\n");
	    // Теперь все элементы переупорядочены в соответствии со значениями front.
	    // и front тоже переупорядочен в порядке возрастания.
	    // dmax есть максимальная длина !

		printf("dmax=%e, \n",dmax);

		ilevel_info_nested_desection=1;

		// Параметры t, fglobal, flow_interior - нужны лишь для графической визуализации которая нужна для отладки.
		separate2(nd.b0,front,front_orig,f,ifrontregulation,dmax,0,f.maxelm+f.maxbound-1,t,fglobal,flow_interior); // два ядра.
        // Визуализация разбиения !!!
		visible_partition(fglobal, f, flow_interior, t, ifrontregulation, nd);
		//QuickSortnd(ifrontregulation,  front, nd.b0.ileft_start,nd.b0.ileft_finish);
		//QuickSortnd(ifrontregulation,  front, nd.b0.iright_start,nd.b0.iright_finish);
		//printf("sort complete...\n");
		// Визуализация критериальной функции.
		// Это очень важно для понимания работы алгоритма.
		//exporttecplotxy360T_3D_part2nd(t.maxelm, t.ncell, fglobal, t, flow_interior, 0, front);
		//system("pause");

		if ((nd.ncore==4)||(nd.ncore==8)) {

			ilevel_info_nested_desection=2;

			// четыре ядра или 8 ядер.
			if (nd.b0.active) {
			    separate2(nd.b00,front,front_orig,f,ifrontregulation,dmax,nd.b0.ileft_start,nd.b0.ileft_finish,t,fglobal,flow_interior);
				// Визуализация разбиения !!!
		        visible_partition(fglobal, f, flow_interior, t, ifrontregulation, nd);
			    separate2(nd.b01,front,front_orig,f,ifrontregulation,dmax,nd.b0.iright_start,nd.b0.iright_finish,t,fglobal,flow_interior);
				// Визуализация разбиения !!!
		        visible_partition(fglobal, f, flow_interior, t, ifrontregulation, nd);
			}
		}
		if (nd.ncore==8) {

			ilevel_info_nested_desection=3;

			// Восемь ядер.
			if (nd.b00.active) {
               separate2(nd.b000,front,front_orig,f,ifrontregulation,dmax,nd.b00.ileft_start,nd.b00.ileft_finish,t,fglobal,flow_interior);
			   // Визуализация разбиения !!!
		        visible_partition(fglobal, f, flow_interior, t, ifrontregulation, nd);
			   separate2(nd.b001,front,front_orig,f,ifrontregulation,dmax,nd.b00.iright_start,nd.b00.iright_finish,t,fglobal,flow_interior);
			   // Визуализация разбиения !!!
		        visible_partition(fglobal, f, flow_interior, t, ifrontregulation, nd);
			}
			if (nd.b01.active) {
			   separate2(nd.b010,front,front_orig,f,ifrontregulation,dmax,nd.b01.ileft_start,nd.b01.ileft_finish,t,fglobal,flow_interior);
			   // Визуализация разбиения !!!
		        visible_partition(fglobal, f, flow_interior, t, ifrontregulation, nd);
			   separate2(nd.b011,front,front_orig,f,ifrontregulation,dmax,nd.b01.iright_start,nd.b01.iright_finish,t,fglobal,flow_interior);
			   // Визуализация разбиения !!!
		        visible_partition(fglobal, f, flow_interior, t, ifrontregulation, nd);
			}
		}

		// Визуализация разбиения !!!
		visible_partition(fglobal, f, flow_interior, t, ifrontregulation, nd);
		

	    printf("construct nested dissection sucsefull finesed. Ok\n");

		f.ifrontregulationgl=new int[f.maxelm+f.maxbound];
	    for (integer i=0; i<f.maxelm+f.maxbound; ++i) {
			// получаем iP по номеру i.
		      f.ifrontregulationgl[i]=ifrontregulation[i]; // Это и есть искомое переупорядочивание nested desection.
	    }

		f.ibackregulationgl=new int[f.maxelm+f.maxbound];

		/*
		// массив не отсортирован.
		
#pragma omp for
		for (integer iP=0; iP<f.maxelm+f.maxbound; iP++) {
			for (integer i=0; i<f.maxelm+f.maxbound; ++i) {
				if (iP==ifrontregulation[i]) {
					// получаем  i по номеру iP.
					ibackregulationgl[iP]=i;
					break;
				}
			}
		}
		*/

		// Есть идея заменить квадратичный поиск Быстрой сортировкой и бинарным поиском.
		int *ifrQ=new int[f.maxelm+f.maxbound];
		int *indx=new int[f.maxelm+f.maxbound];
		// copy
		for (int iP=0; iP<f.maxelm+f.maxbound; ++iP) {
			ifrQ[iP]=ifrontregulation[iP];
			indx[iP]=iP;
		}
		// Быстрая сортировка:
		quickSort(ifrQ,indx,0,f.maxelm+f.maxbound-1);
		// Двоичный поиск:
		for (integer iP = 0; iP < f.maxelm + f.maxbound; iP++) {
			// iP надо найти в ifrontregulation[i]
			integer i = BinarySearchi(ifrQ, iP, f.maxelm + f.maxbound - 1);
			if ((i >= 0) && (i < f.maxelm + f.maxbound)) {
			f.ibackregulationgl[iP] = indx[i];
		}
			else {
				printf("calc front error ibackregulationgl\n");
				//system("pause");
				system("pause");
				exit(1);
			}
		}

		

		delete[] ifrQ;
		delete[] indx;

	    // Освобождение памяти из под фронта.
	    delete[] front;
		delete[] front_orig;
	    delete[] ifrontregulation;

	}
}


// 18.01.2019
// метод получает будущий вектор цвета не инициализированный, но память выделена.
// 1 - первая область, 2 вторая область, 3 - граница областей. Для number_part=2.
// 1 певая часть, 2 часть, 3 - сшивка, 4 третья часть, 5 - четвертая часть, 6 - сшивка, 7 сшивка. Для number_part=4.
// связи хранятся в матрице в координатном формате crs.
// length_separator = 1, 2, 3, ... клетки.
// number_part - число частей, степень двойки 1, 2, 4, 8, 16, 32,...
// Нумерация начинается с нуля.
void nested_desection_crs(integer*& col_ind, integer*& row_ptr, integer n, integer*& color, integer& dist_max
/*integer number_part, integer length_separator*/)
{
	typedef struct TSTACK_NODE_ND {
		integer dist, inumber_node;
	} STACK_NODE_ND;
	integer* dist = new integer[n];
	STACK_NODE_ND* stack = new STACK_NODE_ND[n];// стек.
	/*if (number_part == 1) {
		for (integer i_1 = 0; i_1 < n; ++i_1) color[i_1] = 1;
	}
	else
	{*/
	integer itop = 0; // вершина стека.
	integer ibottom = 0; // Начало очереди.
	for (integer i_1 = 0; i_1 < n; ++i_1) dist[i_1] = 0;
	integer ifirst = 0;
	dist[ifirst] = 1;
	for (integer i_1 = row_ptr[ifirst]; i_1 <= row_ptr[ifirst + 1] - 1; ++i_1) {
		integer ja = col_ind[i_1];
		if (ja != ifirst) {
			// внедиагональный элемент.
			dist[ja] = 2;
			stack[itop].inumber_node = ja;
			stack[itop].dist = 2;
			itop++;
		}
	}
	while ((ibottom<itop)&&(ibottom < n)) {
		integer icandidate = stack[ibottom].inumber_node;
		integer dist_current = stack[ibottom].dist + 1;
		ibottom++;
		for (integer i_1 = row_ptr[icandidate]; i_1 <= row_ptr[icandidate + 1] - 1; ++i_1) {
			integer ja = col_ind[i_1];
			if (ja != icandidate) {
				// внедиагональный элемент.
				if (dist[ja] == 0) {
					// Если узел еще не был посещен.
					dist[ja] = dist_current;
					if (itop < n) {
						stack[itop].inumber_node = ja;
						stack[itop].dist = dist_current;
						itop++;
					}
					else {
						printf("CORE DUMP: STACK OVERFLOW!!!!\n");
						system("pause");
					}
				}
			}
		}
	}
	// поиск максимального расстояния.
	integer max = -1;
	for (integer i_1 = 0; i_1 < n; ++i_1) {
		if ((row_ptr[i_1 + 1] - 1 - row_ptr[i_1]) > 1) {
		    if (dist[i_1] > max) {
				max = dist[i_1];
				ifirst = i_1;
			}
		}
	}
	itop = 0; // вершина стека.
	ibottom = 0; // Начало очереди.
	for (integer i_1 = 0; i_1 < n; ++i_1) dist[i_1] = 0;
	dist[ifirst] = 1;
	if (row_ptr[ifirst + 1] - 1 - row_ptr[ifirst] == 0) {
		printf("error limit\n");
		system("pause");
	}
	for (integer i_1 = row_ptr[ifirst]; i_1 <= row_ptr[ifirst + 1] - 1; ++i_1) {
		integer ja = col_ind[i_1];
		if (ja != ifirst) {
			// внедиагональный элемент.
			dist[ja] = 2;
			stack[itop].inumber_node = ja;
			stack[itop].dist = 2;
			itop++;
		}
	}
	while ((ibottom < itop) && (ibottom < n)) {
		integer icandidate = stack[ibottom].inumber_node;
		integer dist_current = stack[ibottom].dist + 1;
		ibottom++;
		for (integer i_1 = row_ptr[icandidate]; i_1 <= row_ptr[icandidate + 1] - 1; ++i_1) {
			integer ja = col_ind[i_1];
			if (ja != icandidate) {
				// внедиагональный элемент.
				if (dist[ja] == 0) {
					// Если узел еще не был посещен.
					dist[ja] = dist_current;
					if (itop < n) {
						stack[itop].inumber_node = ja;
						stack[itop].dist = dist_current;
						itop++;
					}
					else {
						printf("CORE DUMP: STACK OVERFLOW!!!!\n");
						system("pause");
					}
				}
			}
		}
	}

	max = -1;
	for (integer i_1 = 0; i_1 < n; ++i_1) {
		if (dist[i_1] > max) {
			max = dist[i_1];
			ifirst = i_1;
		}
	}

	dist_max = max;
	for (integer i_1 = 0; i_1 < n; ++i_1) {
		color[i_1] = dist[i_1];
	}

	/*
	if (number_part == 2) {
		// две части.
		integer isep = 0.5 * (max);
		for (integer i_1 = 0; i_1 < n; ++i_1) {
			if (length_separator == 1) {
				if (dist[i_1] < isep) {
					color[i_1] = 1; // первая часть.
				}
				if (dist[i_1] > isep) {
					color[i_1] = 2; // вторая часть.
				}
				if (dist[i_1] == isep) {
					color[i_1] = 3; // разделитель.
				}
			}
			if (length_separator == 2) {
				if (dist[i_1] < isep-1) {
					color[i_1] = 1; // первая часть.
				}
				if (dist[i_1] > isep) {
					color[i_1] = 2; // вторая часть.
				}
				if ((dist[i_1] == isep)||(dist[i_1] == isep-1)) {
					color[i_1] = 3; // разделитель.
				}
			}
			if (length_separator == 3) {
				if (dist[i_1] < isep - 1) {
					color[i_1] = 1; // первая часть.
				}
				if (dist[i_1] > isep+1) {
					color[i_1] = 2; // вторая часть.
				}
				if ((dist[i_1] == isep) || (dist[i_1] == isep - 1) || (dist[i_1] == isep + 1)) {
					color[i_1] = 3; // разделитель.
				}
			}
		}
		*/




//}

delete[] dist;
delete[] stack;
} // nested_desection_crs


#endif