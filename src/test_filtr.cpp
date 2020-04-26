// Файл test_filtr.cpp содержит исходный код
// для получения дважды фильтрованного поля скоростей, а 
// также производных величин от него: 
// 1. градиента дважды фильтрованного поля скорости,
// 2. вихря (модуля ротора скорости) для дважды фильтрованного поля скорости,
// 3. угла между вихрем и фильтрованным вихрем (используется для Selective Smagorinsky SS),
// 4. инварианта тензора скоростей деформаций для дважды фильтрованного поля скоростей.
// begin 19 апреля 2012 года.

#ifndef MY_TEST_FILTR_CPP
#define MY_TEST_FILTR_CPP 1

// Пользователя предлагается на выбор несколько фильтров.
// 1) среднее по объёму, 2) Trapezoidal Filtr, 3) Simpson filtr.

// Самый простой фильтр просто среднее по объёму значение.
// layer TOP
// 1	1	1
// 1	1	1
// 1	1	1
// layer center
// 1	1	1
// 1	1	1
// 1	1	1
// layer Bottom
// 1	1	1
// 1	1	1
// 1	1	1
// множитель 1.0/27.0
#define VOLUME_AVERAGE_FILTR 0
// Фильтры Trapezoidal и SIMPSON пригодны по видимому (в теории) только для равномерной сетки.
// Trapezoidal Filtr
// layer TOP
// 1	2	1
// 2	4	2
// 1	2	1
// layer center
// 2	4	2
// 4	8	4
// 2	4	2
// layer Bottom
// 1	2	1
// 2	4	2
// 1	2	1
// множитель 1.0/64.0.
#define TRAPEZOIDAL_FILTR 1
// Simpson filtr
// layer TOP
// 1	4	1
// 4	16	4
// 1	4	1
// layer center
// 4	16	4
// 16	64	16
// 4	16	4
// layer Bottom
// 1	4	1
// 4	16	4
// 1	4	1
// множитель 1.0/216.0.
#define SIMPSON_FILTR 2

// представляем метод flattener - сглаживатель который
// проектируется с целью убрать не физичный скачок производной или 
// полевой величины зависящей от производных при переходе через плоскость
// смены шага сетки.
// Неоспоримо что данный метод должен подлежать наибольшей критике из-за его необоснованности.
// Внимание ! на полностью равномерной сетке данный метод не изменяет передаваемый вектор potent.
// Априори предполагается что пользователь использует почти равномерную расчётную сетку.
// В общем не обязательно применять этот код. Данный код скорее служит для исследования того как он повлияет на алгоритм.
// Предположительно этот метод требуется применить сразу после вычисления производных (или градиента величины).
void flattener(doublereal* &potent, integer maxelm, integer maxbound, ALICE_PARTITION** neighbors_for_the_internal_node,
	           integer** nvtx, TOCHKA* pa, BOUND* border_neighbor) {

	doublereal* potentcopy = nullptr;
	potentcopy=new doublereal[maxelm + maxbound]; // рабочий вектор сохраняющий передаваемое поле.
	for (integer i=0; i<maxelm+maxbound; i++) {
		potentcopy[i]=potent[i];
	}

	bool** avgonsosed=new bool*[maxelm];
	for (integer i=0; i<maxelm; i++) avgonsosed[i]=new bool[6];
	for (integer i=0; i<maxelm; i++) for (integer j=0; j<6; j++) avgonsosed[i][j]=false;

	for (integer i=0; i<maxelm+maxbound; i++) {
		potent[i]=0.0; // обнуление.
	}

	for (integer iP=0; iP<maxelm; iP++) {

		  doublereal dx, dy, dz;

		  
		  dx = 0.0; dy = 0.0; dz = 0.0; // обнуляющий сброс значений
		
		  volume3D(iP, nvtx, pa, dx, dy, dz);

		// цикл по всем внутренним контрольным объёмам.
		  bool bE, bW, bN, bS, bT, bB;
		  bE=false; bW=false; bN=false; bS=false; bT=false; bB=false;

		  integer iE, iN, iT, iW, iS, iB; // номера соседних контрольных объёмов
	      iE=neighbors_for_the_internal_node[ESIDE][iP].iNODE1; iN=neighbors_for_the_internal_node[NSIDE][iP].iNODE1; iT=neighbors_for_the_internal_node[TSIDE][iP].iNODE1; iW=neighbors_for_the_internal_node[WSIDE][iP].iNODE1; iS=neighbors_for_the_internal_node[SSIDE][iP].iNODE1; iB=neighbors_for_the_internal_node[BSIDE][iP].iNODE1;

          if (iE>=maxelm) bE=true; // если true то узел является граничным.
	      if (iN>=maxelm) bN=true;
	      if (iT>=maxelm) bT=true;
          if (iW>=maxelm) bW=true;
	      if (iS>=maxelm) bS=true;
	      if (iB>=maxelm) bB=true;

		  if (!bE) {
			  doublereal dx1, dy1, dz1;

		      dx1=0.0; dy1=0.0; dz1=0.0; // обнуляющий сброс значений
		      volume3D(iE, nvtx, pa, dx1, dy1, dz1);

			  if (fabs(dx-dx1)>1e-20) avgonsosed[iP][ESIDE]=true;
		  }

		  if (!bW) {
			  doublereal dx1, dy1, dz1;

		      dx1=0.0; dy1=0.0; dz1=0.0; // обнуляющий сброс значений
		      volume3D(iW, nvtx, pa, dx1, dy1, dz1);

			  if (fabs(dx-dx1)>1e-20) avgonsosed[iP][WSIDE]=true;
		  }

		  if (!bN) {
			  doublereal dx1, dy1, dz1;

		      dx1=0.0; dy1=0.0; dz1=0.0; // обнуляющий сброс значений
		      volume3D(iN, nvtx, pa, dx1, dy1, dz1);

			  if (fabs(dy-dy1)>1e-20) avgonsosed[iP][NSIDE]=true;
		  }

		  if (!bS) {
			  doublereal dx1, dy1, dz1;

		      dx1=0.0; dy1=0.0; dz1=0.0; // обнуляющий сброс значений
		      volume3D(iS, nvtx, pa, dx1, dy1, dz1);

			  if (fabs(dy-dy1)>1e-20) avgonsosed[iP][SSIDE]=true;
		  }

		  if (!bT) {
			  doublereal dx1, dy1, dz1;

		      dx1=0.0; dy1=0.0; dz1=0.0; // обнуляющий сброс значений
		      volume3D(iT, nvtx, pa, dx1, dy1, dz1);

			  if (fabs(dz-dz1)>1e-20) avgonsosed[iP][TSIDE]=true;
		  }

		  if (!bB) {
			  doublereal dx1, dy1, dz1;

		      dx1=0.0; dy1=0.0; dz1=0.0; // обнуляющий сброс значений
		      volume3D(iB, nvtx, pa, dx1, dy1, dz1);

			  if (fabs(dz-dz1)>1e-20) avgonsosed[iP][BSIDE]=true;
		  }
	}

	for (integer iP=0; iP<maxelm; iP++) {

		// цикл по всем внутренним контрольным объёмам.
		  bool bE, bW, bN, bS, bT, bB;
		  bE=false; bW=false; bN=false; bS=false; bT=false; bB=false;

		  integer iE, iN, iT, iW, iS, iB; // номера соседних контрольных объёмов
	      iE=neighbors_for_the_internal_node[ESIDE][iP].iNODE1; iN=neighbors_for_the_internal_node[NSIDE][iP].iNODE1; iT=neighbors_for_the_internal_node[TSIDE][iP].iNODE1; iW=neighbors_for_the_internal_node[WSIDE][iP].iNODE1; iS=neighbors_for_the_internal_node[SSIDE][iP].iNODE1; iB=neighbors_for_the_internal_node[BSIDE][iP].iNODE1;

          if (iE>=maxelm) bE=true; // если true то узел является граничным.
	      if (iN>=maxelm) bN=true;
	      if (iT>=maxelm) bT=true;
          if (iW>=maxelm) bW=true;
	      if (iS>=maxelm) bS=true;
	      if (iB>=maxelm) bB=true;

		  doublereal dx, dy, dz;

		  dx=0.0; dy=0.0; dz=0.0; // обнуляющий сброс значений
		  volume3D(iP, nvtx, pa, dx, dy, dz);

		if (!(avgonsosed[iP][ESIDE]||avgonsosed[iP][WSIDE]||avgonsosed[iP][NSIDE]||avgonsosed[iP][SSIDE]||avgonsosed[iP][TSIDE]||avgonsosed[iP][BSIDE])) {
			potent[iP]=potentcopy[iP]; // оставляем значение без изменений (это справедливо для всюду равномерной сетки).
		}
		else {
			// Среднее арифметическое - состоятельная оценка математического ожидания.

			doublereal rn=0.0;
			if (avgonsosed[iP][ESIDE]) {

				doublereal dx1, dy1, dz1;

		        dx1=0.0; dy1=0.0; dz1=0.0; // обнуляющий сброс значений
		        volume3D(iE, nvtx, pa, dx1, dy1, dz1);

				potent[iP]+=0.5*(potentcopy[iP]/(dx*dy*dz)+potentcopy[iE]/(dx1*dy1*dz1));

				rn+=1.0;
			}

            if (avgonsosed[iP][WSIDE]) {

				doublereal dx1, dy1, dz1;

		        dx1=0.0; dy1=0.0; dz1=0.0; // обнуляющий сброс значений
		        volume3D(iW, nvtx, pa, dx1, dy1, dz1);

				potent[iP]+=0.5*(potentcopy[iP]/(dx*dy*dz)+potentcopy[iW]/(dx1*dy1*dz1));

				rn+=1.0;
			}

			if (avgonsosed[iP][NSIDE]) {

				doublereal dx1, dy1, dz1;

		        dx1=0.0; dy1=0.0; dz1=0.0; // обнуляющий сброс значений
		        volume3D(iN, nvtx, pa, dx1, dy1, dz1);

				potent[iP]+=0.5*(potentcopy[iP]/(dx*dy*dz)+potentcopy[iN]/(dx1*dy1*dz1));

				rn+=1.0;
			}

            if (avgonsosed[iP][SSIDE]) {

				doublereal dx1, dy1, dz1;

		        dx1=0.0; dy1=0.0; dz1=0.0; // обнуляющий сброс значений
		        volume3D(iS, nvtx, pa, dx1, dy1, dz1);

				potent[iP]+=0.5*(potentcopy[iP]/(dx*dy*dz)+potentcopy[iS]/(dx1*dy1*dz1));

				rn+=1.0;
			}

			if (avgonsosed[iP][TSIDE]) {

				doublereal dx1, dy1, dz1;

		        dx1=0.0; dy1=0.0; dz1=0.0; // обнуляющий сброс значений
		        volume3D(iT, nvtx, pa, dx1, dy1, dz1);

				potent[iP]+=0.5*(potentcopy[iP]/(dx*dy*dz)+potentcopy[iT]/(dx1*dy1*dz1));

				rn+=1.0;
			}

            if (avgonsosed[iP][BSIDE]) {

				doublereal dx1, dy1, dz1;

		        dx1=0.0; dy1=0.0; dz1=0.0; // обнуляющий сброс значений
		        volume3D(iB, nvtx, pa, dx1, dy1, dz1);

				potent[iP]+=0.5*(potentcopy[iP]/(dx*dy*dz)+potentcopy[iB]/(dx1*dy1*dz1));

				rn+=1.0;
			}

			potent[iP]/=rn; // истинное среднее арифметическое.
		}
	}

	// Для получения значений на границе нужно применить метод квадратичной интерполяции.
	// нужно применить квадратичную интерполяцию для 
	// нахождения фильтрованной величины в граничных контрольных объёмах.
	for (integer iG=0; iG<maxbound; iG++) {
		// цикл по всем граничным контрольным объёмам.

		doublereal dx=0.0, dy=0.0, dz=0.0;// объём текущего контрольного объёма
		// вычисление размеров текущего контрольного объёма:
	    volume3D(border_neighbor[iG].iI, nvtx, pa, dx, dy, dz);
		TOCHKA pp,pb,pbb;

		// сканируем внутреннюю нормаль
		switch (border_neighbor[iG].Norm) {
		  case ESIDE: // внешняя нормаль W
			       // квадратичная интерполяция.
			      
				   center_cord3D(border_neighbor[iG].iI, nvtx, pa, pp,WSIDE);
				   center_cord3D(border_neighbor[iG].iII, nvtx, pa, pb,WSIDE);
			       center_cord3D(neighbors_for_the_internal_node[ESIDE][border_neighbor[iG].iII].iNODE1, nvtx, pa, pbb,ESIDE);

				   potent[border_neighbor[iG].iB]=my_quadratic_interpolation('-', potent[neighbors_for_the_internal_node[ESIDE][border_neighbor[iG].iII].iNODE1], potent[border_neighbor[iG].iII], potent[border_neighbor[iG].iI], pbb.x , pb.x, pp.x, pp.x-0.5*dx);
				   
			       break;
		  case WSIDE: // внешняя нормаль E
			       // квадратичная интерполяция.

				   
		           center_cord3D(border_neighbor[iG].iI, nvtx, pa, pp,ESIDE);
		           center_cord3D(border_neighbor[iG].iII, nvtx, pa, pb,EE);
			       center_cord3D(neighbors_for_the_internal_node[WSIDE][border_neighbor[iG].iII].iNODE1, nvtx, pa, pbb,WSIDE);
					
			       potent[border_neighbor[iG].iB]=my_quadratic_interpolation('+', potent[neighbors_for_the_internal_node[WSIDE][border_neighbor[iG].iII].iNODE1], potent[border_neighbor[iG].iII], potent[border_neighbor[iG].iI], pbb.x , pb.x, pp.x, pp.x+0.5*dx);

			       break;
		  case NSIDE: // внешняя нормаль S
			       // квадратичная интерполяция.
			       
			       
		           center_cord3D(border_neighbor[iG].iI, nvtx, pa, pp,SSIDE);
		           center_cord3D(border_neighbor[iG].iII, nvtx, pa, pb,SS);
			       center_cord3D(neighbors_for_the_internal_node[NSIDE][border_neighbor[iG].iII].iNODE1, nvtx, pa, pbb,NSIDE);

			       potent[border_neighbor[iG].iB]=my_quadratic_interpolation('-', potent[neighbors_for_the_internal_node[NSIDE][border_neighbor[iG].iII].iNODE1], potent[border_neighbor[iG].iII], potent[border_neighbor[iG].iI], pbb.y , pb.y, pp.y, pp.y-0.5*dy);
			  
			       break;
		  case SSIDE:// внешняя нормаль N
			       // квадратичная интерполяция.

			       
		           center_cord3D(border_neighbor[iG].iI, nvtx, pa, pp,NSIDE);
		           center_cord3D(border_neighbor[iG].iII, nvtx, pa, pb,NN);
			       center_cord3D(neighbors_for_the_internal_node[SSIDE][border_neighbor[iG].iII].iNODE1, nvtx, pa, pbb,SSIDE);

			       potent[border_neighbor[iG].iB]=my_quadratic_interpolation('+', potent[neighbors_for_the_internal_node[SSIDE][border_neighbor[iG].iII].iNODE1], potent[border_neighbor[iG].iII], potent[border_neighbor[iG].iI], pbb.y , pb.y, pp.y, pp.y+0.5*dy);
			  
			       break;
		  case TSIDE: // внешняя нормаль B
			       // квадратичная интерполяция.

		           center_cord3D(border_neighbor[iG].iI, nvtx, pa, pp,BSIDE);
		           center_cord3D(border_neighbor[iG].iII, nvtx, pa, pb,BB);
			       center_cord3D(neighbors_for_the_internal_node[TSIDE][border_neighbor[iG].iII].iNODE1, nvtx, pa, pbb,TSIDE);


			       potent[border_neighbor[iG].iB]=my_quadratic_interpolation('-', potent[neighbors_for_the_internal_node[TSIDE][border_neighbor[iG].iII].iNODE1], potent[border_neighbor[iG].iII], potent[border_neighbor[iG].iI], pbb.z , pb.z, pp.z, pp.z-0.5*dz);

			       break;
		  case BSIDE: // внешняя нормаль T
			       // квадратичная интерполяция.
			  
			      
		          center_cord3D(border_neighbor[iG].iI, nvtx, pa, pp,TSIDE);
		          center_cord3D(border_neighbor[iG].iII, nvtx, pa, pb,TTSIDE);
			      center_cord3D(neighbors_for_the_internal_node[BSIDE][border_neighbor[iG].iII].iNODE1, nvtx, pa, pbb,BSIDE);
					
			      potent[border_neighbor[iG].iB]=my_quadratic_interpolation('+', potent[neighbors_for_the_internal_node[BSIDE][border_neighbor[iG].iII].iNODE1], potent[border_neighbor[iG].iII], potent[border_neighbor[iG].iI], pbb.z , pb.z, pp.z, pp.z+0.5*dz); 
				  
			       break;
		}

	}

	if (potentcopy != nullptr) {
		delete[] potentcopy;
	}

	for (integer i = 0; i < maxelm; i++) {
		delete[] avgonsosed[i];
	}
	delete[] avgonsosed;

} // flattener

// на основе 27 значений вокруг данного КО, с помощью квадратичной интерполяции
// находится значение функции в точке с координатами pzvezda.
/*
doublereal rinterpolFinposition(TOCHKA* c27, doublereal* potent27, TOCHKA pzvezda) {
	const integer dirP=0;
	const integer dirE=1;

} // rinterpolFinposition
*/

// данная небольшая функция используется с целью сокращения объёма кода
// в функции double_average_potent
void my_additional_to_calc_average(bool bflag, integer iNODE, doublereal rmultiplyer,
	                               integer** nvtx, TOCHKA* pa, doublereal* potent_in,
								   doublereal &rsvol, doublereal &rsum, doublereal &rcol, 
								   doublereal &rcolnumber) {

	// rmultiplyer - весовой коэффициент связанный с типом фильтра.

	if (!bflag&&(iNODE!=-1)) {
		doublereal dx, dy, dz;
        dx=0.0; dy=0.0; dz=0.0; // обнуляющий сброс значений
		volume3D(iNODE, nvtx, pa, dx, dy, dz);
		//rsvol+=dx*dy*dz;
		rsvol+=rmultiplyer*dx*dy*dz;
		//->//rsvol=1.0;
		rsum+=rmultiplyer*potent_in[iNODE]*dx*dy*dz;
		//-->//rsum+=rmultiplyer*potent_in[iNODE];
		rcol+=rmultiplyer;
		rcolnumber+=1.0;
    }

} // my_additional_to_calc_average

// получение дважды фильтрованного поля potent_test_filtr=test_filtr(potent_in),
// после применения тестового фильтра. Величиной potent_in может являться любая величина:
// компонента скорости, завихрённость и др.
// 16 мая 2012 года сделана квадратичная интерполяция изнутри области на границу.
// Работает только на структурированной прямоугольной ортогональной сетке.
// НЕ работает на АЛИС сетке.
void double_average_potent(doublereal* potent_in, doublereal* &potent_test_filtr, 
	integer maxelm, integer maxbound, ALICE_PARTITION** neighbors_for_the_internal_node,
							 integer** nvtx, TOCHKA* pa, doublereal* &delta_test_filtr,
							 integer itype_filtr, BOUND* border_neighbor, integer iquadraticinterpolboud) {

	// iquadraticinterpolboud - порядок интерполяции доступной для продолжения фильтрованной величины изнутри области на границу.
	// возможные значения 0 - просто снесение величины из ближайшего внутреннего КО  на границу. 2 - квадратичная интерполяция.

	// potent_in - на вход подаётся компонента скорости которая получена в результате 
	// решения уравнений Навье-Стокса. Это как бы один раз фильтрованная скорость
	// Box фильтром для которой является сам контрольный объём.
	// potent_test_filtr - выходная компонента скорости полученная из potent_in применением тестового
	// Box фильтра размеров в 27 контрольных объёмов окружающих данный контрольный объём вместе с данным.

    // delta_test_filtr размер тестового фильтра как корень кубический из
	// объёма ячейки (являющейся объединением 27 контрольных объёмов). 
	// содержит элементы от 0 .. maxelm-1. Вместо него в функцию может передаваться значение nullptr.

	// переменная itype_filtr - отвечает за тип фильтра: 0 - обычное среднее по объёму, 
	// 1 - фильтр на основе формулы трапеций, 2 - фильтр на основе формулы Симпсона.

	// инициализация.
	for (integer i=0; i<maxelm+maxbound; i++) potent_test_filtr[i]=0.0;

	// цикл по всем внутренним контрольным объёмам.
	// в цикле производится операция осреднения (тестовый фильтр).
	for (integer iP=0; iP<maxelm; iP++) {

		 // Алгоритм. Для каждого внутреннего КО может существовать не более 27 (включая его самого) объёмных (не граничных) соседей.
		 // Прежде чем вычислить среднюю величину которая рассчитывается просто как среднее арифметическое взвешенное на объём КО среди
		 // всех объёмных контрольных объёмов окружающих данный нужно определить какие из 26 возможных соседей действительно присутствуют.
		 // Определение присутствия соседей в условиях "произвольной" геометрии объёмная по коду операция. Когда все соседи известны получение 
		 // самой средней величины остаётся делом техники: вычисления суммы и деления её на общий объём. Т.е. это среднее по объёму значение.
		 // Вот смысл фильтрования - получение среднего. Сам метод контрольного объёма является фильтром - для величины хранящейся в центре КО 
		 // подразумевается что её значение распространяется на весь контрольный объём.

          bool bE, bW, bN, bS, bT, bB;
		  bE=false; bW=false; bN=false; bS=false; bT=false; bB=false;

		  integer iE, iN, iT, iW, iS, iB; // номера соседних контрольных объёмов
	      iE=neighbors_for_the_internal_node[ESIDE][iP].iNODE1; iN=neighbors_for_the_internal_node[NSIDE][iP].iNODE1; iT=neighbors_for_the_internal_node[TSIDE][iP].iNODE1; iW=neighbors_for_the_internal_node[WSIDE][iP].iNODE1; iS=neighbors_for_the_internal_node[SSIDE][iP].iNODE1; iB=neighbors_for_the_internal_node[BSIDE][iP].iNODE1;

          if (iE>=maxelm) bE=true; // если true то узел является граничным.
	      if (iN>=maxelm) bN=true;
	      if (iT>=maxelm) bT=true;
          if (iW>=maxelm) bW=true;
	      if (iS>=maxelm) bS=true;
	      if (iB>=maxelm) bB=true;

		  integer iEN=-1, iES=-1, iWS=-1, iWN=-1, iNT=-1, iNB=-1, iET=-1, iEB=-1, iST=-1, iSB=-1, iWT=-1, iWB=-1; // номера соседних Контрольных объёмов (12 шт).
          bool bEN, bES, bWS, bWN, bNT, bNB, bET, bEB, bST, bSB, bWT, bWB; 
		  // false означает что КО внутренний, а true что он граничный либо его вообще не существует.
		  bEN=false; bES=false; bWS=false; bWN=false; bNT=false; bNB=false;
		  bET=false; bEB=false; bST=false; bSB=false; bWT=false; bWB=false;

		  if (!bE && !bN) {
			  // внутренний КО.
			  iEN=neighbors_for_the_internal_node[NSIDE][neighbors_for_the_internal_node[ESIDE][iP].iNODE1].iNODE1;
			  if (iEN>=maxelm) bEN=true; // граничный ко.
		  } else {
			  bEN=true;
			  iEN=maxelm+maxbound; // KO не существует.
		  }

		  if (!bE && !bS) {
			  // внутренний КО.
			  iES=neighbors_for_the_internal_node[SSIDE][neighbors_for_the_internal_node[ESIDE][iP].iNODE1].iNODE1;
			  if (iES>=maxelm) bES=true; // граничный ко.
		  } else {
			  bES=true;
			  iES=maxelm+maxbound; // KO не существует.
		  }

		  if (!bW && !bS) {
			  // внутренний КО.
			  iWS=neighbors_for_the_internal_node[WSIDE][neighbors_for_the_internal_node[SSIDE][iP].iNODE1].iNODE1;
			  if (iWS>=maxelm) bWS=true; // граничный ко.
		  } else {
			  bWS=true;
			  iWS=maxelm+maxbound; // KO не существует.
		  }

		  if (!bW && !bN) {
			  // внутренний КО.
			  iWN=neighbors_for_the_internal_node[WSIDE][neighbors_for_the_internal_node[NSIDE][iP].iNODE1].iNODE1;
			  if (iWN>=maxelm) bWN=true; // граничный ко.
		  } else {
			  bWN=true;
			  iWN=maxelm+maxbound; // KO не существует.
		  }

		  if (!bN && !bT) {
			  // внутренний КО.
			  iNT=neighbors_for_the_internal_node[NSIDE][neighbors_for_the_internal_node[TSIDE][iP].iNODE1].iNODE1;
			  if (iNT>=maxelm) bNT=true; // граничный ко.
		  } else {
			  bNT=true;
			  iNT=maxelm+maxbound; // KO не существует.
		  }

		  if (!bN && !bB) {
			  // внутренний КО.
			  iNB=neighbors_for_the_internal_node[NSIDE][neighbors_for_the_internal_node[BSIDE][iP].iNODE1].iNODE1;
			  if (iNB>=maxelm) bNB=true; // граничный ко.
		  } else {
			  bNB=true;
			  iNB=maxelm+maxbound; // KO не существует.
		  }

		  if (!bE && !bT) {
			  // внутренний КО.
			  iET=neighbors_for_the_internal_node[ESIDE][neighbors_for_the_internal_node[TSIDE][iP].iNODE1].iNODE1;
			  if (iET>=maxelm) bET=true; // граничный ко.
		  } else {
			  bET=true;
			  iET=maxelm+maxbound; // KO не существует.
		  }

		  if (!bE && !bB) {
			  // внутренний КО.
			  iEB=neighbors_for_the_internal_node[ESIDE][neighbors_for_the_internal_node[BSIDE][iP].iNODE1].iNODE1;
			  if (iEB>=maxelm) bEB=true; // граничный ко.
		  } else {
			  bEB=true;
			  iEB=maxelm+maxbound; // KO не существует.
		  }

		  if (!bS && !bT) {
			  // внутренний КО.
			  iST=neighbors_for_the_internal_node[SSIDE][neighbors_for_the_internal_node[TSIDE][iP].iNODE1].iNODE1;
			  if (iST>=maxelm) bST=true; // граничный ко.
		  } else {
			  bST=true;
			  iST=maxelm+maxbound; // KO не существует.
		  }

		  if (!bS && !bB) {
			  // внутренний КО.
			  iSB=neighbors_for_the_internal_node[SSIDE][neighbors_for_the_internal_node[BSIDE][iP].iNODE1].iNODE1;
			  if (iSB>=maxelm) bSB=true; // граничный ко.
		  } else {
			  bSB=true;
			  iSB=maxelm+maxbound; // KO не существует.
		  }

		  if (!bW && !bT) {
			  // внутренний КО.
			  iWT=neighbors_for_the_internal_node[WSIDE][neighbors_for_the_internal_node[TSIDE][iP].iNODE1].iNODE1;
			  if (iWT>=maxelm) bWT=true; // граничный ко.
		  } else {
			  bWT=true;
			  iWT=maxelm+maxbound; // KO не существует.
		  }

		  if (!bW && !bB) {
			  // внутренний КО.
			  iWB=neighbors_for_the_internal_node[WSIDE][neighbors_for_the_internal_node[BSIDE][iP].iNODE1].iNODE1;
			  if (iWB>=maxelm) bWB=true; // граничный ко.
		  } else {
			  bWB=true;
			  iWB=maxelm+maxbound; // KO не существует.
		  }

		  integer iTNE=-1, iTSE=-1, iTNW=-1, iTSW=-1, iBNE=-1, iBSE=-1, iBNW=-1, iBSW=-1; // угловые точки куба (8 шт).
		  bool bTNE, bTSE, bTNW, bTSW, bBNE, bBSE, bBNW, bBSW;
          bTNE=true; bTSE=true; bTNW=true; bTSW=true; bBNE=true; bBSE=true; bBNW=true; bBSW=true; // по умолчанию данный ко не внутренний,
		  // он либо граничный либо его не существует.

		  if (!bEN) {
			  if (bTNE) {
				  // последовательность вызовов очень важна.
				  iTNE=neighbors_for_the_internal_node[TSIDE][neighbors_for_the_internal_node[NSIDE][neighbors_for_the_internal_node[ESIDE][iP].iNODE1].iNODE1].iNODE1;
			      bTNE=false;
			      if (iTNE>=maxelm) bTNE=true;
			  }
			  if (bBNE) {
				  // последовательность вызовов очень важна.
				  iBNE=neighbors_for_the_internal_node[BSIDE][neighbors_for_the_internal_node[NSIDE][neighbors_for_the_internal_node[ESIDE][iP].iNODE1].iNODE1].iNODE1;
				  bBNE=false;
			      if (iBNE>=maxelm) bBNE=true;
			  }
		  }

          if (!bES) {
			  if (bTSE) {
				  // последовательность вызовов очень важна.
				  iTSE=neighbors_for_the_internal_node[TSIDE][neighbors_for_the_internal_node[SSIDE][neighbors_for_the_internal_node[ESIDE][iP].iNODE1].iNODE1].iNODE1;
			      bTSE=false;
			      if (iTSE>=maxelm) bTSE=true;
			  }
			  if (bBSE) {
				  // последовательность вызовов очень важна.
				  iBSE=neighbors_for_the_internal_node[BSIDE][neighbors_for_the_internal_node[SSIDE][neighbors_for_the_internal_node[ESIDE][iP].iNODE1].iNODE1].iNODE1;
				  bBSE=false;
			      if (iBSE>=maxelm) bBSE=true;
			  }
		  }

		  if (!bWN) {
			  if (bTNW) {
				  // последовательность вызовов очень важна.
				  iTNW=neighbors_for_the_internal_node[TSIDE][neighbors_for_the_internal_node[NSIDE][neighbors_for_the_internal_node[WSIDE][iP].iNODE1].iNODE1].iNODE1;
			      bTNW=false;
			      if (iTNW>=maxelm) bTNW=true;
			  }
			  if (bBNW) {
				  // последовательность вызовов очень важна.
				  iBNW=neighbors_for_the_internal_node[BSIDE][neighbors_for_the_internal_node[NSIDE][neighbors_for_the_internal_node[WSIDE][iP].iNODE1].iNODE1].iNODE1;
				  bBNW=false;
			      if (iBNW>=maxelm) bBNW=true;
			  }
		  }

          if (!bWS) {
			  if (bTSW) {
				  // последовательность вызовов очень важна.
				  iTSW=neighbors_for_the_internal_node[TSIDE][neighbors_for_the_internal_node[SSIDE][neighbors_for_the_internal_node[WSIDE][iP].iNODE1].iNODE1].iNODE1;
			      bTSW=false;
			      if (iTSW>=maxelm) bTSW=true;
			  }
			  if (bBSW) {
				  // последовательность вызовов очень важна.
				  iBSW=neighbors_for_the_internal_node[BSIDE][neighbors_for_the_internal_node[SSIDE][neighbors_for_the_internal_node[WSIDE][iP].iNODE1].iNODE1].iNODE1;
				  bBSW=false;
			      if (iBSW>=maxelm) bBSW=true;
			  }
		  }

		  if (!bNT) {
			  if (bTNE) {
				  // последовательность вызовов очень важна.
				  iTNE=neighbors_for_the_internal_node[ESIDE][neighbors_for_the_internal_node[NSIDE][neighbors_for_the_internal_node[TSIDE][iP].iNODE1].iNODE1].iNODE1;
			      bTNE=false;
			      if (iTNE>=maxelm) bTNE=true;
			  }
			  if (bTNW) {
				  // последовательность вызовов очень важна.
				  iTNW=neighbors_for_the_internal_node[WSIDE][neighbors_for_the_internal_node[NSIDE][neighbors_for_the_internal_node[TSIDE][iP].iNODE1].iNODE1].iNODE1;
				  bTNW=false;
			      if (iTNW>=maxelm) bTNW=true;
			  }
		  }

		  if (!bNB) {
			  if (bBNE) {
				  // последовательность вызовов очень важна.
				  iBNE=neighbors_for_the_internal_node[ESIDE][neighbors_for_the_internal_node[NSIDE][neighbors_for_the_internal_node[BSIDE][iP].iNODE1].iNODE1].iNODE1;
			      bBNE=false;
			      if (iBNE>=maxelm) bBNE=true;
			  }
			  if (bBNW) {
				  // последовательность вызовов очень важна.
				  iBNW=neighbors_for_the_internal_node[WSIDE][neighbors_for_the_internal_node[NSIDE][neighbors_for_the_internal_node[BSIDE][iP].iNODE1].iNODE1].iNODE1;
				  bBNW=false;
			      if (iBNW>=maxelm) bBNW=true;
			  }
		  }

		  if (!bET) {
			  if (bTNE) {
				  // последовательность вызовов очень важна.
				  iTNE=neighbors_for_the_internal_node[NSIDE][neighbors_for_the_internal_node[TSIDE][neighbors_for_the_internal_node[ESIDE][iP].iNODE1].iNODE1].iNODE1;
			      bTNE=false;
			      if (iTNE>=maxelm) bTNE=true;
			  }
			  if (bTSE) {
				  // последовательность вызовов очень важна.
				  iTSE=neighbors_for_the_internal_node[SSIDE][neighbors_for_the_internal_node[TSIDE][neighbors_for_the_internal_node[ESIDE][iP].iNODE1].iNODE1].iNODE1;
				  bTSE=false;
			      if (iTSE>=maxelm) bTSE=true;
			  }
		  }

		  if (!bEB) {
			  if (bBNE) {
				  // последовательность вызовов очень важна.
				  iBNE=neighbors_for_the_internal_node[NSIDE][neighbors_for_the_internal_node[ESIDE][neighbors_for_the_internal_node[BSIDE][iP].iNODE1].iNODE1].iNODE1;
			      bBNE=false;
			      if (iBNE>=maxelm) bBNE=true;
			  }
			  if (bBSE) {
				  // последовательность вызовов очень важна.
				  iBSE=neighbors_for_the_internal_node[SSIDE][neighbors_for_the_internal_node[ESIDE][neighbors_for_the_internal_node[BSIDE][iP].iNODE1].iNODE1].iNODE1;
				  bBSE=false;
			      if (iBSE>=maxelm) bBSE=true;
			  }
		  }

          if (!bST) {
			  if (bTSE) {
				  // последовательность вызовов очень важна.
				  iTSE=neighbors_for_the_internal_node[ESIDE][neighbors_for_the_internal_node[SSIDE][neighbors_for_the_internal_node[TSIDE][iP].iNODE1].iNODE1].iNODE1;
			      bTSE=false;
			      if (iTSE>=maxelm) bTSE=true;
			  }
			  if (bTSW) {
				  // последовательность вызовов очень важна.
				  iTSW=neighbors_for_the_internal_node[WSIDE][neighbors_for_the_internal_node[SSIDE][neighbors_for_the_internal_node[TSIDE][iP].iNODE1].iNODE1].iNODE1;
				  bTSW=false;
			      if (iTSW>=maxelm) bTSW=true;
			  }
		  }

		  if (!bSB) {
			  if (bBSE) {
				  // последовательность вызовов очень важна.
				  iBSE=neighbors_for_the_internal_node[ESIDE][neighbors_for_the_internal_node[SSIDE][neighbors_for_the_internal_node[BSIDE][iP].iNODE1].iNODE1].iNODE1;
			      bBSE=false;
			      if (iBSE>=maxelm) bBSE=true;
			  }
			  if (bBSW) {
				  // последовательность вызовов очень важна.
				  iBSW=neighbors_for_the_internal_node[WSIDE][neighbors_for_the_internal_node[SSIDE][neighbors_for_the_internal_node[BSIDE][iP].iNODE1].iNODE1].iNODE1;
				  bBSW=false;
			      if (iBSW>=maxelm) bBSW=true;
			  }
		  }

		  if (!bWT) {
			  if (bTNW) {
				  // последовательность вызовов очень важна.
				  iTNW=neighbors_for_the_internal_node[NSIDE][neighbors_for_the_internal_node[WSIDE][neighbors_for_the_internal_node[TSIDE][iP].iNODE1].iNODE1].iNODE1;
			      bTNW=false;
			      if (iTNW>=maxelm) bTNW=true;
			  }
			  if (bTSW) {
				  // последовательность вызовов очень важна.
				  iTSW=neighbors_for_the_internal_node[SSIDE][neighbors_for_the_internal_node[WSIDE][neighbors_for_the_internal_node[TSIDE][iP].iNODE1].iNODE1].iNODE1;
				  bTSW=false;
			      if (iTSW>=maxelm) bTSW=true;
			  }
		  }

		  if (!bWB) {
			  if (bBNW) {
				  // последовательность вызовов очень важна.
				  iBNW=neighbors_for_the_internal_node[NSIDE][neighbors_for_the_internal_node[WSIDE][neighbors_for_the_internal_node[BSIDE][iP].iNODE1].iNODE1].iNODE1;
			      bBNW=false;
			      if (iBNW>=maxelm) bBNW=true;
			  }
			  if (bBSW) {
				  // последовательность вызовов очень важна.
				  iBSW=neighbors_for_the_internal_node[SSIDE][neighbors_for_the_internal_node[WSIDE][neighbors_for_the_internal_node[BSIDE][iP].iNODE1].iNODE1].iNODE1;
				  bBSW=false;
			      if (iBSW>=maxelm) bBSW=true;
			  }
		  }

          // Всё все возможные 26 кандидатов на соседи опознаны, можно приступать к осреднению.
		  // potent_in

		  doublereal rsum=0.0, rsvol=0.0, rcol=0.0, rcolnumber=0.0;
		  doublereal dx, dy, dz;

		  dx=0.0; dy=0.0; dz=0.0; // обнуляющий сброс значений
		  volume3D(iP, nvtx, pa, dx, dy, dz);
		  		  
		  //rsvol+=dx*dy*dz;	
		  switch (itype_filtr) {
		     case VOLUME_AVERAGE_FILTR: rcol+=1.0; rsum+=1.0*potent_in[iP]*dx*dy*dz; rsvol+=dx*dy*dz; break;
             case TRAPEZOIDAL_FILTR: rcol+=8.0; rsum+=8.0*potent_in[iP]*dx*dy*dz; rsvol+=8.0*dx*dy*dz; break;
	         case SIMPSON_FILTR: rcol+=64.0; rsum+=64.0*potent_in[iP]*dx*dy*dz; rsvol+=64.0*dx*dy*dz; break;
		  }
		  
		  /*
		  rsvol=1.0;
		  switch (itype_filtr) {
		     case VOLUME_AVERAGE_FILTR: rcol+=1.0; rsum+=1.0*potent_in[iP]; break;
             case TRAPEZOIDAL_FILTR: rcol+=8.0; rsum+=8.0*potent_in[iP]; break;
	         case SIMPSON_FILTR: rcol+=64.0; rsum+=64.0*potent_in[iP]; break;
		  }
		  */
		  rcolnumber+=1.0;

		  /*
		  // Этот и аналогичный код 
		  if (!bE) {
              dx=0.0; dy=0.0; dz=0.0; // обнуляющий сброс значений
		      volume3D(iE, nvtx, pa, dx, dy, dz);
		      rsvol+=dx*dy*dz;
		      rsum+=potent_in[iE]*dx*dy*dz; 
		  }
		  заменён на вызов функции
		   my_additional_to_calc_average(bE, iE, nvtx, pa, potent_in, rsvol, rsum);
		   // что выглядит гораздо компактнее.
		  */

		  switch (itype_filtr) {
			  case VOLUME_AVERAGE_FILTR:

				  // простое среднее по объёму.
				  
		           // сначала ближайшие соседи. (6шт)
		           my_additional_to_calc_average(bE, iE, 1.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           my_additional_to_calc_average(bW, iW, 1.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           my_additional_to_calc_average(bN, iN, 1.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           my_additional_to_calc_average(bS, iS, 1.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           my_additional_to_calc_average(bT, iT, 1.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           my_additional_to_calc_average(bB, iB, 1.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           // теперь соседи второго уровня (12шт).
		           my_additional_to_calc_average(bEN, iEN, 1.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           my_additional_to_calc_average(bES, iES, 1.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           my_additional_to_calc_average(bWS, iWS, 1.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           my_additional_to_calc_average(bWN, iWN, 1.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           my_additional_to_calc_average(bNT, iNT, 1.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           my_additional_to_calc_average(bNB, iNB, 1.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           my_additional_to_calc_average(bET, iET, 1.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           my_additional_to_calc_average(bEB, iEB, 1.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           my_additional_to_calc_average(bST, iST, 1.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           my_additional_to_calc_average(bSB, iSB, 1.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           my_additional_to_calc_average(bWT, iWT, 1.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           my_additional_to_calc_average(bWB, iWB, 1.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           // теперь самые дальние соседи (8шт).
		           my_additional_to_calc_average(bTNE, iTNE, 1.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           my_additional_to_calc_average(bTSE, iTSE, 1.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           my_additional_to_calc_average(bTNW, iTNW, 1.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           my_additional_to_calc_average(bTSW, iTSW, 1.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           my_additional_to_calc_average(bBNE, iBNE, 1.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           my_additional_to_calc_average(bBSE, iBSE, 1.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           my_additional_to_calc_average(bBNW, iBNW, 1.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           my_additional_to_calc_average(bBSW, iBSW, 1.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);  

				   break;

		      case TRAPEZOIDAL_FILTR: 
				  
				   // трапецевидный фильтр основанный на формуле трапеций.
				   // весовые коэффициенты получены для равномерной сетки.

				   // сначала ближайшие соседи. (6шт)
		           my_additional_to_calc_average(bE, iE, 4.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           my_additional_to_calc_average(bW, iW, 4.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           my_additional_to_calc_average(bN, iN, 4.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           my_additional_to_calc_average(bS, iS, 4.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           my_additional_to_calc_average(bT, iT, 4.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           my_additional_to_calc_average(bB, iB, 4.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           // теперь соседи второго уровня (12шт).
		           my_additional_to_calc_average(bEN, iEN, 2.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           my_additional_to_calc_average(bES, iES, 2.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           my_additional_to_calc_average(bWS, iWS, 2.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           my_additional_to_calc_average(bWN, iWN, 2.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           my_additional_to_calc_average(bNT, iNT, 2.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           my_additional_to_calc_average(bNB, iNB, 2.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           my_additional_to_calc_average(bET, iET, 2.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           my_additional_to_calc_average(bEB, iEB, 2.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           my_additional_to_calc_average(bST, iST, 2.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           my_additional_to_calc_average(bSB, iSB, 2.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           my_additional_to_calc_average(bWT, iWT, 2.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           my_additional_to_calc_average(bWB, iWB, 2.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           // теперь самые дальние соседи (8шт).
		           my_additional_to_calc_average(bTNE, iTNE, 1.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           my_additional_to_calc_average(bTSE, iTSE, 1.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           my_additional_to_calc_average(bTNW, iTNW, 1.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           my_additional_to_calc_average(bTSW, iTSW, 1.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           my_additional_to_calc_average(bBNE, iBNE, 1.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           my_additional_to_calc_average(bBSE, iBSE, 1.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           my_additional_to_calc_average(bBNW, iBNW, 1.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           my_additional_to_calc_average(bBSW, iBSW, 1.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);

				   break;

              case SIMPSON_FILTR:

				   // фильтр основанный на формуле Симпсона.
				   // весовые коэффициенты получены для равномерной сетки.

				   // сначала ближайшие соседи. (6шт)
		           my_additional_to_calc_average(bE, iE, 16.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           my_additional_to_calc_average(bW, iW, 16.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           my_additional_to_calc_average(bN, iN, 16.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           my_additional_to_calc_average(bS, iS, 16.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           my_additional_to_calc_average(bT, iT, 16.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           my_additional_to_calc_average(bB, iB, 16.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           // теперь соседи второго уровня (12шт).
		           my_additional_to_calc_average(bEN, iEN, 4.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           my_additional_to_calc_average(bES, iES, 4.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           my_additional_to_calc_average(bWS, iWS, 4.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           my_additional_to_calc_average(bWN, iWN, 4.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           my_additional_to_calc_average(bNT, iNT, 4.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           my_additional_to_calc_average(bNB, iNB, 4.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           my_additional_to_calc_average(bET, iET, 4.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           my_additional_to_calc_average(bEB, iEB, 4.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           my_additional_to_calc_average(bST, iST, 4.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           my_additional_to_calc_average(bSB, iSB, 4.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           my_additional_to_calc_average(bWT, iWT, 4.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           my_additional_to_calc_average(bWB, iWB, 4.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           // теперь самые дальние соседи (8шт).
		           my_additional_to_calc_average(bTNE, iTNE, 1.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           my_additional_to_calc_average(bTSE, iTSE, 1.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           my_additional_to_calc_average(bTNW, iTNW, 1.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           my_additional_to_calc_average(bTSW, iTSW, 1.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           my_additional_to_calc_average(bBNE, iBNE, 1.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           my_additional_to_calc_average(bBSE, iBSE, 1.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           my_additional_to_calc_average(bBNW, iBNW, 1.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber);
		           my_additional_to_calc_average(bBSW, iBSW, 1.0, nvtx, pa, potent_in, rsvol, rsum, rcol, rcolnumber); 

				   break;
		  }

		  if (delta_test_filtr!=nullptr) {
             delta_test_filtr[iP]=exp((1.0/3.0)*log(rsvol)); // ширина тестового фильтра в ячейке iP.
		  }

		  // Ради этой строчки писался весь код данной функции.
		  //potent_test_filtr[iP]=(rsum/rcol)*(rcolnumber/rsvol); // среднее по объёму значение.
          potent_test_filtr[iP]=(rsum/rsvol); // среднее по объёму значение.
		  // но для такого среднего значения лучше применять квадратичную интерполяцию на равномерную сетку с неравномерной.
		  //potent_test_filtr[iP]=(rsum/rcol);//*(rcolnumber); // просто среднее без учёта объёма.
	}

	
	if (iquadraticinterpolboud==2) {
		// нужно применить квадратичную интерполяцию для 
	// нахождения фильтрованной величины в граничных контрольных объёмах.
	for (integer iG=0; iG<maxbound; iG++) {
		// цикл по всем граничным контрольным объёмам.

		doublereal dx=0.0, dy=0.0, dz=0.0;// объём текущего контрольного объёма
		// вычисление размеров текущего контрольного объёма:
	    volume3D(border_neighbor[iG].iI, nvtx, pa, dx, dy, dz);
		TOCHKA pp,pb,pbb;

		// сканируем внутреннюю нормаль
		switch (border_neighbor[iG].Norm) {
		  case ESIDE: // внешняя нормаль W
			       // квадратичная интерполяция.
			      
				   center_cord3D(border_neighbor[iG].iI, nvtx, pa, pp, WSIDE);
				   center_cord3D(border_neighbor[iG].iII, nvtx, pa, pb, WW);
			       center_cord3D(neighbors_for_the_internal_node[ESIDE][border_neighbor[iG].iII].iNODE1, nvtx, pa, pbb, ESIDE);

				   potent_test_filtr[border_neighbor[iG].iB]=my_quadratic_interpolation('-', potent_test_filtr[neighbors_for_the_internal_node[ESIDE][border_neighbor[iG].iII].iNODE1], potent_test_filtr[border_neighbor[iG].iII], potent_test_filtr[border_neighbor[iG].iI], pbb.x , pb.x, pp.x, pp.x-0.5*dx);
				   
			       break;
		  case WSIDE: // внешняя нормаль E
			       // квадратичная интерполяция.

				   
		           center_cord3D(border_neighbor[iG].iI, nvtx, pa, pp, ESIDE);
		           center_cord3D(border_neighbor[iG].iII, nvtx, pa, pb, EE);
			       center_cord3D(neighbors_for_the_internal_node[WSIDE][border_neighbor[iG].iII].iNODE1, nvtx, pa, pbb, WSIDE);
					
			       potent_test_filtr[border_neighbor[iG].iB]=my_quadratic_interpolation('+', potent_test_filtr[neighbors_for_the_internal_node[WSIDE][border_neighbor[iG].iII].iNODE1], potent_test_filtr[border_neighbor[iG].iII], potent_test_filtr[border_neighbor[iG].iI], pbb.x , pb.x, pp.x, pp.x+0.5*dx);

			       break;
		  case NSIDE: // внешняя нормаль S
			       // квадратичная интерполяция.
			       
			       
		           center_cord3D(border_neighbor[iG].iI, nvtx, pa, pp, SSIDE);
		           center_cord3D(border_neighbor[iG].iII, nvtx, pa, pb, SS);
			       center_cord3D(neighbors_for_the_internal_node[NSIDE][border_neighbor[iG].iII].iNODE1, nvtx, pa, pbb, NSIDE);

			       potent_test_filtr[border_neighbor[iG].iB]=my_quadratic_interpolation('-', potent_test_filtr[neighbors_for_the_internal_node[NSIDE][border_neighbor[iG].iII].iNODE1], potent_test_filtr[border_neighbor[iG].iII], potent_test_filtr[border_neighbor[iG].iI], pbb.y , pb.y, pp.y, pp.y-0.5*dy);
			  
			       break;
		  case SSIDE:// внешняя нормаль N
			       // квадратичная интерполяция.

			       
		           center_cord3D(border_neighbor[iG].iI, nvtx, pa, pp, NSIDE);
		           center_cord3D(border_neighbor[iG].iII, nvtx, pa, pb,NN);
			       center_cord3D(neighbors_for_the_internal_node[SSIDE][border_neighbor[iG].iII].iNODE1, nvtx, pa, pbb, SSIDE);

			       potent_test_filtr[border_neighbor[iG].iB]=my_quadratic_interpolation('+', potent_test_filtr[neighbors_for_the_internal_node[SSIDE][border_neighbor[iG].iII].iNODE1], potent_test_filtr[border_neighbor[iG].iII], potent_test_filtr[border_neighbor[iG].iI], pbb.y , pb.y, pp.y, pp.y+0.5*dy);
			  
			       break;
		  case TSIDE: // внешняя нормаль B
			       // квадратичная интерполяция.

		           center_cord3D(border_neighbor[iG].iI, nvtx, pa, pp, BSIDE);
		           center_cord3D(border_neighbor[iG].iII, nvtx, pa, pb, BB);
			       center_cord3D(neighbors_for_the_internal_node[TSIDE][border_neighbor[iG].iII].iNODE1, nvtx, pa, pbb, TSIDE);


			       potent_test_filtr[border_neighbor[iG].iB]=my_quadratic_interpolation('-', potent_test_filtr[neighbors_for_the_internal_node[TSIDE][border_neighbor[iG].iII].iNODE1], potent_test_filtr[border_neighbor[iG].iII], potent_test_filtr[border_neighbor[iG].iI], pbb.z , pb.z, pp.z, pp.z-0.5*dz);

			       break;
		  case BSIDE: // внешняя нормаль T
			       // квадратичная интерполяция.
			  
			      
		          center_cord3D(border_neighbor[iG].iI, nvtx, pa, pp, TSIDE);
		          center_cord3D(border_neighbor[iG].iII, nvtx, pa, pb, TTSIDE);
			      center_cord3D(neighbors_for_the_internal_node[BSIDE][border_neighbor[iG].iII].iNODE1, nvtx, pa, pbb, BSIDE);
					
			      potent_test_filtr[border_neighbor[iG].iB]=my_quadratic_interpolation('+', potent_test_filtr[neighbors_for_the_internal_node[BSIDE][border_neighbor[iG].iII].iNODE1], potent_test_filtr[border_neighbor[iG].iII], potent_test_filtr[border_neighbor[iG].iI], pbb.z , pb.z, pp.z, pp.z+0.5*dz); 
				  
			       break;
		}

	}
	}
	else {
		for (integer iG=0; iG<maxbound; iG++) {
		   // цикл по всем граничным контрольным объёмам.
           potent_test_filtr[border_neighbor[iG].iB]=potent_test_filtr[border_neighbor[iG].iI];

	    }
	}


} // double_average_velocity

// Для каждого внутреннего Контрольного объёма возвращает значение true если значение угла между вихрем и 
// фильтрованным вихрем больше чем beta0==15 градусов. Используется в модели турбулентности Selective Smagorinsky (SS).
void calc_selective_smagorinsky(FLOW &f, bool* &bfibeta, integer itype_filtr, doublereal beta0) {

	// Описание данного кода такое.
	// Модель турбулентности Selective Smagorinsky предлагает не учитывать вклад турбулентной динамической вязкости
	// в том случае если угол между вихрем и осреднённым вихрем меньше порогового значения beta0. Значение beta0 установлено
	// в эксперименте и примерно составляет около 15 градусов. Данный приём избирательного учёта вклада турбулентной динамической 
	// вязкости позволяет избирательно уменьшить диссипативность модели Смагоринского.
	// Входные параметры: 
	// f - вся информация о гидродинамической подобласти,
	// beta0 - пороговое значение угла между вихрем и осреднённым вихрем в избирательной модели Смагоринского,
	// itype_filtr - один из трёх фильтров на выбор , например, SIMPSON_FILTR - фильтр типа Симпсона.
	// bfibeta - возвращаемое для всех строго внутренних контрольных объёмов булево значение, говорящее о том 
	// нужно ли учитывать вклад турбулентной динамической вязкости или нет.

	// Вспомогательный материал:
	// Формула векторного произведения:
	// [vectora x vectorb]=(ay*bz-az*by, az*bx-ax*bz, ax*by-ay*bx);

	// bfibeta - принимает значение true (турбулентную динамическую вязкость нужно учитывать) и значение false
	// если турбулентную динамическую вязкость учитывать ненужно. 
	// bfibeta имеет размерность 0..maxelm-1 т.е. имеет значение для всех внутренних контрольных объёмов.

	// Алгоритм такой:
	// 1. Взять три компоненты скорости.
	// 2. Зная частные производные компонент скорости из п.1 вычислить три компоненты вихря по формуле векторного произведения.
	// 3. Отфильтровать и запомнить отфильтрованные компоненты вихря.
	// 4. Вычислить модуль векторного произведения вихря на фильтрованный вихрь как корень из суммы квадратов компонент соответствующего вектора.
	// 5. Найти модуль вихря и модуль фильтрованного вихря. 
	// 6. Вычислить искомый угол в градусах применяя функцию арксинус.

	// VX - 0, VY - 1, VZ - 2

	
	// Компоненты Вихря.
	const integer iCURLX=0;
	const integer iCURLY=1;
	const integer iCURLZ=2;

	doublereal** curl_component = nullptr;
	curl_component = new doublereal*[3];
	for (integer i=0; i<3; i++) curl_component[i]=new doublereal[f.maxelm+f.maxbound];
	for (integer i=0; i<f.maxelm+f.maxbound; i++) {
		// инициализация .
		curl_component[iCURLX][i]=0.0;
		curl_component[iCURLY][i]=0.0;
		curl_component[iCURLZ][i]=0.0;
	}

	// компоненты Вихря.
    for (integer i=0; i<f.maxelm+f.maxbound; i++) {
		curl_component[iCURLX][i]=f.potent[GRADYVZ][i]-f.potent[GRADZVY][i];
		curl_component[iCURLY][i]=f.potent[GRADZVX][i]-f.potent[GRADXVZ][i];
	    curl_component[iCURLZ][i]=f.potent[GRADXVY][i]-f.potent[GRADYVX][i];
	}

	// Выделение оперативной памяти под отфильтрованные компоненты вихря.
	doublereal** curl_component_test_filtr = nullptr;
	curl_component_test_filtr = new doublereal*[3];
	for (integer i=0; i<3; i++) curl_component_test_filtr[i]=new doublereal[f.maxelm+f.maxbound];
	for (integer i=0; i<f.maxelm+f.maxbound; i++) {
		// инициализация .
		curl_component_test_filtr[iCURLX][i]=0.0;
		curl_component_test_filtr[iCURLY][i]=0.0;
		curl_component_test_filtr[iCURLZ][i]=0.0;
	}

	doublereal* delta_test_filtr;
	delta_test_filtr=nullptr;

	// Вычисление фильтрованных компонент вихря:
	// Внимание ! фильтр от производной некоторой величины не равен производной от этой фильтрованной величины.
	// Для вычисления значений на границе области применяется квадратичная интерполяция изнутри области наружу.
	double_average_potent(curl_component[iCURLX], curl_component_test_filtr[iCURLX], f.maxelm, f.maxbound, f.neighbors_for_the_internal_node, 
		                     f.nvtx, f.pa, delta_test_filtr, itype_filtr, f.border_neighbor,2);
	double_average_potent(curl_component[iCURLY], curl_component_test_filtr[iCURLY], f.maxelm, f.maxbound, f.neighbors_for_the_internal_node, 
							 f.nvtx, f.pa, delta_test_filtr, itype_filtr, f.border_neighbor,2);
	double_average_potent(curl_component[iCURLZ], curl_component_test_filtr[iCURLZ], f.maxelm, f.maxbound, f.neighbors_for_the_internal_node, 
							 f.nvtx, f.pa, delta_test_filtr, itype_filtr, f.border_neighbor,2);

	bfibeta=new bool[f.maxelm]; // выделение оперативной памяти под возвращаемое значение.
	for (integer i=0; i<f.maxelm; i++) bfibeta[i]=true; // по умолчанию турбулентную вязкость нужно везде учитывать.
	beta0=fabs(beta0);

	doublereal* curlxtest_curl = nullptr;
	curlxtest_curl = new doublereal[3];
	// вычисление значения возвращаемой булевой функции:
	for (integer i=0; i<f.maxelm; i++) {
		// для всех внутренних контрольных объёмов.
		doublereal module_Curl=0.0;
		// модуль вихря.
		module_Curl=sqrt(curl_component[iCURLX][i]*curl_component[iCURLX][i]+curl_component[iCURLY][i]*curl_component[iCURLY][i]+curl_component[iCURLZ][i]*curl_component[iCURLZ][i]);

        doublereal module_Curl_test_filtr=0.0;
		// модуль фильтрованного вихря.
		module_Curl_test_filtr=sqrt(curl_component_test_filtr[iCURLX][i]*curl_component_test_filtr[iCURLX][i]+curl_component_test_filtr[iCURLY][i]*curl_component_test_filtr[iCURLY][i]+curl_component_test_filtr[iCURLZ][i]*curl_component_test_filtr[iCURLZ][i]);

		
		curlxtest_curl[iCURLX]=curl_component[iCURLY][i]*curl_component_test_filtr[iCURLZ][i]-curl_component[iCURLZ][i]*curl_component_test_filtr[iCURLY][i];
        curlxtest_curl[iCURLY]=curl_component[iCURLZ][i]*curl_component_test_filtr[iCURLX][i]-curl_component[iCURLX][i]*curl_component_test_filtr[iCURLZ][i];
		curlxtest_curl[iCURLZ]=curl_component[iCURLX][i]*curl_component_test_filtr[iCURLY][i]-curl_component[iCURLY][i]*curl_component_test_filtr[iCURLX][i];
		doublereal module_curlxcurl=0.0;
        module_curlxcurl=sqrt(curlxtest_curl[iCURLX]*curlxtest_curl[iCURLX]+curlxtest_curl[iCURLY]*curlxtest_curl[iCURLY]+curlxtest_curl[iCURLZ]*curlxtest_curl[iCURLZ]);
		

		// при вычисление угла может произойти исключительная ситуация,
		// например, деление на ноль и т.п. Поэтому надо сообщить пользователю
		// о такой ошибке и выйти из приложения. 
		try
		{
			if ((module_Curl<1e-30) || (module_Curl_test_filtr<1e-30)) {
				// какие-то ничтожно малые завихрённости надо добавить
				// турбулизации выключив сглаживающую турбулентную вязкость.
                bfibeta[i]=false;
				// Эта ситуация может сложится на первой итерации когда вихрь везде равен нулю.
				//printf("warning little module Curl...\n");
				//getchar();
			}
			else {
				// аргумент арксинуса должен лежать от -1.0 до +1.0.
			    doublereal beta=180.0*fabs(asin(module_curlxcurl/(module_Curl*module_Curl_test_filtr)))/3.141;
			    if (beta>=beta0) {
				   // угол значителен и значит нужно учитывать турбулентную вязкость.
				   bfibeta[i]=true;
			    }
			    else {
				   // угол слишком маленький и турбулентную вязкость 
				   // учитывать ненужно.
				   bfibeta[i]=false;
			    }
			}
		}
		catch (integer a)
		{
			printf("Selective Smagorinsky angle beta exeption...\n");
#if doubleintprecision == 1
			printf("Caught exception number:  %lld\n", a);
#else
			printf("Caught exception number:  %d\n", a);
#endif

			
			printf("Please, press any key to exit...\n");
			//getchar();
			system("pause");
			exit(0);
		}
		
	}

	if (curlxtest_curl != nullptr) {
		delete[] curlxtest_curl;
	}
	// Освобождение оперативной памяти.
	// уничтожение вихря.
	if (curl_component != nullptr) {
		for (integer i = 0; i < 3; i++) {
			if (curl_component[i] != nullptr) {
				delete[] curl_component[i];
			}
		}
		delete[] curl_component;
	}
	// уничтожение отфильтрованного вихря.
	if (curl_component_test_filtr != nullptr) {
		for (integer i = 0; i < 3; i++) {
			if (curl_component_test_filtr[i] != nullptr) {
				delete[] curl_component_test_filtr[i];
			}
		}
		delete[] curl_component_test_filtr;
	}

} // calc_selective_smagorinsky


// begin 28 апреля 2012 года.
// Реализация модели турбулентности Германо.
// Dynamic subgrid-scale model 1991 год.
// Это относительно новая модель турбулентности 1991 год.
// Суть модели состоит в том, что она берёт стандартного классического
// Смагоринского (расчёты по которому проводились ещё в 1960 годах для метерологических
// исследований на сетках 20x20x20) и определяет значение его константы Cs заменяя её 
// квадрат на свою константу Cg - константа Германо. Значение константы Германо вычисляется
// на основе информации содержащейся в прослойке между обычным фильтром и тестовым фильтром.
// Под обычным фильтром понимается неявное фильтрование присущее методу контрольного объёма.
// См. книгу [1] А.А. Юн Теория и практика моделирования турбулентных течений с теплообменом, 
// смешением, химическими реакциями и двухфазных течений. Москва 2009 года.
// [2] Волков К.Н., Емельянов В.Н. Моделирование крупных вихрей в расчётах турбулентных течений.
// ФИЗМАТЛИТ, 2008.
void my_Germano_model(FLOW &f, doublereal* &Cs2, integer itype_filtr) {
	// данный метод возвращает квадрат константы Смагоринского.
	// Внимание значение квадрата может быть и отрицательным это тоже имеет физический
	// смысл как передача энергии по каскаду от мелких вихрей к крупным (т.е. в противоположную
	// сторону от естественного направления передачи энергии).

	// Значение Cs2 - определено только для внутренних контрольных объёмов.
	// размерность 0..f.maxelm-1.

	// внутри этой функции память для Cs2 выделяется автоматически.

	// Внимание нужно применить интерполяцию для расчёта граничных значений
	// фильтрованной величины.


	// VX 0 VY 1 VZ 2
	doublereal** speed_test_filtering=new doublereal*[3];
	for (integer i=0; i<3; i++) {
		speed_test_filtering[i]=new doublereal[f.maxelm+f.maxbound];
	}
	// инициализация
	for (integer i=0; i<3; i++) {
		#pragma omp parallel for
		for (integer j=0; j<f.maxelm+f.maxbound; j++) {
			speed_test_filtering[i][j]=0.0;  // инициализация.
		}
	}

	doublereal* delta_test_filtr;
	delta_test_filtr=nullptr;


	// Вычисление фильтрованных компонент скорости.
	// Внимание ! фильтр от производной некоторой величины не равен производной от этой фильтрованной величины.
	// Для вычисления значений на границе области применяется квадратичная интерполяция изнутри области наружу.
	double_average_potent(f.potent[VX], speed_test_filtering[VX], f.maxelm, f.maxbound, f.neighbors_for_the_internal_node, 
							 f.nvtx, f.pa, delta_test_filtr, itype_filtr, f.border_neighbor,2);
	double_average_potent(f.potent[VY], speed_test_filtering[VY], f.maxelm, f.maxbound, f.neighbors_for_the_internal_node, 
							 f.nvtx, f.pa, delta_test_filtr, itype_filtr, f.border_neighbor,2);
	double_average_potent(f.potent[VZ], speed_test_filtering[VZ], f.maxelm, f.maxbound, f.neighbors_for_the_internal_node, 
							 f.nvtx, f.pa, delta_test_filtr, itype_filtr, f.border_neighbor,2);

	

	// Для вычисления напряжений Леонарда нужны средние от попарных произведений компонент скорости.
	const integer VXVX=0;
	const integer VXVY=1;
	const integer VXVZ=2;
	const integer VYVX=3;
	const integer VYVY=4;
	const integer VYVZ=5;
	const integer VZVX=6;
	const integer VZVY=7;
	const integer VZVZ=8;

	doublereal** speedxspeed_test_filtering = nullptr;
	speedxspeed_test_filtering = new doublereal*[9];
	for (integer i=0; i<9; i++) {
		speedxspeed_test_filtering[i]=new doublereal[f.maxelm+f.maxbound];
	}
	// инициализация
	for (integer i=0; i<9; i++) {
		#pragma omp parallel for
		for (integer j=0; j<f.maxelm+f.maxbound; j++) {
			speedxspeed_test_filtering[i][j]=0.0;  // инициализация.
		}
	}

	doublereal* speedxspeed_buf = nullptr;
	speedxspeed_buf = new doublereal[f.maxelm + f.maxbound];
	// инициализация
	#pragma omp parallel for
	for (integer j=0; j<f.maxelm+f.maxbound; j++) {
		speedxspeed_buf[j]=0.0;
	}

	// Осреднение попарных произведений:
	integer** imarker=new integer*[9];
	for (integer i=0; i<9; i++) imarker[i]=new integer[2];
	imarker[VXVX][0]=VX; imarker[VXVX][1]=VX;
	imarker[VXVY][0]=VX; imarker[VXVY][1]=VY;
	imarker[VXVZ][0]=VX; imarker[VXVZ][1]=VZ;
	imarker[VYVX][0]=VY; imarker[VYVX][1]=VX;
	imarker[VYVY][0]=VY; imarker[VYVY][1]=VY;
	imarker[VYVZ][0]=VY; imarker[VYVZ][1]=VZ;
	imarker[VZVX][0]=VZ; imarker[VZVX][1]=VX;
	imarker[VZVY][0]=VZ; imarker[VZVY][1]=VY;
	imarker[VZVZ][0]=VZ; imarker[VZVZ][1]=VZ;

    for (integer i=0; i<9; i++) {
		integer fmr=imarker[i][0], smr=imarker[i][1]; // first marker, second marker

	    #pragma omp parallel for
		for (integer j=0; j<f.maxelm+f.maxbound; j++) {
		    speedxspeed_buf[j]=f.potent[fmr][j]*f.potent[smr][j];
	    }
		// Вычисление фильтрованных попарных произведений компонент скорости.
	    // Внимание ! фильтр от производной некоторой величины не равен производной от этой фильтрованной величины.
	    // Для вычисления значений на границе области применяется квадратичная интерполяция изнутри области наружу.
	    double_average_potent(speedxspeed_buf, speedxspeed_test_filtering[i], f.maxelm, f.maxbound, f.neighbors_for_the_internal_node, 
							 f.nvtx, f.pa, delta_test_filtr, itype_filtr, f.border_neighbor,2);
	}

	// Вычисление напряжений Леонарда.
	doublereal** rLeonard_stress = nullptr;
	rLeonard_stress = new doublereal*[9];
	for (integer i=0; i<9; i++) {
		rLeonard_stress[i]=new doublereal[f.maxelm+f.maxbound];
	}
	// непосредственное вычисление напряжений Леонарда (1974):
	// см. формулу (2.20) на 43 странице книги А.А. Юна см. [1].
	for (integer i=0; i<9; i++) {
		integer fmr=imarker[i][0], smr=imarker[i][1]; // first marker, second marker

		#pragma omp parallel for 
		for (integer j=0; j<f.maxelm+f.maxbound; j++) {
			// Это определения взято из книги Юна. 
			// Оно согласовано с определение константы Германо (в её формуле определения не нужен минус,
			// который можно встретить в других статьях.
			rLeonard_stress[i][j]=speed_test_filtering[fmr][j]*speed_test_filtering[smr][j]-speedxspeed_test_filtering[i][j];  
		}
	}
	// Освобождение памяти из под фильтрованных
	// попарных произведений скорости:
	if (speedxspeed_test_filtering != nullptr) {
		for (integer i = 0; i < 9; i++) {
			if (speedxspeed_test_filtering[i] != nullptr) {
				delete[] speedxspeed_test_filtering[i];
			}
		}
		delete[] speedxspeed_test_filtering;
	}

	if (speed_test_filtering != nullptr) {
		for (integer i = 0; i < 3; i++) {
			if (speed_test_filtering[i] != nullptr) {
				delete[] speed_test_filtering[i];
			}
		}
		delete[] speed_test_filtering;
	}

    // нужно вычислить компоненты тензора скоростей деформаций на основе 
	// уже вычисленных значений частных производных от компонент скорости по теореме Грина Гаусса.
	// Потом компоненты тензора скоростей деформаций нужно осреднить.
	// Первые 9 значений это компоненты тензора скоростей деформаций.
	// последние значение с индексом 9 это модуль тензора скоростей деформаций.
	doublereal** StRT = nullptr;
	StRT = new doublereal*[10]; // Strain Rate Tensor
	for (integer i=0; i<10; i++) StRT[i]=new doublereal[f.maxelm+f.maxbound];
	const integer MODULEStRt=9; 

	// Вычисление компонент тензора скорости деформации.

	// цикл по всем контрольным объёмам
    #pragma omp parallel for
	for (integer i=0; i<f.maxelm+f.maxbound; i++) {
		StRT[VXVX][i]=f.potent[GRADXVX][i];
		StRT[VXVY][i]=0.5*(f.potent[GRADYVX][i]+f.potent[GRADXVY][i]);
		StRT[VXVZ][i]=0.5*(f.potent[GRADZVX][i]+f.potent[GRADXVZ][i]);

		StRT[VYVX][i]=0.5*(f.potent[GRADXVY][i]+f.potent[GRADYVX][i]);
		StRT[VYVY][i]=f.potent[GRADYVY][i];
		StRT[VYVZ][i]=0.5*(f.potent[GRADZVY][i]+f.potent[GRADYVZ][i]);

		StRT[VZVX][i]=0.5*(f.potent[GRADXVZ][i]+f.potent[GRADZVX][i]);
		StRT[VZVY][i]=0.5*(f.potent[GRADYVZ][i]+f.potent[GRADZVY][i]);
		StRT[VZVZ][i]=f.potent[GRADZVZ][i];

		// Модуль тензора скоростей деформаций:
		StRT[MODULEStRt][i]=sqrt(2.0*(StRT[VXVX][i]*StRT[VXVX][i]+
			                     StRT[VYVY][i]*StRT[VYVY][i]+
								 StRT[VZVZ][i]*StRT[VZVZ][i])+
								 4.0*(StRT[VXVY][i]*StRT[VXVY][i]+
								 StRT[VXVZ][i]*StRT[VXVZ][i]+
								 StRT[VZVY][i]*StRT[VZVY][i]));
	}

	if (speedxspeed_buf != nullptr) {
		delete[] speedxspeed_buf;
		speedxspeed_buf = nullptr;
	}
	doublereal* Mij_buf = nullptr;
	Mij_buf = new doublereal[f.maxelm + f.maxbound];

	doublereal** Mij=new doublereal*[9]; // сначала составляющая Mij которая вычитается справа.
	for (integer i=0; i<9; i++) Mij[i]=new doublereal[f.maxelm+f.maxbound];
	for (integer i=0; i<9; i++) {
		for (integer j=0; j<f.maxelm+f.maxbound; j++) {
			Mij[i][j]=0.0; // инициализация
		}
	}


	for (integer i=0; i<9; i++) {

		for (integer j=0; j<f.maxelm+f.maxbound; j++) {
			Mij_buf[j]=StRT[MODULEStRt][j]*StRT[i][j];
		}

        // Вычисление фильтрованных попарных произведений записанных в Mij_buf.
	    // Внимание ! фильтр от производной некоторой величины не равен производной от этой фильтрованной величины.
	    // Для вычисления значений на границе области применяется квадратичная интерполяция изнутри области наружу.
	    double_average_potent(Mij_buf, Mij[i], f.maxelm, f.maxbound, f.neighbors_for_the_internal_node, 
							 f.nvtx, f.pa, delta_test_filtr, itype_filtr, f.border_neighbor,2);

	}

	if (Mij_buf != nullptr) {
		delete[] Mij_buf;
	}

	doublereal** StRT_filtr = nullptr;
	StRT_filtr = new doublereal*[10]; // дважды фильтрованный Strain Rate Tensor
	for (integer i=0; i<10; i++) StRT_filtr[i]=new doublereal[f.maxelm+f.maxbound];
	for (integer i=0; i<10; i++) for (integer j=0; j<f.maxelm+f.maxbound; j++) StRT_filtr[i][j]=0.0; // инициализация


    for (integer i=0; i<9; i++) {

		// Вычисление дважды отфильтрованных компонент Strain Rate Tensor`а.
	    // Внимание ! фильтр от производной некоторой величины не равен производной от этой фильтрованной величины.
	    // Для вычисления значений на границе области применяется квадратичная интерполяция изнутри области наружу.
	    double_average_potent(StRT[i], StRT_filtr[i], f.maxelm, f.maxbound, f.neighbors_for_the_internal_node, 
							 f.nvtx, f.pa, delta_test_filtr, itype_filtr, f.border_neighbor,2);
	}

	delta_test_filtr=new doublereal[f.maxelm]; // ширина тестового фильтра.
	for (integer i=0; i<f.maxelm; i++) delta_test_filtr[i]=0.0;

    // Вычисление дважды отфильтрованных компонент Strain Rate Tensor`а.
	// Внимание ! фильтр от производной некоторой величины не равен производной от этой фильтрованной величины.
	// Для вычисления значений на границе области применяется квадратичная интерполяция изнутри области наружу.
	double_average_potent(StRT[MODULEStRt], StRT_filtr[MODULEStRt], f.maxelm, f.maxbound, f.neighbors_for_the_internal_node, 
							 f.nvtx, f.pa, delta_test_filtr, itype_filtr, f.border_neighbor,2);

	if (StRT != nullptr) {
		for (integer i = 0; i < 10; i++) {
			if (StRT[i] != nullptr) {
				delete[] StRT[i];
			}
		}
		delete[] StRT;
	}

	// Наконец-то вычисление Mij:
	for (integer i=0; i<9; i++) {
		for (integer j=0; j<f.maxelm; j++) {
			doublereal dx=0.0, dy=0.0, dz=0.0; // обнуляющий сброс значений
		    volume3D(j, f.nvtx, f.pa, dx, dy, dz);
		    doublereal rsvol=dx*dy*dz;
			doublereal delta2=exp((2.0/3.0)*log(rsvol));
			// этот вариант даёт не соответствующее физическому смыслу решение на стыке двух равномерных сеток.
			// для равномерной сетки тестовый фильтр ровно в три раза больше чем базовый фильтр.
			// Mij[i][j]=(delta_test_filtr[j]*delta_test_filtr[j])*StRT_filtr[MODULEStRt][j]*StRT_filtr[i][j]-delta2*Mij[i][j];
			// Квадрат ширины тестового фильтра в 9 раз больше ширины квадрата базового фильтра.
			// Этот вариант в теории должен лучше себя вести на стыке двух равномерных сеток.
			// По результатам численного моделирования не оказалось разницы между двумя вариантами.
			// Опять на стыке двух равномерных сеток наблюдается не соответствующее физическому смыслу 
			// распределение постоянной Смагоринского
            Mij[i][j]=delta2*(9.0*StRT_filtr[MODULEStRt][j]*StRT_filtr[i][j]-Mij[i][j]);
			// Путь выхода из ситуации по видимому состоит в написании более адекватной процедуры фильтрования.
			// Суть более адекватной процедуры фильтрования состоит в том что нужно окружить фильтруемый объект
			// (рассматриваемый текущий контрольный объём) 26 одинаковыми с ним контрольными объёмами. Таким образом
			// будет получено 26 дополнительных геометрических позиций в который нужно восстановить фильтруемую функцию
			// Восстановление функции с неравномерной сетки на равномерную можно сделать с помощью квадратичной интерполяции.
			// Сделав всё это мы будем иметь как бы результат который полностью буде отвечать равномерной сетке и не соответствующего
			// физическому смыслу поведения на стыке равномерных сеток наблюдаться не будет.
			// Резюме: требуется переписать операцию фильтрования.
		}
	}

	// Освобождение оперативной памяти.

	if (StRT_filtr != nullptr) {
		for (integer i = 0; i < 10; i++) {
			if (StRT_filtr[i] != nullptr) {
				delete[] StRT_filtr[i];
			}
		}
		delete[] StRT_filtr;
	}

	if (delta_test_filtr != nullptr) {
		delete[] delta_test_filtr;
	}

	Cs2=new doublereal[f.maxelm];
	for (integer i=0; i<f.maxelm; i++) Cs2[i]=0.0; // инициализация.

	// Cs2 == Cg - константа Германо, квадрат константы Смагоринского.
	for (integer i=0; i<f.maxelm; i++) {
		doublereal rLijMij=0.0;
		for (integer j=0; j<9; j++) rLijMij+=rLeonard_stress[j][i]*Mij[j][i];
		doublereal rMijMij=0.0;
		for (integer j=0; j<9; j++) rMijMij+=Mij[j][i]*Mij[j][i];
		bool bzero=false; // если true то константа Германо равна нулю так всё поле скорости скорее всего нулевое.

		if (fabs(rMijMij)<1e-23) {
			if (fabs(rLijMij)<1e-23) {
				bzero=true;
			}
			else {
				printf("Error in Germano Dynamic Model: Mij is equal 0.0...\n");
			    printf("division by zero! ... \n");
			    printf("Please, press any key to exit...");
			    //getchar();
				system("pause");
			    exit(0);
			}
		}

		// Замечание: значение константы Германо может сильно колебаться
		// в пространстве и времени, что может привести к численной неустойчивости.
		// Предлагается либо применить к ней осреднение, либо применить к ней нижнюю 
		// релаксацию. 
		// В некоторых работах предлагается установить рамки для константы Смагоринского.
		// В [3] предлагается взять 0.06 <= Cs <= 0.25. Это даёт для квадрата константы 
		// Смагоринского значение: 0.0036 <= константа Германо <= 0.0625. Применение 
		// данного ограничителя обеспечивает диссипативность модели турбулентности.
		// [3]. У.С. Абдибеков, Н.Б. Усенбаев, О.Л.Каруна Численное моделирование турбулентного
		// течения в канале. Казахский национальный университет им. аль-Фараби, Алматы.
		// e-mail: uali:kazsu.kz


		// Согласно книге Юна здесь не нужен минус, 
		// смотри определение напряжений Леонарда.
		// Эти две формулы должны быть согласованы.
		if (!bzero) {
			Cs2[i]=rLijMij/(2.0*rMijMij);

			// Обеспечиваем необходимую диссипативность модели.
		    //Cs2[i]=fmin(0.0036, fmax(Cs2[i], 0.0625)); 
		    if (f.smaginfo.bLimiters_Cs) {
				// ограничиваем константу Смагоринского:
			    doublereal rfmax=f.smaginfo.maxCs*f.smaginfo.maxCs;
			    if (f.smaginfo.maxCs<0.0) rfmax*=-1.0;
			    doublereal rfmin=f.smaginfo.minCs*f.smaginfo.minCs;
			    if (f.smaginfo.minCs<0.0) rfmin*=-1.0;
			    if (rfmax>rfmin) {
                    //Cs2[i]=fmax( rfmin, fmin(Cs2[i], rfmax)); 
			    }
			    else {
				    printf("improper restriction of the constant Smagorinsky...\n");
				    printf("Please, press any key to exit is programm...\n");
				   // getchar();
					system("pause");
				    exit(0);
			    }
		    }
		    // если никаких ограничений нет то используется вычисленное значение квадрата постоянной Смагоринского.
		}
		else {
			// и числитель и знаменатель равны нулю, значит
			// имеем нулевое поле скорости, а значит нету никакой
			// турбулентности.
			Cs2[i]=0.0; 
		}
		
	}

	// В результате одной численной проверки было выяснено,
	// что сглаживатель не исправляет положение.
	// Приведём возможную причину ошибки.
	// Дано. Стыкуются в плоскости две равномерные сетки.
	// у одной шаг по нормали к стенке a у другой b.
	// Пусть для определённости b>>a. Тогда имеем отношение 
	// ширины базового фильтра к ширине тестового фильтра слева 
	// направо для четырёх последовательных контрольных объёмов.
	// Пусть для полной определённости b==2.0*a;
	// тогда: a/(3.0*a)==1.0/3.0, // для равномерной сетки ширины a.
	// a/(2.0*a+b)=a/(4.0*a)==1.0/4.0
	// b/(2.0*b+a)=2.0/(4.0+1.0)=2.0/5.0
	// b/(3.0*b)=1.0/3.0 // для равномерной сетки.
	// Итого:
	// a      a  < b   b 
	// 0.333 0.25 0.4 0.333
	// мы видим всплеск на стыке сеток.
	/*
	doublereal* Cs_flat=new doublereal[f.maxelm+f.maxbound];
	for (integer i=0; i<f.maxelm+f.maxbound; i++) {
		if (i<f.maxelm) {
	     	Cs_flat[i]=Cs2[i];
		}
		else Cs_flat[i]=0.0;
	}
	// Применим сглаживатель применение которого необоснованно.
	// Сглаживатель лучше применять до лимитирующих ограничений.
	flattener(Cs_flat, f.maxelm, f.maxbound, f.neighbors_for_the_internal_node, f.nvtx, f.pa, f.border_neighbor);
	for (integer i=0; i<f.maxelm; i++) Cs2[i]=Cs_flat[i];
	delete Cs_flat;
	*/

	// Освобождение оперативной памяти.
	
	// освобождаем память из под напряжений Леонарда.
	if (rLeonard_stress != nullptr) {
		for (integer i = 0; i < 9; i++) {
			if (rLeonard_stress[i] != nullptr) {
				delete[] rLeonard_stress[i];
			}
		}
		delete[] rLeonard_stress;
	}

	if (imarker != nullptr) {
		for (integer i = 0; i < 9; i++) {
			if (imarker[i] != nullptr) {
				delete[] imarker[i];
			}
		}
		delete[] imarker;
	}

	if (Mij != nullptr) {
		for (integer i = 0; i < 9; i++) {
			if (Mij[i] != nullptr) {
				delete[] Mij[i];
			}
		}
		delete[] Mij;
	}

	

} // модель Германо.


#endif