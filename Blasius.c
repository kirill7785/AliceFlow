// файл Blasius.c 
// 1. вычисление толщины гидродинамического пограничного слоя для задачи Блазиуса 1908 года (Boundary layer thickness);
// 2. Вычисление безразмерного локального коэффициента трения на пластине (Wall skin friction distribution);
// 3. Вычисление толщины вытеснения на поверхности пластины (displacement thickness);
// 4. Вычисление толщины температурного пограничного слоя (Thermal boundary layer distribution);
// 5. Вычисление распределеия числа Нуссельта по поверхности пластины.
// begin : 31 января 2012 года.
// end first : 5 февраля 2012 года.

#pragma once
#ifndef MY_BLASIUS_C
#define MY_BLASIUS_C 1

// Постпроцессинг для задачи Блазиуса 1908 года.
void boundarylayer_info(FLOW* &f, TEMPER &t, integer flow_interior_count, WALL* w, integer lw) {

	// Критическое число Рейнольдса при обтекании плоской пластины равно Re_кр=1e+5;

	// Предполагается что задача свелась к двумерной в плоскости YZ, 
	// а скорость VX равна нулю.
	// Т.к. программа оперирует трёхмерными объектами то нас интересует именно плоскость
	// находящаяся по центру оси X : avgX=0.5*(minX+maxX);

	// Пластина совпадает с осью Y при Z==0.0. Пластина начинается при Y==0.0 и длится 1 метр.

	const doublereal length_plate=1.0; // m
	doublereal U_inf=0.0; // скорость набегающего потока на бесконечности.
	doublereal Twall=-1e3; // Температура нагретой пластины
	doublereal Tinf=1e5; // Температура воздуха на входе.

	// определение скорости набегающего потока на бесконечности:
	for (integer i=0; i<lw; i++) {
		if ((!w[i].bsymmetry)&&(!w[i].bpressure)) {
			// стенка не является ни границей симметрии ни выходной границей.
			U_inf=fmax(U_inf,w[i].Vy); // т.к. U_inf совпадает по направлению с положительным направлением оси OY.
			Twall=fmax(Twall,w[i].Tamb); // находим температуру нагретой пластины.
			Tinf=fmin(Tinf,w[i].Tamb); // находим температуру охлаждающего пластину набегающего потока.
		}
	}

	// Пластина предполагается нагретой относительно температуры набегающего воздуха, т.е. Twall>Tinf.

    // Теперь в U_inf - хранится скорость набегающего потока на бесконечности.
	doublereal rho_avg=0.0; // средняя плотность по объёму.
	doublereal mu_avg=0.0; // средняя по объёму динамическая вязкость.
	doublereal lam_avg=0.0; // средняя по объёму теплопроводность
	doublereal cp_avg=0.0; // средняя по объёму теплоёмкость
	doublereal volume_default_interior=0.0; // Объём расчётной области. 
	doublereal VyMAX=-1e30; // максимум по расчётной области компоненты скорости параллельной пластине

	for (integer iP=0; iP<f[0].maxelm; iP++) {
		// вычисление размеров текущего контрольного объёма:
	    doublereal dx=0.0, dy=0.0, dz=0.0;// объём текущего контроольного объёма
	    volume3D(iP, f[0].nvtx, f[0].pa, dx, dy, dz);

		rho_avg+=dx*dy*dz*f[0].prop[RHO][iP];
		mu_avg+=dx*dy*dz*f[0].prop[MU][iP];
		lam_avg+=dx*dy*dz*t.prop[LAM][f[0].ptr[iP]];
		cp_avg+=dx*dy*dz*t.prop[HEAT_CAPACITY][f[0].ptr[iP]];
		volume_default_interior+=dx*dy*dz;

		VyMAX=fmax(VyMAX,f[0].potent[VY][iP]);
	}

	rho_avg=rho_avg/volume_default_interior; // средняя по объёму плотность среды.
	mu_avg=mu_avg/volume_default_interior; // средняя по объёму динамическая вязкость.
	lam_avg=lam_avg/volume_default_interior; // средняя по объёму теплопроводность.
	cp_avg=cp_avg/volume_default_interior; // средняя по объёму теплоёмкость.

	doublereal avg_Re_number=U_inf*length_plate*rho_avg/mu_avg; // Число Рейнольдса по длине пластинки.
	doublereal avg_Pr_number=mu_avg*cp_avg/lam_avg; // Число Прандтля.
	doublereal avg_Pe_number=avg_Re_number*avg_Pr_number; // Число Пекле.

	// нахождение контрольного объёма максимально близкого к точке avgX=0.5*(minX+maxX);
	doublereal dx=0.0, dy=0.0, dz=0.0;// объём текущего контроольного объёма
	TOCHKA p; // координаты центра текущего котрольного объёма.
	integer iP=0;
	doublereal minX, maxX, avgX;
	while (f[0].sosedi[WSIDE][iP].iNODE1<f[0].maxelm) iP=f[0].sosedi[WSIDE][iP].iNODE1;
    // вычисление размеров текущего контрольного объёма:
	volume3D(iP, f[0].nvtx, f[0].pa, dx, dy, dz);
	center_cord3D(iP, f[0].nvtx, f[0].pa, p,100); // вычисление координат центра КО.
	minX=p.x-0.5*dx;
	while (f[0].sosedi[ESIDE][iP].iNODE1<f[0].maxelm) iP=f[0].sosedi[ESIDE][iP].iNODE1;
	volume3D(iP, f[0].nvtx, f[0].pa, dx, dy, dz);
	center_cord3D(iP, f[0].nvtx, f[0].pa, p,100); // вычисление координат центра КО.
	maxX=p.x+0.5*dx;
	avgX=0.5*(minX+maxX); // координаты центральной плоскости в которой и будет происходить измерение.

	center_cord3D(iP, f[0].nvtx, f[0].pa, p,100); // вычисление координат центра КО.
	doublereal mindist=fabs(p.x-avgX);
	integer iPC=iP;
	while (f[0].sosedi[WSIDE][iP].iNODE1<f[0].maxelm) {
		iP=f[0].sosedi[WSIDE][iP].iNODE1;
		center_cord3D(iP, f[0].nvtx, f[0].pa, p,100); // вычисление координат центра КО.
		if (fabs(p.x-avgX)<mindist) {
			iPC=iP;
			mindist=fabs(p.x-avgX);
		}
	}

	// контрольный объём iPC наиболее приближен к центральной плоскости в которой будут производится замеры.
	iP=iPC;

	// Вычисление количества контрольных объёмов расположенных по длине пластины :
	integer iclength=0;
	while (f[0].sosedi[SSIDE][iP].iNODE1<f[0].maxelm) iP=f[0].sosedi[SSIDE][iP].iNODE1;

	while (iP<f[0].maxelm) {
		center_cord3D(iP, f[0].nvtx, f[0].pa, p,100); // вычисление координат центра КО.
		if ((p.y>0.0) && (p.y<1.0)) iclength++; // пластина расположена между 0.0 m и 1.0 m.
		iP=f[0].sosedi[NSIDE][iP].iNODE1;
	}

	// Гидродинамический пограничный слой: (Boundary layer thickness).
	doublereal* delta=new doublereal[iclength];
	// соответствующее расстояние от передней кромки пластины:
	doublereal* yposition=new doublereal[iclength];
	// Безразмерная величина локальный коэффициент трения на поверхности пластины :
	// Cx[i]=SInvariantStrainRateTensor[i]/(0.5*rho_avg*U_inf*U_inf); См. Лыков.
	// Wall skin friction distribution
	doublereal* Cx=new doublereal[iclength];
	// Средний коэффициент трения :
	doublereal avg_Cx=0.0; // integral(Cx[i]*dlength_plate)/length_plate;
	// Толщина вытеснения на поверхности пластины. (displacement thickness) 
	// см. Л.Д.Ландау, Е.М.Лифшиц Теоретическая физика. ГИДРОДИНАМИКА том VI. стр. 228.
	doublereal* displacement_thickness=new doublereal[iclength];
	// температурный пограничный слой (Thermal boundary layer distribution):
	doublereal* deltaT=new doublereal[iclength];
	// локальное число Вильгельма Нуссельта на стенке:
	// Число Нуссельта показывает отношение интенсивности теплообмена за счёт конвекции к
	// к интенсивности теплообмена за счёт теплопроводности.
	doublereal* local_Nusselt_number=new doublereal[iclength];


	// Вычисление 95% толщины гидродинамического пограничного слоя на пластине:
	doublereal deltascal=0.95; // 95% пограничный слой.
	doublereal deltaTscal=0.05; // 5% тепловой пограничный слой на нагретой пластине.
	iP=iPC;
	while (f[0].sosedi[SSIDE][iP].iNODE1<f[0].maxelm) iP=f[0].sosedi[SSIDE][iP].iNODE1; // перемотка в начало.

	integer ilengthcounter=0;
	for (ilengthcounter=0; ilengthcounter<iclength; ilengthcounter++) {
		// инициализация.
		delta[ilengthcounter]=0.0; // гидродинамический пограничный слой 
		yposition[ilengthcounter]=0.0; // абсцисса для графиков
		Cx[ilengthcounter]=0.0; // безразмерный локальный коэффициент трения на стенке
        displacement_thickness[ilengthcounter]=0.0; // толщина вытеснения
		deltaT[ilengthcounter]=0.0; // температурный пограничный слой
		local_Nusselt_number[ilengthcounter]=0.0; // локальное число Нуссельта.
	}
	ilengthcounter=0;

	while (iP<f[0].maxelm) {
		center_cord3D(iP, f[0].nvtx, f[0].pa, p,100); // вычисление координат центра КО.
		if ((p.y>0.0) && (p.y<length_plate)) {
			// пластина расположена между 0.0 m и length_plate m.
			yposition[ilengthcounter]=p.y; // абсцисса для толщины пограничного слоя.
			
			integer iBT=iP;
			while (f[0].sosedi[BSIDE][iBT].iNODE1<f[0].maxelm) iBT=f[0].sosedi[BSIDE][iBT].iNODE1; // перемотка в начало пластины.
			// вычисление распределения безразмерного локального коэффициента трения на стенке:
			// см. А.В.Лыков  Тепломассообмен справочник Москва., "Энергия", 1978г. стр. 184.
			Cx[ilengthcounter]=mu_avg*f[0].SInvariantStrainRateTensor[f[0].sosedi[BSIDE][iBT].iNODE1]/(0.5*rho_avg*U_inf*U_inf);

			// Вычисление толщины гидродинамического пограничного слоя.
			doublereal VYB, VYT, zB, zT;
			VYB=f[0].potent[VY][f[0].sosedi[BSIDE][iBT].iNODE1];
			VYT=f[0].potent[VY][iBT];
			center_cord3D(iBT, f[0].nvtx, f[0].pa, p, BSIDE); // вычисление координат центра КО.
			volume3D(iBT, f[0].nvtx, f[0].pa, dx, dy, dz);
			zB=p.z-0.5*dz;
			zT=p.z;
			while ((f[0].sosedi[TSIDE][iBT].iNODE1<f[0].maxelm) && (VYT<deltascal*U_inf)) {
				iBT=f[0].sosedi[TSIDE][iBT].iNODE1; // удаляемся от пластины перпендикулярно её плоскости.
				VYB=f[0].potent[VY][f[0].sosedi[BSIDE][iBT].iNODE1];
			    VYT=f[0].potent[VY][iBT];
				zB=zT;
				center_cord3D(iBT, f[0].nvtx, f[0].pa, p, TSIDE); // вычисление координат центра КО.
				zT=p.z;
			}
			doublereal a, b;
			a=(VYT-VYB)/(zT-zB);
			b=(zT*VYB-VYT*zB)/(zT-zB);
			delta[ilengthcounter]=(deltascal*U_inf-b)/a; // искомая толщина гидродинамического пограничного слоя.

			while (f[0].sosedi[BSIDE][iBT].iNODE1<f[0].maxelm) iBT=f[0].sosedi[BSIDE][iBT].iNODE1; // перемотка в начало пластины.
			// Вычисление толщины вытеснения по скорости Uoperating (displacement thickness):
			doublereal Uoperating=VyMAX; // U_inf - по скорости набегающего потока на бесконечности. (другой вариант по максимальной скорости - VyMAX)
			while (f[0].sosedi[TSIDE][iBT].iNODE1<f[0].maxelm) {
				volume3D(iBT, f[0].nvtx, f[0].pa, dx, dy, dz);
			    displacement_thickness[ilengthcounter]+=dz*(Uoperating-f[0].potent[VY][iBT]);
                iBT=f[0].sosedi[TSIDE][iBT].iNODE1; // удаляемся от пластины перпендикулярно её плоскости.
			}
			displacement_thickness[ilengthcounter]/=Uoperating; // толщина вытеснения.

			// Температурный пограничный слой.
			while (f[0].sosedi[BSIDE][iBT].iNODE1<f[0].maxelm) iBT=f[0].sosedi[BSIDE][iBT].iNODE1; // перемотка в начало пластины.
			doublereal TempB, TempT;
			iBT=f[0].ptr[iBT];
			TempB=t.potent[t.sosedi[BSIDE][iBT].iNODE1];
			TempT=t.potent[iBT];
			center_cord3D(iBT, t.nvtx, t.pa, p, BSIDE); // вычисление координат центра КО.
			volume3D(iBT, t.nvtx, t.pa, dx, dy, dz);
			zB=p.z-0.5*dz;
			zT=p.z;

			// Вычисляем локальное Число Нуссельта:
			local_Nusselt_number[ilengthcounter]=(t.prop_b[LAM][t.sosedi[BSIDE][iBT].iNODE1-t.maxelm]*((Twall-TempT)/(0.5*dz))*yposition[ilengthcounter])/((Twall-Tinf)*lam_avg);

			// Продолжение вычисления температурного пограничного слоя.
			doublereal Temp_critical=Tinf+deltaTscal*(Twall-Tinf);
			while ((t.sosedi[TSIDE][iBT].iNODE1<t.maxelm)&&(t.potent[iBT]>Temp_critical)) {
				iBT=t.sosedi[TSIDE][iBT].iNODE1; // удаляемся от пластины перпендикулярно её плоскости.
				TempB=t.potent[t.sosedi[BSIDE][iBT].iNODE1];
			    TempT=t.potent[iBT];
				zB=zT;
                center_cord3D(iBT, t.nvtx, t.pa, p, TSIDE); // вычисление координат центра КО.
				zT=p.z;
			}
			a=(TempB-TempT)/(zB-zT);
			b=(TempT*zB-TempB*zT)/(zB-zT);
			deltaT[ilengthcounter]=(Temp_critical-b)/a; // толщина температурного пограничного слоя.

            ilengthcounter++;
		}
		iP=f[0].sosedi[NSIDE][iP].iNODE1;
	}

	// Средний безразмерный коэффициент трения на стенке:
	// интеграл считается по формуле трапеций.
	avg_Cx=0.0;
	doublereal slen=0.0;
	for (integer i=0; i<(iclength-1); i++) {
		slen+=(yposition[i+1]-yposition[i]);
		avg_Cx+=0.5*(yposition[i+1]-yposition[i])*(Cx[i+1]+Cx[i]);
	}
	avg_Cx=avg_Cx/slen; // средний безразмерный коэффициент трения на стенке.

	// Печать подготовленной информации в текстовый файл blasius_1908.txt
	FILE *fpblas; // файл в который будет записываться информация о задаче Блазиуса.
	errno_t err_blas;
	err_blas = fopen_s(&fpblas, "blasius_1908.txt", "w");

	if ((err_blas) != 0) {
		 printf("Create File blasius_1908.txt Error\n");
         //getchar();
		 system("pause");
         //exit(0);
     }
	 else {

		 if (fpblas != NULL) {
			 // Запись информации в текстовый файл:
			 fprintf(fpblas, "Laminar Flow and Heat Transfer over a Flat Plate.\n\n");
			 // свойства материалов:
			 fprintf(fpblas, "Average material properties: \n");
			 fprintf(fpblas, "Density= %+.16f kg/m!3\n", rho_avg);
			 fprintf(fpblas, "Dynamic_Viscosity= %+.16f Pa*s\n", mu_avg);
			 fprintf(fpblas, "Thermal_conductivity= %+.16f W/(m*K)\n", lam_avg);
			 fprintf(fpblas, "Specific_Heat= %+.16f J/(kg*K)\n\n", cp_avg);
			 // геометрические размеры:
			 fprintf(fpblas, "Geometric_parametrs: \n");
			 fprintf(fpblas, "length_plate= %+.16f m\n\n", length_plate);
			 // граничные условия:
			 fprintf(fpblas, "Boundary_conditions: \n");
			 fprintf(fpblas, "Inlet_fluid_velocity= %+.16f m/s\n", U_inf);
			 fprintf(fpblas, "Inlet_fluid_temperature= %+.16f oC\n", Tinf);
			 fprintf(fpblas, "Plate_Wall_temperature= %+.16f oC\n\n", Twall);
			 // остальные скалярные параметры:
			 fprintf(fpblas, "Other_scalar_parametrs: \n");
			 fprintf(fpblas, "Critical_Renolds_number= %+.16f\n", 1e5);
			 fprintf(fpblas, "Average_Renolds_number= %+.16f\n", avg_Re_number);
			 fprintf(fpblas, "Average_Prandtl_number= %+.16f\n", avg_Pr_number);
			 fprintf(fpblas, "Average_Peclet_number= %+.16f\n", avg_Pe_number);
			 fprintf(fpblas, "Local_average_coefficient_of_friction_at_the_wall= %+.16f \n\n", avg_Cx);
			 // XY plot :
			 fprintf(fpblas, "XY Plots: \n");
			 fprintf(fpblas, "ypos - y_position m; \n");
			 fprintf(fpblas, "delta - 95%%_Boundary_layer_thickness m; \n");
			 fprintf(fpblas, "Cx - Wall_skin_friction_distribution; \n");
			 fprintf(fpblas, "delta* - displacement_thickness m; \n");
			 fprintf(fpblas, "deltaT - 5%%_Thermal_boundary_layer_distribution m; \n");
			 fprintf(fpblas, "local_Nusselt_number - Wall_Nusselt_number_distribution; \n\n");
			 // Собственно сами числовые данные для построения графиков:
			 fprintf(fpblas, "ypos	delta	Cx	delta*	deltaT	local_Nusselt_number\n");
			 for (integer i = 0; i < iclength; i++) {
				 fprintf(fpblas, "%+.16f %+.16f %+.16f %+.16f %+.16f %+.16f\n", yposition[i], delta[i], Cx[i], displacement_thickness[i], deltaT[i], local_Nusselt_number[i]);
			 }

			 fclose(fpblas); // закрытие файла для записи постпроцессинга в задаче Блазиуса 1908 года.
		 }
	 }


	// Освобождение оперативной памяти:

	// Гидродинамические характеристики
	delete[] delta; 
	delete[] yposition;
	delete[] Cx;
	delete[] displacement_thickness;
	// Тепловые характеристики
	delete[] deltaT;
	delete[] local_Nusselt_number;

} // boundarylayer_info

#endif