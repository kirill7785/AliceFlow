// Файл icepak_unstructured_mesh.cpp
// Hexa unstructured mesh по мотивам unsys icepak
// Многоблочная структурированная сетка (используется в программе sigmaFlow Гаврилова Андрея.).
// начало 26.06.2021

#pragma once
#ifndef _ICEPAK_UNSTRUCTURED_MESH_CPP_
#define _ICEPAK_UNSTRUCTURED_MESH_CPP_ 1


void write_unstructured_mesh(doublereal* &xpos, doublereal* &ypos, doublereal* &zpos, int** nvtx,
	int maxnod, int maxelm) {


	FILE* fp = NULL;

#ifdef MINGW_COMPILLER
	int err = 0;
	fp = fopen64("ALICEFLOW0_57.PLT", "w");
	if (fp == NULL) err = 1;
#else
	errno_t err = 0;
	err = fopen_s(&fp, "ALICEFLOW0_57.PLT", "w");
#endif


	// создание файла для записи.
	if ((err) != 0) {
		printf("Create File Error\n");
		exit(1);
	}
	else {
		if (fp != NULL) {
			// запись заголовка
			fprintf(fp, "TITLE = \"ALICEFLOW0_57\"\n");

			// запись имён переменных
			fprintf(fp, "VARIABLES = x, y, z\n");

			// запись информации о зонах

			//if (nve==3) fprintf(fp, "ZONE T=\"Rampant\", N=%lld, E=%lld, ET=TRIANGLE, F=FEBLOCK\n\n", maxelm, ncell);
			//if (nve==4) fprintf(fp, "ZONE T=\"Rampant\", N=%lld, E=%lld, ET=QUADRILATERAL, F=FEBLOCK\n\n", maxelm, ncell);
			fprintf(fp, "ZONE T=\"Rampant\", N=%lld, E=%lld, ET=BRICK, F=FEBLOCK\n\n", maxnod, maxelm);

			// запись x
			for (int i = 0; i < maxnod; i++) {
				fprintf(fp, "%e ", xpos[i]));
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			fprintf(fp, "\n");

			// запись y
			for (int i = 0; i < maxnod; i++) {
				fprintf(fp, "%e ", ypos[i]));
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			fprintf(fp, "\n");

			// запись z
			for (int i = 0; i < maxnod; i++) {
				fprintf(fp, "%e ", zpos[i]));
				if (i % 10 == 0) fprintf(fp, "\n");
			}

			fprintf(fp, "\n");

			// запись nvtx
			for (int i = 0; i < maxelm; i++) {
				fprintf(fp, "%d %d %d %d ", nvtx[0][i], nvtx[1][i], nvtx[2][i], nvtx[3][i]);
				fprintf(fp, "%d %d %d %d\n", nvtx[4][i], nvtx[5][i], nvtx[6][i], nvtx[7][i]);
			}

			fprintf(fp, "\n");

			fclose(fp);
		}
	}

	printf("Mesh succsefull writing\n");
	getchar();
	exit();
} // write_mesh



void unstructured_grid(	integer lb, integer ls, integer lw,
	BLOCK*& b, SOURCE*& s, WALL*& w, integer lu, 
	UNION*& my_union, TPROP* matlist,
	integer& iunion_id_p1)
{
	doublereal* xpos=nullptr, *ypos=nullptr,  *zpos=nullptr;
	int** nvtx = nullptr;

	int maxnod = 0;
	int maxelm = 0;


	// Требования 27.06.2021.
	// Граница сетки в точности повторяет форму тел. Сетка скошенная не ортогональная.
	// Стыковка на границах объектов не обязательно узел в узел, хотя если это возможно сделать то желательно узел в узел.

	// Алгоритм. 27.06.2021.
	// A. Разить всю область на непересекающиеся тела достаточно произвольной формы каждое.
	// 	   При этом нужно выяснить с какими телами пересекается исследуемое тело и скоректировать его границу с помощью
	//  разбивки на еще более простые примитивы.
	// 	   Нужно  формировать стек тел.
	// Тело произвольной формы представить набором простых тел:
	// 1. Скошенная шестигранная призма.
	// 2. В простом случае просто прямая прямоугльная шестигранная призма.
	// 3. Цилиндр.
	// 4. Цилиндрическое кольцо.
	// 5. Внешняя часть цилиндра дополненная до прямой прямоугольной призмы. 
	// B. Рабить простые тела структурированной сеткой.

	// C. Решатель температурного поля по аналогии с графовым методом Future.


	write_unstructured_mesh(xpos, ypos, zpos, nvtx,
		 maxnod, maxelm);
	
}

#endif