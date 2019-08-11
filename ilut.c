// Файл ilut.c
// ILUT preconditioner (предобуславливатель)
// начало 3 июля 2011.
// 31 марта 2013 года. Уже переведён первый том книги Ю.Саада на русский.
// Здесь производится попытка интеграции его библиотеки SPARSKIT2 в код AliceFlowv0_07
// с помощью замечательной программы f2c.exe скачанной в интернете.

#include "iluk.c" // ILUK decomposition

// quick-sort split
// 31 марта 2013 года.
void qsplit(doublereal* a, integer* ind, integer n, integer ncut) {

    /* Локальные переменные */
	doublereal tmp, abskey;
	integer itmp, first, last;
	integer mid, j;


    /* ----------------------------------------------------------------------- */
    /*     does a quick-sort split of a real array. */
    /*     on input a(1:n). is a real array */
    /*     on output a(1:n) is permuted such that its elements satisfy: */

    /*     abs(a(i)) .ge. abs(a(ncut)) for i .lt. ncut and */
    /*     abs(a(i)) .le. abs(a(ncut)) for i .gt. ncut */

    /*     ind(1:n) is an integer array which permuted in the same way as a(*). */
    /* ----------------------------------------------------------------------- */
    /* ----- */
    /*
    Примерный перевод вышенаписанного:
    делает быстрое упорядочивание вещественного массива a.
    На входе вещественный массив a содержащий информацию в позициях 1:n включительно.
    На выходе переупорядоченный массив a элементы которого удовлетворяют следующему условию:

    fabs(a[i])>=fabs(a[ncut]) для всех i < ncut а также
    fabs(a[i])<=fabs(a[ncut]) для всех i > ncut.

    ind это целочисленный массив который хранит историю перестановок элементов массива a.
    ind тоже имеет хранение информации в ячейках 1:n включительно.

    Особое внимание следует обратить, что все массивы начинаются с 1 и заканчиваются значением n !!!. 
    */


    // Следующих трёх строк не было в моём коде, они появились после применения f2c.exe.
    /* В языке СИ указатель массива указывает на первый элемент индексация которого начинаетсясо значения ноль.
    Таким образом если мы сдвинем указатель начала на одну позицию влево то первый элемент, который раньше имел индекс ноль
    теперь будет иметь индекс 1 и следовательно последний элемент будет иметь индес n. Значит можно один в один применять код Ю. Саада.
    В конце тока надо вернуть всё на свои места, т.е. опять чтобы индексация начиналась с нуля. Для этого есть обратная операция ++ind; ++a;
    */
    /* Parameter adjustments */
    --ind; // Внимание важно, в этом отличии может заключаться причина неработоспособности.
    --a; // Отличие тока в этих трёх строках, остальное выверено и перепроверено многократно.

	first=1;
	last=n;
	if (((ncut >= first) && (ncut <= last))) {
		do { // внешний цикл пока mid не равен ncut делать.
			mid=first;
			abskey=fabs(a[mid]);
			for (j=first+1; j<=last; j++) {
				if (fabs(a[j])>abskey) {
					mid++;
                    // меняем местами значения в позициях mid и j.
				    tmp=a[mid];
				    itmp=ind[mid];
				    a[mid]=a[j];
				    ind[mid]=ind[j];
				    a[j]=tmp;
				    ind[j]=itmp;
			    }
			}

            // меняем местами позиции mid и first:

			tmp=a[mid];
			a[mid]=a[first];
			a[first]=tmp;

			itmp=ind[mid];
			ind[mid]=ind[first];
			ind[first]=itmp;

            // проверяем условие окончания цикла.


			if (mid > ncut) {
				last=mid-1;
			}
			else if (mid < ncut) {
				first=mid+1;
			}

		} while (mid != ncut);
	}

    // Возвращение на индексацию в стиле СИ.
    // теперь первый элемент снова имеет индекс 0, а последний индекс n-1; !!!
    // А массив переупорядочен и индексы ind тоже !.
    ++ind;
    ++a; 
	
} // qsplit


// Эта версия ilut вручную переведена с фортрана на язык СИ.
// Проблема в том, что я не уверен в правильности перевода.
// Причина сложностей в том, что у Саада всё начинается с единицы в то время как в Си все массивы начинаются с нуля.
// 31 марта 2013 года было принято руководствоваться машинным переводом с помощью утилиты f2c.exe.
void ilut(integer n, doublereal* &a, integer* &ja, integer* &ia, 
		  integer lfil, doublereal droptol, integer &iwk, 
		  doublereal* &alu, integer* &jlu, integer* &ju,
		  doublereal* &w, integer* &jw, integer &ierr)
{
	/*
	* Входные данные:
	* n - размер матрицы СЛАУ или число уравнений.
	* a,ja,ia - матрица СЛАУ в CRS формате. С 
	*          симметричным портретом.
	* lfil - максимальное число внедиагональных
	*        элементов в каждой строке L и U.
	*        lfil должно быть >= 0.
	* droptol - пороговое значение для отсева.
	* iwk - максимальный размер alu и jlu.
	*/

	/*
	* Возвращаемые значения: 
	* alu, jlu - матрица в модифицированном Sparse Row формате (MSR).
	* ju - вектор длины n, указывающий позиции начала каждой строки U
	*      в матрицах alu,jlu.
	* ierr - целочисленный код ошибки. Возвращает 0 при успешном завершении.
	*
	* Модифицированная строчная схема (Modified Sparse Row, MSR) использует 
	* только два массива : массив alu - значений элементов матрицы LU (на 
	* диагонали L стоят единицы) и целочисленный массив jlu. В первых n 
	* позициях alu содержит диагональные элементы матрицы U по порядку. 
	* Элемент alu[n+1] не заполняется или несёт дополнительную информацию о матрице.
	* Начиная с позиции n+2 записываются ненулевые элементы исходной матрицы по строкам,
	* исключая диагональные. Для каждого элемента alu[k] элемент jlu[k] показывает 
	* столбцовый индекс в исходной матрице.На n+1 позициях матрицы jlu размещаются указатели 
	* входа для каждой строки матрицы alu и jlu.
	*/

	/*
	*  Рабочие массивы:
	*  jw - целочисленнный рабочий длины 2*n.
	*  w - вещественный рабочий длины n+1.
	*      1..ii-1 L part, ii..n U part.
	*/

	/*
	* Численный алгоритм вариационной инициализации океанологических полей. 
	* Русаков Александр Сергеевич. инфа о параметрах.
	*
	* 1) Любой элемент L и U величина которого меньше некоторого 
	*    порогового значения (связанного с нормой рассматриваемой строки) отбрасывается. 
	* 2) Сохранение только крупнейших lfil элементов в
	*    i-ой строке L и U (за исключением диагональных элементов).
	* 
	* Гибкость стратегии:
	* можно использовать droptol==0.0 чтобы
	* получить стратегии основанные на хранении
	* наибольших элементов в каждой строке L && U.
	* принимая же droptol!=0.0 но lfil==n получаем обычную
	* стратегию пороовога для отсева элементов.
	*/

	// локальные переменные:
	integer ju0, k, j1, j2, j, ii, i, lenl, lenu, jj, jrow, jpos, len;
	doublereal tnorm, t,  s, fact; // absolute,

    /* Parameter adjustments */
    --jw;
    --w;
    --ju;
    --ia;
    --a;
    --ja;
    --alu;
    --jlu;


	if (lfil < 0) {
		// параметр lfil задан неправильно.
		printf("lfil incorrect detected.../n");
        //getchar();
		system("pause");
	}
	else {
       ju0 = n+2;
	   jlu[1] = ju0;

	   for (j=1; j<=n; j++) jw[n+j]=0; // инициализация массива ненулевых индексов.

	   // начало main цикла
	   for (ii=1; ii<=n; ii++) {
		   j1 = ia[ii]; // начало строки ii
		   j2 = ia[ii+1]-1; // конец строки ii
           tnorm=0.0; 
		   for (k=j1; k<=j2; k++) tnorm+=fabs(a[k]);
		   if (tnorm < 1e-30) {
			   // в строке все элементы равны нулю
#if doubleintprecision == 1
			   printf("%lld is zero string in matrix.../n", ii - 1);
#else
			   printf("%d is zero string in matrix.../n", ii - 1);
#endif
			  
               //getchar();
			   system("pause");
			   break; // досрочный выход из цикла for
		   }
		   else {
			   tnorm/=(doublereal)(j2-j1+1);
               
			   // распаковка частей матриц L и U 
			   // соответствующих строке матрицы A  
			   // в вещественный массив w.
			   lenu=1;
			   lenl=0;
			   jw[ii]=ii;
			   w[ii]=0.0;
			   jw[n+ii]=ii;

			   for (j=j1; j<=j2; j++) {
				   k=ja[j];
				   t=a[j];
				   if (k<ii) {
					   lenl++;
					   jw[lenl]=k;
					   w[lenl]=t;
					   jw[n+k]=lenl;
				   }
				   else if (k==ii) w[ii]=t;
				   else {
					   lenu++;
					   jpos=ii+lenu-1;
					   jw[jpos]=k;
					   w[jpos]=t;
					   jw[n+k]=jpos;
				   }
			   }
			   jj=0;
			   len=0;

			   // Удаление предыдущих строк
			   bool bweShouldbeContinue=true;
			   while (bweShouldbeContinue) {
				   jj++; //FORTRAN 150

			       if (jj<=lenl) {

				       jrow=jw[jj];
				       k=jj;

				       for (j=jj+1; j<=lenl; j++) {
					       if (jw[j] < jrow) {
						       jrow=jw[j];
						       k=j;
					       }
				       }

				       if (k!=jj) {
					       j=jw[jj];
					       jw[jj]=jw[k];
					       jw[k]=j;

					       jw[n+jrow]=jj;
					       jw[n+j]=k;
                           // swap(w[jj],w[k]);
					       s=w[jj];
					       w[jj]=w[k];
					       w[k]=s;
				       }

				       jw[n+jrow]=0;

					   // Здесь возможно нужно не умножать а делить
                       fact=w[jj]*alu[jrow];
				   
				       if (fabs(fact) > droptol) { 

					      for (k=ju[jrow]; k<=jlu[jrow+1]-1; k++) {
						      s=fact*alu[k];
						      j=jlu[k];
						      jpos=jw[n+j];
						      if (j >= ii) {
							     if (jpos == 0) {
								     lenu++;
								     if (lenu > n) {
									      // Матрица задана неправильно.
										  // крушение приложения.
										  printf(" Error!!! Matrix caused by incorrect.../n");
										  //getchar();
										  system("pause");
									      exit(0);
								     }
								     i=ii+lenu-1;
								     jw[i]=j;
								     jw[n+j]=i;
								     w[i]=-s;
							     }
							     else
							     {
								    w[jpos]-= s;
							     }
						       }
						       else 
						       {
							       if (jpos == 0) {
								       lenl++;
								       if (lenl > n) {
									       // Матрица задана неправильно.
										  // крушение приложения.
										  printf(" Error!!! Matrix caused by incorrect.../n");
										  //getchar();
										  system("pause");
									      exit(0);
								       }
								       jw[lenl]=j;
								       jw[n+j]=lenl;
								       w[lenl]= -s;
							       }
							        else
							       {
								       w[jpos]-=s;
							       }
						       }
					      } // end for k

					      len++;
					      w[len]=fact;
					      jw[len]=jrow;
				       }

			       }
			       else {
				   
				      bweShouldbeContinue=false;

				      // установить двойной указатель на 0.
				      for (k=1; k<=lenu; k++) jw[n+jw[ii+k-1]]=0;

				      // обновить матрицу L
				      lenl=len;
				      if (lenl < lfil) len=lenl; 
				      else len=lfil;

				      // sort by quick-split
				      //qsplit(w,jw,lenl,len);
                      qsplit(&w[1],&jw[1],lenl,len);

				      for (k=1; k<=len; k++) {
					      if (ju0 > iwk) {
						       // ошибка. нехватает места для матрицы LU разложения.
						       printf("error! there is no plase for the matrix LU.../n");
						      // getchar();
							   system("pause");
					    	   exit(0);
					      }
					      else {
						      alu[ju0]=w[k];
						      jlu[ju0]=jw[k];
						      ju0++;
					      }
				      }

				      // сохраним указатель на начало строки ii в матрице U
				      ju[ii]=ju0;

				      // обновление матрицы U - первое применение пороговой стратегии.
				      len=0;
				      for (k=1; k<=lenu-1; k++) {
					      if (fabs(w[ii+k])>droptol*tnorm) {
						      len++;
						      w[ii+len]=w[ii+k];
						      jw[ii+len]=jw[ii+k];
					      }
				      }
				      lenu=len+1;
                      // len=min(lenu,lfil);
				      if (lenu<lfil) len=lenu;
				      else len=lfil;

				      // В оригинале непонятный вызов функции
				      //qsplit(w[ii+1], jw[ii+1], lenu-1, len);
				      // Здесь он заменён на 
				      //qsplit(w,jw,lenu-1,len);
                      // Нижеследующая строка получена путём автоматической трансляции с фортрана на Си с помощью f2c.exe.
                      qsplit(&w[ii+1],&jw[ii+1],lenu-1, len);

				      t=fabs(w[ii]);
				      if ((len+ju0)> iwk) {
                           // ошибка. нехватает места для матрицы LU разложения.
					    	printf("error! there is no plase for the matrix LU.../n");
						   // getchar();
							system("pause");
						    exit(0);
				      }
				      for (k=ii+1; k<=ii+len-1; k++) {
					       jlu[ju0]=jw[k];
					       alu[ju0]=w[k];
					       t+=fabs(w[k]);
					       ju0++;
				      }

                      // Если на диагонали ноль, то мы делаем чтобы был не ноль.
                      // величина 1e-30 импирическая и её значение требует проверки.
                      if (fabs(w[ii])<1e-30) {
                          w[ii]=(droptol + 1.0e-4)*tnorm;
					  }

				      // В оригинальной работе присваивается 
				      // элемент обратный диагональному.
				      alu[ii]=1.0/w[ii]; 

				      // Обновим указатель на начало строки в матрице U
				      jlu[ii+1]=ju0;
			       }
			   } // while 
		   }

	   }  // конец main цикла

	}

    /* Parameter adjustments */
    ++jw;
    ++w;
    ++ju;
    ++ia;
    ++a;
    ++ja;
    ++alu;
    ++jlu;

} // ilut

// msrcsr  : converts modified sparse row format to compressed sparse   
//           row format.   
// ilut выдаёт матрицу в msr формате, а рабочий формат crs.
/* ----------------------------------------------------------------------- */
/* Subroutine */ void msrcsr_(integer n, doublereal *a, integer *ja, 
	doublereal *ao, integer *jao, integer *iao, doublereal *wk, integer *iwk)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    integer i__, j, k, ii, iptr;
    bool added;
    integer idiag;

/* ----------------------------------------------------------------------- */
/*       Modified - Sparse Row  to   Compressed Sparse Row */

/* ----------------------------------------------------------------------- */
/* converts a compressed matrix using a separated diagonal */
/* (modified sparse row format) in the Compressed Sparse Row */
/* format. */
/* does not check for zero elements in the diagonal. */


/* on entry : */
/* --------- */
/* n          = row dimension of matrix */
/* a, ja      = sparse matrix in msr sparse storage format */
/*              see routine csrmsr for details on data structure */

/* on return : */
/* ----------- */

/* ao,jao,iao = output matrix in csr format. */

/* work arrays: */
/* ------------ */
/* wk       = real work array of length n */
/* iwk      = integer work array of length n+1 */

/* notes: */
/*   The original version of this was NOT in place, but has */
/*   been modified by adding the vector iwk to be in place. */
/*   The original version had ja instead of iwk everywhere in */
/*   loop 500.  Modified  Sun 29 May 1994 by R. Bramley (Indiana). */

/* ----------------------------------------------------------------------- */
    /* Parameter adjustments */
    --iwk;
    --wk;
    --iao;
    --a;
    --ja;
    --ao;
    --jao;

    /* Function Body */
    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	wk[i__] = a[i__];
	iwk[i__] = ja[i__];
/* L1: */
    }
    iwk[n + 1] = ja[n + 1];
    iao[1] = 1;
    iptr = 1;
/* --------- */
    i__1 = n;
    for (ii = 1; ii <= i__1; ++ii) {
	added = false;
	idiag = iptr + (iwk[ii + 1] - iwk[ii]);
	i__2 = iwk[ii + 1] - 1;
	for (k = iwk[ii]; k <= i__2; ++k) {
	    j = ja[k];
	    if (j < ii) {
		ao[iptr] = a[k];
		jao[iptr] = j;
		++iptr;
	    } else if (added) {
		ao[iptr] = a[k];
		jao[iptr] = j;
		++iptr;
	    } else {
/* add diag element - only reserve a position for it. */
		idiag = iptr;
		++iptr;
		added = true;
/*     then other element */
		ao[iptr] = a[k];
		jao[iptr] = j;
		++iptr;
	    }
/* L100: */
	}
	ao[idiag] = wk[ii];
	jao[idiag] = ii;
	if (! added) {
	    ++iptr;
	}
	iao[ii + 1] = iptr;
/* L500: */
    }
    
    ++iwk;
    ++wk;
    ++iao;
    ++a;
    ++ja;
    ++ao;
    ++jao;


/* ------------ end of subroutine msrcsr --------------------------------- */
/* ----------------------------------------------------------------------- */
} /* msrcsr_ */

// Неполное LU разложение. Возвращает код ошибки или признак успешного завершения.
// На входе матрица в CRS формате : a, ja, ia. 
// На выходе матрица ILU0 разложения в MSR формате alu, jlu.
// ju - указатель на диагональные элементы в MSR.
// ierr - возвращаемый код ошибки. 
// iw - рабочий массив длины n;
/* ---------------------------------------------------------------------- */
/* Subroutine */ integer ilu0_(integer n, doublereal* &a, integer* &ja, integer* &ia,
						   doublereal* &alu, integer* &jlu, integer* &ju, integer* &iw, integer &ierr)
{
    

    /* Local variables */
    integer i__, j, jf, ii, jj, jm, js;
    doublereal tl;
    integer jw, ju0, jcol, jrow;

/* ------------------ right preconditioner ------------------------------* */
/*                    ***   ilu(0) preconditioner.   ***                * */
/* ----------------------------------------------------------------------* */
/* Note that this has been coded in such a way that it can be used */
/* with pgmres. Normally, since the data structure of the L+U matrix is */
/* the same as that the A matrix, savings can be made. In fact with */
/* some definitions (not correct for general sparse matrices) all we */
/* need in addition to a, ja, ia is an additional diagonal. */
/* ILU0 is not recommended for serious problems. It is only provided */
/* here for comparison purposes. */
/* ----------------------------------------------------------------------- */

/* on entry: */
/* --------- */
/* n       = dimension of matrix */
/* a, ja, */
/* ia      = original matrix in compressed sparse row storage. */

/* on return: */
/* ----------- */
/* alu,jlu = matrix stored in Modified Sparse Row (MSR) format containing */
/*           the L and U factors together. The diagonal (stored in */
/*           alu(1:n) ) is inverted. Each i-th row of the alu,jlu matrix */
/*           contains the i-th row of L (excluding the diagonal entry=1) */
/*           followed by the i-th row of U. */

/* ju	  = pointer to the diagonal elements in alu, jlu. */

/* ierr	  = integer indicating error code on return */
/* 	     ierr = 0 --> normal return */
/* 	     ierr = k --> code encountered a zero pivot at step k. */
/* work arrays: */
/* ------------- */
/* iw	    = integer work array of length n. */
/* ------------ */
/* IMPORTANT */
/* ----------- */
/* it is assumed that the the elements in the input matrix are stored */
/*    in such a way that in each row the lower part comes first and */
/*    then the upper part. To get the correct ILU factorization, it is */
/*    also necessary to have the elements of L sorted by increasing */
/*    column number. It may therefore be necessary to sort the */
/*    elements of a, ja, ia prior to calling ilu0. This can be */
/*    achieved by transposing the matrix twice using csrcsc. */

/* ----------------------------------------------------------------------- */
    /* Parameter adjustments */
    --iw;
    --ju;
    --jlu;
    --alu;
    --ia;
    --ja;
    --a;

    /* Function Body */
    ju0 = n + 2;
    jlu[1] = ju0;

/* initialize work vector to zero's */


    for (i__ = 1; i__ <= n; ++i__) {
	    iw[i__] = 0;
/* L31: */
    }

/* main loop */

    for (ii = 1; ii <= n; ++ii) 
	{
	   js = ju0;

       /* generating row number ii of L and U. */

       for (j = ia[ii]; j <= ia[ii + 1] - 1; ++j) 
	   {

/*     copy row ii of a, ja, ia into row ii of alu, jlu (L/U) matrix. */

	       jcol = ja[j];
	       if (jcol == ii) {
		      alu[ii] = a[j];
		      iw[jcol] = ii;
		      ju[ii] = ju0;
	       } else {
		      alu[ju0] = a[j];
		      jlu[ju0] = ja[j];
		      iw[jcol] = ju0;
		      ++ju0;
	       }
/* L100: */
	   }
	   jlu[ii + 1] = ju0;
	   jf = ju0 - 1;
	   jm = ju[ii] - 1;

/*     exit if diagonal element is reached. */


	   for (j = js; j <= jm; ++j) {
	       jrow = jlu[j];
	       tl = alu[j] * alu[jrow];
	       alu[j] = tl;

/*     perform  linear combination */

	       for (jj = ju[jrow]; jj <= jlu[jrow + 1] - 1; ++jj) {
		       jw = iw[jlu[jj]];
		       if (jw != 0) {
		           alu[jw] -= tl * alu[jj];
		       }
/* L140: */
	       }
/* L150: */
	    }

/*     invert  and store diagonal element. */

	if (fabs(alu[ii]) < 1.0e-30) {
	    /*     zero pivot : */
        // Нулевой диагональный элемент в строке ii-1
        ierr = ii;
#if doubleintprecision == 1
		printf("Error !!! Zero diagonal element in string %lld \n", ii);
#else
		printf("Error !!! Zero diagonal element in string %d \n", ii);
#endif
        
        printf("Please, press any key to exit calculation...\n");
       // getchar();
		system("pause");
        /* Parameter adjustments */
        ++iw;
        ++ju;
        ++jlu;
        ++alu;
        ++ia;
        ++ja;
        ++a;
        exit(0);
        return 0;
	}
	else {
	  
		alu[ii] = 1.0 / alu[ii];

         /*     reset pointer iw to zero */

	     iw[ii] = 0;
	    
	     for (i__ = js; i__ <= jf; ++i__) {
             /* L201: */
	         iw[jlu[i__]] = 0;
	     }
	}
/* L500: */
    } // закрывающая скобка main цикла.
    
	ierr = 0;

    /* Parameter adjustments */
    ++iw;
    ++ju;
    ++jlu;
    ++alu;
    ++ia;
    ++ja;
    ++a;

    return 0;  

    
/* ------- end-of-ilu0 --------------------------------------------------- */
/* ----------------------------------------------------------------------- */
} /* ilu0_ */

// Решает ILU систему. (LU)x=y; 
// Матрица ILU декомпозиции подаётся в MSR формате :
// alu, jlu, ju; n - это размерность вектора и размер квадратной матрицы СЛАУ.
/* ----------------------------------------------------------------------- */
/* Subroutine */ integer lusol_1(integer n, doublereal* &y, doublereal* &x, 
	doublereal* &alu, integer* &jlu, integer* &ju, integer maxelm)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    integer i__, k;

/* ----------------------------------------------------------------------- */

/* This routine solves the system (LU) x = y, */
/* given an LU decomposition of a matrix stored in (alu, jlu, ju) */
/* modified sparse row format */

/* ----------------------------------------------------------------------- */
/* on entry: */
/* n   = dimension of system */
/* y   = the right-hand-side vector */
/* alu, jlu, ju */
/*     = the LU matrix as provided from the ILU routines. */

/* on return */
/* x   = solution of LU x = y. */
/* ----------------------------------------------------------------------- */

/* Note: routine is in place: call lusol (n, x, x, alu, jlu, ju) */
/*       will solve the system with rhs x and overwrite the result on x . */

/* ----------------------------------------------------------------------- */
/* local variables */


/* forward solve */

    /* Parameter adjustments */
    --x;
    --y;
    --alu;
    --jlu;
    --ju;

	// Литература Г.Б.Сушко, С.А.Харченко.
	// Многопоточная параллельная реализация итерационного алгоритма
	// решения систем линейных уравнений с динамическим распределением 
	// нагрузки по нитям вычислений.
	/* 9 августа 2015.
	LUx=y;
	L^(-1)*LUx=L^(-1)*y;
	z=L^(-1)*y;
	Ux=z;
	x=U^(-1)*z;
	*/

	/*
#ifdef _OPENMP

	// При решении треугольной системы с L имеют место зависимости от листьев дерева к корням.
	// Сначала в произвольном порядке можно вычислить неизвестные, соответствующие листьям верхнего
	// уровня бинарного дерева зависимостей. Затем можно вычислять неизвестные любого узла дерева 
	// зависимостей как только закончены вычисления с двумя узлами, от которых зависит этот узел.

			if (inumcore==2) {
				if (nd.b0.active) {

					// первый поток
					for (integer iscan_par=nd.b0.ileft_start; iscan_par<=nd.b0.ileft_finish; iscan_par++) {
						integer iPloc=ifrontregulationgl[iscan_par];

						// Сначала обработке подвергнутся лишь граничные КО.
						// А уже во вторую очередь внутренние КО.
						if (iPloc>=maxelm) {
						   // i__1 = n;
                           //for (i__ = 1; i__ <= i__1; ++i__) {
						   i__=iPloc+1;
	                       x[i__] = y[i__];
	                       i__2 = ju[i__] - 1;
	                       for (k = jlu[i__]; k <= i__2; ++k) {
	                           x[i__] -= alu[k] * x[jlu[k]];
                               
	                       }
                           //}
						}

					}
					for (integer iscan_par=nd.b0.ileft_start; iscan_par<=nd.b0.ileft_finish; iscan_par++) {
						integer iPloc=ifrontregulationgl[iscan_par];

						// Сначала обработке подвергнутся лишь граничные КО.
						// А уже во вторую очередь внутренние КО.
						if (iPloc<maxelm) {
						   // i__1 = n;
                           //for (i__ = 1; i__ <= i__1; ++i__) {
						   i__=iPloc+1;
	                       x[i__] = y[i__];
	                       i__2 = ju[i__] - 1;
	                       for (k = jlu[i__]; k <= i__2; ++k) {
	                           x[i__] -= alu[k] * x[jlu[k]];
                               
	                       }
                          
                           //}
						}

					}
					// второй поток
					for (integer iscan_par=nd.b0.iright_start; iscan_par<=nd.b0.iright_finish; iscan_par++) {
						integer iPloc=ifrontregulationgl[iscan_par];
						// Сначала обработке подвергнутся лишь граничные КО.
						// А уже во вторую очередь внутренние КО.
						if (iPloc>=maxelm) {
						// i__1 = n;
                           //for (i__ = 1; i__ <= i__1; ++i__) {
						   i__=iPloc+1;
	                       x[i__] = y[i__];
	                       i__2 = ju[i__] - 1;
	                       for (k = jlu[i__]; k <= i__2; ++k) {
	                           x[i__] -= alu[k] * x[jlu[k]];
                              
	                       }
                           
                           //}
						}
					}
					for (integer iscan_par=nd.b0.iright_start; iscan_par<=nd.b0.iright_finish; iscan_par++) {
						integer iPloc=ifrontregulationgl[iscan_par];
						// Сначала обработке подвергнутся лишь граничные КО.
						// А уже во вторую очередь внутренние КО.
						if (iPloc<maxelm) {
						// i__1 = n;
                           //for (i__ = 1; i__ <= i__1; ++i__) {
						   i__=iPloc+1;
	                       x[i__] = y[i__];
	                       i__2 = ju[i__] - 1;
	                       for (k = jlu[i__]; k <= i__2; ++k) {
	                           x[i__] -= alu[k] * x[jlu[k]];
                             
	                       }
                          
                           //}
						}
					}
					// серийный смыкающий кусок
					for (integer iscan_par=nd.b0.iseparate_start; iscan_par<=nd.b0.iseparate_finish; iscan_par++) {
						integer iPloc=ifrontregulationgl[iscan_par];
						// Сначала обработке подвергнутся лишь граничные КО.
						// А уже во вторую очередь внутренние КО.
						if (iPloc>=maxelm) {
						// i__1 = n;
                           //for (i__ = 1; i__ <= i__1; ++i__) {
						   i__=iPloc+1;
	                       x[i__] = y[i__];
	                       i__2 = ju[i__] - 1;
	                       for (k = jlu[i__]; k <= i__2; ++k) {
	                           x[i__] -= alu[k] * x[jlu[k]];
                              
	                       }
                          
                           //}
						}
					}
					for (integer iscan_par=nd.b0.iseparate_start; iscan_par<=nd.b0.iseparate_finish; iscan_par++) {
						integer iPloc=ifrontregulationgl[iscan_par];
						// Сначала обработке подвергнутся лишь граничные КО.
						// А уже во вторую очередь внутренние КО.
						if (iPloc<maxelm) {
						// i__1 = n;
                           //for (i__ = 1; i__ <= i__1; ++i__) {
						   i__=iPloc+1;
	                       x[i__] = y[i__];
	                       i__2 = ju[i__] - 1;
	                       for (k = jlu[i__]; k <= i__2; ++k) {
	                           x[i__] -= alu[k] * x[jlu[k]];
                              
	                       }
                           
                           //}
						}
					}


				}
			}

#else*/

    /* Function Body */
	// L : z=L^(-1)*y;
    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	    x[i__] = y[i__];
	    i__2 = ju[i__] - 1;
		//doublereal dsum = 0.0;
//#pragma loop(hint_parallel(8))
//#pragma omp parallel for private(k) reduction(+:dsum)
	    for (k = jlu[i__]; k <= i__2; ++k) {
	        x[i__] -= alu[k] * x[jlu[k]];
			//dsum -= alu[k] * x[jlu[k]];
            /* L41: */
	    }
		//x[i__] += dsum;
        /* L40: */
    }

//#endif

	/*
#ifdef _OPENMP

			if (inumcore==2) {
				if (nd.b0.active) {

					// серийный смыкающий кусок
					for (integer iscan_par=nd.b0.iseparate_start; iscan_par<=nd.b0.iseparate_finish; iscan_par++) {
						integer iPloc=ifrontregulationgl[iscan_par];
						// Сначала обработке подвергнутся лишь граничные КО.
						// А уже во вторую очередь внутренние КО.
						if (iPloc>=maxelm) {
							i__=iPloc+1;

						    i__1 = jlu[i__ + 1] - 1;
	                        for (k = ju[i__]; k <= i__1; ++k) {
	                            x[i__] -= alu[k] * x[jlu[k]];
                                
	                        }
	                        x[i__] = alu[i__] * x[i__];
						}
					}
					for (integer iscan_par=nd.b0.iseparate_start; iscan_par<=nd.b0.iseparate_finish; iscan_par++) {
						integer iPloc=ifrontregulationgl[iscan_par];
						// Сначала обработке подвергнутся лишь граничные КО.
						// А уже во вторую очередь внутренние КО.
						if (iPloc<maxelm) {
							i__=iPloc+1;

						    i__1 = jlu[i__ + 1] - 1;
	                        for (k = ju[i__]; k <= i__1; ++k) {
	                        x[i__] -= alu[k] * x[jlu[k]];
                           
	                    }
	                    x[i__] = alu[i__] * x[i__];
						}
					}
					// первый поток
					for (integer iscan_par=nd.b0.ileft_start; iscan_par<=nd.b0.ileft_finish; iscan_par++) {
						integer iPloc=ifrontregulationgl[iscan_par];
						// Сначала обработке подвергнутся лишь граничные КО.
						// А уже во вторую очередь внутренние КО.
						if (iPloc>=maxelm) {
							i__=iPloc+1;

						    i__1 = jlu[i__ + 1] - 1;
	                        for (k = ju[i__]; k <= i__1; ++k) {
	                            x[i__] -= alu[k] * x[jlu[k]];
                               
	                        }
	                        x[i__] = alu[i__] * x[i__];
						}
					}
					for (integer iscan_par=nd.b0.ileft_start; iscan_par<=nd.b0.ileft_finish; iscan_par++) {
						integer iPloc=ifrontregulationgl[iscan_par];
						// Сначала обработке подвергнутся лишь граничные КО.
						// А уже во вторую очередь внутренние КО.
						if (iPloc<maxelm) {
							i__=iPloc+1;

						    i__1 = jlu[i__ + 1] - 1;
	                        for (k = ju[i__]; k <= i__1; ++k) {
	                        x[i__] -= alu[k] * x[jlu[k]];
                           
	                    }
	                    x[i__] = alu[i__] * x[i__];
						}
					}
					// второй поток
					for (integer iscan_par=nd.b0.iright_start; iscan_par<=nd.b0.iright_finish; iscan_par++) {
						integer iPloc=ifrontregulationgl[iscan_par];
						// Сначала обработке подвергнутся лишь граничные КО.
						// А уже во вторую очередь внутренние КО.
						if (iPloc>=maxelm) {
							i__=iPloc+1;

						    i__1 = jlu[i__ + 1] - 1;
	                        for (k = ju[i__]; k <= i__1; ++k) {
	                            x[i__] -= alu[k] * x[jlu[k]];
                               
	                        }
	                        x[i__] = alu[i__] * x[i__];
						}
					}
					for (integer iscan_par=nd.b0.iright_start; iscan_par<=nd.b0.iright_finish; iscan_par++) {
						integer iPloc=ifrontregulationgl[iscan_par];
						// Сначала обработке подвергнутся лишь граничные КО.
						// А уже во вторую очередь внутренние КО.
						if (iPloc<maxelm) {
							i__=iPloc+1;

						    i__1 = jlu[i__ + 1] - 1;
	                        for (k = ju[i__]; k <= i__1; ++k) {
	                        x[i__] -= alu[k] * x[jlu[k]];
                           
	                    }
	                    x[i__] = alu[i__] * x[i__];
						}
					}


				}
			}

#else*/

/*     backward solve. */
	// U : x=U^(-1)*z;
    for (i__ = n; i__ >= 1; --i__) {
	     i__1 = jlu[i__ + 1] - 1;
//#pragma loop(hint_parallel(8))
	     for (k = ju[i__]; k <= i__1; ++k) {
	          x[i__] -= alu[k] * x[jlu[k]];
              /* L91: */
	     }
	     x[i__] = alu[i__] * x[i__];
         /* L90: */
    }

//#endif

    /* Parameter adjustments */
    ++x;
    ++y;
    ++alu;
    ++jlu;
    ++ju;


    return 0;
/* ----------------end of lusol ------------------------------------------ */
/* ----------------------------------------------------------------------- */
} /* lusol_ */

// 9 ноября 2016 года.
// Экономия оперативной памяти для ILU.
// Централизованное хранилище.
typedef struct TLEVEL_ADDITIONAL_DATA_BUFER {
	// lusol_ портит матрицы и вектор правой части
	// поэтому будем работать только на копиях объектов
	doublereal* alu_copy;
	integer* jlu_copy;
	integer* ju_copy;
} LEVEL_ADDITIONAL_DATA_BUFER;

// Централизованное хранилище вспомогательных даных,
// эта вещь используется на всех уровнях как временное хранилище.
// Раньше хранение производилось в каждом уровне индивидуально, 
// что было неприемлемо по памяти на больших задачах.
LEVEL_ADDITIONAL_DATA_BUFER milu_gl_buffer;

// 10 августа 2016 года.
// Дополнительная информация на каждом из уровней,
// она может использоваться и для использования 
// сглаживателей на основе методов подпространства Крылова.
typedef struct TLEVEL_ADDITIONAL_DATA {
	// для хранения матрицы на данном уровне.
	doublereal* val;
	integer* col_ind;
	integer* row_ptr;
	// для хранения ilu2 декомпозиции.
	static const integer lfil = 2;
	integer maxelm_plus_maxbound;
	integer iwk;
	doublereal* alu;
	integer* jlu;
	integer* ju;
	integer* levs;
	doublereal* w;
	integer* jw;
	// lusol_ портит матрицы и вектор правой части
	// поэтому будем работать только на копиях объектов
	//doublereal* alu_copy;
	//integer* jlu_copy;
	//integer* ju_copy;
	doublereal* b_copy;
	doublereal* x_copy;
	doublereal* zbuf;
	doublereal* zbuf2;
} LEVEL_ADDITIONAL_DATA;

// 4 ноября 2016 года.
// Дополнительная информация на каждом из уровней,
// она может использоваться и для использования 
// сглаживателей на основе методов подпространства Крылова.
typedef struct TLEVEL_ADDITIONAL_DATA0 {
	// для хранения матрицы на данном уровне.
	doublereal* val;
	integer* col_ind;
	integer* row_ptr;
	// для хранения ilu0 декомпозиции.

	// Матрица разложения.
	doublereal* alu;
	integer* jlu;
	integer* ju;
	// Вспомогательный вектор.
	integer* iw;
	
	// for lusol_1patchforRUMBA
	integer iwk;
	
	// lusol_ портит матрицы и вектор правой части
	// поэтому будем работать только на копиях объектов
	doublereal* alu_copy;
	integer* jlu_copy;
	integer* ju_copy;
	doublereal* b_copy;
	doublereal* x_copy;
	doublereal* zbuf;
	doublereal* zbuf2;
} LEVEL_ADDITIONAL_DATA0;

// 10.08.2016
// Патч для совместимости с РУМБА Алгоритмом.
// 9 ноября 2016 переписано для централизованного хранения памяти в буффере. 
integer lusol_1patchforRUMBA(integer n, doublereal* &y, doublereal* &x,
	LEVEL_ADDITIONAL_DATA &milu2)
{
	// x - argument
	// y - rthdsd

	integer iret = 0;

	
	// move
	for (integer i_1 = 0; i_1 < n; i_1++) {
		milu2.x_copy[i_1] = x[i_1 + 1];
		milu2.b_copy[i_1] = y[i_1 + 1];
	}

	for (integer i7 = 0; i7<milu2.iwk + 2; i7++) {
		//milu2.alu_copy[i7] = milu2.alu[i7];
		//milu2.jlu_copy[i7] = milu2.jlu[i7];
		// Записываем для использования.
		// При этом диспетчер памяти сам контролирует что памяти должно хватить.
		milu_gl_buffer.alu_copy[i7] = milu2.alu[i7];
		milu_gl_buffer.jlu_copy[i7] = milu2.jlu[i7];
	}
	//for (integer i7 = 0; i7<n + 2; i7++) milu2.ju_copy[i7] = milu2.ju[i7];
	for (integer i7 = 0; i7<n + 2; i7++) milu_gl_buffer.ju_copy[i7] = milu2.ju[i7];

	//iret = lusol_1(n, milu2.b_copy, milu2.x_copy,
		//milu2.alu_copy, milu2.jlu_copy, milu2.ju_copy, 0);

	iret = lusol_1(n, milu2.b_copy, milu2.x_copy,
		milu_gl_buffer.alu_copy, milu_gl_buffer.jlu_copy, milu_gl_buffer.ju_copy, 0);

	// move
	for (integer i_1 = 0; i_1 < n; i_1++) {
		x[i_1 + 1] = milu2.x_copy[i_1];
		milu2.b_copy[i_1] = y[i_1 + 1];
	}

	// неразрушающее восстановление.
	for (integer i7 = 0; i7<milu2.iwk + 2; i7++) {
		//milu2.alu_copy[i7] = milu2.alu[i7];
		//milu2.jlu_copy[i7] = milu2.jlu[i7];
		milu_gl_buffer.alu_copy[i7] = milu2.alu[i7];
		milu_gl_buffer.jlu_copy[i7] = milu2.jlu[i7];
	}
	//for (integer i7 = 0; i7<n + 2; i7++) milu2.ju_copy[i7] = milu2.ju[i7];
	for (integer i7 = 0; i7<n + 2; i7++) milu_gl_buffer.ju_copy[i7] = milu2.ju[i7];

	return iret;

} // lusol_1patchforRUMBA

  // 10.08.2016 04.november.2016
  // Патч для совместимости с РУМБА Алгоритмом.
integer lusol_1patchforRUMBA(integer n, doublereal* &y, doublereal* &x,
	LEVEL_ADDITIONAL_DATA0 &milu0)
{
	// x - argument
	// y - rthdsd

	integer iret = 0;


	// move
	for (integer i_1 = 0; i_1 < n; i_1++) {
		milu0.x_copy[i_1] = x[i_1 + 1];
		milu0.b_copy[i_1] = y[i_1 + 1];
	}

	for (integer i7 = 0; i7<milu0.iwk + 2; i7++) {
		milu0.alu_copy[i7] = milu0.alu[i7];
		milu0.jlu_copy[i7] = milu0.jlu[i7];
	}
	for (integer i7 = 0; i7<n + 2; i7++) milu0.ju_copy[i7] = milu0.ju[i7];


	iret = lusol_1(n, milu0.b_copy, milu0.x_copy,
		milu0.alu_copy, milu0.jlu_copy, milu0.ju_copy, 0);

	// move
	for (integer i_1 = 0; i_1 < n; i_1++) {
		x[i_1 + 1] = milu0.x_copy[i_1];
		milu0.b_copy[i_1] = y[i_1 + 1];
	}

	// неразрушающее восстановление.
	for (integer i7 = 0; i7<milu0.iwk + 2; i7++) {
		milu0.alu_copy[i7] = milu0.alu[i7];
		milu0.jlu_copy[i7] = milu0.jlu[i7];
	}
	for (integer i7 = 0; i7<n + 2; i7++) milu0.ju_copy[i7] = milu0.ju[i7];

	return iret;

} // lusol_1patchforRUMBA

// Решает ILU систему. (LU)x=y; 
// Матрица ILU декомпозиции подаётся в MSR формате :
// alu, jlu, ju; n - это размерность вектора и размер квадратной матрицы СЛАУ.
/* ----------------------------------------------------------------------- */
/* Subroutine */ integer lusol_2(integer n, doublereal* &y, doublereal* &x, 
	doublereal* &alu, integer* &jlu, integer* &ju,  doublereal* &x_copy, integer maxelm)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    integer i__, k;

/* ----------------------------------------------------------------------- */

/* This routine solves the system (LU) x = y, */
/* given an LU decomposition of a matrix stored in (alu, jlu, ju) */
/* modified sparse row format */

/* ----------------------------------------------------------------------- */
/* on entry: */
/* n   = dimension of system */
/* y   = the right-hand-side vector */
/* alu, jlu, ju */
/*     = the LU matrix as provided from the ILU routines. */

/* on return */
/* x   = solution of LU x = y. */
/* ----------------------------------------------------------------------- */

/* Note: routine is in place: call lusol (n, x, x, alu, jlu, ju) */
/*       will solve the system with rhs x and overwrite the result on x . */

/* ----------------------------------------------------------------------- */
/* local variables */


/* forward solve */

    /* Parameter adjustments */
    --x;
    --y;
    --alu;
    --jlu;
    --ju;

	// Литература Г.Б.Сушко, С.А.Харченко.
	// Многопоточная параллельная реализация итерационного алгоритма
	// решения систем линейных уравнений с динамическим распределением 
	// нагрузки по нитям вычислений.
	/* 9 августа 2015.
	LUx=y;
	L^(-1)*LUx=L^(-1)*y;
	z=L^(-1)*y;
	Ux=z;
	x=U^(-1)*z;
	*/

	
#ifdef _OPENMP

	// При решении треугольной системы с L имеют место зависимости от листьев дерева к корням.
	// Сначала в произвольном порядке можно вычислить неизвестные, соответствующие листьям верхнего
	// уровня бинарного дерева зависимостей. Затем можно вычислять неизвестные любого узла дерева 
	// зависимостей как только закончены вычисления с двумя узлами, от которых зависит этот узел.


			if (inumcore==2) {
				if (nd.b0.active) {

#pragma omp parallel private(i__, i__2,k)
	{
#pragma omp sections
		{
#pragma omp section
			{

					// первый поток
					for (integer iscan_par=nd.b0.ileft_start; iscan_par<=nd.b0.ileft_finish; iscan_par++) {
						
						integer iPloc=iscan_par;

						
						   i__=iPloc+1;
	                       x[i__] = y[i__];
	                       i__2 = ju[i__] - 1;
	                       for (k = jlu[i__]; k <= i__2; ++k) {
	                           x[i__] -= alu[k] * x[jlu[k]];
                               
	                       }
                      

					}
					
			}
#pragma omp section
			{

					// второй поток
					for (integer iscan_par=nd.b0.iright_start; iscan_par<=nd.b0.iright_finish; iscan_par++) {
						
						integer iPloc=iscan_par;

						
						   i__=iPloc+1;
	                       x[i__] = y[i__];
	                       i__2 = ju[i__] - 1;
	                       for (k = jlu[i__]; k <= i__2; ++k) {
	                           x[i__] -= alu[k] * x[jlu[k]];
                              
	                       }
                           
                         
					}
			} // pragma omp section
		} // pragma omp sections
	} // pragma omp parallel
					
					// серийный смыкающий кусок
					for (integer iscan_par=nd.b0.iseparate_start; iscan_par<=nd.b0.iseparate_finish; iscan_par++) {
						
						integer iPloc=iscan_par;
						
						
						   i__=iPloc+1;
	                       x[i__] = y[i__];
	                       i__2 = ju[i__] - 1;
	                       for (k = jlu[i__]; k <= i__2; ++k) {
	                           x[i__] -= alu[k] * x[jlu[k]];
                              
	                       }
                          
                         
					}
					


				}
			}

#else

    /* Function Body */
	// L : z=L^(-1)*y;
    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	    x[i__] = y[i__];
	    i__2 = ju[i__] - 1;
	    for (k = jlu[i__]; k <= i__2; ++k) {
	        x[i__] -= alu[k] * x[jlu[k]];
            /* L41: */
	    }
        /* L40: */
    }

#endif

	
#ifdef _OPENMP

			if (inumcore==2) {
				if (nd.b0.active) {

					// серийный смыкающий кусок
					for (integer iscan_par=nd.b0.iseparate_finish; iscan_par>=nd.b0.iseparate_start; iscan_par--) {
						
						integer iPloc=iscan_par;
						
							i__=iPloc+1;

						    i__1 = jlu[i__ + 1] - 1;
	                        for (k = ju[i__]; k <= i__1; ++k) {
	                            x[i__] -= alu[k] * x[jlu[k]];
                                
	                        }
	                        x[i__] = alu[i__] * x[i__];
						
					}
					//doublereal *x_copy=new doublereal[n];
					for (integer icop=1; icop<=n; icop++) x_copy[icop-1]=x[icop];
					
#pragma omp parallel private(i__, i__1,k)
					{
#pragma omp sections
						{
#pragma omp section
							{

					// первый поток
					for (integer iscan_par=nd.b0.ileft_finish; iscan_par>=nd.b0.ileft_start; iscan_par--) {
						
						integer iPloc=iscan_par;
						
						
							i__=iPloc+1;

						    i__1 = jlu[i__ + 1] - 1;
	                        for (k = ju[i__]; k <= i__1; ++k) {
	                            x_copy[i__-1] -= alu[k] * x_copy[jlu[k]-1];
                               
	                        }
	                        x_copy[i__-1] = alu[i__] * x_copy[i__-1];
						
					}
							}
#pragma omp section
							{
					
					// второй поток
					for (integer iscan_par=nd.b0.iright_finish; iscan_par>=nd.b0.iright_start; iscan_par--) {
						
						integer iPloc=iscan_par;
						
						
							i__=iPloc+1;

						    i__1 = jlu[i__ + 1] - 1;
	                        for (k = ju[i__]; k <= i__1; ++k) {
	                            x[i__] -= alu[k] * x[jlu[k]];
                               
	                        }
	                        x[i__] = alu[i__] * x[i__];
						
					}
							} // omp section
						} // omp sections
					} // omp parallel
				
					for (integer iscan_par=nd.b0.ileft_finish; iscan_par>=nd.b0.ileft_start; iscan_par--) {
						x[iscan_par+1]=x_copy[iscan_par];
					}
                   // delete x_copy;

				}
			}

#else

/*     backward solve. */
	// U : x=U^(-1)*z;
    for (i__ = n; i__ >= 1; --i__) {
	     i__1 = jlu[i__ + 1] - 1;
	     for (k = ju[i__]; k <= i__1; ++k) {
	          x[i__] -= alu[k] * x[jlu[k]];
              /* L91: */
	     }
	     x[i__] = alu[i__] * x[i__];
         /* L90: */
    }

#endif

    /* Parameter adjustments */
    ++x;
    ++y;
    ++alu;
    ++jlu;
    ++ju;


    return 0;
/* ----------------end of lusol ------------------------------------------ */
/* ----------------------------------------------------------------------- */
} /* lusol_2 */


// Решает ILU систему. (LU)x=y; 
// Матрица ILU декомпозиции подаётся в MSR формате :
// alu, jlu, ju; n - это размерность вектора и размер квадратной матрицы СЛАУ.
/* ----------------------------------------------------------------------- */
/* Subroutine */ integer lusol_3(integer n, doublereal* &y, doublereal* &x, 
	doublereal* &alu, integer* &jlu, integer* &ju, integer maxelm)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    integer i__, k;

/* ----------------------------------------------------------------------- */

/* This routine solves the system (LU) x = y, */
/* given an LU decomposition of a matrix stored in (alu, jlu, ju) */
/* modified sparse row format */

/* ----------------------------------------------------------------------- */
/* on entry: */
/* n   = dimension of system */
/* y   = the right-hand-side vector */
/* alu, jlu, ju */
/*     = the LU matrix as provided from the ILU routines. */

/* on return */
/* x   = solution of LU x = y. */
/* ----------------------------------------------------------------------- */

/* Note: routine is in place: call lusol (n, x, x, alu, jlu, ju) */
/*       will solve the system with rhs x and overwrite the result on x . */

/* ----------------------------------------------------------------------- */
/* local variables */


/* forward solve */

    /* Parameter adjustments */
    --x;
    --y;
    --alu;
    --jlu;
    --ju;

	// Литература Г.Б.Сушко, С.А.Харченко.
	// Многопоточная параллельная реализация итерационного алгоритма
	// решения систем линейных уравнений с динамическим распределением 
	// нагрузки по нитям вычислений.
	/* 9 августа 2015.
	LUx=y;
	L^(-1)*LUx=L^(-1)*y;
	z=L^(-1)*y;
	Ux=z;
	x=U^(-1)*z;
	*/

	
#ifdef _OPENMP

	// При решении треугольной системы с L имеют место зависимости от листьев дерева к корням.
	// Сначала в произвольном порядке можно вычислить неизвестные, соответствующие листьям верхнего
	// уровня бинарного дерева зависимостей. Затем можно вычислять неизвестные любого узла дерева 
	// зависимостей как только закончены вычисления с двумя узлами, от которых зависит этот узел.


			if (inumcore==2) {
				if (nd.b0.active) {

#pragma omp parallel private(i__, i__2,k)
	{
#pragma omp sections
		{
#pragma omp section
			{

					// первый поток
					for (integer iscan_par=nd.b0.ileft_start; iscan_par<=nd.b0.ileft_finish; iscan_par++) {
						
						integer iPloc=iscan_par;

						
						   i__=iPloc+1;
	                       x[i__] = y[i__];
	                       i__2 = ju[i__] - 1;
	                       for (k = jlu[i__]; k <= i__2; ++k) {
	                           x[i__] -= alu[k] * x[jlu[k]];
                               
	                       }
                      

					}
					
			}
#pragma omp section
			{

					// второй поток
					for (integer iscan_par=nd.b0.iright_start; iscan_par<=nd.b0.iright_finish; iscan_par++) {
						
						integer iPloc=iscan_par;

						
						   i__=iPloc+1;
	                       x[i__] = y[i__];
	                       i__2 = ju[i__] - 1;
	                       for (k = jlu[i__]; k <= i__2; ++k) {
	                           x[i__] -= alu[k] * x[jlu[k]];
                              
	                       }
                           
                         
					}
			} // pragma omp section
		} // pragma omp sections
	} // pragma omp parallel
					
					// серийный смыкающий кусок
					for (integer iscan_par=nd.b0.iseparate_start; iscan_par<=nd.b0.iseparate_finish; iscan_par++) {
						
						integer iPloc=iscan_par;
						
						
						   i__=iPloc+1;
	                       x[i__] = y[i__];
	                       i__2 = ju[i__] - 1;
	                       for (k = jlu[i__]; k <= i__2; ++k) {
	                           x[i__] -= alu[k] * x[jlu[k]];
                              
	                       }
                          
                         
					}
					


				}
			}

#else

    /* Function Body */
	// L : z=L^(-1)*y;
    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	    x[i__] = y[i__];
	    i__2 = ju[i__] - 1;
	    for (k = jlu[i__]; k <= i__2; ++k) {
	        x[i__] -= alu[k] * x[jlu[k]];
            /* L41: */
	    }
        /* L40: */
    }

#endif

	
#ifdef _OPENMP

			if (inumcore==2) {
				if (nd.b0.active) {

					// серийный смыкающий кусок
					for (integer iscan_par=nd.b0.iseparate_finish; iscan_par>=nd.b0.iseparate_start; iscan_par--) {
						
						integer iPloc=iscan_par;
						
							i__=iPloc+1;

						    i__1 = jlu[i__ + 1] - 1;
	                        for (k = ju[i__]; k <= i__1; ++k) {
	                            x[i__] -= alu[k] * x[jlu[k]];
                                
	                        }
	                        x[i__] = alu[i__] * x[i__];
						
					}
					
					
#pragma omp parallel private(i__, i__1,k)
					{
#pragma omp sections
						{
#pragma omp section
							{

					// первый поток
					for (integer iscan_par=nd.b0.ileft_finish; iscan_par>=nd.b0.ileft_start; iscan_par--) {
						
						integer iPloc=iscan_par;
						
						
							i__=iPloc+1;

						    i__1 = jlu[i__ + 1] - 1;
	                        for (k = ju[i__]; k <= i__1; ++k) {
	                            x[i__] -= alu[k] * x[jlu[k]];
                               
	                        }
	                        x[i__] = alu[i__] * x[i__];
						
					}
							}
#pragma omp section
							{
					
					// второй поток
					for (integer iscan_par=nd.b0.iright_finish; iscan_par>=nd.b0.iright_start; iscan_par--) {
						
						integer iPloc=iscan_par;
						
						
							i__=iPloc+1;

						    i__1 = jlu[i__ + 1] - 1;
	                        for (k = ju[i__]; k <= i__1; ++k) {
	                            x[i__] -= alu[k] * x[jlu[k]];
                               
	                        }
	                        x[i__] = alu[i__] * x[i__];
						
					}
							} // omp section
						} // omp sections
					} // omp parallel
				
					

				}
			}

#else

/*     backward solve. */
	// U : x=U^(-1)*z;
    for (i__ = n; i__ >= 1; --i__) {
	     i__1 = jlu[i__ + 1] - 1;
	     for (k = ju[i__]; k <= i__1; ++k) {
	          x[i__] -= alu[k] * x[jlu[k]];
              /* L91: */
	     }
	     x[i__] = alu[i__] * x[i__];
         /* L90: */
    }

#endif

    /* Parameter adjustments */
    ++x;
    ++y;
    ++alu;
    ++jlu;
    ++ju;


    return 0;
/* ----------------end of lusol ------------------------------------------ */
/* ----------------------------------------------------------------------- */
} /* lusol_3 */

integer lusol_(integer n, doublereal* &y, doublereal* &x, 
	doublereal* &alu, integer* &jlu, integer* &ju, integer maxelm)
{
	// Serial
	if (1) {
		// первоначальная версия из  Sparskit2. 
	    lusol_1(n, y, x, alu, jlu, ju, maxelm);
	}
	else {
		// распараллеленная версия 9 августа 2015 года.
		//lusol_2(n, y, x, alu, jlu, ju, maxelm);

		printf("lusol_ -> lusol_2 ustarevshii vjzov. 10 08 2015\n");
	//	getchar();
		system("pause");
		exit(1);
	}

	return 0;
}