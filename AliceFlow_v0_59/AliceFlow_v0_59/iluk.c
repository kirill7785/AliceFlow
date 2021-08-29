/* iluk.f -- translated by f2c (version 20100827).
   SPARSKIT2 Yosef Saad
*/

#include "iluk_quick.cpp"


// Сдержит медленный линейный поиск и из-за этого непригодна.
/* ----------------------------------------------------------------------- */
/* Subroutine */ integer iluk_Saad(integer n, doublereal* &a, integer* &ja, integer* &ia,
						   integer lfil, doublereal* &alu, integer* &jlu, integer* &ju, 
	integer* &levs, integer iwk, doublereal* &w, integer* &jw, integer &ierr)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4;

    /* Local variables */
    integer i__, j, k;
    doublereal s, t;
    integer j1, j2, n2, ii, jj, ju0;
    doublereal fact;
    integer lenl, jlev, lenu, jpos, jrow;

/* ----------------------------------------------------------------------* */
/*     SPARSKIT ROUTINE ILUK -- ILU WITH LEVEL OF FILL-IN OF K (ILU(k)) * */
/* ----------------------------------------------------------------------* */

/* on entry: */
/* ========== */
/* n       = integer. The row dimension of the matrix A. The matrix */

/* a,ja,ia = matrix stored in Compressed Sparse Row format. */

/* lfil    = integer. The fill-in parameter. Each element whose */
/*           leve-of-fill exceeds lfil during the ILU process is dropped. */
/*           lfil must be .ge. 0 */

/* tol     = real*8. Sets the threshold for dropping small terms in the */
/*           factorization. See below for details on dropping strategy. */

/* iwk     = integer. The minimum length of arrays alu, jlu, and levs. */

/* On return: */
/* =========== */

/* alu,jlu = matrix stored in Modified Sparse Row (MSR) format containing */
/*           the L and U factors together. The diagonal (stored in */
/*           alu(1:n) ) is inverted. Each i-th row of the alu,jlu matrix */
/*           contains the i-th row of L (excluding the diagonal entry=1) */
/*           followed by the i-th row of U. */

/* ju      = integer array of length n containing the pointers to */
/*           the beginning of each row of U in the matrix alu,jlu. */

/* levs    = integer (work) array of size iwk -- which contains the */
/*           levels of each element in alu, jlu. */

/* ierr    = integer. Error message with the following meaning. */
/*           ierr  = 0    --> successful return. */
/*           ierr .gt. 0  --> zero pivot encountered at step number ierr. */
/*           ierr  = -1   --> Error. input matrix may be wrong. */
/*                            (The elimination process has generated a */
/*                            row in L or U whose length is .gt.  n.) */
/*           ierr  = -2   --> The matrix L overflows the array al. */
/*           ierr  = -3   --> The matrix U overflows the array alu. */
/*           ierr  = -4   --> Illegal value for lfil. */
/*           ierr  = -5   --> zero row encountered in A or U. */

/* work arrays: */
/* ============= */
/* jw      = integer work array of length 3*n. */
/* w       = real work array of length n */

/* Notes/known bugs: This is not implemented efficiently storage-wise. */
/*       For example: Only the part of the array levs(*) associated with */
/*       the U-matrix is needed in the routine.. So some storage can */
/*       be saved if needed. The levels of fills in the LU matrix are */
/*       output for information only -- they are not needed by LU-solve. */

/* ---------------------------------------------------------------------- */
/* w, ju (1:n) store the working array [1:ii-1 = L-part, ii:n = u] */
/* jw(n+1:2n)  stores the nonzero indicator. */

/* Notes: */
/* ------ */
/* All the diagonal elements of the input matrix must be  nonzero. */

/* ----------------------------------------------------------------------* */
/*     locals */
    /* Parameter adjustments */
    --jw;
    --w;
    --ju;
    --ia;
    --a;
    --ja;
    --alu;
    --jlu;
    --levs;

    /* Function Body */
    if (lfil < 0) {
	goto L998;
    }
/* ----------------------------------------------------------------------- */
/*     initialize ju0 (points to next element to be added to alu,jlu) */
/*     and pointer array. */
/* ----------------------------------------------------------------------- */
    n2 = n + n;
    ju0 = n + 2;
    jlu[1] = ju0;

/*     initialize nonzero indicator array + levs array -- */

    i__1 = n << 1;
    for (j = 1; j <= i__1; ++j) {
	jw[j] = 0;
/* L1: */
    }
/* ----------------------------------------------------------------------- */
/*     beginning of main loop. */
/* ----------------------------------------------------------------------- */
    i__1 = n;
    for (ii = 1; ii <= i__1; ++ii) {
	j1 = ia[ii];
	j2 = ia[ii + 1] - 1;

/*     unpack L-part and U-part of row of A in arrays w */

	lenu = 1;
	lenl = 0;
	jw[ii] = ii;
	w[ii] = 0.f;
	jw[n + ii] = ii;

	i__2 = j2;
	for (j = j1; j <= i__2; ++j) {
	    k = ja[j];
	    t = a[j];
	    if (t == 0.f) {
		goto L170;
	    }
	    if (k < ii) {
		++lenl;
		jw[lenl] = k;
		w[lenl] = t;
		jw[n2 + lenl] = 0;
		jw[n + k] = lenl;
	    } else if (k == ii) {
		w[ii] = t;
		jw[n2 + ii] = 0;
	    } else {
		++lenu;
		jpos = ii + lenu - 1;
		jw[jpos] = k;
		w[jpos] = t;
		jw[n2 + jpos] = 0;
		jw[n + k] = jpos;
	    }
L170:
	    ;
	}

	jj = 0;

/*     eliminate previous rows */

L150:
	++jj;
	if (jj > lenl) {
	    goto L160;
	}
/* ----------------------------------------------------------------------- */
/*     in order to do the elimination in the correct order we must select */
/*     the smallest column index among jw(k), k=jj+1, ..., lenl. */
/* ----------------------------------------------------------------------- */
	jrow = jw[jj];
	k = jj;

/*     determine smallest column index */

	i__2 = lenl;
	// Дьявольски медленный поиск минимума. Это линейный поиск.
	//printf("jj=%d\n",jj);// jj==1 далеко не всегда.
	// Это означает что нужно поддерживать удаление элемента по ключу.
	//getchar();
	for (j = jj + 1; j <= i__2; ++j) {
	    if (jw[j] < jrow) {
		jrow = jw[j];
		k = j;
	    }
/* L151: */
	}

	if (k != jj) {
/*     exchange in jw */
	    j = jw[jj];
	    jw[jj] = jw[k];
	    jw[k] = j;
/*     exchange in jw(n+  (pointers/ nonzero indicator). */
	    jw[n + jrow] = jj;
	    jw[n + j] = k;
/*     exchange in jw(n2+  (levels) */
	    j = jw[n2 + jj];
	    jw[n2 + jj] = jw[n2 + k];
	    jw[n2 + k] = j;
/*     exchange in w */
	    s = w[jj];
	    w[jj] = w[k];
	    w[k] = s;
	}

/*     zero out element in row by resetting jw(n+jrow) to zero. */

	jw[n + jrow] = 0;

/*     get the multiplier for row to be eliminated (jrow) + its level */

	fact = w[jj] * alu[jrow];
	jlev = jw[n2 + jj];
	if (jlev > lfil) {
	    goto L150;
	}

/*     combine current row and row jrow */

	i__2 = jlu[jrow + 1] - 1;
	for (k = ju[jrow]; k <= i__2; ++k) {
	    s = fact * alu[k];
	    j = jlu[k];
	    jpos = jw[n + j];
	    if (j >= ii) {

/*     dealing with upper part. */

		if (jpos == 0) {

/*     this is a fill-in element */

		    ++lenu;
		    if (lenu > n) {
			goto L995;
		    }
		    i__ = ii + lenu - 1;
		    jw[i__] = j;
		    jw[n + j] = i__;
		    w[i__] = -s;
		    jw[n2 + i__] = jlev + levs[k] + 1;
		} else {

/*     this is not a fill-in element */

		    w[jpos] -= s;
/* Computing MIN */
		    i__3 = jw[n2 + jpos], i__4 = jlev + levs[k] + 1;
		    jw[n2 + jpos] = (i__3 < i__4 ? i__3 : i__4);
		}
	    } else {

/*     dealing with lower part. */

		if (jpos == 0) {

/*     this is a fill-in element */

		    ++lenl;
		    if (lenl > n) {
			goto L995;
		    }
		    jw[lenl] = j;
		    jw[n + j] = lenl;
		    w[lenl] = -s;
		    jw[n2 + lenl] = jlev + levs[k] + 1;
		} else {

/*     this is not a fill-in element */

		    w[jpos] -= s;
/* Computing MIN */
		    i__3 = jw[n2 + jpos], i__4 = jlev + levs[k] + 1;
		    jw[n2 + jpos] = (i__3 < i__4 ? i__3 : i__4);
		}
	    }
/* L203: */
	}
	w[jj] = fact;
	jw[jj] = jrow;
	goto L150;
L160:

/*     reset double-pointer to zero (U-part) */

	i__2 = lenu;
	for (k = 1; k <= i__2; ++k) {
	    jw[n + jw[ii + k - 1]] = 0;
/* L308: */
	}

/*     update l-matrix */

	i__2 = lenl;
	for (k = 1; k <= i__2; ++k) {
	    if (ju0 > iwk) {
		goto L996;
	    }
	    if (jw[n2 + k] <= lfil) {
		alu[ju0] = w[k];
		jlu[ju0] = jw[k];
		++ju0;
	    }
/* L204: */
	}

/*     save pointer to beginning of row ii of U */

	ju[ii] = ju0;

/*     update u-matrix */

	i__2 = ii + lenu - 1;
	for (k = ii + 1; k <= i__2; ++k) {
	    if (jw[n2 + k] <= lfil) {
		jlu[ju0] = jw[k];
		alu[ju0] = w[k];
		levs[ju0] = jw[n2 + k];
		++ju0;
	    }
/* L302: */
	}
	if (fabs(w[ii]) < 1.0e-30) {
		printf("w[%lld]=%e\n",ii,w[ii]);
	    goto L999;
	}

	alu[ii] = 1.0 / w[ii];

/*     update pointer to beginning of next row of U. */

	jlu[ii + 1] = ju0;
/* ----------------------------------------------------------------------- */
/*     end main loop */
/* ----------------------------------------------------------------------- */
/* L500: */
    }
    
	++jw;
    ++w;
    ++ju;
    ++ia;
    ++a;
    ++ja;
    ++alu;
    ++jlu;
    ++levs;

	ierr = 0;
    return 0;

/*     incomprehensible error. Matrix must be wrong. */

L995:
	++jw;
    ++w;
    ++ju;
    ++ia;
    ++a;
    ++ja;
    ++alu;
    ++jlu;
    ++levs;

    ierr = -1;
    return 0;

/*     insufficient storage in L. */

L996:
	++jw;
    ++w;
    ++ju;
    ++ia;
    ++a;
    ++ja;
    ++alu;
    ++jlu;
    ++levs;

    ierr = -2;
    return 0;

/*     insufficient storage in U. */

/* L997: */
   // ierr = -3;
   // return 0;

/*     illegal lfil entered. */

L998:
	++jw;
    ++w;
    ++ju;
    ++ia;
    ++a;
    ++ja;
    ++alu;
    ++jlu;
    ++levs;

    ierr = -4;
    return 0;

/*     zero row encountered in A or U. */

L999:
	++jw;
    ++w;
    ++ju;
    ++ia;
    ++a;
    ++ja;
    ++alu;
    ++jlu;
    ++levs;

    ierr = -5;
    return 0;
/* ----------------end-of-iluk-------------------------------------------- */
/* ----------------------------------------------------------------------- */
} /* iluk_ */

integer iluk_(integer n, doublereal* &a, integer* &ja, integer* &ia,
	integer lfil, doublereal* &alu, integer* &jlu, integer* &ju,
	integer* &levs, integer iwk, doublereal* &w, integer* &jw, integer &ierr)
{
	integer ir = 0;
	//ir=iluk_Saad(n,a,ja,ia,lfil,alu,jlu,ju,levs,iwk,w,jw,ierr);
	
	// на основе быстродействующей хеш таблицы и двоичной кучи.
	ir=iluk_quick_stable(n, a, ja, ia, lfil, alu, jlu, ju, levs, iwk, w, jw, ierr);

	return ir;
}


integer iluk_2_serialpart(integer iscan_par, integer n2, integer &ju0, integer n, 
					   doublereal* &a, integer* &ja, integer* &ia,
					   integer lfil, doublereal* &alu, integer* &jlu, integer* &ju, 
	                   integer* &levs, integer iwk, doublereal* &w, integer* &jw, integer &ierr_loc)
{

	ierr_loc=0; // без ошибки.

	integer ii;
	doublereal fact;
	doublereal t;
	/* System generated locals */
    integer /*i__1,*/ i__2, i__3, i__4;

	 /* Local variables */
    integer i__, j, k;
    doublereal s;
    integer j1, j2,  jj;
	integer lenl, jlev, lenu, jpos, jrow;


	ii=iscan_par+1;

	j1 = ia[ii];
	j2 = ia[ii + 1] - 1;

/*     unpack L-part and U-part of row of A in arrays w */

	lenu = 1;
	lenl = 0;
	jw[ii] = ii;
	w[ii] = 0.f;
	jw[n + ii] = ii;

	i__2 = j2;
	for (j = j1; j <= i__2; ++j) {
	    k = ja[j];
	    t = a[j];
	    if (t == 0.f) {
		goto L170;
	    }
	    if (k < ii) {
		++lenl;
		jw[lenl] = k;
		w[lenl] = t;
		jw[n2 + lenl] = 0;
		jw[n + k] = lenl;
	    } else if (k == ii) {
		w[ii] = t;
		jw[n2 + ii] = 0;
	    } else {
		++lenu;
		jpos = ii + lenu - 1;
		jw[jpos] = k;
		w[jpos] = t;
		jw[n2 + jpos] = 0;
		jw[n + k] = jpos;
	    }
L170:
	    ;
	}

	jj = 0;

/*     eliminate previous rows */

L150:
	++jj;
	if (jj > lenl) {
	    goto L160;
	}
/* ----------------------------------------------------------------------- */
/*     in order to do the elimination in the correct order we must select */
/*     the smallest column index among jw(k), k=jj+1, ..., lenl. */
/* ----------------------------------------------------------------------- */
	jrow = jw[jj];
	k = jj;

/*     determine smallest column index */

	i__2 = lenl;
	for (j = jj + 1; j <= i__2; ++j) {
	    if (jw[j] < jrow) {
		jrow = jw[j];
		k = j;
	    }
/* L151: */
	}

	if (k != jj) {
/*     exchange in jw */
	    j = jw[jj];
	    jw[jj] = jw[k];
	    jw[k] = j;
/*     exchange in jw(n+  (pointers/ nonzero indicator). */
	    jw[n + jrow] = jj;
	    jw[n + j] = k;
/*     exchange in jw(n2+  (levels) */
	    j = jw[n2 + jj];
	    jw[n2 + jj] = jw[n2 + k];
	    jw[n2 + k] = j;
/*     exchange in w */
	    s = w[jj];
	    w[jj] = w[k];
	    w[k] = s;
	}

/*     zero out element in row by resetting jw(n+jrow) to zero. */

	jw[n + jrow] = 0;

/*     get the multiplier for row to be eliminated (jrow) + its level */

	fact = w[jj] * alu[jrow];
	jlev = jw[n2 + jj];
	if (jlev > lfil) {
	    goto L150;
	}

/*     combine current row and row jrow */

	i__2 = jlu[jrow + 1] - 1;
	for (k = ju[jrow]; k <= i__2; ++k) {
	    s = fact * alu[k];
	    j = jlu[k];
	    jpos = jw[n + j];
	    if (j >= ii) {

/*     dealing with upper part. */

		if (jpos == 0) {

/*     this is a fill-in element */

		    ++lenu;
		    if (lenu > n) {
			goto L995;
		    }
		    i__ = ii + lenu - 1;
		    jw[i__] = j;
		    jw[n + j] = i__;
		    w[i__] = -s;
		    jw[n2 + i__] = jlev + levs[k] + 1;
		} else {

/*     this is not a fill-in element */

		    w[jpos] -= s;
/* Computing MIN */
		    i__3 = jw[n2 + jpos], i__4 = jlev + levs[k] + 1;
		    jw[n2 + jpos] = (i__3 < i__4 ? i__3 : i__4);
		}
	    } else {

/*     dealing with lower part. */

		if (jpos == 0) {

/*     this is a fill-in element */

		    ++lenl;
		    if (lenl > n) {
			goto L995;
		    }
		    jw[lenl] = j;
		    jw[n + j] = lenl;
		    w[lenl] = -s;
		    jw[n2 + lenl] = jlev + levs[k] + 1;
		} else {

/*     this is not a fill-in element */

		    w[jpos] -= s;
/* Computing MIN */
		    i__3 = jw[n2 + jpos], i__4 = jlev + levs[k] + 1;
		    jw[n2 + jpos] = (i__3 < i__4 ? i__3 : i__4);
		}
	    }
/* L203: */
	}
	w[jj] = fact;
	jw[jj] = jrow;
	goto L150;
L160:

/*     reset double-pointer to zero (U-part) */

	i__2 = lenu;
	for (k = 1; k <= i__2; ++k) {
	    jw[n + jw[ii + k - 1]] = 0;
/* L308: */
	}

/*     update l-matrix */

	i__2 = lenl;
	for (k = 1; k <= i__2; ++k) {
	    if (ju0 > iwk) {
		goto L996;
	    }
	    if (jw[n2 + k] <= lfil) {
		alu[ju0] = w[k];
		jlu[ju0] = jw[k];
		++ju0;
	    }
/* L204: */
	}

/*     save pointer to beginning of row ii of U */

	ju[ii] = ju0;

/*     update u-matrix */

	i__2 = ii + lenu - 1;
	for (k = ii + 1; k <= i__2; ++k) {
	    if (jw[n2 + k] <= lfil) {
		jlu[ju0] = jw[k];
		alu[ju0] = w[k];
		levs[ju0] = jw[n2 + k];
		++ju0;
	    }
/* L302: */
	}
	if (fabs(w[ii]) < 1.0e-30) {
		printf("w[%lld]=%e\n", ii, w[ii]);
	    goto L999;
	}

	alu[ii] = 1.0 / w[ii];

/*     update pointer to beginning of next row of U. */

	jlu[ii + 1] = ju0;
/* ----------------------------------------------------------------------- */
/*     end main loop */
/* ----------------------------------------------------------------------- */

	return 0;

L995:
	// ierr==-1
	ierr_loc=-1;
	return 0;

L996:
	// ierr==-2
	ierr_loc=-2;
	return 0;

L999:
	// ierr==-5
	ierr_loc=-5;
	return 0;
}


bool BiCGStab_internal3_incomming_now = false;

// 13 августа 2015 года распараллелил iluk алгоритм из библиотеки SPARSKIT2
// Й. Саада на два потока. Это дало ускорение одной итерации SIMPLE алгоритма 
// на 18.4%. 
/* ----------------------------------------------------------------------- */
/* Subroutine */ integer iluk_2_work(integer n, doublereal*& a, integer*& ja, integer*& ia,
	integer lfil, doublereal*& alu, integer*& jlu, integer*& ju,
	integer*& levs, integer iwk, doublereal*& w, integer*& jw, doublereal*& w_dubl, integer*& jw_dubl, integer& ierr)
{

	/* Function Body */
	if (lfil < 0) {

		/*     illegal lfil entered. */

		/*
		++jw;
		++w;
		++ju;
		++ia;
		++a;
		++ja;
		++alu;
		++jlu;
		++levs;
		*/

		ierr = -4;
		return 0;
	}
	else {

		/* System generated locals */
		integer i__1/*, i__2,*/ /*i__3,*//* i__4*/;

		/* Local variables */
		integer /*i__,*/ j /*k*/;
		// doublereal /*s,*/ t;
		integer /*j1, j2,*/ n2,/* ii,*/ /*jj,*/ ju0;
		// doublereal fact;
		// integer /*lenl,*/ /*jlev,*/ lenu /*jpos, jrow*/;

	 /* ----------------------------------------------------------------------* */
	 /*     SPARSKIT ROUTINE ILUK -- ILU WITH LEVEL OF FILL-IN OF K (ILU(k)) * */
	 /* ----------------------------------------------------------------------* */

	 /* on entry: */
	 /* ========== */
	 /* n       = integer. The row dimension of the matrix A. The matrix */

	 /* a,ja,ia = matrix stored in Compressed Sparse Row format. */

	 /* lfil    = integer. The fill-in parameter. Each element whose */
	 /*           leve-of-fill exceeds lfil during the ILU process is dropped. */
	 /*           lfil must be .ge. 0 */

	 /* tol     = real*8. Sets the threshold for dropping small terms in the */
	 /*           factorization. See below for details on dropping strategy. */

	 /* iwk     = integer. The minimum length of arrays alu, jlu, and levs. */

	 /* On return: */
	 /* =========== */

	 /* alu,jlu = matrix stored in Modified Sparse Row (MSR) format containing */
	 /*           the L and U factors together. The diagonal (stored in */
	 /*           alu(1:n) ) is inverted. Each i-th row of the alu,jlu matrix */
	 /*           contains the i-th row of L (excluding the diagonal entry=1) */
	 /*           followed by the i-th row of U. */

	 /* ju      = integer array of length n containing the pointers to */
	 /*           the beginning of each row of U in the matrix alu,jlu. */

	 /* levs    = integer (work) array of size iwk -- which contains the */
	 /*           levels of each element in alu, jlu. */

	 /* ierr    = integer. Error message with the following meaning. */
	 /*           ierr  = 0    --> successful return. */
	 /*           ierr .gt. 0  --> zero pivot encountered at step number ierr. */
	 /*           ierr  = -1   --> Error. input matrix may be wrong. */
	 /*                            (The elimination process has generated a */
	 /*                            row in L or U whose length is .gt.  n.) */
	 /*           ierr  = -2   --> The matrix L overflows the array al. */
	 /*           ierr  = -3   --> The matrix U overflows the array alu. */
	 /*           ierr  = -4   --> Illegal value for lfil. */
	 /*           ierr  = -5   --> zero row encountered in A or U. */

	 /* work arrays: */
	 /* ============= */
	 /* jw      = integer work array of length 3*n. */
	 /* w       = real work array of length n */

	 /* Notes/known bugs: This is not implemented efficiently storage-wise. */
	 /*       For example: Only the part of the array levs(*) associated with */
	 /*       the U-matrix is needed in the routine.. So some storage can */
	 /*       be saved if needed. The levels of fills in the LU matrix are */
	 /*       output for information only -- they are not needed by LU-solve. */

	 /* ---------------------------------------------------------------------- */
	 /* w, ju (1:n) store the working array [1:ii-1 = L-part, ii:n = u] */
	 /* jw(n+1:2n)  stores the nonzero indicator. */

	 /* Notes: */
	 /* ------ */
	 /* All the diagonal elements of the input matrix must be  nonzero. */

	 /* ----------------------------------------------------------------------* */
	 /*     locals */
		 /* Parameter adjustments */
		--jw;
		--w;
		--ju;
		--ia;
		--a;
		--ja;
		--alu;
		--jlu;
		--levs;


		/* ----------------------------------------------------------------------- */
		/*     initialize ju0 (points to next element to be added to alu,jlu) */
		/*     and pointer array. */
		/* ----------------------------------------------------------------------- */
		n2 = n + n;
		ju0 = n + 2;
		jlu[1] = ju0;

		/*     initialize nonzero indicator array + levs array -- */

		i__1 = n << 1;
		for (j = 1; j <= i__1; ++j) {
			jw[j] = 0;
			/* L1: */
		}

		//if (bparallelizm_old) {
		if (BiCGStab_internal3_incomming_now&&(number_cores() == 2)&&((my_amg_manager.lfil < 3))) {
			//ILUK 839 +2
			//doublereal *w_dubl=new doublereal[n+2]; // +2 запас по памяти.
			for (integer i35 = 0; i35 < n+2; i35++) w_dubl[i35] = w[i35];

			//integer *jw_dubl=new integer[3*n+2];
			for (integer i35 = 0; i35 < 3 * n+2; i35++) jw_dubl[i35] = jw[i35];
		}

		//integer debug1=0;
		integer ierr_messag = ierr;

#ifdef _OPENMP
		//if (bparallelizm_old) {
			//if (inumcore == 2) 
		// двухпоточный код 18.01.2020
		if (BiCGStab_internal3_incomming_now && (number_cores() == 2)) {
			if (my_amg_manager.lfil<3) {

				int n_thread = omp_get_num_threads();
				omp_set_num_threads(2);
				

				if (nd.b0.active) {
					integer ju0gl;
					integer ju02start;


					// Для второго смыкающего потока.
					integer r87 = (nd.b0.ileft_finish - nd.b0.ileft_start);
					integer r88 = nd.b0.iright_finish - nd.b0.iright_start + nd.b0.ileft_finish - nd.b0.ileft_start;
					ju02start = (integer)(iwk * ((1.0 * r87) / (1.0 * r88)));
					ju0gl = ju02start;


#pragma omp parallel shared(ju0gl,ju02start,ierr_messag)
					{
#pragma omp parallel sections 
						{
#pragma omp section
							{

								integer ierr_loc = 0;


								// первый поток
								for (integer iscan_par = nd.b0.ileft_start; iscan_par <= nd.b0.ileft_finish; iscan_par++) {

									iluk_2_serialpart(iscan_par, n2, ju0, n,
										a, ja, ia,
										lfil, alu, jlu, ju,
										levs, iwk, w, jw, ierr_loc);

									if (ierr_loc < 0) {
										break;
									}

								}



								if (ierr_loc < 0) {

#pragma omp critical
									{

										ierr_messag = ierr_loc;
									}
								}


							}

#pragma omp section
							{
								// второй поток
								integer ierr_loc = 0;





								jlu[nd.b0.iright_start + 1] = ju02start; // обязательно важно.



								for (integer iscan_par = nd.b0.iright_start; iscan_par <= nd.b0.iright_finish; iscan_par++) {


									iluk_2_serialpart(iscan_par, n2, ju0gl, n,
										a, ja, ia,
										lfil, alu, jlu, ju,
										levs, iwk, w_dubl, jw_dubl, ierr_loc);

									if (ierr_loc < 0) {
										break;
									}

								}



								if (ierr_loc < 0) {

#pragma omp critical
									{

										ierr_messag = ierr_loc;
									}
								}

							} // end section
						} // end sections
					} // end parallel




								// серийный смыкающий кусок
					if (ierr_messag == 0) {


						for (integer iscan_par = nd.b0.iright_start; iscan_par <= nd.b0.iright_finish; iscan_par++) {
							jw[iscan_par + 1] = jw_dubl[iscan_par + 1]; // модификация рабочего массива.
							ju[iscan_par + 1] = ju[iscan_par + 1] - ju02start + ju0;
							jlu[iscan_par + 1] -= ju02start - ju0;
						}
						jlu[nd.b0.iright_finish + 2] = ju[nd.b0.iright_finish + 1]; // Вот так наверно совсем правильно.
						jw[nd.b0.iright_finish + 2] = jw[nd.b0.iright_finish + 1];

						for (integer i87 = ju02start; i87 <= ju0gl; i87++) {
							// no move
							//alu[ju0]=alu[i87];
							//jlu[ju0]=jlu[i87];
							//levs[ju0]=levs[i87];

							// need move and swap

							integer i1 = ju0;
							doublereal abub = alu[i1];
							alu[i1] = alu[i87];
							alu[i87] = abub;

							integer ibuf = jlu[i1];
							jlu[i1] = jlu[i87];
							jlu[i87] = ibuf;

							ibuf = levs[i1];
							levs[i1] = levs[i87];
							levs[i87] = ibuf;


							ju0++;
						}


						for (integer iscan_par = nd.b0.iseparate_start; iscan_par <= nd.b0.iseparate_finish; iscan_par++) {

							iluk_2_serialpart(iscan_par, n2, ju0, n,
								a, ja, ia,
								lfil, alu, jlu, ju,
								levs, iwk, w, jw, ierr_messag);


						}
					}

					omp_set_num_threads(n_thread);

				}
			}
		}
		else {
			/* ----------------------------------------------------------------------- */
			/*     beginning of main loop. */
			/* ----------------------------------------------------------------------- */
			i__1 = n;
			for (integer ii = 1; ii <= i__1; ++ii) {

				integer iscan_par = ii - 1;

				iluk_2_serialpart(iscan_par, n2, ju0, n,
					a, ja, ia,
					lfil, alu, jlu, ju,
					levs, iwk, w, jw, ierr_messag);
				if (ierr_messag != 0) {
					break; // досрочный выход из цикла for.
				}

				/* L500: */
			}
		}
#else


		/* ----------------------------------------------------------------------- */
		/*     beginning of main loop. */
		/* ----------------------------------------------------------------------- */
		i__1 = n;
		for (integer ii = 1; ii <= i__1; ++ii) {

			integer iscan_par = ii - 1;

			iluk_2_serialpart(iscan_par, n2, ju0, n,
				a, ja, ia,
				lfil, alu, jlu, ju,
				levs, iwk, w, jw, ierr_messag);
			if (ierr_messag != 0) {
				break; // досрочный выход из цикла for.
			}

			/* L500: */
		}


#endif

		switch (ierr_messag) {
		case -1: goto L995; break;
		case -2: goto L996; break;
		case -5: goto L999; break;
		}

		//delete w_dubl; 
		//delete jw_dubl;

		++jw;
		++w;
		++ju;
		++ia;
		++a;
		++ja;
		++alu;
		++jlu;
		++levs;

		ierr = 0;
		return 0;

		/*     incomprehensible error. Matrix must be wrong. */

	L995:
		++jw;
		++w;
		++ju;
		++ia;
		++a;
		++ja;
		++alu;
		++jlu;
		++levs;

		ierr = -1;
		return 0;

		/*     insufficient storage in L. */

	L996:
		++jw;
		++w;
		++ju;
		++ia;
		++a;
		++ja;
		++alu;
		++jlu;
		++levs;

		ierr = -2;
		return 0;

		/*     insufficient storage in U. */

		/* L997: */
		//ierr = -3;
		//return 0;





		/*     zero row encountered in A or U. */

	L999:
		++jw;
		++w;
		++ju;
		++ia;
		++a;
		++ja;
		++alu;
		++jlu;
		++levs;

		ierr = -5;
		return 0;
		/* ----------------end-of-iluk-------------------------------------------- */
		/* ----------------------------------------------------------------------- */

	}

} /* iluk_2 */

/* Subroutine */ integer iluk_2(integer n, doublereal*& a, integer*& ja, integer*& ia,
	integer lfil, doublereal*& alu, integer*& jlu, integer*& ju,
	integer*& levs, integer iwk, doublereal*& w, integer*& jw, doublereal*& w_dubl, integer*& jw_dubl, integer& ierr)
{
	integer ireturn_state = 0;
	if (1) {
		
		//ireturn_state=iluk_Saad(n,a,ja,ia,lfil,alu,jlu,ju,levs,iwk,w,jw,ierr);

		// на основе быстродействующей хеш таблицы и двоичной кучи.
		ireturn_state = iluk_quick_stable(n, a, ja, ia, lfil, alu, jlu, ju, levs, iwk, w, jw, ierr);

		
	}
	else {
		ireturn_state =iluk_2_work( n,  a, ja, ia, lfil, alu, jlu, ju, levs, iwk, w, jw, w_dubl, jw_dubl, ierr);
	}
	return ireturn_state;
}

// 13 августа 2015 года распараллелил iluk алгоритм из библиотеки SPARSKIT2
// Й. Саада на два потока. Это дало ускорение одной итерации SIMPLE алгоритма 
// на 18.4%. 
/* ----------------------------------------------------------------------- */
/* Subroutine */ integer iluk_2_memory_18_01_2020(integer n, doublereal* &a, integer* &ja, integer* &ia,
						   integer lfil, doublereal* &alu, integer* &jlu, integer* &ju, 
	                       integer* &levs, integer iwk, doublereal* &w, integer* &jw,
	                       doublereal* &w_dubl, integer* &jw_dubl, integer &ierr)
{

	/* Function Body */
	if (lfil < 0) {

		/*     illegal lfil entered. */

		/*
		++jw;
		++w;
		++ju;
		++ia;
		++a;
		++ja;
		++alu;
		++jlu;
		++levs;
		*/

		ierr = -4;
		return 0;
	}
	else {

    /* System generated locals */
    integer i__1/*, i__2,*/ /*i__3,*//* i__4*/;

    /* Local variables */
    integer /*i__,*/ j /*k*/;
   // doublereal /*s,*/ t;
    integer /*j1, j2,*/ n2,/* ii,*/ /*jj,*/ ju0;
   // doublereal fact;
   // integer /*lenl,*/ /*jlev,*/ lenu /*jpos, jrow*/;

/* ----------------------------------------------------------------------* */
/*     SPARSKIT ROUTINE ILUK -- ILU WITH LEVEL OF FILL-IN OF K (ILU(k)) * */
/* ----------------------------------------------------------------------* */

/* on entry: */
/* ========== */
/* n       = integer. The row dimension of the matrix A. The matrix */

/* a,ja,ia = matrix stored in Compressed Sparse Row format. */

/* lfil    = integer. The fill-in parameter. Each element whose */
/*           leve-of-fill exceeds lfil during the ILU process is dropped. */
/*           lfil must be .ge. 0 */

/* tol     = real*8. Sets the threshold for dropping small terms in the */
/*           factorization. See below for details on dropping strategy. */

/* iwk     = integer. The minimum length of arrays alu, jlu, and levs. */

/* On return: */
/* =========== */

/* alu,jlu = matrix stored in Modified Sparse Row (MSR) format containing */
/*           the L and U factors together. The diagonal (stored in */
/*           alu(1:n) ) is inverted. Each i-th row of the alu,jlu matrix */
/*           contains the i-th row of L (excluding the diagonal entry=1) */
/*           followed by the i-th row of U. */

/* ju      = integer array of length n containing the pointers to */
/*           the beginning of each row of U in the matrix alu,jlu. */

/* levs    = integer (work) array of size iwk -- which contains the */
/*           levels of each element in alu, jlu. */

/* ierr    = integer. Error message with the following meaning. */
/*           ierr  = 0    --> successful return. */
/*           ierr .gt. 0  --> zero pivot encountered at step number ierr. */
/*           ierr  = -1   --> Error. input matrix may be wrong. */
/*                            (The elimination process has generated a */
/*                            row in L or U whose length is .gt.  n.) */
/*           ierr  = -2   --> The matrix L overflows the array al. */
/*           ierr  = -3   --> The matrix U overflows the array alu. */
/*           ierr  = -4   --> Illegal value for lfil. */
/*           ierr  = -5   --> zero row encountered in A or U. */

/* work arrays: */
/* ============= */
/* jw      = integer work array of length 3*n. */
/* w       = real work array of length n */

/* Notes/known bugs: This is not implemented efficiently storage-wise. */
/*       For example: Only the part of the array levs(*) associated with */
/*       the U-matrix is needed in the routine.. So some storage can */
/*       be saved if needed. The levels of fills in the LU matrix are */
/*       output for information only -- they are not needed by LU-solve. */

/* ---------------------------------------------------------------------- */
/* w, ju (1:n) store the working array [1:ii-1 = L-part, ii:n = u] */
/* jw(n+1:2n)  stores the nonzero indicator. */

/* Notes: */
/* ------ */
/* All the diagonal elements of the input matrix must be  nonzero. */

/* ----------------------------------------------------------------------* */
/*     locals */
    /* Parameter adjustments */
    --jw;
    --w;
    --ju;
    --ia;
    --a;
    --ja;
    --alu;
    --jlu;
    --levs;

    
/* ----------------------------------------------------------------------- */
/*     initialize ju0 (points to next element to be added to alu,jlu) */
/*     and pointer array. */
/* ----------------------------------------------------------------------- */
    n2 = n + n;
    ju0 = n + 2;
    jlu[1] = ju0;

/*     initialize nonzero indicator array + levs array -- */

    i__1 = n << 1;
    for (j = 1; j <= i__1; ++j) {
	jw[j] = 0;
/* L1: */
    }

	if (bparallelizm_old) {
		//ILUK 839 +2
		//doublereal *w_dubl=new doublereal[n+2]; // +2 запас по памяти.
		for (integer i35 = 0; i35 < n; i35++) w_dubl[i35] = w[i35];

		//integer *jw_dubl=new integer[3*n+2];
		for (integer i35 = 0; i35 < 3 * n; i35++) jw_dubl[i35] = jw[i35];
	}

	//integer debug1=0;
	integer ierr_messag = ierr;

	#ifdef _OPENMP
	if (bparallelizm_old) {
		if (inumcore == 2) {
			if (nd.b0.active) {
				integer ju0gl;
				integer ju02start;


				// Для второго смыкающего потока.
				integer r87 = (nd.b0.ileft_finish - nd.b0.ileft_start);
				integer r88 = nd.b0.iright_finish - nd.b0.iright_start + nd.b0.ileft_finish - nd.b0.ileft_start;
				ju02start = (integer)(iwk*((1.0*r87) / (1.0*r88)));
				ju0gl = ju02start;


#pragma omp parallel shared(ju0gl,ju02start,ierr_messag)
				{
#pragma omp parallel sections 
					{
#pragma omp section
						{

							integer ierr_loc = 0;


							// первый поток
							for (integer iscan_par = nd.b0.ileft_start; iscan_par <= nd.b0.ileft_finish; iscan_par++) {

								iluk_2_serialpart(iscan_par, n2, ju0, n,
									a, ja, ia,
									lfil, alu, jlu, ju,
									levs, iwk, w, jw, ierr_loc);

								if (ierr_loc < 0) {
									break;
								}

							}



							if (ierr_loc < 0) {

#pragma omp critical
								{

									ierr_messag = ierr_loc;
								}
							}


						}

#pragma omp section
						{
							// второй поток
							integer ierr_loc = 0;





							jlu[nd.b0.iright_start + 1] = ju02start; // обязательно важно.



							for (integer iscan_par = nd.b0.iright_start; iscan_par <= nd.b0.iright_finish; iscan_par++) {


								iluk_2_serialpart(iscan_par, n2, ju0gl, n,
									a, ja, ia,
									lfil, alu, jlu, ju,
									levs, iwk, w_dubl, jw_dubl, ierr_loc);

								if (ierr_loc < 0) {
									break;
								}

							}



							if (ierr_loc < 0) {

#pragma omp critical
								{

									ierr_messag = ierr_loc;
								}
							}

						} // end section
					} // end sections
				} // end parallel




							// серийный смыкающий кусок
				if (ierr_messag == 0) {


					for (integer iscan_par = nd.b0.iright_start; iscan_par <= nd.b0.iright_finish; iscan_par++) {
						jw[iscan_par + 1] = jw_dubl[iscan_par + 1]; // модификация рабочего массива.
						ju[iscan_par + 1] = ju[iscan_par + 1] - ju02start + ju0;
						jlu[iscan_par + 1] -= ju02start - ju0;
					}
					jlu[nd.b0.iright_finish + 2] = ju[nd.b0.iright_finish + 1]; // Вот так наверно совсем правильно.
					jw[nd.b0.iright_finish + 2] = jw[nd.b0.iright_finish + 1];

					for (integer i87 = ju02start; i87 <= ju0gl; i87++) {
						// no move
						//alu[ju0]=alu[i87];
						//jlu[ju0]=jlu[i87];
						//levs[ju0]=levs[i87];

						// need move and swap

						integer i1 = ju0;
						doublereal abub = alu[i1];
						alu[i1] = alu[i87];
						alu[i87] = abub;

						integer ibuf = jlu[i1];
						jlu[i1] = jlu[i87];
						jlu[i87] = ibuf;

						ibuf = levs[i1];
						levs[i1] = levs[i87];
						levs[i87] = ibuf;


						ju0++;
					}


					for (integer iscan_par = nd.b0.iseparate_start; iscan_par <= nd.b0.iseparate_finish; iscan_par++) {

						iluk_2_serialpart(iscan_par, n2, ju0, n,
							a, ja, ia,
							lfil, alu, jlu, ju,
							levs, iwk, w, jw, ierr_messag);


					}
				}

			}
		}
	}
	else {
	   /* ----------------------------------------------------------------------- */
       /*     beginning of main loop. */
       /* ----------------------------------------------------------------------- */
       i__1 = n;
       for (integer ii = 1; ii <= i__1; ++ii) {

	        integer iscan_par = ii - 1;

	        iluk_2_serialpart(iscan_par, n2, ju0, n,
	         	a, ja, ia,
	        	lfil, alu, jlu, ju,
	        	levs, iwk, w, jw, ierr_messag);
	         if (ierr_messag != 0) {
	              break; // досрочный выход из цикла for.
	         }

	       /* L500: */
        }
    }
#else


	/* ----------------------------------------------------------------------- */
	/*     beginning of main loop. */
	/* ----------------------------------------------------------------------- */
	i__1 = n;
	for (integer ii = 1; ii <= i__1; ++ii) {

		integer iscan_par = ii - 1;

		iluk_2_serialpart(iscan_par, n2, ju0, n,
			a, ja, ia,
			lfil, alu, jlu, ju,
			levs, iwk, w, jw, ierr_messag);
		if (ierr_messag != 0) {
			break; // досрочный выход из цикла for.
		}

		/* L500: */
						}


#endif

	switch (ierr_messag) {
	case -1: goto L995; break;
	case -2: goto L996; break;
	case -5: goto L999; break;
	}

	//delete w_dubl; 
	//delete jw_dubl;

	++jw;
	++w;
	++ju;
	++ia;
	++a;
	++ja;
	++alu;
	++jlu;
	++levs;

	ierr = 0;
	return 0;

	/*     incomprehensible error. Matrix must be wrong. */

L995:
	++jw;
	++w;
	++ju;
	++ia;
	++a;
	++ja;
	++alu;
	++jlu;
	++levs;

	ierr = -1;
	return 0;

	/*     insufficient storage in L. */

L996:
	++jw;
	++w;
	++ju;
	++ia;
	++a;
	++ja;
	++alu;
	++jlu;
	++levs;

	ierr = -2;
	return 0;

	/*     insufficient storage in U. */

	/* L997: */
	//ierr = -3;
	//return 0;





	/*     zero row encountered in A or U. */

L999:
	++jw;
	++w;
	++ju;
	++ia;
	++a;
	++ja;
	++alu;
	++jlu;
	++levs;

	ierr = -5;
	return 0;
	/* ----------------end-of-iluk-------------------------------------------- */
	/* ----------------------------------------------------------------------- */

					}

} /* iluk_2_memory_18.01.2020 */


/*
/* ----------------------------------------------------------------------- */
/* Subroutine */ /*integer iluk_2(integer n, doublereal* &a, integer* &ja, integer* &ia,
						   integer lfil, doublereal* &alu, integer* &jlu, integer* &ju, 
	integer* &levs, integer iwk, doublereal* &w, integer* &jw, integer &ierr)
{*/
    /* System generated locals */
   // integer i__1/*, i__2,*/ /*i__3,*//* i__4*/;

    /* Local variables */
    //integer /*i__,*/ j /*k*/;
   // doublereal /*s,*/ t;
    //integer /*j1, j2,*/ n2,/* ii,*/ /*jj,*/ ju0;
   // doublereal fact;
   // integer /*lenl,*/ /*jlev,*/ lenu /*jpos, jrow*/;

/* ----------------------------------------------------------------------* */
/*     SPARSKIT ROUTINE ILUK -- ILU WITH LEVEL OF FILL-IN OF K (ILU(k)) * */
/* ----------------------------------------------------------------------* */

/* on entry: */
/* ========== */
/* n       = integer. The row dimension of the matrix A. The matrix */

/* a,ja,ia = matrix stored in Compressed Sparse Row format. */

/* lfil    = integer. The fill-in parameter. Each element whose */
/*           leve-of-fill exceeds lfil during the ILU process is dropped. */
/*           lfil must be .ge. 0 */

/* tol     = real*8. Sets the threshold for dropping small terms in the */
/*           factorization. See below for details on dropping strategy. */

/* iwk     = integer. The minimum length of arrays alu, jlu, and levs. */

/* On return: */
/* =========== */

/* alu,jlu = matrix stored in Modified Sparse Row (MSR) format containing */
/*           the L and U factors together. The diagonal (stored in */
/*           alu(1:n) ) is inverted. Each i-th row of the alu,jlu matrix */
/*           contains the i-th row of L (excluding the diagonal entry=1) */
/*           followed by the i-th row of U. */

/* ju      = integer array of length n containing the pointers to */
/*           the beginning of each row of U in the matrix alu,jlu. */

/* levs    = integer (work) array of size iwk -- which contains the */
/*           levels of each element in alu, jlu. */

/* ierr    = integer. Error message with the following meaning. */
/*           ierr  = 0    --> successful return. */
/*           ierr .gt. 0  --> zero pivot encountered at step number ierr. */
/*           ierr  = -1   --> Error. input matrix may be wrong. */
/*                            (The elimination process has generated a */
/*                            row in L or U whose length is .gt.  n.) */
/*           ierr  = -2   --> The matrix L overflows the array al. */
/*           ierr  = -3   --> The matrix U overflows the array alu. */
/*           ierr  = -4   --> Illegal value for lfil. */
/*           ierr  = -5   --> zero row encountered in A or U. */

/* work arrays: */
/* ============= */
/* jw      = integer work array of length 3*n. */
/* w       = real work array of length n */

/* Notes/known bugs: This is not implemented efficiently storage-wise. */
/*       For example: Only the part of the array levs(*) associated with */
/*       the U-matrix is needed in the routine.. So some storage can */
/*       be saved if needed. The levels of fills in the LU matrix are */
/*       output for information only -- they are not needed by LU-solve. */

/* ---------------------------------------------------------------------- */
/* w, ju (1:n) store the working array [1:ii-1 = L-part, ii:n = u] */
/* jw(n+1:2n)  stores the nonzero indicator. */

/* Notes: */
/* ------ */
/* All the diagonal elements of the input matrix must be  nonzero. */

/* ----------------------------------------------------------------------* */
/*     locals */
    /* Parameter adjustments */
	/*
    --jw;
    --w;
    --ju;
    --ia;
    --a;
    --ja;
    --alu;
    --jlu;
    --levs;

    /* Function Body *//*
    if (lfil < 0) {
	goto L998;
    }*/
/* ----------------------------------------------------------------------- */
/*     initialize ju0 (points to next element to be added to alu,jlu) */
/*     and pointer array. */
/* ----------------------------------------------------------------------- */
  /*  n2 = n + n;
    ju0 = n + 2;
    jlu[1] = ju0;*/

/*     initialize nonzero indicator array + levs array -- */
	/*
    i__1 = n << 1;
    for (j = 1; j <= i__1; ++j) {
	jw[j] = 0;

    }


	integer debug1=0;
	integer ierr_messag=ierr;

	#ifdef _OPENMP

			if (inumcore==2) {
				if (nd.b0.active) {
					 integer ju0gl;
					 integer ju02start;

					 integer r87=(nd.b0.ileft_finish-nd.b0.ileft_start);
					 integer r88=nd.b0.iright_finish-nd.b0.iright_start+nd.b0.ileft_finish-nd.b0.ileft_start;
					 ju02start=(integer)(iwk*((1.0*r87)/(1.0*r88)));
					 ju0gl=ju02start;

					 integer ideb=0;

#pragma omp parallel shared(ju0gl,ju02start)
					{
#pragma omp parallel sections 
						{
#pragma omp section
							{

								integer ierr_loc=0;
							

					// первый поток
					for (integer iscan_par=nd.b0.ileft_start; iscan_par<=nd.b0.ileft_finish; iscan_par++) {
						
                    iluk_2_serialpart(iscan_par, n2, ju0, n, 
					   a, ja, ia,
					    lfil, alu, jlu, ju, 
	                   levs, iwk, w, jw, ierr_loc);
						
					if (ierr_loc<0) {
						break;
					}

					}

					debug1=ju0;
					ideb=ju0;

					if (ierr_loc<0) {
						ierr_messag=ierr_loc;
					}


					}
							*/
							/*
									// debug
	FILE *fdeb=NULL;
					errno_t err87;
	err87 = fopen_s( &fdeb, "lognow.txt", "w");
	//for (integer i54=0; i54<n+2; i54++) {
		//fprintf(fdeb,"%e \n",w[i54]);
//	}
	// Полный лог всех переменных.
	
	
	for (integer i54=0; i54<iwk+1; i54++) {
		if (i54<n+1) {
			fprintf(fdeb,"%e %d %d %d %e %d\n",alu[i54],jlu[i54],levs[i54], ju[i54], w[i54], jw[i54]);
		}
		else if (i54<3*n+1) {
			fprintf(fdeb,"%e %d %d NO NO %d\n",alu[i54],jlu[i54],levs[i54], jw[i54]);
		}
		else {
			fprintf(fdeb,"%e %d %d NO NO NO\n",alu[i54],jlu[i54],levs[i54]);
		}
	}
	fprintf(fdeb,"alu   jlu    levs    ju    w    jw\n");
	fprintf(fdeb,"ileftstart=%d, ileftfinish=%d, irightstart=%d, irightfinish=%d, isepstart=%d, isepfinish=%d\n",nd.b0.ileft_start,nd.b0.ileft_finish,nd.b0.iright_start,nd.b0.iright_finish, nd.b0.iseparate_start,nd.b0.iseparate_finish);
	fprintf(fdeb,"ju0finish1=%d, ju0finish2=%d\n",ideb, ju0);
	printf("hok_pre\n");
	fclose(fdeb);
	getchar();
	*/
/*
							//printf("ju0=%d jluls=%d jlulf=%d jlurs=%d jlurf=%d jluss=%d jlusf=%d\n",ju0,jlu[nd.b0.ileft_start+1],jlu[nd.b0.ileft_finish+1],jlu[nd.b0.iright_start+1],jlu[nd.b0.iright_finish+1],jlu[nd.b0.iseparate_start+1],jlu[nd.b0.iseparate_finish+1]);
							//getchar();
#pragma omp section
					{
					// второй поток
						integer ierr_loc=0;
						//integer r87=(nd.b0.ileft_finish-nd.b0.ileft_start);
						//integer r88=nd.b0.iright_finish-nd.b0.iright_start+nd.b0.ileft_finish-nd.b0.ileft_start;
						//ju02start=(integer)(iwk*((1.0*r87)/(1.0*r88)));
						ju0gl=ju02start;
						//printf("start 2=%d\n",ju0gl);
						//getchar();

						//printf("%d %d %d %d\n",ju0, jlu[ju0-1],jlu[ju0],jlu[ju0+1]);
						//printf("%d %d %d %d\n",ju02start, jlu[ju02start-1],jlu[ju02start],jlu[ju02start+1]);
						//getchar();

						 //ju0 = n + 2;
                        // jlu[nd.b0.iright_start+1] = ju0gl;
						//--->//jlu[ju02start]=ju0gl;
						//jlu[1] = ju0;
						//jlu[ju02start]=n+2;
						jlu[nd.b0.iright_start+1]=ju02start; // обязательно важно.


					//	printf("ju0=%d, ju02start=%d, iwk=%d %d %d\n",ju0,ju02start,iwk,nd.b0.iright_finish-nd.b0.iright_start,nd.b0.ileft_finish-nd.b0.ileft_start);
					//	printf("jlu[ileftstart+1]=%d, %d %d\n ",jlu[nd.b0.ileft_start+1],jlu[nd.b0.ileft_finish+1],jlu[nd.b0.ileft_finish+2]);

					for (integer iscan_par=nd.b0.iright_start; iscan_par<=nd.b0.iright_finish; iscan_par++) {
					

					iluk_2_serialpart(iscan_par, n2, ju0gl, n, 
					   a, ja, ia,
					    lfil, alu, jlu, ju, 
	                   levs, iwk, w, jw, ierr_loc);

					if (ierr_loc<0) {
						break;
					}

					}

					//printf("ju0gl=%d,iwk=%d\n",ju0gl,iwk);
					//getchar();

					if (ierr_loc<0) {
						ierr_messag=ierr_loc;
					}

				} // end section
			} // end sections
		} // end parallel*/
					//printf("ju0=%d jluls=%d jlulf=%d jlurs=%d jlurf=%d jluss=%d jlusf=%d\n",ju0,jlu[nd.b0.ileft_start+1],jlu[nd.b0.ileft_finish+1],jlu[nd.b0.iright_start+1],jlu[nd.b0.iright_finish+1],jlu[nd.b0.iseparate_start+1],jlu[nd.b0.iseparate_finish+1]);

					//printf("alu jlu\n");
					//for (integer i86=debug1-5; i86<debug1+6; i86++) {
				//		printf("%e, %d\n",alu[i86],jlu[i86]);
					//}
					//getchar();
		/*
					// серийный смыкающий кусок
					if (ierr_messag==0) {
					*/
						/*FILE *fdeb=NULL;
					errno_t err87;
	err87 = fopen_s( &fdeb, "lognow.txt", "w");
//	for (integer i54=0; i54<n+2; i54++) {
	//	fprintf(fdeb,"%e \n",w[i54]);
	//}
	for (integer i54=0; i54<iwk+2; i54++) {
		fprintf(fdeb,"%d \n",jlu[i54]);
	}
	fclose(fdeb);
	getchar();*//*

						for (integer iscan_par=nd.b0.iright_start; iscan_par<=nd.b0.iright_finish; iscan_par++) {
							//if (ju[iscan_par+1]<ju02start+1) {
								//printf("min perebor 100: ju0start=%d, ju0=%d, ju_no_preobr=%d, ju_preobr=%d\n",ju02start, ju0, ju[iscan_par+1],ju[iscan_par+1]-ju02start+ju0);
								//getchar();
							//}
							//if (ju[iscan_par+1]>ju0gl-100) {
								//printf("max perebor 100: ju0start=%d, ju0=%d, ju_no_preobr=%d, ju_preobr=%d\n",ju02start, ju0, ju[iscan_par+1],ju[iscan_par+1]-ju02start+ju0);
								//getchar();
						//	}
							ju[iscan_par+1]=ju[iscan_par+1]-ju02start+ju0;
							jlu[iscan_par+1]-=ju02start-ju0;
						}
						//ju[nd.b0.iright_finish+2]=ju[nd.b0.iright_finish+2]-ju02start+ju0;  // обрыв последней строки. // это неверно 
						//ju[nd.b0.iright_finish+1]=0; // Вот обрыв соответствующий даным из правилього файла.
						//jlu[nd.b0.iright_finish+2]=0; // дополнение к правильному обрыву второй части.
						jlu[nd.b0.iright_finish+2]=ju[nd.b0.iright_finish+1]; // Вот так наверно совсем правильно.
						jw[nd.b0.iright_finish+2]=jw[nd.b0.iright_finish+1];
						
						for (integer i87=ju02start; i87<=ju0gl; i87++) {
							//alu[ju0]=alu[i87];
							//jlu[ju0]=jlu[i87];
							//levs[ju0]=levs[i87];

						integer i1=ju0;	
		doublereal abub=alu[i1];
		alu[i1]=alu[i87];
		alu[i87]=abub;
		
		integer ibuf=jlu[i1];
		jlu[i1]=jlu[i87];
		jlu[i87]=ibuf;
							
		ibuf=levs[i1];
		levs[i1]=levs[i87];
		levs[i87]=ibuf;


							ju0++;
						}
						*/
						/*
						// debug
	FILE *fdeb=NULL;
					errno_t err87;
	err87 = fopen_s( &fdeb, "lognow.txt", "w");
	//for (integer i54=0; i54<n+2; i54++) {
		//fprintf(fdeb,"%e \n",w[i54]);
//	}
	// Полный лог всех переменных.
	
	
	for (integer i54=0; i54<iwk+1; i54++) {
		if (i54<n+1) {
			fprintf(fdeb,"%e %d %d %d %e %d\n",alu[i54],jlu[i54],levs[i54], ju[i54], w[i54], jw[i54]);
		}
		else if (i54<3*n+1) {
			fprintf(fdeb,"%e %d %d NO NO %d\n",alu[i54],jlu[i54],levs[i54], jw[i54]);
		}
		else {
			fprintf(fdeb,"%e %d %d NO NO NO\n",alu[i54],jlu[i54],levs[i54]);
		}
	}
	fprintf(fdeb,"alu   jlu    levs    ju    w    jw\n");
	fprintf(fdeb,"ileftstart=%d, ileftfinish=%d, irightstart=%d, irightfinish=%d, isepstart=%d, isepfinish=%d\n",nd.b0.ileft_start,nd.b0.ileft_finish,nd.b0.iright_start,nd.b0.iright_finish, nd.b0.iseparate_start,nd.b0.iseparate_finish);
	fprintf(fdeb,"ju0finish1=%d, ju0finish2=%d\n",ideb, ju0);
	printf("hok\n");
	fclose(fdeb);
	getchar();
	*/
						/*
						
						//jlu[ju0]=ju0;
						//ju0--; // убрать в след тесте
						//printf("ju0=%d\n");
						//getchar();
		
					//	printf("alu jlu\n");
					//for (integer i86=debug1-5; i86<debug1+6; i86++) {
						//printf("%e, %d\n",alu[i86],jlu[i86]);
					//}
					//getchar();

					for (integer iscan_par=nd.b0.iseparate_start; iscan_par<=nd.b0.iseparate_finish; iscan_par++) {
						
					iluk_2_serialpart(iscan_par, n2, ju0, n, 
					   a, ja, ia,
					    lfil, alu, jlu, ju, 
	                   levs, iwk, w, jw,ierr_messag);
						

					}
					}*/
					//printf("ju0=%d jluls=%d jlulf=%d jlurs=%d jlurf=%d jluss=%d jlusf=%d\n",ju0,jlu[nd.b0.ileft_start+1],jlu[nd.b0.ileft_finish+1],jlu[nd.b0.iright_start+1],jlu[nd.b0.iright_finish+1],jlu[nd.b0.iseparate_start+1],jlu[nd.b0.iseparate_finish+1]);
						//	getchar();
					/*
					FILE *fdeb=NULL;
					errno_t err87;
	err87 = fopen_s( &fdeb, "lognow.txt", "w");
//	for (integer i54=0; i54<n+2; i54++) {
	//	fprintf(fdeb,"%e \n",w[i54]);
	//}
	for (integer i54=0; i54<iwk+2; i54++) {
		fprintf(fdeb,"%d \n",jlu[i54]);
	}
	fclose(fdeb);
	getchar();*//*
	
	
				}
			}

#else
integer ideb;  */ 
// ju0=n+2;
//ju0=2*n+2;

/* ----------------------------------------------------------------------- */
/*     beginning of main loop. */
/* ----------------------------------------------------------------------- *//*
    i__1 = n;
    for (integer ii = 1; ii <= i__1; ++ii) {

		

		integer iscan_par=ii-1;

	    iluk_2_serialpart(iscan_par, n2, ju0, n, 
					   a, ja, ia,
					    lfil, alu, jlu, ju, 
	                   levs, iwk, w, jw, ierr_messag);
              	if (ierr_messag!=0) {
		            break; // досрочный выход из цикла for.
	            }

				if (ii==nd.b0.ileft_finish) ideb=ju0;
						

    }*/
	/*
	for (integer i87=2*n+2; i87<=ju0; i87++) {
		integer i1=i87-n;
		doublereal abub=alu[i1];
		alu[i1]=alu[i87];
		alu[i87]=abub;
		
		integer ibuf=jlu[i1];
		jlu[i1]=jlu[i87];
		jlu[i87]=ibuf;
							
		ibuf=levs[i1];
		levs[i1]=levs[i87];
		levs[i87]=ibuf;
							
							
						}
						*//*
	// debug
	FILE *fdeb=NULL;
					errno_t err87;
	err87 = fopen_s( &fdeb, "lognow.txt", "w");
	//for (integer i54=0; i54<n+2; i54++) {
		//fprintf(fdeb,"%e \n",w[i54]);
//	}
	// Полный лог всех переменных.
	
	
	for (integer i54=0; i54<iwk+1; i54++) {
		if (i54<n+1) {
			fprintf(fdeb,"%e %d %d %d %e %d\n",alu[i54],jlu[i54],levs[i54], ju[i54], w[i54], jw[i54]);
		}
		else if (i54<3*n+1) {
			fprintf(fdeb,"%e %d %d NO NO %d\n",alu[i54],jlu[i54],levs[i54], jw[i54]);
		}
		else {
			fprintf(fdeb,"%e %d %d NO NO NO\n",alu[i54],jlu[i54],levs[i54]);
		}
	}
	fprintf(fdeb,"alu   jlu    levs    ju    w    jw\n");
	fprintf(fdeb,"ileftstart=%d, ileftfinish=%d, irightstart=%d, irightfinish=%d, isepstart=%d, isepfinish=%d\n",nd.b0.ileft_start,nd.b0.ileft_finish,nd.b0.iright_start,nd.b0.iright_finish, nd.b0.iseparate_start,nd.b0.iseparate_finish);
	fprintf(fdeb,"ju0finish1=%d, ju0finish2=%d\n",ideb, ju0);
	printf("hok\n");
	fclose(fdeb);
	getchar();*//*
	
#endif

	switch (ierr_messag) {
	   case -1: goto L995; break;
	   case -2: goto L996; break;
	   case -5: goto L999; break;
	}
    
	++jw;
    ++w;
    ++ju;
    ++ia;
    ++a;
    ++ja;
    ++alu;
    ++jlu;
    ++levs;

	ierr = 0;
    return 0;



L995:
	++jw;
    ++w;
    ++ju;
    ++ia;
    ++a;
    ++ja;
    ++alu;
    ++jlu;
    ++levs;

    ierr = -1;
    return 0;



L996:
	++jw;
    ++w;
    ++ju;
    ++ia;
    ++a;
    ++ja;
    ++alu;
    ++jlu;
    ++levs;

    ierr = -2;
    return 0;




    ierr = -3;
    return 0;



L998:
	++jw;
    ++w;
    ++ju;
    ++ia;
    ++a;
    ++ja;
    ++alu;
    ++jlu;
    ++levs;

    ierr = -4;
    return 0;



L999:
	++jw;
    ++w;
    ++ju;
    ++ia;
    ++a;
    ++ja;
    ++alu;
    ++jlu;
    ++levs;

    ierr = -5;
    return 0;

} 
*/