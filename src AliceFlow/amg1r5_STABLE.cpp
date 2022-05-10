/* amg1r5new1.f (amg1r6new1.f) -- translated by f2c (version 20100827).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

//typedef long int integer;
// #define integer int64_t
//typedef double doublereal;

//#define doublereal double 
//#define integer int

bool yes_print_amg = true;




//#include "f2c.h"  // этот заголовочный файл просто ненужен.
#include <time.h>
//#include <math.h>

// недостающий функционал
// эта процедура определена выше по коду в mysolverv0_03.c, поэтому здесь её дефениция излишна.
/*integer min(integer ia, integer ib ) 
{
	integer ir;
	if (ia<ib) ir=ia;
	else ir=ib;
	return ir;
} // min
*/

// недостающий функционал 1 октября 2016.
// эта процедура определена выше по коду в mysolverv0_03.c, поэтому здесь её дефениция излишна.
integer my_imin(integer ia, integer ib )
{
    integer ir;
    if (ia<ib) ir=ia;
    else ir=ib;
    return ir;
} // imin


// 1 october 2016 cuda compiller.


integer myi_max(integer ia, integer ib ) 
{
	integer ir;
	if (ia<ib) ir=ib;
	else ir=ia;
	return ir;
} // max


doublereal myr_max(doublereal da, doublereal db ) 
{
	doublereal dr;
	if (da<db) dr=db;
	else dr=da;
	return dr;
} // max



doublereal d_lg10(doublereal * x) {
    // десятичный логарифм вещественого числа.
	return log10(*x);
}

doublereal pow_dd(doublereal * xa, doublereal * xb) {
	// возведение вещественного числа *xa в вещественную степень *xb.
	return pow(*xa,*xb);
}

integer pow_ii(integer *ia, integer *ib) {
	// Возведение целого числа ia в целую степень ib.
	// на основе вещественной функции класса <math.h>
	doublereal ra = static_cast<doublereal>(ia[0]);
	doublereal rb = static_cast<doublereal>(ib[0]);
	integer ii = static_cast<integer>(pow(ra, rb));
	return ii;
}

integer i_sign(integer *ia, integer *ib) {
	// Возвращает abs(ia)*s, где s=+1 если ib больше либо равно нулю.
	// -1 наоборот.
	integer s=-1;
	if (ib[0]>=0) {
		s=1;
	}
	return (abs(ia[0])*s);
}

/* Table of constant values */

integer c__2 = 2;
integer c__4 = 4;
integer c__25 = 25;
integer c__1 = 1;
integer c__0 = 0;
integer c__9 = 9;
integer c__3 = 3;
integer c__10 = 10;

//#if AMG1R6_LABEL==1

/*     Last change:  K    25 Jul 2002   12:55 pm */
/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */

/*     AMG1R6                                        MAIN SUBROUTINE */

/*     RELEASE 1.6, July 2002 */
/* 1.  changed: value of ntrim in pcol */
/* 2.  dimensioning (1) changed to (*) in some subroutines to avoid subscript */
/*     range checks in sparse solvers */

//#endif


/*     Last change:  ERB  22 Aug 2000   10:31 am */
/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */

/*     AMG1R5                                        MAIN SUBROUTINE */

/*     RELEASE 1.5, OCTOBER 1990 */

/*     CHANGES AGAINST VERSION 1.1, JULY 1985: */

/* 1.  A BUG WAS DETECTED WHICH UNDER CERTAIN CIRCUMSTANCES INFLUENCED */
/*     SLIGHTLY THE CONVERGENCE RATE OF AMG1R1. FOR THAT REASON, THE */
/*     FOLLOWING LINE IN SUBROUTINE RESC: */
/*     IW(IMAXW(KC-1)+1) = IA(IMIN(KC)) */
/*     HAS BEEN CHANGED TO: */
/*     IW(IMAXW(KC-1)+1) = IAUX */

/*
1.  БЫЛА ОБНАРУЖЕНА ОШИБКА, КОТОРАЯ ПРИ ОПРЕДЕЛЕННЫХ ОБСТОЯТЕЛЬСТВАХ 
НЕМНОГО ПОВЛИЯЛА НА СКОРОСТЬ КОНВЕРГЕНЦИИ AMG1R1. ПО ЭТОЙ ПРИЧИНЕ 
СЛЕДУЮЩАЯ СТРОКА В ПОДПРОГРАММЕ RESC:
IW(IMAXW(KC-1)+1) = IA(IMIN(KC))
БЫЛ ИЗМЕНЕН НА:
IW(IMAXW(KC-1)+1) = IAUX
*/

/* 2.  A BUG WAS DETECTED IN SUBROUTINE PWINT. UNDER CERTAIN CIRCUM- */
/*     STANCES AN UNDEFINED VARIABLE WAS USED. ALTHOUGH THIS DID NOT */
/*     AFFECT THE NUMERICAL RESULTS, PROBLEMS CAN OCCUR IF CHECKING */
/*     FOR UNDEFINED VARIABLES IS USED. TO FIX THIS ERROR, IN PWINT */
/*     THE LABEL 1000 WAS MOVED TO THE STATEMENT */
/*     IBLCK1 = IMINW(K). */

/* 3.  A PARAMETER LRATIO HAS BEEN INTRODUCED, DENOTING THE RATIO */
/*     OF SPACE OCCUPIED BY A DOUBLE PRECISION REAL VARIABLE AND */
/*     THAT OF AN INTEGER. FOR THE IBM-VERSION LRATIO HAS BEEN SET */
/*     TO 2. CHANGE THIS VALUE IF NECESSARY. (IF, FOR EXAMPLE, YOU */
/*     WANT TO CHANGE THE DOUBLE PRECISION VECTORS TO SINGLE PRE- */
/*     CISION, LRATIO HAS TO BE SET TO 1. IN THE YALE SMP - ROUTINE */
/*     NDRV THERE IS A PARAMETER LRATIO, TOO. */

/*
3.  ВВЕДЕН ПАРАМЕТР LRATIO, ОБОЗНАЧАЮЩИЙ ОТНОШЕНИЕ
ПРОСТРАНСТВА ЗАНЯТОГО ПЕРЕМЕННОЙ ДВОЙНОЙ ТОЧНОСТИ РЕАЛЬНОЙ И
 ЦЕЛОГО ЧИСЛА. ДЛЯ IBM ЗАДАНО СООТНОШЕНИЕ LRATIO
РАВНОЕ 2. ПРИ НЕОБХОДИМОСТИ ИЗМЕНИТЕ ЭТО ЗНАЧЕНИЕ. (ЕСЛИ, НАПРИМЕР, ВЫ
ЧТОБЫ ИЗМЕНИТЬ ВЕКТОРЫ ДВОЙНОЙ ТОЧНОСТИ НА ОДИНАРНУЮ ТОЧНОСТЬ, 
НЕОБХОДИМО УСТАНОВИТЬ КОЭФФИЦИЕНТ РАВНЫМ 1. В ЙЕЛЬСКОМ SMP-РЕЖИМЕ
NDRV СУЩЕСТВУЕТ ТАКЖЕ ПАРАМЕТР LRATIO.
*/

/* 4.  TYPE DECLARATIONS REAL*4 AND REAL*8 HAVE BEEN CHANGED TO THE */
/*     STANDARD-CONFORMING KEYWORDS REAL AND DOUBLE PRECISION, RESPEC- */
/*     TIVELY. */

/* 5.  CALLS TO THE FOLLOWING INTRINSIC FUNCTIONS HAVE BEEN REPLACED BY */
/*     CALLS USING GENERIC NAMES: DSQRT, MIN0, MAX0, IABS, DABS, FLOAT, */
/*     DFLOAT, DMAX1, ISIGN, IDINT, DLOG10. */

/* 6.  A SAVE STATEMENT HAS BEEN INSERTED IN ALL SUBROUTINES. */

/* 7.  EXTERNAL DECLARATION STATEMENTS HAVE BEEN INSERTED IN ALL SUB- */
/*     ROUTINES FOR ALL EXTERNAL REFERENCES. */

/* ----------------------------------------------------------------------- */

/*     CHANGE AGAINST VERSION 1.3, APRIL 1986: */

/* 1.  A BUG IN SUBROUTINE CHECK HAS BEEN REMOVED. IF THE ORIGINAL MATRIX */
/*     WAS STORED IN AN UNSYMMETRIC WAY, THE SYMMETRIZATION BY AMG1R3 */
/*     COULD FAIL UNDER CERTAIN CIRCUMSTANCES. FOR A FIX, THE FOLLOWING */
/*     STATEMENTS IN SUBROUTINE CHECK HAVE BEEN CHANGED: */
/*
1.  ОШИБКА В ПРОВЕРКЕ ПОДПРОГРАММЫ БЫЛА УДАЛЕНА. ЕСЛИ ИСХОДНАЯ МАТРИЦА
БЫЛА СОХРАНЕНА НЕСИММЕТРИЧНЫМ СПОСОБОМ, СИММЕТРИЗАЦИЯ AMG1R3 МОЖЕТ
ПОТЕРПЕТЬ НЕУДАЧУ ПРИ ОПРЕДЕЛЕННЫХ ОБСТОЯТЕЛЬСТВАХ. ДЛЯ ИСПРАВЛЕНИЯ
БЫЛИ ИЗМЕНЕНЫ СЛЕДУЮЩИЕ ИНСТРУКЦИИ В ФУНКЦИИ CHECK:
*/

/*     DO 450 J=IA(I)+1,IA(I+1)-1 WAS CHANGED TO */
/*     DO 450 J=IA(I)+1,ICG(I)-1 */

/*     DO 430 J1=IA(I1)+1,IA(I1+1)-1 WAS CHANGED TO */
/*     DO 430 J1=IA(I1)+1,ICG(I1)-1 */

/*     DO 550 J=IA(I)+1,IA(I+1)-1 WAS CHANGED TO */
/*     DO 550 J=IA(I)+1,ICG(I)-1 */

/*     DO 530 J1=IA(I1)+1,IA(I1+1)-1 WAS CHANGED TO */
/*     DO 530 J1=IA(I1)+1,ICG(I1)-1 */

/* 2.  THE EXPLANATORY PART IN SUBROUTINE AMG1R5 HAS BEEN ENLARGED TO */
/*     AVOID MISUNDERSTANDINGS IN THE DEFINITION OF THE ARGUMENT LIST. */
/*
2.  ПОЯСНИТЕЛЬНАЯ ЧАСТЬ ПОДПРОГРАММЫ AMG1R5 БЫЛА РАСШИРЕНА
ВО ИЗБЕЖАНИЕ НЕДОРАЗУМЕНИЙ В ОПРЕДЕЛЕНИИ СПИСКА АРГУМЕНТОВ.
*/

/* ----------------------------------------------------------------------- */

/*     CHANGE AGAINST VERSION 1.4, OCTOBER, 1990 (BY JOHN W. RUGE) */

/* 1.  A BUG IN SUBROUTINE CHECK HAS BEEN REMOVED. IF THE ORIGINAL MATRIX */
/*     WAS STORED IN AN UNSYMMETRIC WAY, THE SYMMETRIZATION BY AMG1R3 */
/*     COULD STILL FAIL UNDER CERTAIN CIRCUMSTANCES, AND WAS NOT FIXED */
/*     IN THE PREVIOUS VERSION. IN ADDITION, THE ROUTINE WAS CHANGED */
/*     IN ORDER TO AVOID SOME UNNECESSARY ROW SEARCHES FOR TRANSOSE */
/*     ENTRIES. */
/*
1.  ОШИБКА В ПРОВЕРКЕ ПОДПРОГРАММЫ БЫЛА УДАЛЕНА. ЕСЛИ ИСХОДНАЯ МАТРИЦА
БЫЛА СОХРАНЕНА НЕСИММЕТРИЧНЫМ СПОСОБОМ, СИММЕТРИЗАЦИЯ ПО AMG1R3 ВСЕ ЕЩЕ
МОГЛА ПОТЕРПЕТЬ НЕУДАЧУ ПРИ ОПРЕДЕЛЕННЫХ ОБСТОЯТЕЛЬСТВАХ И НЕ БЫЛА
ЗАФИКСИРОВАНА В ПРЕДЫДУЩЕЙ ВЕРСИИ. КРОМЕ ТОГО, ПРОЦЕДУРА БЫЛА ИЗМЕНЕНА, 
ЧТОБЫ ИЗБЕЖАТЬ НЕКОТОРЫХ НЕНУЖНЫХ ПОИСКОВ СТРОК ДЛЯ ТРАНСПОНИРОВАНИЯ ЗАПИСЕЙ.
*/

/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */

/* Subroutine */ integer amg1r5_(doublereal *a, integer *ia, integer *ja, 
	doublereal *u, doublereal *f, integer *ig, integer *nda, integer *
	ndia, integer *ndja, integer *ndu, integer *ndf, integer *ndig, 
	integer *nnu, integer *matrix, integer *iswtch, integer *iout, 
	integer *iprint, integer *levelx, integer *ifirst, integer *ncyc, 
	doublereal *eps, integer *madapt, integer *nrd, integer *nsolco, 
	integer *nru, doublereal *ecg1, doublereal *ecg2, doublereal *ewt2, 
	integer *nwt, integer *ntr, integer *ierr)
{
    /* Format strings */
   

    /* Builtin functions */
    

    /* Local variables */
    integer mda=0, mdf=0, mdu=0;
    doublereal res=0.0;
    integer ium=0, iup=0;
    doublereal res0=0.0;
    extern /* Subroutine */ integer idec_(integer *, integer *, integer *, 
	    integer *);
	integer mdia = 0, mdja = 0, mdig = 0, iarr[25] = { 0 };
    //static real time[20];
	unsigned int time[20] = {0};
	integer imin[25] = { 0 }, imax[25] = { 0 };
	doublereal resi[25] = { 0.0 };
	integer kout = 0, ncyc0 = 0, irow0 = 0, ndicg = 0, icgst = 0, iminw[25] = { 0 }, imaxw[25] = { 0 };
    extern /* Subroutine */ integer first_(integer *, doublereal *, integer *, 
	    integer *, integer *, integer *), solve_(integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, doublereal 
	    *, doublereal *, doublereal *, integer *, integer *, integer *, 
	    doublereal *, integer *, integer *, integer *, integer *, integer 
	    *, integer *, integer *, integer *, /*real*/ unsigned int *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, doublereal 
	    *, doublereal *, doublereal *), setup_(integer *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, /*real*/ unsigned int *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *);
	integer ndigit = 0, levels = 0, kevelx = 0, nstcol[25] = { 0 }, kswtch=0;
    extern /* Subroutine */ integer wrkcnt_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, /*real*/ unsigned int *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, doublereal *, doublereal *);

    /* Fortran I/O blocks */
  



/*         ----------------------------------------------- */
/*         | AMG-MODULE FOR SOLVING LINEAR SYSTEMS L*U=F | */
/*         ----------------------------------------------- */

/*         -------------------------------------------------------------- */

/*     ASSUMPTIONS ON L: */

/*         THE PROGRAM REQUIRES: */ // Требования к входным данным:

/*             - DIAGONAL ENTRIES ARE ALWAYS POSITIVE (ON ALL GRIDS); */
/*             - L IS A SQUARE MATRIX WHICH IS EITHER REGULAR OR SINGULAR */
/*               WITH ROWSUMS=0. */
	// Матрица L это квадратная матрица с всегда положительными диагональными элементами с нулевой суммой коэффициентов в каждой строке (положительная определённость).

/*         FOR THEORETICAL REASONS THE FOLLOWING SHOULD HOLD: */

/*             - L POSITIVE DEFINITE (OR SEMI-DEFINITE WITH ROWSUM=0) */
/*             - L "ESSENTIALLY" POSITIVE TYPE, I.E., */

/*                  -- DIAGONAL ENTRIES MUST BE > 0 ; */
/*                  -- MOST OF THE OFF-DIAGONAL ENTRIES <= 0 ; */
/*                  -- ROWSUMS SHOULD BE >= 0 . */

// Матрица L положительно определённая: диагональные элементы  строго больше нуля, большинство внедиагональных элементов <= 0.
// Сумма коэффициентов в строке больше либо равна нулю - диагональное преобладание.


/*     THE USER HAS TO PROVIDE THE MATRIX L, THE RIGHT HAND SIDE F AND */
/*     CERTAIN POINTER VECTORS IA AND JA. */

// Пользователь задаёт матрицу коэффициентов L, правую часть F, а также информацию о связях между элементами матрицы в специальных столбцах IA и JA.
// IA - позиция первого элемента в строке. JA - номер столбца для каждого элемента, первым записывается диагональный элемент.
// В фортране нумерация начинается с единицы.

/*         -------------------------------------------------------------- */

/*     STORAGE OF L: */ // Требования к хранению матрицы L.

/*         THE NON-ZERO ENTRIES OF THE MATRIX L ARE STORED IN */
/*         "COMPRESSED" SKY-LINE FASHION IN A 1-D VECTOR A, I.E., ROW */
/*         AFTER ROW, EACH ROW STARTING WITH ITS DIAGONAL ELEMENT. THE */
/*         OTHER NON-ZERO ROW ENTRIES FOLLOW THEIR DIAGONAL ENTRY IN ANY */
/*         ORDER. */


/*         IN ORDER TO IDENTIFY EACH ELEMENT IN A, THE USER HAS TO */
/*         PROVIDE TWO POINTER ARRAYS IA AND JA. IF NNU DENOTES THE TOTAL */
/*         NUMBER OF UNKNOWNS, THE NON-ZERO ENTRIES OF ANY ROW I OF L */
/*         (1.LE.I.LE.NNU) ARE STORED IN A(J) WHERE THE RANGE OF J */
/*         IS GIVEN BY */

/*                     IA(I) .LE. J .LE. IA(I+1)-1. */

/*         THUS, IA(I) POINTS TO THE POSITION OF THE DIAGONAL ENTRY OF */
/*         ROW I WITHIN THE VECTOR A. IN PARTICULAR, */

/*                     IA(1) = 1 ,  IA(NNU+1) = 1 + NNA */

/*         WHERE NNA DENOTES THE TOTAL NUMBER OF MATRIX ENTRIES STORED. */
/*         THE POINTER VECTOR JA HAS TO BE DEFINED SUCH THAT */
/*         ANY ENTRY A(J) CORRESPONDS TO THE UNKNOWN U(JA(J)), I.E., */
/*         JA(J) POINTS TO THE COLUMN INDEX OF A(J). */
/*         IN PARTICULAR, A(IA(I)) IS THE DIAGONAL ENTRY OF ROW I */
/*         AND CORRESPONDS TO THE UNKNOWN U(I): JA(IA(I))=I. */

/*         IN THIS TERMINOLOGY, THE I-TH EQUATION READS AS FOLLOWS */
/*         (FOR ANY I WITH  1.LE.I.LE.NNU): */

/*                  F(I) =        SUM      A(J) * U(JA(J)) */
/*                           J1.LE.J.LE.J2 */

/*         WHERE F(I) DENOTES THE I-TH COMPONENT OF THE RIGHT HAND */
/*         SIDE AND */

/*                     J1 = IA(I) ,  J2 = IA(I+1)-1. */

/*         NOTES: THE ENTRY IA(NNU+1) HAS TO TOCHKA TO THE FIRST FREE */
/*                ENTRY IN VECTORS A AND JA, RESPECTIVELY. OTHERWISE, */
/*                AMG CANNOT KNOW THE LENGTH OF THE LAST MATRIX ROW. */

/*                THE INPUT VECTORS A, IA AND JA ARE CHANGED BY AMG1R5. */
/*                SO, AFTER RETURN FROM AMG1R5, THE PACKAGE MUST NOT */
/*                BE CALLED A SECOND TIME WITHOUT HAVING NEWLY DEFINED */
/*                THE INPUT VECTORS AND USING ISWTCH=4. OTHERWISE, THE */
/*                SETUP PHASE WILL FAIL. */
/*                  ON THE OTHER HAND, RUNNING AMG A SECOND TIME ON THE */
/*                SAME INPUT DATA WITH ISWTCH=4 HAS NO SENSE, BECAUSE */
/*                THE RESULTS OF THE FIRST SETUP PHASE ARE STILL STORED */
/*                AND THUS THIS PHASE CAN BE SKIPPED IN A SECOND CALL. */
/*                IN ORDER TO DO THIS, SET ISWTCH TO 1, 2 OR 3. */

/* ----------------------------------------------------------------------- */

/*         THE FORM OF THE CALLING PROGRAM HAS TO BE AS FOLLOWS: */

/*               PROGRAM DRIVER */
/*         C */
/*               DOUBLE PRECISION A(#NDA),U(#NDU),F(#NDF) */
/*               INTEGER IA(#NDIA),JA(#NDJA),IG(#NDIG) */
/*         C */
/*               NDA  = #NDA */
/*               NDU  = #NDU */
/*               NDF  = #NDF */
/*               NDIA = #NDIA */
/*               NDJA = #NDJA */
/*               NDIG = #NDIG */
/*         C */
/*         C     SET UP A, F, IA, JA AND SPECIFY NECESSARY PARAMETERS */
/*         C */
/*               .... */
/*               .... */
/*         C */
/*               CALL AMG1R5(A,IA,JA,U,F,IG, */
/*        +                  NDA,NDIA,NDJA,NDU,NDF,NDIG,NNU,MATRIX, */
/*        +                  ISWTCH,IOUT,IPRINT, */
/*        +                LEVELX,IFIRST,NCYC,EPS,MADAPT,NRD,NSOLCO,NRU, */
/*        +                  ECG1,ECG2,EWT2,NWT,NTR, */
/*        +                  IERR) */
/*         C */
/*               .... */
/*               .... */
/*         C */
//#if AMG1R6_LABEL==1
/*               CALL USTOP(' ') */
//#else
/*               STOP */
//#endif
/*               END */

/* ----------------------------------------------------------------------- */

/*     INPUT VIA ARRAYS (SEE ABOVE): */

/*     A        -   MATRIX L */

/*     IA       -   POINTER VECTOR */

/*     JA       -   POINTER VECTOR */

/*     U        -   FIRST APPROXIMATION TO SOLUTION */

/*     F        -   RIGHT HAND SIDE */


/* ----------------------------------------------------------------------- */


/*     SCALAR INPUT PARAMETERS OF AMG1R5: */

/*     THE INPUT PARAMETERS OF AMG1R5 IN THE LIST BELOW ARE ARRANGED */
/*     ACCORDING TO THEIR IMPORTANCE TO THE GENERAL USER. THE PARAMETERS */
/*     PRECEEDED BY A * MUST BE SPECIFIED EXPLICITELY. ALL THE OTHER */
/*     PARAMETERS ARE SET TO STANDARD VALUES IF ZERO ON INPUT. */

/*     THERE ARE FOUR CLASSES OF INPUT PARAMETERS WITH DECREASING PRI- */
/*     ORITY: */

/*     1. PARAMETERS DESCRIBING THE USER-DEFINED PROBLEM AND DIMENSIONING */
/*        OF VECTORS IN THE CALLING PROGRAM */

/*     2. PARAMETERS SPECIFYING SOME GENERAL ALGORITHMIC ALTERNATIVES AND */
/*        THE AMOUNT OF OUTPUT DURING SOLUTION */

/*     3. PARAMETERS CONTROLLING THE MULTIGRID CYCLING DURING THE SOLU- */
/*        TION PHASE */

/*     4. PARAMETERS CONTROLLING THE CREATION OF COARSER GRIDS AND INTER- */
/*        POLATION FORMULAS. */

/*     ONLY THE CLASS 1 - PARAMETERS MUST BE SPECIFIED EXPLICITELY BY */
/*     THE USER. CLASS 2 - PARAMETERS CONTROL THE GENERAL PERFORMANCE OF */
/*     AMG1R5. CHANGING THEM DOESN'T REQUIRE UNDERSTANDING THE AMG - */
/*     ALGORITHM. SPECIFYING NON-STANDARD-VALUES FOR CLASS 3 - PARAMETERS */
/*     PRESUPPOSES A GENERAL KNOWLEDGE OF MULTIGRID METHODS, WHEREAS THE */
/*     FUNCTION OF CLASS 4 - PARAMETERS IS ONLY UNDERSTANDABLE AFTER */
/*     STUDYING THE AMG-ALGORITHM IN DETAIL. FORTUNATELY IN MOST CASES */
/*     THE CHOICE OF CLASS 3 AND 4 - PARAMETERS ISN'T CRITICAL AND USING */
/*     THE AMG1R5 - SUPPLIED STANDARD VALUES SHOULD GIVE SATISFACTORY */
/*     RESULTS. */

/*         -------------------------------------------------------------- */

/*     CLASS 1 - PARAMETERS: */

/*  *  NDA      -   DIMENSIONING OF VECTOR A IN CALLING PROGRAM */

/*  *  NDIA     -   DIMENSIONING OF VECTOR IA IN CALLING PROGRAM */

/*  *  NDJA     -   DIMENSIONING OF VECTOR JA IN CALLING PROGRAM */

/*  *  NDU      -   DIMENSIONING OF VECTOR U IN CALLING PROGRAM */

/*  *  NDF      -   DIMENSIONING OF VECTOR F IN CALLING PROGRAM */

/*  *  NDIG     -   DIMENSIONING OF VECTOR IG IN CALLING PROGRAM */

/*  *  NNU      -   NUMBER OF UNKNOWNS */

/*  *  MATRIX   -   INTEGER VALUE CONTAINING INFO ABOUT THE MATRIX L. */

/*                  1ST DIGIT OF MATRIX  --  ISYM: */
/*                    =1: L IS SYMMETRIC; */
/*                    =2: L IS NOT SYMMETRIC. */

/*                  2ND DIGIT OF MATRIX  --  IROW0: */
/*                    =1: L HAS ROWSUM ZERO; */
/*                    =2: L DOES NOT HAVE ROWSUM ZERO. */

/*         -------------------------------------------------------------- */

/*     CLASS 2 - PARAMETERS: */

/*     ISWTCH   -   PARAMETER CONTROLLING WHICH MODULES OF AMG1R5 ARE TO */
/*                  BE USED. */
/*                    =1:   CALL FOR -----, -----, -----, WRKCNT. */
/*                    =2:   CALL FOR -----, -----, SOLVE, WRKCNT. */
/*                    =3:   CALL FOR -----, FIRST, SOLVE, WRKCNT. */
/*                    =4:   CALL FOR SETUP, FIRST, SOLVE, WRKCNT. */
/*                  SETUP DEFINES THE OPERATORS NEEDED IN THE SOLUTION */
/*                         PHASE. */
/*                  FIRST INITIALIZES THE SOLUTION VECTOR (SEE PARAMETER */
/*                         IFIRST). */
/*                  SOLVE COMPUTES THE SOLUTION BY AMG CYCLING (SEE */
/*                         PARAMETER NCYC). */
/*                  WRKCNT PROVIDES THE USER WITH INFORMATION ABOUT */
/*                         RESIDUALS, STORAGE REQUIREMENTS AND CP-TIMES */
/*                         (SEE PARAMETER IOUT). */
/*                  IF AMG1R5 IS CALLED THE FIRST TIME, ISWTCH HAS TO */
/*                  BE =4. INDEPENDENT OF ISWTCH, SINGLE MODULES CAN BE */
/*                  BYPASSED BY A PROPER CHOICE OF THE CORRESPONDING */
/*                  PARAMETER. */

/*     IOUT     -   PARAMETER CONTROLLING THE AMOUNT OF OUTPUT DURING */
/*                  SOLUTION PHASE: */

/*                  1ST DIGIT: NOT USED; HAS TO BE NON-ZERO. */

/*                  2ND DIGIT: */
/*                    =0: NO OUTPUT (EXCEPT FOR MESSAGES) */
/*                    =1: RESIDUAL BEFORE AND AFTER SOLUTION PROCESS */
/*                    =2: ADD.: STATISTICS ON CP-TIMES AND STORAGE REQUI- */
/*                        REMENTS */
/*                    =3: ADD.: RESIDUAL AFTER EACH AMG-CYCLE */

/*     IPRINT   -   PARAMETER SPECIFYING THE FORTRAN UNIT NUMBERS FOR */
/*                  OUTPUT: */

/*                  1ST DIGIT: NOT USED; HAS TO BE NON-ZERO */

/*                  2ND AND 3RD DIGIT  --  IUP: UNIT NUMBER FOR RESULTS */

/*                  4TH AND 5TH DIGIT  --  IUM: UNIT NUMBER FOR MESSAGES */

/*         -------------------------------------------------------------- */

/*     CLASS 3 - PARAMETERS: */

/*     LEVELX   -   MAXIMUM NUMBER OF MG-LEVELS TO BE CREATED (>=1). */

/*     IFIRST   -   PARAMETER FOR FIRST APPROXIMATION. */

/*                  1ST DIGIT OF IFIRST: NOT USED; HAS TO BE NON-ZERO. */

/*                  2ND DIGIT OF IFIRST  --  ITYPU: */
/*                    =0: NO SETTING OF FIRST APPROXIMATION, */
/*                    =1: FIRST APPROXIMATION CONSTANT TO ZERO, */
/*                    =2: FIRST APPROXIMATION CONSTANT TO ONE, */
/*                    =3: FIRST APPROXIMATION IS RANDOM FUNCTION WITH */
/*                        THE CONCRETE RANDOM SEQUENCE BEING DETERMINED */
/*                        BY THE FOLLWING DIGITS. */

/*                  REST OF IFIRST  --  RNDU: */
/*                    DETERMINES THE CONCRETE RANDOM SEQUENCE USED IN */
/*                    THE CASE ITYPU=3. (IFIRST=13 IS EQUIVALENT TO */
/*                    IFIRST=1372815) */

/*     NCYC     -   INTEGER PARAMETER DESCRIBING THE TYPE OF CYCLE TO BE */
/*                  USED AND THE NUMBER OF CYCLES TO BE PERFORMED. */

/*                  1ST DIGIT OF NCYC  --  IGAM: */
/*                    =1: V -CYCLE, */
/*                    =2: V*-CYCLE, */
/*                    =3: F -CYCLE, */
/*                    =4: W -CYCLE. */
/*                  IF NCYC IS NEGATIV, THEN THE APPROXIMATION OF THE */
/*                  PROBLEM ON THE SECOND FINEST GRID IS COMPUTED BY */
/*                  IGAM V-CYCLES ON THAT PARTICULAR GRID. */

/*                  2ND DIGIT OF NCYC  --  ICGR: */
/*                    =0: NO CONJUGATE GRADIENT, */
/*                    =1: CONJUGATE GRADIENT (ONLY FIRST STEP OF CG), */
/*                    =2: CONJUGATE GRADIENT (FULL CG). */

/*                  3RD DIGIT OF NCYC  --  ICONV: */
/*                    CONVERGENCE CRITERION FOR THE USER-DEFINED PROBLEM */
/*                    (FINEST GRID): */
/*                    =1: PERFORM A FIXED NUMBER OF CYCLES AS GIVEN BY */
/*                        NCYCLE (SEE BELOW) */
/*                    =2: STOP, IF  ||RES|| < EPS */
/*                    =3: STOP, IF  ||RES|| < EPS * |F| */
/*                    =4: STOP, IF  ||RES|| < EPS * |U| * |DIAG| */
/*                    WITH ||RES|| = L2-NORM OF RESIDUAL, */
/*                           EPS     (SEE INPUT PARAMETER EPS) */
/*                           |F|   = SUPREMUM NORM OF RIGHT HAND SIDE */
/*                           |U|   = SUPREMUM NORM OF SOLUTION */
/*                         |DIAG|  = MAXIMAL DIAGONAL ENTRY IN MATRIX L */
/*                    NOTE THAT IN ANY CASE THE SOLUTION PROCESS STOPS */
/*                    AFTER AT MOST NCYCLE CYCLES. */

/*                  REST OF NCYC  --  NCYCLE: */
/*                    MAXIMAL NUMBER OF CYCLES TO BE PERFORMED (>0) OR */
/*                    NCYCLE=0: NO CYCLING. */

/*     EPS      -   CONVERGENCE CRITERION FOR SOLUTION PROCESS: (SEE */
/*                  PARAMETER NCYC). NOTE THAT NO MORE THAN NCYCLE CYCLES */
/*                  ARE PERFORMED, REGARDLESS OF EPS. */

/*     MADAPT   -   INTEGER VALUE SPECIFYING THE CHOICE OF COARSEST */
/*                  GRID IN CYCLING: */

/*                  1ST DIGIT OF MADAPT  --  MSEL: */
/*                    =1: IN CYCLING, ALL GRIDS CONSTRUCTED IN THE SETUP */
/*                        PHASE ARE USED WITHOUT CHECK. */
/*                    =2: THE NUMBER OF GRIDS IS AUTOMATICALLY REDUCED */
/*                        IF THE CONVERGENCE FACTOR ON THE COARSER GRIDS */
/*                        IS FOUND TO BE LARGER THAN A GIVEN VALUE FAC */
/*                        (SEE BELOW). */

/*                  REST OF MADAPT  --  FAC */
/*                        THE REST OF MADAPT DEFINES THE FRACTIONAL PART */
/*                        OF A REAL NUMBER FAC BETWEEN 0.1 AND 0.99, E.G. */
/*                        MADAPT=258 MEANS MSEL=2 AND FAC=0.58. IF MADAPT */
/*                        CONSISTS OF ONLY ONE DIGIT, FAC IS SET TO 0.7 */
/*                        BY DEFAULT. */


/*     NRD      -   PARAMETER DESCRIBING RELAXATION (DOWNWARDS): */

/*                  1ST DIGIT OF NRD: NOT USED; HAS TO BE NON-ZERO. */

/*                  2ND DIGIT OF NRD  --  NRDX: */
/*                    ACTUAL NUMBER OF SMOOTHING STEPS TO BE PERFORMED */
/*                    THE TYPE OF WHICH IS GIVEN BY THE FOLLOWING DIGITS */

/*                  FOLLOWING DIGITS  --  ARRAY NRDTYP: */
/*                    =1: RELAXATION OVER THE F-POINTS ONLY */
/*                    =2: FULL GS SWEEP */
/*                    =3: RELAXATION OVER THE C-POINTS ONLY */
/*                    =4: FULL MORE COLOR SWEEP, HIGHEST COLOR FIRST */

/*     NSOLCO   -   PARAMETER CONTROLLING THE SOLUTION ON COARSEST GRID: */

/*                  1ST DIGIT  --  NSC: */
/*                    =1: GAUSS-SEIDEL METHOD */
/*                    =2: DIRECT SOLVER (YALE SMP) */

/*                  REST OF NSOLCO  --  NRCX: (ONLY IF NSC=1) */
/*                  NUMBER OF GS SWEEPS ON COARSEST GRID (>=0). */
/*                  IF NRCX=0, THEN AS MANY GS SWEEPS ARE PERFORMED */
/*                  AS ARE NEEDED TO REDUCE THE RESIDUAL BY TWO ORDERS */
/*                  OF MAGNITUDE. (MAXIMAL 100 RELAXATION SWEEPS) */

/*     NRU      -   PARAMETER FOR RELAXATION (UPWARDS), ANALOGOUS TO NRD. */

/*         -------------------------------------------------------------- */

/*     CLASS 4 - PARAMETERS: */

/*     ECG1,ECG2-   REAL PARAMETERS AFFECTING THE CREATION OF COARSER */
/*     EWT2     -   GRIDS AND/OR THE DEFINITION OF THE INTERPOLATION. */
/*                  THE CHOICE OF THESE PARAMETERS DEPENDS ON */
/*                  THE ACTUAL AMG VERSION (SEE SUBROUTINE CRSNG) */

/*     NWT      -   INTEGER PARAMETER AFFECTING THE CREATION OF COARSER */
/*                  GRIDS AND/OR THE DEFINITION OF THE INTERPOLATION. */
/*                  THE CHOICE OF THIS PARAMETER DEPENDS ON */
/*                  THE ACTUAL AMG VERSION (SEE SUBROUTINE CRSNG) */

/*     NTR      -   PARAMETER CONTROLLING COARSE-GRID OPERATOR TRUNCATION */
/*                    =0: PAIRS OF ZEROES ARE REMOVED FROM COARSE GRID */
/*                        OPERATORS */
/*                    =1: NO COARSE-GRID OPERATOR TRUNCATION */


/* ----------------------------------------------------------------------- */

/*     OUTPUT: */

/*     U        -   CONTAINS THE COMPUTED SOLUTION */


/*     IERR     -   ERROR PARAMETER: */

/*                    >0: FATAL ERROR (ABNORMAL TERMINATION OF AMG1R5) */
/*                    <0: NON-FATAL ERROR (EXECUTION OF AMG1R5 CONTINUES) */

/*                  ERROR CODES IN DETAIL: */

/*                  1. DIMENSIONING TOO SMALL FOR VECTOR */
/*                        A      (IERR = 1) */
/*                        IA     (IERR = 2) */
/*                        JA     (IERR = 3) */
/*                        U      (IERR = 4) */
/*                        F      (IERR = 5) */
/*                        IG     (IERR = 6) */

/*                     NO YALE-SMP BECAUSE OF STORAGE (NDA TOO SMALL): */
/*                               (IERR = -1) */
/*                     NO YALE-SMP BECAUSE OF STORAGE (NDJA TOO SMALL): */
/*                               (IERR = -3) */
/*                     NO CG BECAUSE OF STORAGE (NDU TOO SMALL): */
/*                               (IERR = -4) */
/*                     NO SPACE FOR TRANSPOSE OF INTERPOLATION (NDA OR */
/*                                                     NDJA TOO SMALL): */
/*                               (IERR = -1) */

/*                  2. INPUT DATA ERRONEOUS: */

/*                     A-ENTRY MISSING, ISYM = 1:           (IERR = -11) */
/*                     PARAMETER MATRIX MAY BE ERRONEOUS:   (IERR = -12) */
/*                     DIAGONAL ELEMENT NOT STORED FIRST:   (IERR =  13) */
/*                     DIAGONAL ELEMENT NOT POSITIV:        (IERR =  14) */
/*                     POINTER IA ERRONEOUS:                (IERR =  15) */
/*                     POINTER JA ERRONEOUS:                (IERR =  16) */
/*                     PARAMETER ISWTCH ERRONEOUS:          (IERR =  17) */
/*                     PARAMETER LEVELX ERRONEOUS:          (IERR =  18) */

/*                  3. ERRORS OF THE AMG1R5-SYSTEM (SHOULD NOT OCCUR): */

/*                     TRANSPOSE A-ENTRY MISSING:           (IERR =  21) */
/*                     INTERPOLATION ENTRY MISSING:         (IERR =  22) */

/*                  4. ALGORITHMIC ERRORS: */

/*                     CG-CORRECTION NOT DEFINED:           (IERR =  31) */
/*                     NO YALE-SMP BECAUSE OF ERROR IN */
/*                     FACTORIZATION:                       (IERR = -32) */

/* ----------------------------------------------------------------------- */

/*     WORK SPACE: */

/*     THE INTEGER VECTOR IG HAS TO BE PASSED TO AMG1R5 AS WORK SPACE. */

/* ----------------------------------------------------------------------- */

/*     DIMENSIONING OF INPUT VECTORS AND WORK SPACE: */

/*     IT'S IMPOSSIBLE TO TELL IN ADVANCE THE EXACT STORAGE REQUIREMENTS */
/*     OF AMG. THUS, THE FOLLOWING FORMULAS GIVE ONLY REASONABLE GUESSES */
/*     FOR THE VECTOR LENGTHS WHICH HAVE TO BE DECLARED IN THE CALLING */
/*     PROGRAM. IN THESE FORMULAS NNA DENOTES THE NUMBER OF NON-ZERO */
/*     ENTRIES IN THE INPUT-MATRIX L AND NNU IS THE NUMBER OF UNKNOWNS. */

/*     VECTOR         NEEDED LENGTH (GUESS) */
/*       A               3*NNA + 5*NNU */
/*       JA              3*NNA + 5*NNU */
/*       IA              2.2*NNU */
/*       U               2.2*NNU */
/*       F               2.2*NNU */
/*       IG              5.4*NNU */

/* ----------------------------------------------------------------------- */


/*     STANDARD CHOICES OF PARAMETERS (AS FAR AS MEANINGFUL): */

/*          ISWTCH = 4 */
/*          IOUT   = 12 */
/*          IPRINT = 10606 */

/*          LEVELX = 25 */
/*          IFIRST = 13 */
/*          NCYC   = 10110 */
/*          EPS    = 1.D-12 */
/*          MADAPT = 27 */
/*          NRD    = 1131 */
/*          NSOLCO = 110 */
/*          NRU    = 1131 */

/*          ECG1   = 0. */
/*          ECG2   = 0.25 */
/*          EWT2   = 0.35 */
/*          NWT    = 2 */
/*          NTR    = 0 */

/*     IF ANY ONE OF THESE PARAMETERS IS 0 ON INPUT, ITS CORRESPONDING */
/*     STANDARD VALUE IS USED BY AMG1R5. */

/* ----------------------------------------------------------------------- */

/*     PORTABILITY RESTRICTIONS: */

/*     1. ROUTINE CTIME IS MACHINE DEPENDENT AND HAS TO BE ADAPTED TO */
/*        YOUR COMPUTER INSTALLATION OR REPLACED BY A DUMMY ROUTINE. */

/*     2. MOST INPUT PARAMETERS ARE COMPOSED OF SEVERAL DIGITS, THEIR */
/*        SIGNIFICANCE HAVING BEEN DESCRIBED ABOVE. BE SURE NOT TO ENTER */
/*        MORE DIGITS THAN YOUR COMPUTER CAN STORE ON AN INTEGER VARI- */
/*        ABLE. */

/*     3. APART FROM FORTRAN INTRINSIC FUNCTIONS AND SERVICE ROUTINES, */
/*        THERE IS ONLY ONE EXTERNAL REFERENCE TO A PROGRAM NOT CONTAINED */
/*        IN THE AMG1R5 - SYSTEM, I.E. THE LINEAR SYSTEM SOLVER NDRV OF */
/*        THE YALE SPARSE MATRIX PACKAGE. IF YOU HAVN'T ACCESS TO THIS */
/*        PACKAGE, ENTER A DUMMY ROUTINE NDRV AND AVOID CHOOSING NSC=2 */
/*        (SUBPARAMETER OF NSOLCO). THEN NDRV ISN'T CALLED BY AMG1R5. */
/*        IN THIS CASE, HOWEVER, INDEFINITE PROBLEMS WILL NOT BE SOLV- */
/*        ABLE. */
/*          THE YALE SPARSE MATRIX PACKAGE IS FREELY AVAILABLE FOR NON- */
/*        PROFIT PURPOSES. CONTACT THE DEPARTMENT OF COMPUTER SCIENCE, */
/*        YALE UNITVERSITY. */

/*     4. IN AMG1R5 THERE IS THE PARAMETER LRATIO, DENOTING THE RATIO */
/*        OF SPACE OCCUPIED BY A DOUBLE PRECISION REAL VARIABLE AND */
/*        THAT OF AN INTEGER. FOR THE IBM-VERSION LRATIO HAS BEEN SET */
/*        TO 2. CHANGE THIS VALUE IF NECESSARY. (THE SAME HAS TO BE */
/*        DONE WITH THE YALE SMP-ROUTINE NDRV.) */


/* ----------------------------------------------------------------------- */

/*     AUTHORS: */

/*          JOHN RUGE, FORT COLLINS (USA), */
/*              INSTITUTE FOR COMPUTATIONAL STUDIES AT CSU; */

/*          KLAUS STUEBEN, D-5205 ST. AUGUSTIN (W.-GERMANY), */
/*              GESELLSCHAFT FUER MATHEMATIK UND DATENVERARBEITUNG (GMD). */

/*          ROLF HEMPEL, D-5205 ST. AUGUSTIN (W.-GERMANY), */
/*              GESELLSCHAFT FUER MATHEMATIK UND DATENVERARBEITUNG (GMD). */

/* ----------------------------------------------------------------------- */


/* ===> LRATIO HAS TO BE SET TO THE NUMBER OF INTEGERS OCCUPYING THE SAME */
/*     AMOUNT OF STORAGE AS ONE DOUBLE PRECISION REAL. */


/* ===> MAXGR IS THE MAXIMAL NUMBER OF GRIDS. CHANGING THIS UPPER LIMIT */
/*     JUST REQUIRES CHANGING THE PARAMETER STATEMENT. */


    /* Parameter adjustments */
    --ig;
    --f;
    --u;
    --ja;
    --ia;
    --a;

    /* Function Body */
    *ierr = 0;

/* ===> SET PARAMETERS TO STANDARD VALUES, IF NECCESSARY */


    if (*iout != 0) {
	idec_(iout, &c__2, &ndigit, iarr);
	   kout = iarr[1];
    } else {
		kout = 2;
    }

    if (*iswtch != 0) {
	kswtch = *iswtch;
    } else {
	kswtch = 4;
    }

    if (*levelx > 0) {
		kevelx = my_imin(*levelx, 25);
	} else if (*levelx < 0) {
	goto L70;
    } else {
	kevelx = 25;
    }

    if (*iprint != 0) {
	idec_(iprint, &c__4, &ndigit, iarr);
	iup = iarr[1] * 10 + iarr[2];
	ium = iarr[3];
    } else {
	iup = 6;
	ium = 6;
    }
    icgst = *nnu + 3;
    ndicg = (*ndig - icgst + 1) / 2;
    if (ndicg <= 0) {
	goto L60;
    }

    switch (kswtch) {
	case 1:  goto L10;
	case 2:  goto L20;
	case 3:  goto L30;
	case 4:  goto L40;
    }
    //io___10.ciunit = ium;
    //s_wsfe(&io___10);
    //e_wsfe();

	
	printf("*** ERROR IN AMG1R5: ILLEGAL PARAMETER I.  SWTCH ***\n");
	
    *ierr = 17;
    return 0;

L40:
    setup_(nnu, matrix, &kevelx, ecg1, ecg2, ewt2, nwt, ntr, ierr, &a[1], &u[
	    1], &ia[1], &ja[1], &ig[1], imin, imax, iminw, imaxw, &ig[icgst], 
	    &ig[icgst + ndicg], nstcol, iarr, time, &levels, &irow0, nda, 
	    ndia, ndja, ndu, ndf, &ndicg, &ium, &mda, &mdia, &mdja, &mdu, &
	    mdf, &mdig, &c__25, &c__2);
    if (*ierr > 0) {
	return 0;
    }
L30:
    first_(ifirst, &u[1], imin, imax, iarr, &irow0);
L20:

	/*
	// Для корректности работы надо передавать информацию о модифицированном размере в вызывающие функции вверх по коду.
	if ((a != NULL) && (ja != NULL)) {

		// На данном этапе доступен истинный размер матрицы СЛАУ - хранимое число ненулевых элементов.
		// Вычислим реальное число ненулевых элементов матрицы a. Если номер столбца равен нулю то число ненулевых элементов
		// в матрице заведомо меньше чем позиция этого нуля.
		// С помощью низкоуровневой операции realloc ужимаем размер матрицы a.
		// Последующие выделения оперативной памяти для внешнего Крыловского итерационного процесса BiCGStab не приведут к ещё большему расходу
		// оперативной памяти, т.к. мы её освободили.

		printf("nda=%lld\n", nda[0]);
		integer isize97 = 0;
		for (integer i_95 = 1; i_95 < ndja[0]; i_95++) {
			if (ja[i_95] == 0) {
				isize97 = i_95;
				if (i_95 + 2 < ndja[0] && ja[i_95 + 1] == 0 && ja[i_95 + 2] == 0) {
					break;
				}
			}
		}
		printf("nda_new=%lld\n", isize97);
		++a;
		++ja;
		a = (doublereal*)realloc(a, (static_cast<integer>(isize97)+2) * sizeof(doublereal));
		ja = (integer*)realloc(ja, (static_cast<integer>(isize97)+2) * sizeof(integer));
		--a;
		--ja;
		//system("pause");
	}
	*/

	solve_(madapt, ncyc, nrd, nsolco, nru, &kout, ierr, &a[1], &u[1], &f[1], &
	    ia[1], &ja[1], &ig[1], eps, imin, imax, iminw, imaxw, &ig[icgst], 
	    &ig[icgst + ndicg], nstcol, iarr, time, &ncyc0, &irow0, &levels, 
	    nda, ndja, ndu, ndf, &mda, &mdja, &mdu, &mdf, &iup, &ium, resi, &
	    res0, &res);
	if (*ierr > 0) {
	return 0;
    }
L10:
    wrkcnt_(&kout, &ia[1], &ig[1], imin, imax, iminw, &levels, time, &ncyc0, &
	    iup, &mda, &mdia, &mdja, &mdu, &mdf, &mdig, &res0, &res);
    return 0;
L60:
    //io___29.ciunit = ium;
    //s_wsfe(&io___29);
    //e_wsfe();
	printf("*** ERROR IN AMG1R5: NDIG TOO SMALL ***\n");

    *ierr = 6;
    return 0;
L70:
    //io___30.ciunit = ium;
    //s_wsfe(&io___30);
    //e_wsfe();

	printf("*** ERROR IN AMG1R5: ILLEGAL PARAMETER LEVELX ***\n");

    *ierr = 18;
    return 0;
} 
//#if AMG1R6_LABEL==1
/* amg1r6_ */
//#else
/* amg1r5_ */
//#endif


//#if AMG1R6_LABEL==1

/*     Last change:  K    25 Jul 2002   12:55 pm */
/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */

/*     AMG1R6                                        MAIN SUBROUTINE */

/*     RELEASE 1.6, July 2002 */
/* 1.  changed: value of ntrim in pcol */
/* 2.  dimensioning (1) changed to (*) in some subroutines to avoid subscript */
/*     range checks in sparse solvers */

//#endif

  /*     Last change:  ERB  22 Aug 2000   10:31 am */
  /* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */

  /*     AMG1R5                                        MAIN SUBROUTINE */

  /*     RELEASE 1.5, OCTOBER 1990 */

  /*     CHANGES AGAINST VERSION 1.1, JULY 1985: */

  /* 1.  A BUG WAS DETECTED WHICH UNDER CERTAIN CIRCUMSTANCES INFLUENCED */
  /*     SLIGHTLY THE CONVERGENCE RATE OF AMG1R1. FOR THAT REASON, THE */
  /*     FOLLOWING LINE IN SUBROUTINE RESC: */
  /*     IW(IMAXW(KC-1)+1) = IA(IMIN(KC)) */
  /*     HAS BEEN CHANGED TO: */
  /*     IW(IMAXW(KC-1)+1) = IAUX */

  /* 2.  A BUG WAS DETECTED IN SUBROUTINE PWINT. UNDER CERTAIN CIRCUM- */
  /*     STANCES AN UNDEFINED VARIABLE WAS USED. ALTHOUGH THIS DID NOT */
  /*     AFFECT THE NUMERICAL RESULTS, PROBLEMS CAN OCCUR IF CHECKING */
  /*     FOR UNDEFINED VARIABLES IS USED. TO FIX THIS ERROR, IN PWINT */
  /*     THE LABEL 1000 WAS MOVED TO THE STATEMENT */
  /*     IBLCK1 = IMINW(K). */

  /* 3.  A PARAMETER LRATIO HAS BEEN INTRODUCED, DENOTING THE RATIO */
  /*     OF SPACE OCCUPIED BY A DOUBLE PRECISION REAL VARIABLE AND */
  /*     THAT OF AN INTEGER. FOR THE IBM-VERSION LRATIO HAS BEEN SET */
  /*     TO 2. CHANGE THIS VALUE IF NECESSARY. (IF, FOR EXAMPLE, YOU */
  /*     WANT TO CHANGE THE DOUBLE PRECISION VECTORS TO SINGLE PRE- */
  /*     CISION, LRATIO HAS TO BE SET TO 1. IN THE YALE SMP - ROUTINE */
  /*     NDRV THERE IS A PARAMETER LRATIO, TOO. */

  /* 4.  TYPE DECLARATIONS REAL*4 AND REAL*8 HAVE BEEN CHANGED TO THE */
  /*     STANDARD-CONFORMING KEYWORDS REAL AND DOUBLE PRECISION, RESPEC- */
  /*     TIVELY. */

  /* 5.  CALLS TO THE FOLLOWING INTRINSIC FUNCTIONS HAVE BEEN REPLACED BY */
  /*     CALLS USING GENERIC NAMES: DSQRT, MIN0, MAX0, IABS, DABS, FLOAT, */
  /*     DFLOAT, DMAX1, ISIGN, IDINT, DLOG10. */

  /* 6.  A SAVE STATEMENT HAS BEEN INSERTED IN ALL SUBROUTINES. */

  /* 7.  EXTERNAL DECLARATION STATEMENTS HAVE BEEN INSERTED IN ALL SUB- */
  /*     ROUTINES FOR ALL EXTERNAL REFERENCES. */

  /* ----------------------------------------------------------------------- */

  /*     CHANGE AGAINST VERSION 1.3, APRIL 1986: */

  /* 1.  A BUG IN SUBROUTINE CHECK HAS BEEN REMOVED. IF THE ORIGINAL MATRIX */
  /*     WAS STORED IN AN UNSYMMETRIC WAY, THE SYMMETRIZATION BY AMG1R3 */
  /*     COULD FAIL UNDER CERTAIN CIRCUMSTANCES. FOR A FIX, THE FOLLOWING */
  /*     STATEMENTS IN SUBROUTINE CHECK HAVE BEEN CHANGED: */

  /*     DO 450 J=IA(I)+1,IA(I+1)-1 WAS CHANGED TO */
  /*     DO 450 J=IA(I)+1,ICG(I)-1 */

  /*     DO 430 J1=IA(I1)+1,IA(I1+1)-1 WAS CHANGED TO */
  /*     DO 430 J1=IA(I1)+1,ICG(I1)-1 */

  /*     DO 550 J=IA(I)+1,IA(I+1)-1 WAS CHANGED TO */
  /*     DO 550 J=IA(I)+1,ICG(I)-1 */

  /*     DO 530 J1=IA(I1)+1,IA(I1+1)-1 WAS CHANGED TO */
  /*     DO 530 J1=IA(I1)+1,ICG(I1)-1 */

  /* 2.  THE EXPLANATORY PART IN SUBROUTINE AMG1R5 HAS BEEN ENLARGED TO */
  /*     AVOID MISUNDERSTANDINGS IN THE DEFINITION OF THE ARGUMENT LIST. */

  /* ----------------------------------------------------------------------- */

  /*     CHANGE AGAINST VERSION 1.4, OCTOBER, 1990 (BY JOHN W. RUGE) */

  /* 1.  A BUG IN SUBROUTINE CHECK HAS BEEN REMOVED. IF THE ORIGINAL MATRIX */
  /*     WAS STORED IN AN UNSYMMETRIC WAY, THE SYMMETRIZATION BY AMG1R3 */
  /*     COULD STILL FAIL UNDER CERTAIN CIRCUMSTANCES, AND WAS NOT FIXED */
  /*     IN THE PREVIOUS VERSION. IN ADDITION, THE ROUTINE WAS CHANGED */
  /*     IN ORDER TO AVOID SOME UNNECESSARY ROW SEARCHES FOR TRANSOSE */
  /*     ENTRIES. */

  /* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */



/* Subroutine */ integer amg1r5_fgmres_version(doublereal *a, integer *ia, integer *ja,
	doublereal *u, doublereal *f, integer *ig, integer *nda, integer *
	ndia, integer *ndja, integer *ndu, integer *ndf, integer *ndig,
	integer *nnu, integer *matrix, integer *iswtch, integer *iout,
	integer *iprint, integer *levelx, integer *ifirst, integer *ncyc,
	doublereal *eps, integer *madapt, integer *nrd, integer *nsolco,
	integer *nru, doublereal *ecg1, doublereal *ecg2, doublereal *ewt2,
	integer *nwt, integer *ntr, integer *ierr, integer iVar,
	equation3D* &sl, equation3D_bon* &slb, integer maxelm, integer maxbound,
	bool &bOkfgmres_amg1r5);


/* Subroutine */ integer amg1r5_fgmres_version_matrix_Assemble2(
	doublereal *a, integer *ia, integer *ja,// 16.10.2018
	doublereal *u, doublereal *f, integer *ig, integer *nda, integer *
	ndia, integer *ndja, integer *ndu, integer *ndf, integer *ndig,
	integer *nnu, integer *matrix, integer *iswtch, integer *iout,
	integer *iprint, integer *levelx, integer *ifirst, integer *ncyc,
	doublereal *eps, integer *madapt, integer *nrd, integer *nsolco,
	integer *nru, doublereal *ecg1, doublereal *ecg2, doublereal *ewt2,
	integer *nwt, integer *ntr, integer *ierr, integer iVar,
	SIMPLESPARSE &sparseM, integer n,
	bool &bOkfgmres_amg1r5);


//#if AMG1R6_LABEL==1

/*     Last change:  K    25 Jul 2002   12:55 pm */
/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */

/*     AMG1R6                                        MAIN SUBROUTINE */

/*     RELEASE 1.6, July 2002 */
/* 1.  changed: value of ntrim in pcol */
/* 2.  dimensioning (1) changed to (*) in some subroutines to avoid subscript */
/*     range checks in sparse solvers */

//#endif


  /*     Last change:  ERB  22 Aug 2000   10:31 am */
  /* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */

  /*     AMG1R5                                        MAIN SUBROUTINE */

  /*     RELEASE 1.5, OCTOBER 1990 */

  /*     CHANGES AGAINST VERSION 1.1, JULY 1985: */

  /* 1.  A BUG WAS DETECTED WHICH UNDER CERTAIN CIRCUMSTANCES INFLUENCED */
  /*     SLIGHTLY THE CONVERGENCE RATE OF AMG1R1. FOR THAT REASON, THE */
  /*     FOLLOWING LINE IN SUBROUTINE RESC: */
  /*     IW(IMAXW(KC-1)+1) = IA(IMIN(KC)) */
  /*     HAS BEEN CHANGED TO: */
  /*     IW(IMAXW(KC-1)+1) = IAUX */

  /* 2.  A BUG WAS DETECTED IN SUBROUTINE PWINT. UNDER CERTAIN CIRCUM- */
  /*     STANCES AN UNDEFINED VARIABLE WAS USED. ALTHOUGH THIS DID NOT */
  /*     AFFECT THE NUMERICAL RESULTS, PROBLEMS CAN OCCUR IF CHECKING */
  /*     FOR UNDEFINED VARIABLES IS USED. TO FIX THIS ERROR, IN PWINT */
  /*     THE LABEL 1000 WAS MOVED TO THE STATEMENT */
  /*     IBLCK1 = IMINW(K). */

  /* 3.  A PARAMETER LRATIO HAS BEEN INTRODUCED, DENOTING THE RATIO */
  /*     OF SPACE OCCUPIED BY A DOUBLE PRECISION REAL VARIABLE AND */
  /*     THAT OF AN INTEGER. FOR THE IBM-VERSION LRATIO HAS BEEN SET */
  /*     TO 2. CHANGE THIS VALUE IF NECESSARY. (IF, FOR EXAMPLE, YOU */
  /*     WANT TO CHANGE THE DOUBLE PRECISION VECTORS TO SINGLE PRE- */
  /*     CISION, LRATIO HAS TO BE SET TO 1. IN THE YALE SMP - ROUTINE */
  /*     NDRV THERE IS A PARAMETER LRATIO, TOO. */

  /* 4.  TYPE DECLARATIONS REAL*4 AND REAL*8 HAVE BEEN CHANGED TO THE */
  /*     STANDARD-CONFORMING KEYWORDS REAL AND DOUBLE PRECISION, RESPEC- */
  /*     TIVELY. */

  /* 5.  CALLS TO THE FOLLOWING INTRINSIC FUNCTIONS HAVE BEEN REPLACED BY */
  /*     CALLS USING GENERIC NAMES: DSQRT, MIN0, MAX0, IABS, DABS, FLOAT, */
  /*     DFLOAT, DMAX1, ISIGN, IDINT, DLOG10. */

  /* 6.  A SAVE STATEMENT HAS BEEN INSERTED IN ALL SUBROUTINES. */

  /* 7.  EXTERNAL DECLARATION STATEMENTS HAVE BEEN INSERTED IN ALL SUB- */
  /*     ROUTINES FOR ALL EXTERNAL REFERENCES. */

  /* ----------------------------------------------------------------------- */

  /*     CHANGE AGAINST VERSION 1.3, APRIL 1986: */

  /* 1.  A BUG IN SUBROUTINE CHECK HAS BEEN REMOVED. IF THE ORIGINAL MATRIX */
  /*     WAS STORED IN AN UNSYMMETRIC WAY, THE SYMMETRIZATION BY AMG1R3 */
  /*     COULD FAIL UNDER CERTAIN CIRCUMSTANCES. FOR A FIX, THE FOLLOWING */
  /*     STATEMENTS IN SUBROUTINE CHECK HAVE BEEN CHANGED: */

  /*     DO 450 J=IA(I)+1,IA(I+1)-1 WAS CHANGED TO */
  /*     DO 450 J=IA(I)+1,ICG(I)-1 */

  /*     DO 430 J1=IA(I1)+1,IA(I1+1)-1 WAS CHANGED TO */
  /*     DO 430 J1=IA(I1)+1,ICG(I1)-1 */

  /*     DO 550 J=IA(I)+1,IA(I+1)-1 WAS CHANGED TO */
  /*     DO 550 J=IA(I)+1,ICG(I)-1 */

  /*     DO 530 J1=IA(I1)+1,IA(I1+1)-1 WAS CHANGED TO */
  /*     DO 530 J1=IA(I1)+1,ICG(I1)-1 */

  /* 2.  THE EXPLANATORY PART IN SUBROUTINE AMG1R5 HAS BEEN ENLARGED TO */
  /*     AVOID MISUNDERSTANDINGS IN THE DEFINITION OF THE ARGUMENT LIST. */

  /* ----------------------------------------------------------------------- */

  /*     CHANGE AGAINST VERSION 1.4, OCTOBER, 1990 (BY JOHN W. RUGE) */

  /* 1.  A BUG IN SUBROUTINE CHECK HAS BEEN REMOVED. IF THE ORIGINAL MATRIX */
  /*     WAS STORED IN AN UNSYMMETRIC WAY, THE SYMMETRIZATION BY AMG1R3 */
  /*     COULD STILL FAIL UNDER CERTAIN CIRCUMSTANCES, AND WAS NOT FIXED */
  /*     IN THE PREVIOUS VERSION. IN ADDITION, THE ROUTINE WAS CHANGED */
  /*     IN ORDER TO AVOID SOME UNNECESSARY ROW SEARCHES FOR TRANSOSE */
  /*     ENTRIES. */

  /* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */

/* Subroutine */ integer amg1r5_Vorst_modification(doublereal *a, integer *ia, integer *ja,
	doublereal *u, doublereal *f, integer *ig, integer *nda, integer *
	ndia, integer *ndja, integer *ndu, integer *ndf, integer *ndig,
	integer *nnu, integer *matrix, integer *iswtch, integer *iout,
	integer *iprint, integer *levelx, integer *ifirst, integer *ncyc,
	doublereal *eps, integer *madapt, integer *nrd, integer *nsolco,
	integer *nru, doublereal *ecg1, doublereal *ecg2, doublereal *ewt2,
	integer *nwt, integer *ntr, integer *ierr, integer iVar, 
	equation3D* &sl, equation3D_bon* &slb, integer maxelm, integer maxbound)
{

	// 23-24 декабря 2017.
	// В данной функции amg1r5 алгоритм является предобуславливателем для алгоритма 
	// Хенка Ван дер Ворста BiCGStab[1992].

	// sl, slb - Матрица СЛАУ для алгоритма BiCGStab [1992] Хенка Ван дер Ворста.
	// Матрица СЛАУ (a, ia, ja) модифицируется в процессе работы алгебраического многосеточного метода.
	// априори: nnu[0]==maxelm+maxbound
	// ndu[0]==ndf[0] ...

	/* Format strings */


	/* Builtin functions */


	/* Local variables */
	integer mda = 0, mdf = 0, mdu = 0;
	doublereal res = 0.0;
	integer ium = 0, iup = 0;
	doublereal res0 = 0.0;
	extern /* Subroutine */ integer idec_(integer *, integer *, integer *,
		integer *);
	integer mdia = 0, mdja = 0, mdig = 0, iarr[25] = { 0 };
	//static real time[20];
	unsigned int time[20] = { 0 };
	integer imin[25] = { 0 }, imax[25] = { 0 };
	doublereal resi[25] = { 0.0 };
	integer kout = 0, ncyc0 = 0, irow0 = 0, ndicg = 0, icgst = 0, iminw[25] = { 0 }, imaxw[25] = { 0 };
	extern /* Subroutine */ integer first_(integer *, doublereal *, integer *,
		integer *, integer *, integer *), solve_(integer *, integer *,
			integer *, integer *, integer *, integer *, integer *, doublereal
			*, doublereal *, doublereal *, integer *, integer *, integer *,
			doublereal *, integer *, integer *, integer *, integer *, integer
			*, integer *, integer *, integer *, /*real*/ unsigned int *, integer *, integer *,
			integer *, integer *, integer *, integer *, integer *, integer *,
			integer *, integer *, integer *, integer *, integer *, doublereal
			*, doublereal *, doublereal *), setup_(integer *, integer *,
				integer *, doublereal *, doublereal *, doublereal *, integer *,
				integer *, integer *, doublereal *, doublereal *, integer *,
				integer *, integer *, integer *, integer *, integer *, integer *,
				integer *, integer *, integer *, integer *, /*real*/ unsigned int *, integer *,
				integer *, integer *, integer *, integer *, integer *, integer *,
				integer *, integer *, integer *, integer *, integer *, integer *,
				integer *, integer *, integer *, integer *);
	integer ndigit = 0, levels = 0, kevelx = 0, nstcol[25] = { 0 }, kswtch = 0;
	extern /* Subroutine */ integer wrkcnt_(integer *, integer *, integer *,
		integer *, integer *, integer *, integer *, /*real*/ unsigned int *, integer *,
		integer *, integer *, integer *, integer *, integer *, integer *,
		integer *, doublereal *, doublereal *);

	/* Fortran I/O blocks */

	bool debug_reshime = false;
	integer count_iter_for_film_coef75 = 0;
	doublereal delta_old_iter75 = 1.0e10;

	// Если число расходимостей превысит оговорённую константу то произойдёт выход из алгоритма.
	integer i_signal_break_pam_opening75 = 0;
	// x хорошее значение.
	const integer i_limit_signal_pam_break_opening75 = 4000;//20

	integer maxit75 = 2000;
	integer iN75 = 10;

	integer iflag75 = 1, icount75 = 0;
	doublereal delta075 = 1.0e30, deltai75 = 1.0e30;
	doublereal bet75 = 0.0, roi75 = 0.0;
	doublereal roim175 = 1.0, al75 = 1.0, wi75 = 1.0;

	doublereal epsilon75 = dterminatedTResudual;  // точность вычисления
	integer i75 = 0;
	integer iflag175 = 1;

	// Вектора необходимые для работы BiCGStab.
	doublereal* ri75 = NULL;
	doublereal* roc75 = NULL;
	doublereal* s75 = NULL;
	doublereal* t75 = NULL;
	doublereal* vec75 = NULL;
	doublereal* vi75 = NULL;
	doublereal* pi75 = NULL;
	doublereal* dx75 = NULL;
	doublereal* dax75 = NULL;
	doublereal* y75 = NULL;
	doublereal* z75 = NULL;
	// Первое предобуславливание:
	doublereal* y76 = NULL;
	doublereal* pi76 = NULL;
	// Второе предобуславливание:
	doublereal* z76 = NULL;
	doublereal* s76 = NULL;

	// нумерация векторов начинается с нуля.
	integer n75 = maxelm + maxbound; // число неизвестных на подробном уровне.
	doublereal* val75 = NULL;
	//val75 = new doublereal[nnz];
	integer* col_ind75 = NULL;
	//col_ind75 = new integer[nnz];
	integer* row_ptr75 = NULL;
	//row_ptr75 = new integer[n75 + 1];


	/*         ----------------------------------------------- */
	/*         | AMG-MODULE FOR SOLVING LINEAR SYSTEMS L*U=F | */
	/*         ----------------------------------------------- */

	/*         -------------------------------------------------------------- */

	/*     ASSUMPTIONS ON L: */

	/*         THE PROGRAM REQUIRES: */ // Требования к входным данным:

										/*             - DIAGONAL ENTRIES ARE ALWAYS POSITIVE (ON ALL GRIDS); */
										/*             - L IS A SQUARE MATRIX WHICH IS EITHER REGULAR OR SINGULAR */
										/*               WITH ROWSUMS=0. */
										// Матрица L это квадратная матрица с всегда положительными диагональными элементами с нулевой суммой коэффициентов в каждой строке (положительная определённость).

										/*         FOR THEORETICAL REASONS THE FOLLOWING SHOULD HOLD: */

										/*             - L POSITIVE DEFINITE (OR SEMI-DEFINITE WITH ROWSUM=0) */
										/*             - L "ESSENTIALLY" POSITIVE TYPE, I.E., */

										/*                  -- DIAGONAL ENTRIES MUST BE > 0 ; */
										/*                  -- MOST OF THE OFF-DIAGONAL ENTRIES <= 0 ; */
										/*                  -- ROWSUMS SHOULD BE >= 0 . */

										// Матрица L положительно определённая: диагональные элементы  строго больше нуля, большинство внедиагональных элементов <= 0.
										// Сумма коэффициентов в строке больше либо равна нулю - диагональное преобладание.


										/*     THE USER HAS TO PROVIDE THE MATRIX L, THE RIGHT HAND SIDE F AND */
										/*     CERTAIN POINTER VECTORS IA AND JA. */

										// Пользователь задаёт матрицу коэффициентов L, правую часть F, а также информацию о связях между элементами матрицы в специальных столбцах IA и JA.
										// IA - позиция первого элемента в строке. JA - номер столбца для каждого элемента, первым записывается диагональный элемент.
										// В фортране нумерация начинается с единицы.

										/*         -------------------------------------------------------------- */

										/*     STORAGE OF L: */ // Требования к хранению матрицы L.

																/*         THE NON-ZERO ENTRIES OF THE MATRIX L ARE STORED IN */
																/*         "COMPRESSED" SKY-LINE FASHION IN A 1-D VECTOR A, I.E., ROW */
																/*         AFTER ROW, EACH ROW STARTING WITH ITS DIAGONAL ELEMENT. THE */
																/*         OTHER NON-ZERO ROW ENTRIES FOLLOW THEIR DIAGONAL ENTRY IN ANY */
																/*         ORDER. */


																/*         IN ORDER TO IDENTIFY EACH ELEMENT IN A, THE USER HAS TO */
																/*         PROVIDE TWO POINTER ARRAYS IA AND JA. IF NNU DENOTES THE TOTAL */
																/*         NUMBER OF UNKNOWNS, THE NON-ZERO ENTRIES OF ANY ROW I OF L */
																/*         (1.LE.I.LE.NNU) ARE STORED IN A(J) WHERE THE RANGE OF J */
																/*         IS GIVEN BY */

																/*                     IA(I) .LE. J .LE. IA(I+1)-1. */

																/*         THUS, IA(I) POINTS TO THE POSITION OF THE DIAGONAL ENTRY OF */
																/*         ROW I WITHIN THE VECTOR A. IN PARTICULAR, */

																/*                     IA(1) = 1 ,  IA(NNU+1) = 1 + NNA */

																/*         WHERE NNA DENOTES THE TOTAL NUMBER OF MATRIX ENTRIES STORED. */
																/*         THE POINTER VECTOR JA HAS TO BE DEFINED SUCH THAT */
																/*         ANY ENTRY A(J) CORRESPONDS TO THE UNKNOWN U(JA(J)), I.E., */
																/*         JA(J) POINTS TO THE COLUMN INDEX OF A(J). */
																/*         IN PARTICULAR, A(IA(I)) IS THE DIAGONAL ENTRY OF ROW I */
																/*         AND CORRESPONDS TO THE UNKNOWN U(I): JA(IA(I))=I. */

																/*         IN THIS TERMINOLOGY, THE I-TH EQUATION READS AS FOLLOWS */
																/*         (FOR ANY I WITH  1.LE.I.LE.NNU): */

																/*                  F(I) =        SUM      A(J) * U(JA(J)) */
																/*                           J1.LE.J.LE.J2 */

																/*         WHERE F(I) DENOTES THE I-TH COMPONENT OF THE RIGHT HAND */
																/*         SIDE AND */

																/*                     J1 = IA(I) ,  J2 = IA(I+1)-1. */

																/*         NOTES: THE ENTRY IA(NNU+1) HAS TO TOCHKA TO THE FIRST FREE */
																/*                ENTRY IN VECTORS A AND JA, RESPECTIVELY. OTHERWISE, */
																/*                AMG CANNOT KNOW THE LENGTH OF THE LAST MATRIX ROW. */

																/*                THE INPUT VECTORS A, IA AND JA ARE CHANGED BY AMG1R5. */
																/*                SO, AFTER RETURN FROM AMG1R5, THE PACKAGE MUST NOT */
																/*                BE CALLED A SECOND TIME WITHOUT HAVING NEWLY DEFINED */
																/*                THE INPUT VECTORS AND USING ISWTCH=4. OTHERWISE, THE */
																/*                SETUP PHASE WILL FAIL. */
																/*                  ON THE OTHER HAND, RUNNING AMG A SECOND TIME ON THE */
																/*                SAME INPUT DATA WITH ISWTCH=4 HAS NO SENSE, BECAUSE */
																/*                THE RESULTS OF THE FIRST SETUP PHASE ARE STILL STORED */
																/*                AND THUS THIS PHASE CAN BE SKIPPED IN A SECOND CALL. */
																/*                IN ORDER TO DO THIS, SET ISWTCH TO 1, 2 OR 3. */

																/* ----------------------------------------------------------------------- */

																/*         THE FORM OF THE CALLING PROGRAM HAS TO BE AS FOLLOWS: */

																/*               PROGRAM DRIVER */
																/*         C */
																/*               DOUBLE PRECISION A(#NDA),U(#NDU),F(#NDF) */
																/*               INTEGER IA(#NDIA),JA(#NDJA),IG(#NDIG) */
																/*         C */
																/*               NDA  = #NDA */
																/*               NDU  = #NDU */
																/*               NDF  = #NDF */
																/*               NDIA = #NDIA */
																/*               NDJA = #NDJA */
																/*               NDIG = #NDIG */
																/*         C */
																/*         C     SET UP A, F, IA, JA AND SPECIFY NECESSARY PARAMETERS */
																/*         C */
																/*               .... */
																/*               .... */
																/*         C */
																/*               CALL AMG1R5(A,IA,JA,U,F,IG, */
																/*        +                  NDA,NDIA,NDJA,NDU,NDF,NDIG,NNU,MATRIX, */
																/*        +                  ISWTCH,IOUT,IPRINT, */
																/*        +                LEVELX,IFIRST,NCYC,EPS,MADAPT,NRD,NSOLCO,NRU, */
																/*        +                  ECG1,ECG2,EWT2,NWT,NTR, */
																/*        +                  IERR) */
																/*         C */
																/*               .... */
																/*               .... */
																/*         C */
																/*               STOP */
																/*               END */

																/* ----------------------------------------------------------------------- */

																/*     INPUT VIA ARRAYS (SEE ABOVE): */

																/*     A        -   MATRIX L */

																/*     IA       -   POINTER VECTOR */

																/*     JA       -   POINTER VECTOR */

																/*     U        -   FIRST APPROXIMATION TO SOLUTION */

																/*     F        -   RIGHT HAND SIDE */


																/* ----------------------------------------------------------------------- */


																/*     SCALAR INPUT PARAMETERS OF AMG1R5: */

																/*     THE INPUT PARAMETERS OF AMG1R5 IN THE LIST BELOW ARE ARRANGED */
																/*     ACCORDING TO THEIR IMPORTANCE TO THE GENERAL USER. THE PARAMETERS */
																/*     PRECEEDED BY A * MUST BE SPECIFIED EXPLICITELY. ALL THE OTHER */
																/*     PARAMETERS ARE SET TO STANDARD VALUES IF ZERO ON INPUT. */

																/*     THERE ARE FOUR CLASSES OF INPUT PARAMETERS WITH DECREASING PRI- */
																/*     ORITY: */

																/*     1. PARAMETERS DESCRIBING THE USER-DEFINED PROBLEM AND DIMENSIONING */
																/*        OF VECTORS IN THE CALLING PROGRAM */

																/*     2. PARAMETERS SPECIFYING SOME GENERAL ALGORITHMIC ALTERNATIVES AND */
																/*        THE AMOUNT OF OUTPUT DURING SOLUTION */

																/*     3. PARAMETERS CONTROLLING THE MULTIGRID CYCLING DURING THE SOLU- */
																/*        TION PHASE */

																/*     4. PARAMETERS CONTROLLING THE CREATION OF COARSER GRIDS AND INTER- */
																/*        POLATION FORMULAS. */

																/*     ONLY THE CLASS 1 - PARAMETERS MUST BE SPECIFIED EXPLICITELY BY */
																/*     THE USER. CLASS 2 - PARAMETERS CONTROL THE GENERAL PERFORMANCE OF */
																/*     AMG1R5. CHANGING THEM DOESN'T REQUIRE UNDERSTANDING THE AMG - */
																/*     ALGORITHM. SPECIFYING NON-STANDARD-VALUES FOR CLASS 3 - PARAMETERS */
																/*     PRESUPPOSES A GENERAL KNOWLEDGE OF MULTIGRID METHODS, WHEREAS THE */
																/*     FUNCTION OF CLASS 4 - PARAMETERS IS ONLY UNDERSTANDABLE AFTER */
																/*     STUDYING THE AMG-ALGORITHM IN DETAIL. FORTUNATELY IN MOST CASES */
																/*     THE CHOICE OF CLASS 3 AND 4 - PARAMETERS ISN'T CRITICAL AND USING */
																/*     THE AMG1R5 - SUPPLIED STANDARD VALUES SHOULD GIVE SATISFACTORY */
																/*     RESULTS. */

																/*         -------------------------------------------------------------- */

																/*     CLASS 1 - PARAMETERS: */

																/*  *  NDA      -   DIMENSIONING OF VECTOR A IN CALLING PROGRAM */

																/*  *  NDIA     -   DIMENSIONING OF VECTOR IA IN CALLING PROGRAM */

																/*  *  NDJA     -   DIMENSIONING OF VECTOR JA IN CALLING PROGRAM */

																/*  *  NDU      -   DIMENSIONING OF VECTOR U IN CALLING PROGRAM */

																/*  *  NDF      -   DIMENSIONING OF VECTOR F IN CALLING PROGRAM */

																/*  *  NDIG     -   DIMENSIONING OF VECTOR IG IN CALLING PROGRAM */

																/*  *  NNU      -   NUMBER OF UNKNOWNS */

																/*  *  MATRIX   -   INTEGER VALUE CONTAINING INFO ABOUT THE MATRIX L. */

																/*                  1ST DIGIT OF MATRIX  --  ISYM: */
																/*                    =1: L IS SYMMETRIC; */
																/*                    =2: L IS NOT SYMMETRIC. */

																/*                  2ND DIGIT OF MATRIX  --  IROW0: */
																/*                    =1: L HAS ROWSUM ZERO; */
																/*                    =2: L DOES NOT HAVE ROWSUM ZERO. */

																/*         -------------------------------------------------------------- */

																/*     CLASS 2 - PARAMETERS: */

																/*     ISWTCH   -   PARAMETER CONTROLLING WHICH MODULES OF AMG1R5 ARE TO */
																/*                  BE USED. */
																/*                    =1:   CALL FOR -----, -----, -----, WRKCNT. */
																/*                    =2:   CALL FOR -----, -----, SOLVE, WRKCNT. */
																/*                    =3:   CALL FOR -----, FIRST, SOLVE, WRKCNT. */
																/*                    =4:   CALL FOR SETUP, FIRST, SOLVE, WRKCNT. */
																/*                  SETUP DEFINES THE OPERATORS NEEDED IN THE SOLUTION */
																/*                         PHASE. */
																/*                  FIRST INITIALIZES THE SOLUTION VECTOR (SEE PARAMETER */
																/*                         IFIRST). */
																/*                  SOLVE COMPUTES THE SOLUTION BY AMG CYCLING (SEE */
																/*                         PARAMETER NCYC). */
																/*                  WRKCNT PROVIDES THE USER WITH INFORMATION ABOUT */
																/*                         RESIDUALS, STORAGE REQUIREMENTS AND CP-TIMES */
																/*                         (SEE PARAMETER IOUT). */
																/*                  IF AMG1R5 IS CALLED THE FIRST TIME, ISWTCH HAS TO */
																/*                  BE =4. INDEPENDENT OF ISWTCH, SINGLE MODULES CAN BE */
																/*                  BYPASSED BY A PROPER CHOICE OF THE CORRESPONDING */
																/*                  PARAMETER. */

																/*     IOUT     -   PARAMETER CONTROLLING THE AMOUNT OF OUTPUT DURING */
																/*                  SOLUTION PHASE: */

																/*                  1ST DIGIT: NOT USED; HAS TO BE NON-ZERO. */

																/*                  2ND DIGIT: */
																/*                    =0: NO OUTPUT (EXCEPT FOR MESSAGES) */
																/*                    =1: RESIDUAL BEFORE AND AFTER SOLUTION PROCESS */
																/*                    =2: ADD.: STATISTICS ON CP-TIMES AND STORAGE REQUI- */
																/*                        REMENTS */
																/*                    =3: ADD.: RESIDUAL AFTER EACH AMG-CYCLE */

																/*     IPRINT   -   PARAMETER SPECIFYING THE FORTRAN UNIT NUMBERS FOR */
																/*                  OUTPUT: */

																/*                  1ST DIGIT: NOT USED; HAS TO BE NON-ZERO */

																/*                  2ND AND 3RD DIGIT  --  IUP: UNIT NUMBER FOR RESULTS */

																/*                  4TH AND 5TH DIGIT  --  IUM: UNIT NUMBER FOR MESSAGES */

																/*         -------------------------------------------------------------- */

																/*     CLASS 3 - PARAMETERS: */

																/*     LEVELX   -   MAXIMUM NUMBER OF MG-LEVELS TO BE CREATED (>=1). */

																/*     IFIRST   -   PARAMETER FOR FIRST APPROXIMATION. */

																/*                  1ST DIGIT OF IFIRST: NOT USED; HAS TO BE NON-ZERO. */

																/*                  2ND DIGIT OF IFIRST  --  ITYPU: */
																/*                    =0: NO SETTING OF FIRST APPROXIMATION, */
																/*                    =1: FIRST APPROXIMATION CONSTANT TO ZERO, */
																/*                    =2: FIRST APPROXIMATION CONSTANT TO ONE, */
																/*                    =3: FIRST APPROXIMATION IS RANDOM FUNCTION WITH */
																/*                        THE CONCRETE RANDOM SEQUENCE BEING DETERMINED */
																/*                        BY THE FOLLWING DIGITS. */

																/*                  REST OF IFIRST  --  RNDU: */
																/*                    DETERMINES THE CONCRETE RANDOM SEQUENCE USED IN */
																/*                    THE CASE ITYPU=3. (IFIRST=13 IS EQUIVALENT TO */
																/*                    IFIRST=1372815) */

																/*     NCYC     -   INTEGER PARAMETER DESCRIBING THE TYPE OF CYCLE TO BE */
																/*                  USED AND THE NUMBER OF CYCLES TO BE PERFORMED. */

																/*                  1ST DIGIT OF NCYC  --  IGAM: */
																/*                    =1: V -CYCLE, */
																/*                    =2: V*-CYCLE, */
																/*                    =3: F -CYCLE, */
																/*                    =4: W -CYCLE. */
																/*                  IF NCYC IS NEGATIV, THEN THE APPROXIMATION OF THE */
																/*                  PROBLEM ON THE SECOND FINEST GRID IS COMPUTED BY */
																/*                  IGAM V-CYCLES ON THAT PARTICULAR GRID. */

																/*                  2ND DIGIT OF NCYC  --  ICGR: */
																/*                    =0: NO CONJUGATE GRADIENT, */
																/*                    =1: CONJUGATE GRADIENT (ONLY FIRST STEP OF CG), */
																/*                    =2: CONJUGATE GRADIENT (FULL CG). */

																/*                  3RD DIGIT OF NCYC  --  ICONV: */
																/*                    CONVERGENCE CRITERION FOR THE USER-DEFINED PROBLEM */
																/*                    (FINEST GRID): */
																/*                    =1: PERFORM A FIXED NUMBER OF CYCLES AS GIVEN BY */
																/*                        NCYCLE (SEE BELOW) */
																/*                    =2: STOP, IF  ||RES|| < EPS */
																/*                    =3: STOP, IF  ||RES|| < EPS * |F| */
																/*                    =4: STOP, IF  ||RES|| < EPS * |U| * |DIAG| */
																/*                    WITH ||RES|| = L2-NORM OF RESIDUAL, */
																/*                           EPS     (SEE INPUT PARAMETER EPS) */
																/*                           |F|   = SUPREMUM NORM OF RIGHT HAND SIDE */
																/*                           |U|   = SUPREMUM NORM OF SOLUTION */
																/*                         |DIAG|  = MAXIMAL DIAGONAL ENTRY IN MATRIX L */
																/*                    NOTE THAT IN ANY CASE THE SOLUTION PROCESS STOPS */
																/*                    AFTER AT MOST NCYCLE CYCLES. */

																/*                  REST OF NCYC  --  NCYCLE: */
																/*                    MAXIMAL NUMBER OF CYCLES TO BE PERFORMED (>0) OR */
																/*                    NCYCLE=0: NO CYCLING. */

																/*     EPS      -   CONVERGENCE CRITERION FOR SOLUTION PROCESS: (SEE */
																/*                  PARAMETER NCYC). NOTE THAT NO MORE THAN NCYCLE CYCLES */
																/*                  ARE PERFORMED, REGARDLESS OF EPS. */

																/*     MADAPT   -   INTEGER VALUE SPECIFYING THE CHOICE OF COARSEST */
																/*                  GRID IN CYCLING: */

																/*                  1ST DIGIT OF MADAPT  --  MSEL: */
																/*                    =1: IN CYCLING, ALL GRIDS CONSTRUCTED IN THE SETUP */
																/*                        PHASE ARE USED WITHOUT CHECK. */
																/*                    =2: THE NUMBER OF GRIDS IS AUTOMATICALLY REDUCED */
																/*                        IF THE CONVERGENCE FACTOR ON THE COARSER GRIDS */
																/*                        IS FOUND TO BE LARGER THAN A GIVEN VALUE FAC */
																/*                        (SEE BELOW). */

																/*                  REST OF MADAPT  --  FAC */
																/*                        THE REST OF MADAPT DEFINES THE FRACTIONAL PART */
																/*                        OF A REAL NUMBER FAC BETWEEN 0.1 AND 0.99, E.G. */
																/*                        MADAPT=258 MEANS MSEL=2 AND FAC=0.58. IF MADAPT */
																/*                        CONSISTS OF ONLY ONE DIGIT, FAC IS SET TO 0.7 */
																/*                        BY DEFAULT. */


																/*     NRD      -   PARAMETER DESCRIBING RELAXATION (DOWNWARDS): */

																/*                  1ST DIGIT OF NRD: NOT USED; HAS TO BE NON-ZERO. */

																/*                  2ND DIGIT OF NRD  --  NRDX: */
																/*                    ACTUAL NUMBER OF SMOOTHING STEPS TO BE PERFORMED */
																/*                    THE TYPE OF WHICH IS GIVEN BY THE FOLLOWING DIGITS */

																/*                  FOLLOWING DIGITS  --  ARRAY NRDTYP: */
																/*                    =1: RELAXATION OVER THE F-POINTS ONLY */
																/*                    =2: FULL GS SWEEP */
																/*                    =3: RELAXATION OVER THE C-POINTS ONLY */
																/*                    =4: FULL MORE COLOR SWEEP, HIGHEST COLOR FIRST */

																/*     NSOLCO   -   PARAMETER CONTROLLING THE SOLUTION ON COARSEST GRID: */

																/*                  1ST DIGIT  --  NSC: */
																/*                    =1: GAUSS-SEIDEL METHOD */
																/*                    =2: DIRECT SOLVER (YALE SMP) */

																/*                  REST OF NSOLCO  --  NRCX: (ONLY IF NSC=1) */
																/*                  NUMBER OF GS SWEEPS ON COARSEST GRID (>=0). */
																/*                  IF NRCX=0, THEN AS MANY GS SWEEPS ARE PERFORMED */
																/*                  AS ARE NEEDED TO REDUCE THE RESIDUAL BY TWO ORDERS */
																/*                  OF MAGNITUDE. (MAXIMAL 100 RELAXATION SWEEPS) */

																/*     NRU      -   PARAMETER FOR RELAXATION (UPWARDS), ANALOGOUS TO NRD. */

																/*         -------------------------------------------------------------- */

																/*     CLASS 4 - PARAMETERS: */

																/*     ECG1,ECG2-   REAL PARAMETERS AFFECTING THE CREATION OF COARSER */
																/*     EWT2     -   GRIDS AND/OR THE DEFINITION OF THE INTERPOLATION. */
																/*                  THE CHOICE OF THESE PARAMETERS DEPENDS ON */
																/*                  THE ACTUAL AMG VERSION (SEE SUBROUTINE CRSNG) */

																/*     NWT      -   INTEGER PARAMETER AFFECTING THE CREATION OF COARSER */
																/*                  GRIDS AND/OR THE DEFINITION OF THE INTERPOLATION. */
																/*                  THE CHOICE OF THIS PARAMETER DEPENDS ON */
																/*                  THE ACTUAL AMG VERSION (SEE SUBROUTINE CRSNG) */

																/*     NTR      -   PARAMETER CONTROLLING COARSE-GRID OPERATOR TRUNCATION */
																/*                    =0: PAIRS OF ZEROES ARE REMOVED FROM COARSE GRID */
																/*                        OPERATORS */
																/*                    =1: NO COARSE-GRID OPERATOR TRUNCATION */


																/* ----------------------------------------------------------------------- */

																/*     OUTPUT: */

																/*     U        -   CONTAINS THE COMPUTED SOLUTION */


																/*     IERR     -   ERROR PARAMETER: */

																/*                    >0: FATAL ERROR (ABNORMAL TERMINATION OF AMG1R5) */
																/*                    <0: NON-FATAL ERROR (EXECUTION OF AMG1R5 CONTINUES) */

																/*                  ERROR CODES IN DETAIL: */

																/*                  1. DIMENSIONING TOO SMALL FOR VECTOR */
																/*                        A      (IERR = 1) */
																/*                        IA     (IERR = 2) */
																/*                        JA     (IERR = 3) */
																/*                        U      (IERR = 4) */
																/*                        F      (IERR = 5) */
																/*                        IG     (IERR = 6) */

																/*                     NO YALE-SMP BECAUSE OF STORAGE (NDA TOO SMALL): */
																/*                               (IERR = -1) */
																/*                     NO YALE-SMP BECAUSE OF STORAGE (NDJA TOO SMALL): */
																/*                               (IERR = -3) */
																/*                     NO CG BECAUSE OF STORAGE (NDU TOO SMALL): */
																/*                               (IERR = -4) */
																/*                     NO SPACE FOR TRANSPOSE OF INTERPOLATION (NDA OR */
																/*                                                     NDJA TOO SMALL): */
																/*                               (IERR = -1) */

																/*                  2. INPUT DATA ERRONEOUS: */

																/*                     A-ENTRY MISSING, ISYM = 1:           (IERR = -11) */
																/*                     PARAMETER MATRIX MAY BE ERRONEOUS:   (IERR = -12) */
																/*                     DIAGONAL ELEMENT NOT STORED FIRST:   (IERR =  13) */
																/*                     DIAGONAL ELEMENT NOT POSITIV:        (IERR =  14) */
																/*                     POINTER IA ERRONEOUS:                (IERR =  15) */
																/*                     POINTER JA ERRONEOUS:                (IERR =  16) */
																/*                     PARAMETER ISWTCH ERRONEOUS:          (IERR =  17) */
																/*                     PARAMETER LEVELX ERRONEOUS:          (IERR =  18) */

																/*                  3. ERRORS OF THE AMG1R5-SYSTEM (SHOULD NOT OCCUR): */

																/*                     TRANSPOSE A-ENTRY MISSING:           (IERR =  21) */
																/*                     INTERPOLATION ENTRY MISSING:         (IERR =  22) */

																/*                  4. ALGORITHMIC ERRORS: */

																/*                     CG-CORRECTION NOT DEFINED:           (IERR =  31) */
																/*                     NO YALE-SMP BECAUSE OF ERROR IN */
																/*                     FACTORIZATION:                       (IERR = -32) */

																/* ----------------------------------------------------------------------- */

																/*     WORK SPACE: */

																/*     THE INTEGER VECTOR IG HAS TO BE PASSED TO AMG1R5 AS WORK SPACE. */

																/* ----------------------------------------------------------------------- */

																/*     DIMENSIONING OF INPUT VECTORS AND WORK SPACE: */

																/*     IT'S IMPOSSIBLE TO TELL IN ADVANCE THE EXACT STORAGE REQUIREMENTS */
																/*     OF AMG. THUS, THE FOLLOWING FORMULAS GIVE ONLY REASONABLE GUESSES */
																/*     FOR THE VECTOR LENGTHS WHICH HAVE TO BE DECLARED IN THE CALLING */
																/*     PROGRAM. IN THESE FORMULAS NNA DENOTES THE NUMBER OF NON-ZERO */
																/*     ENTRIES IN THE INPUT-MATRIX L AND NNU IS THE NUMBER OF UNKNOWNS. */

																/*     VECTOR         NEEDED LENGTH (GUESS) */
																/*       A               3*NNA + 5*NNU */
																/*       JA              3*NNA + 5*NNU */
																/*       IA              2.2*NNU */
																/*       U               2.2*NNU */
																/*       F               2.2*NNU */
																/*       IG              5.4*NNU */

																/* ----------------------------------------------------------------------- */


																/*     STANDARD CHOICES OF PARAMETERS (AS FAR AS MEANINGFUL): */

																/*          ISWTCH = 4 */
																/*          IOUT   = 12 */
																/*          IPRINT = 10606 */

																/*          LEVELX = 25 */
																/*          IFIRST = 13 */
																/*          NCYC   = 10110 */
																/*          EPS    = 1.D-12 */
																/*          MADAPT = 27 */
																/*          NRD    = 1131 */
																/*          NSOLCO = 110 */
																/*          NRU    = 1131 */

																/*          ECG1   = 0. */
																/*          ECG2   = 0.25 */
																/*          EWT2   = 0.35 */
																/*          NWT    = 2 */
																/*          NTR    = 0 */

																/*     IF ANY ONE OF THESE PARAMETERS IS 0 ON INPUT, ITS CORRESPONDING */
																/*     STANDARD VALUE IS USED BY AMG1R5. */

																/* ----------------------------------------------------------------------- */

																/*     PORTABILITY RESTRICTIONS: */

																/*     1. ROUTINE CTIME IS MACHINE DEPENDENT AND HAS TO BE ADAPTED TO */
																/*        YOUR COMPUTER INSTALLATION OR REPLACED BY A DUMMY ROUTINE. */

																/*     2. MOST INPUT PARAMETERS ARE COMPOSED OF SEVERAL DIGITS, THEIR */
																/*        SIGNIFICANCE HAVING BEEN DESCRIBED ABOVE. BE SURE NOT TO ENTER */
																/*        MORE DIGITS THAN YOUR COMPUTER CAN STORE ON AN INTEGER VARI- */
																/*        ABLE. */

																/*     3. APART FROM FORTRAN INTRINSIC FUNCTIONS AND SERVICE ROUTINES, */
																/*        THERE IS ONLY ONE EXTERNAL REFERENCE TO A PROGRAM NOT CONTAINED */
																/*        IN THE AMG1R5 - SYSTEM, I.E. THE LINEAR SYSTEM SOLVER NDRV OF */
																/*        THE YALE SPARSE MATRIX PACKAGE. IF YOU HAVN'T ACCESS TO THIS */
																/*        PACKAGE, ENTER A DUMMY ROUTINE NDRV AND AVOID CHOOSING NSC=2 */
																/*        (SUBPARAMETER OF NSOLCO). THEN NDRV ISN'T CALLED BY AMG1R5. */
																/*        IN THIS CASE, HOWEVER, INDEFINITE PROBLEMS WILL NOT BE SOLV- */
																/*        ABLE. */
																/*          THE YALE SPARSE MATRIX PACKAGE IS FREELY AVAILABLE FOR NON- */
																/*        PROFIT PURPOSES. CONTACT THE DEPARTMENT OF COMPUTER SCIENCE, */
																/*        YALE UNITVERSITY. */

																/*     4. IN AMG1R5 THERE IS THE PARAMETER LRATIO, DENOTING THE RATIO */
																/*        OF SPACE OCCUPIED BY A DOUBLE PRECISION REAL VARIABLE AND */
																/*        THAT OF AN INTEGER. FOR THE IBM-VERSION LRATIO HAS BEEN SET */
																/*        TO 2. CHANGE THIS VALUE IF NECESSARY. (THE SAME HAS TO BE */
																/*        DONE WITH THE YALE SMP-ROUTINE NDRV.) */


																/* ----------------------------------------------------------------------- */

																/*     AUTHORS: */

																/*          JOHN RUGE, FORT COLLINS (USA), */
																/*              INSTITUTE FOR COMPUTATIONAL STUDIES AT CSU; */

																/*          KLAUS STUEBEN, D-5205 ST. AUGUSTIN (W.-GERMANY), */
																/*              GESELLSCHAFT FUER MATHEMATIK UND DATENVERARBEITUNG (GMD). */

																/*          ROLF HEMPEL, D-5205 ST. AUGUSTIN (W.-GERMANY), */
																/*              GESELLSCHAFT FUER MATHEMATIK UND DATENVERARBEITUNG (GMD). */

																/* ----------------------------------------------------------------------- */


																/* ===> LRATIO HAS TO BE SET TO THE NUMBER OF INTEGERS OCCUPYING THE SAME */
																/*     AMOUNT OF STORAGE AS ONE DOUBLE PRECISION REAL. */


																/* ===> MAXGR IS THE MAXIMAL NUMBER OF GRIDS. CHANGING THIS UPPER LIMIT */
																/*     JUST REQUIRES CHANGING THE PARAMETER STATEMENT. */


																/* Parameter adjustments */
	--ig;
	--f;
	--u;
	--ja;
	--ia;
	--a;

	/* Function Body */
	*ierr = 0;

	/* ===> SET PARAMETERS TO STANDARD VALUES, IF NECCESSARY */


	if (*iout != 0) {
		idec_(iout, &c__2, &ndigit, iarr);
		kout = iarr[1];
	}
	else {
		kout = 2;
	}

	if (*iswtch != 0) {
		kswtch = *iswtch;
	}
	else {
		kswtch = 4;
	}

	if (*levelx > 0) {
		kevelx = my_imin(*levelx, 25);
	}
	else if (*levelx < 0) {
		goto L70;
	}
	else {
		kevelx = 25;
	}

	if (*iprint != 0) {
		idec_(iprint, &c__4, &ndigit, iarr);
		iup = iarr[1] * 10 + iarr[2];
		ium = iarr[3];
	}
	else {
		iup = 6;
		ium = 6;
	}
	icgst = *nnu + 3;
	ndicg = (*ndig - icgst + 1) / 2;
	if (ndicg <= 0) {
		goto L60;
	}

	

	switch (kswtch) {
	case 1:  goto L10;
	case 2:  goto L20;
	case 3:  goto L30;
	case 4:  goto L40;
	}
	//io___10.ciunit = ium;
	//s_wsfe(&io___10);
	//e_wsfe();


	printf("*** ERROR IN AMG1R5: ILLEGAL PARAMETER I.  SWTCH ***\n");

	*ierr = 17;
	return 0;

L40:
	setup_(nnu, matrix, &kevelx, ecg1, ecg2, ewt2, nwt, ntr, ierr, &a[1], &u[
		1], &ia[1], &ja[1], &ig[1], imin, imax, iminw, imaxw, &ig[icgst],
			&ig[icgst + ndicg], nstcol, iarr, time, &levels, &irow0, nda,
			ndia, ndja, ndu, ndf, &ndicg, &ium, &mda, &mdia, &mdja, &mdu, &
			mdf, &mdig, &c__25, &c__2);
	if (*ierr > 0) {
		return 0;
	}
L30:
	first_(ifirst, &u[1], imin, imax, iarr, &irow0);
L20:

	// Алгебраический Многосеточный Метод как предобуславливатель
	// к алгоритму Крыловского типа Хенка Ван Дер Ворста BiCGStab
	// со стабилизацией.
	// Требует ещё одну память под матрицу А на самом подробном уровне.
	// 5.01.2017 Алгоритм BiCGStab изобретён в 1992 году.
	
	/*
	// Для корректности работы надо передавать информацию о модифицированном размере в вызывающие функции вверх по коду.
	if ((a != NULL)&&(ja!=NULL)) {

		// На данном этапе доступен истинный размер матрицы СЛАУ - хранимое число ненулевых элементов.
		// Вычислим реальное число ненулевых элементов матрицы a. Если номер столбца равен нулю то число ненулевых элементов
		// в матрице заведомо меньше чем позиция этого нуля.
		// С помощью низкоуровневой операции realloc ужимаем размер матрицы a.
		// Последующие выделения оперативной памяти для внешнего Крыловского итерационного процесса BiCGStab не приведут к ещё большему расходу
		// оперативной памяти, т.к. мы её освободили.

		printf("nda=%lld\n", nda[0]);
		integer isize97 = 0;
		for (integer i_95 = 1; i_95 < ndja[0]; i_95++) {
			if (ja[i_95] == 0) {
				isize97 = i_95;
				if (i_95 + 2 < ndja[0] && ja[i_95 + 1] == 0 && ja[i_95 + 2] == 0) {
					break;
				}
			}
		}
		printf("nda_new=%lld\n", isize97);
		++a;
		++ja;
		a = (doublereal*)realloc(a, (static_cast<integer>(isize97)+2) * sizeof(doublereal));
		ja = (integer*)realloc(ja, (static_cast<integer>(isize97)+2) * sizeof(integer));
		--a;
		--ja;
		//system("pause");
	}
	*/

	ncyc[0] = 1011; // Предобуславливание в один V цикл.

	/*
	// для отладки.
	for (integer i_numberV_cycle = 0; i_numberV_cycle < 10; i_numberV_cycle++) {
		ncyc[0] = 1011;
		solve_(madapt, ncyc, nrd, nsolco, nru, &kout, ierr, &a[1], &u[1], &f[1], &
			ia[1], &ja[1], &ig[1], eps, imin, imax, iminw, imaxw, &ig[icgst],
			&ig[icgst + ndicg], nstcol, iarr, time, &ncyc0, &irow0, &levels,
			nda, ndja, ndu, ndf, &mda, &mdja, &mdu, &mdf, &iup, &ium, resi, &
			res0, &res);
		system("pause");
	}
	*/
	printf("sizeof  ndu=%lld nnu=%lld ndf=%lld\n",ndu[0],nnu[0],ndf[0]);
	integer nnz;
	/*
	integer n75 = -1 , nnz;

	for (integer i_83 = 1; i_83 <= nda[0]; i_83++) {
		if (ja[i_83] > n75) n75 = ja[i_83];
	}
	nnz = ia[n75 + 1];
	printf("n_75=%d nnz=%d\n",n75,nnz);
	system("pause");
	*/
	


	// Разреженная матрица СЛАУ
	// в CRS формате.


	if ((iVar == VX) || (iVar == VY) || (iVar == VZ) || (iVar == PAM)) {
		/*
		if (ibackregulationgl != NULL) {
			// nested desection версия алгоритма.
			integer ierr = equation3DtoCRSnd(sl, slb, val75, col_ind75, row_ptr75, maxelm, maxbound, 1.0, true, NULL, NULL);
			if (ierr > 0) {
				switch (iVar) {
				case VX: printf("VX equation problem.\n"); break;
				case VY: printf("VY equation problem.\n"); break;
				case VZ: printf("VZ equation problem.\n"); break;
				case PAM: printf("PAM equation problem.\n"); break;
				}
			}
		}
		else {
		*/
			integer ierr = equation3DtoCRS(sl, slb, val75, col_ind75, row_ptr75, maxelm, maxbound, 1.0, true);
			if (ierr > 0) {
				switch (iVar) {
				case VX: printf("VX equation problem.\n"); break;
				case VY: printf("VY equation problem.\n"); break;
				case VZ: printf("VZ equation problem.\n"); break;
				case PAM: printf("PAM equation problem.\n"); break;
				}
			}
		//}
	}
	if (iVar == TEMP) {
		integer ierr = equation3DtoCRS(sl, slb, val75, col_ind75, row_ptr75, maxelm, maxbound, 1.0, true);
		if (ierr > 0) {
			printf("Temperature equation problem.\n");
		}
	}

	if ((val75 == NULL) || (col_ind75 == NULL) || (row_ptr75 == NULL)) {
		// недостаточно памяти на данном оборудовании.
		printf("Problem: not enough memory on your equipment for val, col_ind or row_ptr: bicgStab + camg...\n");
		printf("Please any key to exit...\n");
		exit(1);
	}

	nnz = row_ptr75[n75];
	printf("n=%lld nnz=%lld ndu=%lld nda=%lld\n", n75, nnz, ndu[0], nda[0]);

	

	/*
	// Так делать нельзя, т.к. матрица СЛАУ (a,ia,ja) 
	// была модифицирована в процессе работы amg1r5.f алгоритма.

	// инициализация матрицы.
#pragma omp parallel for
	for (integer i_91 = 1; i_91 <= n75; i_91++) {

		for (integer i_92 = ia[i_91]; i_92 < ia[i_91 + 1]; i_92++) {
			
			val75[i_92 - 1] = a[i_92];
			col_ind75[i_92 - 1] = ja[i_92]- 1;
			
		}
		row_ptr75[i_91 - 1] = ia[i_91] - 1;
	}
	row_ptr75[n75] = ia[n75 + 1] - 1;
	*/

	//system("PAUSE");

	
	
	

	//y76 = new doublereal[n75 + 1];
	//pi76 = new doublereal[n75 + 1];
	y76 = new doublereal[ndu[0] + 1];
	pi76 = new doublereal[ndu[0] + 1];
	
	//z76 = new doublereal[n75 + 1];
	//s76 = new doublereal[n75 + 1];
	z76 = new doublereal[ndu[0] + 1];
	s76 = new doublereal[ndu[0] + 1];

	ri75 = new doublereal[n75];
	roc75 = new doublereal[n75];
	s75 = new doublereal[n75];
	t75 = new doublereal[n75];
	vec75 = new doublereal[n75];
	vi75 = new doublereal[n75];
	pi75 = new doublereal[n75];
	dx75 = new doublereal[n75];
	dax75 = new doublereal[n75];
	y75 = new doublereal[n75];
	z75 = new doublereal[n75];
	if ((ri75 == NULL) || (roc75 == NULL) || (s75 == NULL) || (t75 == NULL) || (vi75 == NULL) || (pi75 == NULL) || (dx75 == NULL) || (dax75 == NULL) || (y75 == NULL) || (z75 == NULL)) {
		// недостаточно памяти на данном оборудовании.
		printf("Problem: not enough memory on your equipment for: bicgStab + camg...\n");
		printf("Please any key to exit...\n");
		exit(1);
	}


	

	if (iVar == TEMP) {
		epsilon75 *= 1.0e-4; // 1.0e-4
	}
	if (iVar == TOTALDEFORMATIONVAR) {
		epsilon75 *= 1.0e-4; // 1.0e-4
		//epsilon75 *= 1.0e-12;
	}
	

	// initialize
#pragma omp parallel for
	for (i75 = 0; i75<n75; i75++) {
		s75[i75] = 0.0;
		t75[i75] = 0.0;
		vi75[i75] = 0.0;
		pi75[i75] = 0.0;
		// инициализатор массивов для предобуславливания
		y75[i75] = 0.0;
		z75[i75] = 0.0;
		// результат умножения матрицы на вектор.
		dax75[i75] = 0.0;
		// Начальное приближение.
		dx75[i75] = u[i75 + 1];
	}


	// Умножение матрицы на вектор. Нумерации векторов начинаются с нуля.
	MatrixCRSByVector(val75, col_ind75, row_ptr75, dx75, dax75, n75); // результат занесён в  dax75

																	  // Вычисление ri75 и roc75.
#pragma omp parallel for
	for (i75 = 0; i75 < n75; i75++) {
		ri75[i75] = f[i75 + 1] - dax75[i75];
		roc75[i75] = 1.0;
	}
	delta075 = NormaV(ri75, n75);


	// Если решение сразу хорошее то не считать:
	if (iVar == TEMP) {
		if (fabs(delta075)<1.0e-4*dterminatedTResudual) iflag75 = 0;
	}
	else {
		if (fabs(delta075)<dterminatedTResudual) iflag75 = 0;
	}
	
	if (fabs(delta075)<1e-14) iflag175 = 0;
	if ((iVar == TEMP) && (iflag75 == 0) && (iflag175 == 0)) {
#if doubleintprecision == 1
		printf("bicgStab+camg: iflag=%lld, iflag1=%lld, delta0=%e\n", iflag75, iflag175, delta075);
#else
		printf("bicgStab+camg: iflag=%d, iflag1=%d, delta0=%e\n", iflag75, iflag175, delta075);
#endif

		system("PAUSE");
	}

	
	if (n75<30000) {
		// задача очень малой размерности !
		if ((iVar == VX) || (iVar == VY) || (iVar == VZ)) {
			iN75 = 1; // обязательно нужна хотя бы одна итерация.
					  // если этого будет недостаточно то мы всё равно будем итерировать до тех пор пока невязка не станет меньше epsilon.
			if (1.0e-3*fabs(delta075)<epsilon75) {
				epsilon75 = 1.0e-3*fabs(delta075);
			}
			if (iflag175 == 1) {
				iflag75 = 1;
			}
		}
		if (iVar == TEMP) {
			iN75 = 2;
			epsilon75 = fmin(0.1*fabs(delta075), epsilon75);
			if (bSIMPLErun_now_for_temperature  ) {
				//printf("epsilon=%e \n",epsilon);
				//system("pause");
				// Экспериментальным образом обнаружена недоэтерированость по температуре для гидродинамического решателя.
				// поэтому точность было решено увеличить на 5 порядков.
				// 27.07.2016
				epsilon75 *= 1e-10;
				iN75 = 20;
				//epsilon75 *= 1e-16;
				//iN75 = 30;
			}
		}
		if (iVar == PAM) {
			iN75 = 3; // решение для поправки давления должно быть получено точно.
			if (1.0e-3*fabs(delta075)<epsilon75) {
				epsilon75 = 1.0e-3*fabs(delta075);
			}
			if (iflag175 == 1) {
				iflag75 = 1;
			}
			//printf("%e",epsilon75); system("pause");
		}
	}
	else if ((n75 >= 30000) && (n75 < 100000)) {
		// Здесь я немного увеличил число итераций и 
		// скорректировал условие окончания чтобы считало 
		// поточнее, но это не повлияло.
		// Главный вопрос в том что невязка по температуре почему-то не меняется.
		// задача небольшой размерности.
		if ((iVar == VX) || (iVar == VY) || (iVar == VZ)) {
			//iN75 = 3; // обязательно нужна хотя бы одна итерация.
					  // если этого будет недостаточно то мы всё равно будем итерировать до тех пор пока невязка не станет меньше epsilon.
			if (1.0e-3*fabs(delta075)<epsilon75) {
				epsilon75 = 1.0e-3*fabs(delta075);
			}
			if (iflag175 == 1) {
				iflag75 = 1;
			}
			// 27.07.2016
			iN75 = 12;
			epsilon75 *= 1e-2;
		}
		if (iVar == TEMP) {
			iN75 = 4;
			epsilon75 = fmin(0.1*fabs(delta075), epsilon75);
			if (bSIMPLErun_now_for_temperature  ) {
				//printf("epsilon75=%e \n",epsilon75);
				//system("pause");
				// Экспериментальным образом обнаружена недоэтерированость по температуре для гидродинамического решателя.
				// поэтому точность было решено увеличить на 5 порядков.
				// 27.07.2016
				epsilon75 *= 1e-10;
				iN75 = 20;
				//epsilon75 *= 1e-16;
				//iN75 = 30;
			}
		}
		if (iVar == PAM) {
			//iN75 = 6; // решение для поправки давления должно быть получено точно.
			if (1.0e-3*fabs(delta075)<epsilon75) {
				epsilon75 = 1.0e-3*fabs(delta075);
			}
			if (iflag175 == 1) {
				iflag75 = 1;
			}
			//printf("%e",epsilon75); system("pause");
			// 27.07.2016.
			epsilon75 *= 1e-2;
			iN75 = 20;
		}
	}
	else if ((n75 >= 100000) && (n75<300000)) {
		// задача небольшой средней размерности.
		if ((iVar == VX) || (iVar == VY) || (iVar == VZ)) {
			iN75 = 3; // обязательно нужна хотя бы одна итерация.
					  // Вообще говоря невязка для скоростей падает очень быстро поэтому всегда достаточно iN итераций для скорости.
					  // если этого будет недостаточно то мы всё равно будем итерировать до тех пор пока невязка не станет меньше epsilon.
			if (1.0e-3*fabs(delta075)<epsilon75) {
				epsilon75 = 1.0e-3*fabs(delta075);
			}
			if (iflag175 == 1) {
				iflag75 = 1;
			}
		}
		if (iVar == TEMP) {
			iN75 = 4;
			epsilon75 = fmin(0.1*fabs(delta075), epsilon75);
			if (bSIMPLErun_now_for_temperature  ) {
				//printf("epsilon75=%e \n",epsilon75);
				//system("pause");
				// Экспериментальным образом обнаружена недоэтерированость по температуре для гидродинамического решателя.
				// поэтому точность было решено увеличить на 5 порядков.
				// 27.07.2016
				epsilon75 *= 1e-10;
				iN75 = 20;
				//epsilon75 *= 1e-16;
				//iN75 = 30;
			}
		}
		if (iVar == PAM) {
			iN75 = 8; // решение для поправки давления должно быть получено точно.
			if (1.0e-4*fabs(delta075)<epsilon75) {
				epsilon75 = 1.0e-4*fabs(delta075);
			}
			if (iflag175 == 1) {
				iflag75 = 1;
			}
			//printf("%e",epsilon75); system("pause");
		}
	}
	else if ((n75 >= 300000) && (n75<1000000)) {
		// задача истинно средней размерности.
		if ((iVar == VX) || (iVar == VY) || (iVar == VZ)) {
			iN75 = 3; // обязательно нужна хотя бы одна итерация.
					  // если этого будет недостаточно то мы всё равно будем итерировать до тех пор пока невязка не станет меньше epsilon.
			if (1.0e-3*fabs(delta075)<epsilon75) {
				epsilon75 = 1.0e-3*fabs(delta075);
			}
			if (iflag175 == 1) {
				iflag75 = 1;
			}
		}
		if (iVar == TEMP) {
			iN75 = 4;
			epsilon75 = 1e-5*fmin(0.1*fabs(delta075), epsilon75);
			if (bSIMPLErun_now_for_temperature  ) {
				//printf("epsilon75=%e \n",epsilon75);
				//system("pause");
				// Экспериментальным образом обнаружена недоэтерированость по температуре для гидродинамического решателя.
				// поэтому точность было решено увеличить на 5 порядков.
				// 27.07.2016
				epsilon75 *= 1e-8;
				iN75 = 20;
				//epsilon75 *= 1e-16;
				//iN75 = 30;
			}
		}
		if (iVar == PAM) {
			iN75 = 16; // решение для поправки давления должно быть получено точно.
			if (1.0e-4*fabs(delta075)<epsilon75) {
				epsilon75 = 1.0e-4*fabs(delta075);
			}
			if (iflag175 == 1) {
				iflag75 = 1;
			}
			//printf("%e",epsilon75); system("pause");
		}
	}
	else if ((n75 >= 1000000) && (n75<3000000)) {
		// задача достаточно большой размерности.
		if ((iVar == VX) || (iVar == VY) || (iVar == VZ)) {
			iN75 = 6; // обязательно нужна хотя бы одна итерация.
					  // если этого будет недостаточно то мы всё равно будем итерировать до тех пор пока невязка не станет меньше epsilon.
			if (1.0e-3*fabs(delta075)<epsilon75) {
				epsilon75 = 1.0e-3*fabs(delta075);
			}
			if (iflag175 == 1) {
				iflag75 = 1;
			}
		}
		if (iVar == TEMP) {
			iN75 = 8;
			epsilon75 = 1e-5*fmin(0.1*fabs(delta075), epsilon75);
			if (bSIMPLErun_now_for_temperature  ) {
				//printf("epsilon75=%e \n",epsilon75);
				//system("pause");
				// Экспериментальным образом обнаружена недоэтерированость по температуре для гидродинамического решателя.
				// поэтому точность было решено увеличить на 5 порядков.
				// 27.07.2016
				epsilon75 *= 1e-8;
				iN75 = 20;
				//epsilon75 *= 1e-16;
				//iN75 = 30;
			}
		}
		if (iVar == PAM) {
			iN75 = 23; // решение для поправки давления должно быть получено точно.
			if (1.0e-4*fabs(delta075)<epsilon75) {
				epsilon75 = 1.0e-4*fabs(delta075);
			}
			if (iflag175 == 1) {
				iflag75 = 1;
			}
			//printf("%e",epsilon75); system("pause");
		}
	}
	else if (n75 >= 3000000) {
		// задача очень большой размерности.
		if ((iVar == VX) || (iVar == VY) || (iVar == VZ)) {
			iN75 = 6; // обязательно нужна хотя бы одна итерация.
					  // если этого будет недостаточно то мы всё равно будем итерировать до тех пор пока невязка не станет меньше epsilon.
			if (1.0e-3*fabs(delta075)<epsilon75) {
				epsilon75 = 1.0e-3*fabs(delta075);
			}
			if (iflag175 == 1) {
				iflag75 = 1;
			}
		}
		if (iVar == TEMP) {
			iN75 = 8;
			epsilon75 = 1e-10*fmin(0.1*fabs(delta075), epsilon75);
		}
		if (iVar == PAM) {
			iN75 = 36; // решение для поправки давления должно быть получено точно.
			if (1.0e-4*fabs(delta075)<epsilon75) {
				epsilon75 = 1.0e-4*fabs(delta075);
			}
			if (iflag175 == 1) {
				iflag75 = 1;
			}
			//printf("%e",epsilon); system("pause");
		}
	}

	
	if (iVar == TEMP) {
		maxit75 = 2000;
	}
	if (iVar == PAM) {
		maxit75 = 2000; // 2000
	}
	if ((iVar == VX) || (iVar == VY) || (iVar == VZ)) {
		maxit75 = 100;//100
	}
	if (iVar == TOTALDEFORMATIONVAR) {
		maxit75 = 2000; // 2000
		//if (1.0e-4*fabs(delta075) < epsilon75) {
			//epsilon75 = 1.0e-4*fabs(delta075);
		//}
		epsilon75 = 1.0e-16;
		iN75 = 100;
		if (iflag175 == 1) {
			iflag75 = 1;
		}

	}

	
	

	


	// Мы обязательно должны сделать несколько итераций. (не менее 10).
	// Если только решение не удовлетворяет уравнению тождественно.
	while (((icount75 < iN75) && (iflag175 != 0)) || (iflag75 != 0 && icount75 < maxit75)) {

		// 6.01.2017: Body BiCGStab + AMG. (BiCGStab_internal4).


		icount75++;

		count_iter_for_film_coef75++;
		// В случае задачи Ньютона - Рихмана, Стефана-Больцмана и миксового условия не итерируем до конца обрываем, 
		// т.к. нам требуется частая пересборка матрицы. 13 марта 2016.
		//if (((adiabatic_vs_heat_transfer_coeff > ADIABATIC_WALL_BC) || (breakRUMBAcalc_for_nonlinear_boundary_condition)) && (count_iter_for_film_coef75>5)) break;

		roi75 = Scal(roc75, ri75, n75);
		bet75 = (roi75 / roim175)*(al75 / wi75);


		//printf("%e %e %e %e\n",roi75,roim175,al75,wi75);
		//system("pause");

#pragma omp parallel for 
		for (i75 = 0; i75<n75; i75++) {
			doublereal pibuf75 = ri75[i75] + (pi75[i75] - vi75[i75] * wi75)*bet75;
			pi75[i75] = pibuf75;
		}

		// Первое предобуславливание.
		// Ky=pi
#pragma omp parallel for
		for (i75 = 0; i75 < n75; i75++) {
			y75[i75] = 0.0; // Если начинать не с нуля то небудет сходимости для PAM !.
			y76[i75+1 ] = 0.0;//+1
			pi76[i75+1] = pi75[i75];//+1
		}


		// multigrid Ruge and Stuben preconditioning [1986].
		// достаточно одного V цикла.
		// K*y76 = pi76;
		ifirst[0] = 11; // Нулевое начальное приближение
	for (integer i_numberV_cycle = 0; i_numberV_cycle < 1; i_numberV_cycle++) {
		ncyc[0] = 1011; // достаточно одного V цикла.
		solve_(madapt, ncyc, nrd, nsolco, nru, &kout, ierr, &a[1], &y76[1], &pi76[1], &
			ia[1], &ja[1], &ig[1], eps, imin, imax, iminw, imaxw, &ig[icgst],
			&ig[icgst + ndicg], nstcol, iarr, time, &ncyc0, &irow0, &levels,
			nda, ndja, ndu, ndf, &mda, &mdja, &mdu, &mdf, &iup, &ium, resi, &
			res0, &res);
		//system("pause");
	}


	// Возвращение результата.
#pragma omp parallel for
	for (i75 = 0; i75 < n75; i75++) {
		y75[i75] = y76[i75 +1];//+1
	}

	MatrixCRSByVector(val75, col_ind75, row_ptr75, y75, vi75, n75); // vi==A*y;

	if ((fabs(roi75)<1e-30) && (fabs(Scal(roc75, vi75, n75))<1e-30)) {
		al75 = 1.0;
	}
	else if (fabs(roi75)<1e-30) {
		al75 = 0.0;
	}
	else {
		al75 = roi75 / Scal(roc75, vi75, n75);
	}


#pragma omp parallel for
	for (i75 = 0; i75<n75; i75++) {
		s75[i75] = ri75[i75] - al75*vi75[i75];
	}

	// Второе предобуславливание.
	// Kz=s

#pragma omp parallel for
	for (i75 = 0; i75<n75; i75++) z75[i75] = 0.0; // Если начинать не с нуля то небудет сходимости для PAM !.

#pragma omp parallel for
	for (i75 = 0; i75 < n75; i75++) {
		vec75[i75] = s75[i75];
		z76[i75 +1] = 0.0;//+1
		s76[i75 +1] = s75[i75];//+1
	}

	// multigrid Ruge and Stuben preconditioning [1986].
	// достаточно одного V цикла.
	// K*z76 = s76;
	ifirst[0] = 11; // Нулевое начальное приближение
	for (integer i_numberV_cycle = 0; i_numberV_cycle < 1; i_numberV_cycle++) {
		ncyc[0] = 1011; // достаточно одного V цикла.
		solve_(madapt, ncyc, nrd, nsolco, nru, &kout, ierr, &a[1], &z76[1], &s76[1], &
			ia[1], &ja[1], &ig[1], eps, imin, imax, iminw, imaxw, &ig[icgst],
			&ig[icgst + ndicg], nstcol, iarr, time, &ncyc0, &irow0, &levels,
			nda, ndja, ndu, ndf, &mda, &mdja, &mdu, &mdf, &iup, &ium, resi, &
			res0, &res);
		//system("pause");
	}

#pragma omp parallel for
	for (i75 = 0; i75 < n75; i75++) {
		s75[i75] = vec75[i75];
		// Возвращаем результат.
		z75[i75] = z76[i75 +1];//+1
	}

	MatrixCRSByVector(val75, col_ind75, row_ptr75, z75, t75, n75); // t==A*z;

	wi75 = Scal(t75, s75, n75) / Scal(t75, t75, n75);
	// printf("%e %e",Scal(t75,s75,n75),Scal(t75,t75,n75));

#pragma omp parallel for
	for (i75 = 0; i75<n75; i75++) {
		//dx75[i75]+=al75*pi75[i75]+wi75*s75[i75]; // так было без предобуславливателя
		dx75[i75] += al75*y75[i75] + wi75*z75[i75]; // так стало с предобуславливателем
		ri75[i75] = s75[i75] - wi75*t75[i75];
	}
	deltai75 = NormaV(ri75, n75);

	//printf("deltai75=%e\n",deltai75); system("pause");

	// печать невязки на консоль
	bool bprint_mesage_diagnostic = true;
	if (bprint_mesage_diagnostic) {
		if ((icount75 % 10) == 0) {
			std::cout << "iter  residual" << std::endl;
			
		}
		std::cout << icount75 << " " << deltai75 << std::endl; 
	}

	// 28.07.2016.
	//std::cout << icount75 << " " << deltai75 << std::endl;

	//system("pause");
	if (deltai75 > delta_old_iter75) i_signal_break_pam_opening75++;
	delta_old_iter75 = deltai75;
	if (iVar == PAM) {
		if (i_signal_break_pam_opening75 > i_limit_signal_pam_break_opening75) {
			// досрочный выход из цикла.
			std::cout << "icount PAM=" <<  icount75 << std::endl;
			break;
		}
	}

	if (deltai75 <epsilon75) iflag75 = 0; // конец вычисления
	else {
		/*
		// 02.05.2018
		if (bSIMPLErun_now_for_temperature) {
			// Эти значения невязок для CFD задач были 
			// успешно опробованы на задаче теплового расчёта
			// радиатора водяного охлаждения 3л/мин (совместное решение cfd + temperature 
			// + приближение Обербека-Буссинеска.).
			switch (iVar) {
			case VX: if (deltai75/ delta075 < 1e-2) iflag75 = 0;   break; //5e-5
			case VY: if (deltai75/ delta075 < 1e-2) iflag75 = 0;   break; // 5e-5
			case VZ: if (deltai75/ delta075 < 1e-2) iflag75 = 0;   break; // 5e-5
			//case TEMP: tolerance = 1e-8;  break; // 1e-7
			case PAM: if (deltai75/ delta075 < 1e-1) iflag75 = 0;  break; // 1e-6
			}
		}
		else {
			if (iflag75 != 0) {
			*/
			roim175 = roi75;
			//}
		//}
	}

	if (iVar == TEMP) {
#if doubleintprecision == 1
		//printf("epsilon=%e deltai=%e icount=%lld\n",epsilon75,deltai75, icount75);
#else
		//printf("epsilon=%e deltai=%e icount=%d\n",epsilon75,deltai75, icount75);
#endif

		//system("pause");
	}

	//icount_V_cycle = icount75; // количество итераций в BiCGStabP для лога.

	//if (icount75 > 2600) break; // 15.02.2017

	}

	// Возвращение результата вычислений.
#pragma omp parallel for
	for (i75 = 0; i75 < n75; i75++) {
		u[i75 + 1] = dx75[i75];
		//x_best_search[i75 + 1] = dx75[i75];
	}

	// Освобождение оперативной памяти.
	// Первое предобуславливание
	if (pi76 != NULL) {
		delete[] pi76;
		pi76 = NULL;
	}
	if (y76 != NULL) {
		delete[] y76;
		y76 = NULL;
	}
	// Второе предобуславливание
	if (z76 != NULL) {
		delete[] z76;
		z76 = NULL;
	}
	if (s76 != NULL) {
		delete[] s76;
		s76 = NULL;
	}
	if (ri75 != NULL) {
		delete[] ri75;
		ri75 = NULL;
	}
	if (roc75 != NULL) {
		delete[] roc75;
		roc75 = NULL;
	}
	if (s75 != NULL) {
		delete[] s75;
		s75 = NULL;
	}
	if (t75 != NULL) {
		delete[] t75;
		t75 = NULL;
	}
	if (vec75 != NULL) {
		delete[] vec75;
		vec75 = NULL;
	}
	if (vi75 != NULL) {
		delete[] vi75;
		vi75 = NULL;
	}
	if (pi75 != NULL) {
		delete[] pi75;
		pi75 = NULL;
	}
	if (dx75 != NULL) {
		delete[] dx75;
		dx75 = NULL;
	}
	if (dax75 != NULL) {
		delete[] dax75;
		dax75 = NULL;
	}
	if (y75 != NULL) {
		delete[] y75;
		y75 = NULL;
	}
	if (z75 != NULL) {
		delete[] z75;
		z75 = NULL;
	}

	// Освобождение оперативной памяти.
	if (val75 != NULL) {
		delete[] val75;
		val75 = NULL;
	}
	if (col_ind75 != NULL) {
		delete[] col_ind75;
		col_ind75 = NULL;
	}
	if (row_ptr75 != NULL) {
		delete[] row_ptr75;
		row_ptr75 = NULL;
	}

	
	
	if (debug_reshime) system("pause");
	//system("pause");

	if (*ierr > 0) {
		return 0;
	}
L10:
	wrkcnt_(&kout, &ia[1], &ig[1], imin, imax, iminw, &levels, time, &ncyc0, &
		iup, &mda, &mdia, &mdja, &mdu, &mdf, &mdig, &res0, &res);
	return 0;
L60:
	//io___29.ciunit = ium;
	//s_wsfe(&io___29);
	//e_wsfe();
	printf("*** ERROR IN AMG1R5: NDIG TOO SMALL ***\n");

	*ierr = 6;
	return 0;
L70:
	//io___30.ciunit = ium;
	//s_wsfe(&io___30);
	//e_wsfe();

	printf("*** ERROR IN AMG1R5: ILLEGAL PARAMETER LEVELX ***\n");

	*ierr = 18;
	return 0;
} /* amg1r5_Vorst_modification */


  /*     Last change:  ERB  22 Aug 2000   10:31 am */
  /* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */

  /*     AMG1R5                                        MAIN SUBROUTINE */

  /*     RELEASE 1.5, OCTOBER 1990 */

  /*     CHANGES AGAINST VERSION 1.1, JULY 1985: */

  /* 1.  A BUG WAS DETECTED WHICH UNDER CERTAIN CIRCUMSTANCES INFLUENCED */
  /*     SLIGHTLY THE CONVERGENCE RATE OF AMG1R1. FOR THAT REASON, THE */
  /*     FOLLOWING LINE IN SUBROUTINE RESC: */
  /*     IW(IMAXW(KC-1)+1) = IA(IMIN(KC)) */
  /*     HAS BEEN CHANGED TO: */
  /*     IW(IMAXW(KC-1)+1) = IAUX */

  /* 2.  A BUG WAS DETECTED IN SUBROUTINE PWINT. UNDER CERTAIN CIRCUM- */
  /*     STANCES AN UNDEFINED VARIABLE WAS USED. ALTHOUGH THIS DID NOT */
  /*     AFFECT THE NUMERICAL RESULTS, PROBLEMS CAN OCCUR IF CHECKING */
  /*     FOR UNDEFINED VARIABLES IS USED. TO FIX THIS ERROR, IN PWINT */
  /*     THE LABEL 1000 WAS MOVED TO THE STATEMENT */
  /*     IBLCK1 = IMINW(K). */

  /* 3.  A PARAMETER LRATIO HAS BEEN INTRODUCED, DENOTING THE RATIO */
  /*     OF SPACE OCCUPIED BY A DOUBLE PRECISION REAL VARIABLE AND */
  /*     THAT OF AN INTEGER. FOR THE IBM-VERSION LRATIO HAS BEEN SET */
  /*     TO 2. CHANGE THIS VALUE IF NECESSARY. (IF, FOR EXAMPLE, YOU */
  /*     WANT TO CHANGE THE DOUBLE PRECISION VECTORS TO SINGLE PRE- */
  /*     CISION, LRATIO HAS TO BE SET TO 1. IN THE YALE SMP - ROUTINE */
  /*     NDRV THERE IS A PARAMETER LRATIO, TOO. */

  /* 4.  TYPE DECLARATIONS REAL*4 AND REAL*8 HAVE BEEN CHANGED TO THE */
  /*     STANDARD-CONFORMING KEYWORDS REAL AND DOUBLE PRECISION, RESPEC- */
  /*     TIVELY. */

  /* 5.  CALLS TO THE FOLLOWING INTRINSIC FUNCTIONS HAVE BEEN REPLACED BY */
  /*     CALLS USING GENERIC NAMES: DSQRT, MIN0, MAX0, IABS, DABS, FLOAT, */
  /*     DFLOAT, DMAX1, ISIGN, IDINT, DLOG10. */

  /* 6.  A SAVE STATEMENT HAS BEEN INSERTED IN ALL SUBROUTINES. */

  /* 7.  EXTERNAL DECLARATION STATEMENTS HAVE BEEN INSERTED IN ALL SUB- */
  /*     ROUTINES FOR ALL EXTERNAL REFERENCES. */

  /* ----------------------------------------------------------------------- */

  /*     CHANGE AGAINST VERSION 1.3, APRIL 1986: */

  /* 1.  A BUG IN SUBROUTINE CHECK HAS BEEN REMOVED. IF THE ORIGINAL MATRIX */
  /*     WAS STORED IN AN UNSYMMETRIC WAY, THE SYMMETRIZATION BY AMG1R3 */
  /*     COULD FAIL UNDER CERTAIN CIRCUMSTANCES. FOR A FIX, THE FOLLOWING */
  /*     STATEMENTS IN SUBROUTINE CHECK HAVE BEEN CHANGED: */

  /*     DO 450 J=IA(I)+1,IA(I+1)-1 WAS CHANGED TO */
  /*     DO 450 J=IA(I)+1,ICG(I)-1 */

  /*     DO 430 J1=IA(I1)+1,IA(I1+1)-1 WAS CHANGED TO */
  /*     DO 430 J1=IA(I1)+1,ICG(I1)-1 */

  /*     DO 550 J=IA(I)+1,IA(I+1)-1 WAS CHANGED TO */
  /*     DO 550 J=IA(I)+1,ICG(I)-1 */

  /*     DO 530 J1=IA(I1)+1,IA(I1+1)-1 WAS CHANGED TO */
  /*     DO 530 J1=IA(I1)+1,ICG(I1)-1 */

  /* 2.  THE EXPLANATORY PART IN SUBROUTINE AMG1R5 HAS BEEN ENLARGED TO */
  /*     AVOID MISUNDERSTANDINGS IN THE DEFINITION OF THE ARGUMENT LIST. */

  /* ----------------------------------------------------------------------- */

  /*     CHANGE AGAINST VERSION 1.4, OCTOBER, 1990 (BY JOHN W. RUGE) */

  /* 1.  A BUG IN SUBROUTINE CHECK HAS BEEN REMOVED. IF THE ORIGINAL MATRIX */
  /*     WAS STORED IN AN UNSYMMETRIC WAY, THE SYMMETRIZATION BY AMG1R3 */
  /*     COULD STILL FAIL UNDER CERTAIN CIRCUMSTANCES, AND WAS NOT FIXED */
  /*     IN THE PREVIOUS VERSION. IN ADDITION, THE ROUTINE WAS CHANGED */
  /*     IN ORDER TO AVOID SOME UNNECESSARY ROW SEARCHES FOR TRANSOSE */
  /*     ENTRIES. */

  /* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */

/* Subroutine */ integer amg1r5_Vorst_modification_matrix_Assemble2(doublereal *a, integer *ia, integer *ja,
	doublereal *u, doublereal *f, integer *ig, integer *nda, integer *
	ndia, integer *ndja, integer *ndu, integer *ndf, integer *ndig,
	integer *nnu, integer *matrix, integer *iswtch, integer *iout,
	integer *iprint, integer *levelx, integer *ifirst, integer *ncyc,
	doublereal *eps, integer *madapt, integer *nrd, integer *nsolco,
	integer *nru, doublereal *ecg1, doublereal *ecg2, doublereal *ewt2,
	integer *nwt, integer *ntr, integer *ierr, SIMPLESPARSE &sparseM, integer n)
{

	// 23-24 декабря 2017.
	// В данной функции amg1r5 алгоритм является предобуславливателем для алгоритма 
	// Хенка Ван дер Ворста BiCGStab[1992].

	// sl, slb - Матрица СЛАУ для алгоритма BiCGStab [1992] Хенка Ван дер Ворста.
	// Матрица СЛАУ (a, ia, ja) модифицируется в процессе работы алгебраического многосеточного метода.
	// априори: nnu[0]==maxelm+maxbound
	// ndu[0]==ndf[0] ...

	/* Format strings */


	/* Builtin functions */


	/* Local variables */
	integer iVar = TEMP;
	integer mda = 0, mdf = 0, mdu = 0;
	doublereal res = 0.0;
	integer ium = 0, iup = 0;
	doublereal res0 = 0.0;
	extern /* Subroutine */ integer idec_(integer *, integer *, integer *,
		integer *);
	integer mdia = 0, mdja = 0, mdig = 0, iarr[25] = { 0 };
	//static real time[20];
	unsigned int time[20] = { 0 };
	integer imin[25] = { 0 }, imax[25] = { 0 };
	doublereal resi[25] = { 0.0 };
	integer kout = 0, ncyc0 = 0, irow0 = 0, ndicg = 0, icgst = 0, iminw[25] = { 0 }, imaxw[25] = { 0 };
	extern /* Subroutine */ integer first_(integer *, doublereal *, integer *,
		integer *, integer *, integer *), solve_(integer *, integer *,
			integer *, integer *, integer *, integer *, integer *, doublereal
			*, doublereal *, doublereal *, integer *, integer *, integer *,
			doublereal *, integer *, integer *, integer *, integer *, integer
			*, integer *, integer *, integer *, /*real*/ unsigned int *, integer *, integer *,
			integer *, integer *, integer *, integer *, integer *, integer *,
			integer *, integer *, integer *, integer *, integer *, doublereal
			*, doublereal *, doublereal *), setup_(integer *, integer *,
				integer *, doublereal *, doublereal *, doublereal *, integer *,
				integer *, integer *, doublereal *, doublereal *, integer *,
				integer *, integer *, integer *, integer *, integer *, integer *,
				integer *, integer *, integer *, integer *, /*real*/ unsigned int *, integer *,
				integer *, integer *, integer *, integer *, integer *, integer *,
				integer *, integer *, integer *, integer *, integer *, integer *,
				integer *, integer *, integer *, integer *);
	integer ndigit = 0, levels = 0, kevelx = 0, nstcol[25] = { 0 }, kswtch = 0;
	extern /* Subroutine */ integer wrkcnt_(integer *, integer *, integer *,
		integer *, integer *, integer *, integer *, /*real*/ unsigned int *, integer *,
		integer *, integer *, integer *, integer *, integer *, integer *,
		integer *, doublereal *, doublereal *);

	/* Fortran I/O blocks */

	bool debug_reshime = false;
	integer count_iter_for_film_coef75 = 0;
	doublereal delta_old_iter75 = 1.0e10;

	// Если число расходимостей превысит оговорённую константу то произойдёт выход из алгоритма.
	integer i_signal_break_pam_opening75 = 0;
	// x хорошее значение.
	const integer i_limit_signal_pam_break_opening75 = 4000;//20

	integer maxit75 = 2000;
	integer iN75 = 10;

	integer iflag75 = 1, icount75 = 0;
	doublereal delta075 = 1.0e30, deltai75 = 1.0e30;
	doublereal bet75 = 0.0, roi75 = 0.0;
	doublereal roim175 = 1.0, al75 = 1.0, wi75 = 1.0;

	doublereal epsilon75 = dterminatedTResudual;  // точность вычисления
	integer i75 = 0;
	integer iflag175 = 1;

	// Вектора необходимые для работы BiCGStab.
	doublereal* ri75 = NULL;
	doublereal* roc75 = NULL;
	doublereal* s75 = NULL;
	doublereal* t75 = NULL;
	doublereal* vec75 = NULL;
	doublereal* vi75 = NULL;
	doublereal* pi75 = NULL;
	doublereal* dx75 = NULL;
	doublereal* dax75 = NULL;
	doublereal* y75 = NULL;
	doublereal* z75 = NULL;
	// Первое предобуславливание:
	doublereal* y76 = NULL;
	doublereal* pi76 = NULL;
	// Второе предобуславливание:
	doublereal* z76 = NULL;
	doublereal* s76 = NULL;

	// нумерация векторов начинается с нуля.
	integer n75 = n; // число неизвестных на подробном уровне.
	doublereal* val75 = NULL;
	//val75 = new doublereal[nnz];
	integer* col_ind75 = NULL;
	//col_ind75 = new integer[nnz];
	integer* row_ptr75 = NULL;
	//row_ptr75 = new integer[n75 + 1];


	/*         ----------------------------------------------- */
	/*         | AMG-MODULE FOR SOLVING LINEAR SYSTEMS L*U=F | */
	/*         ----------------------------------------------- */

	/*         -------------------------------------------------------------- */

	/*     ASSUMPTIONS ON L: */

	/*         THE PROGRAM REQUIRES: */ // Требования к входным данным:

										/*             - DIAGONAL ENTRIES ARE ALWAYS POSITIVE (ON ALL GRIDS); */
										/*             - L IS A SQUARE MATRIX WHICH IS EITHER REGULAR OR SINGULAR */
										/*               WITH ROWSUMS=0. */
										// Матрица L это квадратная матрица с всегда положительными диагональными элементами с нулевой суммой коэффициентов в каждой строке (положительная определённость).

										/*         FOR THEORETICAL REASONS THE FOLLOWING SHOULD HOLD: */

										/*             - L POSITIVE DEFINITE (OR SEMI-DEFINITE WITH ROWSUM=0) */
										/*             - L "ESSENTIALLY" POSITIVE TYPE, I.E., */

										/*                  -- DIAGONAL ENTRIES MUST BE > 0 ; */
										/*                  -- MOST OF THE OFF-DIAGONAL ENTRIES <= 0 ; */
										/*                  -- ROWSUMS SHOULD BE >= 0 . */

										// Матрица L положительно определённая: диагональные элементы  строго больше нуля, большинство внедиагональных элементов <= 0.
										// Сумма коэффициентов в строке больше либо равна нулю - диагональное преобладание.


										/*     THE USER HAS TO PROVIDE THE MATRIX L, THE RIGHT HAND SIDE F AND */
										/*     CERTAIN POINTER VECTORS IA AND JA. */

										// Пользователь задаёт матрицу коэффициентов L, правую часть F, а также информацию о связях между элементами матрицы в специальных столбцах IA и JA.
										// IA - позиция первого элемента в строке. JA - номер столбца для каждого элемента, первым записывается диагональный элемент.
										// В фортране нумерация начинается с единицы.

										/*         -------------------------------------------------------------- */

										/*     STORAGE OF L: */ // Требования к хранению матрицы L.

																/*         THE NON-ZERO ENTRIES OF THE MATRIX L ARE STORED IN */
																/*         "COMPRESSED" SKY-LINE FASHION IN A 1-D VECTOR A, I.E., ROW */
																/*         AFTER ROW, EACH ROW STARTING WITH ITS DIAGONAL ELEMENT. THE */
																/*         OTHER NON-ZERO ROW ENTRIES FOLLOW THEIR DIAGONAL ENTRY IN ANY */
																/*         ORDER. */


																/*         IN ORDER TO IDENTIFY EACH ELEMENT IN A, THE USER HAS TO */
																/*         PROVIDE TWO POINTER ARRAYS IA AND JA. IF NNU DENOTES THE TOTAL */
																/*         NUMBER OF UNKNOWNS, THE NON-ZERO ENTRIES OF ANY ROW I OF L */
																/*         (1.LE.I.LE.NNU) ARE STORED IN A(J) WHERE THE RANGE OF J */
																/*         IS GIVEN BY */

																/*                     IA(I) .LE. J .LE. IA(I+1)-1. */

																/*         THUS, IA(I) POINTS TO THE POSITION OF THE DIAGONAL ENTRY OF */
																/*         ROW I WITHIN THE VECTOR A. IN PARTICULAR, */

																/*                     IA(1) = 1 ,  IA(NNU+1) = 1 + NNA */

																/*         WHERE NNA DENOTES THE TOTAL NUMBER OF MATRIX ENTRIES STORED. */
																/*         THE POINTER VECTOR JA HAS TO BE DEFINED SUCH THAT */
																/*         ANY ENTRY A(J) CORRESPONDS TO THE UNKNOWN U(JA(J)), I.E., */
																/*         JA(J) POINTS TO THE COLUMN INDEX OF A(J). */
																/*         IN PARTICULAR, A(IA(I)) IS THE DIAGONAL ENTRY OF ROW I */
																/*         AND CORRESPONDS TO THE UNKNOWN U(I): JA(IA(I))=I. */

																/*         IN THIS TERMINOLOGY, THE I-TH EQUATION READS AS FOLLOWS */
																/*         (FOR ANY I WITH  1.LE.I.LE.NNU): */

																/*                  F(I) =        SUM      A(J) * U(JA(J)) */
																/*                           J1.LE.J.LE.J2 */

																/*         WHERE F(I) DENOTES THE I-TH COMPONENT OF THE RIGHT HAND */
																/*         SIDE AND */

																/*                     J1 = IA(I) ,  J2 = IA(I+1)-1. */

																/*         NOTES: THE ENTRY IA(NNU+1) HAS TO TOCHKA TO THE FIRST FREE */
																/*                ENTRY IN VECTORS A AND JA, RESPECTIVELY. OTHERWISE, */
																/*                AMG CANNOT KNOW THE LENGTH OF THE LAST MATRIX ROW. */

																/*                THE INPUT VECTORS A, IA AND JA ARE CHANGED BY AMG1R5. */
																/*                SO, AFTER RETURN FROM AMG1R5, THE PACKAGE MUST NOT */
																/*                BE CALLED A SECOND TIME WITHOUT HAVING NEWLY DEFINED */
																/*                THE INPUT VECTORS AND USING ISWTCH=4. OTHERWISE, THE */
																/*                SETUP PHASE WILL FAIL. */
																/*                  ON THE OTHER HAND, RUNNING AMG A SECOND TIME ON THE */
																/*                SAME INPUT DATA WITH ISWTCH=4 HAS NO SENSE, BECAUSE */
																/*                THE RESULTS OF THE FIRST SETUP PHASE ARE STILL STORED */
																/*                AND THUS THIS PHASE CAN BE SKIPPED IN A SECOND CALL. */
																/*                IN ORDER TO DO THIS, SET ISWTCH TO 1, 2 OR 3. */

																/* ----------------------------------------------------------------------- */

																/*         THE FORM OF THE CALLING PROGRAM HAS TO BE AS FOLLOWS: */

																/*               PROGRAM DRIVER */
																/*         C */
																/*               DOUBLE PRECISION A(#NDA),U(#NDU),F(#NDF) */
																/*               INTEGER IA(#NDIA),JA(#NDJA),IG(#NDIG) */
																/*         C */
																/*               NDA  = #NDA */
																/*               NDU  = #NDU */
																/*               NDF  = #NDF */
																/*               NDIA = #NDIA */
																/*               NDJA = #NDJA */
																/*               NDIG = #NDIG */
																/*         C */
																/*         C     SET UP A, F, IA, JA AND SPECIFY NECESSARY PARAMETERS */
																/*         C */
																/*               .... */
																/*               .... */
																/*         C */
																/*               CALL AMG1R5(A,IA,JA,U,F,IG, */
																/*        +                  NDA,NDIA,NDJA,NDU,NDF,NDIG,NNU,MATRIX, */
																/*        +                  ISWTCH,IOUT,IPRINT, */
																/*        +                LEVELX,IFIRST,NCYC,EPS,MADAPT,NRD,NSOLCO,NRU, */
																/*        +                  ECG1,ECG2,EWT2,NWT,NTR, */
																/*        +                  IERR) */
																/*         C */
																/*               .... */
																/*               .... */
																/*         C */
																/*               STOP */
																/*               END */

																/* ----------------------------------------------------------------------- */

																/*     INPUT VIA ARRAYS (SEE ABOVE): */

																/*     A        -   MATRIX L */

																/*     IA       -   POINTER VECTOR */

																/*     JA       -   POINTER VECTOR */

																/*     U        -   FIRST APPROXIMATION TO SOLUTION */

																/*     F        -   RIGHT HAND SIDE */


																/* ----------------------------------------------------------------------- */


																/*     SCALAR INPUT PARAMETERS OF AMG1R5: */

																/*     THE INPUT PARAMETERS OF AMG1R5 IN THE LIST BELOW ARE ARRANGED */
																/*     ACCORDING TO THEIR IMPORTANCE TO THE GENERAL USER. THE PARAMETERS */
																/*     PRECEEDED BY A * MUST BE SPECIFIED EXPLICITELY. ALL THE OTHER */
																/*     PARAMETERS ARE SET TO STANDARD VALUES IF ZERO ON INPUT. */

																/*     THERE ARE FOUR CLASSES OF INPUT PARAMETERS WITH DECREASING PRI- */
																/*     ORITY: */

																/*     1. PARAMETERS DESCRIBING THE USER-DEFINED PROBLEM AND DIMENSIONING */
																/*        OF VECTORS IN THE CALLING PROGRAM */

																/*     2. PARAMETERS SPECIFYING SOME GENERAL ALGORITHMIC ALTERNATIVES AND */
																/*        THE AMOUNT OF OUTPUT DURING SOLUTION */

																/*     3. PARAMETERS CONTROLLING THE MULTIGRID CYCLING DURING THE SOLU- */
																/*        TION PHASE */

																/*     4. PARAMETERS CONTROLLING THE CREATION OF COARSER GRIDS AND INTER- */
																/*        POLATION FORMULAS. */

																/*     ONLY THE CLASS 1 - PARAMETERS MUST BE SPECIFIED EXPLICITELY BY */
																/*     THE USER. CLASS 2 - PARAMETERS CONTROL THE GENERAL PERFORMANCE OF */
																/*     AMG1R5. CHANGING THEM DOESN'T REQUIRE UNDERSTANDING THE AMG - */
																/*     ALGORITHM. SPECIFYING NON-STANDARD-VALUES FOR CLASS 3 - PARAMETERS */
																/*     PRESUPPOSES A GENERAL KNOWLEDGE OF MULTIGRID METHODS, WHEREAS THE */
																/*     FUNCTION OF CLASS 4 - PARAMETERS IS ONLY UNDERSTANDABLE AFTER */
																/*     STUDYING THE AMG-ALGORITHM IN DETAIL. FORTUNATELY IN MOST CASES */
																/*     THE CHOICE OF CLASS 3 AND 4 - PARAMETERS ISN'T CRITICAL AND USING */
																/*     THE AMG1R5 - SUPPLIED STANDARD VALUES SHOULD GIVE SATISFACTORY */
																/*     RESULTS. */

																/*         -------------------------------------------------------------- */

																/*     CLASS 1 - PARAMETERS: */

																/*  *  NDA      -   DIMENSIONING OF VECTOR A IN CALLING PROGRAM */

																/*  *  NDIA     -   DIMENSIONING OF VECTOR IA IN CALLING PROGRAM */

																/*  *  NDJA     -   DIMENSIONING OF VECTOR JA IN CALLING PROGRAM */

																/*  *  NDU      -   DIMENSIONING OF VECTOR U IN CALLING PROGRAM */

																/*  *  NDF      -   DIMENSIONING OF VECTOR F IN CALLING PROGRAM */

																/*  *  NDIG     -   DIMENSIONING OF VECTOR IG IN CALLING PROGRAM */

																/*  *  NNU      -   NUMBER OF UNKNOWNS */

																/*  *  MATRIX   -   INTEGER VALUE CONTAINING INFO ABOUT THE MATRIX L. */

																/*                  1ST DIGIT OF MATRIX  --  ISYM: */
																/*                    =1: L IS SYMMETRIC; */
																/*                    =2: L IS NOT SYMMETRIC. */

																/*                  2ND DIGIT OF MATRIX  --  IROW0: */
																/*                    =1: L HAS ROWSUM ZERO; */
																/*                    =2: L DOES NOT HAVE ROWSUM ZERO. */

																/*         -------------------------------------------------------------- */

																/*     CLASS 2 - PARAMETERS: */

																/*     ISWTCH   -   PARAMETER CONTROLLING WHICH MODULES OF AMG1R5 ARE TO */
																/*                  BE USED. */
																/*                    =1:   CALL FOR -----, -----, -----, WRKCNT. */
																/*                    =2:   CALL FOR -----, -----, SOLVE, WRKCNT. */
																/*                    =3:   CALL FOR -----, FIRST, SOLVE, WRKCNT. */
																/*                    =4:   CALL FOR SETUP, FIRST, SOLVE, WRKCNT. */
																/*                  SETUP DEFINES THE OPERATORS NEEDED IN THE SOLUTION */
																/*                         PHASE. */
																/*                  FIRST INITIALIZES THE SOLUTION VECTOR (SEE PARAMETER */
																/*                         IFIRST). */
																/*                  SOLVE COMPUTES THE SOLUTION BY AMG CYCLING (SEE */
																/*                         PARAMETER NCYC). */
																/*                  WRKCNT PROVIDES THE USER WITH INFORMATION ABOUT */
																/*                         RESIDUALS, STORAGE REQUIREMENTS AND CP-TIMES */
																/*                         (SEE PARAMETER IOUT). */
																/*                  IF AMG1R5 IS CALLED THE FIRST TIME, ISWTCH HAS TO */
																/*                  BE =4. INDEPENDENT OF ISWTCH, SINGLE MODULES CAN BE */
																/*                  BYPASSED BY A PROPER CHOICE OF THE CORRESPONDING */
																/*                  PARAMETER. */

																/*     IOUT     -   PARAMETER CONTROLLING THE AMOUNT OF OUTPUT DURING */
																/*                  SOLUTION PHASE: */

																/*                  1ST DIGIT: NOT USED; HAS TO BE NON-ZERO. */

																/*                  2ND DIGIT: */
																/*                    =0: NO OUTPUT (EXCEPT FOR MESSAGES) */
																/*                    =1: RESIDUAL BEFORE AND AFTER SOLUTION PROCESS */
																/*                    =2: ADD.: STATISTICS ON CP-TIMES AND STORAGE REQUI- */
																/*                        REMENTS */
																/*                    =3: ADD.: RESIDUAL AFTER EACH AMG-CYCLE */

																/*     IPRINT   -   PARAMETER SPECIFYING THE FORTRAN UNIT NUMBERS FOR */
																/*                  OUTPUT: */

																/*                  1ST DIGIT: NOT USED; HAS TO BE NON-ZERO */

																/*                  2ND AND 3RD DIGIT  --  IUP: UNIT NUMBER FOR RESULTS */

																/*                  4TH AND 5TH DIGIT  --  IUM: UNIT NUMBER FOR MESSAGES */

																/*         -------------------------------------------------------------- */

																/*     CLASS 3 - PARAMETERS: */

																/*     LEVELX   -   MAXIMUM NUMBER OF MG-LEVELS TO BE CREATED (>=1). */

																/*     IFIRST   -   PARAMETER FOR FIRST APPROXIMATION. */

																/*                  1ST DIGIT OF IFIRST: NOT USED; HAS TO BE NON-ZERO. */

																/*                  2ND DIGIT OF IFIRST  --  ITYPU: */
																/*                    =0: NO SETTING OF FIRST APPROXIMATION, */
																/*                    =1: FIRST APPROXIMATION CONSTANT TO ZERO, */
																/*                    =2: FIRST APPROXIMATION CONSTANT TO ONE, */
																/*                    =3: FIRST APPROXIMATION IS RANDOM FUNCTION WITH */
																/*                        THE CONCRETE RANDOM SEQUENCE BEING DETERMINED */
																/*                        BY THE FOLLWING DIGITS. */

																/*                  REST OF IFIRST  --  RNDU: */
																/*                    DETERMINES THE CONCRETE RANDOM SEQUENCE USED IN */
																/*                    THE CASE ITYPU=3. (IFIRST=13 IS EQUIVALENT TO */
																/*                    IFIRST=1372815) */

																/*     NCYC     -   INTEGER PARAMETER DESCRIBING THE TYPE OF CYCLE TO BE */
																/*                  USED AND THE NUMBER OF CYCLES TO BE PERFORMED. */

																/*                  1ST DIGIT OF NCYC  --  IGAM: */
																/*                    =1: V -CYCLE, */
																/*                    =2: V*-CYCLE, */
																/*                    =3: F -CYCLE, */
																/*                    =4: W -CYCLE. */
																/*                  IF NCYC IS NEGATIV, THEN THE APPROXIMATION OF THE */
																/*                  PROBLEM ON THE SECOND FINEST GRID IS COMPUTED BY */
																/*                  IGAM V-CYCLES ON THAT PARTICULAR GRID. */

																/*                  2ND DIGIT OF NCYC  --  ICGR: */
																/*                    =0: NO CONJUGATE GRADIENT, */
																/*                    =1: CONJUGATE GRADIENT (ONLY FIRST STEP OF CG), */
																/*                    =2: CONJUGATE GRADIENT (FULL CG). */

																/*                  3RD DIGIT OF NCYC  --  ICONV: */
																/*                    CONVERGENCE CRITERION FOR THE USER-DEFINED PROBLEM */
																/*                    (FINEST GRID): */
																/*                    =1: PERFORM A FIXED NUMBER OF CYCLES AS GIVEN BY */
																/*                        NCYCLE (SEE BELOW) */
																/*                    =2: STOP, IF  ||RES|| < EPS */
																/*                    =3: STOP, IF  ||RES|| < EPS * |F| */
																/*                    =4: STOP, IF  ||RES|| < EPS * |U| * |DIAG| */
																/*                    WITH ||RES|| = L2-NORM OF RESIDUAL, */
																/*                           EPS     (SEE INPUT PARAMETER EPS) */
																/*                           |F|   = SUPREMUM NORM OF RIGHT HAND SIDE */
																/*                           |U|   = SUPREMUM NORM OF SOLUTION */
																/*                         |DIAG|  = MAXIMAL DIAGONAL ENTRY IN MATRIX L */
																/*                    NOTE THAT IN ANY CASE THE SOLUTION PROCESS STOPS */
																/*                    AFTER AT MOST NCYCLE CYCLES. */

																/*                  REST OF NCYC  --  NCYCLE: */
																/*                    MAXIMAL NUMBER OF CYCLES TO BE PERFORMED (>0) OR */
																/*                    NCYCLE=0: NO CYCLING. */

																/*     EPS      -   CONVERGENCE CRITERION FOR SOLUTION PROCESS: (SEE */
																/*                  PARAMETER NCYC). NOTE THAT NO MORE THAN NCYCLE CYCLES */
																/*                  ARE PERFORMED, REGARDLESS OF EPS. */

																/*     MADAPT   -   INTEGER VALUE SPECIFYING THE CHOICE OF COARSEST */
																/*                  GRID IN CYCLING: */

																/*                  1ST DIGIT OF MADAPT  --  MSEL: */
																/*                    =1: IN CYCLING, ALL GRIDS CONSTRUCTED IN THE SETUP */
																/*                        PHASE ARE USED WITHOUT CHECK. */
																/*                    =2: THE NUMBER OF GRIDS IS AUTOMATICALLY REDUCED */
																/*                        IF THE CONVERGENCE FACTOR ON THE COARSER GRIDS */
																/*                        IS FOUND TO BE LARGER THAN A GIVEN VALUE FAC */
																/*                        (SEE BELOW). */

																/*                  REST OF MADAPT  --  FAC */
																/*                        THE REST OF MADAPT DEFINES THE FRACTIONAL PART */
																/*                        OF A REAL NUMBER FAC BETWEEN 0.1 AND 0.99, E.G. */
																/*                        MADAPT=258 MEANS MSEL=2 AND FAC=0.58. IF MADAPT */
																/*                        CONSISTS OF ONLY ONE DIGIT, FAC IS SET TO 0.7 */
																/*                        BY DEFAULT. */


																/*     NRD      -   PARAMETER DESCRIBING RELAXATION (DOWNWARDS): */

																/*                  1ST DIGIT OF NRD: NOT USED; HAS TO BE NON-ZERO. */

																/*                  2ND DIGIT OF NRD  --  NRDX: */
																/*                    ACTUAL NUMBER OF SMOOTHING STEPS TO BE PERFORMED */
																/*                    THE TYPE OF WHICH IS GIVEN BY THE FOLLOWING DIGITS */

																/*                  FOLLOWING DIGITS  --  ARRAY NRDTYP: */
																/*                    =1: RELAXATION OVER THE F-POINTS ONLY */
																/*                    =2: FULL GS SWEEP */
																/*                    =3: RELAXATION OVER THE C-POINTS ONLY */
																/*                    =4: FULL MORE COLOR SWEEP, HIGHEST COLOR FIRST */

																/*     NSOLCO   -   PARAMETER CONTROLLING THE SOLUTION ON COARSEST GRID: */

																/*                  1ST DIGIT  --  NSC: */
																/*                    =1: GAUSS-SEIDEL METHOD */
																/*                    =2: DIRECT SOLVER (YALE SMP) */

																/*                  REST OF NSOLCO  --  NRCX: (ONLY IF NSC=1) */
																/*                  NUMBER OF GS SWEEPS ON COARSEST GRID (>=0). */
																/*                  IF NRCX=0, THEN AS MANY GS SWEEPS ARE PERFORMED */
																/*                  AS ARE NEEDED TO REDUCE THE RESIDUAL BY TWO ORDERS */
																/*                  OF MAGNITUDE. (MAXIMAL 100 RELAXATION SWEEPS) */

																/*     NRU      -   PARAMETER FOR RELAXATION (UPWARDS), ANALOGOUS TO NRD. */

																/*         -------------------------------------------------------------- */

																/*     CLASS 4 - PARAMETERS: */

																/*     ECG1,ECG2-   REAL PARAMETERS AFFECTING THE CREATION OF COARSER */
																/*     EWT2     -   GRIDS AND/OR THE DEFINITION OF THE INTERPOLATION. */
																/*                  THE CHOICE OF THESE PARAMETERS DEPENDS ON */
																/*                  THE ACTUAL AMG VERSION (SEE SUBROUTINE CRSNG) */

																/*     NWT      -   INTEGER PARAMETER AFFECTING THE CREATION OF COARSER */
																/*                  GRIDS AND/OR THE DEFINITION OF THE INTERPOLATION. */
																/*                  THE CHOICE OF THIS PARAMETER DEPENDS ON */
																/*                  THE ACTUAL AMG VERSION (SEE SUBROUTINE CRSNG) */

																/*     NTR      -   PARAMETER CONTROLLING COARSE-GRID OPERATOR TRUNCATION */
																/*                    =0: PAIRS OF ZEROES ARE REMOVED FROM COARSE GRID */
																/*                        OPERATORS */
																/*                    =1: NO COARSE-GRID OPERATOR TRUNCATION */


																/* ----------------------------------------------------------------------- */

																/*     OUTPUT: */

																/*     U        -   CONTAINS THE COMPUTED SOLUTION */


																/*     IERR     -   ERROR PARAMETER: */

																/*                    >0: FATAL ERROR (ABNORMAL TERMINATION OF AMG1R5) */
																/*                    <0: NON-FATAL ERROR (EXECUTION OF AMG1R5 CONTINUES) */

																/*                  ERROR CODES IN DETAIL: */

																/*                  1. DIMENSIONING TOO SMALL FOR VECTOR */
																/*                        A      (IERR = 1) */
																/*                        IA     (IERR = 2) */
																/*                        JA     (IERR = 3) */
																/*                        U      (IERR = 4) */
																/*                        F      (IERR = 5) */
																/*                        IG     (IERR = 6) */

																/*                     NO YALE-SMP BECAUSE OF STORAGE (NDA TOO SMALL): */
																/*                               (IERR = -1) */
																/*                     NO YALE-SMP BECAUSE OF STORAGE (NDJA TOO SMALL): */
																/*                               (IERR = -3) */
																/*                     NO CG BECAUSE OF STORAGE (NDU TOO SMALL): */
																/*                               (IERR = -4) */
																/*                     NO SPACE FOR TRANSPOSE OF INTERPOLATION (NDA OR */
																/*                                                     NDJA TOO SMALL): */
																/*                               (IERR = -1) */

																/*                  2. INPUT DATA ERRONEOUS: */

																/*                     A-ENTRY MISSING, ISYM = 1:           (IERR = -11) */
																/*                     PARAMETER MATRIX MAY BE ERRONEOUS:   (IERR = -12) */
																/*                     DIAGONAL ELEMENT NOT STORED FIRST:   (IERR =  13) */
																/*                     DIAGONAL ELEMENT NOT POSITIV:        (IERR =  14) */
																/*                     POINTER IA ERRONEOUS:                (IERR =  15) */
																/*                     POINTER JA ERRONEOUS:                (IERR =  16) */
																/*                     PARAMETER ISWTCH ERRONEOUS:          (IERR =  17) */
																/*                     PARAMETER LEVELX ERRONEOUS:          (IERR =  18) */

																/*                  3. ERRORS OF THE AMG1R5-SYSTEM (SHOULD NOT OCCUR): */

																/*                     TRANSPOSE A-ENTRY MISSING:           (IERR =  21) */
																/*                     INTERPOLATION ENTRY MISSING:         (IERR =  22) */

																/*                  4. ALGORITHMIC ERRORS: */

																/*                     CG-CORRECTION NOT DEFINED:           (IERR =  31) */
																/*                     NO YALE-SMP BECAUSE OF ERROR IN */
																/*                     FACTORIZATION:                       (IERR = -32) */

																/* ----------------------------------------------------------------------- */

																/*     WORK SPACE: */

																/*     THE INTEGER VECTOR IG HAS TO BE PASSED TO AMG1R5 AS WORK SPACE. */

																/* ----------------------------------------------------------------------- */

																/*     DIMENSIONING OF INPUT VECTORS AND WORK SPACE: */

																/*     IT'S IMPOSSIBLE TO TELL IN ADVANCE THE EXACT STORAGE REQUIREMENTS */
																/*     OF AMG. THUS, THE FOLLOWING FORMULAS GIVE ONLY REASONABLE GUESSES */
																/*     FOR THE VECTOR LENGTHS WHICH HAVE TO BE DECLARED IN THE CALLING */
																/*     PROGRAM. IN THESE FORMULAS NNA DENOTES THE NUMBER OF NON-ZERO */
																/*     ENTRIES IN THE INPUT-MATRIX L AND NNU IS THE NUMBER OF UNKNOWNS. */

																/*     VECTOR         NEEDED LENGTH (GUESS) */
																/*       A               3*NNA + 5*NNU */
																/*       JA              3*NNA + 5*NNU */
																/*       IA              2.2*NNU */
																/*       U               2.2*NNU */
																/*       F               2.2*NNU */
																/*       IG              5.4*NNU */

																/* ----------------------------------------------------------------------- */


																/*     STANDARD CHOICES OF PARAMETERS (AS FAR AS MEANINGFUL): */

																/*          ISWTCH = 4 */
																/*          IOUT   = 12 */
																/*          IPRINT = 10606 */

																/*          LEVELX = 25 */
																/*          IFIRST = 13 */
																/*          NCYC   = 10110 */
																/*          EPS    = 1.D-12 */
																/*          MADAPT = 27 */
																/*          NRD    = 1131 */
																/*          NSOLCO = 110 */
																/*          NRU    = 1131 */

																/*          ECG1   = 0. */
																/*          ECG2   = 0.25 */
																/*          EWT2   = 0.35 */
																/*          NWT    = 2 */
																/*          NTR    = 0 */

																/*     IF ANY ONE OF THESE PARAMETERS IS 0 ON INPUT, ITS CORRESPONDING */
																/*     STANDARD VALUE IS USED BY AMG1R5. */

																/* ----------------------------------------------------------------------- */

																/*     PORTABILITY RESTRICTIONS: */

																/*     1. ROUTINE CTIME IS MACHINE DEPENDENT AND HAS TO BE ADAPTED TO */
																/*        YOUR COMPUTER INSTALLATION OR REPLACED BY A DUMMY ROUTINE. */

																/*     2. MOST INPUT PARAMETERS ARE COMPOSED OF SEVERAL DIGITS, THEIR */
																/*        SIGNIFICANCE HAVING BEEN DESCRIBED ABOVE. BE SURE NOT TO ENTER */
																/*        MORE DIGITS THAN YOUR COMPUTER CAN STORE ON AN INTEGER VARI- */
																/*        ABLE. */

																/*     3. APART FROM FORTRAN INTRINSIC FUNCTIONS AND SERVICE ROUTINES, */
																/*        THERE IS ONLY ONE EXTERNAL REFERENCE TO A PROGRAM NOT CONTAINED */
																/*        IN THE AMG1R5 - SYSTEM, I.E. THE LINEAR SYSTEM SOLVER NDRV OF */
																/*        THE YALE SPARSE MATRIX PACKAGE. IF YOU HAVN'T ACCESS TO THIS */
																/*        PACKAGE, ENTER A DUMMY ROUTINE NDRV AND AVOID CHOOSING NSC=2 */
																/*        (SUBPARAMETER OF NSOLCO). THEN NDRV ISN'T CALLED BY AMG1R5. */
																/*        IN THIS CASE, HOWEVER, INDEFINITE PROBLEMS WILL NOT BE SOLV- */
																/*        ABLE. */
																/*          THE YALE SPARSE MATRIX PACKAGE IS FREELY AVAILABLE FOR NON- */
																/*        PROFIT PURPOSES. CONTACT THE DEPARTMENT OF COMPUTER SCIENCE, */
																/*        YALE UNITVERSITY. */

																/*     4. IN AMG1R5 THERE IS THE PARAMETER LRATIO, DENOTING THE RATIO */
																/*        OF SPACE OCCUPIED BY A DOUBLE PRECISION REAL VARIABLE AND */
																/*        THAT OF AN INTEGER. FOR THE IBM-VERSION LRATIO HAS BEEN SET */
																/*        TO 2. CHANGE THIS VALUE IF NECESSARY. (THE SAME HAS TO BE */
																/*        DONE WITH THE YALE SMP-ROUTINE NDRV.) */


																/* ----------------------------------------------------------------------- */

																/*     AUTHORS: */

																/*          JOHN RUGE, FORT COLLINS (USA), */
																/*              INSTITUTE FOR COMPUTATIONAL STUDIES AT CSU; */

																/*          KLAUS STUEBEN, D-5205 ST. AUGUSTIN (W.-GERMANY), */
																/*              GESELLSCHAFT FUER MATHEMATIK UND DATENVERARBEITUNG (GMD). */

																/*          ROLF HEMPEL, D-5205 ST. AUGUSTIN (W.-GERMANY), */
																/*              GESELLSCHAFT FUER MATHEMATIK UND DATENVERARBEITUNG (GMD). */

																/* ----------------------------------------------------------------------- */


																/* ===> LRATIO HAS TO BE SET TO THE NUMBER OF INTEGERS OCCUPYING THE SAME */
																/*     AMOUNT OF STORAGE AS ONE DOUBLE PRECISION REAL. */


																/* ===> MAXGR IS THE MAXIMAL NUMBER OF GRIDS. CHANGING THIS UPPER LIMIT */
																/*     JUST REQUIRES CHANGING THE PARAMETER STATEMENT. */


																/* Parameter adjustments */
	--ig;
	--f;
	--u;
	--ja;
	--ia;
	--a;

	/* Function Body */
	*ierr = 0;

	/* ===> SET PARAMETERS TO STANDARD VALUES, IF NECCESSARY */


	if (*iout != 0) {
		idec_(iout, &c__2, &ndigit, iarr);
		kout = iarr[1];
	}
	else {
		kout = 2;
	}

	if (*iswtch != 0) {
		kswtch = *iswtch;
	}
	else {
		kswtch = 4;
	}

	if (*levelx > 0) {
		kevelx = my_imin(*levelx, 25);
	}
	else if (*levelx < 0) {
		goto L70;
	}
	else {
		kevelx = 25;
	}

	if (*iprint != 0) {
		idec_(iprint, &c__4, &ndigit, iarr);
		iup = iarr[1] * 10 + iarr[2];
		ium = iarr[3];
	}
	else {
		iup = 6;
		ium = 6;
	}
	icgst = *nnu + 3;
	ndicg = (*ndig - icgst + 1) / 2;
	if (ndicg <= 0) {
		goto L60;
	}



	switch (kswtch) {
	case 1:  goto L10;
	case 2:  goto L20;
	case 3:  goto L30;
	case 4:  goto L40;
	}
	//io___10.ciunit = ium;
	//s_wsfe(&io___10);
	//e_wsfe();


	printf("*** ERROR IN AMG1R5: ILLEGAL PARAMETER I.  SWTCH ***\n");

	*ierr = 17;
	return 0;

L40:
	setup_(nnu, matrix, &kevelx, ecg1, ecg2, ewt2, nwt, ntr, ierr, &a[1], &u[
		1], &ia[1], &ja[1], &ig[1], imin, imax, iminw, imaxw, &ig[icgst],
			&ig[icgst + ndicg], nstcol, iarr, time, &levels, &irow0, nda,
			ndia, ndja, ndu, ndf, &ndicg, &ium, &mda, &mdia, &mdja, &mdu, &
			mdf, &mdig, &c__25, &c__2);
	if (*ierr > 0) {
		return 0;
	}
L30:
	first_(ifirst, &u[1], imin, imax, iarr, &irow0);
L20:

	// Алгебраический Многосеточный Метод как предобуславливатель
	// к алгоритму Крыловского типа Хенка Ван Дер Ворста BiCGStab
	// со стабилизацией.
	// Требует ещё одну память под матрицу А на самом подробном уровне.
	// 5.01.2017 Алгоритм BiCGStab изобретён в 1992 году.

	/*
	// Для корректности работы надо передавать информацию о модифицированном размере в вызывающие функции вверх по коду.
	if ((a != NULL)&&(ja!=NULL)) {

	// На данном этапе доступен истинный размер матрицы СЛАУ - хранимое число ненулевых элементов.
	// Вычислим реальное число ненулевых элементов матрицы a. Если номер столбца равен нулю то число ненулевых элементов
	// в матрице заведомо меньше чем позиция этого нуля.
	// С помощью низкоуровневой операции realloc ужимаем размер матрицы a.
	// Последующие выделения оперативной памяти для внешнего Крыловского итерационного процесса BiCGStab не приведут к ещё большему расходу
	// оперативной памяти, т.к. мы её освободили.

	printf("nda=%lld\n", nda[0]);
	integer isize97 = 0;
	for (integer i_95 = 1; i_95 < ndja[0]; i_95++) {
	if (ja[i_95] == 0) {
	isize97 = i_95;
	if (i_95 + 2 < ndja[0] && ja[i_95 + 1] == 0 && ja[i_95 + 2] == 0) {
	break;
	}
	}
	}
	printf("nda_new=%lld\n", isize97);
	++a;
	++ja;
	a = (doublereal*)realloc(a, (static_cast<integer>(isize97)+2) * sizeof(doublereal));
	ja = (integer*)realloc(ja, (static_cast<integer>(isize97)+2) * sizeof(integer));
	--a;
	--ja;
	//system("pause");
	}
	*/

	ncyc[0] = 1011; // Предобуславливание в один V цикл.

					/*
					// для отладки.
					for (integer i_numberV_cycle = 0; i_numberV_cycle < 10; i_numberV_cycle++) {
					ncyc[0] = 1011;
					solve_(madapt, ncyc, nrd, nsolco, nru, &kout, ierr, &a[1], &u[1], &f[1], &
					ia[1], &ja[1], &ig[1], eps, imin, imax, iminw, imaxw, &ig[icgst],
					&ig[icgst + ndicg], nstcol, iarr, time, &ncyc0, &irow0, &levels,
					nda, ndja, ndu, ndf, &mda, &mdja, &mdu, &mdf, &iup, &ium, resi, &
					res0, &res);
					system("pause");
					}
					*/
	printf("sizeof  ndu=%lld nnu=%lld ndf=%lld\n", ndu[0], nnu[0], ndf[0]);
	integer nnz;
	/*
	integer n75 = -1 , nnz;

	for (integer i_83 = 1; i_83 <= nda[0]; i_83++) {
	if (ja[i_83] > n75) n75 = ja[i_83];
	}
	nnz = ia[n75 + 1];
	printf("n_75=%d nnz=%d\n",n75,nnz);
	system("pause");
	*/



	// Разреженная матрица СЛАУ
	// в CRS формате.
	simplesparsetoCRS(sparseM, val75, col_ind75, row_ptr75, n); // преобразование матрицы из одного формата хранения в другой.
	

	if ((val75 == NULL) || (col_ind75 == NULL) || (row_ptr75 == NULL)) {
		// недостаточно памяти на данном оборудовании.
		printf("Problem: not enough memory on your equipment for val, col_ind or row_ptr: bicgStab + camg...\n");
		printf("Please any key to exit...\n");
		exit(1);
	}

	nnz = row_ptr75[n75];
	printf("n=%lld nnz=%lld ndu=%lld nda=%lld\n", n75, nnz, ndu[0], nda[0]);



	/*
	// Так делать нельзя, т.к. матрица СЛАУ (a,ia,ja)
	// была модифицирована в процессе работы amg1r5.f алгоритма.

	// инициализация матрицы.
	#pragma omp parallel for
	for (integer i_91 = 1; i_91 <= n75; i_91++) {

	for (integer i_92 = ia[i_91]; i_92 < ia[i_91 + 1]; i_92++) {

	val75[i_92 - 1] = a[i_92];
	col_ind75[i_92 - 1] = ja[i_92]- 1;

	}
	row_ptr75[i_91 - 1] = ia[i_91] - 1;
	}
	row_ptr75[n75] = ia[n75 + 1] - 1;
	*/

	//system("PAUSE");





	//y76 = new doublereal[n75 + 1];
	//pi76 = new doublereal[n75 + 1];
	y76 = new doublereal[ndu[0] + 1];
	pi76 = new doublereal[ndu[0] + 1];

	//z76 = new doublereal[n75 + 1];
	//s76 = new doublereal[n75 + 1];
	z76 = new doublereal[ndu[0] + 1];
	s76 = new doublereal[ndu[0] + 1];

	ri75 = new doublereal[n75];
	roc75 = new doublereal[n75];
	s75 = new doublereal[n75];
	t75 = new doublereal[n75];
	vec75 = new doublereal[n75];
	vi75 = new doublereal[n75];
	pi75 = new doublereal[n75];
	dx75 = new doublereal[n75];
	dax75 = new doublereal[n75];
	y75 = new doublereal[n75];
	z75 = new doublereal[n75];
	if ((ri75 == NULL) || (roc75 == NULL) || (s75 == NULL) || (t75 == NULL) || (vi75 == NULL) || (pi75 == NULL) || (dx75 == NULL) || (dax75 == NULL) || (y75 == NULL) || (z75 == NULL)) {
		// недостаточно памяти на данном оборудовании.
		printf("Problem: not enough memory on your equipment for: bicgStab + camg...\n");
		printf("Please any key to exit...\n");
		exit(1);
	}



	
	if (iVar == TEMP) {
		epsilon75 *= 1.0e-4; // 1.0e-4
	}
	if (iVar == TOTALDEFORMATIONVAR) {
		epsilon75 *= 1.0e-4; // 1.0e-4
		//epsilon75 *= 1.0e-12;
	}


	// initialize
#pragma omp parallel for
	for (i75 = 0; i75<n75; i75++) {
		s75[i75] = 0.0;
		t75[i75] = 0.0;
		vi75[i75] = 0.0;
		pi75[i75] = 0.0;
		// инициализатор массивов для предобуславливания
		y75[i75] = 0.0;
		z75[i75] = 0.0;
		// результат умножения матрицы на вектор.
		dax75[i75] = 0.0;
		// Начальное приближение.
		dx75[i75] = u[i75 + 1];
	}


	// Умножение матрицы на вектор. Нумерации векторов начинаются с нуля.
	MatrixCRSByVector(val75, col_ind75, row_ptr75, dx75, dax75, n75); // результат занесён в  dax75

																	  // Вычисление ri75 и roc75.
#pragma omp parallel for
	for (i75 = 0; i75 < n75; i75++) {
		ri75[i75] = f[i75 + 1] - dax75[i75];
		roc75[i75] = 1.0;
	}
	delta075 = NormaV(ri75, n75);


	// Если решение сразу хорошее то не считать:
	if (iVar == TEMP) {
		if (fabs(delta075)<1.0e-4*dterminatedTResudual) iflag75 = 0;
	}
	else {
		if (fabs(delta075)<dterminatedTResudual) iflag75 = 0;
	}

	if (fabs(delta075)<1e-14) iflag175 = 0;
	if ((iVar == TEMP) && (iflag75 == 0) && (iflag175 == 0)) {
#if doubleintprecision == 1
		printf("bicgStab+camg: iflag=%lld, iflag1=%lld, delta0=%e\n", iflag75, iflag175, delta075);
#else
		printf("bicgStab+camg: iflag=%d, iflag1=%d, delta0=%e\n", iflag75, iflag175, delta075);
#endif

		system("PAUSE");
	}


	if (n75<30000) {
		// задача очень малой размерности !
		if ((iVar == VX) || (iVar == VY) || (iVar == VZ)) {
			iN75 = 1; // обязательно нужна хотя бы одна итерация.
					  // если этого будет недостаточно то мы всё равно будем итерировать до тех пор пока невязка не станет меньше epsilon.
			if (1.0e-3*fabs(delta075)<epsilon75) {
				epsilon75 = 1.0e-3*fabs(delta075);
			}
			if (iflag175 == 1) {
				iflag75 = 1;
			}
		}
		if (iVar == TEMP) {
			iN75 = 2;
			epsilon75 = fmin(0.1*fabs(delta075), epsilon75);
			if (bSIMPLErun_now_for_temperature  ) {
				//printf("epsilon=%e \n",epsilon);
				//system("pause");
				// Экспериментальным образом обнаружена недоэтерированость по температуре для гидродинамического решателя.
				// поэтому точность было решено увеличить на 5 порядков.
				// 27.07.2016
				epsilon75 *= 1e-10;
				iN75 = 20;
				//epsilon75 *= 1e-16;
				//iN75 = 30;
			}
		}
		if (iVar == PAM) {
			iN75 = 3; // решение для поправки давления должно быть получено точно.
			if (1.0e-3*fabs(delta075)<epsilon75) {
				epsilon75 = 1.0e-3*fabs(delta075);
			}
			if (iflag175 == 1) {
				iflag75 = 1;
			}
			//printf("%e",epsilon75); system("pause");
		}
	}
	else if ((n75 >= 30000) && (n75 < 100000)) {
		// Здесь я немного увеличил число итераций и 
		// скорректировал условие окончания чтобы считало 
		// поточнее, но это не повлияло.
		// Главный вопрос в том что невязка по температуре почему-то не меняется.
		// задача небольшой размерности.
		if ((iVar == VX) || (iVar == VY) || (iVar == VZ)) {
			iN75 = 3; // обязательно нужна хотя бы одна итерация.
					  // если этого будет недостаточно то мы всё равно будем итерировать до тех пор пока невязка не станет меньше epsilon.
			if (1.0e-3*fabs(delta075)<epsilon75) {
				epsilon75 = 1.0e-3*fabs(delta075);
			}
			if (iflag175 == 1) {
				iflag75 = 1;
			}
			// 27.07.2016
			iN75 = 12;
			epsilon75 *= 1e-2;
		}
		if (iVar == TEMP) {
			iN75 = 4;
			epsilon75 = fmin(0.1*fabs(delta075), epsilon75);
			if (bSIMPLErun_now_for_temperature  ) {
				//printf("epsilon75=%e \n",epsilon75);
				//system("pause");
				// Экспериментальным образом обнаружена недоэтерированость по температуре для гидродинамического решателя.
				// поэтому точность было решено увеличить на 5 порядков.
				// 27.07.2016
				epsilon75 *= 1e-10;
				iN75 = 20;
				//epsilon75 *= 1e-16;
				//iN75 = 30;
			}
		}
		if (iVar == PAM) {
			iN75 = 6; // решение для поправки давления должно быть получено точно.
			if (1.0e-3*fabs(delta075)<epsilon75) {
				epsilon75 = 1.0e-3*fabs(delta075);
			}
			if (iflag175 == 1) {
				iflag75 = 1;
			}
			//printf("%e",epsilon75); system("pause");
			// 27.07.2016.
			epsilon75 *= 1e-2;
			iN75 = 20;
		}
	}
	else if ((n75 >= 100000) && (n75<300000)) {
		// задача небольшой средней размерности.
		if ((iVar == VX) || (iVar == VY) || (iVar == VZ)) {
			iN75 = 3; // обязательно нужна хотя бы одна итерация.
					  // Вообще говоря невязка для скоростей падает очень быстро поэтому всегда достаточно iN итераций для скорости.
					  // если этого будет недостаточно то мы всё равно будем итерировать до тех пор пока невязка не станет меньше epsilon.
			if (1.0e-3*fabs(delta075)<epsilon75) {
				epsilon75 = 1.0e-3*fabs(delta075);
			}
			if (iflag175 == 1) {
				iflag75 = 1;
			}
		}
		if (iVar == TEMP) {
			iN75 = 4;
			epsilon75 = fmin(0.1*fabs(delta075), epsilon75);
			if (bSIMPLErun_now_for_temperature  ) {
				//printf("epsilon75=%e \n",epsilon75);
				//system("pause");
				// Экспериментальным образом обнаружена недоэтерированость по температуре для гидродинамического решателя.
				// поэтому точность было решено увеличить на 5 порядков.
				// 27.07.2016
				epsilon75 *= 1e-10;
				iN75 = 20;
				//epsilon75 *= 1e-16;
				//iN75 = 30;
			}
		}
		if (iVar == PAM) {
			iN75 = 8; // решение для поправки давления должно быть получено точно.
			if (1.0e-4*fabs(delta075)<epsilon75) {
				epsilon75 = 1.0e-4*fabs(delta075);
			}
			if (iflag175 == 1) {
				iflag75 = 1;
			}
			//printf("%e",epsilon75); system("pause");
		}
	}
	else if ((n75 >= 300000) && (n75<1000000)) {
		// задача истинно средней размерности.
		if ((iVar == VX) || (iVar == VY) || (iVar == VZ)) {
			iN75 = 3; // обязательно нужна хотя бы одна итерация.
					  // если этого будет недостаточно то мы всё равно будем итерировать до тех пор пока невязка не станет меньше epsilon.
			if (1.0e-3*fabs(delta075)<epsilon75) {
				epsilon75 = 1.0e-3*fabs(delta075);
			}
			if (iflag175 == 1) {
				iflag75 = 1;
			}
		}
		if (iVar == TEMP) {
			iN75 = 4;
			epsilon75 = 1e-5*fmin(0.1*fabs(delta075), epsilon75);
			if (bSIMPLErun_now_for_temperature  ) {
				//printf("epsilon75=%e \n",epsilon75);
				//system("pause");
				// Экспериментальным образом обнаружена недоэтерированость по температуре для гидродинамического решателя.
				// поэтому точность было решено увеличить на 5 порядков.
				// 27.07.2016
				epsilon75 *= 1e-8;
				iN75 = 20;
				//epsilon75 *= 1e-16;
				//iN75 = 30;
			}
		}
		if (iVar == PAM) {
			iN75 = 16; // решение для поправки давления должно быть получено точно.
			if (1.0e-4*fabs(delta075)<epsilon75) {
				epsilon75 = 1.0e-4*fabs(delta075);
			}
			if (iflag175 == 1) {
				iflag75 = 1;
			}
			//printf("%e",epsilon75); system("pause");
		}
	}
	else if ((n75 >= 1000000) && (n75<3000000)) {
		// задача достаточно большой размерности.
		if ((iVar == VX) || (iVar == VY) || (iVar == VZ)) {
			iN75 = 6; // обязательно нужна хотя бы одна итерация.
					  // если этого будет недостаточно то мы всё равно будем итерировать до тех пор пока невязка не станет меньше epsilon.
			if (1.0e-3*fabs(delta075)<epsilon75) {
				epsilon75 = 1.0e-3*fabs(delta075);
			}
			if (iflag175 == 1) {
				iflag75 = 1;
			}
		}
		if (iVar == TEMP) {
			iN75 = 8;
			epsilon75 = 1e-5*fmin(0.1*fabs(delta075), epsilon75);
			if (bSIMPLErun_now_for_temperature  ) {
				//printf("epsilon75=%e \n",epsilon75);
				//system("pause");
				// Экспериментальным образом обнаружена недоэтерированость по температуре для гидродинамического решателя.
				// поэтому точность было решено увеличить на 5 порядков.
				// 27.07.2016
				epsilon75 *= 1e-8;
				iN75 = 20;
				//epsilon75 *= 1e-16;
				//iN75 = 30;
			}
		}
		if (iVar == PAM) {
			iN75 = 23; // решение для поправки давления должно быть получено точно.
			if (1.0e-4*fabs(delta075)<epsilon75) {
				epsilon75 = 1.0e-4*fabs(delta075);
			}
			if (iflag175 == 1) {
				iflag75 = 1;
			}
			//printf("%e",epsilon75); system("pause");
		}
	}
	else if (n75 >= 3000000) {
		// задача очень большой размерности.
		if ((iVar == VX) || (iVar == VY) || (iVar == VZ)) {
			iN75 = 6; // обязательно нужна хотя бы одна итерация.
					  // если этого будет недостаточно то мы всё равно будем итерировать до тех пор пока невязка не станет меньше epsilon.
			if (1.0e-3*fabs(delta075)<epsilon75) {
				epsilon75 = 1.0e-3*fabs(delta075);
			}
			if (iflag175 == 1) {
				iflag75 = 1;
			}
		}
		if (iVar == TEMP) {
			iN75 = 8;
			epsilon75 = 1e-10*fmin(0.1*fabs(delta075), epsilon75);
		}
		if (iVar == PAM) {
			iN75 = 36; // решение для поправки давления должно быть получено точно.
			if (1.0e-4*fabs(delta075)<epsilon75) {
				epsilon75 = 1.0e-4*fabs(delta075);
			}
			if (iflag175 == 1) {
				iflag75 = 1;
			}
			//printf("%e",epsilon); system("pause");
		}
	}


	if (iVar == TEMP) {
		maxit75 = 2000;
	}
	if (iVar == PAM) {
		maxit75 = 2000; // 2000
	}
	if ((iVar == VX) || (iVar == VY) || (iVar == VZ)) {
		maxit75 = 100;//100
	}
	if (iVar == TOTALDEFORMATIONVAR) {
		maxit75 = 2000; // 2000
		if (1.0e-4*fabs(delta075) < epsilon75) {
			epsilon75 = 1.0e-4*fabs(delta075);
		}
		epsilon75 = 1.0e-16;
		iN75 = 100;
		if (iflag175 == 1) {
			iflag75 = 1;
		}

	}







	// Мы обязательно должны сделать несколько итераций. (не менее 10).
	// Если только решение не удовлетворяет уравнению тождественно.
	while (((icount75 < iN75) && (iflag175 != 0)) || (iflag75 != 0 && icount75 < maxit75)) {

		// 6.01.2017: Body BiCGStab + AMG. (BiCGStab_internal4).


		icount75++;

		count_iter_for_film_coef75++;
		// В случае задачи Ньютона - Рихмана, Стефана-Больцмана и миксового условия не итерируем до конца обрываем, 
		// т.к. нам требуется частая пересборка матрицы. 13 марта 2016.
		//if (((adiabatic_vs_heat_transfer_coeff > ADIABATIC_WALL_BC) || (breakRUMBAcalc_for_nonlinear_boundary_condition)) && (count_iter_for_film_coef75>5)) break;

		roi75 = Scal(roc75, ri75, n75);
		bet75 = (roi75 / roim175)*(al75 / wi75);


		//printf("%e %e %e %e\n",roi75,roim175,al75,wi75);
		//system("pause");

#pragma omp parallel for 
		for (i75 = 0; i75<n75; i75++) {
			doublereal pibuf75 = ri75[i75] + (pi75[i75] - vi75[i75] * wi75)*bet75;
			pi75[i75] = pibuf75;
		}

		// Первое предобуславливание.
		// Ky=pi
#pragma omp parallel for
		for (i75 = 0; i75 < n75; i75++) {
			y75[i75] = 0.0; // Если начинать не с нуля то небудет сходимости для PAM !.
			y76[i75 + 1] = 0.0;//+1
			pi76[i75 + 1] = pi75[i75];//+1
		}


		// multigrid Ruge and Stuben preconditioning [1986].
		// достаточно одного V цикла.
		// K*y76 = pi76;
		ifirst[0] = 11; // Нулевое начальное приближение
		for (integer i_numberV_cycle = 0; i_numberV_cycle < 1; i_numberV_cycle++) {
			ncyc[0] = 1011; // достаточно одного V цикла.
			solve_(madapt, ncyc, nrd, nsolco, nru, &kout, ierr, &a[1], &y76[1], &pi76[1], &
				ia[1], &ja[1], &ig[1], eps, imin, imax, iminw, imaxw, &ig[icgst],
				&ig[icgst + ndicg], nstcol, iarr, time, &ncyc0, &irow0, &levels,
				nda, ndja, ndu, ndf, &mda, &mdja, &mdu, &mdf, &iup, &ium, resi, &
				res0, &res);
			//system("pause");
		}


		// Возвращение результата.
#pragma omp parallel for
		for (i75 = 0; i75 < n75; i75++) {
			y75[i75] = y76[i75 + 1];//+1
		}

		MatrixCRSByVector(val75, col_ind75, row_ptr75, y75, vi75, n75); // vi==A*y;

		if ((fabs(roi75)<1e-30) && (fabs(Scal(roc75, vi75, n75))<1e-30)) {
			al75 = 1.0;
		}
		else if (fabs(roi75)<1e-30) {
			al75 = 0.0;
		}
		else {
			al75 = roi75 / Scal(roc75, vi75, n75);
		}


#pragma omp parallel for
		for (i75 = 0; i75<n75; i75++) {
			s75[i75] = ri75[i75] - al75 * vi75[i75];
		}

		// Второе предобуславливание.
		// Kz=s

#pragma omp parallel for
		for (i75 = 0; i75<n75; i75++) z75[i75] = 0.0; // Если начинать не с нуля то небудет сходимости для PAM !.

#pragma omp parallel for
		for (i75 = 0; i75 < n75; i75++) {
			vec75[i75] = s75[i75];
			z76[i75 + 1] = 0.0;//+1
			s76[i75 + 1] = s75[i75];//+1
		}

		// multigrid Ruge and Stuben preconditioning [1986].
		// достаточно одного V цикла.
		// K*z76 = s76;
		ifirst[0] = 11; // Нулевое начальное приближение
		for (integer i_numberV_cycle = 0; i_numberV_cycle < 1; i_numberV_cycle++) {
			ncyc[0] = 1011; // достаточно одного V цикла.
			solve_(madapt, ncyc, nrd, nsolco, nru, &kout, ierr, &a[1], &z76[1], &s76[1], &
				ia[1], &ja[1], &ig[1], eps, imin, imax, iminw, imaxw, &ig[icgst],
				&ig[icgst + ndicg], nstcol, iarr, time, &ncyc0, &irow0, &levels,
				nda, ndja, ndu, ndf, &mda, &mdja, &mdu, &mdf, &iup, &ium, resi, &
				res0, &res);
			//system("pause");
		}

#pragma omp parallel for
		for (i75 = 0; i75 < n75; i75++) {
			s75[i75] = vec75[i75];
			// Возвращаем результат.
			z75[i75] = z76[i75 + 1];//+1
		}

		MatrixCRSByVector(val75, col_ind75, row_ptr75, z75, t75, n75); // t==A*z;

		wi75 = Scal(t75, s75, n75) / Scal(t75, t75, n75);
		// printf("%e %e",Scal(t75,s75,n75),Scal(t75,t75,n75));

#pragma omp parallel for
		for (i75 = 0; i75<n75; i75++) {
			//dx75[i75]+=al75*pi75[i75]+wi75*s75[i75]; // так было без предобуславливателя
			dx75[i75] += al75 * y75[i75] + wi75 * z75[i75]; // так стало с предобуславливателем
			ri75[i75] = s75[i75] - wi75 * t75[i75];
		}
		deltai75 = NormaV(ri75, n75);

		//printf("deltai75=%e\n",deltai75); system("pause");

		// печать невязки на консоль
		bool bprint_mesage_diagnostic = true;
		if (bprint_mesage_diagnostic) {
			if ((icount75 % 10) == 0) {
				std::cout << "iter  residual" << std::endl;				
			}
			std::cout << icount75 <<" "<< deltai75 << std::endl;
		}

		// 28.07.2016.
		//std::cout << icount75 <<" "<< deltai75 << std::endl;		

		//system("pause");
		if (deltai75 > delta_old_iter75) i_signal_break_pam_opening75++;
		delta_old_iter75 = deltai75;
		if (iVar == PAM) {
			if (i_signal_break_pam_opening75 > i_limit_signal_pam_break_opening75) {
				// досрочный выход из цикла.
				std::cout << "icount PAM=" << icount75 << std::endl;
				break;
			}
		}

		if (deltai75 <epsilon75) iflag75 = 0; // конец вычисления
		else {
			/*
			// 02.05.2018
			if (bSIMPLErun_now_for_temperature) {
			// Эти значения невязок для CFD задач были
			// успешно опробованы на задаче теплового расчёта
			// радиатора водяного охлаждения 3л/мин (совместное решение cfd + temperature
			// + приближение Обербека-Буссинеска.).
			switch (iVar) {
			case VX: if (deltai75/ delta075 < 1e-2) iflag75 = 0;   break; //5e-5
			case VY: if (deltai75/ delta075 < 1e-2) iflag75 = 0;   break; // 5e-5
			case VZ: if (deltai75/ delta075 < 1e-2) iflag75 = 0;   break; // 5e-5
			//case TEMP: tolerance = 1e-8;  break; // 1e-7
			case PAM: if (deltai75/ delta075 < 1e-1) iflag75 = 0;  break; // 1e-6
			}
			}
			else {
			if (iflag75 != 0) {
			*/
			roim175 = roi75;
			//}
			//}
		}

		if (iVar == TEMP) {
#if doubleintprecision == 1
			//printf("epsilon=%e deltai=%e icount=%lld\n",epsilon75,deltai75, icount75);
#else
			//printf("epsilon=%e deltai=%e icount=%d\n",epsilon75,deltai75, icount75);
#endif

			//system("pause");
		}

		//icount_V_cycle = icount75; // количество итераций в BiCGStabP для лога.

		//if (icount75 > 2600) break; // 15.02.2017

	}

	// Возвращение результата вычислений.
#pragma omp parallel for
	for (i75 = 0; i75 < n75; i75++) {
		u[i75 + 1] = dx75[i75];
		//x_best_search[i75 + 1] = dx75[i75];
	}

	// Освобождение оперативной памяти.
	// Первое предобуславливание
	if (pi76 != NULL) {
		delete[] pi76;
		pi76 = NULL;
	}
	if (y76 != NULL) {
		delete[] y76;
		y76 = NULL;
	}
	// Второе предобуславливание
	if (z76 != NULL) {
		delete[] z76;
		z76 = NULL;
	}
	if (s76 != NULL) {
		delete[] s76;
		s76 = NULL;
	}
	if (ri75 != NULL) {
		delete[] ri75;
		ri75 = NULL;
	}
	if (roc75 != NULL) {
		delete[] roc75;
		roc75 = NULL;
	}
	if (s75 != NULL) {
		delete[] s75;
		s75 = NULL;
	}
	if (t75 != NULL) {
		delete[] t75;
		t75 = NULL;
	}
	if (vec75 != NULL) {
		delete[] vec75;
		vec75 = NULL;
	}
	if (vi75 != NULL) {
		delete[] vi75;
		vi75 = NULL;
	}
	if (pi75 != NULL) {
		delete[] pi75;
		pi75 = NULL;
	}
	if (dx75 != NULL) {
		delete[] dx75;
		dx75 = NULL;
	}
	if (dax75 != NULL) {
		delete[] dax75;
		dax75 = NULL;
	}
	if (y75 != NULL) {
		delete[] y75;
		y75 = NULL;
	}
	if (z75 != NULL) {
		delete[] z75;
		z75 = NULL;
	}

	// Освобождение оперативной памяти.
	if (val75 != NULL) {
		delete[] val75;
		val75 = NULL;
	}
	if (col_ind75 != NULL) {
		delete[] col_ind75;
		col_ind75 = NULL;
	}
	if (row_ptr75 != NULL) {
		delete[] row_ptr75;
		row_ptr75 = NULL;
	}



	if (debug_reshime) system("pause");
	//system("pause");

	if (*ierr > 0) {
		return 0;
	}
L10:
	wrkcnt_(&kout, &ia[1], &ig[1], imin, imax, iminw, &levels, time, &ncyc0, &
		iup, &mda, &mdia, &mdja, &mdu, &mdf, &mdig, &res0, &res);
	return 0;
L60:
	//io___29.ciunit = ium;
	//s_wsfe(&io___29);
	//e_wsfe();
	printf("*** ERROR IN AMG1R5: NDIG TOO SMALL ***\n");

	*ierr = 6;
	return 0;
L70:
	//io___30.ciunit = ium;
	//s_wsfe(&io___30);
	//e_wsfe();

	printf("*** ERROR IN AMG1R5: ILLEGAL PARAMETER LEVELX ***\n");

	*ierr = 18;
	return 0;
} /* amg1r5_Vorst_modification_matrix_Assemble2 */




/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */

/*     GENERAL AMG1R5 SETUP-SUBROUTINES */

/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
/* ....................................................................... */

/*     SETUP                                                  SUBROUTINE */

/* ....................................................................... */

/* Subroutine */ integer setup_(integer *nnu, integer *matrix, integer *levelx, 
	doublereal *ecg1, doublereal *ecg2, doublereal *ewt2, integer *nwt, 
	integer *ntr, integer *ierr, doublereal *a, doublereal *u, integer *
	ia, integer *ja, integer *iw, integer *imin, integer *imax, integer *
	iminw, integer *imaxw, integer *icg, integer *ifg, integer *nstcol, 
	integer *iarr, /*real*/ unsigned int *time, integer *levels, integer *irow0, integer *
	nda, integer *ndia, integer *ndja, integer *ndu, integer *ndf, 
	integer *ndicg, integer *ium, integer *mda, integer *mdia, integer *
	mdja, integer *mdu, integer *mdf, integer *mdig, integer *maxgr, 
	integer *lratio)
{
    integer i__=0;
    extern /* Subroutine */ integer idec_(integer *, integer *, integer *, 
	    integer *);
    integer isym=0;
    extern /* Subroutine */ integer check_(integer *, doublereal *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, /*real*/ unsigned int *, 
	    integer *, integer *, integer *, integer *, integer *, integer *),
	     crsng_(integer *, doublereal *, doublereal *, doublereal *, 
	    integer *, integer *, integer *, doublereal *, doublereal *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    /*real*/ unsigned int *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *);
    integer ndigit=0;


/*     PREPARATION PHASE OF AMG1R5 (GENERAL PART) */


/* ===> DECOMPOSE "MATRIX" */

    /* Parameter adjustments */
    --time;
    --iarr;
    --nstcol;
    --ifg;
    --icg;
    --imaxw;
    --iminw;
    --imax;
    --imin;
    --iw;
    --ja;
    --ia;
    --u;
    --a;

    /* Function Body */
    idec_(matrix, &c__2, &ndigit, &iarr[1]);
    isym = iarr[1];
    *irow0 = iarr[2];

/* ===> PREPARATION (IGNORED IN TIMING) */

    imin[1] = 1;
    imax[1] = *nnu;
    check_(ierr, &a[1], &ia[1], &ja[1], &imin[1], &imax[1], &icg[1], &ifg[1], 
	    &time[1], &isym, irow0, nda, ndicg, ndja, ium);
    if (*ierr > 0) {
	return 0;
    }

/* ===> RESET TIME COUNTERS */

    for (i__ = 1; i__ <= 20; ++i__) {
	//time[i__] = 0.f;
		time[i__] = 0;
/* L30: */
    }

/* ===> DEFINE COARSER GRIDS + OPERATORS. RESET LEVELS IF NECESSARY. */

    crsng_(levelx, ecg1, ecg2, ewt2, nwt, ntr, ierr, &a[1], &u[1], &ia[1], &
	    ja[1], &iw[1], &imin[1], &imax[1], &iminw[1], &imaxw[1], &icg[1], 
	    &ifg[1], &nstcol[1], levels, irow0, nda, ndja, ndia, ndu, ndf, 
	    ndicg, &time[1], ium, mda, mdia, mdja, mdu, mdf, mdig, maxgr, 
	    lratio);
    return 0;
} /* setup_ */


/* ....................................................................... */

/*     CHECK                                                 SUBROUTINE */

/* ....................................................................... */

/* Subroutine */ integer check_(integer *ierr, doublereal *a, integer *ia, 
	integer *ja, integer *imin, integer *imax, integer *icg, integer *ifg,
	 /*real*/ unsigned int *time, integer *isym, integer *irow0, integer *nda, integer *
	ndicg, integer *ndja, integer *ium)
{
    /* Format strings */
    //static char fmt_9000[] = "(\002 CHECK: A PROBABLY SYMMETRIC\002)";
    //static char fmt_9005[] = "(\002 CHECK: A PROBABLY NOT SYMMETRIC. MEASU"
	//    "RE:\002,d11.3)";
    //static char fmt_9010[] = "(\002 CHECK: A PROBABLY NOT POS. TYPE:\002,i6"
	  //  ",\002 OFF-DIAGONAL\002,\002 ELEMENTS POSITIVE\002)";
    //static char fmt_9020[] = "(\002 CHECK: A PROBABLY NOT POS. TYPE:\002,i6"
	  //  ",\002 ROWSUMS\002,\002 NEGATIVE\002)";
    //static char fmt_9025[] = "(\002 CHECK: A PROBABLY SINGULAR - ROWSUMS ARE"
	//    " ZERO\002)";
   // static char fmt_9030[] = "(\002 CHECK: A PROBABLY POSITIVE TYPE\002)";
    //static char fmt_9100[] = "(\002 --- WARNG IN CHECK: PARAM MATRIX MAY BE "
	//    "BAD ---\002)";
    //static char fmt_9600[] = "(\002 CHECK: MATRIX A WAS SYMMETRICALLY STORE"
	//    "D\002)";
    //static char fmt_9500[] = "(\002 CHECK: STORAGE OF A HAS BEEN SYMMETRIZED"
	//    " BY\002,\002 INTRODUCING\002,i6,\002 ZEROES\002)";
   // static char fmt_9220[] = "(\002 *** WARNG IN CHECK:\002,i5,\002 A-ENTRIE"
	//    "S MISSING ***\002/\002     MISSING TRANSPOSE CONNECTIONS WILL BE"
	 //   " FILLED IN\002)";
    //static char fmt_9200[] = "(\002 *** ERROR IN CHECK: NDIG TOO SMALL **"
	//    "*\002)";
    //static char fmt_9210[] = "(\002 *** ERROR IN CHECK: NDA TOO SMALL ***"
	  //  "\002)";
    //static char fmt_9230[] = "(\002 *** ERROR IN CHECK: NDJA TOO SMALL **"
	//    "*\002)";
    //static char fmt_9300[] = "(\002 *** ERROR IN CHECK: POINTER IA ERRONEOUS"
	//    " ***\002)";
    //static char fmt_9310[] = "(\002 *** ERROR IN CHECK: POINTER JA ERRONEOUS"
	//    " ***\002)";
    //static char fmt_9320[] = "(\002 *** ERROR IN CHECK: DIAGONAL IS NOT STOR"
	//    "ED FIRST ***\002)";
    //static char fmt_9330[] = "(\002 *** ERROR IN CHECK: DIAGONAL IS NON-POSI"
	//    "TIVE ***\002)";

    /* System generated locals */
    integer i__1=0, i__2=0, i__3=0;
    doublereal d__1=0.0;

    /* Builtin functions */
    //double sqrt(double);
    //integer s_wsfe(cilist *), e_wsfe(void), do_fio(integer *, char *, ftnlen);

    /* Local variables */
    doublereal d__=0.0;
    integer i__=0, j=0, i1=0, j1=0, j2=0,  new__=0;
	integer nna=0, nnu=0;
    doublereal deps=0.0;
    integer jnew=0;
    doublereal asym=0.0;
    integer naneg=0, naoff=0, nazer=0, napos=0;
    extern /* Subroutine */ integer trunc_(integer *, integer *, doublereal *, 
	    integer *, integer *, integer *, integer *, /*real*/ unsigned int *, integer *, 
	    integer *);
    integer ishift=0;
    doublereal anormm=0.0, anormp=0.0, rowsum=0.0;

    /* Fortran I/O blocks */
    //static cilist io___52 = { 0, 0, 0, fmt_9000, 0 };
    //static cilist io___53 = { 0, 0, 0, fmt_9005, 0 };
    //static cilist io___54 = { 0, 0, 0, fmt_9010, 0 };
    //static cilist io___55 = { 0, 0, 0, fmt_9020, 0 };
    //static cilist io___56 = { 0, 0, 0, fmt_9025, 0 };
    //static cilist io___57 = { 0, 0, 0, fmt_9030, 0 };
    //static cilist io___58 = { 0, 0, 0, fmt_9100, 0 };
    //static cilist io___59 = { 0, 0, 0, fmt_9600, 0 };
    //static cilist io___62 = { 0, 0, 0, fmt_9500, 0 };
    //static cilist io___63 = { 0, 0, 0, fmt_9220, 0 };
    //static cilist io___64 = { 0, 0, 0, fmt_9200, 0 };
    //static cilist io___65 = { 0, 0, 0, fmt_9210, 0 };
    //static cilist io___66 = { 0, 0, 0, fmt_9230, 0 };
    //static cilist io___67 = { 0, 0, 0, fmt_9300, 0 };
    //static cilist io___68 = { 0, 0, 0, fmt_9310, 0 };
    //static cilist io___69 = { 0, 0, 0, fmt_9320, 0 };
    //static cilist io___70 = { 0, 0, 0, fmt_9330, 0 };



/*     CHECKS FOR SEVERAL PROPERTIES OF A, IA, JA. IN PARTICULAR, */
/*     CHECKS FOR SYMMETRIC STORAGE OF GIVEN MATRIX (I.E. L(I,J) IS */
/*     STORED IN A IFF L(J,I) IS STORED IN A). */

/*     IF STORAGE IS SYMMETRIC, PAIRS OF ZEROES ARE REMOVED (IF THERE */
/*     ARE ANY) AND PROGRAM EXECUTION CONTINUES. */

/*     IF, HOWEVER, STORAGE OF A IS NOT SYMMETRIC, IT IS SYMMETRIZED. */
/*     IF ISYM.EQ.1, THE MISSING TRANSPOSE CONNECTIONS ARE COPIED. */
/*     IF ISYM.NE.1, THE STORAGE OF THE MATRIX A IS SYMMETRIZED (BY */
/*     ADDING CERTAIN ZERO ELEMENTS). THEN PAIRS OF ZEROES ARE REMOVED */
/*     (IF THERE ARE ANY) AND PROGRAM EXECUTION CONTINUES. */

/*     ARRAYS USED FOR TEMPORARY STORAGE: ICG, IFG */


/* ===> CHECK IF POINTERS IA, JA ARE REASONABLE */

    /* Parameter adjustments */
    --time;
    --ifg;
    --icg;
    --imax;
    --imin;
    --ja;
    --ia;
    --a;

    /* Function Body */
    if (ia[1] != 1) {
	goto L3300;
    }
    nnu = imax[1];
	printf("nnu=%lld\n",nnu);
    if (nnu >= *ndicg) {
	goto L3200;
    }
    i__1 = nnu;
    for (i__ = 1; i__ <= i__1; ++i__) {
	icg[i__] = 0;
/* L2: */
    }
    i__1 = nnu;
    for (i__ = 1; i__ <= i__1; ++i__) {
	j1 = ia[i__];
	j2 = ia[i__ + 1] - 1;
	if (j2 < j1 || j2 > *nda) {
	    goto L3300;
	}
	// Первый элемент всегда диагональный.
	if (ja[j1] != i__) {
	    goto L3320;
	}
	i__2 = j2;
	for (j = j1; j <= i__2; ++j) {
		// Сканируем строку.
	    i1 = ja[j];//номер столбца.
	    if (i1 < 1 || i1 > nnu || icg[i1] == 1) {
			//icg[i1] == 1 - столбец был посещен.
			printf("nnu=%lld i1=%lld icg=%lld\n",nnu,i1, icg[i1]);
			if (icg[i1] == 1) {
				printf("two identical columns in a row\n");
				for (integer j69 = j1; j69 <= i__2; ++j69) {
					printf("%lld %lld %lld val=%e col_ind=%lld row_ind=%lld\n",j69,j1,i__2,a[j69],ja[j69],i__);
				}
			}
			system("PAUSE");
		goto L3310;
	    }
	    icg[i1] = 1;// столбец был посещен.
/* L5: */
	}
	i__2 = j2;
	for (j = j1; j <= i__2; ++j) {
	    icg[ja[j]] = 0;// сброс посещения.
/* L7: */
	}
/* L10: */
    }
    nna = ia[nnu + 1] - 1;

/* ===> CHECK FOR PROPERTIES OF A. IN PARTICULAR, COUNT */
/* ===> MISSING STORAGE PLACES ("NEW"). RETURN IF NEW=0 */

    anormm = 0.;
    anormp = 0.;
    naoff = 0;
    napos = 0;
    naneg = 0;
    nazer = 0;
    new__ = 0;
    i__1 = nnu;
    for (i__ = 1; i__ <= i__1; ++i__) {
	icg[i__] = ia[i__ + 1] - ia[i__];
/* L100: */
    }

    i__1 = nnu;
    for (i__ = 1; i__ <= i__1; ++i__) {
	d__ = a[ia[i__]];
	if (d__ <= 0.) {
	    goto L3330;
	}
	deps = d__ * 1e-12;
	rowsum = d__;
/* Computing 2nd power */
	d__1 = d__;
	anormp += d__1 * d__1 * 2.;
	i__2 = ia[i__ + 1] - 1;
	for (j = ia[i__] + 1; j <= i__2; ++j) {
	    rowsum += a[j];
	    if (a[j] >= deps) {
		++naoff;
	    }
	    i1 = ja[j];
	    if (i1 > 0) {
		goto L140;
	    }
	    ja[j] = -i1;
	    goto L150;
L140:
	    if (i1 < i__) {
		goto L135;
	    }
	    i__3 = ia[i1 + 1] - 1;
	    for (j1 = ia[i1] + 1; j1 <= i__3; ++j1) {
		if (ja[j1] != i__) {
		    goto L130;
		}
		ja[j1] = -ja[j1];
/* Computing 2nd power */
		d__1 = a[j] - a[j1];
		anormm += d__1 * d__1;
/* Computing 2nd power */
		d__1 = a[j] + a[j1];
		anormp += d__1 * d__1;
		goto L150;
L130:
		;
	    }
L135:
	    ja[j] = -ja[j];
/* Computing 2nd power */
	    d__1 = a[j];
	    anormm += d__1 * d__1;
/* Computing 2nd power */
	    d__1 = a[j];
	    anormp += d__1 * d__1;
	    ++new__;
	    ++icg[i1];
L150:
	    ;
	}
	if (rowsum > deps) {
	    ++napos;
	} else if (rowsum < -deps) {
	    ++naneg;
	} else {
	    ++nazer;
	}
/* L200: */
    }
    anormm = sqrt(anormm);
    anormp = sqrt(anormp);
    asym = anormm / anormp;
    if (asym <= 1e-12) {
	asym = 0.;
    }

/* ===> MESSAGES ON A */

    if (asym == 0.) {
	   if (yes_print_amg) {
	      printf("CHECK: A PROBABLY SYMMETRIC\n");
	   }
    }
    if (asym != 0.) {
		 if (yes_print_amg) {
	         printf("CHECK: A PROBABLY NOT SYMMETRIC. MEASURE %1.5f\n",asym);
		 }
    }
    if (naoff > 0) {
	 if (yes_print_amg) {
#if doubleintprecision == 1
		 printf("CHECK: A PROBABLY NOT POS. TYPE: %lld OFF-DIAGONAL ELEMENTS POSITIVE\n", naoff);
#else
		 printf("CHECK: A PROBABLY NOT POS. TYPE: %d OFF-DIAGONAL ELEMENTS POSITIVE\n", naoff);
#endif
	    
	 }
    }
    if (naneg > 0) {
	 if (yes_print_amg) {
#if doubleintprecision == 1
		 printf("CHECK: A PROBABLY NOT POS. TYPE: %lld ROWSUMS NEGATIVE\n", naneg);
#else
		 printf("CHECK: A PROBABLY NOT POS. TYPE: %d ROWSUMS NEGATIVE\n", naneg);
#endif
	      
	 }
    }
    if (nazer == nnu) {
	 if (yes_print_amg) {
	      printf("CHECK: A PROBABLY SINGULAR - ROWSUMS ARE ZERO\n");
	 }
    }
    if (naoff == 0 && naneg == 0) {
	 if (yes_print_amg) {
	     printf("CHECK: A PROBABLY POSITIVE TYPE\n");
	 }
    }

/* ===> WARNINGS */

    if (*isym == 1 && asym != 0. || *isym > 1 && asym == 0. || *irow0 == 1 && 
	    nazer != nnu || *irow0 > 1 && nazer == nnu) {
	 if (yes_print_amg) {
         printf("--- WARNG IN CHECK: PARAM MATRIX MAY BE BAD ---");
	 }
	*ierr = -12;
    }

    if (new__ > 0) {
	goto L220;
    }
    if (yes_print_amg) {
	     printf("CHECK: MATRIX A WAS SYMMETRICALLY STORED\n");
	}
    goto L600;

/* ===> REPLACE A BY SYMMETRIZED VERSION */

L220:
    if (nna + new__ >= *nda) {
	goto L3210;
    }
    if (nna + new__ >= *ndja) {
	goto L3220;
    }

/* ===> EXTEND MATRIX A IN SITU */

    ifg[1] = 1;
    i__1 = nnu + 1;
    for (i__ = 2; i__ <= i__1; ++i__) {
	ifg[i__] = ifg[i__ - 1] + icg[i__ - 1];
/* L230: */
    }
    for (i__ = nnu; i__ >= 1; --i__) {
	ishift = ifg[i__] - ia[i__];
	i__1 = ia[i__];
	for (j = ia[i__ + 1] - 1; j >= i__1; --j) {
	    a[j + ishift] = a[j];
	    ja[j + ishift] = ja[j];
/* L240: */
	}
	icg[i__] = ia[i__ + 1] + ishift;
	ia[i__ + 1] = ifg[i__ + 1];
/* L250: */
    }

/* ===> SYMMETRIZE MATRIX A:  ISYM=1: COPY MISSING TRANSPOSE ENTRIES */
/*                           ISYM=2: FILL IN ZEROES */

    if (*isym != 1) {
	i__1 = nnu;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    j2 = icg[i__] - 1;
	    i__2 = j2;
	    for (j = ia[i__] + 1; j <= i__2; ++j) {
		i1 = ja[j];
		if (i1 > 0) {
		    goto L450;
		}
		ja[j] = -i1;
		jnew = icg[ja[j]];
		a[jnew] = 0.;
		ja[jnew] = i__;
		icg[ja[j]] = jnew + 1;
L450:
		;
	    }
/* L400: */
	}
	 if (yes_print_amg) {
	     printf("CHECK: STORAGE OF A HAS BEEN \n");
#if doubleintprecision == 1
		 printf("SYMMETRIZED BY INTRODUCING %lld ZEROES\n", new__);
#else
		 printf("SYMMETRIZED BY INTRODUCING %d ZEROES\n", new__);
#endif
	     
	 }
    } else {
	 if (yes_print_amg) {
#if doubleintprecision == 1
		 printf("*** WARNG IN CHECK: %lld A-ENTRIES MISSING ***\n", new__);
#else
		 printf("*** WARNG IN CHECK: %d A-ENTRIES MISSING ***\n", new__);
#endif
	    
	    printf("MISSING TRANSPOSE CONNECTIONS WILL BE FILLED IN\n");
	 }
	*ierr = -11;
	i__1 = nnu;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    j2 = icg[i__] - 1;
	    i__2 = j2;
	    for (j = ia[i__] + 1; j <= i__2; ++j) {
		i1 = ja[j];
		if (i1 > 0) {
		    goto L550;
		}
		ja[j] = -i1;
		jnew = icg[ja[j]];
		a[jnew] = 0.;
		ja[jnew] = i__;
		icg[ja[j]] = jnew + 1;
L550:
		;
	    }
/* L500: */
	}
    }

/* ===> REMOVE PAIRS OF ZEROES */

L600:
    trunc_(&c__1, &c__0, &a[1], &ia[1], &ja[1], &imin[1], &imax[1], &time[1], 
	    ierr, ium);
    return 0;

/* ===> ERROR MESSAGES */

L3200:
    //io___64.ciunit = *ium;
    //s_wsfe(&io___64);
    //e_wsfe();
	printf("*** ERROR IN CHECK: NDIG TOO SMALL ***\n");
    *ierr = 6;
    return 0;

L3210:
    //io___65.ciunit = *ium;
    //s_wsfe(&io___65);
    //e_wsfe();
	printf("*** ERROR IN CHECK: NDA TOO SMALL ***\n");
    *ierr = 1;
    return 0;

L3220:
    //io___66.ciunit = *ium;
    //s_wsfe(&io___66);
    //e_wsfe();
	printf("*** ERROR IN CHECK: NDJA TOO SMALL ***\n");
    *ierr = 3;
    return 0;

L3300:
    //io___67.ciunit = *ium;
    //s_wsfe(&io___67);
    //e_wsfe();
	printf("*** ERROR IN CHECK: POINTER IA ERRONEOUS***\n");
    *ierr = 15;
    return 0;

L3310:
    //io___68.ciunit = *ium;
    //s_wsfe(&io___68);
    //e_wsfe();
	printf("*** ERROR IN CHECK: POINTER JA ERRONEOUS***\n");
    *ierr = 16;
    return 0;

L3320:
    //io___69.ciunit = *ium;
    //s_wsfe(&io___69);
    //e_wsfe();
	printf("*** ERROR IN CHECK: DIAGONAL IS NOT STORED FIRST ***\n");
    *ierr = 13;
    return 0;

L3330:
    //io___70.ciunit = *ium;
    //s_wsfe(&io___70);
    //e_wsfe();
	printf("*** ERROR IN CHECK: DIAGONAL IS NON-POSITIVE ***\n");
    *ierr = 14;
    return 0;


} /* check_ */


/* ....................................................................... */

/*     TRUNC                                                 SUBROUTINE */

/* ....................................................................... */

/* Subroutine */ integer trunc_(integer *k, integer *ntr, doublereal *a, integer *
	ia, integer *ja, integer *imin, integer *imax, /*real*/unsigned int *time, integer *
	ierr, integer *ium)
{
    /* Format strings */
   // static char fmt_9000[] = "(\002 *** ERROR IN TRUNC: TRANSPOSE A-ENTRY MI"
	//    "SSING ON GRID\002,i3,\002 ***\002)";

    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Builtin functions */
   // integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);

    /* Local variables */
    integer i__=0, j=0, i1=0, j1=0, j2=0, j3=0;
    doublereal at=0.0;
    integer jt=0, imn=0, imx=0;
	integer nna=0;
	unsigned int told=0;
    integer jpos=0;
	unsigned int tnew=0;

    /* Fortran I/O blocks */
    //static cilist io___85 = { 0, 0, 0, fmt_9000, 0 };



/*     TRUNCATES OPERATOR ON GRID K CORRESPONDING TO THE VALUE OF NTR. */
/*     NTR HAS TO BE 0 OR 1: */

/*       =0:    PAIRS OF ZEROES ARE REMOVED FROM COARSE GRID OPERATORS; */
/*       =1:    COARSE GRID OPERATORS REMAIN UNCHANGED. */



    /* Parameter adjustments */
    --time;
    --imax;
    --imin;
    --ja;
    --ia;
    --a;

    /* Function Body */
    if (*ntr == 1) {
	return 0;
    }

	told=clock();
    imn = imin[*k];
    imx = imax[*k];
    nna = ia[imx + 1] - ia[imn];
    jpos = ia[imn];

    i__1 = imx;
    for (i__ = imn; i__ <= i__1; ++i__) {
	j1 = ia[i__] + 1;
	j2 = ia[i__ + 1] - 1;
	a[jpos] = a[ia[i__]];
	ja[jpos] = i__;
	ia[i__] = jpos;
	++jpos;
	i__2 = j2;
	for (j = j1; j <= i__2; ++j) {
	    i1 = ja[j];
	    if (i1 < 0) {
		goto L250;
	    }
	    if (i1 < i__) {
		goto L230;
	    }
	    if (a[j] != 0.) {
		goto L230;
	    }
	    i__3 = ia[i1 + 1] - 1;
	    for (j3 = ia[i1]; j3 <= i__3; ++j3) {
		if (ja[j3] != i__) {
		    goto L210;
		}
		jt = j3;
		at = a[j3];
		goto L215;
L210:
		;
	    }
	    goto L1000;
L215:
	    if (at != 0.) {
		goto L230;
	    }
	    ja[jt] = -ja[jt];
	    goto L250;
L230:
	    a[jpos] = a[j];
	    ja[jpos] = ja[j];
	    ++jpos;
L250:
	    ;
	}
/* L205: */
    }
    ia[imx + 1] = jpos;

/* ===> EXIT */

	tnew=clock();
    time[7] = time[7] + tnew - told;
    return 0;

/* ===> ERROR MESSAGE */

L1000:
    //io___85.ciunit = *ium;
    //s_wsfe(&io___85);
    //do_fio(&c__1, (char *)&(*k), (ftnlen)sizeof(integer));
    //e_wsfe();
#if doubleintprecision == 1
	printf("*** ERROR IN TRUNC: TRANSPOSE A-ENTRY MISSING ON GRID %lld ***\n", k[0]);
#else
	printf("*** ERROR IN TRUNC: TRANSPOSE A-ENTRY MISSING ON GRID %d ***\n", k[0]);
#endif
	
    *ierr = 21;
    return 0;

} /* trunc_ */

/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */

/*     AMG1R5 SETUP ROUTINES (VERSION 3) */

/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
/* ....................................................................... */

/*     CRSNG                                           SUBROUTINE */

/* ....................................................................... */

/* Subroutine */ integer crsng_(integer *levelx, doublereal *eecg1, doublereal *
	eecg2, doublereal *eewt2, integer *nnwt, integer *ntr, integer *ierr, 
	doublereal *a, doublereal *u, integer *ia, integer *ja, integer *iw, 
	integer *imin, integer *imax, integer *iminw, integer *imaxw, integer 
	*icg, integer *ifg, integer *nstcol, integer *levels, integer *irow0, 
	integer *nda, integer *ndja, integer *ndia, integer *ndu, integer *
	ndf, integer *ndicg, /*real*/ unsigned int *time, integer *ium, integer *mda, integer *
	mdia, integer *mdja, integer *mdu, integer *mdf, integer *mdig, 
	integer *maxgr, integer *lratio)
{
    /* Format strings */
    //static char fmt_9000[] = "(/\002 **************** SPACE REQUIREMENTS ***"
	//    "*************\002//\002 VECTOR          NEEDED                  "
	//    "            \002/\002 ------------------------------------------"
	//    "----------\002/\002    A \002,i16,\002   ADJUST THE DIMENSIONING"
	//    " OF \002/\002    JA\002,i16,\002   VECTORS A - IG IN THE     "
	//    " \002/\002    IA\002,i16,\002   CALLING PROGRAM ACCORDING  \002"
	//    "/\002    U \002,i16,\002   TO THE CALCULATED SPACE RE-\002/\002 "
	//    "   F \002,i16,\002   QUIREMENTS AND RERUN THE   \002/\002    I"
	//    "G\002,i16,\002   PROGRAM.                   \002/\002 ----------"
	//    "------------------------------------------\002)";
   // static char fmt_9010[] = "(/\002 NOTE: IF YOU WANT TO USE CG-CORRECTIONS"
   //	    " IN THE SOLU-\002/\002       TION PROCESS (NCYC-SUBPARAMETER ICG"
   //	    "R=1 OR =2),\002/\002       PROVIDE FOR ADDITIONAL\002,i6,\002 ST"
   //	    "ORAGE LOCATIONS\002/\002       IN VECTORS U AND F.              "
   //	    "             \002/\002         SIMILARLY, USAGE OF THE YALE-SMP "
   //	    "SOLVER ON  \002/\002       THE COARSEST GRID (NSOLCO=2) WILL REQ"
   //	    "UIRE     \002/\002       ADDITIONAL SPACE IN VECTOR A DURING THE"
   //	    " SOLU- \002/\002       TION PHASE. IN THIS CASE, HOWEVER, ITS EX"
   //	    "ACT  \002/\002       AMOUNT ISN'T PREDICTABLE.                  "
   //	    "  \002/)";

    /* System generated locals */
    integer i__1=0, i__2=0, i__3=0;

    /* Builtin functions */
   // integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);

    /* Local variables */
    integer i__=0, k=0, nwt=0;
    doublereal ecg1=0.0, ecg2=0.0, ewt2=0.0;
    integer ichk=0, iias=0, isia=0, mdir=0;
    extern /* Subroutine */ integer pcol_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, doublereal /*integer*/ *, doublereal /*integer*/ *, /*real*/ unsigned int *, integer *, integer *,
	     integer *, integer *, integer *, integer *);
    integer mdiw=0, mmax=0, kerr=0, ndir=0, iirs=0, irst=0, mdicg=0, iajas=0, isaja=0;
    extern /* Subroutine */ integer opdfn_(integer *, integer *, integer *, 
	    doublereal *, integer *, integer *, integer *, integer *, integer 
	    *, integer *, integer *, integer *, integer *, doublereal /*integer*/ *, 
	    integer *, integer *, /*real*/unsigned int *, integer *, integer *, integer *, 
	    integer *, integer *, integer *);
    integer mdjtr=0, ndjtr=0, ncolx=0;
    extern /* Subroutine */ integer trunc_(integer *, integer *, doublereal *, 
	    integer *, integer *, integer *, integer *, /*real*/unsigned int *, integer *, 
	    integer *), pwint_(integer *, doublereal *, integer *, integer *, 
	    doublereal *, integer *, integer *, integer *, integer *, integer 
	    *, integer *, integer *, integer *, integer *, doublereal /*integer*/ *, /*real*/ unsigned int 
	    *, integer *, integer *, integer *, integer *, integer *, integer 
	    *, integer *, integer *, integer *, integer *, integer *);
    integer jtrst=0;
    extern /* Subroutine */ integer rwsrt_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, integer *, integer 
	    *, integer *, integer *, integer *, doublereal /*integer*/ *, /*real*/ unsigned int *, integer 
	    *, integer *, integer *), setifg_(integer *, integer *, integer *,
	     integer *, integer *, integer *, /*real*/ unsigned int *);
    integer kfirst=0;

    /* Fortran I/O blocks */
   // static cilist io___110 = { 0, 0, 0, fmt_9000, 0 };
   // static cilist io___111 = { 0, 0, 0, fmt_9010, 0 };



/*     PERFORMS COARSENING, DEFINES INTERPOLATIONS AND CGC-OPERATORS */

/*     ============== STANDARD VALUES OF PARAMETERS ===================== */

/*             ECG1=0.,    ECG2=0.25,    EWT2=0.35,   NWT=2 */

/*     ================ DESCRIPTION OF PARAMETERS ======================= */

/*     ECG1 --    DEFINES CRITERION FOR DETERMINING DIAGONAL */
/*                DOMINANCE IN RWSRT. I.E., IF THE ABSOLUTE */
/*                VALUE OF THE SUM OF THE OFF-DIAGONALS OF ROW */
/*                I IS SMALLER THAN ECG1 TIMES THE ABSOLUTE */
/*                VALUE OF THE DIAGONAL ENTRY, THEN TOCHKA I IS */
/*                IMMEDIATELY FORCED TO BE AN F-TOCHKA IN THE */
/*                PRE-COLORING ALGORITHM PCOL. */
/*                IN THE SECOND PART (WINT), NO INTERPOLATION */
/*                IS DEFINED FOR TOCHKA I, AND NO POINTS USE */
/*                I FOR INTERPOLATION. IN ADDITION, THE WEIGHT */
/*                FOR TOCHKA I IS NOT DISTRIBUTED TO OTHER POINTS */
/*                WHEN DEFINING THE INTERPOLATION WEIGHTS FOR */
/*                POINTS WHICH DEPEND ON TOCHKA I. (THIS IS */
/*                EQUIVALENT TO COMPLETELY IGNORING SUCH POINTS */
/*                IN DETERMINING THE COARSE GRID AND INTERPOLATION */
/*                WEIGHTS. */

/*     ECG2 --    DEFINES STRONG CONNECTIONS (ALPHA IN THE PAPER) */

/*     EWT2 --    DEFINES STRONG DEPENDENCE ON A SET (BETA IN THE PAPER) */

/*     NWT  --    PARAMETER CONTROLLING THE DEFINITION OF INTERPOLATION */
/*                FORMULAS: */
/*                  =1 - CHECKING OF INT-FORMULA: OFF */
/*                  =2 - CHECKING OF INT-FORMULA: ON */


/* ===> ASSIGN DEFAULT VALUES IF ZERO */

    /* Parameter adjustments */
    --time;
    --nstcol;
    --ifg;
    --icg;
    --imaxw;
    --iminw;
    --imax;
    --imin;
    --iw;
    --ja;
    --ia;
    --u;
    --a;

    /* Function Body */
    nwt = *nnwt;
    ecg1 = *eecg1;
    ecg2 = *eecg2;
    ewt2 = *eewt2;
    if (nwt == 0) {
	nwt = 2;
    }
    if (ecg2 == 0.) {
	ecg2 = .25;
    }
    if (ewt2 == 0.) {
	ewt2 = .35;
    }

/* ===> DECODE PARAMETER NWT */

/* L5: */
    ichk = nwt;

/* ===> COARSENING */

/* L10: */
    *levels = my_imin(*levelx,*maxgr);
    mmax = *maxgr;

/* ===> INITIALIZE PARAMETERS MDA - MDIW, LATER SET TO ACTUAL STORAGE */
/*     REQUIREMENTS OF CORRESPONDING VECTORS */

    *mda = 0;
    *mdja = 0;
    mdicg = 0;
    mdir = 0;
    mdiw = imax[1];

/* ===> KFIRST IS THE NUMBER OF THE FIRST STORED GRID, IAJAS, IIAS AND */
/*     IIRS ARE THE SHIFTS IN VECTORS A, JA, IA AND IR, RESPECTIVELY, AS */
/*     COMPARED TO FULL STORAGE OF ALL GRIDS */

    kerr = 0;
    kfirst = 1;
    iajas = 0;
    iias = 0;
    iirs = 0;

    i__1 = *levels;
    for (k = 2; k <= i__1; ++k) {

/* ===> JTRST: INITIAL POINTER FOR WORK SPACE IN VECTOR A, TO CONTAIN */
/*            THE STRONG TRANSPOSE CONNECTIONS */
/*     NDJTR: AVAILABLE WORK SPACE */

L20:
	jtrst = ia[imax[k - 1] + 1];
	ndjtr = *lratio * (*nda - jtrst + 1);
	i__2 = k - 1;

	

	rwsrt_(&i__2, &ecg1, &ecg2, ierr, &a[1], &ia[1], &ja[1], &iw[1], &
		imin[1], &imax[1], &imaxw[1], &ifg[1],  &a[jtrst], &time[1], &
		ndjtr, ium, &mdjtr);

	

	if (*ierr > 0) {
	    if (*ierr >= 1 && *ierr <= 6) {
		goto L30;
	    }
	    return 0;
	}
/* Computing MAX */
	i__2 = *mda, i__3 = jtrst + (mdjtr - 1) / *lratio + iajas;
	*mda = myi_max(i__2,i__3);

/* ===>   IRST: INITIAL POINTER FOR WORK SPACE IN VECTOR U, TO CONTAIN */
/*             RESET STACK */
/*       NDIR: AVAILABLE WORK SPACE */

	irst = mdiw + 1 - iirs;
	ndir = *lratio * (*ndu - irst + 1);
	i__2 = k - 1;

	

	

	pcol_(&i__2, ierr, &ia[1], &ja[1], &iw[1], &imin[1], &imax[1], &imaxw[
		1], &icg[1], &ifg[1],  &u[irst] , &a[jtrst], &time[1], ndicg, &
		ndir, ium, &mdicg, &mdir, &iias);

   

	if (*ierr > 0) {
	    if (*ierr >= 1 && *ierr <= 6) {
		goto L30;
	    }
	    return 0;
	}
	i__2 = k - 1;

	

	pwint_(&i__2, &ewt2, &ichk, &mmax, &a[1], &ia[1], &ja[1], &iw[1], &
		imin[1], &imax[1], &iminw[1], &imaxw[1], &ifg[1], &icg[1], &u[irst],
		&time[1], ierr, irow0, &ncolx, nda, ndja, ndicg, ndia, 
		ium, mda, mdja, &iajas);

		

	if (*ierr > 0) {
	    if (*ierr >= 1 && *ierr <= 6) {
		goto L30;
	    }
	    return 0;
	}

	

	opdfn_(&k, ierr, &mmax, &a[1], &ia[1], &ja[1], &iw[1], &imin[1], &
		imax[1], &iminw[1], &imaxw[1], &icg[1], &ifg[1],  &u[irst], &
		nstcol[1], &ncolx, &time[1], nda, ndja, ium, mda, mdja, &
		iajas);

	

	if (*ierr > 0) {
	    if (*ierr >= 1 && *ierr <= 6) {
		goto L30;
	    }
	    return 0;
	}
	if (k > mmax) {
	    goto L100;
	}
	trunc_(&k, ntr, &a[1], &ia[1], &ja[1], &imin[1], &imax[1], &time[1], 
		ierr, ium);
	if (*ierr > 0) {
	    return 0;
	}
	iw[iminw[k]] = ia[imax[k] + 1];
	if (k >= mmax) {
	    goto L100;
	}
	goto L50;
L30:
	if (k <= kfirst + 1 && iirs != 0) {
	    return 0;
	}
	kerr = *ierr;
	*ierr = 0;
	kfirst = k - 1;
	iirs = mdiw;
	isia = imin[k - 1] - 1;
	iias += isia;
	isaja = ia[imin[k - 1]] - 1;
	iajas += isaja;
	i__2 = ia[imax[k - 1] + 1] - 1;
	for (i__ = ia[imin[k - 1]]; i__ <= i__2; ++i__) {
	    ja[i__ - isaja] = ja[i__] - isia;
	    a[i__ - isaja] = a[i__];
/* L35: */
	}
	i__2 = imax[k - 1] + 1;
	for (i__ = imin[k - 1]; i__ <= i__2; ++i__) {
	    ia[i__ - isia] = ia[i__] - isaja;
/* L40: */
	}
	imin[k - 1] = 1;
	imax[k - 1] -= isia;
	iw[iminw[k - 1]] = ia[imax[k - 1] + 1];
	goto L20;
L50:
	;
    }
    goto L200;
L100:
    *levels = mmax;
L200:
    *mdia = imax[*levels] + 1 + iias;
/* Computing MAX */
    i__1 = imax[*levels] + iias, i__2 = mdiw + 1 + (mdir - 1) / *lratio, i__1 
	    = myi_max(i__1,i__2), i__2 = mdiw + 1 + (mdiw - 1) / *lratio;
    *mdu = myi_max(i__1,i__2);
    *mdf = imax[*levels] + iias;
    *mdig = mdiw + 2 + (mdicg << 1);
    if (kerr != 0 || *mdu > *ndu || *mdf > *ndf) {
	//io___110.ciunit = *ium;
	//s_wsfe(&io___110);
	//do_fio(&c__1, (char *)&(*mda), (ftnlen)sizeof(integer));
	//do_fio(&c__1, (char *)&(*mdja), (ftnlen)sizeof(integer));
	//do_fio(&c__1, (char *)&(*mdia), (ftnlen)sizeof(integer));
	//do_fio(&c__1, (char *)&(*mdu), (ftnlen)sizeof(integer));
	//do_fio(&c__1, (char *)&(*mdf), (ftnlen)sizeof(integer));
	//do_fio(&c__1,(char *)&(*mdig) , (ftnlen)sizeof(integer));
	//e_wsfe();

		if (yes_print_amg) 
		{
#if doubleintprecision == 1
			printf("( **************** SPACE REQUIREMENTS ***\n");
			printf("************* VECTOR          NEEDED                  \n");
			printf("             -----------------------------------------\n");
			printf("----------    A ,%lld,   ADJUST THE DIMENSIONING \n", mda[0]);
			printf(" OF     JA ,%lld,   VECTORS A - IG IN THE     \n", mdja[0]);
			printf("     IA,%lld,   CALLING PROGRAM ACCORDING  \n", mdia[0]);
			printf("    U ,%lld,   TO THE CALCULATED SPACE  \n", mdu[0]);
			printf("  REF ,%lld,   QUIREMENTS AND RERUN THE      \n", mdf[0]);
			printf("IG %lld ,   PROGRAM.                    ----------\n", mdig[0]);
			printf("------------------------------------------)\n");
#else
			printf("( **************** SPACE REQUIREMENTS ***\n");
			printf("************* VECTOR          NEEDED                  \n");
			printf("             -----------------------------------------\n");
			printf("----------    A ,%d,   ADJUST THE DIMENSIONING \n", mda[0]);
			printf(" OF     JA ,%d,   VECTORS A - IG IN THE     \n", mdja[0]);
			printf("     IA,%d,   CALLING PROGRAM ACCORDING  \n", mdia[0]);
			printf("    U ,%d,   TO THE CALCULATED SPACE  \n", mdu[0]);
			printf("  REF ,%d,   QUIREMENTS AND RERUN THE      \n", mdf[0]);
			printf("IG %d ,   PROGRAM.                    ----------\n", mdig[0]);
			printf("------------------------------------------)\n");
#endif
	         
		}
	
	//io___111.ciunit = *ium;
	//s_wsfe(&io___111);
	//do_fio(&c__1, (char *)&mdiw, (ftnlen)sizeof(integer));
	//e_wsfe();
		if (yes_print_amg) {
	        printf("NOTE: IF YOU WANT TO USE CG-CORRECTIONS\n");
        	printf(" IN THE SOLUTION PROCESS (NCYC-SUBPARAMETER ICG\n");
#if doubleintprecision == 1
			printf("R=1 OR =2),  PROVIDE FOR ADDITIONAL %lld ST\n", mdiw);
#else
			printf("R=1 OR =2),  PROVIDE FOR ADDITIONAL %d ST\n", mdiw);
#endif
	        
	        printf("ORAGE LOCATIONS IN VECTORS U AND F. \n");
	        printf("SIMILARLY, USAGE OF THE YALE-SMP \n");
	        printf("SOLVER ON  THE COARSEST GRID (NSOLCO=2) WILL REQUIRE\n");
	        printf(" ADDITIONAL SPACE IN VECTOR A DURING THE\n");
        	printf(" SOLUTION PHASE. IN THIS CASE, HOWEVER, ITS \n");
	        printf("EXACT  AMOUNT ISN'T PREDICTABLE.   \n");
		}    

	*ierr = kerr;
	return 0;
    }
    setifg_(&imin[1], &imax[1], &icg[1], &ifg[1], &nstcol[1], levels, &time[1]
	    );
    return 0;
} /* crsng_ */


/* ....................................................................... */

/*     RWSRT                                             SUBROUTINE */

/* ....................................................................... */

/* Subroutine */ integer rwsrt_(integer *k, doublereal *ecg1, doublereal *ecg2, 
	integer *ierr, doublereal *a, integer *ia, integer *ja, integer *iw, 
	integer *imin, integer *imax, integer *imaxw, integer *ifg, /*integer*/ doublereal *
	jtr, /*real*/ unsigned int *time, integer *ndjtr, integer *ium, integer *mdjtr)
{
    /* Format strings */
    //static char fmt_9910[] = "(\002 *** ERROR IN RWSRT: NDA TOO SMALL ***"
	//    "\002)";

    /* System generated locals */
    integer i__1=0, i__2=0;
    doublereal d__1=0.0;

    /* Builtin functions */
   // integer s_wsfe(cilist *), e_wsfe(void);

    /* Local variables */
    integer i__=0, j=0, ii=0;
    doublereal rs=0.0;
    integer ihi=0, jhi=0;
    doublereal amn=0.0;
    integer ilo=0, jlo=0;
    doublereal amx=0.0, ast=0.0;
    integer imx=0, jmx=0, iws=0;
	unsigned int told=0;
    doublereal atmp=0.0;
    integer itmp=0;
	unsigned int tnew=0;

    /* Fortran I/O blocks */
   // static cilist io___130 = { 0, 0, 0, fmt_9910, 0 };



/*     ROW-SORT ALGORITHM FOR ROWS OF A(K). IN DETAIL: */

/*     - SORTS THE ELEMENTS OF EACH ROW OF A(K) SUCH THAT THE STRONG */
/*       CONNECTIONS JUST FOLLOW THE DIAGONAL ELEMENT. ONE OF THE STRONG- */
/*       EST IS ALWAYS FIRST (EVEN IF PARTIAL SORTING IS PERFORMED!). */
/*       NON-NEGATIVE CONNECTIONS ARE ALWAYS DEFINED TO BE WEAK. */
/*       JA(IA(I)) IS RE-DEFINED TO TOCHKA TO THE LAST STRONG CONNECTION */
/*       OF TOCHKA I (OR TO IA(I) IF THERE IS NO STRONG CONNECTION). */

/*     - THE STRONG TRANSPOSE CONNECTIONS ARE LOADED LOGICALLY */
/*       INTO JTR, I.E., ITS POINTERS ARE STORED IN JTR. */

/*     ============== COMMENTS ON WORK SPACE USED ======================= */

/*     - JA(IA(I)) ----- I=IMIN(K),...,IMAX(K) */

/*     IS DEFINED TO TOCHKA TO THE LAST STRONG CONNECTION OF TOCHKA I. */
/*     NOTE THAT THE ORIGINAL CONTENTS OF JA(IA(I)) IS OVERWRITTEN. */
/*     (IT IS PUT BACK IN SUBROUTINE WINT.) */

/*     - JTR(J)--------- J=1,IW(IMAX(K)+IWS+1)-1 */
/*     - IW(I) --------- I=IMIN(K)+IWS,...,IMAX(K)+IWS */
/*     - IFG(I)--------- I=IMIN(K),...,IMAX(K) */

/*     (IWS=0, IF K=1; IWS=IMAXW(K-1)+2-IMIN(K) OTHERWISE) */

/*     JTR IS INITIALIZED TO HAVE SAME FORM AS JA. JTR(J) */
/*     CONTAINS INFORMATION ON STRONG TRANSPOSE CONNECTIONS: */
/*     JTR(J) WITH IW(I+IWS)<=J<=IW(I+IWS+1)-1 POINTS TO THE STRONG */
/*     TRANSPOSE CONNECTIONS OF I. */


/* ===> INITIALIZATION OF WORK SPACE */

    /* Parameter adjustments */
    --time;
    --jtr;
    --ifg;
    --imaxw;
    --imax;
    --imin;
    --iw;
    --ja;
    --ia;
    --a;

    /* Function Body */
	told=clock();
    ilo = imin[*k];
    ihi = imax[*k];
    if (*k != 1) {
	iws = imaxw[*k - 1] + 2 - ilo;
    } else {
	iws = 0;
    }
    i__1 = ihi;
    for (i__ = ilo; i__ <= i__1; ++i__) {
	ifg[i__] = 0;
	ja[ia[i__]] = ia[i__];
/* L5: */
    }
	tnew=clock();
    time[8] = time[8] + tnew - told;
    told = tnew;

/*     ************************* */
/*     * PARTIAL STANDARD SORT * */
/*     ************************* */

/*     NOTE: NON-NEGATIVE CONNECTIONS ARE ALWAYS DEFINED TO BE WEAK! THE */
/*           STRONGEST CONNECTION IS ALWAYS FOLLOWING THE DIAGONAL. */

    i__1 = ihi;
    for (i__ = ilo; i__ <= i__1; ++i__) {
	jlo = ia[i__] + 1;
	jhi = ia[i__ + 1] - 1;
	if (jhi < jlo) {
	    goto L590;
	}

/* ===>   FIND STRONGEST CONNECTION AND SUM OF |OFF-DIAGONALS| OF ROW I */

	amx = a[jlo];
	amn = a[jlo];
	jmx = jlo;
	if (*ecg1 != 0.) {
	    rs = 0.;
	    i__2 = jhi;
	    for (j = jlo + 1; j <= i__2; ++j) {
		rs += (d__1 = a[j], fabs(d__1));
		if (a[j] < amx) {
		    amx = a[j];
		    jmx = j;
		} else if (a[j] > amn) {
		    amn = a[j];
		}
/* L550: */
	    }

/* ===>     TEST FOR POSITIVE OFF-DIAGONALS / DIAGONAL DOMINANCE */

	    if (amx >= 0. || rs <= *ecg1 * a[ia[i__]]) {
		goto L590;
	    }
	} else {
	    i__2 = jhi;
	    for (j = jlo + 1; j <= i__2; ++j) {
		if (a[j] < amx) {
		    amx = a[j];
		    jmx = j;
		} else if (a[j] > amn) {
		    amn = a[j];
		}
/* L555: */
	    }

/* ===>   TEST FOR POSITIVE OFF-DIAGONALS */

	    if (amx >= 0.) {
		goto L590;
	    }
	}

/* ===>   PUT STRONGEST CONNECTION IN FIRST POSITION */

	ast = *ecg2 * amx;
	imx = ja[jmx];
	a[jmx] = a[jlo];
	ja[jmx] = ja[jlo];
	a[jlo] = amx;
	ja[jlo] = imx;
	if (amn <= ast) {
	    goto L580;
	}
	++jhi;

/* ===>   DECREASE JHI UNTIL A STRONG CONNECTION IS FOUND */
/*       (IF JLO >= JHI STOP: ALL CONNECTIONS ARE SORTED) */

L560:
	--jhi;
	if (jlo >= jhi) {
	    goto L580;
	}
	if (a[jhi] > ast) {
	    goto L560;
	}

/* ===>   INCREASE JLO UNTIL A WEAK CONNECTION IS FOUND */
/*       (IF JLO >= JHI STOP: ALL CONNECTIONS ARE SORTED) */

L570:
	++jlo;
	if (jlo >= jhi) {
	    goto L580;
	}
	if (a[jlo] <= ast) {
	    goto L570;
	}

/* ===>   INTERCHANGE A(JHI) AND A(JLO) */

	atmp = a[jhi];
	itmp = ja[jhi];
	a[jhi] = a[jlo];
	ja[jhi] = ja[jlo];
	a[jlo] = atmp;
	ja[jlo] = itmp;
	goto L560;

/* ===>   ROW SORTED --  SET JA(IA(I)) TO LAST STRONG CONNECTION */

L580:
	ja[ia[i__]] = jhi;

/* ===>   COUNT STRONG TRANSPOSE CONNECTIONS IN ROW I */

	i__2 = jhi;
	for (j = ia[i__] + 1; j <= i__2; ++j) {
	    ++ifg[ja[j]];
/* L585: */
	}
L590:
	;
    }

/* ===> INITIALIZATION OF WORK SPACE FOR STRONG TRANSPOSE CONNECTIONS */

    iw[ilo + iws] = 1;
    i__1 = ihi;
    for (i__ = ilo; i__ <= i__1; ++i__) {
	iw[i__ + iws + 1] = iw[i__ + iws] + ifg[i__];
	ifg[i__] = iw[i__ + iws];
/* L5010: */
    }
    *mdjtr = iw[ihi + iws + 1] - 1;
    if (*mdjtr > *ndjtr) {
	goto L9901;
    }

/* ===> LOAD POINTERS TO STRONG TRANSPOSE CONNECTIONS INTO JTR */

    i__1 = ihi;
    for (i__ = ilo; i__ <= i__1; ++i__) {
	i__2 = ja[ia[i__]];
	for (j = ia[i__] + 1; j <= i__2; ++j) {
	    ii = ja[j];
	   // jtr[ifg[ii]] = i__;
		jtr[ifg[ii]]=static_cast<doublereal>(i__);
	    ++ifg[ii];
/* L5020: */
	}
/* L5030: */
    }
	tnew=clock();
    time[1] = time[1] + tnew - told;
    return 0;

/* ===> ERROR MESSAGES */

L9901:
    //io___130.ciunit = *ium;
    //s_wsfe(&io___130);
    //e_wsfe();
	printf("*** ERROR IN RWSRT: NDA TOO SMALL ***\n");
    *ierr = 1;
    return 0;

} /* rwsrt_ */


/* ....................................................................... */

/*     PCOL                                            SUBROUTINE */

/* ....................................................................... */

/* Subroutine */ integer pcol_(integer *k, integer *ierr, integer *ia, integer *
	ja, integer *iw, integer *imin, integer *imax, integer *imaxw, 
	integer *icg, integer *ifg, /*integer*/ doublereal *ir, /*integer*/ doublereal *jtr, /*real*/ unsigned int *time, 
	integer *ndicg, integer *ndir, integer *ium, integer *mdicg, integer *
	mdir, integer *iias)
{
    /* Format strings */
   // static char fmt_9910[] = "(\002 *** ERROR IN PCOL: NDIG TOO SMALL ***"
	//    "\002)";
   // static char fmt_9920[] = "(\002 *** ERROR IN PCOL: NDU TOO SMALL ***\002)"
	    ;

    /* System generated locals */
    integer i__1=0, i__2=0;

//#if AMG1R6_LABEL==1

	integer i__3 = 0;

//#endif

    /* Builtin functions */
    //integer s_wsfe(cilist *), e_wsfe(void);

    /* Local variables */
    integer i__=0, j=0, ic=0, ii=0, jj=0, in=0, jv=0, iic=0, ihi=0, iii=0, iip=0, ilo=0, iws=0, 
	    ilo1=0, iiii=0;
	unsigned int told=0;
//#if AMG1R6_LABEL==1
	integer itop = 0;
//#else
    //integer nscn=0, itop=0;
	integer nscn = 0;
//#endif
	unsigned int tnew=0;
    integer npts=0, jval0=0, npts1=0;
    integer jvalx=0, jcnbhi=0, ntrlim=0;
	integer nscnmx=0;

    /* Fortran I/O blocks */
    //static cilist io___157 = { 0, 0, 0, fmt_9910, 0 };
    //static cilist io___158 = { 0, 0, 0, fmt_9920, 0 };



/*     PRE-COLORING ALGORITHM FOR GRID K. THIS IS THE VERSION AS */
/*     DESCRIBED IN RUGE/STUEBEN (BRISTOL). THE GOAL IS TO OBTAIN QUICKLY */
/*     A TENTATIVE SET OF C-POINTS WITH THE FOLLOWING PROPERTIES: */

/*       - THE C-POINTS ARE ONLY WEAKLY CONNECTED AMONG EACH OTHER; */
/*       - EACH F-TOCHKA HAS (AT LEAST) ONE STRONG CONNECTION TO A C-TOCHKA */
/*         (EXCEPT FOR SOME EXCEPTIONAL F-POINTS, E.G., THOSE WHICH DO */
/*         NOT HAVE A STRONG CONNECTION AT ALL: THE "FORCED F-POINTS"). */

/*     ON EXIT, ICG(I) (IMIN(K)<=I<=IMAX(K)) CONTAINS ALL THE INFORMATION */
/*     ON THE COLORING. IN DETAIL: */

/*        ICG(I) = 1: I IS C-TOCHKA; */
/*        ICG(I) = 0: I IS FORCED F-TOCHKA, I.E., IT IS AN F-TOCHKA WHICH */
/*                    HAS NO STRONG CONNECTION AT ALL; */
/*        ICG(I) =-1: I IS "REGULAR" F-TOCHKA, I.E., I HAS AT LEAST ONE */
/*                    STRONG CONNECTION TO A C-TOCHKA; */
/*        ICG(I) =-2: I IS "EXCEPTIONAL" F-TOCHKA, I.E., I HAS NO STRONG */
/*                    CONNECTION TO A C-TOCHKA. FURTHERMORE, NO REGULAR */
/*                    F-TOCHKA IS STRONGLY CONNECTED TO I AND ALL POINTS */
/*                    WITH ICG=-2 HAVE NO STRONG CONNECTION AMONG EACH */
/*                    OTHER. (NOTE: THESE POINTS DON'T CONTRIBUTE FROM */
/*                    ANY C-TOCHKA. ALSO, IT DOES NOT MAKE SENSE TO MAKE */
/*                    THESE POINTS C-POINTS AS NO F-TOCHKA CONTRIBUTES */
/*                    FROM THEM.) */

/*     ============ COMMENTS ON INPUT =================================== */

/*     - JA(IA(I)) ----- I=IMIN(K),...,IMAX(K) */

/*     USED AS DEFINED IN RWSRT. NOT CHANGED. */

/*     - JTR(J)--------- J=1,IW(IMAX(K)+IWS+1)-1 */
/*     - IW(I) --------- I=IMIN(K)+IWS,...,IMAX(K)+IWS */

/*     (IWS=0, IF K=1; IWS=IMAXW(K-1)+2-IMIN(K) OTHERWISE) */

/*     USED AS DEFINED IN RWSRT. NOT CHANGED. */

/*     ============== COMMENTS ON WORK SPACE USED ======================= */

/*     - ICG(I) -------- I=IMIN(K),...,IMAX(K) */

/*     USED TO DISTINGUISH F-, C- AND U- (UNDECIDED) POINTS. */
/*     ON ENTRY, ICG IS DEFINED TO BE 0 FOR FORCED F-POINTS AND -2 */
/*     FOR U-PNTS. ON EXIT: SEE ABOVE. */

/*     - IFG(I) -------- I=IMIN(K),...,IMAX(K) */

/*     DEFINED TO BE A MEASURE FOR IMPORTANCE OF MAKING THE U-PNT I */
/*     A C-TOCHKA. THIS MEASURE IS A VALUE BETWEEN JVAL0 AND JVALX WITH */

/*       JVAL0=IMAX(K)+NPTS+2,   NPTS=IMAX(K)-IMIN(K)+1; */
/*       JVALX=JVAL0+4*NTRAV+1,  NTRAV=AVVERAGE NUMBER OF STRONG TRANS- */
/*                                     POSE CONNECTIONS */

/*     THE HIGHER THIS VALUE, THE MORE IMPORTANT IS IT TO MAKE U-TOCHKA I */
/*     A C-TOCHKA. ALSO, THE VALUE JV=IFG(I) IS JUST THE "ORIGIN" OF */
/*     A LIST WHICH CONNECTS ALL POINTS HAVING THE SAME MEASURE JV OF */
/*     IMPORTANCE. */

/*     - ICG(II) -------- II=IMAX(K)+2,...,JVAL0-1 */
/*     - IFG(II) -------- II=IMAX(K)+2,...,JVAL0-1 */

/*     LEFT AND RIGHT STACK POINTERS IN A DOUBLY LINKED LIST. THERE */
/*     ARE SEVERAL SUCH LISTS, EACH OF THEM LINKING U-POINTS WHICH */
/*     HAVE THE SAME MEASURE OF IMPORTANCE JV TO BE MADE C-POINTS. */
/*     EACH OF THE VALUES JV MAY BE REGARDED AS THE "ORIGIN" OF SUCH */
/*     A LIST DESCRIBED IN THE FOLLOWING. */

/*     - ICG(JV) -------- JV=JVAL0,...,JVALX */
/*     - IFG(JV) -------- JV=JVAL0,...,JVALX */

/*     USED AS POINTERS INTO THE LISTS TO INDICATE ITS BEGINNING AND ITS */
/*     END: IFG(JV) POINTS TO THE FIRST TOCHKA IN THE LIST AND ICG(JV) */
/*     POINTS TO THE LAST ONE. THE ROUGH PICTURE OF THE LIST WHICH */
/*     CORRESPONDS TO A VALUE JVAL0 <= JV <= JVALX LOOKS AS  FOLLOWS: */


/*                                  IFG */
/*       ---------------------------------------------------------- */
/*       |                                                        | */
/*       |            IFG           IFG           IFG             | */
/*       ---> ------ ----> ------- ----> ------- ----> ------- ---- */
/*            |    |       |     |       |     |       |     | */
/*            | JV |       | II1 |       | II2 |       | II3 | */
/*            |    |  ICG  |     |  ICG  |     |  ICG  |     | */
/*       ---- ------ <---- ------- <---- ------- <---- ------- <--- */
/*       |                                                        | */
/*       |                          ICG                           | */
/*       ---------------------------------------------------------- */


/*     HERE, II1, II2,... LOGICALLY REPRESENT THE "PHYSICAL" POINTS */
/*     I1, I2,... (WHICH ARE VALUES BETWEEN IMIN(K) AND IMAX(K)). */
/*     THE INTERCONNECTION BETWEEN THESE TWO REPRESENTATIONS OF THE */
/*     SAME POINTS IS GIVEN BY THE RELATION */

/*                          II:= I+NPTS+1. */

/*     AN EMPTY LIST IS CHARACTERIZED BY ICG(JV)=JV, IFG(JV)=JV. */
/*     OBVIOUSLY, IT IS QUITE EASY TO REMOVE OR ADD POINTS TO THE LIST. */
/*     SUCH RE-ARRANGEMENTS ARE NECESSARY AS (IN THE COLORING PART OF */
/*     THE COARSENING ALGORITHM) UNDECIDED POINTS BECOME C- OR F-POINTS */
/*     (THEY HAVE TO BE REMOVED FROM THEIR LIST) AND THE MEASURE VALUES */
/*     JV OF SOME POINTS CHANGE DURING THE ALGORITHM. SUCH POINTS HAVE */
/*     TO BE MOVED FROM ONE LIST TO ANOTHER ONE. IN ORDER TO KEEP TRACK */
/*     OF THOSE POINTS WHICH HAVE TO BE RE-ARRANGED, A RESET-STACK IS */
/*     USED WHICH CONTAINS ALL THESE POINTS I (I.E. THEIR "LOGICAL" */
/*     NUMBERS II). THE POINTER IR IS USED TO TOCHKA FROM ONE TOCHKA IN */
/*     THE STACK TO THE PREVIOUS ONE. */

/*     - IR(J) --------- J=1,...,NPTS */

/*     USED BELOW AS STACK-POINTER FOR POINTS TO BE RESET IN PERFORMING */
/*     THE COLORING ALGORITHM. ITOP IS THE TOP-OF-STACK POINTER. */
/*     THE STACK IS EMPTY IF ITOP=-1. */

/*     THE GLOBAL PICTURE OF THE POINTERS IS SKETCHED IN THE FOLLOWING */


/*                                          LIST END */
/*                                       <-------------| */
/*                                         LIST START  | */
/*                                       <------------|| */
/*                                                    || */
/*                  |-- II=I+NPTS+1 ->|            IFG||ICG */
/*                  |                 |               || */
/*                  |                 |               || */
/*       IMIN(K)          IMAX(K)           JVAL0             JVALX */
/*          |       I        |        II      |       JV        | */
/*          |=======*========|========*=======|========*========| */
/*          |                |                |                 | */
/*               GRID K          LOGICAL NUM      LIST ORIGINS */
/*                                              (MEASURE VALUES) */
/*                 ||                                  | */
/*                 ||         IFG                      | */
/*                 ||--------------------------------->| */
/*                 | */
/*                 |     |---> 1 (C) */
/*                 |     | */
/*                 | ICG |---> 0 (FF) */
/*                 |-----| */
/*                       |--->-1 (F) */
/*                       | */
/*                       |--->-2 (U) */


/*          1                NPTS */
/*          |        J        | */
/*          |========*========| */
/*          |                 | */
/*           SHIFTED GRID K */
/*                   | */
/*                   | IR */
/*                   |-----> RESET STACK (ITOP = TOP-OF-STACK) */



/* ===> PREPARATION */

    /* Parameter adjustments */
    --time;
    --jtr;
    --ir;
    --ifg;
    --icg;
    --imaxw;
    --imax;
    --imin;
    --iw;
    --ja;
    --ia;

    /* Function Body */
	told=clock();
    ilo = imin[*k];
    ihi = imax[*k];
    if (*k != 1) {
	iws = imaxw[*k - 1] + 2 - ilo;
    } else {
	iws = 0;
    }
    ilo1 = ilo - 1;
    npts = imax[*k] - imin[*k] + 1;
    npts1 = npts + 1;
	
	if (AMG1R6_LABEL == 1) 
	{

		/* old value of ntrlim has been removed */
		/*     NTRLIM = 2*(IW(IHI+IWS+1)-IW(ILO+IWS))/NPTS */
		/* new value of ntrlim: (modified by Krechel 22.07.02) */
		ntrlim = 0;
		i__1 = ihi;
		for (i__ = ilo; i__ <= i__1; ++i__) {
			/* Computing MAX */
			i__2 = ntrlim, i__3 = iw[i__ + iws + 1] - iw[i__ + iws];
			ntrlim = myi_max(i__2, i__3);
		}
	}
	else {
		ntrlim = ((iw[ihi + iws + 1] - iw[ilo + iws]) << 1) / npts; // скобки поставлены в соответствии с приоритетом.
	}


    nscnmx = 0;
    jval0 = ihi + npts1 + 1;
    jvalx = jval0 + (ntrlim << 1) + 1;
/* Computing MAX */
    i__1 = *mdicg, i__2 = jvalx + *iias;
    *mdicg = myi_max(i__1,i__2);
    *mdir = myi_max(*mdir,npts);
    if (jvalx > *ndicg) {
	goto L9901;
    }
    if (npts > *ndir) {
	goto L9902;
    }

/* ===> PUT INITIAL "MEASURE" FOR EACH TOCHKA I INTO IFG(I). */

    i__1 = ihi;
    for (i__ = ilo; i__ <= i__1; ++i__) {
		if (AMG1R6_LABEL == 1) {
			/*       NSCN = IW(I+IWS+1)-IW(I+IWS) */
			/*       IF (NSCN.LE.NTRLIM) THEN */
			/*         IFG(I) = JVAL0+NSCN */
			/*       ELSE */
			/*         IFG(I) = JVALX */
			/*       ENDIF */
			ifg[i__] = jval0 + iw[i__ + iws + 1] - iw[i__ + iws];
		}
		else {

			nscn = iw[i__ + iws + 1] - iw[i__ + iws];
			if (nscn <= ntrlim) {
				ifg[i__] = jval0 + nscn;
			}
			else {
				ifg[i__] = jvalx;
			}
		}
/* L2: */
    }

/* ===> SET CIRCULARLY LINKED LISTS AND RESET-STACK TO EMPTY */

    itop = -1;
    i__1 = jvalx;
    for (j = jval0; j <= i__1; ++j) {
	icg[j] = j;
	ifg[j] = j;
/* L10: */
    }

/* ===> PUT ALL U-POINTS OF GRID K INTO LISTS (I.E., NO FORCED F-POINTS). */
/*     ADD POINTS ALWAYS TO THE END OF THEIR CORRESPONDING LIST. */
/*     IN THE FOLLOWING, JCNBHI DENOTES THE ACTUAL HIGHEST MEASURE VALUE */
/*     (AMONG U-POINTS). */

    jcnbhi = 0;
    i__1 = ihi;
    for (i__ = ilo; i__ <= i__1; ++i__) {
	if (ja[ia[i__]] > ia[i__]) {
	    goto L15;
	}
	icg[i__] = 0;
	goto L20;
L15:
	icg[i__] = -2;
	ii = i__ + npts1;
	//ir[i__ - ilo1] = 0;
	ir[i__ - ilo1] = 0.0;
	jv = ifg[i__];
	if (jv > jcnbhi) {
	    jcnbhi = jv;
	}
	icg[ii] = icg[jv];
	ifg[ii] = jv;
	icg[jv] = ii;
	ifg[icg[ii]] = ii;
L20:
	;
    }

/*     **************** */
/*     * PRE-COLORING * */
/*     **************** */

/*     PICK A U-TOCHKA IC WITH MAXIMAL MEASURE AS GIVEN BY JCNBHI. */
/*     (TAKE THE FIRST FROM CORRESPONDING LIST: FIRST-IN/FIRST-OUT) */
/*     MAKE THAT TOCHKA A C-TOCHKA AND REMOVE IT FROM LISTS. THEN MAKE */
/*     ALL STRONG TRANSPOSE U-CONNECTIONS F-POINTS AND REMOVE THEM FROM */
/*     THEIR LISTS. UPDATE THE MEASURE OF IMPORTANCE FOR U-POINTS TO */
/*     BECOME C-POINTS AND ADD THESE POINTS TO THE RESET-STACK. FINALLY, */
/*     RE-ARRANGE THE LISTS BY SWEEPING THROUGH THE RESET STACK. THEN */
/*     PICK ANOTHER U-TOCHKA IC. */

/*     IF LIST CORRESPONDING TO JCNBHI IS EMPTY, GO TO NEXT LOWER VALUE. */
/*     PRE-COLOURING IS FINISHED IF JCNBHI<=JVAL0. ALL U-POINTS LEFT AT */
/*     THAT TIME, WILL BE REGARDED AS F-POINTS LATER. */

L30:
    if (jcnbhi <= jval0) {
	goto L100;
    }
    iic = ifg[jcnbhi];
    if (iic != jcnbhi) {
	goto L40;
    }
    --jcnbhi;
    goto L30;

/* ===> CREATE C-TOCHKA */

L40:
    ic = iic - npts1;
    icg[ic] = 1;
    icg[ifg[iic]] = icg[iic];
    ifg[icg[iic]] = ifg[iic];

/* ===> FOR TOCHKA IC WITH ECCESSIVE NUMBER OF STRONG TRANSPOSE */
/*     CONNECTIONS: MAKE IT A C-TOCHKA BUT LET ITS CONNECTED POINTS */
/*     REMAIN UNDECIDED. */

    if (jcnbhi == jvalx) {
	goto L78;
    }

/* ===> CREATE F-POINTS AROUND ABOVE C-TOCHKA */

    i__1 = iw[ic + iws + 1] - 1;
    for (j = iw[ic + iws]; j <= i__1; ++j) {
	//i__ = jtr[j];
    i__ = static_cast<integer>(jtr[j]);
	if (icg[i__] != -2) {
	    goto L77;
	}
	icg[i__] = -1;
	ii = i__ + npts1;
	icg[ifg[ii]] = icg[ii];
	ifg[icg[ii]] = ifg[ii];

/* ===>   INCREMENT MEASURE FOR ALL STRONG U-CONNECTIONS III OF I */
/*       (IF NOT YET MEASURE=JVALX) AND PUT THEM ON RESET STACK (IF NOT */
/*       YET THERE) */

	i__2 = ja[ia[i__]];
	for (jj = ia[i__] + 1; jj <= i__2; ++jj) {
	    iii = ja[jj];
	    if (icg[iii] != -2 || ifg[iii] >= jvalx) {
		goto L76;
	    }
	    ++ifg[iii];
	    iiii = iii - ilo1;
	   // if (ir[iiii] != 0) {
		//goto L76;
	    //}
		if ((static_cast<integer>(ir[iiii])) != 0) {
		goto L76;
	    }
	   // ir[iiii] = itop;
		 ir[iiii] = static_cast<doublereal>(itop);
	    itop = iiii;
L76:
	    ;
	}
L77:
	;
    }

/* ===> DECREMENT MEASURE FOR ALL STRONG U-CONNECTIONS I OF IC */
/*     AND PUT THEM ON RESET-STACK (IF NOT YET THERE) */

L78:
    i__1 = ja[ia[ic]];
    for (j = ia[ic] + 1; j <= i__1; ++j) {
	i__ = ja[j];
	if (icg[i__] != -2) {
	    goto L87;
	}
	--ifg[i__];
	ii = i__ - ilo1;
	//if (ir[ii] != 0) {
	  //  goto L87;
	//}
	if ((static_cast<integer>(ir[ii])) != 0) {
	    goto L87;
	}
	//ir[ii] = itop;
	ir[ii] = static_cast<doublereal>(itop);
	itop = ii;
L87:
	;
    }

/* ===> REARRANGE THE LISTS BY SWEEPING THROUGH RESET-STACK. */
/*     THEN GO BACK TO PICK ANOTHER U-TOCHKA IC. */

    in = itop;
    itop = -1;

L90:
    if (in <= 0) {
	goto L30;
    }
    i__ = in + ilo1;
    ii = i__ + npts1;
    if (icg[i__] != -2) {
	goto L95;
    }
    ifg[icg[ii]] = ifg[ii];
    icg[ifg[ii]] = icg[ii];
    jv = ifg[i__];
    if (jv > jcnbhi) {
	jcnbhi = jv;
    }
    icg[ii] = icg[jv];
    ifg[ii] = jv;
    icg[jv] = ii;
    ifg[icg[ii]] = ii;
L95:
    iip = in;
   // in = ir[in];
	in =static_cast<integer>(ir[in]);
    //ir[iip] = 0;
	ir[iip] = 0.0;
    goto L90;

L100:
    //ctime_(&tnew);
	tnew=clock();
    time[2] = time[2] + tnew - told;
    return 0;

/* ===> ERROR MESSAGES */

L9901:
    //io___157.ciunit = *ium;
    //s_wsfe(&io___157);
    //e_wsfe();
	printf("*** ERROR IN PCOL: NDIG TOO SMALL ***\n");
    *ierr = 6;
    return 0;

L9902:
    //io___158.ciunit = *ium;
    //s_wsfe(&io___158);
    //e_wsfe();
	printf("*** ERROR IN PCOL: NDU TOO SMALL ***\n");
    *ierr = 4;
    return 0;

} /* pcol_ */


/* ....................................................................... */

/*     PWINT                                             SUBROUTINE */

/* ....................................................................... */

/* Subroutine */ integer pwint_(integer *k, doublereal *ewt2, integer *ichk, 
	integer *mmax, doublereal *a, integer *ia, integer *ja, integer *iw, 
	integer *imin, integer *imax, integer *iminw, integer *imaxw, integer 
	*ifg, integer *icg, /*integer*/ doublereal *ncolor, /*real*/ unsigned int *time, integer *ierr, 
	integer *irow0, integer *ncolx, integer *nda, integer *ndja, integer *
	ndicg, integer *ndia, integer *ium, integer *mda, integer *mdja, 
	integer *iajas)
{
    /* Format strings */
   // static char fmt_9000[] = "(\002 INTERPOLATION OPERATOR NO.\002,i3,\002 C"
   //	    "OMPLETED. C-POINTS\002,\002 ADDED IN PWINT:\002,i4)";
   // static char fmt_9030[] = "(\002 *** ERROR IN PWINT: NDIA TOO SMALL **"
	//    "*\002)";
    //static char fmt_9040[] = "(\002 *** ERROR IN PWINT: NDIG TOO SMALL **"
	//    "*\002)";
    //static char fmt_9020[] = "(\002 *** ERROR IN PWINT: NDA TOO SMALL ***"
	 //   "\002)";
   // static char fmt_9050[] = "(\002 *** ERROR IN PWINT: NDJA TOO SMALL **"
	//    "*\002)";

    /* System generated locals */
    integer i__1=0, i__2=0, i__3=0;

    /* Builtin functions */
    //integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);

    /* Local variables */
    integer i__=0, j=0;
    doublereal s=0.0;
    integer j1=0, j2=0, ic=0, nc=0, ii=0, jj=0, ip=0;
    doublereal si=0.0;
    integer is=0;
    doublereal ww=0.0;
    integer jw0=0, ihi=0, jhi=0, ilo=0, jlo=0, jwx=0, icgp=0, jjhi=0, jjlo=0;
	unsigned int told=0, tnew=0;
    integer npts=0;
    doublereal ewt2i=0.0;
    integer ndaja=0, iblck=0;
    doublereal scale=0.0;
    integer nptsc=0, jwpos=0, iblck1=0, ncondc=0, ncount=0;

    /* Fortran I/O blocks */
   // static cilist io___192 = { 0, 0, 0, fmt_9000, 0 };
   // static cilist io___194 = { 0, 0, 0, fmt_9030, 0 };
    //static cilist io___195 = { 0, 0, 0, fmt_9040, 0 };
   // static cilist io___196 = { 0, 0, 0, fmt_9020, 0 };
    //static cilist io___197 = { 0, 0, 0, fmt_9050, 0 };



/*     SET UP FINAL COARSER GRID K+1 AND INTERPOLATION FORMULA FROM */
/*     GRID K+1 TO GRID K. THIS IS THE VERSION AS DESCRIBED IN RUGE/ */
/*     STUEBEN (BRISTOL). PWINT ASSUMES THE GRID TO BE PRE-COLORED */
/*     BY SUBROUTINE PCOL. */

/*     ON EXIT, A, JA AND IW ARE SET TO CONTAIN THE INTERPOLATION */
/*     WEIGHTS AND CORRESPONDING POINTERS AS REQUIRED IN THE SOLUTION */
/*     PHASE OF AMG1R5. ALSO, ICG(I) (IMIN(K)<=I<=IMAX(K)) ARE SET TO */
/*     THEIR FINAL VALUES, EXCEPT FOR THOSE I WITH ICG(I)<0. */

/*     ============ COMMENTS ON INPUT =================================== */

/*     - JA(IA(I)) ----- I=IMIN(K),...,IMAX(K) */

/*     ASSUMED TO TOCHKA TO THE LAST STRONG CONNECTION OF TOCHKA I (OR TO */
/*     IA(I) IF THERE IS NO SUCH CONNECTION). ON EXIT, RESET TO ORIGINAL */
/*     VALUES (I.E. JA(IA(I))=I). */

/*     - ICG(I) -------- I=IMIN(K),...,IMAX(K) */

/*     ASSUMED TO CONTAIN INFORMATION ON PRE-COLORING: */

/*       ICG(I)>0: I IS C-TOCHKA */
/*       ICG(I)=0: I IS FORCED F-TOCHKA (I.E. I HAS NO STRONG CONNNECTION) */
/*       ICG(I)<0: I IS F-TOCHKA WITH (AT LEAST) ONE STRONG CONNECTION. */

/*     ============== COMMENTS ON WORK SPACE USED ======================= */

/*     - IFG(I) -------- I=IMIN(K),...,IMAX(K) */

/*     IS USED FOR SEVERAL PURPOSES. IN PARTICULAR, TO DISTINGUISH */
/*     INTERPOLATORY AND NON-INTERPOLATORY POINTS. */

/*     - NCOLOR(I) ----- I=1,...,#POINTS ON GRID K */

/*     IS SET TO F-TOCHKA-COLORS TO BE USED LATER IN SUBROUTINE OPDFN. */


    /* Parameter adjustments */
    --time;
    --ncolor;
    --icg;
    --ifg;
    --imaxw;
    --iminw;
    --imax;
    --imin;
    --iw;
    --ja;
    --ia;
    --a;

    /* Function Body */
	told=clock();
    ndaja = my_imin(*nda,*ndja);
    ncount = 0;
    if (*k == 1) {
	iminw[1] = 1;
	iw[1] = ia[imax[1] + 1];
    }

    ilo = imin[*k];
    ihi = imax[*k];
    i__1 = ihi;
    for (i__ = ilo; i__ <= i__1; ++i__) {
	ifg[i__] = 0;
/* L12: */
    }

/* ===> SWEEP OVER F-POINTS I WHICH HAVE AT LEAST ONE STRONG CONNECTION */

    iblck = iminw[*k];
    i__1 = ihi;
    for (i__ = ilo; i__ <= i__1; ++i__) {
	if (icg[i__] >= 0) {
	    goto L400;
	}
	jlo = ia[i__] + 1;
	jhi = ja[ia[i__]];
	ewt2i = *ewt2 / a[jlo];
	ncondc = 0;

/* ===>   INITIALIZE "BLOCK" OF INTERPOLATION WEIGHTS OF TOCHKA I */

L30:
	jw0 = iw[iblck];
	jwx = jw0;
	if (jwx > ndaja) {
	    goto L2000;
	}
	i__2 = jhi;
	for (j = jlo; j <= i__2; ++j) {
	    ii = ja[j];
	    if (icg[ii] <= 0) {
		goto L20;
	    }
	    a[jwx] = a[j];
	    ja[jwx] = ii;
	    ifg[ii] = jwx;
	    ++jwx;
	    if (jwx > ndaja) {
		goto L2000;
	    }
L20:
	    ;
	}
	a[jwx] = a[ia[i__]];
	ja[jwx] = i__;
/* Computing MAX */
	i__2 = *mda, i__3 = jwx + *iajas;
	*mda = myi_max(i__2,i__3);
/* Computing MAX */
	i__2 = *mdja, i__3 = jwx + *iajas;
	*mdja = myi_max(i__2,i__3);
	ifg[i__] = jwx;

/* ===>   SWEEP OVER STRONGLY CONNECTED F-POINTS. THESE MUST BE "COVERED" */
/*       BY A TOTAL WEIGHT DEFINED BY EWT2. IF AN F-TOCHKA HAS NO STRONG */
/*       CONNECTIONS, REGARD IT TO BE COVERED, BUT DO NOT DISTRIBUTE */
/*       THE CORRESPONDING WEIGHT (ERROR AT SUCH A TOCHKA CAN BE ASSUMED */
/*       TO BE VERY SMALL!). */

	i__2 = jhi;
	for (j = jlo; j <= i__2; ++j) {
	    ii = ja[j];
	    if (icg[ii] >= 0) {
		goto L150;
	    }

/* ===>     COMPUTE DEPENDENCE ON SET OF INTERPOLATION POINTS */

	    s = 0.;
	    si = 0.;
	    jjlo = ia[ii] + 1;
	    jjhi = ia[ii + 1] - 1;
	    i__3 = jjhi;
	    for (jj = jjlo; jj <= i__3; ++jj) {
		if (ifg[ja[jj]] < jw0) {
		    goto L110;
		}
		if (ja[jj] == i__) {
		    si = a[jj];
		}
		s += a[jj];
L110:
		;
	    }
	    if (*ichk == 2) {
		goto L111;
	    }
	    if (s == 0.) {
		a[jwx] += a[j];
		goto L150;
	    } else {
		goto L135;
	    }

/* ===>     CHECK DEPENDENCE ON SET OF INTERPOLATION POINTS */

L111:
	    if (s - si <= ewt2i * a[j] * a[jjlo]) {
		goto L135;
	    }

/* ===>     DEPENDENCE TOO SMALL: IF THERE IS NOT YET A CONDITIONAL */
/*         C-TOCHKA, MAKE II SUCH A TOCHKA AND RESTART THE */
/*         PROCESS FOR DEFINING INTERPOLATION WEIGHTS FOR TOCHKA I. */
/*         OTHERWISE MAKE I ITSELF A C-TOCHKA AND LEAVE II AN F-TOCHKA. */

	    if (ncondc == 0) {
		++ncount;
		ncondc = 1;
		ip = ii;
		icgp = icg[ii];
		icg[ii] = 1;
		goto L30;
	    } else {
		icg[i__] = 1;
		icg[ip] = icgp;
		i__3 = jwx;
		for (jj = jw0; jj <= i__3; ++jj) {
		    ifg[ja[jj]] = 0;
/* L120: */
		}
		goto L400;
	    }

/* ===>     DISTRIBUTE THE WEIGHT OF TOCHKA II */

L135:
	    ww = a[j] / s;
	    i__3 = jjhi;
	    for (jj = jjlo; jj <= i__3; ++jj) {
		if (ifg[ja[jj]] >= jw0) {
		    a[ifg[ja[jj]]] += a[jj] * ww;
		}
/* L140: */
	    }
L150:
	    ;
	}

/* ===>   ALL NECESSARY POINTS ARE COVERED. NOW DISTRIBUTE WEIGHTS FROM */
/*       WEAK CONNECTIONS OF TOCHKA I (ANALOGOUS AS ABOVE) */

	i__2 = ia[i__ + 1] - 1;


	// 3,03,2019
// Нельзя параллелить портит сходимость
//#pragma omp parallel for
	for (integer j_loc = jhi + 1; j_loc <= i__2; ++j_loc) 
	{
	    integer ii_loc = ja[j_loc];
		if (icg[ii_loc] != 0) {
			doublereal s_loc = 0.;
			integer jjlo_loc = ia[ii_loc] + 1;
			integer jjhi_loc = ia[ii_loc + 1] - 1;
			integer i__3_loc = jjhi_loc;
// Нельзя параллелить замедляет время приложения в двое.
//#pragma omp parallel  for shared(jjlo_loc, i__3_loc, jw0, ifg, ja, a) reduction(+:s_loc)
			for (integer jj_loc = jjlo_loc; jj_loc <= i__3_loc; ++jj_loc) 
			{
				if (ifg[ja[jj_loc]] >= jw0) {
					s_loc += a[jj_loc];
				}
				/* L160: */
			}

			
			if (s_loc == 0.) {
				a[jwx] += a[j_loc];
			}
			else {
				doublereal ww_loc = a[j_loc] / s_loc;
				integer i__3_loc2 = jjhi_loc;
				for (integer jj_loc = jjlo_loc; jj_loc <= i__3_loc2; ++jj_loc) {
					if (ifg[ja[jj_loc]] >= jw0) {
						a[ifg[ja[jj_loc]]] += a[jj_loc] * ww_loc;
					}
					/* L170: */
				}
			}
			
		}
	}

	icg[i__] = -iblck;
	++iblck;
	iw[iblck] = jwx + 1;
L400:
	;
    }

/* ===> SET ICG; RESET JA(IA(I)); CHECK SIZE OF COARSEST GRID */

    ic = ihi;
    imin[*k + 1] = ic + 1;
    i__1 = ihi;
    for (i__ = ilo; i__ <= i__1; ++i__) {
	ja[ia[i__]] = i__;
	if (icg[i__] <= 0) {
	    goto L900;
	}
	++ic;
	icg[i__] = ic;
	if (ic < *ndicg) {
	    icg[ic] = 0;
	}
L900:
	;
    }
    imax[*k + 1] = ic;

    npts = ihi - ilo + 1;
    nptsc = imax[*k + 1] - imin[*k + 1] + 1;
    if (nptsc == 1) {
	*mmax = *k + 1;
    }
    if (nptsc == 1 && *irow0 < 2) {
	*mmax = *k;
    }
    if (nptsc == npts || nptsc == 0) {
	*mmax = *k;
    }
    if (*k >= *mmax) {
	goto L1000;
    }
    if (ic >= *ndia) {
	goto L1700;
    }
    if (ic >= *ndicg) {
	goto L1800;
    }

/* ===> RE-ARRANGE A */

L1000:
    iblck1 = iminw[*k];
    jwpos = iw[iblck1];
    i__1 = ihi;
    for (i__ = ilo; i__ <= i__1; ++i__) {
	if (icg[i__] >= 0) {
	    goto L950;
	}
	iblck = -icg[i__];
	j1 = iw[iblck];
	j2 = iw[iblck + 1] - 1;
	if (j2 <= j1) {
	    icg[i__] = 0;
	} else {
	    icg[i__] = -iblck1;
	    iw[iblck1] = jwpos;
	    scale = -1. / a[j2];
	    i__2 = j2 - 1;
	    for (j = j1; j <= i__2; ++j) {
		a[jwpos] = a[j] * scale;
		ja[jwpos] = icg[ja[j]];
		++jwpos;
/* L920: */
	    }
	    ++iblck1;
	}
L950:
	;
    }
    imaxw[*k] = iblck1 - 1;
    iw[iblck1] = jwpos;

/* ===> STORE TYPE OF POINTS (I.E. C, F OR FF)  ON VECTOR NCOLOR */

    is = 1 - ilo;
    i__1 = ihi;
    for (i__ = ilo; i__ <= i__1; ++i__) {
	nc = icg[i__];
	if (nc < 0) {
	   // ncolor[i__ + is] = 1;
		ncolor[i__ + is]=1.0;
	} else if (nc > 0) {
	   // ncolor[i__ + is] = 2;
		ncolor[i__ + is] = 2.0;
	} else {
	   // ncolor[i__ + is] = 3;
		ncolor[i__ + is] = 3.0;
	}
/* L960: */
    }
    *ncolx = 1;

/* ===> EXIT */

    //io___192.ciunit = *ium;
    //s_wsfe(&io___192);
    //do_fio(&c__1, (char *)&(*k), (ftnlen)sizeof(integer));
    //do_fio(&c__1, (char *)&ncount, (ftnlen)sizeof(integer));
    //e_wsfe();
	if (yes_print_amg) {
#if doubleintprecision == 1
		printf("( INTERPOLATION OPERATOR NO.,%lld, \n", k[0]);
		printf("COMPLETED. C-POINTS, ADDED IN PWINT:,%lld)\n", ncount);
#else
		printf("( INTERPOLATION OPERATOR NO.,%d, \n", k[0]);
		printf("COMPLETED. C-POINTS, ADDED IN PWINT:,%d)\n", ncount);
#endif
	 
	}
    //ctime_(&tnew);
	tnew=clock();
    time[4] = time[4] + tnew - told;
    return 0;

/* ===> ERROR MESSAGES */

L1700:
    //io___194.ciunit = *ium;
    //s_wsfe(&io___194);
    //e_wsfe();
	printf("( *** ERROR IN PWINT: NDIA TOO SMALL ***)\n");
    *ierr = 2;
    return 0;

L1800:
    //io___195.ciunit = *ium;
    //s_wsfe(&io___195);
    //e_wsfe();
	printf( "( *** ERROR IN PWINT: NDIG TOO SMALL ***)\n");
    *ierr = 6;
    return 0;

L2000:
    if (*nda <= *ndja) {
	//io___196.ciunit = *ium;
	//s_wsfe(&io___196);
	//e_wsfe();
	printf("( *** ERROR IN PWINT: NDA TOO SMALL ***)\n");
	*ierr = 1;
    } else {
	//io___197.ciunit = *ium;
	//s_wsfe(&io___197);
	//e_wsfe();
	printf("( *** ERROR IN PWINT: NDJA TOO SMALL ***)\n");
	*ierr = 3;
    }
    return 0;

} /* pwint_ */


/* ....................................................................... */

/*     OPDFN                                             SUBROUTINE */

/* ....................................................................... */

/* Subroutine */ integer opdfn_(integer *k, integer *ierr, integer *mmax, 
	doublereal *a, integer *ia, integer *ja, integer *iw, integer *imin, 
	integer *imax, integer *iminw, integer *imaxw, integer *icg, integer *
	ifg, /*integer*/ doublereal *ncolor, integer *nstcol, integer *ncolx, /*real*/ unsigned int *time, 
	integer *nda, integer *ndja, integer *ium, integer *mda, integer *
	mdja, integer *iajas)
{
    /* Format strings */
    //static char fmt_9940[] = "(\002 *** ERROR IN OPDFN: INTERPOLATION ENTRY "
	//    "MISSING ON GRID\002,i3)";
   // static char fmt_9930[] = "(\002 --- WARNG: UNABLE TO STORE TRANSPOSE OF "
	//    "INTERPOLATION\002,\002 ON GRID \002,i2,\002 DURING \002/\002    "
	//    "        EXECUTION OF OPDFN, BECAUSE NDA OR NDJA \002,\002TOO SMA"
	 //   "LL.\002/\002            SETUP COMPUTATION IS SLOWING DOWN.\002)";
    //static char fmt_9000[] = "(\002 COARSE  GRID  OPERATOR NO.\002,i3,\002 C"
	//    "OMPLETED\002)";
    //static char fmt_9910[] = "(\002 *** ERROR IN OPDFN: NDA TOO SMALL ***"
	  //  "\002)";
   // static char fmt_9920[] = "(\002 *** ERROR IN OPDFN: NDJA TOO SMALL **"
	//    "*\002)";

    /* System generated locals */
    integer i__1=0, i__2=0, i__3=0, i__4=0;

    /* Builtin functions */
    //integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);

    /* Local variables */
    integer i__=0, j=0, k1=0, ic=0, jb=0, if__=0;
    doublereal ww=0.0;
    integer ic1=0, ic2=0, ic3=0, if1=0, jf1=0, jf2=0, jc3=0, if2=0, jf3=0, ibl=0, ihi=0, ilo=0,
	     ist=0, ilo1=0;
    doublereal wjf1=0.0;
    integer icol=0;
	unsigned int told=0;
    integer jpos=0;
	unsigned int tnew=0;
    integer npts=0, ndaja=0;
    integer iadrs=0, ialow=0, istti=0;

    /* Fortran I/O blocks */
   // static cilist io___216 = { 0, 0, 0, fmt_9940, 0 };
   // static cilist io___223 = { 0, 0, 0, fmt_9930, 0 };
    //static cilist io___224 = { 0, 0, 0, fmt_9000, 0 };
    //static cilist io___232 = { 0, 0, 0, fmt_9910, 0 };
    //static cilist io___233 = { 0, 0, 0, fmt_9920, 0 };



/*     THIS SUBROUTINE CONSTRUCTS THE CG-OPERATOR A(K) (K>1). */

/*     ONE ROW IS CONSTRUCTED AT A TIME, AND THE ROW STRUCTURE */
/*     (I.E., WHICH CONNECTION GOES WHERE) MUST BE DETERMINED */
/*     AT THE SAME TIME.  IN ORDER TO AVOID SEARCHES THROUGH */
/*     THE CURRENT ROW TO DETERMINE IF A POSITION FOR A CONNECTION */
/*     HAS ALREADY BEEN DEFINED, ICG (FOR LEVEL K) AND IFG (FOR LEVEL */
/*     K-1) ARE USED AUXILIARILY. FURTHERMORE, TO SPEED UP COMPUTATION, */
/*     THE TRANSPOSE OF INTERPOLATION IS TEMPORARILY STORED ON A (AT */
/*     THE END OF THE AVAILABLE SPACE). IFG (FOR LEVEL K) IS USED */
/*     AS POINTER TO THE CORRESPONDING ROWS. MORE DETAILS: SEE BELOW. */

/*     ============ COMMENTS ON INPUT =================================== */

/*     - IMIN(K), IMAX(K) */

/*     FIRST/LAST NUMBER OF GRID POINTS ON GRID K. */

/*     - IMINW(K-1), IMAXW(K-1) */

/*     FIRST/LAST "BLOCK" NUMBER OF INTERPOLATION WEIGHTS CONTRIBUTING */
/*     IN INTERPOLATION TO GRID K-1. */

/*     - IW(IBL) ----- IBL=IMINW(K-1),IMAXW(K-1)+1 */
/*     - A(J) -------- J=JLO,JHI      (JLO=IW(IMINW(K-1))) */
/*     - JA(J) ------- J=JLO,JHI      (JHI=IW(IMAXW(K-1)+1)-1) */

/*     ARE ASSUMED TO CONTAIN THE WEIGHTS OF INTERPOLATION ALONG WITH */
/*     THE NECESSARY POINTERS. IW(IBL) POINTS TO THE FIRST ENTRY OF */
/*     "BLOCK" NUMBER IBL IN A (CF. BELOW). */

/*     - ICG(IF) ----- IF=IMIN(K-1),IMAX(K-1) */

/*     ARE ASSUMED TO BE SET TO THEIR FINAL VALUES: */

/*       ICG(IF)>0: IF IS C-TOCHKA OF GRID K-1 AND ICG(IF) POINTS JUST */
/*                  TO THE CORRESPONDING TOCHKA ON GRID K; */
/*       ICG(IF)=0: IF IS F-TOCHKA OF GRID K-1 WITHOUT ANY CONTRIBUTION */
/*                  IN INTERPOLATION FROM GRID K; */
/*       ICG(IF)<0: IF IS F-TOCHKA OF GRID K-1; IBL=-ICG(IF) POINTS JUST */
/*                  TO THE "BLOCK" OF INTERPOLATION WEIGHTS. THAT MEANS, */
/*                  JW(J) (IW(IBL)<=J<=IW(IBL+1)-1) POINTS TO THE POINTS */
/*                  (ON GRID K) WHICH CONTRIBUTE IN INTERPOLATION TO IF. */
/*                  THE CORRESPONDING WEIGHTS ARE STORED IN A(J). */

/*     - IFG(IF) ----- IF=IMIN(K-1),IMAX(K-1) */

/*     WORK SPACE (SEE BELOW). IFG(IF) IS ASSUMED TO BE > -IMIN(K). */

/*     - IFG(IC) ----- IC=IMIN(K),IMAX(K)+1 */

/*     WORK SPACE (SEE BELOW). THE CONTENTS OF IFG(IC) IS ARBITRARY. */

/*     - ICG(IC) ----- IC=IMIN(K),IMAX(K) */

/*     WORK SPACE (SEE BELOW). ICG(IC) IS ASSUMED TO BE ZERO. */

/*     ============== COMMENTS ON WORK SPACE USED ======================= */

/*     - A(J) -------- J=JJLO,JJHI */
/*     - JA(J) ------- J=JJLO,JJHI */

/*     (JJLO=MIN(NDA,NDJA)-IW(IMAXW(K-1)+1)+IW(IMINW(K-1)+1), */
/*      JJHI=MIN(NDA,NDJA), I.E. THE SPACE OF LENGTH OF INTERPOLATION */
/*      (K-1) AT THE END OF A AND JA, RESPECTIVELY.) */

/*     IS USED TO STORE THE TRANSPOSE OF INTERPOLATION W(K-1). THIS IS */
/*     DONE IN ORDER TO SPEED UP THE COMPUTATION OF THE OPERATOR A(K). */
/*     IF DURING ASSEMBLAGE OF A(K) THIS WORK SPACE IS REQUIRED BY A(K), */
/*     PROCESSING CONTINUES USING THE NOT YET REDEFINED PART OF THE WORK */
/*     SPACE AND RECALCULATING THE LOST ENTRIES EACH TIME THEY ARE */
/*     NEEDED. THUS THE CALCULATION SLOWS DOWN. */

/*     - IFG(IC) ----- IC=IMIN(K),IMAX(K)+1 */

/*     IS USED AS POINTER FOR THE TRANSPOSE OF INTERPOLATION: THE COARSE- */
/*     GRID TOCHKA IC CONTRIBUTES IN INTERPOLATION TO THE FINE-GRID POINTS */
/*     IF=JA(J) (IFG(IC)<=J<=IFG(IC+1)-1) WITH WEIGHT A(J). THE */
/*     CONTRIBUTION TO ITSELF (BY THE WEIGHT 1.0) IS NOT CONTAINED. */

/*     - ICG(IC) ----- IC=IMIN(K),IMAX(K) */

/*     IS USED FOR SEVERAL PURPOSES. IN ASSEMBLING THE CG-OPERATOR A(K), */
/*     IT SERVES AS POINTER TO POSITIONS IN A(K) WHICH HAVE ALREADY BEEN */
/*     DEFINED: IF THE CURRENT ROW CORRESPONDS TO TOCHKA IC1, AND A */
/*     CONNECTION TO ANOTHER CG-TOCHKA IC2 HAS JUST BEEN FOUND, THEN IF */
/*     ICG(IC2)<IA(IC1), THE CORRESPONDING ENTRY IN ROW IC1 OF A(K) */
/*     HAS NOT YET BEEN DEFINED. OTHERWISE, ICG(IC2) POINTS TO THE */
/*     LOCATION FOR THAT ENTRY. (ALSO SEE IFG BELOW.) */

/*     - IFG(IF) ----- IF=IMIN(K-1),IMAX(K-1) */

/*     IN ASSEMBLING THE CG-OPERATOR A(K), THIS VECTOR CONTAINS INFORMAT- */
/*     ION ON WHETHER THE EXISTENCE OF ENTRIES IN A(K) HAS TO BE CHECKED */
/*     OR NOT. */


/* ===> EXTEND A, JA TO STORE TRANSPOSE OF INTERPOLATION */

    /* Parameter adjustments */
    --time;
    --nstcol;
    --ncolor;
    --ifg;
    --icg;
    --imaxw;
    --iminw;
    --imax;
    --imin;
    --iw;
    --ja;
    --ia;
    --a;

    /* Function Body */
	told=clock();
    ndaja = my_imin(*nda,*ndja);
    iminw[*k] = imaxw[*k - 1] + 1;
    if (*k > *mmax) {
	goto L8000;
    }
    jpos = iw[iminw[*k]];
    i__1 = iw[imaxw[*k - 1] + 1] - 1;
    for (j = iw[iminw[*k - 1]]; j <= i__1; ++j) {
	++icg[ja[j]];
/* L200: */
    }

/* Computing MAX */
    i__1 = *mda, i__2 = jpos + jpos - iw[iminw[*k - 1]] + *iajas;
    *mda = myi_max(i__1,i__2);
/* Computing MAX */
    i__1 = *mdja, i__2 = jpos + jpos - iw[iminw[*k - 1]] + *iajas;
    *mdja = myi_max(i__1,i__2);
    ifg[imin[*k]] = ndaja - jpos + iw[iminw[*k - 1]] + 1;
    if (ifg[imin[*k]] <= jpos) {
	goto L9900;
    }
    i__1 = imax[*k];
    for (ic = imin[*k]; ic <= i__1; ++ic) {
	ifg[ic + 1] = ifg[ic] + icg[ic];
	icg[ic] = ifg[ic];
/* L220: */
    }

    i__1 = imax[*k - 1];
    for (if__ = imin[*k - 1]; if__ <= i__1; ++if__) {
	if (icg[if__] >= 0) {
	    goto L250;
	}
	ibl = -icg[if__];
	i__2 = iw[ibl + 1] - 1;
	for (j = iw[ibl]; j <= i__2; ++j) {
	    ic = ja[j];
	    a[icg[ic]] = a[j];
	    ja[icg[ic]] = if__;
	    ++icg[ic];
/* L260: */
	}
L250:
	;
    }

    i__1 = imax[*k];
    for (ic = imin[*k]; ic <= i__1; ++ic) {
	icg[ic] = 0;
/* L300: */
    }

/* ===> SWEEP OVER ALL CG-POINTS IC TO ASSEMBLE ROWS OF CG MATRIX */

    istti = 1;
    i__1 = imax[*k - 1];
    for (if__ = imin[*k - 1]; if__ <= i__1; ++if__) {
	if (icg[if__] <= 0) {
	    goto L100;
	}
	ic = icg[if__];
	if (jpos > ndaja) {
	    goto L9900;
	}
	ialow = jpos;
	icg[ic] = jpos;
	a[jpos] = a[ia[if__]];
	ja[jpos] = ic;
	++jpos;

/*       ---------------------------------------------- */
/*       | SEARCH FOR C-C-C-C AND C-C-F-C CONNECTIONS | */
/*       ---------------------------------------------- */

	i__2 = ia[if__ + 1] - 1;
	for (jf1 = ia[if__] + 1; jf1 <= i__2; ++jf1) {
	    if1 = ja[jf1];
	    ic1 = icg[if1];
	    if (ic1 < 0) {
		goto L11;
	    } else if (ic1 == 0) {
		goto L25;
	    } else {
		goto L20;
	    }

/* ===>     IF1 IS F-TOCHKA: SWEEP OVER C-C-F-C CONNECTIONS */

L11:
	    ifg[if1] = -ic;
	    i__3 = iw[-ic1 + 1] - 1;
	    for (jf2 = iw[-ic1]; jf2 <= i__3; ++jf2) {
		ic2 = ja[jf2];
		if (icg[ic2] >= ialow) {
		    goto L10;
		}
		if (jpos > ndaja) {
		    goto L9900;
		}
		icg[ic2] = jpos;
		a[jpos] = a[jf1] * a[jf2];
		ja[jpos] = ic2;
		++jpos;
		goto L15;
L10:
		a[icg[ic2]] += a[jf1] * a[jf2];
L15:
		;
	    }
	    goto L25;

/* ===>     IF1 IS C-TOCHKA: C-C-C-C CONNECTION */

L20:
	    if (icg[ic1] >= ialow) {
		goto L23;
	    }
	    if (jpos > ndaja) {
		goto L9900;
	    }
	    icg[ic1] = jpos;
	    a[jpos] = a[jf1];
	    ja[jpos] = ic1;
	    ++jpos;
	    goto L25;
L23:
	    a[icg[ic1]] += a[jf1];
L25:
	    ;
	}
	ist = imin[*k - 1] - 1;

/*        ---------------------------------------------- */
/*        | SEARCH FOR C-F-C-C AND C-F-F-C CONNECTIONS | */
/*        ---------------------------------------------- */

	i__2 = ifg[ic + 1] - 1;
	for (jf1 = ifg[ic]; jf1 <= i__2; ++jf1) {
	    if (jf1 >= jpos) {
		if1 = ja[jf1];
		wjf1 = a[jf1];
		ist = if1;
	    } else {
		istti = 0;
		i__3 = imax[*k - 1];
		for (jf2 = ist + 1; jf2 <= i__3; ++jf2) {
		    if (icg[jf2] >= 0) {
			goto L120;
		    }
		    ibl = -icg[jf2];
		    i__4 = iw[ibl + 1] - 1;
		    for (jb = iw[ibl]; jb <= i__4; ++jb) {
			jc3 = ja[jb];
			if (jc3 != ic) {
			    goto L110;
			}
			if1 = jf2;
			wjf1 = a[jb];
			goto L130;
L110:
			;
		    }
L120:
		    ;
		}

/* ===>       ERROR EXIT */

		//io___216.ciunit = *ium;
		//s_wsfe(&io___216);
		i__3 = *k - 1;
		//do_fio(&c__1, (char *)&i__3, (ftnlen)sizeof(integer));
		//e_wsfe();
#if doubleintprecision == 1
		printf("( *** ERROR IN OPDFN: INTERPOLATION ENTRY MISSING ON GRID,%lld)\n", i__3);
#else
		printf("( *** ERROR IN OPDFN: INTERPOLATION ENTRY MISSING ON GRID,%d)\n", i__3);
#endif
		
		*ierr = 22;
		return 0;
L130:
		ist = if1;
	    }
	    i__3 = ia[if1 + 1] - 1;
	    for (jf2 = ia[if1]; jf2 <= i__3; ++jf2) {
		if2 = ja[jf2];
		ic2 = icg[if2];
		ww = wjf1 * a[jf2];
		if (ic2 < 0) {
		    goto L35;
		} else if (ic2 == 0) {
		    goto L90;
		} else {
		    goto L50;
		}

/* ===>       IF2 IS F-TOCHKA: SWEEP OVER C-F-F-C CONNECTIONS */

L35:
		if (ifg[if2] == -ic) {
		    goto L70;
		}
		ifg[if2] = -ic;
		i__4 = iw[-ic2 + 1] - 1;
		for (jf3 = iw[-ic2]; jf3 <= i__4; ++jf3) {
		    ic3 = ja[jf3];
		    if (icg[ic3] >= ialow) {
			goto L30;
		    }
		    if (jpos > ndaja) {
			goto L9900;
		    }
		    icg[ic3] = jpos;
		    a[jpos] = ww * a[jf3];
		    ja[jpos] = ic3;
		    ++jpos;
		    goto L40;
L30:
		    a[icg[ic3]] += ww * a[jf3];
L40:
		    ;
		}
		goto L90;

/* ===>       IF2 HAS BEEN ENCOUNTERED BEFORE; DO NOT CHECK POSITIONS! */

L70:
		i__4 = iw[-ic2 + 1] - 1;
		for (jf3 = iw[-ic2]; jf3 <= i__4; ++jf3) {
		    iadrs = icg[ja[jf3]];
		    a[iadrs] += ww * a[jf3];
/* L80: */
		}
		goto L90;

/* ===>       IF2 IS C-TOCHKA: C-F-C-C CONNECTION */

L50:
		if (icg[ic2] >= ialow) {
		    goto L60;
		}
		if (jpos > ndaja) {
		    goto L9900;
		}
		icg[ic2] = jpos;
		a[jpos] = ww;
		ja[jpos] = ic2;
		++jpos;
		goto L90;
L60:
		a[icg[ic2]] += ww;
L90:
		;
	    }
/* L95: */
	}
	ia[ic + 1] = jpos;
L100:
	;
    }
/* Computing MAX */
    i__1 = *mda, i__2 = ia[imax[*k] + 1] - 1 + *iajas;
    *mda = myi_max(i__1,i__2);
/* Computing MAX */
    i__1 = *mdja, i__2 = ia[imax[*k] + 1] - 1 + *iajas;
    *mdja = myi_max(i__1,i__2);

/* ===> WARNING FOR STORAGE SHORTAGE */

    if (istti != 1) {
	k1 = *k - 1;
	//io___223.ciunit = *ium;
	//s_wsfe(&io___223);
	//do_fio(&c__1, (char *)&k1, (ftnlen)sizeof(integer));
	//e_wsfe();

	    printf("( --- WARNG: UNABLE TO STORE TRANSPOSE OF \n");
#if doubleintprecision == 1
		printf("INTERPOLATION ON GRID ,%lld, DURING     \n", k1);
#else
		printf("INTERPOLATION ON GRID ,%d, DURING     \n", k1);
#endif
	    
	    printf("  EXECUTION OF OPDFN, BECAUSE NDA OR NDJA TOO SMALL.\n");
	    printf("  SETUP COMPUTATION IS SLOWING DOWN.)\n");

	*ierr = -1;
    }
    //io___224.ciunit = *ium;
    //s_wsfe(&io___224);
    //do_fio(&c__1, (char *)&(*k), (ftnlen)sizeof(integer));
    //e_wsfe();
	if (yes_print_amg) {
#if doubleintprecision == 1
		printf("( COARSE  GRID  OPERATOR NO.,%lld, COMPLETED)\n", k[0]);
#else
		printf("( COARSE  GRID  OPERATOR NO.,%d, COMPLETED)\n", k[0]);
#endif
	   
	}

/* ===> SET UP LINKED LIST FOR RELAXATION ON GRID K-1 */

L8000:
    ia[imin[*k]] = iw[iminw[*k]];
    ilo = imin[*k - 1];
    ihi = imax[*k - 1];
    ilo1 = ilo - 1;
    npts = ihi - ilo1;
#if doubleintprecision == 1
	// 1000000M
	ist = 1000000000000;
#else
	// 100M
	ist = 100000000;
#endif
    
    i__1 = *ncolx;
    for (icol = 1; icol <= i__1; ++icol) {
	for (i__ = npts; i__ >= 1; --i__) {
	    if ((static_cast<integer>(ncolor[i__])) != icol) {
		goto L8040;
	    }
	    icg[i__ + ilo1] = -ist;
	    ist = i__ + ilo1;
L8040:
	    ;
	}
/* L8100: */
    }
    nstcol[*k - 1] = ist;

/* ===> EXIT / ERROR MESSAGES */

    //ctime_(&tnew);
	tnew=clock();
    time[6] = time[6] + tnew - told;
    return 0;

L9900:
    if (*nda <= *ndja) {
	//io___232.ciunit = *ium;
	//s_wsfe(&io___232);
	//e_wsfe();
	printf( "( *** ERROR IN OPDFN: NDA TOO SMALL ***)\n");
	*ierr = 1;
    } else {
	//io___233.ciunit = *ium;
	//s_wsfe(&io___233);
	//e_wsfe();
	printf("( *** ERROR IN OPDFN: NDJA TOO SMALL ***)\n");
	*ierr = 3;
    }
    return 0;

} /* opdfn_ */


/* ....................................................................... */

/*     SETIFG                                              SUBROUTINE */

/* ....................................................................... */

/* Subroutine */ integer setifg_(integer *imin, integer *imax, integer *icg, 
	integer *ifg, integer *nstcol, integer *levels, /*real*/ unsigned int *time)
{
    /* System generated locals */
    integer i__1=0;

    /* Local variables */
    integer i__=0, k=0, ib=0, ist=0;
	unsigned int told=0, tnew=0;


/*     SET "INVERSE" POINTER IFG */


    /* Parameter adjustments */
    --time;
    --nstcol;
    --ifg;
    --icg;
    --imax;
    --imin;

    /* Function Body */
	told=clock();
    i__1 = imax[*levels - 1];
    for (i__ = imin[1]; i__ <= i__1; ++i__) {
	if (icg[i__] > 0) {
	    ifg[icg[i__]] = i__;
	}
/* L10: */
    }
    ib = 1;
    i__1 = *levels - 1;
    for (k = 1; k <= i__1; ++k) {
	ist = nstcol[k];
L20:
#if doubleintprecision == 1
	// 1000000M
	if (ist >= 1000000000000) {
		goto L30;
	}
#else
	// 100M
	if (ist >= 100000000) {
	    goto L30;
	}
#endif
	ifg[ib] = ist;
	ist = -icg[ist];
	++ib;
	goto L20;
L30:
	;
    }
	tnew=clock();
    time[8] = time[8] + tnew - told;
    return 0;
} /* setifg_ */

/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */

/*     AMG1R5 SOLUTION-SUBROUTINES */

/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
/* ....................................................................... */

/*     SOLVE                                                SUBROUTINE */

/* ....................................................................... */

/* Subroutine */ integer solve_(integer *madapt, integer *ncyc, integer *nrd, 
	integer *nsolco, integer *nru, integer *iout, integer *ierr, 
	doublereal *a, doublereal *u, doublereal *f, integer *ia, integer *ja,
	 integer *iw, doublereal *eps, integer *imin, integer *imax, integer *
	iminw, integer *imaxw, integer *icg, integer *ifg, integer *nstcol, 
	integer *iarr, /*real*/unsigned int *time, integer *ncyc0, integer *irow0, integer *
	levels, integer *nda, integer *ndja, integer *ndu, integer *ndf, 
	integer *mda, integer *mdja, integer *mdu, integer *mdf, integer *iup,
	 integer *ium, doublereal *resi, doublereal *res0, doublereal *res)
{
    /* Format strings */
   // static char fmt_9050[] = "(\002 *** ERROR IN SOLVE: NDU TOO SMALL ***"
	//    "\002)";
    //static char fmt_9060[] = "(\002 *** ERROR IN SOLVE: NDF TOO SMALL ***"
	//    "\002)";
    //static char fmt_9005[] = "(/\002 ************* CYCLING..... **********"
	//    "***\002/)";
   // static char fmt_9000[] = "(\002 CYCLE  0:\002,3x,\002RES=\002,d9.3)";
    //static char fmt_9040[] = "(/\002 CYCLING BETWEEN GRIDS 1 AND\002,i3"
	  //  ",\002:\002/)";
   // static char fmt_9010[] = "(\002 CYCLE \002,i2,\002:   RESCG=\002,d9.3"
	//    ",\002   RES=\002,d9.3,\002   CFAC=\002,d9.3)";
    //static char fmt_9020[] = "(\002 CYCLE \002,i2,\002:   RES=\002,d9.3,\002"
	//    "   CFAC=\002,d9.3)";

    /* System generated locals */
    integer i__1=0, i__2=0;
    doublereal d__1=0, d__2=0, d__3=0;

    /* Builtin functions */
   // integer s_wsfe(cilist *), e_wsfe(void), i_sign(integer *, integer *), 
	 //   do_fio(integer *, char *, ftnlen);
	integer i_sign(integer *, integer *);

    /* Local variables */
    integer i__=0,  m=0;
	integer l=0, n=0;
    extern /* Subroutine */ integer cg_(integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, /*real*/unsigned int *, 
	    integer *, integer *);
    doublereal fac=0.0, ama=0.0;
    extern /* Subroutine */ integer cyc_(integer *, integer *, integer *, integer 
	    *, integer *, integer *, integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, /*real*/ unsigned int *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, doublereal *, doublereal *, 
	    integer *);
    integer nsc=0;
    doublereal cfac=0.0;
    extern /* Subroutine */ integer idec_(integer *, integer *, integer *, 
	    integer *);
    integer igam=0, icgr=0;
    doublereal fmax=0.0, epsi=0.0;
    integer msel=0, iter=0, nrcx=0, nrdx=0;
    doublereal umax=0.0;
    integer nrux=0;
    doublereal rescg=0.0;
    extern /* Subroutine */ integer resid_(integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, integer *, integer *, 
	    integer *, integer *, integer *);
    doublereal epsil=0.0;
    integer iconv=0;
    extern /* Subroutine */ integer usave_(integer *, integer *, doublereal *, 
	    integer *, integer *, integer *, integer *, /*real*/unsigned int *, integer *, 
	    integer *, integer *, integer *, integer *);
    integer ncycle=0, ndigit=0, nrdlen=0;
    doublereal resold=0.0;
	integer nrulen = 0, mfirst = 0, nrdtyp[10] = { 0 }, nrutyp[10] = { 0 };

    /* Fortran I/O blocks */
   // static cilist io___240 = { 0, 0, 0, fmt_9050, 0 };
    //static cilist io___241 = { 0, 0, 0, fmt_9060, 0 };
   // static cilist io___264 = { 0, 0, 0, fmt_9005, 0 };
   // static cilist io___265 = { 0, 0, 0, fmt_9000, 0 };
  //  static cilist io___269 = { 0, 0, 0, fmt_9040, 0 };
   // static cilist io___270 = { 0, 0, 0, fmt_9040, 0 };
   // static cilist io___273 = { 0, 0, 0, fmt_9010, 0 };
    //static cilist io___274 = { 0, 0, 0, fmt_9020, 0 };



/*     SOLUTION PHASE OF AMG1R5 */


/* ===> TEST OF AVAILABLE STORAGE */

    /* Parameter adjustments */
    --resi;
    --time;
    --iarr;
    --nstcol;
    --ifg;
    --icg;
    --imaxw;
    --iminw;
    --imax;
    --imin;
    --iw;
    --ja;
    --ia;
    --f;
    --u;
    --a;

    /* Function Body */
    if (*ndu < imax[*levels]) {
	//io___240.ciunit = *ium;
	//s_wsfe(&io___240);
	//e_wsfe();
	printf("( *** ERROR IN SOLVE: NDU TOO SMALL ***)\n");
	*ierr = 4;
	return 0;
    }
    if (*ndf < imax[*levels]) {
	//io___241.ciunit = *ium;
	//s_wsfe(&io___241);
	//e_wsfe();
	printf("( *** ERROR IN SOLVE: NDF TOO SMALL ***)\n");
	*ierr = 5;
	return 0;
    }

    m = *levels;
    *ncyc0 = 0;
    for (n = 11; n <= 20; ++n) {
	//time[n] = 0.f;
		time[n]=0;
/* L5: */
    }
    if (*eps != 0.) {
	epsi = *eps;
    } else {
	epsi = 1e-12;
    }

/* ===> DECOMPOSE MADAPT */

    if (*madapt != 0) {
	idec_(madapt, &c__2, &ndigit, &iarr[1]);
	msel = iarr[1];
	if (msel == 2) {
	    if (iarr[2] != 0) {
		fac = static_cast<doublereal>(iarr[2]);
		for (i__ = 1; i__ <= 100; ++i__) {
		    fac /= 10.;
		    if (fac <= 1.) {
			goto L9;
		    }
/* L8: */
		}
	    } else {
		fac = .7;
	    }
	}
    } else {
	msel = 2;
	fac = .7;
    }

/* ===> DECOMPOSE NCYC */

L9:
    if (*ncyc != 0) {
	i__1 = abs(*ncyc);
	idec_(&i__1, &c__4, &ndigit, &iarr[1]);
	igam = i_sign(&iarr[1], ncyc);
	icgr = iarr[2];
	iconv = iarr[3];
	ncycle = iarr[4];
	if (ncycle == 0) {
	    return 0;
	}
    } else {
	igam = 1;
	icgr = 0;
	iconv = 1;
	ncycle = 10;
    }

/* ===> SET EPSI ACCORDING TO CONVERGENCE CRITERION GIVEN BY ICONV */

    if (iconv != 3) {
	if (iconv == 4) {
	    ama = 0.;
	    i__1 = imax[1];
	    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MAX */
		d__1 = ama, d__2 = a[ia[i__]];
		ama = myr_max(d__1,d__2);
/* L6: */
	    }
	    epsi *= ama;
	}
    } else {
	fmax = 0.;
	i__1 = imax[1];
	for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MAX */
	    d__2 = fmax, d__3 = (d__1 = f[i__], fabs(d__1));
	    fmax = myr_max(d__2,d__3);
/* L7: */
	}
	epsi *= fmax;
    }

/* ===> DECOMPOSE NRD */

    if (*nrd != 0) {
	idec_(nrd, &c__9, &ndigit, nrdtyp);
	nrdx = nrdtyp[1];
	nrdlen = ndigit - 2;
	i__1 = nrdlen;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    nrdtyp[i__ - 1] = nrdtyp[i__ + 1];
/* L10: */
	}
    } else {
	nrdx = 1;
	nrdlen = 2;
	nrdtyp[0] = 3;
	nrdtyp[1] = 1;
    }

/* ===> DECOMPOSE NRU */

    if (*nru != 0) {
	idec_(nru, &c__9, &ndigit, nrutyp);
	nrux = nrutyp[1];
	nrulen = ndigit - 2;
	i__1 = nrulen;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    nrutyp[i__ - 1] = nrutyp[i__ + 1];
/* L40: */
	}
    } else {
	nrux = 1;
	nrulen = 2;
	nrutyp[0] = 3;
	nrutyp[1] = 1;
    }

/* ===> DECOMPOSE NSOLCO */

    if (*nsolco != 0) {
	idec_(nsolco, &c__2, &ndigit, &iarr[1]);
	nsc = iarr[1];
	nrcx = iarr[2];

/* ===> IN CASE OF YALE-SMP COARSE GRID SOLUTION, DON'T USE COARSEST */
/*     GRID WITH LESS THAN 10 POINTS */

	if (nsc == 2) {
	    for (i__ = m; i__ >= 1; --i__) {
		l = i__;
		if (imax[i__] - imin[i__] >= 9) {
		    goto L60;
		}
/* L50: */
	    }
L60:
	    m = i__;
	    *levels = i__;
	}
    } else {
	nsc = 1;
	nrcx = 0;
    }

/* ===> CYCLING */

/* L100: */
    if (*iout != 0) {
	resid_(&c__1, res0, &a[1], &u[1], &f[1], &ia[1], &ja[1], &iw[1], &
		imin[1], &imax[1], &iminw[1]);
	if (*iout == 3) {
		if (yes_print_amg) {
			printf("( ************* CYCLING..... *************)\n");
			printf("( CYCLE  0:,3x,RES=,%1.4f)\n", *res0);
		}
		
	}
	
	resold = *res0;
    }

    i__1 = ncycle;
    for (iter = 1; iter <= i__1; ++iter) {
	usave_(&c__1, &icgr, &u[1], &imin[1], &imax[1], ndu, &m, &time[1], 
		ierr, ium, mdu, ndf, mdf);
	cyc_(&c__1, &nrdx, nrdtyp, &nrdlen, &nrcx, &nrux, nrutyp, &nrulen, &
		igam, &a[1], &u[1], &f[1], &ia[1], &ja[1], &iw[1], &imin[1], &
		imax[1], &iminw[1], &imaxw[1], &ifg[1], &icg[1], &nstcol[1], &
		iarr[1], &time[1], irow0, &m, ium, ierr, &iter, &nsc, nda, 
		ndja, mda, mdja, &msel, &fac, &resi[1], levels);
	cg_(&c__1, &icgr, &iter, &a[1], &u[1], &f[1], &ia[1], &ja[1], &iw[1], 
		&imin[1], &imax[1], &iminw[1], &m, &time[1], ierr, ium);
	if (*ierr > 0) {
	    return 0;
	}
	if (iter == 1) {
	    mfirst = m;
	    if (*iout == 3) {
		
			if (yes_print_amg) {
#if doubleintprecision == 1
				printf("(/ CYCLING BETWEEN GRIDS 1 AND ,%lld,)\n", m);
#else
				printf("(/ CYCLING BETWEEN GRIDS 1 AND ,%d,)\n", m);
#endif
		         
			}
	    }
	} else if (*iout == 3 && m != mfirst) {
	    mfirst = m;
	  
		if (yes_print_amg) {
#if doubleintprecision == 1
			printf("(/ CYCLING BETWEEN GRIDS 1 AND ,%lld,)\n", m);
#else
			printf("(/ CYCLING BETWEEN GRIDS 1 AND ,%d,)\n", m);
#endif
                
		}
		
		
	}
	if (*iout == 3 || iconv != 1) {
	    resid_(&c__1, res, &a[1], &u[1], &f[1], &ia[1], &ja[1], &iw[1], &
		    imin[1], &imax[1], &iminw[1]);
	}
	*ncyc0 = iter;
	if (*iout != 3) {
	    goto L110;
	}
	cfac = *res / (resold + 1e-40);
	resold = *res;
	if (1 == m) {
	    goto L150;
	}
	resid_(&c__2, &rescg, &a[1], &u[1], &f[1], &ia[1], &ja[1], &iw[1], &
		imin[1], &imax[1], &iminw[1]);
	
	if (yes_print_amg) {
#if doubleintprecision == 1
		printf("( CYCLE ,%lld,:   RESCG=,%1.4f\n", iter, rescg);
#else
		printf("( CYCLE ,%d,:   RESCG=,%1.4f\n", iter, rescg);
#endif
	     
	     printf(",   RES=,%1.4f,   CFAC=,%1.4f)\n", *res, cfac);
	}
	goto L110;
L150:
	if (yes_print_amg) {
#if doubleintprecision == 1
		printf("( CYCLE ,%lld,:   RES=,%1.4f,\n", iter, *res);
#else
		printf("( CYCLE ,%d,:   RES=,%1.4f,\n", iter, *res);
#endif
	     
	     printf("   CFAC=,%1.4f)\n",cfac);
	}

L110:
	if (iconv == 1) {
	    goto L120;
	}
	epsil = epsi;
	if (iconv != 4) {
	    goto L115;
	}
	umax = 0.;
	i__2 = imax[1];
	for (i__ = imin[1]; i__ <= i__2; ++i__) {
/* Computing MAX */
	    d__2 = umax, d__3 = (d__1 = u[i__], fabs(d__1));
	    umax = myr_max(d__2,d__3);
/* L160: */
	}
	epsil = epsi * umax;
L115:
	if (*res < epsil) {
	    goto L170;
	}
L120:
	;
    }
L170:
    if (*iout != 3 && *iout != 0) {
	resid_(&c__1, res, &a[1], &u[1], &f[1], &ia[1], &ja[1], &iw[1], &imin[
		1], &imax[1], &iminw[1]);
    }
    return 0;

/* L9030: */
} /* solve_ */


/* ....................................................................... */

/*     CYC                                                    SUBROUTINE */

/* ....................................................................... */

/* Subroutine */ integer cyc_(integer *l, integer *nrdx, integer *nrdtyp, integer 
	*nrdlen, integer *nrcx, integer *nrux, integer *nrutyp, integer *
	nrulen, integer *igam, doublereal *a, doublereal *u, doublereal *f, 
	integer *ia, integer *ja, integer *iw, integer *imin, integer *imax, 
	integer *iminw, integer *imaxw, integer *ifg, integer *icg, integer *
	nstcol, integer *ng, /*real*/ unsigned int *time, integer *irow0, integer *m, integer *
	ium, integer *ierr, integer *iter, integer *nsc, integer *nda, 
	integer *ndja, integer *mda, integer *mdja, integer *msel, doublereal 
	*fac, doublereal *resi, integer *levels)
{
    /* System generated locals */
    integer i__1=0, i__2=0;

    /* Local variables */
    integer k=0, n=0, nl=0, ifi=0;
    doublereal res=0.0;
    integer ifac=0;
    extern /* Subroutine */ integer inta_(integer *, doublereal *, doublereal *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, /*real*/unsigned int *), resc_(integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, /*real*/unsigned int *);
	unsigned int told=0;
    integer mink=0;
    extern /* Subroutine */ integer relx_(integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, /*real*/ unsigned int *);
	unsigned int tnew=0;
    extern /* Subroutine */ integer nrmu_(integer *, doublereal *, integer *, 
	    integer *, integer *), putz_(integer *, doublereal *, integer *, 
	    integer *, /*real*/unsigned int *),/* ctime_(real *),*/ resid_(integer *, doublereal *
	    , doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    integer *, integer *, integer *, integer *);
    integer nptsf=0;
    extern /* Subroutine */ integer coarse_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, integer *, integer *, integer 
	    *, integer *, integer *, integer *, /*real*/ unsigned int *, integer *, integer *, 
	    integer *, integer *, integer *), vscale_(integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, /*real*/unsigned int *);
    integer ivstar=0;


/*     PERFORMS ONE AMG CYCLE WITH GRID L AS FINEST GRID */



/* ===> DURING FIRST CYCLE: INITIALIZE PARAMETERS CONTROLLING YALE-SMP */
/*     FACTORIZATION AND ADAPTIVE DETERMINATION OF COARSEST GRID: */
/*     IFAC=1: ON NEXT CALL OF YALE-SMP FACTORIZE MATRIX */
/*     IFI =1: ON FIRST RETURN TO NEXT TO COARSEST GRID AFTER COARSE */
/*             GRID SOLUTION COMPARE RESIDUAL WITH THE RESIDUAL ON THE */
/*             SAME GRID BEFORE COARSE GRID SOLUTION. IF REDUCTION OF */
/*             RESIDUAL NOT SATISFYING, REDUCE NUMBER OF GRIDS USED IN */
/*             CYCLING BY ONE AND REPETE THE PROCESS WITH THE NOW */
/*             COARSEST GRID. */
/*     NPTSF:  NUMBER OF GRID POINTS ON FINEST GRID USED IN CYCLE, */
/*             DIVIDED BY 10. ONLY GRIDS WITH LESS THEN NPTSF POINTS */
/*             ARE ALLOWED TO BECOME COARSEST GRID DURING ADAPTIVE */
/*             COARSE GRID DETERMINATION. */

    /* Parameter adjustments */
    --resi;
    --time;
    --ng;
    --nstcol;
    --icg;
    --ifg;
    --imaxw;
    --iminw;
    --imax;
    --imin;
    --iw;
    --ja;
    --ia;
    --f;
    --u;
    --a;
    --nrutyp;
    --nrdtyp;

    /* Function Body */
    if (*iter == 1) {
	ifac = 1;
    }
    if (*msel != 2) {
	ifi = 0;
	nptsf = 0;
    } else {
	mink = 1000000;
	ifi = 1;
	nptsf = (imax[*l] - imin[*l] + 1) / 10;
    }
    if (*l < *m) {
	goto L100;
    }

/* ===> ONE GRID ONLY */

    coarse_(m, &ifac, nsc, nrcx, ium, ierr, &a[1], &u[1], &f[1], &ia[1], &ja[
	    1], &iw[1], &imin[1], &imax[1], &iminw[1], &icg[1], &nstcol[1], &
	    time[1], nda, ndja, mda, mdja, irow0);
    goto L1000;

/* ===> MORE THAN ONE GRID */

L100:
    i__1 = *m;
    for (k = *l; k <= i__1; ++k) {
	ng[k] = 0;
/* L110: */
    }
    ivstar = 3 - *igam;
    k = *l;

/* ===> RELAX (DOWNWARDS) */

L150:
    i__1 = *nrdx;
    for (n = 1; n <= i__1; ++n) {
	i__2 = *nrdlen;
	for (nl = 1; nl <= i__2; ++nl) {
	    relx_(&k, &nrdtyp[nl], &a[1], &u[1], &f[1], &ia[1], &ja[1], &iw[1]
		    , &imin[1], &imax[1], &iminw[1], &icg[1], &nstcol[1], &
		    time[1]);
/* L160: */
	}
/* L170: */
    }
    if (ifi != 1 || imax[k] - imin[k] >= nptsf) {
	goto L190;
    }

/* ===> MINK: LOWEST GRID NUMBER FOR WHICH RESIDUAL IS STORED DURING */
/*           FIRST DOWNWARDS RELAXATION */

    if (mink == 1000000) {
	mink = k;
    }
	told=clock();
    resid_(&k, &resi[k], &a[1], &u[1], &f[1], &ia[1], &ja[1], &iw[1], &imin[1]
	    , &imax[1], &iminw[1]);
	tnew=clock();
    time[15] = time[15] + tnew - told;
L190:
    ++ng[k];
    ++k;
L195:
    putz_(&k, &u[1], &imin[1], &imax[1], &time[1]);
    resc_(&k, &a[1], &u[1], &f[1], &ia[1], &ja[1], &iw[1], &imin[1], &imax[1],
	     &iminw[1], &imaxw[1], &ifg[1], &time[1]);
    if (k < *m) {
	goto L150;
    }

/* ===> SOLVE ON COARSEST GRID */

    coarse_(m, &ifac, nsc, nrcx, ium, ierr, &a[1], &u[1], &f[1], &ia[1], &ja[
	    1], &iw[1], &imin[1], &imax[1], &iminw[1], &icg[1], &nstcol[1], &
	    time[1], nda, ndja, mda, mdja, irow0);

/* ===> RELAX (UPWARDS) */

L200:
    vscale_(&k, &ivstar, &a[1], &u[1], &f[1], &ia[1], &ja[1], &iw[1], &imin[1]
	    , &imax[1], &iminw[1], &time[1]);
    --k;
    inta_(&k, &a[1], &u[1], &ia[1], &ja[1], &iw[1], &imin[1], &imax[1], &
	    iminw[1], &imaxw[1], &ifg[1], &time[1]);
    i__1 = *nrux;
    for (n = 1; n <= i__1; ++n) {
	i__2 = *nrulen;
	for (nl = 1; nl <= i__2; ++nl) {
	    relx_(&k, &nrutyp[nl], &a[1], &u[1], &f[1], &ia[1], &ja[1], &iw[1]
		    , &imin[1], &imax[1], &iminw[1], &icg[1], &nstcol[1], &
		    time[1]);
/* L210: */
	}
/* L215: */
    }
    if (ifi != 1 || k < mink) {
	goto L219;
    }

/* ===> ON FIRST RETURN TO NEXT TO COARSEST GRID COMPARE RESIDUAL WITH */
/*     THE PREVIOUS ONE */

	told=clock();
    resid_(&k, &res, &a[1], &u[1], &f[1], &ia[1], &ja[1], &iw[1], &imin[1], &
	    imax[1], &iminw[1]);
	tnew=clock();
    time[15] = time[15] + tnew - told;

/* ===> IF RESIDUAL REDUCTION SATISFYING: COARSE GRID ADAPTION FINISHED */

    if (res < resi[k] * *fac) {
	ifi = 0;
    } else {
	if (*nsc == 2) {
	    *levels = k;
	}
	*m = k;
	ifac = 1;
	goto L195;
    }
L219:
    if (k == *l) {
	goto L1000;
    }

/* ===> GRID SWITCHING CORRESPONDING TO IGAM */

    if (*igam >= 3) {
	goto L220;
    }
    if (*igam >= 0) {
	goto L200;
    }
    if (k == *l + 1 && ng[k] < abs(*igam)) {
	goto L150;
    }
    goto L200;
L220:
    if (ng[k] < 2) {
	goto L150;
    }
    if (*igam == 4) {
	ng[k] = 0;
    }
    goto L200;

/* ===> RETURN */

L1000:
    nrmu_(l, &u[1], &imin[1], &imax[1], irow0);
    return 0;
} /* cyc_ */


/* ....................................................................... */

/*     COARSE                                                SUBROUTINE */

/* ....................................................................... */

/* Subroutine */ integer coarse_(integer *m, integer *ifac, integer *nsc, integer 
	*nrcx, integer *ium, integer *ierr, doublereal *a, doublereal *u, 
	doublereal *f, integer *ia, integer *ja, integer *iw, integer *imin, 
	integer *imax, integer *iminw, integer *icg, integer *nstcol, /*real*/ unsigned int *
	time, integer *nda, integer *ndja, integer *mda, integer *mdja, 
	integer *irow0)
{
    /* Format strings */
    //static char fmt_9020[] = "(\002 --- WARNG IN COARSE: NO YALE-SMP BECAUSE"
	//    " NDJA TOO \002,\002SMALL\002)";
    //static char fmt_9030[] = "(\002 --- WARNG IN COARSE: NO YALE-SMP BECAUSE"
	//    " NDA TOO SMALL\002)";
    //static char fmt_9040[] = "(\002 --- WARNG IN COARSE: NO YALE-SMP BECAUSE"
	//    " OF ERROR IN \002,\002FACTORIZATION,\002/\002     CODE=\002,i8)";

    /* System generated locals */
    integer i__1=0, i__2=0;
    doublereal d__1=0.0, d__2=0.0, d__3=0.0;
	


    /* Builtin functions */
    //integer s_wsfe(cilist *), e_wsfe(void), do_fio(integer *, char *, ftnlen);

    /* Local variables */
    integer i__=0, j=0, ii=0, jj=0, is=0, js=0, np=0, ihi=0, jhi=0, ilo=0, jlo=0, esp=0, nsp=0, 
	    flag__=0;
    doublereal fmax=0.0;
    integer path=0;
    doublereal aaux=0.0;
	unsigned int told=0;
    integer iter=0, iaux=0;
    extern /* Subroutine */ integer ndrv_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, doublereal /*integer*/ *, doublereal *, integer *, 
	    integer *, integer *);
	 extern /* Subroutine */ integer relx_(integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, /*real*/ unsigned int *);
    integer jpos=0;
	unsigned int tnew=0;
    extern /* Subroutine */ integer  resid_(integer *, doublereal *
	    , doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    integer *, integer *, integer *, integer *);
    doublereal resold=0.0, resnew=0.0;
    integer npoint=0;

    /* Fortran I/O blocks */
    //static cilist io___296 = { 0, 0, 0, fmt_9020, 0 };
    //static cilist io___309 = { 0, 0, 0, fmt_9030, 0 };
   // static cilist io___310 = { 0, 0, 0, fmt_9040, 0 };



/*     SOLVES ON COARSEST GRID, EITHER WITH GAUSS-SEIDEL RELAXATION */
/*     (NSC=1) OR WITH THE YALE-SMP DIRECT SOLVER NDRV (NSC=2) */


/* ===> CONV: IF COARSE GRID SOLUTION IS DONE WITH GS-RELAXATION AND */
/*     NRCX=0, AS MANY GS-SWEEPS ARE PERFORMED AS ARE NECESSARY TO RE- */
/*     DUCE THE RESIDUAL BY THE FACTOR CONV */


    /* Parameter adjustments */
    --time;
    --nstcol;
    --icg;
    --iminw;
    --imax;
    --imin;
    --iw;
    --ja;
    --ia;
    --f;
    --u;
    --a;

    /* Function Body */
    if (*nsc != 2) {
	goto L400;
    }

/* ===> SOLUTION WITH YALE-SMP */

	told=clock();
    ilo = imin[*m];
    jlo = ia[ilo];
    if (*ifac == 1) {

/* ===> FIRST CALL ON GRID M, FIRST FACTORIZE MATRIX */

	ihi = imax[*m];
	jhi = iw[iminw[*m]] - 1;
	np = ihi - ilo + 1;
	is = ilo - 1;
	js = jlo - 1;

/* ===>   TEST OF AVAILABLE WORK SPACE */

	if (jhi + np * 3 > *ndja) {
	    *nsc = 1;
	    //io___296.ciunit = *ium;
	    //s_wsfe(&io___296);
	    //e_wsfe();

		printf("( --- WARNG IN COARSE: NO YALE-SMP BECAUSE\n");
	    printf(" NDJA TOO ,SMALL)\n");

	    *ierr = -3;
	    goto L400;
	}

/* ===>   INITIALISATION OF YALE-SMP POINTER VECTORS */

	i__1 = np;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    ja[jhi + i__] = i__;
	    ja[jhi + np + i__] = i__;
	    ja[jhi + (np << 1) + i__] = i__;
/* L310: */
	}
	if (*irow0 != 1) {

/* ===>   COARSE GRID OPERATOR REGULAR, SHIFT CONTENTS OF POINTER */
/*       VECTORS IA AND JA */

	    i__1 = ihi;
	    for (i__ = ilo; i__ <= i__1; ++i__) {
		ia[i__] -= js;
/* L320: */
	    }
	    iaux = ia[ihi + 1];
	    ia[ihi + 1] = iw[iminw[*m]] - js;
	    i__1 = jhi;
	    for (i__ = jlo; i__ <= i__1; ++i__) {
		ja[i__] -= is;
/* L300: */
	    }
	    npoint = np;
	} else {

/* ===>   COARSE GRID OPERATOR HAS ROWSUM ZERO, ELIMINATE LAST SOLUTION */
/*       COMPONENT BY SETTING U(IHI) TO ZERO AND CANCELLING THE COR- */
/*       RESPONDING ENTRIES IN A AND JA, RESPECTIVELY. THE CANCELLED */
/*       ENTRIES ARE STORED BEFORE THE LAST ROW IN A AND JA. */

	    iaux = ia[ihi];

/* ===>     JPOS: POINTER TO POSITION IN A AND JA TO CONTAIN NEXT */
/*               ELIMINATED ENTRY */

	    jpos = iaux - 1;
	    j = ia[ilo];
	    i__1 = ihi - 1;
	    for (i__ = ilo; i__ <= i__1; ++i__) {
		ia[i__] = j - js;
		ja[j] = i__ - is;
		++j;
L13:
		if (j == ia[i__ + 1]) {
		    goto L20;
		}
		if (ja[j] != ihi) {
		    ja[j] -= is;
		    ++j;
		} else {
		    aaux = a[j];
		    i__2 = jpos - 1;
		    for (jj = j; jj <= i__2; ++jj) {
			a[jj] = a[jj + 1];
			ja[jj] = ja[jj + 1];
/* L11: */
		    }
		    i__2 = ihi;
		    for (ii = i__ + 1; ii <= i__2; ++ii) {
			--ia[ii];
/* L12: */
		    }
		    a[jpos] = aaux;
		    ja[jpos] = i__;
		    --jpos;
		}
		goto L13;
L20:
		;
	    }
	    ia[ihi] -= js;

/* ===>     DECREASE NUMBER OF POINTS BY ONE AND SET LAST SOLUTION */
/*         COMPONENT TO ZERO */

	    npoint = np - 1;
	    u[ihi] = 0.;
	}
	nsp = *nda - jhi;
	path = 1;

	


	ndrv_(&npoint, &ja[jhi + 1], &ja[jhi + np + 1], &ja[jhi + (np << 1) + 
		1], &ia[ilo], &ja[jlo], &a[jlo], &f[ilo], &u[ilo], &nsp, 
		&a[jhi + 1], &a[jhi + 1], &esp, &path, &flag__);

		

/* ===>   RESTORE PREVIOUS VALUES FOR IA, JA AND A, IN PARTICULAR PUT */
/*       BACK ELIMINATED MATRIX ENTRIES, IF ROWSUM ZERO */

	i__1 = ihi;
	for (i__ = ilo; i__ <= i__1; ++i__) {
	    ia[i__] += js;
/* L330: */
	}
	if (*irow0 != 1) {
	    ia[ihi + 1] = iaux;
	    i__1 = jhi;
	    for (i__ = jlo; i__ <= i__1; ++i__) {
		ja[i__] += is;
/* L335: */
	    }
	} else {
	    i__1 = jpos;
	    for (i__ = jlo; i__ <= i__1; ++i__) {
		ja[i__] += is;
/* L25: */
	    }
	    i__1 = iaux - 1;
	    for (j = jpos + 1; j <= i__1; ++j) {
		aaux = a[j];
		i__ = ja[j];
		i__2 = ihi;
		for (ii = i__ + 1; ii <= i__2; ++ii) {
		    ++ia[ii];
/* L27: */
		}
		i__2 = ia[i__ + 1];
		for (jj = j; jj >= i__2; --jj) {
		    a[jj] = a[jj - 1];
		    ja[jj] = ja[jj - 1];
/* L30: */
		}
		a[ia[i__ + 1] - 1] = aaux;
		ja[ia[i__ + 1] - 1] = ihi;
/* L40: */
	    }
	}

/* ===>   IF AN ERROR OCCURED DURING EXECUTION OF NDRV, SOLVE WITH */
/*       GAUSS-SEIDEL RELAXATION */

	if (flag__ != 0) {
	    *nsc = 1;
	    if (esp < 0) {
		

		printf("( --- WARNG IN COARSE: NO YALE-SMP BECAUSE\n");
	    printf(" NDA TOO SMALL)\n");

		*ierr = -1;
	    } else {
		

		printf("( --- WARNG IN COARSE: NO YALE-SMP BECAUSE");
#if doubleintprecision == 1
		printf(" OF ERROR IN ,FACTORIZATION,     CODE=,%lld)", flag__);
#else
		printf(" OF ERROR IN ,FACTORIZATION,     CODE=,%d)", flag__);
#endif
	   

		*ierr = -32;
	    }
	    goto L400;
	} else {

/* ===>     FACTORIZATION SUCCESSFULL, UPDATE LIMITS OF USED STORAGE */

/* Computing MAX */
	    i__1 = *mda, i__2 = *nda - esp;
	    *mda = myi_max(i__1,i__2);
/* Computing MAX */
	    i__1 = *mdja, i__2 = jhi + np * 3;
	    *mdja = myi_max(i__1,i__2);
	    *ifac = 0;
	}
    } else {

/* ===>   FACTORIZATION ALLREADY DONE */

	if (*irow0 != 1) {
	    npoint = np;
	} else {
	    npoint = np - 1;
	    u[ihi] = 0.;
	}
	path = 3;

	

	ndrv_(&npoint, &ja[jhi + 1], &ja[jhi + np + 1], &ja[jhi + (np << 1) + 
		1], &ia[ilo], &ja[jlo], &a[jlo], &f[ilo], &u[ilo], &nsp, &a[jhi + 1], 
		&a[jhi + 1], &esp, &path, &flag__);
    }

	

/* ===> UPDATE TIME COUNTER */

	tnew=clock();
    time[17] = time[17] + tnew - told;
    goto L190;

/* ===> SOLUTION WITH GAUSS-SEIDEL RELAXATION */

L400:
    if (*nrcx != 0) {
	i__1 = *nrcx;
	for (iter = 1; iter <= i__1; ++iter) {
	    relx_(m, &c__2, &a[1], &u[1], &f[1], &ia[1], &ja[1], &iw[1], &
		    imin[1], &imax[1], &iminw[1], &icg[1], &nstcol[1], &time[
		    1]);
/* L180: */
	}
    } else {

/* ===> IF NRCX=0: REDUCE RESIDUAL ON COARSEST GRID BY FACTOR CONV */
/*                (IF NOT YET IN THE RANGE OF THE TRUNCATION ERROR) */

	told=clock();

/* ===>   CALCULATE SUPREMUM NORM OF RIGHT HAND SIDE */

	fmax = 0.;
	i__1 = imax[*m];
	for (i__ = imin[*m]; i__ <= i__1; ++i__) {
/* Computing MAX */
	    d__2 = fmax, d__3 = (d__1 = f[i__], fabs(d__1));
	    fmax = myr_max(d__2,d__3);
/* L181: */
	}
	resid_(m, &resold, &a[1], &u[1], &f[1], &ia[1], &ja[1], &iw[1], &imin[
		1], &imax[1], &iminw[1]);
	tnew=clock();
	time[15] = time[15] + tnew - told;
/* Computing MAX */
	d__1 = resold * .01, d__2 = fmax * 1e-12;
	resold = myr_max(d__1,d__2);
	for (i__ = 1; i__ <= 10; ++i__) {
	    for (j = 1; j <= 10; ++j) {
		relx_(m, &c__2, &a[1], &u[1], &f[1], &ia[1], &ja[1], &iw[1], &
			imin[1], &imax[1], &iminw[1], &icg[1], &nstcol[1], &
			time[1]);
/* L185: */
	    }
		told=clock();
	    resid_(m, &resnew, &a[1], &u[1], &f[1], &ia[1], &ja[1], &iw[1], &
		    imin[1], &imax[1], &iminw[1]);
		tnew=clock();
	    time[15] = time[15] + tnew - told;
	    if (resnew <= resold) {
		goto L190;
	    }
/* L187: */
	}
    }

/* ===> COMPUTE RESIDUAL AFTER SOLUTION ON GRID M */

L190:
    return 0;

} /* coarse_ */


/* ....................................................................... */

/*     RELX                                                  SUBROUTINE */

/* ....................................................................... */

/* Subroutine */ integer relx_(integer *k, integer *irel, doublereal *a, 
	doublereal *u, doublereal *f, integer *ia, integer *ja, integer *iw, 
	integer *imin, integer *imax, integer *iminw, integer *icg, integer *
	nstcol, /*real*/ unsigned int *time)
{
    /* System generated locals */
    integer i__1=0, i__2=0, i__S0=0;

    /* Local variables */
    integer i__=0, j=0;
    doublereal s=0.0;
	unsigned int told=0;
    integer iaux=0;
	unsigned int tnew=0;


/*     PERFORMS ONE (PARTIAL) GAUSS-SEIDEL SWEEP ON GRIK K: */

/*     IREL = 1:   PARTIAL GAUSS-SEIDEL SWEEP (ONLY F-POINTS) */
/*          = 2:   FULL GAUSS-SEIDEL SWEEP (ALL POINTS) */
/*          = 3:   PARTIAL GAUSS-SEIDEL SWEEP (ONLY C-POINTS) */
/*          = 4:   FULL SWEEP: FF -- C -- COLORS (HIGHEST FIRST) */


    /* Parameter adjustments */
    --time;
    --nstcol;
    --icg;
    --iminw;
    --imax;
    --imin;
    --iw;
    --ja;
    --ia;
    --f;
    --u;
    --a;

    /* Function Body */
	told=clock();
    iaux = ia[imax[*k] + 1];
    ia[imax[*k] + 1] = iw[iminw[*k]];
    switch (*irel) {
	case 1:  goto L100;
	case 2:  goto L200;
	case 3:  goto L300;
	case 4:  goto L400;
    }
    goto L200;

/* ===> F-RELAXATION */

L100:
    i__1 = imax[*k];
    for (i__ = imin[*k]; i__ <= i__1; ++i__) {
	if (icg[i__] > 0) {
	    goto L120;
	}
	s = f[i__];
	i__2 = ia[i__ + 1] - 1;
	for (j = ia[i__] + 1; j <= i__2; ++j) {
	    s -= a[j] * u[ja[j]];
/* L110: */
	}
	u[i__] = s / a[ia[i__]];
L120:
	;
    }
    goto L1000;

/* ===> FULL GS RELAXATION */

L200:
    i__1 = imax[*k];
    for (i__ = imin[*k]; i__ <= i__1; ++i__) {
	s = f[i__];
	i__2 = ia[i__ + 1] - 1;
	for (j = ia[i__] + 1; j <= i__2; ++j) {
	    s -= a[j] * u[ja[j]];
/* L210: */
	}
	u[i__] = s / a[ia[i__]];
/* L220: */
    }
    goto L1000;

/* ===> C-RELAXATION */

L300:
    i__1 = imax[*k];
	i__S0 = imin[*k];

#pragma omp for 
    for (integer i__loc = i__S0; i__loc <= i__1; ++i__loc) {
		if (icg[i__loc] > 0) {
			doublereal s_loc = f[i__loc];
			integer i__2loc = ia[i__loc + 1] - 1;
			for (integer j_loc = ia[i__loc] + 1; j_loc <= i__2loc; ++j_loc) {
				s_loc -= a[j_loc] * u[ja[j_loc]];
				/* L310: */
			}
			u[i__loc] = s_loc / a[ia[i__loc]];
		}
    }
    goto L1000;

/* ===> FF-RELAXATION */

L400:
    i__1 = imax[*k];
    for (i__ = imin[*k]; i__ <= i__1; ++i__) {
	if (icg[i__] != 0) {
	    goto L420;
	}
	s = f[i__];
	i__2 = ia[i__ + 1] - 1;
	for (j = ia[i__] + 1; j <= i__2; ++j) {
	    s -= a[j] * u[ja[j]];
/* L410: */
	}
	u[i__] = s / a[ia[i__]];
L420:
	;
    }

/* ===> C-RELAXATION */

    i__1 = imax[*k];
    for (i__ = imin[*k]; i__ <= i__1; ++i__) {
	if (icg[i__] <= 0) {
	    goto L440;
	}
	s = f[i__];
	i__2 = ia[i__ + 1] - 1;
	for (j = ia[i__] + 1; j <= i__2; ++j) {
	    s -= a[j] * u[ja[j]];
/* L430: */
	}
	u[i__] = s / a[ia[i__]];
L440:
	;
    }

/* ===> MORE-COLOR RELAXATION */

    i__ = nstcol[*k];
L470:
#if doubleintprecision == 1
	// 1000000M
	if (i__ >= 1000000000000) {
		goto L1000;
	}
#else
	// 100M
	if (i__ >= 100000000) {
		goto L1000;
	}
#endif
    
    s = f[i__];
    i__1 = ia[i__ + 1] - 1;
    for (j = ia[i__] + 1; j <= i__1; ++j) {
	s -= a[j] * u[ja[j]];
/* L480: */
    }
    u[i__] = s / a[ia[i__]];
    i__ = -icg[i__];
    goto L470;

L1000:
	tnew=clock();
    time[13] = time[13] + tnew - told;
    ia[imax[*k] + 1] = iaux;
    return 0;
} /* relx_ */


/* ....................................................................... */

/*     VSCALE                                                SUBROUTINE */

/* ....................................................................... */

/* Subroutine */ integer vscale_(integer *k, integer *ivstar, doublereal *a, 
	doublereal *u, doublereal *f, integer *ia, integer *ja, integer *iw, 
	integer *imin, integer *imax, integer *iminw, /*real*/unsigned int  *time)
{
    /* System generated locals */
    integer i__1=0, i__2=0;

    /* Local variables */
    integer i__=0, j=0;
    doublereal s1=0.0, s2=0.0, sa=0.0, fac=0.0;
	unsigned int told=0;
    integer iaux=0;
	unsigned int tnew=0;
   


/*     SCALES ACTUAL APPROXIMATE SOLUTION ON GRID K (V*-CYCLE); SCALING */
/*     IS DONE SUCH THAT ENERGY NORM BECOMES MINIMAL */

/*     NOTE: THIS SCALING MAKES SENSE ONLY FOR SYMMETRIC PROBLEMS */


/* ===> COMPUTATION OF SCALING FACTOR "FAC" */

    /* Parameter adjustments */
    --time;
    --iminw;
    --imax;
    --imin;
    --iw;
    --ja;
    --ia;
    --f;
    --u;
    --a;

    /* Function Body */
    if (*ivstar != 1) {
	return 0;
    }
   
	told=clock();
    iaux = ia[imax[*k] + 1];
    ia[imax[*k] + 1] = iw[iminw[*k]];
    s1 = 0.;
    s2 = 0.;
    i__1 = imax[*k];
    for (i__ = imin[*k]; i__ <= i__1; ++i__) {
	sa = 0.;
	i__2 = ia[i__ + 1] - 1;
	for (j = ia[i__]; j <= i__2; ++j) {
	    sa += a[j] * u[ja[j]];
/* L20: */
	}
	s1 += u[i__] * f[i__];
	s2 += u[i__] * sa;
/* L10: */
    }
    fac = 1.;
    if (s2 != 0.) {
	fac = s1 / s2;
    }

/* ===> SCALING */

    i__1 = imax[*k];
    for (i__ = imin[*k]; i__ <= i__1; ++i__) {
	u[i__] *= fac;
/* L30: */
    }
    ia[imax[*k] + 1] = iaux;
	tnew=clock();
    time[14] = time[14] + tnew - told;
    return 0;
} /* vscale_ */


/* ....................................................................... */

/*     INTA                                                  SUBROUTINE */

/* ....................................................................... */

/* Subroutine */ integer inta_(integer *kf, doublereal *a, doublereal *u, integer 
	*ia, integer *ja, integer *iw, integer *imin, integer *imax, integer *
	iminw, integer *imaxw, integer *ifg, /*real*/ unsigned int *time)
{
    /* System generated locals */
    integer i__1=0, i__2=0;

    /* Local variables */
    integer i__=0, j=0, ic=0, if__=0;
	unsigned int told=0;
    integer iaux=0;
	unsigned int tnew=0;


/*     INTERPOLATES CORRECTION FROM GRID KF+1 TO GRID KF */


/* ===> C->C CONTRIBUTIONS */

    /* Parameter adjustments */
    --time;
    --ifg;
    --imaxw;
    --iminw;
    --imax;
    --imin;
    --iw;
    --ja;
    --ia;
    --u;
    --a;

    /* Function Body */
	told=clock();
    i__1 = imax[*kf + 1];
    for (ic = imin[*kf + 1]; ic <= i__1; ++ic) {
	if__ = ifg[ic];
	u[if__] += u[ic];
/* L50: */
    }

/* ===> C->F CONTRIBUTIONS */

    iaux = iw[imaxw[*kf] + 1];
    iw[imaxw[*kf] + 1] = ia[imin[*kf + 1]];
    i__1 = imaxw[*kf];
    for (i__ = iminw[*kf]; i__ <= i__1; ++i__) {
	if__ = ifg[i__];
	i__2 = iw[i__ + 1] - 1;
	for (j = iw[i__]; j <= i__2; ++j) {
	    u[if__] += a[j] * u[ja[j]];
/* L150: */
	}
/* L100: */
    }
    iw[imaxw[*kf] + 1] = iaux;
	tnew=clock();
    time[11] = time[11] + tnew - told;
    return 0;
} /* inta_ */


/* ....................................................................... */

/*     RESC                                                  SUBROUTINE */

/* ....................................................................... */

/* Subroutine */ integer resc_(integer *kc, doublereal *a, doublereal *u, 
	doublereal *f, integer *ia, integer *ja, integer *iw, integer *imin, 
	integer *imax, integer *iminw, integer *imaxw, integer *ifg, /*real*/ unsigned int *
	time)
{
    /* System generated locals */
    integer i__1=0, i__2=0;

    /* Local variables */
    doublereal d__=0.0;
    integer i__=0, j=0, ic=0, if__=0;
	unsigned int told=0;
    integer iaux=0;
	unsigned int tnew=0;
    integer iaux1=0;


/*     RESTRICTS RESIDUALS FROM GRID KC-1 TO GRID KC */


/* ===> TRANSFER OF C-TOCHKA DEFECTS */

    /* Parameter adjustments */
    --time;
    --ifg;
    --imaxw;
    --iminw;
    --imax;
    --imin;
    --iw;
    --ja;
    --ia;
    --f;
    --u;
    --a;

    /* Function Body */
	told=clock();
    iaux = ia[imax[*kc - 1] + 1];
    ia[imax[*kc - 1] + 1] = iw[iminw[*kc - 1]];
    iaux1 = iw[imaxw[*kc - 1] + 1];
    iw[imaxw[*kc - 1] + 1] = iaux;
    i__1 = imax[*kc];
    for (ic = imin[*kc]; ic <= i__1; ++ic) {
	if__ = ifg[ic];
	d__ = f[if__];
	i__2 = ia[if__ + 1] - 1;
	for (j = ia[if__]; j <= i__2; ++j) {
	    d__ -= a[j] * u[ja[j]];
/* L80: */
	}
	f[ic] = d__;
/* L100: */
    }

/* ===> TRANSFER OF F-TOCHKA DEFECTS */

    i__1 = imaxw[*kc - 1];
    for (i__ = iminw[*kc - 1]; i__ <= i__1; ++i__) {
	if__ = ifg[i__];
	d__ = f[if__];
	i__2 = ia[if__ + 1] - 1;
	for (j = ia[if__]; j <= i__2; ++j) {
	    d__ -= a[j] * u[ja[j]];
/* L20: */
	}
	i__2 = iw[i__ + 1] - 1;
	for (j = iw[i__]; j <= i__2; ++j) {
	    f[ja[j]] += a[j] * d__;
/* L250: */
	}
/* L200: */
    }
    ia[imax[*kc - 1] + 1] = iaux;
    iw[imaxw[*kc - 1] + 1] = iaux1;
	tnew=clock();
    time[12] = time[12] + tnew - told;
    return 0;
} /* resc_ */


/* ....................................................................... */

/*     PUTZ                                                   SUBROUTINE */

/* ....................................................................... */

/* Subroutine */ integer putz_(integer *k, doublereal *u, integer *imin, integer *
	imax, /*real*/ unsigned int *time)
{
    /* System generated locals */
    integer i__1=0;

    /* Local variables */
    integer i__=0;
	unsigned int told=0, tnew=0;


/*     PUTS ZERO TO U-VALUES OF GRID K */


    /* Parameter adjustments */
    --time;
    --imax;
    --imin;
    --u;

    /* Function Body */
	told=clock();
    i__1 = imax[*k];
    for (i__ = imin[*k]; i__ <= i__1; ++i__) {
	u[i__] = 0.;
/* L10: */
    }
	tnew=clock();
    time[15] = time[15] + tnew - told;
    return 0;
} /* putz_ */


/* ....................................................................... */

/*     FIRST                                                 SUBROUTINE */

/* ....................................................................... */

/* Subroutine */ integer first_(integer *ifirst, doublereal *u, integer *imin, 
	integer *imax, integer *iarr, integer *irow0)
{
    /* System generated locals */
    integer i__1=0;

    /* Local variables */
    integer i__=0;
    doublereal s=0.0;
    extern /* Subroutine */ integer idec_(integer *, integer *, integer *, 
	    integer *);
    integer ifrst=0, ndigit=0;
    extern doublereal random_(doublereal *);


/*     PUTS A FIRST APPROXIMATION TO FINEST GRID */


    /* Parameter adjustments */
    --iarr;
    --imax;
    --imin;
    --u;

    /* Function Body */
    if (*ifirst != 0) {
	idec_(ifirst, &c__3, &ndigit, &iarr[1]);
	ifrst = iarr[2];
    } else {
	ifrst = 3;
    }

    switch (ifrst) {
	case 1:  goto L100;
	case 2:  goto L200;
	case 3:  goto L300;
    }
    return 0;

L100:
    i__1 = imax[1];
    for (i__ = imin[1]; i__ <= i__1; ++i__) {
	u[i__] = 0.;
/* L110: */
    }
    return 0;

L200:
    i__1 = imax[1];
    for (i__ = imin[1]; i__ <= i__1; ++i__) {
	u[i__] = 1.;
/* L210: */
    }
    if (*irow0 < 2) {
	u[imax[1]] = 0.;
    }
    return 0;

L300:
    if (iarr[3] * *ifirst == 0) {
	goto L350;
    }
    s = static_cast<doublereal>(iarr[3]);
    for (i__ = 1; i__ <= 10; ++i__) {
	s *= .1;
	if (s < 1.) {
	    goto L370;
	}
/* L310: */
    }
L350:
    s = .72815;
L370:
    i__1 = imax[1];
    for (i__ = imin[1]; i__ <= i__1; ++i__) {
	u[i__] = random_(&s);
/* L390: */
    }
    if (*irow0 < 2) {
	u[imax[1]] = 0.;
    }
    return 0;
} /* first_ */


/* ....................................................................... */

/*     INJF                                                   SUBROUTINE */

/* ....................................................................... */

/* Subroutine */ integer injf_(integer *kc, doublereal *f, integer *imin, integer 
	*imax, integer *ifg)
{
    /* System generated locals */
    integer i__1=0;

    /* Local variables */
    integer ic=0;


/*     INJECTS F-VALUES FROM GRID KC-1 TO GRID KC */

    /* Parameter adjustments */
    --ifg;
    --imax;
    --imin;
    --f;

    /* Function Body */
    i__1 = imax[*kc];
    for (ic = imin[*kc]; ic <= i__1; ++ic) {
	f[ic] = f[ifg[ic]];
/* L10: */
    }
    return 0;
} /* injf_ */


/* ....................................................................... */

/*     INJU                                                   SUBROUTINE */

/* ....................................................................... */

/* Subroutine */ integer inju_(integer *kc, doublereal *u, integer *imin, integer 
	*imax, integer *ifg)
{
    /* System generated locals */
    integer i__1=0;

    /* Local variables */
    integer ic=0;


/*     INJECTS U-VALUES FROM GRID KC-1 TO GRID KC */

    /* Parameter adjustments */
    --ifg;
    --imax;
    --imin;
    --u;

    /* Function Body */
    i__1 = imax[*kc];
    for (ic = imin[*kc]; ic <= i__1; ++ic) {
	u[ic] = u[ifg[ic]];
/* L10: */
    }
    return 0;
} /* inju_ */


/* ....................................................................... */

/*     NRMU                                                SUBROUTINE */

/* ....................................................................... */

/* Subroutine */ integer nrmu_(integer *k, doublereal *u, integer *imin, integer *
	imax, integer *irow0)
{
    /* System generated locals */
    integer i__1=0;

    /* Local variables */
    integer i__=0;
    doublereal fac=0.0;


/*     NORMALIZES U ON GRID K IF ROWSUM=0 (LAST COMPONENT =0) */


    /* Parameter adjustments */
    --imax;
    --imin;
    --u;

    /* Function Body */
    if (*irow0 > 1) {
	return 0;
    }
    fac = u[imax[*k]];
    i__1 = imax[*k];
    for (i__ = imin[*k]; i__ <= i__1; ++i__) {
	u[i__] -= fac;
/* L10: */
    }
    return 0;
} /* nrmu_ */


/* ....................................................................... */

/*     RESID                                                 SUBROUTINE */

/* ....................................................................... */

/* Subroutine */ integer resid_(integer *k, doublereal *res, doublereal *a, 
	doublereal *u, doublereal *f, integer *ia, integer *ja, integer *iw, 
	integer *imin, integer *imax, integer *iminw)
{
    /* System generated locals */
    integer i__1=0, i__2=0;

    /* Builtin functions */
    //double sqrt(double);

    /* Local variables */
    integer i__=0, j=0;
    doublereal s=0.0;
    integer iaux=0;


/*     COMPUTES L2-NORM OF RESIDUAL ON GRID K */


    /* Parameter adjustments */
    --iminw;
    --imax;
    --imin;
    --iw;
    --ja;
    --ia;
    --f;
    --u;
    --a;

    /* Function Body */
    *res = 0.;
    iaux = ia[imax[*k] + 1];
    ia[imax[*k] + 1] = iw[iminw[*k]];
    i__1 = imax[*k];
    for (i__ = imin[*k]; i__ <= i__1; ++i__) {
	s = f[i__];
	i__2 = ia[i__ + 1] - 1;
	for (j = ia[i__]; j <= i__2; ++j) {
	    s -= a[j] * u[ja[j]];
/* L10: */
	}
	*res += s * s;
/* L20: */
    }
    ia[imax[*k] + 1] = iaux;
    *res = sqrt(*res);
    return 0;
} /* resid_ */


/* ....................................................................... */

/*     CG                                                 SUBROUTINE */

/* ....................................................................... */

/* Subroutine */ integer cg_(integer *k, integer *icgr, integer *iter, doublereal 
	*a, doublereal *u, doublereal *f, integer *ia, integer *ja, integer *
	iw, integer *imin, integer *imax, integer *iminw, integer *m, /*real*/ unsigned int *
	time, integer *ierr, integer *ium)
{
    /* System generated locals */
    integer i__1=0;

    /* Local variables */
    integer i__=0;
    doublereal s2=0.0, alf=0.0, eps=0.0;
    //integer nnu;
	unsigned int told=0, tnew=0;
    extern doublereal cgalf_(integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *), cgeps_(integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     integer *, integer *, integer *, integer *, integer *, integer *,
	     integer *, integer *);
    integer ishift=0;


/*     PERFORMS ONE STEP OF PRECONDITIONED CONJUGATE GRADIENT */


    /* Parameter adjustments */
    --time;
    --iminw;
    --imax;
    --imin;
    --iw;
    --ja;
    --ia;
    --f;
    --u;
    --a;

    /* Function Body */
    if (*icgr == 0) {
	return 0;
    }

/* ===> COMPUTE MOST RECENT MG CORRECTION */

	told=clock();
    //nnu = imax[*k] - imin[*k] + 1;
    ishift = imax[*m] + 1 - imin[*k];
    i__1 = imax[*k];
    for (i__ = imin[*k]; i__ <= i__1; ++i__) {
	u[i__] -= u[i__ + ishift];
/* L10: */
    }
    if (*icgr == 2 && *iter > 1) {
	goto L100;
    }

/* ===> FIRST CG STEP */

    i__1 = imax[*k];
    for (i__ = imin[*k]; i__ <= i__1; ++i__) {
	f[i__ + ishift] = u[i__];
/* L20: */
    }
    goto L200;

/* ===> NEXT CG STEPS (IF ICGR=2 ONLY) */

L100:
    alf = cgalf_(k, &s2, &a[1], &u[1], &f[1], &ia[1], &ja[1], &iw[1], &imin[1]
	    , &imax[1], &iminw[1], m);
    i__1 = imax[*k];
    for (i__ = imin[*k]; i__ <= i__1; ++i__) {
	f[i__ + ishift] = u[i__] + alf * f[i__ + ishift];
/* L110: */
    }
L200:
    eps = cgeps_(k, &s2, &a[1], &u[1], &f[1], &ia[1], &ja[1], &iw[1], &imin[1]
	    , &imax[1], &iminw[1], m, ierr, ium);
    if (*ierr > 0) {
	return 0;
    }
    i__1 = imax[*k];
    for (i__ = imin[*k]; i__ <= i__1; ++i__) {
	u[i__] = u[i__ + ishift] + eps * f[i__ + ishift];
/* L210: */
    }
	tnew=clock();
    time[16] = time[16] + tnew - told;
    return 0;
} /* cg_ */


/* ....................................................................... */

/*     USAVE                                            SUBROUTINE */

/* ....................................................................... */

/* Subroutine */ integer usave_(integer *k, integer *icgr, doublereal *u, integer 
	*imin, integer *imax, integer *ndu, integer *m, /*real*/ unsigned int *time, integer *
	ierr, integer *ium, integer *mdu, integer *ndf, integer *mdf)
{
    /* Format strings */
   // static char fmt_9000[] = "(\002 --- WARNG IN USAVE: NO CG BECAUSE OF STO"
	//    "RAGE ---\002/\002     REQUIRED: NDU =\002,i9/\002               "
	 //   "NDF =\002,i9)";

    /* System generated locals */
    integer i__1=0, i__2=0;

    /* Builtin functions */
    //integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);

    /* Local variables */
    integer i__=0;
	unsigned int told=0, tnew=0;
    integer ishift=0;

    /* Fortran I/O blocks */
    //static cilist io___371 = { 0, 0, 0, fmt_9000, 0 };



/*     MAKES A BACK-UP OF THE CURRENT APPROXIMATION ON LEVEL K */
/*     IF ICGR.NE.0. */


    /* Parameter adjustments */
    --time;
    --imax;
    --imin;
    --u;

    /* Function Body */
    if (*icgr == 0) {
	return 0;
    }
    ishift = imax[*m] + 1 - imin[*k];
/* Computing MAX */
    i__1 = *mdu, i__2 = imax[*k] + ishift;
    *mdu = myi_max(i__1,i__2);
/* Computing MAX */
    i__1 = *mdf, i__2 = imax[*k] + ishift;
    *mdf = myi_max(i__1,i__2);
    if (*mdu <= *ndu && *mdf <= *ndf) {
	goto L10;
    }
   

	printf("( --- WARNG IN USAVE: NO CG BECAUSE OF STORAGE\n");
#if doubleintprecision == 1
	printf(" ---    REQUIRED: NDU =,%lld,      \n", mdu[0]);
	printf("NDF =,%lld)\n", mdf[0]);
#else
	printf(" ---    REQUIRED: NDU =,%d,      \n", mdu[0]);
	printf("NDF =,%d)\n", mdf[0]);
#endif
	
	    
	   

    *ierr = -4;
    *icgr = 0;
    return 0;

L10:
	told=clock();
    i__1 = imax[*k];
    for (i__ = imin[*k]; i__ <= i__1; ++i__) {
	u[i__ + ishift] = u[i__];
/* L20: */
    }
	tnew=clock();
    time[16] = time[16] + tnew - told;
    return 0;

} /* usave_ */


/* ....................................................................... */

/*     CGEPS                                            FUNCTION */

/* ....................................................................... */

doublereal cgeps_(integer *k, doublereal *s2, doublereal *a, doublereal *u, 
	doublereal *f, integer *ia, integer *ja, integer *iw, integer *imin, 
	integer *imax, integer *iminw, integer *m, integer *ierr, integer *
	ium)
{
    /* Format strings */
   // static char fmt_9000[] = "(\002 *** ERROR IN CGEPS: CG CORRECTION NOT DE"
	//    "FINED ***\002)";

    /* System generated locals */
    integer i__1=0, i__2=0;
    doublereal ret_val=0.0; // инициализация была добавлена позже мной её не было изначально.

    /* Builtin functions */
   // integer s_wsfe(cilist *), e_wsfe(void);

    /* Local variables */
    integer i__=0, j=0;
    doublereal s1=0.0, sp=0.0, sr=0.0;
    integer iaux=0, ishift=0;

    /* Fortran I/O blocks */
    //static cilist io___382 = { 0, 0, 0, fmt_9000, 0 };




    /* Parameter adjustments */
    --iminw;
    --imax;
    --imin;
    --iw;
    --ja;
    --ia;
    --f;
    --u;
    --a;

    /* Function Body */
    ishift = imax[*m] + 1 - imin[*k];
    s1 = 0.;
    *s2 = 0.;
    iaux = ia[imax[*k] + 1];
    ia[imax[*k] + 1] = iw[iminw[*k]];
    i__1 = imax[*k];
    for (i__ = imin[*k]; i__ <= i__1; ++i__) {
	sr = f[i__];
	sp = 0.;
	i__2 = ia[i__ + 1] - 1;
	for (j = ia[i__]; j <= i__2; ++j) {
	    sr -= a[j] * u[ja[j] + ishift];
	    sp += a[j] * f[ja[j] + ishift];
/* L40: */
	}
	s1 += sr * f[i__ + ishift];
	*s2 += sp * f[i__ + ishift];
/* L50: */
    }
    ia[imax[*k] + 1] = iaux;
    if (*s2 == 0.) {
	goto L100;
    }
    ret_val = s1 / *s2;
    return ret_val;

/* ===> ERROR EXIT */

L100:
    //io___382.ciunit = *ium;
    //s_wsfe(&io___382);
    //e_wsfe();
	printf("( *** ERROR IN CGEPS: CG CORRECTION NOT DEFINED ***)\n");
    *ierr = 31;
    return ret_val;
} /* cgeps_ */


/* ....................................................................... */

/*     CGALF                                            FUNCTION */

/* ....................................................................... */

doublereal cgalf_(integer *k, doublereal *s2, doublereal *a, doublereal *u, 
	doublereal *f, integer *ia, integer *ja, integer *iw, integer *imin, 
	integer *imax, integer *iminw, integer *m)
{
    /* System generated locals */
    integer i__1=0, i__2=0;
    doublereal ret_val=0.0;

    /* Local variables */
    integer i__=0, j=0;
    doublereal s1=0.0, sr=0.0;
    integer iaux=0, ishift=0;



    /* Parameter adjustments */
    --iminw;
    --imax;
    --imin;
    --iw;
    --ja;
    --ia;
    --f;
    --u;
    --a;

    /* Function Body */
    ishift = imax[*m] + 1 - imin[*k];
    s1 = 0.;
    iaux = ia[imax[*k] + 1];
    ia[imax[*k] + 1] = iw[iminw[*k]];
    i__1 = imax[*k];
    for (i__ = imin[*k]; i__ <= i__1; ++i__) {
	sr = 0.;
	i__2 = ia[i__ + 1] - 1;
	for (j = ia[i__]; j <= i__2; ++j) {
	    sr += a[j] * u[ja[j]];
/* L40: */
	}
	s1 += sr * f[i__ + ishift];
/* L50: */
    }
    ret_val = -s1 / *s2;
    ia[imax[*k] + 1] = iaux;
    return ret_val;
} /* cgalf_ */

/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */

/*     OUTPUT OF STATISTICAL INFORMATION ON CP-TIMES AND DIMENSIONING */

/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
/* ....................................................................... */

/*     WRKCNT                                           SUBROUTINE */

/* ....................................................................... */

/* Subroutine */ integer wrkcnt_(integer *iout, integer *ia, integer *iw, integer 
	*imin, integer *imax, integer *iminw, integer *levels, /*real*/ unsigned int *time, 
	integer *ncyc0, integer *iup, integer *mda, integer *mdia, integer *
	mdja, integer *mdu, integer *mdf, integer *mdig, doublereal *res0, 
	doublereal *res)
{
    /* Format strings */
    //static char fmt_9110[] = "(//\002 **************** CONVERGENCE *********"
	//    "********\002/\002 L2-NORM OF RESIDUAL BEFORE CYCLING =\002,d10.3/"
	//    "\002 L2-NORM OF RESIDUAL AFTER  CYCLING =\002,d10.3/\002 CONVERG"
	//    "ENCE FACTOR                 =\002,d10.3)";
   // static char fmt_9120[] = "(\002 CONVERGENCE FACTOR PER CYCLE       =\002"
	//    ",d10.3)";
    //static char fmt_9000[] = "(//\002 ************** WORK COUNT ************"
	//    "***\002/)";
    //static char fmt_9020[] = "(\002   PREP       SEC       SOL      SEC/CY"
	  //  "CLE\002)";
    //static char fmt_9030[] = "(\002 --------------------------------------"
	//    "---\002)";
   // static char fmt_9040[] = "(\002 1 RWSRT   \002,f7.2,\002   11 INTADD  "
   //	    " \002,f7.2/\002 2 PRE-COL \002,f7.2,\002   12 RESCAL   \002,f7.2/"
   //	    "\002 3 CHK-COL \002,f7.2,\002   13 RELAX    \002,f7.2/\002 4 INT"
   //	    "ERPOL\002,f7.2,\002   14 V-*      \002,f7.2/\002 5 RESTRICT\002,"
//	    "f7.2,\002   15 OTHERS   \002,f7.2/\002 6 OPDFN   \002,f7.2,\002 "
	//    "  16 CONJ-GRAD\002,f7.2/\002 7 TRUNC   \002,f7.2,\002   17 YALE-"
	 //   "SMP \002,f7.2/\002 8 OTHERS  \002,f7.2,\002   18 ------   \002,f"
	 //   "7.2)";
    //static char fmt_9050[] = "(\002   SUM     \002,f7.2,\002      SUM     "
	//    " \002,f7.2)";
    //static char fmt_9100[] = "(//\002 ********* SPACE REQUIREMENTS ********"
	//    "*\002//\002 VECTOR      NEEDED      THEOR. MINIMUM\002/\002 ----"
	//    "----------------------------------\002/\002   A \002,i14,2x,i16"
	//    "/\002   JA\002,i14,2x,i16/\002   IA\002,i14,2x,i16/\002   U \002"
	//    ",i14,2x,i16/\002   F \002,i14,2x,i16/\002   IG\002,i14,2x,i16"
	//    "/\002 --------------------------------------\002/)";
    //static char fmt_9080[] = "(/\002 ******************* COMPLEXITIES ******"
	//    "**************\002/\002 SPACE OCCUPIED BY ALL OPERATORS / SPACE "
	//    "OF OPERATOR  \002/\002 ON THE FINEST GRID   = \002,f8.2,\002   ("
	//    "A-COMPLEXITY)     \002/\002 TOTAL NUMBER OF GRID POINTS / NUMBER"
	//    " OF POINTS IN    \002/\002 THE  FINEST  GRID    = \002,f8.2,\002"
	//    "   (O-COMPLEXITY)     \002/\002 TOTAL SPACE USED BY AMG1R5 / SPA"
	//    "CE OCCUPIED BY USER- \002/\002 DEFINED  PROBLEM     = \002,f8.2"
	//    ",\002   (S-COMPLEXITY)     \002/\002 SPACE USED DURING SOLUTION "
	//    "PHASE / SPACE OCCUPIED BY \002/\002 USER-DEFINED PROBLEM = \002,"
	//    "f8.2/\002 ****************************************************"
	//    "*\002)";

    /* System generated locals */
    integer i__1=0;
    doublereal d__1=0.0;

    /* Builtin functions */
   // integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);
	doublereal pow_dd(doublereal *, doublereal *);

    /* Local variables */
    integer i__=0, k=0;
    //static real t[10];
	unsigned int t[10] = { 0 };
	integer nnu=0;
   // static real sum1, sum2;
	unsigned int sum1=0, sum2=0;
    doublereal cfac=0.0, cfpc=0.0;
    integer mdta=0, mdtf=0, mdtu=0, idima=0, mdtia=0, mdtja=0, mdtig=0;
    doublereal acmplx=0.0, ocmplx=0.0, scmplx=0.0, tcmplx=0.0;

    /* Fortran I/O blocks */
    //static cilist io___390 = { 0, 0, 0, fmt_9110, 0 };
   // static cilist io___392 = { 0, 0, 0, fmt_9120, 0 };
   // static cilist io___393 = { 0, 0, 0, fmt_9000, 0 };
   // static cilist io___399 = { 0, 0, 0, fmt_9020, 0 };
   // static cilist io___400 = { 0, 0, 0, fmt_9030, 0 };
    //static cilist io___401 = { 0, 0, 0, fmt_9040, 0 };
   // static cilist io___402 = { 0, 0, 0, fmt_9030, 0 };
   // static cilist io___403 = { 0, 0, 0, fmt_9050, 0 };
    //static cilist io___404 = { 0, 0, 0, fmt_9030, 0 };
    //static cilist io___413 = { 0, 0, 0, fmt_9100, 0 };
    //static cilist io___418 = { 0, 0, 0, fmt_9080, 0 };



/*     RESIDUALS / CP-TIMES / COMPLEXITY / DIMENSIONING */


    /* Parameter adjustments */
    --time;
    --iminw;
    --imax;
    --imin;
    --iw;
    --ia;

    /* Function Body */
    if (*iout < 1) {
	return 0;
    }

/* ===> RESIDUALS / CONVERGENCE */

    if (*ncyc0 > 0) {
	cfac = *res / (*res0 + 1e-40);
	
	if (yes_print_amg) {
	    printf("( **************** CONVERGENCE *********\n");
	    printf("******** L2-NORM OF RESIDUAL BEFORE CYCLING =,%1.4f\n",res0[0]);
	    printf(" L2-NORM OF RESIDUAL AFTER  CYCLING =,%1.4f, CONVERG\n",res[0]);
	    printf("ENCE FACTOR                 =,%1.4f)\n",cfac);
	}

	d__1 = 1. / static_cast<doublereal> (*ncyc0);
	cfpc = pow_dd(&cfac, &d__1);
	
	if (yes_print_amg) {
	    printf("( CONVERGENCE FACTOR PER CYCLE       =%1.4f\n", cfpc);
	}
	 
    }
    if (*iout <= 1) {
	return 0;
    }
   
	if (yes_print_amg) {
	    printf("( ************** WORK COUNT ***************)\n");
	}

    nnu = imax[1];

/* ===> COMPUTING TIMES */

    //sum1 = 0.f;
    //sum2 = 0.f;
	sum1=0;
	sum2=0;
    for (i__ = 1; i__ <= 10; ++i__) {
	//t[i__ - 1] = 0.f;
		t[i__ - 1]=0;
	if (*ncyc0 > 0) {
	   // t[i__ - 1] = time[i__ + 10] / (real) (*ncyc0);
		 t[i__ - 1] = (unsigned int)( time[i__ + 10] / (*ncyc0));
	}
	sum1 += time[i__];
	sum2 += t[i__ - 1];
/* L10: */
    }

   
	if (yes_print_amg) {
	    printf("(   PREP       SEC       SOL      SEC/CYCLE)\n");
	}

   
	if (yes_print_amg) {
	    printf("( -----------------------------------------)\n");
	}

   
if (yes_print_amg) {

	printf("( 1 RWSRT   ,%d,   11 INTADD  \n", time[1]);
	printf(" ,%d 2 PRE-COL ,%d,   12 RESCAL   ,%d\n", t[1 - 1], time[2], t[2 - 1]);
	printf(" 3 CHK-COL ,%d,   13 RELAX    ,%d 4 INTERPOL\n", time[3], t[3 - 1]);
	printf(",%d,   14 V-*      ,%d 5 RESTRICT,\n", time[4], t[4 - 1]);
	printf("%d,   15 OTHERS   ,%d 6 OPDFN   ,%d, \n", time[5], t[5 - 1], time[6]);
	printf("  16 CONJ-GRAD ,%d 7 TRUNC   ,%d,   17 YALE-\n", t[6 - 1], time[7]);
	printf("SMP ,%d 8 OTHERS  ,%d,   18 ------   ,%d)\n", t[7 - 1], time[8], t[8 - 1]);





	printf("( -----------------------------------------)\n");



	printf("(   SUM     ,%d,      SUM      ,%d)\n", sum1, sum2);



	printf("( -----------------------------------------)\n");

	

}

/* ===> SPACE OCCUPIED BY OPERATORS A(1) - A(LEVELS) */

    idima = 0;
    i__1 = *levels;
    for (k = 1; k <= i__1; ++k) {
	idima = idima + iw[iminw[k]] - ia[imin[k]];
/* L20: */
    }

/* ===> THEORETICAL MINIMAL SPACE REQUIREMENTS */

    if (*levels < 2) {
	return 0;
    }
    mdta = iw[iminw[*levels]] - 1;
    mdtja = mdta;
    mdtia = imax[*levels] + 1;
    mdtu = imax[*levels];
    mdtf = mdtu;
    mdtig = (mdtu << 1) + nnu;
   

	if (yes_print_amg) {
#if doubleintprecision == 1
		printf("( ********* SPACE REQUIREMENTS ********\n");
		printf("* VECTOR      NEEDED      THEOR. MINIMUM ----\n");
		printf("----------------------------------   A ,%lld,2x,%lld\n", mda[0], mdta);
		printf("/   JA,%lld,2x,%lld   IA,%lld,2x,%lld   U \n", mdja[0], mdtja, mdia[0], mdtia);
		printf(",%lld,2x,%lld   F ,%lld,2x,%lld   IG,%lld,2x,%lld\n", mdu[0], mdtu, mdf[0], mdtf, mdig[0], mdtig);
		printf("/ --------------------------------------)\n");
#else
		printf("( ********* SPACE REQUIREMENTS ********\n");
		printf("* VECTOR      NEEDED      THEOR. MINIMUM ----\n");
		printf("----------------------------------   A ,%d,2x,%d\n", mda[0], mdta);
		printf("/   JA,%d,2x,%d   IA,%d,2x,%d   U \n", mdja[0], mdtja, mdia[0], mdtia);
		printf(",%d,2x,%d   F ,%d,2x,%d   IG,%d,2x,%d\n", mdu[0], mdtu, mdf[0], mdtf, mdig[0], mdtig);
		printf("/ --------------------------------------)\n");
#endif
	 
	}

/* ===> COMPLEXITIES */

    scmplx = static_cast<doublereal> (((*mda + *mdu + *mdf) << 1) + *mdja + *mdia + *mdig) 
	    / static_cast<doublereal> (nnu * 5 + 1 + (iw[iminw[1]] - ia[imin[1]]) * 3);
    tcmplx = static_cast<doublereal> (((mdta + mdtu + mdtf) << 1) + mdtja + mdtia + mdtig) // скобки поставлены в соответствии с приоритетом.
	    / static_cast<doublereal> (nnu * 5 + 1 + (iw[iminw[1]] - ia[imin[1]]) * 3);

    acmplx = static_cast<doublereal> (idima) / static_cast<doublereal> (iw[1] - 1);
    ocmplx = static_cast<doublereal> (mdtu) / static_cast<doublereal> (nnu);
   
	if (yes_print_amg) {

     	printf("( ******************* COMPLEXITIES ******\n");
	    printf("************** SPACE OCCUPIED BY ALL OPERATORS / SPACE \n");
	    printf("OF OPERATOR   ON THE FINEST GRID   = ,%1.4f,   (\n",acmplx);
	    printf("A-COMPLEXITY)      TOTAL NUMBER OF GRID POINTS / NUMBER\n");
	    printf(" OF POINTS IN     THE  FINEST  GRID    = ,%1.4f,\n",ocmplx);
	    printf("   (O-COMPLEXITY)      TOTAL SPACE USED BY AMG1R5 / SPA\n");
	    printf("CE OCCUPIED BY USER-  DEFINED  PROBLEM     = ,%1.4f\n",scmplx);
	    printf(",   (S-COMPLEXITY)      SPACE USED DURING SOLUTION \n");
	    printf("PHASE / SPACE OCCUPIED BY  USER-DEFINED PROBLEM = %1.4f,\n",tcmplx);
	    printf(" *****************************************************)\n");

	}

    return 0;

} /* wrkcnt_ */

/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */

/*     AMG1R5 AUXILIARY PROGRAMS */

/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
/* ....................................................................... */

/*     IDEC                                                  SUBROUTINE */

/* ....................................................................... */

/* Subroutine */ integer idec_(integer *into, integer *nnum, integer *ndigit, 
	integer *iarr)
{
    /* Initialized data */

    doublereal eps = .5;

    /* System generated locals */
    integer i__1=0;
    doublereal d__1=0.0;

    /* Builtin functions */
	doublereal d_lg10(doublereal *);
    integer pow_ii(integer *, integer *);

    /* Local variables */
    integer i__=0, iq=0, nrest=0;


/*     DECOMPOSE NON-NEGATIVE INTEGER INTO INTO NNUM INTEGERS */

/*     INPUT:  INTO   - INTEGER (0.LE. INTO .LE.999999999) */
/*             NNUM   - INTEGER (1.LE. NNUM .LE.9); NUMBER OF INTEGERS */
/*                      TO BE RETURNED ON ARRAY IARR (SEE BELOW) */

/*     OUTPUT: NDIGIT - INTEGER; NUMBER OF DIGITS OF INTO */
/*             IARR   - INTEGER-ARRAY OF LENGTH 10: */
/*                      IARR(1)        = FIRST      DIGIT OF INTO, */
/*                      IARR(2)        = SECOND     DIGIT OF INTO, ..... */
/*                      IARR(NNUM-1)   = (NNUM-1)ST DIGIT OF INTO, */
/*                      IARR(NNUM)     = REST OF INTO */
/*                      IF NNUM > NDIGIT, THE CORRESPONDING COMPONENTS */
/*                      OF IARR ARE PUT TO ZERO. */

/*     WARNING: BE SURE THAT YOUT COMPUTER CAN STORE NNUM DIGITS ON AN */
/*              INTEGER VARIABLE. */

    /* Parameter adjustments */
    --iarr;

    /* Function Body */

    nrest = *into;
    d__1 = eps + static_cast<doublereal> (*into);
    *ndigit = static_cast<integer>(d_lg10(&d__1)) + 1;
    if (*nnum >= *ndigit) {
	goto L20;
    }
    i__1 = *ndigit - *nnum + 1;
    iq = pow_ii(&c__10, &i__1);
    iarr[*nnum] = nrest - nrest / iq * iq;
    nrest /= iq;
    for (i__ = *nnum - 1; i__ >= 1; --i__) {
	iarr[i__] = nrest - nrest / 10 * 10;
	nrest /= 10;
/* L10: */
    }
    return 0;

L20:
    for (i__ = *ndigit; i__ >= 1; --i__) {
	iarr[i__] = nrest - nrest / 10 * 10;
	nrest /= 10;
/* L30: */
    }
    i__1 = *nnum;
    for (i__ = *ndigit + 1; i__ <= i__1; ++i__) {
	iarr[i__] = 0;
/* L40: */
    }
    return 0;
} /* idec_ */


/* ....................................................................... */

/*     RANDOM                                                FUNCTION */

/* ....................................................................... */

doublereal random_(doublereal *s)
{
    /* System generated locals */
    doublereal ret_val=0.0;

    /* Builtin functions */
   // double exp(double);


/*     FUNCTION TO CREATE "RANDOM" SEQUENCE OF NUMBERS BETWEEN 0 AND 0.1 */

/*     INPUT:   S      - NUMBER BETWEEN 0 AND 0.1 */

/*     OUTPUT:  RANDOM - NUMBER BETWEEN 0 AND 0.1 */
/*              S      - S=RANDOM */

    ret_val = exp(*s) * 100.;
    ret_val -= static_cast<doublereal> (static_cast<integer> (ret_val));
    *s = ret_val;
    return ret_val;
} /* random_ */

/* *********************************************************************** */
/* $LARGE */
/* $NOFLOATCALLS */
/*                            APPENDIX 1                          1/15/81 */

/*        SUBROUTINES FOR SOLVING SPARSE NONSYMMETRIC SYSTEMS */
/*        OF LINEAR EQUATIONS  (UNCOMPRESSED POINTER STORAGE) */

/*        REAL*8 VERSION. NOTE: THE ORIGINAL SUBROUTINES */

/*            NDRV, NSF, NNF, NNS AND NNT */

/*        HAVE BEEN RENAMED TO */

/*            YALE8, NSF8, NNF8, NNS8 AND NNT8, RESPECTIVELY. */

/* *** SUBROUTINE YALE8 (OLD NAME: NDRV) */
/* *** DRIVER FOR SUBROUTINES FOR SOLVING SPARSE NONSYMMETRIC SYSTEMS OF */
/*       LINEAR EQUATIONS (UNCOMPRESSED POINTER STORAGE) */

/*       SUBROUTINE  NDRV  (= OLD NAME) */
/*       SUBROUTINE  YALE8 (= NEW NAME) */
/* Subroutine */ integer ndrv_(integer *n, integer *r__, integer *c__, integer *
	ic, integer *ia, integer *ja, doublereal *a, doublereal *b, 
	doublereal *z__, integer *nsp, /*integer*/ doublereal *isp, doublereal *rsp, integer 
	*esp, integer *path, integer *flag__)
{
    /* Initialized data */

    integer lratio = 2;

    /* System generated locals */
    integer i__1=0;

    /* Local variables */
    integer d__=0, j=0, l=0, q=0, u=0, il=0, im=0, jl=0, iu=0, ju=0, max__=0, tmp=0, row=0;
    extern /* Subroutine */ integer nnf8_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, /*integer*/ doublereal *, /*integer*/ doublereal *, doublereal *, integer *, 
	    doublereal *, /*integer*/ doublereal *, /*integer*/ doublereal *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *), nsf8_(integer *, integer *
	    , integer *, integer *, integer *, /*integer*/ doublereal *, /*integer*/ doublereal *, integer *
	    , /*integer*/ doublereal *, /*integer*/ doublereal *, integer *, /*integer*/ doublereal *, /*integer*/ doublereal *, integer *
	    ), nns8_(integer *, integer *, integer *, /*integer*/ doublereal *, /*integer*/ doublereal *, 
	    doublereal *, doublereal *, /*integer*/ doublereal *, /*integer*/ doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *), nnt8_(integer *, 
	    integer *, integer *, /*integer*/ doublereal *, /*integer*/ doublereal *, doublereal *, 
	    doublereal *, /*integer*/ doublereal *, /*integer*/ doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *);
    integer lmax=0, umax=0, jlmax=0, jumax=0, jutmp=0;


/*    PARAMETERS */
/*    CLASS ABBREVIATIONS ARE -- */
/*       N - INTEGER VARIABLE */
/*       F - REAL VARIABLE */
/*       V - SUPPLIES A VALUE TO THE DRIVER */
/*       R - RETURNS A RESULT FROM THE DRIVER */
/*       I - USED INTERNALLY BY THE DRIVER */
/*       A - ARRAY */

/* CLASS   PARAMETER */
/* ------+---------- */

/*         THE NONZERO ENTRIES OF THE COEFFICIENT MATRIX M ARE STORED */
/*    ROW-BY-ROW IN THE ARRAY A.  TO IDENTIFY THE INDIVIDUAL NONZERO */
/*    ENTRIES IN EACH ROW, WE NEED TO KNOW IN WHICH COLUMN EACH ENTRY */
/*    LIES.  THE COLUMN INDICES WHICH CORRESPOND TO THE NONZERO ENTRIES */
/*    OF M ARE STORED IN THE ARRAY JA;  I.E., IF  A(K) = M(I,J),  THEN */
/*    JA(K) = J.  IN ADDITION, WE NEED TO KNOW WHERE EACH ROW STARTS AND */
/*    HOW LONG IT IS.  THE INDEX POSITIONS IN JA AND A WHERE THE ROWS OF */
/*    M BEGIN ARE STORED IN THE ARRAY IA;  I.E., IF M(I,J) IS THE FIRST */
/*    NONZERO ENTRY (STORED) IN THE I-TH ROW AND A(K) = M(I,J),  THEN */
/*    IA(I) = K.  MOREOVER, THE INDEX IN JA AND A OF THE FIRST LOCATION */
/*    FOLLOWING THE LAST ELEMENT IN THE LAST ROW IS STORED IN IA(N+1). */
/*    THUS, THE NUMBER OF ENTRIES IN THE I-TH ROW IS GIVEN BY */
/*    IA(I+1) - IA(I),  THE NONZERO ENTRIES OF THE I-TH ROW ARE STORED */
/*    CONSECUTIVELY IN */
/*            A(IA(I)),  A(IA(I)+1),  ..., A(IA(I+1)-1), */
/*    AND THE CORRESPONDING COLUMN INDICES ARE STORED CONSECUTIVELY IN */
/*            JA(IA(I)), JA(IA(I)+1), ..., JA(IA(I+1)-1). */
/*    FOR EXAMPLE, THE 5 BY 5 MATRIX */
/*                ( 1. 0. 2. 0. 0.) */
/*                ( 0. 3. 0. 0. 0.) */
/*            M = ( 0. 4. 5. 6. 0.) */
/*                ( 0. 0. 0. 7. 0.) */
/*                ( 0. 0. 0. 8. 9.) */
/*    WOULD BE STORED AS */
/*                 1  2  3  4  5  6  7  8  9 */
/*            ---+-------------------------- */
/*            IA   1  3  4  7  8 10 */
/*            JA   1  3  2  2  3  4  4  4  5 */
/*             A   1. 2. 3. 4. 5. 6. 7. 8. 9.         . */

/* NV      N     - NUMBER OF VARIABLES/EQUATIONS. */
/* FVA     A     - NONZERO ENTRIES OF THE COEFFICIENT MATRIX M, STORED */
/*                   BY ROWS. */
/*                   SIZE = NUMBER OF NONZERO ENTRIES IN M. */
/* NVA     IA    - POINTERS TO DELIMIT THE ROWS IN A. */
/*                   SIZE = N+1. */
/* NVA     JA    - COLUMN NUMBERS CORRESPONDING TO THE ELEMENTS OF A. */
/*                   SIZE = SIZE OF A. */
/* FVA     B     - RIGHT-HAND SIDE B;  B AND Z CAN THE SAME ARRAY. */
/*                   SIZE = N. */
/* FRA     Z     - SOLUTION X;  B AND Z CAN BE THE SAME ARRAY. */
/*                   SIZE = N. */

/*         THE ROWS AND COLUMNS OF THE ORIGINAL MATRIX M CAN BE */
/*    REORDERED (E.G., TO REDUCE FILLIN OR ENSURE NUMERICAL STABILITY) */
/*    BEFORE CALLING THE DRIVER.  IF NO REORDERING IS DONE, THEN SET */
/*    R(I) = C(I) = IC(I) = I  FOR I=1,...,N.  THE SOLUTION Z IS RETURNED */
/*    IN THE ORIGINAL ORDER. */

/* NVA     R     - ORDERING OF THE ROWS OF M. */
/*                   SIZE = N. */
/* NVA     C     - ORDERING OF THE COLUMNS OF M. */
/*                   SIZE = N. */
/* NVA     IC    - INVERSE OF THE ORDERING OF THE COLUMNS OF M;  I.E., */
/*                   IC(C(I)) = I  FOR I=1,...,N. */
/*                   SIZE = N. */

/*         THE SOLUTION OF THE SYSTEM OF LINEAR EQUATIONS IS DIVIDED INTO */
/*    THREE STAGES -- */
/*      NSF -- THE MATRIX M IS PROCESSED SYMBOLICALLY TO DETERMINE WHERE */
/*              FILLIN WILL OCCUR DURING THE NUMERIC FACTORIZATION. */
/*      NNF -- THE MATRIX M IS FACTORED NUMERICALLY INTO THE PRODUCT LDU */
/*              OF A UNIT LOWER TRIANGULAR MATRIX L, A DIAGONAL MATRIX D, */
/*              AND A UNIT UPPER TRIANGULAR MATRIX U, AND THE SYSTEM */
/*              MX = B  IS SOLVED. */
/*      NNS -- THE LINEAR SYSTEM  MX = B  IS SOLVED USING THE LDU */
/*  OR          FACTORIZATION FROM NNF. */
/*      NNT -- THE TRANSPOSED LINEAR SYSTEM  MT X = B  IS SOLVED USING */
/*              THE LDU FACTORIZATION FROM NNF. */
/*    FOR SEVERAL SYSTEMS WHOSE COEFFICIENT MATRICES HAVE THE SAME */
/*    NONZERO STRUCTURE, NSF NEED BE DONE ONLY ONCE (FOR THE FIRST */
/*    SYSTEM);  THEN NNF IS DONE ONCE FOR EACH ADDITIONAL SYSTEM.  FOR */
/*    SEVERAL SYSTEMS WITH THE SAME COEFFICIENT MATRIX, NSF AND NNF NEED */
/*    BE DONE ONLY ONCE (FOR THE FIRST SYSTEM);  THEN NNS OR NNT IS DONE */
/*    ONCE FOR EACH ADDITIONAL RIGHT-HAND SIDE. */

/* NV      PATH  - PATH SPECIFICATION;  VALUES AND THEIR MEANINGS ARE -- */
/*                   1  PERFORM NSF AND NNF. */
/*                   2  PERFORM NNF ONLY  (NSF IS ASSUMED TO HAVE BEEN */
/*                       DONE IN A MANNER COMPATIBLE WITH THE STORAGE */
/*                       ALLOCATION USED IN THE DRIVER). */
/*                   3  PERFORM NNS ONLY  (NSF AND NNF ARE ASSUMED TO */
/*                       HAVE BEEN DONE IN A MANNER COMPATIBLE WITH THE */
/*                       STORAGE ALLOCATION USED IN THE DRIVER). */
/*                   4  PERFORM NNT ONLY  (NSF AND NNF ARE ASSUMED TO */
/*                       HAVE BEEN DONE IN A MANNER COMPATIBLE WITH THE */
/*                       STORAGE ALLOCATION USED IN THE DRIVER). */
/*                   5  PERFORM NSF ONLY. */

/*         VARIOUS ERRORS ARE DETECTED BY THE DRIVER AND THE INDIVIDUAL */
/*    SUBROUTINES. */

/* NR      FLAG  - ERROR FLAG;  VALUES AND THEIR MEANINGS ARE -- */
/*                     0     NO ERRORS DETECTED */
/*                     N+K   NULL ROW IN A  --  ROW = K */
/*                    2N+K   DUPLICATE ENTRY IN A  --  ROW = K */
/*                    3N+K   INSUFFICIENT STORAGE IN NSF  --  ROW = K */
/*                    4N+1   INSUFFICIENT STORAGE IN NNF */
/*                    5N+K   NULL PIVOT  --  ROW = K */
/*                    6N+K   INSUFFICIENT STORAGE IN NSF  --  ROW = K */
/*                    7N+1   INSUFFICIENT STORAGE IN NNF */
/*                    8N+K   ZERO PIVOT  --  ROW = K */
/*                   10N+1   INSUFFICIENT STORAGE IN NDRV */
/*                   11N+1   ILLEGAL PATH SPECIFICATION */

/*         WORKING STORAGE IS NEEDED FOR THE FACTORED FORM OF THE MATRIX */
/*    M PLUS VARIOUS TEMPORARY VECTORS.  THE ARRAYS ISP AND RSP SHOULD BE */
/*    EQUIVALENCED;  INTEGER STORAGE IS ALLOCATED FROM THE BEGINNING OF */
/*    ISP AND REAL STORAGE FROM THE END OF RSP. */

/* NV      NSP   - DECLARED DIMENSION OF RSP;  NSP GENERALLY MUST */
/*                   BE LARGER THAN  5N+3 + 2K  (WHERE  K = (NUMBER OF */
/*                   NONZERO ENTRIES IN M)). */
/* NVIRA   ISP   - INTEGER WORKING STORAGE DIVIDED UP INTO VARIOUS ARRAYS */
/*                   NEEDED BY THE SUBROUTINES;  ISP AND RSP SHOULD BE */
/*                   EQUIVALENCED. */
/*                   SIZE = LRATIO*NSP */
/* FVIRA   RSP   - REAL WORKING STORAGE DIVIDED UP INTO VARIOUS ARRAYS */
/*                   NEEDED BY THE SUBROUTINES;  ISP AND RSP SHOULD BE */
/*                   EQUIVALENCED. */
/*                   SIZE = NSP. */
/* NR      ESP   - IF SUFFICIENT STORAGE WAS AVAILABLE TO PERFORM THE */
/*                   SYMBOLIC FACTORIZATION (NSF), THEN ESP IS SET TO THE */
/*                   AMOUNT OF EXCESS STORAGE PROVIDED (NEGATIVE IF */
/*                   INSUFFICIENT STORAGE WAS AVAILABLE TO PERFORM THE */
/*                   NUMERIC FACTORIZATION (NNF)). */


/*  CONVERSION TO DOUBLE PRECISION */

/*    TO CONVERT THESE ROUTINES FOR DOUBLE PRECISION ARRAYS, SIMPLY USE */
/*    THE DOUBLE PRECISION DECLARATIONS IN PLACE OF THE REAL DECLARATIONS */
/*    IN EACH SUBPROGRAM;  IN ADDITION, THE DATA VALUE OF THE INTEGER */
/*    VARIABLE LRATIO MUST BE SET AS INDICATED IN SUBROUTINE NDRV */

/*       REAL  A(1),  B(1),  Z(1),  RSP(1) */

/*  SET LRATIO EQUAL TO THE RATIO BETWEEN THE LENGTH OF FLOATING TOCHKA */
/*  AND INTEGER ARRAY DATA;  E. G., LRATIO = 1 FOR (REAL, INTEGER), */
/*  LRATIO = 2 FOR (DOUBLE PRECISION, INTEGER) */

    /* Parameter adjustments */
    --rsp;
    --isp;
    --z__;
    --b;
    --a;
    --ja;
    --ia;
    --ic;
    --c__;
    --r__;

    /* Function Body */

    if (*path < 1 || 5 < *path) {
	goto L111;
    }
/*  ******  INITIALIZE AND DIVIDE UP TEMPORARY STORAGE  ***************** */
    il = 1;
    iu = il + *n + 1;
    jl = iu + *n + 1;

/*  ******  CALL NSF IF FLAG IS SET  ************************************ */
    if ((*path - 1) * (*path - 5) != 0) {
	goto L2;
    }
    max__ = lratio * *nsp + 1 - jl - (*n + 1) - *n;
    jlmax = max__ / 2;
    q = jl + jlmax;
    im = q + (*n + 1);
    jutmp = im + *n;
    jumax = lratio * *nsp + 1 - jutmp;
    *esp = max__ / lratio;
    if (jlmax <= 0 || jumax <= 0) {
	goto L110;
    }
    nsf8_(n, &r__[1], &ic[1], &ia[1], &ja[1], &isp[il], &isp[jl], &jlmax, &
	    isp[iu], &isp[jutmp], &jumax, &isp[q], &isp[im], flag__);
    if (*flag__ != 0) {
	goto L100;
    }
/*  ******  MOVE JU NEXT TO JL  ***************************************** */
    //jlmax = isp[il + *n] - 1;
	jlmax = (static_cast<integer>(isp[il + *n])) - 1;
    ju = jl + jlmax;
    //jumax = isp[iu + *n] - 1;
	jumax = (static_cast<integer>(isp[iu + *n])) - 1;
    if (jumax <= 0) {
	goto L2;
    }
    i__1 = jumax;
    for (j = 1; j <= i__1; ++j) {
/* L1: */
	isp[ju + j - 1] = isp[jutmp + j - 1];
    }

/*  ******  CALL REMAINING SUBROUTINES  ********************************* */
L2:
    //jlmax = isp[il + *n] - 1;
	jlmax = (static_cast<integer>(isp[il + *n])) - 1;
    ju = jl + jlmax;
    //jumax = isp[iu + *n] - 1;
	jumax = (static_cast<integer>(isp[iu + *n])) - 1;
    l = (ju + jumax - 2 + lratio) / lratio + 1;
    lmax = jlmax;
    d__ = l + lmax;
    u = d__ + *n;
    row = *nsp + 1 - *n;
    tmp = row - *n;
    umax = tmp - u;
    *esp = umax - jumax;

    if ((*path - 1) * (*path - 2) != 0) {
	goto L3;
    }
    if (umax <= 0) {
	goto L110;
    }
    nnf8_(n, &r__[1], &c__[1], &ic[1], &ia[1], &ja[1], &a[1], &z__[1], &b[1], 
	    &isp[il], &isp[jl], &rsp[l], &lmax, &rsp[d__], &isp[iu], &isp[ju],
	     &rsp[u], &umax, &rsp[row], &rsp[tmp], flag__);
    if (*flag__ != 0) {
	goto L100;
    }
    return 0;

L3:
    if (*path - 3 != 0) {
	goto L4;
    }
    nns8_(n, &r__[1], &c__[1], &isp[il], &isp[jl], &rsp[l], &rsp[d__], &isp[
	    iu], &isp[ju], &rsp[u], &z__[1], &b[1], &rsp[tmp]);

L4:
    if (*path - 4 != 0) {
	goto L5;
    }
    nnt8_(n, &r__[1], &c__[1], &isp[il], &isp[jl], &rsp[l], &rsp[d__], &isp[
	    iu], &isp[ju], &rsp[u], &z__[1], &b[1], &rsp[tmp]);
L5:
    return 0;

/* ** ERROR:  ERROR DETECTED IN NSF, NNF, NNS, OR NNT */
L100:
    return 0;
/* ** ERROR:  INSUFFICIENT STORAGE */
L110:
    *flag__ = *n * 10 + 1;
    return 0;
/* ** ERROR:  ILLEGAL PATH SPECIFICATION */
L111:
    *flag__ = *n * 11 + 1;
    return 0;
} /* ndrv_ */


/*       ---------------------------------------------------------------- */

/*               YALE SPARSE MATRIX PACKAGE - NONSYMMETRIC CODES */
/*                    SOLVING THE SYSTEM OF EQUATIONS MX = B */
/*                        (UNCOMPRESSED POINTER STORAGE) */

/*    I.   CALLING SEQUENCES */
/*         THE COEFFICIENT MATRIX CAN BE PROCESSED BY AN ORDERING ROUTINE */
/*    (E.G., TO REDUCE FILLIN OR ENSURE NUMERICAL STABILITY) BEFORE USING */
/*    THE REMAINING SUBROUTINES.  IF NO REORDERING IS DONE, THEN SET */
/*    R(I) = C(I) = IC(I) = I  FOR I=1,...,N.  THE CALLING SEQUENCE IS -- */
/*        (      (MATRIX ORDERING)) */
/*         NSF   (SYMBOLIC FACTORIZATION TO DETERMINE WHERE FILLIN WILL */
/*                 OCCUR DURING NUMERIC FACTORIZATION) */
/*         NNF   (NUMERIC FACTORIZATION INTO PRODUCT LDU OF UNIT LOWER */
/*                 TRIANGULAR MATRIX L, DIAGONAL MATRIX D, AND UNIT UPPER */
/*                 TRIANGULAR MATRIX U, AND SOLUTION OF LINEAR SYSTEM) */
/*         NNS   (SOLUTION OF LINEAR SYSTEM FOR ADDITIONAL RIGHT-HAND */
/*     OR          SIDE USING LDU FACTORIZATION FROM NNF) */
/*         NNT   (SOLUTION OF TRANSPOSED LINEAR SYSTEM FOR ADDITIONAL */
/*                 RIGHT-HAND SIDE USING LDU FACTORIZATION FROM NNF) */

/*    II.  STORAGE OF SPARSE MATRICES */
/*         THE NONZERO ENTRIES OF THE COEFFICIENT MATRIX M ARE STORED */
/*    ROW-BY-ROW IN THE ARRAY A.  TO IDENTIFY THE INDIVIDUAL NONZERO */
/*    ENTRIES IN EACH ROW, WE NEED TO KNOW IN WHICH COLUMN EACH ENTRY */
/*    LIES.  THE COLUMN INDICES WHICH CORRESPOND TO THE NONZERO ENTRIES */
/*    OF M ARE STORED IN THE ARRAY JA;  I.E., IF  A(K) = M(I,J),  THEN */
/*    JA(K) = J.  IN ADDITION, WE NEED TO KNOW WHERE EACH ROW STARTS AND */
/*    HOW LONG IT IS.  THE INDEX POSITIONS IN JA AND A WHERE THE ROWS OF */
/*    M BEGIN ARE STORED IN THE ARRAY IA;  I.E., IF M(I,J) IS THE FIRST */
/*    NONZERO ENTRY (STORED) IN THE I-TH ROW AND A(K) = M(I,J),  THEN */
/*    IA(I) = K.  MOREOVER, THE INDEX IN JA AND A OF THE FIRST LOCATION */
/*    FOLLOWING THE LAST ELEMENT IN THE LAST ROW IS STORED IN IA(N+1). */
/*    THUS, THE NUMBER OF ENTRIES IN THE I-TH ROW IS GIVEN BY */
/*    IA(I+1) - IA(I),  THE NONZERO ENTRIES OF THE I-TH ROW ARE STORED */
/*    CONSECUTIVELY IN */
/*            A(IA(I)),  A(IA(I)+1),  ..., A(IA(I+1)-1), */
/*    AND THE CORRESPONDING COLUMN INDICES ARE STORED CONSECUTIVELY IN */
/*            JA(IA(I)), JA(IA(I)+1), ..., JA(IA(I+1)-1). */
/*    FOR EXAMPLE, THE 5 BY 5 MATRIX */
/*                ( 1. 0. 2. 0. 0.) */
/*                ( 0. 3. 0. 0. 0.) */
/*            M = ( 0. 4. 5. 6. 0.) */
/*                ( 0. 0. 0. 7. 0.) */
/*                ( 0. 0. 0. 8. 9.) */
/*    WOULD BE STORED AS */
/*                 1  2  3  4  5  6  7  8  9 */
/*            ---+-------------------------- */
/*            IA   1  3  4  7  8 10 */
/*            JA   1  3  2  2  3  4  4  4  5 */
/*             A   1. 2. 3. 4. 5. 6. 7. 8. 9.         . */

/*         THE STRICT TRIANGULAR PORTIONS OF THE MATRICES L AND U ARE */
/*    STORED IN THE SAME FASHION USING THE ARRAYS  IL, JL, L  AND */
/*    IU, JU, U  RESPECTIVELY.  THE DIAGONAL ENTRIES OF L AND U ARE */
/*    ASSUMED TO BE EQUAL TO ONE AND ARE NOT STORED.  THE ARRAY D */
/*    CONTAINS THE RECIPROCALS OF THE DIAGONAL ENTRIES OF THE MATRIX D. */

/*    III. ADDITIONAL STORAGE SAVINGS */
/*         IN NSF, R AND IC CAN BE THE SAME ARRAY IN THE CALLING */
/*    SEQUENCE IF NO REORDERING OF THE COEFFICIENT MATRIX HAS BEEN DONE. */
/*         IN NNF, R, C AND IC CAN ALL BE THE SAME ARRAY IF NO REORDERING */
/*    HAS BEEN DONE.  IF ONLY THE ROWS HAVE BEEN REORDERED, THEN C AND IC */
/*    CAN BE THE SAME ARRAY.  IF THE ROW AND COLUMN ORDERINGS ARE THE */
/*    SAME, THEN R AND C CAN BE THE SAME ARRAY.  Z AND ROW CAN BE THE */
/*    SAME ARRAY. */
/*         IN NNS OR NNT, R AND C CAN BE THE SAME ARRAY IF NO REORDERING */
/*    HAS BEEN DONE OR IF THE ROW AND COLUMN ORDERINGS ARE THE SAME.  Z */
/*    AND B CAN BE THE SAME ARRAY;  HOWEVER, THEN B WILL BE DESTROYED. */

/*    IV.  PARAMETERS */
/*         FOLLOWING IS A LIST OF PARAMETERS TO THE PROGRAMS.  NAMES ARE */
/*    UNIFORM AMONG THE VARIOUS SUBROUTINES.  CLASS ABBREVIATIONS ARE -- */
/*       N - INTEGER VARIABLE */
/*       F - REAL VARIABLE */
/*       V - SUPPLIES A VALUE TO A SUBROUTINE */
/*       R - RETURNS A RESULT FROM A SUBROUTINE */
/*       I - USED INTERNALLY BY A SUBROUTINE */
/*       A - ARRAY */

/* CLASS   PARAMETER */
/* ------+---------- */
/* FVA     A     - NONZERO ENTRIES OF THE COEFFICIENT MATRIX M, STORED */
/*                   BY ROWS. */
/*                   SIZE = NUMBER OF NONZERO ENTRIES IN M. */
/* FVA     B     - RIGHT-HAND SIDE B. */
/*                   SIZE = N. */
/* NVA     C     - ORDERING OF THE COLUMNS OF M. */
/*                   SIZE = N. */
/* FVRA    D     - RECIPROCALS OF THE DIAGONAL ENTRIES OF THE MATRIX D. */
/*                   SIZE = N. */
/* NR      FLAG  - ERROR FLAG;  VALUES AND THEIR MEANINGS ARE -- */
/*                    0     NO ERRORS DETECTED */
/*                    N+K   NULL ROW IN A  --  ROW = K */
/*                   2N+K   DUPLICATE ENTRY IN A  --  ROW = K */
/*                   3N+K   INSUFFICIENT STORAGE FOR JL  --  ROW = K */
/*                   4N+1   INSUFFICIENT STORAGE FOR L */
/*                   5N+K   NULL PIVOT  --  ROW = K */
/*                   6N+K   INSUFFICIENT STORAGE FOR JU  --  ROW = K */
/*                   7N+1   INSUFFICIENT STORAGE FOR U */
/*                   8N+K   ZERO PIVOT  --  ROW = K */
/* NVA     IA    - POINTERS TO DELIMIT THE ROWS IN A. */
/*                   SIZE = N+1. */
/* NVA     IC    - INVERSE OF THE ORDERING OF THE COLUMNS OF M;  I.E., */
/*                   IC(C(I) = I  FOR I=1,...N. */
/*                   SIZE = N. */
/* NVRA    IL    - POINTERS TO DELIMIT THE ROWS IN L. */
/*                   SIZE = N+1. */
/* NVRA    IU    - POINTERS TO DELIMIT THE ROWS IN U. */
/*                   SIZE = N+1. */
/* NVA     JA    - COLUMN NUMBERS CORRESPONDING TO THE ELEMENTS OF A. */
/*                   SIZE = SIZE OF A. */
/* NVRA    JL    - COLUMN NUMBERS CORRESPONDING TO THE ELEMENTS OF L. */
/*                   SIZE = JLMAX. */
/* NV      JLMAX - DECLARED DIMENSION OF JL;  JLMAX MUST BE LARGER THAN */
/*                   THE NUMBER OF NONZERO ENTRIES IN THE STRICT LOWER */
/*                   TRIANGLE OF M PLUS FILLIN  (IL(N+1)-1 AFTER NSF). */
/* NVRA    JU    - COLUMN NUMBERS CORRESPONDING TO THE ELEMENTS OF U. */
/*                   SIZE = JUMAX. */
/* NV      JUMAX - DECLARED DIMENSION OF JU;  JUMAX MUST BE LARGER THAN */
/*                   THE NUMBER OF NONZERO ENTRIES IN THE STRICT UPPER */
/*                   TRIANGLE OF M PLUS FILLIN  (IU(N+1)-1 AFTER NSF). */
/* FVRA    L     - NONZERO ENTRIES IN THE STRICT LOWER TRIANGULAR PORTION */
/*                   OF THE MATRIX L, STORED BY ROWS. */
/*                   SIZE = LMAX */
/* NV      LMAX  - DECLARED DIMENSION OF L;  LMAX MUST BE LARGER THAN */
/*                   THE NUMBER OF NONZERO ENTRIES IN THE STRICT LOWER */
/*                   TRIANGLE OF M PLUS FILLIN  (IL(N+1)-1 AFTER NSF). */
/* NV      N     - NUMBER OF VARIABLES/EQUATIONS. */
/* NVA     R     - ORDERING OF THE ROWS OF M. */
/*                   SIZE = N. */
/* FVRA    U     - NONZERO ENTRIES IN THE STRICT UPPER TRIANGULAR PORTION */
/*                   OF THE MATRIX U, STORED BY ROWS. */
/*                   SIZE = UMAX. */
/* NV      UMAX  - DECLARED DIMENSION OF U;  UMAX MUST BE LARGER THAN */
/*                   THE NUMBER OF NONZERO ENTRIES IN THE STRICT UPPER */
/*                   TRIANGLE OF M PLUS FILLIN  (IU(N+1)-1 AFTER NSF). */
/* FRA     Z     - SOLUTION X. */
/*                   SIZE = N. */


/*       ---------------------------------------------------------------- */

/* *** SUBROUTINE NSF */
/* *** SYMBOLIC LDU-FACTORIZATION OF A NONSYMMETRIC SPARSE MATRIX */
/*      (UNCOMPRESSED POINTER STORAGE) */

/* Subroutine */ integer nsf8_(integer *n, integer *r__, integer *ic, integer *ia,
	 integer *ja, /*integer*/ doublereal *il, /*integer*/ doublereal *jl, integer *jlmax, /*integer*/ doublereal *iu, 
	/*integer*/ doublereal *ju, integer *jumax, /*integer*/ doublereal *q, /*integer*/ doublereal *im, integer *flag__)
{
    /* System generated locals */
    integer i__1=0, i__2=0;

    /* Local variables */
    integer i__=0, j=0, k=0, m=0, qm=0, vj=0, jmin=0, jmax=0, jlptr=0, juptr=0;


/*       INPUT VARIABLES:   N, R,IC, IA,JA, JLMAX, JUMAX. */
/*       OUTPUT VARIABLES:  IL,JL, IU,JU, FLAG. */

/*       PARAMETERS USED INTERNALLY: */
/* NIA     Q     - SUPPOSE M' IS THE RESULT OF REORDERING M;  IF */
/*                   PROCESSING OF THE KTH ROW OF M' (HENCE THE KTH ROWS */
/*                   OF L AND U) IS BEING DONE, THEN Q(J) IS INITIALLY */
/*                   NONZERO IF M'(K,J) IS NONZERO;  SINCE VALUES NEED */
/*                   NOT BE STORED, EACH ENTRY POINTS TO THE NEXT */
/*                   NONZERO;  FOR EXAMPLE, IF  N=9  AND THE 5TH ROW OF */
/*                   M' IS */
/*                           0 X X 0 X 0 0 X 0, */
/*                   THEN Q WILL INITIALLY BE */
/*                           A 3 5 A 8 A A 10 A 2        (A - ARBITRARY); */
/*                   Q(N+1) POINTS TO THE FIRST NONZERO IN THE ROW AND */
/*                   THE LAST NONZERO POINTS TO  N+1;  AS THE ALGORITHM */
/*                   PROCEEDS, OTHER ELEMENTS OF Q ARE INSERTED IN THE */
/*                   LIST BECAUSE OF FILLIN. */
/*                   SIZE = N+1. */
/* NIA     IM    - AT EACH STEP IN THE FACTORIZATION, IM(I) IS THE LAST */
/*                   ELEMENT IN THE ITH ROW OF U WHICH NEEDS TO BE */
/*                   CONSIDERED IN COMPUTING FILLIN. */
/*                   SIZE = N. */

/*  INTERNAL VARIABLES-- */
/*    JLPTR - POINTS TO THE LAST POSITION USED IN  JL. */
/*    JUPTR - POINTS TO THE LAST POSITION USED IN  JU. */


/*  ******  INITIALIZE POINTERS  **************************************** */
    /* Parameter adjustments */
    --im;
    --q;
    --ju;
    --iu;
    --jl;
    --il;
    --ja;
    --ia;
    --ic;
    --r__;

    /* Function Body */
    jlptr = 0;
   // il[1] = 1;
	 il[1] = 1.0;
    juptr = 0;
   // iu[1] = 1;
	iu[1] = 1.0;

/*  ******  FOR EACH ROW OF L AND U  ************************************ */
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
/*  ******  SET Q TO THE REORDERED ROW OF A  **************************** */
	//q[*n + 1] = *n + 1;
	q[*n + 1] = static_cast<doublereal>(*n + 1);
	jmin = ia[r__[k]];
	jmax = ia[r__[k] + 1] - 1;
	if (jmin > jmax) {
	    goto L101;
	}
	i__2 = jmax;
	for (j = jmin; j <= i__2; ++j) {
	    vj = ic[ja[j]];
	    qm = *n + 1;
L1:
	    m = qm;
	   // qm = q[m];
		 qm = static_cast<integer>(q[m]);
	    if (qm < vj) {
		goto L1;
	    }
	    if (qm == vj) {
		goto L102;
	    }
	   // q[m] = vj;
	   // q[vj] = qm;
		 q[m] = static_cast<doublereal>(vj);
	    q[vj] = static_cast<doublereal>(qm);
/* L2: */
	}

/*  ******  FOR EACH ENTRY IN THE LOWER TRIANGLE  *********************** */
	i__ = *n + 1;
L3:
	//i__ = q[i__];
	i__ = static_cast<integer>(q[i__]);
	if (i__ >= k) {
	    goto L7;
	}
/*  ******  L(K,I) WILL BE NONZERO, SO ADD IT TO JL  ******************** */
	++jlptr;
	if (jlptr > *jlmax) {
	    goto L103;
	}
	//jl[jlptr] = i__;
	jl[jlptr] = static_cast<doublereal>(i__);
	qm = i__;
/*  ******  INSPECT ITH ROW FOR FILLIN, ADJUST IM IF POSSIBLE  ********** */
	//jmin = iu[i__];
	jmin = static_cast<integer>(iu[i__]);
	//jmax = im[i__];
	jmax = static_cast<integer>(im[i__]);
	if (jmin > jmax) {
	    goto L6;
	}
	i__2 = jmax;
	for (j = jmin; j <= i__2; ++j) {
	   // vj = ju[j];
		vj = static_cast<integer>(ju[j]);
	    if (vj == k) {
		//im[i__] = j;
			im[i__] = static_cast<doublereal>(j);
	    }
L4:
	    m = qm;
	    //qm = q[m];
		qm = static_cast<integer>(q[m]);
	    if (qm < vj) {
		goto L4;
	    }
	    if (qm == vj) {
		goto L5;
	    }
	    //q[m] = vj;
	    //q[vj] = qm;
		q[m] = static_cast<doublereal>(vj);
	    q[vj] = static_cast<doublereal>(qm);
	    qm = vj;
L5:
	    ;
	}
L6:
	goto L3;

/*  ******  CHECK FOR NULL PIVOT  *************************************** */
L7:
	if (i__ != k) {
	    goto L105;
	}
/*  ******  REMAINING ELEMENTS OF Q DEFINE STRUCTURE OF U(K, )  ********* */
L8:
	//i__ = q[i__];
	i__ = static_cast<integer>(q[i__]);
	if (i__ > *n) {
	    goto L9;
	}
	++juptr;
	if (juptr > *jumax) {
	    goto L106;
	}
	//ju[juptr] = i__;
	ju[juptr] = static_cast<doublereal>(i__);
	goto L8;
/*  ******  GET READY FOR NEXT ROW  ************************************* */
L9:
	//im[k] = juptr;
	im[k] = static_cast<doublereal>(juptr);
	//il[k + 1] = jlptr + 1;
	il[k + 1] = static_cast<doublereal>(jlptr + 1);
/* L10: */
	//iu[k + 1] = juptr + 1;
	iu[k + 1] = static_cast<doublereal>(juptr + 1);
    }

    *flag__ = 0;
    return 0;

/* ** ERROR:  NULL ROW IN A */
L101:
    *flag__ = *n + r__[k];
    return 0;
/* ** ERROR:  DUPLICATE ENTRY IN A */
L102:
    *flag__ = (*n << 1) + r__[k];
    return 0;
/* ** ERROR:  INSUFFICIENT STORAGE FOR JL */
L103:
    *flag__ = *n * 3 + k;
    return 0;
/* ** ERROR:  NULL PIVOT */
L105:
    *flag__ = *n * 5 + k;
    return 0;
/* ** ERROR:  INSUFFICIENT STORAGE FOR JU */
L106:
    *flag__ = *n * 6 + k;
    return 0;
} /* nsf8_ */


/*       ---------------------------------------------------------------- */

/* *** SUBROUTINE NNF */
/* *** NUMERIC LDU-FACTORIZATION OF SPARSE NONSYMMETRIC MATRIX AND */
/*      SOLUTION OF SYSTEM OF LINEAR EQUATIONS (UNCOMPRESSED POINTER */
/*      STORAGE) */

/* Subroutine */ integer nnf8_(integer *n, integer *r__, integer *c__, integer *
	ic, integer *ia, integer *ja, doublereal *a, doublereal *z__, 
	doublereal *b, /*integer*/ doublereal *il, /*integer*/ doublereal *jl, doublereal *l, integer *lmax,
	 doublereal *d__, /*integer*/ doublereal *iu, /*integer*/ doublereal *ju, doublereal *u, integer *
	umax, doublereal *row, doublereal *tmp, integer *flag__)
{
    /* System generated locals */
    integer i__1=0, i__2=0, i__3=0;

    /* Local variables */
    integer i__=0, j=0, k=0;
    doublereal dk=0.0, li=0.0, sum=0.0;
    integer imin=0, jmin=0, imax=0, jmax=0;


/*       INPUT VARIABLES:   N, R,C,IC, IA,JA,A, B, IL,JL,LMAX, IU,JU,UMAX */
/*       OUTPUT VARIABLES:  Z, L,D,U, FLAG */

/*       PARAMETERS USED INTERNALLY: */
/* FIA     ROW   - HOLDS INTERMEDIATE VALUES IN CALCULATION OF L, D, U. */
/*                   SIZE = N. */
/* FIA     TMP   - HOLDS NEW RIGHT-HAND SIDE B' FOR SOLUTION OF THE */
/*                   EQUATION  UX = B'. */
/*                   SIZE = N. */

/*       REAL  A(1), Z(1), B(1),  L(1), D(1), U(1), */
/*    *     ROW(1), TMP(1),  LI, SUM, DK */

/*  ******  CHECK STORAGE  ********************************************** */
    /* Parameter adjustments */
    --tmp;
    --row;
    --u;
    --ju;
    --iu;
    --d__;
    --l;
    --jl;
    --il;
    --b;
    --z__;
    --a;
    --ja;
    --ia;
    --ic;
    --c__;
    --r__;

    /* Function Body */
    if ((static_cast<integer>(il[*n + 1])) - 1 > *lmax) {
	goto L104;
    }
    if ((static_cast<integer>(iu[*n + 1])) - 1 > *umax) {
	goto L107;
    }

/*  ******  FOR EACH ROW  *********************************************** */
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
/*  ******  SET THE INITIAL STRUCTURE OF ROW  *************************** */
	//jmin = il[k];
	//jmax = il[k + 1] - 1;
	jmin = static_cast<integer>(il[k]);
	jmax = static_cast<integer>(il[k + 1]) - 1;
	if (jmin > jmax) {
	    goto L2;
	}
/*  ******  IF L(K,M) .NE. 0, ROW(M)=0  ********************************* */
	i__2 = jmax;
	for (j = jmin; j <= i__2; ++j) {
/* L1: */
	   // row[jl[j]] = 0.;
		 row[(static_cast<integer>(jl[j]))] = 0.;
	}
L2:
	row[k] = 0.;
	//jmin = iu[k];
	//jmax = iu[k + 1] - 1;
	jmin = static_cast<integer>(iu[k]);
	jmax = static_cast<integer>(iu[k + 1]) - 1;
	if (jmin > jmax) {
	    goto L4;
	}
/*  ******  IF U(K,M) .NE. 0, ROW(M)=0  ********************************* */
	i__2 = jmax;
	for (j = jmin; j <= i__2; ++j) {
/* L3: */
	    //row[ju[j]] = 0.;
		 row[static_cast<integer>(ju[j])] = 0.;
	}
L4:
	jmin = ia[r__[k]];
	jmax = ia[r__[k] + 1] - 1;
/*  ******  SET ROW TO KTH ROW OF REORDERED A  ************************** */
	i__2 = jmax;
	for (j = jmin; j <= i__2; ++j) {
/* L5: */
	    row[ic[ja[j]]] = a[j];
	}
/*  ******  INITIALIZE SUM  ********************************************* */
	sum = b[r__[k]];

/*  ******  ASSIGN THE KTH ROW OF L AND ADJUST ROW, SUM  **************** */
	//imin = il[k];
	//imax = il[k + 1] - 1;
	imin = static_cast<integer>(il[k]);
	imax = static_cast<integer>(il[k + 1]) - 1;
	if (imin > imax) {
	    goto L8;
	}
	i__2 = imax;
	for (i__ = imin; i__ <= i__2; ++i__) {
	   // li = -row[jl[i__]];
/*  ******  IF L IS NOT REQUIRED, THEN COMMENT OUT THE FOLLOWING LINE  ** */
	    //l[i__] = -li;
	    //sum += li * tmp[jl[i__]];
	    //jmin = iu[jl[i__]];
	    //jmax = iu[jl[i__] + 1] - 1;

		integer i__ind=static_cast<integer>(jl[i__]);
		li = -row[i__ind];
/*  ******  IF L IS NOT REQUIRED, THEN COMMENT OUT THE FOLLOWING LINE  ** */
	    l[i__] = -li;
	    sum += li * tmp[i__ind];
	    jmin = static_cast<integer>(iu[i__ind]);
	    jmax = static_cast<integer>(iu[i__ind + 1]) - 1;

	    if (jmin > jmax) {
		goto L7;
	    }
	    i__3 = jmax;
	    for (j = jmin; j <= i__3; ++j) {
/* L6: */
		//row[ju[j]] += li * u[j];
			row[static_cast<integer>(ju[j])] += li * u[j];
	    }
L7:
	    ;
	}

/*  ******  ASSIGN DIAGONAL D AND KTH ROW OF U, SET TMP(K)  ************* */
L8:
	if (row[k] == 0.) {
	    goto L108;
	}
	dk = 1 / row[k];
	d__[k] = dk;
	tmp[k] = sum * dk;
	//jmin = iu[k];
	//jmax = iu[k + 1] - 1;
	jmin = static_cast<integer>(iu[k]);
	jmax = static_cast<integer>(iu[k + 1]) - 1;

	if (jmin > jmax) {
	    goto L10;
	}
	i__2 = jmax;
	for (j = jmin; j <= i__2; ++j) {
/* L9: */
	   // u[j] = row[ju[j]] * dk;
		u[j] = row[static_cast<integer>(ju[j])] * dk;
	}
L10:
	;
    }

/*  ******  SOLVE  UX = TMP  BY BACK SUBSTITUTION  ********************** */
    k = *n;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	sum = tmp[k];
	//jmin = iu[k];
	//jmax = iu[k + 1] - 1;
	jmin = static_cast<integer>(iu[k]);
	jmax = static_cast<integer>(iu[k + 1]) - 1;

	if (jmin > jmax) {
	    goto L12;
	}
	i__2 = jmax;
	for (j = jmin; j <= i__2; ++j) {
/* L11: */
	   // sum -= u[j] * tmp[ju[j]];
		sum -= u[j] * tmp[static_cast<integer>(ju[j])];
	}
L12:
	tmp[k] = sum;
	z__[c__[k]] = sum;
/* L13: */
	--k;
    }

    *flag__ = 0;
    return 0;

/* ** ERROR:  INSUFFICIENT STORAGE FOR L */
L104:
    *flag__ = (*n << 2) + 1;
    return 0;
/* ** ERROR:  INSUFFICIENT STORAGE FOR U */
L107:
    *flag__ = *n * 7 + 1;
    return 0;
/* ** ERROR:  ZERO PIVOT */
L108:
    *flag__ = (*n << 3) + k;
    return 0;
} /* nnf8_ */


/*       ---------------------------------------------------------------- */

/* *** SUBROUTINE NNS 
/* *** NUMERIC SOLUTION OF A SPARSE NONSYMMETRIC SYSTEM OF LINEAR */
/*      EQUATIONS GIVEN LDU-FACTORIZATION (UNCOMPRESSED POINTER STORAGE) */

/* Subroutine */ integer nns8_(integer *n, integer *r__, integer *c__, /*integer*/ doublereal *
	il, /*integer*/ doublereal *jl, doublereal *l, doublereal *d__, /*integer*/ doublereal *iu, /*integer*/ 
	doublereal *ju, doublereal *u, doublereal *z__, doublereal *b, doublereal *tmp)
{
    /* System generated locals */
    integer i__1=0, i__2=0;

    /* Local variables */
    integer i__=0, j=0, k=0;
    doublereal sum=0.0;
    integer jmin=0, jmax=0;


/*       INPUT VARIABLES:   N, R,C, IL,JL,L, D, IU,JU,U, B */
/*       OUTPUT VARIABLES:  Z */

/*       PARAMETERS USED INTERNALLY: */
/* FIA     TMP   - HOLDS NEW RIGHT-HAND SIDE B' FOR SOLUTION OF THE */
/*                   EQUATION UX = B'. */
/*                   SIZE = N. */

/*       REAL  L(1), D(1), U(1),  Z(1), B(1),  TMP(1), SUM */

/*  ******  SOLVE LDY = B  BY FORWARD SUBSTITUTION  ********************* */
    /* Parameter adjustments */
    --tmp;
    --b;
    --z__;
    --u;
    --ju;
    --iu;
    --d__;
    --l;
    --jl;
    --il;
    --c__;
    --r__;

    /* Function Body */
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	sum = b[r__[k]];
	//jmin = il[k];
	//jmax = il[k + 1] - 1;
	jmin = static_cast<integer>(il[k]);
	jmax = static_cast<integer>(il[k + 1]) - 1;
	if (jmin > jmax) {
	    goto L2;
	}
	i__2 = jmax;
	for (j = jmin; j <= i__2; ++j) {
/* L1: */
	    //sum -= l[j] * tmp[jl[j]];
		sum -= l[j] * tmp[static_cast<integer>(jl[j])];
	}
L2:
	tmp[k] = sum * d__[k];
    }

/*  ******  SOLVE  UX = Y  BY BACK SUBSTITUTION  ************************ */
    k = *n;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	sum = tmp[k];
	//jmin = iu[k];
	//jmax = iu[k + 1] - 1;
	jmin = static_cast<integer>(iu[k]);
	jmax = static_cast<integer>(iu[k + 1]) - 1;
	if (jmin > jmax) {
	    goto L4;
	}
	i__2 = jmax;
	for (j = jmin; j <= i__2; ++j) {
/* L3: */
	    //sum -= u[j] * tmp[ju[j]];
		sum -= u[j] * tmp[static_cast<integer>(ju[j])];
	}
L4:
	tmp[k] = sum;
	z__[c__[k]] = sum;
/* L5: */
	--k;
    }
    return 0;
} /* nns8_ */

/*       ---------------------------------------------------------------- */

/* *** SUBROUTINE NNT */
/* *** NUMERIC SOLUTION OF THE TRANSPOSE OF A SPARSE NONSYMMETRIC SYSTEM */
/*      OF LINEAR EQUATIONS GIVEN LDU-FACTORIZATION (UNCOMPRESSED POINTER */
/*      STORAGE) */

/* Subroutine */ integer nnt8_(integer *n, integer *r__, integer *c__, /*integer*/ doublereal *
	il, /*integer*/ doublereal *jl, doublereal *l, doublereal *d__, /*integer*/ doublereal *iu, /*integer*/ doublereal 
	*ju, doublereal *u, doublereal *z__, doublereal *b, doublereal *tmp)
{
    /* System generated locals */
    integer i__1=0, i__2=0;

    /* Local variables */
    integer i__=0, j=0, k=0, jmin=0, jmax=0;
    doublereal tmpk=0.0;


/*       INPUT VARIABLES:   N, R,C, IL,JL,L, D, IU,JU,U, B */
/*       OUTPUT VARIABLES:  Z */

/*       PARAMETERS USED INTERNALLY: */
/* FIA     TMP   - HOLDS NEW RIGHT-HAND SIDE B' FOR SOLUTION OF THE */
/*                   EQUATION LX = B'. */
/*                   SIZE = N. */

/*       REAL  L(1), D(1), U(1),  Z(1), B(1),  TMP(1), TMPK */

/*  ******  SOLVE  UT Y = B  BY FORWARD SUBSTITUTION  ******************* */
    /* Parameter adjustments */
    --tmp;
    --b;
    --z__;
    --u;
    --ju;
    --iu;
    --d__;
    --l;
    --jl;
    --il;
    --c__;
    --r__;

    /* Function Body */
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
/* L1: */
	tmp[k] = b[c__[k]];
    }
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	tmpk = -tmp[k];
	//jmin = iu[k];
	//jmax = iu[k + 1] - 1;
	jmin = static_cast<integer>(iu[k]);
	jmax = static_cast<integer>(iu[k + 1]) - 1;
	if (jmin > jmax) {
	    goto L3;
	}
	i__2 = jmax;
	for (j = jmin; j <= i__2; ++j) {
/* L2: */
	   // tmp[ju[j]] += u[j] * tmpk;
		tmp[static_cast<integer>(ju[j])] += u[j] * tmpk;
	}
L3:
	;
    }

/*  ******  SOLVE  D LT X = Y  BY BACK SUBSTITUTION  ******************** */
    k = *n;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	tmpk = -tmp[k] * d__[k];
	//jmin = il[k];
	//jmax = il[k + 1] - 1;
	jmin = static_cast<integer>(il[k]);
	jmax = static_cast<integer>(il[k + 1]) - 1;
	if (jmin > jmax) {
	    goto L5;
	}
	i__2 = jmax;
	for (j = jmin; j <= i__2; ++j) {
/* L4: */
	    //tmp[jl[j]] += l[j] * tmpk;
		tmp[static_cast<integer>(jl[j])] += l[j] * tmpk;
	}
L5:
	z__[r__[k]] = -tmpk;
/* L6: */
	--k;
    }
    return 0;
} /* nnt8_ */


// Для генерации матрицы СЛАУ требуется в случае реализации
// на динамических массивах переупорядочивание элементов:
// сортировка. Здесь будет реализована быстрая сортировка.
// Брайан Керниган и Денис Ритчи "The C programming language".
// swap: Обмен местами v[i] и v[j]
void swapCSIRamg(integer * &v, doublereal * &dr, integer i, integer j)
{
        integer tempi=0;
		doublereal tempr=0.0;

		// change v[i] <-> v[j]
		tempi = v[i];
		v[i] = v[j];
		v[j] = tempi;
		// change dr[i] <-> dr[j]
		tempr = dr[i];
		dr[i] = dr[j];
		dr[j] = tempr;

} // swapCSIRamg

// Вот алгоритм PivotList
integer PivotListCSIRamg(integer * &jptr, doublereal * &altr, integer first, integer last) {
	// list==jptr and altr обрабатываемый список
	// first номер первого элемента
	// last номер последнего элемента

	integer PivotValue = jptr[first];
	integer PivotPoint = first;

	for (integer index=(first+1); index<=last; index++) {
		if (jptr[index]<PivotValue) {
			PivotPoint++;
			swapCSIRamg(jptr, altr, PivotPoint, index);
		}
	}

	swapCSIRamg(jptr, altr, first, PivotPoint);

	return PivotPoint;
} // PivotListamg


// Быстрая сортировка Хоара.
// Запрограммировано с использованием ДЖ. Макконелл Анализ алгоритмов
// стр. 106.
void QuickSortCSIR_amg(integer * &jptr, doublereal * &altr, integer first, integer last) {
	// list упорядочиваемый список элементов
	// first номер первого элемента в сортируемой части списка
	// last номер последнего элемента в сортируемой части списка

	if (0) {
		// BubbleSort
		integer numberOfPairs=last-first+1;
		bool swappedElements=true;
		while (swappedElements) {
			 numberOfPairs--;
			 swappedElements=false;
			 for (integer i=first; i<=first+numberOfPairs-1; ++i) {
				 if (jptr[i]>jptr[i+1]) {
					 swapCSIRamg(jptr, altr, i, i+1);
					 swappedElements=true;
				 }
			 }
		}
	}
	else
	{
	    integer pivot;

	    if (first < last) {
             pivot = PivotListCSIRamg(jptr, altr, first, last);
             QuickSortCSIR_amg(jptr, altr, first, pivot-1);
		     QuickSortCSIR_amg(jptr, altr, pivot+1, last);
	    }
	}
} // QuickSortCSIR_amg

// переменная с глобальной областью видимости необходима для 
// защиты от холостого рестарта, (перезапуск на сошедшмся решении).
// Это должно существенным образом экономить время пользователя.
// Нехолостой рестарт необходим на нелинейных задачах.
// 23 июля 2015.
doublereal finish_residual=0.0;

/* 
13 апреля 2013 года вся память выделяемая и уничтожаемая внутри
BICGSTAB_internal3 вынесена наружу, с целью препятствовать частым выделениям и уничтожениям памяти.
Это должно положительным образом сказаться на скорости работы BiCGStab.
*/
typedef struct TQuickMemVorst {

	 //doublereal *rthdsd; // правая часть системы уравнений

	 // Исходная матрица в формате CRS.
	 doublereal *val;
     integer *col_ind;
	 integer *row_ptr;

	 // Рабочие вектора.
	 doublereal *ri, *roc, *s, *t, *vi, *pi, *dx, *dax;
	 doublereal *y, *z; // результат предобуславливания

	 doublereal *a; // CRS
	 integer *ja;
	 integer *ia;
	 // ILU предобуславливатель:
	 // Будем сразу хранить в MSR формате, который заведём с запасом по памяти.
	 doublereal *alurc;
	 integer *jlurc;
	 integer *jurc;
	 doublereal *vec;
	 doublereal *alu;
	 integer *jlu;
	 integer *ju;
	 // Копия для распараллеливания lusol_2
	 // Копия матрицы разложения нужна для распараллеливания обратного хода по U матрице.
	 doublereal *x1; // Копия искомого решения.
	 doublereal *alu1; // Копия матрицы 
	 integer *jlu1; // ILU2 разложения.
	 integer *ju1; // это тоже относится к копии матрицы.
	 integer *iw; // для ILU(0)
	 // для ILU(lfil)
	 integer* levs; 
	 doublereal* w; 
	 integer* jw;
	 doublereal* w_dubl; // для распараллеливания iluk.c
	 integer* jw_dubl;

	 integer iwk; // размерность памяти под матрицу предобуславливания.

	 //doublereal *trthdsd; // правая часть системы уравнений

	 // Исходная матрица в формате CRS.
	 doublereal *tval;
     integer *tcol_ind;
	 integer *trow_ptr;

	  // Рабочие вектора.
	 doublereal *tri, *troc, *ts, *tt, *tvi, *tpi, *tdx, *tdax;
	 doublereal *ty, *tz;

	 bool ballocCRScfd;
	 bool bsignalfreeCRScfd; // сигнал к освобождению памяти из под a,ja,ia. если true
	 doublereal *ta; // CRS
	 integer *tja;
	 integer *tia;
	  // ILU предобуславливатель:
	 // Будем сразу хранить в MSR формате, который заведём с запасом по памяти.
	 doublereal *talu;
	 integer *tjlu;
	 integer *tju;
	 integer *tiw; // для ILU(0)
	 // для ILU(lfil)
	 integer* tlevs; // для ILU(lfil)
	 doublereal* tw; 
	 integer* tjw;

	 integer tiwk; // размерность памяти под матрицу предобуславливания.

	 bool ballocCRSt;
	 bool bsignalfreeCRSt; // сигнал к освобождению памяти из под a,ja,ia. если true

	 // Идея запомнить сколько итераций было сделано для компонент скорости 
	 // чтобы использовать такое-же количество итераций для температуры в задачах с естественой конвекцией.
	 integer icount_vel; 
} QuickMemVorst;



// Здесь содержится обвязка вызывающая amg1r5.
// локальное выделение памяти:всё внутри, многократные alloc и free.
void amg_loc_memory(equation3D* &sl, equation3D_bon* &slb,
			   integer maxelm, integer maxbound,
			   doublereal *dV, doublereal* &dX0, 
			   doublereal alpharelax, integer iVar, bool bLRfree, QuickMemVorst& m,
	integer iVorst_version)
{

	// Замер времени.
	unsigned int calculation_main_start_time; // начало счёта мс.
	unsigned int calculation_main_end_time; // окончание счёта мс.

	calculation_main_start_time=clock(); // момент начала счёта.

	

	 // На случай если память не была выделена.
	 if (dX0==NULL) {
	    dX0=new doublereal[maxelm+maxbound];	    
	    for (integer i=0; i<maxelm+maxbound; ++i) {
	        dX0[i]=0.0;
	    }
	 }

	
	const doublereal nonzeroEPS=1e-37; // для отделения вещественного нуля
	doublereal res_sum=0.0;
	res_sum=0.0;
	for (integer i=0; i<maxelm; ++i) {
		// внутренность матрицы.
		doublereal buf=0.0;
		buf=(sl[i].ap*dX0[sl[i].iP]-dV[sl[i].iP]);
		if ((sl[i].iB>-1) && (fabs(sl[i].ab) > nonzeroEPS)) buf-=sl[i].ab*dX0[sl[i].iB];
		if ((sl[i].iE>-1) && (fabs(sl[i].ae) > nonzeroEPS)) buf-=sl[i].ae*dX0[sl[i].iE];
		if ((sl[i].iN>-1) && (fabs(sl[i].an) > nonzeroEPS)) buf-=sl[i].an*dX0[sl[i].iN];
		if ((sl[i].iS>-1) && (fabs(sl[i].as) > nonzeroEPS)) buf-=sl[i].as*dX0[sl[i].iS];
		if ((sl[i].iT>-1) && (fabs(sl[i].at) > nonzeroEPS)) buf-=sl[i].at*dX0[sl[i].iT];
		if ((sl[i].iW>-1) && (fabs(sl[i].aw) > nonzeroEPS)) buf-=sl[i].aw*dX0[sl[i].iW];
		// Модификация на АЛИС сетке:
		if ((sl[i].iB2>-1) && (fabs(sl[i].ab2) > nonzeroEPS)) buf -= sl[i].ab2*dX0[sl[i].iB2];
		if ((sl[i].iE2>-1) && (fabs(sl[i].ae2) > nonzeroEPS)) buf -= sl[i].ae2*dX0[sl[i].iE2];
		if ((sl[i].iN2>-1) && (fabs(sl[i].an2) > nonzeroEPS)) buf -= sl[i].an2*dX0[sl[i].iN2];
		if ((sl[i].iS2>-1) && (fabs(sl[i].as2) > nonzeroEPS)) buf -= sl[i].as2*dX0[sl[i].iS2];
		if ((sl[i].iT2>-1) && (fabs(sl[i].at2) > nonzeroEPS)) buf -= sl[i].at2*dX0[sl[i].iT2];
		if ((sl[i].iW2>-1) && (fabs(sl[i].aw2) > nonzeroEPS)) buf -= sl[i].aw2*dX0[sl[i].iW2];
		// Модификация на АЛИС сетке:
		if ((sl[i].iB3>-1) && (fabs(sl[i].ab3) > nonzeroEPS)) buf -= sl[i].ab3*dX0[sl[i].iB3];
		if ((sl[i].iE3>-1) && (fabs(sl[i].ae3) > nonzeroEPS)) buf -= sl[i].ae3*dX0[sl[i].iE3];
		if ((sl[i].iN3>-1) && (fabs(sl[i].an3) > nonzeroEPS)) buf -= sl[i].an3*dX0[sl[i].iN3];
		if ((sl[i].iS3>-1) && (fabs(sl[i].as3) > nonzeroEPS)) buf -= sl[i].as3*dX0[sl[i].iS3];
		if ((sl[i].iT3>-1) && (fabs(sl[i].at3) > nonzeroEPS)) buf -= sl[i].at3*dX0[sl[i].iT3];
		if ((sl[i].iW3>-1) && (fabs(sl[i].aw3) > nonzeroEPS)) buf -= sl[i].aw3*dX0[sl[i].iW3];
		// Модификация на АЛИС сетке:
		if ((sl[i].iB4>-1) && (fabs(sl[i].ab4) > nonzeroEPS)) buf -= sl[i].ab4*dX0[sl[i].iB4];
		if ((sl[i].iE4>-1) && (fabs(sl[i].ae4) > nonzeroEPS)) buf -= sl[i].ae4*dX0[sl[i].iE4];
		if ((sl[i].iN4>-1) && (fabs(sl[i].an4) > nonzeroEPS)) buf -= sl[i].an4*dX0[sl[i].iN4];
		if ((sl[i].iS4>-1) && (fabs(sl[i].as4) > nonzeroEPS)) buf -= sl[i].as4*dX0[sl[i].iS4];
		if ((sl[i].iT4>-1) && (fabs(sl[i].at4) > nonzeroEPS)) buf -= sl[i].at4*dX0[sl[i].iT4];
		if ((sl[i].iW4>-1) && (fabs(sl[i].aw4) > nonzeroEPS)) buf -= sl[i].aw4*dX0[sl[i].iW4];
		buf*=buf;
		res_sum+=buf;
	}
	for (integer i=0; i<maxbound; ++i) {
		// граничные узлы.
		doublereal buf=0.0;
		buf=slb[i].aw*dX0[slb[i].iW]-dV[slb[i].iW];
		if ((slb[i].iI>-1) && (fabs(slb[i].ai) > nonzeroEPS)) buf-=slb[i].ai*dX0[slb[i].iI];
		buf*=buf;
		res_sum+=buf;
	}
	res_sum=sqrt(res_sum);
	//printf("residual start=%1.4e\n",res_sum);
	//system("pause");
	
	// результаты тестирования
	// задача, начальная невязка , значение евклидовой нормы невязки при которой решение является полученным.
	// tgf01 5.4357e-1 1.0209e-11
	// CGHV1J с метализацией 3.3667e-1 5.0712e-12
	// tgf02 7.6872e-11 1.434e-11
	// tgf05 1.0871e+0  2.2895e-11
	// резистор на 1мм поликоре 5.0e-2 4.9174e-14
	//Diamond ZUb 4 4.0016e-1  4.64444e-11
	// DiamondZUB 4.0016e-1 1.1443e-8
	// NXP100 4.3399e+0  7.8347e-11 (для решения хватило 8Гб ОЗУ.)



	doublereal res_sum_previos = 1.05*finish_residual;
	if (adiabatic_vs_heat_transfer_coeff == NEWTON_RICHMAN_BC) {
		// Работает задача Ньютона Рихмана.
		res_sum_previos = 1.0e-12;
	}


	//if (res_sum>1.0E-10) 
	if (res_sum>res_sum_previos) // защита от повторного холостого запуска экономит время конечного пользователя.
	{

	//yes_print_amg=false;
	yes_print_amg=false;
	 
	

	integer id=0;

	integer ierr=0;
	doublereal eps=1.0e-12;

	ierr=0; // изначальное состояние безошибочное.
	// Порог точности решения СЛАУ. Значение 1.0E-12 достаточно что проверено в ANSYS icepak.
	eps=1.0e-3; // рекомендуемое значение которого достаточно. 

// Требования к оперативной памяти.
/*     VECTOR         NEEDED LENGTH (GUESS) */
/*       A               3*NNA + 5*NNU */
/*       JA              3*NNA + 5*NNU */
/*       IA              2.2*NNU */
/*       U               2.2*NNU */
/*       F               2.2*NNU */
/*       IG              5.4*NNU */

	
	integer nna=0; // количество ненулевых элементов в матрице СЛАУ.
	

	// подсчёт числа ненулевых элементов в матрице.
	nna=0;
	for (integer i=0; i<maxelm; ++i) {
		// внутренность матрицы.
		if ((sl[i].iB>-1) && (fabs(sl[i].ab) > nonzeroEPS)) (nna)++;
		if ((sl[i].iE>-1) && (fabs(sl[i].ae) > nonzeroEPS)) (nna)++;
		if ((sl[i].iN>-1) && (fabs(sl[i].an) > nonzeroEPS)) (nna)++;
		if ((sl[i].iS>-1) && (fabs(sl[i].as) > nonzeroEPS)) (nna)++;
		if ((sl[i].iT>-1) && (fabs(sl[i].at) > nonzeroEPS)) (nna)++;
		if ((sl[i].iW>-1) && (fabs(sl[i].aw) > nonzeroEPS)) (nna)++;
		if ((sl[i].iP>-1) && (fabs(sl[i].ap) > nonzeroEPS)) (nna)++;
		// Дополнение для АЛИС сетки:
		if ((sl[i].iB2>-1) && (fabs(sl[i].ab2) > nonzeroEPS)) (nna)++;
		if ((sl[i].iE2>-1) && (fabs(sl[i].ae2) > nonzeroEPS)) (nna)++;
		if ((sl[i].iN2>-1) && (fabs(sl[i].an2) > nonzeroEPS)) (nna)++;
		if ((sl[i].iS2>-1) && (fabs(sl[i].as2) > nonzeroEPS)) (nna)++;
		if ((sl[i].iT2>-1) && (fabs(sl[i].at2) > nonzeroEPS)) (nna)++;
		if ((sl[i].iW2>-1) && (fabs(sl[i].aw2) > nonzeroEPS)) (nna)++;
		// Дополнение для АЛИС сетки:
		if ((sl[i].iB3>-1) && (fabs(sl[i].ab3) > nonzeroEPS)) (nna)++;
		if ((sl[i].iE3>-1) && (fabs(sl[i].ae3) > nonzeroEPS)) (nna)++;
		if ((sl[i].iN3>-1) && (fabs(sl[i].an3) > nonzeroEPS)) (nna)++;
		if ((sl[i].iS3>-1) && (fabs(sl[i].as3) > nonzeroEPS)) (nna)++;
		if ((sl[i].iT3>-1) && (fabs(sl[i].at3) > nonzeroEPS)) (nna)++;
		if ((sl[i].iW3>-1) && (fabs(sl[i].aw3) > nonzeroEPS)) (nna)++;
		// Дополнение для АЛИС сетки:
		if ((sl[i].iB4>-1) && (fabs(sl[i].ab4) > nonzeroEPS)) (nna)++;
		if ((sl[i].iE4>-1) && (fabs(sl[i].ae4) > nonzeroEPS)) (nna)++;
		if ((sl[i].iN4>-1) && (fabs(sl[i].an4) > nonzeroEPS)) (nna)++;
		if ((sl[i].iS4>-1) && (fabs(sl[i].as4) > nonzeroEPS)) (nna)++;
		if ((sl[i].iT4>-1) && (fabs(sl[i].at4) > nonzeroEPS)) (nna)++;
		if ((sl[i].iW4>-1) && (fabs(sl[i].aw4) > nonzeroEPS)) (nna)++;
	}
	for (integer i=0; i<maxbound; ++i) {
		// граничные узлы.
		if ((slb[i].iW>-1) && (fabs(slb[i].aw) > nonzeroEPS)) (nna)++;
		if ((slb[i].iI>-1) && (fabs(slb[i].ai) > nonzeroEPS)) (nna)++;
	}

	integer nnu=0; // число неизвестных.
	nnu = maxelm + maxbound;

	/*
	// Рекомендуемые по умолчанию параметры.
	integer nda=0; // память под вектор значений матрицы слау.
	nda=3*(nna)+5*(nnu);
	integer ndia=0;
	ndia=static_cast<integer>(2.2*(nnu));
	integer ndja=0;
	ndja=3*(nna)+5*(nnu);
	integer ndu=0;
	ndu=static_cast<integer>(2.2*(nnu));
	integer ndf=0;
	ndf=static_cast<integer>(2.2*(nnu));
	integer ndig=0;
	ndig=static_cast<integer>(5.4*(nnu));
	*/

	/*
	// в двое больше памяти чем рекомендовано.
	integer nda=0; // память под вектор значений матрицы слау.
	nda=6*(nna)+10*(nnu);
	integer ndia=0;
	ndia=static_cast<integer>(4.4*(nnu));
	integer ndja=0;
	ndja=6*(nna)+10*(nnu);
	integer ndu=0;
	ndu=static_cast<integer>(4.4*(nnu));
	integer ndf=0;
	ndf=static_cast<integer>(4.4*(nnu));
	integer ndig=0;
	ndig=static_cast<integer>(10.8*(nnu));
	*/

	// данная константа работоспособна вплоть до размерностей сетки равных 34млн 463тысячи 250узлов.
	//doublereal rsize=1.51; // 1048416
	// Вынужденные течения достаточно 2.5.
	// значения 3.5 недостаточно для 8 модулей Пионер. 
	doublereal rsize=4.5; // на задаче Концевого Ю.А. Электростатика со столбиком в случае сетки со сгущением достаточно 2.0.

	integer nda=0; // память под вектор значений матрицы слау.
	nda=static_cast<integer>(rsize*(3*(nna)+5*(nnu)));
	printf("nda=%lld\n",nda);
	integer ndia=0;
	ndia=static_cast<integer>(rsize*2.2*(nnu));
	integer ndja=0;
	ndja=static_cast<integer>(rsize*(3*(nna)+5*(nnu)));
	integer ndu=0;
	ndu=static_cast<integer>(rsize*2.2*(nnu));
	integer ndf=0;
	ndf=static_cast<integer>(rsize*2.2*(nnu));
	integer ndig=0;
	ndig=static_cast<integer>(rsize*5.4*(nnu));

	/*     CLASS 3 - PARAMETERS: */

/*     LEVELX   -   MAXIMUM NUMBER OF MG-LEVELS TO BE CREATED (>=1). */

/*     IFIRST   -   PARAMETER FOR FIRST APPROXIMATION. */

/*                  1ST DIGIT OF IFIRST: NOT USED; HAS TO BE NON-ZERO. */

/*                  2ND DIGIT OF IFIRST  --  ITYPU: */
/*                    =0: NO SETTING OF FIRST APPROXIMATION, */
/*                    =1: FIRST APPROXIMATION CONSTANT TO ZERO, */
/*                    =2: FIRST APPROXIMATION CONSTANT TO ONE, */
/*                    =3: FIRST APPROXIMATION IS RANDOM FUNCTION WITH */
/*                        THE CONCRETE RANDOM SEQUENCE BEING DETERMINED */
/*                        BY THE FOLLWING DIGITS. */

/*                  REST OF IFIRST  --  RNDU: */
/*                    DETERMINES THE CONCRETE RANDOM SEQUENCE USED IN */
/*                    THE CASE ITYPU=3. (IFIRST=13 IS EQUIVALENT TO */
/*                    IFIRST=1372815) */

/*     NCYC     -   INTEGER PARAMETER DESCRIBING THE TYPE OF CYCLE TO BE */
/*                  USED AND THE NUMBER OF CYCLES TO BE PERFORMED. */

/*                  1ST DIGIT OF NCYC  --  IGAM: */
/*                    =1: V -CYCLE, */
/*                    =2: V*-CYCLE, */
/*                    =3: F -CYCLE, */
/*                    =4: W -CYCLE. */
/*                  IF NCYC IS NEGATIV, THEN THE APPROXIMATION OF THE */
/*                  PROBLEM ON THE SECOND FINEST GRID IS COMPUTED BY */
/*                  IGAM V-CYCLES ON THAT PARTICULAR GRID. */

/*                  2ND DIGIT OF NCYC  --  ICGR: */
/*                    =0: NO CONJUGATE GRADIENT, */
/*                    =1: CONJUGATE GRADIENT (ONLY FIRST STEP OF CG), */
/*                    =2: CONJUGATE GRADIENT (FULL CG). */

/*                  3RD DIGIT OF NCYC  --  ICONV: */
/*                    CONVERGENCE CRITERION FOR THE USER-DEFINED PROBLEM */
/*                    (FINEST GRID): */
/*                    =1: PERFORM A FIXED NUMBER OF CYCLES AS GIVEN BY */
/*                        NCYCLE (SEE BELOW) */
/*                    =2: STOP, IF  ||RES|| < EPS */
/*                    =3: STOP, IF  ||RES|| < EPS * |F| */
/*                    =4: STOP, IF  ||RES|| < EPS * |U| * |DIAG| */
/*                    WITH ||RES|| = L2-NORM OF RESIDUAL, */
/*                           EPS     (SEE INPUT PARAMETER EPS) */
/*                           |F|   = SUPREMUM NORM OF RIGHT HAND SIDE */
/*                           |U|   = SUPREMUM NORM OF SOLUTION */
/*                         |DIAG|  = MAXIMAL DIAGONAL ENTRY IN MATRIX L */
/*                    NOTE THAT IN ANY CASE THE SOLUTION PROCESS STOPS */
/*                    AFTER AT MOST NCYCLE CYCLES. */

/*                  REST OF NCYC  --  NCYCLE: */
/*                    MAXIMAL NUMBER OF CYCLES TO BE PERFORMED (>0) OR */
/*                    NCYCLE=0: NO CYCLING. */

/*     EPS      -   CONVERGENCE CRITERION FOR SOLUTION PROCESS: (SEE */
/*                  PARAMETER NCYC). NOTE THAT NO MORE THAN NCYCLE CYCLES */
/*                  ARE PERFORMED, REGARDLESS OF EPS. */

/*     MADAPT   -   INTEGER VALUE SPECIFYING THE CHOICE OF COARSEST */
/*                  GRID IN CYCLING: */

/*                  1ST DIGIT OF MADAPT  --  MSEL: */
/*                    =1: IN CYCLING, ALL GRIDS CONSTRUCTED IN THE SETUP */
/*                        PHASE ARE USED WITHOUT CHECK. */
/*                    =2: THE NUMBER OF GRIDS IS AUTOMATICALLY REDUCED */
/*                        IF THE CONVERGENCE FACTOR ON THE COARSER GRIDS */
/*                        IS FOUND TO BE LARGER THAN A GIVEN VALUE FAC */
/*                        (SEE BELOW). */

/*                  REST OF MADAPT  --  FAC */
/*                        THE REST OF MADAPT DEFINES THE FRACTIONAL PART */
/*                        OF A REAL NUMBER FAC BETWEEN 0.1 AND 0.99, E.G. */
/*                        MADAPT=258 MEANS MSEL=2 AND FAC=0.58. IF MADAPT */
/*                        CONSISTS OF ONLY ONE DIGIT, FAC IS SET TO 0.7 */
/*                        BY DEFAULT. */


/*     NRD      -   PARAMETER DESCRIBING RELAXATION (DOWNWARDS): */

/*                  1ST DIGIT OF NRD: NOT USED; HAS TO BE NON-ZERO. */

/*                  2ND DIGIT OF NRD  --  NRDX: */
/*                    ACTUAL NUMBER OF SMOOTHING STEPS TO BE PERFORMED */
/*                    THE TYPE OF WHICH IS GIVEN BY THE FOLLOWING DIGITS */

/*                  FOLLOWING DIGITS  --  ARRAY NRDTYP: */
/*                    =1: RELAXATION OVER THE F-POINTS ONLY */
/*                    =2: FULL GS SWEEP */
/*                    =3: RELAXATION OVER THE C-POINTS ONLY */
/*                    =4: FULL MORE COLOR SWEEP, HIGHEST COLOR FIRST */

/*     NSOLCO   -   PARAMETER CONTROLLING THE SOLUTION ON COARSEST GRID: */

/*                  1ST DIGIT  --  NSC: */
/*                    =1: GAUSS-SEIDEL METHOD */
/*                    =2: DIRECT SOLVER (YALE SMP) */

/*                  REST OF NSOLCO  --  NRCX: (ONLY IF NSC=1) */
/*                  NUMBER OF GS SWEEPS ON COARSEST GRID (>=0). */
/*                  IF NRCX=0, THEN AS MANY GS SWEEPS ARE PERFORMED */
/*                  AS ARE NEEDED TO REDUCE THE RESIDUAL BY TWO ORDERS */
/*                  OF MAGNITUDE. (MAXIMAL 100 RELAXATION SWEEPS) */

/*     NRU      -   PARAMETER FOR RELAXATION (UPWARDS), ANALOGOUS TO NRD. */

/*         -------------------------------------------------------------- */

/*     CLASS 4 - PARAMETERS: */

/*     ECG1,ECG2-   REAL PARAMETERS AFFECTING THE CREATION OF COARSER */
/*     EWT2     -   GRIDS AND/OR THE DEFINITION OF THE INTERPOLATION. */
/*                  THE CHOICE OF THESE PARAMETERS DEPENDS ON */
/*                  THE ACTUAL AMG VERSION (SEE SUBROUTINE CRSNG) */

/*     NWT      -   INTEGER PARAMETER AFFECTING THE CREATION OF COARSER */
/*                  GRIDS AND/OR THE DEFINITION OF THE INTERPOLATION. */
/*                  THE CHOICE OF THIS PARAMETER DEPENDS ON */
/*                  THE ACTUAL AMG VERSION (SEE SUBROUTINE CRSNG) */

/*     NTR      -   PARAMETER CONTROLLING COARSE-GRID OPERATOR TRUNCATION */
/*                    =0: PAIRS OF ZEROES ARE REMOVED FROM COARSE GRID */
/*                        OPERATORS */
/*                    =1: NO COARSE-GRID OPERATOR TRUNCATION */



/*     STANDARD CHOICES OF PARAMETERS (AS FAR AS MEANINGFUL): */

/*          ISWTCH = 4 */
/*          IOUT   = 12 */
/*          IPRINT = 10606 */

/*          LEVELX = 25 */
/*          IFIRST = 13 */
/*          NCYC   = 10110 */
/*          EPS    = 1.D-12 */
/*          MADAPT = 27 */
/*          NRD    = 1131 */
/*          NSOLCO = 110 */
/*          NRU    = 1131 */

/*          ECG1   = 0. */
/*          ECG2   = 0.25 */
/*          EWT2   = 0.35 */
/*          NWT    = 2 */
/*          NTR    = 0 */



	// рекомедуемые параметры по дефолту.

	integer iswtch=0;
	iswtch=4;
    integer iout=0;
	iout=13; // 13 обеспечивает печать изменения невязки в процессе счёта.
	integer iprint=0;
	iprint=10606;
	integer levelx=0;
	levelx=25;
	integer ifirst=0;
	// начальное приближение:
	// 0 - используется из вне.
	// 1 - нулевое.
	// 2 - единицы.
	// 3 - случайная последовательность.
	ifirst=13;//13 по умолчанию.
	//ifirst=11; // нулевое начальное приближение.
	//ifirst=10; // вроде как начальное приближение берётся из dX0.
	// но 10 никоим образом не улучшает сходимость.
	integer ncyc=0;
	//ncyc=10110;
	ncyc = 10299; // максимум 99 V циклов
	integer madapt=0;
	madapt=27;
	integer nrd=0;
	nrd=1131;
	integer nsolco=0;
	nsolco=110;
	integer nru=0;
	nru=1131;
	doublereal ecg1=0.0;
	ecg1=0.0;
	doublereal ecg2=0.0;
	ecg2=0.25;
	doublereal ewt2=0.0;
	ewt2=0.35;
	integer nwt=0;
	nwt=2;
	integer ntr=0;
    ntr=0;

	integer matrix=0;
	//matrix=11; // symmetric SPD.
	matrix=22; 

	if ((iVar == TEMP) && (adiabatic_vs_heat_transfer_coeff == NEWTON_RICHMAN_BC)) {
		ifirst = 10;// начальное приближение с предыдущего шага.
		ncyc = 10101; // Всего один V цикл.
		matrix = 11;
	}

	if ((iVar==PAM)&&(bLRfree)) {
		//printf("work amg1r5\n");
		//system("pause");
		// Симметричная положительно определённая матрица это такая матрица
		// которая возникает для поправки давления при решении вязких несжимаемых уравнений Навье-Стокса в 
		// случае задач: каверна, тест Валь-Девиса. Для задач промышленного масштаба это всякие естественные
		// конвекции охлаждающие висящие в воздухе без контакта с теплоотводом греющиеся изделия.
		// Это особый специфический класс задач.
		matrix=11;
	}

	if (iVar!=TEMP)
	{
		//cfd
		//4.08.2018
		if (iVar == PAM) {
			ifirst = 13;
			ncyc = 10110;
			//ncyc = 10399; // максимум 99 V циклов
			//eps = 1.0e-5;
			//5000->10;
			//3720->10
			//255194->
		}
		else {
			ifirst = 10;
			ncyc = 10101;
		}
	}
	//system("pause");
	// allocate memory.
	doublereal *a=NULL;
	//a=new doublereal[nda+1];
	// 15 jan 2016
	a=(doublereal*)malloc((static_cast<integer>(nda) + 1)*sizeof(doublereal));
	if (a==NULL) {
	    // недостаточно памяти на данном оборудовании.
		printf("Problem: not enough memory on your equipment for a matrix in amg1r5 algorithm...\n");
		printf("Please any key to exit...\n");
		//system("pause");
		system("pause");
		exit(1);
	}
	integer *ia=NULL;
	//ia=new integer[ndia+1];
	ia = (integer*)malloc((static_cast<integer>(ndia)+1)*sizeof(integer));
	if (ia==NULL) {
	    // недостаточно памяти на данном оборудовании.
		printf("Problem: not enough memory on your equipment for ia matrix in amg1r5 algorithm...\n");
		printf("Please any key to exit...\n");
		//system("pause");
		system("pause");
		exit(1);
	}
	integer *ja=NULL;
	//ja=new integer[ndja+1];
	ja = (integer*)malloc((static_cast<integer>(ndja)+1)*sizeof(integer));
	if (ja==NULL) {
	    // недостаточно памяти на данном оборудовании.
		printf("Problem: not enough memory on your equipment for ja matrix in amg1r5 algorithm...\n");
		printf("Please any key to exit...\n");
		//system("pause");
		system("pause");
		exit(1);
	}
	doublereal *u=NULL;
	//u = new doublereal[ndu + 1];
	u = (doublereal*)malloc((static_cast<integer>(ndu)+1)*sizeof(doublereal));
	if (u==NULL) {
	    // недостаточно памяти на данном оборудовании.
		printf("Problem: not enough memory on your equipment for u vector in amg1r5 algorithm...\n");
		printf("Please any key to exit...\n");
		//system("pause");
		system("pause");
		exit(1);
	}
	doublereal *f=NULL;
	//f=new doublereal[ndf+1];
	f = (doublereal*)malloc((static_cast<integer>(ndf)+1)*sizeof(doublereal));
	if (f==NULL) {
	    // недостаточно памяти на данном оборудовании.
		printf("Problem: not enough memory on your equipment for f vector in amg1r5 algorithm...\n");
		printf("Please any key to exit...\n");
		//system("pause");
		system("pause");
		exit(1);
	}
	integer *ig=NULL;
	//ig=new integer[ndig+1];
	ig = (integer*)malloc((static_cast<integer>(ndig)+1)*sizeof(integer));
	if (ig==NULL) {
	    // недостаточно памяти на данном оборудовании.
		printf("Problem: not enough memory on your equipment for ig vector in amg1r5 algorithm...\n");
		printf("Please any key to exit...\n");
		//system("pause");
		system("pause");
		exit(1);
	}

	// Блок инициализации нулём, возможно будет работоспособно и без него.

	for (integer k=0; k<=nda; ++k) {
		a[k]=0.0;
	}
	for (integer k=0; k<=ndia; ++k) {
		ia[k]=0;
	}
	for (integer k=0; k<=ndja; ++k) {
		ja[k]=0;
	}
	for (integer k=0; k<=ndu; ++k) {
		u[k]=0.0;
	}
	for (integer k=0; k<=ndf; ++k) {
		f[k]=0.0;
	}
	for (integer k=0; k<=ndig; ++k) {
		ig[k]=0;
	}


	// обязателная инициализация.
	for (integer k=0; k<=nnu+1; ++k) ia[k+id]=nna+1; // инициализация.
	if (id==1) ia[nnu+2]=0;
	

	
	   


	// начальное приближение.
	for (integer i=0; i<=ndu; ++i) {
		u[i]=0.0;
		if (i<maxelm+maxbound) {
			// обязательно нужно проверить была ли выделена оперативная память. 
			u[i+id]=dX0[i]; 
		}
	}

	// правая часть.
    for (integer i=0; i<=ndf; ++i) {
		f[i]=0.0;
		if (i<maxelm+maxbound) {
			// обязательно нужно проверить была ли выделена оперативная память. 
			f[i+id]=dV[i]; 
		}
	}

	// см. equation3DtoCRS.

	    integer ik=0; // счётчик ненулевых элементов СЛАУ
		
		// для внутренних узлов расчётной области:
        for (integer k=0; k<maxelm; ++k) {

			if (fabs(sl[k].ap) > nonzeroEPS) {
                a[ik+id]=sl[k].ap/alpharelax;
				ja[ik+id]=sl[k].iP+1;
                ia[k+id]=my_imin(ik+1,ia[k+id]);
				ik++;
			}
			if ((sl[k].iE>-1) && (fabs(sl[k].ae) > nonzeroEPS)) {
                a[ik+id]=-sl[k].ae;
				ja[ik+id]=sl[k].iE+1;
				ia[k + id] = my_imin(ik + 1, ia[k + id]);
				ik++;
			}
			if ((sl[k].iN>-1) && (fabs(sl[k].an) > nonzeroEPS)) {
                a[ik+id]=-sl[k].an;
				ja[ik+id]=sl[k].iN+1;
				ia[k + id] = my_imin(ik + 1, ia[k + id]);
				ik++;
			}
			if ((sl[k].iT>-1) && (fabs(sl[k].at) > nonzeroEPS)) {
                a[ik+id]=-sl[k].at;
				ja[ik+id]=sl[k].iT+1;
				ia[k + id] = my_imin(ik + 1, ia[k + id]);
				ik++;
			}		
			if ((sl[k].iS>-1) && (fabs(sl[k].as) > nonzeroEPS)) {
                a[ik+id]=-sl[k].as;
				ja[ik+id]=sl[k].iS+1;
				ia[k + id] = my_imin(ik + 1, ia[k + id]);
				ik++;
			}
			if ((sl[k].iW>-1) && (fabs(sl[k].aw) > nonzeroEPS)) {
				a[ik+id]=-sl[k].aw;
				ja[ik+id]=sl[k].iW+1;
				ia[k + id] = my_imin(ik + 1, ia[k + id]);
				ik++;
			}
			if ((sl[k].iB>-1) && (fabs(sl[k].ab) > nonzeroEPS)) {
				a[ik+id]=-sl[k].ab;
				ja[ik+id]=sl[k].iB+1;
				ia[k + id] = my_imin(ik + 1, ia[k + id]);
				ik++;
			}

			// Дополнение на АЛИС сетке:
			if ((sl[k].iE2>-1) && (fabs(sl[k].ae2) > nonzeroEPS)) {
				a[ik + id] = -sl[k].ae2;
				ja[ik + id] = sl[k].iE2 + 1;
				ia[k + id] = my_imin(ik + 1, ia[k + id]);
				ik++;
			}
			if ((sl[k].iN2>-1) && (fabs(sl[k].an2) > nonzeroEPS)) {
				a[ik + id] = -sl[k].an2;
				ja[ik + id] = sl[k].iN2 + 1;
				ia[k + id] = my_imin(ik + 1, ia[k + id]);
				ik++;
			}
			if ((sl[k].iT2>-1) && (fabs(sl[k].at2) > nonzeroEPS)) {
				a[ik + id] = -sl[k].at2;
				ja[ik + id] = sl[k].iT2 + 1;
				ia[k + id] = my_imin(ik + 1, ia[k + id]);
				ik++;
			}
			if ((sl[k].iS2>-1) && (fabs(sl[k].as2) > nonzeroEPS)) {
				a[ik + id] = -sl[k].as2;
				ja[ik + id] = sl[k].iS2 + 1;
				ia[k + id] = my_imin(ik + 1, ia[k + id]);
				ik++;
			}
			if ((sl[k].iW2>-1) && (fabs(sl[k].aw2) > nonzeroEPS)) {
				a[ik + id] = -sl[k].aw2;
				ja[ik + id] = sl[k].iW2 + 1;
				ia[k + id] = my_imin(ik + 1, ia[k + id]);
				ik++;
			}
			if ((sl[k].iB2>-1) && (fabs(sl[k].ab2) > nonzeroEPS)) {
				a[ik + id] = -sl[k].ab2;
				ja[ik + id] = sl[k].iB2 + 1;
				ia[k + id] = my_imin(ik + 1, ia[k + id]);
				ik++;
			}

			// Дополнение на АЛИС сетке:
			if ((sl[k].iE3>-1) && (fabs(sl[k].ae3) > nonzeroEPS)) {
				a[ik + id] = -sl[k].ae3;
				ja[ik + id] = sl[k].iE3 + 1;
				ia[k + id] = my_imin(ik + 1, ia[k + id]);
				ik++;
			}
			if ((sl[k].iN3>-1) && (fabs(sl[k].an3) > nonzeroEPS)) {
				a[ik + id] = -sl[k].an3;
				ja[ik + id] = sl[k].iN3 + 1;
				ia[k + id] = my_imin(ik + 1, ia[k + id]);
				ik++;
			}
			if ((sl[k].iT3>-1) && (fabs(sl[k].at3) > nonzeroEPS)) {
				a[ik + id] = -sl[k].at3;
				ja[ik + id] = sl[k].iT3 + 1;
				ia[k + id] = my_imin(ik + 1, ia[k + id]);
				ik++;
			}
			if ((sl[k].iS3>-1) && (fabs(sl[k].as3) > nonzeroEPS)) {
				a[ik + id] = -sl[k].as3;
				ja[ik + id] = sl[k].iS3 + 1;
				ia[k + id] = my_imin(ik + 1, ia[k + id]);
				ik++;
			}
			if ((sl[k].iW3>-1) && (fabs(sl[k].aw3) > nonzeroEPS)) {
				a[ik + id] = -sl[k].aw3;
				ja[ik + id] = sl[k].iW3 + 1;
				ia[k + id] = my_imin(ik + 1, ia[k + id]);
				ik++;
			}
			if ((sl[k].iB3>-1) && (fabs(sl[k].ab3) > nonzeroEPS)) {
				a[ik + id] = -sl[k].ab3;
				ja[ik + id] = sl[k].iB3 + 1;
				ia[k + id] = my_imin(ik + 1, ia[k + id]);
				ik++;
			}

			// Дополнение на АЛИС сетке:
			if ((sl[k].iE4>-1) && (fabs(sl[k].ae4) > nonzeroEPS)) {
				a[ik + id] = -sl[k].ae4;
				ja[ik + id] = sl[k].iE4 + 1;
				ia[k + id] = my_imin(ik + 1, ia[k + id]);
				ik++;
			}
			if ((sl[k].iN4>-1) && (fabs(sl[k].an4) > nonzeroEPS)) {
				a[ik + id] = -sl[k].an4;
				ja[ik + id] = sl[k].iN4 + 1;
				ia[k + id] = my_imin(ik + 1, ia[k + id]);
				ik++;
			}
			if ((sl[k].iT4>-1) && (fabs(sl[k].at4) > nonzeroEPS)) {
				a[ik + id] = -sl[k].at4;
				ja[ik + id] = sl[k].iT4 + 1;
				ia[k + id] = my_imin(ik + 1, ia[k + id]);
				ik++;
			}
			if ((sl[k].iS4>-1) && (fabs(sl[k].as4) > nonzeroEPS)) {
				a[ik + id] = -sl[k].as4;
				ja[ik + id] = sl[k].iS4 + 1;
				ia[k + id] = my_imin(ik + 1, ia[k + id]);
				ik++;
			}
			if ((sl[k].iW4>-1) && (fabs(sl[k].aw4) > nonzeroEPS)) {
				a[ik + id] = -sl[k].aw4;
				ja[ik + id] = sl[k].iW4 + 1;
				ia[k + id] = my_imin(ik + 1, ia[k + id]);
				ik++;
			}
			if ((sl[k].iB4>-1) && (fabs(sl[k].ab4) > nonzeroEPS)) {
				a[ik + id] = -sl[k].ab4;
				ja[ik + id] = sl[k].iB4 + 1;
				ia[k + id] = my_imin(ik + 1, ia[k + id]);
				ik++;
			}


		}


		// для внутренних узлов расчётной области:
        for (integer k=0; k<maxbound; ++k) {
			if (fabs(slb[k].aw) > nonzeroEPS) {
               // val[ik]=slb[k].aw/alpharelax;
				a[ik+id]=slb[k].aw; // релаксация для граничных узлов не применяется.
				/*if ((slb[k].iI>-1) && (fabs(slb[k].ai) > nonzeroEPS)) {
				     // Внимание !!! было произведено тестирование: один вариант был с нижней релаксацией для граничных узлов,
					 // а второй вариант был без нижней релаксации на граничных узлах. Было выяснено, что для сходимости
					 // более благоприятен вариант без нижней релаксации на граничных узлах.
					 // Данное изменение согласовано с функцией solve.

					 val[ik]/=alpharelax; // Если условия Неймана то нижняя релаксация.
				}*/
				ja[ik+id]=slb[k].iW+1;
				ia[maxelm + k + id] = my_imin(ik + 1, ia[maxelm + k + id]);
				ik++;
			}
			if ((slb[k].iI>-1) && (fabs(slb[k].ai) > nonzeroEPS)) {
				a[ik+id]=-slb[k].ai;
				ja[ik+id]=slb[k].iI+1;
				ia[maxelm + k + id] = my_imin(ik + 1, ia[maxelm + k + id]);
				// Это очень важный вопрос и он требует проверки !
				
				ik++;
			}

		}


		//  
		// нужно акуратно прописать выделения и уничтожения памяти с учётом того что было сделано в BiCGStabP.

        // в каждой строке элементы отсортированы по номерам столбцов:
		// Но диагональный элемент всегда на первом месте в строке матрицы.
		integer imove=0;
		if (id==0) imove=-1;

		// сортировка ненужна порядок следования любой, но главное чтобы первый в строке был имено диагональный элемент.
       //for (integer k=0; k<(maxelm+maxbound); ++k) QuickSortCSIR_amg(ja, a, ia[k+1]+1+imove, ia[k+2]-1+imove); // первый элемент всегда диагональный.
		//for (integer k=0; k<(maxelm+maxbound); ++k) QuickSortCSIR_amg(ja, a, ia[k+1]+imove, ia[k+2]-1+imove); 

		for (integer k=1; k<=nnu; ++k) ig[k+imove]=ia[k+1+imove]; // инициализация.

		
		bool bOkfgmres_amg1r5;
		//printf("getready ...");
		//system("pause");

		if (iVorst_version == 0) {
			// amg - особенно хорош для поправки давления в SIMPLE алгоритме.
			// алгоритм 1985 года.
			amg1r5_(a, ia, ja,
				u, f, ig, &nda, &ndia,
				&ndja, &ndu, &ndf, &ndig,
				&nnu, &matrix, &iswtch, &iout,
				&iprint, &levelx, &ifirst, &ncyc,
				&eps, &madapt, &nrd, &nsolco,
				&nru, &ecg1, &ecg2, &ewt2,
				&nwt, &ntr, &ierr);
		}
		else if (iVorst_version == 1) {

			// 23-24 декабря 2017.

			// В качестве внешнего итерационного процесса используется 
			// алгоритм Хенка Ван Дер Ворста BiCGStab. amg1r5 используется только как
			// многосеточный предобуславливатель.
			amg1r5_Vorst_modification(a, ia, ja,
				u, f, ig, &nda, &ndia,
				&ndja, &ndu, &ndf, &ndig,
				&nnu, &matrix, &iswtch, &iout,
				&iprint, &levelx, &ifirst, &ncyc,
				&eps, &madapt, &nrd, &nsolco,
				&nru, &ecg1, &ecg2, &ewt2,
				&nwt, &ntr, &ierr, iVar, sl, slb, maxelm, maxbound);

		}
		else {
			//31 декабря 2017.

			// В качестве внешнего итерационного процесса используется 
			// алгоритм Ю.Саада и Шульца FGMRes. amg1r5 используется только как
			// многосеточный предобуславливатель.
			amg1r5_fgmres_version(a, ia, ja,
				u, f, ig, &nda, &ndia,
				&ndja, &ndu, &ndf, &ndig,
				&nnu, &matrix, &iswtch, &iout,
				&iprint, &levelx, &ifirst, &ncyc,
				&eps, &madapt, &nrd, &nsolco,
				&nru, &ecg1, &ecg2, &ewt2,
				&nwt, &ntr, &ierr, iVar, sl, slb,
				maxelm, maxbound, bOkfgmres_amg1r5);
		}


		switch (ierr) {
			case 1: printf("dimension A small\n.");
			//system("pause");
				system("pause");
			break;
			case 2: printf("dimension IA small\n.");
			//system("pause");
				system("pause");
			break;
			case 3: printf("dimension JA small\n.");
			//system("pause");
				system("pause");
			break;
			case 4: printf("dimension U small\n.");
			//system("pause");
				system("pause");
			break;
			case 5: printf("dimension F small\n.");
			//system("pause");
				system("pause");
			break;
			case 6: printf("dimension IG small\n.");
			//system("pause");
				system("pause");
			break;
		}

	     // возвращаем решение СЛАУ.
	     for (integer i=0; i<maxelm+maxbound; ++i) {
	        // обратное копирование.
		    dX0[i]=u[i+1+imove]; 
	     }
	

	     // освобождение памяти.
		 if (a!=NULL) {
	       // delete[] a;
			 free(a);
		 }
		 if (ia!=NULL) {
	       // delete[] ia;
			free(ia);
		 }
		 if (ja!=NULL) {
	        //delete[] ja;
			free(ja);
		 }
		 if (u!=NULL) {
			 //delete[] u;
			 free(u);
		 }
		 if (f!=NULL) {
	        //delete[] f;
			 free(f);
		 }
		 if (ig!=NULL) {
	        //delete[] ig;
			free(ig);
		 }

		 
		 res_sum=0.0;
	     for (integer i1=0; i1<maxelm; i1++) {
		     // внутренность матрицы.
		     doublereal buf=0.0;
		     buf=(sl[i1].ap*dX0[sl[i1].iP]-dV[sl[i1].iP]);
		     if ((sl[i1].iB>-1) && (fabs(sl[i1].ab) > nonzeroEPS)) buf-=sl[i1].ab*dX0[sl[i1].iB];
	         if ((sl[i1].iE>-1) && (fabs(sl[i1].ae) > nonzeroEPS)) buf-=sl[i1].ae*dX0[sl[i1].iE];
		     if ((sl[i1].iN>-1) && (fabs(sl[i1].an) > nonzeroEPS)) buf-=sl[i1].an*dX0[sl[i1].iN];
		     if ((sl[i1].iS>-1) && (fabs(sl[i1].as) > nonzeroEPS)) buf-=sl[i1].as*dX0[sl[i1].iS];
		     if ((sl[i1].iT>-1) && (fabs(sl[i1].at) > nonzeroEPS)) buf-=sl[i1].at*dX0[sl[i1].iT];
		     if ((sl[i1].iW>-1) && (fabs(sl[i1].aw) > nonzeroEPS)) buf-=sl[i1].aw*dX0[sl[i1].iW];
			 // дополнение для АЛИС сетки.
			 if ((sl[i1].iB2>-1) && (fabs(sl[i1].ab2) > nonzeroEPS)) buf -= sl[i1].ab2*dX0[sl[i1].iB2];
			 if ((sl[i1].iE2>-1) && (fabs(sl[i1].ae2) > nonzeroEPS)) buf -= sl[i1].ae2*dX0[sl[i1].iE2];
			 if ((sl[i1].iN2>-1) && (fabs(sl[i1].an2) > nonzeroEPS)) buf -= sl[i1].an2*dX0[sl[i1].iN2];
			 if ((sl[i1].iS2>-1) && (fabs(sl[i1].as2) > nonzeroEPS)) buf -= sl[i1].as2*dX0[sl[i1].iS2];
			 if ((sl[i1].iT2>-1) && (fabs(sl[i1].at2) > nonzeroEPS)) buf -= sl[i1].at2*dX0[sl[i1].iT2];
			 if ((sl[i1].iW2>-1) && (fabs(sl[i1].aw2) > nonzeroEPS)) buf -= sl[i1].aw2*dX0[sl[i1].iW2];
			 // дополнение для АЛИС сетки.
			 if ((sl[i1].iB3>-1) && (fabs(sl[i1].ab3) > nonzeroEPS)) buf -= sl[i1].ab3*dX0[sl[i1].iB3];
			 if ((sl[i1].iE3>-1) && (fabs(sl[i1].ae3) > nonzeroEPS)) buf -= sl[i1].ae3*dX0[sl[i1].iE3];
			 if ((sl[i1].iN3>-1) && (fabs(sl[i1].an3) > nonzeroEPS)) buf -= sl[i1].an3*dX0[sl[i1].iN3];
			 if ((sl[i1].iS3>-1) && (fabs(sl[i1].as3) > nonzeroEPS)) buf -= sl[i1].as3*dX0[sl[i1].iS3];
			 if ((sl[i1].iT3>-1) && (fabs(sl[i1].at3) > nonzeroEPS)) buf -= sl[i1].at3*dX0[sl[i1].iT3];
			 if ((sl[i1].iW3>-1) && (fabs(sl[i1].aw3) > nonzeroEPS)) buf -= sl[i1].aw3*dX0[sl[i1].iW3];
			 // дополнение для АЛИС сетки.
			 if ((sl[i1].iB4>-1) && (fabs(sl[i1].ab4) > nonzeroEPS)) buf -= sl[i1].ab4*dX0[sl[i1].iB4];
			 if ((sl[i1].iE4>-1) && (fabs(sl[i1].ae4) > nonzeroEPS)) buf -= sl[i1].ae4*dX0[sl[i1].iE4];
			 if ((sl[i1].iN4>-1) && (fabs(sl[i1].an4) > nonzeroEPS)) buf -= sl[i1].an4*dX0[sl[i1].iN4];
			 if ((sl[i1].iS4>-1) && (fabs(sl[i1].as4) > nonzeroEPS)) buf -= sl[i1].as4*dX0[sl[i1].iS4];
			 if ((sl[i1].iT4>-1) && (fabs(sl[i1].at4) > nonzeroEPS)) buf -= sl[i1].at4*dX0[sl[i1].iT4];
			 if ((sl[i1].iW4>-1) && (fabs(sl[i1].aw4) > nonzeroEPS)) buf -= sl[i1].aw4*dX0[sl[i1].iW4];
	         buf*=buf;
		     res_sum+=buf;
	    }
	    for (integer i1=0; i1<maxbound; i1++) {
	    	// граничные узлы.
		    doublereal buf=0.0;
		    buf=slb[i1].aw*dX0[slb[i1].iW]-dV[slb[i1].iW];
		    if ((slb[i1].iI>-1) && (fabs(slb[i1].ai) > nonzeroEPS)) buf-=slb[i1].ai*dX0[slb[i1].iI];
		    buf*=buf;
		    res_sum+=buf;
	   }
	   res_sum=sqrt(res_sum);
	   //printf("residual finish=%1.4e\n",res_sum);
	   //system("pause");
	   if (bsolid_static_only) {
		   // используется только для теплопередачи в твёрдом теле для ускорения
		   // решения задачи - защита от рестарта.
	       finish_residual=res_sum; // значение невязки решённой задачи.
	   }
	   
		 }

         calculation_main_end_time=clock();
         calculation_vorst_seach_time+=calculation_main_end_time-calculation_main_start_time;

} // amg_loc_memory

  // Здесь содержится обвязка вызывающая amg1r5.
  // локальное выделение памяти:всё внутри, многократные alloc и free.
void amg_loc_memory_for_Matrix_assemble2(SIMPLESPARSE &sparseM, integer n,
	doublereal *dV, doublereal* &dX0,
	integer maxit,
	bool bprintmessage, QuickMemVorst& m)
{

	// Замер времени.
	unsigned int calculation_main_start_time; // начало счёта мс.
	unsigned int calculation_main_end_time; // окончание счёта мс.

	calculation_main_start_time = clock(); // момент начала счёта.

	simplesparsetoCRS(sparseM, m.val, m.col_ind, m.row_ptr, n); // преобразование матрицы из одного формата хранения в другой.

	//simplesparsefree(sparseM, n); // Очистка памяти из под матрицы sparseM.

	//**** apriory matrix check begin ******
	{

	integer i__1 = n, i__;
	integer* icg = new integer[n];
	for (i__ = 0; i__ < i__1; ++i__) {
		icg[i__] = 0;
		// L2: 
	}
	i__1 = n;
	for (i__ = 1; i__ <= i__1; ++i__) {
		i__--;
		integer j1 = m.row_ptr[i__];
		integer j2 = m.row_ptr[i__ + 1] - 1;
		
		
		/*
		// Первый элемент всегда диагональный.
		if (m.col_ind[j1] != i__ ) {
			printf("bug amg1r5 goto L3320\n");
			printf("j1=%d j2=%d, i__=%d ja[j1]=%d\n", j1, j2, i__, m.col_ind[j1]);
			system("pause");
			//goto L3320;
		}
		*/
		integer i__2 = j2;
		for (integer j = j1; j <= i__2; ++j) {
			// Сканируем строку.
			integer i1 = m.col_ind[j];//номер столбца.
			if (i1 < 0 || i1 > n - 1 || icg[i1] == 1) {
				//icg[i1] == 1 - столбец был посещен.
				printf("n=%lld i1=%lld icg=%lld\n", n, i1, icg[i1]);
				if (icg[i1] == 1) {
					printf("POST simplesparsetoCRS in amg_loc_memory_m_ass2 amg two identical columns in a row\n");
					for (integer j69 = j1; j69 <= i__2; ++j69) {
						printf("%lld %lld %lld val=%e col_ind=%lld row_ind=%lld i1=%lld\n", j69, j1, i__2, m.val[j69], m.col_ind[j69], i__, i1);
					}
				}
				system("PAUSE");
			}
			icg[i1] = 1;// столбец был посещен.
						// L5: 
		}
		i__2 = j2;
		for (integer j = j1; j <= i__2; ++j) {
			icg[m.col_ind[j]] = 0;// сброс посещения.
							   // L7: 
		}
		i__++;
	}
	delete[] icg;
}

	//**** apriory matrix check end ******

	// На случай если память не была выделена.
	if (dX0 == NULL) {
		dX0 = new doublereal[n];
		for (integer i = 0; i<n; ++i) {
			dX0[i] = 0.0;
		}
	}

	
	const doublereal nonzeroEPS = 1e-37; // для отделения вещественного нуля
	doublereal res_sum = 0.0;
	res_sum = 0.0;
	res_sum = 1.0;
	
	// результаты тестирования
	// задача, начальная невязка , значение евклидовой нормы невязки при которой решение является полученным.
	// tgf01 5.4357e-1 1.0209e-11
	// CGHV1J с метализацией 3.3667e-1 5.0712e-12
	// tgf02 7.6872e-11 1.434e-11
	// tgf05 1.0871e+0  2.2895e-11
	// резистор на 1мм поликоре 5.0e-2 4.9174e-14
	//Diamond ZUb 4 4.0016e-1  4.64444e-11
	// DiamondZUB 4.0016e-1 1.1443e-8
	// NXP100 4.3399e+0  7.8347e-11 (для решения хватило 8Гб ОЗУ.)



	doublereal res_sum_previos = 1.05*finish_residual;
	if (adiabatic_vs_heat_transfer_coeff == NEWTON_RICHMAN_BC) {
		// Работает задача Ньютона Рихмана.
		res_sum_previos = 1.0e-12;
	}


	//if (res_sum>1.0E-10) 
	if (res_sum>res_sum_previos) // защита от повторного холостого запуска экономит время конечного пользователя.
	{

		//yes_print_amg=false;
		yes_print_amg = true;



		integer id = 0;

		integer ierr = 0;
		doublereal eps = 1.0e-12;

		ierr = 0; // изначальное состояние безошибочное.
				  // Порог точности решения СЛАУ. Значение 1.0E-12 достаточно что проверено в ANSYS icepak.
		eps = 1.0e-3; // рекомендуемое значение которого достаточно. 

					  // Требования к оперативной памяти.
					  /*     VECTOR         NEEDED LENGTH (GUESS) */
					  /*       A               3*NNA + 5*NNU */
					  /*       JA              3*NNA + 5*NNU */
					  /*       IA              2.2*NNU */
					  /*       U               2.2*NNU */
					  /*       F               2.2*NNU */
					  /*       IG              5.4*NNU */


		integer nna = 0; // количество ненулевых элементов в матрице СЛАУ.


						 // подсчёт числа ненулевых элементов в матрице.
		//nna = m.row_ptr[n];
		for (integer k = 0; k < m.row_ptr[n]; ++k) {

			if (fabs(m.val[k]) > nonzeroEPS) {
				nna++;
			}
		}
		
		integer nnu = n; // число неизвестных.
		//nnu = maxelm + maxbound;

		/*
		// Рекомендуемые по умолчанию параметры.
		integer nda=0; // память под вектор значений матрицы слау.
		nda=3*(nna)+5*(nnu);
		integer ndia=0;
		ndia=static_cast<integer>(2.2*(nnu));
		integer ndja=0;
		ndja=3*(nna)+5*(nnu);
		integer ndu=0;
		ndu=static_cast<integer>(2.2*(nnu));
		integer ndf=0;
		ndf=static_cast<integer>(2.2*(nnu));
		integer ndig=0;
		ndig=static_cast<integer>(5.4*(nnu));
		*/

		/*
		// в двое больше памяти чем рекомендовано.
		integer nda=0; // память под вектор значений матрицы слау.
		nda=6*(nna)+10*(nnu);
		integer ndia=0;
		ndia=static_cast<integer>(4.4*(nnu));
		integer ndja=0;
		ndja=6*(nna)+10*(nnu);
		integer ndu=0;
		ndu=static_cast<integer>(4.4*(nnu));
		integer ndf=0;
		ndf=static_cast<integer>(4.4*(nnu));
		integer ndig=0;
		ndig=static_cast<integer>(10.8*(nnu));
		*/

		// данная константа работоспособна вплоть до размерностей сетки равных 34млн 463тысячи 250узлов.
		//doublereal rsize=1.51; // 1048416
		// Вынужденные течения достаточно 2.5.
		// значения 3.5 недостаточно для 8 модулей Пионер. 
		doublereal rsize = 4.5; // на задаче Концевого Ю.А. Электростатика со столбиком в случае сетки со сгущением достаточно 2.0.

		integer nda = 0; // память под вектор значений матрицы слау.
		nda = static_cast<integer>(rsize*(3 * (nna)+5 * (nnu)));
		printf("nda=%lld\n", nda);
		integer ndia = 0;
		ndia = static_cast<integer>(rsize*2.2*(nnu));
		integer ndja = 0;
		ndja = static_cast<integer>(rsize*(3 * (nna)+5 * (nnu)));
		integer ndu = 0;
		ndu = static_cast<integer>(rsize*2.2*(nnu));
		integer ndf = 0;
		ndf = static_cast<integer>(rsize*2.2*(nnu));
		integer ndig = 0;
		ndig = static_cast<integer>(rsize*5.4*(nnu));

		/*     CLASS 3 - PARAMETERS: */

		/*     LEVELX   -   MAXIMUM NUMBER OF MG-LEVELS TO BE CREATED (>=1). */

		/*     IFIRST   -   PARAMETER FOR FIRST APPROXIMATION. */

		/*                  1ST DIGIT OF IFIRST: NOT USED; HAS TO BE NON-ZERO. */

		/*                  2ND DIGIT OF IFIRST  --  ITYPU: */
		/*                    =0: NO SETTING OF FIRST APPROXIMATION, */
		/*                    =1: FIRST APPROXIMATION CONSTANT TO ZERO, */
		/*                    =2: FIRST APPROXIMATION CONSTANT TO ONE, */
		/*                    =3: FIRST APPROXIMATION IS RANDOM FUNCTION WITH */
		/*                        THE CONCRETE RANDOM SEQUENCE BEING DETERMINED */
		/*                        BY THE FOLLWING DIGITS. */

		/*                  REST OF IFIRST  --  RNDU: */
		/*                    DETERMINES THE CONCRETE RANDOM SEQUENCE USED IN */
		/*                    THE CASE ITYPU=3. (IFIRST=13 IS EQUIVALENT TO */
		/*                    IFIRST=1372815) */

		/*     NCYC     -   INTEGER PARAMETER DESCRIBING THE TYPE OF CYCLE TO BE */
		/*                  USED AND THE NUMBER OF CYCLES TO BE PERFORMED. */

		/*                  1ST DIGIT OF NCYC  --  IGAM: */
		/*                    =1: V -CYCLE, */
		/*                    =2: V*-CYCLE, */
		/*                    =3: F -CYCLE, */
		/*                    =4: W -CYCLE. */
		/*                  IF NCYC IS NEGATIV, THEN THE APPROXIMATION OF THE */
		/*                  PROBLEM ON THE SECOND FINEST GRID IS COMPUTED BY */
		/*                  IGAM V-CYCLES ON THAT PARTICULAR GRID. */

		/*                  2ND DIGIT OF NCYC  --  ICGR: */
		/*                    =0: NO CONJUGATE GRADIENT, */
		/*                    =1: CONJUGATE GRADIENT (ONLY FIRST STEP OF CG), */
		/*                    =2: CONJUGATE GRADIENT (FULL CG). */

		/*                  3RD DIGIT OF NCYC  --  ICONV: */
		/*                    CONVERGENCE CRITERION FOR THE USER-DEFINED PROBLEM */
		/*                    (FINEST GRID): */
		/*                    =1: PERFORM A FIXED NUMBER OF CYCLES AS GIVEN BY */
		/*                        NCYCLE (SEE BELOW) */
		/*                    =2: STOP, IF  ||RES|| < EPS */
		/*                    =3: STOP, IF  ||RES|| < EPS * |F| */
		/*                    =4: STOP, IF  ||RES|| < EPS * |U| * |DIAG| */
		/*                    WITH ||RES|| = L2-NORM OF RESIDUAL, */
		/*                           EPS     (SEE INPUT PARAMETER EPS) */
		/*                           |F|   = SUPREMUM NORM OF RIGHT HAND SIDE */
		/*                           |U|   = SUPREMUM NORM OF SOLUTION */
		/*                         |DIAG|  = MAXIMAL DIAGONAL ENTRY IN MATRIX L */
		/*                    NOTE THAT IN ANY CASE THE SOLUTION PROCESS STOPS */
		/*                    AFTER AT MOST NCYCLE CYCLES. */

		/*                  REST OF NCYC  --  NCYCLE: */
		/*                    MAXIMAL NUMBER OF CYCLES TO BE PERFORMED (>0) OR */
		/*                    NCYCLE=0: NO CYCLING. */

		/*     EPS      -   CONVERGENCE CRITERION FOR SOLUTION PROCESS: (SEE */
		/*                  PARAMETER NCYC). NOTE THAT NO MORE THAN NCYCLE CYCLES */
		/*                  ARE PERFORMED, REGARDLESS OF EPS. */

		/*     MADAPT   -   INTEGER VALUE SPECIFYING THE CHOICE OF COARSEST */
		/*                  GRID IN CYCLING: */

		/*                  1ST DIGIT OF MADAPT  --  MSEL: */
		/*                    =1: IN CYCLING, ALL GRIDS CONSTRUCTED IN THE SETUP */
		/*                        PHASE ARE USED WITHOUT CHECK. */
		/*                    =2: THE NUMBER OF GRIDS IS AUTOMATICALLY REDUCED */
		/*                        IF THE CONVERGENCE FACTOR ON THE COARSER GRIDS */
		/*                        IS FOUND TO BE LARGER THAN A GIVEN VALUE FAC */
		/*                        (SEE BELOW). */

		/*                  REST OF MADAPT  --  FAC */
		/*                        THE REST OF MADAPT DEFINES THE FRACTIONAL PART */
		/*                        OF A REAL NUMBER FAC BETWEEN 0.1 AND 0.99, E.G. */
		/*                        MADAPT=258 MEANS MSEL=2 AND FAC=0.58. IF MADAPT */
		/*                        CONSISTS OF ONLY ONE DIGIT, FAC IS SET TO 0.7 */
		/*                        BY DEFAULT. */


		/*     NRD      -   PARAMETER DESCRIBING RELAXATION (DOWNWARDS): */

		/*                  1ST DIGIT OF NRD: NOT USED; HAS TO BE NON-ZERO. */

		/*                  2ND DIGIT OF NRD  --  NRDX: */
		/*                    ACTUAL NUMBER OF SMOOTHING STEPS TO BE PERFORMED */
		/*                    THE TYPE OF WHICH IS GIVEN BY THE FOLLOWING DIGITS */

		/*                  FOLLOWING DIGITS  --  ARRAY NRDTYP: */
		/*                    =1: RELAXATION OVER THE F-POINTS ONLY */
		/*                    =2: FULL GS SWEEP */
		/*                    =3: RELAXATION OVER THE C-POINTS ONLY */
		/*                    =4: FULL MORE COLOR SWEEP, HIGHEST COLOR FIRST */

		/*     NSOLCO   -   PARAMETER CONTROLLING THE SOLUTION ON COARSEST GRID: */

		/*                  1ST DIGIT  --  NSC: */
		/*                    =1: GAUSS-SEIDEL METHOD */
		/*                    =2: DIRECT SOLVER (YALE SMP) */

		/*                  REST OF NSOLCO  --  NRCX: (ONLY IF NSC=1) */
		/*                  NUMBER OF GS SWEEPS ON COARSEST GRID (>=0). */
		/*                  IF NRCX=0, THEN AS MANY GS SWEEPS ARE PERFORMED */
		/*                  AS ARE NEEDED TO REDUCE THE RESIDUAL BY TWO ORDERS */
		/*                  OF MAGNITUDE. (MAXIMAL 100 RELAXATION SWEEPS) */

		/*     NRU      -   PARAMETER FOR RELAXATION (UPWARDS), ANALOGOUS TO NRD. */

		/*         -------------------------------------------------------------- */

		/*     CLASS 4 - PARAMETERS: */

		/*     ECG1,ECG2-   REAL PARAMETERS AFFECTING THE CREATION OF COARSER */
		/*     EWT2     -   GRIDS AND/OR THE DEFINITION OF THE INTERPOLATION. */
		/*                  THE CHOICE OF THESE PARAMETERS DEPENDS ON */
		/*                  THE ACTUAL AMG VERSION (SEE SUBROUTINE CRSNG) */

		/*     NWT      -   INTEGER PARAMETER AFFECTING THE CREATION OF COARSER */
		/*                  GRIDS AND/OR THE DEFINITION OF THE INTERPOLATION. */
		/*                  THE CHOICE OF THIS PARAMETER DEPENDS ON */
		/*                  THE ACTUAL AMG VERSION (SEE SUBROUTINE CRSNG) */

		/*     NTR      -   PARAMETER CONTROLLING COARSE-GRID OPERATOR TRUNCATION */
		/*                    =0: PAIRS OF ZEROES ARE REMOVED FROM COARSE GRID */
		/*                        OPERATORS */
		/*                    =1: NO COARSE-GRID OPERATOR TRUNCATION */



		/*     STANDARD CHOICES OF PARAMETERS (AS FAR AS MEANINGFUL): */

		/*          ISWTCH = 4 */
		/*          IOUT   = 12 */
		/*          IPRINT = 10606 */

		/*          LEVELX = 25 */
		/*          IFIRST = 13 */
		/*          NCYC   = 10110 */
		/*          EPS    = 1.D-12 */
		/*          MADAPT = 27 */
		/*          NRD    = 1131 */
		/*          NSOLCO = 110 */
		/*          NRU    = 1131 */

		/*          ECG1   = 0. */
		/*          ECG2   = 0.25 */
		/*          EWT2   = 0.35 */
		/*          NWT    = 2 */
		/*          NTR    = 0 */



		// рекомедуемые параметры по дефолту.

		integer iswtch = 0;
		iswtch = 4;
		integer iout = 0;
		iout = 13; // 13 обеспечивает печать изменения невязки в процессе счёта.
		integer iprint = 0;
		iprint = 10606;
		integer levelx = 0;
		levelx = 25;
		integer ifirst = 0;
		// начальное приближение:
		// 0 - используется из вне.
		// 1 - нулевое.
		// 2 - единицы.
		// 3 - случайная последовательность.
		ifirst = 13;//13 по умолчанию.
					//ifirst=11; // нулевое начальное приближение.
					//ifirst=10; // вроде как начальное приближение берётся из dX0.
					// но 10 никоим образом не улучшает сходимость.
		integer ncyc = 0;
		//ncyc=10110;
		ncyc = 10299; // максимум 99 V циклов
		integer madapt = 0;
		madapt = 27;
		integer nrd = 0;
		nrd = 1131;
		integer nsolco = 0;
		nsolco = 110;
		integer nru = 0;
		nru = 1131;
		doublereal ecg1 = 0.0;
		ecg1 = 0.0;
		doublereal ecg2 = 0.0;
		ecg2 = 0.25;
		doublereal ewt2 = 0.0;
		ewt2 = 0.35;
		integer nwt = 0;
		nwt = 2;
		integer ntr = 0;
		ntr = 0;

		integer matrix = 0;
		//matrix=11; // symmetric SPD.
		matrix = 22;
		int iVar = TEMP;
		if ((iVar == TEMP) && (adiabatic_vs_heat_transfer_coeff == NEWTON_RICHMAN_BC)) {
			ifirst = 10;// начальное приближение с предыдущего шага.
			ncyc = 10101; // Всего один V цикл.
			matrix = 11;
		}

		

		if (iVar != TEMP)
		{
			//cfd
			//4.08.2018
			if (iVar == PAM) {
				ifirst = 13;
				ncyc = 10110;
				//ncyc = 10399; // максимум 99 V циклов
				//eps = 1.0e-5;
				//5000->10;
				//3720->10
				//255194->
			}
			else {
				ifirst = 10;
				ncyc = 10101;
			}
		}
		//system("pause");
		// allocate memory.
		doublereal *a = NULL;
		//a=new doublereal[nda+1];
		// 15 jan 2016
		a = (doublereal*)malloc((static_cast<integer>(nda)+1) * sizeof(doublereal));
		if (a == NULL) {
			// недостаточно памяти на данном оборудовании.
			printf("Problem: not enough memory on your equipment for a matrix in amg1r5 algorithm...\n");
			printf("Please any key to exit...\n");
			//system("pause");
			system("pause");
			exit(1);
		}
		integer *ia = NULL;
		//ia=new integer[ndia+1];
		ia = (integer*)malloc((static_cast<integer>(ndia)+1) * sizeof(integer));
		if (ia == NULL) {
			// недостаточно памяти на данном оборудовании.
			printf("Problem: not enough memory on your equipment for ia matrix in amg1r5 algorithm...\n");
			printf("Please any key to exit...\n");
			//system("pause");
			system("pause");
			exit(1);
		}
		integer *ja = NULL;
		//ja=new integer[ndja+1];
		ja = (integer*)malloc((static_cast<integer>(ndja)+1) * sizeof(integer));
		if (ja == NULL) {
			// недостаточно памяти на данном оборудовании.
			printf("Problem: not enough memory on your equipment for ja matrix in amg1r5 algorithm...\n");
			printf("Please any key to exit...\n");
			//system("pause");
			system("pause");
			exit(1);
		}
		doublereal *u = NULL;
		//u = new doublereal[ndu + 1];
		u = (doublereal*)malloc((static_cast<integer>(ndu)+1) * sizeof(doublereal));
		if (u == NULL) {
			// недостаточно памяти на данном оборудовании.
			printf("Problem: not enough memory on your equipment for u vector in amg1r5 algorithm...\n");
			printf("Please any key to exit...\n");
			//system("pause");
			system("pause");
			exit(1);
		}
		doublereal *f = NULL;
		//f=new doublereal[ndf+1];
		f = (doublereal*)malloc((static_cast<integer>(ndf)+1) * sizeof(doublereal));
		if (f == NULL) {
			// недостаточно памяти на данном оборудовании.
			printf("Problem: not enough memory on your equipment for f vector in amg1r5 algorithm...\n");
			printf("Please any key to exit...\n");
			//system("pause");
			system("pause");
			exit(1);
		}
		integer *ig = NULL;
		//ig=new integer[ndig+1];
		ig = (integer*)malloc((static_cast<integer>(ndig)+1) * sizeof(integer));
		if (ig == NULL) {
			// недостаточно памяти на данном оборудовании.
			printf("Problem: not enough memory on your equipment for ig vector in amg1r5 algorithm...\n");
			printf("Please any key to exit...\n");
			//system("pause");
			system("pause");
			exit(1);
		}

		// Блок инициализации нулём, возможно будет работоспособно и без него.

		for (integer k = 0; k <= nda; ++k) {
			a[k] = 0.0;
		}
		for (integer k = 0; k <= ndia; ++k) {
			ia[k] = 0;
		}
		for (integer k = 0; k <= ndja; ++k) {
			ja[k] = 0;
		}
		for (integer k = 0; k <= ndu; ++k) {
			u[k] = 0.0;
		}
		for (integer k = 0; k <= ndf; ++k) {
			f[k] = 0.0;
		}
		for (integer k = 0; k <= ndig; ++k) {
			ig[k] = 0;
		}


		// обязателная инициализация.
		for (integer k = 0; k <= nnu + 1; ++k) ia[k + id] = nna+1; // инициализация.//<=
		if (id == 1) ia[nnu + 2] = 0;






		// начальное приближение.
		for (integer i = 0; i <= ndu; ++i) {
			u[i] = 0.0;
			if (i<n) {
				// обязательно нужно проверить была ли выделена оперативная память. 
				u[i + id] = dX0[i];
			}
		}

		// правая часть.
		for (integer i = 0; i <= ndf; ++i) {
			f[i] = 0.0;
			if (i<n) {
				// обязательно нужно проверить была ли выделена оперативная память. 
				f[i + id] = dV[i];
			}
		}

		// см. equation3DtoCRS.

		integer ik = 0; // счётчик ненулевых элементов СЛАУ начинается с нуля.

		{
			integer i__1 = nnu, i__;
			integer* icg = new integer[nnu + 1];
			for ( i__ = 1; i__ <= i__1; ++i__) {
				icg[i__] = 0;
				// L2: 
			}

			// для внутренних узлов расчётной области:
			for (integer k = 0; k < n; ++k) {

				for (integer k1 = m.row_ptr[k]; k1 <= m.row_ptr[k + 1] - 1; k1++) {

					if (m.col_ind[k1] == k) {

						// Диагональ
						if (fabs(m.val[k1]) > nonzeroEPS) {
							// id==0

							a[ik + id] = m.val[k1];
							ja[ik + id] = m.col_ind[k1] + 1;
							ia[k + id] = m.row_ptr[k + id];
							if (icg[m.col_ind[k1]] == 1) {
								printf("error k1=%lld k=%lld\n", k1, k);
								system("PAUSE");
							}
							icg[m.col_ind[k1]] = 1;
							ik++;
						}
					}
				}

				for (integer k1 = m.row_ptr[k]; k1 <= m.row_ptr[k + 1] - 1; k1++) {

					if (m.col_ind[k1] != k) {

						if (fabs(m.val[k1]) > nonzeroEPS) {

							a[ik + id] = m.val[k1];
							ja[ik + id] = m.col_ind[k1] + 1;
							ia[k + id] = m.row_ptr[k + id];
							if (icg[m.col_ind[k1]] == 1) {
								printf("error k1=%lld k=%lld\n",k1,k);
								system("PAUSE");
							}
							icg[m.col_ind[k1]] = 1;
							ik++;
						}
					}
				}

				for (integer k1 = m.row_ptr[k]; k1 <= m.row_ptr[k + 1] - 1; k1++) {
					icg[m.col_ind[k1]] = 0;
				}
				if (0 && k == 11) {
					printf("internal preobraqzovation in amg_loc_memory_for_Matrix_assemble2\n");
					integer ik23 = ik - 1;
					for (integer k1 = m.row_ptr[k]; k1 < m.row_ptr[k + 1]; k1++) {
						printf("val=%e col_ind=%lld row_ind=%lld\n", a[ik23], ja[ik23], k);
						ik23--;
					}
					system("PAUSE");
				}

			}

			if (icg != NULL) {
				delete[] icg;
				icg = NULL;
			}
		}
		for (integer k = 0; k < nnu; ++k) ia[k + id]++;

		//****debug print message*******
		/*
		integer istr_1 = 12;
		printf("POST PRINT bug string\n");
		for (integer ideb = ia[istr_1]; ideb < ia[istr_1 + 1]; ideb++) {
				printf("val=%e col_ind=%d row_ind=%d\n", a[ideb], ja[ideb], istr_1);
		}
		*/
		//for (integer ideb = 0; ideb < 30; ideb++) {
		  // printf("%e %d %d\n", a[ideb], ja[ideb], ia[ideb]);
		//}
		//printf("*********");
		//for (integer ideb = nna-1; ideb >nna- 30; ideb--) {
			//printf("%e %d\n", a[ideb], ja[ideb]);
		//}
		//printf("%d %d\n",nna,ia[n]);
		//system("pause");
		//****debug print message*******
		

		//  
		// нужно акуратно прописать выделения и уничтожения памяти с учётом того что было сделано в BiCGStabP.

		// в каждой строке элементы отсортированы по номерам столбцов:
		// Но диагональный элемент всегда на первом месте в строке матрицы.
		integer imove = 0;
		if (id == 0) imove = -1;

		// сортировка ненужна порядок следования любой, но главное чтобы первый в строке был имено диагональный элемент.
		//for (integer k=0; k<(maxelm+maxbound); ++k) QuickSortCSIR_amg(ja, a, ia[k+1]+1+imove, ia[k+2]-1+imove); // первый элемент всегда диагональный.
		//for (integer k=0; k<(maxelm+maxbound); ++k) QuickSortCSIR_amg(ja, a, ia[k+1]+imove, ia[k+2]-1+imove); 

		for (integer k = 1; k <= nnu; ++k) ig[k + imove] = ia[k + 1 + imove]; // инициализация.

		//**** apriory matrix check begin ******
		
		integer i__1 = nnu, i__;
		integer* icg = new integer[nnu + 1];
		for (i__ = 1; i__ <= i__1; ++i__) {
			icg[i__] = 0;
			// L2: 
		}
		i__1 = nnu;
		for (i__ = 1; i__ <= i__1; ++i__) {
			i__--;
			integer j1 = ia[i__];
			integer j2 = ia[i__ + 1] - 1;
			j1--;
			j2--;
			if (j2 < j1 || j2 > nda) {
				printf("bug amg1r5 goto L3300\n");
				system("PAUSE");
				//goto L3300;
			}
			// Первый элемент всегда диагональный.
			if (ja[j1] != i__+1) {
				printf("bug amg1r5 goto L3320\n");
				printf("j1=%lld j2=%lld, i__=%lld ja[j1]=%lld\n",j1,j2, i__, ja[j1]);
				system("PAUSE");
				//goto L3320;
			}
			integer i__2 = j2;
			for (integer j = j1; j <= i__2; ++j) {
				// Сканируем строку.
				integer i1 = ja[j]-1;//номер столбца.
				if (i1 < 0 || i1 > nnu-1 || icg[i1] == 1) {
					//icg[i1] == 1 - столбец был посещен.
					printf("APRIORY nnu=%lld i1=%lld icg=%lld\n", nnu, i1, icg[i1]);
					if (icg[i1] == 1) {
						printf("APRIORY amg two identical columns in a row\n");
						for (integer j69 = j1; j69 <= i__2; ++j69) {
							printf("%lld %lld %lld val=%e col_ind=%lld row_ind=%lld i1=%lld\n", j69, j1, i__2, a[j69], ja[j69], i__,i1);
						}
					}
					system("PAUSE");
				}
				icg[i1] = 1;// столбец был посещен.
							// L5: 
			}
			i__2 = j2;
			for (integer j = j1; j <= i__2; ++j) {
				icg[ja[j]-1] = 0;// сброс посещения.
							   // L7: 
			}
			i__++;
		}
		if (icg != NULL) {
			delete[] icg;
			icg = NULL;
		}
			
			//**** apriory matrix check end ******


		//printf("getready ...");
		//system("pause");
		if ((AMG1R5_SECOND_T_SOLVER == iswitchsolveramg_vs_BiCGstab_plus_ILU6)&&
		    (NONE_only_amg1r5 == stabilization_amg1r5_algorithm)){
			// amg - особенно хорош для поправки давления в SIMPLE алгоритме.
			// алгоритм 1985 года.
			amg1r5_(a, ia, ja,
				u, f, ig, &nda, &ndia,
				&ndja, &ndu, &ndf, &ndig,
				&nnu, &matrix, &iswtch, &iout,
				&iprint, &levelx, &ifirst, &ncyc,
				&eps, &madapt, &nrd, &nsolco,
				&nru, &ecg1, &ecg2, &ewt2,
				&nwt, &ntr, &ierr);
		}
		else if ((AMG1R5_SECOND_T_SOLVER==iswitchsolveramg_vs_BiCGstab_plus_ILU6)&&
		         (BiCGStab_plus_amg1r5 == stabilization_amg1r5_algorithm)){
			// 23-24 декабря 2017.
			//13.10.2018
			// BiCGStab + amg1r5.

			// В качестве внешнего итерационного процесса используется 
			// алгоритм Хенка Ван Дер Ворста BiCGStab. amg1r5 используется только как
			// многосеточный предобуславливатель.
			amg1r5_Vorst_modification_matrix_Assemble2(a, ia, ja,
				u, f, ig, &nda, &ndia,
				&ndja, &ndu, &ndf, &ndig,
				&nnu, &matrix, &iswtch, &iout,
				&iprint, &levelx, &ifirst, &ncyc,
				&eps, &madapt, &nrd, &nsolco,
				&nru, &ecg1, &ecg2, &ewt2,
				&nwt, &ntr, &ierr,sparseM,n);
		}
		else if ((AMG1R5_SECOND_T_SOLVER==iswitchsolveramg_vs_BiCGstab_plus_ILU6)&&
		(FGMRes_plus_amg1r5 == stabilization_amg1r5_algorithm)){
			// FGMres + amg1r5.
			//31 декабря 2017.

			bool bOkfgmres_amg1r5=false;

			// В качестве внешнего итерационного процесса используется 
			// алгоритм Ю.Саада и Шульца FGMRes. amg1r5 используется только как
			// многосеточный предобуславливатель.
			amg1r5_fgmres_version_matrix_Assemble2(a, ia, ja,
				u, f, ig, &nda, &ndia,
				&ndja, &ndu, &ndf, &ndig,
				&nnu, &matrix, &iswtch, &iout,
				&iprint, &levelx, &ifirst, &ncyc,
				&eps, &madapt, &nrd, &nsolco,
				&nru, &ecg1, &ecg2, &ewt2,
				&nwt, &ntr, &ierr, iVar, sparseM, n, bOkfgmres_amg1r5);
		}

		simplesparsefree(sparseM, n);

		switch (ierr) {
		case 1: printf("dimension A small\n.");
			//system("pause");
			system("pause");
			break;
		case 2: printf("dimension IA small\n.");
			//system("pause");
			system("pause");
			break;
		case 3: printf("dimension JA small\n.");
			//system("pause");
			system("pause");
			break;
		case 4: printf("dimension U small\n.");
			//system("pause");
			system("pause");
			break;
		case 5: printf("dimension F small\n.");
			//system("pause");
			system("pause");
			break;
		case 6: printf("dimension IG small\n.");
			//system("pause");
			system("pause");
			break;
		}

		// возвращаем решение СЛАУ.
		for (integer i = 0; i<n; ++i) {
			// обратное копирование.
			dX0[i] = u[i + 1 + imove];
		}


		// освобождение памяти.
		if (a != NULL) {
			// delete[] a;
			free(a);
		}
		if (ia != NULL) {
			// delete[] ia;
			free(ia);
		}
		if (ja != NULL) {
			//delete[] ja;
			free(ja);
		}
		if (u != NULL) {
			//delete[] u;
			free(u);
		}
		if (f != NULL) {
			//delete[] f;
			free(f);
		}
		if (ig != NULL) {
			//delete[] ig;
			free(ig);
		}

		if (m.val != NULL) {
			delete[] m.val;
			m.val = NULL;
		}
		if (m.col_ind != NULL) {
			delete[] m.col_ind;
			m.col_ind = NULL;
		}
		if (m.row_ptr != NULL) {
			delete[] m.row_ptr;
			m.row_ptr = NULL;
		}

	}

	calculation_main_end_time = clock();
	calculation_vorst_seach_time += calculation_main_end_time - calculation_main_start_time;

} // amg_loc_memory_for_Matrix_assemble2


// глобальная память для amg1r5
typedef struct TamgGlobalMemory {
	doublereal *a;
	integer *ia;
	integer *ja;
	doublereal *u;
	doublereal *f;
	integer *ig;
	integer nda;
	integer ndia;
	integer ndja;
	integer ndu;
	integer ndf;
	integer ndig;
} amgGlobalMemory;

// этот метод показывает значительно более лучшую сходимость, чем простой BiCGStabCRS,
// а также он гораздо лучше (и повидимому правильней) чем Bi_CGStab_internal1.
// Bi_CGStab_internal3 использует предобуславливание из библиотеки Ю.Саада.
// дата написания Bi_CGStab_internal3: 31.03.2013. 
void Bi_CGStab_internal3(equation3D* &sl, equation3D_bon* &slb,
			   integer maxelm, integer maxbound,
			   doublereal *dV, doublereal* &dX0, integer maxit, doublereal alpharelax,
			   bool bprintmessage, integer iVar, QuickMemVorst& m,
	           integer* &ifrontregulationgl, integer* &ibackregulationgl,
	           integer inumber_iteration_SIMPLE);

amgGlobalMemory amgGM;

// Здесь содержится обвязка вызывающая amg1r5.
// Внешняя память, нет выделений и уничтожений памяти.
void amg_global_memory(equation3D* &sl, equation3D_bon* &slb,
			   integer maxelm, integer maxbound,
			   doublereal *dV, doublereal* &dX0, 
			   doublereal alpharelax, integer iVar, bool bLRfree, QuickMemVorst& m,
	           integer* &ifrontregulationgl, integer* &ibackregulationgl,
	           integer iVorst_version)
{


	// iVorst_version == 0 - просто amg1r5 алгоритм.
	// iVorst_version == 1 - amg1r5 алгоритм является предобуславливателем к алгоритму Хенка Ван дер Ворста BiCGStab.

	// Замер времени.
	unsigned int calculation_main_start_time; // начало счёта мс.
	unsigned int calculation_main_end_time; // окончание счёта мс.

	calculation_main_start_time=clock(); // момент начала счёта.



	 // На случай если память не была выделена.
	 if (dX0==NULL) {
	    dX0=new doublereal[maxelm+maxbound];	    
	    for (integer i=0; i<maxelm+maxbound; ++i) {
	        dX0[i]=0.0;
	    }
	 }

	
	const doublereal nonzeroEPS=1e-37; // для отделения вещественного нуля
	doublereal res_sum=0.0;
	res_sum=0.0;
	for (integer i=0; i<maxelm; ++i) {
		// внутренность матрицы.
		doublereal buf=0.0;
		buf=(sl[i].ap*dX0[sl[i].iP]-dV[sl[i].iP]);
		if ((sl[i].iB>-1) && (fabs(sl[i].ab) > nonzeroEPS)) buf-=sl[i].ab*dX0[sl[i].iB];
		if ((sl[i].iE>-1) && (fabs(sl[i].ae) > nonzeroEPS)) buf-=sl[i].ae*dX0[sl[i].iE];
		if ((sl[i].iN>-1) && (fabs(sl[i].an) > nonzeroEPS)) buf-=sl[i].an*dX0[sl[i].iN];
		if ((sl[i].iS>-1) && (fabs(sl[i].as) > nonzeroEPS)) buf-=sl[i].as*dX0[sl[i].iS];
		if ((sl[i].iT>-1) && (fabs(sl[i].at) > nonzeroEPS)) buf-=sl[i].at*dX0[sl[i].iT];
		if ((sl[i].iW>-1) && (fabs(sl[i].aw) > nonzeroEPS)) buf-=sl[i].aw*dX0[sl[i].iW];
		// Дополнение на АЛИС сетке:
		if ((sl[i].iB2>-1) && (fabs(sl[i].ab2) > nonzeroEPS)) buf -= sl[i].ab2*dX0[sl[i].iB2];
		if ((sl[i].iE2>-1) && (fabs(sl[i].ae2) > nonzeroEPS)) buf -= sl[i].ae2*dX0[sl[i].iE2];
		if ((sl[i].iN2>-1) && (fabs(sl[i].an2) > nonzeroEPS)) buf -= sl[i].an2*dX0[sl[i].iN2];
		if ((sl[i].iS2>-1) && (fabs(sl[i].as2) > nonzeroEPS)) buf -= sl[i].as2*dX0[sl[i].iS2];
		if ((sl[i].iT2>-1) && (fabs(sl[i].at2) > nonzeroEPS)) buf -= sl[i].at2*dX0[sl[i].iT2];
		if ((sl[i].iW2>-1) && (fabs(sl[i].aw2) > nonzeroEPS)) buf -= sl[i].aw2*dX0[sl[i].iW2];
		// Дополнение на АЛИС сетке:
		if ((sl[i].iB3>-1) && (fabs(sl[i].ab3) > nonzeroEPS)) buf -= sl[i].ab3*dX0[sl[i].iB3];
		if ((sl[i].iE3>-1) && (fabs(sl[i].ae3) > nonzeroEPS)) buf -= sl[i].ae3*dX0[sl[i].iE3];
		if ((sl[i].iN3>-1) && (fabs(sl[i].an3) > nonzeroEPS)) buf -= sl[i].an3*dX0[sl[i].iN3];
		if ((sl[i].iS3>-1) && (fabs(sl[i].as3) > nonzeroEPS)) buf -= sl[i].as3*dX0[sl[i].iS3];
		if ((sl[i].iT3>-1) && (fabs(sl[i].at3) > nonzeroEPS)) buf -= sl[i].at3*dX0[sl[i].iT3];
		if ((sl[i].iW3>-1) && (fabs(sl[i].aw3) > nonzeroEPS)) buf -= sl[i].aw3*dX0[sl[i].iW3];
		// Дополнение на АЛИС сетке:
		if ((sl[i].iB4>-1) && (fabs(sl[i].ab4) > nonzeroEPS)) buf -= sl[i].ab4*dX0[sl[i].iB4];
		if ((sl[i].iE4>-1) && (fabs(sl[i].ae4) > nonzeroEPS)) buf -= sl[i].ae4*dX0[sl[i].iE4];
		if ((sl[i].iN4>-1) && (fabs(sl[i].an4) > nonzeroEPS)) buf -= sl[i].an4*dX0[sl[i].iN4];
		if ((sl[i].iS4>-1) && (fabs(sl[i].as4) > nonzeroEPS)) buf -= sl[i].as4*dX0[sl[i].iS4];
		if ((sl[i].iT4>-1) && (fabs(sl[i].at4) > nonzeroEPS)) buf -= sl[i].at4*dX0[sl[i].iT4];
		if ((sl[i].iW4>-1) && (fabs(sl[i].aw4) > nonzeroEPS)) buf -= sl[i].aw4*dX0[sl[i].iW4];
		buf*=buf;
		res_sum+=buf;
	}
	for (integer i=0; i<maxbound; ++i) {
		// граничные узлы.
		doublereal buf=0.0;
		buf=slb[i].aw*dX0[slb[i].iW]-dV[slb[i].iW];
		if ((slb[i].iI>-1) && (fabs(slb[i].ai) > nonzeroEPS)) buf-=slb[i].ai*dX0[slb[i].iI];
		buf*=buf;
		res_sum+=buf;
	}
	res_sum=sqrt(res_sum);
	//printf("residual start=%1.4e\n",res_sum);
	//system("pause");
	
	doublereal res_sum_start=res_sum;


	// результаты тестирования
	// задача, начальная невязка , значение евклидовой нормы невязки при которой решение является полученным.
	// tgf01 5.4357e-1 1.0209e-11
	// CGHV1J с метализацией 3.3667e-1 5.0712e-12
	// tgf02 7.6872e-11 1.434e-11
	// tgf05 1.0871e+0  2.2895e-11
	// резистор на 1мм поликоре 5.0e-2 4.9174e-14
	//Diamond ZUb 4 4.0016e-1  4.64444e-11
	// DiamondZUB 4.0016e-1 1.1443e-8
	// NXP100 4.3399e+0  7.8347e-11 (для решения хватило 8Гб ОЗУ.)

	doublereal res_sum_previos = 1.05*finish_residual;
	if (adiabatic_vs_heat_transfer_coeff == NEWTON_RICHMAN_BC) {
		// Работает задача Ньютона Рихмана.
		res_sum_previos = 1.0e-12;
	}

	//if (res_sum>1.0E-10) 
	if (res_sum>res_sum_previos) // защита от повторного холостого запуска экономит время конечного пользователя.
	{

		integer iprogon=0; // в случае расходимости мы будем повторно производить решение.

		//LabelAMGdivergenceDetected:

		LabelReallocMemory:


	//yes_print_amg=false;
		if (bsolid_static_only) {
			yes_print_amg = true;
		}
		else
		{
			yes_print_amg = false;
		}
	 
	

	integer id=0;

	integer ierr=0;
	doublereal eps=1.0e-12;

	ierr=0; // изначальное состояние безошибочное.
	// Порог точности решения СЛАУ. Значение 1.0E-12 достаточно что проверено в ANSYS icepak.
	eps=1.0e-3; // рекомендуемое значение которого достаточно. 

// Требования к оперативной памяти.
/*     VECTOR         NEEDED LENGTH (GUESS) */
/*       A               3*NNA + 5*NNU */
/*       JA              3*NNA + 5*NNU */
/*       IA              2.2*NNU */
/*       U               2.2*NNU */
/*       F               2.2*NNU */
/*       IG              5.4*NNU */

	
	integer nna=0; // количество ненулевых элементов в матрице СЛАУ.
	

	// подсчёт числа ненулевых элементов в матрице.
	nna=0;
	for (integer i=0; i<maxelm; ++i) {
		// внутренность матрицы.
		if ((sl[i].iB>-1) && (fabs(sl[i].ab) > nonzeroEPS)) (nna)++;
		if ((sl[i].iE>-1) && (fabs(sl[i].ae) > nonzeroEPS)) (nna)++;
		if ((sl[i].iN>-1) && (fabs(sl[i].an) > nonzeroEPS)) (nna)++;
		if ((sl[i].iS>-1) && (fabs(sl[i].as) > nonzeroEPS)) (nna)++;
		if ((sl[i].iT>-1) && (fabs(sl[i].at) > nonzeroEPS)) (nna)++;
		if ((sl[i].iW>-1) && (fabs(sl[i].aw) > nonzeroEPS)) (nna)++;
		if ((sl[i].iP>-1) && (fabs(sl[i].ap) > nonzeroEPS)) (nna)++;
		// Дополнение для АЛИС сетки:
		if ((sl[i].iB2>-1) && (fabs(sl[i].ab2) > nonzeroEPS)) (nna)++;
		if ((sl[i].iE2>-1) && (fabs(sl[i].ae2) > nonzeroEPS)) (nna)++;
		if ((sl[i].iN2>-1) && (fabs(sl[i].an2) > nonzeroEPS)) (nna)++;
		if ((sl[i].iS2>-1) && (fabs(sl[i].as2) > nonzeroEPS)) (nna)++;
		if ((sl[i].iT2>-1) && (fabs(sl[i].at2) > nonzeroEPS)) (nna)++;
		if ((sl[i].iW2>-1) && (fabs(sl[i].aw2) > nonzeroEPS)) (nna)++;
		// Дополнение для АЛИС сетки:
		if ((sl[i].iB3>-1) && (fabs(sl[i].ab3) > nonzeroEPS)) (nna)++;
		if ((sl[i].iE3>-1) && (fabs(sl[i].ae3) > nonzeroEPS)) (nna)++;
		if ((sl[i].iN3>-1) && (fabs(sl[i].an3) > nonzeroEPS)) (nna)++;
		if ((sl[i].iS3>-1) && (fabs(sl[i].as3) > nonzeroEPS)) (nna)++;
		if ((sl[i].iT3>-1) && (fabs(sl[i].at3) > nonzeroEPS)) (nna)++;
		if ((sl[i].iW3>-1) && (fabs(sl[i].aw3) > nonzeroEPS)) (nna)++;
		// Дополнение для АЛИС сетки:
		if ((sl[i].iB4>-1) && (fabs(sl[i].ab4) > nonzeroEPS)) (nna)++;
		if ((sl[i].iE4>-1) && (fabs(sl[i].ae4) > nonzeroEPS)) (nna)++;
		if ((sl[i].iN4>-1) && (fabs(sl[i].an4) > nonzeroEPS)) (nna)++;
		if ((sl[i].iS4>-1) && (fabs(sl[i].as4) > nonzeroEPS)) (nna)++;
		if ((sl[i].iT4>-1) && (fabs(sl[i].at4) > nonzeroEPS)) (nna)++;
		if ((sl[i].iW4>-1) && (fabs(sl[i].aw4) > nonzeroEPS)) (nna)++;
	}
	for (integer i=0; i<maxbound; ++i) {
		// граничные узлы.
		if ((slb[i].iW>-1) && (fabs(slb[i].aw) > nonzeroEPS)) (nna)++;
		if ((slb[i].iI>-1) && (fabs(slb[i].ai) > nonzeroEPS)) (nna)++;
	}

	integer nnu=0; // число неизвестных.
	nnu=maxelm+maxbound;

	/*
	// Рекомендуемые по умолчанию параметры.
	integer nda=0; // память под вектор значений матрицы слау.
	nda=3*(nna)+5*(nnu);
	integer ndia=0;
	ndia=static_cast<integer>(2.2*(nnu));
	integer ndja=0;
	ndja=3*(nna)+5*(nnu);
	integer ndu=0;
	ndu=static_cast<integer>(2.2*(nnu));
	integer ndf=0;
	ndf=static_cast<integer>(2.2*(nnu));
	integer ndig=0;
	ndig=static_cast<integer>(5.4*(nnu));
	*/

	/*
	// в двое больше памяти чем рекомендовано.
	integer nda=0; // память под вектор значений матрицы слау.
	nda=6*(nna)+10*(nnu);
	integer ndia=0;
	ndia=static_cast<integer>(4.4*(nnu));
	integer ndja=0;
	ndja=6*(nna)+10*(nnu);
	integer ndu=0;
	ndu=static_cast<integer>(4.4*(nnu));
	integer ndf=0;
	ndf=static_cast<integer>(4.4*(nnu));
	integer ndig=0;
	ndig=static_cast<integer>(10.8*(nnu));
	*/

	// данная константа работоспособна вплоть до размерностей сетки равных 34млн 463тысячи 250узлов.
	//doublereal rsize=1.51; // 1048416
	// Вынужденные течения достаточно 2.5.
	// 3.5 недостаточно для 8 модулей ПИОНЕР.
	doublereal rsize=3.5; // на задаче Концевого Ю.А. Электростатика со столбиком в случае сетки со сгущением достаточно 2.0.
	// 7.05.2017 Opening тест сообщил о нехватке памяти.
	if ((bSIMPLErun_now_for_temperature) && (bSIMPLErun_now_for_natural_convection)) {
		// 8.5 подходит.
		rsize = 8.5; //6.5 неподходит. 10.0 подходит.
	}

	integer nda=0; // память под вектор значений матрицы слау.
	nda=static_cast<integer>(rsize*(3*(nna)+5*(nnu)));
	printf("nda=%lld nna=%lld nnu=%lld maxelm=%lld maxbound=%lld \n",nda, nna, nnu, maxelm, maxbound);
	integer ndia=0;
	ndia=static_cast<integer>(rsize*2.2*(nnu));//2.2
	integer ndja=0;
	ndja=static_cast<integer>(rsize*(3*(nna)+5*(nnu)));
	integer ndu=0;
	ndu=static_cast<integer>(rsize*2.2*(nnu));
	integer ndf=0;
	ndf=static_cast<integer>(rsize*2.2*(nnu));
	integer ndig=0;
	ndig=static_cast<integer>(rsize*5.4*(nnu));

	/*     CLASS 3 - PARAMETERS: */

/*     LEVELX   -   MAXIMUM NUMBER OF MG-LEVELS TO BE CREATED (>=1). */

/*     IFIRST   -   PARAMETER FOR FIRST APPROXIMATION. */

/*                  1ST DIGIT OF IFIRST: NOT USED; HAS TO BE NON-ZERO. */

/*                  2ND DIGIT OF IFIRST  --  ITYPU: */
/*                    =0: NO SETTING OF FIRST APPROXIMATION, */
/*                    =1: FIRST APPROXIMATION CONSTANT TO ZERO, */
/*                    =2: FIRST APPROXIMATION CONSTANT TO ONE, */
/*                    =3: FIRST APPROXIMATION IS RANDOM FUNCTION WITH */
/*                        THE CONCRETE RANDOM SEQUENCE BEING DETERMINED */
/*                        BY THE FOLLWING DIGITS. */

/*                  REST OF IFIRST  --  RNDU: */
/*                    DETERMINES THE CONCRETE RANDOM SEQUENCE USED IN */
/*                    THE CASE ITYPU=3. (IFIRST=13 IS EQUIVALENT TO */
/*                    IFIRST=1372815) */

/*     NCYC     -   INTEGER PARAMETER DESCRIBING THE TYPE OF CYCLE TO BE */
/*                  USED AND THE NUMBER OF CYCLES TO BE PERFORMED. */

/*                  1ST DIGIT OF NCYC  --  IGAM: */
/*                    =1: V -CYCLE, */
/*                    =2: V*-CYCLE, */
/*                    =3: F -CYCLE, */
/*                    =4: W -CYCLE. */
/*                  IF NCYC IS NEGATIV, THEN THE APPROXIMATION OF THE */
/*                  PROBLEM ON THE SECOND FINEST GRID IS COMPUTED BY */
/*                  IGAM V-CYCLES ON THAT PARTICULAR GRID. */

/*                  2ND DIGIT OF NCYC  --  ICGR: */
/*                    =0: NO CONJUGATE GRADIENT, */
/*                    =1: CONJUGATE GRADIENT (ONLY FIRST STEP OF CG), */
/*                    =2: CONJUGATE GRADIENT (FULL CG). */

/*                  3RD DIGIT OF NCYC  --  ICONV: */
/*                    CONVERGENCE CRITERION FOR THE USER-DEFINED PROBLEM */
/*                    (FINEST GRID): */
/*                    =1: PERFORM A FIXED NUMBER OF CYCLES AS GIVEN BY */
/*                        NCYCLE (SEE BELOW) */
/*                    =2: STOP, IF  ||RES|| < EPS */
/*                    =3: STOP, IF  ||RES|| < EPS * |F| */
/*                    =4: STOP, IF  ||RES|| < EPS * |U| * |DIAG| */
/*                    WITH ||RES|| = L2-NORM OF RESIDUAL, */
/*                           EPS     (SEE INPUT PARAMETER EPS) */
/*                           |F|   = SUPREMUM NORM OF RIGHT HAND SIDE */
/*                           |U|   = SUPREMUM NORM OF SOLUTION */
/*                         |DIAG|  = MAXIMAL DIAGONAL ENTRY IN MATRIX L */
/*                    NOTE THAT IN ANY CASE THE SOLUTION PROCESS STOPS */
/*                    AFTER AT MOST NCYCLE CYCLES. */

/*                  REST OF NCYC  --  NCYCLE: */
/*                    MAXIMAL NUMBER OF CYCLES TO BE PERFORMED (>0) OR */
/*                    NCYCLE=0: NO CYCLING. */

/*     EPS      -   CONVERGENCE CRITERION FOR SOLUTION PROCESS: (SEE */
/*                  PARAMETER NCYC). NOTE THAT NO MORE THAN NCYCLE CYCLES */
/*                  ARE PERFORMED, REGARDLESS OF EPS. */

/*     MADAPT   -   INTEGER VALUE SPECIFYING THE CHOICE OF COARSEST */
/*                  GRID IN CYCLING: */

/*                  1ST DIGIT OF MADAPT  --  MSEL: */
/*                    =1: IN CYCLING, ALL GRIDS CONSTRUCTED IN THE SETUP */
/*                        PHASE ARE USED WITHOUT CHECK. */
/*                    =2: THE NUMBER OF GRIDS IS AUTOMATICALLY REDUCED */
/*                        IF THE CONVERGENCE FACTOR ON THE COARSER GRIDS */
/*                        IS FOUND TO BE LARGER THAN A GIVEN VALUE FAC */
/*                        (SEE BELOW). */

/*                  REST OF MADAPT  --  FAC */
/*                        THE REST OF MADAPT DEFINES THE FRACTIONAL PART */
/*                        OF A REAL NUMBER FAC BETWEEN 0.1 AND 0.99, E.G. */
/*                        MADAPT=258 MEANS MSEL=2 AND FAC=0.58. IF MADAPT */
/*                        CONSISTS OF ONLY ONE DIGIT, FAC IS SET TO 0.7 */
/*                        BY DEFAULT. */


/*     NRD      -   PARAMETER DESCRIBING RELAXATION (DOWNWARDS): */

/*                  1ST DIGIT OF NRD: NOT USED; HAS TO BE NON-ZERO. */

/*                  2ND DIGIT OF NRD  --  NRDX: */
/*                    ACTUAL NUMBER OF SMOOTHING STEPS TO BE PERFORMED */
/*                    THE TYPE OF WHICH IS GIVEN BY THE FOLLOWING DIGITS */

/*                  FOLLOWING DIGITS  --  ARRAY NRDTYP: */
/*                    =1: RELAXATION OVER THE F-POINTS ONLY */
/*                    =2: FULL GS SWEEP */
/*                    =3: RELAXATION OVER THE C-POINTS ONLY */
/*                    =4: FULL MORE COLOR SWEEP, HIGHEST COLOR FIRST */

/*     NSOLCO   -   PARAMETER CONTROLLING THE SOLUTION ON COARSEST GRID: */

/*                  1ST DIGIT  --  NSC: */
/*                    =1: GAUSS-SEIDEL METHOD */
/*                    =2: DIRECT SOLVER (YALE SMP) */

/*                  REST OF NSOLCO  --  NRCX: (ONLY IF NSC=1) */
/*                  NUMBER OF GS SWEEPS ON COARSEST GRID (>=0). */
/*                  IF NRCX=0, THEN AS MANY GS SWEEPS ARE PERFORMED */
/*                  AS ARE NEEDED TO REDUCE THE RESIDUAL BY TWO ORDERS */
/*                  OF MAGNITUDE. (MAXIMAL 100 RELAXATION SWEEPS) */

/*     NRU      -   PARAMETER FOR RELAXATION (UPWARDS), ANALOGOUS TO NRD. */

/*         -------------------------------------------------------------- */

/*     CLASS 4 - PARAMETERS: */

/*     ECG1,ECG2-   REAL PARAMETERS AFFECTING THE CREATION OF COARSER */
/*     EWT2     -   GRIDS AND/OR THE DEFINITION OF THE INTERPOLATION. */
/*                  THE CHOICE OF THESE PARAMETERS DEPENDS ON */
/*                  THE ACTUAL AMG VERSION (SEE SUBROUTINE CRSNG) */

/*     NWT      -   INTEGER PARAMETER AFFECTING THE CREATION OF COARSER */
/*                  GRIDS AND/OR THE DEFINITION OF THE INTERPOLATION. */
/*                  THE CHOICE OF THIS PARAMETER DEPENDS ON */
/*                  THE ACTUAL AMG VERSION (SEE SUBROUTINE CRSNG) */

/*     NTR      -   PARAMETER CONTROLLING COARSE-GRID OPERATOR TRUNCATION */
/*                    =0: PAIRS OF ZEROES ARE REMOVED FROM COARSE GRID */
/*                        OPERATORS */
/*                    =1: NO COARSE-GRID OPERATOR TRUNCATION */



/*     STANDARD CHOICES OF PARAMETERS (AS FAR AS MEANINGFUL): */

/*          ISWTCH = 4 */
/*          IOUT   = 12 */
/*          IPRINT = 10606 */

/*          LEVELX = 25 */
/*          IFIRST = 13 */
/*          NCYC   = 10110 */
/*          EPS    = 1.D-12 */
/*          MADAPT = 27 */
/*          NRD    = 1131 */
/*          NSOLCO = 110 */
/*          NRU    = 1131 */

/*          ECG1   = 0. */
/*          ECG2   = 0.25 */
/*          EWT2   = 0.35 */
/*          NWT    = 2 */
/*          NTR    = 0 */



	// рекомедуемые параметры по дефолту.

	integer iswtch=0;
	iswtch=4;
    integer iout=0;
	iout=13; // 13 обеспечивают печать изменения невязки в процессе счёта.
	integer iprint=0;
	iprint=10606;
	integer levelx=0;
	levelx=25;
	integer ifirst=0;
	// начальное приближение:
	// 0 - используется из вне.
	// 1 - нулевое.
	// 2 - единицы.
	// 3 - случайная последовательность.
	ifirst=13;//13 по умолчанию.
	//ifirst=11; // нулевое начальное приближение.
	//ifirst=10; // вроде как начальное приближение берётся из dX0.
	// но 10 никоим образом не улучшает сходимость.
	integer ncyc=0;
	//ncyc=10110; // 10 V циклов
	ncyc = 10299; // максимум 99 V циклов
	if ((iVar == TEMP) && (starting_speed_Vx*starting_speed_Vx + starting_speed_Vy*starting_speed_Vy + starting_speed_Vz*starting_speed_Vz > 1.0e-30)) {
		ncyc = 10120; // 20 / 99 V циклов.
	}
	if ((iVar == TEMP) && (breakRUMBAcalc_for_nonlinear_boundary_condition  )) {
		ncyc = 1012; // 2 - V цикла.
		eps = res_sum*sqrt(static_cast<doublereal>(maxelm + maxbound))*0.0002;
	}
	if (bSIMPLErun_now_for_temperature) {
		//ncyc = 102; // 2 - V цикла.
		//eps = res_sum*sqrt(maxelm + maxbound)*0.001;
	}
	integer madapt=0;
	madapt=27;
	integer nrd=0;
	nrd=1131;
	integer nsolco=0;
	nsolco=110;
	if ((iVar == TEMP) && (breakRUMBAcalc_for_nonlinear_boundary_condition  )) {
		//nsolco = 122;
	}
	integer nru=0;
	nru=1131;
	doublereal ecg1=0.0;
	ecg1=0.0;
	doublereal ecg2=0.0;
	ecg2=0.25;
	doublereal ewt2=0.0;
	ewt2=0.35;
	integer nwt=0;
	nwt=2;
	integer ntr=0;
    ntr=0;

	integer matrix=0;
	//matrix=11; // symmetric SPD.
	matrix=22; 


	// Я так и недобился работоспособности amg1r5 решателя на задаче 
	// с условием Ньютона-Рихмана на всей дефолтной границе. Пока с этой задачей 
	// справляется только my_cl_agl_amg_v0_14 (РУМБА алгоритм).

	if ((iVar==TEMP)&& (adiabatic_vs_heat_transfer_coeff == NEWTON_RICHMAN_BC)) {
		ifirst = 13;// начальное приближение с предыдущего шага.
		ncyc = 1019; // Всего 9 V циклов.
		//if (blocker_Newton_Richman) {
			matrix = 21;
		//}
	}

	if ((iVar==PAM)&&(bLRfree)) {
		//printf("work amg1r5\n");
		//system("pause");
		// Симметричная положительно определённая матрица это такая матрица
		// которая возникает для поправки давления при решении вязких несжимаемых уравнений Навье-Стокса в 
		// случае задач: каверна, тест Валь-Девиса. Для задач промышленного масштаба это всякие естественные
		// конвекции охлаждающие висящие в воздухе без контакта с теплоотводом греющиеся изделия.
		// Это особый специфический класс задач.
		matrix=11;
	}

	if ((iVar==VX)||(iVar==VY)||(iVar==VZ)) {
		// A: Используем начальное приближение.
		// 0 - используется из вне.
	    // 1 - нулевое.
	    // 2 - единицы.
	    // 3 - случайная последовательность.
	    ifirst=10;//13 по умолчанию.
		// B: считаем с точностью до eps.
		// точность на два порядка меньше чем начальное значение невязки.
		//eps=0.01*res_sum*res_sum;
		eps=0.1*res_sum*res_sum;
		// W цикл, Зейдель, с точностью до eps, всего 10 W циклов если надо.
		//ncyc=40210;
		//ncyc=30210;
		//ncyc = 10110; // 10 фиксированных V циклов
		ncyc = 12110; // 10 фиксированных V циклов congruate gradient
	}

	if (0&&(iVar == TEMP) && (bSIMPLErun_now_for_temperature)&&(bSIMPLErun_now_for_natural_convection)) {
		// Моделирование естественной конвекции в рамках SIMPLE алгоритма.
		ifirst = 10;
		eps = 0.0001*res_sum*res_sum;
		//ncyc = 12120; // 20 фиксированных V циклов congruate gradient
		ncyc = 10220; // V seidel 0.1 не более 20 циклов.
	}

	if ((iVar==PAM)) {
		// лучше стартовать с приближения которое было на прошлой итерации так как оно наиболее близкое.
		switch (iprogon) {
			case 0: ifirst=10; break; // начальное приближение с предыдущей итерации.
			case 1: ifirst=11; break; // нулевое начальное приближение.
			case 2: ifirst=12; break; // единичное начальное приближение.
			case 3: ifirst=13; break; // случайное начальное приближение.
			default: exit(1); break;
		}
		
		// порог точности взят из BicgStab+ilu2.
		// Для поправки давления слау должна быть решена с большой точностью.
		if (nnu<100000) {
			if (1.0e-5*fabs(res_sum)<3.0e-8) { //1.0e-3
				eps=sqrt(fabs(1.0e-5*fabs(res_sum))); //1.0e-3
			}
			else {
				eps=sqrt(3.0e-8); // 3.0e-6
			}
		}
		else {
			if (1.0e-6*fabs(res_sum)<3.0e-8) { //1.0e-4
				eps=sqrt(fabs(1.0e-6*fabs(res_sum)));
			}
			else {
				eps=sqrt(3.0e-8);
			}
		}
		//eps=0.00001*res_sum*res_sum;
		//eps = 0.1*res_sum*res_sum;
		//ncyc=12140;
		ncyc = 10140; // фиксированное число V циклов на eps плюём. 40 V циклов.
	}

	// allocate memory.
	if (amgGM.a==NULL) {
	    amgGM.a=new doublereal[nda+1];
		amgGM.nda=nda;
	    if (amgGM.a==NULL) {
	        // недостаточно памяти на данном оборудовании.
		    printf("Problem: not enough memory on your equipment for a matrix in amg1r5 algorithm...\n");
		    printf("Please any key to exit...\n");
		    //system("pause");
			system("pause");
		    exit(1);
	    }
	}
	else {
		if (amgGM.nda<nda) {
			// realloc
			printf("\nPlease wait realloc in amg1r5 for matrix a... ... ...\n");
			delete amgGM.a;
            amgGM.a=NULL;
			amgGM.a=new doublereal[nda+1];
		    amgGM.nda=nda;
	        if (amgGM.a==NULL) {
	            // недостаточно памяти на данном оборудовании.
		        printf("Problem: not enough memory on your equipment for a matrix in amg1r5 algorithm...\n");
		        printf("Please any key to exit...\n");
		        //system("pause");
				system("pause");
		        exit(1);
	        }
		}
	}
	if (amgGM.ia==NULL) {
	   amgGM.ia=new integer[ndia+1];
	   amgGM.ndia=ndia;
	   if (amgGM.ia==NULL) {
	       // недостаточно памяти на данном оборудовании.
		   printf("Problem: not enough memory on your equipment for ia matrix in amg1r5 algorithm...\n");
		   printf("Please any key to exit...\n");
		   //system("pause");
		   system("pause");
		   exit(1);
	   }
	}
	else
	{
		if (amgGM.ndia<ndia) {
		// realloc
			printf("\nPlease wait realloc in amg1r5 for matrix ia... ... ...\n");
			delete amgGM.ia;
            amgGM.ia=NULL;
			amgGM.ia=new integer[ndia+1];
	        amgGM.ndia=ndia;
	        if (amgGM.ia==NULL) {
	           // недостаточно памяти на данном оборудовании.
		       printf("Problem: not enough memory on your equipment for ia matrix in amg1r5 algorithm...\n");
		       printf("Please any key to exit...\n");
		       //system("pause");
			   system("pause");
		       exit(1);
	        }
		}
	}
	if (amgGM.ja==NULL) {
	    amgGM.ja=new integer[ndja+1];
		amgGM.ndja=ndja;
	    if (amgGM.ja==NULL) {
	        // недостаточно памяти на данном оборудовании.
		    printf("Problem: not enough memory on your equipment for ja matrix in amg1r5 algorithm...\n");
		    printf("Please any key to exit...\n");
		    //system("pause");
			system("pause");
		    exit(1);
	    }
	}
	{
		if (amgGM.ndja<ndja) {
		// realloc
			printf("\nPlease wait realloc in amg1r5 for matrix ja... ... ...\n");
			delete amgGM.ja;
            amgGM.ja=NULL;
            amgGM.ja=new integer[ndja+1];
		    amgGM.ndja=ndja;
	        if (amgGM.ja==NULL) {
	             // недостаточно памяти на данном оборудовании.
		         printf("Problem: not enough memory on your equipment for ja matrix in amg1r5 algorithm...\n");
		         printf("Please any key to exit...\n");
		         //system("pause");
				 system("pause");
		         exit(1);
	        }
		}
	}
	if (amgGM.u==NULL) {
	    amgGM.u=new doublereal[ndu+1];
		amgGM.ndu=ndu;
	    if (amgGM.u==NULL) {
	        // недостаточно памяти на данном оборудовании.
		    printf("Problem: not enough memory on your equipment for u vector in amg1r5 algorithm...\n");
	    	printf("Please any key to exit...\n");
		    //system("pause");
			system("pause");
		    exit(1);
	    }
	}
	else {
		if (amgGM.ndu<ndu) {
		// realloc
			printf("\nPlease wait realloc in amg1r5 for matrix u... ... ...\n");
			delete amgGM.u;
            amgGM.u=NULL;
			amgGM.u=new doublereal[ndu+1];
		    amgGM.ndu=ndu;
	        if (amgGM.u==NULL) {
	             // недостаточно памяти на данном оборудовании.
		         printf("Problem: not enough memory on your equipment for u vector in amg1r5 algorithm...\n");
	    	     printf("Please any key to exit...\n");
		        // system("pause");
				 system("pause");
		         exit(1);
	        }
		}

	}
	if (amgGM.f==NULL) {
	    amgGM.f=new doublereal[ndf+1];
		amgGM.ndf=ndf;
	    if (amgGM.f==NULL) {
	        // недостаточно памяти на данном оборудовании.
		    printf("Problem: not enough memory on your equipment for f vector in amg1r5 algorithm...\n");
		    printf("Please any key to exit...\n");
		    //system("pause");
			system("pause");
		    exit(1);
	    }
	}
	else {
		if (amgGM.ndf<ndf) {
			printf("\nPlease wait realloc in amg1r5 for matrix f... ... ...\n");
			delete amgGM.f;
            amgGM.f=NULL;
			amgGM.f=new doublereal[ndf+1];
		    amgGM.ndf=ndf;
	        if (amgGM.f==NULL) {
	            // недостаточно памяти на данном оборудовании.
		        printf("Problem: not enough memory on your equipment for f vector in amg1r5 algorithm...\n");
		        printf("Please any key to exit...\n");
		        //system("pause");
				system("pause");
		        exit(1);
	        }
		}
	}
	if (amgGM.ig==NULL) {
	    amgGM.ig=new integer[ndig+1];
		amgGM.ndig=ndig;
	    if (amgGM.ig==NULL) {
	         // недостаточно памяти на данном оборудовании.
		     printf("Problem: not enough memory on your equipment for ig vector in amg1r5 algorithm...\n");
		     printf("Please any key to exit...\n");
		     //system("pause");
			 system("pause");
		     exit(1);
	    }
	}
	else {
		if (amgGM.ndig<ndig) {
			printf("\nPlease wait realloc in amg1r5 for matrix ig... ... ...\n");
			delete amgGM.ig;
            amgGM.ig=NULL;
			amgGM.ig=new integer[ndig+1];
		    amgGM.ndig=ndig;
	        if (amgGM.ig==NULL) {
	            // недостаточно памяти на данном оборудовании.
		        printf("Problem: not enough memory on your equipment for ig vector in amg1r5 algorithm...\n");
		        printf("Please any key to exit...\n");
		       // system("pause");
				system("pause");
		        exit(1);
	         }
		}
	}

	

	// Блок инициализации нулём, возможно будет работоспособно и без него.

	for (integer k=0; k<=nda; ++k) {
		amgGM.a[k]=0.0;
	}
	for (integer k=0; k<=ndia; ++k) {
		amgGM.ia[k]=0;
	}
	for (integer k=0; k<=ndja; ++k) {
		amgGM.ja[k]=0;
	}
	for (integer k=0; k<=ndu; ++k) {
		amgGM.u[k]=0.0;
	}
	for (integer k=0; k<=ndf; ++k) {
		amgGM.f[k]=0.0;
	}
	for (integer k=0; k<=ndig; ++k) {
		amgGM.ig[k]=0;
	}


	// обязателная инициализация.
	for (integer k=0; k<=nnu+1; ++k) amgGM.ia[k+id]=nna+1; // инициализация.
	if (id==1) amgGM.ia[nnu+2]=0;
	

	
	   


	// начальное приближение.
	for (integer i=0; i<=ndu; ++i) {
		amgGM.u[i]=0.0;
		if (i<maxelm+maxbound) {
			// обязательно нужно проверить была ли выделена оперативная память. 
			amgGM.u[i+id]=dX0[i]; 
		}
	}

	// правая часть.
    for (integer i=0; i<=ndf; ++i) {
		amgGM.f[i]=0.0;
		if (i<maxelm+maxbound) {
			// обязательно нужно проверить была ли выделена оперативная память. 
			amgGM.f[i+id]=dV[i]; 
		}
	}

	// см. equation3DtoCRS.

	    integer ik=0; // счётчик ненулевых элементов СЛАУ
		
		// для внутренних узлов расчётной области:
        for (integer k=0; k<maxelm; ++k) {

			if (fabs(sl[k].ap) > nonzeroEPS) {
				if (iprogon==1) {
					amgGM.a[ik+id]=1.05*sl[k].ap/alpharelax;
				}
				else if (iprogon==2) {
					amgGM.a[ik+id]=1.1*sl[k].ap/alpharelax;
				}
				else {
                    amgGM.a[ik+id]=sl[k].ap/alpharelax;
				}
				amgGM.ja[ik+id]=sl[k].iP+1;
				amgGM.ia[k + id] = my_imin(ik + 1, amgGM.ia[k + id]);
				ik++;
			}
			if ((sl[k].iE>-1) && (fabs(sl[k].ae) > nonzeroEPS)) {
                amgGM.a[ik+id]=-sl[k].ae;
				amgGM.ja[ik+id]=sl[k].iE+1;
				amgGM.ia[k + id] = my_imin(ik + 1, amgGM.ia[k + id]);
				ik++;
			}
			if ((sl[k].iN>-1) && (fabs(sl[k].an) > nonzeroEPS)) {
                amgGM.a[ik+id]=-sl[k].an;
				amgGM.ja[ik+id]=sl[k].iN+1;
				amgGM.ia[k + id] = my_imin(ik + 1, amgGM.ia[k + id]);
				ik++;
			}
			if ((sl[k].iT>-1) && (fabs(sl[k].at) > nonzeroEPS)) {
                amgGM.a[ik+id]=-sl[k].at;
				amgGM.ja[ik+id]=sl[k].iT+1;
				amgGM.ia[k + id] = my_imin(ik + 1, amgGM.ia[k + id]);
				ik++;
			}		
			if ((sl[k].iS>-1) && (fabs(sl[k].as) > nonzeroEPS)) {
                amgGM.a[ik+id]=-sl[k].as;
				amgGM.ja[ik+id]=sl[k].iS+1;
				amgGM.ia[k + id] = my_imin(ik + 1, amgGM.ia[k + id]);
				ik++;
			}
			if ((sl[k].iW>-1) && (fabs(sl[k].aw) > nonzeroEPS)) {
				amgGM.a[ik+id]=-sl[k].aw;
				amgGM.ja[ik+id]=sl[k].iW+1;
				amgGM.ia[k + id] = my_imin(ik + 1, amgGM.ia[k + id]);
				ik++;
			}
			if ((sl[k].iB>-1) && (fabs(sl[k].ab) > nonzeroEPS)) {
				amgGM.a[ik+id]=-sl[k].ab;
				amgGM.ja[ik+id]=sl[k].iB+1;
				amgGM.ia[k + id] = my_imin(ik + 1, amgGM.ia[k + id]);
				ik++;
			}

			if ((sl[k].iE2>-1) && (fabs(sl[k].ae2) > nonzeroEPS)) {
				amgGM.a[ik + id] = -sl[k].ae2;
				amgGM.ja[ik + id] = sl[k].iE2 + 1;
				amgGM.ia[k + id] = my_imin(ik + 1, amgGM.ia[k + id]);
				ik++;
			}
			if ((sl[k].iN2>-1) && (fabs(sl[k].an2) > nonzeroEPS)) {
				amgGM.a[ik + id] = -sl[k].an2;
				amgGM.ja[ik + id] = sl[k].iN2 + 1;
				amgGM.ia[k + id] = my_imin(ik + 1, amgGM.ia[k + id]);
				ik++;
			}
			if ((sl[k].iT2>-1) && (fabs(sl[k].at2) > nonzeroEPS)) {
				amgGM.a[ik + id] = -sl[k].at2;
				amgGM.ja[ik + id] = sl[k].iT2 + 1;
				amgGM.ia[k + id] = my_imin(ik + 1, amgGM.ia[k + id]);
				ik++;
			}
			if ((sl[k].iS2>-1) && (fabs(sl[k].as2) > nonzeroEPS)) {
				amgGM.a[ik + id] = -sl[k].as2;
				amgGM.ja[ik + id] = sl[k].iS2 + 1;
				amgGM.ia[k + id] = my_imin(ik + 1, amgGM.ia[k + id]);
				ik++;
			}
			if ((sl[k].iW2>-1) && (fabs(sl[k].aw2) > nonzeroEPS)) {
				amgGM.a[ik + id] = -sl[k].aw2;
				amgGM.ja[ik + id] = sl[k].iW2 + 1;
				amgGM.ia[k + id] = my_imin(ik + 1, amgGM.ia[k + id]);
				ik++;
			}
			if ((sl[k].iB2>-1) && (fabs(sl[k].ab2) > nonzeroEPS)) {
				amgGM.a[ik + id] = -sl[k].ab2;
				amgGM.ja[ik + id] = sl[k].iB2 + 1;
				amgGM.ia[k + id] = my_imin(ik + 1, amgGM.ia[k + id]);
				ik++;
			}

			if ((sl[k].iE3>-1) && (fabs(sl[k].ae3) > nonzeroEPS)) {
				amgGM.a[ik + id] = -sl[k].ae3;
				amgGM.ja[ik + id] = sl[k].iE3 + 1;
				amgGM.ia[k + id] = my_imin(ik + 1, amgGM.ia[k + id]);
				ik++;
			}
			if ((sl[k].iN3>-1) && (fabs(sl[k].an3) > nonzeroEPS)) {
				amgGM.a[ik + id] = -sl[k].an3;
				amgGM.ja[ik + id] = sl[k].iN3 + 1;
				amgGM.ia[k + id] = my_imin(ik + 1, amgGM.ia[k + id]);
				ik++;
			}
			if ((sl[k].iT3>-1) && (fabs(sl[k].at3) > nonzeroEPS)) {
				amgGM.a[ik + id] = -sl[k].at3;
				amgGM.ja[ik + id] = sl[k].iT3 + 1;
				amgGM.ia[k + id] = my_imin(ik + 1, amgGM.ia[k + id]);
				ik++;
			}
			if ((sl[k].iS3>-1) && (fabs(sl[k].as3) > nonzeroEPS)) {
				amgGM.a[ik + id] = -sl[k].as3;
				amgGM.ja[ik + id] = sl[k].iS3 + 1;
				amgGM.ia[k + id] = my_imin(ik + 1, amgGM.ia[k + id]);
				ik++;
			}
			if ((sl[k].iW3>-1) && (fabs(sl[k].aw3) > nonzeroEPS)) {
				amgGM.a[ik + id] = -sl[k].aw3;
				amgGM.ja[ik + id] = sl[k].iW3 + 1;
				amgGM.ia[k + id] = my_imin(ik + 1, amgGM.ia[k + id]);
				ik++;
			}
			if ((sl[k].iB3>-1) && (fabs(sl[k].ab3) > nonzeroEPS)) {
				amgGM.a[ik + id] = -sl[k].ab3;
				amgGM.ja[ik + id] = sl[k].iB3 + 1;
				amgGM.ia[k + id] = my_imin(ik + 1, amgGM.ia[k + id]);
				ik++;
			}

			if ((sl[k].iE4>-1) && (fabs(sl[k].ae4) > nonzeroEPS)) {
				amgGM.a[ik + id] = -sl[k].ae4;
				amgGM.ja[ik + id] = sl[k].iE4 + 1;
				amgGM.ia[k + id] = my_imin(ik + 1, amgGM.ia[k + id]);
				ik++;
			}
			if ((sl[k].iN4>-1) && (fabs(sl[k].an4) > nonzeroEPS)) {
				amgGM.a[ik + id] = -sl[k].an4;
				amgGM.ja[ik + id] = sl[k].iN4 + 1;
				amgGM.ia[k + id] = my_imin(ik + 1, amgGM.ia[k + id]);
				ik++;
			}
			if ((sl[k].iT4>-1) && (fabs(sl[k].at4) > nonzeroEPS)) {
				amgGM.a[ik + id] = -sl[k].at4;
				amgGM.ja[ik + id] = sl[k].iT4 + 1;
				amgGM.ia[k + id] = my_imin(ik + 1, amgGM.ia[k + id]);
				ik++;
			}
			if ((sl[k].iS4>-1) && (fabs(sl[k].as4) > nonzeroEPS)) {
				amgGM.a[ik + id] = -sl[k].as4;
				amgGM.ja[ik + id] = sl[k].iS4 + 1;
				amgGM.ia[k + id] = my_imin(ik + 1, amgGM.ia[k + id]);
				ik++;
			}
			if ((sl[k].iW4>-1) && (fabs(sl[k].aw4) > nonzeroEPS)) {
				amgGM.a[ik + id] = -sl[k].aw4;
				amgGM.ja[ik + id] = sl[k].iW4 + 1;
				amgGM.ia[k + id] = my_imin(ik + 1, amgGM.ia[k + id]);
				ik++;
			}
			if ((sl[k].iB4>-1) && (fabs(sl[k].ab4) > nonzeroEPS)) {
				amgGM.a[ik + id] = -sl[k].ab4;
				amgGM.ja[ik + id] = sl[k].iB4 + 1;
				amgGM.ia[k + id] = my_imin(ik + 1, amgGM.ia[k + id]);
				ik++;
			}

		}


		// для внутренних узлов расчётной области:
        for (integer k=0; k<maxbound; ++k) {
			if (fabs(slb[k].aw) > nonzeroEPS) {
               // val[ik]=slb[k].aw/alpharelax;
				amgGM.a[ik+id]=slb[k].aw; // релаксация для граничных узлов не применяется.
				/*if ((slb[k].iI>-1) && (fabs(slb[k].ai) > nonzeroEPS)) {
				     // Внимание !!! было произведено тестирование: один вариант был с нижней релаксацией для граничных узлов,
					 // а второй вариант был без нижней релаксации на граничных узлах. Было выяснено, что для сходимости
					 // более благоприятен вариант без нижней релаксации на граничных узлах.
					 // Данное изменение согласовано с функцией solve.

					 val[ik]/=alpharelax; // Если условия Неймана то нижняя релаксация.
				}*/
				amgGM.ja[ik+id]=slb[k].iW+1;
				amgGM.ia[maxelm + k + id] = my_imin(ik + 1, amgGM.ia[maxelm + k + id]);
				ik++;
			}
			if ((slb[k].iI>-1) && (fabs(slb[k].ai) > nonzeroEPS)) {
				amgGM.a[ik+id]=-slb[k].ai;
				amgGM.ja[ik+id]=slb[k].iI+1;
				amgGM.ia[maxelm + k + id] = my_imin(ik + 1, amgGM.ia[maxelm + k + id]);
				// Это очень важный вопрос и он требует проверки !
				
				ik++;
			}

		}

		// debug*****
		//for (integer ideb = 1; ideb < 30; ideb++) {
			//printf("%e %d %d\n", amgGM.a[ideb], amgGM.ja[ideb], amgGM.ia[ideb]);
		//}
		//printf("%d %d\n", nna, amgGM.ia[nnu]);
		//system("pause");
		//debug 13.10.2018*******
		//printf("ik=%d\n",ik);
		//system("pause");

		//  
		// нужно акуратно прописать выделения и уничтожения памяти с учётом того что было сделано в BiCGStabP.

        // в каждой строке элементы отсортированы по номерам столбцов:
		// Но диагональный элемент всегда на первом месте в строке матрицы.
		integer imove=0;
		if (id==0) imove=-1;

		// сортировка ненужна порядок следования любой, но главное чтобы первый в строке был имено диагональный элемент.
       //for (integer k=0; k<(maxelm+maxbound); ++k) QuickSortCSIR_amg(ja, a, ia[k+1]+1+imove, ia[k+2]-1+imove); // первый элемент всегда диагональный.
		//for (integer k=0; k<(maxelm+maxbound); ++k) QuickSortCSIR_amg(ja, a, ia[k+1]+imove, ia[k+2]-1+imove); 

		for (integer k=1; k<=nnu; ++k) amgGM.ig[k+imove]=amgGM.ia[k+1+imove]; // инициализация.

		bool bOkfgmres_amg1r5;

		//printf("getready ...");
		//system("pause");

	    // amg - особенно хорош для поправки давления в SIMPLE алгоритме.
	    // алгоритм Руге и Штубена 1985 года.
		if (iVorst_version == 0) {
			// просто amg1r5 алгоритм.

			amg1r5_(amgGM.a, amgGM.ia, amgGM.ja,
				amgGM.u, amgGM.f, amgGM.ig, &nda, &ndia,
				&ndja, &ndu, &ndf, &ndig,
				&nnu, &matrix, &iswtch, &iout,
				&iprint, &levelx, &ifirst, &ncyc,
				&eps, &madapt, &nrd, &nsolco,
				&nru, &ecg1, &ecg2, &ewt2,
				&nwt, &ntr, &ierr);
		}
		else if (iVorst_version == 1) {

			// 23-24 декабря 2017.

			// В качестве внешнего итерационного процесса используется 
			// алгоритм Хенка Ван Дер Ворста BiCGStab. amg1r5 используется только как
			// многосеточный предобуславливатель.
			amg1r5_Vorst_modification(amgGM.a, amgGM.ia, amgGM.ja,
				amgGM.u, amgGM.f, amgGM.ig, &nda, &ndia,
				&ndja, &ndu, &ndf, &ndig,
				&nnu, &matrix, &iswtch, &iout,
				&iprint, &levelx, &ifirst, &ncyc,
				&eps, &madapt, &nrd, &nsolco,
				&nru, &ecg1, &ecg2, &ewt2,
				&nwt, &ntr, &ierr, iVar, sl, slb, maxelm, maxbound);

		}
		else {
			//31 декабря 2017.

			// В качестве внешнего итерационного процесса используется 
			// алгоритм Ю.Саада и Шульца FGMRes. amg1r5 используется только как
			// многосеточный предобуславливатель.
			amg1r5_fgmres_version(amgGM.a, amgGM.ia, amgGM.ja,
				amgGM.u, amgGM.f, amgGM.ig, &nda, &ndia,
				&ndja, &ndu, &ndf, &ndig,
				&nnu, &matrix, &iswtch, &iout,
				&iprint, &levelx, &ifirst, &ncyc,
				&eps, &madapt, &nrd, &nsolco,
				&nru, &ecg1, &ecg2, &ewt2,
				&nwt, &ntr, &ierr, iVar, sl, slb,
				maxelm, maxbound, bOkfgmres_amg1r5);
		}



		switch (ierr) {
			case 1: printf("dimension A small\n.");
			//system("pause");
				system("pause");
			break;
			case 2: printf("dimension IA small\n.");
			//system("pause");
				system("pause");
			break;
			case 3: printf("dimension JA small\n.");
			//system("pause");
				system("pause");
			break;
			case 4: printf("dimension U small\n.");
			//system("pause");
				system("pause");
			break;
			case 5: printf("dimension F small\n.");
			//system("pause");
				system("pause");
			break;
			case 6: printf("dimension IG small\n.");
			//system("pause");
				system("pause");
			break;
		}

		if ((ierr>0) && (ierr<7)) {
			printf("please, WAIT...  ... ...");
            // освобождение памяти из под amg1r5.
	        if (amgGM.a!=NULL) {
	           delete amgGM.a;
               amgGM.a=NULL;
	        }
	        if (amgGM.ia!=NULL) {
	           delete amgGM.ia;
	           amgGM.ia=NULL;
	        }
	        if (amgGM.ja!=NULL) {
	           delete amgGM.ja;
		       amgGM.ja=NULL;
	        }
	        if (amgGM.u!=NULL) {
	           delete amgGM.u;
		       amgGM.u=NULL;
	        }
	        if (amgGM.f!=NULL) {
	           delete amgGM.f;
		       amgGM.f=NULL;
	        }
	        if (amgGM.ig!=NULL) {
	           delete amgGM.ig;
		       amgGM.ig=NULL;
	        }

	        // повторное выделение большегокуска памяти.
			goto LabelReallocMemory;
		}

	     // возвращаем решение СЛАУ.
	     for (integer i=0; i<maxelm+maxbound; ++i) {
	        // обратное копирование.
		    dX0[i]=amgGM.u[i+1+imove]; 
	     }
	

	    // Освобождение памяти осуществляется единожды на глобальном уровне видимости.
		 
		 
		 res_sum=0.0;
	     for (integer i1=0; i1<maxelm; i1++) {
		     // внутренность матрицы.
		     doublereal buf=0.0;
		     buf=(sl[i1].ap*dX0[sl[i1].iP]-dV[sl[i1].iP]);
		     if ((sl[i1].iB>-1) && (fabs(sl[i1].ab) > nonzeroEPS)) buf-=sl[i1].ab*dX0[sl[i1].iB];
	         if ((sl[i1].iE>-1) && (fabs(sl[i1].ae) > nonzeroEPS)) buf-=sl[i1].ae*dX0[sl[i1].iE];
		     if ((sl[i1].iN>-1) && (fabs(sl[i1].an) > nonzeroEPS)) buf-=sl[i1].an*dX0[sl[i1].iN];
		     if ((sl[i1].iS>-1) && (fabs(sl[i1].as) > nonzeroEPS)) buf-=sl[i1].as*dX0[sl[i1].iS];
		     if ((sl[i1].iT>-1) && (fabs(sl[i1].at) > nonzeroEPS)) buf-=sl[i1].at*dX0[sl[i1].iT];
		     if ((sl[i1].iW>-1) && (fabs(sl[i1].aw) > nonzeroEPS)) buf-=sl[i1].aw*dX0[sl[i1].iW];
			 // С учётом АЛИС сетки:
			 if ((sl[i1].iB2>-1) && (fabs(sl[i1].ab2) > nonzeroEPS)) buf -= sl[i1].ab2*dX0[sl[i1].iB2];
			 if ((sl[i1].iE2>-1) && (fabs(sl[i1].ae2) > nonzeroEPS)) buf -= sl[i1].ae2*dX0[sl[i1].iE2];
			 if ((sl[i1].iN2>-1) && (fabs(sl[i1].an2) > nonzeroEPS)) buf -= sl[i1].an2*dX0[sl[i1].iN2];
			 if ((sl[i1].iS2>-1) && (fabs(sl[i1].as2) > nonzeroEPS)) buf -= sl[i1].as2*dX0[sl[i1].iS2];
			 if ((sl[i1].iT2>-1) && (fabs(sl[i1].at2) > nonzeroEPS)) buf -= sl[i1].at2*dX0[sl[i1].iT2];
			 if ((sl[i1].iW2>-1) && (fabs(sl[i1].aw2) > nonzeroEPS)) buf -= sl[i1].aw2*dX0[sl[i1].iW2];
			 // С учётом АЛИС сетки:
			 if ((sl[i1].iB3>-1) && (fabs(sl[i1].ab3) > nonzeroEPS)) buf -= sl[i1].ab3*dX0[sl[i1].iB3];
			 if ((sl[i1].iE3>-1) && (fabs(sl[i1].ae3) > nonzeroEPS)) buf -= sl[i1].ae3*dX0[sl[i1].iE3];
			 if ((sl[i1].iN3>-1) && (fabs(sl[i1].an3) > nonzeroEPS)) buf -= sl[i1].an3*dX0[sl[i1].iN3];
			 if ((sl[i1].iS3>-1) && (fabs(sl[i1].as3) > nonzeroEPS)) buf -= sl[i1].as3*dX0[sl[i1].iS3];
			 if ((sl[i1].iT3>-1) && (fabs(sl[i1].at3) > nonzeroEPS)) buf -= sl[i1].at3*dX0[sl[i1].iT3];
			 if ((sl[i1].iW3>-1) && (fabs(sl[i1].aw3) > nonzeroEPS)) buf -= sl[i1].aw3*dX0[sl[i1].iW3];
			 // С учётом АЛИС сетки:
			 if ((sl[i1].iB4>-1) && (fabs(sl[i1].ab4) > nonzeroEPS)) buf -= sl[i1].ab4*dX0[sl[i1].iB4];
			 if ((sl[i1].iE4>-1) && (fabs(sl[i1].ae4) > nonzeroEPS)) buf -= sl[i1].ae4*dX0[sl[i1].iE4];
			 if ((sl[i1].iN4>-1) && (fabs(sl[i1].an4) > nonzeroEPS)) buf -= sl[i1].an4*dX0[sl[i1].iN4];
			 if ((sl[i1].iS4>-1) && (fabs(sl[i1].as4) > nonzeroEPS)) buf -= sl[i1].as4*dX0[sl[i1].iS4];
			 if ((sl[i1].iT4>-1) && (fabs(sl[i1].at4) > nonzeroEPS)) buf -= sl[i1].at4*dX0[sl[i1].iT4];
			 if ((sl[i1].iW4>-1) && (fabs(sl[i1].aw4) > nonzeroEPS)) buf -= sl[i1].aw4*dX0[sl[i1].iW4];
	         buf*=buf;
		     res_sum+=buf;
	    }
	    for (integer i1=0; i1<maxbound; i1++) {
	    	// граничные узлы.
		    doublereal buf=0.0;
		    buf=slb[i1].aw*dX0[slb[i1].iW]-dV[slb[i1].iW];
		    if ((slb[i1].iI>-1) && (fabs(slb[i1].ai) > nonzeroEPS)) buf-=slb[i1].ai*dX0[slb[i1].iI];
		    buf*=buf;
		    res_sum+=buf;
	   }
	   res_sum=sqrt(res_sum);

	   if ((iVar == TEMP) && ((adiabatic_vs_heat_transfer_coeff == NEWTON_RICHMAN_BC)||(breakRUMBAcalc_for_nonlinear_boundary_condition==true))) {
		   // Если мы решаем Задачу с условием Ньютона Рихмана то у нас нету ни одного условия Дирихле.
		   // Сходимость определяется на глобальном уровне в solve nonlinear temp.

		   // Если у нас выставлен Стефан - Больцман то мы тоже не запускаем BiCGStab + ILU2.
	   }
	   else {

		   if (res_sum > res_sum_start) {
			   if ((iVorst_version == 2) && (bOkfgmres_amg1r5)) {
				   // Для решения был вызван гибкий вариант обобщенного 
				   // метода минимальных невязок (iVorst_version == 2) и 
				   // он точно сошелся (bOkfgmres_amg1r5) 
				   // (невязка опустилась менее 1E-7 и этого 
				   // в полной мере достаточно для достижения сходимости).
				   // Поэтому вызов вспомогательного решателя отменён, т.к. 
				   // всё впорядке, решение получено.
				   res_sum = res_sum_start;
			   }
			   else {
				   // расходимость нужен перезапуск.
				   printf("\namg solver divergence detected...\n");
				   //printf("\nno panic: because WACEB run. restart solver.\n");
				   printf("iVar=");
				   switch (iVar) {
				   case VX: printf("VX\n"); break;
				   case VY: printf("VY\n"); break;
				   case VZ: printf("VZ\n"); break;
				   case PAM: printf("PAM\n"); break;
				   case TEMP: printf("TEMP\n"); break;
				   }
				   printf("\nplease wait...\n");

				   // смена начального приближения и новый запуск.
				   iprogon++; // в случае расходимости мы будем повторно производить решение.

				   integer maxit = 2000;
				   bool bprintmessage = false;
				   Bi_CGStab_internal3(sl, slb, maxelm, maxbound, dV, dX0, maxit, alpharelax, bprintmessage, iVar, m, ifrontregulationgl, ibackregulationgl);

				   ///goto LabelAMGdivergenceDetected;
			   }
		   }
	   }
	   //printf("residual finish=%1.4e\n",res_sum);
	   //system("pause");
	   if (bsolid_static_only) {
		   // используется только для теплопередачи в твёрдом теле для ускорения
		   // решения задачи - защита от рестарта.
	       finish_residual=res_sum; // значение невязки решённой задачи.
	   }
	   
		 }

         calculation_main_end_time=clock();
         calculation_vorst_seach_time+=calculation_main_end_time-calculation_main_start_time;

} // amg_global_memory

// Здесь содержится обвязка вызывающая amg1r5.
// amg1r5 теперь полностью работоспособен на АЛИС сетке 29.09.2016.
void amg(equation3D* &sl, equation3D_bon* &slb,
			   integer maxelm, integer maxbound,
			   doublereal *dV, doublereal* &dX0, 
			   doublereal alpharelax, integer iVar,
	           bool bLRfree, QuickMemVorst& m,
	           integer* &ifrontregulationgl, integer* &ibackregulationgl,
	           integer iVorst_version) 
{

	              // iVorst_version == 0 - просто amg1r5 алгоритм.
	              // iVorst_version == 1 - amg1r5 алгоритм является предобуславливателем к алгоритму Хенка Ван дер Ворста BiCGStab.

				   bool bmemory_local=false;
				   if ((iVar == VX) || (iVar == VY) || (iVar == VZ) || (iVar == PAM)) {
					   bmemory_local = true;//4.08.2018
				   }
				   if (adiabatic_vs_heat_transfer_coeff == NEWTON_RICHMAN_BC) {
					   // Нелинейное условие Ньютона - Рихмана 
					   // матрица всё равно пересобирается каждый раз.
					 //  bmemory_local = true;
				   }

				   if (bmemory_local) {
					   // локальное выделение памяти , много alloc и free.
					   amg_loc_memory(sl, slb, maxelm,  maxbound, dV, dX0, alpharelax, iVar, bLRfree,m, iVorst_version);
				   }
				   else {
					   // память выделяется лишь единожды.
					   amg_global_memory(sl, slb, maxelm,  maxbound, dV, dX0, alpharelax, iVar, bLRfree,m, ifrontregulationgl, ibackregulationgl, iVorst_version);
				   }
}