/*
 *  ============================================================================= 
 *  ALADDIN Version 1.0 :
 *             matrix.h : Data Structures and Function Declarations for Matrix
 *                                                                     
 *  Copyright (C) 1995 by Mark Austin, Xiaoguang Chen, and Wane-Jang Lin
 *  Institute for Systems Research,                                           
 *  University of Maryland, College Park, MD 20742                                   
 *                                                                     
 *  This software is provided "as is" without express or implied warranty.
 *  Permission is granted to use this software for any on any computer system
 *  and to redistribute it freely, subject to the following restrictions:
 * 
 *  1. The authors are not responsible for the consequences of use of
 *     this software, even if they arise from defects in the software.
 *  2. The origin of this software must not be misrepresented, either
 *     by explicit claim or by omission.
 *  3. Altered versions must be plainly marked as such, and must not
 *     be misrepresented as being the original software.
 *  4. This notice is to remain intact.
 *                                                                    
 *  Written by: Mark Austin and Wane-Jang Lin                       December 1995
 *  ============================================================================= 
 */

#ifndef MATRIX_H
#define MATRIX_H

/* [a] : Data Structures for Matrices of numbers/engineering quantities */

typedef struct {
        double dReal, dImaginary;
} COMPLEX;

typedef enum {
	INTEGER_ARRAY  = 1,
	DOUBLE_ARRAY   = 2,
	COMPLEX_ARRAY  = 3
} DATA_TYPE;

typedef enum {
	SEQUENTIAL = 1,
	INDIRECT   = 2,
	SKYLINE    = 3,
	SPARSE     = 4
} INTERNAL_REP;

typedef struct matrix {
        char      *cpMatrixName;    /*  *name   */
        int             iNoRows;    /*  no_rows   */
        int          iNoColumns;    /*  no_columns   */
	DIMENSIONS  *spRowUnits;    /*  *row_units_buf  */
	DIMENSIONS  *spColUnits;    /*  *col_units_buf  */
        INTERNAL_REP       eRep;
        DATA_TYPE         eType;    /*  type  */
        union {
            int           **iaa;
            double        **daa;    /*  **d  */
            COMPLEX       **caa;
        } uMatrix;
} MATRIX;

/* [b] : Declarations for Root Matrix Functions */

#ifdef __STDC__ 

extern MATRIX *MatrixAllocate( MATRIX * );
extern MATRIX *MatrixDiag( MATRIX * );
extern MATRIX *MatrixZero( MATRIX * );
extern MATRIX *MatrixOne( MATRIX * );
extern MATRIX *MatrixAdd( MATRIX * , MATRIX * );
extern MATRIX *MatrixAddReplace( MATRIX *, MATRIX * );
extern MATRIX *MatrixSub( MATRIX * , MATRIX * );
extern MATRIX *MatrixSubReplace( MATRIX *, MATRIX * );
extern MATRIX *MatrixMult( MATRIX * , MATRIX * );
extern MATRIX *MatrixPower( MATRIX * , QUANTITY * );
extern MATRIX *MatrixNegate( MATRIX * );
extern MATRIX *MatrixNegateReplace( MATRIX *);
extern MATRIX *MatrixTranspose( MATRIX * );
extern MATRIX *MatrixCopy( MATRIX * );
extern MATRIX *MatrixSolve( MATRIX *, MATRIX * );
extern MATRIX *MatrixLU( MATRIX * );
extern MATRIX *MatrixFB( MATRIX *, MATRIX * );
extern MATRIX *MatrixInverse( MATRIX * );
extern MATRIX *MatrixDimension( MATRIX * );
extern MATRIX *MatrixScale( MATRIX *, double );
extern double  MatrixContentScale( MATRIX *, int, int );
extern MATRIX *MatrixQuanMult( QUANTITY *, MATRIX * );
extern MATRIX *MatrixQuanDiv( MATRIX * , QUANTITY * );
extern MATRIX *MatrixZeroUnits( MATRIX *, int, int );
extern MATRIX *MatrixUnitsLess( MATRIX * );
extern MATRIX *MatrixUnitsSimplify( MATRIX * );

extern MATRIX *MatrixPrintVar( MATRIX *, ... );
extern MATRIX *MatrixPrint( MATRIX *, ... );
extern MATRIX *MatrixColumnUnits( MATRIX *, ... );
extern MATRIX *MatrixRowUnits( MATRIX *, ... );
extern MATRIX *MatrixExtract( MATRIX *, ... );
extern MATRIX *MatrixPut( MATRIX *, ... );

extern void      MatrixFree( MATRIX * );
extern QUANTITY *MatrixDet( MATRIX * );
extern QUANTITY *MatrixMax( MATRIX * );
extern QUANTITY *MatrixMin( MATRIX * );
extern QUANTITY *MatrixL2Norm( MATRIX * );
extern QUANTITY *QuantityCast( MATRIX * );

/* [c] : Matrix Functions with INDIRECT Storage Pattern */

extern MATRIX   *MatrixAllocIndirect( char *, DATA_TYPE , int , int );
extern double  **MatrixAllocIndirectDouble( int, int );
extern int     **MatrixAllocIndirectInteger( int, int );
extern void      MatrixFreeIndirectDouble( double **, int );
extern void      MatrixFreeIndirectInteger( int **, int );

extern void    MatrixPrintIndirectDouble( MATRIX * );
extern void    MatrixPrintIndirectInteger( MATRIX * );

extern MATRIX *MatrixCopyIndirectDouble( MATRIX * );
extern MATRIX *MatrixAddIndirectDouble( MATRIX *, MATRIX *);
extern MATRIX *MatrixAddReplaceIndirectDouble( MATRIX *, MATRIX *);
extern MATRIX *MatrixSubIndirectDouble( MATRIX *, MATRIX *);
extern MATRIX *MatrixSubReplaceIndirectDouble( MATRIX *, MATRIX *);
extern MATRIX *MatrixMultIndirectDouble( MATRIX *, MATRIX *);
extern MATRIX *MatrixMultIndirectSkylineDouble( MATRIX *, MATRIX *);
extern MATRIX *MatrixMultSkylineIndirectDouble( MATRIX *, MATRIX *);
extern MATRIX *MatrixNegateIndirectDouble( MATRIX * );
extern MATRIX *MatrixNegateReplaceIndirectDouble( MATRIX * );

extern MATRIX *MatrixTransposeIndirectDouble( MATRIX * );
extern MATRIX *MatrixInverseIndirectDouble( MATRIX * );
extern MATRIX *MatrixScaleIndirectDouble( MATRIX *, double );
extern double  MatrixContentScaleIndirectDouble( MATRIX *, int, int );

/* [d] : Matrix Functions with SKYLINE Storage Pattern */

extern MATRIX *MatrixAllocSkyline( char *, DATA_TYPE, int, int, int *);
extern void    MatrixFreeSkyline( MATRIX * );
extern void    MatrixPrintSkylineDouble( MATRIX * );

extern MATRIX *MatrixReallocSkyline( MATRIX * );
extern MATRIX *MatrixReallocSkylineDouble( MATRIX * );

extern MATRIX *MatrixAddSkyline( MATRIX *, MATRIX * );
extern MATRIX *MatrixSubSkyline( MATRIX *, MATRIX * );
extern MATRIX *MatrixMultSkyline( MATRIX *, MATRIX *);
extern MATRIX *MatrixNegateSkyline( MATRIX *);
extern MATRIX *MatrixNegateReplaceSkyline( MATRIX * );
extern MATRIX *MatrixCopySkyline( MATRIX * );
extern MATRIX *MatrixTransposeSkyline( MATRIX * );
extern MATRIX *MatrixInverseSkyline( MATRIX * );
extern MATRIX *MatrixScaleSkyline( MATRIX *, double );
extern double  MatrixContentScaleSkyline( MATRIX *, int, int );

extern MATRIX *LUDecompositionSkyline( MATRIX *);
extern MATRIX *LUBacksubstitutionSkyline( MATRIX *, MATRIX *);

extern MATRIX *MatrixIndirectToSkyline( MATRIX * );
extern MATRIX *MatrixSkylineToIndirect( MATRIX * );

extern MATRIX *MatrixAssembleSkyline( MATRIX *, MATRIX *, int *, int * );
extern MATRIX *CholeskyDecompositionIndirect( MATRIX * );

extern void    MatrixSolveEigen( MATRIX *, MATRIX *, MATRIX *, MATRIX *, int );
extern MATRIX *Solve_Eigen( MATRIX *, MATRIX *, MATRIX * );
extern MATRIX *Extract_Eigenvalue( MATRIX * );
extern MATRIX *Extract_Eigenvector( MATRIX * );
extern void    Print_Eigen( MATRIX * );

extern void         dMatrixPrint( char *, double **, int, int );
extern double     **dMatrixCopy( double **, int, int );
extern double     **dMatrixCopyRep( double **, double **, int, int );
extern double     **dVmatrixCrossProduct( double **, double **, int, int, double **, int, int );
extern double       dVmatrixInnerProduct( double **, int, int, double **, int, int );
extern double     **dMatrixMult( double **, int, int, double **, int, int );
extern double     **dMatrixMultRep( double **, double **, int, int, double **, int, int );
extern double     **dMatrixTranspose( double **, int, int );
extern double       dVmatrixL2Norm( double **, int, int );
extern double       dMatrixDet( double **, int, int );

#else  /* start case not STDC */

/* [b] : Declarations for Root Matrix Functions */

extern MATRIX *MatrixPrint();
extern MATRIX *MatrixPrintVar();
extern MATRIX *MatrixAllocate();
extern MATRIX *MatrixDiag();
extern MATRIX *MatrixZero();
extern MATRIX *MatrixOne();
extern MATRIX *MatrixScale();
extern double  MatrixContentScale();
extern MATRIX *MatrixCopy();
extern MATRIX *MatrixTranspose();
extern MATRIX *MatrixDimension();
extern MATRIX *MatrixAdd();
extern MATRIX *MatrixAddReplace();
extern MATRIX *MatrixSub();
extern MATRIX *MatrixSubReplace();
extern MATRIX *MatrixMult();
extern MATRIX *MatrixPower();
extern MATRIX *MatrixNegate();
extern MATRIX *MatrixNegateReplace();
extern void    MatrixFree();
extern MATRIX *MatrixSolve();
extern MATRIX *MatrixLU();
extern MATRIX *MatrixFB();
extern void    MatrixSolveEigen();
extern MATRIX *MatrixInverse();
extern QUANTITY  *MatrixDet();
extern QUANTITY  *MatrixL2Norm();
extern QUANTITY  *MatrixMax();
extern QUANTITY  *MatrixMin();

/* [b.1] : Operations between MATRIX and QUANTITY */

extern MATRIX *MatrixQuanMult();
extern MATRIX *MatrixQuanDiv();

/* [b.2] : Declarations for Matrix Functions about Units */

extern MATRIX      *MatrixColumnUnits();
extern MATRIX      *MatrixRowUnits();
extern MATRIX      *MatrixZeroUnits();
extern MATRIX      *MatrixUnitsSimplify();
extern MATRIX      *MatrixUnitsLess();
extern MATRIX      *MatrixExtract();
extern MATRIX      *MatrixPut();
extern QUANTITY    *QuantityCast();


/* [c] : Matrix Functions with INDIRECT Storage Pattern */

extern MATRIX  *MatrixAllocIndirect();
extern double **MatrixAllocIndirectDouble();
extern int    **MatrixAllocIndirectInteger();

extern void     MatrixFreeIndirectDouble();
extern void     MatrixFreeIndirectInteger();

extern MATRIX *MatrixCopyIndirectDouble();
extern MATRIX *MatrixScaleIndirectDouble();
extern double  MatrixContentScaleIndirectDouble();

extern MATRIX *MatrixAddIndirectDouble();
extern MATRIX *MatrixAddReplaceIndirectDouble();
extern MATRIX *MatrixSubIndirectDouble();
extern MATRIX *MatrixSubReplaceIndirectDouble();
extern MATRIX *MatrixNegateIndirectDouble();
extern MATRIX *MatrixNegateReplaceIndirectDouble();

extern MATRIX *MatrixMultIndirectDouble();
extern MATRIX *MatrixTransposeIndirectDouble();

extern MATRIX *LUDecompositionIndirect();
extern MATRIX *LUSubstitutionIndirect();

extern MATRIX *MatrixInverseIndirectDouble();

/* [d] : Matrix Functions with SKYLINE Storage Pattern */

extern MATRIX *MatrixAllocSkyline();
extern void    MatrixPrintSkylineDouble();
extern void    MatrixFreeSkyline();

extern MATRIX *ReSkyMatrix();
extern MATRIX *ReSkyMatrixDouble();

extern MATRIX *MatrixAddSkyline();
extern MATRIX *MatrixSubSkyline();
extern MATRIX *MatrixMultSkyline();
extern MATRIX *MatrixNegateSkyline();
extern MATRIX *MatrixNegateReplaceSkyline();

extern MATRIX *MatrixMultIndirectSkylineDouble();
extern MATRIX *MatrixMultSkylineIndirectDouble();

extern MATRIX *MatrixCopySkyline();
extern MATRIX *MatrixScaleSkyline();
extern double  MatrixContentScaleSkyline();

extern MATRIX *LUDecompositionSkyline();
extern MATRIX *LUBacksubstitutionSkyline();

extern MATRIX *MatrixTransposeSkyline();
extern MATRIX *MatrixInverseSkyline();

extern MATRIX *MatrixIndirectToSkyline();
extern MATRIX *MatrixSkylineToIndirect();

extern MATRIX *MatrixAssembleSkyline();
extern MATRIX *CholeskyDecompositionIndirect();

extern MATRIX *Solve_Eigen();
extern MATRIX *Extract_Eigenvalue();
extern MATRIX *Extract_Eigenvector();
extern void    Print_Eigen();

/* [e] : Declarations for Double Matrix (without units) Functions */

extern void         dMatrixPrint();
extern double     **dMatrixCopy();
extern double     **dMatrixCopyRep();
extern double     **dVmatrixCrossProduct();
extern double       dVmatrixInnerProduct();
extern double     **dMatrixMult();
extern double     **dMatrixMultRep();
extern double     **dMatrixTranspose();
extern double       dVmatrixL2Norm();
extern double       dMatrixDet();

#endif /* end case not STDC */

#endif /* end case MATRIX_H */
