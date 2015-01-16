/*
 *  ============================================================================= 
 *  ALADDIN Version 1.0 :
 *             vector.h : Data Structures for Vector Module                                  
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
 *  Written by: M.A. Austin                                             July 1993
 *  ============================================================================= 
 */

#ifndef VECTOR_H
#define VECTOR_H

/* Vector Data Structure */

typedef struct {
        char *cpVectorName;    /*  *name  */
        int        iLength;    /*  length  */
        DATA_TYPE    eType;    /*  type  */
        union {
            int    *ia;        /*  *i  */
            double *da;        /*  *d  */
        } uVector;             /*  array  */
} VECTOR;

#if (__STDC__ == 1)

/* Function Declarations for Standard ANSI C */

extern VECTOR *VectorAlloc( char * , DATA_TYPE , int );
extern VECTOR *VectorAdd( VECTOR * , VECTOR * );
extern VECTOR *VectorSub( VECTOR * , VECTOR * );
extern void    VectorPrint( VECTOR * );
extern void    VectorFree( VECTOR * );

extern void    PrintVectorInteger( VECTOR * );
extern void    VectorFreeInteger( VECTOR * );
extern VECTOR *VectorAddInteger( VECTOR * , VECTOR * );
extern VECTOR *VectorSubInteger( VECTOR * , VECTOR * );
extern int    *iVectorAlloc( int );

extern void    PrintVectorDouble( VECTOR * );
extern void    VectorFreeDouble( VECTOR * );
extern VECTOR *VectorAddDouble( VECTOR *, VECTOR *);
extern VECTOR *VectorSubDouble( VECTOR *, VECTOR *);
extern double *dVectorAlloc( int );

extern VECTOR *NaiveGaussElimination( MATRIX *, VECTOR *);
extern VECTOR *GaussElimination( char *, MATRIX *, VECTOR *);
extern VECTOR *SetupScaleFactors( MATRIX * );
extern VECTOR *SetupPivotVector( MATRIX * );

extern MATRIX *LUDecompositionIndirect( MATRIX *, VECTOR *);
extern MATRIX *LUSubstitutionIndirect( char *, VECTOR *, MATRIX *, MATRIX *);

#else  /* Start case not STDC */

/* Function Declarations for K&R C */

extern VECTOR *VectorAlloc();
extern VECTOR *VectorAdd();
extern VECTOR *VectorSub();
extern void    VectorPrint();
extern void    VectorFree();

extern void    PrintVectorInteger();
extern void    VectorFreeInteger();
extern VECTOR *VectorAddInteger();
extern VECTOR *VectorSubInteger();
extern int    *iVectorAlloc();

extern void    PrintVectorDouble();
extern void    VectorFreeDouble();
extern VECTOR *VectorAddDouble();
extern VECTOR *VectorSubDouble();
extern double *dVectorAlloc();

extern VECTOR *NaiveGaussElimination();
extern VECTOR *GaussElimination();
extern VECTOR *SetupScaleFactors();
extern VECTOR *SetupPivotVector();

extern MATRIX *LUDecompositionIndirect();
extern MATRIX  *LUSubstitutionIndirect();

#endif /* End case not STDC */

#endif /* end case VECTOR_H */
