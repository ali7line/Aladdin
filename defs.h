/*
 *  ============================================================================= 
 *  ALADDIN Version 1.0 :
 *               defs.h : Standard definitions
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
 *  Written by: Mark Austin, Xiaoguang Chen, and Wane-Jang Lin      December 1995
 *  ============================================================================= 
 */

/*
 *  -----------------
 *  Numerical Values 
 *  -----------------
 */

#define MAXTOKSIZE             80
#define SHORTOKSIZE            16
#define MAX_ITERATIONS         75
#define MAXREAL            1.7e38
#define MINREAL           1.0e-37

#define PI     3.14159265358979323846
#define DEG   57.29577951308232087680

static double ETA = 2.938736e-38;
static double EPS = 1.387779e-16;
static double TOL = 2.117583e-22;

/*
 *  -----------------
 *  Math Macros 
 *  -----------------
 */

#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define SIGN(a)  ( a < 0 ? 1 : -1)
#define ABS(a)   (((a) >= 0) ? (a) : -(a))

#define isletter(c)   (isalpha(c) || (c) == '_')
#define cot(a) ((double)(cos((double)(a)) /sin((double)(a))))
#define square(x) x*x
#define Deg2Rad    (PI / 180.0)

/*
 *  -------------------
 *  Logical Definitions
 *  -------------------
 */

#define ONE          1

#define TRUE         1 
#define FALSE        0 
#define ON           1 
#define OFF          0 
#define NOT_FINISHED 1 
#define FINISHED     0 

#define PASSED       1    
#define FAILED       0   
#define YES          1    
#define NO           0   
#define UP           1
#define DOWN         0
#define COUNTED      1
#define NOTCOUNTED   2

#define FAILURE      0
#define SUCCESS      1

/* --------------------------------------------------------------- */
/* Finite Element : Global Arrays & numbers : THIS WILL BE CHANGED */
/* --------------------------------------------------------------- */

double  **R, **RT;
double  **K,  **M, **C, **L,   *d,    *F, *pe, **KBAC;
double  *Fe,  *B, *Qs, **pr, *nopl;

int   NSTEPS;
float Eq_angle;
float DT;
float BETA;

/* ------------ */
/* Problem Type */
/* ------------ */

#define LUMPED          1
#define CONSISTENT     -1
#define EQNOS           1
#define COLHGT          2
#define STATIC          1
#define DYNAMIC         2
#define IMPLICIT       -1
#define EXPLICIT        1

/* ------------------------- */
/* Cases for Element Library */
/* ------------------------- */

#define PROPTY             1 
#define CHERROR            2   
#define ELMT_S_MAT         3
#define STIFF              3 
#define AXISYM             3
#define STRESS             4
#define MASS_MATRIX        5 
#define LOAD_MATRIX        6
#define PRESSLD            7
#define STRESS_LOAD        8
#define EQUIV_NODAL_LOAD   9
#define ACCELERATION_LOAD  10
#define STRESS_UPDATE      11

#define ATTR               15
#define ELEMENT_LIBRARY    1000

/* --------------------------------------------- */
/* UNIT_QTYs  definitions for INITIAL ALLOCATION */
/* --------------------------------------------- */

#define UNIT_NDM             2
#define UNIT_NDF             3           
#define UNIT_NEN             8          /* max nodes per element */
#define UNIT_INTEG_PTS       2          /* 2 pts for linear, 5 for non-linear analysis   */
                                        /* this needs to be changed to adapt dynamic pts */
                                        /* for different problems and elements           */
#define UNIT_IN_PLANE_INTEG_PTS  4          

#define UNIT_NAD             0
#define UNIT_NODES           10
#define UNIT_ELEMENTS        10
#define UNIT_RIGIDS          10
#define UNIT_NFORCES         10 
#define UNIT_EFORCES         10

#define UNIT_SECTION_ATTR    20   /* must >= no. of variables in data structures SECTION_ATTR defined in fe_database.h */
#define UNIT_MATERIAL_ATTR   20   /* must >= no. of variables in data structures MATERIAL_ATTR defined in fe_database.h */

#define UNIT_ELEMENT_ATTR    17 
#define MAX_NO_MEMBER_LOADS  10

/* -------------------------------------------- */
/* UNIT_QTYs  definitions for Arrays of P_ARRAY */ 
/* -------------------------------------------- */

#define UNIT_ENST            6 
#define UNIT_ENDM            3 
#define UNIT_ENDF            6 
#define UNIT_ENEN            8   /* max nodes per element */

extern  int   max_no_nodes;
extern  int   max_no_eqs;
extern  int   max_no_elements;
extern  int   max_no_rigid;
extern  int   max_no_nforces;
extern  int   max_no_loaded_elmts;
extern  int   TDOF;
extern  int   TNEQ;
extern  int   MDOF;
extern  int   UNITS_TYPE;

int   H_Print;
int   STATE_PLASTIC;
int   NO_OF_ELMTS;
int   LOCAL_DOF;      /* total local dofs per element  */
int   NODES_PER_ELMT; /* not the max nodes per element */

/* -------------------------------------- */
/* Control integers for print statements  */
/* -------------------------------------- */

extern  unsigned int PRINT_PROFILE;
extern  unsigned int PRINT_PLINK;
extern  unsigned int PRINT_EFRAME;
extern  unsigned int PRINT_E_NUM;
extern  unsigned int PRINT_E_STIFF;
extern  unsigned int PRINT_G_STIFF;
extern  unsigned int PRINT_E_MASS;
extern  unsigned int PRINT_R_MASS;
extern  unsigned int PRINT_G_MASS;
extern  unsigned int PRINT_E_DESTIN;
extern  unsigned int PRINT_R_DESTIN;
extern  unsigned int PRINT_E_DESTIN_SIZE;
extern  unsigned int PRINT_R_DESTIN_SIZE;
extern  unsigned int PRINT_LOAD_VECTOR;
extern  unsigned int PRINT_RESPONSE_VECTOR;
extern  unsigned int PRINT_DISP;
extern  unsigned int PRINT_NODAL_DISP;
extern  unsigned int PRINT_RIGIDBODY_DISP;
extern  unsigned int PRINT_STRESS;
extern  unsigned int PRINT_STORY_DRIFT;
extern  unsigned int PRINT_ELEM_LOAD;
extern  unsigned int PRINT_FEF;

/* =================================== */
/* constraints and boundary conditions */
/* =================================== */

#define         FIXED  0
#define      NOTFIXED  1
#define UNCONSTRAINED  0

