/*
 *  ============================================================================= 
 *  ALADDIN Version 1.0 :
 *       elmt_dummy.c : Skeleton code for "dummy" finite element.
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

#include <math.h>
#include "defs.h"
#include "units.h"
#include "matrix.h"
#include "fe_database.h"
#include "symbol.h"
#include "fe_functions.h"
#include "elmt.h"

/* #define DEBUG */


/*  =======================================================================
 *  Main function for finite element matrix formulation.
 *
 *  The variable "isw" controls the purpose of the function call,
 *  with the various cases being as follows: 
 *
 *      case PROPTY:        Compute finite element material and section
 *                          properties.
 *      case CHERROR:       Check for errors (often this section is blank).
 *      case STRESS_LOAD:
 *      case STRESS_UPDATE:
 *      case PRESSLD:
 *      case LOAD_MATRIX:   Compute equivalent nodal load forces.
 *      case STRESS:        Compute and print element stresses.
 *      case STIFF:         Form finite element stiffness matrix.
 *      case MASS_MATRIX:   Form finite element mass matrix.
 *
 *  =======================================================================
 *
 *  Material Properties:
 *
 *        p->work_material[0] = E;
 *        p->work_material[1] = G;
 *        p->work_material[2] = fy;
 *        p->work_material[3] = ET;
 *        p->work_material[4] = nu;
 *        p->work_material[5] = density;
 *        p->work_material[6] = fu;
 *
 *  Section Properties:
 *
 *        p->work_section[0] = Ixx;
 *        p->work_section[1] = Iyy;
 *        p->work_section[2] = Izz;
 *        p->work_section[3] = Ixy;
 *        p->work_section[4] = Ixz;
 *        p->work_section[5] = Iyz;
 *        p->work_section[6] = weight;
 *        p->work_section[7] = bf;
 *        p->work_section[8] = tf;
 *        p->work_section[9] = depth;                                  
 *        p->work_section[10] = area;
 *        p->work_section[11] = plate_thickness;
 *        p->work_section[12] = J;
 *        p->work_section[13] = rT;
 *        p->work_section[14] = width;
 *        p->work_section[15] = tw;                                 
 *
 *  =======================================================================
 */

#ifdef __STDC__
ARRAY *elmt_frame_2d(ARRAY *p, int isw)
#else
ARRAY *elmt_frame_2d(p, isw)
ARRAY *p;
int     isw;
#endif
{

    UNITS_SWITCH = CheckUnits();

    switch(isw) {
        case PROPTY:
             printf("Assign material properties\n");
             printf("Assign section  properties\n");
             break;
        case CHERROR:
             break;
        case STRESS_LOAD:
             break;
        case STRESS_UPDATE:
             break;
	case PRESSLD:
             break;
	case LOAD_MATRIX:
        case STRESS:
             printf("Compute and print element level nodal forces and/or stresses\n");
             break;
        case STIFF: 
             printf("Form finite element stiffness matrix\n");
             break;
        case MASS_MATRIX:
             printf("Form finite element mass matrix\n");
             break;
        default:
             break;
    }

    return(p);
}


/* 
 *  ====================================================================
 *  Form Finite Element Stiffness Matrix
 * 
 *  Typical Input/Output Parameters : 
 * 
 *  Input  : ARRAY  *p     = p array for transferring properties.
 *         : MATRIX *s     = stiffness matrix.
 *         : double EA     = Young's modulus times Area.
 *         : double EI     = Flexural stiffness
 *         : double length = length of element.
 *         : double cs     = cos() of element orientation.
 *         : double sn     = sin() of element orientation.
 *         : double size_of_stiff = size of stiffness matrix.
 *         : double n_dof         = no degrees of freedom.
 *  Output : MATRIX *s     = stiffness matrix.
 *  ====================================================================
 */

#ifdef __STDC__
MATRIX *beamst(ARRAY *p, MATRIX *s, double EA, double EI, double length,
               double cs, double sn, int size_of_stiff, int no_dof)
#else
MATRIX *beamst(p, s, EA, EI, length, cs, sn, size_of_stiff, no_dof)
ARRAY  *p;
MATRIX *s;
double  EA, EI, length, cs, sn;
int     size_of_stiff, no_dof ;
#endif
{

    /* [a] : Use material/section properties to form elements of stiffness matrix */

    /* [b] : Add to stiffness matrix computation */

    /* [c] : Rotate local coordinates to global coordinate frame */

    return(s);
}

/* 
 *  ====================================================================
 *  Form Finite Element Mass Matrix
 * 
 *  Typical Input/Output Parameters : 
 * 
 *  Input  : ARRAY  *p     = p array for transferring properties.
 *         : MATRIX *s     = stiffness matrix.
 *         : double mtype  = mass matrix type (CONSISTENT/LUMPED).
 *         : double mass   = mass per unit length.
 *         : double xl     = element length.
 *         : double cs     = cos() of element orientation.
 *         : double sn     = sin() of element orientation.
 *         : double nst    = size of stiffness matrix.
 *         : double ndf    = no degrees of freedom.
 *  Output : MATRIX *s     = mass matrix.
 *  ====================================================================
 */

#ifdef __STDC__
MATRIX *beamms(ARRAY *p, MATRIX *s, int mtype, double mass, double xl,
              double cs, double sn,   int nst, int ndf )
#else
MATRIX *beamms(p, s, mtype, mass, xl, cs, sn, nst, ndf)
ARRAY  *p;
MATRIX *s;
double mass, xl, cs, sn;  /* mass, length, cosine, ans sine */
int nst, ndf,   mtype;
#endif
{

    /* [a] : Form lumped/consistent mass matrices */

    switch (mtype) {
	case LUMPED:
	     break;
	case CONSISTENT:
	     break;
	default:
             FatalError("In elmt_frame2d() : beamms() : Type of Mass Matrix Undefined",
                       (char*) NULL);
	     break;
    }

    /* [b] : Assign Units to mass matrix */

    /* [c] : Rotate mass matrix */

    return(s);
}


/* 
 *  ===========================================================
 *  Rotate element from local coordinates to global coordinates               
 *  
 *  Input  : double **s = pointer to mass/stiffness matrix.
 *         : double cs  = cos() of element orientation.
 *         : double sn  = sin() of element orientation.
 *         : double nst = size of stiffness matrix.
 *         : double ndf = no degrees of freedom.
 *  Output : double **s = pointer to mass/stiffness matrix.
 *  ===========================================================
 */ 

#ifdef __STDC__
double **rotate(double **s, double cs, double sn, int nst, int ndf)
#else
double **rotate(s, cs, sn, nst, ndf)
double **s;
int  nst, ndf;
double cs, sn;
#endif
{

}


/* 
 *  ===================================================
 *  Equivalent Loading Procedure for 2 D frame element 
 * 
 *  Input : 
 *  Output: 
 *  ===================================================
 */ 

#ifdef __STDC__
ARRAY *sld07(ARRAY *p, int task)
#else
ARRAY *sld07(p,task)
ARRAY *p;
int task;
#endif
{

    /* [a] : Initialize total load */

    /* [b] : Form load matrix for different loading types */

    switch(task){
       case PRESSLD:
       case STRESS:
            break;
       case STRESS_LOAD:
            break;
       default:
            break;
    }

    /* [c] : Rotate local forces to global forces */
 
    return(p);
}
