/*
 *  ============================================================================= 
 *  ALADDIN Version 1.0 :
 *         elmt_fiber.c : Linear/Nonlinear Fiber Finite Element
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
 *  Written by: Wane-Jang Lin                                            May 1996
 *  ============================================================================= 
 */

#include "defs.h"
#include "units.h"
#include "matrix.h"
#include "fe_database.h"
#include "fe_functions.h"
#include "symbol.h"
#include "elmt.h"

/*
#define DEBUG 
*/


/* ============================================================== */
/*   Element FIBER                                                */
/*   2D   Frame Element: Beam_Col Elmt                            */
/*        Frame element :   material properties array             */
/*        Input Properties:                                       */
/* ============================================================== */

/*    p->work_material[0] = E;
      p->work_material[1] = G;
      p->work_material[2] = fy;
      p->work_material[3] = ET;
      p->work_material[4] = nu;
      p->work_material[5] = density;
      p->work_material[6] = fu;

      p->work_section[0] = Ixx;
      p->work_section[1] = Iyy;
      p->work_section[2] = Izz;
      p->work_section[3] = Ixy;
      p->work_section[4] = Ixz;
      p->work_section[5] = Iyz;
      p->work_section[6] = weight;
      p->work_section[7] = bf;
      p->work_section[8] = tf;
      p->work_section[9] = depth;                                  
      p->work_section[10] = area;
      p->work_section[11] = plate_thickness;
      p->work_section[12] = J;
      p->work_section[13] = rT;
      p->work_section[14] = width;
      p->work_section[15] = tw;                                   */
/* ============================================================== */

typedef struct load_history {
        int     elmt_no;
        MATRIX  *Dx, *dx, *rx;
        MATRIX  *Q, *q;
        MATRIX  *sr, *er, *s0, *e0, *sx, *ex;
        int     **yielding, **pre_load, **pre_range, **loading;
} HISTORY_DATA;

typedef struct fiber_respond {
        int   total_fiber_elmt;
        HISTORY_DATA  *history;
} FIBER_RESPOND;

static  FIBER_RESPOND  *FiberRespondBuffer;

#ifdef __STDC__
ARRAY *elmt_fiber(ARRAY *p, int isw)
#else
ARRAY *elmt_fiber(p, isw)
ARRAY *p;
int     isw;
#endif
{

    return(p);
}

extern  EFRAME  *frame;

/* 
 * ==================================================================   
 * Setup the static flags, stresses and strains to store load history
 * ==================================================================   
 */ 

#ifdef  __STDC__
void SetUpFiberRespondBuff( int total_fiber_elmt, EFRAME *frp )
#else
void SetUpFiberRespondBuff( total_fiber_elmt, frp )
int  total_fiber_elmt;
EFRAME *frp;
#endif
{
HISTORY_DATA   *hp;
ELEMENT_ATTR   *eap;
FIBER_ELMT     *fep;
int   ii, jj, elmt_no;

    FiberRespondBuffer = (FIBER_RESPOND *)MyCalloc(1, sizeof(FIBER_RESPOND));
    FiberRespondBuffer->total_fiber_elmt = total_fiber_elmt;
    FiberRespondBuffer->history = (HISTORY_DATA *)MyCalloc(total_fiber_elmt, sizeof(HISTORY_DATA));

    jj = 0;
    for( ii=1 ; ii <= frame->no_elements; ++ii )
    {
      eap = lookup(frame->element[ii-1].elmt_attr_name)->u.eap;
      if( !(strcmp(eap->elmt_type, "FIBER")) )
      {
         fep = lookup(eap->fiber_attr_name)->u.fep;
	 elmt_no = ii;
	 jj++;
         hp = &FiberRespondBuffer->history[jj-1];
	 hp->elmt_no = elmt_no;
      }
    }
}

#ifdef __STDC__
ARRAY *sld02(ARRAY *p, int isw)
#else
ARRAY *sld02(p, isw)
ARRAY *p;
int isw;
#endif
{
    printf("ERROR >> In sld02() : elmt no =%3d : isw= %3d\n",p->elmt_no, isw);
    return(p);
}

