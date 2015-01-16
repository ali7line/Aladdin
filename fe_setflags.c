/*
 *  ============================================================================= 
 *  ALADDIN Version 1.0 :
 *        fe_setflags.c : Finite Element Preprocessor & Base Module
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

#include "units.h"
#include "matrix.h"
#include "fe_database.h"
#include "symbol.h"
#include "fe_functions.h"
#include "defs.h"

unsigned int PRINT_PROFILE;
unsigned int PRINT_PLINK;
unsigned int PRINT_EFRAME;
unsigned int PRINT_E_NUM;
unsigned int PRINT_E_STIFF;
unsigned int PRINT_G_STIFF;
unsigned int PRINT_E_MASS;
unsigned int PRINT_R_MASS;
unsigned int PRINT_G_MASS;
unsigned int PRINT_E_DESTIN;
unsigned int PRINT_R_DESTIN;
unsigned int PRINT_E_DESTIN_SIZE;
unsigned int PRINT_R_DESTIN_SIZE;
unsigned int PRINT_LOAD_VECTOR;
unsigned int PRINT_RESPONSE_VECTOR;
unsigned int PRINT_DISP;
unsigned int PRINT_NODAL_DISP;
unsigned int PRINT_RIGIDBODY_DISP;
unsigned int PRINT_STRESS;
unsigned int PRINT_STORY_DRIFT;
unsigned int PRINT_ELEM_LOAD;
unsigned int PRINT_FEF;
int mtype;


/* -------------------------------------------------------- */
/* set_default_values();                                    */
/* function to set default values for problem/analysis      */
/* called in fera_preprocessor()                            */
/*                                                          */
/* These parameters are set as default values for analysis  */
/* and can be overridden by  values given in input file     */
/* -------------------------------------------------------- */

set_default_values() 
{
      mtype = LUMPED; 
}

/* -------------------- */
/* set_print_output()   */
/* -------------------- */

set_print_output()
{

     /* set default values */

     PRINT_PROFILE          = YES;
     PRINT_PLINK            = YES;
     PRINT_EFRAME           = YES;
     PRINT_E_NUM            = 0;
     PRINT_E_STIFF          = NO;
     PRINT_G_STIFF          = YES;
     PRINT_E_MASS           = NO;
     PRINT_R_MASS           = YES;
     PRINT_G_MASS           = YES;
     PRINT_E_DESTIN         = NO;
     PRINT_R_DESTIN         = NO;
     PRINT_E_DESTIN_SIZE    = NO;
     PRINT_R_DESTIN_SIZE    = NO;
     PRINT_LOAD_VECTOR      = YES;
     PRINT_RESPONSE_VECTOR  = NO;

     PRINT_DISP             = OFF;
     PRINT_NODAL_DISP       = YES;
     PRINT_RIGIDBODY_DISP   = YES;
     PRINT_STRESS           = OFF;
     PRINT_STORY_DRIFT      = NO;
     PRINT_ELEM_LOAD        = NO;
     PRINT_FEF              = NO;

     return(1);
}
