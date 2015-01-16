/*
 *  ============================================================================= 
 *  ALADDIN Version 1.0 :
 *               elmt.h : Definitions for finite element library.
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

#ifndef ELMT_H
#define ELMT_H

/* element library functions */

extern ARRAY     *elmlib();
extern ARRAY     *elmt_fiber();
extern ARRAY     *elmt_frame_2d();
extern ARRAY     *elmt_frame_3d();
extern ARRAY     *elmt_psps();
extern ARRAY     *elmt_plate();
extern ARRAY     *elmt_shell_4nodes_q();   /* Shell Element with Drill DOF */
extern ARRAY     *elmt_shell();            /* X.G. shell element, 4 node   */
extern ARRAY     *elmt_shell_8n();         /* X.G. shell element, 8 node   */

/* Equivalent Loading Procedure for elements */

extern ARRAY     *sld01();   /* elmt_psps.c : PLANE_STRAIN and PLANE_STRESS */
extern ARRAY     *sld02();   /* elmt_fiber.c : FIBER */
extern ARRAY     *sld04();   /* elmt_shell_4n.c : SHELL_4N and SHELL_4NQ */
extern ARRAY     *sld05();   /* elmt_frame3d.c : FRAME_3D */
extern ARRAY     *sld07();   /* elmt_frame2d.c : FRAME_2D */
extern ARRAY     *sld08();   /* elmt_plate.c   : DKT_PLATE */
extern ARRAY     *sld108();  /* elmt_shell_8n.c : SHELL_8N */

/* Generic Template for Item in Finite Element Library */

static struct {
	char                 *name;         /* name of elment type                      */
        ARRAY  *(*elmt_lib_func)();         /* pointer to elmt library function         */
        ARRAY  *(*elmt_sld_func)();         /* pointer to elmt library sld0*() function */
        int                 no_dof;         /* No dof per node                          */
        int       no_node_per_elmt;         /* No nodes per element                     */
        int               no_dimen;         /* No dimension of problem                  */
	} elmt_library[] = {
		"FRAME_2D",          elmt_frame_2d,          sld07,  3, 2, 2,
		"FRAME_3D",          elmt_frame_3d,          sld05,  6, 2, 3,
                "SHELL_4N",          elmt_shell,             sld04,  5, 4, 3,
                "SHELL_4NQ",         elmt_shell_4nodes_q,    sld04,  6, 4, 3,
                "SHELL_8N",          elmt_shell_8n,          sld108, 5, 8, 3,
		"PLANE_STRAIN",      elmt_psps,              sld01,  2, 4, 2,
		"PLANE_STRESS",      elmt_psps,              sld01,  2, 4, 2,
        	"DKT_PLATE",         elmt_plate,             sld08,  3, 4, 3,
        	"FIBER",             elmt_fiber,             sld02,  3, 2, 2,
	};

#define NO_ELEMENTS_IN_LIBRARY (sizeof(elmt_library)/sizeof(elmt_library[0]))

/* ------------------------- */
/* Cases for Element Library */
/* ------------------------- */

#define PROPTY           1 
#define CHERROR          2   
#define ELMT_S_MAT       3
#define STIFF            3 
#define AXISYM           3
#define STRESS           4
                                        /* #define MASS           5  */
#define LOAD             6
#define PRESSLD          7
#define STRESS_LOAD      8
#define EQUIV_NODAL_LOAD 9

#define ATTR           15

/* ======================== */
/* element_type definitions */
/* ======================== */

#define NO_ELEMENT_TYPES 6

#define COLUMN         1
#define BEARING        2
#define GIRDER         3
#define SHEARGIRDER    4
#define DISSIPATOR     5
#define BRACE          6
#define GRID           9 
#define STOREY        10 

#endif /* end case ELMT_H */
