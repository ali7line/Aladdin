/*
 *  ============================================================================= 
 *  ALADDIN Version 1.0 :
 *       fe_functions.h : External Function Declarations for Finite Elements
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

#ifndef FE_FUNCTIONS_H
#define FE_FUNCTIONS_H

#include <stdio.h>
#include <malloc.h>
#include <math.h>
#include <string.h>

extern EFRAME    *FrameAlloc();
extern EFRAME    *profile();
extern EFRAME    *Set_Elmt_Attrs();
extern EFRAME    *assign_properties();
extern EFRAME    *plink();
extern EFRAME    *rplink();
extern void       print_property();

extern ARRAY     *Alloc_p_Array();
extern ARRAY     *Assign_p_Array();
extern ARRAY     *Element_Property();
extern ARRAY     *Mate_Property();
extern ARRAY     *Eload_Property();

extern MATRIX    *Element_Matrix();
extern MATRIX    *Element_Equiv();
extern QUANTITY  *Element_Vector();

extern MATRIX    *Assemble_Global();
extern MATRIX    *Assemble_Global_Load();

extern int       *Destination_Array();
extern int        Bound_Disp();
extern double   **Boundary_Conditions();
extern double   **MATER_MAT_PLANE();    /* function for elmt_psps.c only */

extern QUANTITY  *pload();
extern QUANTITY  *Modify_Load();
extern QUANTITY  *Addload_Vector();
extern QUANTITY  *Transform_Force();
extern QUANTITY  *Nodal_Forces();

extern QUANTITY  *Assemble_Nodal_Load();
extern QUANTITY  *Assemble_Gravity_Load();
extern QUANTITY  *Assemble_Ctrfgl_Load();
extern double    *Assemble_Equiv_Load(); 

extern int       *Destination_Array_for_Rigid_Body();
extern ARRAY     *Assign_p_Array_for_Rigid_Body();
extern int        Rigidbody_conn();
extern MATRIX    *Transform_Stiff_Matrix();
extern MATRIX    *Transform_Rigid_Body_Mass_Matrix();
extern double   **Transformation_Matrix();
extern double   **Modify_T_Matrix();

/* functions called in code.c for finite element solution procedures */

extern  void      Start_Mesh();
extern  void      End_Mesh();
extern  void      Print_Mesh();
extern  void      Add_Node();
extern  void      Fix_Node();
extern  void      Node_Load();
extern  void      Link_Node();
extern  void      Add_Elmt();

extern  void      Print_Displ();
#ifdef __STDC__
extern MATRIX    *Print_Stress( MATRIX *, ... );
#else
extern MATRIX    *Print_Stress();
#endif

extern MATRIX    *Form_Stiffness();
extern MATRIX    *Form_Mass();
extern MATRIX    *Form_External_Load();
extern MATRIX    *Form_Equiv_Nodal_Load();
#ifdef __STDC__
extern MATRIX    *Form_Internal_Load( MATRIX *, ... );
#else
extern MATRIX    *Form_Internal_Load();
#endif

extern MATRIX    *Solve_Eigen();
extern MATRIX    *Velocity_Extract();
extern MATRIX    *Displacement_Extract();
extern void       Ldof_to_gdof();

/* Finite Element Allocation Routines */

extern ELEMENT_ATTR     *Alloc_Element_Attr_Item();
extern SECTION_ATTR     *Alloc_Section_Attr_Item();
extern MATERIAL_ATTR    *Alloc_Material_Attr_Item();
extern FIBER_ELMT       *Alloc_Fiber_Elmt_Attr_Item();

/* functions declarations for Add_Elmt, */

extern EFRAME    *CheckMaterialsSpace(); 
extern EFRAME    *CheckElementSpace();
extern EFRAME    *CheckRigidSpace();
extern EFRAME    *CheckJdiagSpace();   
extern EFRAME    *CheckNodeSpace();   
extern EFRAME    *CheckNforcesSpace(); 
extern EFRAME    *CheckEforcesSpace(); 
extern EFRAME    *CheckEtypesSpace(); 

/* functions for rule checking */

extern MATRIX    *Get_Coord();
extern MATRIX    *Get_Node();
extern MATRIX    *Get_Displ();
extern MATRIX    *Get_Stress();
extern MATRIX    *Get_Dof();
extern MATRIX    *Get_Section();
extern MATRIX    *Get_Material();

/* functions used in elmt_*.c */

extern MATRIX    *beamst();
extern MATRIX    *beamms();
extern MATRIX    *beamms3d();
extern MATRIX    *beamst3d();
extern double   **tmat();
extern double   **rotate();
extern double   **rotate3d();
extern int        pstres();
extern            gauss();
extern int        pgauss();
extern int        shape();
extern int        shp0();
extern double   **qushp8();
extern double   **dktqbm();
extern int        jacqud();
extern void       dktb06();
extern void       dktq06();
extern void       hshp06();
extern void       jacq06();
extern void       jtri06();
extern void       proj06();
extern void       rots06();
extern void       rshp06();
extern void       stre06();
extern void       tran06();
extern double    *pstres06();
extern void       shp_prt();

/* functions about 4 node shell elmt */

extern ARRAY     *elmt_shell_4nodes_implicit();
extern void       Lamina_Sys();
extern void       Lamina_Sys_Implicit();
extern void       elmt_shell_shape_4node();
extern double   **Hourglass_Stress_Rate();
extern void       Shell_4Node_Mass();
extern double   **B_MATRIX_4Node();
extern void       Shell_Stiff_Plane_4node();
extern double   **Shell_Nodal_Load_Plane();

/* functions about 8 node shell elmt */

extern ARRAY     *elmt_shell_8nodes_implicit();
extern void       Lamina_Sys_8node();
extern void       elmt_shell_shape_8node();
extern void       Shell_8Node_Mass();
extern void       Stress_Update_8Node();
extern double   **B_MATRIX_8Node();
extern void       Shell_Stiff_Plane_8node();
extern double   **Shell_Nodal_Load_8Node();

/* functions declarations for 4 Node and 8 Node shell elements */

extern void       MATER_SHELL_UPDATE();
extern double   **STRAIN_RATE_SHELL();
extern double   **Hourglass_Stress_Rate();
extern double   **Hourglass_stiff();
extern double   **Strain_Displ_Matrix();
extern double   **MATER_MAT_SHELL();
extern void       Load_Curve();
extern void       Plastic_Deform();

extern void       DISPL_UPDATE();
extern void       Stress_Update();
extern void       BB_Vector();

/* functions for non-linear analysis */

extern void       UpdateResponse();
extern void       save_action();
extern void       SetUpRespondBuffer();
extern void       SaveRespondBuffer();
extern void       SetUpFiberRespondBuff();

#endif /* end case FE_FUNCTIONS_H */
