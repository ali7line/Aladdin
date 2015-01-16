/*
 *  ============================================================================= 
 *  ALADDIN Version 1.0 :
 *             symbol.h : Declarations for Symbol Table
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
 *  Written by: Mark Austin                                         December 1995
 *  ============================================================================= 
 */

#ifndef SYMBOL_H
#define SYMBOL_H

#define NO_LINES_IN_HEADER           10 
#define MAX_NO_CHAR_PER_LINE        100

/*
 * ===============================
 * Data Structure for Symbol Table
 * ===============================
 */
 
typedef struct Symbtab_element {
	char    *cpSymName;
	short	type;
	union	{
		double                  value;   /* Store a number             */
		char                     *str;   /* String                     */
		int                 (*defn)();   /* Function, Procedure        */
		void             (*voidptr)();   /* Builtin Math/FE Functions  */
		double         (*doubleptr)();   /* Builtin Math/FE Functions  */
		MATRIX        *(*matrixptr)();   /* Builtin Matrix Functions   */
		QUANTITY    *(*quantityptr)();   /* Builtin Quantity Functions */
                ARRAY          *(*elmt_ptr)();   /* Elmt Library ptr to func   */
		QUANTITY                   *q;   /* Engineering Quantity       */
		MATRIX                     *m;   /* Matrix                     */
                DIMENSIONS             *dimen;   /* Basic Dimension            */
                ELEMENT_ATTR             *eap;   /* Element Attribute ptr      */
                SECTION_ATTR             *sap;   /* Section Attribute ptr      */
                MATERIAL_ATTR            *map;   /* Materials Attribute ptr    */
		FIBER_ELMT               *fep;   /* Fiber Element ptr          */
		} u;
	struct  Symbtab_element	*next;
} SYMBOL;

extern SYMBOL     *build_table();
extern SYMBOL     *install();
extern SYMBOL     *lookup();
extern void        print_symtable();

#endif /* end case SYMBOL_H */
