/*
 *  ============================================================================= 
 *  ALADDIN Version 1.0 :
 *              units.h : Data Structures for Engineering Quantities
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
 *  Written by: Mark Austin and Xiaoguang Chen                      December 1995
 *  ============================================================================= 
 */

#ifndef UNITS_H
#define UNITS_H

/* Units Type */

#define  SI     100
#define  US     200
#define  SI_US  300

/* Variables Defined for Units Functions */

#define COLUMN       1
#define ROW          2
#define COL_UNITS    3
#define ROW_UNITS    4

/* Engineering Quantities */

typedef struct dimensional_exponents {
        char         *units_name; /* units name              */
        double      scale_factor; /* scale/conversion factor */
        double      length_expnt;
        double        mass_expnt;
        double        time_expnt;
        double        temp_expnt; /* temperature             */
        int           units_type; /* US or SI units          */
} DIMENSIONS;

typedef struct engineering_quantity {
        double       value;
        DIMENSIONS  *dimen;
} QUANTITY;

/* Units functions */

#ifdef __STDC__

extern int            SameUnits( DIMENSIONS *, DIMENSIONS * );
extern DIMENSIONS    *UnitsMult( DIMENSIONS *, DIMENSIONS * );
extern DIMENSIONS    *UnitsMultRep( DIMENSIONS *, DIMENSIONS *, DIMENSIONS * );
extern DIMENSIONS    *UnitsDiv( DIMENSIONS *, DIMENSIONS *, int );
extern DIMENSIONS    *UnitsDivRep( DIMENSIONS *, DIMENSIONS *, DIMENSIONS *, int );
extern DIMENSIONS    *UnitsPower( DIMENSIONS *, double , int );
extern DIMENSIONS    *UnitsPowerRep( DIMENSIONS *, DIMENSIONS *, double , int );
extern DIMENSIONS    *UnitsCopy( DIMENSIONS *, DIMENSIONS * );
extern DIMENSIONS    *ZeroUnits( DIMENSIONS * );
extern DIMENSIONS    *DefaultUnits( char * );
extern DIMENSIONS    *UnitsSimplify( DIMENSIONS * );
extern DIMENSIONS    *RadUnitsSimplify( DIMENSIONS * );

extern void           UnitsPrint( DIMENSIONS * );
extern int            UnitsLength( char *, ... );

extern DIMENSIONS    *UnitsTypeConvert( DIMENSIONS * , int );
extern DIMENSIONS     UnitsScaleConvert( DIMENSIONS , int );
extern double         ConvertTempUnits( char *, double, int );

extern DIMENSIONS    *BufferInit( int );
extern void           BufferPrint( char *, DIMENSIONS *, int );

int       SetUnitsOn();
int       SetUnitsOff();
int       CheckUnits();
int       CheckUnitsType();

/*----------------------------------------------*/

extern QUANTITY          *QuantityUnitsLess( QUANTITY * );

/*----------------------------------------------*/

#else  /* start case not STDC */

extern int            SameUnits();
extern DIMENSIONS    *UnitsMult();
extern DIMENSIONS    *UnitsMultRep();
extern DIMENSIONS    *UnitsDiv();
extern DIMENSIONS    *UnitsDivRep();
extern DIMENSIONS    *UnitsPower();
extern DIMENSIONS    *UnitsPowerRep();
extern DIMENSIONS    *UnitsCopy();
extern DIMENSIONS    *ZeroUnits();
extern DIMENSIONS    *DefaultUnits();
extern DIMENSIONS    *UnitsSimplify();
extern DIMENSIONS    *RadUnitsSimplify();

extern void           UnitsPrint();
extern int            UnitsLength();

extern DIMENSIONS    *UnitsTypeConvert();
extern DIMENSIONS     UnitsScaleConvert();
extern double         ConvertTempUnits();

extern DIMENSIONS    *BufferInit();
extern void           BufferPrint();

int       SetUnitsOn();
int       SetUnitsOff();
int       CheckUnits();
int       CheckUnitsType();

/*----------------------------------------------*/

extern QUANTITY          *QuantityUnitsLess();

/*----------------------------------------------*/

#endif /* end case not STDC */
#endif /* end case UNITS_H */
