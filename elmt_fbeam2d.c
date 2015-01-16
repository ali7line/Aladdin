/*
 *  ============================================================================= 
 *  ALADDIN Version 2.0 :
 *                                                                     
 *  elmt_fbeam2d.c : Two-dimensional Flexible Beam-Column Element
 *                                                                     
 *  Copyright (C) 1995 by Mark Austin, Xiaoguang Chen, and Wane-Jang Lin
 *  Institute for Systems Research,                                           
 *  University of Maryland, College Park, MD 20742                                   
 *                                                                     
 *  This software is provided "as is" without express or implied warranty.
 *  Permission is granted to use this software for any purpose, and on any
 *  computer system, and to redistribute it freely, subject to the following
 *  restrictions:
 * 
 *  1. The authors are not responsible for the consequences of use of
 *     this software, even if they arise from defects in the software.
 *  2. The origin of this software must not be misrepresented, either
 *     by explicit claim or by omission.
 *  3. Altered versions must be plainly marked as such, and must not
 *     be misrepresented as being the original software.
 *  4. This notice is to remain intact.
 *                                                                    
 *  -------------------------------------------------------------------
 *  Convention for Nodal Forces
 *             +ve M     -  anticlockwise(RHT Rule)
 *             +ve X,Y,Z -  along +ve axis
 *  Convention for Member End  Forces
 *             +ve M     -  Sagging Moments
 *             +ve SF    -  LHS upwards
 *             +ve AF    -  Tension(LHS outwards)
 *  -------------------------------------------------------------------
 *                                                                    
 *  Written by: Ji Zhai                                             April, 1998
 *  ============================================================================= 
 */

#include <math.h>
#include "defs.h"
#include "miscellaneous.h"
#include "units.h"
#include "matrix.h"
#include "fe_database.h"
#include "symbol.h"
#include "fe_functions.h"
#include "elmt.h"



/* 
 *  ============================================================== 
 *  Input Properties:                                      
 *
 *       p->work_material[0] = E;
 *       p->work_material[1] = G;
 *       p->work_material[2] = fy;
 *       p->work_material[3] = ET;
 *       p->work_material[4] = nu;
 *       p->work_material[5] = density;
 *       p->work_material[6] = fu;
 *
 *       p->work_section[0] = Ixx;
 *       p->work_section[1] = Iyy;
 *       p->work_section[2] = Izz;
 *       p->work_section[3] = Ixy;
 *       p->work_section[4] = Ixz;
 *       p->work_section[5] = Iyz;
 *       p->work_section[6] = weight;
 *       p->work_section[7] = bf;
 *       p->work_section[8] = tf;
 *       p->work_section[9] = depth;                                  
 *       p->work_section[10] = area;
 *       p->work_section[11] = plate_thickness;
 *       p->work_section[12] = J;
 *       p->work_section[13] = rT;
 *       p->work_section[14] = width;
 *       p->work_section[15] = tw;                                  
 * 
 *  Input  :  ARRAY *p  -- pointer to working ARRAY data structure
 *         :  int isw   -- flag for task to be computed.
 *  Output :  ARRAY *p  -- pointer to working ARRAY data structure
 *  ============================================================== 
 */

extern EFRAME             *frame;

#ifdef __STDC__
ARRAY *elmt_fbeam_2d(ARRAY *p, int isw)
#else
ARRAY *elmt_fbeam_2d(p, isw)
ARRAY *p;
int     isw;
#endif
{
static double nu; 
static QUANTITY fy, G, E, ET, density; 
static double  Ixx, Iyy, Izz, Ixy, Ixz, Iyz, bf, tf, A, depth, weight, EA, EIzz;

static int       no_integ_pt, elmt_no;
static QUANTITY  elmt_length;

double cs, sn, tn, xl, xx, yy, eps, chi, gam, mzc, mass;
double vl1,vl2,tl1,tl2, fx1,fx2,fy1,fy2,mz1,mz2,e6,e2;
double  eplas , alp, sig;
int NS_Sub_Incr;

DIMENSIONS *dp_length, *dp_force, *dp_moment;               /*  delete */
DIMENSIONS *dp_stress, *dp_degree, *dp_temperature;         /* delete  */


double sum, temp;
int    i, j, k, ii, jj, iGaussP;
int    UNITS_SWITCH, UnitsType;
int node_no;

double      *abscissas, *weights;
double      xi;
double      scale;

ELEMENT     *ep;

MATRIX      *Q, *q;
MATRIX      *F, *K, *Ke;
MATRIX      *L, *Ltrans, *R, *Rtrans, *temp_m1, *temp_m2, *temp_m3;
MATRIX      *kx1, *kx2, *kt, *ft;
MATRIX      *bx, *bxtrans;
    /* [a] : Jump to task case */

    UNITS_SWITCH = CheckUnits();
    UnitsType = CheckUnitsType();
    switch(isw) {
        case PROPTY: /* beam element :   material properties  */
             E.value       =  p->work_material[0].value;
             nu            =  p->work_material[4].value;
             density.value =  p->work_material[5].value;
             if( UNITS_SWITCH == ON ) {
                 E.dimen       =  p->work_material[0].dimen;
                 density.dimen =  p->work_material[5].dimen;
             }

             /* [a] : check  poi_ratio value */

             if( nu == 0.0 || nu > 0.5 ) {
                 printf("WARNING >> ... In 2d frame elmt () \n");
                 printf("WARNING >> ... nu value = %9.4f,reset to 0.3 !\n", nu);
                 nu = 0.3;    /* default poi_ratio value */
             }

             /* [b] : calculate  G value */
            
             G.value = p->work_material[1].value = E.value/(1.0 -2.0 *  nu) ;
             if( UNITS_SWITCH == ON )  G.dimen = E.dimen;

             if(E.value/((1.0 - 2.0*nu)) != p->work_material[1].value) {
                 printf("WARNING >> G is not equal to E/(1-2nu), check G for homogeneous material \n");
                 printf("WARNING >> Ignore this message for non-homogeneous materials \n");
             }

             Ixx    = p->work_section[0].value;
             Iyy    = p->work_section[1].value;
             Izz    = p->work_section[2].value;
             Ixy    = p->work_section[3].value;
             Ixz    = p->work_section[4].value;
             Iyz    = p->work_section[5].value;
             weight = p->work_section[6].value;
             bf     = p->work_section[7].value;
             tf     = p->work_section[8].value;
             depth  = p->work_section[9].value;
             A      = p->work_section[10].value;

             EA     = E.value*A;
             EIzz   = E.value*Izz;

	     /* Some element qualtity */
	      no_integ_pt = p->integ_ptr->integ_pts;
	      elmt_no     = p->elmt_no;
	      cs = p->coord[0][1].value - p->coord[0][0].value;
	      sn = p->coord[1][1].value - p->coord[1][0].value;
	      elmt_length.value = sqrt( cs*cs + sn*sn );
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
#ifdef DEBUGZJ
	  printf("***** Entering STRESS or LOAD in fbeam ******\n");
#endif

             if( UNITS_SWITCH == ON )
                 SetUnitsOff();

	     ep = &frame->element[elmt_no - 1];

	     if( p->elmt_state == 0 ) {

	       cs = p->coord[0][1].value - p->coord[0][0].value;
	       sn = p->coord[1][1].value - p->coord[1][0].value;
	       xl = sqrt(cs * cs + sn * sn);
	       cs = cs/xl; 
	       sn = sn/xl;
	       p->length.value = xl;
	       elmt_length.value = xl;
          
	       abscissas = (double *)MyCalloc( no_integ_pt + 1, sizeof(double) );   /* bad design of gauss ==> +1 */
	       weights   = (double *)MyCalloc( no_integ_pt + 1, sizeof(double) );

	       gauss( abscissas, weights, no_integ_pt ); 

	       kx1 = MatrixAllocIndirect( "kx1", DOUBLE_ARRAY, 6, 6);

	       /* Gauss integeration for the tangent stiffess --- Eq.(A.10a) */
 
	       for( iGaussP = 1; iGaussP < no_integ_pt + 1; iGaussP++)
		 {  
		   scale = 0.5 * elmt_length.value * weights[iGaussP];
		   T_Stiffness_2d(p, kx1, abscissas[iGaussP], elmt_no );       /* for K at integeration points
                                                                           --- Eq.(A.10a) */ 
		   if( iGaussP == 1 )
		     K = MatrixScale(kx1, scale);
		   else 
		     {
		       temp_m1 = MatrixScale(kx1, scale);
		       MatrixAddReplace( K, temp_m1);
		       MatrixFree(temp_m1);
		     }
		 }
#ifdef DEBUGZJ
	     for( ii = 0; ii < 6; ii++)
		 for( jj = 0; jj < 6; jj++)
		   printf("Kt[%d][%d] = %g\n", ii, jj, K->uMatrix.daa[ii][jj] );
#endif 
	     MatrixFree( kx1 );
	     free((char *) abscissas);
	     free((char *) weights);
	     
	       
	       /* put displacement matrix to one column form */
	       /* p->displ=[px1,px2; py1,py2; pz1,pz2]3x2    */
	       /* => temp_m1=pe=[px1;py1;pz1;px2;py2;pz2]6x1 */

	       temp_m1  = MatrixAllocIndirect( (char *)NULL, DOUBLE_ARRAY, 6, 1 );
	       for( ii=0 ; ii < p->dof_per_node ; ii++ )
		 for( jj=0 ; jj < p->nodes_per_elmt ; jj++ ){
		   temp_m1->uMatrix.daa[ii+jj*p->dof_per_node][0] = 
		                                    p->displ->uMatrix.daa[ii][jj];

		 }
	       L = Local_To_Global_2d( p->coord, elmt_length );
	       q = MatrixMult( L, temp_m1 );            
	       Q = MatrixMult( K, q );           /* Nodal force in local coordiante system */
	       Q = MatrixMult( K, temp_m1 );
	       MatrixFree( temp_m1 );
	     }
	     else if( p->elmt_state == 1 ) {
	       L = Local_To_Global_2d( p->coord, elmt_length );
	       Q = p->Q_saved;
	     }

	     if( isw == LOAD_MATRIX ) {
	       Ltrans = MatrixTranspose( L );       
	       temp_m1 = MatrixMult( Ltrans, Q );   /* element nodal force in global coord. */
	       
	  
#ifdef DEBUGZJ1
	  
	       temp_m2 = MatrixMult( Ltrans, K );   
	       temp_m3 = MatrixMult( temp_m2, L );

	       /*  for( ii = 0; ii < 6; ii++)
		 for( jj = 0; jj < 6; jj++)
		   printf("L[%d][%d] = %g\n", ii, jj, L->uMatrix.daa[ii][jj] );
		   */
	       for( ii = 0; ii < 6; ii++)
		 for( jj = 0; jj < 6; jj++)
		   printf("K[%d][%d] = %g\n", ii, jj, temp_m3->uMatrix.daa[ii][jj] );

	       MatrixFree( temp_m2 );
	       MatrixFree( temp_m3 );
	   
	       for( ii =0; ii < 6; ii++)		     
		 printf("Q[%d][0] = %g   temp_m1[%d][0] = %g\n", 
		     ii, Q->uMatrix.daa[ii][0],
		     ii, temp_m1->uMatrix.daa[ii][0]);	     
#endif
	       MatrixFree( Ltrans );
	     }
          
	     if( isw == STRESS ) {
	       R = Rigid_Body_Rotation_2d( elmt_length );
	       Rtrans = MatrixTranspose( R );
	       temp_m1 = MatrixMult( Rtrans, Q );   /* element nodal force in local coord. */
	       MatrixFree( R );
	       MatrixFree( Rtrans );
	     }
	     MatrixFree( L );

	     /* Assign force values */
	     for( ii=0 ; ii < p->size_of_stiff ; ii++ )
	       /*  p->nodal_loads[ii].value = temp_m1->uMatrix.daa[ii][0]; */
	       p->nodal_loads[ii].value = Q->uMatrix.daa[ii][0];

	     if( p->elmt_state == 0 ) {
	       MatrixFree( K );
	       MatrixFree( Q );
	       MatrixFree( q );
	     }
	       MatrixFree( temp_m1 );
	     
	     /* Assign force units */

	     if( UNITS_SWITCH == ON ) {
	       SetUnitsOn();
	       switch( UnitsType ) {
                case SI:
                case SI_US:
                   dp_force  = DefaultUnits("N");
                   dp_length = DefaultUnits("m");
                   break;
                case US:
                   dp_force  = DefaultUnits("lbf");
                   dp_length = DefaultUnits("in");
                   break;
	       }

	       /* node 1 */

	       UnitsCopy( p->nodal_loads[0].dimen, dp_force );
	       UnitsCopy( p->nodal_loads[1].dimen, dp_force );
	       UnitsMultRep( p->nodal_loads[2].dimen, dp_force, dp_length );
	       
	       /* node 2 */
	       
	       UnitsCopy( p->nodal_loads[3].dimen, p->nodal_loads[0].dimen );
	       UnitsCopy( p->nodal_loads[4].dimen, p->nodal_loads[1].dimen );
	       UnitsCopy( p->nodal_loads[5].dimen, p->nodal_loads[2].dimen );

	       free((char *) dp_force->units_name);
	       free((char *) dp_length->units_name);
	       free((char *) dp_force);
	       free((char *) dp_length);
	     }

          if(isw == STRESS ) {
             xx  = 0.5*(p->coord[0][0].value + p->coord[0][1].value);   /* xx = 0.5(x1+x2) */
             yy  = 0.5*(p->coord[1][0].value + p->coord[1][1].value);   /* yy = 0.5(y1+y2) */
             printf("\n");
             printf("Elmt No %3d : \n", p->elmt_no);
             switch( UNITS_SWITCH ) {
                case ON :
                   printf("Coords (X,Y) = (%10.3f %s, %10.3f %s)\n", 
                       xx/elmt_length.dimen->scale_factor,elmt_length.dimen->units_name,
                       yy/elmt_length.dimen->scale_factor,elmt_length.dimen->units_name);
                   printf("\n");

                   /* node i */
                   printf(" Fx1 = %13.5e %s  Fy1 = %13.5e %s  Mz1 = %13.5e %s\n",
                       p->nodal_loads[0].value/p->nodal_loads[0].dimen->scale_factor,
                       p->nodal_loads[0].dimen->units_name,
                       p->nodal_loads[1].value/p->nodal_loads[1].dimen->scale_factor,
                       p->nodal_loads[1].dimen->units_name,
                       p->nodal_loads[2].value/p->nodal_loads[2].dimen->scale_factor,
                       p->nodal_loads[2].dimen->units_name);

                   /* node j */
                   printf(" Fx2 = %13.5e %s  Fy2 = %13.5e %s  Mz2 = %13.5e %s\n",
                       p->nodal_loads[3].value/p->nodal_loads[3].dimen->scale_factor,
                       p->nodal_loads[3].dimen->units_name,
                       p->nodal_loads[4].value/p->nodal_loads[4].dimen->scale_factor,
                       p->nodal_loads[4].dimen->units_name,
                       p->nodal_loads[5].value/p->nodal_loads[5].dimen->scale_factor,
                       p->nodal_loads[5].dimen->units_name);
                       printf("\n");
                   /* Member Forces */
                   printf(" Axial Force : x-direction = %13.5e %s\n",
                       -p->nodal_loads[0].value/p->nodal_loads[0].dimen->scale_factor,
                        p->nodal_loads[0].dimen->units_name);
                   printf(" Shear Force : y-direction = %13.5e %s\n",
                        p->nodal_loads[1].value/p->nodal_loads[1].dimen->scale_factor,
                        p->nodal_loads[1].dimen->units_name);
                   printf("\n");
                   break;
                case OFF :
                   printf("Coords (X,Y) = (%10.3f , %10.3f )\n", xx, yy);
                   printf("\n");

                   /* node i */
                   printf(" Fx1 = %13.5e   Fy1 = %13.5e   Mz1 = %13.5e \n",
                       p->nodal_loads[0].value, 
                       p->nodal_loads[1].value, 
                       p->nodal_loads[2].value);
                   /* node j */
                   printf(" Fx2 = %13.5e   Fy2 = %13.5e   Mz2 = %13.5e \n",
                       p->nodal_loads[3].value,
                       p->nodal_loads[4].value,
                       p->nodal_loads[5].value);
                   printf("\n");
                   /* Member Forces */
                   printf(" Axial Force : x-direction = %13.5e \n", -p->nodal_loads[0].value);
                   printf(" Shear Force : y-direction = %13.5e \n",  p->nodal_loads[1].value);
                   printf("\n");
                   break;
                default:
                   break;
             }
          }
#ifdef DEBUGZJ
	  printf("***** Ending STRESS or LOAD in fbeam ******\n");
#endif
          break;  /* end of case STRESS, LOAD_MATRIX */

  
        case STIFF: /* form element stiffness */
#ifdef DEBUGZJ
	  printf("***** Entering STIFF in fbeam ******\n");
#endif
	     if( UNITS_SWITCH == ON )
                 SetUnitsOff();

	     ep = &frame->element[elmt_no - 1];

#ifdef DEBUGZJ
	     printf("elmt no=%d   x1=%g  x2=%g  y1=%g  y2=%g \n", 
		    p->elmt_no,
		    p->coord[0][0].value,
		    p->coord[0][1].value,
		    p->coord[1][0].value,
		    p->coord[1][1].value );
#endif
	     cs = p->coord[0][1].value - p->coord[0][0].value;
             sn = p->coord[1][1].value - p->coord[1][0].value;
             xl = sqrt(cs * cs + sn * sn);
             cs = cs/xl; 
             sn = sn/xl;
             p->length.value = xl;
             elmt_length.value = xl;
          
             abscissas = (double *)MyCalloc( no_integ_pt + 1, sizeof(double) );   /* bad design of gauss ==> +1 */
             weights   = (double *)MyCalloc( no_integ_pt + 1, sizeof(double) );

             gauss( abscissas, weights, no_integ_pt ); 

             kx1 = MatrixAllocIndirect( "kx1", DOUBLE_ARRAY, 6, 6);

           /* Gauss integeration for the tangent stiffness --- Eq.(A.10a) */
 
 	   for( iGaussP = 1; iGaussP < no_integ_pt + 1; iGaussP++)
	    {  
             scale = 0.5 * elmt_length.value * weights[iGaussP];
             T_Stiffness_2d(p, kx1, abscissas[iGaussP], elmt_no );       /* for K at integeration points
                                                                           --- Eq.(A.10a) */ 
	     if( iGaussP == 1 )
                K = MatrixScale(kx1, scale);
             else 
               {
                temp_m1 = MatrixScale(kx1, scale);
                MatrixAddReplace( K, temp_m1);
                MatrixFree(temp_m1);
               }
            }
#ifdef DEBUGZJ
	     for( ii = 0; ii < 6; ii++)
		 for( jj = 0; jj < 6; jj++)
		   printf("Kt[%d][%d] = %g\n", ii, jj, K->uMatrix.daa[ii][jj] );
#endif
         
          MatrixFree( kx1 );
         
          free((char *) abscissas);
          free((char *) weights);
         
         
          /* Transform local coordinate to global */
	  L = Local_To_Global_2d( p->coord, elmt_length );   
          Ltrans = MatrixTranspose( L );
          temp_m1 = MatrixMult( Ltrans, K );   
          Ke = MatrixMult( temp_m1, L );
	  
          /* Copy stiffness matrix to p array */

          for(ii = 1; ii <= p->stiff->iNoRows; ii++)
	    for(jj = 1; jj <= p->stiff->iNoColumns; jj++){
                p->stiff->uMatrix.daa[ii-1][jj-1] = K->uMatrix.daa[ii-1][jj-1];
#ifdef DEBUGZJ1
		printf("Ke[%d][%d]= %g\n", ii-1, jj-1, Ke->uMatrix.daa[ii-1][jj-1] );
#endif
	    }

          MatrixFree( K );
          MatrixFree( temp_m1 );
	  MatrixFree( L );
          MatrixFree( Ltrans ); 
          MatrixFree( Ke );
	 
          /* Assign units to p array stiffness */
          if( UNITS_SWITCH == ON ) {
             SetUnitsOn();
             switch( UnitsType ) {
                case SI:
                case SI_US:
                   dp_force  = DefaultUnits("N");
                   dp_length = DefaultUnits("m");
                   break;
                case US:
                   dp_force  = DefaultUnits("lbf");
                   dp_length = DefaultUnits("in");
                   break;
             }

             ZeroUnits( &(p->stiff->spRowUnits[0]) );
             ZeroUnits( &(p->stiff->spRowUnits[1]) );
             UnitsCopy( &(p->stiff->spRowUnits[2]), dp_length );
             UnitsCopy( &(p->stiff->spRowUnits[3]), &(p->stiff->spRowUnits[0]) );
             UnitsCopy( &(p->stiff->spRowUnits[4]), &(p->stiff->spRowUnits[1]) );
             UnitsCopy( &(p->stiff->spRowUnits[5]), &(p->stiff->spRowUnits[2]) );

             UnitsDivRep( &(p->stiff->spColUnits[0]), dp_force, dp_length, NO );
             UnitsCopy(   &(p->stiff->spColUnits[1]), &(p->stiff->spColUnits[0]) );
             UnitsCopy(   &(p->stiff->spColUnits[2]), dp_force );
             UnitsCopy( &(p->stiff->spColUnits[3]), &(p->stiff->spColUnits[0]) );
             UnitsCopy( &(p->stiff->spColUnits[4]), &(p->stiff->spColUnits[1]) );
             UnitsCopy( &(p->stiff->spColUnits[5]), &(p->stiff->spColUnits[2]) );

             free((char *) dp_force->units_name);
             free((char *) dp_length->units_name);
             free((char *) dp_force);
             free((char *) dp_length);
          
         } /* end of units on/off for case STIFF */
#ifdef DEBUGZJ
	  printf("***** Ending STIFF in fbeam ******\n");
#endif            
             break;
        case MASS_MATRIX:

             cs = p->coord[0][1].value - p->coord[0][0].value;
             sn = p->coord[1][1].value - p->coord[1][0].value;
             xl = sqrt(cs * cs + sn * sn);
             cs = cs/xl; 
             sn = sn/xl;
             p->length.value = xl;

             /* Calculate mass = m_bar = mass/length                   */
             /* in units of (kg/m) or ((lbf-sec^2/in)/in)=(lb/in)      */
             /* if no units, then assume gravity g = 9.80665 m/sec^2   */

             if( weight != 0.0 )
               mass = weight/9.80665;
             else
               if( density.value > 0 )  mass = A * density.value ;
             else {
                  printf("ERROR >> In input file \n");
                  printf("ERROR >> You need a density value to calculate mass matrix\n");
                  exit(1);
             }

             if( UNITS_SWITCH == ON ) 
                p->stiff->spColUnits[0].units_type = CheckUnitsType();

             p->stiff = beamms(p, p->stiff, p->type, mass, xl, cs, sn,
                               p->size_of_stiff, p->dof_per_node);
             break;
        default:
             break;
    }

    return(p);
}

/*
 *  ==============================================================================
 *  T_Stiffness_2d(): Compute flexible beam stiffness at integer point 
 *  ==============================================================================
 */

#ifdef  __STDC__
void T_Stiffness_2d(ARRAY *p, MATRIX *kx, double xi, int elmt_no )
#else
void T_Stiffness_2d(p, kx, xi, elmt_no )
ARRAY * p;
MATRIX *kx;
double xi;
int elmt_no;
#endif
{
static QUANTITY  elmt_length;

 ELEMENT  *ep;
 int	  i, j;
 int      total_fbeam;
 double   y, z, A, E, G;
 double   du1, du2, dtheta, dN1, dN2;
 double   cs, sn;
 double   theta;
 MATRIX   *B, *B2, *N, *C, *T,*n, *Gmatrix;
 MATRIX   *temp_1, *temp_2, *temp_3, *temp_4, *temp_5, *temp_6, *temp_7, 
          *temp_8, *temp_9, *temp_10, *temp_11, *temp_12, *temp_13;

/* Assign the node displacements by using frame ( global to local) */
 ep = &frame->element[p->elmt_no - 1];
 
 for( i=0 ; i < kx->iNoRows ; i++ )
   for( j=0 ; j < kx->iNoColumns ; j++ )
     kx->uMatrix.daa[i][j] = 0.0;                   /* initiate the Matrice */
    
 cs = p->coord[0][1].value - p->coord[0][0].value;
 sn = p->coord[1][1].value - p->coord[1][0].value;
 elmt_length.value = sqrt( cs*cs + sn*sn );
 
 B = MatrixAllocIndirect(   "B",  DOUBLE_ARRAY, 3, 6);
 B2 = MatrixAllocIndirect(  "B2",  DOUBLE_ARRAY, 3, 6);
 N = MatrixAllocIndirect(   "N",  DOUBLE_ARRAY, 2, 1);
 C = MatrixAllocIndirect(   "C",  DOUBLE_ARRAY, 3, 3);
 T = MatrixAllocIndirect(   "T",  DOUBLE_ARRAY, 3, 3);
 Gmatrix = MatrixAllocIndirect(   "Gmatrix",  DOUBLE_ARRAY, 3, 3); 
 temp_1 = MatrixAllocIndirect(   "temp_1",  DOUBLE_ARRAY, 3, 1); 
 
 for( i=0 ; i < T->iNoRows ; i++ )
   for( j=0 ; j < T->iNoColumns ; j++ )
     T->uMatrix.daa[i][j] = 0.0;
 
 for( i=0 ; i < C->iNoRows ; i++ )
   for( j=0 ; j < C->iNoColumns ; j++ )
     C->uMatrix.daa[i][j] = 0.0;
 
 for( i=0 ; i < B->iNoRows ; i++ )
   for( j=0 ; j < B->iNoColumns ; j++ )
     B->uMatrix.daa[i][j] = 0.0;
 
 for( i=0 ; i < B2->iNoRows ; i++ )
   for( j=0 ; j < B2->iNoColumns ; j++ )
     B2->uMatrix.daa[i][j] = 0.0;
 
 for( i=0 ; i < Gmatrix ->iNoRows ; i++ )
   for( j=0 ; j < Gmatrix ->iNoColumns ; j++ )
     Gmatrix ->uMatrix.daa[i][j] = 0.0;
 
 N -> uMatrix.daa[0][0] = -0.5 * ( xi - 1.0 );             /* shap function */
 N -> uMatrix.daa[1][0] =  0.5 * ( xi + 1.0 );
 
 /*theta = theta1*N1 + theta2*N2, here thetai is ith node theta*/
 theta = ( p -> displ -> uMatrix.daa[2][0] ) * ( N -> uMatrix.daa[0][0] )  + 
         ( p -> displ -> uMatrix.daa[2][1] ) * ( N -> uMatrix.daa[1][0] ) ;                
 
 /* Tansform matrix, Eq.(A.5b) */
 T -> uMatrix.daa[0][0] = cos( theta );                 
 T -> uMatrix.daa[0][1] = - sin( theta );
 T -> uMatrix.daa[1][0] = sin( theta );
 T -> uMatrix.daa[1][1] = cos( theta );
 T -> uMatrix.daa[2][2] = 1.0;
 
 /* elastic matrix, Eq.(A.5b) */
 
 C -> uMatrix.daa[0][0] = p -> work_material[0].value * p -> work_section[10].value;     /* EA */
 C -> uMatrix.daa[1][1] = p -> work_material[1].value * p -> work_section[10].value * 5.0 / 6.0 ;     /* GAs ?  Here is for retangular section ! ---- CHANGE it !!*/
 C -> uMatrix.daa[2][2] = p -> work_material[0].value * p -> work_section[2].value;     /* EI  */
 
 /* according to strain at xi, get du1/dx, du2/dx , dtheta/dx for the fbeam (isoparamater) */
 
 du1 = (p -> displ -> uMatrix.daa[0][1] - p -> displ -> uMatrix.daa[0][0]) / elmt_length.value;
 du2 = (p -> displ -> uMatrix.daa[1][1] - p -> displ -> uMatrix.daa[1][0]) / elmt_length.value;
 dtheta = (p -> displ -> uMatrix.daa[2][1] - p -> displ -> uMatrix.daa[2][0]) / elmt_length.value;
 
#ifdef DEBUGZJ
 printf("elm = %d  u11 = %g  u12 = %g  theta1 = %g  u21 = %g  u22 = %g  theta2 = %g\n", elmt_no, 
	p -> displ -> uMatrix.daa[0][0], 
	p -> displ -> uMatrix.daa[1][0], 
	p -> displ -> uMatrix.daa[2][0], 
	p -> displ -> uMatrix.daa[0][1], 
	p -> displ -> uMatrix.daa[1][1], 
	p -> displ -> uMatrix.daa[2][1] );
#endif
 /* form relative matrix */

 dN1 = - 1.0 / elmt_length.value;    /* dN1/dx */
 dN2 =   1.0 / elmt_length.value;    /* dN2/dx */
 
 /* form B = [D1][N] */
 
 B -> uMatrix.daa[0][0] = dN1;
 B -> uMatrix.daa[1][1] = dN1;
 B -> uMatrix.daa[2][2] = dN1;
 B -> uMatrix.daa[0][2] = du2 * N -> uMatrix.daa[0][0];
 B -> uMatrix.daa[1][2] = - ( 1.0 + du1 ) * N -> uMatrix.daa[0][0];
 B -> uMatrix.daa[0][3] = dN2;
 B -> uMatrix.daa[1][4] = dN2;
 B -> uMatrix.daa[2][5] = dN2;
 B -> uMatrix.daa[0][5] = du2 * N -> uMatrix.daa[1][0];
 B -> uMatrix.daa[1][5] = - ( 1.0 + du1 ) * N -> uMatrix.daa[1][0];
       
 /* form B2 = [D2][N] */

 B2 -> uMatrix.daa[0][0] = dN1;
 B2 -> uMatrix.daa[1][1] = dN1;
 B2 -> uMatrix.daa[2][2] = N -> uMatrix.daa[0][0];
 B2 -> uMatrix.daa[0][3] = dN2;
 B2 -> uMatrix.daa[1][4] = dN2;
 B2 -> uMatrix.daa[2][5] = N -> uMatrix.daa[1][0];

 /* form n matrix, Eq.(A.5a) */      
 
 temp_1 -> uMatrix.daa[0][0] = 1.0 + du1;
 temp_1 -> uMatrix.daa[1][0] = du2;
 temp_1 -> uMatrix.daa[2][0] = dtheta;
       
 temp_2 = MatrixTranspose( T );                  
 temp_9 = MatrixMult(temp_2, temp_1);
       
 temp_1 -> uMatrix.daa[0][0] = 1.0;
 temp_1 -> uMatrix.daa[1][0] = 0.0;
 temp_1 -> uMatrix.daa[2][0] = 0.0;

 MatrixSubReplace( temp_9, temp_1);
 temp_12 = MatrixMult(C, temp_9);
 n      = MatrixMult(T, temp_12);
       
 /* form the integrating function of K, Eq.(A.10b) */

 temp_3 = MatrixTranspose( B );
 temp_4 = MatrixMult(temp_3, T);
 temp_5 = MatrixMult(temp_4, C); 
 temp_6 = MatrixTranspose( T );
 temp_13 = MatrixMult(temp_5, temp_6);
 temp_10 = MatrixMult(temp_13, B);                   

 /* form G matrix, Eq.(A.8) */
      
 Gmatrix->uMatrix.daa[0][2] = - n -> uMatrix.daa[1][0];
 Gmatrix->uMatrix.daa[1][2] =   n -> uMatrix.daa[0][0];
 Gmatrix->uMatrix.daa[2][0] = - n -> uMatrix.daa[1][0];
 Gmatrix->uMatrix.daa[2][1] =   n -> uMatrix.daa[0][0];
 Gmatrix->uMatrix.daa[2][2] = - ( (1.0 + du1) * n -> uMatrix.daa[0][0] + 
				  du2 * n -> uMatrix.daa[1][0] );
      
 /* form the integrating function of Kg, Eq.(A.10c) */

 temp_5 = MatrixTranspose( B2 );
 temp_11 = MatrixMult(temp_5, Gmatrix);
 temp_7 = MatrixMult(temp_11, B2);                

 /* 
#ifdef DEBUG
 for( i=0 ; i < kx ->iNoRows ; ++i )
   for( j=0 ; j < kx ->iNoColumns ; ++j )
      printf("kg[%d][%d] = %g\n ", i, j, temp_7->uMatrix.daa[i][j] );
#endif
*/ 
 temp_8 = MatrixAdd(temp_10, temp_7);            /* Eq.A.10a */

 for( i=0 ; i < kx ->iNoRows ; ++i )
   for( j=0 ; j < kx ->iNoColumns ; ++j )
     kx->uMatrix.daa[i][j] = temp_8->uMatrix.daa[i][j];
 
 MatrixFree( B );
 MatrixFree( B2 );
 MatrixFree( n );
 MatrixFree( N );
 MatrixFree( C );
 MatrixFree( T );
 MatrixFree( Gmatrix);
 MatrixFree( temp_1);
 MatrixFree( temp_2);
 MatrixFree( temp_3);
 MatrixFree( temp_4);
 MatrixFree( temp_5);
 MatrixFree( temp_6);
 MatrixFree( temp_7);
 MatrixFree( temp_8);
 MatrixFree( temp_9);
 MatrixFree( temp_10);
 MatrixFree( temp_11);
 MatrixFree( temp_12);
 MatrixFree( temp_13);
}


/*
 *  ===============================================================
 *  print_property_fbeam_2d() : print Flexible Beam Element Properties
 *
 *  Input  : EFRAME *frp  --
 *         : int i        -- 
 *  Output : void 
 *  ===============================================================
 */

#ifdef __STDC__
void print_property_fbeam_2d(EFRAME *frp, int i)
#else
void print_property_fbeam_2d(frp, i)
EFRAME    *frp;
int          i;                 /* elmt_attr_no */
#endif
{
int     UNITS_SWITCH;
ELEMENT_ATTR    *eap;

#ifdef DEBUG
       printf("*** Enter print_property_fbeam_2d()\n");
#endif

     UNITS_SWITCH = CheckUnits();
     eap = &frp->eattr[i-1];

     if( PRINT_MAP_DOF == ON ) {
        if(frp->no_dof == 3 || frp->no_dof == 2) { 
           printf("             ");
           printf("         : gdof [0] = %4d : gdof[1] = %4d : gdof[2] = %4d\n",
                           eap->map_ldof_to_gdof[0],
                           eap->map_ldof_to_gdof[1],
                           eap->map_ldof_to_gdof[2]);
        }

        if(frp->no_dof == 6) { /* 3d analysis */
           printf("             ");
           printf("         : dof-mapping : gdof[0] = %4d : gdof[1] = %4d : gdof[2] = %4d\n",
                           eap->map_ldof_to_gdof[0],
                           eap->map_ldof_to_gdof[1],
                           eap->map_ldof_to_gdof[2]);
           printf("             ");
           printf("                         gdof[3] = %4d : gdof[4] = %4d : gdof[5] = %4d\n",
                           eap->map_ldof_to_gdof[3],
                           eap->map_ldof_to_gdof[4],
                           eap->map_ldof_to_gdof[5]);
        } 
     }

     switch(UNITS_SWITCH) {
       case ON:
        UnitsSimplify( eap->work_material[0].dimen );
        UnitsSimplify( eap->work_material[2].dimen );
        UnitsSimplify( eap->work_material[5].dimen );
        UnitsSimplify( eap->work_section[2].dimen );
        UnitsSimplify( eap->work_section[10].dimen );
        if( eap->work_material[0].dimen->units_name != NULL ) {
           printf("             ");
           printf("         : Young's Modulus =  E = %16.3e %s\n",
                           eap->work_material[0].value/eap->work_material[0].dimen->scale_factor,
                           eap->work_material[0].dimen->units_name);
        }
        if( eap->work_material[4].value != 0.0 ) {
           printf("             ");
           printf("         : Poisson's ratio = nu = %16.3e   \n", eap->work_material[4].value);
        }
        if( eap->work_material[2].dimen->units_name != NULL ) {
           printf("             ");
           printf("         : Yielding Stress = fy = %16.3e %s\n",
                           eap->work_material[2].value/eap->work_material[2].dimen->scale_factor,
                           eap->work_material[2].dimen->units_name);
        }
	if( eap->work_material[5].dimen->units_name != NULL ) {
          printf("             ");
          printf("         : Density         = %16.3e %s\n",
                           eap->work_material[5].value/eap->work_material[5].dimen->scale_factor,
                           eap->work_material[5].dimen->units_name);
	}
	if( eap->work_section[2].dimen->units_name != NULL ) {
          printf("             ");
          printf("         : Inertia Izz     = %16.3e %s\n",
                           eap->work_section[2].value/eap->work_section[2].dimen->scale_factor,
                           eap->work_section[2].dimen->units_name);
	}
	if( eap->work_section[10].dimen->units_name != NULL ) {
          printf("             ");
          printf("         : Area            = %16.3e %s\n",
                           eap->work_section[10].value/eap->work_section[10].dimen->scale_factor,
                           eap->work_section[10].dimen->units_name);
	}
       break;
       case OFF:
        if( eap->work_material[0].value != 0.0 ) {
           printf("             ");
           printf("         : Young's Modulus =  E = %16.3e\n",
                            eap->work_material[0].value);
        }
        if( eap->work_material[2].value != 0.0 ) {
           printf("             ");
           printf("         : Yielding Stress = fy = %16.3e\n",
                            eap->work_material[2].value);
        }
        if( eap->work_material[4].value != 0.0 ) {
           printf("             ");
           printf("         : Poisson's ratio = nu = %16.3e   \n", eap->work_material[4].value);
        }
        if( eap->work_material[0].value != 0.0 ) {
           printf("             ");
           printf("         : Density         = %16.3e\n",
                            eap->work_material[5].value);
        }
        if( eap->work_section[2].value != 0.0 ) {
           printf("             ");
           printf("         : Inertia Izz     = %16.3e\n",
                            eap->work_section[2].value);
        }
        if( eap->work_section[10].value != 0.0 ) {
           printf("             ");
           printf("         : Area            = %16.3e\n",
                            eap->work_section[10].value);
        }
        break;
        default:
        break;
     }
#ifdef DEBUG
       printf("*** Leave print_property_fbeam_2d()\n");
#endif
}

/*
 *============================================================================
 * Local_To_Global_2d(): Ratation transformation matrix from local to global
 *                       coordinate system for FBEAM
 * Input: QUANTITY **coord ---------  coordates of nodes
 *        QUANTITY length  ---------  length of element
 *
 *============================================================================
 */
#ifdef  __STDC__
MATRIX *Local_To_Global_2d( QUANTITY **coord, QUANTITY length )
#else
MATRIX *Local_To_Global_2d( coord, length )
QUANTITY **coord;
QUANTITY length;
#endif
{
MATRIX *Lele;
double cs, sn;
double temp;
int    i, j, k;

 Lele  = MatrixAllocIndirect( (char *)NULL, DOUBLE_ARRAY, 6, 6 );

for( i=0 ; i < 6; i++ ) { 
  for( j=0 ; j < 6 ; j++ ) {
    Lele -> uMatrix.daa[i][j] = 0.0;
  }
}

Lele->uMatrix.daa[0][0] =  cs;
Lele->uMatrix.daa[1][1] =  cs;
Lele->uMatrix.daa[3][3] =  cs;
Lele->uMatrix.daa[4][4] =  cs;
Lele->uMatrix.daa[0][1] =  sn;
Lele->uMatrix.daa[1][0] = -sn;
Lele->uMatrix.daa[3][4] =  sn;
Lele->uMatrix.daa[4][3] = -sn;
Lele->uMatrix.daa[2][2] = 1.0; 
Lele->uMatrix.daa[5][5] = 1.0;         

return ( Lele );
}

