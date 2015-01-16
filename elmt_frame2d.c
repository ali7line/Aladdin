/*
 *  ============================================================================= 
 *  ALADDIN Version 1.0 :
 *       elmt_frame2d.c : Two-dimensional Beam-Column Element
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

/*
#define DEBUG 
*/


/* ============================================================== */
/*   Element FRAME-2D                                             */
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

#ifdef __STDC__
ARRAY *elmt_frame_2d(ARRAY *p, int isw)
#else
ARRAY *elmt_frame_2d(p, isw)
ARRAY *p;
int     isw;
#endif
{
static double nu; 
static QUANTITY fy, G, E, ET, density; 
static double  Ixx, Iyy, Izz, Ixy, Ixz, Iyz, bf, tf, A, depth, weight, EA, EIzz;

double cs, sn, xl,xx,yy,  eps,chi,gam, mzc,mass;
double vl1,vl2,tl1,tl2, fx1,fx2,fy1,fy2,mz1,mz2,e6,e2;
double  eplas , alp, sig;
int NS_Sub_Incr;
DIMENSIONS *dp_length, *dp_force, *dp_moment;
DIMENSIONS *dp_stress, *dp_degree, *dp_temperature;

double sum, temp;
int    i, j, k, ii;
int    UNITS_SWITCH;

#ifdef DEBUG
       printf("*** Enter elmt_frame_2d() : isw = %4d\n", isw);
#endif

    UNITS_SWITCH = CheckUnits();

    switch(isw) {
        case PROPTY:
             /* beam element :   material properties  */

             E.value       =  p->work_material[0].value;
             fy.value      =  p->work_material[2].value;
             ET.value      =  p->work_material[3].value;
             nu            =  p->work_material[4].value;
             density.value =  p->work_material[5].value;
             if( UNITS_SWITCH == ON ) {
                 E.dimen       =  p->work_material[0].dimen;
                 fy.dimen      =  p->work_material[2].dimen;
                 ET.dimen      =  p->work_material[3].dimen;
                 density.dimen =  p->work_material[5].dimen;
             }

           /* (1)   check  poi_ratio value */

            if( nu == 0.0 || nu > 0.5 ) {
                printf("WARNING >> ... In 2d beam element() - frame_2d -  nu value = %9.4f,reset to 0.3 !\n",
                       nu);
                nu = 0.3;    /* default poi_ratio value */
            }

            /* (2)   calculate  G value */
            
            G.value = p->work_material[1].value = E.value/(1.0 - 2.0*nu) ;
            if( UNITS_SWITCH == ON )  G.dimen = E.dimen;

            if(E.value/((1.0 - 2.0*nu)) != p->work_material[1].value) {
                    printf(" elmt_frame_2d(): WARNING: G is not equal to E/(1-2nu), check G for homogeneous material \n");
                    printf("                : ignore this message for non-homogeneous materials \n");
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

             cs = p->coord[0][1].value - p->coord[0][0].value;
             sn = p->coord[1][1].value - p->coord[1][0].value;

             xl = sqrt(cs * cs + sn * sn);
             cs = cs/xl; 
             sn = sn/xl;
             p->length.value = xl;

             xx  = 0.5*(p->coord[0][0].value + p->coord[0][1].value);            /* xx = 0.5(x1+x2) */
             yy  = 0.5*(p->coord[1][0].value + p->coord[1][1].value);            /* yy = 0.5(y1+y2) */
             eps = (cs*(p->displ->uMatrix.daa[0][1] - p->displ->uMatrix.daa[0][0])
                    +sn*(p->displ->uMatrix.daa[1][1] - p->displ->uMatrix.daa[1][0])) / xl;

                    /* eps = (u2-u1)cos(phi)+(v2-v1)sin(phi) */

             chi = (p->displ->uMatrix.daa[2][1] - p->displ->uMatrix.daa[2][0]) / xl;
             gam = -12 *(-sn*(p->displ->uMatrix.daa[0][0] - p->displ->uMatrix.daa[0][1])
                   + cs*(p->displ->uMatrix.daa[1][0]- p->displ->uMatrix.daa[1][1]))
                   / (xl*xl*xl)-6.*(p->displ->uMatrix.daa[2][0] + p->displ->uMatrix.daa[2][1])/xl/xl;

             /* local deflections */

             vl1 = -sn * p->displ->uMatrix.daa[0][0] + cs * p->displ->uMatrix.daa[1][0];
             vl2 = -sn * p->displ->uMatrix.daa[0][1] + cs * p->displ->uMatrix.daa[1][1];
             tl1 = p->displ->uMatrix.daa[2][0];
             tl2 = p->displ->uMatrix.daa[2][1];

             /* computer axial forces fx, transverse forces fy, and moment mz */

             e6 = 6* EIzz /xl/xl;
             e2 = 2* EIzz /xl;

             /* Elasto_Plastic load case                            */
             /* Considering uniaxial element with x_displ, Fx force */

             fx1 = - EA * eps;
             fy1 = 2*e6/xl*(vl1-vl2) + e6*(tl1 + tl2);
             fx2 = - fx1;
             fy2 = - fy1;
             mz1 = e6 * ( vl1 - vl2) + e2 * ( 2 * tl1 + tl2 );
             mz2 = e6 * ( vl1 - vl2) + e2 * (     tl1 + 2 * tl2 );
             mzc = (mz1 - mz2)/2;

             /* Add FEF if elmt loaded */
 
             if( p->elmt_load_ptr != NULL) { 
                 p = sld07(p, STRESS);

                 /* Add FEF to joint forces */

                 fx1 = fx1  - p->nodal_loads[0].value;
                 fy1 = fy1  - p->nodal_loads[1].value;
                 mz1 = mz1  - p->nodal_loads[2].value;
                 fx2 = fx2  - p->nodal_loads[3].value;
                 fy2 = fy2  - p->nodal_loads[4].value;
                 mz2 = mz2  - p->nodal_loads[5].value;
             } 

             /* Assign force values */

             if( UNITS_SWITCH == ON ) {
                if(UNITS_TYPE == SI) {
                    dp_length = DefaultUnits("m");
                    dp_force  = DefaultUnits("N");
                }
                else if(UNITS_TYPE == US) {
                    dp_length = DefaultUnits("in");
                    dp_force  = DefaultUnits("lbf");
                }
            
                /* node no 1 */
                UnitsCopy( p->nodal_loads[0].dimen, dp_force );
                UnitsCopy( p->nodal_loads[1].dimen, dp_force );
                UnitsMultRep( p->nodal_loads[2].dimen, dp_force, dp_length );

                /* node no = 2 */
                UnitsCopy( p->nodal_loads[3].dimen, p->nodal_loads[0].dimen );
                UnitsCopy( p->nodal_loads[4].dimen, p->nodal_loads[1].dimen );
                UnitsCopy( p->nodal_loads[5].dimen, p->nodal_loads[2].dimen );
             }

             p->nodal_loads[0].value = fx1; 
             p->nodal_loads[1].value = fy1; 
             p->nodal_loads[2].value = mz1;
             p->nodal_loads[3].value = fx2; 
             p->nodal_loads[4].value = fy2; 
             p->nodal_loads[5].value = mz2; 

             if(isw == LOAD_MATRIX ) {
                p->nodal_loads[0].value = fx1*cs - fy1*sn; 
                p->nodal_loads[1].value = fx1*sn + fy1*cs; 
                p->nodal_loads[2].value = mz1;
                p->nodal_loads[3].value = fx2*cs - fy2*sn; 
                p->nodal_loads[4].value = fx2*sn + fy2*cs; 
                p->nodal_loads[5].value = mz2; 
             }

             if(isw == STRESS && PRINT_STRESS == ON) {
                printf("\n");
                printf("Elmt No %3d : \n", p->elmt_no);
                if( UNITS_SWITCH == ON ) {
                    printf("Coords (X,Y) = (%10.3f %s, %10.3f %s)\n", 
                        xx/dp_length->scale_factor,dp_length->units_name,
                        yy/dp_length->scale_factor,dp_length->units_name);

                    printf("exx = %13.5e , curva = %13.5g , gamma = %13.5e\n", eps,chi, gam);
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

                    free((char *) dp_length->units_name);
                    free((char *) dp_length);
                    free((char *) dp_force->units_name);
                    free((char *) dp_force);
                } else {
                    printf("Coords (X,Y) = (%10.3f , %10.3f )\n", xx, yy);
                    printf("exx = %13.5e , curva = %13.5g , gamma = %13.5e\n", eps,chi, gam);
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
                }
             }
             break;
        case STIFF: /* form element stiffness */

#ifdef DEBUG
       printf("*** In elmt_frame_2d() : start case STIFF\n", isw);
#endif

             cs = p->coord[0][1].value - p->coord[0][0].value;
             sn = p->coord[1][1].value - p->coord[1][0].value;
             xl = sqrt(cs * cs + sn * sn);
             cs = cs/xl; 
             sn = sn/xl;
             p->length.value = xl;

             if( UNITS_SWITCH==ON )
                  p->stiff->spColUnits[0].units_type = UNITS_TYPE;

             if(p->nodes_per_elmt == 2) /* elastic elment use 2-node elmt */
                p->stiff = beamst(p, p->stiff, EA, EIzz, xl, cs, sn,
                                  p->size_of_stiff, p->dof_per_node);

#ifdef DEBUG
      for(i = 1; i <= p->size_of_stiff; i++) {
         sum = 0.0;
         temp = 0.0;
         for(j = 1; j <= p->dof_per_node; j++) {
            for(k = 1; k <= p->nodes_per_elmt; k++) {
                ii = p->dof_per_node*(k-1)+j;
                if(j <= 2)
                  sum += p->stiff->uMatrix.daa[i-1][ii-1];
                if(j > 1)
                  temp += p->stiff->uMatrix.daa[i-1][ii-1];
            }
         }
         printf(" Force[%d] = %lf\n", i, sum);
         printf(" Moment[%d] = %lf\n", i, temp);
      }
#endif

             

#ifdef DEBUG
       MatrixPrintIndirectDouble(p->stiff);
       printf("*** In elmt_frame_2d() : end case STIFF\n", isw);
#endif

             break;
        case MASS_MATRIX:

#ifdef DEBUG
       printf("*** In elmt_frame_2d() : start case MASS\n", isw);
       printf("                : EA      = %8.2f\n", EA);
       printf("                : EI      = %8.2f\n", EIzz);
       printf("                : Density = %8.2f\n", p->work_material[5].value);
#endif


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
                  printf("\nError in input: Need density value to calculate mass matrix\n");
                  exit(1);
             }

             if( UNITS_SWITCH == ON ) 
                p->stiff->spColUnits[0].units_type = UNITS_TYPE;

             p->stiff = beamms(p, p->stiff, p->type, mass, xl, cs, sn,
                               p->size_of_stiff, p->dof_per_node);

#ifdef DEBUG
       printf("mass\n");
       MatrixPrintIndirectDouble( p->stiff);
       printf("*** In elmt_frame_2d() : end case MASS\n", isw);
#endif
             break;
        default:
             break;
    }
#ifdef DEBUG
       printf("*** leaving elmt_frame_2d() \n");
#endif

    return(p);
}


/* ==================== */
/* Beam Stiffness       */
/* ==================== */

#ifdef __STDC__
MATRIX *beamst(ARRAY *p, MATRIX *s, double EA, double EI, double length, double cs, double sn, int size_of_stiff, int no_dof)
#else
MATRIX *beamst(p, s, EA, EI, length, cs, sn, size_of_stiff, no_dof)
ARRAY  *p;
MATRIX *s;
double  EA, EI, length, cs, sn;
int     size_of_stiff, no_dof ;
#endif
{
int             i, j, k;
int    length1, length2;
double                t;
DIMENSIONS     *d1, *d2;

/***********************************************/
/* Elastic Stiffness Matrix with two node at   */
/* each element                                */
/***********************************************/
    i = no_dof + 1;
    j = no_dof + 2; 
    k = no_dof + 3;

#ifdef DEBUG
 printf(" in beamst(): EA = %lf length = %lf \n", EA, length);
 printf(" in beamst(): EI = %lf cs= %lf sn = %lf \n", EI, cs, sn);
#endif


    t = EA/length; 
 
    s->uMatrix.daa[0][0]     =  t;
    s->uMatrix.daa[i-1][i-1] =  t;
    s->uMatrix.daa[0][i-1]   = -t;
    s->uMatrix.daa[i-1][0]   = -t;

    t = 12 * EI /(length*length*length);

    s->uMatrix.daa[1][1]     =  t;
    s->uMatrix.daa[j-1][j-1] =  t;
    s->uMatrix.daa[1][j-1]   = -t;
    s->uMatrix.daa[j-1][1]   = -t;

    t = (EI+ EI) / length;

    s->uMatrix.daa[2][2]   = s->uMatrix.daa[k-1][k-1] = t+t;
    s->uMatrix.daa[2][k-1] = s->uMatrix.daa[k-1][2]   = t;

    t = 6 * EI/length/length;
    s->uMatrix.daa[1][2]   = s->uMatrix.daa[2][1]       = t;
    s->uMatrix.daa[1][k-1] = s->uMatrix.daa[k-1][1]     = t;
    s->uMatrix.daa[2][j-1] = s->uMatrix.daa[j-1][2]     = -t;
    s->uMatrix.daa[j-1][k-1] = s->uMatrix.daa[k-1][j-1] = -t;

    /* ==================== */
    /* units buffer         */
    /* ==================== */

    if( CheckUnits() == ON ) {
       if(UNITS_TYPE == SI) {
          d1 = DefaultUnits("Pa");
          d2 = DefaultUnits("m");
       }
       else {
          d1 = DefaultUnits("psi");
          d2 = DefaultUnits("in");
       }

       UnitsMultRep( &(s->spColUnits[0]), d1, d2 );
       UnitsCopy( &(s->spColUnits[1]), &(s->spColUnits[0]) );
       UnitsMultRep( &(s->spColUnits[2]), &(s->spColUnits[0]), d2 );

       ZeroUnits( &(s->spRowUnits[0]) );
       ZeroUnits( &(s->spRowUnits[1]) );
       UnitsCopy( &(s->spRowUnits[2]), d2 ); 

       UnitsCopy( &(s->spColUnits[3]), &(s->spColUnits[0]) ); 
       UnitsCopy( &(s->spColUnits[4]), &(s->spColUnits[1]) ); 
       UnitsCopy( &(s->spColUnits[5]), &(s->spColUnits[2]) ); 

       UnitsCopy( &(s->spRowUnits[3]), &(s->spRowUnits[0]) ); 
       UnitsCopy( &(s->spRowUnits[4]), &(s->spRowUnits[1]) ); 
       UnitsCopy( &(s->spRowUnits[5]), &(s->spRowUnits[2]) ); 

       free((char *) d1->units_name);
       free((char *) d1);
       free((char *) d2->units_name);
       free((char *) d2);
    }

    s->uMatrix.daa = (double **) rotate(s->uMatrix.daa, cs, sn, size_of_stiff, no_dof);

#ifdef DEBUG
  printf("flag: in beamst() : STIFF after rotation \n");
  MatrixPrintIndirectDouble( s );
#endif
    return(s);

}

/* =================== */
/* Beam Mass Matrix    */
/* =================== */

#ifdef __STDC__
MATRIX *beamms(ARRAY *p, MATRIX *s, int mtype, double mass, double xl, double cs, double sn, int nst, int ndf)
#else
MATRIX *beamms(p, s, mtype, mass, xl, cs, sn, nst, ndf)
ARRAY  *p;
MATRIX *s;
double mass, xl, cs, sn;  /* mass, length, cosine, ans sine */
int nst, ndf,   mtype;
#endif
{
int    i, j, k, n, i1, j1;
int      length1, length2;
double  t, s1, s2, s3, dv;
double sg[5], ag[5], ba[3], rba[3], bb[5], rbb[5];
DIMENSIONS  *d1, *d2, *d3;

#ifdef DEBUG
       printf("*** Enter beamms() : matrix type = %10d  \n", mtype);
       printf("                   : t           = %10.3e\n", (mass * xl/2));
#endif

    t = mass * xl/2;
    switch(mtype) {
	case LUMPED:
             s->uMatrix.daa[0][0] = s->uMatrix.daa[1][1]  =  t;
             s->uMatrix.daa[2][2] =  t*xl*xl/12.0;   /* Using 16.0 in original version */

             s->uMatrix.daa[ndf][ndf] = t;
             s->uMatrix.daa[ndf+1][ndf+1] = t;
             s->uMatrix.daa[ndf+2][ndf+2] = t*xl*xl/12.0;   /* Using 16.0 in original version */

	     break;
	case CONSISTENT:

             gauss(sg, ag, 4);

             for(n =1; n<=4; n++) {

                 dv = 0.5 * xl * mass * ag[n];

                 s1 = (1 + sg[n])/2;
                 s2 = s1 * s1;
                 s3 = s1 * s2;

                 ba[1] = 1 - s1;
                 ba[2] = s1;
                 bb[1] = 1 - 3*s2 + s3 + s3;
                 bb[2] = xl* (s1- s2 - s2 + s3);
                 bb[3] = 3*s2 - s3 - s3;
                 bb[4] = xl  *  (-s2 + s3 );

                 rba[1] = ba[1] * dv;
                 rba[2] = ba[2] * dv;
                 for(i=1; i<= 4; i++)
                     rbb[i] = bb[i] * dv;

                 i = ndf;
                 s->uMatrix.daa[0][0] = s->uMatrix.daa[0][0] + ba[1] * rba[1];
                 s->uMatrix.daa[0][i] = s->uMatrix.daa[0][i] + ba[1] * rba[2];
                 s->uMatrix.daa[i][i] = s->uMatrix.daa[i][i] + ba[2] * rba[2];

                 for(i=1; i<= 4; i++){
                     i1 = i + 1;
                     if(i1 > 3 ) i1 = i1 + 1;
                     for(j=1; j<= 4; j++){
                         j1 = j  + 1;
                         if(j1 > 3) j1 = j1 +1;
                         s->uMatrix.daa[i1-1][j1-1] = s->uMatrix.daa[i1-1][j1-1] + bb[i]* rbb[j];
                     }
                 }
             }

             for(i = 0; i < 4; i++)
                 for(j = 0; j < i; j++)
                     s->uMatrix.daa[i][j] = s->uMatrix.daa[j][i];

	     break;
	default:
             FatalError("In elmt_frame2d() : beamms() : Type of Mass Matrix Undefined",(char*)NULL);
	     break;
    }

    /* ==================== */
    /* Initial units buffer */
    /* ==================== */

    if( CheckUnits() == ON ) {
       if(UNITS_TYPE == SI) {
          d1 = DefaultUnits("Pa");
          d2 = DefaultUnits("m");
       }
       else {
          d1 = DefaultUnits("psi");
          d2 = DefaultUnits("in");
       }
       d3 = DefaultUnits("sec");

       UnitsMultRep( &(s->spColUnits[0]), d1, d2 );
       UnitsCopy( &(s->spColUnits[1]), &(s->spColUnits[0]) );
       UnitsMultRep( &(s->spColUnits[2]), &(s->spColUnits[0]), d2 );

       UnitsPowerRep( &(s->spRowUnits[0]), d3, 2.0, NO );
       UnitsCopy( &(s->spRowUnits[1]), &(s->spRowUnits[0]) );
       UnitsMultRep( &(s->spRowUnits[2]), d2, &(s->spRowUnits[0]) );

       UnitsCopy( &(s->spColUnits[3]), &(s->spColUnits[0]) );
       UnitsCopy( &(s->spColUnits[4]), &(s->spColUnits[1]) );
       UnitsCopy( &(s->spColUnits[5]), &(s->spColUnits[2]) );

       UnitsCopy( &(s->spRowUnits[3]), &(s->spRowUnits[0]) );
       UnitsCopy( &(s->spRowUnits[4]), &(s->spRowUnits[1]) );
       UnitsCopy( &(s->spRowUnits[5]), &(s->spRowUnits[2]) );

       free((char *) d1->units_name);
       free((char *) d1);
       free((char *) d2->units_name);
       free((char *) d2);
       free((char *) d3->units_name);
       free((char *) d3);
    }

    /* ==================== */
    /* Rotate mass matrix   */
    /* ==================== */

    s->uMatrix.daa = (double **) rotate(s->uMatrix.daa,cs,sn,nst,ndf);

#ifdef DEBUG
       MatrixPrint(s);
       printf("*** leaving beamms()  \n");
#endif

    return(s);
}

/* ====================== */
/* Rotate                 */
/* ====================== */

#ifdef __STDC__
double **rotate(double **s, double cs, double sn, int nst, int ndf)
#else
double **rotate(s, cs, sn, nst, ndf)
double **s;
int  nst, ndf;
double cs, sn;
#endif
{
int    i,j,n;
double t;

   if(cs == 1.0)
      return(s);

   for(i = 0; i < nst; i = i + ndf){
       j = i+1;
       for(n = 0; n < nst;n++){
           t = s[n][i] * cs - s[n][j] * sn;
           s[n][j] = s[n][i] * sn + s[n][j] *cs;
           s[n][i] = t;
       }
   }

   for(i = 0; i < nst; i=i+ndf){
       j = i+1;
       for(n = 0; n < nst; n++){
           t = cs*s[i][n] - sn*s[j][n];
           s[j][n] = sn * s[i][n] + cs * s[j][n];
           s[i][n] = t;
       }
   }
 
   return(s);
}


/* =================================================== */
/* Equivalent Loading Procedure for 2 D frame element  */
/* =================================================== */

#ifdef __STDC__
ARRAY *sld07(ARRAY *p, int task)
#else
ARRAY *sld07(p,task)
ARRAY *p;
int task;
#endif
{
ELEMENT_LOADS *elsptr;
ELOAD_LIB        *elp;
double  P, a ,b ;
double  L, load[8];
double  px,py,pz, mx,my,mz, bx,by,bz, ze,ze2,ze3;
double  X1,Y1,MZ1,
        X2,Y2,MZ2, 
        MCZ, MCZT;
double  f1,f2,f3,f4,f5,f6, Z1, Z2,
        df1x, df3x, df5x,df2x, df6x;
double  **rmat,**rt;
int     *da;
int     inc,j, i;

    /* Initialize total load */

    for(inc=1; inc<=7;inc++)
        load[inc-1] = 0.0;
    MCZT = 0.0;

    switch(task){
        case PRESSLD:
        case STRESS:
             L       = (double) p->length.value;  
             elsptr  =  p->elmt_load_ptr;
             for (j=1; j<= elsptr->no_loads_faces;j++) {
                  elp = &elsptr->elib_ptr[j-1];

             /* Pt loads */

             P = elp->P.value;
	     a = elp->a.value;
	     b = elp->b.value;
             px = elp->px.value;
	     py = elp->py.value;
       
             /* moments */

             mz = elp->mz.value;

             /* distributed loading */

             bx =  elp->bx.value;
             by =  elp->by.value;

             if(a > L)
                printf(">>ERROR in sld; Elmt Load dist. 'a' > Elmt Length; El_no= %d\n",p->elmt_no);  

             if (elp->type == -1) { /* Distributed loading Condition */
                inc = 0;
                /* set default values */
                if(b == 0.0) b = L; /* dist loading acts on entire length */

                /* first calc f(b) */

                ze = b/L;    
SHP_START:
                ze2 = ze * ze; ze3 = ze2 * ze;
                f1 = 1   -  ze2/2;
                f2 = ze3*ze/2  - ze3 + ze ;
                f3 = (ze3*ze/4 - 2*ze3/3 + ze2/2) * L;
                f4 = ze2/2 ;
                f5 = -ze *ze3/2 +  ze3;
                f6 = (ze3*ze/4  - ze3/3) * L;
                inc++;
                 
                if( inc == 1) {
                   /* temp hold of values f(b) */
                   X1 = f1; Y1 = f2;  Z1 = f3; 
                   X2 = f4; Y2 = f5;  Z2 = f6; 
                   ze = a/L;
                   goto SHP_START;
                }
                else{
                   /* f() = f(b) - f(a)  */
                   f1 = X1 - f1; f2 =  Y1 - f2;f3 = Z1 - f3; 
                   f4 = X2 - f4; f5 =  Y2 - f5;f6 = Z2 - f6; 
                }

                X1 = bx * f1 * L;
                Y1 = by * f2  * L;
                MZ1 = by * f3 * L;

                X2  = bx * f4 * L;
                Y2  = by * f5 * L;
                MZ2 = by * f6 * L;

                /* +ve  simply support moment at center */

                if (task == STRESS){

                   if(b==L && a== 0.0)  /* udl acting on entire length */
                      MCZ = -by * (L * L)/8;   
                   else   /* approximate mom at center */
                    MCZ = -by *(b-a)* (L - (a+b)/2)/2;   
                }

            } /* end of dist loading */

	    /* Concentrated Loading Condition */
            /* distributed body force         */

            else {
              /* shape functions */
              /* =========================================*/
              /* according to :                           */
              /* "Structural Dynamics by Finite Elements" */
              /* W. Weaver, Jr. and P. R. Johnson         */
              /* Chapter 6: Space Frame                   */
              /* =========================================*/


              ze = a/L;      ze2 = ze * ze;        ze3 = ze2 * ze;

              f1 =     1  -  ze;
              f2 =     2 * ze3  -  3 * ze2 + 1;
              f3 =    (ze3 - 2 * ze2 + ze) * L;
              f4 =     ze ;
              f5 =    -2 * ze3  +  3 * ze2;
              f6 =    (ze3 - ze2 ) * L;

              /* derivatives of shape function */  

              if(mz != 0.0){
                 df2x =  6 *( ze2  -  ze) / L;
                 df3x =  3 *ze2  - 4*ze  + 1;
                 df5x = - df2x; 
                 df6x =  3 *ze2  - 2*ze;
              }

              X1 = px * f1 + 0;
              Y1 = py * f2 + mz *  df2x;
              MZ1 = py * f3 + mz *  df3x;

              X2 = px * f4 + 0;
              Y2 = py * f5 + mz *  df5x;
              MZ2 = py * f6 + mz *  df6x;

              /* +ve  simply support moment at center */
              if(task == STRESS) 
                 MCZ = -py * (L -a)/2 ;   
                  
            }

            /* Add Contributation to Total Equivalent Load */

               load[0] = load[0] + X1;
               load[1] = load[1] + Y1;
               load[2] = load[2] + MZ1;
               load[3] = load[3] + X2;
               load[4] = load[4] + Y2;
               load[5] = load[5] + MZ2;
	       
               if(task == STRESS)
	          MCZT = MCZT+ MCZ;
            }
            break;
       case STRESS_LOAD:
            for(i = 1; i<= 6; i++)
            load[i-1] = p->nodal_loads[i-1].value ;
            break;
       default:
            break;
    }

    /* Rotate Local forces to Global forces */
 
    rmat = (double **) MatrixAllocIndirectDouble(p->size_of_stiff, p->size_of_stiff);
    rmat = (double **) tmat(rmat,4,p);
    rt   = (double **) MatrixAllocIndirectDouble(p->size_of_stiff, p->size_of_stiff);
    for( i=1 ; i<=p->size_of_stiff ; i++ )
       for( j=1 ; j<=p->size_of_stiff ; j++ )
          rt[j-1][i-1] = rmat[i-1][j-1];

    for (inc=1; inc<=p->size_of_stiff; inc++){
         p->nodal_loads[inc-1].value = 0.0;
         for (j=1; j<=p->size_of_stiff; j++) {
              p->nodal_loads[inc-1].value += rt[inc-1][j-1]* (double) load[j-1];
         }
    }
    /* mid pt moment */
/*
    if (task == STRESS) {
        p->nodal_loads[7].value = MCZT; 
    }
*/
 
    MatrixFreeIndirectDouble(rmat, p->size_of_stiff);
    MatrixFreeIndirectDouble(rt, p->size_of_stiff);

    return(p);
}

