/*
 *  ============================================================================= 
 *  ALADDIN Version 1.0 :
 *    elmt_shell_4n_q.c : 4 nodes shell quadratic element
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
 *   Element SHELL_FOUR_NODES
 *        Shell Element:
 *        material properties array
 *        Input Properties:
 *
 *    p->work_material[0] = E;
 *    p->work_material[1] = G;
 *    p->work_material[2] = fy;
 *    p->work_material[3] = ET;
 *    p->work_material[4] = nu;
 *    p->work_material[5] = density;
 *    p->work_material[6] = fu;
 *
 *    p->work_section[0] = Ixx;
 *    p->work_section[1] = Iyy;
 *    p->work_section[2] = Izz;
 *    p->work_section[3] = Ixy;
 *    p->work_section[4] = Ixz;
 *    p->work_section[5] = Iyz;
 *    p->work_section[6] = weight;
 *    p->work_section[7] = bf;
 *    p->work_section[8] = tf;
 *    p->work_section[9] = depth;                                  
 *    p->work_section[10] = area;                                 
 *    p->work_section[11] = thickness; 
 *
 *    p->no_dimen          = ndm;    ( = 3 )
 *    p->elmt_type         = &iel;   ( = 16 )
 *    p->nodes_per_elmt    = nen;    ( = 4 )
 *    p->nodes_per_elmt    = nel;    ( = 4 )    = nen ???
 *    p->dof_per_node      = ndf;    ( = 6 )
 *    p->size_of_stiff     = nst;
 *    p->type              = d[16];   ???
 *    p->coord[ndm][nen]   = xl(ndm,*);
 *    p->node_connect[nen] = ix(nen);

 *    p->LC_ptr->ialph     = ialph;     Rotation constant
 *    p->LC_ptr->pen       = pen;       Penalty constant
 *    p->LC_ptr->load[0-5] = d[6-11];   Loading parameters
 *
 *  Written by: Lanheng Jin                                         December 1995
 *  ============================================================================= 
 */


#include <stdio.h>
#include <string.h>
#include <math.h>

#include "defs.h"
#include "units.h"
#include "matrix.h"
#include "vector.h"
#include "fe_database.h"
#include "symbol.h"
#include "fe_functions.h"
#include "elmt.h"
#include "miscellaneous.h"

/*
#define DEBUG 
*/

static double      b[4], c[4], aa[4], bb[4], cc[4], dd[4], ee[4];
static double      bm[3][6][9];
static int         ii1, ii2;
static double      shp[3][8][9], shp1[3][4][9], shp2[3][4][9];
static int         nlla, nllb, iel, ia[2][26], ir[2][26];
static int         false=0, true=1;
static int         mct = 0;

ARRAY *elmt_shell_4nodes_q(p, isw)
ARRAY *p;
int    isw;
{

static double                   nu, h; 
static QUANTITY                 G, E, density; 

static double *sg = NULL;
static double *tg = NULL;
static double *wg = NULL;
static double **yl = NULL;
static double **tr = NULL;
static double *eps = NULL;
static double *sigi = NULL;
static double *sigb = NULL;
static double **gshp1 = NULL;
static double **gshp2 = NULL;
static double **btd   = NULL;
static double **vl    = NULL;
static double *dvl    = NULL;
static double **s     = NULL;
static double *pp     = NULL;
static double *gg     = NULL;
static double *bl     = NULL;
static double *d      = NULL;

double        pen;
int           ialph;
int           l, lint;
int           qdflg;
int           ndm, ndf, nst, nen, nel;
int           i1, i2, j1, j2, ii, jj, kk;
int           i, j, k;
int           length1, length2, length;
double        a11, a12, a13, a21, a22, a23, a31, a32, a33, xx, yy, zz;
double        dv, dv1, dv2, dv3, dv4, dv5, dv6, ggi, ggv, tc, xsj, xsjt;
double        xyn, x1n, x2n, y1n, y2n, xn, yn, sn;
double        shp1i, shp2i, shp3i, shp11, shp12, shp13, shp21, shp22, shp23;
double        **ul;
DIMENSIONS    *dimen1, *dimen2, *dimen3;
double        coord, stress, moment, curvtr;

FILE *fopen(), *fp;

#ifdef DEBUG
    printf("*** Enter elmt_shell_4nodes_q() : isw = %4d\n", isw);
#endif

    ndm = p->no_dimen;        /*    = 3      */ 
    nel = p->nodes_per_elmt;  /*    = 4      */
    nen = p->nodes_per_elmt;  /*    = 4      */
    ndf = p->dof_per_node;    /*    = 6      */
    nst = p->size_of_stiff;   /*    = nst    */

    if( isw != PROPTY && isw != CHERROR ) {

       l = d[14];
       if( l*l != lint )
          pgauss( l, &lint, sg, tg, wg );
       
/* -- Calculation of transformation array and surface coordiantes -- */

       if(yl == NULL)  yl = MatrixAllocIndirectDouble(3,4);
       if(tr == NULL)  tr = MatrixAllocIndirectDouble(3,3);

       tran06(p, yl, tr);

/*-----Test for triangular element-----*/
 
       if( p->node_connect[0] == p->node_connect[1] ||
           p->node_connect[1] == p->node_connect[2] ||
           p->node_connect[2] == p->node_connect[3] ||
           p->node_connect[3] == p->node_connect[0] ) {

          qdflg = false;
          for( i=1; i<=3; i++ ) {
             for( j=1; j<=6; j++ ) {
                for( k=1; k<=9; k++ ) {
                   bm[i-1][j-1][k-1] = 0.0;
                }
             }
          }

          shp11 = 0.0;
          shp12 = 0.0;
          shp13 = 0.0;
          shp21 = 0.0;
          shp22 = 0.0;
          shp23 = 0.0;
          a13   = 0.0;
          a23   = 0.0;
          a31   = 0.0;
          a32   = 0.0;
          a33   = 0.0;
          x1n   = 0.0;
          y1n   = 0.0;
          x2n   = 0.0;
          y2n   = 0.0;
          xyn   = 0.0;

          jtri06( yl, &xsjt );
       }
       else {
          qdflg = true;
          jacq06( yl );
       }

       if( qdflg && yl[2][0] == 0.0 ) {
          ii1 = 3;
          ii2 = 5;
       }
       else {
          ii1 = 1;
          ii2 = 6;
       }

/* Construct the integrals of the drilling shape functions */

       for( i=1; i<=3; i++ ) {
          for( j=1; j<=4; j++ ) {
             gshp1[i-1][j-1] = 0.0;
             gshp2[i-1][j-1] = 0.0;
          }
       }
       
       for( i=1; i<=24; i++ ) 
          gg[i-1] = 0.0;

       dv = 0.0;
     
       for( k=1; k<=lint; k++ ) {

          rshp06( k-1, sg[k-1], tg[k-1], yl, &xsj, ndm);

          dvl[k-1] = xsj * wg[k-1];
          dv      += dvl[k-1];
          for( i=1; i<=3; i++ ) {
             for( j=1; j<=4; j++ ) {
                gshp1[i-1][j-1]  += shp1[i-1][j-1][k-1]*dvl[k-1];
                gshp2[i-1][j-1]  += shp2[i-1][j-1][k-1]*dvl[k-1];
             }
          }
       }

       ggv = pen * d[2]/dv;
    }

/**********************************************/
/**    Begining of the main switch           **/
/**********************************************/

    switch(isw) {

        case PROPTY: 

#ifdef DEBUG
    printf(" In elmt_shell_4nodes_q(): \n");
    printf("     : enter case of PROPTY\n");
#endif

        if(sg == NULL)     sg = dVectorAlloc( 9 );
        if(tg == NULL)     tg = dVectorAlloc( 9 );
        if(wg == NULL)     wg = dVectorAlloc( 9 );
        if(gg == NULL)     gg = dVectorAlloc( 24 );
        if(dvl == NULL)    dvl = dVectorAlloc( 9 );
        if(pp == NULL)     pp = dVectorAlloc(nst);
        if(gshp1 == NULL)  gshp1 = MatrixAllocIndirectDouble( 3, 4 );
        if(gshp2 == NULL)  gshp2 = MatrixAllocIndirectDouble( 3, 4 );
        if(btd == NULL)    btd = MatrixAllocIndirectDouble(3,6);
        if(s == NULL)      s = MatrixAllocIndirectDouble(nst, nst);

/* ---- material properties:  elastic ------------- */

        E.value       =  p->work_material[0].value;
        nu            =  p->work_material[4].value;
        density.value =  p->work_material[5].value;

/* ---- check possion_ratio value ------------------ */

        if( nu == 0.0 || nu > 0.5 ) {
           printf("WARNING >> ... In elmt_shell_4nodes_q() : nu value = %9.4f, reset to 0.3 !\n", nu);
           nu = 0.3;                /* default poi_ratio value */
        }

/* ---- calculate  G value ------------------------- */

        G.value = p->work_material[1].value = E.value/(1.0 - 2.0*nu) ;

	if( CheckUnits() == ON ) {
            E.dimen       =  p->work_material[0].dimen;
            density.dimen =  p->work_material[5].dimen;
            G.dimen = E.dimen;
	}

        if(E.value/((1.0 - 2.0*nu)) != p->work_material[1].value) {
           printf("WARNING >> ... In elmt_shell_4nodes_q(): G is not equal to E/(1-2nu), check G for homogeneous material \n");
           printf("               Ignore this message for non-homogeneous materials \n");
        }

        /* Get the thickness of the shell  */

        h = p->work_section[11].value;  

        if(d == NULL)  d = dVectorAlloc( 20 );

/*
        fp = fopen("data", "r");

        if(fp) fscanf(fp, "%d%lf%lf%lf%lf%lf%lf%lf%lf", 
                      &ialph, &pen, d+6, d+7, d+8, d+9, d+10, d+11, d+16);
        else {
           printf("\nError for opening a file\n");
           exit(1);
        }

        fclose (fp);
*/

	ialph = p->LC_ptr->ialph;
	pen   = p->LC_ptr->pen;
	d[6]  = p->LC_ptr->load[0];
	d[7]  = p->LC_ptr->load[1];
	d[8]  = p->LC_ptr->load[2];
	d[9]  = p->LC_ptr->load[3];
	d[10] = p->LC_ptr->load[4];
	d[11] = p->LC_ptr->load[5];
        d[0]  = E.value/(1.0-nu*nu)*h;
        d[1]  = nu*d[0];
        d[2]  = (d[0] - d[1])/2.0;

        xx    = h*h/12.0;

        d[3]  = d[0]*xx;
        d[4]  = d[1]*xx;
        d[5]  = d[2]*xx;
        d[12] = density.value*h;
        d[13] = d[12]*xx;

/* ---- Quadrature order --------------------------- */

        d[14] = 2.0;
        d[15] = 1.0;
        lint  = 0;

/* ---- Construct rotation parameters : u-x = 1, u-y = 2 */

        iel = 16;       /* set up the element type 16   */

        ia[0][iel-1] = 1;
        ia[1][iel-1] = 2;

/* ---- Construct rotation parameters : theta-x = 4, theta-y = 5 */

        ir[0][iel-1] = 4;
        ir[1][iel-1] = 5;

/* ---- Print out some parameters ---------------------- */
/*
        printf("Output of material properties, section parameters, load and some constants\n");
        printf("    : Young's Modulus =   %l3.5e\n", E.value);
        printf("    : Poisson Ratio   = %13.5e\n", nu);
        printf("    : Density         = %13.5e\n", density.value);
        printf("    : Shell thickness = %13.5e\n", h);
        printf("    : 1-gravity       = %13.5e\n", d[6]);
        printf("    : 2-gravity       = %13.5e\n", d[7]);
        printf("    : pressure        = %13.5e\n", d[8]);
        printf("    : x-gravity       = %13.5e\n", d[9]);
        printf("    : y-gravity       = %13.5e\n", d[10]);
        printf("    : z-gravity       = %13.5e\n", d[11]);
        printf("    : Load Constant   = %13.5e   (0=lump, 1=cons.)\n", d[16]);
        printf("    : Rotation Const. = %3d             (0=none, 1=hughes)\n", ialph);
        printf("    : Penalty Const.  = %13.5e\n\n", pen);
*/
#ifdef DEBUG
    printf(" In elmt_shell_4nodes_q(): \n");
    printf("    : leaving case of PROPTY\n");
#endif
        break;

        case CHERROR:
        break;

        case STIFF:              /* form element stiffness */
                                 /* Compute the element tangent array  */
        case LOAD_MATRIX:        /* form load matrix */

#ifdef DEBUG
    printf("*** elmt_shell_4nodes_q() : In case STIFF\n");
#endif

/* ---- Construct the modified drilling shape functions -- */

        for( i=1; i<=3; i++ ) {
           for( j=1; j<=4; j++ ) {
              dv1 = gshp1[i-1][j-1]/dv;
              dv2 = gshp2[i-1][j-1]/dv;
              for( k=1; k<=lint; k++ ) {
                 shp1[i-1][j-1][k-1] -= dv1;
                 shp2[i-1][j-1][k-1] -= dv2;
              }
           }
        }

/* ---- Transform the element load vector ---------- */

        if(bl == NULL)  bl = dVectorAlloc(3);

        for( i=1; i<=3; i++ ) {
           bl[i-1] = d[i+5];                        /*  set dm=1 here  */
           for( j=1; j<=3; j++ )
              bl[i-1]  += tr[i-1][j-1]*d[i+8];      /*  set dm=1 here  */
        }

        tc = d[17-1];

/* ---- Compute membrane/load parts ---------------- */
 
        for( k=1; k<=lint; k++ ) {

        /* Membrane and bending stiffness parts */

           if( qdflg ) 
              dktq06( k-1 );
           
           else 
              dktb06( sg[k-1], tg[k-1], xsjt );
         
           dv  = dvl[k-1];
           dv1 = d[1-1] * dv;
           dv2 = d[2-1] * dv;
           dv3 = d[3-1] * dv;
           dv4 = d[4-1] * dv;
           dv5 = d[5-1] * dv;
           dv6 = d[6-1] * dv;
 
           i1 = 1;
           for( i=1; i<=4; i++ ) {

           /* Recover the previously computed shape functions */

              shp1i = shp[0][i-1][k-1];
              shp2i = shp[1][i-1][k-1];
              shp3i = shp[2][i-1][k-1];

              if( qdflg ) {
                shp11 = shp1[0][i-1][k-1];
                shp12 = shp1[1][i-1][k-1];
                shp13 = shp1[2][i-1][k-1];
                shp21 = shp2[0][i-1][k-1];
                shp22 = shp2[1][i-1][k-1];
                shp23 = shp2[2][i-1][k-1];
                a13   = -dv1*shp11 - dv2*shp22;
                a23   = -dv2*shp11 - dv1*shp22;
                a31   =  dv3*shp2i;
                a32   =  dv3*shp1i;
                a33   = -dv3*(shp12+shp21);
              }

              i2 = i1-1;

              /* Compute the loading term */

              pp[i1-1] +=              shp3i*bl[0] *dv;
              pp[i1  ] +=              shp3i*bl[1] *dv;
              pp[i1+1] +=              shp3i*bl[2] *dv;
              pp[i1+2] +=              shp13*bl[2] *dv*tc;
              pp[i1+3] +=              shp23*bl[2] *dv*tc;
              pp[i1+4] -= (shp13*bl[0]+shp23*bl[1])*dv*tc;

              /* Form the stress-displacement matrix (Bi-trans 8 D) */

              a11 = dv1*shp1i;
              a12 = dv2*shp2i;
              a21 = dv2*shp1i;
              a22 = dv1*shp2i;

              /* Form the plate stress-displacement matrix */

              for( ii=ii1; ii<=ii2; ii++ ) {
                 btd[0][ii-1] = dv4*bm[0][ii-1][i-1] + dv5*bm[1][ii-1][i-1];
                 btd[1][ii-1] = dv5*bm[0][ii-1][i-1] + dv4*bm[1][ii-1][i-1];
                 btd[2][ii-1] = dv6*bm[2][ii-1][i-1];
              }

              /* Loop on columns */

              j1 = i1;
              for( j=i; j<=4; j++ ) {

                 j2 = j1 - 1;
                 xn = shp[0][j-1][k-1];
                 yn = shp[1][j-1][k-1];
                 sn = shp[2][j-1][k-1];

                 if( qdflg ) {
                   x1n = - shp1[0][j-1][k-1];
                   y2n = - shp2[1][j-1][k-1];
                   x2n = - shp1[1][j-1][k-1];
                   y1n = - shp2[0][j-1][k-1];
                   xyn =   x2n + y1n;
                 }

                 /* Compute the membrane part */
 
                 s[i1-1][j1-1] += a11*xn + a31*yn;
                 s[i1  ][j1-1] += a12*xn + a32*yn;
                 s[i1-1][j1  ] += a21*yn + a31*xn;
                 s[i1  ][j1  ] += a22*yn + a32*xn;
                 s[i1+4][j1-1] += a13*xn + a33*yn;
                 s[i1+4][j1  ] += a23*yn + a33*xn;
                 s[i1-1][j1+4] += a11*x1n + a21*y2n + a31*xyn;
                 s[i1  ][j1+4] += a12*x1n + a22*y2n + a32*xyn;
                 s[i1+4][j1+4] += a13*x1n + a23*y2n + a33*xyn;

                 /* Compute the bending part */

                 for( ii=ii1; ii<=ii2; ii++ ) {
                    for( jj=ii1; jj<=ii2; jj++ ) {
                       for( kk=1; kk<=3; kk++ ) 
                          s[i2+ii-1][j2+jj-1] += btd[kk-1][ii-1]*bm[kk-1][jj-1][j-1];
                    }
                 }
                 j1 = j1 + ndf;
              }
              i1 = i1 + ndf;
           }

           /* Compute the hughes/brezzi rotation matrix */

           if( ialph != 0 && qdflg ) {

              j1 = 0;

              for( j=1; j<=4; j++ ) {
                 gg[j1  ] -= shp[1][j-1][k-1]*dv;
                 gg[j1+1] += shp[0][j-1][k-1]*dv;
                 gg[j1+5] -= 2.0*shp[2][j-1][k-1]*dv;

                 if( ialph > 1 ) 
                    gg[j1+5] += (shp1[1][j-1][k-1] - shp2[0][j-1][k-1])*dv;
                 
                 j1 = j1 + ndf;
              }
           }
        }

/* -- Perform the rank one update for the hughes/brezzi term -- */

        if( ialph > 0 && qdflg ) {

           for( i=1; i<=nst; i++ ) {

              ggi = ggv * gg[i-1];
              for( j=1; j<=nst; j++ ) 
                 s[i-1][j-1] += ggi * gg[j-1];
              
           }
        }

        i1 = 1;
        if( yl[2][0] != 0.0 ) {
           for( i=1; i<=4; i++ ) {
              pp[i1+2] += yl[2][i-1]*pp[i1  ];
              pp[i1+3] -= yl[2][i-1]*pp[i1-1];
              j1 = i1;
              for( j=i; j<=4; j++ ) {
                 proj06( i1-1, j1-1, s, yl[2][i-1], yl[2][j-1], nst);
                 j1 = j1 + ndf;
              }
              i1 = i1 + ndf;
           }
        }

        rots06( s, pp, tr, nst, ndf );

/*
#ifdef DEBUG
    dMatrixPrint("s -- Stiffness Matrix", s, nst, nst);
    printf("pp vector\n");
    for(i=1; i<=nst; i++)
       printf("%5d%15.8e\n", i, pp[i-1]);
    dMatrixPrint("p->displ->uMatrix.daa", p->displ->uMatrix.daa, p->displ->iNoRows, p->displ->iNoColumns);
#endif
*/

/* ---- Form residual ------------------------ */

        ul = p->displ->uMatrix.daa;

        for( i=1; i<=nst; i++ ) {
           for( j=1; j<=nst; j++ ) 
              pp[i-1] -= s[i-1][j-1]*ul[(j-1)/4][(j-1)%4];
        }
/*
    printf("pp vector\n");
    for(i=1; i<=nst; i++)
    printf("%5d%15.8e\n", i, pp[i-1]);
*/
            

/* Transfer stiff matrix and equiv. nodal load vector to ARRAY *p */

    switch(isw) {
      case STIFF:
        p->stiff->uMatrix.daa = s;

    /**************************************************/
    /* Assign Units to Stiffness Matrix               */
    /**************************************************/

       if( CheckUnits() == ON ) {
          if(UNITS_TYPE == SI || UNITS_TYPE == SI_US ) {
             dimen1 = DefaultUnits("Pa");
             dimen2 = DefaultUnits("m");
          }
          else {
             dimen1 = DefaultUnits("psi");
             dimen2 = DefaultUnits("in");
          }

          /* node 1 */
          UnitsMultRep( &(p->stiff->spColUnits[0]), dimen1, dimen2 );
          UnitsCopy( &(p->stiff->spColUnits[1]), &(p->stiff->spColUnits[0]) );
          UnitsCopy( &(p->stiff->spColUnits[2]), &(p->stiff->spColUnits[0]) );
          UnitsMultRep( &(p->stiff->spColUnits[3]), &(p->stiff->spColUnits[0]), dimen2 );
          UnitsCopy( &(p->stiff->spColUnits[4]), &(p->stiff->spColUnits[3]) );
          UnitsCopy( &(p->stiff->spColUnits[5]), &(p->stiff->spColUnits[3]) );

          ZeroUnits( &(p->stiff->spRowUnits[0]) );
          ZeroUnits( &(p->stiff->spRowUnits[1]) );
          ZeroUnits( &(p->stiff->spRowUnits[2]) );
          UnitsCopy( &(p->stiff->spRowUnits[3]), dimen2 );
          UnitsCopy( &(p->stiff->spRowUnits[4]), dimen2 );
          UnitsCopy( &(p->stiff->spRowUnits[5]), dimen2 );
     
          /* node i  i > 1*/
        
          for ( i = 2; i <= p->nodes_per_elmt; i++) {
                kk = p->dof_per_node*(i-1) + 3; 
                for( j = 1; j <= p->dof_per_node; j++) {
                     k  = p->dof_per_node*(i-1) + j;
                     if( k <= kk) {
                        UnitsCopy( &(p->stiff->spColUnits[k-1]), &(p->stiff->spColUnits[0]) );
                        UnitsCopy( &(p->stiff->spRowUnits[k-1]), &(p->stiff->spRowUnits[0]) );
                     }
                     if(k > kk) {
                        UnitsCopy( &(p->stiff->spColUnits[k-1]), &(p->stiff->spColUnits[3]) );
                        UnitsCopy( &(p->stiff->spRowUnits[k-1]), &(p->stiff->spRowUnits[3]) );
                     }
                 }
           }
           free((char *) dimen1->units_name);
           free((char *) dimen1);
           free((char *) dimen2->units_name);
           free((char *) dimen2);
      }

      break;
      case LOAD_MATRIX:
         for( i=1; i<=nst; i++ ) 
           p->equiv_nodal_load->uMatrix.daa[i-1][0] = pp[i-1];
      break;
      default:
      break;
    }

/*
#ifdef DEBUG
    dMatrixPrint("p->equiv_nodal_load->uMatrix.daa", p->equiv_nodal_load->uMatrix.daa, 
                   p->equiv_nodal_load->iNoRows, p->equiv_nodal_load->iNoColumns);
#endif
*/

#ifdef DEBUG
    printf(" In elmt_shell_4nodes_q(): \n");
    printf("    : leaving case of STIFF\n");
#endif
        break;

        case STRESS:     /* Compute and output the element variables */

#ifdef DEBUG
        printf("*** elmt_shell_4nodes_q() : In case STRESS\n");
#endif

        if(vl == NULL)    vl = MatrixAllocIndirectDouble( 6, 4 );
        if(eps == NULL)   eps  = dVectorAlloc(6);
        if(sigi == NULL)  sigi = dVectorAlloc(6);
        if(sigb == NULL)  sigb = dVectorAlloc(6);

        ul = p->displ->uMatrix.daa;
        for( i=1; i<=4; i++ ) {
           for( j=1; j<=3; j++ ) {
              vl[j-1][i-1] = 0.0; 
              vl[j+2][i-1] = 0.0; 
              for( k=1; k<=3; k++ ) {
                 vl[j-1][i-1] += tr[j-1][k-1]*ul[k-1][i-1]; 
                 vl[j+2][i-1] += tr[j-1][k-1]*ul[k+2][i-1]; 
              }
           }
        }

        l = d[15];

        if( l*l != lint ) 
           pgauss ( l, &lint, sg, tg, wg );

        for( k=1; k<=lint; k++ ) {

           rshp06( k-1, sg[k-1], tg[k-1], yl, &xsj, ndm);

           /* modify the rotational shape functions */

           for( i=1; i<=3; i++ ) {
              for( j=1; j<=4; j++ ) {
                 shp1[i-1][j-1][0] -= gshp1[i-1][j-1]/dv;
                 shp2[i-1][j-1][0] -= gshp2[i-1][j-1]/dv;
              }
           }

           stre06(p,d,vl, ndm, nel,1, &xx, &yy, &zz, eps, sigi, sigb);

	   if( UNITS_TYPE == SI || UNITS_TYPE == SI_US ) {
	      dimen1 = DefaultUnits("N");
	      dimen2 = DefaultUnits("m");
	   }
	   else {
	      dimen1 = DefaultUnits("lbf");
	      dimen2 = DefaultUnits("in");
	   }

	   coord  = dimen2->scale_factor;
	   stress = dimen1->scale_factor/dimen2->scale_factor;
	   moment = dimen1->scale_factor;
	   curvtr = 1.0/dimen2->scale_factor;

           free((char *) dimen1->units_name);
           free((char *) dimen1);
           free((char *) dimen2->units_name);
           free((char *) dimen2);

           printf("                       S h e l l   S t r e s s e s");
           printf("\n elmt,x-coord  xx-stress  xy-stress  yy-stress   1-stress   2-stress   angle");
           printf("\n matl,y-coord  xx-moment  xy-moment  yy-moment   1-moment   2-moment   angle");
           printf("\n      z-coord  xx-strain  xy-strain  yy-strain  xx-curvtr  xy-curvtr  yy-curvtr");
           if( UNITS_TYPE == SI )
              printf("\n   coord = m, stress = N/m, moment = N*m/m, curvtr = 1/m, angle = rad");
           else
              printf("\n  coord = in, stress = lbf/in, moment = lbf*in/in, curvtr = 1/in, angle = rad");
           printf("\n - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -");
           printf("\n     %8.3lf%11.3e%11.3e%11.3e%11.3e%11.3e%8.2lf",
                   xx/coord,sigi[0]/stress,sigi[1]/stress,sigi[2]/stress,sigi[3]/stress,sigi[4]/stress,sigi[5]);
           printf("\n     %8.3lf%11.3e%11.3e%11.3e%11.3e%11.3e%8.2lf",
                   yy/coord,sigb[0]/moment,sigb[1]/moment,sigb[2]/moment,sigb[3]/moment,sigb[4]/moment,sigb[5]);
           printf("\n     %8.3lf%11.3e%11.3e%11.3e%11.3e%11.3e%11.3e\n\n",
                   zz/coord,eps[0],eps[1],eps[2],eps[3]/curvtr,eps[4]/curvtr,eps[5]/curvtr);
        }
        
#ifdef DEBUG
    printf(" In elmt_shell_4nodes_q(): \n");
    printf("    : leaving case of STRESS\n");
#endif
        break;

        case MASS_MATRIX:                   /* Compute element mass array */

#ifdef DEBUG
    printf(" In elmt_shell_4nodes_q(): \n");
    printf("    : entering case of MASS_MATRIX\n");
#endif
 
        for(i=1; i<=nst; i++)
           pp[i-1] = 0.0;

        for( k=1; k<=lint; k++ ) {
           dv1 = dvl[k-1]*d[12];
           i1 = 0;
           for( j=1; j<=4; j++ ) {
              for( i=1; i<=3; i++ ) {
                 pp[i1+i-1] += shp[2][j-1][k-1]*dv1;
                 s[i1+i-1][i1+i-1] = pp[i1+i-1];
              }
              i1 = i1 + ndf;
           }
        }
	p->stiff->uMatrix.daa = s;

 /* ------------ MASS UNITS ---------------------------- */
 /* The units type is determined by the SetUnitsType()   */
 /* ---------------------------------------------------- */

 /* Initiation of Mass Units Buffer                      */

      switch( CheckUnits() ) {
        case ON:
           if(UNITS_TYPE == SI || UNITS_TYPE == SI_US ) {
              dimen1 = DefaultUnits("Pa");
              dimen2 = DefaultUnits("m");
           }
           else {
              dimen1 = DefaultUnits("psi");
              dimen2 = DefaultUnits("in");
           }
           dimen3 = DefaultUnits("sec");

          /* node no 1 */
           UnitsMultRep( &(p->stiff->spColUnits[0]), dimen1, dimen2 );
           UnitsCopy( &(p->stiff->spColUnits[1]), &(p->stiff->spColUnits[0]) );
           UnitsCopy( &(p->stiff->spColUnits[2]), &(p->stiff->spColUnits[0]) );
           UnitsMultRep( &(p->stiff->spColUnits[3]), &(p->stiff->spColUnits[0]), dimen2 );
           UnitsCopy( &(p->stiff->spColUnits[4]), &(p->stiff->spColUnits[3]) );
           UnitsCopy( &(p->stiff->spColUnits[5]), &(p->stiff->spColUnits[3]) );

           UnitsPowerRep( &(p->stiff->spRowUnits[0]), dimen3, 2.0, NO );
           UnitsCopy( &(p->stiff->spRowUnits[1]), &(p->stiff->spRowUnits[0]) );
           UnitsCopy( &(p->stiff->spRowUnits[2]), &(p->stiff->spRowUnits[0]) );
           UnitsMultRep( &(p->stiff->spRowUnits[3]), dimen2, &(p->stiff->spRowUnits[0]) );
           UnitsCopy( &(p->stiff->spRowUnits[4]), &(p->stiff->spRowUnits[3]) );
           UnitsCopy( &(p->stiff->spRowUnits[5]), &(p->stiff->spRowUnits[3]) );

           /* node no > 1 */
           for(i = 2; i <= p->nodes_per_elmt; i++) {
               for(j = 1; j <= p->dof_per_node; j++) {
                   k = p->dof_per_node*(i-1)+j;
                   if(j <= 3) {
                      UnitsCopy( &(p->stiff->spColUnits[k-1]), &(p->stiff->spColUnits[0]) );
                      UnitsCopy( &(p->stiff->spRowUnits[k-1]), &(p->stiff->spRowUnits[0]) );
                   }
                   else {
                      UnitsCopy( &(p->stiff->spColUnits[k-1]), &(p->stiff->spColUnits[3]) );
                      UnitsCopy( &(p->stiff->spRowUnits[k-1]), &(p->stiff->spRowUnits[3]) );
                   }
               }
           }
           free((char *) dimen1->units_name);
           free((char *) dimen1);
           free((char *) dimen2->units_name);
           free((char *) dimen2);
           free((char *) dimen3->units_name);
           free((char *) dimen3);

           break;
        case OFF:
           break;
        default:
           break;
     }


#ifdef DEBUG
    printf(" In elmt_shell_4nodes_q(): \n");
    printf("    : leaving case of MASS_MATRIX\n");
#endif
        
        break;
           
        default:
        break;
    }

/*
        MatrixFreeIndirectDouble(yl, 3);
        MatrixFreeIndirectDouble(tr, 3);
        MatrixFreeIndirectDouble(gshp1, 3);
        MatrixFreeIndirectDouble(gshp2, 3);
        MatrixFreeIndirectDouble(vl, 6);
        MatrixFreeIndirectDouble(btd, 3);
        MatrixFreeIndirectDouble(s, nst);
 
        free ((char *) eps);
        free ((char *) sigi);
        free ((char *) sigb);
        free ((char *) sg);
        free ((char *) tg);
        free ((char *) wg);
        free ((char *) gg);
        free ((char *) dvl);
        free ((char *) pp);
        free ((char *) d);
        free ((char *) bl);
*/

#ifdef DEBUG
    printf("*** leaving elmt_shell_4nodes_q() \n");
#endif

    return(p);
}


/*===================================================================*/
/* Be transfered from dktb06 by L. Jin,  Jan 16                      */
/* ==================================================================*/
     
void dktb06( sg, tg,xsj )
double sg, tg,xsj;
{
double      *el, **shn, *bd, *cd;
double      shm1, shm2;
double      *a1, *b1, *c1, *d1, *e1;
double      *a2, *b2, *c2, *d2, *e2;
int         i, j, k;
  
#ifdef DEBUG
    printf("*** Entering dktb06() ***\n");
#endif

       el  = dVectorAlloc(3);
       bd  = dVectorAlloc(3);
       cd  = dVectorAlloc(3);
       a1  = dVectorAlloc(3);
       a2  = dVectorAlloc(3);
       b1  = dVectorAlloc(3);
       b2  = dVectorAlloc(3);
       c1  = dVectorAlloc(3);
       c2  = dVectorAlloc(3);
       d1  = dVectorAlloc(3);
       d2  = dVectorAlloc(3);
       e1  = dVectorAlloc(3);
       e2  = dVectorAlloc(3);
       shn = MatrixAllocIndirectDouble( 2, 3 );

/* -------    Compute the area coordinates for the point    -------- */

       el[0] = 0.25 * (1.0 - sg) * (1.0 - tg);
       el[1] = 0.25 * (1.0 + sg) * (1.0 - tg);
       el[2] = 0.5  * (1.0 + tg);

/* -------     Form the shape function terms for trains     -------- */

       for( i=1; i<=3; i++ ) {
          bd[i-1] = b[i-1]/xsj;
          cd[i-1] = c[i-1]/xsj;
       }
       for( i=1; i<=3; i++ ) {
          j = i%3 + 1;
          k = j%3 + 1;
          shn[0][i-1] = bd[i-1] * (el[i-1] - 1.0);
          shn[1][i-1] = cd[i-1] * (el[i-1] - 1.0);
          shm1 = bd[j-1]*el[k-1] + bd[k-1]*el[j-1];
          shm2 = cd[j-1]*el[k-1] + cd[k-1]*el[j-1];
          a1[i-1] = aa[i-1] * shm1;
          a2[i-1] = aa[i-1] * shm2;
          b1[i-1] = bb[i-1] * shm1;
          b2[i-1] = bb[i-1] * shm2;
          c1[i-1] = cc[i-1] * shm1;
          c2[i-1] = cc[i-1] * shm2;
          d1[i-1] = dd[i-1] * shm1;
          d2[i-1] = dd[i-1] * shm2;
          e1[i-1] = ee[i-1] * shm1;
          e2[i-1] = ee[i-1] * shm2;
       }

/* ---------      Plate train displacement matrix       ------------ */

       for( i=1; i<=3; i++ ) {

          j = i%3 + 1;
          k = j%3 + 1;
          bm[0][2][i-1] = a1[k-1] - a1[j-1];
          bm[0][3][i-1] = b1[k-1] + b1[j-1];
          bm[0][4][i-1] = c1[k-1] + c1[j-1] - shn[0][i-1];
          bm[1][2][i-1] = d2[k-1] - d2[j-1];
          bm[1][3][i-1] =-e2[k-1] + e2[j-1] + shn[1][i-1];
          bm[1][4][i-1] =-b2[k-1] + b2[j-1];
          bm[2][2][i-1] = a2[k-1] - a2[j-1] + d1[k-1] - d1[j-1];
          bm[2][3][i-1] =-e1[k-1] - e1[j-1] + shn[0][i-1] - bm[1][4][i-1];
          bm[2][4][i-1] = c2[k-1] + c2[j-1] - shn[1][i-1] - bm[0][3][i-1];
       }
    
       free ((char *) el);
       free ((char *) bd);
       free ((char *) cd);
       MatrixFreeIndirectDouble(shn,2);
       free ((char *) a1);
       free ((char *) a2);
       free ((char *) b1);
       free ((char *) b2);
       free ((char *) c1);
       free ((char *) c2);
       free ((char *) d1);
       free ((char *) d2);
       free ((char *) e1);
       free ((char *) e2);

#ifdef DEBUG
    printf("*** leaving dktb06() ***\n");
#endif

}


/*===================================================================*/
/* Be transfered from dktq06 by L. Jin,  Jan 16                      */
/* ==================================================================*/
     
void dktq06( si )
int         si;
{
int         i, j, k;
  
#ifdef DEBUG
    printf("*** Entering dktq06() ***\n");
#endif

    /* Form strain-displacement array for DKQ element */

    for( i=1; i<=4; i++ ) {

       j = (i+2) % 4 + 1;

       bm[0][2][i-1] = aa[i-1]*shp[0][i+3][si] - aa[j-1]*shp[0][j+3][si];
       bm[0][3][i-1] = bb[i-1]*shp[0][i+3][si] + bb[j-1]*shp[0][j+3][si];
       bm[0][4][i-1] = cc[i-1]*shp[0][i+3][si] + cc[j-1]*shp[0][j+3][si] - shp[0][i-1][si];
       bm[1][2][i-1] = dd[i-1]*shp[1][i+3][si] - dd[j-1]*shp[1][j+3][si];
       bm[1][3][i-1] =-ee[i-1]*shp[1][i+3][si] - ee[j-1]*shp[1][j+3][si] + shp[1][i-1][si];
       bm[1][4][i-1] =-bb[i-1]*shp[1][i+3][si] - bb[j-1]*shp[1][j+3][si];
       bm[2][2][i-1] = aa[i-1]*shp[1][i+3][si] - aa[j-1]*shp[1][j+3][si] +
                       dd[i-1]*shp[0][i+3][si] - dd[j-1]*shp[0][j+3][si];
       bm[2][3][i-1] =-ee[i-1]*shp[0][i+3][si] - ee[j-1]*shp[0][j+3][si] +
                               shp[0][i-1][si] - bm[1][4][i-1];
       bm[2][4][i-1] = cc[i-1]*shp[1][i+3][si] + cc[j-1]*shp[1][j+3][si] -
                               shp[1][i-1][si] - bm[0][3][i-1];
    }
    
#ifdef DEBUG
    printf("*** leaving dktq06() ***\n");
#endif

}


/*===================================================================*/
/* Be transfered from hshp06 by L. Jin,  Jan 16                      */
/* ==================================================================*/
     
void hshp06( si )
int         si;
{
double      *shx, *shy;
int         i, j, k;
  
#ifdef DEBUG
    printf("*** Entering hshp06() ***\n");
#endif

/* -----    Form the shape function    ----------------------------- */

       shx = dVectorAlloc( 4 );
       shy = dVectorAlloc( 4 );

       for( k=1; k<=3; k++ ) {
          for( i=1; i<=4; i++ ) {
             shx[i-1] = shp[k-1][i+3][si]*c[i-1]; 
             shy[i-1] = shp[k-1][i+3][si]*b[i-1]; 
          }

          j = 4;

          for( i=1; i<=4; i++ ) {
             shp1[k-1][i-1][si] = shy[i-1] - shy[j-1];
             shp2[k-1][i-1][si] = shx[i-1] - shx[j-1];
             j = i;
          }
       }
    
       free ((char *) shx );
       free ((char *) shy );
       
#ifdef DEBUG
    printf("*** leaving hshp06() ***\n");
#endif

}


/*===================================================================*/
/* Be transfered from jacq06 by L. Jin,  Jan 16                      */
/* ==================================================================*/
     
void jacq06( xl )
double      **xl;
{
double      cxx, cyy, cxy, dz, ad, a4, b2, c2, sql, x31, y31, x42, y42;
int         i, j, k;
  
#ifdef DEBUG
    printf("*** Entering jacq06() ***\n");
#endif

/* -----    Form geometric constants for DKQ element   ------- */

       cxx = 0.0;
       cyy = 0.0;     
       cxy = 0.0;     

       for( i=1; i<=4; i++ ) {
          k = i % 4 + 1;
          b[i-1] = xl[1][k-1] - xl[1][i-1];
          c[i-1] = xl[0][i-1] - xl[0][k-1];
          dz     = xl[2][k-1] - xl[2][i-1];
          b2     = b[i-1]*b[i-1];
          c2     = c[i-1]*c[i-1];
          cxx   += dz*b2/(b2 + c2);
          cyy   += dz*c2/(b2 + c2);
          cxy   += (dz+dz)*b[i-1]*c[i-1]/(b2 + c2);
          sql    = (b2 + c2)/0.75;
         aa[i-1] = (c[i-1] + c[i-1])/sql;
         bb[i-1] =  b[i-1] * c[i-1] /sql;
         cc[i-1] =  c2/sql;
         dd[i-1] =-(b[i-1] + b[i-1])/sql;
         ee[i-1] =  b2/sql;
          b[i-1] =  b[i-1]/8.0;
          c[i-1] =  c[i-1]/8.0;
       }

       if( xl[2][0] != 0.0 ) {
          x31 = xl[0][2] - xl[0][0];
          y31 = xl[1][2] - xl[1][0];
          x42 = xl[0][3] - xl[0][1];
          y42 = xl[1][3] - xl[1][1];
          a4  = (x31*y42 - x42*y31)*2.0;
          ad  = a4*a4/4.0;
          x31 = x31/ad;
          y31 = y31/ad;
          x42 = x42/ad;
          y42 = y42/ad;
          bm[0][0][0] = -cxx*x42;
          bm[1][0][0] = -cyy*x42;
          bm[2][0][0] = -cxy*x42;

          bm[0][1][0] = -cxx*y42;
          bm[1][1][0] = -cyy*y42;
          bm[2][1][0] = -cxy*y42;

          bm[0][5][0] = -cxx/a4;
          bm[1][5][0] = -cyy/a4;
          bm[2][5][0] = -cxy/a4;

          bm[0][0][1] =  cxx*x31;
          bm[1][0][1] =  cyy*x31;
          bm[2][0][1] =  cxy*x31;

          bm[0][1][1] =  cxx*y31;
          bm[1][1][1] =  cyy*y31;
          bm[2][1][1] =  cxy*y31;

          bm[0][5][1] = -cxx/a4;
          bm[1][5][1] = -cyy/a4;
          bm[2][5][1] = -cxy/a4;

          for( i=1; i<=3; i++ ) {
             bm[i-1][0][2] = -bm[i-1][0][0];
             bm[i-1][1][2] = -bm[i-1][1][0];
             bm[i-1][5][2] =  bm[i-1][5][0];

             bm[i-1][0][3] = -bm[i-1][0][1];
             bm[i-1][1][3] = -bm[i-1][1][1];
             bm[i-1][5][3] =  bm[i-1][5][1];
          }
       }
    
#ifdef DEBUG
    printf("*** leaving jacq06() ***\n");
#endif

}


/*===================================================================*/
/* Be transfered from jtri06 by L. Jin,  Jan 16                      */
/* ==================================================================*/
     
void jtri06( xl, xsj )
double      **xl, *xsj;
{
double      sql, cs, bs;
int         i, j, k;
  
#ifdef DEBUG
    printf("*** Entering jtri06() ***\n");
#endif

/* -----    Form the terms for the Jacobian of a triangle   ------- */

       for( i=1; i<=3; i++ ) {
          j = i%3 + 1;
          k = j%3 + 1;
          b[i-1] = xl[1][k-1] - xl[1][j-1];
          c[i-1] = xl[0][j-1] - xl[0][k-1];
          sql    = 4.0*(b[i-1]*b[i-1] + c[i-1]*c[i-1]);
          cs     = c[i-1]/sql;
          bs     = b[i-1]/sql;
         aa[i-1] = 6.0*cs;
         bb[i-1] = 3.0*bs*c[i-1];
         cc[i-1] = c[i-1]*cs - b[i-1]*(bs+bs);
         dd[i-1] =-6.0*bs;
         ee[i-1] = b[i-1]*bs - c[i-1]*(cs+cs);
       }

       b[3]  = 0.0;
       c[3]  = 0.0;
       aa[3] = 0.0;
       bb[3] = 0.0;
       cc[3] = 0.0;
       dd[3] = 0.0;
       ee[3] = 0.0;

/* -----   Store jacobian   ----------------------------------------*/

       *xsj   = - xl[0][0]*b[0] - xl[0][1]*b[1] - xl[0][2]*b[2];

#ifdef DEBUG
    printf("*** leaving jtri06() ***\n");
#endif

}


/*===================================================================*/
/* Be transfered from proj06 by L. Jin,  Jan 16                      */
/* ==================================================================*/
     
void proj06( i1, j1, s, zi, zj, nst )
double      **s;
double      zi, zj;
int         i1, j1, nst; 
{
int         i;
  
#ifdef DEBUG
    printf("*** Entering proj06() ***\n");
#endif

/* -----    Modify the stiffness for offset projections      ------- */
/* -----    Postmultiply by transformation      ------------ ------- */

       for( i=1; i<=6; i++ ) {
          s[i1+i-1][j1+3] += zj*s[i1+i-1][j1+1];
          s[i1+i-1][j1+4] -= zj*s[i1+i-1][j1  ];
       }

/* --- Premultiply using modified terms from prosmultiplication  --- */

       for( i=1; i<=6; i++ ) {
          s[i1+3][j1+i-1] += zi*s[i1+1][j1+i-1];
          s[i1+4][j1+i-1] -= zi*s[i1  ][j1+i-1];
       }

#ifdef DEBUG
    printf("*** leaving proj06() ***\n");
#endif

}


/*===================================================================*/
/* Be transfered from rots06 by L. Jin,  Jan 16                      */
/* ==================================================================*/
     
void rots06( s, p, t, nst, ndf )
double      **s, *p, **t;
int         nst, ndf; 
{
double      **a, *b;
int         i, i0, i1, ir, ii, j, j0, j1, jc, jj, k;

#ifdef DEBUG
    printf("*** Entering rots06() ***\n");
#endif

/* -----  Transform the loads and stiffness to global coords.  ----- */

       a = MatrixAllocIndirectDouble(3,3);
       b = dVectorAlloc(6);

       i0 = 0;
       for( ir=1; ir<=4; ir++ ) {
          for( ii=1; ii<=3; ii++ ) {
             b[ii-1] = 0.0;
             b[ii+2] = 0.0;
             for( k=1; k<=3; k++ ) {
                b[ii-1] += t[k-1][ii-1]*p[i0+k-1];
                b[ii+2] += t[k-1][ii-1]*p[i0+k+2];
             }
          }
          for( ii=1; ii<=6; ii++ ) {
             p[i0+ii-1] = b[ii-1];
          }
          j0 = i0;
          for( jc=ir; jc<=4; jc++ ) {
             i1 = i0;
             for( i=1; i<=2; i++ ) {
                j1 = j0;
                for( j=1; j<=2; j++ ) {
                   for( ii=1; ii<=3; ii++ ) {
                      for( jj=1; jj<=3; jj++ ) {
                         a[jj-1][ii-1] = 0.0;
                         for( k=1; k<=3; k++ ) {
                            a[jj-1][ii-1] += t[k-1][ii-1]*s[i1+k-1][jj+j1-1];
                         }
                      }
                   }
                   for( ii=1; ii<=3; ii++ ) {
                      for( jj=1; jj<=3; jj++ ) {
                         s[ii+i1-1][jj+j1-1] = 0.0;
                         for( k=1; k<=3; k++ ) {
                            s[ii+i1-1][jj+j1-1] += a[k-1][ii-1]*t[k-1][jj-1];
                         }
                      }
                   }
                j1 = j1 + 3;
                }
             i1 = i1 + 3;
             }

/* -----   Compute the symmetric block    -------------------------- */
 
             if( ir != jc ) {
                for( i=1; i<=6; i++ ) {
                   for( j=1; j<=6; j++ ) {
                      s[j0+j-1][i0+i-1] = s[i0+i-1][j0+j-1];
                   }
                }
             }
             j0 = j0 + ndf;
          }
          i0 = i0 + ndf;
       }

       MatrixFreeIndirectDouble(a,3);
       free ((char *) b);

#ifdef DEBUG
    printf("*** leaving rots06() ***\n");
#endif

}


/*===================================================================*/
/* Be transfered from rshp06 by L. Jin,  Jan 16                      */
/* ==================================================================*/
     
void rshp06( si, ss, tt, x, xsj, ndm )
double      **x;
double      ss, tt, *xsj;
int         si, ndm; 
{
/*  
double      s[4] = {-0.5, 0.5, 0.5, -0.5};
double      t[4] = {-0.5, -0.5, 0.5, 0.5};   
*/
double      **sx;
double      s2, t2, tp, xsj1;
int         i;
        
#ifdef DEBUG
    printf("*** Entering rshp06() ***\n");
#endif

/* -----  Form 4-node quatrilateral shape functions    ------------- */

       sx = MatrixAllocIndirectDouble( 2, 2 );

       shp[0][1][si] = 0.25*(1.0-tt);
       shp[0][2][si] = 0.25*(1.0+tt);
       shp[0][0][si] = -shp[0][1][si];
       shp[0][3][si] = -shp[0][2][si];
       shp[1][3][si] = 0.25*(1.0-ss);
       shp[1][2][si] = 0.25*(1.0+ss);
       shp[1][1][si] = -shp[1][2][si];
       shp[1][0][si] = -shp[1][3][si];
       shp[2][0][si] = shp[0][1][si]*(1.0-ss);
       shp[2][1][si] = shp[0][1][si]*(1.0+ss);
       shp[2][2][si] = shp[0][2][si]*(1.0+ss);
       shp[2][3][si] = shp[0][2][si]*(1.0-ss);
       s2        = (1.0-ss*ss)*0.5;
       t2        = (1.0-tt*tt)*0.5;
       shp[2][4][si] = s2*(1.0-tt);
       shp[2][5][si] = t2*(1.0+ss);
       shp[2][6][si] = s2*(1.0+tt);
       shp[2][7][si] = t2*(1.0-ss);
       shp[0][4][si] =-ss*(1.0-tt);
       shp[0][5][si] = t2;
       shp[0][6][si] =-ss*(1.0+tt);
       shp[0][7][si] =-t2;
       shp[1][4][si] =-s2;
       shp[1][5][si] =-tt*(1.0+ss);
       shp[1][6][si] = s2;
       shp[1][7][si] =-tt*(1.0-ss);

/* -----  Construct the jacobian and its inverse    ---------------- */

       sx[1][1] = x[0][0]*shp[0][0][si] + x[0][1]*shp[0][1][si] + 
                  x[0][2]*shp[0][2][si] + x[0][3]*shp[0][3][si];
       sx[0][1] = x[0][0]*shp[1][0][si] + x[0][1]*shp[1][1][si] + 
                  x[0][2]*shp[1][2][si] + x[0][3]*shp[1][3][si];
       sx[1][0] = x[1][0]*shp[0][0][si] + x[1][1]*shp[0][1][si] + 
                  x[1][2]*shp[0][2][si] + x[1][3]*shp[0][3][si];
       sx[0][0] = x[1][0]*shp[1][0][si] + x[1][1]*shp[1][1][si] + 
                  x[1][2]*shp[1][2][si] + x[1][3]*shp[1][3][si];

       *xsj     = sx[0][0]*sx[1][1] - sx[0][1]*sx[1][0];
       xsj1     = *xsj;

       if(xsj1 == 0.0 )
          xsj1 = 1.0;

       sx[1][1] = sx[1][1]/xsj1;
       sx[0][0] = sx[0][0]/xsj1;
       sx[0][1] =-sx[0][1]/xsj1;
       sx[1][0] =-sx[1][0]/xsj1;

/* -----  Form global derivatives    ------------------------------ */

       for( i=1; i<=8; i++ ) {
          tp              = shp[0][i-1][si]*sx[0][0] + shp[1][i-1][si]*sx[1][0];
          shp[1][i-1][si] = shp[0][i-1][si]*sx[0][1] + shp[1][i-1][si]*sx[1][1];
          shp[0][i-1][si] = tp; 
       }

/* -----  Form the rotational and 5th shape functions  ------------- */

       hshp06 ( si );

        MatrixFreeIndirectDouble(sx,2);

#ifdef DEBUG
    printf("*** leaving rshp06() ***\n");
#endif

}


/*===================================================================*/
/* Be transfered from stre06 by L. Jin,  Jan 16                      */
/* ==================================================================*/
     
void stre06( p, d, vl, ndm, nel, k, xx, yy, zz, eps, sigi, sigb )
ARRAY       *p;
double      *d, **vl, *eps, *sigi, *sigb;
double      *xx, *yy, *zz;
int         ndm, nel, k; 
{
int         i, j;

#ifdef DEBUG
    printf("*** Entering stre06() ***\n");
#endif

/*
#ifdef DEBUG
    printf("d vector\n");
    for(i=1; i<=20; i++) printf("%5d%15.8e\n", i, d[i-1]);
    printf("p->coord array\n");
           printf("                   1              2              3               4\n");
    for(i=1; i<=3; i++){ 
        printf("%5d", i);
        for(j=1; j<=4; j++) 
           printf("%15.4e", p->coord[i-1][j-1].value);
        printf("\n");
    }

    dMatrixPrint("vl", vl, 6, 4);
#endif
*/

/* -----  Compute the membrane and bending strains    -------------- */

      dktq06( k-1 );

      for( i=1; i<=6; i++ ) 
         eps[i-1] = 0.0;

      *xx = 0.0;
      *yy = 0.0;
      *zz = 0.0;

      for( j=1; j<=nel; j++ ) {

         *xx += shp[2][j-1][k-1]*p->coord[0][j-1].value;
         *yy += shp[2][j-1][k-1]*p->coord[1][j-1].value;
         *zz += shp[2][j-1][k-1]*p->coord[2][j-1].value;
         eps[0] += shp[0][j-1][k-1]*vl[0][j-1] - shp1[0][j-1][k-1]*vl[5][j-1];
         eps[2] += shp[1][j-1][k-1]*vl[1][j-1] - shp2[1][j-1][k-1]*vl[5][j-1];
         eps[1] += shp[0][j-1][k-1]*vl[1][j-1] + shp[1][j-1][k-1]*vl[0][j-1] 
                   - (shp1[1][j-1][k-1] + shp2[0][j-1][k-1])*vl[5][j-1];

         for( i=ii1; i<=ii2; i++ ) {
            eps[3] += bm[0][i-1][j-1]*vl[i-1][j-1];
            eps[5] += bm[1][i-1][j-1]*vl[i-1][j-1];
            eps[4] += bm[2][i-1][j-1]*vl[i-1][j-1];
         }
      }

      sigi[0] = d[0]*eps[0] + d[1]*eps[2];
      sigi[2] = d[0]*eps[2] + d[1]*eps[0];
      sigi[1] = d[2]*eps[1];

      sigi = pstres06(sigi);

      sigb[0] = d[3]*eps[3] + d[4]*eps[5];
      sigb[2] = d[4]*eps[3] + d[3]*eps[5];
      sigb[1] = d[5]*eps[4];

      sigb = pstres06(sigb);

#ifdef DEBUG
    printf("*** leaving stre06() ***\n");
#endif

}
      

/* ----------------------------------------*/
/*     shp_prt                             */
/* ----------------------------------------*/

void shp_prt(nel)
int nel;
{
int i, j, k;

    printf("shp[3][8][9] array\n");
    for(k=1; k<=3; k++) {
    printf("k = %2d\n           1         2         3         4         5         6         7         8\n", k);
       for(i=1; i<=9; i++){ 
          printf("%2d", i);
          for(j=1; j<=8; j++) 
             printf("%10.3e", shp[k-1][j-1][i-1]);
          printf("\n");
       }
    }
    printf("shp1[3][4][9] array\n");
    for(k=1; k<=3; k++) {
    printf("k = %2d\n            1          2          3          4\n", k);
       for(i=1; i<=9; i++){ 
          printf("%2d", i);
          for(j=1; j<=4; j++) 
             printf("%11.3e", shp1[k-1][j-1][i-1]);
          printf("\n");
       }
    }
    printf("shp2[3][4][9] array\n");
    for(k=1; k<=3; k++) {
    printf("k = %2d\n            1          2          3          4\n", k);
       for(i=1; i<=9; i++){ 
          printf("%2d", i);
          for(j=1; j<=4; j++) 
             printf("%11.3e", shp2[k-1][j-1][i-1]);
          printf("\n");
       }
    }
    printf("bm[3][6][9] array\n");
    for(k=1; k<=3; k++) {
    printf("k = %2d\n            1          2          3          4          5          6\n", k);
       for(i=1; i<=9; i++){ 
          printf("%2d", i);
          for(j=1; j<=6; j++) 
             printf("%11.3e", bm[k-1][j-1][i-1]);
          printf("\n");
       }
    }

    printf("nel = %3d, ii1 = %3d, ii2 = %3d\n",nel,ii1,ii2);

}


/*===================================================================*/
/* Calculation  of  Transformation  Array  and  Surface  Coordiantes */
/* Be transfered from tran06 by L. Jin,  Jan 14                      */
/* ==================================================================*/
     
void tran06(p, yl, t)
ARRAY      *p;
double     **yl, **t;
{
int      i, j, k;
double   v1, v2, htol, d11, d12, *x0;

#ifdef DEBUG
    printf("*** Entering tran06() ***\n");
#endif

/* --- Compute the inplane direction cosines (bisect diagonals) --- */

       x0 = dVectorAlloc(3);

       for( i=1; i<=3; i++ ) {
           t[0][i-1] = p->coord[i-1][2].value - p->coord[i-1][0].value;
           t[1][i-1] = p->coord[i-1][1].value - p->coord[i-1][3].value;
       }

       d11 = sqrt(t[0][0]*t[0][0] + t[0][1]*t[0][1] + t[0][2]*t[0][2]);
       d12 = sqrt(t[1][0]*t[1][0] + t[1][1]*t[1][1] + t[1][2]*t[1][2]);
                   
       for( i=1; i<=3; i++ ) {
           v1 = t[0][i-1]/d11;
           v2 = t[1][i-1]/d12;
           t[0][i-1] = v1 + v2;
           t[1][i-1] = v1 - v2;
       }

       d11 = sqrt(t[0][0]*t[0][0] + t[0][1]*t[0][1] + t[0][2]*t[0][2]);
       d12 = sqrt(t[1][0]*t[1][0] + t[1][1]*t[1][1] + t[1][2]*t[1][2]);

       for( i=1; i<=3; i++ ) {
           t[0][i-1] = t[0][i-1]/d11;
           t[1][i-1] = t[1][i-1]/d12;

/* ------------ Compute the center (0,0) displacement ------------- */

           x0[i-1] = 0.25*(p->coord[i-1][0].value + p->coord[i-1][1].value +
                           p->coord[i-1][2].value + p->coord[i-1][3].value);
       }

/* -------------- Compute the normal to the surface --------------- */

       t[2][0] = t[0][1]*t[1][2] - t[1][1]*t[0][2] ;
       t[2][1] = t[0][2]*t[1][0] - t[1][2]*t[0][0] ;
       t[2][2] = t[0][0]*t[1][1] - t[1][0]*t[0][1] ;

/* ------- Compute the projected middle surface coordinates ------- */

       for( i=1; i<=4; i++ ) {
          for( j=1; j<=3; j++ ) {   
             yl[j-1][i-1] = 0.0;
             for( k=1; k<=3; k++ ) {
                yl[j-1][i-1] += t[j-1][k-1]*(p->coord[k-1][i-1].value - x0[k-1]);
             }
          }
       }

/* - Set offset coordinates to zero if small compared to plan size - */

       htol = 0.0;

       for( i=1; i<=4; i++ ) {
          htol = MAX( htol, ABS(yl[0][i-1]));
          htol = MAX( htol, ABS(yl[1][i-1]));
       }

       htol = htol * 1.0e-7;

       for( i=1; i<=4; i++ ) {
          if( ABS(yl[2][i-1]) <= htol )
             yl[2][i-1] = 0.0;
       }

       free ((char *) x0);

#ifdef DEBUG
    printf("*** leaving tran06() ***\n");
#endif

}


/* ============================================= */
/*          pstres06                             */
/* ============================================= */

double *pstres06(sig)
double *sig;
{
double xi1, xi2, rho;

#ifdef DEBUG
     printf("*** entering pstres06() ***\n");
#endif

     xi1 = 0.5*(sig[0] + sig[2]);
     xi2 = 0.5*(sig[0] - sig[2]);
     rho = sqrt(xi2*xi2 + sig[1]*sig[1]);
     sig[3] = xi1 + rho;
     sig[4] = xi1 - rho;
     sig[5] = 45.0;
     if(xi2 != 0.0) sig[5] = 22.5*atan2(sig[1],xi2)/atan(1.0);
     
#ifdef DEBUG
    printf("*** leaving pstres06() ***\n");
#endif

     return (sig);

}
