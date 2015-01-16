/*
 *  ============================================================================= 
 *  ALADDIN Version 1.0 :
 *          elmt_psps.c : Plane Stress Plane Strain Element
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
#include "vector.h"
#include "fe_database.h"
#include "symbol.h"
#include "fe_functions.h"

/*
#define DEBUG
*/
#define Streq(s1, s2) (strcmp(s1, s2) == 0)


/* ============================================================== */
/*   Element01                                                    */
/*   2D   Actions on Plane Stress/Plane Strain Element            */
/*        For Properties, see note below.                         */
/* ============================================================== */

save_actions_elmt01(ep,p)
ELEMENT *ep;
ARRAY    *p;
{
   /* ADD DETAILS LATER */
}


/* ================================================ */
/*   Plane-stress and plane strain Element          */
/*   elements: PLANE_STRAIN; PLANE_STRESS           */
/*   Plane Linear Elastic Element                   */
/* ================================================ */
/*    p->work_material[0] = E;
      p->work_material[1] = G;
      p->work_material[2] = fy;
      p->work_material[3] = ET;
      p->work_material[4] = nu;
      p->work_material[5] = density;
      p->work_material[6] = fu;
      p->work_material[7] = alpha_thermal[0];
      p->work_material[8] = alpha_thermal[1];

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
      p->work_section[10] = area;                                  */
/* ============================================================== */


ARRAY *elmt_psps(p,isw)
ARRAY *p;
int isw;
{
static QUANTITY fy, G, E, ET, density, lambda;
double nu;
double temperature, *alpha_thermal;

char *elmt_type;
int no_integ_pts, no_stress_pts;
int l, k, i, j,j1, k1, lint, kk;

double sg[17],tg[17],wg[17],sig[7],
       jacobain,w11,w12,w21,w22, xx, yy, dv, wd[3];
double **B_matrix;
double **B_Transpose;
double **stiff;
double **Mater_matrix;
double **body_force;
double **strain;
double **stress;
double **temp_change;
double **load;
double **displ;
double **nodal_load;
double **shp;
double **temp_matrix;

int             dof_per_elmt;
int          length1, length;
int             UNITS_SWITCH;
double         *sum_row, sum;
DIMENSIONS        *dp_stress;
DIMENSIONS     *d1, *d2, *d3;

#ifdef DEBUG
       printf("*** enter elmt_psps()\n");
#endif

    UNITS_SWITCH = CheckUnits();
    dof_per_elmt  = p->nodes_per_elmt*p->dof_per_node;  

   /* Allocation for matrices */

    alpha_thermal = dVectorAlloc(2);
    strain        = MatrixAllocIndirectDouble(3, 1);
    stress        = MatrixAllocIndirectDouble(3, 1);
    displ         = MatrixAllocIndirectDouble(dof_per_elmt, 1);
    stiff         = MatrixAllocIndirectDouble(dof_per_elmt, dof_per_elmt);
    load          = MatrixAllocIndirectDouble(dof_per_elmt, 1);
    nodal_load    = MatrixAllocIndirectDouble(dof_per_elmt, 1);
    B_matrix      = MatrixAllocIndirectDouble(3, dof_per_elmt);
    shp           = MatrixAllocIndirectDouble(3, p->nodes_per_elmt);
    sum_row       = dVectorAlloc(dof_per_elmt);

    B_Transpose   = MatrixAllocIndirectDouble(dof_per_elmt, 3);
    temp_matrix   = MatrixAllocIndirectDouble(dof_per_elmt, 3);
    body_force    = MatrixAllocIndirectDouble(3, 1);
/*
   Mater_matrix  = MatrixAllocIndirectDouble(3, 3);
   temp_change   = MatrixAllocIndirectDouble(3, 1);
*/

            E.value       =  p->work_material[0].value; 
            fy.value      =  p->work_material[2].value;
            ET.value      =  p->work_material[3].value;
            nu            =  p->work_material[4].value; 
            density.value =  p->work_material[5].value;

            if(UNITS_SWITCH == ON) {
               E.dimen       =  p->work_material[0].dimen; 
               fy.dimen      =  p->work_material[2].dimen;
               ET.dimen      =  p->work_material[3].dimen;
               density.dimen =  p->work_material[5].dimen;
            }

            alpha_thermal[0] = p->work_material[7].value; /* thermal expansion coeff. in x-dirct */
            alpha_thermal[1] = p->work_material[8].value; /* thermal expansion coeff. in y-direc */

            if(p->nodes_per_elmt == 4) no_integ_pts  =  2;/* 4-node element */
            if(p->nodes_per_elmt == 8) no_integ_pts  =  3;/* 8-node element */

            no_stress_pts = no_integ_pts;                 /* stress points, subject to change later */

            elmt_type  = p->elmt_type;          

            if(Streq("PLANE_STRAIN", elmt_type)) {
               E.value  = p->work_material[0].value = E.value /(1+nu)/(1-nu);
               nu = p->work_material[4].value = nu * E.value;
               G.value  = p->work_material[1].value = E.value/2/(1+nu);
               for (i = 1; i <= 2; i++)
               alpha_thermal[i-1]  =  alpha_thermal[i-1]*(1.0 + nu);
            }

   switch(isw) {
       case PROPTY:  /* input material properties */

            lint = 0;
            break;
       case CHERROR:
            break;
       case STIFF:

#ifdef DEBUG
       printf("*** in elmt_psps() : start case STIFF \n");
#endif

            for(i = 1; i <= dof_per_elmt; i++) 
                for(j = 1; j<= dof_per_elmt; j++)
                    stiff[i-1][j-1] = 0.0;

            if(no_integ_pts*no_integ_pts != lint)
               pgauss(no_integ_pts,&lint,sg,tg,wg);

            /* start gaussian integration  */

            for(l = 1; l<= lint; l++) { 
               shape(sg[l-1],tg[l-1],p->coord,shp,&jacobain,p->no_dimen,p->nodes_per_elmt, p->node_connect,FALSE);
            /* output:                                    */
            /*     shp[0][i-1] = dN_i/dx                  */
            /*     shp[1][i-1] = dN_i/dy                  */
            /*     compute [B] matrix at each Gaussian    */
            /*     integration point                      */

            /*****************************************************/
            /*  Derivative matrix                                */
            /*     (B_i)3x2; B_i[0][0] = dNi/dx, B_i[0][1] = 0   */
            /*     B_i[1][0] = 0,      B_i[1][1] = dNi/dy        */
            /*     B_i[2][0] = dNi/dy, B_i[2][2] = dNi/dx        */
            /*     [B] = [B_1, B_2, B_3, B_4, ..., B_n]          */
            /*       n = no of node                              */ 
            /*****************************************************/
            
            /* material matrix */
            Mater_matrix = MATER_MAT_PLANE(E.value, nu);

            /* Form [B] matrix */
            
            for(j = 1; j <= p->nodes_per_elmt; j++) { 
               k = 2*(j-1);
               B_matrix[0][k]   = shp[0][j-1];
               B_matrix[1][k+1] = shp[1][j-1];
               B_matrix[2][k]   = shp[1][j-1];
               B_matrix[2][k+1] = shp[0][j-1];
            }
           
           /* muti. Jacobain determinant with weight coefficents and mater matrix */
          
            jacobain = jacobain* wg[l-1];

            for( i = 1; i <=3 ; i++ ) 
              for (j = 1; j <= 3; j++)  
                 Mater_matrix[i-1][j-1]  *= jacobain;

           /* Transpose [B] matrix */

            for(i = 1; i <= 3; i++)
                for(j = 1; j <= dof_per_elmt; j++) 
                    B_Transpose[j-1][i-1] = B_matrix[i-1][j-1];

           /* [a] [B]^T * [Mater]  and save in [B]^T    */

           temp_matrix = dMatrixMultRep(temp_matrix, B_Transpose, dof_per_elmt, 3, Mater_matrix, 3, 3);

           /* [b] calculate stiffness : integral of [B]^T * [Mater] * [B] */

           stiff = dMatrixMultRep(stiff, temp_matrix, dof_per_elmt, 3, B_matrix, 3, dof_per_elmt);

           for (i = 1; i <= dof_per_elmt; i++) 
             for (j = 1; j <= dof_per_elmt; j++) 
                 p->stiff->uMatrix.daa[i-1][j-1]  += stiff[i-1][j-1];
           }

#ifdef DEBUG
       dMatrixPrint("element stiffness matrix", p->stiff->uMatrix.daa, dof_per_elmt, dof_per_elmt);

      /* check element stiffness : sum of row of K = 0 */

       for (i = 1; i <= dof_per_elmt; i++) 
          sum_row[i-1] =  0.0;

       for (i = 1; i <= dof_per_elmt; i++) 
          for (j = 1; j <= dof_per_elmt; j++) 
            sum_row[i-1] +=  p->stiff->uMatrix.daa[i-1][j-1]; 

       for (i = 1; i <= dof_per_elmt; i++) 
         printf(" Stiffness K_e : sum of row[%d] = %lf \n", i, sum_row[i-1]); 

#endif

     /**************************************************/
     /* Assign Units to Stiffness Matrix               */
     /**************************************************/

       /* Initiation of Stiffness Units Buffer                      */
     
       switch(UNITS_SWITCH) {
         case ON:
           if(UNITS_TYPE == SI) {
              d1 = DefaultUnits("Pa");
              d2 = DefaultUnits("m");
           }
           else {
              d1 = DefaultUnits("psi");
              d2 = DefaultUnits("in");
           }

           /* node 1 */
           UnitsMultRep( &(p->stiff->spColUnits[0]), d1, d2 );
           UnitsCopy( &(p->stiff->spColUnits[1]), &(p->stiff->spColUnits[0]) );

           /* node i  i > 1*/
           for(i = 2; i <= p->nodes_per_elmt; i++) {
                for(j = 1; j <= p->dof_per_node; j++) {
                    k  = p->dof_per_node*(i-1) + j;
                    UnitsCopy( &(p->stiff->spColUnits[k-1]), &(p->stiff->spColUnits[0]) );
                }
           }
           free((char *) d1->units_name);
           free((char *) d1);
           free((char *) d2->units_name);
           free((char *) d2);

           for( i=1 ; i<=p->size_of_stiff ; i++ )
               ZeroUnits( &(p->stiff->spRowUnits[i-1]) );

         break;
         case OFF:
         break;
         default:
         break;
       }
     
#ifdef DEBUG
       printf("*** in elmt_psps() : leaving case STIFF\n");
#endif
       
      break;
      case EQUIV_NODAL_LOAD:
      /* calculate the equivalent nodal load */
      /* due to distributed loading          */
           
#ifdef DEBUG
       printf("*** in elmt_psps() : start case EQUIV_NODAL_LOAD: \n");
#endif
      /* initilize nodal_load */
      
      for(i = 1; i <= dof_per_elmt; i++) {
          p->equiv_nodal_load->uMatrix.daa[i-1][0] = 0.0;
          nodal_load[i-1][0] = 0.0;
          load[i-1][0]= 0.0;
      }
 
     for(i = 1; i<= no_integ_pts*no_integ_pts; i++) {
         sg[i-1] = 0.0;
         tg[i-1] = 0.0;
         wg[i-1] = 0.0;
     }
      
     if(no_integ_pts*no_integ_pts != lint)
        pgauss(no_integ_pts,&lint,sg,tg,wg); 

     /* start gaussian integration */

     for(l = 1; l <= lint; l++) { 
         shape(sg[l-1],tg[l-1],p->coord,shp,&jacobain,p->no_dimen,p->nodes_per_elmt, p->node_connect,0);

            /* output: shp[0][i-1] = dN_i/dx          */
            /*         shp[1][i-1] = dN_i/dy          */
            /* compute [B] matrix at each Gaussian    */
            /* integration point                      */

            /*****************************************************/
            /*  Derivative matrix (B_i)3x2;                      */
            /*     B_i[0][0] = dNi/dx, B_i[0][1] = 0             */
            /*     B_i[1][0] = 0,      B_i[1][1] = dNi/dy        */
            /*     B_i[2][0] = dNi/dy, B_i[2][2] = dNi/dx        */
            /*     [B] = [B_1, B_2, B_3, B_4, ..., B_n]          */
            /*      n  = no of node                              */ 
            /*****************************************************/
            
            /* material matrix */
            Mater_matrix = MATER_MAT_PLANE(E.value, nu);      /* need to be changed */
    
            /* Form [B] matrix */
            
            for(j = 1; j <= p->nodes_per_elmt; j++) { 
               k = 2*(j-1);
               B_matrix[0][k]   = shp[0][j-1];
               B_matrix[1][k+1] = shp[1][j-1];
               B_matrix[2][k]   = shp[1][j-1];
               B_matrix[2][k+1] = shp[0][j-1];
            }
            
           /* muti. Jacobain determinant with weight coefficents and mater matrix */
            jacobain = jacobain* wg[l-1];
           
           /* [a] CALCULATE EQUIVALENT NODAL FORCE DUE TO INITIAL STRAIN     */

            if(p->nodal_init_strain != NULL) {

                for( i = 1; i <=3 ; i++ ) 
                  for (j = 1; j <= 3; j++)  
                    Mater_matrix[i-1][j-1]  *= jacobain;

            /* Transpose [B] matrix */

                for(i = 1; i <= 3; i++) 
                    for(j = 1; j <= dof_per_elmt; j++)
                        B_Transpose[j-1][i-1] = B_matrix[i-1][j-1];

            /* calculate strain[] at gaussian points : strain = sum Ni*nodal_strain */

                for (i = 1; i <= 3; i++)
                  strain[i-1][0]  =  0.0;
                
                for (i = 1; i <= 3; i++) {
                  for( j = 1; j <= p->nodes_per_elmt; j++) {
                    strain[i-1][0]  += shp[2][j-1] * p->nodal_init_strain[i-1][j-1];
                  }
                }
           /* [Mater]_3x3 * [strain]_3x1 and save in [stress]_3x1    */

                stress = dMatrixMultRep(stress, Mater_matrix, 3, 3, strain, 3, 1);
           
           /* mutiply [B]^T * [Mater]* [strain] */

              if(nodal_load == NULL ) { 
                  nodal_load = dMatrixMultRep(nodal_load, B_Transpose, dof_per_elmt, 3, stress, 3, 1);
              } else {
                  load = dMatrixMultRep(load, B_Transpose, dof_per_elmt, 3, stress, 3, 1);
                  for (i = 1; i<= dof_per_elmt; i++)  {
                       nodal_load[i-1][0] += load[i-1][0];
                  }
              }
            }

           /* [b] CALCULATE EQUIVALENT NODAL FORCE DUE TO INITIAL STRESS     */

            if(p->nodal_init_stress != NULL) {

            /* Transpose [B] matrix */

                for(i = 1; i <= 3; i++)
                    for(j = 1; j <= dof_per_elmt; j++)
                        B_Transpose[j-1][i-1] = B_matrix[i-1][j-1];

            /* calculate stress[] at gaussian points : stress = sum Ni*nodal_stress */

                for (i = 1; i <= 3; i++)
                  stress[i-1][0]  =  0.0;
                
                for (i = 1; i <= 3; i++) {
                  for( j = 1; j <= p->nodes_per_elmt; j++) {
                    stress[i-1][0]  += shp[2][j-1] * p->nodal_init_stress[i-1][j-1].value;
                  }
                  stress[i-1][0]  *=  jacobain;
                }

           /* mutiply [B]^T * [stress] */
               
                if(nodal_load == NULL) {
                  nodal_load  = dMatrixMultRep(nodal_load, B_Transpose, dof_per_elmt, 3, stress, 3, 1);
                } else {
                  load        = dMatrixMultRep(load, B_Transpose, dof_per_elmt, 3, stress, 3, 1);
                  for (i = 1; i<= dof_per_elmt; i++) 
                    nodal_load[i-1][0] -= load[i-1][0];
                  }
            }
           /* [c] CALCULATE EQUIVALENT NODAL FORCE DUE TO BODY FORCE         */

            if(p->nodal_body_force != NULL) {

            /* calculate body_force[] at gaussian points : body_force = sum Ni*body_force */

                for (i = 1; i <= p->no_dimen; i++)
                  body_force[i-1][0]  =  0.0;
                
                for(i = 1; i <= p->no_dimen; i++) {
                    for( j = 1; j <= p->nodes_per_elmt; j++) {
                        body_force[i-1][0] 
                            += shp[2][j-1] * p->nodal_body_force[i-1][j-1].value*jacobain;

                    /* mutiply [N]^T * [body_force] */

                        k = (j-1)*p->no_dimen + i;
                        if(nodal_load == NULL) {
                           nodal_load[k-1][0] =  shp[2][j-1]*body_force[i-1][0];
                        } else
                           nodal_load[k-1][0] += shp[2][j-1]*body_force[i-1][0];
                    }
                }
            }


           /* [d] CALCULATE EQUIVALENT NODAL FORCE DUE TO TEMPERATURE CHANGE */ 

            if(p->nodal_init_strain != NULL) {
               for(i = 1; i <= 3 ; i++ ) 
                   for(j = 1; j <= 3; j++)  
                       Mater_matrix[i-1][j-1] *= jacobain;

            /* Transpose [B] matrix */

                for(i = 1; i <= 3; i++)
                    for(j = 1; j <= dof_per_elmt; j++)
                        B_Transpose[j-1][i-1] = B_matrix[i-1][j-1];

           /*  [B]^T * [Mater]  and save in [B]^T    */

                temp_matrix =
                   dMatrixMultRep(temp_matrix, B_Transpose, dof_per_elmt, 3, Mater_matrix, 3, 3);
               
            /* calculate strain[] at gaussian points : strain = sum Ni*nodal_strain */

                for(i = 1; i <= 3; i++)
                    strain[i-1][0]  =  0.0;
                
                for(j = 1; j <= p->nodes_per_elmt; j++) {
                    temperature  += shp[2][j-1] * p->nodal_temp[j-1].value;
                }
                for(i = 1; i <= 2; i++) {
                    strain[i-1][0]  = temperature *alpha_thermal[i-1];
                }
                strain[2][0]    = 0.0;
                    
           /* mutiply [B]^T * [Mater]* [strain] */
               
                if(nodal_load == NULL) {
                  nodal_load  = dMatrixMultRep(nodal_load, temp_matrix, dof_per_elmt, 3, strain, 3, 1);
                } else {
                  load        = dMatrixMultRep(load, temp_matrix, dof_per_elmt, 3, strain, 3, 1);
                  for(i = 1; i<= dof_per_elmt; i++) 
                      nodal_load[i-1][0] += load[i-1][0];
                }
            }

           /* [e] CALCULATE EQUIVALENT NODAL FORCE DUE TO SURFACE TRACTION   */

           /* add code later */

          /* Transfer nodal_load to p->equiv_nodal_load */

                 for(i = 1; i <= dof_per_elmt; i++) {
                     p->equiv_nodal_load->uMatrix.daa[i-1][0] +=  nodal_load[i-1][0];
                 }
          }
            
   /* ------------ EQUIVALENT NODAL LOAD UNITS ------------*/
   /* Young's Modulus E is in SI then Use SI, else use US  */
   /* ---------------------------------------------------- */
    switch(UNITS_SWITCH) {
      case ON:
       /* Initiation of Equivalent nodal load Units Buffer */

       if(UNITS_TYPE == SI)
          d1 = DefaultUnits("N");
       else
          d1 = DefaultUnits("lbf");

       /* node 1 */
       UnitsCopy( &(p->equiv_nodal_load->spRowUnits[0]), d1 ); 
       UnitsCopy( &(p->equiv_nodal_load->spRowUnits[1]), d1 );

       /* node i  i > 1*/
       for ( i = 2; i <= p->nodes_per_elmt; i++) {
             for( j = 1; j <= p->dof_per_node; j++) {
                  k  = p->dof_per_node*(i-1) + j;
                  UnitsCopy( &(p->equiv_nodal_load->spRowUnits[k-1]), d1 ); 
             }
       }

       ZeroUnits( &(p->equiv_nodal_load->spColUnits[0]) );

       free((char *) d1->units_name);
       free((char *) d1);
      break;
      case OFF:
      break;
      default:
      break;
    }
     
#ifdef DEBUG
       printf("*** in elmt_psps() : leaving case EQUIV_NODAL_LOAD: \n");
#endif
           break;

      case STRESS_UPDATE:
           break;
      case STRESS:
      case LOAD_MATRIX:
	
#ifdef DEBUG
       printf("*** in elmt_psps() : start case STRESS: or case LOAD_MATRIX: \n");
#endif

           lint = (int ) 0;
           if(isw == STRESS)      l=  no_stress_pts; /* stress pts         */
           if(isw == LOAD_MATRIX) l=  no_integ_pts;  /* guassian integ pts */

           /* initilize nodal_load */
      
           for (i = 1; i <= dof_per_elmt; i++) {
              nodal_load[i-1][0] = 0.0;
              load[i-1][0]= 0.0;
           }

           for (i = 1; i<= no_integ_pts*no_integ_pts; i++) {
              sg[i-1] = 0.0;
              tg[i-1] = 0.0;
              wg[i-1] = 0.0;
           }

          if(l*l != lint)
              pgauss(l,&lint,sg,tg,wg);

           /* element stress, strains and forces */

            if(p->no_dimen == 2 && isw == STRESS) {

                printf("\n STRESS in  Element No  %d \n", p->elmt_no);
                printf(" ================================================================================================== \n");
                printf(" Gaussion    xi        eta         x         y          Stress-11       Stress-22        Stress-12 \n");
                if(UNITS_SWITCH == OFF)
                   printf("  Points \n");
                if(UNITS_SWITCH == ON){
                   if(UNITS_TYPE == SI)
                      dp_stress = DefaultUnits("Pa");
                   else
                      dp_stress = DefaultUnits("psi");
                   printf("  Points                           %s         %s             %s             %s               %s\n", 
                             p->coord[0][0].dimen->units_name,
                             p->coord[1][0].dimen->units_name,
                             dp_stress->units_name,
                             dp_stress->units_name,
                             dp_stress->units_name);
                   free((char *) dp_stress->units_name);
                   free((char *) dp_stress);
                }
                printf(" --------------------------------------------------------------------------------------------------\n \n");
            }
#ifdef DEBUG
            if(p->no_dimen == 2 && isw == LOAD_MATRIX) {
                printf("\n NODAL LOAD in  Element No  %d \n", p->elmt_no);
                printf(" ===================================================================================\n");
                printf(" Gaussion    xi        eta         x         y          nodal_no        nodal_load  \n");
                printf("  Points \n");
                printf(" ----------------------------------------------------------------------------------- \n");
            }
#endif

           for( l= 1; l<= lint; l++) {
             /* element shape function and their derivatives */
                shape(sg[l-1], tg[l-1],p->coord,shp,&jacobain,p->no_dimen,p->nodes_per_elmt, p->node_connect,0);

            /* Form [B] matrix */
            
            for(j = 1; j <= p->nodes_per_elmt; j++) { 
               k = 2*(j-1);
               B_matrix[0][k]   = shp[0][j-1];
               B_matrix[1][k+1] = shp[1][j-1];
               B_matrix[2][k]   = shp[1][j-1];
               B_matrix[2][k+1] = shp[0][j-1];
            }
            
            Mater_matrix = MATER_MAT_PLANE(E.value, nu);
            
             /* calculate strains at guassian integretion pts */

                for(i = 1;i <= 3; i++)
                    strain[i-1][0] = 0.0;

                xx = 0.0;
                yy = 0.0;

                for(j = 1; j <= p->nodes_per_elmt; j++) {
                    xx = xx + shp[2][j-1] * p->coord[0][j-1].value;
                    yy = yy + shp[2][j-1] * p->coord[1][j-1].value;

                /*  converting p->displ into a array */
                    
                   for ( k = 1; k <= p->dof_per_node; k++) {
                      j1 = p->dof_per_node*(j-1) + k; 
                      displ[j1-1][0] = p->displ->uMatrix.daa[k-1][j-1];
                   }
                }

                strain = dMatrixMultRep(strain, B_matrix, 3, dof_per_elmt, displ, dof_per_elmt, 1);
                
             /* compute stress */

                stress = dMatrixMultRep(stress, Mater_matrix, 3, 3, strain, 3, 1);

                if(isw == LOAD_MATRIX) {       
                   for(i = 1; i <= 3; i++)
                       for(j = 1; j <= dof_per_elmt; j++)
                           B_Transpose[j-1][i-1] = B_matrix[i-1][j-1];

                  /* compute eqivalent node forces due to stress */
                  /* f_equiv = int [B]^T stress dV               */ 

                   dv = jacobain*wg[l-1];

                   for (i = 1; i<= 3; i++)
                     stress[i-1][0] *= dv;

                     nodal_load = dMatrixMultRep(nodal_load, B_Transpose, dof_per_elmt, 3, stress, 3, 1);
                }
                else { /* case STRESS */

                  /* output stresses  */
                  printf("   %d  %10.4f %10.4f %10.4f %10.4f", l, sg[l-1], tg[l-1], xx, yy);
                  printf("\t%10.4f\t%10.4f\t%10.4f\n", stress[0][0], stress[1][0], stress[2][0]);  
                }
               
                for(i = 1; i <= dof_per_elmt; i++) 
                   p->nodal_loads[i-1].value += nodal_load[i-1][0];
            }

   /* ------------NODAL LOAD UNITS ------------------------*/
   /* The units type is determined by the SetUnitsType()   */
   /* ---------------------------------------------------- */
    switch(UNITS_SWITCH) {
      case ON:
        if(UNITS_TYPE == SI)
           d1 = DefaultUnits("N");
        else
           d1 = DefaultUnits("lbf");

        /* node no 1 */
        UnitsCopy( p->nodal_loads[0].dimen, d1 );
        UnitsCopy( p->nodal_loads[1].dimen, d1 );
        /* node no > 1 */
        for(i = 2; i <= p->nodes_per_elmt; i++) {    
            for(j = 1; j <= p->dof_per_node; j++) {
                k = p->dof_per_node*(i-1)+j;
                UnitsCopy( p->nodal_loads[k-1].dimen, d1 );
            }
        }
        free((char *) d1->units_name);
        free((char *) d1);
      break;
      case OFF:
      break;
      default:
      break;
    }

 /* ------------====================-------------------- */
            break;
       case MASS_MATRIX:
             
            /* compute consistent mass matrix      */
            /* p->type should be -1 for consistent */

            l = no_integ_pts; 

           printf(" no_integ_pts = %d \n", l);
           for (i = 1; i<= no_integ_pts*no_integ_pts; i++) {
              sg[i-1] = 0.0;
              tg[i-1] = 0.0;
              wg[i-1] = 0.0;
           }

            if(l*l != lint)
               pgauss(l,&lint,sg,tg,wg);

            for(l=1; l<= lint; l++) {

             /* compute shape functions */

                shape(sg[l-1],tg[l-1],p->coord,shp,&jacobain,p->no_dimen, p->nodes_per_elmt,p->node_connect,0);

                dv = density.value * wg[l-1] * jacobain;

             /* for each node compute db = shape * dv  */
                j1 = 1;
                for(j = 1; j<= p->nodes_per_elmt; j++){
                    w11 = shp[2][j-1] * dv;

                 /* compute lumped mass */
                 /* ?? store lumped mass in p->nodal_loads ?? */
                    p->nodal_loads[j1-1].value = p->nodal_loads[j1-1].value + w11; 

                 /* for each node k compute mass matrix ( upper triangular part ) */
                    k1 = j1;
                    for(k = j; k <= p->nodes_per_elmt; k++) {
                        stiff[j1-1][k1-1] += shp[2][k-1] * w11;
                        k1 = k1 + p->dof_per_node;
                    }
                    j1 = j1 + p->dof_per_node;
                } 
                      
                for (i = 1; i <= dof_per_elmt; i++) {
                   for (j = 1; j <= dof_per_elmt; j++) {
                      p->stiff->uMatrix.daa[i-1][j-1] += stiff[i-1][j-1];
                   }
                }
                 
            }

            /* compute missing parts and lower part by symmetries */

            dof_per_elmt = p->nodes_per_elmt* p->dof_per_node;

            for(j = 1; j <= dof_per_elmt; j++){
                /* ??????? Store lumped mass in p->nodal_loads ?? */
                p->nodal_loads[j].value = p->nodal_loads[j-1].value;
                for(k = j; k <= p->dof_per_node; k = k + p->dof_per_node) {
                    p->stiff->uMatrix.daa[j][k]      = p->stiff->uMatrix.daa[j-1][k-1];
                    p->stiff->uMatrix.daa[k-1][j-1]  = p->stiff->uMatrix.daa[j-1][k-1];
                    p->stiff->uMatrix.daa[k][j]      = p->stiff->uMatrix.daa[j-1][k-1];
                }  
            }

#ifdef DEBUG
      /* check element mass : sum of row of Mass */

       for (i = 1; i <= dof_per_elmt; i++)
          sum_row[i-1] =  0.0;

       for (i = 1; i <= dof_per_elmt; i++)
          for (j = 1; j <= dof_per_elmt; j++)
            sum_row[i-1] +=  p->stiff->uMatrix.daa[i-1][j-1];

       sum = 0;
       for (i = 1; i <= dof_per_elmt; i++){
         printf(" Mass M_e : sum of row[%d] = %lf \n", i, sum_row[i-1]);
         sum += sum_row[i-1];
       }
         printf(" Mass M_e : sum = %lf \n",sum);


#endif

 /* ------------ MASS UNITS ---------------------------- */
 /* Initiation of Mass Units Buffer                      */

    switch(UNITS_SWITCH) {
      case ON:
        if(UNITS_TYPE == SI || UNITS_TYPE == SI_US ) {
           d1 = DefaultUnits("Pa");
           d1 = DefaultUnits("m");
        }
        else {
           d1 = DefaultUnits("psi");
           d1 = DefaultUnits("in");
        }
        d3 = DefaultUnits("sec");

        /* node no 1 */
        UnitsMultRep( &(p->stiff->spColUnits[0]), d1, d2 );
        UnitsCopy( &(p->stiff->spColUnits[1]), &(p->stiff->spColUnits[0]) );

        UnitsPowerRep( &(p->stiff->spRowUnits[0]), d3, 2.0, NO );
        UnitsCopy( &(p->stiff->spRowUnits[1]), &(p->stiff->spRowUnits[0]) );

        /* node no > 1 */
        for(i = 2; i <= p->nodes_per_elmt; i++) {    
            for(j = 1; j <= p->dof_per_node; j++) {
                k = p->dof_per_node*(i-1)+j;
                UnitsCopy( &(p->stiff->spColUnits[k-1]), &(p->stiff->spColUnits[0]) );
                UnitsCopy( &(p->stiff->spRowUnits[k-1]), &(p->stiff->spRowUnits[0]) );
            }
        }
        free((char *) d1->units_name);
        free((char *) d1);
        free((char *) d2->units_name);
        free((char *) d2);
        free((char *) d3->units_name);
        free((char *) d3);
      break;
      case OFF:
      break;
      default:
      break;
    }
     
 /* ------------====================-------------------- */

#ifdef DEBUG
       printf("*** leaving elmt_psps() case : MASS_MATRIX  \n");
#endif
            break;
            default:
            break;
    }
   
    free(alpha_thermal);
    free(alpha_thermal);
    free(sum_row);
    MatrixFreeIndirectDouble(strain, 3);
    MatrixFreeIndirectDouble(stress, 3);
    MatrixFreeIndirectDouble(body_force, 3);
    MatrixFreeIndirectDouble(load, dof_per_elmt);
    MatrixFreeIndirectDouble(nodal_load, dof_per_elmt);
    MatrixFreeIndirectDouble(displ, dof_per_elmt);
    MatrixFreeIndirectDouble(stiff, dof_per_elmt);
    MatrixFreeIndirectDouble(B_matrix, 3);
    MatrixFreeIndirectDouble(B_Transpose, dof_per_elmt);
    MatrixFreeIndirectDouble(temp_matrix, dof_per_elmt);
    MatrixFreeIndirectDouble(shp, 3);
    

    return(p);
}


/* ================================================== */
/* shape function                                     */
/* ================================================== */

int shape(ss,tt,coord,shp,jacobian,no_dimen,nodes_per_elmt,node_connect,flg)
double  ss, tt, **shp,*jacobian, *node_connect;
QUANTITY                          **coord;
int              no_dimen, nodes_per_elmt, flg;
{
int i,j,k;
double  s[4], t[4], xs[2][2], sx[2][2], tp;
double  **shp_temp;
 
    shp_temp = MatrixAllocIndirectDouble(2, 4);

    s[0] =  0.5; s[1] = -0.5;
    s[2] = -0.5; s[3] =  0.5;

    t[0] =  0.5; t[1] =  0.5;
    t[2] = -0.5; t[3] = -0.5;

  switch(nodes_per_elmt) { 
    case 3: 
    case 4:

    /* form 4-node quadrilateral shape function                  */
    /* shape function:                  Ni = shape[2][i]         */
    /*                                  node no. i = 1, .. 4     */
    /* derivatives of shape functions:  dNi/d(ss) = shape[0][i]  */
    /*                                  dNi/d(tt) = shape[1][i]  */

    for(i = 1; i <= 4; i++){
        shp[2][i-1]      = (0.5 + s[i-1] * ss) * ( 0.5 + t[i-1] * tt); 
        shp_temp[0][i-1] = s[i-1] * (0.5 + t[i-1] * tt);                
        shp_temp[1][i-1] = t[i-1] * (0.5 + s[i-1] * ss);
    }

    /* form triangle by adding third and fourth node together */

    if(nodes_per_elmt == 3) { 
        shp[2][2] = shp[2][2] + shp[2][3];
        for(i = 0; i <= 1; i++)
            shp_temp[i][2] = shp_temp[i][2] + shp_temp[i][3];
    }

     /* construct jacobian matrix and its determinant */

     for(i = 1; i <= no_dimen; i++)      /* no_dimen = 2 */
         for(j = 1; j <= no_dimen; j++) {
             xs[i-1][j-1] = 0.0;
             for(k = 1; k <= nodes_per_elmt; k++)
                 xs[i-1][j-1] = xs[i-1][j-1] + coord[i-1][k-1].value*shp_temp[j-1][k-1];
         }

     *jacobian = xs[0][0] * xs[1][1] - xs[0][1] *xs[1][0];

    if(flg == TRUE) return;

    /* compute Jacobain inverse matrix */

    sx[0][0] = xs[1][1]/ *jacobian;
    sx[1][1] = xs[0][0]/ *jacobian;
    sx[0][1] = - xs[0][1]/ *jacobian;
    sx[1][0] = - xs[1][0]/ *jacobian;


    /* form global derivatives */

    /* save dNi/dx, dNi/dy into shp[2][node] */
  
    for(i = 1; i <= nodes_per_elmt; i++){
        shp[0][i-1] = shp_temp[0][i-1]*sx[0][0] + shp_temp[1][i-1]*sx[0][1];
        shp[1][i-1] = shp_temp[0][i-1]*sx[1][0] + shp_temp[1][i-1]*sx[1][1];
    
    }
    break;

    case 5: case 6: case 7: case 8:
       shp0(ss, tt, shp, node_connect, nodes_per_elmt);
    break;

    default:
    break;
  }
#ifdef DEBUG

   printf(" in shape () \n");
    for(j = 1; j <= nodes_per_elmt; j++){
            printf("   dN/dx[%d] = %lf ", j , shp[0][j-1]);
            printf("   dN/dy[%d] = %lf ", j , shp[1][j-1]);
            printf("       N[%d] = %lf \n", j , shp[2][j-1]);
   }
   printf(" leaving shape () \n");
#endif

   MatrixFreeIndirectDouble(shp_temp, 2);
    return;
}

/* ================================ */
/* shap0                            */
/* ================================ */

int shp0(s,t,shp,ix,nel)
double s,t,shp[4][10];
int  *ix;
int nel;
{
double s2, t2;
int i,j,k,l;

    s2 = (1- s * s)/2;
    t2 = (1 - t * t)/2;

    for(i = 5; i<= 9 ;i++)
    for(j =1;j<= 3 ; j++)
        shp[j][i] = 0.0;

    /* midside nodes (serendipty)  */

    if(ix[5] == 0) goto NEXT1;
        shp[1][5] = -s*(1-t);
        shp[2][5] = -s2;
        shp[3][5] =  s2*(1-t);

    NEXT1:
        if(nel < 6) goto NEXT7;
        if(ix[6] == 0) goto NEXT2;

        shp[1][6] = t2;
        shp[2][6] = - t *(1+s);
        shp[3][6] =   t2*(1+s);

   NEXT2:
        if(nel < 7) goto NEXT7;
        if(ix[7] == 0) goto NEXT3;

        shp[1][7] = -s*(1 + t);
        shp[2][7] = s2;
        shp[3][7] = s2*(1+t);

   NEXT3:
        if(nel < 8) goto NEXT7;
        if(ix[8] == 0) goto NEXT4;

        shp[1][8] = -t2;
        shp[2][8] = - t *(1-s);
        shp[3][8] = t2*(1-s);

     /* interior node (langragian) */
   NEXT4:
        if(nel < 9) goto NEXT7;
        if(ix[9] == 0) goto NEXT7;

        shp[1][9] = -4 * s *  t2;
        shp[2][9] = -4 * t * s2;
        shp[3][9] =  4 * s2 * t2;

     /* correct edge nodes for interior node(langrangian) */

        for(j=1;j<=3;j++) {
            for(i=1;i<=4;i++)
                shp[j][i] = shp[j][i] - 0.25 * shp[j][9];
                for(i=5;i<=8;i++)
                if(ix[i] != 0)
	           shp[j][i] = shp[j][i] - .5 *shp[j][9];
        }

     /* correct corner nodes for presence of mideside nodes */
   NEXT7:
        k = 8;
        for(i=1;i<=4;i++){
            l = i +4;
            for(j=1;j<=3;j++)
                shp[j][i] = shp[j][i] - 0.50 *(shp[j][k] + shp[j][l]);
                k = 1;
        }

        return(1);

}


/* ======================== */
/* gauss  integration       */
/* ======================== */

#ifdef __STDC__
gauss(double *sg, double *ag, int lt)
#else
gauss(sg, ag, lt)
double *sg, *ag;
int lt;
#endif
{
double t;
int    i;

#ifdef DEBUG
   printf(" ******enter gauss() \n");
   printf("       no_integ_pts = %d \n", lt);
#endif

switch(lt) {
  case 1:
         sg[1] = 0.0;
         ag[1] = 2.0;
         break;
  case 2:
         sg[1] = -1/sqrt(3.);
         sg[2] =  -sg[1];
         ag[1] = 1.0; 
         ag[2] = 1.0; 
         break;
  case 3:
         sg[1] = -sqrt(0.6);
         sg[2] = 0.0;
         sg[3] = -sg[1];

         ag[1] = 5.0/9.0; 
         ag[2] = 8.0/9.0; 
         ag[3] = ag[1]; 
         break;
  case 4:
         t = sqrt(4.8);
         sg[1] =  sqrt((3+t)/7);
         sg[2] =  sqrt((3-t)/7);
         sg[3] = -sg[2];
         sg[4] = -sg[1];
         t  = 1/3/t;
         ag[1] = 0.5-t; 
         ag[2] = 0.5+t; 
         ag[3] = ag[2]; 
         ag[4] = ag[1]; 
         break;
  default:
         break;
  }

#ifdef DEBUG
      printf(" In gauss(): \n");
   for (i = 1; i <= lt; i++) {
      printf("           : integ_coord[%d] = %lf \n",i, sg[i]);
      printf("           : weight[%d]      = %lf \n",i, ag[i]);
   }
   printf(" ******leaving gauss() \n");
#endif

  return (1);
}

/*====================================================*/
/* pgauss(no_integ_pts,lint,sg,tg,wg)                 */
/* input:                                             */
/*        no_integ_pts   no of gaussian integ pts     */
/*                       in one-direction             */
/* output:                                            */
/*        lint           no_integ_pts*no_integ_pts    */
/*        sg             coordinates in xi direction  */
/*        tg             coordinates in eta direction */
/*        wg             weight coeffcient            */
/*====================================================*/

int pgauss(l, lint, r, z, w)
int l, *lint;
double r[], z[], w[];
{
int i, j, k;
double g, h;
double g4[5], h4[5], lr[10], lz[10], lw[10];

lr[1] = lr[4] =lr[8] = -1;
lr[2] = lr[3] =lr[6] =  1;
lr[5] = lr[7] =lr[9] =  0;

lz[1] = lz[2] =lz[5] = -1;
lz[3] = lz[4] =lz[7] =  1;
lz[6] = lz[8] =lz[9] =  0;

lw[3] = lw[4] = lw[1] = lw[2] = 25;
lw[5] = lw[6] = lw[7] = lw[8] = 40;
lw[9] = 64;

    if(l < 0) {
       *lint = 3;
       g = sqrt((double) 1.0/3.0);
       h = sqrt((double) 2.0/3.0);
       r[1] = -h; r[2] =  h; r[3] = 0;
       z[1] = -g; z[2] = -g; z[3] = g;
       w[1] =  1; w[2] =  1; w[3] = 2;

       return (1);
    }

    *lint = l*l;
    switch (l) {
	case 1: /*  1 x 1 integration */
	     r[0] = 0.0;
	     z[0] = 0.0;
	     w[0] = 4.0;
	     break;
        case 2: /*  2 x 2 integration */
             g = 1.0/sqrt((double) 3.0);
             for(i = 1; i<= 4; i++) {
                 r[i-1] = g * lr[i];
                 z[i-1] = g * lz[i];
                 w[i-1] = 1.0;
             }
             break;
        case 3: /* 3 x 3 integration */
             g = sqrt((double) 0.6);
             h = 1.0/81.0;
             for(i = 1; i<= 9; i++) {
                 r[i-1] = g * lr[i];
                 z[i-1] = g * lz[i];
                 w[i-1] = h * lw[i];
             }
             break;
        case 4: /* 4 x 4 integration */
             g = sqrt((double) 4.8);
             h = sqrt((double) 30.0)/36;
             g4[1] = sqrt((double) (3+g)/7.0);
             g4[4] = -g4[1];
             g4[2] = sqrt((double) (3-g)/7.0);
             g4[3] = -g4[2];
             h4[1] = 0.5 - h;
             h4[2] = 0.5 + h;
             h4[3] = 0.5 + h;
             h4[4] = 0.5 - h;
             i = 0;
             for(j = 1; j<= 4; j++) {
                 for(k = 1; k<= 4; k++) {
                     i = i +1;
                     r[i-1] = g4[k];
                     z[i-1] = g4[j];
                     w[i-1] = h4[j]* h4[k];
                 }
             }
             break;
    }

    return(1);
}


/* ========================================== */
/* pstres                                     */
/* ========================================== */

int pstres(sig,sg4,sg5,sg6)
double *sig;
double sg4,sg5,sg6;
{
    printf("******In function pstres\n ");
    return(1);
}

/*=============================================*/
/* material matrix                             */
/*=============================================*/
double **MATER_MAT_PLANE(E, nu)
double E, nu;
{
 
double **m1;
double temp;

    m1 = MatrixAllocIndirectDouble(3, 3);
    
    temp =  E/(1.0-nu*nu);

    m1[0][0]= m1[1][1] = temp;
    m1[0][1]= m1[1][0] = nu*temp; 
    m1[2][2]= E/2.0/(1.0+nu); 

    m1[0][2] = m1[2][0] = m1[1][2] = m1[2][1] = 0.0;

    return(m1);
}


/* =================================================================== */
/* Calculate Equivalent Nodal Loads                                    */
/*                                                                     */
/* These element level procedures determine Eqvivalent Nodal Loads due */
/* to element loadings. The loads are added directly to nodal loads.   */
/*                                                                     */
/*  Negative Eqv (= Fixed end forces) added to local elmt forces       */
/*   Sign Convention:  P, by,bx,etc positive along Y,X axis            */
/*   sld01 Equivalent Loading Procedure for 2 D finite element         */
/* =================================================================== */

ARRAY *sld01( p,task)
ARRAY *p;
int task;
{
ELEMENT         *el;  
ELEMENT_LOADS  *ell; 
ELOAD_LIB      *elp; 

double shp[4][10],sg[17],tg[17],wg[17],pq[4],dj[4],xyp[4],xye[4][5];
double h, dv,dv1,dv2,ax,ay,jacobain;
int mash,gash,gish,gosh,nash,l,lint,nsfr,lne,ncord,nln;
int i,j,lnew, k,igauss;
int no_integ_pts;

    switch(task){
        case PRESSLD:
             nsfr  = 2;/* conter for variation type; para =2 */
             lne   = 4;/* no of element face */
             ncord = 3;
             nln = nsfr + 1;
             ell  =  p->elmt_load_ptr;
             for(lnew=1;lnew<=lne; lnew++){
                 h = -1;
                 for (j=1; j<=nln; j++){
                      elp  =  &ell->elib_ptr[lnew-1];
                      mash =   elp->nopl[j];
                      for(k=1; k<=ncord; k++)
                          xye[k][j] = p->coord[mash][k].value;
                 }

            if(p->nodes_per_elmt == 4) no_integ_pts  =  2;    /* 4-node element */
            if(p->nodes_per_elmt == 8) no_integ_pts  =  3;    /* 8-node element */
            l = no_integ_pts;

                 pgauss(l,lint,sg,tg,wg);

                 for(igauss=1;igauss<=lint;igauss++){
                     shape(sg[igauss],h,p->coord,shp,jacobain,p->no_dimen,p->nodes_per_elmt, p->node_connect,0);
                     for(i = 1; i <= ncord; i++){
                         gash = 0;
                         gosh = 0;
                         gish = 0;
                         for (k=1; k<=nln; k++){
                              if(i<=2){;
                                 gosh = gosh +   elp->pr->uMatrix.daa[k][i] * shp[3][k];
                                 gish = gish + xye[i][k] * shp[1][k];
                              }
                         }
                         xyp[i] = gash;
                         if(i<=2){
                            pq[i] = gosh;
                            dj[i] = gish;
                         }

                     }
                     dv = wg[igauss];
                     if(ncord == 3) dv = dv + xyp[3];
                     if(p->nodes_per_elmt == AXISYM) dv = 6.283185 * xyp[1];
                        dv1 = dv + dj[1];
                        dv2 = dv + dj[2];
                        ay = dv1*pq[1] + dv2*pq[2];
                        ax = dv1 * pq[2] - dv2 * pq[1];
                        for (i=1; i<=nln; i++){
                             nash = p->dof_per_node * ( ell->elib_ptr[lnew -1].nopl[i] - 1)  + 2; 
                             mash = nash -1;
                             F[mash] = F[mash] + shp[3][i] * ax;
                             F[nash] = F[nash] + shp[3][i] * ay;
                        }
                 }
            }
            break;
       default:
            break;
    }

    return(p);

}
