/*
 * bn256-20080525/eta_tate.c 
 * Peter Schwabe
 * Public domain
 */

#include "fp.h"
#include "fp2e.h"
#include "fp6e.h"
#include "fp12e.h"
#include "curvepoint_fp.h"
#include "twistpoint_fp2.h"
#include "init.h"
#include "points.h"
#include "final_expo.h"
#include "eta_tate.h"

#include <stdio.h>
#include <gmp.h>
#include <cpucycles.h>

static void linefunction_add_eta(
    fp12e_t rop1, 
    curvepoint_fp_t rop2, 
    const curvepoint_fp_t op1, 
    const curvepoint_fp_t op2, 
    const twistpoint_fp2_t op3,
    const fpe_t r2
    )
{
  fpe_t tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10; // Temporary variables needed for intermediary results
  fp2e_t tfp2e1, tfp2e2, tfp2e3, tfp2e4;
  fp6e_t tfp6e1, tfp6e2;

  fpe_mul(tmp0, op2->m_x, op1->m_t); /* tmp0 = B = x2 * T1  = x2z1^2*/

  fpe_add(tmp1, op2->m_y, op1->m_z);
  fpe_square(tmp1, tmp1);
  fpe_sub(tmp1, tmp1, r2);
  fpe_sub(tmp1, tmp1, op1->m_t);
  fpe_mul(tmp1, tmp1, op1->m_t); /* tmp1 = D = ((y2 + Z1)^2 - R2 - T1)T1  = 2y2z1^3 */

  fpe_sub(tmp2, tmp0, op1->m_x); /* tmp2 = H = B - X1  = x2z1^2 - x1*/

  fpe_square(tmp3, tmp2); /* tmp3 = I = H^2  = (x2z1^2 - x1)^2*/

  fpe_double(tmp4, tmp3); 
  fpe_double(tmp4, tmp4); /* tmp4 = E = 4I = 4(x2z1^2 - x1)^2*/

  fpe_mul(tmp5, tmp2, tmp4); /* tmp5 = J = HE =  4(x2z1^2 - x1)(x2z1^2 - x1)^2*/

  fpe_sub(tmp6, tmp1, op1->m_y); 
  fpe_sub(tmp6, tmp6, op1->m_y); /* tmp6 = r = 2(D - 2Y1) = (2y2z1^3 - 2y1)*/

  fpe_mul(tmp9, tmp6, op2->m_x); /* Needed later: tmp9 = x2(2y2z1^3 - 2y1)*/

  fpe_mul(tmp7, op1->m_x, tmp4); /* tmp7 = V = X1*E = 4x1(x2z1^2 - x1)^2*/

  fpe_square(rop2->m_x, tmp6);
  fpe_sub(rop2->m_x, rop2->m_x, tmp5);
  fpe_sub(rop2->m_x, rop2->m_x, tmp7);
  fpe_sub(rop2->m_x, rop2->m_x, tmp7); /* X3 = r^2 - J - 2V = (2y2z1^3 - 2y1)^2 - 4(x2z1^2 - x1)(x2z1^2 - x1)^2 - 8x1(x2z1^2 - x1)^2*/

  fpe_add(rop2->m_z, op1->m_z, tmp2);
  fpe_square(rop2->m_z, rop2->m_z);
  fpe_sub(rop2->m_z, rop2->m_z, op1->m_t);
  fpe_sub(rop2->m_z, rop2->m_z, tmp3); /* Z3 = (z1 + H)^2 - T1 - I  = 2z1(x2z1^2 - x1) */

  fpe_add(tmp10, op2->m_y, rop2->m_z); /* Needed later: tmp10 = y2 + z3*/

  fpe_sub(tmp8, tmp7, rop2->m_x);
  fpe_mul(tmp8, tmp8, tmp6);
  fpe_mul(tmp0, op1->m_y, tmp5);
  fpe_double(tmp0, tmp0);
  fpe_sub(rop2->m_y, tmp8, tmp0); /* Y3 = r(V - X3) - 2Y1*J = (2y2z1^3 - 2y1)(4x1(x2z1^2 - x1)^2 - x3) - 8y1(x2z1^2 - x1)(x2z1^2 - x1)^2*/


  fpe_square(rop2->m_t, rop2->m_z); /* T3 = Z3^2 */

  fpe_square(tmp10, tmp10); /* tmp10 = (y2 + z3)^2 */
  fpe_sub(tmp10, tmp10, r2);
  fpe_sub(tmp10, tmp10, rop2->m_t); 
  fpe_double(tmp9, tmp9);
  fpe_sub(tmp9, tmp9, tmp10); /* tmp9 = 4x2(y2z1^3 - y1) - 2z3y2 */

  fp2e_mul_fpe(tfp2e3, op3->m_y, rop2->m_z); /* tmp10 = z3y_Q */
  fp2e_double(tfp2e3, tfp2e3);

  fpe_neg(tmp6, tmp6);
  fp2e_mul_fpe(tfp2e2, op3->m_x, tmp6);
  fp2e_double(tfp2e2, tfp2e2);

  fp2e_set_fpe(tfp2e4, tmp9);

  fp2e_setzero(tfp2e1);

  fp6e_set_fp2e(tfp6e1, tfp2e1, tfp2e3, tfp2e1);
  fp6e_set_fp2e(tfp6e2, tfp2e1, tfp2e2, tfp2e4);

  // tfp2e3 = - yQz1^3 * (x2z1^2 - x1)
  // tfp2e2 = xQz1^2 * (y2z1^3 - y1)
  // tfp2e4 = y1(x2z1^2 - x1) - x1(y2z1^3 - y1)

  fp12e_set_fp6e(rop1, tfp6e1, tfp6e2);
}

static void linefunction_double_eta(fp12e_t rop1, curvepoint_fp_t rop2, const curvepoint_fp_t op1, const twistpoint_fp2_t op3)
{
  fpe_t tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp7, dummy; // Temporary variables needed for intermediary results
  fp2e_t tfp2e1, tfp2e2, tfp2e3, tfp2e4;
  fp6e_t tfp6e1, tfp6e2;

  fpe_square(tmp0, op1->m_x); /* tmp0 = A = X1^2 = x1^2 */
  fpe_square(tmp1, op1->m_y); /* tmp1 = B = Y1^2 = y1^2 */

  fpe_square(tmp2, tmp1); /* tmp2 = C = B^2 = y1^4 */

  fpe_add(tmp3, op1->m_x, tmp1);
  fpe_square(tmp3, tmp3);
  fpe_sub(tmp3, tmp3, tmp0);
  fpe_sub(tmp3, tmp3, tmp2);
  fpe_double(tmp3, tmp3); /* tmp3 = D = 2(X1 + B)^2 - A - C) = 4x1y1^2 */

  fpe_triple(tmp4, tmp0); /* tmp4 = E = 3A = 3x1^2 */

  fpe_add(tmp7, tmp4, op1->m_x); /* Needed later */

  fpe_square(tmp5, tmp4); /* tmp5 = G = E^2 = 9x1^4 */

  fpe_sub(rop2->m_x, tmp5, tmp3);
  fpe_sub(rop2->m_x, rop2->m_x, tmp3); /* X3 = G - 2D = 9x1^4 - 8x1y1^2 */

  fpe_add(rop2->m_z, op1->m_y, op1->m_z);
  fpe_square(rop2->m_z, rop2->m_z);
  fpe_sub(rop2->m_z, rop2->m_z, tmp1);
  fpe_sub(rop2->m_z, rop2->m_z, op1->m_t); /* Z3 = (Y1 + Z1)^2 - B - T1 = 2y1z1; */

  fpe_sub(rop2->m_y, tmp3, rop2->m_x);
  fpe_mul(rop2->m_y, rop2->m_y, tmp4); 
  fpe_double(dummy, tmp2);
  fpe_double(dummy, dummy);
  fpe_double(dummy, dummy);
  fpe_sub(rop2->m_y, rop2->m_y, dummy); /* Y3 = E(D - X3) - 8C = 3x1^2(4x1y1^2 - X3) - 8y1^4 */

  fpe_mul(tmp3, tmp4, op1->m_t);
  fpe_double(tmp3, tmp3);
  fpe_neg(tmp3, tmp3);
  fp2e_mul_fpe(tfp2e2, op3->m_x, tmp3); /* tfp2e2 = -6x1^2z1^2 * x_Q */

  fpe_square(tmp7, tmp7);
  fpe_sub(tmp7, tmp7, tmp0);
  fpe_sub(tmp7, tmp7, tmp5);
  fpe_double(dummy, tmp1);
  fpe_double(dummy, dummy);
  fpe_sub(tmp7, tmp7, dummy); /* tmp7 = 6x1^3 - 4y1^2 */
  fp2e_set_fpe(tfp2e4, tmp7);

  fpe_mul(tmp0, rop2->m_z, op1->m_t);
  fpe_double(tmp0, tmp0);
  fp2e_mul_fpe(tfp2e3, op3->m_y, tmp0);

  fpe_square(rop2->m_t, rop2->m_z); 

  fp2e_setzero(tfp2e1);

  fp6e_set_fp2e(tfp6e1, tfp2e1, tfp2e3, tfp2e1);
  fp6e_set_fp2e(tfp6e2, tfp2e1, tfp2e2, tfp2e4);

  fp12e_set_fp6e(rop1, tfp6e1, tfp6e2);
}


void eta(fp12e_t rop, const curvepoint_fp_t op1, const twistpoint_fp2_t op2)
{
  fp12e_setone(rop);

  fp12e_t dummy;
  fpe_t r2;

  curvepoint_fp_t r;
  curvepoint_fp_set(r, op1);
  fpe_setone(r->m_t); /* As r has to be in affine coordinates this is ok */
  
  fpe_set(r2, op1->m_y);
  fpe_square(r2, r2);

  size_t i;

  unsigned long long int t1, t2, t3;
  t1 = cpucycles();
  for(i = mpz_sizeinbase(looplength_eta, 2) - 1; i > 0; i--)
  {
    linefunction_double_eta(dummy, r, r, op2);

    fp12e_square(rop, rop);
    fp12e_mul(rop, rop, dummy);

    if (mpz_tstbit(looplength_eta, i - 1)) 
    {
      linefunction_add_eta(dummy, r, r, op1, op2, r2);
      fp12e_mul(rop, rop, dummy);
    }
  }
  t2 = cpucycles();
  final_expo(rop);
  t3 = cpucycles();
  printf("Miller: %llu\nFinal expo: %llu\n", t2-t1, t3-t2);
}

void tate(fp12e_t rop, const curvepoint_fp_t op1, const twistpoint_fp2_t op2)
{
  fp12e_setone(rop);

  fp12e_t dummy;
  fpe_t r2;

  curvepoint_fp_t r;
  curvepoint_fp_set(r, op1);
  fpe_setone(r->m_t); /* As r has to be in affine coordinates this is ok */

  fpe_set(r2, op1->m_y);
  fpe_square(r2, r2);

  size_t i;

  unsigned long long int t1, t2, t3;
  t1 = cpucycles();
  for(i = mpz_sizeinbase(n, 2) - 1; i > 0; i--)
  {
    linefunction_double_eta(dummy, r, r, op2);
    fp12e_square(rop, rop);
    fp12e_mul(rop, rop, dummy);

    if (mpz_tstbit(n, i - 1)) 
    {
      linefunction_add_eta(dummy, r, r, op1, op2, r2);
      fp12e_mul(rop, rop, dummy);
    }
  }

  t2 = cpucycles();
  final_expo(rop);
  t3 = cpucycles();
  printf("Miller: %llu\nFinal expo: %llu\n", t2-t1, t3-t2);

}

