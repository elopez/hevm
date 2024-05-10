/*
 * bn256-20080525/cometa_comtate.c 
 * Peter Schwabe
 * Public domain
*/

#include <stdio.h>
#include <gmp.h>
#include <cpucycles.h>

#include "fp.h"
#include "fp2.h"
#include "fp2e.h"
#include "fp6.h"
#include "fp6e.h"
#include "comfp12e.h"
#include "curvepoint_fp.h"
#include "twistpoint_fp2.h"
#include "init.h"
#include "points.h"
#include "final_expo.h"
#include "cometa_comtate.h"

static void linefunction_evaluate_cometa(comfp12e_t rop,
		const fpe_t slope_numerator, 	// numerator of the slope of the line function
		const fpe_t slope_denominator, 	// denominator of the slope of the line function
		const fpe_t zapow3,			// z_A^3
		const curvepoint_fp_t op1,		// One of the two points, from which the line function is constructed
		const twistpoint_fp2_t op2)		// The point on the twist, at which the linefunction is evaluated
{
	fpe_t tfpe1, tfpe2;	// Temporary variables
	fp2e_t tfp2e1, tfp2e2;	// Temporary variables

	fpe_mul(tfpe2, slope_numerator, op1->m_x);
	fpe_mul(tfpe2, tfpe2, op1->m_z);
	fpe_mul(tfpe1, slope_denominator, op1->m_y);
	fpe_sub(tfpe1, tfpe1, tfpe2);
	fpe_mul(tfpe2, cometa_c0_const, tfpe1); // C0 computed

	fpe_mul(tfpe1, slope_numerator, zapow3);
	fpe_mul(tfpe1, tfpe1, cometa_c1_const);
	fp2e_mul_fpe(tfp2e1, op2->m_x, tfpe1); // C1 computed

	fpe_mul(tfpe1, slope_denominator, zapow3);
	fp2e_mul_fpe(tfp2e2, op2->m_y, tfpe1); // C computed

	//Negate the coefficient m_a of c:
	fpe_neg(tfp2e2->m_a, tfp2e2->m_a);

	fp2e_mul_fpe(rop->m_a0, tfp2e2, tfpe2);
	fp2e_mul(rop->m_a1, tfp2e1, tfp2e2);

	fpe_square(tfpe1, tfp2e2->m_a);
	fpe_square(tfpe2, tfp2e2->m_b);
#if(ALPHA == 2)
	fpe_double(tfpe1, tfpe1);
	fpe_sub(rop->m_a, tfpe1, tfpe2);
#elif(ALPHA == -2)
	fpe_double(tfpe1, tfpe1);
	fpe_add(rop->m_a, tfpe1, tfpe2);
#elif(ALPHA == -1)
	fpe_add(rop->m_a, tfpe1, tfpe2);
#else
#error "ALPHA must be -1, 2 or -2"
#endif
}


static void linefunction_add_cometa(comfp12e_t rop1, curvepoint_fp_t rop2, const curvepoint_fp_t op1, const curvepoint_fp_t op2, const twistpoint_fp2_t op3)
{
	// op2 is assumed to be in affine coordinates!
	fpe_t tfpe1, tfpe2, tfpe3, tfpe4, tfpe5;	// Temporary variables

	// Computation of the slope (numerator and denominator), see Guide to ECC, page 91f.
	fpe_square(tfpe1, op1->m_z);
	fpe_mul(tfpe2, tfpe1, op1->m_z); // dummy2 = z_A^3

	// Store z_1^3 for use in the evaluation:
	fpe_set(tfpe3, tfpe2);	

	fpe_mul(tfpe1, tfpe1, op2->m_x); // dummy1 = z_A^2x_B
	fpe_mul(tfpe2, tfpe2, op2->m_y); // dummy2 = y_Bz_A^3
	fpe_sub(tfpe1, tfpe1, op1->m_x); // dummy1 = z_A^2x_B - x_A
	fpe_sub(tfpe2, tfpe2, op1->m_y); // dummy2 = y_Bz_A^3 - y_A
	fpe_mul(tfpe5, tfpe1, op1->m_z); // z-coordinate computed, stored in tfpe5, dummy0 = (z_A^2x_B - x_A)z_A

	// Compute the linefunction:
	linefunction_evaluate_cometa(rop1, tfpe2, tfpe5, tfpe3, op1, op3);

	// Copy z-coordinate to rop:
	fpe_set(rop2->m_z, tfpe5);

	// Continue computation of point coordinates:
	fpe_square(tfpe3, tfpe1);
	fpe_mul(tfpe4, tfpe3, tfpe1);
	fpe_mul(tfpe3, tfpe3, op1->m_x);
	fpe_double(tfpe1, tfpe3);
	fpe_square(rop2->m_x, tfpe2);
	fpe_sub(rop2->m_x, rop2->m_x, tfpe1);
	fpe_sub(rop2->m_x, rop2->m_x, tfpe4);
	fpe_sub(tfpe3, tfpe3, rop2->m_x);
	fpe_mul(tfpe3, tfpe3, tfpe2);
	fpe_mul(tfpe4, tfpe4, op1->m_y);
	fpe_sub(rop2->m_y, tfpe3, tfpe4);
}

static void linefunction_double_cometa(comfp12e_t rop1, curvepoint_fp_t rop2, const curvepoint_fp_t op1, const twistpoint_fp2_t op3)
{
	fpe_t tfpe1, tfpe2, tfpe3, tfpe4, tfpe5;	// Temporary variables

	// Computation of the slope (numerator and denominator), see DA Peter Schwabe, page 42
	fpe_square(tfpe1, op1->m_y);
	fpe_mul(tfpe2, op1->m_x, tfpe1);
	fpe_double(tfpe2, tfpe2);
	fpe_double(tfpe2, tfpe2);

	fpe_square(tfpe3, tfpe1);
	fpe_double(tfpe3, tfpe3);
	fpe_double(tfpe3, tfpe3);
	fpe_double(tfpe3, tfpe3);

	fpe_square(tfpe1, op1->m_x);
	fpe_triple(tfpe1, tfpe1);

	// Compute z_1^3 for evaluation
	fpe_square(tfpe4, op1->m_z);
	fpe_mul(tfpe4, tfpe4, op1->m_z);

	fpe_mul(tfpe5, op1->m_y, op1->m_z);
	fpe_double(tfpe5, tfpe5);

	// Compute the linefunction:
	linefunction_evaluate_cometa(rop1, tfpe1, tfpe5, tfpe4, op1, op3);

	// Copy z-coordinate
	fpe_set(rop2->m_z, tfpe5);

	fpe_square(rop2->m_x, tfpe1);
	fpe_double(tfpe5, tfpe2);
	fpe_sub(rop2->m_x, rop2->m_x, tfpe5);
	fpe_sub(rop2->m_y, tfpe2, rop2->m_x);
	fpe_mul(rop2->m_y, rop2->m_y, tfpe1);
	fpe_sub(rop2->m_y, rop2->m_y, tfpe3);

}

void cometa(comfp12e_t rop, const curvepoint_fp_t op1, const twistpoint_fp2_t op2)
{
	comfp12e_setone(rop);
	
	comfp12e_t dummy;

	curvepoint_fp_t r;
	curvepoint_fp_set(r, op1);

	size_t i;
    unsigned long long int t1, t2, t3;
    t1 = cpucycles();

	for(i = mpz_sizeinbase(looplength_eta, 2) - 1; i > 0; i--)
	{
		linefunction_double_cometa(dummy, r, r, op2);

		comfp12e_square(rop, rop);
		comfp12e_mul(rop, rop, dummy);

		if (mpz_tstbit(looplength_eta, i - 1)) 
		{
			linefunction_add_cometa(dummy, r, r, op1, op2);
			comfp12e_mul(rop, rop, dummy);
		}
	}

    t2 = cpucycles();
	final_expo_com(rop);
    t3 = cpucycles();
    printf("Miller: %llu\nFinal expo: %llu\n", t2-t1, t3-t2);
}

void comtate(comfp12e_t rop, const curvepoint_fp_t op1, const twistpoint_fp2_t op2)
{
	comfp12e_setone(rop);
	
	comfp12e_t dummy;

	curvepoint_fp_t r;
	curvepoint_fp_init_set(r, op1);

	size_t i;
    unsigned long long int t1, t2, t3;
    t1 = cpucycles();

	for(i = mpz_sizeinbase(n, 2) - 1; i > 0; i--)
	{
		linefunction_double_cometa(dummy, r, r, op2);

		comfp12e_square(rop, rop);
		comfp12e_mul(rop, rop, dummy);

		if (mpz_tstbit(n, i - 1)) 
		{
			linefunction_add_cometa(dummy, r, r, op1, op2);
			comfp12e_mul(rop, rop, dummy);
		}
	}

    t2 = cpucycles();
	final_expo_com(rop);
    t3 = cpucycles();
    printf("Miller: %llu\nFinal expo: %llu\n", t2-t1, t3-t2);

}


