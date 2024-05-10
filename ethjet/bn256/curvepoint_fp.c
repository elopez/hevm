/*
 * bn256-20080525/curvepoint_fp.c 
 * Peter Schwabe
 * Public domain
*/

#include <stdio.h>
#include <gmp.h>
#include <stdlib.h>
#include <string.h>
#include "fp.h"
#include "curvepoint_fp.h"
#include "curve.h"
// #include "cpucycles.h"

//////////////////////////////////////////////////////////////////////////////////////////////////////////
//            Point initialization and deletion functions
//////////////////////////////////////////////////////////////////////////////////////////////////////////

// Global dummies usable by all curvepoints:
fpe_t curvepoint_dummy_fpe1;

// Initialize a curvepoint_fp:
void curvepoint_fp_init(curvepoint_fp_t rop)
{
	fpe_set_ui(rop->m_x, 0);
	fpe_set_ui(rop->m_y, 1);
	fpe_set_ui(rop->m_z, 0);
	fpe_set_ui(rop->m_t, 0);
}

void curvepoint_fp_init_set_str(curvepoint_fp_t rop, const char* x, const char* y, const char* z)
{
	fpe_set_str(rop->m_x, x);
	fpe_set_str(rop->m_y, y);
	fpe_set_str(rop->m_z, z);
	fpe_set_ui(rop->m_t, 0);
}

void curvepoint_fp_init_set_mpz(curvepoint_fp_t rop, const mpz_t x, const mpz_t y, const mpz_t z)
{
	fpe_set_mpz(rop->m_x, x);
	fpe_set_mpz(rop->m_y, y);
	fpe_set_mpz(rop->m_z, z);
	fpe_set_ui(rop->m_t, 0);
}

void curvepoint_fp_get_mpz(curvepoint_fp_t rop, mpz_t x, mpz_t y, mpz_t z)
{
	fpe_get_mpz(rop->m_x, x);
	fpe_get_mpz(rop->m_y, y);
	fpe_get_mpz(rop->m_z, z);
}

void curvepoint_fp_init_set(curvepoint_fp_t rop, const curvepoint_fp_t op)
{
	fpe_set(rop->m_x, op->m_x);
	fpe_set(rop->m_y, op->m_y);
	fpe_set(rop->m_z, op->m_z);
	fpe_set_ui(rop->m_t, 0);
}

// Set the coordinates of a curvepoint_fp:
void curvepoint_fp_set_str(curvepoint_fp_t rop, const char* x, const char* y, const char* z)
{
	fpe_set_str(rop->m_x, x);
	fpe_set_str(rop->m_y, y);
	fpe_set_str(rop->m_z, z);
	fpe_set_ui(rop->m_t, 0);
}

// Set the coordinates of a curvepoint_fp_t by copying the coordinates from another curvepoint_fp
void curvepoint_fp_set(curvepoint_fp_t rop, const curvepoint_fp_t op)
{
	fpe_set(rop->m_x, op->m_x);
	fpe_set(rop->m_y, op->m_y);
	fpe_set(rop->m_z, op->m_z);
	fpe_set_ui(rop->m_t, 0);
}

// Addition of two points, op2 is assumed to be in affine coordinates 
// For the algorithm see e.g. DA Peter Schwabe
void curvepoint_fp_mixadd(curvepoint_fp_t rop, const curvepoint_fp_t op1, const curvepoint_fp_t op2)
{
	if (!memcmp(op1, op2, sizeof(*op1))) {
		return curvepoint_fp_double(rop, op1);
	}
	fpe_t tfpe1, tfpe2, tfpe3, tfpe4, tfpe5, tfpe6, tfpe7, tfpe8, tfpe9; // Temporary variables needed for intermediary results
	fpe_square(tfpe1, op1->m_z);
	fpe_mul(tfpe2, op1->m_z, tfpe1);
	fpe_mul(tfpe3, op2->m_x, tfpe1);
	fpe_mul(tfpe4, op2->m_y, tfpe2);
	fpe_sub(tfpe5, tfpe3, op1->m_x);
	fpe_sub(tfpe6, tfpe4, op1->m_y);
	fpe_square(tfpe7, tfpe5);
	fpe_mul(tfpe8, tfpe7, tfpe5);
	fpe_mul(tfpe9, op1->m_x, tfpe7);

	fpe_double(tfpe1, tfpe9);
	fpe_add(tfpe1, tfpe1, tfpe8);
	fpe_square(rop->m_x, tfpe6);
	fpe_sub(rop->m_x, rop->m_x, tfpe1);
	fpe_sub(tfpe1, tfpe9, rop->m_x);
	fpe_mul(tfpe2, tfpe1, tfpe6);
	fpe_mul(tfpe3, op1->m_y, tfpe8);
	fpe_sub(rop->m_y, tfpe2, tfpe3);
	fpe_mul(rop->m_z, op1->m_z, tfpe5);
}

void curvepoint_fp_double(curvepoint_fp_t rop, const curvepoint_fp_t op)
{
	fpe_t tfpe1, tfpe2, tfpe3, tfpe4; // Temporary variables needed for intermediary results
	fpe_square(tfpe1, op->m_y);
	fpe_mul(tfpe2, tfpe1, op->m_x);
	fpe_double(tfpe2, tfpe2);
	fpe_double(tfpe2, tfpe2);
	fpe_square(tfpe3, tfpe1);
	fpe_double(tfpe3, tfpe3);
	fpe_double(tfpe3, tfpe3);
	fpe_double(tfpe3, tfpe3);
	fpe_square(tfpe4, op->m_x);
	fpe_triple(tfpe4, tfpe4);
	fpe_square(rop->m_x, tfpe4);
	fpe_double(tfpe1, tfpe2);
	fpe_sub(rop->m_x, rop->m_x, tfpe1);
	fpe_sub(tfpe1, tfpe2, rop->m_x);
	fpe_mul(rop->m_z, op->m_y, op->m_z);
	fpe_double(rop->m_z, rop->m_z);
	fpe_mul(rop->m_y, tfpe4, tfpe1);
	fpe_sub(rop->m_y, rop->m_y, tfpe3);
}

void curvepoint_fp_mul(curvepoint_fp_t rop, const curvepoint_fp_t op, const mpz_t scalar)
{
	size_t i;
	curvepoint_fp_t r;
	curvepoint_fp_set(r, op);
	for(i = mpz_sizeinbase(scalar, 2) - 1; i > 0; i--)
	{
		curvepoint_fp_double(r, r);
		if(mpz_tstbit(scalar, i - 1)) 
			curvepoint_fp_mixadd(r, r, op);
	}
	curvepoint_fp_set(rop, r);
}



// Negate a point, store in rop:
void curvepoint_fp_neg(curvepoint_fp_t rop, const curvepoint_fp_t op)
{
	fpe_neg(curvepoint_dummy_fpe1, op->m_y);
	fpe_set(rop->m_x, op->m_x);
	fpe_set(rop->m_y, curvepoint_dummy_fpe1);
	fpe_set(rop->m_z, op->m_z);
}

// Transform to Affine Coordinates (z=1)
void curvepoint_fp_makeaffine(curvepoint_fp_t point)
{
	if(fpe_iszero(point->m_z))
	{
		fpe_setzero(point->m_x);
		fpe_setone(point->m_y);
		fpe_setzero(point->m_z);
	}
	else
	{
		fpe_invert(curvepoint_dummy_fpe1, point->m_z);
		fpe_mul(point->m_x, point->m_x, curvepoint_dummy_fpe1);
		fpe_mul(point->m_x, point->m_x, curvepoint_dummy_fpe1);

		fpe_mul(point->m_y, point->m_y, curvepoint_dummy_fpe1);
		fpe_mul(point->m_y, point->m_y, curvepoint_dummy_fpe1);
		fpe_mul(point->m_y, point->m_y, curvepoint_dummy_fpe1);

		fpe_setone(point->m_z);
	}
}

// Print a point:
void curvepoint_fp_print(FILE *outfile, const curvepoint_fp_t point)
{
	fprintf(outfile, "[");
	fpe_print(outfile, point->m_x);
	fprintf(outfile, ", ");
	fpe_print(outfile, point->m_y);
	fprintf(outfile, ", ");
	fpe_print(outfile, point->m_z);
	fprintf(outfile, "]");
}

int curvepoint_fp_well_formed(const curvepoint_fp_t op)
{
    if (fpe_iszero(op->m_x) && fpe_isone(op->m_y) && fpe_iszero(op->m_z))
    {
        return 1;
    }
    /*
        y^2 = x^3 + b

        We are using Jacobian coordinates, so equation we need to check is actually

        (y/z^3)^2 = (x/z^2)^3 + b
        y^2 / z^6 = x^3 / z^6 + b
        y^2 = x^3 + b z^6
    */

	fpe_t X2, Y2, Z2, X3, Z3, Z6, bcoef, rhs; // Temporary variables needed for intermediary results
	fpe_square(X2, op->m_x);
	fpe_square(Y2, op->m_y);
	fpe_square(Z2, op->m_z);

	fpe_mul(X3, X2, op->m_x);
	fpe_mul(Z3, Z2, op->m_z);
	fpe_square(Z6, Z3);

	fpe_set_mpz(bcoef, b);
	fpe_mul(rhs, Z6, bcoef);
	fpe_add(rhs, rhs, X3);
	return !memcmp(rhs, Y2, sizeof(*rhs));
}