/*
 * bn256-20080525/twistpoint_fp2.c 
 * Peter Schwabe
 * Public domain
*/

#include "twistpoint_fp2.h"
#include "curve.h"
#include "fp.h"
#include "fp6.h"

int twistpoint_fp2_zero(twistpoint_fp2_t op)
{
	// twistpoint_fp2_t zero;
	// twistpoint_fp2_init(zero);
	// return !memcmp(zero, op, sizeof(*op));
	return fp2e_iszero(op->m_z);
}

void twistpoint_fp2_init(twistpoint_fp2_t rop)
{
	fp2e_setzero(rop->m_x);
	fp2e_setone(rop->m_y);
	fp2e_setzero(rop->m_z);
  fp2e_setzero(rop->m_t);
}

void twistpoint_fp2_set(twistpoint_fp2_t rop, const twistpoint_fp2_t op)
{
	fp2e_set(rop->m_x, op->m_x);
	fp2e_set(rop->m_y, op->m_y);
	fp2e_set(rop->m_z, op->m_z);
  fp2e_setzero(rop->m_t);
}

void twistpoint_fp2_set_fp2e(twistpoint_fp2_t rop, const fp2e_t x, const fp2e_t y, const fp2e_t z)
{
	fp2e_set(rop->m_x, x);
	fp2e_set(rop->m_y, y);
	fp2e_set(rop->m_z, z);
  fp2e_setzero(rop->m_t);
}

void twistpoint_fp2_affineset_fp2e(twistpoint_fp2_t rop, const fp2e_t x, const fp2e_t y)
{
	fp2e_set(rop->m_x, x);
	fp2e_set(rop->m_y, y);
	fp2e_setone(rop->m_z);
  fp2e_setzero(rop->m_t);
}

void twistpoint_fp2_affineset_mpz(twistpoint_fp2_t rop, const mpz_t x0, const mpz_t x1, const mpz_t y0, const mpz_t y1)
{
	fp2e_set_mpz(rop->m_x, x0, x1);
	fp2e_set_mpz(rop->m_y, y0, y1);
	fp2e_setone(rop->m_z);
	fp2e_setzero(rop->m_t);
}

void twistpoint_fp2_mixadd(twistpoint_fp2_t rop, const twistpoint_fp2_t op1, const twistpoint_fp2_t op2)
{
	fp2e_t tfp2e1, tfp2e2, tfp2e3, tfp2e4, tfp2e5, tfp2e6, tfp2e7, tfp2e8, tfp2e9; // Temporary variables needed for intermediary results
	fp2e_square(tfp2e1, op1->m_z);
	fp2e_mul(tfp2e2, op1->m_z, tfp2e1);
	fp2e_mul(tfp2e3, op2->m_x, tfp2e1);
	fp2e_mul(tfp2e4, op2->m_y, tfp2e2);
	fp2e_sub(tfp2e5, tfp2e3, op1->m_x);
	fp2e_sub(tfp2e6, tfp2e4, op1->m_y);
	fp2e_square(tfp2e7, tfp2e5);
	fp2e_mul(tfp2e8, tfp2e7, tfp2e5);
	fp2e_mul(tfp2e9, op1->m_x, tfp2e7);

	fp2e_double(tfp2e1, tfp2e9);
	fp2e_add(tfp2e1, tfp2e1, tfp2e8);
	fp2e_square(rop->m_x, tfp2e6);
	fp2e_sub(rop->m_x, rop->m_x, tfp2e1);
	fp2e_sub(tfp2e1, tfp2e9, rop->m_x);
	fp2e_mul(tfp2e2, tfp2e1, tfp2e6);
	fp2e_mul(tfp2e3, op1->m_y, tfp2e8);
	fp2e_sub(rop->m_y, tfp2e2, tfp2e3);
	fp2e_mul(rop->m_z, op1->m_z, tfp2e5);
}

void twistpoint_fp2_double(twistpoint_fp2_t rop, const twistpoint_fp2_t op)
{
	fp2e_t tfp2e1, tfp2e2, tfp2e3, tfp2e4; // Temporary variables needed for intermediary results
	fp2e_square(tfp2e1, op->m_y);
	fp2e_mul(tfp2e2, tfp2e1, op->m_x);
	fp2e_double(tfp2e2, tfp2e2);
	fp2e_double(tfp2e2, tfp2e2);
	fp2e_square(tfp2e3, tfp2e1);
	fp2e_double(tfp2e3, tfp2e3);
	fp2e_double(tfp2e3, tfp2e3);
	fp2e_double(tfp2e3, tfp2e3);
	fp2e_square(tfp2e4, op->m_x);
	fp2e_triple(tfp2e4, tfp2e4);
	fp2e_square(rop->m_x, tfp2e4);
	fp2e_double(tfp2e1, tfp2e2);
	fp2e_sub(rop->m_x, rop->m_x, tfp2e1);
	fp2e_sub(tfp2e1, tfp2e2, rop->m_x);
	fp2e_mul(rop->m_z, op->m_y, op->m_z);
	fp2e_double(rop->m_z, rop->m_z);
	fp2e_mul(rop->m_y, tfp2e4, tfp2e1);
	fp2e_sub(rop->m_y, rop->m_y, tfp2e3);
}

void twistpoint_fp2_mul(twistpoint_fp2_t rop, const twistpoint_fp2_t op, const mpz_t scalar)
{
	// TODO: Test...
	size_t i;
	twistpoint_fp2_t r;
	twistpoint_fp2_set(r, op);
	for(i = mpz_sizeinbase(scalar, 2) - 1; i > 0; i--)
	{
		twistpoint_fp2_double(r, r);
		if(mpz_tstbit(scalar, i - 1)) 
			twistpoint_fp2_mixadd(r, r, op);
	}
	twistpoint_fp2_set(rop, r);
}


void twistpoint_fp2_makeaffine(twistpoint_fp2_t op)
{
	fp2e_invert(op->m_z, op->m_z);
	fp2e_mul(op->m_y, op->m_y, op->m_z);
	fp2e_square(op->m_z, op->m_z);
	fp2e_mul(op->m_x, op->m_x, op->m_z);
	fp2e_mul(op->m_y, op->m_y, op->m_z);
	fp2e_setone(op->m_z);
}

void twistpoint_fp2_print(FILE *outfile, const twistpoint_fp2_t op)
{
	fprintf(outfile, "[");
	fp2e_print(outfile, op->m_x);
	fprintf(outfile, ", ");
	fp2e_print(outfile, op->m_y);
	fprintf(outfile, ", ");
	fp2e_print(outfile, op->m_z);
	fprintf(outfile, "]");
}

int twistpoint_fp2_well_formed(const twistpoint_fp2_t op)
{
    if (fp2e_iszero(op->m_x) && fp2e_isone(op->m_y) && fp2e_iszero(op->m_z))
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

	fp2e_t X2, Y2, Z2, X3, Z3, Z6, bcoef, rhs; // Temporary variables needed for intermediary results
	fp2e_square(X2, op->m_x);
	fp2e_square(Y2, op->m_y);
	fp2e_square(Z2, op->m_z);

	fp2e_mul(X3, X2, op->m_x);
	fp2e_mul(Z3, Z2, op->m_z);
	fp2e_square(Z6, Z3);

	fpe_t bfpe;
	fpe_set_mpz(bfpe, b);
	fp2e_invert(bcoef, xi);
	fp2e_mul_fpe(bcoef, bcoef, bfpe);
	fp2e_mul(rhs, Z6, bcoef);
	fp2e_add(rhs, rhs, X3);
	return !memcmp(rhs, Y2, sizeof(*rhs));
}