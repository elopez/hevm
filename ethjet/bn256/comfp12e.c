/*
 * bn256-20080525/comfp12e.c 
 * Peter Schwabe
 * Public domain
*/

#include <stdio.h>
#include <gmp.h>

#include "comfp12e.h"
#include "fp2.h"
#include "fp2e.h"
#include "fp6.h"
#include "fp12.h"


void comfp12e_set(comfp12e_t rop, const comfp12e_t op)
{
	fp2e_set(rop->m_a0, op->m_a0);
	fp2e_set(rop->m_a1, op->m_a1);
	fpe_set(rop->m_a, op->m_a);
}

void comfp12e_setone(comfp12e_t rop)
{
	fp2e_setone(rop->m_a0);
	fp2e_setzero(rop->m_a1);
	fpe_setone(rop->m_a);
}

void comfp12e_set_coeffs(comfp12e_t rop, const fp2e_t a0, const fp2e_t a1, const fpe_t a)
{
	fp2e_set(rop->m_a0, a0);
	fp2e_set(rop->m_a1, a1);
	fpe_set(rop->m_a, a);
}

int comfp12e_isone(const comfp12e_t op)
{
	return(fpe_iszero((op->m_a1)->m_a) && fpe_iszero((op->m_a1)->m_b));
}

// Computes rop *= op
void comfp12e_mul(comfp12e_t rop, const comfp12e_t op1, const comfp12e_t op2)
{
	if(comfp12e_isone(op1))
	{
		comfp12e_set(rop, op2);
		return;
	}
	if(comfp12e_isone(op2))
	{
		comfp12e_set(rop, op1);
		return;
	}

	// PART I 
	fp2e_t s0, s1, s2, t0, t1, t2; // Temporary variables, see "On Compressed Pairings and their Computation"
	fp2e_t tmp1, tmp2, tmp3, tmp4, tmp5, tmp6; // More temporary variables, needed for intermediary results
	fpe_t tfp1, tfp2; // Also temporary variables 

	fp2e_square(tmp1, op1->m_a0);
	fpe_square(tfp1, op1->m_a);
	fp2e_mul_fpe(tmp3, _1o3xi, tfp1);
	fp2e_add(tmp1, tmp1, tmp3);
	fp2e_square(tmp2, op2->m_a0);
	fpe_square(tfp1, op2->m_a);
	fp2e_mul_fpe(tmp3, _1o3xi, tfp1);
	fp2e_add(tmp2, tmp2, tmp3);
	fp2e_mul(tmp3, op1->m_a1, op2->m_a1);
	fp2e_mul(tmp4, op1->m_a0, op2->m_a0);
	fpe_mul(tfp1, op1->m_a, op2->m_a);
    fp2e_mulxi_fpe(tmp5, tfp1);
	fp2e_add(tmp4, tmp4, tmp5);
	fp2e_mul(s0, tmp3, tmp4);
	fp2e_square(s2, tmp3);
	fp2e_mulxi(s2, s2);

	fp2e_mul(tmp4, op1->m_a0, op2->m_a1);
	fp2e_mul(tmp5, op1->m_a1, op2->m_a0);
	fp2e_add(tmp4, tmp4, tmp5);
	fp2e_mulxi(tmp6, tmp3);
	fp2e_mul(s1, tmp6, tmp4);
	fp2e_mul(tmp4, tmp1, tmp2);
	fp2e_add(s1, s1, tmp4);

	fp2e_mul(tmp4, op1->m_a1, tmp2);
	fp2e_mul(tmp5, tmp4, op1->m_a0);
	fp2e_add(s2, s2, tmp5);
	fp2e_mul_fpe(t2, tmp4, op1->m_a);
	fp2e_mul(tmp4, tmp4, op1->m_a1);
	fp2e_add(s0, s0, tmp4);
	fp2e_mul(tmp4, op2->m_a1, tmp1);

	fp2e_mul_fpe(tmp5, tmp4, op2->m_a);
	fp2e_add(t2, t2, tmp5);
	fp2e_mul(tmp5, tmp4, op2->m_a0);
	fp2e_add(s2, s2, tmp5);
	fp2e_mul(tmp4, tmp4, op2->m_a1);
	fp2e_add(s0, s0, tmp4);
	fp2e_mulxi(s0, s0);

	fp2e_mul_fpe(t0, op1->m_a0, op2->m_a);
	fp2e_mul_fpe(tmp4, op2->m_a0, op1->m_a);
	fp2e_add(t0, t0, tmp4);
	fp2e_mul(t0, t0, tmp6);
	fp2e_mul_fpe(t1, op1->m_a1, op2->m_a);
	fp2e_mul_fpe(tmp4, op2->m_a1, op1->m_a);
	fp2e_add(t1, t1, tmp4);
	fp2e_mul(t1, t1, tmp6);

	// PART II:
	fp2e_t u0, u1, u2, t, v0, v1; // And again: Temporary variables, see "On Compressed Pairings and their Computation"

	fp2e_square(tmp1, t0);
	fp2e_square(tmp2, t1);
	fp2e_square(tmp3, t2);
	fp2e_mul(t, tmp1, t0);
	fp2e_mul(tmp4, tmp2, t1);
	fp2e_mulxi(tmp4, tmp4);
	fp2e_add(t, t, tmp4);
	fp2e_mul(tmp4, tmp3, t2);
	fp2e_mul(tmp4, tmp4, xi2);

	fp2e_add(t, t, tmp4);
	fp2e_mul(tmp4, t1, t2);
	fp2e_mulxi(tmp4, tmp4);
	fp2e_sub(u0, tmp1, tmp4);
	fp2e_mul(tmp4, tmp4, t0);
	fp2e_triple(tmp4, tmp4);
	fp2e_sub(t, t, tmp4);
	fp2e_mul(tmp4, t0, t1);
	fp2e_mulxi(u1, tmp3);
	fp2e_sub(u1, u1, tmp4);

	fp2e_mul(tmp4, t0, t2);
	fp2e_sub(u2, tmp2, tmp4);
	fp2e_mul(v0, s0, u0);
	fp2e_mul(tmp1, s1, u2);
	fp2e_mul(tmp2, s2, u1);
	fp2e_add(tmp1, tmp1, tmp2);
	fp2e_mulxi(tmp1, tmp1);
	fp2e_add(v0, v0, tmp1);

	fp2e_mul(v1, s0, u1);
	fp2e_mul(tmp1, s1, u0);
	fp2e_add(v1, v1, tmp1);
	fp2e_mul(tmp1, s2, u2);
	fp2e_mulxi(tmp1, tmp1);
	fp2e_add(v1, v1, tmp1);

	//Negate the coefficient m_a of T:
	fpe_neg(t->m_a, t->m_a);

	fp2e_mul(rop->m_a0, v0, t); // C0 
	fp2e_mul(rop->m_a1, v1, t); // C1

	fpe_square(tfp2, t->m_b);
	fpe_square(tfp1, t->m_a);
#if(ALPHA == 2)
	fpe_double(tfp1, tfp1);
	fpe_sub(rop->m_a, tfp2, tfp1);
#elif(ALPHA == -2)
	fpe_double(tfp1, tfp1);
	fpe_add(rop->m_a, tfp2, tfp1);
#elif(ALPHA == -1)
	fpe_add(rop->m_a, tfp2, tfp1);
#else
#error "ALPHA must be -1, 2 or -2"
#endif
}

// Computes rop *= rop
void comfp12e_square(comfp12e_t rop, const comfp12e_t op)
{
	if(comfp12e_isone(op))
		return;

	fp2e_t s0, s1, s; // Temporary variables, see "On Compressible Pairings and their Computation"
	fp2e_t tmp1, tmp2, tmp3, tmp4, tmp5, tmp6; // Some more temporary variables, needed for intermediary results
	fpe_t tfp1, tfp2; // Even more temporary variables

	fp2e_square(tmp2, op->m_a0);
	fp2e_mul(tmp3, op->m_a0, tmp2);
	fp2e_mul(s0, tmp2, tmp3);
	fpe_square(tfp2, op->m_a);
	fp2e_mul_fpe(tmp5, tmp3, tfp2);
	fp2e_square(tmp6, op->m_a1);
	fp2e_mul(tmp6, tmp6, op->m_a1);
	fp2e_mul(tmp4, tmp2, tmp6);
	fp2e_sub(tmp5, tmp5, tmp4);
	fp2e_mulxi(tmp1, tmp5);
	fp2e_double(tmp1, tmp1);
	fp2e_add(s1, s0, tmp1);
	
	fp2e_sub(tmp5, tmp5, tmp4);
	fp2e_mulxi(tmp5, tmp5);
	fp2e_add(s0, s0, tmp5);
	fpe_square(tfp1, tfp2);
	fp2e_mul_fpe(tmp5, op->m_a0, tfp1);
	fp2e_mul_fpe(tmp1, tmp5, _1o3modp);

	fp2e_mul_fpe(tmp2, tmp6, tfp2);
	fp2e_sub(tmp1, tmp1, tmp2);
	fp2e_double(tmp2, tmp2);
	fp2e_sub(tmp5, tmp5, tmp2);
	fp2e_mul(tmp1, tmp1, xi2);
	fp2e_mul(tmp5, tmp5, xi2);
	fp2e_add(s0, s0, tmp1);
	fp2e_mul(s0, s0, op->m_a0);
	fp2e_add(s1, s1, tmp5);
	fp2e_mul(s1, s1, op->m_a0);
	fp2e_square(tmp3, tmp6);
	fp2e_mul(tmp3, tmp3, xi2);
	fpe_mul(tfp2, tfp1, tfp2);
	fp2e_mul_fpe(tmp5, _1o27xi3, tfp2);
	fp2e_add(tmp2, tmp3, tmp5);
	fp2e_add(s0, s0, tmp2);
	fp2e_mul_fpe(s, s0, op->m_a);
	fp2e_mul(s0, s0, op->m_a0);
	fp2e_double(s, s);
	fp2e_double(tmp5, tmp5);
	fp2e_double(tmp5, tmp5);
	fp2e_add(tmp2, tmp3, tmp5);
	fp2e_add(s1, s1, tmp2);
	fp2e_mul(s1, s1, op->m_a1);
	fp2e_mul_fpe(tmp2, tmp6, tfp1);
	fp2e_mul(tmp2, tmp2, _1o3xi3);
	fp2e_sub(s0, s0, tmp2);
	
	//Negate the coefficient m_a of T:
	fpe_neg(s->m_a, s->m_a);

	fp2e_mul(rop->m_a0, s0, s); // C0 
	fp2e_mul(rop->m_a1, s1, s); // C1

	fpe_square(tfp2, s->m_b);
	fpe_square(tfp1, s->m_a);
#if(ALPHA == 2)
	fpe_double(tfp1, tfp1);
	fpe_sub(rop->m_a, tfp2, tfp1);
#elif(ALPHA == -2)
	fpe_double(tfp1, tfp1);
	fpe_add(rop->m_a, tfp2, tfp1);
#elif(ALPHA == -1)
	fpe_add(rop->m_a, tfp2, tfp1);
#else
#error "ALPHA must be -1, 2 or -2"
#endif
}

void comfp12e_pow(comfp12e_t rop, const comfp12e_t op, const mpz_t exp)
{
	// TODO: Implement Square-and-multiply generically using function pointers
	comfp12e_t dummy;
	comfp12e_set(dummy, op);
	comfp12e_set(rop, op);
	int i;
	for(i = mpz_sizeinbase(exp, 2) - 1; i > 0; i--)
	{
		comfp12e_square(rop, rop);
		if(mpz_tstbit(exp, i - 1)) 
			comfp12e_mul(rop, rop, dummy);
	}
}


void comfp12e_invert(comfp12e_t rop, const comfp12e_t op)
{
	comfp12e_set(rop, op);
	fpe_neg(rop->m_a, rop->m_a);
}

void comfp12e_makeaffine(comfp12e_t rop, const comfp12e_t op)
{
	fpe_t tmp;
	fpe_invert(tmp, op->m_a);
	fp2e_mul_fpe(rop->m_a0, op->m_a0, tmp);
	fp2e_mul_fpe(rop->m_a1, op->m_a1, tmp);
	fpe_setone(rop->m_a);
}

void comfp12e_frobenius_p(comfp12e_t rop, const comfp12e_t op)
{
	fp2e_t tmp;

	comfp12e_set(rop, op);
	fpe_neg((rop->m_a0)->m_a, (rop->m_a0)->m_a);
	fpe_neg((rop->m_a1)->m_a, (rop->m_a1)->m_a);

	fp2e_mul(rop->m_a1, rop->m_a1, zpminus1inv);
	fp2e_square(tmp, zpminus1inv);
	fp2e_mul(tmp, tmp, zpminus1inv);
	fp2e_mul(rop->m_a0, rop->m_a0, tmp);
}

void comfp12e_frobenius_p2(comfp12e_t rop, const comfp12e_t op)
{
	fpe_t tmp;

	fpe_square(tmp, zeta);
	fp2e_set(rop->m_a0, op->m_a0);
	fp2e_mul_fpe(rop->m_a1, op->m_a1, tmp);
	fpe_neg(rop->m_a, op->m_a);
}

void comfp12e_print(FILE *outfile, const comfp12e_t op)
{
	fprintf(outfile, "(");
	fp2e_print(outfile, op->m_a0);
	fprintf(outfile, ", ");
	fp2e_print(outfile, op->m_a1);
	fprintf(outfile, ", ");
	fpe_print(outfile, op->m_a);
	fprintf(outfile, ")");
}
