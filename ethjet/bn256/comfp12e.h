/*
 * bn256-20080525/comfp12e.h 
 * Peter Schwabe
 * Public domain
*/

#ifndef COMFP12E_H
#define COMFP12E_H

#include <gmp.h>
#include <stdio.h>

#include "fp12e.h"
#include "fp2e.h"

typedef struct comfp12e_struct comfp12e_struct_t;

struct comfp12e_struct
{
	fp2e_t m_a0;
	fp2e_t m_a1;
	fpe_t m_a;
};

typedef comfp12e_struct_t comfp12e_t[1];


void comfp12e_set(comfp12e_t rop, const comfp12e_t op);

void comfp12e_setone(comfp12e_t rop);

void comfp12e_set(comfp12e_t rop, const comfp12e_t op);

void comfp12e_set_coeffs(comfp12e_t rop, const fp2e_t a0, const fp2e_t a1, const fpe_t a);

int comfp12e_isone(const comfp12e_t op);

// Computes rop *= op
void comfp12e_mul(comfp12e_t rop, const comfp12e_t op1, const comfp12e_t op2);

// Computes rop *= rop
void comfp12e_square(comfp12e_t rop, const comfp12e_t op);

void comfp12e_pow(comfp12e_t rop, const comfp12e_t op, const mpz_t exp);

void comfp12e_invert(comfp12e_t rop, const comfp12e_t op);

void comfp12e_makeaffine(comfp12e_t rop, const comfp12e_t op);

void comfp12e_frobenius_p(comfp12e_t rop, const comfp12e_t op);

void comfp12e_frobenius_p2(comfp12e_t rop, const comfp12e_t op);

void comfp12e_print(FILE *outfile, const comfp12e_t op);

#endif // ifdef COMFP12E_H
