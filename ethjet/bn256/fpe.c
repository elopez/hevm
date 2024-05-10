/*
 * bn256-20080525/fpe.c 
 * Peter Schwabe
 * Public domain
*/

#include <gmp.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
//#include <cpucycles.h>
#include <math.h>

#include "fpe.h"
#include "fp.h"

extern unsigned long long multpcycles; extern unsigned long long nummultp;
extern unsigned long long invpcycles; extern unsigned long long numinvp;
extern unsigned long long multp2cycles; extern unsigned long long nummultp2;
extern unsigned long long invp2cycles; extern unsigned long long numinvp2;
extern unsigned long long multp12cycles; extern unsigned long long nummultp12;
extern unsigned long long invp12cycles; extern unsigned long long numinvp12;


////////////////////////////////////////////////////////////////
//     "Private functions"    //////////////////////////////////
////////////////////////////////////////////////////////////////

// Montgomery reduction, see HAC, alg 14.32, size of op must be 2 * N_LIMBS:
static void fpe_montgomery_reduce(fpe_t rop, mp_limb_t *op)
{
	mp_limb_t u, t, c, d; // dummies
	int i;
	c = 0;

	for(i = 0; i < N_LIMBS; i++)
	{
		u = op[i] * p_inv;
		d = mpn_addmul_1(op + i, p->_mp_d, N_LIMBS, u);
		t = op[N_LIMBS + i] + d;
		op[N_LIMBS + i] = t + c;
		if(t < d || (c == 1 && t + c == 0))
			c = 1;
		else
			c = 0;
	}
	if(c || mpn_cmp(op + N_LIMBS, p->_mp_d, N_LIMBS) >= 0)
		// Result is larger than p, subtract p from the result:
		mpn_sub_n(rop->m_value, op + N_LIMBS, p->_mp_d, N_LIMBS);
	else
	{
		for(i = 0; i < N_LIMBS; i++)
			rop->m_value[i] = op[N_LIMBS + i];
	}
}

// Montgomery transformation:
static void fpe_montgomery_trans(fpe_t rop)
{
	mp_limb_t tmp1[2 * N_LIMBS], tmp2[2 * N_LIMBS]; // needed for intermediate results
	int i;

	for(i = 0; i < N_LIMBS; i++)
		tmp1[i] = 0;
	for(i = 0; i < N_LIMBS; i++)
		tmp1[N_LIMBS + i] = rop->m_value[i];
	mpn_tdiv_qr(tmp2, rop->m_value, 0, tmp1, 2 * (N_LIMBS), p->_mp_d, N_LIMBS);
}

// Montgomery retransformation:
static void fpe_montgomery_retrans(fpe_t rop, const fpe_t op)
{
	mp_limb_t tmp[2 * N_LIMBS]; // needed for intermediate results
	int i;

	for(i = 0; i < N_LIMBS; i++)
		tmp[N_LIMBS + i] = 0;
	for(i = 0; i < N_LIMBS; i++)
		tmp[i] = op->m_value[i];
	fpe_montgomery_reduce(rop, tmp);
}


////////////////////////////////////////////////////////////////
//     "public functinos"    /////////////////////////////////
///////////////////////////////////////////////////////////////


// Set fpe_t rop to given value:
void fpe_set(fpe_t rop, const fpe_t op)
{
	int i;
	for(i = 0; i < N_LIMBS; i++)
		rop->m_value[i] = op->m_value[i];
}

// Set fpe_t rop to given value:
void fpe_set_ui(fpe_t rop, const unsigned int op)
{
	int i;
	for(i = 1; i < N_LIMBS; i++)
		rop->m_value[i] = 0;
	rop->m_value[0] = op;

	fpe_montgomery_trans(rop);
}

void fpe_set_bytearray(fpe_t rop, const unsigned char *op, size_t oplen)
{
	int i, j = 0;
	
	for(i = 0; i < N_LIMBS; i++)
		rop->m_value[i] = 0;

	size_t gmp_limb_bytes = GMP_LIMB_BITS/8;
	assert(oplen <= gmp_limb_bytes * N_LIMBS);
	for(i = 0; i < N_LIMBS; i++)
	{
		for(j = 0; j < gmp_limb_bytes && oplen > i * gmp_limb_bytes + j; j++)
		{
			rop->m_value[i] += (mp_limb_t) op[i * gmp_limb_bytes + j] << 8*j;
		}
	}
}

// Set fpe_t rop to given value given as mpz
void fpe_set_mpz(fpe_t rop, const mpz_t op)
{
	for(int i = 0; i < N_LIMBS; i++)
	{
		rop->m_value[i] = mpz_getlimbn(op, i);
	}
	fpe_montgomery_trans(rop);
}

// Set mpz_t rop to given value given as fpe_t
void fpe_get_mpz(fpe_t op, mpz_t rop)
{
	fpe_t dummy;
	fpe_montgomery_retrans(dummy, op);

	mp_limb_t* xp = mpz_limbs_modify(rop, N_LIMBS);
	for(int i = 0; i < N_LIMBS; i++)
	{
		xp[i] = dummy->m_value[i];
	}
	mpz_limbs_finish(rop, N_LIMBS);
}

// Set fpe_t rop to given value given as (ASCII) string
void fpe_set_str(fpe_t rop, const char* op)
{

	// Initialize all limbs with 0:
	int i;
	for(i = 0; i < N_LIMBS; i++)
	{
		rop->m_value[i] = 0;
	}
	// Determine the length of op:
	const char *scan = op;
	int size = 0;
	while(*scan != 0)
	{
		++scan;
		++size;
	}

	unsigned char str[size];
	memcpy(str, op, size);

	// Convert from ASCII:
	for(i = 0; i < size; ++i)
		str[i] -= '0';

	mpn_set_str(rop->m_value, str, size, 10);

	fpe_montgomery_trans(rop);
}

// Set rop to one
void fpe_setone(fpe_t rop)
{
	int i;
	for(i = 1; i < N_LIMBS; i++)
		rop->m_value[i] = 0;
	rop->m_value[0] = 1;
	fpe_montgomery_trans(rop);
}

// Set rop to zero
void fpe_setzero(fpe_t rop)
{
	int i;
	for(i = 0; i < N_LIMBS; i++)
		rop->m_value[i] = 0;
}

// Return 1 if op is zero, 0 otherwise
int fpe_iszero(const fpe_t op)
{
	int i;
	for(i = 0; i < N_LIMBS; i++)
		if(op->m_value[i]) return 0;
	return 1;
}

// Return 1 if op is one, 0 otherwise
int fpe_isone(const fpe_t op)
{
	fpe_t dummy;
	fpe_montgomery_retrans(dummy, op);

	int i;
	for(i = 1; i < N_LIMBS; i++)
		if(dummy->m_value[i]) return 0;
	if(dummy->m_value[0] != 1) return 0;
	return 1;
}

// Compute the negative of an fpe
void fpe_neg(fpe_t rop, const fpe_t op)
{
	mpn_sub_n(rop->m_value, p->_mp_d, op->m_value, N_LIMBS);
}

// Double an fpe:
void fpe_double(fpe_t rop, const fpe_t op)
{
	mp_limb_t c;
	c = mpn_lshift(rop->m_value, op->m_value, N_LIMBS, 1);
	
	// Reduce if result is larger than p: 
	if(c || mpn_cmp(rop->m_value, p->_mp_d, N_LIMBS) > 0)
		mpn_sub_n(rop->m_value, rop->m_value, p->_mp_d, N_LIMBS);
}

// Halve an fpe:
void fpe_halve(fpe_t rop, const fpe_t op)
{
	if((op->m_value[0] % 2) == 0)
		mpn_rshift(rop->m_value, op->m_value, N_LIMBS, 1);
	else
	{
		int c;
		c = (mpn_add_n(rop->m_value, op->m_value, p->_mp_d, N_LIMBS)) << (GMP_LIMB_BITS - 1);
		mpn_rshift(rop->m_value, rop->m_value, N_LIMBS, 1);
		rop->m_value[N_LIMBS - 1] ^= c;
	}
}

// Triple an fpe:
void fpe_triple(fpe_t rop, const fpe_t op)
{
	mp_limb_t tmp[2 * N_LIMBS]; // needed for intermediate results
	int i, c;

	for(i = 2 * N_LIMBS - 1; i >= N_LIMBS; i--)
		tmp[i] = 0;

	// Double
	c = mpn_lshift(tmp, op->m_value, N_LIMBS, 1);
	if(c || mpn_cmp(tmp, p->_mp_d, N_LIMBS) >= 0)
		mpn_sub_n(tmp, tmp, p->_mp_d, N_LIMBS);

	// Add
	c = mpn_add_n(rop->m_value, tmp, op->m_value, N_LIMBS);	
	if(c || mpn_cmp(rop->m_value, p->_mp_d, N_LIMBS) >= 0)
		mpn_sub_n(rop->m_value, rop->m_value, p->_mp_d, N_LIMBS);
}

// Add two fpe, store result in rop:
void fpe_add(fpe_t rop, const fpe_t op1, const fpe_t op2)
{
	mp_limb_t c;
	c = mpn_add_n(rop->m_value, op1->m_value, op2->m_value, N_LIMBS);

	// Reduce if result is larger than p: 
	if(c || mpn_cmp(rop->m_value, p->_mp_d, N_LIMBS) >= 0)
		mpn_sub_n(rop->m_value, rop->m_value, p->_mp_d, N_LIMBS);
}

// Subtract op2 from op1, store result in rop:
void fpe_sub(fpe_t rop, const fpe_t op1, const fpe_t op2)
{
	mp_limb_t b;
	b = mpn_sub_n(rop->m_value, op1->m_value, op2->m_value, N_LIMBS);
	
	if(b)
		mpn_add_n(rop->m_value, rop->m_value, p->_mp_d, N_LIMBS);
}

// Multiply two fpe, store result in rop:
void fpe_mul(fpe_t rop, const fpe_t op1, const fpe_t op2)
{
#ifdef BENCH
    nummultp++;
    multpcycles -= cpucycles();
#endif

	mp_limb_t tmp[2 * N_LIMBS]; // needed for intermediate results


	if(fpe_iszero(op1) || fpe_iszero(op2))
		fpe_setzero(rop);
	else
	{
		mpn_mul_n(tmp, op1->m_value, op2->m_value, N_LIMBS);
		fpe_montgomery_reduce(rop, tmp);
	}

#ifdef BENCH
    multpcycles += cpucycles();
#endif
}



/*
// Square an fpe, store result in rop:
void fpe_square(fpe_t rop, const fpe_t op)
{
#ifdef BENCH
    numsqp++;
    sqpcycles -= cpucycles();
#endif
	fpe_mul(rop, op, op);
#ifdef BENCH
    sqpcycles += cpucycles();
#endif
}
*/

// Compute inverse of an fpe, store result in rop:
void fpe_invert(fpe_t rop, const fpe_t op)
{
#ifdef BENCH
    numinvp++;
    invpcycles -= cpucycles();
#endif
	/*
	 * FIXME: This code doesn't work they way it should
	// Using mpn_gcdext:
	mp_limb_t mp_limb_t d1[N_LIMBS + 1], mp_limb_t d2[N_LIMBS + 1], tmp1[2 * N_LIMBS], tmp2[2 * N_LIMBS]; // needed for intermediate results
	mp_size_t n1, n2;
	
	int i;
	for(i = 0; i < N_LIMBS; i++)
	{
		d1[i] = p->_mp_d[i];
		d2[i] = op->m_value[i];
	}

	n2 = mpn_gcdext(tmp1, tmp2, &n1, d2, N_LIMBS, d1, N_LIMBS);
	

	for(i = 0; i < N_LIMBS; i++)
		rop->m_value[i] = tmp2[i];

	fpe_montgomery_trans(rop);
	fpe_montgomery_trans(rop);
	
	*/


	// Partial Montgomery inversion, see Guide to ECC, alg. 2.23
	mp_limb_t u[N_LIMBS], v[N_LIMBS], s[N_LIMBS];
	unsigned long i, k, c;
	c = k = 0; 
	for(i = 0; i < N_LIMBS; i++)
	{
		u[i] = p->_mp_d[i];
		v[i] = op->m_value[i];
		rop->m_value[i] = 0;
		s[i] = 0;
	}
	s[0] = 1;

	int v_zero = 1;
	for(i = 0; i < N_LIMBS && v_zero; i++)
		if(v[i]) v_zero = 0;
	// Don't try to invert 0:
	assert(!v_zero);

	while(!v_zero)
	{
		if(u[0] % 2 == 0)
		{
			mpn_rshift(u, u, N_LIMBS, 1);
			mpn_lshift(s, s, N_LIMBS, 1);
		}
		else
		{
			if(v[0] % 2 == 0)
			{
				mpn_rshift(v, v, N_LIMBS, 1);
				mpn_lshift(rop->m_value, rop->m_value, N_LIMBS, 1);
			}
			else
			{
				if(mpn_cmp(u, v, N_LIMBS) > 0)
				{
					mpn_sub_n(u, u, v, N_LIMBS);
					mpn_rshift(u, u, N_LIMBS, 1);
					c = mpn_add_n(rop->m_value, rop->m_value, s, N_LIMBS);
					mpn_lshift(s, s, N_LIMBS, 1);
				}
				else
				{
					mpn_sub_n(v, v, u, N_LIMBS);
					mpn_rshift(v, v, N_LIMBS, 1);
					mpn_add_n(s, s, rop->m_value, N_LIMBS);
					c = mpn_lshift(rop->m_value, rop->m_value, N_LIMBS, 1);
				}
			}
		}
		++k;
		// Check, whether v is zero:
		v_zero = 1;
		for(i = 0; i < N_LIMBS; i++)
		{
			if(v[i]) v_zero = 0;
		}
	}
	
	if(c || mpn_cmp(rop->m_value, p->_mp_d, N_LIMBS) >= 0)
		mpn_sub_n(rop->m_value, rop->m_value, p->_mp_d, N_LIMBS);

	mpn_sub_n(rop->m_value, p->_mp_d, rop->m_value, N_LIMBS);
	
	// Make the Montgomery Inversion complete:
	for(; k < N_LIMBS * sizeof(mp_limb_t) * 16; k++)
	{
		c = mpn_lshift(rop->m_value, rop->m_value, N_LIMBS, 1);
		if(c || mpn_cmp(rop->m_value, p->_mp_d, N_LIMBS) >= 0)
			mpn_sub_n(rop->m_value, rop->m_value, p->_mp_d, N_LIMBS);
	}
#ifdef BENCH
    invpcycles += cpucycles();
#endif
}

void fpe_sqrt(fpe_t rop, const fpe_t op)
{
	
}

// Print the element to stdout:
void fpe_print(FILE *outfile, const fpe_t op)
{
	int i;

	// Retransform from Montgomery domain:
	fpe_t dummy;
	fpe_montgomery_retrans(dummy, op);

	i = 0;
	
	while(i < N_LIMBS && !dummy->m_value[N_LIMBS - 1 - i])
		i++;

	// Print '0' if op is zero:
	if(i == N_LIMBS && !dummy->m_value[0])
	{
		fputc('0', outfile);
		fputc(0, outfile);
	}
	else
	{
		unsigned char *str;
		double log2 = M_LOG10E * M_LN2; // log_{10}2
		
		size_t str_size = log2 * (GMP_LIMB_BITS * (N_LIMBS)) + 2;
		str = malloc(str_size);

		str_size = mpn_get_str(str, 10, dummy->m_value, N_LIMBS - i);
	
		// Strip leading zeros:
		while(*str == 0 && str_size > 1)
		{
			str++;
			str_size--;
		}
		
		for(i = 0; i < str_size; i++)
		{
			fputc(*(str + i) + '0', outfile);
		}
		fputc(0, outfile);
		free(str);
	}
}

// Convert fpe into a bytearray, the destination must have space for N_LIMBS mpz_limb_t
void fpe_to_bytearray(unsigned char * rop, const fpe_t op)
{
	int i,j;
	mp_limb_t dummy;
	size_t gmp_limb_bytes = GMP_LIMB_BITS/8;
	for(i = 0; i < N_LIMBS; i++)
	{
		dummy = op->m_value[i];
		for(j = 0; j < gmp_limb_bytes; j++)
		{
			rop[gmp_limb_bytes * i + j] = (unsigned char) (dummy & 0xff);
			dummy >>= 8;
		}
	}
}
