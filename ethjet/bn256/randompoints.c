/*
 * bn256-20080525/randompoints.c 
 * Peter Schwabe
 * Public domain
*/

#include <gmp.h>
#include <time.h>

#include "randompoints.h"
#include "curvepoint_fp.h"
#include "twistpoint_fp2.h"
#include "points.h"

void get_random_curvepoint_fp(curvepoint_fp_t rop)
{
	gmp_randstate_t state;
	gmp_randinit_default(state);
	unsigned long int seed;
	seed = (unsigned long int)time(NULL);
	gmp_randseed_ui(state, seed);

	mpz_t scalar;
	mpz_init(scalar);
	
	do
		mpz_urandomm(scalar, state, p);
	while(!mpz_cmp_ui(scalar, 0));

	curvepoint_fp_mul(rop, curve_gen, scalar);

	mpz_clear(scalar);
	gmp_randclear(state);
}

void get_random_twistpoint_fp2(twistpoint_fp2_t rop)
{
	gmp_randstate_t state;
	gmp_randinit_default(state);
	unsigned long int seed;
	seed = (unsigned long int)time(NULL);
	gmp_randseed_ui(state, seed);

	mpz_t scalar;
	mpz_init(scalar);
	
	do
		mpz_urandomm(scalar, state, p);
	while(!mpz_cmp_ui(scalar, 0));

	twistpoint_fp2_mul(rop, twist_gen, scalar);

	mpz_clear(scalar);
	gmp_randclear(state);
}
