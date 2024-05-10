/*
 * bn256-20080525/randompoints.h 
 * Peter Schwabe
 * Public domain
*/

#ifndef RANDOMPOINTS_H
#define RANDOMPOINTS_H

#include "curvepoint_fp.h"
#include "twistpoint_fp2.h"

void get_random_curvepoint_fp(curvepoint_fp_t rop);

void get_random_twistpoint_fp2(twistpoint_fp2_t rop);

#endif // ifdef RANDOMPOINTS_H
