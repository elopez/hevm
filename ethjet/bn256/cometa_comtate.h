/*
 * bn256-20080525/cometa_comtate.h 
 * Peter Schwabe
 * Public domain
*/

#ifndef COMETA_COMTATE_H
#define COMETA_COMTATE_H

#include "curvepoint_fp.h"
#include "twistpoint_fp2.h"
#include "comfp12e.h"


void cometa(comfp12e_t rop, const curvepoint_fp_t op1, const twistpoint_fp2_t op2);
void comtate(comfp12e_t rop, const curvepoint_fp_t op1, const twistpoint_fp2_t op2);

#endif // ifdef COMETA_COMTATE_H
