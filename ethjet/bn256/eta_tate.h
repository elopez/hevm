/*
 * bn256-20080525/eta_tate.h 
 * Peter Schwabe
 * Public domain
*/

#ifndef ETA_TATE_H
#define ETA_TATE_H

#include "curvepoint_fp.h"
#include "twistpoint_fp2.h"
#include "fp12e.h"

void eta(fp12e_t rop, const curvepoint_fp_t op1, const twistpoint_fp2_t op2);
void tate(fp12e_t rop, const curvepoint_fp_t op1, const twistpoint_fp2_t op2);

#endif // ifdef ETA_TATE_H
