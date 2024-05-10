/*
 * bn256-20080525/final_expo.h 
 * Peter Schwabe
 * Public domain
*/

#ifndef FINAL_EXPO_H
#define FINAL_EXPO_H

#include "fp12e.h"
#include "comfp12e.h"

void final_expo(fp12e_t rop);
void final_expo_com(comfp12e_t rop);

#endif // ifdef FINAL_EXPO_H
