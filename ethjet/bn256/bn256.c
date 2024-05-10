/*
 * bn256-20080525/bnpairings.c 
 * Peter Schwabe
 * Public domain
 */

#include <stdio.h>

#include "init.h"
#include "curvepoint_fp.h"
#include "twistpoint_fp2.h"
#include "points.h"
#include "fp12e.h"
#include "comfp12e.h"
#include "randompoints.h"

#include "ate_optate.h"
#include "eta_tate.h"
#include "cometa_comtate.h"
#include "cpucycles.h"

#ifdef BENCH
#include "fp.h"
#include "fp2.h"
#include "fp12.h"
#endif

static void print_help(void)
{
  printf("t: Tate Pairing\n");
  printf("e: Generalized Eta Pairing\n");
  printf("a: Generalized Ate Pairing\n");
  printf("o: Optimal Ate Pairing\n");
  printf("u: Compressed generalized Tate Pairing\n");
  printf("f: Compressed generalized Eta Pairing\n");
  printf("n: Generate two new random points\n");
  printf("s: Set P and Q to predefined generators\n");
  printf("q: Quit the program");
}

#ifdef BENCH
static void set_ctr_zero(void)
{
      nummultp=0; multpcycles = 0;
      nummultp2=0; multp2cycles = 0;
      nummultp12=0; multp12cycles = 0;
      numsqp=0; sqpcycles = 0;
      numsqp2=0; sqp2cycles = 0;
      numsqp12=0; sqp12cycles = 0;
      numinvp=0; invpcycles = 0;
      numinvp2=0; invp2cycles = 0;
      numinvp12=0; invp12cycles = 0;
}

static void print_cycles(void)
{
  if(nummultp > 0)
    printf("Number of multiplications F_p: %llu, cycles per multiplication: %llu\n", nummultp, multpcycles/nummultp);
  if(numsqp > 0)
    printf("Number of squarings in F_p: %llu, cycles per squaring: %llu\n", numsqp, sqpcycles/numsqp);
  if(numinvp > 0)
    printf("Number of inversions in F_p: %llu, cycles per inversion: %llu\n", numinvp, invpcycles/numinvp);
  if(nummultp2 > 0)
    printf("Number of multiplications F_p^2: %llu, cycles per multiplication: %llu\n", nummultp2, multp2cycles/nummultp2);
  if(numsqp2 > 0)
    printf("Number of squarings in F_p^2: %llu, cycles per squaring: %llu\n", numsqp2, sqp2cycles/numsqp2);
  if(numinvp2 > 0)
    printf("Number of inversions in F_p^2: %llu, cycles per inversion: %llu\n", numinvp2, invp2cycles/numinvp2);
  if(nummultp12 > 0)
    printf("Number of multiplications F_p^12: %llu, cycles per multiplication: %llu\n", nummultp12, multp12cycles/nummultp12);
  if(numsqp12 > 0)
    printf("Number of squarings in F_p^12: %llu, cycles per squaring: %llu\n", numsqp12, sqp12cycles/numsqp12);
  if(numinvp12 > 0)
    printf("Number of inversions in F_p^12: %llu, cycles per inversion: %llu\n", numinvp12, invp12cycles/numinvp12);
}
#endif


int main(int argc, char* argv[])
{
  init_globals();

  curvepoint_fp_t point1;
  twistpoint_fp2_t point2;

  long long int start, stop;

  get_random_curvepoint_fp(point1);
  get_random_twistpoint_fp2(point2);

  printf("Pairing computation of the points\nP =\n");
  curvepoint_fp_makeaffine(point1);
  curvepoint_fp_print(stdout, point1);
  printf("\n\nand Q =\n");
  twistpoint_fp2_makeaffine(point2);
  twistpoint_fp2_print(stdout, point2);
  printf("\n\n");

  int c1, c2, ctr;
  char choice;
  for(;;)	
  {
    do
    {
      ctr = 0;
      printf("\nCommand (h for help)> ");
      c1 = getchar();
      if(c1 != '\n' && c1 != EOF)
      {
        ++ctr;
        while((c2 = getchar()) != EOF && c2 != '\n')
          ++ctr;
        if(ctr != 1)
          printf("Unknown command, press 'h' for a list of available commands");
      }
      choice = (char) c1;
    }
    while(ctr != 1);

    if (choice == 't')
    {	
#ifdef BENCH
      set_ctr_zero();
#endif
      fp12e_t rop;
      start = cpucycles();
      tate(rop, point1, point2);
      stop = cpucycles();
      fp12e_print(stdout, rop);
      printf("\n\nCycles needed for computation: %llu\n", stop - start);
      printf("This corresponds to %f seconds\n", (double) (stop - start) / cpucycles_persecond());
#ifdef BENCH
      print_cycles();
#endif
    }
    else if (choice == 'e')
    {	
#ifdef BENCH
      set_ctr_zero();
#endif
      fp12e_t rop;
      start = cpucycles();
      eta(rop, point1, point2);
      stop = cpucycles();
      fp12e_print(stdout, rop);
      printf("\n\nCycles needed for computation: %llu\n", stop - start);
      printf("This corresponds to %f seconds\n", (double) (stop - start) / cpucycles_persecond());
#ifdef BENCH 
      print_cycles();
#endif
    }
    else if (choice == 'a')
    {	
#ifdef BENCH
      set_ctr_zero();
#endif
      fp12e_t rop;
      start = cpucycles();
      ate(rop, point2, point1);
      stop = cpucycles();
      fp12e_print(stdout, rop);
      printf("\n\nCycles needed for computation: %llu\n", stop - start);
      printf("This corresponds to %f seconds\n", (double) (stop - start) / cpucycles_persecond());
#ifdef BENCH 
      print_cycles();
#endif
    }
    else if (choice == 'o')
    {	
#ifdef BENCH 
      set_ctr_zero();
#endif
      fp12e_t rop;
      start = cpucycles();
      optate(rop, point2, point1);
      stop = cpucycles();
      fp12e_print(stdout, rop);
      printf("\n\nCycles needed for computation: %llu\n", stop - start);
      printf("This corresponds to %f seconds\n", (double) (stop - start) / cpucycles_persecond());
#ifdef BENCH 
      print_cycles();
#endif
    }
    else if (choice == 'u')
    {	
#ifdef BENCH 
      set_ctr_zero();
#endif
      comfp12e_t rop;
      start = cpucycles();
      comtate(rop, point1, point2);
      comfp12e_makeaffine(rop, rop);
      stop = cpucycles();
      comfp12e_print(stdout, rop);
      printf("\n\nCycles needed for computation: %llu\n", stop - start);
      printf("This corresponds to %f seconds\n", (double) (stop - start) / cpucycles_persecond());
#ifdef BENCH 
      print_cycles();
#endif
    }
    else if (choice == 'f')
    {	
#ifdef BENCH 
      set_ctr_zero();
#endif
      comfp12e_t rop;
      start = cpucycles();
      cometa(rop, point1, point2);
#ifndef BENCH 
      comfp12e_makeaffine(rop, rop);
#endif
      stop = cpucycles();
      comfp12e_print(stdout, rop);
      printf("\n\nCycles needed for computation: %llu\n", stop - start);
      printf("This corresponds to %f seconds\n", (double) (stop - start) / cpucycles_persecond());
#ifdef BENCH 
      print_cycles();
#endif
    }
    else if (choice == 'n')
    {
      get_random_curvepoint_fp(point1);
      get_random_twistpoint_fp2(point2);
      printf("The new points are\nP =\n");
      curvepoint_fp_print(stdout, point1);
      printf("\n\nand Q =\n");
      twistpoint_fp2_print(stdout, point2);
      printf("\n\n");

    }
    else if (choice == 's')
    {
      curvepoint_fp_set(point1, curve_gen);
      twistpoint_fp2_set(point2, twist_gen);
      printf("The new points are\nP =\n");
      curvepoint_fp_print(stdout, point1);
      printf("\n\nand Q =\n");
      twistpoint_fp2_print(stdout, point2);
      printf("\n\n");

    }
    else if (choice == 'h')
    {
      print_help();
    }
    else if (choice == 'q')
    {
      clear_globals();
      printf("Exiting gracefully...\n");
      return 0;
    }
    else if (choice != '\n')
    {
      printf("Unknown command, press 'h' for a list of available commands");
    }
  }
}


