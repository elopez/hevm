#include "ethjet.h"
#include "ethjet-ff.h"
#include "tinykeccak.h"
#include "blake2.h"

#include "bn256/init.h"
#include "bn256/curvepoint_fp.h"
#include "bn256/twistpoint_fp2.h"

#include <gmp.h>
#include <secp256k1_recovery.h>

#include <stdlib.h>
#include <stdint.h>
#include <string.h>

struct ethjet_context *
ethjet_init ()
{
  struct ethjet_context *ctx;
  ctx = malloc (sizeof *ctx);
  if (!ctx) return NULL;

  ctx->ec = secp256k1_context_create (SECP256K1_CONTEXT_VERIFY);

  init_globals();

  return ctx;
}

void
ethjet_free (struct ethjet_context *ctx)
{
  secp256k1_context_destroy (ctx->ec);
  free (ctx);
}

/* 
 * The example contract at 0xdeadbeef just reverses its input.
 */
int
ethjet_example (struct ethjet_context *ctx,
                uint8_t *in, size_t in_size,
                uint8_t *out, size_t out_size)
{
  if (out_size != in_size)
    return 0;

  for (int i = 0; i < in_size; i++)
    out[i] = in[in_size - i - 1];

  return 1;
}

int
ethjet_ecrecover (secp256k1_context *ctx,
                  uint8_t *in, size_t in_size,
                  uint8_t *out, size_t out_size)
{
  /* Input: H V R S, all 32 bytes. */

  secp256k1_pubkey pubkey;
  secp256k1_ecdsa_recoverable_signature rsig;

  uint8_t *input64;
  uint8_t pubkey_hex[65];
  size_t hexlen = 65;

  int recid;

  if (in_size != 128)
    return 0;

  if (out_size != 32)
    return 0;

  input64 = in + 64;
  recid = in[63] - 27;

  /* higher bytes of V should be zero  */
  static const char z31 [31];
  if (memcmp (z31, in + 32, 31))
    return 0;

  if (recid < 0 || recid > 3)
    return 0;

  if (!secp256k1_ecdsa_recoverable_signature_parse_compact
      (ctx, &rsig, input64, recid))
    return 0;

  if (!secp256k1_ecdsa_recover (ctx, &pubkey, &rsig, in))
    return 0;

  if (!secp256k1_ec_pubkey_serialize
      (ctx, pubkey_hex, &hexlen, &pubkey, SECP256K1_EC_UNCOMPRESSED))
    return 0;

  if (sha3_256 (out, 32, pubkey_hex + 1, 64))
    return 0;

  memset (out, 0, 12);

  return 1;
}

int ethjet_blake2(uint8_t *in, size_t in_size,
                  uint8_t *out, size_t out_size) {
  uint32_t rounds = in[0] << 24 | in[1] << 16 | in[2] << 8 | in[3];
  unsigned char f = in[212];
  uint64_t *h = (uint64_t *)&in[4];
  uint64_t *m = (uint64_t *)&in[68];
  uint64_t *t = (uint64_t *)&in[196];

  blake2b_compress(h, m, t, f, rounds);

  memcpy(out, h, out_size);

  return 1;
}

  // for loading an element of F_q (a coordinate of G_1)
  // consumes 32 bytes
  int read_Fq_element_v2 (uint8_t *in, mpz_t x_data) {
    mpz_init(x_data);
    mpz_import(x_data, 32, 1, sizeof(in[0]), 1, 0, in);

    // mpz_t q;
    // mpz_init(q);
    // alt_bn128_modulus_q.to_mpz(q);
    // const mp_size_t limbs = alt_bn128_q_limbs;

    // if (mpz_cmp(x_data, q) >= 0)
    //   throw 0;
    return 0;
  }

  // for loading an element of F_{q^2} (a coordinate of G_2)
  // consumes 64 bytes
  int read_Fq2_element_v2 (uint8_t *in, mpz_t x0_data, mpz_t x1_data) {
    // suprising "big-endian" encoding
    read_Fq_element_v2(in+32, x0_data);
    read_Fq_element_v2(in, x1_data);

    return 0;
  }

  // for loading an element of F_r (a scalar for G_1)
  // consumes 32 bytes
  int read_Fr_element (uint8_t *in, mpz_t x_data) {
    mpz_init(x_data);
    mpz_import(x_data, 32, 1, sizeof(in[0]), 1, 0, in);
    return 0;
  }

  // for loading an element of F_q (a coordinate of G_1)
  // consumes 32 bytes
  int write_Fq_element_v2 (uint8_t *out, mpz_t x_data) {
    size_t x_size;
    uint8_t *x_arr = (uint8_t *)mpz_export(NULL, &x_size, 1, 1, 1, 0, x_data);
    if (x_size > 32) {
      puts("ERR!!!");
      return 0;
    }
    // copy the result to the output buffer
    // with padding
    for (size_t i = 1; i <= 32; i++) {
      if (i <= x_size)
        out[32-i] = x_arr[x_size-i];
      else
        out[32-i] = 0;
    }
    return 0; 
  }

  // for loading a point in G_1
  // consumes 64 bytes
  void read_G1_point_v2 (uint8_t *in, curvepoint_fp_t p) {
    mpz_t ax, ay, az;
    read_Fq_element_v2(in, ax);
    read_Fq_element_v2(in+32, ay);

    // create curve point from affine coordinates
    // the point at infinity (0,0) is a special case
    if (mpz_sgn(ax) == 0 && mpz_sgn(ay) == 0) {
        //curvepoint_fp_init_set_str(p, "0", "1", "0");
        curvepoint_fp_init(p); // initialize to infinity
        puts("INF POINT:");
        curvepoint_fp_print(stdout,p);
    }
    else {
      mpz_init_set_ui(az, 1);
      gmp_printf("READ X: %Zd \n", ax);
      gmp_printf("READ Y: %Zd \n", ay);
      gmp_printf("READ Z: %Zd \n", az);

      curvepoint_fp_init_set_mpz(p, ax, ay, az);
      puts("POINT:");
      curvepoint_fp_print(stdout,p);
    }
  }

  // for loading a point in G_2
  // consumes 128 bytes
  void read_G2_point (uint8_t *in, twistpoint_fp2_t a) {
    mpz_t ax0, ax1, ay0, ay1;
    read_Fq2_element_v2(in, ax0, ax1);
    read_Fq2_element_v2(in+64, ay0, ay1);

    // create curve point from affine coordinates
    // the point at infinity (0,0) is a special case
    if (mpz_sgn(ax0) == 0 && mpz_sgn(ax1) == 0 && mpz_sgn(ay0) == 0 && mpz_sgn(ay1) == 0) {
      twistpoint_fp2_init(a);
      puts("INF TWIST");
    } else {
      twistpoint_fp2_affineset_mpz(a, ax0, ax1, ay0, ay1);
    }

    // a = alt_bn128_G2(ax, ay, alt_bn128_Fq2::one());
    // if (! a.is_well_formed()) {
    //   throw 0;
    // }
    // // additionally check that the element has the right order
    // if (-alt_bn128_Fr::one() * a + a != alt_bn128_G2::G2_zero) {
    //   throw 0;
    // }
  }

  // writes a point of G1
  // produces 64 bytes
  void write_G1_point_v2(uint8_t *out, curvepoint_fp_t p) {
    mpz_t x, y, z;
    mpz_init(x);
    mpz_init(y);
    mpz_init(z);
    
    //curvepoint_fp_makeaffine(p);
    puts("WRITING");
    curvepoint_fp_print(stdout,p);
    curvepoint_fp_get_mpz(p, x, y, z);

    gmp_printf("WRITE X: %Zd \n", x);
    gmp_printf("WRITE Y: %Zd \n", y);
    gmp_printf("WRITE Z: %Zd \n", z);

    if (!mpz_cmp_ui(x, 0) && !mpz_cmp_ui(y, 1) && !mpz_cmp_ui(z, 0)) {
      puts("INF!");
      mpz_set_ui(y, 0);
    }

    write_Fq_element_v2(out,    x);
    write_Fq_element_v2(out+32, y);
    return;
  }

int ethjet_ecadd_v2 (uint8_t *in, size_t in_size,
                     uint8_t *out, size_t out_size) {
  if (in_size != 128) {
    return 0;
  }
  if (out_size != 64) {
    return 0;
  }

  curvepoint_fp_t a;
  curvepoint_fp_t b, c;
  read_G1_point_v2(in, a);
  read_G1_point_v2(in+64, b);

  //curvepoint_fp_makeaffine(a); 
  //curvepoint_fp_makeaffine(b);

  if (!curvepoint_fp_well_formed(a) || !curvepoint_fp_well_formed(b))
    return 0;

  if (fpe_iszero(a->m_x) && fpe_isone(a->m_y) && fpe_iszero(a->m_z)) {
    write_G1_point_v2(out, b);
    return 1;
  }

  if (fpe_iszero(b->m_x) && fpe_isone(b->m_y) && fpe_iszero(b->m_z)) {
    write_G1_point_v2(out, a);
    return 1;
  }

  //puts("BEFORE");
  //curvepoint_fp_print(stdout,c);

// curvepoint_fp_makeaffine(a); 
//   curvepoint_fp_makeaffine(b); 
  curvepoint_fp_mixadd(b, a, b);
  curvepoint_fp_makeaffine(b); 

  write_G1_point_v2(out, b);

  return 1;
}

    int
  ethjet_ecmul_v2 (uint8_t *in, size_t in_size,
                uint8_t *out, size_t out_size) {

    if (in_size != 96) {
      return 0;
    }
    if (out_size != 64) {
      return 0;
    }

    curvepoint_fp_t a;
    mpz_t n;
   // curvepoint_fp_t ;
    read_G1_point_v2(in, a);
    read_Fr_element(in+64, n);
    //curvepoint_fp_makeaffine(a);
    //curvepoint_fp_makeaffine(b);


  if (!curvepoint_fp_well_formed(a))
    return 0;
  // if (fpe_iszero(a->m_x) && fpe_isone(a->m_y) && fpe_iszero(a->m_z)) {
  //   write_G1_point_v2(out, a);
  //   return 1;
  // }

    if (mpz_sgn(n) == 0) {
      curvepoint_fp_init(a);
    } else {
      curvepoint_fp_mul(a, a, n);
      curvepoint_fp_makeaffine(a);
      }
     write_G1_point_v2(out, a);

      // alt_bn128_G1 a = read_G1_point(in);
      // alt_bn128_Fr n = read_Fr_element(in+64);
      // alt_bn128_G1 na = n * a;

      // write_G1_point(out, na);

    return 1;
  }

int
ethjet (struct ethjet_context *ctx,
        enum ethjet_operation op,
        uint8_t *in, size_t in_size,
        uint8_t *out, size_t out_size)
{
  switch (op) {
  case ETHJET_ECRECOVER:
    return ethjet_ecrecover (ctx->ec, in, in_size, out, out_size);
    break;

  case ETHJET_EXAMPLE:
    return ethjet_example (ctx, in, in_size, out, out_size);

  case ETHJET_ECADD:
    return ethjet_ecadd_v2 (in, in_size, out, out_size);

  case ETHJET_ECMUL:
    return ethjet_ecmul_v2 (in, in_size, out, out_size);

  case ETHJET_ECPAIRING:
    return ethjet_ecpairing (in, in_size, out, out_size);

  case ETHJET_BLAKE2:
    return ethjet_blake2 (in, in_size, out, out_size);

  default:
    return 0;
  }
}
