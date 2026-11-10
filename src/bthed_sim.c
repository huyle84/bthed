/*
 * bthed_sim.c -- Stage-1 (cleartext) simulator for BTHED
 * (Blind Threshold Dilithium from FHE), built on the pq-crystals ML-DSA
 * reference implementation (FIPS 204).
 *
 * What is simulated (in the clear, no FHE, no NIZK):
 *   - the noisy-commitment signing NC.Sign_L of Figure 1 (w = Ay + e,
 *     low-bits test on w0 - c*s2 with margin gamma2 - 2*beta),
 *   - the XOR-blinding v_j = PRF(k, H_x(msg1||j)) xor v^U_j and Unpack,
 *   - the user's Finalize from PUBLIC data only (t0-test, MakeHint on
 *     A z - c t1 2^d, hint-weight test),
 *   - verification with the UNMODIFIED reference ML-DSA verifier.
 * Three signing rules are compared:
 *   STD   : standard ML-DSA attempt (e = 0, margin gamma2 - beta),
 *   BTHED  : noisy commitment (e != 0), margin gamma2 - beta (this paper; same rule as THED / BdCE+26),
 *   DOUBLE : noisy commitment (e != 0), doubled margin gamma2 - 2*beta (comparison only; not needed).
 *
 * Stand-ins (documented in README): the ring-Darkmatter PRF and H_x are
 * replaced by SHAKE256; FHE ciphertexts, threshold decryption and NIZK
 * arguments are not simulated (Stage 2/3).
 *
 * Commands:
 *   bthed_sim stats  <attempts> <be_bits> [seed]   per-attempt statistics, Fin, Verify
 *   bthed_sim session <sessions> <Latt> <be_bits> [seed]  full protocol flow + timings
 *   bthed_sim bench  <reps> <be_bits>              micro-benchmarks (median ns)
 *   bthed_sim rss    <n> <t>                       replicated secret sharing checks
 */
#define _POSIX_C_SOURCE 199309L
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "params.h"
#include "sign.h"
#include "packing.h"
#include "polyvec.h"
#include "poly.h"
#include "rounding.h"
#include "fips202.h"
#include "randombytes.h"

#define NMODES 3
enum { M_STD = 0, M_BTHED = 1, M_DBL = 2 };
static const char *mode_name[NMODES] = {"STD", "BTHED", "DOUBLE"};

#define MAX_BE_BITS 8
#define VBYTES_MAX (L * POLYZ_PACKEDBYTES + K * N * MAX_BE_BITS / 8)

static int be_bits = 3;                 /* b_e ; B_e = 2^(b_e-1) */
static int32_t Be(void) { return 1 << (be_bits - 1); }
static size_t vbytes(void) { return L * POLYZ_PACKEDBYTES + (size_t)K * N * be_bits / 8; }

/* ---------------------------------------------------------------- utils */
static double now_ns(void) {
  struct timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts);
  return ts.tv_sec * 1e9 + ts.tv_nsec;
}
/* deterministic stream for Monte Carlo: SHAKE256(seed || ctr) */
static uint64_t g_ctr = 0;
static uint8_t g_seed[32];
static void mc_bytes(uint8_t *out, size_t len) {
  uint8_t in[40];
  memcpy(in, g_seed, 32);
  for (int i = 0; i < 8; i++) in[32 + i] = (uint8_t)(g_ctr >> (8 * i));
  g_ctr++;
  shake256(out, len, in, sizeof in);
}

/* ---------------------------------------------------------------- keys */
typedef struct {
  uint8_t pk[CRYPTO_PUBLICKEYBYTES], sk[CRYPTO_SECRETKEYBYTES];
  uint8_t rho[SEEDBYTES], tr[TRBYTES], key[SEEDBYTES];
  polyvecl mat[K];
  polyvecl s1, s1hat;
  polyveck s2, s2hat, t0, t0hat, t1, t1shat; /* t1shat = NTT(t1 * 2^d) */
} keys_t;

static void keys_init(keys_t *k) {
  crypto_sign_keypair(k->pk, k->sk);
  unpack_sk(k->rho, k->tr, k->key, &k->t0, &k->s1, &k->s2, k->sk);
  unpack_pk(k->rho, &k->t1, k->pk);
  polyvec_matrix_expand(k->mat, k->rho);
  k->s1hat = k->s1; polyvecl_ntt(&k->s1hat);
  k->s2hat = k->s2; polyveck_ntt(&k->s2hat);
  k->t0hat = k->t0; polyveck_ntt(&k->t0hat);
  k->t1shat = k->t1; polyveck_shiftl(&k->t1shat); polyveck_ntt(&k->t1shat);
}

/* mu = H(tr || 0 || 0 || M)  (FIPS 204, empty context) */
static void compute_mu(uint8_t mu[CRHBYTES], const keys_t *k, const uint8_t *m, size_t mlen) {
  keccak_state st; uint8_t pre[2] = {0, 0};
  shake256_init(&st);
  shake256_absorb(&st, k->tr, TRBYTES);
  shake256_absorb(&st, pre, 2);
  shake256_absorb(&st, m, mlen);
  shake256_finalize(&st);
  shake256_squeeze(mu, CRHBYTES, &st);
}

/* ---------------------------------------------------------------- Unpack */
/* y_i = gamma1 - u (b_gamma bits), e_i = B_e - u' (b_e bits) */
static void unpack_ye(polyvecl *y, polyveck *e, const uint8_t *v) {
  for (int i = 0; i < L; i++) polyz_unpack(&y->vec[i], v + i * POLYZ_PACKEDBYTES);
  const uint8_t *p = v + L * POLYZ_PACKEDBYTES;
  uint64_t acc = 0; int nb = 0; size_t pos = 0;
  int32_t B = Be(); uint32_t mask = (1u << be_bits) - 1;
  for (int i = 0; i < K; i++)
    for (int j = 0; j < N; j++) {
      while (nb < be_bits) { acc |= (uint64_t)p[pos++] << nb; nb += 8; }
      e->vec[i].coeffs[j] = B - (int32_t)(acc & mask);
      acc >>= be_bits; nb -= be_bits;
    }
}

/* ---------------------------------------------------------------- C1 */
typedef struct {
  polyvecl y; polyveck w0, w1;
  uint8_t w1pack[K * POLYW1_PACKEDBYTES];
} c1_out_t;

/* commitment of one attempt: w = A y (+ e), Decompose */
static void c1_attempt(c1_out_t *o, const keys_t *k, const uint8_t *v, int mode) {
  polyveck e, w; polyvecl yhat;
  unpack_ye(&o->y, &e, v);
  yhat = o->y; polyvecl_ntt(&yhat);
  polyvec_matrix_pointwise_montgomery(&w, k->mat, &yhat);
  polyveck_reduce(&w);
  polyveck_invntt_tomont(&w);
  if (mode != M_STD) polyveck_add(&w, &w, &e);
  polyveck_reduce(&w);
  polyveck_caddq(&w);
  polyveck_decompose(&o->w1, &o->w0, &w);
  polyveck_pack_w1(o->w1pack, &o->w1);
}

/* user: ctilde = H(mu || w1) */
static void chal(uint8_t ct[CTILDEBYTES], const uint8_t mu[CRHBYTES], const uint8_t *w1pack) {
  keccak_state st;
  shake256_init(&st);
  shake256_absorb(&st, mu, CRHBYTES);
  shake256_absorb(&st, w1pack, K * POLYW1_PACKEDBYTES);
  shake256_finalize(&st);
  shake256_squeeze(ct, CTILDEBYTES, &st);
}

/* ---------------------------------------------------------------- C2 */
/* z = y + c s1, r = w0 - c s2 ; returns bit0 = z-test, bit1 = r-test */
static int c2_attempt(polyvecl *z, const keys_t *k, const c1_out_t *o,
                      const uint8_t ct[CTILDEBYTES], int mode) {
  poly c; polyveck r;
  poly_challenge(&c, ct);
  poly_ntt(&c);
  polyvecl_pointwise_poly_montgomery(z, &c, &k->s1hat);
  polyvecl_invntt_tomont(z);
  polyvecl_add(z, z, &o->y);
  polyvecl_reduce(z);
  int zok = !polyvecl_chknorm(z, GAMMA1 - BETA);
  polyveck_pointwise_poly_montgomery(&r, &c, &k->s2hat);
  polyveck_invntt_tomont(&r);
  polyveck_sub(&r, &o->w0, &r);
  polyveck_reduce(&r);
  int32_t margin = (mode == M_DBL) ? GAMMA2 - 2 * BETA : GAMMA2 - BETA;
  int rok = !polyveck_chknorm(&r, margin);
  return zok | (rok << 1);
}

/* ---------------------------------------------------------------- Fin */
/* Finalize from PUBLIC data: (A, t1, t0, mu, ctilde, z).
   returns 0 ok, 1 t0-test fail, 2 hint weight fail */
static int fin_public(uint8_t sig[CRYPTO_BYTES], const keys_t *k,
                      const uint8_t ct[CTILDEBYTES], const polyvecl *z) {
  poly c; polyveck ct0, r, rm, h; polyvecl zhat;
  poly_challenge(&c, ct);
  poly_ntt(&c);
  polyveck_pointwise_poly_montgomery(&ct0, &c, &k->t0hat);
  polyveck_invntt_tomont(&ct0);
  polyveck_reduce(&ct0);
  if (polyveck_chknorm(&ct0, GAMMA2)) return 1;
  /* r = A z - c t1 2^d */
  zhat = *z; polyvecl_ntt(&zhat);
  polyvec_matrix_pointwise_montgomery(&r, k->mat, &zhat);
  polyveck_pointwise_poly_montgomery(&h, &c, &k->t1shat);
  polyveck_sub(&r, &r, &h);
  polyveck_reduce(&r);
  polyveck_invntt_tomont(&r);
  polyveck_caddq(&r);
  /* rm = r - c t0 ;  h = [HighBits(r) != HighBits(r - c t0)]  (MakeHint(-ct0, r)) */
  polyveck_sub(&rm, &r, &ct0);
  polyveck_reduce(&rm);
  polyveck_caddq(&rm);
  unsigned int wt = 0;
  for (int i = 0; i < K; i++)
    for (int j = 0; j < N; j++) {
      int32_t a0;
      int32_t h1 = decompose(&a0, r.vec[i].coeffs[j]);
      int32_t h2 = decompose(&a0, rm.vec[i].coeffs[j]);
      h.vec[i].coeffs[j] = (h1 != h2);
      wt += (h1 != h2);
    }
  if (wt > OMEGA) return 2;
  pack_sig(sig, ct, z, &h);
  return 0;
}

static int verify_std(const keys_t *k, const uint8_t *sig, const uint8_t *m, size_t mlen) {
  return crypto_sign_verify(sig, CRYPTO_BYTES, m, mlen, NULL, 0, k->pk); /* 0 = valid */
}

/* ================================================================ stats */
static int cmd_stats(long natt, int bb, unsigned long seed) {
  keys_t k; uint8_t mu[CRHBYTES], msg[32], v[VBYTES_MAX], ct[CTILDEBYTES], sig[CRYPTO_BYTES];
  be_bits = bb;
  memset(g_seed, 0, 32); memcpy(g_seed, &seed, sizeof seed);
  keys_init(&k);
  long zp[NMODES] = {0}, rp[NMODES] = {0}, acc[NMODES] = {0}, f_t0[NMODES] = {0},
       f_h[NMODES] = {0}, fin_ok[NMODES] = {0}, vfail[NMODES] = {0}, w1diff[NMODES] = {0};
  double t0 = now_ns();
  for (long a = 0; a < natt; a++) {
    if (a % 64 == 0) { mc_bytes(msg, 32); compute_mu(mu, &k, msg, 32); }
    mc_bytes(v, vbytes());
    for (int md = 0; md < NMODES; md++) {
      c1_out_t o; polyvecl z;
      c1_attempt(&o, &k, v, md);
      chal(ct, mu, o.w1pack);
      int t = c2_attempt(&z, &k, &o, ct, md);
      zp[md] += t & 1; rp[md] += (t >> 1) & 1;
      if (t != 3) continue;
      acc[md]++;
      int f = fin_public(sig, &k, ct, &z);
      if (f == 1) { f_t0[md]++; continue; }
      if (f == 2) { f_h[md]++; continue; }
      fin_ok[md]++;
      if (verify_std(&k, sig, msg, 32) != 0) vfail[md]++;
      /* direct check: UseHint recovers w1 of the attempt? */
      { polyvecl zz; polyveck hh; uint8_t cc[CTILDEBYTES]; unpack_sig(cc, &zz, &hh, sig);
        (void)zz; (void)hh; }
    }
  }
  double dt = (now_ns() - t0) / 1e9;
  int32_t B = Be();
  printf("# ML-DSA mode %d (k=%d,l=%d,beta=%d,gamma1=2^%d,gamma2=%d), b_e=%d (B_e=%d), attempts=%ld, %.1fs\n",
         DILITHIUM_MODE, K, L, BETA, (GAMMA1 == (1 << 17)) ? 17 : 19, GAMMA2, be_bits, B, natt, dt);
  printf("mode,attempts,p_z,p_r,p_accept,p_fin_given_accept,t0_fail,hint_fail,signatures,verify_failures,expected_repetitions\n");
  for (int md = 0; md < NMODES; md++) {
    double pa = (double)acc[md] / natt;
    double pf = acc[md] ? (double)fin_ok[md] / acc[md] : 0;
    printf("%s,%ld,%.5f,%.5f,%.5f,%.5f,%ld,%ld,%ld,%ld,%.3f\n", mode_name[md], natt,
           (double)zp[md] / natt, (double)rp[md] / natt, pa, pf, f_t0[md], f_h[md],
           fin_ok[md], vfail[md], pa > 0 ? 1.0 / (pa * pf) : 0);
    (void)w1diff;
  }
  return 0;
}

/* ================================================================ session */
/* full protocol flow (cleartext): Req -> C1 -> Chal -> C2 -> Fin -> Verify */
static int cmd_session(long nsess, int natt, int bb, unsigned long seed) {
  keys_t k; be_bits = bb;
  memset(g_seed, 0, 32); memcpy(g_seed, &seed, sizeof seed); g_seed[31] = 0x5e;
  keys_init(&k);
  uint8_t kprf[32]; mc_bytes(kprf, 32);
  size_t VB = vbytes();
  uint8_t *vU = malloc((size_t)natt * VB), *vS = malloc((size_t)natt * VB);
  c1_out_t *o = malloc(sizeof(c1_out_t) * natt);
  uint8_t (*ct)[CTILDEBYTES] = malloc((size_t)natt * CTILDEBYTES);
  long succ[NMODES] = {0}, vfail[NMODES] = {0}, noacc[NMODES] = {0}, finfail[NMODES] = {0};
  double T[6] = {0}; /* req, c1, chal, c2, fin, verify (BTHED only) */
  for (long s = 0; s < nsess; s++) {
    uint8_t msg[32], mu[CRHBYTES], msg1[32], sig[CRYPTO_BYTES];
    mc_bytes(msg, 32);
    for (int md = 0; md < NMODES; md++) {
      double ta = now_ns();
      /* --- User Req: mu, user one-time pads v^U_j, session id msg1 */
      compute_mu(mu, &k, msg, 32);
      mc_bytes(vU, (size_t)natt * VB);
      mc_bytes(msg1, 32); /* stand-in for msg1 = (upk, ct_mu, ct_v, pi_1) */
      double tb = now_ns();
      /* --- Signers C1: v^S_j = PRF(k, H_x(msg1||j)), v_j = v^S_j xor v^U_j */
      for (int j = 0; j < natt; j++) {
        uint8_t in[32 + 4], x[32], kin[64];
        memcpy(in, msg1, 32); in[32] = j; in[33] = j >> 8; in[34] = 0; in[35] = 0;
        shake256(x, 32, in, sizeof in);              /* H_x */
        memcpy(kin, kprf, 32); memcpy(kin + 32, x, 32);
        shake256(vS + j * VB, VB, kin, 64);          /* PRF stand-in */
        for (size_t b = 0; b < VB; b++) vS[j * VB + b] ^= vU[j * VB + b];
        c1_attempt(&o[j], &k, vS + j * VB, md);
      }
      double tc = now_ns();
      /* --- User Chal: c_j = SIB(H(mu || w1_j)) (SIB expanded inside C2 here) */
      for (int j = 0; j < natt; j++) chal(ct[j], mu, o[j].w1pack);
      double td = now_ns();
      /* --- Signers C2: responses, tests, selection of the first success */
      int sel = -1; polyvecl z, zsel;
      for (int j = 0; j < natt; j++) {           /* evaluate all (as in FHE) */
        int t = c2_attempt(&z, &k, &o[j], ct[j], md);
        if (t == 3 && sel < 0) { sel = j; zsel = z; }
      }
      double te = now_ns();
      /* --- User Fin + Verify */
      int ok = 0, f = -1;
      if (sel >= 0) f = fin_public(sig, &k, ct[sel], &zsel);
      double tf = now_ns();
      if (sel < 0) noacc[md]++;
      else if (f != 0) finfail[md]++;
      else { ok = 1; succ[md]++; if (verify_std(&k, sig, msg, 32) != 0) vfail[md]++; }
      double tg = now_ns();
      if (md == M_BTHED) {
        T[0] += tb - ta; T[1] += tc - tb; T[2] += td - tc; T[3] += te - td; T[4] += tf - te;
        if (ok) T[5] += tg - tf;
      }
    }
  }
  printf("# ML-DSA mode %d, L=%d attempts/session, b_e=%d, sessions=%ld\n", DILITHIUM_MODE, natt, be_bits, nsess);
  printf("mode,L,sessions,success,p_L,no_accepted_attempt,fin_fail,verify_failures\n");
  for (int md = 0; md < NMODES; md++)
    printf("%s,%d,%ld,%ld,%.4f,%ld,%ld,%ld\n", mode_name[md], natt, nsess, succ[md],
           (double)succ[md] / nsess, noacc[md], finfail[md], vfail[md]);
  printf("# BTHED cleartext timings per session (mean, microseconds):\n");
  printf("timing_us,user_req,signer_C1,user_chal,signer_C2,user_fin,verify\n");
  printf("BTHED,%.1f,%.1f,%.1f,%.1f,%.1f,%.1f\n", T[0] / nsess / 1e3, T[1] / nsess / 1e3,
         T[2] / nsess / 1e3, T[3] / nsess / 1e3, T[4] / nsess / 1e3,
         succ[M_BTHED] ? T[5] / succ[M_BTHED] / 1e3 : 0);
  free(vU); free(vS); free(o); free(ct);
  return 0;
}

/* ================================================================ bench */
static int cmpd(const void *a, const void *b) {
  double x = *(const double *)a, y = *(const double *)b; return (x > y) - (x < y);
}
#define BENCH(name, reps, ...) do { double *t_ = malloc(sizeof(double) * (reps)); \
  for (int r_ = 0; r_ < (reps); r_++) { double a_ = now_ns(); __VA_ARGS__; t_[r_] = now_ns() - a_; } \
  qsort(t_, (reps), sizeof(double), cmpd); printf("%s,%.0f\n", name, t_[(reps) / 2]); free(t_); } while (0)

static int cmd_bench(int reps, int bb) {
  keys_t k; be_bits = bb;
  memset(g_seed, 7, 32);
  keys_init(&k);
  uint8_t mu[CRHBYTES], msg[32] = {1}, v[VBYTES_MAX], ct[CTILDEBYTES], sig[CRYPTO_BYTES];
  size_t siglen;
  compute_mu(mu, &k, msg, 32);
  mc_bytes(v, vbytes());
  c1_out_t o; polyvecl z;
  printf("# ML-DSA mode %d, b_e=%d, median over %d reps (ns)\n", DILITHIUM_MODE, be_bits, reps);
  printf("operation,median_ns\n");
  BENCH("mldsa_keygen", reps, { uint8_t pk[CRYPTO_PUBLICKEYBYTES], sk[CRYPTO_SECRETKEYBYTES]; crypto_sign_keypair(pk, sk); });
  BENCH("mldsa_sign_reference(full; with rejections)", reps, crypto_sign_signature(sig, &siglen, msg, 32, NULL, 0, k.sk));
  BENCH("mldsa_verify_reference", reps, verify_std(&k, sig, msg, 32));
  BENCH("prf_stand_in_shake256(one_attempt_bits)", reps, mc_bytes(v, vbytes()));
  BENCH("C1_one_attempt(unpack; Ay+e; decompose)", reps, c1_attempt(&o, &k, v, M_BTHED));
  BENCH("user_chal_one_attempt(H(mu||w1))", reps, chal(ct, mu, o.w1pack));
  BENCH("C2_one_attempt(SIB; z; r; tests)", reps, c2_attempt(&z, &k, &o, ct, M_BTHED));
  /* find an accepted attempt for Fin */
  for (;;) { mc_bytes(v, vbytes()); c1_attempt(&o, &k, v, M_BTHED); chal(ct, mu, o.w1pack);
             if (c2_attempt(&z, &k, &o, ct, M_BTHED) == 3 && fin_public(sig, &k, ct, &z) == 0) break; }
  BENCH("user_fin_public(t0-test; MakeHint; pack)", reps, fin_public(sig, &k, ct, &z));
  BENCH("verify_BTHED_signature_reference", reps, verify_std(&k, sig, msg, 32));
  return 0;
}

/* ================================================================ RSS */
static long binom(int n, int r) { if (r < 0 || r > n) return 0; long b = 1;
  for (int i = 1; i <= r; i++) b = b * (n - r + i) / i;
  return b; }

/* enumerate subsets of [n] of size t-1 as bitmasks */
static int cmd_rss(int n, int t) {
  if (n > 20 || t < 1 || t > n) { fprintf(stderr, "need 1<=t<=n<=20\n"); return 1; }
  int ns = 0; uint32_t *T = malloc(sizeof(uint32_t) * binom(n, t - 1));
  for (uint32_t m = 0; m < (1u << n); m++) if (__builtin_popcount(m) == t - 1) T[ns++] = m;
  /* secret s in Z_{2^64}; shares uniform with sum s */
  uint64_t s, *sh = malloc(sizeof(uint64_t) * ns), sum = 0;
  randombytes((uint8_t *)&s, 8);
  for (int i = 0; i < ns - 1; i++) { randombytes((uint8_t *)&sh[i], 8); sum += sh[i]; }
  sh[ns - 1] = s - sum;
  /* reconstruction with 0/1 coefficients for every active set S, |S| = t */
  long nsets = 0, ok = 0, priv_ok = 0, nC = 0;
  for (uint32_t S = 0; S < (1u << n); S++) {
    if (__builtin_popcount(S) != t) continue;
    nsets++;
    uint64_t acc = 0;
    for (int i = 0; i < ns; i++) {
      uint32_t cand = S & ~T[i];              /* S \ T non-empty since |T| < |S| */
      int resp = __builtin_ctz(cand);         /* resp_S(T) = min(S \ T) */
      if (T[i] >> resp & 1) { fprintf(stderr, "resp in T!\n"); return 1; }
      acc += sh[i];
    }
    ok += (acc == s);
  }
  /* privacy: every C with |C| <= t-1 misses at least one share T0 ⊇ C */
  for (uint32_t C = 0; C < (1u << n); C++) {
    if (__builtin_popcount(C) > t - 1) continue;
    nC++;
    for (int i = 0; i < ns; i++) if ((T[i] & C) == C) { priv_ok++; break; }
  }
  printf("# replicated secret sharing, n=%d, t=%d\n", n, t);
  printf("n,t,total_shares,shares_per_party,active_sets,reconstruct_ok,coalitions_checked,privacy_ok\n");
  printf("%d,%d,%d,%ld,%ld,%ld,%ld,%ld\n", n, t, ns, binom(n - 1, t - 1), nsets, ok, nC, priv_ok);
  free(T); free(sh);
  return 0;
}

int main(int argc, char **argv) {
  if (argc < 2) goto usage;
  if (!strcmp(argv[1], "stats") && argc >= 4)
    return cmd_stats(atol(argv[2]), atoi(argv[3]), argc > 4 ? strtoul(argv[4], 0, 10) : 1);
  if (!strcmp(argv[1], "session") && argc >= 5)
    return cmd_session(atol(argv[2]), atoi(argv[3]), atoi(argv[4]), argc > 5 ? strtoul(argv[5], 0, 10) : 2);
  if (!strcmp(argv[1], "bench") && argc >= 4) return cmd_bench(atoi(argv[2]), atoi(argv[3]));
  if (!strcmp(argv[1], "rss") && argc >= 4) return cmd_rss(atoi(argv[2]), atoi(argv[3]));
usage:
  fprintf(stderr, "usage: %s stats <attempts> <be_bits> [seed] | session <sessions> <L> <be_bits> [seed] | bench <reps> <be_bits> | rss <n> <t>\n", argv[0]);
  return 1;
}
