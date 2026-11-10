/*
 * BTHED - Stage 2a: two-key threshold RLWE layer for the arithmetic core.
 *
 * What is real here:
 *   - a two-key (user U + committee S) BGV-type scheme over R_Q = Z_Q[X]/(X^NN+1), RNS + NTT,
 *     plaintext modulus t = q_Dilithium = 8380417, Dilithium ring embedded by Y = X^(NN/256);
 *   - CDKS19-style multi-key relinearization from individually generated keys (common random a);
 *   - replicated secret sharing of the committee key, partial decryptions with Gaussian flooding,
 *     decryption towards the user only (the user adds its own key component);
 *   - homomorphic evaluation of the arithmetic part of C1/C2 on real ML-DSA data:
 *       H1: w = A y + e      H2: z = y + c s1      H3: u = w - c s2      H4: f * z
 *   - the user's Fin on public data and the unmodified reference ML-DSA verifier.
 * What is emulated in the clear (marked ORACLE): the non-arithmetic steps of C1/C2
 *   (PRF, XOR, Unpack, Decompose/HighBits, the two norm tests, selection of the first accepted
 *   attempt). They need bootstrapped comparisons (THED's CKKS-BFV machinery) and are Stage 2b/3.
 * Randomness of the FHE layer comes from xoshiro256** (NOT cryptographic): prototype only.
 */
#define _POSIX_C_SOURCE 199309L
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "params.h"
#include "sign.h"
#include "packing.h"
#include "polyvec.h"
#include "poly.h"
#include "rounding.h"
#include "reduce.h"
#include "fips202.h"
#include "randombytes.h"
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define NN 8192
#define LOGNN 13
#define NP 6
#define DEMB (NN / N) /* 32 Dilithium polynomials fit in one plaintext */
#define TT ((uint32_t)Q)
#define QDEC_LOG2 64.0 /* decryption-query budget used for the flooding parameter */
#define BE_BITS 3      /* commitment noise: B_e = 4, as THED */

typedef struct { uint32_t c[NP][NN]; } rp;  /* ring element; NTT form unless stated */
typedef struct { rp c[3]; } ct_t;           /* c[0] + c[1]*s_S + c[2]*s_U = m + t*e   */
enum { PS = 1, PU = 2 };
typedef struct { rp s, *b, *d0, *d1, *d2; } party_t;

static uint32_t P[NP], tmod[NP], ninv[NP], ninvs[NP], ginv[NP][NP], ppref_t[NP], Qmodt;
static uint32_t *psi[NP], *psis[NP], *ipsi[NP], *ipsis[NP];
static double lgpref[NP + 1];
static rp *CRS, *DG, *DG2, *MONO;
static party_t party[3];

static double now_ms(void) { struct timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts); return ts.tv_sec * 1e3 + ts.tv_nsec * 1e-6; }
static void *xalloc(size_t n) { void *p = calloc(1, n); if (!p) { fprintf(stderr, "out of memory\n"); exit(1); } return p; }

/* ------------------------------------------------------------ PRNG (prototype only) */
static uint64_t xs[4];
static inline uint64_t rotl(uint64_t x, int k) { return (x << k) | (x >> (64 - k)); }
static uint64_t rng(void) { uint64_t r = rotl(xs[1] * 5, 7) * 9, t = xs[1] << 17; xs[2] ^= xs[0]; xs[3] ^= xs[1]; xs[1] ^= xs[2]; xs[0] ^= xs[3]; xs[2] ^= t; xs[3] = rotl(xs[3], 45); return r; }
static void rng_seed(uint64_t s) { for (int i = 0; i < 4; i++) { s += 0x9e3779b97f4a7c15ULL; uint64_t z = s; z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL; z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL; xs[i] = z ^ (z >> 31); } }
static double rng_u(void) { return ((rng() >> 11) + 0.5) * (1.0 / 9007199254740992.0); }
static double gauss(void) { static int have = 0; static double spare; if (have) { have = 0; return spare; }
  double u = rng_u(), v = rng_u(), r = sqrt(-2 * log(u)); spare = r * sin(2 * M_PI * v); have = 1; return r * cos(2 * M_PI * v); }

/* ------------------------------------------------------------ modular arithmetic, NTT */
static uint32_t powmod(uint32_t a, uint64_t e, uint32_t p) { uint64_t r = 1, b = a % p; while (e) { if (e & 1) r = r * b % p; b = b * b % p; e >>= 1; } return (uint32_t)r; }
static int is_prime(uint32_t n) { if (n < 4) return n > 1; if (!(n & 1)) return 0; uint32_t d = n - 1; int s = 0; while (!(d & 1)) { d >>= 1; s++; }
  const uint32_t bases[3] = {2, 7, 61};
  for (int i = 0; i < 3; i++) { uint32_t a = bases[i] % n; if (!a) continue; uint64_t x = powmod(a, d, n); if (x == 1 || x == n - 1) continue; int ok = 0;
    for (int r = 1; r < s; r++) { x = x * x % n; if (x == n - 1) { ok = 1; break; } } if (!ok) return 0; } return 1; }
static inline uint32_t shoup(uint32_t w, uint32_t p) { return (uint32_t)(((uint64_t)w << 32) / p); }
static inline uint32_t mulsh(uint32_t x, uint32_t w, uint32_t ws, uint32_t p) { uint32_t q = (uint32_t)(((uint64_t)x * ws) >> 32); uint32_t r = x * w - q * p; return r >= p ? r - p : r; }
static uint32_t brv(uint32_t x) { uint32_t r = 0; for (int i = 0; i < LOGNN; i++) { r = (r << 1) | (x & 1); x >>= 1; } return r; }

static void fntt(uint32_t *a, int i) { uint32_t p = P[i]; int t = NN;
  for (int m = 1; m < NN; m <<= 1) { t >>= 1;
    for (int k = 0; k < m; k++) { int j1 = 2 * k * t, j2 = j1 + t; uint32_t w = psi[i][m + k], ws = psis[i][m + k];
      for (int j = j1; j < j2; j++) { uint32_t u = a[j], v = mulsh(a[j + t], w, ws, p), s = u + v, d = u + p - v; a[j] = s >= p ? s - p : s; a[j + t] = d >= p ? d - p : d; } } } }
static void fintt(uint32_t *a, int i) { uint32_t p = P[i]; int t = 1;
  for (int m = NN; m > 1; m >>= 1) { int j1 = 0, h = m >> 1;
    for (int k = 0; k < h; k++) { int j2 = j1 + t; uint32_t w = ipsi[i][h + k], ws = ipsis[i][h + k];
      for (int j = j1; j < j2; j++) { uint32_t u = a[j], v = a[j + t], s = u + v, d = u + p - v; a[j] = s >= p ? s - p : s; a[j + t] = mulsh(d >= p ? d - p : d, w, ws, p); }
      j1 += 2 * t; } t <<= 1; }
  for (int j = 0; j < NN; j++) a[j] = mulsh(a[j], ninv[i], ninvs[i], p); }

static void rp_ntt(rp *a) { for (int i = 0; i < NP; i++) fntt(a->c[i], i); }
static void rp_intt(rp *a) { for (int i = 0; i < NP; i++) fintt(a->c[i], i); }
static void rp_add(rp *r, const rp *a, const rp *b) { for (int i = 0; i < NP; i++) { uint32_t p = P[i]; for (int j = 0; j < NN; j++) { uint32_t s = a->c[i][j] + b->c[i][j]; r->c[i][j] = s >= p ? s - p : s; } } }
static void rp_sub(rp *r, const rp *a, const rp *b) { for (int i = 0; i < NP; i++) { uint32_t p = P[i]; for (int j = 0; j < NN; j++) { uint32_t s = a->c[i][j] + p - b->c[i][j]; r->c[i][j] = s >= p ? s - p : s; } } }
static void rp_mul(rp *r, const rp *a, const rp *b) { for (int i = 0; i < NP; i++) { uint32_t p = P[i]; for (int j = 0; j < NN; j++) r->c[i][j] = (uint32_t)((uint64_t)a->c[i][j] * b->c[i][j] % p); } }
static void rp_muladd(rp *r, const rp *a, const rp *b) { for (int i = 0; i < NP; i++) { uint32_t p = P[i]; for (int j = 0; j < NN; j++) r->c[i][j] = (uint32_t)(((uint64_t)a->c[i][j] * b->c[i][j] + r->c[i][j]) % p); } }
static void rp_mult(rp *r) { for (int i = 0; i < NP; i++) { uint32_t p = P[i]; for (int j = 0; j < NN; j++) r->c[i][j] = (uint32_t)((uint64_t)r->c[i][j] * tmod[i] % p); } }
static void rp_from_i32(rp *r, const int32_t *v) { for (int i = 0; i < NP; i++) for (int j = 0; j < NN; j++) r->c[i][j] = v[j] >= 0 ? (uint32_t)v[j] : P[i] - (uint32_t)(-v[j]); rp_ntt(r); }
static void rp_uniform(rp *r) { for (int i = 0; i < NP; i++) for (int j = 0; j < NN; j++) r->c[i][j] = (uint32_t)(rng() % P[i]); }
static void rp_small(rp *r, int ternary) { static int32_t v[NN];
  for (int j = 0; j < NN; j++) { uint64_t x = rng();
    if (ternary) { while ((x & 3) == 3) x = rng(); v[j] = (int32_t)(x & 3) - 1; }
    else v[j] = __builtin_popcountll(x & 0x1FFFFF) - __builtin_popcountll((x >> 21) & 0x1FFFFF); /* sigma = 3.24 */ }
  rp_from_i32(r, v); }

static void setup(void) {
  uint32_t cand = (uint32_t)((((1ull << 31) - 1) / (2 * NN)) * (2 * NN) + 1); int n = 0;
  while (n < NP) { if (is_prime(cand) && cand != TT) P[n++] = cand; cand -= 2 * NN; }
  lgpref[0] = 0;
  for (int i = 0; i < NP; i++) { uint32_t p = P[i], g = 2, x;
    for (;; g++) { x = powmod(g, (p - 1) / (2 * NN), p); if (powmod(x, NN, p) == p - 1) break; }
    uint32_t xi = powmod(x, p - 2, p), *pw = xalloc(4 * NN), *pwi = xalloc(4 * NN);
    pw[0] = pwi[0] = 1; for (int k = 1; k < NN; k++) { pw[k] = (uint32_t)((uint64_t)pw[k - 1] * x % p); pwi[k] = (uint32_t)((uint64_t)pwi[k - 1] * xi % p); }
    psi[i] = xalloc(4 * NN); psis[i] = xalloc(4 * NN); ipsi[i] = xalloc(4 * NN); ipsis[i] = xalloc(4 * NN);
    for (uint32_t k = 0; k < NN; k++) { psi[i][k] = pw[brv(k)]; psis[i][k] = shoup(psi[i][k], p); ipsi[i][k] = pwi[brv(k)]; ipsis[i][k] = shoup(ipsi[i][k], p); }
    free(pw); free(pwi);
    ninv[i] = powmod(NN, p - 2, p); ninvs[i] = shoup(ninv[i], p); tmod[i] = TT % p;
    for (int k = 0; k < NP; k++) if (k != i) ginv[k][i] = powmod(P[k] % p, p - 2, p);
    lgpref[i + 1] = lgpref[i] + log2((double)p); }
  uint64_t pr = 1; for (int i = 0; i < NP; i++) { ppref_t[i] = (uint32_t)pr; pr = pr * (P[i] % TT) % TT; } Qmodt = (uint32_t)pr;
  /* self-test: X^(NN-1) * X = -1 */
  rp *a = xalloc(sizeof(rp)), *b = xalloc(sizeof(rp)); static int32_t v[NN]; memset(v, 0, sizeof v); v[NN - 1] = 1; rp_from_i32(a, v); v[NN - 1] = 0; v[1] = 1; rp_from_i32(b, v);
  rp_mul(a, a, b); rp_intt(a); for (int i = 0; i < NP; i++) if (a->c[i][0] != P[i] - 1 || a->c[i][1] != 0) { fprintf(stderr, "NTT self-test failed\n"); exit(1); }
  free(a); free(b);
  MONO = xalloc(DEMB * sizeof(rp)); for (int s = 0; s < DEMB; s++) { memset(v, 0, sizeof v); v[s] = 1; rp_from_i32(&MONO[s], v); }
}

/* ------------------------------------------------------------ keys, encryption */
static void keygen(party_t *p) { rp *e = xalloc(sizeof(rp)), *r = xalloc(sizeof(rp));
  p->b = xalloc(NP * sizeof(rp)); p->d0 = xalloc(NP * sizeof(rp)); p->d1 = xalloc(NP * sizeof(rp)); p->d2 = xalloc(NP * sizeof(rp));
  rp_small(&p->s, 1); rp_small(r, 1);
  for (int i = 0; i < NP; i++) {
    rp_small(e, 0); rp_mult(e); rp_mul(&p->b[i], &p->s, &CRS[i]); rp_sub(&p->b[i], e, &p->b[i]);            /* b = -s a + t e      */
    rp_uniform(&p->d1[i]); rp_small(e, 0); rp_mult(e); rp_mul(&p->d0[i], &p->s, &p->d1[i]); rp_sub(&p->d0[i], e, &p->d0[i]);
    for (int j = 0; j < NN; j++) { uint32_t s = p->d0[i].c[i][j] + r->c[i][j]; p->d0[i].c[i][j] = s >= P[i] ? s - P[i] : s; } /* + r g_i */
    rp_small(e, 0); rp_mult(e); rp_mul(&p->d2[i], r, &CRS[i]); rp_add(&p->d2[i], &p->d2[i], e);
    for (int j = 0; j < NN; j++) { uint32_t s = p->d2[i].c[i][j] + p->s.c[i][j]; p->d2[i].c[i][j] = s >= P[i] ? s - P[i] : s; } /* + s g_i */ }
  free(e); free(r); }

/* single-key encryption of m (NTT form) under party 'who' */
static void enc1(ct_t *ct, int who, const rp *m) { static rp u, e; memset(ct, 0, sizeof *ct);
  rp_small(&u, 1); rp_small(&e, 0); rp_mult(&e); rp_mul(&ct->c[0], &u, &party[who].b[0]); rp_add(&ct->c[0], &ct->c[0], &e); rp_add(&ct->c[0], &ct->c[0], m);
  rp_small(&e, 0); rp_mult(&e); rp_mul(&ct->c[who], &u, &CRS[0]); rp_add(&ct->c[who], &ct->c[who], &e); }
/* fresh two-key ciphertext of m (coefficients mod t): Enc_S(mA) + Enc_U(m - mA) */
static void enc2(ct_t *ct, const int32_t *m) { static int32_t a[NN], b[NN]; static rp ma; static ct_t tmp;
  for (int j = 0; j < NN; j++) { a[j] = (int32_t)(rng() % TT) - (int32_t)(TT / 2); int64_t d = ((int64_t)m[j] - a[j]) % (int64_t)TT; if (d > (int64_t)(TT / 2)) d -= TT; if (d < -(int64_t)(TT / 2)) d += TT; b[j] = (int32_t)d; }
  rp_from_i32(&ma, a); enc1(ct, PS, &ma); rp_from_i32(&ma, b); enc1(&tmp, PU, &ma);
  rp_add(&ct->c[0], &ct->c[0], &tmp.c[0]); memcpy(&ct->c[PU], &tmp.c[PU], sizeof(rp)); }
static void emb(int32_t *out, const int32_t *coef) { memset(out, 0, NN * sizeof(int32_t));
  for (int k = 0; k < N; k++) { int32_t c = coef[k] % Q; if (c > Q / 2) c -= Q; if (c < -(Q / 2)) c += Q; out[DEMB * k] = c; } }

/* ------------------------------------------------------------ homomorphic operations */
static void decomp(rp *dg, const rp *c) { static rp tmp; memcpy(&tmp, c, sizeof tmp); rp_intt(&tmp);
  for (int i = 0; i < NP; i++) { for (int k = 0; k < NP; k++) { uint32_t p = P[k]; for (int j = 0; j < NN; j++) { uint32_t v = tmp.c[i][j]; dg[i].c[k][j] = v >= p ? v - p : v; } } rp_ntt(&dg[i]); } }
/* replace term * s_i * s_j by a ciphertext part (CDKS19, Sect. 3) */
static void relin(ct_t *ct, const rp *term, int i, int j) { static rp cp; memset(&cp, 0, sizeof cp);
  decomp(DG, term);
  for (int d = 0; d < NP; d++) { rp_muladd(&cp, &DG[d], &party[j].b[d]); rp_muladd(&ct->c[j], &DG[d], &party[i].d2[d]); }
  decomp(DG2, &cp);
  for (int d = 0; d < NP; d++) { rp_muladd(&ct->c[0], &DG2[d], &party[i].d0[d]); rp_muladd(&ct->c[i], &DG2[d], &party[i].d1[d]); } }
/* product of a U-only and an S-only ciphertext */
static void mul_cross(ct_t *o, const ct_t *cu, const ct_t *cs) { static rp term;
  rp_mul(&o->c[0], &cu->c[0], &cs->c[0]); rp_mul(&o->c[PS], &cu->c[0], &cs->c[PS]); rp_mul(&o->c[PU], &cu->c[PU], &cs->c[0]);
  rp_mul(&term, &cs->c[PS], &cu->c[PU]); relin(o, &term, PS, PU); }
/* product of two two-key ciphertexts */
static void mul_full(ct_t *o, const ct_t *a, const ct_t *b) { static rp tss, tuu, tsu;
  rp_mul(&o->c[0], &a->c[0], &b->c[0]);
  rp_mul(&o->c[PS], &a->c[0], &b->c[PS]); rp_muladd(&o->c[PS], &a->c[PS], &b->c[0]);
  rp_mul(&o->c[PU], &a->c[0], &b->c[PU]); rp_muladd(&o->c[PU], &a->c[PU], &b->c[0]);
  rp_mul(&tss, &a->c[PS], &b->c[PS]); rp_mul(&tuu, &a->c[PU], &b->c[PU]); rp_mul(&tsu, &a->c[PS], &b->c[PU]); rp_muladd(&tsu, &a->c[PU], &b->c[PS]);
  relin(o, &tss, PS, PS); relin(o, &tuu, PU, PU); relin(o, &tsu, PS, PU); }
static void ct_mulpt_add(ct_t *o, const ct_t *a, const rp *pt) { for (int k = 0; k < 3; k++) rp_muladd(&o->c[k], &a->c[k], pt); }
static void ct_add(ct_t *o, const ct_t *a, const ct_t *b) { for (int k = 0; k < 3; k++) rp_add(&o->c[k], &a->c[k], &b->c[k]); }
static void ct_sub(ct_t *o, const ct_t *a, const ct_t *b) { for (int k = 0; k < 3; k++) rp_sub(&o->c[k], &a->c[k], &b->c[k]); }

/* ------------------------------------------------------------ decoding (Garner, no bignum) */
/* x in coefficient form; dec = centered (x mod Q) mod t; returns log2 of the largest |x| */
static double garner(int32_t *dec, const rp *x) { double mx = 0;
  for (int j = 0; j < NN; j++) { uint32_t v[NP];
    for (int i = 0; i < NP; i++) { uint64_t tmp = x->c[i][j]; uint32_t p = P[i]; for (int k = 0; k < i; k++) tmp = (tmp + p - v[k] % p) % p * ginv[k][i] % p; v[i] = (uint32_t)tmp; }
    int neg = 0; for (int i = NP - 1; i >= 0; i--) { uint32_t h = (P[i] - 1) / 2; if (v[i] > h) { neg = 1; break; } if (v[i] < h) break; }
    uint64_t r = 0; for (int i = 0; i < NP; i++) r = (r + (uint64_t)(v[i] % TT) * ppref_t[i]) % TT;
    if (neg) r = (r + TT - Qmodt) % TT;
    int32_t c = (int32_t)r; if (c > (int32_t)(TT / 2)) c -= TT; dec[j] = c;
    for (int i = NP - 1; i >= 0; i--) { uint32_t d = neg ? P[i] - 1 - v[i] : v[i]; if (d) { double b = log2((double)d + 1) + lgpref[i]; if (b > mx) mx = b; break; } } }
  return mx; }
static double debug_dec(int32_t *dec, const ct_t *ct) { static rp x, y; rp_mul(&x, &ct->c[PS], &party[PS].s); rp_mul(&y, &ct->c[PU], &party[PU].s);
  rp_add(&x, &x, &y); rp_add(&x, &x, &ct->c[0]); rp_intt(&x); return garner(dec, &x); }

/* ------------------------------------------------------------ RSS and threshold decryption */
#define MAXSH 512
static int nsh, n_sig, t_thr; static uint32_t shmask[MAXSH]; static rp *sh;
static void rss_share(int n, int t) { n_sig = n; t_thr = t; nsh = 0;
  for (uint32_t m = 0; m < (1u << n); m++) if (__builtin_popcount(m) == t - 1) { if (nsh == MAXSH) { fprintf(stderr, "too many shares\n"); exit(1); } shmask[nsh++] = m; }
  sh = xalloc((size_t)nsh * sizeof(rp)); memcpy(&sh[nsh - 1], &party[PS].s, sizeof(rp));
  for (int k = 0; k < nsh - 1; k++) { rp_uniform(&sh[k]); rp_sub(&sh[nsh - 1], &sh[nsh - 1], &sh[k]); } }
/* E ~ rounded Gaussian with parameter nu = 2^lognu (standard deviation nu / sqrt(2 pi)), 40-bit precision */
static void flood(rp *E, double lognu) { double lstd = lognu - 0.5 * log2(2 * M_PI); int s = (int)floor(lstd) - 40; if (s < 0) s = 0; if (s > 127) { fprintf(stderr, "flooding too large\n"); exit(1); }
  double stdhi = pow(2.0, lstd - s); uint32_t p2s[NP], half[NP], p64[NP];
  for (int i = 0; i < NP; i++) { p2s[i] = powmod(2, s, P[i]); half[i] = s ? powmod(2, s - 1, P[i]) : 0; p64[i] = powmod(2, 64, P[i]); }
  for (int j = 0; j < NN; j++) { int64_t hi = llround(gauss() * stdhi); uint64_t lo = 0, up = 0;
    if (s > 0) { lo = rng(); if (s < 64) lo &= (1ull << s) - 1; } if (s > 64) { up = rng(); if (s < 128) up &= (1ull << (s - 64)) - 1; }
    for (int i = 0; i < NP; i++) { uint64_t p = P[i], h = hi >= 0 ? (uint64_t)hi % p : (p - (uint64_t)(-hi) % p) % p;
      uint64_t l = ((up % p) * p64[i] + lo % p) % p; E->c[i][j] = (uint32_t)((h * p2s[i] + l + p - half[i]) % p); } }
  rp_ntt(E); }
static void partdec(rp *out, int i, uint32_t Smask, const rp *cS, double lognu) { static rp sum, E; memset(&sum, 0, sizeof sum);
  for (int k = 0; k < nsh; k++) { uint32_t rest = Smask & ~shmask[k]; if (rest && __builtin_ctz(rest) == i) rp_add(&sum, &sum, &sh[k]); }
  rp_mul(out, cS, &sum); flood(&E, lognu); rp_mult(&E); rp_add(out, out, &E); }

static double T_part, T_user, noise_fin; static long n_part;
/* threshold decryption towards the user of cnt ciphertexts, each holding one Dilithium polynomial;
   they are first packed DEMB per ciphertext (free: multiplication by X^slot) */
static void thr_decrypt(int32_t (*out)[N], ct_t **cts, int cnt, double logB, uint32_t Smask) { static ct_t pk; static rp sum, pd; static int32_t dec[NN];
  int ng = (cnt + DEMB - 1) / DEMB; double lognu = logB - log2((double)TT) + 0.5 * (log2(2 * M_PI / log(2.0)) + log2((double)NN * ng) + QDEC_LOG2);
  for (int g = 0; g < ng; g++) { memset(&pk, 0, sizeof pk);
    for (int s = 0; s < DEMB && g * DEMB + s < cnt; s++) ct_mulpt_add(&pk, cts[g * DEMB + s], &MONO[s]);
    memcpy(&sum, &pk.c[0], sizeof sum);
    for (int i = 0; i < n_sig; i++) if (Smask >> i & 1) { double t0 = now_ms(); partdec(&pd, i, Smask, &pk.c[PS], lognu); T_part += now_ms() - t0; n_part++; rp_add(&sum, &sum, &pd); }
    double t0 = now_ms(); rp_mul(&pd, &pk.c[PU], &party[PU].s); rp_add(&sum, &sum, &pd); rp_intt(&sum); double nb = garner(dec, &sum); T_user += now_ms() - t0; if (nb > noise_fin) noise_fin = nb;
    for (int s = 0; s < DEMB && g * DEMB + s < cnt; s++) for (int k = 0; k < N; k++) out[g * DEMB + s][k] = dec[DEMB * k + s]; } }

/* ------------------------------------------------------------ ML-DSA side (as in Stage 1) */
typedef struct { uint8_t pk[CRYPTO_PUBLICKEYBYTES], sk[CRYPTO_SECRETKEYBYTES], rho[SEEDBYTES], tr[TRBYTES], key[SEEDBYTES];
  polyvecl mat[K], s1, s1hat; polyveck s2, s2hat, t0, t0hat, t1, t1shat; } keys_t;
static void keys_init(keys_t *k) { crypto_sign_keypair(k->pk, k->sk); unpack_sk(k->rho, k->tr, k->key, &k->t0, &k->s1, &k->s2, k->sk); unpack_pk(k->rho, &k->t1, k->pk);
  polyvec_matrix_expand(k->mat, k->rho); k->s1hat = k->s1; polyvecl_ntt(&k->s1hat); k->s2hat = k->s2; polyveck_ntt(&k->s2hat);
  k->t0hat = k->t0; polyveck_ntt(&k->t0hat); k->t1shat = k->t1; polyveck_shiftl(&k->t1shat); polyveck_ntt(&k->t1shat); }
static void compute_mu(uint8_t mu[CRHBYTES], const keys_t *k, const uint8_t *m, size_t mlen) { keccak_state st; uint8_t pre[2] = {0, 0};
  shake256_init(&st); shake256_absorb(&st, k->tr, TRBYTES); shake256_absorb(&st, pre, 2); shake256_absorb(&st, m, mlen); shake256_finalize(&st); shake256_squeeze(mu, CRHBYTES, &st); }
static void unpack_ye(polyvecl *y, polyveck *e, const uint8_t *v) { for (int i = 0; i < L; i++) polyz_unpack(&y->vec[i], v + i * POLYZ_PACKEDBYTES);
  const uint8_t *p = v + L * POLYZ_PACKEDBYTES; uint64_t acc = 0; int nb = 0; size_t pos = 0; int32_t B = 1 << (BE_BITS - 1); uint32_t mask = (1u << BE_BITS) - 1;
  for (int i = 0; i < K; i++) for (int j = 0; j < N; j++) { while (nb < BE_BITS) { acc |= (uint64_t)p[pos++] << nb; nb += 8; } e->vec[i].coeffs[j] = B - (int32_t)(acc & mask); acc >>= BE_BITS; nb -= BE_BITS; } }
static void chal(uint8_t ct[CTILDEBYTES], const uint8_t mu[CRHBYTES], const uint8_t *w1pack) { keccak_state st;
  shake256_init(&st); shake256_absorb(&st, mu, CRHBYTES); shake256_absorb(&st, w1pack, K * POLYW1_PACKEDBYTES); shake256_finalize(&st); shake256_squeeze(ct, CTILDEBYTES, &st); }
/* Fin from public data (A, t1, t0, ctilde, z): 0 ok, 1 t0-test, 2 hint weight */
static int fin_public(uint8_t sig[CRYPTO_BYTES], const keys_t *k, const uint8_t ct[CTILDEBYTES], const polyvecl *z) { poly c; polyveck ct0, r, rm, h; polyvecl zhat;
  poly_challenge(&c, ct); poly_ntt(&c); polyveck_pointwise_poly_montgomery(&ct0, &c, &k->t0hat); polyveck_invntt_tomont(&ct0); polyveck_reduce(&ct0);
  if (polyveck_chknorm(&ct0, GAMMA2)) return 1;
  zhat = *z; polyvecl_ntt(&zhat); polyvec_matrix_pointwise_montgomery(&r, k->mat, &zhat); polyveck_pointwise_poly_montgomery(&h, &c, &k->t1shat);
  polyveck_sub(&r, &r, &h); polyveck_reduce(&r); polyveck_invntt_tomont(&r); polyveck_caddq(&r);
  polyveck_sub(&rm, &r, &ct0); polyveck_reduce(&rm); polyveck_caddq(&rm); unsigned wt = 0;
  for (int i = 0; i < K; i++) for (int j = 0; j < N; j++) { int32_t a0, h1 = decompose(&a0, r.vec[i].coeffs[j]), h2 = decompose(&a0, rm.vec[i].coeffs[j]); h.vec[i].coeffs[j] = (h1 != h2); wt += (h1 != h2); }
  if (wt > OMEGA) return 2;
  pack_sig(sig, ct, z, &h); return 0; }

typedef struct { polyvecl y, z; polyveck e, w, w1, w0, u; uint8_t w1pack[K * POLYW1_PACKEDBYTES], ct[CTILDEBYTES]; poly c; int acc; } att_t;
static int eqm(int32_t a, int32_t b) { return ((int64_t)a - b) % Q == 0; }

int main(int argc, char **argv) {
  if (argc < 5) { fprintf(stderr, "usage: %s <sessions> <attempts L> <n> <t> [seed] [logB_W logB_Z]\n", argv[0]); return 1; }
  int nsess = atoi(argv[1]), natt = atoi(argv[2]), n = atoi(argv[3]), t = atoi(argv[4]); uint64_t seed = argc > 5 ? strtoull(argv[5], 0, 10) : 1;
  double logBW = argc > 7 ? atof(argv[6]) : 40, logBZ = argc > 7 ? atof(argv[7]) : 117; /* measured maxima (34.2, 110.7) + 6 bits */
  if (natt < 1 || natt > DEMB || t < 1 || t > n || n > 16) { fprintf(stderr, "bad parameters\n"); return 1; }
  rng_seed(seed); double t0 = now_ms(); setup(); double t_setup = now_ms() - t0;
  CRS = xalloc(NP * sizeof(rp)); for (int i = 0; i < NP; i++) rp_uniform(&CRS[i]);
  DG = xalloc(NP * sizeof(rp)); DG2 = xalloc(NP * sizeof(rp));
  t0 = now_ms(); keygen(&party[PS]); double t_kgS = now_ms() - t0; t0 = now_ms(); keygen(&party[PU]); double t_kgU = now_ms() - t0;
  t0 = now_ms(); rss_share(n, t); double t_rss = now_ms() - t0;
  uint32_t Smask = (1u << t) - 1; /* active set: the first t signers */
  keys_t *key = xalloc(sizeof *key); keys_init(key);
  static int32_t m[NN], dec[NN], tmpc[N]; rp *pt = xalloc(sizeof(rp));
  /* public data: A as plaintexts; committee data: Enc_S(s1), Enc_S(s2) */
  rp *Ahat = xalloc((size_t)K * L * sizeof(rp));
  for (int i = 0; i < K; i++) for (int l = 0; l < L; l++) { poly a = key->mat[i].vec[l]; poly_invntt_tomont(&a); for (int k = 0; k < N; k++) tmpc[k] = montgomery_reduce((int64_t)a.coeffs[k]); emb(m, tmpc); rp_from_i32(&Ahat[i * L + l], m); }
  ct_t *ctS1 = xalloc(L * sizeof(ct_t)), *ctS2 = xalloc(K * sizeof(ct_t));
  for (int l = 0; l < L; l++) { emb(m, key->s1.vec[l].coeffs); rp_from_i32(pt, m); enc1(&ctS1[l], PS, pt); }
  for (int i = 0; i < K; i++) { emb(m, key->s2.vec[i].coeffs); rp_from_i32(pt, m); enc1(&ctS2[i], PS, pt); }
  ct_t *ctY = xalloc((size_t)natt * L * sizeof(ct_t)), *ctE = xalloc((size_t)natt * K * sizeof(ct_t)), *ctw = xalloc((size_t)natt * K * sizeof(ct_t)), *ctW1 = xalloc((size_t)natt * K * sizeof(ct_t));
  ct_t *ctC = xalloc((size_t)natt * sizeof(ct_t)), *ctZ = xalloc((size_t)natt * L * sizeof(ct_t)), *ctU = xalloc((size_t)natt * K * sizeof(ct_t)), *ctF = xalloc((size_t)natt * sizeof(ct_t)), *ctSel = xalloc((size_t)natt * L * sizeof(ct_t));
  ct_t *prod = xalloc(sizeof(ct_t)); ct_t **list = xalloc((size_t)natt * (K + L + 1) * sizeof(ct_t *));
  int32_t (*outW)[N] = xalloc((size_t)natt * K * sizeof *outW), (*outZ)[N] = xalloc((size_t)natt * (L + 1) * sizeof *outZ);
  att_t *at = xalloc((size_t)natt * sizeof *at);
  size_t vb = L * POLYZ_PACKEDBYTES + (size_t)K * N * BE_BITS / 8; uint8_t *v = xalloc(vb), seedb[40] = {0}, mu[CRHBYTES], sig[CRYPTO_BYTES], sig2[CRYPTO_BYTES]; uint64_t ctr = 0;
  memcpy(seedb, &seed, 8);
  double T_in = 0, T_h1 = 0, T_uc = 0, T_h2 = 0, T_h3 = 0, T_h4 = 0, nb_fresh = 0, nb_w = 0, nb_z = 0, nb_u = 0, nb_sel = 0, nfW = 0, nfZ = 0;
  long n_sig_ok = 0, n_fail_none = 0, n_fail_fin = 0, bad_eval = 0, bad_dec = 0, bad_sig = 0, bad_verify = 0, cor_mismatch = 0;
  for (int ss = 0; ss < nsess; ss++) { uint8_t msg[32]; memcpy(msg, &ss, sizeof ss); memset(msg + sizeof ss, 0x5a, sizeof msg - sizeof ss); compute_mu(mu, key, msg, sizeof msg);
    /* ---- ORACLE: bit-level prefix of C1 in the clear (PRF xor user bits, Unpack); result given as fresh two-key ciphertexts */
    for (int j = 0; j < natt; j++) { memcpy(seedb + 32, &ctr, 8); ctr++; shake256(v, vb, seedb, sizeof seedb); unpack_ye(&at[j].y, &at[j].e, v);
      t0 = now_ms();
      for (int l = 0; l < L; l++) { emb(m, at[j].y.vec[l].coeffs); enc2(&ctY[j * L + l], m); }
      for (int i = 0; i < K; i++) { emb(m, at[j].e.vec[i].coeffs); enc2(&ctE[j * K + i], m); }
      T_in += now_ms() - t0;
      /* cleartext reference: w = A y + e, Decompose */
      polyvecl yh = at[j].y; polyvecl_ntt(&yh); polyvec_matrix_pointwise_montgomery(&at[j].w, key->mat, &yh); polyveck_reduce(&at[j].w); polyveck_invntt_tomont(&at[j].w);
      polyveck_add(&at[j].w, &at[j].w, &at[j].e); polyveck_reduce(&at[j].w); polyveck_caddq(&at[j].w); polyveck_decompose(&at[j].w1, &at[j].w0, &at[j].w); polyveck_pack_w1(at[j].w1pack, &at[j].w1); }
    nb_fresh = fmax(nb_fresh, debug_dec(dec, &ctY[0]));
    /* ---- H1 (signers): ct_w = A ct_y + ct_e */
    t0 = now_ms();
    for (int j = 0; j < natt; j++) for (int i = 0; i < K; i++) { ct_t *o = &ctw[j * K + i]; memcpy(o, &ctE[j * K + i], sizeof *o); for (int l = 0; l < L; l++) ct_mulpt_add(o, &ctY[j * L + l], &Ahat[i * L + l]); }
    T_h1 += now_ms() - t0;
    /* ---- ORACLE: HighBits. Check H1 against the reference, then hand W = w1 over as fresh two-key ciphertexts */
    for (int j = 0; j < natt; j++) for (int i = 0; i < K; i++) { nb_w = fmax(nb_w, debug_dec(dec, &ctw[j * K + i]));
      for (int k = 0; k < N; k++) if (!eqm(dec[DEMB * k], at[j].w.vec[i].coeffs[k])) bad_eval++;
      emb(m, at[j].w1.vec[i].coeffs); enc2(&ctW1[j * K + i], m); list[j * K + i] = &ctW1[j * K + i]; }
    /* ---- round 1: threshold decryption of W towards the user */
    noise_fin = 0; thr_decrypt(outW, list, natt * K, logBW, Smask); nfW = fmax(nfW, noise_fin);
    for (int j = 0; j < natt; j++) for (int i = 0; i < K; i++) for (int k = 0; k < N; k++) if (outW[j * K + i][k] != at[j].w1.vec[i].coeffs[k]) bad_dec++;
    /* ---- user: challenges from the decrypted commitments, encrypted under the user key */
    t0 = now_ms();
    for (int j = 0; j < natt; j++) { polyveck w1d; uint8_t pk1[K * POLYW1_PACKEDBYTES]; for (int i = 0; i < K; i++) for (int k = 0; k < N; k++) w1d.vec[i].coeffs[k] = outW[j * K + i][k];
      polyveck_pack_w1(pk1, &w1d); chal(at[j].ct, mu, pk1); poly_challenge(&at[j].c, at[j].ct); emb(m, at[j].c.coeffs); rp_from_i32(pt, m); enc1(&ctC[j], PU, pt); }
    T_uc += now_ms() - t0;
    /* ---- H2, H3 (signers): ct_z = ct_y + ct_c * ct_s1 ; ct_u = ct_w - ct_c * ct_s2 */
    t0 = now_ms(); for (int j = 0; j < natt; j++) for (int l = 0; l < L; l++) { mul_cross(prod, &ctC[j], &ctS1[l]); ct_add(&ctZ[j * L + l], &ctY[j * L + l], prod); } T_h2 += now_ms() - t0;
    t0 = now_ms(); for (int j = 0; j < natt; j++) for (int i = 0; i < K; i++) { mul_cross(prod, &ctC[j], &ctS2[i]); ct_sub(&ctU[j * K + i], &ctw[j * K + i], prod); } T_h3 += now_ms() - t0;
    /* ---- ORACLE: the two norm tests and the selection of the first accepted attempt */
    int j0 = -1;
    for (int j = 0; j < natt; j++) { poly ch = at[j].c; polyveck r, r2; poly_ntt(&ch);
      polyvecl_pointwise_poly_montgomery(&at[j].z, &ch, &key->s1hat); polyvecl_invntt_tomont(&at[j].z); polyvecl_add(&at[j].z, &at[j].z, &at[j].y); polyvecl_reduce(&at[j].z);
      polyveck_pointwise_poly_montgomery(&r, &ch, &key->s2hat); polyveck_invntt_tomont(&r); polyveck_sub(&at[j].u, &at[j].w, &r); polyveck_reduce(&at[j].u); polyveck_caddq(&at[j].u);
      polyveck_sub(&r2, &at[j].w0, &r); polyveck_reduce(&r2); int rok_ref = !polyveck_chknorm(&r2, GAMMA2 - BETA);
      int zok = 1, rok = 1;
      for (int l = 0; l < L; l++) { nb_z = fmax(nb_z, debug_dec(dec, &ctZ[j * L + l])); for (int k = 0; k < N; k++) { if (!eqm(dec[DEMB * k], at[j].z.vec[l].coeffs[k])) bad_eval++; if (abs(dec[DEMB * k]) >= GAMMA1 - BETA) zok = 0; } }
      for (int i = 0; i < K; i++) { nb_u = fmax(nb_u, debug_dec(dec, &ctU[j * K + i])); for (int k = 0; k < N; k++) { if (!eqm(dec[DEMB * k], at[j].u.vec[i].coeffs[k])) bad_eval++;
          int32_t a0, x = dec[DEMB * k]; if (x < 0) x += Q; (void)decompose(&a0, x); if (abs(a0) >= GAMMA2 - BETA) rok = 0; } }
      if (rok != rok_ref) cor_mismatch++; /* Corollary 1: test on LowBits(w - c s2) == test on LowBits(w) - c s2 */
      at[j].acc = zok && rok; int f = at[j].acc && j0 < 0; if (f) j0 = j;
      memset(m, 0, sizeof m); m[0] = f; enc2(&ctF[j], m); }
    /* ---- H4 (signers): ct_Z = ct_f * ct_z */
    t0 = now_ms(); for (int j = 0; j < natt; j++) for (int l = 0; l < L; l++) mul_full(&ctSel[j * L + l], &ctF[j], &ctZ[j * L + l]); T_h4 += now_ms() - t0;
    for (int j = 0; j < natt; j++) for (int l = 0; l < L; l++) { nb_sel = fmax(nb_sel, debug_dec(dec, &ctSel[j * L + l])); list[j * L + l] = &ctSel[j * L + l]; }
    for (int j = 0; j < natt; j++) list[natt * L + j] = &ctF[j];
    /* ---- round 2: threshold decryption of Z = (f_j z_j, f_j) towards the user */
    noise_fin = 0; thr_decrypt(outZ, list, natt * (L + 1), logBZ, Smask); nfZ = fmax(nfZ, noise_fin);
    int ju = -1; for (int j = 0; j < natt; j++) { int f = outZ[natt * L + j][0]; if (f != (j == j0)) bad_dec++; if (f == 1 && ju < 0) ju = j;
      for (int l = 0; l < L; l++) for (int k = 0; k < N; k++) if (outZ[j * L + l][k] != (j == j0 ? at[j].z.vec[l].coeffs[k] : 0)) bad_dec++; }
    /* ---- user: Fin on public data, reference verifier; compare with the cleartext pipeline */
    if (ju < 0) { n_fail_none++; continue; }
    polyvecl zd; for (int l = 0; l < L; l++) for (int k = 0; k < N; k++) zd.vec[l].coeffs[k] = outZ[ju * L + l][k];
    int fr = fin_public(sig, key, at[ju].ct, &zd); if (fr) { n_fail_fin++; continue; }
    if (fin_public(sig2, key, at[j0].ct, &at[j0].z) || memcmp(sig, sig2, CRYPTO_BYTES)) bad_sig++;
    if (crypto_sign_verify(sig, CRYPTO_BYTES, msg, sizeof msg, NULL, 0, key->pk)) bad_verify++; else n_sig_ok++;
  }
  double lgQ = lgpref[NP], polyB = NN * lgQ / 8.0, ns = nsess; int ngW = (natt * K + DEMB - 1) / DEMB, ngZ = (natt * (L + 1) + DEMB - 1) / DEMB;
  printf("# BTHED stage 2a  mode=%d (K=%d,L=%d)  NN=%d  primes=%d  log2Q=%.1f  t=%u  sigma=3.24  (n,t)=(%d,%d)  attempts=%d  sessions=%d  seed=%llu\n", DILITHIUM_MODE, K, L, NN, NP, lgQ, TT, n, t, natt, nsess, (unsigned long long)seed);
  printf("section,item,value\n");
  printf("correctness,signatures_verified,%ld\ncorrectness,sessions_no_accepted_attempt,%ld\ncorrectness,sessions_fin_failed,%ld\n", n_sig_ok, n_fail_none, n_fail_fin);
  printf("correctness,eval_mismatches,%ld\ncorrectness,decryption_mismatches,%ld\ncorrectness,signature_differs_from_cleartext,%ld\ncorrectness,verify_failures,%ld\ncorrectness,corollary1_test_mismatches,%ld\n", bad_eval, bad_dec, bad_sig, bad_verify, cor_mismatch);
  printf("noise_log2,fresh_two_key,%.1f\nnoise_log2,after_H1_w,%.1f\nnoise_log2,after_H2_z,%.1f\nnoise_log2,after_H3_u,%.1f\nnoise_log2,after_H4_fz,%.1f\n", nb_fresh, nb_w, nb_z, nb_u, nb_sel);
  printf("noise_log2,bound_used_W,%.1f\nnoise_log2,bound_used_Z,%.1f\nnoise_log2,after_flooding_W,%.1f\nnoise_log2,after_flooding_Z,%.1f\nnoise_log2,decryption_limit_Q_over_2,%.1f\n", logBW, logBZ, nfW, nfZ, lgQ - 1);
  printf("size_bytes,ring_element,%.0f\nsize_bytes,crs_a,%.0f\nsize_bytes,public_key_b,%.0f\nsize_bytes,evaluation_key_per_party,%.0f\nsize_bytes,fresh_single_key_ct,%.0f\nsize_bytes,two_key_ct,%.0f\n", polyB, NP * polyB, NP * polyB, 3 * NP * polyB, 2 * polyB, 3 * polyB);
  printf("size_bytes,rss_shares_per_signer,%.0f\nsize_bytes,user_to_signers_challenges_per_session,%.0f\nsize_bytes,signer_to_user_round1_per_session,%.0f\nsize_bytes,signer_to_user_round2_per_session,%.0f\n",
         polyB * nsh * (double)(n - t + 1) / n, natt * 2 * polyB, ngW * polyB, ngZ * polyB);
  printf("time_ms,setup_ntt_tables,%.1f\ntime_ms,keygen_committee_incl_evk,%.1f\ntime_ms,keygen_user_incl_evk,%.1f\ntime_ms,rss_share,%.1f\n", t_setup, t_kgS, t_kgU, t_rss);
  printf("time_ms_per_session,oracle_input_encryption,%.1f\ntime_ms_per_session,H1_w_eq_Ay_plus_e,%.1f\ntime_ms_per_session,user_challenges_and_encryption,%.1f\n", T_in / ns, T_h1 / ns, T_uc / ns);
  printf("time_ms_per_session,H2_z_eq_y_plus_c_s1,%.1f\ntime_ms_per_session,H3_u_eq_w_minus_c_s2,%.1f\ntime_ms_per_session,H4_select_f_times_z,%.1f\n", T_h2 / ns, T_h3 / ns, T_h4 / ns);
  printf("time_ms_per_session,partial_decryptions_one_signer,%.1f\ntime_ms_per_session,user_final_decryption,%.1f\ntime_ms,one_partial_decryption,%.2f\n", T_part / ns / t, T_user / ns, T_part / n_part);
  return (bad_eval || bad_dec || bad_sig || bad_verify) ? 2 : 0;
}
