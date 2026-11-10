/*
 * BTHED - Stage 2b: cost of one homomorphic comparison (CKKS, composite polynomials of
 * Cheon-Kim-Kim, ASIACRYPT 2020), with a single key and with two keys (user + committee).
 *
 * Leveled RNS-CKKS written from scratch: canonical-embedding encoding (N/2 real slots),
 * 31-bit NTT primes, one special prime for key switching, exact scale tracking, rescaling.
 * Single key: standard relinearization. Two keys: CDKS19 relinearization from individually
 * generated keys (common random a), as in Stage 2a, here with the special prime.
 * sign(x) is approximated by f3^(df) o g3^(dg), each a degree-7 odd polynomial (depth 3).
 * No bootstrapping: the number of levels fixes the modulus. Randomness: xoshiro256** (NOT
 * cryptographic); secrets are sparse ternary (Hamming weight 256, as in THED). Prototype only.
 */
#define _POSIX_C_SOURCE 199309L
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include <time.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#ifndef LOGNN
#define LOGNN 15
#endif
#define NN (1 << LOGNN)
#define NH (NN / 2)
#define NPMAX 48
#define HW 256

typedef uint32_t *rp;                        /* (NPQ+1) x NN residues, NTT form unless stated */
#define C(r, i) ((r) + (size_t)(i) * NN)
typedef struct { rp c[3]; int nk, lvl; double scale; } ct_t; /* c0 + c1 s1 (+ c2 s2) */
typedef struct { rp s, *b, *d0, *d1, *d2; } party_t;

static int NPQ; /* Q-primes 0..NPQ-1, special prime at index NPQ */
static uint32_t P[NPMAX + 1], ninv[NPMAX + 1], ninvs[NPMAX + 1], *psi[NPMAX + 1], *psis[NPMAX + 1], *ipsi[NPMAX + 1], *ipsis[NPMAX + 1];
static rp sk1, *ksa, *ksb, *crs, *DG, *DG2, skv[3], T0, T1, T2, T3;
static party_t party[3];
static uint32_t *tmpc, *tmpu;
static double complex *ksi; static int *rot;
static double T_target; static long n_mul;

static double now_ms(void) { struct timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts); return ts.tv_sec * 1e3 + ts.tv_nsec * 1e-6; }
static void *xalloc(size_t n) { void *p = calloc(1, n); if (!p) { fprintf(stderr, "out of memory\n"); exit(1); } return p; }
static rp rp_new(void) { return xalloc((size_t)(NPQ + 1) * NN * 4); }
static uint64_t xs[4];
static inline uint64_t rotl(uint64_t x, int k) { return (x << k) | (x >> (64 - k)); }
static uint64_t rng(void) { uint64_t r = rotl(xs[1] * 5, 7) * 9, t = xs[1] << 17; xs[2] ^= xs[0]; xs[3] ^= xs[1]; xs[1] ^= xs[2]; xs[0] ^= xs[3]; xs[2] ^= t; xs[3] = rotl(xs[3], 45); return r; }
static void rng_seed(uint64_t s) { for (int i = 0; i < 4; i++) { s += 0x9e3779b97f4a7c15ULL; uint64_t z = s; z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL; z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL; xs[i] = z ^ (z >> 31); } }

/* ------------------------------------------------------------ arithmetic */
static uint32_t powmod(uint32_t a, uint64_t e, uint32_t p) { uint64_t r = 1, b = a % p; while (e) { if (e & 1) r = r * b % p; b = b * b % p; e >>= 1; } return (uint32_t)r; }
static int is_prime(uint32_t n) { if (!(n & 1)) return 0; uint32_t d = n - 1; int s = 0; while (!(d & 1)) { d >>= 1; s++; }
  const uint32_t bases[3] = {2, 7, 61};
  for (int i = 0; i < 3; i++) { uint64_t x = powmod(bases[i], d, n); if (x == 1 || x == n - 1) continue; int ok = 0;
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
/* loops over primes 0..lvl and, if ext, the special prime */
#define FORP(lvl, ext) for (int ii_ = 0, i; ii_ <= (lvl) + (ext) && ((i = ii_ <= (lvl) ? ii_ : NPQ), 1); ii_++)
static void p_add(rp r, rp a, rp b, int lvl, int ext) { FORP(lvl, ext) { uint32_t p = P[i], *x = C(r, i), *y = C(a, i), *z = C(b, i); for (int j = 0; j < NN; j++) { uint32_t s = y[j] + z[j]; x[j] = s >= p ? s - p : s; } } }
static void p_sub(rp r, rp a, rp b, int lvl, int ext) { FORP(lvl, ext) { uint32_t p = P[i], *x = C(r, i), *y = C(a, i), *z = C(b, i); for (int j = 0; j < NN; j++) { uint32_t s = y[j] + p - z[j]; x[j] = s >= p ? s - p : s; } } }
static void p_mul(rp r, rp a, rp b, int lvl, int ext) { FORP(lvl, ext) { uint32_t p = P[i], *x = C(r, i), *y = C(a, i), *z = C(b, i); for (int j = 0; j < NN; j++) x[j] = (uint32_t)((uint64_t)y[j] * z[j] % p); } }
static void p_muladd(rp r, rp a, rp b, int lvl, int ext) { FORP(lvl, ext) { uint32_t p = P[i], *x = C(r, i), *y = C(a, i), *z = C(b, i); for (int j = 0; j < NN; j++) x[j] = (uint32_t)(((uint64_t)y[j] * z[j] + x[j]) % p); } }
static void p_zero(rp r) { memset(r, 0, (size_t)(NPQ + 1) * NN * 4); }
static void p_from_i64(rp r, const int64_t *v, int lvl, int ext) { FORP(lvl, ext) { int64_t p = P[i]; uint32_t *x = C(r, i); for (int j = 0; j < NN; j++) { int64_t t = v[j] % p; x[j] = (uint32_t)(t < 0 ? t + p : t); } fntt(x, i); } }
static void p_uniform(rp r, int lvl, int ext) { FORP(lvl, ext) { uint32_t *x = C(r, i); for (int j = 0; j < NN; j++) x[j] = (uint32_t)(rng() % P[i]); } }
static void p_err(rp r, int lvl, int ext) { static int64_t v[NN]; for (int j = 0; j < NN; j++) { uint64_t x = rng(); v[j] = __builtin_popcountll(x & 0x1FFFFF) - __builtin_popcountll((x >> 21) & 0x1FFFFF); } p_from_i64(r, v, lvl, ext); }
static void p_sparse(rp r) { static int64_t v[NN]; memset(v, 0, sizeof v); for (int k = 0; k < HW;) { int j = (int)(rng() % NN); if (!v[j]) { v[j] = (rng() & 1) ? 1 : -1; k++; } } p_from_i64(r, v, NPQ - 1, 1); }
/* x (primes 0..lvl and src, NTT form) -> round(x / P[src]) on primes 0..lvl */
static void div_by_prime(rp x, int src, int lvl) { uint32_t ps = P[src], half = ps / 2; memcpy(tmpc, C(x, src), (size_t)NN * 4); fintt(tmpc, src);
  for (int k = 0; k <= lvl; k++) { uint32_t pk = P[k], pm = ps % pk, inv = powmod(pm, pk - 2, pk), invs = shoup(inv, pk), *y = C(x, k);
    for (int j = 0; j < NN; j++) { uint32_t v = tmpc[j], u = v >= pk ? v - pk : v; if (v > half) { u += pk - pm; if (u >= pk) u -= pk; } tmpu[j] = u; }
    fntt(tmpu, k); for (int j = 0; j < NN; j++) { uint32_t d = y[j] + pk - tmpu[j]; y[j] = mulsh(d >= pk ? d - pk : d, inv, invs, pk); } } }
/* RNS digits of c (level lvl), each extended to primes 0..lvl and the special prime */
static void decomp(rp *dg, rp c, int lvl) { for (int i = 0; i <= lvl; i++) { memcpy(tmpc, C(c, i), (size_t)NN * 4); fintt(tmpc, i);
    for (int kk = 0; kk <= lvl + 1; kk++) { int k = kk <= lvl ? kk : NPQ; uint32_t pk = P[k], *y = C(dg[i], k); for (int j = 0; j < NN; j++) { uint32_t v = tmpc[j]; y[j] = v >= pk ? v - pk : v; } fntt(y, k); } } }

static void setup(int npq) { NPQ = npq; uint32_t cand = (uint32_t)((((1ull << 31) - 1) / (2 * NN)) * (2 * NN) + 1); int n = 0, order[NPMAX + 1];
  order[0] = npq; for (int i = 0; i < npq; i++) order[i + 1] = i;
  while (n <= npq) { if (is_prime(cand)) P[order[n++]] = cand; cand -= 2 * NN; }
  for (int i = 0; i <= npq; i++) { uint32_t p = P[i], g = 2, x; for (;; g++) { x = powmod(g, (p - 1) / (2 * NN), p); if (powmod(x, NN, p) == p - 1) break; }
    uint32_t xi = powmod(x, p - 2, p), *pw = xalloc(4 * NN), *pwi = xalloc(4 * NN); pw[0] = pwi[0] = 1;
    for (int k = 1; k < NN; k++) { pw[k] = (uint32_t)((uint64_t)pw[k - 1] * x % p); pwi[k] = (uint32_t)((uint64_t)pwi[k - 1] * xi % p); }
    psi[i] = xalloc(4 * NN); psis[i] = xalloc(4 * NN); ipsi[i] = xalloc(4 * NN); ipsis[i] = xalloc(4 * NN);
    for (uint32_t k = 0; k < NN; k++) { psi[i][k] = pw[brv(k)]; psis[i][k] = shoup(psi[i][k], p); ipsi[i][k] = pwi[brv(k)]; ipsis[i][k] = shoup(ipsi[i][k], p); }
    free(pw); free(pwi); ninv[i] = powmod(NN, p - 2, p); ninvs[i] = shoup(ninv[i], p); }
  tmpc = xalloc(4 * NN); tmpu = xalloc(4 * NN); T0 = rp_new(); T1 = rp_new(); T2 = rp_new(); T3 = rp_new();
  ksi = xalloc((2 * NN + 1) * sizeof *ksi); rot = xalloc(NH * sizeof *rot);
  for (int k = 0; k <= 2 * NN; k++) ksi[k] = cexp(2 * M_PI * I * k / (2.0 * NN));
  long f = 1; for (int j = 0; j < NH; j++) { rot[j] = (int)f; f = f * 5 % (2 * NN); }
  DG = xalloc(npq * sizeof(rp)); DG2 = xalloc(npq * sizeof(rp)); for (int i = 0; i < npq; i++) { DG[i] = rp_new(); DG2[i] = rp_new(); } }

/* ------------------------------------------------------------ keys */
static void addP(rp x, int i, rp y) { uint32_t p = P[i], pm = P[NPQ] % p, *a = C(x, i), *b = C(y, i); for (int j = 0; j < NN; j++) a[j] = (uint32_t)(((uint64_t)b[j] * pm + a[j]) % p); } /* x += P * g_i * y */
static void keygen1(void) { rp e = rp_new(), s2 = rp_new(); sk1 = rp_new(); p_sparse(sk1); p_mul(s2, sk1, sk1, NPQ - 1, 1);
  ksa = xalloc(NPQ * sizeof(rp)); ksb = xalloc(NPQ * sizeof(rp));
  for (int i = 0; i < NPQ; i++) { ksa[i] = rp_new(); ksb[i] = rp_new(); p_uniform(ksa[i], NPQ - 1, 1); p_err(e, NPQ - 1, 1); p_mul(ksb[i], ksa[i], sk1, NPQ - 1, 1); p_sub(ksb[i], e, ksb[i], NPQ - 1, 1); addP(ksb[i], i, s2); }
  free(e); free(s2); }
static void keygen2(party_t *p) { rp e = rp_new(), r = rp_new(); p->s = rp_new(); p_sparse(p->s); p_sparse(r);
  p->b = xalloc(NPQ * sizeof(rp)); p->d0 = xalloc(NPQ * sizeof(rp)); p->d1 = xalloc(NPQ * sizeof(rp)); p->d2 = xalloc(NPQ * sizeof(rp));
  for (int i = 0; i < NPQ; i++) { p->b[i] = rp_new(); p->d0[i] = rp_new(); p->d1[i] = rp_new(); p->d2[i] = rp_new();
    p_err(e, NPQ - 1, 1); p_mul(p->b[i], p->s, crs[i], NPQ - 1, 1); p_sub(p->b[i], e, p->b[i], NPQ - 1, 1);                       /* b  = -s a + e          */
    p_uniform(p->d1[i], NPQ - 1, 1); p_err(e, NPQ - 1, 1); p_mul(p->d0[i], p->s, p->d1[i], NPQ - 1, 1); p_sub(p->d0[i], e, p->d0[i], NPQ - 1, 1); addP(p->d0[i], i, r); /* d0 = -s d1 + e + P r g */
    p_err(e, NPQ - 1, 1); p_mul(p->d2[i], r, crs[i], NPQ - 1, 1); p_add(p->d2[i], p->d2[i], e, NPQ - 1, 1); addP(p->d2[i], i, p->s); }  /* d2 = r a + e + P s g   */
  free(e); free(r); }

/* ------------------------------------------------------------ ciphertext operations */
static ct_t *ct_new(int nk) { ct_t *c = xalloc(sizeof *c); c->nk = nk; for (int k = 0; k < nk; k++) c->c[k] = rp_new(); return c; }
static void ct_add(ct_t *o, const ct_t *a, const ct_t *b) { if (a->lvl != b->lvl || fabs(a->scale / b->scale - 1) > 1e-6) { fprintf(stderr, "ct_add: level/scale mismatch (%d,%d,%g)\n", a->lvl, b->lvl, a->scale / b->scale); exit(1); }
  for (int k = 0; k < a->nk; k++) p_add(o->c[k], a->c[k], b->c[k], a->lvl, 0);
  o->lvl = a->lvl; o->scale = a->scale; }
static void ct_rescale(ct_t *a) { for (int k = 0; k < a->nk; k++) div_by_prime(a->c[k], a->lvl, a->lvl - 1); a->scale /= P[a->lvl]; a->lvl--; }
static void ct_mulconst(ct_t *o, const ct_t *a, double cst, double ks) { long long K = llround(cst * ks);
  for (int k = 0; k < a->nk; k++) for (int i = 0; i <= a->lvl; i++) { long long p = P[i], m = K % p; if (m < 0) m += p; uint32_t *x = C(o->c[k], i), *y = C(a->c[k], i); for (int j = 0; j < NN; j++) x[j] = (uint32_t)((uint64_t)y[j] * (uint64_t)m % (uint64_t)p); }
  o->lvl = a->lvl; o->scale = a->scale * ks; }
static void relin1(ct_t *o, rp term) { int l = o->lvl; decomp(DG, term, l); p_zero(T0); p_zero(T1);
  for (int d = 0; d <= l; d++) { p_muladd(T0, DG[d], ksb[d], l, 1); p_muladd(T1, DG[d], ksa[d], l, 1); }
  div_by_prime(T0, NPQ, l); div_by_prime(T1, NPQ, l); p_add(o->c[0], o->c[0], T0, l, 0); p_add(o->c[1], o->c[1], T1, l, 0); }
/* two keys: replace term * s_i * s_j (CDKS19, with the special prime) */
static void relin2(ct_t *o, rp term, int i, int j) { int l = o->lvl; decomp(DG, term, l); p_zero(T0); p_zero(T1);
  for (int d = 0; d <= l; d++) { p_muladd(T0, DG[d], party[j].b[d], l, 1); p_muladd(T1, DG[d], party[i].d2[d], l, 1); }
  div_by_prime(T1, NPQ, l); p_add(o->c[j], o->c[j], T1, l, 0);
  div_by_prime(T0, NPQ, l); decomp(DG2, T0, l); p_zero(T0); p_zero(T1);
  for (int d = 0; d <= l; d++) { p_muladd(T0, DG2[d], party[i].d0[d], l, 1); p_muladd(T1, DG2[d], party[i].d1[d], l, 1); }
  div_by_prime(T0, NPQ, l); div_by_prime(T1, NPQ, l); p_add(o->c[0], o->c[0], T0, l, 0); p_add(o->c[i], o->c[i], T1, l, 0); }
static void ct_mul(ct_t *o, const ct_t *a, const ct_t *b) { int l = a->lvl; if (b->lvl != l) { fprintf(stderr, "ct_mul: level mismatch\n"); exit(1); }
  o->lvl = l; o->scale = a->scale * b->scale; n_mul++;
  if (a->nk == 2) { p_mul(o->c[0], a->c[0], b->c[0], l, 0); p_mul(o->c[1], a->c[0], b->c[1], l, 0); p_muladd(o->c[1], a->c[1], b->c[0], l, 0); p_mul(T2, a->c[1], b->c[1], l, 0); relin1(o, T2); }
  else { p_mul(o->c[0], a->c[0], b->c[0], l, 0);
    p_mul(o->c[1], a->c[0], b->c[1], l, 0); p_muladd(o->c[1], a->c[1], b->c[0], l, 0); p_mul(o->c[2], a->c[0], b->c[2], l, 0); p_muladd(o->c[2], a->c[2], b->c[0], l, 0);
    p_mul(T2, a->c[1], b->c[1], l, 0); relin2(o, T2, 1, 1); p_mul(T2, a->c[2], b->c[2], l, 0); relin2(o, T2, 2, 2);
    p_mul(T2, a->c[1], b->c[2], l, 0); p_muladd(T2, a->c[2], b->c[1], l, 0); relin2(o, T2, 1, 2); }
  ct_rescale(o); }

/* ------------------------------------------------------------ encoding, encryption */
static void bitrev_c(double complex *v, int n) { for (int i = 1, j = 0; i < n; i++) { int bit = n >> 1; for (; j >= bit; bit >>= 1) j -= bit; j += bit; if (i < j) { double complex t = v[i]; v[i] = v[j]; v[j] = t; } } }
static void fft_special(double complex *v, int n) { bitrev_c(v, n);
  for (int len = 2; len <= n; len <<= 1) { int lenh = len >> 1, lenq = len << 2, gap = 2 * NN / lenq;
    for (int i = 0; i < n; i += len) for (int j = 0; j < lenh; j++) { double complex u = v[i + j], w = v[i + j + lenh] * ksi[(rot[j] % lenq) * gap]; v[i + j] = u + w; v[i + j + lenh] = u - w; } } }
static void fft_special_inv(double complex *v, int n) {
  for (int len = n; len >= 2; len >>= 1) { int lenh = len >> 1, lenq = len << 2, gap = 2 * NN / lenq;
    for (int i = 0; i < n; i += len) for (int j = 0; j < lenh; j++) { double complex u = v[i + j] + v[i + j + lenh], w = (v[i + j] - v[i + j + lenh]) * ksi[(lenq - rot[j] % lenq) * gap]; v[i + j] = u; v[i + j + lenh] = w; } }
  bitrev_c(v, n); for (int i = 0; i < n; i++) v[i] /= n; }
static void encrypt(ct_t *ct, const double *vals, int lvl, double scale) { static double complex u[NH]; static int64_t m[NN];
  for (int i = 0; i < NH; i++) u[i] = vals[i];
  fft_special_inv(u, NH);
  for (int i = 0; i < NH; i++) { m[i] = llround(creal(u[i]) * scale); m[i + NH] = llround(cimag(u[i]) * scale); }
  ct->lvl = lvl; ct->scale = scale; p_from_i64(ct->c[0], m, lvl, 0); p_err(T0, lvl, 0); p_add(ct->c[0], ct->c[0], T0, lvl, 0);
  for (int k = 1; k < ct->nk; k++) { p_uniform(ct->c[k], lvl, 0); p_mul(T0, ct->c[k], skv[k], lvl, 0); p_sub(ct->c[0], ct->c[0], T0, lvl, 0); } }
static void decrypt(double *vals, const ct_t *ct) { static double complex u[NH]; static double m[NN];
  memcpy(T0, ct->c[0], (size_t)2 * NN * 4); for (int k = 1; k < ct->nk; k++) p_muladd(T0, ct->c[k], skv[k], 1, 0);
  fintt(C(T0, 0), 0); fintt(C(T0, 1), 1); uint64_t q0 = P[0], q1 = P[1], inv = powmod(P[0] % P[1], P[1] - 2, P[1]);
  for (int j = 0; j < NN; j++) { uint64_t v0 = C(T0, 0)[j], v1 = C(T0, 1)[j], t = (v1 + q1 - v0 % q1) % q1 * inv % q1; int64_t x = (int64_t)(v0 + q0 * t); if ((uint64_t)x > q0 * q1 / 2) x -= (int64_t)(q0 * q1); m[j] = (double)x; }
  for (int i = 0; i < NH; i++) u[i] = (m[i] + I * m[i + NH]) / ct->scale;
  fft_special(u, NH);
  for (int i = 0; i < NH; i++) vals[i] = creal(u[i]); }

/* ------------------------------------------------------------ degree-7 odd polynomial, depth 3 */
static ct_t *W[8];
static void eval_odd7(ct_t *x, const double c[4]) { int l = x->lvl; double S = x->scale, q0 = P[l], q1 = P[l - 1], q2 = P[l - 2];
  ct_t *x2 = W[0], *x4 = W[1], *t = W[2], *A = W[3], *B = W[4], *u = W[5];
  ct_mul(x2, x, x); ct_mul(x4, x2, x2); double S2 = x2->scale, S4 = x4->scale, SA = T_target, SB = T_target * q2 / S4;
  ct_mulconst(t, x, c[1], SA * q0 * q1 / (S * S2)); ct_rescale(t); ct_mul(A, t, x2);          /* c3 x^3 */
  ct_mulconst(t, x, c[0], SA * q0 / S); ct_rescale(t); t->lvl--; ct_add(A, A, t);            /* + c1 x  */
  ct_mulconst(t, x, c[3], SB * q0 * q1 / (S * S2)); ct_rescale(t); ct_mul(B, t, x2);          /* c7 x^3 */
  ct_mulconst(t, x, c[2], SB * q0 / S); ct_rescale(t); t->lvl--; ct_add(B, B, t);            /* + c5 x  */
  ct_mul(u, x4, B); A->lvl--; ct_add(x, u, A); }                                              /* A + x^4 B */
static double poly7(double x, const double c[4]) { double y = x * x; return x * (c[0] + y * (c[1] + y * (c[2] + y * c[3]))); }

int main(int argc, char **argv) {
  if (argc < 4) { fprintf(stderr, "usage: %s <alpha: inputs are multiples of 2^-alpha> <dg> <df> [seed]\n", argv[0]); return 1; }
  int alpha = atoi(argv[1]), dg = atoi(argv[2]), df = atoi(argv[3]); uint64_t seed = argc > 4 ? strtoull(argv[4], 0, 10) : 1; int ncomp = dg + df, npq = 2 + 3 * ncomp;
  if (npq > NPMAX || alpha < 1 || alpha > 30) { fprintf(stderr, "bad parameters\n"); return 1; }
  const double F3[4] = {35.0 / 16, -35.0 / 16, 21.0 / 16, -5.0 / 16}, G3[4] = {4589.0 / 1024, -16577.0 / 1024, 25614.0 / 1024, -12860.0 / 1024};
  rng_seed(seed); double t0 = now_ms(); setup(npq); double t_setup = now_ms() - t0; T_target = (double)P[npq / 2];
  t0 = now_ms(); keygen1(); double t_kg1 = now_ms() - t0;
  crs = xalloc(npq * sizeof(rp)); for (int i = 0; i < npq; i++) { crs[i] = rp_new(); p_uniform(crs[i], npq - 1, 1); }
  t0 = now_ms(); keygen2(&party[1]); keygen2(&party[2]); double t_kg2 = (now_ms() - t0) / 2;
  static double x[NH], ref[NH], out[NH]; double lgQ = 0, lgP = log2((double)P[npq]); for (int i = 0; i < npq; i++) lgQ += log2((double)P[i]);
  for (int i = 0; i < NH; i++) { long d = (long)(rng() % ((1ul << alpha))) + 1; if (rng() & 1) d = -d; if (i < 2) d = i ? 1 : -1; x[i] = ldexp((double)d, -alpha);
    double v = x[i]; for (int k = 0; k < dg; k++) v = poly7(v, G3); for (int k = 0; k < df; k++) v = poly7(v, F3); ref[i] = v; }
  double errp = 0; for (int i = 0; i < NH; i++) errp = fmax(errp, fabs(ref[i] - (x[i] > 0 ? 1 : -1)));
  printf("# BTHED stage 2b  N=%d  slots=%d  Q-primes=%d (+1 special)  log2Q=%.1f  log2PQ=%.1f  scale=2^%.2f  secret HW=%d  sign ~ f3^%d o g3^%d  inputs: multiples of 2^-%d  seed=%llu\n",
         NN, NH, npq, lgQ, lgQ + lgP, log2(T_target), HW, df, dg, alpha, (unsigned long long)seed);
  printf("variant,item,value\n");
  printf("cleartext,max_error_of_the_polynomial_vs_sign,%.3g\n", errp);
  double polyB = NN * (lgQ + lgP) / 8, ttot[4] = {0}, tmul[4] = {0};
  for (int nk = 2; nk <= 3; nk++) { const char *nm = nk == 2 ? "single_key" : "two_key";
    if (nk == 2) skv[1] = sk1; else { skv[1] = party[1].s; skv[2] = party[2].s; }
    ct_t *ct = ct_new(nk); for (int k = 0; k < 6; k++) { if (W[k]) { for (int c = 0; c < W[k]->nk; c++) free(W[k]->c[c]); free(W[k]); } W[k] = ct_new(nk); }
    encrypt(ct, x, npq - 1, T_target); decrypt(out, ct); double e0 = 0; for (int i = 0; i < NH; i++) e0 = fmax(e0, fabs(out[i] - x[i]));
    /* one multiplication at the top level, for the ratio */
    t0 = now_ms(); ct_mul(W[0], ct, ct); tmul[nk] = now_ms() - t0;
    n_mul = 0; t0 = now_ms();
    for (int k = 0; k < dg; k++) eval_odd7(ct, G3);
    for (int k = 0; k < df; k++) eval_odd7(ct, F3);
    ttot[nk] = now_ms() - t0; decrypt(out, ct);
    double es = 0, er = 0; long wrong = 0; for (int i = 0; i < NH; i++) { double sg = x[i] > 0 ? 1 : -1; es = fmax(es, fabs(out[i] - sg)); er = fmax(er, fabs(out[i] - ref[i])); if (out[i] * sg <= 0.5) wrong++; }
    printf("%s,fresh_encryption_error_log2,%.1f\n%s,levels_consumed,%d\n%s,final_level,%d\n%s,ciphertext_multiplications,%ld\n", nm, log2(e0), nm, 3 * ncomp, nm, ct->lvl, nm, n_mul);
    printf("%s,max_error_vs_sign,%.3g\n%s,max_error_vs_cleartext_polynomial_log2,%.1f\n%s,wrong_signs_out_of_%d,%ld\n", nm, es, nm, log2(er), nm, NH, wrong);
    printf("%s,time_one_multiplication_top_level_ms,%.1f\n%s,time_comparison_total_ms,%.1f\n%s,time_per_comparison_amortized_us,%.2f\n", nm, tmul[nk], nm, ttot[nk], nm, ttot[nk] * 1e3 / NH);
    printf("%s,ciphertext_top_level_bytes,%.0f\n%s,evaluation_key_bytes,%.0f\n", nm, nk * NN * lgQ / 8, nm, (nk == 2 ? 2.0 : 4.0) * npq * polyB);
    for (int c = 0; c < nk; c++) free(ct->c[c]);
    free(ct); }
  printf("ratio,two_key_over_single_key_one_multiplication,%.2f\nratio,two_key_over_single_key_comparison,%.2f\n", tmul[3] / tmul[2], ttot[3] / ttot[2]);
  printf("time_ms,setup,%.0f\ntime_ms,keygen_single_key,%.0f\ntime_ms,keygen_per_party_two_key,%.0f\n", t_setup, t_kg1, t_kg2);
  return 0;
}
