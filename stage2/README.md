# BTHED – Stage 2: homomorphic layer (C)

Stage 2a (`bthed_fhe.c`): two-key threshold RLWE layer and arithmetic core. Stage 2b (`bthed_cmp.c`): cost of one
homomorphic comparison, with one key and with two keys.

## Stage 2a

Stage 2a runs the **arithmetic core** of the BTHED signing protocol under a real two-key
(user + committee) homomorphic encryption scheme with threshold decryption towards the user,
on real ML-DSA data (pq-crystals reference code, FIPS 204). Every signature is checked with
the **unmodified** reference verifier and compared byte by byte with the cleartext pipeline.

## What is real

| Component (paper) | Function | Notes |
|---|---|---|
| Two-key scheme, `ct = (c0, cS, cU)`, `c0 + cS·sS + cU·sU = m + t·e` | `keygen`, `enc1`, `enc2` | BGV-type (LSB) encoding, RNS + NTT, plaintext modulus `t = q = 8380417`; the ML-DSA ring is embedded by `Y = X^(N/256)` |
| Multi-key relinearization from individually generated keys | `relin`, `mul_cross`, `mul_full` | CDKS19, Sect. 3 (common random `a`, RNS gadget) |
| Replicated secret sharing of `dkS`, 0/1 reconstruction | `rss_share`, `partdec` | `resp_S(T) = min(S \ T)` |
| Partial decryption with Gaussian flooding, decryption towards the user | `partdec`, `thr_decrypt` | `ν = B·sqrt(2π·N_max·Q_dec/ln 2)`, `Q_dec = 2^64` |
| H1: `w = A·y + e` (plaintext × ciphertext) | `main` | all `L` attempts |
| H2: `z = y + c·s1`, H3: `u = w − c·s2` (user-key × committee-key products) | `mul_cross` | `c` is encrypted under the user key |
| H4: `f·z` (two-key × two-key, depth 2) | `mul_full` | only the selected response decrypts to a non-zero value |
| User: challenges, `Fin` from public data, verification | `chal`, `fin_public`, `crypto_sign_verify` | as in Stage 1 |

## What is emulated in the clear (`ORACLE` in the code)

The non-arithmetic steps of the circuits: the PRF, the XOR with the user's bits, `Unpack`,
`Decompose`/`HighBits`, the two norm tests and the selection of the first accepted attempt.
The harness decrypts the intermediate values with both keys, computes these steps in the clear
and hands the results back as fresh two-key ciphertexts. Under FHE these steps need bootstrapped
comparisons (the CKKS–BFV machinery of THED, ring degree 2^16, GPU); this is Stage 2b/3.
**The timings and sizes below are therefore those of the arithmetic core and of the threshold
layer, not of the full protocol.**

## Parameters

`N = 8192`, `Q` = product of six 31-bit NTT primes (`log2 Q = 186`), ternary secrets, errors
with `σ = 3.24` (within the 128-bit range of the HomomorphicEncryption.org table for `N = 8192`).
No bootstrapping and no modulus switching: the circuit has multiplicative depth 2.

## Limitations

* The noise bounds used for flooding (`2^40` for `W`, `2^117` for `Z`) are the measured maxima
  plus 6 bits, not worst-case bounds as required by the definition in the paper.
* Flooding noise is a rounded Gaussian with 40 bits of precision, not an exact discrete Gaussian.
* The FHE layer takes its randomness from xoshiro256** (not cryptographic). Prototype only.
* One ML-DSA polynomial per ciphertext during evaluation (32 would fit when the same operation is
  applied to all of them); packing is only used before threshold decryption.
* Not constant time, single-threaded, no vector instructions beyond what the compiler emits.

## Build and run

```bash
cd src && make                       # needs the pq-crystals reference code (see Makefile)
./bthed_fhe2 20 8 5 3 1              # <sessions> <attempts L> <n> <t> [seed] [logB_W logB_Z]
./run_all.sh && python3 analyze.py   # -> ../results/*.csv and ../results/RESULTS.md
```

The program returns a non-zero exit code if any homomorphic value, decrypted value or signature
differs from the cleartext reference.

## Stage 2b: one homomorphic comparison

`bthed_cmp.c` is a leveled RNS-CKKS written from scratch (no dependency): canonical-embedding
encoding with `N/2` real slots, 31-bit NTT primes, one special prime for key switching, exact scale
tracking, rescaling, sparse ternary secrets of Hamming weight 256 (as THED). It evaluates
`sign(x) ≈ f3^(df) ∘ g3^(dg)(x)` (Cheon–Kim–Kim, ASIACRYPT 2020; each polynomial has degree 7 and
depth 3, five ciphertext multiplications) on `N/2 = 16384` values at once, and checks every slot.
It runs the same circuit twice: with a single key (standard relinearization) and with two keys
(CDKS19 relinearization from individually generated keys, as in Stage 2a).

```bash
./bthed_cmp 10 5 2 1        # <alpha> <dg> <df> [seed]: inputs are non-zero multiples of 2^-alpha
```

What it does **not** contain: bootstrapping, the "modulus-thrifty" precision schedule of THED,
rotations, and the CKKS–BFV conversions. With a 2^31 scale the fresh encryption error is about
2^-20, so comparisons of 20-bit precision (the low-bits test of ML-DSA) cannot be run here: they
need about 12 compositions (36 levels), a larger scale at the first levels and ring degree 2^16
or a bootstrap. The measured times are single-core CPU times of an unoptimized implementation.
