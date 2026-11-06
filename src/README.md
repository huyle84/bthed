# BTHED – Stage 1: cleartext simulator (C)

Stage 1 of the BTHED implementation plan. It runs the **plaintext logic** of the
BTHED signing protocol on top of the pq-crystals ML-DSA reference implementation
(FIPS 204) and checks every produced signature with the **unmodified** reference
ML-DSA verifier.

## What is implemented

| Component (paper) | File / function | Status |
|---|---|---|
| Noisy-commitment signing `NC.Sign_L` (Fig. 1), `w = Ay + e`, low-bits test on `w0 − c s2` | `c1_attempt`, `c2_attempt` | exact |
| `Unpack` (bits → `y ∈ S̃_{γ1}^ℓ`, `e ∈ S̃_{B_e}^k`) | `unpack_ye` | exact |
| XOR-blinding `v_j = PRF(k, H_x(msg1‖j)) ⊕ v^U_j` | `cmd_session` | PRF and `H_x` are SHAKE256 stand-ins |
| User challenge `c̃_j = H(μ‖w1_j)`, `c_j = SampleInBall(c̃_j)` | `chal`, `c2_attempt` | exact |
| Selection of the first successful attempt | `cmd_session` | exact (all L attempts evaluated, as under FHE) |
| User `Finalize` from **public data only** (t0-test, MakeHint on `Az − c t1 2^d`, hint weight) | `fin_public` | exact |
| Verification | `crypto_sign_verify` (reference) | unmodified |
| Replicated secret sharing (reconstruction with 0/1 coefficients, privacy) | `cmd_rss` | exhaustive check |

Three signing rules are compared on the **same random bits** (paired comparison):

* `STD`   – standard ML-DSA attempt (`e = 0`, margin `γ2 − β`);
* `THED`  – noisy commitment, margin `γ2 − β` (THED / Quorus rule);
* `BTHED` – noisy commitment, margin `γ2 − 2β` (this paper).

**Not simulated (Stages 2–3):** FHE encryption and homomorphic evaluation, threshold
decryption with flooding noise, NIZK arguments. The ring-Darkmatter PRF is replaced by
SHAKE256 (it only changes the cost of C1, not the statistics).

## Build and run

```bash
git clone https://github.com/pq-crystals/dilithium.git   # next to this folder
cd src && make                     # builds bthed_sim2 / 3 / 5 (ML-DSA-44 / 65 / 87)
./bthed_sim2 stats 1000000 3       # per-attempt statistics, b_e = 3 (B_e = 4)
./bthed_sim2 session 20000 20 3    # full protocol flow, L = 20 attempts
./bthed_sim2 bench 2001 3          # micro-benchmarks (median ns)
./bthed_sim2 rss 7 4               # RSS checks for n = 7, t = 4
./run_all.sh                       # full batch -> ../results/*.csv
python3 analyze.py                 # -> ../results/RESULTS.md
```

Arguments: `stats <attempts> <b_e> [seed]`, `session <sessions> <L> <b_e> [seed]`,
`bench <reps> <b_e>`, `rss <n> <t>` (n ≤ 20). `b_e ∈ [1, 8]`, `B_e = 2^(b_e − 1)`.
Runs are deterministic given the seed (SHAKE256 stream), except for key generation.

## Files

* `bthed_sim.c` – simulator; `Makefile`; `run_all.sh` – experiment batch;
  `analyze.py` – produces `../results/RESULTS.md` from the CSV files.
