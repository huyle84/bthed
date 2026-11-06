# bthed-stage1

Stage-1 (cleartext) simulator of **BTHED – Blind Threshold Dilithium from FHE**:
the plaintext logic of the BTHED signing protocol (noisy-commitment ML-DSA signing
with doubled low-bits margin, XOR-blinded signing randomness, user-side
finalization from public data, replicated secret sharing), built on the
[pq-crystals ML-DSA reference implementation](https://github.com/pq-crystals/dilithium)
(FIPS 204). Every produced signature is checked with the unmodified reference verifier.

FHE evaluation, threshold decryption and NIZK arguments are **not** included (Stages 2–3).

```bash
git clone https://github.com/pq-crystals/dilithium.git
cd src && make && ./run_all.sh && python3 analyze.py   # ~10 min, writes ../results/
```

* `src/` – simulator (`bthed_sim.c`), `Makefile`, experiment batch (`run_all.sh`),
  analysis (`analyze.py`); see `src/README.md` for details and command-line options.
* `results/` – raw CSV outputs and the generated report `results/RESULTS.md`.

The pq-crystals code is not redistributed here; it is fetched by the user (see above).
License: MIT (this repository only).
