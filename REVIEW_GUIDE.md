# External review guide

The most useful review separates machine-checkable finite arithmetic from the
mathematical interfaces surrounding it.

## 1. Artifact integrity

```bash
./scripts/check_artifacts.sh
```

Confirm both release archives, the paper, the transparent JSON certificate,
and the public moment files match `ARTIFACT_SHA256SUMS`.

## 2. Independent exact replay

```bash
python3 -m pip install -r requirements-replay.txt
./scripts/run_full_replay.sh
```

The script extracts the self-contained bundle into a fresh temporary
directory, verifies its internal manifest, reconstructs the moment vector,
and replays both JSON and pickle certificates.

## 3. Finite-certificate audit

Check independently that:

- descriptor parsing rejects malformed and duplicate rows;
- all multipliers and both envelope legs have the required signs;
- the two legs sum exactly to the common denominator;
- every one of 107 K7 and 410 K8 roots satisfies its multiplier lower bound;
- K7, K8, and Horn functionals are reconstructed from exact combinatorial
  data rather than trusted floating matrices;
- the digit-split sparse product cannot overflow signed 64-bit arithmetic;
- representative columns agree with direct Python-integer dot products;
- all 12,172 rational residuals are nonnegative;
- the stored objective equals the reconstructed exact objective and is
  strictly negative.

## 4. Mathematical bridge

Audit the following claims against the manuscript and cited sources:

1. The normalization is
   \(d_{\rm edge}=2|E(G)|/N^2\) and
   \(d_{\rm mono}=2\operatorname{bip}(G)/N^2\).
2. The certificate covers the *closed* band
   \([0.2486,0.3197]\).
3. Every selected per-root K7/K8 row is valid for every genuine
   triangle-free graphon.
4. The rooted-Horn rows follow from the stated copositivity /
   Motzkin-Straus inequality.
5. The fixed Gram functional is manifestly positive semidefinite.
6. The dual sign convention gives
   \(d_{\rm mono}(W)\le2/25+\delta_\star\).
7. Ferudun's blow-up identity is exact:
   \(\operatorname{bip}(G[t])=t^2\operatorname{bip}(G)\).
8. The Balogh-Clemen-Lidický tail theorem transfers through blow-ups for
   densities strictly outside the band, while the certificate covers both
   boundary values.

## 5. Reporting

Please open a GitHub issue with:

- operating system and Python version;
- the exact commit SHA;
- whether artifact checks, moment replay, JSON replay, and pickle replay pass;
- the first failing command and complete terminal output for any failure;
- any mathematical objection identified by section and equation number.
