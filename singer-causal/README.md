# A lossless causal schedule for the PG(2,9) construction

This folder gives a causal realization of the 166-session, 744-edge parameter point in §3 of [Braverman–He](https://arxiv.org/abs/2610.10108v1). **The ratio 415/372 already appears in their paper.** The addition is a causal schedule with no zero controls, together with exact routing certificates. The local rules are exactly those of their §3–§4; this is a small supplement to their construction.

Label points P_i and lines L_i modulo 91, with P_i adjacent to L_(i+d) for d in {3,4,15,18,28,51,60,68,88,90}. This is a perfect difference set, giving PG(2,9). Orient the first five offsets point→line and the last five line→point, and use the order P_0, L_0, P_1, L_1, …, P_90, L_90. Delete the edges running backward: their counts are 3+4+15+18+28+40+31+23+3+1 = 166. This leaves 744 physical edges. Each deleted edge becomes a session from its earlier endpoint to its later endpoint; there are no zero-control or ignored edges.

For causality, restrict the paper's compatible assignment space W to the deleted edges. This restriction is injective: if those values vanish, induction through the order gives five known zero values at each vertex; its quartic interpolation forces every remaining value to vanish. Since §3 gives dim W ≥ 166, restriction to the 166 deleted values is an isomorphism. Arbitrary source symbols therefore extend uniquely, and interpolation through the order computes and decodes that extension causally.

All rates below are common rates per independent session, with unit edge capacity shared between directions in the first row and unit capacity in each direction in the other rows. Coding has exact rate 1, with a matching source-cut bound.

| Capacity and traffic | Exact routing optimum | Coding/routing |
| --- | --- | --- |
| Shared capacity, 166 sessions | 372/415 | 415/372 |
| Per-direction capacity, one-way traffic | 3025/3068 | 3068/3025 |
| Per-direction capacity, paired with reversal: 332 sessions | 372/415 | 415/372 |

The one-way result is unchanged with unlimited capacity on arcs against the causal order: all such arcs have length zero in its routing certificate. In contrast, the published 157-session code has gap 1 under per-direction capacity. Pairing the Singer code with its reversal is an instance of known linear-code reciprocity; the verifier checks both directions explicitly and uses a shared-capacity routing flow and its reversal to attain the paired optimum.

Before deletion, every vertex has in = out = 5, and the absence of controls makes the retained-arc and demand divergences equal. The potential ψ(P_i) = floor(i/4), ψ(L_i) = floor((i+1)/4) − 1/2 gives lengths max(0, 1/2 + ψ(v) − ψ(u)); it adds 2,653 to both the half-hop cost 372 and distance total 415, giving 3025 and 3068 while conserving the margin 43.

Every local map satisfies MᵀM = −I over F9 (local self-duality); the paper's orthogonality identity is also the residue theorem applied to the product of two local quartics.

[Zhang–Li–Li, arXiv:2610.09367v1](https://arxiv.org/abs/2610.09367v1), give a different construction; the reviewed version does not contain this compact Singer certificate.

From the repository root, run:

```sh
python singer-causal/verify.py
```

Only the Python standard library is required. The checker rebuilds the geometry and local rules, verifies every input symbolically, checks the exact rational flows and integer metrics in [certificates.json](certificates.json), and rejects two deliberately corrupted controls. It completes in under one minute; run without Python's optimization flag. The JSON stores the order and deleted-edge IDs in the paper's canonical projective-coordinate indexing; the checker derives their Singer labels independently.

I'm not a mathematician; I'm currently studying the shape of AI intelligence, and this popped out along the way. Found with AI research agents; checked by exact computation and a separately written checker, but not checked by a proof assistant. Apologies if you already have this.

This folder is offered under the repository's existing license terms, if any; no separate license is added.
