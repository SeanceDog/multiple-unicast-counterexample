# A 182-vertex causal network coding counterexample

This repository contains the complete finite specification, an exact symbolic
verifier, and an executable simulator for a causal linear network code over
$\mathbb F_9$ on a simple undirected graph.

| Quantity | Value |
| --- | ---: |
| Vertices | 182 |
| Physical edges | 735 |
| Independent unicast messages | 157 field symbols |
| Transmissions per execution | 735 |
| Distance of every source–sink pair | 5 |
| Shortest-path routing cost | 785 |
| Transmission saving | 50 |

Each physical edge carries exactly one symbol, including when that symbol is
zero. Each encoder uses only its vertex's source messages and previously received
packets. Different sessions may share terminal vertices.

With unit capacity per undirected edge, shared between the two directions,
pipelining achieves common coding rate at least $1$. Any fractional routing
has common rate at most $735/785=147/157<1$. This is an upper bound on flow
throughput, not a claim that its exact optimum has been computed.

The [construction](CONSTRUCTION.md) defines the field, graph, local codes,
source–sink pairs, chronological schedule, and the cost and throughput comparison.

## Quick start

Requires **Python 3.10 or later**, using only its standard library. No installation,
third-party packages, network access, or search procedure is required.
Run the following commands from this repository's root directory.

### Verify the complete mathematical certificate

```sh
python3 verify.py
```

Expected output:

```json
{
  "vertices": 182,
  "edges": 735,
  "messages": 157,
  "transmissions": 735,
  "distance_histogram": {
    "5": 157
  },
  "routing": 785,
  "saving": 50,
  "all_decoders_verified": true,
  "all_encoders_causal": true
}
```

The verifier reconstructs the geometry, field arithmetic, local interpolation
matrices, and communication schedule from the four lists in
[`data/182-vertex-code.json`](data/182-vertex-code.json). It tracks each packet as
an exact vector of coefficients in the 157 independent source variables and
checks that every sink recovers the corresponding unit vector. Thus it verifies
**all $9^{157}$ possible source inputs by linearity**; it does not enumerate or
sample those inputs. It also checks causal availability and all terminal distances.
It does not read the stored transcript or assume an analytic dimension bound.
Checks remain active under `python3 -O`.

To regenerate the complete transcript of local encoding and decoding recipes:

```sh
python3 verify.py --trace regenerated-transcript.json
```

The verifier exits with status zero on success and a nonzero status on failure.
An optional positional argument selects another parameter file.

### Execute the code on actual messages

```sh
python3 causal_code.py examples/example_messages.json
```

Expected output:

```text
ACCEPT: all 157 sink symbols are correct (735 transmissions).
```

The input is a JSON list of 157 integers in `0,...,8`. Entry `i` is message `i`.
The integers label field elements, as specified below. Message locations are in
[`data/message_locations.json`](data/message_locations.json).

An alternative input format groups messages explicitly by their source vertex:

```sh
python3 causal_code.py examples/example_sources.json --decoded decoded.json --trace packets.json
```

Its format is `{"sources": {"vertex": {"message_id": symbol, ...}, ...}}`.
Both supplied example files contain the same messages. Incorrect owners, missing
messages, and invalid symbols are rejected. `decoded.json` records the 157 sink
outputs; `packets.json` records all 735 actual transmitted values in order.
Runtime exit codes are 0 for ACCEPT, 1 for incorrect sink outputs, and 2 for
input or execution errors.

The simulator prepares the fixed public interpolation matrices once. During
execution, each vertex reads only its own source inputs, public zero ports, and
its previously received packets. It never reads the reference transcript or
precomputed global coefficient vectors. Only the testing harness compares sink
outputs with the expected input vector, after execution has finished.

### Random testing

```sh
python3 random_tester.py --seed 182 --cases 10000
```

The tester generates fresh independent uniform field symbols and executes the
whole protocol for each input, printing a counter every 1,000 successful cases.
Run `python3 random_tester.py` without arguments to continue until Ctrl+C.
`--report-every` changes the reporting interval. A failing input is saved as
`failure.json` and can be replayed with `python3 causal_code.py failure.json`.
Random testing exercises the implementation; the symbolic verifier supplies the
universal correctness check.

### Run the tests

```sh
python3 -m unittest discover -s tests -v
```

The tests cover actual execution, a separate replay of the reference recipes,
all 1,256 nonzero single-source inputs, random and constant inputs, field
arithmetic, terminal distances, ownership, unavailable packets, and corruption
of an encoder. Verifier tests also reject malformed parameters and check that
verification remains enabled in optimized Python.

The included GitHub Actions workflow runs the symbolic verifier and tests on
Python 3.10 and 3.12 when pushed or invoked manually.

## Field and indexing conventions

The field is

$$
\mathbb F_9=\mathbb F_3[\theta]/(\theta^2+1).
$$

Integer `a + 3*b` denotes $a+b\theta$, for $a,b\in\{0,1,2\}$.
For example, the square of label `3` has label `2`.
**These integers are labels; arithmetic is not performed modulo nine.**
All calculations use exact finite-field operations without floating point.

Vertex, edge, message, and packet identifiers start at zero. Candidate edges
are indexed `0,...,909`; only 735 of these are physical edges. The data field
`backward` is the set $R$ in the mathematical description. The other fields
are the vertex permutation `order` ($\pi$) and the two nine-element sets `C`
and `J`. Message `i` corresponds to entry `i` of the sorted list
$B=R\setminus(C\cup J)$, with its source at the earlier endpoint in $\pi$.

## Python interface

```python
from causal_code import CausalCode

code = CausalCode()
messages = [i % 9 for i in range(code.message_count)]
result = code.run(messages, capture_trace=True)

print("ACCEPT" if result.accepted else "REJECT")
print(result.decoded)
print(result.transmissions)  # 735
```

`code.pairs[i]` gives the source and sink of message `i`. For distributed input,
use `code.run_at_sources({vertex: {message_id: symbol, ...}, ...})` with integer
dictionary keys. Every call allocates fresh communication state.

## Repository contents

| File | Purpose |
| --- | --- |
| `CONSTRUCTION.md` | Mathematical definition and justification |
| `verify.py` | Independent exact symbolic verification |
| `causal_code.py` | Actual-symbol causal simulator |
| `random_tester.py` | Reproducible random-input testing |
| `data/182-vertex-code.json` | Complete finite parameters; the essential certificate |
| `data/reference-transcript.json` | Reconstructed local recipes and coefficient vectors |
| `data/verification.json` | Expected verification summary |
| `data/message_locations.json` | Message-to-terminal mapping |
| `examples/` | Inputs in flat and source-grouped formats |
| `tests/` | Runtime and verifier tests |
| `.github/workflows/verify.yml` | Automated verification on GitHub |

