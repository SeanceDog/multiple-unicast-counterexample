#!/usr/bin/env python3
"""Execute the 182-vertex causal code on actual F9 source symbols.

Python 3.10+, standard library only. No precomputed global message vectors
are used. Public interpolation matrices are prepared once; each execution
uses only a vertex's local inputs, public zeros, and received packets.
"""

import argparse
import json
import sys
from dataclasses import asdict, dataclass
from pathlib import Path


# Label a + b*theta by a + 3*b, where theta^2 = -1 in F3.
def add(x, y):
    return (x % 3 + y % 3) % 3 + 3 * ((x // 3 + y // 3) % 3)


def neg(x):
    return (-x % 3) + 3 * ((-(x // 3)) % 3)


def mul(x, y):
    a, b, c, d = x % 3, x // 3, y % 3, y // 3
    return (a * c - b * d) % 3 + 3 * ((a * d + b * c) % 3)


ADD = tuple(tuple(add(x, y) for y in range(9)) for x in range(9))
MUL = tuple(tuple(mul(x, y) for y in range(9)) for x in range(9))
INV = (0,) + tuple(next(y for y in range(1, 9) if MUL[x][y] == 1)
                   for x in range(1, 9))


def power(x, n):
    result = 1
    for _ in range(n):
        result = MUL[result][x]
    return result


def dot(row, values):
    result = 0
    for coefficient, value in zip(row, values):
        result = ADD[result][MUL[coefficient][value]]
    return result


def inverse(matrix):
    """Gauss-Jordan inversion over F9; no floating point arithmetic."""
    n = len(matrix)
    rows = [list(row) + [int(i == j) for j in range(n)]
            for i, row in enumerate(matrix)]
    for j in range(n):
        pivot = next((i for i in range(j, n) if rows[i][j]), None)
        if pivot is None:
            raise ValueError("Singular local interpolation matrix")
        rows[j], rows[pivot] = rows[pivot], rows[j]
        scale = INV[rows[j][j]]
        rows[j] = [MUL[scale][x] for x in rows[j]]
        for i in range(n):
            if i != j and rows[i][j]:
                scale = neg(rows[i][j])
                rows[i] = [ADD[x][MUL[scale][y]]
                           for x, y in zip(rows[i], rows[j])]
    return tuple(tuple(row[n:]) for row in rows)


def require(condition, message):
    # Deliberately not "assert": checks remain enabled with python -O.
    if not condition:
        raise ValueError(message)


@dataclass(frozen=True)
class Port:
    kind: str
    index: int


@dataclass(frozen=True)
class Send:
    edge: int
    receiver: int
    feature: tuple


@dataclass(frozen=True)
class Decode:
    message: int
    feature: tuple


@dataclass(frozen=True)
class Step:
    vertex: int
    known: tuple
    inverse: tuple
    sends: tuple
    decodes: tuple


@dataclass(frozen=True)
class Packet:
    edge: int
    sender: int
    receiver: int
    symbol: int


@dataclass(frozen=True)
class RunResult:
    accepted: bool
    decoded: tuple
    transmissions: int
    mismatches: tuple
    trace: tuple


class CausalCode:
    """Fixed public code; create once and reuse for many independent runs."""

    vertex_count = 182
    message_count = 157

    def __init__(self, parameter_path=None):
        path = (Path(parameter_path) if parameter_path is not None else
                Path(__file__).resolve().parent / "data" / "182-vertex-code.json")
        data = json.loads(path.read_text())
        points = [(0, 0, 1)] + [(0, 1, b) for b in range(9)]
        points += [(1, a, b) for a in range(9) for b in range(9)]
        edges = tuple((i, 91 + j) for i, p in enumerate(points)
                      for j, line in enumerate(points) if dot(p, line) == 0)
        require(len(edges) == 910, "Geometry must have 910 edges")
        incident = [[] for _ in range(182)]
        for e, (u, v) in enumerate(edges):
            incident[u].append(e)
            incident[v].append(e)
        require(all(len(es) == 10 for es in incident), "Incorrect geometry degree")
        order = tuple(data["order"])
        require(sorted(order) == list(range(182)), "Invalid vertex order")
        position = {v: i for i, v in enumerate(order)}
        backward, zeros, unused = (set(data[key]) for key in ("backward", "C", "J"))
        require(len(backward) == 175 and backward <= set(range(910)),
                "Invalid backward edges")
        require(len(zeros) == len(unused) == 9 and zeros.isdisjoint(unused)
                and zeros | unused <= backward, "Invalid C or J edges")
        information = tuple(sorted(backward - zeros - unused))
        message_index = {e: i for i, e in enumerate(information)}
        directions = []
        for e, (u, v) in enumerate(edges):
            earlier, later = sorted((u, v), key=position.__getitem__)
            directions.append((later, earlier) if e in backward else (earlier, later))
        require(all(sum(directions[e][1] == v for e in incident[v]) == 5
                    for v in range(182)), "Orientation is not balanced")
        require(all(sum(directions[e][0] == v for e in zeros)
                    == sum(directions[e][1] == v for e in unused)
                    for v in range(182)), "C/J balance condition fails")
        pairs = tuple((directions[e][1], directions[e][0]) for e in information)

        steps = []
        sent = set()
        decoded = set()
        for v in order:
            # Choose the same two kernel representatives as the proof/verifier.
            kernel = [p for p in points if dot(points[v % 91], p) == 0]
            b0, b1 = kernel[:2]
            coordinates = {
                tuple(ADD[MUL[a][x]][MUL[b][y]] for x, y in zip(b0, b1)): (a, b)
                for a in range(9) for b in range(9)
            }
            feature = {}
            for e in incident[v]:
                u, w = edges[e]
                neighbor = points[w - 91 if v < 91 else u]
                a, b = coordinates[neighbor]
                feature[e] = tuple(power(MUL[power(a, 4 - j)][power(b, j)],
                                         1 if v < 91 else 3) for j in range(5))
            known_edges = [e for e in incident[v] if e in zeros or
                           (directions[e][1] == v and
                            (e in message_index or e not in backward))]
            require(len(known_edges) == 5, f"Vertex {v} does not know five ports")
            known = []
            for e in known_edges:
                if e in zeros:
                    known.append(Port("zero", e))
                elif e in message_index:
                    i = message_index[e]
                    require(pairs[i][0] == v, "Nonlocal source reference")
                    known.append(Port("source", i))
                else:
                    require(e in sent and directions[e][1] == v,
                            "Packet is not available at this vertex")
                    known.append(Port("packet", e))
            sends, decodes = [], []
            for e in incident[v]:
                if e not in backward and directions[e][0] == v:
                    require(e not in sent, "Physical edge sent twice")
                    sends.append(Send(e, directions[e][1], feature[e]))
                    sent.add(e)
                elif e in message_index and directions[e][0] == v:
                    i = message_index[e]
                    require(i not in decoded and pairs[i][1] == v, "Invalid sink")
                    decodes.append(Decode(i, feature[e]))
                    decoded.add(i)
            steps.append(Step(v, tuple(known),
                              inverse([feature[e] for e in known_edges]),
                              tuple(sends), tuple(decodes)))
        require(len(sent) == 735 and decoded == set(range(157)), "Incomplete schedule")
        self.geometry_edges = edges
        self.physical_edges = tuple(sorted(sent))
        self.information_edges = information
        self.pairs = pairs
        self.steps = tuple(steps)

    def run(self, messages, *, capture_trace=False):
        """messages[i] is initially placed only at self.pairs[i][0]."""
        require(isinstance(messages, (list, tuple)) and len(messages) == 157,
                "Input must be a list or tuple of exactly 157 symbols")
        sources = {}
        for i, symbol in enumerate(messages):
            sources.setdefault(self.pairs[i][0], {})[i] = symbol
        return self.run_at_sources(sources, capture_trace=capture_trace)

    def run_at_sources(self, sources, *, capture_trace=False):
        """sources[vertex][message_id] = F9 symbol, with integer keys."""
        require(isinstance(sources, dict), "Sources must be a dictionary")
        source_boxes = [{} for _ in range(182)]
        expected = [None] * 157
        for vertex, inputs in sources.items():
            require(type(vertex) is int and 0 <= vertex < 182, "Invalid source vertex")
            require(isinstance(inputs, dict), "Each source must contain a dictionary")
            for i, symbol in inputs.items():
                require(type(i) is int and 0 <= i < 157, "Invalid message ID")
                require(type(symbol) is int and 0 <= symbol < 9,
                        f"Message {i}: symbol must be an integer from 0 to 8")
                require(self.pairs[i][0] == vertex,
                        f"Message {i} belongs at vertex {self.pairs[i][0]}, not {vertex}")
                require(expected[i] is None, f"Duplicate message {i}")
                source_boxes[vertex][i] = symbol
                expected[i] = symbol
        require(all(x is not None for x in expected), "Some source messages are missing")

        # The network execution never receives the expected output vector.
        decoded, count, trace = self._execute(source_boxes, capture_trace)
        mismatches = tuple(i for i in range(157) if decoded[i] != expected[i])
        return RunResult(not mismatches, decoded, count, mismatches, trace)

    def _execute(self, source_boxes, capture_trace):
        inboxes = [{} for _ in range(182)]
        decoded = [None] * 157
        trace = []
        count = 0
        for step in self.steps:
            v = step.vertex
            values = []
            for port in step.known:
                if port.kind == "zero":
                    values.append(0)
                elif port.kind == "source":
                    values.append(source_boxes[v][port.index])
                else:
                    # Only this vertex's previously received packets are readable.
                    if port.index not in inboxes[v]:
                        raise RuntimeError(f"Vertex {v}: packet {port.index} not received")
                    values.append(inboxes[v].pop(port.index))

            # Interpolate the local quartic from its five known evaluations.
            coefficients = tuple(dot(row, values) for row in step.inverse)
            for operation in step.sends:
                symbol = dot(operation.feature, coefficients)
                inboxes[operation.receiver][operation.edge] = symbol
                count += 1  # Zero-valued packets are charged too.
                if capture_trace:
                    trace.append(Packet(operation.edge, v, operation.receiver, symbol))
            for operation in step.decodes:
                decoded[operation.message] = dot(operation.feature, coefficients)
        if count != 735 or any(x is None for x in decoded) or any(inboxes):
            raise RuntimeError("Incomplete causal execution")
        return tuple(decoded), count, tuple(trace)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", help="JSON list of 157 symbols, or {\"sources\": {...}}")
    parser.add_argument("--decoded", type=Path, help="Save decoded symbols as JSON")
    parser.add_argument("--trace", type=Path, help="Save all 735 actual packet values")
    args = parser.parse_args()
    try:
        data = json.loads(Path(args.input).read_text())
        code = CausalCode()
        if isinstance(data, dict) and "sources" in data:
            sources = {int(v): {int(i): x for i, x in inputs.items()}
                       for v, inputs in data["sources"].items()}
            result = code.run_at_sources(sources, capture_trace=args.trace is not None)
        else:
            messages = data["messages"] if isinstance(data, dict) else data
            result = code.run(messages, capture_trace=args.trace is not None)
        if args.decoded:
            args.decoded.write_text(json.dumps(result.decoded) + "\n")
        if args.trace:
            args.trace.write_text(json.dumps([asdict(p) for p in result.trace], indent=2) + "\n")
    except (OSError, ValueError, TypeError, KeyError, AttributeError, RuntimeError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2
    if result.accepted:
        print(f"ACCEPT: all 157 sink symbols are correct ({result.transmissions} transmissions).")
        return 0
    print(f"REJECT: wrong sink outputs for message IDs {result.mismatches}")
    return 1


if __name__ == "__main__":
    sys.exit(main())
