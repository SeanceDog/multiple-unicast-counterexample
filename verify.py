#!/usr/bin/env python3
"""Verify the 182-vertex causal code using exact, symbolic arithmetic.

Only the Python standard library is needed.  The certificate provides a vertex
order and three edge sets; this verifier independently rebuilds the geometry,
field, local interpolation rules, every packet, and every source--sink distance.
Each source is a distinct formal variable, so a successful run proves correct
decoding for every one of the 9**157 possible input vectors, without sampling.

Usage:
    python3 verify.py
    python3 verify.py path/to/certificate.json --trace transcript.json
"""

import argparse
import collections
import itertools
import json
import sys
from pathlib import Path


class VerificationError(ValueError):
    """A certificate or one of its required mathematical checks is invalid."""


def require(condition, message):
    """Enforce a check even when Python is run with the -O option."""
    if not condition:
        raise VerificationError(message)


# Field elements a + b*theta are encoded by a + 3*b, where a,b are in {0,1,2}
# and theta**2 = -1 = 2.  These are not integers with arithmetic modulo nine.
def add(x, y):
    return (x % 3 + y % 3) % 3 + 3 * ((x // 3 + y // 3) % 3)


def neg(x):
    return (-x % 3) + 3 * ((-(x // 3)) % 3)


def mul(x, y):
    a, b = x % 3, x // 3
    c, d = y % 3, y // 3
    return (a * c + 2 * b * d) % 3 + 3 * ((a * d + b * c) % 3)


ADD = [[add(x, y) for y in range(9)] for x in range(9)]
MUL = [[mul(x, y) for y in range(9)] for x in range(9)]
INV = [0] + [next(y for y in range(1, 9) if MUL[x][y] == 1)
             for x in range(1, 9)]


def power(x, exponent):
    result = 1
    for _ in range(exponent):
        result = MUL[result][x]
    return result


def dot(x, y):
    require(len(x) == len(y), "Internal error: dot-product dimension mismatch")
    result = 0
    for a, b in zip(x, y):
        result = ADD[result][MUL[a][b]]
    return result


def normalized(triple):
    """Choose the projective representative with first nonzero entry one."""
    first = next(x for x in triple if x)
    return tuple(MUL[x][INV[first]] for x in triple)


def inverse(matrix):
    """Invert a square matrix by Gaussian elimination over F_9."""
    size = len(matrix)
    require(size > 0 and all(len(row) == size for row in matrix),
            "Interpolation matrix must be square and nonempty")
    augmented = [row[:] + [int(i == j) for j in range(size)]
                 for i, row in enumerate(matrix)]
    for column in range(size):
        pivot = next((i for i in range(column, size)
                      if augmented[i][column]), None)
        require(pivot is not None, "A local interpolation matrix is singular")
        augmented[column], augmented[pivot] = augmented[pivot], augmented[column]
        scale = INV[augmented[column][column]]
        augmented[column] = [MUL[scale][x] for x in augmented[column]]
        for row in range(size):
            if row != column and augmented[row][column]:
                scale = neg(augmented[row][column])
                augmented[row] = [ADD[x][MUL[scale][y]] for x, y in
                                  zip(augmented[row], augmented[column])]
    return [row[size:] for row in augmented]


def combine(coefficients, rows):
    """Form a linear combination of formal message-coefficient vectors."""
    require(len(coefficients) == len(rows) and len(rows) > 0,
            "Internal error: invalid linear-combination dimensions")
    width = len(rows[0])
    require(all(len(row) == width for row in rows),
            "Internal error: inconsistent vector widths")
    result = [0] * width
    for coefficient, row in zip(coefficients, rows):
        if coefficient:
            result = [ADD[x][MUL[coefficient][y]] for x, y in zip(result, row)]
    return result


def validate_parameters(data):
    """Validate the certificate before converting its lists to sets."""
    require(type(data) is dict, "Certificate must be a JSON object")
    require(set(data) == {"order", "backward", "C", "J"},
            "Certificate keys must be exactly order, backward, C, and J")
    for name, count, bound in (("order", 182, 182), ("backward", 175, 910),
                               ("C", 9, 910), ("J", 9, 910)):
        values = data[name]
        require(type(values) is list, name + " must be a JSON array")
        require(len(values) == count,
                "{} must contain exactly {} entries".format(name, count))
        # bool is a subclass of int in Python; use exact types to reject it.
        require(all(type(value) is int for value in values),
                name + " must contain integers (not booleans or floats)")
        require(all(0 <= value < bound for value in values),
                "{} entries must lie in [0, {})".format(name, bound))
        require(len(set(values)) == count, name + " contains duplicate entries")
    require(set(data["order"]) == set(range(182)),
            "order must be a permutation of the 182 vertex indices")
    backward, zeros, unused = map(set, (data["backward"], data["C"], data["J"]))
    require(zeros.isdisjoint(unused), "C and J must be disjoint")
    require((zeros | unused) <= backward, "C and J must be subsets of backward")


def verify(data, trace_file=None):
    """Verify a certificate and return the summary; optionally save a trace."""
    validate_parameters(data)

    # Lexicographic integer ordering fixes every vertex and edge index.
    points = sorted({normalized(triple)
                     for triple in itertools.product(range(9), repeat=3)
                     if any(triple)})
    require(len(points) == 91, "Geometry must have 91 projective points")
    edges = [(i, 91 + j) for i, point in enumerate(points)
             for j, line in enumerate(points) if dot(point, line) == 0]
    require(len(edges) == 910 and len(set(edges)) == 910,
            "Geometry must have 910 distinct incidence edges")
    incident = [[] for _ in range(182)]
    for edge, (u, v) in enumerate(edges):
        incident[u].append(edge)
        incident[v].append(edge)
    require(all(len(port_list) == 10 for port_list in incident),
            "Each geometry vertex must have degree ten")

    order = data["order"]
    position = {vertex: rank for rank, vertex in enumerate(order)}
    backward, zeros, unused = map(set, (data["backward"], data["C"], data["J"]))
    information = sorted(backward - zeros - unused)
    require(len(information) == 157, "There must be 157 independent messages")

    # The auxiliary orientation points backward on deleted edges, forward on
    # physical edges.  A backward information edge specifies a demand in the
    # opposite direction: the earlier endpoint is the source.
    orientation = []
    for edge, (u, v) in enumerate(edges):
        if (position[u] > position[v]) != (edge in backward):
            u, v = v, u
        orientation.append((u, v))
    for vertex in range(182):
        require(sum(orientation[e][1] == vertex for e in incident[vertex]) == 5,
                "Auxiliary indegree is not five at vertex {}".format(vertex))
        require(sum(orientation[e][0] == vertex for e in zeros) ==
                sum(orientation[e][1] == vertex for e in unused),
                "C/J balance fails at vertex {}".format(vertex))
    pairs = [(orientation[e][1], orientation[e][0]) for e in information]

    # Each local basis consists of the first two normalized neighboring triples.
    # Use their exact coordinates; do not projectively normalize them again.
    features = {}
    for vertex in range(182):
        neighbors = [point for point in points if dot(points[vertex % 91], point) == 0]
        basis0, basis1 = neighbors[:2]
        lookup = {tuple(ADD[MUL[a][x]][MUL[b][y]] for x, y in zip(basis0, basis1)):
                  (a, b) for a, b in itertools.product(range(9), repeat=2)}
        require(len(lookup) == 81, "A local basis is linearly dependent")
        for edge in incident[vertex]:
            other = points[edges[edge][1] - 91 if vertex < 91 else edges[edge][0]]
            require(other in lookup, "A neighboring triple is outside its local plane")
            a, b = lookup[other]
            # Points evaluate Q(a,b); lines evaluate Q(a**3,b**3).
            frobenius = 1 if vertex < 91 else 3
            features[vertex, edge] = [power(MUL[power(a, 4 - j)][power(b, j)],
                                            frobenius) for j in range(5)]

    unit = [[int(i == j) for j in range(157)] for i in range(157)]
    labels = {edge: unit[i] for i, edge in enumerate(information)}
    labels.update({edge: [0] * 157 for edge in zeros})
    source_index = {edge: i for i, edge in enumerate(information)}
    local, packet_index = {}, {}
    packets, decoders = [], []

    # Process vertices chronologically.  Formal vectors encode every packet's
    # dependence on the independent inputs, not specific random test inputs.
    for vertex in order:
        known = [edge for edge in incident[vertex]
                 if edge in zeros or (orientation[edge][1] == vertex and
                                      (edge in information or edge not in backward))]
        require(len(known) == 5 and len(set(known)) == 5,
                "Vertex {} does not have exactly five known ports".format(vertex))
        require(all(edge in labels for edge in known),
                "Vertex {} uses a value before it is available".format(vertex))
        interpolation = inverse([features[vertex, edge] for edge in known])
        values = [labels[edge] for edge in known]
        for edge in incident[vertex]:
            weights = [dot(features[vertex, edge],
                           [interpolation[row][column] for row in range(5)])
                       for column in range(5)]
            local[vertex, edge] = combine(weights, values)
            references = []
            for coefficient, known_edge in zip(weights, known):
                if coefficient == 0 or known_edge in zeros:
                    continue
                if known_edge in source_index:
                    index = source_index[known_edge]
                    require(pairs[index][0] == vertex,
                            "An encoder uses a source message at the wrong vertex")
                    references.append(["source", index, coefficient])
                else:
                    require(known_edge in packet_index and
                            orientation[known_edge][1] == vertex,
                            "An encoder uses a packet not previously received")
                    references.append(["packet", packet_index[known_edge], coefficient])
            if edge not in backward and orientation[edge][0] == vertex:
                require(edge not in labels, "A physical edge is transmitted more than once")
                labels[edge] = local[vertex, edge]
                packet_index[edge] = len(packets)
                packets.append(dict(edge=edge, sender=vertex,
                                    receiver=orientation[edge][1], encoder=references,
                                    vector="".join(map(str, labels[edge]))))
            if edge in source_index and orientation[edge][0] == vertex:
                index = source_index[edge]
                require(local[vertex, edge] == unit[index],
                        "Decoder {} at vertex {} is incorrect".format(index, vertex))
                decoders.append(dict(message=index, vertex=vertex, decoder=references))

    physical = [edge for edge in range(910) if edge not in backward]
    require(len(physical) == 735 and set(packet_index) == set(physical) and len(packets) == 735,
            "Exactly the 735 physical edges must each transmit once")
    require(len(decoders) == 157 and
            {decoder["message"] for decoder in decoders} == set(range(157)),
            "Each of the 157 messages must have exactly one verified decoder")
    # This additionally checks every deleted port, including all nine J ports.
    for edge, (u, v) in enumerate(edges):
        require(local[u, edge] == local[v, edge],
                "The two local values disagree on geometry edge {}".format(edge))

    adjacency = [[] for _ in range(182)]
    for edge in physical:
        u, v = edges[edge]
        adjacency[u].append(v)
        adjacency[v].append(u)
    distances = []
    for index, (source, sink) in enumerate(pairs):
        require(source != sink, "A demand has coincident terminals")
        distance = {source: 0}
        queue = collections.deque([source])
        while queue:
            u = queue.popleft()
            for v in adjacency[u]:
                if v not in distance:
                    distance[v] = distance[u] + 1
                    queue.append(v)
        require(sink in distance, "A demand is disconnected")
        require(distance[sink] == 5,
                "Demand {} does not have distance exactly five".format(index))
        distances.append(distance[sink])
    routing = sum(distances)
    require(routing == 785 and routing - len(packets) == 50,
            "The routing and coding costs are not 785 and 735")
    result = dict(vertices=182, edges=735, messages=157, transmissions=735,
                  distance_histogram=dict(sorted(collections.Counter(distances).items())),
                  routing=routing, saving=routing - len(packets),
                  all_decoders_verified=True, all_encoders_causal=True)
    if trace_file is not None:
        trace = dict(field="F3[z]/(z^2+1); a+3b means a+b*z", projective_points=points,
                     geometry_edges=edges, physical_edge_indices=physical,
                     information_edges=information, pairs=pairs, transmissions=packets,
                     decoders=sorted(decoders, key=lambda item: item["message"]),
                     distances=distances, verification=result)
        Path(trace_file).write_text(json.dumps(trace, separators=(",", ":")) + "\n",
                                    encoding="utf-8")
    return result


def unique_object(pairs):
    """Reject repeated JSON keys rather than silently accepting the last one."""
    result = {}
    for key, value in pairs:
        require(key not in result, "Duplicate JSON key: " + key)
        result[key] = value
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("certificate", nargs="?", type=Path,
                        default=Path(__file__).resolve().parent / "data" / "182-vertex-code.json",
                        help="certificate JSON (default: data/182-vertex-code.json beside this script)")
    parser.add_argument("--trace", type=Path, help="write the complete symbolic transcript to this file")
    args = parser.parse_args()
    try:
        data = json.loads(args.certificate.read_text(encoding="utf-8"),
                          object_pairs_hook=unique_object)
        result = verify(data, args.trace)
    except (VerificationError, OSError, UnicodeError, json.JSONDecodeError) as error:
        print("Verification failed: {}".format(error), file=sys.stderr)
        return 1
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
