"""Exact standard-library verification of the finite Singer causal construction.

Run: python verify.py
The graph, field arithmetic, local quartics, order and linear maps are rebuilt.
The certificate file supplies rational paths and consistency checks for derived
orders, edge sets, metrics and the published appendix parameters.
No files are written and no external packages are imported.
"""
from collections import Counter, deque
from fractions import Fraction
from heapq import heappop, heappush
from itertools import product
from pathlib import Path
import copy
import json
import time

ADD = [[(a % 3 + b % 3) % 3 + 3 * ((a // 3 + b // 3) % 3)
        for b in range(9)] for a in range(9)]
MUL = [[((a % 3) * (b % 3) + 2 * (a // 3) * (b // 3)) % 3
        + 3 * (((a % 3) * (b // 3) + (a // 3) * (b % 3)) % 3)
        for b in range(9)] for a in range(9)]
NEG = [(-a % 3) % 3 + 3 * (-(a // 3) % 3) for a in range(9)]
INV = [0] + [next(b for b in range(1, 9) if MUL[a][b] == 1)
             for a in range(1, 9)]
BH = {"order":[93,23,16,143,107,54,88,38,157,40,171,124,123,118,90,32,15,146,114,179,181,39,128,0,119,83,21,6,132,160,17,79,22,42,41,174,5,167,26,96,8,61,168,101,75,43,102,127,99,12,30,57,147,92,138,159,113,31,172,7,116,78,28,126,19,120,89,117,9,14,122,56,72,115,153,47,65,103,74,137,97,10,110,154,135,95,4,125,86,139,34,162,80,87,129,27,134,140,36,13,148,149,82,37,176,91,73,1,51,141,105,109,35,44,130,66,11,98,121,25,76,169,49,53,150,2,104,29,158,69,180,81,52,106,112,131,164,63,71,170,155,24,64,156,62,48,46,145,163,45,77,175,33,142,18,178,166,144,152,55,67,70,177,85,111,20,136,84,100,58,151,68,50,133,165,94,60,173,3,59,161,108],"deleted":[7,8,9,30,31,38,39,56,57,60,68,75,79,84,88,94,102,103,117,122,127,128,141,151,154,155,157,163,164,166,167,168,176,179,183,187,203,207,208,212,215,221,224,225,231,234,236,237,238,242,243,251,253,266,282,290,300,309,318,321,324,325,330,333,338,341,349,380,386,387,389,396,398,400,401,406,407,409,411,412,418,425,426,428,439,444,454,458,464,465,468,470,476,480,481,492,502,503,504,513,517,521,522,527,539,541,544,548,552,557,559,569,579,583,585,588,592,593,598,602,603,606,607,612,619,620,632,634,643,646,670,672,673,674,682,685,688,694,695,701,702,703,709,716,718,737,746,754,757,759,773,790,821,835,837,838,839,841,842,843,846,850,853,856,857,861,862,882,884,885,887,900,901,904,907],"constants":[102,163,164,341,386,406,476,513,862],"ignored":[103,349,470,517,603,606,674,716,861]}

def total(values):
    out = 0
    for value in values:
        out = ADD[out][value]
    return out

def power(a, n):
    out = 1
    for _ in range(n):
        out = MUL[out][a]
    return out

def dot(a, b):
    return total(MUL[x][y] for x, y in zip(a, b))

def matmul(a, b):
    columns = list(zip(*b))
    return [[dot(row, col) for col in columns] for row in a]

def inverse(a):
    n = len(a)
    rows = [list(row) + [int(i == j) for j in range(n)]
            for i, row in enumerate(a)]
    for col in range(n):
        pivot = next((i for i in range(col, n) if rows[i][col]), None)
        assert pivot is not None, "Singular interpolation matrix"
        rows[col], rows[pivot] = rows[pivot], rows[col]
        scale = INV[rows[col][col]]
        rows[col] = [MUL[scale][x] for x in rows[col]]
        for i in range(n):
            if i != col and rows[i][col]:
                scale = NEG[rows[i][col]]
                rows[i] = [ADD[x][MUL[scale][y]]
                           for x, y in zip(rows[i], rows[col])]
    return [row[n:] for row in rows]

def normalize(v):
    scale = INV[next(x for x in v if x)]
    return tuple(MUL[scale][x] for x in v)

def reconstruct():
    triples = [(0, 0, 1)] + [(0, 1, a) for a in range(9)]
    triples += [(1, a, b) for a in range(9) for b in range(9)]
    edges = [(i, 91 + j) for i, p in enumerate(triples)
             for j, line in enumerate(triples) if dot(p, line) == 0]
    incident = [[] for _ in range(182)]
    for q, (u, v) in enumerate(edges):
        incident[u].append(q)
        incident[v].append(q)
    assert len(edges) == 910 and all(len(x) == 10 for x in incident)
    # Derive coordinates in the first two neighboring normalized triples.
    rows = []
    for v, qs in enumerate(incident):
        neighbors = [triples[(sum(edges[q]) - v) % 91] for q in qs]
        u, w = neighbors[:2]
        coordinates = {
            tuple(ADD[MUL[a][x]][MUL[b][y]] for x, y in zip(u, w)): (a, b)
            for a, b in product(range(9), repeat=2)
        }
        local = []
        for y in neighbors:
            a, b = coordinates[y]
            if v >= 91:
                a, b = power(a, 3), power(b, 3)
            local.append([MUL[power(a, 4 - r)][power(b, r)]
                          for r in range(5)])
        assert matmul(list(zip(*local)), local) == [[0] * 5 for _ in range(5)]
        rows.append(local)
    # Companion x^3 - x - 4 and its inverse-transpose action.
    assert all(ADD[ADD[power(x, 3)][NEG[x]]][NEG[4]] for x in range(9))
    points, lines = [], []
    p = line = (1, 0, 0)
    for _ in range(91):
        points.append(p)
        lines.append(line)
        x, y, z = p
        p = normalize((MUL[4][z], ADD[x][z], y))
        x, y, z = line
        line = normalize((MUL[INV[4]][ADD[z][NEG[x]]], x, y))
    assert p == line == (1, 0, 0)
    assert len(set(points)) == len(set(lines)) == 91
    lookup = {v: i for i, v in enumerate(triples)}
    ps = [lookup[points[5 * i % 91]] for i in range(91)]
    ls = [91 + lookup[lines[5 * (i + 21) % 91]] for i in range(91)]
    order = [v for pair in zip(ps, ls) for v in pair]
    position = {v: i for i, v in enumerate(order)}
    edge_id = {edge: q for q, edge in enumerate(edges)}
    shifts = [3, 4, 15, 18, 28, 51, 60, 68, 88, 90]
    oriented = {}
    deleted_by_shift = {shift: 0 for shift in shifts}
    assert Counter((a - b) % 91 for a in shifts for b in shifts if a != b) == Counter(range(1, 91))
    incoming = [0] * 182
    outgoing = [0] * 182
    for i in range(91):
        for shift in shifts:
            u, v = ps[i], ls[(i + shift) % 91]
            q = edge_id[u, v]
            tail, head = (u, v) if shift in shifts[:5] else (v, u)
            oriented[q] = tail, head
            if position[tail] > position[head]:
                deleted_by_shift[shift] += 1
            outgoing[tail] += 1
            incoming[head] += 1
    assert len(oriented) == 910
    assert incoming == outgoing == [5] * 182
    deleted = sorted(q for q, (u, v) in oriented.items()
                     if position[u] > position[v])
    assert len(deleted) == 166
    assert [deleted_by_shift[s] for s in shifts] == [3, 4, 15, 18, 28, 40, 31, 23, 3, 1]
    potential = {}
    for i in range(91):
        potential[ps[i]] = 2 * (i // 4)
        potential[ls[i]] = 2 * ((i + 1) // 4) - 1
    return edges, incident, rows, order, deleted, potential, deleted_by_shift

def code_rules(edges, incident, rows, order, deleted, constants=(), ignored=()):
    assert sorted(order) == list(range(182))
    deleted, constants, ignored = set(deleted), set(constants), set(ignored)
    assert constants <= deleted and ignored <= deleted
    assert not constants & ignored
    sessions = sorted(deleted - constants - ignored)
    position = {v: i for i, v in enumerate(order)}
    early = [min(edge, key=position.get) for edge in edges]
    late = [max(edge, key=position.get) for edge in edges]
    physical = [q for q in range(910) if q not in deleted]
    session_set = set(sessions)
    rules = []
    for v, qs in enumerate(incident):
        known = [q for q in qs if q in constants
                 or q not in deleted and late[q] == v
                 or q in session_set and early[q] == v]
        assert len(known) == 5
        matrix = matmul(rows[v], inverse([rows[v][qs.index(q)] for q in known]))
        assert [matrix[qs.index(q)] for q in known] == [
            [int(i == j) for j in range(5)] for i in range(5)]
        unknown = [matrix[i] for i, q in enumerate(qs) if q not in known]
        assert matmul(list(zip(*unknown)), unknown) == [
            [2 * int(i == j) for j in range(5)] for i in range(5)]
        rules.append((known, matrix))
    return dict(order=order, deleted=deleted, constants=constants,
                sessions=sessions, early=early, late=late,
                physical=physical, rules=rules)

def propagate(edges, incident, code, mutation=None):
    k = len(code["sessions"])
    inputs = {q: [int(i == j) for i in range(k)]
              for j, q in enumerate(code["sessions"])}
    values = {q: [0] * k for q in code["constants"]}
    decoded = {}
    local = {}
    for v in code["order"]:
        known, matrix = code["rules"][v]
        available = []
        for q in known:
            if q in inputs and code["early"][q] == v:
                available.append(inputs[q])
            else:
                assert q in values, "Noncausal dependency"
                assert q in code["constants"] or code["late"][q] == v
                available.append(values[q])
        local[v] = {}
        for q, coefficients in zip(incident[v], matrix):
            coefficients = list(coefficients)
            if mutation and (v, q) == mutation[:2]:
                coefficients[mutation[2]] = ADD[coefficients[mutation[2]]][1]
            form = [0] * k
            for a, source in zip(coefficients, available):
                if a:
                    form = [ADD[x][MUL[a][y]] for x, y in zip(form, source)]
            local[v][q] = form
            if q not in code["deleted"] and code["early"][q] == v:
                assert q not in values
                values[q] = form
            if q in inputs and code["late"][q] == v:
                decoded[q] = form
    errors = [q for q in inputs if decoded.get(q) != inputs[q]]
    consistent = all(local[u][q] == local[v][q]
                     for q, (u, v) in enumerate(edges))
    nonzero = all(any(values[q]) for q in code["physical"])
    return values, errors, consistent, nonzero

def check_paths(paths, pairs, physical_edges, rate, directed,
                common_denominator=None, common_numerator=None):
    edge_index = {tuple(sorted(edge)): i for i, edge in enumerate(physical_edges)}
    if common_denominator is None:
        loads = [Fraction(0)] * (2 * len(physical_edges) if directed
                                else len(physical_edges))
        totals = [Fraction(0)] * len(pairs)
        limit, demand = Fraction(1), rate
    else:
        assert type(common_denominator) is int and common_denominator > 0
        assert type(common_numerator) is int and common_numerator > 0
        assert Fraction(common_numerator, common_denominator) == rate
        loads = [0] * (2 * len(physical_edges) if directed else len(physical_edges))
        totals = [0] * len(pairs)
        limit, demand = common_denominator, common_numerator
    for row in paths:
        if common_denominator is None:
            j, numerator, denominator, vertices = row
            assert type(numerator) is type(denominator) is int
            assert numerator > 0 and denominator > 0
            amount = Fraction(numerator, denominator)
        else:
            j, amount, vertices = row
            assert type(amount) is int and amount > 0
        assert type(j) is int and 0 <= j < len(pairs)
        assert (vertices[0], vertices[-1]) == pairs[j], "Path endpoint mismatch"
        assert len(vertices) >= 2 and len(set(vertices)) == len(vertices)
        for u, v in zip(vertices, vertices[1:]):
            edge = tuple(sorted((u, v)))
            assert edge in edge_index, "Missing physical edge"
            q = edge_index[edge]
            index = 2 * q + int(u != physical_edges[q][0]) if directed else q
            loads[index] += amount
        totals[j] += amount
    assert totals == [demand] * len(pairs), "Flow demand mismatch"
    assert all(x <= limit for x in loads), "Capacity violation"
    return loads, limit

def shortest(physical_edges, weights, pairs):
    adjacency = [[] for _ in range(182)]
    for i, (u, v) in enumerate(physical_edges):
        adjacency[u].append((v, weights[2 * i]))
        adjacency[v].append((u, weights[2 * i + 1]))
    cache = {}
    result = []
    for s, t in pairs:
        if s not in cache:
            dist = [None] * 182
            dist[s] = 0
            queue = [(0, s)]
            while queue:
                du, u = heappop(queue)
                if du != dist[u]:
                    continue
                for v, weight in adjacency[u]:
                    value = du + weight
                    if dist[v] is None or value < dist[v]:
                        dist[v] = value
                        heappush(queue, (value, v))
            cache[s] = dist
        assert cache[s][t] is not None
        result.append(cache[s][t])
    return result

def coding_cut(code, edges, pairs, required=None):
    candidates = []
    for v in range(182):
        degree = sum(v in edges[q] for q in code["physical"])
        count = sum(s == v and t != v for s, t in pairs)
        if count and count == degree:
            candidates.append(v)
    assert candidates, "No rate-one source cut"
    if required is not None:
        assert required in candidates
        return required
    return candidates[0]

def main():
    if not __debug__:
        raise RuntimeError("Run without optimized mode so all exact checks remain active")
    start = time.monotonic()
    cert = json.loads(Path(__file__).with_name("certificates.json").read_text(
        encoding="utf-8"))
    edges, incident, rows, order, deleted, potential, shift_counts = reconstruct()
    assert cert["singer"]["order"] == order
    assert cert["singer"]["deleted"] == deleted
    forward = code_rules(edges, incident, rows, order, deleted)
    reverse = code_rules(edges, incident, rows, order[::-1], deleted)
    fv, fe, fc, fn = propagate(edges, incident, forward)
    rv, re, rc, rn = propagate(edges, incident, reverse)
    assert not fe and not re and fc and rc and fn and rn
    assert all(fv[q] == rv[q] for q in forward["physical"])
    assert all(set(reverse["rules"][v][0]) == set(incident[v])
               - set(forward["rules"][v][0]) for v in range(182))
    physical = [edges[q] for q in forward["physical"]]
    pairs = [(forward["early"][q], forward["late"][q])
             for q in forward["sessions"]]
    paired_pairs = pairs + [(t, s) for s, t in pairs]
    cut = coding_cut(forward, edges, pairs, required=10)
    assert order[0] == cut
    shared = cert["shared"]
    assert shared["metric"] == [1] * 744
    sl, sd = check_paths(shared["paths"], pairs, physical, Fraction(372, 415),
                         False, shared["denominator"], shared["session_numerator"])
    assert sl == [sd] * 744
    hops = shortest(physical, [1] * 1488, pairs)
    assert hops == [5] * 166 and Fraction(744, sum(hops)) == Fraction(372, 415)
    # Reconstruct the zero-return metric from the closed-form Singer potential.
    weights = []
    for u, v in physical:
        for tail, head in [(u, v), (v, u)]:
            numerator = 1 + potential[head] - potential[tail]
            assert numerator % 2 == 0
            weights.append(max(0, numerator // 2))
    assert cert["oneway"]["metric"] == weights
    position = {v: i for i, v in enumerate(order)}
    for q, (u, v) in enumerate(physical):
        forward_arc = 2 * q + int(position[u] > position[v])
        assert weights[forward_arc] > 0 and weights[forward_arc ^ 1] == 0
    check_paths(cert["oneway"]["paths"], pairs, physical, Fraction(3025, 3068), True)
    distances = shortest(physical, weights, pairs)
    assert sum(weights) == 3025 and sum(distances) == 3068
    assert Fraction(sum(weights), sum(distances)) == Fraction(3025, 3068)
    half_hop_cost = Fraction(len(physical), 2)
    half_hop_demand = Fraction(sum(hops), 2)
    physical_increment = Fraction(sum(
        potential[forward["late"][q]] - potential[forward["early"][q]]
        for q in forward["physical"]), 2)
    demand_increment = Fraction(sum(potential[t] - potential[s]
                                    for s, t in pairs), 2)
    assert (half_hop_cost, half_hop_demand) == (372, 415)
    assert physical_increment == demand_increment == 2653
    assert sum(weights) == half_hop_cost + physical_increment
    assert sum(distances) == half_hop_demand + demand_increment
    assert all(2 * distance == 5 + potential[t] - potential[s]
               for distance, (s, t) in zip(distances, pairs))
    assert sum(distances) - sum(weights) == half_hop_demand - half_hop_cost == 43
    # Explicit paired flow must be exactly shared flow plus its reversed copy.
    paired = cert["paired"]
    assert paired["denominator"] == shared["denominator"]
    assert paired["session_numerator"] == shared["session_numerator"]
    assert paired["metric"] == [1] * 1488
    expected = Counter()
    for j, count, path in shared["paths"]:
        expected[j, count, tuple(path)] += 1
        expected[j + 166, count, tuple(reversed(path))] += 1
    assert Counter((j, count, tuple(path)) for j, count, path
                   in paired["paths"]) == expected
    pl, pd = check_paths(paired["paths"], paired_pairs, physical, Fraction(372, 415),
                         True, paired["denominator"], paired["session_numerator"])
    assert pl == [pd] * 1488
    assert shortest(physical, [1] * 1488, paired_pairs) == [5] * 332
    coding_cut(forward, edges, paired_pairs, required=10)
    # Appendix-A published-instance contrast, reconstructed with its controls.
    old = cert["bh157"]
    for key in BH:
        assert old[key] == BH[key], "Published appendix data mismatch"
    bh = code_rules(edges, incident, rows, old["order"], old["deleted"],
                    old["constants"], old["ignored"])
    bv, be, bc, bn = propagate(edges, incident, bh)
    assert len(bh["sessions"]) == 157 and len(bh["physical"]) == 735
    assert not be and bc and bn
    bp = [(bh["early"][q], bh["late"][q]) for q in bh["sessions"]]
    bg = [edges[q] for q in bh["physical"]]
    check_paths(old["paths"], bp, bg, Fraction(1), True)
    bh_cut = coding_cut(bh, edges, bp, required=23)
    assert old["cut_vertex"] == bh_cut
    bh_metric = [int(tail == bh_cut) for u, v in bg
                 for tail, head in [(u, v), (v, u)]]
    assert old["metric"] == bh_metric
    assert sum(bh_metric) == sum(shortest(bg, bh_metric, bp)) == 5
    # Two deliberate corruptions: one local coefficient, one flow amount.
    _, bad_errors, _, _ = propagate(edges, incident, forward, (10, 101, 0))
    assert len(bad_errors) == 5
    bad_paths = copy.deepcopy(shared["paths"])
    bad_paths[0][1] += shared["denominator"]
    rejected = False
    try:
        check_paths(bad_paths, pairs, physical, Fraction(372, 415), False,
                    shared["denominator"], shared["session_numerator"])
    except AssertionError:
        rejected = True
    assert rejected, "Corrupted flow was accepted"
    result = {
        "status": "PASS",
        "reconstructed": {"vertices": 182, "candidate_edges": 910,
                          "singer_sessions": 166, "physical_edges": 744,
                          "constants": 0, "ignored_edges": 0,
                          "candidate_degree": 10,
                          "original_indegree": 5, "original_outdegree": 5,
                          "singer_shifts": list(shift_counts),
                          "deleted_counts_by_shift": list(shift_counts.values()),
                          "perfect_difference_set": {"modulus": 91, "size": 10,
                                                     "nonzero_difference_multiplicity": 1}},
        "local_matrix_checks": {"forward": 182, "reverse": 182, "bh157": 182,
                                "identity": "M transpose M = -I over F9"},
        "decoded_independent_inputs": {"forward": 166, "reverse": 166, "bh157": 157},
        "equal_forward_reverse_edge_forms": True,
        "potential_identity": {"half_hop_cost": int(half_hop_cost),
                               "half_hop_demand": int(half_hop_demand),
                               "physical_increment": int(physical_increment),
                               "demand_increment": int(demand_increment),
                               "preserved_margin": 43},
        "shared": {"routing": "372/415", "coding": "1", "gap": "415/372"},
        "oneway_full_duplex": {"routing": "3025/3068", "coding": "1",
                              "gap": "3068/3025", "metric_cost": 3025,
                              "distance_sum": 3068, "reverse_metric_weights": 0,
                              "unlimited_reverse_capacity_bound": True},
        "paired_full_duplex": {"routing": "372/415", "coding": "1",
                               "gap": "415/372", "independent_inputs": 332},
        "published_bh157_full_duplex": {"routing": "1", "coding": "1", "gap": "1",
                                       "source_cut_vertex": bh_cut},
        "negative_controls": {"coefficient_corruption_rejected": True,
                              "flow_corruption_rejected": True},
        "elapsed_seconds": round(time.monotonic() - start, 3)
    }
    print(json.dumps(result, indent=2))

if __name__ == "__main__":
    main()
