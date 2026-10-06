#!/usr/bin/env python3
"""Independent standard-library tests for the actual-input GF9 simulator.

The oracle replays the previously audited transcript, using field arithmetic
implemented here without importing any arithmetic or interpolation helper from
the simulator.  Run: python3 -m unittest discover -s tests -v
"""
from __future__ import annotations

import collections
import copy
from dataclasses import replace
import json
import random
import unittest
from pathlib import Path

from causal_code import CausalCode

ROOT = Path(__file__).resolve().parents[1]
TRANSCRIPT = ROOT / 'data' / 'reference-transcript.json'


def plus(x: int, y: int) -> int:
    a, b = x % 3, x // 3
    c, d = y % 3, y // 3
    return (a + c) % 3 + 3 * ((b + d) % 3)


def times(x: int, y: int) -> int:
    a, b = x % 3, x // 3
    c, d = y % 3, y // 3
    # z^2 = -1 = 2 in F3.
    return (a * c - b * d) % 3 + 3 * ((a * d + b * c) % 3)


def oracle(messages, transcript):
    """Evaluate each recipe using only its sender's actual local observations."""
    pairs = transcript['pairs']
    inbox = [dict() for _ in range(182)]
    sources = [dict() for _ in range(182)]
    for i, value in enumerate(messages):
        sources[pairs[i][0]][i] = value
    packets = []

    def evaluate(recipe, vertex):
        value = 0
        for kind, index, coefficient in recipe:
            if kind == 'source':
                assert pairs[index][0] == vertex
                symbol = sources[vertex][index]
            else:
                assert kind == 'packet'
                # The dictionary lookup checks that this vertex has already
                # received this particular packet; global values are not read.
                symbol = inbox[vertex][index]
            value = plus(value, times(coefficient, symbol))
        return value

    for j, packet in enumerate(transcript['transmissions']):
        value = evaluate(packet['encoder'], packet['sender'])
        packets.append(value)
        inbox[packet['receiver']][j] = value
    decoded = [None] * len(messages)
    for decoder in transcript['decoders']:
        i = decoder['message']
        assert decoder['vertex'] == pairs[i][1]
        decoded[i] = evaluate(decoder['decoder'], decoder['vertex'])
    return tuple(decoded), tuple(packets)


class CausalCodeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.code = CausalCode()
        cls.transcript = json.loads(TRANSCRIPT.read_text())

    def assert_success(self, messages, *, replay=False):
        before = list(messages)
        result = self.code.run(messages, capture_trace=replay)
        self.assertIs(result.accepted, True)
        self.assertEqual(tuple(result.decoded), tuple(messages))
        self.assertEqual(result.transmissions, 735)
        self.assertEqual(tuple(result.mismatches), ())
        self.assertEqual(list(messages), before, 'run mutated caller input')
        if replay:
            decoded, packets = oracle(messages, self.transcript)
            self.assertEqual(tuple(result.decoded), decoded)
            self.assertEqual(tuple(packet.symbol for packet in result.trace), packets)
            self.assertEqual(
                tuple((packet.edge, packet.sender, packet.receiver) for packet in result.trace),
                tuple((packet['edge'], packet['sender'], packet['receiver'])
                      for packet in self.transcript['transmissions']))
        else:
            self.assertEqual(tuple(result.trace), ())
        return result

    def test_all_nonzero_single_source_assignments(self):
        # All eight scalar multiples, rather than only the prime-field basis.
        for i in range(157):
            for value in range(1, 9):
                with self.subTest(message=i, value=value):
                    messages = [0] * 157
                    messages[i] = value
                    self.assert_success(messages, replay=(value == 1))

    def test_constants_zero_and_random_against_local_replay(self):
        for value in range(9):
            self.assert_success([value] * 157, replay=True)
        rng = random.Random(0x182735)
        for trial in range(64):
            with self.subTest(trial=trial):
                self.assert_success([rng.randrange(9) for _ in range(157)], replay=True)

    def test_sources_are_owned_at_correct_vertices(self):
        rng = random.Random(157)
        messages = [rng.randrange(9) for _ in range(157)]
        sources = collections.defaultdict(dict)
        for i, (source, _) in enumerate(self.code.pairs):
            sources[source][i] = messages[i]
        saved = {v: dict(values) for v, values in sources.items()}
        result = self.code.run_at_sources(sources)
        self.assertIs(result.accepted, True)
        self.assertEqual(tuple(result.decoded), tuple(messages))
        self.assertEqual(result.transmissions, 735)
        self.assertEqual(dict(sources), saved)

    def test_repeated_runs_do_not_reuse_previous_state(self):
        first = [8] * 157
        self.assert_success(first)
        zero = self.assert_success([0] * 157, replay=True)
        self.assertEqual(tuple(packet.symbol for packet in zero.trace), (0,) * 735)
        self.assert_success(first, replay=True)

    def test_malformed_messages_are_rejected(self):
        bad_values = (-1, 9, 10, True, False, 1.0, '1', None)
        for bad in bad_values:
            messages = [0] * 157
            messages[73] = bad
            with self.subTest(value=repr(bad)):
                with self.assertRaises((TypeError, ValueError)):
                    self.code.run(messages)
        for messages in ([], [0] * 156, [0] * 158, '0' * 157, None):
            with self.subTest(input=repr(messages)[:30]):
                with self.assertRaises((TypeError, ValueError)):
                    self.code.run(messages)

    def test_source_mapping_cannot_supply_foreign_or_missing_inputs(self):
        sources = collections.defaultdict(dict)
        for i, (source, _) in enumerate(self.code.pairs):
            sources[source][i] = 0
        source, sink = self.code.pairs[0]
        missing = {v: dict(values) for v, values in sources.items()}
        del missing[source][0]
        with self.assertRaises((TypeError, ValueError)):
            self.code.run_at_sources(missing)
        foreign = {v: dict(values) for v, values in missing.items()}
        foreign.setdefault(sink, {})[0] = 0
        with self.assertRaises((TypeError, ValueError)):
            self.code.run_at_sources(foreign)
        extra = {v: dict(values) for v, values in sources.items()}
        extra[source][157] = 0
        with self.assertRaises((TypeError, ValueError)):
            self.code.run_at_sources(extra)
        boolean = {v: dict(values) for v, values in sources.items()}
        boolean[source][0] = True
        with self.assertRaises((TypeError, ValueError)):
            self.code.run_at_sources(boolean)

    def test_runtime_schedule_uses_only_local_available_ports(self):
        received = [set() for _ in range(182)]
        sent = set()
        seen_vertices = set()
        decoded = set()
        for step in self.code.steps:
            vertex = step.vertex
            self.assertNotIn(vertex, seen_vertices)
            seen_vertices.add(vertex)
            self.assertEqual(len(step.known), 5)
            self.assertEqual(len(step.inverse), 5)
            self.assertTrue(all(len(row) == 5 for row in step.inverse))
            for port in step.known:
                if port.kind == 'source':
                    self.assertEqual(self.code.pairs[port.index][0], vertex)
                elif port.kind == 'packet':
                    self.assertIn(port.index, received[vertex])
                else:
                    self.assertEqual(port.kind, 'zero')
            for send in step.sends:
                self.assertNotIn(send.edge, sent)
                self.assertIn(send.edge, self.code.physical_edges)
                self.assertEqual(set(self.code.geometry_edges[send.edge]),
                                 {vertex, send.receiver})
                sent.add(send.edge)
                received[send.receiver].add(send.edge)
            for decoder in step.decodes:
                self.assertNotIn(decoder.message, decoded)
                self.assertEqual(self.code.pairs[decoder.message][1], vertex)
                decoded.add(decoder.message)
        self.assertEqual(seen_vertices, set(range(182)))
        self.assertEqual(sent, set(self.code.physical_edges))
        self.assertEqual(decoded, set(range(157)))

    def test_corrupted_local_encoder_is_rejected(self):
        # Alter a private copy of the schedule; the production code is untouched.
        mutant = copy.copy(self.code)
        first = mutant.steps[0]
        source = next(port.index for port in first.known if port.kind == 'source')
        broken = tuple(replace(send, feature=(0,) * 5) for send in first.sends)
        mutant.steps = (replace(first, sends=broken),) + mutant.steps[1:]
        messages = [0] * 157
        messages[source] = 1
        result = mutant.run(messages, capture_trace=True)
        self.assertIs(result.accepted, False)
        self.assertIn(source, result.mismatches)
        self.assertNotEqual(tuple(result.decoded), tuple(messages))
        self.assertEqual(result.transmissions, 735)
        # A later healthy run must not inherit this private mutation.
        self.assert_success(messages, replay=True)

    def test_future_packet_reference_cannot_be_read(self):
        mutant = copy.copy(self.code)
        first = mutant.steps[0]
        unavailable = replace(first.known[0], kind='packet', index=first.sends[0].edge)
        mutant.steps = (replace(first, known=(unavailable,) + first.known[1:]),) + mutant.steps[1:]
        with self.assertRaisesRegex(RuntimeError, 'not received'):
            mutant.run([0] * 157)

    def test_graph_and_routing_cost_independently(self):
        transcript = self.transcript
        self.assertEqual(tuple(map(tuple, transcript['pairs'])), tuple(self.code.pairs))
        self.assertEqual(len(transcript['transmissions']), 735)
        edges = [transcript['geometry_edges'][e] for e in transcript['physical_edge_indices']]
        self.assertEqual(len(edges), 735)
        self.assertEqual(len(set(map(tuple, edges))), 735)
        adjacency = [set() for _ in range(182)]
        for u, v in edges:
            self.assertNotEqual(u, v)
            adjacency[u].add(v)
            adjacency[v].add(u)
        distances = []
        for source, target in self.code.pairs:
            self.assertNotEqual(source, target)
            distance = {source: 0}
            queue = collections.deque([source])
            while queue:
                u = queue.popleft()
                for v in adjacency[u]:
                    if v not in distance:
                        distance[v] = distance[u] + 1
                        queue.append(v)
            distances.append(distance[target])
        self.assertEqual(distances, [5] * 157)
        self.assertEqual(sum(distances), 785)
        self.assertLess(735, sum(distances))

    def test_independent_field_reference(self):
        for x in range(9):
            self.assertEqual(plus(x, 0), x)
            self.assertEqual(times(x, 1), x)
            self.assertEqual(times(x, 0), 0)
            self.assertEqual(plus(plus(x, x), x), 0)
            if x:
                self.assertTrue(any(times(x, y) == 1 for y in range(1, 9)))
            for y in range(9):
                for z in range(9):
                    self.assertEqual(times(x, plus(y, z)), plus(times(x, y), times(x, z)))
        self.assertEqual(times(3, 3), 2)


if __name__ == '__main__':
    unittest.main(verbosity=2)
