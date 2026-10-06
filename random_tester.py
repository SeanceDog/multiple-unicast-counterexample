#!/usr/bin/env python3
"""Keep testing independent random F9 inputs; stop at the first failure."""

import argparse
import json
import random
import secrets
import sys
import time
from pathlib import Path

from causal_code import CausalCode


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cases", type=int, default=0, help="0 = continue until Ctrl+C (default)")
    parser.add_argument("--report-every", type=int, default=1000)
    parser.add_argument("--seed", type=int, help="Optional reproducible random seed")
    parser.add_argument("--failure-file", type=Path, default=Path("failure.json"))
    args = parser.parse_args()
    if args.cases < 0 or args.report_every < 1:
        parser.error("--cases must be nonnegative and --report-every must be positive")

    seed = args.seed if args.seed is not None else secrets.randbits(64)
    rng = random.Random(seed)
    code = CausalCode()
    passed = 0
    start = time.monotonic()
    print(f"Testing 157 independent F9 symbols per case; seed={seed}.", flush=True)
    print("735 causal transmissions per case. Press Ctrl+C to stop.", flush=True)
    try:
        while args.cases == 0 or passed < args.cases:
            messages = [rng.randrange(9) for _ in range(157)]
            try:
                result = code.run(messages)
            except Exception as error:
                # Preserve the exact trial even if execution fails before decoding.
                args.failure_file.write_text(json.dumps({
                    "messages": messages,
                    "seed": seed,
                    "case": passed + 1,
                    "error": f"{type(error).__name__}: {error}",
                }, indent=2) + "\n")
                print(f"ERROR at case {passed + 1:,}; {passed:,} cases passed before error.", flush=True)
                print(f"{type(error).__name__}: {error}", flush=True)
                print(f"Replay: python3 causal_code.py {args.failure_file}", flush=True)
                return 1
            if not result.accepted:
                args.failure_file.write_text(json.dumps({
                    "messages": messages,
                    "decoded": result.decoded,
                    "mismatches": result.mismatches,
                    "seed": seed,
                    "case": passed + 1,
                }, indent=2) + "\n")
                print(f"FAIL at case {passed + 1:,}; {passed:,} cases passed before failure.", flush=True)
                print(f"Replay: python3 causal_code.py {args.failure_file}", flush=True)
                return 1
            passed += 1
            if passed % args.report_every == 0:
                elapsed = time.monotonic() - start
                print(f"{passed // args.report_every:,} x {args.report_every:,} cases passed "
                      f"({passed:,} total; {passed / max(elapsed, 1e-9):,.0f} cases/s)", flush=True)
    except KeyboardInterrupt:
        print(f"\nStopped: {passed:,} cases passed; no failed case.", flush=True)
        return 0
    print(f"PASS: {passed:,} cases passed; no failed case.", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
