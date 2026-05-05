#!/usr/bin/env python3
"""
order-domains.py

Robust parser + overlap resolver that picks non-overlapping hits in
order of best confidence (i-evalue, then longer length, then higher score).

Usage:
  # Read from file:
  python3 order-domains.py domtbl_tabbed.txt [--debug]

  # Read from stdin (piped from grep or other command):
  grep pattern file.txt | python3 order-domains.py [--debug]
  cat domtbl_tabbed.txt | python3 order-domains.py [--debug]
"""

import sys, re, argparse
from collections import defaultdict

sci_re = re.compile(r"^[0-9]*\.?[0-9]+[eE][-+]?\d+$")
float_re = re.compile(r"^[0-9]*\.?[0-9]+$")


def as_int(x):
    try:
        return int(x)
    except:
        return None


def as_float(x):
    try:
        return float(x)
    except:
        return None


def parse_line(tokens):
    # tokens expected from tab-split; fallback to whitespace split if tiny
    domain = tokens[0] if len(tokens) > 0 else None
    query = tokens[3] if len(tokens) > 3 else None

    # try env coords at canonical positions if present
    env_from = env_to = None
    if len(tokens) > 20:
        env_from = as_int(tokens[19])
        env_to = as_int(tokens[20])

    # fallback to ali -> hmm if env missing
    if env_from is None or env_to is None:
        if len(tokens) > 18:
            env_from = env_from or as_int(tokens[17])
            env_to = env_to or as_int(tokens[18])
    if env_from is None or env_to is None:
        if len(tokens) > 16:
            env_from = env_from or as_int(tokens[15])
            env_to = env_to or as_int(tokens[16])

    # detect acc token (0..1) by scanning from right and get env_from/env_to as two tokens before it if needed
    acc = None
    if env_from is None or env_to is None or len(tokens) <= 20:
        for i in range(len(tokens) - 1, 1, -1):
            v = as_float(tokens[i])
            if v is not None and 0.0 <= v <= 1.0:
                acc = v
                # tokens i-2, i-1 expected to be env_from/env_to
                a = as_int(tokens[i - 2])
                b = as_int(tokens[i - 1])
                if a is not None and b is not None and a <= b:
                    env_from, env_to = a, b
                    break

    # i-evalue: prefer standard column if present else smallest sci token
    i_eval = None
    if len(tokens) > 12:
        i_eval = as_float(tokens[12])
    if i_eval is None:
        evs = [as_float(t) for t in tokens if t and sci_re.match(t)]
        evs = [e for e in evs if e is not None]
        i_eval = min(evs) if evs else float("inf")

    # score: prefer standard col else best plausible float
    score = None
    if len(tokens) > 13:
        score = as_float(tokens[13])
    if score is None:
        scores = [
            as_float(t)
            for t in tokens
            if t and float_re.match(t) and not sci_re.match(t)
        ]
        scores = [s for s in scores if s is not None and 0 <= s <= 5000]
        score = max(scores) if scores else 0.0

    length = (
        (env_to - env_from + 1) if (env_from is not None and env_to is not None) else 0
    )

    return {
        "domain": domain,
        "query": query,
        "env_from": env_from,
        "env_to": env_to,
        "i_eval": i_eval,
        "score": score,
        "acc": acc,
        "length": length,
        "tokens": tokens,
    }


def parse_file(path=None):
    """
    Parse domain data from file or stdin.
    If path is None, reads from stdin (useful for piped input from grep).
    If path is provided, reads from that file.
    """
    by_query = defaultdict(list)

    if path is None:
        # Read from stdin
        fh = sys.stdin
    else:
        # Read from file
        fh = open(path)

    try:
        for ln in fh:
            line = ln.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            tokens = re.split(r"\t+", line)
            if len(tokens) < 6:
                tokens = line.split()
            parsed = parse_line(tokens)
            if parsed["query"] is None:
                continue
            by_query[parsed["query"]].append(parsed)
    finally:
        # Only close if we opened a file (not stdin)
        if path is not None:
            fh.close()

    return by_query


def select_nonoverlapping_bestfirst(hits):
    """
    Best-first selection:
    - Sort hits by preference tuple: (i_eval asc, -length desc, -score desc)
    - Iterate, adding a hit if it does NOT overlap any already chosen hit.
    """
    # discard hits without valid coords
    hits = [h for h in hits if h["env_from"] is not None and h["env_to"] is not None]
    # sort by preference
    hits.sort(key=lambda h: (h["i_eval"], -h["length"], -h["score"]))
    chosen = []
    for h in hits:
        overlaps = False
        for c in chosen:
            if not (h["env_to"] < c["env_from"] or h["env_from"] > c["env_to"]):
                overlaps = True
                break
        if not overlaps:
            chosen.append(h)
    # final sort by coordinate to get N->C order
    chosen.sort(key=lambda x: x["env_from"])
    return chosen


def main():
    p = argparse.ArgumentParser(
        description="Parse and select non-overlapping domain hits with best-first strategy",
        epilog="If DOMTBL is omitted, reads from stdin (useful for piped input: grep ... | order-domain.py)",
    )
    p.add_argument(
        "domtbl",
        nargs="?",
        default=None,
        help="tabbed domtbl file (optional; reads from stdin if omitted)",
    )
    p.add_argument("--debug", action="store_true", help="show debug output")
    args = p.parse_args()

    # Detect if input is piped or file is provided
    input_file = args.domtbl
    if input_file is None and sys.stdin.isatty():
        # No file provided and no stdin piped
        p.print_help()
        sys.exit(1)

    by_query = parse_file(input_file)
    for q in sorted(by_query):
        hits = by_query[q]
        chosen = select_nonoverlapping_bestfirst(hits)
        arch = "-".join(h["domain"] for h in chosen)
        if args.debug:
            print(f"=== QUERY: {q} ===")
            print("All parsed hits (unsorted):")
            for h in hits:
                print(
                    f"  {h['domain']:20s} {h['env_from']:4} - {h['env_to']:4} iE={h['i_eval']:.2e} score={h['score']:.1f} len={h['length']}"
                )
            print("Chosen (best-first, non-overlap):")
            for h in chosen:
                print(
                    f"  -> {h['domain']:20s} {h['env_from']:4} - {h['env_to']:4} iE={h['i_eval']:.2e} score={h['score']:.1f} len={h['length']}"
                )
            print()
        print(f"{q}\t{arch}")


if __name__ == "__main__":
    main()
