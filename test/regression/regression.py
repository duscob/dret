#!/usr/bin/env python3
"""Regression harness for the dret redesign: answers and sizes of every index, against a baseline.

Builds every cell of sweep.json from scratch on a small collection, queries it, and reduces the
run to a manifest: per benchmark, the md5 of its answers (bm_doc_list -print_results_dir, one
sorted line of document ids per pattern), its resident and serialized sizes, and its query time.

  make-data  cut a small collection out of a large one (whole documents, in order)
  run        build, query, write the manifest; --record GOLDEN stores it as the baseline,
             --check GOLDEN compares with it and exits non-zero on any difference

Checks, in order:
  1. every index except the count-only one gives the same answers (a cross-index check, no baseline needed);
  2. against the baseline: the same benchmarks, the same answers, the same sizes (exact; pass
     --allow-size-change when a change is meant to move sizes);
  3. query times are only reported (median and worst ratio against the baseline); pass
     --time-tolerance to gate on them, on a quiet machine.

A refactor that renames benchmarks passes --names NAMES.json, a {"old name": "new name"} map
applied to the baseline. Python standard library only.
"""

import argparse
import hashlib
import json
import os
import shutil
import statistics
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
COUNT_ONLY_MARK = "Count"  # the count-only brute kind reports no documents

# Sizes that differ between two builds of the same code, found by comparing two fresh builds.
# (substring of the benchmark name, size fields not compared, why). Answers are still compared.
# Empty: two builds agree on every size (README.md, "Sizes that are not reproducible").
UNSTABLE_SIZES = []


def unstable_fields(name):
    return {f for mark, fields, _ in UNSTABLE_SIZES if mark in name for f in fields}


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def md5(path):
    return hashlib.md5(Path(path).read_bytes()).hexdigest()


def make_data(args):
    docs = Path(args.source).read_bytes().split(b"\0")
    if docs and docs[-1] == b"":
        docs.pop()
    out, total = [], 0
    for d in docs:
        if args.max_doc_bytes and len(d) > args.max_doc_bytes:
            continue
        out.append(d)
        total += len(d) + 1
        if total >= args.max_bytes:
            break
    dest = Path(args.out)
    dest.mkdir(parents=True, exist_ok=True)
    (dest / "data").write_bytes(b"\0".join(out) + b"\0")
    print(f"{dest / 'data'}: {len(out)} documents, {total} bytes, sha256 {sha256(dest / 'data')}")


def write_spec(path, mode, dataset, sweep, reps, version):
    spec = {
        "spec_version": version,
        "mode": mode,
        "dataset": dataset,
        "workload": {"reps": reps, "min_time": 0.0},
        "sweep": sweep,
    }
    Path(path).write_text(json.dumps(spec, indent=2))


def run_bm(binary, spec, cwd, extra, log):
    cmd = [str(binary), "-spec", str(spec)] + extra
    with open(log, "w") as f:
        rc = subprocess.run(cmd, cwd=cwd, stdout=f, stderr=subprocess.STDOUT).returncode
    if rc != 0:
        sys.exit(f"bm_doc_list failed ({rc}): {' '.join(cmd)}\n  see {log}")


def build_manifest(args):
    binary = Path(args.bin).resolve()
    data = Path(args.data).resolve()
    patterns = Path(args.patterns).resolve()
    sweep_spec = json.loads(Path(args.sweep).read_text())
    sweep, version = sweep_spec["sweep"], sweep_spec.get("spec_version", 1)

    work = Path(args.work).resolve() if args.work else Path(tempfile.mkdtemp(prefix="dret-regression-"))
    if work.exists() and any(work.iterdir()):
        if not args.reuse:
            sys.exit(f"{work} is not empty; dret reuses cached files, so start from an empty directory "
                     "(or pass --reuse to query an existing build)")
    index_dir, results_dir = work / "index", work / "results"
    index_dir.mkdir(parents=True, exist_ok=True)
    if results_dir.exists():
        shutil.rmtree(results_dir)
    results_dir.mkdir()

    dataset = {"dir": str(data.parent), "name": data.name, "data_width": 8}
    if not args.reuse:
        write_spec(work / "construct.json", "construct", dataset, sweep, 1, version)
        print(f"construct -> {index_dir} (log {work / 'construct.log'})", flush=True)
        run_bm(binary, work / "construct.json", index_dir, [], work / "construct.log")

    dataset = {"dir": str(index_dir), "name": data.name, "data_width": 8, "patterns": str(patterns)}
    write_spec(work / "query.json", "query", dataset, sweep, args.reps, version)
    print(f"query -> {work / 'query.out.json'} (log {work / 'query.log'})", flush=True)
    run_bm(binary, work / "query.json", index_dir,
           ["-print_results_dir", str(results_dir),
            f"--benchmark_out={work / 'query.out.json'}", "--benchmark_out_format=json"],
           work / "query.log")

    out = json.loads((work / "query.out.json").read_text())
    entries = {}
    for b in out["benchmarks"]:
        if b.get("run_type") == "aggregate":
            continue
        name = b.get("run_name", b["name"])
        e = entries.setdefault(name, {"times": []})
        e["times"].append(b["Time_x_Pattern"])
        e["size_bytes"] = int(b["Size(bytes)"])
        e["serialized_bytes"] = int(b["serialized_bytes"])

    sanitized = {p.stem: p for p in results_dir.glob("*.txt")}
    manifest = {}
    for name, e in sorted(entries.items()):
        key = "".join(c if c.isalnum() or c in "-." else "_" for c in name)
        if key not in sanitized:
            sys.exit(f"no answers written for {name}")
        manifest[name] = {
            "answers_md5": md5(sanitized[key]),
            "size_bytes": e["size_bytes"],
            "serialized_bytes": e["serialized_bytes"],
            "time_x_pattern": statistics.median(e["times"]),
        }

    with open(work / "query.log") as f:
        validate = [line.rstrip() for line in f if line.startswith("[VALIDATE]")]

    return {
        "data": {"file": data.name, "sha256": sha256(data)},
        "patterns": {"file": patterns.name, "sha256": sha256(patterns)},
        "sweep_sha256": sha256(args.sweep),
        "commit": subprocess.run(["git", "-C", str(HERE), "rev-parse", "--short", "HEAD"],
                                 capture_output=True, text=True).stdout.strip(),
        "validate": validate,
        "indexes": manifest,
    }, work


def check_cross_index(m):
    md5s = {}
    for name, e in m["indexes"].items():
        if COUNT_ONLY_MARK in name:
            continue
        md5s.setdefault(e["answers_md5"], []).append(name)
    problems = list(m["validate"])
    if len(md5s) > 1:
        groups = sorted(md5s.values(), key=len, reverse=True)
        problems.append(f"indexes disagree: {len(groups[0])} agree, the others give different answers:")
        problems += [f"  {n}" for g in groups[1:] for n in g]
    return problems


def compare(golden, got, args):
    problems = []
    for key in ("data", "patterns"):
        if golden[key]["sha256"] != got[key]["sha256"]:
            problems.append(f"{key} differs from the baseline's ({golden[key]['file']}): not the same input")
    if golden["sweep_sha256"] != got["sweep_sha256"]:
        print("note: sweep.json changed since the baseline; only the common benchmarks are compared")

    names = json.loads(Path(args.names).read_text()) if args.names else {}
    base = {names.get(n, n): e for n, e in golden["indexes"].items()}
    now = got["indexes"]

    missing = sorted(set(base) - set(now))
    added = sorted(set(now) - set(base))
    problems += [f"missing: {n}" for n in missing]
    if added:
        print(f"note: {len(added)} benchmark(s) not in the baseline: " + ", ".join(added))

    ratios = []
    for n in sorted(set(base) & set(now)):
        b, g = base[n], now[n]
        if b["answers_md5"] != g["answers_md5"]:
            problems.append(f"answers changed: {n}")
        if not args.allow_size_change:
            for k in ("size_bytes", "serialized_bytes"):
                if k in unstable_fields(n):
                    continue
                if b[k] != g[k]:
                    problems.append(f"{k} changed: {n}: {b[k]} -> {g[k]}")
        if b["time_x_pattern"] > 0:
            ratios.append((g["time_x_pattern"] / b["time_x_pattern"], n))

    if ratios:
        ratios.sort()
        med = statistics.median(r for r, _ in ratios)
        print(f"query time vs baseline: median ratio {med:.3f}, "
              f"fastest {ratios[0][0]:.3f} ({ratios[0][1]}), slowest {ratios[-1][0]:.3f} ({ratios[-1][1]})")
        if args.time_tolerance is not None:
            problems += [f"slower by more than {args.time_tolerance:.0%}: {n} ({r:.3f})"
                         for r, n in ratios if r > 1 + args.time_tolerance]
    return problems


def run(args):
    manifest, work = build_manifest(args)
    (work / "manifest.json").write_text(json.dumps(manifest, indent=1, sort_keys=True))
    print(f"{len(manifest['indexes'])} benchmarks, manifest {work / 'manifest.json'}")

    problems = check_cross_index(manifest)
    if args.record:
        if problems:
            print("\n".join(problems))
            sys.exit("not recording a baseline from a run whose indexes disagree")
        Path(args.record).write_text(json.dumps(manifest, indent=1, sort_keys=True) + "\n")
        print(f"recorded baseline {args.record}")
    elif args.check:
        problems += compare(json.loads(Path(args.check).read_text()), manifest, args)

    if problems:
        print("FAIL\n" + "\n".join(problems))
        sys.exit(1)
    print("OK")


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)

    d = sub.add_parser("make-data", help="cut a small collection out of a large one")
    d.add_argument("--source", required=True, help="collection file, documents ended by \\0")
    d.add_argument("--out", required=True, help="directory to write <out>/data into")
    d.add_argument("--max-bytes", type=int, required=True, help="stop once this many bytes are taken")
    d.add_argument("--max-doc-bytes", type=int, default=0, help="skip documents longer than this")
    d.set_defaults(func=make_data)

    r = sub.add_parser("run", help="build, query and compare")
    r.add_argument("--bin", required=True, help="bm_doc_list binary (a Release build)")
    r.add_argument("--data", required=True, help="collection file")
    r.add_argument("--patterns", required=True, help="patterns file, one per line")
    r.add_argument("--sweep", default=str(HERE / "sweep.json"))
    r.add_argument("--work", help="empty working directory (default: a new temporary one)")
    r.add_argument("--reuse", action="store_true", help="query the build already in --work")
    r.add_argument("--reps", type=int, default=3, help="query repetitions (time = median)")
    mode = r.add_mutually_exclusive_group()
    mode.add_argument("--record", metavar="GOLDEN", help="store this run as the baseline")
    mode.add_argument("--check", metavar="GOLDEN", help="compare this run with the baseline")
    r.add_argument("--names", help="JSON map of renamed benchmarks, old name -> new name")
    r.add_argument("--allow-size-change", action="store_true")
    r.add_argument("--time-tolerance", type=float, help="fail if a query is slower by more than this fraction")
    r.set_defaults(func=run)

    args = p.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
