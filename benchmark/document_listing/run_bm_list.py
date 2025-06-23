import argparse
from pathlib import Path
import json
import subprocess


def main():
    # Parsing args
    parser = argparse.ArgumentParser(
        description="Run benchmarks for list operation",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    script_dir = Path(__file__).resolve().parent
    parser.add_argument("-b", "--benchmarks", default=script_dir / "benchmarks.json", help="Benchmarks config file")
    parser.add_argument("-c", "--cmd_path", default=script_dir / "build", help="Benchmarks command path")
    parser.add_argument("-o", "--output_path", default="./", help="Output path")
    parser.add_argument("-s", "--stats", action="store_true", help="Report statistics (mean, median, ...)")
    parser.add_argument("-d", "--details", action="store_true", help="Print detailed result")
    parser.add_argument("--patterns_filename", default="patterns", help="Name of files with patterns")
    parser.add_argument("--patterns_path", default="", help="Root path of files with patterns")
    parser.add_argument("--pattern_code", default="PLAIN", const="PLAIN", nargs="?", choices=["PLAIN", "BASE64"],
                        help="Patterns encoding")
    parser.add_argument("-g", "--group", action="store_true", help="Group of collections")
    parser.add_argument("collection_path", help="Collection path")
    parser.add_argument("index_path", help="Index path")
    args = parser.parse_args()

    # Reading benchmark configurations
    benchmarks = json.load(open(args.benchmarks))

    collection_path = Path(args.collection_path).resolve()
    index_path = Path(args.index_path).resolve()
    if args.group:
        # Processing collections group
        for child in collection_path.iterdir():
            if child.is_dir():
                process_collection(benchmarks, child, index_path / child.name, args)
    else:
        # Processing single collection
        process_collection(benchmarks, collection_path, index_path, args)


def esc(code):
    return f'\033[{code}m'


def process_collection(global_benchmarks, collection_path, index_path, args):
    collection_path = Path(collection_path).resolve()
    index_path = Path(index_path).resolve()
    output_path = Path(args.output_path).resolve()
    bm_cmd_path = Path(args.cmd_path).resolve()

    # Creating output directory for the given collection
    collection_name = collection_path.name
    print(f"{esc('1;42;34')}Collection '{collection_name}'{esc(0)}")
    output_path = output_path / collection_name
    output_path.mkdir(parents=True, exist_ok=True)

    patterns_path = Path(args.patterns_path).resolve() / collection_name if args.patterns_path else collection_path
    patterns_path /= args.patterns_filename

    local_benchmarks_path = collection_path / Path(args.benchmarks).resolve().name
    benchmarks = global_benchmarks if not local_benchmarks_path.exists() else json.load(open(local_benchmarks_path))

    # Running each benchmark
    for key, value in benchmarks.items():
        cmd_list = value.get("cmd_list")
        if cmd_list is None:
            continue

        print(f" {esc('34')}Index '{key}'{esc(0)}")

        cmd = str(bm_cmd_path / cmd_list)
        cmd += " --benchmark_counters_tabular=true" \
               " --benchmark_out_format=json" \
               " --benchmark_out=" + collection_name + "-" + key + "-list.json"
        cmd += " --data_dir=" + str(index_path / value.get("index_dir", ""))
        cmd += " --patterns=" + str(patterns_path)
        cmd += " --pattern_code=" + args.pattern_code
        if args.stats:
            cmd += " --report_stats"
        if args.details:
            cmd += " --print_result"
        cmd += " 2>" + key + "_list-error.txt"

        print(f"  {esc('38;5;22')}Running '{cmd}'{esc(0)}")

        subprocess.run(cmd, shell=True, cwd=output_path)


if __name__ == "__main__":
    main()
