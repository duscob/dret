import argparse
import json
from pathlib import Path
import subprocess


def main():
    # Parsing args
    parser = argparse.ArgumentParser(
        description="Run benchmarks for indexes construction",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    script_dir = Path(__file__).resolve().parent
    parser.add_argument("-b", "--benchmarks", default=script_dir / "benchmarks_build.json",
                        help="Benchmarks config file")
    parser.add_argument("-c", "--cmd_path", default="./build", help="Benchmarks command path")
    parser.add_argument("-o", "--output_path", default="./", help="Output path")
    parser.add_argument("-g", "--group", action="store_true", help="Group of collections")
    parser.add_argument("-e", "--end_of_text", default="\\x01", help="End of text char")
    parser.add_argument("collection", help="Collection path")
    args = parser.parse_args()

    # Reading benchmark configurations
    benchmarks = json.load(open(args.benchmarks))

    collection_path = args.collection
    if args.group:
        # Processing a collections group
        for child in collection_path.iterdir():
            if child.is_dir():
                process_collection(benchmarks, child, args.output_path, args.cmd_path, args.end_of_text)
    else:
        # Processing a single collection
        process_collection(benchmarks, collection_path, args.output_path, args.cmd_path, args.end_of_text)


def esc(code):
    return f'\033[{code}m'


def process_collection(benchmarks, collection_path, output_path, bm_cmd_path, end_of_text):
    collection_path = Path(collection_path).resolve()
    output_path = Path(output_path).resolve()
    bm_cmd_path = Path(bm_cmd_path).resolve()

    # Creating output directory for the given collection
    collection_name = collection_path.name
    print(f"{esc('1;42;34')}Collection '{collection_name}'{esc(0)}")
    output_path = output_path / collection_name

    # Running each benchmark
    for key, value in benchmarks.items():
        print(f"\n {esc('34')}Index '{key}'{esc(0)}")

        # Creating output directory for the given collection
        index_output_path = output_path / value.get("subdir", "")
        index_output_path.mkdir(parents=True, exist_ok=True)

        if not create_output_data(collection_path, index_output_path, end_of_text):
            return

        cmd = Path(value["cmd"])  # .resolve()

        if not cmd.is_absolute():
            cmd = bm_cmd_path / cmd
        cmd = str(cmd)

        if value.get("benchmark_args", True):
            cmd += " --benchmark_counters_tabular=true" \
                   " --benchmark_dry_run" \
                   " --benchmark_repetitions=1" \
                   " --benchmark_out_format=json" \
                   " --benchmark_out=" + collection_name + "-" + key + "-build.json"

        # cmd += " --data=./data"
        cmd += " " + value["args"]

        cmd += " 2>" + key + "_build-error.txt"

        must_run = True
        if value.get("output_files", []):
            must_run = False
            for file in value["output_files"]:
                if Path(file).exists():
                    must_run = True
                    break

        if not must_run:
            print(f"  {esc('93')}WARNING: Skipping '{cmd}'{esc(0)}")
            return

        print(f"  {esc('38;5;22')}Running '{cmd}'{esc(0)}")

        subprocess.run(cmd, shell=True, cwd=index_output_path)


def create_output_data(input_path, output_path, end_of_text):
    output_data_path = output_path / "data"
    if output_data_path.exists():
        return True

    input_data_path = input_path / "data"

    cmd = f"grep -c -Pa '\\x00' {input_data_path}"
    n_null_chars = int(subprocess.run(cmd, shell=True, cwd=output_path, capture_output=True).stdout.decode())

    cmd = f"grep -c -Pa '{end_of_text}' {input_data_path}"
    n_etx_chars = int(subprocess.run(cmd, shell=True, cwd=output_path, capture_output=True).stdout.decode())

    if n_null_chars < 1 and n_etx_chars < 1:
        print(f"  {esc('31')}ERROR: Collection does not contain document delimiter.{esc(0)}")
        return False

    if 1 < n_null_chars and 1 <= n_etx_chars:
        print(f"  {esc('31')}ERROR: Collection contains mixed document delimiters.{esc(0)}")
        return False

    if n_null_chars <= 1 <= n_etx_chars:
        # Creating symbolic link to data file
        print(f"  # docs: {n_etx_chars}")
        print(f"  {esc('32')}Creating symbolic link to data file.{esc(0)}")
        output_data_path.symlink_to(input_path / "data")
    # elif n_etx_chars < 1 <= n_null_chars:
    else:
        # Creating new data replacing end of document delimiter
        print(f"  # docs: {n_null_chars}")
        print(f"  {esc('32')}Creating new data file replacing end of document delimiters.{esc(0)}")
        cmd = f"sed 's/\\x00/{end_of_text}/g' {input_data_path} > {output_data_path}"
        subprocess.run(cmd, shell=True, cwd=output_path)

    return True


if __name__ == "__main__":
    main()
