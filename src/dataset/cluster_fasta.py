import subprocess
import argparse
import logging


def cd_hit(
    executable_location,
    input_fasta,
    output_fasta,
    log_file,
    identity_threshold,
    n_threads=4,
    memory=4096,
    verbose=True,
):
    log = logging.getLogger("CDHIT")
    logging.basicConfig(
        format="[%(levelname)s] %(filename)s: %(message)s",
        level=logging.DEBUG if verbose else logging.INFO,
        filename=log_file,
    )
    log.debug("#" * 60)
    log.debug(f"clustering file {input_fasta} with threshold {identity_threshold}...")
    try:
        word_length = {40: 2, 50: 3, 60: 4, 70: 5, 80: 5, 90: 5, 100: 5}[
            identity_threshold
        ]
    except KeyError:
        raise ValueError(f"Incorrect identity threshold: {identity_threshold}")
    execution = [
        executable_location,
        "-i",
        input_fasta,
        "-o",
        output_fasta,
        "-c",
        str(identity_threshold / 100.0),
        "-n",
        str(word_length),
        "-T",
        str(n_threads),
        "-M",
        str(memory),
    ]
    result = subprocess.run(
        execution, check=True, stdout=subprocess.PIPE, universal_newlines=True
    )
    for line in result.stdout.split("\n"):
        if "finished" in line:
            line = line.split()
            log.debug(
                f"clustered {line[0]} sequences into {line[2]} clusters at threshold {identity_threshold}"
            )
            break


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--executable",
        type=str,
        required=True,
    )
    parser.add_argument(
        "--input-fasta",
        type=str,
        required=True,
    )
    parser.add_argument(
        "--output-log",
        type=str,
    )
    parser.add_argument(
        "--output-fasta",
        type=str,
        required=True,
    )
    parser.add_argument(
        "--identity-threshold",
        type=int,
        choices=[40, 50, 60, 70, 80, 90, 100],
        required=True,
    )

    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--memory", type=int, default=4096)
    parser.add_argument(
        "--verbose",
        action="store_true",
    )
    args = parser.parse_args()
    cd_hit(
        executable_location=args.executable,
        input_fasta=args.input_fasta,
        output_fasta=args.output_fasta,
        log_file=args.output_log,
        identity_threshold=args.identity_threshold,
        n_threads=args.threads,
        memory=args.memory,
        verbose=args.verbose,
    )
