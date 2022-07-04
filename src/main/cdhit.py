import tempfile
import subprocess
import os

from main.fasta import read_fasta, write_fasta


def __auto_word_length(identity_threshold: int):
    if identity_threshold >= 70:
        return 5
    if identity_threshold >= 60:
        return 4
    if identity_threshold >= 50:
        return 3
    if identity_threshold >= 40:
        return 2
    raise ValueError("Invalid identity threshold: ", identity_threshold)


def __flatten_kwargs(**kwargs):
    kwargs_list = list()
    for k, v in kwargs.items():
        kwargs_list.append(k)
        kwargs_list.append(v)
    return kwargs_list


def cd_hit(
    fasta_data,
    identity_threshold,
    verbose=True,
    pandas_series=True,
    executable_location="cd-hit",
    n_threads=4,
    memory=4096,
    **kwargs,
):
    if pandas_series:
        fasta_data = list(
            zip(
                [">" + ac for ac in fasta_data.index.tolist()],
                fasta_data.values.tolist(),
            )
        )
    tmp_fasta_in = tempfile.NamedTemporaryFile(suffix=".fasta")
    write_fasta(fasta_file_name=tmp_fasta_in.name, fasta_data=fasta_data)

    tmp_fasta_out = tempfile.NamedTemporaryFile(suffix=".fasta")
    word_length = __auto_word_length(identity_threshold)
    execution = [
        executable_location,
        "-i",
        tmp_fasta_in.name,
        "-o",
        tmp_fasta_out.name,
        "-c",
        str(identity_threshold / 100.0),
        "-n",
        str(word_length),
        "-T",
        str(n_threads),
        "-M",
        str(memory),
    ] + __flatten_kwargs(**kwargs)
    result = subprocess.run(
        execution, check=True, stdout=subprocess.PIPE, universal_newlines=True
    )
    if verbose:
        for line in result.stdout.split("\n"):
            if "finished" in line:
                line = line.split()
                print(
                    f"cd-hit: clustered {line[0]} sequences into {line[2]} clusters at threshold {identity_threshold}"
                )
                break
    fasta_data_clustered = read_fasta(tmp_fasta_out.name)
    os.remove(tmp_fasta_out.name + ".clstr")
    tmp_fasta_in.close()
    tmp_fasta_out.close()
    return [ac[1:] for ac, _ in fasta_data_clustered]

