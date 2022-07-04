import gzip


def read_soft(input_file: str, delimit: str = "\t") -> list:
    # Doc: https://www.ncbi.nlm.nih.gov/geo/info/soft.html
    if input_file.endswith(".gz"):
        with gzip.open(input_file, "rb") as file:
            lines = [line.decode("utf-8") for line in file.readlines()]
    elif input_file.endswith(".soft"):
        with open(input_file) as file:
            lines = file.readlines()
    else:
        raise ValueError(f"unknown file type: {input_file}")
    lines = [
        line.rstrip("\n").split(delimit)
        for line in lines
        if len(line.strip()) > 0 and not line[0] in {"^", "#", "!"}
    ]
    return lines
