def read_fasta(fasta_file_name):
    fasta_data = []
    current_sequence = ""
    last_header = ""

    with open(fasta_file_name) as fasta_file:
        for line in fasta_file:
            line = line.strip()
            if line.startswith(">"):
                if current_sequence == "":
                    last_header = line
                    continue
                else:
                    fasta_data.append((last_header, current_sequence))
                    current_sequence = ""
                    # last_header = ""
                    last_header = line
            else:
                current_sequence += line
        if current_sequence != "":
            fasta_data.append((last_header, current_sequence))
            # current_sequence = ""
            # last_header = ""
    return fasta_data


def write_fasta(fasta_file_name, fasta_data):
    with open(fasta_file_name, "w") as outfile:
        for fasta_entry in fasta_data:
            outfile.write(fasta_entry[0] + "\n")
            sequence = fasta_entry[1]
            for i in range(0, len(sequence), 60):
                outfile.write(sequence[i : min(i + 60, len(sequence))] + "\n")
