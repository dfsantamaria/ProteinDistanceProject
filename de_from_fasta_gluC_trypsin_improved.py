import csv
import os


class ProteinAnalyzer:
    def __init__(self, fasta_file, output_file):
        self.fasta_file = fasta_file
        self.output_file = output_file
        self.sequences = {}

    # ---------- CLUSTERING ----------
    def clustering_distance(self, sequence, residues):
        residues = set(residues)

        positions = [i + 1 for i, aa in enumerate(sequence) if aa in residues]

        if len(positions) < 2:
            return None

        def adjusted_distance(a, b):
            return max(abs(a - b) - 1, 0)

        distances = []

        for i, pos in enumerate(positions):
            candidates = []

            if i > 0:
                candidates.append(
                    adjusted_distance(pos, positions[i - 1])
                )

            if i < len(positions) - 1:
                candidates.append(
                    adjusted_distance(pos, positions[i + 1])
                )

            distances.append(min(candidates))

        return sum(distances) / len(distances)

    # ---------- PEPTIDE LENGTHS ----------
    def peptide_lengths(self, sequence, residues):
        residues = set(residues)

        cleavage_sites = [
            i for i, aa in enumerate(sequence)
            if aa in residues
        ]

        lengths = []
        prev_cut = 0

        for site in cleavage_sites:
            peptide_length = site + 1 - prev_cut

            if peptide_length > 0:
                lengths.append(peptide_length)

            prev_cut = site + 1

        last_length = len(sequence) - prev_cut

        if last_length > 0:
            lengths.append(last_length)

        if len(lengths) == 0:
            return None, None, None, None

        avg_len = sum(lengths) / len(lengths)
        min_len = min(lengths)
        max_len = max(lengths)

        ms_peptides = [
            l for l in lengths
            if 7 <= l <= 25
        ]

        ms_fraction = len(ms_peptides) / len(lengths)

        return avg_len, min_len, max_len, ms_fraction

    # ---------- FASTA READER ----------
    def read_fasta(self):

        if not os.path.exists(self.fasta_file):
            raise FileNotFoundError(
                f"FASTA file not found: {self.fasta_file}"
            )

        sequences = {}
        name = None
        seq = []

        with open(self.fasta_file, "r") as file:

            for line in file:
                line = line.strip()

                if not line:
                    continue

                if line.startswith(">"):

                    if name:
                        sequences[name] = "".join(seq).upper()

                    name = line[1:]
                    seq = []

                else:
                    seq.append(line)

            if name:
                sequences[name] = "".join(seq).upper()

        self.sequences = sequences

    # ---------- RUN ANALYSIS ----------
    def run(self):

        print("Working directory:", os.getcwd())

        self.read_fasta()

        print("Sequences loaded:", len(self.sequences))

        with open(self.output_file, "w", newline="") as f:

            writer = csv.writer(f)

            writer.writerow([
                "Protein",
                "Length",
                "D_count", "E_count", "%DE",
                "R_count", "K_count", "%RK",
                "acidic_clustering",
                "basic_clustering",
                "acidic_avg_peptide_length",
                "acidic_min_peptide_length",
                "acidic_max_peptide_length",
                "acidic_MS_fraction",
                "basic_avg_peptide_length",
                "basic_min_peptide_length",
                "basic_max_peptide_length",
                "basic_MS_fraction"
            ])

            for name, seq in self.sequences.items():

                print("Processing:", name)

                D_count = seq.count("D")
                E_count = seq.count("E")
                R_count = seq.count("R")
                K_count = seq.count("K")

                length = len(seq)

                DE_percent = (
                    ((D_count + E_count) / length) * 100
                    if length > 0 else 0
                )

                RK_percent = (
                    ((R_count + K_count) / length) * 100
                    if length > 0 else 0
                )

                acidic_score = self.clustering_distance(
                    seq,
                    ["D", "E"]
                )

                basic_score = self.clustering_distance(
                    seq,
                    ["R", "K"]
                )

                acid_avg, acid_min, acid_max, acid_ms = (
                    self.peptide_lengths(seq, ["D", "E"])
                )

                base_avg, base_min, base_max, base_ms = (
                    self.peptide_lengths(seq, ["R", "K"])
                )

                writer.writerow([
                    name,
                    length,
                    D_count,
                    E_count,
                    round(DE_percent, 2),
                    R_count,
                    K_count,
                    round(RK_percent, 2),
                    round(acidic_score, 3)
                    if acidic_score is not None else "NA",
                    round(basic_score, 3)
                    if basic_score is not None else "NA",
                    round(acid_avg, 2)
                    if acid_avg is not None else "NA",
                    acid_min if acid_min is not None else "NA",
                    acid_max if acid_max is not None else "NA",
                    round(acid_ms, 3)
                    if acid_ms is not None else "NA",
                    round(base_avg, 2)
                    if base_avg is not None else "NA",
                    base_min if base_min is not None else "NA",
                    base_max if base_max is not None else "NA",
                    round(base_ms, 3)
                    if base_ms is not None else "NA"
                ])

        print(f"Results saved to {self.output_file}")


def main():
    fasta_file = "proteins.fasta"
    output_file = "global_results.csv"

    analyzer = ProteinAnalyzer(
        fasta_file,
        output_file
    )

    analyzer.run()


if __name__ == "__main__":
    main()