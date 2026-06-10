import csv
import os


# =========================
# BIO CONFIG
# =========================

RESIDUE_GROUPS = {
    "acidic": {"D", "E"},
    "basic": {"R", "K"}
}

MS_MIN_LENGTH = 7
MS_MAX_LENGTH = 25


# =========================
# FEATURE REGISTRY
# =========================

FEATURES = []


def feature(name):
    """Decorator to register a feature."""
    def wrapper(func):
        FEATURES.append((name, func))
        return func
    return wrapper


# =========================
# ANALYZER
# =========================

class ProteinAnalyzer:

    def __init__(self, fasta_file, output_file):
        self.fasta_file = fasta_file
        self.output_file = output_file
        self.sequences = {}

    # ---------- UTILS ----------
    @staticmethod
    def format_value(value, digits):
        return round(value, digits) if value is not None else "NA"

    # ---------- CORE METHODS ----------
    def clustering_distance(self, sequence, residues):
        residues = set(residues)

        positions = [
            i + 1 for i, aa in enumerate(sequence)
            if aa in residues
        ]

        if len(positions) < 2:
            return None

        def dist(a, b):
            return max(abs(a - b) - 1, 0)

        values = []

        for i, pos in enumerate(positions):
            candidates = []

            if i > 0:
                candidates.append(dist(pos, positions[i - 1]))

            if i < len(positions) - 1:
                candidates.append(dist(pos, positions[i + 1]))

            values.append(min(candidates))

        return sum(values) / len(values)

    def peptide_lengths(self, sequence, residues):
        residues = set(residues)

        cuts = [i for i, aa in enumerate(sequence) if aa in residues]

        lengths = []
        prev = 0

        for c in cuts:
            l = c + 1 - prev
            if l > 0:
                lengths.append(l)
            prev = c + 1

        last = len(sequence) - prev
        if last > 0:
            lengths.append(last)

        if not lengths:
            return None, None, None, None

        avg = sum(lengths) / len(lengths)
        mn = min(lengths)
        mx = max(lengths)

        ms = [
            l for l in lengths
            if MS_MIN_LENGTH <= l <= MS_MAX_LENGTH
        ]

        frac = len(ms) / len(lengths)

        return avg, mn, mx, frac

    # ---------- FASTA ----------
    def read_fasta(self):

        if not os.path.exists(self.fasta_file):
            raise FileNotFoundError(
                f"FASTA not found: {self.fasta_file}"
            )

        data = {}
        name = None
        seq = []

        with open(self.fasta_file) as f:

            for line in f:
                line = line.strip()
                if not line:
                    continue

                if line.startswith(">"):

                    if name:
                        data[name] = "".join(seq).upper()

                    name = line[1:]
                    seq = []

                else:
                    seq.append(line)

            if name:
                data[name] = "".join(seq).upper()

        self.sequences = data

    # =========================
    # FEATURE DEFINITIONS
    # =========================

    @feature("Protein")
    def f_name(self, seq, name):
        return name

    @feature("Length")
    def f_length(self, seq, name):
        return len(seq)

    @feature("D_count")
    def f_d(self, seq, name):
        return seq.count("D")

    @feature("E_count")
    def f_e(self, seq, name):
        return seq.count("E")

    @feature("%DE")
    def f_de(self, seq, name):
        l = len(seq)
        return ((seq.count("D") + seq.count("E")) / l) * 100 if l else 0

    @feature("R_count")
    def f_r(self, seq, name):
        return seq.count("R")

    @feature("K_count")
    def f_k(self, seq, name):
        return seq.count("K")

    @feature("%RK")
    def f_rk(self, seq, name):
        l = len(seq)
        return ((seq.count("R") + seq.count("K")) / l) * 100 if l else 0

    @feature("acidic_clustering")
    def f_acid_cluster(self, seq, name):
        return self.clustering_distance(seq, RESIDUE_GROUPS["acidic"])

    @feature("basic_clustering")
    def f_basic_cluster(self, seq, name):
        return self.clustering_distance(seq, RESIDUE_GROUPS["basic"])

    @feature("acidic_avg_peptide_length")
    def f_acid_avg(self, seq, name):
        return self.peptide_lengths(seq, RESIDUE_GROUPS["acidic"])[0]

    @feature("acidic_min_peptide_length")
    def f_acid_min(self, seq, name):
        return self.peptide_lengths(seq, RESIDUE_GROUPS["acidic"])[1]

    @feature("acidic_max_peptide_length")
    def f_acid_max(self, seq, name):
        return self.peptide_lengths(seq, RESIDUE_GROUPS["acidic"])[2]

    @feature("acidic_MS_fraction")
    def f_acid_ms(self, seq, name):
        return self.peptide_lengths(seq, RESIDUE_GROUPS["acidic"])[3]

    @feature("basic_avg_peptide_length")
    def f_basic_avg(self, seq, name):
        return self.peptide_lengths(seq, RESIDUE_GROUPS["basic"])[0]

    @feature("basic_min_peptide_length")
    def f_basic_min(self, seq, name):
        return self.peptide_lengths(seq, RESIDUE_GROUPS["basic"])[1]

    @feature("basic_max_peptide_length")
    def f_basic_max(self, seq, name):
        return self.peptide_lengths(seq, RESIDUE_GROUPS["basic"])[2]

    @feature("basic_MS_fraction")
    def f_basic_ms(self, seq, name):
        return self.peptide_lengths(seq, RESIDUE_GROUPS["basic"])[3]

    # =========================
    # RUN
    # =========================

    def run(self):

        self.read_fasta()

        header = [name for name, _ in FEATURES]

        with open(self.output_file, "w", newline="") as f:

            writer = csv.writer(f)
            writer.writerow(header)

            for name, seq in self.sequences.items():

                row = [
                    func(self, seq, name)
                    for _, func in FEATURES
                ]

                writer.writerow(row)


# =========================
# MAIN
# =========================

def main():
    analyzer = ProteinAnalyzer(
        "proteins.fasta",
        "global_results.csv"
    )
    analyzer.run()


if __name__ == "__main__":
    main()