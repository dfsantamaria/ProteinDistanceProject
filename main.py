import csv
import os
from collections import Counter

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


def reset_features():
    FEATURES.clear()


# =========================
# ANALYZER
# =========================

class ProteinAnalyzer:

    def __init__(self, fasta_file, output_file):
        self.fasta_file = fasta_file
        self.output_file = output_file
        self.sequences = {}
        self._cache = {}   # FIX caching


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
                line = line.strip().replace("\r", "")
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
    # CACHE HELPERS
    # =========================

    def _acidic_cache(self, seq):
        key = ("acidic", seq)
        if key not in self._cache:
            self._cache[key] = self.peptide_lengths(
                seq, RESIDUE_GROUPS["acidic"]
            )
        return self._cache[key]


    def _basic_cache(self, seq):
        key = ("basic", seq)
        if key not in self._cache:
            self._cache[key] = self.peptide_lengths(
                seq, RESIDUE_GROUPS["basic"]
            )
        return self._cache[key]


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
    def f_d(self, seq, name, aa):
        return aa.get("D", 0)


    @feature("E_count")
    def f_e(self, seq, name, aa):
        return aa.get("E", 0)


    @feature("%DE")
    def f_de(self, seq, name, aa):
        l = len(seq)
        if l == 0:
            return None
        return ((aa.get("D", 0) + aa.get("E", 0)) / l) * 100


    @feature("R_count")
    def f_r(self, seq, name, aa):
        return aa.get("R", 0)


    @feature("K_count")
    def f_k(self, seq, name, aa):
        return aa.get("K", 0)


    @feature("%RK")
    def f_rk(self, seq, name, aa):
        l = len(seq)
        if l == 0:
            return None
        return ((aa.get("R", 0) + aa.get("K", 0)) / l) * 100


    # =========================
    # CACHED PEPTIDE FEATURES
    # =========================

    @feature("acidic_clustering")
    def f_acid_cluster(self, seq, name):
        return self.clustering_distance(seq, RESIDUE_GROUPS["acidic"])


    @feature("basic_clustering")
    def f_basic_cluster(self, seq, name):
        return self.clustering_distance(seq, RESIDUE_GROUPS["basic"])


    @feature("acidic_avg_peptide_length")
    def f_acid_avg(self, seq, name):
        return self._acidic_cache(seq)[0]


    @feature("acidic_min_peptide_length")
    def f_acid_min(self, seq, name):
        return self._acidic_cache(seq)[1]


    @feature("acidic_max_peptide_length")
    def f_acid_max(self, seq, name):
        return self._acidic_cache(seq)[2]


    @feature("acidic_MS_fraction")
    def f_acid_ms(self, seq, name):
        return self._acidic_cache(seq)[3]


    @feature("basic_avg_peptide_length")
    def f_basic_avg(self, seq, name):
        return self._basic_cache(seq)[0]


    @feature("basic_min_peptide_length")
    def f_basic_min(self, seq, name):
        return self._basic_cache(seq)[1]


    @feature("basic_max_peptide_length")
    def f_basic_max(self, seq, name):
        return self._basic_cache(seq)[2]


    @feature("basic_MS_fraction")
    def f_basic_ms(self, seq, name):
        return self._basic_cache(seq)[3]


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

                aa = Counter(seq)   # OPTIMIZATION FIX

                row = [
                    func(self, seq, name, aa) if "aa" in func.__code__.co_varnames
                    else func(self, seq, name)
                    for _, func in FEATURES
                ]

                writer.writerow(row)


# =========================
# MAIN
# =========================

def main():
    reset_features()

    analyzer = ProteinAnalyzer(
        "proteins.fasta",
        "global_results.csv"
    )
    analyzer.run()


if __name__ == "__main__":
    main()