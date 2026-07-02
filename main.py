import csv
import os
from typing import Dict, List, Tuple, Optional


# =========================
# BIO CONFIG
# =========================

RESIDUE_GROUPS = {
    "acidic": {"D", "E"},
    "basic": {"R", "K"}
}

MIN_PEPTIDE_LENGTH = 7
MAX_PEPTIDE_LENGTH = 25


# =========================
# ANALYZER
# =========================

class ProteinAnalyzer:

    def __init__(self, fasta_file: str, output_file: str):
        self.fasta_file = fasta_file
        self.output_file = output_file
        self.sequences: Dict[str, str] = {}

        self.features: List[Tuple[str, callable]] = []

    # =========================
    # FEATURE REGISTRY
    # =========================

    def feature(self, name: str):
        def wrapper(func):
            self.features.append((name, func))
            return func
        return wrapper

    # =========================
    # CORE METHODS
    # =========================

    def nearest_neighbor_distance(
        self,
        sequence: str,
        residues: set[str]
    ) -> Optional[float]:

        positions = [
            i + 1 for i, aa in enumerate(sequence)
            if aa in residues
        ]

        if len(positions) < 2:
            return None

        def dist(a: int, b: int) -> int:
            return max(abs(a - b) - 1, 0)

        values: List[int] = []

        for i, pos in enumerate(positions):
            candidates = []

            if i > 0:
                candidates.append(dist(pos, positions[i - 1]))

            if i < len(positions) - 1:
                candidates.append(dist(pos, positions[i + 1]))

            values.append(min(candidates))

        return sum(values) / len(values)

    # =========================
    # PEPTIDE LENGTHS (CACHE SHARED)
    # =========================

    def peptide_lengths(
        self,
        sequence: str,
        residues: set[str]
    ) -> Tuple[Optional[float], Optional[int], Optional[int], Optional[float]]:

        cuts = [i for i, aa in enumerate(sequence) if aa in residues]

        lengths: List[int] = []
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
            if MIN_PEPTIDE_LENGTH <= l <= MAX_PEPTIDE_LENGTH
        ]

        frac = len(ms) / len(lengths)

        return avg, mn, mx, frac

    def _peptide_cache(self, seq: str, group: str):
        return self.peptide_lengths(seq, RESIDUE_GROUPS[group])

    # =========================
    # FASTA
    # =========================

    def read_fasta(self) -> None:
        if not os.path.exists(self.fasta_file):
            raise FileNotFoundError(f"FASTA not found: {self.fasta_file}")

        data: Dict[str, str] = {}
        name = None
        seq: List[str] = []

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
    # FEATURES
    # =========================

    def f_name(self, seq: str, name: str) -> str:
        return name

    def f_length(self, seq: str, name: str) -> int:
        return len(seq)

    def f_d(self, seq: str, name: str) -> int:
        return seq.count("D")

    def f_e(self, seq: str, name: str) -> int:
        return seq.count("E")

    def f_de(self, seq: str, name: str) -> Optional[float]:
        if len(seq) == 0:
            return None
        return ((seq.count("D") + seq.count("E")) / len(seq)) * 100

    def f_r(self, seq: str, name: str) -> int:
        return seq.count("R")

    def f_k(self, seq: str, name: str) -> int:
        return seq.count("K")

    def f_rk(self, seq: str, name: str) -> Optional[float]:
        if len(seq) == 0:
            return None
        return ((seq.count("R") + seq.count("K")) / len(seq)) * 100

    # =========================
    # NEAREST NEIGHBOR FEATURES
    # =========================

    def f_acidic_nearest_neighbor_distance(self, seq: str, name: str):
        return self.nearest_neighbor_distance(seq, RESIDUE_GROUPS["acidic"])

    def f_basic_nearest_neighbor_distance(self, seq: str, name: str):
        return self.nearest_neighbor_distance(seq, RESIDUE_GROUPS["basic"])

    # =========================
    # PEPTIDE FEATURES
    # =========================

    def f_acidic_avg(self, seq: str, name: str):
        return self._peptide_cache(seq, "acidic")[0]

    def f_acidic_min(self, seq: str, name: str):
        return self._peptide_cache(seq, "acidic")[1]

    def f_acidic_max(self, seq: str, name: str):
        return self._peptide_cache(seq, "acidic")[2]

    def f_acidic_fraction(self, seq: str, name: str):
        return self._peptide_cache(seq, "acidic")[3]

    def f_basic_avg(self, seq: str, name: str):
        return self._peptide_cache(seq, "basic")[0]

    def f_basic_min(self, seq: str, name: str):
        return self._peptide_cache(seq, "basic")[1]

    def f_basic_max(self, seq: str, name: str):
        return self._peptide_cache(seq, "basic")[2]

    def f_basic_fraction(self, seq: str, name: str):
        return self._peptide_cache(seq, "basic")[3]

    # =========================
    # RUN
    # =========================

    def run(self) -> None:
        self.read_fasta()

        self.feature("Protein")(self.f_name)
        self.feature("Length")(self.f_length)

        self.feature("D_count")(self.f_d)
        self.feature("E_count")(self.f_e)
        self.feature("%DE")(self.f_de)

        self.feature("R_count")(self.f_r)
        self.feature("K_count")(self.f_k)
        self.feature("%RK")(self.f_rk)

        self.feature("acidic_nearest_neighbor_distance")(
            self.f_acidic_nearest_neighbor_distance
        )
        self.feature("basic_nearest_neighbor_distance")(
            self.f_basic_nearest_neighbor_distance
        )

        self.feature("acidic_avg")(self.f_acidic_avg)
        self.feature("acidic_min")(self.f_acidic_min)
        self.feature("acidic_max")(self.f_acidic_max)
        self.feature("acidic_fraction")(self.f_acidic_fraction)

        self.feature("basic_avg")(self.f_basic_avg)
        self.feature("basic_min")(self.f_basic_min)
        self.feature("basic_max")(self.f_basic_max)
        self.feature("basic_fraction")(self.f_basic_fraction)

        header = [name for name, _ in self.features]

        with open(self.output_file, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(header)

            for name, seq in self.sequences.items():
                row = [
                    func(seq, name)
                    for _, func in self.features
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