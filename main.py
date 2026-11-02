
import csv
import os
from collections import Counter
from functools import lru_cache
from dataclasses import dataclass
from typing import Dict, Optional


# ==========================================================
# CONFIGURATION
# ==========================================================

RESIDUE_GROUPS = {
    "acidic": {"D", "E"},
    "basic": {"R", "K"},
}

MIN_PEPTIDE_LENGTH = 7
MAX_PEPTIDE_LENGTH = 25


# Peptide statistics limits
MIN_STATS_LENGTH = 7
MAX_STATS_LENGTH = 25
MAX_MISSED_CLEAVAGES = 2

# ==========================================================
# SEQUENCE STATISTICS
# ==========================================================

@dataclass(frozen=True)
class SeqStats:
    counts: Counter
    acidic_nnd: Optional[float]
    basic_nnd: Optional[float]
    acidic: tuple
    basic: tuple


# ==========================================================
# PROTEIN ANALYZER
# ==========================================================

class ProteinAnalyzer:

    def __init__(
        self,
        fasta_file: str,
        output_file: str,
        max_missed_cleavages: int = MAX_MISSED_CLEAVAGES,
    ):
        if max_missed_cleavages < 0:
            raise ValueError("Maximum missed cleavages must be >= 0.")

        self.fasta_file = fasta_file
        self.output_file = output_file
        self.max_missed_cleavages = max_missed_cleavages
        self.sequences: Dict[str, str] = {}

    # ======================================================
    # READ FASTA FILE
    # ======================================================

    def read_fasta(self):
        if not os.path.exists(self.fasta_file):
            raise FileNotFoundError(self.fasta_file)

        self.sequences.clear()
        name = None
        seq = []

        with open(self.fasta_file, encoding="utf-8") as f:
            for line in f:
                line = line.strip()

                if not line:
                    continue

                if line.startswith(">"):
                    if name is not None:
                        self.sequences[name] = "".join(seq).upper()

                    name = line[1:].strip()
                    seq = []

                else:
                    if name is None:
                        raise ValueError(
                             "The FASTA file must start with a header line beginning with '>'."
                        )

                    seq.append(line)

            if name is not None:
                self.sequences[name] = "".join(seq).upper()

    # ======================================================
    # TRYPSIN DIGEST
    # Cleaves after K or R, except before P
    # Coordinates: zero-based, end-exclusive
    # ======================================================

    @staticmethod
    def digest_trypsin(sequence: str) -> list:
        peptides = []
        start = 0

        for i in range(len(sequence) - 1):
            if (
                sequence[i] in {"K", "R"}
                and sequence[i + 1] != "P"
            ):
                peptides.append(
                    (start, i + 1, sequence[start:i + 1])
                )
                start = i + 1

        if start < len(sequence):
            peptides.append(
                (start, len(sequence), sequence[start:])
            )

        return peptides

    # ======================================================
    # GLU-C DIGEST
    # Cleaves after D or E
    # ======================================================

    @staticmethod
    def digest_gluc(sequence: str) -> list:
        peptides = []
        start = 0

        for i, aa in enumerate(sequence):
            if aa in {"D", "E"}:
                peptides.append(
                    (start, i + 1, sequence[start:i + 1])
                )
                start = i + 1

        if start < len(sequence):
            peptides.append(
                (start, len(sequence), sequence[start:])
            )

        return peptides

    # ======================================================
    # GENERATE PEPTIDES WITH MISSED CLEAVAGES
    # ======================================================

    @staticmethod
    def apply_missed_cleavages(
        peptides: list,
        mc: int,
    ) -> list:

        if mc < 0:
            raise ValueError("Missed cleavages must be >= 0.")

        result = []
        n = len(peptides)

        for i in range(n):
            start = peptides[i][0]
            current_end = peptides[i][1]
            current_seq = peptides[i][2]

            # Zero missed cleavages
            result.append(
                (start, current_end, current_seq)
            )

            # Merge adjacent peptides
            for j in range(1, mc + 1):
                if i + j >= n:
                    break

                next_peptide = peptides[i + j]
                current_seq += next_peptide[2]
                current_end = next_peptide[1]

                result.append((start, current_end, current_seq))

        return result

    # ======================================================
    # CALCULATE PEPTIDE COVERAGE
    # ======================================================

    @staticmethod
    def calculate_coverage(
        sequence_length: int,
        peptides: list,
        min_length: int = MIN_PEPTIDE_LENGTH,
        max_length: int = MAX_PEPTIDE_LENGTH,
    ) -> Optional[float]:

        if sequence_length == 0:
            return None

        covered = set()

        for start, end, peptide in peptides:
            peptide_length = len(peptide)

            if min_length <= peptide_length <= max_length:
                covered.update(range(start, end))

        return round(
            100 * len(covered) / sequence_length,
            2
        )

    # ======================================================
    # CALCULATE ORIGINAL PEPTIDE STATISTICS
    # ======================================================

    @staticmethod
    def calculate_peptide_statistics(
            seq: str,
            residues: set,
    ) -> tuple:

        cuts = [
            i for i, aa in enumerate(seq[:-1])
            if aa in residues and seq[i + 1] != "P"
        ]

        # Include the last residue if it belongs to the cleavage group
        if seq and seq[-1] in residues:
            cuts.append(len(seq) - 1)

        prev = 0
        lengths = []

        for c in cuts:
            lengths.append(c + 1 - prev)
            prev = c + 1

        if prev < len(seq):
            lengths.append(len(seq) - prev)

        if not lengths:
            return (None, None, None, 0.0)

        # Average, minimum and maximum - no length restrictions
        avg = sum(lengths) / len(lengths)
        minimum = min(lengths)
        maximum = max(lengths)

        # Fraction restricted to 7-25 amino acids
        valid_lengths = [
            length for length in lengths
            if MIN_STATS_LENGTH <= length <= MAX_STATS_LENGTH
        ]

        fraction = len(valid_lengths) / len(lengths)

        return avg, minimum, maximum, fraction



    # ======================================================
    # NEAREST-NEIGHBOR DISTANCE
    # ======================================================

    @staticmethod
    def calculate_nnd(
        seq: str,
        residues: set,
    ) -> Optional[float]:

        positions = [
            i + 1
            for i, aa in enumerate(seq)
            if aa in residues
        ]

        if len(positions) < 2:
            return None

        distances = []

        for i, position in enumerate(positions):
            neighbors = []

            if i > 0:
                neighbors.append(
                    position - positions[i - 1] - 1
                )

            if i < len(positions) - 1:
                neighbors.append(
                    positions[i + 1] - position - 1
                )

            distances.append(min(neighbors))

        return sum(distances) / len(distances)

    # ======================================================
    # CACHE SEQUENCE STATISTICS
    # ======================================================

    @lru_cache(maxsize=None)
    def _stats(self, seq: str) -> SeqStats:
        return SeqStats(
            counts=Counter(seq),

            acidic_nnd=self.calculate_nnd(
                seq, RESIDUE_GROUPS["acidic"]
            ),

            basic_nnd=self.calculate_nnd(
                seq, RESIDUE_GROUPS["basic"]
            ),

            acidic=self.calculate_peptide_statistics(
                seq, RESIDUE_GROUPS["acidic"]
            ),

            basic=self.calculate_peptide_statistics(
                seq, RESIDUE_GROUPS["basic"]
            )
        )

    # ======================================================
    # CALCULATE COVERAGE FOR EACH MISSED-CLEAVAGE LEVEL
    # ======================================================

    def coverage_by_missed_cleavages(
        self,
        sequence: str,
        base_peptides: list,
    ) -> list:

        coverages = []

        for mc in range(self.max_missed_cleavages + 1):
            peptides = self.apply_missed_cleavages(
                base_peptides, mc
            )

            coverage = self.calculate_coverage(
                len(sequence), peptides
            )

            coverages.append(coverage)

        return coverages

    # ======================================================
    # RUN ANALYSIS AND WRITE CSV
    # ======================================================

    def run(self):
        self.read_fasta()

        header = [
            "Protein",
            "Length",
            "D_count",
            "E_count",
            "%DE",
            "R_count",
            "K_count",
            "%RK",
            "acidic_nearest_neighbor_distance",
            "basic_nearest_neighbor_distance",
            "acidic_avg",
            "acidic_min",
            "acidic_max",
            "acidic_fraction",
            "basic_avg",
            "basic_min",
            "basic_max",
            "basic_fraction",

            *[
                f"Trypsin_{mc}MC_%"
                for mc in range(self.max_missed_cleavages + 1)
            ],

            *[
                f"GluC_{mc}MC_%"
                for mc in range(self.max_missed_cleavages + 1)
            ]
        ]

        with open(
            self.output_file,
            "w",
            newline="",
            encoding="utf-8",
        ) as f:

            writer = csv.writer(f)
            writer.writerow(header)

            for name, seq in self.sequences.items():
                stats = self._stats(seq)
                counts = stats.counts
                length = len(seq)

                D = counts["D"]
                E = counts["E"]
                R = counts["R"]
                K = counts["K"]

                percent_DE = (
                    (D + E) / length * 100
                    if length else None
                )

                percent_RK = (
                    (R + K) / length * 100
                    if length else None
                )

                # Digest protein
                trypsin_base = self.digest_trypsin(seq)
                gluc_base = self.digest_gluc(seq)

                # Coverage for 0, 1, ..., max missed cleavages
                trypsin_coverages = (
                    self.coverage_by_missed_cleavages(
                        seq, trypsin_base
                    )
                )

                gluc_coverages = (
                    self.coverage_by_missed_cleavages(
                        seq, gluc_base
                    )
                )

                writer.writerow([
                    name,
                    length,
                    D,
                    E,
                    percent_DE,
                    R,
                    K,
                    percent_RK,

                    stats.acidic_nnd,
                    stats.basic_nnd,

                    *stats.acidic,
                    *stats.basic,

                    *trypsin_coverages,
                    *gluc_coverages
                ])


# ==========================================================
# MAIN
# ==========================================================

if __name__ == "__main__":
    analyzer = ProteinAnalyzer("proteins.fasta", "global_results.csv",  max_missed_cleavages=2)

    analyzer.run()
