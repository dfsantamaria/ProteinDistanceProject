import csv
import os
from collections import Counter
from functools import lru_cache
from dataclasses import dataclass
from typing import Dict, Optional


# ==========================================================
# CONFIGURATION
# ==========================================================

RESIDUE_GROUPS = {"acidic": {"D", "E"}, "basic": {"R", "K"},}

MIN_PEPTIDE_LENGTH = 7
MAX_PEPTIDE_LENGTH = 25
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
# TRYPSIN DIGEST
# Cleaves after K or R, except before P
# Coordinates: zero-based, end-exclusive
# ==========================================================

def digest_trypsin(sequence: str) -> list:
    peptides = []
    start = 0

    for i in range(len(sequence) - 1):
        if sequence[i] in {"K", "R"} and sequence[i + 1] != "P":
            peptides.append((start, i + 1, sequence[start:i + 1]))
            start = i + 1
    if start < len(sequence):
        peptides.append((start, len(sequence), sequence[start:]))
    return peptides


# ==========================================================
# GLU-C DIGEST
# Cleaves after D or E, following the supplied rule
# Coordinates: zero-based, end-exclusive
# ==========================================================

def digest_gluc(sequence: str) -> list:
    peptides = []
    start = 0

    for i, aa in enumerate(sequence):
        if aa in {"D", "E"}:
            peptides.append((start, i + 1, sequence[start:i + 1]))
            start = i + 1

    if start < len(sequence):
        peptides.append((start, len(sequence), sequence[start:]))
    return peptides


# ==========================================================
# GENERATE MISSED CLEAVAGES
# ==========================================================

def apply_missed_cleavages(peptides: list, mc: int) -> list:
    if mc < 0:
        raise ValueError("Missed cleavages must be >= 0.")
    result = []
    n = len(peptides)

    for i in range(n):
        start = peptides[i][0]
        current_end = peptides[i][1]
        current_seq = peptides[i][2]

        # Peptide with zero missed cleavages
        result.append((start, current_end, current_seq))

        # Merge adjacent peptides
        for j in range(1, mc + 1):
            if i + j >= n:
                break
            next_peptide = peptides[i + j]
            current_seq += next_peptide[2]
            current_end = next_peptide[1]
            result.append((start, current_end, current_seq))
    return result


# ==========================================================
# CALCULATE PEPTIDE COVERAGE
# Only peptides within the configured length range count
# ==========================================================

def calculate_coverage(sequence_length: int, peptides: list, min_length: int = MIN_PEPTIDE_LENGTH, max_length: int = MAX_PEPTIDE_LENGTH,) -> Optional[float]:

    if sequence_length == 0:
        return None

    covered = set()

    for start, end, peptide in peptides:
        peptide_length = len(peptide)

        if min_length <= peptide_length <= max_length:
            covered.update(range(start, end))

    return round( 100 * len(covered) / sequence_length, 2 )


# ==========================================================
# CALCULATE ORIGINAL PEPTIDE STATISTICS
# ==========================================================

def _peptides(seq: str, residues: set) -> tuple:
    cuts = [i for i, aa in enumerate(seq)  if aa in residues  ]

    prev = 0
    lengths = []

    for c in cuts:
        lengths.append(c + 1 - prev)
        prev = c + 1

    if prev < len(seq):
        lengths.append(len(seq) - prev)

    if not lengths:
        return (None, None, None, None)

    avg = sum(lengths) / len(lengths)
    minimum = min(lengths)
    maximum = max(lengths)

    fraction = sum(
        MIN_PEPTIDE_LENGTH <= length <= MAX_PEPTIDE_LENGTH
        for length in lengths
    ) / len(lengths)

    return avg, minimum, maximum, fraction


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
            raise ValueError(
                "max_missed_cleavages must be >= 0."
            )

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

        name = None
        seq = []

        with open(
            self.fasta_file,
            encoding="utf-8"
        ) as f:

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
                            "FASTA sequence data found before a header."
                        )

                    seq.append(line)

            if name is not None:
                self.sequences[name] = "".join(seq).upper()

    # ======================================================
    # NEAREST-NEIGHBOR DISTANCE
    # ======================================================

    @staticmethod
    def _nnd(seq: str, residues: set) -> Optional[float]:

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

            acidic_nnd=self._nnd(
                seq, RESIDUE_GROUPS["acidic"]
            ),

            basic_nnd=self._nnd(
                seq, RESIDUE_GROUPS["basic"]
            ),

            acidic=_peptides(
                seq, RESIDUE_GROUPS["acidic"]
            ),

            basic=_peptides(
                seq, RESIDUE_GROUPS["basic"]
            ),
        )

    # ======================================================
    # RUN ANALYSIS AND WRITE CSV
    # ======================================================

    def run(self):
        self.read_fasta()

        header = [
            # Original columns
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

            # New digestion and coverage columns
            "trypsin_peptide_count",
            "trypsin_coverage_0mc",
            "trypsin_coverage_max_mc",
            "gluc_peptide_count",
            "gluc_coverage_0mc",
            "gluc_coverage_max_mc",
        ]

        with open(
            self.output_file,
            "w",
            newline="",
            encoding="utf-8"
        ) as f:

            writer = csv.writer(f)
            writer.writerow(header)

            for name, seq in self.sequences.items():

                stats = self._stats(seq)
                counts = stats.counts
                length = len(seq)

                # Original residue counts
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

                # Digest with zero missed cleavages
                trypsin_base = digest_trypsin(seq)
                gluc_base = digest_gluc(seq)

                # Generate peptides with missed cleavages
                trypsin_max = apply_missed_cleavages(
                    trypsin_base,
                    self.max_missed_cleavages
                )

                gluc_max = apply_missed_cleavages(
                    gluc_base,
                    self.max_missed_cleavages
                )

                # Calculate coverage
                trypsin_coverage_0mc = calculate_coverage(
                    length, trypsin_base
                )

                trypsin_coverage_max_mc = calculate_coverage(
                    length, trypsin_max
                )

                gluc_coverage_0mc = calculate_coverage(
                    length, gluc_base
                )

                gluc_coverage_max_mc = calculate_coverage(
                    length, gluc_max
                )

                # Write original and new CSV fields
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

                    len(trypsin_base),
                    trypsin_coverage_0mc,
                    trypsin_coverage_max_mc,

                    len(gluc_base),
                    gluc_coverage_0mc,
                    gluc_coverage_max_mc,
                ])


# ==========================================================
# MAIN
# ==========================================================

if __name__ == "__main__":
    analyzer = ProteinAnalyzer("proteins.fasta","global_results.csv", max_missed_cleavages=2,)
    analyzer.run()

