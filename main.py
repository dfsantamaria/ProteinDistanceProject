import csv
from collections import Counter
from functools import lru_cache
from dataclasses import dataclass
from typing import Dict, Optional, List, Tuple


# Amino acid groups
RESIDUE_GROUPS = {
    "acidic": {"D", "E"},
    "basic": {"R", "K"},
}

# Original peptide statistics
MIN_PEPTIDE_LENGTH = 7
MAX_PEPTIDE_LENGTH = 25

# Digestion coverage peptide-length limits
MIN_COVERAGE_LENGTH = 7
MAX_COVERAGE_LENGTH = 25


@dataclass
class SeqStats:
    counts: Counter
    acidic_nnd: Optional[float]
    basic_nnd: Optional[float]
    acidic: tuple
    basic: tuple


class ProteinAnalyzer:
    def __init__(self, fasta_file: str, output_file: str):
        self.fasta_file = fasta_file
        self.output_file = output_file
        self.sequences: Dict[str, str] = {}

    def read_fasta(self):
        """Read protein sequences from a local FASTA file."""
        self.sequences = {}
        current_name = None
        current_sequence = []

        with open(self.fasta_file, "r", encoding="utf-8") as f:
            for line in f:
                line = line.strip()

                if not line:
                    continue

                if line.startswith(">"):
                    if current_name is not None:
                        self.sequences[current_name] = "".join(
                            current_sequence
                        ).upper()

                    current_name = line[1:].strip()
                    current_sequence = []
                else:
                    if current_name is None:
                        raise ValueError(
                            "The FASTA file must start with a header line "
                            "beginning with '>'."
                        )

                    current_sequence.append(line)

        # Save the final sequence
        if current_name is not None:
            self.sequences[current_name] = "".join(
                current_sequence
            ).upper()

        if not self.sequences:
            raise ValueError(
                f"No protein sequences were found in {self.fasta_file}"
            )

    @staticmethod
    def _nnd(seq: str, residues: set) -> Optional[float]:
        """
        Calculate the mean distance between consecutive occurrences
        of residues in the specified group.
        """
        positions = [
            i for i, aa in enumerate(seq)
            if aa in residues
        ]

        if len(positions) < 2:
            return None

        distances = [
            positions[i + 1] - positions[i]
            for i in range(len(positions) - 1)
        ]

        return sum(distances) / len(distances)

    @staticmethod
    def _peptides(seq: str, residues: set) -> tuple:
        """
        Split a sequence after each residue in the specified group.

        Return average length, minimum length, maximum length, and
        fraction of resulting peptides within the original MIN_PEPTIDE_LENGTH and MAX_PEPTIDE_LENGTH
        amino-acid length range.
        """
        peptide_lengths = []
        start = 0

        for i, aa in enumerate(seq):
            if aa in residues:
                peptide_lengths.append(i + 1 - start)
                start = i + 1

        # Include any remaining sequence after the last cleavage
        if start < len(seq):
            peptide_lengths.append(len(seq) - start)

        valid_lengths = [
            length for length in peptide_lengths
            if MIN_PEPTIDE_LENGTH <= length <= MAX_PEPTIDE_LENGTH
        ]

        if not peptide_lengths or not valid_lengths:
            return None, None, None, 0.0

        average = sum(valid_lengths) / len(valid_lengths)
        minimum = min(valid_lengths)
        maximum = max(valid_lengths)
        fraction = len(valid_lengths) / len(peptide_lengths)

        return average, minimum, maximum, fraction

    @lru_cache(maxsize=100_000)
    def _stats(self, seq: str) -> SeqStats:
        """Calculate the original sequence statistics."""
        counts = Counter(seq)

        acidic_nnd = self._nnd(
            seq, RESIDUE_GROUPS["acidic"]
        )
        basic_nnd = self._nnd(
            seq, RESIDUE_GROUPS["basic"]
        )

        acidic = self._peptides(
            seq, RESIDUE_GROUPS["acidic"]
        )
        basic = self._peptides(
            seq, RESIDUE_GROUPS["basic"]
        )

        return SeqStats(
            counts=counts,
            acidic_nnd=acidic_nnd,
            basic_nnd=basic_nnd,
            acidic=acidic,
            basic=basic,
        )

    @staticmethod
    def _digest_trypsin(seq: str) -> List[Tuple[int, int, str]]:
        """
        Digest with Trypsin.

        Cleaves after K or R unless the next amino acid is P.

        Each fragment is represented as:
            (start_index, end_index, peptide_sequence)

        Indices are zero-based, with end_index exclusive.
        """
        fragments = []
        start = 0

        for i, aa in enumerate(seq):
            if aa in {"K", "R"}:
                next_is_proline = (
                    i + 1 < len(seq) and seq[i + 1] == "P"
                )

                if not next_is_proline:
                    end = i + 1
                    fragments.append((start, end, seq[start:end]))
                    start = end

        if start < len(seq):
            fragments.append((start, len(seq), seq[start:]))

        return fragments

    @staticmethod
    def _digest_gluc(seq: str) -> List[Tuple[int, int, str]]:
        """
        Digest with Glu-C using the requested cleavage rule:
        cleaves after D or E.
        """
        fragments = []
        start = 0

        for i, aa in enumerate(seq):
            if aa in {"D", "E"}:
                end = i + 1
                fragments.append((start, end, seq[start:end]))
                start = end

        if start < len(seq):
            fragments.append((start, len(seq), seq[start:]))

        return fragments

    @staticmethod
    def _apply_missed_cleavages(
        fragments: List[Tuple[int, int, str]],
        missed_cleavages: int,
    ) -> List[Tuple[int, int, str]]:
        """
        Generate peptides allowing the specified number of missed
        cleavages.

        0 missed cleavages: individual fragments
        1 missed cleavage: combine up to 2 adjacent fragments
        2 missed cleavages: combine up to 3 adjacent fragments
        """
        if missed_cleavages < 0:
            raise ValueError(
                "missed_cleavages must be zero or greater."
            )

        peptides = []
        max_fragments = missed_cleavages + 1

        for i in range(len(fragments)):
            combined_sequence = ""

            for j in range(i, min(i + max_fragments, len(fragments))):
                combined_sequence += fragments[j][2]

                peptides.append(
                    (
                        fragments[i][0],
                        fragments[j][1],
                        combined_sequence,
                    )
                )

        return peptides

    @staticmethod
    def _calculate_coverage(
        peptides: List[Tuple[int, int, str]],
        sequence_length: int,
    ) -> float:
        """
        Calculate percentage of sequence positions covered by at least
        one peptide whose length is between 7 and 50 amino acids.
        """
        if sequence_length == 0:
            return 0.0

        covered = bytearray(sequence_length)

        for start, end, peptide in peptides:
            peptide_length = len(peptide)

            if (
                MIN_COVERAGE_LENGTH
                <= peptide_length
                <= MAX_COVERAGE_LENGTH
            ):
                for position in range(start, end):
                    covered[position] = 1

        covered_residues = sum(covered)

        return round(
            covered_residues / sequence_length * 100,
            2,
        )

    def _digest_coverage(
        self,
        seq: str,
        enzyme: str,
    ) -> Dict[str, float]:
        """Calculate coverage for 0, 1, and 2 missed cleavages."""
        enzyme_name = enzyme.strip().lower()

        if enzyme_name == "trypsin":
            fragments = self._digest_trypsin(seq)
            output_prefix = "Trypsin"
        elif enzyme_name in {"gluc", "glu-c", "gluc"}:
            fragments = self._digest_gluc(seq)
            output_prefix = "GluC"
        else:
            raise ValueError(
                f"Unsupported enzyme: {enzyme}. "
                "Use 'Trypsin' or 'GluC'."
            )

        results = {}

        for missed_cleavages in range(3):
            peptides = self._apply_missed_cleavages(
                fragments,
                missed_cleavages,
            )

            coverage = self._calculate_coverage(
                peptides,
                len(seq),
            )

            column_name = (
                f"{output_prefix}_{missed_cleavages}MC_%"
            )
            results[column_name] = coverage

        return results

    def run(self):
        """Analyze all FASTA sequences and write one combined CSV."""
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
            #digestion coverage columns
            "Trypsin_0MC_%",
            "Trypsin_1MC_%",
            "Trypsin_2MC_%",
            "GluC_0MC_%",
            "GluC_1MC_%",
            "GluC_2MC_%",
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
                sequence_length = len(seq)

                D = stats.counts["D"]
                E = stats.counts["E"]
                R = stats.counts["R"]
                K = stats.counts["K"]

                acidic_percentage = (
                    (D + E) / sequence_length * 100
                    if sequence_length else None
                )
                basic_percentage = (
                    (R + K) / sequence_length * 100
                    if sequence_length else None
                )

                trypsin = self._digest_coverage(
                    seq, "Trypsin"
                )
                gluc = self._digest_coverage(
                    seq, "GluC"
                )

                writer.writerow([
                    name,
                    sequence_length,
                    D,
                    E,
                    acidic_percentage,
                    R,
                    K,
                    basic_percentage,
                    stats.acidic_nnd,
                    stats.basic_nnd,
                    *stats.acidic,
                    *stats.basic,
                    trypsin["Trypsin_0MC_%"],
                    trypsin["Trypsin_1MC_%"],
                    trypsin["Trypsin_2MC_%"],
                    gluc["GluC_0MC_%"],
                    gluc["GluC_1MC_%"],
                    gluc["GluC_2MC_%"],
                ])

                #print(f"Processed {name}")

        #print(f"Results written to: {self.output_file}")


if __name__ == "__main__":
    analyzer = ProteinAnalyzer(fasta_file="proteins.fasta", output_file="global_results.csv")
    analyzer.run()
