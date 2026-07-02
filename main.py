import csv
import os
from collections import Counter
from functools import lru_cache
from dataclasses import dataclass
from typing import Dict, Optional

RESIDUE_GROUPS={"acidic":{"D","E"},"basic":{"R","K"}}
MIN_PEPTIDE_LENGTH=7
MAX_PEPTIDE_LENGTH=25

@dataclass(frozen=True)
class SeqStats:
    counts: Counter
    acidic_nnd: Optional[float]
    basic_nnd: Optional[float]
    acidic: tuple
    basic: tuple

class ProteinAnalyzer:
    def __init__(self,fasta_file:str,output_file:str):
        self.fasta_file=fasta_file
        self.output_file=output_file
        self.sequences:Dict[str,str]={}

    def read_fasta(self):
        if not os.path.exists(self.fasta_file):
            raise FileNotFoundError(self.fasta_file)
        name=None;seq=[]
        with open(self.fasta_file) as f:
            for line in f:
                line=line.strip()
                if not line: continue
                if line.startswith(">"):
                    if name:self.sequences[name]="".join(seq).upper()
                    name=line[1:];seq=[]
                else: seq.append(line)
            if name:self.sequences[name]="".join(seq).upper()

    @staticmethod
    def _nnd(seq,res):
        p=[i+1 for i,a in enumerate(seq) if a in res]
        if len(p)<2:return None
        d=[]
        for i,x in enumerate(p):
            c=[]
            if i:c.append(max(x-p[i-1]-1,0))
            if i<len(p)-1:c.append(max(p[i+1]-x-1,0))
            d.append(min(c))
        return sum(d)/len(d)

    @lru_cache(maxsize=None)
    def _stats(self,seq):
        return SeqStats(
            counts=Counter(seq),
            acidic_nnd=self._nnd(seq,RESIDUE_GROUPS["acidic"]),
            basic_nnd=self._nnd(seq,RESIDUE_GROUPS["basic"]),
            acidic=_peptides(seq, RESIDUE_GROUPS["acidic"]),
            basic=_peptides(seq, RESIDUE_GROUPS["basic"])
        )

    def run(self):
        self.read_fasta()
        header=["Protein","Length","D_count","E_count","%DE","R_count","K_count","%RK",
                "acidic_nearest_neighbor_distance","basic_nearest_neighbor_distance",
                "acidic_avg","acidic_min","acidic_max","acidic_fraction",
                "basic_avg","basic_min","basic_max","basic_fraction"]
        with open(self.output_file,"w",newline="") as f:
            w=csv.writer(f);w.writerow(header)
            for name,seq in self.sequences.items():
                s=self._stats(seq)
                D,E,R,K=s.counts["D"],s.counts["E"],s.counts["R"],s.counts["K"]
                w.writerow([
                    name,len(seq),D,E,(D+E)/len(seq)*100 if seq else None,
                    R,K,(R+K)/len(seq)*100 if seq else None,
                    s.acidic_nnd,s.basic_nnd,
                    *s.acidic,*s.basic
                ])


def _peptides(seq,res):
    cuts=[i for i,a in enumerate(seq) if a in res]
    prev=0;L=[]
    for c in cuts:
        L.append(c+1-prev);prev=c+1
    if prev<len(seq):L.append(len(seq)-prev)
    avg=sum(L)/len(L)
    return (avg,min(L),max(L),sum(MIN_PEPTIDE_LENGTH<=x<=MAX_PEPTIDE_LENGTH for x in L)/len(L))


if __name__=="__main__":
    ProteinAnalyzer("proteins.fasta","global_results.csv").run()
