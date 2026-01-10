from Bio import SeqIO
import sys
import os
from collections import Counter
import matplotlib.pyplot as plt 
import re


class SequenceAnalyzer:
    def __init__(self, file_path):
        """
        Initializes the analyzer and loads sequences from a FASTA file.

        Args:
            file_path (str): The path to the DNA sequence file (FASTA format).
        """
        self.records = list(SeqIO.parse(file_path, "fasta"))
        self.file_name = file_path.split('.')[0]

    def get_codons(self, sequence):
        """
        Slices a DNA string into a list of 3-character codons.

        Args:
            sequence (str): The DNA sequence string to be parsed.

        Returns:
            List[str]: A list of triplet codons.
        """
        #List comprehension for dummies: Slice the string (Sequence) from including i up to but not including i+3 then the i is iterated from 0 up to the length of the string
        #The range operator also has an optional parameter "step" that will change how much the i jumps per round. I.E. round 1 i=0 slice [0:4] then next round i=3 slice [3:7]
        # len(sequence) - 2 to avoid a codon < 3
        return [sequence[i:i+3] for i in range(0, len(sequence) - 2 ,3)]

    def validate_sequence(self, sequence):
        """
        Checks if a sequence contains only valid nucleotides (A, T, C, G).

        Args:
            sequence (str): The sequence string to validate.

        Returns:
            bool: True if valid, False if it contains unknown characters.
        """
        sequence = str(sequence)
        matches = re.search(r'[^ATCG]', sequence)
        if matches: return False
        else: return True

    def plot_gc_content(self):
        """
        Calculates GC content for valid sequences and displays a bar chart.
        """
        gc_data = {}
        for record in self.records:
           if self.validate_sequence(record.seq):
            guanine_count = record.seq.count('G')
            cytosine_count = record.seq.count('C')
            total_count = self.sequence_count(record)
            gc_content = (guanine_count + cytosine_count)/total_count * 100
            gc_data[record.id] = gc_content
           else: print(f"Skipping {record.id}: Invalid sequence")
        plt.bar(x=gc_data.keys(), height=gc_data.values())
        plt.ylabel('GC Content (%)')
        plt.xlabel('Sequence ID')
        plt.title('GC content')
        plt.show()

    def sequence_count(self, record):
      return len(record.seq)

    def codon_usage(self):
       """
        Analyzes the frequency of each codon per sequence.

        Returns:
            Dict[str, Counter]: A dictionary mapping Sequence IDs to their codon counts.
        """
       codon_count_dict = {}
       for record in self.records:
          if self.validate_sequence(record.seq):
             codon_count_dict[record.id] = Counter(self.get_codons(record.seq))
          else: print(f"Skipping {record.id}: Invalid sequence")
       return codon_count_dict



    def dna_to_protein(self):
      """
        Translates all valid DNA sequences in the file into protein strings.

        Returns:
            Dict[str, str]: A dictionary mapping Sequence IDs to protein sequences.
        """
      translations = {}
      for record in self.records:
         if self.validate_sequence(record.seq):
            translations[record.id] = record.seq.translate()
         else: print(f"Skipping {record.id}: Invalid sequence")
      return translations
    
       

# -------------------------------------------------------------------------
# EXECUTION BLOCK
# -------------------------------------------------------------------------
if __name__ == "__main__":
    # Once refactored, your main execution will look like this:
    analyzer = SequenceAnalyzer("TP53seq.fasta")
    analyzer.plot_gc_content()
    print(analyzer.dna_to_protein())
    print(analyzer.codon_usage())