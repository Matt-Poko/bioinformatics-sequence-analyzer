from Bio import SeqIO
from Bio.Seq import Seq
from collections import Counter
import matplotlib.pyplot as plt 
import random

class SequenceAnalyzer:
    def __init__(self, file_path):
        """
        Initializes the analyzer and loads sequences from a FASTA file.

        Args:
            file_path (str): The path to the DNA sequence file (FASTA format).
        """
        self.records = list(SeqIO.parse(file_path, "fasta"))
        self.file_name = file_path.split('.')[0]

    def get_codons(self, sequence_id):
        """
        Slices a DNA string into a list of 3-character codons.

        Args:
            sequence (str): The DNA sequence string to be parsed.

        Returns:
            List[str]: A list of triplet codons.
        """
        target_record = next(record for record in self.records if record.id == sequence_id)
        sequence = target_record.seq
        if self.validate_sequence(sequence):
            return [sequence[i:i+3] for i in range(0, len(sequence) - 2 ,3)]
        else: print(f"Skipping {target_record.id}: Invalid sequence")
        

    def validate_sequence(self, sequence):
        """
        Checks if a sequence contains only valid nucleotides (A, T, C, G).

        Args:
            sequence (str): The sequence string to validate.

        Returns:
            bool: True if valid, False if it contains unknown characters.
        """
        sequence = str(sequence)
        sequence = sequence.strip().upper()
        allowed = set("ATCG")
        if not set(sequence) <= allowed: return False
        else: return True

    def plot_gc_content(self, file_name=None):
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
        if file_name:
           plt.savefig(file_name)
           plt.close()
        else: plt.show()

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
    
    def point_mutation(self, sequence_id, position, new_base):
        """Manually changes one base at a specific index."""
        target_record = next(record for record in self.records if record.id == sequence_id)
        base_lst = list(target_record.seq)
        base_lst[position] = new_base
        new_seq = "".join(base_lst)
        target_record.seq = Seq(new_seq)
        return new_seq
    

    def simulate_random_mutation(self, sequence_id):
        """Randomly swaps one base to see the effect on the protein."""
        target_record = next(record for record in self.records if record.id == sequence_id)
        position = random.randrange(0, len(target_record.seq))
        base_lst = list(target_record.seq)
        new_base = ""
        bases = ['A','G','T','C']
        while new_base == "":
           temp_base = random.choice(bases)
           if temp_base != base_lst[position]:
              new_base = temp_base
           else: continue
        base_lst[position] = new_base
        new_seq = "".join(base_lst)
        target_record.seq = Seq(new_seq)
        return new_seq
    
        
