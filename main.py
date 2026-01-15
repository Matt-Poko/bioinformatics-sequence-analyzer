import argparse
import sys
from src import SequenceAnalyzer

def run():
    parser = argparse.ArgumentParser(prog="Sequence Analyzer", description="Analyzes gene sequences from fasta files")
    parser.add_argument('--input', required=True, help="Path to desired fasta file")
    parser.add_argument('--operation', required=True, 
                        help='Operation to perform on the fasta file data. Available options[CODON-USEAGE, GET-CODONS, TRANSLATE, MUTATE-POINT, MUTATE-RANDOM, GRAPH-GC]')
    args = parser.parse_args()
    analyzer = SequenceAnalyzer(args.input)
    operation = args.operation
    match operation:
       case 'CODON-USAGE':
          print(analyzer.codon_usage())  
       case 'GET-CODONS':
          print('Enter a {Sequence ID}')
          sequence_id = input()
          print(analyzer.get_codons(sequence_id))
       case 'TRANSLATE':
          print(analyzer.dna_to_protein())
       case 'MUTATE-POINT':
          print('Enter a {Sequence ID}, {Position}, and {New Base}')
          user_input = input()
          func_params = user_input.split(" ")
          print(analyzer.point_mutation(func_params[0], int(func_params[1]), func_params[2]))
       case 'MUTATE-RANDOM':
          print('Enter a {Sequence ID}')
          user_input = input()
          func_params = user_input.split(" ")
          print(analyzer.simulate_random_mutation(func_params[0]))
       case 'GRAPH-GC':
          print('Enter a {File Name}')
          file_name = input()
          analyzer.plot_gc_content(file_name)
       case _:
          sys.exit("Invalid operation Please ensure operation matches casing. Terminating")
if __name__ == "__main__":
   run()