# SequenceAnalyzer: A Python Tool for Bioinformatics

SequenceAnalyzer is a robust, object-oriented Python library designed to automate the analysis of DNA sequences. It provides tools for validating genetic data, calculating GC content, analyzing codon usage, and performing DNA-to-protein translation.

## 🌟 Features
- **Data Validation:** Uses Regular Expressions to ensure sequences contain only valid (A, T, C, G) nucleotides.
- **GC Content Analysis:** Calculates and visualizes GC percentage across multiple sequences using Matplotlib.
- **Codon Usage Mapping:** Generates frequency dictionaries for DNA triplets using `collections.Counter`.
- **Protein Translation:** Leverages Biopython's high-performance translation engine to convert DNA to Amino Acid sequences.
- **Type-Safe & Documented:** Fully implemented with Python type hints and descriptive docstrings for professional-grade readability.

## 🧬 Data Sources
The example data provided in the `data/` directory is sourced from the **NCBI (National Center for Biotechnology Information)** Nucleotide database. 

* **NCBI Reference:** [NCBI Nucleotide Database](https://www.ncbi.nlm.nih.gov/nucleotide/)
* **Included Genes:**
  * **TP53:** A critical tumor suppressor gene often referred to as the "guardian of the genome."
  * **BRCA1:** A gene involved in DNA repair; mutations are significantly linked to hereditary breast and ovarian cancer risks.
* **Format:** Data is stored in standard FASTA format, containing raw nucleotide sequences for analysis.

## 🛠️ Installation
To use this tool, you will need Python 3.x and the following libraries:

bash
pip install biopython matplotlib


## 🚀 Quick Start
```python
from sequence_analyzer import SequenceAnalyzer

# Initialize the analyzer with a FASTA file
analyzer = SequenceAnalyzer("data/your_data.fasta")

# Generate a GC content bar chart
analyzer.plot_gc_content()

# Get protein translations as a dictionary
proteins = analyzer.dna_to_protein()

# View codon frequencies
usage = analyzer.codon_usage()
