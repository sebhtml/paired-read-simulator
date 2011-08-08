# Description

This is a simple paired-end read simulator.
Dependencies are boost, a C++ 1998 compiler and the C++ standard library.


# Arguments

VirtualNextGenSequencer GENOME.FASTA SUBSTITUTION_RATE AVERAGE_OUTER_DISTANCE STANDARD_DEVIATION PAIRS READ_LENGTH OUT1.fasta OUT2.fasta


# Example

VirtualNextGenSequencer Ecoli.fasta 0.005 250 25 2000000 75 1.fasta 2.fasta

This generates 2000000 pairs of sequences of length 75 nucleotides. They contain 0.5% substitution errors and are the ends of 
DNA fragments of average length 250 with standard deviation 25. Note that the source strand is random. Finally,
Sequences are written to 1.fasta and 2.fasta.

# Compilation with the GNU C++ compiler

make

# Compilation with the Intel C++ compiler

make CXX=icpc
