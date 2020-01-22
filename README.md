# Multi-purpose and transferrable k-meric enrichment/depletion analysis software.

User inputs:
1. A table of genomic coordinates pointing to damage patterns in csv format.
   - A bedfile format should be accepted as well. Extension end with .bed or .bedfile
   - A genomic coordinate should have these columns: chromosome, start, end, strand
2. Damage pattern
3. Location of control regions
4. Size of kmer
5. Selection of genome
   - If genome of interest is not available from the selection, provide a genome sequence in a single fasta file where each header must contain only chromosome name just like in the genomic coordinate table.

In-store data:
1. Human genome hg19 and hg38