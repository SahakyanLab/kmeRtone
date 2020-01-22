# Multi-purpose and transferrable k-meric enrichment/depletion analysis software.

User inputs:
1. A table of genomic coordinates pointing to damage patterns in csv format.
   - A bedfile format should be accepted as well. Extension end with .bed or .bedfile
   - A genomic coordinate should have these columns: chromosome, start, end, strand
```
dt <- fread(dt)
print(dt)
         chromosome start   end strand
      1:       chrI     1     1      +
      2:       chrI     2     2      +
      3:       chrI     4     4      +
      4:       chrI     5     5      +
      5:       chrI     6     6      +
     ---                                                                           
 682620:      chrmt 80012 80012      -
```
2. Damage pattern: a string e.g. "G"
3. Mode of directionality of strands: "sensitive" or "insensitive"
3. Location of control regions: a vector e.g c(80, 500)
4. Size of kmer
5. Selection of genome
   - If genome of interest is not available from the selection, provide a genome sequence in a single fasta file where each header must contain only chromosome name just like in the genomic coordinate table.

In-store data:
1. Human genome hg19 and hg38
2. Data for analysis e.g. location of transcription start site, G4, etc.

## Kmertone Workflow
1. Error checking: make sure the genomic coordinates point to the damage pattern.