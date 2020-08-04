# prepare kmers for tardigrade

library(Biostrings)
library(data.table)

fasta.path = "data/genome/tardigrade/YOKOZUNA-1.scaffolds.fa"
tardigrade = readDNAStringSet(fasta.path)

for (i in 1:10) {
  kmers <- extractGenomeKmers(tardigrade, k=i)
  fwrite(kmers, paste0("data/kmers/tardigrade/kmers_k-", i, ".csv"))
}


UVizs <- fread("../kmertone/")




# save kmers
fwrite(data.table(kmer = names(kmers), count = kmers), "data/tardigrade/kmers_8.csv")

# load data
kmers <- fread("data/tardigrade/kmers.csv")

# includes only central DNA pattern
kmers[substring(kmer, pattern.pos[1], pattern.pos[length(pattern.pos)]) %in% DNA.pattern]


# apply UViz score