.onLoad <- function(libname, pkgname) {
  if (Sys.getenv("R_kmeRtone_DATA") == "")
    Sys.setenv(R_kmeRtone_DATA = "~/kmeRtone_data")
  assign("kmeRtone.data.path", Sys.getenv("R_kmeRtone_DATA"), envir = topenv())
  
  # declare global variables 
  utils::globalVariables(c(
    "chromStart", ".", "kmer_skew", "pos_strand", "rc.seq",
    "proportion", "coor", "kmer_part1", "kmer_part2", 
    "organism_name", "assembly_accession", "asm_name",
    "is_overlap", "i.pos_strand", "i.neg_strand",
    "txStart", "cdsStart", "exonStarts", 
    "exonEnds", "ExonCount", "cdsEnd", "txEnd", "exonCount",
    "chromosome", "element", "V1", "V2", "response", "vcf.file",
    "case", "control", "TUNE", "STUDY_G4_SUSCEPTIBILITY", "G",
    "..col.nums", "INFO", "assembly_level", "genome_rep", "refseq_category", 
    "temp", "species_taxid", "organism_group",
    "cdsStartStat", "cdsEndStat", "name", "chrom",
    "group", "rel_end", "i.start", "i.end", "prob",
    "neg_strand", "strand", "sense", "antisense", "N", 
    "CHROM", "POS", "ALT", "FILTER", "position2", "ref_base",
    "kmer", "len", "z", "Genome.composition", "idx", "Virus.name(s)", "Species",
    "Sequence-Role", "Assigned-Molecule-Location/Type", 
    "RefSeq-Accn", "GenBank-Accn", "susceptible_kmer_count",
    "total_kmer_count", "top_kmer_count", "susceptible k-mer content",
    "highly susceptible k-mer content", "label", "FC_median",
    "Gene Symbol", "selected", "cancer", "susceptible_kmer", "total_kmer", 
    "top_kmer", "transcriptClass", "log2_FC"
  ))
}