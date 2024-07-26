# KmeRtone: multi-purpose and transferrable k-meric enrichment/depletion analysis software.

![alt text](man/figures/KmeRtone_logo.png)

## Installation

``` r
# The easiest way is to install from CRAN:
install.packages('kmeRtone')

# Otherwise, please download the latest release, then install with
R CMD INSTALL kmeRtone_1.0.tar.gz
```

Alternatively, download and install using the [latest release files from here.](https://github.com/SahakyanLab/kmertone/releases/)

# Overview of `kmeRtone` operations

`KmeRtone` contains many modules. The core module (SCORE) calculates the z-score of k-meric enrichment and depletion. Briefly, the input source are case coordinates for the DNA-related phenomenon under study (e.g. DNA damage, DNA binding, DNA breakage, etc.) and a reference to the chromosome-separated FASTA files. `KmeRtone` calculates the k-mer z-score for every k-mer sequence and generates a table of all k-mer sequences and their associated z-scores. Here, the resulting z-scores indicate how enriched ($z \gg 1$) or depleted ($z \ll 1$) a given k-mer sequence is under the studied phenomenon. 

## `kmeRtone` Input Flags - Overview

Here, we highlight some of the key arguments as input to the `kmeRtone` function. Please refer to the documentation of the function for further details on the required and optional arguments.

1.  **Case coordinate**

    | Flag           | Class                  | Description                                                                                                                                                   |
    |-----------------------------------|-------------------------------------------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------|
    | case.coor.path | `<character>`          | A path to a **folder** containing chromosome-separated genomic coordinates or chromosome-combined BED files. This flag is ignored when case.coor is not NULL. |
    | case      | `<genomic.coordinate>` | A pre-loaded `<genomic.coordinate>` class object..                                                                                          |

2.  **Genome**

    | Flag        | Class         | Description                                                                                                                                 |
    |-------------------------|---------------------------|---------------------------------------------------------------------------------------------------------------------------------------------|
    | genome.name | `<character>` | Available: "hg19" or "hg38". User's own genome name.                                                                                        |
    | genome.path | `<character>` | A path to a user's **folder** containing chromosome-separated fasta files. Default is `NULL`. The file name must be the name of chromosome. |
    | genome      | `<genome>`    | Pre-loaded `<genome>` class object. Default is `NULL`. The two flags above are ignored when this is used.                                   |

3.  **Case characteristics**

    | Flag             | Class         | Description                                      |
    |------------------|---------------|--------------------------------------------------|
    | strand.sensitive | `<bool>`      | Does strand polarity matter?                     |
    | single.case.length      | `<int>`       | Default is `NULL` for unspecified/varied length. |
    | case.pattern     | `<character>` | Default is `NULL` for no pattern.                |

4.  **Case coordinate operation**

    | Flag                   | Class         | Description                                                                                                      |
    |----------------------------------------|-------------------------------|------------------------------------------------------------------------------------------------------------------|
    | rm.case.kmer.overlaps  | `<bool>`      | Default is `TRUE`. This is important to remove neighbouring effect.                                              |
    | merge.replicates       | `<bool>`      | Default is `TRUE`. When merging replicates, duplicated coordinates coming from different replicates are removed. |
    | k                      | `<int>`       | Length of k-mer                                                                                                  |
    | ctrl.rel.pos | `<character>` | Position of control regions relative to the case positions. Input is a vector of length two: `c(from, to)`       |

5.  **Other module flags**

    | Flag       | Class          | Description                                                      |
    |------------|----------------|------------------------------------------------------------------|
    | kmer.table | `<data.table>` | Pre-loaded k-mer table with calculated score. Default is `NULL`. |

6.  **kmeRtone module**

    | Flag   | Class         | Description                                                                                |
    |--------|---------------|--------------------------------------------------------------------------------------------|
    | module | `<character>` | Available module: "score", "tune", "explore", "evolution", "genic element", "cancer", etc. |

7.  **Other**

    | Flag        | Class         | Description                                        |
    |-------------|---------------|----------------------------------------------------|
    | ncpu        | `<int>`       | Number of CPU cores. Default is 1.                 |
    | output.path | `<character>` | A path to an output **folder**. Default is "data/" |

## kmeRtone Input Flags - Additional Description

| Flag           | Description                                                                                                                                                                                                                                                                                                                                                                                                  |
|---------------------------------------------------------------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| single.case.length    | The case length unit is number of nucleotide. In an event where case happens in between two nucleotide e.g. DNA breakage, the case.length is 2 nt.                                                                                                                                                                                                                                                           |
| case.coor.path | Three situations can happen. (1) A folder containing a BED file. A second or more BED files indicates a presence of replicates. (2) A folder containing chromosome-separated files. The file name must be the name of chromosome. (3) A folder containing sub-folders of chromosome-separated files, indicating a presence of replicates. In situation (2) and (3), the coordinates must be a 1-based index due to `R` language conventions. Alternatively, user can specify this with the `case.coor.1st.idx` argument |

## kmeRtone Objects

kmeRtone introduce two class objects: `<genome>` and `<genomic.coordinate>`

1. `<genome>`

kmeRtone comes with two pre-built `<genome>`: hg19 and hg38. The `<genome>`s are saved as uncompressed RDS binary object for fast loading. `print.genome` function is built to print the `<genome>` object. It will show the genome name (e.g. hg19) and genome length by chromosome. The default `base::print` showing the very long sequence will crash the R console.

`<genome>` is an S3-class object with the following contents:

```{bash, include = TRUE}
$seq
 named <character> vector
   chr1              chr2             ...
 c(AACTCGTACC......, ACGTTGGTTC....)  ...

$chr.names
 <character> vector
 c(chr1, chr2, ...) 

$length
 <character> vector
 c(2947924, 2093123, ...)

$name
 <character>
 hg19
```

2. `<genomic.coordinate>`

`<genomic.coordinate>` is an S3-class object. The reason for building this class is to reduce data redundancy in genomic coordinate table (e.g. repeated number of chromosome name and unnecessary column end when case length is fixed). It also helps with organisation of kmeRtone configuration (e.g. k-mer size, case length, etc.) as the `<genomic.coordinate>` object will carry and contain those information. It utilises `<data.table>` to use its inherent feature to update by reference (instead of memory copy) for genomic coordinate table and coordinate status (case vs. k-mer coordinate). This will help to reduce memory (RAM) consumption and keep track what the coordinates refer to (whether the case itself or k-mer). The contents of the `<genomic.coordinate>` object are as follow:

```{bash, include = TRUE}
$chr1
 <data.table>
   start strand ...
 1:   12      +
 2:   16      +
 3:  499      -
 ...

$chr2 .__C__.externalptr

$chr3 ...

$chr... ...

$chr.names
 <character> vector
 c(chr1, chr2, ...)

$status
 <data.table> single row
   is.kmer
 1:   TRUE

$case.length
 <character>
 2

$case.pattern
 <character> vector
 c(CT, TT, ...)
```

## Code Convention

-   Table column name is written in lowercase and snake_case.

-   Function name is written in camelCase. The function filename if it is saved will be the same like the function name except for workflow functions which begin with capital case corresponds to their module letter.

-   Module workflow code begins with a function calling (left-aligned) and ends with variable assignment (right-aligned).

-   Workflow boolean is designed to make it natural to read in English e.g. `if(coor$status$is.kmer)` or `if(coor$is.strand.sensitive)`.

-   Looping uses singular and plural as variable name i.e. `for (chr.name in chr.names)`.

-   The code finish at a standard column number 80 for better viewing.

-   This symbol \<\> refers to R class object e.g. `<character>`

## Quick example

Below is an example code that generates random genomic coordinates and runs the default kmeRtone `SCORE` function to quantify the k-meric enrichment and depletion.

For a detailed explanation, please refer to the `kmeRtone.pdf` in the `vignettes` folder.

```R
library(data.table)
library(kmeRtone)
temp_dir <- tempdir()

#' 1. Randomly generate genomic positions and save results
dir.create("./data", showWarnings = FALSE)

set.seed(1234)
temp_files <- character(22)
for(chr in 1:22){
    genomic_coor <- data.table(
        seqnames = paste0("chr", chr),
        start = sample(
            x = 10000:100000, 
            size = 1000, 
            replace = FALSE
        ),
        width = 2
    )

    f <- file.path(temp_dir, paste0("chr", chr, ".csv"))
    fwrite(genomic_coor, f)
    temp_files[chr] <- f
}

#' 2. Run kmeRtone `score` function
temp_dir_genome <- tempdir()
kmeRtone::kmeRtone(
    case.coor.path = temp_dir, 
    genome.name = "hg19", 
    genome.path = temp_dir_genome,
    strand.sensitive = FALSE, 
    k = 2,
    ctrl.rel.pos = c(80, 500),
    case.pattern = NULL,
    single.case.len = 2,
    output.dir = temp_dir,
    module = "score",
    rm.case.kmer.overlaps = FALSE,
    merge.replicate = TRUE, 
    kmer.table = NULL,
    verbose = TRUE
)
```

The above should generate the below output

```bash
------------------------------------------------------------
                 Extraction of Case K-mers                 
------------------------------------------------------------
Extracting 2-mers from chr1.............DONE! -- 3.23 secs
Extracting 2-mers from chr2.............DONE! -- 3.28 secs
Extracting 2-mers from chr3.............DONE! -- 2.64 secs
Extracting 2-mers from chr4.............DONE! -- 2.56 secs
Extracting 2-mers from chr5.............DONE! -- 2.31 secs
Extracting 2-mers from chr6.............DONE! -- 2.33 secs
Extracting 2-mers from chr7.............DONE! -- 2.04 secs
Extracting 2-mers from chr8.............DONE! -- 1.97 secs
Extracting 2-mers from chr9.............DONE! -- 1.75 secs
Extracting 2-mers from chr10............DONE! -- 1.82 secs
Extracting 2-mers from chr11............DONE! -- 1.75 secs
Extracting 2-mers from chr12............DONE! -- 1.8 secs
Extracting 2-mers from chr13............DONE! -- 1.35 secs
Extracting 2-mers from chr14............DONE! -- 1.25 secs
Extracting 2-mers from chr15............DONE! -- 1.22 secs
Extracting 2-mers from chr16............DONE! -- 1.04 secs
Extracting 2-mers from chr17............DONE! -- 0.97 secs
Extracting 2-mers from chr18............DONE! -- 0.97 secs
Extracting 2-mers from chr19............DONE! -- 0.74 secs
Extracting 2-mers from chr20............DONE! -- 0.71 secs
Extracting 2-mers from chr21............DONE! -- 0.53 secs
Extracting 2-mers from chr22............DONE! -- 0.55 secs

Total time taken: 37.2 secs 
------------------------------------------------------------
                Extraction of Control K-mers                
------------------------------------------------------------
Building control regions of chr1........DONE! -- 3.14 secs
Building control regions of chr2........DONE! -- 3.14 secs
Building control regions of chr3........DONE! -- 2.57 secs
Building control regions of chr4........DONE! -- 2.56 secs
Building control regions of chr5........DONE! -- 2.31 secs
Building control regions of chr6........DONE! -- 2.28 secs
Building control regions of chr7........DONE! -- 2.06 secs
Building control regions of chr8........DONE! -- 1.98 secs
Building control regions of chr9........DONE! -- 1.77 secs
Building control regions of chr10.......DONE! -- 1.78 secs
Building control regions of chr11.......DONE! -- 1.9 secs
Building control regions of chr12.......DONE! -- 1.79 secs
Building control regions of chr13.......DONE! -- 1.36 secs
Building control regions of chr14.......DONE! -- 1.32 secs
Building control regions of chr15.......DONE! -- 1.16 secs
Building control regions of chr16.......DONE! -- 1.06 secs
Building control regions of chr17.......DONE! -- 1.06 secs
Building control regions of chr18.......DONE! -- 0.95 secs
Building control regions of chr19.......DONE! -- 0.71 secs
Building control regions of chr20.......DONE! -- 0.76 secs
Building control regions of chr21.......DONE! -- 0.53 secs
Building control regions of chr22.......DONE! -- 0.63 secs

Total time taken: 36.97 secs 
Extracting 2-mers from chr1.............DONE! -- 3.11 secs
Extracting 2-mers from chr2.............DONE! -- 3.13 secs
Extracting 2-mers from chr3.............DONE! -- 2.62 secs
Extracting 2-mers from chr4.............DONE! -- 2.43 secs
Extracting 2-mers from chr5.............DONE! -- 2.41 secs
Extracting 2-mers from chr6.............DONE! -- 2.22 secs
Extracting 2-mers from chr7.............DONE! -- 2.17 secs
Extracting 2-mers from chr8.............DONE! -- 1.93 secs
Extracting 2-mers from chr9.............DONE! -- 1.86 secs
Extracting 2-mers from chr10............DONE! -- 1.78 secs
Extracting 2-mers from chr11............DONE! -- 1.85 secs
Extracting 2-mers from chr12............DONE! -- 1.76 secs
Extracting 2-mers from chr13............DONE! -- 1.31 secs
Extracting 2-mers from chr14............DONE! -- 1.32 secs
Extracting 2-mers from chr15............DONE! -- 1.18 secs
Extracting 2-mers from chr16............DONE! -- 1.09 secs
Extracting 2-mers from chr17............DONE! -- 1.02 secs
Extracting 2-mers from chr18............DONE! -- 1.07 secs
Extracting 2-mers from chr19............DONE! -- 0.68 secs
Extracting 2-mers from chr20............DONE! -- 0.72 secs
Extracting 2-mers from chr21............DONE! -- 0.53 secs
Extracting 2-mers from chr22............DONE! -- 0.55 secs

Total time taken: 36.97 secs 
------------------------------------------------------------
            Calculation of K-mer Susceptibility            
------------------------------------------------------------
The 2-mer scores are saved at {temp_dir}/score_2-mer.csv

FINISH! Total time taken: 1.85 mins 
```