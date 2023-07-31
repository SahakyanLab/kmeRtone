# Multi-purpose and transferrable k-meric enrichment/depletion analysis software.

## Kmertone Operation

Kmertone contains many modules. The core module (SCORE) is calculating z-score of k-meric enrichment and depletion. Draw the flowchart.

## Kmertone Input Flags - Overview

1.  Case coordinate

    | Flag           | Class                  | Description                                                                                                                                                   |
    |----------------|------------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------|
    | case.coor.path | `<character>`          | A path to a **folder** containing chromosome-separated genomic coordinates or chromosome-combined BED files. This flag is ignored when case.coor is not NULL. |
    | case.coor      | `<genomic.coordinate>` | A pre-loaded `<genomic.coordinate>` class object..                                                                                          |

2.  Genome

    | Flag        | Class         | Description                                                                                                                                 |
    |-------------|---------------|---------------------------------------------------------------------------------------------------------------------------------------------|
    | genome.name | `<character>` | Available: "hg19" or "hg38". User's own genome name.                                                                                        |
    | genome.path | `<character>` | A path to a user's **folder** containing chromosome-separated fasta files. Default is `NULL`. The file name must be the name of chromosome. |
    | genome      | `<genome>`    | Pre-loaded `<genome>` class object. Default is `NULL`. The two flags above are ignored when this is used.                                   |

3.  Case characteristics

    | Flag             | Class         | Description                                      |
    |------------------|---------------|--------------------------------------------------|
    | strand.sensitive | `<bool>`      | Does strand polarity matter?                     |
    | case.length      | `<int>`       | Default is `NULL` for unspecified/varied length. |
    | case.pattern     | `<character>` | Default is `NULL` for no pattern.                |

4.  Case coordinate operation

    | Flag                   | Class         | Description                                                                                                      |
    |------------------------|---------------|------------------------------------------------------------------------------------------------------------------|
    | rm.case.kmer.overlaps  | `<bool>`      | Default is `TRUE`. This is important to remove neighbouring effect.                                              |
    | merge.replicates       | `<bool>`      | Default is `TRUE`. When merging replicates, duplicated coordinates coming from different replicates are removed. |
    | rm.diff.pattern        | `<bool>`      | Default is `FALSE`. This removes false positive coordinate (different DNA sequence) as a result of sequencing.   |
    | k                      | `<int>`       | Length of k-mer                                                                                                  |
    | ctrl.rel.pos | `<character>` | Position of control regions relative to the case positions. Input is a vector of length two: `c(from, to)`       |

5.  Other module flags

    | Flag       | Class          | Description                                                      |
    |------------|----------------|------------------------------------------------------------------|
    | kmer.table | `<data.table>` | Pre-loaded k-mer table with calculated score. Default is `NULL`. |

6.  Kmertone module

    | Flag   | Class         | Description                                                                                |
    |--------|---------------|--------------------------------------------------------------------------------------------|
    | module | `<character>` | Available module: "score", "tune", "explore", "evolution", "genic element", "cancer", etc. |

7.  Other

    | Flag        | Class         | Description                                        |
    |-------------|---------------|----------------------------------------------------|
    | ncpu        | `<int>`       | Number of CPU cores. Default is 1.                 |
    | output.path | `<character>` | A path to an output **folder**. Default is "data/" |

## Kmertone Input Flags - Additional Description

| Flag           | Description                                                                                                                                                                                                                                                                                                                                                                                                  |
|----------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| case.length    | The case length unit is number of nucleotide. In an event where case happens in between two nucleotide e.g. DNA breakage, the case.length is 2 nt.                                                                                                                                                                                                                                                           |
| case.coor.path | Three situations can happen. (1) A folder containing a BED file. A second or more BED files indicates a presence of replicates. (2) A folder containing chromosome-separated files. The file name must be the name of chromosome. (3) A folder containing sub-folders of chromosome-separated files, indicating a presence of replicates. In situation (2) and (3), the coordinates must be a 1-based index. |

## Kmertone Objects

Kmertone introduce two class objects: `<genome>` and `<genomic.coordinate>`

### `<genome>`

Kmertone comes with two pre-built `<genome>`: hg19 and hg38. The `<genome>`s are saved as uncompressed RDS binary object for fast loading. `print.genome` function is built to print the `<genome>` object. It will show the genome name (e.g. hg19) and genome length by chromosome. The default `base::print` showing the very long sequence will crash the R console.

`<genome>` is an S3-class object with the following contents:

``` r
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

### `<genomic.coordinate>`

`<genomic.coordinate>` is an S3-class object. The reason for building this class is to reduce data redundancy in genomic coordinate table (e.g. repeated number of chromosome name and unnecessary column end when case length is fixed). It also helps with organisation of kmertone configuration (e.g. k-mer size, case length, etc.) as the `<genomic.coordinate>` object will carry and contain those information. It utilises `<data.table>` to use its inherent feature to update by reference (instead of memory copy) for genomic coordinate table and coordinate status (case vs. k-mer coordinate). This will help to reduce memory (RAM) consumption and keep track what the coordinates refer to (whether the case itself or k-mer). The contents of the `<genomic.coordinate>` object are as follow:

``` r
$chr1
 <data.table>
   start strand ...
 1:   12      +
 2:   16      +
 3:  499      -
 ...

$chr2 ...

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
