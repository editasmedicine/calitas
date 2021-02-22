# CALITAS

This repository is home to **CALITAS**, a _CRISPR/Cas-aware ALigner for In silico off-TArget Search_.  CALITAS implements a customized gapped alignment of guide sequences to genomes and other reference sequences, returning consistent and non-redundant alignments.

## Overview

CALITAS is a suite of bioinformatic tools for enumerating candidate off-target sites for CRISPR guide sequences and standardizing their alignment.  It features tools for searching an entire genome for candidate off-target sites, as well as for aligning guides to specific sequences or locations in a genome.

Key features of CALITAS include:

1. Detection of all candidate off-target sites up to the requested number of mismatches and gaps
2. Elimination of redundant alignments resulting in a single best or canonical alignment per locus
3. Searches with multiple PAM sequences and/or PAM-less searching
4. Integration of known variants in VCF format into genome-wide off-target searches
5. Customizable scoring system for weighting mismatches vs. gaps and differences in the protospacer vs. PAM

CALITAS is not intended to _predict_ active off-target sites, but rather to _enumerate_ candidate off-target sites for further investigation.

## Getting CALITAS

### Releases

Binary releases of CALITAS are available from the [releases](releases) page on GitHub.  The downloadable JAR files contain CALITAS and all dependencies, and need only a [Java RunTime version 8 or higher](https://openjdk.java.net/install/) installed.  Once downloaded CALITAS can be run as follows to produce usage information:

```
java -Xmx8g -jar calitas.jar
```

### Bioconda

Relases of CALITAS are also available via [Bioconda](https://bioconda.github.io/).  To use these releases you will first need to install the conda package manager - [Miniconda](https://docs.conda.io/en/latest/miniconda.html) is recommended if you do not already use conda.  CALITAS can then be installed with:

```
conda install -c bioconda calitas
```

The conda release comes with a small helper script and can be run simply as `calitas`.


### Building from Source

Both release and development versions of the code can be built from source.  Builds are performed with [sbt](http://www.scala-sbt.org/).  Once sbt is installed run:

```
sbt clean assembly
```

This will build a JAR equivalent to a release JAR at `calitas/target/scala-2.12/calitas.jar`

## Publication

A CALITAS manuscript is currently in review.  This section will be updated with a reference to the publication when it is available.

## Usage

CALITAS has four available sub-commands:

* **AlignToReference** performs glocal alignment of query sequence to a window on the reference
* **PairwiseAlignSequences** performs pairwise alignment of sequences
* **PrepareVcf** prepares a VCF for optimal use by SearchReference
* **SearchReference** searches a FASTA file for alignments of a guide+PAM(s)

It should be noted that commands which use a genome or reference FASTA file require that the FASTA file have both an index and a sequence dictionary.  These can be generated using [samtools](https://github.com/samtools/samtools) as follows:

```bash
samtools faidx ref.fa
samtools dict -a <assembly-name> -s <species> -o ref.dict ref.fa
```

The following is an example of invoking `SearchReference` to find candidate off target sequences in the HG38 genome for a single guide and PAM (note that location and sequence of the PAM is indicated by providing guide sequence in upper case and PAM sequence in lower case):

```
calitas SearchReference \
  -i CTTGCCCCACAGGGCAGTAAnrg \
  -I myguide \
  -r hg38.fa \
  -o myguide.hits.txt \
  --max-guide-diffs 5 \
  --max-pam-mismatches 1 \
  --max-gaps-between-guide-and-pam 3
```

The last three parameters are optional and replicate the defaults.  Additional options are available; detailed usage including all available parameters can be obtained by running `calitas SearchReference`.

The following is an example of running `AlignToReference` to produce standarized alignments at locations where guide(s) are known to align.  The invocation will produce the single best alignment per query sequence and target location:

```
calitas AlignToReference -i input.txt -r hg38.fa -o output.txt --window-size 60
```

With the following being an example of the tab-delimited input file for `AlignToReference`:

```
id	query	chrom	position
1	CTTGCCCCACAGGGCAGTAAnrg	chr1	13358
2	CTTGCCCCACAGGGCAGTAAnrg	chr1	510578
3	CTTGCCCCACAGGGCAGTAAnrg	chr1	844033
```

The output of both `SearchReference` and `AlignToReference` is a tab-delimited text file with one row per candidate off-target site including the following columns:

|column name|description|
|:---|:---|
|`guide_id`|Name/ID of guide.|
|`unpadded_guide_sequence`|The sequence of the guide used, unpadded.|
|`genome_build`|The assembly name of the searched genome (e.g. `HG38`).|
|`chromosome`|Chromosome for target sequence alignment (eg: `chr3`).|
|`coordinate_start`|Start of the unpadded target sequence in the genome, 0-based open ended, excluding PAM.|
|`coordinate_end`|End of the unpadded target sequence in the genome, 0-based open ended, excluding PAM.|
|`strand`|Either `+` or `-`. The reported strand is the strand of the genome which matches the guide sequence.  E.g. if strand is reported as `+` this means the guide resembles the sequence on the top strand of the genome, and will bind to the bottom strand of the genome.|
|`unpadded_target_sequence`|The unpadded target sequence (as DNA) as found in the genome, without gaps/bulges, excluding PAM. Reported sequence matches the reported `strand` (i.e. `-` strand hits will report the reverse complement of the genomic sequence). |
|`ten_bases_5_prime`|The 10 bases from the reference genome immediately 5' of the off-target location (`coordinate_start`/`coordinate_end`), respecting `strand`.|
|`ten_bases_3_prime`|The 10 bases from the reference genome immediately 3' of the off-target location (`coordinate_start`/`coordinate_end`), respecting `strand`.|
|`pam_used`|PAM used in the alignment (eg: `nrg`).|
|`variant_id`|When searching using a VCF, a semi-colon separated list of variant IDs (e.g. `rs1234;rs2345`) that have non-reference alleles present in the off-target alignment. May be empty when no variants are present.|
|`variant_description`|When searching using a VCF, a semi-colon separated list of variant descriptions in the format `id:pos:ref>alt:af` where pos is the position within the alignment and af is the allele frequency of the alternate allele.|
|`variant_vcf`|When searching using a VCF, a string composed of filname of the VCF followed by a colon (`:`) and then the MD5 of the VCF.|
|`allele_frequency`|When searching using a VCF, the minimum allele frequency of any variant included in the target alignment.|
|`score`|Alignment score (including PAM).|
|`guide_mm`|Mismatches in the guide region (excluding PAM).|
|`guide_gaps`|Gaps in the guide region (excluding PAM).|
|`guide_mm_plus_gaps`|Total gaps and mismatches in the guide region (excluding PAM).|
|`pam_mm`|Mismatches in the PAM region.|
|`total_mm_plus_gaps`|Total count of mismatches and gaps across the both guide and PAM regions.|
|`padded_guide`|Guide + PAM sequence with padding for mismatches and bulges/gaps.|
|`padded_alignment`|Visual representation of the guide-target alignment: `|` for matches, `.` for mismatches, `~` for gaps.|
|`padded_target`|Target sequence, including PAM, with padding for mismatches and bulges/gaps.|
|`padded_extra_8_bases_5_prime`|an additional 8 bases on the 5' side of the padded_target.|
|`padded_extra_8_bases_3_prime`|an additional 8 bases on the 3' side of the padded_target.|
|`cigar`|Cigar representation of guide sequence alignment.|
|`unpadded_guide_sequence_length`|Length of guide sequence not including PAM, nor gaps/bulges.|
|`unpadded_target_sequence_length`|Length of target sequence not including PAM, nor gaps/bulges.|
|`aligner`|Aligner name, e.g. `CALITAS:SearchReference`.|
|`aligner_version`|Version number of the CALITAS software used to produce the alignments.|
|`aligner_search_pam`|Comma-separated list of pams used during the search.|
|`aligner_other_parameters`|A semicolon-separated list of parameters provided when running CALITAS.|
|`time_stamp`|A date and time stamp for when the alignment run was started (UTC) in this format: Wed Jan 6 16:58:29 UTC 2021.|

### Note on coordinates produced by CALITAS

All genome coordinates in CALITAS output files are 0-based open ended, the same systemm used in BED files.  I.e. the first base of a chromosome or sequence is represented by 0, and when describing a region on chromosome we specify the first base _included_ in the interval as the start and the base _after_ the interval as the end.  E.g. the first 10bp on a chromosome would be represented as `0-10`.
