# ContigChecker
A tool for checking circularity and correctness of contigs of especially metagenomes

## Motivation

Recently, by virtue of the development of long-read sequencing technology, it has become possible to assemble complete DNA sequences of bacterial chromosomes, plasmids and phages. In terms of a sequence of nucleotides, a DNA sequence has two possible shapes: **linear or circular**. In genome assembly, we additionally have to consider about the completeness of contigs: **partial or complete**. That is to say, in microbial genome assembly we would like to classify each assembled contig into one of the following three categories:

1. (complete) circular
   * In general, we cannot conclude a contig as "partial circular" (unless we already know the true circular seuquence).
2. complete linear
3. partial linear
   * In the era of short reads, only this type of contigs was usually obtained.
   * (Maybe we should split this category into "left-partial linear" and "right-partial linear"?)

However, current genome assemblers usually do not care about the categories above and sometimes output contigs with wrong shape; especially, when there exists a pair of terminal repeats in a linear contig, usual assembly strategies would recognize it as circular because the 5'-end and 3'-end of the contigs have overlap. Another minor point is the distinguishment between complete linear and partial linear, which current assemblers have not reported as well so far.

Our main focus here is to detect and revise misclassification of complete/partial linear as circular just like above, but we also try to identify other types of misassembly. Although for now we do not go into assembly strategy itself, we care about the guarantee on the generality of our method regardless of assembly strategy.

## What triggers misassembly?

Misassembly, including misclassification of a linear contig as circular, is caused by false overlaps (no false overlaps lead no misassembly) stemming from

1. repeats, or
2. non-repeats.

In Case 1, it is expected that we have no (or a few) reads spanning the repeat which wrongly connects its upstream unique sequence and downstream unique sequence, because the original (meta)genome actually does not include such a region corresponding to the spanning reads. In other words, for the unambiguity of a contig, the maximal repeat of every position of the contig must be spanned by (sufficient number of) reads. The problem for this case is to detect repeats in contigs *de novo*.

Note that, in our PacBio metagenome assembly, we used FALCON's unitigs instead of contigs, which implies that the assembled sequences are expected to be stopped at any ambiguous sites in the graph and basically do not contain non-spanned repeats. However, it does not hold for other assemblies.

Case 2 is due to sequencing artifacts or noise. In this case, supporting reads around the misassembly points should be very few. Conversely speaking, most reads end at the breakpoint. Using orthogonal datasets (e.g. mapping Illumina reads to PacBio contigs) would be effective as well for detecting technology-specific errors.

In both cases, we perform analysis basically similar to DAMASKER and DASCRUBBER. They leverage read pileup data for repeat annotation and adapter/chimer/low-quality region detection, respectively.



### Case 1: How to detect repeats?

​    To determine repeat-spanning reads, we have to annotate repeat regions in the contigs. If every repeat properly appears multiple times in the assembled contigs, we can simply avoid repeats by using only uniquely mapped reads. However, in practice, assemblers can output only a single, collapsed repeat instance, and then reads can be uniquely mapped even on repeats.

​    Therefore, we have to find repeat regions only by the mapped reads. DAZZLER's REPmask module determines repeat regions from raw reads by simply setting a single, user-specified threshold for depth; however, we cannot directly use it for metagenome due to the heterogeneity of abundance. Moreover, just looking at read depth would sometimes fail to detect repeats. By the nature of a repeat, a significant number of (clipped) reads should start from the 5'-end of the repeat and end at the 3'-end with sufficient coverage. We call these as start/end spikes.

​    Accordingly, we handle each contig independently, and also use start/end spikes by clipped reads for initial screening of the candidates of repeat intervals. In order to achieve this, we must admit multiple mapping (instead of random assignment to one of the multiple mappings) and clipped alignment (or non-proper overlap in assembly terminology) during the mapping step. Otherwise, we would miss increase of read depth and also start/end spikes.



### Case 2: How to detect potential assembly breakpoints (chimera)?

​    Assembly breakpoints would show rapid drop of read depth and also an end spike by clipped reads just before the breakpoint and a start spike just after the breakpoint. Both a repeat and a breakpoint should have a pair of a start spike and an end spike. Most of the remaining "orphan spikes" would be just noise, but we should check whether the read depth around them is anormal or not.



## Characteristics of complete linear

​    The difference between complete linear and partial linear is start/end spikes by only properly mapped (i.e. non-clipped) reads at the terminals of a complete linear contig.

​    The story gets a little complicated when repeats are involved in. For the above example of crAssphage, if a contig $c$ has a "collapsed repeat" $c[s,t]$ of a pair of terminal repeats, we would not have spanning reads, but at the same time observe a start spike at $s$ and an end spike at $t$. Then we should output a modified sequence $c[s,t]+c[t,|c|]+c[0,s]+c[s,t]$ as complete linear. Internal repeats (local loop having incoming edge(s) and outgoing edge(s) in an assembly graph) would let us observe neither spanning reads nor start/end spikes.



## Considering the cut point

​    So far we have basically ignored about the cut point of a contig (when assuming as circular). Checking reads bridging the cut point requires a special care, and moreover, reads are less likely to be mapped around the cut point. To avoid this, we can simply prepare two sequences for each contig $c$; one is $c$ itself, and the other one is $c[|c|/2,|c|]+c[0,|c|/2]$ (namely, concatenating two sequences generated by cutting $c$ at the furthest point from the original cut point of $c$) if the concatenated sequence is neither redundant nor gapped around the original cut point.

​    However, contig sequences around the cut point depend on the assembler used (some assemblers output redundant sequences? FALCON seems to output non-redundant contigs though; but anyway it depends on the assembler or assembly graph used). If the sequences around the cut point are not appropriate, there exist no spanning reads around the seuquences and also no start/end spikes, thus the contig would be wrongly categorized as partial linear. Therefore, it is worth looking at reads mapped to both (original) terminals of the contig and modifying the terminal sequences before we concatenate the two half sequences.



## In summary,

|    Category     |     Repeat-spanning reads      |    Spike by properly mapped reads    |
| :-------------: | :----------------------------: | :----------------------------------: |
|    Circular     | Everywhere including cut point |                 None                 |
| Complete linear |    Except terminal repeats     | Start [end] spike at 5'-end [3'-end] |
| Partial linear  |   (Except terminal repeats)    |                 None                 |

Therefore, important data we have to determine are essentially:

- Regional annotation of repeats
- Positional start/end frequencies of reads
  - the number of reads starting from [ending at] each position

In most steps, start/end spikes are used for a good and robust indicator of repeats, breakpoints or complete linear, because normally multiple reads rarely start/end at the same position although read depth is frequently fluctuated.

|      Type       | Spike by properly mapped reads |     Spike by clipped reads      |
| :-------------: | :----------------------------: | :-----------------------------: |
|     Repeat      |               -                | start + end with depth increase |
|     Chimera     |               -                |   end + start with depth drop   |
| Complete linear |          start + end           |                -                |



## Prerequirements

* DAZZLER modules
  * DAZZ_DB
  * daligner
  * damapper
* python libraries
  * logzero

## How to install

```bash
$ python3 setup.py install
```

## Input data

* `<contigs_fasta>`
  * A fasta file of the contigs you assembled using any assembler.
* `<reads_fasta>`
  * A fasta file of the reads you used for the assembly.

## How to run

```bash
$ contig_checker.py map <contigs_fasta> <reads_fasta>
$ contig_checker.py plot
$ contig_checker.py categorize
$ contig_checker.py modify
```