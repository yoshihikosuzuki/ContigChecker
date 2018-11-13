# ContigChecker
A tool for checking circularity and correctness of contigs of especially metagenomes

## Quick summary

This tool offers a function to classify each given contig into one of the following three categories:

|    Category     |     Repeat-spanning reads      |    Spike by properly mapped reads    |
| :-------------: | :----------------------------: | :----------------------------------: |
|    Circular     | Everywhere including cut point |                 None                 |
| Complete linear |    Except terminal repeats     | Start [end] spike at 5'-end [3'-end] |
| Partial linear  |   (Except terminal repeats)    |                 None                 |

In addition, following types of regions in the given contigs are annotated and (optionally) modified:

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

```
$ python3 setup.py install
```


## Input data

* `<contigs_fasta>`
  * A fasta file of the contigs you assembled using any assembler.
* `<reads_fasta>`
  * A fasta file of the reads you used for the assembly.


## How to run

```
$ contig_checker.py map <contigs_fasta> <reads_fasta>
$ contig_checker.py plot
$ contig_checker.py categorize
$ contig_checker.py modify
```