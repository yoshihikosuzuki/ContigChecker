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


## Requirements

* DAZZLER modules
  * DAZZ_DB
  * daligner
  * damapper
* python libraries
  * logzero


## How to install and run

```
$ git clone https://github.com/yoshihikosuzuki/ContigChecker
$ cd ContigChecker
$ python3 setup.py install
```

After that, you need to copy `input.cfg` to somewhere and edit it, and run as follows.

```
$ ctcc.py input.cfg
```


## Config file format and input data

* `<contigs_fasta>`
  * A fasta file of the contigs you assembled using any assembler.
* `<reads_fasta>`
  * A fasta file of the reads you used for the assembly.
