# cyp51A promoter fixture (real sequence)

`c87_published_cyp51A.fa` / `.gbk` — a 3.6-kb slice of the published *Aspergillus fumigatus* C87
assembly (ENA GCA_949125185, contig CASBLV010000001, Hemmings, Rhodes & Fisher 2023), 2 kb
upstream of the cyp51A ATG to 1.6 kb into the gene, with a minimal `gene` feature. The promoter
carries the genuine TR34 tandem duplication; `bin/cyp51a_TR.py` must report `TR34`, 2 copies.

Note (2026-09-15): the reads deposited for this isolate (ONT ERR10820709, Illumina ERR9791656)
do NOT carry the duplication (0/84 ONT and 1/76 Illumina reads spanning the site show an
insertion), so those accessions cannot serve as a read-level positive control; this slice of the
published assembly is the assembly-level truth used instead.
