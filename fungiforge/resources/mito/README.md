# Mitochondrial reference sets (W2.7)

`mito_core_proteins.faa` — the 15 core mitochondrial protein-coding genes (cox1-3, cob, atp6/8/9,
nad1-6, nad4L, rps3) from three mitogenomes spanning the fungal kingdom: *Aspergillus fumigatus*
Af293 (NC_017016), *Saccharomyces cerevisiae* S288c (NC_001224; no nad genes), *Cryptococcus
neoformans* JEC21 (NC_004336). Headers `gene|protein_id|species`. Used by `bin/mito_extract.py`
(tblastn: which assembly contigs are mitochondrial) and `bin/mito_annotate.py` (gene set, copies,
introns from HSP gaps).

`afum_mito_rrna.fna` — the Af293 rnl (large subunit) and rns (small subunit) mitochondrial rRNAs
for blastn location of the rRNA genes.

Intron-encoded ORFs (LAGLIDADG/GIY-YIG endonucleases, maturases) are deliberately excluded: they
are mobile elements, counted by stage 10.
