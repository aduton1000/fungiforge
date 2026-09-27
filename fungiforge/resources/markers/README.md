# Secondary-locus reference proteins (W2.3)

`marker_reference_proteins.faa` — query proteins used by `bin/extract_markers.py` (tblastn) to
locate the secondary identification loci in an assembly:

| Marker | Protein | Source |
|---|---|---|
| CaM | calmodulin (cmdA) | *A. fumigatus* Af293, NCBI XP_751821.2 |
| BenA | β-tubulin (benA) | *A. fumigatus* Af293, NCBI XP_752456.1 |
| TEF1 | translation elongation factor 1-α (tef1) | *A. fumigatus* Af293, NCBI XP_750388.2 |
| RPB2 | RNA polymerase II second largest subunit | *S. cerevisiae*, UniProt P08518 |

These only locate the loci; the extracted nucleotide regions are then compared with the
type-material reference sets staged by `bin/fetch_references.sh markers` (NCBI "sequence from
type" records per locus) for the species-level call. LSU (D1/D2) comes from barrnap's 28S
coordinates, ITS from ITSx as before.
