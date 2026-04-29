# collector_plasmids_pipelines

Bioinformatic Analysis of CRISPR Spacer Acquisition

This repository contains a suite of custom scripts developed for the analysis of CRISPR spacer acquisition. The pipeline covers the entire workflow from processing raw sequencing data to host identification, evolutionary trajectory reconstruction, and functional annotation of acquired spacers.

1. Spacer Extraction and Pre-processing
spacer_extraction.py
  This script processes raw FASTQ files from CRISPR array sequencing to build a non-redundant library of spacers.

    Array Identification: Scans reads in both orientations to identify CRISPR repeats.
    Spacer Isolation: Extracts sequences between repeats. It is configured to include the last nucleotide of the upstream repeat (a PAM-derived nucleotide) to         extend the spacer by 1 bp, enhancing alignment resolution.
    Filtering: Automatically filters out sequences exceeding 38 bp.
    Output: Generates deduplicated FASTA/FASTQ files and detailed tables including global metrics (counts, lengths) and array-specific labeling.

2. Spacer Assignment in Defined Pairs
For experiments involving known bacterial strains, the following scripts ensure high-precision mapping:

    cross_mapping.py: Parsers BAM files using the pysam library to inspect the XA tag (alternative hits). It identifies "ambiguous" spacers that align to the          primary target but also have secondary alignments to other genomes with comparable scores, filtering out potential sequence homology artifacts.

    extract_PAMs.py: Cross-references alignment coordinates with genomic FASTA files to retrieve the 3 nucleotides immediately upstream of the protospacer on the      target strand.

3. Host Identification and Taxonomy
  LCA_algorithm.py
    When analyzing acquisitions from complex communities (e.g., wastewater), this script resolves multi-species hits from BLASTn searches.

    LCA Logic: It evaluates hits sharing the highest bitscore and assigns the taxonomy to the Lowest Common Ancestor (LCA).

4. Trajectory Determination
  chord_diagram_trajectories.R
    This R script reconstructs the movement trajectories of the Collector Plasmid between different bacterial hosts.

    Chronological Logic: Since CRISPR systems add new spacers at the leader end, the script reverses the physical sequence order to reflect the chronological          acquisition (e.g., Spacer 2 was acquired before Spacer 1).
    Reconstruction: It evaluates taxonomic transitions between consecutive valid spacers at the genus level.
    Visualization: Aggregated trajectories are represented via a Chord Diagram.

5. Functional Classification of Targets
  get_annotation.py
      Uses the NCBI Entrez API to query the genomic window defined by BLAST hit coordinates. It identifies overlapping Coding Sequences (CDS) and extracts product       names. Regions without a CDS are classified as intergenic.
  categories.py
Takes the output from the annotation script and classifies entries into 29 functional categories using a custom keyword dictionary.
  Priority System: Specific terms take priority over broader terms.
  Majority-Vote: For spacers with multiple annotation hits, a majority-vote rule is applied to determine the most representative category.

Requirements
Bioinformatics Tools: BWA-MEM (v0.7.18), Samtools (v1.21), BLASTn (v2.16.0).
Python Dependencies: pysam, Biopython (for Entrez API).

R Dependencies: circlize (or similar for Chord Diagrams)
