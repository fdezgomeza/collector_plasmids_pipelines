#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Cross-Mapping Alignment Analyzer

Description:
    This script identifies and analyzes "cross-mapping" events where CRISPR 
    spacers align to more than one of the provided reference genomes. It 
    evaluates primary alignments and searches for alternative alignments 
    (via the BWA 'XA' tag) to determine if a single spacer sequence matches 
    multiple biological sources.

Author: Andrea Fernández
Date: April 2026

Usage:
    python cross_mapping.py <input_file.sam or input_file.bam>
"""

import pysam
import sys
import os
import csv
from collections import Counter

# --- REFERENCE CONFIGURATION ---
# Mapping keywords in reference names to biological sources
SOURCE_MAPPING = {
    "E. coli": ["Escherichia", "coli"],
    "Klebsiella": ["Klebsiella", "oxytoca", "pneumoniae", "variicola", "AQCG"], 
    "pMTX1001": ["pMtx1001"],
    "pMTX1002": ["pMtx1002"]
}

def identify_source(reference_name):
    """
    Determines the biological source based on the reference sequence name
    using the SOURCE_MAPPING dictionary.
    """
    if not reference_name:
        return "Unknown"
        
    for source, keywords in SOURCE_MAPPING.items():
        for key in keywords:
            if key.lower() in reference_name.lower():
                return source
    return "Unknown"

def analyze_multimapping(input_file):
    """
    Parses a SAM/BAM file to detect reads that align to multiple sources.
    Outputs a CSV file with detailed cross-mapping data.
    """
    if not os.path.exists(input_file):
        print(f"Error: The file '{input_file}' was not found.")
        return

    base_name = os.path.basename(input_file)
    file_root, _ = os.path.splitext(base_name)
    output_csv = f"multi_{file_root}.csv"

    print(f"--- Alignment Analysis Initiated ---")
    print(f"Input File: {input_file}")
    
    total_aligned_reads = 0
    multimapped_reads_count = 0
    cross_mapping_types = Counter()

    # Determine file mode (BAM binary vs SAM text)
    mode = "rb" if input_file.endswith(".bam") else "r"

    try:
        with pysam.AlignmentFile(input_file, mode) as samfile, \
             open(output_csv, mode='w', newline='') as csvfile:

            csv_writer = csv.writer(csvfile)
            csv_writer.writerow([
                "Read_ID", 
                "Cross_Mapping_Type", 
                "Primary_Source", 
                "Alternative_Sources", 
                "Sequence_Length",
                "Sequence_Content"
            ])

            # Iterate through alignments
            for read in samfile.fetch(until_eof=True):
                # FILTERING:
                # 1. Ignore unmapped reads
                if read.is_unmapped:
                    continue
                
                # 2. Ignore Secondary (flag 256) and Supplementary (flag 2048) alignments.
                # This prevents double-counting biological spacers. 
                # Alternative alignment information is extracted from the 'XA' tag.
                if read.is_secondary or read.is_supplementary:
                    continue

                total_aligned_reads += 1
                detected_sources = set()

                # A. Analyze Primary Alignment
                primary_ref = read.reference_name
                primary_source = identify_source(primary_ref)
                if primary_source != "Unknown":
                    detected_sources.add(primary_source)

                # B. Analyze Alternative Alignments (XA Tag)
                alt_sources_list = []
                if read.has_tag("XA"):
                    xa_tag = read.get_tag("XA")
                    alternatives = xa_tag.split(';')
                    
                    for alt in alternatives:
                        if not alt: continue
                        ref_alt_name = alt.split(',')[0]
                        source = identify_source(ref_alt_name)
                        if source != "Unknown":
                            detected_sources.add(source)
                            alt_sources_list.append(source)

                # C. Detect and Record Cross-Mapping
                if len(detected_sources) > 1:
                    multimapped_reads_count += 1
                    mapping_combination = " + ".join(sorted(detected_sources))
                    cross_mapping_types[mapping_combination] += 1
                    
                    csv_writer.writerow([
                        read.query_name,
                        mapping_combination,
                        primary_source,
                        ", ".join(set(alt_sources_list)),
                        len(read.query_sequence),
                        read.query_sequence
                    ])

        # Final Report Generation
        print("--- Analysis Complete ---")
        print(f"Total Unique Spacers Aligned: {total_aligned_reads}")
        
        if multimapped_reads_count > 0:
            percentage = (multimapped_reads_count / total_aligned_reads) * 100
            print(f"Cross-mapping Detected: {multimapped_reads_count} reads ({percentage:.2f}%)")
            
            print("\nDistribution of Cross-Mapping Types:")
            print(f"{'Source Combination':<40} | {'Count':<10}")
            print("-" * 55)
            for combo, count in cross_mapping_types.most_common():
                print(f"{combo:<40} | {count:<10}")
            
            print(f"\nDetailed data saved to: {output_csv}")
        else:
            print("No cross-species mapping detected.")

    except Exception as e:
        print(f"An error occurred during processing: {e}")

if __name__ == "__main__":
    if len(sys.argv) > 1:
        analyze_multimapping(sys.argv[1])
    else:
        print("Usage: python cross_mapping.py <input_file.sam or input_file.bam>")