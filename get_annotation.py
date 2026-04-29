#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
BLAST Hit Annotator (Coordinate-based)

Description:
    This script retrieves functional annotations (gene products) for CRISPR 
    spacer alignments. Instead of downloading complete reference genomes, 
    it efficiently uses the NCBI Entrez API to fetch only the exact genomic 
    region defined by the BLAST hit coordinates. By analyzing this specific 
    window, it identifies if the spacer targets a Coding Sequence (CDS) 
    and extracts its product, or classifies it as an intergenic region.

Author: Andrea Fernández
Date: April 2026

Usage:
    python Anotaciones_v2.py
"""

import pandas as pd
from Bio import Entrez, SeqIO
import time
from tqdm import tqdm
import warnings

# --- CONFIGURATION ---
Entrez.email = "*"
# Entrez.api_key = "*" 
blast_input_file = "*.tsv"
output_filename = "*.tsv"
# ---------------------

# A dictionary to cache API requests. If multiple spacers hit the exact same 
# region (or if the input has duplicate rows), we prevent redundant NCBI queries.
api_cache = {}

def get_gene_by_coordinates(accession, start, end, max_retries=3):
    """
    Queries NCBI for the specific genomic window where the spacer aligned.
    This approach is highly efficient and avoids downloading entire genomes.
    """
    # BLAST coordinates can be inverted if the hit is on the reverse strand
    seq_start = min(start, end)
    seq_stop = max(start, end)
    
    cache_key = (accession, seq_start, seq_stop)
    if cache_key in api_cache:
        return api_cache[cache_key]
        
    for attempt in range(max_retries):
        try:
            # Respect NCBI rate limits (change to 0.11 if using an API key)
            time.sleep(0.35) 
            
            # CORE LOGIC: seq_start and seq_stop limit the download to the hit area
            handle = Entrez.efetch(
                db="nucleotide", 
                id=accession, 
                rettype="gb", 
                retmode="text",
                seq_start=seq_start,
                seq_stop=seq_stop
            )
            
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                record = SeqIO.read(handle, "genbank")
            handle.close()

            # Since the downloaded region is restricted to the spacer's target area, 
            # any CDS feature present automatically overlaps the spacer.
            hits_products = []
            for feature in record.features:
                if feature.type == "CDS":
                    # Extract the product name or mark as unknown
                    product = feature.qualifiers.get('product', ['unknown protein'])[0]
                    hits_products.append(product)

            # Assign category based on findings
            if not hits_products:
                result = "Intergenic"
            else:
                # Remove duplicates and join if the hit falls between overlapping ORFs
                result = "; ".join(list(set(hits_products)))
                
            api_cache[cache_key] = result
            return result
            
        except Exception as e:
            if attempt == max_retries - 1:
                return "Error/Not Found"
            # Wait a bit before retrying if connection fails
            time.sleep(1)

def process_annotations():
    """
    Main function to parse the BLAST output and annotate each alignment.
    """
    # Define column names based on the custom BLAST output format
    blast_cols = ["qseqid", "sseqid", "pident", "length", "mismatch",
                  "gapopen", "evalue", "bitscore", "staxids", "stitle",
                  "sstart", "send", "sstrand"]

    print(f"Loading BLAST file: {blast_input_file} ...")
    df_blast = pd.read_csv(blast_input_file, sep="\t", names=blast_cols, low_memory=False)

    print(f"Annotating {len(df_blast)} alignments by fetching exact coordinates...")
    
    # Apply the annotation function row by row with a progress bar
    tqdm.pandas(desc="Mapping ORFs")
    df_blast['Gene_Product'] = df_blast.progress_apply(
        lambda row: get_gene_by_coordinates(row['sseqid'], row['sstart'], row['send']), 
        axis=1
    )

    # Save to file
    df_blast.to_csv(output_filename, sep="\t", index=False)
    print(f"\nProcess complete! Annotated file saved as '{output_filename}'.")

if __name__ == "__main__":
    process_annotations()