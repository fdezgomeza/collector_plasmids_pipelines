#!/usr/bin/env python3
import pysam
from pyfaidx import Fasta
import argparse
import csv
import os

def reverse_complement(dna_seq):
    """Computes the reverse complement of a DNA sequence."""
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    return "".join(complement.get(base, base) for base in reversed(dna_seq.upper()))

def get_pam_from_alignment(alignment, reference_fasta, pam_length=3):
    """
    Extracts the PAM sequence for a given pysam alignment object.
    The PAM is defined as `pam_length` nucleotides 5' to the spacer sequence
    on the strand the spacer aligns to.

    Args:
        alignment (pysam.AlignedSegment): The alignment object.
        reference_fasta (pyfaidx.Fasta): The reference genome object.
        pam_length (int): Length of the PAM sequence.

    Returns:
        str: The PAM sequence, or 'N/A' if it cannot be determined.
    """
    if alignment.is_unmapped:
        return "N/A"

    chrom = alignment.reference_name
    ref_start = alignment.reference_start  # 0-based start of alignment on reference
    ref_end = alignment.reference_end      # 0-based, exclusive end of alignment on reference

    pam_seq_on_fwd_strand = ""

    try:
        contig = reference_fasta[chrom]
    except KeyError:
        # print(f"Warning: Contig {chrom} not found in reference FASTA.")
        return "N/A (Contig_Not_Found)"


    if not alignment.is_reverse:
        # Spacer aligns to the FORWARD strand of the reference
        # PAM is 5' to the spacer, so it's upstream on the forward strand
        # Ref: ---PAM--->[Spacer_Match]--->
        #      ^pam_s   ^pam_e=ref_start
        pam_start_coord = ref_start - pam_length
        pam_end_coord = ref_start

        if pam_start_coord < 0:
            # print(f"Warning: PAM extends before start of contig {chrom} for {alignment.query_name}")
            return "N/A (Boundary)"
        
        pam_seq_on_fwd_strand = contig[pam_start_coord:pam_end_coord].seq.upper()
        return pam_seq_on_fwd_strand # PAM is already on the correct strand

    else:
        # Spacer aligns to the REVERSE strand of the reference
        # The spacer sequence itself would be the reverse complement of what's on the fwd strand.
        # PAM is 5' to the spacer *on the reverse strand*.
        # This means it's 3' to the aligned region on the FORWARD strand.
        # Ref Fwd: <---[Spacer_Match_RevComp]<---PAM---
        #                     ^ref_end=pam_s ^pam_e
        pam_start_coord = ref_end
        pam_end_coord = ref_end + pam_length

        if pam_end_coord > len(contig):
            # print(f"Warning: PAM extends after end of contig {chrom} for {alignment.query_name}")
            return "N/A (Boundary)"
            
        pam_seq_on_fwd_strand = contig[pam_start_coord:pam_end_coord].seq.upper()
        # Since the spacer is on the reverse strand, the PAM we want is also on the reverse strand.
        # So, we take the reverse complement of the sequence fetched from the forward reference.
        return reverse_complement(pam_seq_on_fwd_strand)

    return "N/A (Error)" # Should not reach here

def main():
    parser = argparse.ArgumentParser(description="Extract PAM sequences for aligned spacers from a BAM file.")
    parser.add_argument("-b", "--bam", required=True, help="Input BAM file (sorted and indexed).")
    parser.add_argument("-f", "--fasta", required=True, help="Reference FASTA file (must be indexed with .fai).")
    parser.add_argument("-o", "--output", required=True, help="Output CSV file for PAMs.")
    parser.add_argument("--pam_length", type=int, default=3, help="Length of the PAM sequence (default: 3).")

    args = parser.parse_args()

    if not os.path.exists(args.bam + ".bai") and not os.path.exists(args.bam + ".csi"):
         print(f"Error: BAM index (.bai or .csi) not found for {args.bam}. Please index the BAM file.")
         return
    if not os.path.exists(args.fasta + ".fai"):
        print(f"Error: FASTA index (.fai) not found for {args.fasta}. Please index using 'samtools faidx'.")
        return

    try:
        reference_genomes = Fasta(args.fasta, sequence_always_upper=True)
        bamfile = pysam.AlignmentFile(args.bam, "rb")
    except Exception as e:
        print(f"Error opening files: {e}")
        return

    print(f"Processing BAM file: {args.bam}")
    print(f"Using reference FASTA: {args.fasta}")
    print(f"Outputting to: {args.output}")

    results = []
    processed_alignments = 0
    pams_found = 0

    for alignment in bamfile.fetch(): # Iterate over all mapped reads
        if alignment.is_unmapped:
            continue
        
        processed_alignments += 1
        if processed_alignments % 10000 == 0:
            print(f"  Processed {processed_alignments} alignments...")

        spacer_id = alignment.query_name
        pam = get_pam_from_alignment(alignment, reference_genomes, args.pam_length)

        if pam not in ["N/A", "N/A (Boundary)", "N/A (Contig_Not_Found)", "N/A (Error)"]:
             pams_found +=1

        results.append({
            "Spacer_ID": spacer_id,
            "Reference_Contig": alignment.reference_name,
            "Alignment_Start_on_Ref": alignment.reference_start, # 0-based
            "Alignment_End_on_Ref": alignment.reference_end,     # 0-based, exclusive
            "Strand": "-" if alignment.is_reverse else "+",
            "Spacer_Sequence_Aligned": alignment.query_sequence, # Sequence from BAM (may be RC if is_reverse)
            "CIGAR": alignment.cigarstring,
            "PAM_Length": args.pam_length,
            "PAM_Sequence": pam
        })

    bamfile.close()
    reference_genomes.close()

    # Write results to CSV
    if results:
        with open(args.output, 'w', newline='') as csvfile:
            fieldnames = results[0].keys()
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(results)
        print(f"\nSuccessfully processed {processed_alignments} alignments.")
        print(f"Found {pams_found} valid PAM sequences.")
        print(f"Results saved to {args.output}")
    else:
        print("No alignments processed or no results to write.")

if __name__ == "__main__":
    main()