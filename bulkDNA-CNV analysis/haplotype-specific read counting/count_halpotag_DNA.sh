#!/bin/bash

# Check if correct number of arguments is provided
if [[ $# -ne 3 ]]; then
    echo "Usage: $0 <BAM_FILE> <BINSIZE> <OUTPUT_FILE>"
    exit 1
fi

# Define input BAM file, bin size, and output file
BAM_FILE=$1
BINSIZE=$2
OUTPUT_FILE=$3

# Check if BAM file exists
if [[ ! -f "$BAM_FILE" ]]; then
    echo "Error: BAM file '$BAM_FILE' not found!"
    exit 1
fi

# Ensure binsize is a valid integer
if ! [[ "$BINSIZE" =~ ^[0-9]+$ ]]; then
    echo "Error: BINSIZE must be a positive integer!"
    exit 1
fi

# Process BAM file and cache bin contributions
samtools view -F 3852 -q 20 "$BAM_FILE" | awk -v binsize="$BINSIZE" -v output_file="$OUTPUT_FILE" '
BEGIN {
    OFS="\t";
    prev_chrom = "";
    min_bin = -1;
}

{
    chrom = $3;                 # Chromosome
    start = $4;                 # Start position
    cigar = $6;                 # CIGAR string
    has_hp = 0;                 # Flag for HP tag
    has_ps = 0;                 # Flag for PS tag
    ps_tag = "NA";              # Default PS tag value
    hp_tag = "NA";              # Default HP tag value

    # Scan auxiliary fields (from column 11 onwards) for HP and PS tags
    for (i=11; i<=NF; i++) {
        if ($i ~ /^HP:i:/) { has_hp = 1; hp_tag = substr($i, 6); }
        if ($i ~ /^PS:i:/) { has_ps = 1; ps_tag = substr($i, 6); }
    }

    # Process only if both HP and PS tags are present
    if (has_hp && has_ps) {
        match_length = 0;
        while (match(cigar, /[0-9]+M/)) {
            match_length += substr(cigar, RSTART, RLENGTH-1);
            cigar = substr(cigar, RSTART + RLENGTH);
        }

        end = start + match_length;  # Compute read end position

        # Compute bin indices for start and end positions
        bin_start = int(start / binsize);
        bin_end = int(end / binsize);

        # If chromosome changes or new reads are at least 2 bins ahead, flush the cache
        if (chrom != prev_chrom || bin_start >= min_bin + 2) {
            for (key in cache) {
                split(key, fields, SUBSEP);
                print fields[1], fields[2], fields[3], fields[4], cache[key] >> output_file;
                delete cache[key];
            }
            min_bin = bin_start;
            prev_chrom = chrom;
        }

        # Case 1: Read is entirely within a single bin
        if (bin_start == bin_end) {
            cache[chrom, ps_tag, hp_tag, bin_start] += 1;
        }
        # Case 2: Read spans exactly two adjacent bins
        else if (bin_end == bin_start + 1) {
            frac1 = ((bin_start + 1) * binsize - start) / match_length;
            frac2 = (end - bin_end * binsize) / match_length;
            cache[chrom, ps_tag, hp_tag, bin_start] += frac1;
            cache[chrom, ps_tag, hp_tag, bin_end] += frac2;
        }
        # Case 3: Read spans more than two bins → Ignore
    }
}

END {
    # Flush remaining cache after processing all reads
    for (key in cache) {
        split(key, fields, SUBSEP);
        print fields[1], fields[2], fields[3], fields[4], cache[key] >> output_file;
    }
}'

