#!/bin/bash

# Check if correct number of arguments is provided
if [[ $# -ne 2 ]]; then
    echo "Usage: $0 <input_file> <output_file>"
    exit 1
fi

input_file=$1
output_file=$2

awk '
BEGIN {
    # Initialize arrays
    delete hp1_counts
    delete hp2_counts
    delete positions
}
{
    chr = $1
    pos = $2
    hp = $3
    read_count = $5

    # Create a unique key for each position
    pos_key = chr "\t" pos

    # Sum read counts by haplotype
    if (hp == 1) {
        hp1_counts[pos_key] += read_count
    } else if (hp == 2) {
        hp2_counts[pos_key] += read_count
    }

    # Track all positions weve seen
    positions[pos_key] = 1
}
END {
    # Print header
    print "chr", "PS", "HP1counts", "HP2counts"

    # Output counts for each position
    for (p in positions) {
        # Default to 0 if no counts exist for a haplotype
        hp1 = (p in hp1_counts) ? hp1_counts[p] : 0
        hp2 = (p in hp2_counts) ? hp2_counts[p] : 0
        print p, hp1, hp2
    }
}' OFS="\t" "$input_file" | sort -k1,1 -k2,2n > "$output_file"

echo "Processing complete. Output saved to $output_file"
