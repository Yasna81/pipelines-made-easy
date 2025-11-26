import sys
from collections import Counter
import pysam

def main():
    if len(sys.argv) < 2:
        sys.exit("usage: python aligncheck.py in.bam")
    
    fname = sys.argv[1]

    try:
        bamfile = pysam.AlignmentFile(fname, "rb")
    except (IOError, ValueError) as e:
        sys.exit(f"Error while opening file: {e}")
    
    stats = Counter()
    mapq_scores = []  # Store mapping qualities for distribution
    read_lengths = []  # Store read lengths
    insert_sizes = []  # insert size for paired end 

    for read in bamfile:
        stats["total"] += 1
        stats["qcfail"] += int(read.is_qcfail)
        stats["duplicate"] += int(read.is_duplicate)
        stats["secondary"] += int(read.is_secondary)
        stats["supplementary"] += int(read.is_supplementary)  

        # paired
        stats["paired"] += int(read.is_paired)
        stats["read1"] += int(read.is_read1)
        stats["read2"] += int(read.is_read2)

        read_len = read.query_length
        if read_len > 0:
            read_lengths.append(read_len)
        
        if read.is_unmapped:
            stats["unmapped"] += 1
            continue
        
        # only for mapped reads
        stats["mapped"] += 1
        stats["proper pair"] += int(read.is_proper_pair)
        
        # mapping quality of reads 
        mapq = read.mapping_quality
        mapq_scores.append(mapq)  

        stats['mapq=0'] += int(mapq == 0)
        stats['mapq<=10'] += int(mapq <= 10)
        stats['mapq<=30'] += int(mapq <= 30)
        stats['mapq>30'] += int(mapq > 30)

        # strand info
        stats['forward strand'] += int(not read.is_reverse)
        stats['reverse strand'] += int(read.is_reverse)

        # Insert size for paired-end reads
        if read.is_paired and read.is_proper_pair and read.template_length > 0 and abs(read.template_length) < 10000:
            insert_sizes.append(abs(read.template_length))
        
        # CIGAR operations analysis
        if read.cigartuples:
            for operation, length in read.cigartuples:
                if operation == 0:  # M
                    stats['cigar_M'] += length
                elif operation == 1:  # I
                    stats['cigar_I'] += length
                elif operation == 2:  # D
                    stats['cigar_D'] += length
                elif operation == 4:  # S
                    stats['cigar_S'] += length
                elif operation == 5:  # H
                    stats['cigar_H'] += length
    
    bamfile.close()
    
    # Calculate additional statistics
    if stats['mapped'] > 0:
        stats['mapping_rate'] = 100 * stats['mapped'] / stats['total']
    
    if mapq_scores:
        stats['avg_mapq'] = sum(mapq_scores) / len(mapq_scores)
        stats['median_mapq'] = sorted(mapq_scores)[len(mapq_scores)//2]
    
    if read_lengths:
        stats['avg_read_length'] = sum(read_lengths) / len(read_lengths)
        stats['median_read_length'] = sorted(read_lengths)[len(read_lengths)//2]
    
    if insert_sizes:
        stats['avg_insert_size'] = sum(insert_sizes) / len(insert_sizes)
        stats['median_insert_size'] = sorted(insert_sizes)[len(insert_sizes)//2]
    
    # Output all statistics
    output_order = (
        "total", "mapped", "unmapped", "mapping_rate",
        "paired", "read1", "read2", "proper pair",
        "qcfail", "duplicate", "secondary", "supplementary",
        "mapq=0", "mapq<=10", "mapq<=30", "mapq>30", "avg_mapq", "median_mapq",
        "forward strand", "reverse strand",
        "avg_read_length", "median_read_length",
        "avg_insert_size", "median_insert_size",
        "cigar_M", "cigar_I", "cigar_D", "cigar_S", "cigar_H"
    )

    print(f"Detailed alignment statistics for {fname}:")
    print("=" * 50)
    
    for key in output_order:
        if key in stats:
            value = stats[key]
            if key in ['mapping_rate', 'avg_mapq', 'median_mapq', 
                      'avg_read_length', 'median_read_length',
                      'avg_insert_size', 'median_insert_size']:
                print(f"{key}: {value:.2f}")
            else:
                percentage = 100 * value / stats["total"] if stats["total"] > 0 else 0
                print(f"{key}: {value} ({percentage:.2f}%)")

if __name__ == "__main__":
    main()