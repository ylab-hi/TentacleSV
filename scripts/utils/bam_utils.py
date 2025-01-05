#!/usr/bin/env python3
"""
Utility functions for BAM file processing in TentacleSV pipeline.
"""

import os
import sys
import pysam
import logging
from typing import Dict, List, Tuple, Optional

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def check_bam_file(bam_path: str) -> Tuple[bool, Optional[str]]:
    """
    Check if BAM file exists and is properly formatted/indexed.
    
    Args:
        bam_path: Path to BAM file
        
    Returns:
        Tuple of (is_valid, error_message)
    """
    try:
        if not os.path.exists(bam_path):
            return False, f"BAM file not found: {bam_path}"
            
        # Check if BAM is indexed
        index_path = bam_path + ".bai"
        if not os.path.exists(index_path):
            return False, f"BAM index not found: {index_path}"
            
        # Try opening BAM file
        with pysam.AlignmentFile(bam_path, "rb") as bam:
            # Check if BAM is sorted
            if not bam.header.get('HD', {}).get('SO') == 'coordinate':
                return False, f"BAM file not coordinate sorted: {bam_path}"
                
        return True, None
        
    except Exception as e:
        return False, f"Error checking BAM file: {str(e)}"

def get_read_groups(bam_path: str) -> List[Dict[str, str]]:
    """
    Extract read group information from BAM file.
    
    Args:
        bam_path: Path to BAM file
        
    Returns:
        List of read group dictionaries
    """
    try:
        with pysam.AlignmentFile(bam_path, "rb") as bam:
            return bam.header.get('RG', [])
    except Exception as e:
        logger.error(f"Error getting read groups: {str(e)}")
        return []

def get_sequencing_type(bam_path: str) -> str:
    """
    Determine sequencing type (short/long) from BAM read groups and read lengths.
    
    Args:
        bam_path: Path to BAM file
        
    Returns:
        String indicating sequencing type ('short', 'long_ont', or 'long_pacbio')
    """
    try:
        with pysam.AlignmentFile(bam_path, "rb") as bam:
            # Check read groups first
            rg_info = bam.header.get('RG', [])
            for rg in rg_info:
                platform = rg.get('PL', '').lower()
                if platform == 'illumina':
                    return 'short'
                elif platform == 'ont':
                    return 'long_ont'
                elif platform in ['pacbio', 'pb']:
                    return 'long_pacbio'
            
            # If no platform info in read groups, check read lengths
            total_reads = 0
            long_reads = 0
            for read in bam.fetch():
                total_reads += 1
                if read.query_length > 1000:  # Typical threshold for long reads
                    long_reads += 1
                if total_reads >= 1000:  # Sample first 1000 reads
                    break
            
            if long_reads / total_reads > 0.5:
                # Default to ONT if we can't determine specific long-read technology
                return 'long_ont'
            return 'short'
            
    except Exception as e:
        logger.error(f"Error determining sequencing type: {str(e)}")
        return 'short'  # Default to short reads if we can't determine

def calculate_coverage(bam_path: str, reference_lengths: Dict[str, int]) -> float:
    """
    Calculate average coverage across the genome.
    
    Args:
        bam_path: Path to BAM file
        reference_lengths: Dictionary of reference sequence names and lengths
        
    Returns:
        Float indicating average coverage
    """
    try:
        total_bases = 0
        total_length = sum(reference_lengths.values())
        
        with pysam.AlignmentFile(bam_path, "rb") as bam:
            for read in bam.fetch():
                if not read.is_secondary and not read.is_supplementary:
                    total_bases += read.query_length
                    
        return total_bases / total_length
        
    except Exception as e:
        logger.error(f"Error calculating coverage: {str(e)}")
        return 0.0

def filter_supplementary_alignments(
    input_bam: str,
    output_bam: str,
    min_mapping_quality: int = 20
) -> bool:
    """
    Filter BAM file to remove supplementary alignments and low mapping quality reads.
    
    Args:
        input_bam: Path to input BAM file
        output_bam: Path to output filtered BAM file
        min_mapping_quality: Minimum mapping quality threshold
        
    Returns:
        Boolean indicating success
    """
    try:
        with pysam.AlignmentFile(input_bam, "rb") as in_bam:
            with pysam.AlignmentFile(output_bam, "wb", template=in_bam) as out_bam:
                for read in in_bam.fetch():
                    if (not read.is_supplementary and 
                        not read.is_secondary and 
                        read.mapping_quality >= min_mapping_quality):
                        out_bam.write(read)
                        
        # Index output BAM
        pysam.index(output_bam)
        return True
        
    except Exception as e:
        logger.error(f"Error filtering BAM file: {str(e)}")
        return False

if __name__ == "__main__":
    # Example usage
    if len(sys.argv) != 2:
        print("Usage: python bam_utils.py <bam_file>")
        sys.exit(1)
        
    bam_file = sys.argv[1]
    is_valid, error = check_bam_file(bam_file)
    
    if not is_valid:
        print(f"BAM file validation failed: {error}")
        sys.exit(1)
        
    read_groups = get_read_groups(bam_file)
    print(f"Read groups: {read_groups}")
    
    seq_type = get_sequencing_type(bam_file)
    print(f"Sequencing type: {seq_type}")
