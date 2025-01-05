#!/usr/bin/env python3
"""
Utility functions for VCF file processing in TentacleSV pipeline.
"""

import os
import sys
import pysam
import logging
from typing import Dict, List, Tuple, Optional, Set
from collections import defaultdict

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

class SVRecord:
    """Class to represent a structural variant record."""
    def __init__(self, chrom: str, pos: int, id: str, ref: str, alt: str,
                 qual: float, filter_status: str, info: Dict[str, str],
                 format_keys: List[str], samples: Dict[str, Dict[str, str]]):
        self.chrom = chrom
        self.pos = pos
        self.id = id
        self.ref = ref
        self.alt = alt
        self.qual = qual
        self.filter = filter_status
        self.info = info
        self.format = format_keys
        self.samples = samples
        
    def get_sv_type(self) -> str:
        """Extract SV type from INFO field or ALT allele."""
        if 'SVTYPE' in self.info:
            return self.info['SVTYPE']
        elif '<' in self.alt and '>' in self.alt:
            return self.alt.split('<')[1].split('>')[0]
        return 'BND'  # Default to breakend if type cannot be determined

def parse_vcf(vcf_path: str) -> List[SVRecord]:
    """
    Parse VCF file and return list of SV records.
    
    Args:
        vcf_path: Path to VCF file
        
    Returns:
        List of SVRecord objects
    """
    records = []
    try:
        with pysam.VariantFile(vcf_path) as vcf:
            for record in vcf:
                # Extract INFO fields
                info = {}
                for key, value in record.info.items():
                    if isinstance(value, tuple):
                        info[key] = ','.join(map(str, value))
                    else:
                        info[key] = str(value)
                
                # Extract FORMAT fields and sample data
                format_keys = list(record.format.keys())
                samples = {}
                for sample in record.samples:
                    samples[sample] = {
                        key: str(record.samples[sample][key])
                        for key in format_keys
                    }
                
                sv_record = SVRecord(
                    record.chrom,
                    record.pos,
                    record.id or '.',
                    record.ref,
                    str(record.alts[0]),
                    record.qual or 0,
                    ','.join(record.filter.keys()) or 'PASS',
                    info,
                    format_keys,
                    samples
                )
                records.append(sv_record)
                
        return records
        
    except Exception as e:
        logger.error(f"Error parsing VCF file {vcf_path}: {str(e)}")
        return []

def analyze_sv_types(records: List[SVRecord]) -> Dict[str, int]:
    """
    Count frequency of each SV type in records.
    
    Args:
        records: List of SVRecord objects
        
    Returns:
        Dictionary mapping SV types to counts
    """
    type_counts = defaultdict(int)
    for record in records:
        sv_type = record.get_sv_type()
        type_counts[sv_type] += 1
    return dict(type_counts)

def filter_sv_records(
    records: List[SVRecord],
    min_size: Optional[int] = None,
    max_size: Optional[int] = None,
    min_qual: Optional[float] = None,
    sv_types: Optional[Set[str]] = None,
    filters: Optional[Set[str]] = None
) -> List[SVRecord]:
    """
    Filter SV records based on various criteria.
    
    Args:
        records: List of SVRecord objects
        min_size: Minimum SV size
        max_size: Maximum SV size
        min_qual: Minimum quality score
        sv_types: Set of allowed SV types
        filters: Set of allowed FILTER values
        
    Returns:
        Filtered list of SVRecord objects
    """
    filtered = []
    for record in records:
        # Check SV size if applicable
        sv_size = int(record.info.get('SVLEN', 0))
        if min_size and sv_size < min_size:
            continue
        if max_size and sv_size > max_size:
            continue
            
        # Check quality score
        if min_qual and record.qual < min_qual:
            continue
            
        # Check SV type
        if sv_types and record.get_sv_type() not in sv_types:
            continue
            
        # Check filter status
        if filters and record.filter not in filters:
            continue
            
        filtered.append(record)
        
    return filtered

def merge_overlapping_svs(
    records: List[SVRecord],
    max_distance: int = 50,
    min_overlap: float = 0.8
) -> List[SVRecord]:
    """
    Merge overlapping SV records based on position and overlap.
    
    Args:
        records: List of SVRecord objects
        max_distance: Maximum distance between breakpoints to consider merging
        min_overlap: Minimum reciprocal overlap required for merging
        
    Returns:
        List of merged SVRecord objects
    """
    # Sort records by chromosome and position
    sorted_records = sorted(records, key=lambda x: (x.chrom, x.pos))
    merged = []
    current = None
    
    for record in sorted_records:
        if current is None:
            current = record
            continue
            
        # Check if records can be merged
        if (current.chrom == record.chrom and
            abs(current.pos - record.pos) <= max_distance and
            current.get_sv_type() == record.get_sv_type()):
            
            # Calculate overlap for non-BND events
            if current.get_sv_type() != 'BND':
                curr_end = current.pos + abs(int(current.info.get('SVLEN', 0)))
                rec_end = record.pos + abs(int(record.info.get('SVLEN', 0)))
                
                overlap_start = max(current.pos, record.pos)
                overlap_end = min(curr_end, rec_end)
                overlap_length = max(0, overlap_end - overlap_start)
                
                curr_length = curr_end - current.pos
                rec_length = rec_end - record.pos
                
                # Check reciprocal overlap
                if (overlap_length / curr_length >= min_overlap and
                    overlap_length / rec_length >= min_overlap):
                    # Merge records (keep the one with higher quality)
                    current = record if record.qual > current.qual else current
                    continue
                    
        merged.append(current)
        current = record
        
    if current is not None:
        merged.append(current)
        
    return merged

if __name__ == "__main__":
    # Example usage
    if len(sys.argv) != 2:
        print("Usage: python vcf_utils.py <vcf_file>")
        sys.exit(1)
        
    vcf_file = sys.argv[1]
    records = parse_vcf(vcf_file)
    
    # Print summary of SV types
    type_counts = analyze_sv_types(records)
    print("\nSV Type Counts:")
    for sv_type, count in type_counts.items():
        print(f"{sv_type}: {count}")
        
    # Example filtering
    filtered_records = filter_sv_records(
        records,
        min_size=50,
        min_qual=20,
        sv_types={'DEL', 'INS', 'DUP', 'INV'}
    )
    print(f"\nFiltered records: {len(filtered_records)}/{len(records)}")
