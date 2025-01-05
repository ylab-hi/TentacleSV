#!/usr/bin/env python3
"""
Quality control utilities for read mapping in TentacleSV pipeline.
"""

import os
import sys
import pysam
import logging
import numpy as np
from typing import Dict, List, Tuple, Optional
from collections import defaultdict
import matplotlib.pyplot as plt

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

class MappingQC:
    """Class for performing mapping quality control."""
    
    def __init__(self, bam_path: str):
        """
        Initialize MappingQC object.
        
        Args:
            bam_path: Path to BAM file
        """
        self.bam_path = bam_path
        self.stats = {}
        self._calculate_basic_stats()
        
    def _calculate_basic_stats(self):
        """Calculate basic mapping statistics."""
        try:
            with pysam.AlignmentFile(self.bam_path, "rb") as bam:
                # Initialize counters
                total_reads = 0
                mapped_reads = 0
                proper_pairs = 0
                supplementary = 0
                duplicate = 0
                mapq_stats = []
                insert_sizes = []
                read_lengths = []
                
                # Collect statistics
                for read in bam.fetch():
                    total_reads += 1
                    
                    if read.is_mapped:
                        mapped_reads += 1
                        mapq_stats.append(read.mapping_quality)
                        read_lengths.append(len(read.query_sequence))
                        
                    if read.is_proper_pair:
                        proper_pairs += 1
                        if read.template_length > 0:  # Only count positive template lengths
                            insert_sizes.append(abs(read.template_length))
                            
                    if read.is_supplementary:
                        supplementary += 1
                    if read.is_duplicate:
                        duplicate += 1
                
                # Calculate percentages and statistics
                self.stats = {
                    'total_reads': total_reads,
                    'mapped_reads': mapped_reads,
                    'mapping_rate': mapped_reads / total_reads if total_reads > 0 else 0,
                    'proper_pairs': proper_pairs,
                    'proper_pair_rate': proper_pairs / (total_reads/2) if total_reads > 0 else 0,
                    'supplementary_alignments': supplementary,
                    'duplicate_rate': duplicate / total_reads if total_reads > 0 else 0,
                    'mean_mapq': np.mean(mapq_stats) if mapq_stats else 0,
                    'median_mapq': np.median(mapq_stats) if mapq_stats else 0,
                    'mean_read_length': np.mean(read_lengths) if read_lengths else 0,
                    'median_read_length': np.median(read_lengths) if read_lengths else 0,
                    'mean_insert_size': np.mean(insert_sizes) if insert_sizes else 0,
                    'median_insert_size': np.median(insert_sizes) if insert_sizes else 0,
                    'insert_size_std': np.std(insert_sizes) if len(insert_sizes) > 1 else 0
                }
                
        except Exception as e:
            logger.error(f"Error calculating mapping statistics: {str(e)}")
            self.stats = {}
            
    def get_coverage_stats(self, min_mapq: int = 20) -> Dict[str, float]:
        """
        Calculate coverage statistics.
        
        Args:
            min_mapq: Minimum mapping quality to include in coverage calculation
            
        Returns:
            Dictionary of coverage statistics
        """
        try:
            with pysam.AlignmentFile(self.bam_path, "rb") as bam:
                coverage_stats = {}
                for ref in bam.references:
                    # Calculate coverage array
                    coverage = bam.count_coverage(
                        ref,
                        quality_threshold=min_mapq
                    )
                    
                    # Convert to numpy array and calculate statistics
                    cov_array = np.sum(coverage, axis=0)
                    coverage_stats[ref] = {
                        'mean': np.mean(cov_array),
                        'median': np.median(cov_array),
                        'std': np.std(cov_array),
                        'zero_coverage_bases': np.sum(cov_array == 0),
                        'bases_above_10x': np.sum(cov_array >= 10) / len(cov_array)
                    }
                    
                return coverage_stats
                
        except Exception as e:
            logger.error(f"Error calculating coverage statistics: {str(e)}")
            return {}
            
    def plot_insert_size_distribution(self, output_path: str):
        """
        Generate insert size distribution plot.
        
        Args:
            output_path: Path to save the plot
        """
        try:
            insert_sizes = []
            with pysam.AlignmentFile(self.bam_path, "rb") as bam:
                for read in bam.fetch():
                    if read.is_proper_pair and read.template_length > 0:
                        insert_sizes.append(abs(read.template_length))
            
            if insert_sizes:
                plt.figure(figsize=(10, 6))
                plt.hist(insert_sizes, bins=50, alpha=0.75)
                plt.title('Insert Size Distribution')
                plt.xlabel('Insert Size (bp)')
                plt.ylabel('Count')
                plt.axvline(np.median(insert_sizes), color='r', linestyle='dashed',
                          label=f'Median: {np.median(insert_sizes):.0f}')
                plt.legend()
                plt.savefig(output_path)
                plt.close()
                
        except Exception as e:
            logger.error(f"Error plotting insert size distribution: {str(e)}")
            
    def plot_mapping_quality_distribution(self, output_path: str):
        """
        Generate mapping quality distribution plot.
        
        Args:
            output_path: Path to save the plot
        """
        try:
            mapq_values = []
            with pysam.AlignmentFile(self.bam_path, "rb") as bam:
                for read in bam.fetch():
                    if read.is_mapped:
                        mapq_values.
