# draft_genome_simulator_with_visualization_red_circle_fixed.py
import argparse
import random
import os
import sys
import tempfile
import shutil
from pathlib import Path
from typing import List, Tuple, Dict
from Bio import SeqIO, SeqRecord
from Bio.Seq import Seq
import numpy as np
import subprocess
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt
from pycirclize import Circos

# Force unbuffered output
sys.stdout = sys.__stdout__
sys.stderr = sys.__stderr__

# Debug print function with immediate flush
def debug_print(message):
    print(f"{message}", flush=True)
    sys.stdout.flush()

# Utility to compute GC% in windows
def compute_gc_windows(seq: str, window: int) -> List[Tuple[int, int, float]]:
    gc_windows = []
    for i in range(0, len(seq) - window + 1, window):
        fragment = seq[i:i + window]
        gc = (fragment.count("G") + fragment.count("C")) / len(fragment)
        gc_windows.append((i, i + window, gc))
    return gc_windows

# Fragment large contigs for manageable Red processing
def fragment_large_contigs(records: List[SeqRecord.SeqRecord], max_contig_size: int = 1_000_000, 
                          fragment_size: int = 500_000, overlap: int = 50_000) -> Tuple[List[SeqRecord.SeqRecord], List[Tuple[int, int, int]]]:
    """
    Fragment contigs larger than max_contig_size into smaller overlapping fragments.
    
    Returns:
        fragmented_records: List of SeqRecord objects (original small contigs + fragments of large contigs)
        fragment_mapping: List of (original_start, original_end, fragment_index) for coordinate mapping
    """
    fragmented_records = []
    fragment_mapping = []
    current_pos = 0
    
    for record in records:
        if len(record.seq) <= max_contig_size:
            # Keep small contigs as-is
            fragmented_records.append(record)
            fragment_mapping.append((current_pos, current_pos + len(record.seq), len(fragmented_records) - 1))
            current_pos += len(record.seq)
        else:
            # Fragment large contigs
            print(f"Fragmenting large contig {record.id}: {len(record.seq):,} bp -> {fragment_size:,} bp fragments with {overlap:,} bp overlap")
            
            seq_str = str(record.seq)
            contig_start = current_pos
            fragment_idx = 0
            
            for start in range(0, len(seq_str), fragment_size - overlap):
                end = min(start + fragment_size, len(seq_str))
                
                # Create fragment record
                fragment_seq = seq_str[start:end]
                fragment_record = SeqRecord.SeqRecord(
                    Seq(fragment_seq),
                    id=f"{record.id}_frag_{fragment_idx}",
                    description=f"Fragment {fragment_idx} of {record.id} ({start}-{end})"
                )
                fragmented_records.append(fragment_record)
                
                # Map fragment coordinates to original genome coordinates
                fragment_mapping.append((contig_start + start, contig_start + end, len(fragmented_records) - 1))
                
                fragment_idx += 1
                
                # Break if we've reached the end
                if end >= len(seq_str):
                    break
            
            current_pos += len(record.seq)
    
    print(f"Fragmentation complete: {len(records)} original contigs -> {len(fragmented_records)} fragments")
    return fragmented_records, fragment_mapping

# Map fragment coordinates back to original genome coordinates
def map_fragment_coords_to_original(repeat_regions: List[Tuple[int, int]], 
                                   fragment_mapping: List[Tuple[int, int, int]],
                                   fragmented_records: List[SeqRecord.SeqRecord]) -> List[Tuple[int, int]]:
    """
    Map coordinates from fragmented sequences back to original genome coordinates.
    """
    original_repeats = []
    
    for frag_start, frag_end in repeat_regions:
        # Find which fragment this coordinate belongs to
        cumulative_pos = 0
        for i, record in enumerate(fragmented_records):
            if cumulative_pos <= frag_start < cumulative_pos + len(record.seq):
                # Found the fragment
                # Get the original coordinates from fragment_mapping
                for orig_start, orig_end, frag_idx in fragment_mapping:
                    if frag_idx == i:
                        # Map fragment coordinates to original coordinates
                        relative_start = frag_start - cumulative_pos
                        relative_end = frag_end - cumulative_pos
                        
                        # Ensure we don't exceed the original fragment boundaries
                        mapped_start = orig_start + relative_start
                        mapped_end = min(orig_start + relative_end, orig_end)
                        
                        if mapped_start < mapped_end:  # Valid mapping
                            original_repeats.append((mapped_start, mapped_end))
                        break
                break
            cumulative_pos += len(record.seq)
    
    return original_repeats

# FIXED: Parse Red output file and map coordinates properly
def parse_red_output_with_mapping(red_output_file: Path, contig_mapping: Dict[str, int]) -> List[Tuple[int, int]]:
    """
    Parse Red repeat output file and map contig-relative coordinates to concatenated genome coordinates.
    Red BED format: chrName\tstart\tend (tab-separated, 0-based coordinates)
    
    Args:
        red_output_file: Path to Red BED output file
        contig_mapping: Dictionary mapping contig names to their start positions in concatenated genome
    
    Returns:
        List of repeat regions with coordinates mapped to concatenated genome
    """
    repeat_regions = []
    
    if not red_output_file.exists():
        print(f"Warning: Red output file {red_output_file} does not exist")
        return repeat_regions
    
    debug_print(f"Parsing Red output: {red_output_file}")
    debug_print(f"Available contig mappings: {list(contig_mapping.keys())[:5]}...")  # Show first 5
    
    try:
        with open(red_output_file, 'r') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                
                # Red BED format: chrName\tstart\tend
                parts = line.split('\t')
                if len(parts) >= 3:
                    try:
                        contig_full_name = parts[0]
                        # Extract just the contig ID (part before first space) since Red outputs full headers
                        contig_name = contig_full_name.split()[0]
                        contig_start = int(parts[1])  # 0-based coordinates from Red
                        contig_end = int(parts[2])    # 0-based coordinates from Red
                        
                        # Map contig-relative coordinates to concatenated genome coordinates
                        if contig_name in contig_mapping:
                            # Get the start position of this contig in the concatenated genome
                            genome_offset = contig_mapping[contig_name]
                            
                            # Calculate concatenated genome coordinates
                            concat_start = genome_offset + contig_start
                            concat_end = genome_offset + contig_end
                            
                            repeat_regions.append((concat_start, concat_end))
                            
                            debug_print(f"  Line {line_num}: {contig_name}:{contig_start}-{contig_end} -> concatenated {concat_start:,}-{concat_end:,}")
                        else:
                            print(f"Warning: Contig '{contig_name}' (from '{contig_full_name}') not found in mapping. Available contigs: {list(contig_mapping.keys())[:3]}...")
                    except ValueError as e:
                        print(f"Warning: Could not parse line {line_num}: {line} - {e}")
                        continue
    except Exception as e:
        print(f"Error parsing Red output file {red_output_file}: {e}")
    
    debug_print(f"Parsed {len(repeat_regions)} repeat regions from {red_output_file}")
    return repeat_regions

# LEGACY: Keep the old function for backwards compatibility (but it's broken)
def parse_red_output(red_output_file: Path) -> List[Tuple[int, int]]:
    """
    LEGACY FUNCTION - BROKEN COORDINATE MAPPING
    Parse Red repeat output file and return list of repeat regions.
    Red BED format: chrName\tstart\tend (tab-separated, 0-based coordinates)
    
    WARNING: This function has the coordinate mapping bug!
    """
    repeat_regions = []
    
    if not red_output_file.exists():
        print(f"Warning: Red output file {red_output_file} does not exist")
        return repeat_regions
    
    try:
        with open(red_output_file, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                
                # Red BED format: chrName\tstart\tend
                parts = line.split('\t')
                if len(parts) >= 3:
                    try:
                        start = int(parts[1])
                        end = int(parts[2])
                        # BED format is 0-based, keep as-is since we need 0-based for Python
                        repeat_regions.append((start, end))
                    except ValueError:
                        continue
    except Exception as e:
        print(f"Error parsing Red output file {red_output_file}: {e}")
    
    return repeat_regions

# Create contig mapping for coordinate conversion
def create_contig_mapping(records: List[SeqRecord.SeqRecord]) -> Dict[str, int]:
    """
    Create a mapping from contig names to their start positions in the concatenated genome.
    
    Args:
        records: List of SeqRecord objects in the order they appear in the concatenated genome
        
    Returns:
        Dictionary mapping contig IDs to their start positions in concatenated genome
    """
    contig_mapping = {}
    cumulative_pos = 0
    
    for record in records:
        contig_mapping[record.id] = cumulative_pos
        debug_print(f"  Contig mapping: {record.id} -> start position {cumulative_pos:,}")
        cumulative_pos += len(record.seq)
    
    debug_print(f"Created contig mapping for {len(contig_mapping)} contigs")
    return contig_mapping

# FIXED: Simple Red processing with proper coordinate mapping
def find_repeats_red_simple_fixed(fasta_path: Path, records: List[SeqRecord.SeqRecord], num_threads: int = 8) -> List[Tuple[int, int]]:
    """Fixed version - Fast path for small genomes with proper coordinate mapping"""
    
    # Create contig mapping for coordinate conversion
    contig_mapping = create_contig_mapping(records)
    
    with tempfile.TemporaryDirectory() as temp_dir:
        temp_path = Path(temp_dir)
        genome_dir = temp_path / "genome"
        output_dir = temp_path / "output"
        genome_dir.mkdir()
        output_dir.mkdir()
        
        # Simple copy with extension change
        temp_fasta = genome_dir / f"{fasta_path.stem}.fa"
        shutil.copy2(fasta_path, temp_fasta)
        
        print(f"Running Red directly on {temp_fasta.name}...")
        
        try:
            red_cmd = [
                "Red",
                "-gnm", str(genome_dir),
                "-rpt", str(output_dir), 
                "-frm", "2",
                "-cor", str(num_threads)
            ]
            
            result = subprocess.run(red_cmd, capture_output=True, text=True, timeout=1800)
            
            if result.returncode != 0:
                print(f"Red failed: {result.stderr}")
                return []
            
            # Parse output files with proper coordinate mapping
            red_output_files = list(output_dir.glob("*.bed"))
            all_repeats = []
            for red_file in red_output_files:
                file_repeats = parse_red_output_with_mapping(red_file, contig_mapping)
                all_repeats.extend(file_repeats)
            
            print(f"Found {len(all_repeats)} repeat regions with proper coordinate mapping")
            return all_repeats
            
        except Exception as e:
            print(f"Error running Red: {e}")
            return []

# FIXED: Find repeat regions using Red with proper coordinate mapping
def find_repeats_red_fixed(fasta_path: Path, num_threads: int = 8, 
                          max_contig_size: int = 1_000_000, fragment_size: int = 500_000, 
                          overlap: int = 50_000) -> List[Tuple[int, int]]:
    
    # Read original records for coordinate mapping
    records = list(SeqIO.parse(fasta_path, "fasta"))
    
    # Check if any contig needs fragmentation (regardless of total file size)
    needs_fragmentation = any(len(rec.seq) > max_contig_size for rec in records)
    
    if not needs_fragmentation:
        # Only use fast path if no contigs exceed max_contig_size
        print(f"No large contigs detected (all <{max_contig_size:,} bp), using fast path with coordinate mapping...")
        return find_repeats_red_simple_fixed(fasta_path, records, num_threads)
    
    total_length = sum(len(rec.seq) for rec in records)
    max_contig_length = max(len(rec.seq) for rec in records) if records else 0
    
    print(f"Processing {fasta_path.name}: {len(records)} contigs, max={max_contig_length:,} bp, total={total_length:,} bp")
    print(f"Fragmentation needed: {needs_fragmentation} (max contig: {max_contig_length:,} bp > {max_contig_size:,} bp threshold)")
    
    # Create temporary directories for Red processing
    with tempfile.TemporaryDirectory() as temp_dir:
        temp_path = Path(temp_dir)
        genome_dir = temp_path / "genome"
        output_dir = temp_path / "output"
        genome_dir.mkdir()
        output_dir.mkdir()
        
        if needs_fragmentation:
            print(f"Large contigs detected (>{max_contig_size:,} bp), fragmenting for Red processing...")
            
            # Fragment large contigs
            fragmented_records, fragment_mapping = fragment_large_contigs(
                records, max_contig_size, fragment_size, overlap
            )
            
            # Write fragmented sequences to temporary file
            temp_fasta = genome_dir / f"{fasta_path.stem}_fragmented.fa"
            with open(temp_fasta, 'w') as f:
                SeqIO.write(fragmented_records, f, "fasta")
            print(f"Written fragmented genome to {temp_fasta} ({os.path.getsize(temp_fasta)} bytes)")
            
            # Create mapping for fragmented sequences
            contig_mapping = create_contig_mapping(fragmented_records)
            
        else:
            # Use original approach for genomes without large contigs
            print(f"No large contigs detected, using standard Red processing with coordinate mapping...")
            temp_fasta = genome_dir / f"{fasta_path.stem}.fa"
            with open(temp_fasta, 'w') as f:
                SeqIO.write(records, f, "fasta")
            print(f"Written genome to {temp_fasta} ({os.path.getsize(temp_fasta)} bytes)")
            fragment_mapping = []
            fragmented_records = records
            contig_mapping = create_contig_mapping(records)
        
        try:
            # Verify files exist before running Red
            fa_files = list(genome_dir.glob("*.fa"))
            print(f"Found {len(fa_files)} .fa files in {genome_dir}: {[f.name for f in fa_files]}")
            
            if not fa_files:
                print(f"ERROR: No .fa files found in {genome_dir}")
                return []
            
            # Run Red to detect repeats
            print(f"Running Red with {num_threads} threads...")
            
            red_cmd = [
                "Red",
                "-gnm", str(genome_dir),
                "-rpt", str(output_dir),
                "-frm", "2",  # Use format 2: chrName\tstart\tend
                "-cor", str(num_threads)
            ]
            
            print(f"Red command: {' '.join(red_cmd)}")
            
            # Ensure Red runs with a clean environment and proper error handling
            result = subprocess.run(
                red_cmd,
                capture_output=True,
                text=True,
                timeout=1800,  # 30 minute timeout
                cwd=temp_dir,
                env=dict(os.environ, TMPDIR=temp_dir)  # Force temp operations to use local dir
            )
            
            print(f"Red completed with return code {result.returncode}")
            print(f"Red stdout: {result.stdout}")
            
            if result.returncode != 0:
                print(f"Red failed with return code {result.returncode}")
                print(f"Red stderr: {result.stderr}")
                
                # Check for specific known Red issues
                stderr_lower = result.stderr.lower()
                if "the size of the genome cannot be zero" in stderr_lower or "invalid state exception" in stderr_lower:
                    print(f"Red encountered genome size calculation error - this is a known Red issue with certain genomes")
                
                return []
            
            # Find Red output file (Red creates .bed files, not .rpt)
            red_output_files = list(output_dir.glob("*.bed"))
            if not red_output_files:
                print(f"No Red output files found in {output_dir}")
                return []
            
            # Parse all Red output files with proper coordinate mapping
            all_repeats = []
            for red_file in red_output_files:
                file_repeats = parse_red_output_with_mapping(red_file, contig_mapping)
                all_repeats.extend(file_repeats)
            
            if needs_fragmentation and fragment_mapping:
                # Map fragment coordinates back to original genome coordinates
                original_repeats = map_fragment_coords_to_original(
                    all_repeats, fragment_mapping, fragmented_records
                )
                print(f"Found {len(all_repeats)} repeat regions in fragments -> {len(original_repeats)} in original coordinates")
                return original_repeats
            else:
                print(f"Found {len(all_repeats)} repeat regions using Red with proper coordinate mapping")
                return all_repeats
                
        except subprocess.TimeoutExpired:
            print(f"Red timeout for {fasta_path.name} after 30 minutes - skipping repeat detection")
            return []
        except Exception as e:
            print(f"Error running Red: {e}")
            return []

# Merge overlapping regions
def merge_overlapping_regions(regions: List[Tuple[int, int]]) -> List[Tuple[int, int]]:
    if not regions:
        return []
    
    # Sort regions by start position
    sorted_regions = sorted(regions)
    merged = [sorted_regions[0]]
    
    for current_start, current_end in sorted_regions[1:]:
        last_start, last_end = merged[-1]
        
        # If current region overlaps with the last merged region
        if current_start <= last_end:
            # Merge by extending the end position
            merged[-1] = (last_start, max(last_end, current_end))
        else:
            # No overlap, add as new region
            merged.append((current_start, current_end))
    
    return merged

# Select candidate regions based on GC% and repeats
def select_candidate_regions(gc_windows, gc_thresh, repeats) -> Tuple[List[Tuple[int, int]], List[Tuple[int, int]], List[Tuple[int, int]]]:
    gc_only = []
    repeat_only = []
    both = []

    # Get high GC regions
    high_gc_regions = [(start, end) for start, end, gc in gc_windows if gc >= gc_thresh]
    
    # Helper function to check if two regions overlap
    def regions_overlap(region1, region2):
        start1, end1 = region1
        start2, end2 = region2
        return start1 < end2 and start2 < end1
    
    # Track which regions have been classified
    gc_classified = set()
    repeat_classified = set()
    
    # Find overlapping regions (both GC and repeat)
    for gc_idx, gc_region in enumerate(high_gc_regions):
        for repeat_idx, repeat_region in enumerate(repeats):
            if regions_overlap(gc_region, repeat_region):
                # These regions overlap - add to "both" category
                if gc_region not in both:
                    both.append(gc_region)
                if repeat_region not in both:
                    both.append(repeat_region)
                gc_classified.add(gc_idx)
                repeat_classified.add(repeat_idx)
    
    # Add remaining high GC regions that don't overlap with repeats
    for gc_idx, gc_region in enumerate(high_gc_regions):
        if gc_idx not in gc_classified:
            gc_only.append(gc_region)
    
    # Add remaining repeat regions that don't overlap with high GC
    for repeat_idx, repeat_region in enumerate(repeats):
        if repeat_idx not in repeat_classified:
            repeat_only.append(repeat_region)
    
    # Remove duplicates and sort
    both = sorted(list(set(both)))
    
    return sorted(gc_only), sorted(repeat_only), both

# Remove blocks from genome with fallback to random non-overlapping blocks - FIXED VERSION
def remove_blocks(seq: str, preferred_blocks: List[Tuple[int, int]], block_categories: dict, completeness: float, seed: int, debug: bool = False) -> Tuple[List[str], dict, List[dict]]:
    random.seed(seed)
    genome_len = len(seq)
    to_remove = int((1 - completeness) * genome_len)
    removed = 0
    mask = np.ones(len(seq), dtype=bool)
    used_blocks = set()
    
    if debug:
        print(f"DEBUG: Initial setup - genome_len={genome_len}, to_remove={to_remove}, completeness={completeness}")
        print(f"DEBUG: Initial mask True count: {np.sum(mask)}")
    
    # Track removal statistics
    removal_stats = {
        'gc_only': 0,
        'repeat_only': 0,
        'both': 0,
        'random': 0,
        'total_removed': 0,
        'total_bases': genome_len,
        'target_removal': to_remove
    }
    
    # Track removed regions with their details
    removed_regions = []

    # Primary removal: from preferred blocks
    for start, end in preferred_blocks:
        if removed >= to_remove:
            break
        if (start, end) in used_blocks:
            continue
        
        # Check if this block would exceed our target
        block_size = end - start
        if removed + block_size > to_remove:
            # Skip this block if it would cause us to exceed target
            continue
            
        mask[start:end] = False
        removed += block_size
        
        if debug:
            print(f"DEBUG: Removed targeted block ({start}, {end}), size={block_size}, total_removed={removed}")
            print(f"DEBUG: Mask False count after removal: {np.sum(~mask)}")
        
        # Track which category this block belongs to
        category = block_categories.get((start, end), 'unknown')
        removal_stats[category] += block_size
        
        # Record the removed region
        removed_regions.append({
            'start': start,
            'end': end,
            'length': block_size,
            'category': category,
            'type': 'targeted'
        })
        
        used_blocks.add((start, end))

    # Secondary removal: fill in remaining with random blocks adapted to available fragments
    max_attempts = 10000  # Prevent infinite loop
    attempts = 0
    
    while removed < to_remove and attempts < max_attempts:
        attempts += 1
        
        # Calculate remaining bases to remove
        remaining_to_remove = to_remove - removed
        
        # Find all available contiguous fragments
        available_fragments = []
        frag_start = None
        for i, keep in enumerate(mask):
            if keep and frag_start is None:
                frag_start = i
            elif not keep and frag_start is not None:
                frag_length = i - frag_start
                if frag_length > 0:
                    available_fragments.append((frag_start, i, frag_length))
                frag_start = None
        if frag_start is not None:
            frag_length = genome_len - frag_start
            if frag_length > 0:
                available_fragments.append((frag_start, genome_len, frag_length))
        
        # If no fragments available, break
        if not available_fragments:
            if debug:
                print(f"DEBUG: No available fragments remaining. Breaking at removed={removed}, target={to_remove}")
            break
        
        # Sort fragments by size (largest first)
        available_fragments.sort(key=lambda x: x[2], reverse=True)
        
        # Adaptively choose block size based on available fragments
        largest_available = available_fragments[0][2]
        
        # If the largest available fragment is very small, consider removing entire small fragments
        if largest_available < 1000:
            # Remove small fragments entirely to reach target
            for frag_start, frag_end, frag_length in available_fragments:
                if removed >= to_remove:
                    break
                if frag_length <= remaining_to_remove:
                    mask[frag_start:frag_end] = False
                    removed += frag_length
                    removal_stats['random'] += frag_length
                    
                    if debug:
                        print(f"DEBUG: Removed entire small fragment ({frag_start}, {frag_end}), size={frag_length}, total_removed={removed}")
                    
                    removed_regions.append({
                        'start': frag_start,
                        'end': frag_end,
                        'length': frag_length,
                        'category': 'random',
                        'type': 'random'
                    })
            continue
        
        # For larger fragments, choose an appropriate block size
        max_block_size = min(20000, remaining_to_remove, largest_available)
        min_block_size = min(1000, max_block_size)  # Reduced minimum to 1000
        
        # If we need less than min_block_size, just take what we need
        if remaining_to_remove < min_block_size:
            block_size = remaining_to_remove
        else:
            block_size = random.randint(min_block_size, max_block_size)
        
        # Try to find a suitable position within available fragments
        suitable_fragments = [f for f in available_fragments if f[2] >= block_size]
        
        if not suitable_fragments:
            # If no fragment is large enough, use the largest available
            frag_start, frag_end, frag_length = available_fragments[0]
            if frag_length <= remaining_to_remove:
                mask[frag_start:frag_end] = False
                removed += frag_length
                removal_stats['random'] += frag_length
                
                if debug:
                    print(f"DEBUG: Removed largest available fragment ({frag_start}, {frag_end}), size={frag_length}, total_removed={removed}")
                
                removed_regions.append({
                    'start': frag_start,
                    'end': frag_end,
                    'length': frag_length,
                    'category': 'random',
                    'type': 'random'
                })
            else:
                # Fragment is too large, but we need to remove something
                # Remove exactly what we need from the beginning of this fragment
                mask[frag_start:frag_start + remaining_to_remove] = False
                removed += remaining_to_remove
                removal_stats['random'] += remaining_to_remove
                
                if debug:
                    print(f"DEBUG: Partially removed from fragment ({frag_start}, {frag_start + remaining_to_remove}), size={remaining_to_remove}, total_removed={removed}")
                
                removed_regions.append({
                    'start': frag_start,
                    'end': frag_start + remaining_to_remove,
                    'length': remaining_to_remove,
                    'category': 'random',
                    'type': 'random'
                })
            continue
        
        # Choose a random fragment from suitable ones
        chosen_frag = random.choice(suitable_fragments)
        frag_start, frag_end, frag_length = chosen_frag
        
        # Choose a random position within this fragment
        max_start_pos = frag_end - block_size
        start = random.randint(frag_start, max_start_pos)
        end = start + block_size
        
        # Remove the block
        mask[start:end] = False
        removed += block_size
        removal_stats['random'] += block_size
        
        if debug:
            print(f"DEBUG: Removed random block ({start}, {end}), size={block_size}, total_removed={removed}")
            print(f"DEBUG: Mask False count after random removal: {np.sum(~mask)}")
        
        # Record the random removal
        removed_regions.append({
            'start': start,
            'end': end,
            'length': block_size,
            'category': 'random',
            'type': 'random'
        })
    
    if attempts >= max_attempts:
        print(f"WARNING: Reached maximum attempts ({max_attempts}) in random block removal. Removed {removed}/{to_remove} bases.")

    removal_stats['total_removed'] = removed
    
    if debug:
        print(f"DEBUG: Final mask status - True count: {np.sum(mask)}, False count: {np.sum(~mask)}")
        print(f"DEBUG: Expected to keep: {genome_len - to_remove}, Expected to remove: {to_remove}")
    
    remaining_seq = []
    frag_start = None
    for i, keep in enumerate(mask):
        if keep and frag_start is None:
            frag_start = i
        elif not keep and frag_start is not None:
            remaining_seq.append(seq[frag_start:i])
            frag_start = None
    if frag_start is not None:
        remaining_seq.append(seq[frag_start:])

    if debug:
        actual_kept = sum(len(frag) for frag in remaining_seq)
        print(f"DEBUG: Fragment reconstruction - kept {actual_kept} bases in {len(remaining_seq)} fragments")
        print(f"DEBUG: Fragments lengths: {[len(f) for f in remaining_seq]}")

    return remaining_seq, removal_stats, removed_regions


# Create circular genome visualization using pyCirclize
def create_circular_genome_plot(genome_length: int, removed_regions: List[dict], removal_stats: dict, 
                              output_pdf: str, genome_name: str = "Genome"):
    """
    Create a circular genome plot showing removed regions using pyCirclize.
    """
    # Define colors for each category
    colors = {
        'gc_only': '#FF6B6B',      # Red
        'repeat_only': '#7e57c2',   # Purple
        'both': '#9ccc65',          # Light green
        'random': '#FFE66D'         # Yellow
    }
    
    # Initialize Circos plot
    sectors = {genome_name: genome_length}
    circos = Circos(sectors, space=2)
    
    # Add sector track for the genome
    sector = circos.sectors[0]
    
    # Add genome track (backbone)
    track = sector.add_track((90, 100))
    track.axis(fc="lightgray", ec="none", lw=0)
    
    # Add removed regions as colored bars
    for region in removed_regions:
        # Calculate angular positions
        start_angle = (region['start'] / genome_length) * 360
        end_angle = (region['end'] / genome_length) * 360
        
        # Add colored arc for removed region
        track.rect(region['start'], region['end'], fc=colors[region['category']], 
                  ec="none", lw=0, alpha=0.8)
    
    # Add tick marks and labels
    major_ticks = []
    tick_interval = 10 ** (len(str(genome_length)) - 1)  # Dynamic tick interval based on genome size
    for i in range(0, genome_length, tick_interval):
        major_ticks.append(i)
    
    # Add axis with ticks
    track.xticks(major_ticks, labels=[f"{int(t/1000)}kb" if t > 0 else "0" for t in major_ticks], 
                label_size=8, label_orientation="vertical")
    
    # Add title and statistics
    title_text = f"{genome_name}\n"
    title_text += f"Genome size: {genome_length:,} bp | "
    title_text += f"Removed: {removal_stats['total_removed']:,} bp ({removal_stats['total_removed']/genome_length*100:.1f}%)"
    
    # Create figure with title
    fig = circos.plotfig()
    fig.suptitle(title_text, fontsize=14, y=0.98)
    
    # Add legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor=colors['gc_only'], label="High GC only"),
        Patch(facecolor=colors['repeat_only'], label="Repeat only"),
        Patch(facecolor=colors['both'], label="Both GC+Repeat"),
        Patch(facecolor=colors['random'], label="Random")
    ]
    fig.legend(handles=legend_elements, loc='lower center', ncol=4, 
              bbox_to_anchor=(0.5, -0.05), fontsize=10)
    
    # Save figure
    plt.tight_layout()
    plt.savefig(output_pdf, format='pdf', dpi=300, bbox_inches='tight')
    plt.close()

# Main function - FIXED VERSION
def simulate_draft(input_fasta: str, output_fasta: str, completeness: float, seed: int, num_threads: int = 8, 
                  max_contig_size: int = 1_000_000, fragment_size: int = 500_000, overlap: int = 50_000,
                  log_file: str = None, create_visualization: bool = False):
    debug_print("Starting simulate_draft function (FIXED VERSION)")
    debug_print(f"Reading input FASTA: {input_fasta}")
    
    # Check if file exists and is readable
    if not os.path.exists(input_fasta):
        raise FileNotFoundError(f"Input file not found: {input_fasta}")
    
    file_size = os.path.getsize(input_fasta)
    debug_print(f"Input file size: {file_size:,} bytes")
    
    # Parse FASTA with progress
    debug_print("Parsing FASTA file...")
    try:
        records = []
        for i, record in enumerate(SeqIO.parse(input_fasta, "fasta")):
            records.append(record)
            debug_print(f"  Read record {i+1}: {record.id} ({len(record.seq):,} bp)")
            if i > 0 and i % 100 == 0:
                debug_print(f"  ... processed {i} records so far")
    except Exception as e:
        debug_print(f"ERROR parsing FASTA: {e}")
        raise
    
    debug_print(f"Loaded {len(records)} records from FASTA")
    num_original_contigs = len(records)
    
    # More efficient concatenation with progress tracking
    debug_print("Concatenating sequences...")
    seq_parts = []
    total_length = 0
    for i, rec in enumerate(records):
        rec_seq = str(rec.seq)
        seq_parts.append(rec_seq)
        total_length += len(rec_seq)
        debug_print(f"  Processed record {i+1}/{len(records)}: {rec.id} ({len(rec_seq):,} bp)")
    
    debug_print(f"Joining {len(seq_parts)} sequences...")
    seq = ''.join(seq_parts)
    debug_print(f"Total sequence length: {len(seq):,} bp")
    genome_name = Path(input_fasta).stem

    # Adjust window size for small genomes
    debug_print("Calculating GC windows...")
    window_size = min(10000, len(seq) // 10)  # Use 1/10 of genome or 10kb, whichever is smaller
    window_size = max(window_size, 100)  # But at least 100 bases
    debug_print(f"Using window size: {window_size}")
    
    gc_windows = compute_gc_windows(seq, window_size)
    debug_print(f"Computed {len(gc_windows)} GC windows")
    gc_vals = [gc for _, _, gc in gc_windows]
    
    # Handle empty gc_vals case
    debug_print("Computing GC statistics...")
    if gc_vals:
        mean_gc = np.mean(gc_vals)
        std_gc = np.std(gc_vals)
        gc_thresh = mean_gc + 2 * std_gc
        debug_print(f"GC stats - mean: {mean_gc:.3f}, std: {std_gc:.3f}, threshold: {gc_thresh:.3f}")
    else:
        # If genome is too small for any windows, use whole genome GC
        gc_count = seq.count("G") + seq.count("C")
        mean_gc = gc_count / len(seq) if len(seq) > 0 else 0
        std_gc = 0
        gc_thresh = 1.0  # Set high threshold so no regions are selected by GC
        debug_print(f"Small genome - whole GC: {mean_gc:.3f}, threshold: {gc_thresh:.3f}")

    debug_print(f"Starting Red analysis with {num_threads} threads (FIXED VERSION)...")
    repeats = find_repeats_red_fixed(Path(input_fasta), num_threads=num_threads, 
                              max_contig_size=max_contig_size, fragment_size=fragment_size, overlap=overlap)
    debug_print(f"Red analysis completed, found {len(repeats)} initial repeat regions")
    
    # Show distribution of repeat regions
    if repeats:
        halfway_point = len(seq) // 2
        first_half_repeats = sum(1 for start, end in repeats if start < halfway_point)
        second_half_repeats = len(repeats) - first_half_repeats
        debug_print(f"Repeat distribution: {first_half_repeats} in first half, {second_half_repeats} in second half")
    
    # Merge overlapping repeat regions to ensure non-overlapping requirement
    debug_print("Merging overlapping repeat regions...")
    repeats = merge_overlapping_regions(repeats)
    debug_print(f"After merging: {len(repeats)} repeat regions")
    
    debug_print("Selecting candidate regions...")
    gc_only, repeat_only, both = select_candidate_regions(gc_windows, gc_thresh, repeats)
    debug_print(f"Candidate regions - GC only: {len(gc_only)}, repeat only: {len(repeat_only)}, both: {len(both)}")
    
    # Merge overlapping regions within each category
    debug_print("Merging overlapping regions within categories...")
    gc_only = merge_overlapping_regions(gc_only)
    repeat_only = merge_overlapping_regions(repeat_only)
    both = merge_overlapping_regions(list(both))
    debug_print(f"After category merging - GC only: {len(gc_only)}, repeat only: {len(repeat_only)}, both: {len(both)}")

    # Create a mapping of blocks to their categories
    debug_print("Creating block categories and prioritization...")
    block_categories = {}
    for block in both:
        block_categories[block] = 'both'
    for block in gc_only:
        block_categories[block] = 'gc_only'
    for block in repeat_only:
        block_categories[block] = 'repeat_only'

    prioritized_blocks = both + gc_only + repeat_only
    debug_print(f"Total prioritized blocks: {len(prioritized_blocks)}")

    debug_print("Starting block removal...")
    fragments, removal_stats, removed_regions = remove_blocks(seq, prioritized_blocks, block_categories, completeness, seed, debug=False)
    debug_print(f"Block removal completed, generated {len(fragments)} fragments")

    # Calculate actual output genome size
    debug_print("Calculating output statistics...")
    actual_output_size = sum(len(frag) for frag in fragments)
    actual_removed_size = len(seq) - actual_output_size
    debug_print(f"Output size: {actual_output_size:,} bp, removed: {actual_removed_size:,} bp")

    # Write log file and create visualization if requested
    debug_print("Writing log file and creating output...")
    if log_file:
        with open(log_file, 'w') as f:
            f.write("=== Draft Genome Simulation Report (FIXED VERSION with Red) ===\n\n")
            f.write(f"Input file: {input_fasta}\n")
            f.write(f"Output file: {output_fasta}\n")
            f.write(f"Target completeness: {completeness:.2%}\n")
            f.write(f"Random seed: {seed}\n")
            f.write(f"Number of original contigs: {num_original_contigs}\n")
            f.write(f"Red threads used: {num_threads}\n")
            f.write(f"Fragmentation settings:\n")
            f.write(f"  - Max contig size: {max_contig_size:,} bp\n")
            f.write(f"  - Fragment size: {fragment_size:,} bp\n")
            f.write(f"  - Overlap: {overlap:,} bp\n\n")
            
            f.write("=== Removal Statistics ===\n")
            f.write(f"Total genome size: {removal_stats['total_bases']:,} bases\n")
            f.write(f"Target removal: {removal_stats['target_removal']:,} bases\n")
            f.write(f"Actual removal: {removal_stats['total_removed']:,} bases\n")
            f.write(f"Actual completeness: {(removal_stats['total_bases'] - removal_stats['total_removed']) / removal_stats['total_bases']:.2%}\n\n")
            
            f.write("=== Actual Output Statistics ===\n")
            f.write(f"Actual output genome size: {actual_output_size:,} bases\n")
            f.write(f"Actual removed region size: {actual_removed_size:,} bases\n")
            f.write(f"True completeness (output/original): {actual_output_size / removal_stats['total_bases']:.2%}\n\n")
            
            f.write("=== Bases Removed by Category ===\n")
            if removal_stats['total_removed'] > 0:
                f.write(f"High GC regions only: {removal_stats['gc_only']:,} bases ({removal_stats['gc_only']/removal_stats['total_removed']*100:.1f}% of removed)\n")
                f.write(f"Repeat regions only: {removal_stats['repeat_only']:,} bases ({removal_stats['repeat_only']/removal_stats['total_removed']*100:.1f}% of removed)\n")
                f.write(f"Both high GC and repeat: {removal_stats['both']:,} bases ({removal_stats['both']/removal_stats['total_removed']*100:.1f}% of removed)\n")
                f.write(f"Random regions: {removal_stats['random']:,} bases ({removal_stats['random']/removal_stats['total_removed']*100:.1f}% of removed)\n\n")
            else:
                f.write("No bases were removed.\n\n")
            
            # Show repeat region distribution
            if repeats:
                halfway_point = len(seq) // 2
                first_half_repeats = sum(1 for start, end in repeats if start < halfway_point)
                second_half_repeats = len(repeats) - first_half_repeats
                f.write("=== Repeat Region Distribution ===\n")
                f.write(f"Total repeat regions found: {len(repeats)}\n")
                f.write(f"Repeat regions in first half (0-{halfway_point:,}): {first_half_repeats}\n")
                f.write(f"Repeat regions in second half ({halfway_point:,}-{len(seq):,}): {second_half_repeats}\n\n")
            
            f.write("=== Region Analysis ===\n")
            f.write(f"GC threshold (mean + 2*std): {gc_thresh:.2%}\n")
            f.write(f"Number of high GC-only regions: {len(gc_only)}\n")
            f.write(f"Number of repeat-only regions: {len(repeat_only)}\n")
            f.write(f"Number of both GC+repeat regions: {len(both)}\n")
            f.write(f"Number of fragments in output: {len(fragments)}\n\n")
            
            f.write("=== Detailed List of Removed Regions ===\n")
            f.write(f"Total regions removed: {len(removed_regions)}\n\n")
            
            # Sort regions by start position
            removed_regions.sort(key=lambda x: x['start'])
            
            f.write("Start\tEnd\tLength\tCategory\tType\n")
            f.write("-" * 60 + "\n")
            for region in removed_regions:
                f.write(f"{region['start']:,}\t{region['end']:,}\t{region['length']:,}\t{region['category']}\t{region['type']}\n")
            
            # Summary by category
            f.write("\n=== Removed Regions Summary by Category ===\n")
            category_counts = {'gc_only': 0, 'repeat_only': 0, 'both': 0, 'random': 0}
            for region in removed_regions:
                category_counts[region['category']] += 1
            
            f.write(f"High GC only regions: {category_counts['gc_only']}\n")
            f.write(f"Repeat only regions: {category_counts['repeat_only']}\n")
            f.write(f"Both GC+repeat regions: {category_counts['both']}\n")
            f.write(f"Random regions: {category_counts['random']}\n")
        
        # Create PDF visualization only if requested
        if create_visualization:
            if log_file.endswith('.txt'):
                pdf_path = log_file.replace('.txt', '_visualization.pdf')
            else:
                pdf_path = log_file + '_visualization.pdf'
            
            print(f"Creating circular visualization: {pdf_path}")
            try:
                create_circular_genome_plot(len(seq), removed_regions, removal_stats, pdf_path, genome_name)
                debug_print(f"Visualization saved to: {pdf_path}")
            except Exception as e:
                print(f"Warning: Failed to create visualization: {e}")
                debug_print(f"Visualization error: {e}")

    debug_print("Creating output FASTA records...")
    frag_records = []
    for i, frag in enumerate(fragments):
        new_record = SeqRecord.SeqRecord(
            Seq(frag),
            id=f"{genome_name}_frag_{i+1}",
            description=""
        )
        frag_records.append(new_record)

    debug_print(f"Writing {len(frag_records)} fragments to {output_fasta}")
    SeqIO.write(frag_records, output_fasta, "fasta")
    debug_print("simulate_draft function completed successfully! (FIXED VERSION)")

if __name__ == "__main__":
    debug_print("Script started, parsing arguments... (FIXED VERSION)")
    parser = argparse.ArgumentParser(description="FIXED: Simulate realistic draft bacterial genome using Red for repeat detection with circular visualization")
    parser.add_argument("--input", required=True, help="Path to complete genome FASTA")
    parser.add_argument("--output", required=True, help="Output FASTA of simulated draft")
    parser.add_argument("--completeness", type=float, default=0.5, help="Target completeness (0-100 for percentage, or 0-1 for fraction)")
    parser.add_argument("--seed", type=int, default=42, help="Random seed for reproducibility")
    parser.add_argument("--num_threads", type=int, default=8, help="Number of threads for Red (default: 8)")
    parser.add_argument("--max_contig_size", type=int, default=1_000_000, help="Maximum contig size before fragmentation (default: 1MB)")
    parser.add_argument("--fragment_size", type=int, default=500_000, help="Fragment size for large contigs (default: 500kb)")
    parser.add_argument("--overlap", type=int, default=50_000, help="Overlap between fragments (default: 50kb)")
    parser.add_argument("--log", help="Path to log file (text report only)")
    parser.add_argument("--vis_log", help="Path to log file with circular PDF visualization")
    
    debug_print("Parsing command line arguments...")
    args = parser.parse_args()
    debug_print("Arguments parsed successfully!")
    debug_print(f"Input: {args.input}")
    debug_print(f"Output: {args.output}")  
    debug_print(f"Completeness: {args.completeness}")
    debug_print(f"Threads: {args.num_threads}")

    # Convert completeness to fraction if given as percentage
    debug_print("Processing completeness value...")
    completeness = args.completeness
    if completeness > 1:
        completeness = completeness / 100.0
    debug_print(f"Final completeness: {completeness}")
    
    # Validate completeness is in valid range
    if completeness < 0 or completeness > 1:
        raise ValueError(f"Completeness must be between 0 and 1 (or 0-100 for percentage), got {args.completeness}")
    
    # Determine which log file to use and whether to create visualization
    debug_print("Setting up log file and visualization options...")
    log_file = args.log or args.vis_log
    create_vis = bool(args.vis_log)
    debug_print(f"Log file: {log_file}, Create visualization: {create_vis}")
    
    # Get number of threads from environment or use provided value
    debug_print("Determining thread count...")
    num_threads = args.num_threads
    if 'SLURM_CPUS_PER_TASK' in os.environ:
        num_threads = int(os.environ['SLURM_CPUS_PER_TASK'])
        debug_print(f"Using {num_threads} threads from SLURM allocation")
    else:
        debug_print(f"Using {num_threads} threads from argument")
    
    debug_print("About to call simulate_draft function... (FIXED VERSION)")
    try:
        simulate_draft(
            input_fasta=args.input,
            output_fasta=args.output,
            completeness=completeness,
            seed=args.seed,
            num_threads=num_threads,
            max_contig_size=args.max_contig_size,
            fragment_size=args.fragment_size,
            overlap=args.overlap,
            log_file=log_file,
            create_visualization=create_vis
        )
        debug_print("Script completed successfully! (FIXED VERSION)")
    except Exception as e:
        debug_print(f"Script failed with error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)