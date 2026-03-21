# draft_genome_simulator_with_visualization_red_circle_fixed.py
import argparse
from bisect import bisect_left
import json
import random
import os
import sys
import tempfile
import shutil
from pathlib import Path
from typing import Dict, List, Optional, Tuple
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
    if shutil.which("Red") is None:
        print("Warning: Red is not available in PATH; skipping Red repeat detection")
        return []
    
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

def find_self_homologies_nucmer(fasta_path: Path, min_identity: float = 95.0,
                                min_length: int = 500) -> List[Tuple[int, int]]:
    """
    Detect exact and near-exact self-homologies with nucmer.

    Coordinates are mapped onto the same concatenated coordinate system used by
    the Red repeat detector.
    """
    missing_tools = [tool for tool in ("nucmer", "show-coords") if shutil.which(tool) is None]
    if missing_tools:
        missing = ", ".join(missing_tools)
        raise RuntimeError(f"--use-nucmer requires the following tools in PATH: {missing}")

    records = list(SeqIO.parse(fasta_path, "fasta"))
    contig_mapping = create_contig_mapping(records)
    repeat_regions = set()

    with tempfile.TemporaryDirectory() as temp_dir:
        temp_path = Path(temp_dir)
        prefix = temp_path / "nucmer_self"
        delta_path = prefix.with_suffix(".delta")

        nucmer_cmd = [
            "nucmer",
            "--maxmatch",
            "-p", str(prefix),
            str(fasta_path),
            str(fasta_path)
        ]

        debug_print(f"Running nucmer self-alignment: {' '.join(nucmer_cmd)}")
        nucmer_result = subprocess.run(
            nucmer_cmd,
            capture_output=True,
            text=True,
            timeout=1800
        )

        if nucmer_result.returncode != 0:
            raise RuntimeError(
                "nucmer failed with return code "
                f"{nucmer_result.returncode}: {nucmer_result.stderr.strip()}"
            )

        coords_cmd = ["show-coords", "-THrcl", str(delta_path)]
        debug_print(f"Parsing nucmer alignments: {' '.join(coords_cmd)}")
        coords_result = subprocess.run(
            coords_cmd,
            capture_output=True,
            text=True,
            timeout=1800
        )

        if coords_result.returncode != 0:
            raise RuntimeError(
                "show-coords failed with return code "
                f"{coords_result.returncode}: {coords_result.stderr.strip()}"
            )

        for line_num, line in enumerate(coords_result.stdout.splitlines(), 1):
            if not line.strip():
                continue

            parts = line.split("\t")
            if len(parts) < 9:
                continue

            try:
                ref_start, ref_end, qry_start, qry_end = map(int, parts[:4])
                aln_len_ref = int(parts[4])
                aln_len_qry = int(parts[5])
                identity = float(parts[6])
                ref_name = parts[-2].split()[0]
                qry_name = parts[-1].split()[0]
            except ValueError:
                debug_print(f"Skipping unparsable nucmer line {line_num}: {line}")
                continue

            if min(aln_len_ref, aln_len_qry) < min_length or identity < min_identity:
                continue

            ref_interval = (min(ref_start, ref_end) - 1, max(ref_start, ref_end))
            qry_interval = (min(qry_start, qry_end) - 1, max(qry_start, qry_end))

            # Skip the trivial self-diagonal alignment for a sequence against itself.
            if ref_name == qry_name and ref_interval == qry_interval:
                continue

            if ref_name not in contig_mapping or qry_name not in contig_mapping:
                debug_print(
                    f"Skipping nucmer alignment with unmapped contigs: {ref_name}, {qry_name}"
                )
                continue

            repeat_regions.add(
                (contig_mapping[ref_name] + ref_interval[0], contig_mapping[ref_name] + ref_interval[1])
            )
            repeat_regions.add(
                (contig_mapping[qry_name] + qry_interval[0], contig_mapping[qry_name] + qry_interval[1])
            )

    merged_regions = merge_overlapping_regions(sorted(repeat_regions))
    debug_print(f"nucmer detected {len(merged_regions)} self-homology regions after merging")
    return merged_regions

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

def load_input_records(input_fasta: str) -> List[SeqRecord.SeqRecord]:
    debug_print(f"Reading input FASTA: {input_fasta}")

    if not os.path.exists(input_fasta):
        raise FileNotFoundError(f"Input file not found: {input_fasta}")

    file_size = os.path.getsize(input_fasta)
    debug_print(f"Input file size: {file_size:,} bytes")

    debug_print("Parsing FASTA file...")
    records = []
    try:
        for i, record in enumerate(SeqIO.parse(input_fasta, "fasta")):
            records.append(record)
            debug_print(f"  Read record {i+1}: {record.id} ({len(record.seq):,} bp)")
            if i > 0 and i % 100 == 0:
                debug_print(f"  ... processed {i} records so far")
    except Exception as exc:
        debug_print(f"ERROR parsing FASTA: {exc}")
        raise

    debug_print(f"Loaded {len(records)} records from FASTA")
    return records

def concatenate_records(records: List[SeqRecord.SeqRecord]) -> Tuple[str, List[Dict[str, int]]]:
    debug_print("Concatenating sequences...")
    seq_parts = []
    record_layout = []
    cumulative_pos = 0

    for i, record in enumerate(records):
        rec_seq = str(record.seq)
        start = cumulative_pos
        end = start + len(rec_seq)
        seq_parts.append(rec_seq)
        record_layout.append({
            "id": record.id,
            "start": start,
            "end": end,
            "length": len(rec_seq),
        })
        cumulative_pos = end
        debug_print(f"  Processed record {i+1}/{len(records)}: {record.id} ({len(rec_seq):,} bp)")

    debug_print(f"Joining {len(seq_parts)} sequences...")
    seq = ''.join(seq_parts)
    debug_print(f"Total sequence length: {len(seq):,} bp")
    return seq, record_layout

def compute_gc_statistics(seq: str) -> Tuple[int, List[Tuple[int, int, float]], float, float, float]:
    debug_print("Calculating GC windows...")
    window_size = min(10000, len(seq) // 10)
    window_size = max(window_size, 100)
    debug_print(f"Using window size: {window_size}")

    gc_windows = compute_gc_windows(seq, window_size)
    debug_print(f"Computed {len(gc_windows)} GC windows")
    gc_vals = [gc for _, _, gc in gc_windows]

    debug_print("Computing GC statistics...")
    if gc_vals:
        mean_gc = float(np.mean(gc_vals))
        std_gc = float(np.std(gc_vals))
        gc_thresh = mean_gc + 2 * std_gc
        debug_print(f"GC stats - mean: {mean_gc:.3f}, std: {std_gc:.3f}, threshold: {gc_thresh:.3f}")
    else:
        gc_count = seq.count("G") + seq.count("C")
        mean_gc = gc_count / len(seq) if len(seq) > 0 else 0.0
        std_gc = 0.0
        gc_thresh = 1.0
        debug_print(f"Small genome - whole GC: {mean_gc:.3f}, threshold: {gc_thresh:.3f}")

    return window_size, gc_windows, mean_gc, std_gc, gc_thresh

def calculate_n50(lengths: List[int]) -> Tuple[int, int]:
    if not lengths:
        return 0, 0

    sorted_lengths = sorted(lengths, reverse=True)
    half_total = sum(sorted_lengths) / 2
    cumulative = 0

    for idx, length in enumerate(sorted_lengths, 1):
        cumulative += length
        if cumulative >= half_total:
            return length, idx

    return 0, 0

def calculate_nx(lengths: List[int], fraction: float) -> Tuple[int, int]:
    if not lengths:
        return 0, 0

    sorted_lengths = sorted(lengths, reverse=True)
    target = sum(sorted_lengths) * fraction
    cumulative = 0

    for idx, length in enumerate(sorted_lengths, 1):
        cumulative += length
        if cumulative >= target:
            return length, idx

    return 0, 0

def summarize_contig_lengths(lengths: List[int]) -> Dict[str, object]:
    if not lengths:
        return {
            "contig_count": 0,
            "total_bases": 0,
            "min_contig_length": 0,
            "max_contig_length": 0,
            "mean_contig_length": 0.0,
            "median_contig_length": 0.0,
            "n50": 0,
            "l50": 0,
            "n90": 0,
            "l90": 0,
            "contigs_ge_500bp": 0,
            "contigs_ge_1kb": 0,
            "contigs_ge_10kb": 0,
            "contigs_ge_50kb": 0,
            "top_10_lengths": [],
        }

    sorted_lengths = sorted(lengths, reverse=True)
    n50, l50 = calculate_n50(sorted_lengths)
    n90, l90 = calculate_nx(sorted_lengths, 0.9)

    return {
        "contig_count": len(sorted_lengths),
        "total_bases": int(sum(sorted_lengths)),
        "min_contig_length": int(sorted_lengths[-1]),
        "max_contig_length": int(sorted_lengths[0]),
        "mean_contig_length": float(np.mean(sorted_lengths)),
        "median_contig_length": float(np.median(sorted_lengths)),
        "n50": int(n50),
        "l50": int(l50),
        "n90": int(n90),
        "l90": int(l90),
        "contigs_ge_500bp": int(sum(length >= 500 for length in sorted_lengths)),
        "contigs_ge_1kb": int(sum(length >= 1_000 for length in sorted_lengths)),
        "contigs_ge_10kb": int(sum(length >= 10_000 for length in sorted_lengths)),
        "contigs_ge_50kb": int(sum(length >= 50_000 for length in sorted_lengths)),
        "top_10_lengths": [int(length) for length in sorted_lengths[:10]],
    }

def default_stats_json_path(output_fasta: str) -> str:
    return str(Path(output_fasta).with_suffix('.stats.json'))

def write_contig_stats_json(stats_json_path: str, payload: Dict[str, object]):
    with open(stats_json_path, 'w') as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")

def prepare_simulation_context(input_fasta: str, num_threads: int = 8,
                               max_contig_size: int = 1_000_000,
                               fragment_size: int = 500_000,
                               overlap: int = 50_000,
                               use_nucmer: bool = False,
                               nucmer_min_identity: float = 95.0,
                               nucmer_min_length: int = 500) -> Dict[str, object]:
    records = load_input_records(input_fasta)
    num_original_contigs = len(records)
    seq, record_layout = concatenate_records(records)
    genome_name = Path(input_fasta).stem

    window_size, gc_windows, mean_gc, std_gc, gc_thresh = compute_gc_statistics(seq)

    debug_print(f"Starting Red analysis with {num_threads} threads (FIXED VERSION)...")
    red_repeats = find_repeats_red_fixed(
        Path(input_fasta),
        num_threads=num_threads,
        max_contig_size=max_contig_size,
        fragment_size=fragment_size,
        overlap=overlap
    )
    debug_print(f"Red analysis completed, found {len(red_repeats)} initial repeat regions")

    nucmer_repeats = []
    if use_nucmer:
        debug_print("Starting nucmer self-alignment analysis...")
        nucmer_repeats = find_self_homologies_nucmer(
            Path(input_fasta),
            min_identity=nucmer_min_identity,
            min_length=nucmer_min_length
        )
        debug_print(f"nucmer analysis completed, found {len(nucmer_repeats)} self-homology regions")

    repeats = merge_overlapping_regions(red_repeats + nucmer_repeats)
    debug_print(f"After merging repeat evidence sources: {len(repeats)} repeat regions")

    if repeats:
        halfway_point = len(seq) // 2
        first_half_repeats = sum(1 for start, _ in repeats if start < halfway_point)
        second_half_repeats = len(repeats) - first_half_repeats
        debug_print(f"Repeat distribution: {first_half_repeats} in first half, {second_half_repeats} in second half")

    debug_print("Selecting candidate regions...")
    gc_only, repeat_only, both = select_candidate_regions(gc_windows, gc_thresh, repeats)
    debug_print(f"Candidate regions - GC only: {len(gc_only)}, repeat only: {len(repeat_only)}, both: {len(both)}")

    debug_print("Merging overlapping regions within categories...")
    gc_only = merge_overlapping_regions(gc_only)
    repeat_only = merge_overlapping_regions(repeat_only)
    both = merge_overlapping_regions(list(both))
    debug_print(f"After category merging - GC only: {len(gc_only)}, repeat only: {len(repeat_only)}, both: {len(both)}")

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

    initial_lengths = [len(record.seq) for record in records]
    initial_n50, initial_l50 = calculate_n50(initial_lengths)
    record_boundaries = [layout["end"] for layout in record_layout[:-1]]

    return {
        "records": records,
        "seq": seq,
        "genome_name": genome_name,
        "num_original_contigs": num_original_contigs,
        "record_layout": record_layout,
        "record_boundaries": record_boundaries,
        "window_size": window_size,
        "gc_windows": gc_windows,
        "mean_gc": mean_gc,
        "std_gc": std_gc,
        "gc_thresh": gc_thresh,
        "red_repeats": red_repeats,
        "nucmer_repeats": nucmer_repeats,
        "repeats": repeats,
        "gc_only": gc_only,
        "repeat_only": repeat_only,
        "both": both,
        "block_categories": block_categories,
        "prioritized_blocks": prioritized_blocks,
        "initial_n50": initial_n50,
        "initial_l50": initial_l50,
    }

def build_contig_intervals(genome_length: int, breakpoints: List[int]) -> List[Tuple[int, int]]:
    cleaned_breakpoints = sorted({point for point in breakpoints if 0 < point < genome_length})
    boundaries = [0] + cleaned_breakpoints + [genome_length]
    return [
        (boundaries[i], boundaries[i + 1])
        for i in range(len(boundaries) - 1)
        if boundaries[i] < boundaries[i + 1]
    ]

def summarize_fragmentation_state(genome_length: int, breakpoints: List[int],
                                  events: List[Dict[str, object]]) -> Dict[str, object]:
    intervals = build_contig_intervals(genome_length, breakpoints)
    lengths = [end - start for start, end in intervals]
    n50, l50 = calculate_n50(lengths)
    random_breakpoints = sum(1 for event in events if event["category"] == "random")

    return {
        "breakpoints": sorted(set(breakpoints)),
        "intervals": intervals,
        "lengths": lengths,
        "n50": n50,
        "l50": l50,
        "contig_count": len(intervals),
        "events": [dict(event) for event in events],
        "random_breakpoints": random_breakpoints,
    }

def can_place_breakpoint(position: int, existing_breakpoints: List[int], genome_length: int,
                         min_contig_size: int) -> bool:
    if position <= 0 or position >= genome_length:
        return False

    idx = bisect_left(existing_breakpoints, position)
    if idx < len(existing_breakpoints) and existing_breakpoints[idx] == position:
        return False

    left_boundary = 0 if idx == 0 else existing_breakpoints[idx - 1]
    right_boundary = genome_length if idx == len(existing_breakpoints) else existing_breakpoints[idx]
    return (position - left_boundary) >= min_contig_size and (right_boundary - position) >= min_contig_size

def sample_jittered_breakpoint(anchor: int, existing_breakpoints: List[int], genome_length: int,
                               min_contig_size: int, rng: random.Random,
                               breakpoint_jitter: int) -> Optional[int]:
    candidate_positions = [anchor]

    if breakpoint_jitter > 0:
        jitter_std = max(1.0, breakpoint_jitter / 2)
        for _ in range(8):
            jitter = int(round(rng.gauss(0, jitter_std)))
            jitter = max(-breakpoint_jitter, min(breakpoint_jitter, jitter))
            candidate_positions.append(anchor + jitter)

    seen = set()
    for position in candidate_positions:
        position = int(round(position))
        position = max(1, min(genome_length - 1, position))
        if position in seen:
            continue
        seen.add(position)

        if can_place_breakpoint(position, existing_breakpoints, genome_length, min_contig_size):
            return position

    return None

def collect_breakpoint_candidates(prioritized_blocks: List[Tuple[int, int]],
                                  block_categories: Dict[Tuple[int, int], str],
                                  genome_length: int) -> List[Dict[str, object]]:
    candidates = []
    seen_anchors = set()

    for start, end in prioritized_blocks:
        category = block_categories.get((start, end), 'unknown')
        for edge_name, anchor in (("start", start), ("end", end)):
            if anchor <= 0 or anchor >= genome_length or anchor in seen_anchors:
                continue

            seen_anchors.add(anchor)
            candidates.append({
                "anchor": anchor,
                "category": category,
                "type": "targeted",
                "edge": edge_name,
                "source_start": start,
                "source_end": end,
            })

    return candidates

def choose_random_breakpoint(genome_length: int, existing_breakpoints: List[int],
                             min_contig_size: int, rng: random.Random) -> Optional[int]:
    intervals = build_contig_intervals(genome_length, existing_breakpoints)
    eligible_intervals = [
        (start, end, end - start)
        for start, end in intervals
        if (end - start) >= (2 * min_contig_size)
    ]

    if not eligible_intervals:
        return None

    weights = [length * length for _, _, length in eligible_intervals]
    selected_start, selected_end, _ = rng.choices(eligible_intervals, weights=weights, k=1)[0]
    low = selected_start + min_contig_size
    high = selected_end - min_contig_size

    if low > high:
        return None
    if low == high:
        return low

    return rng.randint(low, high)

def is_better_n50_state(candidate_state: Dict[str, object], best_state: Dict[str, object],
                        target_n50: int) -> bool:
    candidate_diff = abs(candidate_state["n50"] - target_n50)
    best_diff = abs(best_state["n50"] - target_n50)

    if candidate_diff != best_diff:
        return candidate_diff < best_diff

    candidate_above_target = candidate_state["n50"] >= target_n50
    best_above_target = best_state["n50"] >= target_n50
    if candidate_above_target != best_above_target:
        return candidate_above_target

    if candidate_state["random_breakpoints"] != best_state["random_breakpoints"]:
        return candidate_state["random_breakpoints"] < best_state["random_breakpoints"]

    return candidate_state["contig_count"] < best_state["contig_count"]

def select_breakpoints_for_target_n50(seq: str, prioritized_blocks: List[Tuple[int, int]],
                                      block_categories: Dict[Tuple[int, int], str],
                                      target_n50: int, seed: int,
                                      mandatory_breakpoints: List[int],
                                      min_contig_size: int = 500,
                                      breakpoint_jitter: int = 100) -> Dict[str, object]:
    genome_length = len(seq)
    rng = random.Random(seed)
    current_breakpoints = sorted({point for point in mandatory_breakpoints if 0 < point < genome_length})
    current_events: List[Dict[str, object]] = []
    current_state = summarize_fragmentation_state(genome_length, current_breakpoints, current_events)
    best_state = summarize_fragmentation_state(genome_length, current_breakpoints, current_events)

    if target_n50 <= 0:
        raise ValueError(f"Target N50 must be a positive integer, got {target_n50}")

    if target_n50 > genome_length:
        debug_print(
            f"Requested target N50 ({target_n50:,} bp) exceeds the total input genome "
            f"length ({genome_length:,} bp); returning the input assembly unchanged"
        )
        return current_state

    if target_n50 > current_state["n50"]:
        raise ValueError(
            f"Target N50 ({target_n50:,} bp) is larger than the input assembly N50 "
            f"({current_state['n50']:,} bp). This mode can fragment, but it cannot scaffold."
        )

    if min_contig_size <= 0:
        raise ValueError(f"Minimum contig size must be positive, got {min_contig_size}")

    debug_print(
        f"Initial N50 fragmentation state: N50={current_state['n50']:,} bp, "
        f"contigs={current_state['contig_count']}"
    )

    candidates = collect_breakpoint_candidates(prioritized_blocks, block_categories, genome_length)
    debug_print(f"Collected {len(candidates)} targeted breakpoint candidates")

    for candidate in candidates:
        if current_state["n50"] <= target_n50:
            break

        position = sample_jittered_breakpoint(
            candidate["anchor"],
            current_state["breakpoints"],
            genome_length,
            min_contig_size,
            rng,
            breakpoint_jitter
        )
        if position is None:
            continue

        current_breakpoints = current_state["breakpoints"] + [position]
        current_events = current_state["events"] + [{
            **candidate,
            "position": position,
            "jitter": position - candidate["anchor"],
        }]
        current_state = summarize_fragmentation_state(genome_length, current_breakpoints, current_events)

        if is_better_n50_state(current_state, best_state, target_n50):
            best_state = summarize_fragmentation_state(genome_length, current_breakpoints, current_events)

    random_attempts = 0
    max_random_attempts = max(1000, genome_length // max(min_contig_size, 1))
    while current_state["n50"] > target_n50 and random_attempts < max_random_attempts:
        random_attempts += 1
        position = choose_random_breakpoint(
            genome_length,
            current_state["breakpoints"],
            min_contig_size,
            rng
        )
        if position is None:
            break

        current_breakpoints = current_state["breakpoints"] + [position]
        current_events = current_state["events"] + [{
            "anchor": position,
            "category": "random",
            "type": "random",
            "edge": "random",
            "source_start": None,
            "source_end": None,
            "position": position,
            "jitter": 0,
        }]
        current_state = summarize_fragmentation_state(genome_length, current_breakpoints, current_events)

        if is_better_n50_state(current_state, best_state, target_n50):
            best_state = summarize_fragmentation_state(genome_length, current_breakpoints, current_events)

    if best_state["n50"] > target_n50:
        debug_print(
            f"WARNING: Could not reach target N50 of {target_n50:,} bp; "
            f"best achievable state was {best_state['n50']:,} bp"
        )
    else:
        debug_print(f"Target N50 reached or crossed; best state N50={best_state['n50']:,} bp")

    return best_state

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

def create_circular_breakpoint_plot(genome_length: int, breakpoint_events: List[Dict[str, object]],
                                    target_n50: int, actual_n50: int, output_pdf: str,
                                    genome_name: str = "Genome"):
    """
    Create a circular genome plot showing breakpoint positions used for N50 fragmentation.
    """
    colors = {
        'gc_only': '#FF6B6B',
        'repeat_only': '#7e57c2',
        'both': '#9ccc65',
        'random': '#FFE66D'
    }

    sectors = {genome_name: genome_length}
    circos = Circos(sectors, space=2)
    sector = circos.sectors[0]
    track = sector.add_track((90, 100))
    track.axis(fc="lightgray", ec="none", lw=0)

    marker_half_width = max(25, min(1000, genome_length // 2000 if genome_length else 25))
    for event in breakpoint_events:
        start = max(0, event["position"] - marker_half_width)
        end = min(genome_length, event["position"] + marker_half_width)
        if start == end:
            end = min(genome_length, start + 1)
        track.rect(start, end, fc=colors[event["category"]], ec="none", lw=0, alpha=0.8)

    major_ticks = []
    tick_interval = 10 ** (len(str(genome_length)) - 1)
    for i in range(0, genome_length, tick_interval):
        major_ticks.append(i)

    track.xticks(
        major_ticks,
        labels=[f"{int(t / 1000)}kb" if t > 0 else "0" for t in major_ticks],
        label_size=8,
        label_orientation="vertical"
    )

    title_text = f"{genome_name}\n"
    title_text += f"Genome size: {genome_length:,} bp | "
    title_text += f"Target N50: {target_n50:,} bp | Actual N50: {actual_n50:,} bp"

    fig = circos.plotfig()
    fig.suptitle(title_text, fontsize=14, y=0.98)

    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor=colors['gc_only'], label="High GC breakpoint"),
        Patch(facecolor=colors['repeat_only'], label="Repeat breakpoint"),
        Patch(facecolor=colors['both'], label="GC+Repeat breakpoint"),
        Patch(facecolor=colors['random'], label="Random breakpoint")
    ]
    fig.legend(handles=legend_elements, loc='lower center', ncol=4,
               bbox_to_anchor=(0.5, -0.05), fontsize=10)

    plt.tight_layout()
    plt.savefig(output_pdf, format='pdf', dpi=300, bbox_inches='tight')
    plt.close()

def create_n50_contig_length_plot(contig_lengths: List[int], target_n50: int, actual_n50: int,
                                  l50: int, output_path: str, genome_name: str = "Genome"):
    """
    Create a ranked contig-length figure for N50 mode.

    The figure contains:
    - A linear-scale zoom over the longest contigs so the N50-defining contig is readable
    - A full ranked overview on a log scale so all contigs remain visible
    """
    if not contig_lengths:
        raise ValueError("Cannot create an N50 plot without contig lengths")

    sorted_lengths = sorted(contig_lengths, reverse=True)
    contig_count = len(sorted_lengths)
    n50_index = max(0, min(l50 - 1, contig_count - 1))

    default_color = "#B8B8B8"
    highlight_color = "#D95F02"
    target_color = "#000000"
    colors = [default_color] * contig_count
    colors[n50_index] = highlight_color

    zoom_count = min(contig_count, max(25, l50 + 15))
    zoom_lengths = sorted_lengths[:zoom_count]
    zoom_colors = colors[:zoom_count]
    zoom_positions = np.arange(zoom_count)

    overview_positions = np.arange(contig_count)
    overview_xmin = max(1, min(sorted_lengths) * 0.5)

    fig = plt.figure(figsize=(13, 10), constrained_layout=True)
    grid = fig.add_gridspec(2, 1, height_ratios=[1.1, 1.6], hspace=0.25)
    ax_zoom = fig.add_subplot(grid[0])
    ax_overview = fig.add_subplot(grid[1])

    ax_zoom.barh(zoom_positions, zoom_lengths, color=zoom_colors, edgecolor="none", height=0.85)
    ax_zoom.invert_yaxis()
    n50_line_width = 0.5
    ax_zoom.axvline(actual_n50, color=highlight_color, linestyle="--", linewidth=n50_line_width, label=f"Actual N50: {actual_n50:,} bp")
    ax_zoom.axvline(target_n50, color=target_color, linestyle=":", linewidth=n50_line_width, label=f"Target N50: {target_n50:,} bp")
    ax_zoom.set_title(f"{genome_name} Contig Lengths", fontsize=14)
    ax_zoom.set_xlabel("Contig length (bp)")
    ax_zoom.set_ylabel("Top contig ranks")
    ax_zoom.grid(axis="x", linestyle=":", linewidth=0.5, alpha=0.6)

    if zoom_count <= 40:
        ax_zoom.set_yticks(zoom_positions)
        ax_zoom.set_yticklabels([str(i + 1) for i in zoom_positions], fontsize=8)
    else:
        tick_positions = np.linspace(0, zoom_count - 1, num=12, dtype=int)
        ax_zoom.set_yticks(tick_positions)
        ax_zoom.set_yticklabels([str(pos + 1) for pos in tick_positions], fontsize=8)

    if n50_index < zoom_count:
        n50_length = sorted_lengths[n50_index]
        ax_zoom.annotate(
            f"N50 contig: rank {n50_index + 1}, {n50_length:,} bp",
            xy=(n50_length, n50_index),
            xytext=(8, 0),
            textcoords="offset points",
            va="center",
            ha="left",
            fontsize=9,
            color=highlight_color,
        )

    line_widths = np.full(contig_count, 0.6)
    line_widths[n50_index] = 2.2
    for y, contig_length, color, width in zip(overview_positions, sorted_lengths, colors, line_widths):
        ax_overview.hlines(y, overview_xmin, contig_length, colors=color, linewidth=width)

    ax_overview.invert_yaxis()
    ax_overview.set_xscale("log")
    ax_overview.axvline(actual_n50, color=highlight_color, linestyle="--", linewidth=n50_line_width)
    ax_overview.axvline(target_n50, color=target_color, linestyle=":", linewidth=n50_line_width)
    ax_overview.set_xlabel("Contig length (bp, log scale)")
    ax_overview.set_ylabel("All contig ranks")
    ax_overview.grid(axis="x", linestyle=":", linewidth=0.5, alpha=0.6)

    overview_tick_positions = np.linspace(0, contig_count - 1, num=min(12, contig_count), dtype=int)
    ax_overview.set_yticks(overview_tick_positions)
    ax_overview.set_yticklabels([str(pos + 1) for pos in overview_tick_positions], fontsize=8)

    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], color=default_color, lw=2, label="Other contigs"),
        Line2D([0], [0], color=highlight_color, lw=2.5, label=f"N50 contig (rank {n50_index + 1})"),
        Line2D([0], [0], color=highlight_color, lw=n50_line_width, linestyle="--", label=f"Actual N50: {actual_n50:,} bp"),
        Line2D([0], [0], color=target_color, lw=n50_line_width, linestyle=":", label=f"Target N50: {target_n50:,} bp"),
    ]
    ax_zoom.legend(handles=legend_elements, loc="lower right")

    fig.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close(fig)

def build_fragment_records_from_intervals(seq: str, intervals: List[Tuple[int, int]],
                                          genome_name: str) -> List[SeqRecord.SeqRecord]:
    frag_records = []
    for i, (start, end) in enumerate(intervals):
        frag_records.append(
            SeqRecord.SeqRecord(
                Seq(seq[start:end]),
                id=f"{genome_name}_frag_{i+1}",
                description=f"coords={start + 1}-{end}"
            )
        )
    return frag_records

def write_n50_fragmentation_report(log_file: str, input_fasta: str, output_fasta: str,
                                   target_n50: int, state: Dict[str, object],
                                   context: Dict[str, object], seed: int,
                                   num_threads: int, max_contig_size: int,
                                   fragment_size: int, overlap: int,
                                   min_contig_size: int, breakpoint_jitter: int,
                                   use_nucmer: bool):
    category_counts = {'gc_only': 0, 'repeat_only': 0, 'both': 0, 'random': 0}
    for event in state["events"]:
        category_counts[event["category"]] += 1

    contig_lengths = sorted(state["lengths"], reverse=True)

    with open(log_file, 'w') as f:
        f.write("=== Draft Genome Fragmentation Report (N50 MODE) ===\n\n")
        f.write(f"Input file: {input_fasta}\n")
        f.write(f"Output file: {output_fasta}\n")
        f.write(f"Target N50: {target_n50:,} bp\n")
        f.write(f"Actual N50: {state['n50']:,} bp\n")
        f.write(f"Input assembly N50: {context['initial_n50']:,} bp\n")
        f.write(f"L50: {state['l50']}\n")
        f.write(f"Random seed: {seed}\n")
        f.write(f"Number of original contigs: {context['num_original_contigs']}\n")
        f.write(f"Output contig count: {state['contig_count']}\n")
        f.write("Completeness preserved: 100.00%\n")
        f.write(f"Red threads used: {num_threads}\n")
        f.write(f"Use nucmer: {use_nucmer}\n")
        f.write(f"Breakpoint jitter (max absolute): {breakpoint_jitter} bp\n")
        f.write(f"Minimum contig size constraint: {min_contig_size} bp\n")
        f.write("Fragmentation settings:\n")
        f.write(f"  - Max contig size: {max_contig_size:,} bp\n")
        f.write(f"  - Fragment size: {fragment_size:,} bp\n")
        f.write(f"  - Overlap: {overlap:,} bp\n\n")

        f.write("=== Repeat Evidence ===\n")
        f.write(f"Red regions: {len(context['red_repeats'])}\n")
        f.write(f"nucmer regions: {len(context['nucmer_repeats'])}\n")
        f.write(f"Combined repeat regions after merging: {len(context['repeats'])}\n\n")

        f.write("=== Region Analysis ===\n")
        f.write(f"GC threshold (mean + 2*std): {context['gc_thresh']:.2%}\n")
        f.write(f"High GC-only regions: {len(context['gc_only'])}\n")
        f.write(f"Repeat-only regions: {len(context['repeat_only'])}\n")
        f.write(f"Both GC+repeat regions: {len(context['both'])}\n")
        f.write(f"Targeted breakpoints used: {len(state['events']) - state['random_breakpoints']}\n")
        f.write(f"Random breakpoints used: {state['random_breakpoints']}\n\n")

        f.write("=== Breakpoint Summary by Category ===\n")
        f.write(f"High GC only breakpoints: {category_counts['gc_only']}\n")
        f.write(f"Repeat only breakpoints: {category_counts['repeat_only']}\n")
        f.write(f"Both GC+repeat breakpoints: {category_counts['both']}\n")
        f.write(f"Random breakpoints: {category_counts['random']}\n\n")

        f.write("=== Output Contig Lengths ===\n")
        for idx, length in enumerate(contig_lengths, 1):
            f.write(f"{idx}\t{length:,}\n")

        f.write("\n=== Detailed Breakpoint List ===\n")
        f.write("Position\tCategory\tType\tAnchor\tJitter\tSourceStart\tSourceEnd\n")
        f.write("-" * 90 + "\n")
        for event in sorted(state["events"], key=lambda item: item["position"]):
            source_start = "" if event["source_start"] is None else f"{event['source_start']:,}"
            source_end = "" if event["source_end"] is None else f"{event['source_end']:,}"
            f.write(
                f"{event['position']:,}\t{event['category']}\t{event['type']}\t"
                f"{event['anchor']:,}\t{event['jitter']:+,}\t{source_start}\t{source_end}\n"
            )

def simulate_draft(input_fasta: str, output_fasta: str, completeness: float, seed: int, num_threads: int = 8,
                  max_contig_size: int = 1_000_000, fragment_size: int = 500_000, overlap: int = 50_000,
                  log_file: str = None, create_visualization: bool = False,
                  stats_json_path: str = None,
                  use_nucmer: bool = False, nucmer_min_identity: float = 95.0,
                  nucmer_min_length: int = 500):
    debug_print("Starting simulate_draft function (FIXED VERSION)")
    context = prepare_simulation_context(
        input_fasta=input_fasta,
        num_threads=num_threads,
        max_contig_size=max_contig_size,
        fragment_size=fragment_size,
        overlap=overlap,
        use_nucmer=use_nucmer,
        nucmer_min_identity=nucmer_min_identity,
        nucmer_min_length=nucmer_min_length
    )
    seq = context["seq"]
    genome_name = context["genome_name"]
    repeats = context["repeats"]
    gc_only = context["gc_only"]
    repeat_only = context["repeat_only"]
    both = context["both"]
    block_categories = context["block_categories"]
    prioritized_blocks = context["prioritized_blocks"]

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
            f.write("=== Draft Genome Simulation Report (FIXED VERSION) ===\n\n")
            f.write(f"Input file: {input_fasta}\n")
            f.write(f"Output file: {output_fasta}\n")
            f.write(f"Target completeness: {completeness:.2%}\n")
            f.write(f"Random seed: {seed}\n")
            f.write(f"Number of original contigs: {context['num_original_contigs']}\n")
            f.write(f"Red threads used: {num_threads}\n")
            f.write(f"Use nucmer: {use_nucmer}\n")
            f.write(f"Fragmentation settings:\n")
            f.write(f"  - Max contig size: {max_contig_size:,} bp\n")
            f.write(f"  - Fragment size: {fragment_size:,} bp\n")
            f.write(f"  - Overlap: {overlap:,} bp\n\n")

            f.write("=== Repeat Evidence ===\n")
            f.write(f"Red regions: {len(context['red_repeats'])}\n")
            f.write(f"nucmer regions: {len(context['nucmer_repeats'])}\n")
            f.write(f"Combined repeat regions after merging: {len(repeats)}\n\n")
            
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
            f.write(f"GC threshold (mean + 2*std): {context['gc_thresh']:.2%}\n")
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

    contig_lengths = [len(frag) for frag in fragments]
    contig_stats = summarize_contig_lengths(contig_lengths)
    stats_json_path = stats_json_path or default_stats_json_path(output_fasta)
    write_contig_stats_json(
        stats_json_path,
        {
            "mode": "completeness",
            "input_fasta": input_fasta,
            "output_fasta": output_fasta,
            "seed": seed,
            "parameters": {
                "completeness": completeness,
                "num_threads": num_threads,
                "max_contig_size": max_contig_size,
                "fragment_size": fragment_size,
                "overlap": overlap,
                "use_nucmer": use_nucmer,
                "nucmer_min_identity": nucmer_min_identity,
                "nucmer_min_length": nucmer_min_length,
            },
            "mode_summary": {
                "target_completeness": completeness,
                "actual_completeness": actual_output_size / removal_stats["total_bases"] if removal_stats["total_bases"] else 0.0,
                "bases_removed": int(removal_stats["total_removed"]),
                "bases_kept": int(actual_output_size),
            },
            "contig_stats": contig_stats,
        }
    )

    debug_print(f"Writing {len(frag_records)} fragments to {output_fasta}")
    SeqIO.write(frag_records, output_fasta, "fasta")
    debug_print(f"Saved contig statistics JSON to {stats_json_path}")
    debug_print("simulate_draft function completed successfully! (FIXED VERSION)")

def simulate_target_n50(input_fasta: str, output_fasta: str, target_n50: int, seed: int,
                        num_threads: int = 8, max_contig_size: int = 1_000_000,
                        fragment_size: int = 500_000, overlap: int = 50_000,
                        log_file: str = None, create_visualization: bool = False,
                        n50_plot_path: str = None, stats_json_path: str = None,
                        min_contig_size: int = 500, breakpoint_jitter: int = 100,
                        use_nucmer: bool = False, nucmer_min_identity: float = 95.0,
                        nucmer_min_length: int = 500):
    debug_print("Starting simulate_target_n50 function")
    context = prepare_simulation_context(
        input_fasta=input_fasta,
        num_threads=num_threads,
        max_contig_size=max_contig_size,
        fragment_size=fragment_size,
        overlap=overlap,
        use_nucmer=use_nucmer,
        nucmer_min_identity=nucmer_min_identity,
        nucmer_min_length=nucmer_min_length
    )

    short_input_records = [layout for layout in context["record_layout"] if layout["length"] < min_contig_size]
    if short_input_records:
        debug_print(
            f"Warning: {len(short_input_records)} input records are shorter than the "
            f"minimum contig size of {min_contig_size} bp and will be preserved as-is"
        )

    fragmentation_state = select_breakpoints_for_target_n50(
        seq=context["seq"],
        prioritized_blocks=context["prioritized_blocks"],
        block_categories=context["block_categories"],
        target_n50=target_n50,
        seed=seed,
        mandatory_breakpoints=context["record_boundaries"],
        min_contig_size=min_contig_size,
        breakpoint_jitter=breakpoint_jitter
    )

    frag_records = build_fragment_records_from_intervals(
        context["seq"],
        fragmentation_state["intervals"],
        context["genome_name"]
    )

    debug_print(
        f"N50 fragmentation completed: target={target_n50:,} bp, "
        f"actual={fragmentation_state['n50']:,} bp, contigs={fragmentation_state['contig_count']}"
    )

    if n50_plot_path:
        debug_print(f"Creating N50 contig-length plot: {n50_plot_path}")
        create_n50_contig_length_plot(
            contig_lengths=fragmentation_state["lengths"],
            target_n50=target_n50,
            actual_n50=fragmentation_state["n50"],
            l50=fragmentation_state["l50"],
            output_path=n50_plot_path,
            genome_name=context["genome_name"]
        )

    if log_file:
        write_n50_fragmentation_report(
            log_file=log_file,
            input_fasta=input_fasta,
            output_fasta=output_fasta,
            target_n50=target_n50,
            state=fragmentation_state,
            context=context,
            seed=seed,
            num_threads=num_threads,
            max_contig_size=max_contig_size,
            fragment_size=fragment_size,
            overlap=overlap,
            min_contig_size=min_contig_size,
            breakpoint_jitter=breakpoint_jitter,
            use_nucmer=use_nucmer
        )

        if create_visualization:
            if log_file.endswith('.txt'):
                pdf_path = log_file.replace('.txt', '_visualization.pdf')
            else:
                pdf_path = log_file + '_visualization.pdf'

            print(f"Creating circular visualization: {pdf_path}")
            try:
                create_circular_breakpoint_plot(
                    genome_length=len(context["seq"]),
                    breakpoint_events=fragmentation_state["events"],
                    target_n50=target_n50,
                    actual_n50=fragmentation_state["n50"],
                    output_pdf=pdf_path,
                    genome_name=context["genome_name"]
                )
                debug_print(f"Visualization saved to: {pdf_path}")
            except Exception as exc:
                print(f"Warning: Failed to create visualization: {exc}")
                debug_print(f"Visualization error: {exc}")

    contig_stats = summarize_contig_lengths(fragmentation_state["lengths"])
    stats_json_path = stats_json_path or default_stats_json_path(output_fasta)
    write_contig_stats_json(
        stats_json_path,
        {
            "mode": "target_n50",
            "input_fasta": input_fasta,
            "output_fasta": output_fasta,
            "seed": seed,
            "parameters": {
                "target_n50": target_n50,
                "num_threads": num_threads,
                "max_contig_size": max_contig_size,
                "fragment_size": fragment_size,
                "overlap": overlap,
                "min_contig_size": min_contig_size,
                "breakpoint_jitter": breakpoint_jitter,
                "use_nucmer": use_nucmer,
                "nucmer_min_identity": nucmer_min_identity,
                "nucmer_min_length": nucmer_min_length,
            },
            "mode_summary": {
                "input_assembly_n50": int(context["initial_n50"]),
                "input_assembly_l50": int(context["initial_l50"]),
                "target_n50": int(target_n50),
                "actual_n50": int(fragmentation_state["n50"]),
                "actual_l50": int(fragmentation_state["l50"]),
                "breakpoint_count": int(len(fragmentation_state["events"])),
                "random_breakpoints": int(fragmentation_state["random_breakpoints"]),
                "completeness_preserved": 1.0,
            },
            "contig_stats": contig_stats,
        }
    )

    debug_print(f"Writing {len(frag_records)} fragments to {output_fasta}")
    SeqIO.write(frag_records, output_fasta, "fasta")
    debug_print(f"Saved contig statistics JSON to {stats_json_path}")
    debug_print("simulate_target_n50 function completed successfully!")

if __name__ == "__main__":
    debug_print("Script started, parsing arguments... (FIXED VERSION)")
    parser = argparse.ArgumentParser(
        description=(
            "Simulate realistic draft bacterial genomes either by reducing completeness "
            "or by fragmenting to a target N50 using GC/repeat-aware breakpoint selection"
        )
    )
    parser.add_argument("--input", required=True, help="Path to complete genome FASTA")
    parser.add_argument("--output", required=True, help="Output FASTA of simulated draft")
    mode_group = parser.add_mutually_exclusive_group(required=True)
    mode_group.add_argument(
        "--completeness",
        type=float,
        help="Target completeness (0-100 for percentage, or 0-1 for fraction)"
    )
    mode_group.add_argument(
        "--target-n50",
        type=int,
        dest="target_n50",
        help="Target output contig N50 in bp while preserving 100%% completeness"
    )
    parser.add_argument("--seed", type=int, default=42, help="Random seed for reproducibility")
    parser.add_argument("--num_threads", type=int, default=8, help="Number of threads for Red (default: 8)")
    parser.add_argument("--max_contig_size", type=int, default=1_000_000, help="Maximum contig size before fragmentation (default: 1MB)")
    parser.add_argument("--fragment_size", type=int, default=500_000, help="Fragment size for large contigs (default: 500kb)")
    parser.add_argument("--overlap", type=int, default=50_000, help="Overlap between fragments (default: 50kb)")
    parser.add_argument("--min-contig-size", type=int, default=500, help="Minimum contig size enforced in --target-n50 mode (default: 500)")
    parser.add_argument("--breakpoint-jitter", type=int, default=100, help="Maximum absolute Gaussian breakpoint jitter in bp for --target-n50 mode (default: 100)")
    parser.add_argument("--use-nucmer", action="store_true", help="Augment Red repeat detection with nucmer self-alignment evidence")
    parser.add_argument("--nucmer-min-identity", type=float, default=95.0, help="Minimum percent identity for nucmer self-homology intervals (default: 95.0)")
    parser.add_argument("--nucmer-min-length", type=int, default=500, help="Minimum alignment length for nucmer self-homology intervals (default: 500)")
    parser.add_argument("--log", help="Path to log file (text report only)")
    parser.add_argument("--vis_log", help="Path to log file with circular PDF visualization")
    parser.add_argument("--n50-plot", dest="n50_plot", help="Path to contig-length plot for --target-n50 mode (PDF/PNG/etc.)")
    parser.add_argument("--stats-json", dest="stats_json", help="Path to JSON file with output contig statistics (default: <output>.stats.json)")
    
    debug_print("Parsing command line arguments...")
    args = parser.parse_args()
    debug_print("Arguments parsed successfully!")
    debug_print(f"Input: {args.input}")
    debug_print(f"Output: {args.output}")  
    if args.completeness is not None:
        debug_print(f"Completeness: {args.completeness}")
    else:
        debug_print(f"Target N50: {args.target_n50}")
    debug_print(f"Threads: {args.num_threads}")
    debug_print(f"Use nucmer: {args.use_nucmer}")
    
    completeness = None
    if args.completeness is not None:
        debug_print("Processing completeness value...")
        completeness = args.completeness
        if completeness > 1:
            completeness = completeness / 100.0
        debug_print(f"Final completeness: {completeness}")
        if completeness < 0 or completeness > 1:
            raise ValueError(
                f"Completeness must be between 0 and 1 (or 0-100 for percentage), got {args.completeness}"
            )
    else:
        if args.target_n50 <= 0:
            raise ValueError(f"Target N50 must be a positive integer, got {args.target_n50}")
        if args.min_contig_size <= 0:
            raise ValueError(f"Minimum contig size must be positive, got {args.min_contig_size}")
        if args.breakpoint_jitter < 0:
            raise ValueError(f"Breakpoint jitter must be non-negative, got {args.breakpoint_jitter}")
        if args.nucmer_min_identity <= 0 or args.nucmer_min_identity > 100:
            raise ValueError(
                f"nucmer minimum identity must be in the range (0, 100], got {args.nucmer_min_identity}"
            )
        if args.nucmer_min_length <= 0:
            raise ValueError(
                f"nucmer minimum alignment length must be positive, got {args.nucmer_min_length}"
            )
    if completeness is not None and args.n50_plot:
        raise ValueError("--n50-plot can only be used together with --target-n50")
    
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
    
    try:
        if completeness is not None:
            debug_print("About to call simulate_draft function... (FIXED VERSION)")
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
                create_visualization=create_vis,
                stats_json_path=args.stats_json,
                use_nucmer=args.use_nucmer,
                nucmer_min_identity=args.nucmer_min_identity,
                nucmer_min_length=args.nucmer_min_length
            )
        else:
            debug_print("About to call simulate_target_n50 function...")
            simulate_target_n50(
                input_fasta=args.input,
                output_fasta=args.output,
                target_n50=args.target_n50,
                seed=args.seed,
                num_threads=num_threads,
                max_contig_size=args.max_contig_size,
                fragment_size=args.fragment_size,
                overlap=args.overlap,
                log_file=log_file,
                create_visualization=create_vis,
                n50_plot_path=args.n50_plot,
                stats_json_path=args.stats_json,
                min_contig_size=args.min_contig_size,
                breakpoint_jitter=args.breakpoint_jitter,
                use_nucmer=args.use_nucmer,
                nucmer_min_identity=args.nucmer_min_identity,
                nucmer_min_length=args.nucmer_min_length
            )
        debug_print("Script completed successfully! (FIXED VERSION)")
    except Exception as e:
        debug_print(f"Script failed with error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
