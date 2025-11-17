# Genome Completeness Simulation

A Python-based tool for simulating realistic draft bacterial genomes by selectively removing regions based on their genomic characteristics. This simulator preferentially removes high GC regions, repetitive sequences, and randomly selected fragments to mimic the natural biases observed in draft genome assemblies.

## Overview

Draft genome assemblies are often incomplete due to various technical limitations in sequencing and assembly processes. This tool simulates incomplete genomes by intelligently removing specific genomic regions that are typically challenging to assemble:

- **High GC content regions**: Difficult to sequence and assemble
- **Repetitive sequences**: Detected using [Red](https://github.com/BioinformaticsToolsmith/Red) (REpeat Detector)
- **Random fragments**: To achieve target completeness levels

The tool generates both a simulated draft genome and detailed visualizations showing which regions were removed and why.

## Features

- **Intelligent region selection**: Prioritizes removal of biologically challenging regions
- **Repeat detection**: Uses Red for accurate identification of repetitive elements
- **Large genome support**: Automatically fragments large contigs for efficient processing
- **Circular visualization**: Generates publication-quality circular genome plots
- **Detailed reporting**: Comprehensive statistics and region-by-region analysis
- **Reproducible**: Seed-based random number generation for consistent results

## Installation

### Install Pixi

```bash
# Linux/macOS
curl -fsSL https://pixi.sh/install.sh | bash

# Windows PowerShell
iwr -useb https://pixi.sh/install.ps1 | iex
```

### Set Up Environment

```bash
# Install all dependencies (Python, BioPython, NumPy, Matplotlib, pyCirclize, Red)
pixi install

# Activate the environment
pixi shell
```

## Usage

### Basic Usage

```bash
# Using pixi (recommended)
pixi run simulate --input GCA_000157015_1.fna \
  --output example_50_output.fna \
  --completeness 0.5 \
  --seed 42

# Or after activating the environment
pixi shell
python simulate_dna_completeness.py \
  --input GCA_000157015_1.fna \
  --output example_50_output.fna \
  --completeness 0.5 \
  --seed 42
```

### With Visualization

```bash
pixi run simulate --input GCA_000157015_1.fna \
  --output example_50_output.fna \
  --completeness 0.5 \
  --vis_log simulation_report.txt \
  --seed 42
```

This example uses the provided `GCA_000157015_1.fna` as input and generates `example_50_output.fna` as the output with 50% target completeness.

### Command-Line Options

| Option | Description | Default |
|--------|-------------|---------|
| `--input` | Path to complete genome FASTA file (required) | - |
| `--output` | Path for output draft genome FASTA (required) | - |
| `--completeness` | Target completeness (0-1 or 0-100 for percentage) | 0.5 |
| `--seed` | Random seed for reproducibility | 42 |
| `--num_threads` | Number of threads for Red | 8 |
| `--max_contig_size` | Maximum contig size before fragmentation (bp) | 1,000,000 |
| `--fragment_size` | Fragment size for large contigs (bp) | 500,000 |
| `--overlap` | Overlap between fragments (bp) | 50,000 |
| `--log` | Path to text log file (report only) | - |
| `--vis_log` | Path to log file with PDF visualization | - |

## Example Output

The tool generates a circular visualization showing the genomic regions that were removed during simulation:

![Example visualization](example_output.png)

**Legend:**
- **Red**: High GC content regions only
- **Purple**: Repetitive regions only
- **Green**: Regions with both high GC and repeats
- **Yellow**: Randomly selected regions

In this example, a 5.5 Mbp genome was reduced to 50% completeness by removing 2.76 Mbp distributed across the genome based on the region characteristics.

## How It Works

1. **GC Content Analysis**: The genome is divided into windows, and regions with unusually high GC content (mean + 2σ) are identified

2. **Repeat Detection**: Red analyzes the genome to identify repetitive sequences
   - Large contigs are automatically fragmented for efficient processing
   - Coordinates are properly mapped back to the original genome

3. **Region Prioritization**: Regions are categorized and prioritized:
   - Highest priority: Regions with both high GC and repeats
   - Medium priority: High GC regions only
   - Medium priority: Repeat regions only
   - Lowest priority: Random regions (to reach target completeness)

4. **Selective Removal**: Regions are removed according to priority until target completeness is achieved

5. **Output Generation**:
   - Draft genome FASTA with remaining fragments
   - Detailed report with statistics
   - Circular visualization (if requested)

## Output Files

### Draft Genome FASTA

Contains the simulated draft genome as multiple contigs (fragments):

```
>genome_frag_1
ATGCGATCG...
>genome_frag_2
GCTAGCTA...
```

### Report File

Detailed statistics including:
- Removal statistics by category
- Region distribution analysis
- List of all removed regions with coordinates
- Actual vs. target completeness

### Visualization PDF

Circular genome plot showing:
- Genome structure
- Removed regions color-coded by category
- Statistics and legend

## Algorithm Details

### Coordinate Mapping

The tool implements proper coordinate mapping when processing fragmented genomes:
- Maintains mapping between fragment coordinates and original genome positions
- Ensures accurate repeat region identification across fragment boundaries
- Properly handles overlapping fragments

### Region Merging

Overlapping regions within each category are merged to prevent:
- Double-counting of removed bases
- Coordinate conflicts
- Inaccurate statistics

## Performance Considerations

- **Small genomes** (< 10 Mbp): Fast processing without fragmentation
- **Large genomes**: Automatic contig fragmentation for efficient Red processing
- **Memory usage**: Proportional to genome size
- **Threading**: Red supports multi-threading for faster repeat detection

## Citation

If you use this tool in your research, please cite:

- **Red**: Girgis, H. Z. (2015). Red: an intelligent, rapid, accurate tool for detecting repeats de-novo on the genomic scale. *BMC Bioinformatics*, 16(1), 227.

## License

This project is provided as-is for academic and research purposes.

## Author

This tool was developed for genome assembly simulation and benchmarking studies.

## Troubleshooting

### Red Installation Issues

If Red is not found, ensure it's installed and in your PATH:
```bash
which Red
```

### Memory Errors

For very large genomes, adjust fragmentation parameters:
```bash
--max_contig_size 500000 --fragment_size 250000
```

### Visualization Errors

Ensure all required Python packages are installed:
```bash
pip install matplotlib pycirclize
```

## Acknowledgments

This tool uses the Red repeat detector developed by the Bioinformatics Toolsmith team.
