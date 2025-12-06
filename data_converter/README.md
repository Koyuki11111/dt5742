# CAEN DT5742 Waveform Processing Pipeline

A three-stage pipeline for converting and analyzing CAEN DT5742 digitizer data.

## Quick Start

```bash
# 1. Setup ROOT environment
# Make sure root-config is in PATH; pick the setup for your system
# macOS (Homebrew): source /opt/homebrew/bin/thisroot.sh
# Ubuntu/Debian:    source /usr/local/bin/thisroot.sh   (or /usr/share/root/bin/thisroot.sh)

# 2. Build
make

# 3. Edit configuration
# Edit converter_config.json with your input paths and settings

# 4. Run
./run_full_pipeline.sh
```

## Pipeline Overview

```
Binary/ASCII → [waveform_converter] → ROOT → [waveform_analyzer] → ROOT + Features → [hdf5_exporter] → HDF5
                convert_to_root         analyze_waveforms               export_to_hdf5
```

**waveform_converter**: Converts binary/ASCII waveform files to ROOT with pedestal subtraction  
**waveform_analyzer**: Extracts timing and amplitude features (CFD, leading edge, charge, etc.)  
**hdf5_exporter**: Exports to HDF5 format for Python/NumPy/pandas analysis

## Configuration

Edit `converter_config.json`:

```json
{
  "common": {
    "output_dir": "output",
    "n_channels": 16,
    "max_cores": 8              // Set > 1 to enable parallel processing
  },
  "waveform_converter": {
    "input_dir": "../AC_LGAD_TEST/",
    "input_pattern": "wave_%d.dat",
    "input_is_ascii": false
  },
  "waveform_analyzer": {
    "signal_region_min": [0.0, 0.0, ...],
    "signal_region_max": [190.0, 190.0, ...]
    // ... per-channel analysis settings
  }
}
```

## Usage

### Run Full Pipeline
```bash
./run_full_pipeline.sh                # All stages
./run_full_pipeline.sh --stage1-only  # Convert only (waveform_converter)
./run_full_pipeline.sh --stage2-only  # Analyze only (waveform_analyzer)
```

### Run Individual Stages
```bash
./convert_to_root                     # waveform_converter
./analyze_waveforms                   # waveform_analyzer
./export_to_hdf5 --mode analysis      # hdf5_exporter
```

## Output Files

All outputs in `output/` directory:

- `output/root/waveforms.root` - Raw waveforms (~46 MB for 523 events)
- `output/root/waveforms_analyzed.root` - Extracted features (~250 KB)
- `output/hdf5/waveforms_analyzed.hdf5` - Features in HDF5 format (~400 KB)

## Python Analysis

```python
import h5py
import pandas as pd

# Load data
with h5py.File('output/hdf5/waveforms_analyzed.hdf5', 'r') as f:
    data = f['AnalysisFeatures'][:]

# Analyze
df = pd.DataFrame(data)
ch0 = df[df['channel'] == 0]
print(f"Mean amplitude: {ch0['ampMax'].mean():.3f} V")
```

## Requirements

- ROOT 6.x (`root-config` available in PATH)
- HDF5 development libraries (pkg-config `hdf5`; e.g., `sudo apt install libhdf5-dev` on Ubuntu or `brew install hdf5` on macOS)
- simdjson development libraries (pkg-config `simdjson`; e.g., `sudo apt install libsimdjson-dev` on Ubuntu or `brew install simdjson` on macOS)
- C++17 compiler

## Troubleshooting

**ROOT not found**: source the appropriate `thisroot.sh` for your install (see Quick Start)
**HDF5 not found**: ensure `pkg-config hdf5` works; otherwise set `HDF5_PREFIX` (e.g., `/usr` on Ubuntu or `/opt/homebrew/opt/hdf5` on macOS)
**simdjson not found**: ensure `pkg-config simdjson` works; otherwise set `JSON_INCLUDES`/`JSON_LIBS` in `Makefile`
**Input files not found**: Check `input_dir` in `converter_config.json`

## Documentation

- **CLAUDE.md** - Detailed documentation (architecture, algorithms, development guide)
- **AGENTS.md** - Agent-based workflow automation
- Each executable has `--help` option

## Project Structure

```
data_converter/
├── include/          # Headers (analysis, config, utils modules)
├── src/              # Implementation files
├── output/           # Generated outputs (root/, hdf5/, plots/)
├── Makefile          # Build system
└── converter_config.json  # Unified configuration
```

Total: ~5,600 lines of modular C++ code across 20 files.
