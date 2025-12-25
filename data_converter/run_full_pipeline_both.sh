#!/bin/bash

# Process both DAQ01 and DAQ02 data directories
# This script runs the full waveform processing pipeline for both daq01 and daq02 folders

set -e  # Exit on error

# Default configuration
BASE_DATA_DIR="/home/blim/epic/data/000004"
DAQ_NAMES=("daq00" "daq01")
RUN_PARALLEL=false
VERBOSE=false
WITH_QA_COMPARISON=false
NUM_QA_EVENTS=5
POSITIONAL_ARGS=()

# Resolve script directory to locate helper scripts even when invoked elsewhere
SCRIPT_DIR="$(cd -- "$(dirname "$0")" && pwd)"
BASE_CONFIG="${SCRIPT_DIR}/converter_config.json"
RUN_FULL_PIPELINE="${SCRIPT_DIR}/run_full_pipeline.sh"

# Parse command line arguments
print_usage() {
    cat << EOF
Process both DAQ01 and DAQ02 data directories

Usage: $0 [options] [RUN_ID]

Options:
    --base-config FILE       Base configuration file (default: converter_config.json)
    --base-dir DIR           Base data directory containing daq01 and daq02 (default: /data/test07)
    --daq-names NAME1,NAME2  Comma-separated DAQ directory names (default: daq01,daq02)
    --parallel               Run both DAQs in parallel (instead of sequential)
    --with-qa-comparison     Run QA comparison after processing (requires merged HDF5 directory)
    --num-qa-events N        Number of events for QA comparison (default: 5)
    --verbose                Verbose output
    -h, --help               Show this help message

Positional:
    RUN_ID                  Optional; overrides base dir to /data/RUN_ID
                           If RUN_ID is all digits, it is zero-padded to 6 digits

Examples:
    # Process both daq01 and daq02 sequentially
    $0

    # Process both DAQs in parallel
    $0 --parallel

    # Override base dir via RUN_ID (uses /data/test08)
    $0 test08

    # Override base dir via numeric RUN_ID (uses /data/000123)
    $0 123

    # Custom base directory
    $0 --base-dir /data/test08

    # Custom DAQ names
    $0 --daq-names daqA,daqB

Output Structure:
    Each DAQ will have its own output directory:
      /data/test07/daq01/output/
      /data/test07/daq02/output/

EOF
}

while [[ $# -gt 0 ]]; do
    case $1 in
        --base-config)
            BASE_CONFIG="$2"
            shift 2
            ;;
        --base-dir)
            BASE_DATA_DIR="$2"
            shift 2
            ;;
        --daq-names)
            IFS=',' read -ra DAQ_NAMES <<< "$2"
            shift 2
            ;;
        --parallel)
            RUN_PARALLEL=true
            shift
            ;;
        --with-qa-comparison)
            WITH_QA_COMPARISON=true
            shift
            ;;
        --num-qa-events)
            NUM_QA_EVENTS="$2"
            shift 2
            ;;
        --verbose)
            VERBOSE=true
            shift
            ;;
        -h|--help)
            print_usage
            exit 0
            ;;
        *)
            if [[ $1 == -* ]]; then
                echo "ERROR: Unknown option $1"
                print_usage
                exit 1
            fi
            POSITIONAL_ARGS+=("$1")
            shift
            ;;
    esac
done

# Handle positional base-dir override
if [ ${#POSITIONAL_ARGS[@]} -gt 1 ]; then
    echo "ERROR: Too many positional arguments. Provide at most one RUN_ID."
    print_usage
    exit 1
fi

if [ ${#POSITIONAL_ARGS[@]} -eq 1 ]; then
    RUN_ID="${POSITIONAL_ARGS[0]}"
    if [[ $RUN_ID =~ ^[0-9]+$ ]]; then
        printf -v RUN_ID "%06d" "$RUN_ID"
    fi
    BASE_DATA_DIR="/data/${RUN_ID}"
fi

# Check if base config exists
if [ ! -f "$BASE_CONFIG" ]; then
    echo "ERROR: Base configuration file '$BASE_CONFIG' not found"
    exit 1
fi

# Check if Python3 is available
if ! command -v python3 &> /dev/null; then
    echo "ERROR: python3 not found. Required for JSON config manipulation"
    exit 1
fi

# Check if run_full_pipeline.sh exists
if [ ! -f "$RUN_FULL_PIPELINE" ]; then
    echo "ERROR: run_full_pipeline.sh not found at '$RUN_FULL_PIPELINE'"
    exit 1
fi

echo "=========================================="
echo "Dual DAQ Processing Pipeline"
echo "=========================================="
echo "Base directory: $BASE_DATA_DIR"
echo "DAQ names: ${DAQ_NAMES[@]}"
echo "Base config: $BASE_CONFIG"
echo "Processing mode: $([ "$RUN_PARALLEL" = true ] && echo "PARALLEL" || echo "SEQUENTIAL")"
echo ""

# Function to process a single DAQ
process_daq() {
    local daq_name=$1
    local daq_dir="${BASE_DATA_DIR}/${daq_name}"
    local config_file="converter_config_${daq_name}_temp.json"

    echo "=========================================="
    echo "Processing: $daq_name"
    echo "=========================================="
    echo "DAQ directory: $daq_dir"
    echo "Config file: $config_file"
    echo ""

    # Check if DAQ directory exists
    if [ ! -d "$daq_dir" ]; then
        echo "WARNING: DAQ directory '$daq_dir' not found. Skipping."
        return 1
    fi

    # Check if fixed config file exists, otherwise create from base config
    local fixed_config="${SCRIPT_DIR}/converter_config_${daq_name}.json"

    if [ -f "$fixed_config" ]; then
        echo "Using existing config: $fixed_config"
        python3 - "$fixed_config" "$config_file" "$daq_dir" "$daq_name" <<'EOF'
import json
import sys
import os

fixed_config_file = sys.argv[1]
output_config_file = sys.argv[2]
daq_dir = sys.argv[3]
daq_name = sys.argv[4]

# Read fixed configuration
with open(fixed_config_file, 'r') as f:
    config = json.load(f)

# Extract runnumber from daq_dir path
# daq_dir format: /home/blim/epic/data/000004/daq00
parts = daq_dir.rstrip('/').split('/')
run_str = parts[-2] if len(parts) >= 2 else "000000"
runnumber = int(run_str)
base_data_dir = '/'.join(parts[:-2])

# Update only runnumber (keep other settings from fixed config)
config['common']['runnumber'] = runnumber

# Ensure paths are set correctly
if 'output_dir' not in config['common'] or not config['common']['output_dir']:
    config['common']['output_dir'] = base_data_dir
if 'daq_name' not in config['common'] or not config['common']['daq_name']:
    config['common']['daq_name'] = daq_name
if 'input_dir' not in config['waveform_converter'] or not config['waveform_converter']['input_dir']:
    config['waveform_converter']['input_dir'] = base_data_dir

# Write configuration for this run
with open(output_config_file, 'w') as f:
    json.dump(config, f, indent=2)

print(f"  Loaded from: {os.path.basename(fixed_config_file)}")
print(f"  runnumber: {config['common']['runnumber']}")
print(f"  sensor_ids: {config.get('waveform_analyzer', {}).get('sensor_mapping', {}).get('sensor_ids', 'N/A')[:4]}...")
EOF
    else
        echo "Creating new config from base: $BASE_CONFIG"
        python3 - "$BASE_CONFIG" "$config_file" "$daq_dir" "$daq_name" "$fixed_config" <<'EOF'
import json
import sys

base_config_file = sys.argv[1]
output_config_file = sys.argv[2]
daq_dir = sys.argv[3]
daq_name = sys.argv[4]
fixed_config_file = sys.argv[5]

# Read base configuration
with open(base_config_file, 'r') as f:
    config = json.load(f)

# Extract runnumber from daq_dir path
# daq_dir format: /home/blim/epic/data/000004/daq00
parts = daq_dir.rstrip('/').split('/')
run_str = parts[-2] if len(parts) >= 2 else "000000"
runnumber = int(run_str)
base_data_dir = '/'.join(parts[:-2])

# Update paths using runnumber/daq_name pattern
config['common']['output_dir'] = base_data_dir
config['common']['runnumber'] = runnumber
config['common']['daq_name'] = daq_name
config['waveform_converter']['input_dir'] = base_data_dir

# Update sensor mapping based on DAQ name
# DAQ00: sensor 0, columns 0-1, rows 0-7
# DAQ01: sensor 1, columns 0-1, rows 0-7
if 'waveform_analyzer' in config and 'sensor_mapping' in config['waveform_analyzer']:
    if 'daq00' in daq_name.lower():
        config['waveform_analyzer']['sensor_mapping']['sensor_ids'] = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        config['waveform_analyzer']['sensor_mapping']['column_ids'] = [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1]
        config['waveform_analyzer']['sensor_mapping']['strip_ids'] = [0, 1, 2, 3, 4, 5, 6, 7, 0, 1, 2, 3, 4, 5, 6, 7]
        print(f"  Sensor 0: Channels 0-7 -> Col 0 (rows 0-7), Channels 8-15 -> Col 1 (rows 0-7)")
    elif 'daq01' in daq_name.lower():
        config['waveform_analyzer']['sensor_mapping']['sensor_ids'] = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
        config['waveform_analyzer']['sensor_mapping']['column_ids'] = [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1]
        config['waveform_analyzer']['sensor_mapping']['strip_ids'] = [0, 1, 2, 3, 4, 5, 6, 7, 0, 1, 2, 3, 4, 5, 6, 7]
        print(f"  Sensor 1: Channels 0-7 -> Col 0 (rows 0-7), Channels 8-15 -> Col 1 (rows 0-7)")
    else:
        print(f"  WARNING: Unknown DAQ name '{daq_name}', using default sensor_ids")

# Save as fixed config file for future use
with open(fixed_config_file, 'w') as f:
    json.dump(config, f, indent=2)
print(f"  Created fixed config: {fixed_config_file}")

# Also write as temporary config for this run
with open(output_config_file, 'w') as f:
    json.dump(config, f, indent=2)

print(f"  runnumber: {runnumber}")
print(f"  -> Will create: {base_data_dir}/{run_str}/{daq_name}/output/")
EOF
    fi

    if [ $? -ne 0 ]; then
        echo "ERROR: Failed to create config for $daq_name"
        return 1
    fi

    echo ""

    # Run the pipeline for this DAQ
    echo "Running pipeline for $daq_name..."
    if [ "$VERBOSE" = true ]; then
        "$RUN_FULL_PIPELINE" --config "$config_file" --verbose
    else
        "$RUN_FULL_PIPELINE" --config "$config_file"
    fi

    if [ $? -eq 0 ]; then
        echo ""
        echo "SUCCESS: $daq_name processing completed"
        echo ""

        # Clean up temporary config file
        rm -f "$config_file"
        return 0
    else
        echo ""
        echo "ERROR: $daq_name processing failed"
        echo ""
        return 1
    fi
}

# Process DAQs
FAILED_DAQS=()

if [ "$RUN_PARALLEL" = true ]; then
    echo "Starting parallel processing of all DAQs..."
    echo ""

    # Array to store background PIDs
    PIDS=()

    # Launch all DAQs in background
    for daq_name in "${DAQ_NAMES[@]}"; do
        process_daq "$daq_name" &
        PIDS+=($!)
    done

    # Wait for all background jobs
    for i in "${!PIDS[@]}"; do
        daq_name="${DAQ_NAMES[$i]}"
        pid="${PIDS[$i]}"

        if wait "$pid"; then
            echo "PARALLEL: $daq_name completed successfully"
        else
            echo "PARALLEL: $daq_name failed"
            FAILED_DAQS+=("$daq_name")
        fi
    done
else
    echo "Starting sequential processing of DAQs..."
    echo ""

    # Process DAQs sequentially
    for daq_name in "${DAQ_NAMES[@]}"; do
        if ! process_daq "$daq_name"; then
            FAILED_DAQS+=("$daq_name")
        fi
    done
fi

echo "=========================================="
echo "Dual DAQ Processing Summary"
echo "=========================================="
echo ""

if [ ${#FAILED_DAQS[@]} -eq 0 ]; then
    echo "SUCCESS: All DAQs processed successfully!"
    echo ""

    for daq_name in "${DAQ_NAMES[@]}"; do
        daq_dir="${BASE_DATA_DIR}/${daq_name}"
        echo "$daq_name output:"
        echo "  ${daq_dir}/output/"

        # List output files if they exist
        if [ -d "${daq_dir}/output/root" ]; then
            echo "  ROOT files:"
            ls -lh "${daq_dir}/output/root"/*.root 2>/dev/null | awk '{print "    " $9 " (" $5 ")"}'
        fi

        if [ -d "${daq_dir}/output/hdf5" ]; then
            echo "  HDF5 files:"
            ls -lh "${daq_dir}/output/hdf5"/*.h5 2>/dev/null | awk '{print "    " $9 " (" $5 ")"}'
        fi
        echo ""
    done

    # Merge quality check files into shared QA directory
    echo "=========================================="
    echo "Merging Quality Check Files"
    echo "=========================================="

    SHARED_QA_DIR="${BASE_DATA_DIR}/qa"
    mkdir -p "$SHARED_QA_DIR"

    # Collect quality check files from all DAQs
    QC_FILES=()
    for daq_name in "${DAQ_NAMES[@]}"; do
        daq_qc_file="${BASE_DATA_DIR}/${daq_name}/output/quality_check/quality_check.root"
        if [ -f "$daq_qc_file" ]; then
            QC_FILES+=("$daq_qc_file")
            echo "Found: $daq_qc_file"
        else
            echo "WARNING: No quality check file for $daq_name"
        fi
    done

    # Merge if we have quality check files
    if [ ${#QC_FILES[@]} -gt 0 ]; then
        MERGED_QC_FILE="${SHARED_QA_DIR}/quality_check.root"

        echo ""
        echo "Merging ${#QC_FILES[@]} file(s) into: $MERGED_QC_FILE"

        hadd -f "$MERGED_QC_FILE" "${QC_FILES[@]}" > /dev/null 2>&1

        if [ $? -eq 0 ]; then
            echo "SUCCESS: Quality check files merged"

            MERGED_SIZE=$(ls -lh "$MERGED_QC_FILE" | awk '{print $5}')
            echo "  Merged file size: $MERGED_SIZE"
            echo ""

            # Delete individual DAQ quality check files
            echo "Cleaning up individual DAQ files..."
            for qc_file in "${QC_FILES[@]}"; do
                rm -f "$qc_file"
                echo "  Deleted: $qc_file"
            done
            echo ""

            echo "All 4 sensors available in: $MERGED_QC_FILE"
        else
            echo "ERROR: Failed to merge quality check files"
            echo "  Individual files preserved"
        fi
    else
        echo "WARNING: No quality check files found"
    fi

    echo ""
    echo "=========================================="
    echo ""

    # Run QA comparison if requested
    if [ "$WITH_QA_COMPARISON" = true ]; then
        echo "=========================================="
        echo "Running QA Comparison"
        echo "=========================================="

        # Extract run number from BASE_DATA_DIR
        RUN_ID=$(basename "$BASE_DATA_DIR")

        # Check if merged HDF5 directory exists
        MERGED_HDF5_DIR="${BASE_DATA_DIR}/merged/hdf5"
        if [ ! -d "$MERGED_HDF5_DIR" ]; then
            echo "WARNING: Merged HDF5 directory not found: $MERGED_HDF5_DIR"
            echo "Skipping QA comparison. Please merge HDF5 files first."
            echo ""
        else
            # Check if quality check files exist
            DAQ00_QC="${BASE_DATA_DIR}/daq00/output/quality_check/quality_check.root"
            DAQ01_QC="${BASE_DATA_DIR}/daq01/output/quality_check/quality_check.root"

            if [ ! -f "$DAQ00_QC" ] || [ ! -f "$DAQ01_QC" ]; then
                echo "WARNING: Quality check files not found for both DAQs"
                echo "  DAQ00: $DAQ00_QC"
                echo "  DAQ01: $DAQ01_QC"
                echo "Skipping QA comparison."
                echo ""
            else
                # Extract base directory (parent of run directory)
                BASE_DIR=$(dirname "$BASE_DATA_DIR")

                echo "Running QA comparison for run $RUN_ID..."
                echo "  Base dir: $BASE_DIR"
                echo "  Events:   $NUM_QA_EVENTS"
                echo ""

                # Run qa_comparison
                if "${SCRIPT_DIR}/src/qa_comparison" "$RUN_ID" --base-dir "$BASE_DIR" --num-events "$NUM_QA_EVENTS"; then
                    echo ""
                    echo "SUCCESS: QA comparison completed"
                    QA_OUTPUT="${BASE_DATA_DIR}/merged/qa"
                    echo "  Output: $QA_OUTPUT"
                    echo ""
                else
                    echo ""
                    echo "WARNING: QA comparison failed"
                    echo ""
                fi
            fi
        fi

        echo "=========================================="
        echo ""
    fi

    exit 0
else
    echo "WARNING: Some DAQs failed to process:"
    for failed in "${FAILED_DAQS[@]}"; do
        echo "  - $failed"
    done
    echo ""
    exit 1
fi
