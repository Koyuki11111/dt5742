#!/bin/bash

# Process both DAQ01 and DAQ02 data directories
# This script runs the full waveform processing pipeline for both daq01 and daq02 folders

set -e  # Exit on error

# Default configuration
BASE_DATA_DIR="/data/test07"
DAQ_NAMES=("daq01" "daq02")
RUN_PARALLEL=false
VERBOSE=false
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
    local config_file="converter_config_${daq_name}.json"

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

    # Create temporary config file for this DAQ
    echo "Creating configuration for $daq_name..."
    python3 - "$BASE_CONFIG" "$config_file" "$daq_dir" "$daq_name" <<'EOF'
import json
import sys

base_config_file = sys.argv[1]
output_config_file = sys.argv[2]
daq_dir = sys.argv[3]
daq_name = sys.argv[4]

# Read base configuration
with open(base_config_file, 'r') as f:
    config = json.load(f)

# Update paths for this DAQ
config['common']['output_dir'] = f"{daq_dir}/output"
config['waveform_converter']['input_dir'] = daq_dir

# Update sensor_ids based on DAQ name
# DAQ01 (daq01): sensors 1,2 (channels 0-7: sensor 1, channels 8-15: sensor 2)
# DAQ02 (daq02): sensors 3,4 (channels 0-7: sensor 3, channels 8-15: sensor 4)
if 'waveform_analyzer' in config and 'sensor_mapping' in config['waveform_analyzer']:
    if 'daq01' in daq_name.lower():
        # DAQ01: sensors 1 and 2
        config['waveform_analyzer']['sensor_mapping']['sensor_ids'] = [1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2]
        print(f"  sensor_ids set to: Channels 0-7 -> Sensor 1, Channels 8-15 -> Sensor 2")
    elif 'daq02' in daq_name.lower():
        # DAQ02: sensors 3 and 4
        config['waveform_analyzer']['sensor_mapping']['sensor_ids'] = [3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4]
        print(f"  sensor_ids set to: Channels 0-7 -> Sensor 3, Channels 8-15 -> Sensor 4")
    else:
        print(f"  WARNING: Unknown DAQ name '{daq_name}', using default sensor_ids")

# Write updated configuration
with open(output_config_file, 'w') as f:
    json.dump(config, f, indent=2)

print(f"Configuration saved to {output_config_file}")
print(f"  input_dir: {config['waveform_converter']['input_dir']}")
print(f"  output_dir: {config['common']['output_dir']}")
EOF

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

    exit 0
else
    echo "WARNING: Some DAQs failed to process:"
    for failed in "${FAILED_DAQS[@]}"; do
        echo "  - $failed"
    done
    echo ""
    exit 1
fi
