#!/bin/bash

# Full waveform processing pipeline
# Stage 1 (waveform_converter): Convert binary/ASCII to ROOT
# Stage 2 (waveform_analyzer): Analyze waveforms
# Stage 3 (hdf5_exporter): Export to HDF5

set -e  # Exit on error

# Resolve script directory so we can find configs/binaries when invoked elsewhere
SCRIPT_DIR="$(cd -- "$(dirname "$0")" && pwd)"

# Default configuration file (unified)
PIPELINE_CONFIG="${SCRIPT_DIR}/converter_config.json"

# Output directory will be loaded from config (common.output_dir)
OUTPUT_DIR="output"

# Default intermediate/output files (relative names). For a minimal smoke test,
# drop a few tiny mixed-length inputs under ./temp/smoke/ (see README) and run:
#   ./run_full_pipeline.sh --config converter_config.json
RAW_ROOT="waveforms.root"
ANALYSIS_ROOT="waveforms_analyzed.root"
ANALYSIS_HDF5="waveforms_analyzed.h5"

# Parse command line arguments
RUN_STAGE1=true
RUN_STAGE2=true
RUN_STAGE3=true
VERBOSE=false
USE_PARALLEL=false

print_usage() {
    cat << EOF
Full waveform processing pipeline

Usage: $0 [options]

Options:
    --config FILE            Pipeline configuration (default: converter_config.json)
    --raw-root FILE          Raw ROOT filename (default: waveforms.root)
    --analysis-root FILE     Analysis ROOT filename (default: waveforms_analyzed.root)
    --analysis-hdf5 FILE     Analysis HDF5 filename (default: waveforms_analyzed.h5)
    --stage1-only            Run only stage 1 (waveform_converter)
    --stage2-only            Run only stage 2 (waveform_analyzer)
    --stage3-only            Run only stage 3 (hdf5_exporter)
    --skip-stage1            Skip stage 1 (waveform_converter)
    --skip-stage2            Skip stage 2 (waveform_analyzer)
    --skip-stage3            Skip stage 3 (hdf5_exporter)
    --parallel               Force parallel processing (overrides config auto-detection)
    --verbose                Verbose output
    -h, --help               Show this help message

Parallel Processing:
    By default, the pipeline auto-detects parallel mode from config files:
      - analysis_config.json: parallel_max_cores > 1 enables parallel Stage 2
      - wave_converter_config.json: max_threads > 1 enables parallel Stage 1
    Use --parallel to force parallel mode regardless of config settings.

Stages:
    Stage 1: waveform_converter (convert_to_root)   - Convert binary/ASCII waveforms to ROOT format
    Stage 2: waveform_analyzer (analyze_waveforms)  - Extract timing and amplitude features
    Stage 3: hdf5_exporter (export_to_hdf5)         - Export analyzed ROOT data to HDF5 format

Output Organization:
    All outputs are organized in subdirectories (default: output/):
      output/root/            - ROOT files (waveforms.root, waveforms_analyzed.root)
      output/hdf5/            - HDF5 files (waveforms_analyzed.h5 and sensor splits)
      output/waveform_plots/  - Waveform plots debug output (if enabled)
      output/temp/            - Temporary files for parallel processing

    To change the output directory, edit the "output_dir" field in:
      - wave_converter_config.json
      - analysis_config.json

Examples:
    # Run full pipeline (auto-detects parallel mode from config files)
    $0

    # Force parallel mode regardless of config settings
    $0 --parallel

    # Run only conversion step
    $0 --stage1-only

    # Skip analysis, only convert (no analysis HDF5 will be produced)
    $0 --skip-stage2

    # Custom configuration files
    $0 --converter-config my_config.json --analysis-config my_analysis.json

Note:
    To control parallel processing, edit the config files:
      - analysis_config.json: set "parallel_max_cores" (e.g., 8 for 8 cores, 1 for sequential)
      - wave_converter_config.json: set "max_threads" (e.g., 16 for 16 threads, 1 for sequential)

EOF
}

while [[ $# -gt 0 ]]; do
    case $1 in
        --config)
            PIPELINE_CONFIG="$2"
            shift 2
            ;;
        --raw-root)
            RAW_ROOT="$2"
            shift 2
            ;;
        --analysis-root)
            ANALYSIS_ROOT="$2"
            shift 2
            ;;
        --analysis-hdf5)
            ANALYSIS_HDF5="$2"
            shift 2
            ;;
        --stage1-only)
            RUN_STAGE2=false
            RUN_STAGE3=false
            shift
            ;;
        --stage2-only)
            RUN_STAGE1=false
            RUN_STAGE3=false
            shift
            ;;
        --stage3-only)
            RUN_STAGE1=false
            RUN_STAGE2=false
            shift
            ;;
        --skip-stage1)
            RUN_STAGE1=false
            shift
            ;;
        --skip-stage2)
            RUN_STAGE2=false
            shift
            ;;
        --skip-stage3)
            RUN_STAGE3=false
            shift
            ;;
        --parallel)
            USE_PARALLEL=false
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
            echo "ERROR: Unknown option $1"
            print_usage
            exit 1
            ;;
    esac
done

# Normalize PIPELINE_CONFIG to an absolute path.
# Prefer the user’s working directory; if not found, fall back to the script directory.
if [[ "$PIPELINE_CONFIG" != /* ]]; then
    CWD_CONFIG="$(cd -- "$(dirname "$PIPELINE_CONFIG")" && pwd)/$(basename "$PIPELINE_CONFIG")"
    SCRIPT_CONFIG="$(cd -- "$SCRIPT_DIR" && cd -- "$(dirname "$PIPELINE_CONFIG")" && pwd)/$(basename "$PIPELINE_CONFIG")"

    if [ -f "$CWD_CONFIG" ]; then
        PIPELINE_CONFIG="$CWD_CONFIG"
    elif [ -f "$SCRIPT_CONFIG" ]; then
        PIPELINE_CONFIG="$SCRIPT_CONFIG"
    else
        PIPELINE_CONFIG="$CWD_CONFIG"
    fi
fi

# Ensure the config file exists before continuing
if [ ! -f "$PIPELINE_CONFIG" ]; then
    echo "ERROR: Pipeline configuration file '$PIPELINE_CONFIG' not found"
    exit 1
fi

# Require python3 for config parsing
if ! command -v python3 &> /dev/null; then
    echo "ERROR: python3 not found. Required for config parsing"
    exit 1
fi

# Load output_dir from config if present
if [ -f "$PIPELINE_CONFIG" ]; then
    RAW_OUTPUT_DIR=$(python3 - "$PIPELINE_CONFIG" <<'EOF'
import json, sys
cfg_path = sys.argv[1]
try:
    with open(cfg_path) as f:
        cfg = json.load(f)
    print(cfg.get("common", {}).get("output_dir", "output"))
except Exception:
    print("output")
EOF
)
    CONFIG_DIR="$(cd -- "$(dirname "$PIPELINE_CONFIG")" && pwd)"
    # If the config's output_dir is relative, resolve it relative to the config file directory
    case "$RAW_OUTPUT_DIR" in
        /*) OUTPUT_DIR="$RAW_OUTPUT_DIR" ;;
        *) OUTPUT_DIR="${CONFIG_DIR}/${RAW_OUTPUT_DIR}" ;;
    esac
fi


daq_name=$(sed -nE 's/^[[:space:]]*"daq_name"[[:space:]]*:[[:space:]]*"([^"]+)".*/\1/p' "$PIPELINE_CONFIG" | head -n1)
runnumber=$(sed -nE 's/^[[:space:]]*"runnumber"[[:space:]]*:[[:space:]]*([0-9]+).*/\1/p' "$PIPELINE_CONFIG" | head -n1)
output_dir=$(sed -nE 's/^[[:space:]]*"output_dir"[[:space:]]*:[[:space:]]*"([^"]+)".*/\1/p' "$PIPELINE_CONFIG" | head -n1)

run_str=$(printf "%06d" "$runnumber")

BASE_DIR="${output_dir}/${run_str}/${daq_name}"
OUTPUT_DIR=$BASE_DIR/output

# Check ROOT environment
if ! command -v root-config &> /dev/null; then
    echo "ERROR: ROOT not found. Please source ROOT environment:"
    echo "  source /opt/root6/root_install/bin/thisroot.sh"
    exit 1
fi

# Auto-detect parallel mode from config file (unless explicitly set by --parallel flag)
if [ "$USE_PARALLEL" = false ]; then
    if [ -f "$PIPELINE_CONFIG" ]; then
        MAX_CORES=$(python3 -c "
import json, sys
try:
    with open('$PIPELINE_CONFIG') as f:
        config = json.load(f)
        common = config.get('common', {})
        print(common.get('max_cores', 1))
except:
    print(1)
" 2>/dev/null)

        if [ "$MAX_CORES" -gt 1 ]; then
            USE_PARALLEL=true
            echo "INFO: Auto-detected parallel mode from $PIPELINE_CONFIG (max_cores=$MAX_CORES)"
        fi
    fi
fi

echo "=========================================="
echo "Waveform Processing Pipeline"
echo "=========================================="
echo ""

# Stage 1: Convert to ROOT
if [ "$RUN_STAGE1" = true ]; then
    echo "Stage 1: Converting waveforms to ROOT format..."
    echo "  Config: $PIPELINE_CONFIG"
    echo "  Output: $OUTPUT_DIR/root/$RAW_ROOT"
    
    if [ "$USE_PARALLEL" = true ]; then
        echo "  Mode:   PARALLEL"
    else
        echo "  Mode:   SEQUENTIAL"
    fi
    echo ""

    if [ ! -f "${SCRIPT_DIR}/convert_to_root" ]; then
        echo "ERROR: convert_to_root executable not found at ${SCRIPT_DIR}/convert_to_root"
        echo "Please run 'make' to build the executables"
        exit 1
    fi

    # Note: convert_to_root reads output_dir from config and automatically creates output_dir/root/
    if [ "$USE_PARALLEL" = true ]; then
        "${SCRIPT_DIR}/convert_to_root" --config "$PIPELINE_CONFIG" --root "$RAW_ROOT" --parallel
    else
        "${SCRIPT_DIR}/convert_to_root" --config "$PIPELINE_CONFIG" --root "$RAW_ROOT"
    fi

    if [ $? -ne 0 ]; then
        echo "ERROR: Stage 1 failed"
        exit 1
    fi
    echo ""
fi

# Stage 2: Analyze waveforms
if [ "$RUN_STAGE2" = true ]; then
    echo "Stage 2: Analyzing waveforms..."
    echo "  Config: $PIPELINE_CONFIG"
    echo "  Input:  $OUTPUT_DIR/root/$RAW_ROOT"
    echo "  Output: $OUTPUT_DIR/root/$ANALYSIS_ROOT"
    
    # Check for input file (it should be in output_dir/root/)
    if [ ! -f "$OUTPUT_DIR/root/$RAW_ROOT" ]; then
        echo "ERROR: Input file $OUTPUT_DIR/root/$RAW_ROOT not found"
        echo "Please run stage 1 first or specify correct input file"
        exit 1
    fi

    if [ "$USE_PARALLEL" = true ]; then
        echo "  Mode:   PARALLEL"
        echo ""

        if [ ! -x "${SCRIPT_DIR}/parallel_analyze.sh" ]; then
            echo "ERROR: parallel_analyze.sh not found or not executable at ${SCRIPT_DIR}/parallel_analyze.sh"
            echo "Please run 'make' to build the executables"
            exit 1
        fi

        "${SCRIPT_DIR}/parallel_analyze.sh" --config "$PIPELINE_CONFIG" --input "$RAW_ROOT" --output "$ANALYSIS_ROOT"

	if [ $? -ne 0 ]; then
            echo "ERROR: Stage 2 (parallel) failed"
            exit 1
        fi
    else
        echo "  Mode:   SEQUENTIAL"
        echo ""

        if [ ! -f "${SCRIPT_DIR}/analyze_waveforms" ]; then
            echo "ERROR: analyze_waveforms executable not found at ${SCRIPT_DIR}/analyze_waveforms"
            echo "Please run 'make' to build the executables"
            exit 1
        fi

        "${SCRIPT_DIR}/analyze_waveforms" --config "$PIPELINE_CONFIG" --input "$RAW_ROOT" --output "$ANALYSIS_ROOT"

        if [ $? -ne 0 ]; then
            echo "ERROR: Stage 2 failed"
            exit 1
        fi
    fi
    echo ""

    # Stage 2.5: Generate quality check plots
    echo "Stage 2.5: Generating quality check plots..."
    echo "  Config: $PIPELINE_CONFIG"
    echo "  Input:  $OUTPUT_DIR/root/$ANALYSIS_ROOT"
    echo "  Output: $OUTPUT_DIR/quality_check/quality_check.root"
    echo ""

    if [ ! -f "${SCRIPT_DIR}/fast_qa" ]; then
        echo "WARNING: fast_qa executable not found at ${SCRIPT_DIR}/fast_qa"
        echo "         Skipping quality check generation..."
    else
        "${SCRIPT_DIR}/fast_qa" --config "$PIPELINE_CONFIG"

        if [ $? -ne 0 ]; then
            echo "WARNING: Fast QA failed (continuing anyway)"
        fi
    fi
    echo ""
fi

# Stage 3: Export to HDF5
if [ "$RUN_STAGE3" = true ]; then
    echo "Stage 3: Exporting analyzed data to HDF5 format..."

    if [ ! -f "${SCRIPT_DIR}/export_to_hdf5" ]; then
        echo "ERROR: export_to_hdf5 executable not found at ${SCRIPT_DIR}/export_to_hdf5"
        echo "Please run 'make' to build the executables"
        exit 1
    fi

    # Extract unique sensor IDs from pipeline config
    if [ -f "$PIPELINE_CONFIG" ]; then
        SENSOR_IDS=$(python3 -c "
import json, sys
try:
    with open('$PIPELINE_CONFIG') as f:
        config = json.load(f)
        stage2 = config.get('waveform_analyzer', {})
        sensor_mapping = stage2.get('sensor_mapping', {})
        sensor_ids = sensor_mapping.get('sensor_ids', [])
        unique_ids = sorted(set(sensor_ids))
        print(' '.join(map(str, unique_ids)))
except:
    print('')
" 2>/dev/null)

        if [ -z "$SENSOR_IDS" ]; then
            echo "WARNING: Could not extract sensor IDs from $PIPELINE_CONFIG"
            echo "         Exporting all channels to single file"
            SENSOR_IDS=""
        else
            echo "  Found sensors: $SENSOR_IDS"
        fi
    else
        SENSOR_IDS=""
    fi

    # Export to Corryvreckan HDF5 format if analysis exists
    if [ -f "$OUTPUT_DIR/root/$ANALYSIS_ROOT" ] && [ "$RUN_STAGE2" = true ] || [ "$RUN_STAGE2" = false ]; then
        if [ -n "$SENSOR_IDS" ]; then
            # Export per sensor in Corryvreckan format
            for SENSOR_ID in $SENSOR_IDS; do
                OUTPUT_FILE="waveforms_corry_sensor$(printf '%02d' $SENSOR_ID).h5"
                echo "  Exporting Corryvreckan Hits for sensor $SENSOR_ID..."
                echo "    Input:  $OUTPUT_DIR/root/$ANALYSIS_ROOT"
                echo "    Output: $OUTPUT_DIR/hdf5/$OUTPUT_FILE"
		
		echo "OUTPUT_DIR=" $OUTPUT_DIR
		echo "ANALYSIS_ROOT=" $ANALYSIS_ROOT
		echo "OUTPUT_FILE=" $OUTPUT_FILE
		
                "${SCRIPT_DIR}/export_to_hdf5" --mode corry --input "$ANALYSIS_ROOT" --tree Analysis \
                    --output "$OUTPUT_FILE" --output-dir "$OUTPUT_DIR" \
                    --sensor-id "$SENSOR_ID" --sensor-mapping "$PIPELINE_CONFIG" \
                    --column-id 1

                if [ $? -ne 0 ]; then
                    echo "WARNING: Corryvreckan export failed for sensor $SENSOR_ID (continuing anyway)"
                fi
            done
        else
            # Export all channels to single file in Corryvreckan format
            echo "  Exporting Corryvreckan Hits (all channels)..."
            echo "    Input:  $OUTPUT_DIR/root/$ANALYSIS_ROOT"
            echo "    Output: $OUTPUT_DIR/hdf5/$ANALYSIS_HDF5"

            "${SCRIPT_DIR}/export_to_hdf5" --mode corry --input "$ANALYSIS_ROOT" --tree Analysis \
                --output "$ANALYSIS_HDF5" --output-dir "$OUTPUT_DIR" \
                --sensor-mapping "$PIPELINE_CONFIG" --column-id 1

            if [ $? -ne 0 ]; then
                echo "WARNING: Corryvreckan export failed (continuing anyway)"
            fi
        fi
        echo ""
    fi
fi

echo "=========================================="
echo "Pipeline completed successfully!"
echo "=========================================="
echo ""
echo "Output files:"
if [ -f "$OUTPUT_DIR/root/$RAW_ROOT" ]; then
    echo "  Raw ROOT:        $OUTPUT_DIR/root/$RAW_ROOT"
fi
if [ -f "$OUTPUT_DIR/root/$ANALYSIS_ROOT" ]; then
    echo "  Analysis ROOT:   $OUTPUT_DIR/root/$ANALYSIS_ROOT"
fi

# List HDF5 files (sensor-specific or single file)
HDF5_FILES=$(ls "$OUTPUT_DIR/hdf5/"*.h5 2>/dev/null)
if [ -n "$HDF5_FILES" ]; then
    echo "  HDF5 files:"
    for FILE in $HDF5_FILES; do
        echo "    - $(basename $FILE)"
    done
fi
echo ""
echo "Output directory structure:"
echo "  $OUTPUT_DIR/"
echo "    ├── root/            (ROOT files)"
echo "    ├── hdf5/            (HDF5 files)"
echo "    ├── waveform_plots/  (debug output, if enabled)"
echo "    └── temp/            (temporary files, if parallel processing)"
echo ""
