#!/bin/bash

# Full waveform processing pipeline
# Stage 1: Convert binary/ASCII to ROOT
# Stage 2: Analyze waveforms
# Stage 3: Export to HDF5

set -e  # Exit on error

# Default configuration file (unified)
PIPELINE_CONFIG="converter_config.json"

# Output directory is read from JSON configs (default: "output")
# Edit converter_config.json (common.output_dir) to change output_dir
OUTPUT_DIR="output"

# Default intermediate/output files (relative names). For a minimal smoke test,
# drop a few tiny mixed-length inputs under ./temp/smoke/ (see README) and run:
#   ./run_full_pipeline.sh --config converter_config.json
RAW_ROOT="waveforms.root"
ANALYSIS_ROOT="waveforms_analyzed.root"
ANALYSIS_HDF5="waveforms_analyzed.hdf5"

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
    --analysis-hdf5 FILE     Analysis HDF5 filename (default: waveforms_analyzed.hdf5)
    --stage1-only            Run only stage 1 (convert to ROOT)
    --stage2-only            Run only stage 2 (analysis)
    --stage3-only            Run only stage 3 (export to HDF5)
    --skip-stage1            Skip stage 1
    --skip-stage2            Skip stage 2
    --skip-stage3            Skip stage 3
    --parallel               Force parallel processing (overrides config auto-detection)
    --verbose                Verbose output
    -h, --help               Show this help message

Parallel Processing:
    By default, the pipeline auto-detects parallel mode from config files:
      - analysis_config.json: parallel_max_cores > 1 enables parallel Stage 2
      - wave_converter_config.json: max_threads > 1 enables parallel Stage 1
    Use --parallel to force parallel mode regardless of config settings.

Stages:
    Stage 1: convert_to_root     - Convert binary/ASCII waveforms to ROOT format
    Stage 2: analyze_waveforms   - Extract timing and amplitude features
    Stage 3: export_to_hdf5      - Export analyzed ROOT data to HDF5 format

Output Organization:
    All outputs are organized in subdirectories (default: output/):
      output/root/            - ROOT files (waveforms.root, waveforms_analyzed.root)
      output/hdf5/            - HDF5 files (waveforms_analyzed.hdf5 and sensor splits)
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
            USE_PARALLEL=true
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

    if [ ! -f "./convert_to_root" ]; then
        echo "ERROR: convert_to_root executable not found"
        echo "Please run 'make' to build the executables"
        exit 1
    fi

    # Note: convert_to_root reads output_dir from config and automatically creates output_dir/root/
    if [ "$USE_PARALLEL" = true ]; then
        ./convert_to_root --config "$PIPELINE_CONFIG" --root "$RAW_ROOT" --parallel
    else
        ./convert_to_root --config "$PIPELINE_CONFIG" --root "$RAW_ROOT"
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

        if [ ! -x "./parallel_analyze.sh" ]; then
            echo "ERROR: parallel_analyze.sh not found or not executable"
            echo "Please run 'make' to build the executables"
            exit 1
        fi

        ./parallel_analyze.sh --config "$PIPELINE_CONFIG" --input "$RAW_ROOT" --output "$ANALYSIS_ROOT"

        if [ $? -ne 0 ]; then
            echo "ERROR: Stage 2 (parallel) failed"
            exit 1
        fi
    else
        echo "  Mode:   SEQUENTIAL"
        echo ""

        if [ ! -f "./analyze_waveforms" ]; then
            echo "ERROR: analyze_waveforms executable not found"
            echo "Please run 'make' to build the executables"
            exit 1
        fi

        ./analyze_waveforms --config "$PIPELINE_CONFIG" --input "$RAW_ROOT" --output "$ANALYSIS_ROOT"

        if [ $? -ne 0 ]; then
            echo "ERROR: Stage 2 failed"
            exit 1
        fi
    fi
    echo ""
fi

# Stage 3: Export to HDF5
if [ "$RUN_STAGE3" = true ]; then
    echo "Stage 3: Exporting analyzed data to HDF5 format..."

    if [ ! -f "./export_to_hdf5" ]; then
        echo "ERROR: export_to_hdf5 executable not found"
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
        stage2 = config.get('stage2', {})
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

    # Export analysis features if they exist
    if [ -f "$OUTPUT_DIR/root/$ANALYSIS_ROOT" ] && [ "$RUN_STAGE2" = true ] || [ "$RUN_STAGE2" = false ]; then
        if [ -n "$SENSOR_IDS" ]; then
            # Export per sensor
            for SENSOR_ID in $SENSOR_IDS; do
                OUTPUT_FILE="waveforms_analyzed_sensor$(printf '%02d' $SENSOR_ID).hdf5"
                echo "  Exporting analysis features for sensor $SENSOR_ID..."
                echo "    Input:  $OUTPUT_DIR/root/$ANALYSIS_ROOT"
                echo "    Output: $OUTPUT_DIR/hdf5/$OUTPUT_FILE"

                ./export_to_hdf5 --mode analysis --input "$ANALYSIS_ROOT" --tree Analysis \
                    --output "$OUTPUT_FILE" --output-dir "$OUTPUT_DIR" \
                    --sensor-id "$SENSOR_ID" --sensor-mapping "$PIPELINE_CONFIG"

                if [ $? -ne 0 ]; then
                    echo "WARNING: Analysis feature export failed for sensor $SENSOR_ID (continuing anyway)"
                fi
            done
        else
            # Export all channels to single file
            echo "  Exporting analysis features..."
            echo "    Input:  $OUTPUT_DIR/root/$ANALYSIS_ROOT"
            echo "    Output: $OUTPUT_DIR/hdf5/$ANALYSIS_HDF5"

            ./export_to_hdf5 --mode analysis --input "$ANALYSIS_ROOT" --tree Analysis \
                --output "$ANALYSIS_HDF5" --output-dir "$OUTPUT_DIR"

            if [ $? -ne 0 ]; then
                echo "WARNING: Analysis feature export failed (continuing anyway)"
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
HDF5_FILES=$(ls "$OUTPUT_DIR/hdf5/"*.hdf5 2>/dev/null)
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
