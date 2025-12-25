#!/bin/bash

# Parallel waveform analysis script
# Splits analysis into chunks and processes them in parallel

set -e

# Resolve script directory so we can find binaries regardless of caller cwd
SCRIPT_DIR="$(cd -- "$(dirname "$0")" && pwd)"

# Default values (can be overridden by config or command line)
DEFAULT_CONFIG="converter_config.json"
DEFAULT_CHUNK_SIZE=500
DEFAULT_MAX_CORES=8
DEFAULT_TEMP_DIR="./temp_analysis"

# Parse command line arguments
CONFIG="$DEFAULT_CONFIG"
INPUT_ROOT=""
OUTPUT_ROOT=""
CHUNK_SIZE=$DEFAULT_CHUNK_SIZE
MAX_CORES=$DEFAULT_MAX_CORES
TEMP_DIR=$DEFAULT_TEMP_DIR

print_usage() {
    cat << EOF
Parallel waveform analysis

Usage: $0 [options]

Options:
    --config FILE       Pipeline configuration file (default: converter_config.json)
    --input FILE        Input ROOT file name (relative to output_dir/root/)
    --output FILE       Output ROOT file name (relative to output_dir/root/)
    --chunk-size N      Events per chunk (default: 100)
    --max-cores N       Maximum parallel processes (default: 8)
    --temp-dir DIR      Temporary directory for chunks (default: ./temp_analysis)
    -h, --help          Show this help

Example:
    $0 --config converter_config.json --input waveforms.root --output waveforms_analyzed.root

EOF
}

while [[ $# -gt 0 ]]; do
    case $1 in
        --config)
            CONFIG="$2"
            shift 2
            ;;
        --input)
            INPUT_ROOT="$2"
            shift 2
            ;;
        --output)
            OUTPUT_ROOT="$2"
            shift 2
            ;;
        --chunk-size)
            CHUNK_SIZE="$2"
            shift 2
            ;;
        --max-cores)
            MAX_CORES="$2"
            shift 2
            ;;
        --temp-dir)
            TEMP_DIR="$2"
            shift 2
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

# Validate arguments
if [ ! -f "$CONFIG" ]; then
    echo "ERROR: Config file not found: $CONFIG"
    exit 1
fi

ANALYZE_BIN="${SCRIPT_DIR}/analyze_waveforms"

if [ ! -x "$ANALYZE_BIN" ]; then
    echo "ERROR: analyze_waveforms executable not found at $ANALYZE_BIN"
    echo "Please run 'make' to build the executables"
    exit 1
fi

daq_name=$(sed -nE 's/^[[:space:]]*"daq_name"[[:space:]]*:[[:space:]]*"([^"]+)".*/\1/p' "$CONFIG" | head -n1)
runnumber=$(sed -nE 's/^[[:space:]]*"runnumber"[[:space:]]*:[[:space:]]*([0-9]+).*/\1/p' "$CONFIG" | head -n1)
output_dir=$(sed -nE 's/^[[:space:]]*"output_dir"[[:space:]]*:[[:space:]]*"([^"]+)".*/\1/p' "$CONFIG" | head -n1)
run_str=$(printf "%06d" "$runnumber")
OUTPUT_DIR="${output_dir}/${run_str}/${daq_name}"

# If inputs not specified, try to read from config
if [ -z "$INPUT_ROOT" ]; then
    INPUT_ROOT=$(sed -nE 's/^[[:space:]]*"waveforms_root"[[:space:]]*:[[:space:]]*"([^"]+)".*/\1/p' "$CONFIG" | head -n1)
fi
if [ -z "$OUTPUT_ROOT" ]; then
    OUTPUT_ROOT=$(sed -nE 's/^[[:space:]]*"analysis_root"[[:space:]]*:[[:space:]]*"([^"]+)".*/\1/p' "$CONFIG" | head -n1)
fi

if [ -z "$INPUT_ROOT" ] || [ -z "$OUTPUT_ROOT" ]; then
    echo "ERROR: --input and --output are required (or must be in config)"
    print_usage
    exit 1
fi

# Build full input path
INPUT_PATH="$OUTPUT_DIR/output/root/$INPUT_ROOT"

if [ ! -f "$INPUT_PATH" ]; then
    echo "ERROR: Input file not found: $INPUT_PATH"
    exit 1
fi

echo "=========================================="
echo "Parallel Waveform Analysis"
echo "=========================================="
echo "Config:      $CONFIG"
echo "Input:       $INPUT_PATH"
echo "Output:      $OUTPUT_DIR/output/root/$OUTPUT_ROOT"
echo "Chunk size:  $CHUNK_SIZE events"
echo "Max cores:   $MAX_CORES"
echo "Temp dir:    $TEMP_DIR"
echo ""

# Get total number of events from ROOT file using temporary macro
echo "Counting events in input file..."

NUM_EVENTS=$(
INPUT_PATH="$INPUT_PATH" root -l -b -q -e '
    TFile *f = TFile::Open(gSystem->Getenv("INPUT_PATH"));
    if (f && !f->IsZombie()) {
        TTree *t = nullptr;
        f->GetObject("Waveforms", t);
        if (t) {
            std::cout << "NUM_ENTRIES=" << t->GetEntries() << std::endl;
        }
        f->Close();
    }
    gSystem->Exit(0);
' 2>&1 | grep "NUM_ENTRIES=" | cut -d= -f2
)

if [ -z "$NUM_EVENTS" ] || [ "$NUM_EVENTS" -le 0 ]; then
    echo "ERROR: Could not determine number of events in input file"
    echo "       Tried to read: $INPUT_PATH"
    exit 1
fi

echo "Total events: $NUM_EVENTS"
echo ""

# Calculate number of chunks
NUM_CHUNKS=$(( (NUM_EVENTS + CHUNK_SIZE - 1) / CHUNK_SIZE ))
echo "Processing in $NUM_CHUNKS chunks..."
echo ""

# Create temporary directory
rm -rf "$TEMP_DIR"
mkdir -p "$TEMP_DIR"

# Function to process a single chunk
process_chunk() {
    local CHUNK_ID=$1
    local START_EVENT=$2
    local END_EVENT=$3

    local CHUNK_OUTPUT="$TEMP_DIR/chunk_${CHUNK_ID}.root"
    local CHUNK_PLOTS="waveform_plots_chunk_${CHUNK_ID}"

    echo "  Chunk $CHUNK_ID: events [$START_EVENT, $END_EVENT)"

    "$ANALYZE_BIN" \
        --config "$CONFIG" \
        --input "$INPUT_ROOT" \
        --output "$(basename $CHUNK_OUTPUT)" \
        --event-range "$START_EVENT:$END_EVENT" \
        --waveform-plots-file "$CHUNK_PLOTS" \
        > "$TEMP_DIR/chunk_${CHUNK_ID}.log" 2>&1    
    
    # Move output(analysis result) to temp directory
    if [ -f "$OUTPUT_DIR/output/root/$(basename $CHUNK_OUTPUT)" ]; then	
	    mv "$OUTPUT_DIR/output/root/$(basename $CHUNK_OUTPUT)" "$CHUNK_OUTPUT"
    else
	    echo "ERROR: Missing $OUTPUT_DIR/$(basename $CHUNK_OUTPUT)"
        return 1
    fi
    
    # Move waveform plots to temp directory if it exists
    # If plots override output is enabled, it should be in waveform_plots dir
    # Actually analyze_waveforms puts it in output_dir/output/waveform_plots/[name].root
    
    if [ -f "$OUTPUT_DIR/output/waveform_plots/${CHUNK_PLOTS}.root" ]; then
        mv "$OUTPUT_DIR/output/waveform_plots/${CHUNK_PLOTS}.root" "$TEMP_DIR/${CHUNK_PLOTS}.root"
    fi
    # Also handle possibility of split files (.root, _001.root ...) if large
    # But for chunks it's unlikely to be > 4GB unless chunk is huge
    
    # Move quality check files to temp directory if they exist
    local CHUNK_QC=$(echo "$CHUNK_PLOTS" | sed 's/waveform_plots/quality_check/')
    if [ -f "$OUTPUT_DIR/output/quality_check/${CHUNK_QC}.root" ]; then
        mv "$OUTPUT_DIR/output/quality_check/${CHUNK_QC}.root" "$TEMP_DIR/${CHUNK_QC}.root"
    fi

    if [ $? -eq 0 ]; then
        echo "  Chunk $CHUNK_ID: DONE"
    else
        echo "  Chunk $CHUNK_ID: FAILED (see $TEMP_DIR/chunk_${CHUNK_ID}.log)"
        return 1
    fi
}

export -f process_chunk
export CONFIG INPUT_ROOT OUTPUT_DIR TEMP_DIR

# Process chunks in parallel
PIDS=()
CHUNK_ID=0

for ((START=0; START<NUM_EVENTS; START+=CHUNK_SIZE)); do
    END=$((START + CHUNK_SIZE))
    if [ $END -gt $NUM_EVENTS ]; then
        END=$NUM_EVENTS
    fi

    # Wait if we've reached max parallel processes
    while [ ${#PIDS[@]} -ge $MAX_CORES ]; do
        # Check if any process has finished
        for i in "${!PIDS[@]}"; do
            if ! kill -0 "${PIDS[$i]}" 2>/dev/null; then
                wait "${PIDS[$i]}"
                unset 'PIDS[$i]'
            fi
        done
        PIDS=("${PIDS[@]}")  # Re-index array
        sleep 0.1
    done

    # Start new chunk processing in background
    process_chunk $CHUNK_ID $START $END &
    PIDS+=($!)

    CHUNK_ID=$((CHUNK_ID + 1))
done

# Wait for all remaining processes
echo ""
echo "Waiting for all chunks to complete..."
for PID in "${PIDS[@]}"; do
    wait $PID
done

echo "All chunks processed."
echo ""

# Merge chunk results
echo "Merging results..."

# Check that all chunk files exist
MISSING_CHUNKS=0
for ((i=0; i<NUM_CHUNKS; i++)); do
    if [ ! -f "$TEMP_DIR/chunk_${i}.root" ]; then
        echo "ERROR: Missing chunk file: chunk_${i}.root"
        MISSING_CHUNKS=$((MISSING_CHUNKS + 1))
    fi
done

if [ $MISSING_CHUNKS -gt 0 ]; then
    echo "ERROR: $MISSING_CHUNKS chunk(s) failed. Cannot merge."
    exit 1
fi

# Use hadd to merge analysis ROOT files
OUTPUT_PATH="$OUTPUT_DIR/output/root/$OUTPUT_ROOT"
hadd -f "$OUTPUT_PATH" "$TEMP_DIR"/chunk_*.root > "$TEMP_DIR/merge.log" 2>&1

if [ $? -ne 0 ]; then
    echo "ERROR: Failed to merge chunk files (see $TEMP_DIR/merge.log)"
    exit 1
fi

echo "Merged $NUM_CHUNKS analysis chunks into $OUTPUT_PATH"
echo ""

# Handle waveform plots (SKIP MERGE by default, copy instead)
# Check if --merge-plots is passed? No, user asked generally.
# I'll implement a flag variable $MERGE_PLOTS and add it to args.
if [ "$MERGE_PLOTS" = "true" ]; then
    PLOTS_FILES=("$TEMP_DIR"/waveform_plots_chunk_*.root)
    if [ -f "${PLOTS_FILES[0]}" ]; then
        echo "Merging waveform plots files (this may take time)..."
        PLOTS_OUTPUT="$OUTPUT_DIR/output/waveform_plots/waveform_plots.root"
        mkdir -p "$OUTPUT_DIR/output/waveform_plots"
        hadd -f "$PLOTS_OUTPUT" "$TEMP_DIR"/waveform_plots_chunk_*.root > "$TEMP_DIR/merge_plots.log" 2>&1
        if [ $? -eq 0 ]; then
            echo "Merged waveform plots into $PLOTS_OUTPUT"
        else
            echo "WARNING: Failed to merge waveform plots (see $TEMP_DIR/merge_plots.log)"
        fi
    fi
else
    # Copy instead of merge
    PLOTS_FILES=("$TEMP_DIR"/waveform_plots_chunk_*.root)
    if [ -f "${PLOTS_FILES[0]}" ]; then
        echo "Copying waveform plots chunks (skipping merge)..."
        mkdir -p "$OUTPUT_DIR/output/waveform_plots"
        # We rename them to be meaningful if possible, or just chunk_N
        # Maybe use the start/end event if we tracked it, but chunk_ID is simpler
        cp "$TEMP_DIR"/waveform_plots_chunk_*.root "$OUTPUT_DIR/output/waveform_plots/"
        echo "Copied plot chunks to $OUTPUT_DIR/output/waveform_plots/"
    fi
fi

# Merge quality check files (usually small so ok to merge)
QC_FILES=("$TEMP_DIR"/quality_check_chunk_*.root)
if [ -f "${QC_FILES[0]}" ]; then
    echo "Merging quality check files..."
    QC_OUTPUT="$OUTPUT_DIR/output/quality_check/quality_check.root"
    mkdir -p "$OUTPUT_DIR/output/quality_check"
    hadd -f "$QC_OUTPUT" "$TEMP_DIR"/quality_check_chunk_*.root > "$TEMP_DIR/merge_qc.log" 2>&1
    if [ $? -eq 0 ]; then
        echo "Merged quality check files into $QC_OUTPUT"
    else
        echo "WARNING: Failed to merge quality check files (see $TEMP_DIR/merge_qc.log)"
    fi
fi
echo ""

# Clean up temporary files (optional)
echo "Cleaning up temporary files..."
# Uncomment to remove temp directory:
# rm -rf "$TEMP_DIR"
echo "Temporary files kept in: $TEMP_DIR"
echo ""

echo "=========================================="
echo "Parallel analysis complete!"
echo "=========================================="
echo "Output: $OUTPUT_PATH"
echo ""
