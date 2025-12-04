#!/bin/bash

# Parallel waveform analysis script
# Splits analysis into chunks and processes them in parallel

set -e

# Default values (can be overridden by config or command line)
DEFAULT_CONFIG="converter_config.json"
DEFAULT_CHUNK_SIZE=100
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

# Validate required arguments
if [ -z "$INPUT_ROOT" ] || [ -z "$OUTPUT_ROOT" ]; then
    echo "ERROR: --input and --output are required"
    print_usage
    exit 1
fi

if [ ! -f "$CONFIG" ]; then
    echo "ERROR: Config file not found: $CONFIG"
    exit 1
fi

if [ ! -f "./analyze_waveforms" ]; then
    echo "ERROR: analyze_waveforms executable not found"
    echo "Please run 'make' to build the executables"
    exit 1
fi

# Extract output_dir from config
OUTPUT_DIR=$(python3 -c "
import json, sys
try:
    with open('$CONFIG') as f:
        config = json.load(f)
        print(config.get('output_dir', 'output'))
except Exception:
    print('output')
" 2>/dev/null)

if [ -z "$OUTPUT_DIR" ]; then
    OUTPUT_DIR="output"
fi

# Build full input path
INPUT_PATH="$OUTPUT_DIR/root/$INPUT_ROOT"

if [ ! -f "$INPUT_PATH" ]; then
    echo "ERROR: Input file not found: $INPUT_PATH"
    exit 1
fi

echo "=========================================="
echo "Parallel Waveform Analysis"
echo "=========================================="
echo "Config:      $CONFIG"
echo "Input:       $INPUT_PATH"
echo "Output:      $OUTPUT_DIR/root/$OUTPUT_ROOT"
echo "Chunk size:  $CHUNK_SIZE events"
echo "Max cores:   $MAX_CORES"
echo "Temp dir:    $TEMP_DIR"
echo ""

# Get total number of events from ROOT file using temporary macro
echo "Counting events in input file..."
# cat > /tmp/count_events_$$.C <<MACRO_EOF
# void count_events_$$() {
#     TFile *f = TFile::Open("$INPUT_PATH");
#     if (f && !f->IsZombie()) {
#         TTree *t = (TTree*)f->Get("Waveforms");
#         if (t) {
#             std::cout << "NUM_ENTRIES=" << t->GetEntries() << std::endl;
#         }
#         f->Close();
#     }
# }
# MACRO_EOF

# NUM_EVENTS=$(root -l -b -q "/tmp/count_events_$$.C" 2>&1 | grep "NUM_ENTRIES=" | cut -d= -f2)
# rm -f /tmp/count_events_$$.C
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

    echo "  Chunk $CHUNK_ID: events [$START_EVENT, $END_EVENT)"

    ./analyze_waveforms \
        --config "$CONFIG" \
        --input "$INPUT_ROOT" \
        --output "$(basename $CHUNK_OUTPUT)" \
        --event-range "$START_EVENT:$END_EVENT" \
        > "$TEMP_DIR/chunk_${CHUNK_ID}.log" 2>&1

    # Move output to temp directory
    if [ -f "$OUTPUT_DIR/root/$(basename $CHUNK_OUTPUT)" ]; then
        mv "$OUTPUT_DIR/root/$(basename $CHUNK_OUTPUT)" "$CHUNK_OUTPUT"
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

# Use hadd to merge ROOT files
OUTPUT_PATH="$OUTPUT_DIR/root/$OUTPUT_ROOT"
hadd -f "$OUTPUT_PATH" "$TEMP_DIR"/chunk_*.root > "$TEMP_DIR/merge.log" 2>&1

if [ $? -ne 0 ]; then
    echo "ERROR: Failed to merge chunk files (see $TEMP_DIR/merge.log)"
    exit 1
fi

echo "Merged $NUM_CHUNKS chunks into $OUTPUT_PATH"
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
