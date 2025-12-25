import os
import sys
import glob
from PIL import Image

def stitch_images(input_dir, output_dir, num_events=5):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Sensors 0, 1, 2, 3
    sensors = [0, 1, 2, 3]

    for i in range(num_events):
        event_id = i
        images = []
        
        # Find images for this event for each sensor
        for sid in sensors:
            # Pattern: event_000000_sensor0_sensor0_merged_analysis.png
            pattern = os.path.join(input_dir, f"event_{event_id:06d}_sensor{sid}_*.png")
            files = glob.glob(pattern)
            if files:
                images.append(Image.open(files[0]))
            else:
                print(f"Warning: Missing plot for Event {event_id}, Sensor {sid}")
                images.append(None)
        
        # Stitch them horizontally
        if not any(images):
            continue

        # Filter None
        valid_images = [img for img in images if img is not None]
        if not valid_images:
            continue

        # Assume all same size
        w, h = valid_images[0].size
        total_w = w * len(valid_images)
        
        stitched = Image.new('RGB', (total_w, h), (255, 255, 255))
        
        x_offset = 0
        for img in valid_images:
            stitched.paste(img, (x_offset, 0))
            x_offset += w
            
        out_name = os.path.join(output_dir, f"stitched_hdf5_event_{event_id:06d}.png")
        stitched.save(out_name)
        print(f"Saved {out_name}")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python3 stitch_hdf5_plots.py <input_dir> <output_dir>")
        sys.exit(1)
        
    stitch_images(sys.argv[1], sys.argv[2])
