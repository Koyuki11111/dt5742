import os
import sys
import glob
from PIL import Image, ImageDraw, ImageFont

def stack_images(input_dir, output_dir, num_events=5):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    for i in range(num_events):
        event_id = i
        
        # Define file paths
        root_img_path = os.path.join(input_dir, f"merged_event_{event_id:06d}_quality_check.png")
        hdf5_img_path = os.path.join(input_dir, f"stitched_hdf5_event_{event_id:06d}.png")
        
        # Check if both exist
        if not os.path.exists(root_img_path):
            print(f"Missing ROOT image for event {event_id}: {root_img_path}")
            continue
        if not os.path.exists(hdf5_img_path):
            print(f"Missing HDF5 image for event {event_id}: {hdf5_img_path}")
            continue
            
        try:
            img_top = Image.open(root_img_path)
            img_bot = Image.open(hdf5_img_path)
            
            # Target width is the max of the two
            target_w = max(img_top.width, img_bot.width)
            
            # Resize bottom to match top width? Or just center them?
            # Let's center them on a canvas.
            
            total_h = img_top.height + img_bot.height + 60 # 60px for labels/spacing
            
            canvas = Image.new('RGB', (target_w, total_h), (255, 255, 255))
            
            # Draw labels
            draw = ImageDraw.Draw(canvas)
            # Try to load a default font, otherwise use default
            try:
                font = ImageFont.truetype("DejaVuSans-Bold.ttf", 20)
            except:
                font = ImageFont.load_default()

            # Paste Top (ROOT)
            x_top = (target_w - img_top.width) // 2
            y_top = 30
            canvas.paste(img_top, (x_top, y_top))
            draw.text((10, 5), "ROOT Merged (Original)", fill="black", font=font)
            
            # Paste Bottom (HDF5)
            x_bot = (target_w - img_bot.width) // 2
            y_bot = y_top + img_top.height + 30
            canvas.paste(img_bot, (x_bot, y_bot))
            draw.text((10, y_top + img_top.height + 5), "HDF5 Stitched (Reconstructed)", fill="black", font=font)
            
            out_name = os.path.join(output_dir, f"comparison_event_{event_id:06d}.png")
            canvas.save(out_name)
            print(f"Saved {out_name}")
            
        except Exception as e:
            print(f"Error processing event {event_id}: {e}")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python3 stack_comparison_images.py <input_dir> <output_dir>")
        sys.exit(1)
        
    stack_images(sys.argv[1], sys.argv[2])
