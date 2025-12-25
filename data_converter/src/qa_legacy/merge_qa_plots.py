import ROOT
import os
import sys

def merge_and_compare(file0_path, file1_path, output_dir, num_events=5):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    f0 = ROOT.TFile.Open(file0_path, "READ")
    f1 = ROOT.TFile.Open(file1_path, "READ")

    if not f0 or f0.IsZombie():
        print(f"Error opening {file0_path}")
        return
    if not f1 or f1.IsZombie():
        print(f"Error opening {file1_path}")
        return

    # Get events directory
    dir0 = f0.Get("events")
    dir1 = f1.Get("events")

    if not dir0 or not dir1:
        print("Error: 'events' directory not found in one of the files")
        return

    # Assuming keys are "event_XXXXXX_quality_check"
    keys0 = [k.GetName() for k in dir0.GetListOfKeys() if k.GetName().startswith("event_")]
    keys0.sort()

    print(f"Found {len(keys0)} events in file 0. Processing first {num_events}...")

    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetPalette(ROOT.kRainBow)

    for key_name in keys0[:num_events]:
        print(f"Processing {key_name}...")
        
        c0 = dir0.Get(key_name)
        c1 = dir1.Get(key_name)

        if not c0 or not c1:
            print(f"  Missing canvas for {key_name} in one file. Skipping.")
            continue

        # Extract histograms
        hists = {} # Map sensor_id -> combined histogram

        # Helper to extract and add
        def extract_from_canvas(canvas):
            primitives = canvas.GetListOfPrimitives()
            for obj in primitives:
                # The canvas is divided into pads. We need to go into pads.
                if obj.InheritsFrom("TPad"):
                    pad_prims = obj.GetListOfPrimitives()
                    for prim in pad_prims:
                        if prim.InheritsFrom("TH2"):
                            name = prim.GetName()
                            # name is "sensorXX_amplitude_map"
                            if "sensor" in name:
                                sensor_id = int(name.split("sensor")[1].split("_")[0])
                                if sensor_id not in hists:
                                    hists[sensor_id] = prim.Clone(f"merged_sensor{sensor_id}_{key_name}")
                                    hists[sensor_id].SetDirectory(0)
                                else:
                                    hists[sensor_id].Add(prim)
        
        extract_from_canvas(c0)
        extract_from_canvas(c1)

        # Now draw merged histograms
        num_sensors = len(hists)
        if num_sensors == 0:
            continue

        c_out = ROOT.TCanvas(f"merged_{key_name}", f"Merged {key_name}", 600 * num_sensors, 600)
        c_out.Divide(num_sensors, 1)

        sorted_ids = sorted(hists.keys())
        for i, sid in enumerate(sorted_ids):
            c_out.cd(i + 1)
            h = hists[sid]
            h.SetTitle(f"Merged Sensor {sid} - {key_name}")
            # Ensure proper range if needed, or let auto-scale
            h.SetMinimum(0)
            h.SetMaximum(300) # Set to 300 for better color contrast
            h.Draw("COLZ TEXT")

        out_name = os.path.join(output_dir, f"merged_{key_name}.png")
        c_out.SaveAs(out_name)
        print(f"  Saved {out_name}")

    f0.Close()
    f1.Close()

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Usage: python3 merge_qa_plots.py <daq00_qc.root> <daq01_qc.root> <output_dir>")
        sys.exit(1)
    
    merge_and_compare(sys.argv[1], sys.argv[2], sys.argv[3])
