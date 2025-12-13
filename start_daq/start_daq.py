#!/usr/bin/env python3
"""
DAQ Startup Script
Starts DAQ acquisition with tmux session management
"""

import os
import sys
import subprocess
import json
import argparse
import shutil
import time


def check_tmux_session(session_name):
    """Check if tmux session exists"""
    result = subprocess.run(
        ['tmux', 'has-session', '-t', session_name],
        capture_output=True,
        text=True
    )
    return result.returncode == 0


def kill_tmux_session(session_name):
    """Kill tmux session"""
    subprocess.run(['tmux', 'kill-session', '-t', session_name])


def create_folder_structure(run_identifier):
    """Create folder structure for DAQ data

    Args:
        run_identifier: Either a number (will be formatted as 6-digit) or a string

    Returns:
        Path to created folder
    """
    # Determine folder name
    try:
        run_number = int(run_identifier)
        folder_name = f"{run_number:06d}"
    except ValueError:
        folder_name = run_identifier

    # Create base folder
    base_path = f"/data/{folder_name}"

    # Check if folder already exists
    if os.path.exists(base_path):
        print(f"Folder '{base_path}' already exists.")
        response = input("Delete existing folder and continue? (y/n): ").strip().lower()
        if response == 'y':
            print(f"Deleting '{base_path}'...")
            shutil.rmtree(base_path)
        else:
            print("Aborting.")
            sys.exit(1)

    daq00_path = os.path.join(base_path, "daq00")
    daq01_path = os.path.join(base_path, "daq01")

    # Create directories
    os.makedirs(daq00_path, exist_ok=True)
    os.makedirs(daq01_path, exist_ok=True)

    print(f"Created folder structure: {base_path}")
    print(f"  - {daq00_path}")
    print(f"  - {daq01_path}")

    return base_path


def update_monitor_config(config_path, daq_path, daq_number):
    """Update monitor config file with correct input path

    Args:
        config_path: Path to config file
        daq_path: Path to daq folder
        daq_number: DAQ number (01 or 02)
    """
    with open(config_path, 'r') as f:
        config = json.load(f)

    # Update input_file path in monitor section
    if 'monitor' not in config:
        config['monitor'] = {}

    config['monitor']['input_file'] = f"{daq_path}/wave_0.dat"
    config['monitor']['log_file'] = f"{daq_path}/monitor.log"

    with open(config_path, 'w') as f:
        json.dump(config, f, indent=2)

    print(f"Updated config file: {config_path}")


def setup_tmux_session(base_path, args):
    """Setup tmux session with panes for DAQ operations

    Layout:
    ┌─────────────┬─────────────┐
    │             │             │
    │ DAQ00       │ DAQ01       │
    │             │             │
    ├─────┬───────┤             │
    │Mon1 │       │             │
    ├─────┤       │ glances     │
    │Mon2 │       │             │
    └─────┴───────┴─────────────┘

    Args:
        base_path: Base path to data folder
    """
    session_name = "caen_daq"
    daq00_path = os.path.join(base_path, "daq00")
    daq01_path = os.path.join(base_path, "daq01")

    # Copy config files using Python (to avoid race conditions with tmux)
    # Copy WaveDump configs
    shutil.copy(
        args.wavedump_usb0,
        os.path.join(daq00_path, 'WaveDumpConfig_USB0.txt')
    )

    shutil.copy(
        args.wavedump_usb1,
        os.path.join(daq01_path, 'WaveDumpConfig_USB1.txt')
    )    
    #shutil.copy('/opt/WaveDumpConfig_USB0.txt', os.path.join(daq00_path, 'WaveDumpConfig_USB0.txt'))
    #shutil.copy('/opt/WaveDumpConfig_USB1.txt', os.path.join(daq01_path, 'WaveDumpConfig_USB1.txt'))
    print(f"Copied WaveDump config files")

    # Copy and update monitor configs
    config_path_01 = os.path.join(daq00_path, 'monitor_config.json')
    config_path_02 = os.path.join(daq01_path, 'monitor_config.json')

    shutil.copy('/opt/dt5742/daq_monitor/monitor_config.json', config_path_01)
    update_monitor_config(config_path_01, daq00_path, '00')

    shutil.copy('/opt/dt5742/daq_monitor/monitor_config.json', config_path_02)
    update_monitor_config(config_path_02, daq01_path, '01')

    # Create new tmux session and get the initial pane ID
    result = subprocess.run(
        ['tmux', 'new-session', '-d', '-s', session_name, '-P', '-F', '#{pane_id}'],
        capture_output=True, text=True
    )
    pane_left = result.stdout.strip()  # Initial pane (will be left side)

    # Split 50:50 horizontally (left | right) - creates new pane on the right
    result = subprocess.run(
        ['tmux', 'split-window', '-h', '-t', pane_left, '-P', '-F', '#{pane_id}'],
        capture_output=True, text=True
    )
    pane_right = result.stdout.strip()  # Right side pane
    # pane_left stays as left side

    # Now split LEFT side vertically (top | bottom)
    result = subprocess.run(
        ['tmux', 'split-window', '-v', '-t', pane_left, '-P', '-F', '#{pane_id}'],
        capture_output=True, text=True
    )
    pane_left_bottom = result.stdout.strip()  # Left bottom (monitor area)
    pane_left_top = pane_left  # Left top (DAQ00)

    # Split left bottom into two monitors (Monitor1 | Monitor2)
    result = subprocess.run(
        ['tmux', 'split-window', '-v', '-t', pane_left_bottom, '-P', '-F', '#{pane_id}'],
        capture_output=True, text=True
    )
    pane_monitor2 = result.stdout.strip()  # Monitor2 (bottom)
    pane_monitor1 = pane_left_bottom  # Monitor1 (top)

    # Now split RIGHT side vertically (DAQ01 top | glances bottom)
    result = subprocess.run(
        ['tmux', 'split-window', '-v', '-t', pane_right, '-P', '-F', '#{pane_id}'],
        capture_output=True, text=True
    )
    pane_glances = result.stdout.strip()  # glances (right bottom)
    pane_right_top = pane_right  # DAQ01 (right top)

    # Resize panes to match heights - DAQ panes should be larger
    # Resize left top (DAQ00) to 60% of left column
    subprocess.run(['tmux', 'resize-pane', '-t', pane_left_top, '-y', '60%'])
    # Resize right top (DAQ01) to 60% of right column
    subprocess.run(['tmux', 'resize-pane', '-t', pane_right_top, '-y', '60%'])

    # Configure DAQ00 pane (top-left)
    subprocess.run(['tmux', 'send-keys', '-t', pane_left_top,
                   f'cd {daq00_path} && wavedump WaveDumpConfig_USB0.txt', 'C-m'])

    # Wait for DAQ00 to initialize
    print("Starting DAQ00...")
    time.sleep(1)

    # Configure DAQ01 pane (top-right)
    subprocess.run(['tmux', 'send-keys', '-t', pane_right_top,
                   f'cd {daq01_path} && wavedump WaveDumpConfig_USB1.txt', 'C-m'])

    # Wait for DAQ01 to initialize
    print("Starting DAQ01...")
    time.sleep(1)

    # Configure Monitor1 pane (left-bottom-top)
    subprocess.run(['tmux', 'send-keys', '-t', pane_monitor1,
                   f'cd {daq00_path} && /opt/dt5742/daq_monitor/monitor_realtime --config monitor_config.json', 'C-m'])

    # Configure Monitor2 pane (left-bottom-bottom)
    subprocess.run(['tmux', 'send-keys', '-t', pane_monitor2,
                   f'cd {daq01_path} && /opt/dt5742/daq_monitor/monitor_realtime --config monitor_config.json', 'C-m'])

    # Configure glances pane (right-bottom)
    subprocess.run(['tmux', 'send-keys', '-t', pane_glances, 'glances', 'C-m'])

    # Return session name
    return session_name


def main():
    parser = argparse.ArgumentParser(
        description='Start DAQ acquisition with tmux session management'
    )
    parser.add_argument(
        'run_identifier',
        help='Run number (will be formatted as 6-digit) or folder name (string)'
    )
    parser.add_argument("--wavedump-usb0", default="/opt/WaveDumpConfig_USB0.txt")
    parser.add_argument("--wavedump-usb1", default="/opt/WaveDumpConfig_USB1.txt")
    #parser.add_argument("--monitor-config", default="/opt/dt5742/daq_monitor/monitor_config.json")
    
    args = parser.parse_args()

    session_name = "caen_daq"

    # Check if tmux session already exists
    if check_tmux_session(session_name):
        print(f"Tmux session '{session_name}' already exists.")
        print("Options:")
        print("  1) Kill and restart")
        print("  2) Attach to existing session")
        print("  3) Cancel")
        response = input("Select option (1/2/3): ").strip()

        if response == '1':
            print(f"Killing existing session '{session_name}'...")
            kill_tmux_session(session_name)
        elif response == '2':
            print(f"Attaching to existing session '{session_name}'...")
            subprocess.run(['tmux', 'attach', '-t', session_name])
            sys.exit(0)
        else:
            print("Aborting.")
            sys.exit(1)

    # Create folder structure
    base_path = create_folder_structure(args.run_identifier)
    
    # Setup tmux session
    session_name = setup_tmux_session(base_path,args)

    # Wait for processes to start
    time.sleep(1)

    # Display session info
    daq00_path = os.path.join(base_path, "daq00")
    daq01_path = os.path.join(base_path, "daq01")

    print("\n" + "="*63)
    print(f"║         Tmux Session: {session_name} (5 panes)            ║")
    print(f"║         Working Directory: {base_path:<28}║")
    print("╠═════════════════════════════╦═════════════════════════════╣")
    print("║  DAQ00 WaveDump             ║  DAQ01 WaveDump             ║")
    print("║  Digitizer USB0             ║  Digitizer USB1             ║")
    print(f"║  {daq00_path:<27}║  {daq01_path:<27}║")
    print("╠═════════════════════════════╬═════════════════════════════╣")
    print("║  DAQ00 Monitor              ║                             ║")
    print("║  Real-time Analysis         ║  System Monitor             ║")
    print("╟─────────────────────────────╢  glances                    ║")
    print("║  DAQ01 Monitor              ║                             ║")
    print("║  Real-time Analysis         ║                             ║")
    print("="*63)
    print("\nPress ENTER to attach to tmux session...")
    input()

    # Attach to tmux session
    subprocess.run(['tmux', 'attach', '-t', session_name])


if __name__ == '__main__':
    main()
