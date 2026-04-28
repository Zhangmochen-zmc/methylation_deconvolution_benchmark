import os
import glob
import subprocess
import sys
import argparse
import time


def main():
    # 1. Define parameters
    parser = argparse.ArgumentParser(
        description="MEnet predict",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # config
    parser.add_argument('-m', '--model', required=True, 
                        help="Path to the trained model file (.pkl)")

    # Optional parameters
    parser.add_argument('-i', '--input_dir', default='test_data', 
                        help="Enter the folder path (including the CSV file).")
    
    parser.add_argument('-o', '--output_dir', default='menet_result', 
                        help="Path to the output folder")
    
    parser.add_argument('--bedtools', 
                        help="Specify the absolute path to bedtools")

    args = parser.parse_args()

    # 2. check
    if not os.path.exists(args.model):
        print(f"[Error] Model file not found: {args.model}")
        sys.exit(1)

    # Check the input directory
    if not os.path.exists(args.input_dir):
        print(f"[Error] The input directory does not exist.: {args.input_dir}")
        sys.exit(1)

    # Create output directory
    if not os.path.exists(args.output_dir):
        print(f"[-] Create output directory: {args.output_dir}")
        os.makedirs(args.output_dir)

    # 3. Get a list of CSV files
    search_pattern = os.path.join(args.input_dir, "*.csv")
    csv_files = glob.glob(search_pattern)

    if not csv_files:
        print(f"[Warning]  No .csv file was found in '{args.input_dir}'.")
        sys.exit(0)

    print(f"==========================================")
    print(f"{len(csv_files)} files were detected.")
    print(f"model: {args.model}")
    print(f"output: {args.output_dir}")
    print(f"==========================================\n")

    # 4. run
    success_count = 0
    fail_count = 0

    for idx, csv_file in enumerate(csv_files):
        file_name = os.path.basename(csv_file)
        print(f"[{idx+1}/{len(csv_files)}] processing: {file_name} ...")
        
        start_time = time.time()

        # Build command
        # MEnet predict -i input.csv -m model.pkl -o output_dir --input_type table
        cmd = [
            "MEnet", "predict",
            "-i", csv_file,
            "-m", args.model,
            "-o", args.output_dir,
            "--input_type", "table"  
        ]

        # If the user specifies the bedtools path
        if args.bedtools:
            cmd.extend(["--bedtools", args.bedtools])

        try:
            # run

            env = os.environ.copy()
            env["CUDA_VISIBLE_DEVICES"] = ""
            
            subprocess.run(cmd, check=True, env=env)
            
            elapsed = time.time() - start_time
            print(f"  -> done (time {elapsed:.2f}s)")
            success_count += 1

        except subprocess.CalledProcessError as e:
            print(f"  >>> [Error] Processing failure: {file_name}")
            fail_count += 1
        except FileNotFoundError:
            print(f"  >>> [Error] The 'MEnet' command could not be found. Please ensure that the environment is installed and activated.")
            sys.exit(1)
        except Exception as e:
            print(f"  >>> [Error] Unknown error: {e}")
            fail_count += 1

    print("\n" + "="*30)
    print(f"success: {success_count}, failure: {fail_count}")
    print("="*30)

if __name__ == "__main__":
    main()
