#!/usr/bin/env python3
import os
import re
import argparse

def update_files(directory, new_range):
    # Compile a regular expression pattern to find the restraintmask with numbers
    pattern = re.compile(r"(restraintmask=':)(\d+)-(\d+)(&!@H=')")
    
    # Convert new_range from string to tuple of integers
    new_start, new_end = map(int, new_range.split('-'))
    
    # Variables for feedback
    files_checked = 0
    files_updated = 0

    # Walk through all files in the directory
    for root, dirs, files in os.walk(directory):
        for file in files:
            if file.endswith(".in"):
                file_path = os.path.join(root, file)
                files_checked += 1
                
                try:
                    with open(file_path, 'r') as f:
                        contents = f.read()

                    # Replace old range with new range in the file contents
                    new_contents, count = pattern.subn(
                        lambda m: f"{m.group(1)}{new_start}-{new_end}{m.group(4)}", 
                        contents
                    )

                    # Only write back if a substitution was made
                    if count > 0:
                        try:
                            with open(file_path, 'w') as f:
                                f.write(new_contents)
                            print(f"Updated {file_path}")
                            files_updated += 1
                        except IOError as e:
                            print(f"Error writing to file {file_path}: {e}")
                except IOError as e:
                    print(f"Error reading file {file_path}: {e}")

    print(f"Total files checked: {files_checked}")
    print(f"Total files updated: {files_updated}")

def main():
    parser = argparse.ArgumentParser(description='Update restraintmask ranges in .in files.')
    parser.add_argument('new_range', help='New range in the format X-Y.')
    parser.add_argument('-d', '--directory', default=os.getcwd(), help='Directory to process (default: current directory).')
    args = parser.parse_args()

    update_files(args.directory, args.new_range)

if __name__ == "__main__":
    main()

