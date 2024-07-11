import os
import sys
import pdb

def rename_files_in_directory(directory_path, ext):
    for filename in os.listdir(directory_path):
        if filename.endswith(ext):
            old_filepath = os.path.join(directory_path, filename)
            
            # Replacing dashes with underscores in the file basename
            new_filename = f"{filename.split('hmmalign')[0]}hmmalign_dashmergeseed{filename.split('hmmalign')[1]}"
            
            new_filepath = os.path.join(directory_path, new_filename)
            
            os.rename(old_filepath, new_filepath)

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python rename_files.py <directory_path>")
        sys.exit(1)

    directory_path = sys.argv[1]
    ext = ".a3m"

    # Calling the function to rename files
    rename_files_in_directory(directory_path, ext)
