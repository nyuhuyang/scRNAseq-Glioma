import argparse
parser = argparse.ArgumentParser(description='select a integer, start from 1')
parser.add_argument("ID", help="Select file to process according to SGE_TASK_ID integer",
                    type=int)
args = parser.parse_args()
print(args.ID)
import os
import loompy

path="/athena/elementolab/scratch/yah2014/Projects/scRNAseq-Glioma/data"
os.chdir(path) # change current path
print(os.getcwd())
# List all filer folder's names.
file_folders=os.listdir(os.getcwd())  # list files
file=file_folders[args.ID-1]
print(file)

# Select the file according to SGE_TASK_ID
file_path=os.path.join(path, file, "velocyto",file+'.loom')
print("file_path= "+file_path)


