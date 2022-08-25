import os
import pysam
import pickle
from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser(description="Parse BAM files to pickle files for tombo input.")
parser.add_argument("-b", "--bam_input", metavar="path", type=str, help="path to bam input file", required=True)
parser.add_argument("-f", "--reference_file", metavar="path", type=str, help="path to bam input file", required=True)
parser.add_argument("-d", "--data_output_dir", metavar="path", type=str, help="path to bam input file", required=True)
parser.add_argument("-c", "--continue_work", action="store_true", help="allows resuming unfinished parsing")
parser.add_argument("-ow", "--overwrite", action="store_true", help="overwrites/deletes specificed dicitionary before starting")



def parse_BAM_to_pickle(bam_input, reference_file, data_output_dir, continue_opt = False, overwrite_opt=False):

    alignment_dir = paths["project_dir"] + "/" + paths["analysis_prefix"] + "/data/" + sample + "/bamtx/"
    dict_output_dir = data_output_dir + "/tombo_read_align_dicts/"
    os.mkdirs(dict_output_dir, exist_ok=True)

    filenames_present = (filename.split(".pickle")[0] for filename in os.listdir(dict_outputdir))
    num_files_in_dict = len(filenames_present)
    if num_files_in_dict>0:
        print("Already files present at output location")
        if not continue_opt and not overwrite_opt:
            print("Neither [continue_opt] nor [overwrite_opt] options are set, exiting without generating dicts")
            return 1
        elif overwrite_opt:
            print("Overwrite option set, deleting old and overwriting all data")
            os.system(f"rm -r {dict_output_dir}")
            os.mkdirs(dict_output_dir, exist_ok=True)
        elif continue_opt:
            print(f"Continuing from previous parsing, already {num_files_in_dict} files present")
        
    reference = SeqIO.read(reference_file, "fasta")
    
    start = 0
    end = len(reference.seq)
    chrm = reference.id
    
    ref_sequence_list = {}

    current_alignment = pysam.AlignmentFile(bam_input, "rb")
    
    for read in current_alignment.fetch(chrm, start, end):
        if continue_opt:
            if read.query_name in filenames_present:
                continue
        
        orientation = 1
        cigartuples = []
        if read.is_reverse:
            orientation = 0
        for operation, op_len in read.cigartuples:
            cigartuples.append((op_len, operation))
        attributes = {
            "ctg"    : read.reference_name,
            "r_st"   : read.reference_start, 
            "r_en"   : read.reference_end,
            "strand" : orientation,
            "mlen"   : read.query_alignment_length,
            "cigar"  : cigartuples,
            "q_st"   : read.query_alignment_start,
            "q_en"   : read.query_alignment_end,
            "ref_name": read.reference_name
        }
        
        with open(dict_output_dir + "/" + read.query_name + ".pickle", "wb") as outfile:
            pickle.dump(attributes, outfile)
        
    with open(f"{dict_output_dir}/{reference.id}.pickle", "wb") as outfile:
        pickle.dump(reference.seq.upper().__str__(), outfile)
    print(f"Wrote reference file {dict_output_dir}/{reference.id}.pickle")

    return dict_output_dir


if name == "__main__":
    args = parser.parse_args()
    parse_BAM_to_pickle(args.bam_input, args.reference_file, args.data_output_dir, args.continue_work, args.overwrite)