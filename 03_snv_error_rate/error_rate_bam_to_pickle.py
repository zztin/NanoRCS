import pysam
import random
import os
import pandas as pd
import argparse

BATCH_SIZE = 100_000
SET_BATCH_COUNT = 10
# pysam.index(filename)

def get_error_df(filename, outpath, BATCH_SIZE, SET_BATCH_COUNT):
    filename_prefix = filename.split(".bam")[0].split("/")[-1]
    inpath = filename.split(".bam")[0].rsplit("/",1)[0]
    all_reads = []
    # Read all reads in bam file into python memory
    in_bam = pysam.AlignmentFile(filename, "rb")
    for read in in_bam:
        # Only record the main chromosomes
        if read.reference_name not in [str(x) for x in range(1,23)]:
            continue
        if read.mapping_quality < 60:
            continue
        else:
            all_reads.append(read)
    print(len(all_reads), "has mapQ >= 60.")

    # Choose batch size
    batch_count = len(all_reads) // BATCH_SIZE
    if batch_count > SET_BATCH_COUNT:
        batch_count = SET_BATCH_COUNT

    # Subsample a list of reads containing no replacement:
    read_list = random.sample(all_reads, batch_count * BATCH_SIZE)

    read_features = []
    for file_i in range(batch_count):
        with pysam.AlignmentFile(f"{outpath}/{filename_prefix}_subsample_{file_i+1}.bam", "wb", header=in_bam.header) as out_bam:
            for i in range(BATCH_SIZE):
                readX = read_list[i + file_i * BATCH_SIZE]
                out_bam.write(readX)
                # Get attributes per read
                soft_clip = sum([y for (x, y) in readX.cigartuples if x == 4])
                nm_edit_distance = readX.get_tag("NM")
                indels = sum([y for (x, y) in readX.cigartuples if x == 1 or x == 2])
                mismatches = nm_edit_distance - indels
                read_length = readX.inferred_length - soft_clip  # exclude soft clip

                read_features.append(
                    [
                        readX.qname,
                        f"{readX.reference_name}",
                        f"{readX.reference_name}_{readX.reference_start}_{readX.reference_end}",
                        readX.mapping_quality,
                        read_length,
                        nm_edit_distance,
                        indels,
                        mismatches,
                        readX.cigartuples,
                        readX.cigarstring,
                        f"batch{file_i}",
                    ]
                )

        pysam.sort("-o", f"{outpath}/{filename_prefix}_subsample_{file_i+1}.sorted.bam", f"{outpath}/{filename_prefix}_subsample_{file_i+1}.bam", )
        pysam.index( f"{outpath}/{filename_prefix}_subsample_{file_i+1}.sorted.bam")
        os.remove(f"{outpath}/{filename_prefix}_subsample_{file_i+1}.bam")
    in_bam.close()
    df = pd.DataFrame(read_features, columns = ['readname',
                                                'chromosome',
                                                'map_pos',
                                                'mapping_quality',
                                                'inferred_read_length',
                                                'nm_tag',
                                                'indels',
                                                'mismatches',
                                                'cigartuples',
                                                'cigarstring',
                                                'batch',
                                                ]
                      )

    df.to_pickle(f"{outpath}/{filename_prefix}.pickle.gz")
    df.to_csv(f"{outpath}/{filename_prefix}.csv")
    print( df.shape[0], "read collected")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-bam-folder", "-i", help="path to bam files to process")
    parser.add_argument("--output-folder", "-o", help="path to write output files")
    args = parser.parse_args()
    outpath = args.output_folder
    filepath = args.input_bam_folder
    # Create folder if not exist
    if not os.path.exists(outpath):
        os.makedirs(outpath)
    for file in os.listdir(filepath):
        if os.path.isdir(file):
            continue
        else:
            # Create error rate per read dataframe
            filename = os.path.join(filepath, file)
            print("Execute:", filename)
            get_error_df(filename, outpath, BATCH_SIZE, SET_BATCH_COUNT)
            print("Output at:", outpath)
