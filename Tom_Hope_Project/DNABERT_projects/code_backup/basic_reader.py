import os

REFERENCE_FILE = "C:/Users/USER/Downloads/ncbi_dataset/data/GCF_000001405.40/GCF_000001405.40_GRCh38.p14_genomic.fna"
CREATED_FILES_PATH = "C:/Users/USER/Dropbox/study/tom_lab/GWAS_REP"
POS_DICT ={
    1:-67,
    4:-67,
    6:-67
}




def read_range(chrom,position,length):
    counter = -67
    if chrom in POS_DICT:
        counter = POS_DICT[chrom]
    sequences = []
    end = position + length
    file = os.path.join(CREATED_FILES_PATH,str(chrom) +".fna")
    with open(file) as f:
        for l in f:
            line = l[:-1]
            if counter + len(line)< position:
                counter += len(line)
                continue
            if counter < position:
                to_advance = position - counter
                if end < counter + len(line):
                    # in case all fits in first line
                    return line[to_advance:to_advance+length]
                else:
                    # in case it ends in a future line
                    sequences.append(line[to_advance:])
                    counter += len(line)
                continue
            if end > counter + len(line):
                sequences.append(line)
                counter += len(line)
                continue
#             reached end
            sequences.append(line[: end - counter + 1])
            return "".join(sequences)









if __name__ == "__main__":
    print(read_range(5,318000,100))
