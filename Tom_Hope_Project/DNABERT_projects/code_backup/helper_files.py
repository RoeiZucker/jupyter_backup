import os

REFERENCE_FILE = "C:/Users/USER/Downloads/ncbi_dataset/data/GCF_000001405.40/GCF_000001405.40_GRCh38.p14_genomic.fna"
DESTINATION_PATH = "C:/Users/USER/Dropbox/study/tom_lab/GWAS_REP"


def split_file_to_chromosomes(file):
    with open(file) as f:
        count = 0
        current_chr = ""
        in_correct_zone = False
        curr_lines = []
        for line in f:
            if line.startswith(">") and line.endswith("Primary Assembly\n"):
                if current_chr != "":
                    write_chr_file(curr_lines, current_chr)
                # start_new_file
                current_chr = get_current_chr(line)
                print(current_chr)
                curr_lines = []
                in_correct_zone = True
            if line.startswith(">") and not line.endswith("Primary Assembly\n"):
                in_correct_zone = False
            if in_correct_zone:
                curr_lines.append(line)
            count += 1
    write_chr_file(curr_lines, current_chr)
    print(count)
    pass


def get_current_chr(line):
    val = (line.find("Homo sapiens chromosome "))
    current_chr = (
        line[val + len("Homo sapiens chromosome "):val + 2 + len("Homo sapiens chromosome ")].replace(",", ""))
    return current_chr


def write_chr_file(curr_lines, current_chr):
    with open(os.path.join(DESTINATION_PATH, current_chr + ".fna"), "w") as chr_file:
        chr_file.writelines(curr_lines)


if __name__ == "__main__":
    split_file_to_chromosomes(REFERENCE_FILE)
