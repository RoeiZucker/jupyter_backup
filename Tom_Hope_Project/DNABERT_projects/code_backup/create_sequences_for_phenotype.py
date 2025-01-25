import os.path

import pandas as pd
from basic_reader import read_range

def apply_read_range(row):
    chrom = row["GRCh38Chromosome"]
    # if "-" in row["GRCh38Location"]:
    #     pos = int(row["GRCh38Location"].split("-")[0])
    # else:
    pos = int(row["GRCh38Location"])
    return read_range(chrom,pos - BATCH_SIZE // 2, pos + BATCH_SIZE // 2)

def join_neighbors(df):
    dfs = []
    for chrom in df["GRCh38Chromosome"].unique():
        curr_df = df[df["GRCh38Chromosome"] == chrom]
        lst = (curr_df.values.tolist())
        counter = 0
        for i in range(len(lst)):
            if "-" in lst[i][1]:
                splited = lst[i][1].split(" - ")
                lst[i][1] = str((int(splited[0]) + int(splited[1]) )//2)
        while counter < len(lst) - 1:
            if int(lst[counter + 1][1]) - int(lst[counter][1]) < COMBINING_DISTANCE:
                lst[counter][1] = (int(lst[counter][1]) + int(lst[counter + 1][1])) // 2
                lst.remove(lst[counter + 1])
                # counter = 0
                continue

            counter += 1
        new_df = pd.DataFrame(lst,columns=["GRCh38Chromosome","GRCh38Location","key"])
        dfs.append(new_df)
    return pd.concat(dfs)

    # for chrom in df["GRCh38Chromosome"].unique():
    #     curr_df =
#     split to different chromosomes
#     run iteretavly on each df, if one less than COMBINING_DISTANCE than the other, remove them both and add the avarage
#     start the loop again
    pass

PHEN = "Amyotrophic lateral sclerosis"
CLINVAR_RESULTS_PATH = "C:\\Users\\USER\\Downloads\\clinvar_result_Multiple sclerosis.txt"
RESULTS_PATH = "C:/Users/USER/Dropbox/study/tom_lab/GWAS_REP/seq_files"

BATCH_SIZE = 512
COMBINING_DISTANCE = 206
df = pd.read_csv(CLINVAR_RESULTS_PATH,sep="\t")
# print(df.columns)
df = df[df["Condition(s)"].str.contains(PHEN)]
df_benign = df[(df["Germline classification"] == "Benign") | (df["Germline classification"] == "Likely benign")]
df_pathogenic = df[(df["Germline classification"] == "Pathogenic") | (df["Germline classification"] == "Likely pathogenic")]
# print(df_pathogenic[["GRCh38Chromosome","GRCh38Location"]])
# print(df_pathogenic[["GRCh38Chromosome","GRCh38Location"]])
df_pathogenic = df_pathogenic[["GRCh38Chromosome","GRCh38Location"]]
df_benign = df_benign[["GRCh38Chromosome","GRCh38Location"]]
df_pathogenic["val"] = 1
df_benign["val"] = 0
df_pathogenic = join_neighbors(df_pathogenic)
df_benign = join_neighbors(df_benign)
df_pathogenic["seq"] = (df_pathogenic[["GRCh38Chromosome","GRCh38Location"]].apply(apply_read_range,axis=1))
df_benign["seq"] = (df_benign[["GRCh38Chromosome","GRCh38Location"]].apply(apply_read_range,axis=1))

df_result = pd.concat([df_benign,df_pathogenic])
df_result = df_result.dropna()
print(df_result)
df_result.to_csv(os.path.join(RESULTS_PATH,PHEN + "_reduced.csv"))