import os
import re
import configparser
import pandas as pd


def parseFlagstat(config):
    flagstat_dir = config.get("FLAG", "bin")
    ff = [fi for fi in os.listdir(flagstat_dir) if "flagstat" in fi]
    print(ff)
    name = re.compile("HM.*_U.*_")
    keys = [re.search(name, pr).group() for pr in ff if re.search(name, pr)]
    print(keys)
    parsed = {key:{} for key in keys}
    print(parsed)
    args = ["genome", "sample", "attr", "total_reads", "mapped"]


    for fi in ff:
        if re.search(name, fi):
            prefix = re.search(name, fi).group()
            parsed[prefix]["genome"] = fi.split("_")[0]
            parsed[prefix]["sample"] = fi.split("_")[1]
            if fi.split("_")[2] == "trimmed":
                parsed[prefix]["attr"] = "NaN"
            else:
                parsed[prefix]["attr"] = fi.split("_")[2]
            with open(os.path.join(flagstat_dir,fi), "r") as fh:
                line1 = fh.readline()
                parsed[prefix]["total_reads"] = int(line1.split()[0])
                fh.readline()
                fh.readline()
                fh.readline()
                line5 = fh.readline()
                parsed[prefix]["mapped"] = int(line5.split()[0])
                parsed[prefix]["percent_mapped"] = parsed[prefix]["mapped"]/parsed[prefix]["total_reads"]

    parsed_df = pd.DataFrame.from_dict(parsed, orient="index")
    path = os.path.join(flagstat_dir, "flagstat_summary.txt")
    parsed_df.to_csv(path)

if __name__ == "__main__":
    config = configparser.ConfigParser()
    config.read("/Users/annasintsova/git_repos/HUTI-RNAseq/analysis/DE/2018-01-16-sync-genomes/lib/config")
    parseFlagstat(config)