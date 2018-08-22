import os


def process_flagstat(flagstat_directory):

    fstats = [os.path.join(flagstat_directory, f) for f in os.listdir(flagstat_directory)]
    flag_summary = {}
    for fstat in fstats:
        if "flagstat.txt" in fstat:
            with open(fstat, "r") as fh:
                sample = os.path.basename(fstat).split('_')[0]
                summary = fh.readlines()
                flag_summary[sample] = (summary[0].split()[0], summary[4].split()[0])
    return flag_summary

