from itertools import product

ALL = ["HM01", "HM03","HM06", "HM07",
              "HM14", "HM17", "HM43", "HM54", "HM56",
              "HM57", "HM60", 'HM68',"HM66", "HM86"]
COND = ["UR", "UTI"]

def generateSampleNames(genomes = ALL, sample = COND, reseqed = True):
    reseq = ["HM06", "HM07", "HM43", "HM57", "HM60", 'HM68']

    sample_names = list(product(genomes, sample))
    if reseqed:
        attr = ["UTI_seq1", "UTI_seq2"]
        additional_samples = [ge for ge in genomes if ge in reseq]
        additional_names = list(product(additional_samples, attr))
        sample_names += additional_names
    names = ["_".join(s) for s in sample_names]
    return names