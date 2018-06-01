from corda import test_model
from corda import CORDA
import pandas as pd
import cobra
from corda import reaction_confidence

data = "/Users/annasintsova/git_repos/"\
        "HUTI-RNAseq/data/get_homologs_output/C50_S90_e0_/"\
        "run_C50_S90_e0__pan_C50_S90/"\
        "2018-02-26_pangenome_matrix_t0_crossRef.csv"
model = "/Users/annasintsova/git_repos/"\
        "HUTI-RNAseq/analysis/fba/2017-12-16-model-2/"\
        "data/HM54_homologues/HM54_model/HM54_iML1515.xml"

counts = pd.read_csv(data)
HM54_rpkm = counts[["MG1655", "HM54_UR_RPKM", "HM54_UTI_RPKM"]]
HM54_rpkm = HM54_rpkm.dropna(subset=["MG1655"])
HM54_rpkm = HM54_rpkm[HM54_rpkm.MG1655 != "PARALOGS"]
HM54_rpkm.set_index('MG1655', inplace=True)
HM54_rpkm = HM54_rpkm.dropna()


# Let's play
def assignConfidence2(x):
    if x > 1000.0:
        return 2
    # elif x > 500.0:
    #     return 2
    # elif x > 100.0:
    #     return 1
    elif x <= 100.0:
        return -1
    else:
        return 0

def runCORDA(model, conf_dict, file_out, metabolite = ""):
    reaction_conf = {}
    for rx in model.reactions:
        if rx.id == "BIOMASS_Ec_iML1515_core_75p37M":
            reaction_conf[rx.id] = 3
        else:
            reaction_conf[rx.id] = reaction_confidence(rx.gene_reaction_rule, conf_dict)
    if not metabolite:
        opt = CORDA(model, reaction_conf)
        opt.build()
    else:
        opt = CORDA(model, reaction_conf, met_prod=metabolite)
        opt.build()
    print(opt)
    cobra.io.write_sbml_model(opt.model, (file_out+".sbml"))
    with open((file_out+".tab"), "w") as fh:
        for k, used in opt.included.items():
            if used:
                fh.write(opt.model.reactions.get_by_id(k).id +
                         "\t" + opt.model.reactions.get_by_id(k).name + "\n")

    return opt

confidences = {'UR_conf':list(map(assignConfidence2, HM54_rpkm.HM54_UR_RPKM.values)),
               'UTI_conf':list(map(assignConfidence2, HM54_rpkm.HM54_UTI_RPKM.values))}

HM54 = pd.DataFrame(confidences,index=HM54_rpkm.index)
ur_conf = HM54.UR_conf.to_dict()
uti_conf = HM54.UTI_conf.to_dict()

base = cobra.io.read_sbml_model(model)

runCORDA(base, ur_conf, "data/2018-03-26-HM54-UR-BM-CORDA")
#runCORDA(base, uti_conf, "data/2018-03-26-HM54-UTI-BM-CORDA")
#runCORDA(base, ur_conf, "data/2018-03-26-HM54-UR-ATPM-CORDA", metabolite)
#runCORDA(base, uti_conf, "data/2018-03-26-HM54-UTI-ATPM-CORDA", metabolite)
