
import pandas as pd
import os
import numpy as np
import scipy.stats
import datetime as dt
import cobra
from corda import reaction_confidence
from corda import CORDA
import sys

#########################################################################################################
# Original functions, working with RPKMs directly

def assignConfidence(x, y = 500, z = 10):

    if x > y:
        return 2
    elif x < z:
        return - 1
    else:
        return 0


def binRPKMs(data, strain):
    """
    Take in counts, select rpkms for correct strain, assign them confidence,
    return confidence for ur and uti conditions (2 dictionaries)
    """

    counts = pd.read_csv(data)
    try:
        rpkm = counts[["MG1655", "{}_UR_RPKM".format(strain), "{}_UTI_RPKM".format(strain)]]
    except:
        strain = "HM0" + strain.split("HM")[1]
        rpkm = counts[["MG1655", "{}_UR_RPKM".format(strain), "{}_UTI_RPKM".format(strain)]]
    finally:
        rpkm = rpkm.dropna(subset=["MG1655"])
        rpkm = rpkm[rpkm.MG1655 != "PARALOGS"]
        rpkm.set_index('MG1655', inplace=True)
        rpkm = rpkm.dropna()
        print("UTI")
        print(rpkm[rpkm["{}_UTI_RPKM".format(strain)] > 1].shape)
        print(rpkm[rpkm["{}_UTI_RPKM".format(strain)] < -1.5].shape)
        print("URINE")
        print(rpkm[rpkm["{}_UR_RPKM".format(strain)] > 1].shape)
        print(rpkm[rpkm["{}_UR_RPKM".format(strain)] < -1.5].shape)
        confidences = {'UR_conf': list(map(assignConfidence, rpkm["{}_UR_RPKM".format(strain)].values)),
                       'UTI_conf': list(map(assignConfidence, rpkm["{}_UTI_RPKM".format(strain)].values))}
        confidence_df = pd.DataFrame(confidences, index=rpkm.index)
        ur_conf = confidence_df.UR_conf.to_dict()
        uti_conf = confidence_df.UTI_conf.to_dict()
        return ur_conf, uti_conf
###########################################################################################################


def invnorm(x):
    """
    Inverse rank normalize a list of numbers

    :param x: a list of numbers
    :return: a normalized list of numbers

    """
    return scipy.stats.norm.ppf((x.rank() -0.5)/x.count())


def assignNormConfidence(x):

    if x > 1.5:
        return 2
    elif x < - 1.5:
        return - 1
    else:
        return 0


def getModel(strain, model_dir):

    model_path = os.path.join(model_dir,
                              "{}_homologues/{}_model/{}_iML1515.xml".format(strain, strain, strain))

    return cobra.io.read_sbml_model(model_path)


def binNormRPKMs(data, strain):

    """
    Take in counts, select rpkms for correct strain, normalaize them, assign them confidence,
    return confidence for ur and uti conditions (2 dictionaries)

    :param data: path to counts file
    :return: two dictionaries, one for urine, one for uti

    """
    zero_strains = ["HM1", "HM3", "HM6", "HM7"]

    new_strain = "HM0" + strain.split("HM")[-1] if strain in zero_strains else strain

    counts = pd.read_csv(data)
    rpkm = counts[["MG1655", "{}_UR_RPKM".format(new_strain), "{}_UTI_RPKM".format(new_strain)]]
    rpkm = rpkm.dropna(subset=["MG1655"])
    rpkm = rpkm[rpkm.MG1655 != "PARALOGS"]
    rpkm.set_index('MG1655', inplace=True)
    rpkm = rpkm.dropna()
    rpkm["{}_UR".format(strain)] = invnorm(rpkm["{}_UR_RPKM".format(new_strain)])
    rpkm["{}_UTI".format(strain)] = invnorm(rpkm["{}_UTI_RPKM".format(new_strain)])
    confidences = {'UR_conf': list(map(assignNormConfidence,
                                       rpkm["{}_UR".format(strain)].values)),
                   'UTI_conf': list(map(assignNormConfidence,
                                        rpkm["{}_UTI".format(strain)].values))}
    confidence_df = pd.DataFrame(confidences, index=rpkm.index)
    ur_conf = confidence_df.UR_conf.to_dict()
    uti_conf = confidence_df.UTI_conf.to_dict()

    print("In URINE:")
    print("Confidence of 2: {}".format(list(ur_conf.values()).count(2)))
    print("Confidence of -1: {}".format(list(ur_conf.values()).count(-1)))
    print("In UTI:")
    print("Confidence of 2: {}".format(list(uti_conf.values()).count(2)))
    print("Confidence of -1: {}".format(list(uti_conf.values()).count(-1)))

    return ur_conf, uti_conf





def runCORDA(model, conf_dict, file_out, metabolite = ""):
    """
    
    """
    print(model.id)
    print(len(model.reactions))
    reaction_conf = {}
    for rx in model.reactions:
        if rx.id == "BIOMASS_Ec_iML1515_core_75p37M":
            reaction_conf[rx.id] = 3
        else:
            reaction_conf[rx.id] = reaction_confidence(rx.gene_reaction_rule, conf_dict)
    with open(file_out + ".csv", "w") as fd:
        for k, v in reaction_conf.items():
            fd.write(str(k) + ',' + str(v) + "\n")

    if not metabolite:
        opt = CORDA(model, reaction_conf)
        opt.build()
    else:
        opt = CORDA(model, reaction_conf, met_prod=metabolite)
        opt.build()
    print(opt)
   
    with open(file_out+".tab", "w") as fh:
        for k, used in opt.included.items():
            if used:
                fh.write(opt.model.reactions.get_by_id(k).id +
                         "\t" + opt.model.reactions.get_by_id(k).name + "\t" +
                         opt.model.reactions.get_by_id(k).gene_reaction_rule + "\n")

    return opt


def buildEnvModel(included_dict, model, strain, condition, output_dir):

    today = dt.datetime.today().strftime("%Y-%m-%d")
    print(strain)
    print(condition)
    print(len(model.reactions))

    all_reactions = list(model.reactions)

    for rx in all_reactions:
        if not included_dict[rx.id]:
            model.reactions.get_by_id(rx.id).remove_from_model()
        
    print(len(model.reactions))
    name_suffix = model.id.split("_")[-1]
    model.id = "{}_{}_{}".format(strain, condition, name_suffix)
    cobra.io.write_sbml_model(model,
                              os.path.join(output_dir, (today + "_"+ strain +
                                                        "_" + condition +".sbml")))

    return model



def cordaAnalysis(strains, data, output_dir, model_dir):

    today = dt.datetime.today().strftime("%Y-%m-%d")
    all_envs ={}
    all_models = []
    for strain in strains:
        print(strain)
        # bin
        print("Binning data")
        UR, UTI = binNormRPKMs(data, strain)
        # run CORDA
        # for UTI
        print("Running CORDA for UTI")
        UTI_key = "{}_UTI".format(strain)
        uti_env = runCORDA(getModel(strain, model_dir),
                           UTI, os.path.join(output_dir, today+"_"+UTI_key))
        # for urine
        print("Running CORDA for urine")
        UR_key = "{}_UR".format(strain)
        urine_env = runCORDA(getModel(strain, model_dir),
                             UR, os.path.join(output_dir, today+"_"+UR_key))

        all_envs[UR_key] = urine_env
        all_envs[UTI_key] = uti_env

    for key, val in all_envs.items():
        st = key.split("_")[0]
        condition = key.split("_")[1]
        print("Building Model")
        m = buildEnvModel(val.included, getModel(st, model_dir), st, condition, output_dir)
        all_models.append(m)
        print("done")

    return all_envs, all_models


if __name__ == "__main__":

    # Best Models:
    # HM56, HM14, HM43, HM54, HM86

    data = "/Users/annasintsova/git_repos/" \
           "HUTI-RNAseq/data/get_homologs_output/C50_S90_e0_/" \
           "run_C50_S90_e0__pan_C50_S90/" \
           "2018-02-26_pangenome_matrix_t0_crossRef.csv"

    model_dir = "/Users/annasintsova/git_repos/"   \
                "HUTI-RNAseq/analysis/fba/2017-12-16-model-2/data"

    strains = ["HM1", "HM6", "HM14", "HM43", "HM54", "HM56", "HM68", "HM86", "HM3"]

    #strains = ["HM14","HM43", "HM54", "HM56", "HM86"]
    #strains = ["HM3"]
    output_dir = "/Users/annasintsova/git_repos/HUTI-RNAseq/analysis/fba/2018-03-21-clinical-strain-metabolism/data/2018-04-11-CORDA-analysis"

    all_envs, all_models = cordaAnalysis(strains, data, output_dir, model_dir)


