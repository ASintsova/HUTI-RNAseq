
# coding: utf-8


# Run *in vitro* growth simulations, as described in Monk et al. 10.1073/pnas.1307797110

import cobra
import os
import libsbml
import pandas as pd
import datetime as dt
import scipy
import itertools


# ## For specific growth conditions simulations:
# ### For anaerobic condtions:
# EX_o2_e (0, 1000)
# ### For carbon source profilling:
# EX_glc__D_e (0, 1000)
# ### For nitrogen source profilling:
# EX_nh4_e (0, 1000)
# ### For phosphorous source profilling:
# EX_pi_e (0, 1000)
# ### For sulfur source profilling:
# EX_so4_e (0, 1000)


#################################################################################


# Conditions and settings taken from Monk PNAS paper

def defineM9():
    """
    Reads in M9 conditions from file, returns dataframe of M9 conditons, reactions as index
    """
    
    m9_cond_file = "/Users/annasintsova/git_repos/HUTI-RNAseq/"                "analysis/fba/2018-03-21-clinical-strain-metabolism/data/"                "M9_conditions.csv"
    m9 = pd.read_csv(m9_cond_file, index_col=0)
    m9_reactions = list(m9.index)
    m9_reactions[19] = "EX_glc__D_e"
    m9.index = m9_reactions
    return m9


def setToM9(model, m9_settings):
    m9_reactions = list(m9_settings.index)
    original_settings = {k: [model.reactions.get_by_id(k).name, model.reactions.get_by_id(k).lower_bound,
                             model.reactions.get_by_id(k).upper_bound] for k in m9_reactions}

    # print(original_settings)

    for rx in original_settings.keys():
        new_lower_bound = m9_settings.loc[rx]["lower_bound"]
        new_upper_bound = m9_settings.loc[rx]["upper_bound"]
        model.reactions.get_by_id(rx).lower_bound = new_lower_bound
        model.reactions.get_by_id(rx).uppper_bound = new_upper_bound
    new_settings = {k: [model.reactions.get_by_id(k).name, model.reactions.get_by_id(k).lower_bound,
                        model.reactions.get_by_id(k).upper_bound] for k in m9_reactions}
    # print(new_settings)
    return model


def getReactionIds(model):
    conditions = []
    r_ids = {}
    conditions_with_growth_file = "/Users/annasintsova/git_repos/HUTI-RNAseq/"                    "analysis/fba/2018-03-21-clinical-strain-metabolism/data/"                    "conditionsWithGrowthMonk.csv"
    growth_cond = pd.read_csv(conditions_with_growth_file, index_col=0)

    for i in growth_cond.index:
        rnames = set([c.split("_")[2].strip() for c in growth_cond.index])
        for rx in model.reactions:
            if rx.name in rnames:
                r_ids[rx.name] = rx.id
                rnames.remove(rx.name)

    rx_info_file = "/Users/annasintsova/git_repos/HUTI-RNAseq/"                    "analysis/fba/2018-03-21-clinical-strain-metabolism/data/"                    "reactionInformationMonk.csv"
    rx_info = pd.read_csv(rx_info_file, index_col=1)

    df = rx_info.loc[rnames]
    for n, r in zip(df.index, df["Reaction Abbreviation"]):

        try:
            if model.reactions.get_by_id(r):
                r_ids[n] = r
                rnames.remove(n)
        except:
            pass
    # Not ideal, but I don't think this matters a huge deal, got 207 out of 221
    return r_ids


def setConditions(model, reaction_ids):
    conditions_with_growth_file = "/Users/annasintsova/git_repos/HUTI-RNAseq/"                    "analysis/fba/2018-03-21-clinical-strain-metabolism/data/"                    "conditionsWithGrowthMonk.csv"

    growth_cond = pd.read_csv(conditions_with_growth_file, index_col=0)

    conditions = []
    for c in growth_cond.index:
        if c.split("_")[2] in reaction_ids.keys():
            element = c.split("_")[0]
            ox = True if c.split("_")[1] == "AERO" else False
            reaction = reaction_ids[c.split("_")[2]]
            conditions.append((element, reaction, ox))
    return conditions


def growthSimulation2(models, out_dir):
    """
    Return dataframe with all models growth under all conditions
    For each model run FBA with all specified conditions: seperate function

    """
    all_growth = []
    model_index = []
    for model in models:
        conditions = setConditions(model, getReactionIds(model))
        growth = doesItGrow(model, conditions)
        all_growth.append(growth)
        model_index.append(model.id)
    growth_simulation = pd.DataFrame(all_growth)
    growth_simulation.index = model_index
    growth_simulation.to_csv(os.path.join(out_dir,
                                          "{}_growth_simulations.csv".format(dt.datetime.today().strftime("%Y-%m-%d"))))
    return growth_simulation
###################################################################################



def growthSimulation(models, conditions, out_dir):
    """
    Return dataframe with all models growth under all conditions
    For each model run FBA with all specified conditions: seperate function
    
    """
    all_growth = []
    model_index = []
    for model in models:
        print("Working on: {}".format(model.id))
        growth = doesItGrow(model, conditions)
        all_growth.append(growth)
        model_index.append(model.id)
    print("Putting data together into a dataframe")
    growth_simulation = pd.DataFrame(all_growth)
    growth_simulation.index = model_index
    
    print("Removing 'no growth' conditions")
    growing = growth_simulation.loc[:,(growth_simulation.sum(axis=0) != 0)]
    
    growing.to_csv(os.path.join(out_dir, 
                    "{}_growth_simulations.csv".format(dt.datetime.today().strftime("%Y-%m-%d"))))
    print("Done")
    return growing


def doesItGrow(model, conditions):
    """
    Return dictionary of {condition:(no)growth} for a model
    
    1. Set model to M9 conditions
    2. Go through list of conditions run fba, collect data on growth
    
        simulation_condition: list of tuples: [(element (C, N, P, S), reaction.id, oxygen (True/False))
        
    3. restore original settings
    
    """
    # 1
    m9_settings = defineM9()
    print("Setting M9 conditions")
    model = setToM9(model, m9_settings)
    
  
    growth = {}
    for condition in conditions:
        
        air_status = "AERO" if condition[2] else "ANAERO"
        cname = "{}_{}_{}".format(condition[0], air_status, condition[1])
        print("Testing: {}".format(cname))
        
        #2
        if condition[0] == "C":
            element = model.reactions.EX_glc__D_e
        elif condition[0] == "N":
            element = model.reactions.EX_nh4_e
        elif condition[0] == "P":
            element = model.reactions.EX_pi_e
        elif condition[0] == "S":
            element = model.reactions.EX_so4_e
        else:
            print("wrong condition format")
            
        original_element_lb = element.lower_bound
        original_element_up = element.uppper_bound
        element.lower_bound = 0
        element.upper_bound = 1000
        
        r = model.reactions.get_by_id(condition[1])
        original_lower_bound = r.lower_bound
        r.lower_bound = -1000
            
        
        ox = model.reactions.EX_o2_e
        if not condition[2]:
            ox.lower_bound = 0
        
        obj = round(model.optimize().objective_value,4)
        growth[cname] = obj if obj > 0 else 0
        
        #3 
        element.lower_bound = original_element_lb
        element.upper_bound = original_element_up
        r.lower_bound = original_lower_bound
        ox.lower_bound = -1000.0
            
    return growth


def allCond(model):

    """
    generate all possible conditions using all exchange reactions from model

    :param model: metabolic model
    :return: list of tuples (element, reaction_id, oxygen(true/flase))

    """
    
    elements = ["C", "N", "P", "S"]
    reactions = [rx.id for rx in model.reactions if "EX_" in rx.id]
    ox = [True, False]
    conditions = list(itertools.product(elements, reactions, ox))
    return conditions
    

if __name__ == "__main__":

    # Run Simulation
    model_path = "/Users/annasintsova/git_repos/HUTI-RNAseq/analysis/fba/2018-03-21-clinical-strain-metabolism/data/models"
    models = [cobra.io.read_sbml_model(os.path.join(model_path, mp)) for mp in os.listdir(model_path)]
    out_dir = "/Users/annasintsova/git_repos/HUTI-RNAseq/analysis/fba/2018-03-21-clinical-strain-metabolism/data/"
    conditions = allCond(models[0])
    df = growthSimulation(models, conditions, out_dir)




