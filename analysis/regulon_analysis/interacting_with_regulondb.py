import sqlite3
import pandas as pd
import os
import sys
sys.path.append('/Users/annasintsova/git_repos/HUTI-RNAseq/analysis/methods')
import helpers


def get_all_regulators(db):
    """Returns a list of regulators"""
    conn = sqlite3.connect(db)
    cur = conn.cursor()
    cur.execute("SELECT transcription_factor_id "
                "FROM TRANSCRIPTION_FACTOR ")
    rows = cur.fetchall()
    conn.close()
    return [r[0] for r in rows]


def get_all_regulators_iterator(db):
    """Returns one regulator at a time """
    conn = sqlite3.connect(db)
    cur = conn.cursor()
    for row in cur.execute("SELECT transcription_factor_id "
                           "FROM TRANSCRIPTION_FACTOR "):
        yield row[0]
    conn.close()


def get_bnum(db, object_id=""):
    """given id try to find bnum"""
    conn = sqlite3.connect(db)
    cur = conn.cursor()
    cur.execute("SELECT ext_reference_id "
                "FROM OBJECT_EXTERNAL_DB_LINK "
                "WHERE object_id=? ",
                (object_id,))
    rows = cur.fetchall()
    conn.close()
    if rows:
        return rows[0][0]
    return object_id


def get_regulon(db, gene_name_regulator="", tf_id_regulator=""):
    """Returns list of tuples (regulator_name,
    regulated_name, regulated_id, repressed/activated)"""
    conn = sqlite3.connect(db)
    cur = conn.cursor()
    cur.execute("SELECT gene_name_regulator, gene_name_regulated, "
                "gene_id_regulated, generegulation_function "
                "FROM GENEREGULATION_TMP "
                "WHERE gene_name_regulator=? OR tf_id_regulator=?",
                (gene_name_regulator, tf_id_regulator))
    regulon = cur.fetchall()
    conn.close()
    reg_network = []
    for reg in regulon:
        eck = reg[2]
        bnum = get_bnum(db, object_id=eck)
        reg_network.append((reg[0], reg[1], bnum, reg[3]))  # Should I make this a dictionary?
    return reg_network


def get_all_regulons(db, filename):
    regulators = get_all_regulators(db)
    pd_list = []
    for regulator in regulators:
        regulon = get_regulon(db, tf_id_regulator=regulator)
        pd_list.append(pd.DataFrame(regulon))
    df = pd.concat(pd_list).reset_index(drop=True)
    df.to_csv(filename)
    return df


def create_regulon_csv():
    config = os.path.join(os.path.dirname(os.path.realpath(__file__)), "config")
    config_dict = helpers.process_config(config)
    filename = config_dict["out_dir"]["regulon_csv"]
    rdb = config_dict["db"]["path"]
    rn = get_all_regulons(rdb, filename)
    return rn


if __name__ == "__main__":
    print(create_regulon_csv().head())


