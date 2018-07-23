import sqlite3
import pandas as pd


def view(db):
    conn = sqlite3.connect(db)
    cur = conn.cursor()
    cur.execute("SELECT * FROM SIGMA_TMP")
    rows = cur.fetchall()
    conn.close()
    return rows


def get_all_regulators(db):
    conn = sqlite3.connect(db)
    cur = conn.cursor()
    cur.execute("SELECT transcription_factor_id "
                "FROM TRANSCRIPTION_FACTOR ")
    rows = cur.fetchall()
    conn.close()
    return [r[0] for r in rows]  # returns list of regulators


def get_bnum(db, object_id=""):  # given id try to find bnum
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
        reg_network.append((reg[0], reg[1], bnum, reg[3]))

    return reg_network  # returns list of tuples (regulated_name, id, reprssed/activated)


def get_all_regulons(db):
    regulators = get_all_regulators(db)
    # t_r = [regulators[10], regulators[22]]
    pd_list = []
    for regulator in regulators:
        regulon = get_regulon(db, tf_id_regulator=regulator)
        pd_list.append(pd.DataFrame(regulon))

    df = pd.concat(pd_list).reset_index(drop=True)
    df.to_csv("/Users/annasintsova/git_repos/HUTI-RNAseq/results/regulon_analysis/regulatory_network.csv")
    return df


if __name__ == "__main__":
    rdb = "/Users/annasintsova/git_repos/HUTI-RNAseq/data/regulondb9.4.db"
    rn = get_all_regulons(rdb)
    print(rn.head())
