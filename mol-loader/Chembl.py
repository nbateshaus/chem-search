"""
Encapsulate interaction with ChEMBLdb.
"""

import psycopg2

from rdkit_utils import rdkit_standardize, rdkit_descriptors, rdkit_smiles, rdkit_mol_from_smiles, rdkit_fps_from_mol


class Chembl:
    def __init__(self, dsn, limit=None):
        """
        Encapsulate interaction with ChEMBLdb.

        Keyword arguments:
        dsn - Connection string for PostgreSQL, following http://www.postgresql.org/docs/current/static/libpq-connect.html#LIBPQ-CONNSTRING
        limit - integer limit on number of molecules to generate
        Other keyword arguments will be passed to psycopg2 for database connection.
        """
        self.dsn = dsn
        self.limit = limit

    def __iter__(self):
        """Yield molecules"""
        q = self.__base_query
        if self.limit is not None:
            q += "LIMIT {0}".format(self.limit) # TODO: use SQL parameters
        cur = psycopg2.connect(self.dsn).cursor()
        cur.execute(q)
        cols = ['ChEMBL_' + col[0] for col in cur.description]
        row = cur.fetchone()
        while row is not None:
            d = dict(zip(cols, row))
            d["id"] = d["ChEMBL_chembl_id"]
            d["URL"] = "https://www.ebi.ac.uk/chembl/compound/inspect/" + d["ChEMBL_chembl_id"]
            s = d["ChEMBL_smiles"]
            if s is not None:
                smiles = [s]
                mol = rdkit_standardize(rdkit_mol_from_smiles(s))
                if mol is not None:
                    fps, fps_bits = rdkit_fps_from_mol(mol)
                    d['RDKit_Fingerprint'] = fps
                    d['RDKit_Fingerprint_bits'] = fps_bits

                    rs = rdkit_smiles(mol)
                    if rs is not None:
                        d['RDKit_SMILES'] = rs
                        smiles.append(rs)
                    descs = rdkit_descriptors(mol)
                    for name in descs:
                        d['RDKit_' + name] = descs[name]
                d["SMILES"] = list(set(smiles))
            yield d
            row = cur.fetchone()

    __base_query = '''
WITH synonyms_agg AS (
    SELECT
        molregno,
        array_agg(synonyms) AS synonyms
    FROM molecule_synonyms
    GROUP BY molregno
  ), alerts_agg AS (
    SELECT
        molregno,
        array_agg(alert_name) AS alerts
    FROM compound_structural_alerts csa
      JOIN structural_alerts sa ON (sa.alert_id = csa.alert_id)
    GROUP BY molregno
  )
SELECT
  md.molregno                       AS "molregno",
  CASE md.availability_type
    -- -1 appears in the data; let it go to NULL
    WHEN 0 THEN 'discontinued'
    WHEN 1 THEN 'prescription only'
    WHEN 2 THEN 'over the counter'
    END                             AS "availability_type",
  (md.black_box_warning <> 0)       AS "black_box_warning",
  md.chebi_par_id                   AS "chebi_id",
  md.chembl_id                      AS "chembl_id",
  CASE md.chirality
    -- -1 appears in the data; let it go to NULL
    WHEN 0 THEN 'racemic mixture'
    WHEN 1 THEN 'single stereoisomer'
    WHEN 2 THEN 'achiral molecule'
    END                             AS "chirality",
  (md.dosed_ingredient <> 0)        AS "dosed_ingredient",
  md.first_approval                 AS "first_approval",
  (md.first_in_class <> 0)          AS "first_in_class",
  md.indication_class               AS "indication_class",
  (md.inorganic_flag <> 0)          AS "inorganic",
  md.max_phase                      AS "max_phase",
  md.molecule_type                  AS "molecule_type",
  (md.natural_product <> 0)         AS "natural_product",
  (md.oral <> 0)                    AS "oral",
  (md.parenteral <> 0)              AS "parenteral",
  (md.polymer_flag <> 0)            AS "polymer",
  md.pref_name                      AS "preferred_name",
  (md.prodrug <> 0)                 AS "prodrug",
  md.structure_type                 AS "structure_type",
  (md.therapeutic_flag <> 0)        AS "therapeutic",
  (md.topical <> 0)                 AS "topical",
  md.usan_stem_definition           AS "usan_stem_definition",
  md.usan_stem                      AS "usan_stem",
  md.usan_substem                   AS "usan_substem",
  md.usan_year                      AS "usan_year",
  cp.mw_freebase                    AS "mw_freebase",
  cp.alogp                          AS "alogp",
  cp.hba                            AS "hba",
  cp.hbd                            AS "hbd",
  cp.psa                            AS "psa",
  cp.rtb                            AS "rtb",
  CASE cp.ro3_pass
    WHEN 'Y' THEN true
    WHEN 'N' THEN false
  END                               AS "ro3_pass",
  cp.num_ro5_violations             AS "num_ro5_violations",
  CASE cp.med_chem_friendly
    WHEN 'Y' THEN true
    WHEN 'N' THEN false
  END                               AS "med_chem_friendly",
  cp.acd_most_apka                  AS "acd_most_apka",
  cp.acd_most_bpka                  AS "acd_most_bpka",
  cp.acd_logp                       AS "acd_logp",
  cp.acd_logd                       AS "acd_logd",
  cp.molecular_species              AS "molecular_species",
  cp.full_mwt                       AS "full_mwt",
  cp.aromatic_rings                 AS "aromatic_rings",
  cp.heavy_atoms                    AS "heavy_atom_count",
  cp.num_alerts                     AS "num_alerts",
  cp.qed_weighted                   AS "qed_weighted",
  cp.mw_monoisotopic                AS "monoisotopic_weight",
  cp.full_molformula                AS "full_molformula",
  cp.hba_lipinski                   AS "hba_lipinski",
  cp.hbd_lipinski                   AS "hbd_lipinski",
  cp.num_lipinski_ro5_violations    AS "num_lipinski_ro5_violations",
  cs.standard_inchi                 AS "inchi",
  cs.standard_inchi_key             AS "inchi_key",
  cs.canonical_smiles               AS "smiles",
  bt.helm_notation                  AS "helm_notation",
  syn.synonyms                      AS "synonyms",
  alrt.alerts                       AS "alerts"
FROM molecule_dictionary md
LEFT JOIN compound_properties cp ON md.molregno = cp.molregno
LEFT JOIN compound_structures cs ON md.molregno = cs.molregno
LEFT JOIN biotherapeutics bt     ON md.molregno = bt.molregno
LEFT JOIN synonyms_agg syn       ON md.molregno = syn.molregno
LEFT JOIN alerts_agg alrt        ON md.molregno = alrt.molregno
'''


def test_iter():
    c = Chembl("dbname=chembl_20", limit=10)
    for m in c:
        print(m)

if __name__ == '__main__':
    test_iter()
