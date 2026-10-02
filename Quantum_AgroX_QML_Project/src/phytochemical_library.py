"""
Quantum_AgroX – phytochemical_library.py
Curated SMILES + resolver shared by AgroPhytoX -> AgroDockX -> AgroReport.

Strategy:
  1. Curated RDKit-validated SMILES for common benchmark phytochemicals.
  2. Case-insensitive / synonym lookup.
  3. Live PubChem PUG-REST fallback (optional, needs internet).
  4. Manual SMILES entry fallback (for obscure GC-MS names).
"""

from __future__ import annotations

import re
import urllib.parse

# ---------------------------------------------------------------------------
# Curated SMILES (all validated with RDKit MolFromSmiles).
# For complex steroids the scaffold SMILES is representative; exact
# stereochemistry of the GC-MS-reported derivative should be curated
# from PubChem / Supplementary File 1 before production docking.
# ---------------------------------------------------------------------------
CURATED_SMILES: dict[str, str] = {
    # Sesquiterpenes / monoterpenoids
    "alpha-caryophyllene": "CC1=CCCC(=CCCC1)C",  # humulene (= alpha-caryophyllene)
    "humulene": "CC1=CCCC(=CCCC1)C",
    "alpha-cedrol": "CC1CCC2C1C1CCCC1(C)C2O",
    "cedrol": "CC1CCC2C1C1CCCC1(C)C2O",
    "cis-farnesol": "CC(=CCCC(=CCCC(=CCO)C)C)C",
    "farnesol": "CC(=CCCC(=CCCC(=CCO)C)C)C",
    "valencene": "CC1=CCC2CCC(C)C2C1=C(C)C",
    "copaene": "CC1=CCC2C3CC(C2C1)C3(C)C",
    "alpha-copaene": "CC1=CCC2C3CC(C2C1)C3(C)C",
    "curzerene": "C=CC1=C(C)C2=C(O1)C(C)CC2",
    "dendrolasin": "CC(=CCCC(=CCOC1=CC=CO1)C)C",
    "torreyol": "CC1CCC(C)C2CCC(C)(O)C2C1",
    "caryophyllene": "CC1=CCCC(C)(C)C2CC1C2C=C",  # beta-caryophyllene scaffold
    # Diterpenoids / fatty acids / esters
    "phytol": "CC(C)CCCC(C)CCCC(C)CCCC(=CCO)C",
    "hexadecanoic acid": "CCCCCCCCCCCCCCCC(=O)O",
    "palmitic acid": "CCCCCCCCCCCCCCCC(=O)O",
    "ethyl pentadecanoate": "CCCCCCCCCCCCCCC(=O)OCC",
    "16-hentriacontanol": "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCO",
    # Sterols / steroids (representative scaffolds)
    "lathosterol": "CC(C)CCCC(C)C1CCC2C3CC=C4CC(O)CCC4(C)C3CCC12C",
    "7alpha-methylcholesterol": "CC(C)CCCC(C)C1CCC2C3CCC4CC(O)CCC4(C)C3CCC12C",
    "cholesterol scaffold": "CC(C)CCCC(C)C1CCC2C3CCC4CC(O)CCC4(C)C3CCC12C",
    "cholest-7-en-3-ol (3alpha)": "CC(C)CCCC(C)C1CCC2C3CC=C4CC(O)CCC4(C)C3CCC12C",
    "cholest-8(14)-en-3-ol (3alpha,5alpha)": "CC(C)CCCC(C)C1CCC2C3CCC4CC(O)CCC4(C)C3CCC12C",
    # Aromatics / misc
    "2-phenylaniline": "NC1=CC=CC=C1C1=CC=CC=C1",
    "o-aminobiphenyl": "NC1=CC=CC=C1C1=CC=CC=C1",
    # Reference acaricide
    "chlorfenvinphos": "CCOP(=O)(OCC)OC(=CCl)C1=CC=C(Cl)C=C1Cl",
}

# Extra synonyms mapping -> canonical key in CURATED_SMILES
SYNONYMS: dict[str, str] = {
    "alpha caryophyllene": "alpha-caryophyllene",
    "α-caryophyllene": "alpha-caryophyllene",
    "alpha cedrol": "alpha-cedrol",
    "alpha copaene": "alpha-copaene",
    "α-copaene": "alpha-copaene",
    "cis farnesol": "cis-farnesol",
    "2 phenylaniline": "2-phenylaniline",
    "o-aminobiphenyl": "2-phenylaniline",
    "hexadecanoic acid": "hexadecanoic acid",
    "palmitic acid": "hexadecanoic acid",
    "lathosterol": "lathosterol",
    "7alpha methylcholesterol": "7alpha-methylcholesterol",
    "7α-methylcholesterol": "7alpha-methylcholesterol",
    "chlorfenvinphos": "chlorfenvinphos",
}


def normalize_name(name: str) -> str:
    """Lowercase, strip, collapse whitespace, unify dashes/quotes."""
    n = name.strip().lower()
    n = n.replace("–", "-").replace("—", "-")
    n = re.sub(r"\s+", " ", n)
    return n


def get_curated_smiles(name: str) -> tuple[str, str]:
    """
    Return (smiles, source) from the local curated library.
    source is 'curated' or '' if not found.
    """
    if not name or not name.strip():
        return "", ""
    key = normalize_name(name)
    # direct hit
    if key in CURATED_SMILES:
        return CURATED_SMILES[key], "curated"
    # synonym hit
    if key in SYNONYMS:
        canon = SYNONYMS[key]
        return CURATED_SMILES.get(canon, ""), "curated (synonym)"
    # fuzzy: strip parenthetical stereochemistry, e.g. "Cholest-7-en-3-ol (3alpha)"
    no_paren = re.sub(r"\s*\(.*?\)\s*", "", key).strip()
    if no_paren in CURATED_SMILES:
        return CURATED_SMILES[no_paren], "curated (scaffold)"
    return "", ""


def fetch_pubchem_smiles(name: str, timeout: int = 8) -> tuple[str, str]:
    """Live PubChem PUG-REST lookup. Returns (smiles, source) or ('','')."""
    try:
        import requests  # lazy import so offline use still works

        q = urllib.parse.quote(name.strip())
        url = (
            f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{q}"
            "/property/CanonicalSMILES/JSON"
        )
        r = requests.get(url, timeout=timeout)
        if r.status_code != 200:
            return "", ""
        data = r.json()
        props = data.get("PropertyTable", {}).get("Properties", [])
        if props and props[0].get("CanonicalSMILES"):
            return props[0]["CanonicalSMILES"], "pubchem"
        return "", ""
    except Exception:
        return "", ""


def resolve_smiles(name: str, try_pubchem: bool = False) -> tuple[str, str]:
    """
    Resolve a phytochemical name -> (SMILES, source).
    Order: curated -> PubChem (if try_pubchem) -> ('','needs-manual-entry').
    """
    smiles, src = get_curated_smiles(name)
    if smiles:
        return smiles, src
    if try_pubchem and name.strip():
        smiles, src = fetch_pubchem_smiles(name)
        if smiles:
            return smiles, src
    return "", "needs-manual-entry"


def rdkit_descriptors(smiles: str) -> dict:
    """Compute basic RDKit descriptors; returns {} if invalid/empty."""
    try:
        from rdkit import Chem
        from rdkit.Chem import Descriptors, Lipinski

        mol = Chem.MolFromSmiles(smiles) if smiles else None
        if mol is None:
            return {"valid": False}
        return {
            "valid": True,
            "formula": Descriptors.rdMolDescriptors.CalcMolFormula(mol),
            "mw": round(Descriptors.MolWt(mol), 2),
            "logp": round(Descriptors.MolLogP(mol), 2),
            "tpsa": round(Descriptors.TPSA(mol), 2),
            "hbd": Lipinski.NumHDonors(mol),
            "hba": Lipinski.NumHAcceptors(mol),
            "rot_bonds": Lipinski.NumRotatableBonds(mol),
        }
    except Exception:
        return {"valid": False}
