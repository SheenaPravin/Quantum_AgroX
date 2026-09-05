window.AGROX = window.AGROX || {};

window.AGROX.MODEL = {
 "label": "Quantum_AgroX — interactive engineering layer (demo schema)",
 "fingerprint": "Morgan-like topological fingerprint (RDKit.js, 256-bit)",
 "affinity": "Surrogate ridge QSAR fitted live in-browser on 20 reference ligands; per-target calibration via literature anchor affinity. Pre-screening hypothesis only — not docking, not wet-lab.",
 "window": "Guidance mapped to the study's tested topical window 5–80 mg/mL.",
 "caveat": "Treat every number on this app as a hypothesis. Confirm with AutoDock Vina docking and wet-lab assays before any recommendation.",
 "train": [
  {
   "name": "donepezil",
   "smiles": "O=C(c1ccc(cc1)CCN2CCC(CC2)Cc3ccccc3)OC",
   "dG": -9.8,
   "mw": 337.2,
   "clogp": 3.97,
   "tpsa": 29.5,
   "hba": 3,
   "hbd": 0,
   "rotb": 6,
   "rings": 3,
   "arom": 2,
   "heavy": 25,
   "fsp3": 0.41
  },
  {
   "name": "galantamine",
   "smiles": "CN1CCC23C4C5=CC(O)=CC=C5C(=O)C2C1OC3C=C4O",
   "dG": -8.4,
   "mw": 299.1,
   "clogp": 1.79,
   "tpsa": 70,
   "hba": 5,
   "hbd": 2,
   "rotb": 0,
   "rings": 5,
   "arom": 1,
   "heavy": 22,
   "fsp3": 0.47
  },
  {
   "name": "physostigmine",
   "smiles": "CNC(=O)OC1C2CCN(C)C(C1)C3=C2C=CC=C3N(C)C",
   "dG": -8.1,
   "mw": 303.2,
   "clogp": 2.34,
   "tpsa": 44.8,
   "hba": 5,
   "hbd": 1,
   "rotb": 2,
   "rings": 4,
   "arom": 1,
   "heavy": 22,
   "fsp3": 0.59
  },
  {
   "name": "alpha_pinene",
   "smiles": "CC1=CCC2CC1C2(C)C",
   "dG": -6.9,
   "mw": 136.1,
   "clogp": 3,
   "tpsa": 0,
   "hba": 0,
   "hbd": 0,
   "rotb": 0,
   "rings": 3,
   "arom": 0,
   "heavy": 10,
   "fsp3": 0.8
  },
  {
   "name": "beta_pinene",
   "smiles": "CC1(C)C2CCC1C(=C)C2",
   "dG": -6.8,
   "mw": 136.1,
   "clogp": 3,
   "tpsa": 0,
   "hba": 0,
   "hbd": 0,
   "rotb": 0,
   "rings": 2,
   "arom": 0,
   "heavy": 10,
   "fsp3": 0.8
  },
  {
   "name": "linalool",
   "smiles": "CC(=CCCC(C)(O)C=C)C",
   "dG": -6.4,
   "mw": 154.1,
   "clogp": 2.67,
   "tpsa": 20.2,
   "hba": 1,
   "hbd": 1,
   "rotb": 4,
   "rings": 0,
   "arom": 0,
   "heavy": 11,
   "fsp3": 0.6
  },
  {
   "name": "geraniol",
   "smiles": "CC(=CCCC(=CCO)C)C",
   "dG": -6.6,
   "mw": 154.1,
   "clogp": 2.67,
   "tpsa": 20.2,
   "hba": 1,
   "hbd": 1,
   "rotb": 4,
   "rings": 0,
   "arom": 0,
   "heavy": 11,
   "fsp3": 0.6
  },
  {
   "name": "limonene",
   "smiles": "CC1=CCC(CC1)C(=C)C",
   "dG": -6.5,
   "mw": 136.1,
   "clogp": 3.31,
   "tpsa": 0,
   "hba": 0,
   "hbd": 0,
   "rotb": 1,
   "rings": 1,
   "arom": 0,
   "heavy": 10,
   "fsp3": 0.6
  },
  {
   "name": "thymol",
   "smiles": "Cc1ccc(C(C)C)c(O)c1",
   "dG": -7,
   "mw": 150.1,
   "clogp": 2.82,
   "tpsa": 20.2,
   "hba": 1,
   "hbd": 1,
   "rotb": 1,
   "rings": 1,
   "arom": 1,
   "heavy": 11,
   "fsp3": 0.4
  },
  {
   "name": "carvacrol",
   "smiles": "Cc1ccc(C(C)C)cc1O",
   "dG": -7.1,
   "mw": 150.1,
   "clogp": 2.82,
   "tpsa": 20.2,
   "hba": 1,
   "hbd": 1,
   "rotb": 1,
   "rings": 1,
   "arom": 1,
   "heavy": 11,
   "fsp3": 0.4
  },
  {
   "name": "menthol",
   "smiles": "CC1CCC(C(C)C)C(O)C1",
   "dG": -6.6,
   "mw": 156.2,
   "clogp": 2.44,
   "tpsa": 20.2,
   "hba": 1,
   "hbd": 1,
   "rotb": 1,
   "rings": 1,
   "arom": 0,
   "heavy": 11,
   "fsp3": 1
  },
  {
   "name": "camphor",
   "smiles": "CC1(C2CCC1(C(=O)C2)C)C",
   "dG": -6.7,
   "mw": 152.1,
   "clogp": 2.4,
   "tpsa": 17.1,
   "hba": 1,
   "hbd": 0,
   "rotb": 0,
   "rings": 2,
   "arom": 0,
   "heavy": 11,
   "fsp3": 0.9
  },
  {
   "name": "eugenol",
   "smiles": "COc1cc(CC=C)ccc1O",
   "dG": -6.8,
   "mw": 164.1,
   "clogp": 2.13,
   "tpsa": 29.5,
   "hba": 2,
   "hbd": 1,
   "rotb": 3,
   "rings": 1,
   "arom": 1,
   "heavy": 12,
   "fsp3": 0.2
  },
  {
   "name": "caffeine",
   "smiles": "Cn1c(=O)n(C)c2ncn(C)c2n1C",
   "dG": -6.5,
   "mw": 195.1,
   "clogp": 0.27,
   "tpsa": 44.6,
   "hba": 6,
   "hbd": 0,
   "rotb": 0,
   "rings": 2,
   "arom": 1,
   "heavy": 14,
   "fsp3": 0.5
  },
  {
   "name": "nicotine",
   "smiles": "CN1CCC[C@@H]1c2cccnc2",
   "dG": -7.6,
   "mw": 162.1,
   "clogp": 1.85,
   "tpsa": 16.1,
   "hba": 2,
   "hbd": 0,
   "rotb": 1,
   "rings": 2,
   "arom": 1,
   "heavy": 12,
   "fsp3": 0.5
  },
  {
   "name": "chlorfenvinphos",
   "smiles": "CCOP(=O)(O/C=C(/Cl)c1ccc(Cl)cc1)OCCCl",
   "dG": -6.5,
   "mw": 358,
   "clogp": 5.29,
   "tpsa": 44.8,
   "hba": 4,
   "hbd": 0,
   "rotb": 8,
   "rings": 1,
   "arom": 1,
   "heavy": 20,
   "fsp3": 0.33
  },
  {
   "name": "chlorpyrifos",
   "smiles": "CCOP(=S)(OCC)Oc1nc(Cl)nc(Cl)c1Cl",
   "dG": -7.7,
   "mw": 349.9,
   "clogp": 4.11,
   "tpsa": 53.5,
   "hba": 5,
   "hbd": 0,
   "rotb": 6,
   "rings": 1,
   "arom": 1,
   "heavy": 18,
   "fsp3": 0.5
  },
  {
   "name": "propoxur",
   "smiles": "CNC(=O)OC1=CC=CC=C1OC(C)C",
   "dG": -7.2,
   "mw": 209.1,
   "clogp": 2.19,
   "tpsa": 47.6,
   "hba": 4,
   "hbd": 1,
   "rotb": 3,
   "rings": 1,
   "arom": 1,
   "heavy": 15,
   "fsp3": 0.36
  },
  {
   "name": "carbaryl",
   "smiles": "CNC(=O)Oc1c2ccccc2ccc1",
   "dG": -7.3,
   "mw": 201.1,
   "clogp": 2.56,
   "tpsa": 38.3,
   "hba": 3,
   "hbd": 1,
   "rotb": 1,
   "rings": 2,
   "arom": 2,
   "heavy": 15,
   "fsp3": 0.08
  },
  {
   "name": "malathion",
   "smiles": "CCOC(=O)CC(SP(=S)(OCC)OCC)C(=O)OCC",
   "dG": -6.9,
   "mw": 358.1,
   "clogp": 2.9,
   "tpsa": 71.1,
   "hba": 6,
   "hbd": 0,
   "rotb": 11,
   "rings": 0,
   "arom": 0,
   "heavy": 21,
   "fsp3": 0.83
  }
 ],
 "features": [
  "mw",
  "clogp",
  "tpsa",
  "hba",
  "hbd",
  "rotb",
  "rings",
  "arom",
  "heavy",
  "fsp3"
 ],
 "targets": [
  {
   "key": "ache_microp",
   "name": "Acetylcholinesterase",
   "organism": "Rhipicephalus microplus (cattle tick)",
   "cls": "mite",
   "protein": "AChE",
   "accession": "CAA11702",
   "template": "6ARX · chain A",
   "anchor": "Chlorfenvinphos",
   "anchorKey": "chlorfenvinphos",
   "focus": "Tick acaricide target — OP/carbamate gorge",
   "note": "Benchmark target from the supplied paper.",
   "box": "35.4 × 37.3 × 38.8 Å",
   "center": "−62.28, 57.63, 13.27",
   "anchorDG": -6.5,
   "anchorFeat": [
    358,
    5.29,
    44.8,
    4,
    0,
    8,
    1,
    1,
    20,
    0.33
   ]
  },
  {
   "key": "ache_decol",
   "name": "Acetylcholinesterase",
   "organism": "Rhipicephalus decoloratus (cattle tick)",
   "cls": "mite",
   "protein": "AChE",
   "accession": "CAA06980",
   "template": "5YDH · chain A",
   "anchor": "Chlorfenvinphos",
   "anchorKey": "chlorfenvinphos",
   "focus": "Tick acaricide target — OP/carbamate gorge",
   "note": "Benchmark target from the supplied paper.",
   "box": "41.7 × 32.5 × 39.8 Å",
   "center": "56.19, −51.49, 9.45",
   "anchorDG": -6.5,
   "anchorFeat": [
    358,
    5.29,
    44.8,
    4,
    0,
    8,
    1,
    1,
    20,
    0.33
   ]
  },
  {
   "key": "ache_insect",
   "name": "Acetylcholinesterase",
   "organism": "Insect (general)",
   "cls": "insect",
   "protein": "AChE",
   "accession": "—",
   "template": "insect-2 paralog",
   "anchor": "Propoxur",
   "anchorKey": "propoxur",
   "focus": "OP/carbamate neurotoxicity",
   "note": "Generic insect AChE anchor.",
   "anchorDG": -7.2,
   "anchorFeat": [
    209.1,
    2.19,
    47.6,
    4,
    1,
    3,
    1,
    1,
    15,
    0.36
   ]
  },
  {
   "key": "nAChR",
   "name": "Nicotinic acetylcholine receptor",
   "organism": "Insect (sucking pests)",
   "cls": "insect",
   "protein": "nAChR α1β2",
   "accession": "—",
   "template": "α1 subunit",
   "anchor": "Imidacloprid",
   "anchorKey": "imidacloprid",
   "focus": "Neonicotinoid agonist site",
   "note": "Imidacloprid anchor.",
   "anchorDG": -8.9,
   "anchorFeat": [
    256,
    -0.07,
    97.7,
    8,
    1,
    0,
    2,
    1,
    17,
    0.38
   ]
  },
  {
   "key": "VGSC",
   "name": "Voltage-gated sodium channel",
   "organism": "Insect (nerve)",
   "cls": "insect",
   "protein": "VGSC",
   "accession": "kdr site",
   "template": "NaV, pore/II–III",
   "anchor": "Permethrin",
   "anchorKey": "permethrin",
   "focus": "Pyrethroid/DDT site (kdr)",
   "note": "Permethrin anchor.",
   "anchorDG": -8.1,
   "anchorFeat": [
    390.1,
    6.11,
    35.5,
    3,
    0,
    6,
    3,
    2,
    26,
    0.29
   ]
  },
  {
   "key": "GABA",
   "name": "GABA-gated chloride channel",
   "organism": "Insect (nerve)",
   "cls": "insect",
   "protein": "Rdl subunit",
   "accession": "Rdl",
   "template": "A302S hot spot",
   "anchor": null,
   "focus": "Cyclodiene / fipronil site",
   "note": "No anchor bundled — low-confidence calibration.",
   "anchorDG": null,
   "anchorFeat": null
  },
  {
   "key": "ryanodine",
   "name": "Ryanodine receptor",
   "organism": "Insect (muscle)",
   "cls": "insect",
   "protein": "RyR",
   "accession": "—",
   "template": "N-terminal domain",
   "anchor": null,
   "focus": "Diamide insecticide site",
   "note": "No anchor bundled — low-confidence calibration.",
   "anchorDG": null,
   "anchorFeat": null
  },
  {
   "key": "chitin",
   "name": "Chitin synthase",
   "organism": "Insect / mite (cuticle)",
   "cls": "insect",
   "protein": "CHS1",
   "accession": "—",
   "template": "catalytic domain",
   "anchor": "Diflubenzuron",
   "anchorKey": "diflubenzuron",
   "focus": "IGR site (benzoylureas)",
   "note": "Diflubenzuron anchor.",
   "anchorDG": -7.8,
   "anchorFeat": [
    376,
    4.97,
    58.2,
    4,
    2,
    2,
    2,
    2,
    24,
    0.07
   ]
  },
  {
   "key": "GST",
   "name": "Glutathione S-transferase",
   "organism": "Detox enzyme (all pest classes)",
   "cls": "defense",
   "protein": "GST",
   "accession": "—",
   "template": "class Delta",
   "anchor": "Mefenamic acid",
   "anchorKey": "mefenamic_acid",
   "focus": "Detoxification blockers / synergists",
   "note": "Mefenamic acid anchor.",
   "anchorDG": -6.8,
   "anchorFeat": [
    227.1,
    3.44,
    49.3,
    3,
    2,
    3,
    2,
    2,
    17,
    0.07
   ]
  },
  {
   "key": "amylase",
   "name": "Alpha-amylase",
   "organism": "Insect (gut, stored-grain)",
   "cls": "insect",
   "protein": "AMY",
   "accession": "—",
   "template": "catalytic site",
   "anchor": "Limonene",
   "anchorKey": "limonene",
   "focus": "Digestion inhibitor",
   "note": "Limonene (monoterpene) anchor.",
   "anchorDG": -6.2,
   "anchorFeat": [
    136.1,
    3.31,
    0,
    0,
    0,
    1,
    1,
    0,
    10,
    0.6
   ]
  },
  {
   "key": "SDH",
   "name": "Succinate dehydrogenase",
   "organism": "Fungus (respiration, complex II)",
   "cls": "fungus",
   "protein": "SDH",
   "accession": "—",
   "template": "SDHB catalytic",
   "anchor": "Boscalid",
   "anchorKey": "boscalid",
   "focus": "SDHI fungicide site",
   "note": "Boscalid anchor.",
   "anchorDG": -7.4,
   "anchorFeat": [
    342,
    5.31,
    42,
    3,
    1,
    3,
    3,
    3,
    23,
    0
   ]
  },
  {
   "key": "CYP51",
   "name": "CYP51 14α-demethylase",
   "organism": "Fungus (sterol biosynthesis)",
   "cls": "fungus",
   "protein": "CYP51",
   "accession": "—",
   "template": "haem pocket",
   "anchor": "Tebuconazole",
   "anchorKey": "tebuconazole",
   "focus": "DMI / triazole fungicide site",
   "note": "Tebuconazole anchor.",
   "anchorDG": -8,
   "anchorFeat": [
    320.2,
    4.51,
    38,
    3,
    1,
    6,
    2,
    2,
    22,
    0.5
   ]
  },
  {
   "key": "tubulin",
   "name": "β-tubulin (MBC site)",
   "organism": "Fungus (mitosis)",
   "cls": "fungus",
   "protein": "β-tubulin",
   "accession": "—",
   "template": "colchicine-site region",
   "anchor": "Carbendazim",
   "anchorKey": "carbendazim",
   "focus": "Benzimidazole fungicide site",
   "note": "Carbendazim anchor.",
   "anchorDG": -7.6,
   "anchorFeat": [
    191.1,
    1.74,
    67,
    5,
    2,
    1,
    2,
    2,
    14,
    0.11
   ]
  },
  {
   "key": "PSII",
   "name": "Photosystem II subunit D1",
   "organism": "Weed (photosynthesis)",
   "cls": "weed",
   "protein": "D1 (PsbA)",
   "accession": "—",
   "template": "QB pocket",
   "anchor": "Atrazine",
   "anchorKey": "atrazine",
   "focus": "Triazine / urea herbicide site",
   "note": "Atrazine anchor.",
   "anchorDG": -6.1,
   "anchorFeat": [
    215.1,
    1.78,
    62.7,
    5,
    2,
    4,
    1,
    1,
    14,
    0.63
   ]
  },
  {
   "key": "EPSPS",
   "name": "EPSP synthase",
   "organism": "Weed (shikimate pathway)",
   "cls": "weed",
   "protein": "EPSPS",
   "accession": "—",
   "template": "glyphosate pocket",
   "anchor": "Glyphosate",
   "anchorKey": "glyphosate",
   "focus": "Glyphosate site",
   "note": "Glyphosate anchor.",
   "anchorDG": -4.9,
   "anchorFeat": [
    183,
    -0.86,
    98.1,
    6,
    3,
    4,
    0,
    0,
    11,
    0.75
   ]
  },
  {
   "key": "nAChR_nema",
   "name": "Nicotinic AChR",
   "organism": "Nematode (muscle)",
   "cls": "nematode",
   "protein": "AChR",
   "accession": "—",
   "template": "α-subunit",
   "anchor": "Levamisole",
   "anchorKey": "levamisole",
   "focus": "Nematicide (levamisole site)",
   "note": "Levamisole anchor.",
   "anchorDG": -6.7,
   "anchorFeat": [
    170.1,
    0.38,
    27.6,
    3,
    1,
    0,
    3,
    0,
    11,
    0.86
   ]
  },
  {
  "key": "ache2_microp",
  "name": "Acetylcholinesterase 2 (AChE2 variant)",
  "organism": "Rhipicephalus microplus (cattle tick)",
  "cls": "mite",
  "protein": "AChE2",
  "accession": "—",
  "template": "—",
  "anchor": "—",
  "anchorKey": null,
  "focus": "Tick acaricide target — OP/carbamate gorge, BME-2 paralog",
  "note": "Unanchored paralog — generic AChE calibration applies.",
  "box": "—",
  "center": "—",
  "anchorDG": null,
  "anchorFeat": null
 },
  {
  "key": "gaba_mite",
  "name": "GABA-gated chloride channel",
  "organism": "Rhipicephalus microplus (cattle tick)",
  "cls": "mite",
  "protein": "GABA-ClC (Rdl)",
  "accession": "—",
  "template": "—",
  "anchor": "—",
  "anchorKey": null,
  "focus": "Tick acaricide target — fipronil-type site",
  "note": "Unanchored — generic insect calibration applies.",
  "box": "—",
  "center": "—",
  "anchorDG": null,
  "anchorFeat": null
 },
  {
  "key": "vgsc_mite",
  "name": "Voltage-gated sodium channel (kdr site)",
  "organism": "Rhipicephalus microplus (cattle tick)",
  "cls": "mite",
  "protein": "VGSC",
  "accession": "—",
  "template": "—",
  "anchor": "—",
  "anchorKey": null,
  "focus": "Tick acaricide target — pyrethroid site I",
  "note": "Unanchored — generic insect VGSC calibration applies.",
  "box": "—",
  "center": "—",
  "anchorDG": null,
  "anchorFeat": null
 },
  {
  "key": "glucl_mite",
  "name": "Glutamate-gated chloride channel",
  "organism": "Rhipicephalus microplus (cattle tick)",
  "cls": "mite",
  "protein": "GluCl",
  "accession": "—",
  "template": "—",
  "anchor": "—",
  "anchorKey": null,
  "focus": "Tick acaricide target — abamectin/milbemycin site",
  "note": "Unanchored — generic calibration applies.",
  "box": "—",
  "center": "—",
  "anchorDG": null,
  "anchorFeat": null
 },
  {
  "key": "glucl_insect",
  "name": "Glutamate-gated chloride channel",
  "organism": "Musca domestica (house fly)",
  "cls": "insect",
  "protein": "GluCl",
  "accession": "—",
  "template": "—",
  "anchor": "—",
  "anchorKey": null,
  "focus": "Insecticide target — macrocyclic lactones",
  "note": "Unanchored — generic calibration applies.",
  "box": "—",
  "center": "—",
  "anchorDG": null,
  "anchorFeat": null
 },
  {
  "key": "octopamine",
  "name": "Octopamine receptor",
  "organism": "Rhipicephalus microplus (cattle tick)",
  "cls": "insect",
  "protein": "OctR",
  "accession": "—",
  "template": "—",
  "anchor": "—",
  "anchorKey": null,
  "focus": "Formamidine (amitraz) site — acaricide",
  "note": "Unanchored — generic calibration applies.",
  "box": "—",
  "center": "—",
  "anchorDG": null,
  "anchorFeat": null
 },
  {
  "key": "tyramine",
  "name": "Tyramine receptor",
  "organism": "Apis mellifera (honeybee)",
  "cls": "insect",
  "protein": "TyrR",
  "accession": "—",
  "template": "—",
  "anchor": "—",
  "anchorKey": null,
  "focus": "Biogenic-amine GPCR — selectivity comparator",
  "note": "Unanchored — generic calibration applies.",
  "box": "—",
  "center": "—",
  "anchorDG": null,
  "anchorFeat": null
 },
  {
  "key": "jh_receptor",
  "name": "Juvenile hormone receptor (Met)",
  "organism": "Tribolium castaneum (red flour beetle)",
  "cls": "insect",
  "protein": "Met",
  "accession": "—",
  "template": "—",
  "anchor": "—",
  "anchorKey": null,
  "focus": "IGR target — methoprene-type mimics",
  "note": "Unanchored — generic calibration applies.",
  "box": "—",
  "center": "—",
  "anchorDG": null,
  "anchorFeat": null
 },
  {
  "key": "ecdysone",
  "name": "Ecdysone receptor (EcR/USP)",
  "organism": "Helicoverpa armigera (cotton bollworm)",
  "cls": "insect",
  "protein": "EcR",
  "accession": "—",
  "template": "—",
  "anchor": "—",
  "anchorKey": null,
  "focus": "IGR target — diacylhydrazine moulting disruptors",
  "note": "Unanchored — generic calibration applies.",
  "box": "—",
  "center": "—",
  "anchorDG": null,
  "anchorFeat": null
 },
  {
  "key": "mt_ci",
  "name": "Mitochondrial complex I (NADH dehydrogenase)",
  "organism": "Blattella germanica (German cockroach)",
  "cls": "insect",
  "protein": "ND1",
  "accession": "—",
  "template": "—",
  "anchor": "—",
  "anchorKey": null,
  "focus": "Respiratory inhibitor target — rotenone site",
  "note": "Unanchored — generic calibration applies.",
  "box": "—",
  "center": "—",
  "anchorDG": null,
  "anchorFeat": null
 },
  {
  "key": "mt_atp",
  "name": "Mitochondrial ATP synthase",
  "organism": "Aedes aegypti (yellow fever mosquito)",
  "cls": "insect",
  "protein": "ATPase",
  "accession": "—",
  "template": "—",
  "anchor": "—",
  "anchorKey": null,
  "focus": "Respiratory inhibitor target",
  "note": "Unanchored — generic calibration applies.",
  "box": "—",
  "center": "—",
  "anchorDG": null,
  "anchorFeat": null
 },
  {
  "key": "acc_insect",
  "name": "Acetyl-CoA carboxylase (lipid)",
  "organism": "Bemisia tabaci (whitefly)",
  "cls": "insect",
  "protein": "ACC",
  "accession": "—",
  "template": "—",
  "anchor": "—",
  "anchorKey": null,
  "focus": "Lipid-synthesis inhibitor (spirotetramat-like)",
  "note": "Unanchored — generic calibration applies.",
  "box": "—",
  "center": "—",
  "anchorDG": null,
  "anchorFeat": null
 },
  {
  "key": "machr_insect",
  "name": "Muscarinic acetylcholine receptor",
  "organism": "Helicoverpa zea (corn earworm)",
  "cls": "insect",
  "protein": "mAChR",
  "accession": "—",
  "template": "—",
  "anchor": "—",
  "anchorKey": null,
  "focus": "Insect neuroreceptor target",
  "note": "Unanchored — generic calibration applies.",
  "box": "—",
  "center": "—",
  "anchorDG": null,
  "anchorFeat": null
 },
  {
  "key": "ache_bee",
  "name": "Acetylcholinesterase (synaptic)",
  "organism": "Apis mellifera (honeybee) — selectivity comparator",
  "cls": "insect",
  "protein": "AChE",
  "accession": "—",
  "template": "—",
  "anchor": "—",
  "anchorKey": null,
  "focus": "Pollinator AAChE selectivity comparator",
  "note": "Unanchored — used for selectivity balance, not efficacy.",
  "box": "—",
  "center": "—",
  "anchorDG": null,
  "anchorFeat": null
 },
  {
  "key": "nachr_bee",
  "name": "Nicotinic acetylcholine receptor (α-subunit)",
  "organism": "Apis mellifera (honeybee) — selectivity comparator",
  "cls": "insect",
  "protein": "nAChR",
  "accession": "—",
  "template": "—",
  "anchor": "—",
  "anchorKey": null,
  "focus": "Pollinator nAChR selectivity comparator (neonicotinoid site)",
  "note": "Unanchored — used for selectivity balance.",
  "box": "—",
  "center": "—",
  "anchorDG": null,
  "anchorFeat": null
 },
  {
  "key": "carboxyles",
  "name": "Carboxylesterase (detox)",
  "organism": "Rhipicephalus microplus (cattle tick)",
  "cls": "defense",
  "protein": "CaE",
  "accession": "—",
  "template": "—",
  "anchor": "—",
  "anchorKey": null,
  "focus": "Detoxification enzyme — synergist target",
  "note": "Unanchored — generic calibration applies.",
  "box": "—",
  "center": "—",
  "anchorDG": null,
  "anchorFeat": null
 },
  {
  "key": "catalase",
  "name": "Catalase",
  "organism": "Rhipicephalus microplus (cattle tick)",
  "cls": "defense",
  "protein": "CAT",
  "accession": "—",
  "template": "—",
  "anchor": "—",
  "anchorKey": null,
  "focus": "Antioxidant defense enzyme",
  "note": "Unanchored — generic calibration applies.",
  "box": "—",
  "center": "—",
  "anchorDG": null,
  "anchorFeat": null
 },
  {
  "key": "sod",
  "name": "Superoxide dismutase",
  "organism": "Rhipicephalus microplus (cattle tick)",
  "cls": "defense",
  "protein": "SOD",
  "accession": "—",
  "template": "—",
  "anchor": "—",
  "anchorKey": null,
  "focus": "Antioxidant defense enzyme",
  "note": "Unanchored — generic calibration applies.",
  "box": "—",
  "center": "—",
  "anchorDG": null,
  "anchorFeat": null
 },
  {
  "key": "cyp450",
  "name": "Cytochromes P450 monooxygenase",
  "organism": "Rhipicephalus microplus (cattle tick)",
  "cls": "defense",
  "protein": "CYP",
  "accession": "—",
  "template": "—",
  "anchor": "—",
  "anchorKey": null,
  "focus": "Metabolic detoxification — synergist target",
  "note": "Unanchored — generic calibration applies.",
  "box": "—",
  "center": "—",
  "anchorDG": null,
  "anchorFeat": null
 },
  {
  "key": "ugt",
  "name": "UDP-glycosyltransferase",
  "organism": "Rhipicephalus microplus (cattle tick)",
  "cls": "defense",
  "protein": "UGT",
  "accession": "—",
  "template": "—",
  "anchor": "—",
  "anchorKey": null,
  "focus": "Phase II detoxification — synergist target",
  "note": "Unanchored — generic calibration applies.",
  "box": "—",
  "center": "—",
  "anchorDG": null,
  "anchorFeat": null
 },
  {
  "key": "bc1_fungi",
  "name": "Cytochrome bc1 complex (QoI site)",
  "organism": "Botrytis cinerea (grey mould)",
  "cls": "fungus",
  "protein": "bc1",
  "accession": "—",
  "template": "—",
  "anchor": "—",
  "anchorKey": null,
  "focus": "Fungicide target — strobilurin QoI site",
  "note": "Unanchored — generic calibration applies.",
  "box": "—",
  "center": "—",
  "anchorDG": null,
  "anchorFeat": null
 },
  {
  "key": "glcn6p_fungi",
  "name": "Glucosamine-6-phosphate synthase",
  "organism": "Fusarium graminearum (wheat blight)",
  "cls": "fungus",
  "protein": "GlmS",
  "accession": "—",
  "template": "—",
  "anchor": "—",
  "anchorKey": null,
  "focus": "Fungicide target — chitin/ergosterol precursor pathway",
  "note": "Unanchored — generic calibration applies.",
  "box": "—",
  "center": "—",
  "anchorDG": null,
  "anchorFeat": null
 },
  {
  "key": "glucanase_fungi",
  "name": "β-1,3-Glucan synthase",
  "organism": "Magnaporthe oryzae (rice blast)",
  "cls": "fungus",
  "protein": "FKS",
  "accession": "—",
  "template": "—",
  "anchor": "—",
  "anchorKey": null,
  "focus": "Fungicide target — cell-wall glucan synthesis",
  "note": "Unanchored — generic calibration applies.",
  "box": "—",
  "center": "—",
  "anchorDG": null,
  "anchorFeat": null
 },
  {
  "key": "chitin_fungi",
  "name": "Chitin synthase",
  "organism": "Botrytis cinerea (grey mould)",
  "cls": "fungus",
  "protein": "Chs",
  "accession": "—",
  "template": "—",
  "anchor": "—",
  "anchorKey": null,
  "focus": "Fungicide target — cell-wall chitin synthesis",
  "note": "Unanchored — generic calibration applies.",
  "box": "—",
  "center": "—",
  "anchorDG": null,
  "anchorFeat": null
 },
  {
  "key": "als_weed",
  "name": "Acetolactate synthase (ALS/AHAS)",
  "organism": "Amaranthus palmeri (Palmer amaranth)",
  "cls": "weed",
  "protein": "ALS",
  "accession": "—",
  "template": "—",
  "anchor": "—",
  "anchorKey": null,
  "focus": "Herbicide target — branched-chain amino-acid biosynthesis",
  "note": "Unanchored — generic calibration applies.",
  "box": "—",
  "center": "—",
  "anchorDG": null,
  "anchorFeat": null
 },
  {
  "key": "psi_weed",
  "name": "Photosystem I acceptors",
  "organism": "Chenopodium album (fat-hen)",
  "cls": "weed",
  "protein": "PSI",
  "accession": "—",
  "template": "—",
  "anchor": "—",
  "anchorKey": null,
  "focus": "Herbicide target — bipyridylium (paraquat) site",
  "note": "Unanchored — generic calibration applies.",
  "box": "—",
  "center": "—",
  "anchorDG": null,
  "anchorFeat": null
 },
  {
  "key": "hppd_weed",
  "name": "4-Hydroxyphenylpyruvate dioxygenase",
  "organism": "Avena fatua (wild oat)",
  "cls": "weed",
  "protein": "HPPD",
  "accession": "—",
  "template": "—",
  "anchor": "—",
  "anchorKey": null,
  "focus": "Herbicide target — carotenoid/plastoquinone biosynthesis",
  "note": "Unanchored — generic calibration applies.",
  "box": "—",
  "center": "—",
  "anchorDG": null,
  "anchorFeat": null
 },
  {
  "key": "ppo_weed",
  "name": "Protoporphyrinogen oxidase (PPO)",
  "organism": "Ipomoea hederacea (ivyleaf morningglory)",
  "cls": "weed",
  "protein": "PPO",
  "accession": "—",
  "template": "—",
  "anchor": "—",
  "anchorKey": null,
  "focus": "Herbicide target — chlorophyll biosynthesis",
  "note": "Unanchored — generic calibration applies.",
  "box": "—",
  "center": "—",
  "anchorDG": null,
  "anchorFeat": null
 },
  {
  "key": "pds_weed",
  "name": "Phytoene desaturase",
  "organism": "Avena fatua (wild oat)",
  "cls": "weed",
  "protein": "PDS",
  "accession": "—",
  "template": "—",
  "anchor": "—",
  "anchorKey": null,
  "focus": "Herbicide target — carotenoid bleaching",
  "note": "Unanchored — generic calibration applies.",
  "box": "—",
  "center": "—",
  "anchorDG": null,
  "anchorFeat": null
 },
  {
  "key": "accase_weed",
  "name": "Acetyl-CoA carboxylase (ACCase, plastid)",
  "organism": "Lolium rigidum (annual ryegrass)",
  "cls": "weed",
  "protein": "ACCase",
  "accession": "—",
  "template": "—",
  "anchor": "—",
  "anchorKey": null,
  "focus": "Herbicide target — graminicide lipid synthesis",
  "note": "Unanchored — generic calibration applies.",
  "box": "—",
  "center": "—",
  "anchorDG": null,
  "anchorFeat": null
 },
  {
  "key": "gs_weed",
  "name": "Glutamine synthetase",
  "organism": "Sorghum halepense (johnsongrass)",
  "cls": "weed",
  "protein": "GS",
  "accession": "—",
  "template": "—",
  "anchor": "—",
  "anchorKey": null,
  "focus": "Herbicide target — amino-acid assimilation (glufosinate-like)",
  "note": "Unanchored — generic calibration applies.",
  "box": "—",
  "center": "—",
  "anchorDG": null,
  "anchorFeat": null
 },
  {
  "key": "slo1_nema",
  "name": "SLO-1 Ca²⁺-activated K⁺ channel",
  "organism": "Haemonchus contortus (barber's pole worm)",
  "cls": "nematode",
  "protein": "SLO-1",
  "accession": "—",
  "template": "—",
  "anchor": "—",
  "anchorKey": null,
  "focus": "Anthelmintic target — emodepside site",
  "note": "Unanchored — generic calibration applies.",
  "box": "—",
  "center": "—",
  "anchorDG": null,
  "anchorFeat": null
 },
  {
  "key": "glucl_nema",
  "name": "Glutamate-gated chloride channel",
  "organism": "Haemonchus contortus (barber's pole worm)",
  "cls": "nematode",
  "protein": "GluCl",
  "accession": "—",
  "template": "—",
  "anchor": "—",
  "anchorKey": null,
  "focus": "Anthelmintic target — avermectin site",
  "note": "Unanchored — generic calibration applies.",
  "box": "—",
  "center": "—",
  "anchorDG": null,
  "anchorFeat": null
 },
  {
  "key": "tubulin_nema",
  "name": "β-Tubulin (benzimidazole site)",
  "organism": "Haemonchus contortus (barber's pole worm)",
  "cls": "nematode",
  "protein": "TUB",
  "accession": "—",
  "template": "—",
  "anchor": "—",
  "anchorKey": null,
  "focus": "Anthelmintic/fungicide target — microtubule disruptors",
  "note": "Unanchored — generic calibration applies.",
  "box": "—",
  "center": "—",
  "anchorDG": null,
  "anchorFeat": null
 }
  ],
  "smilesExamples": [
   {
   "name": "α-Pinene (monoterpene)",
   "smi": "CC1=CCC2CC1C2(C)C",
   "cls": "botanical"
  },
   {
   "name": "β-Pinene (monoterpene)",
   "smi": "CC1(C)C2CCC1C(=C)C2",
   "cls": "botanical"
  },
   {
   "name": "Limonene (monoterpene)",
   "smi": "CC1=CCC(CC1)=C(C)C",
   "cls": "botanical"
  },
   {
   "name": "Linalool (monoterpene)",
   "smi": "CC(C)=CCCC(C)(O)C=C",
   "cls": "botanical"
  },
   {
   "name": "Geraniol (monoterpene)",
   "smi": "CC(C)=CCCC(C)=CCO",
   "cls": "botanical"
  },
   {
   "name": "Thymol (monoterpene phenol)",
   "smi": "CC(C)c1cc(C)ccc1O",
   "cls": "botanical"
  },
   {
   "name": "Carvacrol (monoterpene phenol)",
   "smi": "CC(C)c1ccc(O)cc1C",
   "cls": "botanical"
  },
   {
   "name": "Menthol (monoterpene)",
   "smi": "CC(C)C1CCC(C)CC1O",
   "cls": "botanical"
  },
   {
   "name": "Eugenol (phenylpropanoid)",
   "smi": "COc1cc(CC=C)ccc1O",
   "cls": "botanical"
  },
   {
   "name": "Anethole (phenylpropanoid)",
   "smi": "COc1ccc(C=CC)cc1",
   "cls": "botanical"
  },
   {
   "name": "Cinnamaldehyde (phenylpropanoid)",
   "smi": "O=C/C=C/c1ccccc1",
   "cls": "botanical"
  },
   {
   "name": "Camphor (bicyclic ketone)",
   "smi": "CC1(C)C2CCC1(C)C(=O)C2",
   "cls": "botanical"
  },
   {
   "name": "Phytol (diterpenoid)",
   "smi": "CC(C)CCCC(C)CCCC(C)CCCC(C)=CCO",
   "cls": "botanical"
  },
   {
   "name": "Squalene (triterpene)",
   "smi": "CC(=CCCC(=CCCC(=CCCC=C(C)C)C)C)C",
   "cls": "botanical"
  },
   {
   "name": "Quercetin (flavonoid)",
   "smi": "O=C1c2c(O)cc(O)cc2OC(c2cc(O)c(O)c(O)c2)=C1O",
   "cls": "botanical"
  },
   {
   "name": "Kaempferol (flavonoid)",
   "smi": "O=C1c2c(O)cc(O)cc2OC(c2ccc(O)cc2)=C1O",
   "cls": "botanical"
  },
   {
   "name": "Nicotine (alkaloid)",
   "smi": "CN1CCC[C@H]1c1cccnc1",
   "cls": "botanical"
  },
   {
   "name": "Donepezil (AChE drug)",
   "smi": "COc1ccc2CC(CC(=O)N3CCN(CC3)Cc3ccccc3)Cc2c1OC",
   "cls": "drug"
  },
   {
   "name": "Galantamine (AChE drug)",
   "smi": "COc1cc2c(cc1OC)O[C@H]1C[C@@H](O)[C@@H](CN2CC=C)C1",
   "cls": "drug"
  },
   {
   "name": "Physostigmine (AChE drug)",
   "smi": "CNC(=O)OC1=CC2=C(C=C1)N(C)C1CC(N(C)C)C2CC1",
   "cls": "drug"
  },
   {
   "name": "Chlorpyrifos (OP insecticide)",
   "smi": "CCOP(=S)(OCC)Oc1nc(Cl)nc(Cl)c1Cl",
   "cls": "insecticide"
  },
   {
   "name": "Chlorfenvinphos (OP acaricide)",
   "smi": "Cl/C(Cl)=C(OP(=O)(OCC)OCC)c1ccc(Cl)cc1Cl",
   "cls": "insecticide"
  },
   {
   "name": "Propoxur (carbamate)",
   "smi": "CC(C)OC(=O)NC1=CC=CC=C1OC",
   "cls": "insecticide"
  },
   {
   "name": "Carbaryl (carbamate)",
   "smi": "CNC(=O)Oc1cccc2ccccc12",
   "cls": "insecticide"
  },
   {
   "name": "Malathion (OP)",
   "smi": "CCOC(=O)CC(SP(=S)(OC)OC)C(=O)OCC",
   "cls": "insecticide"
  },
   {
   "name": "Permethrin (pyrethroid)",
   "smi": "CC1(C)C(C=C(Cl)Cl)C1C(=O)OCCc1cccc(OCc2ccccc2)c1",
   "cls": "insecticide"
  },
   {
   "name": "Imidacloprid (neonicotinoid)",
   "smi": "Clc1ccc(CN2CCNC2=N[N+](=O)[O-])cn1",
   "cls": "insecticide"
  },
   {
   "name": "Acetamiprid (neonicotinoid)",
   "smi": "CN=C(C#N)NCC1=CC=C(N=C1)Cl",
   "cls": "insecticide"
  },
   {
   "name": "Amitraz (formamidine acaricide)",
   "smi": "CN(C=Nc1ccc(C)cc1C)C=Nc1ccc(C)cc1C",
   "cls": "insecticide"
  },
   {
   "name": "Diflubenzuron (IGR)",
   "smi": "Fc1ccc(F)cc1NC(=O)NC(=O)c1ccccc1Cl",
   "cls": "insecticide"
  },
   {
   "name": "Methomyl (carbamate)",
   "smi": "CNC(=O)O/N=C(\\C)SC",
   "cls": "insecticide"
  },
   {
   "name": "Cyromazine (triazine IGR)",
   "smi": "Nc1nc(N)nc(NC2CC2)n1",
   "cls": "insecticide"
  },
   {
   "name": "Levamisole (anthelmintic)",
   "smi": "C1CS[C@H]2CN=C(N2C1)c3ccccc3",
   "cls": "nematicide"
  },
   {
   "name": "Albendazole (anthelmintic)",
   "smi": "CCCSc1ccc2nc(NC(=O)OC)[nH]c2c1",
   "cls": "nematicide"
  },
   {
   "name": "Fenbendazole (anthelmintic)",
   "smi": "COC(=O)Nc1nc2ccc(Sc3ccccc3)cc2[nH]1",
   "cls": "nematicide"
  },
   {
   "name": "Atrazine (triazine herbicide)",
   "smi": "CCNc1nc(Cl)nc(NC(C)C)n1",
   "cls": "herbicide"
  },
   {
   "name": "Paraquat (bipyridylium herbicide)",
   "smi": "C[n+]1ccc(-c2cc[n+](C)cc2)cc1",
   "cls": "herbicide"
  },
   {
   "name": "Diuron (urea herbicide)",
   "smi": "CC(C)N(C)C(=O)Nc1ccc(Cl)cc1Cl",
   "cls": "herbicide"
  },
   {
   "name": "2,4-D (auxin herbicide)",
   "smi": "OC(=O)COc1ccc(Cl)cc1Cl",
   "cls": "herbicide"
  },
   {
   "name": "Dicamba (benzoic herbicide)",
   "smi": "COc1c(Cl)cc(Cl)cc1C(=O)O",
   "cls": "herbicide"
  },
   {
   "name": "Metolachlor (chloroacetanilide)",
   "smi": "CCC1=CC=CC=C1N(C(C)COC)C(=O)CCl",
   "cls": "herbicide"
  },
   {
   "name": "Pretilachlor (chloroacetanilide)",
   "smi": "CCC1=CC=CC=C1N(CCCOCC)C(=O)CCl",
   "cls": "herbicide"
  },
   {
   "name": "Pendimethalin (dinitroaniline)",
   "smi": "CCC(CC)Nc1c(C)c(C)c([N+](=O)[O-])cc1[N+](=O)[O-]",
   "cls": "herbicide"
  },
   {
   "name": "Glyphosate (amino-acid herbicide)",
   "smi": "OC(=O)CNCP(=O)(O)O",
   "cls": "herbicide"
  },
   {
   "name": "Bromoxynil (nitrile herbicide)",
   "smi": "N#Cc1cc(Br)c(O)c(Br)c1C#N",
   "cls": "herbicide"
  },
   {
   "name": "Carbendazim (benzimidazole)",
   "smi": "COC(=O)Nc1nc2ccccc2[nH]1",
   "cls": "fungicide"
  },
   {
   "name": "Thiabendazole (benzimidazole)",
   "smi": "C1=CSC(=N1)c1nc2ccccc2[nH]1",
   "cls": "fungicide"
  },
   {
   "name": "Chlorothalonil (nitrile)",
   "smi": "N#Cc1c(Cl)c(Cl)c(C#N)c(Cl)c1Cl",
   "cls": "fungicide"
  },
   {
   "name": "Tebuconazole (triazole)",
   "smi": "CC(C)(C)CC(CC1=CNN=N1)c1ccc(Cl)cc1",
   "cls": "fungicide"
  },
   {
   "name": "Boscalid (pyridine carboxamide)",
   "smi": "O=C(Nc1ccccc1-c1ccc(Cl)cc1)c1cccnc1Cl",
   "cls": "fungicide"
  }
  ],
  "mutations": [
  {
   "name": "AChE G119S",
   "target": [
    "ache_microp",
    "ache_decol",
    "ache_insect"
   ],
   "gene": "AChE",
   "delta": 0.6,
   "effect": "Reduces organophosphate/carbamate sensitivity (resistance)",
   "note": "Classic OP-resistance substitution; analogous form in tick BME-2 AChE suspected.",
   "lit": "Weill et al., EMBO J 2003"
  },
  {
   "name": "AChE F331W / W279R",
   "target": [
    "ache_microp",
    "ache_decol"
   ],
   "gene": "AChE",
   "delta": 0.5,
   "effect": "Broad OP/carbamate resistance under acaricide selection",
   "note": "Reported in Rhipicephalus (Boophilus) AChE gene screens.",
   "lit": "Chen et al., acaricide-resistance screens"
  },
  {
   "name": "VGSC kdr T929I + L1014F",
   "target": [
    "VGSC"
   ],
   "gene": "VGSC",
   "delta": 0.5,
   "effect": "Pyrethroid knock-down resistance (kdr)",
   "note": "The 'diagonal' mutation pair reduces pyrethroid binding.",
   "lit": "Williamson et al. 1996"
  },
  {
   "name": "GABA-Rdl A302S",
   "target": [
    "GABA"
   ],
   "gene": "Rdl",
   "delta": 0.6,
   "effect": "Cyclodiene / fipronil resistance",
   "note": "Canonical Rdl mutation across arthropods.",
   "lit": "ffrench-Constant 1993 (Drosophila Rdl)"
  },
  {
   "name": "nAChR Y151S",
   "target": [
    "nAChR"
   ],
   "gene": "nAChR α",
   "delta": 0.4,
   "effect": "Imidacloprid binding loss",
   "note": "Neonicotinoid resistance in aphids.",
   "lit": "Liu & Han 2006"
  },
  {
   "name": "nAChR R81T",
   "target": [
    "nAChR"
   ],
   "gene": "nAChR β1",
   "delta": 0.5,
   "effect": "Neonicotinoid resistance (sucking pests)",
   "note": "Field-selected R81T in aphids.",
   "lit": "Bass et al. 2011"
  },
  {
   "name": "β-tubulin E198K",
   "target": [
    "tubulin"
   ],
   "gene": "β-tubulin",
   "delta": 0.3,
   "effect": "Benzimidazole (MBC) fungicide resistance",
   "note": "Classic MBC-resistance substitution.",
   "lit": "Ko et al. 1986 (A. nidulans)"
  },
  {
   "name": "CYP P450 over-expression",
   "target": null,
   "gene": "CYP",
   "delta": null,
   "effect": "Metabolic resistance (not structural)",
   "note": "Over-expression of P450 / esterase metabolisers — needs expression data, not docking.",
   "lit": "Phase-6 roadmap"
  }
 ],
 "modules": [
  {
   "id": "datahub",
   "name": "AgroDataHub",
   "func": "Import & standardize",
   "desc": "Ingest and standardize GC-MS, SDF/SMILES, FASTA, VCF, bioassay and docking datasets with an auditable provenance step.",
   "icon": "M13 3h7v7M16 8l4-4M6 3H3v17h17v-3M3 11h6v10M14 11h3M14 16h3M14 21h3"
  },
  {
   "id": "phytox",
   "name": "AgroPhytoX",
   "func": "Phytochemical library",
   "desc": "Natural-product and phytochemical library; plant part, extraction solvent, abundance and structure metadata.",
   "icon": "M12 3c1 3 3 4 3 7a3 3 0 0 1-6 0c0-3 2-4 3-7ZM12 13v8M8 21l4-4M16 21l-4-4"
  },
  {
   "id": "syntheticx",
   "name": "AgroSyntheticX",
   "func": "Draw → SMILES",
   "desc": "Sketch a synthetic candidate on the in-dashboard editor; export a canonical, RDKit-validated SMILES into the Molecule Lab and AgroDockX pair.",
   "icon": "M12 19l7-7a4.95 4.95 0 0 0-7-7l-7 7v7h7ZM14 6l4 4M16 3h5v5"
  },
  {
   "id": "targetx",
   "name": "AgroTargetX",
   "func": "Target discovery",
   "desc": "Target discovery for pests, pathogens, fungi, bacteria, nematodes and weeds.",
   "icon": "M12 3l2.5 5.5 6 .5-4.5 4 1.2 5.9L12 16l-5.2 2.9L8 13l-4.5-4 6-.5L12 3Z"
  },
  {
   "id": "sitemap",
   "name": "AgroSiteMap",
   "func": "Binding-site analysis",
   "desc": "Binding-site, pocket, cryptic-pocket and druggability analysis.",
   "icon": "M4 5h16M4 5l3-2M4 5l3 2M20 5l-3-2M20 5l-3 2M12 5v14M12 19l-3-2M12 19l3-2"
  },
  {
   "id": "dockx",
   "name": "AgroDockX",
   "func": "Docking & affinity ranking",
   "desc": "Protein–ligand docking and affinity ranking.",
   "icon": "M5 4h14v6H5V4ZM5 14h14v6H5v-6Z"
  },
  {
   "id": "dosex",
   "name": "AgroDoseX",
   "func": "Dose–response modelling",
   "desc": "Dose-response modelling, LC50/LC90 and time-dependent efficacy prediction.",
   "icon": "M3 12h4l2-5 2 10 2-6 2 8h10"
  },
  {
   "id": "qml",
   "name": "AgroQML",
   "func": "Quantum learning",
   "desc": "Quantum kernels, QSVM, VQC, QNN and hybrid quantum-classical prediction.",
   "icon": "M12 2l3 5 5 1-3 5 3 5-5 1-3 5-3-5-5-1 3-5-3-5 5-1 3-5Z"
  },
  {
   "id": "md",
   "name": "AgroMD",
   "func": "Molecular dynamics",
   "desc": "MD, RMSD/RMSF, interaction persistence and MM/PBSA thermodynamic estimates.",
   "icon": "M3 12c3-5 6-5 9 0s6 5 9 0M3 17c3-5 6-5 9 0s6 5 9 0"
  },
  {
   "id": "resistance",
   "name": "AgroResistanceScan",
   "func": "Resistance hypotheses",
   "desc": "Mutation → structure → affinity-change workflow for resistance hypotheses.",
   "icon": "M12 2v20M12 2l5 4M12 6L7 2M12 22l-5-4M12 18l5 4"
  },
  {
   "id": "select",
   "name": "AgroSelect",
   "func": "Selectivity scoring",
   "desc": "Pest/pathogen versus crop/beneficial-organism selectivity scoring.",
   "icon": "M12 2l6 3v6c0 4.5-2.6 8-6 11-3.4-3-6-6.5-6-11V5l6-3ZM9 11l2 2 4-4"
  },
  {
   "id": "ecorisk",
   "name": "AgroEcoRisk",
   "func": "Environmental risk",
   "desc": "Environmental and non-target risk prioritization.",
   "icon": "M12 3l9 16H3l9-16ZM12 10v4M12 17.2v.4"
  },
  {
   "id": "synergyx",
   "name": "AgroSynergyX",
   "func": "Combination analysis",
   "desc": "Combination/synergy analysis for botanical mixtures and pesticide combinations.",
   "icon": "M4 7h16M4 12h16M4 17h16M7 7V3M7 21v-4M17 7V3M17 21v-4M7 12v5M17 7v5"
  },
  {
   "id": "optimize",
   "name": "AgroOptimize",
   "func": "Multi-objective optimization",
   "desc": "Multi-objective candidate and dose/formulation optimization.",
   "icon": "M5 20V10M10 20V4M15 20v-8M20 20V7"
  },
  {
   "id": "report",
   "name": "AgroReport",
   "func": "Reproducible reporting",
   "desc": "Reproducible reports, rankings, figures and audit trails.",
   "icon": "M6 2h9l4 4v16H6V2ZM15 2v5h4M9 11h6M9 15h6M9 7h2"
  }
 ],
 "refparams": {
  "R. decoloratus AChE": {
   "accession": "CAA06980",
   "template": "5YDH, chain A",
   "center": "56.1947, −51.4858, 9.4480",
   "box": "41.65 × 32.54 × 39.77 Å",
   "exhaust": "8",
   "modes": "9"
  },
  "R. microplus AChE": {
   "accession": "CAA11702",
   "template": "6ARX, chain A",
   "center": "−62.2842, 57.6277, 13.2699",
   "box": "35.36 × 37.33 × 38.82 Å",
   "exhaust": "8",
   "modes": "9"
  }
 },
 "phases": [
  "Phase 1 — Benchmark: run the supplied paper-derived dataset through classical and QML models.",
  "Phase 2 — Chemical intelligence: canonical structures, RDKit descriptors and full compound metadata.",
  "Phase 3 — Target intelligence: automated target retrieval, binding-site prediction and docking.",
  "Phase 4 — Experimental learning: raw replicate-level bioassays and measured AChE inhibition.",
  "Phase 5 — Advanced QML: QNN/QGNN, quantum metric learning and quantum Bayesian optimization.",
  "Phase 6 — Safety/resistance: non-target prediction, ecological endpoints and mutation-resistance modelling.",
  "Phase 7 — AgroX decision engine: auditable multi-objective candidate score."
 ],
 "selectivity": [
  {
   "family": "Acetylcholinesterase blockers (OP/carbamate-type)",
   "chemistry": "ester + thio/aryl, HBD≈0, logP 2–5",
   "beneficial": "HIGH caution — pollinator & parasitoid acute toxicity",
   "flag": "HIGH",
   "note": "Neuroactives hit bees at low doses."
  },
  {
   "family": "Nicotinic actives (neonicotinoid-type)",
   "chemistry": "heteroaryl-N, HBA high, low logP",
   "beneficial": "HIGH caution — systemic bee exposure",
   "flag": "HIGH",
   "note": "Systemic translocation increases blossom exposure."
  },
  {
   "family": "Terpenoid botanicals (mono/sesquiterpenes)",
   "chemistry": "lipophilic volatiles, logP 2.5–3.5",
   "beneficial": "LOW–MOD — repellent, short residual",
   "flag": "MOD",
   "note": "Fast photodegradation limits non-target burden."
  },
  {
   "family": "Phenol / coumarin botanicals",
   "chemistry": "phenolic OH, logP 2–3.5",
   "beneficial": "LOW — feeding deterrent, biodegradable",
   "flag": "LOW",
   "note": "Time-limited, rain-sensitive."
  },
  {
   "family": "Lipophilic IGRs (benzoylurea-type)",
   "chemistry": "aromatic amide, logP 4.5–5.5, low fsp3",
   "beneficial": "LOW–MOD for adult bees; larval beneficials affected",
   "flag": "MOD",
   "note": "Moult inhibitors spare adults."
  }
 ],
 "ecorisk": [
  {
   "rule": "logP ≥ 4",
   "flag": "Bioconcentration / soil sorption",
   "level": "HIGH"
  },
  {
   "rule": "logP ≤ 2 and soluble",
   "flag": "Aquatic mobility / runoff",
   "level": "MOD"
  },
  {
   "rule": "Neuroactive pharmacophore",
   "flag": "Non-target neurotoxicity (bees, aquatic arthropods)",
   "level": "HIGH"
  },
  {
   "rule": "Repellent monoterpene tip",
   "flag": "Short DT50, low residue",
   "level": "LOW"
  },
  {
   "rule": "Phenolic / acidic group",
   "flag": "Phytotoxic at high dose on young tissue",
   "level": "MOD"
  },
  {
   "rule": "≥3 aromatic rings",
   "flag": "Persistence / metabolism burden",
   "level": "MOD"
  }
 ],
 "sitepanels": [
  {
   "site": "AChE gorge",
   "note": "Catalytic triad at base of ~20 Å deep gorge; carbamates block gorge entry, OPs phosphorylate the triad Ser."
  },
  {
   "site": "VGSC (kdr area)",
   "note": "Hydrophobic fenvalerate/pyrethroid site at IIS4–S5; T929I / L1014F reduce binding."
  },
  {
   "site": "nAChR agonist pocket",
   "note": "Orthosteric α/β interface; nitroimine/neonic anchors maximise electrostatics."
  },
  {
   "site": "CYP51 haem pocket",
   "note": "Triazole N coordinates haem Fe; pocket volume gates sterol access."
  },
  {
   "site": "PSII QB pocket",
   "note": "Triazine/urea occlude plastoquinone QB, blocking electron transfer."
  },
  {
   "site": "Note on cryptic pockets",
   "note": "Cryptic/induced-fit pockets need MD sampling (AgroMD) — static maps may miss them."
  }
 ],
 "targetClasses": [
  {
   "id": "mite",
   "label": "Pests — mites & ticks",
   "hint": "AChE · VGSC · chitin synthase"
  },
  {
   "id": "insect",
   "label": "Pests — insects",
   "hint": "nAChR · VGSC · GABA-Rdl · RyR"
  },
  {
   "id": "nematode",
   "label": "Nematodes",
   "hint": "nicotinic AChR (levamisole site)"
  },
  {
   "id": "fungus",
   "label": "Pathogens — fungi",
   "hint": "SDH · CYP51 · β-tubulin (MBC)"
  },
  {
   "id": "weed",
   "label": "Weeds",
   "hint": "PSII-D1 · EPSPS"
  },
  {
   "id": "defense",
   "label": "Detox / defense enzymes",
   "hint": "GST — synergist targets"
  }
 ]
};

window.AGROX.LIT = [
  {
    "focus": "Mortality bioassay",
    "status": "Reported means on file",
    "detail": "240 aggregate bioassay records (larval + adult) across 4 extracts, 2 Rhipicephalus spp., 5 concentrations, 2 exposure regimes. Replicate-level data pending.",
    "verify": "See Evidence tables; duplicate wet-lab runs recommended."
  },
  {
    "focus": "Docking (top compounds)",
    "status": "Reported re-docking values",
    "detail": "10 top-scoring candidates per extract with reported affinity and interaction notes; validated against 7 PDB complexes (RMSD).",
    "verify": "Re-run AutoDock Vina + consensus scoring."
  },
  {
    "focus": "AChE homology models",
    "status": "Reported QA scores",
    "detail": "Model generators, Verify3D/ERRAT/ProSA-Z/Ramachandran reported for both tick species.",
    "verify": "Publish coordinate files; deposit on ModelArchive."
  },
  {
    "focus": "LC50 / LC90",
    "status": "Reported probit outputs",
    "detail": "Probit slope, 95% CI and chi-square df reported per extract/species.",
    "verify": "Replicate probit fits on pooled raw counts."
  },
  {
    "focus": "Interactive surrogates",
    "status": "Demo-grade",
    "detail": "In-browser affinity/synergy/eco surrogates are hypotheses for pipeline triage — not validated endpoints.",
    "verify": "Replace with real docking/MD/experiments in production."
  }
];
