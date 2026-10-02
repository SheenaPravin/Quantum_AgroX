"""
Quantum_AgroX™ – Streamlit app
Modules: AgroPhytoX (enter 2 phytochemicals) -> AgroSyntheticX (draw molecule
         -> SMILES) -> AgroDockX (SMILES handoff) -> AgroReport (parameters
         grouped under sub-modules only).

Run:  streamlit run app.py   (from Quantum_AgroX_QML_Project/)
"""

from pathlib import Path

import pandas as pd
import streamlit as st
import yaml

try:
    from streamlit_ketcher import st_ketcher

    HAS_KETCHER = True
except Exception:
    st_ketcher = None
    HAS_KETCHER = False

BASE_DIR = Path(__file__).resolve().parent
RAW_DIR = BASE_DIR / "data" / "raw"

import sys

sys.path.insert(0, str(BASE_DIR / "src"))
from phytochemical_library import rdkit_descriptors, resolve_smiles  # noqa: E402

st.set_page_config(page_title="Quantum_AgroX", layout="wide")

# ---------------------------------------------------------------- session ---
DEFAULT_A = "Phytol"
DEFAULT_B = "alpha-Caryophyllene"
for k, v in {
    "phyto_a_name": DEFAULT_A,
    "phyto_b_name": DEFAULT_B,
    "smiles_a": "",
    "smiles_b": "",
    "smiles_a_src": "",
    "smiles_b_src": "",
    "dock_label_a": "",
    "dock_label_b": "",
    "synth_pending": "",
    "synth_smiles_a": "",
    "synth_smiles_b": "",
    "synth_name_a": "Synthetic-A",
    "synth_name_b": "Synthetic-B",
    "target_species": "R. microplus",
    "try_pubchem": False,
}.items():
    st.session_state.setdefault(k, v)


def dock_labels() -> tuple[str, str]:
    """Names shown in AgroDockX: explicit handoff labels, else PhytoX pair."""
    la = st.session_state.get("dock_label_a") or st.session_state["phyto_a_name"]
    lb = st.session_state.get("dock_label_b") or st.session_state["phyto_b_name"]
    return la, lb


# ---------------------------------------------------------------- loaders ---
@st.cache_data
def load_csv(name: str) -> pd.DataFrame:
    p = RAW_DIR / name
    return pd.read_csv(p) if p.exists() else pd.DataFrame()


@st.cache_data
def load_config() -> dict:
    p = BASE_DIR / "config" / "project_config.yaml"
    if p.exists():
        with open(p) as f:
            return yaml.safe_load(f)
    return {}


phyto_df = load_csv("01_major_phytochemicals_reported.csv")
dock_df = load_csv("05_docking_top_compounds_reported.csv")
params_df = load_csv("08_AChE_model_and_docking_parameters.csv")
config = load_config()

_phyto_names = phyto_df["Compound"].dropna().tolist() if "Compound" in phyto_df else []
_dock_names = dock_df["Compound"].dropna().tolist() if "Compound" in dock_df else []
ALL_COMPOUNDS = sorted(set(_phyto_names + _dock_names)) or [DEFAULT_A, DEFAULT_B]


def phyto_metadata(name: str) -> pd.DataFrame:
    if phyto_df.empty or not name:
        return pd.DataFrame()
    return phyto_df[phyto_df["Compound"].str.lower() == name.lower()]


def docking_hit(name: str) -> pd.DataFrame:
    if dock_df.empty or not name:
        return pd.DataFrame()
    # exact match first, else short-name containment for truncated paper names
    exact = dock_df[dock_df["Compound"].str.lower() == name.lower()]
    if not exact.empty:
        return exact
    short = name.split(",")[0].split("(")[0].strip().lower()
    if len(short) > 8:
        return dock_df[dock_df["Compound"].str.lower().str.contains(short[:20])]
    return pd.DataFrame()


def ensure_smiles(slot: str, name: str) -> tuple[str, str]:
    """Resolve + cache SMILES for slot 'a'/'b'; returns (smiles, source)."""
    smi_key, src_key = f"smiles_{slot}", f"smiles_{slot}_src"
    # re-resolve when the name changed or nothing cached yet
    if not st.session_state[smi_key] or st.session_state.get(f"resolved_for_{slot}") != name:
        smi, src = resolve_smiles(name, try_pubchem=st.session_state["try_pubchem"])
        st.session_state[smi_key], st.session_state[src_key] = smi, src
        st.session_state[f"resolved_for_{slot}"] = name
    return st.session_state[smi_key], st.session_state[src_key]


def mol_image(smiles: str, size=(380, 260)):
    try:
        from rdkit import Chem
        from rdkit.Chem import Draw

        mol = Chem.MolFromSmiles(smiles) if smiles else None
        return Draw.MolToImage(mol, size=size) if mol is not None else None
    except Exception:
        return None


# ================================================================== AgroPhytoX
def page_phytosx():
    st.header("AgroPhytoX™ — Phytochemical entry (2 compounds)")
    st.caption(
        "Enter two phytochemical names. They are stored in the session and "
        "their SMILES are automatically carried over to AgroDockX."
    )
    c1, c2 = st.columns(2)
    with c1:
        st.subheader("Compound A")
        a = st.selectbox(
            "Select / search Compound A",
            ALL_COMPOUNDS,
            index=ALL_COMPOUNDS.index(st.session_state["phyto_a_name"])
            if st.session_state["phyto_a_name"] in ALL_COMPOUNDS
            else 0,
            key="sel_a",
        )
        a_custom = st.text_input(
            "…or type a custom name for A", value="", placeholder="e.g. Phytol"
        )
        if a_custom.strip():
            a = a_custom.strip()
        st.session_state["phyto_a_name"] = a
    with c2:
        st.subheader("Compound B")
        b = st.selectbox(
            "Select / search Compound B",
            ALL_COMPOUNDS,
            index=ALL_COMPOUNDS.index(st.session_state["phyto_b_name"])
            if st.session_state["phyto_b_name"] in ALL_COMPOUNDS
            else 1,
            key="sel_b",
        )
        b_custom = st.text_input(
            "…or type a custom name for B",
            value="",
            placeholder="e.g. alpha-Caryophyllene",
        )
        if b_custom.strip():
            b = b_custom.strip()
        st.session_state["phyto_b_name"] = b

    st.checkbox(
        "Try live PubChem lookup for names missing from the curated library (needs internet)",
        key="try_pubchem",
    )
    # pre-resolve so AgroDockX opens with SMILES ready
    smi_a, src_a = ensure_smiles("a", a)
    smi_b, src_b = ensure_smiles("b", b)

    for label, name, smi, src in [("A", a, smi_a, src_a), ("B", b, smi_b, src_b)]:
        st.markdown(f"### Compound {label}: {name}")
        meta = phyto_metadata(name)
        if not meta.empty:
            st.dataframe(meta, use_container_width=True)
        else:
            st.info("No GC-MS metadata row for this exact name (custom entry).")
        if smi:
            st.success(f"SMILES resolved ({src}): `{smi}`")
        else:
            st.warning(
                "No curated/PubChem SMILES found — paste it manually in AgroDockX."
            )
    if st.button("Send both to AgroDockX →", type="primary"):
        st.session_state["dock_label_a"] = a
        st.session_state["dock_label_b"] = b
        st.session_state["nav"] = "AgroDockX"
        st.rerun()


# ============================================================ AgroSyntheticX
def show_synthetic_result(mol) -> str:
    """Render canonical SMILES + 2D image + descriptors; stage as pending."""
    from rdkit import Chem
    from rdkit.Chem import Draw

    smiles = Chem.MolToSmiles(mol, canonical=True)
    st.session_state["synth_pending"] = smiles
    st.success(f"**SMILES:** `{smiles}`")
    st.subheader("2D Structure")
    st.image(Draw.MolToImage(mol, size=(300, 300)))
    desc = rdkit_descriptors(smiles)
    if desc.get("valid"):
        st.table(pd.DataFrame([desc]).T.rename(columns={0: "value"}))
    return smiles


def page_syntheticx():
    st.header("AgroSyntheticX™ — Draw molecule → SMILES")
    st.caption(
        "Draw a synthetic candidate below. Its SMILES can be sent to AgroDockX "
        "as Molecule A or B (drawn with Ketcher, canonicalised with RDKit)."
    )
    if not HAS_KETCHER:
        st.error(
            "The `streamlit-ketcher` package is not installed. "
            "Run `pip install streamlit-ketcher` to enable drawing."
        )
        return
    molecule = st_ketcher(key="molecule_drawing")
    if molecule:
        try:
            from rdkit import Chem

            mol = None
            if isinstance(molecule, str):
                # Case 1: Ketcher returns SMILES (most common)
                mol = Chem.MolFromSmiles(molecule)
            else:
                # Case 2: Ketcher returns MolBlock
                mol = Chem.MolFromMolBlock(molecule)
            if mol is not None:
                show_synthetic_result(mol)
            else:
                st.error("Invalid molecule structure. Please redraw.")
        except Exception as e:
            st.error(f"Error processing molecule: {e}")
    else:
        st.info("Draw a molecule in the editor above ⬆️")

    st.markdown("**…or paste a SMILES directly** (fallback if the board fails to load)")
    pasted = st.text_area(
        "SMILES input",
        value="",
        placeholder="e.g. CC(=O)Oc1ccccc1C(=O)O",
        key="synth_paste",
    )
    if st.button("Convert → canonical SMILES"):
        try:
            from rdkit import Chem

            mol = Chem.MolFromSmiles(pasted.strip()) if pasted.strip() else None
            if mol is not None:
                show_synthetic_result(mol)
            else:
                st.error("Invalid SMILES — RDKit cannot parse it.")
        except Exception as e:
            st.error(f"Error processing SMILES: {e}")
    smiles = st.session_state.get("synth_pending", "")
    if smiles:
        st.markdown("### Store this drawing separately as Synthetic A / B")
        c1, c2 = st.columns(2)
        with c1:
            st.text_input("Name for Synthetic A", key="synth_name_a")
            if st.button("Store as Synthetic A"):
                st.session_state["synth_smiles_a"] = smiles
                st.success(f"Stored Synthetic A: `{smiles}`")
        with c2:
            st.text_input("Name for Synthetic B", key="synth_name_b")
            if st.button("Store as Synthetic B"):
                st.session_state["synth_smiles_b"] = smiles
                st.success(f"Stored Synthetic B: `{smiles}`")
    sa, sb = st.session_state.get("synth_smiles_a", ""), st.session_state.get(
        "synth_smiles_b", ""
    )
    if sa or sb:
        st.markdown("### Stored synthetics")
        st.table(
            pd.DataFrame(
                [
                    {
                        "slot": "A",
                        "name": st.session_state.get("synth_name_a", ""),
                        "SMILES": sa or "—",
                    },
                    {
                        "slot": "B",
                        "name": st.session_state.get("synth_name_b", ""),
                        "SMILES": sb or "—",
                    },
                ]
            )
        )
        if st.button("Send both to AgroDockX →", type="primary"):
            for s, smi, name in [
                ("a", sa, st.session_state.get("synth_name_a", "")),
                ("b", sb, st.session_state.get("synth_name_b", "")),
            ]:
                if smi:
                    st.session_state[f"smiles_{s}"] = smi
                    st.session_state[f"smiles_{s}_src"] = "drawn (AgroSyntheticX)"
                    st.session_state[f"resolved_for_{s}"] = name
                    st.session_state[f"dock_label_{s}"] = name
            st.session_state["nav"] = "AgroDockX"
            st.rerun()
    st.markdown("---")
    st.caption("Powered by [RDKit](https://www.rdkit.org/) + [Ketcher](https://lifescience.opensource.epam.com/ketcher/)")


# ================================================================== AgroDockX
def take_synthetic(slot: str):
    """Pull stored Synthetic A/B from AgroSyntheticX into an AgroDockX slot."""
    smi = st.session_state.get(f"synth_smiles_{slot}", "")
    name = st.session_state.get(f"synth_name_{slot}", "")
    if smi:
        st.session_state[f"smiles_{slot}"] = smi
        st.session_state[f"smiles_{slot}_src"] = "drawn (AgroSyntheticX)"
        st.session_state[f"resolved_for_{slot}"] = name
        st.session_state[f"dock_label_{slot}"] = name
        st.success(f"Synthetic {slot.upper()} taken into slot {slot.upper()} ✓")
    else:
        st.warning(f"No Synthetic {slot.upper()} stored yet — draw it first.")


def compound_dock_card(slot: str, name: str):
    st.markdown(f"### Compound {slot.upper()}: {name}")
    smi, src = ensure_smiles(slot, name)
    smi_key = f"smiles_{slot}"
    edited = st.text_input(
        f"SMILES for {slot.upper()} (auto-filled from AgroPhytoX, editable)",
        value=smi,
        key=f"smi_edit_{slot}",
    )
    if edited != smi:
        st.session_state[smi_key] = edited
        st.session_state[f"smiles_{slot}_src"] = "manual-edit"
        smi = edited
        src = st.session_state[f"smiles_{slot}_src"]
    st.caption(f"SMILES source: {src or '—'}")
    if not smi:
        st.warning("Enter a SMILES to enable structure + descriptors.")
        return smi
    img = mol_image(smi)
    if img is not None:
        st.image(img, caption=f"{name} — 2D structure")
    else:
        st.error("Invalid SMILES — RDKit cannot parse it.")
    desc = rdkit_descriptors(smi)
    if desc.get("valid"):
        st.table(pd.DataFrame([desc]).T.rename(columns={0: "value"}))
    hit = docking_hit(name)
    if not hit.empty:
        st.markdown("**Reported docking hit(s)**")
        st.dataframe(hit, use_container_width=True)
    return smi


def page_dockx():
    st.header("AgroDockX™ — SMILES handoff + docking setup")
    st.caption(
        "SMILES below are auto-filled from the AgroPhytoX pair and/or the "
        "AgroSyntheticX drawing (edit in place if needed)."
    )
    a, b = dock_labels()
    pa, pb = st.session_state["phyto_a_name"], st.session_state["phyto_b_name"]
    st.info(f"AgroPhytoX selection → **A:** {pa} · **B:** {pb}")
    if st.session_state.get("synth_smiles_a") or st.session_state.get(
        "synth_smiles_b"
    ):
        st.info(
            f"AgroSyntheticX stored → **A:** `{st.session_state.get('synth_smiles_a', '') or '—'}` "
            f"({st.session_state.get('synth_name_a', '')}) · **B:** `{st.session_state.get('synth_smiles_b', '') or '—'}` "
            f"({st.session_state.get('synth_name_b', '')})"
        )
        t1, t2 = st.columns(2)
        with t1:
            if st.button("Take Synthetic A → slot A"):
                take_synthetic("a")
        with t2:
            if st.button("Take Synthetic B → slot B"):
                take_synthetic("b")
    if st.session_state.get("mix_note"):
        st.info(st.session_state.pop("mix_note"))
    if st.button("Reset both slots to AgroPhytoX pair"):
        for s, name in [("a", pa), ("b", pb)]:
            st.session_state[f"smiles_{s}"] = ""
            st.session_state[f"smiles_{s}_src"] = ""
            st.session_state[f"resolved_for_{s}"] = ""
            st.session_state[f"dock_label_{s}"] = name
        st.rerun()
    st.markdown("**Mixed pair — e.g. phytochemical as A + synthetic as B**")
    m1, m2 = st.columns(2)
    with m1:
        mix_a = st.selectbox(
            "Molecule A from",
            ["Keep current", "AgroPhytoX A", "Synthetic A"],
            key="mix_a",
        )
    with m2:
        mix_b = st.selectbox(
            "Molecule B from",
            ["Keep current", "AgroPhytoX B", "Synthetic B"],
            key="mix_b",
        )
    if st.button("Build mixed pair →"):
        notes = []
        for slot, choice in [("a", mix_a), ("b", mix_b)]:
            if choice == "Keep current":
                notes.append(f"{slot.upper()}: kept")
                continue
            if choice.startswith("AgroPhytoX"):
                name = pa if slot == "a" else pb
                smi, _ = ensure_smiles(slot, name)
                if smi:
                    st.session_state[f"dock_label_{slot}"] = name
                    notes.append(f"{slot.upper()}: {name} ✓")
                else:
                    notes.append(f"{slot.upper()}: {name} has no SMILES")
                continue
            smi = st.session_state.get(f"synth_smiles_{slot}", "")
            name = st.session_state.get(f"synth_name_{slot}", "")
            if smi:
                st.session_state[f"smiles_{slot}"] = smi
                st.session_state[f"smiles_{slot}_src"] = "drawn (AgroSyntheticX)"
                st.session_state[f"resolved_for_{slot}"] = name
                st.session_state[f"dock_label_{slot}"] = name
                notes.append(f"{slot.upper()}: {name} ✓")
            else:
                notes.append(f"{slot.upper()}: no Synthetic {slot.upper()} stored")
        st.session_state["mix_note"] = "Mixed pair — " + " · ".join(notes)
        st.rerun()
    st.selectbox(
        "Target species (AChE model)",
        ["R. decoloratus", "R. microplus"],
        key="target_species",
    )
    c1, c2 = st.columns(2)
    with c1:
        compound_dock_card("a", a)
    with c2:
        compound_dock_card("b", b)

    st.subheader("Docking parameters (from 08_AChE_model_and_docking_parameters.csv)")
    if not params_df.empty:
        st.dataframe(
            params_df[params_df["Species"] == st.session_state["target_species"]],
            use_container_width=True,
        )
    else:
        st.warning("Docking parameter file not found.")
    if st.button("Generate AgroReport →", type="primary"):
        st.session_state["nav"] = "AgroReport"
        st.rerun()


# ================================================================ AgroReport
def kv_table(d: dict) -> pd.DataFrame:
    return pd.DataFrame([{"parameter": k, "value": v} for k, v in d.items()])


def page_report():
    st.header("AgroReport™ — parameters grouped by sub-module")
    a, b = st.session_state["phyto_a_name"], st.session_state["phyto_b_name"]
    la, lb = dock_labels()
    smi_a = st.session_state.get("smiles_a", "")
    smi_b = st.session_state.get("smiles_b", "")
    tgt = st.session_state.get("target_species", "")
    q = config.get("quantum", {})
    ev = config.get("evaluation", {})

    sections: dict[str, pd.DataFrame] = {}
    dock_used = bool(smi_a or smi_b)
    synth_used = bool(
        st.session_state.get("synth_smiles_a") or st.session_state.get("synth_smiles_b")
    )

    # Summary — brief
    sections["Summary — brief"] = kv_table(
        {
            "Pair": f"{a} + {b}",
            "Synthetic_A": f"{st.session_state.get('synth_name_a', '')}: "
            f"{st.session_state.get('synth_smiles_a', '') or '—'}",
            "Synthetic_B": f"{st.session_state.get('synth_name_b', '')}: "
            f"{st.session_state.get('synth_smiles_b', '') or '—'}",
            "AgroDockX": f"A: {smi_a or '—'} · B: {smi_b or '—'}",
            "Target": tgt or "—",
            "Caveat": "Pre-screening hypotheses — docking scores are not proof "
            "of inhibition.",
        }
    )
    if not dock_used and not synth_used:
        sections["Note"] = kv_table(
            {
                "next": "Pick a pair in AgroPhytoX, draw in AgroSyntheticX, or "
                "fill AgroDockX to populate the sections below."
            }
        )
    synth_a_smi = st.session_state.get("synth_smiles_a", "") or "—"
    synth_b_smi = st.session_state.get("synth_smiles_b", "") or "—"
    sections["1 · AgroPhytoX — inputs"] = kv_table(
        {
            "compound_A": a,
            "compound_B": b,
            "synthetic_A": f"{st.session_state.get('synth_name_a', '')}: "
            f"{synth_a_smi}",
            "synthetic_B": f"{st.session_state.get('synth_name_b', '')}: "
            f"{synth_b_smi}",
        }
    )

    # 2. AgroTargetX / AgroSiteMap — only when docking is in play
    if dock_used:
        if not params_df.empty:
            row = params_df[params_df["Species"] == tgt]
            if not row.empty:
                r0 = row.iloc[0].to_dict()
                sections["2 · AgroTargetX / AgroSiteMap — target & pocket"] = kv_table(r0)
            else:
                sections["2 · AgroTargetX / AgroSiteMap — target & pocket"] = kv_table(
                    {"species": tgt, "note": "no parameter row for species"}
                )
        else:
            sections["2 · AgroTargetX / AgroSiteMap — target & pocket"] = kv_table(
                {"species": tgt}
            )

    # 3. AgroDockX — only slots holding a SMILES
    for slot, name, smi in [("A", la, smi_a), ("B", lb, smi_b)]:
        if not smi:
            continue
        d = rdkit_descriptors(smi)
        hit = docking_hit(name)
        aff = (
            hit.iloc[0]["Binding_affinity_kcal_mol"]
            if not hit.empty and "Binding_affinity_kcal_mol" in hit
            else "—"
        )
        sections[f"3 · AgroDockX — compound {slot} ({name})"] = kv_table(
            {
                "SMILES": smi or "—",
                "SMILES_source": st.session_state.get(f"smiles_{slot.lower()}_src", ""),
                "MW": d.get("mw", "—"),
                "LogP": d.get("logp", "—"),
                "reported_affinity_kcal/mol": aff,
            }
        )
    if dock_used:
        sections["3 · AgroDockX — run parameters"] = kv_table(
        {"target_species": tgt, "vina_exhaustiveness": 8, "vina_modes": 9}
    )

    # 4. AgroDoseX / AgroQML — only when docking is in play
    if dock_used:
        sections["4 · AgroDoseX / AgroQML — model parameters"] = kv_table(
            {
                "n_qubits": q.get("n_qubits", "—"),
                "feature_map": q.get("feature_map", "—"),
                "shots": q.get("shots", "—"),
                "cv_folds": ev.get("cv_folds", "—"),
                "classification_threshold": config.get("classification_threshold", {}).get(
                    "high_efficacy", "—"
                ),
            }
        )

    for title, df in sections.items():
        st.subheader(title)
        st.table(df)

    # Markdown export
    md = ["# Quantum_AgroX — AgroReport", ""]
    for title, df in sections.items():
        md.append(f"## {title}")
        for _, r in df.iterrows():
            md.append(f"- **{r['parameter']}**: {r['value']}")
        md.append("")
    st.download_button(
        "Download report (.md)",
        "\n".join(md),
        file_name="agroreport.md",
        mime="text/markdown",
    )


# ------------------------------------------------------------------ nav -----
st.sidebar.title("Quantum_AgroX™")
NAV_MODULES = ["AgroPhytoX", "AgroSyntheticX", "AgroDockX", "AgroReport"]
nav = st.sidebar.radio(
    "Module",
    NAV_MODULES,
    index=NAV_MODULES.index(st.session_state.get("nav", "AgroPhytoX")),
    key="nav",
)
st.sidebar.divider()
st.sidebar.caption(
    f"PhytoX: A={st.session_state['phyto_a_name']} · "
    f"B={st.session_state['phyto_b_name']}"
)
if st.session_state.get("synth_smiles_a") or st.session_state.get("synth_smiles_b"):
    st.sidebar.caption(
        f"Synthetic A: `{st.session_state.get('synth_smiles_a', '') or '—'}` · "
        f"B: `{st.session_state.get('synth_smiles_b', '') or '—'}`"
    )

if nav == "AgroPhytoX":
    page_phytosx()
elif nav == "AgroSyntheticX":
    page_syntheticx()
elif nav == "AgroDockX":
    page_dockx()
else:
    page_report()
