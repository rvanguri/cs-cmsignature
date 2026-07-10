"""Shared helpers for the HeartMap-CS pipeline scripts.

Keeps logging, gene panels, and small IO utilities in one place so the step scripts stay short.
All step scripts import from here: `from _common import log, PANELS, ...`.
"""
from __future__ import annotations
import logging
import os
import sys

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(name)s | %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)


def log(name: str) -> logging.Logger:
    return logging.getLogger(name)


def ensure_dir(path: str) -> str:
    os.makedirs(path, exist_ok=True)
    return path


def require(path: str, what: str = "input") -> str:
    """Fail fast with a clear message if a required input is missing."""
    if not os.path.exists(path):
        sys.exit(f"ERROR: required {what} not found: {path}")
    return path


# ---- the harmonized cell-type marker set (Plan Step 3; B split from T/NK) ----
CELLTYPE_MARKERS = {
    "Cardiomyocyte": ["MYH7", "TNNT2", "ACTC1", "MYL2", "TNNI3", "MYBPC3", "TTN", "RYR2"],
    "Fibroblast": ["PDGFRA", "COL1A2", "DCN"],
    "Endothelial": ["PECAM1", "FLT1", "VWF"],
    "Pericyte": ["RGS5", "PDGFRB"],
    "SmoothMuscle": ["MYH11", "TAGLN", "ACTA2"],
    "Macrophage": ["CD163", "MRC1", "C1QA"],
    "Bcell": ["EBF1", "MS4A1", "CD79A"],
    "TNK": ["CD3E", "IL7R", "THEMIS", "NKG7"],
    "Adipocyte": ["GPAM", "ADIPOQ"],
    "Neuronal": ["NRXN1", "PLP1"],
}

# ---- panels used by the two gates and the DE filters (Plan Steps 6b/6c/7) ----
PANELS = {
    # 6b decontamination QC gate: non-CM lineage panel (should be ~0 in CS-CMs post-decontam)
    "lineage_ambient": {
        "Bcell": ["EBF1", "MS4A1"],
        "TNK": ["CD3E", "IL7R", "THEMIS"],
        "Macrophage": ["CD163", "MRC1", "C1QA"],
        "Endothelial": ["PECAM1", "FLT1"],
        "Fibroblast": ["PDGFRA", "COL1A2"],
    },
    # Canonical CS-cardiomyocyte panel (v3). TNNI3K/GJB7 were the original two-gene core; MLIP and
    # PANK1 are RE-ADMITTED because they carry the batch-clean support the core genes lack in one
    # respect or another: MLIP/PANK1 replicate in the within-study CS-vs-ICM (Liu, 5') contrast (both
    # q<0.05), whereas GJB7/TNNI3K are untestable there (<1% 5' capture). All four are recovered in the
    # independent Foong spatial data and clear every artifact gate below (not IEG/HSP, not
    # lineage-ambient, dissociation-stress program not elevated). The cohort-confounded discovery
    # contrasts are NOT the basis for membership.
    "cs_cm_expected_up": ["TNNI3K", "GJB7", "MLIP", "PANK1"],
    # NOT expected to validate (artifacts of decontam-inflated CS-vs-DCM / mismatched arms)
    "cs_cm_not_expected": ["SPECC1L", "CALM2", "TJP2"],
    # CM-dedifferentiation panel = TECHNICAL CONTROL ONLY, not CS biology
    "cm_dediff_control": ["TNNT2", "MYL2", "MB", "CKB"],
    # ischemia/agonal IEG/HSP panel to flag/exclude (surgical handling artifact)
    "ieg_hsp": ["FOS", "FOSB", "JUN", "JUNB", "EGR1", "ATF3",
                "HSPA1A", "HSPA1B", "DNAJB1", "NR4A1"],
    # HeartMap cross-study DCM (generic-HF / failing-axis) signature
    "generic_hf_up": ["NPPA", "NPPB", "C5AR1", "PLCE1", "UNC80"],
    "generic_hf_down": ["BMP7"],
    # activated-fibroblast niche markers (HeartMap)
    "fib_niche": {"COL22A1_DCM": ["COL22A1"], "TNC_ICM": ["TNC"]},
}

# the five 3' datasets that form the integrated arm
THREE_PRIME_DATASETS = ["neyazi", "larson", "chin2022", "chin2021", "reichart"]

# ---- disease vocabulary harmonization (one shared label set across all datasets) ----
# CELLxGENE (Reichart) gives ontology strings; GEO sample_meta gives short codes. Normalize both to:
CANONICAL_DISEASES = {"CS", "DCM", "ICM", "ARVC", "HCM", "NF", "NCC"}
DISEASE_ALIASES = {
    "cardiac sarcoidosis": "CS", "sarcoidosis": "CS", "cs": "CS",
    "dilated cardiomyopathy": "DCM", "dcm": "DCM", "idiopathic dilated cardiomyopathy": "DCM",
    "ischemic cardiomyopathy": "ICM", "ischaemic cardiomyopathy": "ICM",
    "ischemic heart disease": "ICM", "icm": "ICM",
    "arrhythmogenic right ventricular cardiomyopathy": "ARVC",
    "arrhythmogenic cardiomyopathy": "ARVC", "acm": "ARVC", "arvc": "ARVC",
    "hypertrophic cardiomyopathy": "HCM", "hcm": "HCM",
    "normal": "NF", "non-failing": "NF", "nonfailing": "NF", "control": "NF",
    "healthy": "NF", "donor": "NF", "nf": "NF",
    "non-compaction cardiomyopathy": "NCC", "left ventricular non-compaction": "NCC", "ncc": "NCC",
}


def normalize_disease(val):
    """Map a disease label (ontology string or short code) to the canonical vocab.
    Unrecognized values are returned unchanged (caller should warn so they get caught)."""
    if val is None:
        return val
    return DISEASE_ALIASES.get(str(val).strip().lower(), str(val))


# ---- cell-type harmonization: CELLxGENE Cell-Ontology labels -> the CELLTYPE_MARKERS vocab ----
# so the scANVI reference (Reichart) uses the same names the gates/DE key on ("Cardiomyocyte" etc.).
CELLTYPE_ALIASES = {
    "cardiac muscle cell": "Cardiomyocyte", "cardiomyocyte": "Cardiomyocyte",
    "fibroblast of cardiac tissue": "Fibroblast", "fibroblast": "Fibroblast",
    "endothelial cell": "Endothelial",
    "mural cell": "Pericyte", "pericyte": "Pericyte",
    "smooth muscle cell": "SmoothMuscle",
    "myeloid cell": "Macrophage", "macrophage": "Macrophage", "monocyte": "Macrophage",
    "adipocyte": "Adipocyte",
    "cardiac neuron": "Neuronal", "neuron": "Neuronal", "glial cell": "Neuronal",
    "lymphocyte": "TNK", "t cell": "TNK", "nk cell": "TNK", "natural killer cell": "TNK",
    "b cell": "Bcell", "plasma cell": "Bcell",
    "mast cell": "Mast",
    "unknown": "Unknown", "nan": "Unknown",
}


def normalize_celltype(val):
    """Map a cell-type label to the CELLTYPE_MARKERS vocab (exact, then substring). Fallback: unchanged."""
    if val is None:
        return "Unknown"
    k = str(val).strip().lower()
    if k in CELLTYPE_ALIASES:
        return CELLTYPE_ALIASES[k]
    for key, v in CELLTYPE_ALIASES.items():
        if key in k:                      # e.g. "endothelial cell of coronary artery" -> Endothelial
            return v
    return str(val)
