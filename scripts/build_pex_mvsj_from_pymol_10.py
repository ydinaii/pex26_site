"""
Convert the PyMOL PEX1/PEX6/PEX26 coloring + inset script into a
two-snapshot MolViewSpec (MVSJ) file.

Adjust STRUCTURE_URL to point at wherever the CIF is hosted on your
Quarto site (relative path, e.g. "data/hs_7549_db_2026...cif"), since
the MVS builder loads structures via download(url=...) rather than a
local filesystem path.
"""

from molviewspec import create_builder, ComponentExpression
from molviewspec.builder import States, GlobalMetadata

# ---------------- CONFIG ----------------
# Absolute URL (matches the pattern used by the working pex26_complex.mvsj on this site —
# relative paths resolve against the MVSJ file's own location, which caused issues before)
STRUCTURE_URL = "https://ydinaii.github.io/pex26_site/data/hs_7549_db_2026_force_0_permuted_minimized_run_1.cif"

PEX1_CHAINS  = ["A", "B", "C"]
PEX6_CHAINS  = ["D", "E", "F"]
PEX26_CHAINS = ["G", "H", "I"]

COLOR_PEX1  = "#0072B2"   # blue
COLOR_PEX6  = "#E69F00"   # gold
COLOR_PEX26 = "#009E73"   # bluish-green

# Residues labeled in the inset (chain, resi, one-letter aa code)
INSET_LABEL_RESIDUES = [
    ("G", 51, "F"),   # Phe
    ("G", 49, "L"),   # Leu
    ("G", 95, "E"),   # Glu
    ("D", 218, "F"),  # Phe
    ("A", 229, "N"),  # Asn
]
# -----------------------------------------


def chain_selector(chains):
    """Build a union selector (list of ComponentExpression) for a list of chain IDs."""
    return [ComponentExpression(auth_asym_id=c) for c in chains]


def build_overview_snapshot():
    builder = create_builder()
    structure = (
        builder.download(url=STRUCTURE_URL)
        .parse(format="mmcif")
        .model_structure()
    )

    # --- PEX1 / PEX6 / PEX26 colored by chain group ---
    structure.component(selector=chain_selector(PEX1_CHAINS)) \
        .representation(type="cartoon").color(color=COLOR_PEX1)

    structure.component(selector=chain_selector(PEX6_CHAINS)) \
        .representation(type="cartoon").color(color=COLOR_PEX6)

    structure.component(selector=chain_selector(PEX26_CHAINS)) \
        .representation(type="cartoon").color(color=COLOR_PEX26)

    # --- Ligands / ions (ATP, ADP, PO4, Mg2+, etc.) as sticks ---
    # Explicit label_comp_id list rather than the generic "ligand"/"ion" selector —
    # AF3-generated CIFs often don't populate the non-polymer categorization fields
    # that the generic MVS keyword relies on, so it silently matches nothing.
    # Matches the proven pattern in build_pex_mvsj.py.
    LIGAND_RESN = ["ATP", "ADP", "PO4", "NA", "K", "MG", "CA", "CL"]
    ligand_selector = [ComponentExpression(label_comp_id=r) for r in LIGAND_RESN]

    ligands = structure.component(selector=ligand_selector)
    ligands.representation(type="ball_and_stick").color(custom={"molstar_color_theme_name": "element-symbol"})

    # --- F51 label, visible only in this (overview) snapshot ---
    f51 = structure.component(selector=ComponentExpression(auth_asym_id="G", auth_seq_id=51))
    f51.representation(type="ball_and_stick").color(color="red")
    f51.label(text="F51")

    # Overview camera — fit the whole assembly in frame
    whole = structure.component(selector="all")
    whole.focus()

    return builder.get_snapshot(
        key="cartoon_view",
        title="Overview",
        description="PEX1 (blue) / PEX6 (gold) / PEX26 (green) hexameric complex, F51 labeled.",
    )


def build_inset_snapshot():
    builder = create_builder()
    structure = (
        builder.download(url=STRUCTURE_URL)
        .parse(format="mmcif")
        .model_structure()
    )

    # Keep the same base coloring so the region reads consistently with the overview
    structure.component(selector=chain_selector(PEX1_CHAINS)) \
        .representation(type="cartoon").color(color=COLOR_PEX1)
    structure.component(selector=chain_selector(PEX6_CHAINS)) \
        .representation(type="cartoon").color(color=COLOR_PEX6)
    structure.component(selector=chain_selector(PEX26_CHAINS)) \
        .representation(type="cartoon").color(color=COLOR_PEX26)

    # F51 site + neighboring residues (approximate "byres around 5" with an explicit list;
    # MVS selectors don't support a live spatial-proximity expression, so the neighbor set
    # needs to be enumerated up front — reuse whatever residue list your PyMOL "around 5"
    # selection produced, if you want an exact match).
    site_selector = [
        ComponentExpression(auth_asym_id=chain_id, auth_seq_id=resi)
        for chain_id, resi, _ in INSET_LABEL_RESIDUES
    ]
    site = structure.component(selector=site_selector)
    site.representation(type="ball_and_stick").color(custom={"molstar_color_theme_name": "element-symbol"})

    # Mg2+ near the site, if relevant to this pocket
    mg = structure.component(selector=ComponentExpression(label_comp_id="MG"))
    mg.representation(type="spacefill").color(color="#8FE388")  # pale green, matches PyMOL sphere_color

    # Labels for all five residues (single-letter aa + resi, e.g. "F51")
    for chain_id, resi, aa in INSET_LABEL_RESIDUES:
        res = structure.component(selector=ComponentExpression(auth_asym_id=chain_id, auth_seq_id=resi))
        res.label(text=f"{aa}{resi}")

    # Zoomed camera on the F51 site
    site.focus(radius=10, radius_extent=5)

    return builder.get_snapshot(
        key="inset_view",
        title="F51 site",
        description="Zoomed view of PEX26 F51 and neighboring residues (G49, G95, D218, A229).",
    )


if __name__ == "__main__":
    snapshot_overview = build_overview_snapshot()
    snapshot_inset = build_inset_snapshot()

    states = States(
        snapshots=[snapshot_overview, snapshot_inset],
        metadata=GlobalMetadata(description="PEX1/PEX6/PEX26 complex — overview and F51 inset"),
    )

    with open("pex_complex.mvsj", "w") as f:
        f.write(states.dumps(indent=2))

    print("Wrote pex_complex.mvsj")
