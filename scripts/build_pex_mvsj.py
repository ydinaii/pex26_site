"""
Convert the PEX1/PEX6/PEX26/PEX5 PyMOL coloring script into a MolViewSpec (.mvsj)
file with 3 navigable snapshots: cartoon_view, surface_view, pLDDT_view.

Run:
    pip install molviewspec --break-system-packages
    python build_pex_mvsj.py

Output:
    pex26_complex.mvsj   (upload alongside your structure .cif to pex26_site repo)

IMPORTANT — fill this in before running:
    STRUCTURE_URL must point to a *hosted* .cif file (e.g. raw.githubusercontent.com
    URL to your AF3 model_7549_0 structure, once it's in the repo).
"""

import molviewspec as mvs

# ---------------------------------------------------------------------------
# Config — EDIT THIS
# ---------------------------------------------------------------------------
STRUCTURE_URL = "https://ydinaii.github.io/pex26_site/data/hs_7549_1283x3_980x3_305x1_31x10_27x2_1x10_5x0_4x3_db_2025_12_21_seed-1_sample-0_model.cif"
STRUCTURE_FORMAT = "mmcif"
OUTPUT_PATH = "pex26_complex.mvsj"

# ---------------------------------------------------------------------------
# Chain groups (from your PyMOL script)
# ---------------------------------------------------------------------------
PEX1_CHAINS = ["A", "B", "C"]
PEX6_CHAINS = ["D", "E", "F"]
PEX26_CHAINS = ["G", "H", "I"]
PEX5_CHAINS = ["J"]

# First residue number of each protein IN THIS MODEL'S STRUCTURE FILE.
# This is the equivalent of your `first_resi_of_pex1/6/26/5` variables.
# Domain/motif ranges below are defined in CANONICAL numbering (starting at 1);
# an offset of (first_resi - 1) is added when building each selector, exactly
# like `cmd.alter(..., f"resi=str(int(resi)+{first_resi}-1)")` did in PyMOL.
# For model_7549_0 these are all 1 (no shift). Change per-model as needed.
FIRST_RESI = {
    "PEX1": 1,
    "PEX6": 1,
    "PEX26": 1,
    "PEX5": 1,
}


LIGAND_RESN = ["ATP", "ADP", "PO4", "NA", "K", "MG", "CA", "CL"]
PHOSPHO_RESN = ["SEP", "TPO", "PTR"]

# ---------------------------------------------------------------------------
# Domain / motif residue ranges (verbatim from your PyMOL dicts)
# ---------------------------------------------------------------------------
PEX1_RANGES = {
    "PEX1_N1": "1-180",
    "PEX1_N2": "266-532",
    "PEX1_D1": "560-819",
    "PEX1_D2": "846-1084",
    "PEX1_D1_Walker_A": "599-606",
    "PEX1_D1_alpha_insert": "629-632",
    "PEX1_D1_Walker_B": "658-663",
    "PEX1_D1_Sensor_1": "712-712",
    "PEX1_D2_Walker_A": "881-888",
    "PEX1_D2_alpha_insert": "907-910",
    "PEX1_D2_Pore_loop_1": "914-914",
    "PEX1_D2_Walker_B": "936-941",
    "PEX1_D2_DE": "940-941",
    "PEX1_D2_ISS": "969-971",
    "PEX1_D2_Sensor_1": "983-983",
    "PEX1_D2_R_fingers": "995-998",
}

PEX6_RANGES = {
    "PEX6_N1": "1-184",
    "PEX6_N2": "185-426",
    "PEX6_D1": "439-694",
    "PEX6_D2": "711-957",
    "PEX6_D1_Walker_A": "470-477",
    "PEX6_D1_alpha_insert": "496-499",
    "PEX6_D1_Walker_B": "525-530",
    "PEX6_D1_Sensor_1": "573-573",
    "PEX6_D2_Walker_A": "744-751",
    "PEX6_D2_alpha_insert": "771-773",
    "PEX6_D2_Pore_loop_1": "777-777",
    "PEX6_D2_Walker_B": "799-804",
    "PEX6_D2_DE": "803-804",
    "PEX6_D2_ISS": "834-836",
    "PEX6_D2_Sensor_1": "848-848",
    "PEX6_D2_R_fingers": "860-863",
}

# ---------------------------------------------------------------------------
# Colors (converted from your set_color RGB triplets -> hex)
# ---------------------------------------------------------------------------
COLOR_GRAY = "#BEBEBE"          # pymol "gray"
COLOR_FOREST = "#339933"        # pymol "forest" (PEX26)
COLOR_YELLOW = "#FFFF00"        # PEX5, Sensor_1, phospho show
COLOR_ORANGE = "#FF8000"        # phospho color
COLOR_CYAN = "#00FFFF"          # ISS motif
COLOR_DARK_PURPLE = "#660066"   # Walker A
COLOR_DARK_GREEN = "#008000"    # Walker B
COLOR_AMBER = "#E6B450"         # alpha insert
COLOR_WHITE = "#FFFFFF"         # R fingers carbons (util.cbaw approx)

DOMAIN_COLORS = {
    "PEX1_N1": "#989898",
    "PEX1_N2": "#DBC5AE",
    "PEX1_D1": "#FEBDAF",
    "PEX1_D2": "#DB798A",
    "PEX6_N1": "#B0D1BE",
    "PEX6_N2": "#97DFEB",
    "PEX6_D1": "#57A0CD",
    "PEX6_D2": "#8290F3",
}

MOTIF_COLORS = {
    "PEX1_D1_Walker_A": COLOR_DARK_PURPLE,
    "PEX1_D2_Walker_A": COLOR_DARK_PURPLE,
    "PEX6_D1_Walker_A": COLOR_DARK_PURPLE,
    "PEX6_D2_Walker_A": COLOR_DARK_PURPLE,
    "PEX1_D1_Walker_B": COLOR_DARK_GREEN,
    "PEX1_D2_Walker_B": COLOR_DARK_GREEN,
    "PEX6_D1_Walker_B": COLOR_DARK_GREEN,
    "PEX6_D2_Walker_B": COLOR_DARK_GREEN,
    "PEX1_D1_Sensor_1": COLOR_YELLOW,
    "PEX1_D2_Sensor_1": COLOR_YELLOW,
    "PEX6_D1_Sensor_1": COLOR_YELLOW,
    "PEX6_D2_Sensor_1": COLOR_YELLOW,
    "PEX1_D1_alpha_insert": COLOR_AMBER,
    "PEX1_D2_alpha_insert": COLOR_AMBER,
    "PEX6_D1_alpha_insert": COLOR_AMBER,
    "PEX6_D2_alpha_insert": COLOR_AMBER,
    "PEX1_D2_ISS": COLOR_CYAN,
    "PEX6_D2_ISS": COLOR_CYAN,
}

STICK_MOTIFS = [
    "PEX1_D2_DE", "PEX6_D2_DE",
    "PEX1_D2_Pore_loop_1", "PEX6_D2_Pore_loop_1",
    "PEX1_D2_R_fingers", "PEX6_D2_R_fingers",
]


def chain_range_selector(chains, resi_range, offset=0):
    """
    Build a list of ComponentExpressions for one residue range across a chain group.
    `offset` = FIRST_RESI[protein] - 1, shifting canonical ranges to actual structure numbering.
    """
    lo, hi = (int(x) for x in resi_range.split("-"))
    lo, hi = lo + offset, hi + offset
    return [
        mvs.ComponentExpression(auth_asym_id=c, beg_auth_seq_id=lo, end_auth_seq_id=hi)
        for c in chains
    ]


def chain_selector(chains):
    return [mvs.ComponentExpression(auth_asym_id=c) for c in chains]


def resn_selector(resn_list):
    return [mvs.ComponentExpression(label_comp_id=r) for r in resn_list]


def build_scene(mode):
    """
    Build a fully independent MVS tree for one scene.
    mode: "cartoon" | "surface" | "plddt"

    NOTE: molviewspec's builder.get_snapshot() does NOT deep-copy the state tree
    (as of molviewspec 1.8.1 -- confirmed in source, has a `# TODO create deep
    copy of node` comment). Reusing one builder across snapshots means later
    mutations silently leak into earlier snapshots. To avoid that, each scene
    gets its own from-scratch builder here.
    """
    builder = mvs.create_builder()
    structure = (
        builder.download(url=STRUCTURE_URL)
        .parse(format=STRUCTURE_FORMAT)
        .model_structure()
    )

    polymer = structure.component(selector="polymer")
    cartoon_repr = polymer.representation(type="cartoon", ref="cartoon")

    if mode == "plddt":
        cartoon_repr.color_from_source(
            schema="atom",
            category_name="atom_site",
            field_name="B_iso_or_equiv",
            palette={
                "kind": "continuous",
                "colors": ["#E68C3A", "#FFFF00", "#5C88DA", "#1E3A8A"],  # orange -> yellow -> marine -> blue
                "mode": "normalized",
                "value_domain": [20, 100],
                "underflow_color": "auto",  # clip values < 20 to the lowest checkpoint color (matches PyMOL `spectrum` clipping)
                "overflow_color": "auto",   # clip values > 100 to the highest checkpoint color
            },
        )
    elif mode == "plddt_bins":
        # Official AlphaFold DB / EBI discrete 4-bin scheme (hard thresholds, not a gradient):
        # >90 dark blue (very high), 70-90 light blue (confident), 50-70 yellow (low), <50 orange (very low)
        cartoon_repr.color_from_source(
            schema="atom",
            category_name="atom_site",
            field_name="B_iso_or_equiv",
            palette={
                "kind": "discrete",
                "mode": "absolute",
                "colors": [
                    ["#FF8000", None, 50.0],
                    ["#FFFF00", 50.0, 70.0],
                    ["#87CEFA", 70.0, 90.0],
                    ["#00008B", 90.0, None],
                ],
            },
        )
    else:
        cartoon_repr.color(color=COLOR_GRAY)

        all_ranges = {
            **{f"PEX1::{k}": (PEX1_CHAINS, v, FIRST_RESI["PEX1"] - 1) for k, v in PEX1_RANGES.items()},
            **{f"PEX6::{k}": (PEX6_CHAINS, v, FIRST_RESI["PEX6"] - 1) for k, v in PEX6_RANGES.items()},
        }
        for key, (chains, resi_range, offset) in all_ranges.items():
            domain_name = key.split("::")[1]
            color = DOMAIN_COLORS.get(domain_name) or MOTIF_COLORS.get(domain_name)
            if color is None:
                continue
            for sel in chain_range_selector(chains, resi_range, offset):
                cartoon_repr.color(selector=sel, color=color)

        for sel in chain_selector(PEX26_CHAINS):
            cartoon_repr.color(selector=sel, color=COLOR_FOREST)

        phospho_component = structure.component(selector=resn_selector(PHOSPHO_RESN))
        phospho_component.representation(type="ball_and_stick").color(color=COLOR_ORANGE)

        ligand_component = structure.component(selector=resn_selector(LIGAND_RESN))
        ligand_component.representation(type="ball_and_stick")

        for motif in STICK_MOTIFS:
            protein = "PEX1" if motif.startswith("PEX1") else "PEX6"
            chains = PEX1_CHAINS if protein == "PEX1" else PEX6_CHAINS
            resi_range = PEX1_RANGES.get(motif) or PEX6_RANGES.get(motif)
            offset = FIRST_RESI[protein] - 1
            if "R_fingers" in motif:
                color = COLOR_WHITE
            elif "Pore_loop" in motif:
                color = COLOR_ORANGE
            else:
                # DE motif: PyMOL never gives it its own `color` command -- it just
                # inherits whatever the D2 domain color already was. Match that here
                # instead of a hardcoded color.
                color = DOMAIN_COLORS[f"{protein}_D2"]
            sticks = structure.component(selector=chain_range_selector(chains, resi_range, offset))
            sticks.representation(type="ball_and_stick").color(color=color)

        for sel in chain_selector(PEX5_CHAINS):
            cartoon_repr.color(selector=sel, color=COLOR_YELLOW)

    pex5_component = structure.component(selector=chain_selector(PEX5_CHAINS))

    if mode == "surface":
        pex5_surface = pex5_component.representation(type="surface", ref="pex5_surface")
        pex5_surface.color(color=COLOR_YELLOW)
        pex5_surface.opacity(opacity=0.7)  # PyMOL transparency=0.3 -> opacity=0.7

    titles = {
        "cartoon": ("Cartoon view", "Domain-colored cartoon (PEX1/PEX6 domains, Walker motifs, ISS, PEX26, PEX5)."),
        "surface": ("Surface view", "Same as cartoon view, with PEX5 (if applicable) shown as a semi-transparent surface."),
        "plddt": ("pLDDT view", "Colored by predicted confidence (B-factor), continuous gradient matching "
                                 "`spectrum b, tv_orange tv_yellow marine blue, minimum=20, maximum=100`."),
        "plddt_bins": ("pLDDT view (AF DB bins)", "Colored by predicted confidence (B-factor) using the official "
                                                    "AlphaFold DB discrete scheme: dark blue >90 (very high), "
                                                    "light blue 70-90 (confident), yellow 50-70 (low), orange <50 (very low)."),
    }
    title, description = titles[mode]
    return builder.get_snapshot(title=title, description=description)


snapshot_cartoon = build_scene("cartoon")
snapshot_surface = build_scene("surface")
snapshot_plddt = build_scene("plddt")
snapshot_plddt_bins = build_scene("plddt_bins")

# ---------------------------------------------------------------------------
# Combine into a navigable multi-snapshot state
# ---------------------------------------------------------------------------
states = mvs.States(
    snapshots=[snapshot_cartoon, snapshot_surface, snapshot_plddt, snapshot_plddt_bins],
    metadata=mvs.GlobalMetadata(description="PEX1-PEX6-PEX26-PEX5 complex — scene navigation"),
)

with open(OUTPUT_PATH, "w") as f:
    f.write(states.model_dump_json(exclude_none=True, indent=2))

print(f"Wrote {OUTPUT_PATH}")
