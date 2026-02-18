import os
from plasmid_maker import (
    read_fasta,
    parse_design,
    load_markers,
    build_plasmid,
    PlasmidDesign,
    RE_SITE,
)


def test_puc19_ecoRI_removed():
    
    input_fa = "pUC19.fa"
    design_txt = "Design_pUC19.txt"
    markers_tab = "markers.tab"

    assert os.path.exists(input_fa), f"Missing test file: {input_fa}"
    assert os.path.exists(design_txt), f"Missing test file: {design_txt}"
    assert os.path.exists(markers_tab), f"Missing test file: {markers_tab}"

    _, original_seq = read_fasta(input_fa)
    eco_motif = RE_SITE["EcoRI"]

    # Ensure EcoRI site exists in the original pUC19 input
    assert eco_motif in original_seq, "Original pUC19.fa must contain EcoRI site GAATTC"

    design: PlasmidDesign = parse_design(design_txt)
    markers = load_markers(markers_tab)

    # For test, always force EcoRI removal, independent of infer_sites_to_remove().
    out_seq = build_plasmid(
        original_seq,
        design=design,
        markers=markers,
        remove_sites=["EcoRI"],
    )

    assert eco_motif not in out_seq, "Output plasmid still contains EcoRI site GAATTC"
