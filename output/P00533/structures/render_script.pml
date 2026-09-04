# ProteinScope Automated PyMOL Visualization Script
reinitialize
load C:/Users/Sree Ganesh/Dropbox/computational_protein_pipeline/output/P00533/structures/experimental/8a27.pdb, protein
bg_color white
hide everything, all

# Color Scheme: Secondary Structure
show cartoon, protein
color cyan, ss h        # Alpha helices
color yellow, ss s      # Beta sheets
color white, ss l+''    # Loops / Coils
set cartoon_fancy_helices, 1
set ray_shadows, 0
zoom protein
orient protein
png C:/Users/Sree Ganesh/Dropbox/computational_protein_pipeline/output/P00533/figures/structure_full.png, width=1600, height=1200, dpi=300, ray=1

# Surface View
show surface, protein
set transparency, 0.2
png C:/Users/Sree Ganesh/Dropbox/computational_protein_pipeline/output/P00533/figures/structure_surface.png, width=1600, height=1200, dpi=300, ray=1
hide surface, protein

# Domain View
color spectrum, protein
png C:/Users/Sree Ganesh/Dropbox/computational_protein_pipeline/output/P00533/figures/structure_domains.png, width=1600, height=1200, dpi=300, ray=1

save C:/Users/Sree Ganesh/Dropbox/computational_protein_pipeline/output/P00533/structures/session.pse
quit