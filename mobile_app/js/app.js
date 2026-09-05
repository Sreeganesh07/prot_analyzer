/**
 * ProteinScope Mobile Application Controller
 * Manages 3Dmol WebGL viewer, offline database, client-side biophysics, and online search.
 */


let currentProtein = null;
let viewer3D = null;
let renderMode = 'cartoon';
let isDegradedVisual = false;
let isSpinning = false;
let blinkInterval = null;

// ── AMINO ACID ENCYCLOPEDIA (PROPERTIES, CHARGES, PKA, & BIOLOGICAL ROLES) ──
const AA_DICT = {
    'A': { name: 'Alanine', code3: 'ALA', category: 'Aliphatic Hydrophobic', charge: 0.0, pKa: 'None', hydropathy: 1.80, role: 'Stabilizes hydrophobic core; high helix-forming propensity.' },
    'R': { name: 'Arginine', code3: 'ARG', category: 'Basic Positively Charged', charge: +1.0, pKa: '12.48 (Guanidinium)', hydropathy: -4.50, role: 'Forms bidentate salt bridges with acidic residues and clamps DNA phosphate backbones.' },
    'N': { name: 'Asparagine', code3: 'ASN', category: 'Polar Uncharged', charge: 0.0, pKa: 'None', hydropathy: -3.50, role: 'Stabilizes turns via amide H-bonds; consensus N-glycosylation acceptor (Asn-X-Ser/Thr).' },
    'D': { name: 'Aspartate', code3: 'ASP', category: 'Acidic Negatively Charged', charge: -1.0, pKa: '3.90 (β-Carboxyl)', hydropathy: -3.50, role: 'Electrostatic salt bridges; coordinates catalytic Mg²⁺/Ca²⁺ ions; acid-base catalysis.' },
    'C': { name: 'Cysteine', code3: 'CYS', category: 'Special / Thiol Reactive', charge: 0.0, pKa: '8.33 (Thiol -SH)', hydropathy: 2.50, role: 'Forms covalent disulfide bonds (-S-S-) locking tertiary fold; coordinates catalytic Zn²⁺ ions.' },
    'E': { name: 'Glutamate', code3: 'GLU', category: 'Acidic Negatively Charged', charge: -1.0, pKa: '4.07 (γ-Carboxyl)', hydropathy: -3.50, role: 'Surface hydration shell stabilizer; acts as general acid/base nucleophile in active sites.' },
    'Q': { name: 'Glutamine', code3: 'GLN', category: 'Polar Uncharged', charge: 0.0, pKa: 'None', hydropathy: -3.50, role: 'Flexible polar side chain promoting tertiary H-bonding networks and metabolic nitrogen shuttle.' },
    'G': { name: 'Glycine', code3: 'GLY', category: 'Special / Achiral & Flexible', charge: 0.0, pKa: 'None', hydropathy: -0.40, role: 'Lacks side chain; confers conformational freedom to tight beta-hairpin turns and hinges.' },
    'H': { name: 'Histidine', code3: 'HIS', category: 'Basic / Imidazole Buffer', charge: +0.10, pKa: '6.00 (Imidazole)', hydropathy: -3.20, role: 'Near-neutral physiological pKa allows rapid reversible proton shuttling; key biological fluid buffer.' },
    'I': { name: 'Isoleucine', code3: 'ILE', category: 'Aliphatic Hydrophobic', charge: 0.0, pKa: 'None', hydropathy: 4.50, role: 'Bulky β-branched hydrophobic core stabilizer; strongly drives hydrophobic collapse.' },
    'L': { name: 'Leucine', code3: 'LEU', category: 'Aliphatic Hydrophobic', charge: 0.0, pKa: 'None', hydropathy: 3.80, role: 'Major hydrophobic packaging residue; forms leucine zipper coiled-coil dimerization motifs.' },
    'K': { name: 'Lysine', code3: 'LYS', category: 'Basic Positively Charged', charge: +1.0, pKa: '10.53 (ε-Amino)', hydropathy: -3.90, role: 'Solvation surface salt bridges; primary substrate for ubiquitination, acetylation, and SUMOylation.' },
    'M': { name: 'Methionine', code3: 'MET', category: 'Hydrophobic / Thioether', charge: 0.0, pKa: 'None', hydropathy: 1.90, role: 'Universal translation initiator; hydrophobic packing; acts as reversible reactive oxygen sensor.' },
    'F': { name: 'Phenylalanine', code3: 'PHE', category: 'Aromatic Hydrophobic', charge: 0.0, pKa: 'None', hydropathy: 2.80, role: 'Nonpolar interior core packaging; stabilizes tertiary fold through aromatic π-π stacking.' },
    'P': { name: 'Proline', code3: 'PRO', category: 'Special / Cyclic Imino Acid', charge: 0.0, pKa: 'None', hydropathy: -1.60, role: 'Rigid pyrrolidine ring introduces strict conformational kinks; canonical α-helix breaker and turn inducer.' },
    'S': { name: 'Serine', code3: 'SER', category: 'Polar Hydroxyl', charge: 0.0, pKa: '13.00 (Hydroxyl)', hydropathy: -0.80, role: 'Forms active site catalytic triads (e.g. serine proteases); major target for regulatory phosphorylation.' },
    'T': { name: 'Threonine', code3: 'THR', category: 'Polar Hydroxyl / β-branched', charge: 0.0, pKa: '13.00 (Hydroxyl)', hydropathy: -0.70, role: 'Secondary alcohol side chain; regulatory phosphorylation and O-linked glycosylation site.' },
    'W': { name: 'Tryptophan', code3: 'TRP', category: 'Aromatic Indole', charge: 0.0, pKa: 'None', hydropathy: -0.90, role: 'Bulkiest sidechain; dominant contributor to intrinsic protein UV absorbance (280 nm) & fluorescence.' },
    'Y': { name: 'Tyrosine', code3: 'TYR', category: 'Aromatic Hydroxyl', charge: 0.0, pKa: '10.07 (Phenolic OH)', hydropathy: -1.30, role: 'Amphipathic aromatic; target for receptor tyrosine kinase (RTK) phosphorylation and H-bonding.' },
    'V': { name: 'Valine', code3: 'VAL', category: 'Aliphatic Hydrophobic', charge: 0.0, pKa: 'None', hydropathy: 4.20, role: 'Sterically rigid β-branched nonpolar side chain favoring β-sheet conformations.' }
};

// ── APP INITIALIZATION ──
document.addEventListener('DOMContentLoaded', () => {
    initNetworkStatus();
    initCanisterTrack();
    setupEventListeners();
    setupSearchModal();

    // Load initial default protein (p53)
    const db = window.PROTEINSCOPE_CACHED_DB || {};
    const firstAcc = Object.keys(db)[0] || 'P04637';
    loadProtein(firstAcc);

    // Dismiss splash screen with Instagram/YouTube style smooth transition
    setTimeout(() => {
        const splash = document.getElementById('splash-screen');
        if (splash) {
            splash.classList.add('splash-hidden');
            setTimeout(() => {
                splash.style.display = 'none';
            }, 600);
        }
    }, 1300);
});

// ── NETWORK STATUS ──
function initNetworkStatus() {
    function update() {
        const badge = document.getElementById('net-status-badge');
        const text = document.getElementById('net-status-text');
        if (!badge || !text) return;

        if (navigator.onLine) {
            badge.classList.remove('offline-mode');
            text.textContent = 'ONLINE READY';
        } else {
            badge.classList.add('offline-mode');
            text.textContent = 'OFFLINE (9 CACHED)';
        }
    }

    window.addEventListener('online', update);
    window.addEventListener('offline', update);
    update();
}

// ── SPECIMEN CAROUSEL ──
function initCanisterTrack() {
    const track = document.getElementById('canister-track');
    if (!track) return;
    track.innerHTML = '';

    const db = window.PROTEINSCOPE_CACHED_DB || {};
    const keys = Object.keys(db);

    const countEl = document.getElementById('canister-count');
    if (countEl) countEl.textContent = `${keys.length} Cached`;

    keys.forEach(acc => {
        const p = db[acc];
        const card = document.createElement('div');
        card.className = 'can-card';
        card.id = `canister-${acc}`;
        card.innerHTML = `
            <span class="can-acc">${p.accession}</span>
            <span class="can-name">${(p.name || acc).split('(')[0].trim()}</span>
            <span class="can-tm">${(p.estimated_tm || 65.0).toFixed(1)}°C</span>
        `;
        card.onclick = () => loadProtein(acc);
        track.appendChild(card);
    });
}

// ── TAB NAVIGATION ──
window.switchTab = function(tabName) {
    const tabs = ['structure', 'biophysics', 'properties'];
    tabs.forEach(t => {
        const el = document.getElementById(`tab-${t}`);
        const nav = document.getElementById(`nav-${t}`);
        if (el) el.classList.remove('active-tab');
        if (nav) nav.classList.remove('nav-active');
    });

    const activeEl = document.getElementById(`tab-${tabName}`);
    const activeNav = document.getElementById(`nav-${tabName}`);
    if (activeEl) activeEl.classList.add('active-tab');
    if (activeNav) activeNav.classList.add('nav-active');

    // If returning to 3D view, ensure WebGL canvas is properly sized and rendered
    if (tabName === 'structure' && viewer3D) {
        setTimeout(() => {
            viewer3D.resize();
            viewer3D.render();
        }, 50);
    }
};

// ── EVENT LISTENERS ──
function setupEventListeners() {
    const sliderTemp = document.getElementById('slider-temp');
    const sliderPH = document.getElementById('slider-ph');

    if (sliderTemp) {
        sliderTemp.addEventListener('input', (e) => {
            const val = parseFloat(e.target.value);
            const valEl = document.getElementById('val-temp');
            if (valEl) {
                valEl.innerHTML = `${val.toFixed(1)}°C <span style="font-size:10px; color:var(--text-muted);">(${(val + 273.15).toFixed(2)} K)</span>`;
            }
            recalculateStability();
        });
    }

    if (sliderPH) {
        sliderPH.addEventListener('input', (e) => {
            const val = parseFloat(e.target.value);
            const valEl = document.getElementById('val-ph');
            if (valEl) {
                valEl.textContent = `pH ${val.toFixed(2)}`;
            }
            recalculateStability();
        });
    }

    // Window resize handler for 3Dmol canvas
    window.addEventListener('resize', () => {
        if (viewer3D) {
            viewer3D.resize();
            viewer3D.render();
        }
    });
}

// ── GLOBAL SLIDER PRESET SHORTCUTS ──
window.setTemp = function(temp) {
    const slider = document.getElementById('slider-temp');
    if (slider) {
        slider.value = temp;
        const valEl = document.getElementById('val-temp');
        if (valEl) {
            valEl.innerHTML = `${temp.toFixed(1)}°C <span style="font-size:10px; color:var(--text-muted);">(${(temp + 273.15).toFixed(2)} K)</span>`;
        }
        recalculateStability();
    }
};

window.setPH = function(ph) {
    const slider = document.getElementById('slider-ph');
    if (slider) {
        slider.value = ph;
        const valEl = document.getElementById('val-ph');
        if (valEl) {
            valEl.textContent = `pH ${ph.toFixed(2)}`;
        }
        recalculateStability();
    }
};

// ── 3D MOLECULAR VIEWER CONTROLS ──
window.setRenderMode = function(mode) {
    renderMode = mode;
    document.querySelectorAll('.style-btn').forEach(btn => {
        if (btn.textContent.toLowerCase().trim() === mode.toLowerCase().trim()) {
            btn.classList.add('active');
        } else {
            btn.classList.remove('active');
        }
    });
    apply3DStyle();
};

window.toggleSpin = function() {
    if (!viewer3D) return;
    isSpinning = !isSpinning;
    const spinBtn = document.getElementById('spin-btn');
    if (spinBtn) {
        spinBtn.classList.toggle('active-toggle', isSpinning);
    }
    if (isSpinning) {
        viewer3D.spin('y', 0.8);
    } else {
        viewer3D.spin(false);
    }
};

window.reset3DCamera = function() {
    if (viewer3D) {
        apply3DStyle();
        viewer3D.zoomTo();
        viewer3D.render();
    }
};

window.highlightDomainIn3D = function(domainIdx) {
    if (!currentProtein || !currentProtein.domains || !currentProtein.domains[domainIdx]) return;
    const dom = currentProtein.domains[domainIdx];
    
    // Switch to 3D viewer tab so user immediately sees the visual focus
    switchTab('structure');

    const resiRange = Array.from({ length: dom.end - dom.start + 1 }, (_, i) => dom.start + i);
    const infoChip = document.getElementById('viewer-info-chip');

    // Scroll sequence strip smoothly to start of domain
    const firstPill = document.getElementById(`seq-pill-${dom.start}`);
    if (firstPill) {
        firstPill.scrollIntoView({ behavior: 'smooth', inline: 'center', block: 'nearest' });
    }

    if (viewer3D) {
        // Remove any surface meshes and reset
        viewer3D.removeAllSurfaces();
        viewer3D.removeAllShapes();

        // Safe atom check to prevent camera NaN when domain residues are unresolved in crystal
        const selectedAtoms = viewer3D.selectedAtoms({ resi: resiRange });
        if (selectedAtoms && selectedAtoms.length > 0) {
            // Mute non-domain residues into translucent ghost slate
            viewer3D.setStyle({}, {
                cartoon: { color: '#334155', opacity: 0.28, thickness: 0.40 }
            });

            // Highlight selected domain with bold glowing color
            viewer3D.setStyle(
                { resi: resiRange },
                { cartoon: { color: dom.color, thickness: 0.85, opacity: 1.0 } }
            );

            // Safe zoom without second duration parameter (prevents camera matrix corruption)
            viewer3D.zoomTo({ resi: resiRange });
            viewer3D.render();

            if (infoChip) {
                infoChip.textContent = `FOCUS // ${dom.name.toUpperCase()} (${dom.start}-${dom.end})`;
            }
        } else {
            // Residues outside the experimental crystal boundary: keep whole model safely in view
            viewer3D.zoomTo();
            viewer3D.render();
            if (infoChip) {
                infoChip.textContent = `DOMAIN // ${dom.name.toUpperCase()} (${dom.start}-${dom.end}) [UNRESOLVED IN CRYSTAL]`;
            }
        }
    }

    // Trigger synchronized 3-second blinking pulse across 3D canvas and sequence strip
    startCoordinatedBlink(resiRange);
};

window.toggleStructureComparison = async function(mode) {
    if (!currentProtein) return;
    const infoChip = document.getElementById('viewer-info-chip');
    const comp = currentProtein.structural_comparison;

    if (mode === 'experimental') {
        const expId = comp?.experimental_id || (currentProtein.pdb_cross_references && currentProtein.pdb_cross_references[0]);
        if (expId && expId !== "None Solved") {
            await loadExperimentalPDB(expId);
        } else {
            alert("No solved experimental co-crystal registered for this sequence. Rendering comparative predicted model.");
        }
    } else if (mode === 'predicted') {
        if (currentProtein.pdb_content) {
            init3DViewer(currentProtein.pdb_content);
            switchTab('structure');
            if (infoChip) {
                infoChip.textContent = `PREDICTED MODEL // ${comp?.predicted_model || 'AlphaFold AI'}`;
            }
        }
    }
};

// ── LOAD PROTEIN ──
function loadProtein(acc) {
    const db = window.PROTEINSCOPE_CACHED_DB || {};
    const p = db[acc];
    if (!p) return;

    currentProtein = p;

    // Update active canister tab in carousel
    document.querySelectorAll('.can-card').forEach(b => b.classList.remove('active'));
    const activeBtn = document.getElementById(`canister-${acc}`);
    if (activeBtn) {
        activeBtn.classList.add('active');
        activeBtn.scrollIntoView({ behavior: 'smooth', inline: 'center', block: 'nearest' });
    }

    // Update Hero Banner
    document.getElementById('hero-acc').textContent = p.accession;
    document.getElementById('hero-name').textContent = p.name;
    document.getElementById('hero-org').textContent = p.organism;

    // Update Properties & Biology Tab Values
    const len = p.sequence_length || (p.sequence ? p.sequence.length : 300);
    const mw = p.molecular_weight_kda || (p.molecular_weight ? (p.molecular_weight / 1000).toFixed(1) : (len * 0.11).toFixed(1));
    const pi = p.isoelectric_point ? p.isoelectric_point.toFixed(2) : '6.50';
    const tm = (p.estimated_tm || 65.0).toFixed(1);

    document.getElementById('stat-len').textContent = `${len} aa`;
    document.getElementById('stat-mw').textContent = `${mw} kDa`;
    document.getElementById('stat-pi-val').textContent = pi;
    document.getElementById('disp-pi').textContent = `pI: ${pi}`;
    document.getElementById('stat-thresh-tm').textContent = `${tm}°C`;
    document.getElementById('stat-tm').textContent = `${tm}°C`;

    // Populate Rich Biological Fields
    const bioOrg = document.getElementById('bio-org-chip');
    if (bioOrg) bioOrg.textContent = p.organism || 'Biological Specimen';

    const bioGene = document.getElementById('bio-gene');
    if (bioGene) bioGene.textContent = `GENE: ${p.gene_name || p.accession}`;

    const bioSyn = document.getElementById('bio-synonyms');
    if (bioSyn) {
        const synList = p.gene_synonyms || [];
        bioSyn.textContent = synList.length > 0 ? `(${synList.slice(0, 3).join(', ')})` : '';
    }

    const subcellRow = document.getElementById('bio-subcell-row');
    if (subcellRow) {
        subcellRow.innerHTML = '';
        const locs = p.subcellular_location || ['Intracellular'];
        locs.forEach(loc => {
            const pill = document.createElement('span');
            pill.className = 'subcell-pill';
            pill.innerHTML = `📍 ${loc}`;
            subcellRow.appendChild(pill);
        });
    }

    const diseaseBanner = document.getElementById('bio-disease-banner');
    const diseaseText = document.getElementById('bio-disease-text');
    if (diseaseBanner && diseaseText) {
        if (p.disease_associations && p.disease_associations !== 'No direct clinical pathology registered.') {
            diseaseBanner.style.display = 'flex';
            diseaseText.textContent = p.disease_associations;
        } else {
            diseaseBanner.style.display = 'none';
        }
    }

    const functionText = document.getElementById('bio-function-text');
    if (functionText) {
        functionText.textContent = p.function_summary || 'Biological macromolecule investigated under physiological and thermal stress.';
    }

    // GRAVY & Extinction Coefficient
    const gravyEl = document.getElementById('stat-gravy');
    if (gravyEl) {
        gravyEl.textContent = p.gravy_score !== undefined ? (p.gravy_score > 0 ? '+' : '') + p.gravy_score.toFixed(3) : '-0.420';
    }

    const extEl = document.getElementById('stat-extcoeff');
    if (extEl) {
        extEl.textContent = p.extinction_coefficient ? p.extinction_coefficient.toLocaleString() : '35,870';
    }

    // Structure source chip
    const structType = document.getElementById('bio-struct-type');
    if (structType) {
        structType.textContent = p.pdb_content ? (p.accession.length === 4 ? 'RCSB PDB Experimental' : 'AlphaFold AI Structure') : 'Sequence Only';
    }

    // Populate Clickable Experimental PDB Cross-References
    const pdbRow = document.getElementById('bio-pdb-chips');
    if (pdbRow) {
        pdbRow.innerHTML = '';
        const xrefs = p.pdb_cross_references || [];
        if (xrefs.length > 0) {
            xrefs.forEach(pid => {
                const btn = document.createElement('button');
                btn.className = 'pdb-chip-btn';
                btn.innerHTML = `<span>🏛️</span><span>${pid}</span>`;
                btn.title = `Load experimental crystal ${pid}`;
                btn.onclick = () => loadExperimentalPDB(pid);
                pdbRow.appendChild(btn);
            });
        } else {
            pdbRow.innerHTML = '<span style="font-size:11px; color:var(--text-muted); font-style:italic;">No experimental PDB co-crystals cataloged for this entry.</span>';
        }
    }

    // ── POPULATE PROTEIN DOMAINS & ARCHITECTURE ──
    const domainsBar = document.getElementById('bio-domains-bar');
    const domainsList = document.getElementById('bio-domains-list');
    const domainsCountEl = document.getElementById('bio-domains-count');
    const doms = p.domains || [];

    if (domainsCountEl) {
        domainsCountEl.textContent = `${doms.length} Functional Modules`;
    }

    if (domainsBar) {
        domainsBar.innerHTML = '';
        if (doms.length > 0 && len > 0) {
            doms.forEach((dom, idx) => {
                const span = dom.end - dom.start + 1;
                const pct = Math.max(3, (span / len) * 100);
                const seg = document.createElement('div');
                seg.className = 'domain-bar-segment';
                seg.style.width = `${pct}%`;
                seg.style.backgroundColor = dom.color;
                seg.title = `${dom.name} (${dom.start}-${dom.end}): ${dom.purpose}`;
                seg.onclick = () => highlightDomainIn3D(idx);
                domainsBar.appendChild(seg);
            });
        } else {
            domainsBar.innerHTML = '<div style="width:100%; height:100%; background:var(--accent-primary); opacity:0.5; border-radius:4px;"></div>';
        }
    }

    if (domainsList) {
        domainsList.innerHTML = '';
        if (doms.length > 0) {
            doms.forEach((dom, idx) => {
                const card = document.createElement('div');
                card.className = 'domain-detail-card';
                card.innerHTML = `
                    <div class="domain-card-header">
                        <div class="domain-title-group">
                            <span class="domain-color-dot" style="background-color: ${dom.color}; box-shadow: 0 0 8px ${dom.color}88;"></span>
                            <span class="domain-card-name">${dom.name}</span>
                        </div>
                        <span class="domain-coords-badge">[${dom.start} – ${dom.end}]</span>
                    </div>
                    <div class="domain-purpose-text">${dom.purpose}</div>
                    <button class="domain-inspect-btn" onclick="highlightDomainIn3D(${idx})">
                        <span>🔍</span> Highlight in 3D
                    </button>
                `;
                domainsList.appendChild(card);
            });
        } else {
            domainsList.innerHTML = '<div style="font-size:11px; color:var(--text-muted); font-style:italic; padding:8px 0;">No segmented functional domains identified for this primary sequence.</div>';
        }
    }

    // ── POPULATE COMPARATIVE STRUCTURE PREDICTION & DELTA ANALYSIS ──
    const compBadge = document.getElementById('bio-comp-badge');
    const compModel = document.getElementById('bio-comp-model');
    const compRmsd = document.getElementById('bio-comp-rmsd');
    const compIdentity = document.getElementById('bio-comp-identity');
    const compNotes = document.getElementById('bio-comp-notes');
    const compLoops = document.getElementById('bio-comp-loops');
    const comp = p.structural_comparison;

    if (compBadge) {
        const isPred = p.is_predicted || (p.accession.length !== 4 && !comp?.has_experimental);
        if (isPred) {
            compBadge.className = 'comp-badge comp-badge-pred';
            compBadge.innerHTML = '<span>🔮</span> PREDICTED STRUCTURE';
        } else {
            compBadge.className = 'comp-badge comp-badge-exp';
            compBadge.innerHTML = '<span>🏛️</span> EXPERIMENTAL CO-CRYSTAL';
        }
    }

    if (compModel && comp) compModel.textContent = comp.predicted_model || 'AlphaFold AI v4';
    if (compRmsd && comp) compRmsd.textContent = `${comp.rmsd_angstroms.toFixed(2)} Å`;
    if (compIdentity && comp) compIdentity.textContent = `${comp.sequence_identity_percent.toFixed(1)}%`;
    if (compNotes && comp) compNotes.textContent = comp.conformational_deltas || 'Core tertiary fold exhibits high structural concordance.';
    if (compLoops && comp) compLoops.textContent = comp.flexible_loops || 'Terminal disordered tails';

    // External DB Links
    const linkUni = document.getElementById('link-uniprot');
    if (linkUni) linkUni.href = `https://www.uniprot.org/uniprotkb/${p.accession}`;

    const linkRcsb = document.getElementById('link-rcsb');
    if (linkRcsb) {
        const firstPdb = (p.pdb_cross_references && p.pdb_cross_references[0]) || (p.accession.length === 4 ? p.accession : '');
        linkRcsb.href = firstPdb ? `https://www.rcsb.org/structure/${firstPdb}` : `https://www.rcsb.org/search?request=%7B%22query%22%3A%7B%22type%22%3A%22terminal%22%2C%22service%22%3A%22text%22%2C%22parameters%22%3A%7B%22value%22%3A%22${p.accession}%22%7D%7D%7D`;
    }

    const linkAf = document.getElementById('link-alphafold');
    if (linkAf) linkAf.href = `https://alphafold.ebi.ac.uk/entry/${p.accession}`;

    // Render 3D molecular structure
    init3DViewer(p.pdb_content);

    // Recalculate biophysics stability
    recalculateStability();

    // Render Coordinated Sequence Strip & 3D Inspector
    renderSequenceStrip(p);

    // Render Comprehensive 5-Pillar Architecture Summary
    renderComprehensiveSummary(p);
}

window.loadExperimentalPDB = async function(pdbId) {
    const infoChip = document.getElementById('viewer-info-chip');
    if (infoChip) infoChip.textContent = `FETCHING RCSB PDB // ${pdbId}...`;

    try {
        const res = await fetch(`https://files.rcsb.org/download/${pdbId}.pdb`);
        if (!res.ok) throw new Error(`HTTP ${res.status}`);
        const pdbText = await res.text();
        init3DViewer(pdbText);
        switchTab('structure');
        if (infoChip) infoChip.textContent = `RCSB PDB // ${pdbId} CO-CRYSTAL ACTIVE`;
    } catch (e) {
        alert(`Could not download experimental crystal ${pdbId}: ${e.message}`);
    }
};

// ── INITIALIZE 3D VIEWER (OFFLINE LOCAL 3DMOL WITH MOBILE GPU OPTIMIZATION) ──
function init3DViewer(pdbText) {
    const container = document.getElementById('mol-canvas');
    const infoChip = document.getElementById('viewer-info-chip');
    if (!container) return;

    container.innerHTML = '';

    if (typeof $3Dmol === 'undefined') {
        container.innerHTML = `
            <div style="color:var(--text-secondary); padding:40px 20px; text-align:center; font-family:var(--font-mono); font-size:12px;">
                WebGL 3Dmol Library Initializing...
            </div>`;
        return;
    }

    try {
        // antialias: false eliminates severe GPU fill-rate bottlenecks on high-DPI mobile devices
        viewer3D = $3Dmol.createViewer(container, {
            backgroundColor: '#050811',
            defaultcolors: $3Dmol.rasmolElementColors,
            antialias: false,
            id: 'viewer3D_canvas'
        });

        // ── TOUCH ROTATION OPTIMIZATION ──
        // Pausing auto-spin during user touch prevents RAF spin loop and touch events from colliding
        let isTouching = false;
        container.addEventListener('touchstart', () => {
            isTouching = true;
            if (viewer3D && isSpinning) {
                viewer3D.spin(false);
            }
        }, { passive: true });

        container.addEventListener('touchend', () => {
            isTouching = false;
            if (viewer3D && isSpinning) {
                setTimeout(() => {
                    if (!isTouching && isSpinning && viewer3D) {
                        viewer3D.spin('y', 0.8);
                    }
                }, 400);
            }
        }, { passive: true });

        if (pdbText && pdbText.trim().length > 0) {
            viewer3D.addModel(pdbText, 'pdb');
            apply3DStyle();
            viewer3D.zoomTo();
            viewer3D.render();

            if (infoChip) {
                const badge = currentProtein.is_predicted ? 'PREDICTED' : 'ACTIVE';
                infoChip.textContent = `3DMOL // ${currentProtein.accession} ${badge}`;
            }
        } else {
            container.innerHTML = `
                <div style="display:flex; flex-direction:column; align-items:center; justify-content:center; height:100%; padding:24px; text-align:center;">
                    <div style="font-size:24px; margin-bottom:8px;">🧬</div>
                    <div style="color:#ffffff; font-weight:700; font-size:13px; margin-bottom:4px;">No 3D Coordinates Available</div>
                    <div style="color:var(--text-muted); font-size:11px; max-width:260px;">Biophysical stability and titration curves computed from primary sequence.</div>
                </div>`;
            if (infoChip) {
                infoChip.textContent = `BIOPHYSICS // SEQUENCE COMPUTED`;
            }
        }

        // Maintain auto-rotation state
        if (isSpinning) {
            viewer3D.spin('y', 0.8);
        }
    } catch (e) {
        console.error('Error in 3D viewer initialization:', e);
        container.innerHTML = `
            <div style="color:var(--accent-crimson); padding:30px; text-align:center; font-family:var(--font-mono); font-size:11px;">
                WebGL Initialization Error: ${e.message}
            </div>`;
    }
}

// ── APPLY 3D VISUAL STYLES (WITH CRITICAL SURFACE REMOVAL FIX) ──
function apply3DStyle() {
    if (!viewer3D) return;

    // CRITICAL: Must purge surface meshes and shapes so switching to cartoon/spheres removes them
    viewer3D.removeAllSurfaces();
    viewer3D.removeAllShapes();
    viewer3D.setStyle({}, {}); // Clear previous atom styles

    if (isDegradedVisual) {
        // Disorganized, aggregated appearance in red/crimson
        viewer3D.setStyle({}, {
            stick: { radius: 0.16, color: '#ef4444' },
            sphere: { radius: 0.42, color: '#991b1b' }
        });
    } else if (renderMode === 'surface') {
        viewer3D.setStyle({}, {
            cartoon: { color: 'spectrum', opacity: 0.35, thickness: 0.45 }
        });
        viewer3D.addSurface($3Dmol.SurfaceType.VDW, {
            opacity: 0.62,
            colorscheme: 'hydrophobicity'
        });
    } else if (renderMode === 'domains') {
        // Color-code distinct functional regions and domains
        const doms = currentProtein?.domains || [];
        if (doms.length > 0) {
            viewer3D.setStyle({}, { cartoon: { color: '#475569', opacity: 0.30, thickness: 0.40 } });
            doms.forEach(d => {
                const resiArr = Array.from({ length: d.end - d.start + 1 }, (_, i) => d.start + i);
                viewer3D.setStyle(
                    { resi: resiArr },
                    { cartoon: { color: d.color, thickness: 0.65 } }
                );
            });
        } else {
            viewer3D.setStyle({}, { cartoon: { color: 'spectrum', thickness: 0.55 } });
        }
    } else if (renderMode === 'spheres') {
        viewer3D.setStyle({}, {
            sphere: { colorscheme: 'amino', scale: 0.28 }
        });
    } else {
        // Default: High-contrast Secondary Structure Cartoon Ribbon
        viewer3D.setStyle({ ss: 'h' }, { cartoon: { color: '#0ea5e9', thickness: 0.55 } }); // Helices (Cyan)
        viewer3D.setStyle({ ss: 's' }, { cartoon: { color: '#f59e0b', thickness: 0.55 } }); // Sheets (Gold)
        viewer3D.setStyle({ ss: 'c' }, { cartoon: { color: '#64748b', thickness: 0.35 } }); // Coils (Slate)
    }

    viewer3D.render();
}

// ── BIOPHYSICS & THERMODYNAMIC STABILITY RECALCULATION ──
function recalculateStability() {
    if (!currentProtein) return;

    const temp = parseFloat(document.getElementById('slider-temp').value) || 37.0;
    const ph = parseFloat(document.getElementById('slider-ph').value) || 7.4;

    const evalResult = BiophysicsEngine.evaluateState(
        temp,
        ph,
        currentProtein.estimated_tm || 65.0,
        currentProtein.titratable_counts,
        currentProtein.thresholds
    );

    // Update Quick Telemetry Under 3D Viewer
    const foldedEl = document.getElementById('stat-folded');
    foldedEl.textContent = `${evalResult.fraction_folded_percent.toFixed(1)}%`;
    foldedEl.style.color = evalResult.fraction_folded_percent > 80 ? 'var(--accent-emerald)' : 
                          (evalResult.fraction_folded_percent > 40 ? 'var(--accent-amber)' : 'var(--accent-crimson)');

    const chargeEl = document.getElementById('stat-charge');
    chargeEl.textContent = `${evalResult.net_charge > 0 ? '+' : ''}${evalResult.net_charge.toFixed(2)} e`;

    document.getElementById('stat-dg').textContent = `${evalResult.delta_g_folding_kcal_mol.toFixed(1)} kcal`;
    document.getElementById('stat-salt').textContent = `${evalResult.salt_bridge_retention_percent.toFixed(1)}%`;

    // Update Hero State Pill
    const statePill = document.getElementById('hero-state-pill');
    const stateText = document.getElementById('hero-state-text');
    statePill.style.color = evalResult.badge_color;
    stateText.textContent = evalResult.badge_label;

    // Update Critical Degradation Alert
    const banner = document.getElementById('deg-banner');
    const viewerCard = document.getElementById('viewer-card');
    if (evalResult.is_degraded) {
        banner.style.display = 'block';
        document.getElementById('alarm-title').textContent = evalResult.badge_label;
        document.getElementById('alarm-desc').textContent = evalResult.reasons.join(' • ') + ' — ' + evalResult.mechanism;
        viewerCard.classList.add('degraded-glow');

        if (!isDegradedVisual) {
            isDegradedVisual = true;
            apply3DStyle();
        }
    } else {
        banner.style.display = 'none';
        viewerCard.classList.remove('degraded-glow');

        if (isDegradedVisual) {
            isDegradedVisual = false;
            apply3DStyle();
        }
    }

    // Render Henderson-Hasselbalch Titration Curve
    renderTitrationCurve(ph, evalResult.net_charge);
}

// ── SVG CHARGE TITRATION CURVE ──
function renderTitrationCurve(currentPH, currentQ) {
    if (!currentProtein) return;
    const path = document.getElementById('svg-curve');
    const dot = document.getElementById('svg-dot');
    if (!path || !dot) return;

    let d = '';
    const maxQ = Math.max(20, Math.abs(currentQ) * 1.5, 40);

    for (let ph = 0; ph <= 14; ph += 0.5) {
        const q = BiophysicsEngine.calculateNetCharge(currentProtein.titratable_counts, ph);
        const x = (ph / 14.0) * 320.0;
        const y = 50.0 - (q / maxQ) * 45.0;
        d += (ph === 0 ? 'M ' : ' L ') + `${x.toFixed(1)},${Math.max(2, Math.min(98, y)).toFixed(1)}`;
    }
    path.setAttribute('d', d);

    // Position dot
    const dotX = (currentPH / 14.0) * 320.0;
    const dotY = 50.0 - (currentQ / maxQ) * 45.0;
    dot.setAttribute('cx', dotX.toFixed(1));
    dot.setAttribute('cy', Math.max(4, Math.min(96, dotY)).toFixed(1));
}

// ── SEARCH MODAL & RETRIEVAL WORKFLOW ──
function setupSearchModal() {
    const modal = document.getElementById('search-modal');
    const headerBtn = document.getElementById('header-search-btn');
    const submitBtn = document.getElementById('search-submit-btn');
    const input = document.getElementById('search-input');

    if (headerBtn) {
        headerBtn.onclick = openSearchModal;
    }

    if (modal) {
        modal.onclick = (e) => {
            if (e.target === modal) closeSearchModal();
        };
    }

    if (submitBtn) {
        submitBtn.onclick = handleSearch;
    }

    if (input) {
        input.onkeypress = (e) => {
            if (e.key === 'Enter') handleSearch();
        };
    }
}

window.openSearchModal = function() {
    const modal = document.getElementById('search-modal');
    const feedback = document.getElementById('search-feedback');
    if (feedback) feedback.style.display = 'none';

    if (modal) {
        modal.classList.add('modal-open');
        setTimeout(() => {
            const input = document.getElementById('search-input');
            if (input) input.focus();
        }, 150);
    }
};

window.closeSearchModal = function() {
    const modal = document.getElementById('search-modal');
    if (modal) modal.classList.remove('modal-open');
};

window.quickSearch = function(query) {
    const input = document.getElementById('search-input');
    if (input) input.value = query;
    handleSearch();
};

async function handleSearch() {
    const input = document.getElementById('search-input');
    const feedback = document.getElementById('search-feedback');
    const query = input.value.trim();

    if (!query) {
        feedback.style.display = 'block';
        feedback.style.background = 'rgba(245, 158, 11, 0.15)';
        feedback.style.color = 'var(--accent-amber)';
        feedback.textContent = 'Please enter a PDB ID, UniProt accession, or protein name.';
        return;
    }

    feedback.style.display = 'block';
    feedback.style.background = 'rgba(56, 189, 248, 0.12)';
    feedback.style.color = 'var(--accent-cyan)';
    feedback.textContent = `Connecting & retrieving '${query}' (PDB / UniProt / AlphaFold)...`;

    try {
        const fetchedData = await BiophysicsEngine.fetchOnlineProtein(query);
        feedback.style.background = 'rgba(16, 185, 129, 0.15)';
        const has3D = fetchedData.pdb_content && fetchedData.pdb_content.length > 0;
        const structNote = has3D ? ` • 3D Structure Ready (${(fetchedData.pdb_content.length / 1024).toFixed(0)} KB)` : '';
        feedback.textContent = `✓ Retrieved: ${fetchedData.name} (${fetchedData.sequence_length} aa)${structNote}`;

        // Cache in memory database
        window.PROTEINSCOPE_CACHED_DB[fetchedData.accession] = fetchedData;

        // Add to canister carousel if not already present
        let card = document.getElementById(`canister-${fetchedData.accession}`);
        if (!card) {
            const track = document.getElementById('canister-track');
            card = document.createElement('div');
            card.className = 'can-card';
            card.id = `canister-${fetchedData.accession}`;
            card.innerHTML = `
                <span class="can-acc">${fetchedData.accession}</span>
                <span class="can-name">${(fetchedData.name || fetchedData.accession).split('(')[0].trim()}</span>
                <span class="can-tm">${(fetchedData.estimated_tm || 65.0).toFixed(1)}°C</span>
            `;
            card.onclick = () => loadProtein(fetchedData.accession);
            track.appendChild(card);

            const countEl = document.getElementById('canister-count');
            if (countEl) countEl.textContent = `${Object.keys(window.PROTEINSCOPE_CACHED_DB).length} Cached`;
        }

        // Load new protein immediately and switch to 3D structure tab
        loadProtein(fetchedData.accession);
        switchTab('structure');

        setTimeout(() => {
            closeSearchModal();
        }, 1200);
    } catch (err) {
        feedback.style.background = 'rgba(239, 68, 68, 0.15)';
        feedback.style.color = 'var(--accent-crimson)';
        feedback.textContent = `Search Error: ${err.message}`;
    }
}

// ── COORDINATED 3-SECOND BLINKING ENGINE (SEQUENCE STRIP <-> 3D VIEWPORT) ──
function startCoordinatedBlink(targetResidues) {
    if (blinkInterval) {
        clearInterval(blinkInterval);
        blinkInterval = null;
    }

    // Clear previous blinking animation classes
    document.querySelectorAll('.blinking-pill').forEach(el => el.classList.remove('blinking-pill'));

    if (!targetResidues || targetResidues.length === 0) return;

    // Apply pulsing CSS keyframe animation to corresponding sequence pills
    targetResidues.forEach(resNum => {
        const pill = document.getElementById(`seq-pill-${resNum}`);
        if (pill) pill.classList.add('blinking-pill');
    });

    if (!viewer3D) return;

    // Verify atoms exist in 3D canvas
    const selectedAtoms = viewer3D.selectedAtoms({ resi: targetResidues });
    if (!selectedAtoms || selectedAtoms.length === 0) return;

    let flash = false;
    let count = 0;
    const maxFlashes = 10; // 10 pulses * 300ms = 3000ms (3.0 seconds)

    blinkInterval = setInterval(() => {
        count++;
        flash = !flash;

        if (viewer3D) {
            if (flash) {
                viewer3D.setStyle({ resi: targetResidues }, {
                    cartoon: { color: '#fbbf24', thickness: 1.15, opacity: 1.0 },
                    stick: { color: '#fbbf24', radius: 0.36 },
                    sphere: { color: '#fbbf24', scale: 0.45 }
                });
            } else {
                viewer3D.setStyle({ resi: targetResidues }, {
                    cartoon: { color: '#0ea5e9', thickness: 0.65, opacity: 0.45 },
                    stick: { color: '#0ea5e9', radius: 0.16 },
                    sphere: { color: '#0ea5e9', scale: 0.22 }
                });
            }
            viewer3D.render();
        }

        if (count >= maxFlashes) {
            clearInterval(blinkInterval);
            blinkInterval = null;
            document.querySelectorAll('.blinking-pill').forEach(el => el.classList.remove('blinking-pill'));
            apply3DStyle();
        }
    }, 300);
}

// ── SINGLE RESIDUE INSPECTOR (DISPLAY ONE AT A TIME & TRIGGER 3D BLINK) ──
window.inspectResidue = function(resNum, aaChar) {
    if (!currentProtein) return;
    const aa = AA_DICT[aaChar] || {
        name: `Residue ${aaChar}`,
        code3: aaChar,
        category: 'Constituent Amino Acid',
        charge: 0.0,
        pKa: 'None',
        hydropathy: 0.0,
        role: 'Polypeptide constituent building block.'
    };

    // 1. Update Single Residue Card (Ensuring clean one-at-a-time presentation)
    const singleCard = document.getElementById('single-residue-card');
    if (singleCard) singleCard.style.display = 'block';

    const nameEl = document.getElementById('sres-name');
    if (nameEl) nameEl.textContent = `${aa.name} (${aa.code3}) — #${resNum}`;

    const catEl = document.getElementById('sres-cat');
    if (catEl) catEl.textContent = aa.category;

    const pkaEl = document.getElementById('sres-pka');
    if (pkaEl) pkaEl.textContent = aa.pKa;

    const chargeEl = document.getElementById('sres-charge');
    if (chargeEl) chargeEl.textContent = `${aa.charge > 0 ? '+' : ''}${aa.charge.toFixed(1)} e`;

    const hydroEl = document.getElementById('sres-hydro');
    if (hydroEl) hydroEl.textContent = `${aa.hydropathy > 0 ? '+' : ''}${aa.hydropathy.toFixed(2)}`;

    // Check if residue is contained within any annotated domain or functional region
    let domainContext = '';
    if (currentProtein.domains && currentProtein.domains.length > 0) {
        const dom = currentProtein.domains.find(d => resNum >= d.start && resNum <= d.end);
        if (dom) {
            domainContext = `[Located in ${dom.name} (${dom.start}–${dom.end})] `;
        }
    }

    const roleEl = document.getElementById('sres-role');
    if (roleEl) roleEl.textContent = domainContext + aa.role;

    // 2. Update sequence strip pill selection state
    document.querySelectorAll('.seq-residue-pill').forEach(p => p.classList.remove('active-pill'));
    const pill = document.getElementById(`seq-pill-${resNum}`);
    if (pill) {
        pill.classList.add('active-pill');
        pill.scrollIntoView({ behavior: 'smooth', inline: 'center', block: 'nearest' });
    }

    // 3. Safe 3D coordinate camera pinpointing
    if (viewer3D) {
        const selected = viewer3D.selectedAtoms({ resi: resNum });
        const infoChip = document.getElementById('viewer-info-chip');
        if (selected && selected.length > 0) {
            // Safe zoom to individual residue without second duration argument
            viewer3D.zoomTo({ resi: resNum });
            viewer3D.render();
            if (infoChip) {
                infoChip.textContent = `PINPOINT // ${aa.code3}${resNum} [3D COORDINATES LINKED]`;
            }
        } else {
            if (infoChip) {
                infoChip.textContent = `RESIDUE // ${aa.code3}${resNum} (UNRESOLVED IN SOLVED PDB)`;
            }
        }
    }

    // 4. Trigger synchronized 3-second blinking pulse
    startCoordinatedBlink([resNum]);
};

// ── RENDER SEQUENCE STRIP ──
function renderSequenceStrip(protein) {
    const container = document.getElementById('seq-strip-container');
    if (!container) return;
    container.innerHTML = '';

    const seq = protein.sequence || '';
    if (!seq || seq.length === 0) {
        container.innerHTML = '<div style="color:var(--text-muted); font-size:11px; padding:10px; font-style:italic;">Sequence data unavailable for this record.</div>';
        return;
    }

    const fragment = document.createDocumentFragment();
    for (let i = 0; i < seq.length; i++) {
        const char = seq[i];
        const resNum = i + 1;
        const aa = AA_DICT[char];

        let colorClass = 'res-polar';
        if (['K', 'R', 'H'].includes(char)) colorClass = 'res-pos';
        else if (['D', 'E'].includes(char)) colorClass = 'res-neg';
        else if (['A', 'V', 'I', 'L', 'M'].includes(char)) colorClass = 'res-hydro';
        else if (['F', 'W', 'Y'].includes(char)) colorClass = 'res-aro';
        else if (char === 'C') colorClass = 'res-cys';
        else if (char === 'P' || char === 'G') colorClass = 'res-special';

        const pill = document.createElement('div');
        pill.className = `seq-residue-pill ${colorClass}`;
        pill.id = `seq-pill-${resNum}`;
        pill.innerHTML = `
            <span class="pill-num">${resNum}</span>
            <span class="pill-letter">${char}</span>
        `;
        pill.title = `${aa ? aa.name : char} (#${resNum})`;
        pill.onclick = () => inspectResidue(resNum, char);
        fragment.appendChild(pill);
    }
    container.appendChild(fragment);

    // Initial inspection preview of the first residue
    if (seq.length > 0) {
        inspectResidue(1, seq[0]);
    }
}

// ── RENDER COMPREHENSIVE 5-PILLAR ARCHITECTURE SUMMARY ──
function renderComprehensiveSummary(protein) {
    if (!protein) return;

    // Synthesize 5-pillar biological architecture summary
    const s = protein.comprehensive_summary ||
              (window.BiophysicsEngine && window.BiophysicsEngine.computeComprehensiveSummary(protein)) || {};

    // ── 1. Concordance Card (Exact Percentage & Deltas) ──
    const comp = protein.structural_comparison || {};
    const concPct = s.concordance_percent !== undefined ? s.concordance_percent : (comp.concordance_percent || 90.3);
    const concRmsd = comp.rmsd_angstroms !== undefined ? comp.rmsd_angstroms : 0.82;
    const concIdent = comp.sequence_identity_percent !== undefined ? comp.sequence_identity_percent : 100.0;
    const modelRef = `${comp.predicted_model || 'AlphaFold'} / ${comp.experimental_id || 'Experimental'}`;

    const pctValEl = document.getElementById('concordance-pct-val');
    if (pctValEl) pctValEl.textContent = `${concPct.toFixed(1)}%`;

    const fillEl = document.getElementById('concordance-progress-fill');
    if (fillEl) fillEl.style.width = `${Math.min(100, Math.max(5, concPct))}%`;

    const rmsdEl = document.getElementById('concordance-rmsd-val');
    if (rmsdEl) rmsdEl.textContent = `${concRmsd.toFixed(2)} Å`;

    const identEl = document.getElementById('concordance-identity-val');
    if (identEl) identEl.textContent = `${concIdent.toFixed(1)}%`;

    const refEl = document.getElementById('concordance-ref-val');
    if (refEl) refEl.textContent = modelRef;

    // ── 2. Pillar 1: Structural Properties ──
    const p1 = s.structural_pillar || {};
    const primaryDesc = document.getElementById('sum-primary-desc');
    if (primaryDesc && p1.primary_structure) {
        primaryDesc.textContent = p1.primary_structure;
    }

    const sec = p1.secondary_structure || {};
    const helixVal = document.getElementById('sum-sec-helix');
    if (helixVal && sec.alpha_helix_percent !== undefined) {
        helixVal.textContent = `${sec.alpha_helix_percent.toFixed(1)}%`;
    }
    const sheetVal = document.getElementById('sum-sec-sheet');
    if (sheetVal && sec.beta_sheet_percent !== undefined) {
        sheetVal.textContent = `${sec.beta_sheet_percent.toFixed(1)}%`;
    }
    const loopsVal = document.getElementById('sum-sec-loops');
    if (loopsVal && sec.turns_loops_percent !== undefined) {
        loopsVal.textContent = `${sec.turns_loops_percent.toFixed(1)}%`;
    }

    const tertDesc = document.getElementById('sum-tertiary-desc');
    if (tertDesc && p1.tertiary_structure) {
        const t = p1.tertiary_structure;
        tertDesc.textContent = `${t.bonding_mechanisms} (${t.disulfide_bonds} disulfide bridges, ~${t.salt_bridges_estimated} salt bridges)`;
    }

    const quatDesc = document.getElementById('sum-quaternary-desc');
    if (quatDesc && p1.quaternary_structure) {
        quatDesc.textContent = `${p1.quaternary_structure.oligomeric_state} (${p1.quaternary_structure.subunits})`;
    }

    // ── 3. Pillar 2: Physicochemical Properties ──
    const p2 = s.physicochemical_pillar || {};
    const amphoDesc = document.getElementById('sum-ampho-desc');
    if (amphoDesc && p2.amphoteric_behavior) {
        amphoDesc.textContent = `${p2.amphoteric_behavior} (pI: ${p2.isoelectric_point?.toFixed(2) || '6.50'}, Net Charge at pH 7.4: ${p2.net_charge_ph74 > 0 ? '+' : ''}${p2.net_charge_ph74?.toFixed(2) || '0.00'} e)`;
    }

    const solDesc = document.getElementById('sum-solubility-desc');
    if (solDesc && p2.solubility_profile) {
        solDesc.textContent = `GRAVY: ${p2.gravy_index > 0 ? '+' : ''}${p2.gravy_index?.toFixed(3) || '-0.400'} — ${p2.solubility_profile} (MW: ${p2.molecular_weight_kda?.toFixed(1) || '30.0'} kDa)`;
    }

    const denatDesc = document.getElementById('sum-denat-desc');
    if (denatDesc && p2.denaturation_mechanism) {
        denatDesc.textContent = `${p2.denaturation_mechanism} (Estimated Tm: ${p2.estimated_tm?.toFixed(1) || '65.0'}°C)`;
    }

    // ── 4. Pillar 3: Functional Properties ──
    const p3 = s.functional_pillar || {};
    const specDesc = document.getElementById('sum-spec-desc');
    if (specDesc && p3.specificity_binding) {
        specDesc.textContent = p3.specificity_binding;
    }

    const catDesc = document.getElementById('sum-catalysis-desc');
    if (catDesc && p3.catalytic_active_site) {
        catDesc.textContent = p3.catalytic_active_site;
    }

    const ptmDesc = document.getElementById('sum-ptm-desc');
    if (ptmDesc && p3.ptm_capacity) {
        ptmDesc.textContent = `${p3.ptm_capacity} ${p3.conformational_flexibility ? '• ' + p3.conformational_flexibility : ''}`;
    }

    // ── 5. Pillar 4: Buffering Capacity ──
    const p4 = s.buffering_pillar || {};
    const bufDesc = document.getElementById('sum-buffer-desc');
    if (bufDesc && p4.fluid_buffering_role) {
        bufDesc.textContent = `${p4.fluid_buffering_role} (${p4.histidine_count || 0} Histidine residues, Buffering Capacity Score: ${p4.buffering_capacity_score?.toFixed(1) || '7.5'}/10.0)`;
    }
}

