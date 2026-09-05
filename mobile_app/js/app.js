/**
 * ProteinScope Mobile Application Controller
 * Manages 3Dmol WebGL viewer, offline database, client-side biophysics, and online search.
 */


let currentProtein = null;
let viewer3D = null;
let renderMode = 'cartoon';
let isDegradedVisual = false;
let isSpinning = false;

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
    if (!viewer3D || !currentProtein || !currentProtein.domains || !currentProtein.domains[domainIdx]) return;
    const dom = currentProtein.domains[domainIdx];
    
    // Switch to 3D viewer tab so user immediately sees the visual focus
    switchTab('structure');

    // Remove any surface meshes and reset
    viewer3D.removeAllSurfaces();
    viewer3D.removeAllShapes();

    // Mute non-domain residues into translucent ghost slate
    viewer3D.setStyle({}, {
        cartoon: { color: '#334155', opacity: 0.28, thickness: 0.40 }
    });

    // Highlight selected domain with bold glowing color
    const resiRange = Array.from({ length: dom.end - dom.start + 1 }, (_, i) => dom.start + i);
    viewer3D.setStyle(
        { resi: resiRange },
        { cartoon: { color: dom.color, thickness: 0.85, opacity: 1.0 } }
    );

    // Zoom and center directly onto this domain
    viewer3D.zoomTo({ resi: resiRange }, 600);
    viewer3D.render();

    const infoChip = document.getElementById('viewer-info-chip');
    if (infoChip) {
        infoChip.textContent = `FOCUS // ${dom.name.toUpperCase()} (${dom.start}-${dom.end})`;
    }
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
