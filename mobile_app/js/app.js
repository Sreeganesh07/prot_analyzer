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
        if (btn.textContent.toLowerCase() === mode.toLowerCase()) {
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
        viewer3D.spin('y', 1.0);
    } else {
        viewer3D.spin(false);
    }
};

window.reset3DCamera = function() {
    if (viewer3D) {
        viewer3D.zoomTo();
        viewer3D.render();
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

    // Update Properties Tab Values
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

    // Render 3D molecular structure
    init3DViewer(p.pdb_content);

    // Recalculate biophysics stability
    recalculateStability();
}

// ── INITIALIZE 3D VIEWER (OFFLINE LOCAL 3DMOL) ──
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
        viewer3D = $3Dmol.createViewer(container, {
            backgroundColor: '#050811',
            defaultcolors: $3Dmol.rasmolElementColors
        });

        if (pdbText && pdbText.trim().length > 0) {
            viewer3D.addModel(pdbText, 'pdb');
            apply3DStyle();
            viewer3D.zoomTo();
            viewer3D.render();

            if (infoChip) {
                infoChip.textContent = `3DMOL // ${currentProtein.accession} PDB ACTIVE`;
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
            viewer3D.spin('y', 1.0);
        }
    } catch (e) {
        console.error('Error in 3D viewer initialization:', e);
        container.innerHTML = `
            <div style="color:var(--accent-crimson); padding:30px; text-align:center; font-family:var(--font-mono); font-size:11px;">
                WebGL Initialization Error: ${e.message}
            </div>`;
    }
}

// ── APPLY 3D VISUAL STYLES ──
function apply3DStyle() {
    if (!viewer3D) return;
    viewer3D.setStyle({}, {}); // Clear previous styles

    if (isDegradedVisual) {
        // Disorganized, aggregated appearance in red/crimson
        viewer3D.setStyle({}, {
            stick: { radius: 0.16, color: '#ef4444' },
            sphere: { radius: 0.42, color: '#991b1b' }
        });
    } else if (renderMode === 'surface') {
        viewer3D.setStyle({}, {
            cartoon: { color: 'spectrum', opacity: 0.35 }
        });
        viewer3D.addSurface($3Dmol.SurfaceType.VDW, {
            opacity: 0.65,
            colorscheme: 'hydrophobicity'
        });
    } else if (renderMode === 'spheres') {
        viewer3D.setStyle({}, {
            sphere: { colorscheme: 'amino', scale: 0.32 }
        });
    } else {
        // Default: High-contrast Secondary Structure Cartoon Ribbon
        viewer3D.setStyle({ ss: 'h' }, { cartoon: { color: '#0ea5e9', thickness: 0.65 } }); // Helices (Cyan/Blue)
        viewer3D.setStyle({ ss: 's' }, { cartoon: { color: '#f59e0b', thickness: 0.65 } }); // Sheets (Gold)
        viewer3D.setStyle({ ss: 'c' }, { cartoon: { color: '#64748b', thickness: 0.40 } }); // Coils (Slate)
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
        feedback.style.color = 'var(--accent-emerald)';
        feedback.textContent = `✓ Retrieved: ${fetchedData.name} (${fetchedData.sequence_length} aa). Estimated Tm: ${fetchedData.estimated_tm}°C`;

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
