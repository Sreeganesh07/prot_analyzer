/**
 * ProteinScope Mobile Application Controller
 * Manages 3Dmol viewer, offline protein database, client-side biophysics, and NCBI search.
 */

// BiophysicsEngine is attached to window
const BiophysicsEngine = window.BiophysicsEngine;

let currentProtein = null;
let viewer3D = null;
let renderMode = 'cartoon';
let isDegradedVisual = false;

// Initialize on DOM ready
document.addEventListener('DOMContentLoaded', () => {
    initCanisterTrack();
    setupEventListeners();

    // Load default protein (p53)
    const db = window.PROTEINSCOPE_CACHED_DB || {};
    const firstAcc = Object.keys(db)[0] || 'P04637';
    loadProtein(firstAcc);
});

function initCanisterTrack() {
    const track = document.getElementById('canister-track');
    if (!track) return;
    track.innerHTML = '';

    const db = window.PROTEINSCOPE_CACHED_DB || {};
    Object.keys(db).forEach(acc => {
        const p = db[acc];
        const btn = document.createElement('button');
        btn.className = 'canister-btn';
        btn.id = `canister-${acc}`;
        btn.innerHTML = `
            <span class="can-acc">${p.accession}</span>
            <span class="can-name">${p.name.split('(')[0].trim()}</span>
            <span class="can-tm">${p.estimated_tm.toFixed(1)}°C</span>
        `;
        btn.onclick = () => loadProtein(acc);
        track.appendChild(btn);
    });

    // Add NCBI Search toggle button at end of track
    const ncbiBtn = document.createElement('button');
    ncbiBtn.className = 'canister-btn ncbi-tab-btn';
    ncbiBtn.innerHTML = `
        <span class="can-acc">+ NCBI</span>
        <span class="can-name">Online Query</span>
        <span class="can-tm">ENTREZ DB</span>
    `;
    ncbiBtn.onclick = toggleNCBIDock;
    track.appendChild(ncbiBtn);
}

function toggleNCBIDock() {
    const dock = document.getElementById('ncbi-dock');
    if (!dock) return;
    dock.classList.toggle('active');
    if (dock.classList.contains('active')) {
        document.getElementById('ncbi-input').focus();
    }
}

function setupEventListeners() {
    const sliderTemp = document.getElementById('slider-temp');
    const sliderPH = document.getElementById('slider-ph');

    if (sliderTemp) {
        sliderTemp.addEventListener('input', (e) => {
            const val = parseFloat(e.target.value);
            document.getElementById('val-temp').innerHTML = `${val.toFixed(1)}°C <span style="font-size:10px; color:var(--text-muted);">(${(val + 273.15).toFixed(2)} K)</span>`;
            recalculateStability();
        });
    }

    if (sliderPH) {
        sliderPH.addEventListener('input', (e) => {
            const val = parseFloat(e.target.value);
            document.getElementById('val-ph').textContent = `pH ${val.toFixed(2)}`;
            recalculateStability();
        });
    }

    const ncbiFetchBtn = document.getElementById('ncbi-fetch-btn');
    if (ncbiFetchBtn) {
        ncbiFetchBtn.addEventListener('click', handleNCBIFetch);
    }
    const ncbiInput = document.getElementById('ncbi-input');
    if (ncbiInput) {
        ncbiInput.addEventListener('keypress', (e) => {
            if (e.key === 'Enter') handleNCBIFetch();
        });
    }
}

// Global functions for HTML onclick
window.setTemp = function(temp) {
    const slider = document.getElementById('slider-temp');
    if (slider) {
        slider.value = temp;
        document.getElementById('val-temp').innerHTML = `${temp.toFixed(1)}°C <span style="font-size:10px; color:var(--text-muted);">(${(temp + 273.15).toFixed(2)} K)</span>`;
        recalculateStability();
    }
};

window.setPH = function(ph) {
    const slider = document.getElementById('slider-ph');
    if (slider) {
        slider.value = ph;
        document.getElementById('val-ph').textContent = `pH ${ph.toFixed(2)}`;
        recalculateStability();
    }
};

window.setRenderMode = function(mode) {
    renderMode = mode;
    document.querySelectorAll('.tool-btn').forEach(b => {
        if (b.textContent.toLowerCase() === mode) b.classList.add('active');
        else if (b.textContent.toLowerCase() === 'cartoon' || b.textContent.toLowerCase() === 'surface') b.classList.remove('active');
    });
    apply3DStyle();
};

window.reset3DView = function() {
    if (viewer3D) {
        viewer3D.zoomTo();
        viewer3D.render();
    }
};

function loadProtein(acc) {
    const db = window.PROTEINSCOPE_CACHED_DB || {};
    const p = db[acc];
    if (!p) return;

    currentProtein = p;

    // Update active canister highlight
    document.querySelectorAll('.canister-btn').forEach(b => b.classList.remove('active'));
    const activeBtn = document.getElementById(`canister-${acc}`);
    if (activeBtn) activeBtn.classList.add('active');

    // Update specimen header
    document.getElementById('disp-acc').textContent = p.accession;
    document.getElementById('disp-name').textContent = p.name;
    document.getElementById('disp-org').textContent = p.organism;
    document.getElementById('disp-pi').textContent = `pI: ${p.isoelectric_point}`;
    document.getElementById('tele-tm').textContent = `${p.estimated_tm.toFixed(1)}°C`;

    // Render 3D structure
    init3DViewer(p.pdb_content);

    // Initial recalculate
    recalculateStability();
}

function init3DViewer(pdbText) {
    const container = document.getElementById('mol-canvas');
    if (!container) return;
    container.innerHTML = '';

    try {
        if (typeof $3Dmol === 'undefined') {
            container.innerHTML = '<div style="color:var(--text-muted); padding:30px; text-align:center; font-family:var(--font-mono); font-size:11px;">WebGL 3Dmol Library Initializing...</div>';
            return;
        }

        viewer3D = $3Dmol.createViewer(container, {
            backgroundColor: '#07090e',
            defaultcolors: $3Dmol.rasmolElementColors
        });

        if (pdbText && pdbText.trim().length > 0) {
            viewer3D.addModel(pdbText, 'pdb');
            apply3DStyle();
            viewer3D.zoomTo();
            viewer3D.render();
            document.getElementById('viewer-status-text').textContent = `3DMOL // ${currentProtein.accession} PDB LOADED`;
        } else {
            container.innerHTML = '<div style="color:var(--text-muted); padding:40px 20px; text-align:center; font-family:var(--font-mono); font-size:11px;">NO EXPERIMENTAL PDB CO-CRYSTAL AVAILABLE<br><span style="color:var(--accent-cyan);">SEQUENCE BIOPHYSICAL ANALYSIS ACTIVE</span></div>';
            document.getElementById('viewer-status-text').textContent = `BIOPHYSICS COMPUTED // NO 3D COORDINATES`;
        }
    } catch (e) {
        console.error('Error in 3D viewer initialization:', e);
    }
}

function apply3DStyle() {
    if (!viewer3D) return;
    viewer3D.setStyle({}, {}); // Clear

    if (isDegradedVisual) {
        // Disorganized, aggregated appearance
        viewer3D.setStyle({}, {
            stick: { radius: 0.15, color: '#ef4444' },
            sphere: { radius: 0.4, color: '#991b1b' }
        });
    } else if (renderMode === 'surface') {
        viewer3D.setStyle({}, {
            cartoon: { color: 'spectrum', opacity: 0.4 }
        });
        viewer3D.addSurface($3Dmol.SurfaceType.VDW, {
            opacity: 0.65,
            colorscheme: 'hydrophobicity'
        });
    } else {
        // High-contrast secondary structure coloring
        viewer3D.setStyle({ ss: 'h' }, { cartoon: { color: '#0284c7', thickness: 0.6 } }); // Helices
        viewer3D.setStyle({ ss: 's' }, { cartoon: { color: '#f59e0b', thickness: 0.6 } }); // Sheets
        viewer3D.setStyle({ ss: 'c' }, { cartoon: { color: '#64748b', thickness: 0.4 } }); // Coils
    }
    viewer3D.render();
}

function recalculateStability() {
    if (!currentProtein) return;

    const temp = parseFloat(document.getElementById('slider-temp').value) || 37.0;
    const ph = parseFloat(document.getElementById('slider-ph').value) || 7.4;

    const evalResult = BiophysicsEngine.evaluateState(
        temp,
        ph,
        currentProtein.estimated_tm,
        currentProtein.titratable_counts,
        currentProtein.thresholds
    );

    // Update Telemetry Display
    document.getElementById('tele-folded').textContent = `${evalResult.fraction_folded_percent.toFixed(1)}%`;
    document.getElementById('tele-folded').style.color = evalResult.fraction_folded_percent > 80 ? 'var(--accent-emerald)' : (evalResult.fraction_folded_percent > 40 ? 'var(--accent-amber)' : 'var(--accent-crimson)');

    document.getElementById('tele-charge').textContent = `${evalResult.net_charge > 0 ? '+' : ''}${evalResult.net_charge.toFixed(2)} e`;
    document.getElementById('tele-salt').textContent = `${evalResult.salt_bridge_retention_percent.toFixed(1)}%`;
    document.getElementById('tele-dg').textContent = `${evalResult.delta_g_folding_kcal_mol.toFixed(1)} kcal`;

    // Update State Badge
    const stateBadge = document.getElementById('disp-state-badge');
    const stateDot = document.getElementById('disp-state-dot');
    const stateText = document.getElementById('disp-state-text');
    stateBadge.style.color = evalResult.badge_color;
    stateDot.style.background = evalResult.badge_color;
    stateText.textContent = evalResult.badge_label;

    // Update Degradation Alarm
    const banner = document.getElementById('deg-banner');
    const viewerCard = document.getElementById('viewer-container');
    if (evalResult.is_degraded) {
        banner.style.display = 'block';
        document.getElementById('deg-badge').textContent = evalResult.badge_label;
        document.getElementById('deg-desc').textContent = evalResult.reasons.join(' • ') + ' — ' + evalResult.mechanism;
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

    // Draw Henderson-Hasselbalch Titration Curve
    renderTitrationCurve(ph, evalResult.net_charge);
}

function renderTitrationCurve(currentPH, currentQ) {
    if (!currentProtein) return;
    const path = document.getElementById('svg-curve');
    const dot = document.getElementById('svg-dot');
    if (!path || !dot) return;

    // Compute curve points from pH 0 to 14
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

async function handleNCBIFetch() {
    const input = document.getElementById('ncbi-input');
    const feedback = document.getElementById('ncbi-feedback');
    const query = input.value.trim();

    if (!query) {
        feedback.style.display = 'block';
        feedback.style.color = 'var(--accent-amber)';
        feedback.textContent = 'Please enter an NCBI accession code (e.g. NP_000537) or protein name.';
        return;
    }

    feedback.style.display = 'block';
    feedback.style.color = 'var(--accent-cyan)';
    feedback.textContent = `Connecting to NCBI Entrez (eutils.ncbi.nlm.nih.gov)... Fetching '${query}'`;

    try {
        const fetchedData = await BiophysicsEngine.fetchNCBIProtein(query);
        feedback.style.color = 'var(--accent-emerald)';
        feedback.textContent = `Retrieved: ${fetchedData.name} (${fetchedData.sequence_length} aa). Estimated Tm: ${fetchedData.estimated_tm}°C`;

        // Cache in memory
        window.PROTEINSCOPE_CACHED_DB[fetchedData.accession] = fetchedData;

        // Add canister tab if not already present
        let canBtn = document.getElementById(`canister-${fetchedData.accession}`);
        if (!canBtn) {
            const track = document.getElementById('canister-track');
            canBtn = document.createElement('button');
            canBtn.className = 'canister-btn';
            canBtn.id = `canister-${fetchedData.accession}`;
            canBtn.innerHTML = `
                <span class="can-acc">${fetchedData.accession}</span>
                <span class="can-name">${fetchedData.name.split('(')[0].trim()}</span>
                <span class="can-tm">${fetchedData.estimated_tm.toFixed(1)}°C</span>
            `;
            canBtn.onclick = () => loadProtein(fetchedData.accession);
            track.insertBefore(canBtn, track.lastElementChild);
        }

        // Load immediately
        loadProtein(fetchedData.accession);

        setTimeout(() => {
            document.getElementById('ncbi-dock').classList.remove('active');
        }, 1800);
    } catch (err) {
        feedback.style.color = 'var(--accent-crimson)';
        feedback.textContent = `NCBI Error: ${err.message}`;
    }
}
