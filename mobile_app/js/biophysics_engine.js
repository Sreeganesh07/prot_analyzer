/**
 * ProteinScope Mobile Biophysics & Thermodynamics Engine
 * Implements Alberts et al. (Molecular Biology of the Cell, Chapter 3: NBK26830)
 * Client-side high-performance execution (0 latency, 60-120 FPS).
 */

const PKA_VALUES = {
    N_term: 9.69,
    C_term: 2.34,
    D: 3.86,
    E: 4.25,
    H: 6.00,
    C: 8.33,
    Y: 10.07,
    K: 10.53,
    R: 12.48
};

class BiophysicsEngine {
    /**
     * Compute Henderson-Hasselbalch net molecular charge at given pH.
     */
    static calculateNetCharge(titratableCounts, pH) {
        const counts = titratableCounts || {};
        const nTerm = 1.0 / (1.0 + Math.pow(10, pH - PKA_VALUES.N_term));
        const kPos = (counts.K || 0) / (1.0 + Math.pow(10, pH - PKA_VALUES.K));
        const rPos = (counts.R || 0) / (1.0 + Math.pow(10, pH - PKA_VALUES.R));
        const hPos = (counts.H || 0) / (1.0 + Math.pow(10, pH - PKA_VALUES.H));
        const posCharge = nTerm + kPos + rPos + hPos;

        const cTerm = 1.0 / (1.0 + Math.pow(10, PKA_VALUES.C_term - pH));
        const dNeg = (counts.D || 0) / (1.0 + Math.pow(10, PKA_VALUES.D - pH));
        const eNeg = (counts.E || 0) / (1.0 + Math.pow(10, PKA_VALUES.E - pH));
        const cNeg = (counts.C || 0) / (1.0 + Math.pow(10, PKA_VALUES.C - pH));
        const yNeg = (counts.Y || 0) / (1.0 + Math.pow(10, PKA_VALUES.Y - pH));
        const negCharge = cTerm + dNeg + eNeg + cNeg + yNeg;

        return Number((posCharge - negCharge).toFixed(2));
    }

    /**
     * Calculate Isoelectric Point (pI) where net charge Q(pH) = 0.
     */
    static calculateIsoelectricPoint(titratableCounts) {
        let low = 2.0;
        let high = 13.0;
        for (let i = 0; i < 20; i++) {
            const mid = (low + high) / 2.0;
            const q = this.calculateNetCharge(titratableCounts, mid);
            if (q > 0) {
                low = mid;
            } else {
                high = mid;
            }
        }
        return Number(((low + high) / 2.0).toFixed(2));
    }

    /**
     * Sequence-specific Melting Temperature Tm estimation.
     */
    static estimateMeltingTemperature(sequence, properties) {
        if (!sequence) return 65.0;
        const seq = sequence.toUpperCase();
        const L = seq.length || 300;

        const aCount = (seq.match(/A/g) || []).length;
        const vCount = (seq.match(/V/g) || []).length;
        const iCount = (seq.match(/I/g) || []).length;
        const lCount = (seq.match(/L/g) || []).length;
        const pCount = (seq.match(/P/g) || []).length;
        const cCount = (seq.match(/C/g) || []).length;

        const aliphIndex = ((aCount * 1.0 + vCount * 2.9 + (iCount + lCount) * 3.9) / L) * 100.0;
        const proPct = (pCount / L) * 100.0;
        const cysPct = (cCount / L) * 100.0;

        let tm = 45.0 + 0.15 * aliphIndex + 0.3 * proPct + 0.4 * cysPct;
        if (properties && properties.instability_index) {
            tm -= (properties.instability_index - 40.0) * 0.1;
        }
        return Number(Math.max(45.0, Math.min(95.0, tm)).toFixed(1));
    }

    /**
     * Calculate Salt Bridge Retention % as a function of pH.
     */
    static calculateSaltBridgeRetention(pH) {
        const pKaAcid = 4.0;
        const pKaBase = 11.0;
        const fAcid = 1.0 / (1.0 + Math.pow(10, pKaAcid - pH));
        const fBase = 1.0 / (1.0 + Math.pow(10, pH - pKaBase));
        return Number((fAcid * fBase * 100.0).toFixed(1));
    }

    /**
     * Two-State Folding Thermodynamics (Alberts NBK26830).
     */
    static evaluateState(tempC, pH, tm, titratableCounts, thresholds) {
        const T = tempC + 273.15;
        const Tm_K = tm + 273.15;
        const R = 0.001987; // kcal/(mol·K)

        // Enthalpy & Heat Capacity
        const dH_Tm = 100.0; // kcal/mol baseline
        const dCp = 1.8;     // kcal/(mol·K)

        const dH_T = dH_Tm + dCp * (T - Tm_K);
        const dS_T = (dH_Tm / Tm_K) + dCp * Math.log(T / Tm_K);
        let dG_thermal = -(dH_T - T * dS_T);

        // Electrostatic repulsion penalty
        const q = this.calculateNetCharge(titratableCounts, pH);
        const dG_elec = 0.008 * Math.pow(q, 2);
        const dG_total = dG_thermal + dG_elec;

        // Fraction Folded
        const expArg = Math.max(-40.0, Math.min(40.0, -dG_total / (R * T)));
        const K_fold = Math.exp(expArg);
        const fractionFolded = Number(((K_fold / (1.0 + K_fold)) * 100.0).toFixed(1));
        const saltRetention = this.calculateSaltBridgeRetention(pH);

        const th = thresholds || {
            temp_min_stable: 15.0,
            temp_max_stable: tm - 8.0,
            temp_degradation: Math.max(70.0, tm + 10.0),
            ph_min_stable: 5.5,
            ph_max_stable: 8.5,
            ph_acid_degradation: 3.0,
            ph_alkaline_degradation: 11.5
        };

        // Determine Degradation State
        let state = "NATIVE";
        let badge = "Native Folded State";
        let color = "#10b981";
        let isDegraded = false;
        let reasons = [];
        let mechanism = "";

        if (tempC >= th.temp_degradation) {
            state = "DEGRADED";
            badge = "Thermal Degradation / Amyloid Aggregation";
            color = "#ef4444";
            isDegraded = true;
            reasons.push(`Temperature ${tempC}°C exceeds thermal degradation boundary (${th.temp_degradation}°C)`);
            mechanism = "High thermal kinetic energy breaks backbone hydrogen bonds and noncovalent interactions, exposing hydrophobic cores leading to irreversible precipitation.";
        } else if (pH <= th.ph_acid_degradation) {
            state = "DEGRADED";
            badge = "Acid Hydrolysis / Unfolded";
            color = "#ef4444";
            isDegraded = true;
            reasons.push(`pH ${pH} reaches acid degradation threshold (≤ ${th.ph_acid_degradation})`);
            mechanism = "Extreme acidity protonates aspartic and glutamic acid carboxyl groups, abolishing stabilizing salt bridges and creating massive electrostatic repulsion.";
        } else if (pH >= th.ph_alkaline_degradation) {
            state = "DEGRADED";
            badge = "Alkaline Denaturation / Aggregation";
            color = "#ef4444";
            isDegraded = true;
            reasons.push(`pH ${pH} exceeds alkaline threshold (≥ ${th.ph_alkaline_degradation})`);
            mechanism = "Severe alkalinity deprotonates lysine and arginine basic side chains, stripping positive charges and unraveling tertiary folding.";
        } else if (fractionFolded < 35.0 || tempC >= tm) {
            state = "DENATURED";
            badge = "Denatured Random Coil";
            color = "#f97316";
            mechanism = "The polypeptide chain has lost its specific tertiary structure due to thermal disruption of weak bonds, existing as a dynamic random coil.";
        } else if (fractionFolded < 80.0 || tempC > th.temp_max_stable || pH < th.ph_min_stable || pH > th.ph_max_stable) {
            state = "PERTURBED";
            badge = "Perturbed / Molten Globule";
            color = "#f59e0b";
            mechanism = "Sub-optimal physiological condition induces conformation breathing; native contacts loosen while maintaining overall secondary structure.";
        } else {
            state = "NATIVE";
            badge = "Native Folded State";
            color = "#10b981";
            mechanism = "At standard/physiological conditions, noncovalent bonds (H-bonds, hydrophobic clustering, ionic salt bridges) hold the protein in its minimum free energy conformation.";
        }

        return {
            temperature_celsius: tempC,
            ph: pH,
            state: state,
            badge_label: badge,
            badge_color: color,
            is_degraded: isDegraded,
            reasons: reasons,
            mechanism: mechanism,
            net_charge: q,
            fraction_folded_percent: fractionFolded,
            salt_bridge_retention_percent: saltRetention,
            delta_g_folding_kcal_mol: Number(dG_total.toFixed(2)),
            thresholds: th
        };
    }

    /**
     * Live NCBI Entrez Protein database search and sequence retrieval.
     */
    static async fetchNCBIProtein(query) {
        const cleanQuery = query.trim();
        if (!cleanQuery) throw new Error("Query cannot be empty.");

        // 1. Search NCBI Protein DB
        const searchUrl = `https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=protein&term=${encodeURIComponent(cleanQuery)}&retmode=json&retmax=1`;
        const searchRes = await fetch(searchUrl);
        if (!searchRes.ok) throw new Error(`NCBI Search failed (HTTP ${searchRes.status})`);

        const searchData = await searchRes.json();
        const idList = searchData?.esearchresult?.idlist || [];
        if (idList.length === 0) {
            throw new Error(`No NCBI protein records found for '${cleanQuery}'`);
        }
        const ncbiId = idList[0];

        // 2. Fetch FASTA
        const fetchUrl = `https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=${ncbiId}&rettype=fasta&retmode=text`;
        const fastaRes = await fetch(fetchUrl);
        if (!fastaRes.ok) throw new Error("Failed to fetch FASTA sequence from NCBI.");

        const fastaText = (await fastaRes.text()).trim();
        const lines = fastaText.split("\n");
        const header = lines[0] || "";
        const sequence = lines.slice(1).join("").replace(/\s+/g, "").toUpperCase();

        // 3. Fetch Summary for title/organism
        let proteinName = cleanQuery;
        let organism = "Unknown";
        try {
            const sumUrl = `https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=protein&id=${ncbiId}&retmode=json`;
            const sumRes = await fetch(sumUrl);
            if (sumRes.ok) {
                const sumData = await sumRes.json();
                const obj = sumData?.result?.[ncbiId];
                if (obj && obj.title) {
                    proteinName = obj.title.split("[")[0].trim();
                    if (obj.title.includes("[") && obj.title.includes("]")) {
                        organism = obj.title.split("[").pop().split("]")[0].trim();
                    }
                }
            }
        } catch (e) {
            console.warn("Could not fetch NCBI summary:", e);
        }

        // Count titratable residues
        const titratableCounts = {
            D: (sequence.match(/D/g) || []).length,
            E: (sequence.match(/E/g) || []).length,
            H: (sequence.match(/H/g) || []).length,
            C: (sequence.match(/C/g) || []).length,
            Y: (sequence.match(/Y/g) || []).length,
            K: (sequence.match(/K/g) || []).length,
            R: (sequence.match(/R/g) || []).length,
            len: sequence.length
        };

        const tmEst = this.estimateMeltingTemperature(sequence, {});
        const piEst = this.calculateIsoelectricPoint(titratableCounts);

        // Check if AlphaFold model exists for query
        let pdbContent = "";
        try {
            const afUrl = `https://alphafold.ebi.ac.uk/files/AF-${cleanQuery.toUpperCase()}-F1-model_v4.pdb`;
            const afRes = await fetch(afUrl);
            if (afRes.ok) {
                pdbContent = await afRes.text();
            }
        } catch (e) {
            console.debug("AlphaFold DB check skipped:", e);
        }

        return {
            accession: cleanQuery.toUpperCase(),
            name: proteinName,
            organism: organism,
            sequence: sequence,
            sequence_length: sequence.length,
            titratable_counts: titratableCounts,
            isoelectric_point: piEst,
            estimated_tm: tmEst,
            thresholds: {
                tm_celsius: tmEst,
                temp_min_stable: 15.0,
                temp_max_stable: Number((tmEst - 8.0).toFixed(1)),
                temp_degradation: Number(Math.max(70.0, tmEst + 10.0).toFixed(1)),
                ph_min_stable: 5.5,
                ph_max_stable: 8.5,
                ph_acid_degradation: 3.0,
                ph_alkaline_degradation: 11.5
            },
            pdb_content: pdbContent
        };
    }
}

window.BiophysicsEngine = BiophysicsEngine;
