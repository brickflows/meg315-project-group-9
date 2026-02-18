/* ═══════════════════════════════════════════════════════════════
   AD-HTC Fuel-Enhanced Power Gas Cycle — Thermodynamic Engine
   ═══════════════════════════════════════════════════════════════
   AD  = Anaerobic Digestion (processes moisture-rich biomass → biogas)
   HTC = Hydrothermal Carbonization (processes moisture-lean biomass)
   Gas Cycle = Brayton cycle fueled by biogas from AD
   Steam Cycle = Rankine cycle driven by HTC process heat
   ═══════════════════════════════════════════════════════════════ */

// ────────────────────────────────────────────────
//  1. GAS PROPERTY TABLE API (Air / Combustion Gas)
// ────────────────────────────────────────────────

const GasTable = (() => {
    const R_air = 0.287; // kJ/(kg·K)

    // Polynomial cp(T) for air [kJ/(kg·K)] — valid 250–2000 K
    function cp(T) {
        const t = T / 1000;
        return 1.0483 - 0.3717 * t + 0.9483 * t * t
            - 0.6271 * t * t * t + 0.1507 * t * t * t * t;
    }

    // Enthalpy h(T) relative to 0 K [kJ/kg] — Simpson integration
    function h(T) {
        const n = 200;
        const dT = T / n;
        let sum = cp(1); // avoid T=0
        for (let i = 1; i < n; i++) {
            sum += (i % 2 === 0 ? 2 : 4) * cp(i * dT);
        }
        sum += cp(T);
        return (dT / 3) * sum;
    }

    // Entropy s(T, P) relative to T_ref=298.15K, P_ref=101.325 kPa [kJ/(kg·K)]
    function s(T, P) {
        const T_ref = 298.15, P_ref = 101.325;
        const n = 200;
        const dT = (T - T_ref) / n;
        let sum = cp(T_ref) / T_ref;
        for (let i = 1; i < n; i++) {
            const Ti = T_ref + i * dT;
            sum += (i % 2 === 0 ? 2 : 4) * (cp(Ti) / Ti);
        }
        sum += cp(T) / T;
        return (dT / 3) * sum - R_air * Math.log(P / P_ref);
    }

    // Isentropic T2 for given T1, P1→P2 (bisection)
    function T_isentropic(T1, P1, P2) {
        const s_target = s(T1, P1);
        let lo = T1 * 0.4, hi = T1 * 3.5;
        for (let i = 0; i < 80; i++) {
            const mid = (lo + hi) / 2;
            if (s(mid, P2) < s_target) lo = mid; else hi = mid;
        }
        return (lo + hi) / 2;
    }

    function gamma(T) { const c = cp(T); return c / (c - R_air); }

    return { cp, h, s, T_isentropic, gamma, R: R_air };
})();


// ────────────────────────────────────────────────
//  2. STEAM TABLE API (IAPWS-IF97 Approximations)
// ────────────────────────────────────────────────

const SteamTable = (() => {
    // Saturation temperature from pressure [kPa → K]
    function T_sat(P_kPa) {
        const P = P_kPa / 1000; // MPa
        const n1 = 1167.0521452767, n3 = -17.073846940092, n6 = 14.91510861353;
        const n4 = 12020.82470247, n7 = -4823.2657361591;
        const n2 = -724213.16703206, n5 = -3232555.0322333, n8 = 405113.40542057;
        const n9 = -0.23855557567849, n10 = 650.17534844798;
        const beta = Math.pow(P, 0.25);
        const E = beta * beta + n3 * beta + n6;
        const F = n1 * beta * beta + n4 * beta + n7;
        const G = n2 * beta * beta + n5 * beta + n8;
        const D = 2 * G / (-F - Math.sqrt(F * F - 4 * E * G));
        return (n10 + D - Math.sqrt((n10 + D) * (n10 + D) - 4 * (n9 + n10 * D))) / 2;
    }

    function P_sat(T) {
        let lo = 0.5, hi = 25000;
        for (let i = 0; i < 60; i++) {
            const mid = (lo + hi) / 2;
            if (T_sat(mid) < T) lo = mid; else hi = mid;
        }
        return (lo + hi) / 2;
    }

    function hf(P) { const Tc = T_sat(P) - 273.15; return 4.18 * Tc + 0.00088 * Tc * Tc; }
    function hg(P) { return 2501.3 + 1.86 * (T_sat(P) - 273.15); }
    function hfg(P) { return hg(P) - hf(P); }
    function sf(P) { return 4.18 * Math.log(T_sat(P) / 273.15); }
    function sg(P) { return sf(P) + hfg(P) / T_sat(P); }
    function vf() { return 0.001; }

    function h_super(P, T) {
        const cp_s = 1.87 + 0.00036 * (T - 373.15);
        return hg(P) + cp_s * (T - T_sat(P));
    }
    function s_super(P, T) {
        const cp_s = 1.87 + 0.00036 * (T - 373.15);
        return sg(P) + cp_s * Math.log(T / T_sat(P));
    }
    function T_from_s_super(P, s_target) {
        let lo = T_sat(P), hi = 1200;
        for (let i = 0; i < 80; i++) {
            const mid = (lo + hi) / 2;
            if (s_super(P, mid) < s_target) lo = mid; else hi = mid;
        }
        return (lo + hi) / 2;
    }
    function x_from_s(P, sv) {
        const sfv = sf(P), sgv = sg(P);
        if (sv <= sfv) return 0;
        if (sv >= sgv) return 1;
        return (sv - sfv) / (sgv - sfv);
    }
    function h_from_x(P, x) { return hf(P) + x * hfg(P); }

    function saturationDome(n = 50) {
        const dome = [];
        for (let i = 0; i <= n; i++) {
            const P = 5 + (18000 - 5) * (i / n);
            dome.push({ sf: sf(P), hf: hf(P), sg: sg(P), hg: hg(P), P, T: T_sat(P) });
        }
        return dome;
    }

    return { T_sat, P_sat, hf, hg, hfg, sf, sg, vf, h_super, s_super, T_from_s_super, x_from_s, h_from_x, saturationDome };
})();


// ────────────────────────────────────────────────
//  3. CYCLE CALCULATIONS
// ────────────────────────────────────────────────

function getParams() {
    return {
        T1: +document.getElementById('T1').value,
        P1: +document.getElementById('P1').value,
        rp: +document.getElementById('rp').value,
        TIT: +document.getElementById('TIT').value,
        eta_c: +document.getElementById('eta_c').value,
        eta_t: +document.getElementById('eta_t').value,
        LHV: +document.getElementById('LHV').value,
        m_air: +document.getElementById('m_air').value,
        P_boiler: +document.getElementById('P_boiler').value * 1000,
        T_steam: +document.getElementById('T_steam').value,
        P_cond: +document.getElementById('P_cond').value,
        eta_st: +document.getElementById('eta_st').value,
        eta_fp: +document.getElementById('eta_fp').value,
        m_biomass: +document.getElementById('m_biomass').value,
        moisture_split: +document.getElementById('moisture_split').value,
        ad_yield: +document.getElementById('ad_yield').value,
        htc_temp: +document.getElementById('htc_temp').value,
    };
}


/* ── Gas (Brayton) Cycle ── */
function calculateGasCycle(p) {
    const states = [];

    // State 1: Air inlet (ambient)
    const T1 = p.T1, P1 = p.P1;
    states.push({ id: '①', label: 'Air Inlet', T: T1, P: P1, h: GasTable.h(T1), s: GasTable.s(T1, P1) });

    // State 2: After Compressor (isentropic + η_c)
    const P2 = P1 * p.rp;
    const T2s = GasTable.T_isentropic(T1, P1, P2);
    const h1 = GasTable.h(T1), h2s = GasTable.h(T2s);
    const h2 = h1 + (h2s - h1) / p.eta_c;
    const T2 = T2s + (T2s - T1) * (1 / p.eta_c - 1);
    states.push({ id: '②', label: 'Compressor Exit', T: T2, P: P2, h: h2, s: GasTable.s(T2, P2) });

    // State 3: After Biogas Combustion Chamber (heat addition to TIT)
    const T3 = p.TIT;
    const P3 = P2 * 0.96; // ~4% pressure drop in combustor
    states.push({ id: '③', label: 'Combustor Inlet', T: T2, P: P2, h: h2, s: GasTable.s(T2, P2) });

    // State 4: Combustor exit = TIT
    const h4 = GasTable.h(T3);
    states.push({ id: '④', label: 'Combustor Exit (TIT)', T: T3, P: P3, h: h4, s: GasTable.s(T3, P3) });

    // State 5: After Turbine expansion to ~P1
    const P5 = P1 * 1.02;
    const T5s = GasTable.T_isentropic(T3, P3, P5);
    const h5s = GasTable.h(T5s);
    const h5 = h4 - p.eta_t * (h4 - h5s);
    const T5 = T5s + (T3 - T5s) * (1 - p.eta_t);
    states.push({ id: '⑤', label: 'Turbine Exhaust', T: T5, P: P5, h: h5, s: GasTable.s(T5, P5) });

    // Work & heat
    const w_c = h2 - h1;      // compressor work (input)
    const w_t = h4 - h5;      // turbine work (output)
    const q_in = h4 - h2;      // heat from biogas combustion
    const w_net = w_t - w_c;

    // Fuel (biogas) mass flow needed
    const m_fuel = (q_in * p.m_air) / p.LHV;

    return { states, w_c, w_t, q_in, w_net, m_fuel, T_exhaust: T5 };
}


/* ── HTC Steam (Rankine) Cycle ── */
function calculateSteamCycle(p) {
    const states = [];
    const Pa = p.P_boiler, Ta = p.T_steam, Pb = p.P_cond;

    // State a: Superheated steam leaving Boiler
    const ha = SteamTable.h_super(Pa, Ta);
    const sa = SteamTable.s_super(Pa, Ta);
    states.push({ id: 'ⓐ', label: 'Boiler Exit (Superheated)', T: Ta, P: Pa, h: ha, s: sa, x: '-' });

    // State b: After Steam Turbine expansion
    const sb_s = sa;
    const sg_c = SteamTable.sg(Pb), sf_c = SteamTable.sf(Pb);
    let hb_s, xb_s, Tb;
    if (sb_s < sg_c) {
        xb_s = SteamTable.x_from_s(Pb, sb_s);
        hb_s = SteamTable.h_from_x(Pb, xb_s);
        Tb = SteamTable.T_sat(Pb);
    } else {
        xb_s = 1;
        Tb = SteamTable.T_from_s_super(Pb, sb_s);
        hb_s = SteamTable.h_super(Pb, Tb);
    }
    const hb = ha - p.eta_st * (ha - hb_s);
    let xb = (hb - SteamTable.hf(Pb)) / SteamTable.hfg(Pb);
    let sb;
    if (xb <= 1) {
        sb = sf_c + xb * (sg_c - sf_c);
        Tb = SteamTable.T_sat(Pb);
    } else {
        sb = SteamTable.s_super(Pb, Tb);
        xb = '-';
    }
    states.push({ id: 'ⓑ', label: 'ST Exit', T: Tb, P: Pb, h: hb, s: sb, x: typeof xb === 'number' ? xb.toFixed(3) : xb });

    // State c: Condenser exit (saturated liquid)
    const Pc = Pb, Tc = SteamTable.T_sat(Pc);
    const hc = SteamTable.hf(Pc), sc = SteamTable.sf(Pc);
    states.push({ id: 'ⓒ', label: 'Condenser Exit (Sat. Liq.)', T: Tc, P: Pc, h: hc, s: sc, x: '0.000' });

    // State d: After Feed Pump
    const Pd = Pa;
    const wp_s = SteamTable.vf() * (Pd - Pc);
    const wp = wp_s / p.eta_fp;
    const hd = hc + wp;
    const sd = sc + 0.001;
    const Td = Tc + wp / 4.18;
    states.push({ id: 'ⓓ', label: 'Feed Pump Exit', T: Td, P: Pd, h: hd, s: sd, x: '-' });

    const w_st = ha - hb;
    const w_fp = wp;
    const q_boiler = ha - hd;
    const q_cond = hb - hc;
    const w_net = w_st - w_fp;

    return { states, w_st, w_fp, q_boiler, q_cond, w_net };
}


/* ── AD-HTC Mass & Energy Balance ── */
function calculateADHTC(p) {
    const m_total = p.m_biomass;
    const m_rich = m_total * p.moisture_split;
    const m_lean = m_total * (1 - p.moisture_split);
    const biogas_vol = m_rich * p.ad_yield;         // m³/s
    const biogas_density = 1.15;                     // kg/m³ (approx)
    const m_biogas = biogas_vol * biogas_density;    // kg/s
    const htc_energy = m_lean * 1.5 * (p.htc_temp - 298); // kJ/s (approx)

    return { m_total, m_rich, m_lean, biogas_vol, m_biogas, htc_energy };
}


// ────────────────────────────────────────────────
//  4. CHART RENDERING
// ────────────────────────────────────────────────

let hsChart = null, thChart = null;

function interpPts(s1, h1, s2, h2, n = 25) {
    const pts = [];
    for (let i = 0; i <= n; i++) {
        const t = i / n;
        pts.push({ x: s1 + t * (s2 - s1), y: h1 + t * (h2 - h1) });
    }
    return pts;
}

function renderHSChart(steamStates) {
    const ctx = document.getElementById('hs-chart').getContext('2d');
    if (hsChart) hsChart.destroy();

    const sm = {};
    steamStates.forEach(s => { sm[s.id] = s; });

    const processes = [
        { name: 'Boiler (ⓓ→ⓐ, Heat addition)', from: 'ⓓ', to: 'ⓐ', color: '#f97316' },
        { name: 'Steam Turbine (ⓐ→ⓑ, Expansion)', from: 'ⓐ', to: 'ⓑ', color: '#a78bfa' },
        { name: 'Condenser (ⓑ→ⓒ, Heat rejection)', from: 'ⓑ', to: 'ⓒ', color: '#22d3ee' },
        { name: 'Feed Pump (ⓒ→ⓓ, Compression)', from: 'ⓒ', to: 'ⓓ', color: '#3b82f6' },
    ];

    const datasets = processes.map(p => ({
        label: p.name,
        data: interpPts(sm[p.from].s, sm[p.from].h, sm[p.to].s, sm[p.to].h),
        borderColor: p.color, backgroundColor: p.color + '20',
        borderWidth: 3, pointRadius: 0, fill: false, showLine: true, tension: 0.3,
    }));

    // State markers
    const pts = ['ⓓ', 'ⓐ', 'ⓑ', 'ⓒ'].map(id => ({ x: sm[id].s, y: sm[id].h }));
    datasets.push({
        label: 'State Points', data: pts,
        borderColor: '#f1f5f9', backgroundColor: '#38bdf8',
        pointRadius: 8, pointHoverRadius: 11, pointBorderWidth: 2,
        showLine: false, pointStyle: 'circle',
    });

    // Saturation dome
    const dome = SteamTable.saturationDome(50);
    const liq = dome.map(d => ({ x: d.sf, y: d.hf }));
    const vap = dome.map(d => ({ x: d.sg, y: d.hg })).reverse();
    datasets.unshift({
        label: 'Saturation Dome', data: [...liq, ...vap, liq[0]],
        borderColor: 'rgba(148,163,184,0.35)', backgroundColor: 'rgba(148,163,184,0.04)',
        borderWidth: 1.5, borderDash: [6, 3], pointRadius: 0, fill: true, showLine: true, tension: 0.4,
    });

    hsChart = new Chart(ctx, {
        type: 'scatter', data: { datasets },
        options: {
            responsive: true, maintainAspectRatio: false,
            plugins: {
                legend: { position: 'top', labels: { color: '#94a3b8', font: { family: 'Inter', size: 11 }, usePointStyle: true, padding: 15 } },
                tooltip: {
                    backgroundColor: 'rgba(10,10,30,0.95)', borderColor: '#38bdf8', borderWidth: 1,
                    titleFont: { family: 'JetBrains Mono', size: 12 }, bodyFont: { family: 'JetBrains Mono', size: 11 },
                    callbacks: { label: ctx => `h=${ctx.parsed.y.toFixed(1)} kJ/kg, s=${ctx.parsed.x.toFixed(4)} kJ/(kg·K)` }
                },
            },
            scales: {
                x: {
                    title: { display: true, text: 'Entropy s [kJ/(kg·K)]', color: '#94a3b8', font: { family: 'Inter', size: 13, weight: 600 } },
                    grid: { color: 'rgba(56,189,248,0.06)' }, ticks: { color: '#64748b', font: { family: 'JetBrains Mono', size: 10 } }
                },
                y: {
                    title: { display: true, text: 'Enthalpy h [kJ/kg]', color: '#94a3b8', font: { family: 'Inter', size: 13, weight: 600 } },
                    grid: { color: 'rgba(56,189,248,0.06)' }, ticks: { color: '#64748b', font: { family: 'JetBrains Mono', size: 10 } }
                },
            }
        }
    });
}


function renderTHdotChart(gasStates) {
    const ctx = document.getElementById('th-chart').getContext('2d');
    if (thChart) thChart.destroy();

    // Gas cycle processes: 1→2 Compression, 2→4 Combustion (skip ③ duplicate), 4→5 Expansion
    const procs = [
        { name: 'Compression (①→②)', from: 0, to: 1, color: '#3b82f6' },
        { name: 'Biogas Combustion (③→④)', from: 2, to: 3, color: '#f97316' },
        { name: 'Turbine Expansion (④→⑤)', from: 3, to: 4, color: '#ef4444' },
    ];

    const datasets = [];
    let cumH = 0;

    procs.forEach(proc => {
        const s1 = gasStates[proc.from], s2 = gasStates[proc.to];
        const dh = Math.abs(s2.h - s1.h);
        const data = [];
        for (let i = 0; i <= 25; i++) {
            const t = i / 25;
            data.push({ x: cumH + t * dh, y: s1.T + t * (s2.T - s1.T) });
        }
        cumH += dh;
        datasets.push({
            label: proc.name, data,
            borderColor: proc.color, backgroundColor: proc.color + '15',
            borderWidth: 3, pointRadius: 0, fill: false, showLine: true, tension: 0.3,
        });
    });

    // State markers
    let mCum = 0;
    const markers = [{ x: 0, y: gasStates[0].T }];
    procs.forEach(proc => {
        const s1 = gasStates[proc.from], s2 = gasStates[proc.to];
        mCum += Math.abs(s2.h - s1.h);
        markers.push({ x: mCum, y: s2.T });
    });
    datasets.push({
        label: 'State Points', data: markers,
        borderColor: '#f1f5f9', backgroundColor: '#fb923c',
        pointRadius: 7, pointHoverRadius: 10, pointBorderWidth: 2,
        showLine: false, pointStyle: 'rectRot',
    });

    thChart = new Chart(ctx, {
        type: 'scatter', data: { datasets },
        options: {
            responsive: true, maintainAspectRatio: false,
            plugins: {
                legend: { position: 'top', labels: { color: '#94a3b8', font: { family: 'Inter', size: 11 }, usePointStyle: true, padding: 12 } },
                tooltip: {
                    backgroundColor: 'rgba(10,10,30,0.95)', borderColor: '#fb923c', borderWidth: 1,
                    titleFont: { family: 'JetBrains Mono', size: 12 }, bodyFont: { family: 'JetBrains Mono', size: 11 },
                    callbacks: { label: ctx => `T=${ctx.parsed.y.toFixed(1)} K, Ḣ=${ctx.parsed.x.toFixed(1)} kJ/kg` }
                },
            },
            scales: {
                x: {
                    title: { display: true, text: 'Cumulative Heat / Work Rate Ḣ [kJ/kg]', color: '#94a3b8', font: { family: 'Inter', size: 13, weight: 600 } },
                    grid: { color: 'rgba(251,146,60,0.06)' }, ticks: { color: '#64748b', font: { family: 'JetBrains Mono', size: 10 } }
                },
                y: {
                    title: { display: true, text: 'Temperature T [K]', color: '#94a3b8', font: { family: 'Inter', size: 13, weight: 600 } },
                    grid: { color: 'rgba(251,146,60,0.06)' }, ticks: { color: '#64748b', font: { family: 'JetBrains Mono', size: 10 } }
                },
            }
        }
    });
}


// ────────────────────────────────────────────────
//  5. REPORT POPULATION
// ────────────────────────────────────────────────

function populatePerformance(gas, steam, p) {
    const container = document.getElementById('perf-metrics');
    const eta_gas = ((gas.w_net / gas.q_in) * 100);
    const eta_steam = steam.q_boiler > 0 ? ((steam.w_net / steam.q_boiler) * 100) : 0;
    const bwr = (gas.w_c / gas.w_t * 100);
    const W_gas_total = gas.w_net * p.m_air;     // kW
    const W_steam_total = steam.w_net * 2;        // approx steam mass flow

    const metrics = [
        { value: gas.w_net.toFixed(1), label: 'Net Gas Cycle Work', unit: 'kJ/kg' },
        { value: steam.w_net.toFixed(1), label: 'Net Steam Cycle Work', unit: 'kJ/kg' },
        { value: eta_gas.toFixed(1), label: 'Gas Cycle Efficiency', unit: '%' },
        { value: eta_steam.toFixed(1), label: 'Steam Cycle Efficiency', unit: '%' },
        { value: bwr.toFixed(1), label: 'Back Work Ratio', unit: '%' },
        { value: gas.q_in.toFixed(1), label: 'Biogas Heat Input', unit: 'kJ/kg' },
        { value: (W_gas_total / 1000).toFixed(2), label: 'Gas Cycle Power', unit: 'MW' },
        { value: gas.states[4].T.toFixed(0), label: 'GT Exhaust Temp', unit: 'K' },
    ];
    container.innerHTML = metrics.map(m => `
        <div class="metric-item">
            <div class="metric-value">${m.value}</div>
            <div class="metric-unit">${m.unit}</div>
            <div class="metric-label">${m.label}</div>
        </div>
    `).join('');
}

function populateADMetrics(ad, gas, p) {
    const container = document.getElementById('ad-metrics');
    const metrics = [
        { value: ad.m_total.toFixed(2), label: 'Total Biomass Feed', unit: 'kg/s' },
        { value: ad.m_rich.toFixed(2), label: 'Moisture-Rich (→AD)', unit: 'kg/s' },
        { value: ad.m_lean.toFixed(2), label: 'Moisture-Lean (→HTC)', unit: 'kg/s' },
        { value: ad.biogas_vol.toFixed(2), label: 'Biogas Production', unit: 'm³/s' },
        { value: ad.m_biogas.toFixed(2), label: 'Biogas Mass Flow', unit: 'kg/s' },
        { value: gas.m_fuel.toFixed(2), label: 'Fuel Demand (Combustor)', unit: 'kg/s' },
        { value: (ad.htc_energy / 1000).toFixed(2), label: 'HTC Thermal Energy', unit: 'MW' },
        { value: p.htc_temp.toFixed(0), label: 'HTC Reactor Temp', unit: 'K' },
    ];
    container.innerHTML = metrics.map(m => `
        <div class="metric-item">
            <div class="metric-value">${m.value}</div>
            <div class="metric-unit">${m.unit}</div>
            <div class="metric-label">${m.label}</div>
        </div>
    `).join('');
}

function populateTable(tbodyId, states) {
    const tbody = document.getElementById(tbodyId);
    tbody.innerHTML = states.map(st => `
        <tr>
            <td style="color:#38bdf8;font-weight:600">${st.id}</td>
            <td>${st.label}</td>
            <td>${st.T.toFixed(1)}</td>
            <td>${st.P.toFixed(1)}</td>
            <td>${st.h.toFixed(1)}</td>
            <td>${st.s.toFixed(4)}</td>
            ${st.x !== undefined ? `<td>${st.x}</td>` : ''}
        </tr>
    `).join('');
}


// ────────────────────────────────────────────────
//  6. UI EVENTS
// ────────────────────────────────────────────────

const modal = document.getElementById('report-modal');
const analyseBtn = document.getElementById('analyse-btn');
const closeBtn = document.getElementById('modal-close');

analyseBtn.addEventListener('click', () => {
    const p = getParams();
    const gasResult = calculateGasCycle(p);
    const steamResult = calculateSteamCycle(p);
    const adResult = calculateADHTC(p);

    modal.classList.add('active');
    document.body.style.overflow = 'hidden';

    populatePerformance(gasResult, steamResult, p);
    populateADMetrics(adResult, gasResult, p);

    setTimeout(() => {
        renderHSChart(steamResult.states);
        renderTHdotChart(gasResult.states);
    }, 100);

    populateTable('gas-table-body', gasResult.states);
    populateTable('steam-table-body', steamResult.states);
});

closeBtn.addEventListener('click', () => { modal.classList.remove('active'); document.body.style.overflow = ''; });
modal.addEventListener('click', e => { if (e.target === modal) { modal.classList.remove('active'); document.body.style.overflow = ''; } });
document.addEventListener('keydown', e => { if (e.key === 'Escape' && modal.classList.contains('active')) { modal.classList.remove('active'); document.body.style.overflow = ''; } });


// ────────────────────────────────────────────────
//  7. SVG TOOLTIPS
// ────────────────────────────────────────────────

const tooltip = document.getElementById('svg-tooltip');
document.querySelectorAll('.component').forEach(comp => {
    comp.addEventListener('mouseenter', () => {
        tooltip.textContent = comp.getAttribute('data-tooltip');
        tooltip.classList.add('visible');
    });
    comp.addEventListener('mousemove', e => {
        tooltip.style.left = (e.clientX + 15) + 'px';
        tooltip.style.top = (e.clientY - 10) + 'px';
    });
    comp.addEventListener('mouseleave', () => tooltip.classList.remove('visible'));
});

console.log('AD-HTC Fuel-Enhanced Power Gas Cycle loaded.');
console.log('GasTable API:', Object.keys(GasTable));
console.log('SteamTable API:', Object.keys(SteamTable));
