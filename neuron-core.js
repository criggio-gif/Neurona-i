// neuron-core.js — modello biofisico multicompartimentale di un neurone piramidale.
// Unità: V in mV, t in ms, capacità in pF, conduttanze in nS, correnti in pA, [Ca] in µM.
// Il voltaggio è integrato con Eulero implicito sull'albero dei compartimenti (algoritmo di Hines),
// le variabili di gating con Eulero esponenziale. Nessuna dipendenza dal DOM: gira anche in Node.
(function (root) {
  'use strict';

  const DT = 0.025; // ms
  const E = { Na: 50, K: -90, Ca: 120, h: -45, AMPA: 0, NMDA: 0, GABAA: -75, GABAB: -95 };
  const RA = 100;       // resistività assiale, Ω·cm
  const VT0 = -56.2;    // soglia di riferimento dei canali Na/K (Pospischil et al. 2008)
  const CA_REST = 0.05; // µM
  const TAU_CA = 80;    // ms, estrusione/buffering del calcio

  function efun(u) { return Math.abs(u) < 1e-6 ? 1 - u / 2 : u / (Math.exp(u) - 1); }

  // ---------------------------------------------------------------------------------------------
  // Morfologia: ogni compartimento è un cilindro (d, L in µm) con coordinate di disegno (x0,y0)->(x1,y1)
  // nello spazio logico 600×760. L'ordine garantisce parent[i] < i (necessario per Hines).
  function buildMorphology() {
    const comps = [];
    function add(c) { c.i = comps.length; comps.push(c); return c.i; }
    function branch(parent, n, x0, y0, x1, y1, L, d0, d1, type, name) {
      let p = parent;
      for (let k = 0; k < n; k++) {
        const f0 = k / n, f1 = (k + 1) / n;
        const d = d0 + (d1 - d0) * (k + 0.5) / n;
        const dist = (p >= 0 ? comps[p].dist : 0) + L;
        p = add({
          parent: p, type, name, L, d, dist, seg: k, nseg: n,
          x0: x0 + (x1 - x0) * f0, y0: y0 + (y1 - y0) * f0,
          x1: x0 + (x1 - x0) * f1, y1: y0 + (y1 - y0) * f1,
        });
      }
      return p;
    }
    const soma = add({ parent: -1, type: 'soma', name: 'Soma', L: 20, d: 20, dist: 0, x0: 300, y0: 348, x1: 300, y1: 312, seg: 0, nseg: 1 });
    const trunk = branch(soma, 4, 300, 312, 300, 172, 60, 3.5, 2.4, 'apic', 'Tronco apicale');
    // obliquo dal secondo compartimento del tronco
    branch(2, 3, 300, 242, 398, 214, 45, 1.2, 0.9, 'apic', 'Dendrite obliquo');
    branch(trunk, 4, 300, 172, 200, 52, 50, 1.6, 0.9, 'apic', 'Ciuffo apicale sinistro');
    branch(trunk, 4, 300, 172, 410, 58, 50, 1.6, 0.9, 'apic', 'Ciuffo apicale destro');
    branch(soma, 4, 288, 344, 170, 428, 40, 1.6, 0.8, 'dend', 'Dendrite basale sinistro');
    branch(soma, 4, 312, 344, 430, 428, 40, 1.6, 0.8, 'dend', 'Dendrite basale destro');
    let ax = branch(soma, 2, 300, 348, 300, 392, 25, 1.6, 1.3, 'ais', 'Segmento iniziale (AIS)');
    let y = 392;
    for (let k = 0; k < 5; k++) {
      ax = branch(ax, 1, 300, y, 300, y + 46, 150, 1.0, 1.0, 'internode', 'Internodo mielinizzato ' + (k + 1));
      y += 46;
      ax = branch(ax, 1, 300, y, 300, y + 5, 1, 1.0, 1.0, 'node', 'Nodo di Ranvier ' + (k + 1));
      y += 5;
    }
    branch(ax, 1, 300, y, 300, y + 16, 3, 2.5, 2.5, 'term', 'Bottone presinaptico');
    return comps;
  }

  // Densità dei canali in mS/cm² (e cm in µF/cm²) per tipo di compartimento.
  function densities(c, myelin) {
    const z = { Na: 0, K: 0, M: 0, Ca: 0, H: 0, SK: 0, L: 0.1, cm: 1, vt: VT0 };
    switch (c.type) {
      case 'soma': Object.assign(z, { Na: 50, K: 6, M: 0.08, Ca: 0.2, H: 0.02, SK: 0.15 }); break;
      case 'dend': Object.assign(z, { Na: 5, K: 3, M: 0.02, Ca: 0.05, H: 0.03, SK: 0.1 }); break;
      case 'apic': { // gradiente di HCN lungo il dendrite apicale (Magee 1998)
        const f = Math.min(1, c.dist / 350);
        Object.assign(z, { Na: 6, K: 3, M: 0.02, Ca: 0.06, H: 0.03 + 0.2 * f, SK: 0.1 });
        break;
      }
      case 'ais': Object.assign(z, { Na: 300, K: 150, M: 5, vt: VT0 - 6 }); break;
      case 'internode':
        // Kv1 juxtaparanodali sotto la mielina stabilizzano l'internodo dopo ogni spike
        if (myelin) Object.assign(z, { L: 0.01, K: 2, cm: 0.04 });
        else Object.assign(z, { Na: 0, K: 0, L: 0.1, cm: 1 });
        break;
      case 'node': Object.assign(z, { Na: 2000, K: 800, M: 1, L: 5 }); break;
      case 'term': Object.assign(z, { Na: 200, K: 120, Ca: 0.15, SK: 0 }); break;
    }
    return z;
  }

  // ---------------------------------------------------------------------------------------------
  function createNeuron(opts) {
    opts = opts || {};
    let rand = opts.random || Math.random;
    const comps = buildMorphology();
    const N = comps.length;
    const parent = new Int32Array(N);
    const area = new Float64Array(N);
    const C = new Float64Array(N);
    const gax = new Float64Array(N); // conduttanza verso il genitore (nS)
    const g = { Na: new Float64Array(N), K: new Float64Array(N), M: new Float64Array(N), Ca: new Float64Array(N),
      H: new Float64Array(N), SK: new Float64Array(N), L: new Float64Array(N) };
    const vt = new Float64Array(N);
    const depth = new Float64Array(N), tauCa = new Float64Array(N); // spessore dello shell sottomembrana per il Ca²⁺ (µm)
    const EL = -70;

    const V = new Float64Array(N).fill(-70);
    const m = new Float64Array(N), h = new Float64Array(N), n = new Float64Array(N), p = new Float64Array(N);
    const q = new Float64Array(N), r = new Float64Array(N), hc = new Float64Array(N), zk = new Float64Array(N);
    const ca = new Float64Array(N).fill(CA_REST);
    const iinj = new Float64Array(N);
    const Iion = { Na: new Float64Array(N), K: new Float64Array(N), M: new Float64Array(N), Ca: new Float64Array(N),
      H: new Float64Array(N), SK: new Float64Array(N), L: new Float64Array(N), syn: new Float64Array(N) };

    const dg = new Float64Array(N), bb = new Float64Array(N), Gs = new Float64Array(N), GEs = new Float64Array(N);
    const caIn = new Float64Array(N); // influsso extra di Ca (es. da NMDA), µM/ms

    const params = {
      temperature: 36, mg: 1.0, myelin: true,
      block: { Na: 1, K: 1, M: 1, Ca: 1, H: 1, SK: 1, AMPA: 1, NMDA: 1, GABAA: 1, GABAB: 1 },
      excRate: 8, inhRate: 8, // Hz per sinapsi
    };

    function setupPassive() {
      for (let i = 0; i < N; i++) {
        const c = comps[i];
        parent[i] = c.parent;
        area[i] = Math.PI * c.d * c.L;
        const z = densities(c, params.myelin);
        c.dens = z;
        C[i] = z.cm * area[i] * 1e-2;
        g.Na[i] = z.Na * area[i] * 1e-2; g.K[i] = z.K * area[i] * 1e-2; g.M[i] = z.M * area[i] * 1e-2;
        g.Ca[i] = z.Ca * area[i] * 1e-2; g.H[i] = z.H * area[i] * 1e-2; g.SK[i] = z.SK * area[i] * 1e-2;
        g.L[i] = z.L * area[i] * 1e-2;
        vt[i] = z.vt;
        depth[i] = c.type === 'term' ? 0.1 : 1.0;
        tauCa[i] = c.type === 'term' ? 6 : TAU_CA; // nel bottone il Ca²⁺ del nanodominio si dissipa in pochi ms
      }
      for (let i = 1; i < N; i++) {
        const a = comps[i], b = comps[parent[i]];
        const Ra_ = RA * (a.L / 2) * 4e4 / (Math.PI * a.d * a.d);
        const Rb = RA * (b.L / 2) * 4e4 / (Math.PI * b.d * b.d);
        gax[i] = 1e9 / (Ra_ + Rb);
      }
    }
    setupPassive();

    // ---- Sinapsi in ingresso -------------------------------------------------------------------
    const synapses = [];
    // Ogni punto sinaptico rappresenta un fascio di afferenze indipendenti (aff), ciascuna con più
    // contatti (contacts); ogni contatto ha una vescicola pronta al rilascio con probabilità U e si
    // ricarica con costante trec (depressione a breve termine). Conduttanze in nS per quanto.
    const EXC = { gA: 0.25, gN: 0.15, U: 0.5, trec: 450, contacts: 6, aff: 20 };
    const INH = { gA: 1.0, gB: 1.5, U: 0.35, trec: 250, contacts: 5, aff: 10 };
    const kA = 1 - Math.exp(-1.1 * 1 * 1), kN = 1 - Math.exp(-0.072 * 1 * 1), kGA = 1 - Math.exp(-5 * 1 * 1), kGB = 1 - Math.exp(-0.09 * 0.5 * 0.3); // GABA-B: impulso 0.5 mM × 0.3 ms (Destexhe 1994)
    const dA = Math.exp(-0.19 * DT), dN = Math.exp(-0.0066 * DT), dGA = Math.exp(-0.18 * DT), dGB = Math.exp(-0.0012 * DT);
    let side = 1;
    function addSyn(ci, kind, frac) {
      side = -side;
      const P = kind === 'exc' ? EXC : INH;
      const S = P.contacts * P.aff;
      synapses.push({ comp: ci, kind, frac: frac == null ? 0.5 : frac, side, sites: S, aff: P.aff,
        rA: 0, rN: 0, rGA: 0, rGB: 0, sGB: 0, docked: S, lastRel: -1e9, lastSpike: -1e9, lastQ: 0, gNow: 0 });
    }
    for (let i = 0; i < N; i++) {
      const c = comps[i];
      if (c.type === 'apic' || c.type === 'dend') {
        const nS = (c.type === 'apic' && c.dist > 150) || (c.type === 'dend' && c.seg > 0) ? 2 : 1;
        for (let k = 0; k < nS; k++) addSyn(i, 'exc', (k + 0.5) / nS);
        if (c.seg === 0 && c.name !== 'Tronco apicale') addSyn(i, 'inh', 0.3);
      }
    }
    for (let k = 0; k < 4; k++) addSyn(0, 'inh', 0.2 + 0.2 * k);
    comps.forEach((c, i) => { if (c.type === 'ais') addSyn(i, 'inh', 0.5); });

    // ---- Terminale presinaptico in uscita --------------------------------------------------------
    const term = comps.findIndex(c => c.type === 'term');
    const aisEnd = comps.reduce((a, c, i) => (c.type === 'ais' ? i : a), -1);
    const NVES = 12;
    const vesicles = new Uint8Array(NVES).fill(1);
    const releases = []; // tempi di rilascio (ms) recenti
    let totalReleased = 0;

    // ---- Stato globale, spike detection -----------------------------------------------------------
    let t = 0;
    const det = { soma: { above: false, times: [] }, ais: { above: false, times: [] }, term: { above: false, times: [] } };
    let lastDelay = null;

    function steadyState() {
      phi = Math.pow(3, (params.temperature - 36) / 10);
      for (let i = 0; i < N; i++) gatesStep(i, V[i], 1e9);
    }

    let phi = 1;
    function gatesStep(i, v, dt) {
      const u = v - vt[i];
      // Na (m³h) e K delayed rectifier (n⁴): Traub/Pospischil
      let a = 1.28 * efun(-(u - 13) / 4), b = 1.4 * efun((u - 40) / 5);
      let tau = 1 / ((a + b) * phi); m[i] = a / (a + b) + (m[i] - a / (a + b)) * Math.exp(-dt / tau);
      a = 0.128 * Math.exp(-(u - 17) / 18); b = 4 / (1 + Math.exp(-(u - 40) / 5));
      tau = 1 / ((a + b) * phi); h[i] = a / (a + b) + (h[i] - a / (a + b)) * Math.exp(-dt / tau);
      a = 0.16 * efun(-(u - 15) / 5); b = 0.5 * Math.exp(-(u - 10) / 40);
      tau = 1 / ((a + b) * phi); n[i] = a / (a + b) + (n[i] - a / (a + b)) * Math.exp(-dt / tau);
      // M (Kv7), lenta, responsabile dell'adattamento
      let inf = 1 / (1 + Math.exp(-(v + 35) / 10));
      tau = 1000 / (3.3 * Math.exp((v + 35) / 20) + Math.exp(-(v + 35) / 20)) / phi;
      p[i] = inf + (p[i] - inf) * Math.exp(-dt / tau);
      // Ca alta soglia (L / P/Q nel terminale): q²r
      a = 0.209 * efun((-27 - v) / 3.8); b = 0.94 * Math.exp((-75 - v) / 17);
      tau = 1 / ((a + b) * phi); q[i] = a / (a + b) + (q[i] - a / (a + b)) * Math.exp(-dt / tau);
      a = 0.000457 * Math.exp((-13 - v) / 50); b = 0.0065 / (Math.exp((-15 - v) / 28) + 1);
      tau = 1 / ((a + b) * phi); r[i] = a / (a + b) + (r[i] - a / (a + b)) * Math.exp(-dt / tau);
      // HCN (Ih), attivata dall'iperpolarizzazione
      inf = 1 / (1 + Math.exp((v + 82) / 7));
      a = 0.00643 * (v + 154.9) / (Math.exp((v + 154.9) / 11.9) - 1); b = 0.193 * Math.exp(v / 33.1);
      tau = Math.min(1000, 1 / (a + b)) / phi; hc[i] = inf + (hc[i] - inf) * Math.exp(-dt / tau);
      // SK, attivata dal calcio intracellulare
      const c4 = Math.pow(ca[i], 4); inf = c4 / (c4 + 0.0625);
      zk[i] = inf + (zk[i] - inf) * Math.exp(-dt / (15 / phi));
    }

    function synStep() {
      const pe = params.excRate * DT / 1000, pi = params.inhRate * DT / 1000;
      for (let k = 0; k < synapses.length; k++) {
        const s = synapses[k];
        s.rA *= dA; s.rN *= dN; s.rGA *= dGA; s.rGB *= dGB;
        s.sGB += DT * (0.18 * s.rGB - 0.034 * s.sGB);
        const P = s.kind === 'exc' ? EXC : INH;
        s.docked += (s.sites - s.docked) * DT / P.trec;
        if (rand() < (s.kind === 'exc' ? pe : pi) * s.aff) presynSpike(s, 1);
      }
    }
    // nAff afferenti del fascio scaricano insieme; ognuna dispone in media di docked/aff vescicole pronte
    function presynSpike(s, nAff, pRel) {
      s.lastSpike = t;
      const P = s.kind === 'exc' ? EXC : INH;
      const avail = s.docked * nAff / s.aff;
      const whole = Math.floor(avail), n = whole + (rand() < avail - whole ? 1 : 0);
      const pr = pRel || P.U;
      let nq = 0;
      for (let j = 0; j < n; j++) if (rand() < pr) nq++;
      if (!nq) return 0;
      s.docked = Math.max(0, s.docked - nq); s.lastRel = t; s.lastQ = nq;
      const S = s.sites;
      // ogni quanto apre i recettori del proprio contatto; r ∈ [0, S] (saturazione per contatto)
      if (s.kind === 'exc') { s.rA += nq * kA * (1 - s.rA / S); s.rN += nq * kN * (1 - s.rN / S); }
      else { s.rGA += nq * kGA * (1 - s.rGA / S); s.rGB += nq * kGB * (1 - s.rGB / S); }
      return nq;
    }

    function step() {
      synStep();
      Gs.fill(0); GEs.fill(0); caIn.fill(0); Iion.syn.fill(0);
      const B = params.block;
      const mgf = params.mg / 3.57;
      for (let k = 0; k < synapses.length; k++) {
        const s = synapses[k], i = s.comp, v = V[i];
        let gsum = 0;
        if (s.kind === 'exc') {
          const gA = EXC.gA * s.rA * B.AMPA;
          const gN = EXC.gN * s.rN * B.NMDA / (1 + Math.exp(-0.062 * v) * mgf);
          Gs[i] += gA + gN; GEs[i] += gA * E.AMPA + gN * E.NMDA; gsum = gA + gN;
          Iion.syn[i] += gA * (v - E.AMPA) + gN * (v - E.NMDA);
          caIn[i] += -0.1 * gN * (v - E.NMDA); // ~10% della corrente NMDA è portata dal Ca²⁺ (pA)
        } else {
          const gA = INH.gA * s.rGA * B.GABAA;
          const s4 = Math.pow(s.sGB / s.sites * INH.contacts, 4);
          const gB = INH.gB * s.sites * B.GABAB * s4 / (s4 + 100);
          Gs[i] += gA + gB; GEs[i] += gA * E.GABAA + gB * E.GABAB; gsum = gA + gB;
          Iion.syn[i] += gA * (v - E.GABAA) + gB * (v - E.GABAB);
        }
        s.gNow = gsum;
      }
      const dtC = 1 / DT;
      phi = Math.pow(3, (params.temperature - 36) / 10);
      for (let i = 0; i < N; i++) {
        const v = V[i];
        gatesStep(i, v, DT);
        const gNa = g.Na[i] * B.Na * m[i] * m[i] * m[i] * h[i];
        const n2 = n[i] * n[i];
        const gK = g.K[i] * B.K * n2 * n2;
        const gM = g.M[i] * B.M * p[i];
        const gCa = g.Ca[i] * B.Ca * q[i] * q[i] * r[i];
        const gH = g.H[i] * B.H * hc[i];
        const gSK = g.SK[i] * B.SK * zk[i];
        const gL = g.L[i];
        Iion.Na[i] = gNa * (v - E.Na); Iion.K[i] = gK * (v - E.K); Iion.M[i] = gM * (v - E.K);
        Iion.Ca[i] = gCa * (v - E.Ca); Iion.H[i] = gH * (v - E.h); Iion.SK[i] = gSK * (v - E.K); Iion.L[i] = gL * (v - EL);
        const G = gNa + gK + gM + gCa + gH + gSK + gL + Gs[i];
        const GE = gNa * E.Na + (gK + gM + gSK) * E.K + gCa * E.Ca + gH * E.h + gL * EL + GEs[i];
        dg[i] = C[i] * dtC + G;
        bb[i] = C[i] * dtC * v + GE + iinj[i];
        // dinamica del calcio (shell sottomembrana di 0.1 µm)
        const iCaDens = (Iion.Ca[i] - caIn[i]) * 100 / area[i]; // µA/cm²
        const drive = Math.max(0, -0.0518 * iCaDens / depth[i]);
        ca[i] += DT * (drive - (ca[i] - CA_REST) / tauCa[i]);
      }
      // Hines
      for (let i = 1; i < N; i++) { dg[i] += gax[i]; dg[parent[i]] += gax[i]; }
      for (let i = N - 1; i >= 1; i--) {
        const pi = parent[i], f = -gax[i] / dg[i];
        dg[pi] -= f * -gax[i];
        bb[pi] -= f * bb[i];
      }
      V[0] = bb[0] / dg[0];
      for (let i = 1; i < N; i++) V[i] = (bb[i] + gax[i] * V[parent[i]]) / dg[i];
      t += DT;

      // rilascio vescicolare nel terminale: tasso ∝ [Ca]^4 (cooperatività del sensore sinaptotagmina)
      const c4 = Math.pow(ca[term], 4), rate = 1.0 * c4 / (c4 + Math.pow(5, 4));
      for (let k = 0; k < NVES; k++) {
        if (vesicles[k]) { if (rand() < rate * DT) { vesicles[k] = 0; releases.push(t); totalReleased++; } }
        else if (rand() < DT / 700) vesicles[k] = 1;
      }
      while (releases.length && releases[0] < t - 2000) releases.shift();

      detect(det.soma, V[0], 'soma');
      detect(det.ais, V[aisEnd], 'ais');
      detect(det.term, V[term], 'term');
    }

    function detect(d, v, which) {
      if (!d.above && v > -10) {
        d.above = true; d.times.push(t);
        if (d.times.length > 400) d.times.shift();
        if (which === 'term' && det.ais.times.length) {
          const ts = det.ais.times[det.ais.times.length - 1];
          if (t - ts < 5) lastDelay = t - ts;
        }
      } else if (d.above && v < -30) d.above = false;
    }

    function firingRate(which, win) {
      win = win || 1000;
      const arr = det[which].times; let c = 0;
      for (let k = arr.length - 1; k >= 0 && arr[k] > t - win; k--) c++;
      return c * 1000 / Math.min(win, Math.max(t, 1));
    }

    // Stimolo extracellulare: recluta in sincronia una frazione delle afferenze di ogni fascio selezionato
    function stimulateSynapses(filterFn, recruit) {
      let c = 0;
      for (const s of synapses) if (filterFn(s)) c += presynSpike(s, Math.max(1, Math.round(s.aff * (recruit || 0.12))));
      return c;
    }

    function setMyelin(on) {
      params.myelin = on;
      setupPassive();
    }

    function reset() {
      t = 0; V.fill(-70); ca.fill(CA_REST); vesicles.fill(1); releases.length = 0; totalReleased = 0;
      for (const k in det) { det[k].times.length = 0; det[k].above = false; } lastDelay = null;
      for (const s of synapses) Object.assign(s, { rA: 0, rN: 0, rGA: 0, rGB: 0, sGB: 0, docked: s.sites, lastRel: -1e9, lastSpike: -1e9 });
      steadyState();
      // 1 s di assestamento con gli ingressi attivi: canali lenti (M, HCN), depressione sinaptica e GABA-B
      // raggiungono il loro stato stazionario prima che la registrazione cominci.
      const saveInj = iinj.slice(); iinj.fill(0);
      for (let k = 0; k < 40000; k++) step();
      iinj.set(saveInj);
      t = 0; for (const k in det) { det[k].times.length = 0; } releases.length = 0; totalReleased = 0; lastDelay = null;
      for (const s of synapses) { s.lastRel = -1e9; s.lastSpike = -1e9; }
    }

    // Conduttanze sinaptiche totali (nS) per tipo di recettore
    function synConductances() {
      const out = { AMPA: 0, NMDA: 0, GABAA: 0, GABAB: 0 }, B = params.block, mgf = params.mg / 3.57;
      for (const s of synapses) {
        if (s.kind === 'exc') {
          out.AMPA += EXC.gA * s.rA * B.AMPA;
          out.NMDA += EXC.gN * s.rN * B.NMDA / (1 + Math.exp(-0.062 * V[s.comp]) * mgf);
        } else {
          out.GABAA += INH.gA * s.rGA * B.GABAA;
          const s4 = Math.pow(s.sGB / s.sites * INH.contacts, 4); out.GABAB += INH.gB * s.sites * B.GABAB * s4 / (s4 + 100);
        }
      }
      return out;
    }

    function channelState(i) {
      const B = params.block;
      return [
        { key: 'Na', label: 'Na⁺ voltaggio-dip. (Nav)', gbar: g.Na[i], open: m[i] ** 3 * h[i] * B.Na, I: Iion.Na[i] },
        { key: 'K', label: 'K⁺ rettificatore ritardato (Kv)', gbar: g.K[i], open: n[i] ** 4 * B.K, I: Iion.K[i] },
        { key: 'M', label: 'K⁺ M (Kv7)', gbar: g.M[i], open: p[i] * B.M, I: Iion.M[i] },
        { key: 'Ca', label: 'Ca²⁺ alta soglia (Cav)', gbar: g.Ca[i], open: q[i] * q[i] * r[i] * B.Ca, I: Iion.Ca[i] },
        { key: 'H', label: 'HCN (Ih)', gbar: g.H[i], open: hc[i] * B.H, I: Iion.H[i] },
        { key: 'SK', label: 'K⁺ attivato da Ca²⁺ (SK)', gbar: g.SK[i], open: zk[i] * B.SK, I: Iion.SK[i] },
        { key: 'L', label: 'Leak', gbar: g.L[i], open: 1, I: Iion.L[i] },
      ].filter(c => c.gbar > 0);
    }

    reset();

    return {
      DT, E, comps, N, synapses, params, V, ca, iinj, area, C, m, h, n, term, aisEnd, vesicles, releases,
      step, reset, setMyelin, stimulateSynapses, synConductances, channelState, firingRate,
      get t() { return t; }, get lastDelay() { return lastDelay; }, get totalReleased() { return totalReleased; },
      det, setRandom(fn) { rand = fn; },
    };
  }

  root.NeuronCore = { createNeuron, DT };
  if (typeof module !== 'undefined') module.exports = root.NeuronCore;
})(typeof window !== 'undefined' ? window : globalThis);
