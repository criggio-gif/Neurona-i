// snn-core.js — rete di neuroni spiking che impara a riconoscere cifre senza supervisione.
// Architettura di Diehl & Cook (2015), ridotta:
//   retina (256 pixel → treni di spike di Poisson)
//   → strato eccitatorio di neuroni integrate-and-fire a conduttanza (sinapsi AMPA-like)
//   ↔ inibizione laterale (interneuroni GABAergici impliciti: chi scarica zittisce gli altri)
// Apprendimento: STDP sulle sinapsi retina→strato, soglia adattativa (omeostasi intrinseca),
// normalizzazione dei pesi (scaling sinaptico). Nessuna etichetta durante l'apprendimento:
// le etichette servono solo dopo, per dare un nome a ciò che ogni neurone ha imparato.
(function (root) {
  'use strict';

  function createNetwork(opts) {
    opts = opts || {};
    const rand = opts.random || Math.random;
    const nIn = opts.nIn || 256;
    const nE = opts.nE || 49;
    const P = Object.assign({
      dt: 0.5,          // ms
      tPresent: 350,    // ms di presentazione per immagine
      maxRate: 64,      // Hz per un pixel pieno
      vRest: -65, vReset: -65, vThresh: -52, tRef: 5, tauM: 100,
      eExc: 0, eInh: -100, tauGe: 1, tauGi: 2,
      wInh: 8,          // forza dell'inibizione laterale
      thetaPlus: 0.1,   // mV aggiunti alla soglia a ogni spike (omeostasi)
      tauTheta: 1e6,    // ms: lento ritorno della soglia verso il valore di base
      tauPre: 20, tauPost: 20,
      etaPre: 0.0002, etaPost: 0.01, xTar: 0.4, wMax: 1,
      wSum: 0.1 * nIn,  // somma dei pesi in ingresso a ogni neurone
      minSpikes: 5,
    }, opts.params || {});

    const W = new Float32Array(nE * nIn);
    for (let k = 0; k < W.length; k++) W[k] = 0.3 * rand();
    const v = new Float32Array(nE), ge = new Float32Array(nE), gi = new Float32Array(nE);
    const theta = new Float32Array(nE), refr = new Float32Array(nE);
    const xPre = new Float32Array(nIn), xPost = new Float32Array(nE);
    const labels = new Int8Array(nE).fill(-1);
    const resp = new Float32Array(nE * 10), seen = new Float32Array(10);
    normalize();

    function normalize() {
      for (let j = 0; j < nE; j++) {
        let s = 0; const o = j * nIn;
        for (let i = 0; i < nIn; i++) s += W[o + i];
        const f = P.wSum / (s || 1);
        for (let i = 0; i < nIn; i++) W[o + i] = Math.min(P.wMax, W[o + i] * f);
      }
    }

    // Presenta un'immagine. Restituisce gli spike di ogni neurone.
    // onSpike(tipo, indice, tempo) permette all'interfaccia di disegnare il raster.
    function present(img, o) {
      o = o || {};
      const learn = !!o.learn;
      let rateBoost = 0, counts;
      for (let attempt = 0; attempt < 5; attempt++) {
        counts = run(img, learn, P.maxRate + rateBoost, o.onSpike);
        let tot = 0; for (let j = 0; j < nE; j++) tot += counts[j];
        if (tot >= P.minSpikes) break;
        rateBoost += 32; // come in Diehl & Cook: se la rete tace, si alza il contrasto
      }
      if (learn) {
        normalize();
        const f = Math.exp(-P.tPresent / P.tauTheta);
        for (let j = 0; j < nE; j++) theta[j] *= f;
      }
      return counts;
    }

    function run(img, learn, maxRate, onSpike) {
      const dt = P.dt, steps = Math.round(P.tPresent / dt);
      v.fill(P.vRest); ge.fill(0); gi.fill(0); refr.fill(0); xPre.fill(0); xPost.fill(0);
      const counts = new Uint16Array(nE);
      const pIn = new Float32Array(nIn);
      for (let i = 0; i < nIn; i++) pIn[i] = img[i] * maxRate * dt / 1000;
      const dGe = Math.exp(-dt / P.tauGe), dGi = Math.exp(-dt / P.tauGi);
      const dPre = Math.exp(-dt / P.tauPre), dPost = Math.exp(-dt / P.tauPost);
      for (let s = 0; s < steps; s++) {
        const t = s * dt;
        if (learn) { for (let i = 0; i < nIn; i++) xPre[i] *= dPre; for (let j = 0; j < nE; j++) xPost[j] *= dPost; }
        // spike in ingresso
        for (let i = 0; i < nIn; i++) {
          if (pIn[i] === 0 || rand() >= pIn[i]) continue;
          if (onSpike) onSpike('in', i, t);
          for (let j = 0, o = i; j < nE; j++, o += nIn) {
            ge[j] += W[o];
            if (learn) { const w = W[o] - P.etaPre * xPost[j]; W[o] = w < 0 ? 0 : w; } // pre dopo post → depressione
          }
          if (learn) xPre[i] = 1;
        }
        // neuroni eccitatori
        for (let j = 0; j < nE; j++) {
          ge[j] *= dGe; gi[j] *= dGi;
          if (refr[j] > 0) { refr[j] -= dt; continue; }
          v[j] += dt * ((P.vRest - v[j]) + ge[j] * (P.eExc - v[j]) + gi[j] * (P.eInh - v[j])) / P.tauM;
          if (v[j] > P.vThresh + theta[j]) {
            v[j] = P.vReset; refr[j] = P.tRef; counts[j]++;
            if (onSpike) onSpike('exc', j, t);
            for (let k = 0; k < nE; k++) if (k !== j) gi[k] += P.wInh; // inibizione laterale
            if (learn) {
              theta[j] += P.thetaPlus;
              const o = j * nIn;
              for (let i = 0; i < nIn; i++) { // pre prima di post → potenziamento (dipendente dal peso)
                const w = W[o + i] + P.etaPost * (xPre[i] - P.xTar) * (P.wMax - W[o + i]);
                W[o + i] = w < 0 ? 0 : (w > P.wMax ? P.wMax : w);
              }
              xPost[j] = 1;
            }
          }
        }
      }
      return counts;
    }

    // Etichettatura: ogni neurone prende il nome della cifra a cui risponde di più (in media).
    function observe(counts, label, decay) {
      decay = decay == null ? 0.995 : decay;
      for (let k = 0; k < resp.length; k++) resp[k] *= decay;
      for (let c = 0; c < 10; c++) seen[c] *= decay;
      for (let j = 0; j < nE; j++) resp[j * 10 + label] += counts[j];
      seen[label] += 1;
    }
    function assignLabels() {
      for (let j = 0; j < nE; j++) {
        let best = -1, bv = 0;
        for (let c = 0; c < 10; c++) { const r = seen[c] ? resp[j * 10 + c] / seen[c] : 0; if (r > bv) { bv = r; best = c; } }
        labels[j] = best;
      }
    }
    // Voto: per ogni cifra, media degli spike dei neuroni che la rappresentano.
    function classify(counts) {
      const sum = new Float32Array(10), n = new Float32Array(10);
      for (let j = 0; j < nE; j++) if (labels[j] >= 0) { sum[labels[j]] += counts[j]; n[labels[j]]++; }
      const score = Array.from(sum, (s, c) => (n[c] ? s / n[c] : 0));
      let best = -1, bv = 0;
      for (let c = 0; c < 10; c++) if (score[c] > bv) { bv = score[c]; best = c; }
      return { digit: best, score };
    }

    function exportState() { return { nIn, nE, P, W: Array.from(W), theta: Array.from(theta), labels: Array.from(labels) }; }

    return { nIn, nE, P, W, theta, labels, present, observe, assignLabels, classify, exportState };
  }

  root.SNN = { createNetwork };
  if (typeof module !== 'undefined') module.exports = root.SNN;
})(typeof window !== 'undefined' ? window : globalThis);
