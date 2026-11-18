// digits.js — generatore di cifre "scritte a mano" 16×16 e normalizzazione delle immagini.
// Ogni cifra è un insieme di tratti (polilinee nel quadrato unitario) deformati a caso:
// rotazione, scala, inclinazione, spostamento, tremolio dei punti e spessore del tratto.
(function (root) {
  'use strict';
  const S = 16;

  function ellipse(cx, cy, rx, ry, a0, a1, n) {
    const pts = [];
    for (let k = 0; k <= n; k++) { const a = a0 + (a1 - a0) * k / n; pts.push([cx + rx * Math.cos(a), cy + ry * Math.sin(a)]); }
    return pts;
  }
  const TAU = Math.PI * 2;
  const STROKES = [
    [ellipse(0.5, 0.5, 0.27, 0.4, 0, TAU, 24)],
    [[[0.36, 0.24], [0.52, 0.1], [0.52, 0.9]]],
    [[[0.25, 0.3], [0.33, 0.15], [0.52, 0.09], [0.7, 0.18], [0.72, 0.36], [0.5, 0.6], [0.24, 0.9], [0.78, 0.9]]],
    [[[0.25, 0.15], [0.55, 0.09], [0.73, 0.22], [0.7, 0.38], [0.46, 0.48], [0.72, 0.6], [0.75, 0.78], [0.55, 0.91], [0.24, 0.85]]],
    [[[0.66, 0.9], [0.66, 0.1], [0.22, 0.64], [0.8, 0.64]]],
    [[[0.76, 0.1], [0.32, 0.1], [0.28, 0.46], [0.55, 0.4], [0.74, 0.55], [0.73, 0.78], [0.55, 0.91], [0.25, 0.85]]],
    [[[0.7, 0.11], [0.45, 0.16], [0.3, 0.4], [0.27, 0.7], [0.4, 0.9], [0.62, 0.89], [0.74, 0.7], [0.62, 0.52], [0.4, 0.52], [0.28, 0.66]]],
    [[[0.22, 0.1], [0.78, 0.1], [0.42, 0.9]]],
    [ellipse(0.5, 0.29, 0.2, 0.19, 0, TAU, 18), ellipse(0.5, 0.7, 0.25, 0.21, 0, TAU, 20)],
    [ellipse(0.5, 0.32, 0.22, 0.21, 0, TAU, 18), [[0.72, 0.32], [0.69, 0.6], [0.56, 0.91]]],
  ];

  function segDist(px, py, ax, ay, bx, by) {
    const dx = bx - ax, dy = by - ay, l2 = dx * dx + dy * dy || 1e-9;
    const f = Math.max(0, Math.min(1, ((px - ax) * dx + (py - ay) * dy) / l2));
    return Math.hypot(px - ax - f * dx, py - ay - f * dy);
  }

  function render(digit, rand) {
    rand = rand || Math.random;
    const r = (a, b) => a + (b - a) * rand();
    const rot = r(-0.22, 0.22), sx = r(0.8, 1.12), sy = r(0.85, 1.1), sh = r(-0.25, 0.25), tx = r(-0.05, 0.05), ty = r(-0.05, 0.05);
    const c = Math.cos(rot), s = Math.sin(rot), jit = 0.03;
    const strokes = STROKES[digit].map(pl => pl.map(([x, y]) => {
      let X = (x - 0.5) * sx + sh * (y - 0.5), Y = (y - 0.5) * sy;
      X += r(-jit, jit); Y += r(-jit, jit);
      return [(c * X - s * Y + 0.5 + tx) * S, (s * X + c * Y + 0.5 + ty) * S];
    }));
    const th = r(0.55, 1.05); // semispessore in pixel
    const img = new Float32Array(S * S);
    for (let y = 0; y < S; y++) for (let x = 0; x < S; x++) {
      let d = 1e9;
      for (const pl of strokes) for (let k = 1; k < pl.length; k++) d = Math.min(d, segDist(x + 0.5, y + 0.5, pl[k - 1][0], pl[k - 1][1], pl[k][0], pl[k][1]));
      img[y * S + x] = Math.max(0, Math.min(1, th + 0.5 - d));
    }
    return normalize(img, S);
  }

  // Ritaglia l'inchiostro e lo riscala in un riquadro di 12 px centrato (stessa procedura per
  // cifre generate e disegnate), come si fa con MNIST.
  function normalize(src, n) {
    let x0 = n, y0 = n, x1 = -1, y1 = -1, mass = 0;
    for (let y = 0; y < n; y++) for (let x = 0; x < n; x++) if (src[y * n + x] > 0.1) { x0 = Math.min(x0, x); y0 = Math.min(y0, y); x1 = Math.max(x1, x); y1 = Math.max(y1, y); mass++; }
    const out = new Float32Array(S * S);
    if (!mass) return out;
    const bw = x1 - x0 + 1, bh = y1 - y0 + 1, box = 12 * n / S;
    const scale = Math.max(bw, bh) / box; // pixel sorgente per pixel di destinazione (in unità di n)
    const cx = (x0 + x1 + 1) / 2, cy = (y0 + y1 + 1) / 2, k = n / S;
    for (let y = 0; y < S; y++) for (let x = 0; x < S; x++) {
      // media su una piccola finestra per evitare aliasing quando si rimpicciolisce
      let acc = 0, cnt = 0;
      for (let sy = 0; sy < 2; sy++) for (let sx = 0; sx < 2; sx++) {
        const fx = cx + ((x + 0.25 + sx * 0.5) - S / 2) * k * scale - 0.5;
        const fy = cy + ((y + 0.25 + sy * 0.5) - S / 2) * k * scale - 0.5;
        acc += bilinear(src, n, fx, fy); cnt++;
      }
      out[y * S + x] = Math.min(1, acc / cnt);
    }
    return out;
  }
  function bilinear(a, n, x, y) {
    const xi = Math.floor(x), yi = Math.floor(y), fx = x - xi, fy = y - yi;
    const g = (X, Y) => (X < 0 || Y < 0 || X >= n || Y >= n) ? 0 : a[Y * n + X];
    return g(xi, yi) * (1 - fx) * (1 - fy) + g(xi + 1, yi) * fx * (1 - fy) + g(xi, yi + 1) * (1 - fx) * fy + g(xi + 1, yi + 1) * fx * fy;
  }

  root.Digits = { S, render, normalize };
  if (typeof module !== 'undefined') module.exports = root.Digits;
})(typeof window !== 'undefined' ? window : globalThis);
