# Neurona-i

Simulatore biofisico di un neurone piramidale corticale, che gira interamente nel browser.

Apri `index.html` (oppure servi la cartella con `python3 -m http.server`).

## Il modello (`neuron-core.js`)

- **Struttura**: 37 compartimenti cilindrici (soma, dendriti apicali e basali, segmento iniziale dell'assone,
  5 internodi mielinizzati e nodi di Ranvier, bottone presinaptico). Equazione del cavo risolta con
  Eulero implicito e algoritmo di Hines, dt = 25 µs.
- **Canali**: Nav, Kv (Traub/Pospischil 2008), Kv7/M, Cav ad alta soglia, HCN con gradiente apicale,
  SK Ca²⁺-dipendenti, leak. Densità specifiche per regione (AIS e nodi ricchi di Nav).
- **Recettori**: AMPA, NMDA (blocco Mg²⁺ di Jahr e Stevens), GABA-A, GABA-B con proteina G
  (cinetiche di Destexhe 1994). 53 fasci di afferenze con rilascio stocastico e depressione a breve termine.
- **Calcio e rilascio**: Ca²⁺ sottomembrana per compartimento; nel bottone 12 siti di rilascio con
  tasso ∝ [Ca²⁺]⁴.
- **Farmaci**: TTX, TEA, XE991, Cd²⁺, ZD7288, apamina, CNQX, D-AP5, bicucullina, CGP 55845.
  Si possono variare anche [Mg²⁺], la temperatura (Q10 = 3) e la demielinizzazione.

Il nucleo non dipende dal DOM e si usa anche in Node:

```js
const { createNeuron } = require('./neuron-core.js');
const n = createNeuron();
n.iinj[0] = 200;                 // pA nel soma
for (let k = 0; k < 40000; k++) n.step(); // 1 s
console.log(n.det.soma.times.length, 'spike');
```

## Limiti

Il modello lavora a livello di popolazioni di canali. Non include la geometria 3D reale con le spine,
la dinamica delle concentrazioni di Na⁺/K⁺/Cl⁻, il gating stocastico del singolo canale, le cascate
di secondi messaggeri (oltre GABA-B), la plasticità a lungo termine, la neuromodulazione e la glia.
