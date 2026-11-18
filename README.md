# Neurona-i

Simulatore biofisico di un neurone piramidale corticale, che gira interamente nel browser.

Due pagine (aprile direttamente, oppure servi la cartella con `python3 -m http.server`):

- `index.html`: un singolo neurone biofisico dettagliato.
- `rete/index.html`: una rete di 49 neuroni spiking che impara da sola a riconoscere le cifre.

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

## La rete (`rete/`)

Un'intelligenza artificiale rudimentale fatta di neuroni spiking, sull'architettura di Diehl e Cook (2015):

- **Retina**: immagini 16×16 trasformate in treni di spike di Poisson (fino a 64 Hz per pixel).
- **Strato eccitatorio**: 49 neuroni integrate-and-fire a conduttanza con inibizione laterale
  (competizione winner-take-all) e soglia adattativa (omeostasi).
- **Apprendimento non supervisionato**: STDP sulle sinapsi retina→neurone e normalizzazione dei pesi.
  Le etichette servono solo a dare un nome a ciò che ogni neurone ha imparato.
- **Dati**: cifre generate al volo (`digits.js`) con deformazioni casuali; sulla pagina si può
  disegnare una cifra a mano per metterla alla prova.

Con circa 2000 esempi la rete riconosce l'82–91% delle cifre nuove (prove in Node su 3 semi casuali).

```js
const Digits = require('./rete/digits.js'), SNN = require('./rete/snn-core.js');
const net = SNN.createNetwork();
for (let k = 0; k < 2000; k++) {
  const d = Math.floor(Math.random() * 10);
  net.observe(net.present(Digits.render(d), { learn: true }), d);
}
net.assignLabels();
console.log(net.classify(net.present(Digits.render(3))).digit);
```
