# Esame finale
L'esame finale consiste in una simulazione di analisi dati ispirata ad un
caso reale (http://resonaances.blogspot.com/2015/12/a-new-boson-at-750-gev.html).

In un ipotetico esperimento su acceleratore viene richiesto di cercare l'esistenza
di una risonanza ad alta massa (~750 GeV) che decade in due fotoni.

L'esercizio consiste in
- Stimare la presenza o meno di un segnale alla massa di ~750 GeV e riportarne
la significanza statistica e p-value.
  - In caso la significanza sia troppo bassa, porre un upper-limit al numero di eventi di segnale.
  - In caso la significanza sia sufficiente, stimare il numero di eventi di segnale e riportarne l’incertezza.

La scelta del metodo con cui stimare il numero di eventi e' lasciata a voi

## Come ottenere i dati
Potete scaricare i dati da https://github.com/valerioformato/StatFinalExam

I dataset presenti sono:
- Dataset 1 (10% statistica): `data/data_low.root`
- Dataset 2 (100% statistica): `data/data_high.root`
- Simulazione Montecarlo (solo background, ~1e6 eventi): `data/mc.root`

## Come leggere i dati
Tutti i file root contengono un `TTree` chiamato `Data`
```
root [1] .ls
TFile**		data/data_low.root
 TFile*		data/data_low.root
  KEY: TTree	Data;1	Data
```

Il tree ha due branch, chiamate `photon1` e `photon2`
```
root [2] Data->Print()
******************************************************************************
*Tree    :Data      : Data                                                   *
*Entries :     9953 : Total =          320756 bytes  File  Size =     298269 *
*        :          : Tree compression factor =   1.07                       *
******************************************************************************
*Br    0 :photon1   : ph1[4]/F                                               *
*Entries :     9953 : Total  Size=     160211 bytes  File Size  =     119298 *
*Baskets :        4 : Basket Size=      32000 bytes  Compression=   1.07     *
*............................................................................*
*Br    1 :photon2   : ph2[4]/F                                               *
*Entries :     9953 : Total  Size=     160211 bytes  File Size  =     119321 *
*Baskets :        4 : Basket Size=      32000 bytes  Compression=   1.07     *
*............................................................................*
```
entrambe sono degli array con dimensione 4 e le componenti del quadrimpulso
sono organizzate come
```cpp
float p[4] = {px, py, pz, E};
```

Per iterare sugli eventi nel tree
```cpp
void ReadData() {

  // 4-vectors are stored with the convention
  // p = (px, py, pz, E)
  float photon1[4], photon2[4];

  auto inFile = TFile::Open("data/data_high.root");
  auto tree = (TTree *)inFile->Get("Data");

  tree->SetBranchAddress("photon1", &photon1);
  tree->SetBranchAddress("photon2", &photon2);

  for (Long64_t iEv = 0; iEv < tree->GetEntries(); iEv++) {
    tree->GetEntry(iEv);

    // do stuff here...
  }
}
```
