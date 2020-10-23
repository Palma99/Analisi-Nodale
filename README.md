# Analisi-Nodale

## Funzionamento
Data una netlist, calcola i valori di tensione e corrente del circuito associato a tale netlist.

La netlist deve essere contenuta in un file dedicato `.net` che deve essere nella stessa cartella degli altri file.

# Librerie
Assicurarsi di avere `python3` e `pip3` installato. Nel caso consultare [qui](https://phoenixnap.com/kb/install-pip-windows).

Installare le seguenti librerie per il corretto funzionamento
- sympy
- numpy
- pandas

Per installare le librerie usare:

`
pip3 install <nome_libreria>
`

# Run
Per eseguire il programma da terminale:

`
  python3 main.py
`

Una volta avviato il programma chiede di inserire il nome della netlist (**Il file deve avere estensione .net e deve essere nella stessa cartella degli altri file**).

Digitare il nome del file in cui è contenuta la netlist e premere `Invio`.

Se non ci sono stati errori nella lettura del file .net, il programma stampa una lista con i valori delle tensioni e correnti incognite.

## Netlist di esempio

```
R1 2 3 10000

R2 3 0 22000

V1 1 0 10

R3 2 1 30000
```

Verranno ignorate tutti i commenti o direttive all'interno della netlist.

**Non utilizzare** prefissi come k, m, u per indicare le potenze del 10.
