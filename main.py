# Main file

import parser
from sympy import solve

cifre_mostrate = 7

# Equazioni del circuito
equ = parser.get_equation()

# Variabili note (Esempio valori di resistenze e generatori)
content = parser.get_variable_values()

# Incognite
X = parser.X

# Dizionario contente variabile e valore
# 'var': '1000' ad esempio
variabili = {}


# Assegno alle variabili note i propri 
# valori specificati nella netlist
for line in content:
    val = line.split(' ');
    variabili[val[0]] = val[3]
    

# Risolvo equazioni e in forma simbolica
# solutions è un array che contiente le equazioni in forma di stringa
solutions = str(solve(equ, X)).split(',')


# Sostituisco i valori nelle equazioni e ne calcolo il risultato
# Infine stampa gli output
for s in solutions:

    # Elimino eventuali parentesi graffe derivanti dalla conversione
    # Delle equazioni in stringhe
    if s[0] == "{":
        s = s[1:]
    if s[-1:] == "}":
        s = s[:-1]

    # Rimpiazzo le variabili con i propri valori
    for val in variabili:
        s = s.replace(val, variabili[val])
        # Per i componenti reattivi si ha s = jw
        # Siccome in regime di corrente continua si ha 
        # w = 2*π*f, dove f = 0 => w = s = 0
        s = s.replace('s', '0')

    # Elimino eventuali spazi dalla stringa
    s = s.strip()

    # Separo il nome dell'incognita dal suo valore
    #   res[0] -> nome dell'incognita
    #   res[1] -> espressione algebrica con simboli
    res = s.split(":")

    # Unità di misura
    unit = ' '

    # Calcolo il valore dell'incognita
    value = float(str(eval(res[1])))

    # Nel caso i valori siano piccoli, aggiungo milli
    if value < 1:
        value *= 1000;
        unit += 'm'

    # Imposto l'unità di misura in funzione della prima lettera
    # del nome della variabile
    if res[0][0].lower() == 'i':
        unit += 'A'
    else:
        unit += 'V'

    if '_' in res[0]:
        for val in variabili:
            if variabili[val] in res[0]:
                res[0] = res[0].replace(variabili[val], val)
                break

    # Sommo 1 alle cifre mostrate perche viene contato anche il .
    print(res[0] + " = " + str(value)[:cifre_mostrate + 1] + unit)