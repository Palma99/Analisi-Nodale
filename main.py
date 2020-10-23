# Main file

import parser
from sympy import solve

# Equazioni del circuito
equ = parser.get_equation()

# Variabili note (Esempio valori di resistenze e generatori)
content = parser.get_variable_values()

# Incognite
X = parser.X

# Dizionario contente variabile e valore
# 'var': '1000' ad esempio
variabili = {}

# Risolvo equazioni e in forma simbolica
# solutions è un array che contiente le equazioni in forma di stringa
solutions = str(solve(equ, X)).split(',')

# Assegno alle variabili note i propri 
# valori specificati nella netlist
for line in content:
    val = line.split(' ');
    variabili[val[0]] = val[3]
    

# Sostituisco i valori nelle equazioni e ne calcolo il risultato
# Infine stampa gli output
for s in solutions:
    if s[0] == "{":
        s = s[1:]
    if s[-1:] == "}":
        s = s[:-1]
    for val in variabili:
        s = s.replace(val, variabili[val])
    s = s.strip()
    res = s.split(":")
    print(res[0] + " = " + str(eval(res[1])))