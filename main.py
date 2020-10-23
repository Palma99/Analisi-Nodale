import parser
from sympy import solve


equ = parser.get_equation()
content = parser.get_variable_values()
X = parser.X

variabili = {}

# Risolvo equazioni
solutions = str(solve(equ, X)).split(',')
variabili = {}

# Assegno alle variabili note i valori
for line in content:
    val = line.split(' ');
    variabili[val[0]] = val[3]
    
# Sostituisco i valori nelle equazioni e le stampo
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