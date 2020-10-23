#!/usr/bin/env python
# coding: utf-8

# # Symbolic modified nodal analysis
# Last update: 12/9/2017
# 
# **Abstract:** The python code in this [notebook](https://cocalc.com/projects/715a4699-f882-4848-af86-6e5c14f24be2/files/node%20analysis.ipynb?session=default) will read in a spice like circuit netlist file and formulate a set of network equations in symbolic form using sympy. These equations can then be copied to a different notebook where the node voltages can be numerically solved using sympy or numpy.  Linear resistors, capacitors, inductors, independent sources and controlled sources are supported.
# 
# **Introduction:** This node analysis code started as a translation from some C code to generate a nodal admittance matrix that I had written in 1988.  I wrote this code for two reasons.  Free versions of Spice for the PC didn't exist at the time and I wanted to use some of the code from the Numerical Recipes in C [[1]](#ref1) Book.  The original C code worked well and calculated numeric solutions.  I then started writing some C code to generate the matrices with symbolic values and then intended to use LISP to symbolically solve the equations.  I didn’t get too far with this effort.  The LISP code would generate huge symbolic strings with no simplification.  The output was a big pile of trash that was not in the least bit useful or decipherable.
# 
# In 2014, I started to use python for my little coding projects and engineering calculations.  There are some nice python libraries for numeric and symbolic calculations (such as numpy and sympy), so I decided to try writing a python script to generate the node equations based on the old C code I had written many years before.  Part way into this project I discovered that there is a new nodal analysis technique being taught today in engineering school called the modified nodal analysis [[2]](#ref2)[[3]](#ref3).  My motivation for reviving this coding project is my continued interest in circuit analysis and synthesis.  
# 
# **Description:** The modified nodal analysis provides an algorithmic method for generating systems of independent equations for linear circuit analysis.  Some of my younger colleagues at work were taught this method, but I never heard of it until a short time ago.  These days, I never really analyze a circuit by hand, unless it’s so simple that you can almost do it by inspection.  Most problems that an electrical engineer encounters on the job are complex enough that they use computers to analyze the circuits.  LTspice [[4]](#ref4) is the version of spice that I use, since it’s free and does a good job converging when analyzing switching circuits.
# 
# My code started initially by following Erik Cheever's Analysis of Resistive Circuits, reference [[5]](#ref5) MATLAB code, to generate modified nodal equations. I somewhat followed his MATLAB file for resistors, capacitors, opamps and independent sources.  The naming of the matrices follows his convention.  The preprocessor and parser code was converted from my old C code.  The use of pandas for a data frame is new and sympy [[6]](#ref6) is used to do the math and the use of element stamps is from reference [[7]](#ref7).
# 
# Inductors are being addressed in the D matrix.  Erik's code puts inductors into the G matrix as 1/s/L.  My code puts the inductor contribution into the D matrix and the unknown current from the unductor into the B and C matricies.  Coupled inductors also affect the D matrix, so it makes sense to allow the inductors to be in the D matrix rather than the G matrix.
# 
# **Network equations:** The network equations are a set of independent equations expressed in this code in matrix form.  There is an equation for each node based on Kirchhoff's current law (KCL) [[8]](#ref8) and an equation for each current unknown.   The current unknowns are the currents from the voltages sources, op amps, voltage controlled voltage sources, current controlled voltage sources, current controlled current sources and inductors.
# 
# Equation 1 is the network equations in matrix form.  
# 
# $A\cdot X = Z$
# 
# The A matrix describes the connectivity of the resistors, capacitors and G type (VCCS) circuit elements.  The column vector X are the unknown node voltages and unknown currents terms from the voltage sources and inductors.  The column vector Z is made of the known voltages and currents.  The A is formed by four sub matrices, G, B, C and D, which are described below.
# 
# $A = \begin{bmatrix}G B\\C D\end{bmatrix}$
# 
# The matrix G is formed from the coefficients representing the KCL equations for each node.
# The positive diagonal of G$_{k,k}$ are the conductance terms of the resistor and capacitor elements connected to node k.  The off diagonal terms of G$_{k,j}$ are the resistors and capacitor conductances connecting node k to node j.  G type elements (VCCS) have input to the G matrix at the connection and controlling node positions.
# 
# The B matrix describes the connectivity of the unknown branch currents.  Independent voltage sources, opamps, H, F and E type elements as well as inductors have inputs to the B matrix.
# 
# The C matrix describes the connectivity of the unknown branch currents and is mainly the transpose of B matrix, with the exception of the F type elements (CCCS) and includes the E type value. 
# 
# The D matrix describes also connectivity of the unknown currents.  The D matrix is composed of zeros unless there are controlled sources and inductors in the network.  
# 
# The X vector is comprised of the V and J vectors as shown below.   
# $X = \begin{bmatrix}V\\J\end{bmatrix}$  
# The V vector contains the node voltages which are the voltage unknowns to be solved for.  The J vector contains the unknown currents from each voltage source.
# 
# The Z vector is comprised of the I and Ev vectors as shown below.  
# $Z = \begin{bmatrix}I\\Ev\end{bmatrix}$  
# The I vector contains the known currents and the Ev vector contains the known voltages.  Ev is used as the variable because sympy uses e and E sometimes for the constant e=2.71, sometimes called Euler's number [[9]](#ref9).  The use of E or e as a symbol was causing some errors when the code was run.  
# 
# Putting all the parts together:
# 
# $\begin{bmatrix}G B\\C D\end{bmatrix} \cdot \begin{bmatrix}V\\J\end{bmatrix} = \begin{bmatrix}I\\Ev\end{bmatrix}$
# 
# **Stamps:** Stamps are templates for modifying the B, C and D matrices and facilitate the construction of the matrices. The stamps used in this implementation of the MNA follow the stamps of reference [[7]](#ref7).  
# 
# **Code description:**  The code is divided in the following sections.  
# Preprocessor:  The preprocessor reads in the netlist text file and removes comments, extra spaces and blank lines.  The first letter of the element type is capitalized to make subsequent parsing of the file easier.  The number of lines are counted and the number of entries on each line are checked to make sure the count is consistent with the element type.
# 
# Parser:  The parser code loads the preprocessed netlist into a data frame.  A report is generated which consists of a count of the element types in the netlist.  
# 
# Matrix formulation: Each of the matrices and vectors are generated.  
# 
# Circuit equation generation:  The circuit equations are generated in a for loop.  Sympy automatically does some simplification according to its default settings.  Two for loops perform the matrix multiplication on equation 1.  The laplace variable s is used when inductors and capacitors are included in the circuit[[7]](#ref7).  
# 
# **Code validation:**  The python code was verified by analyzing test circuits and comparing the results to LTspice.  A collection of worked circuits can be found in reference [[8]](#ref8).  See the revision history below for an indication of validation performed so far.  Other test circuits can be found in reference [[11]](#ref11).  [Test vectors](https://github.com/Tiburonboy/Symbolic-modified-nodal-analysis/tree/master/test%20circuits) can be found in the github repository for this project.  See the readme file in the test circuits folder for a discription of the validation tests.  Code validation is on going at this time.  
# 
# **Usage:**  The input file is a text file called the net list.  It can generated by using a text editor using the format listed below or by drawing the schematic and exporting the net list.  LTspice can be used to draw the schematic of the circuit to be analyzed.  The bit mapped image of the schematic can be copied and included in a document.  The net list can be imported into the python code and the circuit equations can be generated in symbolic form.  The [User’s guide](https://github.com/Tiburonboy/Symbolic-modified-nodal-analysis/blob/master/user_guide.md) can be found on github repository.
# 
# **Change log:**  The code development change log can be found on github [here](https://github.com/Tiburonboy/Symbolic-modified-nodal-analysis/blob/master/Change%20Log.md).
# 
# The backup history can also be found on github.
# [Backups](https://github.com/Tiburonboy/Symbolic-modified-nodal-analysis/tree/master/backup)
# 
# **Survey of other symbolic circuit analysis code:**  The python code presented in this notebook is somewhat unique since python is open source, free and runs on a variety of platforms.  The code presented in this ipython notebook is portable.  As described below, this code is made available under a public domain licence and archived in a github [repository](https://github.com/Tiburonboy/Symbolic-modified-nodal-analysis).  
# 
# There are other symbolic circuit analysis codes available and some of these are described here.  Some of these codes are based on commercial software such as MATLAB [[12]](#ref12), TINA  [[13]](#ref13) and Maple [[14]](#ref14).  
# 
# [SLiCAP](https://www.analog-electronics.eu/slicap/slicap.html) is a symbolic linear analysis tool.  SLiCAP runs in MATLAB.  
# 
# [TINA](https://www.tina.com) is an acronym of Toolkit for Interactive Network Analysis.  The TINA design suite is a circuit simulator and PCB design software package for analyzing, designing, and real time testing of analog, digital, HDL, MCU, and mixed electronic circuits and their PCB layouts. TINA has some [symbolic analysis capability](https://www.tina.com/symbolic-analysis).
# 
# Maple is a mathematical package and there is an application [note](https://www.maplesoft.com/applications/view.aspx?SID=1427) available describing it use in symbolic circuit analysis. The application note presents an method for evaluating, solving and designing a common, but not so simple pulse-mode high-gain transimpedance amplifier or TIA circuit. 
# 
# [Symbolic Circuit Analysis](https://rodanski.net/ben/work/symbolic/index.htm) is a web page devoted to symbolic circuit analysis.  
# 
# [SAPWIN](http://www.ewh.ieee.org/soc/es/May2001/12/Begin.htm) is a windows program package for symbolic and numerical simulation of analog circuits.
# 
# [Lcapy](https://github.com/mph-/lcapy) is an experimental Python package for teaching linear circuit analysis. It uses SymPy for symbolic mathematics.  
# 
# **License:**  This work (includes python code, documentation, test circuits, etc.) is licensed under a Creative Commons Attribution-ShareAlike 4.0 International License.  
# Share — copy and redistribute the material in any medium or format  
# Adapt — remix, transform, and build upon the material for any purpose, even commercially.  
# https://creativecommons.org/licenses/by-sa/4.0/  
# 
# <img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by-sa/4.0/88x31.png" />
# 
# ---
# **References:**  
# <a id='ref1'></a>
# 1. Numerical Recipes in C: The Art of Scientific Computing, William H. Press, Brian P. Flannery, Saul A. Teukolsky, William T. Vetterling, Cambridge University Press; 1988
# <a id='ref2'></a>
# 1. The modified nodal approach to network analysis, Chung-Wen Ho, A. Ruehli, P. Brennan, IEEE Transactions on Circuits and Systems ( Volume: 22, Issue: 6, Jun 1975 )
# <a id='ref3'></a>
# 1. [Modified nodal analysis](https://en.wikipedia.org/wiki/Modified_nodal_analysis), wikipedia.org, retrieved October 6, 2017
# <a id='ref4'></a>
# 1. [LTspice](http://www.linear.com/solutions/ltspice), Linear Technology Corporation, retrieved October 6, 2017
# <a id='ref5'></a>
# 1. [Analysis of  Resistive Circuits](http://www.swarthmore.edu/NatSci/echeeve1/Ref/mna/MNA1.html), retrieved October 6, 2017
# <a id='ref6'></a>
# 1. [Sympy](https://www.scipy.org/), Scipy.org, retrieved October 8, 2017
# <a id='ref7'></a>
# 1. [Laplace transform](https://en.wikipedia.org/wiki/Laplace_transform), wikipedia.org, retrieved December 3, 2017
# <a id='ref8'></a>
# 1. ECE 570 Session 3, Computer Aided Engineering for Integrated Circuits, http://www2.engr.arizona.edu/~ece570/session3.pdf
# <a id='ref9'></a>
# 1. [Kirchhoff's circuit laws](https://en.wikipedia.org/wiki/Kirchhoff%27s_circuit_laws), Wikipedia.com, retrieved October 8, 2017
# <a id='ref10'></a>
# 1. [e (mathematical constant)](https://en.wikipedia.org/wiki/E_(mathematical_constant)), Wikipedia.com, retrieved October 8, 2017
# <a id='ref11'></a>
# 1. Solved Problems, A Source of Free Solved Problems,[Category Archives: Electrical Circuits](http://www.solved-problems.com/circuits/electrical-circuits-problems/716/supernode-dependent-voltage-source/), retrieved October 6, 2017
# <a id='ref12'></a>
# 1. [MATLAB](https://www.mathworks.com/products/matlab.html), retrieved October 6, 2017
# <a id='ref13'></a>
# 1. [TINA](https://www.tina.com/), retrieved October 6, 2017
# <a id='ref14'></a>
# 1. [Maple](https://www.maplesoft.com/), retrieved October 6, 2017
# 

import os
from sympy import *
import numpy as np
import pandas as pd
init_printing()


# initialize variables
num_rlc = 0 # number of passive elements
num_ind = 0 # number of inductors
num_v = 0    # number of independent voltage sources
num_i = 0    # number of independent current sources
i_unk = 0  # number of current unknowns
num_opamps = 0   # number of op amps
num_vcvs = 0     # number of controlled sources of various types
num_vccs = 0
num_cccs = 0
num_ccvs = 0
num_cpld_ind = 0 # number of coupled inductors


# ## Open net list and preprocess it
# The following steps are performed:  
# 1. file name extenstion is defaulted to .net
# 2. remove blank lines and comments
# 3. convert first letter of element name to upper case
# 4. removes extra spaces between entries
# 5. count number of entries on each line, make sure the count is correct, count each element type


# Richedo all'utente il nome della netlist
fn = input("Nome della netlist (Il file .net deve essere nella cartella corrente): ")
fn = fn.replace('.net', '')
try:
    fd1 = open(fn + '.net', 'r')
except:
    print('\nErrore nell\'apertura del file. Assicurati che il nome sia giusto e che sia nella cartella corrente.')
    exit(-1)
content = fd1.readlines()
content = [x.strip() for x in content]  #remove leading and trailing white space
# remove empty lines
while '' in content:
    content.pop(content.index(''))

# remove comment lines, these start with a asterisk *
content = [n for n in content if not n.startswith('*')]
# remove other comment lines, these start with a semicolon ;
content = [n for n in content if not n.startswith(';')]
# remove spice directives, these start with a period, .
content = [n for n in content if not n.startswith('.')]
# converts 1st letter to upper case
#content = [x.upper() for x in content] <- this converts all to upper case
content = [x.capitalize() for x in content]
# removes extra spaces between entries
content = [' '.join(x.split()) for x in content]
# print(content)


line_cnt = len(content) # number of lines in the netlist
branch_cnt = 0  # number of branches in the netlist
# check number of entries on each line, count each element type
for i in range(line_cnt):
    x = content[i][0]
    tk_cnt = len(content[i].split()) # split the line into a list of words

    if (x == 'R') or (x == 'L') or (x == 'C'):
        if tk_cnt != 4:
            print("branch {:d} not formatted correctly, {:s}".format(i,content[i]))
            print("had {:d} items and should only be 4".format(tk_cnt))
        num_rlc += 1
        branch_cnt += 1
        if x == 'L':
            num_ind += 1
    elif x == 'V':
        if tk_cnt != 4:
            print("branch {:d} not formatted correctly, {:s}".format(i,content[i]))
            print("had {:d} items and should only be 4".format(tk_cnt))
        num_v += 1
        branch_cnt += 1
    elif x == 'I':
        if tk_cnt != 4:
            print("branch {:d} not formatted correctly, {:s}".format(i,content[i]))
            print("had {:d} items and should only be 4".format(tk_cnt))
        num_i += 1
        branch_cnt += 1
    elif x == 'O':
        if tk_cnt != 4:
            print("branch {:d} not formatted correctly, {:s}".format(i,content[i]))
            print("had {:d} items and should only be 4".format(tk_cnt))
        num_opamps += 1
    elif x == 'E':
        if (tk_cnt != 6):
            print("branch {:d} not formatted correctly, {:s}".format(i,content[i]))
            print("had {:d} items and should only be 6".format(tk_cnt))
        num_vcvs += 1
        branch_cnt += 1
    elif x == 'G':
        if (tk_cnt != 6):
            print("branch {:d} not formatted correctly, {:s}".format(i,content[i]))
            print("had {:d} items and should only be 6".format(tk_cnt))
        num_vccs += 1
        branch_cnt += 1
    elif x == 'F':
        if (tk_cnt != 5):
            print("branch {:d} not formatted correctly, {:s}".format(i,content[i]))
            print("had {:d} items and should only be 5".format(tk_cnt))
        num_cccs += 1
        branch_cnt += 1
    elif x == 'H':
        if (tk_cnt != 5):
            print("branch {:d} not formatted correctly, {:s}".format(i,content[i]))
            print("had {:d} items and should only be 5".format(tk_cnt))
        num_ccvs += 1
        branch_cnt += 1
    elif x == 'K':
        if (tk_cnt != 4):
            print("branch {:d} not formatted correctly, {:s}".format(i,content[i]))
            print("had {:d} items and should only be 4".format(tk_cnt))
        num_cpld_ind += 1
    else:
        print("unknown element type in branch {:d}, {:s}".format(i,content[i]))


# ## Parser
# The parser performs the following operations.
# 1. puts branch elements into data frame
# 2. counts number of nodes
# 
# data frame lables:
# - element: type of element
# - p node: positive node
# - n node: negitive node, for a current source, the arrow point terminal, LTspice puts the inductor phasing dot on this terminal
# - cp node: controlling positive node of branch
# - cn node: controlling negitive node of branch
# - Vout: opamp output node
# - value: value of element or voltage
# - Vname: voltage source through which the controlling current flows. Need to add a zero volt voltage source to the controlling branch.
# - Lname1: name of coupled inductor 1
# - Lname2: name of coupled inductor 2




# build the pandas data frame
df = pd.DataFrame(columns=['element','p node','n node','cp node','cn node',
    'Vout','value','Vname','Lname1','Lname2'])

# this data frame is for branches with unknown currents
df2 = pd.DataFrame(columns=['element','p node','n node'])


# ### Functions to load branch elements into data frame and check for gaps in node numbering

# loads voltage or current sources into branch structure
def indep_source(line_nu):
    tk = content[line_nu].split()
    df.loc[line_nu,'element'] = tk[0]
    df.loc[line_nu,'p node'] = int(tk[1])
    df.loc[line_nu,'n node'] = int(tk[2])
    df.loc[line_nu,'value'] = float(tk[3])

# loads passive elements into branch structure
def rlc_element(line_nu):
    tk = content[line_nu].split()
    df.loc[line_nu,'element'] = tk[0]
    df.loc[line_nu,'p node'] = int(tk[1])
    df.loc[line_nu,'n node'] = int(tk[2])
    df.loc[line_nu,'value'] = float(tk[3])

# loads multi-terminal elements into branch structure
# O - Op Amps
def opamp_sub_network(line_nu):
    tk = content[line_nu].split()
    df.loc[line_nu,'element'] = tk[0]
    df.loc[line_nu,'p node'] = int(tk[1])
    df.loc[line_nu,'n node'] = int(tk[2])
    df.loc[line_nu,'Vout'] = int(tk[3])

# G - VCCS
def vccs_sub_network(line_nu):
    tk = content[line_nu].split()
    df.loc[line_nu,'element'] = tk[0]
    df.loc[line_nu,'p node'] = int(tk[1])
    df.loc[line_nu,'n node'] = int(tk[2])
    df.loc[line_nu,'cp node'] = int(tk[3])
    df.loc[line_nu,'cn node'] = int(tk[4])
    df.loc[line_nu,'value'] = float(tk[5])

# E - VCVS
# in sympy E is the number 2.718, replacing E with Ea otherwise, sympify() errors out
def vcvs_sub_network(line_nu):
    tk = content[line_nu].split()
    df.loc[line_nu,'element'] = tk[0].replace('E', 'Ea')
    df.loc[line_nu,'p node'] = int(tk[1])
    df.loc[line_nu,'n node'] = int(tk[2])
    df.loc[line_nu,'cp node'] = int(tk[3])
    df.loc[line_nu,'cn node'] = int(tk[4])
    df.loc[line_nu,'value'] = float(tk[5])

# F - CCCS
def cccs_sub_network(line_nu):
    tk = content[line_nu].split()
    df.loc[line_nu,'element'] = tk[0]
    df.loc[line_nu,'p node'] = int(tk[1])
    df.loc[line_nu,'n node'] = int(tk[2])
    df.loc[line_nu,'Vname'] = tk[3].capitalize()
    df.loc[line_nu,'value'] = float(tk[4])

# H - CCVS
def ccvs_sub_network(line_nu):
    tk = content[line_nu].split()
    df.loc[line_nu,'element'] = tk[0]
    df.loc[line_nu,'p node'] = int(tk[1])
    df.loc[line_nu,'n node'] = int(tk[2])
    df.loc[line_nu,'Vname'] = tk[3].capitalize()
    df.loc[line_nu,'value'] = float(tk[4])

# K - Coupled inductors
def cpld_ind_sub_network(line_nu):
    tk = content[line_nu].split()
    df.loc[line_nu,'element'] = tk[0]
    df.loc[line_nu,'Lname1'] = tk[1].capitalize()
    df.loc[line_nu,'Lname2'] = tk[2].capitalize()
    df.loc[line_nu,'value'] = float(tk[3])

# function to scan df and get largest node number
def count_nodes():
    # need to check that nodes are consecutive
    # fill array with node numbers
    p = np.zeros(line_cnt+1)
    for i in range(line_cnt):
        # need to skip coupled inductor 'K' statements
        if df.loc[i,'element'][0] != 'K': #get 1st letter of element name
            p[df['p node'][i]] = df['p node'][i]
            p[df['n node'][i]] = df['n node'][i]

    # find the largest node number
    if df['n node'].max() > df['p node'].max():
        largest = df['n node'].max()
    else:
        largest =  df['p node'].max()

    largest = int(largest)
    # check for unfilled elements, skip node 0
    for i in range(1,largest):
        if p[i] == 0:
            print('nodes not in continuous order, node {:.0f} is missing'.format(p[i-1]+1))

    return largest


# ### Load circuit net list into the data frames




# load branch info into data frame
for i in range(line_cnt):
    x = content[i][0]

    if (x == 'R') or (x == 'L') or (x == 'C'):
        rlc_element(i)
    elif (x == 'V') or (x == 'I'):
        indep_source(i)
    elif x == 'O':
        opamp_sub_network(i)
    elif x == 'E':
        vcvs_sub_network(i)
    elif x == 'G':
        vccs_sub_network(i)
    elif x == 'F':
        cccs_sub_network(i)
    elif x == 'H':
        ccvs_sub_network(i)
    elif x == 'K':
        cpld_ind_sub_network(i)
    else:
        print("unknown element type in branch {:d}, {:s}".format(i,content[i]))

# count number of nodes
num_nodes = count_nodes()

# Build df2: consists of branches with current unknowns, used for C & D matrices
# walk through data frame and find these parameters
count = 0
for i in range(len(df)):
    # process all the elements creating unknown currents
    x = df.loc[i,'element'][0]   #get 1st letter of element name
    if (x == 'L') or (x == 'V') or (x == 'O') or (x == 'E') or (x == 'H') or (x == 'F'):
        df2.loc[count,'element'] = df.loc[i,'element']
        df2.loc[count,'p node'] = df.loc[i,'p node']
        df2.loc[count,'n node'] = df.loc[i,'n node']
        count += 1


# ## Print net list report

# 


# # print a report
# print('Net list report')
# print('number of lines in netlist: {:d}'.format(line_cnt))
# print('number of branches: {:d}'.format(branch_cnt))
# print('number of nodes: {:d}'.format(num_nodes))
# # count the number of element types that affect the size of the B, C, D, E and J arrays
# # these are current unknows
# i_unk = num_v+num_opamps+num_vcvs+num_ccvs+num_cccs+num_ind
# print('number of unknown currents: {:d}'.format(i_unk))
# print('number of RLC (passive components): {:d}'.format(num_rlc))
# print('number of inductors: {:d}'.format(num_ind))
# print('number of independent voltage sources: {:d}'.format(num_v))
# print('number of independent current sources: {:d}'.format(num_i))
# print('number of op amps: {:d}'.format(num_opamps))
# print('number of E - VCVS: {:d}'.format(num_vcvs))
# print('number of G - VCCS: {:d}'.format(num_vccs))
# print('number of F - CCCS: {:d}'.format(num_cccs))
# print('number of H - CCVS: {:d}'.format(num_ccvs))
# print('number of K - Coupled inductors: {:d}'.format(num_cpld_ind))



# store the data frame as a pickle file
# df.to_pickle(fn+'.pkl')  # <- uncomment if needed


# initialize some symbolic matrix with zeros
# A is formed by [[G, C] [B, D]]
# Z = [I,E]
# X = [V, J]
V = zeros(num_nodes,1)
I = zeros(num_nodes,1)
G = zeros(num_nodes,num_nodes)  # also called Yr, the reduced nodal matrix
s = Symbol('s')  # the Laplace variable

# count the number of element types that affect the size of the B, C, D, E and J arrays
# these are element types that have unknown currents
i_unk = num_v+num_opamps+num_vcvs+num_ccvs+num_ind+num_cccs
# if i_unk == 0, just generate empty arrays
B = zeros(num_nodes,i_unk)
C = zeros(i_unk,num_nodes)
D = zeros(i_unk,i_unk)
Ev = zeros(i_unk,1)
J = zeros(i_unk,1)


# #### some debugging notes:
# Is is possible to have i_unk == 0 ?, what about a network with only current sources?  This would make B = 0 for example. Did one test, need to run others  
# Is there a valid op amp case where B is n by 1?

# ## G matrix
# The G matrix is n by n, where n is the number of nodes. The matrix is formed by the interconnections between the resistors, capacitors and VCCS type elements.  In the original paper G is called Yr, where Yr, is a reduced form of the nodal matrix excluding the contributions due to voltage sources, current controlling elements, etc.  In python row and columns are: G[row, column]
# 




# G matrix
for i in range(len(df)):  # process each row in the data frame
    n1 = df.loc[i,'p node']
    n2 = df.loc[i,'n node']
    cn1 = df.loc[i,'cp node']
    cn2 = df.loc[i,'cn node']
    # process all the passive elements, save conductance to temp value
    x = df.loc[i,'element'][0]   #get 1st letter of element name
    if x == 'R':
        g = 1/sympify(df.loc[i,'element'])
    if x == 'C':
        g = s*sympify(df.loc[i,'element'])
    if x == 'G':   #vccs type element
        g = sympify(df.loc[i,'element'].lower())  # use a symbol for gain value

    if (x == 'R') or (x == 'C'):
        # If neither side of the element is connected to ground
        # then subtract it from appropriate location in matrix.
        if (n1 != 0) and (n2 != 0):
            G[n1-1,n2-1] += -g
            G[n2-1,n1-1] += -g

        # If node 1 is connected to ground, add element to diagonal of matrix
        if n1 != 0:
            G[n1-1,n1-1] += g

        # same for for node 2
        if n2 != 0:
            G[n2-1,n2-1] += g

    if x == 'G':    #vccs type element
        # check to see if any terminal is grounded
        # then stamp the matrix
        if n1 != 0 and cn1 != 0:
            G[n1-1,cn1-1] += g

        if n2 != 0 and cn2 != 0:
            G[n2-1,cn2-1] += g

        if n1 != 0 and cn2 != 0:
            G[n1-1,cn2-1] -= g

        if n2 != 0 and cn1 != 0:
            G[n2-1,cn1-1] -= g

G  # display the G matrix


# ## B Matrix
# The B matrix is an n by m matrix with only 0, 1 and -1 elements, where n = number of nodes and m is the number of current unknowns, i_unk. There is one column for each unknown current. The code loop through all the branches and process elements that have stamps for the B matrix:  
# - Voltage sources (V)
# - Opamps (O)
# - CCVS (H)
# - CCCS (F)
# - VCVS (E)
# - Inductors (L)  
# 
# The order of the columns is as they appear in the netlist.  CCCS (F) does not get its own column because the controlling current is through a zero volt voltage source, called Vname and is already in the net list.




# generate the B Matrix
sn = 0   # count source number as code walks through the data frame
for i in range(len(df)):
    n1 = df.loc[i,'p node']
    n2 = df.loc[i,'n node']
    n_vout = df.loc[i,'Vout'] # node connected to op amp output

    # process elements with input to B matrix
    x = df.loc[i,'element'][0]   #get 1st letter of element name
    if x == 'V':
        if i_unk > 1:  #is B greater than 1 by n?, V
            if n1 != 0:
                B[n1-1,sn] = 1
            if n2 != 0:
                B[n2-1,sn] = -1
        else:
            if n1 != 0:
                B[n1-1] = 1
            if n2 != 0:
                B[n2-1] = -1
        sn += 1   #increment source count
    if x == 'O':  # op amp type, output connection of the opamg goes in the B matrix
        B[n_vout-1,sn] = 1
        sn += 1   # increment source count
    if (x == 'H') or (x == 'F'):  # H: ccvs, F: cccs,
        if i_unk > 1:  #is B greater than 1 by n?, H, F
            # check to see if any terminal is grounded
            # then stamp the matrix
            if n1 != 0:
                B[n1-1,sn] = 1
            if n2 != 0:
                B[n2-1,sn] = -1
        else:
            if n1 != 0:
                B[n1-1] = 1
            if n2 != 0:
                B[n2-1] = -1
        sn += 1   #increment source count
    if x == 'E':   # vcvs type, only ik column is altered at n1 and n2
        if i_unk > 1:  #is B greater than 1 by n?, E
            if n1 != 0:
                B[n1-1,sn] = 1
            if n2 != 0:
                B[n2-1,sn] = -1
        else:
            if n1 != 0:
                B[n1-1] = 1
            if n2 != 0:
                B[n2-1] = -1
        sn += 1   #increment source count
    if x == 'L':
        if i_unk > 1:  #is B greater than 1 by n?, L
            if n1 != 0:
                B[n1-1,sn] = 1
            if n2 != 0:
                B[n2-1,sn] = -1
        else:
            if n1 != 0:
                B[n1-1] = 1
            if n2 != 0:
                B[n2-1] = -1
        sn += 1   #increment source count

# check source count
if sn != i_unk:
    print('source number, sn={:d} not equal to i_unk={:d} in matrix B'.format(sn,i_unk))

B   # display the B matrix


# ## C matrix
# The C matrix is an m by n matrix with only 0, 1 and -1 elements (except for controlled sources).  The code is similar to the B matrix code, except the indices are swapped.   The code loops through all the branches and process elements that have stamps for the C matrix:  
# - Voltage sources (V)
# - Opamps (O)
# - CCVS (H)
# - CCCS (F)
# - VCVS (E)
# - Inductors (L)
# 
# References use in the debugging of the opamp stamp:  
# Design of Analog Circuits Through Symbolic Analysis
# edited by Mourad Fakhfakh, Esteban Tlelo-Cuautle, Francisco V. Fernández  
# Computer Aided Design and Design Automation
# edited by Wai-Kai Chen  
# http://users.ecs.soton.ac.uk/mz/CctSim/chap1_4.htm
# 


# find the the column position in the C and D matrix for controlled sources
# needs to return the node numbers and branch number of controlling branch
def find_vname(name):
    # need to walk through data frame and find these parameters
    for i in range(len(df2)):
        # process all the elements creating unknown currents
        if name == df2.loc[i,'element']:
            n1 = df2.loc[i,'p node']
            n2 = df2.loc[i,'n node']
            return n1, n2, i  # n1, n2 & col_num are from the branch of the controlling element

    print('failed to find matching branch element in find_vname')


# generate the C Matrix
sn = 0   # count source number as code walks through the data frame
for i in range(len(df)):
    n1 = df.loc[i,'p node']
    n2 = df.loc[i,'n node']
    cn1 = df.loc[i,'cp node'] # nodes for controlled sources
    cn2 = df.loc[i,'cn node']
    n_vout = df.loc[i,'Vout'] # node connected to op amp output

    # process elements with input to B matrix
    x = df.loc[i,'element'][0]   #get 1st letter of element name
    if x == 'V':
        if i_unk > 1:  #is B greater than 1 by n?, V
            if n1 != 0:
                C[sn,n1-1] = 1
            if n2 != 0:
                C[sn,n2-1] = -1
        else:
            if n1 != 0:
                C[n1-1] = 1
            if n2 != 0:
                C[n2-1] = -1
        sn += 1   #increment source count

    if x == 'O':  # op amp type, input connections of the opamp go into the C matrix
        # C[sn,n_vout-1] = 1
        if i_unk > 1:  #is B greater than 1 by n?, O
            # check to see if any terminal is grounded
            # then stamp the matrix
            if n1 != 0:
                C[sn,n1-1] = 1
            if n2 != 0:
                C[sn,n2-1] = -1
        else:
            if n1 != 0:
                C[n1-1] = 1
            if n2 != 0:
                C[n2-1] = -1
        sn += 1   # increment source count

    if x == 'F':  # need to count F (cccs) types
        sn += 1   #increment source count
    if x == 'H':  # H: ccvs
        if i_unk > 1:  #is B greater than 1 by n?, H
            # check to see if any terminal is grounded
            # then stamp the matrix
            if n1 != 0:
                C[sn,n1-1] = 1
            if n2 != 0:
                C[sn,n2-1] = -1
        else:
            if n1 != 0:
                C[n1-1] = 1
            if n2 != 0:
                C[n2-1] = -1
        sn += 1   #increment source count
    if x == 'E':   # vcvs type, ik column is altered at n1 and n2, cn1 & cn2 get value
        if i_unk > 1:  #is B greater than 1 by n?, E
            if n1 != 0:
                C[sn,n1-1] = 1
            if n2 != 0:
                C[sn,n2-1] = -1
            # add entry for cp and cn of the controlling voltage
            if cn1 != 0:
                C[sn,cn1-1] = -sympify(df.loc[i,'element'].lower())
            if cn2 != 0:
                C[sn,cn2-1] = sympify(df.loc[i,'element'].lower())
        else:
            if n1 != 0:
                C[n1-1] = 1
            if n2 != 0:
                C[n2-1] = -1
            vn1, vn2, df2_index = find_vname(df.loc[i,'Vname'])
            if vn1 != 0:
                C[vn1-1] = -sympify(df.loc[i,'element'].lower())
            if vn2 != 0:
                C[vn2-1] = sympify(df.loc[i,'element'].lower())
        sn += 1   #increment source count

    if x == 'L':
        if i_unk > 1:  #is B greater than 1 by n?, L
            if n1 != 0:
                C[sn,n1-1] = 1
            if n2 != 0:
                C[sn,n2-1] = -1
        else:
            if n1 != 0:
                C[n1-1] = 1
            if n2 != 0:
                C[n2-1] = -1
        sn += 1   #increment source count

# check source count
if sn != i_unk:
    print('source number, sn={:d} not equal to i_unk={:d} in matrix C'.format(sn,i_unk))


# ## D matrix
# The D matrix is an m by m matrix, where m is the number of unknown currents.
# > m = i_unk = num_v+num_opamps+num_vcvs+num_ccvs+num_ind+num_cccs
# 
# Stamps that affect the D matrix are: inductor, ccvs and cccs  
# inductors: minus sign added to keep current flow convention consistent  
# 
# Coupled inductors notes:  
# Can the K statement be anywhere in the net list, even before Lx and Ly?  
# 12/6/2017 doing some debugging on with coupled inductors  
# LTspice seems to put the phasing dot on the neg node when it generates the netlist  
# This code uses M for mutual inductance, LTspice uses k for the coupling coefficient.

# generate the D Matrix
sn = 0   # count source number as code walks through the data frame
for i in range(len(df)):
    n1 = df.loc[i,'p node']
    n2 = df.loc[i,'n node']
    #cn1 = df.loc[i,'cp node'] # nodes for controlled sources
    #cn2 = df.loc[i,'cn node']
    #n_vout = df.loc[i,'Vout'] # node connected to op amp output

    # process elements with input to D matrix
    x = df.loc[i,'element'][0]   #get 1st letter of element name
    if (x == 'V') or (x == 'O') or (x == 'E'):  # need to count V, E & O types
        sn += 1   #increment source count

    if x == 'L':
        if i_unk > 1:  #is D greater than 1 by 1?
            D[sn,sn] += -s*sympify(df.loc[i,'element'])
        else:
            D[sn] += -s*sympify(df.loc[i,'element'])
        sn += 1   #increment source count

    if x == 'H':  # H: ccvs
        # if there is a H type, D is m by m
        # need to find the vn for Vname
        # then stamp the matrix
        vn1, vn2, df2_index = find_vname(df.loc[i,'Vname'])
        D[sn,df2_index] += -sympify(df.loc[i,'element'].lower())
        sn += 1   #increment source count

    if x == 'F':  # F: cccs
        # if there is a F type, D is m by m
        # need to find the vn for Vname
        # then stamp the matrix
        vn1, vn2, df2_index = find_vname(df.loc[i,'Vname'])
        D[sn,df2_index] += -sympify(df.loc[i,'element'].lower())
        D[sn,sn] = 1
        sn += 1   #increment source count

    if x == 'K':  # K: coupled inductors, KXX LYY LZZ value
        # if there is a K type, D is m by m
        vn1, vn2, ind1_index = find_vname(df.loc[i,'Lname1'])  # get i_unk position for Lx
        vn1, vn2, ind2_index = find_vname(df.loc[i,'Lname2'])  # get i_unk position for Ly
        # enter sM on diagonals = value*sqrt(LXX*LZZ)

        D[ind1_index,ind2_index] += -s*sympify('M{:s}'.format(df.loc[i,'element'].lower()[1:]))  # s*Mxx
        D[ind2_index,ind1_index] += -s*sympify('M{:s}'.format(df.loc[i,'element'].lower()[1:]))  # -s*Mxx


# ## V matrix
# The V matrix is an n by 1 matrix formed of the node voltages, where n is the number of nodes. Each element in V corresponds to the voltage at the node.  
# 
# Maybe make small v's v_1 so as not to confuse v1 with V1.

# generate the V matrix
for i in range(num_nodes):
    V[i] = sympify('v{:d}'.format(i+1))


# ## J matrix
# The J matrix is an m by 1 matrix, where m is the number of unknown currents.
# >i_unk = num_v+num_opamps+num_vcvs+num_ccvs+num_ind+num_cccs


# The J matrix is an mx1 matrix, with one entry for each i_unk from a source
#sn = 0   # count i_unk source number
#oan = 0   #count op amp number
for i in range(len(df2)):
    # process all the unknown currents
    J[i] = sympify('I_{:s}'.format(df2.loc[i,'element']))



# ## I matrix
# The I matrix is an n by 1 matrix, where n is the number of nodes. The value of each element of I is determined by the sum of current sources into the corresponding node. If there are no current sources connected to the node, the value is zero.


# generate the I matrix, current sources have n2 = arrow end of the element
for i in range(len(df)):
    n1 = df.loc[i,'p node']
    n2 = df.loc[i,'n node']
    # process all the passive elements, save conductance to temp value
    x = df.loc[i,'element'][0]   #get 1st letter of element name
    if x == 'I':
        g = sympify(df.loc[i,'element'])
        # sum the current into each node
        if n1 != 0:
            I[n1-1] -= g
        if n2 != 0:
            I[n2-1] += g


# ## Ev matrix
# The Ev matrix is mx1 and holds the values of the independent voltage sources.

# generate the E matrix
sn = 0   # count source number
for i in range(len(df)):
    # process all the passive elements
    x = df.loc[i,'element'][0]   #get 1st letter of element name
    if x == 'V':
        Ev[sn] = sympify(df.loc[i,'element'])
        sn += 1


# ## Z matrix
# The Z matrix holds the independent voltage and current sources and is the combination of 2 smaller matrices I and Ev. The Z matrix is (m+n) by 1, n is the number of nodes, and m is the number of independent voltage sources. The I matrix is n by 1 and contains the sum of the currents through the passive elements into the corresponding node (either zero, or the sum of independent current sources). The Ev matrix is m by 1 and holds the values of the independent voltage sources.

Z = I[:] + Ev[:]  # the + operator in python concatinates the lists


# ## X matrix
# The X matrix is an (n+m) by 1 vector that holds the unknown quantities (node voltages and the currents through the independent voltage sources). The top n elements are the n node voltages. The bottom m elements represent the currents through the m independent voltage sources in the circuit. The V matrix is n by 1 and holds the unknown voltages. The J matrix is m by 1 and holds the unknown currents through the voltage sources

X = V[:] + J[:]  # the + operator in python concatinates the lists


# ## A matrix
# The A matrix is (m+n) by (m+n) and will be developed as the combination of 4 smaller matrices, G, B, C, and D.



n = num_nodes
m = i_unk
A = zeros(m+n,m+n)
for i in range(n):
    for j in range(n):
        A[i,j] = G[i,j]

if i_unk > 1:
    for i in range(n):
        for j in range(m):
            A[i,n+j] = B[i,j]
            A[n+j,i] = C[j,i]

    for i in range(m):
        for j in range(m):
            A[n+i,n+j] = D[i,j]

if i_unk == 1:
    for i in range(n):
        A[i,n] = B[i]
        A[n,i] = C[i]

# ## generate the circuit equations


# Funzioni utilizzate nel main.py

# Funzione che ritorna la matrice X, contente le incognite
def get_X():
    return X

# Funzione che ritorna i valori delle variabili note 
def get_variable_values():
    return content

# Definisco una funzione che ritorna le equazioni da risolvere
def get_equation():
    n = num_nodes
    m = i_unk
    eq_temp = 0  # temporary equation used to build up the equation
    equ = zeros(m+n,1)  #initialize the array to hold the equations
    for i in range(n+m):
        for j in range(n+m):
            eq_temp += A[i,j]*X[j]
        equ[i] = Eq(eq_temp,Z[i])
        eq_temp = 0
    return equ

