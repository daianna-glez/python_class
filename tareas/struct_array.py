'''
NAME
    struct_array.py
VERSION
    [1.0]
AUTHOR
    Daianna Gonzalez Padilla <daianna@lcg.unam.mx>
DESCRIPTION
    This programs prints the structred arrays of an array with data of the production of metabolite of 4 genes
    at two different temperatures, an array with the cost of induction for each gene and an array with the cost
    per g/L.
    
CATEGORY
     Gene expression analysis
USAGE
    None
ARGUMENTS
    This program doesn't take arguments
INPUT
    The input arrays are already given
OUTPUT
    Structred arrays of the input and output arrays 
EXAMPLES
    Example 1: with the following arrays:   production = np.array([ [5,3], [11,7], [4,9], [2,6]])
                                            costs = np.array([3.5, 5, 7, 4.3])

            prints the structured arrays:   Production for each gene at diferente temperatures: 
                                             [('Gen 1',  5, 3) ('Gen 2', 11, 7) ('Gen 3',  4, 9) ('Gen 4',  2, 6)]
                                            Cost of induction for each gene: 
                                             [('Gen 1', 3.5) ('Gen 2', 5. ) ('Gen 3', 7. ) ('Gen 4', 4.3)]
                                            Cost of induction per g/L for each gene: 
                                             [('Gen 1', 0.7   , 1.1666) ('Gen 2', 0.4545, 0.7142)
                                             ('Gen 3', 1.75  , 0.777 ) ('Gen 4', 2.15  , 0.7166)]

               
SOURCE
    https://github.com/daianna21/python_class/blob/master/tareas/struct_array.py
'''

import numpy as np

# Input data
production=np.array([[5,3], [11,7], [4,9], [2,6]])
costs = np.array([3.5, 5, 7, 4.3])
# Cost per g/L is given by cost of each gene/its production at each temperature
# Adjust the dimensions of the matrices by using their transposed to operate with them
costs_per_gL=(costs/production.T).T


# Create the structured arrays: rows are genes and columns different temperatures
production= np.array([('Gen 1', 5, 3), ('Gen 2', 11,7),('Gen 3',4,9),('Gen 4',2,6)],
       dtype=[('gen', (np.str_, 6)), ('Temp1', np.int32), ('Temp2', np.int32)])

print('Production for each gene at diferente temperatures:','\n',production)

# Costs are floats
costs= np.array([('Gen 1', 3.5), ('Gen 2', 5),('Gen 3',7),('Gen 4',4.3)],
       dtype=[('gen', (np.str_, 6)), ('Cost', np.float64)])

print('Cost of induction for each gene:','\n',costs)

# Costs per g/L are also floats
costs_per_gL= np.array([('Gen 1', 0.7, 1.1666), ('Gen 2', 0.4545,0.7142),('Gen 3',1.75,0.777),('Gen 4',2.15,0.7166)],
       dtype=[('gen', (np.str_, 6)), ('Temp1', np.float64), ('Temp2', np.float64)])

print('Cost of induction per g/L for each gene:','\n',costs_per_gL)
