#Clase 9:

# Data frames de Pandas
import pandas as pd
import numpy as np
pd_DF = pd.DataFrame(np.random.rand(3, 2),
       columns=["columna_1", "columna_2"],
                     index=['a','b','c'])


produccion = pd.Series([5, 11, 4, 7, 2],
    index= ['gen1', 'gen2', 'gen3','gen4', 'gen5'],
    name='production')
costos = pd.Series([ 5, 4.3, 7, 3.5],
    index=['gen1', 'gen2', 'gen3', 'gen5'],
    name='costos')

costo_benecio = pd.DataFrame({'costos':costos, 'produccion':produccion})

##Ejercicio 1:

produccion = pd.Series([5, 11, 4, 7, 2],
        index=['gen1', 'gen2', 'gen3', 'gen4', 'gen5'],
        name='produccion')
costos = pd.Series([ 3.5, 5, 7, 4.3],
        index=['gen1', 'gen2', 'gen3', 'gen5'],
                    name='costos')
costo_beneficio = pd.DataFrame({'costos': costos,
                'produccion': produccion})
costo_beneficio['costo unitario']=costo_beneficio.costos/costo_beneficio.produccion
print(costo_beneficio)
