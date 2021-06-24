'''
NAME
    bio_analysis

VERSION
    [1.0]

AUTHOR
    Daianna Gonzalez Padilla <daianna@lcg.unam.mx>

DESCRIPTION
    This package works with biological data: dna sequences and regulatory networks to obtain general characteristics

CATEGORY
     Biological data analysis

USAGE
    import bio_analysis

MODULES
    DNA.py
    network.py

'''



print(f'Invoking __init__.py for {__name__}')
import bio_analysis.DNA
import bio_analysis.network

