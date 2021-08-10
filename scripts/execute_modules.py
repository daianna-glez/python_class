import bio_analysis.DNA


#print(bio_analysis.DNA.__doc__)

dna=input("Insert dna sequence:\n")
#TTTAATGGCCATGGCGCCCAGAACTGAGATCAATAGTACCCGTATTAACGGGTGAAAAAAAATA
print(bio_analysis.DNA.get_orf(dna))
print(bio_analysis.DNA.get_protein(dna))

print(bio_analysis.DNA.atcg_content(dna))
print(bio_analysis.DNA.at_content(dna))
print(bio_analysis.DNA.gc_content(dna))
print(bio_analysis.DNA.at_regions(dna))
print(bio_analysis.DNA.revcomp_seq(dna))


