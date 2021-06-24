import bio_analysis.DNA
import bio_analysis.network

print(bio_analysis.DNA.__doc__)
print(print.__doc__)
dna=input("Insert dna sequence:\n")

print(bio_analysis.DNA.get_orf(dna))
print(bio_analysis.DNA.get_protein(dna))

print(bio_analysis.DNA.atcg_content(dna))
print(bio_analysis.DNA.at_content(dna))
print(bio_analysis.DNA.gc_content(dna))
print(bio_analysis.DNA.at_regions(dna))
print(bio_analysis.DNA.revcomp_seq(dna))



f = bio_analysis.network.net("../docs/NewGeneRegAbr03.trn")
print(f.get_edges("narG"))
print(f.neighbors("narG"))
print(f.grades("narG"))
print(f.clustering("narG"))
print(f.P_k(3))
print(f.C_k(3))