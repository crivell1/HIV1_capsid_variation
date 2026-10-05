# HIV-1 CAPSID PROTEIN VARIANT INVESTIGATION

The HIV-1 genome codes for 9 genes which translate into 15 different viral proteins. One of the major structural genes known as Gag encodes the core structual matrix and capsid proteins (https://pmc.ncbi.nlm.nih.gov/articles/PMC7471995/). 
Investigation into each viral protein is paramount for full understanding of the genome and how to best target HIV-1 on an individualized basis.
p24 in the HIV-1 GAG gene is the major capsid protein and is the target of the capsid-inhibitor lenacapavir https://pubmed.ncbi.nlm.nih.gov/40623458/. 
Understanding the conservation of the amino acid sequence and the effects on structure between the different HIV-1 variants would allow for better targeting techniques to arise and for drugs to target conserved regions that are exposed for proper drug binding.
p24/CA provides a nice target due to the prevalance of HIV-1 variant sequences and the depth of analysis done into the HIV1 viral genome. Databases such as the Los Alamos National Laboratory contains a large amassment of sequencing data for each HIV gene. 
p24/CA itself has been highly studied and its structure regionally annotated with specific interfaces such as the N-terminus domain, C-terminus domain, linker region, and the corresponding amino acids to alpha helices and B-sheets.


This project is based heavily upon the work of Paloma Troyano-Hernáez et al., 2022
link: https://pmc.ncbi.nlm.nih.gov/articles/PMC9039614/#S5

Figure to recreate:
![<img width="700" height="892" alt="HIV paper fig 1" src="https://github.com/user-attachments/assets/9f694fd2-0013-49a6-ab85-6b71991a4033" />]


This paper utilizes the LANL in order to do a conservation analysis for each variant of HIV1 in the database. 
LANL HIV1 database link: https://www.hiv.lanl.gov/components/sequence/HIV/search/search.html




Our project will follow a similar pipeline with a couple goals in mind.

1a. Align the sequences of the different HIV1 variants for p24/CA \n
1b. Calculate frequency rate of mutation for each variant (Shannon Entropy) \n

2a. Plot Shannon entropy by amino acid in the p24 sequence \n
2b. Visualize variability within structural categories of p24 : NTD vs CTD vs linker, helix vs loop, exposed vs buried residues (SASA calculated from PDB file)

3a. (STRETCH) compare with experimental values of directed mutations and show differences in unexpected high or low variability in strucutal categories \n
3b. (STRETCH) compare other variables within the LANL database variants that encode for 


Software:

Bash/CLI \n
	MUSCLE v5 https://drive5.com/muscle5/ \n
	If Nextstrain is used - intallation instructions: (https://docs.nextstrain.org/en/latest/install.html) \n

Python v3.12.3 \n
	Biopython v1.88 https://biopython.org/docs/latest/Tutorial/chapter_msa.html#muscle \n
	Biotite v1.71 https://pypi.org/project/biotite/ \n
	FreeSASA v2.2.1 https://pypi.org/project/freesasa/ \n
 	pyMSAviz v0.5.0 https://pypi.org/project/pyMSAviz/ \n 
	MatPlotLib v3.9.2 https://pypi.org/project/matplotlib/ \n
	NumPy v1.26.4 https://pypi.org/project/numpy/ \n

