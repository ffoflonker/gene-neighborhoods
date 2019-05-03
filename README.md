# gene-neighborhoods
README gene_neighborhoods.pl

Usage:
```
"USAGE: perl gene_neighborhoods.pl <OrthoMCL_output_file> <gff_file> <window_size_integer> <min_ortholog_species_number_integer>  OPTIONS: [--tabular] [--stats] [--method default/cooccurring/select/default_clean/select_clean] 
```

Finds gene neighborhoods of proximal orthologous gene pairs conserved in a minimum number of species

---------------------------------------------------------------------------------

Inputs:

1. Orthofinder output in legacy format or OrthoMCL output

Example:

OG0005365: Bathy_XP_007510238.1 ChlNC64A_e_gw1.9.175.1 Cre_06.g276050.t1 Csub_18799 Czo_Cz09g24110.t1 Knitens_00150_0030_v1.1 MpusCCMP1545_44192 Vocar_0002s0367.1

OG0005366: Bathy_XP_007510256.1 ChlNC64A_IGS.gm_4_00538 Cre_12.g497101.t1 Czo_Cz02g17020.t1 MpusCCMP1545_196356 Ostta4221_3_69920_OT_ostta10g02300T0 Pse_000006F_consensus.g874.t1 Vocar_0006s0387.1

2. gff file of all species. tab delimited. chromosome identifier (format: specieschromosome#), gene identifier (format: species_gene#), start coordinate, end coordinate). 

Example:

Cr9	Cre_09.g388750.t1	2824325	2830821

Cr9	Cre_09.g388763.t1	3373489	3377999

Vc9	Vocar_0009s0028.1	193952	208105

Vc9	Vocar_0009s0029.1	209754	212303


3. window size (integer) -window size in number of genes (Recommended <10)

4. minimum number of orthologs (integer) -filter for gene neighborhoods with proximal orthologous gene pairs from a minimum number of species 


OPTIONS:

--methods

default: find clusters based on minimum number of orthologous proximal genes

cooccur:find genes that always clustered in the species in which they are present (can lower minimum number of orthologs).

select: select for neighborhoods containing certain species. Useful for comparing to outgroup species, can lower minimum number of orthologs.

default_clean: default method of cluster finding, with a cleanup step to remove genes under the minimum ortholog number threshold. Helps remove background synteny. 

select_clean: 'select' method of cluster finding, with a cleanup step to remove genes under the minimum ortholog number threshold and not selected for. Helps remove background synteny.

--stats

  reports cluster statistics. 
  

---------------------------------------------------------------------------------

Outputs:

1. File.txt -full output of sliding windows
2. File.merged.txt -merged windows
3. File.ranked.txt -gene neighborhoods ranked giving weight to genes from smaller orthologous groups
4. File.ranked.uniq.txt

Optional outputs

--tabular

File.tabular.txt
File.html

tabular text and html output showing gene IDs for clustered genes, Present= orthologous gene present in genome, but not clustered, Absent= orthologous gene absent from genome. 

File.heatmap.html

heatmap output: light red= present, dark red = clustering, white=absent. Species and gene neighborhoods are clustered in hierarchical order. 

--stats

File.stats.txt -Statistics on cluster size, homolog number, number of genes per species clustering. 
