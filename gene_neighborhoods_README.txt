README gene_neighborhoods.pl

usage: perl gene_neighborhoods.pl <OrthoMCL_output> <gff> <window_size> <min_num_orthologs>

Finds gene neighborhoods of proximal orthologous gene pairs conserved in a minimum number of species

---------------------------------------------------------------------------------

Inputs:

1. Orthofinder output in legacy format or OrthoMCL output

Example:

OG0005365: Bathy_XP_007510238.1 ChlNC64A_e_gw1.9.175.1 Cre_06.g276050.t1 Csub_18799 Czo_Cz09g24110.t1 Knitens_00150_0030_v1.1 MpusCCMP1545_44192 Vocar_0002s0367.1
OG0005366: Bathy_XP_007510256.1 ChlNC64A_IGS.gm_4_00538 Cre_12.g497101.t1 Czo_Cz02g17020.t1 MpusCCMP1545_196356 Ostta4221_3_69920_OT_ostta10g02300T0 Pse_000006F_consensus.g874.t1 Vocar_0006s0387.1
OG0005367: Bathy_XP_007510307.1 Cre_10.g448950.t1 Csub_66649 Czo_Cz13g08110.t1 Knitens_00185_0060_v1.1 MpusCCMP1545_163162 Ostta4221_3_72635_OT_ostta02g01180T0 Vocar_0070s0008.1

2. gff file of all species. tab delimited. chromosome identifier (format: specieschromosome#), gene identifier (format: species_gene#), start coordinate, end coordinate). 

Example:

Cr9	Cre_09.g388750.t1	2824325	2830821
Cr9	Cre_09.g388763.t1	3373489	3377999

3. window size (integer) -window size in number of genes 
4. minimum number of orthologs (integer) -filer for gene neighborhoods with proximal orthologous gene pairs from a minimum number of species 

---------------------------------------------------------------------------------

Outputs:

1. File.txt -full output of sliding windows
2. File.merged.txt -merged windows
3. File.ranked.txt -gene neighborhoods ranked giving weight to genes from smaller orthologous groups
4. File.ranked.uniq.txt 