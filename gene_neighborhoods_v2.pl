#! /usr/bin/perl 
use strict; use warnings;
use Data::Dumper;
use List::MoreUtils 'none';
use List::MoreUtils 'uniq';
use PerlX::Window;
use Time::HiRes qw( time );
use Getopt::Long;
use Array::Compare;


#finds gene clusters based on a set window size in number of genes and returns clusters with a set minimum number of homologs. 
# average base pair/gene is calculated for each chromosome, multiplied by window size and used to find clustered homologs on a sliding window basis.
#Usage: perl gene_neighborhoods.pl <OrthoMCL_output_file> <gff_file> <window_size_integer> <min_ortholog_species_number_integer>  OPTIONS: [--tabular][--stats] [--method default/cooccurring/select/default_clean/select_clean]


 ## collect error messages
my %WARNS;
local $SIG{__WARN__} = sub {
  my $message = shift;
  return if $WARNS{$message}++;
  logger('warning', $message);
};
sub logger {
    my ($level, $msg) = @_;
    if (open my $out, '>>', 'log.txt') {
        chomp $msg;
        print $out "$level - $msg\n";
    }
  }

if (scalar(@ARGV) < 4) {die "Usage: perl gene_neighborhoods.pl <OrthoMCL_output_file> <gff_file> <window_size_integer> <min_ortholog_species_number_integer>  OPTIONS: [--tabular][--stats] [--method default/cooccurring/select/default_clean/select_clean] \n";}

my $file= $ARGV[0]; # OrthoFinder legacy format or Orthomcl,input defaults  ******make sure no | in headers. headers must have species name before _.  eval 1e-5
my $gff_file= $ARGV[1]; # must be ordered list by gene coordinates: species-contig (must be unique from all species), gene start, gene end. 
my $window_size=$ARGV[2]; # change window size in number of genes
my $ortholog_number=$ARGV[3]; # min ortholog number of species
my $stats= 0;
my $tabular = 0;
my $method= "default";
GetOptions ("stats" => \$stats, "method=s" => \$method,"tabular" => \$tabular);

if (($method ne "default") && ($method ne "cooccur") && ($method ne "select") && ($method ne "default_clean" ) && ($method ne "select_clean")) {die "Method entered $method. Method options available: default, cooccur, select, or default_clean"};

my $select_spec;
if (($method eq "select")|| ($method eq "select_clean")){
	print "Enter species to select: ";
	$select_spec = <STDIN>;
	chomp ($select_spec);
	print "Selecting for $select_spec species.\n";
}

my $starttime=time();

#####read gff file and create gene-chromosome-coordinate hash pairs and calculate gene density per chromosome ######################

my %gff=();
my @gene_list;
my $chr_prev="none";
my $min=0;
my $max=0;
my $j=0;
my %gene_density_chr; # gene density per chromosome
my @good_chr;


open (GFF, $gff_file); 
#check for chromosomes larger than window size
	while (<GFF>){
		chomp;
		my @field = split ( /\t/, $_);
		if ($field[0] eq $chr_prev) {
			$j++;
			$max= $field[3];
		}
		else{
			if ($j > $window_size){
				push @good_chr, $chr_prev;
				$gene_density_chr{$chr_prev}=($max-$min)/$j*($window_size); # kb_chr/gene * window size
				#print "$j $min $max $gene_density_chr{$chr_prev} $chr_prev\n";
			}
			else {
				print "$chr_prev less than window size $window_size, omitted\n";
			}
			$j=1;
			$min=$field[2];	
		}
		if (eof){
			if ($j > $window_size){
				push @good_chr, $chr_prev;
				$gene_density_chr{$chr_prev}=($max-$min)/$j*($window_size);
			}
			else {
				print "$chr_prev less than window size $window_size, omitted\n";
			}
			$j=1;
			$min=$field[2];	
		}
		$chr_prev = $field[0];
	}

# read only chromosomes larger than window size

seek GFF, 0,0;
while (<GFF>){
	chomp;
	my @field = split ( /\t/, $_);
	if ($field[0] ~~ @good_chr) {
		$gff{$field[1]}= [$field[0],$field[2],$field[3]]; #id=chr,start,end
		push @gene_list, $field[1];	
	}
}


print "\nDone loading " . keys(%gff). " gene coordinates from ". scalar @good_chr ." chromosomes \n\n";

#######read orthogroups	#################################

my %ortholog_list;
my %OG;
my $OG_name;
my %score;
my %OG_spec;

open (IN, $file) or die "error ortholog file";
while (<IN>) { 
     chomp;
     my @field1 = (split /: /,$_,2);
     $OG_name= $field1[0];
     if (defined $field1[1]){
		 my @orthologs = split (/ /, $field1[1]);
		 if (scalar @orthologs > 1) { #only save orthogroups with at least two genes
			my $OG_score= 1/(scalar @orthologs);
			foreach my $ortholog (@orthologs) {
				if (length $ortholog > 1){ #ignore whitespaces
					my $spec = (split /_/,$ortholog)[0];
					$OG{$ortholog} = $OG_name;
					push @{$ortholog_list{$OG_name}}, $ortholog;
					$score{$ortholog} = $OG_score;
					if (($method eq "cooccur") || ($tabular == 1)) { #creates species list per OG for cooccurring comparison
						if (none {$_ eq $spec} @{$OG_spec{$OG_name}}){
							push @{$OG_spec{$OG_name}}, $spec;
						}
					}
				}
			}
		 }
	}
}

#print Dumper (%ortholog_list);
print "\nDone loading ". keys(%ortholog_list). " Ortholog groups \n\n";

############ Find proximal orthologous gene pairs in sliding window ################################

my %homolog_match=(); # actually same as window [#]/ref query
my %close_homologs=(); # array of close homologs to ref query
my @close_homolog_list=();
my @homolog_spec=();
my @uniq_homolog_spec=();
my $gene;
my @window_chr=();
my $k=0;
my $comp = Array::Compare->new;
my $header;
my $selected ="no";
my %spec_cluster_keep;
my @uniq_homolog_spec2;
my @homolog_spec2;

if ($method eq "default"){
	$header = "Win.$window_size.min.$ortholog_number";
}
if ($method eq "cooccur"){
	$header = "Win.$window_size.min.$ortholog_number.cooccur";
}
if ($method eq "select"){
	$header = "Win.$window_size.min.$ortholog_number.select.$select_spec";
}
if ($method eq "default_clean"){
	$header = "Win.$window_size.min.$ortholog_number.default-clean";
	open (TMP, ">$header.tmp");
}
if ($method eq "select_clean"){
	$header = "Win.$window_size.min.$ortholog_number.select-clean";
	open (TMP, ">$header.tmp");
}

open (OUT, ">$header.txt"); 

while (defined window @gene_list, $window_size){
	foreach $gene (@window) {
		push @window_chr, $gff{$gene}[0];
	}
	if (@window_chr == grep {$_ eq $gff{$window[0]}[0] }@window_chr){ #check if all in window are in same chr
		for(my $i = 0; $i < @window; $i++) {
			if (exists $ortholog_list{$OG{$window[$i]}}){
				foreach my $homolog (@{$ortholog_list{$OG{$window[$i]}}}){
					$homolog_match{$homolog}=$window[$i];
				}			
			}	
		}
		for my $key1 (keys %homolog_match){ #check references in window for pairs within range. key1 and key2 are genes from the same chromosome and same species
			my $key1_spec = (split /_/,$key1)[0]; 
			for my $key2 (keys %homolog_match){
				if ((exists $gff{$key1}) &&(exists $gff{$key2})) { #check chromosome data exists
					 if ($gff{$key1}[0] eq $gff{$key2}[0]) { #make sure on same chr 
						if (($gff{$key1}[1] < ($gff{$key2}[1] + $gene_density_chr{$gff{$key1}[0]})) && ($gff{$key1}[1] > ($gff{$key2}[1] - $gene_density_chr{$gff{$key1}[0]})) && ($homolog_match{$key2} !~m/$homolog_match{$key1}/)) { #make sure references within a range and not in the same homologous group
							my $homologmatch_key1_spec = (split /_/, $homolog_match{$key1}) [0];
							if ($homologmatch_key1_spec !~m/$key1_spec/){ # make sure homologs are not from the same species 
								if (none {$_ eq $key1} @close_homolog_list){
									push @{$close_homologs{$homolog_match{$key1}}}, $key1;
									push @close_homolog_list, $key1;
								}
							}
							
						}
					}
				}
			}
		}
		if ($method eq "default"){
			method_default();
		}
		if ($method eq "cooccur"){
			method_cooccur();
		}
		if ($method eq "select"){
			method_select();
		}
		if ($method eq "default_clean"){
			method_default_clean();
		}
		if ($method eq "select_clean"){
			method_select_clean();
		}
		$k=0;	
		@close_homolog_list=(); # reset window
		%close_homologs=();
		%homolog_match=();	
	}
	@window_chr = ();
}



print OUT "----------------\n";


print "Done finding gene pairs\n";


###sub method default

sub method_default {
	OUTER: for (my $j = 0; $j < @window; $j++) { #for each of the genes in window with homologs check if greater than minimum number of homolog species set then print window
		my $spec1 = (split /_/,$window[$j])[0];
		push @homolog_spec, $spec1;
		foreach (@{$close_homologs{$window[$j]}}) {
			my $spec = (split /_/,$_)[0];
			push @homolog_spec, $spec;
		}
		@uniq_homolog_spec = uniq @homolog_spec;
		if (scalar(@uniq_homolog_spec) >= $ortholog_number){
			$k++;
			if ($k==2){ #minimum two genes >min ortholog number
				print OUT "----------------\n";
				for (my $i = 0; $i < @window; $i++) {
					if (exists $close_homologs{$window[$i]}){
						print OUT "$window[$i]\t@{$close_homologs{$window[$i]}}\n";
					}
					else{
						print OUT "$window[$i]\n";
					}	
				}
				@homolog_spec=();
				last OUTER;
			}
	
		}
		@homolog_spec=();
	}
}

### sub method coocccurring 

sub method_cooccur {
	OUTER: for (my $j = 0; $j < @window; $j++) { #for each of the genes in window with homologs check if greater than minimum number of homolog species set then print window
		my $spec1 = (split /_/,$window[$j])[0];
		push @homolog_spec, $spec1;
		foreach (@{$close_homologs{$window[$j]}}) {
			my $spec = (split /_/,$_)[0];
			push @homolog_spec, $spec;
		}
		@uniq_homolog_spec = uniq @homolog_spec;

		if (scalar(@uniq_homolog_spec) >= $ortholog_number){
			if ($comp->perm(\@uniq_homolog_spec, \@{$OG_spec{$OG{$window[$j]}}}) ){ #count if all orthologs per OG are also cooccurring
				$k++;
			}
			if ($k==2) { #min two genes cooccurring and > min ortholog number
				print OUT "----------------\n";
				for (my $i = 0; $i < @window; $i++) {
					if (exists $close_homologs{$window[$i]}){
						print OUT "$window[$i]\t@{$close_homologs{$window[$i]}}\n";
					}
					else{
						print OUT "$window[$i]\n";
					}	
				}
				@homolog_spec=();
				last OUTER;
			}
		
		}
		@homolog_spec=();
	}
}

###sub method select

sub method_select {
	OUTER: for (my $j = 0; $j < @window; $j++) { #for each of the genes in window with homologs check if greater than minimum number of homolog species set then print window
		my $spec1 = (split /_/,$window[$j])[0];
		push @homolog_spec, $spec1;
		foreach (@{$close_homologs{$window[$j]}}) {
			my $spec = (split /_/,$_)[0];
			push @homolog_spec, $spec;
		}
		@uniq_homolog_spec = uniq @homolog_spec;

		if (scalar(@uniq_homolog_spec) >= $ortholog_number){
			if (grep {$_ eq $select_spec} @uniq_homolog_spec){
				$k++;
				if ($k==2) { #min two genes have selected species and > min ortholog number
					print OUT "----------------\n";
					for (my $i = 0; $i < @window; $i++) {
						if (exists $close_homologs{$window[$i]}){
							print OUT "$window[$i]\t@{$close_homologs{$window[$i]}}\n";
						}
						else{
							print OUT "$window[$i]\n";
						}	
					}
					@homolog_spec=();
					$selected="no";
					last OUTER;
				}
			}
		}
		@homolog_spec=();
	}
}


####### Clean-up methods
sub method_default_clean {
	for (my $j = 0; $j < @window; $j++) { #for each of the genes in window with homologs check if greater than minimum number of homolog species set then print window
		my $spec1 = (split /_/,$window[$j])[0];
		push @homolog_spec, $spec1;
		foreach (@{$close_homologs{$window[$j]}}) {
			my $spec = (split /_/,$_)[0];
			push @homolog_spec, $spec;
		}
		@uniq_homolog_spec = uniq @homolog_spec;
		
		if (scalar @uniq_homolog_spec >= $ortholog_number){
			$k++;
			foreach my $spec2 (@uniq_homolog_spec) {
				$spec_cluster_keep{$spec2}++; #count species per line with enough homologs
			}
		}
		@homolog_spec=();	
	}
	if ($k>=2){ #minimum two genes >min ortholog number
	print TMP "----------------\n";
		for (my $i = 0; $i < @window; $i++) {
			my $ref_spec = (split /_/,$window[$i])[0];
			push @homolog_spec2, $ref_spec;
			print TMP "$window[$i]\t";
			if (exists $close_homologs{$window[$i]}) {
				foreach (@{$close_homologs{$window[$i]}}) {
					my $spec4 = (split /_/,$_)[0];
					push @homolog_spec2, $spec4;
				}
				my @uniq_homolog_spec2 = uniq @homolog_spec2;
				if (scalar @uniq_homolog_spec2 >= $ortholog_number) { 
					#print "length good $window[$i]\t@{$close_homologs{$window[$i]}} \t @uniq_homolog_spec2\n";
					foreach my $gene (@{$close_homologs{$window[$i]}}) {
						my $spec3 = (split /_/,$gene)[0];
						#print "spec count $spec_cluster_keep{$spec3} $gene $spec3 k= $k\n";
						if ($spec_cluster_keep{$spec3} >= 2) { #remove genes that are not present in at least pairs in the highly conserved genes 
							print TMP "$gene ";
						}
					}
				}
			}
			print TMP "\n";	
			@homolog_spec2=();
		}
	}
	%spec_cluster_keep =();	
}

sub method_select_clean {
	for (my $j = 0; $j < @window; $j++) { #for each of the genes in window with homologs check if greater than minimum number of homolog species set then print window
		my $spec1 = (split /_/,$window[$j])[0];
		push @homolog_spec, $spec1;
		foreach (@{$close_homologs{$window[$j]}}) {
			my $spec = (split /_/,$_)[0];
			push @homolog_spec, $spec;
		}
		@uniq_homolog_spec = uniq @homolog_spec;
		
		if (scalar @uniq_homolog_spec >= $ortholog_number){
			if (grep {$_ eq $select_spec} @uniq_homolog_spec){
				$k++;
				foreach my $spec2 (@uniq_homolog_spec) {
					$spec_cluster_keep{$spec2}++; #count species per line with enough homologs
				}
			}
		}
		@homolog_spec=();	
	}
	if ($k>=2){ #minimum two genes >min ortholog number
	print TMP "----------------\n";
		for (my $i = 0; $i < @window; $i++) {
			my $ref_spec = (split /_/,$window[$i])[0];
			push @homolog_spec2, $ref_spec;
			print TMP "$window[$i]\t";
			if (exists $close_homologs{$window[$i]}) {
				foreach (@{$close_homologs{$window[$i]}}) {
					my $spec4 = (split /_/,$_)[0];
					push @homolog_spec2, $spec4;
				}
				my @uniq_homolog_spec2 = uniq @homolog_spec2;
				if (scalar @uniq_homolog_spec2 >= $ortholog_number) { 
					#print "length good $window[$i]\t@{$close_homologs{$window[$i]}} \t @uniq_homolog_spec2\n";
					foreach my $gene (@{$close_homologs{$window[$i]}}) {
						my $spec3 = (split /_/,$gene)[0];
						#print "spec count $spec_cluster_keep{$spec3} $gene $spec3 k= $k\n";
						if ($spec_cluster_keep{$spec3} >= 2) { #remove genes that are not present in at least pairs in the highly conserved genes 
							if (grep {$_ eq $select_spec} @uniq_homolog_spec2) {
								print TMP "$gene ";
							}
						}
					}
				}
			}
			print TMP "\n";	
			@homolog_spec2=();
			$selected="no";
		}
	}
	%spec_cluster_keep =();	
	
}



####################### Final Clean-up Step #############################
if (($method eq "select_clean") || ($method eq "default_clean")){
	my $m=0;
	my @line_ortho;
	my $cluster2;
	
	close TMP;
	open (TMPED, "<$header.tmp"); 
	while (<TMPED>){
		chomp;
		$cluster2 .= "$_\n";
		my @field1 = (split /\t/,$_,2);
		my $ref_spec2 = (split /_/,$field1[0])[0];
		push @line_ortho, $ref_spec2 ;
		if (defined $field1[1]){
			my @field2 = split (/ /, $field1[1]);
			foreach my $genes (@field2){
				my $spec5 = (split /_/,$genes)[0];
				push @line_ortho, $spec5;
			}
		}
		my @uniq_line = uniq @line_ortho;
		
		if (scalar @uniq_line >= $ortholog_number) { 
			$m++;
		}
		if ($_ eq '----------------') {
			if ($m >= 2){
			print OUT "$cluster2";
			}
			$m = 0;
			$cluster2 = "";
		}
	@line_ortho = ();
	
	
	}
	close TMPED;
	close OUT;
	
}


########################### Merge neighborhoods ############################
close OUT;

my $m=1;
my %all = ();
my $group=$m;
my %rest=();

open (MERGE, "<$header.txt");

while (<MERGE>) {
	chomp;
	if (!/^$/) { #if line not blank
		my @field = split (/\t/, $_);
		if ($_ eq '----------------'){
			$m++;
			$group=$m;
		}
		else {
			push @{$all{$group}},$field[0];
			if (defined $field[1]){
				my @results = split (/ /, $field[1]);
				foreach my $result (@results){
					push @{$rest{$field[0]}},$result;
				}
			}
			@{$rest{$field[0]}} = uniq @{$rest{$field[0]}};	
		}
	}
	
}

			
my %found;

foreach my $key ( sort keys %all ) {
    foreach my $element ( sort @{ $all{$key} } ) {
        if ( $found{$element} ) {
            #this element is present in another list, we merge.
            #print "$element found in $found{$element}\n";
            merge_and_delete( $found{$element}, $key );
            last;
        }
        else {
            #add this unique match to our index
            #print "$element -> $key\n";
            $found{$element} = $key;
        }
    }
}
	
open (MERGED, ">$header.merged.txt");
for my $key (sort keys %all){
	print MERGED "----------------\n";
	for my $ref (@{$all{$key}}) {
		print MERGED "$ref\t@{$rest{$ref}}\n"
	}
}
print MERGED "----------------\n";
close MERGED;
close MERGE;




sub merge_and_delete {
    my ( $first_key, $second_key ) = @_;
    #print "Merging $first_key with $second_key\n";

    #use hash to remove dupes.
    my %elements;
    foreach my $element ( @{ $all{$first_key} }, @{ $all{$second_key} } )
    {
        $elements{$element}++;
        #update index - don't want to point it to an array we're deleting
        $found{$element} = $first_key;
    }
    #sorting for neatness - you might want to do a numeric sort instead, 
    #as by default %all contains text elements. 
    $all{$first_key} = [ sort keys %elements ];
    delete $all{$second_key};

}

################################### Rank and Unique ##################################
# rank clusters giving weight to genes from smaller orthologous groups and get  unique clusters
open (RANK, "<$header.merged.txt");
open (RANKED, ">$header.ranked.txt");

my $cluster_score= 0;
my $cluster_size=0;
my $final_score =0;
my $cluster="";
my %ranked_cluster;
my @done;

while (<RANK>){
	chomp;
	my @field1 = (split /\t/,$_,2);
	if (defined $field1[1]){
		$cluster_size++;
		my @field2 = split (/ /, $field1[1]);
		foreach my $close_gene (@field2){
			if ($close_gene ne ""){
		 		$cluster_score= $cluster_score + $score{$close_gene};
		 	}
		}
	}
	if (($_ =~m/---/) && ($cluster_size >0)){
		$final_score= $cluster_score/$cluster_size*10;
		if( defined $ranked_cluster{$final_score}){ #make sure no repeated score
			my $taken= 1;
			until ($taken <1){
				$final_score= $final_score + .00000000000001;
				if( defined $ranked_cluster{$final_score}){
					$taken=1;
				}
				else {
					$taken=0;
				}
			}
		}
		$ranked_cluster{$final_score}=$cluster;
		$cluster_score =0;
		$cluster_size=0;
		$final_score=0;
		$cluster="";
	}
	$cluster .= "$_\n";
}

#print Dumper(%ranked_cluster);
close RANK;


for my $key (sort {$b <=> $a} keys %ranked_cluster){
	print RANKED "$ranked_cluster{$key}\ncluster score:$key\n";
}
close RANKED;

#find unique clusters
open (RANKED, "<$header.ranked.txt");
open (UNIQ, ">$header.uniq.txt");

my $print='yes';

while (<RANKED>){
	chomp;
	$cluster .= "$_\n";
	my @field1 = (split /\t/,$_,2);
	if (defined $field1[1]){
		my @field2 = split (/ /, $field1[1]);
		if ($field1[0] ~~ @done){
				$print ='no';
		}
		foreach my $close_gene (@field2){	
			push @done, $close_gene;
		}
		
	}
	if ($_ =~m/---/) {
		if ($print =~m/yes/){
			print UNIQ "$cluster\n";
		}
		$print = 'yes';
		$cluster="";
	}
	
}

close RANKED;
close UNIQ;
################Options###########################

################ Stats #################

#gets stats
if ($stats ==1) {


	open (UNIQED, "<$header.uniq.txt"); 
	open(STATS, ">$header.stats.txt");

	my @item_spec=();
	my @all_spec=();
	my %count_spec=();
	my $score=0;
	my @spec_order=();
	my $cluster_size =0;
	my $num_ortho=0;


 	#find all species for columns
	while (<UNIQED>){
		chomp;
		my @field1 = (split /\t/,$_,2);
		my @field2 = (split / /, $field1[1]);
		foreach my $items (@field2){
			my $species = (split /_/,$items)[0];
			if (none {$_ eq $species} @all_spec){
				push @all_spec, $species;
			}
		}
	}

	my @uniq_all_spec = uniq @all_spec;

	print STATS "Score\tMax_Number_orthologs_per_gene\tCluster_size\tTotal_number_species";
	foreach my $spec2 (@uniq_all_spec){
		print STATS "\t$spec2";
		push @spec_order, $spec2;
	}

	#read file and print only columns with >min_ortholog_number and get stats
	seek UNIQED, 0,0;
	while (<UNIQED>){
		chomp;
		my @field1 = (split /\t/,$_,2);
		my $species2 = (split /_/,$field1[0])[0];
		push @item_spec, $species2;
		if ($_ =~m/score/){
			$score = (split /:/,$_) [1];
		}
		if (length $field1[1] > 1){ #make sure not count just ref chromosome
			my @field2 = (split / /, $field1[1]);
			foreach my $items (@field2){
				my $species = (split /_/,$items)[0];
				push @item_spec, $species;
			}
			my @uniq_spec = uniq @item_spec;
			$cluster_size++;
			if (scalar @uniq_spec > $num_ortho){
				$num_ortho= scalar @uniq_spec;
			}
			foreach my $item2 (@item_spec){
				$count_spec{$item2}++;
			}
		
		}

		if ($_ eq '----------------'){
			my $cluster_spec = keys %count_spec;
			print STATS "\n$score\t$num_ortho\t$cluster_size\t$cluster_spec";
			for (my $j = 0; $j < @spec_order; $j++) {
				if (exists $count_spec{$spec_order[$j]}){
					print STATS "\t$count_spec{$spec_order[$j]}";
				}
				else {
					print STATS "\t0";
				}	
			}
			$num_ortho=0;
			%count_spec=();
			$score=0;
			$cluster_size=0;	
		}
		@item_spec =();
	}

	close UNIQED;
	close STATS;
}

#################Tabular output ################

if ($tabular == 1){
	open (UNIQED, "<$header.uniq.txt"); 
	open(TAB, ">$header.tabular.txt");

	
	my @all_spec;
	my @spec_order;
	my @line_genes;
	my %species_genes;

 	#find all species for column names
	while (<UNIQED>){
		chomp;
		my @field1 = (split /\t/,$_,2);
		if (defined $field1[1]){
			my @field2 = (split / /, $field1[1]);
			foreach my $items (@field2){
				my $species = (split /_/,$items)[0];
				if (none {$_ eq $species} @all_spec){
					push @all_spec, $species;
				}
			}
		}
	}

	my @uniq_all_spec = uniq @all_spec;
	print TAB "Ref_chromosome";
	foreach my $spec2 (@uniq_all_spec){
		print TAB "\t$spec2";
		push @spec_order, $spec2;
	}
	print TAB "\n";
	
	seek UNIQED, 0,0;
	while (<UNIQED>){
		chomp;
		my @field1 = (split /\t/,$_,2);
		push @line_genes, $field1[0];
		if (($_ eq '----------------') || (eof ==1)){
			print TAB "__________\t";
			for (my $j = 0; $j < @spec_order; $j++) {
				print TAB "__________\t";
			}
		}
		if (defined $field1[1]){
			print TAB "$field1[0]\t";
			my @field2 = (split / /, $field1[1]);
			push @line_genes, @field2;
			foreach my $gene (@line_genes){
				my $species2 = (split /_/,$gene)[0];
				push @{$species_genes{$species2}}, $gene;
			}
			for (my $j = 0; $j < @spec_order; $j++) {
				if (exists $species_genes{$spec_order[$j]}) {
					print TAB "@{$species_genes{$spec_order[$j]}}\t";
				}	
				elsif (grep {$_ eq $spec_order[$j]} @{$OG_spec{$OG{$field1[0]}}}) {
				
					print TAB "Present\t";
				}
				else {
					print TAB "Absent\t";	
				}
			}
		}
		%species_genes=();
		@line_genes=();
		print TAB "\n";
	}
	close TAB;
	
		## HTML table
	open (TABBED, "<$header.tabular.txt");
	open (HTML, ">$header.html");
	my $i=0;
	
	print HTML ' <style>.mytable{border-collapse:collapse; background-color: lightcoral;} .mytable td {border-right:2px solid white}.mytable td:nth-child(1) { background: gold; border-right: 5px solid black; }</style><table class= "mytable"><tbody>';
	while (<TABBED>){
		chomp;
		$i++;
		my @field = (split /\t/,$_);
		if ($i == 1){
			print HTML '<tr style = "background-color: gold;">';
			foreach my $col_name (@field){
				print HTML '<th style = "border-right:2px solid white";>'; print HTML "$col_name</th>";
			}
			print HTML "</tr>\n";
		}
		elsif (defined $field[1]){
			if ($field[1] =~ m/______/){
				print HTML '<tr style="border-top: 5px solid black;">';
			}
			else {
				foreach my $item (@field){
					chomp ($item);
					if ($item eq "Absent") {
						print HTML '<td style = "background-color: ghostwhite;">';print HTML "$item</td>";
					}
					elsif ($item eq "Present") {
						print HTML '<td style="background-color:mistyrose;">'; print HTML "$item</td>";
					}
					else{
						print HTML "<td>$item</td>";
					}
				}
				print HTML "</tr>\n";
			}
		}
		
	}
	print HTML "</tbody>\n</table>";

		use Data::Dumper;
	## HTML table
	open (TABBED, "<$header.tabular.txt");
	open (HTML, ">$header.html");
	open (HEAT, ">$header.heatmap.html");
	#open (TMP, ">$header.heatmap.tmp.txt");
	my $a=0;
	my @array=();
	my $tot_gene; #includes lines
	
	print HTML ' <style>.mytable{border-collapse:collapse; background-color: lightcoral;} .mytable td {border-right:4px solid white}.mytable td:nth-child(1) { background: powderblue; border-right: 5px solid black; }</style><table class= "mytable"><tbody>';
	while (<TABBED>){
		chomp;
		$a++;
		my @field = (split /\t/,$_);
		push @array, \@field;
		if ($a == 1){
			print HTML '<tr style = "background-color: powderblue;">';
			foreach my $col_name (@field){
				print HTML '<th style = "border-right:2px solid powderblue";>'; print HTML "$col_name</th>";
			}
			print HTML "</tr>\n";
		}
		elsif (defined $field[1]){
			$tot_gene++;
			if ($field[1] =~ m/______/){
				print HTML '<tr style="border-top: 5px solid black;">';
			}
			else {
				foreach my $item (@field){
					chomp ($item);
					if ($item eq "Absent") {
						print HTML '<td style = "background-color: white;">';print HTML "$item</td>";
					}
					elsif ($item eq "Present") {
						print HTML '<td style="background-color:mistyrose;">'; print HTML "$item</td>";
					}
					else{
						print HTML "<td>$item</td>";
					}
				}
				print HTML "</tr>\n";
			}
		}
		
	}
	print HTML "</tbody>\n</table>";
	

	
	
	########## for Heatmap
	use Array::Transpose;
	use Text::Table;
	
	my $col=0;
	my $row=0;	
	@array=transpose(\@array); 

	for my $ref (@array) {
		$row++;
		if ($row >= 2){
			for (my $j = 0; $j < scalar (@$ref); $j++) {
				$col++;
				if ($col >=2){
					if (@$ref[$j] =~ m/______/){
						@$ref[$j] =0;
						
					}
					elsif (@$ref[$j] eq "Absent"){
						@$ref[$j]=1;
					}
					elsif (@$ref[$j] eq "Present"){
						@$ref[$j]=2;
					}
					elsif (@$ref[$j] eq "") {
						delete @$ref[$j];
					}
					else {
						@$ref[$j]=3;
					}
				}
				
			}
			$col=0;
		}
	}

	
	
	####### For Clustering input
	open (TMP2, ">input_hierarchy.tmp.txt");
	open (TMP3, ">input_hierarchy2.tmp.txt");
	
	my @array2; #sum of presence/absence data per neighborhood
	my @newrow;
	my $sum=0;
	my $row=0;
	my $col=0;
	my %groups;
	my $group_counter=0;
	my %spec;
	my $spec_num=0;
	my %colgene;
	my $col1=0;

	
	for my $ref (@array) {
		if (defined $ref){
			$row++;
			if ($row == 1){
				for (my $c = 0; $c < scalar (@$ref); $c++) {
					if (defined @$ref[$c]){
						if (@$ref[$c] =~/_______/){
							$col1++;
						}
						elsif ( @$ref[$c] !~ /Ref_chromosome/) {
							push @{$colgene{$col1}}, @$ref[$c];
						}
					}
				}
			}
			if ($row >= 2){
				for (my $c = 0; $c < scalar (@$ref); $c++) {
					if (defined @$ref[$c]){
						$col++;
						if ($col==1){
							$spec{$spec_num}= @$ref[$c];
							$spec_num++;
						}
						if ($col >=2){
							$sum= $sum+@$ref[$c];
							if (@$ref[$c] == 0) {
								push @newrow, $sum;
								$sum=0;
								print TMP3 "$group_counter\t";
								$group_counter++;
							}
							else{
								push @{$groups{$group_counter}},@$ref[$c]; 
							}
						}
					}
				}
			}
		}			
		push @array2, [@newrow];
		@newrow=();
		$col=0;
		$sum=0;
		print TMP3 "\n";
	}
					
	my $tb = Text::Table->new;
	$tb-> load(@array2);
	print TMP2 $tb;
	
	close TMP2;
	close TMP3;
	
	####### Do Python script for hierarchical clustering


	system ("python hierarchical_cluster.py input_hierarchy.tmp.txt input_hierarchy2.tmp.txt");



	###reorder matrix

	open (ORDERSP, "<order_spec.tmp.txt");
	open (ORDER1, "<ordered.tmp.txt");
	open (ORDER2, ">ordered_genes.tmp.txt");
	open (ORDERCOL, "<order_col.tmp.txt");
	my @spec_order2;
	my @col_reorder;
	my $row=0;
	
	while (<ORDERSP>){
		chomp;
		push @spec_order2, $_;
	}
	
	while (<ORDERCOL>){
		chomp;
		push @col_reorder, $_;
	}
	
	print ORDER2 "Ref_chromosome\t";
	foreach my $colgene_id (@col_reorder){
		foreach my $colgene_name (@{$colgene{$colgene_id}}){
			print ORDER2 "$colgene_name\t";
		}
		print ORDER2 "______________________\t";
	}
	print ORDER2 "\n";
	while (<ORDER1>){
		chomp;
		my @field = (split /\t/,$_);
		print ORDER2 "$spec{$spec_order2[$row]}\t\0";
		foreach my $groupid (@field){
			foreach my $gene_ex ( @{$groups{$groupid}} ){
				print ORDER2 "$gene_ex\t";
			}
			print ORDER2 "0\t";
		}
		print ORDER2 "\n";
		$row++;
	}
	close ORDERSP;
	close ORDER1;
	close ORDER2;
	
	## Make HTML tabular output heatmap
	open (ORDER3, "<ordered_genes.tmp.txt");
	my @final_array;
	
	while (<ORDER3>) {
		chomp;
		my @field = (split /\t/,$_);
		push @final_array, (\@field);
	}
	print HEAT ' <style>.mytable{border-collapse:collapse; background-color: lightcoral;} .mytable th {border-bottom: 5px solid; background-color:powderblue} .mytable td {height: 34px ; border-top: 5px solid white; border-right: 5px solid white}.mytable td:nth-child(1) { background: powderblue; } th.rotate {height: 160px; white-space: nowrap;} th.rotate > div { transform: translate(25px, 51px) rotate(315deg);width: 34px;} th.rotate > div > span { padding: 0px 0px;} .wrapper {position: relative;} .scroller {margin-left: 141px;overflow-x: scroll;overflow-y: visible;padding-bottom: 5px;width: 2000px;} .mytable .headcol {left: 0;position: absolute;top: auto;width: 120px;} </style><div class="wrapper"> <div class="scroller"><table class= "mytable"><tbody>';
	print HEAT "\n";
	
	my $row=0;
	my $col=0;
	for my $ref (@final_array) {
		$row++;
		print HEAT "<tr>";
		if ($row == 1){
			print HEAT '<tr style = "background-color: white;">'; 
			foreach my $inner (@$ref) {
			
				if (length $inner >2){
					if ($inner =~ m/______/){
						print HEAT '<th class="rotate"><div><span>'; print HEAT "$inner</span></div></th>";
					}
					elsif ($inner eq "Ref_chromosome"){
						print HEAT '<th class = "headcol" class="rotate"><div><span>'; print HEAT "$inner</span></div></th>";
					}
					else{
						print HEAT '<th class="rotate"><div><span>';print HEAT  "$inner</span></div></th>";
					}
				}
				
			}
			print HEAT "</tr>\n";
			$col=0;
		}
		
		else{
			for my $inner (@$ref) {
				if (defined $inner){
					$col++;
					if ($col ==1){
						print HEAT '<td class="headcol" style=" border-top: 5px solid">'; print HEAT "$inner</td>";
					}
					elsif ($inner ==0){
							print HEAT '<td style="background-color:black; border-top: 5px solid; border-right: 5px solid"></td>';
					}
					else{
						if ($inner ==1){
							print HEAT '<td style = "background-color: white;"></td>';
						}
						elsif ($inner ==2){
							print HEAT '<td style="background-color:mistyrose;"></td>';
						}
						else{
							print HEAT "<td></td>";
						}
					}
				}
			}
			$col=0;
			print HEAT "</tr>\n";
		}
	}
	print HEAT "</tbody></table></div></div>";	
	#my $tb = Text::Table->new;
	$tb-> load(@final_array);
	#print TMP $tb;
	system ("rm *tmp.txt");
	
}

###############################################
my $end = time();
printf("runtime %.2f\n", $end - $starttime);
