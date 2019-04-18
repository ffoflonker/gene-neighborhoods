#! /usr/bin/perl 
use strict; use warnings;
use Data::Dumper;
use List::MoreUtils 'none';
use List::MoreUtils 'uniq';
use PerlX::Window;
use Time::HiRes qw( time );
use Getopt::Long;
use Array::Compare;

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

#finds gene clusters based on a set window size in number of genes and returns clusters with a set minimum number of homologs. 
# average base pair/gene is calculated for each chromosome, multiplied by window size and used to find clustered homologs on a sliding window basis.
#
# usage: perl gene_neighborhoods.pl Orthofinder_file.txt file.gff window_size min_ortholog_species_number 
# next step: perl rank_gene_neighborhoods.pl gene_neighborhoods.out.txt Orthofinder_file.txt



if (scalar(@ARGV) < 4) {die "Usage: perl gene_neighborhoods.pl <OrthoMCL_output_file> <gff_file> <window_size_integer> <min_ortholog_species_number_integer>  OPTIONS: [--stats] [--method default/cooccurring/select] \n";}

my $file= $ARGV[0]; # OrthoFinder legacy format or Orthomcl,input defaults  ******make sure no | in headers. headers must have species name before _.  eval 1e-5
my $gff_file= $ARGV[1]; # must be ordered list by gene coordinates: species-contig (must be unique from all species), gene start, gene end. 
my $window_size=$ARGV[2]; # change window size in number of genes
my $ortholog_number=$ARGV[3]; # min ortholog number of species
my $stats= 0;
my $method= "default";
GetOptions ("stats" => \$stats, "method=s" => \$method,);

if (($method ne "default") && ($method ne "cooccur") && ($method ne "select") ) {die "Method entered $method. Method options available: default, cooccur, or select"};

my $select_spec;
if ($method eq "select"){
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
     if (exists $field1[1]){
		 my @orthologs = split (/ /, $field1[1]);
		 if (scalar @orthologs > 1) { #only save orthogroups with at least two genes
			my $OG_score= 1/(scalar @orthologs);
			foreach my $ortholog (@orthologs) {
				if (length $ortholog > 1){ #ignore whitespaces
					my $spec = (split /_/,$ortholog)[0];
					$OG{$ortholog} = $OG_name;
					push @{$ortholog_list{$OG_name}}, $ortholog;
					$score{$ortholog} = $OG_score;
					if ($method eq "cooccur"){ #creates species list per OG for cooccurring comparison
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

if ($method eq "default"){
	$header = "Win.$window_size.min.$ortholog_number";
}
if ($method eq "cooccur"){
	$header = "Win.$window_size.min.$ortholog_number.cooccur";
}
if ($method eq "select"){
	$header = "Win.$window_size.min.$ortholog_number.select";
}
open (TMP, ">$header.txt"); 

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
		OUTER: for (my $j = 0; $j < @window; $j++) { #for each of the genes in window with homologs check if greater than minimum number of homolog species set then print window
			my $spec1 = (split /_/,$window[$j])[0];
			push @homolog_spec, $spec1;
			foreach (@{$close_homologs{$window[$j]}}) {
				my $spec = (split /_/,$_)[0];
				push @homolog_spec, $spec;
			}
			@uniq_homolog_spec = uniq @homolog_spec;
			if ($method eq "default"){
				method_default();
			}
			if ($method eq "cooccur"){
				method_cooccur();
			}
			if ($method eq "select"){
				method_select();
			}
			@homolog_spec=();
		}
		$k=0;	
		@close_homolog_list=(); # reset window
		%close_homologs=();
		%homolog_match=();	
	}
	@window_chr = ();
}



print TMP "----------------\n";
close TMP;

print "Done finding gene pairs\n";


###sub method default

sub method_default {

	if (scalar(@uniq_homolog_spec) >= $ortholog_number){
		$k++;
		if ($k==2){ #minimum two genes >min ortholog number
			print TMP "----------------\n";
			for (my $i = 0; $i < @window; $i++) {
				if (exists $close_homologs{$window[$i]}){
					print TMP "$window[$i]\t@{$close_homologs{$window[$i]}}\n";
				}
				else{
					print TMP "$window[$i]\n";
				}	
			}
			@homolog_spec=();
			last OUTER;
		}
	
	}
}

### sub method coocccurring 

sub method_cooccur {
	

	if (scalar(@uniq_homolog_spec) >= $ortholog_number){
		if ($comp->perm(\@uniq_homolog_spec, \@{$OG_spec{$OG{$window[$j]}}}) ){ #count if all orthologs per OG are also cooccurring
			$k++;
		}
		if ($k==2) { #min two genes cooccurring and > min ortholog number
			print TMP "----------------\n";
			for (my $i = 0; $i < @window; $i++) {
				if (exists $close_homologs{$window[$i]}){
					print TMP "$window[$i]\t@{$close_homologs{$window[$i]}}\n";
				}
				else{
					print TMP "$window[$i]\n";
				}	
			}
			@homolog_spec=();
			last OUTER;
		}
		
	}

}

###sub method select

sub method_select {
	if (scalar(@uniq_homolog_spec) >= $ortholog_number){
		if (grep {$_ eq $select_spec} @uniq_homolog_spec){
			$k++;
			if ($k==2) { #min two genes have selected species and > min ortholog number
				print TMP "----------------\n";
				for (my $i = 0; $i < @window; $i++) {
					if (exists $close_homologs{$window[$i]}){
						print TMP "$window[$i]\t@{$close_homologs{$window[$i]}}\n";
					}
					else{
						print TMP "$window[$i]\n";
					}	
				}
				@homolog_spec=();
				$selected="no";
				last OUTER;
			}
		}
	}
}

############## Merge neighborhoods ############################


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
			if (exists $field[1]){
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

################Options###########################

if ($stats == 1){
	clean_stats();
}

################Cleanup and Stats#################

#filters for highly conserved genes >min_ortholog number and gets stats
sub clean_stats {

	open (CLEAN, ">$header.clean.txt");
	open (UNIQED, "<$header.ranked.uniq.txt"); 
	open(STATS, ">$header.stats.txt");

	my @item_spec;
	my @all_spec;
	my %count_spec;
	my $score=0;
	my @spec_order;
	my $cluster_size =0;
	my $num_ortho;

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

	print STATS "Score\tNumber_orthologs\tCluster_size";
	foreach my $spec2 (@uniq_all_spec){
		print STATS "\t$spec2";
		push @spec_order, $spec2;
	}

	#read file and print only columns with >min_ortholog_number and get stats
	seek UNIQED, 0,0;
	while (<UNIQED>){
		chomp;
		my @field1 = (split /\t/,$_,2);
		push @item_spec, $field1[0];
		if ($_ =~m/score/){
			$score = (split /:/,$_) [1];
		}
		if (exists $field1[1]){
			my @field2 = (split / /, $field1[1]);
			foreach my $items (@field2){
				my $species = (split /_/,$items)[0];
				push @item_spec, $species;
			}
			my @uniq_spec = uniq @item_spec;
		
			if (scalar @uniq_spec>= $ortholog_number) {
				print CLEAN "$_\n";
				$cluster_size++;
				if (scalar @uniq_spec > $num_ortho){
					$num_ortho= scalar @uniq_spec;
				}
				foreach my $item2 (@item_spec){
					$count_spec{$item2}++;
				}
			}
			else {
				print CLEAN "$field1[0]\n";
			}
		}
		else {
			print CLEAN "$field1[0]\n";
		}
		if ($_ eq '----------------'){
			print STATS "\n$score\t$num_ortho\t$cluster_size";
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
	close CLEAN;
	close UNIQED;
	close STATS;
}

###############################################
my $end = time();
printf("runtime %.2f\n", $end - $starttime);
