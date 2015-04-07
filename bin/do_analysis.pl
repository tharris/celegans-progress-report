#!/usr/bin/perl

# Author: Todd Harris (todd@wormbase.org).
# You may freely modify and distribute this code
# provided you retain this citation.

use strict;
use Getopt::Long;
use Ace;

$|++;

my ($acedb,$host,$port,$help,$debug,);
GetOptions ('acedb=s'       => \$acedb,
            'host=s'        => \$host,
            'port=s'        => \$port,
            'help=s'        => \$help,
            'debug=s'       => \$debug,
           );

if ($help) {
  die <<END;
 Usage:  [ [options]

  Options:
  -acedb    Full path to acedb database to use (/usr/local/acedb/elegans)
    OR
  -host     Hostname of aceserver
  -port     Port of aceserver

  -debug    Boolean true to log debugging information

END
}


# Connect to ace and fetch the database version
print STDERR "\nConnecting to acedb: " . (($host) ? "$host:$port\n" : "$acedb\n");
$acedb        ||= '/usr/local/acedb/elegans';
my $DB = ($host)
  ? Ace->connect(-host=>$host,-port=>$port)
  : Ace->connect(-path=>$acedb);

my $version = $DB->status->{database}{version};

open OUT,">c_elegans.$version.gene_report.txt";

# How many genes have been genetically defined and cloned
#my @genes = $DB->fetch(-query=>qq{find Live_genes},-fill=>1);
my @genes = $DB->fetch(Gene => '*');

my %global_stats;
#my %permitted_classes = map { $_ => 1 } qw/let/;
foreach my $gene (@genes) {
    #  next unless defined $permitted_classes{$gene->Gene_class};
    next unless $gene->Species eq 'Caenorhabditis elegans';
    next unless $gene->Live(0);

    # Protein coding genes only
    # Uncloned genes or protein coding genes ONLY
    next if ($gene->Corresponding_transcript && !$gene->Corresponding_CDS);
    next unless ($gene->CGC_name || $gene->Sequence_name);

    $global_stats{total_genes}++;    

    my %gene_stats;
    
    my @variations = $gene->Variation;
    foreach (@variations) {
	next unless $_->Variation_type eq 'Allele';
	$gene_stats{number_of_alleles}++;
    }

    $gene_stats{genetically_defined} = defined $gene_stats{number_of_alleles} ? 'yes' : 'no';   
    $gene_stats{number_of_alleles}   ||= 0;

    $global_stats{genetically_defined}++                if ($gene_stats{number_of_alleles} > 0);
    $global_stats{genetically_defined_but_uncloned}++   if ($gene_stats{number_of_alleles} > 0 && !$gene->Corresponding_CDS);
    $global_stats{genetically_defined_and_cloned}++     if ($gene_stats{number_of_alleles} > 0 && $gene->Corresponding_CDS);
    

    my @rnai = $gene->RNAi_result;
    $gene_stats{total_rnai_experiments} = scalar @rnai;

    if (@rnai) {
	$global_stats{with_rnai_result}++;
    } else {
	$global_stats{lacking_rnai_result}++;
    }

    my (%unique_phenes);
    foreach (@rnai) {
	my $phenotypes = get_phenotypes($_);
	
	if (@$phenotypes == 0) {
	    $gene_stats{rnai_wildtype}++;
	} else {
	    $gene_stats{rnai_mutant}++;
	}
	
	foreach (@$phenotypes) {
	    $unique_phenes{$_}++;
	}
    }
    
    $global_stats{with_mutant_rnai}++ if defined $gene_stats{rnai_mutant};
    $global_stats{with_wt_rnai}++     if defined $gene_stats{rnai_wildtype};
    
    $gene_stats{phene_string} = join(';',keys %unique_phenes);
    $gene_stats{rnai_mutant}   ||= 0;
    
    # Does this gene have a: human ortholog?
    # Is the ortholog in OMIM?
    my @orthologs = $gene->Ortholog_other;
    foreach (@orthologs) {

	next unless $_ =~ /ENSEMBL:ENSP\d/;  # Human ortholog IDs.
	my @databases = $_->Database;
	my %dbs;
	foreach my $db (@databases) {
	    
	    my $name            = $db->Name || "$db";
	    my $description     = $db->Description;
	    my $url             = $db->URL;
	    my $url_constructor = $db->URL_constructor;
#	    my $email           = $db->Email;
	    my $remote_text     = $db->right(1);
	    
	    # Possibly multiple entries for a single DB                                                     
	    my @ids = map {
		my @types = $_->col;
		@types ? map { "$_" } @types : $_->right->name;
	    } $db->col;
	    if ($db =~ /omim/i || $name =~ /omim/i) {
		$gene_stats{omim} = join(";",@ids);
	    }
	}
    }
    
    if (defined $gene_stats{omim}) {
	$global_stats{with_omim}++;
	if ($gene_stats{genetically_defined} eq 'yes') {
	    $global_stats{with_omim_and_mutation}++;
	}
    }	

    # count unique phenes observed
    my $unique_phenes = scalar keys %unique_phenes;
    

#  if ($mut > 0 && $rnai > 0) {
    print OUT join("\t",$gene,$gene->CGC_name,$gene->Sequence_name || 'unknown',
		   $gene_stats{total_rnai_experiments},$gene_stats{rnai_mutant},
		   $gene_stats{genetically_defined},$gene_stats{number_of_alleles},
		   $gene_stats{omim}),"\n";

#		   $rnai_wildtype,$rnai_mutant,$phene_string,$unique_phenes,

#  }

#  print join("\t",$gene,$gene->Public_name,(scalar @rnai),$wt,$mut),"\n";

}


print "Total genes                            : $global_stats{total_genes}\n";
print "Genetically-defined                    : $global_stats{genetically_defined}\n";
print "Genetically-defined and cloned genes   : $global_stats{genetically_defined_and_cloned}\n";
print "Genetically-defined but uncloned genes : $global_stats{genetically_defined_but_uncloned}\n";
print "With OMIM orthologs  : $gobal_stats{with_omim}\n";
print "With OMIM orthologs and mutation : $global_stats{with_omim_and_mutation}\n";    
print "Lacking RNAi results : $global_stats{lacking_rnai_result}\n";
print "   With RNAi results : $global_stats{with_rnai_result}\n";
print "   With mutant phene : $global_stats{with_mutant_rnai}\n";   # at least one mutant RNAi phenotype


#print "Proportion of genetically defined genes with a WT RNAi phenotype: ",
#  (($stats{with_wt_rnai_only} / $stats{rnai_result} ) * 100),"\n";

#print "Proportion of genetically defined genes with a mutant RNAi phenotype: ",
#  ((($stats{with_mutant_rnai_only} + $stats{with_both}) / $stats{rnai_result} ) * 100), "\n";




sub get_phenotypes {
    my $rnai = shift;

    my @phenotypes = $rnai->Phenotype;
    my @phenes;
    foreach my $phenotype (@phenotypes) {
	my $phenotype_desc = $phenotype->Primary_name;
	my $string = "$phenotype:$phenotype_desc";
	push @phenes,$string;
    }
    return \@phenes;
}
