#!/usr/bin/perl -w

# Modules and packages
use strict;
use warnings;
use Getopt::Long;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::ApiVersion;
use Bio::EnsEMBL::Utils::Exception qw(throw);
use PatternModule::CleaningVariants;
use PatternModule::GeneAnno;
use PatternModule::RegAnno;
use List::Util qw(min);
use List::Util qw(max);

#########################
##### Documentation #####
#########################

=head1 NAME

Transcription Factor Binding Site Annotator

=cut

=head1 DESCRIPTION

This script is designed to annotate the effects of variants at gene level, and the position
and effect of variants in transcription factor binding motifs. This version of the script
scans all adjoining sequences for alternative TFBSs for both SNVs and INDELs

=cut

=head1 LICENSE

TFBS Annotator, Version 1.0
Copyright(c) 2016 Esben Eickhardt & Francesco Lescai
All Rights Reserved.

The authors are with Aarhus University, Department of Biomedicine.

Permission is hereby granted, free of charge, to any person obtaining a copy of this 
software and associated documentation files (the "Software"), to deal in the Software 
without restriction for educational and research purposes only, including without 
limitation the rights to use, copy, modify, merge, publish and distribute the Software, 
and to permit persons to whom the Software is furnished to do so, provided that this 
copyright notice and the original authors names appear on all copies and supporting 
documentation, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

This software shall not be used, rewritten, or adapted as the basis of a commercial
software or hardware product without first obtaining permission of the
authors. The authors make no representations about the suitability of
this software for any purpose.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

=cut

=head1 SYNOPSIS

=head3 USAGE

perl TFBSAnnotator.pl -in vcf_file_input.vcf -out vcf_file_output.vcf

=head3 INPUT

The script takes as input a VCF file with the same annotations a the 1000 genomes project
phase 3. For variants that are multi allelic all consequences on TFBSs are returned together,
and are not separated by allele. When the alternative alleles are "complex" (i.e. SNPs, INS 
or DEL in the formats AA/AG, AT/ATT, ATT/AT), the variants are further "cleaned" by a module 
in the script, in order to correct for sequence overlaps and other format issues (See 
CleaningVariants.pm).

=head3 RETURN		
												
The output is a txt-file with effects of variants at gene-level and TFBS-level. The script 
fixes the format of those variants present in the "complex" format (i.e. resulting from 
multisample calling) as explained previously. See detailed documentation on the cleaning 
methods (CleaningVariants.pm).

=head3 OTHER		
												
The script requires that the ENSEMBL API is installed on your system, as it is from this
database the annotations are pulled. The script has been tested only on API version 75 with
the four following databases:

homo_sapiens_core_75_37
homo_sapiens_funcgen_75_37
homo_sapiens_variation_75_37
ensembl_compara_75

For details on the installation of databases and the API see ENSEMBL's official website

=head1 FEEDBACK

Esben Eickhardt esbeneickhardt@biomed.au.dk

=head1 AUTHORS

Esben Eickhardt - Email E<lt>esbeneickhardt@biomed.au.dkE<gt>

Francesco Lescai - Email E<lt>lescai@biomed.au.dkE<gt>

=cut

##############################################
##### Connecting to the ENSEMBL database #####
##############################################

my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org',
    -user => 'anonymous'
);

# Getting adaptors
my $cs_adaptor = $registry->get_adaptor( 'Human', 'Core', 'CoordSystem' );
my $slice_adaptor = $registry->get_adaptor( 'human', 'core', 'slice');
my $variant_feature_adaptor = $registry->get_adaptor('human', 'variation', 'variationfeature');
my $transcript_adaptor = $registry->get_adaptor( 'Human', 'Core', 'Transcript' );
my $gene_adaptor  = $registry->get_adaptor( 'Human', 'Core', 'Gene' );
my $regfeat_adaptor = $registry->get_adaptor('Human', 'Funcgen', 'RegulatoryFeature');
my $motifFeature_adaptor = $registry->get_adaptor('Human', 'funcgen', 'motiffeature');
my $mlss_adaptor = $registry->get_adaptor("Multi", "compara", "MethodLinkSpeciesSet");
my $csscore_adaptor = $registry->get_adaptor("Multi", 'compara', 'ConservationScore');

# Get method_link_species_set object for gerp conservation scores for mammals
my $mlss = $mlss_adaptor->fetch_by_method_link_type_species_set_name("GERP_CONSERVATION_SCORE", "mammals");
throw("Unable to find method_link_species_set") if (!defined($mlss));

###################
##### Options #####
###################

# Creating options system allowing one to decide what needs to be annotated
my ($geneint, $tfbsint, $inputFile, $outputFile) = (0,0,0,0);

# Setting usage message
my $usage="\nUsage: 
	PatternAnnotator.pl -in vcf_file_input.vcf -out vcf_file_output.vcf
optional parameters:
	-gene_interval integer for chosing a custom interval for finding nearest gene
	
	\nPlease try again.\n\n\n";

# Input and output are the only mandatory options
die $usage unless GetOptions(
    'gene_interval:i' 	=> \$geneint,
    'in:s' => \$inputFile,
    'out:s' => \$outputFile
	)
	&& $inputFile && $outputFile;

$geneint ||= 5000;
$tfbsint ||= 5000;

# Opening file and creating output file
open (FILE, $inputFile) or die $!;
open (OUTPUT, "> $outputFile");

# Printing out what coordinate version is used:
my $cs = $cs_adaptor->fetch_by_name('chromosome');
print OUTPUT "##Coordinate System:", "\n";
printf OUTPUT "##%s %s\n", $cs->name(), $cs->version();

# Printing our what API version is used:
my $soft_version = software_version(); 
print OUTPUT "##ENSEMBL API Version:", "\n";
print OUTPUT "##".$soft_version."\n";

# Printing what option is used
print OUTPUT "##Intervals Around Variants for Finding Closest Features:", "\n";
print OUTPUT "##"."Genes: ".$geneint." bases"."\n";
print OUTPUT "##\n";
print OUTPUT "#Chromosome	Position	rsNumber	Ref	Alt	VariantType	MAF	NearestGene MostSevereConsequence	GeneLength	TSS	TFBS1_NAME	TFBS1_MATRIX	TFBS1_POS	TFBS1_LENGTH	TFBS1_DIST_TO_START	TFBS1_SIZE_OF_INDEL	TFBS1_ORIG	TFBS1_NEW	TFBS1_INFO_CONTENT	TFBS1_CONS_SCORE	TFBS2_NAME TFBS2_MATRIX	TFBS2_POS TFBS2_LENGTH	TFBS2_DIST_TO_START	TFBS2_SIZE_OF_INDEL	TFBS2_ORIG	TFBS2_NEW	TFBS2_INFO_CONTENT	TFBS2_CONS_SCORE	TFBS3_NAME	TFBS3_MATRIX	TFBS3_POS TFBS3_LENGTH	TFBS3_DIST_TO_START	TFBS3_SIZE_OF_INDEL	TFBS3_ORIG	TFBS3_NEW	TFBS3_INFO_CONTENT	TFBS3_CONS_SCORE	TFBS4_NAME	TFBS4_MATRIX TFBS4_POS TFBS4_LENGTH	TFBS4_DIST_TO_START	TFBS4_SIZE_OF_INDEL	TFBS4_ORIG	TFBS4_NEW	TFBS4_INFO_CONTENT	TFBS4_CONS_SCORE TFBS5_NAME	TFBS5_MATRIX	TFBS5_POS TFBS5_LENGTH	TFBS5_DIST_TO_START	TFBS5_SIZE_OF_INDEL	TFBS5_ORIG	TFBS5_NEW	TFBS5_INFO_CONTENT	TFBS5_CONS_SCORE	TFBS6_NAME	TFBS6_MATRIX	TFBS6_POS TFBS6_LENGTH	TFBS6_DIST_TO_START	TFBS6_SIZE_OF_INDEL	TFBS6_ORIG	TFBS6_NEW	TFBS6_INFO_CONTENT	TFBS6_CONS_SCORE	TFBS7_NAME	TFBS7_MATRIX	TFBS7_POS TFBS7_LENGTH	TFBS7_DIST_TO_START	TFBS7_SIZE_OF_INDEL	TFBS7_ORIG	TFBS7_NEW	TFBS7_INFO_CONTENT	TFBS7_CONS_SCORE	TFBS8_NAME	TFBS8_MATRIX	TFBS8_POS TFBS8_LENGTH	TFBS8_DIST_TO_START	TFBS8_SIZE_OF_INDEL	TFBS8_ORIG	TFBS8_NEW	TFBS8_INFO_CONTENT	TFBS8_CONS_SCORE TFBS9_NAME	TFBS9_MATRIX	TFBS7_POS TFBS9_LENGTH	TFBS9_DIST_TO_START	TFBS9_SIZE_OF_INDEL	TFBS9_ORIG	TFBS9_NEW	TFBS9_INFO_CONTENT	TFBS9_CONS_SCORE TFBS10_NAME	TFBS10_MATRIX	TFBS10_POS TFBS10_LENGTH	TFBS10_DIST_TO_START	TFBS10_SIZE_OF_INDEL	TFBS10_ORIG	TFBS10_NEW	TFBS10_INFO_CONTENT TFBS10_CONS_SCORE";
print OUTPUT "\n";

###############################################
##### Going through the file line by line #####
###############################################

while (<FILE>) {
    chomp;
    my $line = $_;    
    my $annotationString ="";
    my $vf = ();
    
    if ($line !~ /^#/) {
		my @col = split("\t", $line);
		my $slice = $slice_adaptor->fetch_by_region('chromosome', $col[0]);
				
		# A VariationFeature is created for variants with a single alternative allele
    	if ($col[4] !~ /,/) {
    	
			############################################################
			##### Cleaning variant and creating a VariationFeature #####
			############################################################
    		
    		# Variants in the form <.....> are not cleaned
    		if ($col[3] !~ "<.*>" && $col[4] !~ "<.*>") {
    		
				# Clean SNPs
				my @cleaningInfoSNP = PatternModule::CleaningVariants::cleanSNP($col[3], $col[4]);
				if (scalar @cleaningInfoSNP == 3) {
					print STDOUT "Cleaning complex variant:", "\n";
					print STDOUT "\tFrom: " . $col[3] . "/" . $col[4] . "\n";
					print STDOUT "\tStart position: " . $col[1] . "\n";
					$col[3] = $cleaningInfoSNP[0];
					$col[4] = $cleaningInfoSNP[1];
					$col[1] = $col[1] + $cleaningInfoSNP[2];
					print STDOUT "\tTo: " . $col[3] . "/" . $col[4] . "\n";
					print STDOUT "\t" . "Start position: " . $col[1] . "\n";
				}
			
				# Clean Deletions
				my @cleaningInfoDELETION = PatternModule::CleaningVariants::cleanDELETIONS($col[3], $col[4]);
				if (scalar @cleaningInfoDELETION == 3) {
					print STDOUT "Cleaning deletion variant:", "\n";
					print STDOUT "\tFrom: " . $col[3] . "/" . $col[4] . "\n";
					print STDOUT "\tStart position: " . $col[1] . "\n";
					$col[3] = $cleaningInfoDELETION[0];
					$col[4] = $cleaningInfoDELETION[1]; 
					$col[1] = $col[1] + $cleaningInfoDELETION[2];	
					print STDOUT "\tTo: " . $col[3] . "/" . $col[4] . "\n";	
					print STDOUT "\t" . "Start position: " . $col[1] . "\n";
				}
			
				# Clean Insertions
				my @cleaningInfoINSERTION = PatternModule::CleaningVariants::cleanINSERTIONS($col[3], $col[4]);
				if (scalar @cleaningInfoINSERTION == 3) {
					print STDOUT "Cleaning insertion variant:", "\n";
					print STDOUT "\tFrom: " . $col[3] . "/" . $col[4] . "\n";
					print STDOUT "\tStart position: " . $col[1] . "\n";
					$col[3] = $cleaningInfoINSERTION[0];
					$col[4] = $cleaningInfoINSERTION[1]; 
					$col[1] = $col[1] + $cleaningInfoINSERTION[2];	
					print STDOUT "\tTo: " . $col[3] . "/" . $col[4] . "\n";	
					print STDOUT "\t" . "Start position: " . $col[1] . "\n";
				}
				
				# Creates variationfeature
				$vf = PatternModule::GeneAnno::createVariationFeature($col[1], $col[3], $col[4], $slice, $variant_feature_adaptor);
			}
		}
		
		# Variation features not created for variants in the form <....>
    	if ($col[3] !~ "<.*>" && $col[4] !~ "<.*>") {
    	
			# Creating one VariationFeature for each multi allelic variant
			if ($col[4] =~ /,/) {
				$vf = PatternModule::GeneAnno::createMultiAllelicVariationFeature($col[1], $col[3], $col[4], $slice, $variant_feature_adaptor);
			}
			
		}

			#############################################
			##### Annotation of variant information #####
			#############################################
			
			# Annotating chromosome, position, rs-number, reference and alternative sequences directly from the vcf-file
			$annotationString .= 
				$col[0]."\t".
				$col[1]."\t".
				$col[2]."\t".
				$col[3]."\t";
			
			# Assigning the variant type for variants that are not of in the format "<.*>"
			if ($col[3] !~ "<.*>" && $col[4] !~ "<.*>" && $col[4] !~ /,/) {
				$annotationString .= $col[4]."\t"; 
			
				# Testing variant type using length of reference and alternative alleles
				if (length($col[3]) == 1 && length($col[4]) == 1) {
					$annotationString .= "SNP"."\t";
				}
			
				if (length($col[3]) == length($col[4]) && length($col[3]) > 1) {
					$annotationString .= "SUBSTITUTION"."\t";
				}
				
				if (length($col[3]) > length($col[4])) {
					$annotationString .= "DELETION"."\t";
				}
			
				if (length($col[3]) < length($col[4])) {
					$annotationString .= "INSERTION"."\t";
				}
				
			}
			
			# Variants where the alternative allele contains a comma are assigned as multiallelic
			if ($col[4] =~ /,/)  {
				my $str = $col[4];
				$str =~ s/<+//g;
				$str =~ s/>+//g;
				$annotationString .= $str."\t"."MULTIALLELIC"."\t";
			}

			# Variants of format "<.*>" are assigned as the variant type contained within the brackets
			if ($col[4] =~ "<.*>" && $col[4] !~ /,/) {
				my $str = $col[4];
				$str =~ s/<+//g;
				$str =~ s/>+//g;
				$annotationString .= $str."\t".$str."\t";
			}

			#############################
			##### Annotation of MAF #####
			#############################
			
			my $MAF = 0;
			
			# The allele frequencies are pulled from that info-section of the vcf-file, 
			# and if the variant is multi allelic the allele frequencies are added together
			
			# VCF-file for Y-chromosome has a different structure than that for the other chromosomes
			if ($col[0] eq "Y") {
				my @infocol = split(";", $col[7]);
				my @infocolsplit1 = split("=", $infocol[2]);
				my @infocolsplit2 = split(",", $infocolsplit1[1]);
			
				# If only a single allele
				if (scalar @infocolsplit2 == 1) {
					$MAF = $infocolsplit2[0];
				}
			
				# If multiallelic variant
				if (scalar @infocolsplit2 > 1) {
					$MAF += $_ for @infocolsplit2;
				}				
			}
			
			# VCF-file for Y-chromosome has a different structure than that for the other chromosomes
			if ($col[0] ne "Y") {
				my @infocol = split(";", $col[7]);
				my @infocolsplit1 = split("=", $infocol[1]);
				my @infocolsplit2 = split(",", $infocolsplit1[1]);
			
				# If only a single allele
				if (scalar @infocolsplit2 == 1) {
					$MAF = $infocolsplit2[0];
				}
			
				# If multiallelic variant
				if (scalar @infocolsplit2 > 1) {
					$MAF += $_ for @infocolsplit2;
				}
			}
			
			$annotationString .= $MAF."\t";

			#############################################################
			##### Annotation of gene regions related to the variant #####
			#############################################################
			
			# Annotating variants not of the type <.....>
			if ($col[3] !~ "<.*>" && $col[4] !~ "<.*>") {
			
				my %variantInformation = PatternModule::GeneAnno::createVariantProteinCodingHash();
				$variantInformation{'ENSPCANNO_GENE_LENGTH'} = "NA";
				$variantInformation{'ENSPCANNO_TRANSCRIPT_START_SITE'} = "NA";
			
				if ($vf) {
					# If overlaps are found with genes, the first one is annotated here
					my ($one, $two) = PatternModule::GeneAnno::doesVariantOverlapWithAnyGene($vf);
					$variantInformation{'ENSPCANNO_VARIANT_OVERLAPS_WITH_GENE'} = $one;
					$variantInformation{'ENSPCANNO_GENE'} = $two;
				
					# Most severely hit transcript is found together with its consequence
					my ($three, $four) = PatternModule::GeneAnno::findMostSeverelyHitTranscript($vf,$transcript_adaptor, "ANY");
					$variantInformation{'ENSPCANNO_TRANSCRIPT'} = $three;
					$variantInformation{'ENSPCANNO_TRANSCRIPT_CONSEQUENCE'} = $four;
				
					# If a transcript is found, the gene belonging to the most severely hit transcript is annotated
					if ($three ne "NA") {
						$variantInformation{'ENSPCANNO_GENE'} = $transcript_adaptor->fetch_by_stable_id($three)->get_Gene()->display_xref->display_id();
					}
				
					# If no overlapping genes are found, the nearest gene is found within 5000 bases as default
					if ($variantInformation{'ENSPCANNO_GENE'} eq "NA") {
						my ($five, $six) = PatternModule::GeneAnno::getNearestGene($vf, $slice_adaptor, "ANY", $geneint);
						$variantInformation{'ENSPCANNO_GENE'} = $five;
						$variantInformation{'ENSPCANNO_TRANSCRIPT_CONSEQUENCE'} = $six;
					}
										
					# Adding length of gene
					if ($variantInformation{'ENSPCANNO_GENE'} ne "NA") {
						my $gene = $gene_adaptor->fetch_by_display_label($variantInformation{'ENSPCANNO_GENE'});
						$variantInformation{'ENSPCANNO_GENE_LENGTH'} = $gene->length();
						$variantInformation{'ENSPCANNO_TRANSCRIPT_START_SITE'} = $gene->canonical_transcript->seq_region_start();	
					}
				
					# Adding information to the annotation
					$annotationString .= 
						$variantInformation{'ENSPCANNO_GENE'}."\t".
						$variantInformation{'ENSPCANNO_TRANSCRIPT_CONSEQUENCE'}."\t".
						$variantInformation{'ENSPCANNO_GENE_LENGTH'}."\t".
						$variantInformation{'ENSPCANNO_TRANSCRIPT_START_SITE'}."\t";
				}
			
			}
			
			# If variants are of type <....> they are annotated whether they overlap with a gene and what gene
			if ($col[3] =~ "<.*>" | $col[4] =~ "<.*>") {
				# Creates a "fake" variation feature to access variationfeature methods
				$vf = PatternModule::GeneAnno::createVariationFeature($col[1], $col[3], $col[3], $slice, $variant_feature_adaptor);
								
				# Finds overlapping and closest genes
				my %variantInformation = PatternModule::GeneAnno::createVariantProteinCodingHash();
				$variantInformation{'ENSPCANNO_GENE_LENGTH'} = "NA";
				$variantInformation{'ENSPCANNO_TRANSCRIPT_START_SITE'} = "NA";
				
				# Overlapping gene
				my ($one, $two) = PatternModule::GeneAnno::doesVariantOverlapWithAnyGene($vf);
				$variantInformation{'ENSPCANNO_VARIANT_OVERLAPS_WITH_GENE'} = $one;
				$variantInformation{'ENSPCANNO_GENE'} = $two;
			
				# Most severely hit transcript is found together with its consequence
				my ($three, $four) = PatternModule::GeneAnno::findMostSeverelyHitTranscript($vf,$transcript_adaptor, "ANY");
				$variantInformation{'ENSPCANNO_TRANSCRIPT'} = $three;
				$variantInformation{'ENSPCANNO_TRANSCRIPT_CONSEQUENCE'} = $four;
			
				# If a transcript is found, the gene belonging to the most severely hit transcript is annotated
				if ($three ne "NA") {
					$variantInformation{'ENSPCANNO_GENE'} = $transcript_adaptor->fetch_by_stable_id($three)->get_Gene()->display_xref->display_id();
				}
			
				# If no overlapping genes are found, the nearest gene is found within 5000 bases as default
				if ($variantInformation{'ENSPCANNO_GENE'} eq "NA") {
					my ($five, $six) = PatternModule::GeneAnno::getNearestGene($vf, $slice_adaptor, "ANY", $geneint);
					$variantInformation{'ENSPCANNO_GENE'} = $five;
					$variantInformation{'ENSPCANNO_TRANSCRIPT_CONSEQUENCE'} = $six;
				}
				
				# Adding length of gene
					if ($variantInformation{'ENSPCANNO_GENE'} ne "NA") {
						my $gene = $gene_adaptor->fetch_by_display_label($variantInformation{'ENSPCANNO_GENE'});
						$variantInformation{'ENSPCANNO_GENE_LENGTH'} = $gene->length();
						$variantInformation{'ENSPCANNO_TRANSCRIPT_START_SITE'} = $gene->canonical_transcript->seq_region_start();	
				}
				
				# Adding information to the annotation
				$annotationString .= 
					$variantInformation{'ENSPCANNO_GENE'}."\t".
					$variantInformation{'ENSPCANNO_TRANSCRIPT_CONSEQUENCE'}."\t".
					$variantInformation{'ENSPCANNO_GENE_LENGTH'}."\t".
					$variantInformation{'ENSPCANNO_TRANSCRIPT_START_SITE'}."\t";			
			}
			
			#############################################################
			##### Annotation of transcription factor binding motifs #####
			#############################################################
			
			#######################################################
			##### Annotating variants not of the type <.....> #####
			#######################################################
			
			if ($col[3] !~ "<.*>" && $col[4] !~ "<.*>" && $col[4] !~ /,/) {

				# Getting information on transcription factor binding modules that are hit by variants
				my $motifFeature_adaptor = $registry->get_adaptor('Human', 'funcgen', 'motiffeature');
				my $variantSlice = $slice_adaptor->fetch_by_region('chromosome',$col[0],$col[1],$col[1]-1+length($col[3]));
				my @motif_features = @{$motifFeature_adaptor->fetch_all_by_Slice($variantSlice)};
				my @keys = ();			
			
				if (scalar @motif_features > 0) {

				my %motifScoreHash = ();
				my %motifScoreSortingHash = ();
				
					for my $motif (@motif_features) {
						if ($vf->start() <= $motif->seq_region_end() && $motif->seq_region_start() <= $vf->end()) {
							my $slice = $slice_adaptor->fetch_by_region('chromosome', $vf->slice->seq_region_name(), $motif->seq_region_start(), $motif->seq_region_end());				
												
							# For SNPs
							if ($vf->allele_string() !~ /-/ && length($vf->allele_string()) == 3) {
								if ($motif->strand() == 1) {
									# Calculating binding scores
									my $refString = $slice->seq();
									my $refScore = $motif->binding_matrix->relative_affinity($refString);
									
									# Creating a sequence including the SNV
									my $SNVSlice = $slice_adaptor->fetch_by_region('chromosome', $vf->slice->seq_region_name(), $motif->seq_region_start() - 1000, $motif->seq_region_end() + 1000);
									my $SNVString = $SNVSlice->seq();
									substr($SNVString, $vf->start() - $motif->seq_region_start() + 1000, 1) = substr($vf->allele_string(), 2, 1);

									my @newBindingScores = ();	
													
									# Scoring all adjoining sequences and taking the highest scoring sequence as the new TFBS 					
									for (my $i = 0; $i < $motif->length + 2; $i++) {
										my $altString = substr($SNVString, 1000 + $vf->start() - $motif->seq_region_start() - $motif->length() + $i, $motif->length());
										push @newBindingScores, $motif->binding_matrix->relative_affinity($altString);
									}
									my $totalScore = max(@newBindingScores)-$refScore;
									
									# Motif start
									my $StartPos = $motif->seq_region_start();
									
									# Length of motif
									my $motiflength = $motif->binding_matrix->length();
								
									# Distance from variant to start of motif
									my $DistToStart = $vf->start() - $motif->seq_region_start();
								
									# Size of indel
									my $SizeOfIndel = "NA";
							
									# Name of motif
									my @motifName = $motif->display_label();

									# Conservation score calculation
									my $scores = $csscore_adaptor->fetch_all_by_MethodLinkSpeciesSet_Slice($mlss, $slice, $motiflength);
									my ($consScore) = calculateConservationScore($scores);
									
									# TFBS information content
									my $infoCont = ();
									for (my $i = 0; $i < $motif->length(); $i++) {
										if ($i != $motif->length() - 1) {
											$infoCont .= $motif->binding_matrix->{'ic'}->[$i] . ",";
										}
										if ($i == $motif->length() - 1) {
											$infoCont .= $motif->binding_matrix->{'ic'}->[$i];
										}
									}
								
									# Placing information into Hash. If the same motif has two scores for a variant, the most severe score is stored
									unless (exists $motifScoreHash{$motifName[0]} and $motifScoreHash{$motifName[0]}{Difference} <= $totalScore) {
										$motifScoreHash{$motifName[0]}{Start} = $StartPos;
										$motifScoreHash{$motifName[0]}{Length} = $motiflength;
										$motifScoreHash{$motifName[0]}{DistToStart} = $DistToStart;
										$motifScoreHash{$motifName[0]}{SizeOfIndel} = $SizeOfIndel;
										$motifScoreHash{$motifName[0]}{Original} = $refScore;
										$motifScoreHash{$motifName[0]}{New} = max(@newBindingScores);
										$motifScoreHash{$motifName[0]}{Difference} = $totalScore;
										$motifScoreHash{$motifName[0]}{InfoContent} = $infoCont;
										$motifScoreHash{$motifName[0]}{ConsScore} = $consScore;
									}
								}
						
								if ($motif->strand() == -1) {
									# Calculating binding scores
									my $refString = PatternModule::RegAnno::createComplementarySequence($slice->seq());
									my $refScore = $motif->binding_matrix->relative_affinity($refString);
									
									# Creating a sequence including the SNV
									my $SNVSlice = $slice_adaptor->fetch_by_region('chromosome', $vf->slice->seq_region_name(), $motif->seq_region_start() - 1000, $motif->seq_region_end() + 1000);
									my $SNVString = $SNVSlice->seq();
									substr($SNVString, $vf->start() - $motif->seq_region_start() + 1000, 1) = substr($vf->allele_string(), 2, 1);
									my @newBindingScores = ();	
													
									# Scoring all adjoining sequences and taking the highest scoring sequence as the new TFBS 					
									for (my $i = 0; $i < $motif->length + 2; $i++) {
										my $altString = PatternModule::RegAnno::createComplementarySequence(substr($SNVString, 1000 + $vf->start() - $motif->seq_region_start() - $motif->length() + $i, $motif->length()));
										push @newBindingScores, $motif->binding_matrix->relative_affinity($altString);
									}
									my $totalScore = max(@newBindingScores)-$refScore;	
									
									# Motif start
									my $StartPos = $motif->seq_region_start();
									
									# Length of motif
									my $motiflength = $motif->binding_matrix->length();
							
									# Distance from variant to start of motif
									my $DistToStart = $motif->seq_region_end() - $vf->start();
								
									# Size of indel
									my $SizeOfIndel = "NA";
								
									# Name of motif
									my @motifName = $motif->display_label();
									
									# Conservation score calculation
									my $scores = $csscore_adaptor->fetch_all_by_MethodLinkSpeciesSet_Slice($mlss, $slice, $motiflength);
									my ($consScore) = calculateConservationScore($scores);
									
									# TFBS information content
									my $infoCont = ();
									for (my $i = 0; $i < $motif->length(); $i++) {
										if ($i != $motif->length() - 1) {
											$infoCont .= $motif->binding_matrix->{'ic'}->[$i] . ",";
										}
										if ($i == $motif->length() - 1) {
											$infoCont .= $motif->binding_matrix->{'ic'}->[$i];
										}
									}
								
									unless (exists $motifScoreHash{$motifName[0]} and $motifScoreHash{$motifName[0]}{Difference} <= $totalScore) {
										$motifScoreHash{$motifName[0]}{Start} = $StartPos;
										$motifScoreHash{$motifName[0]}{Length} = $motiflength;
										$motifScoreHash{$motifName[0]}{DistToStart} = $DistToStart;
										$motifScoreHash{$motifName[0]}{SizeOfIndel} = $SizeOfIndel;
										$motifScoreHash{$motifName[0]}{Original} = $refScore;
										$motifScoreHash{$motifName[0]}{New} = max(@newBindingScores);
										$motifScoreHash{$motifName[0]}{Difference} = $totalScore;
										$motifScoreHash{$motifName[0]}{InfoContent} = $infoCont;
										$motifScoreHash{$motifName[0]}{ConsScore} = $consScore;
									}
								}	
							}
				
							# For Substitutions
							if ($vf->allele_string() !~ /-/ && length($vf->allele_string()) > 3) {
								if ($motif->strand() == 1) {
								
									# Calculating binding scores
									my $refString = $slice->seq();
									my $refScore = $motif->binding_matrix->relative_affinity($refString);
						
									# Creating a sequence including the substitution	
									# A slice is made with 1000 bases in each direction of the motif
									my $deletionSlice = $slice_adaptor->fetch_by_region('chromosome', $vf->slice->seq_region_name(), $motif->seq_region_start() - 1000, $motif->seq_region_end() + 1000);
									my $deletionString = $deletionSlice->seq();
									my @varrefalt = split("/", $vf->allele_string());
									substr($deletionString, $vf->start() - $motif->seq_region_start() + 1000, $vf->end() - $vf->start() + 1) = $varrefalt[1];
									my @newBindingScores = ();						
									for (my $i = 0; $i < $motif->length + 1; $i++) {
										my $altString = substr($deletionString, 1000 + $vf->start() - $motif->seq_region_start() - $motif->length() + $i, $motif->length());
										push @newBindingScores, $motif->binding_matrix->relative_affinity($altString);
									}
									my $totalScore = max(@newBindingScores)-$refScore;
									
									# Motif start
									my $StartPos = $motif->seq_region_start();
									
									# Length of motif
									my $motiflength = $motif->binding_matrix->length();
					
									# Distance from variant to start of motif
									my $DistToStart = -1;
									if ($vf->start() - $motif->seq_region_start() > -1) {
										$DistToStart = $vf->start() - $motif->seq_region_start();
									}
								
									# Size of substitution
									my @splitallele = split("/", $vf->allele_string());
									my $SizeOfsubstitution = length($splitallele[0]);
								
									# If the substitution is of the entire transcription factor binding motif
									if ($vf->start() - $motif->seq_region_start() <= 0 && $motif->seq_region_end() - $vf->end() <= 0) {
										$SizeOfsubstitution = $motif->binding_matrix->length();
									}
								
									# If the substitution start/end both lie within the transcription factor binding motif
									if ($vf->start() - $motif->seq_region_start() >= 0 && $motif->seq_region_end() - $vf->end() >= 0) {
										$SizeOfsubstitution = length($splitallele[0]);
									}
								
									# If the substitution start lies before the motif start, and the substitution end lies within the motif
									if ($vf->start() - $motif->seq_region_start() <= 0 && $motif->seq_region_end() - $vf->end() >= 0) {
										$SizeOfsubstitution = $vf->end() - $motif->seq_region_start() + 1;
									}
								
									# If the substitution start lies within motif, and the substitution end lies outside motif
									if ($vf->start() - $motif->seq_region_start() >= 0 && $motif->seq_region_end() - $vf->end() <= 0) {
										$SizeOfsubstitution = $motif->seq_region_end() - $vf->start() + 1;
									}								
								
									# Name of motif
									my @motifName = $motif->display_label();
									
									# Conservation score calculation
									my $scores = $csscore_adaptor->fetch_all_by_MethodLinkSpeciesSet_Slice($mlss, $slice, $motiflength);
									my ($consScore) = calculateConservationScore($scores);
									
									# TFBS information content
									my $infoCont = ();
									for (my $i = 0; $i < $motif->length(); $i++) {
										if ($i != $motif->length() - 1) {
											$infoCont .= $motif->binding_matrix->{'ic'}->[$i] . ",";
										}
										if ($i == $motif->length() - 1) {
											$infoCont .= $motif->binding_matrix->{'ic'}->[$i];
										}
									}
								
									unless (exists $motifScoreHash{$motifName[0]} and $motifScoreHash{$motifName[0]}{Difference} <= $totalScore) {
										$motifScoreHash{$motifName[0]}{Start} = $StartPos;
										$motifScoreHash{$motifName[0]}{Length} = $motiflength;
										$motifScoreHash{$motifName[0]}{DistToStart} = $DistToStart;
										$motifScoreHash{$motifName[0]}{SizeOfsubstitution} = $SizeOfsubstitution;
										$motifScoreHash{$motifName[0]}{Original} = $refScore;
										$motifScoreHash{$motifName[0]}{New} = max(@newBindingScores);
										$motifScoreHash{$motifName[0]}{Difference} = $totalScore;
										$motifScoreHash{$motifName[0]}{InfoContent} = $infoCont;
										$motifScoreHash{$motifName[0]}{ConsScore} = $consScore;
									}
								}
								if ($motif->strand == -1) {
								
									# Calculating binding scores
									my $refString = PatternModule::RegAnno::createComplementarySequence($slice->seq());
									my $refScore = $motif->binding_matrix->relative_affinity($refString);
						
									# Creating a sequence including the substitution
									# A slice is made with 1000 bases in each direction of the motif
									my $deletionSlice = $slice_adaptor->fetch_by_region('chromosome', $vf->slice->seq_region_name(), $motif->seq_region_start() - 1000, $motif->seq_region_end() + 1000);
									my $deletionString = $deletionSlice->seq();
									my @varrefalt = split("/", $vf->allele_string());
									substr($deletionString, $vf->start - $motif->seq_region_start() + 1000, $vf->end() - $vf->start() + 1) = $varrefalt[1];
									my @newBindingScores = ();						
									for (my $i = 0; $i < $motif->length() + 1; $i++) {
										my $altString = PatternModule::RegAnno::createComplementarySequence(substr($deletionString, 1000 + $vf->start() - $motif->seq_region_start() - $motif->length() + $i, $motif->length()));
										push @newBindingScores, $motif->binding_matrix->relative_affinity($altString);
									}
									my $totalScore = max(@newBindingScores)-$refScore;
									
									# Motif start
									my $StartPos = $motif->seq_region_start();
									
									# Length of motif
									my $motiflength = $motif->binding_matrix->length();

									# Distance from variant to start of motif
									my $DistToStart = -1;
									if ($motif->seq_region_end() - $vf->end() > -1) {
										$DistToStart = $motif->seq_region_end() - $vf->end();
									}
								
									# Size of substitution
									my @splitallele = split("/", $vf->allele_string());
									my $SizeOfsubstitution = length($splitallele[0]);
								
									# If the substitution is of the entire transcription factor binding motif
									if ($vf->start() - $motif->seq_region_start() <= 0 && $motif->seq_region_end() - $vf->end() <= 0) {
										$SizeOfsubstitution = $motif->binding_matrix->length();
									}
								
									# If the substitution start/end both lie within the transcription factor binding motif
									if ($vf->start() - $motif->seq_region_start() >= 0 && $motif->seq_region_end() - $vf->end() >= 0) {
										$SizeOfsubstitution = length($splitallele[0]);
									}
								
									# If the substitution start lies before the motif start, and the substitution end lies within the motif
									if ($vf->start() - $motif->seq_region_start() <= 0 && $motif->seq_region_end() - $vf->end() >= 0) {
										$SizeOfsubstitution = $vf->end() - $motif->seq_region_start() + 1;
									}
								
									# If the substitution start lies within motif, and the substitution end lies outside motif
									if ($vf->start() - $motif->seq_region_start() >= 0 && $motif->seq_region_end() - $vf->end() <= 0) {
										$SizeOfsubstitution = $motif->seq_region_end() - $vf->start() + 1;
									}		
								
									# Name of motif
									my @motifName = $motif->display_label();
									
									# Conservation score calculation
									my $scores = $csscore_adaptor->fetch_all_by_MethodLinkSpeciesSet_Slice($mlss, $slice, $motiflength);
									my ($consScore) = calculateConservationScore($scores);
									
									# TFBS information content
									my $infoCont = ();
									for (my $i = 0; $i < $motif->length(); $i++) {
										if ($i != $motif->length() - 1) {
											$infoCont .= $motif->binding_matrix->{'ic'}->[$i] . ",";
										}
										if ($i == $motif->length() - 1) {
											$infoCont .= $motif->binding_matrix->{'ic'}->[$i];
										}
									}
								
									unless (exists $motifScoreHash{$motifName[0]} and $motifScoreHash{$motifName[0]}{Difference} <= $totalScore) {
										$motifScoreHash{$motifName[0]}{Start} = $StartPos;
										$motifScoreHash{$motifName[0]}{Length} = $motiflength;
										$motifScoreHash{$motifName[0]}{DistToStart} = $DistToStart;
										$motifScoreHash{$motifName[0]}{SizeOfsubstitution} = $SizeOfsubstitution;
										$motifScoreHash{$motifName[0]}{Original} = $refScore;
										$motifScoreHash{$motifName[0]}{New} = max(@newBindingScores);
										$motifScoreHash{$motifName[0]}{Difference} = $totalScore;
										$motifScoreHash{$motifName[0]}{InfoContent} = $infoCont;
										$motifScoreHash{$motifName[0]}{ConsScore} = $consScore;
									}
								}
							}
				
							# For Insertions
							if (substr($vf->allele_string(), 0, 1) eq "-") {	
								if ($motif->strand() == 1) {
									# Calculating binding scores
									my $refString = $slice->seq();
									my $refScore = $motif->binding_matrix->relative_affinity($refString);
							
									# Creating a sequence including the insertion
									# A slice is made with 1000 bases in each direction of the motif
									my $insertionSlice = $slice_adaptor->fetch_by_region('chromosome', $vf->slice->seq_region_name(), $motif->seq_region_start() - 1000, $motif->seq_region_end() + 1000);
									my $insertionString = $insertionSlice->seq();
									my $insertion = substr($vf->allele_string(), 2);
									$insertionString = substr($insertionString, 0, 1000 + $vf->start() - $motif->seq_region_start()) . $insertion . substr($insertionString, 1000 + $vf->start() - $motif->seq_region_start());							
									my @newBindingScores = ();
									for (my $i = 0; $i < $motif->length + length($insertion) + 1; $i++) {
										my $altString = substr($insertionString, 1000 + $vf->start() - $motif->seq_region_start() - $motif->length() + $i, $motif->seq_region_end() - $motif->seq_region_start() + 1);
										push @newBindingScores, $motif->binding_matrix->relative_affinity($altString);
									}
									my $totalScore = max(@newBindingScores)-$refScore;
									
									# Motif start
									my $StartPos = $motif->seq_region_start();
								
									# Length of motif
									my $motiflength = $motif->binding_matrix->length();
								
									# Distance from variant to start of motif
									my $DistToStart = $vf->start() - $motif->seq_region_start();
								
									# Size of indel
									my @splitallele = split("/", $vf->allele_string());
									my $SizeOfIndel = length($splitallele[1]);
								
									# Name of motif
									my @motifName = $motif->display_label();
									
									# Conservation score calculation
									my $scores = $csscore_adaptor->fetch_all_by_MethodLinkSpeciesSet_Slice($mlss, $slice, $motiflength);
									my ($consScore) = calculateConservationScore($scores);
									
									# TFBS information content
									my $infoCont = ();
									for (my $i = 0; $i < $motif->length(); $i++) {
										if ($i != $motif->length() - 1) {
											$infoCont .= $motif->binding_matrix->{'ic'}->[$i] . ",";
										}
										if ($i == $motif->length() - 1) {
											$infoCont .= $motif->binding_matrix->{'ic'}->[$i];
										}
									}
								
									unless (exists $motifScoreHash{$motifName[0]} and $motifScoreHash{$motifName[0]}{Difference} <= $totalScore) {
										$motifScoreHash{$motifName[0]}{Start} = $StartPos;
										$motifScoreHash{$motifName[0]}{Length} = $motiflength;
										$motifScoreHash{$motifName[0]}{DistToStart} = $DistToStart;
										$motifScoreHash{$motifName[0]}{SizeOfIndel} = $SizeOfIndel;
										$motifScoreHash{$motifName[0]}{Original} = $refScore;
										$motifScoreHash{$motifName[0]}{New} = max(@newBindingScores);
										$motifScoreHash{$motifName[0]}{Difference} = $totalScore;
										$motifScoreHash{$motifName[0]}{InfoContent} = $infoCont;
										$motifScoreHash{$motifName[0]}{ConsScore} = $consScore;
									}
								}
						
								if ($motif->strand == -1) {
									# Calculating binding scores
									my $refString = PatternModule::RegAnno::createComplementarySequence($slice->seq);
									my $refScore = $motif->binding_matrix->relative_affinity($refString);
							
									# Creating a sequence including the insertion
									# A slice is made with 1000 bases in each direction of the motif
									my $insertionSlice = $slice_adaptor->fetch_by_region('chromosome', $vf->slice->seq_region_name(), $motif->seq_region_start() - 1000, $motif->seq_region_end() + 1000);
									my $insertionString = $insertionSlice->seq();
									my $insertion = substr($vf->allele_string(), 2);
									$insertionString = substr($insertionString, 0, 1000 + $vf->start() - $motif->seq_region_start()) . $insertion . substr($insertionString, 1000 + $vf->start() - $motif->seq_region_start());
									my @newBindingScores = ();
									for (my $i = 0; $i < $motif->length + length($insertion) + 1; $i++) {
										my $altString = PatternModule::RegAnno::createComplementarySequence(substr($insertionString, 1000 + $vf->start() - $motif->seq_region_start() - $motif->length() + $i, $motif->seq_region_end() - $motif->seq_region_start() + 1));
										push @newBindingScores, $motif->binding_matrix->relative_affinity($altString);
									}
									my $totalScore = max(@newBindingScores)-$refScore;
									
									# Motif start
									my $StartPos = $motif->seq_region_start();
									
									# Length of motif
									my $motiflength = $motif->binding_matrix->length();
								
									# Distance from variant to start of motif
									my $DistToStart = $motif->seq_region_end() - $vf->start();
								
									# Size of indel
									my @splitallele = split("/", $vf->allele_string());
									my $SizeOfIndel = length($splitallele[1]);
								
									# Name of motif
									my @motifName = $motif->display_label();
									
									# Conservation score calculation
									my $scores = $csscore_adaptor->fetch_all_by_MethodLinkSpeciesSet_Slice($mlss, $slice, $motiflength);
									my ($consScore) = calculateConservationScore($scores);
									
									# TFBS information content
									my $infoCont = ();
									for (my $i = 0; $i < $motif->length(); $i++) {
										if ($i != $motif->length() - 1) {
											$infoCont .= $motif->binding_matrix->{'ic'}->[$i] . ",";
										}
										if ($i == $motif->length() - 1) {
											$infoCont .= $motif->binding_matrix->{'ic'}->[$i];
										}
									}
								
									unless (exists $motifScoreHash{$motifName[0]} and $motifScoreHash{$motifName[0]}{Difference} <= $totalScore) {
										$motifScoreHash{$motifName[0]}{Start} = $StartPos;
										$motifScoreHash{$motifName[0]}{Length} = $motiflength;
										$motifScoreHash{$motifName[0]}{DistToStart} = $DistToStart;
										$motifScoreHash{$motifName[0]}{SizeOfIndel} = $SizeOfIndel;
										$motifScoreHash{$motifName[0]}{Original} = $refScore;
										$motifScoreHash{$motifName[0]}{New} = max(@newBindingScores);
										$motifScoreHash{$motifName[0]}{Difference} = $totalScore;
										$motifScoreHash{$motifName[0]}{InfoContent} = $infoCont;
										$motifScoreHash{$motifName[0]}{ConsScore} = $consScore;
									}
								}
							}
				
							# For Deletions
							if (substr($vf->allele_string(), -1) eq "-") {
								if ($motif->strand() == 1) {
								
									# Calculating binding scores
									my $refString = $slice->seq();
									my $refScore = $motif->binding_matrix->relative_affinity($refString);
						
									# Creating a sequence including the Deletion	
									# A slice is made with 1000 bases in each direction of the motif
									my $deletionSlice = $slice_adaptor->fetch_by_region('chromosome', $vf->slice->seq_region_name(), $motif->seq_region_start() - 1000, $motif->seq_region_end() + 1000);
									my $deletionString = $deletionSlice->seq();
									substr($deletionString, $vf->start() - $motif->seq_region_start() + 1000, $vf->end() - $vf->start() + 1) = "";
									my @newBindingScores = ();						
									for (my $i = 0; $i < $motif->length + 1; $i++) {
										my $altString = substr($deletionString, 1000 + $vf->start() - $motif->seq_region_start() - $motif->length() + $i, $motif->length());
										push @newBindingScores, $motif->binding_matrix->relative_affinity($altString);
									}
									my $totalScore = max(@newBindingScores)-$refScore;
									
									# Motif start
									my $StartPos = $motif->seq_region_start();
								
									# Length of motif
									my $motiflength = $motif->binding_matrix->length();
					
									# Distance from variant to start of motif
									my $DistToStart = -1;
									if ($vf->start() - $motif->seq_region_start() > -1) {
										$DistToStart = $vf->start() - $motif->seq_region_start();
									}
								
									# Size of indel
									my @splitallele = split("/", $vf->allele_string());
									my $SizeOfIndel = length($splitallele[0]);
								
									# If the deletion is of the entire transcription factor binding motif
									if ($vf->start() - $motif->seq_region_start() <= 0 && $motif->seq_region_end() - $vf->end() <= 0) {
										$SizeOfIndel = $motif->binding_matrix->length();
									}
								
									# If the deletion start/end both lie within the transcription factor binding motif
									if ($vf->start() - $motif->seq_region_start() >= 0 && $motif->seq_region_end() - $vf->end() >= 0) {
										$SizeOfIndel = length($splitallele[0]);
									}
								
									# If the deletion start lies before the motif start, and the deletion end lies within the motif
									if ($vf->start() - $motif->seq_region_start() <= 0 && $motif->seq_region_end() - $vf->end() >= 0) {
										$SizeOfIndel = $vf->end() - $motif->seq_region_start() + 1;
									}
								
									# If the deletion start lies within motif, and the deletion end lies outside motif
									if ($vf->start() - $motif->seq_region_start() >= 0 && $motif->seq_region_end() - $vf->end() <= 0) {
										$SizeOfIndel = $motif->seq_region_end() - $vf->start() + 1;
									}								
								
									# Name of motif
									my @motifName = $motif->display_label();
									
									# Conservation score calculation
									my $scores = $csscore_adaptor->fetch_all_by_MethodLinkSpeciesSet_Slice($mlss, $slice, $motiflength);
									my ($consScore) = calculateConservationScore($scores);
									
									# TFBS information content
									my $infoCont = ();
									for (my $i = 0; $i < $motif->length(); $i++) {
										if ($i != $motif->length() - 1) {
											$infoCont .= $motif->binding_matrix->{'ic'}->[$i] . ",";
										}
										if ($i == $motif->length() - 1) {
											$infoCont .= $motif->binding_matrix->{'ic'}->[$i];
										}
									}
								
									unless (exists $motifScoreHash{$motifName[0]} and $motifScoreHash{$motifName[0]}{Difference} <= $totalScore) {
										$motifScoreHash{$motifName[0]}{Start} = $StartPos;
										$motifScoreHash{$motifName[0]}{Length} = $motiflength;
										$motifScoreHash{$motifName[0]}{DistToStart} = $DistToStart;
										$motifScoreHash{$motifName[0]}{SizeOfIndel} = $SizeOfIndel;
										$motifScoreHash{$motifName[0]}{Original} = $refScore;
										$motifScoreHash{$motifName[0]}{New} = max(@newBindingScores);
										$motifScoreHash{$motifName[0]}{Difference} = $totalScore;
										$motifScoreHash{$motifName[0]}{InfoContent} = $infoCont;
										$motifScoreHash{$motifName[0]}{ConsScore} = $consScore;
									}
								}
								if ($motif->strand == -1) {
								
									# Calculating binding scores
									my $refString = PatternModule::RegAnno::createComplementarySequence($slice->seq());
									my $refScore = $motif->binding_matrix->relative_affinity($refString);
						
									# Creating a sequence including the Deletion
									# A slice is made with 1000 bases in each direction of the motif
									my $deletionSlice = $slice_adaptor->fetch_by_region('chromosome', $vf->slice->seq_region_name(), $motif->seq_region_start() - 1000, $motif->seq_region_end() + 1000);
									my $deletionString = $deletionSlice->seq();
									substr($deletionString, $vf->start - $motif->seq_region_start() + 1000, $vf->end() - $vf->start() + 1) = "";
									my @newBindingScores = ();						
									for (my $i = 0; $i < $motif->length() + 1; $i++) {
										my $altString = PatternModule::RegAnno::createComplementarySequence(substr($deletionString, 1000 + $vf->start() - $motif->seq_region_start() - $motif->length() + $i, $motif->length()));
										push @newBindingScores, $motif->binding_matrix->relative_affinity($altString);
									}
									my $totalScore = max(@newBindingScores)-$refScore;
									
									# Motif start
									my $StartPos = $motif->seq_region_start();
									
									# Length of motif
									my $motiflength = $motif->binding_matrix->length();

									# Distance from variant to start of motif
									my $DistToStart = -1;
									if ($motif->seq_region_end() - $vf->end() > -1) {
										$DistToStart = $motif->seq_region_end() - $vf->end();
									}
								
									# Size of indel
									my @splitallele = split("/", $vf->allele_string());
									my $SizeOfIndel = length($splitallele[0]);
								
									# If the deletion is of the entire transcription factor binding motif
									if ($vf->start() - $motif->seq_region_start() <= 0 && $motif->seq_region_end() - $vf->end() <= 0) {
										$SizeOfIndel = $motif->binding_matrix->length();
									}
								
									# If the deletion start/end both lie within the transcription factor binding motif
									if ($vf->start() - $motif->seq_region_start() >= 0 && $motif->seq_region_end() - $vf->end() >= 0) {
										$SizeOfIndel = length($splitallele[0]);
									}
								
									# If the deletion start lies before the motif start, and the deletion end lies within the motif
									if ($vf->start() - $motif->seq_region_start() <= 0 && $motif->seq_region_end() - $vf->end() >= 0) {
										$SizeOfIndel = $vf->end() - $motif->seq_region_start() + 1;
									}
								
									# If the deletion start lies within motif, and the deletion end lies outside motif
									if ($vf->start() - $motif->seq_region_start() >= 0 && $motif->seq_region_end() - $vf->end() <= 0) {
										$SizeOfIndel = $motif->seq_region_end() - $vf->start() + 1;
									}		
								
									# Name of motif
									my @motifName = $motif->display_label();
									
									# Conservation score calculation
									my $scores = $csscore_adaptor->fetch_all_by_MethodLinkSpeciesSet_Slice($mlss, $slice, $motiflength);
									my ($consScore) = calculateConservationScore($scores);
									
									# TFBS information content
									my $infoCont = ();
									for (my $i = 0; $i < $motif->length(); $i++) {
										if ($i != $motif->length() - 1) {
											$infoCont .= $motif->binding_matrix->{'ic'}->[$i] . ",";
										}
										if ($i == $motif->length() - 1) {
											$infoCont .= $motif->binding_matrix->{'ic'}->[$i];
										}
									}
								
									unless (exists $motifScoreHash{$motifName[0]} and $motifScoreHash{$motifName[0]}{Difference} <= $totalScore) {
										$motifScoreHash{$motifName[0]}{Start} = $StartPos;
										$motifScoreHash{$motifName[0]}{Length} = $motiflength;
										$motifScoreHash{$motifName[0]}{DistToStart} = $DistToStart;
										$motifScoreHash{$motifName[0]}{SizeOfIndel} = $SizeOfIndel;
										$motifScoreHash{$motifName[0]}{Original} = $refScore;
										$motifScoreHash{$motifName[0]}{New} = max(@newBindingScores);
										$motifScoreHash{$motifName[0]}{Difference} = $totalScore;
										$motifScoreHash{$motifName[0]}{InfoContent} = $infoCont;
										$motifScoreHash{$motifName[0]}{ConsScore} = $consScore;
									}
								}
							}
						}
					}
				
					# TFBS names are extracted
					@keys = sort keys %motifScoreHash;

					# A new simple hash with only TFBS name, and Difference between old and new binding scores are present
					for my $key (@keys) {
						$motifScoreSortingHash{$key} = $motifScoreHash{$key}{Difference};
					}
				
					# Keys are sorted such that the lowest Difference-score is first
					my @DifferenceSortedKeys = sort {$motifScoreSortingHash{$a} <=> $motifScoreSortingHash{$b}} keys %motifScoreSortingHash;
				
					# The TFBS information is added to the annotation string
					for my $Rkey (@DifferenceSortedKeys) {
						# Splitting TFBS name and matrix into two
						my @TFBSnames = split(":", $Rkey);
						
						$annotationString .= $TFBSnames[0] . "\t" . $TFBSnames[-1] . "\t" . $motifScoreHash{$Rkey}{Start} . "\t" . $motifScoreHash{$Rkey}{Length} . "\t" . $motifScoreHash{$Rkey}{DistToStart} . "\t" . $motifScoreHash{$Rkey}{SizeOfIndel} . "\t" . $motifScoreHash{$Rkey}{Original} . "\t" . $motifScoreHash{$Rkey}{New} . "\t" . $motifScoreHash{$Rkey}{InfoContent} . "\t" . $motifScoreHash{$Rkey}{ConsScore}  . "\t";
					}
				
					# The empty spaces in the annotation string are filled out with NA
					if (scalar @keys != 0) {
						for (my $i = 0; $i < 10 - scalar @keys; $i++) {
							$annotationString .= "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t";
						}
					}
				}	
			
				# If the variant overlaps with no TFBS, NA is filled in into the annotation string
				if (scalar @keys == 0) {
					$annotationString .= "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" ."NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t". "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t". "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t";
				}
			}
						
			####################################################################
			##### Annotating multiallelic variants not of the type <.....> #####
			####################################################################

			if ($col[3] !~ "<.*>" && $col[4] !~ "<.*>" && $col[4] =~ /,/) {				
				# Getting information on transcription factor binding modules that are hit by variants
				my @variationFeatureList = PatternModule::GeneAnno::splitMultiAllelicVariationFeature($vf, $slice, $variant_feature_adaptor);							
				my $motifFeature_adaptor = $registry->get_adaptor('Human', 'funcgen', 'motiffeature');
				my $variantSlice = $slice_adaptor->fetch_by_region('chromosome',$col[0],$col[1],$col[1]-1+length($col[3]));
				my @motif_features = @{$motifFeature_adaptor->fetch_all_by_Slice($variantSlice)};
				my @keys = ();
			
				if (scalar @motif_features > 0) {
				my %motifScoreHash = ();
				my %motifScoreSortingHash = ();
				
					for my $varfeat (@variationFeatureList) {
						for my $motif (@motif_features) {
							if ($varfeat->start() <= $motif->seq_region_end() && $motif->seq_region_start() <= $varfeat->end()) {
								my $slice = $slice_adaptor->fetch_by_region('chromosome', $varfeat->slice->seq_region_name(), $motif->seq_region_start(), $motif->seq_region_end());				
												
								# For SNPs
								if ($vf->allele_string() !~ /-/ && length($vf->allele_string()) == 3) {
									if ($motif->strand() == 1) {
								
										# Calculating binding scores
										my $refString = $slice->seq();
										my $refScore = $motif->binding_matrix->relative_affinity($refString);
										my $snpIndex = $varfeat->start() - $slice->start();
										my $altString = $refString;
										substr($altString, $snpIndex, 1) = substr($varfeat->allele_string(), 2, 1);
										my $altScore = $motif->binding_matrix->relative_affinity($altString);
										my $totalScore = $altScore-$refScore;
										
										# Motif start
										my $StartPos = $motif->seq_region_start();
										
										# Length of motif
										my $motiflength = $motif->binding_matrix->length();
								
										# Distance from variant to start of motif
										my $DistToStart = $varfeat->start() - $motif->seq_region_start();
								
										# Size of indel
										my $SizeOfIndel = "NA";
							
										# Name of motif
										my @motifName = $motif->display_label();
										
										# Conservation score calculation
										my $scores = $csscore_adaptor->fetch_all_by_MethodLinkSpeciesSet_Slice($mlss, $slice, $motiflength);
										my ($consScore) = calculateConservationScore($scores);
										
										# TFBS information content
										my $infoCont = ();
										for (my $i = 0; $i < $motif->length(); $i++) {
											if ($i != $motif->length() - 1) {
												$infoCont .= $motif->binding_matrix->{'ic'}->[$i] . ",";
											}
											if ($i == $motif->length() - 1) {
												$infoCont .= $motif->binding_matrix->{'ic'}->[$i];
											}
										}
								
										# Placing information into Hash. If the same motif has two scores for a variant, the most severe score is stored
										unless (exists $motifScoreHash{$motifName[0]} and $motifScoreHash{$motifName[0]}{Difference} <= $totalScore) {
											$motifScoreHash{$motifName[0]}{Start} = $StartPos;
											$motifScoreHash{$motifName[0]}{Length} = $motiflength;
											$motifScoreHash{$motifName[0]}{DistToStart} = $DistToStart;
											$motifScoreHash{$motifName[0]}{SizeOfIndel} = $SizeOfIndel;
											$motifScoreHash{$motifName[0]}{Original} = $refScore;
											$motifScoreHash{$motifName[0]}{New} = $altScore;
											$motifScoreHash{$motifName[0]}{Difference} = $totalScore;
											$motifScoreHash{$motifName[0]}{InfoContent} = $infoCont;
											$motifScoreHash{$motifName[0]}{ConsScore} = $consScore;
										}
									}
						
									if ($motif->strand() == -1) {
								
										# Calculating binding scores
										my $refString = PatternModule::RegAnno::createComplementarySequence($slice->seq());
										my $refScore = $motif->binding_matrix->relative_affinity($refString);
										my $snpIndex = $varfeat->start() - $slice->start();
										my $altString = $slice->seq();
										substr($altString, $snpIndex, 1) = substr($varfeat->allele_string(), 2, 1);
										$altString = PatternModule::RegAnno::createComplementarySequence($altString);
										my $altScore = $motif->binding_matrix->relative_affinity($altString);
										my $totalScore = $altScore-$refScore;
										
										# Motif start
										my $StartPos = $motif->seq_region_start();
										
										# Length of motif
										my $motiflength = $motif->binding_matrix->length();
							
										# Distance from variant to start of motif
										my $DistToStart = $motif->seq_region_end() - $varfeat->start();
								
										# Size of indel
										my $SizeOfIndel = "NA";
								
										# Name of motif
										my @motifName = $motif->display_label();
										
										# Conservation score calculation
										my $scores = $csscore_adaptor->fetch_all_by_MethodLinkSpeciesSet_Slice($mlss, $slice, $motiflength);
										my ($consScore) = calculateConservationScore($scores);
										
										# TFBS information content
										my $infoCont = ();
										for (my $i = 0; $i < $motif->length(); $i++) {
											if ($i != $motif->length() - 1) {
												$infoCont .= $motif->binding_matrix->{'ic'}->[$i] . ",";
											}
											if ($i == $motif->length() - 1) {
												$infoCont .= $motif->binding_matrix->{'ic'}->[$i];
											}
										}
								
										unless (exists $motifScoreHash{$motifName[0]} and $motifScoreHash{$motifName[0]}{Difference} <= $totalScore) {
											$motifScoreHash{$motifName[0]}{Start} = $StartPos;
											$motifScoreHash{$motifName[0]}{Length} = $motiflength;
											$motifScoreHash{$motifName[0]}{DistToStart} = $DistToStart;
											$motifScoreHash{$motifName[0]}{SizeOfIndel} = $SizeOfIndel;
											$motifScoreHash{$motifName[0]}{Original} = $refScore;
											$motifScoreHash{$motifName[0]}{New} = $altScore;
											$motifScoreHash{$motifName[0]}{Difference} = $totalScore;
											$motifScoreHash{$motifName[0]}{InfoContent} = $infoCont;
											$motifScoreHash{$motifName[0]}{ConsScore} = $consScore;
										}
									}	
								}
				
							# For Substitutions
							if ($varfeat->allele_string() !~ /-/ && length($varfeat->allele_string()) > 3) {
								if ($motif->strand() == 1) {
								
									# Calculating binding scores
									my $refString = $slice->seq();
									my $refScore = $motif->binding_matrix->relative_affinity($refString);
						
									# Creating a sequence including the substitution	
									# A slice is made with 1000 bases in each direction of the motif
									my $deletionSlice = $slice_adaptor->fetch_by_region('chromosome', $varfeat->slice->seq_region_name(), $motif->seq_region_start() - 1000, $motif->seq_region_end() + 1000);
									my $deletionString = $deletionSlice->seq();
									my @varrefalt = split("/", $varfeat->allele_string());
									substr($deletionString, $varfeat->start() - $motif->seq_region_start() + 1000, $varfeat->end() - $varfeat->start() + 1) = $varrefalt[1];
									my @newBindingScores = ();						
									for (my $i = 0; $i < $motif->length + 1; $i++) {
										my $altString = substr($deletionString, 1000 + $varfeat->start() - $motif->seq_region_start() - $motif->length() + $i, $motif->length());
										push @newBindingScores, $motif->binding_matrix->relative_affinity($altString);
									}
									my $totalScore = max(@newBindingScores)-$refScore;
									
									# Motif start
									my $StartPos = $motif->seq_region_start();
									
									# Length of motif
									my $motiflength = $motif->binding_matrix->length();
					
									# Distance from variant to start of motif
									my $DistToStart = -1;
									if ($varfeat->start() - $motif->seq_region_start() > -1) {
										$DistToStart = $varfeat->start() - $motif->seq_region_start();
									}
								
									# Size of substitution
									my @splitallele = split("/", $varfeat->allele_string());
									my $SizeOfsubstitution = length($splitallele[0]);
								
									# If the substitution is of the entire transcription factor binding motif
									if ($varfeat->start() - $motif->seq_region_start() <= 0 && $motif->seq_region_end() - $varfeat->end() <= 0) {
										$SizeOfsubstitution = $motif->binding_matrix->length();
									}
								
									# If the substitution start/end both lie within the transcription factor binding motif
									if ($varfeat->start() - $motif->seq_region_start() >= 0 && $motif->seq_region_end() - $varfeat->end() >= 0) {
										$SizeOfsubstitution = length($splitallele[0]);
									}
								
									# If the substitution start lies before the motif start, and the substitution end lies within the motif
									if ($varfeat->start() - $motif->seq_region_start() <= 0 && $motif->seq_region_end() - $varfeat->end() >= 0) {
										$SizeOfsubstitution = $varfeat->end() - $motif->seq_region_start() + 1;
									}
								
									# If the substitution start lies within motif, and the substitution end lies outside motif
									if ($varfeat->start() - $motif->seq_region_start() >= 0 && $motif->seq_region_end() - $varfeat->end() <= 0) {
										$SizeOfsubstitution = $motif->seq_region_end() - $varfeat->start() + 1;
									}								
								
									# Name of motif
									my @motifName = $motif->display_label();
									
									# Conservation score calculation
									my $scores = $csscore_adaptor->fetch_all_by_MethodLinkSpeciesSet_Slice($mlss, $slice, $motiflength);
									my ($consScore) = calculateConservationScore($scores);
								
									# TFBS information content
									my $infoCont = ();
									for (my $i = 0; $i < $motif->length(); $i++) {
										if ($i != $motif->length() - 1) {
											$infoCont .= $motif->binding_matrix->{'ic'}->[$i] . ",";
										}
										if ($i == $motif->length() - 1) {
											$infoCont .= $motif->binding_matrix->{'ic'}->[$i];
										}
									}
								
									unless (exists $motifScoreHash{$motifName[0]} and $motifScoreHash{$motifName[0]}{Difference} <= $totalScore) {
										$motifScoreHash{$motifName[0]}{Start} = $StartPos;
										$motifScoreHash{$motifName[0]}{Length} = $motiflength;
										$motifScoreHash{$motifName[0]}{DistToStart} = $DistToStart;
										$motifScoreHash{$motifName[0]}{SizeOfsubstitution} = $SizeOfsubstitution;
										$motifScoreHash{$motifName[0]}{Original} = $refScore;
										$motifScoreHash{$motifName[0]}{New} = max(@newBindingScores);
										$motifScoreHash{$motifName[0]}{Difference} = $totalScore;
										$motifScoreHash{$motifName[0]}{InfoContent} = $infoCont;
										$motifScoreHash{$motifName[0]}{ConsScore} = $consScore;
									}
								}
								if ($motif->strand == -1) {
								
									# Calculating binding scores
									my $refString = PatternModule::RegAnno::createComplementarySequence($slice->seq());
									my $refScore = $motif->binding_matrix->relative_affinity($refString);
						
									# Creating a sequence including the substitution
									# A slice is made with 1000 bases in each direction of the motif
									my $deletionSlice = $slice_adaptor->fetch_by_region('chromosome', $varfeat->slice->seq_region_name(), $motif->seq_region_start() - 1000, $motif->seq_region_end() + 1000);
									my $deletionString = $deletionSlice->seq();
									my @varrefalt = split("/", $varfeat->allele_string());
									substr($deletionString, $varfeat->start - $motif->seq_region_start() + 1000, $varfeat->end() - $varfeat->start() + 1) = $varrefalt[1];
									my @newBindingScores = ();						
									for (my $i = 0; $i < $motif->length() + 1; $i++) {
										my $altString = PatternModule::RegAnno::createComplementarySequence(substr($deletionString, 1000 + $varfeat->start() - $motif->seq_region_start() - $motif->length() + $i, $motif->length()));
										push @newBindingScores, $motif->binding_matrix->relative_affinity($altString);
									}
									my $totalScore = max(@newBindingScores)-$refScore;
									
									# Motif start
									my $StartPos = $motif->seq_region_start();
									
									# Length of motif
									my $motiflength = $motif->binding_matrix->length();

									# Distance from variant to start of motif
									my $DistToStart = -1;
									if ($motif->seq_region_end() - $varfeat->end() > -1) {
										$DistToStart = $motif->seq_region_end() - $varfeat->end();
									}
								
									# Size of substitution
									my @splitallele = split("/", $varfeat->allele_string());
									my $SizeOfsubstitution = length($splitallele[0]);
								
									# If the substitution is of the entire transcription factor binding motif
									if ($varfeat->start() - $motif->seq_region_start() <= 0 && $motif->seq_region_end() - $varfeat->end() <= 0) {
										$SizeOfsubstitution = $motif->binding_matrix->length();
									}
								
									# If the substitution start/end both lie within the transcription factor binding motif
									if ($varfeat->start() - $motif->seq_region_start() >= 0 && $motif->seq_region_end() - $varfeat->end() >= 0) {
										$SizeOfsubstitution = length($splitallele[0]);
									}
								
									# If the substitution start lies before the motif start, and the substitution end lies within the motif
									if ($varfeat->start() - $motif->seq_region_start() <= 0 && $motif->seq_region_end() - $varfeat->end() >= 0) {
										$SizeOfsubstitution = $varfeat->end() - $motif->seq_region_start() + 1;
									}
								
									# If the substitution start lies within motif, and the substitution end lies outside motif
									if ($varfeat->start() - $motif->seq_region_start() >= 0 && $motif->seq_region_end() - $varfeat->end() <= 0) {
										$SizeOfsubstitution = $motif->seq_region_end() - $varfeat->start() + 1;
									}		
								
									# Name of motif
									my @motifName = $motif->display_label();
									
									# Conservation score calculation
									my $scores = $csscore_adaptor->fetch_all_by_MethodLinkSpeciesSet_Slice($mlss, $slice, $motiflength);
									my ($consScore) = calculateConservationScore($scores);
									
									# TFBS information content
									my $infoCont = ();
									for (my $i = 0; $i < $motif->length(); $i++) {
										if ($i != $motif->length() - 1) {
											$infoCont .= $motif->binding_matrix->{'ic'}->[$i] . ",";
										}
										if ($i == $motif->length() - 1) {
											$infoCont .= $motif->binding_matrix->{'ic'}->[$i];
										}
									}
								
									unless (exists $motifScoreHash{$motifName[0]} and $motifScoreHash{$motifName[0]}{Difference} <= $totalScore) {
										$motifScoreHash{$motifName[0]}{Start} = $StartPos;
										$motifScoreHash{$motifName[0]}{Length} = $motiflength;
										$motifScoreHash{$motifName[0]}{DistToStart} = $DistToStart;
										$motifScoreHash{$motifName[0]}{SizeOfsubstitution} = $SizeOfsubstitution;
										$motifScoreHash{$motifName[0]}{Original} = $refScore;
										$motifScoreHash{$motifName[0]}{New} = max(@newBindingScores);
										$motifScoreHash{$motifName[0]}{Difference} = $totalScore;
										$motifScoreHash{$motifName[0]}{InfoContent} = $infoCont;
										$motifScoreHash{$motifName[0]}{ConsScore} = $consScore;
									}
								}
							}				
				
								# For Insertions
								if (substr($varfeat->allele_string(), 0, 1) eq "-") {											
									if ($motif->strand() == 1) {
								
										# Calculating binding scores
										my $refString = $slice->seq();
										my $refScore = $motif->binding_matrix->relative_affinity($refString);
							
										# Creating a sequence including the insertion
										# A slice is made with 1000 bases in each direction of the motif
										my $insertionSlice = $slice_adaptor->fetch_by_region('chromosome', $varfeat->slice->seq_region_name(), $motif->seq_region_start() - 1000, $motif->seq_region_end() + 1000);
										my $insertionString = $insertionSlice->seq();
										my $insertion = substr($varfeat->allele_string(), 2);
										$insertionString = substr($insertionString, 0, 1000 + $varfeat->start() - $motif->seq_region_start()) . $insertion . substr($insertionString, 1000 + $varfeat->start() - $motif->seq_region_start());							
										my @newBindingScores = ();
										for (my $i = 0; $i < $motif->length + length($insertion) + 1; $i++) {
											my $altString = substr($insertionString, 1000 + $varfeat->start() - $motif->seq_region_start() - $motif->length() + $i, $motif->seq_region_end() - $motif->seq_region_start() + 1);
											push @newBindingScores, $motif->binding_matrix->relative_affinity($altString);
										}
										my $totalScore = max(@newBindingScores)-$refScore;
										
										# Motif start
										my $StartPos = $motif->seq_region_start();
										
										# Length of motif
										my $motiflength = $motif->binding_matrix->length();
								
										# Distance from variant to start of motif
										my $DistToStart = $varfeat->start() - $motif->seq_region_start();
								
										# Size of indel
										my @splitallele = split("/", $varfeat->allele_string());
										my $SizeOfIndel = length($splitallele[1]);
								
										# Name of motif
										my @motifName = $motif->display_label();
										
										# Conservation score calculation
										my $scores = $csscore_adaptor->fetch_all_by_MethodLinkSpeciesSet_Slice($mlss, $slice, $motiflength);
										my ($consScore) = calculateConservationScore($scores);
										
										# TFBS information content
										my $infoCont = ();
										for (my $i = 0; $i < $motif->length(); $i++) {
											if ($i != $motif->length() - 1) {
												$infoCont .= $motif->binding_matrix->{'ic'}->[$i] . ",";
											}
											if ($i == $motif->length() - 1) {
												$infoCont .= $motif->binding_matrix->{'ic'}->[$i];
											}
										}
								
										unless (exists $motifScoreHash{$motifName[0]} and $motifScoreHash{$motifName[0]}{Difference} <= $totalScore) {
											$motifScoreHash{$motifName[0]}{Start} = $StartPos;
											$motifScoreHash{$motifName[0]}{Length} = $motiflength;
											$motifScoreHash{$motifName[0]}{DistToStart} = $DistToStart;
											$motifScoreHash{$motifName[0]}{SizeOfIndel} = $SizeOfIndel;
											$motifScoreHash{$motifName[0]}{Original} = $refScore;
											$motifScoreHash{$motifName[0]}{New} = max(@newBindingScores);
											$motifScoreHash{$motifName[0]}{Difference} = $totalScore;
											$motifScoreHash{$motifName[0]}{InfoContent} = $infoCont;
											$motifScoreHash{$motifName[0]}{ConsScore} = $consScore;
										}
									}
						
									if ($motif->strand == -1) {
								
										# Calculating binding scores
										my $refString = PatternModule::RegAnno::createComplementarySequence($slice->seq);
										my $refScore = $motif->binding_matrix->relative_affinity($refString);
							
										# Creating a sequence including the insertion
										# A slice is made with 1000 bases in each direction of the motif
										my $insertionSlice = $slice_adaptor->fetch_by_region('chromosome', $varfeat->slice->seq_region_name(), $motif->seq_region_start() - 1000, $motif->seq_region_end() + 1000);
										my $insertionString = $insertionSlice->seq();
										my $insertion = substr($varfeat->allele_string(), 2);
										$insertionString = substr($insertionString, 0, 1000 + $varfeat->start() - $motif->seq_region_start()) . $insertion . substr($insertionString, 1000 + $varfeat->start() - $motif->seq_region_start());
										my @newBindingScores = ();
										for (my $i = 0; $i < $motif->length + length($insertion) + 1; $i++) {
											my $altString = PatternModule::RegAnno::createComplementarySequence(substr($insertionString, 1000 + $varfeat->start() - $motif->seq_region_start() - $motif->length() + $i, $motif->seq_region_end() - $motif->seq_region_start() + 1));
											push @newBindingScores, $motif->binding_matrix->relative_affinity($altString);
										}
										my $totalScore = max(@newBindingScores)-$refScore;
										
										# Motif start
										my $StartPos = $motif->seq_region_start();
										
										# Length of motif
										my $motiflength = $motif->binding_matrix->length();
								
										# Distance from variant to start of motif
										my $DistToStart = $motif->seq_region_end() - $varfeat->start();
								
										# Size of indel
										my @splitallele = split("/", $varfeat->allele_string());
										my $SizeOfIndel = length($splitallele[1]);
								
										# Name of motif
										my @motifName = $motif->display_label();
										
										# Conservation score calculation
										my $scores = $csscore_adaptor->fetch_all_by_MethodLinkSpeciesSet_Slice($mlss, $slice, $motiflength);
										my ($consScore) = calculateConservationScore($scores);
										
										# TFBS information content
										my $infoCont = ();
										for (my $i = 0; $i < $motif->length(); $i++) {
											if ($i != $motif->length() - 1) {
												$infoCont .= $motif->binding_matrix->{'ic'}->[$i] . ",";
											}
											if ($i == $motif->length() - 1) {
												$infoCont .= $motif->binding_matrix->{'ic'}->[$i];
											}
										}
								
										unless (exists $motifScoreHash{$motifName[0]} and $motifScoreHash{$motifName[0]}{Difference} <= $totalScore) {
											$motifScoreHash{$motifName[0]}{Start} = $StartPos;
											$motifScoreHash{$motifName[0]}{Length} = $motiflength;
											$motifScoreHash{$motifName[0]}{DistToStart} = $DistToStart;
											$motifScoreHash{$motifName[0]}{SizeOfIndel} = $SizeOfIndel;
											$motifScoreHash{$motifName[0]}{Original} = $refScore;
											$motifScoreHash{$motifName[0]}{New} = max(@newBindingScores);
											$motifScoreHash{$motifName[0]}{Difference} = $totalScore;
											$motifScoreHash{$motifName[0]}{InfoContent} = $infoCont;
											$motifScoreHash{$motifName[0]}{ConsScore} = $consScore;
										}
									}
								}
				
								# For Deletions
								if (substr($varfeat->allele_string(), -1) eq "-") {
									if ($motif->strand() == 1) {
								
										# Calculating binding scores
										my $refString = $slice->seq();
										my $refScore = $motif->binding_matrix->relative_affinity($refString);
						
										# Creating a sequence including the Deletion	
										# A slice is made with 1000 bases in each direction of the motif
										my $deletionSlice = $slice_adaptor->fetch_by_region('chromosome', $varfeat->slice->seq_region_name(), $motif->seq_region_start() - 1000, $motif->seq_region_end() + 1000);
										my $deletionString = $deletionSlice->seq();
										substr($deletionString, $varfeat->start() - $motif->seq_region_start() + 1000, $varfeat->end() - $varfeat->start() + 1) = "";
										my @newBindingScores = ();						
										for (my $i = 0; $i < $motif->length + 1; $i++) {
											my $altString = substr($deletionString, 1000 + $varfeat->start() - $motif->seq_region_start() - $motif->length() + $i, $motif->length());
											push @newBindingScores, $motif->binding_matrix->relative_affinity($altString);
										}
										my $totalScore = max(@newBindingScores)-$refScore;
										
										# Motif start
										my $StartPos = $motif->seq_region_start();
										
										# Length of motif
										my $motiflength = $motif->binding_matrix->length();
					
										# Distance from variant to start of motif
										my $DistToStart = -1;
										if ($varfeat->start() - $motif->seq_region_start() > -1) {
											$DistToStart = $varfeat->start() - $motif->seq_region_start();
										}
								
										# Size of indel
										my @splitallele = split("/", $varfeat->allele_string());
										my $SizeOfIndel = length($splitallele[0]);
								
										# If the deletion is of the entire transcription factor binding motif
										if ($varfeat->start() - $motif->seq_region_start() <= 0 && $motif->seq_region_end() - $varfeat->end() <= 0) {
											$SizeOfIndel = $motif->binding_matrix->length();
										}
								
										# If the deletion start/end both lie within the transcription factor binding motif
										if ($varfeat->start() - $motif->seq_region_start() >= 0 && $motif->seq_region_end() - $varfeat->end() >= 0) {
											$SizeOfIndel = length($splitallele[0]);
										}
								
										# If the deletion start lies before the motif start, and the deletion end lies within the motif
										if ($varfeat->start() - $motif->seq_region_start() <= 0 && $motif->seq_region_end() - $varfeat->end() >= 0) {
											$SizeOfIndel = $varfeat->end() - $motif->seq_region_start() + 1;
										}
								
										# If the deletion start lies within motif, and the deletion end lies outside motif
										if ($varfeat->start() - $motif->seq_region_start() >= 0 && $motif->seq_region_end() - $varfeat->end() <= 0) {
											$SizeOfIndel = $motif->seq_region_end() - $varfeat->start() + 1;
										}								
								
										# Name of motif
										my @motifName = $motif->display_label();
										
										# Conservation score calculation
										my $scores = $csscore_adaptor->fetch_all_by_MethodLinkSpeciesSet_Slice($mlss, $slice, $motiflength);
										my ($consScore) = calculateConservationScore($scores);
										
										# TFBS information content
										my $infoCont = ();
										for (my $i = 0; $i < $motif->length(); $i++) {
											if ($i != $motif->length() - 1) {
												$infoCont .= $motif->binding_matrix->{'ic'}->[$i] . ",";
											}
											if ($i == $motif->length() - 1) {
												$infoCont .= $motif->binding_matrix->{'ic'}->[$i];
											}
										}
								
										unless (exists $motifScoreHash{$motifName[0]} and $motifScoreHash{$motifName[0]}{Difference} <= $totalScore) {
											$motifScoreHash{$motifName[0]}{Start} = $StartPos;
											$motifScoreHash{$motifName[0]}{Length} = $motiflength;
											$motifScoreHash{$motifName[0]}{DistToStart} = $DistToStart;
											$motifScoreHash{$motifName[0]}{SizeOfIndel} = $SizeOfIndel;
											$motifScoreHash{$motifName[0]}{Original} = $refScore;
											$motifScoreHash{$motifName[0]}{New} = max(@newBindingScores);
											$motifScoreHash{$motifName[0]}{Difference} = $totalScore;
											$motifScoreHash{$motifName[0]}{InfoContent} = $infoCont;
											$motifScoreHash{$motifName[0]}{ConsScore} = $consScore;
										}
									}
									if ($motif->strand == -1) {
								
										# Calculating binding scores
										my $refString = PatternModule::RegAnno::createComplementarySequence($slice->seq());
										my $refScore = $motif->binding_matrix->relative_affinity($refString);
						
										# Creating a sequence including the Deletion
										# A slice is made with 1000 bases in each direction of the motif
										my $deletionSlice = $slice_adaptor->fetch_by_region('chromosome', $varfeat->slice->seq_region_name(), $motif->seq_region_start() - 1000, $motif->seq_region_end() + 1000);
										my $deletionString = $deletionSlice->seq();
										substr($deletionString, $varfeat->start - $motif->seq_region_start() + 1000, $varfeat->end() - $varfeat->start() + 1) = "";
										my @newBindingScores = ();						
										for (my $i = 0; $i < $motif->length() + 1; $i++) {
											my $altString = PatternModule::RegAnno::createComplementarySequence(substr($deletionString, 1000 + $varfeat->start() - $motif->seq_region_start() - $motif->length() + $i, $motif->length()));
											push @newBindingScores, $motif->binding_matrix->relative_affinity($altString);
										}
										my $totalScore = max(@newBindingScores)-$refScore;
										
										# Motif start
										my $StartPos = $motif->seq_region_start();
										
										# Length of motif
										my $motiflength = $motif->binding_matrix->length();

										# Distance from variant to start of motif
										my $DistToStart = -1;
										if ($motif->seq_region_end() - $varfeat->end() > -1) {
											$DistToStart = $motif->seq_region_end() - $varfeat->end();
										}
								
										# Size of indel
										my @splitallele = split("/", $varfeat->allele_string());
										my $SizeOfIndel = length($splitallele[0]);
								
										# If the deletion is of the entire transcription factor binding motif
										if ($varfeat->start() - $motif->seq_region_start() <= 0 && $motif->seq_region_end() - $varfeat->end() <= 0) {
											$SizeOfIndel = $motif->binding_matrix->length();
										}
								
										# If the deletion start/end both lie within the transcription factor binding motif
										if ($varfeat->start() - $motif->seq_region_start() >= 0 && $motif->seq_region_end() - $varfeat->end() >= 0) {
											$SizeOfIndel = length($splitallele[0]);
										}
								
										# If the deletion start lies before the motif start, and the deletion end lies within the motif
										if ($varfeat->start() - $motif->seq_region_start() <= 0 && $motif->seq_region_end() - $varfeat->end() >= 0) {
											$SizeOfIndel = $varfeat->end() - $motif->seq_region_start() + 1;
										}
								
										# If the deletion start lies within motif, and the deletion end lies outside motif
										if ($varfeat->start() - $motif->seq_region_start() >= 0 && $motif->seq_region_end() - $varfeat->end() <= 0) {
											$SizeOfIndel = $motif->seq_region_end() - $varfeat->start() + 1;
										}		
								
										# Name of motif
										my @motifName = $motif->display_label();
										
										# Conservation score calculation
										my $scores = $csscore_adaptor->fetch_all_by_MethodLinkSpeciesSet_Slice($mlss, $slice, $motiflength);
										my ($consScore) = calculateConservationScore($scores);
										
										# TFBS information content
										my $infoCont = ();
										for (my $i = 0; $i < $motif->length(); $i++) {
											if ($i != $motif->length() - 1) {
												$infoCont .= $motif->binding_matrix->{'ic'}->[$i] . ",";
											}
											if ($i == $motif->length() - 1) {
												$infoCont .= $motif->binding_matrix->{'ic'}->[$i];
											}
										}
								
										unless (exists $motifScoreHash{$motifName[0]} and $motifScoreHash{$motifName[0]}{Difference} <= $totalScore) {
											$motifScoreHash{$motifName[0]}{Start} = $StartPos;
											$motifScoreHash{$motifName[0]}{Length} = $motiflength;
											$motifScoreHash{$motifName[0]}{DistToStart} = $DistToStart;
											$motifScoreHash{$motifName[0]}{SizeOfIndel} = $SizeOfIndel;
											$motifScoreHash{$motifName[0]}{Original} = $refScore;
											$motifScoreHash{$motifName[0]}{New} = max(@newBindingScores);
											$motifScoreHash{$motifName[0]}{Difference} = $totalScore;
											$motifScoreHash{$motifName[0]}{InfoContent} = $infoCont;
											$motifScoreHash{$motifName[0]}{ConsScore} = $consScore;
										}
									}
								}
							}
						}
					}				
				
					# TFBS names are extracted
					@keys = sort keys %motifScoreHash;
				
					# A new simple hash with only TFBS name, and Difference between old and new binding scores are present
					for my $key (@keys) {
						$motifScoreSortingHash{$key} = $motifScoreHash{$key}{Difference};
					}
				
					# Keys are sorted such that the lowest Difference-score is first
					my @DifferenceSortedKeys = sort {$motifScoreSortingHash{$a} <=> $motifScoreSortingHash{$b}} keys %motifScoreSortingHash;
				
					# The TFBS information is added to the annotation string
					for my $Rkey (@DifferenceSortedKeys) {
						# Splitting TFBS name and matrix into two
						my @TFBSnames = split(":", $Rkey);
						
						$annotationString .= $TFBSnames[0] . "\t" . $TFBSnames[-1] . "\t" . $motifScoreHash{$Rkey}{Start} . "\t" . $motifScoreHash{$Rkey}{Length} . "\t" . $motifScoreHash{$Rkey}{DistToStart} . "\t" . $motifScoreHash{$Rkey}{SizeOfIndel} . "\t" . $motifScoreHash{$Rkey}{Original} . "\t" . $motifScoreHash{$Rkey}{New} . "\t" . $motifScoreHash{$Rkey}{InfoContent} . "\t" . $motifScoreHash{$Rkey}{ConsScore}  . "\t";
					}
				
					# The empty spaces in the annotation string are filled out with NA
					if (scalar @keys != 0) {
						for (my $i = 0; $i < 10 - scalar @keys; $i++) {
							$annotationString .= "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t";
						}
					}
				}	
			
				# If the variant overlaps with no TFBS, NA is filled in into the annotation string
				if (scalar @keys == 0) {
					$annotationString .= "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" ."NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t". "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t". "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t";
				}
			}
			
			###################################################
			##### Annotating variants of the type <.....> #####
			###################################################

			# If variants are of type <....> overlapping TFBS are annotated, but not their consequences
			
			if ($col[3] =~ "<.*>" | $col[4] =~ "<.*>") {
			
				# Creates a "fake" variation feature to access variationfeature methods
				my $slice = $slice_adaptor->fetch_by_region('chromosome', $col[0]);
				$vf = PatternModule::GeneAnno::createVariationFeature($col[1], $col[3], $col[3], $slice, $variant_feature_adaptor);
				
				# Getting information on transcription factor binding motifs that are hit by variants
				my $motifFeature_adaptor = $registry->get_adaptor('Human', 'funcgen', 'motiffeature');
				my $variantSlice = $slice_adaptor->fetch_by_region('chromosome',$col[0],$col[1],$col[1]);
				my @motif_features = @{$motifFeature_adaptor->fetch_all_by_Slice($variantSlice)};
				my @keys = ();
				
				if (scalar @motif_features > 0) {
					my %motifScoreHash = ();
				
					for my $motif (@motif_features) {
						if ($vf->start() <= $motif->seq_region_end() && $motif->seq_region_start() <= $vf->end()) {
							
							# Length of motif
							my $motiflength = $motif->binding_matrix->length();
							
							# Conservation score calculation
							my $slice = $slice_adaptor->fetch_by_region('chromosome', $variantSlice->seq_region_name(), $motif->seq_region_start(), $motif->seq_region_end());
							my $scores = $csscore_adaptor->fetch_all_by_MethodLinkSpeciesSet_Slice($mlss, $slice, $motiflength);
							my ($consScore) = calculateConservationScore($scores);
							
							# Name of motif
							my @motifName = $motif->display_label();
							
							# Motif start
							my $StartPos = $motif->seq_region_start();
							
							# Placing motif names in hash
							unless (exists $motifScoreHash{$motifName[0]}) {
								$motifScoreHash{$motifName[0]}{Start} = $StartPos;
								$motifScoreHash{$motifName[0]}{ConsScore} = $consScore;								
							}
						}
					}
					
					# TFBS names are extracted
					@keys = sort keys %motifScoreHash;
					
					# The TFBS information is added to the annotation string
					for my $key (@keys) {
						# Splitting TFBS name and matrix into two
						my @TFBSnames = split(":", $key);
						
						$annotationString .= $TFBSnames[0] . "\t" . $TFBSnames[-1] . "\t" . $motifScoreHash{$key}{Start} . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . $motifScoreHash{$key}{ConsScore} . "\t";
					}
				
					# The empty spaces in the annotation string are filled out with NA
					if (scalar @keys != 0) {
						for (my $i = 0; $i < 10 - scalar @keys; $i++) {
							$annotationString .= "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t";
						}
					}	
				}
				
				# If the variant overlaps with no TFBS, NA is filled in into the annotation string
				if (scalar @keys == 0) {
					$annotationString .= "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" ."NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t". "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t". "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t";
				}		
			}
			
			print OUTPUT $annotationString;
			print OUTPUT "\n";	
    	}
	}

# Calculates average conservation score for a slice given an array of scores
sub calculateConservationScore {
	my ($scores) = @_;
	my $sum = 0;
	my $count = 0;
	my $diff_score;
	foreach my $score (@$scores) {
		if (defined $score->diff_score) {
			$sum += $score->diff_score;
			$count++;
		}
	}	
	if ($count > 0){
		$diff_score = $sum / $count;
	}
	if (defined $diff_score){
		return $diff_score;
	}
	else {
		return "NA";
	}
}

