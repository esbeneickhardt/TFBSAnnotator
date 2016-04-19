package PatternModule::RegAnno;

use strict;
use warnings;
use List::Util qw(min);
use List::Util qw(max);

# Regulatory consequences
my %consequecetable = (
	"5_prime_UTR_variant" => 1,
	"3_prime_UTR_variant" => 2,
	"TFBS_ablation" => 3,
	"TFBS_amplification" => 4,
	"TF_binding_site_variant" => 5,
	"regulatory_region_ablation" => 6,
	"regulatory_region_amplification" => 7,
	"regulatory_region_variant" => 8,
);
my %rconsequecetable = reverse %consequecetable;

=head2 createVariantRegulatoryHash

  Args       : None
  Example    : my %variantInformation = createVariantRegulatoryHash()
  Description: Creates a hash that stores information on the variant's effects on 
               regulatory elements
  Returntype : Hash

=cut

sub createVariantRegulatoryHash {
	my %data = ('ENSREGANNO_VARIANT_OVERLAPS_WITH_REGULATORY'	=> "NA",
				'ENSREGANNO_REGULATORY_FEATURE'					=> "NA",
				'ENSREGANNO_REGULATORY_FEATURE_TYPE'			=> "NA",
				'ENSREGANNO_DISTANCE_TO_START' 					=> ".",
				'ENSREGANNO_DISTANCE_TO_END' 					=> ".",
				'ENSREGANNO_CONSEQUENCE' 						=> "NA",
				'ENSREGANNO_MULTICELL_TFBS_HIT' 				=> "NA",
				'ENSREGANNO_MULTICELL_TFBS' 					=> "NA",
				'ENSREGANNO_MULTICELL_MULTITFBS'				=> "NA",
				'ENSREGANNO_MULTICELL_TFBS_SCORE' 				=> ".",
				'ENSREGANNO_MULTICELL_MULTITFBS_SCORE'			=> ".",
				'ENSREGANNO_MULTICELL_EVIDENCE' 				=> "NA",
				'ENSREGANNO_MULTICELL_MULTIEVIDENCE'			=> "NA",
				'ENSREGANNO_NEUROCELL_TFBS_HIT' 				=> "NA",
				'ENSREGANNO_NEUROCELL_TFBS' 					=> "NA",
				'ENSREGANNO_NEUROCELL_MULTITFBS'				=> "NA",
				'ENSREGANNO_NEUROCELL_TFBS_SCORE' 				=> ".",
				'ENSREGANNO_NEUROCELL_MULTITFBS_SCORE'			=> ".",
				'ENSREGANNO_NEUROCELL_EVIDENCE' 				=> "NA",
				'ENSREGANNO_NEUROCELL_MULTIEVIDENCE'			=> "NA",
				'ENSREGANNO_IMMUNOCELL_TFBS_HIT' 				=> "NA",
				'ENSREGANNO_IMMUNOCELL_TFBS' 					=> "NA",
				'ENSREGANNO_IMMUNOCELL_MULTITFBS'				=> "NA",
				'ENSREGANNO_IMMUNOCELL_TFBS_SCORE' 				=> ".",
				'ENSREGANNO_IMMUNOCELL_MULTITFBS_SCORE'			=> ".",
				'ENSREGANNO_IMMUNOCELL_EVIDENCE' 				=> "NA",
				'ENSREGANNO_IMMUNOCELL_MULTIEVIDENCE'			=> "NA",
	);	
	return %data;
}

=head2 doesVariantOverlapWithRegulatoryFeature

  Args       : VariationFeature
  Example    : my ($answer, $reg_feature, $reg_feature_consequence) = doesVariantOverlapWithRegulatoryFeature($vf)
  Description: Tells whether a VariationFeature overlaps with a RegulatoryFeature and 
               the regulatory feature is returned together with the consequence of the 
               variant
  Returntype : Array with YES/NO, name of the overlapping regulatory feature and 
               consequence of the variant

=cut

sub doesVariantOverlapWithRegulatoryFeature {
	my ($vf) = @_;	
	my $answer = "NO";
	my $reg_feature = "NA";
	my $reg_feature_consequence = "NA";
	
	if (scalar @{$vf->get_all_RegulatoryFeatureVariations()} > 0) {
		$answer = "YES";
		my @consequenceValues = ();
		
		# Looking for specific regulatory consequence
		for my $consequence (@{$vf->consequence_type()}) {
			if (exists $consequecetable{$consequence}) {
				push @consequenceValues, $consequecetable{$consequence};
			} 
		}
		if (scalar @consequenceValues > 0) {
			my $consequenceValue = min(@consequenceValues);
			$reg_feature_consequence = $rconsequecetable{$consequenceValue};
		}
		
		# If none of the consequences from the table are annotated, the standard annotation is used
		if ($reg_feature_consequence eq "NA") {
			$reg_feature_consequence = ${$vf->get_all_RegulatoryFeatureVariations()}[0]->most_severe_OverlapConsequence->SO_term();
		}
		
		$reg_feature = ${$vf->get_all_RegulatoryFeatureVariations()}[0]->regulatory_feature_stable_id();
	}
	return ($answer, $reg_feature, $reg_feature_consequence);
}

=head2 findClosestRegulatoryFeature

  Args       : "VariationFeature", "RegulatoryFeature Adaptor", "Slice Adaptor", "Interval"
  Example    : my ($closestRegFeat, $vartype) = findClosestRegulatoryFeature($vf, $regfeat_adaptor, $slice_adaptor, $regint)
  Description: Finds closest regulatory feature within the specified number of bases 
               up-/downstream and tells whether it is upstream or downstream
  Returntype : Array with RegulatoryFeature stable id and the consequence of the variant

=cut

sub findClosestRegulatoryFeature {
	my ($vf, $regfeat_adaptor, $slice_adaptor, $regint) = @_;	
	my $slice = $slice_adaptor->fetch_by_region('chromosome', $vf->slice->seq_region_name, $vf->start - $regint, $vf->start + $regint);
	my @features = @{$regfeat_adaptor->fetch_all_by_Slice($slice)};
	my %featureHash = ();
	my %upstreamDownstreamHash = ();
	
	for my $feature (@features) {
		my $distToRegFeat = $vf->start() - $feature->seq_region_start();
		my $regFeatID = $feature->stable_id();
		$featureHash{$regFeatID} = abs($distToRegFeat);
		if ($distToRegFeat > 0) {
			$upstreamDownstreamHash{$regFeatID} = "upstream_regulatory_region_variant";
		}
		if ($distToRegFeat < 0) {
			$upstreamDownstreamHash{$regFeatID} = "downstream_regulatory_region_variant";
		}
	}
	
	my $closestRegFeat = (sort {$featureHash{$a} <=> $featureHash{$b}} keys %featureHash)[0];
	
	if ($closestRegFeat) {
		return ($closestRegFeat, $upstreamDownstreamHash{$closestRegFeat});
	}
	else {
		return ("NA","NA");
	}
}

=head2 findDistanceToStartAndEndOfRegFeature

  Args       : "RegulatoryFeature", "VariationFeature"
  Example    : my ($distToStart, $distToEnd) = findDistanceToStartAndEndOfRegFeature($reg_feature, $vf)
  Description: Finds distance from variant to start and end of RegulatoryFeature          
  Returntype : Array with distance from variant to start and end of RegulatoryFeature

=cut

sub findDistanceToStartAndEndOfRegFeature {
	my ($reg_feature, $vf) = @_;	
	my $distToStart = ".";
	my $distToEnd = ".";
	
	# All regulatory features are placed on strand 0
	
	# Insertions and SNPs
	if ($vf->end() - $vf->start() < 1) {
		$distToStart = $vf->start() - $reg_feature->seq_region_start();
		$distToEnd = $reg_feature->seq_region_end() - $vf->start();
	}
	
	# Deletions
	if ($vf->end() - $vf->start() > 0) {
		$distToStart = $vf->start() - $reg_feature->seq_region_start();
		$distToEnd = $reg_feature->seq_region_end() - $vf->end();
	}
	return ($distToStart, $distToEnd);	
}

# Hash of tissues that each cell type comes from:
# http://www.ensembl.org/Homo_sapiens/Experiment/Sources?db=funcgen;ex=all;fdb=funcgen;r=17:46617590-46621119
my %tissuetable = (
	"HepG2" => "Liver",
	"IMR90" => "Lung",
	"K562" => "Bone_Marrow",
	"GM12878" => "Immuno",
	"HeLa-S3" => "Cervix",
	"HUVEC" => "Umbilical_Cord",
	"NHEK" => "Skin",
	"H1ESC" => "Embryo",
	"A549" => "Lung",
	"HMEC" => "Skin",
	"HSMM" => "Muscle",
	"HSMMtube" => "Muscle",
	"K562b" => "Bone_Marrow",
	"NH-A" => "Neuro",
	"NHDF-AD" => "Skin",
	"NHLF" => "Lung",
	"Osteobl" => "Bone",
	"GM06990" => "Immuno",
	"CD4" => "Immuno",
	"MultiCell" => "Multiple_Tissues",
	"NA" => "NA",
);	

# Motif Features represent short genomic regions where the Transcription Factor (TF) protein is thought to be directly interacting with DNA (TF binding sites). 
# These are defined by combining binding matrices and annotated features. More information on how they are defined can be found on the RegulatoryBuild page.
# sub getTranscriptionFactorBindingMotifs {
# 	my ($per_cell_reg_features, $vf, $tissue) = @_;
# 	my @TFBS = ();
# 	my $MultiTFBS = "NA";
# 	
# 	for my $rf (@{$per_cell_reg_features}) {
# 		my @motif_features = @{$rf->regulatory_attributes('motif')};
# 		
# 		if ($tissuetable{$rf->cell_type->name} eq $tissue) {
# 			for my $motif (@motif_features) {
# 				if ($vf->start <= $motif->end && $motif->start <= $vf->end) {
# 					push @TFBS, $motif->display_label();
# 				}
# 			}
# 		}
# 	}
# 	
# 	@TFBS = @{stripText(\@TFBS, ":")};
# 	@TFBS = removeDuplicatesFromList(@TFBS);
# 
# 	if (scalar @TFBS > 1) {
# 		$MultiTFBS = join(",", @TFBS);
# 		@TFBS = "Multi";
# 	}
# 
# 	if (scalar @TFBS == 0) {
# 		$MultiTFBS = "NA";
# 		@TFBS = "NA";
# 	}
# 	
# 	return($TFBS[0], $MultiTFBS);
# }

=head2 getNearestTranscriptionFactorBindingMotif

  Args       : "VariationFeature", "MotifFeature Adaptor", "Slice Adaptor", "Tissue", "Interval"
  Example    : my ($closestMotifFeat) = getNearestTranscriptionFactorBindingMotif($vf, $motifFeature_adaptor, $slice_adaptor, $tissue, $tfbsint)
  Description: Finds closest transcription factor binding motif within the specified number of bases 
               up-/downstream within the chosen tissue type
  Returntype : Motif

=cut

sub getNearestTranscriptionFactorBindingMotif {
	my ($vf, $motifFeature_adaptor, $slice_adaptor, $tissue, $tfbsint) = @_;
	my $slice = $slice_adaptor->fetch_by_region('chromosome', $vf->slice->seq_region_name(), $vf->start() - $tfbsint, $vf->start() + $tfbsint);
	my @motif_features = @{$motifFeature_adaptor->fetch_all_by_Slice($slice)};
	my %featureHash = ();
	
	for my $feature (@motif_features) {
		my $distToRegFeat = $vf->start() - $feature->seq_region_start();	
		my $regFeatID = $feature->display_label();
		my @associatedAnnoFeat = @{$feature->associated_annotated_features()};			
		
		if ($tissue eq "Multiple_Tissues") {
			$featureHash{$regFeatID} = abs($distToRegFeat);
		}
		
		if ($tissue ne "Multiple_Tissues") {
			for my $assAnno (@associatedAnnoFeat) {
				if ($tissuetable{$assAnno->cell_type->name()} eq $tissue) {
					$featureHash{$regFeatID} = abs($distToRegFeat);
				}
			}
		}		
	}

	my $closestMotifFeat = (sort {$featureHash{$a} <=> $featureHash{$b}} keys %featureHash)[0];

	if ($closestMotifFeat) {
		my @splitClostestMotifFeat = split(":", $closestMotifFeat);
		return $splitClostestMotifFeat[0];
	}
	else {
		return ("NA");
	}
}

=head2 getTranscriptionFactorBindingScore

  Args       : "Per Cell RegulatoryFeatures", "VariationFeature", "Slice Adaptor", "Tissue"
  Example    : my ($TFBS, $TFBSs, $Score, $Scores) = getTranscriptionFactorBindingScore($per_cell_reg_features, $vf, $slice_adaptor, $tissue)
  Description: Finds transcription factor binding motifs that the variant overlaps with.
               If there are multiple overlaps, the first output will be "Multi", and the second a 
               list of Motifs. A relative binding score is calculated for each motif that 
               is hit, and is the difference in affinity scores between the reference sequence
               and the alternative sequence. A sliding window is examined around the variant 
               for SNVs, deletions and insertions for alternative binding sites that may bind 
               better than the initially interrupted motif. If a TF has multiple binding 
               motifs the most severely hit motif is chosen and scored
  Returntype : Array with four elements, single TFBS, list of TFBSs, single Score, list of Scores

=cut


sub getTranscriptionFactorBindingScore {
	my ($per_cell_reg_features, $vf, $slice_adaptor, $tissue) = @_;
	my %motifScoreHash = ();	
	for my $rf (@{$per_cell_reg_features}) {
		my @motif_features = @{$rf->regulatory_attributes('motif')};
		
		if ($tissuetable{$rf->cell_type->name()} eq $tissue) {
			for my $motif (@motif_features) {
				if ($vf->start() <= $motif->seq_region_end() && $motif->seq_region_start() <= $vf->end()) {
					my $slice = $slice_adaptor->fetch_by_region('chromosome', $vf->slice->seq_region_name(), $motif->start(), $motif->end());					
														
					# For SNPs
					if ($vf->allele_string() !~ /-/) {
						if ($motif->strand() == 1) {
							# Calculating binding score
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
														
							# Name of motif
							my @motifName = split(":", $motif->display_label());
							
							# Placing information into Hash. If the same motif has two scores for a variant, the most severe score is stored
							unless (exists $motifScoreHash{$motifName[0]} and $motifScoreHash{$motifName[0]} <= $totalScore) {
								$motifScoreHash{$motifName[0]} = $totalScore;
							}
						}
						
						if ($motif->strand() == -1) {
							# Calculating bindin score
							my $refString = createComplementarySequence($slice->seq());
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
							
							# Name of motif
							my @motifName = split(":", $motif->display_label());
							
							# Placing information into Hash. If the same motif has two scores for a variant, the most severe score is stored
							unless (exists $motifScoreHash{$motifName[0]} and $motifScoreHash{$motifName[0]} <= $totalScore) {
								$motifScoreHash{$motifName[0]} = $totalScore;
							}
						}	
					}
				
					# For Insertions
					if (substr($vf->allele_string(), 0, 1) eq "-") {
						if ($motif->strand() == 1) {
							my $refString = $slice->seq();
							my $refScore = $motif->binding_matrix->relative_affinity($refString);
							
							# Creating a sequence including the insertion
							
							# A slice is made with 1000 bases in each direction of the motif
							my $insertionSlice = $slice_adaptor->fetch_by_region('chromosome', $vf->slice->seq_region_name(), $motif->start() - 1000, $motif->end() + 1000);
							my $insertionString = $insertionSlice->seq();
							my $insertion = substr($vf->allele_string(), 2);
							$insertionString = substr($insertionString, 0, 1000 + $vf->start() - $motif->start()) . $insertion . substr($insertionString, 1000 + $vf->start() - $motif->start());							
							my @newBindingScores = ();
							for (my $i = 0; $i < $motif->length + length($insertion) + 1; $i++) {
								my $altString = substr($insertionString, 1000 + $vf->start() - $motif->start() - $motif->length() + $i, $motif->end() - $motif->start() + 1);
								push @newBindingScores, $motif->binding_matrix->relative_affinity($altString);
							}

							my $totalScore = max(@newBindingScores)-$refScore;
							
							my @motifName = split(":", $motif->display_label());
							unless (exists $motifScoreHash{$motifName[0]} and $motifScoreHash{$motifName[0]} <= $totalScore) {
								$motifScoreHash{$motifName[0]} = $totalScore;
							}
						}
						
						if ($motif->strand == -1) {
							my $refString = createComplementarySequence($slice->seq);
							my $refScore = $motif->binding_matrix->relative_affinity($refString);
							
							# Creating a sequence including the insertion
							
							# A slice is made with 1000 bases in each direction of the motif
							my $insertionSlice = $slice_adaptor->fetch_by_region('chromosome', $vf->slice->seq_region_name(), $motif->start() - 1000, $motif->end() + 1000);
							my $insertionString = $insertionSlice->seq();
							my $insertion = substr($vf->allele_string(), 2);
							$insertionString = substr($insertionString, 0, 1000 + $vf->start() - $motif->start()) . $insertion . substr($insertionString, 1000 + $vf->start() - $motif->start());
							
							my @newBindingScores = ();
							for (my $i = 0; $i < $motif->length + length($insertion) + 1; $i++) {
								my $altString = createComplementarySequence(substr($insertionString, 1000 + $vf->start() - $motif->start() - $motif->length() + $i, $motif->end() - $motif->start() + 1));
								push @newBindingScores, $motif->binding_matrix->relative_affinity($altString);
							}

							my $totalScore = max(@newBindingScores)-$refScore;
							
							my @motifName = split(":", $motif->display_label());
							unless (exists $motifScoreHash{$motifName[0]} and $motifScoreHash{$motifName[0]} <= $totalScore) {
								$motifScoreHash{$motifName[0]} = $totalScore;
							}
						}
					}
				
					# For Deletions
					if (substr($vf->allele_string(), -1) eq "-") {
						if ($motif->strand() == 1) {
							my $refString = $slice->seq();
							my $refScore = $motif->binding_matrix->relative_affinity($refString);
						
							# Creating a sequence including the Deletion
													
							# A slice is made with 1000 bases in each direction of the motif
							my $deletionSlice = $slice_adaptor->fetch_by_region('chromosome', $vf->slice->seq_region_name(), $motif->start() - 1000, $motif->end() + 1000);
							my $deletionString = $deletionSlice->seq();
							substr($deletionString, $vf->start() - $motif->start() + 1000, $vf->end() - $vf->start() + 1) = "";
							
							my @newBindingScores = ();						
							for (my $i = 0; $i < $motif->length + 1; $i++) {
								my $altString = substr($deletionString, 1000 + $vf->start() - $motif->start() - $motif->length() + $i, $motif->length());
								push @newBindingScores, $motif->binding_matrix->relative_affinity($altString);
							}
							
							my $totalScore = max(@newBindingScores)-$refScore;
							
							my @motifName = split(":", $motif->display_label());
							unless (exists $motifScoreHash{$motifName[0]} and $motifScoreHash{$motifName[0]} <= $totalScore) {
								$motifScoreHash{$motifName[0]} = $totalScore;
							}
						}
						if ($motif->strand == -1) {
							my $refString = createComplementarySequence($slice->seq());
							my $refScore = $motif->binding_matrix->relative_affinity($refString);
						
							# Creating a sequence including the Deletion
													
							# A slice is made with 1000 bases in each direction of the motif
							my $deletionSlice = $slice_adaptor->fetch_by_region('chromosome', $vf->slice->seq_region_name(), $motif->start() - 1000, $motif->end() + 1000);
							my $deletionString = $deletionSlice->seq();
							substr($deletionString, $vf->start - $motif->start() + 1000, $vf->end() - $vf->start() + 1) = "";
							
							my @newBindingScores = ();						
							for (my $i = 0; $i < $motif->length() + 1; $i++) {
								my $altString = createComplementarySequence(substr($deletionString, 1000 + $vf->start() - $motif->start() - $motif->length() + $i, $motif->length()));
								push @newBindingScores, $motif->binding_matrix->relative_affinity($altString);
							}
							
							my $totalScore = max(@newBindingScores)-$refScore;
							
							my @motifName = split(":", $motif->display_label());
							unless (exists $motifScoreHash{$motifName[0]} and $motifScoreHash{$motifName[0]} <= $totalScore) {
								$motifScoreHash{$motifName[0]} = $totalScore;
							}
						}
					}
				}
			}
		}
	}
	
	my @keys = sort keys %motifScoreHash;
	
	if (scalar @keys == 1) {
		return ($keys[0], "NA", $motifScoreHash{$keys[0]}, ".");
	}
	
	if (scalar @keys > 1) {
		my @scores = ();
		foreach my $key (@keys) {
			push @scores, $motifScoreHash{$key};
		}
		return ("Multi", join(",", @keys), ".", join(",", @scores));
	}
	
	else {
		return ("NA", "NA", ".", ".");
	}
}

=head2 getEvidenceType

  Args       : "Per Cell RegulatoryFeatures", "VariationFeature", "Tissue"
  Example    : my ($evidence, $multiEvidence) = getEvidenceType($per_cell_reg_features, $vf, $tissue)
  Description: Finds evidence for the RegulatoryFeatures. If there is more that one type of
               evidence, the first output is set as "Multi", and a list of evidence is 
               returned in the second output as a string with a comma separator
  Returntype : Array with singleEvidence and multiEvidence

=cut

# Regulatory Features are built based on results from experiments like Dnase1 sensitivity assays (Dnase-Seq) to detect regions of open chromatin, or TF binding assays, 
# like Chromatin ImmunoPrecipitation (ChIP) coupled with high throughput sequencing (ChIP-Seq). Results from these experiments are stored as Annotated Features.
sub getEvidenceType {
	my ($per_cell_reg_features, $vf, $tissue) = @_;
	my @Evidence = ();
	my $MultiEvidence = "NA";
	
	for my $rf (@{$per_cell_reg_features}) {
		my @annotated_features = @{$rf->regulatory_attributes('annotated')};
		
		if ($tissuetable{$rf->cell_type->name()} eq $tissue) {
			for my $annotated (@annotated_features) {
				if ($vf->start() <= $annotated->end() && $annotated->start() <= $vf->end()) {
					push @Evidence, $annotated->display_label();
				}
			}
		}
	}
	
	@Evidence = @{stripText(\@Evidence, " - ")};
	@Evidence = removeDuplicatesFromList(@Evidence);
		
	if (scalar @Evidence > 1) {
		$MultiEvidence = join(",", @Evidence);
		@Evidence = "Multi";
	}

	if (scalar @Evidence == 0) {
		$MultiEvidence = "NA";
		@Evidence = "NA";
	}

	return($Evidence[0], $MultiEvidence);	
}

=head2 createComplementarySequence

  Args       : "String"
  Example    : my ($complimentaryString) = createComplementarySequence($string)
  Description: Creates a complimentary string
  Returntype : String that is complimentary

=cut

sub createComplementarySequence {
	my ($string) = @_;
	my $length = length $string;
	
	for (my $i = 0; $i < $length; $i++) {
		if (substr($string, $i, 1) eq "A") {
			substr($string, $i, 1) = "T";
			next;
		}
		if (substr($string, $i, 1) eq "T") {
			substr($string, $i, 1) = "A";
			next;
		}
		if (substr($string, $i, 1) eq "G") {
			substr($string, $i, 1) = "C";
			next;
		}
		if (substr($string, $i, 1) eq "C") {
			substr($string, $i, 1) = "G";
			next;
		}
	}
	
	$string = reverse $string;
	
	return $string;
}

=head2 stripText

  Args       : "String", "Separator"
  Example    : my ($strippedString) = stripText($string)
  Description: Creates a list from the string and takes out first element
  Returntype : String that is simplified

=cut

sub stripText {
	my ($textArray, $divider) = @_;
	my @splitArray = ();
	
	for my $text (@{$textArray}) {
		my @splitText = split($divider, $text);
		push @splitArray, $splitText[0];
	}
	
	return \@splitArray;
}

=head2 removeDuplicatesFromList

  Args       : "Array"
  Example    : my (@noDuplicateArray) = removeDuplicatesFromList(@Array)
  Description: Removes duplicates from array
  Returntype : Array with no duplicates

=cut

sub removeDuplicatesFromList {
    my %seen;
    grep !$seen{$_}++, @_;
}

return 1;