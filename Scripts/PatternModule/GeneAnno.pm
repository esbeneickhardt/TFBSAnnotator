package PatternModule::GeneAnno;

use strict;
use warnings;
use Bio::EnsEMBL::Variation::VariationFeature;

# Hash table with ranks of all the consequence types in the SO-format from 
# http://www.ensembl.org/info/genome/variation/predicted_data.html#consequences
my %consequecetable = (
	"transcript_ablation" => 1,
	"splice_donor_variant" => 2,
	"splice_acceptor_variant" => 3,
	"stop_gained" => 4,
	"frameshift_variant" => 5,
	"stop_lost" => 6,
	"initiator_codon_variant" => 7,
	"transcript_amplification" => 8,
	"inframe_insertion" => 9,
	"inframe_deletion" => 10,
	"missense_variant" => 11,
	"splice_region_variant" => 12,
	"incomplete_terminal_codon_variant" => 13,
	"stop_retained_variant" => 14,
	"synonymous_variant" => 15,
	"coding_sequence_variant" => 16,
	"mature_miRNA_variant" => 17,
	"5_prime_UTR_variant" => 18,
	"3_prime_UTR_variant" => 19,
	"non_coding_transcript_exon_variant" => 20,
	"non_coding_exon_variant" => 20, 
	"intron_variant" => 21,
	"NMD_transcript_variant" => 22,
	"non_coding_transcript_variant" => 23,
	"nc_transcript_variant" => 23,
	"upstream_gene_variant" => 24,
	"downstream_gene_variant" => 25,
	"TFBS_ablation" => 26,
	"TFBS_amplification" => 27,
	"TF_binding_site_variant" => 28,
	"regulatory_region_ablation" => 29,
	"regulatory_region_amplification" => 30,
	"regulatory_region_variant" => 31,
	"feature_elongation" => 32,
	"feature_truncation" => 33,
	"intergenic_variant" => 34,
);

# Ensembl biotypes are found on: http://www.ensembl.org/Help/Faq?id=468
# Protein coding: IG_C_gene, IG_D_gene, IG_gene, IG_J_gene, IG_LV_gene, IG_M_gene, IG_V_gene, IG_Z_gene, nonsense_mediated_decay, nontranslating_CDS, non_stop_decay, polymorphic, polymorphic_pseudogene, protein_coding, TR_C_gene, TR_D_gene, TR_gene, TR_J_gene, TR_V_gene
# Pseudogene: disrupted_domain, IG_C_pseudogene, IG_J_pseudogene, IG_pseudogene, IG_V_pseudogene, processed_pseudogene, pseudogene, transcribed_processed_pseudogene, transcribed_unitary_pseudogene, transcribed_unprocessed_pseudogene, translated_processed_pseudogene, TR_J_pseudogene, TR_pseudogene, TR_V_pseudogene, unitary_pseudogene, unprocessed_pseudogene
# Long noncoding: 3prime_overlapping_ncrna, ambiguous_orf, antisense, antisense_RNA, lincRNA, ncrna_host, non_coding, processed_transcript, retained_intron, sense_intronic, sense_overlapping
# Short noncoding: miRNA, miRNA_pseudogene, misc_RNA, misc_RNA_pseudogene, Mt_rRNA, Mt_tRNA, Mt_tRNA_pseudogene, ncRNA, ncRNA_pseudogene, rRNA, rRNA_pseudogene, scRNA, scRNA_pseudogene, snlRNA, snoRNA, snoRNA_pseudogene, snRNA, snRNA_pseudogene, tRNA, tRNA_pseudogene

# Two hash tables with long/short noncoding RNA are created to look up transcript biotypes in
my $shortnoncoding = "miRNA,miRNA_pseudogene,misc_RNA,misc_RNA_pseudogene,Mt_rRNA,Mt_tRNA,Mt_tRNA_pseudogene,ncRNA,ncRNA_pseudogene,rRNA,rRNA_pseudogene,scRNA,scRNA_pseudogene,snlRNA,snoRNA,snoRNA_pseudogene,snRNA,snRNA_pseudogene,tRNA,tRNA_pseudogene";
my $longnoncoding = "3prime_overlapping_ncrna,ambiguous_orf,antisense,antisense_RNA,lincRNA,ncrna_host,non_coding,processed_transcript,retained_intron,sense_intronic,sense_overlapping";

my %shortNoncodingTable;
for my $item (split(",", $shortnoncoding)) {
	$shortNoncodingTable{$item} = 1;
}

my %longNoncodingTable;
for my $item (split(",", $longnoncoding)) {
	$longNoncodingTable{$item} = 1;
}

=head2 createVariationFeature

  Args       : "Position", "Reference Sequence", "Alternative Sequence", 
               "Slice with Relevant Chromosome", "VariationFeature Adaptor"
  Example    : my $vf = createVariationFeature($pos, $ref, $alt, $slice, $vfa)
  Description: Creates a VariationFeature for SNPs and INDELs using a trimmed version of 
               allelic strings (e.g. start=100 string=G/GTT is converted into start=101
               string=-/TT)
  Returntype : VariationFeature

=cut

sub createVariationFeature {
	my ($pos, $ref, $alt, $slice, $vfa) = @_;
	my $variationfeature = ();
	
	# Note: ENSEMBL uses a trimmed version of allelic strings compared to vcf
	# VCF: G/GTT ENSEMBL: -/TT
	
	# Catching SNPs
	if (length($ref) == 1 && length($alt) == 1) {
		my $string = $ref . "/" . $alt;
		$variationfeature = Bio::EnsEMBL::Variation::VariationFeature->new(
			-start   => $pos,
			-end     => $pos,
			-strand  => 1,
			-allele_string => $string,
			-slice   => $slice,
			-adaptor => $vfa
		);
	}

	# Catching substitutions (e.g. AA/GC)
	if (length($ref) == length($alt) && length($ref) > 1) {
		my $string = $ref . "/" . $alt;
		$variationfeature = Bio::EnsEMBL::Variation::VariationFeature->new(
			-start   => $pos,
			-end     => $pos + length($alt) - 1,
			-strand  => 1,
			-allele_string => $string,
			-slice   => $slice,
			-adaptor => $vfa
		);
	}

	# Catching Insertions
	if (length($ref) < length($alt)) {	
	
		# Creating trimmed string in ENSEMBL format -/GTTT
		my $string = "-" . "/" . substr($alt, 2-1, length($alt) - 1);
		$variationfeature = Bio::EnsEMBL::Variation::VariationFeature->new(
			-start   => $pos + 1,
			-end     => $pos,
			-strand  => 1,
			-allele_string => $string,
			-slice   => $slice,
			-adaptor => $vfa
		);
	}
	

	# Catching Deletions
	if (length($ref) > length($alt)) {
		
		# Creating trimmed string in ENSEMBL format GTTT/-
		my $string = substr($ref, 2-1, length($ref) - 1) . "/" . "-";
		$variationfeature = Bio::EnsEMBL::Variation::VariationFeature->new(
			-start   => $pos + 1,
			-end     => $pos + length($ref) - length($alt),
			-strand  => 1,
			-allele_string => $string,
			-slice   => $slice,
			-adaptor => $vfa
		);
	}
	return $variationfeature;
}

=head2 createMultiAllelicVariationFeature

  Args       : "Position", "Reference Sequence", "Alternative Sequence", 
               "Slice with Relevant Chromosome", "VariationFeature Adaptor"
  Example    : my $vf = createMultiAllelicVariationFeature($pos, $ref, $alt, $slice, $vfa)
  Description: Creates a VariationFeature for multi allelic variants. The allelic string is
               in the form Ref/Alt1/Alt2/Alt3 (e.g. ATT/A/ATTT/AGT). The start and end 
               positions of the VariationFeature are relative to the reference sequence, though 
               the precise position of a specific variant may differ (e.g. start=100 
               string=ATT/A/ATTT/AGT for Alt3 the positions would be start=101 string T/G).
               
               NB: in this example Alt3 is formatted in this way because it is represented as
               part of a multisample-called variant, i.e. with multiple alleles. But the variant
               itself is just a SNP.
               
  Returntype : VariationFeature

=cut

sub createMultiAllelicVariationFeature {
	my ($pos, $ref, $alt, $slice, $vfa) = @_;
	my $variationfeature = ();
	my @altSplit = split(",", $alt);
	my $string = $ref."/".join("/", @altSplit);
	$variationfeature = Bio::EnsEMBL::Variation::VariationFeature->new(
		-start   => $pos,
		-end     => $pos + length($ref) - 1,
		-strand  => 1,
		-allele_string => $string,
		-slice   => $slice,
		-adaptor => $vfa
	);
	return $variationfeature;
}

=head2 splitMultiAllelicVariationFeature

  Args       : "Multi Allelic VariationFeature", "Slice with Relevant Chromosome", 
               "Variation Feature Adaptor"
  Example    : my @variationFeatures = splitMultiAllelicVariationFeature($vf, $slice, $vfa)
  Description: Splits a multi allelic VariationFeature into several simple VariationFeatures
               where the start/end positions have been adjusted and where the allelic 
               strings have been simplified (e.g. Multi allelic: start=100 end=103 
               string=ATT/A/ATTT/AGT is now VF1: start=101 end=102 string=TT/-, VF2: 
               start=103 end=102 string=-/T, VF3: start= 101 end=101 string=T/G)
  Returntype : Array with VariationFeatures

=cut

sub splitMultiAllelicVariationFeature {
	my ($vf, $slice, $vfa) = @_;
	my $allele_string = $vf->allele_string();
	my @alleleSplit = split("/", $allele_string);
	my @variationFeatureList = ();
	
	for (my $i = 0; $i < scalar @alleleSplit - 1; $i++) {
		my $ref = $alleleSplit[0];
		my $alt = $alleleSplit[$i + 1];
		my $pos = $vf->start();
	
		my @cleaningInfoSNP = PatternModule::CleaningVariants::cleanSNP($alleleSplit[0], $alleleSplit[$i + 1]);
		my @cleaningInfoDEL = PatternModule::CleaningVariants::cleanDELETIONS($alleleSplit[0], $alleleSplit[$i + 1]);
		my @cleaningInfoINS = PatternModule::CleaningVariants::cleanINSERTIONS($alleleSplit[0], $alleleSplit[$i + 1]);
		
		if (scalar @cleaningInfoSNP == 3) {
			$ref = $cleaningInfoSNP[0];
			$alt = $cleaningInfoSNP[1];
			$pos = $vf->start + $cleaningInfoSNP[2];
		}
		
		if (scalar @cleaningInfoDEL == 3) {
			$ref = $cleaningInfoDEL[0];
			$alt = $cleaningInfoDEL[1];
			$pos = $vf->start + $cleaningInfoDEL[2];
		}
		
		if (scalar @cleaningInfoINS == 3) {
			$ref = $cleaningInfoINS[0];
			$alt = $cleaningInfoINS[1];
			$pos = $vf->start + $cleaningInfoINS[2];
		}
		
		push @variationFeatureList, createVariationFeature($pos, $ref, $alt, $slice, $vfa);
	}
	
	return @variationFeatureList;
}

=head2 createVariantProteinCodingHash

  Args       : None
  Example    : my %variantInformation = createVariantProteinCodingHash()
  Description: Creates a hash that stores annotation information from protein coding
               variants
  Returntype : Hash

=cut

sub createVariantProteinCodingHash {
	my %data = ('ENSPCANNO_VARIANT_OVERLAPS_WITH_GENE' 		=> "NA", 
				'ENSPCANNO_GENE' 							=> "NA",
				'ENSPCANNO_GENE_LENGTH' 					=> ".",
				'ENSPCANNO_GENE_NUMBER_OF_TRANSCRIPTS' 		=> ".",
				'ENSPCANNO_GENE_DISTANCE_TO_START' 			=> ".",
				'ENSPCANNO_GENE_DISTANCE_TO_END' 			=> ".",
				'ENSPCANNO_TRANSCRIPT' 						=> "NA",
				'ENSPCANNO_TRANSCRIPT_BIOTYPE'				=> "NA",
				'ENSPCANNO_TRANSCRIPT_LENGTH' 				=> ".",
				'ENSPCANNO_TRANSCRIPT_START_SITE'			=> ".",
				'ENSPCANNO_TRANSCRIPT_DISTANCE_TO_START' 	=> ".",
				'ENSPCANNO_TRANSCRIPT_DISTANCE_TO_END' 		=> ".",
				'ENSPCANNO_TRANSCRIPT_CONSEQUENCE' 			=> "NA",
	);	
	return %data;
}

=head2 createVariantLongNonCodingHash

  Args       : None
  Example    : my %variantInformation = createVariantLongNonCodingHash()
  Description: Creates a hash that stores annotation information from long non-coding
               variants
  Returntype : Hash

=cut

sub createVariantLongNonCodingHash {
	my %data = ('ENSLNCANNO_GENE' 								=> "NA",
				'ENSLNCANNO_GENE_LENGTH' 						=> ".",
				'ENSLNCANNO_GENE_NUMBER_OF_TRANSCRIPTS' 		=> ".",
				'ENSLNCANNO_GENE_DISTANCE_TO_START' 			=> ".",
				'ENSLNCANNO_GENE_DISTANCE_TO_END' 				=> ".",
				'ENSLNCANNO_TRANSCRIPT' 						=> "NA",
				'ENSLNCANNO_TRANSCRIPT_BIOTYPE'					=> "NA",
				'ENSLNCANNO_TRANSCRIPT_LENGTH' 					=> ".",
				'ENSLNCANNO_TRANSCRIPT_DISTANCE_TO_START' 		=> ".",
				'ENSLNCANNO_TRANSCRIPT_DISTANCE_TO_END' 		=> ".",
				'ENSLNCANNO_TRANSCRIPT_CONSEQUENCE' 			=> "NA",
	);	
	return %data;
}

=head2 createVariantShortNonCodingHash

  Args       : None
  Example    : my %variantInformation = createVariantShortNonCodingHash()
  Description: Creates a hash that stores annotation information from short non-coding
               variants
  Returntype : Hash

=cut

sub createVariantShortNonCodingHash {
	my %data = ('ENSSNCANNO_GENE' 								=> "NA",
				'ENSSNCANNO_GENE_LENGTH' 						=> ".",
				'ENSSNCANNO_GENE_NUMBER_OF_TRANSCRIPTS' 		=> ".",
				'ENSSNCANNO_GENE_DISTANCE_TO_START' 			=> ".",
				'ENSSNCANNO_GENE_DISTANCE_TO_END' 				=> ".",
				'ENSSNCANNO_TRANSCRIPT' 						=> "NA",
				'ENSSNCANNO_TRANSCRIPT_BIOTYPE'					=> "NA",
				'ENSSNCANNO_TRANSCRIPT_LENGTH' 					=> ".",
				'ENSSNCANNO_TRANSCRIPT_DISTANCE_TO_START' 		=> ".",
				'ENSSNCANNO_TRANSCRIPT_DISTANCE_TO_END' 		=> ".",
				'ENSSNCANNO_TRANSCRIPT_CONSEQUENCE' 			=> "NA",
	);	
	return %data;
}

=head2 doesVariantOverlapWithAnyGene

  Args       : VariationFeature
  Example    : my ($answer, $overlapGene) = doesVariantOverlapWithProteinCodingGene($vf)
  Description: Tells whether a VariationFeature overlaps with a gene, and 
               the first of the overlapping genes is returned
  Returntype : Array with YES/NO and name of the overlapping gene

=cut

sub doesVariantOverlapWithAnyGene {
	my ($vf) = @_;
	my $answer = "NO";
	my $overlapGene = "NA";
	
	if (scalar @{$vf->get_overlapping_Genes()} > 0) {
		$answer = "YES";
		$overlapGene = @{$vf->get_overlapping_Genes()}[0]->display_xref->display_id();
	}
	
	return ($answer, $overlapGene);
}

=head2 doesVariantOverlapWithProteinCodingGene

  Args       : VariationFeature
  Example    : my ($answer, $overlapGene) = doesVariantOverlapWithProteinCodingGene($vf)
  Description: Tells whether a VariationFeature overlaps with a protein coding gene, and 
               the first of the overlapping genes is returned
  Returntype : Array with YES/NO and name of the overlapping gene

=cut

sub doesVariantOverlapWithProteinCodingGene {
	my ($vf) = @_;
	my $answer = "NO";
	my $overlapGene = "NA";
	foreach my $gene ( @{$vf->get_overlapping_Genes()} ) {
		if ($gene->biotype() eq "protein_coding") {
			$answer = "YES";
			if (scalar @{$vf->get_overlapping_Genes()} == 1) {
				$overlapGene = @{$vf->get_overlapping_Genes()}[0]->display_xref->display_id();	
			}
		}
	}
	return ($answer, $overlapGene);
}

=head2 doesVariantOverlapWithLongNonCodingGene

  Args       : VariationFeature
  Example    : my ($answer, $overlapGene) = doesVariantOverlapWithLongNonCodingGene($vf)
  Description: Tells whether a VariationFeature overlaps with a long non-coding gene, and 
               the first of the overlapping genes is returned
  Returntype : Array with YES/NO and name of the overlapping gene

=cut

sub doesVariantOverlapWithLongNonCodingGene {
	my ($vf) = @_;
	my $overlapGene = "NA";
	foreach my $gene ( @{$vf->get_overlapping_Genes()} ) {
		if (exists $longNoncodingTable{$gene->biotype()}) {
			if (scalar @{$vf->get_overlapping_Genes()} == 1) {
				$overlapGene = @{$vf->get_overlapping_Genes()}[0]->display_xref->display_id();	
			}
		}
	}
	return ($overlapGene);
}

=head2 doesVariantOverlapWithShortNonCodingGene

  Args       : VariationFeature
  Example    : my ($answer, $overlapGene) = doesVariantOverlapWithShortNonCodingGene($vf)
  Description: Tells whether a VariationFeature overlaps with a short non-coding gene, and 
               the first of the overlapping genes is returned
  Returntype : Array with YES/NO and name of the overlapping gene

=cut

sub doesVariantOverlapWithShortNonCodingGene {
	my ($vf) = @_;
	my $overlapGene = "NA";
	foreach my $gene ( @{$vf->get_overlapping_Genes()} ) {
		if (exists $shortNoncodingTable{$gene->biotype()}) {
			if (scalar @{$vf->get_overlapping_Genes()} == 1) {
				$overlapGene = @{$vf->get_overlapping_Genes()}[0]->display_xref->display_id();	
			}
		}
	}
	return ($overlapGene);
}

=head2 getNearestGene

  Args       : "VariationFeature", "Slice Adaptor", "Type", "Interval"
  Example    : my ($gene, $vartype) = getNearestGene($vf, $slice_adaptor, $type, $int)
  Description: Finds closest gene of the specified type any, protein coding, long non-coding 
               or short non-coding, within the specified number of bases up-/downstream 
               and tells whether it is upstream or downstream
  Returntype : Array with gene name and variant position

=cut

sub getNearestGene {
	my ($vf, $slice_adaptor, $type, $int) = @_;
	my $slice = $slice_adaptor->fetch_by_region('chromosome', $vf->slice->seq_region_name(), $vf->start() - $int, $vf->start() + $int);
	my %transcriptHash = ();
	my %upstreamDownstreamHash = ();
	
	foreach my $transcript ( @{$slice->get_all_Transcripts()} ) {
		if ($type eq "ANY") {
			my $distToGene = $vf->start() - $transcript->seq_region_start();				
			my $geneID = $transcript->get_Gene->display_xref->display_id();
			
			$transcriptHash{$geneID} = abs($distToGene);
			
			if ($transcript->strand() == 1) {
				if ($distToGene > 0) {
					$upstreamDownstreamHash{$geneID} = "upstream_gene_variant";
				}
				if ($distToGene < 0) {
					$upstreamDownstreamHash{$geneID} = "downstream_gene_variant";
				}
			}
			
			if ($transcript->strand() == -1) {
				if ($distToGene < 0) {
					$upstreamDownstreamHash{$geneID} = "upstream_gene_variant";
				}
				if ($distToGene > 0) {
					$upstreamDownstreamHash{$geneID} = "downstream_gene_variant";
				}
			}
		}
		if ($type eq "PC") {
			if ($transcript->biotype() eq "protein_coding") {
				my $distToGene = $vf->start() - $transcript->seq_region_start();				
				my $geneID = $transcript->get_Gene->display_xref->display_id();
				
				$transcriptHash{$geneID} = abs($distToGene);
				
				if ($transcript->strand() == 1) {
					if ($distToGene > 0) {
						$upstreamDownstreamHash{$geneID} = "upstream_gene_variant";
					}
					if ($distToGene < 0) {
						$upstreamDownstreamHash{$geneID} = "downstream_gene_variant";
					}
				}
				
				if ($transcript->strand() == -1) {
					if ($distToGene < 0) {
						$upstreamDownstreamHash{$geneID} = "upstream_gene_variant";
					}
					if ($distToGene > 0) {
						$upstreamDownstreamHash{$geneID} = "downstream_gene_variant";
					}
				}
			}
		}
		if ($type eq "LNC") {
			if (exists $longNoncodingTable{$transcript->biotype()}) {
				my $distToGene = $vf->start() - $transcript->seq_region_start();
				my $geneID = $transcript->get_Gene->display_xref->display_id();
				
				$transcriptHash{$geneID} = abs($distToGene);
				
				if ($transcript->strand() == 1) {
					if ($distToGene > 0) {
						$upstreamDownstreamHash{$geneID} = "upstream_gene_variant";
					}
					if ($distToGene < 0) {
						$upstreamDownstreamHash{$geneID} = "downstream_gene_variant";
					}
				}
				
				if ($transcript->strand == -1) {
					if ($distToGene < 0) {
						$upstreamDownstreamHash{$geneID} = "upstream_gene_variant";
					}
					if ($distToGene > 0) {
						$upstreamDownstreamHash{$geneID} = "downstream_gene_variant";
					}
				}
			}
		}
		if ($type eq "SNC") {
			if (exists $shortNoncodingTable{$transcript->biotype()}) {
				my $distToGene = $vf->start() - $transcript->seq_region_start();
				my $geneID = $transcript->get_Gene->display_xref->display_id();
				
				$transcriptHash{$geneID} = abs($distToGene);
				
				if ($transcript->strand() == 1) {
					if ($distToGene > 0) {
						$upstreamDownstreamHash{$geneID} = "upstream_gene_variant";
					}
					if ($distToGene < 0) {
						$upstreamDownstreamHash{$geneID} = "downstream_gene_variant";
					}
				}
				
				if ($transcript->strand() == -1) {
					if ($distToGene < 0) {
						$upstreamDownstreamHash{$geneID} = "upstream_gene_variant";
					}
					if ($distToGene > 0) {
						$upstreamDownstreamHash{$geneID} = "downstream_gene_variant";
					}
				}
			}
		}
	}
	
	my $closestGene = (sort {$transcriptHash{$a} <=> $transcriptHash{$b}} keys %transcriptHash)[0];
	
	if ($closestGene) {
		return ($closestGene, $upstreamDownstreamHash{$closestGene});
	}
	else {
		return ("NA","NA");
	}
}

=head2 findMostSeverelyHitTranscript

  Args       : "VariationFeature", "Transcript Adaptor", "Type"
  Example    : my ($msctranscript, $msvcons) = findMostSeverelyHitTranscript($vf, $ta, $type)
  Description: Finds the most severely hit transcript of the types "protein coding", 
               "long non-coding", "short non-coding" or "any". The transcripts are found in the
               prioritised order: 1. ENSEMBL most severe consequence 2. Is transcript canonical 
               3. Longest transcript            
  Returntype : Array with transcript stable id and most severe consequence

=cut

sub findMostSeverelyHitTranscript {
	my ($vf, $ta, $type) = @_;
	
	# Creates lists with transcripts and their consequences
	my @transcripts = ();
	my @consequences = ();
	foreach my $tv (@{$vf->get_all_TranscriptVariations()}) {
		foreach my $ct (@{$tv->consequence_type()}) {
			if ($type eq "PC") {
				if ($tv->transcript->biotype() eq "protein_coding") {
					push @transcripts, $tv->transcript->stable_id();
					push @consequences, $ct;
				}
			}
			if ($type eq "LNC") {
				if (exists $longNoncodingTable{$tv->transcript->biotype()}) {
					push @transcripts, $tv->transcript->stable_id();
					push @consequences, $ct;
				}
			}
			if ($type eq "SNC") {
				if (exists $shortNoncodingTable{$tv->transcript->biotype()}) {
					push @transcripts, $tv->transcript->stable_id();
					push @consequences, $ct;
				}
			}
			if ($type eq "ANY") {
					push @transcripts, $tv->transcript->stable_id();
					push @consequences, $ct;
			}			
		}
	}

	# Orders the transcripts by most severe consequence (MSC) using ENSEMBLs prioritised list
	my %table = ();
	for(my $i = 0; $i < scalar @consequences; $i++) {
		unless (exists $table{$transcripts[$i]} and $table{$transcripts[$i]} <= $consequecetable{$consequences[$i]}) {
			$table{$transcripts[$i]} = $consequecetable{$consequences[$i]};
		}
	}			
	my @MSCTranscripts = sort { $table{$a} <=> $table{$b} } keys(%table);
	
	# If multiple transcripts have the same MSC, the canonical transcript is chosen. If 
	# the canonical transcript is not among the transcripts, the longest transcript is 
	# returned.	
	my $rank = ();
	my $msctranscript = "NA";
	if (scalar @MSCTranscripts > 0) {
		my $rank = $table{$MSCTranscripts[0]};
		my %tablebylength = ();
		foreach my $MSCTranscript (@MSCTranscripts) {
			if (scalar @MSCTranscripts == 1) {
				$msctranscript = $MSCTranscript;
			}
			
			# Creates table with most severely hit transcripts and their lengths
			else {
				if ($table{$MSCTranscript} == $rank) {
					$tablebylength{$MSCTranscript} = $ta->fetch_by_stable_id($MSCTranscript)->length();
				}
			}
		}
		
		# Transcripts are sorted by length
		my @MSCTranscriptsSortedByLength = sort { $tablebylength{$b} <=> $tablebylength{$a} } keys(%tablebylength);
		
		# If there is a canonical transcript it is returned, else longest MSC transcript is
		if (scalar @MSCTranscripts > 1) {
			my $canonical = ();
			foreach my $lengthSorted (@MSCTranscriptsSortedByLength) {
				if ($ta->fetch_by_stable_id($lengthSorted)->is_canonical()) {
					$canonical = $lengthSorted;
				}
			}
			if ($canonical) {
				$msctranscript = $canonical;
			}
			if (!$canonical) {
				$msctranscript = $MSCTranscriptsSortedByLength[0];
			}
		}
	}
	
	# Finding consequence on the transcript
	my $msvcons = "NA";
	if ($msctranscript ne "NA") {
		my $msckey = $table{$msctranscript};
		($msvcons) = grep { $consequecetable{$_} eq $msckey } keys %consequecetable;
	}	
	return ($msctranscript, $msvcons);
}

=head2 findDistanceToStartAndEndOfGene

  Args       : "Gene", "VariationFeature"
  Example    : my ($distToStart, $distToEnd) = findDistanceToStartAndEndOfGene($gene, $vf)
  Description: Finds distance from variant to start and end of gene          
  Returntype : Array with distance from variant to start and end of gene

=cut

# Multiallelic variants will be positioned properly after comparison with the reference sequence
# see methods explanation (XXX) for details

sub findDistanceToStartAndEndOfGene {
	my ($gene, $vf) = @_;	
	my $distToStart = ".";
	my $distToEnd = ".";
		
	# Insertions and SNPs
	if ($vf->end() - $vf->start() < 1) {
		if ($gene->strand() > 0) {
			$distToStart = $vf->start() - $gene->seq_region_start();
			$distToEnd = $gene->seq_region_end() - $vf->start;
		}
		if ($gene->strand() < 0) {
			$distToStart = $gene->seq_region_end() - $vf->start();
			$distToEnd = $vf->start() - $gene->seq_region_start();
		}
	}
	
	# Deletions
	if ($vf->end() - $vf->start() > 0) {
		if ($gene->strand() > 0) {
			$distToStart = $vf->start() - $gene->seq_region_start();
			$distToEnd = $gene->seq_region_end() - $vf->end();
		}
		if ($gene->strand() < 0) {
			$distToStart = $gene->seq_region_end() - $vf->start();
			$distToEnd = $vf->end() - $gene->seq_region_start();
		}
	}
	return ($distToStart, $distToEnd);	
}

=head2 findDistanceToStartAndEndOfTranscript

  Args       : "Gene", "VariationFeature"
  Example    : my ($distToStart, $distToEnd) = findDistanceToStartAndEndOfTranscript($transcript, $vf) 
  Description: Finds distance from variant to start and end of transcript          
  Returntype : Array with distance from variant to start and end of transcript

=cut

sub findDistanceToStartAndEndOfTranscript {
	my ($transcript, $vf) = @_;	
	my $distToStart = ".";
	my $distToEnd = ".";
	
	# Insertions and SNPs
	if ($vf->end() - $vf->start() < 1) {
		if ($transcript->strand() > 0) {
			$distToStart = $vf->start() - $transcript->seq_region_start();
			$distToEnd = $transcript->seq_region_end() - $vf->start();
		}
		if ($transcript->strand() < 0) {
			$distToStart = $transcript->seq_region_end() - $vf->start();
			$distToEnd = $vf->start() - $transcript->seq_region_start();
		}
	}
	
	# Deletions
	if ($vf->end - $vf->start > 0) {
		if ($transcript->strand() > 0) {
			$distToStart = $vf->start() - $transcript->seq_region_start();
			$distToEnd = $transcript->seq_region_end() - $vf->end();
		}
		if ($transcript->strand() < 0) {
			$distToStart = $transcript->seq_region_end() - $vf->start();
			$distToEnd = $vf->end() - $transcript->seq_region_start();
		}
	}
	
	return ($distToStart, $distToEnd);
}

						# sub findDistanceToStartAndEndOfTranscript {
						# 	# Note: The length of the transcript is of the spliced transcript
						# 	
						# 	my ($transcript, $vf, $ta) = @_;
						# 	my $strand = $ta->fetch_by_stable_id($transcript)->strand;
						# 	my @exons = @{$ta->fetch_by_stable_id($transcript)->get_all_Exons()};
						# 	my $transciptLength = $ta->fetch_by_stable_id($transcript)->length;
						# 	my $distToStart = "NA";
						# 	my $distToEnd = "NA";
						# 	
						# 	# Insertions and SNPs	
						# 	if ($vf->end - $vf->start < 1) {
						# 		my $isInsideExon = "NO";
						# 		my $currentDist = 0;
						# 		my $totalDist = "NA";	
						# 		
						# 		# For transcripts on the 1 strand
						# 		if ($strand > 0) {
						# 			foreach my $exon (@exons) {
						# 				if ($vf->start >= $exon->start() && $vf->start <= $exon->end()) {
						# 					$totalDist = $currentDist + $vf->start - $exon->start();
						# 					$isInsideExon = "YES";
						# 				}
						# 				if ($totalDist eq "NA") {
						# 					$currentDist = $currentDist + $exon->end() - $exon->start() + 1;
						# 				}
						# 			}
						# 		}
						# 	
						# 		# For transcripts on the -1 strand
						# 		# Note: Exons start with high indexes and go towards lower
						# 		if ($strand < 0) {
						# 			foreach my $exon (@exons) {
						# 				if ($vf->start >= $exon->start() && $vf->start <= $exon->end()) {
						# 					$totalDist = $currentDist +  $exon->end() - $vf->start;
						# 					$isInsideExon = "YES";
						# 				}
						# 				if ($totalDist eq "NA") {
						# 					$currentDist = $currentDist + $exon->end() - $exon->start() + 1;
						# 				}
						# 			}
						# 		}
						# 		
						# 		# Calculating distances to start and end
						# 		if ($totalDist ne "NA") {
						# 			$distToStart = 	$totalDist;
						# 			$distToEnd = $transciptLength - $totalDist - 1;
						# 		}
						# 	}
						# 	
						# 	# Deletions
						# 	if ($vf->end - $vf->start > 0) {
						# 		my $lengthOfDeletion = $vf->end - $vf->start;
						# 		my $transcriptStart = $ta->fetch_by_stable_id($transcript)->start;
						# 		my $transcriptEnd = $ta->fetch_by_stable_id($transcript)->end;
						# 		
						# 		# Testing whether deletion start/end overlap with exons
						# 		my $isStartInsideExon = "NO";
						# 		my $isEndInsideExon = "NO";
						# 		foreach my $exon (@exons) {	
						# 			# Testing if start and end falls inside exons
						# 			if ($vf->start >= $exon->start() && $vf->start <= $exon->end()) {
						# 				$isStartInsideExon = "YES";
						# 			}
						# 			if ($vf->end >= $exon->start() && $vf->end <= $exon->end()) {
						# 				$isEndInsideExon = "YES";
						# 			}
						# 		}
						# 		
						# 		# For transcripts on the 1 strand
						# 		my $beforeDeletionDist = 0;
						# 		my $afterDeletionDist = 0;
						# 		my $tempDist = 0;
						# 		my $tempDist2 = 0;
						# 		my $afterDeletionStarter = "NO";
						# 		my $isExonDeleted = "NO";
						# 		
						# 		if ($strand > 0) {
						# 			foreach my $exon (@exons) {	
						# 		
						# 				# Deletions where start and end of deletion fall outside exons
						# 				if ($isStartInsideExon eq "NO" && $isEndInsideExon eq "NO") {
						# 				
						# 					# Testing whether an exon is found between the deletion start and end
						# 					if ($vf->start <= $exon->start() && $vf->end >= $exon->end()) {
						# 						$isExonDeleted = "YES";
						# 					}
						# 				
						# 					# Getting the length of all exons that come before the deletion-start
						# 					if ($exon->end() < $vf->start) {
						# 						$beforeDeletionDist = $exon->end() - $exon->start() + 1 + $beforeDeletionDist;
						# 					}
						# 				
						# 					# Getting the length of all exons that come after the deletion-end
						# 					if ($exon->start() > $vf->end) {
						# 						$afterDeletionDist = $exon->end() - $exon->start() + 1 + $afterDeletionDist;
						# 					}
						# 				
						# 					# Result is only returned if an exon is deleted
						# 					if ($isExonDeleted eq "YES") {
						# 						$distToStart = $beforeDeletionDist;
						# 						$distToEnd = $afterDeletionDist;
						# 					}
						# 				}
						# 			
						# 				# Deletions where both start and end falls inside exons
						# 				if ($isStartInsideExon eq "YES" && $isEndInsideExon eq "YES") {
						# 				
						# 					# Getting distance up to start of deletion
						# 					if ($vf->start >= $exon->start() && $vf->start <= $exon->end()) {
						# 						$beforeDeletionDist = $tempDist + $vf->start - $exon->start();
						# 					}
						# 					if ($beforeDeletionDist == 0) {
						# 						$tempDist = $tempDist + $exon->length();
						# 					}
						# 				
						# 					# Getting distance after end of deletion
						# 					if ($afterDeletionStarter eq "YES") {
						# 						$afterDeletionDist = $afterDeletionDist + $exon->end() - $exon->start() + 1;
						# 					}
						# 					if ($vf->end >= $exon->start() && $vf->end <= $exon->end()) {
						# 						$afterDeletionStarter = "YES";
						# 						$afterDeletionDist = $exon->end() - $vf->end;
						# 					}					
						# 					$distToStart = $beforeDeletionDist;
						# 					$distToEnd = $afterDeletionDist;
						# 				}
						# 			
						# 				# Testing all scenarios where start falls outside exons and the end inside
						# 				if ($isStartInsideExon eq "NO" && $isEndInsideExon eq "YES") {
						# 				
						# 					# Getting length of all exons before start of deletion
						# 					if ($vf->start > $exon->end()) {
						# 						$beforeDeletionDist = $beforeDeletionDist + $exon->end() - $exon->start() + 1;
						# 					}
						# 					
						# 					# Getting distance after end of deletion
						# 					if ($afterDeletionStarter eq "YES") {
						# 						$afterDeletionDist = $afterDeletionDist + $exon->end() - $exon->start() + 1;
						# 					}
						# 					if ($vf->start >= $exon->start() && $vf->start <= $exon->end()) {
						# 						$afterDeletionStarter = "YES";
						# 						$afterDeletionDist = $exon->end() - $vf->start;
						# 					}					
						# 				
						# 					$distToStart = $beforeDeletionDist;
						# 					$distToEnd = $afterDeletionDist;
						# 				}
						# 			
						# 				# Testing all scenarios where start falls inside exons and the end outside
						# 				if ($isStartInsideExon eq "YES" && $isEndInsideExon eq "NO") {	
						# 
						# 					# Getting distance up to start of deletion
						# 					if ($vf->start >= $exon->start() && $vf->start <= $exon->end()) {
						# 						$beforeDeletionDist = $tempDist + $vf->start - $exon->start();
						# 					}
						# 					if ($beforeDeletionDist == 0) {
						# 						$tempDist = $tempDist + $exon->end() - $exon->start() + 1;
						# 					}
						# 				
						# 					# Getting length of exons after the end of deletion
						# 					if ($vf->end < $exon->start()) {
						# 						$afterDeletionDist = $afterDeletionDist + $exon->end() - $exon->start() + 1;
						# 					}
						# 					$distToStart = $beforeDeletionDist;
						# 					$distToEnd = $afterDeletionDist;
						# 				}
						# 			}
						# 		}
						# 		
						# 		# For transcripts on the -1 strand
						# 		# Note: Exons start with high indexes and go towards lower
						# 		if ($strand < 0) {
						# 			foreach my $exon (@exons) {	
						# 		
						# 				# Deletions where start and end of deletion fall outside exons
						# 				if ($isStartInsideExon eq "NO" && $isEndInsideExon eq "NO") {
						# 				
						# 					# Testing whether an exon is found between the deletion start and end
						# 					if ($vf->start <= $exon->start() && $vf->end >= $exon->end()) {
						# 						$isExonDeleted = "YES";
						# 					}
						# 				
						# 					# Getting the length of all exons that come before the deletion-start
						# 					if ($exon->start() > $vf->end) {
						# 						$beforeDeletionDist = $exon->end() - $exon->start() + 1 + $beforeDeletionDist;
						# 					}
						# 				
						# 					# Getting the length of all exons that come after the deletion-end
						# 					if ($exon->end() < $vf->start) {
						# 						$afterDeletionDist = $exon->end() - $exon->start() + 1 + $afterDeletionDist;
						# 					}
						# 				
						# 					# Result is only returned if an exon is deleted
						# 					if ($isExonDeleted eq "YES") {
						# 						$distToStart = $beforeDeletionDist;
						# 						$distToEnd = $afterDeletionDist;
						# 					}
						# 				}
						# 			
						# 				# Deletions where both start and end falls inside exons
						# 				if ($isStartInsideExon eq "YES" && $isEndInsideExon eq "YES") {
						# 				
						# 					# Getting distance up to start of deletion
						# 					if ($vf->end >= $exon->start() && $vf->end <= $exon->end()) {
						# 						$beforeDeletionDist = $tempDist + $exon->end() - $vf->end;
						# 					}
						# 					if ($beforeDeletionDist == 0) {
						# 						$tempDist = $tempDist + $exon->end() - $exon->start() + 1;
						# 					}
						# 				
						# 					# Getting distance after end of deletion
						# 					if ($afterDeletionStarter eq "YES") {
						# 						$afterDeletionDist = $afterDeletionDist + $exon->end() - $exon->start() + 1;
						# 					}
						# 					if ($vf->start >= $exon->start() && $vf->start <= $exon->end()) {
						# 						$afterDeletionStarter = "YES";
						# 						$afterDeletionDist = $vf->start - $exon->start();
						# 					}					
						# 					$distToStart = $beforeDeletionDist;
						# 					$distToEnd = $afterDeletionDist;
						# 				}
						# 			
						#  				# Testing all scenarios where start falls outside exons and the end inside
						# 				if ($isStartInsideExon eq "NO" && $isEndInsideExon eq "YES") {
						# 				
						# 					# Getting distance from trans-start to deletion
						# 					if ($vf->end >= $exon->start() && $vf->end <= $exon->end()) {
						# 						$beforeDeletionDist = $tempDist + $exon->end() - $vf->end;
						# 					}
						# 					if ($beforeDeletionDist == 0) {
						# 						$tempDist = $tempDist + $exon->end() - $exon->start() + 1;
						# 					}
						# 					
						# 					# Getting distance from trans-end to deletion
						# 					if ($vf->start > $exon->end()) {
						# 						$afterDeletionDist = $afterDeletionDist + $exon->length;
						# 					}
						# 				}
						# 
						# 				# Testing all scenarios where start falls inside exons and the end outside
						# 				if ($isStartInsideExon eq "YES" && $isEndInsideExon eq "NO") {	
						# 				
						# 					# Getting distance from trans start to deletion 					
						# 					if ($vf->end < $exon->start()) {
						# 						$beforeDeletionDist = $beforeDeletionDist + $exon->end() - $exon->start() + 1;
						# 					}
						# 					
						# 					# Getting distance from trans end to deletion
						# 					if ($afterDeletionStarter eq "YES") {
						# 						$afterDeletionDist = $afterDeletionDist + $exon->end() - $vf->start + 1;
						# 					}
						# 					if ($vf->start >= $exon->start() && $vf->start <= $exon->end()) {
						# 						$afterDeletionStarter = "YES";
						# 						$afterDeletionDist = $exon->end() - $vf->start;
						# 					}
						# 					$distToStart = $beforeDeletionDist;
						# 					$distToEnd = $afterDeletionDist;
						# 				}
						# 			}
						# 		}
						# 	}
						# 	return ($distToStart, $distToEnd);
						# }
return 1;