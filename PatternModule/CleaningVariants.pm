package PatternModule::CleaningVariants;

use strict;
use warnings;
use List::MoreUtils qw(firstidx);

=head2 cleanSNP

  Args       : "Reference Sequence" and "Alternative Sequence"
  Example    : my ($ref, $alt, $snpIndex) = cleanSNP($ref, $alt)
  Description: If a SNP is represented as a complex variant, this method simplifies the 
               representation to the simplest possible sequences, and also returns the 
               positional difference compared to the original start position.
  Returntype : Array with three elements

=cut

sub cleanSNP {
	my ($ref, $alt) = @_;
	my @varIndex = ();
	
	# Ensuring that it is a SNP
	if (length($ref) == length($alt) && length($ref) != 1) {
		
		# Comparing the sequences of the ref and alt
		my @refSplit = split('',$ref);
		my @altSplit = split('',$alt);
		
		# Finding the index of the difference
		for(my $i = 0; $i < scalar @refSplit; $i++) {
			if ($refSplit[$i] eq $altSplit[$i]) {
				push @varIndex, "S";
			}
			if ($refSplit[$i] ne $altSplit[$i]) {
				push @varIndex, "D";
			}
		}
 		
		# Counting number of differences
		my $count = grep (/D/, @varIndex);
		
		# If there is only one difference between ref and alt, the sequences are reduced 
		# to the simplest SNP, and the index of the difference is found 
		my $snpIndex = 0;
		if ($count == 1) {
			$snpIndex = firstidx { $_ eq 'D' } @varIndex;
			$ref = $refSplit[$snpIndex];
			$alt = $altSplit[$snpIndex];
		}	
		return ($ref, $alt, $snpIndex);
	}
}

=head2 cleanDELETIONS

  Args       : "Reference Sequence" and "Alternative Sequence"
  Example    : my ($ref, $alt, $delIndex) = cleanDELETIONS($ref, $alt)
  Description: This method simplifies the reference and alternative sequences to the 
               simplest possible sequences, and also returns the positional difference 
               compared to the original start position.
  Returntype : Array with three elements

=cut
	
sub cleanDELETIONS {
	my ($ref, $alt) = @_;
	my @varIndexForward = ();
	
	# Ensuring that we catch only deletions
	if (length($ref) > length($alt) && length($alt) > 1) {
		
		# Comparing the sequences of the ref and alt
		my @refSplit = split('',$ref);
		my @altSplit = split('',$alt);
		
		# This section finds the first difference in the string beginning from the start
		# and sliding the comparison of the array forward
		for(my $i = 0; $i < scalar @altSplit; $i++) {
			if ($refSplit[$i] eq $altSplit[$i]) {
				push @varIndexForward, "S";
			}
			if ($refSplit[$i] ne $altSplit[$i]) {
				push @varIndexForward, "D";
			}
		}
				
		# Finding number of differences from start to end
		my $count = grep (/D/, @varIndexForward);
		
		# If simple deletion no differences are found
		if ($count == 0) {
			# Determining how much of the sequences can be simplified
			my $simplifyable = length($alt) - 1;
		
			# Simplifying
			$ref = join("", @refSplit[$simplifyable..(length($ref)-1)]);
			$alt = $altSplit[$simplifyable];

			return ($ref, $alt, $simplifyable);
		}
		
		# If complex deletion differences are found. Indexes of the start of the difference
		# and the end of the difference is found, and this is where the deletion is. All
		# excess information is stripped.
		if ($count > 0) {
			my $ForwardIndex = firstidx { $_ eq 'D' } @varIndexForward;
			
			# Simplifying
			$ref = join("", @refSplit[$ForwardIndex - 1 ..(length($ref)-length($alt)+$ForwardIndex - 1)]);
			$alt = $refSplit[$ForwardIndex - 1];
			
			return ($ref, $alt, ($ForwardIndex - 1));
		}
	}
} 

=head2 cleanINSERTIONS

  Args       : "Reference Sequence" and "Alternative Sequence"
  Example    : my ($ref, $alt, $delIndex) = cleanINSERTIONS($ref, $alt)
  Description: This method simplifies the reference and alternative sequences to the 
               simplest possible sequences, and also returns the positional difference 
               compared to the original start position.
  Returntype : Array with three elements

=cut

sub cleanINSERTIONS {
	my ($ref, $alt) = @_;
	my @varIndexForward = ();
	
	# Ensuring that we catch only insertions	
	if (length($ref) < length($alt) && length($ref) > 1) {
		
		# Splitting ref and alt sequences
		my @refSplit = split('',$ref);
		my @altSplit = split('',$alt);
		
		# Looking for differences between reference and insertion sequences from start to end
		for(my $i = 0; $i < scalar @refSplit; $i++) {
			if ($refSplit[$i] eq $altSplit[$i]) {
				push @varIndexForward, "S";
			}
			if ($refSplit[$i] ne $altSplit[$i]) {
				push @varIndexForward, "D";
			}
		}
		
		# Finding number of differences from start to end
		my $count = grep (/D/, @varIndexForward);
		
		# If simple insertion no differences are found
		if ($count == 0) {
			# Determining how much of the sequences can be simplified
			my $simplifyable = length($ref) - 1;
		
			# Simplifying
			$ref = $refSplit[$simplifyable];
			$alt = join("", @altSplit[$simplifyable..(length($alt)-1)]);

			return ($ref, $alt, $simplifyable);
		}
		
		# If complex insertion (Reference sequence does not match start of alternative sequence)
		# Position of insertion start is found, and everything before and after insertion is stripped
		if ($count > 0) {
			my $ForwardIndex = firstidx { $_ eq 'D' } @varIndexForward;
			
			# Simplifying
			$alt = join("", @altSplit[$ForwardIndex - 1 ..(length($alt)-length($ref)+$ForwardIndex - 1)]);
			$ref = $refSplit[$ForwardIndex - 1];
			
			return ($ref, $alt, ($ForwardIndex - 1));
		}	
	}
}	

return 1;