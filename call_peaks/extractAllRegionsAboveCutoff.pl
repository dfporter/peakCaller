use strict;

my $pileupFile = shift;
my $strand = shift;

my $c = 9; #cutoff
my %regions;

my $lastLoci = 0;
my $lastChr = '';

my $intervalChr = '';
my $intervalStart = 0;
my $intervalEnd = 0;
my $minValueInInterval = 0;
my $maxValueInInterval = 0;
my $positionOfmaxValueInInterval = 0;
my $in_interval = 0;

my $throughFirstLoop = 0;
my $debugl = "";
open(FH, "<$pileupFile") or die($!);
while (<FH>) {
	chomp $_;

	my @s = split(/\t/, $_);
	#was the last loci in an inteval?
	if($in_interval) {
		#are we on a new chromosome?
		if( $s[0] ne $intervalChr ) {
			#did we just enter this loop?
			if($throughFirstLoop) {
				#print last interval
				#my $l = "$intervalChr\t$intervalStart\t$intervalEnd\tmin=$minValueInInterval\tmax=$maxValueInInterval\n";
				my $l = "$intervalChr\t$intervalStart\t$intervalEnd\t.\t$maxValueInInterval\t$strand\n";
				print $l;
				$debugl = $debugl . "hit new chromosome.";
				#new interval
				if($s[3] > $c) {
					$debugl = $_;
				#new interval
					$intervalChr = $s[0];
					$intervalStart = $s[1];
					$intervalEnd = $s[1];
					$in_interval = 1;
				}
				else {
					$in_interval = 0;
				}
			}
			else { $throughFirstLoop = 1; } #we've been through while loop once, it's ok to output now
		}
		else { #same chr
			if( $s[1] < int($lastLoci) ){
				print "error! file not sorted\n";
				print $s[0] . $s[1] . "\n";
				print "the chromosome we are looking for is -" . $intervalChr . "-\n";
				print $debugl;
				exit;
			}
			if( $s[1] > int($lastLoci + 1) ) { #did we jump ahead?
				#print last interval
				my $l = "$intervalChr\t$intervalStart\t$intervalEnd\t.\t$maxValueInInterval\t$strand\n";
				print $l;
				
				$debugl = $debugl . "same chromosome, but more than one nt advance.";
				#new interval
				if($s[3] > $c) {
					$debugl = $_;
					$intervalChr = $s[0];
					$intervalStart = $s[1];
					$intervalEnd = $s[1];
					$in_interval = 1;
					$maxValueInInterval = $s[3];
					$minValueInInterval = $s[3];
				}
				else {
					$in_interval = 0;
				}
			}
			else { #same chr, one step ahead
				$debugl = $debugl . "same chr, one step ahead";
				if ($s[3] > $c) { #are we above the depth cutoff?
					#extend the interval
					$debugl = $debugl . "extending interval:" . $_;
					$intervalEnd = $s[1];
					if($s[3] > $maxValueInInterval) { $maxValueInInterval = $s[3]; }
					if($s[3] < $minValueInInterval) { $minValueInInterval = $s[3]; }
				}
				else {	
					$debugl =  $debugl . "dropped below cutoff";
					#print last interval
					my $l = "$intervalChr\t$intervalStart\t$intervalEnd\t.\t$maxValueInInterval\t$strand\n";
					print $l;
					#new interval
					if($s[3] > $c) {
						$debugl = $_;
						$intervalChr = $s[0];
						$intervalStart = $s[1];
						$intervalEnd = $s[1];
						$in_interval = 1;
						$maxValueInInterval = $s[3];
	                                	$minValueInInterval = $s[3];
					}
					else {
						$in_interval = 0;
					}
				}
			} #block of same chr, one step ahead
		} #block of same chr	
	}#block for last loci being in an interval
	else { #last loci not in an interval
		
		#new interval if above depth
		if($s[3] > $c) {
			$debugl = $_;
			$intervalChr = $s[0];
			$intervalStart = $s[1];
			$intervalEnd = $s[1];
			$in_interval = 1;
			$maxValueInInterval = $s[3];
                        $minValueInInterval = $s[3];
		}
		else {
			$in_interval = 0;
		}
	}

	$lastLoci = $s[1];
	$lastChr = $s[0];
}

close FH;
