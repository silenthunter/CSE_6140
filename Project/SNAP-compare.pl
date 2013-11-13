use strict;

open(SNAP, $ARGV[0]);
open(CUDA, $ARGV[1]);

my %snapBC;
my %cudaBC;

my $e = .07;

while(<SNAP>)
{
	my $line = $_;
	chomp $line;
	#print $line."\n";

	if($line =~ m/^BC\s+(\d+)\s+([\d\.]+)/)
	{
		#print $1.":".$2."\n";
		$snapBC{$1} = $2;
	}
}


my $cnt = 0;
while(<CUDA>)
{
	my $line = $_;
	chomp $line;

	if($line =~ m/^\d/)
	{
		$cnt++;
		$cudaBC{$cnt} = $line;
	}
}

#Check to make sure the numbers are the same
foreach my $key (keys(%snapBC))
{
	my $snapNum = $snapBC{$key};
	my $cudaNum = $cudaBC{$key};
	
	my $result = $snapNum - $cudaNum;
	if($result > $snapNum * $e)
	{
		print $key." = ". $snapBC{$key}." : ".$cudaBC{$key}."\n";
	}
	#print $key." = ". $snapBC{$key}." : ".$cudaBC{$key}."($result)\n";
}
print "Done!\n";
