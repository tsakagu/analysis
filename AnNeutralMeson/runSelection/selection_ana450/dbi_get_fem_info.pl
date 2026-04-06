#!/usr/bin/perl

use strict;
use DBI; # Database interface 

# Read file
my $filename = 'RunList_golden.txt';
open(my $fh,'<',$filename) or die "Could not open file '$filename' $!";

# connect to the database Production_read
my $dbh = DBI->connect("dbi:ODBC:RawdataCatalog_read","") || die $DBI::error;

my @hostnames = ("seb00","seb01","seb02","seb03","seb04","seb05","seb06","seb07","seb08",
                 "seb09","seb10","seb11","seb12","seb13","seb14","seb15","seb16","seb17", "seb18", "seb19", "seb20"); # All SEBS (EMCal, IHCAL, OHCAL, MBD)

while (my $line = <$fh>) {
    chomp($line);

    my $file_path = "../FEM_folder/fem_${line}.txt";
    if -f $file_path {
        next;
    }

    my $fem_matching = 2;
    # select the number of events for each hostname in the given runnumber

    #print "line = $line\n";
    my $gettab = $dbh->prepare("SELECT filename,hostname FROM filelist WHERE runnumber=$line AND filename LIKE '%0000.prdf';");
    $gettab->execute();

    # Fetch hostname, filename and store them in a hash
    my %values;
    while (my ($rawfilename, $hostname) = $gettab->fetchrow_array) {
        $values{$hostname} = $rawfilename;
        #print "values($hostname) = $rawfilename\n";
    }
    
    my $clock_value = -1;
    open(outfh, ">", "fem_$line.txt");
    foreach my $host(@hostnames) {
        unless (exists($values{$host})) {
            $fem_matching = 0;
        }

        my $rawfilename = $values{$host};

        # By default, the bufferbox content is manually transferred to /sphenix/lustre01/sphnxpro.
        # RawdataCatalog does not account for that
        $rawfilename =~ s{^/bbox/bbox\d+/}{/sphenix/lustre01/sphnxpro/physics/};
        
        my $output = `ddump -i -e 10 $rawfilename | grep "FEM Clock"`;
        
        print "rawfilename = $rawfilename\n";
        print outfh "$host\n";
        print outfh "$output\n";
    }
    close(outfh);
}

close($fh); # Close the file handle.

$dbh->disconnect;
