#!/usr/bin/perl

use strict;
use DBI; # Database interface 

# Read file
my $filename = 'RunList_golden.txt';
open(my $fh,'<',$filename) or die "Could not open file '$filename' $!";

# connect to the database Production_read
my $dbh = DBI->connect("dbi:ODBC:RawdataCatalog_read","") || die $DBI::error;

my @hostnames = ("seb00","seb01","seb02","seb03","seb04","seb05","seb06","seb07","seb08",
                 "seb09","seb10","seb11","seb12","seb13","seb14","seb15"); # Only check SEBs for the EMCal

while (my $line = <$fh>) {
    chomp($line);
    my $feminfoname = "FEM_folder/fem_${line}.txt";
    open(my $feminfofile,'<',$feminfoname) or die "Could not open '$feminfoname' $!";

    # Determine whether a run is good or not wrt FEM Clock matching
    my $fem_matching = 1;
    
    # Common clock value (one per file)
    my $clock_value = -1;

    # Current scanned Sub Event Buffer
    my $current_seb = "sebXX";
    
    # equals to 1 if the current SEB is considered
    # equals to 0 if the current SEB is ignored
    my $checkClock = 0;

    while (my $rawline = <$feminfofile>) {
        if ($rawline =~ /^(seb\S+)/) {
            $current_seb = $1;
            # If seb is considered, check FEM Clock matching
            if ($current_seb ~~ @hostnames) {
                # print "$current_seb\n";
                $checkClock = 1;
            }
            # If SEB is not considered, ignore it until the next SEB or the end of the file;
            else {
                $checkClock = 0;
            }
        }

        # Check Clock matching
        if ($checkClock) {
            # Check whether all lines agree between them
            if ($rawline =~ /^FEM Clock: \s+(\S+)./) {
                # the first entry defines the clock value
                if ($clock_value == -1) {
                    if (($current_seb eq "seb16") || ($current_seb eq "seb17") || ($current_seb eq "seb20")) {
                        $clock_value = $1 - 1;
                    }
                    else {
                        $clock_value = $1;
                    }
                }
                else {
                    if ((($current_seb eq "seb16") || ($current_seb eq "seb17") || ($current_seb eq "seb20")) && ($1 ne $clock_value + 1)) {
                        $fem_matching = 0;
                        last;
                    }
                    elsif ((($current_seb ne "seb16") && ($current_seb ne "seb17") || ($current_seb eq "seb20")) && ($1 ne $clock_value)) {
                        $fem_matching = 0;
                        last;
                    }
                }
            }

            # Check whether all clock numbers on the same line agree.
            if ($rawline =~ /^FEM Clock:\s+(\S+)\s+(\S+)\s+(\S+)/) {
                if (($1 ne $2) || ($1 ne $3) || ($2 ne $3)) {
                    $fem_matching = 0;
                    last;
                }
            } 
            elsif ($rawline =~ /^FEM Clock: \s+(\S+)\s+(\S+)/) {
                if ($1 ne $2) {
                    $fem_matching = 0;
                    last;
                }
            }
        }
    }
    close($feminfofile);
    
    if ($fem_matching) {
        print "$line\n";
    }
}

close($fh); # Close the file handle.

$dbh->disconnect;
