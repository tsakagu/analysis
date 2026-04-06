#!/usr/bin/perl

use strict;
use DBI; # Database interface 

# Read file
my $filename = 'RunList_ana450_2024p009_DSTs_available.txt';
open(my $fh,'<',$filename) or die "Could not open file '$filename' $!";

# connect to the database daq
my $dbh = DBI->connect("dbi:ODBC:FileCatalog","") || die $DBI::error;

while (my $line = <$fh>) {
    chomp($line);
    my $gettab = $dbh->prepare("SELECT runnumber FROM datasets WHERE runnumber=$line AND dataset LIKE 'ana450_2024p009' AND dsttype LIKE 'DST_CALO_run2pp'GROUP BY runnumber HAVING SUM(events) >= 500000;");
    $gettab->execute();
# print the first entries:
    my @row;
    while (@row = $gettab->fetchrow_array) {
        print join(" ", @row), "\n";
    }
}

close($fh); # Close the file handle.

$dbh->disconnect;
