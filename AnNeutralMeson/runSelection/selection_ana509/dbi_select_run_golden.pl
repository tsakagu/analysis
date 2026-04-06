#!/usr/bin/env perl

use strict;
use DBI; # Database interface 

# Read file
my $filename = 'RunList_500k.txt';
open(my $fh,'<',$filename) or die "Could not open file '$filename' $!";

# connect to the database Production_read
my $dbh = DBI->connect("dbi:ODBC:Production_read","") || die $DBI::error;

while (my $line = <$fh>) {
    chomp($line);
    my $gettab = $dbh->prepare("SELECT runnumber FROM goodruns WHERE runnumber=$line AND CAST(emcal_auto AS TEXT) LIKE '%GOLDEN%';");
    $gettab->execute();
# print the first entries:
    my @row;
    while (@row = $gettab->fetchrow_array) {
        print join(" ", @row), "\n";
    }
}

close($fh); # Close the file handle.

$dbh->disconnect;
