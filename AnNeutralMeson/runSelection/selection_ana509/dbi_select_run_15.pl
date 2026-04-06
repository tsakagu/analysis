#!/usr/bin/env perl

use strict;
use DBI; # Database interface 

# Read file
my $filename = 'RunList_spinqa.txt';
open(my $fh,'<',$filename) or die "Could not open file '$filename' $!";

# connect to the database Production_read
my $dbh = DBI->connect("dbi:ODBC:SpinDB","") || die $DBI::error;

while (my $line = <$fh>) {
    chomp($line);
    my $gettab = $dbh->prepare("SELECT runnumber FROM spin WHERE runnumber=$line AND is_default=true AND crossingangle BETWEEN -1.7 AND -1.3;");
    $gettab->execute();
# print the first entries:
    my @row;
    while (@row = $gettab->fetchrow_array) {
        print join(" ", @row), "\n";
    }
}

close($fh); # Close the file handle.

$dbh->disconnect;
