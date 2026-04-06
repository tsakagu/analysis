#!/usr/bin/perl

use strict;
use DBI; # Database interface 

# connect to the database FileCatalog
my $dbh = DBI->connect("dbi:ODBC:FileCatalog","") || die $DBI::error;

my $gettab = $dbh->prepare("SELECT DISTINCT runnumber FROM datasets WHERE dataset LIKE 'ana450_2024p009%' AND dsttype LIKE 'DST_CALO_run2pp' ORDER BY runnumber;");
$gettab->execute();
# print the first entries:
my @row;
while (@row = $gettab->fetchrow_array) {
    print join(", ", @row), "\n";
}

$dbh->disconnect;
