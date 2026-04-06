#!/usr/bin/env perl

use strict;
use DBI; # Database interface 

# connect to the database FileCatalog
my $dbh = DBI->connect("dbi:ODBC:FileCatalog","") || die $DBI::error;

my $gettab = $dbh->prepare("SELECT runnumber FROM datasets WHERE dsttype ='DST_CALOFITTING' AND dataset = 'run2pp' AND tag='ana509_2024p022_v001' GROUP BY runnumber HAVING SUM(events) >= 500000 ORDER BY runnumber");
$gettab->execute();
my @row;
while (@row = $gettab->fetchrow_array) {
    print join(", ", @row), "\n";
}

$dbh->disconnect;
