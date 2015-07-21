#!/usr/bin/perl
use strict;
use warnings;

my $inputfile = shift @ARGV or die "No inputfile Detector-* file provided\n";
open(my $fh, "<$inputfile") or die "Cannot open $inputfile: $!\n";
my $n = 0;
while(my $line = <$fh>) {
  chomp $line;
  if($line =~ /Element:/) {
    #print $line,"\n";
    my @cols = split(/\s+/,$line);
    my $ADC = $cols[1];
    my $TDC = $cols[6];
    $ADC =~ s/M\d$//;
    $TDC =~ s/M\d$//;
    #my @xyz = @cols[11..13];
    my $electronenergy = $cols[14];
    my $scaler = $cols[16];
    printf("{%3d, %4d, %4d, %4d, %7.3f}, \n",
            $n, $TDC, $scaler, $ADC, $electronenergy);
    $n++;
  }
}

close $fh;
