#!/usr/bin/perl
use strict;
use warnings;

my $inputfile = shift @ARGV or die "No inputfile Detector-* file provided\n";
open(my $fh, "<$inputfile") or die "Cannot open $inputfile: $!\n";
my $n = 0;
my $m = 0;
my @lines;
while(my $line = <$fh>) {
  chomp $line;
  if($line =~ /Element:/) {
    #print $line,"\n";
    my @cols = split(/\s+/,$line);
    my $ADC = $cols[1];
    my $TDC = $cols[6];
    my $gain = $cols[5];
    $ADC =~ s/M\d$//;
    $TDC =~ s/M\d$//;
    my @xyz = @cols[11..13];
    printf("%3d %f \n",
           $n, $gain);
    $n++;
  }
}

close $fh;
