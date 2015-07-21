#!/usr/bin/perl
use strict;
use warnings;

my $inputfile = shift @ARGV or die "No inputfile Detector-PID-* file provided\n";
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
    my @xyz = @cols[11..13];
    printf("{%2d, %8.3f, %03d, %04d}, \n",
                 $n, $xyz[2], $ADC, $TDC);
    $n++;
  }
}

close $fh;
