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
    $ADC =~ s/M\d$//;
    $TDC =~ s/M\d$//;
    my @xyz = @cols[11..13];
    push(@lines,
         sprintf("{%3d, {%10.6f, %10.6f, %10.6f}, %04d, %04d, {%%s}}, \n",
                 $n, @xyz, $ADC, $TDC));
    $n++;
  }
  elsif($line =~ /Next-Neighbour:/) {
    #print $line, "\n";
    my @cols = split(/\s+/,$line);
    my @neighbours = @cols[3..$#cols];
    printf($lines[$m], join(",", @neighbours));
    $m++;
  }
}

close $fh;
