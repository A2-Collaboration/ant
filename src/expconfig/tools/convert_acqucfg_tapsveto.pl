#!/usr/bin/perl
use strict;
use warnings;
use Data::Dumper;


my $inputfile   = $ARGV[0];
open(my $fh, "<$inputfile") or die "Cannot open $inputfile: $!\n";

sub is_pbwo4 {
  my $n = shift;
  return $n % 64 < 3;
}

sub extract_tdc_adc {
  my @cols = @_;
  my $ADC = $cols[1];
  my $TDC = $cols[6];
  $ADC =~ s/M\d$//;
  $TDC =~ s/M\d$//;
  return [$TDC, $ADC]; # put TDC first, since the rest is LG, LGS and so on
}

my $n=0;

while (my $line = <$fh>) {
  chomp $line;
  if ($line =~ /Element:/) {
    #print $line,"\n";
    my @cols = split(/\s+/,$line);
    my $ADC = $cols[1];
    my $TDC = $cols[6];
    $ADC =~ s/M\d$//;
    $TDC =~ s/M\d$//;
    my @xy = @cols[11..12];
    if(is_pbwo4($n)) {
      printf("{%3d, {%7.3f, %7.3f}, %04d, %04d, %04d}, \n",
             $n, @xy, $TDC, $ADC, $ADC+100);
    }
    # if(is_pbwo4($n)) {
    #   printf("{%3d, {%7.3f, %7.3f}, %04d, %04d, %04d}, \n",
    #          $n, @xy, $TDC, $ADC, $ADC+1);
    # }
    $n++;
  }
}

close $fh;
