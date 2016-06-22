#!/usr/bin/perl
use strict;
use warnings;
use Data::Dumper;

my $inputfile   = $ARGV[0];
open(my $fh, "<$inputfile") or die "Cannot open $inputfile: $!\n";

sub extract_tdc_adc {
  my @cols = @_;
  my $ADC = $cols[1];
  my $TDC = $cols[6];
  $ADC =~ s/M\d$//;
  $TDC =~ s/M\d$//;
  return [$TDC, $ADC]; # put TDC first, since the rest is LG, LGS and so on
}

my $data;
my $n=0;
my $m=0;
my $l=0;

while (my $line = <$fh>) {
  chomp $line;
  if ($line =~ /^Element:/) {
    #print $line,"\n";
    my @cols = split(/\s+/,$line);
    my @xy = @cols[11..12];
    $data->[$n]->{position} = \@xy;
    $data->[$n]->{channels} = extract_tdc_adc(@cols);
    $n++;
  } elsif ($line =~ /^Next-TAPS:/) {
    
    my @cols = split(/\s+/,$line);
    my @neighbours = @cols[3..$#cols];
    $data->[$m]->{neighbours} = \@neighbours;
    $m++;
  } elsif ($line =~ /^TAPSSG:/) {
     my @cols = split(/\s+/,$line);
     my $tdcadc = extract_tdc_adc(@cols);
     push(@{$data->[$l]->{channels}}, $tdcadc->[1]);
    $l++;
  }
}

#print Dumper($data->[0]);
#print Dumper($data->[16]);
#print join("", @lines_pbwo4);

for($n=0;$n<@{$data};$n++) {
  my @xy = @{$data->[$n]->{position}};
  my @channels = @{$data->[$n]->{channels}};
  my $neighbours = join(",", @{$data->[$n]->{neighbours}});
  printf("{%3d, {%7.3f, %7.3f}, %04d, %04d, %04d, %04d, %04d, {%s}}, \n", $n, @xy, $channels[0], $channels[1], $channels[1]+100, $channels[1]-200, $channels[1]-100, $neighbours);
}

close $fh;
