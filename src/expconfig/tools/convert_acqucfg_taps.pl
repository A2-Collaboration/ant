#!/usr/bin/perl
use strict;
use warnings;
use Data::Dumper;

if (@ARGV != 2) {
  die "Need non-sensitive and sensitive detector config file for TAPS\n";
}

my $inputfile   = $ARGV[0];
my $inputfile_S = $ARGV[1];
open(my $fh, "<$inputfile") or die "Cannot open $inputfile: $!\n";
open(my $fh_S, "<$inputfile_S") or die "Cannot open $inputfile_S: $!\n";

sub is_pbwo4 {
  my $n = shift;
  return $n % 73 < 12;
}

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
  if ($line =~ /Element:/) {
    #print $line,"\n";
    my @cols = split(/\s+/,$line);
    my @xy = @cols[11..12];
    $data->[$n]->{position} = \@xy;
    $data->[$n]->{channels} = extract_tdc_adc(@cols);
    $n++;
  } elsif ($line =~ /Next-Neighbour:/) {
    
    my @cols = split(/\s+/,$line);
    my @neighbours = @cols[3..$#cols];
    $data->[$m]->{neighbours} = \@neighbours;
    $m++;
  } elsif ($line =~ /TAPSSG:/) {
    # skip pbwo4 dummy entries
    if (!is_pbwo4($l)) {
      my @cols = split(/\s+/,$line);
      my $tdcadc = extract_tdc_adc(@cols);
      push(@{$data->[$l]->{channels}}, $tdcadc->[1]);
    }
    $l++;
  }
}

# add the sensitive channels,
# two channels for BaF2 (long/short gate)
# but only one for PbWO4
$n=0;
$l=0;
while (my $line = <$fh_S>) {
  chomp $line;
  if ($line =~ /Element:/) {
    #print $line,"\n";
    my @cols = split(/\s+/,$line);
    my $tdcadc = extract_tdc_adc(@cols);
    push(@{$data->[$n]->{channels}}, $tdcadc->[1]);
    $n++;
  } elsif ($line =~ /TAPSSG:/) {
    # skip pbwo4 dummy entries
    if (!is_pbwo4($l)) {
      my @cols = split(/\s+/,$line);
      my $tdcadc = extract_tdc_adc(@cols);
      push(@{$data->[$l]->{channels}}, $tdcadc->[1]);
    }
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
  if(is_pbwo4($n)) {
    printf("{%3d, {%7.3f, %7.3f}, %04d, %04d, %04d, {%s}}, \n", $n, @xy, @channels, $neighbours);
  }
  else {
    #printf("{%3d, {%7.3f, %7.3f}, %04d, %04d, %04d, %04d, %04d, {%s}}, \n", $n, @xy, @channels, $neighbours);
  }
}

close $fh;
close $fh_S;
