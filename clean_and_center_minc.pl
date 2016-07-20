#!/usr/bin/env perl

use strict;
use warnings "all";
use File::Temp qw/ tempdir /;

# Directory for temporary files.
my $tmpdir = &tempdir( "clean-XXXXXX", TMPDIR => 1, CLEANUP => 1 );

my $file = $ARGV[0];
my $output = $ARGV[1];

&make_regular( $file, $output );

sub make_regular {

    my $input = shift;
    my $output = shift;

    # bypass symbolic link, if any
    &run( 'cp', '-f', $input, "${tmpdir}/regular.mnc" );
    # convert to minc2 since minc_modify_header can crash on
    # minc1 files (don't know why).
    my $ret = `file ${tmpdir}/regular.mnc`;
    chomp( $ret ) ;
    if( $ret =~ m/NetCDF/ ) {
      &run( 'mincconvert', '-2', "${tmpdir}/regular.mnc",
            "${tmpdir}/regular_v2.mnc" );
      &run( 'mv', '-f', "${tmpdir}/regular_v2.mnc", "${tmpdir}/regular.mnc" );
    }

    &run( 'minc_modify_header', '-sinsert', 'xspace:spacing=regular__',
          '-sinsert', 'yspace:spacing=regular__', '-sinsert',
          'zspace:spacing=regular__', "${tmpdir}/regular.mnc" );
    &run( 'mincreshape', '-quiet', '-clobber', '+direction', '-dimorder',
          'zspace,yspace,xspace', '-dimsize', 'xspace=-1', '-dimsize', 'yspace=-1',
          '-dimsize', 'zspace=-1', "${tmpdir}/regular.mnc", $output );
    &run( 'volcentre', $output );

    unlink( "${tmpdir}/regular.mnc" );
  }

  # Execute a system call.

  sub run {
    print "@_\n";
    system(@_)==0 or die "Command @_ failed with status: $?";
  }
