#!/usr/bin/perl
# open GFF annotation file ...  see http://plasmodb.org/common/downloads/Current_Release/Pfalciparum3D7/gff/data/PlasmoDB-33_Pfalciparum3D7.gff

open (AFILE, $ARGV[0]);

while (<AFILE>) {

    ($Arm, $t, $type, $Start, $End) = split (' ');  # read in feature

    if ($type eq "gene") {   # gene contains all introns and exons... OK for us

	($t, $t2) = split (/\=/);
	($Key) = split (/\;/, $t2);

	$Hash{$Arm} .= "$Start  $End  ";   # store start + end
	$IDs{$Arm} .= "$Key  $Key  ";      # store gene name parsed from above

    }

}

open (BFILE, $ARGV[1]);   # open VCF/Connor-formatted file.  Only requirement is <chr/seq> <position>

while (<BFILE>) {

    ($Arm, $Pos) = split (' ');  # get required info

    @Genes = split (' ', $Hash{$Arm});    # parse all features on this chr/seq
    @IDs = split (' ', $IDs{$Arm});

    for ($i = 0; $i < @Genes; $i += 2) {   # check to see if this SNP is in that location

	if ($Pos >= $Genes[$i] && $Pos <= $Genes[$i + 1]) {

	    print "$IDs[$i]\n"

	}

    }

}
