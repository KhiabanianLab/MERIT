#!/usr/bin/perl -w

$MERITFile = $ARGV[0];
open (IN, "$MERITFile") || die "Cannot open $MERITFile.\n";
<IN>;

printf STDOUT "##fileformat=VCFv4.2\n";
printf STDOUT "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele frequency\">\n";
printf STDOUT "##INFO=<ID=AP1,Number=1,Type=Integer,Description=\"Allele depth forward\">\n";
printf STDOUT "##INFO=<ID=AP2,Number=1,Type=Integer,Description=\"Allele depth reverse\">\n";
printf STDOUT "##INFO=<ID=DP1,Number=1,Type=Integer,Description=\"Total depth forward\">\n";
printf STDOUT "##INFO=<ID=DP2,Number=1,Type=Integer,Description=\"Total depth reverse\">\n";
printf STDOUT "##INFO=<ID=AQ1,Number=1,Type=Integer,Description=\"Average allele quality forward\">\n";
printf STDOUT "##INFO=<ID=AQ2,Number=1,Type=Integer,Description=\"Average allele quality reverse\">\n";
printf STDOUT "##INFO=<ID=INP1,Number=1,Type=Integer,Description=\"Allele position in forward reads\">\n";
printf STDOUT "##INFO=<ID=INP2,Number=1,Type=Integer,Description=\"Allele position in reverse reads\">\n";
printf STDOUT "##INFO=<ID=RCTX,Number=.,Type=String,Description=\"Reference sequence context\">\n";
printf STDOUT "##INFO=<ID=ACTX,Number=.,Type=String,Description=\"Allele sequence context\">\n";
printf STDOUT "##INFO=<ID=DB,Number=.,Type=String,Description=\"Bases in indels\">\n";
printf STDOUT "##INFO=<ID=NHP,Number=1,Type=Integer,Description=\"Number of homopolymers in indel's context\">\n";
printf STDOUT "##INFO=<ID=ANN,Number=.,Type=String,Description=\"Annotation\">\n";

printf STDOUT "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

while (<IN>) {
  @Line = split (/\t/, substr($_, 0, -1));

  
  $chr = $Line[0];
  $pos = $Line[1];
  $ref = $Line[2];
  $ref_ctx = $Line[3];
  $alt = $Line[4];
  $alt_ctx = $Line[5];
  $total_depth_frw = $Line[6];
  $total_depth_rev = $Line[7];
  $alt_depth_frw = $Line[8];
  $alt_depth_rev = $Line[9];
  $freq = $Line[10];
  $ave_alt_qual_frw = $Line[11];
  $ave_alt_qual_rev = $Line[12];
  $pos_inread_alt_frw = $Line[13];
  $pos_inread_alt_rev = $Line[14];
  $indel_bases = $Line[15];
  $hp_repeats = $Line[16];
  $snp = ".";
  $effect = $Line[17];

  printf STDOUT "$chr\t$pos\t$snp\t$ref\t$alt\t.\t.\t";
  printf STDOUT "AF=$freq";
  printf STDOUT ";";
  printf STDOUT "AP1=$alt_depth_frw";
  printf STDOUT ";";
  printf STDOUT "AP2=$alt_depth_rev";
  printf STDOUT ";";
  printf STDOUT "DP1=$total_depth_frw";
  printf STDOUT ";";
  printf STDOUT "DP2=$total_depth_rev";
  printf STDOUT ";";
  printf STDOUT "AQ1=$ave_alt_qual_frw";
  printf STDOUT ";";
  printf STDOUT "AQ2=$ave_alt_qual_rev";
  printf STDOUT ";";
  printf STDOUT "INP1=$pos_inread_alt_frw";
  printf STDOUT ";";
  printf STDOUT "INP2=$pos_inread_alt_rev";
  printf STDOUT ";";
  printf STDOUT "RCTX=$ref_ctx";
  printf STDOUT ";";
  printf STDOUT "ACTX=$alt_ctx";
  printf STDOUT ";";
  printf STDOUT "DB=$indel_bases";
  printf STDOUT ";";
  printf STDOUT "NHP=$hp_repeats";
  printf STDOUT ";";
  printf STDOUT "ANN=$effect";
  printf STDOUT "\n";
}
