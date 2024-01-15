#!/usr/bin/perl
use strict;

# Expected input from stdin:
#II  6122    7608    intergenic_II_6125_7605 +   II  7605    7733    YBL108C-A   -
#II  6122    7608    intergenic_II_6125_7605 +   II  5790    6125    YBL109W +
#II  7730    8180    intergenic_II_7733_8177 +   II  7605    7733    YBL108C-A   -
#II  7730    8180    intergenic_II_7733_8177 +   II  8177    8482    YBL108W +
#II  8479    8851    intergenic_II_8482_8848 +   II  8177    8482    YBL108W +
#II  9089    9096    intergenic_II_9092_9093 +   .   -1  -1  .   .
#II  9421    9428    intergenic_II_9424_9425 +   .   -1  -1  .   .
#II  9515    9586    intergenic_II_9518_9583 +   II  9583    9666    tL(UAA)B1   +
#II  9663    9964    intergenic_II_9666_9961 +   II  9961    10551   YBL107C -

my %NEW_ID;
while(<>){
    chomp;
    my ($chr, $e5, $e3, $id, $str, $b_chr, $b_e5, $b_e3, $b_id, $b_str) = split /\t/;
    unless (exists($NEW_ID{$id})){
        $NEW_ID{$id}->{"downstream"} = 'NA';
        $NEW_ID{$id}->{"downstream_str"} = 'NA';
        $NEW_ID{$id}->{"upstream"} = 'NA';
        $NEW_ID{$id}->{"upstream_str"} = 'NA';
    }

    if ($b_e3 > $e3 && $b_e3 > 0){
        $NEW_ID{$id}->{"downstream"} = $b_id;
        $NEW_ID{$id}->{"downstream_str"} = $b_str;
    } elsif ($b_e5 < $e5 && $b_e5 > 0){
        $NEW_ID{$id}->{"upstream"} = $b_id;
        $NEW_ID{$id}->{"upstream_str"} = $b_str;
    } 
}

# Print header
print "id\tup_gene\tup_strand\tdown_gene\tdown_strand\n";

# Print genes up and downstream
foreach my $id (keys %NEW_ID){
    my ($dn, $up, $dn_str, $up_str) = ('NA','NA','NA','NA');
    if($NEW_ID{$id}->{"downstream"} ne 'NA'){
        $dn = $NEW_ID{$id}->{"downstream"};
        $dn_str = $NEW_ID{$id}->{"downstream_str"};
    }
    if($NEW_ID{$id}->{"upstream"} ne 'NA'){
        $up = $NEW_ID{$id}->{"upstream"};
        $up_str = $NEW_ID{$id}->{"upstream_str"};
    }
    print join("\t",($id,$up,$up_str,$dn,$dn_str))."\n";
}



