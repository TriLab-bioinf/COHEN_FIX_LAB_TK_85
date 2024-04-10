#!/usr/bin/perl
use strict;

my $usage = "$0 -g <annot gtf file> -d <samtools depth file> -w <window size [200]> -n <number of samples [12]> -D <read-depth cutoff [1]>\n\n";
my %arg = @ARGV;
die $usage unless $arg{-g} && $arg{-d};

my $NUM_SAMPLES=$arg{-n} || 12;
my $WINDOW = $arg{-w} || 200;
my $DEPTH_CUTOFF = $arg{-D} || 1;

#                 <----------- direction of read-depth reading
#   key = chr + e5
#   new_coord    e5  gene_1      e3 
#       | WINDOW |-------------->|     
# -------------------------------------------------------------------------------------------- chromosome
#                                      |<-----------| WINDOW |             
#                                      e5  gene_2   e3       new_coord
#                                                   key = chr + e3
#
#  direction of read-depth reading  ---------->


# Upload gene coordinates from annotation
open (my $gtf, "<$arg{-g}") || die "ERROR, I cannot open $arg{-g}: $!\n\n";

my (%genes);
while(<$gtf>){

    my ($chr, undef, $feat, $e5, $e3, undef, $strand, undef, $comm) = split(m/\t/);
    # skip all feats except gene
    next unless $feat eq "gene";
    
    # Extranct gene_id info from $comm
    my ($gene_id, $gene_name, $biotype) = &parse_comm($comm);

    
    # Only keep protein_coding genes
    next unless $biotype eq "protein_coding";
    
    # Record coords of 5' region (gene 5' end + window_size) within %genes using key = chr : lowest coord
    if ($strand eq "+"){
        my $new_coord = $e5 - $WINDOW;
        my $key = "$chr:$new_coord";
        $genes{$key} = [$gene_id, $new_coord, $e5, $strand];
    } elsif ($strand eq "-"){
        my $new_coord = $e3 + $WINDOW;
        my $key = "$chr:$e3";
        $genes{$key} = [$gene_id, $e3, $new_coord, $strand];

        # Debug line
        #if($gene_id eq "YLR406C-A"){
        #print STDERR "$gene_id, $gene_name, $chr, $feat, $e5, $e3, $strand, new coord = $new_coord, e3 = $e3\n";
        #}

    } else {next}

    

}
close $gtf;

# Open depth file
open (my $depth, "<$arg{-d}") || die "ERROR, I cannot open $arg{-g}: $!\n\n";
    # If current position within any 5' gene coords then track length of
    # non-zero runs before reaching the end for each track (sample)
my $c = 0;
while(<$depth>){
    chomp;
    my ($chr, $pos, @val) = split /\t/;
    next unless $pos > 0;

    my $start_key = "$chr:$pos";
    if ($genes{$start_key}){
        $c++;
        my $gene_id = $genes{$start_key}->[0];
        my $strand = $genes{$start_key}->[3];
        my $end_key = $genes{$start_key}->[2];
        #print "$gene_id, $strand, $end_key\n";
        my @upstream_region;

        # Load depth for region upstream gene
        while($pos < $end_key){
            ($chr, $pos, @val) = split(m/\t/,  <$depth>);
            push @upstream_region, "$chr\t$pos\t". join("\t", @val);

        }
        if ($strand eq "+"){
            # Reverse upstream_region array so we track depth starting from the 5'end of the CDS
            @upstream_region = reverse(@upstream_region);
        }

        # Count lenght of UTR for each sample
        my @utr_length = &count_utr_length(@upstream_region);
        print "$gene_id\t$start_key\tLengths\t".join("\t",@utr_length)."\n";

        #print "$gene_id, $strand, $start_key, $end_key\n@upstream_region\n";
        #print "@utr_length\n\n";
        #exit(0) if $c ==  1;
    }
    
}
close $depth;


##################################################################################
sub count_utr_length {
    my @region = @_;
    #print "@region\n";
    my @len = split(//, "0" x $NUM_SAMPLES);     # (0,0,0,0,0,0,0,0);
    my @flag = split(//, "1" x $NUM_SAMPLES);     #(1,1,1,1,1,1,1,1);
    #print "FLAG = @flag\n";
    foreach my $row (@region){
        chomp($row);
        my @row = split(/\t/, $row);
        @row = @row[2..scalar(@row)-1];
        #print "\n\n\@row = ".join("-",@row)."\n";
        #print "flag=".join("-",@flag)."\n";
        my $p = 0;
        foreach my $sample (@row){
            # Set flag = 0 if read-depth falls below depth cutoff 
            $flag[$p] = 0 if $sample < $DEPTH_CUTOFF;
            $len[$p]++ if ($flag[$p] == 1);
            $p++;
        }
    }
    return @len;
}

sub parse_comm {
    my $desc = shift;
    my $gene_id = $1 if m/gene_id\s"(\S+?)";/ || "none";
    my $gene_name = $1 if m/gene_name\s"(\S+?)";/ || "none";
    my $biotype = $1 if m/gene_biotype\s"(\S+?)";/ || "none";
    return ($gene_id, $gene_name, $biotype);
}
