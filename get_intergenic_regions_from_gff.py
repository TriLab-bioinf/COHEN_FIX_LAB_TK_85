import re

# grep -w '^#\|gene' Saccharomyces_cerevisiae.R64-1-1.106.gtf |sort -k1,1 -k4n,4 > Saccharomyces_cerevisiae.R64-1-1.106_sorted.gtf
gtf = open("./saccharomyces_cerevisiae_R64-4-1_20230830.gff", "r")

(old_chr, old5, old3) = ('',0,0)
counter = 0
for line in gtf.readlines():

    # Comment line
    if re.search(pattern=r"^#", string=line):
        print(line, end="")

    # Feature line
    elif re.search(pattern=r"\t(gene|tRNA_gene|ncRNA_gene|pseudogene|rRNA_gene|snRNA_gene|snoRNA_gene|tRNA_gene|telomerase_RNA_gene|transposable_element_gene|long_terminal_repeat)\t", string=line):
        gene_list = re.split(pattern="\t", string=line)
        feat_id = gene_list[2]
        end5 = int(gene_list[3])
        end3 = int(gene_list[4])
        chr = gene_list[0]

        # Reset old coords if it is a new chromosome
        if chr != old_chr:
            (old5, old3) = (0,0)

        #print(f"{feat_id}")

        #if feat_id != "gene":
        #    continue 

        #print(f"if {chr} == {old_chr} and {old3} < {end5}:")

        # Check if new and old genes are on same chromosome
        if chr == old_chr and old3 < end5:
            #print(f"{chr}")
            gene_id = f"intergenic_{chr}_{old3}_{end5}"
            print(f"{chr}\tsgd\tintergenic\t{old3}\t{end5}\t.\t+\t.\tgene_id {gene_id}")

        old_chr = chr
        old5 = min(end5,old5)
        old3 = max(end3,old3)

gtf.close()
