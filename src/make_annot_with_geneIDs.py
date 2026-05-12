"""
the enriched genes are gene IDs, but the eggnogmapper decorated gff is transcript IDs
read the gff file and replace all the transcript IDs with the corresponding geneIDs
"""


import parse_gff as gff


if __name__ == "__main__":
    
    annot_dir = "/Users/miltr339/work/c_maculatus"
    funcannot_filepath = f"{annot_dir}/C_mac_eggnog_diamond.emapper.annotations"
    funcannot_filepath_out = f"{annot_dir}/C_mac_eggnog_diamond.emapper.annotations_geneIDs"
    gff_filepath = f"{annot_dir}/Cmac_Lome_yes_yTor.gff"

    gff_dict = gff.parse_gff3_general(gff_filepath)

    with open(funcannot_filepath, "r") as infile_annot, open(funcannot_filepath_out, "w") as outfile_annot:
        for func_line in infile_annot.readlines():
            if func_line[0] == "#":
                outfile_annot.write(func_line)
            else:
                lines_list = func_line.strip().split("\t")
                transcript = lines_list[0]
                annots = "\t".join(lines_list[1:])
                geneID = gff_dict[transcript].parent_id
                newline = f"{geneID}\t{annots}\n"
                outfile_annot.write(newline)

