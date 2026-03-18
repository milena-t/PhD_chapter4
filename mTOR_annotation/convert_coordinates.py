"""
Transform the coordinates of genes from the Kaufmann2023 assembly onto the superscaffolded assembly with the agp file
I will not write a general script, since I am only interested in three genes from utg000322l_1 which has a + orientation so I don't need to worry about that.
the agp file with the coordinates is here: "PhD_chapter4/mTOR_annotation/data/SALSA_superscaffolding_contig_coordinates.agp"
"""


def transform_coords(in_path:str, out_path:str, agp_new, agp_old):
    contig_ind = 0
    start_ind = 3
    end_ind = 4
    with open(in_path, "r") as tor_annot, open(out_path, "w") as tor_scaffold:
        for line in tor_annot.readlines():
            line_new = line.strip().split("\t")
            line_old = line.strip().split("\t")
            if line_old[contig_ind] == agp_old["cont"]:
                line_new[contig_ind] = agp_new["cont"]
                line_new[start_ind] = str(int(line_old[start_ind])-agp_old["start"])
                line_new[end_ind] = str(int(line_old[end_ind])-agp_old["start"])
            write_line_new = "\t".join(line_new)
            tor_scaffold.write(f"{write_line_new}\n")
    print(f"transformed coords written to {out_path}")



if __name__ == "__main__":

    username ="miltr339"
    tor_annot_dir = f"/Users/{username}/work/chapter4/mTor_annot/"
    tor_annots = [f"{tor_annot_dir}yTor-{letter}.gff" for letter in ["A", "B", "C"]]

    agp_str = """scaffold_26	1	5935189	1	W	utg000322l_1	2496969	8432157	+
                 scaffold_48	1	2496968	1	W	utg000322l_1	1	2496968	+"""
    agp_new = { "cont" : "scaffold_26",
                "start" : 1,
                "end" : 5935189}
    agp_old = { "cont" : "utg000322l_1",
                "start" : 2496969,
                "end" : 8432157}


    for tor_annot_path in tor_annots:
        tor_scaffold_annot_path = tor_annot_path.replace(".gff", "_superscaffolds.gff")
        transform_coords(in_path=tor_annot_path, out_path=tor_scaffold_annot_path, agp_new=agp_new, agp_old=agp_old)
                