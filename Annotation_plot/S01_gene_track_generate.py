import os
from configparser import ConfigParser
import pandas as pd

PATH = "01_genomic_tracks"
os.makedirs(PATH, exist_ok=True)


config = ConfigParser()

# Universal kws
bigwig_kws = {
	"height": 2, "min_value": 0, "number_of_bins": 700, "nans_to_zeros": True, 
	"summary_method": "mean", "show_data_range": True, "file_type": "bigwig"
}

bed_kws = {
    "display": "collapsed", "border_color": "none",
    "labels": False, "file_type": "bed"
}

link_kws = {
    "height":2, "links_type": "arcs", "line_width": 1, "line_style": "solid", "color": "YlGnBu",
    "compact_arcs_level": 2, "use_middle": True, "file_type": "links"
}

space_index = 0
# Make config
#config.add_section("spacer")

## ReMap
config["ReMap"] = {
	"file": "/home/sdc1/Shify/Data_public/ENCODE/ReMap-2022/reMapDensity2022.bw",
        "Title": "ReMap-2022", "height": 4, "color": "orange", "min_value": 0
}
space_index += 1
config[f"spacer{space_index}"] = {"height": 0.5, "file_type":"spacer"}
## Bigwig


metadata = pd.read_table("/home/sdc1/Shify/Data_public/ENCODE/Immune-cells/CTCF_GRCh38_immunecells_ChIP/metadata_bigwig.tsv")
metadata.index = metadata.loc[:,"File accession"]


#cell_ids = ["ENCFF680XUD","ENCFF461QRZ","ENCFF877OBO","ENCFF609TMO","ENCFF843PZR","ENCFF044YVA","ENCFF557ZTZ","ENCFF113QRU","ENCFF655FYI","ENCFF451AYU","ENCFF395HSD","ENCFF663DNP"]
cell_ids = ["ENCFF680XUD","ENCFF843PZR","ENCFF044YVA","ENCFF557ZTZ","ENCFF113QRU","ENCFF655FYI","ENCFF451AYU"]
#for cell_id in metadata.index:
for cell_id in cell_ids: # Reordered the cell type
    cell = metadata.loc[cell_id,"Biosample term name"]
    config[f"CTCF {cell} ChIP"] = {
        "file": f"/home/sdc1/Shify/Data_public/ENCODE/Immune-cells/CTCF_GRCh38_immunecells_ChIP/{cell_id}.bw",
        "title": f"CTCF ChIP\n{cell}", "color": "lightcoral", **bigwig_kws
    }
    space_index += 1
    config[f"spacer{space_index}"] = {"height": 0.5, "file_type":"spacer"}


## Links
links = ['GM12878','T-cell','T-cell-activated','B-cell']
for link_id in links:
    config[f"{link_id} links"] = {
    "file": f"/home/sdc1/Shify/Data_public/ENCODE/Immune-cells/CTCF_GRCh38_immunecells_ChIP/{link_id}-5to30.links",
    #"file": f"/home/sdc1/Shify/Data_public/ENCODE/Immune-cells/{link_id}.links",
    "title": f"{link_id}", **link_kws
    }
    space_index += 1
    config[f"spacer{space_index}"] = {"height": 0.5, "file_type":"spacer"}


## SNP
config["SNP"] = {
	"file": "/home/sdc1/Shify/Data_public/GWAS/GWAS_Catalog/20220517_GWAS_Catalog_GRCh38/gwas_catalog_download_from_UCSC_modify.bed",
	"title": "GWAS Catalog SNPs", "height": 4, "fontsize": 8, "file_type": "bed"
}

space_index += 1
config[f"spacer{space_index}"] = {"height": 0.5, "file_type":"spacer"}

## Gene
config["Genes"] = {
    "file": "../../GENCODE/GRCh38/gencode.v42.chr_patch_hapl_scaff.annotation.gtf",
    "title": "Genes", "prefered_name": "gene_name", "merge_transcripts": True,
    "fontsize": 7, "height": 5, "labels": True, "max_labels": 100,
    "all_labels_inside": True, "style": "UCSC", "file_type": "gtf"
}


config["x-axis"] = {"fontsize": 10}

with open(f"{PATH}/tracks-use.ini", "w") as f:
    config.write(f)

#chrom, chromStart, chromEnd = tss.loc[gene, ["chrom", "chromStart", "chromEnd"]]
#chromStart -= 200000
##chromEnd += 200000
#suffix = "-additional" if additional else ""
# !pyGenomeTracks --tracks {PATH}/tracks-{gene}.ini --region {chrom}:{chromStart}-{chromEnd} \
#      -t 'Target gene: {gene}' --dpi 600 --width 22 --fontSize 10 \
#      --outFileName {PATH}/tracks-{gene}{suffix}.pdf 2> /dev/null
# !rm {PATH}/vis_peaks.bed {PATH}/*.links
