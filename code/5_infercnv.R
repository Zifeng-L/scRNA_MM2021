library(infercnv)
setwd("../MM_analysis/infercnv/output")

#MM01
infercnv_obj = CreateInfercnvObject(raw_counts_matrix="../MM_analysis/infercnv/input/infercnv_MM01_counts.txt",
                                    annotations_file="../MM_analysis/infercnv/input/infercnv_MM01_meta.txt",
                                    delim="\t",
                                    gene_order_file="../MM_analysis/infercnv/input/gencode_v38_gene_pos.txt",
                                    ref_group_names=c("normal"))
									
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                             out_dir="infercnv_MM01",  # dir is auto-created for storing outputs
                             cluster_by_groups=T,   # cluster
                             denoise=T,
                             HMM=T)
							 
seurat_obj = infercnv::add_to_seurat(infercnv_output_path="infercnv_MM01",
                                     top_n=10)
									 
#MM02
infercnv_obj = CreateInfercnvObject(raw_counts_matrix="../MM_analysis/infercnv/input/infercnv_MM02_counts.txt",
                                    annotations_file="../MM_analysis/infercnv/input/infercnv_MM02_meta.txt",
                                    delim="\t",
                                    gene_order_file="../MM_analysis/infercnv/input/gencode_v38_gene_pos.txt",
                                    ref_group_names=c("normal"))
									
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                             out_dir="infercnv_MM02",  # dir is auto-created for storing outputs
                             cluster_by_groups=T,   # cluster
                             denoise=T,
                             HMM=T)
							 
seurat_obj = infercnv::add_to_seurat(infercnv_output_path="infercnv_MM02",
                                     top_n=10)
									 
#MM03
infercnv_obj = CreateInfercnvObject(raw_counts_matrix="../MM_analysis/infercnv/input/infercnv_MM03_counts.txt",
                                    annotations_file="../MM_analysis/infercnv/input/infercnv_MM03_meta.txt",
                                    delim="\t",
                                    gene_order_file="../MM_analysis/infercnv/input/gencode_v38_gene_pos.txt",
                                    ref_group_names=c("normal"))
									
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                             out_dir="infercnv_MM03",  # dir is auto-created for storing outputs
                             cluster_by_groups=T,   # cluster
                             denoise=T,
                             HMM=T)
							 
seurat_obj = infercnv::add_to_seurat(infercnv_output_path="infercnv_MM03",
                                     top_n=10)
									 
#MM04
infercnv_obj = CreateInfercnvObject(raw_counts_matrix="../MM_analysis/infercnv/input/infercnv_MM04_counts.txt",
                                    annotations_file="../MM_analysis/infercnv/input/infercnv_MM04_meta.txt",
                                    delim="\t",
                                    gene_order_file="../MM_analysis/infercnv/input/gencode_v38_gene_pos.txt",
                                    ref_group_names=c("normal"))
									
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                             out_dir="infercnv_MM04",  # dir is auto-created for storing outputs
                             cluster_by_groups=T,   # cluster
                             denoise=T,
                             HMM=T)
							 
seurat_obj = infercnv::add_to_seurat(infercnv_output_path="infercnv_MM04",
                                     top_n=10)
									 
#MM05
infercnv_obj = CreateInfercnvObject(raw_counts_matrix="../MM_analysis/infercnv/input/infercnv_MM05_counts.txt",
                                    annotations_file="../MM_analysis/infercnv/input/infercnv_MM05_meta.txt",
                                    delim="\t",
                                    gene_order_file="../MM_analysis/infercnv/input/gencode_v38_gene_pos.txt",
                                    ref_group_names=c("normal"))
									
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                             out_dir="infercnv_MM05",  # dir is auto-created for storing outputs
                             cluster_by_groups=T,   # cluster
                             denoise=T,
                             HMM=T)
							 
seurat_obj = infercnv::add_to_seurat(infercnv_output_path="infercnv_MM05",
                                     top_n=10)
									 
#MM06
infercnv_obj = CreateInfercnvObject(raw_counts_matrix="../MM_analysis/infercnv/input/infercnv_MM06_counts.txt",
                                    annotations_file="../MM_analysis/infercnv/input/infercnv_MM06_meta.txt",
                                    delim="\t",
                                    gene_order_file="../MM_analysis/infercnv/input/gencode_v38_gene_pos.txt",
                                    ref_group_names=c("normal"))
									
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                             out_dir="infercnv_MM06",  # dir is auto-created for storing outputs
                             cluster_by_groups=T,   # cluster
                             denoise=T,
                             HMM=T)
							 
seurat_obj = infercnv::add_to_seurat(infercnv_output_path="infercnv_MM06",
                                     top_n=10)
									 
#MM07
infercnv_obj = CreateInfercnvObject(raw_counts_matrix="../MM_analysis/infercnv/input/infercnv_MM07_counts.txt",
                                    annotations_file="../MM_analysis/infercnv/input/infercnv_MM07_meta.txt",
                                    delim="\t",
                                    gene_order_file="../MM_analysis/infercnv/input/gencode_v38_gene_pos.txt",
                                    ref_group_names=c("normal"))
									
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                             out_dir="infercnv_MM07",  # dir is auto-created for storing outputs
                             cluster_by_groups=T,   # cluster
                             denoise=T,
                             HMM=T)
							 
seurat_obj = infercnv::add_to_seurat(infercnv_output_path="infercnv_MM07",
                                     top_n=10)
									 
#MM08
infercnv_obj = CreateInfercnvObject(raw_counts_matrix="../MM_analysis/infercnv/input/infercnv_MM08_counts.txt",
                                    annotations_file="../MM_analysis/infercnv/input/infercnv_MM08_meta.txt",
                                    delim="\t",
                                    gene_order_file="../MM_analysis/infercnv/input/gencode_v38_gene_pos.txt",
                                    ref_group_names=c("normal"))
									
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                             out_dir="infercnv_MM08",  # dir is auto-created for storing outputs
                             cluster_by_groups=T,   # cluster
                             denoise=T,
                             HMM=T)
							 
seurat_obj = infercnv::add_to_seurat(infercnv_output_path="infercnv_MM08",
                                     top_n=10)
									 
#MM09
infercnv_obj = CreateInfercnvObject(raw_counts_matrix="../MM_analysis/infercnv/input/infercnv_MM09_counts.txt",
                                    annotations_file="../MM_analysis/infercnv/input/infercnv_MM09_meta.txt",
                                    delim="\t",
                                    gene_order_file="../MM_analysis/infercnv/input/gencode_v38_gene_pos.txt",
                                    ref_group_names=c("normal"))
									
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                             out_dir="infercnv_MM09",  # dir is auto-created for storing outputs
                             cluster_by_groups=T,   # cluster
                             denoise=T,
                             HMM=T)
							 
seurat_obj = infercnv::add_to_seurat(infercnv_output_path="infercnv_MM09",
                                     top_n=10)
									 
#MM10
infercnv_obj = CreateInfercnvObject(raw_counts_matrix="../MM_analysis/infercnv/input/infercnv_MM10_counts.txt",
                                    annotations_file="../MM_analysis/infercnv/input/infercnv_MM10_meta.txt",
                                    delim="\t",
                                    gene_order_file="../MM_analysis/infercnv/input/gencode_v38_gene_pos.txt",
                                    ref_group_names=c("normal"))
									
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                             out_dir="infercnv_MM10",  # dir is auto-created for storing outputs
                             cluster_by_groups=T,   # cluster
                             denoise=T,
                             HMM=T)
							 
seurat_obj = infercnv::add_to_seurat(infercnv_output_path="infercnv_MM10",
                                     top_n=10)
									 
#MM11
infercnv_obj = CreateInfercnvObject(raw_counts_matrix="../MM_analysis/infercnv/input/infercnv_MM11_counts.txt",
                                    annotations_file="../MM_analysis/infercnv/input/infercnv_MM11_meta.txt",
                                    delim="\t",
                                    gene_order_file="../MM_analysis/infercnv/input/gencode_v38_gene_pos.txt",
                                    ref_group_names=c("normal"))
									
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                             out_dir="infercnv_MM11",  # dir is auto-created for storing outputs
                             cluster_by_groups=T,   # cluster
                             denoise=T,
                             HMM=T)
							 
seurat_obj = infercnv::add_to_seurat(infercnv_output_path="infercnv_MM11",
                                     top_n=10)
									 
#MM12
infercnv_obj = CreateInfercnvObject(raw_counts_matrix="../MM_analysis/infercnv/input/infercnv_MM12_counts.txt",
                                    annotations_file="../MM_analysis/infercnv/input/infercnv_MM12_meta.txt",
                                    delim="\t",
                                    gene_order_file="../MM_analysis/infercnv/input/gencode_v38_gene_pos.txt",
                                    ref_group_names=c("normal"))
									
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                             out_dir="infercnv_MM12",  # dir is auto-created for storing outputs
                             cluster_by_groups=T,   # cluster
                             denoise=T,
                             HMM=T)
							 
seurat_obj = infercnv::add_to_seurat(infercnv_output_path="infercnv_MM12",
                                     top_n=10)
									 
#MM13
infercnv_obj = CreateInfercnvObject(raw_counts_matrix="../MM_analysis/infercnv/input/infercnv_MM13_counts.txt",
                                    annotations_file="../MM_analysis/infercnv/input/infercnv_MM13_meta.txt",
                                    delim="\t",
                                    gene_order_file="../MM_analysis/infercnv/input/gencode_v38_gene_pos.txt",
                                    ref_group_names=c("normal"))
									
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                             out_dir="infercnv_MM13",  # dir is auto-created for storing outputs
                             cluster_by_groups=T,   # cluster
                             denoise=T,
                             HMM=T)
							 
seurat_obj = infercnv::add_to_seurat(infercnv_output_path="infercnv_MM13",
                                     top_n=10)
									 
#MM14
infercnv_obj = CreateInfercnvObject(raw_counts_matrix="../MM_analysis/infercnv/input/infercnv_MM14_counts.txt",
                                    annotations_file="../MM_analysis/infercnv/input/infercnv_MM14_meta.txt",
                                    delim="\t",
                                    gene_order_file="../MM_analysis/infercnv/input/gencode_v38_gene_pos.txt",
                                    ref_group_names=c("normal"))
									
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                             out_dir="infercnv_MM14",  # dir is auto-created for storing outputs
                             cluster_by_groups=T,   # cluster
                             denoise=T,
                             HMM=T)
							 
seurat_obj = infercnv::add_to_seurat(infercnv_output_path="infercnv_MM14",
                                     top_n=10)
									 
#MM15
infercnv_obj = CreateInfercnvObject(raw_counts_matrix="../MM_analysis/infercnv/input/infercnv_MM15_counts.txt",
                                    annotations_file="../MM_analysis/infercnv/input/infercnv_MM15_meta.txt",
                                    delim="\t",
                                    gene_order_file="../MM_analysis/infercnv/input/gencode_v38_gene_pos.txt",
                                    ref_group_names=c("normal"))
									
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                             out_dir="infercnv_MM15",  # dir is auto-created for storing outputs
                             cluster_by_groups=T,   # cluster
                             denoise=T,
                             HMM=T)
							 
seurat_obj = infercnv::add_to_seurat(infercnv_output_path="infercnv_MM15",
                                     top_n=10)
									 
#MM16
infercnv_obj = CreateInfercnvObject(raw_counts_matrix="../MM_analysis/infercnv/input/infercnv_MM16_counts.txt",
                                    annotations_file="../MM_analysis/infercnv/input/infercnv_MM16_meta.txt",
                                    delim="\t",
                                    gene_order_file="../MM_analysis/infercnv/input/gencode_v38_gene_pos.txt",
                                    ref_group_names=c("normal"))
									
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                             out_dir="infercnv_MM16",  # dir is auto-created for storing outputs
                             cluster_by_groups=T,   # cluster
                             denoise=T,
                             HMM=T)
							 
seurat_obj = infercnv::add_to_seurat(infercnv_output_path="infercnv_MM16",
                                     top_n=10)
									 
#MM17
infercnv_obj = CreateInfercnvObject(raw_counts_matrix="../MM_analysis/infercnv/input/infercnv_MM17_counts.txt",
                                    annotations_file="../MM_analysis/infercnv/input/infercnv_MM17_meta.txt",
                                    delim="\t",
                                    gene_order_file="../MM_analysis/infercnv/input/gencode_v38_gene_pos.txt",
                                    ref_group_names=c("normal"))
									
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                             out_dir="infercnv_MM17",  # dir is auto-created for storing outputs
                             cluster_by_groups=T,   # cluster
                             denoise=T,
                             HMM=T)
							 
seurat_obj = infercnv::add_to_seurat(infercnv_output_path="infercnv_MM17",
                                     top_n=10)
									 
#MM18
infercnv_obj = CreateInfercnvObject(raw_counts_matrix="../MM_analysis/infercnv/input/infercnv_MM18_counts.txt",
                                    annotations_file="../MM_analysis/infercnv/input/infercnv_MM18_meta.txt",
                                    delim="\t",
                                    gene_order_file="../MM_analysis/infercnv/input/gencode_v38_gene_pos.txt",
                                    ref_group_names=c("normal"))
									
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                             out_dir="infercnv_MM18",  # dir is auto-created for storing outputs
                             cluster_by_groups=T,   # cluster
                             denoise=T,
                             HMM=T)
							 
seurat_obj = infercnv::add_to_seurat(infercnv_output_path="infercnv_MM18",
                                     top_n=10)
									 
