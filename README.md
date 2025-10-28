# LncRNAs-in-Dryophytes-arenicolor
This repository includes code for the tools used in the assembly of brain transcriptomes of three different life cycle stages of _Dryophytes arenicolor_ (Hylidae), and both known and de novo lncRNA identifcation. The clean sequencing files and the assembled transcriptomes can be found in the NCBI BioProject Repository **PRJNA1295574**. Each sample is formed by pooling three independent tissues. The stages of the samples are Pre-metamorphosis (G26), Metamorphic Climax (G42), and Adults. 

This repository is the code for an original manuscript by Herrera-Orozco et al currently under peer-review by G3:Genes|Genomes|Genetics. 


# The order of the scripts should go

1. quality_and_trimming
2. trinity_assembly
3. cd_hit
4. annotated_lncRNAs
5. g_profiler
6. coding_filter
7. CPAT
8. pfam
9. dryophytes_filter
10. Salmon
11. wgcna
12. wgcna_enrichment
13. enrichment_plots

We also included a comprehensive list of software used, their version, and the links where they can be obtained (Software_version.txt) and some specifications of the PC used (pc_specs.txt). 

**The authors do not own any of the tools, software, or databases used in this pipeline. Please refer to each repository for authorship** 


