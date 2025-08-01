# LEI_LAB_TK_208

### Check potential restriction enzyme used
```
./scrips/check_potential_restriction_enzyme.sh ./data_reference/RumpKD_RP1_Biotin_R1_001.fastq.gz 4 1000000
   2830 AATT
   3192 AATA
   3564 ATTA
   4040 AAAT
   4063 AAAA
   4742 ATAA
   4936 ATTT
   5495 ATAT
  20639 GATC <<<<<<<< DpnII
```

## Juicer HiC pipeline

### Make restriction site file for custom dme assembly
### script downloaded and adapted from https://github.com/aidenlab/juicer/blob/main/misc/generate_site_positions.py to include path to Dme reference fasta file
```
python ./scripts/generate_site_positions.py DpnII dme
```

### Make bwa database
```
bwa index ./data_reference/dmel-all-chromosome-r6.52.fasta
```

### Make chromosome size file
```
seqtk comp ./dmel-only-chromosomes-r6.52.fasta |cut -f 1,2 > ./dmel-only-chrom-sizes.txt
```

### QC and Trim reads with fastp using snakemake hic_pipeline.smk 
```
sbatch ./run_snakemake.sh ./workflow/hic_pipeline.smk
```

### Run juicer pipeline
```
module load juicer/1.6
juicer.sh -y ./data_reference/dme_DpnII.txt -z ./data_reference/dmel-all-chromosome-r6.52.fasta -p ./data_reference/dmel_chrom_sizes.txt
```
### Virtual 4C analysis with R

All analysis was don in R (v4.3.3) with the following R quarto script [virtual-4C-analysis.qmd]("LEI_LAB_TK_208/R/render/virtual-4C-analysis.qmd"). Auxiliary R scripts are located in [R/scripts]("LEI_LAB_TK_208/R/scripts").

R environment with the packages and verions used for the analysis can be regenerated using the `renv.lock` file by running the following command in R.

If you have not initialized `renv` in your project
```
renv::init()
```   
This will set up the project library and, if an `renv.lock` file exists, it will prompt you to restore the environment.

If `renv` is already initialized and active in the project (i.e., you see the `renv/` directory and `.Rprofile` file), simply opening the project should automatically load the `renv` environment based on the `renv.lock` file. You can confirm this by checking the active library paths by typing `.libPaths()` in your R console.

