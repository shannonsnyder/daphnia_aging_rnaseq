#!/bin/bash

#SBATCH -A bgmp       ### Account                                                                                                              
#SBATCH --partition=memory   ### Partition                                                                        \                            
#SBATCH --job-name=sort_bams         ### Job Name                                                                                              
#SBATCH --time=24:00:00        ### WallTime                                                                                                    
#SBATCH --nodes=1              ### Number of Nodes                                                                                             
#SBATCH --ntasks=1             ### Number of tasks per array job                                                                               
##SBATCH --array=0-3           ### Array index                                                                                                 
#SBATCH --mail-type=END,FAIL              ### Mail events (NONE, BEGIN, END, FA$                                                               
##SBATCH --mail-user=ssnyder3@uoregon.edu  ### Where to send mail
#SBATCH --cpus-per-task=12  ### Number of CPU cores per task 

# Path to the QoRTs jar file
QORTS_JAR="/home/ssnyder3/nereus/aging_rnaseq/star_alignment/pulex/QoRTs/QoRTs.jar"

# Input directory containing sorted BAM files
INPUT_DIR="/home/ssnyder3/nereus/aging_rnaseq/star_alignment/pulex/aligned_reads/sorted_bams"

# Output directory
OUTPUT_DIR="/home/ssnyder3/nereus/aging_rnaseq/star_alignment/pulex/QoRTs"

# Annotation directory                                                                                               
GTF="/home/ssnyder3/nereus/aging_rnaseq/genomes/dpulex/ncbi-genomes-2023-05-22/Fixed_GCF_021134715.1_ASM2113471v1_genomic.gtf"

# Loop through each BAM file in the input directory
for bam_file in "$INPUT_DIR"/*.bam; do
    # Extract the file name without the extension
    filename=$(basename -- "$bam_file")
    filename_noext="${filename%.*}"
    
    # Create the output directory for this sample
    output_subdir="$OUTPUT_DIR/$filename_noext"
    mkdir -p "$output_subdir"
    
    # Run QoRTs on the input BAM file
    java -Xmx4G -jar "$QORTS_JAR" QC --generatePlots "$bam_file" "$GTF" "$output_subdir"
    
    # Rename the QoRTs output files based on the input file name
   # mv "$output_subdir/QC_BAMSummary.txt" "$output_subdir/${filename_noext}_QC_BAMSummary.txt"
   # mv "$output_subdir/QC_Metrics.txt" "$output_subdir/${filename_noext}_QC_Metrics.txt"
    
    echo "Processed $filename"
done

echo "QoRTs processing and renaming complete."
