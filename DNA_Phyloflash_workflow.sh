# Phyloflash workflow for DNA (Metagenomes)

# Starting in the Metagenomes data folder
~/Working/IS19_Meta/Genomics


# # Create Filtered folder (assuming you're still in /home/ctrivedi/Working/Meta/Genomics)
# mkdir 02_Filtered


# # Generate a .txt file with sample names and locations. Do this from inside the raw data folder.
# # Download these files, and put into columns "sample", "r1", and "r2", and save as .txt then reupload.

# ls *R1* > ../R1s.txt # for R1s
# ls *R2* > ../R2s.txt # for R2s
# for f in *_R1_001.fastq.gz;         do SAMPLE=`basename ${f%%_S*_*_R1_001*}`;                 echo "$SAMPLE"; done; # for sample names
# # Opened R1s and R2s, added heading, and copy/pasted names from terminal and saved over R1 file then renamed to 'sample_names_for_filtering.txt'


# # Create config files before running the filtering program.
# # HAVE TO RUN THIS FROM THE DATA FOLDER OR IT WON'T CAPTURE THE CORRECT DATA INPUT FOLDER.
# iu-gen-configs ../sample_names_for_filtering.txt -o ../02_Filtered/

# # Filter one at a time
# # iu-filter-quality-minoche 02_Filtered/10_F.ini # OR

# # Filter using a for loop (run from /Genomics/)
# for ini in ../02_Filtered/*.ini; do iu-filter-quality-minoche $ini; done

# # Check sample stats after running (besides looking at individual stat files)
# grep -E 'number of pairs analyzed|total pairs passed' ../02_Filtered/*STATS.txt > ../All_passed_stats.txt


# Create a new sample list for phyloflash (from inside the data folder.
cd ~/Working/IS19_Meta/Genomics/02_Filtered/

for f in *_R1.fastq;         
do SAMPLE=`basename ${f%%-QUALITY*}`;                 
echo "$SAMPLE"; 
done;
# copy and paste this inside nano and save 'pf_names.txt'


# Will use a loop to run through the filtered files
# Use phyloFlash_DNA_loop.sh, or

for SAMPLE in `awk '{print $1}' ~/Working/IS19_Meta/Genomics/pf_names.txt`
do
    if [ "$SAMPLE" == "sample" ]; then continue; fi
    
    # you need to make sure your "ls ~/Working/IS19_Meta/Genomics/02_Filtered/*QUALITY_PASSED_R1*" returns R1 files for all your samples in samples.txt
    R1s=`ls ~/Working/IS19_Meta/Genomics/02_Filtered/$SAMPLE*QUALITY_PASSED_R1* | python -c 'import sys; print ",".join([x.strip() for x in sys.stdin.readlines()])'`
    R2s=`ls ~/Working/IS19_Meta/Genomics/02_Filtered/$SAMPLE*QUALITY_PASSED_R2* | python -c 'import sys; print ",".join([x.strip() for x in sys.stdin.readlines()])'`
    
    /home/ctrivedi/miniconda3/envs/pf/bin/phyloFlash.pl -lib $SAMPLE -read1 $R1s -read2 $R2s -CPUs 16 -readlength 250 -taxlevel 20 -id 80 -tophit -zip
done



# Unzip all the tarballs after!
for a in *.tar.gz
do
    a_dir=`expr $a : '\(.*\).tar.gz'`
    mkdir $a_dir 2>/dev/null
    tar -xvzf $a -C $a_dir
done


# Find NTU files and copy to new directory
mkdir ~/Working/IS19_Meta/Genomics/phyloflash_NTUabundance_files

cp */*phyloFlash.NTUabundance.csv ../phyloflash_NTUabundance_files/ # From within ~/Working/IS19_Meta/Genomics/phyloflash/


# Take analysis into R
# NTU files copied to C:\Users\ctrivedi\Google Drive\PostDoc\1_Projects\3_RNA_Preservation_Methods\Analyses\Chris_results


