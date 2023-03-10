# PhyloPhind

Command line tool for filtering FASTA files in the EXPAT_BENCH curation pipeline. 

# Install instructions 

1. Clone this repository to your desktop

```
git clone https://github.com/SebPorras/PhyloPhind.git
```

2. Create a conda environment 

```
conda create -n phylo_phind python=3.9
```

3. Activate environment 

```
conda activate phylo_phind 
```

4. Install the required Python packages 

```
pip install -r ./PhyloPhind/requirements.txt
```

5. Install binfpy 

http://bioinf1.scmb.uq.edu.au/opensource/binfpy

6. Make binfpy accessible to your Python PATH via terminal config 

```
export PYTHONPATH=/home/seb-porras/binfpy
```

# Running the script 

1. Place PhyloPhind.py in expat_bench directory 

Now you can run the script like so 

```
python ./expat_bench/PhyloPhind.py -e 3_5_2_6 -t 0.15 -m 0 -r 0
```

2. Filtered FASTA will be placed in ./expat_bench/workflows/EC_NUM/files/EC_NUM_filt.fasta

3. For help, run

```
python ./expat_bench/PhyloPhind.py --help
```

Output: 

```
  - EC_NUM: -e , --ec_num
  -> Enter as 3_5_2_6 etc.
  
  - THRESHOLD: -t, --threshold
  -> Jaccard threshold between 0 and 1
  
  - MIN_SEQS: -m, --min_seqs
  -> minimum number of sequences required for a InterPro tag group to be considered, defaults to 0
  
  - ROW_NUM: -r , --row_num 
  -> Choose group to filter on, sorted based on largest InterPro tag groups indexed at 0, defaults to zero(largest group)
```

4. A log of the process will be saved to ./expat_bench/filtering_logs/EC_NUM_filtering_log.txt

```
EC Group: 3_5_2_6 

Threshold: 0.15 
Min seqs in group: 0 
Row number of group: 0

-Inital number of sequences 
>129 
-Inital families present 
>Metallo-beta-lactamase superfamily, Class-B beta-lactamase family 
>Class-A beta-lactamase family 
>Class-C beta-lactamase family 
>Metallo-beta-lactamase superfamily 
>Class-D beta-lactamase family 

-Number of sequences after filtering 
>80 
-Families present after filtering 
>Class-A beta-lactamase family 
>Class-C beta-lactamase family 
>Class-D beta-lactamase family 

-Number of sequences after family sampling 
>81
```


# Changing location of output

All paths can by changing the path constants at the bottom of PhyloPhind.py

Important paths are:

- WORKDIR: defines path from current directory to files in expat_bench 
- FASTA: defines path for where filtered fasta file will be saved 


