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
pip install -r requirements.txt
```

5. Install binfpy 

http://bioinf1.scmb.uq.edu.au/opensource/binfpy

6. Make binfpy accessible to your Python PATH via terminal config 

```
export PYTHONPATH=/home/seb-porras/binfpy
```

# Running the script 

Place PhyloPhind.py in expat_bench directory 

Now you can run the script like so 

```
python PhyloPhind.py -e 3_5_2_6 -t 0.2 -m 1 -r 1
```

Filtered FASTA will be placed in ./workflows/EC_NUM/files/EC_NUM_filt.fasta


