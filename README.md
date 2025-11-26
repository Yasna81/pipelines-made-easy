# pipelines-made-easy
 A collection of python/R scripts that make life easier.


-----


### 1- DEG-analysis :

#### ⬇️ Ever needed a quick look at your DEGs without writing any code?


  - get the number of significant genes with their level of significance
  - find number of duplicates
  - find number of up/down-regulated genes
  - check the distribution of p-values
  - check the distribution of genes in the volcano plot

#### download libraries :

```bash
pip install pandas numpy matplotlib seaborn
```

##### run the Python script with your .csv file of DEGS:

``` bash
python DEG-analysis.py data.csv
```
#### Note 1:  Make sure you have columns named gene, PValue, and Foldchange. It is case sensitive!
#### Note 2:  You can manually adjust cutoffs in the script.

#### Sample output :
<img width="1237" height="900" alt="2025-10-25_22-23" src="https://github.com/user-attachments/assets/4849b7a2-8cdf-4aa7-be7b-1c3684f38138" />


----
-----


### 2- aligncheck.py :

#### ⬇️ Automizing samtool commands by giving a short, complete overview of your aligned reads.


- Comprehensive mapping quality analysis - not just averages, but full distribution across thresholds (0, ≤10, ≤30, >30)

- Complete CIGAR operation breakdown - see exact counts for matches, inserts, deletions, soft/hard clips

- Smart insert size calculation - automatically filters outliers and provides both mean & median

- Strand-specific mapping stats - forward vs reverse strand distribution at a glance

- All quality flags decoded - QC failures, duplicates, secondary/supplementary alignments in one view

#### download libraries :

```bash
conda install -c condaforge pysam
```

##### run the Python script with your BAM file :

``` bash
python aligncheck.py x.bam
```

##### Example output :
<img width="1614" height="866" alt="2025-11-26_15-5226" src="https://github.com/user-attachments/assets/bf2e0a02-7809-49f0-b951-5a6849dc1bd5" />





Inspired by the book: Bioinformatics Data Skills by Vince Buffalo


----
-----


