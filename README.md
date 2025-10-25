# pipelines-made-easy
 A collection of python/R scripts that make life easier.


-----


### 1- DEG-analysis :

#### ⬇️ Ever needed a quick look at your DEGs with no code?


  - get the number of significant genes with their level of significance
  - find number of duplicates
  - find number of up/down down-regulated genes
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
#### Note:  make sure you have columns named gene, PValue, and Foldchange. It is case sensitive !

#### Sample output :
<img width="1237" height="900" alt="2025-10-25_22-23" src="https://github.com/user-attachments/assets/4849b7a2-8cdf-4aa7-be7b-1c3684f38138" />


----




