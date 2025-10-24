import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt 
import warnings
import os
import seaborn as sns
class GeneExpresionAnalyzer :
    """
    Analysis of Expresion for your Data 

    note : make sure your cols are : gene, PValue and Foldchange
    """
    print(">>>  Analysis of Expresion for your Data  <<<<    note : make sure your cols are : gene, PValue and Foldchange")
    def __init__(self,df,gene_col = "gene",pval_col = "PValue",lfc_col = "Foldchange") :
        self.df = df.copy()
        self.gene_col = gene_col
        self.pval_col = pval_col
        self.lfc_col = lfc_col
        self.analysis_results = {}
    
    def pvalue_analysis(self):
        """analysis of pvlue"""
        if self.pval_col not in self.df.columns:
            print(f" P-value column '{self.pval_col}' not found")
            return
        
        print("\n" + "=" * 70)
        print(" P-VALUE ANALYSIS")
        print("=" * 70)


        pvals = self.df[self.pval_col].dropna()
        significance_levels = {
            'p < 0.001': (pvals < 0.001).sum(),
            'p < 0.01': (pvals < 0.01).sum(),
            'p < 0.05': (pvals < 0.05).sum(),
            'p < 0.1': (pvals < 0.1).sum(),
            'Not significant': (pvals >= 0.05).sum()
        }
        print("  Significance Levels:")

        for level, count in significance_levels.items():
            percentage = (count / len(pvals)) * 100
            print(f"    {level:<15}: {count:>6} genes ({percentage:5.1f}%)")
    



    def foldchange_analysis(self):
        """Analyze log2FoldChange distribution"""
        if self.lfc_col not in self.df.columns:
            print(f" Log2FoldChange column '{self.lfc_col}' not found")
            return
        
        #print("= " * 70)
        print(" FOLD CHANGE ANALYSIS")
        #print=("=" * 70)
        lfc = self.df[self.lfc_col].dropna()
        fc_categories = {
            'Strong Down (LFC < -2)': (lfc < -2).sum(),
            'Moderate Down (-2 ≤ LFC < -1)': ((lfc >= -2) & (lfc < -1)).sum(),
            'Mild Down (-1 ≤ LFC < -0.5)': ((lfc >= -1) & (lfc < -0.5)).sum(),
            'Stable (-0.5 ≤ LFC ≤ 0.5)': ((lfc >= -0.5) & (lfc <= 0.5)).sum(),
            'Mild Up (0.5 < LFC ≤ 1)': ((lfc > 0.5) & (lfc <= 1)).sum(),
            'Moderate Up (1 < LFC ≤ 2)': ((lfc > 1) & (lfc <= 2)).sum(),
            'Strong Up (LFC > 2)': (lfc > 2).sum()
        }
        print("  Fold Change Distribution:")
        for category, count in fc_categories.items():
            percentage = (count / len(lfc)) * 100
            print(f"    {category:<30}: {count:>6} genes ({percentage:5.1f}%)")
        
        #print("\n" + "=" * 70)
        print(" DIFFERENTIAL EXPRESSION SUMMARY")
        #print("=" * 70)

        de_genes = self.df[
            (self.df[self.pval_col] < 0.05) &
            (self.df[self.lfc_col].abs() > 1)
        ]
        up_regulated = de_genes[de_genes[self.lfc_col] > 0]
        down_regulated = de_genes[de_genes[self.lfc_col] < 0 ]

        de_stats = {
            'Total DE genes': len(de_genes),
            'Up-regulated': len(up_regulated),
            'Down-regulated': len(down_regulated),
            'Up/Down ratio': f"{len(up_regulated)/len(down_regulated):.2f}" if len(down_regulated) > 0 else "N/A"
        }

        for stat, value in de_stats.items():
            print(f" {stat:<20}: {value}")   

    def quality_metrics(self):
        """Calculate data quality metrics"""
        print("\n" + "=" * 70)
        print("DATA QUALITY METRICS")
        print("=" * 70)
        
        quality_metrics = {
            'Total missing values': self.df.isnull().sum().sum(),
            'Columns with missing data': (self.df.isnull().sum() > 0).sum(),
            'Duplicate rows': self.df.duplicated().sum(),
            'Constant columns': (self.df.nunique() == 1).sum()
        }
        for metric, value in quality_metrics.items():
            print(f"    {metric:<30}: {value}")
    def show_all_plots(self):
        """Show all available plots in a grid"""
        available_plots = []
        
        if self.pval_col in self.df.columns:
            available_plots.append('pvalue')
        if self.lfc_col in self.df.columns:
            available_plots.append('foldchange')
        if self.pval_col in self.df.columns and self.lfc_col in self.df.columns:
            available_plots.append('volcano')
        if self.df.isnull().sum().sum() > 0:
            available_plots.append('missing')
        
        n_plots = len(available_plots)
        
        if n_plots == 0:
            print(" No data available for plotting")
            return
        
        # Create subplot grid
        if n_plots == 1:
            fig, ax = plt.subplots(1, 1, figsize=(10, 8))
            axes = [ax]
        elif n_plots == 2:
            fig, axes = plt.subplots(1, 2, figsize=(15, 6))
        elif n_plots == 3:
            fig, axes = plt.subplots(2, 2, figsize=(15, 10))
            axes = axes.flatten()[:3]  # Use only first 3
        else:  # 4 plots
            fig, axes = plt.subplots(2, 2, figsize=(15, 10))
            axes = axes.flatten()
        
        fig.suptitle('Gene Expression Analysis Plots', fontsize=16, fontweight='bold')
        
        plot_index = 0
        
        # P-value distribution
        if 'pvalue' in available_plots:
            axes[plot_index].hist(self.df[self.pval_col].dropna(), bins=50, alpha=0.7, color='skyblue')
            axes[plot_index].axvline(0.05, color='red', linestyle='--', label='p=0.05')
            axes[plot_index].set_xlabel('P-value')
            axes[plot_index].set_ylabel('Frequency')
            axes[plot_index].set_title('P-value Distribution')
            axes[plot_index].legend()
            axes[plot_index].grid(True, alpha=0.3)
            plot_index += 1
        
        # Fold change distribution
        if 'foldchange' in available_plots:
            axes[plot_index].hist(self.df[self.lfc_col].dropna(), bins=50, alpha=0.7, color='lightgreen')
            axes[plot_index].axvline(0, color='black', linestyle='-')
            axes[plot_index].axvline(1, color='red', linestyle='--', label='LFC=±1')
            axes[plot_index].axvline(-1, color='red', linestyle='--')
            axes[plot_index].set_xlabel('Log2FoldChange')
            axes[plot_index].set_ylabel('Frequency')
            axes[plot_index].set_title('Fold Change Distribution')
            axes[plot_index].legend()
            axes[plot_index].grid(True, alpha=0.3)
            plot_index += 1
        
        # Volcano plot
        if 'volcano' in available_plots:
            df_clean = self.df.dropna(subset=[self.pval_col, self.lfc_col])
    
        # Create separate datasets
        up_reg = df_clean[(df_clean[self.pval_col] < 0.05) & (df_clean[self.lfc_col] > 1)]
        down_reg = df_clean[(df_clean[self.pval_col] < 0.05) & (df_clean[self.lfc_col] < -1)]
        not_sig = df_clean[(df_clean[self.pval_col] >= 0.05) | (df_clean[self.lfc_col].abs() <= 1)]
    
        # Plot each group separately with labels
        axes[plot_index].scatter(not_sig[self.lfc_col], -np.log10(not_sig[self.pval_col]), 
                           alpha=0.6, s=20, c='gray', label='Not significant')
        axes[plot_index].scatter(up_reg[self.lfc_col], -np.log10(up_reg[self.pval_col]), 
                           alpha=0.7, s=20, c='red', label=f'Up-regulated ({len(up_reg)})')
        axes[plot_index].scatter(down_reg[self.lfc_col], -np.log10(down_reg[self.pval_col]), 
                           alpha=0.7, s=20, c='blue', label=f'Down-regulated ({len(down_reg)})')
    
        # Add thresholds
        axes[plot_index].axhline(-np.log10(0.05), color='black', linestyle='--', label='p=0.05')
        axes[plot_index].axvline(1, color='orange', linestyle='--')
        axes[plot_index].axvline(-1, color='orange', linestyle='--')
    
        axes[plot_index].set_xlabel('Log2FoldChange')
        axes[plot_index].set_ylabel('-log10(P-value)')
        axes[plot_index].set_title('Volcano Plot')
        axes[plot_index].legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        axes[plot_index].grid(True, alpha=0.3)
        plot_index += 1
        
        # Missing data heatmap
        if 'missing' in available_plots:
            missing_data = self.df.isnull()
            sns.heatmap(missing_data, cbar=True, cmap='viridis', ax=axes[plot_index])
            axes[plot_index].set_title('Missing Data Heatmap')
        
        plt.tight_layout()
        plt.show()


    def generate_report(self):
        """Generate complete analysis report"""
        print(" GENERATING COMPREHENSIVE GENE EXPRESSION ANALYSIS REPORT")
        print("=" * 70)
        
        self.pvalue_analysis()
        self.foldchange_analysis() 
        self.quality_metrics()
        self.show_all_plots()
        
        print("\n" + "=" * 70)
        print(" ANALYSIS COMPLETE!")
        print("=" * 70)
 

#------------------------------------------------------------------------ usage :
def analyze_gene_expression_data(file_path, **kwargs):
    """
    Main function to analyze gene expression data
    
    Parameters:
    file_path (str): Path to CSV file
    **kwargs: Additional arguments for pd.read_csv()
    """
    # Load data
    print(f"📁 Loading data from: {file_path}")
    try:
        df = pd.read_csv(file_path, **kwargs)
        print(f" Successfully loaded data with {df.shape[0]:,} rows and {df.shape[1]} columns")
    except Exception as e:
        print(f" Error loading data: {e}")
        return
    
    # Initialize analyzer
    analyzer = GeneExpresionAnalyzer(df)
    # Generate report
    analyzer.generate_report()
    
    return analyzer




if __name__ == "__main__":
    # Example usage
    file_path = "~/practice/project-1/CrePosHyperoxiavsCreNegHyperoxia_all_detectable_genes.csv"  # Change this to your file path
    
    # Analyze the data
    analyzer = analyze_gene_expression_data(
        file_path,
        # Add any additional pandas read_csv parameters here
        # sep='\t',  # for tab-separated files
        # comment='#',  # skip lines starting with #
    )
