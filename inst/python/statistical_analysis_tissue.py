#!usr/bin/env/python3
import pandas as pd
import numpy as np
from collections import defaultdict
import json
import re
from tqdm import tqdm
from glob import glob
from scipy.stats import pearsonr
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.font_manager


class Downstream_SHAP_Analysis:
    def __init__(self, pred_target, moa_pair = None, n_top = 20):
        """
        * 1. Summarize SHAP values of all features of machine learning models on specific mode-of-action dataset.
            This function will also generate a .tsv file for plotting figures in R.
        * 2. Subtract SHAP values of molecular markers of ddr genes.
        * 3. Analysis tissue specificity of top genes
        * 4. Analysis the interactions between top genes 
        
        -------------------------------------------------------------------------------------------
        Parameters:
            pred_target: prediction target ('aoc' or 'bliss')
                a string

            moa_pair: the mode-of-action(MOA) subset to analysis (e.g. 'ATMi_ATRi')
                if None, carry our analysis on overall dataset.
                None or a string

            n_top: number of top genes for interaction analysis
                an integer

        Yields:
            self.pred_target: pred_target
            
            self.moa_pair: moa_pair

            self.n_top:

            self.fpaths: input path of SHAP analysis results from model_training.py
                a list of string

            self.processed_fpaths: input path of SHAP analysis results, with SHAP value of every gene calculated.
                a list of string

            self.cat_outpath: output path of categoricalized SHAP analysis results.
                a string

            self.mol_outpath: output path of molecular marker's SHAP analysis results.
                a string

            self.top_gene_outpath: output path of top gene's SHAP analysis results.
                a string

            self.fig_outpath: output path of interaction heatmap of top features.
                a string

            self.category_df: categoricalized SHAP analysis results.
                a Pandas dataframe

            self.molecular_df: molecular marker's SHAP analysis results.
                a Pandas dataframe

            self.top_gene_df: top gene's SHAP analysis results.
                a Pandas dataframe
            self.top_genes: top genes
                a list of string
        """
        self.pred_target = pred_target
        self.moa_pair = moa_pair
        self.n_top = n_top
        if self.moa_pair == None:
            self.outpath = ''
        else:
            self.outpath = 'tissue/'+moa_pair+'/'

        # Check if SHAP analysis results are in the right place
        self.SHAPpath = self.outpath+'SHAP_'+self.pred_target+'_*.csv'
        assert len(glob(self.SHAPpath)) > 0, 'No SHAP Analysis results of machine learning models found!'

        self.fpaths = glob(self.SHAPpath)
        self.gene_SHAPpath = self.outpath+'gene_SHAP_'+self.pred_target+'_*.csv'
        #if len(glob(self.gene_SHAPpath)) < 5:
        if len(glob(self.gene_SHAPpath)) < 1:
            print("Preprocess SHAP features ... ")
            self.Preprocess_and_dissect_SHAP(self.fpaths)

        self.gene_SHAPpath = self.outpath+'gene_SHAP_'+self.pred_target+'_*.csv'
        self.cat_SHAPpath = self.outpath+'cat_SHAP_'+self.pred_target+'_*.csv'
        self.chem_SHAPpath = self.outpath+'chem_SHAP_'+self.pred_target+'_*.csv'
        self.mol_SHAPpath = self.outpath+'mol_SHAP_'+self.pred_target+'_*.csv'

        self.gene_fpaths = glob(self.gene_SHAPpath)
        self.cat_fpaths = glob(self.cat_SHAPpath)
        self.chem_fpaths = glob(self.chem_SHAPpath)
        self.mol_fpaths = glob(self.mol_SHAPpath)
        
        # Sum SHAP of datasets and gives overall SHAP importance on the whole dataset
        self.totalSHAP_outpath = self.outpath+self.pred_target+'_totalSHAP.tsv'
        self.gene_totalSHAP_outpath = self.outpath+self.pred_target+'_gene_totalSHAP.tsv'
        self.cat_totalSHAP_outpath = self.outpath+self.pred_target+'_cat_totalSHAP.tsv'
        self.chem_totalSHAP_outpath = self.outpath+self.pred_target+'_chem_totalSHAP.tsv'
        self.mol_totalSHAP_outpath = self.outpath+self.pred_target+'_mol_totalSHAP.tsv'
        self.gene_each_mol_totalSHAP_outpath = self.outpath+self.pred_target+'_gene_each_mol_totalSHAP.tsv'
        
        ##
        self.top_gene_outpath = self.outpath+self.pred_target+'_top_gene.tsv'
        self.fig_outpath = self.outpath+self.pred_target+'_interaction_SHAP.pdf'
        ##

        if len(glob(self.totalSHAP_outpath)) <= 0:
            self.totalSHAP_df = self.Sum_SHAP_on_dataset(self.fpaths, annote_category=True)
            self.totalSHAP_df.to_csv(self.totalSHAP_outpath, sep = '\t', index = False)
        else:
            self.totalSHAP_df = pd.read_csv(self.totalSHAP_outpath, sep = '\t')

        if len(glob(self.gene_totalSHAP_outpath)) <= 0:
            self.gene_totalSHAP_df = self.Sum_SHAP_on_dataset(self.gene_fpaths, annote_category=False)
            self.gene_totalSHAP_df.to_csv(self.gene_totalSHAP_outpath, sep = '\t', index = False)
        else:
            self.gene_totalSHAP_df = pd.read_csv(self.gene_totalSHAP_outpath, sep = '\t')
        
        if len(glob(self.cat_totalSHAP_outpath)) <= 0:
            self.cat_totalSHAP_df = self.Sum_SHAP_on_dataset(self.cat_fpaths, annote_category=False)
            self.cat_totalSHAP_df.to_csv(self.cat_totalSHAP_outpath, sep = '\t', index = False)
        else:
            self.cat_totalSHAP_df = pd.read_csv(self.cat_totalSHAP_outpath, sep = '\t')

        if len(glob(self.chem_totalSHAP_outpath)) <= 0:
            self.chem_totalSHAP_df = self.Sum_SHAP_on_dataset(self.chem_fpaths, annote_category=False)
            self.chem_totalSHAP_df.to_csv(self.chem_totalSHAP_outpath, sep = '\t', index = False)
        else:
            self.chem_totalSHAP_df = pd.read_csv(self.chem_totalSHAP_outpath, sep = '\t')
        
        if len(glob(self.mol_totalSHAP_outpath)) <= 0:
            self.mol_totalSHAP_df = self.Sum_SHAP_on_dataset(self.mol_fpaths, annote_category=False)
            self.mol_totalSHAP_df.to_csv(self.mol_totalSHAP_outpath, sep = '\t', index = False)
        else:
            self.mol_totalSHAP_df = pd.read_csv(self.mol_totalSHAP_outpath, sep = '\t')
        
        if len(glob(self.gene_each_mol_totalSHAP_outpath)) <= 0:
            self.gene_each_mol_totalSHAP_df = self.Annote_each_molecular_features_of_genes()
            self.gene_each_mol_totalSHAP_df.to_csv(self.gene_each_mol_totalSHAP_outpath, sep = '\t', index = False)
        else:
            self.gene_each_mol_totalSHAP_df = pd.read_csv(self.gene_each_mol_totalSHAP_outpath, sep = '\t')

        if len(glob(self.top_gene_outpath)) <= 0:
            self.top_gene_df = self.Get_top_n_genes()
            self.top_gene_df.to_csv(self.top_gene_outpath, sep = '\t', index = False)
        else:
            self.top_gene_df = pd.read_csv(self.top_gene_outpath, sep = '\t')

        self.top_genes = self.top_gene_df['gene'].to_list()
        #self.top_genes = [re.sub('\*$', '', gene) for gene in self.top_gene_df['gene']]

        #df = self.Analyse_tissue_specificity_of_top_features()
        #self.tissue_specificity_outpath = self.outpath+self.pred_target+'-SHAP-tissue-specificity_of_top_genes.tsv'
        #df.to_csv(self.tissue_specificity_outpath, sep = '\t')
        
        self.Analyse_top_feature_interation_network()

    def feature_to_category(self,f):
        """
        * Annote feature by category.

        Parameters:
            f: name of input feature.
                a string
    
        Yields:
            cat: category of the input feature.
                a string
        """

        cat = 'Unknown'

        if (f == 'Cell_line'):
            cat = 'Cell_line'
        elif (f == 'Synergy_batch'):
            cat = 'Synergy batch'
        elif (f == 'Cancer_type'):
            cat = 'Cancer_type'
        elif (f == 'Cancer_subtype'):
            cat = 'Cancer_subtype'
        elif (f == "Concentration"):
            cat = 'Concentration of Static Compound'
        elif re.match(r'Treatment_[1|2]_Oncolead_[0-9]{3}', f) is not None:
            cat = 'Monotherapy drug efficacy (AOC)'
        elif re.match(r'Treatment_[1|2]_ave', f) is not None:
            cat = 'Monotherapy drug efficacy (AOC)'
        elif re.match(r'Treatment_[1|2]_moa', f) is not None:
            cat = 'Mode-of-action'
        elif re.match(r'Treatment_[1|2]_name', f) is not None:
            cat = 'Drug name'
        elif re.match(r'Treatment_[1|2]_RDK_[0-9]+', f) is not None:
            cat = 'Chemical_Structure_RDK'
        elif re.match(r'Treatment_[1|2]_MACCS_[0-9]+', f) is not None:
            cat = 'Chemical_Structure_MACCS'
        elif re.match(r'Treatment_[1|2]_Morgan_[0-9]+', f) is not None:
            cat = 'Chemical_Structure_Morgan'
        elif re.match(r'Treatment_[1|2]_FP2_[0-9]+', f) is not None:
            cat = 'Chemical_Structure_FP2'
        elif re.match(r'Treatment_[1|2]_FP3_[0-9]+', f) is not None:
            cat = 'Chemical_Structure_FP3'
        elif re.match(r'Treatment_[1|2]_FP4_[0-9]+', f) is not None:
            cat = 'Chemical_Structure_FP4'
        elif  f.endswith('_snv'):
            cat = 'Molecular_snv'
        elif f.endswith('_cnv'):
            cat = 'Molecular_cnv'
        elif f.endswith('_exp'):
            cat = 'Molecular_exp'
        elif f.endswith('_lof'):
            cat = 'Molecular_lof'
        elif f.endswith('__coh_pat'):
            cat = 'Molecular_coh_pat'
        elif f.endswith('__lof_pat'):
            cat = 'Molecular_lof_pat'
        elif f.endswith('_ddr'):
            cat = 'Molecular_ddr'
        elif f.endswith('_all'):
            cat = 'Molecular_all'
        elif f.startswith('Geneset_'):
            cat = 'Geneset'
        elif f.startswith('CRISPR_cancer_dependency_'):
            cat = 'CRISPR Cancer Dependency'
        else:
            print("Unknown feature found: "+f)
    
        return cat

    def Preprocess_and_dissect_SHAP(self, fpaths):
        """
        * Preprocess SHAP value results from models of each fold.

        Parameters:
            self.fpaths

        Yields:
            processed SHAP values, saved to <processed>_fpath
        """
        def Sum_all_mol_of_gene(df):
            """
            * Sum all molecular markers (cnv, snv, exp and lof) of this gene
            
            Parameters:
                df: SHAP value table
                    a Pandas dataframe
            Yields:
                df: processed SHAP value table with extra columns: <gene>_all_mol
                    a Pandas dataframe
            """
            target_genes = open('../../feature/target_gene/target_genes.txt', 'r').read().rstrip().split('\n')
            genes = []
            all_col = df.columns
            for col in all_col:
                if self.feature_to_category(col) in ['Molecular_exp', 'Molecular_cnv', 'Molecular_snv', 'Molecular_lof']:
                    genes.append(col.split('_')[0])
            
            df_new = {}
            genes = sorted(set(genes))
            for g in genes:
                gene_col = [col for col in all_col if col.split('_')[0] == g]
                if g in target_genes:
                    g_new = g+"*"
                else:
                    g_new = g
                df_new.update({g_new:df[gene_col].sum(axis = 1)})
            
            df_new = pd.DataFrame.from_dict(df_new)
            return df_new
        
        def Sum_category(df):
            """ SHAP of each category
                
                Mode-of-action: 'Mode-of-action'
                Drug-name: 'Drug name'
                Monotherapy: 'Monotherapy drug efficacy (AOC)'
                Molecular Biomarkders: 'Molecular_*'
                Geneset Annotations: 'Geneset'
                Chemical Structure Fingerprints: 'Chemical_Structure_*'
            
            Params
            ------
            df:
            
            Yields
            ------
            df_new: SHAP by category
            """
            categories = {"Mode-of-action": r'Mode-of-action',"Drug name": r'Drug*',
                "Monotherapy": r'Monotherapy*',
                "Molecular Biomarkders": r'Molecular_*',
                "Geneset Annotations": r'Geneset',
                "Chemical Structure": r'Chemical_Structure_*'}
            df_new = {}
            all_col = df.columns
            for k,v in categories.items():
                cat_cols = [col for col in all_col if re.match(v,self.feature_to_category(col)) is not None]
                df_new.update({k:df[cat_cols].sum(axis = 1)})

            df_new = pd.DataFrame.from_dict(df_new)
            return df_new

        def Sum_each_chemical(df):
            """ SHAP if each chemical strucure fingerprints
                6: MACCS, Morgan, RDK, FP2, FP3 FP4
            """
            df_new = {}
            all_col = df.columns

            # SHAP of six types of fingerprints respectively
            for chem_fp in ['Chemical_Structure_RDK', 'Chemical_Structure_MACCS', 'Chemical_Structure_Morgan', 'Chemical_Structure_FP2', 'Chemical_Structure_FP3', 'Chemical_Structure_FP4']:
                chem_col = [col for col in all_col if self.feature_to_category(col) == chem_fp]
                df_new.update({chem_fp.split("_")[-1]: df[chem_col].sum(axis = 1)})
            
            df_new = pd.DataFrame.from_dict(df_new)
            return df_new

        def Sum_each_molecular(df):
            """ SHAP of each mol biomarkers
                4+2+1
            """
            df_new = {}
            all_col = df.columns
            # Sum of each types of molecular biomarkers respectively
            for mol_marker in ['Molecular_exp', 'Molecular_cnv', 'Molecular_snv', 'Molecular_lof', 'Molecular_coh_pat', 'Molecular_lof_pat', 'Molecular_ddr']:
                mol_col = [col for col in all_col if self.feature_to_category(col) == mol_marker]
                df_new.update({"_".join(mol_marker.split("_")[1:]): df[mol_col].sum(axis = 1)})

            df_new = pd.DataFrame.from_dict(df_new)
            return df_new

        for path in fpaths:
            idx = path.split('_')[-1].split('.')[0]
            data = pd.read_csv(path, header = 0)
            # summed by gene
            data_gene = Sum_all_mol_of_gene(data)
            data_gene.to_csv(self.outpath+'gene_SHAP_'+self.pred_target+'_'+idx+'.csv', index = False)
            # summed by category
            data_cat = Sum_category(data)
            data_cat.to_csv(self.outpath+'cat_SHAP_'+self.pred_target+'_'+idx+'.csv', index = False)
            # chemical-sum of each fingerprints
            data_chem = Sum_each_chemical(data)
            data_chem.to_csv(self.outpath+'chem_SHAP_'+self.pred_target+'_'+idx+'.csv', index = False)
            # molecular-sum of each type of molecular (4+2+1)
            data_mol = Sum_each_molecular(data)
            data_mol.to_csv(self.outpath+'mol_SHAP_'+self.pred_target+'_'+idx+'.csv', index = False)
            #outpath = self.outpath+'processed_SHAP_'+self.pred_target+'_'+idx+'.csv'
            #data.to_csv(outpath, index = False)

    def Get_top_n_genes(self):
        """ 
        * Get n top genes
        
        Parameters:
            self.molecular_df
            self.n_top
            self.top_gene_outpath

        Yields:
            top_genes: top n genes selected and their mean SHAP values across all folds
                a Pandas dataframe
        """
        df = self.gene_totalSHAP_df[['feature', 'SHAP_val']].rename(columns = {'feature':'gene'})
        top_gene_df = df.groupby('gene', as_index=False).agg('mean').nlargest(self.n_top, 'SHAP_val')
        return top_gene_df


    
    def Sum_SHAP_on_dataset(self, fpaths, annote_category = True):
        """ 
        * Annotate SHAP feature set by feature categories.
            This function will generate a .tsv file for plotting figures in R.  
        
        Parameters:
            self.fpaths
            annote_category: boolean, True then annothe category

        Yields:
            category_df: Summarized SHAP contributions by category.
                a Pandas dataframe
                written to <pred_target>+"_SHAP.tsv"
        """

        print('Start prepare SHAP dataframe for '+self.pred_target+' score on overall dataset...')
        if annote_category:
            df = {'feature':[], 'SHAP_val':[], 'rep':[], 'feature_type': []} 
        else:
            df = {'feature':[], 'SHAP_val':[], 'rep':[]}
        
        for path in fpaths:
            data = pd.read_csv(path, header = 0)
            rep = str(path.split('_')[-1].split('.')[0]) #replcate index
            
            for col in tqdm(data.columns, total = len(data.columns.to_list())):
                feature = col #feature type
                shap_val = np.mean([abs(i) for i in data[col]]) #SHAP value: average of absolue shape value across all samples
                df['feature'].append(feature)
                df['SHAP_val'].append(shap_val)
                df['rep'].append(rep)
                if annote_category:
                    cat = self.feature_to_category(feature)
                    df['feature_type'].append(cat)

        df = pd.DataFrame.from_dict(df)
        return df

    def Annote_each_molecular_features_of_genes(self):
        """ 
        * Dissect the SHAP results of four typesf molecular markers of ddr genes.
        
        Parameters:
            self.pred_target
            self.totalSHAP_df

        Yields:
            molecular_df: SHAP analysis of ddr genes' molecular markers
                a Pandas dataframe
                written to <pred_target>+"_molecular_SHAP.tsv"
        """
        print('Annote molecular biomarker types of all genes for '+self.pred_target+' score prediction models ...')

        target_genes = open('../../feature/target_gene/target_genes.txt', 'r').read().rstrip().split('\n')
        data_mol = self.totalSHAP_df[self.totalSHAP_df['feature_type'].isin(['Molecular_exp', 'Molecular_cnv', 'Molecular_snv', 'Molecular_lof'])]
        df = {'gene':[], 'mol_type':[], 'SHAP_val':[], 'rep':[], 'target':[]}

        for _,r in data_mol.iterrows():
            gene, mol_type = r['feature'].split('_')
            shap_val = r['SHAP_val']
            rep = r['rep']
            df['mol_type'].append(mol_type)
            df['SHAP_val'].append(shap_val)
            df['rep'].append(rep)
            if gene in target_genes: # whether is a direct target of drug
                df['target'].append('True')
                df['gene'].append(gene+'*')
            else:
                df['target'].append('False')
                df['gene'].append(gene)
        
        df = pd.DataFrame.from_dict(df)
        return df

    def Analyse_tissue_specificity_of_top_features(self):
        """
        * Show tissue-specificity of top genes.
            1. to be determined
            (cancer_type and cancer_subtype)

            2. show which tissue type correlates well with which top genes
        
        Parameters:
            self.gene_fpaths
            self.top_genes
            self.pred_target

        Yields:
            df: tissue specificity based on SHAP value correlation
            
        """
        print('Start analysing tissue-specificity of top genes ...')
        all_data = []
        for path in glob(self.gene_SHAPpath):
            idx = path.split('_')[-1].split('.')[0]
            
            # corresponding tissue
            if self.moa_pair == None:
                testpath = '../../test_by_cell_line/fold_'+idx+'/Test.tsv'
            else:
                testpath = '../../test_by_cell_line/fold_'+idx+'/'+self.moa_pair+'.tsv'
            Test = pd.read_csv(testpath, sep = '\t')
            Test = Test[Test['.response_'+self.pred_target].notna()]
            
            data = pd.read_csv(path, header = 0)
            data = data[self.top_genes]
            data['cancer_type'] = Test['.metadata_cancer_type']
            data['cancer_subtype'] = Test['.metadata_cancer_subtype']
            all_data.append(data)
        data = pd.concat(all_data)
        df = {'gene':[], 'shap_val':[],'|shap_val|':[], 'cancer_type':[], 'cancer_subtype':[]}
        for g in self.top_genes:
            df['shap_val'].extend(data[g].to_list())
            df['|shap_val|'].extend(abs(data[g]).to_list())
            df['cancer_type'].extend(data['cancer_type'].to_list())
            df['cancer_subtype'].extend(data['cancer_subtype'].to_list())
            df['gene'].extend([g for i in range(len(data[g].to_list()))])
        
        df = pd.DataFrame.from_dict(df)
        print(df)
        return df


    def Analyse_top_feature_interation_network(self):
        """
        * Draw interaction map between top genes

        Parameters:
            self.pred_target
            self.propcessed_fpaths
            self.category_df
            self.top_genes

        Yields:
            correlation heatmaps of top genes
        """

        #combine SHAP Analysis results from all cross-validation models
        all_data = []
        for path in glob(self.gene_SHAPpath):
            data = pd.read_csv(path, header = 0)
            all_data.append(data) 
        
        data = pd.concat(all_data)
        cols = list(self.top_genes)
        data_gene = data[cols]
        #data_gene.rename(columns = {gene+'_all':gene for gene in self.top_genes}, inplace = True)
    
        cor_matrix = data_gene.corr()
        # output correlation matrix of top features 
        outpath = self.pred_target+'_top_gene_interaction.tsv'
        cor_matrix.to_csv(outpath, sep = '\t', header = True, index = True)
        
        # plot correlation heatmap between features

        fig = plt.figure(dpi=300, figsize = (20,20))
        sns.set_style({'font.family':'serif', 'font.serif':'Helvetica'})
        p = sns.heatmap(cor_matrix, square = True,cmap = 'seismic', xticklabels=1, yticklabels=1,vmin = -1, vmax = 1, center = 0, cbar_kws={'label': "Pearson's corrlation", 'shrink':0.5})
        p.tick_params(labelsize = 8)
        #p.figure.axes[-1].yaxis.label.set_size(8)
        p.axes.set_title("Correlation Heatmap of Top Molecular Biomarkers", fontdict={'fontsize':20, 'fontweight': 'medium'})
        fig = p.get_figure()
        #fig.savefig(rep+'_'+self.fig_outpath)
        fig.savefig(self.fig_outpath)
        plt.clf()


def main():
        
    # 1. Yield dataframe of SHAP contribution by categories 
    # and
    # 2. Yield dataframe of SHAP contributions of ddr genes
    # and
    # 3. Check tissue-specific ddr gene features.
    # and
    # 4. interaction of top genes.
    Downstream_SHAP_Analysis(pred_target = 'aoc')
    Downstream_SHAP_Analysis(pred_target = 'bliss')


    """
    moa_pairs = glob('./tissue/*')
    for m in moa_pairs:
        moa_pair = m.split('/')[-1]
        print(moa_pair)
        Downstream_SHAP_Analysis(pred_target = 'aoc', moa_pair = moa_pair)
        Downstream_SHAP_Analysis(pred_target = 'bliss', moa_pair = moa_pair)
    """

if __name__ == '__main__':
    main()





