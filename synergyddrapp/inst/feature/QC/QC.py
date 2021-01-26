#!usr/bin/python3
import pandas as pd
import numpy as np
import collections
import os, json, pickle
from tqdm import tqdm
from scipy import stats
from glob import glob
from itertools import combinations
from venn import venn
from matplotlib import pyplot 
from utils import *

# chemical structure required modules
from rdkit import Chem,DataStructs
from rdkit.Chem import MACCSkeys, AllChem
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit.Chem import rdmolops
from pubchempy import *
import openbabel
import pybel

def geneset_qc():
    """ QC on geneset information; Construct geneset feature matrix 
    
    The geneset matrix will be laid out like:
         |geneset1|geneset2|geneset3|...
    -------------------------------------
    gene1|   0    |    1   |   1    |
    -------------------------------------
    gene2|   0    |    0   |   1    |
    -------------------------------------
    gene3|   1    |    1   |   1    |
    -------------------------------------
    ...  |

    Params
    -------

    Yields
    ------
    """
    df = pd.read_csv('../data/merck_confidential_genesets_20200414.tsv', sep = '\t', header =0, index_col = 0)
    print("QC on geneset information ...")
    
    outpath = '../geneset_features/'
    os.makedirs(outpath, exist_ok = True)
    
    #print geneset total
    geneset_all = list(set(df['geneset_name'].to_list()))
    print("Total genesets:",len(geneset_all)) #9450
    f = open(outpath+'all_geneset.txt','w')
    f.write('\n'.join(geneset_all))
    f.close()

    #print gene total
    gene_all = sorted(list(set(df['gene_symbol'].to_list())))
    print("Total genes:",len(gene_all)) #18681
    df_mol = pd.read_csv('../molecular_features/exp_merck_confidential_features_20200221.tsv', sep = '\t', header =0, index_col = 0)
    gene = [x.split('_')[0] for x in df_mol.columns[3:]]
    overlap_gene = list(set(gene_all) &set(gene))

    # binray feature: 0 or 1
    features_df = {geneset:[] for geneset in geneset_all}
    for geneset in tqdm(geneset_all):
        gene_geneset = df.loc[df['geneset_name'] == geneset]['gene_symbol'].to_list()
        features_df[geneset] = [1 if x in gene_geneset else 0 for x in overlap_gene]
    features_df = pd.DataFrame.from_dict(features_df, orient = 'index', columns = overlap_gene)
    print(features_df)
    features_df.to_csv(outpath+"geneset_features.csv")

def molecular_qc():
    """
    Data types:

    1. Single gene molecular features: <GENESYMBOL>_snv/cnv/exp/lof;
    2. Geneset features: <GENESET>__coh_pat/lof_pat (geneset sep by '_');
    3. DDR features: *_ddr (78);
    
    Params
    ------
    Yields
    ------
    """

    # Read molecular raw data
    path = '../data/merck_confidential_features_20200221.tsv'
    data = pd.read_csv(path,sep = '\t',index_col=0 ,header = 0)
    
    # Create output path
    out_path = '../molecular_features/'
    print("Saving processed molecular features files to "+ out_path)
    os.makedirs(out_path, exist_ok = True)
    
    # Quantile normalize gene expression level
    print("Quantile normalize gene expression levels of molecular features ...")
    data = quantile_normalize_expression_levels(data)
    data = data.drop(columns = ['.metadata_cancer_subtype','.metadata_cancer_type'])
    data.to_csv(out_path+'quantiled_merck_confidential_features_20200221.tsv', index = False, sep = '\t')
    
    # Split datasets
    print("Split molecular data sets ...")
    columns = list(data.columns)
    snv = []
    cnv = []
    exp = []
    lof = []

    coh_pat = []
    lof_pat = []

    ddr = []
    uk = []
    while len(columns)>0:
        col = columns.pop()
        if col.endswith('_snv'):
            snv.append(col)
        elif col.endswith('_cnv'):
            cnv.append(col)
        elif col.endswith('_exp'):
            exp.append(col)
        elif col.endswith('_lof'):
            lof.append(col)
        elif col.endswith('__coh_pat'):
            coh_pat.append(col)
        elif col.endswith('__lof_pat'):
            lof_pat.append(col)
        elif col.endswith('_ddr'):
            ddr.append(col)
        else:
            uk.append(col)

    print("snv:", len(snv))
    print("cnv:", len(cnv))
    print("exp:", len(exp))
    print("lof:", len(lof))
    print("coh_pat:", len(coh_pat))
    print("lof_pat:", len(lof_pat))
    print("ddr:", len(ddr))
    print("uk columns:", len(uk), uk)
    celltype = ['.identifier_sample_name']
    all_mols = {"snv":set([x.split("_")[0] for x in snv]), "cnv":set([x.split("_")[0] for x in cnv]), "exp":set([x.split("_")[0] for x in exp]), "lof":set([x.split("_")[0] for x in lof])}
    #print(all_mols)
    plt = venn(all_mols)
    pyplot.savefig("mol.pdf")

    #exp: QC of gene expression level
    exp_cols = celltype+exp
    exp_data = data[exp_cols]
    #exp_data = exp_data.drop(columns = ['.metadata_cancer_subtype','.metadata_cancer_type'])
    exp_data.to_csv(out_path+'exp_merck_confidential_features_20200221.tsv', index = False, sep = '\t')

    #snv: QC of number of alleles affected by pathigenic mutation or indel
    snv_cols = celltype+snv
    snv_data = data[snv_cols]
    #snv_data = snv_data.drop(columns = ['.metadata_cancer_subtype','.metadata_cancer_type'])
    snv_data.to_csv(out_path+'snv_merck_confidential_features_20200221.tsv', index = False, sep = '\t')
    
    #cnv: QC of gene copy number variation
    cnv_cols = celltype+cnv
    cnv_data = data[cnv_cols]
    #cnv_data = cnv_data.drop(columns = ['.metadata_cancer_subtype','.metadata_cancer_type'])
    cnv_data.to_csv(out_path+'cnv_merck_confidential_features_20200221.tsv', index = False, sep = '\t')

    #lof: QC of predicted gene loss of function
    lof_cols = celltype+lof
    lof_data = data[lof_cols]
    #lof_data = lof_data.drop(columns = ['.metadata_cancer_subtype','.metadata_cancer_type'])
    lof_data.to_csv(out_path+'lof_merck_confidential_features_20200221.tsv', index = False, sep = '\t')

    # coh_pat
    coh_pat_cols = celltype+coh_pat
    coh_pat_data = data[coh_pat_cols]
    #coh_pat_data = coh_pat_data.drop(columns = ['.metadata_cancer_subtype','.metadata_cancer_type'])
    coh_pat_data.to_csv(out_path+'coh_pat_merck_confidential_features_20200221.tsv', index = False, sep = '\t')
    
    # lof_pat
    lof_pat_cols = celltype+lof_pat
    lof_pat_data = data[lof_pat_cols]
    #lof_pat_data = lof_pat_data.drop(columns = ['.metadata_cancer_subtype','.metadata_cancer_type'])
    lof_pat_data.to_csv(out_path+'lof_pat_merck_confidential_features_20200221.tsv', index = False, sep = '\t')
    
    # ddr
    ddr_cols = celltype + ddr
    ddr_data = data[ddr_cols]
    #ddr_data = ddr_data.drop(columns = ['.metadata_cancer_subtype','.metadata_cancer_type'])
    ddr_data.to_csv(out_path+'ddr_merck_confidential_features_20200221.tsv', index = False, sep = '\t')

def response_qc():
    """ QC of drug response data

    * publication 1 version: use synergy that must with monotherapt  response data
    
    Params
    ------
    Yields
    ------
    """
    def clean_missing_monotherapy(mono_features, synergy):
        """ Remove synergy responses with missing monotherapy responses data
        
        Params
        ------
        mono_features: Pandas dataframe
        synergy: Pandas dataframe
        
        Yields
        ------
        synergy: synergy after removed samples without monotherapy response features
        """
        idx = []
        for i, r in synergy.iterrows():
            if ((r['.identifier_sample_name']+'.'+r['.metadata_treatment_1'] in list(mono_features['cell_line.treatment'])) and (r['.identifier_sample_name']+'.'+r['.metadata_treatment_2'] in list(mono_features['cell_line.treatment']))):
                idx.append(i)
        synergy = synergy.loc[idx]
        return synergy
    # preprocess: remove used drugs
    path = '../data/'
    df = pd.read_csv(path+'merck_confidential_responses_training_20200221.tsv', sep = '\t', header = 0,index_col = 0)
    print(df)
    df_new = df.loc[(df['.metadata_treatment_1'] != 'MSC2531331')&(df['.metadata_treatment_1'] != 'MSC2691088_SRA737')&(df['.metadata_treatment_2'] != 'MSC2531331')&(df['.metadata_treatment_2'] != 'MSC2691088_SRA737')].copy().reset_index(drop = True)
    
    # replace moa subclasses with parent class
    responses =  df_new.replace('Cytostatic_Antimetabolite_Pyrimidineanalogue', 'Cytostatic_Antimetabolite').replace('Cytostatic_Alkylator', 'Cytostatic_Intercalator').replace('MEKi', 'CDKi').replace('CDK4i_6i', 'CDKi')
    # save processed responses
    responses.to_csv(path+'merck_confidential_responses_training_20200221.tsv', sep = '\t')

    # split monotherapy and synergy data
    responses, monotherapy, synergy = split_monotherapy_synergy(responses)
    
    # save cleaned/splited datasets
    responses.to_csv(path+'cleaned_merck_confidential_responses_training_20200221.tsv', sep = '\t', index = False,  na_rep = 'NA')
    monotherapy.to_csv(path+'monotherapy_responses.tsv', sep ='\t', index = False,  na_rep = 'NA')
    synergy.to_csv(path+'synergy_responses.tsv', sep ='\t', index = False,  na_rep = 'NA')

    # construct monotherapy response features
    out_path = '../monotherapy_features/'
    print("Saving processed monotherapy feature file to "+ out_path+"monotherapy_features.tsv")
    os.makedirs(out_path, exist_ok = True)
    mono_features = construct_monotherapy_response_features(monotherapy)
    mono_features.to_csv(out_path+'monotherapy_features.tsv', sep = '\t', index = False, na_rep = 'NA')

    # clean synergy(remove samples missing monotherapy responses)
    print("Before cleaning:", synergy.shape[0])
    synergy = clean_missing_monotherapy(mono_features, synergy)
    print("After cleaning:", synergy.shape[0])
    synergy.to_csv(path+'synergy_responses_with_monotherapy.tsv', sep ='\t', index = False,  na_rep = 'NA')

    # check #treatment in Monotherapy
    treat = list(set(monotherapy['.metadata_treatment']))
    print("Total treatment in monotherapy:", len(treat))
    #print(treat)

    # check #treatment in dual-therapy
    treat = list(set(list(synergy['.metadata_treatment_1'])+list(synergy['.metadata_treatment_2'])))
    print("Total treatment in dual-therapy:", len(treat))
    #print(treat)
    
    # check #cell line
    cl = list(set(monotherapy['.identifier_sample_name']))
    print("Total cell lines in monotherapy:", len(cl))
    cl = list(set(synergy['.identifier_sample_name']))
    print("Total cell lines in synergy:", len(cl))
    
    # check #tissue   
    tissue = list(set(monotherapy['.metadata_cancer_type']))
    print("Total cell lines in monotherapy:", len(tissue))
    tissue = list(set(synergy['.metadata_cancer_type']))
    print("Total cell lines in synergy:", len(tissue))

    # check overall preprodicibility
    check_reproducibility(responses, 'synergy', 'aoc')
    check_reproducibility(responses, 'synergy', 'bliss')
    # check synergy reprodicibility
    check_reproducibility(synergy, 'synergy', 'aoc')
    check_reproducibility(synergy, 'synergy', 'bliss')
    # check in-batch:
    """
    for batch in sorted(set(synergy['.identifier_batch'])):
        print(batch)
        synergy_new = synergy[synergy['.identifier_batch'] == batch]
        check_reproducibility(synergy_new, 'synergy', 'aoc')
        check_reproducibility(synergy_new, 'synergy', 'bliss')
    """
    # check monotherapy reprodicibility
    check_reproducibility(monotherapy, 'monotherapy', 'aoc')
    """
    # check in-batch
    for batch in sorted(set(monotherapy['.identifier_batch'])):
        print(batch)
        monotherapy_new = monotherapy[monotherapy['.identifier_batch'] == batch]
        check_reproducibility(monotherapy_new, 'monotherapy', 'aoc')
    """

    # check drug_pair tested in cell lines
    synergy['drug_pair'] = ['.'.join(sorted(list(r))) for _,r in synergy[['.metadata_treatment_1','.metadata_treatment_2']].iterrows()]
    cellline_summary = synergy.groupby('drug_pair')['.identifier_sample_name'].nunique().sort_values()
    mean = cellline_summary.mean()
    print("Averge # cellines tested for drug pairs:", mean)
    print(cellline_summary)

    # check combinations
    synergy['combination'] = synergy['drug_pair']+'.'+synergy['.identifier_sample_name']
    combination_summary = synergy.groupby('combination')['combination'].count().sort_values()
    print(combination_summary)

    # get all drugs tested
    drugs_summary = monotherapy[['.metadata_treatment', '.metadata_moa']].drop_duplicates().sort_values(by = '.metadata_treatment').rename(columns = {'.metadata_treatment':"drug_name",  '.metadata_moa':"mode-of-action"})
    drugs_summary.to_csv("all_drugs_summary.csv", index = False)


def chemical_structure_qc():
    """ QC of chemical structure information and feature preprocess

    Transform SMILE formatted chemical structure into 6 types of Fingerprints
    """

    def fp_to_feature(fp, max_bit):
        """
        """
        feature = [0 for i in range(max_bit)]
        for n in fp.bits:
            feature[n] = 1
        return feature

    df_map = pd.read_csv('../chemical_structure_features/merck_confidential_molecular_structures_20200420/smiles/smiles_and_ms_no.csv', header = 0)
    
    # make new dirs
    os.makedirs("../chemical_structure_features/MACCS_features/", exist_ok = True)
    os.makedirs("../chemical_structure_features/Morgan_features/", exist_ok = True)
    os.makedirs("../chemical_structure_features/RDK_features/", exist_ok = True)
    os.makedirs("../chemical_structure_features/FP2_features/", exist_ok = True)
    os.makedirs("../chemical_structure_features/FP3_features/", exist_ok = True)
    os.makedirs("../chemical_structure_features/FP4_features/", exist_ok = True)

    max_bit_fp2 = 1024
    max_bit_fp3 = 56
    max_bit_fp4 = 308

    all_features = []
    for _,r in df_map.iterrows():
        smile = r['treatment_smiles']
        ms = Chem.MolFromSmiles(smile)
        mol = pybel.readstring("smi", smile)
        
        # MACCS features (167*1)
        fp = MACCSkeys.GenMACCSKeys(ms)
        tmp = fp.ToBitString()
        feature_1 = list(map(int, tmp))
        np.savetxt('../chemical_structure_features/MACCS_features/'+r['treatment_name_orig'],np.array(feature_1))
        feature_1 = {'MACCS_'+str(idx+1):x for idx,x in enumerate(feature_1)}

        # Morgan Fingerprints (1024*1)
        fp = AllChem.GetMorganFingerprintAsBitVect(ms,2,nBits=1024)
        tmp = fp.ToBitString()
        feature_2 = list(map(int, tmp))
        np.savetxt('../chemical_structure_features/Morgan_features/'+r['treatment_name_orig'],np.array(feature_2))
        feature_2 = {'Morgan_'+str(idx+1):x for idx,x in enumerate(feature_2)}
        
        # FP2 (1024*1)
        fp = mol.calcfp('FP2')
        feature_3 = fp_to_feature(fp, max_bit_fp2)
        np.savetxt('../chemical_structure_features/FP2_features/'+r['treatment_name_orig'],np.array(feature_3))
        feature_3 = {'FP2_'+str(idx+1):x for idx,x in enumerate(feature_3)}

        # FP3 (56*1)
        fp = mol.calcfp('FP3')
        feature_4 = fp_to_feature(fp, max_bit_fp3)
        np.savetxt('../chemical_structure_features/FP3_features/'+r['treatment_name_orig'],np.array(feature_4))
        feature_4 = {'FP3_'+str(idx+1):x for idx,x in enumerate(feature_4)}

        # FP4 (308*1)
        fp = mol.calcfp('FP4')
        feature_5 = fp_to_feature(fp, max_bit_fp4)
        np.savetxt('../chemical_structure_features/FP4_features/'+r['treatment_name_orig'],np.array(feature_5))
        feature_5 = {'FP4_'+str(idx+1):x for idx,x in enumerate(feature_5)}

        # RDK Fingerprints (2048*1)
        fp = rdmolops.RDKFingerprint(ms)
        tmp = fp.ToBitString()
        feature_6 = list(map(int, tmp))
        np.savetxt('../chemical_structure_features/RDK_features/'+r['treatment_name_orig'],np.array(feature_6))
        feature_6 = {'RDK_'+str(idx+1):x for idx,x in enumerate(feature_6)}

        drug_name = {'treatment':r['treatment_name_orig']}
        the_features = {**drug_name,**feature_1,**feature_2,**feature_3,**feature_4,**feature_5,**feature_6}
        all_features.append(the_features)

    all_feature=pd.DataFrame.from_dict(all_features)
    all_feature.to_csv('../chemical_structure_features/all_chemical_structure.csv', index = False)


def cancer_dependency_qc():
    """ QC and feature construction of cancer dependency information
    """
    out_path = '../cancer_dependency_features/'
    print("Saving processed cancer dependency feature files to "+ out_path)
    os.makedirs(out_path, exist_ok = True)
    
    # CERES dataset
    dep = pd.read_csv(out_path+'integrated_Sanger_Broad_essentiality_matrices_20200402/CERES_FC.txt', sep = '\t', index_col = 0)
    moa_dep = build_moa_cencer_dependency(dep)
    p_out = out_path+'integrated_Sanger_Broad_essentiality_matrices_20200402/CERES_FC_dep.json'
    json.dump(moa_dep, open(p_out, 'w'))

    # CRISPRcleanR datase
    dep = pd.read_csv(out_path+'integrated_Sanger_Broad_essentiality_matrices_20200402/CRISPRcleanR_FC.txt', sep = '\t', index_col = 0)
    moa_dep = build_moa_cencer_dependency(dep)
    p_out = out_path+'integrated_Sanger_Broad_essentiality_matrices_20200402/CRISPRcleanR_FC_dep.json'
    json.dump(moa_dep, open(p_out, 'w'))

def main():
    #geneset_qc()
    #response_qc()
    molecular_qc()
    #chemical_structure_qc()
    #cancer_dependency_qc()

if __name__ == "__main__":
    main()
