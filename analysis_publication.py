import numpy as np # scientific computing
import pandas as pd # data loading and processing
import os # os operations
import matplotlib.pyplot as plt # for generating figures
import math
import matplotlib.dates as mdates
import seaborn as sns # for generating visualizations, better support with pandas than matplotlib
from scipy import stats
from sklearn.impute import SimpleImputer

def get_gene_signature_file():
    df = pd.read_csv("./data/oxstress genes.csv", index_col=None, header= 0)
    return df

def extract_rows_by_type(data, hccdb=None, db='PANCAN'):
    if db.startswith("HCCDB"):
        df = hccdb.T
        df = df[df["ptype"] == db]
        df = df.T
        df.drop(["ptype"], inplace = True)
    elif db == "PANCAN":
        df = data
        df = df.T
        df.drop(["ptype","sample_type_id", "sample_type", "_primary_disease"], inplace = True)
    elif db == "Aggregated": 
        ls = ['STAD', 'HNSC', 'SARC', 'UCS','LUSC', 'BRCA']
        df = data[data["ptype"].isin(ls)]
        df = df.T
        df.drop(["ptype","sample_type_id", "sample_type", "_primary_disease"], inplace = True)	
    else:
        df = data[data["ptype"] == db]
        df = df.T
        df.drop(["ptype","sample_type_id", "sample_type", "_primary_disease"], inplace = True)	
    return df

def get_raw_data():

    # script to consolidate all HCCDB data into one dataframe
    hccdb_names = ["1", "3", "4",  "8", "9", "11", "12", "13", "14", "16", "17", "18"]
    hccdb = pd.DataFrame()

    for i in range(len(hccdb_names)):
        n1, n2 = construct_hccdb_filename(hccdb_names[i])
        hccdb_temp = get_hccdb_data(n1)
        hccdb_temp = hccdb_temp.loc[~hccdb_temp.index.duplicated(),:].copy()
        hccdb_temp.loc["ptype",:] = "HCCDB-" + hccdb_names[i]
        hccdb = pd.concat([hccdb, hccdb_temp], axis = 1) # patients x genes

    # load pancan data
    tcga = pd.read_csv("./data/EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp (1).xena", index_col = 0, sep = "\t") # gene x patient
    pheno = pd.read_csv("./data/TCGA_phenotype_denseDataOnlyDownload (1).tsv", index_col = 0, sep = "\t") # patient x phenotype

    # attach cancer type to each patient
    data = tcga.T
    data = pd.concat([data, pheno], axis = 1, join = "inner") # patients x genes
        
    print(data.shape)
    print(tcga.T.shape)

    # attach abbeviations for each cancer type
    ls = data["_primary_disease"].unique().tolist()

    conditions = [
        data['_primary_disease'] == 'adrenocortical cancer',
        data['_primary_disease'] == 'bladder urothelial carcinoma',
        data['_primary_disease'] == 'breast invasive carcinoma',
        data['_primary_disease'] == 'cervical & endocervical cancer',
        data['_primary_disease'] == 'cholangiocarcinoma', 
        data['_primary_disease'] == 'colon adenocarcinoma',
        data['_primary_disease'] == 'diffuse large B-cell lymphoma',
        data['_primary_disease'] == 'esophageal carcinoma',
        data['_primary_disease'] == 'glioblastoma multiforme',
        data['_primary_disease'] == 'head & neck squamous cell carcinoma',
        data['_primary_disease'] == 'kidney chromophobe',
        data['_primary_disease'] == 'kidney clear cell carcinoma',
        data['_primary_disease'] == 'kidney papillary cell carcinoma',
        data['_primary_disease'] == 'acute myeloid leukemia',
        data['_primary_disease'] == 'brain lower grade glioma',
        data['_primary_disease'] == 'liver hepatocellular carcinoma',
        data['_primary_disease'] == 'lung adenocarcinoma',
        data['_primary_disease'] == 'lung squamous cell carcinoma',
        data['_primary_disease'] == 'mesothelioma',
        data['_primary_disease'] == 'ovarian serous cystadenocarcinoma',
        data['_primary_disease'] == 'pancreatic adenocarcinoma',
        data['_primary_disease'] == 'pheochromocytoma & paraganglioma',
        data['_primary_disease'] == 'prostate adenocarcinoma',
        data['_primary_disease'] == 'rectum adenocarcinoma',
        data['_primary_disease'] == 'sarcoma',
        data['_primary_disease'] == 'skin cutaneous melanoma',
        data['_primary_disease'] == 'stomach adenocarcinoma',
        data['_primary_disease'] == 'testicular germ cell tumor',
        data['_primary_disease'] == 'thyroid carcinoma',
        data['_primary_disease'] == 'thymoma',
        data['_primary_disease'] == 'uterine corpus endometrioid carcinoma',
        data['_primary_disease'] == 'uterine carcinosarcoma',
        data['_primary_disease'] == 'uveal melanoma'    
    ]

    choices = ["ACC",
            "BLCA",
            "BRCA",
            "CESC",
            "CHOL",
            "COAD",
            "DBLC",
            "ESCA",
            "GBM",
            "HNSC",
            "KICH",
            "KIRC",
            "KIRP",
            "LAML",
            "LGG",
            "LIHC",
            "LUAD",
            "LUSC",
            "MESO",
            "OV",
            "PAAD",
            "PCPG",
            "PRAD",
            "READ",
            "SARC",
            "SKCM",
            "STAD",
            "TGCT",
            "THCA",
            "THYM",
            "UCEC",
            "UCS",
            "UVM"
            ]

    data["ptype"] = np.select(conditions, choices, default = "null")
    
    return data, hccdb

def get_gene_names(filename,col=None):
    file = pd.read_csv(filename, index_col=None, header= 0).T
    names = file[col].dropna().tolist()
    return names

def construct_filename(c, db):
    if db =="xena":
        n1 = "./data/" + "TCGA." + c + ".sampleMap_HiSeqV2"
        n2 = "./data/" + "TCGA." + c + ".sampleMap_"+ c +"_clinicalMatrix"
    elif db == "cbio":
        n1 = "./data/" + c+"_data_mrna_seq_v2_rsem.txt"
        n2 = "./data/" + c + "_data_clinical_sample.txt"
    else:
        print("db must be either xena or cbio")
    return n1,n2

def construct_hccdb_filename(n):
    n1 = "./data/HCCDB/HCCDB" + n + "_mRNA_level3.txt"
    n2 = "./data/HCCDB/HCCDB" + n  + ".sample.txt"
    return n1,n2

def process_data(df, targets, x_var_names = None, y_var_names = None, pheno_filtered=None, outlier_corrected = False):

    # df is inputted as gene x patient 
    df = df.T # patients x genes

    # subset to get relevant genes
    df_filtered = df[targets]
    
    # impute missing values with mean
    imp_mean = SimpleImputer(missing_values=np.nan, strategy='mean')
    output = imp_mean.fit_transform(df_filtered)
    df_filtered = pd.DataFrame(output, columns = imp_mean.get_feature_names_out().tolist(), index = df_filtered.index)

    # scale numerical data
    df_filtered = df_filtered.astype(np.float64)

    # for each sequenced gene were rescaled to set the median equal to 1
    df_filtered=(df_filtered-df_filtered.median())/(df_filtered.std()+1)
    data = df_filtered  

    print(list(set(x_var_names).symmetric_difference(set(x_var_names))))
    print(list(set(y_var_names).symmetric_difference(set(y_var_names))))

    x_var_names = list(set(x_var_names).intersection(set(data.columns)))
    y_var_names = list(set(y_var_names).intersection(set(data.columns)))

    if x_var_names != None:
        x_var_gene_set = data[x_var_names]
        x_var_gene_set["x_composite_score"] = x_var_gene_set.mean(axis = 1)
        data = pd.concat([data, x_var_gene_set[["x_composite_score"]]], axis = 1) 
        if outlier_corrected == True:
            iqr = data["x_composite_score"].describe()
            data = data.loc[(data["x_composite_score"] > iqr['min']) & (data["x_composite_score"] < iqr['max']), :] 

    if y_var_names != None:
        # take only nrf2 target genes
        y_var_gene_set = data[y_var_names]
        y_var_gene_set["y_composite_score"]= y_var_gene_set.mean(axis = 1)
        data = pd.concat([data, y_var_gene_set[["y_composite_score"] ]], axis = 1) # patients x genes 
    
    return data


def analyse(data, fig, db, ax, fn, x_label, y_label, x_target = "x_composite_score", y_target = "y_composite_score", dataset_screen = False):
    
    #find line of best fit
    y, x = data[y_target].to_numpy(), data[x_target].to_numpy()
    a, b = np.polyfit(x, y, 1)

    # bin the patients into quartiles based on G6PD expression
    iqr = data[x_target].T.describe()
    data["RRM2B levels"] = pd.cut(data["RRM2B"],
                    bins=[ iqr["min"], iqr["25%"], iqr["75%"], iqr["max"]],
                    labels=["Bottom 25%", "-", "Top 25%"])

    # get r sq val
    r = np.corrcoef(x, y)[0, 1]

    #find p-value
    n = data.shape[0]
    t = (r-math.sqrt(n-2))/math.sqrt(1-(r**2))
    p = stats.t.sf(abs(t), df=n)*2 
    if p < 0.0001:
        pval = "< 0.0001"
    elif p <0.001:
        pval = "< 0.001"
    elif p<0.01:
        pval = "< 0.01"
    elif p<0.05:
        pval = "< 0.05"
    else:
        pval = "= N.S."

    # plot the data
    # scatter plot for RRM2B against NRF2 activity
    sns.set_style("ticks")
    sns.scatterplot(data=data, x=x_target, y=y_target, ax= ax)
    ax.plot(x, a*x+b, color="black")
    ax.set_ylabel(y_label,fontsize = 28)
    ax.set_xlabel(x_label + " \n (r = " + str(round(r, 4)) + "," + " p " + pval +")",fontsize = 25)
    name = db + " (n = " + str(data.shape[0]) + ")"
    if dataset_screen == True:
        name = db + " (n = " + str(data.shape[0]) + ")"
    else:
        name = db
    ax.set_title(name, fontsize = 30)
    ax.tick_params(axis='both', which='major', labelsize=25)
    sns.despine()
    plt.show()

    # save the figure 
    fig.savefig(fn)
    return r, p

def get_targets_present(data, targets):
    idx = data.index.to_list()
    # print(idx)
    # print(targets)
    targets_present = list(set(idx).intersection(set(targets)))
    # print(targets_present)
    return targets_present

def get_xena_data(n1):
    df = pd.read_csv(n1, index_col = 0, sep = "\t") # gene x patient
    return df

def get_cbio_data(n1):
    df = pd.read_csv(n1, index_col = 0, sep = "\t").drop(["Entrez_Gene_Id"], axis=1) # gene x patient
    return df

def get_hccdb_data(n1):
    df = pd.read_csv(n1, index_col = 1, sep = "\t").drop(["Entrez_ID"], axis=1) # gene x patient
    return df

def get_xena_pheno(n2):
    pheno = pd.read_csv(n2, index_col=0, sep = "\t")
    pheno = pheno[["sample_type"]]
    pheno_filtered = pheno.dropna()
    return pheno_filtered

def get_cbio_pheno(n2):
    pheno = pd.read_csv(n2, index_col=1,header = 4, sep = "\t")
    pheno = pheno[["SAMPLE_TYPE"]]
    pheno_filtered = pheno.dropna()
    return pheno_filtered

def get_hccdb_pheno(n2):
    pheno = pd.read_csv(n2, index_col=0, sep = "\t").T
    pheno = pheno[["TYPE"]]
    pheno_filtered = pheno.dropna()
    return pheno_filtered

