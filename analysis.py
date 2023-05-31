import numpy as np # scientific computing
import pandas as pd # data loading and processing
import os # os operations
import matplotlib.pyplot as plt # for generating figures
import math
import matplotlib.dates as mdates
import seaborn as sns # for generating visualizations, better support with pandas than matplotlib
from scipy import stats
from sklearn.impute import SimpleImputer

# function to dynamically generate HCCDC file name
def construct_hccdb_filename(n):
    expression = "./data/HCCDB/HCCDB" + n + "_mRNA_level3.txt"
    pheno = "./data/HCCDB/HCCDB" + n  + ".sample.txt"
    return expression,pheno

# function to get HCCDB data
def get_hccdb_data(n):
    df = pd.read_csv(n, index_col = 1, sep = "\t").drop(["Entrez_ID"], axis=1) # gene x patient
    return df

# function to fetch all gene signatures from source file
def get_gene_signature_file():
    df = pd.read_csv("./data/oxstress genes.csv", index_col=None, header= 0)
    return df

# function to get relevant signatures from gene signature dataframe
def get_xy_set(gene_set, xvar=None,yvar=None):
    if xvar != "RRM2B"and xvar != None:
        x_set = gene_set[xvar].dropna().tolist()
    else:
        x_set = ["RRM2B"]

    if yvar != "G6PD" and yvar != None:
        y_set = gene_set[yvar].dropna().tolist()
    else:
        y_set = ["G6PD"]
    targets = list(set(x_set + y_set))
    return x_set, y_set, targets

# extract the specific cancer type from the data
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

# returns 60 x 40 figure
def generate_subplots(title, x, y):
    tot = len(y)
    cols = 8
    rows = tot//8
    if tot % cols != 0:
        rows += 1
    fig, axs = plt.subplots(rows, cols, figsize=(60, 40), sharey=True)
    plt.subplots_adjust(hspace=0.6)
    fig.suptitle(title,fontsize = 40)
    return fig, axs
    
# get raw data from HCCDB and TCGA files
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

# function to subset relevant genes, impute missing data, normalise data, and get composite scores
def process_data(raw_data, targets, x_var_names, y_var_names, pheno_filtered=None, outlier_corrected = False):
    print("entering process data")
    # df is inputted as gene x patient 
    df = raw_data.T # patients x genes
    print("transposed")

    # subset to get relevant genes
    df_filtered = df[targets]
    
    # impute missing values with mean
    imp_mean = SimpleImputer(missing_values=np.nan, strategy='mean')
    output = imp_mean.fit_transform(df_filtered)
    df_filtered = pd.DataFrame(output, columns = imp_mean.get_feature_names_out().tolist(), index = df_filtered.index)

    # normalise data
    df_filtered = df_filtered.astype(np.float64)
    df_filtered=(df_filtered-df_filtered.mean())/(df_filtered.std()+1)

    # handle dropped columns due to imputation
    x_var_names = list(set(x_var_names).intersection(set(df_filtered.columns)))
    y_var_names = list(set(y_var_names).intersection(set(df_filtered.columns)))

    # average the genes in each gene set to get a composite score
    if len(x_var_names) != 0:
        x_var_gene_set = df_filtered.loc[:,x_var_names]
        x_var_gene_set.loc[:,"x_composite_score"] = x_var_gene_set.mean(axis = 1)
        df_filtered = pd.concat([df_filtered, x_var_gene_set[["x_composite_score"]]], axis = 1) 

    else:
        raise KeyError("no x var gene set")

    if len(y_var_names) != 0:
        y_var_gene_set = df_filtered.loc[:,y_var_names]
        y_var_gene_set.loc[:,"y_composite_score"]= y_var_gene_set.mean(axis = 1)
        df_filtered = pd.concat([df_filtered, y_var_gene_set[["y_composite_score"] ]], axis = 1) # patients x genes 
    else:
        raise KeyError("no y var gene set")

    # outlier correction
    if outlier_corrected == True:
        q25, q75 = df_filtered["x_composite_score"].describe()["25%"], df_filtered["x_composite_score"].describe()["75%"]
        iqr = q75 - q25
        cut_off = iqr * 1.5
        lower, upper = q25 - cut_off, q75 + cut_off
        df_filtered = df_filtered.loc[(df_filtered["x_composite_score"] > lower) & (df_filtered["x_composite_score"] < upper), :] 
    
    return df_filtered

# function to  calculate r and p value, and plot the data
def analyse(data, fig, db, ax, fn, x_label, y_label, x_target = "x_composite_score", y_target = "y_composite_score", dataset_screen = False, plotter = True):
    # get x and y data
    y, x = data[y_target].to_numpy(), data[x_target].to_numpy()

    # get r value
    r = np.corrcoef(x, y)[0, 1]

    # find p-value
    n = data.shape[0]
    t = (r - math.sqrt(n - 2)) / math.sqrt(1 - (r ** 2))
    p = stats.t.sf(abs(t), df=n-2)*2  # Survival function (also defined as 1 - cdf, but sf is sometimes more accurate).
    # our t-test is two-sided, we need to multiply this value by 2
    
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
    if plotter == True:
        sns.set_style("ticks")
        sns.scatterplot(data=data, x=x_target, y=y_target, ax= ax)
        
        a, b = np.polyfit(x, y, 1)
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

def impute_nan(df):
    # df is inputted as gene x patient 
    print("imputing data")
    print("transpose")
    df = df.T # patients x genes

    print("impute")
    imp_mean = SimpleImputer(missing_values=np.nan, strategy='mean')
    output = imp_mean.fit_transform(df)
    df = pd.DataFrame(output, columns = df.columns, index = df.index)
    
    df = df.astype(np.float64)
    print("done imputing")
    return df

def impute_nan_general(df):
    # df is inputted as gene x patient 
    print("imputing data")
    print("transpose")

    print("impute")
    imp_mean = SimpleImputer(missing_values=np.nan, strategy='mean')
    output = imp_mean.fit_transform(df)
    df = pd.DataFrame(output, columns = df.columns, index = df.index)
    
    print("done imputing")
    return df

def get_gene_names(filename,col=None):
    file = pd.read_csv(filename, index_col=None, header= 0).T
    names = file[col].dropna().tolist()
    return names
