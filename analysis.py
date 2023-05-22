import numpy as np # scientific computing
import pandas as pd # data loading and processing
import os # os operations
import matplotlib.pyplot as plt # for generating figures
import math
import matplotlib.dates as mdates
import seaborn as sns # for generating visualizations, better support with pandas than matplotlib
from scipy import stats
from sklearn.impute import SimpleImputer

def get_data(data, hccdb, db):
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

    if x_var_names != None and y_var_names!=None:
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


def analyse(data, fig, db, ax, fn, x_label, y_label, x_target = "x_composite_score", y_target = "y_composite_score"):
    
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
        pval = "<0.0001"
    elif p <0.001:
        pval = "<0.001"
    elif p<0.01:
        pval = "<0.01"
    elif p<0.05:
        pval = "<0.05"
    else:
        pval = "N.S."

    # plot the data
    # scatter plot for RRM2B against NRF2 activity
    sns.set_style("whitegrid")
    sns.set()
    sns.scatterplot(data=data, x=x_target, y=y_target, ax= ax)
    ax.plot(x, a*x+b, color="black")
    ax.set_ylabel(y_label,fontsize = 28)
    ax.set_xlabel(x_label + " \n (r = " + str(round(r, 4)) + "," + " p = " + pval +")",fontsize = 25)
    name = db + " (n = " + str(data.shape[0]) + ")"
    ax.set_title(name, fontsize = 30)
    ax.tick_params(axis='both', which='major', labelsize=25)
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

