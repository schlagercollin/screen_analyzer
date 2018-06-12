import numpy as np
import pandas as pd
from scipy import stats
import math

def ratio_test(bot5_comb, top5_comb, output_file):
    # Load in the mageck_input files
    # bot5_comb = pd.read_csv("/Users/collinschlager/Documents/Rohatgi_Lab/R Scripts/Kinome_hiOsmo_rep1_Bot5_mageck_input.txt", sep="\t")
    # top5_comb = pd.read_csv("/Users/collinschlager/Documents/Rohatgi_Lab/R Scripts/Kinome_hiOsmo_rep1_Top5_mageck_input.txt", sep="\t")
    bot5_comb = pd.read_csv(bot5_comb, sep="\t")
    top5_comb = pd.read_csv(top5_comb, sep="\t")

    # Remove control columns
    bot5_comb.drop(["Control1"], axis=1, inplace=True)
    top5_comb.drop(["Control1"], axis=1, inplace=True)

    # Rename counts to Rep1.Reads
    bot5_comb.rename(columns={'Sorted1': 'bot.Rep1.Reads'}, inplace=True)
    top5_comb.rename(columns={'Sorted1': 'top.Rep1.Reads'}, inplace=True)

    # Add 1 to reads
    bot5_comb['bot.Rep1.Reads'] = bot5_comb['bot.Rep1.Reads']+1
    top5_comb['top.Rep1.Reads'] = top5_comb['top.Rep1.Reads']+1

    def quantileNormalize(df_input):
        df = df_input.copy()
        #compute rank
        dic = {}
        for col in df:
            dic.update({col : sorted(df[col])})
        sorted_df = pd.DataFrame(dic)
        rank = sorted_df.mean(axis = 1).tolist()
        #sort
        for col in df:
            t = np.searchsorted(np.sort(df[col]), df[col])
            df[col] = [rank[i] for i in t]
        return df

    forQuant = pd.concat([bot5_comb['bot.Rep1.Reads'], top5_comb['top.Rep1.Reads']], axis=1)
    quantNorm = quantileNormalize(forQuant)
    # rank_mean = forQuant.stack().groupby(forQuant.rank(method='first').stack().astype(int)).mean()
    # quantNorm = forQuant.rank(method='min').stack().astype(int).map(rank_mean).unstack()


    bot5_comb['bot.Rep1.Quant'] = quantNorm.loc[:,'bot.Rep1.Reads']
    top5_comb['top.Rep1.Quant'] = quantNorm.loc[:,'top.Rep1.Reads']

    ratios = bot5_comb[['sgRNA', 'Gene', 'bot.Rep1.Quant']]
    ratios.rename(columns = {'bot.Rep1.Quant':'Bot5.Geom.Mean'}, inplace=True)
    ratios['Top5.Geom.Mean'] = top5_comb['top.Rep1.Quant']
    ratios['ratio'] = ratios['Bot5.Geom.Mean'] / ratios['Top5.Geom.Mean']
    ratios['log_ratio'] = ratios['ratio'].apply(np.log2)

    def geomMean(x):
        product = np.prod(x)
        numElem = len(x)
        result = product**(1/numElem)
        return result

    ratios['mean_abundance'] = forQuant.apply(geomMean, axis=1)
    ratios['log_MA'] = ratios['mean_abundance'].apply(np.log2)

    ratios.sort_values(by=['log_MA'], ascending=True, inplace=True)
    ratios['ZScore'] = np.zeros(len(ratios['log_MA']))

    bins=3

    guidesPerBin = math.ceil(len(ratios)/bins)
    ratios.reset_index(inplace=True)

    for i in range(bins):
        start = (i)*guidesPerBin
        end = start+guidesPerBin
        if end > len(ratios):
            end = len(ratios)
        guides = ratios[start:end]
        guides['ZScore'] = stats.zscore(guides.loc[:,'log_ratio'])
        # ratios['ZScore'][start:end] = guides.loc[:,'ZScore']
        ratios.loc[start:end,'ZScore'] = guides.loc[:,'ZScore']

    ratios.sort_values(by=['ZScore'], ascending=False, inplace=True)
    ratios.to_csv(output_file)
    return ratios

if __name__ == "__main__":
    bot5_test_file = "/Users/collinschlager/Documents/Rohatgi_Lab/R Scripts/Kinome_hiOsmo_rep1_Bot5_mageck_input.txt"
    top5_test_file - "/Users/collinschlager/Documents/Rohatgi_Lab/R Scripts/Kinome_hiOsmo_rep1_Top5_mageck_input.txt"
    output = "/Users/collinschlager/Documents/Rohatgi_Lab/R Scripts/Kinome_hiOsmo_rep1_Top5_Bot5_ratio_test.csv"
    ratios_df = ratio_test(bot5_test_file, top5_test_file, output)
