import pandas as pd
import numpy as np

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

if __name__ == "__main__":
    df = pd.DataFrame([[5,4,3],[2,1,4],[3,4,6],[4,2,8]])
    result = quantileNormalize(df)
    print(result)
