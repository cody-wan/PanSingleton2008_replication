import numpy as np 
import pandas as pd



if __name__ == '__main__':
    """ test code """
    # check own implementation with authors'
    CDS_df = pd.read_csv("OUTC0_CDS_test.csv")
    OUTC11BRAZIL = pd.read_csv("OUTC11BRAZIL", delim_whitespace=True, header=None, dtype=float).iloc[:-1] # exclude last row (non CDS data)
    # total element-wise abs differences between two tables
    print(f"{np.abs((CDS_df.values - OUTC11BRAZIL.values)).sum():.5f}")