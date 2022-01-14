import pandas as pd
import prince
import string
from decimal import Decimal
import numpy as np

## Raw data
data = pd.read_csv("H_DAD_w1a3.csv")
data_filled = pd.read_csv("FAMD_filled.csv")
col_type = pd.read_excel("All_Tables.xlsx", sheet_name="Sheet1")
col_type = col_type[["Variable","Type"]]
categ = col_type.loc[col_type["Type"] == "Categ", "Variable"].tolist()
categ = list(map(lambda x: x.lower(), categ))
data_filled = data_filled.loc[:,data_filled.columns!="Unnamed: 0"]
data_filled = data_filled.loc[:,data_filled.columns!="prim_key"]
data_filled["id"] = data_filled.index + 1

data_filled[categ] = data_filled[categ].astype(str)

famd = prince.FAMD(n_components=3, copy= True, check_input = True, engine = "fbpca", random_state = 42)
famd = famd.fit(data_filled.drop("r1cdr_final",axis="columns"))
res = famd.row_coordinates(data_filled)
print(res)
res.to_csv("FAMD_Coordinates.csv",header=True,index=True)

res_col = famd.column_correlations(data_filled)
print(res_col)
res_col.to_csv("FAMD_column_correlations.csv",header=True,index=True)

## Processed data
data = pd.read_csv("processed_filled.csv")

catdata = pd.read_csv("processed_categorical_data.csv")
col_type = pd.read_excel("All_Tables.xlsx", sheet_name="Sheet1")
col_type = col_type[["Variable","Type"]]
catdata = catdata.loc[:,catdata.columns!="Unnamed: 0"]
categ = catdata.columns.tolist()
data = data.loc[:,data.columns!="Unnamed: 0"]

data[categ] = data[categ].astype(str)


famd = prince.FAMD(n_components=3, copy= True, check_input = True, engine = "fbpca", random_state = 42)
famd = famd.fit(data)
res = famd.row_coordinates(data)
print(res)
res.to_csv("FAMD_Processed_Coordinates.csv",header=True,index=True)

res_col = famd.column_correlations(data)
print(res_col)
res_col.to_csv("FAMD_Processed_column_correlations.csv",header=True,index=True)

