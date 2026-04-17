import numpy as np
import pandas as pd

def resized(array:np.array,new_size:int,defaultvalue:any=0):
    if new_size<=len(array): 
        return array[:new_size]
    tmp = np.full(new_size,defaultvalue,array.dtype)
    tmp[:len(array)] = array
    return tmp

def read_tsv(dir:str)->pd.DataFrame:
    with open(dir,'r') as file:
        firstline = True
        dic_df = {}
        for line in file.readlines():
            line_list_raw = line.split(' ') # Separates the values and deletes the new line char
            line_list = [x for x in line_list_raw if x not in ('','\n')] # Removes empty strings nad new line for safe measure
            if firstline:
                firstline = False
                cols = line_list.copy()
                for col in cols: dic_df[col] = []
            else:
                assert len(cols) == len(line_list), (len(dic_df[cols[0]]),line_list)
                for i,entry in enumerate(line_list): dic_df[cols[i]].append(entry)
    return pd.DataFrame(dic_df)