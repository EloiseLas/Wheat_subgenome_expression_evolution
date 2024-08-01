import pandas as pd
import numpy as np
import math

def category(scoreAB,scoreBD,scoreAD):
    cat="Unknown"
    scores=[scoreAB,scoreBD,scoreAD]
    if sum(math.isnan(x) for x in scores)>0:
        cat="Unknown"
    elif scoreAB>0.9 and scoreBD>0.9 and scoreAD>0.9:
        cat="Same_GO"
    elif scoreAB<0.9 and scoreBD>0.9 and scoreAD<0.9:
        cat="A_diff"
    elif scoreAB<0.9 and scoreBD<0.9 and scoreAD>0.9:
        cat="B_diff"
    elif scoreAB>0.9 and scoreBD<0.9 and scoreAD<0.9:
        cat="D_diff"
    elif scoreAB<0.9 and scoreBD<0.9 and scoreAD<0.9:
        cat="All_diff"
    return cat

homeo=pd.read_csv('Output/homeo_annotation_GO_cat2.csv', sep=",", header=0)
homeo["category"]=homeo.apply(lambda row : category(row["scores_AB"],row["scores_BD"],row["scores_AD"]), axis = 1)
print(homeo)

homeo.to_csv('Output/homeo_annotation_GO_category_precise2.csv')

counts=homeo.groupby(["category"]).count()
print(counts)
counts.to_csv('Output/annotation_category_counts2.csv')
