import pandas as pd
import math

def category(corAB,corBD,corAD):
    cat="Unknown"
    cors=[corAB,corBD,corAD]
    thresh=0.1
    signif_count=sum(1 for i in cors if i>thresh)
    if signif_count>=2:
        cat="S"
    elif corAB<=thresh and corBD>thresh and corAD<=thresh:
        cat="A"
    elif corAB<=thresh and corBD<=thresh and corAD>thresh:
        cat="B"
    elif corAB>thresh and corBD<=thresh and corAD<=thresh:
        cat="D"
    elif corAB<=thresh and corBD<=thresh and corAD<=thresh:
        cat="F"
    return cat

homeo_cor=pd.read_csv('Output/homeo_expression_profiles_cor_new.csv', sep=",", header=0)
homeo_cor=homeo_cor.dropna()
homeo_cor["category"]=homeo_cor.apply(lambda row : category(row["cor_AB"],row["cor_BD"],row["cor_AD"]), axis = 1)
print(homeo_cor)

homeo_cor_clean=pd.concat([homeo_cor.iloc[:,1:8],homeo_cor.iloc[:,12]],axis=1)
print(homeo_cor_clean)

homeo_cor.to_csv('Output/homeo_cor_category_new.csv')

counts=homeo_cor.groupby(["category"]).count()
print(counts)
counts.to_csv('Output/cor_category_counts_new.csv')