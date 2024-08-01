import pandas as pd
import math

def category(mrAB,mrBD,mrAD):
    cat="Unknown"
    mrs=[mrAB,mrBD,mrAD]
    thresh=2
    signif_count=sum(1 for i in mrs if i<thresh)
    if signif_count>=2:
        cat="S"
    elif mrAB>=thresh and mrBD<thresh and mrAD>=thresh:
        cat="A"
    elif mrAB>=thresh and mrBD>=thresh and mrAD<thresh:
        cat="B"
    elif mrAB<thresh and mrBD>=thresh and mrAD>=thresh:
        cat="D"
    elif mrAB>=thresh and mrBD>=thresh and mrAD>=thresh:
        cat="F"
    return cat

homeo_mr=pd.read_csv('Output/homeo_expression_profiles_log_mutualrank_new.csv', sep=",", header=0)
homeo_mr=homeo_mr.dropna()
homeo_mr["category"]=homeo_mr.apply(lambda row : category(row["AB_log_mr"],row["BD_log_mr"],row["AD_log_mr"]), axis = 1)
print(homeo_mr)

homeo_mr_clean=pd.concat([homeo_mr.iloc[:,1:8],homeo_mr.iloc[:,12]],axis=1)
print(homeo_mr_clean)

homeo_mr.to_csv('Output/homeo_logmr_category_new.csv')

counts=homeo_mr.groupby(["category"]).count()
print(counts)
counts.to_csv('Output/logmr_category_counts_new.csv')