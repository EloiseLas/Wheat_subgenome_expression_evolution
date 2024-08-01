import pandas as pd
import math

def category(mrAB,mrBD,mrAD,corAB,corBD,corAD):
    cat="Unknown"
    scores=[(mrAB,corAB),(mrBD,corBD),(mrAD,corAD)]
    thresh_mr=1
    thresh_cor=0.5
    signif_count=sum(1 for i in scores if i[0]<thresh_mr and i[1]>thresh_cor)
    if signif_count>=2:
        cat="S"
    elif (mrAB>=thresh_mr or corAB<=thresh_cor) and (mrBD<thresh_mr and corBD>thresh_cor) and (mrAD>=thresh_mr or corAD<=thresh_cor):
        cat="A"
    elif (mrAB>=thresh_mr or corAB<=thresh_cor) and (mrBD>=thresh_mr or corBD<=thresh_cor) and (mrAD<thresh_mr and corAD>thresh_cor):
        cat="B"
    elif (mrAB<thresh_mr and corAB>thresh_cor) and (mrBD>=thresh_mr or corBD<=thresh_cor) and (mrAD>=thresh_mr or corAD<=thresh_cor):
        cat="D"
    elif (mrAB>=thresh_mr or corAB<=thresh_cor) and (mrBD>=thresh_mr or corBD<=thresh_cor) and (mrAD>=thresh_mr or corAD<=thresh_cor):
        cat="F"
    return cat

homeo=pd.read_csv("Output/homeo_expression_profiles_cor_mutualrank.csv", sep=",", header=0)
homeo=homeo.dropna()
homeo["category"]=homeo.apply(lambda row : category(row["AB_log_mr"],row["BD_log_mr"],row["AD_log_mr"], row["cor_AB"],row["cor_BD"],row["cor_AD"]), axis = 1)
print(homeo)

# homeo_clean=pd.concat([homeo.iloc[:,1:8],homeo.iloc[:,12]],axis=1)
# print(homeo_mr_clean)

homeo.to_csv('Output/homeo_profile_both_category.csv')

counts=homeo.groupby(["category"]).count()
print(counts)
counts.to_csv('Output/profile_both_category_counts.csv')