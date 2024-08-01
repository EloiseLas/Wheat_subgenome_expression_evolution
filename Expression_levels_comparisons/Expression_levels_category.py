import pandas as pd

def category(AB_pv,AD_pv,BD_pv,means_A,means_B,means_D):
    cat="Unknown"
    signif_count=sum(1 for i in [AB_pv,AD_pv,BD_pv] if i<0.05)
    if AB_pv<0.05 and AD_pv<0.05 and BD_pv>=0.05:
        cat="A"
        if means_A>(means_B+means_D)/2:
            cat+="+"
        else:
            cat+="-"
    elif signif_count==3:
        cat="F"
        if means_A>means_B and means_B>means_D:
            cat+="_ABD"
        elif means_D>means_A and means_A>means_B:
            cat+="_DAB"
        elif means_B>means_A and means_A>means_D:
            cat+="_BAD"
        elif means_A>means_D and means_D>means_B:
            cat+="_ADB"
        elif means_D>means_B and means_B>means_A:
            cat+="_DBA"
        elif means_B>means_D and means_D>means_A:
            cat+="_BDA"
    elif AB_pv<0.05 and AD_pv>=0.05 and BD_pv<0.05:
        cat="B"
        if means_B>(means_A+means_D)/2:
            cat+="+"
        else:
            cat+="-"
    elif AB_pv>=0.05 and AD_pv<0.05 and BD_pv<0.05:
        cat="D"
        if means_D>(means_A+means_B)/2:
            cat+="+"
        else:
            cat+="-"
    elif signif_count==0 or signif_count==1:
        cat="S"
    return cat

homeo_pvalue=pd.read_csv('Output/homeologs_pvalue_new.csv', sep=",", header=0)
homeo_pvalue["category"]=homeo_pvalue.apply(lambda row : category(row["AB_pv"],row["AD_pv"],row["BD_pv"],row["means_A"],row["means_B"],row["means_D"]), axis = 1)
print(homeo_pvalue)

categories=["A+","A-","F_ABD","F_DAB","F_BAD","F_ADB","F_DBA","F_BDA","B+","B-","D+","D-","S"]

for c in categories:
    res1 = homeo_pvalue[homeo_pvalue["category"]==c][["A","B","D","ID","length_diff"]]
    res2 = homeo_pvalue[homeo_pvalue["category"]==c][["EGI_LOC_A","EGI_LOC_B","EGI_LOC_D","ID","length_diff"]]
    res1.to_csv('Output/category' + c + '.csv')
    res2.to_csv('Output/category' + c + '_EGI.csv')

counts=homeo_pvalue.groupby(["category"]).count()
print(counts)
counts.to_csv('Output/homeo_pvalue_category_counts_new.csv')
counts=homeo_pvalue[homeo_pvalue["length_diff"]==True].groupby(["category"]).count()
counts.to_csv('Output/homeo_pvalue_category_counts_lengthdiff_new.csv')
counts=homeo_pvalue[homeo_pvalue["length_diff"]==False].groupby(["category"]).count()
counts.to_csv('Output/homeo_pvalue_category_counts_nolengthdiff_new.csv')

homeo_cat=homeo_pvalue[["ID","A","B","D","EGI_LOC_A","EGI_LOC_B","EGI_LOC_D","category","length_diff"]]

homeo_cat.to_csv('Output/homeo_pvalue_category_wilcox_new.csv')
homeo_pvalue.to_csv('Output/homeo_pvalue_category_detail_wilcox_new.csv')