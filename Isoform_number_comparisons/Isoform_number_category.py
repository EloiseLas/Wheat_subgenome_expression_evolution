import pandas as pd

def category(isoA,isoB,isoD):
    cat="Unknown"
    diffAB=abs(isoA-isoB)
    diffBD=abs(isoB-isoD)
    diffAD=abs(isoA-isoD)
    diffs=[diffAB,diffBD,diffAD]
    count=sum(1 for i in diffs if i==0)
    if count==3:
        cat="S"
    elif diffAB>0 and diffBD==0 and diffAD>0:
        cat="A"
        if isoA>(isoB+isoD)/2:
            cat+="+"
        else:
            cat+="-"
    elif diffAB>0 and diffBD>0 and diffAD==0:
        cat="B"
        if isoB>(isoA+isoD)/2:
            cat+="+"
        else:
            cat+="-"
    elif diffAB==0 and diffBD>0 and diffAD>0:
        cat="D"
        if isoD>(isoB+isoA)/2:
            cat+="+"
        else:
            cat+="-"
    elif diffAB>0 and diffBD>0 and diffAD>0:
        cat="F"
        if isoA>isoB and isoB>isoD:
            cat+="_ABD"
        elif isoD>isoA and isoA>isoB:
            cat+="_DAB"
        elif isoB>isoA and isoA>isoD:
            cat+="_BAD"
        elif isoA>isoD and isoD>isoB:
            cat+="_ADB"
        elif isoD>isoB and isoB>isoA:
            cat+="_DBA"
        elif isoB>isoD and isoD>isoA:
            cat+="_BDA"
    return cat

homeo_iso=pd.read_csv('Output/homoeolog_isoforms_detailed.csv', sep=",", header=0)
homeo_iso["category"]=homeo_iso.apply(lambda row : category(row["Isoform_number_A"],row["Isoform_number_B"],row["Isoform_number_D"]), axis = 1)
homeo_iso["AB_iso_diff"]=homeo_iso.apply(lambda row : abs(row["Isoform_number_A"]-row["Isoform_number_B"]), axis = 1)
homeo_iso["BD_iso_diff"]=homeo_iso.apply(lambda row : abs(row["Isoform_number_B"]-row["Isoform_number_D"]), axis = 1)
homeo_iso["AD_iso_diff"]=homeo_iso.apply(lambda row : abs(row["Isoform_number_A"]-row["Isoform_number_D"]), axis = 1)
print(homeo_iso)

homeo_iso_cat=pd.concat([homeo_iso.iloc[:,1:8],homeo_iso.iloc[:,13]],axis=1)
print(homeo_iso_cat)
homeo_iso_cat.to_csv('Output/homeo_iso_category.csv')

homeo_iso_cat_det=pd.concat([homeo_iso.iloc[:,1:11],homeo_iso.iloc[:,13]],axis=1)
print(homeo_iso_cat_det)
homeo_iso_cat_det.to_csv('Output/homeo_iso_category_detail.csv')

homeo_iso_diff=pd.concat([homeo_iso.iloc[:,1:8],homeo_iso.iloc[:,12:17]],axis=1)
print(homeo_iso_diff)
homeo_iso_diff.to_csv('Output/homeo_iso_diff.csv')

counts=homeo_iso.groupby(["category"]).count()
print(counts)
counts.to_csv('Output/isoform_category_counts.csv')