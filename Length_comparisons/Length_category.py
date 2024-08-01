import pandas as pd
import math

def category(lengthA,lengthB,lengthD, length_diff):
    cat="Unknown"
    if not math.isnan(length_diff):
        diffAB=abs(lengthA-lengthB)
        diffBD=abs(lengthB-lengthD)
        diffAD=abs(lengthA-lengthD)
        diffs=[diffAB,diffBD,diffAD]
        count=sum(1 for i in diffs if i>200)
        if count<2:
            cat="S"
        elif diffAB>200 and diffBD<200 and diffAD>200:
            cat="A"
            if lengthA>(lengthB+lengthD)/2:
                cat+="+"
            else:
                cat+="-"
        elif diffAB>200 and diffBD>200 and diffAD<200:
            cat="B"
            if lengthB>(lengthA+lengthD)/2:
                cat+="+"
            else:
                cat+="-"
        elif diffAB<200 and diffBD>200 and diffAD>200:
            cat="D"
            if lengthD>(lengthB+lengthA)/2:
                cat+="+"
            else:
                cat+="-"
        elif diffAB>200 and diffBD>200 and diffAD>200:
            cat="F"
            if lengthA>lengthB and lengthB>lengthD:
                cat+="_ABD"
            elif lengthD>lengthA and lengthA>lengthB:
                cat+="_DAB"
            elif lengthB>lengthA and lengthA>lengthD:
                cat+="_BAD"
            elif lengthA>lengthD and lengthD>lengthB:
                cat+="_ADB"
            elif lengthD>lengthB and lengthB>lengthA:
                cat+="_DBA"
            elif lengthB>lengthD and lengthD>lengthA:
                cat+="_BDA"
    return cat

def flag(High_std_A,High_std_B,High_std_D):
    return sum([High_std_A,High_std_B,High_std_D])>0


homeo_length=pd.read_csv('Output/homoeolog_length_detailed.csv', sep=",", header=0)
homeo_length["category"]=homeo_length.apply(lambda row : category(row["length_A"],row["length_B"],row["length_D"],row["length_diff"]), axis = 1)
homeo_length["High_std_flag"]=homeo_length.apply(lambda row : flag(row["High_std_A"],row["High_std_B"],row["High_std_D"]), axis = 1)
homeo_length["AB_length_diff"]=homeo_length.apply(lambda row : abs(row["length_A"]-row["length_B"]), axis = 1)
homeo_length["BD_length_diff"]=homeo_length.apply(lambda row : abs(row["length_B"]-row["length_D"]), axis = 1)
homeo_length["AD_length_diff"]=homeo_length.apply(lambda row : abs(row["length_A"]-row["length_D"]), axis = 1)
print(homeo_length)

homeo_length_cat=pd.concat([homeo_length.iloc[:,1:8],homeo_length.iloc[:,15:18]],axis=1)
homeo_length_cat.to_csv('Output/homeo_length_category.csv')

homeo_length_diff=pd.concat([homeo_length.iloc[:,1:8],homeo_length.iloc[:,15],homeo_length.iloc[:,17],homeo_length.iloc[:,18:21]],axis=1)
homeo_length_diff.to_csv('Output/homeo_length_diff.csv')


counts=homeo_length.groupby(["category"]).count()
print(counts)
counts.to_csv('Output/length_category_counts.csv')