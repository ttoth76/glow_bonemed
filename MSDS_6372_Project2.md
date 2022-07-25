MSDS_6372_Project2
================
Tamas Toth
2022-07-25

#### Loading the necessary R libraries for the analysis

``` r
# Load the necessary libraries
library(knitr)
library(rmarkdown)
# From class
library(GGally)
library(epitools)
library(MASS)
library(tidyverse)
library(car)
library(caret)
library(glmnet)
library(ROCR)
#Dasta set
library(aplore3)

# Other useful libraries
library(dplyr)
library(ggplot2)
library(gplots)

#library(ggpubr)
#library(tidyr)
#library(plyr)
#library(ggthemes)
#library(e1071)
#library(class)
#library(stringr)
#library(sjPlot)
#library(data.table)
#library(reshape2)
#library(corrplot)
#library(naivebayes)
#library(egg)
#library(rworldmap)
#library(Hmisc)
#library(DataExplorer)
#library(selectiveInference)
#library(dlookr)
```

``` r
# Turn off scientific notation
options(scipen = 100, digits = 4)
```

## The Global Longitudinal Study of Osteoporosis in Women

### Objective: Assessing risk factors and predicting if a woman with osteoperosis will have a bone fracture within the first year after joining the study.

#### The study has enrolled over 60,000 women aged 55 and older in ten countries. The major goals of the study are to use the data to provide insights into the management of fracture risk, patient experience with prevention and treatment of fractures and distribution of risk factors among older women on an international scale over the follow up period. The outcome variable is any fracture in the ﬁrst year of follow up. www.outcomes-umassmed.org/glow

## Objective 1 methodology:

1.  Understand the data
2.  EDA
3.  Feature Selection (Penalized Logistic Regression,
    Stepwise/forward/backward, Manual)
4.  Split the data to Training and Test set
5.  Fit Logistic Regression model
6.  Interpret the model, including hypothesis testing and confidence
    intervals
7.  Conclusion

#### Read the data

``` r
#Read the data
bonemed_df = glow_bonemed
attach(bonemed_df)
bonemed_df_sample = sample_n(bonemed_df, 5)
knitr::kable(bonemed_df_sample, "html")
```

<table>
<thead>
<tr>
<th style="text-align:right;">
sub_id
</th>
<th style="text-align:right;">
site_id
</th>
<th style="text-align:right;">
phy_id
</th>
<th style="text-align:left;">
priorfrac
</th>
<th style="text-align:right;">
age
</th>
<th style="text-align:right;">
weight
</th>
<th style="text-align:right;">
height
</th>
<th style="text-align:right;">
bmi
</th>
<th style="text-align:left;">
premeno
</th>
<th style="text-align:left;">
momfrac
</th>
<th style="text-align:left;">
armassist
</th>
<th style="text-align:left;">
smoke
</th>
<th style="text-align:left;">
raterisk
</th>
<th style="text-align:right;">
fracscore
</th>
<th style="text-align:left;">
fracture
</th>
<th style="text-align:left;">
bonemed
</th>
<th style="text-align:left;">
bonemed_fu
</th>
<th style="text-align:left;">
bonetreat
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:right;">
198
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
21
</td>
<td style="text-align:left;">
No
</td>
<td style="text-align:right;">
60
</td>
<td style="text-align:right;">
79.4
</td>
<td style="text-align:right;">
168
</td>
<td style="text-align:right;">
28.13
</td>
<td style="text-align:left;">
Yes
</td>
<td style="text-align:left;">
Yes
</td>
<td style="text-align:left;">
No
</td>
<td style="text-align:left;">
No
</td>
<td style="text-align:left;">
Less
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:left;">
No
</td>
<td style="text-align:left;">
No
</td>
<td style="text-align:left;">
No
</td>
<td style="text-align:left;">
No
</td>
</tr>
<tr>
<td style="text-align:right;">
74
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
63
</td>
<td style="text-align:left;">
No
</td>
<td style="text-align:right;">
63
</td>
<td style="text-align:right;">
63.5
</td>
<td style="text-align:right;">
165
</td>
<td style="text-align:right;">
23.32
</td>
<td style="text-align:left;">
Yes
</td>
<td style="text-align:left;">
No
</td>
<td style="text-align:left;">
No
</td>
<td style="text-align:left;">
No
</td>
<td style="text-align:left;">
Less
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:left;">
No
</td>
<td style="text-align:left;">
No
</td>
<td style="text-align:left;">
No
</td>
<td style="text-align:left;">
No
</td>
</tr>
<tr>
<td style="text-align:right;">
62
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
287
</td>
<td style="text-align:left;">
No
</td>
<td style="text-align:right;">
60
</td>
<td style="text-align:right;">
80.7
</td>
<td style="text-align:right;">
160
</td>
<td style="text-align:right;">
31.52
</td>
<td style="text-align:left;">
No
</td>
<td style="text-align:left;">
No
</td>
<td style="text-align:left;">
No
</td>
<td style="text-align:left;">
No
</td>
<td style="text-align:left;">
Less
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:left;">
No
</td>
<td style="text-align:left;">
No
</td>
<td style="text-align:left;">
No
</td>
<td style="text-align:left;">
No
</td>
</tr>
<tr>
<td style="text-align:right;">
282
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
137
</td>
<td style="text-align:left;">
No
</td>
<td style="text-align:right;">
58
</td>
<td style="text-align:right;">
61.2
</td>
<td style="text-align:right;">
168
</td>
<td style="text-align:right;">
21.68
</td>
<td style="text-align:left;">
No
</td>
<td style="text-align:left;">
No
</td>
<td style="text-align:left;">
No
</td>
<td style="text-align:left;">
No
</td>
<td style="text-align:left;">
Greater
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:left;">
No
</td>
<td style="text-align:left;">
No
</td>
<td style="text-align:left;">
No
</td>
<td style="text-align:left;">
No
</td>
</tr>
<tr>
<td style="text-align:right;">
152
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
292
</td>
<td style="text-align:left;">
Yes
</td>
<td style="text-align:right;">
61
</td>
<td style="text-align:right;">
87.5
</td>
<td style="text-align:right;">
167
</td>
<td style="text-align:right;">
31.37
</td>
<td style="text-align:left;">
Yes
</td>
<td style="text-align:left;">
No
</td>
<td style="text-align:left;">
Yes
</td>
<td style="text-align:left;">
No
</td>
<td style="text-align:left;">
Greater
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:left;">
No
</td>
<td style="text-align:left;">
No
</td>
<td style="text-align:left;">
No
</td>
<td style="text-align:left;">
No
</td>
</tr>
</tbody>
</table>

``` r
dim(bonemed_df)
```

    ## [1] 500  18

### Data Description

-   **sub_id**: Identification Code (1 - n)

-   **site_id**: Study Site (1 - 6)

-   **phy_id**: Physician ID code (128 unique codes)

-   **priorfrac**: History of Prior Fracture (1: No, 2: Yes)

-   **age**: Age at Enrollment (Years)

-   **weight**: Weight at enrollment (Kilograms)

-   **height**: Height at enrollment (Centimeters)

-   **bmi**: Body Mass Index (Kg/m^2)

-   **premeno**: Menopause before age 45 (1: No, 2: Yes)

-   **momfrac**: Mother had hip fracture (1: No, 2: Yes)

-   **armassist**: Arms are needed to stand from a chair (1: No, 2: Yes)

-   **smoke**: Former or current smoker (1: No, 2: Yes)

-   **raterisk**: Self-reported risk of fracture (1: Less than others of
    the same age, 2: Same as others of the same age, 3: Greater than
    others of the same age)

-   **fracscore**: Fracture Risk Score (Composite Risk Score)

-   **fracture**: Any fracture in first year (1: No, 2: Yes)

-   **bonemed**: Bone medications at enrollment (1: No, 2: Yes)

-   **bonemed_fu**: Bone medications at follow-up (1: No, 2: Yes)

-   **bonetreat**: Bone medications both at enrollment and follow-up (1:
    No, 2: Yes)

``` r
#set random seed
set.seed(329)
```

#### Address the missing values in each column (NA as well as empty strings).

``` r
# Address the missing values in each column (NA as well as empty strings).
missing_df = as.data.frame(sapply(bonemed_df, function(x) sum(is.na(x))))
colnames(missing_df) = c("missing values")
knitr::kable(missing_df, "html")
```

<table>
<thead>
<tr>
<th style="text-align:left;">
</th>
<th style="text-align:right;">
missing values
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
sub_id
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
site_id
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
phy_id
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
priorfrac
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
age
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
weight
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
height
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
bmi
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
premeno
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
momfrac
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
armassist
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
smoke
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
raterisk
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
fracscore
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
fracture
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
bonemed
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
bonemed_fu
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
bonetreat
</td>
<td style="text-align:right;">
0
</td>
</tr>
</tbody>
</table>

``` r
empty_string_df = as.data.frame(sapply(bonemed_df, function(x) sum(x == "")))
colnames(empty_string_df) = c("empty string")
knitr::kable(empty_string_df, "html")
```

<table>
<thead>
<tr>
<th style="text-align:left;">
</th>
<th style="text-align:right;">
empty string
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
sub_id
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
site_id
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
phy_id
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
priorfrac
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
age
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
weight
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
height
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
bmi
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
premeno
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
momfrac
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
armassist
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
smoke
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
raterisk
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
fracscore
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
fracture
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
bonemed
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
bonemed_fu
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
bonetreat
</td>
<td style="text-align:right;">
0
</td>
</tr>
</tbody>
</table>

``` r
# Function to Identify different characteristics of the data frame 
# Getting a concise summary of the dataframe: str()
# Listing the column labels of the dataframe: colnames()
# Size of the dataset: dim()
# # Verify if there is any negative values in the dataset
dfinfo = function(df_name)
  {
  df_structure = str(df_name)
  df_colnames = colnames(df_name)
  df_dimensions = dim(df_name)
  
  num_cols = bonemed_df %>% dplyr::select(where(is.numeric)) %>% colnames()
  df_neg = print(paste("Negative values in the variable:",  
                       sapply(bonemed_df[,num_cols], function(x) sum(x < 0))))
  
  outparam = list(df_structure, df_colnames, df_dimensions, df_neg)
  return (outparam)
}
```

``` r
dfinfo(bonemed_df)
```

    ## 'data.frame':    500 obs. of  18 variables:
    ##  $ sub_id    : int  1 2 3 4 5 6 7 8 9 10 ...
    ##  $ site_id   : int  1 4 6 6 1 5 5 1 1 4 ...
    ##  $ phy_id    : int  14 284 305 309 37 299 302 36 8 282 ...
    ##  $ priorfrac : Factor w/ 2 levels "No","Yes": 1 1 2 1 1 2 1 2 2 1 ...
    ##  $ age       : int  62 65 88 82 61 67 84 82 86 58 ...
    ##  $ weight    : num  70.3 87.1 50.8 62.1 68 68 50.8 40.8 62.6 63.5 ...
    ##  $ height    : int  158 160 157 160 152 161 150 153 156 166 ...
    ##  $ bmi       : num  28.2 34 20.6 24.3 29.4 ...
    ##  $ premeno   : Factor w/ 2 levels "No","Yes": 1 1 1 1 1 1 1 1 1 1 ...
    ##  $ momfrac   : Factor w/ 2 levels "No","Yes": 1 1 2 1 1 1 1 1 1 1 ...
    ##  $ armassist : Factor w/ 2 levels "No","Yes": 1 1 2 1 1 1 1 1 1 1 ...
    ##  $ smoke     : Factor w/ 2 levels "No","Yes": 1 1 1 1 1 2 1 1 1 1 ...
    ##  $ raterisk  : Factor w/ 3 levels "Less","Same",..: 2 2 1 1 2 2 1 2 2 1 ...
    ##  $ fracscore : int  1 2 11 5 1 4 6 7 7 0 ...
    ##  $ fracture  : Factor w/ 2 levels "No","Yes": 1 1 1 1 1 1 1 1 1 1 ...
    ##  $ bonemed   : Factor w/ 2 levels "No","Yes": 1 1 1 1 1 1 1 2 1 1 ...
    ##  $ bonemed_fu: Factor w/ 2 levels "No","Yes": 1 1 1 1 1 1 1 2 1 1 ...
    ##  $ bonetreat : Factor w/ 2 levels "No","Yes": 1 1 1 1 1 1 1 2 1 1 ...
    ## [1] "Negative values in the variable: 0" "Negative values in the variable: 0"
    ## [3] "Negative values in the variable: 0" "Negative values in the variable: 0"
    ## [5] "Negative values in the variable: 0" "Negative values in the variable: 0"
    ## [7] "Negative values in the variable: 0" "Negative values in the variable: 0"

    ## [[1]]
    ## NULL
    ## 
    ## [[2]]
    ##  [1] "sub_id"     "site_id"    "phy_id"     "priorfrac"  "age"       
    ##  [6] "weight"     "height"     "bmi"        "premeno"    "momfrac"   
    ## [11] "armassist"  "smoke"      "raterisk"   "fracscore"  "fracture"  
    ## [16] "bonemed"    "bonemed_fu" "bonetreat" 
    ## 
    ## [[3]]
    ## [1] 500  18
    ## 
    ## [[4]]
    ## [1] "Negative values in the variable: 0" "Negative values in the variable: 0"
    ## [3] "Negative values in the variable: 0" "Negative values in the variable: 0"
    ## [5] "Negative values in the variable: 0" "Negative values in the variable: 0"
    ## [7] "Negative values in the variable: 0" "Negative values in the variable: 0"

### Observations:

-   The data set is comprised of 500 observations and 18 variables
-   There are numerical and categorical variables in the data set
-   There are no missing values or empty strings in the data set
-   No negative values in the data set
-   ??????No duplicated records
-   ‘fracture’ is the dependent variable

##################################################################################### 

# Exploratory Data Analysis

##################################################################################### 

``` r
num_cols = bonemed_df %>% dplyr::select(where(is.numeric)) %>% colnames()
pair_plot = c(num_cols, 'fracture')
pair_plot = pair_plot[-1]
ggpairs(bonemed_df[,pair_plot],aes(color=fracture, alpha = 0.5))
```

<img src="MSDS_6372_Project2_files/figure-gfm/pair plots-1.png" angle=90 style="display: block; margin: auto;" />

### Observations:

-   Height for fracture ‘No’ and ‘Yes’ levels seem to be normally
    distributed.
-   Weight and bmi are strongly positively linearly correlated for both
    factor levels.
-   Age and fracscore are strongly positively linearly correlated for
    both factor levels.
-   There are more “No” fracture in the first year observations than
    “Yes”

``` r
# Summary statistics
t(aggregate(.~ fracture,data=bonemed_df,summary))
```

    ##                    [,1]     [,2]    
    ## fracture           "No"     "Yes"   
    ## sub_id.Min.        "  1.0"  "376.0" 
    ## sub_id.1st Qu.     " 94.5"  "407.0" 
    ## sub_id.Median      "188.0"  "438.0" 
    ## sub_id.Mean        "188.0"  "438.0" 
    ## sub_id.3rd Qu.     "281.5"  "469.0" 
    ## sub_id.Max.        "375.0"  "500.0" 
    ## site_id.Min.       "1.000"  "1.000" 
    ## site_id.1st Qu.    "2.000"  "2.000" 
    ## site_id.Median     "3.000"  "4.000" 
    ## site_id.Mean       "3.363"  "3.656" 
    ## site_id.3rd Qu.    "5.000"  "5.000" 
    ## site_id.Max.       "6.000"  "6.000" 
    ## phy_id.Min.        "  1.0"  "  7.0" 
    ## phy_id.1st Qu.     " 50.0"  " 68.0" 
    ## phy_id.Median      "172.0"  "197.0" 
    ## phy_id.Mean        "173.9"  "192.5" 
    ## phy_id.3rd Qu.     "296.0"  "300.0" 
    ## phy_id.Max.        "325.0"  "325.0" 
    ## priorfrac.Min.     "1.000"  "1.000" 
    ## priorfrac.1st Qu.  "1.000"  "1.000" 
    ## priorfrac.Median   "1.000"  "1.000" 
    ## priorfrac.Mean     "1.197"  "1.416" 
    ## priorfrac.3rd Qu.  "1.000"  "2.000" 
    ## priorfrac.Max.     "2.000"  "2.000" 
    ## age.Min.           "55.00"  "56.00" 
    ## age.1st Qu.        "60.00"  "65.00" 
    ## age.Median         "66.00"  "72.00" 
    ## age.Mean           "67.49"  "71.79" 
    ## age.3rd Qu.        "74.00"  "79.00" 
    ## age.Max.           "90.00"  "89.00" 
    ## weight.Min.        " 39.90" " 45.80"
    ## weight.1st Qu.     " 60.30" " 59.90"
    ## weight.Median      " 68.00" " 68.00"
    ## weight.Mean        " 72.17" " 70.79"
    ## weight.3rd Qu.     " 81.60" " 79.40"
    ## weight.Max.        "127.00" "124.70"
    ## height.Min.        "142.0"  "134.0" 
    ## height.1st Qu.     "158.0"  "155.0" 
    ## height.Median      "162.0"  "160.0" 
    ## height.Mean        "161.9"  "159.9" 
    ## height.3rd Qu.     "166.0"  "164.0" 
    ## height.Max.        "199.0"  "178.0" 
    ## bmi.Min.           "14.88"  "17.04" 
    ## bmi.1st Qu.        "23.32"  "23.05" 
    ## bmi.Median         "26.37"  "26.43" 
    ## bmi.Mean           "27.50"  "27.71" 
    ## bmi.3rd Qu.        "30.62"  "31.09" 
    ## bmi.Max.           "49.08"  "44.04" 
    ## premeno.Min.       "1.000"  "1.000" 
    ## premeno.1st Qu.    "1.000"  "1.000" 
    ## premeno.Median     "1.000"  "1.000" 
    ## premeno.Mean       "1.192"  "1.200" 
    ## premeno.3rd Qu.    "1.000"  "1.000" 
    ## premeno.Max.       "2.000"  "2.000" 
    ## momfrac.Min.       "1.000"  "1.000" 
    ## momfrac.1st Qu.    "1.000"  "1.000" 
    ## momfrac.Median     "1.000"  "1.000" 
    ## momfrac.Mean       "1.109"  "1.192" 
    ## momfrac.3rd Qu.    "1.000"  "1.000" 
    ## momfrac.Max.       "2.000"  "2.000" 
    ## armassist.Min.     "1.000"  "1.000" 
    ## armassist.1st Qu.  "1.000"  "1.000" 
    ## armassist.Median   "1.000"  "2.000" 
    ## armassist.Mean     "1.333"  "1.504" 
    ## armassist.3rd Qu.  "2.000"  "2.000" 
    ## armassist.Max.     "2.000"  "2.000" 
    ## smoke.Min.         "1.000"  "1.000" 
    ## smoke.1st Qu.      "1.000"  "1.000" 
    ## smoke.Median       "1.000"  "1.000" 
    ## smoke.Mean         "1.075"  "1.056" 
    ## smoke.3rd Qu.      "1.000"  "1.000" 
    ## smoke.Max.         "2.000"  "2.000" 
    ## raterisk.Min.      "1.000"  "1.000" 
    ## raterisk.1st Qu.   "1.000"  "2.000" 
    ## raterisk.Median    "2.000"  "2.000" 
    ## raterisk.Mean      "1.891"  "2.168" 
    ## raterisk.3rd Qu.   "3.000"  "3.000" 
    ## raterisk.Max.      "3.000"  "3.000" 
    ## fracscore.Min.     " 0.000" " 0.000"
    ## fracscore.1st Qu.  " 1.500" " 3.000"
    ## fracscore.Median   " 3.000" " 5.000"
    ## fracscore.Mean     " 3.317" " 4.840"
    ## fracscore.3rd Qu.  " 5.000" " 7.000"
    ## fracscore.Max.     "11.000" " 9.000"
    ## bonemed.Min.       "1.000"  "1.000" 
    ## bonemed.1st Qu.    "1.000"  "1.000" 
    ## bonemed.Median     "1.000"  "1.000" 
    ## bonemed.Mean       "1.221"  "1.368" 
    ## bonemed.3rd Qu.    "1.000"  "2.000" 
    ## bonemed.Max.       "2.000"  "2.000" 
    ## bonemed_fu.Min.    "1.000"  "1.000" 
    ## bonemed_fu.1st Qu. "1.000"  "1.000" 
    ## bonemed_fu.Median  "1.000"  "1.000" 
    ## bonemed_fu.Mean    "1.229"  "1.424" 
    ## bonemed_fu.3rd Qu. "1.000"  "2.000" 
    ## bonemed_fu.Max.    "2.000"  "2.000" 
    ## bonetreat.Min.     "1.000"  "1.000" 
    ## bonetreat.1st Qu.  "1.000"  "1.000" 
    ## bonetreat.Median   "1.000"  "1.000" 
    ## bonetreat.Mean     "1.208"  "1.320" 
    ## bonetreat.3rd Qu.  "1.000"  "2.000" 
    ## bonetreat.Max.     "2.000"  "2.000"

### Observations:

-   The minimum age for “No” fracture in the first year is 55. The
    minimum age for “Yes” fracture in the first year is 56.

-   The maximum age for “No” fracture in the first year is 99. The
    maximum age for “Yes” fracture in the first year is 89.

-   The median age for “No” fracture in the first year is 66. The median
    age for “Yes” fracture in the first year is 72.

-   The mean age for “No” fracture in the first year is 67.49. The mean
    age for “Yes” fracture in the first year is 71.79.

-   Age is a likely is good predictor of fracture in the first year
    since there is a observable difference in median and mean age of
    having a fracture.

-   The minimum weight for “No” fracture in the first year is 39.9. The
    minimum weight for “Yes” fracture in the first year is 45.8.

-   The maximum weight for “No” fracture in the first year is 127. The
    maximum weight for “Yes” fracture in the first year is 124.7.

-   The median weight for “No” fracture in the first year is 68. The
    median weight for “Yes” fracture in the first year is 68.

-   The mean weight for “No” fracture in the first year is 72.17. The
    mean weight for “Yes” fracture in the first year is 70.79.

-   Weight does not seem to be a significant variable for predicting
    fracture in the first year as the median and mean weight values are
    very similar to the factor levels.

-   The minimum height for “No” fracture in the first year is 142. The
    minimum height for “Yes” fracture in the first year is 134.

-   The maximum height for “No” fracture in the first year is 199. The
    maximum height for “Yes” fracture in the first year is 178.

-   The median height for “No” fracture in the first year is 162. The
    median height for “Yes” fracture in the first year is 160.

-   The mean height for “No” fracture in the first year is 161.9. The
    mean height for “Yes” fracture in the first year is 159.9.

-   Height seems to be a good predictor for fracture in the first year.

-   The minimum bmi for “No” fracture in the first year is 14.88. The
    minimum bmi for “Yes” fracture in the first year is 17.04.

-   The maximum bmi for “No” fracture in the first year is 49.08 The
    maximum bmi for “Yes” fracture in the first year is 44.04.

-   The median bmi for “No” fracture in the first year is 26.37. The
    median bmi for “Yes” fracture in the first year is 26.43.

-   The mean bmi for “No” fracture in the first year is 27.50. The mean
    bmi for “Yes” fracture in the first year is 27.71.

-   BMI does not seem to be a significant predictor for fracture in the
    first year.

-   The minimum Fracture Risk Score for “No” fracture in the first year
    is 0. The minimum Fracture Risk Score for “Yes” fracture in the
    first year is 0.

-   The maximum Fracture Risk Score for “No” fracture in the first year
    is 11. The maximum Fracture Risk Score for “Yes” fracture in the
    first year is 9.

-   The median Fracture Risk Score for “No” fracture in the first year
    is 3. The median Fracture Risk Score for “Yes” fracture in the first
    year is 5.

-   The mean Fracture Risk Score for “No” fracture in the first year is
    3.317. The mean Fracture Risk Score for “Yes” fracture in the first
    year is 4.840.

-   Fracture Risk Score seems to be a significant predictor for fracture
    in the first year.

##################################################################################### 

# Categorical data plots

##################################################################################### 

``` r
cat_cols = bonemed_df %>% dplyr::select(where(is.factor)) %>% colnames()

# Plot all categorical variables
for (c in cat_cols)
{
  cat_plot = bonemed_df %>% ggplot(aes(x= .data[[c]], group = 1)) + 
    geom_bar(aes(y = ..prop.., fill = factor(..x..)), stat="count") +
    geom_text(aes( label = scales::percent(..prop..),
                   y= ..prop.. ), stat= "count", vjust = -.5) +
    labs(y = "Percent") +
    scale_y_continuous(labels = scales::percent) + theme(legend.position = "none") +
    ggtitle(paste(c, "Categorical Analysis")) + 
    theme(plot.title = element_text(hjust = 0.5)) + 
    theme(axis.text.x = element_text(vjust = 0.5, hjust=1)) # +
    egg::ggarrange(cat_plot, ncol=2) 
}
```

<img src="MSDS_6372_Project2_files/figure-gfm/unnamed-chunk-9-1.png" angle=90 style="display: block; margin: auto;" /><img src="MSDS_6372_Project2_files/figure-gfm/unnamed-chunk-9-2.png" angle=90 style="display: block; margin: auto;" /><img src="MSDS_6372_Project2_files/figure-gfm/unnamed-chunk-9-3.png" angle=90 style="display: block; margin: auto;" /><img src="MSDS_6372_Project2_files/figure-gfm/unnamed-chunk-9-4.png" angle=90 style="display: block; margin: auto;" /><img src="MSDS_6372_Project2_files/figure-gfm/unnamed-chunk-9-5.png" angle=90 style="display: block; margin: auto;" /><img src="MSDS_6372_Project2_files/figure-gfm/unnamed-chunk-9-6.png" angle=90 style="display: block; margin: auto;" /><img src="MSDS_6372_Project2_files/figure-gfm/unnamed-chunk-9-7.png" angle=90 style="display: block; margin: auto;" /><img src="MSDS_6372_Project2_files/figure-gfm/unnamed-chunk-9-8.png" angle=90 style="display: block; margin: auto;" /><img src="MSDS_6372_Project2_files/figure-gfm/unnamed-chunk-9-9.png" angle=90 style="display: block; margin: auto;" /><img src="MSDS_6372_Project2_files/figure-gfm/unnamed-chunk-9-10.png" angle=90 style="display: block; margin: auto;" />

### Observations:

-   75% of the observations have no prior history of fracture. 25% of
    the observations have prior history of fracture.
-   81% of the observations have no menopause before age 45. 19% of the
    observations have menopause before age 45.
-   87% of the subjects’ mother had no hip fracture. 13% of the
    subjects’ mother had hip fracture.
-   62% of the subjects do not need arms to stand from a chair. 38% of
    the subjects need arms to stand from a chair.
-   93% of the subjects don’t smoke. 7% of the subjects smoke are former
    or current smoker.
-   33.4% of the subjects reported less risk of fractures than others of
    the same age. 37.2% of the subjects reported the same risk of
    fractures than others of the same age. 29.4% of the subjects
    reported greater risk of fractures than others of the same age.
-   75% of the patients do not have any fractures in first year. 25% of
    the patients have any fractures in first year.
-   74% of the subjects have enrolled to Bone medications. 26% of the
    subjects have not enrolled to Bone medications.
-   72% of the subjects did not require bone medications at follow-up.
    28% of the subjects required bone medications at follow-up.
-   76% of the patients did not receive bone medications at either
    enrollment or follow-up. 24% of the patients received bone
    medications both at enrollment and follow-up.

##################################################################################### 

# Bi-variate analysis with Fracture variable

##################################################################################### 

``` r
num_cols = bonemed_df %>% dplyr::select(where(is.numeric)) %>% colnames()
bivar_plot = num_cols[c(-1, -2, -3)]


for (i in bivar_plot)
{
multibox = bonemed_df %>%
  ggplot(aes(x=fracture, y = .data[[i]])) +
  geom_boxplot(fill = "sandybrown", color = "black") + 
  xlab("Fracture") +
  ylab(i) + stat_summary(fun=mean, geom="point", shape=20, size=7, color="red", fill="red") +
  ggtitle(paste(i, "vs Fracture bi-variate analysis")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette = "Oranges")  
egg::ggarrange(multibox, ncol=2)
}
```

<img src="MSDS_6372_Project2_files/figure-gfm/unnamed-chunk-10-1.png" angle=90 style="display: block; margin: auto;" /><img src="MSDS_6372_Project2_files/figure-gfm/unnamed-chunk-10-2.png" angle=90 style="display: block; margin: auto;" /><img src="MSDS_6372_Project2_files/figure-gfm/unnamed-chunk-10-3.png" angle=90 style="display: block; margin: auto;" /><img src="MSDS_6372_Project2_files/figure-gfm/unnamed-chunk-10-4.png" angle=90 style="display: block; margin: auto;" /><img src="MSDS_6372_Project2_files/figure-gfm/unnamed-chunk-10-5.png" angle=90 style="display: block; margin: auto;" />

### Observations:

-   We can see an increasing in median and mean age for those who had
    any fractures in the first year.
-   Median and mean weight for both fracture levels are very similar.
-   The median and mean height are lower for those who had fracture in
    the first year. This is likely because the bone density is reduced
    and the subjects shrunk.
-   Median and mean BMI is very similar for both factor levels.
-   The mean and median fracture risk score is higher for those who had
    any fracture in the first year.

##################################################################################### 

# Correlation plot

##################################################################################### 

``` r
corr_df = bonemed_df[,c('age', 'weight', 'height', 'bmi', 'fracscore')]
cont_var.cor = cor(corr_df)
heatmap.2(cont_var.cor,col=redgreen(75), 
          density.info="none", trace="none", dendrogram=c("row"), 
          symm=F,symkey=T,symbreaks=T, scale="none")
```

![](MSDS_6372_Project2_files/figure-gfm/correlation%20plot-1.png)<!-- -->

### Observations:

-   weight and bmi are more similar to each other therefore they would
    form a cluster
-   height, weight and bmi are also similar but more distant from each
    other. Still can form a cluster.
-   age and fracscore are also similar to each other and can form a
    cluster

``` r
# Label encoding Yes=1; No=0
bonemed_df$fracture.num<-ifelse(bonemed_df$fracture=="Yes",1,0)
```

##################################################################################### 

# Loess plots

##################################################################################### 

``` r
num_cols = bonemed_df %>% dplyr::select(where(is.numeric)) %>% colnames()
loess_plot = num_cols[c(-1, -2, -3, -9)]

for (i in loess_plot)
{
loess = bonemed_df %>% 
ggplot(aes(x=.data[[i]],y=fracture.num))+
geom_point()+
geom_smooth(formula = y ~ x, method="loess")+
theme(plot.title = element_text(hjust = 0.5)) +
ggtitle(paste(i, "vs fracture loess smoothing"))
egg::ggarrange(loess, ncol=2)
}
```

![](MSDS_6372_Project2_files/figure-gfm/Loess%20plot-1.png)<!-- -->![](MSDS_6372_Project2_files/figure-gfm/Loess%20plot-2.png)<!-- -->![](MSDS_6372_Project2_files/figure-gfm/Loess%20plot-3.png)<!-- -->![](MSDS_6372_Project2_files/figure-gfm/Loess%20plot-4.png)<!-- -->![](MSDS_6372_Project2_files/figure-gfm/Loess%20plot-5.png)<!-- -->

### Observations:

-   The Loess plots show the ‘height’ has an S curve like logistic
    model. The curve is trending down.
-   The other continuous predictors don’t seem to be important
    predicting fracture in the first year.

##################################################################################### 

# Loess plots to investigate interactions

##################################################################################### 

``` r
num_cols = bonemed_df %>% dplyr::select(where(is.numeric)) %>% colnames()
loess_plot = num_cols[c(-1, -2, -3, -9)]

for (j in cat_cols)
{
  for (i in loess_plot)
  {
  plot1 = ggplot(bonemed_df,aes(x=.data[[i]],y=fracture.num,colour=.data[[j]]))+geom_point()+
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle(paste(i, "&", j, " interaction with fracture"))+
  geom_smooth(formula = y ~ x, method="loess",size=1,span=1.5)+facet_wrap(~.data[[j]])
  ylim(-.2,1.2)
  show(plot1)
  }
}
```

![](MSDS_6372_Project2_files/figure-gfm/Loess%20plot%20interaction-1.png)<!-- -->![](MSDS_6372_Project2_files/figure-gfm/Loess%20plot%20interaction-2.png)<!-- -->![](MSDS_6372_Project2_files/figure-gfm/Loess%20plot%20interaction-3.png)<!-- -->![](MSDS_6372_Project2_files/figure-gfm/Loess%20plot%20interaction-4.png)<!-- -->![](MSDS_6372_Project2_files/figure-gfm/Loess%20plot%20interaction-5.png)<!-- -->![](MSDS_6372_Project2_files/figure-gfm/Loess%20plot%20interaction-6.png)<!-- -->![](MSDS_6372_Project2_files/figure-gfm/Loess%20plot%20interaction-7.png)<!-- -->![](MSDS_6372_Project2_files/figure-gfm/Loess%20plot%20interaction-8.png)<!-- -->![](MSDS_6372_Project2_files/figure-gfm/Loess%20plot%20interaction-9.png)<!-- -->![](MSDS_6372_Project2_files/figure-gfm/Loess%20plot%20interaction-10.png)<!-- -->![](MSDS_6372_Project2_files/figure-gfm/Loess%20plot%20interaction-11.png)<!-- -->![](MSDS_6372_Project2_files/figure-gfm/Loess%20plot%20interaction-12.png)<!-- -->![](MSDS_6372_Project2_files/figure-gfm/Loess%20plot%20interaction-13.png)<!-- -->![](MSDS_6372_Project2_files/figure-gfm/Loess%20plot%20interaction-14.png)<!-- -->![](MSDS_6372_Project2_files/figure-gfm/Loess%20plot%20interaction-15.png)<!-- -->![](MSDS_6372_Project2_files/figure-gfm/Loess%20plot%20interaction-16.png)<!-- -->![](MSDS_6372_Project2_files/figure-gfm/Loess%20plot%20interaction-17.png)<!-- -->![](MSDS_6372_Project2_files/figure-gfm/Loess%20plot%20interaction-18.png)<!-- -->![](MSDS_6372_Project2_files/figure-gfm/Loess%20plot%20interaction-19.png)<!-- -->![](MSDS_6372_Project2_files/figure-gfm/Loess%20plot%20interaction-20.png)<!-- -->![](MSDS_6372_Project2_files/figure-gfm/Loess%20plot%20interaction-21.png)<!-- -->![](MSDS_6372_Project2_files/figure-gfm/Loess%20plot%20interaction-22.png)<!-- -->![](MSDS_6372_Project2_files/figure-gfm/Loess%20plot%20interaction-23.png)<!-- -->![](MSDS_6372_Project2_files/figure-gfm/Loess%20plot%20interaction-24.png)<!-- -->![](MSDS_6372_Project2_files/figure-gfm/Loess%20plot%20interaction-25.png)<!-- -->![](MSDS_6372_Project2_files/figure-gfm/Loess%20plot%20interaction-26.png)<!-- -->![](MSDS_6372_Project2_files/figure-gfm/Loess%20plot%20interaction-27.png)<!-- -->![](MSDS_6372_Project2_files/figure-gfm/Loess%20plot%20interaction-28.png)<!-- -->![](MSDS_6372_Project2_files/figure-gfm/Loess%20plot%20interaction-29.png)<!-- -->![](MSDS_6372_Project2_files/figure-gfm/Loess%20plot%20interaction-30.png)<!-- -->![](MSDS_6372_Project2_files/figure-gfm/Loess%20plot%20interaction-31.png)<!-- -->![](MSDS_6372_Project2_files/figure-gfm/Loess%20plot%20interaction-32.png)<!-- -->![](MSDS_6372_Project2_files/figure-gfm/Loess%20plot%20interaction-33.png)<!-- -->![](MSDS_6372_Project2_files/figure-gfm/Loess%20plot%20interaction-34.png)<!-- -->![](MSDS_6372_Project2_files/figure-gfm/Loess%20plot%20interaction-35.png)<!-- -->![](MSDS_6372_Project2_files/figure-gfm/Loess%20plot%20interaction-36.png)<!-- -->![](MSDS_6372_Project2_files/figure-gfm/Loess%20plot%20interaction-37.png)<!-- -->![](MSDS_6372_Project2_files/figure-gfm/Loess%20plot%20interaction-38.png)<!-- -->![](MSDS_6372_Project2_files/figure-gfm/Loess%20plot%20interaction-39.png)<!-- -->![](MSDS_6372_Project2_files/figure-gfm/Loess%20plot%20interaction-40.png)<!-- -->![](MSDS_6372_Project2_files/figure-gfm/Loess%20plot%20interaction-41.png)<!-- -->![](MSDS_6372_Project2_files/figure-gfm/Loess%20plot%20interaction-42.png)<!-- -->![](MSDS_6372_Project2_files/figure-gfm/Loess%20plot%20interaction-43.png)<!-- -->![](MSDS_6372_Project2_files/figure-gfm/Loess%20plot%20interaction-44.png)<!-- -->![](MSDS_6372_Project2_files/figure-gfm/Loess%20plot%20interaction-45.png)<!-- -->![](MSDS_6372_Project2_files/figure-gfm/Loess%20plot%20interaction-46.png)<!-- -->![](MSDS_6372_Project2_files/figure-gfm/Loess%20plot%20interaction-47.png)<!-- -->![](MSDS_6372_Project2_files/figure-gfm/Loess%20plot%20interaction-48.png)<!-- -->![](MSDS_6372_Project2_files/figure-gfm/Loess%20plot%20interaction-49.png)<!-- -->![](MSDS_6372_Project2_files/figure-gfm/Loess%20plot%20interaction-50.png)<!-- -->

## Observations:

-   There is a interaction with weight and priorfrac
-   There is a interaction with bmi and priorfrac
-   There is a interaction with fracscore and priorfrac
-   There is a interaction with weight and raterisk
-   There is a interaction with bmi and raterisk
-   There is a interaction with weight and bonemed
-   There is a interaction with bmi and bonemed
-   There is a interaction with weight and bonemed_fu
-   There is a interaction with bmi and bonemed_fu
-   There is a interaction with weight and bonetreat
-   There is a interaction with bmi and bonetreat

##################################################################################### 

# Split the Data to Train and Test sets (85%-15%)

##################################################################################### 

``` r
bonmed_df_split = bonemed_df[,-c(1, 2, 3)]
index<-sample(1:dim(bonmed_df_split)[1],round(dim(bonmed_df_split)[1]*0.85),replace=F)
train = bonmed_df_split[index,]
test = bonmed_df_split[-index,]
```

##################################################################################### 

# Feature Selection

##################################################################################### 

## Manual / Intuition selection

``` r
# Fit each variable separately and check the p-value for significance

# priorfrac
simple.log.priorfrac<-glm(fracture~priorfrac,family=binomial(link='logit'),data=train)
simple.log.priorfrac.sum = summary(simple.log.priorfrac)

# age
simple.log.age<-glm(fracture~age,family=binomial(link='logit'),data=train)
simple.log.age.sum = summary(simple.log.age)

# weight
simple.log.weight<-glm(fracture~weight,family=binomial(link='logit'),data=train)
simple.log.weight.sum = summary(simple.log.weight)

# height
simple.log.height<-glm(fracture~height,family=binomial(link='logit'),data=train)
simple.log.height.sum = summary(simple.log.height)

# bmi
simple.log.bmi<-glm(fracture~bmi,family=binomial(link='logit'),data=train)
simple.log.bmi.sum = summary(simple.log.bmi)

# premeno
simple.log.premeno<-glm(fracture~premeno,family=binomial(link='logit'),data=train)
simple.log.premeno.sum = summary(simple.log.premeno)

# momfrac
simple.log.momfrac<-glm(fracture~momfrac,family=binomial(link='logit'),data=train)
simple.log.momfrac.sum = summary(simple.log.momfrac)

# armassist
simple.log.armassist<-glm(fracture~armassist,family=binomial(link='logit'),data=train)
simple.log.armassist.sum = summary(simple.log.armassist)

# smoke
simple.log.smoke<-glm(fracture~smoke,family=binomial(link='logit'),data=train)
simple.log.smoke.sum = summary(simple.log.smoke)

# raterisk
simple.log.raterisk<-glm(fracture~raterisk,family=binomial(link='logit'),data=train)
simple.log.raterisk.sum = summary(simple.log.raterisk)

# fracscore
simple.log.fracscore<-glm(fracture~fracscore,family=binomial(link='logit'),data=train)
simple.log.fracscore.sum = summary(simple.log.fracscore)

# bonemed
simple.log.bonemed<-glm(fracture~bonemed,family=binomial(link='logit'),data=train)
simple.log.bonemed.sum = summary(simple.log.bonemed)

# bonemed_fu
simple.log.bonemed_fu<-glm(fracture~bonemed_fu,family=binomial(link='logit'),data=train)
simple.log.bonemed_fu.sum = summary(simple.log.bonemed_fu)

# bonetreat
simple.log.bonetreat<-glm(fracture~bonetreat,family=binomial(link='logit'),data=train)
simple.log.bonetreat.sum = summary(simple.log.bonetreat)
```

``` r
# priorfrac
simple.log.priorfrac.sum$coefficients
```

    ##              Estimate Std. Error z value                      Pr(>|z|)
    ## (Intercept)    -1.426     0.1414 -10.081 0.000000000000000000000006714
    ## priorfracYes    1.060     0.2437   4.349 0.000013666179481153440988849

``` r
# age
simple.log.age.sum$coefficients
```

    ##             Estimate Std. Error z value    Pr(>|z|)
    ## (Intercept) -4.05002    0.89623  -4.519 0.000006215
    ## age          0.04231    0.01267   3.340 0.000838153

``` r
# weight
simple.log.weight.sum$coefficients
```

    ##              Estimate Std. Error z value Pr(>|z|)
    ## (Intercept) -0.683674   0.509143 -1.3428   0.1793
    ## weight      -0.006016   0.006976 -0.8623   0.3885

``` r
# height
simple.log.height.sum$coefficients
```

    ##             Estimate Std. Error z value Pr(>|z|)
    ## (Intercept)  4.71146    2.99392   1.574  0.11556
    ## height      -0.03618    0.01862  -1.943  0.05201

``` r
# bmi
simple.log.bmi.sum$coefficients
```

    ##              Estimate Std. Error z value Pr(>|z|)
    ## (Intercept) -1.028066    0.53698 -1.9145  0.05555
    ## bmi         -0.003128    0.01905 -0.1642  0.86959

``` r
# premeno
simple.log.premeno.sum$coefficients
```

    ##             Estimate Std. Error z value                  Pr(>|z|)
    ## (Intercept)  -1.1585     0.1274 -9.0960 0.00000000000000000009368
    ## premenoYes    0.2094     0.2721  0.7698 0.44142454587949042643658

``` r
# momfrac
simple.log.momfrac.sum$coefficients
```

    ##             Estimate Std. Error z value                    Pr(>|z|)
    ## (Intercept)  -1.1901     0.1224  -9.720 0.0000000000000000000002478
    ## momfracYes    0.5541     0.3162   1.753 0.0796758824434630613664510

``` r
# armassist
simple.log.armassist.sum$coefficients
```

    ##              Estimate Std. Error z value                  Pr(>|z|)
    ## (Intercept)   -1.4100     0.1547  -9.116 0.00000000000000000007797
    ## armassistYes   0.7075     0.2283   3.098 0.00194505752636372152529

``` r
# smoke
simple.log.smoke.sum$coefficients
```

    ##             Estimate Std. Error z value                   Pr(>|z|)
    ## (Intercept)  -1.1054     0.1165 -9.4848 0.000000000000000000002429
    ## smokeYes     -0.1268     0.4451 -0.2848 0.775814948711020968552532

``` r
# raterisk
simple.log.raterisk.sum$coefficients
```

    ##                 Estimate Std. Error z value         Pr(>|z|)
    ## (Intercept)      -1.6603     0.2275  -7.299 0.00000000000029
    ## rateriskSame      0.5443     0.2944   1.849 0.06449595109024
    ## rateriskGreater   1.0256     0.2941   3.487 0.00048859486001

``` r
# fracscore
simple.log.fracscore.sum$coefficients
```

    ##             Estimate Std. Error z value               Pr(>|z|)
    ## (Intercept)  -2.0098    0.23573  -8.526 0.00000000000000001518
    ## fracscore     0.2237    0.04786   4.674 0.00000294956575793599

``` r
# bonemed
simple.log.bonemed.sum$coefficients
```

    ##             Estimate Std. Error z value                  Pr(>|z|)
    ## (Intercept)  -1.2609     0.1354  -9.312 0.00000000000000000001254
    ## bonemedYes    0.5258     0.2462   2.136 0.03269847370748939913154

``` r
# bonemed_fu
simple.log.bonemed_fu.sum$coefficients
```

    ##               Estimate Std. Error z value                   Pr(>|z|)
    ## (Intercept)    -1.3383     0.1404  -9.529 0.000000000000000000001584
    ## bonemed_fuYes   0.7211     0.2393   3.013 0.002584333250422946558211

``` r
# bonetreat
simple.log.bonetreat.sum$coefficients
```

    ##              Estimate Std. Error z value                  Pr(>|z|)
    ## (Intercept)   -1.2119     0.1315  -9.214 0.00000000000000000003148
    ## bonetreatYes   0.3936     0.2556   1.540 0.12358095157835313759875

### Observations:

-   We can see that the following variables are statistically
    significant since p-value \< 0.05:

1.  priorfracYes
2.  age
3.  height
4.  armassistYes
5.  rateriskSame
6.  rateriskGreater
7.  fracscore
8.  bonemedYes
9.  bonemed_fuYes
10. bonetreatYes

This result is in-line with what we have observed through EDA for the
continuous variables. Next, let’s fit all the variables and observe the
effect and see how it is changing the significance of the predictors.

## Fit all variables at the same time to check effects

``` r
multi_var.log<-glm(fracture~priorfrac+age+weight+height+bmi+premeno+momfrac+armassist+smoke+raterisk+fracscore+bonemed+
                     bonemed_fu+bonetreat,family=binomial(link='logit'),data=train)
multi_var.log.sum = summary(multi_var.log)
multi_var.log.sum$coefficients
```

    ##                   Estimate Std. Error  z value Pr(>|z|)
    ## (Intercept)     -18.133928   13.63463 -1.32999 0.183522
    ## priorfracYes      0.693635    0.43776  1.58451 0.113078
    ## age               0.024911    0.06172  0.40360 0.686508
    ## weight           -0.137226    0.09619 -1.42657 0.153703
    ## height            0.087262    0.08622  1.01210 0.311488
    ## bmi               0.362430    0.25094  1.44428 0.148661
    ## premenoYes        0.253407    0.30653  0.82669 0.408416
    ## momfracYes        0.607017    0.46222  1.31327 0.189093
    ## armassistYes      0.436793    0.66987  0.65205 0.514366
    ## smokeYes         -0.183912    0.58058 -0.31677 0.751416
    ## rateriskSame      0.363463    0.31503  1.15373 0.248609
    ## rateriskGreater   0.674468    0.33738  1.99911 0.045597
    ## fracscore         0.008364    0.30805  0.02715 0.978340
    ## bonemedYes        1.332603    0.68927  1.93336 0.053192
    ## bonemed_fuYes     1.292946    0.52379  2.46845 0.013570
    ## bonetreatYes     -2.399374    0.87892 -2.72990 0.006335

### Check VIF

``` r
vif(multi_var.log)
```

    ##               GVIF Df GVIF^(1/(2*Df))
    ## priorfrac    2.909  1           1.705
    ## age         21.094  1           4.593
    ## weight     170.899  1          13.073
    ## height      20.555  1           4.534
    ## bmi        157.762  1          12.560
    ## premeno      1.132  1           1.064
    ## momfrac      1.845  1           1.358
    ## armassist    7.611  1           2.759
    ## smoke        1.516  1           1.231
    ## raterisk     1.247  2           1.057
    ## fracscore   38.554  1           6.209
    ## bonemed      6.926  1           2.632
    ## bonemed_fu   4.188  1           2.047
    ## bonetreat   10.521  1           3.244

### Remove the multicollinearity

``` r
multi_var.log.vif<-glm(fracture~priorfrac+age+weight+height+premeno+momfrac+armassist+smoke+raterisk+fracscore+bonemed+
                     bonemed_fu+bonetreat,family=binomial(link='logit'),data=train)
multi_var.log.vif.sum = summary(multi_var.log.vif)
multi_var.log.vif.sum$coefficients
```

    ##                 Estimate Std. Error  z value Pr(>|z|)
    ## (Intercept)      0.61028   4.423815  0.13795 0.890278
    ## priorfracYes     0.76995   0.432076  1.78198 0.074752
    ## age              0.04351   0.060262  0.72198 0.470306
    ## weight           0.00061   0.009787  0.06233 0.950303
    ## height          -0.03463   0.021349 -1.62233 0.104733
    ## premenoYes       0.29569   0.304591  0.97078 0.331660
    ## momfracYes       0.70836   0.457781  1.54737 0.121774
    ## armassistYes     0.62487   0.654863  0.95419 0.339985
    ## smokeYes        -0.10380   0.577891 -0.17963 0.857447
    ## rateriskSame     0.34465   0.314249  1.09673 0.272759
    ## rateriskGreater  0.65612   0.335660  1.95471 0.050617
    ## fracscore       -0.09076   0.299061 -0.30347 0.761531
    ## bonemedYes       1.34997   0.684311  1.97275 0.048524
    ## bonemed_fuYes    1.28949   0.524784  2.45717 0.014003
    ## bonetreatYes    -2.41404   0.875252 -2.75811 0.005814

### re-Check VIF

``` r
vif(multi_var.log.vif)
```

    ##              GVIF Df GVIF^(1/(2*Df))
    ## priorfrac   2.850  1           1.688
    ## age        20.252  1           4.500
    ## weight      1.842  1           1.357
    ## height      1.222  1           1.106
    ## premeno     1.122  1           1.059
    ## momfrac     1.795  1           1.340
    ## armassist   7.322  1           2.706
    ## smoke       1.506  1           1.227
    ## raterisk    1.241  2           1.055
    ## fracscore  36.550  1           6.046
    ## bonemed     6.860  1           2.619
    ## bonemed_fu  4.219  1           2.054
    ## bonetreat  10.474  1           3.236

### Observations:

-   Fitting all variables to the logistic regression model, it shows
    that only ‘bonmed’, ‘bonemed_fu’ and ‘bonetreat’ are statistically
    significant.
-   We observed multicollinearity between ‘weight’ and ‘bmi’. We removed
    the ‘bmi’ variable and re-run the logistic regression model.
-   As a result the following variables are statistically significant:

1.  height
2.  bonemed
3.  bonemed_fu
4.  bonetreat

Let’s test other feature selection methods as well.

## Stepwise selection

``` r
bonemed_df.step = train[,c('priorfrac', 'age', 'weight', 'height', 'bmi', 'premeno', 'momfrac', 'armassist', 'smoke', 'raterisk', 'fracscore', 
                                'bonemed', 'bonemed_fu', 'bonetreat', 'fracture')]
step.full.log = glm(fracture~.,family=binomial(link='logit'),data=bonemed_df.step)
step.log = step.full.log %>% stepAIC(trace=FALSE)
```

``` r
summary(step.log)
```

    ## 
    ## Call:
    ## glm(formula = fracture ~ priorfrac + age + weight + bmi + momfrac + 
    ##     armassist + raterisk + bonemed + bonemed_fu + bonetreat, 
    ##     family = binomial(link = "logit"), data = bonemed_df.step)
    ## 
    ## Deviance Residuals: 
    ##    Min      1Q  Median      3Q     Max  
    ## -1.526  -0.726  -0.555  -0.346   2.254  
    ## 
    ## Coefficients:
    ##                 Estimate Std. Error z value Pr(>|z|)   
    ## (Intercept)      -3.9663     1.4369   -2.76   0.0058 **
    ## priorfracYes      0.6735     0.2721    2.48   0.0133 * 
    ## age               0.0238     0.0155    1.54   0.1245   
    ## weight           -0.0442     0.0223   -1.98   0.0479 * 
    ## bmi               0.1217     0.0606    2.01   0.0444 * 
    ## momfracYes        0.6283     0.3486    1.80   0.0715 . 
    ## armassistYes      0.4734     0.2741    1.73   0.0842 . 
    ## rateriskSame      0.3758     0.3123    1.20   0.2289   
    ## rateriskGreater   0.6845     0.3360    2.04   0.0416 * 
    ## bonemedYes        1.3057     0.6881    1.90   0.0578 . 
    ## bonemed_fuYes     1.2774     0.5208    2.45   0.0142 * 
    ## bonetreatYes     -2.3827     0.8787   -2.71   0.0067 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for binomial family taken to be 1)
    ## 
    ##     Null deviance: 475.22  on 424  degrees of freedom
    ## Residual deviance: 423.79  on 413  degrees of freedom
    ## AIC: 447.8
    ## 
    ## Number of Fisher Scoring iterations: 4

``` r
exp(cbind("Odds ratio" = coef(step.log), confint.default(step.log, level = 0.95)))
```

    ##                 Odds ratio    2.5 %  97.5 %
    ## (Intercept)        0.01894 0.001133  0.3167
    ## priorfracYes       1.96112 1.150618  3.3425
    ## age                1.02410 0.993453  1.0557
    ## weight             0.95677 0.915787  0.9996
    ## bmi                1.12946 1.003051  1.2718
    ## momfracYes         1.87435 0.946554  3.7116
    ## armassistYes       1.60538 0.938085  2.7473
    ## rateriskSame       1.45613 0.789516  2.6856
    ## rateriskGreater    1.98272 1.026292  3.8305
    ## bonemedYes         3.69013 0.957843 14.2164
    ## bonemed_fuYes      3.58715 1.292440  9.9561
    ## bonetreatYes       0.09230 0.016492  0.5166

``` r
vif(step.log)
```

    ##              GVIF Df GVIF^(1/(2*Df))
    ## priorfrac   1.129  1           1.062
    ## age         1.351  1           1.162
    ## weight      9.658  1           3.108
    ## bmi         9.252  1           3.042
    ## momfrac     1.042  1           1.021
    ## armassist   1.282  1           1.132
    ## raterisk    1.229  2           1.053
    ## bonemed     6.924  1           2.631
    ## bonemed_fu  4.155  1           2.038
    ## bonetreat  10.550  1           3.248

## Obseravtions:

-   Running a stepwise selection to identify the predictors the process
    is selecting:
-   priorfrac
-   age
-   weight
-   bmi
-   momfrac
-   armassist
-   raterisk
-   bonemed
-   bonemed_fu
-   bonetreat

## PCA

``` r
# Let's use PCA to see if the continuous variables separate or not

num_cols = train %>% dplyr::select(where(is.numeric)) %>% colnames()
pca_var = num_cols[c(-6)]
pca_df = train[pca_var]

pc.result=prcomp(pca_df,scale.=TRUE)
pc.scores=pc.result$x
pc.scores=data.frame(pc.scores)
pc.scores$fracture=train$fracture

#plot the first few pc's
ggplot(data = pc.scores, aes(x = PC1, y = PC2)) +
  geom_point(aes(col=fracture), size=1)+
  theme(plot.title = element_text(hjust = 0.5))+
  ggtitle("PC1 and PC2 on Fracture")
```

![](MSDS_6372_Project2_files/figure-gfm/PCA-1.png)<!-- -->

``` r
# Let's check lower PCs
ggplot(data = pc.scores, aes(x = PC1, y = PC3)) +
  geom_point(aes(col=fracture), size=1)+
  theme(plot.title = element_text(hjust = 0.5))+
  ggtitle("PC1 and PC3 on Fracture")
```

![](MSDS_6372_Project2_files/figure-gfm/PCA-2.png)<!-- -->

``` r
ggplot(data = pc.scores, aes(x = PC1, y = PC4)) +
  geom_point(aes(col=fracture), size=1)+
  theme(plot.title = element_text(hjust = 0.5))+
  ggtitle("PC1 and PC4 on Fracture")
```

![](MSDS_6372_Project2_files/figure-gfm/PCA-3.png)<!-- -->

``` r
ggplot(data = pc.scores, aes(x = PC1, y = PC5)) +
  geom_point(aes(col=fracture), size=1)+
  theme(plot.title = element_text(hjust = 0.5))+
  ggtitle("PC1 and PC5 on Fracture")
```

![](MSDS_6372_Project2_files/figure-gfm/PCA-4.png)<!-- -->

``` r
ggplot(data = pc.scores, aes(x = PC2, y = PC3)) +
  geom_point(aes(col=fracture), size=1)+
  theme(plot.title = element_text(hjust = 0.5))+
  ggtitle("PC2 and PC3 on Fracture")
```

![](MSDS_6372_Project2_files/figure-gfm/PCA-5.png)<!-- -->

``` r
ggplot(data = pc.scores, aes(x = PC2, y = PC4)) +
  geom_point(aes(col=fracture), size=1)+
  theme(plot.title = element_text(hjust = 0.5))+
  ggtitle("PC2 and PC4 on Fracture")
```

![](MSDS_6372_Project2_files/figure-gfm/PCA-6.png)<!-- -->

``` r
ggplot(data = pc.scores, aes(x = PC2, y = PC5)) +
  geom_point(aes(col=fracture), size=1)+
  theme(plot.title = element_text(hjust = 0.5))+
  ggtitle("PC2 and PC5 on Fracture")
```

![](MSDS_6372_Project2_files/figure-gfm/PCA-7.png)<!-- -->

``` r
ggplot(data = pc.scores, aes(x = PC3, y = PC4)) +
  geom_point(aes(col=fracture), size=1)+
  theme(plot.title = element_text(hjust = 0.5))+
  ggtitle("PC3 and PC4 on Fracture")
```

![](MSDS_6372_Project2_files/figure-gfm/PCA-8.png)<!-- -->

``` r
ggplot(data = pc.scores, aes(x = PC3, y = PC5)) +
  geom_point(aes(col=fracture), size=1)+
  theme(plot.title = element_text(hjust = 0.5))+
  ggtitle("PC3 and PC5 on Fracture")
```

![](MSDS_6372_Project2_files/figure-gfm/PCA-9.png)<!-- -->

``` r
ggplot(data = pc.scores, aes(x = PC4, y = PC5)) +
  geom_point(aes(col=fracture), size=1)+
  theme(plot.title = element_text(hjust = 0.5))+
  ggtitle("PC4 and PC5 on Fracture")
```

![](MSDS_6372_Project2_files/figure-gfm/PCA-10.png)<!-- -->

``` r
par(mfrow=c(1,2))
eigenvals = (pc.result$sdev)^2
plot(1:5,eigenvals/sum(eigenvals),type="l",main="Scree Plot",ylab="Prop. Var. Explained")
cumulative.prop = cumsum(eigenvals/sum(eigenvals))
plot(1:5,cumulative.prop,type="l",main="Cumulative proportion",ylim=c(0,1))
```

![](MSDS_6372_Project2_files/figure-gfm/PCA-11.png)<!-- -->

``` r
for (i in colnames(pc.scores))
{
loess_pca = pc.scores %>% 
ggplot(aes(x=.data[[i]],y=train$fracture.num))+
geom_point()+
geom_smooth(formula = y ~ x, method="loess")+
theme(plot.title = element_text(hjust = 0.5)) +
ggtitle(paste(i, " test"))
egg::ggarrange(loess_pca, ncol=2)
}
```

![](MSDS_6372_Project2_files/figure-gfm/loess%20for%20PCI-1.png)<!-- -->![](MSDS_6372_Project2_files/figure-gfm/loess%20for%20PCI-2.png)<!-- -->![](MSDS_6372_Project2_files/figure-gfm/loess%20for%20PCI-3.png)<!-- -->![](MSDS_6372_Project2_files/figure-gfm/loess%20for%20PCI-4.png)<!-- -->![](MSDS_6372_Project2_files/figure-gfm/loess%20for%20PCI-5.png)<!-- -->![](MSDS_6372_Project2_files/figure-gfm/loess%20for%20PCI-6.png)<!-- -->

## Observations:

-   The levels in the PCA graph are completely intermixed. There is no
    obvious separation in the PCA.
-   The Scree Plot shows that if we include PC1 that will explain \~50%
    of the variation of the dataset
-   The cumulative proportion plot helps us see that PC1, PC2, PC3
    combined explains the variation \~98% but if we add PC4 than 100% of
    the variation is explained.
-   This lets us conclude that we need ‘age’, ‘weight’, ‘height’ and
    ‘bmi’ to explain 100% of the variation in the data.
-   The loess plot of PCA shows that PC3 is the only continous variable
    which could be important for predicting fracture in the first year.

## Penalized logistic regression (LASSO)

``` r
# Dummy code categorical predictor variables
x = model.matrix(fracture~., train)[,-c(1,17)]

# Convert the outcome (class) to a numerical variable
y <- train$fracture.num
```

``` r
grid=10^seq(10,-2, length =100)
lasso.mod = glmnet(x, y, family = binomial(link='logit'), alpha=1, lambda=grid)
cv.out=cv.glmnet(x,y,family = "binomial", alpha=1)
plot(cv.out)
```

![](MSDS_6372_Project2_files/figure-gfm/lasso%20feature%20selection-1.png)<!-- -->

``` r
bestlambda = cv.out$lambda.min
coef(lasso.mod,s=bestlambda)
```

    ## 16 x 1 sparse Matrix of class "dgCMatrix"
    ##                        s1
    ## (Intercept)     -1.391038
    ## priorfracYes     0.385785
    ## age              .       
    ## weight           .       
    ## height          -0.002597
    ## bmi              .       
    ## premenoYes       .       
    ## momfracYes       0.049417
    ## armassistYes     0.027311
    ## smokeYes         .       
    ## rateriskSame     .       
    ## rateriskGreater  0.254864
    ## fracscore        0.111860
    ## bonemedYes       .       
    ## bonemed_fuYes    0.169981
    ## bonetreatYes     .

``` r
# Final LASSO model
lasso.mod.final=glmnet(x,y,family = binomial(link='logit'), alpha=1, type.measure = "class", lambda=bestlambda)
```

##################################################################################### 

# Predictions

##################################################################################### 

``` r
# VIF score function
vif_score = function(model_name)
{
  vif_values = vif(model_name)[,3]
  barplot(vif_values, main = 'VIF Values', horiz = TRUE, col="blue", xlim = c(0,12))
  abline(v=10, col="red")
  return (vif_values)
}
```

``` r
y.test = model.matrix(fracture~., test)[,-c(1,17)]
```

##################################################################################### 

# LASSO

##################################################################################### 

``` r
cutoff.lasso<-0.43
# Predicting on the final LASSO model TRAINing data
lasso.mod.final.pred.train <- predict(lasso.mod.final, newx = x, type = "response")

# Confusion matrix
class.lasso.final.train<-factor(ifelse(lasso.mod.final.pred.train>cutoff.lasso,"Yes","No"),levels=c("No","Yes"))

#Confusion Matrix for LASSO
print("Confusion matrix for LASSO TRAINING with 0.5 cutoff")
```

    ## [1] "Confusion matrix for LASSO TRAINING with 0.5 cutoff"

``` r
confusionMatrix(table(class.lasso.final.train,train$fracture))
```

    ## Confusion Matrix and Statistics
    ## 
    ##                        
    ## class.lasso.final.train  No Yes
    ##                     No  309  95
    ##                     Yes  11  10
    ##                                               
    ##                Accuracy : 0.751               
    ##                  95% CI : (0.707, 0.791)      
    ##     No Information Rate : 0.753               
    ##     P-Value [Acc > NIR] : 0.571               
    ##                                               
    ##                   Kappa : 0.083               
    ##                                               
    ##  Mcnemar's Test P-Value : 0.000000000000000753
    ##                                               
    ##             Sensitivity : 0.9656              
    ##             Specificity : 0.0952              
    ##          Pos Pred Value : 0.7649              
    ##          Neg Pred Value : 0.4762              
    ##              Prevalence : 0.7529              
    ##          Detection Rate : 0.7271              
    ##    Detection Prevalence : 0.9506              
    ##       Balanced Accuracy : 0.5304              
    ##                                               
    ##        'Positive' Class : No                  
    ## 

``` r
################ Misclassification rate train ######################################
#cross.table.lasso.train = table(table(class.lasso.final.train,train$fracture))
#MCR_lasso.train = (cross.table.lasso.train[2]+cross.table.lasso.train[3])/dim(train)[1]
misClasificError.lasso.train = mean(class.lasso.final.train != train$fracture)
print(paste('Misclassification Rate for lasso on training set: ', misClasificError.lasso.train))
```

    ## [1] "Misclassification Rate for lasso on training set:  0.249411764705882"

``` r
#####################################################################################################
# Predicting on the final LASSO model TEST data
lasso.mod.final.pred.test <- predict(lasso.mod.final, newx = y.test, type = "response")

# Confusion matrix
class.lasso.final.test<-factor(ifelse(lasso.mod.final.pred.test>cutoff.lasso,"Yes","No"),levels=c("No","Yes"))

#Confusion Matrix for LASSO
print("Confusion matrix for LASSO TEST with 0.5 cutoff")
```

    ## [1] "Confusion matrix for LASSO TEST with 0.5 cutoff"

``` r
confusionMatrix(table(class.lasso.final.test,test$fracture))
```

    ## Confusion Matrix and Statistics
    ## 
    ##                       
    ## class.lasso.final.test No Yes
    ##                    No  53  19
    ##                    Yes  2   1
    ##                                         
    ##                Accuracy : 0.72          
    ##                  95% CI : (0.604, 0.818)
    ##     No Information Rate : 0.733         
    ##     P-Value [Acc > NIR] : 0.65856       
    ##                                         
    ##                   Kappa : 0.019         
    ##                                         
    ##  Mcnemar's Test P-Value : 0.00048       
    ##                                         
    ##             Sensitivity : 0.964         
    ##             Specificity : 0.050         
    ##          Pos Pred Value : 0.736         
    ##          Neg Pred Value : 0.333         
    ##              Prevalence : 0.733         
    ##          Detection Rate : 0.707         
    ##    Detection Prevalence : 0.960         
    ##       Balanced Accuracy : 0.507         
    ##                                         
    ##        'Positive' Class : No            
    ## 

``` r
################ Misclassification rate test ######################################
#cross.table.lasso.test = table(class.lasso.final.test,test$fracture)
#MCR_lasso.test = (cross.table.lasso.test[2]+cross.table.lasso.test[3])/dim(test)[1]
#print(paste('Misclassification Rate for LASSO selection on test set: ', MCR_lasso.test))
misClasificError.lasso.test = mean(class.lasso.final.test != test$fracture)
print(paste('Misclassification Rate for lasso on test set: ', misClasificError.lasso.test))
```

    ## [1] "Misclassification Rate for lasso on test set:  0.28"

``` r
################ LASSO ROC Curve ######################################
library(ROCR)
lasso.roc.test = prediction(lasso.mod.final.pred.test, test$fracture,label.ordering=c("No","Yes"))
roc.lasso.test = performance(lasso.roc.test, measure = "tpr", x.measure = "fpr")
plot(roc.lasso.test,colorize = TRUE)
abline(a=0, b= 1)
```

![](MSDS_6372_Project2_files/figure-gfm/LASSO%20Predict-1.png)<!-- -->

## Assumprions

``` r
# Checking logistic regression model assumptions

# The Observations are Independent
# We can assume that the observations are not coming from repeated measures

#There is No Multicollinearity Among Explanatory Variables
### Visualize VIF
#vif_score(lasso.mod)


# There are No Extreme Outliers (cooks D)
#influenceIndexPlot(lasso.mod)

# There is a Linear Relationship Between Explanatory Variables and the Logit of the Response Variable

# Influential point analysis and residual plots
#residualPlots(lasso.mod)
#influenceIndexPlot(lasso.mod)
#influencePlot(lasso.mod)
```

## Observations:

## All assumptions are met:

-   **The Response Variable is Binary**: ‘Fracture’ as a response
    variable is a factor with binary levels (Yes/No)
-   **Independence**: We can assume that the observations are
    independent
-   **Multicolliearity**: The is no multicollinearity among the
    explanatory variables (VIF values show no multicollinearity)
-   **Outliers**: The largest Cooks D value is 0.015 which indicates
    that there is no extreme outlier
-   **The Sample Size is Sufficiently Large**: ?

## Model scoring

``` r
## Add misclassification rate
```

##################################################################################### 

# Stepwise

##################################################################################### 

``` r
cutoff.step = 0.42
# Predicting on the final Stepwise model TRAINing data
stepwise.mod.final.pred.train <- predict(step.log, newdata = train, type = "response")

# Confusion matrix
class.stepwise.final.train<-factor(ifelse(stepwise.mod.final.pred.train>cutoff.step,"Yes","No"),levels=c("No","Yes"))

#Confusion Matrix for Stepwise
print("Confusion matrix for Stepwise TRAINING with 0.5 cutoff")
```

    ## [1] "Confusion matrix for Stepwise TRAINING with 0.5 cutoff"

``` r
confusionMatrix(table(class.stepwise.final.train,train$fracture))
```

    ## Confusion Matrix and Statistics
    ## 
    ##                           
    ## class.stepwise.final.train  No Yes
    ##                        No  293  73
    ##                        Yes  27  32
    ##                                         
    ##                Accuracy : 0.765         
    ##                  95% CI : (0.721, 0.804)
    ##     No Information Rate : 0.753         
    ##     P-Value [Acc > NIR] : 0.309         
    ##                                         
    ##                   Kappa : 0.258         
    ##                                         
    ##  Mcnemar's Test P-Value : 0.0000068     
    ##                                         
    ##             Sensitivity : 0.916         
    ##             Specificity : 0.305         
    ##          Pos Pred Value : 0.801         
    ##          Neg Pred Value : 0.542         
    ##              Prevalence : 0.753         
    ##          Detection Rate : 0.689         
    ##    Detection Prevalence : 0.861         
    ##       Balanced Accuracy : 0.610         
    ##                                         
    ##        'Positive' Class : No            
    ## 

``` r
################ Misclassification rate train ######################################
#cross.table.stepwise.train = table(class.stepwise.final.train,train$fracture)
#MCR_stepwise.train = (cross.table.stepwise.train[2]+cross.table.stepwise.train[3])/dim(train)[1]
#print(paste('Misclassification Rate for stepwise selection on training set: ', MCR_stepwise.train))

misClasificError.stepwise.train = mean(class.stepwise.final.train != train$fracture)
print(paste('Misclassification Rate for stepwise on train set: ', misClasificError.stepwise.train))
```

    ## [1] "Misclassification Rate for stepwise on train set:  0.235294117647059"

``` r
#####################################################################################################
# Predicting on the final Stepwise model TEST data
stepwise.mod.final.pred.test <- predict(step.log, newdata = test, type = "response")

# Confusion matrix
class.stepwise.final.test<-factor(ifelse(stepwise.mod.final.pred.test>cutoff.step,"Yes","No"),levels=c("No","Yes"))

#Confusion Matrix for Stepwise
print("Confusion matrix for Stepwise TEST with 0.5 cutoff")
```

    ## [1] "Confusion matrix for Stepwise TEST with 0.5 cutoff"

``` r
confusionMatrix(table(class.stepwise.final.test,test$fracture))
```

    ## Confusion Matrix and Statistics
    ## 
    ##                          
    ## class.stepwise.final.test No Yes
    ##                       No  52   9
    ##                       Yes  3  11
    ##                                         
    ##                Accuracy : 0.84          
    ##                  95% CI : (0.737, 0.914)
    ##     No Information Rate : 0.733         
    ##     P-Value [Acc > NIR] : 0.0211        
    ##                                         
    ##                   Kappa : 0.548         
    ##                                         
    ##  Mcnemar's Test P-Value : 0.1489        
    ##                                         
    ##             Sensitivity : 0.945         
    ##             Specificity : 0.550         
    ##          Pos Pred Value : 0.852         
    ##          Neg Pred Value : 0.786         
    ##              Prevalence : 0.733         
    ##          Detection Rate : 0.693         
    ##    Detection Prevalence : 0.813         
    ##       Balanced Accuracy : 0.748         
    ##                                         
    ##        'Positive' Class : No            
    ## 

``` r
################ Misclassification rate test ######################################
#cross.table.stepwise.test = table(class.stepwise.final.test,test$fracture)
#MCR_stepwise.test = (cross.table.stepwise.test[2]+cross.table.stepwise.test[3])/dim(test)[1]
#print(paste('Misclassification Rate for stepwise selection on test set: ', MCR_stepwise.test))

misClasificError.stepwise.test = mean(class.stepwise.final.test != test$fracture)
print(paste('Misclassification Rate for stepwise on test set: ', misClasificError.stepwise.test))
```

    ## [1] "Misclassification Rate for stepwise on test set:  0.16"

``` r
# ODD ratios for interpretation
summary(step.log)
```

    ## 
    ## Call:
    ## glm(formula = fracture ~ priorfrac + age + weight + bmi + momfrac + 
    ##     armassist + raterisk + bonemed + bonemed_fu + bonetreat, 
    ##     family = binomial(link = "logit"), data = bonemed_df.step)
    ## 
    ## Deviance Residuals: 
    ##    Min      1Q  Median      3Q     Max  
    ## -1.526  -0.726  -0.555  -0.346   2.254  
    ## 
    ## Coefficients:
    ##                 Estimate Std. Error z value Pr(>|z|)   
    ## (Intercept)      -3.9663     1.4369   -2.76   0.0058 **
    ## priorfracYes      0.6735     0.2721    2.48   0.0133 * 
    ## age               0.0238     0.0155    1.54   0.1245   
    ## weight           -0.0442     0.0223   -1.98   0.0479 * 
    ## bmi               0.1217     0.0606    2.01   0.0444 * 
    ## momfracYes        0.6283     0.3486    1.80   0.0715 . 
    ## armassistYes      0.4734     0.2741    1.73   0.0842 . 
    ## rateriskSame      0.3758     0.3123    1.20   0.2289   
    ## rateriskGreater   0.6845     0.3360    2.04   0.0416 * 
    ## bonemedYes        1.3057     0.6881    1.90   0.0578 . 
    ## bonemed_fuYes     1.2774     0.5208    2.45   0.0142 * 
    ## bonetreatYes     -2.3827     0.8787   -2.71   0.0067 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for binomial family taken to be 1)
    ## 
    ##     Null deviance: 475.22  on 424  degrees of freedom
    ## Residual deviance: 423.79  on 413  degrees of freedom
    ## AIC: 447.8
    ## 
    ## Number of Fisher Scoring iterations: 4

``` r
exp(cbind("Odds ratio" = coef(step.log), confint.default(step.log, level = 0.95)))
```

    ##                 Odds ratio    2.5 %  97.5 %
    ## (Intercept)        0.01894 0.001133  0.3167
    ## priorfracYes       1.96112 1.150618  3.3425
    ## age                1.02410 0.993453  1.0557
    ## weight             0.95677 0.915787  0.9996
    ## bmi                1.12946 1.003051  1.2718
    ## momfracYes         1.87435 0.946554  3.7116
    ## armassistYes       1.60538 0.938085  2.7473
    ## rateriskSame       1.45613 0.789516  2.6856
    ## rateriskGreater    1.98272 1.026292  3.8305
    ## bonemedYes         3.69013 0.957843 14.2164
    ## bonemed_fuYes      3.58715 1.292440  9.9561
    ## bonetreatYes       0.09230 0.016492  0.5166

``` r
################ Stepwise ROC Curve ######################################
library(ROCR)
stepwise.roc.test = prediction(stepwise.mod.final.pred.test, test$fracture,label.ordering=c("No","Yes"))
roc.stepwise.test = performance(stepwise.roc.test, measure = "tpr", x.measure = "fpr")
plot(roc.stepwise.test,colorize = TRUE)
abline(a=0, b= 1)
```

![](MSDS_6372_Project2_files/figure-gfm/Stepwise%20Predict-1.png)<!-- -->

##################################################################################### 

# Manual/Intuition

##################################################################################### 

``` r
# 5-fold cross validation
cv <- trainControl(
  method = "repeatedcv", 
  number = 5,
  repeats = 10,
  savePredictions = TRUE,
  summaryFunction=mnLogLoss,
  classProbs = TRUE
)

MLogReg = train(
  fracture ~ age + height + bmi + priorfrac,
  data = train,
  method = "glm",
  family = "binomial",
  trControl = cv,
  metric = "logLoss")

### Visualize VIF
MLogReg_VIF = vif(MLogReg$finalModel)
barplot(MLogReg_VIF, main = 'VIF Values (Custom Logistic Regression', horiz = TRUE, col="blue", xlim = c(0,12))
abline(v=10, col="red")
```

![](MSDS_6372_Project2_files/figure-gfm/manual%20selection%20prediction-1.png)<!-- -->

``` r
### Hypothesis testing
summary(MLogReg$finalModel)
```

    ## 
    ## Call:
    ## NULL
    ## 
    ## Deviance Residuals: 
    ##    Min      1Q  Median      3Q     Max  
    ## -1.272  -0.723  -0.622  -0.477   2.063  
    ## 
    ## Coefficients:
    ##              Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)    0.1120     3.4043    0.03  0.97376    
    ## age            0.0301     0.0140    2.16  0.03086 *  
    ## height        -0.0242     0.0190   -1.27  0.20418    
    ## bmi            0.0116     0.0203    0.57  0.56767    
    ## priorfracYes   0.8789     0.2546    3.45  0.00056 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for binomial family taken to be 1)
    ## 
    ##     Null deviance: 475.22  on 424  degrees of freedom
    ## Residual deviance: 449.62  on 420  degrees of freedom
    ## AIC: 459.6
    ## 
    ## Number of Fisher Scoring iterations: 4

``` r
anova(MLogReg$finalModel, test="Chisq")
```

    ## Analysis of Deviance Table
    ## 
    ## Model: binomial, link: logit
    ## 
    ## Response: .outcome
    ## 
    ## Terms added sequentially (first to last)
    ## 
    ## 
    ##              Df Deviance Resid. Df Resid. Dev Pr(>Chi)    
    ## NULL                           424        475             
    ## age           1    11.30       423        464  0.00078 ***
    ## height        1     2.13       422        462  0.14471    
    ## bmi           1     0.53       421        461  0.46852    
    ## priorfracYes  1    11.65       420        450  0.00064 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
library(ResourceSelection) 
```

    ## ResourceSelection 0.3-5   2019-07-22

``` r
hoslem.test(MLogReg$finalModel$y,fitted(MLogReg))
```

    ## 
    ##  Hosmer and Lemeshow goodness of fit (GOF) test
    ## 
    ## data:  MLogReg$finalModel$y, fitted(MLogReg)
    ## X-squared = 12, df = 8, p-value = 0.1

``` r
# Predicting
MLogReg.pred.train = predict(MLogReg, train)
MLogReg.pred.test = predict(MLogReg, test)

#Confusion Matrix for manual
confusionMatrix(table(MLogReg.pred.train,train$fracture))
```

    ## Confusion Matrix and Statistics
    ## 
    ##                   
    ## MLogReg.pred.train  No Yes
    ##                No  311  98
    ##                Yes   9   7
    ##                                              
    ##                Accuracy : 0.748              
    ##                  95% CI : (0.704, 0.789)     
    ##     No Information Rate : 0.753              
    ##     P-Value [Acc > NIR] : 0.614              
    ##                                              
    ##                   Kappa : 0.054              
    ##                                              
    ##  Mcnemar's Test P-Value : <0.0000000000000002
    ##                                              
    ##             Sensitivity : 0.9719             
    ##             Specificity : 0.0667             
    ##          Pos Pred Value : 0.7604             
    ##          Neg Pred Value : 0.4375             
    ##              Prevalence : 0.7529             
    ##          Detection Rate : 0.7318             
    ##    Detection Prevalence : 0.9624             
    ##       Balanced Accuracy : 0.5193             
    ##                                              
    ##        'Positive' Class : No                 
    ## 

``` r
confusionMatrix(table(MLogReg.pred.test,test$fracture))
```

    ## Confusion Matrix and Statistics
    ## 
    ##                  
    ## MLogReg.pred.test No Yes
    ##               No  54  16
    ##               Yes  1   4
    ##                                         
    ##                Accuracy : 0.773         
    ##                  95% CI : (0.662, 0.862)
    ##     No Information Rate : 0.733         
    ##     P-Value [Acc > NIR] : 0.260910      
    ##                                         
    ##                   Kappa : 0.239         
    ##                                         
    ##  Mcnemar's Test P-Value : 0.000685      
    ##                                         
    ##             Sensitivity : 0.982         
    ##             Specificity : 0.200         
    ##          Pos Pred Value : 0.771         
    ##          Neg Pred Value : 0.800         
    ##              Prevalence : 0.733         
    ##          Detection Rate : 0.720         
    ##    Detection Prevalence : 0.933         
    ##       Balanced Accuracy : 0.591         
    ##                                         
    ##        'Positive' Class : No            
    ## 

``` r
################ Misclassification rate train ######################################
cross.table.manual.train = table(MLogReg.pred.train,train$fracture)
MCR_manual.train = (cross.table.manual.train[2]+cross.table.manual.train[3])/dim(train)[1]
print(paste('Misclassification Rate for manual selection on training set: ', MCR_manual.train))
```

    ## [1] "Misclassification Rate for manual selection on training set:  0.251764705882353"

``` r
################ Misclassification rate test ######################################
cross.table.manual.test = table(MLogReg.pred.test,test$fracture)
MCR_manual.test = (cross.table.manual.test[2]+cross.table.manual.test[3])/dim(test)[1]
print(paste('Misclassification Rate for manual selection on test set: ', MCR_manual.test))
```

    ## [1] "Misclassification Rate for manual selection on test set:  0.226666666666667"

``` r
################ Manual ROC Curve ######################################
#library(ROCR)
#manual.roc.test = prediction(as.numeric(MLogReg.pred.test), as.numeric(test$fracture), label.ordering=c("No","Yes"))
#roc.manual.test = performance(manual.roc.test, measure = "tpr", x.measure = "fpr")
#plot(roc.manual.test,colorize = TRUE)
#abline(a=0, b= 1)
```

##################################################################################### 

# Objective II

##################################################################################### 

``` r
cutoff.manual = 0.5
LogReg.Interaction<-glm(fracture~age + height + bmi + priorfrac, family=binomial(link='logit'))
LogReg.Interaction.pred.train <- predict(LogReg.Interaction, newdata = train, type = "response")

Anova(LogReg.Interaction,type=3)
```

    ## Analysis of Deviance Table (Type III tests)
    ## 
    ## Response: fracture
    ##           LR Chisq Df Pr(>Chisq)    
    ## age          10.19  1     0.0014 ** 
    ## height        4.28  1     0.0386 *  
    ## bmi           1.30  1     0.2543    
    ## priorfrac    11.50  1     0.0007 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Confusion matrix
class.LogReg.Interaction.train<-factor(ifelse(LogReg.Interaction.pred.train>cutoff.manual,"Yes","No"),levels=c("No","Yes"))

#Confusion Matrix for interaction
print("Confusion matrix for Interaction TRAINING with 0.5 cutoff")
```

    ## [1] "Confusion matrix for Interaction TRAINING with 0.5 cutoff"

``` r
confusionMatrix(table(class.LogReg.Interaction.train,train$fracture))
```

    ## Confusion Matrix and Statistics
    ## 
    ##                               
    ## class.LogReg.Interaction.train  No Yes
    ##                            No  305  96
    ##                            Yes  15   9
    ##                                             
    ##                Accuracy : 0.739             
    ##                  95% CI : (0.694, 0.78)     
    ##     No Information Rate : 0.753             
    ##     P-Value [Acc > NIR] : 0.769             
    ##                                             
    ##                   Kappa : 0.052             
    ##                                             
    ##  Mcnemar's Test P-Value : 0.0000000000000312
    ##                                             
    ##             Sensitivity : 0.9531            
    ##             Specificity : 0.0857            
    ##          Pos Pred Value : 0.7606            
    ##          Neg Pred Value : 0.3750            
    ##              Prevalence : 0.7529            
    ##          Detection Rate : 0.7176            
    ##    Detection Prevalence : 0.9435            
    ##       Balanced Accuracy : 0.5194            
    ##                                             
    ##        'Positive' Class : No                
    ## 

``` r
################ Misclassification rate train ######################################
cross.table.interaction.train = table(class.LogReg.Interaction.train,train$fracture)
MCR_interaction.train = (cross.table.interaction.train[2]+cross.table.interaction.train[3])/dim(train)[1]
print(paste('Misclassification Rate for stepwise selection on training set: ', MCR_interaction.train))
```

    ## [1] "Misclassification Rate for stepwise selection on training set:  0.261176470588235"
