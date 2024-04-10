This code was documented by Ryo Niwa.

```R=
# library loading
library(tidyverse)
library(data.table)

# Load the raw data
clinvar = fread('gunzip -c ~/Downloads/variant_summary.txt.gz')

# Identify the variant types in both GRCh37 and GRCh38
var_type_version = clinvar %>%
    group_by(Type, Assembly) %>%
    summarize(count = n(), .groups = 'drop') %>%
    spread(key = Assembly, value = count)

var_type_version
# A tibble: 15 × 5
   Type                       GRCh37  GRCh38    na NCBI36
   <chr>                       <int>   <int> <int>  <int>
 1 Complex                        64       3    12      3
 2 copy number gain            21337    6907     1   2898
 3 copy number loss            19916    6655    NA   1749
 4 Deletion                   134069  123853  1423     71
 5 Duplication                 62415   56048   315     50
 6 fusion                          1      NA     5     NA
 7 Indel                       13962   13968    54     NA
 8 Insertion                   11567   11561   282     NA
 9 Inversion                    1258    1249    61     NA
10 Microsatellite              32562   32574    52     NA
11 protein only                   NA      NA    94     NA
12 single nucleotide variant 2492762 2492724   605     NA
13 Tandem duplication              1      NA    NA     NA
14 Translocation                 105      22   210     NA
15 Variation                     354     359   142     NA

# Detecting differences between GRCh37 and GRCh38
GRCh37_alleleIDs <- clinvar %>%
    filter(Assembly == "GRCh37") %>%
    distinct(`#AlleleID`)

GRCh38_alleleIDs <- clinvar %>%
    filter(Assembly == "GRCh38") %>%
    distinct(`#AlleleID`)

# Find the intersection of AlleleIDs between the two versions
overlap_alleleIDs <- intersect(GRCh37_alleleIDs$`#AlleleID`, GRCh38_alleleIDs$`#AlleleID`)

length(overlap_alleleIDs)
274386

gap_GRCh37 = count(GRCh37_alleleIDs) - length(overlap_alleleIDs)
1404

gap_GRCh38 = count(GRCh38_alleleIDs) - length(overlap_alleleIDs)
45793

# Alternatively you can use VennDiagram package
library(VennDiagram)

# Prepare data for Venn Diagram
list_of_alleles <- list(
    GRCh37 = GRCh37_alleleIDs$`#AlleleID`,
    GRCh38 = GRCh38_alleleIDs$`#AlleleID`
)

# Generate Venn Diagram
venn.plot <- venn.diagram(
    x = list_of_alleles,
    category.names = c("GRCh37", "GRCh38"),
    filename = "AlleleID.png"
)

# We concluded GRCh38 mostly represented the current ClinVar entries
# Filter the data by GRCh38 assembly
clinvar_GRCh38 <- clinvar %>% 
filter(Assembly == "GRCh38")

var_type = clinvar_GRCh38 %>% 
group_by(Type) %>%
summarize(count = n())

var_type
# A tibble: 12 × 2
   Type                        count
   <chr>                       <int>
 1 Complex                         3
 2 Deletion                   123853
 3 Duplication                 56048
 4 Indel                       13968
 5 Insertion                   11561
 6 Inversion                    1249
 7 Microsatellite              32574
 8 Translocation                  22
 9 Variation                     359
10 copy number gain             6907
11 copy number loss             6655
12 single nucleotide variant 2492724

# Use the output for drawing the stacked barplots
# Next let's focus on clinical significance

clin_type = clinvar_GRCh38 %>% 
group_by(ClinicalSignificance) %>%
summarize(count = n())

clin_type
# A tibble: 92 × 2
   ClinicalSignificance                        count
   <chr>                                       <int>
 1 -                                               1
 2 Affects                                       139
 3 Affects; risk factor                            1
 4 Benign                                     206064
 5 Benign/Likely benign                        47413
 6 Benign/Likely benign; association               1
 7 Benign/Likely benign; drug response             2
 8 Benign/Likely benign; drug response; other      3
 9 Benign/Likely benign; other                     9
10 Benign/Likely benign; other; risk factor        1
# ℹ 82 more rows
# ℹ Use `print(n = ...)` to see more rows

# Filter the GRCh38 dataset by Pathogenic
var_type_path = clinvar_GRCh38 %>% 
    dplyr::filter(ClinicalSignificance == "Pathogenic") %>%
    group_by(Type) %>%
    summarize(count = n())
```