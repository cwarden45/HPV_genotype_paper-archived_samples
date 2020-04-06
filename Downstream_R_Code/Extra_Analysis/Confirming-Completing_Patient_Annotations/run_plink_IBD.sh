#!/bin/sh
Combined_Prefix="combined"
PRUNE_PREFIX="pruned"

plink --lfile $Combined_Prefix --make-bed --out $Combined_Prefix --noweb
plink --bfile $Combined_Prefix --recode --tab --out $Combined_Prefix --noweb
plink --bfile $Combined_Prefix --noweb
plink --bfile $Combined_Prefix  --indep-pairwise 50 5 0.2 --out $Combined_Prefix --noweb
plink --bfile $Combined_Prefix --extract $Combined_Prefix.prune.in --recode --out $PRUNE_PREFIX
plink --file $PRUNE_PREFIX --genome --out $PRUNE_PREFIX --noweb