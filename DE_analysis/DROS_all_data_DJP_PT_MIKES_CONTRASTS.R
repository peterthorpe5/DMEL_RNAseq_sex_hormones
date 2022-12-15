

library(plotly)
library(ggplot2)
library(edgeR)


setwd("C:/Users/pjt6/OneDrive - University of St Andrews/2022Sep05_RNAseq_DROS_Mike_R_Tanya/RNAseq/gene")

cts <- read.delim("DMEL.genes.isoform.counts - Copy.matrix", 
                  header=T, row.names=1)

coldata = read.table("metadata_2 - Copy.txt", header=T, com='', 
                      sep="\t", check.names=F, row.names=1)


head(cts)


head(coldata)


coldata <- coldata[,c("condition", "treatment", "body_part", "batch", "reps", "strain")]
coldata$condition <- factor(coldata$condition)
coldata$strain <- factor(coldata$strain)
coldata$side <- factor(coldata$body_part)
coldata$Age_of_bat <- factor(coldata$batch)
coldata$reps <- factor(coldata$reps)
strain <- factor(coldata$strain)
condition <- factor(coldata$condition)
treatment <- factor(coldata$treatment)
batch <- factor(coldata$batch)
body_part <- factor(coldata$body_part)



all(rownames(coldata) %in% colnames(cts))
# TRUE

rownames(coldata) <-  rownames(coldata)

all(rownames(coldata) == colnames(cts))
# FALSE

cts <- cts[, rownames(coldata)]

all(rownames(coldata) == colnames(cts))
# TRUE
cts <- cts[, rownames(coldata)]
all(rownames(coldata) == colnames(cts))





# explicitly set up the group (yes, this is contained in the table!)

# these get reordered to this - they are not in this order in the matrix,
# but this is reorderd as above
group=c('CSHA','CSHA','CSHA','CSHA','CSHA',
        'CSHT','CSHT','CSHT','CSHT','CSHT',
        'CSCA','CSCA','CSCA','CSCA','CSCA',
        'CSCT','CSCT','CSCT','CSCT','CSCT',
        'CSFA','CSFA','CSFA','CSFA','CSFA',
        'CSFT','CSFT','CSFT','CSFT','CSFT',
        'ITHA','ITHA','ITHA','ITHA','ITHT',
        'ITHT','ITHT','ITHT','ITCA','ITCA',
        'ITCA','ITCA','ITCT','ITCT','ITCT',
        'ITCT','ITFA','ITFA','ITFA','ITFA',
        'ITFT','ITFT','ITFT','ITFT')





#0.3 R Packages now to start using edgeR to do some analysis: https://bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeR.pdf



y <- DGEList(counts=cts, group=group)

condition<-coldata$condition # just so I dont have to keep using dollars :)



#filter low expressing genes. Stats perfomred in a larger number, when lots are zero to low expressing negatively impacts on the data. 5,640 genes dropped from the analysis. 


keep <- filterByExpr(y)

table(keep)



# calculate nomrlaisation factors and normalise the data for lib depth differences



y <- y[keep, , keep.lib.sizes=FALSE]
#The TMM normalization is applied to account for the compositional biases:
y <- calcNormFactors(y)

y$samples




#following this https://support.bioconductor.org/p/56637/ with a problem with y, see the fix at the webpage. Reasign to d:


d <- DGEList(counts=y,group=group)

keep <- filterByExpr(d)

table(keep)

d <- d[keep, , keep.lib.sizes=FALSE]
d <- calcNormFactors(d)

# 
# 
# Before we fit GLMs, we need to define our design matrix based on the experimental design. 
# We want to test for differential expressions between our conditions within batches, i.e. adjusting for differences between batches. In statistical terms,
# this is an additive linear model. So the design matrix is created as:

#design_i <- model.matrix(~0 + treatment*strain + batch)
#d_i <- estimateDisp(d, design_i)
#fit_i <- glmQLFit(d_i, design_i)
##plotQLDisp(fit_i)
##d_i$common.dispersion
#sqrt(d_i$common.dispersion) ## this is the (common) BCV ## this value (1.162545) is very high! 

## we are removing batch in th emodel

design_i <- model.matrix(~0 + treatment*strain)
d_i <- estimateDisp(d, design_i)
fit_i <- glmQLFit(d_i, design_i)
plotQLDisp(fit_i)
d_i$common.dispersion
#sqrt(d_i$common.dispersion) ## this is the (common) BCV ## this value (1.162545) is very high! 
# now with batch removed - still very high.  1.122826
# P.s I notice that your biological coeff of var is massive (>1, normally I have values ~0.3 
#for my own noisy data (for nice exps it should be in the 0.1 area)) - are there multiple tissues or something else 
# with a huge effect here? You might want to consider splitting the analysis into different models if possible.

sqrt(d_i$common.dispersion)



### what coefs do I want?
colnames(design_i)


# [1] "treatmentblank_dummy"               
# [2] "treatmentdead_female"               
# [3] "treatmentHD_dummy"                  
# [4] "strainIT_IV_69"                     
# [5] "batch2"                             
# [6] "batch3"                             
# [7] "batch4"                             
# [8] "batch5"                             
# [9] "treatmentdead_female:strainIT_IV_69"
# [10] "treatmentHD_dummy:strainIT_IV_69"  

#lrt_inter <- glmQLFTest(fit_i,coef=9:10) ### these are the interaction ones
#topTags(lrt_inter) 

### this is the best way to get genes showing an interaction
### you can extract 'main effect terms' in a similar way:

#lrt_batch <- glmQLFTest(fit_i,coef=5:8) ### batch
#topTags(lrt_batch) 

# now we have removed batch

#[1] "treatmentblank_dummy"                "treatmentdead_female"                "treatmentHD_dummy"                   "strainIT_IV_69"                     
#[5] "treatmentdead_female:strainIT_IV_69" "treatmentHD_dummy:strainIT_IV_69"  

lrt_inter <- glmQLFTest(fit_i,coef=5:6) ### these are the interaction ones
topTags(lrt_inter) 




#######################
### make contrasts
### not sure what contrasts you want, but this should give you the idea

make_Group_matrix <- function(y, fact1,fact2){
  a_samp = data.frame(Sample=colnames(y),fact1,fact2)
  Group <- factor(paste(a_samp$fact1,a_samp$fact2, sep="."))
  print(Group)
  cbind(a_samp,Group=Group)
  cat(Group)
  G_design <- model.matrix(~0+Group)
  colnames(G_design) <- levels(Group)
  print(G_design)
  return(G_design)
}

mat_test <- make_Group_matrix(d,substr(colnames(d) , start = 1 , stop = 2),substr(colnames(d) , start = 3 , stop = 4))
mat_test

# level: Levels: CS.CA CS.CT CS.FA CS.FT CS.HA CS.HT IT.CA IT.CT IT.FA IT.FT IT.HA IT.HT

fit_glmQ_test <-  glmQLFit(d_i, mat_test) ### I think using the full model is OK here.




######################################################################################
######################################################################################

# REMOVE BATCH FROM MODEL ############


### but... this can be a problem. It is often best to extract genes showing a main effect from an additive model
### (doing it this way gives results as you would expect from a classical anova) - it is a v. small effect normally.
### e.g. :

design_a <- model.matrix(~0 + treatment + strain)
d_a <- estimateDisp(d, design_a)
fit_a <- glmQLFit(d_a, design_a)
plotQLDisp(fit_a)
colnames(design_a)
lrt_batch_a <- glmQLFTest(fit_a,coef=1:4) ### batch removed from additive model 
topTags(lrt_batch_a) 


#######################
### make contrasts
### not sure what contrasts you want, but this should give you the idea

make_Group_matrix <- function(y, fact1,fact2){
  a_samp = data.frame(Sample=colnames(y),fact1,fact2)
  Group <- factor(paste(a_samp$fact1,a_samp$fact2, sep="."))
  print(Group)
  cbind(a_samp,Group=Group)
  cat(Group)
  G_design <- model.matrix(~0+Group)
  colnames(G_design) <- levels(Group)
  print(G_design)
  return(G_design)
}

mat_test <- make_Group_matrix(d,substr(colnames(d) , start = 1 , stop = 2),substr(colnames(d) , start = 3 , stop = 4))
mat_test

fit_glmQ_test <-  glmQLFit(d_i, mat_test) ### I think using the full model is OK here.


##### Question 1.	Do  h/t respond differently to treatments?


# 1. For T tissue only (1.	Do  h/t respond differently to treatments?)
test.contrasts <- makeContrasts(CSTC = CS.CT - ((CS.FT + CS.HT)/2), levels=mat_test) 
CSTC <- glmQLFTest(fit_glmQ_test, contrast=test.contrasts[,"CSTC"])

test.contrasts <- makeContrasts(CSTF = CS.FT - ((CS.CT + CS.HT)/2), levels=mat_test) 
CSTF <- glmQLFTest(fit_glmQ_test, contrast=test.contrasts[,"CSTF"])

test.contrasts <- makeContrasts(CSTH = CS.HT - ((CS.CT + CS.FT)/2), levels=mat_test) 
CSTH <- glmQLFTest(fit_glmQ_test, contrast=test.contrasts[,"CSTH"])

test.contrasts <- makeContrasts(ITTC = IT.CT - ((IT.FT + IT.HT)/2), levels=mat_test) 
ITTC <- glmQLFTest(fit_glmQ_test, contrast=test.contrasts[,"ITTC"])

test.contrasts <- makeContrasts(ITTF = IT.FT - ((IT.CT + IT.HT)/2), levels=mat_test) 
ITTF <- glmQLFTest(fit_glmQ_test, contrast=test.contrasts[,"ITTF"])

test.contrasts <- makeContrasts(ITTH = IT.HT - ((IT.CT + IT.FT)/2), levels=mat_test) 
ITTH <- glmQLFTest(fit_glmQ_test, contrast=test.contrasts[,"ITTH"])


####
tTags = topTags(CSTC,n=NULL)

result_table = tTags$table

result_table = data.frame(sampleA="CSTC", sampleB="average_of_others", result_table)
result_table = merge(result_table, cts,by="row.names",all.x=TRUE)
result_table <- result_table[order(result_table$logFC),]


write.table(result_table, file='./question_1_Mike_contrasts/CSTC.ANOVA_LIKE.edgeR.DE_results.txt', sep='	', quote=F, row.names=T, col.names = NA)
###

tTags = topTags(CSTF,n=NULL)


result_table = tTags$table

result_table = data.frame(sampleA="CSTF", sampleB="average_of_others", result_table)
result_table = merge(result_table, cts,by="row.names",all.x=TRUE)
result_table <- result_table[order(result_table$logFC),]


write.table(result_table, file='./question_1_Mike_contrasts/CSTF.ANOVA_LIKE.edgeR.DE_results.txt', sep='	', quote=F, row.names=T, col.names = NA)

###


tTags = topTags(CSTH,n=NULL)

result_table = tTags$table

result_table = data.frame(sampleA="CSTH", sampleB="average_of_others", result_table)
result_table = merge(result_table, cts,by="row.names",all.x=TRUE)
result_table <- result_table[order(result_table$logFC),]


write.table(result_table, file='./question_1_Mike_contrasts/CSTH.ANOVA_LIKE.edgeR.DE_results.txt', sep='	', quote=F, row.names=T, col.names = NA)

####

tTags = topTags(ITTC,n=NULL)

result_table = tTags$table

result_table = data.frame(sampleA="ITTC", sampleB="average_of_others", result_table)
result_table = merge(result_table, cts,by="row.names",all.x=TRUE)
result_table <- result_table[order(result_table$logFC),]


write.table(result_table, file='./question_1_Mike_contrasts/ITTC.ANOVA_LIKE.edgeR.DE_results.txt', sep='	', quote=F, row.names=T, col.names = NA)




#####
tTags = topTags(ITTF,n=NULL)

result_table = tTags$table

result_table = data.frame(sampleA="ITTF", sampleB="average_of_others", result_table)
result_table = merge(result_table, cts,by="row.names",all.x=TRUE)
result_table <- result_table[order(result_table$logFC),]


write.table(result_table, file='./question_1_Mike_contrasts/ITTF.ANOVA_LIKE.edgeR.DE_results.txt', sep='	', quote=F, row.names=T, col.names = NA)

###


tTags = topTags(ITTH,n=NULL)

result_table = tTags$table

result_table = data.frame(sampleA="ITTH", sampleB="average_of_others", result_table)
result_table = merge(result_table, cts,by="row.names",all.x=TRUE)
result_table <- result_table[order(result_table$logFC),]


write.table(result_table, file='./question_1_Mike_contrasts/ITTH.ANOVA_LIKE.edgeR.DE_results.txt', sep='	', quote=F, row.names=T, col.names = NA)


# In words, that will give us a list of DE genes in CantonSThorax in the Control, Female and Hormone treatments, then the same for Italy.

################################################################

#### QUESTION 2 : For A tissue only

test.contrasts <- makeContrasts(CSAC = CS.CA - ((CS.FA + CS.HA)/2) , levels=mat_test) 
CSAC <- glmQLFTest(fit_glmQ_test, contrast=test.contrasts[,"CSAC"])

test.contrasts <- makeContrasts(CSAF = CS.FA - ((CS.CA + CS.HA)/2) , levels=mat_test) 
CSAF <- glmQLFTest(fit_glmQ_test, contrast=test.contrasts[,"CSAF"])

test.contrasts <- makeContrasts(CSAH = CS.HA - ((CS.CA + CS.FA)/2) , levels=mat_test) 
CSAH <- glmQLFTest(fit_glmQ_test, contrast=test.contrasts[,"CSAH"])

test.contrasts <- makeContrasts(ITAC = IT.CA - ((IT.FA + IT.HA)/2) , levels=mat_test) 
ITAC <- glmQLFTest(fit_glmQ_test, contrast=test.contrasts[,"ITAC"])

test.contrasts <- makeContrasts(ITAF = IT.FA - ((IT.CA + IT.HA)/2) , levels=mat_test) 
ITAF <- glmQLFTest(fit_glmQ_test, contrast=test.contrasts[,"ITAF"])

test.contrasts <- makeContrasts(ITAH = IT.HA - ((IT.CA + IT.FA)/2) , levels=mat_test) 
ITAH <- glmQLFTest(fit_glmQ_test, contrast=test.contrasts[,"ITAH"])


###
tTags = topTags(CSAC,n=NULL)

result_table = tTags$table

result_table = data.frame(sampleA='CSAC', sampleB='average_of_others', result_table)
result_table = merge(result_table, cts,by='row.names',all.x=TRUE)
result_table <- result_table[order(result_table$logFC),]


write.table(result_table, file='./question_2_Mike_contrasts/CSAC.ANOVA_LIKE.edgeR.DE_results.txt', sep='	', quote=F, row.names=T, col.names = NA)
######


tTags = topTags(CSAF,n=NULL)

result_table = tTags$table

result_table = data.frame(sampleA='CSAF', sampleB='average_of_others', result_table)
result_table = merge(result_table, cts,by='row.names',all.x=TRUE)
result_table <- result_table[order(result_table$logFC),]


write.table(result_table, file='./question_2_Mike_contrasts/CSAF.ANOVA_LIKE.edgeR.DE_results.txt', sep='	', quote=F, row.names=T, col.names = NA)
######


tTags = topTags(CSAH,n=NULL)

result_table = tTags$table

result_table = data.frame(sampleA='CSAH', sampleB='average_of_others', result_table)
result_table = merge(result_table, cts,by='row.names',all.x=TRUE)
result_table <- result_table[order(result_table$logFC),]


write.table(result_table, file='./question_2_Mike_contrasts/CSAH.ANOVA_LIKE.edgeR.DE_results.txt', sep='	', quote=F, row.names=T, col.names = NA)
######


tTags = topTags(ITAC,n=NULL)

result_table = tTags$table

result_table = data.frame(sampleA='ITAC', sampleB='average_of_others', result_table)
result_table = merge(result_table, cts,by='row.names',all.x=TRUE)
result_table <- result_table[order(result_table$logFC),]


write.table(result_table, file='./question_2_Mike_contrasts/ITAC.ANOVA_LIKE.edgeR.DE_results.txt', sep='	', quote=F, row.names=T, col.names = NA)
######


tTags = topTags(ITAF,n=NULL)

result_table = tTags$table

result_table = data.frame(sampleA='ITAF', sampleB='average_of_others', result_table)
result_table = merge(result_table, cts,by='row.names',all.x=TRUE)
result_table <- result_table[order(result_table$logFC),]


write.table(result_table, file='./question_2_Mike_contrasts/ITAF.ANOVA_LIKE.edgeR.DE_results.txt', sep='	', quote=F, row.names=T, col.names = NA)
######


tTags = topTags(ITAH,n=NULL)

result_table = tTags$table

result_table = data.frame(sampleA='ITAH', sampleB='average_of_others', result_table)
result_table = merge(result_table, cts,by='row.names',all.x=TRUE)
result_table <- result_table[order(result_table$logFC),]


write.table(result_table, file='./question_2_Mike_contrasts/ITAH.ANOVA_LIKE.edgeR.DE_results.txt', sep='	', quote=F, row.names=T, col.names = NA)
######


####################################################################################
#########
# 3 ) For abdomens

test.contrasts <- makeContrasts(CS.CA_vs_IT.CA = CS.CA - IT.CA , levels=mat_test)
CS.CA_vs_IT.CA <- glmQLFTest(fit_glmQ_test, contrast=test.contrasts[,"CS.CA_vs_IT.CA"])

test.contrasts <- makeContrasts(CS.FA_vs_IT.FA = CS.FA - IT.FA , levels=mat_test)
CS.FA_vs_IT.FA <- glmQLFTest(fit_glmQ_test, contrast=test.contrasts[,"CS.FA_vs_IT.FA"])

test.contrasts <- makeContrasts(CS.HA_vs_IT.HA = CS.HA - IT.HA , levels=mat_test)
CS.HA_vs_IT.HA <- glmQLFTest(fit_glmQ_test, contrast=test.contrasts[,"CS.HA_vs_IT.HA"])


###
tTags = topTags(CS.CA_vs_IT.CA,n=NULL)

result_table = tTags$table

result_table = data.frame(sampleA='CS.CA', sampleB='IT.CA', result_table)
result_table = merge(result_table, cts,by='row.names',all.x=TRUE)
result_table <- result_table[order(result_table$logFC),]


write.table(result_table, file='./question_3_Mike_contrasts/CS.CA_vs_IT.CA.ANOVA_LIKE.abdomens.edgeR.DE_results.txt', sep='	', quote=F, row.names=T, col.names = NA)
######

###
tTags = topTags(CS.FA_vs_IT.FA,n=NULL)

result_table = tTags$table

result_table = data.frame(sampleA='CS.FA', sampleB='IT.FA', result_table)
result_table = merge(result_table, cts,by='row.names',all.x=TRUE)
result_table <- result_table[order(result_table$logFC),]


write.table(result_table, file='./question_3_Mike_contrasts/CS.FA_vs_IT.FA.ANOVA_LIKE.abdomens.edgeR.DE_results.txt', sep='	', quote=F, row.names=T, col.names = NA)
######

###
tTags = topTags(CS.HA_vs_IT.HA,n=NULL)

result_table = tTags$table

result_table = data.frame(sampleA='CS.HA', sampleB='IT.HA', result_table)
result_table = merge(result_table, cts,by='row.names',all.x=TRUE)
result_table <- result_table[order(result_table$logFC),]


write.table(result_table, file='./question_3_Mike_contrasts/CS.HA_vs_IT.HA.ANOVA_LIKE.abdomens.edgeR.DE_results.txt', sep='	', quote=F, row.names=T, col.names = NA)
######



#########
# 3 ) For Thorax

test.contrasts <- makeContrasts(CS.CT_vs_IT.CT = CS.CT - IT.CT, levels=mat_test)
CS.CT_vs_IT.CT <- glmQLFTest(fit_glmQ_test, contrast=test.contrasts[,"CS.CT_vs_IT.CT"])

test.contrasts <- makeContrasts(CS.FT_vs_IT.FT = CS.FT - IT.FT, levels=mat_test)
CS.FT_vs_IT.FT <- glmQLFTest(fit_glmQ_test, contrast=test.contrasts[,"CS.FT_vs_IT.FT"])

test.contrasts <- makeContrasts(CS.HT_vs_IT.HT = CS.HT - IT.HT, levels=mat_test)
CS.HT_vs_IT.HT <- glmQLFTest(fit_glmQ_test, contrast=test.contrasts[,"CS.HT_vs_IT.HT"])


###
tTags = topTags(CS.CT_vs_IT.CT,n=NULL)

result_table = tTags$table

result_table = data.frame(sampleA='CS.CT', sampleB='IT.CT', result_table)
result_table = merge(result_table, cts,by='row.names',all.x=TRUE)
result_table <- result_table[order(result_table$logFC),]


write.table(result_table, file='./question_3_Mike_contrasts/CS.CT_vs_IT.CT.ANOVA_LIKE.THorax.edgeR.DE_results.txt', sep='	', quote=F, row.names=T, col.names = NA)
######

###
tTags = topTags(CS.FT_vs_IT.FT,n=NULL)

result_table = tTags$table

result_table = data.frame(sampleA='CS.FT', sampleB='IT.FT', result_table)
result_table = merge(result_table, cts,by='row.names',all.x=TRUE)
result_table <- result_table[order(result_table$logFC),]


write.table(result_table, file='./question_3_Mike_contrasts/CS.FT_vs_IT.FT.ANOVA_LIKE.HTorax.edgeR.DE_results.txt', sep='	', quote=F, row.names=T, col.names = NA)
######

###
tTags = topTags(CS.HT_vs_IT.HT,n=NULL)

result_table = tTags$table

result_table = data.frame(sampleA='CS.HT', sampleB='IT.HT', result_table)
result_table = merge(result_table, cts,by='row.names',all.x=TRUE)
result_table <- result_table[order(result_table$logFC),]


write.table(result_table, file='./question_3_Mike_contrasts/CS.HT_vs_IT.HT.ANOVA_LIKE.Thorax.edgeR.DE_results.txt', sep='	', quote=F, row.names=T, col.names = NA)
######

# q 3 more complex


#((CS.CA - (CS.FA + CS.HA)/2) - ((IT.CA - (IT.FA + IT.HA)/2)


test.contrasts <- makeContrasts(CS_effect_treat_differ = (CS.CA - ((CS.FA + CS.HA)/2)) - ((IT.CA - (IT.FA + IT.HA)/2)), levels=mat_test)
CS_effect_treat_differ <- glmQLFTest(fit_glmQ_test, contrast=test.contrasts[,"CS_effect_treat_differ"])


###
tTags = topTags(CS_effect_treat_differ,n=NULL)

result_table = tTags$table

result_table = data.frame(sampleA='CS_effect', sampleB='IT_effect', result_table)
result_table = merge(result_table, cts,by='row.names',all.x=TRUE)
result_table <- result_table[order(result_table$logFC),]


write.table(result_table, file='./question_3_Mike_contrasts/CS_effect_treat_differ_vs_IT.ANOVA_LIKE.THorax.edgeR.DE_results.txt', sep='	', quote=F, row.names=T, col.names = NA)
######


#######################################################
# EXTRAS:





test.contrasts <- makeContrasts(CS.CA = (CS.CA-(CS.FA+CS.HA)/2) - ((IT.CA-(IT.FA+IT.HA)/2)),  levels=mat_test)

CS.CA <- glmQLFTest(fit_glmQ_test, contrast=test.contrasts[,"CS.CA"])
topTags(CS.CA) 


tTags = topTags(CS.CA,n=NULL)

result_table = tTags$table

result_table = data.frame(sampleA="CS.CA", sampleB="equation", result_table)
result_table = merge(result_table, cts,by="row.names",all.x=TRUE)
result_table <- result_table[order(result_table$logFC),]


write.table(result_table, file='./extras/CS.CA.ANOVA_LIKE.edgeR.DE_results.txt', sep='	', quote=F, row.names=T, col.names = NA)

###############################

test.contrasts <- makeContrasts(CS.FA = (CS.FA-(CS.HA+CS.CA)/2) - ((IT.FA-(IT.HA+IT.CA)/2)),  levels=mat_test)

CS.FA <- glmQLFTest(fit_glmQ_test, contrast=test.contrasts[,"CS.FA"])
topTags(CS.FA) 


tTags = topTags(CS.FA,n=NULL)

result_table = tTags$table

result_table = data.frame(sampleA="CS.FA", sampleB="equation", result_table)
result_table = merge(result_table, cts,by="row.names",all.x=TRUE)
result_table <- result_table[order(result_table$logFC),]


write.table(result_table, file='./extras/CS.FA.ANOVA_LIKE.edgeR.DE_results.txt', sep='	', quote=F, row.names=T, col.names = NA)


############################################
#

test.contrasts <- makeContrasts(CS.HA = (CS.HA-(CS.FA+CS.CA)/2) - ((IT.HA-(IT.FA+IT.CA)/2)),  levels=mat_test)

CS.HA <- glmQLFTest(fit_glmQ_test, contrast=test.contrasts[,"CS.HA"])
topTags(CS.HA) 


tTags = topTags(CS.HA,n=NULL)

result_table = tTags$table

result_table = data.frame(sampleA="CS.HA", sampleB="equation", result_table)
result_table = merge(result_table, cts,by="row.names",all.x=TRUE)
result_table <- result_table[order(result_table$logFC),]


write.table(result_table, file='./extras/CS.HA.ANOVA_LIKE.edgeR.DE_results.txt', sep='	', quote=F, row.names=T, col.names = NA)

##################################
#

test.contrasts <- makeContrasts(CS.CT = (CS.CT-(CS.FT+CS.HT)/2) - ((IT.CT-(IT.FT+IT.HT)/2)),  levels=mat_test)

CS.CT <- glmQLFTest(fit_glmQ_test, contrast=test.contrasts[,"CS.CT"])
topTags(CS.CT) 


tTags = topTags(CS.CT,n=NULL)

result_table = tTags$table

result_table = data.frame(sampleA="CS.CT", sampleB="equation", result_table)
result_table = merge(result_table, cts,by="row.names",all.x=TRUE)
result_table <- result_table[order(result_table$logFC),]


write.table(result_table, file='./extras/CS.CT.ANOVA_LIKE.edgeR.DE_results.txt', sep='	', quote=F, row.names=T, col.names = NA)

#########################################
#

test.contrasts <- makeContrasts(CS.FT = (CS.FT-(CS.HT+CS.CT)/2) - ((IT.FT-(IT.HT+IT.CT)/2)),  levels=mat_test)

CS.FT <- glmQLFTest(fit_glmQ_test, contrast=test.contrasts[,"CS.FT"])
topTags(CS.FT) 


tTags = topTags(CS.FT,n=NULL)

result_table = tTags$table

result_table = data.frame(sampleA="CS.FT", sampleB="equation", result_table)
result_table = merge(result_table, cts,by="row.names",all.x=TRUE)
result_table <- result_table[order(result_table$logFC),]


write.table(result_table, file='./extras/CS.FT.ANOVA_LIKE.edgeR.DE_results.txt', sep='	', quote=F, row.names=T, col.names = NA)


#############################
#

test.contrasts <- makeContrasts(CS.HT = (CS.HT-(CS.FT+CS.CT)/2) - ((IT.HT-(IT.FT+IT.CT)/2)),  levels=mat_test)


CS.HT <- glmQLFTest(fit_glmQ_test, contrast=test.contrasts[,"CS.HT"])
topTags(CS.HT) 


tTags = topTags(CS.HT,n=NULL)

result_table = tTags$table

result_table = data.frame(sampleA="CS.HT", sampleB="equation", result_table)
result_table = merge(result_table, cts,by="row.names",all.x=TRUE)
result_table <- result_table[order(result_table$logFC),]


write.table(result_table, file='./extras/CS.HT.ANOVA_LIKE.edgeR.DE_results.txt', sep='	', quote=F, row.names=T, col.names = NA)

########################


test.contrasts <- makeContrasts(CS.CT2 = (CS.CT+IT.CT)-((CS.FT+CS.HT+IT.FT+IT.HT)/2),  levels=mat_test)


CS.CT2 <- glmQLFTest(fit_glmQ_test, contrast=test.contrasts[,"CS.CT2"])
topTags(CS.CT2) 


tTags = topTags(CS.CT2,n=NULL)

result_table = tTags$table

result_table = data.frame(sampleA="CS.CT2", sampleB="equation", result_table)
result_table = merge(result_table, cts,by="row.names",all.x=TRUE)
result_table <- result_table[order(result_table$logFC),]


write.table(result_table, file='./extras/CS.CT2.ANOVA_LIKE.edgeR.DE_results.txt', sep='	', quote=F, row.names=T, col.names = NA)

##########################
#
test.contrasts <- makeContrasts(CS.FT2 = (CS.FT+IT.FT)-((CS.CT+CS.HT+IT.CT+IT.HT)/2),  levels=mat_test)


CS.FT2 <- glmQLFTest(fit_glmQ_test, contrast=test.contrasts[,"CS.FT2"])
topTags(CS.FT2) 


tTags = topTags(CS.FT2,n=NULL)

result_table = tTags$table

result_table = data.frame(sampleA="CS.FT2", sampleB="equation", result_table)
result_table = merge(result_table, cts,by="row.names",all.x=TRUE)
result_table <- result_table[order(result_table$logFC),]


write.table(result_table, file='./extras/CS.FT2.ANOVA_LIKE.edgeR.DE_results.txt', sep='	', quote=F, row.names=T, col.names = NA)


#################
#
test.contrasts <- makeContrasts(CS.HT2 =   (CS.HT+IT.HT)-((CS.CT+CS.FT+IT.CT+IT.FT)/2),  levels=mat_test)


CS.HT2 <- glmQLFTest(fit_glmQ_test, contrast=test.contrasts[,"CS.HT2"])
topTags(CS.HT2) 


tTags = topTags(CS.HT2,n=NULL)

result_table = tTags$table

result_table = data.frame(sampleA="CS.HT2", sampleB="equation", result_table)
result_table = merge(result_table, cts,by="row.names",all.x=TRUE)
result_table <- result_table[order(result_table$logFC),]


write.table(result_table, file='./extras/CS.HT2.ANOVA_LIKE.edgeR.DE_results.txt', sep='	', quote=F, row.names=T, col.names = NA)


########################
#

test.contrasts <- makeContrasts(CS.CA2 = (CS.CA+IT.CA)-((CS.FA+CS.HA+IT.FA+IT.HA)/2),  levels=mat_test)


CS.CA2 <- glmQLFTest(fit_glmQ_test, contrast=test.contrasts[,"CS.CA2"])
topTags(CS.CA2) 


tTags = topTags(CS.CA2,n=NULL)

result_table = tTags$table

result_table = data.frame(sampleA="CS.CA2", sampleB="equation", result_table)
result_table = merge(result_table, cts,by="row.names",all.x=TRUE)
result_table <- result_table[order(result_table$logFC),]


write.table(result_table, file='./extras/CS.CA2.ANOVA_LIKE.edgeR.DE_results.txt', sep='	', quote=F, row.names=T, col.names = NA)

########################################
#
test.contrasts <- makeContrasts(CS.FA2 = (CS.FA+IT.FA)-((CS.CA+CS.HA+IT.CA+IT.HA)/2),  levels=mat_test)


CS.FA2 <- glmQLFTest(fit_glmQ_test, contrast=test.contrasts[,"CS.FA2"])
topTags(CS.FA2) 


tTags = topTags(CS.FA2,n=NULL)

result_table = tTags$table

result_table = data.frame(sampleA="CS.FA2", sampleB="equation", result_table)
result_table = merge(result_table, cts,by="row.names",all.x=TRUE)
result_table <- result_table[order(result_table$logFC),]


write.table(result_table, file='./extras/CS.FA2.ANOVA_LIKE.edgeR.DE_results.txt', sep='	', quote=F, row.names=T, col.names = NA)

#####################################
#
test.contrasts <- makeContrasts(CS.HA2 =   (CS.HA2 =   (CS.HA+IT.HA)-((CS.CA+CS.FA+IT.CA+IT.FA)/2)),  levels=mat_test)

CS.HA2 <- glmQLFTest(fit_glmQ_test, contrast=test.contrasts[,"CS.HA2"])
topTags(CS.HA2) 


tTags = topTags(CS.HA2,n=NULL)

result_table = tTags$table

result_table = data.frame(sampleA="CS.HA2", sampleB="equation", result_table)
result_table = merge(result_table, cts,by="row.names",all.x=TRUE)
result_table <- result_table[order(result_table$logFC),]


write.table(result_table, file='./extras/CS.HA2.ANOVA_LIKE.edgeR.DE_results.txt', sep='	', quote=F, row.names=T, col.names = NA)


