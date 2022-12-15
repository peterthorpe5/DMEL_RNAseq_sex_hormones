

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

design_i <- model.matrix(~0 + treatment*strain + batch)
d_i <- estimateDisp(d, design_i)
fit_i <- glmQLFit(d_i, design_i)
plotQLDisp(fit_i)
d_i$common.dispersion
sqrt(d_i$common.dispersion) ## this is the (common) BCV ## this value (1.162545) is very high! 



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

lrt_inter <- glmQLFTest(fit_i,coef=9:10) ### these are the interaction ones
topTags(lrt_inter) 

### this is the best way to get genes showing an interaction
### you can extract 'main effect terms' in a similar way:

lrt_batch <- glmQLFTest(fit_i,coef=5:8) ### batch
topTags(lrt_batch) 

### but... this can be a problem. It is often best to extract genes showing a main effect from an additive model
### (doing it this way gives results as you would expect from a classical anova) - it is a v. small effect normally.
### e.g. :

design_a <- model.matrix(~0 + treatment + strain + batch)
d_a <- estimateDisp(d, design_a)
fit_a <- glmQLFit(d_a, design_a)
plotQLDisp(fit_a)
colnames(design_a)
lrt_batch_a <- glmQLFTest(fit_a,coef=5:8) ### batch from additive model 
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

test.contrasts <- makeContrasts(
  CSHAvsCSHT	=	CS.HA - CS.HT,
  levels=mat_test)

CSHAvsCSHT <- glmQLFTest(fit_glmQ_test, contrast=test.contrasts[,"CSHAvsCSHT"])
topTags(CSHAvsCSHT) 

