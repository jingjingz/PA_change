# supp material to manuscript 
# "A Riemann Manifold Model Framework for Longitudinal Changes in Physical Activity"
# Example R code for data analysis 


#--------------------------------------------------------
# load NEW momenta_all 
#--------------------------------------------------------
momenta_all <- matrix(0, 177, 2160)
for (id in 1: 177){
  filename <- paste0("/results/MENU_defor_est_", id, ".csv")
  momenta <- read.csv(file = filename, header = F)
  dim(momenta)
  momenta_all[id, ] <- c(momenta[, 1], momenta[, 2])
}

# save combined momenta
write.table(momenta_all, file = "/results/momenta_all.csv", row.names = F, col.names = F, sep = ",")


#--------------------------------------------------------
# FPCA on new momenta with package fdapace
#--------------------------------------------------------
library(fdapace)

s <- c(1: dim(momenta_all)[2])
L3 <- MakeFPCAInputs(IDs = rep(1:dim(momenta_all)[1], each=length(s)), tVec=rep(s,dim(momenta_all)[1]), t(momenta_all))
FPCAdense <- FPCA(L3$Ly, L3$Lt, optns = list(usergrid = T))
plot(FPCAdense)

FPCAdense$sigma2
FPCAdense$cumFVE 
plot(FPCAdense$cumFVE, type = 'l')

#-------------------------------------------------
# get top PCs
PC.est <- FPCAdense$phi # PCs
PC_select <- PC.est[, 1:30]

# calculate PC scores 
PCscores <- FPCAdense$xiEst # PC scores
PCscores <- PCscores[, 1:30]

#-----------------------------------------------
# standardize scores of each PC 
PCscores.sd <- apply(PCscores, 2, sd)
PC.rescale <- PC_select * (rep(1, dim(PC_select)[1]) %*% t(PCscores.sd))
scores.rescale <- PCscores / (rep(1, dim(PCscores)[1]) %*% t(PCscores.sd))

#-----------------------------------------------------------------------------------
# combine with covariates in data

# follow-up periods
A1 <- blood.l.sub[blood.l.sub$Timepoint == "A1", ]
A2 <- blood.l.sub[blood.l.sub$Timepoint == "A2", ]
A3 <- blood.l.sub[blood.l.sub$Timepoint == "A3", ]

# merge by ID
A13 <- merge(A1, A3, by.x = "PtID", by.y = "PtID", sort = F, suffixes = c(".1",".3"))
bmi_diff <- A13$bmi.3 - A13$bmi.1 # BMI (used in paper)
glu_diff <- A13$glu.3 - A13$glu.1 # Glucose
ins_diff <- A13$ins.3 - A13$ins.1 # Insulin

#-----------------------------------------------------------------------------------
# combine with PC scores
data <- data.frame(dat.base.sub, bmi_diff, glu_diff, ins_diff, scores.rescale) 

# only covariates (not including outcome)
data.x <- data.frame(dat.base.sub, scores.rescale) 


#-------------------------------------------------------------------
# summary table of variables
#-------------------------------------------------------------------
# check type of variables
library(gtsummary)
library(gt)

table.all <- data.x %>% select(HOMA, agerand, Diet, race, Edyears, Marital, EverSmoke, AnyCancer, HeartProbsQ19, PDFDEPRESSION, PDFHeartPROBLEM, PDFHYPERTENSION, pain, workmets, FSH, p2)

# format factors
levels(table.all$race) <- c("White", "Black", "Asian", "Pacific Islander", "Native American", "Mixed Race", "Other Race")
levels(table.all$Marital) <- c("Married or living together", "Single-never married", "Widowed", "Divorced", "Separated")
levels(table.all$HeartProbsQ19) <- c("None"," High blood pressure"," High blood pressure and High cholesterol"," High blood pressure and Other"," High cholesterol"," Other"
)
levels(table.all$EverSmoke) <- c("No", "Yes")


# summarize the data with our package
table1 <- 
  tbl_summary(table.all, 
              type = list(HOMA ~ "continuous", agerand ~ "continuous", Diet~"categorical", Edyears~"continuous", Marital~"categorical", EverSmoke~"dichotomous", AnyCancer~"dichotomous", HeartProbsQ19~"categorical", PDFDEPRESSION~"dichotomous",PDFHeartPROBLEM~"dichotomous", PDFHYPERTENSION~"dichotomous", pain~"continuous", workmets~ "continuous", FSH~"continuous"), 
              label = list(HOMA ~ "HOMA", agerand ~ "Age", Diet ~ "Diet", race ~ "Race", Edyears ~ "Education Years", Marital ~ "Marital status", EverSmoke ~ "Ever Smoked", AnyCancer ~ "Any Cancer", HeartProbsQ19 ~ "Heart Related Problems", PDFDEPRESSION ~ "Depression Meds", PDFHeartPROBLEM ~ "Heart Problem Meds", PDFHYPERTENSION ~ "Hypertension Meds", pain ~ "Pain Scale", workmets ~ "Weekly Met-Minutes", FSH ~ "Follicle stimulating hormone IU/L", p2 ~ "Vigorous activity at work days")) %>% 
  as_flex_table() %>% 
  flextable::save_as_docx(., path = "/Users/izoujj/Desktop/Projects/Accelerometer/Tuo/draft/tables/MENU_summ_new.docx")

table1




#---------------------------------------------------
# check distribution of outcome
#--------------------------------------------------
typeof(bmi_diff)
hist(bmi_diff, breaks = 30)

qqnorm(bmi_diff)
qqline(bmi_diff)


#--------------------------------------------------
# LASSO variable selection
#--------------------------------------------------
cv.lasso <- cv.glmnet(x, as.matrix(bmi_diff), nfolds = 10)

c = coef(cv.lasso, cv.lasso$lambda.min)

par(mfrow = c(1,1))
plot(cv.lasso)

vars.select <- cbind(rownames(c)[which(c!=0)], format(c[which(c!=0)], digits = 2))
vars.select

# example: re-run of regression with Lasso-selected vars only
y <- bmi_diff
data.x.lm <- data.x[, c("Diet","race","Marital","HeartProbsQ19","pain","FSH","X1","X10","X20","X25")]
tst <- model.matrix( ~ ., data.x.lm)
l.sig.2 <- lm(y ~ tst+0)
summary(l.sig.2)
plot(l.sig.2) # diagnostics













