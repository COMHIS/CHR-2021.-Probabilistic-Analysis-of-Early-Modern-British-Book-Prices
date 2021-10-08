#WHAT: This script produces the results of the CHR 2021 Short Paper Probabilistic Analysis of Early Modern British Book Prices.
#Authors: Iiro Tiihonen, Mikko Tolonen, Leo Lahti.
#Unfortunately the script requires data that can not be shared, but it gives an overview of the process.
#For further information, contact iiro.tiihonen@helsinki.fi

#Call the relevant libraries
setwd("~/Git_root")
library(dplyr)
library(magrittr)
library(ggplot2)
library(XLConnect)
library(rstan)
library(gghsci)
library(gridExtra)

#Define functions used in the script.
#Name variables in a data frame.
name_data_frame <- function(X,nameVec){
  
  colnames(X) <- nameVec  
  
  return(X)  
}

#Convert years to a decade level accuracy. E.g 1634 -> 1630.
decade_of_year <- function(x){
  
  res <-  trunc(x/10)*10  
  return(res)
}

#The probabilistic robust regression model used in the paper. The parameter nu was fixed as it is, for example, in the STAN user's
#guide for implementing robust regression (https://mc-stan.org/docs/2_18/stan-users-guide/robust-noise-models.html). The value
#is slighly over 2 so as to allow very high but still finite variance with all parametrisations.
#Note that the model also generates the predictions.
stan_linear_model <-
"data {
  int<lower=0> N;
  int<lower=0> N_new;
  vector[N] x1;
  vector[N] x2;
  vector[N_new] x1_new;
  vector[N_new] x2_new;
  vector[N] y;
}
parameters {
  real beta0;
  real beta1;
  real beta2;
  real<lower=0> sigma;
  
}
model {
  y ~  student_t(2.02,beta0 + beta1 * x1 + beta2 * x2, sigma);
  beta0 ~ student_t(1,0,1);
  beta1 ~ student_t(1,0,1);
  beta2 ~ student_t(1,0,1);
  sigma ~ inv_gamma(1,1);
}
generated quantities {
  vector[N_new] y_new;
  for (n in 1:N_new)
    y_new[n] = student_t_rng(2.02,beta0 + x1_new[n] * beta1 + x2_new[n] * beta2, sigma);
}"





#Read relevant columns of price index information. Only price_index_pence is used in this analysis.
price_index <- readWorksheetFromFile(".../data/london.xlsx",sheet=1)[,c(1,117,119,121)]
colnames(price_index) <- c("year","price_index_silver","price_index_silver_extended","price_index_pence")

#Select the relevant rows of the price index and convert from factors to numbers.
price_index_lim <- price_index[13:nrow(price_index),]
price_index_lim$year <- as.integer(price_index_lim$year)
price_index_lim$price_index_silver <- as.numeric(price_index_lim$price_index_silver)
price_index_lim$price_index_silver_extended <- as.numeric(price_index_lim$price_index_silver_extended)
price_index_lim$price_index_pence <- as.numeric(price_index_lim$price_index_pence)

#The word price is mentioned in other contexts in the 500a than in the declaration of the print prodcuts price. Thereofore we
#made a table that only contains ids of print products for which the field 500a has a generic price declaration phrase
#we will filter the prices with this table
estc_price_reliable_ids <- read.csv(".../data_output/robust_prices.csv",stringsAsFactors = FALSE)
estc_price_reliable_ids$info <- NULL
#We will filter prices related to advertisements,prospectuses and pilot guides out, as these were more often about the price of
#something else than of the record itself (e.g about the price of the product being advertised)
estc_prices_unreliable_ids <- read.csv(".../data_output/unreliable_prices.csv",stringsAsFactors = FALSE)

#We will also filter out public records and periodicals, as it was often unclear whether the price was related
#to a single issue or part of the record or to the entirety of the record (both public records and periodicals are often coerced e.g. to annual collections)
estc_public_records <- read.csv(".../data_output/public_print.csv",stringsAsFactors = FALSE)

estc_periodicals <- read.csv(".../periodicals.csv",stringsAsFactors = FALSE)
colnames(estc_periodicals) <- "system_control_number"


#Read in the information about prices, physical properties and publication place of records from the local github repositories.
#Correct some erroneous prices with handchecked prices (this sometimes led to marking NA as in some cases the price was not that of the record itself or there was a serious error in e.g the plate estimate of the print product)
#Standardise variable names and formating between three handchecked files.

estc_prices_handchecked_all <-  read.csv(".../data/prices_handchecked_unified.csv")
estc_prices_handchecked_all$comment <- NULL

#Join price data with price index data. And select only values without obvious parsing errors (prices above 0)
estc_prices <- read.csv(".../estc_prices_extended.csv") %>% .[,c(1:8)] %>% left_join(.,price_index_lim) %>% subset(.,!is.na(total_price)) %>% subset(.,total_price>0)
colnames(estc_prices)[1] <- "system_control_number"

estc_prices_ids <- estc_prices$system_control_number %>% as.data.frame(.) %>% name_data_frame(.,"system_control_number")

estc_prices_handchecked_all <- inner_join(estc_prices_ids,estc_prices_handchecked_all)
#Match the corrected observations with the main table and replace the automatically parsed values with manually checked values.
estc_prices$total_price[match(estc_prices_handchecked_all$system_control_number,estc_prices$system_control_number)] <- estc_prices_handchecked_all$total_price


#Standardise the prices in the manner described in the article. Normalise the price in decimals of pence (original price in pences) with the price index
#This conversion results in prices that could also be used and were used in experimentation with e.g poisson models (they are counts > 1), but the final
#model presented here is a continuous regression model.
estc_prices$total_price <- estc_prices$total_price/estc_prices$price_index_pence
estc_prices$total_price <- estc_prices$total_price*10 
estc_prices$total_price <- round(estc_prices$total_price)
estc_prices$total_price <- as.integer(estc_prices$total_price)



#Read in the relevant bibliographic data sets about physical characteristics of observations.
estc_physical_extent <- read.csv(".../physicalextent.csv",sep="\t")

estc_physical_dimension <- read.csv(".../physical_dimension.csv",sep="\t")

#set the sheet size used to normalise paper to sheets.

sheet_size <- 5760

#Read in the place information about our observations.
estc_publication_place <- read.csv(".../Estc752Mappings.csv")
colnames(estc_publication_place)[1] <- "system_control_number"
estc_publication_place$system_control_number <- paste0("(CU-RivES)",estc_publication_place$system_control_number)
estc_publication_place <- estc_publication_place[,c("system_control_number","publication_place_752")]
colnames(estc_publication_place)[2] <- "publication_place"

#Select observations without critical missing data and limit to relevant variables. Filter the mentioned unreliable records out.
estc_prices_filtered <- left_join(estc_prices,estc_physical_extent) %>% left_join(.,estc_physical_dimension) %>% left_join(.,estc_publication_place) %>% subset(.,!is.na(gatherings.original) & !is.na(pagecount.orig) & !is.na(year) & !is.na(total_price)) %>% .[,c("system_control_number","total_price","pagecount","area","pagecount.plate","publication_place","year")] %>% unique.data.frame(.) %>% inner_join(.,estc_price_reliable_ids) %>% anti_join(.,estc_periodicals) %>% anti_join(.,estc_prices_unreliable_ids) %>% anti_join(.,estc_public_records)

                                                                                                                                                                                                          
#convert NA's to 0's. If plates are not found it's most likely there were no plate pages in the print product.
estc_prices_filtered$pagecount.plate[is.na(estc_prices_filtered$pagecount.plate)] <- 0

#Make final estimates for the amount of paper and plates in each print product. We remove plate pagecount from the total pagecount and assume that
#that the size of plates in a given book corresponds with the size of the books paper pages.
total_paper <- (estc_prices_filtered$pagecount-estc_prices_filtered$pagecount.plate)*(1/2)*estc_prices_filtered$area*(1/sheet_size)
total_plate <- estc_prices_filtered$pagecount.plate*(1/2)*estc_prices_filtered$area*(1/sheet_size)

#Bring the relevant transformed variables to same frame with the original data. Limit to observations with all relevant information. Presence of plates was not required.
#For the sake of possible errors in other fields, sensible value was also required from the amount of paper. 
estc_prices_filtered_2 <- cbind.data.frame(estc_prices_filtered,total_paper,total_plate) %>% subset(.,!is.na(total_paper) & !is.na(total_price) & !is.na(year) & !is.na(publication_place)) %>% subset(.,total_paper>0)


#Plot the Figure 1 of the article.
estc_prices_ag_by_year <- aggregate(estc_prices_filtered$system_control_number,by=list(estc_prices_filtered$year),FUN=length) %>% name_data_frame(.,c("Year","N"))
ggplot(estc_prices_ag_by_year)+geom_line(aes(x=Year,y=N))+theme_hsci(base_size=20)+ylab("Number of Records")+xlim(c(1680,1800))



#Limit to the wanted time frame.
estc_prices_filtered_enriched <- estc_prices_filtered_2 %>% subset(.,year>=1680 & year<=1800)

#Set seed for replication purposes.
set.seed(101)
#Split the data without to training and test data.
all_data_big_sample <- estc_prices_filtered_enriched[sample(x=1:nrow(estc_prices_filtered_enriched),size=round(nrow(estc_prices_filtered_enriched)/2)),] 

all_data_big_sample_ids <- all_data_big_sample[,c("system_control_number")] %>% as.data.frame(.) %>% name_data_frame(.,"system_control_number")

#Create training and test data sets. Standardise variables 
data_all <- estc_prices_filtered_enriched
data_all$total_paper <- scale(data_all$total_paper) %>% as.numeric(.)
data_all$total_plate <- scale(data_all$total_plate) %>% as.numeric(.)

#Create training and test data sets 

data_sample <- inner_join(data_all,all_data_big_sample_ids)

data_test_data <- anti_join(data_all,all_data_big_sample_ids)

#Model fitting and analysis 

stanData <- list(N=nrow(data_sample),N_new=nrow(data_test_data),x1=data_sample$total_paper,x2=data_sample$total_plate,x1_new=data_test_data$total_paper,x2_new=data_test_data$total_plate,y=data_sample$total_price)

stanDataFit <- stan(model_code = stan_linear_model,iter=10000,chains=2,seed=101,data=stanData)

#Get the marginal distributions of table 1 and save the Stan data to data frame.
stanData <- as.data.frame(stanDataFit)
summary(stanData$beta0)
summary(stanData$beta1)
summary(stanData$beta2)
summary(stanData$sigma)

#Produce the point estimate predictions by taking the mean over the marginal posterior distribution of each prediction.
lm_test_data_all_pred <- apply(stanData[,5:(ncol(stanData)-1)],MARGIN = 2,FUN=mean)
#lm_test_data_all_pred[lm_test_data_all_pred>300] <- 300
#Collect the prediction and standardised residuals to the same table
res <- data_test_data$total_price-lm_test_data_all_pred
pred_data <- cbind.data.frame(lm_test_data_all_pred,data_test_data$total_price,res/sd(data_test_data$total_price))
colnames(pred_data) <- c("price_predicted","price","residual")

#Compute R^2 and Relative Mean Absolute Error 
R2_all <- 1-sum((pred_data$price_predicted-pred_data$price)^2)/sum((pred_data$price-mean(pred_data$price))^2)
mae_all <- 1-mean(abs(pred_data$price_predicted-pred_data$price))/mean(abs(pred_data$price-median(pred_data$price)))

#Plots for figure 2 in the article
ggplot(data=pred_data)+geom_point(aes(x=price_predicted,y=residual))+xlab("Predicted Price")+ylab("Standardised Residual")+theme_hsci(base_size=20)
ggplot(data=pred_data)+geom_point(aes(x=price_predicted,y=residual))+xlim(c(0,200))+ylim(c(-5,5))+xlab("Predicted Price")+ylab("Standardised Residual")+theme_hsci(base_size=20)
fig2_1_1 <-  ggplot(data=pred_data)+geom_point(aes(x=price_predicted,y=residual))+xlim(c(0,100))+ylim(c(-2.5,2.5))+xlab("Predicted Price")+ylab("Standardised Residual")+theme_hsci(base_size=20)
fig2_1_2 <- ggplot(data=pred_data)+geom_point(aes(x=price_predicted,y=residual))+xlim(c(0,200))+ylim(c(-2.5,2.5))+xlab("Predicted Price")+ylab("Standardised Residual")+theme_hsci(base_size=20)
fig2_1_3 <- ggplot(data=pred_data)+geom_point(aes(x=price_predicted,y=residual))+xlab("Predicted Price")+ylab("Standardised Residual")+theme_hsci(base_size=20)

Fig2 <- grid.arrange(fig2_1_1,fig2_1_2,fig2_1_3,ncol=1)

data_test_data_w_residules <- cbind.data.frame(data_test_data,pred_data) 

#Limit the analysis to London and print products of average size.
data_test_data_w_residules_limited <- subset(data_test_data_w_residules,(publication_place=="London") & total_paper>=-0.5 & total_paper<=0.5 & total_plate>=-0.5 & total_plate<=0.5 & year>=1680 & year<=1809)

#Years to decades.  
decade <- decade_of_year(data_test_data_w_residules_limited$year)

#Add the decade variable to the residual data.
data_test_data_w_residules_limited_2 <- cbind.data.frame(data_test_data_w_residules_limited,decade)


#Production of the Figure 3

#Aggregate residuals by time and place to get the median and 1st and 3rd quantiles.
data_test_data_w_residules_limited_ag_by_place_and_time_median <- aggregate(data_test_data_w_residules_limited_2$residual,by=list(data_test_data_w_residules_limited_2$publication_place,data_test_data_w_residules_limited_2$decade),FUN=median) %>% name_data_frame(.,c("publication_place","decade","res"))
data_test_data_w_residules_limited_ag_by_place_and_time_upper <- aggregate(data_test_data_w_residules_limited_2$residual,by=list(data_test_data_w_residules_limited_2$publication_place,data_test_data_w_residules_limited_2$decade),FUN=quantile,probs=0.75) %>% name_data_frame(.,c("publication_place","decade","res_upper"))
data_test_data_w_residules_limited_ag_by_place_and_time_lower <- aggregate(data_test_data_w_residules_limited_2$residual,by=list(data_test_data_w_residules_limited_2$publication_place,data_test_data_w_residules_limited_2$decade),FUN=quantile,probs=0.25) %>% name_data_frame(.,c("publication_place","decade","res_lower"))


#Combine the residual computations to a single data set.
data_test_data_w_residules_limited_ag_by_place_and_time_full <- left_join(data_test_data_w_residules_limited_ag_by_place_and_time_median,data_test_data_w_residules_limited_ag_by_place_and_time_upper) %>% left_join(.,data_test_data_w_residules_limited_ag_by_place_and_time_lower)

#Somewhat obsolete, as the final plot was done only with London.
data_test_data_w_residules_limited_ag_by_place_and_time_full_final <- data_test_data_w_residules_limited_ag_by_place_and_time_full
data_test_data_w_residules_limited_ag_by_place_and_time_full_final$publication_place <- as.factor(data_test_data_w_residules_limited_ag_by_place_and_time_full_final$publication_place)
colnames(data_test_data_w_residules_limited_ag_by_place_and_time_full_final)[1] <- "Place"

#Visualise residuals by decade in London.
#Figure 3 of the article.
ggplot(data=data_test_data_w_residules_limited_ag_by_place_and_time_full_final)+geom_pointrange(aes(x=decade,y=res,ymin=res_lower,ymax=res_upper),alpha=0.5)+xlab("Decade")+ylab("Standardised Residual")+theme_hsci(base_size=20)+ylim(c(-0.5,0.5))
