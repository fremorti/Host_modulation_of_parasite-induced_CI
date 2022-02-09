setwd("D:/OneDrive - UGent/Postdoc Terec/Mijten Nicky/")

############################
##### Libraries ############
############################

library(ggplot2)
library(car)
library(brms)
library(tidyr)
library(dplyr)
library(grid)
library(gridExtra)
library(tibble)
library(egg)
options('max.print' =  200000)

################
# loading data #
################

data__<-read.table("Wybouwetal_CI_panel", header = 1)
colnames(data__)[c(3,9,10)]<-c('inf', 'fgeno', 'mgeno')                 #rename some columns
data__$CIobs<-(data__$adult_females/data__$eggs)                        #proportion of eggs developing into females (related to CI)
data__$MDobs<-(data__$adult_males/data__$eggs)                          #proportion of eggs developing into males (related to MD)
data__$deadeggs <- data__$eggs-data__$adult_males-data__$adult_females  #eggs that did not produce an adult female or male
data__$fdeadeggs <- data__$eggs-data__$adult_males                      #female and dead eggs
data__$FMobs<-((data__$deadeggs)/(data__$fdeadeggs))                    #proportion of dead eggs, without considering (succesful adult) males (this is a proxy for the mortality of females related to FM)
data__$intrasp <- data__$fgeno != 'Bch'                                 #Boolean indicator whether a cross was intraspecific (TRUE) or interspecific (FALSE: only Bch crosses)
levels(data__$inf) <- list("0" = "uninfected", "1" = "infected")        #renames infection levels
data__ <- data__[data__$eggs > 15,]                                         #filter control and incopatible crosses that yielded 15 or less eggs due to bigger stochasticity on proportional results and non-representativity

data_<-subset(data__, Type %in% c('I', 'C'))                            #data set without the rescue test data (Type = 'R' or 'Rc')     
data = data_[data_$intrasp,]                                            #data set with only intraspecific crosses
data <- droplevels(data)


#################
# Effects in CI #
#################

#Simplest data model. Not discussed.
#model with only Day group-level and infection effect
form0 <- brmsformula(adult_females | trials(eggs)  ~ 1 + inf + (1+inf|Day), family = binomial)
get_prior(form0, data = data)
fite0 <- brm(data = data, family = binomial,
             formula = form0,
             prior = c(prior(normal(0,5), class = Intercept),
                       prior(normal(0,5), class = b),
                       prior(cauchy(0, 1), class = sd),
                       prior(lkj(2), class = cor)),
             iter = 5000, warmup = 1000, chains = 2, cores = 2)
fite0 <- add_criterion(fite0, 'waic')
fite0$fit

nd0 <- tibble(unique(data[c('inf','mgeno', 'fgeno', 'Day')]))%>%mutate(eggs = 1) 
poste0 <- fitted(fite0, newdata = nd0, probs = c(0.09, 0.91))%>%
  as_tibble()%>%mutate(nd0[c('mgeno', 'fgeno', 'inf')])
ggplot(poste0, aes(x = fgeno, y = Estimate))+
  geom_boxplot(data = data, aes(x = fgeno, y = CIobs, color = inf), position = position_dodge(0.8), width = 0.7)+
  geom_point(aes(class = inf), position = position_dodge(0.8))+
  geom_errorbar(aes(ymin=Q9, ymax=Q91, class = inf), width=.2, position=position_dodge(0.8))+
  facet_wrap('mgeno')


# full model with FxM interactions
# all models hereafter receive a similar procedure that we will anotate for this model but which applies to all

#model formula
form2 <- brmsformula(adult_females | trials(eggs)  ~ 1 + inf*fgeno*mgeno + (1|Day), family = binomial)
get_prior(form2, data = data)
#MCMC model estimation
fite2 <- brm(data = data, family = binomial,
             formula = form2,
             #priors are defined
             prior = c(prior(normal(0,10), class = Intercept),
                       prior(normal(0,10), class = b),
                       prior(cauchy(0,2), class = sd)),
             #chains are defined
             iter = 5000, warmup = 1000, chains = 2, cores = 2)
#add waic information
fite2 <- add_criterion(fite2, 'waic')
#base table of the posterior of all estimated coefficients
fite2$fit

#define a dummy dataset of all relevant explaining variables for the model to calculate expected values from the posterior
nd <- tibble(unique(data[c('inf','mgeno', 'fgeno')]))%>%mutate(eggs = 1)
#calculate expected values from the posterior
poste2 <- fitted(fite2, newdata = nd, re_formula = NA, probs = c(0.09, 0.91))%>%
  as_tibble()%>%mutate(nd[c('mgeno', 'fgeno', 'inf')])

#plot expected values with observed values. This is mainly (but not exclusively) to check how the model fits
ggplot(poste2, aes(x = fgeno, y = Estimate))+
  geom_boxplot(data = data, aes(x = fgeno, y = CIobs, color = inf), position = position_dodge(0.8), width = 0.7)+
  geom_point(aes(class = inf), position = position_dodge(0.8))+
  geom_errorbar(aes(ymin=Q9, ymax=Q91, class = inf), width=.2, position=position_dodge(0.8))+
  labs(title = 'proportion eggs succesfully developing into adult females', y = 'proportion female', x = 'female genotype')+
  scale_color_discrete(name = "", labels = c('cured', expression(paste(italic("Wolbachia"),"-infected"))))+
  facet_wrap('mgeno')+
  theme(legend.position = c(0.8, 0.3))

#sample expected values from the posterior
poste2_ <- fitted(fite2, newdata=nd, re_formula = NA, summary = 0)%>%as_tibble()
#calculate the corrected CI metric (MD or FM in other models) from posterior distribution
postCI2__ <- 1-poste2_[seq(2,nrow(nd),2)]/poste2_[seq(1,nrow(nd),2)]%>%as_tibble()
postCI2_wide <- postCI2__%>%mutate(across(everything(), function(x) {pmax(x, 0)}))
postCI2_ <- pivot_longer(postCI2_wide, cols = everything(), names_to = 'group', names_prefix = 'V', values_to = 'CI')%>%
  mutate(group = as.numeric(group))%>%mutate(mgeno = rep(nd$mgeno[seq(1,nrow(nd),2)], nrow(postCI2__)), fgeno = rep(nd$fgeno[seq(1,nrow(nd),2)], nrow(postCI2__)))
postCI2 <- group_by(postCI2_, group)%>%summarise(CImean = mean(CI), CIq09 = quantile(CI, probs = 0.09), CIq91 = quantile(CI, probs = 0.91))%>%
  bind_cols(tibble(unique(data[c('mgeno', 'fgeno')])))

#plot estimated CI for the relevant groups in different panels for male genotype
ggplot(postCI2_, aes(x = fgeno, y = CI))+
  facet_wrap('mgeno')+
  theme_light()+
  geom_violin(draw_quantiles = c(0.09, 0.5, 0.91), alpha = 0.5, aes(fill = mgeno), scale = 'width')+
  scale_fill_brewer(palette = 'Dark2')+
  labs(title = 'Cytoplasmatic incompatibility', y = 'CI', x = 'female genotype')+
  theme(legend.position = "none")

#plot CI in one panel and oordered according to mean estimated CI 
postCI2_$fxmgeno <- paste(postCI2_$fgeno, 'X', postCI2_$mgeno)
postCI2_$fxmgeno <- factor(postCI2_$fxmgeno, levels=unique(postCI2_$fxmgeno[order(postCI2_$CI)]), ordered=TRUE)
(fig1b <- ggplot(postCI2_, aes(x = fxmgeno, y = CI))+
  theme_light()+
  geom_violin(draw_quantiles = c(0.09, 0.5, 0.91), alpha = 0.5, aes(fill = mgeno), scale = 'width')+
  coord_flip()+
  scale_fill_brewer(palette = 'Dark2')+
  labs(title = 'Cytoplasmatic incompatibility', y = 'CI', x = 'cross')+
  theme(legend.position = "none"))




#model only female in the model - model comparison
formF <- brmsformula(adult_females | trials(eggs)  ~ 1 + inf*fgeno + (1|Day), family = binomial)
get_prior(formF, data = data)
fitF <- brm(data = data, family = binomial,
             formula = formF,
             prior = c(prior(normal(0,10), class = Intercept),
                       prior(normal(0,10), class = b),
                       prior(cauchy(0, 2), class = sd)),
             iter = 5000, warmup = 1000, chains = 2, cores = 2)
fitF <- add_criterion(fitF, 'waic')
fitF$fit

ndF <- tibble(unique(data[c('inf', 'fgeno')]))%>%mutate(eggs = 1) 
postF <- fitted(fitF, re_formula = adult_females | trials(eggs)  ~ 1 + inf*fgeno, newdata = ndF, probs = c(0.09, 0.91))%>%
  as_tibble()%>%mutate(ndF[c('fgeno', 'inf')])

ggplot(postF, aes(x = fgeno, y = Estimate))+
  geom_boxplot(data = data, aes(x = fgeno, y = CIobs, color = inf), position = position_dodge(0.8), width = 0.7)+
  geom_point(aes(class = inf), position = position_dodge(0.8))+
  geom_errorbar(aes(ymin=Q9, ymax=Q91, class = inf), width=.2, position=position_dodge(0.8))+
  labs(title = 'proportion eggs succesfully developing into adult females', y = 'proportion female', x = 'female genotype')+
  scale_color_discrete(name = expression(italic("Wolbachia")), labels = c('no', 'yes'))+
  theme(legend.position = 'bottom')


#model only male  - model comparison
formM <- brmsformula(adult_females | trials(eggs)  ~ 1 + inf*mgeno + (1|Day), family = binomial)
get_prior(formM, data = data)
fitM <- brm(data = data, family = binomial,
             formula = formM,
             prior = c(prior(normal(0,10), class = Intercept),
                       prior(normal(0,10), class = b),
                       prior(cauchy(0, 2), class = sd)),
             iter = 5000, warmup = 1000, chains = 2, cores = 2)
fitM <- add_criterion(fitM, 'waic')
fitM$fit

ndM <- tibble(unique(data[c('inf', 'mgeno')]))%>%mutate(eggs = 1) 
postM <- fitted(fitM, re_formula = adult_females | trials(eggs)  ~ 1 + inf*mgeno, newdata = ndM, probs = c(0.09, 0.91))%>%
  as_tibble()%>%mutate(ndM[c('mgeno', 'inf')])

ggplot(postM, aes(x = mgeno, y = Estimate))+
  geom_boxplot(data = data, aes(x = mgeno, y = CIobs, color = inf), position = position_dodge(0.8), width = 0.7)+
  geom_point(aes(class = inf), position = position_dodge(0.8))+
  geom_errorbar(aes(ymin=Q9, ymax=Q91, class = inf), width=.2, position=position_dodge(0.8))+
  labs(title = 'proportion eggs succesfully developing into adult females', y = 'proportion female', x = 'male genotype')+
  scale_color_discrete(name = expression(italic("Wolbachia")), labels = c('no', 'yes'))+
  theme(legend.position = 'bottom')

#model female and male but no interaction  - model comparison
formMF <- brmsformula(adult_females | trials(eggs)  ~ 1 + inf*mgeno + inf*fgeno + (1|Day), family = binomial)
get_prior(formMF, data = data)
fitMF <- brm(data = data, family = binomial,
            formula = formMF,
            prior = c(prior(normal(0,10), class = Intercept),
                      prior(normal(0,10), class = b),
                      prior(cauchy(0, 2), class = sd)),
            iter = 5000, warmup = 1000, chains = 2, cores = 2)
fitMF <- add_criterion(fitMF, 'waic')
fitMF$fit

ndMF <- tibble(unique(data[c('inf', 'mgeno', 'fgeno')]))%>%mutate(eggs = 1) 
postMF <- fitted(fitMF, re_formula = adult_females | trials(eggs)  ~ 1 + inf*mgeno+inf*fgeno, newdata = ndMF, probs = c(0.09, 0.91))%>%
  as_tibble()%>%mutate(ndMF[c('fgeno', 'mgeno', 'inf')])

ggplot(postMF, aes(x = fgeno, y = Estimate))+
  geom_boxplot(data = data, aes(x = fgeno, y = CIobs, color = inf), position = position_dodge(0.8), width = 0.7)+
  geom_point(aes(class = inf), position = position_dodge(0.8))+
  geom_errorbar(aes(ymin=Q9, ymax=Q91, class = inf), width=.2, position=position_dodge(0.8))+
  facet_wrap('mgeno')+
  labs(title = 'proportion eggs succesfully developing into adult females', y = 'proportion female', x = 'female genotype')+
  scale_color_discrete(name = expression(italic("Wolbachia")), labels = c('no', 'yes'))+
  theme(legend.position = c(0.8, 0.3))

#model comparison CI
w <- loo_compare(fitF, fitM, fite2, fitMF, criterion = 'waic')
print(w, simplify = 0)
w[, 7:8] %>% 
  data.frame() %>% 
  rownames_to_column(var = "model_name") %>% 
  ggplot(aes(x    = model_name, y    = waic, ymin = waic - se_waic, ymax = waic + se_waic)) +
  geom_pointrange(shape = 21) +
  coord_flip() +
  scale_x_discrete(labels = c('full', 'F', 'M', 'M+F'))+
  labs(x = NULL, y = 'WAIC', title = "WAIC plot") +
  theme_classic()



#male, female genotype and interaction as variable effect
form2_ <- brmsformula(adult_females | trials(eggs)  ~ 1 + inf + (1+inf|mgeno*fgeno) + (1|Day), family = binomial)
get_prior(form2_, data = data)
fite2_ <- brm(data = data, family = binomial,
              formula = form2_,
              prior = c(prior(normal(0,5), class = Intercept),
                        prior(normal(0,5), class = b),
                        prior(exponential(1), class = sd),
                        prior(lkj(2), class = cor)),
              iter = 5000, warmup = 1000, chains = 2, cores = 2)
fite2_ <- add_criterion(fite2_, 'waic')
fite2_$fit

nd2_ <- tibble(unique(data[c('inf','mgeno', 'fgeno')]))%>%mutate(eggs = 1) 
poste2_ <- fitted(fite2_, re_formula = adult_females | trials(eggs)  ~ 1 + inf*mgeno + (1+inf*mgeno|fgeno), newdata = nd2_, probs = c(0.09, 0.91))%>%
  as_tibble()%>%mutate(nd2_[c('mgeno', 'fgeno', 'inf')])

ggplot(poste2_, aes(x = fgeno, y = Estimate))+
  geom_boxplot(data = data, aes(x = fgeno, y = CIobs, color = inf), position = position_dodge(0.8), width = 0.7)+
  geom_point(aes(class = inf), position = position_dodge(0.8))+
  geom_errorbar(aes(ymin=Q9, ymax=Q91, class = inf), width=.2, position=position_dodge(0.8))+
  labs(title = 'proportion eggs succesfully developing into adult females', y = 'proportion female', x = 'female genotype')+
  scale_color_discrete(name = expression(italic("Wolbachia")), labels = c('no', 'yes'))+
  facet_wrap('mgeno')+
  theme(legend.position = c(0.8, 0.3))

posta_ <- fitted(fite2_, newdata=nd2_, re_formula = adult_females | trials(eggs)  ~ 1 + inf + (1+inf|mgeno*fgeno), summary = 0)%>%as_tibble()
postCIa__ <- 1-posta_[seq(2,nrow(nd2_),2)]/posta_[seq(1,nrow(nd2_),2)]%>%as_tibble()
postCIa_wide <- postCIa__%>%mutate(across(everything(), function(x) {pmax(x, 0)}))
postCIa_ <- pivot_longer(postCIa_wide, cols = everything(), names_to = 'group', names_prefix = 'V', values_to = 'CI')%>%
  mutate(group = as.numeric(group))%>%mutate(mgeno = rep(nd2_$mgeno[seq(1,nrow(nd2_),2)], nrow(postCIa__)), fgeno = rep(nd2_$fgeno[seq(1,nrow(nd2_),2)], nrow(postCIa__)))
postCIa <- group_by(postCIa_, group)%>%summarise(mean = mean(CI), q09 = quantile(CI, probs = 0.09), q91 = quantile(CI, probs = 0.91))



ggplot(postCIa_, aes(x = fgeno, y = CI))+
  geom_violin(draw_quantiles = c(0.09, 0.5, 0.91), alpha = 0.3, aes(fill = fgeno), scale = 'width')+
  facet_wrap('mgeno')+
  theme_light()+
  theme(legend.position = "none")


#visualize infection x female x male interaction
paramsinfint <- fite2_%>%posterior_samples()%>%select(starts_with('r_mgeno:fgeno[') & ends_with('inf1]'))%>%
  rename_with(function(x) substr(x, 15, nchar(x)-6))%>%
  pivot_longer(cols = everything(), names_to = 'fmgeno_interaction', values_to = 'estimate')%>%
  mutate(mgeno = matrix(unlist(strsplit(.$fmgeno_interaction, '_')), ncol = 2, byrow = 1)[,1], 
         fgeno = matrix(unlist(strsplit(.$fmgeno_interaction, '_')), ncol = 2, byrow = 1)[,2])

ggplot(paramsinfint, aes(x = fgeno, y = estimate, fill = mgeno))+
  geom_violin(draw_quantiles = c(0.09, 0.5, 0.91), alpha = 0.3, scale = 'width')+
  facet_wrap('mgeno')+
  theme_light()+
  theme(legend.position = "none")+
  labs(title = 'Infection effect in different male x female crosses' , y = 'infection effect on log-odds', x = 'female genotype')+
  scale_fill_brewer(palette = 'Dark2')


#variance analysis
#Expected values on log-odds scale
fitted1 <- fitted(fite2_, summary = 0, scale = 'linear')
#Residuals on log-odds scale
observed_logit <- logit(data$CIobs, adjust = 0.0025)                   
resid <- sweep(as.matrix(fitted1), 2, observed_logit, '-')%>%as_tibble()%>%rowwise()%>%transmute(residuals = sd(c_across(everything())))

#calculate finite-population standard deviation of interactions with infection
pi_ <- as_tibble(fite2_)%>%rowwise()%>%
  transmute(sd_mgenoi = sd(c_across(starts_with('r_mgeno[') & ends_with('inf1]'))),
            sd_fgenoi = sd(c_across(starts_with('r_fgeno[') & ends_with('inf1]'))),
            sd_mfi = sd(c_across(starts_with('r_mgeno:fgeno[') & ends_with('inf1]'))))%>%
  bind_cols(resid$residuals)%>%rename(resid = ...4)%>%
  mutate(dm_fi = sd_mgenoi-sd_fgenoi, dm_mfi = sd_mgenoi-sd_mfi, df_mfi = sd_mfi-sd_fgenoi)

#calculate finite-population standard deviation  of levels of intercepts (no infection effect)
p_ <- as_tibble(fite2_)%>%rowwise()%>%
  transmute(sd_mgeno = sd(c_across(starts_with('r_mgeno[') & ends_with('Intercept]'))),
            sd_fgeno = sd(c_across(starts_with('r_fgeno[') & ends_with('Intercept]'))),
            sd_mf = sd(c_across(starts_with('r_mgeno:fgeno[') & ends_with('Intercept]'))))%>%
  mutate(dm_f = sd_mgeno-sd_fgeno, dm_mf = sd_mgeno-sd_mf, df_mf = sd_mf-sd_fgeno)
ppi <- cbind(p_[c('sd_mgeno', 'sd_fgeno', 'sd_mf')], pi_[c('sd_mgenoi', 'sd_fgenoi', 'sd_mfi', 'resid')])%>%
  rename(fgeno = sd_fgeno, mgeno = sd_mgeno, 'fgeno:inf' = sd_fgenoi, 'mgeno:inf' = sd_mgenoi, 'm:f' = sd_mf, 'm:f:inf' = sd_mfi)%>%
  pivot_longer(.[c('mgeno:inf', 'fgeno:inf', 'm:f:inf', 'resid', 'mgeno', 'fgeno', 'm:f')], cols = everything(), names_to = 'level', values_to = 'sd')
ppi$level <- factor(ppi$level, levels = c('fgeno', 'fgeno:inf','mgeno','mgeno:inf',  'm:f' , 'm:f:inf', 'resid'))

#plot sd of all levels
ggplot(ppi, aes(x = level, y = sd))+
  geom_violin(draw_quantiles = c(0.09, 0.5, 0.91), alpha = 0.3, aes(fill = level), scale = 'width')+
  scale_fill_brewer(palette = 'Dark2')+
  labs(title = 'sd of model level estimates', y = 'sd', x = 'level')+
  theme(legend.position = 'none')



###############################
# CI with Bch as different sp #
###############################

data_i <- data_
#model with intraspecific cross as explanatory variable, but no female genotype (essentially making the distinction between crosses involving females of Bch and females of other lines)
form3_ <- brmsformula(adult_females | trials(eggs)  ~ 1 + inf*mgeno*intrasp + (1|Day), family = binomial)
get_prior(form3_, data = data_i)
fite3_ <- brm(data = data_i, family = binomial,
              formula = form3_,
              prior = c(prior(normal(0,10), class = Intercept),
                        prior(normal(0,10), class = b),
                        prior(cauchy(0, 2), class = sd)),
              iter = 5000, warmup = 1000, chains = 2, cores = 2)
fite3_ <- add_criterion(fite3_, 'waic')
fite3_$fit

nd3_ <- tibble(unique(data_i[c('inf','mgeno', 'intrasp')]))%>%mutate(eggs = 1) 
poste3_ <- fitted(fite3_, re_formula = adult_females | trials(eggs)  ~ 1 + inf*mgeno*intrasp, newdata = nd3_, probs = c(0.09, 0.91))%>%
  as_tibble()%>%mutate(nd3_[c('mgeno', 'intrasp', 'inf')])

ggplot(poste3_, aes(x = intrasp, y = Estimate))+
  geom_boxplot(data = data_, aes(x = intrasp, y = CIobs, color = inf), position = position_dodge(0.8), width = 0.7)+
  geom_point(aes(class = inf), position = position_dodge(0.8))+
  geom_errorbar(aes(ymin=Q9, ymax=Q91, class = inf), width=.2, position=position_dodge(0.8))+
  labs(title = 'proportion eggs succesfully developing into adult females', y = 'proportion female', x = 'intraspecific cross')+
  scale_color_discrete(name = expression(italic("Wolbachia")), labels = c('no', 'yes'))+
  facet_wrap('mgeno')+
  theme(legend.position = c(0.8, 0.3))

#CI with interactions
post3_ <- fitted(fite3_, newdata=nd3_, re_formula = adult_females | trials(eggs)  ~ 1 + inf*mgeno*intrasp, summary = 0)%>%as_tibble()
postCI3__ <- 1-post3_[seq(2,nrow(nd3_),2)]/post3_[seq(1,nrow(nd3_),2)]%>%as_tibble()
postCI3_wide <- postCI3__%>%mutate(across(everything(), function(x) {pmax(x, 0)}))
postCI3_ <- pivot_longer(postCI3_wide, cols = everything(), names_to = 'group', names_prefix = 'V', values_to = 'CI')%>%
  mutate(group = as.numeric(group))%>%mutate(mgeno = rep(nd3_$mgeno[seq(1,nrow(nd3_),2)], nrow(postCI3__)), intrasp = rep(nd3_$intrasp[seq(1,nrow(nd3_),2)], nrow(postCI3__)))
postCI3 <- group_by(postCI3_, group)%>%summarise(mean = mean(CI), q09 = quantile(CI, probs = 0.09), q91 = quantile(CI, probs = 0.91))%>%mutate(intrasp = nd3_$intrasp[seq(1,nrow(nd3_),2)])



postCI3_$mgeno <- factor(postCI3_$mgeno, levels=unique(postCI3_$mgeno[order(postCI3[postCI3$intrasp,]$mean)]), ordered=TRUE)
(fig1c <-ggplot(postCI3_, aes(x = mgeno, y = CI))+
  geom_violin(draw_quantiles = c(0.09, 0.5, 0.91), alpha = 0.3, aes(fill = intrasp), scale = 'width', position = position_dodge(0.8))+
  theme_light()+
  coord_flip()+
  labs(title = 'Effect of genetic distance on CI', y = 'CI', x = 'male genotype', fill = 'intraspecific cross')+
  theme(legend.position = 'bottom'))
#define fixed plot size
fig1bset <- set_panel_size(fig1b,
                           width  = unit(8, "cm"),
                           height = unit(20, "cm"))
fig1cset <- set_panel_size(fig1c,
                           width  = unit(8, "cm"),
                           height = unit(12, "cm"))
grid.newpage()
grid.draw(fig1bset) #rescale png save to height 825, keeping aspect ratio
grid.newpage()
grid.draw(fig1cset)

#model to generate the distribution of T. urticae fgeno effects to see whether Bch falls within the estimated range of effects or not
#female genotype and interaction as variable effect to simulate the population of female genotypes
form2_b <- brmsformula(adult_females | trials(eggs)  ~ 1 + inf*mgeno + (1+inf*mgeno|fgeno) + (1|Day), family = binomial)
get_prior(form2_b, data = data)
fite2_b <- brm(data = data, family = binomial,
               formula = form2_b,
               prior = c(prior(normal(0,5), class = Intercept),
                         prior(normal(0,5), class = b),
                         prior(exponential(1), class = sd),
                         prior(lkj(2), class = cor)),
               iter = 5000, warmup = 1000, chains = 2, cores = 2)
fite2_b <- add_criterion(fite2_b, 'waic')
fite2_b$fit

nd2b_ <- tibble(unique(data[c('inf','mgeno')]))%>%mutate(eggs = 1) 
#allow_new_levels argument asks the fitted function to simulate posterior predictions of new group levels (here: fgeno)
poste2b_ <- fitted(fite2_b, probs = c(0.09, 0.91), re_formula = adult_females | trials(eggs)  ~ 1 + inf*mgeno + (1+inf*mgeno|fgeno), newdata=nd2b_, allow_new_levels = T,
                   sample_new_levels = "gaussian")%>%as_tibble()%>%mutate(nd2b_[c('inf', 'mgeno')])

ggplot(poste2b_, aes(x = mgeno, y = Estimate))+
  geom_boxplot(data = data, aes(x = mgeno, y = CIobs, color = inf), position = position_dodge(0.8), width = 0.7)+
  geom_point(aes(class = inf), position = position_dodge(0.8))+
  geom_errorbar(aes(ymin=Q9, ymax=Q91, class = inf), width=.2, position=position_dodge(0.8))+
  labs(title = 'proportion eggs succesfully developing into adult females', y = 'proportion female', x = 'female genotype')+
  scale_color_discrete(name = expression(italic("Wolbachia")), labels = c('no', 'yes'))+
  theme(legend.position = c(0.8, 0.3))

postb_ <- fitted(fite2_b,re_formula = adult_females | trials(eggs)  ~ 1 + inf*mgeno + (1+inf*mgeno|fgeno), newdata=nd2b_, summary = 0, allow_new_levels = T,
                 sample_new_levels = "gaussian")%>%as_tibble()
postCIb__ <- 1-postb_[seq(2,nrow(nd2b_),2)]/postb_[seq(1,nrow(nd2b_),2)]%>%as_tibble()
postCIb_wide <- postCIb__%>%mutate(across(everything(), function(x) {pmax(x, 0)}))
postCIb_ <- pivot_longer(postCIb_wide, cols = everything(), names_to = 'group', names_prefix = 'V', values_to = 'CI')%>%
  mutate(group = as.numeric(group))%>%mutate(mgeno = rep(nd2b_$mgeno[seq(1,nrow(nd2b_),2)], nrow(postCIb__)))
postCIb <- group_by(postCIb_, group)%>%summarise(mean = mean(CI), q09 = quantile(CI, probs = 0.09), q91 = quantile(CI, probs = 0.91))

ggplot(postCIb_, aes(x = mgeno, y = CI))+
  theme_light()+
  geom_violin(draw_quantiles = c(0.09, 0.5, 0.91), alpha = 0.5, aes(fill = mgeno), scale = 'width')+
  scale_fill_brewer(palette = 'Dark2')+
  labs(title = 'Cytoplasmatic incompatibility', y = 'CI', x = 'male genotype')+
  theme(legend.position = "none")



simfgeno <- postCIb_%>%mutate(intrasp = TRUE)%>%bind_rows(postCI3_[!postCI3_$intrasp, c('group', 'CI', 'intrasp', 'mgeno')]%>%mutate(mgeno = factor(mgeno, ordered = 0)))
ggplot(simgeno, aes(x = mgeno, y = CI))+
  geom_violin(draw_quantiles = c(0.09, 0.5, 0.91), alpha = 0.3, aes(fill = intrasp), scale = 'width')+
  theme_light()+
  coord_flip()+
  labs(title = 'Bch effect on CI compared to ditribution of T. urticae effects', y = 'CI', x = 'male genotype', fill = 'intraspecific cross')+
  theme(legend.position = 'bottom')


######################
#  Male development  #
######################

#Full model for MD
formMD <- brmsformula(adult_males | trials(eggs)  ~ 1 + inf*fgeno*mgeno + (1|Day), family = binomial)
get_prior(formMD, data = data)
fiteMD <- brm(data = data, family = binomial,
             formula = formMD,
             prior = c(prior(normal(0,10), class = Intercept),
                       prior(normal(0,10), class = b),
                       prior(cauchy(0, 2), class = sd)),
             iter = 5000, warmup = 1000, chains = 2, cores = 2)
fiteMD <- add_criterion(fiteMD, 'waic')
fiteMD$fit

ndMD <- tibble(unique(data[c('inf','mgeno', 'fgeno')]))%>%mutate(eggs = 1) 
postMDish2 <- fitted(fiteMD, newdata = ndMD, re_formula = NA, probs = c(0.09, 0.91))%>%
  as_tibble()%>%mutate(ndMD[c('mgeno', 'fgeno', 'inf')])

ggplot(postMDish2, aes(x = fgeno, y = Estimate))+
  geom_boxplot(data = data, aes(x = fgeno, y = MDobs, color = inf), position = position_dodge(0.8), width = 0.7)+
  geom_point(aes(class = inf), position = position_dodge(0.8))+
  geom_errorbar(aes(ymin=Q9, ymax=Q91, class = inf), width=.2, position=position_dodge(0.8))+
  labs(title = 'proportion eggs succesfully developing into adult males', y = 'proportion male', x = 'female genotype')+
  scale_color_discrete(name = expression(italic("Wolbachia")), labels = c('no', 'yes'))+
  facet_wrap('mgeno')+
  theme(legend.position = c(0.8, 0.3))


postMDf_ <- fitted(fiteMD, newdata=ndMD, re_formula = NA, summary = 0)%>%as_tibble()
postMD__ <- (postMDf_[seq(2,nrow(ndMD),2)]-postMDf_[seq(1,nrow(ndMD),2)])/(1-postMDf_[seq(1,nrow(ndMD),2)])%>%as_tibble()
postMD_wide <- postMD__%>%mutate(across(everything(), function(x) {pmax(x, 0)}))
postMD_ <- pivot_longer(postMD_wide, cols = everything(), names_to = 'group', names_prefix = 'V', values_to = 'MD')%>%
  mutate(group = as.numeric(group))%>%mutate(mgeno = rep(ndMD$mgeno[seq(1,nrow(ndMD),2)], nrow(postMD__)), fgeno = rep(ndMD$fgeno[seq(1,nrow(ndMD),2)], nrow(postMD__)))
postMD <- group_by(postMD_, group)%>%summarise(MDmean = mean(MD), MDq09 = quantile(MD, probs = 0.09), MDq91 = quantile(MD, probs = 0.91))%>%
  bind_cols(tibble(unique(data[c('mgeno', 'fgeno')])))

ggplot(postMD_, aes(x = fgeno, y = MD))+
  geom_violin(draw_quantiles = c(0.09, 0.5, 0.91), alpha = 0.3, aes(fill = mgeno), scale = 'width')+
  facet_wrap('mgeno')+
  theme_light()+
  scale_fill_brewer(palette = 'Dark2')+
  labs(title = 'Male development', y = 'MD', x = 'female genotype')+
  theme(legend.position = "none")

#plot MD all together, ordered according to mean CI
postMD_$fxmgeno <- paste(postMD_$fgeno, 'X', postMD_$mgeno)
postMD_$fxmgeno <- factor(postMD_$fxmgeno, levels=unique(postCI2_$fxmgeno[order(postCI2$CImean)]), ordered=TRUE)
ggplot(postMD_, aes(x = fxmgeno, y = MD))+
  geom_violin(draw_quantiles = c(0.09, 0.5, 0.91), alpha = 0.3, aes(fill = mgeno), scale = 'width')+
  coord_flip()+
  theme_light()+
  scale_fill_brewer(palette = 'Dark2')+
  labs(title = 'Male development', y = 'MD', x = 'female genotype')+
  theme(legend.position = "none")

#plot correlation of MD and CI
corrCIMD <- full_join(postCI2, postMD)
ggplot(corrCIMD, aes(x = CImean, y = MDmean))+
  geom_point()+
  geom_errorbar(aes(ymin=MDq09, ymax=MDq91), width=.01)+
  geom_errorbarh(aes(xmin=CIq09, xmax=CIq91), height=.01)+
  theme_light()+
  lims(x = c(0,1), y = c(0,1))+
  scale_fill_brewer(palette = 'Dark2')+
  labs(title = 'CI-MD correlation', y = 'MD', x = 'CI')+
  theme(legend.position = "none")




#interspecific comparison of MD (Bch - rest)
formMD2 <- brmsformula(adult_males | trials(eggs)  ~ 1 + inf*intrasp*mgeno + (1|Day), family = binomial)
get_prior(formMD2, data = data_i)
fiteMD2 <- brm(data = data_i, family = binomial,
              formula = formMD2,
              prior = c(prior(normal(0,10), class = Intercept),
                        prior(normal(0,10), class = b),
                        prior(cauchy(0, 2), class = sd)),
              iter = 5000, warmup = 1000, chains = 2, cores = 2)
fiteMD2 <- add_criterion(fiteMD2, 'waic')
fiteMD2$fit

nd <- tibble(unique(data_i[c('inf','mgeno', 'intrasp')]))%>%mutate(eggs = 1) 
postMD22 <- fitted(fiteMD2, newdata = nd, re_formula = NA, probs = c(0.09, 0.91))%>%
  as_tibble()%>%mutate(nd[c('mgeno', 'intrasp', 'inf')])

ggplot(postMD22, aes(x = intrasp, y = Estimate))+
  geom_boxplot(data = data_i, aes(x = intrasp, y = MDobs, color = inf), position = position_dodge(0.8), width = 0.7)+
  geom_point(aes(class = inf), position = position_dodge(0.8))+
  geom_errorbar(aes(ymin=Q9, ymax=Q91, class = inf), width=.2, position=position_dodge(0.8))+
  labs(title = 'proportion eggs succesfully developing into adult males', y = 'proportion male', x = 'intraspecific cross')+
  scale_color_discrete(name = expression(italic("Wolbachia")), labels = c('no', 'yes'))+
  facet_wrap('mgeno')+
  theme(legend.position = c(0.8, 0.3))


postMD22_ <- fitted(fiteMD2, newdata=nd, re_formula = NA, summary = 0)%>%as_tibble()
postMD2__ <- (postMD22_[seq(2,nrow(nd),2)]-postMD22_[seq(1,nrow(nd),2)])/(1-postMD22_[seq(1,nrow(nd),2)])%>%as_tibble()
postMD2_wide <- postMD2__%>%mutate(across(everything(), function(x) {pmax(x, 0)}))
postMD2_ <- pivot_longer(postMD2_wide, cols = everything(), names_to = 'group', names_prefix = 'V', values_to = 'MD')%>%
  mutate(group = as.numeric(group))%>%mutate(mgeno = rep(nd$mgeno[seq(1,nrow(nd),2)], nrow(postMD2__)), intrasp = rep(nd$intrasp[seq(1,nrow(nd),2)], nrow(postMD2__)))
postMD2 <- group_by(postMD2_, group)%>%summarise(mean = mean(MD), q09 = quantile(MD, probs = 0.09), q91 = quantile(MD, probs = 0.91))

#plot intrasp-intersp comparison
(figMDb <- ggplot(postMD2_, aes(x = mgeno, y = MD, fill = intrasp))+
  geom_violin(draw_quantiles = c(0.09, 0.5, 0.91), alpha = 0.3, aes(fill = intrasp), scale = 'width')+
  theme_light()+
  coord_flip()+
  scale_fill_brewer(palette = 'Dark2')+
  labs(title = 'Male development', y = 'MD', x = 'female genotype', fill ='intrasp. cross'))


#model comparison
#only males
formMDM <- brmsformula(adult_males | trials(eggs)  ~ 1 + inf*mgeno + (1|Day), family = binomial)
get_prior(formMDM, data = data)
fiteMDM <- brm(data = data, family = binomial,
               formula = formMDM,
               prior = c(prior(normal(0,10), class = Intercept),
                         prior(normal(0,10), class = b),
                         prior(cauchy(0, 2), class = sd)),
               iter = 5000, warmup = 1000, chains = 2, cores = 2)
fiteMDM <- add_criterion(fiteMDM, 'waic')
fiteMDM$fit

ndMDM <- tibble(unique(data[c('inf', 'mgeno')]))%>%mutate(eggs = 1) 
postMDM <- fitted(fiteMDM, re_formula = adult_males | trials(eggs)  ~ 1 + inf*mgeno, newdata = ndMDM, probs = c(0.09, 0.91))%>%
  as_tibble()%>%mutate(ndMDM[c('mgeno', 'inf')])

ggplot(postMDM, aes(x = mgeno, y = Estimate))+
  geom_boxplot(data = data, aes(x = mgeno, y = MDobs, color = inf), position = position_dodge(0.8), width = 0.7)+
  geom_point(aes(class = inf), position = position_dodge(0.8))+
  geom_errorbar(aes(ymin=Q9, ymax=Q91, class = inf), width=.2, position=position_dodge(0.8))+
  labs(title = 'proportion eggs succesfully developing into adult males', y = 'proportion male', x = 'male genotype')+
  scale_color_discrete(name = expression(italic("Wolbachia")), labels = c('no', 'yes'))+
  theme(legend.position = 'bottom')


#only female
formMDF <- brmsformula(adult_males | trials(eggs)  ~ 1 + inf*fgeno + (1|Day), family = binomial)
get_prior(formMDF, data = data)
fiteMDF <- brm(data = data, family = binomial,
               formula = formMDF,
               prior = c(prior(normal(0,10), class = Intercept),
                         prior(normal(0,10), class = b),
                         prior(cauchy(0, 2), class = sd)),
               iter = 5000, warmup = 1000, chains = 2, cores = 2)
fiteMDF <- add_criterion(fiteMDF, 'waic')
fiteMDF$fit
ndMDF <- tibble(unique(data[c('inf', 'fgeno')]))%>%mutate(eggs = 1) 
postMDF <- fitted(fiteMDF, re_formula = adult_males | trials(eggs)  ~ 1 + inf*fgeno, newdata = ndMDF, probs = c(0.09, 0.91))%>%
  as_tibble()%>%mutate(ndMDF[c('fgeno', 'inf')])

ggplot(postMDF, aes(x = fgeno, y = Estimate))+
  geom_boxplot(data = data, aes(x = fgeno, y = MDobs, color = inf), position = position_dodge(0.8), width = 0.7)+
  geom_point(aes(class = inf), position = position_dodge(0.8))+
  geom_errorbar(aes(ymin=Q9, ymax=Q91, class = inf), width=.2, position=position_dodge(0.8))+
  labs(title = 'proportion eggs succesfully developing into adult males', y = 'proportion male', x = 'female genotype')+
  scale_color_discrete(name = expression(italic("Wolbachia")), labels = c('no', 'yes'))+
  theme(legend.position = 'bottom')

#males and females, no interaction
formMDFM <- brmsformula(adult_males | trials(eggs)  ~ 1 + inf*mgeno + inf*fgeno + (1|Day), family = binomial)
get_prior(formMDFM, data = data)
fiteMDFM <- brm(data = data, family = binomial,
               formula = formMDFM,
               prior = c(prior(normal(0,10), class = Intercept),
                         prior(normal(0,10), class = b),
                         prior(cauchy(0, 2), class = sd)),
               iter = 5000, warmup = 1000, chains = 2, cores = 2)
fiteMDFM <- add_criterion(fiteMDFM, 'waic')
fiteMDFM$fit

ndMDFM <- tibble(unique(data[c('inf', 'fgeno', 'mgeno')]))%>%mutate(eggs = 1) 
postMDFM <- fitted(fiteMDFM, re_formula = adult_males | trials(eggs)  ~ 1 + inf*fgeno + inf*mgeno, newdata = ndMDFM, probs = c(0.09, 0.91))%>%
  as_tibble()%>%mutate(ndMDFM[c('fgeno', 'mgeno', 'inf')])

ggplot(postMDFM, aes(x = fgeno, y = Estimate))+
  geom_boxplot(data = data, aes(x = fgeno, y = MDobs, color = inf), position = position_dodge(0.8), width = 0.7)+
  geom_point(aes(class = inf), position = position_dodge(0.8))+
  geom_errorbar(aes(ymin=Q9, ymax=Q91, class = inf), width=.2, position=position_dodge(0.8))+
  labs(title = 'proportion eggs succesfully developing into adult males', y = 'proportion male', x = 'female genotype')+
  facet_wrap('mgeno')+
  scale_color_discrete(name = expression(italic("Wolbachia")), labels = c('no', 'yes'))+
  theme(legend.position = c(0.8, 0.3))

#compare all WAICs
w <- loo_compare(fiteMDM, fiteMDF, fiteMDFM, fiteMD, criterion = 'waic')
print(w, simplify = 0)
w[, 7:8] %>% 
  data.frame() %>% 
  rownames_to_column(var = "model_name") %>% 
  ggplot(aes(x    = model_name, y    = waic, ymin = waic - se_waic, ymax = waic + se_waic)) +
  geom_pointrange(shape = 21) +
  coord_flip() +
  scale_x_discrete(labels = c('full', 'F', 'M+F', 'M'))+
  labs(x = NULL, y = 'WAIC', title = "WAIC plot") +
  theme_classic()



#variance analysis MD
#full model but with mgeno and fgeno as variable effects
formMD_ <- brmsformula(adult_males | trials(eggs)  ~ 1 + inf + (1+inf|fgeno*mgeno) + (1|Day), family = binomial)
get_prior(formMD_, data = data)
fiteMD_ <- brm(data = data, family = binomial,
              formula = formMD_,
              prior = c(prior(normal(0,5), class = Intercept),
                        prior(normal(0,5), class = b),
                        prior(cauchy(0, 1), class = sd),
                        prior(lkj(2), class = cor)),
              iter = 5000, warmup = 1000, chains = 2, cores = 2)
fiteMD_ <- add_criterion(fiteMD_, 'waic')
fiteMD_$fit


ndMD_ <- tibble(unique(data[c('inf','mgeno', 'fgeno', 'Day')]))%>%mutate(eggs = 1) 
postMD_f <- fitted(fiteMD_, re_formula = adult_males | trials(eggs)  ~ 1 + inf + (1+inf|mgeno*fgeno), newdata = ndMD_, probs = c(0.09, 0.91))%>%
  as_tibble()%>%mutate(ndMD_[c('mgeno', 'fgeno', 'inf')])

ggplot(postMD_f, aes(x = fgeno, y = Estimate))+
  geom_boxplot(data = data, aes(x = fgeno, y = MDobs, color = inf), position = position_dodge(0.8), width = 0.7)+
  geom_point(aes(class = inf), position = position_dodge(0.8))+
  geom_errorbar(aes(ymin=Q9, ymax=Q91, class = inf), width=.2, position=position_dodge(0.8))+
  labs(title = 'proportion eggs succesfully developing into adult males', y = 'proportion male', x = 'female genotype')+
  scale_color_discrete(name = expression(italic("Wolbachia")), labels = c('no', 'yes'))+
  facet_wrap('mgeno')+
  theme(legend.position = c(0.8, 0.3))


fittedMD1 <- fitted(fiteMD_, summary = 0, scale = 'linear')
observed_logit <- logit(data$MDobs, adjust = 0.0025)                   
residMD <- sweep(as.matrix(fittedMD1), 2, observed_logit, '-')%>%as_tibble()%>%rowwise()%>%transmute(residuals = sd(c_across(everything())))

#calculate finite-population standard deviation of interactions with infection
pi_MD <- as_tibble(fiteMD_)%>%rowwise()%>%
  transmute(sd_mgenoi = sd(c_across(starts_with('r_mgeno[') & ends_with('inf1]'))),
            sd_fgenoi = sd(c_across(starts_with('r_fgeno[') & ends_with('inf1]'))),
            sd_mfi = sd(c_across(starts_with('r_fgeno:mgeno[') & ends_with('inf1]'))))%>%
  bind_cols(residMD$residuals)%>%rename(resid = ...4)%>%
  mutate(dm_fi = sd_mgenoi-sd_fgenoi, dm_mfi = sd_mgenoi-sd_mfi, df_mfi = sd_mfi-sd_fgenoi)

#finite-population standard deviation of intercepts (no infection effect)
p_MD <- as_tibble(fiteMD_)%>%rowwise()%>%
  transmute(sd_mgeno = sd(c_across(starts_with('r_mgeno[') & ends_with('Intercept]'))),
            sd_fgeno = sd(c_across(starts_with('r_fgeno[') & ends_with('Intercept]'))),
            sd_mf = sd(c_across(starts_with('r_fgeno:mgeno[') & ends_with('Intercept]'))))%>%
  mutate(dm_f = sd_mgeno-sd_fgeno, dm_mf = sd_mgeno-sd_mf, df_mf = sd_mf-sd_fgeno)
ppiMD <- cbind(p_MD[c('sd_mgeno', 'sd_fgeno', 'sd_mf')], pi_MD[c('sd_mgenoi', 'sd_fgenoi', 'sd_mfi', 'resid')])%>%
  rename(fgeno = sd_fgeno, mgeno = sd_mgeno, 'fgeno:inf' = sd_fgenoi, 'mgeno:inf' = sd_mgenoi, 'm:f' = sd_mf, 'm:f:inf' = sd_mfi)%>%
  pivot_longer(.[c('mgeno:inf', 'fgeno:inf', 'm:f:inf', 'resid', 'mgeno', 'fgeno', 'm:f')], cols = everything(), names_to = 'level', values_to = 'sd')
ppiMD$level <- factor(ppiMD$level, levels = c('fgeno', 'fgeno:inf','mgeno','mgeno:inf',  'm:f' , 'm:f:inf', 'resid'))

#plot sd of intercepts and interactions with infection together
ggplot(ppiMD, aes(x = level, y = sd))+
  geom_violin(draw_quantiles = c(0.09, 0.5, 0.91), alpha = 0.3, aes(fill = level), scale = 'width')+
  scale_fill_brewer(palette = 'Dark2')+
  labs(title = 'sd of model level estimates', y = 'sd', x = 'level')+
  theme(legend.position = 'none')



######################
#  female mortality  #
######################
#full model for FM
formFM <- brmsformula(deadeggs | trials(fdeadeggs)  ~ 1 + inf*fgeno*mgeno + (1|Day), family = binomial)
get_prior(formFM, data = data)
fiteFM <- brm(data = data, family = binomial,
              formula = formFM,
              prior = c(prior(normal(0,10), class = Intercept),
                        prior(normal(0,10), class = b),
                        prior(cauchy(0, 2), class = sd)),
              iter = 5000, warmup = 1000, chains = 2, cores = 2)
fiteFM <- add_criterion(fiteFM, 'waic')
fiteFM$fit

ndFM <- tibble(unique(data[c('inf','mgeno', 'fgeno')]))%>%mutate(fdeadeggs = 1) 
postFM <- fitted(fiteFM, newdata = ndFM, re_formula = NA, probs = c(0.09, 0.91))%>%
  as_tibble()%>%mutate(ndFM[c('mgeno', 'fgeno', 'inf')])

ggplot(postFM, aes(x = fgeno, y = Estimate))+
  geom_boxplot(data = data, aes(x = fgeno, y = FMobs, color = inf), position = position_dodge(0.8), width = 0.7)+
  geom_point(aes(class = inf), position = position_dodge(0.8))+
  geom_errorbar(aes(ymin=Q9, ymax=Q91, class = inf), width=.2, position=position_dodge(0.8))+
  labs(title = 'proportion non-male eggs failing to develop into adult', y = 'proportion failed eggs', x = 'female genotype')+
  facet_wrap('mgeno')+
  scale_color_discrete(name = expression(italic("Wolbachia")), labels = c('no', 'yes'))+
  theme(legend.position = c(0.8, 0.3))

fitFM_ <- fitted(fiteFM, newdata=ndFM, re_formula = NA, summary = 0)%>%as_tibble()
postFM2__ <- (fitFM_[seq(2,nrow(ndFM),2)]*(1-postMD_wide)-fitFM_[seq(1,nrow(ndFM),2)])/(1-fitFM_[seq(1,nrow(ndFM),2)])%>%as_tibble()
postFM2_wide <- postFM2__%>%mutate(across(everything(), function(x) {pmax(x, 0)}))
postFM2_ <- pivot_longer(postFM2_wide, cols = everything(), names_to = 'group', names_prefix = 'V', values_to = 'FM')%>%
  mutate(group = as.numeric(group))%>%mutate(mgeno = rep(ndFM$mgeno[seq(1,nrow(ndFM),2)], nrow(postFM2__)), fgeno = rep(ndFM$fgeno[seq(1,nrow(ndFM),2)], nrow(postFM2__)))
postFM2 <- group_by(postFM2_, group)%>%summarise(FMmean = mean(FM), FMq09 = quantile(FM, probs = 0.09), FMq91 = quantile(FM, probs = 0.91))%>%
  bind_cols(tibble(unique(data[c('mgeno', 'fgeno')])))

ggplot(postFM2_, aes(x = fgeno, y = FM))+
  geom_violin(draw_quantiles = c(0.09, 0.5, 0.91), alpha = 0.3, aes(fill = fgeno), scale = 'width')+
  facet_wrap('mgeno')+
  theme_light()+
  theme(legend.position = "none")

#plot FM of all crosses together, ordered according to CI
postFM2_$fxmgeno <- paste(postFM2_$fgeno, 'X', postFM2_$mgeno)
postFM2_$fxmgeno <- factor(postFM2_$fxmgeno, levels=unique(postCI2_$fxmgeno[order(postCI2$CImean)]), ordered=TRUE)
ggplot(postFM2_, aes(x = fxmgeno, y = FM))+
  geom_violin(draw_quantiles = c(0.09, 0.5, 0.91), alpha = 0.3, aes(fill = mgeno), scale = 'width')+
  coord_flip()+
  theme_light()+
  scale_fill_brewer(palette = 'Dark2')+
  labs(title = 'Female mortality', y = 'FM', x = 'female genotype')+
  theme(legend.position = "none")

#plot FM CI correlation
corrCIFM <- full_join(postCI2, postFM2)
ggplot(corrCIFM, aes(x = CImean, y = FMmean))+
  geom_point()+
  geom_errorbar(aes(ymin=FMq09, ymax=FMq91), width=.01)+
  geom_errorbarh(aes(xmin=CIq09, xmax=CIq91), height=.01)+
  theme_light()+
  lims(x = c(0,1), y = c(0,1))+
  scale_fill_brewer(palette = 'Dark2')+
  labs(title = 'CI-FM correlation', y = 'FM', x = 'CI')+
  theme(legend.position = "none")

#interspecific comparison (Bch - rest), intrasp as explanatory variable instead of fgeno
formFM2 <- brmsformula(deadeggs | trials(fdeadeggs)  ~ 1 + inf*intrasp*mgeno + (1|Day), family = binomial)
get_prior(formFM2, data = data_)
fiteFM2 <- brm(data = data_i, family = binomial,
              formula = formFM2,
              prior = c(prior(normal(0,10), class = Intercept),
                        prior(normal(0,10), class = b),
                        prior(cauchy(0, 2), class = sd)),
              iter = 5000, warmup = 1000, chains = 2, cores = 2)
fiteFM2 <- add_criterion(fiteFM2, 'waic')
fiteFM2$fit

ndFM2 <- tibble(unique(data_i[c('inf','mgeno', 'intrasp')]))%>%mutate(fdeadeggs = 1) 
postFM2i <- fitted(fiteFM2, newdata = ndFM2, re_formula = NA, probs = c(0.09, 0.91))%>%
  as_tibble()%>%mutate(ndFM2[c('mgeno', 'intrasp', 'inf')])

ggplot(postFM2i, aes(x = intrasp, y = Estimate))+
  geom_boxplot(data = data_, aes(x = intrasp, y = FMobs, color = inf), position = position_dodge(0.8), width = 0.7)+
  geom_point(aes(class = inf), position = position_dodge(0.8))+
  geom_errorbar(aes(ymin=Q9, ymax=Q91, class = inf), width=.2, position=position_dodge(0.8))+
  facet_wrap('mgeno')

postFM22_ <- fitted(fiteFM2, newdata=ndFM2, re_formula = NA, summary = 0)%>%as_tibble()
postFM22__ <- (postFM22_[seq(2,nrow(ndFM2),2)]*(1-postMD2_wide)-postFM22_[seq(1,nrow(ndFM2),2)])/(1-postFM22_[seq(1,nrow(ndFM2),2)])%>%as_tibble()
postFM22_wide <- postFM22__%>%mutate(across(everything(), function(x) {pmax(x, 0)}))
postFM22_ <- pivot_longer(postFM22_wide, cols = everything(), names_to = 'group', names_prefix = 'V', values_to = 'FM')%>%
  mutate(group = as.numeric(group))%>%mutate(mgeno = rep(ndFM2$mgeno[seq(1,nrow(ndFM2),2)], nrow(postFM22__)), intrasp = rep(ndFM2$intrasp[seq(1,nrow(ndFM2),2)], nrow(postFM22__)))
postFM22 <- group_by(postFM22_, group)%>%summarise(mean = mean(FM), q09 = quantile(FM, probs = 0.09), q91 = quantile(FM, probs = 0.91))

ggplot(postFM22_, aes(x = mgeno, y = FM, fill = instrasp))+
  geom_violin(draw_quantiles = c(0.09, 0.5, 0.91), alpha = 0.3, aes(fill = intrasp), scale = 'width')+
  theme_light()+
  ylim(0, 1)+
  labs(title = 'Female mortality', y = 'FM', x = 'female genotype', fill ='intrasp. cross')


#model comparison
#only male
formFMM <- brmsformula(deadeggs | trials(fdeadeggs)  ~ 1 + inf*mgeno + (1|Day), family = binomial)
get_prior(formFMM, data = data)
fiteFMM <- brm(data = data, family = binomial,
               formula = formFMM,
               prior = c(prior(normal(0,10), class = Intercept),
                         prior(normal(0,10), class = b),
                         prior(cauchy(0, 2), class = sd)),
               iter = 5000, warmup = 1000, chains = 2, cores = 2)
fiteFMM <- add_criterion(fiteFMM, 'waic')
fiteFMM$fit

ndFMM <- tibble(unique(data[c('inf', 'mgeno')]))%>%mutate(fdeadeggs = 1)#round(mean(data$fdeadeggs)))
postFMM <- fitted(fiteFMM, re_formula = deadeggs | trials(fdeadeggs)  ~ 1 + inf*mgeno, newdata = ndFMM, probs = c(0.09, 0.91))%>%
  as_tibble()%>%mutate(ndFMM[c('mgeno', 'inf')])

ggplot(postFMM, aes(x = mgeno, y = Estimate))+
  geom_boxplot(data = data, aes(x = mgeno, y = FMobs, color = inf), position = position_dodge(0.8), width = 0.7)+
  geom_point(aes(class = inf), position = position_dodge(0.8))+
  geom_errorbar(aes(ymin=Q9, ymax=Q91, class = inf), width=.2, position=position_dodge(0.8))+
  labs(title = 'proportion non-male eggs failing to develop into adult', y = 'proportion failed eggs', x = 'male genotype')+
  scale_color_discrete(name = expression(italic("Wolbachia")), labels = c('no', 'yes'))+
  theme(legend.position = 'bottom')


#only females
formFMF <- brmsformula(deadeggs | trials(fdeadeggs)  ~ 1 + inf*fgeno + (1|Day), family = binomial)
get_prior(formFMF, data = data)
fiteFMF <- brm(data = data, family = binomial,
               formula = formFMF,
               prior = c(prior(normal(0,10), class = Intercept),
                         prior(normal(0,10), class = b),
                         prior(cauchy(0, 2), class = sd)),
               iter = 5000, warmup = 1000, chains = 2, cores = 2)
fiteFMF <- add_criterion(fiteFMF, 'waic')
fiteFMF$fit
  
ndFMF <- tibble(unique(data[c('inf', 'fgeno')]))%>%mutate(fdeadeggs = 1)#round(mean(data$fdeadeggs)))
postFMF <- fitted(fiteFMF, re_formula = deadeggs | trials(fdeadeggs)  ~ 1 + inf*fgeno, newdata = ndFMF, probs = c(0.09, 0.91))%>%
  as_tibble()%>%mutate(ndFMF[c('fgeno', 'inf')])

ggplot(postFMF, aes(x = fgeno, y = Estimate))+
  geom_boxplot(data = data, aes(x = fgeno, y = FMobs, color = inf), position = position_dodge(0.8), width = 0.7)+
  geom_point(aes(class = inf), position = position_dodge(0.8))+
  geom_errorbar(aes(ymin=Q9, ymax=Q91, class = inf), width=.2, position=position_dodge(0.8))+
  labs(title = 'proportion non-male eggs failing to develop into adult', y = 'proportion failed eggs', x = 'female genotype')+
  scale_color_discrete(name = expression(italic("Wolbachia")), labels = c('no', 'yes'))+
  theme(legend.position = 'bottom')


#males and females, no interaction
formFMFM <- brmsformula(deadeggs | trials(fdeadeggs)  ~ 1 + inf*mgeno+ inf*fgeno + (1|Day), family = binomial)
get_prior(formFMFM, data = data)
fiteFMFM <- brm(data = data, family = binomial,
               formula = formFMFM,
               prior = c(prior(normal(0,10), class = Intercept),
                         prior(normal(0,10), class = b),
                         prior(cauchy(0, 2), class = sd)),
               iter = 5000, warmup = 1000, chains = 2, cores = 2)
fiteFMFM <- add_criterion(fiteFMFM, 'waic')
fiteFMFM$fit

ndFMFM <- tibble(unique(data[c('inf', 'fgeno', 'mgeno')]))%>%mutate(fdeadeggs = 1)#round(mean(data$fdeadeggs)))
postFMFM <- fitted(fiteFMFM, re_formula = deadeggs | trials(fdeadeggs)  ~ 1 + inf*fgeno + inf*mgeno, newdata = ndFMFM, probs = c(0.09, 0.91))%>%
  as_tibble()%>%mutate(ndFMFM[c('fgeno', 'mgeno', 'inf')])

ggplot(postFMFM, aes(x = fgeno, y = Estimate))+
  geom_boxplot(data = data, aes(x = fgeno, y = FMobs, color = inf), position = position_dodge(0.8), width = 0.7)+
  geom_point(aes(class = inf), position = position_dodge(0.8))+
  geom_errorbar(aes(ymin=Q9, ymax=Q91, class = inf), width=.2, position=position_dodge(0.8))+
  facet_wrap('mgeno')+
  labs(title = 'proportion non-male eggs failing to develop into adult', y = 'proportion failed eggs', x = 'female genotype')+
  scale_color_discrete(name = expression(italic("Wolbachia")), labels = c('no', 'yes'))+
  theme(legend.position = c(0.8, 0.3))

#compare all WAICs
FMcomp <- loo_compare(fiteFM, fiteFMM, fiteFMF, fiteFMFM, criterion = 'waic')
print(FMcomp, simplify = 0)
FMcomp[, 7:8] %>% 
  data.frame() %>% 
  rownames_to_column(var = "model_name") %>% 
  ggplot(aes(x    = model_name, y    = waic, ymin = waic - se_waic, ymax = waic + se_waic)) +
  geom_pointrange(shape = 21) +
  coord_flip() +
  scale_x_discrete(labels = c('full', 'F', 'M+F', 'M'))+
  labs(x = NULL, y = 'WAIC', title = "WAIC plot") +
  theme_classic()

#variance analysis FM
#full model with all fgeno, mgeno and interactions as variable effects
formFM_ <- brmsformula(deadeggs | trials(fdeadeggs)  ~ 1 + inf + (1+inf|fgeno*mgeno) + (1|Day), family = binomial)
get_prior(formFM_, data = data)
fiteFM_ <- brm(data = data, family = binomial,
               formula = formFM_,
               prior = c(prior(normal(0,5), class = Intercept),
                         prior(normal(0,5), class = b),
                         prior(cauchy(0, 1), class = sd),
                         prior(lkj(2), class = cor)),
               iter = 5000, warmup = 1000, chains = 2, cores = 2)
fiteFM_ <- add_criterion(fiteFM_, 'waic')
fiteFM_$fit

ndFM_ <- tibble(unique(data[c('inf','mgeno', 'fgeno', 'Day')]))%>%mutate(fdeadeggs = 1) 
postFM_ <- fitted(fiteFM_, re_formula = deadeggs | trials(fdeadeggs) ~ 1 + inf + (1+inf|fgeno*mgeno), newdata = ndFM_, probs = c(0.09, 0.91))%>%
  as_tibble()%>%mutate(ndFM_[c('mgeno', 'fgeno', 'inf')])

ggplot(postFM_, aes(x = fgeno, y = Estimate))+
  geom_boxplot(data = data, aes(x = fgeno, y = FMobs, color = inf), position = position_dodge(0.8), width = 0.7)+
  geom_point(aes(class = inf), position = position_dodge(0.8))+
  geom_errorbar(aes(ymin=Q9, ymax=Q91, class = inf), width=.2, position=position_dodge(0.8))+
  labs(title = 'proportion non-male eggs failing to develop into adult', y = 'proportion failed eggs', x = 'female genotype')+
  scale_color_discrete(name = expression(italic("Wolbachia")), labels = c('no', 'yes'))+
  facet_wrap('mgeno')+
  theme(legend.position = c(0.8, 0.3))

#expected values on log-odds scale
fittedFM1 <- fitted(fiteFM_, summary = 0, scale = 'linear')
#calculate residuals on log-odds scale
observed_logit <- logit(data_$FMobs, adjust = 0.0025)                   
residFM <- sweep(as.matrix(fittedFM1), 2, observed_logit, '-')%>%as_tibble()%>%rowwise()%>%transmute(residuals = sd(c_across(everything())))

#calculate finite-population standard deviation of interactions with infection
pi_FM <- as_tibble(fiteFM_)%>%rowwise()%>%
  transmute(sd_mgenoi = sd(c_across(starts_with('r_mgeno[') & ends_with('inf1]'))),
            sd_fgenoi = sd(c_across(starts_with('r_fgeno[') & ends_with('inf1]'))),
            sd_mfi = sd(c_across(starts_with('r_fgeno:mgeno[') & ends_with('inf1]'))))%>%
  bind_cols(residFM$residuals)%>%rename(resid = ...4)%>%
  mutate(dm_fi = sd_mgenoi-sd_fgenoi, dm_mfi = sd_mgenoi-sd_mfi, df_mfi = sd_mfi-sd_fgenoi)

#calculate finite-population standard deviation of intercepts (no infection effect)
p_FM <- as_tibble(fiteFM_)%>%rowwise()%>%
  transmute(sd_mgeno = sd(c_across(starts_with('r_mgeno[') & ends_with('Intercept]'))),
            sd_fgeno = sd(c_across(starts_with('r_fgeno[') & ends_with('Intercept]'))),
            sd_mf = sd(c_across(starts_with('r_fgeno:mgeno[') & ends_with('Intercept]'))))%>%
  mutate(dm_f = sd_mgeno-sd_fgeno, dm_mf = sd_mgeno-sd_mf, df_mf = sd_mf-sd_fgeno)
ppiFM <- cbind(p_FM[c('sd_mgeno', 'sd_fgeno', 'sd_mf')], pi_FM[c('sd_mgenoi', 'sd_fgenoi', 'sd_mfi', 'resid')])%>%
  rename(fgeno = sd_fgeno, mgeno = sd_mgeno, 'fgeno:inf' = sd_fgenoi, 'mgeno:inf' = sd_mgenoi, 'm:f' = sd_mf, 'm:f:inf' = sd_mfi)%>%
  pivot_longer(.[c('mgeno:inf', 'fgeno:inf', 'm:f:inf', 'resid', 'mgeno', 'fgeno', 'm:f')], cols = everything(), names_to = 'level', values_to = 'sd')
ppiFM$level <- factor(ppiFM$level, levels = c('fgeno', 'fgeno:inf','mgeno','mgeno:inf',  'm:f' , 'm:f:inf', 'resid'))

#plot sd of all effects and residuals
ggplot(ppiFM, aes(x = level, y = sd))+
  geom_violin(draw_quantiles = c(0.09, 0.5, 0.91), alpha = 0.3, aes(fill = level), scale = 'width')+
  scale_fill_brewer(palette = 'Dark2')+
  labs(title = 'sd of model level estimates', y = 'sd', x = 'level')+
  theme(legend.position = 'none')

#FM and MD in one plot
postFMxMD <- cbind(postFM2_, MD = postMD_$MD)%>%pivot_longer(cols = c('MD', 'FM'), names_to = 'metric', values_to = 'estimate')

(fig2 <- ggplot(postFMxMD, aes(x = fxmgeno, y = estimate, scale = metric))+
  geom_violin(draw_quantiles = c(0.09, 0.5, 0.91), alpha = 0.3, aes(fill = metric), scale = 'width', position = position_dodge(0))+
  coord_flip()+
  labs(title = 'Cytoplasmic Incompatibility mechanisms', x = 'cross', y = 'mechanism strength')+
  theme_light())

#FM and MD in one plor for Bch comparison
postFMxMDii <- cbind(postFM22_[!postFM22_$intrasp,], MD = postMD2_[!postMD2_$intrasp,]$MD)%>%pivot_longer(cols = c('MD', 'FM'), names_to = 'metric', values_to = 'estimate')
postFMxMDii$mgeno <- factor(postFMxMDii$mgeno, levels=unique(postCI3_$mgeno[order(postCI3[postCI3$intrasp,]$mean)]), ordered=TRUE)
(fig2b_ <- ggplot(postFMxMDii, aes(x = mgeno, y = estimate, scale = metric))+
    geom_violin(draw_quantiles = c(0.09, 0.5, 0.91), alpha = 0.3, aes(fill = metric), scale = 'width', position = position_dodge(0))+
    coord_flip()+
    ylim(0,1)+
    labs(title = 'Cytoplasmic Incompatibility mechanisms', x = 'cross', y = 'mechanism strength')+
    theme_light()+
    theme(legend.position = 'none'))

#define fixed plot sizes (of panel)
fig2aset <- set_panel_size(fig2,
                          width  = unit(8, "cm"),
                          height = unit(20, "cm"))
fig2bset <- set_panel_size(fig2b_,
                          width  = unit(8, "cm"),
                          height = unit(4, "cm"))
grid.newpage()
grid.draw(fig2aset) #rescale png save to height 825, jeeping aspect ration
grid.newpage()
grid.draw(fig2bset)



####################
# F1 & F2 analysis #
####################
F1 <- read.table('Wybouwetal_CI_F1', header = 1)
F2 <- read.table('Wybouwetal_CI_F2', header = 1)
colnames(F1)[3]<-'inf'
colnames(F2)[3]<-'inf'
F1$Fprop <- F1$adult_females/F1$eggs
F2$Fprop <- F2$adult_females/F2$eggs
F1$Mprop <- F1$adult_males/F1$eggs
F2$Mprop <- F2$adult_males/F2$eggs

str(F1)
F1$deadeggs <- F1$eggs - F1$adult_females - F1$adult_males
F1$fdeadeggs <- F1$eggs - F1$adult_males
F1$FMprop <- F1$deadeggs/F1$fdeadeggs
F2$deadeggs <- F2$eggs - F2$adult_females - F2$adult_males
F2$fdeadeggs <- F2$eggs - F2$adult_males
F2$FMprop <- F2$deadeggs/F2$fdeadeggs

levels(F1$inf) <- list("0" = "uninfected", "1" = "infected")
levels(F2$inf) <- list("0" = "uninfected", "1" = "infected")
ggplot(F1, aes(x = Cross, y = Fprop))+
  geom_boxplot(aes(color = inf), position = position_dodge(0.8), width = 0.7)
str(F2)
ggplot(F2, aes(x = Cross, y = Fprop))+
  geom_boxplot(aes(color = inf), position = position_dodge(0.8), width = 0.7)


#F1 MD
forMDF1 <- brmsformula(adult_males | trials(eggs)  ~ 1 + inf*Cross, family = binomial)
get_prior(forMDF1, data = F1)
fitMDF1 <- brm(data = F1, family = binomial,
             formula = forMDF1,
             prior = c(prior(normal(0,10), class = Intercept),
                       prior(normal(0,10), class = b)),
             iter = 5000, warmup = 1000, chains = 2, cores = 2)
fitMDF1 <- add_criterion(fitMDF1, 'waic')
fitMDF1$fit

ndF1 <- tibble(unique(F1[c('inf','Cross')]))%>%mutate(eggs = 1) 
MDpostF1 <- fitted(fitMDF1, newdata = ndF1, re_formula = NA, probs = c(0.09, 0.91))%>%
  as_tibble()%>%mutate(ndF1[c('Cross', 'inf')])

ggplot(MDpostF1, aes(x = Cross, y = Estimate))+
  geom_boxplot(data = F1, aes(x = Cross, y = Mprop, color = inf), position = position_dodge(0.8), width = 0.7)+
  geom_point(aes(class = inf), position = position_dodge(0.8))+
  geom_errorbar(aes(ymin=Q9, ymax=Q91, class = inf), width=.2, position=position_dodge(0.8))+
  labs(title = 'proportion eggs succesfully developing into adult males', y = 'proportion male', x = 'cross')+
  scale_color_discrete(name = expression(italic("Wolbachia")), labels = c('no', 'yes'))

MDpostF1_ <- fitted(fitMDF1, newdata=ndF1, re_formula = NA, summary = 0)%>%as_tibble()
postMDF1__ <- (MDpostF1_[seq(2,nrow(ndF1),2)]-MDpostF1_[seq(1,nrow(ndF1),2)])/(1-MDpostF1_[seq(1,nrow(ndF1),2)])%>%as_tibble()
postMDF1_wide <- postMDF1__%>%mutate(across(everything(), function(x) {pmax(x, 0)}))
postMDF1_ <- pivot_longer(postMDF1_wide, cols = everything(), names_to = 'group', names_prefix = 'V', values_to = 'MD')%>%
  mutate(group = as.numeric(group))%>%mutate(Cross = rep(ndF1$Cross[seq(1,nrow(ndF1),2)], nrow(postMDF1__)))
postMDF1 <- group_by(postMDF1_, group)%>%summarise(mean = mean(MD), q09 = quantile(MD, probs = 0.09), q91 = quantile(MD, probs = 0.91))


ggplot(postMDF1_, aes(x = Cross, y = MD))+
  theme_light()+
  geom_violin(draw_quantiles = c(0.09, 0.5, 0.91), alpha = 0.5, aes(fill = Cross), scale = 'width')+
  scale_fill_brewer(palette = 'Dark2')+
  labs(title = 'Male development', y = 'MD', x = 'Cross')+
  theme(legend.position = "none")

ff <- postMD_[postMD_$fxmgeno %in% c('LonX X Scp', 'Scp X Scp'),]%>%transmute('group'=group, 'MD' = MD, 'Cross' = fxmgeno)
ff$Cross <- factor(ff$Cross, ordered = 0)
postMDF0F1 <- bind_rows(postMDF1_, ff)
postMDF0F1$Cross <- factor(postMDF0F1$Cross, levels = c('LonX X Scp', 'LonX|Scp_x_Scp', 'Scp|LonX_x_Scp', 'Scp X Scp'))
ggplot(postMDF0F1, aes(x = Cross, y = MD))+
  theme_light()+
  geom_violin(draw_quantiles = c(0.09, 0.5, 0.91), alpha = 0.3, fill = "#00BFC4", scale = 'width')+
  coord_flip()+
  labs(title = 'Male development in F1s and parental crosses', y = 'MD', x = 'Cross')+
  theme(legend.position = "none")

#F1 FM
forFMF1 <- brmsformula(deadeggs | trials(fdeadeggs)  ~ 1 + inf*Cross, family = binomial)
get_prior(forFMF1, data = FMF1)
fitFMF1 <- brm(data = F1, family = binomial,
             formula = forFMF1,
             prior = c(prior(normal(0,10), class = Intercept),
                       prior(normal(0,10), class = b)),
             iter = 5000, warmup = 1000, chains = 2, cores = 2)
fitFMF1 <- add_criterion(fitFMF1, 'waic')
fitFMF1$fit

ndFMF1 <- tibble(unique(F1[c('inf','Cross')]))%>%mutate(fdeadeggs = 1) 
FMpostF1 <- fitted(fitFMF1, newdata = ndFMF1, re_formula = NA, probs = c(0.09, 0.91))%>%
  as_tibble()%>%mutate(ndFMF1[c('Cross', 'inf')])

ggplot(FMpostF1, aes(x = Cross, y = Estimate))+
  geom_boxplot(data = F1, aes(x = Cross, y = FMprop, color = inf), position = position_dodge(0.8), width = 0.7)+
  geom_point(aes(class = inf), position = position_dodge(0.8))+
  geom_errorbar(aes(ymin=Q9, ymax=Q91, class = inf), width=.2, position=position_dodge(0.8))+
  labs(title = 'proportion non-male eggs failing to develop into adult', y = 'proportion failed eggs', x = 'cross')+
  scale_color_discrete(name = expression(italic("Wolbachia")), labels = c('no', 'yes'))

FMpostF1_ <- fitted(fitFMF1, newdata=ndFMF1, re_formula = NA, summary = 0)%>%as_tibble()
postFMF1__ <- (FMpostF1_[seq(2,nrow(ndF1),2)]*(1-postMDF1_wide)-FMpostF1_[seq(1,nrow(ndFMF1),2)])/(1-FMpostF1_[seq(1,nrow(ndFMF1),2)])%>%as_tibble()
postFMF1_wide <- postFMF1__%>%mutate(across(everything(), function(x) {pmax(x, 0)}))
postFMF1_ <- pivot_longer(postFMF1_wide, cols = everything(), names_to = 'group', names_prefix = 'V', values_to = 'FM')%>%
  mutate(group = as.numeric(group))%>%mutate(Cross = rep(ndF1$Cross[seq(1,nrow(ndF1),2)], nrow(postFMF1__)))
postFMF1 <- group_by(postFMF1_, group)%>%summarise(mean = mean(FM), q09 = quantile(FM, probs = 0.09), q91 = quantile(FM, probs = 0.91))


ggplot(postFMF1_, aes(x = Cross, y = FM))+
  theme_light()+
  geom_violin(draw_quantiles = c(0.09, 0.5, 0.91), alpha = 0.5, aes(fill = Cross), scale = 'width')+
  scale_fill_brewer(palette = 'Dark2')+
  labs(title = 'Male development', y = 'FM', x = 'Cross')+
  theme(legend.position = "none")

ffFM <- postFM2_[postFM2_$fxmgeno %in% c('LonX X Scp', 'Scp X Scp'),]%>%transmute('group'=group, 'FM' = FM, 'Cross' = fxmgeno)
ffFM$Cross <- factor(ffFM$Cross, ordered = 0)
postFMF0F1 <- bind_rows(postFMF1_, ffFM)
postFMF0F1$Cross <- factor(postFMF0F1$Cross, levels = c('LonX X Scp', 'LonX|Scp_x_Scp', 'Scp|LonX_x_Scp', 'Scp X Scp'))
ggplot(postFMF0F1, aes(x = Cross, y = FM))+
  theme_light()+
  geom_violin(draw_quantiles = c(0.09, 0.5, 0.91), alpha = 0.3, fill = "#00BFC4", scale = 'width')+
  coord_flip()+
  labs(title = 'Male development in F1s and parental crosses', y = 'FM', x = 'Cross')+
  theme(legend.position = "none")

postF1FMxMD <- cbind(postFMF0F1, MD = postMDF0F1$MD)%>%pivot_longer(cols = c('MD', 'FM'), names_to = 'metric', values_to = 'estimate')
(fig3a <- ggplot(postF1FMxMD, aes(x = Cross, y = estimate, scale = metric))+
    geom_violin(draw_quantiles = c(0.09, 0.5, 0.91), alpha = 0.3, aes(fill = metric), scale = 'width', position = position_dodge(0))+
    coord_flip()+
    ylim(0,1)+
    labs(title = 'Cytoplasmic Incompatibility mechanisms', x = 'cross', y = 'mechanism strength')+
    theme_light())

fig3aset <- set_panel_size(fig3a,
                           width  = unit(8, "cm"),
                           height = unit(3.2, "cm"))
grid.newpage()
grid.draw(fig3aset) #rescale png save to height 825, jeeping aspect ration


#F1 CI
forF1CI <- brmsformula(adult_females | trials(eggs)  ~ 1 + inf*Cross, family = binomial)
get_prior(forF1CI, data = F1)
fitF1CI <- brm(data = F1, family = binomial,
             formula = forF1CI,
             prior = c(prior(normal(0,10), class = Intercept),
                       prior(normal(0,10), class = b)),
             iter = 5000, warmup = 1000, chains = 2, cores = 2)
fitF1CI <- add_criterion(fitF1CI, 'waic')
fitF1CI$fit

ndF1CI <- tibble(unique(F1[c('inf','Cross')]))%>%mutate(eggs = 1) 
postF1CI <- fitted(fitF1CI, newdata = ndF1CI, re_formula = NA, probs = c(0.09, 0.91))%>%
  as_tibble()%>%mutate(ndF1CI[c('Cross', 'inf')])

ggplot(postF1CI, aes(x = Cross, y = Estimate))+
  geom_boxplot(data = F1, aes(x = Cross, y = Fprop, color = inf), position = position_dodge(0.8), width = 0.7)+
  geom_point(aes(class = inf), position = position_dodge(0.8))+
  geom_errorbar(aes(ymin=Q9, ymax=Q91, class = inf), width=.2, position=position_dodge(0.8))+
  labs(title = 'proportion eggs succesfully developing into adult females', y = 'proportion female', x = 'cross')+
  scale_color_discrete(name = expression(italic("Wolbachia")), labels = c('no', 'yes'))

postF1CI_ <- fitted(fitF1CI, newdata=ndF1CI, re_formula = NA, summary = 0)%>%as_tibble()
postCIF1__ <- 1-(postF1CI_[seq(2,nrow(ndF1CI),2)]/postF1CI_[seq(1,nrow(ndF1CI),2)])%>%as_tibble()
postCIF1_wide <- postCIF1__%>%mutate(across(everything(), function(x) {pmax(x, 0)}))
postCIF1_ <- pivot_longer(postCIF1_wide, cols = everything(), names_to = 'group', names_prefix = 'V', values_to = 'CI')%>%
  mutate(group = as.numeric(group))%>%mutate(Cross = rep(ndF1CI$Cross[seq(1,nrow(ndF1CI),2)], nrow(postCIF1__)))
postCIF1 <- group_by(postCIF1_, group)%>%summarise(mean = mean(CI), q09 = quantile(CI, probs = 0.09), q91 = quantile(CI, probs = 0.91))


ggplot(postCIF1_, aes(x = Cross, y = CI))+
  theme_light()+
  geom_violin(draw_quantiles = c(0.09, 0.5, 0.91), alpha = 0.5, aes(fill = Cross), scale = 'width')+
  scale_fill_brewer(palette = 'Dark2')+
  labs(title = 'Cytoplasmatic Incompatibility', y = 'CI', x = 'Cross')+
  theme(legend.position = "none")

ffCI <- postCI2_[postCI2_$fxmgeno %in% c('LonX X Scp', 'Scp X Scp'),]%>%transmute('group'=group, 'CI' = CI, 'Cross' = fxmgeno)
ffCI$Cross <- factor(ffCI$Cross, ordered = 0)
postMDF0F1CI <- bind_rows(postCIF1_, ffCI)
postMDF0F1CI$Cross <- factor(postMDF0F1CI$Cross, levels = c('LonX X Scp', 'LonX|Scp_x_Scp', 'Scp|LonX_x_Scp', 'Scp X Scp'))
ggplot(postMDF0F1CI, aes(x = Cross, y = CI))+
  theme_light()+
  geom_violin(draw_quantiles = c(0.09, 0.5, 0.91), alpha = 0.3, scale = 'width')+
  coord_flip()+
  ylim(0,1)+
  labs(title = 'Cytoplasmatic incompatibility in F1s and parental crosses', y = 'CI', x = 'Cross')+
  theme(legend.position = "none")


#F2
#MD in F2
ggplot(F2[F2$Cross == '(LonX|Scp)|LonX_x_Scp',], aes(x = Mprop, fill = inf))+
  geom_histogram(bins = 20, alpha = 0.3, position = 'dodge')+
  geom_density(alpha=.2)+
  labs(title = 'male development in F2 back-crosses', y = 'count', x = 'proportion male')+
  scale_fill_brewer(palette = 'Dark2', name = "", labels = c('cured', expression(paste(italic("Wolbachia"),"-infected"))))+
  xlim(0, 1)+
  theme(legend.position = 'top')

forF2 <- brmsformula(adult_males | trials(eggs)  ~ 1 + inf*Cross, family = binomial)
get_prior(forF2, data = F2)
fitF2 <- brm(data = F2, family = binomial,
             formula = forF2,
             prior = c(prior(normal(0,10), class = Intercept),
                       prior(normal(0,10), class = b)),
             iter = 5000, warmup = 1000, chains = 2, cores = 2)
fitF2 <- add_criterion(fitF2, 'waic')
fitF2$fit

ndF2 <- tibble(unique(F2[c('inf','Cross')]))%>%mutate(eggs = 1) 
postF2 <- fitted(fitF2, newdata = ndF2, re_formula = NA, probs = c(0.09, 0.91))%>%
  as_tibble()%>%mutate(ndF2[c('Cross', 'inf')])

ggplot(postF2, aes(x = Cross, y = Estimate))+
  geom_boxplot(data = F2, aes(x = Cross, y = Mprop, color = inf), position = position_dodge(0.8), width = 0.7)+
  geom_point(aes(class = inf), position = position_dodge(0.8))+
  geom_errorbar(aes(ymin=Q9, ymax=Q91, class = inf), width=.2, position=position_dodge(0.8))+
  labs(title = 'proportion eggs succesfully developing into adult males', y = 'proportion male', x = 'female genotype')+
  scale_color_discrete(name = expression(italic("Wolbachia")), labels = c('no', 'yes'))

postF2_ <- fitted(fitF2, newdata=ndF2, re_formula = NA, summary = 0)%>%as_tibble()
postMDF2__ <- (postF2_[seq(2,nrow(ndF2),2)]-postF2_[seq(1,nrow(ndF2),2)])/(1-postF2_[seq(1,nrow(ndF2),2)])%>%as_tibble()
postMDF2_wide <- postMDF2__%>%mutate(across(everything(), function(x) {pmax(x, 0)}))
postMDF2_ <- pivot_longer(postMDF2_wide, cols = everything(), names_to = 'group', names_prefix = 'V', values_to = 'MD')%>%
  mutate(group = as.numeric(group))%>%mutate(Cross = rep(ndF2$Cross[seq(1,nrow(ndF2),2)], nrow(postMDF2__)))
postMDF2 <- group_by(postMDF2_, group)%>%summarise(mean = mean(MD), q09 = quantile(MD, probs = 0.09), q91 = quantile(MD, probs = 0.91))


ggplot(postMDF2_, aes(x = Cross, y = MD))+
  theme_light()+
  geom_violin(draw_quantiles = c(0.09, 0.5, 0.91), alpha = 0.5, aes(fill = Cross), scale = 'width')+
  scale_fill_brewer(palette = 'Dark2')+
  labs(title = 'Male development in F2 back-crosses', y = 'MD', x = 'Cross')+
  theme(legend.position = "none")


#mixture model for MD in F2 without the other cross
mix <- mixture(binomial, binomial)
mixF2__ <- brmsformula(adult_males | trials(eggs)  ~ 1 + inf, family = mix)
get_prior(mixF2, data = F2[F2$Cross == '(LonX|Scp)|LonX_x_Scp',])
fitmixF2__ <- brm(data = F2[F2$Cross == '(LonX|Scp)|LonX_x_Scp',], family = binomial,
                formula = mixF2__,
                prior = c(prior(normal(0,5), class = Intercept, dpar = mu1),
                          prior(normal(0,5), class = Intercept, dpar = mu2),
                          prior(normal(0,5), class = b, dpar = mu1),
                          prior(normal(0,5), class = b, dpar = mu2)),
                iter = 5000, warmup = 1000, chains = 2, cores = 2)
fitmixF2__ <- add_criterion(fitmixF2__, 'waic')
fitmixF2__$fit


ndF2mix__ <- tibble(unique(F2[c('inf')]))%>%mutate(eggs = 1) 
postF2mix <- fitted(fitmixF2__, newdata = ndF2mix__, re_formula = NA, probs = c(0.09, 0.91))%>%
  as_tibble()%>%mutate(ndF2mix__[c('inf')])

ggplot(postF2mix, aes(x = 1, y = Estimate))+
  geom_boxplot(data = F2[F2$Cross == '(LonX|Scp)|LonX_x_Scp',], aes(x = Cross, y = Mprop, color = inf), position = position_dodge(0.8), width = 0.7)+
  geom_point(aes(class = inf), position = position_dodge(0.8))+
  geom_errorbar(aes(ymin=Q9, ymax=Q91, class = inf), width=.2, position=position_dodge(0.8))+
  labs(title = 'proportion eggs succesfully developing into adult males', y = 'proportion male', x = 'female genotype')+
  scale_color_discrete(name = expression(italic("Wolbachia")), labels = c('no', 'yes'))

postF2mixd_ <- posterior_samples(fitmixF2__)%>%
  mutate(distnumber = rbinom(nrow(.),1,prob = .$theta2))%>%
  mutate('lo_uninfected' = (b_mu1_Intercept*!distnumber)+(b_mu2_Intercept*distnumber),
         'lo_infected' = (b_mu1_Intercept*!distnumber)+(b_mu2_Intercept*distnumber) + (b_mu1_inf1*!distnumber) + (b_mu2_inf1*distnumber))%>%
  mutate('uninfected' = inv_logit_scaled(lo_uninfected),
         'infected' = inv_logit_scaled(lo_infected))%>%
  select(distnumber, uninfected, infected)%>%
  pivot_longer(c(uninfected, infected), names_to = 'inf', values_to = 'Mprop')%>%
  mutate(inf = as.factor(as.numeric(inf == 'infected')), distnumber = as.factor(distnumber))

ggplot(posterior_samples(fitmixF2__), aes(x = 'theta: mixing coefficient', y = theta2))+
  geom_violin(draw_quantiles = c(0.09, 0.5, 0.91))+
  labs(title = 'proportion of high MD phenotypes', y = 'weight of highest estimated distribution', x = '')+
  ylim(c(0,1))+
  geom_hline(yintercept = 0.5, lty = 2, color = 'blue', lwd = 1)


simeggs <- 100
postF2mixsim <- predict(fitmixF2__, newdata = ndF2mix__%>%mutate(eggs = simeggs), re_formula = NA, probs = c(0.09, 0.91), summary = 0)%>%
  as_tibble()%>%
  rename('0' = V1, '1' = V2)%>%
  pivot_longer(everything(), names_to = 'inf', values_to = 'Mprop')%>%
  mutate(Mprop = Mprop/simeggs, inf = factor(inf))

(fig3d <- ggplot(F2[F2$Cross == '(LonX|Scp)|LonX_x_Scp',], aes(x = Mprop, fill = inf))+
  geom_histogram(bins = 20, alpha = 0.3, position = 'dodge')+
  geom_density(data = postF2mixsim, aes(x = Mprop, fill = inf),alpha=.2, bw = 0.05)+
  labs(title = 'male development in F2 back-crosses', y = 'count', x = 'proportion male')+
  scale_fill_brewer(palette = 'Dark2', name = "", labels = c('cured', expression(paste(italic("Wolbachia"),"-infected"))))+
  xlim(0, 1)+
  theme(legend.position = 'top'))

fig3dset <- set_panel_size(fig3d,
                           width  = unit(8, "cm"),
                           height = unit(4, "cm"))
grid.newpage()
grid.draw(fig3dset) #rescale png save to height 825, jeeping aspect ration


#Fm in F2
forF2FM <- brmsformula(deadeggs | trials(fdeadeggs)  ~ 1 + inf*Cross, family = binomial)
get_prior(forF2FM, data = F2)
fitF2FM <- brm(data = F2, family = binomial,
             formula = forF2FM,
             prior = c(prior(normal(0,10), class = Intercept),
                       prior(normal(0,10), class = b)),
             iter = 5000, warmup = 1000, chains = 2, cores = 2)
fitF2FM <- add_criterion(fitF2FM, 'waic')
fitF2FM$fit

ndF2FM <- tibble(unique(F2[c('inf','Cross')]))%>%mutate(fdeadeggs = 1) 
postF2FM <- fitted(fitF2FM, newdata = ndF2FM, re_formula = NA, probs = c(0.09, 0.91))%>%
  as_tibble()%>%mutate(ndF2FM[c('Cross', 'inf')])

ggplot(postF2FM, aes(x = Cross, y = Estimate))+
  geom_boxplot(data = F2, aes(x = Cross, y = FMprop, color = inf), position = position_dodge(0.8), width = 0.7)+
  geom_point(aes(class = inf), position = position_dodge(0.8))+
  geom_errorbar(aes(ymin=Q9, ymax=Q91, class = inf), width=.2, position=position_dodge(0.8))+
  labs(title = 'proportion non-male eggs failing to develop into adult', y = 'proportion failed eggs', x = 'cross')+
  scale_color_discrete(name = expression(italic("Wolbachia")), labels = c('no', 'yes'))

FMpostF2_ <- fitted(fitF2FM, newdata=ndF2FM, re_formula = NA, summary = 0)%>%as_tibble()
postFMF2__ <- (FMpostF2_[seq(2,nrow(ndF2FM),2)]*(1-postMDF2_wide)-FMpostF2_[seq(1,nrow(ndF2FM),2)])/(1-FMpostF2_[seq(1,nrow(ndF2FM),2)])%>%as_tibble()
postFMF2_wide <- postFMF2__%>%mutate(across(everything(), function(x) {pmax(x, 0)}))
postFMF2_ <- pivot_longer(postFMF2_wide, cols = everything(), names_to = 'group', names_prefix = 'V', values_to = 'FM')%>%
  mutate(group = as.numeric(group))%>%mutate(Cross = rep(ndF2FM$Cross[seq(1,nrow(ndF2FM),2)], nrow(postFMF2__)))
postFMF2 <- group_by(postFMF2_, group)%>%summarise(mean = mean(FM), q09 = quantile(FM, probs = 0.09), q91 = quantile(FM, probs = 0.91))

ggplot(postFMF2_, aes(x = Cross, y = FM))+
  theme_light()+
  geom_violin(draw_quantiles = c(0.09, 0.5, 0.91), alpha = 0.5, aes(fill = Cross), scale = 'width')+
  scale_fill_brewer(palette = 'Dark2')+
  labs(title = 'Male development', y = 'FM', x = 'Cross')+
  theme(legend.position = "none")


postF2FMxMD <- cbind(postFMF2_, MD = postMDF2_$MD)%>%pivot_longer(cols = c('MD', 'FM'), names_to = 'metric', values_to = 'estimate')
(figsF2 <- ggplot(postF2FMxMD, aes(x = Cross, y = estimate, scale = metric))+
    geom_violin(draw_quantiles = c(0.09, 0.5, 0.91), alpha = 0.3, aes(fill = metric), scale = 'width', position = position_dodge(0))+
    coord_flip()+
    ylim(0,1)+
    labs(title = 'Cytoplasmic Incompatibility mechanisms', x = 'cross', y = 'mechanism strength')+
    theme_light())

figsF2set <- set_panel_size(figsF2,
                           width  = unit(8, "cm"),
                           height = unit(1.6, "cm"))
grid.newpage()
grid.draw(figsF2set) #rescale png save to height 825, jeeping aspect ration


#F2 CI
forF2CI <- brmsformula(adult_females | trials(eggs)  ~ 1 + inf*Cross, family = binomial)
get_prior(forF2CI, data = F2)
fitF2CI <- brm(data = F2, family = binomial,
               formula = forF2CI,
               prior = c(prior(normal(0,10), class = Intercept),
                         prior(normal(0,10), class = b)),
               iter = 5000, warmup = 1000, chains = 2, cores = 2)
fitF2CI <- add_criterion(fitF2CI, 'waic')
fitF2CI$fit

ndF2CI <- tibble(unique(F2[c('inf','Cross')]))%>%mutate(eggs = 1) 
postF2CI <- fitted(fitF2CI, newdata = ndF2CI, re_formula = NA, probs = c(0.09, 0.91))%>%
  as_tibble()%>%mutate(ndF2CI[c('Cross', 'inf')])

ggplot(postF2CI, aes(x = Cross, y = Estimate))+
  geom_boxplot(data = F2, aes(x = Cross, y = Fprop, color = inf), position = position_dodge(0.8), width = 0.7)+
  geom_point(aes(class = inf), position = position_dodge(0.8))+
  geom_errorbar(aes(ymin=Q9, ymax=Q91, class = inf), width=.2, position=position_dodge(0.8))+
  labs(title = 'proportion eggs succesfully developing into adult females', y = 'proportion female', x = 'cross')+
  scale_color_discrete(name = expression(italic("Wolbachia")), labels = c('no', 'yes'))

postF2CI_ <- fitted(fitF2CI, newdata=ndF2CI, re_formula = NA, summary = 0)%>%as_tibble()
postCIF2__ <- 1-(postF2CI_[seq(2,nrow(ndF2CI),2)]/postF2CI_[seq(1,nrow(ndF2CI),2)])%>%as_tibble()
postCIF2_wide <- postCIF2__%>%mutate(across(everything(), function(x) {pmax(x, 0)}))
postCIF2_ <- pivot_longer(postCIF2_wide, cols = everything(), names_to = 'group', names_prefix = 'V', values_to = 'CI')%>%
  mutate(group = as.numeric(group))%>%mutate(Cross = rep(ndF2CI$Cross[seq(1,nrow(ndF2CI),2)], nrow(postCIF2__)))
postCIF2 <- group_by(postCIF2_, group)%>%summarise(mean = mean(CI), q09 = quantile(CI, probs = 0.09), q91 = quantile(CI, probs = 0.91))


ggplot(postCIF2_, aes(x = Cross, y = CI))+
  theme_light()+
  geom_violin(draw_quantiles = c(0.09, 0.5, 0.91), alpha = 0.5, aes(fill = Cross), scale = 'width')+
  scale_fill_brewer(palette = 'Dark2')+
  labs(title = 'Male development', y = 'MD', x = 'Cross')+
  theme(legend.position = "none")

ffCI <- postCI2_[postCI2_$fxmgeno %in% c('LonX X Scp', 'Scp X Scp'),]%>%transmute('group'=group, 'CI' = CI, 'Cross' = fxmgeno)
ffCI$Cross <- factor(ffCI$Cross, ordered = 0)
postMDF0F2CI <- bind_rows(postCIF2_[postCIF2_$Cross == '(LonX|Scp)|LonX_x_Scp',], ffCI)
postMDF0F2CI$Cross <- factor(postMDF0F2CI$Cross, levels = c('LonX X Scp', '(LonX|Scp)|LonX_x_Scp', 'Scp X Scp'))
ggplot(postMDF0F2CI, aes(x = Cross, y = CI))+
  theme_light()+
  geom_violin(draw_quantiles = c(0.09, 0.5, 0.91), alpha = 0.3, scale = 'width')+
  coord_flip()+
  ylim(0,1)+
  labs(title = 'Cytoplasmatic incompatibility in F2s and parental crosses', y = 'CI', x = 'Cross')+
  theme(legend.position = "none")

#only rescue
#load same data, but filter to data of rescue crosses (Type = [R, Rc])
datares <- data__[data__$Type %in% c('R', 'Rc'),]
datares <- droplevels(datares)

formres <- brmsformula(adult_females | trials(eggs)  ~ 1 + inf*Cross, family = binomial)
get_prior(formres, data = datares)
fitres <- brm(data = datares, family = binomial,
              formula = formres,
              prior = c(prior(normal(0,10), class = Intercept),
                        prior(normal(0,10), class = b)),
              iter = 5000, warmup = 2000, chains = 2, cores = 2)
fitres <- add_criterion(fitres, 'waic')
fitres$fit

ndres <- tibble(unique(datares[c('inf','Cross')]))%>%mutate(eggs = 1) 
postres <- fitted(fitres, newdata = ndres, re_formula = NA, probs = c(0.09, 0.91))%>%
  as_tibble()%>%mutate(ndres[c('Cross', 'inf')])

ggplot(postres, aes(x = Cross, y = Estimate))+
  geom_boxplot(data = datares, aes(x = Cross, y = CIobs, color = inf), position = position_dodge(0.8), width = 0.7)+
  geom_point(aes(class = inf), position = position_dodge(0.8))+
  geom_errorbar(aes(ymin=Q9, ymax=Q91, class = inf), width=.2, position=position_dodge(0.8))+
  labs(title = 'proportion eggs succesfully developing in adult females', y = 'proportion female', x = 'cross')+
  ylim(0,1)+
  scale_color_discrete(name = "", labels = c('cured', expression(paste(italic("Wolbachia"),"-infected"))))+
  theme(legend.position = 'bottom')

postres_ <- fitted(fitres, newdata=ndres, re_formula = NA, summary = 0)%>%as_tibble()
postresCI2__ <- 1-postres_[seq(2,nrow(ndres),2)]/postres_[seq(1,nrow(ndres),2)]%>%as_tibble()
postresCI2_wide <- postresCI2__%>%mutate(across(everything(), function(x) {pmax(x, 0)}))
postresCI2_ <- pivot_longer(postresCI2_wide, cols = everything(), names_to = 'group', names_prefix = 'V', values_to = 'CI')%>%
  mutate(group = as.numeric(group))%>%mutate(Cross = rep(ndres$Cross[seq(1,nrow(ndres),2)], nrow(postresCI2__)))
postresCI2 <- group_by(postresCI2_, group)%>%summarise(CImean = mean(CI), CIq09 = quantile(CI, probs = 0.09), CIq91 = quantile(CI, probs = 0.91))%>%
  bind_cols(tibble(unique(datares[c('Cross')])))


ggplot(postresCI2_, aes(x = Cross, y = CI))+
  theme_light()+
  geom_violin(draw_quantiles = c(0.09, 0.5, 0.91), alpha = 0.5, scale = 'width')+
  scale_fill_brewer(palette = 'Dark2')+
  ylim(0,1)+
  labs(title = 'Cytoplasmatic incompatibility', y = 'CI', x = 'Cross')


#rescued lines with unrescued
datafres <- data__[data__$Cross %in% c('Beis_x_Beis', 'Scp_x_Scp'),]
datafres <- droplevels(datafres)
datafres$finf <- datafres$Type %in% c('R', 'Rc')

formfres <- brmsformula(adult_females | trials(eggs)  ~ 1 + inf*finf*Cross, family = binomial)
get_prior(formfres, data = datafres)
fitfres <- brm(data = datafres, family = binomial,
              formula = formfres,
              prior = c(prior(normal(0,10), class = Intercept),
                        prior(normal(0,10), class = b)),
              iter = 5000, warmup = 2000, chains = 2, cores = 2)
fitfres <- add_criterion(fitfres, 'waic')
fitfres$fit

ndfres <- tibble(unique(datafres[c('inf', 'finf','Cross')]))%>%mutate(eggs = 1) 
postfres <- fitted(fitfres, newdata = ndfres, re_formula = NA, probs = c(0.09, 0.91))%>%
  as_tibble()%>%mutate(ndfres[c('Cross', 'inf', 'finf')])

ggplot(postfres, aes(x = Cross, y = Estimate))+
  geom_boxplot(data = datafres, aes(x = Cross, y = CIobs, color = inf), position = position_dodge(0.8), width = 0.7)+
  geom_point(aes(class = inf), position = position_dodge(0.8))+
  geom_errorbar(aes(ymin=Q9, ymax=Q91, class = inf), width=.2, position=position_dodge(0.8))+
  facet_wrap('finf')+
  labs(title = 'proportion eggs succesfully developing in adult females', y = 'proportion female', x = 'cross')+
  ylim(0,1)+
  scale_color_discrete(name = "", labels = c('cured', expression(paste(italic("Wolbachia"),"-infected"))))+
  theme(legend.position = 'bottom')

postfres_ <- fitted(fitfres, newdata=ndfres, re_formula = NA, summary = 0)%>%as_tibble()
postfresCI2__ <- 1-postfres_[seq(2,nrow(ndfres),2)]/postfres_[seq(1,nrow(ndfres),2)]%>%as_tibble()
postfresCI2_wide <- postfresCI2__%>%mutate(across(everything(), function(x) {pmax(x, 0)}))
postfresCI2_ <- pivot_longer(postfresCI2_wide, cols = everything(), names_to = 'group', names_prefix = 'V', values_to = 'CI')%>%
  mutate(group = as.numeric(group))%>%mutate(Cross = rep(ndfres$Cross[seq(1,nrow(ndfres),2)], nrow(postfresCI2__)), finf = rep(ndfres$finf[seq(1,nrow(ndfres),2)], nrow(postfresCI2__)))
postfresCI2 <- group_by(postfresCI2_, group)%>%summarise(CImean = mean(CI), CIq09 = quantile(CI, probs = 0.09), CIq91 = quantile(CI, probs = 0.91))%>%
  bind_cols(tibble(unique(datafres[c('Cross', 'finf')])))


ggplot(postfresCI2_, aes(x = Cross, y = CI))+
  theme_light()+
  geom_violin(aes(fill = finf), draw_quantiles = c(0.09, 0.5, 0.91), alpha = 0.5, scale = 'width', position = position_dodge(0.8))+
  scale_fill_brewer(palette = 'Dark2', name = "", labels = c('incompatible', 'cured'))+
  ylim(0,1)+
  coord_flip()+
  labs(title = 'Cytoplasmatic incompatibility', y = 'CI', x = 'Cross')+
  theme(legend.position = 'bottom')