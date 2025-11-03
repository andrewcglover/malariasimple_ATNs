## DECLARE TIME STEPS
n_days <- parameter()

dim(days) <- n_days + 1
days <- parameter()

dim(daily_rain_input) <- n_days +1
daily_rain_input <- parameter()
rain_input <- interpolate(days, daily_rain_input, "linear")

dim(daily_ft) <- n_days + 1
daily_ft <- parameter()
ft <- interpolate(days, daily_ft, "linear")

## MODEL VARIABLES
##------------------------------------------------------------------------------

na <- parameter()
nh <- parameter()

##------------------------------------------------------------------------------
##################
## HUMAN STATES ##
##################
##------------------------------------------------------------------------------
# Human states as specified in full transmission model
# http://journals.plos.org/plosmedicine/article?id=10.1371%2Fjournal.pmed.1000324
# http://www.nature.com/articles/ncomms4136

#--------- DEFINE PARAMETERS FOR HUMAN TRANSITIONS ---------
eta <- parameter()
age_rate <- parameter()
dim(age_rate) <- na
het_wt <- parameter()
dim(het_wt) <- nh
rA <- parameter()
rT <- parameter()
rD <- parameter()
rU <- parameter()
rP <- parameter()

#----------- INITIALISE HUMAN COMPARTMENTS -----------------
# S - SUSCEPTIBLE
init_S <- parameter()
dim(init_S) <- c(na,nh,num_int)
initial(S[,,]) <- init_S[i,j,k]
dim(S) <- c(na,nh,num_int)

# T- SUCCESSFULLY TREATED
init_T <- parameter()
dim(init_T) <- c(na,nh,num_int)
initial(T[,,]) <- init_T[i,j,k]
dim(T) <- c(na,nh,num_int)

# D - CLEAR DISEASE
init_D <- parameter()
dim(init_D) <- c(na,nh,num_int)
initial(D[,,]) <- init_D[i,j,k]
dim(D) <- c(na,nh,num_int)

# A - ASYMPTOMATIC DISEASE
init_A <- parameter()
dim(init_A) <- c(na,nh,num_int)
initial(A[,,]) <- init_A[i,j,k]
dim(A) <- c(na,nh,num_int)

# U - SUBPATENT DISEASE
init_U <- parameter()
dim(init_U) <- c(na,nh,num_int)
initial(U[,,]) <- init_U[i,j,k]
dim(U) <- c(na,nh,num_int)

# P - PROPHYLAXIS
init_P <- parameter()
dim(init_P) <- c(na,nh,num_int)
initial(P[,,]) <- init_P[i,j,k]
dim(P) <- c(na,nh,num_int)

#-------------------- DEFINE HUMAN TRANSITION EVENTS -------------------------
#### BIRTHS ####
dim(births) <- c(1,nh,num_int)
dim(birth_rate) <- c(1,nh,num_int)
dim(birth_prob) <- c(1,nh,num_int)
birth_rate[1,,] <- cov[k]*eta*het_wt[j]
birth_prob[1,,] <- 1-exp(-(birth_rate[1,j,k]*dt))
births[1,,] <- Binomial(H,birth_prob[1,j,k])

#
# #### TRANSITIONS FROM S #####
dim(S_rates) <- c(na,nh,num_int,5)
S_rates[,,,1] <- (1-phi[i,j,k])*FOI_smc[i,j,k] #S -> A
S_rates[,,,2] <- age_rate[i] #Aging
S_rates[,,,3] <- ft*phi[i,j,k]*FOI_smc[i,j,k] #S -> T
S_rates[,,,4] <- (1-ft)*phi[i,j,k]*FOI_smc[i,j,k] #S -> D
S_rates[,,,5] <- eta #Death

#Calculate the total number of individuals leaving the 'S' compartment each time step
dim(S_leave_rate) <- c(na,nh,num_int)
dim(S_leave_prob) <- c(na,nh,num_int)
dim(S_leave) <- c(na,nh,num_int)
S_leave_rate[,,] <- sum(S_rates[i,j,k,1:5]) #Rate of individuals leaving 'S' compartment is equal to the sum of each individual transition rate
S_leave_prob[,,] <- 1 - exp(-(S_leave_rate[i,j,k]*dt)) #Per capita probability of an individual leaving the 'S' compartment per time step.
S_leave[, , ] <- Binomial(S[i, j, k], S_leave_prob[i, j, k])

##Individual transitions are stored in an array. Multinomials are not supported by odin_dust, hence a sequence of binomials is used
dim(S_trans) <- c(na,nh,num_int,5)
S_trans[,,,1] <- Binomial(S_leave[i,j,k], S_rates[i,j,k,1] / S_leave_rate[i,j,k])
S_trans[,,,2:5] <- Binomial(S_leave[i,j,k] - sum(S_trans[i, j, k, 1:(l - 1)]), S_rates[i, j, k, l]/sum(S_rates[i, j, k, l:5]))

##Individual transitions are named for interpretability.
dim(S_death) <- c(na,nh,num_int)
dim(S_age) <- c(na,nh,num_int)
dim(ST_trans) <- c(na,nh,num_int)
dim(SD_trans) <- c(na,nh,num_int)
dim(SA_trans) <- c(na,nh,num_int)

SA_trans[,,] <- S_trans[i,j,k,1]
S_age[,,] <- S_trans[i,j,k,2]
ST_trans[,,] <- S_trans[i,j,k,3]
SD_trans[,,] <- S_trans[i,j,k,4]
S_death[,,] <- S_trans[i,j,k,5]

#### TRANSITIONS FROM T #####
dim(T_rates) <- c(na,nh,num_int,3)
T_rates[,,,1] <- rT #T -> P
T_rates[,,,2] <- age_rate[i] #Aging
T_rates[,,,3] <- eta #Death

dim(T_leave_rate) <- c(na,nh,num_int)
dim(T_leave_prob) <- c(na,nh,num_int)
dim(T_leave) <- c(na,nh,num_int)
T_leave_rate[,,] <- sum(T_rates[i,j,k,1:3])
T_leave_prob[,,] <- 1 - exp(-(T_leave_rate[i,j,k]*dt))
T_leave[, , ] <- Binomial(T[i, j, k], T_leave_prob[i, j, k])

dim(T_trans) <- c(na,nh,num_int,3)
T_trans[,,,1] <- Binomial(T_leave[i,j,k], T_rates[i,j,k,1] / T_leave_rate[i,j,k])
T_trans[, , , 2:3] <- Binomial(T_leave[i, j, k] - sum(T_trans[i, j, k, 1:(l - 1)]), T_rates[i, j, k, l]/sum(T_rates[i, j, k, l:3]))

dim(T_death) <- c(na,nh,num_int)
dim(T_age) <- c(na,nh,num_int)
dim(TP_trans) <- c(na,nh,num_int)

TP_trans[,,] <- T_trans[i,j,k,1]
T_age[,,] <- T_trans[i,j,k,2]
T_death[,,] <- T_trans[i,j,k,3]

#### TRANSITIONS FROM D #####
dim(D_rates) <- c(na,nh,num_int,3)
D_rates[,,,1] <- age_rate[i] * dt#Aging
D_rates[,,,2] <- rD * dt#D -> A
D_rates[,,,3] <- eta * dt#Death

dim(D_smc_cleared) <- c(na,nh,num_int)
D_smc_cleared[,,] <- Binomial(D[i,j,k], alpha_smc_array[i,j,k])

dim(D_leave_rate) <- c(na,nh,num_int)
dim(D_leave_prob) <- c(na,nh,num_int)
dim(D_leave) <- c(na,nh,num_int)
D_leave_rate[,,] <- sum(D_rates[i,j,k,1:3])
D_leave_prob[,,] <- 1 - exp(-(D_leave_rate[i,j,k]))
D_leave[, , ] <- Binomial(D[i, j, k] - D_smc_cleared[i,j,k], D_leave_prob[i, j, k])

dim(D_trans) <- c(na,nh,num_int,3)
D_trans[,,,1] <- Binomial(D_leave[i,j,k], D_rates[i,j,k,1] / D_leave_rate[i,j,k])
D_trans[, , , 2:3] <- Binomial(D_leave[i, j, k] - sum(D_trans[i, j, k, 1:(l - 1)]), D_rates[i, j, k, l]/sum(D_rates[i, j, k, l:3]))

dim(D_death) <- c(na,nh,num_int)
dim(D_age) <- c(na,nh,num_int)
dim(DA_trans) <- c(na,nh,num_int)

D_age[,,] <- D_trans[i,j,k,1]
DA_trans[,,] <- D_trans[i,j,k,2]
D_death[,,] <- D_trans[i,j,k,3]

#### TRANSITIONS FROM A #####
dim(A_rates) <- c(na,nh,num_int,5)
A_rates[,,,1] <- age_rate[i] * dt #Aging
A_rates[,,,2] <- rA * dt#A -> U
A_rates[,,,3] <- phi[i,j,k] * FOI_smc[i,j,k] * ft * dt #A -> T
A_rates[,,,4] <- phi[i,j,k] * FOI_smc[i,j,k] * (1-ft) * dt #A -> D
A_rates[,,,5] <- eta * dt#Death

dim(A_smc_cleared) <- c(na,nh,num_int)
A_smc_cleared[,,] <- Binomial(A[i,j,k], alpha_smc_array[i,j,k])

dim(A_leave) <- c(na,nh,num_int)
dim(A_leave_rate) <- c(na,nh,num_int)
dim(A_leave_prob) <- c(na,nh,num_int)

A_leave_rate[,,] <- sum(A_rates[i,j,k,1:5])
A_leave_prob[,,] <- 1 - exp(-(A_leave_rate[i,j,k]))
A_leave[, , ] <- Binomial(A[i, j, k] - A_smc_cleared[i, j, k], A_leave_prob[i, j, k])

dim(A_trans) <- c(na,nh,num_int,5)
A_trans[,,,1] <- Binomial(A_leave[i,j,k], A_rates[i,j,k,1] / A_leave_rate[i,j,k])
A_trans[, , , 2:5] <- Binomial(A_leave[i, j, k] - sum(A_trans[i, j, k, 1:(l - 1)]), A_rates[i, j, k, l]/sum(A_rates[i, j, k, l:5]))

dim(A_death) <- c(na,nh,num_int)
dim(A_age) <- c(na,nh,num_int)
dim(AU_trans) <- c(na,nh,num_int)
dim(AT_trans) <- c(na,nh,num_int)
dim(AD_trans) <- c(na,nh,num_int)

A_age[,,] <- A_trans[i,j,k,1]
AU_trans[,,] <- A_trans[i,j,k,2]
AT_trans[,,] <- A_trans[i,j,k,3]
AD_trans[,,] <- A_trans[i,j,k,4]
A_death[,,] <- A_trans[i,j,k,5]

#### TRANSITIONS FROM U #####
dim(U_rates) <- c(na,nh,num_int,7)
U_rates[,,,1] <- age_rate[i] * dt #Aging
U_rates[,,,2] <- rU * dt #U -> S
U_rates[,,,3] <- (1-phi[i,j,k]) * FOI_smc[i,j,k] * dt #U -> A
U_rates[,,,4] <- phi[i,j,k] * (1-ft) * FOI_smc[i,j,k] * dt #U -> D
U_rates[,,,5] <- phi[i,j,k] * ft * FOI_smc[i,j,k] * dt #U -> T
U_rates[,,,6] <- eta * dt #Death

#Infection clearance from SMC
dim(U_smc_cleared) <- c(na,nh,num_int)
U_smc_cleared[,,] <- Binomial(U[i,j,k], alpha_smc_array[i,j,k])

dim(U_leave_rate) <- c(na,nh,num_int)
dim(U_leave_prob) <- c(na,nh,num_int)
dim(U_leave) <- c(na,nh,num_int)

U_leave_rate[,,] <- sum(U_rates[i,j,k,1:6])
U_leave_prob[,,] <- 1 - exp(-(U_leave_rate[i,j,k]))
U_leave[, , ] <- Binomial(U[i, j, k] - U_smc_cleared[i,j,k], U_leave_prob[i, j, k])

dim(U_trans) <- c(na,nh,num_int,6)
U_trans[,,,1] <- Binomial(U_leave[i,j,k], U_rates[i,j,k,1] / U_leave_rate[i,j,k])
U_trans[, , , 2:6] <- Binomial(U_leave[i, j, k] - sum(U_trans[i, j, k, 1:(l - 1)]), U_rates[i, j, k, l]/sum(U_rates[i, j, k, l:6]))

dim(U_death) <- c(na,nh,num_int)
dim(U_age) <- c(na,nh,num_int)
dim(US_trans) <- c(na,nh,num_int)
dim(UA_trans) <- c(na,nh,num_int)
dim(UD_trans) <- c(na,nh,num_int)
dim(UT_trans) <- c(na,nh,num_int)

U_age[,,]    <- U_trans[i,j,k,1]
US_trans[,,] <- U_trans[i,j,k,2]
UA_trans[,,] <- U_trans[i,j,k,3]
UD_trans[,,] <- U_trans[i,j,k,4]
UT_trans[,,] <- U_trans[i,j,k,5]
U_death[,,]  <- U_trans[i,j,k,6]

#### TRANSITIONS FROM P #####
dim(P_rates) <- c(na,nh,num_int,3)
P_rates[,,,1] <- rP * dt#P -> S
P_rates[,,,2] <- age_rate[i] * dt #Aging
P_rates[,,,3] <- eta * dt #Death

dim(P_leave_rate) <- c(na,nh,num_int)
dim(P_leave_prob) <- c(na,nh,num_int)
dim(P_leave) <- c(na,nh,num_int)
P_leave_rate[,,] <- sum(P_rates[i,j,k,1:3])
P_leave_prob[,,] <- 1 - exp(-(P_leave_rate[i,j,k]))
P_leave[, , ] <- Binomial(P[i, j, k], P_leave_prob[i, j, k])

dim(P_trans) <- c(na,nh,num_int,3)
P_trans[,,,1] <- Binomial(P_leave[i,j,k], P_rates[i,j,k,1] / P_leave_rate[i,j,k])
P_trans[, , , 2:3] <- Binomial(P_leave[i, j, k] - sum(P_trans[i, j, k, 1:(l - 1)]), P_rates[i, j, k, l]/sum(P_rates[i, j, k, l:3]))

dim(PS_trans) <- c(na,nh,num_int)
dim(P_death) <- c(na,nh,num_int)
dim(P_age) <- c(na,nh,num_int)

PS_trans[,,] <- P_trans[i,j,k,1]
P_age[,,] <- P_trans[i,j,k,2]
P_death[,,] <- P_trans[i,j,k,3]

#------------------------------ TRANSITIONS ----------------------------------
update(S[1,1:nh,1:num_int]) <- S[i,j,k] + PS_trans[i,j,k] + US_trans[i,j,k] + U_smc_cleared[i,j,k] + A_smc_cleared[i,j,k] + D_smc_cleared[i,j,k] - ST_trans[i,j,k] - SD_trans[i,j,k] - SA_trans[i,j,k] - S_death[i,j,k] - S_age[i,j,k] + births[1,j,k]
update(S[2:na,1:nh,1:num_int]) <- S[i,j,k] + PS_trans[i,j,k] + US_trans[i,j,k] + U_smc_cleared[i,j,k] + A_smc_cleared[i,j,k] + D_smc_cleared[i,j,k] - ST_trans[i,j,k] - SD_trans[i,j,k] - SA_trans[i,j,k] - S_death[i,j,k] - S_age[i,j,k] + S_age[i-1,j,k]

update(T[1,1:nh,1:num_int]) <- T[i,j,k] + AT_trans[i,j,k] + ST_trans[i,j,k] + UT_trans[i,j,k] - TP_trans[i,j,k] - T_age[i,j,k] - T_death[i,j,k]
update(T[2:na,1:nh,1:num_int]) <- T[i,j,k] + AT_trans[i,j,k] + ST_trans[i,j,k] + UT_trans[i,j,k] - TP_trans[i,j,k] - T_age[i,j,k] - T_death[i,j,k] + T_age[i-1,j,k]

update(D[1,1:nh,1:num_int]) <- D[i,j,k] + SD_trans[i,j,k] + AD_trans[i,j,k] + UD_trans[i,j,k] - DA_trans[i,j,k] - D_smc_cleared[i,j,k] - D_death[i,j,k] - D_age[i,j,k]
update(D[2:na,1:nh,1:num_int]) <- D[i,j,k] + SD_trans[i,j,k] + AD_trans[i,j,k] + UD_trans[i,j,k] - DA_trans[i,j,k] - D_smc_cleared[i,j,k] - D_death[i,j,k] - D_age[i,j,k] + D_age[i-1,j,k]

update(A[1,1:nh,1:num_int]) <- A[i,j,k] + SA_trans[i,j,k] + DA_trans[i,j,k] + UA_trans[i,j,k] - AT_trans[i,j,k] - AD_trans[i,j,k] - AU_trans[i,j,k] - A_smc_cleared[i,j,k] - A_death[i,j,k] - A_age[i,j,k]
update(A[2:na,1:nh,1:num_int]) <- A[i,j,k] + SA_trans[i,j,k] + DA_trans[i,j,k] + UA_trans[i,j,k] - AT_trans[i,j,k] - AD_trans[i,j,k] - AU_trans[i,j,k] - A_smc_cleared[i,j,k] - A_death[i,j,k] - A_age[i,j,k] + A_age[i-1,j,k]

update(U[1,1:nh,1:num_int]) <- U[i,j,k] + AU_trans[i,j,k] - UD_trans[i,j,k] - UT_trans[i,j,k] - US_trans[i,j,k] - U_smc_cleared[i,j,k] - UA_trans[i,j,k] - U_age[i,j,k] - U_death[i,j,k]
update(U[2:na,1:nh,1:num_int]) <- U[i,j,k] + AU_trans[i,j,k] - UD_trans[i,j,k] - UT_trans[i,j,k] - US_trans[i,j,k] - U_smc_cleared[i,j,k] - UA_trans[i,j,k] - U_age[i,j,k] - U_death[i,j,k] + U_age[i-1,j,k]

update(P[1,1:nh,1:num_int]) <- P[i,j,k] + TP_trans[i,j,k] - PS_trans[i,j,k] - P_death[i,j,k] - P_age[i,j,k]
update(P[2:na,1:nh,1:num_int]) <- P[i,j,k] + TP_trans[i,j,k] - PS_trans[i,j,k] - P_death[i,j,k] - P_age[i,j,k] + P_age[i-1,j,k]

##------------------------------------------------------------------------------
##########################################
## KEY INFECTION PARAMETERS (FOI & EIR) ##
##########################################
##------------------------------------------------------------------------------

#----------FORCE OF INFECTION-----------------
# Force of infection, depends on level of infection blocking immunity
dim(FOI_lag) <- c(na,nh,num_int)
FOI_lag[1:na, 1:nh, 1:num_int] <- EIR[i,j,k] * (if(IB[i,j,k]==0) b0 else b[i,j,k])

# Current FOI depends on humans that have been through the latent period
dE <- parameter()
lag_rates <- parameter(type = "integer")

dim(FOI_eq) <- c(na,nh)
FOI_eq <- parameter()

dim(FOI_XL) <- c(na,nh,num_int,lag_rates)
initial(FOI_XL[,,,]) <- FOI_eq[i,j]

update(FOI_XL[,,,1]) <- FOI_XL[i,j,k,1] + dt*((lag_rates/dE)*FOI_lag[i,j,k] - (lag_rates/dE)*FOI_XL[i,j,k,1])
update(FOI_XL[,,,2:lag_rates]) <- FOI_XL[i,j,k,l] + dt*((lag_rates/dE)*FOI_XL[i,j,k,l-1] - (lag_rates/dE)*FOI_XL[i,j,k,l])

dim(FOI) <- c(na,nh,num_int)
FOI[,,] <- FOI_XL[i,j,k,lag_rates]


#------------ EIR -----------------
# rate for age group * rate for biting category * FOI for age group * prop of
# infectious mosquitoes
dim(foi_age) <- na
foi_age <- parameter()
dim(rel_foi) <- nh
rel_foi <- parameter()
dim(EIR) <- c(na,nh,num_int)
EIR[,,] <- av_mosq[k] * rel_foi[j] * foi_age[i] * Iv/omega

##------------------------------------------------------------------------------
#####################
## IMMUNITY STATES ##
#####################
##------------------------------------------------------------------------------

# See supplementary materials S1 from http://journals.plos.org/plosmedicine/article?id=10.1371/journal.pmed.1000324#s6

# ICM - Maternally acquired immunity acquired by babies from mothers (assumed to be proportional to the immunity of a 15 to 30 year old woman)
# ICA - Immunity acquired due to exposure to previous infection, increases with age
# IC - Clinical immunity. Upon infection, immunity against clinical case. IC = ICA + ICM
# IB - Infection blocking immunity, chances of preventing infection upon receiving infectious bite
# ID - Detection immunity, when immunity suppresses parasite densities this makes it less likely that diagnostics will detect parasite infection

#----------------- DEFINE IMMUNITY PARAMETERS -------------------
dCM <- parameter()
uCA <- parameter()
dCA <- parameter()
dB <- parameter()
uB <- parameter()
dID <- parameter()
uD <- parameter()
x_I <- parameter()
dim(x_I) <- na
age20l <- parameter(type = "integer")
age20u <- parameter(type = "integer")
age_20_factor <- parameter()
PM <- parameter()

#---------- UPDATE IMMUNITY STATES FOR EACH HUMAN COMPARTMENTS ---------

# ICM - maternally acquired immunity
init_ICM <- parameter()
dim(init_ICM) <- c(na,nh,num_int)
initial(ICM[,,]) <- init_ICM[i,j,k]
dim(ICM) <- c(na,nh,num_int)
dim(init_ICM_pre) <- c(nh,num_int)
init_ICM_pre[1:nh,1:num_int] <- PM*(ICA[age20l,i,j] + age_20_factor*(ICA[age20u,i,j]-ICA[age20l,i,j]))

update(ICM[1, 1:nh, 1:num_int]) <- ICM[i,j,k] + dt*(-1/dCM*ICM[i,j,k] + (init_ICM_pre[j,k]-ICM[i,j,k])/x_I[i])
update(ICM[2:na, 1:nh, 1:num_int]) <- ICM[i,j,k] + dt*(-1/dCM*ICM[i,j,k] - (ICM[i,j,k]-ICM[i-1,j,k])/x_I[i])

# ICA - exposure driven immunity
init_ICA <- parameter()
dim(init_ICA) <- c(na,nh,num_int)
initial(ICA[,,]) <- init_ICA[i,j,k]
dim(ICA) <- c(na,nh,num_int)

update(ICA[1, 1:nh, 1:num_int]) <- ICA[i,j,k] + dt*(FOI_smc[i,j,k]/(FOI_smc[i,j,k] * uCA + 1) - 1/dCA*ICA[i,j,k] -ICA[i,j,k]/x_I[i])
update(ICA[2:na, 1:nh, 1:num_int]) <- ICA[i,j,k] + dt*(FOI_smc[i,j,k]/(FOI_smc[i,j,k] * uCA + 1) - 1/dCA*ICA[i,j,k] - (ICA[i,j,k]-ICA[i-1,j,k])/x_I[i])

dim(IC) <- c(na,nh,num_int)
IC[,,] <- ICM[i,j,k] + ICA[i,j,k]

# IB - infection blocking immunity
init_IB <- parameter()
dim(init_IB) <- c(na,nh,num_int)
initial(IB[,,]) <- init_IB[i,j,k]
dim(IB) <- c(na,nh,num_int)

update(IB[1, 1:nh, 1:num_int]) <- IB[i,j,k] + dt*(EIR[i,j,k]/(EIR[i,j,k]* uB + 1) - IB[i,j,k]/dB - IB[i,j,k]/x_I[i])
update(IB[2:na, 1:nh, 1:num_int]) <- IB[i,j,k] + dt*(EIR[i,j,k]/(EIR[i,j,k]* uB + 1) - IB[i,j,k]/dB - (IB[i,j,k]-IB[i-1,j,k])/x_I[i])

# ID - detection immunity
init_ID <- parameter()
dim(init_ID) <- c(na,nh,num_int)
initial(ID[,,]) <- init_ID[i,j,k]
dim(ID) <- c(na,nh,num_int)

update(ID[1, 1:nh, 1:num_int]) <- ID[i,j,k] + dt*(FOI_smc[i,j,k]/(FOI_smc[i,j,k]*uD + 1) - ID[i,j,k]/dID - ID[i,j,k]/x_I[i])
update(ID[2:na, 1:nh, 1:num_int]) <- ID[i,j,k] + dt*(FOI_smc[i,j,k]/(FOI_smc[i,j,k]*uD + 1) - ID[i,j,k]/dID - (ID[i,j,k]-ID[i-1,j,k])/x_I[i])

#----------- PROBABILITIES INFERRED FROM IMMUNITY --------------------
# phi - probability of clinical disease, dependent on clinical immunity
phi0 <- parameter()
phi1 <- parameter()
IC0 <- parameter()
kC <- parameter()
dim(phi) <- c(na,nh,num_int)
phi[1:na,1:nh,1:num_int] <- phi0*((1-phi1)/(1+(IC[i,j,k]/IC0)^kC) + phi1)

# b - probability of disease from infectious bite, depends on infection blocking immunity
b0 <- parameter()
b1 <- parameter()
kB <- parameter()
IB0 <- parameter()
dim(b) <- c(na,nh,num_int)
b[1:na, 1:nh, 1:num_int] <- b0 * ((1-b1)/(1+(IB[i,j,k]/IB0)^kB)+b1)


# p_det - probability of detection by microscopy, immunity decreases chances of
# infection because it pushes parasite density down
aD <- parameter()
fD0 <- parameter()
gammaD <- parameter()
d1 <- parameter()
ID0 <- parameter()
kD <- parameter()
dim(age_vector) <- na
age_vector <- parameter()

dim(fd) <- na
fd[1:na] <- 1-(1-fD0)/(1+(age_vector[i]/aD)^gammaD)
dim(p_det) <- c(na,nh,num_int)
p_det[,,] <- d1 + (1-d1)/(1 + fd[i]*(ID[i,j,k]/ID0)^kD)

##------------------------------------------------------------------------------
#####################
## MOSQUITO STATES ##
#####################
##------------------------------------------------------------------------------

# See supplementary materials S1 from http://journals.plos.org/plosmedicine/article?id=10.1371/journal.pmed.1000324#s6

# Sv - Susceptible mosquitoes
# Pv - Latently infected (Pre-Infectious?) mosquitoes. Number of compartments used to simulate delay in becoming infectious
# Iv - Infectious mosquitoes

# initial state values:
init_Sv <- parameter()
init_Pv <- parameter()
init_Iv <- parameter()
initial(Sv) <- init_Sv * mv0
initial(Pv) <- init_Pv * mv0
initial(Iv) <- init_Iv * mv0

# cA is the infectiousness to mosquitoes of humans in the asmyptomatic compartment broken down
# by age/het/int category, infectiousness depends on p_det which depends on detection immunity
cU <- parameter()
cD <- parameter()
cT <- parameter()
gamma1 <- parameter()
dim(cA) <- c(na,nh,num_int)
cA[,,] <- cU + (cD-cU)*p_det[i,j,k]^gamma1

# Force of infection from humans to mosquitoes
lag_ratesMos <- parameter(type = "integer")

FOIv_eq <- parameter()
dim(FOIv) <- lag_ratesMos
initial(FOIv[]) <- FOIv_eq*delayGam/lag_ratesMos

update(FOIv[1]) <- FOIv[1] + dt*(lag_FOIv - (lag_ratesMos/delayGam)*FOIv[1])
update(FOIv[2:lag_ratesMos]) <- FOIv[i] + dt*((lag_ratesMos/delayGam)*FOIv[i-1] - (lag_ratesMos/delayGam)*FOIv[i])


dim(FOIvijk) <- c(na,nh,num_int)
omega <- parameter()
FOIvijk[1:na, 1:nh, 1:num_int] <- ((cT*smc_rel_c_mask[i,j,k]*T[i,j,k] + cD*smc_rel_c_mask[i,j,k]*D[i,j,k] + cA[i,j,k]*smc_rel_c_mask[i,j,k]*A[i,j,k] + cU*smc_rel_c_mask[i,j,k]*U[i,j,k])/H) *
  rel_foi[j] * av_mosq[k]*foi_age[i]/omega ## For discrete human compartments
lag_FOIv=sum(FOIvijk)

ince <- FOIv[lag_ratesMos] * lag_ratesMos/delayGam * Sv

initial(ince_delay[]) <- FOIv_eq*init_Sv*mv0*delayMos_use/lag_ratesMos
dim(ince_delay) <- lag_ratesMos

update(ince_delay[1]) <- ince_delay[1] + dt*(ince - (lag_ratesMos/delayMos_use)*ince_delay[1])
update(ince_delay[2:lag_ratesMos]) <- ince_delay[i] + dt*((lag_ratesMos/delayMos_use)*ince_delay[i-1] -
                                                            (lag_ratesMos/delayMos_use)*ince_delay[i])

incv <- ince_delay[lag_ratesMos]*lag_ratesMos/delayMos_use * surv

# Current hum->mos FOI depends on the number of individuals now producing gametocytes (12 day lag)
delayGam <- parameter()
delayMos <- parameter()
delayMos_use <- delayMos
# Number of mosquitoes that become infected at each time point
surv <- exp(-mu*delayMos_use)

# Number of mosquitoes born (depends on PL, number of larvae), or is constant outside of seasonality
betaa <- 0.5*PL/dPL

update(Sv) <- if(Sv + dt*(-ince - mu*Sv + betaa) < 0) 0 else Sv + dt*(-ince - mu*Sv + betaa)
update(Pv) <- if(Pv + dt*(ince - incv - mu*Pv) < 0) 0 else Pv + dt*(ince - incv - mu*Pv)
update(Iv) <- if(Iv + dt*(incv - mu*Iv) < 0) 0 else Iv + dt*(incv - mu*Iv)

# Total mosquito population
initial(mv) <- 0
update(mv) <- Sv+Pv+Iv

##------------------------------------------------------------------------------
###################
## LARVAL STATES ##
###################
##------------------------------------------------------------------------------

# Model by White et al.
# (https://parasitesandvectors.biomedcentral.com/articles/10.1186/1756-3305-4-153)

# EL - early larval instar stage
# LL - late larval instar stage
# PL - pupal stage

# mean carrying capacity from initial mosquito density:
dLL <- parameter()
dPL <- parameter()
dEL <- parameter()
muLL <- parameter()
muEL <- parameter()
muPL <- parameter()
gammaL <- parameter()

# fitted entomological parameters:
mv0 <- parameter()
mum <- parameter()
foraging_time <- parameter()
gonotrophic_cycle <- parameter()
betaL <- parameter()

# Entomological variables:
mum_use <- mum
p10 <- exp(-mum_use * foraging_time)  # probability of surviving one feeding cycle
p2 <- exp(-mum_use * gonotrophic_cycle)  # probability of surviving one resting cycle
eov <- betaL/mu*(exp(mu/fv)-1)
beta_larval <- eov*mu*exp(-mu/fv)/(1-exp(-mu/fv)) # Number of eggs laid per day
b_lambda <- (gammaL*muLL/muEL-dEL/dLL+(gammaL-1)*muLL*dEL)
lambda <- -0.5*b_lambda + sqrt(0.25*b_lambda^2 + gammaL*beta_larval*muLL*dEL/(2*muEL*mum_use*dLL*(1+dPL*muPL)))
K0 <- 2*mv0*dLL*mum_use*(1+dPL*muPL)*gammaL*(lambda+1)/(lambda/(muLL*dEL)-1/(muLL*dLL)-1)

# Seasonal carrying capacity KL = base carrying capacity K0 * effect for time of year:
KL <- K0*rain_input
fv <- 1/( foraging_time/(1-zbar) + gonotrophic_cycle ) # mosquito feeding rate (zbar from intervention param.)
mu <- -fv*log(p1*p2) # mosquito death rate

# finding equilibrium and initial values for EL, LL & PL
init_PL <- parameter()
initial(PL) <- init_PL
init_LL <- parameter()
initial(LL) <- init_LL
init_EL <- parameter()
initial(EL) <- init_EL

# (beta_larval (egg rate) * total mosquito) - den. dep. egg mortality - egg hatching
update(EL) <- if(EL + dt*(beta_larval*mv-muEL*(1+(EL+LL)/KL)*EL - EL/dEL) < 0) 0 else (EL + dt*(beta_larval*mv-muEL*(1+(EL+LL)/KL)*EL - EL/dEL))
# egg hatching - den. dep. mortality - maturing larvae
update(LL) <- if(LL + dt*(EL/dEL - muLL*(1+gammaL*(EL + LL)/KL)*LL - LL/dLL) < 0) 0 else (LL + dt*(EL/dEL - muLL*(1+gammaL*(EL + LL)/KL)*LL - LL/dLL))
# pupae - mortality - fully developed pupae
update(PL) <- if(PL + dt*(LL/dLL - muPL*PL - PL/dPL) < 0) 0 else (PL + dt*(LL/dLL - muPL*PL - PL/dPL))


##------------------------------------------------------------------------------
########################
## INTERVENTION MODEL ##
########################
##------------------------------------------------------------------------------

################## SMC #######################
#See supplementary materials of [Thompson, 2022] - https://doi.org/10.1016/S2214-109X(22)00416-8
max_smc_cov <- parameter()

dim(alpha_smc_times, alpha_smc_set) <- parameter(rank = 1)
alpha_smc_times <- parameter(constant = TRUE)
alpha_smc_set <- parameter(constant = TRUE)

alpha_smc <- interpolate(alpha_smc_times, alpha_smc_set, mode = "constant")

dim(alpha_smc_array) <-  c(na,nh,num_int)
alpha_smc_array[,,] <- smc_mask[i,j,k] * alpha_smc #Multiply by proportion of individuals in SMC compartment currently receiving SMC

##Parameters relating to prophylaxis effect of SMC
dim(P_smc_daily) <- n_days + 1
P_smc_daily <- parameter()
P_smc <- interpolate(days, P_smc_daily, "linear")

dim(smc_mask) <- c(na,nh,num_int)
smc_mask <- parameter()
dim(FOI_smc) <- c(na,nh,num_int)
FOI_smc[,,] <- FOI[i,j,k] * (1 - (smc_mask[i,j,k] * P_smc))

##Parameters relating to reduced infectivity of existing infections upon SMC dose
dim(rel_c_days) <- n_days + 1
rel_c_days <- parameter()

rel_c <- interpolate(days, rel_c_days, "constant")
dim(smc_rel_c_mask) <- c(na,nh,num_int)
smc_rel_c_mask[,,] <- 1 - (smc_mask[i,j,k] * (1-rel_c)) #Value = SMC_rel_c for those scheduled to receive SMC on rel_c_days. Else = 1.


################## ITN #######################
# See supplementary materials S2 from http://journals.plos.org/plosmedicine/article?id=10.1371/journal.pmed.1000324#s6
max_itn_cov <- parameter()


Q0 <- parameter()
phi_bednets <- parameter()

# General intervention model terminology:
# r - probability of trying to repeat feed after hitting ITN
# d - probability of dying after hitting ITN
# s - probability of successful feed after hitting ITN

dim(r_itn_daily) <- n_days + 1 #Bednet insecticide decay
dim(s_itn_daily) <- n_days + 1 #Proportion of individuals in the ITN compartment who currently have ITNs.= daily_ITN / max(daily_ITN)
r_itn_daily <- parameter()
s_itn_daily <- parameter()

r_itn <- interpolate(days, r_itn_daily, "linear")
s_itn <- interpolate(days, s_itn_daily, "linear")

################## GENERAL INTERVENTION PARAMETERS #######################
num_int <- parameter()

# cov is a vector of coverages for each intervention category:
dim(cov_) <- 4
cov_[1] <- (1-max_itn_cov)*(1-max_smc_cov)  # {No intervention}
cov_[2] <- max_itn_cov*(1-max_smc_cov) # 	   {ITN only}
cov_[3] <- (1-max_itn_cov)*max_smc_cov	#      {SMC only}
cov_[4] <- max_itn_cov*max_smc_cov #	   {Both itn and SMC}
cov[] <- cov_[i]
dim(cov) <- num_int

# probability that mosquito bites and survives for each intervention category
dim(w_) <- 4
w_[1] <- 1
w_[2] <- 1 - phi_bednets + phi_bednets*s_itn
w_[3] <- 1
w_[4] <- 1 - phi_bednets + phi_bednets*s_itn
w[] <- w_[i]
dim(w) <- num_int

# probability that mosquito is repelled during a single attempt for each int. cat.
dim(z_) <- 4
z_[1] <- 0
z_[2] <- phi_bednets*r_itn
z_[3] <- 0
z_[4] <- phi_bednets*r_itn
z[] <- z_[i]
dim(z) <- num_int

# Calculating Z (zbar) and W (wbar) - see Supplementary materials 2 for details
dim(zhi) <- num_int
dim(whi) <- num_int
zhi[1:num_int] <- cov[i]*z[i]
whi[1:num_int] <- cov[i]*w[i]
zh <- sum(zhi)
wh <- sum(whi)
# Z (zbar) - average probability of mosquito trying again during single feeding attempt
zbar <- Q0*zh
# W (wbar) - average probability of mosquito successfully feeding during single attempt
wbar <- 1 - Q0 + Q0*wh

# p1 is the updated p10 given that interventions are now in place:
p1 <- wbar*p10/(1-zbar*p10)
Q <- 1-(1-Q0)/wbar # updated anthropophagy given interventions
av <- fv*Q # biting rate on humans
dim(av_mosq) <- num_int
av_mosq[1:num_int] <- av*w[i]/wh # rate at which mosquitoes bite each int. cat.


##------------------------------------------------------------------------------
###################
## MODEL OUTPUTS ##
###################
##------------------------------------------------------------------------------
dim(clin_inc) <- c(na,nh,num_int)
initial(clin_inc[,,]) <- init_T[i,j,k]
update(clin_inc[,,]) <- ST_trans[i,j,k] + SD_trans[i,j,k] + AT_trans[i,j,k] + AD_trans[i,j,k] + UD_trans[i,j,k]  + UT_trans[i,j,k]

prev_dim <- parameter()
dim(min_age_prev) <- prev_dim
dim(max_age_prev) <- prev_dim
min_age_prev <- parameter(type = "integer")
max_age_prev <- parameter(type = "integer")

dim(n_prev) <- prev_dim

n_prev[1:prev_dim] <- sum(S[min_age_prev[i]:max_age_prev[i],,]) + sum(T[min_age_prev[i]:max_age_prev[i],,]) + sum(D[min_age_prev[i]:max_age_prev[i],,]) +
  sum(A[min_age_prev[i]:max_age_prev[i],,]) + sum(U[min_age_prev[i]:max_age_prev[i],,]) + sum(P[min_age_prev[i]:max_age_prev[i],,])

dim(detect_prev_full) <- c(na,nh,num_int)
detect_prev_full[,,] <- T[i,j,k] + D[i,j,k]  + Binomial(A[i,j,k],p_det[i,j,k])

dim(detect_prev) <- prev_dim
detect_prev[1:prev_dim] <- sum(detect_prev_full[min_age_prev[i]:max_age_prev[i],,])

##Proportion of population in user defined prevalence age groups
dim(n_ud_prev) <- prev_dim
initial(n_ud_prev[]) <- min_age_prev[i] #user defined prevalence
update(n_ud_prev[]) <- n_prev[i]

##Proportion of population in user defined prevalence age groups with detectable malaria
dim(n_ud_detect_prev) <- prev_dim
initial(n_ud_detect_prev[]) <- min_age_prev[i]
update(n_ud_detect_prev[]) <- detect_prev[i]

# #Incidence
inc_dim <- parameter()
dim(min_age_inc) <- inc_dim
dim(max_age_inc) <- inc_dim
min_age_inc <- parameter(type = "integer")
max_age_inc <- parameter(type = "integer")

dim(n_ud_inc) <- inc_dim
initial(n_ud_inc[]) <- min_age_inc[i]
update(n_ud_inc[]) <- sum(clin_inc[min_age_inc[i]:max_age_inc[i],,]) / dt

Sh <- sum(S[,,])
Th <- sum(T[,,])
Dh <- sum(D[,,])
Ah <- sum(A[,,])
Uh <- sum(U[,,])
Ph <- sum(P[,,])
H <- Sh + Th + Dh + Ah + Uh + Ph

initial(S_count) <- 0
initial(T_count) <- 0
initial(D_count) <- 0
initial(A_count) <- 0
initial(U_count) <- 0
initial(P_count) <- 0

update(S_count) <- Sh
update(T_count) <- Th
update(D_count) <- Dh
update(A_count) <- Ah
update(U_count) <- Uh
update(P_count) <- Ph

dim(all) <- c(na,nh,num_int)
all[,,] <- S[i,j,k] + T[i,j,k] + D[i,j,k] + A[i,j,k] + U[i,j,k] + P[i,j,k]
##-------------------- IMMUNITY OUTPUTS ----------------------------------
dim(icm_pop) <- c(na,nh,num_int)
icm_pop[,,] <- all[i,j,k] * ICM[i,j,k]
initial(icm_mean) <- 0
update(icm_mean) <- sum(icm_pop[,,]) / H

dim(ica_pop) <- c(na,nh,num_int)
ica_pop[,,] <- all[i,j,k] * ICA[i,j,k]
initial(ica_mean) <- 0
update(ica_mean) <- sum(ica_pop[,,]) / H

dim(id_pop) <- c(na,nh,num_int)
id_pop[,,] <- all[i,j,k] * ID[i,j,k]
initial(id_mean) <- 0
update(id_mean) <- sum(id_pop[,,]) / H

dim(ib_pop) <- c(na,nh,num_int)
ib_pop[,,] <- all[i,j,k] * IB[i,j,k]
initial(ib_mean) <- 0
update(ib_mean) <- sum(ib_pop[,,]) / H

dim(ic_pop) <- c(na,nh,num_int)
ic_pop[,,] <- all[i,j,k] * IC[i,j,k]
initial(ic_mean) <- 0
update(ic_mean) <- sum(ic_pop[,,]) / H

##-------------------- OTHER OUTPUTS ----------------------------------
dim(all_deaths) <- c(na,nh,num_int)
all_deaths[,,] <- S_death[i,j,k] + T_death[i,j,k] + D_death[i,j,k] + A_death[i,j,k] + U_death[i,j,k] + P_death[i,j,k]
initial(natural_deaths) <- 0
update(natural_deaths) <- sum(all_deaths[,,])

dim(epsilon_0) <- c(na,nh,num_int)
epsilon_0[,,] <- (all[i,j,k] * EIR[i,j,k]) / foi_age[i]
initial(EIR_mean) <- 0
update(EIR_mean) <- sum(epsilon_0[,,])

initial(mu_mosq) <- 0
update(mu_mosq) <- mu

##------------- FOR INTERACTION WITH MONTY. ALLOWS FITTING TO PREVALENCE DATA --------------
prevalence <- n_ud_detect_prev[1]/ n_ud_prev[1]
tests <- data()
positive <- data()
positive ~ Binomial(tests, prevalence)

##----------------------- FULL DETECT OUTPUT ----------------------------
dim(detect) <- c(na,nh,num_int)
initial(detect[,,]) <- init_D[i,j,k] + init_T[i,j,k]
update(detect[,,]) <- T[i,j,k] + D[i,j,k]  + A[i,j,k]*p_det[i,j,k]

dim(n) <- c(na,nh,num_int)
initial(n[,,]) <- init_A[i,j,k] + init_T[i,j,k] + init_D[i,j,k] + init_P[i,j,k] + init_U[i,j,k] + init_S[i,j,k]
update(n[,,]) <- all[i,j,k]
