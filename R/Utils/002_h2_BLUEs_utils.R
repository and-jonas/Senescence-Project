#HELPER FUNCTIONS ----

f_spats <- function(data, response = "value", random = as.formula("~Yf + Xf"),
                    genotype.as.random = TRUE, genotype = "Gen_Name") {
  SpATS(response = response, random = ~ Yf + Xf,
        spatial = ~PSANOVA(RangeBL, RowBL, nseg = c(20,20), nest.div=c(2,2)), 
        genotype = genotype, genotype.as.random = genotype.as.random, data = data,
        control = list(maxit = 100, tolerance = 1e-03, monitoring = 0))
}

make_design_matrix <- function(geno, names) {
  Nrow <- length(geno)
  Ncol <- length(names)
  col <- match(geno, names)
  frame <- data.frame(i = c(1:Nrow), j = col, v = rep(1,Nrow))
  frame <- subset(frame, is.na(col) == FALSE)
  L <- as.list(frame)
  X <- spam::spam(L, nrow = Nrow, ncol = Ncol)
  return(X)
}

construct_genotype_prediction_matrix <- function(object, newdata) {
  Z_geno = make_design_matrix(newdata[,object$model$geno$genotype], object$terms$geno$geno_names)
  if(object$model$geno$as.random)
    Z_geno <- Z_geno[,object$terms$geno$ndx]
  else
    Z_geno <- Z_geno[, object$terms$geno$ndx[2:length(object$terms$geno$ndx)]]
  Z_geno	
}

#extract BLUE from SpATS
extract_BLUE <- function(object){
  intercept <- object$coeff['Intercept']
  gen_mod_mat <- construct_genotype_prediction_matrix(object, object$data)
  gen_coeff <- object$coeff[1:ncol(gen_mod_mat)]
  BLUE <- as.data.frame(intercept + gen_coeff) %>% 
    rownames_to_column() %>% 
    as_tibble() %>% 
    rename(Gen_Name = 1, BLUE = 2)
}

#calculate heritability using BLUEs

get_h2_years <- function(data, fixed = fixed, random = "Gen_Name"){
  #check if random is factor and convert if required
  if(!is.factor(data$Gen_Name)){
    print("Gen_Name converted to factor")
    data$Gen_Name <- as.factor(data$Gen_Name)
  }
  #fit linear mixed model
  mod <- asreml(fixed = as.formula(paste("BLUE ~", fixed)),
                random = as.formula(paste("~", random)),
                data = data)
  #extract variance components
  GenVar <- summary(mod)$varcomp["Gen_Name", "component"]
  ErrVar <- summary(mod)$varcomp["units!R","component"]
  #calculate broad-sense heritability
  H2 <- round(GenVar/(GenVar + ErrVar/2), 2)
}

