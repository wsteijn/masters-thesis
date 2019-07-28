library(rsm)
library(boot)

#function for inputting data table or matrix and making it into an output that can fit into the format for making coded data
#data is the data frame that contains the data that will be changed into coded data for the RSM
#x1 is the name of the column in the data frame that will be x1 in the RSM
#x2 is the name of the column in the data frame that will be x2 in the RSM
#x1 and x2 variables must be put in quotes
RSM.function = function(data = "", x1 = "", x2 = ""){
  #define x variables as the columns in the data frame that have column names x1 and x2 where x1 and x2 are the defined in the arguments of the function
  x_variables = data[,c(x1,x2)]
  #define the y variables as all the collumns in the data frame other than the x_variable columns
  y_variables = data[,!names(data) %in% c(x1,x2)] 
  #create a data frame with the x_variables as the first columns, and the y_variables as the rest of the columns
  RSM.matrix = as.data.frame(cbind(x_variables,y_variables))
  return(RSM.matrix)
}
RSM.matrix = RSM.function(data = QMP_final_matrix_subset, x1 = "Carbt1", x2 = "Proteint1")
RSM.matrix
  
rsm(Akkermansia ~FO(Carbt1,Proteint1), data = RSM.matrix)
summary(rsm(Akkermansia ~FO(Carbt1,Proteint1), data = RSM.matrix))
#create the coded data dataframe
middle.values = NULL
range.values = NULL
for(i in colnames(RSM.matrix)[1:2]){
  middle.values[i] = apply(RSM.matrix[i],2,function(x)(min(x) + max(x))/2)
  range.values[i] = apply(RSM.matrix[i],2, function(x)(max(x) - min(x))/2)
}
middle.values
range.values
RSM.matrix = coded.data(RSM.matrix, x1 ~ (Carbt1 - middle.values[[1]])/(range.values[[1]]), x2 ~ (Proteint1 - middle.values[[2]])/(range.values[[2]]))
print(RSM.matrix, decode = FALSE)

#############################################################################################################
#Coded data problem
#RSM.matrix = coded.data(RSM.matrix, x1 ~ (Carbt1 - middle.values[[1]])/(range.values[[1]]), x2 ~ (Proteint1 - middle.values[[2]])/(range.values[[2]]))
#the above line of code works correctly, and creates the coded data, because "Carbt1" and "Proteint1" are explicitly named

#RSM.matrix = coded.data(RSM.matrix, x1 ~ (i - middle.values[[1]])/(range.values[[1]]), x2 ~ (i2 - middle.values[[2]])/(range.values[[2]]))
#the above line of code does not work correctly - it runs, but values for Carbt1 and Proteint1 are not changed to coded data

#to check if the coded data is created, each column, x1 and x2, should have a range from -1 to 1
#so the problem is getting the coded.data command to accept variables for column names of x1 and x2 so that the coded.data command can work as part of the pipeline,
#and the user does not have to manually type in the column names of the x1 and x2 variables every time they need to make coded data
###############################################################################################################




#An important aspect of response-surface analysis is using an appropriate coding transformation
#of the data. The way the data are coded affects the results of canonical analysis (see
# Section 4) and steepest-ascent analysis (see Section 6); for example, unless the scaling factors
#are all equal, the path of steepest ascent obtained by fitting a model to the raw predictor values
#will differ from the path obtained in the coded units, decoded to the original scale. Using
#a coding method that makes all coded variables in the experiment vary over the same range
#is a way of giving each predictor an equal share in potentially determining the steepest-ascent
#path. Thus, coding is an important step in response-surface analysis.




#define data frame created by the loop
RSM.output.matrix = NULL

#start with the case where there are two x variables and then multiple y variables and each y variable will be modeled by a FO and SO rsm, 
#and based off certain criteria (p-value, adjusted r-squared, AIC, and 10-fold cross validation), one model will be chosen as the best. If neither the FO or SO 
#model fits the criteria, then nothing is added to the RSM.output.matrix for that 'y-variable' iteration. If one of the three criteria is met, the y variable 
#(bacteria species), which order of model was best, p-value, and adjusted r-squared values for that iteration are added to the RSM.output.matrix and the graph 
#of this 'best' RSM is saved to the specified path as a pdf.

#First two columns are x1 and x2 - start iterating over the 'y-variables' starting with column 3
#i is the column name of the 3rd to the final column - these are the 'y-variables'
for(i in colnames(RSM.matrix)[3:length(RSM.matrix)]){
  #make the first order response surface model
  FO_response_surface_model= rsm(RSM.matrix[[i]] ~ FO(x1,x2), data = RSM.matrix)
  #make the second order model
  SO_response_surface_model = rsm(RSM.matrix[[i]] ~ SO(x1,x2), data = RSM.matrix)
  #define a variable for the summary of the first order and second order models
  FO_summary_response_surface = summary(FO_response_surface_model)
  SO_summary_response_surface = summary(SO_response_surface_model)
  #extract the p-values for the FO and SO model so they can be compared
  FO.p.value = FO_summary_response_surface$lof$`Pr(>F)`[1]
  SO.p.value = pf(SO_summary_response_surface$fstatistic[[1]], SO_summary_response_surface$fstatistic[[2]],SO_summary_response_surface$fstatistic[[3]],lower.tail = FALSE)
  
  #first case: the FO p-value is below a certain threshold, while the p-value for the SO model is not -- the FO model is significant, while the SO model is not
  #in this case, the FO model is clearly the best
  if(FO.p.value<.3 & SO.p.value > .3){
    #extract the adjusted r-squared value from the summary of the model
    adjusted.r.squared = FO_summary_response_surface$adj.r.squared
    #extract the x1 and x2 coefficients of the direction of steepest ascent, and then code them back into the original units
    x1_SA_coded = FO_summary_response_surface$sa[[1]]
    x2_SA_coded = FO_summary_response_surface$sa[[2]]
    a = code2val(data.frame(x1=FO_summary_response_surface$sa[[1]],x2 = FO_summary_response_surface$sa[[2]]),codings(RSM.matrix))
    x1_SA_noncoded = a[[1]] - middle.values[[1]]
    x2_SA_noncoded = a[[2]] - middle.values[[2]]
    #set the outputs that would be for Second order models as "NA"
    x1_SP_noncoded = NA
    x2_SP_noncoded = NA
    SO_eigen = NA
    #set the output as a row to be added to the output matrix
    output = (c(i,"FO",FO.p.value,adjusted.r.squared,x1_SA_coded,x2_SA_coded, x1_SA_noncoded, x2_SA_noncoded, x1_SP_noncoded, x2_SP_noncoded, SO_eigen))
    #output the species of bacteria (the column name), the order of the model (FO),  the p-value, and the adjusted r-squared value
    output = (c(i,"FO",FO.p.value,adjusted.r.squared,x1_SA_coded,x2_SA_coded, x1_SA_noncoded, x2_SA_noncoded))
    #specify a path and file name for the pdf of the graph of the FO model
    mypath = file.path("C:/Users/Will/Documents/thesis stuff/RSM Pipeline", paste("FO_RSM_", i, ".pdf", sep = ""))
    pdf(file = mypath)
    contour(FO_response_surface_model, ~x1  + x2 , image = TRUE)
    title(paste("FO",i, sep = ""))
    dev.off()
    #print the output to the console and save the output to the previously defined matrix
    print(output)
    RSM.output.matrix = as.matrix(rbind(RSM.output.matrix, output))
  }
  #second case: the SO p-value is below a certain threshold, while the p-value for the FO model is not -- the SO model is significant, while the FO model is not
  #in this case, the SO model is clearly the best
  if(SO.p.value<.3 & FO.p.value > .3){
    #extract the adjusted r-squared value from the summary of the model
    adjusted.r.squared = SO_summary_response_surface$adj.r.squared
    #extract the x1 and x2 coefficients of the stationary point, and then code them back into the original units
    a = code2val(data.frame(x1=SO_summary_response_surface$canonical$xs[[1]],x2 = SO_summary_response_surface$canonical$xs[[2]]),codings(RSM.matrix))
    x1_SP_noncoded = (a[[1]])
    x2_SP_noncoded = (a[[2]])
    #extract the eigenvectors so it is known if the stationary point is a min, max, or saddle point
    SO_eigen = SO_summary_response_surface$canonical$eigen[[1]]
    #set the FO-specific output variables to NA
    x1_SA_coded = NA
    x2_SA_coded = NA
    x1_SA_noncoded = NA
    x2_SA_noncoded = NA
    #output the species of bacteria (the column name), the order of the model (SO),  the p-value, and the adjusted r-squared value, and other model information 
    output = (c(i,"SO",FO.p.value,adjusted.r.squared,x1_SA_coded,x2_SA_coded, x1_SA_noncoded, x2_SA_noncoded, x1_SP_noncoded, x2_SP_noncoded, SO_eigen))
    #specify a path and file name for the pdf of the graph of the SO model
    mypath = file.path("C:/Users/Will/Documents/thesis stuff/RSM Pipeline", paste("SO_RSM_", i, ".pdf", sep = ""))
    pdf(file = mypath)
    contour(SO_response_surface_model, ~x1  + x2 , image = TRUE,at = summary(SO_response_surface_model)$canonical$xs)
    title(paste("SO",i, sep = ""))
    dev.off()
    #print the output to the console and save the output to the previously defined matrix
    print(output)
    RSM.output.matrix = as.matrix(rbind(RSM.output.matrix, output))
  }
  
  #third case: both the p-values for the FO and SO models are below a certain threshold
  #in this case, both the SO and FO model are significant, so not clear yet which one is the best
  if(FO.p.value<.3 & SO.p.value < .3){
    #to compare the models, consider adjusted r-squared of the rsm and AIC and cross validation error of the generalized linear model which has
    #the same parameters as the linear model used to create the response surface model
    #extract the adjusted r squared values from the FO and SO model summaries
    FO.adjusted.r.squared = FO_summary_response_surface$adj.r.squared
    SO.adjusted.r.squared = SO_summary_response_surface$adj.r.squared
    #make generalized linear models for the FO and SO case
    formula.FO = formula(bquote(.(as.name(i)) ~ x1 + x2))
    formula.SO = formula(bquote(.(as.name(i)) ~ x1+x2+I(x1*x2) + I(x1*x1) + I(x2*x2)))
    FO.glm = glm(formula.FO, data = RSM.matrix)
    SO.glm = glm(formula.SO, data = RSM.matrix)
    #extract the AIC values for the FO and SO generalized linear models
    FO.AIC = FO.glm$aic
    SO.AIC = SO.glm$aic
    #find the cross validation error from a 10 fold cross validation 
    #delta is a vector of two, where the second component is the adjusted cross validation estimate - the adjustment is designed to compensate for the bias 
    #introduced by not using leave-one-out cross-validation
    FO.cv.error = (cv.glm(RSM.matrix,FO.glm,K=10))$delta[1]
    SO.cv.error = (cv.glm(RSM.matrix,SO.glm,K=10))$delta[1]
    
      #first case - adjusted r-squared for SO is larger than the adjusted r-squared for the FO model AND SO AIC is smaller than FO AIC AND SO cross validation error is smaller than FO cv error
      #SO model is better by these three criteria - SO model is chosen as best and it is treated like the second case described above
      if(SO.adjusted.r.squared>FO.adjusted.r.squared & SO.AIC<FO.AIC &SO.cv.error<FO.cv.error ){ #& SO.cv.error<FO.cv.error
        adjusted.r.squared = SO_summary_response_surface$adj.r.squared
        #extract the x1 and x2 coefficients of the stationary point, and then code them back into the original units
        a = code2val(data.frame(x1=SO_summary_response_surface$canonical$xs[[1]],x2 = SO_summary_response_surface$canonical$xs[[2]]),codings(RSM.matrix))
        x1_SP_noncoded = (a[[1]])
        x2_SP_noncoded = (a[[2]])
        #extract the eigenvectors so it is known if the stationary point is a min, max, or saddle point
        SO_eigen = SO_summary_response_surface$canonical$eigen[[1]]
        #set the FO-specific output variables to NA
        x1_SA_coded = NA
        x2_SA_coded = NA
        x1_SA_noncoded = NA
        x2_SA_noncoded = NA
        #output the species of bacteria (the column name), the order of the model (SO),  the p-value, and the adjusted r-squared value, and other model information 
        output = (c(i,"SO",FO.p.value,adjusted.r.squared,x1_SA_coded,x2_SA_coded, x1_SA_noncoded, x2_SA_noncoded, x1_SP_noncoded, x2_SP_noncoded, SO_eigen))
        mypath = file.path("C:/Users/Will/Documents/thesis stuff/RSM Pipeline", paste("SO_RSM_", i, ".pdf", sep = ""))
        pdf(file = mypath)
        contour(SO_response_surface_model, ~x1  + x2 , image = TRUE,at = summary(SO_response_surface_model)$canonical$xs)
        title(paste("SO",i, sep = ""))
        dev.off()
        print(output)
        RSM.output.matrix = as.matrix(rbind(RSM.output.matrix, output))
      }
    
      #if SO model is not better by all three criteria, FO model is chosen as best due to simplicity 
      #FO model is best - treated like the first case as decribed above
      else{
        adjusted.r.squared = FO_summary_response_surface$adj.r.squared
        x1_SA_coded = FO_summary_response_surface$sa[[1]]
        x2_SA_coded = FO_summary_response_surface$sa[[2]]
        a = code2val(data.frame(x1=FO_summary_response_surface$sa[[1]],x2 = FO_summary_response_surface$sa[[2]]),codings(RSM.matrix))
        x1_SA_noncoded = a[[1]]- middle.values[[1]]
        x2_SA_noncoded = a[[2]] - middle.values[[2]]
        x1_SP_noncoded = NA
        x2_SP_noncoded = NA
        SO_eigen = NA
        output = (c(i,"FO",FO.p.value,adjusted.r.squared,x1_SA_coded,x2_SA_coded, x1_SA_noncoded, x2_SA_noncoded, x1_SP_noncoded, x2_SP_noncoded, SO_eigen))
        mypath = file.path("C:/Users/Will/Documents/thesis stuff/RSM Pipeline", paste("FO_RSM_", i, ".pdf", sep = ""))
        pdf(file = mypath)
        contour(FO_response_surface_model, ~x1  + x2 , image = TRUE)
        title(paste("FO",i, sep = ""))
        dev.off()
        print(output)
        RSM.output.matrix = as.matrix(rbind(RSM.output.matrix, output))
    }
  }
  #define the column names for the output matrix built by the for loop
  #if the model type is FO: columns 5 and 6 are the values for x1 and x2 in coded data of the direction of steepest ascent of a radius of distance 1 (in coded units) from the center of the data, and 
  #                         columns 7 and 8 are the values for x1 and x2 in non-coded of the direction of steepest ascent 
  #                         - because the data in columns 5 and 6 is coded and thus has the same range of variance, this can tell us which variable, x1 or x2, has more of an influence on the maximum
  #                 
  #if the model type is SO: columns 9 and 10 are the values for x1 and x2 in coded data of the direction of steepest ascent, and 
  #                         the last two columns are the eigenvalues - if both are negative -> the stationary point is a max
  #                                                                 - if one is negative and one is positive -> the stationary point is a saddle point
  #                                                                 - if both are positive -> the stationary point is a minimum  
  #
  
  colnames(RSM.output.matrix) = c("bacteria species","Model type","p-value","adj R-squared","x1 S.A. coded","x2 coded S.A.", "x1 S.A. non-coded", "x2 S.A. non-coded", "SO x1 stationary point", "SO x2 stationary point", "SO 1st eigenvalue", "SO 2nd eigenvalue")
  #write the output matrix to a csv file
  write.csv(RSM.output.matrix, file = "C:/Users/Will/Documents/thesis stuff/RSM Pipeline/RSM_Output.csv")
}
