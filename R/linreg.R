#Linear regression Class implementation using RC Class

#' this class is calculating Regressions coefficients, fitted values, residuals, 
#' degrees of freedom, residual variance, variance of the regression coefficients,
#' t-values for each coefficient and p-value for each coefficient.
#'
#' @param formula this is the formula that need to be passed into contructor of the class, it will contain variable names from data frame.
#' @param data this will be the data frame that we are working on.
#' @return this will return class object, which will then be used to access multiple function of the class.
#' @examples
#' \dontrun{
#' data("iris")
#'mat_obj<-linreg$new()
#'mat_obj<-mat_obj$linreg(formula=Petal.Length~Species,data = iris)
#'mat_obj$print()
#'mat_obj$plot()
#'mat_obj$resid()
#'mat_obj$summary()

#'}

linreg<-setRefClass("linreg", fields = list(formula="formula",data="data.frame",
                                            beta_cap ="matrix",y = "numeric",
                                            model_mat ="matrix",trans_my_mat ="matrix",
                                            my_matrix ="matrix", fitted_val ="matrix",
                                            resi_val ="matrix",df ="numeric",
                                            residual_variance ="matrix",variance_reg_coef ="numeric",
                                            t_val ="matrix", p_val ="matrix"),
                    methods = list(
                    linreg=function(formula=as.formula(),data=as.data.frame(),...){
                      #  stopifnot(length(formula())==3)
                      y<<-data[,all.vars(formula)[1]]
                      my_matrix<<-model.matrix(formula,data)
                      my_matrix<<-cbind(my_matrix,y)

                      #----------- Regression Corfficient Beta Cap -----------------
                      
                      y<<-my_matrix[,ncol(my_matrix)]
                      my_matrix<<-my_matrix[,-ncol(my_matrix)]
                      trans_my_mat<<-t(my_matrix)
                      beta_cap<<-solve((trans_my_mat%*%my_matrix))%*%trans_my_mat%*%y
                      
                      #----------- The fitted values fitted_val -----------------
                      
                      fitted_val<<-my_matrix%*%beta_cap
                      
                      #----------- The residuals: resi_val -----------------
                      
                      resi_val<<- y - fitted_val
                      
                      #----------- degree of freedom: df_val -----------------
                      
                      p<-ncol(my_matrix)
                      n<-nrow(my_matrix)
                      df<<-n-p
                      
                      #----------- residual Variance  residual_variace -----------------
                      
                      residual_variance<<-(t(resi_val)%*%resi_val)/df
                      
                      #----------- Variance of regression coefficient variance_reg_coef -----------------
                      
                      variance_reg_coef<<- diag(c(residual_variance)*(solve(t(my_matrix)%*%my_matrix)))
                     
                      #----------- t-value : t_val -----------------
                      
                      t_val<<-beta_cap %/% sqrt(variance_reg_coef)
                      
                      #----------- p-value : p_val -----------------
                      
                      p_val<<- 2*pt(-abs(t_val),df)
                      
                      return(.self)
                    },
                    
                    print=function(){
                     cat("Call:\n")
                     cat("linreg(formula = Petal.Length ~ Species, data = iris)\n\n")
                     cat("Coefficients: \n")
                     disp<-c(beta_cap)
                     names(disp)<-row.names(beta_cap)
                     return(disp)
                    },
                    
                    plot=function(){
                      library(ggplot2)
                      install.packages("gridExtra")
                      library(gridExtra)
                      my_theme<-theme(axis.title=element_text(face="bold.italic",size="12", color="brown"), 
                                      legend.position="top",plot.background = element_rect(fill = "black"),
                                      panel.grid.major = element_line(colour = "black"),
                                      panel.grid.minor.y = element_blank(),
                                      plot.title = element_text(colour = "white",size = 14,face = "bold.italic",
                                                                hjust = .5))
                      p1<-qplot(fitted_val,resi_val,main="Residuals vs Fitted", xlab="Fitted Values",
                               ylab="Residuals")
                      p1<-p1+my_theme
                      std_res<-sqrt(abs(resi_val/sd(resi_val)))
                      p2<-qplot(fitted_val,std_res,main="Scale-Location", xlab="Fitted Values",
                            ylab="|Standardized Residual|")
                      p2<-p2+my_theme
                      figure<-grid.arrange(p1,p1,ncol=1,nrow=2)
                      figure
                     
                    },
                    
                    resid=function(){
                      cat("Residual Vector is : \n")
                      return(as.vector(resi_val))
                    },
                    
                    pred=function(){
                      cat("predicted values are : \n")
                      return(fitted_val)
                    },
                   
                    coef=function(){
                      cat("Coefficients: \n")
                      disp<-c(beta_cap)
                      names(disp)<-row.names(beta_cap)
                      return(disp)
                    },
                    
                    summary=function(){
                      cat("Coefficients: \n")
                      disp<-c(beta_cap)
                      names(disp)<-row.names(beta_cap)
                      cat(as.vector(disp))
                      cat("\n\nCoefficients standard error: \n")
                      
                      cat("\n\nT Value: \n")
                      cat(t_val)
                      cat("\n\nP Value: \n")
                      cat(p_val)
                      cat("\n\nEstimate of Sigma: \n")
                      cat(variance_reg_coef)
                      cat("\n\nDegree of Freedom: \n")
                      cat(df)
                    }
                    
                    
                    ))


