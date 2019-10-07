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
#'@export linreg

linreg<-setRefClass("linreg", fields = list(formula="formula",data="data.frame",
                                            beta_cap ="matrix",y = "numeric",
                                            model_mat ="matrix",trans_my_mat ="matrix",
                                            my_matrix ="matrix", fitted_val ="matrix",
                                            resi_val ="matrix",df ="numeric",
                                            residual_variance ="numeric",variance_reg_coef ="matrix",
                                            t_val ="matrix", p_val ="matrix", f_formula = "character",
                                            data_name="character"),
                    methods = list(
                    initialize=function(formula=as.formula(),data=as.data.frame(),...){
                      check<-all.vars(formula)
                      f_formula<<-Reduce(paste, deparse(formula))
                      data_name<<-deparse(substitute(data))
                      stopifnot(is.data.frame(data)| !is.null(check)|length(check)!=0)
                      
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
                      
                      residual_variance<<-as.numeric((t(resi_val)%*%resi_val)/df)
                      
                      #----------- Variance of regression coefficient variance_reg_coef -----------------
                      
                      #variance_reg_coef<<- diag(c(residual_variance)*(solve(t(my_matrix)%*%my_matrix)))
                      variance_reg_coef<<- residual_variance*solve(t(my_matrix)%*%my_matrix)
                     
                      #----------- t-value : t_val -----------------
                      
                      t_val<<-beta_cap / as.double(sqrt(diag(variance_reg_coef)))
                      
                      #----------- p-value : p_val -----------------
                      
                      p_val<<- 2*pt(-abs(t_val),df)
                      
                      return(.self)
                    },
                    
                    print=function(){
                     cat("Call:\n")
                     text_str<-paste("linreg(formula = ",f_formula,", data = ",data_name,")",sep = "")
                     cat(text_str,"\n")
                     cat("Coefficients: \n\n")
                     disp<-c(round(beta_cap,digits = 2))
                     #names(disp)<-row.names(beta_cap)
                     cat(row.names(beta_cap),"\n")
                     cat(disp)
                    },
                    
                    plot=function(){
                      library(ggplot2)
                      #install.packages("gridExtra")
                      #library(gridExtra)
                      liu_blue <- "#54D8E0"
                      theme_liu <- theme(plot.margin = unit(c(1,1,1,1), "cm"), 
                                         panel.background = element_rect(fill="white"),
                                         panel.grid.major.y = element_blank(),
                                         panel.grid.minor.y = element_blank(),
                                         panel.grid.major.x = element_blank(),
                                         panel.grid.minor.x = element_blank(),
                                         axis.line = element_line(color= "#58585b", size=0.1),
                                         axis.text.x = element_text(color="Black", size="10"),
                                         axis.text.y = element_text(color="Black", size="10"),
                                         axis.title.x = element_text(color="Black", size="10", face="bold"),
                                         axis.title.y = element_text(color="Black", size="10", face="bold"),
                                         axis.ticks.y = element_blank(),
                                         axis.ticks.x = element_line(color = "#58585b", size = 0.3),
                                         plot.title = element_text(color="Black", face="bold", size="14", hjust = 0.5),
                                         legend.position="bottom", legend.title = element_blank(),
                                         legend.key = element_blank(),
                                         legend.text = element_text(color="Black", size="10"))
                     
                      plot1<-ggplot(as.data.frame(fitted_val,resi_val),aes(fitted_val,resi_val))+
                        geom_point(colour = liu_blue) +
                        geom_smooth(color = "red") +
                        geom_abline(slope = 0,
                                    intercept = 0,
                                    linetype = "dotted") +
                        ggtitle("Residual vs Fitted") +
                        ylab("Residuals") +
                        xlab("Fitted Values") +
                        theme_liu
                      
                      #std_res<<-sqrt(abs((resi_val - mean(resi_val)) / sqrt(residual_variance)))
                      std_res<-sqrt(abs(resi_val/sd(resi_val)))
                              
                      plot2<-ggplot(as.data.frame(fitted_val,resi_val),aes(fitted_val,std_res))+
                        geom_point(colour = liu_blue) +
                        geom_smooth(color = "red") +
                        geom_abline(slope = 0,
                                    intercept = 0,
                                    linetype = "dotted") +
                        ggtitle("Scale-Location") +
                        ylab("|Standardized Residuals|") +
                        xlab("Fitted Values") +
                        theme_liu
                      # p1<-qplot(fitted_val,resi_val,main="Residuals vs Fitted", xlab="Fitted Values",
                      #         ylab="Residuals")
                      #p1<-p1+my_theme
                      #std_res<<-sqrt(abs(resi_val/sd(resi_val)))
                      #p2<-qplot(fitted_val,std_res,main="Scale-Location", xlab="Fitted Values",
                      #      ylab="|Standardized Residual|")
                      #p2<-p2+my_theme
                      figure<-grid.arrange(plot1,plot2,ncol=2,nrow=1)
                      figure
                     
                    },
                    
                    resid=function(){
                      #cat("Residual Vector is : \n")
                      return(as.vector(resi_val))
                    },
                    
                    pred=function(){
                      #cat("predicted values are : \n")
                      return(round(fitted_val,digits = 2))
                    },
                   
                    coef=function(){
                      #cat("Coefficients: \n")
                      disp<-c(round(beta_cap,digits = 2))
                      names(disp)<-row.names(beta_cap)
                      return(disp)
                    },
                    
                    summary=function(){
                      #cat("Coefficients: \n")
                     # disp<-1
                      #disp<-c(round(beta_cap,digits = 2))
                      #names(disp)<-row.names(beta_cap)
                      #disp;
                   #   for(i in 1:length(beta_cap))
                    #  {
                     #   text_str<-c(row.names(beta_cap)[i],round(beta_cap[i],2))
                      #  cat(text_str,"\n")#,quote = FALSE,row.names=FALSE)
                      #}
                      #--------------
                      coef_mx <- data.frame(
                             #var = rownames(beta_cap),
                             estimate = round(beta_cap, 2),
                             std.error = round(sqrt(diag(variance_reg_coef)), 2),
                             t_value = round(t_val,2),
                             p_value = round(p_val,4)
                             #sigma_cap<-residual_variance,
                             #deg_fre<-df
                         )
                      for(i in 1:nrow(coef_mx))
                      {
                        if(coef_mx$p_value[i]==0){
                          coef_mx$p_value[i] = "***"
                        }else if(coef_mx$p_value[i] > 0 & coef_mx$p_value[i]<= 0.001){
                          coef_mx$p_value[i] <- "**"
                        } else if(coef_mx$p_value[i] > 0.001 & coef_mx$p_value[i] <= 0.01){
                          coef_mx$p_value[i] <- "*"
                        } else if(coef_mx$p_value[i] > 0.01 & coef_mx$p_value[i] <= 0.05){
                          coef_mx$p_value[i] <- "."
                        } else if(coef_mx$p_value[i] > 0.05 & coef_mx$p_value[i] <= 0.1){
                          coef_mx$p_value[i] <- " "
                        } else if(coef_mx$p_value[i] > 0.1){
                          coef_mx$p_value[i] <- " "
                        }
                      }
                      names(coef_mx)<-NULL
                      #cat(as.vector(coef_mx[1,]))
                      #my_str<-paste(c(coef_mx[1,]),sep = "")
                      #my_print(my_str)
                      my_print(coef_mx)
                      #cat("\n")
                      cat("Residual standard error:",round(sqrt(diag(variance_reg_coef)[1]), 1),"on",df,"degrees of freedom")
                      #my_print(my_str)
                      #--------------
                      
                      #cat("\n\nCoefficients standard error: \n")
                      
                      #cat("\n\nT Value: \n")
                      #print.default(t_val)
                      #cat("\n\nP Value: \n")
                      #print.default(p_val)
                      #cat("\n\nEstimate of Sigma: \n")
                      #print.default(variance_reg_coef)
                      #cat("\n\nDegree of Freedom: \n")
                      #print.default(df)
                    }
                    
                    
                    ))


my_print=function(x){
  
  print(x)
  
}

