soilStrength <- 
function(clay.content, matric.suction = NULL, water.content = NULL)
{
   if (!is.null(matric.suction) || !is.null(water.content)) {
        if (is.null(clay.content)) 
            warning("To estimate soil strength, please inform water.content or matric.suction")

   if (is.numeric(matric.suction) & is.numeric(water.content)) {
            warning("To estimate soil strength, please inform only one of them: water.content or matric.suction")
   }

   pre.cons.water <- function(clay.content, water.content) {
       mh <- c()
       for (j in 1:length(clay.content)) {
               if (clay.content[j] <= 20) {
                  mh[j] <- ((((0.42)/(water.content[j] - 0.049356))^(1/0.42)) - 
                    1)^(1/1.72) * (1/0.79)
                }
                else if (clay.content[j] > 20 & clay.content[j] <= 
                  31) {
                  mh[j] <- ((((0.45)/(water.content[j] - 0.08689))^(1/0.36)) - 
                    1)^(1/1.56) * (1/0.72)
                }
                else if (clay.content[j] > 31 & clay.content[j] <= 
                  37) {
                  mh[j] <- ((((0.46)/(water.content[j] - 0.10696))^(1/0.34)) - 
                    1)^(1/1.52) * (1/1.66)
                }
                else if (clay.content[j] > 37 & clay.content[j] <= 
                  52) {
                  mh[j] <- ((((0.5)/(water.content[j] - 0.125941))^(1/0.33)) - 
                    1)^(1/1.47) * (1/2.04)
                }
                else {
                  mh[j] <- ((((0.51)/(water.content[j] - 0.139358))^(1/0.28)) - 
                    1)^(1/1.38) * (1/2.27)
                }
          }
       return(round(mh, 0))
   }
   if (length(matric.suction) > 0) {
            matric.suction <- matric.suction
   }
      else {
            matric.suction <- pre.cons.water(clay.content = clay.content, 
                water.content = water.content)
   }
   pcs <- c()
   for (j in 1:length(clay.content)) {
            if (clay.content[j] < 20) {
                pcs[j] <- round(129 * matric.suction[j]^(0.15),0)
            }
            else if (clay.content[j] >= 20 & clay.content[j] <= 
                31) {
                pcs[j] <- round(123.3 * matric.suction[j]^(0.13),0)
            }
            else if (clay.content[j] > 31 & clay.content[j] <= 
                37) {
                pcs[j] <- round(85 * matric.suction[j]^(0.17),0)
            }
            else if (clay.content[j] > 37 & clay.content[j] <= 
                52) {
                pcs[j] <- round(70.1 * matric.suction[j]^(0.16),0)
            }
            else if (clay.content[j] > 52) {
                pcs[j] <- round(62.7 * matric.suction[j]^(0.15),0)
            }
     }
     pcs05 <- pcs*0.5
     pcs11 <- pcs*1.1
     soil.strength <- data.frame(pcs,pcs05,pcs11)
     colnames(soil.strength) <- c("Pc","LL.Pc","UL.Pc")
     return(soil.strength)
     }
}
