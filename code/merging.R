strain_nums = c(1,6,54,57, 86)
#strain_nums = c(1)


createDF <- function(num, var){
    x <- read.csv(paste0("HM", num,"iut_counts"), header = TRUE)
    var <- x[,-1]
    row.names(var)<- x[,1]
}
i <- 1:8
print (i)
print (strain_nums[i])
f_list <- paste0("HM", strain_nums[i], "_X_counts")
print (f_list)
HM1 <- read.csv(f_list[1], header = TRUE)
hm1 <- HM1[,-1]
row.names(hm1)<- HM1[,1]

HM6 <- read.csv(f_list[2], header = TRUE)
hm6 <- HM6[,-1]
row.names(hm6)<- HM6[,1]

HM17 <- read.csv(f_list[3], header = TRUE)
hm17 <- HM17[,-1]
row.names(hm17)<- HM17[,1]



HM54 <- read.csv(f_list[4], header = TRUE)
hm54 <- HM54[,-1]
row.names(hm54)<- HM54[,1]

HM57 <- read.csv(f_list[5], header = TRUE)
hm57 <- HM57[,-1]
row.names(hm57)<- HM57[,1]

HM66 <- read.csv(f_list[6], header = TRUE)
hm66 <- HM66[,-1]
row.names(hm66)<- HM66[,1]

HM68 <- read.csv(f_list[7], header = TRUE)
hm68 <- HM68[,-1]
row.names(hm68)<- HM68[,1]

HM86 <- read.csv(f_list[8], header = TRUE)
hm86 <- HM86[,-1]
row.names(hm86)<- HM86[,1]

###merging multiple files

x <- list(hm17, hm54, hm57, hm66, hm68, hm86)
cmb <- merge(hm1, hm6, by = "row.names", all = TRUE)
cmb1 <- cmb[,-1]
row.names(cmb1)<- cmb[,1]

for (m in x){
    print (m)
    cmb <- merge(cmb1, x, by = "row.names", all = TRUE)
    cmb1 <- cmb[,-1]
    row.names(cmb1)<- cmb[,1]
}


#cmb <- Reduce(function(x,y) merge(x,y, by = "row.names", all = TRUE), list(hm1, hm6, hm54, hm57, hm86))
#merge(hm1, hm6, hm54, hm57, hm86, by = "row.names", all = TRUE)


cmb<- rbind(hm1, hm6, hm17)


