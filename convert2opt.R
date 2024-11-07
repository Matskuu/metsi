rm(list=ls())

setwd("C:/MyTemp/code/UTOPIA/metsi/sim_outdir/P6111A-tree-vol")

count <- 1

data <- read.table("data.xda", header=F)

# index 8 is where stock values and after harvest values start
# index 32 is where they end
# index 33- are stock volumes of different tree species:
#   - 33-36 are pine volumes from years 0, 5, 10 and 20
#   - 37-40 are spruce volumes from years 0, 5, 10 and 20
#   - 41-44 are other tree's volumes from years 0, 5, 10 and 20
data[,8:184] <- data[,2] * data[,8:184] # before data[,8:32] then [,8:44]

names(data) <- c("identifier",
"area",
"npv_1_percent",
"npv_2_percent",
"npv_3_percent",
"npv_4_percent",
"npv_5_percent",
"stock_0",
"stock_5",
"stock_10",
"stock_20",
"harvest_5",
"harvest_10",
"harvest_20",
"harvest_value_5",
"harvest_value_10",
"harvest_value_20",
"harvest_value_first_5",
"harvest_value_below_5",
"harvest_value_above_5",
"harvest_value_even_5",
"harvest_value_clearcut_5",
"harvest_value_first_10",
"harvest_value_below_10",
"harvest_value_above_10",
"harvest_value_even_10",
"harvest_value_clearcut_10",
"harvest_value_first_20",
"harvest_value_below_20",
"harvest_value_above_20",
"harvest_value_even_20",
"harvest_value_clearcut_20",
"stock_1_0",
"stock_1_5",
"stock_1_10",
"stock_1_20",
"stock_2_0",
"stock_2_5",
"stock_2_10",
"stock_2_20",
"stock_3_0",
"stock_3_5",
"stock_3_10",
"stock_3_20",
"stock_4_0",
"stock_4_5",
"stock_4_10",
"stock_4_20",
"stock_5_0",
"stock_5_5",
"stock_5_10",
"stock_5_20",
"stock_6_0",
"stock_6_5",
"stock_6_10",
"stock_6_20",
"stock_7_0",
"stock_7_5",
"stock_7_10",
"stock_7_20",
"stock_8_0",
"stock_8_5",
"stock_8_10",
"stock_8_20",
"stock_9_0",
"stock_9_5",
"stock_9_10",
"stock_9_20",
"stock_10_0",
"stock_10_5",
"stock_10_10",
"stock_10_20",
"stock_11_0",
"stock_11_5",
"stock_11_10",
"stock_11_20",
"stock_12_0",
"stock_12_5",
"stock_12_10",
"stock_12_20",
"stock_13_0",
"stock_13_5",
"stock_13_10",
"stock_13_20",
"stock_14_0",
"stock_14_5",
"stock_14_10",
"stock_14_20",
"stock_15_0",
"stock_15_5",
"stock_15_10",
"stock_15_20",
"stock_16_0",
"stock_16_5",
"stock_16_10",
"stock_16_20",
"stock_17_0",
"stock_17_5",
"stock_17_10",
"stock_17_20",
"stock_18_0",
"stock_18_5",
"stock_18_10",
"stock_18_20",
"stock_19_0",
"stock_19_5",
"stock_19_10",
"stock_19_20",
"stock_20_0",
"stock_20_5",
"stock_20_10",
"stock_20_20",
"stock_21_0",
"stock_21_5",
"stock_21_10",
"stock_21_20",
"stock_22_0",
"stock_22_5",
"stock_22_10",
"stock_22_20",
"stock_23_0",
"stock_23_5",
"stock_23_10",
"stock_23_20",
"stock_24_0",
"stock_24_5",
"stock_24_10",
"stock_24_20",
"stock_25_0",
"stock_25_5",
"stock_25_10",
"stock_25_20",
"stock_26_0",
"stock_26_5",
"stock_26_10",
"stock_26_20",
"stock_27_0",
"stock_27_5",
"stock_27_10",
"stock_27_20",
"stock_28_0",
"stock_28_5",
"stock_28_10",
"stock_28_20",
"stock_29_0",
"stock_29_5",
"stock_29_10",
"stock_29_20",
"stock_30_0",
"stock_30_5",
"stock_30_10",
"stock_30_20",
"stock_31_0",
"stock_31_5",
"stock_31_10",
"stock_31_20",
"stock_32_0",
"stock_32_5",
"stock_32_10",
"stock_32_20",
"stock_33_0",
"stock_33_5",
"stock_33_10",
"stock_33_20",
"stock_34_0",
"stock_34_5",
"stock_34_10",
"stock_34_20",
"stock_35_0",
"stock_35_5",
"stock_35_10",
"stock_35_20",
"stock_36_0",
"stock_36_5",
"stock_36_10",
"stock_36_20",
"stock_37_0",
"stock_37_5",
"stock_37_10",
"stock_37_20",
"stock_38_0",
"stock_38_5",
"stock_38_10",
"stock_38_20")


schs <- sepsss <- invalids <- c()
for (i in unique(data$identifier)) {
   schs <- c(schs, 0:(length(which(data$identifier==i))-1))

# index 18 is where harvests that do something start
# index 32 is where harvests that do something end
iiin <- names(which(rowSums(data[which(data$identifier==i),18:32])==0))
if (length(iiin)>1) {
  invalids <- c(invalids, iiin[-1])
}

# schedule 0 no mgmt
sepss <- c("donothing")
for (di in 2:length(which(data$identifier==i))) {
  cur.nmn <- names(data)[18:32][which(data[which(data$identifier==i)[di],18:32]>0)]
  seps <- NULL
  for (nmn in cur.nmn) {  
    if (is.null(seps)) {
      seps <- strsplit(nmn,"harvest_value_")[[1]][2]
    } else {
      seps <- paste(seps, strsplit(nmn,"harvest_value_")[[1]][2], sep=" + ")
    }
  }
  sepss <- rbind(sepss,seps)
}
if (is.null(sepsss)) {
  sepsss <- sepss
} else {
  sepsss <- rbind(sepsss,sepss)
}
}

if(length(invalids)>0) {
bdata <- cbind(count, data[-as.numeric(invalids),c(1:2)], sepsss)
} else {
bdata <- cbind(count, data[,c(1:2)], sepsss)
}

if(length(invalids)>0) {
bdata[,3] <- schs[-as.numeric(invalids)]
} else {
bdata[,3] <- schs
}
names(bdata) <- c("holding", "unit", "schedule", "treatment")

# indexes from 1-11 (also 33-44)
# indexes from 15-17 are total harvest values from the cutting years (2, 7, 17)
sdata <- cbind(count, data[,c(1:11, 15:17, 33:184)])
sdata[,3] <- schs
names(sdata)[c(1:3)] <- c("holding", "unit", "schedule")

if (count==1) { 
  sdatas <- sdata
  bdatas <- bdata
} else {
  sdatas <- rbind(sdatas, sdata)
  bdatas <- rbind(bdatas, bdata)
}

count <- count + 1

write.table(bdatas, "alternatives_key.csv", quote=F, row.names=F, sep=",")
write.table(sdatas, "alternatives.csv", quote=F, row.names=F, sep=",")
 
