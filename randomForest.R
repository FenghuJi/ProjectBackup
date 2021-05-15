library(randomForest)
library(dplyr)

asv_abundance <- read.table("asv_realtive_abundance.txt",header=T,sep="\t")

## 过滤ASV
count <- function(abun){
	num <- sum(abun>0.1)
	return(num)
}

filter_count <- apply(asv_abundance,1,count)
index <- which(filter_count>2)
filter_abundance <- asv_abundance[index,]

## 读入metadata数据并合并数据
metadata <- read.table("healthy_metadata.txt",header=T,sep="\t",row.names=1)
otu <- data.frame(t(filter_abundance))
otu$age_day <- metadata[rownames(otu),]$age_day


# 训练集和测试集

train_sample <- otu[rownames(metadata)[which(metadata$Type=="Training")],]
valid_sample <- otu[rownames(metadata)[which(metadata$Type=="Validation")],]

set.seed(123)
otu_train.forest <- randomForest(age_day~., data = train_sample, importance = TRUE)

otu_train.forest

# Call:
#  randomForest(formula = age_day ~ ., data = train_sample, importance = TRUE) 
#                Type of random forest: regression
#                      Number of trees: 500
# No. of variables tried at each split: 201

#           Mean of squared residuals: 16052.08
#                     % Var explained: 64.26

#查看该模型的预测性能，可以看到具有较高的精度。
importance_otu <- data.frame(otu_train.forest$importance)

varImpPlot(otu_train.forest, n.var = min(30, nrow(otu_train.forest$importance)),
    main = 'Top 30 - variable importance')

#“%IncMSE”即increase in mean squared error，通过对每一个预测变量随机赋值，如果该预测变量更为重要，那么其值被随机替换后模型预测的误差会增大。因此，该值越大表示该变量的重要性越大；

#“IncNodePurity”即increase in node purity，通过残差平方和来度量，代表了每个变量对分类树每个节点上观测值的异质性的影响，从而比较变量的重要性。该值越大表示该变量的重要性越大。
importance_otu <- importance_otu[order(importance_otu$IncNodePurity, decreasing = TRUE), ]
write.table(importance_otu, 'importance_otu.txt', sep = '\t', col.names = NA, quote = FALSE)

# 十折交叉验证
otu_train <- train_sample

set.seed(123)
otu_train.cv <- replicate(10, rfcv(otu_train[-ncol(otu_train)], otu_train$age_day, cv.fold = 10, step = 1.5), simplify = FALSE)
otu_train.cv
otu_train.cv <- data.frame(sapply(otu_train.cv, '[[', 'error.cv'))
otu_train.cv$otus <- rownames(otu_train.cv)
otu_train.cv <- reshape2::melt(otu_train.cv, id = 'otus')
otu_train.cv$otus <- as.numeric(as.character(otu_train.cv$otus))
 
otu_train.cv.mean <- aggregate(otu_train.cv$value, by = list(otu_train.cv$otus), FUN = mean)
head(otu_train.cv.mean, 10)
 
#拟合线图
library(ggplot2)
 
ggplot(otu_train.cv.mean, aes(Group.1, x)) +
geom_line() +
theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent')) +
labs(title = '',x = 'Number of OTUs', y = 'Cross-validation error')


## 还是24个变量
#    Group.1        x
# 1        1 31429.30
# 2        2 26456.25
# 3        3 25148.69
# 4        5 23241.56
# 5        7 21539.77
# 6       10 19436.76
# 7       16 17416.52
# 8       24 16545.88
# 9       35 16572.65
# 10      53 16623.18


importance_otu.select <- importance_otu[1:24, ]
otu_id.select <- rownames(importance_otu.select)
write.table(importance_otu.select, 'importance_otu.select.txt', sep = '\t', col.names = NA, quote = FALSE)

##只包含 24 个重要预测变量的简约回归
otu_train.select <- train_sample[,c(otu_id.select,"age_day")]
otu_test.select <- valid_sample[,c(otu_id.select,"age_day")]

#随机森林计算（默认生成 500 棵决策树），详情 ?randomForest
set.seed(123)
otu_train.select.forest <- randomForest(age_day~., data = otu_train.select, importance = TRUE)
otu_train.select.forest

# Call:
#  randomForest(formula = age_day ~ ., data = otu_train.select,      importance = TRUE) 
#                Type of random forest: regression
#                      Number of trees: 500
# No. of variables tried at each split: 8

#           Mean of squared residuals: 15775.65
#                     % Var explained: 64.88
#使用训练集，查看预测精度
# age_predict <- predict(otu_train.select.forest, otu_train.select)

# plot(otu_train.select$age_day, age_predict, main = '训练集',
#     xlab = 'age (days)', ylab = 'Predict')
# abline(1, 1)

#使用测试集，评估预测性能
# 将测试集分为signle和多胞胎
names <- rownames(otu_test.select)
single_test <- otu_test.select[which(grepl(names,pattern="Bgs")==TRUE),]
twin_test <- otu_test.select[which(grepl(names,pattern="Bgt")==TRUE),]

age_predict_single <- predict(otu_train.select.forest, single_test)
 
plot(single_test$age_day, age_predict_single, main = 'Single_Birth_Cohort_Test',
    xlab = 'age (days)', ylab = 'Predict')
abline(1, 1)
Single_data <- data.frame(single_test$age_day, age_predict_single)
Single_Fit <- lm(age_predict_single~single_test.age_day,data=Single_data)
summary(Single_Fit)
# lm(formula = age_predict_single ~ single_test.age_day, data = Single_data)

# Residuals:
#     Min      1Q  Median      3Q     Max 
# -340.61  -66.55   -1.79   64.17  308.49 

# Coefficients:
#                      Estimate Std. Error t value Pr(>|t|)    
# (Intercept)         115.47437   11.21453   10.30   <2e-16 ***
# single_test.age_day   0.71214    0.02805   25.39   <2e-16 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Residual standard error: 97.4 on 274 degrees of freedom
# Multiple R-squared:  0.7017,	Adjusted R-squared:  0.7006 
# F-statistic: 644.5 on 1 and 274 DF,  p-value: < 2.2e-16

age_predict_twin <- predict(otu_train.select.forest, twin_test)
 
plot(twin_test$age_day, age_predict_twin, main = 'Twin_Birth_Cohort_Test',
    xlab = 'age (days)', ylab = 'Predict')
abline(1, 1)
Twin_data <- data.frame(twin_test$age_day, age_predict_twin)
Twin_Fit <- lm(age_predict_twin~twin_test.age_day,data=Twin_data)
summary(Twin_Fit)
# Call:
# lm(formula = age_predict_twin ~ twin_test.age_day, data = Twin_data)

# Residuals:
#     Min      1Q  Median      3Q     Max 
# -327.22  -69.77   -8.07   63.99  276.87 

# Coefficients:
#                    Estimate Std. Error t value Pr(>|t|)    
# (Intercept)       138.39419    8.50596   16.27   <2e-16 ***
# twin_test.age_day   0.75011    0.02635   28.47   <2e-16 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Residual standard error: 100.3 on 445 degrees of freedom
# Multiple R-squared:  0.6455,	Adjusted R-squared:  0.6447 
# F-statistic: 810.4 on 1 and 445 DF,  p-value: < 2.2e-16


## 获取ASV信息
taxa <- read.table("taxonomy.txt",sep="\t",header=T)
ASV <- read.table("Selected_OTUs.txt",header=T,sep="\t")
ASV_infomation <- inner_join(taxa,ASV,by="ASV")
show_info <- ASV_infomation[,c(1,2)]
show_info