# NMF-for-Identification-of-Human-Disease-Associated-Gene-Module
# NMF for Identification of Human Disease-Associated Gene Modules through Multi-label Classification
 #find unique gene 
# load all_gene dataset And intersect(using MATLAB)
 [a,b] = find(ismember(allgenedata,unique_ppi));
 result=x([a],:);
 # in R studio, create a semantic similarity matrix.
 install.packages(csgl.go);
 library("csbl.go"); 
set.prob.table(organism=9606 , type="similarity")
ent <- entities.from.text("Godata_nonprogressor.txt")
matrix<-entity.sim.many(ent,"BP","Relevance")
write.table(matrix,"simantic_sim_matrix_nonprogressor",sep="\t")
# After that find the top gene sets using the quantile function.
quantile(sim_mat,c(0.75,0.85,0.95),na.rm = TRUE)
#create an adjacency matrix using the PPI dataset and the semantic similarity matrix using GO_dataset.(Adj_mat.m)
# find optimal number of cluster using silhouette method(for both GO and PPI dataset)
library(factoextra)
fviz_nbclust(sim_matrix, kmeans, method = "silhouette")+
      labs(subtitle = "Silhouette method")
fviz_nbclust(adj_matrix, kmeans, method = "silhouette")+
      labs(subtitle = "Silhouette method")
# The two networks (PPI Network and semantic similarity network) are partitioned using a simple k-means clustering.
sim_matrix_cluster=kmeans(sim_matrix,centers = 6,nstart = 25)
adj_matrix_cluster=kmeans(adj_matrix,centers = 4,nstart = 25)
#The two categories of modules are then integrated by using a matrix factorization-based framework.
#create NMF object
library(NMF)
res <- nmf(total_mat, 7,"lee",nrun=40,maxIter=20,seed=1234)
# Seven meta-modules identified.
cluster=list()
for (i in 1:7){
      cluster[[i]]=gene_name[which(predict(res)==i)]
}

# Calculate Z − Score of each meta module mi and for each disease class dj as : Z Scoremi,dj = (Gmi,dj − Emi )/SQRT(Emi ) ,
# based on Z_score create a plot using ggplot2 package 
library(reshape)
library(ggplot2)
z_score.new=melt(final_module_Z_score);
#plot geom_linear
row_name=c("Module1","Module2","Module3","Module4","Module5","Module6","Module7");
final_module_Z_score=cbind(row_name,final_module_Z_score);
z_score.new=melt(final_module_Z_score,id=c("row_name"));
ggplot(z_score.new, aes(Disease_class, Z_score, ymin = 0,ymax = z_score.new)) + geom_linerange(data = z_score.new, aes(col = ifelse(Z_score <0, 'yellow', 'blue')),stat = "identity", position = "identity",size=0.5)+facet_grid(Modules~.)+theme_bw()+ theme(axis.text.x = element_text(angle = 90, hjust = 1))+theme(legend.position="none")+geom_hline(yintercept=1.5, linetype="dashed", color = "grey")+geom_hline(yintercept=-1.5, linetype="dashed", color = "grey");
ggsave(filename = "zscore.png",height=700,width=500,units = "in");
#Now a threshold is selected to label a meta module to one of the disease classes.
#For example,meta-module-2 is predicted to be included in seven of the disease classes (Neurological, Ophthalmological, Skeletal, #Cardiovascular, Cancer, and Endocrine)
# Heatmap of meta-modules;
#Figure 3, shows a heatmap representing the membership of the meta-modules within each of the disease classes.
library(RColorBrewer) 
library(ComplexHeatmap)
mat_hat=as.matrix.data.frame(b)
rownames(b)=gene_name
ht <- Heatmap(
      mat_hat,
     col = structure(brewer.pal(9, "RdYlBu")),
     cluster_rows = TRUE,
     cluster_columns =  TRUE,
     show_column_names = TRUE,
     row_names_gp = gpar(fontsize = 5.5),
     column_names_gp = gpar(fontsize = 1.5),
     row_dend_width = unit(4, "cm"),
     column_dend_height = unit(4, "cm"),
     rect_gp = gpar(col = "white", lwd = 0.5)
     )
# The labeled network was then trained using multi-label classification, here we used Utiml packages for ML classification.
install.packages("utiml")
library("utiml")
library(mldr)
##for SVM learner
library(e1071)
##for random forest learner
library(randomForest)
#for decission tree learner
library(rpart)
##for gradient boosting machine
library(gbm)
##for fuzzy rule based classifier
library(frbs)

library(kknn)
library("C50")

class_num=class_logical*1
my_data=cbind(sim_adj_mat,class_num)
#create mldr object from data frame

mymldr <- mldr_from_dataframe(my_data, labelIndices = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22), name = "testMLDR")

# Create two partitions (train and test) of my multi-label dataset
ds <- create_holdout_partition(mymldr, c(train=0.65, test=0.35))


# Using SVM method

#Create a Binary Relevance Model using the e1071::SVM method
brmodel <- br(ds$train, "SVM", seed=123)
prediction <- predict(brmodel, ds$test)
prediction

head(as.bipartition(prediction))
head(as.ranking(prediction))
newpred <- rcut_threshold(prediction, 2)
result <- multilabel_evaluate(ds$tes, prediction, "bipartition")
thresres <- multilabel_evaluate(ds$tes, newpred, "bipartition")
print(round(cbind(Default=result, RCUT=thresres), 3))
my_result=cbind(result,thresres)
library("csv")
write.csv(my_result,"my_result.csv")

# Using RandomForest method

#Create three partitions (train, val, test) of emotions dataset
partitions <- c(train = 0.6, val = 0.2, test = 0.2)
ds <- create_holdout_partition(mymldr, partitions, method="iterative")

#Create an Ensemble of Classifier Chains using Random Forest (randomForest package)
eccmodel <- ecc(ds$train, "RF", m=3, cores=parallel::detectCores(), seed=123)
val <- predict(eccmodel, ds$val, cores=parallel::detectCores())
test <- predict(eccmodel, ds$test, cores=parallel::detectCores())
thresholds <- scut_threshold(val, ds$val, cores=parallel::detectCores())
new.val <- fixed_threshold(val, thresholds)
new.test <- fixed_threshold(test, thresholds)
measures <- c("subset-accuracy", "F1", "hamming-loss", "macro-based") 

result <- cbind(
  Test = multilabel_evaluate(ds$tes, test, measures),
  TestWithThreshold = multilabel_evaluate(ds$tes, new.test, measures),
  Validation = multilabel_evaluate(ds$val, val, measures),
  ValidationWithThreshold = multilabel_evaluate(ds$val, new.val, measures)
)

print(round(result, 3))
write.csv(result,"my_result_randomForest.csv")

# KNN classifier
partitions <- c(train=0.7, test=0.2, val=0.1)
strat <- create_holdout_partition(mymldr, partitions, "iterative")
#Using KNN with k = 5 and changing ECC parameters
model1 <- ecc(strat$train, "KNN", m=7, subsample=0.8, k=5)
val <- predict(model1, strat$val)
test <- predict(model1, strat$test)
thresholds <- scut_threshold(val, strat$val)
new.val <- fixed_threshold(val, thresholds)
new.test <- fixed_threshold(test, thresholds)

measures <- c("subset-accuracy", "F1", "hamming-loss", "macro-based")
result <- cbind(
  Test = multilabel_evaluate(strat$tes, test, measures),
  TestWithThreshold = multilabel_evaluate(strat$tes, new.test, measures),
  Validation = multilabel_evaluate(strat$val, val, measures),
  ValidationWithThreshold = multilabel_evaluate(strat$val, new.val, measures)
)
write.csv(result,"my_result_knn.csv")


# Decision tree classifier
partitions <- c(train=0.7, test=0.2, val=0.1)
strat <- create_holdout_partition(mymldr, partitions, "iterative")
#using C5.0 and changing ECC parameters
model1=ecc(strat$train,"C5.0",subsample = 0.6,attr.space = 1)
val <- predict(model1, strat$val)
test <- predict(model1, strat$test)
thresholds <- scut_threshold(val, strat$val)
new.val <- fixed_threshold(val, thresholds)
new.test <- fixed_threshold(test, thresholds)

measures <- c("subset-accuracy", "F1", "hamming-loss", "macro-based")
result <- cbind(
  Test = multilabel_evaluate(strat$tes, test, measures),
  TestWithThreshold = multilabel_evaluate(strat$tes, new.test, measures),
  Validation = multilabel_evaluate(strat$val, val, measures),
  ValidationWithThreshold = multilabel_evaluate(strat$val, new.val, measures)
)
write.csv(result,"my_result_decision_tree.csv")


# Naive bayes classifier

partitions <- c(train=0.7, test=0.2, val=0.1)
strat <- create_holdout_partition(mymldr, partitions, "iterative")
#using C5.0 and changing ECC parameters
model1=ecc(strat$train,"NB",subsample = 0.6,attr.space = 1)
val <- predict(model1, strat$val)
test <- predict(model1, strat$test)
thresholds <- scut_threshold(val, strat$val)
new.val <- fixed_threshold(val, thresholds)
new.test <- fixed_threshold(test, thresholds)

measures <- c("subset-accuracy", "F1", "hamming-loss", "macro-based")
result <- cbind(
  Test = multilabel_evaluate(strat$tes, test, measures),
  TestWithThreshold = multilabel_evaluate(strat$tes, new.test, measures),
  Validation = multilabel_evaluate(strat$val, val, measures),
  ValidationWithThreshold = multilabel_evaluate(strat$val, new.val, measures)
)
write.csv(result,"my_result_Naive_bayes.csv");
# The table shows the multi-label Hamming loss of the classification results.
#The Hamming Loss ranges from 0 to 1, where lower values indicate better performance. A Hamming Loss of 0 means perfect predictions, while a value of 1 indicates that all labels are incorrectly predicted.
#The trained model is utilized to predict the class label of unknown genes. We have obtained 3131 genes associated with multiple disease classes.
# Table shows some of the predicted gene-disease associations along with their supporting PubMed IDs, demonstrating the validity of these associations based on existing literature.

 
