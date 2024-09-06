# NMF-for-Identification-of-Human-Disease-Associated-Gene-Module
 NMF for Identification of Human Disease-Associated Gene Modules through Multi-label Classification
 #find unique gene 
#load all_gene dataset And intersect(using MATLAB)
 [a,b] = find(ismember(allgenedata,unique_ppi));
 result=x([a],:);
