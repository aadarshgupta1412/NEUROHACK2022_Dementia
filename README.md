# Challenge-2-London-Team-E

## Team Name: NaNs

### Team Members:
1. **Maitreyee Wairagkar**
2. **Nabila Rahman**
3. **Ana Lawry Aguila** 
4. **Aadarsh Gupta**
5. **Jordan Moore**
6. **Pooja Sarin** 
7. **Winnie (Cheng Wai) Lei** 

### Contents:

## 1. Complete Pipeline for Multiclass PCA+SVD based Model for Clinical Dementia Rating Prediction

File 'Multiclass-Prediction-Model-LASI-DAD' contains the entire pipeline for Machine Learning prediction model for classifying 5 CDR scores. This includes data preprocessing and cleaning, data visualisation and transformation using T-SNE, PCA, Histogram, training PCA+SVM based multiclass model, and detailed evaluation of this model. Average accuracy of CDR prediction on unseen testing data using this model is 95.72%   
**Author: Maitreyee Wairagkar**

## 2. NN Model for using (sub)feature set to predict dementia onset

Implementation of 2 layer model done in Pytorch. Data input and vectorisation needs completing. NN tried due to wide dataset, but issues on vectorisation due to unusual encoding of the features that needs to be examined
**Author: Jordan Moore**

## 3. Multi-modal Ensemble Method for dementia rating prediction
 
Implementation of Ensemble Method (based on Voting Classifier) over 7 different Machine learning algorithms after extensive pre-processing of data and choice of important features to be considered based on demographic, health and cognitive measures. The algorithms shown are: Logistic Regression, LDA, K-Neighbors Classifier, Decision Tree Classifier, Random Forest Classifier, Gaussian Naive Bayes and SVM. Implementation also emphasizes on impact caused by label encoded features over raw data values for each classifier. 
**Author: Aadarsh Gupta**

## Others: Preprocessing -- Missing data imputation using the regularised iterative FAMD algorithm

File "Preprocessing_FAMD_Impute" is an R script that requires the input mixed data (the LASI-DAD continuous and categorical data) to be separated into to continuous ("processed_continuous_data.csv") and categorical ("processed_categorical_data.csv") csv files. The missing values are imputed using the regularised iterative FAMD algorithm with 5 components. 
**Author: Winnie (Cheng Wai) Lei**
