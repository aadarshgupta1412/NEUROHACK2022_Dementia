# Challenge-2-London-Team-E

## Team Name: NaNs

### Team Members:
1. **Maitreyee Wairagkar** (Multiclass PCA_SVM based prediction model for CDR - Final Solution - accuracy 95.72%)
2. **Nabila Rahman** (Association of different factors with standardised cognitive score + correlations of features)
3. **Ana Lawry Aguila** (Associations between health, pollution and demographic measures with cognitive measures)
4. **Aadarsh Gupta** (Exploring multiclass ensamble methods for classification)
5. **Jordan Moore** (Exploring deep neural network for classification)
6. **Pooja Sarin** (Social media analysis and scientific literature review of features critical for dementia)
7. **Winnie (Cheng Wai) Lei** (Preprocessing data - imputing missing value features)

Details of our solution and description of individual contributions are given below in contents.

### Contents:

## 1. Complete Pipeline for Multiclass PCA+SVD based Model for Clinical Dementia Rating Prediction

File 'Multiclass-Prediction-Model-LASI-DAD' contains the entire pipeline for Machine Learning prediction model for classifying 5 CDR scores. This includes data preprocessing and cleaning, data visualisation and transformation using T-SNE, PCA, Histogram, training PCA+SVM based multiclass model, and detailed evaluation of this model. Average accuracy of CDR prediction on unseen testing data using this model is 95.72%   
**Author: Maitreyee Wairagkar**

## 2. NN Model for using (sub)feature set to predict dementia onset

Implementation of 2 layer model done in Pytorch. Data input and vectorisation needs completing. NN tried due to wide dataset, but issues on vectorisation due to unusual encoding of the features that needs to be examined
**Author: Jordan Moore**

## Others: Preprocessing -- Missing data imputation using the regularised iterative FAMD algorithm

File "Preprocessing_FAMD_Impute" is an R script that requires the input mixed data (the LASI-DAD continuous and categorical data) to be separated into to continuous ("processed_continuous_data.csv") and categorical ("processed_categorical_data.csv") csv files. The missing values are imputed using the regularised iterative FAMD algorithm with 5 components. 
**Author: Winnie (Cheng Wai) Lei**
