# Challenge-2-London-Team-E

## Team Name: NaNs

### Team Members:
1. **Maitreyee Wairagkar** (Multiclass PCA_SVM based prediction model for CDR - Final Solution - accuracy 95.72%)
2. **Nabila Rahman** (Association of different factors with standardised cognitive score + correlations of features)
3. **Ana Lawry Aguila** (Associations between health, pollution and demographic measures with cognitive measures)
4. **Aadarsh Gupta** (Exploring multiclass ensemble methods for classification)
5. **Jordan Moore** (Exploring deep neural network for classification)
6. **Pooja Sarin** (Social media analysis and scientific literature review of features critical for dementia)
7. **Winnie (Cheng Wai) Lei** (Preprocessing data - imputing missing value features + FAMD preliminary data exploration)

## Project Title: A predictive model for dementia that overcomes the need for cohort harmonisation

Details of our solution and description of individual contributions are given below in contents.

### Contents:

## 1. Complete Pipeline for Multiclass PCA+SVD based Model for Clinical Dementia Rating Prediction
 
**Author: Maitreyee Wairagkar**

File **'Multiclass-Prediction-Model-LASI-DAD'** contains the **final solution** with the entire pipeline for Machine Learning prediction model for classifying 5 CDR scores. This includes data preprocessing and cleaning, data visualisation and transformation using T-SNE, PCA, Histogram, training PCA+SVM based multiclass model, and detailed evaluation of this model. Average accuracy of CDR prediction on unseen testing data using this model is **95.72%**

Initially data was preprocessed to remove non-relevant features and features with missing values. The cleaned dataset was visualised using PCA as shown in Fig1. PCA components clearly show separate clusters for 5 CDR scores as shown in Fig 1a (colour coded). Fig 1b shows that the variance of these features has a fat tail distribution which means all principal components contribute to the variance and hence we used all principal components for classification of CDR. Since PCA showed clear separable clusters with different CDR scores, SVM model was deemed to be suitable for classification of principal component based featyres. Fig 1c shows the histogram of original data with imbalanced classes which was upsampled to get equal number of samples in each class before training the model. 

![PCA of LASI-DAD features showing clear clusters according to CDR](https://github.com/DEMON-NEUROHACK/Challenge-2-London-Team-E/blob/main/PCA.png)

Fig 2a shows the pipeline of our final solution - a PCA based multiclass SVM and detailed analysis of its performance. We tuned hyperparameters for SVM using Grid Search method with cross validation. Then we used 10-fold cross validation scheme to train and evaluate our multiclass SVM model. We used strong regularisation in SVM to avoid overfitting and weighted penalty in case of imbalanced classes. PCA was trained on training set within wach fold of cross-validation. The normalised principal components were used as the training data. Same PCA weights and normalisation weights that were trained on training set were used to obtain principal components for the testing set. Testing data was then used to predict CDR scores and evaluate the model performance. Fig 2b shows the confusion matrix for all five classes. We can see that all classes were discriminated with high accuracy. CDR of 0 and 0.5 (no dementia and questionable impairment) were confused with each other the most, however, confirmed dementia diagnosis were classified almost perfectly. Fig 2c shows precision, recall and f1-score for all classes. Precision, recall and f1-score are all similar indicating that the model shows equally high performance on both negative and positive samples of each class. Overall accuracy of the model on training data was 98.81% and final **accuracy for test data was 95.72%**. This suggests that PCA+SVM based model has a good potential for use in dementia diagnosis.       
![FinalSolution - Multiclass Model Pipeline and Results ](https://github.com/DEMON-NEUROHACK/Challenge-2-London-Team-E/blob/main/FinalSolution-Multiclass_ModelSVM.png)

## 2. NN Model for using (sub)feature set to predict dementia onset

**Author: Jordan Moore**

Implementation of 2 layer model done in Pytorch. Data input and vectorisation needs completing. NN tried due to wide dataset, but issues on vectorisation due to unusual encoding of the features that needs to be examined


## 3. Multi-modal Ensemble method for Dementia Level classification

**Author: Aadarsh Gupta**

Implementation of Ensemble method (based on Voting Classifier) over 7 different Machine learning algorithms, after extensive pre-processing of data and choice of important features based on demographic, health and cognitive measures. The classifiers explored: Logistic Regression, LDA, K-Neighbors Classifier, Decision Tree Classifier, Random Forest Classifier, Gaussian Naïve Bayes and SVM; Emphasis on impact of label encoded features over raw values for each algorithm. Ensemble method accuracy : 78.9%.
The confusion matrix for the ensemble method obtained is :

![Confusion-matrix](https://github.com/DEMON-NEUROHACK/Challenge-2-London-Team-E/blob/main/Confusion-matrix.png)

The important features and their relative composition in determining clinical rating for:
 - Decision Tree Classifier: 

![Feature importance based on Decision-Tree classifier](https://github.com/DEMON-NEUROHACK/Challenge-2-London-Team-E/blob/main/Decision-tree-features.png)

 - Random Forest Classifier: 

![Feature importance based on Random-Forest classifier](https://github.com/DEMON-NEUROHACK/Challenge-2-London-Team-E/blob/main/Random-forest-features.png)

## 4. Preprocessing - Missing data imputation using the regularised iterative FAMD algorithm

**Author: Winnie (Cheng Wai) Lei**

File "Preprocessing_FAMD_Impute" is an R script that requires the input mixed data (the LASI-DAD continuous and categorical data) to be separated into to continuous ("processed_continuous_data.csv") and categorical ("processed_categorical_data.csv") csv files. The missing values are imputed using the regularised iterative FAMD algorithm with 5 components. 

## 5. Multiview analysis of Cognitive features with demographics, health and pollution features

**Author: Ana Lawry Aguila**

The notebook 'PLS_analysis.ipynb' contains a partial least squares (PLS) analysis considering a subset of the LASI-DAD data as two data views: a cognitive dataset, and a health, demographics and pollution dataset. We explore the associative effect between these two data views and provide an analysis of feature importance. We show that by viewing the dataset in this way, we find PLS projections that generalise to test data and reflect disease severity. This initial exploration of PLS could be extended in future work. For example, using the PLS components rather than PCA components as input for the classification model discussed in 1). 

Figure 1 shows the first PLS components projected onto the test data. We can see some stratification by disease severity suggesting that the strongest association between the two datasets is in part driven by disease effects.

![PLS projections](https://github.com/DEMON-NEUROHACK/Challenge-2-London-Team-E/blob/main/PLS_projections.png)

Figures 2 and 3 show the first set of PLS weights for the two datasets. These weight vectors can be directly interpreted as feature importance and give us an idea of which features contribute to an associative effect.
![xweights](https://github.com/DEMON-NEUROHACK/Challenge-2-London-Team-E/blob/main/xweights.png)
![yweights](https://github.com/DEMON-NEUROHACK/Challenge-2-London-Team-E/blob/main/yweights.png)

## 6. A Systematic Literature Review on "Dementia" and "LASI-DAD" through Scopus Database along with Twitter Analytics for feature enahancement of ML-based model

**Author: Pooja Sarin**

184 final stage journal articles identified out of more than 2 lakh total documents in the repository with major focus on Machine Learning, Deep Learning, Forecasting, Convolution Neural Network based techniques for Dementia Prediction. Also, there are total of 10 published papers on LASI-DAD on Scopus till date. Further, aprroximately 65K recent tweets through Twitter API are extracted and analysed to understand the discussions surrounding "Dementia" on Social Media. 


## 7. Preliminary data exploration using factor analysis of mixed data

**Author: Winnie Lei**

'FAMD_Data_Exploration.py' is a FAMD based analysis using 3 factors to identify key features in the raw data and the processed data (both imputed with the regularised iterative FAMD algorithmn). 'FAMD_Processed_column_correlation.csv' shows similar key features have been identified in the PCA and PLS analysis. 'FAMD_Processed_coordinates.csv' shows the calculated principal components that maybe used instead of the PCA input in the final multi-class model.

## 8. Feature Selection and Exploratory Analyis of LASI-DAD
 
**Author: Nabila**

File **'Nabila-Exploratory_Analysis_LASI-DAD'** contains the **exploratory analysis and visualition** performed using the LASI-DAD data

First I have selected and annotated 66 features that each summarises the individual tests performed in LASI-DAD. This includes cognitive scores, mini nutritional assessment, etc.

In the Exploratory analysis, I used a typical method used for gene analysis of variance to figure out which features are important. Here least squared regression is used to model each cognitive feature against each non-cognitive feature. Then the median R2 value is used to weight the importance of each non-cognitive feature.

For Correlation plot, this was a simple Pearson's correlation test between each cognitive feature and non cognitive feature.
