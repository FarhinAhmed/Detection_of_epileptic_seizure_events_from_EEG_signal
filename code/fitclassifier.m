clc
clear all
close all
load 'band_features01_03_ch1.mat'
Xtrain=[band_features_n(1:100,1:2);band_features_s(1:30,1:2)];
labels=[ones(100,1); 2*ones(30,1)];
Xtest=[band_features_s(31:41,1:2);band_features_n(3200:3210,1:2)];
mdl = fitcknn(Xtrain,labels); 
[label,score,cost] = predict(mdl,Xtest);