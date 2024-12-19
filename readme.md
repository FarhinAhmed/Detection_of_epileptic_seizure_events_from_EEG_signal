# Detection of epileptic seizure events from EEG signal

## Tools used:
## MATLAB - Signal Processing Toolbox, Machine Learning Toolbox

Epilepsy is a neurological disorder that generates abnormal electrical activity in the brain and causes seizures. Electroencephalogram (EEG) is a widely used tool for detecting epileptic seizures. However, identifying subtle differences between epileptic and non-epileptic EEG signals through visual inspection is highly challenging due to individual variability, recording location, and electrode configurations. Manual inspection by neurologists is time-intensive and prone to errors. Automated detection methods can significantly reduce this burden. 
In this project, conducted at the University of Rochester, I have proposed an automated seizure detection approach using the delta and theta band power of EEG signals, employing support vector machine (SVM) and k-nearest neighbor (k-NN) classifiers. Ten-fold cross-validation on individual subjects yielded an average accuracy of 77.86%. Additionally, a composite training set from five subjects was tested on data from six independent subjects, achieving an average accuracy of 76.92%.
