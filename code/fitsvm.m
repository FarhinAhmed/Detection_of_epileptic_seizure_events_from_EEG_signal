clc
clear all
close all
load All_mat1.mat
mat1=All_mat;
load All_mat2
mat2=All_mat;
total_mat=[mat1;mat2];
labels=str2double(total_mat(:,3));

feature_mat=str2double([total_mat(:,1) total_mat(:,2)]); 
sandn=[feature_mat labels];
ind_seizure=find((sandn(:,3)==1));
ind_normal=find((sandn(:,3)==0));
seizure_mat=sandn(ind_seizure,1:2);
normal_mat=sandn(ind_normal,1:2);
[rs,cs]=size(seizure_mat);
trs=round(0.2*rs);

[rn,cn]=size(normal_mat);
trn=round(0.2*rn);
indtest_s=randsample(rs,trs);
indtest_n=randsample(rn,trn);
test_seizure_w_label=[seizure_mat(indtest_s,1:2) ones(size(indtest_s,1),1)];
test_normal_w_label=[normal_mat(indtest_n,1:2) zeros(size(indtest_n,1),1)];
normal_mat_all=normal_mat(indtest_n,1:2);

test_set=[seizure_mat(indtest_s,1:2);normal_mat_all(1:130,:)];
true_label=[test_seizure_w_label(:,3);test_normal_w_label(1:130,3)];
% [trainingset]= intersect(a,randsample(a,8))      % gives training set with 8 random samples from a ;you can set what size your trainign set needs to be
% testset=a(~ismember(a,trainingset))    
% true_label=[
train_seizure=seizure_mat(~ismember(1:rs,indtest_s),1:2);
train_normal_all=normal_mat(~ismember(1:rn,indtest_n),1:2);
train_normal=train_normal_all(1:500,:);
train_label=[ones(size(train_seizure,1),1); zeros(size(train_normal,1),1)];

train_set=[train_seizure;train_normal];
% plot(train_seizure,'.','r')

mysvm=fitcsvm(train_set,train_label,'KernelFunction','RBF');
label = predict(mysvm,test_set);
accu_s=sum(true_label(1:125,1)==label(1:125,1))/length(test_seizure_w_label);
accu_n=sum(true_label(126:255,1)==label(126:255,1))/length(normal_mat_all(1:130,:));
%easy way
% classperf(true_label,label)
% mysvm=fitcsvm(train_set,train_label,'Standardize',true,'KernelFunction','RBF',...
%     'KernelScale','auto');
% CVSVMModel = crossval(mysvm);
% res_label=kfoldPredict(CVSVMModel);
 %label = predict(mysvm,test_set);
 figure(1)
plot(train_seizure(:,1),'linewidth',2); hold on
plot(train_normal(:,1),'linewidth',2); legend('seizure','normal')
title('delta energy of a training set')
%show how distinct the feature is between normal and seizure
figure(2)
plot(train_seizure(:,2),'linewidth',2); hold on
plot(train_normal(:,2),'linewidth',2); legend('seizure','normal')
title('theta energy of a training set')

% figure(3);
% hgscatter = gscatter(train_set(:,1),train_set(:,2),train_label);
% hold on;
% h_sv=plot(mysvm.SupportVectors(:,1),mysvm.SupportVectors(:,2),'ko','markersize',8);

% Confusionmat
c=confusionmat(true_label,label);
figure(3)
plotConfMat(c)
