clc
clear all
close all
filelist=dir('*.mat');
seizure_time=zeros(length(filelist),2);
%input the seizure times for this subject (available in the dataset)
seizure_time(6,:)=[298 320];
seizure_time(7,:)=[2695 2727];
seizure_time(8,:)=[1454 2206];
All_mat=[];
for i=1:length(filelist)
    i
    [new_feature,r,c]=feature_extraction(filelist(i).name);
%     new_feature= total_data(filelist(i),seizure_time(i,:));
    
    label=zeros(r,1);
    if (seizure_time(i,1))& (seizure_time(i,2))~=0
       label(seizure_time(i,1):seizure_time(i,2))=1;
    end
       name=filelist(i).name;
       new_mat=[string(new_feature) string(label) string((repmat(name,r,1))) string(ones(r,1)*seizure_time(i,:))];
       [All_mat]=[All_mat;new_mat];
end

   