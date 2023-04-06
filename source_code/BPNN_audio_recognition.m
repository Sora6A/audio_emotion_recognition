function [Accuracy,Precison,Recall,ps,train_data_mean,...
tran_matrix,weight1,weight2,bias1,bias2]=BPNN_audio_recognition(positive_folder,...
    negetive_folder,test_ratio,PCA_dimension,hiden_layer_neuro_num,epoch_num,lr)                             %使用svm算法进行语音分类，正例文件夹，反例文件夹，随机抽取作为测试数据的比例。
%     clear;
%     positive_folder='..\audiofile\positive';
%     negetive_folder='..\audiofile\negative';
%     test_ratio=0.1;
%     PCA_dimension=100;
%     hiden_layer_neuro_num=100;
%     epoch_num=1000;
%     lr=0.05;
    input_format='wav';
    input_format=['*.',input_format];
    positive_files = dir(fullfile(positive_folder,input_format));
    negetive_files = dir(fullfile(negetive_folder,input_format));
    counter_dir='.\counter';
    if ~exist(counter_dir,'dir')
	    mkdir(counter_dir);
    end
    system('start python ./feature_extracting_waitbar.py');
    positive_L=length(positive_files);
    negetive_L=length(negetive_files);
    positive_test_num=fix(positive_L*test_ratio);
    negetive_test_num=fix(negetive_L*test_ratio);
    positive_trainning_num=positive_L-positive_test_num;
    negetive_trainning_num=negetive_L-negetive_test_num;
%     full_sample_num=positive_L+negetive_L;
    feature_num_per_sample=288;
    positive_features=zeros(positive_L,feature_num_per_sample);
    negetive_features=zeros(negetive_L,feature_num_per_sample);
%     core=14;
%     p=parpool('Threads',core);
    disp('Feature Extracting...');
    parfor i=1:positive_L
        positive_features(i,:)=featurevector...
            (fullfile(positive_files(i).folder,...
            positive_files(i).name),640,160);
        counter_filename=strcat('positive_',num2str(i),'.txt');
        counter_filename=fullfile(counter_dir,counter_filename);
        fp=fopen(counter_filename,'a');
        fprintf(fp,'%d ',1);
        fclose(fp);
    end
    parfor i=1:negetive_L
        negetive_features(i,:)=featurevector...
            (fullfile(negetive_files(i).folder,...
            negetive_files(i).name),640,160);
        counter_filename=strcat('negetive_',num2str(i),'.txt');
        counter_filename=fullfile(counter_dir,counter_filename);
        fp=fopen(counter_filename,'a');
        fprintf(fp,'%d ',1);
        fclose(fp);
    end
    disp('Feature Extracting Complete!');
    disp('Normalization...');
    positive_test_features_extract_num=randperm(positive_L,positive_test_num);
    negetive_test_features_extract_num=randperm(negetive_L,negetive_test_num);
    positive_test_features=positive_features(positive_test_features_extract_num,:);
    negetive_test_features=negetive_features(negetive_test_features_extract_num,:);
    positive_trainning_features=positive_features;
    positive_trainning_features(positive_test_features_extract_num,:)=[];
    negetive_trainning_features=negetive_features;
    negetive_trainning_features(negetive_test_features_extract_num,:)=[];
    full_tranning_vector=zeros(positive_trainning_num+...
        negetive_trainning_num,feature_num_per_sample);
    full_test_vector=zeros(positive_test_num+negetive_test_num,feature_num_per_sample);
    full_tranning_vector(1:positive_trainning_num,:)=positive_trainning_features;
    full_tranning_vector((positive_trainning_num+1):(positive_trainning_num+...
        negetive_trainning_num),:)=negetive_trainning_features;
    full_test_vector(1:positive_test_num,:)=positive_test_features;
    full_test_vector((positive_test_num+1):(positive_test_num+negetive_test_num),:)=negetive_test_features;
    [full_tranning_vector_normalized,ps]=mapstd(full_tranning_vector.',0,1);
    full_tranning_vector_normalized=full_tranning_vector_normalized.';
    full_test_vector_normalized=mapstd('apply',full_test_vector.',ps);
    full_test_vector_normalized=full_test_vector_normalized.';
%     [full_tranning_vector_normalized,ps]=mapminmax(full_tranning_vector.');
%     full_tranning_vector_normalized=full_tranning_vector_normalized.';
%     full_test_vector_normalized=mapminmax('apply',full_test_vector.',ps);
%     full_test_vector_normalized=full_test_vector_normalized.';
%     [full_tranning_vector_normalized,ps]=mapminmax(full_tranning_vector_normalized.');
%     full_tranning_vector_normalized=full_tranning_vector_normalized.';
%     full_test_vector_normalized=mapminmax('apply',full_test_vector_normalized.',ps);
%     full_test_vector_normalized=full_test_vector_normalized.';
    test_set_labels=zeros(positive_test_num+negetive_test_num,1);
    train_set_labels=zeros(positive_trainning_num+negetive_trainning_num,1);
    train_set_labels(1:positive_trainning_num)=1;
    train_set_labels((positive_trainning_num+1):(positive_trainning_num+negetive_trainning_num))=-1;
    test_set_labels(1:positive_test_num)=1;
    test_set_labels((positive_test_num+1):(positive_test_num+negetive_test_num))=-1;
    disp('Normalization Complete!');
    disp('PCA...');
    [coeff,~] = pca(full_tranning_vector_normalized);
    tran_matrix = coeff(:,1:PCA_dimension);
    train_data0 = bsxfun(@minus,full_tranning_vector_normalized,mean(full_tranning_vector_normalized,1));
    low_train_data = train_data0 * tran_matrix;
    train_data_mean=mean(full_tranning_vector_normalized,1);
    test_data0 = bsxfun(@minus,full_test_vector_normalized,train_data_mean); 
    low_test_data = test_data0 * tran_matrix;
    disp('PCA Complete!');
    disp('BPNN Tranning...');
    weight1=normrnd(0,0.01,[PCA_dimension hiden_layer_neuro_num]);
    weight2=normrnd(0,0.01,[hiden_layer_neuro_num 2]);
    bias1=normrnd(0,0.01,[1 hiden_layer_neuro_num]);
    bias2=normrnd(0,0.01,[1 2]);
    train_label=zeros(length(train_set_labels),2);
    for i=1:length(train_set_labels)
        if train_set_labels(i)==1
            train_label(i,1)=1;
        else
            train_label(i,2)=1;
        end
    end
    for k=1:epoch_num
        [weight1,weight2,bias1,bias2,loss]=Train_BP_Neural_Network...
            (low_train_data,weight1,weight2,bias1,bias2,train_label,lr);
        disp(['epoch:',num2str(k),'/',num2str(epoch_num),'  loss:',num2str(loss)]);
    end
    disp('trainning complete!');
    test_label=BP_Neural_Network_Predict(weight1,weight2,bias1,...
        bias2,low_test_data);
    predict_label=zeros(length(test_set_labels),1);
    for i=1:length(test_set_labels)
        if test_label(i,1)>=test_label(i,2)
            predict_label(i)=1;
        else
            predict_label(i)=-1;
        end
    end
    TP=0;
    FP=0;
    FN=0;
    TN=0;
    for i=1:length(test_set_labels)
        if (predict_label(i)==1)&&(test_set_labels(i)==1)
            TP=TP+1;
        end
    end
    for i=1:length(test_set_labels)
        if (predict_label(i)==1)&&(test_set_labels(i)==-1)
            FP=FP+1;
        end
    end 
    for i=1:length(test_set_labels)
        if (predict_label(i)==-1)&&(test_set_labels(i)==1)
            FN=FN+1;
        end
    end
    for i=1:length(test_set_labels)
        if (predict_label(i)==-1)&&(test_set_labels(i)==-1)
            TN=TN+1;
        end
    end
    Accuracy=(TP+TN)/length(test_set_labels);
    Precison=TP/(TP+FP);
    Recall=TP/(TP+FN);
    disp('BPNN Tranning and predict Complete!');
    disp(['Accuracy=',num2str(Accuracy)]);
    disp(['Precison=',num2str(Precison)]);
    disp(['Recall=',num2str(Recall)]);
end

% function output=d_sigmoid(input)
%     output=sigmoid(input).*(1-sigmoid(input));
% end
% 
% function output=sigmoid(input)
%     output=1./(1.+exp(-1.*input));
% end