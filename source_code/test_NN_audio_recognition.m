positive_folder='..\audiofile\positive';
negetive_folder='..\audiofile\negative';
test_ratio=0.25;
PCA_dimension=66;
hiden_layer_neuro_num=100;
epoch_num=100;
lr=0.05;
[Accuracy,~]=BPNN_audio_recognition(positive_folder,...
    negetive_folder,test_ratio,PCA_dimension,...
    hiden_layer_neuro_num,epoch_num,lr);   