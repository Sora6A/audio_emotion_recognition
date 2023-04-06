function result=audio_segment_extarct(folder,input_file_name,...               %提取某音频文件中从指定处开始的指定长度一段，（文件夹路径，输入文件名，输出格式，提取开始时刻，提取时长）
    output_format,start,duration)
% folder='C:\Users\90317\Desktop\audio_work_dir\audiofile';
% input_file_name='蒋奉儒1.mp3';
% output_format='wav';
% start=60;
% duration=10;
output_format=['.',output_format];
    output_dir=[folder,'\output'];
    if ~exist(output_dir,'dir')
	    mkdir(folder, 'output');
    end
    file = dir(fullfile(folder,input_file_name)); % list of mp3 files in folder
    [y,Fs]=audioread(fullfile(file.folder,file.name));
    y=y(start*Fs:(start+duration)*Fs,:);
    str1=file.name; 
    filename=strcat(str1(1:end-4),'_start_',num2str(start),'_duration',num2str(duration),output_format);
    audiowrite(fullfile(output_dir,filename),y,Fs); %转写成.wav格式文件
    result=1;
end