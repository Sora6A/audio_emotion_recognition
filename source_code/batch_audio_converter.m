function result=batch_audio_converter(folder,input_format,...               %批量音频转换分割，（文件夹路径，输入格式，输出格式，分割时长， 输出采样率）
    output_format,Time,out_fs)
% folder='C:\Users\90317\Desktop\audio_work_dir\audiofile';
% input_format='mp3';
% output_format='wav';
% Time=60;
% out_fs=16000;
    input_format=['*.',input_format];
    output_format=['.',output_format];
    output_dir=[folder,'\output'];
    if ~exist(output_dir,'dir')
	    mkdir(folder, 'output');
    end
    files = dir(fullfile(folder,input_format)); % list of mp3 files in folder
    L=length(files);
    bar=my_waitbar_init(1);
    ll=zeros(L,1);
    Fs=zeros(L,1);
    tunnel=zeros(L,1);
    for i=1:L
        [y,Fs(i)]=audioread(fullfile(files(i).folder,files(i).name));
        tunnel(i)=size(y,2);
        ll(i)=fix(length(y)/(Time*Fs(i)))-1;
    end
    LoopNum=sum(ll);
    count=0;
    my_waitbar(LoopNum,0,bar);
    for i=1:L  
        samples=zeros(ll(i),2);
        y=zeros(Time*Fs(i),tunnel(i),ll(i));
        y_to_convert=audioread(fullfile(files(i).folder,files(i).name));
        for j=1:ll(i)
            samples(j,:)=[((j-1)*Time*Fs(i))+1,(j*Time*Fs(i))];
            for k=1:tunnel(i)
                y(:,k,j)=y_to_convert(samples(j,1):samples(j,2),k);
            end
%             y(:,:,j)=audioread(files{i},samples(j,:));
            str1=files(i).name; 
            filename=strcat(str1(1:end-4),'_',num2str(j),output_format);
            filename=fullfile(output_dir,filename);
            y_toWrite=y(:,:,j);
            if out_fs==0
                audiowrite(filename,y_toWrite,Fs(i)); %转写成.wav格式文件
            else
                if Fs(i)==out_fs
                    audiowrite(filename,y_toWrite,out_fs); %转写成.wav格式文件
                else
                    [P,Q] = rat(out_fs/Fs(i));
                    new_y_toWrite=resample(y_toWrite,P,Q);
                    audiowrite(filename,new_y_toWrite,out_fs);
                end
            end
            count=count+1;
            my_waitbar(LoopNum,count,bar);
        end
    end
    close(bar);
    result=1;
end

function bar=my_waitbar_init(i)
    bar=waitbar(i,'Loading your data');
end

function []=my_waitbar(loopNum,i,bar)
tic;
    pause(1)  % replace your real code
    currentProgress = roundn((i/loopNum)*100,-1);
    remainingTime = roundn((loopNum-i)*toc/60,-1);
    barString = ['Current Progress:',num2str(currentProgress),...
        '%, Remaining Time:',num2str(remainingTime),'min :)'];
    waitbar(i/loopNum,bar,barString);
toc;
end