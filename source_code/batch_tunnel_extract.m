function result=batch_tunnel_extract(folder,input_format,...
    output_format,out_tunnel)
%     folder='C:\Users\90317\Desktop\audio_work_dir\audiofile\output';
%     input_format='wav';
%     output_format='wav';
%     out_tunnel=1;
    input_format=['*.',input_format];
    output_format=['.',output_format];
    output_dir=[folder,'\output'];
    if ~exist(output_dir,'dir')
	    mkdir(folder, 'output');
    end
    files = dir(fullfile(folder,input_format)); % list of mp3 files in folder
    L=length(files);
    bar=my_waitbar_init(1);
    my_waitbar(L,0,bar);
    count=0;
    for i=1:L  
            [y,Fs]=audioread(fullfile(files(i).folder,files(i).name));
            str1=files(i).name; 
            filename=strcat(str1(1:end-4),'_tunnel_',num2str(out_tunnel),output_format);
            filename=fullfile(output_dir,filename);
            y_toWrite=y(:,out_tunnel);
            audiowrite(filename,y_toWrite,Fs); %转写成.wav格式文件
            count=count+1;
            my_waitbar(L,count,bar);
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