function feature=featurevector(filename,framesize,step)                     %filename是读取的文件名,framesize是一帧的长度,step是滑动的点数
% filename='C:\Users\90317\Desktop\audio_work_dir\audiofile\negetive\蒋奉儒1_4_tunnel_1.wav';
% framesize=640;
% step=160;
overlap=framesize-step;                                                     %overlap是窗口重叠的点数
[y,fs] = audioread(filename);                                               %y是读取的双声道原始数据，fs是采样频率
y1=y(:,1);                                                                  %y1是获取的单声道原始数据
y1_pre=filter([1 -0.97],1,y1);                                              %y1_pre是预加重后的单声道数据
s1=enframe(y1,framesize,step);                                              %S1是单声道数据直接分帧
s1_hamming=enframe(y1,hamming(framesize),step);                             %s1_hamming是原始数据分帧加汉明窗
% s1_pre=enframe(y1_pre,framesize,step);                                    %s1_pre是单声道预加重数据直接分帧
s1_pre_hamming=enframe(y1_pre,hamming(framesize),step);                     %s1_pre_hamming是单声道预加重数据分帧加汉明窗
% yFFT=fft(y1);                                                             %yFFT是单声道原始数据进行FFT变换 
% yDFT=abs(yFFT);                                                           %yDFT是单声道原始数据的DFT变换结果 
yLen=length(y1);                                                            %yLen是单声道原始数据的点数  
k=0:yLen-1;  
Fk1=k.*fs/yLen;                                                             %Fk1是FFT变换每个点对应的频率
[nframe,~]=size(s1);                                                        %nframe是按帧长分帧后的帧数
T1=0.05;                                                                    %T1设置基音端点检测的参数
dctlength=12;                                                               %设置DCT变换的维数
number_of_mel_filter=24;                                                    %设置Mel滤波器个数
x=y1-mean(y1);                                                              
x=x/max(abs(x));                                                            %x是消除直流分量并幅值归一化的单声道原始数据
y_nod=enframe(x,framesize,step);                                            %y_nod消除直流分量并分帧的单声道归一化数据
% time = (0 : length(x)-1)/fs;                                                %time是时间坐标
% frameTime = frame2time(nframe, framesize, step, fs);                        %frame_time是各帧对应的时间坐标
[voiceseg,vosl,SF,~]=pitch_vad1(y_nod,nframe,T1);                           %voiceseg中有vosl个结构体，vosl是这一段的有语音段数
                                                                            %其中每个结构体的begin参数指明了此语音段的起始帧位置
                                                                            %每个结构体的end参数指出此语音段的结束帧位置
                                                                            %每个结构体的duration参数指出此语音段含有的帧数
                                                                            %SF(i)为1时第i个数据点在语音段中
E=zeros(1,nframe);  
for i=1:nframe
    E(i)=sum(s1_hamming(i,:).*s1_hamming(i,:));                             %E是短时能量
end

x=0;
for i=1:(nframe-1)
    t=abs(E(i)-E(i+1))/(nframe-1);
    x=x+t;
end
E_shimmer=x/mean(E);                                                        %E_shimmer是短时抖动能量

x1=0;x2=0;x3=0;x4=0;
for i=1:nframe
    t1=i*mean(E);
    t2=i*E(i);
    t3=i*i;
    t4=i;
    x1=x1+t1;
    x2=x2+t2;
    x3=x3+t3;
    x4=x4+t4;
end
x4=x4*x4/nframe;
z1=x2-x1;
z2=x3-x4;
E_Reg_coff=z1/z2;                                                           %E_Reg_coff短时能量的线性回归系数

x=0;
for i=1:nframe
    t=E(i)-(mean(E)-z1/z2*x4/nframe)-z1/z2*i;
    x=x+t^2/nframe;
end
E_Sqr_Err=x;                                                                %E_Sqr_Err短时能量的线性回归系数均方误差

[z2,f1,~]=specgram(y1,framesize,fs);
sn=20*log10(abs(z2)+eps);
sn1=sn+min(sn(:));
n=round(length(f1)*250/max(f1(:)));
Eratio=sum(sum(sn1(1:n,:)))/sum(sn1(:));                                    %Eratio是250Hz以下短时能量与全部短时能量的比

S=zeros(nframe,framesize);
for i=1:nframe  
	YPow=s1_hamming(i,:);
	YPowFft=fft(YPow);  
	YPowDFT=abs(YPowFft);  
    PowerY=YPowDFT.^2;     
    S(i,:)=PowerY;                                                          %S是功率谱
end  

Bright=zeros(nframe,1);
for i=1:nframe
    x1=0;
    x2=0;
    for j=1:framesize/2
        x1=x1+Fk1(j)*S(i,j);
        x2=x2+S(i,j);
    end
    if x2~=0
        Bright(i)=x1/x2;                                                   
    end
end
Bright(Bright==0)=[];                                                       %Bright是明亮度           

S_TEO=zeros(nframe,framesize/2+1);
for i=1:nframe  
    for j=2:framesize/2+1
        S_TEO(i,j-1)=S(i,j)^2-S(i,j+1)*S(i,j-1);                            %S_TEO是TEO变换后的功率谱值
    end
end  
S_TEO=abs(S_TEO);

NFD_Mel=zeros(nframe,dctlength);
dctcoef=zeros(dctlength,number_of_mel_filter);
bank=melbankm(number_of_mel_filter, framesize, fs, 0, 0.5, 't');
bank=full(bank);
bank=bank/max(bank(:));
for k=1:dctlength
    n=0:number_of_mel_filter-1;
    dctcoef(k,:)=cos((2*n+1)*k*pi/(2*24));
end
for i=1:nframe
    if ~all(S_TEO(i,:)==0)
        x_tolog=bank*(S_TEO(i,:)).';
        xd=log(x_tolog);
        NFD_Mel(i,:) = dctcoef*xd;                                          %计算特征NFD_Mel
    end
end
NFD_Mel(all(NFD_Mel==0,2),:) = [];

AF_Mel=zeros(nframe,dctlength);
AF=zeros(nframe,framesize/2+1);
for i=1:nframe
    for j=1:framesize/2+1
        AF(i,j)=Fk1(j)^2*S(i,j)^2;
    end
    if ~all(AF(i,:)==0)
        AF_Mel(i,:) = dctcoef*log(bank*AF(i,:).');                          %计算特征AF_Mel
    end
end
AF_Mel(all(AF_Mel==0,2),:) = [];

[mfccs,delta_mfcc,delta_delta_mfcc]= ...
    mfcc(y1_pre,fs,'WindowLength',framesize,...
    'OverlapLength',overlap , ...
    'NumCoeffs',dctlength,'LogEnergy','Ignore'); 
                                                                            %计算特征mfcc及其一阶和二阶差分

lmin=fix(fs/500);                                                           %基音周期的最小值
lmax=fix(fs/60);                                                            %基音周期的最大值              
period=zeros(1,nframe);                                                     %基音周期初始化
R=zeros(nframe,1);
Lc=zeros(nframe,1);
for k=1:nframe 
    if SF(k)==1                                                             %是否在有话帧中
        yb=y_nod(k,:).*hamming(framesize).';                                  %取来一帧数据加窗函数
        xx=fft(yb);                                                         %FFT
        a=2*log(abs(xx)+eps);                                               %取模值和对数
        b=ifft(a);                                                          %求取倒谱 
        [R(k),Lc(k)]=max(b(lmin:lmax));                                     %在lmin和lmax区间中寻找最大值
        period(k)=Lc(k)+lmin-1;                                             %给出基音周期
    end
end
F0=zeros(1,nframe);
for i=1:nframe
    if period(i) ~= 0
        F0(i)=fs/period(i);                                                 %计算基音频率
    end
end                                                                         

% FF=zeros(numel(F0)-1);
% k=1;
% for i=2:numel(F0)
%     if(F0(i)==F0(1))
%         continue;
%     end
%     FF(k)= F0(i);
%     k=k+1;
% end
FF=F0;
if ~all(FF==0)
    FF(FF==0)=[];
end

x=0;
for i=1:(nframe-1)
    t=abs(F0(i)-F0(i+1));
    x=x+t;
end
if ~mean(F0)*(nframe-1)==0
    F_Jitter1=100*x/(mean(F0)*(nframe-1));                                      %计算一阶基音频率抖动
else
    F_Jitter1=0;
end

x=0;
for i=2:(nframe-1)
    t=abs(2*F0(i)-F0(i+1)-F0(i-1));
    x=x+t;
end
if ~(mean(F0)*(nframe-2))==0
    F_Jitter2=100*x/(mean(F0)*(nframe-2));                                      %计算二阶基音频率抖动
else
    F_Jitter2=0;
end

nf=1;
dF=zeros(nframe,1);
for i=1:(nframe-1)
    if(F0(i)*F0(i+1)~=0)
        dF(nf)=F0(i)-F0(i+1);                                               %计算浊音帧差分基音
    end
    nf=nf+1;
end

interval=zeros(1,vosl-1);
for i=2:vosl
    interval(i-1)=voiceseg(i).begin-voiceseg(i-1).begin;                    %计算间隔
end
if isempty(interval)
    interval=0;
end

duration=zeros(vosl,1);
for i=1:vosl
    duration(i)=voiceseg(i).duration;                                       %提取每个有话段长度
end
if isempty(duration)
    duration=0;
end

Fm1=zeros(vosl,1);
Fm2=zeros(vosl,1);
Fm3=zeros(vosl,1);
df=fs/512;
if ~(length(voiceseg)==1 && voiceseg.duration(1)==0)
    for i=1:vosl
            u=s1_pre_hamming(voiceseg(i).begin:voiceseg(i).end,:);
            u=u(:);
            a=lpc(u,12);
            U=lpcar2pf(a,255);
            [Val,Loc]=findpeaks(U);
            ll=length(Loc);
            FPeak=zeros(ll);
            Bw=zeros(ll);
            for k=1 : ll
                m=Loc(k);                                                       %设置m-1,m和m+1
                m1=m-1; m2=m+1;
                p=Val(k);                                                       %设置P(m-1),P(m)和P(m+1)
                p1=U(m1); p2=U(m2);
                aa=(p1+p2)/2-p;                                                 %按式(9-3-4)计算
                bb=(p2-p1)/2;
                cc=p;
                dm=-bb/2/aa;                                                    %按式(9-3-6)计算
                pp=-bb*bb/4/aa+cc;                                              %按式(9-3-8)计算
                m_new=m+dm;
                bf=-sqrt(bb*bb-4*aa*(cc-pp/2))/aa;                              %按式(9-3-13)计算
                FPeak(k)=(m_new-1)*df;                                          %按式(9-3-7)计算
                Bw(k)=bf*df;                                                    %按式(9-3-14)计算
            end
            smooth(FPeak,Bw);
            if ll>=3                                                 %此处存疑
                Fm1(i)=FPeak(1);
                Fm2(i)=FPeak(2);
                Fm3(i)=FPeak(3);
            end
    end
end
if ~all(Fm1==0)
    Fm1(Fm1==0)=[];
end
if ~all(Fm2==0)
    Fm2(Fm2==0)=[];
end
if ~all(Fm3==0)
    Fm3(Fm3==0)=[];
end

x1=0;
x2=0;
x3=0;
for i=1:(numel(Fm1)-1)
    t1=abs(Fm1(i)-Fm1(i+1));
    t2=abs(Fm2(i)-Fm2(i+1));
    t3=abs(Fm3(i)-Fm3(i+1));
    x1=x1+t1;x2=x2+t2;x3=x3+t3;
end

if ~((mean(Fm1)*(numel(Fm1)-1))==0)
    Fm1_Jitter1=100*x1/(mean(Fm1)*(numel(Fm1)-1));                              %前三个共振峰的一阶抖动
else
    Fm1_Jitter1=0;
end
if ~((mean(Fm2)*(numel(Fm1)-1))==0)
    Fm2_Jitter1=100*x2/(mean(Fm2)*(numel(Fm1)-1));
else
    Fm2_Jitter1=0;
end
if ~((mean(Fm3)*(numel(Fm1)-1))==0)
    Fm3_Jitter1=100*x3/(mean(Fm3)*(numel(Fm1)-1));
else
    Fm3_Jitter1=0;
end
if ~((Fm2-Fm1)==0)
    Fm2R=Fm2./(Fm2-Fm1);
else
    Fm2R=0;
end

nFm=[max(Fm1);min(Fm1);mean(Fm1);var(Fm1);Fm1_Jitter1;max(Fm2);min(Fm2);...
    mean(Fm2);var(Fm2);Fm2_Jitter1;max(Fm3);min(Fm3);mean(Fm3);var(Fm3);...
    Fm3_Jitter1;max(Fm2R);min(Fm2R);mean(Fm2R)];                            %提取共振峰相关特征

feature(1:7,1)=[max(E);min(E);mean(E);var(E);E_shimmer;E_Reg_coff;E_Sqr_Err];
feature(8,1)=Eratio;
feature(9:14,1)=[max(F0);min(FF);mean(F0);var(F0);F_Jitter1;F_Jitter2];     %基音频率相关特征
feature(15:18,1)=[max(dF);min(dF);mean(dF);var(dF)];                        %浊音帧差分基音
feature(19:(size(nFm)+18),1)=nFm;                                           %20-37
for i=1:size(mfccs,2)
    feature(37+4*(i-1):37+4*i-1,1)=[max(mfccs(:,i));min(mfccs(:,i));...
        mean(mfccs(:,i));var(mfccs(:,i))];                                  %mel倒谱系数及其一阶差分相关特征
end
for i=1:size(delta_mfcc,2)
    feature(85+4*(i-1):85+4*i-1,1)=[max(delta_mfcc(:,i));...
        min(delta_mfcc(:,i));mean(delta_mfcc(:,i));...
        var(delta_mfcc(:,i))];                        
end
for i=1:size(delta_delta_mfcc,2)
    feature(133+4*(i-1):133+4*i-1,1)=[max(delta_delta_mfcc(:,i));...
        min(delta_delta_mfcc(:,i));mean(delta_delta_mfcc(:,i));...
        var(delta_delta_mfcc(:,i))];                        
end
for i=1:size(NFD_Mel,2)
    feature(181+4*(i-1):181+4*i-1,1)=[max(NFD_Mel(:,i));...
        min(NFD_Mel(:,i));mean(NFD_Mel(:,i));...
        var(NFD_Mel(:,i))];                        
end
for i=1:size(AF_Mel,2)
    feature(229+4*(i-1):229+4*i-1,1)=[max(AF_Mel(:,i));...
        min(AF_Mel(:,i));mean(AF_Mel(:,i));...
        var(AF_Mel(:,i))];                        
end
feature(277:280,1)=[max(interval);min(interval);...
    mean(interval);var(interval)];                       
feature(281:284,1)=[max(duration);min(duration);...
    mean(duration);var(duration)]; 
feature(285:288,1)=[max(Bright);min(Bright);...
    mean(Bright);var(Bright)]; 

end






% function frameTime=frame2time(frameNum,framelen,inc,fs)
%     frameTime=(((1:frameNum)-1)*inc+framelen/2)/fs;
% end

function [voiceseg,vosl,SF,Ef]=pitch_vad1(y,fn,T1,miniL)
if nargin<4, miniL=10; end
if size(y,2)~=fn, y=y'; end                                                 %把y转换为每列数据表示一帧语音信号
wlen=size(y,1);                                                             %取得帧长
Esum=zeros(fn,1);
H=zeros(fn,1);
for i=1:fn
    Sp = abs(fft(y(:,i)));                                                  %FFT取幅值
    Sp = Sp(1:wlen/2+1);	                                                %只取正频率部分
    Esum(i) = sum(Sp.*Sp);                                                  %计算能量值
    prob = Sp/(sum(Sp));	                                                %计算概率
    H(i) = -sum(prob.*log(prob+eps));                                       %求谱熵值
end
hindex=H<0.1;
H(hindex)=max(H);
Ef=sqrt(1 + abs(Esum./H));                                                  %计算能熵比
Ef=Ef/max(Ef);                                                              %归一化
zindex=find(Ef>=T1);                                                        %寻找Ef中大于T1的部分
zseg=findSegment(zindex);                                                   %给出端点检测各段的信息
zsl=length(zseg);                                                           %给出段数
j=0;
SF=zeros(1,fn);
voiceseg = struct('begin',zeros(zsl,1),'end',zeros(zsl,1), ...
    'duration',zeros(zsl,1));
for k=1 : zsl                                                               %在大于T1中剔除小于miniL的部分
    if zseg(k).duration>=miniL
        j=j+1;
        in1=zseg(k).begin;
        in2=zseg(k).end;
        voiceseg(j).begin=in1;
        voiceseg(j).end=in2;
        voiceseg(j).duration=zseg(k).duration;
        SF(in1:in2)=1;                                                      %设置SF
    end
end
vosl=length(voiceseg);                                                      %有话段的段数 
end

function soundSegment=findSegment(express)
if express(1)==0
    voicedIndex=find(express);                                              %寻找express中为1的位置
else
    voicedIndex=express;
end
soundSegment = [];
k = 1;
soundSegment(k).begin = voicedIndex(1);                                     %设置第一组有话段的起始位置
for i=1:length(voicedIndex)-1
	if voicedIndex(i+1)-voicedIndex(i)>1                                    %本组有话段结束
		soundSegment(k).end = voicedIndex(i);                               %设置本组有话段的结束位置
		soundSegment(k+1).begin = voicedIndex(i+1);                         %设置下一组有话段的起始位置  
		k = k+1;
	end
end
soundSegment(k).end = voicedIndex(end);                                     %最后一组有话段的结束位置
                                                                            %计算每组有话段的长度
for i=1:k
    soundSegment(i).duration=soundSegment(i).end-soundSegment(i).begin+1;
end
end