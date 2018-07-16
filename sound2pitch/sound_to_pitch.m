function [pitch] = sound_to_pitch(arr,timeStep,minimumPitch,maximumPitch,ccc,ttt,datafile)
%UNTITLED 此处显示有关此函数的摘要
%声音转pitch
%arr：音频数据，44100，单声道，16bit
%timeStep：时间戳
%minimumPitch：最小pitch
%maximumPitch：最大pitch，不会超过 44100/2

debuglog=1; %是否打log到文件
debugplot = 0;%是否画图

%Step.1
% af = fft(arr);%fft
% l=length(af);
% start1 = floor(l *0.10);
% saf=fftshift(af);%翻转，使0hz为中心，两端是高频
%   
% osaf=zeros(1,length(saf));
% osaf(start1+1:l-start1)=saf(start1+1:l-start1);%只取0-95%部分，其余部分清零
%   
% arr1=ifft(ifftshift(osaf));%逆变换
% arr = arr1;
%以上放低通滤波器,过滤掉 Nyquistfrequency 的 95%-100%的点
%Step.1 end

%000一些参数
dx = 1.0 / 44100;       %每个采样点时间
nx = length(arr);       %采样数
duration = dx * nx;     %总时长
x1 = 0.5 / 44100;       %半个采样点时间
periodsPerWindow=3.0;   %每个window三个周期，至少
%minimumPitch = 75;      %最小pitch 上一层传过来
%maximumPitch = 44100/4; %最大pitch 上一层传过来

brent_depth = 3;
silenceThreshold = 0.03;
voicingThreshold = 0.45;
octaveCost = 0.01;
octaveJumpCost = 0.35;
voicedUnvoicedCost = 0.14;
            

nsamp_window1 = (2048/2 + 1)*2;
dt_window1 = nsamp_window1 * dx;
minimumPitch = periodsPerWindow / dt_window1 ;

nsamp_period = floor (1.0 / dx / minimumPitch);%最大周期采样点数
halfnsamp_period = nsamp_period / 2 + 1;        %用于计算两个方向上的平均值和局部峰值
if(maximumPitch > 0.5 / dx)
    maximumPitch = 0.5 / dx;     %最大pitch不能超过采样率的一半
end

maxnCandidates = floor (maximumPitch / minimumPitch);    %最大候选数量,
%   调试
maxnCandidates = 15;
%   调试
if(timeStep == 0)
    timeStep = periodsPerWindow / minimumPitch / 4.0;   % e.g. 3 periods, 75 Hz: 10 milliseconds 决定帧数
end

dt_window = periodsPerWindow / minimumPitch;    %窗长度，根据最大周期计算，要包含三个最大周期
nsamp_window = floor (dt_window / dx);      %每个窗采样点
halfnsamp_window = nsamp_window / 2 - 1;    %半个窗采样点，为啥-1？
nsamp_window = halfnsamp_window * 2;        %一个窗采样点
interpolation_depth = 0.5;                  %hanning
brent_ixmax = floor (nsamp_window * interpolation_depth);   %半个窗？
minimumLag = floor (1.0 / dx / maximumPitch);    %最小时间
if (minimumLag < 2) 
    minimumLag = 2; %maxmumPitch最大值不会超过采样率的一半，所以为2
end         
minimumLag
maximumLag = floor (nsamp_window / periodsPerWindow) + 2;   %最大时间延时
if (maximumLag > nsamp_window) 
    maximumLag = nsamp_window;  %最大时间延时不会大于窗的采样点数，
end 
%Determine the number of frames. Fit as many frames as possible symmetrically in the total duration.
%We do this even for the forward cross-correlation method,because that allows us to compare the two methods.
if (dt_window > duration) 
    return; %err 窗时间小于文件时长
end 
numberOfFrames = floor ((duration - dt_window) / timeStep) + 1; %帧数，可以看到是富余了一个窗长度
%   调试
%numberOfFrames = 5;
%   调试
if(numberOfFrames <1)
    return ;
end %err 

thyDuration = numberOfFrames * timeStep;
firstTime = 0.5 * duration - 0.5 * thyDuration + 0.5 * timeStep;
pitch = zeros(3,numberOfFrames);    %  3 * numberOfFrames 矩阵 ， 用以保存pitch信息


%Step.2    计算全局峰值
globalPeak = 0.0;
sum = 0.0;
for i=1:nx
    sum = sum + arr(i);
end
mean = sum/nx;
for i=1:nx
    value = abs(arr(i)-mean);
    if(value > globalPeak)
       globalPeak = value; 
    end
end

if(globalPeak==0)
   return; 
end
%Step.2  end


%窗的自相关
window=zeros(1,nsamp_window);
for j=1:nsamp_window
	window(j) = 0.5 - 0.5 * cos(j * 2 * pi / (nsamp_window + 1));
end				
windowR=xcorr(window);
windowR=windowR(length(windowR)-length(window)+1:length(windowR));
windowR = windowR / windowR(1);
if debugplot ~= 0
    figure('NumberTitle', 'off', 'Name', '窗函数及其归一化自相关函数');
    title('窗函数及其归一化自相关函数');
    xlabel('samples数');  %x轴
    ylabel('值float');%y轴
    hold on
    plot(window,'-r'); %窗图形
    hold on
    plot(windowR,'-b'); %窗自相关图形
end
%Step.3   对每一帧进行处理
%Step 3.1 计算窗长度<->MinimumPitch,至少包含三个周期

data=zeros(maxnCandidates*2,numberOfFrames);% 一帧有maxnCandidates个候选，每个候选包含2个元素，freq、strength
Ofreq=1;
Ostrength=2;
nCandidate=zeros(1,numberOfFrames);         % 用于记录每一帧计算出来的满足要求候选者要求的数量
imax=zeros(1,numberOfFrames);               %
intensity = zeros(1,numberOfFrames);        %记录每一帧的intensity （localpeak / globalpeak）

%打开文件log
if debuglog ~= 0
    diary('testlog3.txt');
    diary on;
end

% 打印关键参数
if debuglog ~= 0

    str=['maxnCandidates=' num2str(maxnCandidates)];
    disp(str);

    str=['dt=' num2str(timeStep)];
    disp(str);

    str=['duration=' num2str(duration)];
    disp(str);

    str=['brent_depth=' num2str(brent_depth)];
    disp(str);

    str=['interpolation_depth=' num2str(interpolation_depth)];
    disp(str);

    str=['nsamp_period=' num2str(nsamp_period)];
    disp(str);

    str=['halfnsamp_period=' num2str(halfnsamp_period)];
    disp(str);

    str=['ceiling=' num2str(maximumPitch)];
    disp(str);

    str=['dt_window=' num2str(dt_window)];
    disp(str);

    str=['nsamp_window=' num2str(nsamp_window)];
    disp(str);

    str=['halfnsamp_window=' num2str(halfnsamp_window)];
    disp(str);

    str=['minimumLag=' num2str(minimumLag)];
    disp(str);

    str=['maximumLag=' num2str(maximumLag)];
    disp(str);

    str=['numberOfFrames=' num2str(numberOfFrames)];
    disp(str);

    str=['t1=' num2str(firstTime)];
    disp(str);

    str=['gloable=' num2str(globalPeak)];
    disp(str);

    str=['brent_ixmax=' num2str(brent_ixmax)];
    disp(str);

    %打印window和windowR

    disp('printf window and windowR');
    % for i=1:nsamp_window
    %     str=['window[' num2str(i) ']=' num2str(window(i))];
    %     disp(str);
    % end
    % 
    % for i=1:nsamp_window
    %     str=['windowR[' num2str(i) ']=' num2str(windowR(i))];
    %     disp(str);   
    % end

end

% 逐帧进行处理
for i=1:numberOfFrames
    %计算局部平均值，是按照一个周期采样点数nsamp_period计算的
    localMean = 0.0;
    localPeak = 0.0;
%     x = firstTime + (i - 1) * timeStep; 
    x = (i+1) * timeStep; 
%     leftSample = floor(((x - x1) / dx + 1.0));
    leftSample = floor(((x - x1) / dx + 1.0))-1;

    rightSample = leftSample +1;
    startSample = rightSample - nsamp_period;
	endSample = leftSample + nsamp_period;
    
    for m=startSample:endSample
        localMean = localMean + arr(m);
    end
    localMean = localMean / (2 * nsamp_period);

    if debuglog ~= 0
%         str=['startSample=' num2str(startSample) 'endSample=' num2str(endSample) 'localMean=' num2str(localMean)];
%         disp(str);
    end
    %Step 3.2 减去平均值按照一个窗计算
    startSample = rightSample - halfnsamp_window;
	endSample = leftSample + halfnsamp_window;
    
    frame=arr(startSample:endSample);
    
    if debuglog ~= 0
%          disp('printf frame');
%          for xx=1:nsamp_window
%              str=['frame[' num2str(xx) ']=' num2str(frame(xx))];
%              disp(str);
%          end
    end
    
    for j=1:nsamp_window
        frame(j) = (frame(j) - localMean);
    end

    if debuglog ~= 0
%         str=['startSample=' num2str(startSample) 'endSample=' num2str(endSample) 'localMean=' num2str(localMean)];
%         disp(str);
    end
    
    %Step 3.4 乘以窗函数
    for j=1:nsamp_window
        frame(j) = frame(j)*window(j);
    end
    
    %Step 3.3 计算unvoiced candidate，这里是计算localpeak
    startSample = halfnsamp_window + 1 - halfnsamp_period;
    endSample = halfnsamp_window + halfnsamp_period;
    
    if (startSample < 1)
        startSample = 1;
    end
    
	if (endSample > nsamp_window)
        endSample = nsamp_window;
	end
    
    if debuglog ~= 0
%         str=['startSample=' num2str(startSample) 'endSample=' num2str(endSample)];
%         disp(str);
    end
    
    if debuglog ~= 0
    %     disp('printf frames');
    %     for xx=startSample:endSample
    %         str=['frame[' num2str(xx) ']=' num2str(frame(xx))];
    %         disp(str);
    %     end
    end
    
    for j = 1:nsamp_window
        value = abs (frame(j));
        if (value > localPeak)
            localPeak = value;
        end
    end
    
    if(localPeak > 1.0)
        intensity(i) = 1.0;
    else
        intensity(i) = localPeak/1.0;
    end
    
    if debuglog ~= 0
%         str=['localPeak=' num2str(localPeak) 'intensity=' num2str(intensity(i))];
%         disp(str);
    %     disp('printf frames-mean');
    %     for xx=startSample:endSample
    %         str=['frame[' num2str(xx) ']=' num2str(frame(xx))];
    %         disp(str);
    %     end
    end

    %Step 3.5 - 3.9 计算ra(τ)
    %Step 3.10 Divide by the autocorrelation of the window, which was computed once with steps 3.5 through 3.9 (equation 9). This gives a sampled version of rx(τ).
    r = xcorr(frame);
    rx= r(length(r)-length(frame)+1:length(r));
    for j=2:nsamp_window
        rx(j) = rx(j) / (rx(1) * windowR(j));
    end
    rx(1) = 1;
    rx=rx(2:nsamp_window);
    
    if debuglog ~=0 
        %打印rx
%          for j=1:brent_ixmax
%                   str=['rx(' num2str(j) ']=' num2str(rx(j))];
%                   disp(str);
%          end
    end
     
    %Candidate计算
    nCandidate(i) = 1;  %第i帧拥有的备选数量设置为1，因为必有一个是unvoiced

    data(Ofreq,i)=0;            %第i帧的，Ofreq=1也就是freq，置为0
    data(Ostrength,i)=0;        %第i帧的，Ostrength=2也就是strength，置为0
    if(localPeak == 0)  %全是静音数据，保留unvoiced备选，直接退出
        return;
    end
    %找到这个帧的自相关的最强的极大值，注册他们作为候选
    result = 0.0;
    strengthOfMaximum = 0.0;
    %以半个窗或者最大时间延时为界限
    if(brent_ixmax < maximumLag)
        limit = brent_ixmax;
    else
        limit = maximumLag;
    end
    %从2起
    n_samples=zeros(1,maxnCandidates);
    n_r=zeros(1,maxnCandidates);

    if debuglog ~=0 
%          disp(num2str(nsamp_window));
%          str=['xx=' num2str(i*timeStep)];
%          disp(str);
    end
    for j=minimumLag:limit-1   
            if( rx(j)>rx(j-1)&&rx(j)>rx(j+1)) %% maximum?
                %11Find the strongest maxima of the correlation of this frame, and register them as candidates.
                
                if(rx(j) <= 0.5*voicingThreshold )%可下调rx阀值由0.5->0.4,增加candidate
                    continue;
                end
                
                place = 0;
                dr = 0.5 * (rx(j+1) - rx(j-1));
                d2r = 2.0 * rx(j) - rx(j-1) - rx(j+1);
                frequencyOfMaximum = 1.0 / dx / (j + dr / d2r);
                
                offset = -brent_ixmax - 1;
    
                y=zeros(1,brent_ixmax*2+1);
                y(brent_ixmax+1) = 1.0;%0点
                for m=1:brent_ixmax
                y(m+brent_ixmax+1) = rx(m);
                y(-m+brent_ixmax+1) = rx(m);
                end
                
                strengthOfMaximum = NUM_interpolate_sinc(y,brent_ixmax - offset,1 / dx / frequencyOfMaximum - offset,30);
                if (strengthOfMaximum > 1.0) 
                    strengthOfMaximum = 1.0 / strengthOfMaximum;
                end
                
                if debuglog ~=0 
%                     str=['frequencyOfMaximum=' num2str(frequencyOfMaximum) 'strengthOfMaximum=' num2str(strengthOfMaximum) 'rx [' num2str(j) ']' num2str(rx(j)) ];
%                     disp(str);
                end

                if (nCandidate(i) < maxnCandidates) %is there still a free place?
                    nCandidate(i) = nCandidate(i) + 1;
                    place = nCandidate(i);
                else 
                    % 找到最弱的candidate的位置
                    weakest = 2.0;
                    for iweak = 2:maxnCandidates
                        % High frequencies are to be favoured */
                        % if we want to analyze a perfectly periodic signal correctly. */
                        localStrength = data((iweak-1)*2+Ostrength,i) - octaveCost * log2 (minimumPitch / data((iweak-1)*2+Ofreq,i) );
                        if (localStrength < weakest) 
                            weakest = localStrength;
                            place = iweak;
                        end
                    end
                    %如果这个比最弱的还弱，就什么也不做，置place为0，否则，记录place
                    if (strengthOfMaximum - octaveCost * log2 (minimumPitch / frequencyOfMaximum) <= weakest)
                        place = 0;
                    end
                end
                
                if (place) %  当前找到的替换最弱的canditate
                    data((place-1)*2+Ofreq,i) = frequencyOfMaximum;  %i帧的第place个candidate
                    data((place-1)*2+Ostrength,i) = strengthOfMaximum;
                    imax(place)= i;%这个记录每个candidate对应的帧数
                    n_samples(place)=j;
                    n_r(place)=rx(j);
                end
                
                % Second pass: for extra precision, maximize sin(x)/x interpolation ('sinc').
                %这里对每个candidate的频率和强度有个修正，这步还不知道怎么做。代码还原比较复杂
%                 	for ii = 2:nCandidate(i)
%                         if ( data((ii-1)*2+1,i)> 0.0/timeStep) 
%                             offset = - brent_ixmax - 1;
%                             y=zeros(1,brent_ixmax*2+1);
%                             y(brent_ixmax+1) = 1.0;%0点
%                             for m=1:brent_ixmax
%                             y(m+brent_ixmax+1) = rx(m+1);
%                             y(-m+brent_ixmax+1) = rx(m+1);
%                             end
%                             if(data((ii-1)*2+1,i) > 0.3 / timeStep )
%                                 [ymid,xmid] = NUMimproveMaximum (y, brent_ixmax - offset, imax(ii) - offset,4);
%                             else
%                                 [ymid,xmid] = NUMimproveMaximum (y, brent_ixmax - offset, imax(ii) - offset,brent_depth);
%                             end
%                             xmid = xmid + offset;
%                             data((ii-1)*2+1,i) = 1.0 / timeStep / xmid;
%                             if (ymid > 1.0) 
%                                 ymid = 1.0 / ymid;
%                             end
%                             data((ii-1)*2+2,i) = ymid;
              
%                         end
%                     end
                %11 End
            end
            
%             %这里调整place存放顺序，把强度最大的，放前面
%             best = 0.0;
%             ad = 2;
%             for ii=2 : nCandidate(i)
%                 if (data((ii-1)*2+Ostrength,i)>best)
%                     best = data((ii-1)*2+Ostrength,i);
%                     ad = ii;
%                 end
%             end
%             %交换
%             f = data((ad-1)*2+Ofreq,i);
%             l = data((ad-1)*2+Ostrength,i);
%             data((ad-1)*2+Ofreq,i)      = data(2+Ofreq,i);
%             data((ad-1)*2+Ostrength,i)  = data(2+Ostrength,i);
%             data(2+Ofreq,i) = f;
%             data(2+Ostrength,i) = l;
             
    end %forj=2：limit-1
end

% for i=1:length(data(2+Ofreq,:))
%     fprintf(datafile, '%f\t%f\n',((ccc-1)*ttt/44100+i*timeStep),data(2+Ofreq,i));
% end
% return;

if debuglog ~=0
%      for xx=1:numberOfFrames  
%          str=['xx=' num2str(xx*timeStep)];
%          disp(str);
%          for xxx=1:getMaxnCandidates(nCandidate)
%              str=['xxx=' num2str(data((xxx-1)*2+Ofreq,xx))];
%              disp(str);
%          end
%      end
     diary off;
end

%动态规划部分
%找到最大的candidate数
maxnc = getMaxnCandidates(nCandidate);
maxnc
ceiling2 = maximumPitch;
timeStepCorrection = 0.01/timeStep;%这

octaveJumpCost = octaveJumpCost * timeStepCorrection;
voicedUnvoicedCost = voicedUnvoicedCost * timeStepCorrection;

delta = zeros(numberOfFrames, maxnc);
psi = zeros(numberOfFrames, maxnc);
unvoicedStrength = 0.0;

for i=1:numberOfFrames
    %R
    if(silenceThreshold <= 0)
        unvoicedStrength =  0.0;
    else
        unvoicedStrength = 2.0 - intensity(i) / (silenceThreshold / (1.0 + voicingThreshold));
    end
    if(unvoicedStrength > 0.0)
        unvoicedStrength = voicingThreshold + unvoicedStrength;
    else
        unvoicedStrength = voicingThreshold + 0.0;
    end
    %voice和非voice区分开来
	for icand=1:nCandidate(i)
		voiceless = ~(data((icand-1)*2+Ofreq,i) > 0.0 && data((icand-1)*2+Ofreq,i) < ceiling2);
        if(voiceless)
            delta(i,icand) =  unvoicedStrength;
        else
            delta(i,icand) = data((icand-1)*2+Ostrength,i) - octaveCost * log2(maximumPitch / data((icand-1)*2+Ofreq,i));
        end
    end
end


%寻找一个最优的路径
transitionCost=0.0;
for i=2:numberOfFrames
    for icand2=1:nCandidate(i)          %i帧候选数量
        f2 = data((icand2-1)*2+Ofreq,i);    %i帧第icand2候选者的频率
        place = 0;
        maximum = -1e30;
        %计算前一帧每个candidate到当前icand2的cost
        for icand1=1:nCandidate(i-1)    %i前一帧候选数量
            f1 = data((icand1-1)*2+Ofreq,i-1);  %i前一帧icand1候选者的频率
            previousVoiceless = ~(f1 > 0.0 && f1 < ceiling2);
            currentVoiceless = ~(f2 > 0.0 && f2 < ceiling2);
            if(currentVoiceless)
                if(previousVoiceless)
                    transitionCost = 0.0;   %当前非语音，之前非语音
                else
                    transitionCost = voicedUnvoicedCost;    %当前非语音，之前语音
                end
            else
                if(previousVoiceless)
                    transitionCost = voicedUnvoicedCost;    %当前语音，之前非语音
                else
                    transitionCost = octaveJumpCost * abs(log2 (f1 / f2));    %当前语音，之前语音
                end
            end
            
            value = delta(i-1,icand1) - transitionCost + delta(i,icand2);
            if(value > maximum)
                maximum = value;
                place = icand1;
            elseif(value == maximum)
                
            
            end
        end %计算前一帧每个candidate到当前icand2的cost   结束
        delta(i,icand2) = maximum; %更新delta地i帧第icand2候选者的花费为preallcost
        psi(i,icand2) = place;      %psi记录前一帧最小cost位置
    end
end
%找到最优路径
place = 1;
maximum = delta(numberOfFrames,place);
for icand=2:nCandidate(numberOfFrames)
    if(delta(numberOfFrames,icand) > maximum)
        place = icand;
        maximum = delta(numberOfFrames,place);
    end
end
%follow the path backwards
for ll=1:numberOfFrames
    iframe = numberOfFrames-ll+1;
    %将place位置和第一个candidate进行交换, 
    f = data((place-1)*2+Ofreq,iframe);
    l = data((place-1)*2+Ostrength,iframe);
    data((place-1)*2+Ofreq,iframe)      = data(Ofreq,iframe);
    data((place-1)*2+Ostrength,iframe)  = data(Ostrength,iframe);
    data(Ofreq,iframe) = f;
    data(Ostrength,iframe) = l;
    place = psi(iframe,place);
end

if debugplot ~= 0
    figure('NumberTitle', 'off', 'Name', '最优路径');
    hold on
     title('最优路径');
     xlabel('帧数');  %x轴
     ylabel('频率');%y轴
    plot(data(Ofreq,:),'b*');
end

for i=1:length(data(Ofreq,:))
    fprintf(datafile, '%f\t%f\n',((ccc-1)*ttt/44100+i*timeStep),data(Ofreq,i));
end

end

