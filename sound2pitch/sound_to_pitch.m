function [pitch] = sound_to_pitch(arr,timeStep,minimumPitch,maximumPitch,ccc,ttt,datafile)
%UNTITLED �˴���ʾ�йش˺�����ժҪ
%����תpitch
%arr����Ƶ���ݣ�44100����������16bit
%timeStep��ʱ���
%minimumPitch����Сpitch
%maximumPitch�����pitch�����ᳬ�� 44100/2

debuglog=1; %�Ƿ��log���ļ�
debugplot = 0;%�Ƿ�ͼ

%Step.1
% af = fft(arr);%fft
% l=length(af);
% start1 = floor(l *0.10);
% saf=fftshift(af);%��ת��ʹ0hzΪ���ģ������Ǹ�Ƶ
%   
% osaf=zeros(1,length(saf));
% osaf(start1+1:l-start1)=saf(start1+1:l-start1);%ֻȡ0-95%���֣����ಿ������
%   
% arr1=ifft(ifftshift(osaf));%��任
% arr = arr1;
%���Ϸŵ�ͨ�˲���,���˵� Nyquistfrequency �� 95%-100%�ĵ�
%Step.1 end

%000һЩ����
dx = 1.0 / 44100;       %ÿ��������ʱ��
nx = length(arr);       %������
duration = dx * nx;     %��ʱ��
x1 = 0.5 / 44100;       %���������ʱ��
periodsPerWindow=3.0;   %ÿ��window�������ڣ�����
%minimumPitch = 75;      %��Сpitch ��һ�㴫����
%maximumPitch = 44100/4; %���pitch ��һ�㴫����

brent_depth = 3;
silenceThreshold = 0.03;
voicingThreshold = 0.45;
octaveCost = 0.01;
octaveJumpCost = 0.35;
voicedUnvoicedCost = 0.14;
            

nsamp_window1 = (2048/2 + 1)*2;
dt_window1 = nsamp_window1 * dx;
minimumPitch = periodsPerWindow / dt_window1 ;

nsamp_period = floor (1.0 / dx / minimumPitch);%������ڲ�������
halfnsamp_period = nsamp_period / 2 + 1;        %���ڼ������������ϵ�ƽ��ֵ�;ֲ���ֵ
if(maximumPitch > 0.5 / dx)
    maximumPitch = 0.5 / dx;     %���pitch���ܳ��������ʵ�һ��
end

maxnCandidates = floor (maximumPitch / minimumPitch);    %����ѡ����,
%   ����
maxnCandidates = 15;
%   ����
if(timeStep == 0)
    timeStep = periodsPerWindow / minimumPitch / 4.0;   % e.g. 3 periods, 75 Hz: 10 milliseconds ����֡��
end

dt_window = periodsPerWindow / minimumPitch;    %�����ȣ�����������ڼ��㣬Ҫ���������������
nsamp_window = floor (dt_window / dx);      %ÿ����������
halfnsamp_window = nsamp_window / 2 - 1;    %����������㣬Ϊɶ-1��
nsamp_window = halfnsamp_window * 2;        %һ����������
interpolation_depth = 0.5;                  %hanning
brent_ixmax = floor (nsamp_window * interpolation_depth);   %�������
minimumLag = floor (1.0 / dx / maximumPitch);    %��Сʱ��
if (minimumLag < 2) 
    minimumLag = 2; %maxmumPitch���ֵ���ᳬ�������ʵ�һ�룬����Ϊ2
end         
minimumLag
maximumLag = floor (nsamp_window / periodsPerWindow) + 2;   %���ʱ����ʱ
if (maximumLag > nsamp_window) 
    maximumLag = nsamp_window;  %���ʱ����ʱ������ڴ��Ĳ���������
end 
%Determine the number of frames. Fit as many frames as possible symmetrically in the total duration.
%We do this even for the forward cross-correlation method,because that allows us to compare the two methods.
if (dt_window > duration) 
    return; %err ��ʱ��С���ļ�ʱ��
end 
numberOfFrames = floor ((duration - dt_window) / timeStep) + 1; %֡�������Կ����Ǹ�����һ��������
%   ����
%numberOfFrames = 5;
%   ����
if(numberOfFrames <1)
    return ;
end %err 

thyDuration = numberOfFrames * timeStep;
firstTime = 0.5 * duration - 0.5 * thyDuration + 0.5 * timeStep;
pitch = zeros(3,numberOfFrames);    %  3 * numberOfFrames ���� �� ���Ա���pitch��Ϣ


%Step.2    ����ȫ�ַ�ֵ
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


%���������
window=zeros(1,nsamp_window);
for j=1:nsamp_window
	window(j) = 0.5 - 0.5 * cos(j * 2 * pi / (nsamp_window + 1));
end				
windowR=xcorr(window);
windowR=windowR(length(windowR)-length(window)+1:length(windowR));
windowR = windowR / windowR(1);
if debugplot ~= 0
    figure('NumberTitle', 'off', 'Name', '�����������һ������غ���');
    title('�����������һ������غ���');
    xlabel('samples��');  %x��
    ylabel('ֵfloat');%y��
    hold on
    plot(window,'-r'); %��ͼ��
    hold on
    plot(windowR,'-b'); %�������ͼ��
end
%Step.3   ��ÿһ֡���д���
%Step 3.1 ���㴰����<->MinimumPitch,���ٰ�����������

data=zeros(maxnCandidates*2,numberOfFrames);% һ֡��maxnCandidates����ѡ��ÿ����ѡ����2��Ԫ�أ�freq��strength
Ofreq=1;
Ostrength=2;
nCandidate=zeros(1,numberOfFrames);         % ���ڼ�¼ÿһ֡�������������Ҫ���ѡ��Ҫ�������
imax=zeros(1,numberOfFrames);               %
intensity = zeros(1,numberOfFrames);        %��¼ÿһ֡��intensity ��localpeak / globalpeak��

%���ļ�log
if debuglog ~= 0
    diary('testlog3.txt');
    diary on;
end

% ��ӡ�ؼ�����
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

    %��ӡwindow��windowR

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

% ��֡���д���
for i=1:numberOfFrames
    %����ֲ�ƽ��ֵ���ǰ���һ�����ڲ�������nsamp_period�����
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
    %Step 3.2 ��ȥƽ��ֵ����һ��������
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
    
    %Step 3.4 ���Դ�����
    for j=1:nsamp_window
        frame(j) = frame(j)*window(j);
    end
    
    %Step 3.3 ����unvoiced candidate�������Ǽ���localpeak
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

    %Step 3.5 - 3.9 ����ra(��)
    %Step 3.10 Divide by the autocorrelation of the window, which was computed once with steps 3.5 through 3.9 (equation 9). This gives a sampled version of rx(��).
    r = xcorr(frame);
    rx= r(length(r)-length(frame)+1:length(r));
    for j=2:nsamp_window
        rx(j) = rx(j) / (rx(1) * windowR(j));
    end
    rx(1) = 1;
    rx=rx(2:nsamp_window);
    
    if debuglog ~=0 
        %��ӡrx
%          for j=1:brent_ixmax
%                   str=['rx(' num2str(j) ']=' num2str(rx(j))];
%                   disp(str);
%          end
    end
     
    %Candidate����
    nCandidate(i) = 1;  %��i֡ӵ�еı�ѡ��������Ϊ1����Ϊ����һ����unvoiced

    data(Ofreq,i)=0;            %��i֡�ģ�Ofreq=1Ҳ����freq����Ϊ0
    data(Ostrength,i)=0;        %��i֡�ģ�Ostrength=2Ҳ����strength����Ϊ0
    if(localPeak == 0)  %ȫ�Ǿ������ݣ�����unvoiced��ѡ��ֱ���˳�
        return;
    end
    %�ҵ����֡������ص���ǿ�ļ���ֵ��ע��������Ϊ��ѡ
    result = 0.0;
    strengthOfMaximum = 0.0;
    %�԰�����������ʱ����ʱΪ����
    if(brent_ixmax < maximumLag)
        limit = brent_ixmax;
    else
        limit = maximumLag;
    end
    %��2��
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
                
                if(rx(j) <= 0.5*voicingThreshold )%���µ�rx��ֵ��0.5->0.4,����candidate
                    continue;
                end
                
                place = 0;
                dr = 0.5 * (rx(j+1) - rx(j-1));
                d2r = 2.0 * rx(j) - rx(j-1) - rx(j+1);
                frequencyOfMaximum = 1.0 / dx / (j + dr / d2r);
                
                offset = -brent_ixmax - 1;
    
                y=zeros(1,brent_ixmax*2+1);
                y(brent_ixmax+1) = 1.0;%0��
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
                    % �ҵ�������candidate��λ��
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
                    %�������������Ļ�������ʲôҲ��������placeΪ0�����򣬼�¼place
                    if (strengthOfMaximum - octaveCost * log2 (minimumPitch / frequencyOfMaximum) <= weakest)
                        place = 0;
                    end
                end
                
                if (place) %  ��ǰ�ҵ����滻������canditate
                    data((place-1)*2+Ofreq,i) = frequencyOfMaximum;  %i֡�ĵ�place��candidate
                    data((place-1)*2+Ostrength,i) = strengthOfMaximum;
                    imax(place)= i;%�����¼ÿ��candidate��Ӧ��֡��
                    n_samples(place)=j;
                    n_r(place)=rx(j);
                end
                
                % Second pass: for extra precision, maximize sin(x)/x interpolation ('sinc').
                %�����ÿ��candidate��Ƶ�ʺ�ǿ���и��������ⲽ����֪����ô�������뻹ԭ�Ƚϸ���
%                 	for ii = 2:nCandidate(i)
%                         if ( data((ii-1)*2+1,i)> 0.0/timeStep) 
%                             offset = - brent_ixmax - 1;
%                             y=zeros(1,brent_ixmax*2+1);
%                             y(brent_ixmax+1) = 1.0;%0��
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
            
%             %�������place���˳�򣬰�ǿ�����ģ���ǰ��
%             best = 0.0;
%             ad = 2;
%             for ii=2 : nCandidate(i)
%                 if (data((ii-1)*2+Ostrength,i)>best)
%                     best = data((ii-1)*2+Ostrength,i);
%                     ad = ii;
%                 end
%             end
%             %����
%             f = data((ad-1)*2+Ofreq,i);
%             l = data((ad-1)*2+Ostrength,i);
%             data((ad-1)*2+Ofreq,i)      = data(2+Ofreq,i);
%             data((ad-1)*2+Ostrength,i)  = data(2+Ostrength,i);
%             data(2+Ofreq,i) = f;
%             data(2+Ostrength,i) = l;
             
    end %forj=2��limit-1
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

%��̬�滮����
%�ҵ�����candidate��
maxnc = getMaxnCandidates(nCandidate);
maxnc
ceiling2 = maximumPitch;
timeStepCorrection = 0.01/timeStep;%��

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
    %voice�ͷ�voice���ֿ���
	for icand=1:nCandidate(i)
		voiceless = ~(data((icand-1)*2+Ofreq,i) > 0.0 && data((icand-1)*2+Ofreq,i) < ceiling2);
        if(voiceless)
            delta(i,icand) =  unvoicedStrength;
        else
            delta(i,icand) = data((icand-1)*2+Ostrength,i) - octaveCost * log2(maximumPitch / data((icand-1)*2+Ofreq,i));
        end
    end
end


%Ѱ��һ�����ŵ�·��
transitionCost=0.0;
for i=2:numberOfFrames
    for icand2=1:nCandidate(i)          %i֡��ѡ����
        f2 = data((icand2-1)*2+Ofreq,i);    %i֡��icand2��ѡ�ߵ�Ƶ��
        place = 0;
        maximum = -1e30;
        %����ǰһ֡ÿ��candidate����ǰicand2��cost
        for icand1=1:nCandidate(i-1)    %iǰһ֡��ѡ����
            f1 = data((icand1-1)*2+Ofreq,i-1);  %iǰһ֡icand1��ѡ�ߵ�Ƶ��
            previousVoiceless = ~(f1 > 0.0 && f1 < ceiling2);
            currentVoiceless = ~(f2 > 0.0 && f2 < ceiling2);
            if(currentVoiceless)
                if(previousVoiceless)
                    transitionCost = 0.0;   %��ǰ��������֮ǰ������
                else
                    transitionCost = voicedUnvoicedCost;    %��ǰ��������֮ǰ����
                end
            else
                if(previousVoiceless)
                    transitionCost = voicedUnvoicedCost;    %��ǰ������֮ǰ������
                else
                    transitionCost = octaveJumpCost * abs(log2 (f1 / f2));    %��ǰ������֮ǰ����
                end
            end
            
            value = delta(i-1,icand1) - transitionCost + delta(i,icand2);
            if(value > maximum)
                maximum = value;
                place = icand1;
            elseif(value == maximum)
                
            
            end
        end %����ǰһ֡ÿ��candidate����ǰicand2��cost   ����
        delta(i,icand2) = maximum; %����delta��i֡��icand2��ѡ�ߵĻ���Ϊpreallcost
        psi(i,icand2) = place;      %psi��¼ǰһ֡��Сcostλ��
    end
end
%�ҵ�����·��
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
    %��placeλ�ú͵�һ��candidate���н���, 
    f = data((place-1)*2+Ofreq,iframe);
    l = data((place-1)*2+Ostrength,iframe);
    data((place-1)*2+Ofreq,iframe)      = data(Ofreq,iframe);
    data((place-1)*2+Ostrength,iframe)  = data(Ostrength,iframe);
    data(Ofreq,iframe) = f;
    data(Ostrength,iframe) = l;
    place = psi(iframe,place);
end

if debugplot ~= 0
    figure('NumberTitle', 'off', 'Name', '����·��');
    hold on
     title('����·��');
     xlabel('֡��');  %x��
     ylabel('Ƶ��');%y��
    plot(data(Ofreq,:),'b*');
end

for i=1:length(data(Ofreq,:))
    fprintf(datafile, '%f\t%f\n',((ccc-1)*ttt/44100+i*timeStep),data(Ofreq,i));
end

end

