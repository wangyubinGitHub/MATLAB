function  testwav3()
%UNTITLED4 此处显示有关此函数的摘要
%   此处显示详细说明

filename = '../sourcefile/source1.wav'

[y, fs] = audioread(filename);
r_y=y;

word_time=[1,633;%你
    668,1075;%的
    1093,1539;%城
    1540,2651];%市
word_count=size(word_time,1);

%bi=1.059463%半音
bi=1.33484%五个半音
%目标频率
pitch_desire=[110,0,18749-18263;
    110,18749-18263,19196-18263;
    110,19196-18263,19737-18263;
    139,19737-18263,20742-18263];

%求取基频
steps= 1536*2;
section=1024;
sample_rate=fs;
pitch=zeros(4,floor((length(y)-steps)/section));%1，时间尺度，2，基频值，3，组id，4，目标频率
length(y)%samples
for m=1:floor((length(y)-steps)/section) %包含多少个section采样点，并保证最后一个的数据量，有steps长度，steps是后面做自相关的长度
    pitch(1,m) = m*section*1000/sample_rate;%pitch(1,:)保留了time值(ms)，方便作图
    t=m*section/sample_rate;
%    if(t<4)%这里求取时间，对时长进行控制
        x=y((m-1)*section+1:(m-1)*section+steps);%从y的第一个元素计算起，取steps长度数据，做自相关
        a=myxcorr(x);
        [y_peaks,t_peaks]=findpeaks(a);%查找峰值，查找最大值
        [weizhi,flag]=max(y_peaks);
        flag1=t_peaks(flag);%峰值的采样点数
        if((sample_rate/flag1) < 1000 && sample_rate/flag1 >50 )%这里对基频做一个限制，要大于50，小于1000
            pitch(2,m+1)=sample_rate/flag1;
        end
%    end
end
% olo 对基频的连续性，准确性做个预估，将连续的跨度不大的分为一组，分配组id，其余的置0
group=zeros(length(pitch),2);%
group_count=1;
status=0;%0，组外，1，组中
for i=1:size(pitch,2)
    if(pitch(2,i) > 0 &&i < size(pitch,2))
        if(status==0)
           group(group_count,1)=i;
           status=1;%组中
        elseif(pitch(2,i)/pitch(2,i-1) > 1.5 || pitch(2,i)/pitch(2,i-1) < 0.5)%跨度大，认为不连续
            group(group_count,2)=i-1;
            group_count = group_count + 1;
            group(group_count,1)=i;
        end
    else
        if(status == 1)
            if( i==length(pitch) )
                group(group_count,2)=i;
            else
                group(group_count,2)=i-1;
            end
           group_count = group_count + 1;
           status = 0;
        end
    end
end
plot(pitch(1,:),pitch(2,:),'go');%绘制最初基频点图
hold on
%遍历每个组，删除连续数量不足5的野组,再删除STDs过大的组
se=0;%每个组用红蓝间隔开来好分辨
for i=1:group_count-1
    if(group(i,1) > 0 && group(i,2)>0 && group(i,2)-group(i,1) < 3)
       pitch(2,group(i,1):group(i,2)) = 0; 
    else
        arrA=pitch(2,group(i,1):group(i,2)-1);
        arrB=pitch(2,group(i,1)+1:group(i,2));
        arrV=arrB-arrA;%前一个值减后一个值的差
        stds=var(arrV);
        if(stds > 1000)
            pitch(2,group(i,1):group(i,2)) = 0;
        else
            i
            if(se==0)
                plot(pitch(1,group(i,1):group(i,2)),pitch(2,group(i,1):group(i,2)),'r*');%绘制基频点图
                se =1;
            else
                se=0;
                plot(pitch(1,group(i,1):group(i,2)),pitch(2,group(i,1):group(i,2)),'b*');%绘制基频点图
            end
            pitch(3,group(i,1):group(i,2))=i;%记录下group id

        end

    end

end
% olo 
%mark 根据基频做mark
%begin_samples和end_samples分别是一个字的开始和结束位置，坐标是位于y中的采样点位置
begin_samples=0;
end_samples=0;

%因为每个字的基频可能有较大的跨度被分成了两个段，所以每个字也需要在进行细分成更小的组
%每个组要包含 开始、结束位置，基频（可以通过pitch求），desire频率(可以从字里拆解出来,通过下面这个 步骤就?〈?pitch(4,)里面取)
group_count=0;
mark_group=[];

for i=1:word_count
    begin_samples=floor(word_time(i,1)*sample_rate/1000);
    end_samples=floor(word_time(i,2)*sample_rate/1000);
    begin_m=floor(begin_samples/section);
    if(begin_m == 0)
        begin_m=1;
    end
    end_m=floor(end_samples/section);
    pitch(4, floor(begin_m):floor(end_m))=pitch_desire(i,1)*bi;
    n(1:end_m-begin_m+1)=0;
    for m=begin_m:end_m
       if(pitch(2,m)>0)
           n(m-begin_m+1)=pitch(3,m);%记录组number
       end
    end
    o=unique(n);
    b=find(o~=0);
    group_count=group_count+length(b);
    for m=1:length(b)
       if(group(o(b(m)),1)<=begin_m)
          t_begain=begin_m;
       else
          t_begain=group(o(b(m)),1);
       end
       if(group(o(b(m)),2)>=end_m)
           t_end=end_m;
       else
           t_end=group(o(b(m)),2);
       end
       if(t_end-t_begain>2)
           g=[t_begain,t_end,pitch_desire(i,1)*bi];
           mark_group=[mark_group;g];
           begin_m=t_end+1;
       end
    end
end
mark_group

for cl=1:size(mark_group,1)
    begin_samples=mark_group(cl,1)*section;
    end_samples=floor(mark_group(cl,2)*section);
    begin_samples
    end_samples
    if(begin_samples==0)
       begin_samples=1; 
    end
    if(end_samples >= length(y)-steps)
        end_samples=length(y)-steps;
    end
    %length(y)
    mark_addr=zeros(1,length(y));%第一步mark，根据基频周期mark每个基频周期的峰值点
    mark_z_addr=zeros(1,length(y));%第二步mark，记录目标需要overlap的目标点；如果目标点比第一步mark的多，即提高音调
    %首先在这里对基频的准确性，连续性做个筛选
    for i=begin_samples:end_samples
        
    end
    
    %a. mark找到第一个峰值点
    [mark_y,mark_r] = findpeaks(y(begin_samples:end_samples));%峰值
    [smark_max,s_index]=sort(mark_y);%峰值升序排列
    for oo=0:length(mark_y)-1%找能检测出基频的最大峰值
       cur_mark=begin_samples+mark_r(s_index(length(mark_y)-oo));
       if(pitch(2,floor( cur_mark /section)+1) >50 && pitch(2,floor( cur_mark /section)+1)<1000)%这边再次对基频做了判断
            mark_addr(1,cur_mark)= smark_max(length(mark_y)-oo);
            curT=sample_rate/pitch(2,floor( cur_mark /section)+1);
            break;
       end
    end
    tmp_cur_mark=cur_mark;
    tmp_curT=curT;
    count=0;
    b_s=cur_mark;
    e_s=cur_mark;
    %b. mark根据a找到的峰值点，继续完成mark
    while cur_mark<end_samples
       [cur_y,cur_r] = findpeaks( y(floor(cur_mark+0.7*curT):floor(cur_mark+1.3*curT)));
       [cur_max,cur_i] = max(cur_y); 
       mark_addr(1,floor(cur_mark+0.7*curT)+cur_r(cur_i))=cur_max;
       cur_mark=floor(cur_mark+0.7*curT)+cur_r(cur_i);
       count=count+1;
       e_s=cur_mark;
       if(pitch(2,floor( cur_mark /1024))>0)
       curT = sample_rate/pitch(2,floor( cur_mark /section)+1);
       end
    end
    
    cur_mark = tmp_cur_mark;
    curT=tmp_curT;
    while cur_mark>begin_samples
       [cur_y,cur_r] = findpeaks( y(floor(cur_mark-1.3*curT):floor(cur_mark-0.7*curT)));
       [cur_max,cur_i] = max(cur_y); 
       if(cur_max > 0.01)
       mark_addr(1,floor(cur_mark-1.3*curT)+cur_r(cur_i))=cur_max;
       cur_mark=floor(cur_mark-1.3*curT)+cur_r(cur_i);
       b_s = cur_mark;
       count=count+1;
       else
           break;
       end
       if(pitch(2,floor( cur_mark /section))>0)
       curT = sample_rate/pitch(2,floor( cur_mark /section)+1);
       end
    end
    
%第一步mark完成    

b_s
e_s
%r_y(b_s:e_s)=0;

%
curF = 1/((e_s-b_s)/count/sample_rate);
curF
desir_count=0;
%desir_count=floor((e_s-b_s)/(sample_rate/pitch_desire(cl,1)));
desir_count=floor((e_s-b_s)/(sample_rate/mark_group(cl,3)));
desir_samples=floor((e_s-b_s)/desir_count);

addr2=zeros(1,count);%保存了mark后的一个索引，这个段第n个mark点对应音频的samples数
tmp_count=1;
for i=1:length(y)
    if(mark_addr(1,i) > 0.01)
        addr2(1,tmp_count) = i;
        tmp_count = tmp_count+1;
    end
end


%z mark
z_begin=b_s;
z_end = b_s+(desir_count-1)*desir_samples;

%加窗清数据
r_y(b_s+2049:e_s-2049)=0;
haha=windows(r_y(b_s-2048:b_s+2048));
r_y(b_s:b_s+2048)=haha(2049:4097);
haha=windows(r_y(e_s-2048:e_s+2048));
r_y(e_s-2048:e_s)=haha(1:2049);


ttt =0;
for i=0:desir_count-1
    %记录起来
    mark_z_addr(1,b_s+i*desir_samples)=0.1;
    %conv
    %1找临近的mark点
    lalal=0;
    left=0;
    right=0;
    while 1
        ttt = b_s+i*desir_samples + lalal;
        if(mark_addr(1,ttt) > 0.01)
           %有临近的 
           break;
        end
        ttt = b_s+i*desir_samples - lalal;
        if(mark_addr(1,ttt) > 0.01)
           %有临近的 
           break;
        end
        lalal= lalal + 1;
    end
    %2加窗
    for x=2:count-1
        if(addr2(1,x)==ttt)
            left = addr2(1,x-1);
            right=addr2(1,x+1);
            mm=y(left:right);
            nn=windows(mm);
            r_y(b_s+i*desir_samples-floor((right-left)/2):b_s+i*desir_samples-floor((right-left)/2)+right-left)=r_y(b_s+i*desir_samples-floor((right-left)/2):b_s+i*desir_samples-floor((right-left)/2)+right-left)+nn;
        end
    end

end

filename = 'handel.wav';
audiowrite(filename,r_y,sample_rate);

figure;
plot(y,'-b');
hold on
for i=1:length(y)
    if(mark_addr(1,i) > 0.01)
        plot(i,mark_addr(1,i),'*r');
    end
    if(mark_z_addr(1,i) > 0.01)
       plot(i,mark_z_addr(1,i),'*b'); 
    end
end

end
end

