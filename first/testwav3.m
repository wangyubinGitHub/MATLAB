function  testwav3()
%UNTITLED4 �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��

filename = '../sourcefile/source1.wav'

[y, fs] = audioread(filename);
r_y=y;

word_time=[1,633;%��
    668,1075;%��
    1093,1539;%��
    1540,2651];%��
word_count=size(word_time,1);

%bi=1.059463%����
bi=1.33484%�������
%Ŀ��Ƶ��
pitch_desire=[110,0,18749-18263;
    110,18749-18263,19196-18263;
    110,19196-18263,19737-18263;
    139,19737-18263,20742-18263];

%��ȡ��Ƶ
steps= 1536*2;
section=1024;
sample_rate=fs;
pitch=zeros(4,floor((length(y)-steps)/section));%1��ʱ��߶ȣ�2����Ƶֵ��3����id��4��Ŀ��Ƶ��
length(y)%samples
for m=1:floor((length(y)-steps)/section) %�������ٸ�section�����㣬����֤���һ��������������steps���ȣ�steps�Ǻ���������صĳ���
    pitch(1,m) = m*section*1000/sample_rate;%pitch(1,:)������timeֵ(ms)��������ͼ
    t=m*section/sample_rate;
%    if(t<4)%������ȡʱ�䣬��ʱ�����п���
        x=y((m-1)*section+1:(m-1)*section+steps);%��y�ĵ�һ��Ԫ�ؼ�����ȡsteps�������ݣ��������
        a=myxcorr(x);
        [y_peaks,t_peaks]=findpeaks(a);%���ҷ�ֵ���������ֵ
        [weizhi,flag]=max(y_peaks);
        flag1=t_peaks(flag);%��ֵ�Ĳ�������
        if((sample_rate/flag1) < 1000 && sample_rate/flag1 >50 )%����Ի�Ƶ��һ�����ƣ�Ҫ����50��С��1000
            pitch(2,m+1)=sample_rate/flag1;
        end
%    end
end
% olo �Ի�Ƶ�������ԣ�׼ȷ������Ԥ�����������Ŀ�Ȳ���ķ�Ϊһ�飬������id���������0
group=zeros(length(pitch),2);%
group_count=1;
status=0;%0�����⣬1������
for i=1:size(pitch,2)
    if(pitch(2,i) > 0 &&i < size(pitch,2))
        if(status==0)
           group(group_count,1)=i;
           status=1;%����
        elseif(pitch(2,i)/pitch(2,i-1) > 1.5 || pitch(2,i)/pitch(2,i-1) < 0.5)%��ȴ���Ϊ������
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
plot(pitch(1,:),pitch(2,:),'go');%���������Ƶ��ͼ
hold on
%����ÿ���飬ɾ��������������5��Ұ��,��ɾ��STDs�������
se=0;%ÿ�����ú�����������÷ֱ�
for i=1:group_count-1
    if(group(i,1) > 0 && group(i,2)>0 && group(i,2)-group(i,1) < 3)
       pitch(2,group(i,1):group(i,2)) = 0; 
    else
        arrA=pitch(2,group(i,1):group(i,2)-1);
        arrB=pitch(2,group(i,1)+1:group(i,2));
        arrV=arrB-arrA;%ǰһ��ֵ����һ��ֵ�Ĳ�
        stds=var(arrV);
        if(stds > 1000)
            pitch(2,group(i,1):group(i,2)) = 0;
        else
            i
            if(se==0)
                plot(pitch(1,group(i,1):group(i,2)),pitch(2,group(i,1):group(i,2)),'r*');%���ƻ�Ƶ��ͼ
                se =1;
            else
                se=0;
                plot(pitch(1,group(i,1):group(i,2)),pitch(2,group(i,1):group(i,2)),'b*');%���ƻ�Ƶ��ͼ
            end
            pitch(3,group(i,1):group(i,2))=i;%��¼��group id

        end

    end

end
% olo 
%mark ���ݻ�Ƶ��mark
%begin_samples��end_samples�ֱ���һ���ֵĿ�ʼ�ͽ���λ�ã�������λ��y�еĲ�����λ��
begin_samples=0;
end_samples=0;

%��Ϊÿ���ֵĻ�Ƶ�����нϴ�Ŀ�ȱ��ֳ��������Σ�����ÿ����Ҳ��Ҫ�ڽ���ϸ�ֳɸ�С����
%ÿ����Ҫ���� ��ʼ������λ�ã���Ƶ������ͨ��pitch�󣩣�desireƵ��(���Դ����������,ͨ��������� �����?��?pitch(4,)����ȡ)
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
           n(m-begin_m+1)=pitch(3,m);%��¼��number
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
    mark_addr=zeros(1,length(y));%��һ��mark�����ݻ�Ƶ����markÿ����Ƶ���ڵķ�ֵ��
    mark_z_addr=zeros(1,length(y));%�ڶ���mark����¼Ŀ����Ҫoverlap��Ŀ��㣻���Ŀ���ȵ�һ��mark�Ķ࣬���������
    %����������Ի�Ƶ��׼ȷ�ԣ�����������ɸѡ
    for i=begin_samples:end_samples
        
    end
    
    %a. mark�ҵ���һ����ֵ��
    [mark_y,mark_r] = findpeaks(y(begin_samples:end_samples));%��ֵ
    [smark_max,s_index]=sort(mark_y);%��ֵ��������
    for oo=0:length(mark_y)-1%���ܼ�����Ƶ������ֵ
       cur_mark=begin_samples+mark_r(s_index(length(mark_y)-oo));
       if(pitch(2,floor( cur_mark /section)+1) >50 && pitch(2,floor( cur_mark /section)+1)<1000)%����ٴζԻ�Ƶ�����ж�
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
    %b. mark����a�ҵ��ķ�ֵ�㣬�������mark
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
    
%��һ��mark���    

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

addr2=zeros(1,count);%������mark���һ������������ε�n��mark���Ӧ��Ƶ��samples��
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

%�Ӵ�������
r_y(b_s+2049:e_s-2049)=0;
haha=windows(r_y(b_s-2048:b_s+2048));
r_y(b_s:b_s+2048)=haha(2049:4097);
haha=windows(r_y(e_s-2048:e_s+2048));
r_y(e_s-2048:e_s)=haha(1:2049);


ttt =0;
for i=0:desir_count-1
    %��¼����
    mark_z_addr(1,b_s+i*desir_samples)=0.1;
    %conv
    %1���ٽ���mark��
    lalal=0;
    left=0;
    right=0;
    while 1
        ttt = b_s+i*desir_samples + lalal;
        if(mark_addr(1,ttt) > 0.01)
           %���ٽ��� 
           break;
        end
        ttt = b_s+i*desir_samples - lalal;
        if(mark_addr(1,ttt) > 0.01)
           %���ٽ��� 
           break;
        end
        lalal= lalal + 1;
    end
    %2�Ӵ�
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

