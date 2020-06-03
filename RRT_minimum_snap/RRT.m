%***************************************
%Author: Chaoqun Wang
%Date: 2019-10-15
%***************************************
%% ���̳�ʼ��
clear all; close all;
x_I=1; y_I=1;           % ���ó�ʼ��
x_G=700; y_G=700;       % ����Ŀ���
Thr=50;                 %����Ŀ�����ֵ
Delta= 30;              % ������չ����
%% ������ʼ��
T.v(1).x = x_I;         % T������Ҫ��������v�ǽڵ㣬�����Ȱ���ʼ����뵽T������
T.v(1).y = y_I; 
T.v(1).xPrev = x_I;     % ��ʼ�ڵ�ĸ��ڵ���Ȼ���䱾��
T.v(1).yPrev = y_I;
T.v(1).dist=0;          %�Ӹ��ڵ㵽�ýڵ�ľ��룬�����ȡŷ�Ͼ���
T.v(1).indPrev = 0;     %
%% ��ʼ������������ҵ����
figure(1);
ImpRgb=imread('newmap.png');
Imp=rgb2gray(ImpRgb);
imshow(Imp)
xL=size(Imp,1);%��ͼx�᳤��
yL=size(Imp,2);%��ͼy�᳤��
hold on
plot(x_I, y_I, 'ro', 'MarkerSize',10, 'MarkerFaceColor','r');
plot(x_G, y_G, 'go', 'MarkerSize',10, 'MarkerFaceColor','g');% ��������Ŀ���
count=1;
for iter = 1:3000
    x_rand=[];
    %Step 1: �ڵ�ͼ���������һ����x_rand
    %��ʾ���ã�x_rand(1),x_rand(2)����ʾ�����в����������
    x_rand(1) = randi(xL);
    x_rand(2) = randi(yL);
    
    x_near=[];
    %Step 2: ���������������ҵ�����ڽ���x_near 
    %��ʾ��x_near�Ѿ�����T��
    distance = sqrt((T.v(1).x-x_rand(1))^2+(T.v(1).y-x_rand(2))^2);
    x_near(1) = T.v(1).x;
    x_near(2) = T.v(1).y;
    temp = 1;
    for i=1:count
        dis_temp = sqrt((T.v(i).x-x_rand(1))^2+(T.v(i).y-x_rand(2))^2);
        if(distance>dis_temp)
            distance = dis_temp;
            x_near(1) = T.v(i).x;
            x_near(2) = T.v(i).y;
            temp = i;
        end
        
    end
    
    x_new=[];
    %Step 3: ��չ�õ�x_new�ڵ�
    %��ʾ��ע��ʹ����չ����Delta
    x_new(1) = x_near(1)+(x_rand(1) - x_near(1)) * Delta / distance;
    x_new(2) = x_near(2)+(x_rand(2) - x_near(2)) * Delta / distance;
    
    %���ڵ��Ƿ���collision-free
    if ~collisionChecking(x_near,x_new,Imp) 
        continue;
    end
    count=count+1;
    
    %Step 4: ��x_new������T 
    %��ʾ���½ڵ�x_new�ĸ��ڵ���x_near
    
    T.v(count).x = x_new(1);
    T.v(count).y = x_new(2);
    T.v(count).xPrev = x_near(1);
    T.v(count).yPrev = x_near(2);
    T.v(count).dist = Delta;
    T.v(count).indPrev = temp;
    
    %Step 5:����Ƿ񵽴�Ŀ��㸽�� 
    %��ʾ��ע��ʹ��Ŀ�����ֵThr������ǰ�ڵ���յ��ŷʽ����С��Thr����������ǰforѭ��
    if(sqrt((x_new(1)-x_G)^2+(x_new(2)-y_G)^2) < Thr)
        break;
    end
    
   %Step 6:��x_near��x_new֮���·��������
   %��ʾ 1��ʹ��plot���ƣ���ΪҪ�����ͬһ��ͼ�ϻ����߶Σ�����ÿ��ʹ��plot����Ҫ����hold on����
   %��ʾ 2�����ж��յ���������forѭ��ǰ���ǵð�x_near��x_new֮���·��������
   plot([x_new(1); x_near(1);], [x_new(2); x_near(2)], 'r', 'Linewidth', 3);
   hold on
   pause(0.1); %��ͣ0.1s��ʹ��RRT��չ�������׹۲�
end
%% ·���Ѿ��ҵ��������ѯ
if iter < 2000
    path = [];
    path(1,1) = x_G; path(1,2) = y_G;
    path(2,1) = T.v(end).x; path(2,2) = T.v(end).y;
    
    pathIndex = T.v(end).indPrev; % �յ����·��
    fprintf('%d\n',T.v(end).indPrev);
    fprintf('%d\n',pathIndex);
    j=0;
    while 1
        path((j+3),1) = T.v(pathIndex).x;
        path((j+3),2) = T.v(pathIndex).y;
        pathIndex = T.v(pathIndex).indPrev;
        if pathIndex == 1
            break
        end
        j=j+1;
    end  % ���յ���ݵ����
    path((end+1),1) = x_I; path((end),2) = y_I; % ������·��
    for j = 2:length(path)
        plot([path((j),1); path((j-1),1);], [path((j),2); path((j-1),2)], 'b', 'Linewidth', 3);
    end
    hw1_1(path);
else
    disp('Error, no path found!');
end

