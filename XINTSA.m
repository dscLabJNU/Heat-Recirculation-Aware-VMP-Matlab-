function totalpower=XINTSA(M)
% clc;
% clear;
% close all;
%%
tic
T0=1200;   % 初始温度
Tend= 1.8739277e-11;  % 终止温度  1.873927703884808e-11
L=100;    % 各温度下的迭代次数（链长）
q=0.9;    %降温速率

%% 加载数据
load D;



%%
N=size(D,1);    %节点的个数
%% 初始解
  %随机产生一个初始任务分配
%  M=800%M任务总量
 S1 =  randi([1,10],1,N); %生成1*N,1-10的随机整数矩阵（N=50），因为有50个节点。而矩阵D表示50个节点相互之间的热力学作用。
 S1 = round(M*S1/(sum(S1)));%round四舍五入成整数，将900个任务随机分配进矩阵中，此时S1表示分配情况。
 
  %去掉四舍五入多余的任务
    flag=1;
    while flag==1
        index=randperm(length(S1(:)),1);
        if (S1(index)+M-sum(S1))>=0
            S1(index)= S1(index)+M-sum(S1);
            flag=0;
        else
            S1(index)=0;
        end
    end
 
 %不让单个节点任务数大于节点最大负载数
 
 [a,i]=max(S1);
 [b,j]=min(S1);
 while (a>40)
     c=a-40;
     S1(i)=40;
     S1(j)=S1(j)+c;
     [a,i]=max(S1);
     [b,j]=min(S1);
 end
 

%% 输出随机解的节点入口气流温度和数据中心各节点的任务数
disp('初始种群中的一个随机值:')
OutputPath(S1);%将随机分配任务S1输出 
yy=sum(S1) %yy=任务总量，也就是m
Rlength=PathLength(D,S1); %Rlength函数依据热力学矩阵D以及任务分配情况S1计算出节点的入口气流温度
disp(['初始制冷系统提供气流温度：',num2str(Rlength)]);



%% 计算迭代的次数Time
syms x
Time=ceil(double(solve(1000*(0.9)^x==Tend,x)));
% Time=ceil(double(solve(['1000*(0.9)^x=',num2str(Tend)])));  %当没有指定变量的时候matlab默认求解的是关于x的一元二次方程的解
count=0;                   %迭代计数
Obj=zeros(Time,1);         %目标值矩阵初始化
track=zeros(Time,N);       %每代的最优路线矩阵初始化
%% 迭代
while T0>Tend
    count=count+1;     %更新迭代次数
    temp=zeros(L,N+1);
    for k=1:L
        %% 产生新解
        S2=NewAnswer(S1);
        %% Metropolis法则判断是否接受新解
        [S1,R]=Metropolis(S1,S2,D,T0);  %Metropolis 抽样算法
        temp(k,:)=[S1 R];          %记录下一路线的及其路程
    end
    %% 记录每次迭代过程的最优路线
%     [d0,index]=min(temp(:,end)); %??找出当前温度下最优路线 [Y,U]=min(A)：返回行向量Y和U，Y向量记录A的每列的最小值，U向量记录每列最小值的行号。
%     if count==1 || d0<Obj(count-1)
%         Obj(count)=d0;           %如果当前温度下最优路程小于上一路程则记录当前路程
%     else
%         Obj(count)=Obj(count-1);%如果当前温度下最优路程大于上一路程则记录上一路程
%     end
    [d0,index]=max(temp(:,end));
    if count==1 || d0>Obj(count-1)
        Obj(count)=d0;           %如果当前温度下最优路程小于上一路程则记录当前路程
    else
        Obj(count)=Obj(count-1);%如果当前温度下最优路程大于上一路程则记录上一路程
    end
    track(count,:)=temp(index,1:end-1);  %记录当前温度的最优路线
    T0=q*T0;     %降温
%     fprintf(1,'%d\n',count)  %输出当前迭代次数
end
%% 优化过程迭代图
% figure
% plot(1:count,Obj)
% xlabel('迭代次数')
% ylabel('距离')
% title('优化过程')


%% 输出最优解的节点入口气流温度和数据中心各节点的任务数
disp('优化后数据中心各节点分配的任务：')
S=track(end,:);
%qq=sum(S)
O=S;
p=OutputPath(S);
B = (S~=0);
B=sum(B(:))
disp(['优化后制冷系统提供气流温度：',num2str(PathLength(D,S))]);

disp('-------------------------------------------------------------')
%% 计算优化后的单位时间
% time1=ceil(0.05*S);
% % for rr=1:1
% %     time2(ZZ(rr))=0; 
% % end
% time1;
% cs1=max(time1)
% T1=zeros(1,cs1);
% for ii1=1:cs1
% T1(ii1)=numel(find(time1==ii1));
% end
% T1
serverpower=B*2100
%服务器能耗，启动一台服务器为2100w
coolingpower=serverpower/(0.0068*PathLength(D,S)*PathLength(D,S)+0.0008*PathLength(D,S)+0.458)
%制冷能耗
totalpower=serverpower+coolingpower
toc