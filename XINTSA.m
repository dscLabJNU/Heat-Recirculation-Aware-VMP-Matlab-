function totalpower=XINTSA(M)
% clc;
% clear;
% close all;
%%
tic
T0=1200;   % ��ʼ�¶�
Tend= 1.8739277e-11;  % ��ֹ�¶�  1.873927703884808e-11
L=100;    % ���¶��µĵ���������������
q=0.9;    %��������

%% ��������
load D;



%%
N=size(D,1);    %�ڵ�ĸ���
%% ��ʼ��
  %�������һ����ʼ�������
%  M=800%M��������
 S1 =  randi([1,10],1,N); %����1*N,1-10�������������N=50������Ϊ��50���ڵ㡣������D��ʾ50���ڵ��໥֮�������ѧ���á�
 S1 = round(M*S1/(sum(S1)));%round�����������������900�������������������У���ʱS1��ʾ���������
 
  %ȥ������������������
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
 
 %���õ����ڵ����������ڽڵ��������
 
 [a,i]=max(S1);
 [b,j]=min(S1);
 while (a>40)
     c=a-40;
     S1(i)=40;
     S1(j)=S1(j)+c;
     [a,i]=max(S1);
     [b,j]=min(S1);
 end
 

%% ��������Ľڵ���������¶Ⱥ��������ĸ��ڵ��������
disp('��ʼ��Ⱥ�е�һ�����ֵ:')
OutputPath(S1);%�������������S1��� 
yy=sum(S1) %yy=����������Ҳ����m
Rlength=PathLength(D,S1); %Rlength������������ѧ����D�Լ�����������S1������ڵ����������¶�
disp(['��ʼ����ϵͳ�ṩ�����¶ȣ�',num2str(Rlength)]);



%% ��������Ĵ���Time
syms x
Time=ceil(double(solve(1000*(0.9)^x==Tend,x)));
% Time=ceil(double(solve(['1000*(0.9)^x=',num2str(Tend)])));  %��û��ָ��������ʱ��matlabĬ�������ǹ���x��һԪ���η��̵Ľ�
count=0;                   %��������
Obj=zeros(Time,1);         %Ŀ��ֵ�����ʼ��
track=zeros(Time,N);       %ÿ��������·�߾����ʼ��
%% ����
while T0>Tend
    count=count+1;     %���µ�������
    temp=zeros(L,N+1);
    for k=1:L
        %% �����½�
        S2=NewAnswer(S1);
        %% Metropolis�����ж��Ƿ�����½�
        [S1,R]=Metropolis(S1,S2,D,T0);  %Metropolis �����㷨
        temp(k,:)=[S1 R];          %��¼��һ·�ߵļ���·��
    end
    %% ��¼ÿ�ε������̵�����·��
%     [d0,index]=min(temp(:,end)); %??�ҳ���ǰ�¶�������·�� [Y,U]=min(A)������������Y��U��Y������¼A��ÿ�е���Сֵ��U������¼ÿ����Сֵ���кš�
%     if count==1 || d0<Obj(count-1)
%         Obj(count)=d0;           %�����ǰ�¶�������·��С����һ·�����¼��ǰ·��
%     else
%         Obj(count)=Obj(count-1);%�����ǰ�¶�������·�̴�����һ·�����¼��һ·��
%     end
    [d0,index]=max(temp(:,end));
    if count==1 || d0>Obj(count-1)
        Obj(count)=d0;           %�����ǰ�¶�������·��С����һ·�����¼��ǰ·��
    else
        Obj(count)=Obj(count-1);%�����ǰ�¶�������·�̴�����һ·�����¼��һ·��
    end
    track(count,:)=temp(index,1:end-1);  %��¼��ǰ�¶ȵ�����·��
    T0=q*T0;     %����
%     fprintf(1,'%d\n',count)  %�����ǰ��������
end
%% �Ż����̵���ͼ
% figure
% plot(1:count,Obj)
% xlabel('��������')
% ylabel('����')
% title('�Ż�����')


%% ������Ž�Ľڵ���������¶Ⱥ��������ĸ��ڵ��������
disp('�Ż����������ĸ��ڵ���������')
S=track(end,:);
%qq=sum(S)
O=S;
p=OutputPath(S);
B = (S~=0);
B=sum(B(:))
disp(['�Ż�������ϵͳ�ṩ�����¶ȣ�',num2str(PathLength(D,S))]);

disp('-------------------------------------------------------------')
%% �����Ż���ĵ�λʱ��
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
%�������ܺģ�����һ̨������Ϊ2100w
coolingpower=serverpower/(0.0068*PathLength(D,S)*PathLength(D,S)+0.0008*PathLength(D,S)+0.458)
%�����ܺ�
totalpower=serverpower+coolingpower
toc