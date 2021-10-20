clc;
clear;
close all;
%%
tic


%% ��������
load D;

%% ��ȡxls�ļ�
A=xlsread('Trace2.xlsx');
time1=A(:,1);
tasks=A(:,2);
result1=zeros(1,length(tasks(:)));
result2=zeros(1,length(tasks(:)));


%%
N=size(D,1);    %�ڵ�ĸ���

for n=1:length(tasks(:))
    disp(['n��ֵ---------����������',num2str(n)]);
    
    T0=1000;   % ��ʼ�¶�
    Tend= 1.8739277e-11;  % ��ֹ�¶�  1.873927703884808e-11
    L=100;    % ���¶��µĵ���������������
    q=0.9;    %��������
    
    %% ��ʼ��
    %�������һ����ʼ�������
    M=tasks(n); %�����������ⶥ2000
    if M<=1300
        SERVER=ceil(M/0.65/40);%�������õķ�����������0.65��ʾ������ƽ��ʹ����
        disp('���õķ�����������')
        disp(SERVER)
        S1 =  randi([1,10],1,SERVER);
        S0 = zeros(1,N-SERVER);
        S1 = [S1,S0];
        S1 = round(M*S1/(sum(S1)));%round�����������������900�������������������У���ʱS1��ʾ���������
    else
        S1 =  randi([1,10],1,N); %����1*N,1-10�������������N=50������Ϊ��50���ڵ㡣������D��ʾ50���ڵ��໥֮�������ѧ���á�
        S1 = round(M*S1/(sum(S1)));%round�����������������900�������������������У���ʱS1��ʾ���������
        SERVER=(S1~=0);
        SERVER=sum(SERVER(:));
        disp('���õķ�����������')
        disp(SERVER)
    end
   %ȥ������������������
    flag=1;
    while flag==1
        index=randperm(length(S1(S1~=0)),1);
        if (S1(index)+M-sum(S1))>0
            S1(index)= S1(index)+(M-sum(S1));
            flag=0;
        end
    end
    
    %���õ����ڵ����������ڽڵ��������
    
    [a,i]=max(S1);
    [b,j]=min(S1(S1~=0));
    while (a>40)
        c=a-40;
        S1(i)=40;
        S1(j)=S1(j)+c;
        [a,i]=max(S1);
        [b,j]=min(S1(S1~=0));
    end
    
    %% ��������Ľڵ���������¶Ⱥ��������ĸ��ڵ��������
    disp('��ʼ��Ⱥ�е�һ�����ֵ:')
    OutputPath(S1);%�������������S1���
    yy=sum(S1) %yy=����������Ҳ����m
    Rlength=PathLength(D,S1); %Rlength������������ѧ����D�Լ�����������S1������ڵ����������¶�
    disp(['��ʼ����ϵͳ�ṩ�����¶ȣ�',num2str(Rlength)]);
    
%     initserverpower=SERVER*2100%�������ܺģ�����һ̨������Ϊ2100w
%     initcoolingpower=initserverpower/(0.0068*PathLength(D,S1)*PathLength(D,S1)+0.0008*PathLength(D,S1)+0.458)
%     %�����ܺ�
%     inittotalpower=initserverpower+initcoolingpower
%     % ��¼��ʼ�ܺ�
%     result1(n)=initcoolingpower;
    
    %% ���������г����Ż�
    S2=S1;
    for k=5:5:50
        [x,i]=max(S2);
        y=S1(k);
        S1(k)=S1(i);
        S1(i)=y;
        S2(i)=0;  %��S2�����ֵ��0����Ѱ����һ�����ֵ��
        S2(k)=0;
    end
    disp('�����ֵ���г����Ż�:')
    OutputPath(S1);
    
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
            % if(S2~=S1)
            %% Metropolis�����ж��Ƿ�����½�
            [S1,R]=Metropolis(S1,S2,D,T0);  %Metropolis �����㷨
            temp(k,:)=[S1 R];          %��¼��һ·�ߵļ���·��
            %else
            %   break
            %end
        end
        %% ��¼ÿ�ε������̵�����·��
        [d0,index]=min(temp(:,end)); %�ҳ���ǰ�¶�������·�� [Y,U]=min(A)������������Y��U��Y������¼A��ÿ�е���Сֵ��U������¼ÿ����Сֵ���кš�
        if count==1 || d0<Obj(count-1)
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
    disp('�Ż��������ķ�����������')
    disp(B);
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
    %��¼�Ż�����ܺ�
    result2(n)=coolingpower;
    
    %�������з�����
    result1(n)=XINTSA(M);
    toc
end

hold on
x=time1';
y1=result1(x);
y2=result2(x);
values1=spcrv([[x(1) x x(end)];[y1(1) y1 y1(end)]],3);
values2=spcrv([[x(1) x x(end)];[y2(1) y2 y2(end)]],3);
plot(values1(1,:),values1(2,:),'sc-');
plot(values2(1,:),values2(2,:),'dr-');
legend('SA','SABA');
xlabel('Time');
ylabel('totalpower');
set(gca,'XTick',0:10:time1(end));
hold off