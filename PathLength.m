%% 计算Rlength tsup最小值
% 输入：
% D     热循环矩阵
% Chrom 初始化种群
function len=PathLength(D,Chrom)
%% 计算空调提供温度

[row,col]=size(D);
NIND=size(Chrom,1); %该语句返回的时矩阵Chrom的行数
len=zeros(NIND,1);
% for i=1:NIND
%     p=[Chrom(i,:) Chrom(i,1)];
%     i1=p(1:end-1);
%     i2=p(2:end);
%     len(i,1)=sum(D((i1-1)*col+i2));
% end
% %% 失效矩阵X
% x=ones(1,50);
% r1=randsrc(1,1,[1:50]);
% x(r1)=0
% x
% X=diag(x)
%% 求出最小空调提供温度值
a=50;

Db=zeros(NIND,row);
for i=1:NIND
    Db(i,:)=(2020*ones(1,row))*D';
end

Tin=26*ones(NIND,row);  % 节点入口温度
 

Tsup=Tin-Db-a*Chrom*D';  % Tin-(Tsup+Db)
% T=a*Chrom*D';
l=(min(Tsup,[],2));
for i=1:NIND
    len(i,1)=l(i);
end
