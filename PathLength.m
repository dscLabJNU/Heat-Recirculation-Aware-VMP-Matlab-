%% ����Rlength tsup��Сֵ
% ���룺
% D     ��ѭ������
% Chrom ��ʼ����Ⱥ
function len=PathLength(D,Chrom)
%% ����յ��ṩ�¶�

[row,col]=size(D);
NIND=size(Chrom,1); %����䷵�ص�ʱ����Chrom������
len=zeros(NIND,1);
% for i=1:NIND
%     p=[Chrom(i,:) Chrom(i,1)];
%     i1=p(1:end-1);
%     i2=p(2:end);
%     len(i,1)=sum(D((i1-1)*col+i2));
% end
% %% ʧЧ����X
% x=ones(1,50);
% r1=randsrc(1,1,[1:50]);
% x(r1)=0
% x
% X=diag(x)
%% �����С�յ��ṩ�¶�ֵ
a=50;

Db=zeros(NIND,row);
for i=1:NIND
    Db(i,:)=(2020*ones(1,row))*D';
end

Tin=26*ones(NIND,row);  % �ڵ�����¶�
 

Tsup=Tin-Db-a*Chrom*D';  % Tin-(Tsup+Db)
% T=a*Chrom*D';
l=(min(Tsup,[],2));
for i=1:NIND
    len(i,1)=l(i);
end
