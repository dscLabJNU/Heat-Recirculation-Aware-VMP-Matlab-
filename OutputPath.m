function p=OutputPath(R)
%% ���·������
%���룺R ·��
R=[R,R(1)]
N=length(R);
p=['��>',num2str(R(1))];%����ֵת��Ϊ�ַ����������
for i=2:N-1
    if mod(i,5)~=0
    p=[p,'��>',num2str(R(i))];
    else
            p=[p,'��>',num2str(R(i))];
            p=[p,char(13)]; %����
    end
end
disp(p);