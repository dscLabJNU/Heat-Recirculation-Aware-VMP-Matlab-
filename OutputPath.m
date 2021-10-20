function p=OutputPath(R)
%% 输出路径函数
%输入：R 路径
R=[R,R(1)]
N=length(R);
p=['—>',num2str(R(1))];%把数值转换为字符串进行输出
for i=2:N-1
    if mod(i,5)~=0
    p=[p,'—>',num2str(R(i))];
    else
            p=[p,'—>',num2str(R(i))];
            p=[p,char(13)]; %换行
    end
end
disp(p);