# Heat-Recirculation-Aware-Virtual-Machine-Placement
# Matlab虚拟机放置SABA算法已发表于Towards Heat-Recirculation-Aware Virtual Machine Placement in Data Centers论文（Hao Feng, Y. Deng, Y. Zhou and G. Min, "Towards Heat-Recirculation-Aware Virtual Machine Placement in Data Centers," in IEEE Transactions on Network and Service Management, online, doi: 10.1109/TNSM.2021.3120295. ）。
## 1. 文件列表
1. D.mat
2. Trace2.xlsx
3. XINTSA.m
4. SABA.m
5. PathLength.m
6. NewAnswer.m
7. Metropolis.m
8.  OutputPath.m

## 2. 文件介绍
1. D.mat
热循环干扰系数矩阵
2. Trace2.xlsx
有两列，第一列是时刻，第二列是任务数量。
3. XINTSA.m
实现SA算法，计算数据中心能耗，并打印优化过程。
4. SABA.m
从Trace2.xlsx文件中分别读取所有时刻的所有任务数量，计算各个时刻的能耗，打印出SA算法和SABA算法的总能耗对比图。
5. PathLength.m
计算出当前任务分配下CRAC需要提供的温度。
6. NewAnswer.m
产生新解，与原来的解对比。
7. Metropolis.m
新解取得较好的效果，接受新解，新解未取得较好的效果，以一定概率接受新解，避免产生局部最优解。
8.  OutputPath.m
输出任务分配情况。
## 3. 运行
运行SABA文件即可运行程序。
