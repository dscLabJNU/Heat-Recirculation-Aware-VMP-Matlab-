# The SABA algorithm was proposed in the paper "Towards Heat-Recirculation-Aware Virtual Machine Placement in Data Centers"（Hao Feng, Y. Deng, Y. Zhou and G. Min, "Towards Heat-Recirculation-Aware Virtual Machine Placement in Data Centers," in IEEE Transactions on Network and Service Management, online, doi: 10.1109/TNSM.2021.3120295. ）.
## 1. Files
1. D.mat
2. Trace2.xlsx
3. XINTSA.m
4. SABA.m
5. PathLength.m
6. NewAnswer.m
7. Metropolis.m
8.  OutputPath.m

## 2. Introduction
1. D.mat
The heat-recirculation matrix
2. Trace2.xlsx
The numbers of first column of this table represent different time periods, and the numbers of second column represent the number of VMs needed.
3. XINTSA.m
The SA algorithm. This procedure can calculate the power of data center by using SA, and output the optimization process.
4. SABA.m
The SABA algorithm. Using Trace2.xlsx file to calculate, and output the comparisons results of the power of data center with SA algorithm.
5. PathLength.m
Calculate the Tsup.
6. NewAnswer.m
Generating a new solution.
7. Metropolis.m
The Metropolis Criterion.
8.  OutputPath.m
Output the VMP results.
## 3. Execution
Run the SABA.m.
