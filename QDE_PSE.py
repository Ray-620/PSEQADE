# coding:utf-8
import numpy as np
import time
from numpy import random
import math
import matplotlib.pyplot as plt
# from function_cec2017 import func
# from function_test_v2 import func
from function_cec2013 import func


epsilon = np.finfo(np.float64).tiny

def Dis(NP,GP,value,XG,gbest,gbest_x,fai,gfaimin,xmax,xmin): # 分散操作
    span = xmax - xmin
    span2 = xmax + xmin
    for i in range(NP):
        IS = (value[i] - gbest) / (max(value) - gbest + epsilon)    # 表示第i个个体的适应度状态 与最佳值越像->值越小
        RP = 0.5*(1.0 - GP + IS)    # 再生概率
        for j in range(D):      # 更新种群
            if random.random() < RP or (j == random.randint(1, D)):
                DP = (xmax-gbest_x[j]) / (xmax-xmin)    # 距离上界的百分比
                if random.random() < DP:
                    XG[i][j] = gbest_x[j] + random.random()* GP *(xmax-gbest_x[j])
                    tempsin = math.asin((2*XG[i][j] - span2) / span)    # 根据量子种群反向生成实数种群
                    tempcos = math.acos((2*XG[i][j] - span2) / span)
                    if tempcos > tempsin:
                        fai[i][j] = tempcos
                    else:
                        fai[i][j] = tempsin
                else:
                    XG[i][j] = gbest_x[j] + random.random()* GP *(xmin-gbest_x[j])
                    tempsin = math.asin((2*XG[i][j] - span2) / span)
                    tempcos = math.acos((2*XG[i][j] - span2) / span)
                    if tempcos > tempsin:
                        fai[i][j] = tempcos
                    else:
                        fai[i][j] = tempsin
    fitness = []
    for i in range(NP):
        fitness.append(func(func_num, XG[i], 'value'))

    value_min = min(fitness)
    if value_min < gbest:
        # 更新
        gbest = value_min
        pos_min = fitness.index(value_min)
        gbest_x = XG[pos_min]
        gfaimin = fai[pos_min]

    return XG,fai,gbest,gbest_x,gfaimin

    
def Agg(NP,GP,AP,value,XG,fai,gfaimin,cbest,cbest_x,xmax,xmin,gbest,gbest_x):
    span = xmax - xmin
    span2 = xmax + xmin
    for i in range(NP):
        IS = (max(value)-value[i]) / (max(value) - cbest + epsilon)    
        RP = AP * (GP + IS)    # 再生概率
        for j in range(D):
            R = abs(XG[i][j]-cbest_x[j]) / (xmax-xmin)
            if RP < R or (j == random.randint(1, D)):
                XG[i][j] = cbest_x[j] + random.random()* GP *(XG[i][j]-cbest_x[j])
                tempsin = math.asin((2*XG[i][j] - span2) / span)
                tempcos = math.acos((2*XG[i][j] - span2) / span)
                if tempcos > tempsin:
                    fai[i][j] = tempcos
                else:
                    fai[i][j] = tempsin
            else:
                continue
    
    fitness = []
    for i in range(NP):
        fitness.append(func(func_num, XG[i], 'value'))

    value_min = min(fitness)
    if value_min < gbest:
        # 更新
        gbest = value_min
        pos_min = fitness.index(value_min)
        gbest_x = XG[pos_min]
        gfaimin = fai[pos_min]
    return XG,fai,gbest,gbest_x,gfaimin



def QDE_PSE(func_num,Gm,NP,CR,D,Clock,IT):
    # 计时
    time_start = time.time()
    # Gm 最大迭代次数
    Gmin = np.zeros(Gm)  # 各代最优值
    best_x = np.zeros((Gm, D))  # 各代最优解
    gbest = 0   # 历史最优值
    gbest_x = np.zeros(D)   #历史最优个体
    gfaimin = np.zeros(D)
    value = np.zeros((NP))
    G = 0

    fai = math.pi * np.random.rand(NP, D)
    Pm = np.ones((NP, 1)) * 0.05
    X = np.zeros((2 * NP, D))
    trace = []
    for i in range(0, 2 * NP, 2):
        for j in range(D):
            X[i, j] = math.cos(fai[int(0.5 * (i + 1)), j])
            X[i + 1, j] = math.sin(fai[int(0.5 * (i + 1)), j])

    ## 量子空间变换
    # 读取函数边界
    XBounds = func(func_num, [], 'initial')
    xmin = XBounds[0]
    xmax = XBounds[1]
    # print(xmin,xmax)

    span = xmax - xmin
    span2 = xmax + xmin
    X = (span * X + span2) / 2

    ## 量子测量
    XG = np.zeros((NP, D))
    for i in range(NP):
        for j in range(D):
            pick = random.random()
            if pick > ((X[2 * i - 1, j]) ** 2):
                XG[i, j] = X[2 * i - 1, j]
            else:
                XG[i, j] = X[2 * i, j]

    ## 计算初始适应度值
    fitness = []
    for i in range(NP):
        fitness.append(func(func_num, XG[i], 'value'))

    value_min = min(fitness)
    pos_min = fitness.index(value_min)

    gbest = value_min
    gbest_x = XG[pos_min]

    faimin = np.zeros((Gm,D))
    faimin[0] = fai[pos_min]
    # print(faimin)

    value_min = 0
    pos_min = 0

    while G < Gm:
        if G % IT==0 and G!=0:
            fp = min(Gmin)
            fcbest = Gmin[G]

            GP = (Gm-G+1.0)/Gm
            ATinit = 0.02

            IP = (fp-fcbest) / (fp+epsilon)     # 相对适应度改进水平
            DT = (IT/10000)*GP      # 判断是否需要干涉的阈值

            AP = 0
            if IP < DT: # 需要干扰
                for i in range(NP):
                    for j in range(D):
                        AP = AP + np.abs(XG[i][j]-gbest_x[j]) / (xmax-xmin)
                AP = AP / (NP*D)
                AT = ATinit * GP

                disturb_flag = 1

                if AP < AT:     # 使用历史最优进行分散
                    XG,fai,gbest,gbest_x,gfaimin = Dis(NP,GP,value,XG,gbest,gbest_x,fai,gfaimin,xmax,xmin)
                else:      # NP,GP,AP,value,XG,fai,cbest,cbest_x,xmax,xmin,gbest,gbest_x
                    XG,fai,gbest,gbest_x,gfaimin = Agg(NP,GP,AP,value,XG,fai,gfaimin,Gmin[G],best_x[G],xmax,xmin,gbest,gbest_x)


        fainew1 = np.zeros((NP, D))

        for i in range(NP):
            ########################----变异操作----#################################
            # 产生j,k,p,p1,p2五个不同的数
            a = 1
            b = NP
            li = list(range(NP))
            random.shuffle(li)
            dx = li
            j = dx[0]
            k = dx[1]
            p = dx[2]
            p1 = dx[3]
            p2 = dx[4]
            if j == i:
                j = dx[5]
            elif k == i:
                k = dx[5]
            elif p == i:
                p = dx[5]
            elif p1 == i:
                p1 = dx[5]
            elif p2 == i:
                p2 = dx[5]

            ## 1
            # F0 = 0.5

            ## 2
            lf0 = 0.5
            F0 = lf0 * (math.exp(1 - Gm / (Gm + 1 - G)))
            # F0 = lf0 * (2 ** (math.exp(1-Gm/(Gm+1-G))))


            ## 3
            # x00 = random.random()
            # F0 = 2 / math.sqrt(3) * (math.pi ** (-1/4)) * (1-x00 ** 2) * math.exp(-x00 ** 2/2)
            # if not disturb_flag:
            #     fainew1[i] = gfaimin + F0 * (fai[j] - fai[k])  # F0:变异因子
            # else:
            #     fainew1[i] = fai[p] + F0 * (fai[j] - fai[k])

            # if G<400:
            #     fainew1[i] = gfaimin + F0 * (fai[j] - fai[k]) 
            # else:
            #     fainew1[i] = faimin[G] + F0 * (fai[j] - fai[k])
            fainew1[i] = gfaimin + F0 * (fai[j] - fai[k]) 

        for i in range(NP):
            for j in range(D):
                if fainew1[i, j] > (- math.pi) and fainew1[i, j] < (math.pi):
                    fainew1[i, j] = fainew1[i, j]
                else:
                    fainew1[i, j] = math.pi * random.random()

        ## 量子交叉
        fainew2 = np.zeros((NP, D))
        for i in range(NP):
            for j in range(D):
                if random.random() < CR or (j == random.randint(1, D)):
                    fainew2[i, j] = fainew1[i, j]
                else:
                    fainew2[i, j] = fai[i, j]

        Xnew2 = np.zeros((2 * NP, D))
        for i in range(0, 2 * NP, 2):
            for j in range(D):
                Xnew2[i, j] = math.cos(fainew2[int(0.5 * (i + 1)), j])
                Xnew2[i + 1, j] = math.sin(fainew2[int(0.5 * (i + 1)), j])
        ## 量子空间变换
        Xnew2 = (span * Xnew2 + span2) / 2

        ## 量子测量
        XGnew2 = np.zeros((NP, D))
        for i in range(NP):
            for j in range(D):
                pick = random.random()
                if pick > ((Xnew2[2 * i - 1, j]) ** 2):
                    XGnew2[i, j] = Xnew2[2 * i - 1, j]
                else:
                    XGnew2[i, j] = Xnew2[2 * i, j]

        ##################----选择操作----####################
        fainew3 = np.zeros((NP, D))
        for i in range(NP):
            if func(func_num, XGnew2[i], 'value') < func(func_num, XG[i], 'value'):
                fainew3[i] = fainew2[i]
            else:
                fainew3[i] = fai[i]

        ## 量子染色体非门
        Pm_rand = np.random.rand(NP, 1)

        fainew4 = np.zeros((NP, D))
        for i in range(NP):
            for j in range(D):
                if (Pm[i] > Pm_rand[i]):
                    fainew4[i, j] = 0.5 * np.pi - fainew3[i, j]
                else:
                    fainew4[i, j] = fainew3[i, j]
        ##
        Xnew4 = np.zeros((2 * NP, D))
        for i in range(0, 2 * NP, 2):
            for j in range(D):
                Xnew4[i, j] = math.cos(fainew4[int(0.5 * (i + 1)), j])
                Xnew4[i + 1, j] = math.sin(fainew4[int(0.5 * (i + 1)), j])
        # 量子空间变换
        Xnew4 = (span * Xnew4 + span2) / 2
        # 量子测量
        XGnew4 = np.zeros((NP, D))
        for i in range(NP):
            for j in range(D):
                pick = random.random()
                if pick > (Xnew4[2 * i - 1, j]) ** 2:
                    XGnew4[i, j] = Xnew4[2 * i - 1, j]
                else:
                    XGnew4[i, j] = Xnew4[2 * i, j]
        ##
        
        # 找出最小值
        for i in range(NP):
            value[i] = func(func_num, XGnew4[i], 'value')

        value_min = min(value)
        pos_min = value.tolist().index(value_min)

        if value_min < gbest:
            gbest = value_min
            gbest_x = XGnew4[pos_min]
            gfaimin = fainew4[pos_min]


        # 第G代中的目标函数的最小值
        Gmin[G] = value_min
        # 保存最优的个体
        # print(G)
        faimin[G] = fainew4[pos_min]
        best_x[G] = XGnew4[pos_min]
        XG = XGnew4
        fai = fainew4

        # print(value_min)
        trace.append(value_min)
        # trace(G,2) = value_min;
        # print('正在计算函数，第',Clock,'次','第',G,'代')

        if G % 100 == 0:
            print('正在计算函数',func_num,'，第', Clock, '次', '第', G, '代')
            print(value_min)

        G += 1

        # print('正在计算函数%s , 运行第%d次 ， 第%d代\n' ,char(func) , Clock , G)

    time_end = time.time()
    a = time_end - time_start

    return gbest,trace, a

import pandas as pd
Gm = 10000
NP = 50
CR = 0.9
D = 1000
IT = 100000

func_nums = [1, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 20, 21, 22, 24, 25, 26, 28, 30]
hyfunc_nums = [21, 22, 24, 25, 26, 28, 30]

# for func_num in func_nums:
#     res = []
#     for i in range(25):
#         best,trace,a = QDE_PSE(func_num,Gm,NP,CR,D,i+1,IT)
#         res.append(trace)
#         df = pd.DataFrame(res)
#         print('最优值为：',best)
#         x = list(range(Gm))
#         y = trace
#     df.to_csv('res/QDE_PSE/'+str(D)+'/QDE_PSE_fun(d='+str(D)+')' + str(func_num) + '.csv', index=False)


func_num = 12
res = []
best_res = []
for i in range(1):
    best,trace,a = QDE_PSE(func_num,Gm,NP,CR,D,i+1,IT)
    res.append(trace)
    best_res.append(best)
    df = pd.DataFrame(res)
    print('最优值为：',best)
    x = list(range(Gm))
    y = trace
    # plt.plot(x, y)
    # plt.show()
print(df)


