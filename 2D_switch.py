import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import ArtistAnimation # アニメーション作成のためのメソッドをインポート
import matplotlib.animation as animation
import os
from numba.decorators import jit

#初期状態読む込み先
#init_file="/Users/sakailab/Desktop/nagata/1.27/4-elec_6V_iter=100"
init_file="/Users/zenyanagata/Desktop/test2"
#読み込むNdファイル名
filename="Nd_iter=10.txt"
#保存先
#active_file="/Users/sakailab/Desktop/nagata/1.27/4-elec_6V_switching_iter=25"
active_file = "/Users/zenyanagata/Desktop/test3"

resolution = 2/3.5

delta_Lx=1/resolution# 1 meshgridの長さ
delta_Ly=1/resolution

Lx = int(384*resolution - 1)#分解能をあげるときは、384*○ - 1の形にする
Ly = int(384*resolution - 1)#必ずLx = Lyで。

V = 6.0 # 電圧
end_iter = 25
n_iter = 0
turn = 1 #V_sweep用の変数
v_sweep_list = np.array([])
n_iter_list = np.array([])
def V_sweep(start_voltage, end_voltage):#sweep_rate = V/iter
    global V, turn, v_sweep_list, n_iter, n_iter_list
    sweep_rate = (end_voltage - start_voltage)/(end_iter/2)
    if n_iter == 0:
        V = start_voltage - sweep_rate*turn
    if V >= end_voltage:
        turn = -1
    V += sweep_rate*turn
    v_sweep_list = np.append(v_sweep_list, V)
    n_iter_list = np.append(n_iter_list, n_iter)

k = 1.38064852*10**(-23) #J/K
e = 1.60217662*10**(-19) #C
T = 300 #K
D = 4.6*10**-10*10**8 #μm^2/s
u = 2*e*D/(k*T) #0.01μm^2/Vs
eps = 85*10**(-6) #F/μm

phi_Bound= np.ones([Lx+1,Ly+1])
phi = np.zeros([Lx+1,Ly+1])
Ex = np.zeros([Lx+1,Ly+1])
Ey = np.zeros([Lx+1,Ly+1])
Nd = np.ones([Lx+1, Ly+1])*10**6

os.chdir(init_file)
def import_data():
    global Nd
    Nd = np.genfromtxt(filename, delimiter=" ")
import_data()

os.chdir(active_file)

# 電極の形状
def two_electrode_p():
    global phi_Bound
    #右上辺
    for i in range(int((Lx+1)/2-58*(Ly+1)/384), int((Lx+1)/2+1)):
        for j in range(int(108*(Lx+1)/384), int(108*(Lx+1)/384+i-((Lx+1)/2-58*(Ly+1)/384-1))):
            phi_Bound[i,j] = V
        for j in range(int(275*(Lx+1)/384), int(275*(Lx+1)/384+i-((Lx+1)/2-58*(Ly+1)/384-1))):
            phi_Bound[i,j] = 0
    #右下辺
    for i in range(int((Lx+1)/2), int((Lx+1)/2+59*(Ly+1)/384)):
        for j in range(int(108*(Lx+1)/384), int(167*(Lx+1)/384+((Lx+1)/2-i))):
            phi_Bound[i,j] = V
        for j in range(int(275*(Lx+1)/384), int(334*(Lx+1)/384+((Lx+1)/2-i))):
            phi_Bound[i,j] = 0
    #左上辺
    for i in range(int((Lx+1)/2-58*(Ly+1)/384), int((Lx+1)/2+1)):
        for j in range(int(109*(Lx+1)/384-i+(Lx+1)/2-58*(Ly+1)/384-1), int(109*(Lx+1)/384)):
            phi_Bound[i,j] = V
        for j in range(int(276*(Lx+1)/384-i+(Lx+1)/2-58*(Ly+1)/384-1), int(276*(Lx+1)/384)):
            phi_Bound[i,j] = 0
    #左下辺
    for i in range(int((Lx+1)/2), int((Lx+1)/2+59*(Ly+1)/384)):
        for j in range(int(50*(Lx+1)/384-((Lx+1)/2-i)), int(109*(Lx+1)/384)):
            phi_Bound[i,j] = V
        for j in range(int(217*(Lx+1)/384-((Lx+1)/2-i)), int(276*(Lx+1)/384)):
            phi_Bound[i,j] = 0

    np.savetxt("phi_Bound.txt", phi_Bound, delimiter=" ", fmt="%1.8f")#parallel electrode 4-2
#two_electrode_p()

def two_electrode_v():
    global phi_Bound
    #右上辺
    for i in range(int((Lx+1)/2-58*(Ly+1)/384), int((Lx+1)/2+1)):
        for j in range(int(108*(Lx+1)/384), int(108*(Lx+1)/384+i-((Lx+1)/2-58*(Ly+1)/384-1))):
            phi_Bound[j,i] = V
        for j in range(int(275*(Lx+1)/384), int(275*(Lx+1)/384+i-((Lx+1)/2-58*(Ly+1)/384-1))):
            phi_Bound[j,i] = 0
    #右下辺
    for i in range(int((Lx+1)/2), int((Lx+1)/2+59*(Ly+1)/384)):
        for j in range(int(108*(Lx+1)/384), int(167*(Lx+1)/384+((Lx+1)/2-i))):
            phi_Bound[j,i] = V
        for j in range(int(275*(Lx+1)/384), int(334*(Lx+1)/384+((Lx+1)/2-i))):
            phi_Bound[j,i] = 0
    #左上辺
    for i in range(int((Lx+1)/2-58*(Ly+1)/384), int((Lx+1)/2+1)):
        for j in range(int(109*(Lx+1)/384-i+(Lx+1)/2-58*(Ly+1)/384-1), int(109*(Lx+1)/384)):
            phi_Bound[j,i] = V
        for j in range(int(276*(Lx+1)/384-i+(Lx+1)/2-58*(Ly+1)/384-1), int(276*(Lx+1)/384)):
            phi_Bound[j,i] = 0
    #左下辺
    for i in range(int((Lx+1)/2), int((Lx+1)/2+59*(Ly+1)/384)):
        for j in range(int(50*(Lx+1)/384-((Lx+1)/2-i)), int(109*(Lx+1)/384)):
            phi_Bound[j,i] = V
        for j in range(int(217*(Lx+1)/384-((Lx+1)/2-i)), int(276*(Lx+1)/384)):
            phi_Bound[j,i] = 0

    np.savetxt("phi_Bound.txt", phi_Bound, delimiter=" ", fmt="%1.8f")#vertical electrode 1-3
two_electrode_v()

def four_electrode():
    global phi_Bound
    #右上辺
    for i in range(int((Lx+1)/2-58*(Ly+1)/384), int((Lx+1)/2+1)):
        for j in range(int(108*(Lx+1)/384), int(108*(Lx+1)/384+i-((Lx+1)/2-58*(Ly+1)/384-1))):
            phi_Bound[i,j] = V
        for j in range(int(275*(Lx+1)/384), int(275*(Lx+1)/384+i-((Lx+1)/2-58*(Ly+1)/384-1))):
            phi_Bound[i,j] = V
    #右下辺
    for i in range(int((Lx+1)/2), int((Lx+1)/2+59*(Ly+1)/384)):
        for j in range(int(108*(Lx+1)/384), int(167*(Lx+1)/384+((Lx+1)/2-i))):
            phi_Bound[i,j] = V
        for j in range(int(275*(Lx+1)/384), int(334*(Lx+1)/384+((Lx+1)/2-i))):
            phi_Bound[i,j] = V
    #左上辺
    for i in range(int((Lx+1)/2-58*(Ly+1)/384), int((Lx+1)/2+1)):
        for j in range(int(109*(Lx+1)/384-i+(Lx+1)/2-58*(Ly+1)/384-1), int(109*(Lx+1)/384)):
            phi_Bound[i,j] = V
        for j in range(int(276*(Lx+1)/384-i+(Lx+1)/2-58*(Ly+1)/384-1), int(276*(Lx+1)/384)):
            phi_Bound[i,j] = V
    #左下辺
    for i in range(int((Lx+1)/2), int((Lx+1)/2+59*(Ly+1)/384)):
        for j in range(int(50*(Lx+1)/384-((Lx+1)/2-i)), int(109*(Lx+1)/384)):
            phi_Bound[i,j] = V
        for j in range(int(217*(Lx+1)/384-((Lx+1)/2-i)), int(276*(Lx+1)/384)):
            phi_Bound[i,j] = V
    #右上
    for i in range(int((Lx+1)/2-58*(Ly+1)/384), int((Lx+1)/2+1)):
        for j in range(int(108*(Lx+1)/384), int(108*(Lx+1)/384+i-((Lx+1)/2-58*(Ly+1)/384-1))):
            phi_Bound[j,i] = 0
        for j in range(int(275*(Lx+1)/384), int(275*(Lx+1)/384+i-((Lx+1)/2-58*(Ly+1)/384-1))):
            phi_Bound[j,i] = 0
    #右下辺
    for i in range(int((Lx+1)/2), int((Lx+1)/2+59*(Ly+1)/384)):
        for j in range(int(108*(Lx+1)/384), int(167*(Lx+1)/384+((Lx+1)/2-i))):
            phi_Bound[j,i] = 0
        for j in range(int(275*(Lx+1)/384), int(334*(Lx+1)/384+((Lx+1)/2-i))):
            phi_Bound[j,i] = 0
    #左上辺
    for i in range(int((Lx+1)/2-58*(Ly+1)/384), int((Lx+1)/2+1)):
        for j in range(int(109*(Lx+1)/384-i+(Lx+1)/2-58*(Ly+1)/384-1), int(109*(Lx+1)/384)):
            phi_Bound[j,i] = 0
        for j in range(int(276*(Lx+1)/384-i+(Lx+1)/2-58*(Ly+1)/384-1), int(276*(Lx+1)/384)):
            phi_Bound[j,i] = 0
    #左下辺
    for i in range(int((Lx+1)/2), int((Lx+1)/2+59*(Ly+1)/384)):
        for j in range(int(50*(Lx+1)/384-((Lx+1)/2-i)), int(109*(Lx+1)/384)):
            phi_Bound[j,i] = 0
        for j in range(int(217*(Lx+1)/384-((Lx+1)/2-i)), int(276*(Lx+1)/384)):
            phi_Bound[j,i] = 0

    np.savetxt("phi_Bound.txt", phi_Bound, delimiter=" ", fmt="%1.8f")
#four_electrode()

#静電ポテンシャル分布の計算
@jit
def phi_culc(phi):
    convegence_criterion = 10**-5
    # for SOR method
    aa_recta=0.5*(np.cos(np.pi/Lx)+np.cos(np.pi/Ly))
    omega_SOR_recta = 2/(1+np.sqrt(1-aa_recta**2)) # SOR法における加速パラメータ
    #print("omega_SOR_rect=", omega_SOR_recta)

    #main culculation
    delta = 1.0
    n_iter=0
    while delta > convegence_criterion:
        phi_in=phi.copy()
        if n_iter % 100 ==0:  # 収束状況のモニタリング
            print("iteration No =", n_iter, "delta=",delta)
        for i in range(Lx+1):
            for j in range(Ly+1):
                if phi_Bound[i, j] == 0 or phi_Bound[i, j] == V:
                    phi[i,j] = phi_Bound[i,j]
                elif i == 0:
                    if j == 0:              #角の場合、隣り合うグリッドは2点
                        phi[i,j] = phi[i,j] + omega_SOR_recta * ((phi[i+1,j] + phi[i,j+1])/2 - phi[i,j] + 2*e*Nd[i,j]*delta_Lx*delta_Ly/(2*eps))
                    elif j == Ly:
                        phi[i,j] = phi[i,j] + omega_SOR_recta * ((phi[i+1,j] + phi[i,j-1])/2 - phi[i,j] + 2*e*Nd[i,j]*delta_Lx*delta_Ly/(2*eps))
                    else:                   #辺の場合、隣り合うグリッドは3点
                        phi[i,j] = phi[i,j] + omega_SOR_recta * ((phi[i+1,j]+phi[i,j-1]+phi[i,j+1])/3 - phi[i,j] + 2*e*Nd[i,j]*delta_Lx*delta_Ly/(3*eps))
                elif i == Lx:
                    if j==0:                #角の場合、隣り合うグリッドは2点
                        phi[i,j] = phi[i,j] + omega_SOR_recta * ((phi[i-1,j] + phi[i,j+1])/2 - phi[i,j] + 2*e*Nd[i,j]*delta_Lx*delta_Ly/(2*eps))
                    elif j == Ly:
                        phi[i,j] = phi[i,j] + omega_SOR_recta * ((phi[i-1,j] + phi[i,j-1])/2 - phi[i,j] + 2*e*Nd[i,j]*delta_Lx*delta_Ly/(2*eps))
                    else:                   #辺の場合、隣り合うグリッドは3点
                        phi[i,j] = phi[i,j] + omega_SOR_recta * ((phi[i-1,j]+phi[i,j-1]+phi[i,j+1])/3 - phi[i,j] + 2*e*Nd[i,j]*delta_Lx*delta_Ly/(3*eps))
                elif j == 0:
                    phi[i,j] = phi[i,j] + omega_SOR_recta * ((phi[i+1,j]+phi[i-1,j]+phi[i,j+1])/3 - phi[i,j] + 2*e*Nd[i,j]*delta_Lx*delta_Ly/(3*eps))
                elif j == Lx:
                    phi[i,j] = phi[i,j] + omega_SOR_recta * ((phi[i-1,j]+phi[i+1,j]+phi[i,j-1])/3 - phi[i,j] + 2*e*Nd[i,j]*delta_Lx*delta_Ly/(3*eps))
                else:
                    phi[i,j] = phi[i,j]+omega_SOR_recta *((phi[i+1,j] + phi[i-1,j] + phi[i,j+1] + phi[i,j-1])/4-phi[i,j] + 2*e*Nd[i,j]*delta_Lx*delta_Ly/(4*eps)) # SOR = ガウス-ザイデル + OR

        m=0
        diff_phi = np.empty((Lx+1)*(Ly+1), dtype=np.float64)
        for i in range(Ly+1):
            for j in range(Lx+1):
                diff_phi[m] = abs(phi[i,j]-phi_in[i,j])
                m+=1
        delta = np.max(diff_phi)
        n_iter+=1
    return phi
phi = phi_culc(phi)

#静電場の計算
@jit
def field_culc(Ex, Ey):
    for i in range(Lx+1):
        for j in range(Ly+1):
            if i==0:
                if j==0:
                    Ey[i,j] = (phi[i+1,j]-phi[i,j])/delta_Ly
                    Ex[i,j] = -(phi[i,j+1]-phi[i,j])/delta_Lx
                elif j==Lx:
                    Ey[i,j] = (phi[i+1,j]-phi[i,j])/delta_Ly
                    Ex[i,j] = -(phi[i,j]-phi[i,j-1])/delta_Lx
                else:
                    Ey[i,j] = (phi[i+1,j]-phi[i,j])/delta_Ly
                    Ex[i,j] = (-(phi[i,j+1]-phi[i,j])/delta_Lx-(phi[i,j]-phi[i,j-1])/delta_Lx)/2
            elif i==Ly:
                if j==0:
                    Ey[i,j] = (phi[i,j]-phi[i-1,j])/delta_Ly
                    Ex[i,j] = -(phi[i,j+1]-phi[i,j])/delta_Lx
                elif j==Lx:
                    Ey[i,j] = (phi[i,j]-phi[i-1,j])/delta_Ly
                    Ex[i,j] = -(phi[i,j]-phi[i,j-1])/delta_Lx
                else:
                    Ey[i,j] = (phi[i,j]-phi[i-1,j])/delta_Ly
                    Ex[i,j] = (-(phi[i,j+1]-phi[i,j])/delta_Lx-(phi[i,j]-phi[i,j-1])/delta_Lx)/2
            elif j==0:
                Ey[i,j] = ((phi[i+1,j]-phi[i,j])/delta_Ly+(phi[i,j]-phi[i-1,j])/delta_Ly)/2
                Ex[i,j] = -(phi[i,j+1]-phi[i,j])/delta_Lx
            elif j==Lx:
                Ey[i,j] = ((phi[i+1,j]-phi[i,j])/delta_Ly+(phi[i,j]-phi[i-1,j])/delta_Ly)/2
                Ex[i,j] = -(phi[i,j]-phi[i,j-1])/delta_Lx
            else:
                Ey[i,j] = ((phi[i+1,j]-phi[i,j])/delta_Ly+(phi[i,j]-phi[i-1,j])/delta_Ly)/2 #前後のgridで微分*1/2
                Ex[i,j] = (-(phi[i,j+1]-phi[i,j])/delta_Lx-(phi[i,j]-phi[i,j-1])/delta_Lx)/2 #前後のgridで微分*1/2
    return Ex, Ey
Ex, Ey = field_culc(Ex,Ey)

#ドリフト・拡散シミュレーション
delta_Ndxp = np.zeros([Lx+1, Ly+1])
delta_Ndxn = np.zeros([Lx+1, Ly+1])
delta_Ndyp = np.zeros([Lx+1, Ly+1])
delta_Ndyn = np.zeros([Lx+1, Ly+1])
Drift_x = np.zeros([Lx+1, Ly+1])
Drift_y = np.zeros([Lx+1, Ly+1])
Diff_xp = np.zeros([Lx+1, Ly+1])
Diff_xn = np.zeros([Lx+1, Ly+1])
Diff_yp = np.zeros([Lx+1, Ly+1])
Diff_yn = np.zeros([Lx+1, Ly+1])

#初期状態
@jit
def init_culc(Nd, Drift_x, Drift_y, Diff_xn, Diff_xp, Diff_yn, Diff_yp):
    for i in range(Lx+1):
        for j in range(Ly+1):
            if phi_Bound[i, j] == 0 or phi_Bound[i, j] == V:
                Nd[i,j] = 0

    for i in range(Lx+1):
        for j in range(Ly+1):
            Drift_x[i,j] = u * Nd[i,j] * Ex[i,j]
            Drift_y[i,j] = u * Nd[i,j] * Ey[i,j]

            Diff_xp[i,j] = D * delta_Ndxp[i,j]
            Diff_xn[i,j] = D * delta_Ndxn[i,j]
            Diff_yp[i,j] = D * delta_Ndyp[i,j]
            Diff_yn[i,j] = D * delta_Ndyn[i,j]

    return Nd, Drift_x, Drift_y, Diff_xn, Diff_xp, Diff_yn, Diff_yp
Nd, Drift_x, Drift_y, Diff_xn, Diff_xp, Diff_yn, Diff_yp = init_culc(Nd, Drift_x, Drift_y, Diff_xn, Diff_xp, Diff_yn, Diff_yp)


#for plot
fig2 = plt.figure()
anim2 = []
im2 = plt.imshow(Nd,cmap='rainbow', vmax=0.5*10**7, vmin=0)
anim2.append([im2])
np.savetxt("Nd_iter=" + np.str(n_iter) + ".txt", Nd, delimiter=" ", fmt="%1.8f")

@jit
def drift_culc(Nd):
    for i in range(0, Lx+1):
        for j in range(0, Ly+1):
            move_x, move_y = abs(Drift_x[i,j]), abs(Drift_y[i,j])
            if Nd[i,j] == 0:
                move_x, move_y = 0, 0
            elif 0 < Nd[i,j] < move_x + move_y:
                move_x, move_y = Nd[i,j]*move_x/(move_x+move_y), Nd[i,j]*move_y/(move_x+move_y)
            elif Nd[i,j] < 0:
                print("Nd[",i,j,"] < 0")
            #電極内部でのキャリア移動はなし
            if phi_Bound[i,j] == 0 or phi_Bound[i,j] == V:
                pass

            #界面
            elif i==0:
                if j==0:
                    if Drift_x[i,j] > 0:
                        Nd[i,j], Nd[i,j+1] = Nd[i,j] - move_x, Nd[i,j+1] + move_x
                    if Drift_y[i,j] < 0:
                        Nd[i,j], Nd[i+1,j] = Nd[i,j] - move_y, Nd[i+1,j] + move_y
                elif j==Lx:
                    if Drift_x[i,j] < 0:
                        Nd[i,j], Nd[i,j-1] = Nd[i,j] - move_x, Nd[i,j-1] + move_x
                    if Drift_y[i,j] < 0:
                        Nd[i,j], Nd[i+1,j] = Nd[i,j] - move_y, Nd[i+1,j] + move_y
                else:
                    if Drift_x[i,j] > 0:
                        Nd[i,j], Nd[i,j+1] = Nd[i,j] - move_x, Nd[i,j+1] + move_x
                    elif Drift_x[i,j] < 0:
                        Nd[i,j], Nd[i,j-1] = Nd[i,j] - move_x, Nd[i,j-1] + move_x
                    if Drift_y[i,j] < 0:
                        Nd[i,j], Nd[i+1,j] = Nd[i,j] - move_y, Nd[i+1,j] + move_y
            elif i==Ly:
                if j==0:
                    if Drift_x[i,j] > 0:
                        Nd[i,j], Nd[i,j+1] = Nd[i,j] - move_x, Nd[i,j+1] + move_x
                    if Drift_y[i,j] > 0:
                        Nd[i,j], Nd[i-1,j] = Nd[i,j] - move_y, Nd[i-1,j] + move_y
                elif j==Lx:
                    if Drift_x[i,j] < 0:
                        Nd[i,j], Nd[i,j-1] = Nd[i,j] - move_x, Nd[i,j-1] + move_x
                    if Drift_y[i,j] > 0:
                        Nd[i,j], Nd[i-1,j] = Nd[i,j] - move_y, Nd[i-1,j] + move_y
                else:
                    if Drift_x[i,j] > 0:
                        Nd[i,j], Nd[i,j+1] = Nd[i,j] - move_x, Nd[i,j+1] + move_x
                    elif Drift_x[i,j] < 0:
                        Nd[i,j], Nd[i,j-1] = Nd[i,j] - move_x, Nd[i,j-1] + move_x
                    if Drift_y[i,j] > 0:
                        Nd[i,j], Nd[i-1,j] = Nd[i,j] - move_y, Nd[i-1,j] + move_y
            elif j==0:
                if Drift_x[i,j] > 0:
                    Nd[i,j], Nd[i,j+1] = Nd[i,j] - move_x, Nd[i,j+1] + move_x
                if Drift_y[i,j] > 0:
                    Nd[i,j], Nd[i-1,j] = Nd[i,j] - move_y, Nd[i-1,j] + move_y
                elif Drift_y[i,j] < 0:
                    Nd[i,j], Nd[i+1,j] = Nd[i,j] - move_y, Nd[i+1,j] + move_y
            elif j==Lx:
                if Drift_x[i,j] < 0:
                    Nd[i,j], Nd[i,j-1] = Nd[i,j] - move_x, Nd[i,j-1] + move_x
                if Drift_y[i,j] > 0:
                    Nd[i,j], Nd[i-1,j] = Nd[i,j] - move_y, Nd[i-1,j] + move_y
                elif Drift_y[i,j] < 0:
                    Nd[i,j], Nd[i+1,j] = Nd[i,j] - move_y, Nd[i+1,j] + move_y

            #Bias側での境界条件（電極内にキャリアは入れない）
            elif phi_Bound[i,j-1] == V and phi_Bound[i+1,j] == V:#Bias側右上辺
                if Drift_x[i,j] > 0:
                    Nd[i,j], Nd[i,j+1] = Nd[i,j] - move_x, Nd[i,j+1] + move_x
                if Drift_y[i,j] > 0:
                    Nd[i,j], Nd[i-1,j] = Nd[i,j] - move_y, Nd[i-1,j] + move_y
            elif phi_Bound[i,j-1] == V and phi_Bound[i-1,j] == V:#Bias側右下辺
                if Drift_x[i,j] > 0:
                    Nd[i,j], Nd[i,j+1] = Nd[i,j] - move_x, Nd[i,j+1] + move_x
                if Drift_y[i,j] < 0:
                    Nd[i,j], Nd[i+1,j] = Nd[i,j] - move_y, Nd[i+1,j] + move_y
            elif phi_Bound[i,j+1] == V and phi_Bound[i+1,j] == V:#Bias側左上辺
                if Drift_x[i,j] < 0:
                    Nd[i,j], Nd[i,j-1] = Nd[i,j] - move_x, Nd[i,j-1] + move_x
                if Drift_y[i,j] > 0:
                    Nd[i,j], Nd[i-1,j] = Nd[i,j] - move_y, Nd[i-1,j] + move_y
            elif phi_Bound[i,j+1] == V and phi_Bound[i-1,j] == V:#Bias側左下辺
                if Drift_x[i,j] < 0:
                    Nd[i,j], Nd[i,j-1] = Nd[i,j] - move_x, Nd[i,j-1] + move_x
                if Drift_y[i,j] < 0:
                    Nd[i,j], Nd[i+1,j] = Nd[i,j] - move_y, Nd[i+1,j] + move_y
            #GND側での境界条件（電極内にキャリアは入れない）
            elif phi_Bound[i,j+1] == 0 and phi_Bound[i+1,j] == 0:#GND側の左上辺
                if Drift_x[i,j] < 0:
                    Nd[i,j], Nd[i,j-1] = Nd[i,j] - move_x, Nd[i,j-1] + move_x
                if Drift_y[i,j] > 0:
                    Nd[i,j], Nd[i-1,j] = Nd[i,j] - move_y, Nd[i-1,j] + move_y
            elif phi_Bound[i,j+1] == 0 and phi_Bound[i-1,j] == 0:#GND側の左下辺
                if Drift_x[i,j] < 0:
                    Nd[i,j], Nd[i,j-1] = Nd[i,j] - move_x, Nd[i,j-1] + move_x
                if Drift_y[i,j] < 0:
                    Nd[i,j], Nd[i+1,j] = Nd[i,j] - move_y, Nd[i+1,j] + move_y
            elif phi_Bound[i,j-1] == 0 and phi_Bound[i+1,j] == 0:#GND側の右上辺
                if Drift_x[i,j] > 0:
                    Nd[i,j], Nd[i,j+1] = Nd[i,j] - move_x, Nd[i,j+1] + move_x
                if Drift_y[i,j] > 0:
                    Nd[i,j], Nd[i-1,j] = Nd[i,j] - move_y, Nd[i-1,j] + move_y
            elif phi_Bound[i,j-1] == 0 and phi_Bound[i-1,j] == 0:#GND側右下辺
                if Drift_x[i,j] > 0:
                    Nd[i,j], Nd[i,j+1] = Nd[i,j] - move_x, Nd[i,j+1] + move_x
                if Drift_y[i,j] < 0:
                    Nd[i,j], Nd[i+1,j] = Nd[i,j] - move_y, Nd[i+1,j] + move_y
            #Tipでの境界条件
            elif phi_Bound[i,j+1] == V or phi_Bound[i,j+1] == 0:#Bias,GNDの左Tip
                if Drift_x[i,j] < 0:
                    Nd[i,j], Nd[i,j-1] = Nd[i,j] - move_x, Nd[i,j-1] + move_x
                if Drift_y[i,j] > 0:
                    Nd[i,j], Nd[i-1,j] = Nd[i,j] - move_y, Nd[i-1,j] + move_y
                elif Drift_y[i,j] < 0:
                    Nd[i,j], Nd[i+1,j] = Nd[i,j] - move_y, Nd[i+1,j] + move_y
            elif phi_Bound[i,j-1] == V or phi_Bound[i,j-1] == 0:#Bias,GNDの右Tip
                if Drift_x[i,j] > 0:
                    Nd[i,j], Nd[i,j+1] = Nd[i,j] - move_x, Nd[i,j+1] + move_x
                if Drift_y[i,j] > 0:
                    Nd[i,j], Nd[i-1,j] = Nd[i,j] - move_y, Nd[i-1,j] + move_y
                elif Drift_y[i,j] < 0:
                    Nd[i,j], Nd[i+1,j] = Nd[i,j] - move_y, Nd[i+1,j] + move_y
            elif phi_Bound[i+1,j] == V or phi_Bound[i+1,j] == 0:#Bias,GNDの上Tip
                if Drift_x[i,j] > 0:
                    Nd[i,j], Nd[i,j+1] = Nd[i,j] - move_x, Nd[i,j+1] + move_x
                elif Drift_x[i,j] < 0:
                    Nd[i,j], Nd[i,j-1] = Nd[i,j] - move_x, Nd[i,j-1] + move_x
                if Drift_y[i,j] > 0:
                    Nd[i,j], Nd[i-1,j] = Nd[i,j] - move_y, Nd[i-1,j] + move_y
            elif phi_Bound[i-1,j] == V or phi_Bound[i-1,j] == 0:#Bias,GNDの下Tip
                if Drift_x[i,j] > 0:
                    Nd[i,j], Nd[i,j+1] = Nd[i,j] - move_x, Nd[i,j+1] + move_x
                elif Drift_x[i,j] < 0:
                    Nd[i,j], Nd[i,j-1] = Nd[i,j] - move_x, Nd[i,j-1] + move_x
                if Drift_y[i,j] < 0:
                    Nd[i,j], Nd[i+1,j] = Nd[i,j] - move_y, Nd[i+1,j] + move_y
            #境界以外
            else:
                if Drift_x[i,j] > 0:
                    Nd[i,j], Nd[i,j+1] = Nd[i,j] - move_x, Nd[i,j+1] + move_x
                elif Drift_x[i,j] < 0:
                    Nd[i,j], Nd[i,j-1] = Nd[i,j] - move_x, Nd[i,j-1] + move_x
                if Drift_y[i,j] > 0:
                    Nd[i,j], Nd[i-1,j] = Nd[i,j] - move_y, Nd[i-1,j] + move_y
                elif Drift_y[i,j] < 0:
                    Nd[i,j], Nd[i+1,j] = Nd[i,j] - move_y, Nd[i+1,j] + move_y
    return Nd

@jit
def diff_culc(Diff_xn, Diff_xp, Diff_yn, Diff_yp, Nd):
    for i in range(0, Lx+1):
        for j in range(0, Ly+1):
            if Nd[i,j] <= 0:
                if Diff_xp[i,j] > 0:
                    Diff_xp[i,j] = 0
                if Diff_xn[i,j] < 0:
                    Diff_xn[i,j] = 0
                if Diff_yp[i,j] > 0:
                    Diff_yp[i,j] = 0
                if Diff_yn[i,j] < 0:
                    Diff_yn[i,j] = 0
            elif Nd[i,j] <= Diff_xp[i,j] - Diff_xn[i,j] + Diff_yp[i,j] - Diff_yn[i,j]:#流出 - 流入
                if Diff_xp[i,j] > 0:
                    Diff_xp[i,j] = 0
                if Diff_xn[i,j] < 0:
                    Diff_xn[i,j] = 0
                if Diff_yp[i,j] > 0:
                    Diff_yp[i,j] = 0
                if Diff_yn[i,j] < 0:
                    Diff_yn[i,j] = 0
            if Diff_xn[i,j]>0 and Nd[i,j-1]<=abs(Diff_xn[i,j]):
                Diff_xn[i,j] = 0
            if Diff_xp[i,j]<0 and Nd[i,j+1]<=abs(Diff_xp[i,j]):
                Diff_xp[i,j] = 0
            if Diff_yn[i,j]>0 and Nd[i+1,j]<=abs(Diff_yn[i,j]):
                Diff_yn[i,j] = 0
            if Diff_yp[i,j]<0 and Nd[i-1,j]<=abs(Diff_yp[i,j]):
                Diff_yp[i,j] = 0
            #電極内部でのキャリア移動はなし
            if phi_Bound[i,j] == 0 or phi_Bound[i,j] == V:
                pass

            #界面
            elif i==0:
                if j==0:
                    Nd[i,j], Nd[i,j+1], Nd[i+1,j] = Nd[i,j] - Diff_xp[i,j] + Diff_yn[i,j], Nd[i,j+1] + Diff_xp[i,j], Nd[i+1,j] - Diff_yn[i,j]
                elif j==Lx:
                    Nd[i,j], Nd[i,j-1], Nd[i+1,j] = Nd[i,j] + Diff_xn[i,j] + Diff_yn[i,j], Nd[i,j-1] - Diff_xn[i,j], Nd[i+1,j] - Diff_yn[i,j]
                else:
                    Nd[i,j], Nd[i,j+1], Nd[i,j-1], Nd[i+1,j] = Nd[i,j] - Diff_xp[i,j] + Diff_xn[i,j] + Diff_yn[i,j], Nd[i,j+1] + Diff_xp[i,j],\
                    Nd[i,j-1] - Diff_xn[i,j], Nd[i+1,j] - Diff_yn[i,j]
            elif i==Ly:
                if j==0:
                    Nd[i,j], Nd[i,j+1], Nd[i-1,j] = Nd[i,j] - Diff_xp[i,j] - Diff_yp[i,j], Nd[i,j+1] + Diff_xp[i,j], Nd[i-1,j] + Diff_yp[i,j]
                elif j==Lx:
                    Nd[i,j], Nd[i,j-1], Nd[i-1,j] = Nd[i,j] + Diff_xn[i,j] - Diff_yp[i,j], Nd[i,j-1] - Diff_xn[i,j], Nd[i-1,j] + Diff_yp[i,j]
                else:
                    Nd[i,j], Nd[i,j+1], Nd[i,j-1], Nd[i-1,j] = Nd[i,j] - Diff_xp[i,j] + Diff_xn[i,j] - Diff_yp[i,j], Nd[i,j+1] + Diff_xp[i,j],\
                    Nd[i,j-1] - Diff_xn[i,j], Nd[i-1,j] + Diff_yp[i,j]
            elif j==0:
                Nd[i,j], Nd[i,j+1], Nd[i-1,j], Nd[i+1,j] = Nd[i,j] - Diff_xp[i,j] - Diff_yp[i,j] + Diff_yn[i,j], Nd[i,j+1] + Diff_xp[i,j],\
                Nd[i-1,j] + Diff_yp[i,j], Nd[i+1,j] - Diff_yn[i,j]
            elif j==Lx:
                Nd[i,j], Nd[i,j-1], Nd[i-1,j], Nd[i+1,j] = Nd[i,j] + Diff_xn[i,j] - Diff_yp[i,j] + Diff_yn[i,j], Nd[i,j-1] - Diff_xn[i,j],\
                Nd[i-1,j] + Diff_yp[i,j], Nd[i+1,j] - Diff_yn[i,j]

            #Bias側での境界条件（電極内にキャリアは入れない）
            elif phi_Bound[i,j-1] == V and phi_Bound[i+1,j] == V:#Bias側右上辺
                Nd[i,j], Nd[i,j+1], Nd[i-1,j] = Nd[i,j] - Diff_xp[i,j] - Diff_yp[i,j], Nd[i,j+1] + Diff_xp[i,j], Nd[i-1,j] + Diff_yp[i,j]
            elif phi_Bound[i,j-1] == V and phi_Bound[i-1,j] == V:#Bias側右下辺
                Nd[i,j], Nd[i,j+1], Nd[i+1,j] = Nd[i,j] - Diff_xp[i,j] + Diff_yn[i,j], Nd[i,j+1] + Diff_xp[i,j], Nd[i+1,j] - Diff_yn[i,j]
            elif phi_Bound[i,j+1] == V and phi_Bound[i+1,j] == V:#Bias側左上辺
                Nd[i,j], Nd[i,j-1], Nd[i-1,j] = Nd[i,j] + Diff_xn[i,j] - Diff_yp[i,j], Nd[i,j-1] - Diff_xn[i,j], Nd[i-1,j] + Diff_yp[i,j]
            elif phi_Bound[i,j+1] == V and phi_Bound[i-1,j] == V:#Bias側左下辺
                Nd[i,j], Nd[i,j-1], Nd[i+1,j] = Nd[i,j] + Diff_xn[i,j] + Diff_yn[i,j], Nd[i,j-1] - Diff_xn[i,j], Nd[i+1,j] - Diff_yn[i,j]
            #GND側での境界条件（電極内にキャリアは入れない）
            elif phi_Bound[i,j-1] == 0 and phi_Bound[i+1,j] == 0:#Bias側右上辺
                Nd[i,j], Nd[i,j+1], Nd[i-1,j] = Nd[i,j] - Diff_xp[i,j] - Diff_yp[i,j], Nd[i,j+1] + Diff_xp[i,j], Nd[i-1,j] + Diff_yp[i,j]
            elif phi_Bound[i,j-1] == 0 and phi_Bound[i-1,j] == 0:#Bias側右下辺
                Nd[i,j], Nd[i,j+1], Nd[i+1,j] = Nd[i,j] - Diff_xp[i,j] + Diff_yn[i,j], Nd[i,j+1] + Diff_xp[i,j], Nd[i+1,j] - Diff_yn[i,j]
            elif phi_Bound[i,j+1] == 0 and phi_Bound[i+1,j] == 0:#Bias側左上辺
                Nd[i,j], Nd[i,j-1], Nd[i-1,j] = Nd[i,j] + Diff_xn[i,j] - Diff_yp[i,j], Nd[i,j-1] - Diff_xn[i,j], Nd[i-1,j] + Diff_yp[i,j]
            elif phi_Bound[i,j+1] == 0 and phi_Bound[i-1,j] == 0:#Bias側左下辺
                Nd[i,j], Nd[i,j-1], Nd[i+1,j] = Nd[i,j] + Diff_xn[i,j] + Diff_yn[i,j], Nd[i,j-1] - Diff_xn[i,j], Nd[i+1,j] - Diff_yn[i,j]
            #Tipでの境界条件
            elif phi_Bound[i,j+1] == V or phi_Bound[i,j+1] == 0:#Bias,GNDの左Tip
                Nd[i,j], Nd[i,j-1], Nd[i-1,j], Nd[i+1,j] = Nd[i,j] + Diff_xn[i,j] - Diff_yp[i,j] + Diff_yn[i,j], Nd[i,j-1] - Diff_xn[i,j],\
                Nd[i-1,j] + Diff_yp[i,j], Nd[i+1,j] - Diff_yn[i,j]
            elif phi_Bound[i,j-1] == V or phi_Bound[i,j-1] == 0:#Bias,GNDの右Tip
                Nd[i,j], Nd[i,j+1], Nd[i-1,j], Nd[i+1,j] = Nd[i,j] - Diff_xp[i,j] - Diff_yp[i,j] + Diff_yn[i,j], Nd[i,j+1] + Diff_xp[i,j],\
                Nd[i-1,j] + Diff_yp[i,j], Nd[i+1,j] - Diff_yn[i,j]
            elif phi_Bound[i+1,j] == V or phi_Bound[i+1,j] == 0:#Bias,GNDの上Tip
                Nd[i,j], Nd[i,j+1], Nd[i,j-1], Nd[i-1,j] = Nd[i,j] - Diff_xp[i,j] + Diff_xn[i,j] - Diff_yp[i,j], Nd[i,j+1] + Diff_xp[i,j],\
                Nd[i,j-1] - Diff_xn[i,j], Nd[i-1,j] + Diff_yp[i,j]
            elif phi_Bound[i-1,j] == V or phi_Bound[i-1,j] == 0:#Bias,GNDの下Tip
                Nd[i,j], Nd[i,j+1], Nd[i,j-1], Nd[i+1,j] = Nd[i,j] - Diff_xp[i,j] + Diff_xn[i,j] + Diff_yn[i,j], Nd[i,j+1] + Diff_xp[i,j],\
                Nd[i,j-1] - Diff_xn[i,j], Nd[i+1,j] - Diff_yn[i,j]
            #境界以外
            else:
                Nd[i,j], Nd[i,j+1], Nd[i,j-1], Nd[i-1,j], Nd[i+1,j] = Nd[i,j] - Diff_xp[i,j] - Diff_yp[i,j] + Diff_xn[i,j] + Diff_yn[i,j],\
                Nd[i,j+1] + Diff_xp[i,j], Nd[i,j-1] - Diff_xn[i,j], Nd[i-1,j] + Diff_yp[i,j], Nd[i+1,j] - Diff_yn[i,j]
    return Diff_xn, Diff_xp, Diff_yn, Diff_yp, Nd

@jit
def next_drift_diff_culc(Drift_x, Drift_y, Diff_xn, Diff_xp, Diff_yn, Diff_yp, delta_Ndxn, delta_Ndxp, delta_Ndyn, delta_Ndyp):
    for i in range(0, Lx+1):
        for j in range(0, Ly+1):
            if phi_Bound[i,j] == 0 or phi_Bound[i,j] == V:
                pass
            elif i==0:
                if j==0:
                    delta_Ndxp[i,j] = -(Nd[i,j+1]-Nd[i,j])/delta_Lx
                    delta_Ndxn[i,j] = 0
                    delta_Ndyp[i,j] = 0
                    delta_Ndyn[i,j] = (Nd[i+1,j]-Nd[i,j])/delta_Ly
                elif j==Lx:
                    delta_Ndxp[i,j] = 0
                    delta_Ndxn[i,j] = -(Nd[i,j]-Nd[i,j-1])/delta_Lx
                    delta_Ndyp[i,j] = 0
                    delta_Ndyn[i,j] = (Nd[i+1,j]-Nd[i,j])/delta_Ly
                else:
                    delta_Ndxp[i,j] = -(Nd[i,j+1]-Nd[i,j])/delta_Lx
                    delta_Ndxn[i,j] = -(Nd[i,j]-Nd[i,j-1])/delta_Lx
                    delta_Ndyp[i,j] = 0
                    delta_Ndyn[i,j] = (Nd[i+1,j]-Nd[i,j])/delta_Ly
            elif i==Ly:
                if j==0:
                    delta_Ndxp[i,j] = -(Nd[i,j+1]-Nd[i,j])/delta_Lx
                    delta_Ndxn[i,j] = 0
                    delta_Ndyp[i,j] = (Nd[i,j]-Nd[i-1,j])/delta_Ly
                    delta_Ndyn[i,j] = 0
                elif j==Lx:
                    delta_Ndxp[i,j] = 0
                    delta_Ndxn[i,j] = -(Nd[i,j]-Nd[i,j-1])/delta_Lx
                    delta_Ndyp[i,j] = (Nd[i,j]-Nd[i-1,j])/delta_Ly
                    delta_Ndyn[i,j] = 0
                else:
                    delta_Ndxp[i,j] = -(Nd[i,j+1]-Nd[i,j])/delta_Lx
                    delta_Ndxn[i,j] = -(Nd[i,j]-Nd[i,j-1])/delta_Lx
                    delta_Ndyp[i,j] = (Nd[i,j]-Nd[i-1,j])/delta_Ly
                    delta_Ndyn[i,j] = 0
            elif j==0:
                delta_Ndxp[i,j] = -(Nd[i,j+1]-Nd[i,j])/delta_Lx
                delta_Ndxn[i,j] = 0
                delta_Ndyp[i,j] = (Nd[i,j]-Nd[i-1,j])/delta_Ly
                delta_Ndyn[i,j] = (Nd[i+1,j]-Nd[i,j])/delta_Ly
            elif j==Lx:
                delta_Ndxp[i,j] = 0
                delta_Ndxn[i,j] = -(Nd[i,j]-Nd[i,j-1])/delta_Lx
                delta_Ndyp[i,j] = (Nd[i,j]-Nd[i-1,j])/delta_Ly
                delta_Ndyn[i,j] = (Nd[i+1,j]-Nd[i,j])/delta_Ly
            #電極境界
            elif phi_Bound[i,j-1] == V and phi_Bound[i+1,j] == V:#Bias側右上辺
                delta_Ndxp[i,j] = -(Nd[i,j+1]-Nd[i,j])/delta_Lx
                delta_Ndxn[i,j] = 0
                delta_Ndyp[i,j] = (Nd[i,j]-Nd[i-1,j])/delta_Ly
                delta_Ndyn[i,j] = 0
            elif phi_Bound[i,j-1] == V and phi_Bound[i-1,j] == V:#Bias側右下辺
                delta_Ndxp[i,j] = -(Nd[i,j+1]-Nd[i,j])/delta_Lx
                delta_Ndxn[i,j] = 0
                delta_Ndyp[i,j] = 0
                delta_Ndyn[i,j] = (Nd[i+1,j]-Nd[i,j])/delta_Ly
            elif phi_Bound[i,j+1] == V and phi_Bound[i+1,j] == V:#Bias側左上辺
                delta_Ndxp[i,j] = 0
                delta_Ndxn[i,j] = -(Nd[i,j]-Nd[i,j-1])/delta_Lx
                delta_Ndyp[i,j] = (Nd[i,j]-Nd[i-1,j])/delta_Ly
                delta_Ndyn[i,j] = 0
            elif phi_Bound[i,j+1] == V and phi_Bound[i-1,j] == V:#Bias側左下辺
                delta_Ndxp[i,j] = 0
                delta_Ndxn[i,j] = -(Nd[i,j]-Nd[i,j-1])/delta_Lx
                delta_Ndyp[i,j] = 0
                delta_Ndyn[i,j] = (Nd[i+1,j]-Nd[i,j])/delta_Ly
            elif phi_Bound[i,j+1] == 0 and phi_Bound[i+1,j] == 0:#GND側の左上辺
                delta_Ndxp[i,j] = 0
                delta_Ndxn[i,j] = -(Nd[i,j]-Nd[i,j-1])/delta_Lx
                delta_Ndyp[i,j] = (Nd[i,j]-Nd[i-1,j])/delta_Ly
                delta_Ndyn[i,j] = 0
            elif phi_Bound[i,j+1] == 0 and phi_Bound[i-1,j] == 0:#GND側の左下辺
                delta_Ndxp[i,j] = 0
                delta_Ndxn[i,j] = -(Nd[i,j]-Nd[i,j-1])/delta_Lx
                delta_Ndyp[i,j] = 0
                delta_Ndyn[i,j] = (Nd[i+1,j]-Nd[i,j])/delta_Ly
            elif phi_Bound[i,j-1] == 0 and phi_Bound[i+1,j] == 0:#GND側の右上辺
                delta_Ndxp[i,j] = -(Nd[i,j+1]-Nd[i,j])/delta_Lx
                delta_Ndxn[i,j] = 0
                delta_Ndyp[i,j] = (Nd[i,j]-Nd[i-1,j])/delta_Ly
                delta_Ndyn[i,j] = 0
            elif phi_Bound[i,j-1] == 0 and phi_Bound[i-1,j] == 0:#GND側右下辺
                delta_Ndxp[i,j] = -(Nd[i,j+1]-Nd[i,j])/delta_Lx
                delta_Ndxn[i,j] = 0
                delta_Ndyp[i,j] = 0
                delta_Ndyn[i,j] = (Nd[i+1,j]-Nd[i,j])/delta_Ly
            elif phi_Bound[i,j+1] == V or phi_Bound[i,j+1] == 0:#Bias,GNDの左Tip
                delta_Ndxp[i,j] = 0
                delta_Ndxn[i,j] = -(Nd[i,j]-Nd[i,j-1])/delta_Lx
                delta_Ndyp[i,j] = (Nd[i,j]-Nd[i-1,j])/delta_Ly
                delta_Ndyn[i,j] = (Nd[i+1,j]-Nd[i,j])/delta_Ly
            elif phi_Bound[i,j-1] == V or phi_Bound[i,j-1] == 0:#Bias,GNDの右Tip
                delta_Ndxp[i,j] = -(Nd[i,j+1]-Nd[i,j])/delta_Lx
                delta_Ndxn[i,j] = 0
                delta_Ndyp[i,j] = (Nd[i,j]-Nd[i-1,j])/delta_Ly
                delta_Ndyn[i,j] = (Nd[i+1,j]-Nd[i,j])/delta_Ly
            elif phi_Bound[i+1,j] == V or phi_Bound[i+1,j] == 0:#Bias,GNDの上Tip
                delta_Ndxp[i,j] = -(Nd[i,j+1]-Nd[i,j])/delta_Lx
                delta_Ndxn[i,j] = -(Nd[i,j]-Nd[i,j-1])/delta_Lx
                delta_Ndyp[i,j] = (Nd[i,j]-Nd[i-1,j])/delta_Ly
                delta_Ndyn[i,j] = 0
            elif phi_Bound[i-1,j] == V or phi_Bound[i-1,j] == 0:#Bias,GNDの下Tip
                delta_Ndxp[i,j] = -(Nd[i,j+1]-Nd[i,j])/delta_Lx
                delta_Ndxn[i,j] = -(Nd[i,j]-Nd[i,j-1])/delta_Lx
                delta_Ndyp[i,j] = 0
                delta_Ndyn[i,j] = (Nd[i+1,j]-Nd[i,j])/delta_Ly
            else:
                delta_Ndxp[i,j] = -(Nd[i,j+1]-Nd[i,j])/delta_Lx
                delta_Ndxn[i,j] = -(Nd[i,j]-Nd[i,j-1])/delta_Lx
                delta_Ndyp[i,j] = (Nd[i,j]-Nd[i-1,j])/delta_Ly
                delta_Ndyn[i,j] = (Nd[i+1,j]-Nd[i,j])/delta_Ly

            Drift_x[i,j] = u * Nd[i,j] * Ex[i,j]
            Drift_y[i,j] = u * Nd[i,j] * Ey[i,j]
            Diff_xp[i,j] = D * delta_Ndxp[i,j]/2
            Diff_xn[i,j] = D * delta_Ndxn[i,j]/2
            Diff_yp[i,j] = D * delta_Ndyp[i,j]/2
            Diff_yn[i,j] = D * delta_Ndyn[i,j]/2

    return Drift_x, Drift_y, Diff_xn, Diff_xp, Diff_yn, Diff_yp, delta_Ndxn, delta_Ndxp, delta_Ndyn, delta_Ndyp

while n_iter <= end_iter:
    V_sweep(0, 6)
    two_electrode_v()

    phi = phi_culc(phi)
    if n_iter%1==0:
        np.savetxt("phi_"+np.str(n_iter)+".txt", phi, delimiter=" ", fmt="%1.8f")
    Ex, Ey = field_culc(Ex,Ey)
    if n_iter%1==0:
        np.savetxt("Ex_"+np.str(n_iter)+".txt", Ex, delimiter=" ", fmt="%1.8f")
        np.savetxt("Ey_"+np.str(n_iter)+".txt", Ey, delimiter=" ", fmt="%1.8f")

    Nd = drift_culc(Nd)
    Diff_xn, Diff_xp, Diff_yn, Diff_yp, Nd = diff_culc(Diff_xn, Diff_xp, Diff_yn, Diff_yp, Nd)
    Drift_x, Drift_y, Diff_xn, Diff_xp, Diff_yn, Diff_yp, delta_Ndxn, delta_Ndxp, delta_Ndyn, delta_Ndyp = \
        next_drift_diff_culc(Drift_x, Drift_y, Diff_xn, Diff_xp, Diff_yn, Diff_yp, delta_Ndxn, delta_Ndxp, delta_Ndyn, delta_Ndyp)

    print("Sum of Nd", Nd.sum())
    n_iter += 1
    print("iteration:", n_iter)

    if n_iter%1==0:
        im2=plt.imshow(Nd,cmap='rainbow', vmax=0.5*10**7, vmin=10**2)
        anim2.append([im2])
        name = "iter="+np.str(n_iter)+".png"
        plt.savefig(name, format = 'png')
        name2 = "Nd_iter="+np.str(n_iter)+".txt"
        np.savetxt(name2, Nd, delimiter=" ", fmt="%1.8f")


pp2=plt.colorbar(orientation="vertical") # カラーバーの表示
pp2.set_label("Carrier density", fontname="Ariel", fontsize=20)
plt.xlabel('X',fontname="Ariel", fontsize=20)
plt.ylabel('Y',fontname="Ariel",  fontsize=20)
plt.title('Drift and Diffusion of Carrier Density')
ani2 = animation.ArtistAnimation(fig2, anim2, interval=200, blit=True,repeat_delay=1000)
ani2.save("Nd.mp4", writer='ffmpeg')#, dpi=300追加で高画質 ffmpeg
