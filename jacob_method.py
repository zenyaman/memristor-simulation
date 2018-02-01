import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import ArtistAnimation # アニメーション作成のためのメソッドをインポート
import matplotlib.animation as animation

delta_Lx=0.01
delta_Ly=0.01
LLx = 1 # 正方形の幅
LLy= 1
Lx = int(LLx/delta_Lx)
Ly = int(LLy/delta_Ly)

V = 5.0 # 電圧
convegence_criterion = 10**-4
phi_Bound= np.ones([Lx+1,Ly+1])
phi = np.zeros([Lx+1,Ly+1])
phi_in = np.zeros([Lx+1,Ly+1])
Ex = np.zeros([Lx,Ly])
Ey = np.zeros([Lx,Ly])

# 電極の形状
def electrode_reigion():
    global phi_Bound

    for i in range(15, 51):
        for ii in range(i-14):
            phi_Bound[i, ii] = V #[15, 0] [16, 1] [17, 2] ... [50, 35]
    for i in range(51, 86):
        for ii in range(86-i):
            phi_Bound[i, ii] = V #[51, 34] [52, 33] [53, 32] ... [85, 0]

    for i in range(15, 51):
        for ii in range(115-i, 101):
            phi_Bound[i, ii] = 0 #[15, 100] [16, 99] [17, 98] ... [50, 65]
    for i in range(51, 86):
        for ii in range(i+15, 101):
            phi_Bound[i, ii] = 0 #[51, 66] [52, 67] [53, 68] ... [85, 100]
electrode_reigion()

#静電ポテンシャルの計算
def phi_culc():
    global phi, phi_in

    #main culculation
    delta = 1.0
    n_iter=0
    conv_check=[]
    while delta > convegence_criterion:
        if n_iter % 100 ==0:  # 収束状況のモニタリング
            print("iteration No =", n_iter, "delta=",delta)
        conv_check.append([n_iter, delta])
        for i in range(Lx+1):
            for j in range(Ly+1):
                if phi_Bound[i, j] == 0 or phi_Bound[i, j] == V:
                    phi[i,j] = phi_Bound[i,j]
                elif i == 0:
                    if j == 0:              #角の場合、隣り合うグリッドは2点
                        phi[i,j] = (phi_in[i+1,j] + phi_in[i,j+1])/2
                    elif j == Ly:
                        phi[i,j] = (phi_in[i+1,j] + phi_in[i,j-1])/2
                    else:                   #辺の場合、隣り合うグリッドは3点
                        phi[i,j] = (phi_in[i+1,j]+phi_in[i,j-1]+phi_in[i,j+1])/3
                elif i == Lx:
                    if j==0:                #角の場合、隣り合うグリッドは2点
                        phi[i,j] = (phi_in[i-1,j] + phi_in[i,j+1])/2
                    elif j == Ly:
                        phi[i,j] = (phi_in[i-1,j] + phi_in[i,j-1])/2
                    else:                   #辺の場合、隣り合うグリッドは3点
                        phi[i,j] = (phi_in[i-1,j]+phi_in[i,j-1]+phi_in[i,j+1])/3
                elif j == 0:
                    phi[i,j] = (phi_in[i+1,j]+phi_in[i-1,j]+phi_in[i,j+1])/3
                elif j == Lx:
                    phi[i,j] = (phi_in[i-1,j]+phi_in[i+1,j]+phi_in[i,j-1])/3
                else:
                    phi[i,j] = (phi_in[i+1,j] + phi_in[i-1,j] + phi_in[i,j+1] + phi_in[i,j-1])/4
        delta = np.max(abs(phi-phi_in))
        phi_in=phi.copy()
        n_iter+=1

    #for plot
    plt.imshow(phi,cmap='rainbow')
    pp=plt.colorbar(orientation="vertical") # カラーバーの表示
    pp.set_label("phi", fontname="Ariel", fontsize=20)
    plt.xlabel('X',fontname="Ariel", fontsize=20)
    plt.ylabel('Y',fontname="Ariel",  fontsize=20)
    plt.title('Solution of the Poisson equation for electrostatic potential')
    plt.savefig("Potential.png", format = "png")

    plt.show()
phi_culc()
