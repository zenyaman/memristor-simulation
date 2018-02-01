import numpy as np
import matplotlib.pyplot as plt
import os

load_Nd_file ="/Users/sakailab/Desktop/nagata/1.27/4-elec_4V_iter=100"
active_file ="/Users/sakailab/Desktop/nagata/1.27/4-elec_4V_iter=100_R_1V"

resolution = 2/3.5
delta_Lx=1/resolution# 1 meshgridの長さ
delta_Ly=1/resolution
Lx = int(384*resolution - 1)#分解能をあげるときは、384*○ - 1の形にする
Ly = int(384*resolution - 1)#必ずLx = Lyで。

V = 1 # 電圧
end_iter = 50
n_iter = 0
n_iter_list = np.array([])
R = np.array([])
R_mesh = np.zeros([Lx+1,Ly+1])
phi_Bound = np.ones([Lx+1,Ly+1])*100
phi = np.ones([Lx+1,Ly+1])
Ex = np.zeros([Lx+1, Ly+1])
Ey =  np.zeros([Lx+1, Ly+1])
u = 0.4*10**8#シリコン中の電子移動度 μm^2/Vs
e = 1.60217662*10**(-19) #C
eps = 85*10**(-6) #F/μm

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

def phi_culc():
    global phi
    convegence_criterion = 10**-5
    # for SOR method
    aa_recta=0.5*(np.cos(np.pi/Lx)+np.cos(np.pi/Ly))
    omega_SOR_recta = 2/(1+np.sqrt(1-aa_recta**2)) # SOR法における加速パラメータ
    print("omega_SOR_rect=", omega_SOR_recta)

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
        delta = np.max(abs(phi-phi_in))
        n_iter+=1

    #np.savetxt("phi.txt", phi, delimiter=" ", fmt="%1.8f")
    #plt.imshow(phi,cmap='rainbow')
    #pp=plt.colorbar(orientation="vertical") # カラーバーの表示
    #pp.set_label("phi", fontname="Ariel", fontsize=20)
    #plt.xlabel('X',fontname="Ariel", fontsize=20)
    #plt.ylabel('Y',fontname="Ariel",  fontsize=20)
    #plt.title('Solution of the Poisson equation for electrostatic potential ')
    #plt.savefig("Potential.png", format = 'png', dpi=300)

def field_culc():
    global Ex, Ey

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

    #np.savetxt("Ex.txt", Ex, delimiter=" ", fmt="%1.8f")
    #np.savetxt("Ey.txt", Ey, delimiter=" ", fmt="%1.8f")
    #plt.figure()
    #X, Y = np.meshgrid(np.arange(0, Lx+1, 1), np.arange(Ly+1, 0, -1))  # メッシュ生成
    #plt.quiver(X,Y,Ex,Ey,color='red',angles='xy',scale_units='xy', scale=10) # ベクトル場をプロット
    #plt.xlim([0,Lx+1]) # 描くXの範囲
    #plt.ylim([0,Ly+1]) # 描くyの範囲
    #plt.grid()
    #plt.xlabel('X',fontname="Ariel", fontsize=20)
    #plt.ylabel('Y',fontname="Ariel",  fontsize=20)
    #plt.title('Solution of the Poisson equation for electrostatic field')
    #plt.savefig("電場分布.png", format = 'png', dpi=300)


flag = False
while n_iter <= end_iter:
    root = phi_Bound.copy()
    os.chdir(load_Nd_file)
    Nd = np.genfromtxt("Nd_iter=" + str(n_iter) + ".txt", delimiter=" ")

    os.chdir(active_file)
    phi_culc()
    if n_iter%1==0:
        np.savetxt("phi_"+str(n_iter)+".txt", phi, delimiter=" ", fmt="%1.8f")
    field_culc()
    if n_iter%1==0:
        np.savetxt("Ex_"+str(n_iter)+".txt", Ex, delimiter=" ", fmt="%1.8f")
        np.savetxt("Ey_"+str(n_iter)+".txt", Ey, delimiter=" ", fmt="%1.8f")

    theta = np.arctan(Ey/Ex)*180/np.pi
    R_line = np.array([])
    R_ = np.array([])

    for i in range(Ly+1):
        for j in range(Lx+1):
            if phi_Bound[i,j] == 0 or phi_Bound[i,j] == V:
                pass
            elif Nd[i,j] < 0:
                print("warning: Nd[",i,j,"] = ",Nd[i,j], "n_iter=", n_iter)
                phi_Bound[i,j] = 20
                plt.figure()
                plt.imshow(phi_Bound, cmap="rainbow")
                plt.show()
                flag = True
                break
            else:
                R_mesh[i,j] = 1/(2*e*u*(Nd[i,j]+10**2)) #抵抗/mesh
        if flag:
            break
    if flag:
        break

    for i in range(Ly):
        for j in range(Lx):
            if phi_Bound[i,j] == 0 or phi_Bound[i,j] == V:
                pass
            elif phi_Bound[i,j+1] == 0 and phi_Bound[i+1,j] == 0:#GND側左上辺
                R_line = np.append(R_line, R_mesh[i,j])
                theta_line = abs(theta[i,j])
                k,l,m,n = i,j,i,j
                root[k,l] = 10
                edge = False
                while edge == False:
                    if Ex[m,n]<=0 and Ey[m,n]<=0:
                        if 0<=theta_line<11.25:#1
                            k,l,m,n = k,l+1,m,n+2
                        elif 11.25<=theta_line<22.5:#2
                            k,l,m,n = k,l+1,m-1,n+2
                        elif 22.5<=theta_line<33.75:#3
                            k,l,m,n = k-1,l+1,m-1,n+2
                        elif 33.75<=theta_line<45:#4
                            k,l,m,n = k-1,l+1,m-2,n+2
                        elif 45<=theta_line<56.25:#5
                            k,l,m,n = k-1,l+1,m-2,n+2
                        elif 56.25<=theta_line<67.5:#6
                            k,l,m,n = k-1,l+1,m-2,n+1
                        elif 67.5<=theta_line<78.75:#7
                            k,l,m,n = k-1,l,m-2,n+1
                        elif 78.75<=theta_line<=90:#8
                            k,l,m,n = k-1,l,m-2,n
                    elif Ex[m,n]>=0 and Ey[m,n]>=0:
                        if 0<=theta_line<11.25:#1
                            k,l,m,n = k,l-1,m,n-2
                        elif 11.25<=theta_line<22.5:#2
                            k,l,m,n = k,l-1,m+1,n-2
                        elif 22.5<=theta_line<33.75:#3
                            k,l,m,n = k+1,l-1,m+1,n-2
                        elif 33.75<=theta_line<45:#4
                            k,l,m,n = k+1,l-1,m+2,n-2
                        elif 45<=theta_line<56.25:#5
                            k,l,m,n = k+1,l-1,m+2,n-2
                        elif 56.25<=theta_line<67.5:#6
                            k,l,m,n = k+1,l-1,m+2,n-1
                        elif 67.5<=theta_line<78.75:#7
                            k,l,m,n = k+1,l,m+2,n-1
                        elif 78.75<=theta_line<=90:#8
                            k,l,m,n = k+1,l,m+2,n
                    elif Ex[m,n]<=0 and Ey[m,n]>=0:
                        if 0<=theta_line<11.25:#1
                            k,l,m,n = k,l+1,m,n+2
                        elif 11.25<=theta_line<22.5:#2
                            k,l,m,n = k,l+1,m+1,n+2
                        elif 22.5<=theta_line<33.75:#3
                            k,l,m,n = k+1,l+1,m+1,n+2
                        elif 33.75<=theta_line<45:#4
                            k,l,m,n = k+1,l+1,m+2,n+2
                        elif 45<=theta_line<56.25:#5
                            k,l,m,n = k+1,l+1,m+2,n+2
                        elif 56.25<=theta_line<67.5:#6
                            k,l,m,n = k+1,l+1,m+2,n+1
                        elif 67.5<=theta_line<78.75:#7
                            k,l,m,n = k+1,l,m+2,n+1
                        elif 78.75<=theta_line<=90:#8
                            k,l,m,n = k+1,l,m+2,n
                    elif Ex[m,n]>=0 and Ey[m,n]<=0:
                        if 0<=theta_line<11.25:#1
                            k,l,m,n = k,l-1,m,n-2
                        elif 11.25<=theta_line<22.5:#2
                            k,l,m,n = k,l-1,m-1,n-2
                        elif 22.5<=theta_line<33.75:#3
                            k,l,m,n = k-1,l-1,m-1,n-2
                        elif 33.75<=theta_line<45:#4
                            k,l,m,n = k-1,l-1,m-2,n-2
                        elif 45<=theta_line<56.25:#5
                            k,l,m,n = k-1,l-1,m-2,n-2
                        elif 56.25<=theta_line<67.5:#6
                            k,l,m,n = k-1,l-1,m-2,n-1
                        elif 67.5<=theta_line<78.75:#7
                            k,l,m,n = k-1,l,m-2,n-1
                        elif 78.75<=theta_line<=90:#8
                            k,l,m,n = k-1,l,m-2,n
                    R_line = np.append(R_line, R_mesh[k,l])
                    root[k,l] = 10
                    if phi_Bound[k,l-1] == V and phi_Bound[k+1,l] == V:#Bias側右上辺
                        edge = True
                    elif phi_Bound[k,l-1] == V and phi_Bound[k-1,l] == V:#Bias側右下辺
                        edge = True
                    elif phi_Bound[k,l+1] == V and phi_Bound[k+1,l] == V:#Bias側左上辺
                        edge = True
                    elif phi_Bound[k,l+1] == V and phi_Bound[k-1,l] == V:#Bias側左下辺
                        edge = True
                    elif phi_Bound[k-1,l+1] == V or phi_Bound[k+1,l+1] == V or phi_Bound[k+1,l-1] == V or phi_Bound[k-1,l-1] == V:#斜め
                        edge = True
                    if not edge:
                        R_line = np.append(R_line, R_mesh[m,n])
                        root[m,n] = 10

                    if phi_Bound[m,n-1] == V and phi_Bound[m+1,n] == V:#Bias側右上辺
                        edge = True
                    elif phi_Bound[m,n-1] == V and phi_Bound[m-1,n] == V:#Bias側右下辺
                        edge = True
                    elif phi_Bound[m,n+1] == V and phi_Bound[m+1,n] == V:#Bias側左上辺
                        edge = True
                    elif phi_Bound[m,n+1] == V and phi_Bound[m-1,n] == V:#Bias側左下辺
                        edge = True
                    elif phi_Bound[m-1,n+1] == V or phi_Bound[m+1,n+1] == V or phi_Bound[m+1,n-1] == V or phi_Bound[m-1,n-1] == V:#斜め
                        edge = True

                    theta_line = abs(theta[m,n])
                    k,l = m,n
                R_ = np.append(R_, R_line.sum())
                R_line = np.array([])


            elif phi_Bound[i,j-1] == 0 and phi_Bound[i+1,j] == 0:#GND側右上辺
                R_line = np.append(R_line, R_mesh[i,j])
                theta_line = abs(theta[i,j])
                k,l,m,n = i,j,i,j
                root[k,l] = 10
                edge = False
                while edge == False:
                    if Ex[m,n]<=0 and Ey[m,n]<=0:
                        if 0<=theta_line<11.25:#1
                            k,l,m,n = k,l+1,m,n+2
                        elif 11.25<=theta_line<22.5:#2
                            k,l,m,n = k,l+1,m-1,n+2
                        elif 22.5<=theta_line<33.75:#3
                            k,l,m,n = k-1,l+1,m-1,n+2
                        elif 33.75<=theta_line<45:#4
                            k,l,m,n = k-1,l+1,m-2,n+2
                        elif 45<=theta_line<56.25:#5
                            k,l,m,n = k-1,l+1,m-2,n+2
                        elif 56.25<=theta_line<67.5:#6
                            k,l,m,n = k-1,l+1,m-2,n+1
                        elif 67.5<=theta_line<78.75:#7
                            k,l,m,n = k-1,l,m-2,n+1
                        elif 78.75<=theta_line<=90:#8
                            k,l,m,n = k-1,l,m-2,n
                    elif Ex[m,n]>=0 and Ey[m,n]>=0:
                        if 0<=theta_line<11.25:#1
                            k,l,m,n = k,l-1,m,n-2
                        elif 11.25<=theta_line<22.5:#2
                            k,l,m,n = k,l-1,m+1,n-2
                        elif 22.5<=theta_line<33.75:#3
                            k,l,m,n = k+1,l-1,m+1,n-2
                        elif 33.75<=theta_line<45:#4
                            k,l,m,n = k+1,l-1,m+2,n-2
                        elif 45<=theta_line<56.25:#5
                            k,l,m,n = k+1,l-1,m+2,n-2
                        elif 56.25<=theta_line<67.5:#6
                            k,l,m,n = k+1,l-1,m+2,n-1
                        elif 67.5<=theta_line<78.75:#7
                            k,l,m,n = k+1,l,m+2,n-1
                        elif 78.75<=theta_line<=90:#8
                            k,l,m,n = k+1,l,m+2,n
                    elif Ex[m,n]<=0 and Ey[m,n]>=0:
                        if 0<=theta_line<11.25:#1
                            k,l,m,n = k,l+1,m,n+2
                        elif 11.25<=theta_line<22.5:#2
                            k,l,m,n = k,l+1,m+1,n+2
                        elif 22.5<=theta_line<33.75:#3
                            k,l,m,n = k+1,l+1,m+1,n+2
                        elif 33.75<=theta_line<45:#4
                            k,l,m,n = k+1,l+1,m+2,n+2
                        elif 45<=theta_line<56.25:#5
                            k,l,m,n = k+1,l+1,m+2,n+2
                        elif 56.25<=theta_line<67.5:#6
                            k,l,m,n = k+1,l+1,m+2,n+1
                        elif 67.5<=theta_line<78.75:#7
                            k,l,m,n = k+1,l,m+2,n+1
                        elif 78.75<=theta_line<=90:#8
                            k,l,m,n = k+1,l,m+2,n
                    elif Ex[m,n]>=0 and Ey[m,n]<=0:
                        if 0<=theta_line<11.25:#1
                            k,l,m,n = k,l-1,m,n-2
                        elif 11.25<=theta_line<22.5:#2
                            k,l,m,n = k,l-1,m-1,n-2
                        elif 22.5<=theta_line<33.75:#3
                            k,l,m,n = k-1,l-1,m-1,n-2
                        elif 33.75<=theta_line<45:#4
                            k,l,m,n = k-1,l-1,m-2,n-2
                        elif 45<=theta_line<56.25:#5
                            k,l,m,n = k-1,l-1,m-2,n-2
                        elif 56.25<=theta_line<67.5:#6
                            k,l,m,n = k-1,l-1,m-2,n-1
                        elif 67.5<=theta_line<78.75:#7
                            k,l,m,n = k-1,l,m-2,n-1
                        elif 78.75<=theta_line<=90:#8
                            k,l,m,n = k-1,l,m-2,n
                    R_line = np.append(R_line, R_mesh[k,l])
                    root[k,l] = 10
                    if phi_Bound[k,l-1] == V and phi_Bound[k+1,l] == V:#Bias側右上辺
                        edge = True
                    elif phi_Bound[k,l-1] == V and phi_Bound[k-1,l] == V:#Bias側右下辺
                        edge = True
                    elif phi_Bound[k,l+1] == V and phi_Bound[k+1,l] == V:#Bias側左上辺
                        edge = True
                    elif phi_Bound[k,l+1] == V and phi_Bound[k-1,l] == V:#Bias側左下辺
                        edge = True
                    elif phi_Bound[k-1,l+1] == V or phi_Bound[k+1,l+1] == V or phi_Bound[k+1,l-1] == V or phi_Bound[k-1,l-1] == V:#斜め
                        edge = True
                    if not edge:
                        R_line = np.append(R_line, R_mesh[m,n])
                        root[m,n] = 10

                    if phi_Bound[m,n-1] == V and phi_Bound[m+1,n] == V:#Bias側右上辺
                        edge = True
                    elif phi_Bound[m,n-1] == V and phi_Bound[m-1,n] == V:#Bias側右下辺
                        edge = True
                    elif phi_Bound[m,n+1] == V and phi_Bound[m+1,n] == V:#Bias側左上辺
                        edge = True
                    elif phi_Bound[m,n+1] == V and phi_Bound[m-1,n] == V:#Bias側左下辺
                        edge = True
                    elif phi_Bound[m-1,n+1] == V or phi_Bound[m+1,n+1] == V or phi_Bound[m+1,n-1] == V or phi_Bound[m-1,n-1] == V:#斜め
                        edge = True

                    theta_line = abs(theta[m,n])
                    k,l = m,n
                R_ = np.append(R_, R_line.sum())
                R_line = np.array([])

            elif phi_Bound[i+1,j] == 0:#GND側上辺
                R_line = np.append(R_line, R_mesh[i,j])
                theta_line = abs(theta[i,j])
                k,l = i,j
                root[k,l] = 10
                edge = False
                while edge == False:
                    k = k-1
                    R_line = np.append(R_line, R_mesh[k,l])
                    root[k,l] = 10
                    if phi_Bound[k-1,l] == V:
                        edge = True
                R_ = np.append(R_, R_line.sum())
                R_line = np.array([])
                #np.savetxt("R_"+str(n_iter)+".txt", R_, delimiter=" ", fmt="%1.8f")

    R_ = [1/r for r in R_]
    R = np.append(R, 1/(np.array(R_).sum()))
    n_iter_list = np.append(n_iter_list, n_iter)

    #for plot
    plt.figure()
    plt.imshow(root,cmap='rainbow')
    pp=plt.colorbar(orientation="vertical") # カラーバーの表示
    pp.set_label("Current root", fontname="Ariel", fontsize=20)
    plt.xlabel('X',fontname="Ariel", fontsize=20)
    plt.ylabel('Y',fontname="Ariel",  fontsize=20)
    plt.title('Root of Current')
    plt.savefig("R_root"+str(n_iter)+".png", format = 'png', dpi=300)

    n_iter +=1
    print("n_iter:",n_iter)

R_data = np.zeros([len(R),2])
for i in range(len(R)):
    R_data[i] = [n_iter_list[i], R[i]]

np.savetxt("R.txt", R_data, delimiter=" ", fmt="%1.8f")

plt.figure()
plt.plot(n_iter_list, R)
plt.savefig("抵抗値遷移.png", format = 'png')
