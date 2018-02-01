import numpy as np
import matplotlib.pyplot as plt
import os

os.chdir("/Users/zenyanagata/Desktop/test")

Lx = 384 - 1#分解能をあげるときは、384*○ - 1の形にする
Ly = 384 - 1#必ずLx = Lyで。

end_iter = 10
n_iter = 0
turn = 1 #V_sweep用の変数
v_sweep_list = np.array([])
n_iter_list = np.array([])
def V_sweep(start_voltage, end_voltage):#sweep_rate = V/iter
    global V, turn, v_sweep_list, n_iter, n_iter_list
    sweep_rate = (end_voltage - start_voltage)/(end_iter/2)
    if n_iter == 0:
        V = start_voltage - sweep_rate*turn
    if V == end_voltage:
        turn = -1
    V += sweep_rate*turn
    v_sweep_list = np.append(v_sweep_list, V)
    n_iter_list = np.append(n_iter_list, n_iter)

R = np.array([])
n_iter_list = np.array([])
R_mesh = np.zeros([Lx+1,Ly+1])
phi_Bound = np.genfromtxt("phi_Bound.txt", delimiter=" ")
Ex = np.genfromtxt("Ex.txt", delimiter=" ")
Ey = np.genfromtxt("Ey.txt", delimiter=" ")

u = 1.5*10**4*10**8#電子移動度 μm^2/Vs
e = 1.60217662*10**(-19) #C

root = phi_Bound.copy()
theta = np.arctan(Ey/Ex)*180/np.pi
n_iter = 0
flag = False
while n_iter <= 1000:
    name = "Nd_iter=" + str(n_iter) + ".txt"
    Nd = np.genfromtxt(name, delimiter=" ")
    R_line = np.array([])
    R_ = np.array([])

    for i in range(Ly+1):
        for j in range(Lx+1):
            if phi_Bound[i,j] == 0 or phi_Bound[i,j] == V:
                pass
            elif Nd[i,j] <= 0:
                print("warning: R_mesh[",i,j,"], <= 0(ohm)", "n_iter=", n_iter)
                flag = True
                break
            else:
                R_mesh[i,j] = 1/(2*e*u*Nd[i,j]) #抵抗/mesh
        if flag:
            break
    if flag:
        break


    for i in range(Ly+1):
        for j in range(Lx+1):
            if phi_Bound[i,j] == 0 or phi_Bound[i,j] == V:
                pass
            elif i==95 and j==109:#電極1の下tip
                R_line = np.append(R_line, R_mesh[i,j])
                k,l = i,j
                root[k,l] = 10
                edge = False
                while edge == False:
                    k,l = k+1,l
                    R_line = np.append(R_line, R_mesh[k,l])
                    root[k,l] = 10
                    if phi_Bound[k,l-1] == 0 and phi_Bound[k+1,l] == 0:#GND側右上辺
                        edge = True
                    elif phi_Bound[k,l-1] == 0 and phi_Bound[k-1,l] == 0:#GND側右下辺
                        edge = True
                    elif phi_Bound[k,l+1] == 0 and phi_Bound[k+1,l] == 0:#GND側左上辺
                        edge = True
                    elif phi_Bound[k,l+1] == 0 and phi_Bound[k-1,l] == 0:#GND側左下辺
                        edge = True
                    elif phi_Bound[k-1,l+1] == 0 or phi_Bound[k+1,l+1] == 0 or phi_Bound[k+1,l-1] == 0 or phi_Bound[k-1,l-1] == 0:#斜め
                        edge = True

                R_ = np.append(R_, R_line.sum())
                R_line = np.array([])
                np.savetxt("R_"+str(n_iter)+".txt", R_, delimiter=" ", fmt="%1.8f")

    R_ = [1/r for r in R_]
    R = np.append(R, 1/(np.array(R_).sum()))
    n_iter_list = np.append(n_iter_list, n_iter)
    n_iter +=10

R_data = np.zeros([len(R),2])
for i in range(len(R)):
    R_data[i] = [n_iter_list[i], R[i]]

np.savetxt("R.txt", R_data, delimiter=" ", fmt="%1.8f")
np.savetxt("root.txt", root, delimiter=" ", fmt="%1.8f")
#for plot
plt.imshow(root,cmap='rainbow')
pp=plt.colorbar(orientation="vertical") # カラーバーの表示
pp.set_label("Current root", fontname="Ariel", fontsize=20)
plt.xlabel('X',fontname="Ariel", fontsize=20)
plt.ylabel('Y',fontname="Ariel",  fontsize=20)
plt.title('Root of Current')
plt.savefig("R_root.png", format = 'png', dpi=300)
#plt.show()

plt.show()
plt.figure()
plt.plot(n_iter_list, R)
plt.savefig("抵抗値遷移.png", format = 'png')
plt.show()
