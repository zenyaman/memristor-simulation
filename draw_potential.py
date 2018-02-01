import numpy as np
import matplotlib.pyplot as plt
import os

os.chdir("/Users/zenyanagata/Desktop/STO_simulation/modified_sim2")




n_iter = 0
while n_iter <= 2000:
    phi = np.genfromtxt("phi_"+str(n_iter)+".txt", delimiter=" ")

    plt.imshow(phi,cmap='rainbow')
    pp=plt.colorbar(orientation="vertical") # カラーバーの表示
    pp.set_label("phi", fontname="Ariel", fontsize=20)
    plt.xlabel('X',fontname="Ariel", fontsize=20)
    plt.ylabel('Y',fontname="Ariel",  fontsize=20)
    plt.title('Solution of the Poisson equation for electrostatic potential ')
    plt.savefig("Potential_"+str(n_iter)+".png", format = 'png', dpi=300)

    n_iter += 100
