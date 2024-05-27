import numpy as np
import matplotlib.pyplot as plt
import sys

filenames = ["b1", "b2", "b3", "b4"]
directs = ["wyniki\\b-a-1000-5-4_o50\\", "wyniki\\b-a-1000-5-4_o100\\", "wyniki\\b-a-1000-5-4_o150\\", "wyniki\\e-r-1000-008_o50\\", "wyniki\\e-r-1000-008_o100\\", "wyniki\\e-r-1000-008_o150\\"]
k_width = 50
label_dict = {
    "b1": "beta_vect = [1.0, 1.0, 1.0, 1.0, 1.0]",
    "b2": "beta_vect = [0.5, 0.5, 0.5, 0.5, 0.5]",
    "b3": "beta_vect = [0.2, 0.2, 0.2, 0.2, 0.2]",
    "b4": "beta_vect = [1.0, 0.8, 0.6, 0.4, 0.2]",
}

#plt.figure(figsize=(45,70))
plt.figure(figsize=(1.2 * 10, 1.4 * 10), dpi=400)
#plt.title("Precyzja", fontsize = 30)
#plt.xlabel("gamma", fontsize = 25)
#plt.ylabel("prec", fontsize = 25)

fig, axs = plt.subplots(3, 2)
for i, dir in enumerate(directs):
    for name in filenames:
        file = open(dir + name + ".txt")
        
        gamma = []
        prec_kor = []
        rank_kor = []
        prec_err = []
        rank_err = []
        sizes = []

        for line in file:
            gamma.append(float(line.split(" ")[0]))
            prec_kor.append(float(line.split(" ")[1]))
            rank_kor.append(float(line.split(" ")[2]))
            if i % 3 == 2:
                prec_err.append(float(line.split(" ")[3]) * 0.3 * k_width)
                rank_err.append(float(line.split(" ")[4]) * 0.3 * k_width)
            else:
                prec_err.append(float(line.split(" ")[3]) * 0.3 * k_width)
                rank_err.append(float(line.split(" ")[4]) * 0.3 * k_width)   
            sizes.append(5)
        file.close()
    #print(str(i % 3) + " " + str(i  // 3))
        axs[i % 3, i // 3].scatter(gamma, prec_kor, sizes, label = label_dict[name], linestyle="None")
        axs[i % 3, i // 3].errorbar(gamma, prec_kor, prec_err, ls='none')
        #axs[i % 3, i // 3].xticks(fontsize = 20)
        #axs[i % 3, i // 3].yticks(fontsize = 20)

        #axs[i % 3, i // 3].legend(fontsize = 25)

for ax in axs.flat:
    ax.set(xlabel="gamma", ylabel="precyzja")
    ax.set_xticks(np.arange(0, 4, step=0.5))
    ax.set_yticks(np.arange(0, 1, step=0.2))
    ax.grid()
for ax in axs.flat:
    ax.label_outer()
axs[0, 1].legend(loc='lower left', fontsize="5", bbox_to_anchor=(0.55, 0.7))
axs[0, 1].set_title(label="E-R")
axs[0, 0].set_title(label="B-A")
text1 = '10%'
text2 = '15%'
text3 = '20%'
axs[0, 1].text(4.33, 0.5, text1)
axs[1, 1].text(4.33, 0.5, text2)
axs[2, 1].text(4.33, 0.5, text3)         
plt.savefig(f"precyzja_od_gammy_{sys.argv[1]}.png")