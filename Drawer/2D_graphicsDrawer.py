import matplotlib.pyplot as plt
import numpy as np
import sys

#A_data_path = sys.argv[1]
#E_data_path = sys.argv[2]
#output_path = sys.argv[3]

A_data_path = "D:\\CodeRepos\\Diplom\\Data\\Output\\ToDraw\\3_dim\\Receivers\\A.txt"
E_data_path = "D:\\CodeRepos\\Diplom\\Data\\Output\\ToDraw\\3_dim\\Receivers\\E.txt"
output_path = "D:\\CodeRepos\\Diplom\\Drawer\\Graphics\\"


#A_data_path = "D:\\CodeRepos\\Diplom\\Data\\Output\\ToDraw\\3_dim\\Graphics\\A_With_anomaly_1.05.txt"
#E_data_path = "D:\\CodeRepos\\Diplom\\Data\\Output\\ToDraw\\3_dim\\Graphics\\E_With_anomaly_1.05.txt"
#output_path = "D:\\CodeRepos\\Diplom\\Drawer\\Graphics\\"


a0 = []
a1 = []
a2 = []
a3 = []
e0 = []
e1 = []
e2 = []
e3 = []
time = []

text_a = open(A_data_path, 'r').read().split("\n")

for text in text_a:
    if text == '':
        break
    string = text.split(" ")
    time.append(float(string[0]))
    a0.append(float(string[1]))
    a1.append(float(string[2]))
    a2.append(float(string[3]))
    a3.append(float(string[4]))

#time = time[50:200]
#a0 = a0[50:200]
#a1 = a1[50:200]
#a2 = a2[50:200]
#a3 = a3[50:200]


plt.plot(time, a0, color='green', label="A(0, y, -130, 1.0)")
plt.plot(time, a1, color='blue', label="(2500.0, 0.0, -200.0)")
plt.plot(time, a2, color='red', label="(10.0, 0.0, -700.0)")
plt.plot(time, a3, color='pink', label="(1000.0, 0.0, -1250.0)")
#plt.ylim(-1 * max(a0) - 5 * max(a0), max(a0) + 5 * max(a0))
plt.xlabel("Y [м]")
plt.ylabel("A [Тл * м]")
plt.legend()
plt.savefig(output_path + 'Normal_A_obj2.png')
plt.clf()

plt.semilogy(time, a0, color='green', label="A(0, y, -130, 1.0)")
plt.semilogy(time, a1, color='blue', label="(2500.0, 0.0, -200.0)")
plt.semilogy(time, a2, color='red', label="(10.0, 0.0, -700.0)")
plt.semilogy(time, a3, color='pink', label="(1000.0, 0.0, -1250.0)")
#plt.ylim(min(a0), -1 * min(a0))
plt.xlabel("Y [м]")
plt.ylabel("A [Тл * м]")
plt.legend()
plt.savefig(output_path + 'Log_A_obj2.png')
plt.clf()

text_e = open(E_data_path, 'r').read().split("\n")

for text in text_e:
    if text == '':
        break
    string = text.split(" ")
    e0.append(abs(float(string[1])))
    e1.append(abs(float(string[2])))
    e2.append(abs(float(string[3])))
    e3.append(abs(float(string[4])))

#time = time[6:200]
#e0 = e0[50:200]
#e1 = e1[50:200]
#e2 = e2[50:200]
#e3 = e3[50:200]

plt.plot(time, e0, color='green', label="E(0, y, -130, 1.0)")
plt.plot(time, e1, color='blue', label="(2500.0, 0.0, -200.0)")
plt.plot(time, e2, color='red', label="(10.0, 0.0, -700.0)")
plt.plot(time, e3, color='pink', label="(1000.0, 0.0, -1250.0)")
#plt.ylim(-1 * max(e0) - 5 * max(e0), max(e0) + 5 * max(e0))
plt.xlabel("Y [м]")
plt.ylabel("E [В / м]")
plt.legend()
plt.savefig(output_path + 'Normal_E_obj2.png')
plt.clf()

plt.semilogy(time, e0, color='green', label="E(0, y, -130, 1.0)")
plt.semilogy(time, e1, color='blue', label="(2500.0, 0.0, -200.0)")
plt.semilogy(time, e2, color='red', label="(10.0, 0.0, -700.0)")
plt.semilogy(time, e3, color='pink', label="(1000.0, 0.0, -1250.0)")
#plt.ylim(min(e0), -1 * min(e0))
plt.xlabel("Y [м]")
plt.ylabel("E [В / м]")
plt.legend()
plt.savefig(output_path + 'Log_E_obj2.png')
plt.clf()
#plt.show()