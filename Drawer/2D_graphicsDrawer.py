import matplotlib.pyplot as plt
import numpy as np
import sys

#A_data_path = sys.argv[1]
#E_data_path = sys.argv[2]
#output_path = sys.argv[3]

A_data_path = "D:\\CodeRepos\\Diplom\\Data\\Output\\ToDraw\\2_dim\\Receivers\\A.txt"
E_data_path = "D:\\CodeRepos\\Diplom\\Data\\Output\\ToDraw\\2_dim\\Receivers\\E.txt"
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
time_iter = []

text_e = open(E_data_path, 'r').read().split("\n")

i = 0
for text in text_e:
    if text == '':
        break
    string = text.split(" ")
    time_iter.append(i)
    i+=1
    time.append(float(string[0]))
    e0.append(abs(float(string[1])))
    e1.append(abs(float(string[2])))
    e2.append(abs(float(string[3])))
    e3.append(abs(float(string[4])))

time = time[2:]
e0 = e0[2:]
e1 = e1[2:]
e2 = e2[2:]
e3 = e3[2:]

plt.loglog(time, e0, color='green', label="E(0, y, -130, 1.0)")
plt.loglog(time, e1, color='blue', label="(2500.0, 0.0, -200.0)")
plt.loglog(time, e2, color='red', label="(10.0, 0.0, -700.0)")
plt.loglog(time, e3, color='pink', label="(1000.0, 0.0, -1250.0)")
#plt.ylim(min(e0), -1 * min(e0))
plt.xlabel("T [c]")
plt.ylabel("E [В / м]")
plt.legend()
plt.savefig(output_path + 'Log_E.png')
plt.show()
plt.clf()