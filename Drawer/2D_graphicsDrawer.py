import matplotlib.pyplot as plt
import numpy as np
import sys

#A_data_path = sys.argv[1]
#E_data_path = sys.argv[2]
#output_path = sys.argv[3]

A_data_path = "D:\\CodeRepos\\Diplom\\Data\\Output\\ToDraw\\2_dim\\Receivers\\A.txt"
E_data_path = "D:\\CodeRepos\\Diplom\\Data\\Output\\ToDraw\\2_dim\\Receivers\\E.txt"
output_path = "D:\\CodeRepos\\Diplom\\Drawer\\Graphics\\"
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

time = time[2:100]
a0 = a0[2:100]
a1 = a1[2:100]
a2 = a2[2:100]
a3 = a3[2:100]


#plt.plot(time, a0, color='green', label="(699.997, 0.0, -3.0)")
#plt.plot(time, a1, color='blue', label="(697.737, 0.0, -30.0)")
#plt.plot(time, a2, color='red', label="(673.205, 0.0, -100.0)")
plt.plot(time, a3, color='pink', label="(605.357, 0.0, -170.0)")
plt.xlabel("T [c]")
plt.ylabel("A [Тл * м]")
plt.legend()
plt.savefig(output_path + 'Normal_A.png')
plt.clf()

#plt.semilogy(time, a0, color='green', label="(699.997, 0.0, -3.0)")
#plt.semilogy(time, a1, color='blue', label="(697.737, 0.0, -30.0)")
#plt.semilogy(time, a2, color='red', label="(673.205, 0.0, -100.0)")
plt.semilogy(time, a3, color='pink', label="(605.357, 0.0, -170.0)")
plt.xlabel("T [c]")
plt.ylabel("A [Тл * м]")
plt.legend()
plt.savefig(output_path + 'Log_A.png')
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

e0 = e0[2:100]
e1 = e1[2:100]
e2 = e2[2:100]
e3 = e3[2:100]

#plt.plot(time, e0, color='green', label="(699.997, 0.0, -3.0)")
#plt.plot(time, e1, color='blue', label="(697.737, 0.0, -30.0)")
#plt.plot(time, e2, color='red', label="(673.205, 0.0, -100.0)")
plt.plot(time, e3, color='pink', label="(605.357, 0.0, -170.0)")
plt.xlabel("T [c]")
plt.ylabel("E [В / м]")
plt.legend()
plt.savefig(output_path + 'Normal_E.png')
plt.clf()
#plt.semilogy(time, e0, color='green', label="(699.997, 0.0, -3.0)")
#plt.semilogy(time, e1, color='blue', label="(697.737, 0.0, -30.0)")
#plt.semilogy(time, e2, color='red', label="(673.205, 0.0, -100.0)")
plt.semilogy(time, e3, color='pink', label="(605.357, 0.0, -170.0)")
plt.xlabel("T [c]")
plt.ylabel("E [В / м]")
plt.legend()
plt.savefig(output_path + 'Log_E.png')
plt.clf()
#plt.show()