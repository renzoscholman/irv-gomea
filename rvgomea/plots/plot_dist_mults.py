import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("statistics.dat", delim_whitespace=True, skiprows=1, header=None)
fig, axs = plt.subplots(1, 3, figsize=(10,6))
# df.plot(x=0, y=6, ax=axs[0])
# df.plot(x=0, y=7, ax=axs[1])
df.plot(x=0, y=9, ax=axs[0])
df.plot(x=0, y=10, ax=axs[1])
df.plot(x=0, y=11, ax=axs[2])
plt.show()	