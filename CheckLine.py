import numpy as np
import matplotlib.pyplot as plt

x = np.arange(0,60,1)

y = np.exp(x)^-0.2

plt.plot(y)
plt.show()

if __name__ == "__main__":
    print("plotting")