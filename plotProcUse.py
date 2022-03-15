import numpy as np
import matplotlib.pyplot as plt

nx,ny=1024,3072

tx=[]
txg=[]

tx1,tx2,tx3,tx4,tx5,tx6,tx7,tx8=0.406,0.308,0.227,0.236,0.178,0.271,0.729,1.09
txg8,txg7,txg6,txg5,txg4,txg3,txg2,txg1=1.09,0.545,0.2725,0.13625,0.068125,0.0340625,0.01703125,0.008515625

ty=[]

tx.append(tx1),tx.append(tx2)
tx.append(tx3),tx.append(tx4)
tx.append(tx5),tx.append(tx6)
tx.append(tx7),tx.append(tx8)

ty.append(1024),ty.append(512)
ty.append(256),ty.append(128)
ty.append(64),ty.append(32)
ty.append(16),ty.append(8)

txg.append(txg1),txg.append(txg2)
txg.append(txg3),txg.append(txg4)
txg.append(txg5),txg.append(txg6)
txg.append(txg7),txg.append(txg8)

t1px,t1py=16,64
t2px,t2py=8,64
t3px,t3py=8,32
t4px,t4py=8,16
t5px,t5py=4,16
t6px,t6py=4,8
t7px,t7py=4,4

#print('===========')
#print(np.sqrt((t1px*t1py)))
#print(np.sqrt(t2px*t2py))
#print(np.sqrt(t3px*t3py))
#print(np.sqrt(t4px*t4py))
#print(np.sqrt(t5px*t5py))
#print('===========')
#print(np.sqrt((nx/t1px)*(ny/t1py)))
#print(np.sqrt((nx/t2px)*(ny/t2py)))
#print(np.sqrt((nx/t3px)*(ny/t3py)))
#print(np.sqrt((nx/t4px)*(ny/t4py)))
#print('===========')

#fit=np.polyfit(tx,ty,deg=2)

plt.scatter(ty,tx)
plt.plot(ty,txg,color='orange',ls=':',alpha=0.5)
#plt.plot(ty,tx,color='orange',ls=':',alpha=0.5)
#plt.plot(fit,color='orange')
plt.ylabel(r'$\mu$'+'[s]')
plt.ylim(1e-1)
plt.yscale('log')
plt.xlabel('proc')
plt.xscale('log')
#plt.show()
plt.savefig('pltprocuse.png')
