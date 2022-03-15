import numpy as np
import pandas as pd

nprocx=np.arange(0,35,1)
nprocy=nprocx*3

#nxgrid=60*nprocx
#nygrid=60*nprocy

nxgrid=256
nygrid=nxgrid*3

ncpus=nprocx*nprocy

checkx=nxgrid/ncpus
checky=nygrid/ncpus

data={'nxgrid':nxgrid,
      'nygrid':nygrid,
      'nprocx':nprocx,
      'nprocy':nprocy,
      'checkx':checkx,
      'checky':checky,
      'ncpus':ncpus}

df = pd.DataFrame(data=data)

print('')
print(checkx)
print('')
print(checky)
print('')

print('')
print(df)
print('')
