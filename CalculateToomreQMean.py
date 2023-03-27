T_max = [
    0.007031, # AGNRun2072
    0.007529, # AGNRun2073
    0.007655, # AGNRun2074
    0.008906, # AGNRun2075
    0.006175, # AGNRun2076
    0.009786, # AGNRun2077
    0.007696, # AGNRun2078
    0.01186, # AGNRun2079
    0.009042, # AGNRun2080
    0.007938, # AGNRun2081
          ]

T1=[0.004464,
    0.004464,
    0.004464,
    0.004464,
    0.004464,
    0.004464,
    0.004464,
    0.004464,
    0.004464,
    0.004464]

DeltaT = []

for i in range(len(T1)):
    DeltaT.append((T_max[i] - T1[i])/T1[i])
    print(T_max[i] - T1[i])

meanToomreq = []
meanToomreq = sum(DeltaT)/len(DeltaT)

print("meanToomreq",meanToomreq)