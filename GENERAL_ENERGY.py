# disk and Orbital Energy
# ======================================

# get the disk and Orbital energy

OE = -0.5*par1/ts.xq2

OE_Sum = []

 i = 0
  di = 1

   while i <= len(Orbit_Len)-2:
        OE_Sum.append(
            np.sum(
                OE[
                    Orbit_Len[i]:Orbit_Len[i+1]
                ]
            )
        )
        i = i+di

    Int = 0
    dInt = 1

    Orbit = Max_Orbits-2
    # Orbit=1

    KE_Sum = []
    UE_Sum = []
    UINT_Sum = []

    while Int <= Orbit:

        ff = pc.read_var(trimall=True, ivar=Int, magic=["TT"])
        TT = ff.TT
        rad = ff.x
        phi = ff.y
        uu = ff.uu
        ux = ff.ux
        uy = ff.uy
        rho = ff.rho
        i = 0
        di = 1
        j = 0
        dj = 1
        KE = []
        UE = []
        UINT = []

    while j <= len(phi)-1:
        while i <= len(rad)-1:
            try:
                KE.append(
                    0.5*np.abs(phi[j])*((rad[i+1]**2)-(rad[i]**2))
                    * (ux[j][i]**2+uy[j][j]**2)
                    * 0.5*rho[i][j]
                )
                UE.append(
                    rho[i][j]
                    / rad[i]
                )
                UINT.append(
                    (rho[i][j]
                     * TT[i][j])
                    / gamma
                )
                i = i+di
            except:

                # out of bound

                i = i+di
        i = 0
        j = j+dj

    KE_Sum.append(np.sum(KE[:]))
    UE_Sum.append(np.sum(UE[:]))
    UINT_Sum.append(np.sum(UINT[:]))
    Int = Int+dInt

    i = 0
    di = 1
    Total = []

    while i <= len(KE_Sum)-1:
        Total.append(KE_Sum[i]+UINT_Sum[i]-UE_Sum[i])
        i = i+di
