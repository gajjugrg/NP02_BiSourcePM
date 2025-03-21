import numpy as np
import math

# --------------------------
# Parameter Initialization
# --------------------------
drlong   = 205.0         # PM drift length in mm
drshort  = 45.0
rainn    = 15.0          # inner anode radius
raout    = 30.0          # outer anode radius
resINNL  = 0.10         # resolution for inner events (long drift)
resOUTL  = 0.10          # resolution for outer events (long drift)
resINNS  = 0.10          # resolution for inner events (short drift)
resOUTS  = 0.12          # resolution for outer events (short drift)
dtave    = 1000.0 / (18.0 + 1.4)  # average time between events in μs
srcfrac  = 1.4 / (18.0 + 1.4)
drfiv    = 1.0 / 1.623     # For E = μs/mm
attdist  = 15000.0 / drfiv  # attenuation distance in mm (e lifetime in μs * v_d)
thr      = 0.26

# A-branch parameters (electrons and gamma)
Aprob1   = 0.0015
Aprob2   = 0.0044 + Aprob1
Aprob3   = 0.0152 + Aprob2
EleEne1  = 0.482
EleEne2  = 0.556
EleEne3  = 0.566
GamEne1  = 0.570
GamInt1  = 80.0

# B-branch parameters (electrons and gamma)
Bprob1   = 0.0054
Bprob2   = 0.0184 + Bprob1
Bprob3   = 0.0703 + Bprob2
Bprob4   = 0.7458 + Bprob3
Bprob5   = 0.0687 + Bprob4
Bprob6   = 0.0002 + Bprob5
Bprob7   = 0.0013 + Bprob6
EleEne4  = 1.060
EleEne5  = 1.049
EleEne6  = 0.976
GamEne2  = 1.064
GamInt2  = 120.0
GamEne3  = 1.77
GamInt3  = 160.0
EleEne7  = 1.682
GamEne4  = 1.440
GamInt4  = 140.0

pi2 = 2.0 * 3.1415927
timtot = 0.0

# --------------------------
# Simulation Loop & Output
# --------------------------
with open("bi207stream.txt", "w") as outfile:
    nEvents = 1000000
    #nEvents = 1000000
    for i in range(1, nEvents + 1):
        # Choose drift length and set flag (ilosh)
        if np.random.random() > srcfrac:
            drmax = drlong
            ilosh = 1
        else:
            drmax = drshort
            ilosh = 0

        # Time increment (exponential distribution)
        dtim = -dtave * math.log(np.random.random())
        timtot += dtim

        # --------------------
        # A-Branch Processing
        # --------------------
        Aprob = np.random.random()
        if Aprob <= Aprob1:
            DInterA = 0.0
            EnergyA = EleEne1
            ieleA = 1
        elif Aprob <= Aprob2:
            DInterA = 0.0
            EnergyA = EleEne2
            ieleA = 1
        elif Aprob <= Aprob3:
            DInterA = 0.0
            EnergyA = EleEne3
            ieleA = 1
        else:
            DInterA = GamInt1
            EnergyA = GamEne1
            ieleA = 0

        # Generate a random starting position inside a circle (radius 2.5)
        while True:
            x0 = 5.0 * np.random.random() - 2.5
            y0 = 5.0 * np.random.random() - 2.5
            if x0 * x0 + y0 * y0 < 6.25:
                break
        z0 = 0.0

        # Propagation for A-branch
        dA      = -DInterA * math.log(np.random.random())
        cthetaA = np.random.random()
        phiA    = np.random.random() * pi2
        sthetaA = math.sqrt(1.0 - cthetaA * cthetaA)
        x1A     = x0 + dA * sthetaA * math.sin(phiA)
        y1A     = y0 + dA * sthetaA * math.cos(phiA)
        z1A     = z0 + dA * cthetaA
        r1A     = math.sqrt(x1A * x1A + y1A * y1A)

        IokA = 0
        ElectEA = 0.0
        rrA = 0.0
        if EnergyA > 0 and (0.0 <= z1A < drmax and r1A < raout):
            if ieleA == 0:  # Gamma scattering simulation
                scattered = False
                gg = EnergyA / 0.511
                while not scattered:
                    ct = 2.0 * np.random.random() - 1.0
                    eps = 1.0 / (1.0 + gg * (1.0 - ct))
                    sctprob = 0.5 * (eps ** 2) * (eps + 1.0/eps - (1.0 - ct * ct))
                    if np.random.random() <= sctprob:
                        scattered = True
                        GammaE = EnergyA * eps
                        ElectEA = EnergyA - GammaE
            else:
                ElectEA = EnergyA
            # Attenuate energy based on propagation distance
            ElectEA *= math.exp((z1A - drmax) / attdist)
            # Choose resolution based on radial position
            if r1A <= rainn:
                IokA = 2
                resol = resINNL if ilosh == 1 else resINNS
            else:
                IokA = 1
                resol = resOUTL if ilosh == 1 else resOUTS
            # Add resolution noise (sum of 12 random numbers, shifted by 6)
            res_noise = sum(np.random.random() for _ in range(12))
            res_noise = resol * (res_noise - 6.0)
            rrA = ElectEA + res_noise

        # --------------------
        # B-Branch Processing
        # --------------------
        Bprob = np.random.random()
        if Bprob <= Bprob1:
            DInterB = 0.0
            EnergyB = EleEne4
            ieleB = 1
        elif Bprob <= Bprob2:
            DInterB = 0.0
            EnergyB = EleEne5
            ieleB = 1
        elif Bprob <= Bprob3:
            DInterB = 0.0
            EnergyB = EleEne6
            ieleB = 1
        elif Bprob <= Bprob4:
            DInterB = GamInt2
            EnergyB = GamEne2
            ieleB = 0
        elif Bprob <= Bprob5:
            DInterB = GamInt3
            EnergyB = GamEne3
            ieleB = 0
        elif Bprob <= Bprob6:
            DInterB = 0.0
            EnergyB = EleEne7
            ieleB = 1
        elif Bprob <= Bprob7:
            DInterB = GamInt4
            EnergyB = GamEne4
            ieleB = 0
        else:
            DInterB = 0.0
            EnergyB = 0.0
            ieleB = 0

        # Generate a starting position for B-branch
        while True:
            x0 = 5.0 * np.random.random() - 2.5
            y0 = 5.0 * np.random.random() - 2.5
            if x0 * x0 + y0 * y0 < 6.25:
                break
        z0 = 0.0

        dB      = -DInterB * math.log(np.random.random())
        cthetaB = np.random.random()
        phiB    = np.random.random() * pi2
        sthetaB = math.sqrt(1.0 - cthetaB * cthetaB)
        x1B     = x0 + dB * sthetaB * math.sin(phiB)
        y1B     = y0 + dB * sthetaB * math.cos(phiB)
        z1B     = z0 + dB * cthetaB
        r1B     = math.sqrt(x1B * x1B + y1B * y1B)

        IokB = 0
        ElectEB = 0.0
        rrB = 0.0
        if EnergyB > 0 and (0.0 <= z1B < drmax and r1B < raout):
            if ieleB == 0:  # Gamma scattering simulation
                scattered = False
                gg = EnergyB / 0.511
                while not scattered:
                    ct = 2.0 * np.random.random() - 1.0
                    eps = 1.0 / (1.0 + gg * (1.0 - ct))
                    sctprob = 0.5 * (eps ** 2) * (eps + 1.0/eps - (1.0 - ct * ct))
                    if np.random.random() <= sctprob:
                        scattered = True
                        GammaE = EnergyB * eps
                        ElectEB = EnergyB - GammaE
            else:
                ElectEB = EnergyB
            ElectEB *= math.exp((z1B - drmax) / attdist)
            if r1B <= rainn:
                IokB = 2
                resol = resINNL if ilosh == 1 else resINNS
            else:
                IokB = 1
                resol = resOUTL if ilosh == 1 else resOUTS
            res_noise = sum(np.random.random() for _ in range(12))
            res_noise = resol * (res_noise - 6.0)
            rrB = ElectEB + res_noise

        # --------------------
        # Final Checks & Output
        # --------------------
        if rrA < thr:
            IokA = 0
        if rrB < thr:
            IokB = 0

        dtimA = timtot + (drmax - z1A) * drfiv
        dtimB = timtot + (drmax - z1B) * drfiv

        # Write events in order depending on the z-coordinate
        if IokA != 0 and IokB != 0:
            if z1B > z1A:
                outfile.write(f"{i} {ilosh} {IokB} {dtimB} {ieleB} {r1B} {z1B} {EnergyB} {ElectEB} {rrB}\n")
                outfile.write(f"{i} {ilosh} {IokA} {dtimA} {ieleA} {r1A} {z1A} {EnergyA} {ElectEA} {rrA}\n")
            else:
                outfile.write(f"{i} {ilosh} {IokA} {dtimA} {ieleA} {r1A} {z1A} {EnergyA} {ElectEA} {rrA}\n")
                outfile.write(f"{i} {ilosh} {IokB} {dtimB} {ieleB} {r1B} {z1B} {EnergyB} {ElectEB} {rrB}\n")
        elif IokA != 0:
            outfile.write(f"{i} {ilosh} {IokA} {dtimA} {ieleA} {r1A} {z1A} {EnergyA} {ElectEA} {rrA}\n")
        elif IokB != 0:
            outfile.write(f"{i} {ilosh} {IokB} {dtimB} {ieleB} {r1B} {z1B} {EnergyB} {ElectEB} {rrB}\n")

    print("Simulation finished. Total time:", timtot)
