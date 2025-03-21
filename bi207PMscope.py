import numpy as np
import math

# --------------------------
# Initialization
# --------------------------

# Allocate arrays
trace1 = np.zeros(4000)
trace2 = np.zeros(4000)
trace3 = np.zeros(4000)
trace4 = np.zeros(4000)

spcOUTL = np.zeros(400, dtype=int)
spcEINL = np.zeros(400, dtype=int)
spcGINL = np.zeros(400, dtype=int)
spcOUTS = np.zeros(400, dtype=int)
spcEINS = np.zeros(400, dtype=int)
spcGINS = np.zeros(400, dtype=int)

# Compute signal array (500 points)
signal = np.zeros(500)
for k in range(500):  
    x = (k + 1) * 0.1
    signal[k] = math.exp(-(x - 10.) / 5.) / (0.56988 * (1 + math.exp(-(x - 10.) / 1.25)))

# Open files
fin   = open("bi207stream.txt", "r")
fOUTL = open("outer_anode_long.txt", "w")
fEINL = open("inner_anode_long.txt", "w")
fOUTS = open("outer_anode_short.txt", "w")
fEINS = open("inner_anode_short.txt", "w")
fSPEC = open("bi207spectra.txt", "w")

# --------------------------
# Initialize variables
# --------------------------
thr = 0.28
ioffset = 1500
timtrig = 0.0
phmaxOUTL = 0.0
phmaxINNL = 0.0
phmaxOUTS = 0.0
phmaxINNS = 0.0
ieleOUTL = -1
ieleINNL = -1
ieleOUTS = -1
ieleINNS = -1

# --------------------------
# Main Processing Loop
# --------------------------
# We'll loop until we run out of lines or i exceeds 999999.
# Here, each input line should contain: i, ilosh, iok, dtim, iele, r1, z1, energy, electe, rr
while True:
    line = fin.readline()
    if not line:
        break  # end-of-file
    try:
        parts = line.strip().split()
        # Convert values appropriately
        i_val = int(parts[0])
        ilosh = int(parts[1])
        iok = int(parts[2])
        dtim = float(parts[3])
        iele = int(parts[4])
        r1 = float(parts[5])
        z1 = float(parts[6])
        energy = float(parts[7])
        electe = float(parts[8])
        rr = float(parts[9])
    except Exception as e:
        print("Error parsing line:", line, e)
        continue

    # Compute starting index for traces (kstart)
    kstart = int(10 * (dtim - timtrig) + ioffset)

    # If within time window (dtim - timtrig <= 200)
    if (dtim - timtrig) <= 200.:
        if ilosh == 1:  # long drift
            if iok == 1:
                if rr > phmaxOUTL:
                    phmaxOUTL = rr
                    ieleOUTL = iele
                kk = kstart
                for k in range(500):
                    kk += 1
                    if kk < len(trace1):
                        trace1[kk] += signal[k] * rr
            elif iok == 2:
                if rr > phmaxINNL:
                    phmaxINNL = rr
                    ieleINNL = iele
                kk = kstart
                for k in range(500):
                    kk += 1
                    if kk < len(trace2):
                        trace2[kk] += signal[k] * rr
        elif ilosh == 0:  # short drift
            if iok == 1:
                if rr > phmaxOUTS:
                    phmaxOUTS = rr
                    ieleOUTS = iele
                kk = kstart
                for k in range(500):
                    kk += 1
                    if kk < len(trace3):
                        trace3[kk] += signal[k] * rr
            elif iok == 2:
                if rr > phmaxINNS:
                    phmaxINNS = rr
                    ieleINNS = iele
                kk = kstart
                for k in range(500):
                    kk += 1
                    if kk < len(trace4):
                        trace4[kk] += signal[k] * rr

    # If dtim - timtrig is negative or zero, update trigger parameters
    if (dtim - timtrig) <= 0.:
        ioffset = kstart
        timtrig = dtim

    # If time window exceeds 200 units, process and reset traces
    if (dtim - timtrig) > 200.:
        # Find maximum pulse values in the traces over [ioffset, ioffset+2200]
        phOUTL = np.max(trace1[ioffset:ioffset+2200])
        phINNL = np.max(trace2[ioffset:ioffset+2200])
        phOUTS = np.max(trace3[ioffset:ioffset+2200])
        phINNS = np.max(trace4[ioffset:ioffset+2200])

        # Optionally print if ioffset is out of expected range
        if ioffset > 1800 or ioffset <= 0:
            print(">>>>>>>>", ioffset, "<<<<<<<<")

        # Write out events if maximum pulse exceeds threshold
        if phmaxOUTL > thr:
            fOUTL.write(f"{i_val} {ieleOUTL} {timtrig:.3f} {phmaxOUTL:.3f} {phOUTL:.3f}\n")
            k_spec = int(phOUTL * 200)
            if k_spec < 400:
                spcOUTL[k_spec] += 1

        if phmaxINNL > thr:
            fEINL.write(f"{i_val} {ieleINNL} {timtrig:.3f} {phmaxINNL:.3f} {phINNL:.3f}\n")
            k_spec = int(phINNL * 200)
            if k_spec < 400:
                if ieleINNL == 1:
                    spcEINL[k_spec] += 1
                else:
                    spcGINL[k_spec] += 1

        if phmaxOUTS > thr:
            fOUTS.write(f"{i_val} {ieleOUTS} {timtrig:.3f} {phmaxOUTS:.3f} {phOUTS:.3f}\n")
            k_spec = int(phOUTS * 200)
            if k_spec < 400:
                spcOUTS[k_spec] += 1

        if phmaxINNS > thr:
            fEINS.write(f"{i_val} {ieleINNS} {timtrig:.3f} {phmaxINNS:.3f} {phINNS:.3f}\n")
            k_spec = int(phINNS * 200)
            if k_spec < 400:
                if ieleINNS == 1:
                    spcEINS[k_spec] += 1
                else:
                    spcGINS[k_spec] += 1

        # Reset trigger and trace parameters
        timtrig = dtim
        ioffset = 1500
        kstart = ioffset
        phmaxOUTL = phmaxINNL = phmaxOUTS = phmaxINNS = 0.0
        trace1.fill(0.0)
        trace2.fill(0.0)
        trace3.fill(0.0)
        trace4.fill(0.0)
        ieleOUTL = ieleINNL = ieleOUTS = ieleINNS = -1

        # Immediately accumulate the current event into traces after reset
        kk = ioffset
        if ilosh == 1:
            if iok == 1:
                phmaxOUTL = rr
                ieleOUTL = iele
                for k in range(500):
                    kk += 1
                    if kk < len(trace1):
                        trace1[kk] += signal[k] * rr
            elif iok == 2:
                phmaxINNL = rr
                ieleINNL = iele
                for k in range(500):
                    kk += 1
                    if kk < len(trace2):
                        trace2[kk] += signal[k] * rr
        elif ilosh == 0:
            if iok == 1:
                phmaxOUTS = rr
                ieleOUTS = iele
                for k in range(500):
                    kk += 1
                    if kk < len(trace3):
                        trace3[kk] += signal[k] * rr
            elif iok == 2:
                phmaxINNS = rr
                ieleINNS = iele
                for k in range(500):
                    kk += 1
                    if kk < len(trace4):
                        trace4[kk] += signal[k] * rr

# End of main loop

# --------------------------
# Write Spectra Data
# --------------------------
for k in range(400):
    val = k * 0.005 - 0.0025
    sum_inner = spcEINL[k] + spcGINL[k]
    sum_short = spcEINS[k] + spcGINS[k]
    fSPEC.write(f"{val:.5f} {spcOUTL[k]} {spcEINL[k]} {spcGINL[k]} {sum_inner} "
                f"{spcOUTS[k]} {spcEINS[k]} {spcGINS[k]} {sum_short}\n")

# Close all files
fin.close()
fOUTL.close()
fEINL.close()
fOUTS.close()
fEINS.close()
fSPEC.close()
print("Program Ended")
