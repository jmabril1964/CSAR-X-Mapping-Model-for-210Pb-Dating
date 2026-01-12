
"""
CSAR_map.py , a X-mapping model for 210Pb-dating assumig a constant SAR

==> @ Created by J.M. Abril-Hernández, 2025.
This code reads a file with experimental 210Pb_exc data in sediment core layers,
a file with random numbers from library (\aleat_S1>) and the configuration.json file.

Outputs:
- `Map3D.txt`: contains the four input values used by each solver and the corresponding X_df.
- `Absolute_min.txt`: contains the absolute minimum found across all solvers.
These outputs can be used by CSAR_cronos.py.

Read the  bibliografic support for understanding fundamentals 
and the Readme.txt file to customize the configuration.json file.
"""

import numpy as np
import random
import sys
import json

# Read the configuration file JSON
with open("configuration.json", "r") as archivo:
    datos = json.load(archivo)

# Assign the values to variables
Core_data = datos["Core_data"]
A0_central = datos["A0_central"]
w_central = datos["w_central"]
sA_central = datos["sA_central"]  
OP_LN_A0 = datos["OP_LN_A0"]  
f_A = datos["f_A"]
f_w = datos["f_w"]
f_sA = datos["f_sA"]
NR = datos["NR"]
kr = int(datos["kr"])
Tmr = datos["Tmr"]
sgt = datos["sgt"]
peso = datos["peso"]



# *** UPLOAD INPUT FILES (core data and file from library)

m_i = []       # mass depth at the bottom of layer i
Aexp_i = []    # 210Pb_exc (Bq/kg) in layer i
sgAexp_i = []  # error in 210Pb_exc (Bq/kg) in layer i

with open(Core_data, 'r') as file:
    for line in file:
        if line.strip():
            parts = line.strip().split()
            m_i.append(float(parts[0]))
            Aexp_i.append(float(parts[1]))
            sgAexp_i.append(float(parts[2]))

# Read the randomly sorted canonical representative sample of size N
# from a standardized Normal distribution (see folder aleat/).
N = len(m_i)
print(f"Number of sediment slices : {N} \n")
file2 = f"./aleat_S1/aleat_S1_{N}.txt"
print(file2)

z_1 = []
with open(file2, 'r') as file_aleat:
    for line in file_aleat:
        if line.strip():
            z_1.append(float(line.strip()))
z_1 = np.array(z_1)

ldPb = 0.03118  # 1/yr. 210Pb decay constant (http://www.lnhb.fr/home/nuclear-data/nuclear-data-table/)



DoF=(N+peso) - 3 # Degees of freedom for reduced X_df

# Main part: for this program we only need to return the chi-value.

def sorting(A0_solver, w_solver, sA_solver):
    """
    Given (A0_solver, w_solver, sA_solver), generates N pairs of (A0i, wi)
    values and find their best ordering so that minimizes the distance to the experimental profile, 
    including radioactive decay and an optional time-mark constraint.

    Returns
    -------
    chi_df : float
        Square root of the reduced chi: (chi/DoF)**0.5 for the best ordering found.
    """
    Sol_T = []    # Solution chronology (cumulative age by slice)
    Sol_A = []    # Ordered initial activities
    Sol_w = []    # Ordered SARs
    Sol_Ath = []  # Theoretical profile (solution including decay)

    slices = N
    chif = 0.0
    tant = 0.0
    
    # --- A0 ---
    if OP_LN_A0:
        # Parameters already in log-space; sA_solver is relative in log-space
        A0i = A0_solver * (1.0 + sA_solver * z_1)
        A0i = np.exp(A0i)
    else:
        # Standard (physical) space with relative deviation
        A0i = A0_solver * (1.0 + sA_solver * z_1)

# --- w ---
    wi = np.full(N, w_solver)
    # A0i  = np.clip(A0i, 0.1 * A0_central, None) # guard against non-physical small values. [Activate this line when needed]
    for k in range(N):
        chiant = 1.0E20
        jmin = 0

        if k == 0:
            dm2 = m_i[0]
        else:
            dm2 = m_i[k] - m_i[k - 1]

        for s in range(slices):
            Dt = dm2 / wi[s]
            Act = A0i[s] * np.exp(-ldPb * tant) * (1 - np.exp(-ldPb * Dt)) / (ldPb * Dt)
            chi = (Act - Aexp_i[k])**2
            if chi < chiant:
                jmin = s
                chiant = chi

        Sol_A.append(A0i[jmin])
        Sol_w.append(wi[jmin])
        Sol_T.append(tant + dm2 / wi[jmin])
        tant = tant + dm2 / wi[jmin]

        A0i = np.delete(A0i, jmin)
        wi  = np.delete(wi,  jmin)
        slices -= 1

        if k == 0:
            Dt = Sol_T[0]
            tanterior = 0
        else:
            Dt = Sol_T[k] - Sol_T[k - 1]
            tanterior = Sol_T[k - 1]

        Act = Sol_A[k] * np.exp(-ldPb * tanterior) * (1 - np.exp(-ldPb * Dt)) / (ldPb * Dt)
        chif += ((Act - Aexp_i[k]) / sgAexp_i[k])**2
        Sol_Ath.append(Act)
    # End k-loop

    # Objective function definition (optional time-mark)
    chif = chif + peso * ((Tmr - Sol_T[kr]) / sgt)**2  # see more in publications
    chi_df = (chif / DoF)**0.5
    return chi_df



# Sweep resolutions
d_A  = f_A  * A0_central * 2 / NR
d_w  = f_w  * w_central  * 2 / NR
d_as = f_sA * sA_central * 2 / NR

Chimin = 2.0E20

with open('Map3D.txt', 'w', encoding='utf-8') as file:
    for ia in range(NR):
        A0_solver = A0_central * (1 - f_A) + ia * d_A
        print(f"Counter ia from 0 to NR = {ia}")
        for iw in range(NR):
            w_solver = w_central * (1 - f_w) + iw * d_w
            for isa in range(NR):
                sA_solver = sA_central * (1 - f_sA) + isa * d_as
                valor = sorting(A0_solver, w_solver, sA_solver)
                if valor < Chimin:
                    Chimin  = valor
                    min_A0  = A0_solver
                    min_w   = w_solver
                    min_sA  = sA_solver
                file.write(f"{A0_solver:.2f}\t{w_solver:.4f}\t{sA_solver:.4f}\t{valor:.3f}\n")
        file.write(f"\n")

print(
    f"Minimum A= {min_A0:.2f}, w = {min_w:.4f}, sa={min_sA:.4f} "
    f"chi_df ={Chimin:.3f}, slices= {N}\n"
)

with open('Absolute_min.txt', 'w', encoding='utf-8') as file4:
    file4.write(f"{min_A0:.2f}\t{min_w:.4f}\t{min_sA:.4f}\t{Chimin:.3f}\t{DoF}\n")
    file4.write(f"{d_A/2:.2f}\t{d_w/2:.5f}\t{d_as/2:.5f}\n")
    file4.write(f"{A0_central:.2f}\t{w_central:.5f}\t{sA_central:.5f}\n")
    file4.write(f"{NR}\t{f_A:.3f}\t{f_w:.3f}\t{f_sA:.3f}\n")
    file4.write(f"{kr}\t{Tmr:.1f}\t{sgt:.2f}\t{peso:.3f}\n")
    file4.write(f"# Line 1: min_A0, min_w, min_sA, X_df, DoF\n")
    file4.write(f"# Line 2: half-resolution intervals for A0, w, sA \n")
    file4.write(f"# Line 3: Input parameters A0_central, w_central, sA_central \n")
    file4.write(f"# Line 4: NR, f_A, f_w, f_sA \n")
    file4.write(f"# Line 5: Time-mark information: kr, Tmr, sgt, peso \n")

# end task
print("Bye!!")
sys.exit()

