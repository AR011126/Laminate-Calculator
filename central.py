"""
central.py
------------
Main user interface for composite laminate calculator

Handles:
    - user input
    - laminate definition (four input modes)
    - optional strength and load analysis
    - mid plane strains and ply-by-ply stresses
    - optional first ply failure
    - optional failure index plotting

All numerical mechanics implemented inside mechanics.py


"""

from mechanics import (
    calculate_abd,
    calculate_transformed_plies,
    compute_ply_strains_stresses,
    laminate_effective_properties,
    first_ply_failure
)
import numpy as np
import matplotlib.pyplot as plt

def get_user_laminate():
    """
    Collects laminate definition from the user

    Supports four input modes

    Returns
    ---------
    laminate : list of dicts
        Ply-by-ply definition including angle, thickness, stiffness properties and strength properties

    strength_enabled: bool
        True if every ply has valid strength properties
    
    
    """

# Ask user for laminate input mode
    print("\nChoose laminate input mode:")
    print("1 - Use same material properties and thickness for all plies")
    print("2 - Use same thickness but different material properties for each ply")
    print("3 - Use different thicknesses but the same material properties for each ply")
    print("4 - Use different thicknesses and different material properties for each ply")

    mode = int(input("\n Select mode (1-4):  "))
    if mode not in [1, 2, 3, 4]:
        raise ValueError("Mode must be an integer between 1 and 4.")
    
    laminate = []
    n = int(input("Enter number of plies: "))

#=============
# if properties are shared
#============
    if mode in [1,3]:
        print("\n Enter shared material properties (MPa):")
        E1_shared = float(input("  E1:  ")) * 1e6
        E2_shared = float(input("  E2:  ")) * 1e6
        G12_shared = float(input("  G12:  ")) * 1e6
        nu12_shared = float(input("  nu12:  "))

        print("\nEnter shared strength properties (MPa):")
        Xt_shared = float(input("  Xt:  ")) * 1e6
        Xc_shared = float(input("  Xc:  ")) * 1e6
        Yt_shared = float(input("  Yt:  ")) * 1e6
        Yc_shared = float(input("  Yc:  ")) * 1e6
        S_shared = float(input("  S:  ")) * 1e6

#===============
# if thickness is chared
#==================
    if mode in [1,2]:
        t_shared_mm = float(input("\n Enter shared thickness (mm): "))
        t_shared = t_shared_mm * 1e-3

#=============
# Ply per ply input
#=============
    for i in range(n):
        print(f"\n--- Ply {i+1} ---")
        angle = float(input("  Angle (deg):  "))

    #thickness selection
        if mode in [3,4]:
            t_mm = float(input("  Thickness (mm):  "))
            t = t_mm * 1e-3
        else:
            t = t_shared

    #material selection
        if mode in [2,4]:
            print("  Enter material properties (MPa):")
            E1 = float(input(  "E1:  ")) * 1e6
            E2 = float(input(  "E2:  ")) * 1e6
            G12 = float(input(  "g12:  ")) * 1e6
            nu12 = float(input(  "nu12:  "))

            print(" Enter strength properties (MPa):")
            Xt = float(input(  "Xt:  ")) * 1e6
            Xc = float(input(  "Xc:  ")) * 1e6
            Yt = float(input(  "Yt:  ")) * 1e6
            Yc = float(input(  "Yc:  ")) * 1e6
            S = float(input(  "S:  ")) * 1e6

        else:
            E1 = E1_shared
            E2 = E2_shared
            G12 = G12_shared
            nu12 = nu12_shared
            Xt = Xt_shared
            Xc = Xc_shared
            Yt = Yt_shared
            Yc = Yc_shared
            S  = S_shared
        laminate.append({
            "angle": angle,
            "thickness": t,
            "E1": E1,
            "E2": E2,
            "G12": G12,
            "nu12": nu12,
            "Xt": Xt,
            "Xc": Xc,
            "Yt": Yt,
            "Yc": Yc,
            "S": S
        })

    strength_enabled = all(
        ply["Xt"] > 0 and ply["Xc"] > 0 and ply["Yt"] > 0 and ply["Yc"] > 0 and ply["S"] > 0
        for ply in laminate
    )
    return laminate, strength_enabled

##MAIN PROGRAM
def main():
    """
    Main execution function
    1) Reads laminate definition
    2) Computes ABD matrices
    3) Computes effective laminate properties
    4) Optionally computes stresses and strains under load
    5) Optionally performs first ply failure analysis
    6) Optionally plots failure indices
    
    """
    print("=== Composite Laminate Calculator ===")

    laminate, strength_enabled = get_user_laminate()

    transformed_plies = calculate_transformed_plies(laminate)
    
    for i, ply in enumerate(transformed_plies):
        print(f"\nPly {i+1}: Angle {ply['angle']} deg, Thickness {ply['thickness']*1000:.2f} mm")
        print("Q_bar:")
        for row in ply['Q_bar']:
            print([round(val,2) for val in row])



    A, B, D = calculate_abd(transformed_plies)

    print("\n=== A Matrix (Extensional) ===")
    for row in A:
        print([round(float(v),2) for v in row])

    print("\n=== B Matrix (Coupling) ===")
    for row in B:
        print([round(float(v),2) for v in row])

    print("\n=== D Matrix (Bending) ===")
    for row in D:
        print([round(float(v),2) for v in row])

    #Effective laminate properties
    total_thickness = sum(p['thickness'] for p in transformed_plies)
    Ex, Ey, vxy, Gxy = laminate_effective_properties(A, total_thickness)

    print("\n=== Effective Laminate Properties ===")
    print(f"  Ex  = {Ex/1e9:.3f} GPa")
    print(f"  Ey  = {Ey/1e9:.3f} GPa")
    print(f"  vxy  = {vxy:.3f}")
    print(f"  Gxy  = {Gxy/1e9:.3f} GPa")

    if strength_enabled:
        use_loads = input("\n Apply external loads? (y/n): ").strip().lower()
    else:
        print("\n Strength properties not provided therefore load analysis disabled.")
        use_loads = "n"

    if use_loads == "y":
        #Ask for loads
        print("\nEnter applied loads (Press enter for 0):")
        Nx = float(input("  Nx (N/m):  ") or 0.0)
        Ny = float(input("  Ny (N/m):  ") or 0.0)
        Nxy = float(input("  Nxy (N/m):  ") or 0.0)
        Mx = float(input("  Mx (Nm/m):  ") or 0.0)
        My = float(input("  My (Nm/m):  ") or 0.0)
        Mxy = float(input("  Mxy (Nm/m):  ") or 0.0)

        N_vec = np.array([Nx, Ny, Nxy])
        M_vec = np.array([Mx, My, Mxy])

        #compute strains/curavtures and ply stresses
        eps0, kappa, ply_results = compute_ply_strains_stresses(transformed_plies, A, B, D, N_vec, M_vec)

        print("\nMidplane strains (eps0):", [round(float(x), 6) for x in eps0])
        print("\nCurvatures (kappa):", [round(float(x), 6) for x in kappa])

        #print per ply result
        for p in ply_results:
            print(f"\nPly {p['ply_index']}   |   theta = {p['angle']} degrees, z_mid = {p['z_mid']*1000:.2f} mm")
            print("  eps global:", [round(float(x), 6) for x in p['eps_global']])
            print("  Sigma global  (Pa):", [round(float(x), 2,) for x in p['sigma_global']])
            print("  Sigma Local (Pa):", [round(float(x), 2) for x in p['sigma_local']])

            if strength_enabled:
                print(f"  FI(Max Stress):  {p['FI_max']:.4f}")
                print(f"  FI(Tsai-Hill):   {p['FI_th']:.4f}")
                print(f"  FI(Tsai-Wu):     {p['FI_tw']:.4f}")
                print(f"  FI(Hashin):      {p['FI_hashin']:.4f}")
                print(f"  Dominant failure mode (Hashin): {p['Failure_mode']}")
                print(f"  SF(Max): {p['SF_max']:.2f},  SF(TH): {p['SF_th']:.2f},  SF(TW): {p['SF_tw']:.2f}, SF(Hashin): {p['SF_hashin']:.2f}")

#============FPF ANALYSIS===========

        if strength_enabled:
            print("\n=== First ply Failure (FPF) Analysis ===")

            fpf_result = first_ply_failure(transformed_plies, A, B, D, N_vec, M_vec)

            if fpf_result["load_factor"] is not None:
                print(f"\n--> First ply failure predicted at load factor: {fpf_result['load_factor']:.2f}")
                print(f"    Failed Ply: {fpf_result['ply_index']}")
                print(f"    Dominant Mode: {fpf_result['mode']}")
                print(f"    Failure Indices at failure load:")
                for name, value in fpf_result["FI_values"].items():
                    print(f"      {name}: {value:.3f}")
            else:
                print("\nNo ply failed up to the maximum load factor limit.")
    
            plot_failure_indices(ply_results)

        else:
            print("\n Strength properties not provided therefore skipping FPF Analysis")
    else:
        print("\n Loads disabled - skipping strains, stress and failure calculations")
    



def plot_failure_indices(ply_results):

    plies = [p['ply_index'] for p in ply_results]
    FI_max = [p['FI_max'] for p in ply_results]
    FI_th = [p['FI_th'] for p in ply_results]
    FI_tw = [p['FI_tw'] for p in ply_results]
    FI_hashin = [p['FI_hashin'] for p in ply_results]

    plt.figure(figsize=(8,5))
    width = 0.2
    x = np.arange(len(plies))

    plt.bar(x - 1.5*width, FI_max, width, label='Max Stress')
    plt.bar(x - 0.5*width, FI_th, width, label='Tsai-Hill')
    plt.bar(x + 0.5*width, FI_tw, width, label='Tsai-Wu')
    plt.bar(x + 1.5*width, FI_hashin, width, label='Hashin')

    plt.axhline(1.0, color='r', linestyle='--', label='Failure Limit (FI = 1)')
    plt.xticks(x, [f'Ply {p["ply_index"]}\n({p["angle"]}(Deg))' for p in ply_results])
    plt.ylabel("Failure Index (FI)")
    plt.xlabel("Ply Number / Angle")
    plt.title("Failure Indices per Ply")
    plt.legend(
        loc='center left',
        bbox_to_anchor=(1.02, 0.5),
        frameon = True
    )
    plt.grid(True, linestyle='--', alpha=0.4)
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()

