import math
import numpy as np
"""
mechanics.py - Composite Laminate Mechanics Library
---------------------------------------------------

This module implements classical laminate theory and failure analysis methods for composite materials and includes:
- Q matrix and transformed Q bar matrix calculations
- ABD laminate stiffness matrices
- Mid plane strains and curvatures
- Ply-level stresses in global and local axes
- Maximum stress, tsai-hill, tsai-wu and hashin failure criteria
- first ply failure (FPF) load prediction

"""
#==============================================================
                           #Q-MATRIX
#=============================================================
def calculate_q_matrix(E1, E2, G12, nu12):
    """
    Computes the reduced stiffness matrix (Q) for a single lamina
    in its principle material axes (1-2 coordinate system)

    Parameters
    ---------
    E1, E2 : float (Youngs moduli in fibre and transverse directions)
    G12: float (In-plane shear modulus)
    nu12: float (Major Poisson's ratio)


    Returns
    --------
    Q: list[list[float]] (3x3 reduced stiffness matric for plane-stress)
    """
# Compute the reduced stiffness matrix [Q] for a unidirectional lamina
# Inputs: 
# - E1, E2 : float (Young's moduli along fibre (1) and transverse (2) directions.)
# - G12 : float (In-plane shear modulus.)
# - nu12 : float (Major Poisson’s ratio.)
# Returns:
# list[list[float]]
# 3x3 Q matrix   
    nu21 = (E2/E1) * nu12
#denominator definition 
    denom = 1 - nu12 * nu21
#standard orthotropic stiffness matrix
    Q11 = E1 / denom
    Q22 = E2 / denom
    Q12 = (nu12 * E2) / denom
    Q66 = G12
# [Q] matrix is built following classical laminate theory, Q16 and Q26 are zero in local coordinates for orthotropic materials
    return [[Q11, Q12, 0],
            [Q12, Q22, 0],
            [0, 0, Q66]]
#==============================================================
                           #ROTATE Q-MATRIX TO Q BAR
#=============================================================
def rotate_q_matrix(Q, theta_deg):
    """
    Transforms the lamina stiffness matrix Q into the laminate coordinates for angle theta
    Produces the transformed stiffness Q_bar

    Step explenation:
    Q_bar represents how stiffness behaves once the lamina is rotated in the laminate. It is derived from tensor transformation relations.

    Parameters
    ----------
    Q: 3x3 lamina stiffness matrix in principle material axes
    theta: float (ply angle in degrees)

    Returns
    ---------
    Q_bar: 3x3 matrix (transformed lamina stiffness matrix for use in ABD matrices)
    
    """
    #Transforms the lamina stiffness matrix [Q] into global coordinates [Q bar]
    #  rotates q to q bar for ply at angle theta 
    ##returns a 3 by 3 Q bar matrix
    theta = math.radians(theta_deg)
    m = math.cos(theta)
    n = math.sin(theta)
# The individual values within the [Q] matrix are then extracted from the list of lists built before to be used in calculations
    Q11, Q12, Q16 = Q[0]
    Q12, Q22, Q26 = Q[1]  #Q12 and Q16 repeated intentionally (symmetry)
    Q16, Q26, Q66 = Q[2]
#CLASSICAL LAMINATE THEORY stiffness transformation
# these equations relate local (1-2) behaviour to laminate (x-y)
    Q_bar11 = Q11 * m**4 + 2*(Q12 +2*Q66)*m**2*n**2 + Q22 *n**4
    Q_bar12 = (Q11 + Q22 -4*Q66)*m**2*n**2 + Q12*(m**4+ n**4)
    Q_bar22 = Q11 * n**4 + 2*(Q12 + 2*Q66)*m**2*n**2 + Q22 * m**4
    Q_bar16 = (Q11 - Q12 -2*Q66)*m**3*n - (Q22 - Q12 - 2*Q66)*m*n**3
    Q_bar26 = (Q11 - Q12 - 2*Q66)*m*n**3 - (Q22 - Q12 - 2*Q66)*m**3*n
    Q_bar66 = (Q11 + Q22 - 2*Q12 - 2*Q66)*m**2*n**2 + Q66*(m**4 + n**4)
# A list of list (like the [Q]) matrix is built, this is just the transformed [Q] matrix
    return [
        [Q_bar11, Q_bar12, Q_bar16],
        [Q_bar12, Q_bar22, Q_bar26],
        [Q_bar16, Q_bar26, Q_bar66]
    ]

#==============================================================
                           # BUILD TRANSFORMED PLIES
#=============================================================

def calculate_transformed_plies(laminate):
    """
    For each ply in the laminate definition, compute and store:
    Q_bar
    thickness
    fibre angle
    strength properties

    Returns list of structured ply dictionaries
    
    """
#attaches q bar matrices and strength property to each ply
    transformed = []
#this will call up the relevant material properties from the ply dictionary which we asked the user to define
# computes q in local coordinates
    for ply in laminate:
        Q = calculate_q_matrix(
            ply['E1'],
            ply['E2'],
            ply['G12'],
            ply['nu12']
        )
# rotate q to global coordinates
        Q_bar = rotate_q_matrix(Q, ply["angle"])

#store all relevant ply data
        transformed.append({
            "angle": ply["angle"],
            "thickness": ply["thickness"],
            "Q_bar": Q_bar,
            "Xt": ply['Xt'],
            "Xc": ply['Xc'],
            "Yt": ply['Yt'],
            "Yc": ply['Yc'],
            "S": ply['S'],


        })
        
        

    return transformed

#==============================================================
                           #A B D MATRIX CALCULATIONS
#=============================================================

def calculate_abd(transformed_plies):
    """
    Constructs the laminate ABD matrices using classical laminate theory

    Steps:
    1) Build z coordinates (top and bottom of each ply)
    2) Integrate Q_bar across thickness to obtain ABD

    Returns
    --------
    A, B, D: ndarray
    
    """
    total_t = sum(ply["thickness"] for ply in transformed_plies)
    z_bottom = -total_t / 2

    z = [z_bottom]
    for ply in transformed_plies:
        z.append(z[-1]+ply["thickness"]) #cumulative heights
    
    #Initialise ABD Matrices as zero matrices
    A = np.zeros((3,3,))
    B = np.zeros((3,3,))
    D = np.zeros((3,3,))

    #2 numerical integration per ply
    for k, ply in enumerate(transformed_plies):
        Qb = np.array(ply["Q_bar"])
        z_k = z[k+1]
        z_km1 = z[k]

        A += Qb * (z_k - z_km1)
        B += 0.5 * Qb * (z_k**2 - z_km1**2)
        D += (1/3) * Qb * (z_k**3 - z_km1**3)

    return A, B, D


#==============================================================
#                    EFFECTIVE LAMINATE PROPERTIES
#==============================================================
def laminate_effective_properties(A, total_thickness):
    """
    Computes homogenised laminate-level engineering constants:
    Ex, Ey, vxy, Gxy

    These are derived using the inverse of the normalised A-matrix

    Returns
    --------
    Ex, Ey, vxy, Gxy : float
    
    """

    a_bar = np.linalg.inv(A / total_thickness)
    Ex = 1 / (a_bar[0,0])
    Ey = 1 / (a_bar[1,1])
    vxy = -a_bar[0,1] * Ex
    Gxy = 1 / (a_bar[2,2])

    return Ex, Ey, vxy, Gxy



#==============================================================
                           #SOLVE STRAIN +_ CURVATURE
#=============================================================

def solve_midplane_strains (A,B,D, N_vec, M_vec):
    """
    Soves the 6x6 ABD system to obtain:
    - mid plane strains (eps0)
    - curvatures (kappa)

    Returns
    ---------
    eps0, kappa
    
    """

    top = np.hstack((A, B))
    bot = np.hstack((B, D))
    ABD6 = np.vstack((top, bot))

    rhs = np.concatenate([N_vec, M_vec])

    #solve 
    try:
        sol = np.linalg.solve(ABD6, rhs)
    except np.linalg.LinAlgError:
        raise ValueError("ABD Matrix is singular (cannot invert). Check laminate symmetry/properties")
    
    eps0 = sol[:3]
    kappa = sol[3:]
    return eps0,kappa

#==============================================================
                           #LOCAL STRES TRANSFORMATION
#=============================================================

def transform_stress_global_to_local(sigma_global, theta_deg):
    """"
    Converts global stress components (sx, sy, txy) into ply-local stress components (s1, s2, t12)

    Uses standard 2d rotation relations
    
    """

    theta = math.radians(theta_deg)
    m = math.cos(theta)
    n = math.sin(theta)
    sx, sy, txy = sigma_global

    s1 = m*m*sx + n*n*sy + 2*m*n*txy
    s2 = n*n*sx + m*m*sy - 2*m*n*txy
    t12 = -m*n*sx + m*n*sy + (m*m - n*n)*txy

    return np.array([s1, s2, t12])

#==============================================================
                           #FAILURE CRITERIA
#=============================================================

def failure_max_stress(sig, Xt, Xc, Yt, Yc, S):
    s1, s2, t12 = sig

    fi_ft = s1 / Xt if s1 >= 0 else 0
    fi_fc = -s1 / Xc if s1 < 0 else 0

    fi_mt = s2 / Yt if s2 >= 0 else 0
    fi_mc = -s2 / Yc if s2 < 0 else 0

    fi_s = abs(t12) / S
    
    return max(fi_ft, fi_fc, fi_mt, fi_mc, fi_s)

def failure_tsai_hill(sig, Xt, Xc, Yt, Yc, S):
    s1, s2, t12 = sig
    x = Xt if s1 >= 0 else Xc
    y = Yt if s2 >= 0 else Yc
    fi = (s1/x)**2 - (s1*s2)/(x**2) + (s2/y)**2 + (t12/S)**2
    return fi


def failure_tsai_wu(sig, Xt, Xc, Yt, Yc, S):
    s1, s2, t12 = sig

    f1 = 1/Xt - 1/Xc
    f2 = 1/Yt - 1/Yc
    f11 = 1/(Xt*Xc)
    f22 = 1/(Yt*Yc)
    f66 = 1/(S**2)

    f12 = -0.5 * np.sqrt(f11 * f22)

    fi = (f1*s1 + f2*s2 + f11*s1**2 + f22*s2**2 + 2*f12*s1*s2 + f66*t12**2)

    return fi

def failure_hashin(sig, Xt, Xc, Yt, Yc, S):
    s1, s2, t12 = sig

    #fibre tension
    if s1 >= 0:
        FI_ft = (s1/ Xt)**2 + (t12 / S)**2
    else:
        FI_ft = 0.0
    
    if s1 < 0:
        FI_fc = (s1 / Xc)**2
    else:
        FI_fc = 0.0

    if s2 >= 0:
        FI_mt = (s2 / Yt)**2 + (t12 / S) ** 2
    else:
        FI_mt = 0.0

    if s2 < 0:
        FI_mc = ((Yc / (2*S))**2 - 1) * (s2 / Yc) + (s2 / (2*S))**2 + (t12 / S)**2
    else:
        FI_mc = 0.0

    FI_values = {
        "Fibre Tension": FI_ft,
        "Fibre Compression": FI_fc,
        "Matrix Tension": FI_mt,
        "Matrix Compression": FI_mc
    }
    mode = max(FI_values, key=FI_values.get)
    return FI_values[mode], mode

#==============================================================
                           #FULL STRESS/STRAIN/FAILURE EXtRACTOR
#=============================================================

def compute_ply_strains_stresses(transformed_plies, A, B, D, N_vec, M_vec):
    """
    Computes:
    - mid plane strain and curvature
    - ply strains at mid-surface
    - global stresses
    - local stresses
    - all failure indices


    Returns
    --------
    eps0, kappa, ply_results: dict list
    
    """

    #get z coordinates
    total_t = sum(p['thickness'] for p in transformed_plies)
    z = [-total_t / 2.0]
    for p in transformed_plies:
        z.append(z[-1] + p['thickness'])    #z[i+1] is the top of ply i

#solve laminate response
    eps0, kappa = solve_midplane_strains(A, B, D, N_vec, M_vec)

    ply_results = []
    for k, ply in enumerate(transformed_plies):
        z_bot = z[k]
        z_top = z[k+1]
        z_mid = 0.5 * (z_bot + z_top)

        #strain at ply mid-surface (global coordinates)
        eps_mid = eps0 + z_mid * kappa
        Qb = np.array(ply['Q_bar'])

        sigma_global = Qb.dot(eps_mid)

        sigma_local = transform_stress_global_to_local(sigma_global, ply['angle'])

        #Failure criteria
        FI_max = failure_max_stress(
            sigma_local, ply["Xt"], ply["Xc"], ply["Yt"], ply["Yc"], ply["S"]
        )
        FI_th = failure_tsai_hill(
            sigma_local, ply["Xt"], ply["Xc"], ply["Yt"], ply["Yc"], ply["S"]
        )
        FI_tw = failure_tsai_wu(
            sigma_local, ply["Xt"], ply["Xc"], ply["Yt"], ply["Yc"], ply["S"]
        )
        FI_hashin, mode = failure_hashin(
            sigma_local, ply["Xt"], ply["Xc"], ply["Yt"], ply["Yc"], ply["S"]
        )

        SF_max = 1 / FI_max if FI_max > 0 else np.inf
        SF_th = 1 / FI_th if FI_th > 0 else np.inf
        SF_tw = 1 / FI_tw if FI_tw > 0 else np.inf
        SF_hashin = 1 / FI_hashin if FI_hashin > 0 else np.inf



        ply_results.append({
            "ply_index": k+1,
            "angle": ply['angle'],
            "thickness": ply['thickness'],
            "z_mid": z_mid,
            "eps_global": eps_mid,
            "sigma_global": sigma_global,
            "sigma_local": sigma_local,
            "FI_max": FI_max,
            "FI_th": FI_th,
            "FI_tw": FI_tw,
            "FI_hashin": FI_hashin,
            "Failure_mode": mode,
            "SF_max": SF_max,
            "SF_th": SF_th,
            "SF_tw": SF_tw,
            "SF_hashin": SF_hashin
        })

    return eps0, kappa, ply_results

#============================================
#            FIRST PLY FAILURE ANALYSIS
#============================================
 
def first_ply_failure(transformed_plies, A, B, D, N_vec, M_vec, max_load_factor = 10.0, step=0.01):
    """
    Incrementally scales applied loads until any ply reaches FI >= 1
    Searches up to 'max_load_factor' in increments of 'step'

    Returns
    ----------
    dict with:
        load_factor
        ply_index
        mode (failure mode)
        FI_values : dict



    """

    load_factor = 0.0
    while load_factor <= max_load_factor:
        scaled_N = N_vec * load_factor
        scaled_M = M_vec * load_factor

        eps0, kappa, ply_results = compute_ply_strains_stresses(
            transformed_plies, A, B, D, scaled_N, scaled_M
        )


        for p in ply_results:
            if (
                p["FI_max"] >=  1.0
                or p["FI_th"] >= 1.0
                or p["FI_tw"] >= 1.0
                or p["FI_hashin"] >= 1.0
            ):
                return {
                    "load_factor": load_factor,
                    "ply_index": p["ply_index"],
                    "mode": p["Failure_mode"],
                    "FI_values": {
                        "Max Stress": p["FI_max"],
                        "Tsai-Hill": p["FI_th"],
                        "Tsai-Wu": p["FI_tw"],
                        "Hashin": p["FI_hashin"]
                    }

                    
                }
            
        load_factor += step
    # if no failure occurred
    return {
        "load_factor": None,
        "ply_index": None,
        "mode": None,
        "FI_values": None
    }




