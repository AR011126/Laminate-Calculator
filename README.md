Composite Laminate Analysis Toolkit

A Python-based engineering tool for analysing composite laminates using Classical Laminate Theory (CLT).
The toolkit calculates:

A, B, D stiffness matrices

Effective laminate engineering constants (Ex, Ey, vxy, Gxy)

Ply-level stresses, strains, and local stress transformations

Multiple failure criteria (Max Stress, Tsai-Hill, Tsai-Wu, Hashin)

First Ply Failure (FPF) predictions

Optional plotting of Failure Indices per ply

Flexible laminate definition modes (shared/different thickness & material properties)

This project is structured as an engineering portfolio piece and demonstrates implementation of composite mechanics from first principles.

Project Sturcture
Composite Calculator/
|
|--mechanics.py
|--central.py
|--README.md


Theory Summary
1) Ply stiffness
Each ply is modeled as an orthotropic lamina with properties (E1, E2, G12, ν12).
The ply stiffness matrix Q is computed using the plane-stress constitutive law.

2) Ply Rotation
For off-axis plies, the Q-matrix is rotated using standard tensor transformation.
The rotated Q̅ controls how each ply contributes to laminate stiffness.

3) Laminte Stiffness Matrix
A, B, and D matrices are assembled using thickness-wise integration:

A: In-plane stiffness

B: Coupling stiffness (bending ↔ extension)

D: Bending stiffness

4) Global Strains and Curvatures
The system is solved for midplane strain (epsilon0) and curvature (kappa):
[A  B] [epsilon] = [N]
[B  D] [kappa]   = [M]

5) Local Ply Stresses
Global stresses are rotated to local ply stresses (s1, s2, t12)

6) Failure Criteria
Four failure theories implemented:

Max Stress

Tsai-Hill

Tsai-Wu

Hashin (mode-specific)

7) First Ply Failure

Loads are scaled until any ply reaches FI ≥ 1.0.

