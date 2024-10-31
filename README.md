# Simulation_of_Le-vy_driven_OU_processes

This repository was developed for the final project of the course of Financial Engineering at Politecnico di Milano

The aim of this project is to explore the simulation of two broad categories of mean-reverting stochastic
processes: Ornstein-Uhlenbeck (OU-Lévy) and Lévy-Ornstein-Uhlenbeck (Lévy-OU). Specifically, we
focus on simulating Ornstein-Uhlenbeck Tempered Stable (OU-TS) and Tempered Stable Ornstein-
Uhlenbeck (TS-OU) processes. To achieve this, we employ both the Exact Decomposition algorithm
basedonthepropertyofself-decomposability(asdone in Sabino[4] and SabinoCufaro[5])andtheFast
General Monte Carlo method which requires only the characteristic function of the process (Baviera
& Manzoni [2]). The ultimate goal of our work is to utilize these simulations to price energy European
and American call options, since these kind of processes are able to capture the dynamics of energy
markets.

We have developed a library in MATLAB and in Python in order to meet the requirements for
the processes under study. For a better comprehension, the main script does not include all the
tests conducted, but it does fulfill the assignment’s requirements. We focused
our attention on the optimization in the MATLAB code. The results reported in the PDF document [Download the PDF](Simulation of Lévy-driven OU processes Report.pdf) are
obtained using MATLAB; in Python the results are available in the Jupyter notebook.

[1] Baviera, R. & Manzoni, P. (2024). Fast and General Simulation of Lévy-driven OU processes for Energy
Derivatives. Preprint arXiv: :2401.15483.
[2] Longstaff, F. A., & Schwartz, E. S. (2001). Valuing American options by simulation: a simple least-squares
approach. The review of financial studies, 14(1), 113-147.
[3] Sabino, P. (2022). Pricing energy derivatives in markets driven by tempered stable and CGMY processes of
Ornstein–Uhlenbeck type. Risks, 10(8):148.
[4] Sabino, P. & Cufaro Petroni, N. (2022). Fast simulation of tempered stable Ornstein–Uhlenbeck processes.
Computational Statistics, 37(5):2517–2551.
