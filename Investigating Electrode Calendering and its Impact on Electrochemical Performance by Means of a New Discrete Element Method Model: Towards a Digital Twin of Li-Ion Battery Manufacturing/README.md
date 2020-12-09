In here you will find the codes for modelling in 3D (x, y, z) at the mesoscale (tens of µms) the electrode slurry based on NMC as active material, the associated dried electrode and the calendering process. 
To launch all these code it will be needed the open source software LAMMPS (initial random structure, slurry and dried electrode) or LIGGGHTS (calendering).

- The code "in_structure.run" allows to build the initial randomly generated structure to match the volume fraction of active material (AM) and carbon binder domain (CBD) and the particle size distributiof of the AM used experimentally;
- The code "slurry.run" allows to model the slurry phase and to calculate its density (associated to the force fields parameters values used);
- The code "viscosity_Shear-rate.run" allows to calculate the slurry viscosity for a specific shear-rate and temperature (associated to the FFs parameters values used);
- The code "dried_electrode.run" allows to calculate the electrode mesostructure after the drying step. In addition, it can output the associated electrode density and an estimation of the associated electrode porosity (associated to the FFs parameters values used);
- The code "calendering.run" allows to mimic the electrode calendering and to control the degree of calendering desired.

These codes are offered for smaller (Small_electrodes folder) and bigger (Thick_electrodes folder) electrodes, for the sake of examples. However, the dimension of the final electrode can be modulated as a function of the numbar of particles used through the "in_structure.run" code.

For all the details about these simulations, please check the following scientific articles:

- Investigating Electrode Calendering and its Impact on Electrochemical Performance by Means of a New Discrete Element Method Model: Towards a Digital Twin of Li-Ion Battery Manufacturing by Alain C. Ngandjong, Teo Lombardo, Emiliano N. Primo, Mehdi Chouchane, Abbos Shodiev, Oier Arcelus, Alejandro A. Franco.
- Accelerated Optimization Methods for Force‐Field Parametrization in Battery Electrode Manufacturing Modeling by Teo Lombardo,  Dr. Jean‐Baptiste Hoock,  Dr. Emiliano N. Primo,  Dr. Alain C. Ngandjong,  Marc Duquesnoy and Prof. Alejandro A. Franco.  https://doi.org/10.1002/batt.202000049
- Multiscale Simulation Platform Linking Lithium Ion Battery Electrode Fabrication Process with Performance at the Cell Level by Alain C. Ngandjong, Alexis Rucci, Mariem Maiza, Garima Shukla, Jorge Vazquez-Arenas, and Alejandro A. Franco. https://doi.org/10.1021/acs.jpclett.7b02647
- Tracking variabilities in the simulation of Lithium Ion Battery electrode fabrication and its impact on electrochemical performance by AlexisRucci, Alain C.Ngandjong, Emiliano N. Primo, Mariem Maiza, Alejandro A.Franco. https://doi.org/10.1016/j.electacta.2019.04.110

Please, to reproduce exactly the results reported in the afroementioned articles you should use the force fields parameters values reported in the associated scientific article and the procedure described in there.


For any doubt do not hesitate to contact us (we will be happy to answer and help you!) :
Teo Lombardo --> teo.lombardo@u-picardie.fr *OR* teo.lombardo3@gmail.com
Alain Ngandjong --> alain.ngandjong@u-picardie.fr
Alejandro A. Franco --> alejandro.franco@u-picardie.fr