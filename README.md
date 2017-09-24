# Numerical Modeling of 3D Porosity Waves

This is the repository for the project involving numerical modeling of 3D porosity waves. The simulation was visualized with VisIt software and can be viewed [here](https://drive.google.com/file/d/0B5JHF3TQfmV7NHRqR291Q29TVUk/view). The pressure wave is visualized on the left, with an underpressure and overpressure region, and the porosity wave is visualized on the right. They migrate upwards due to buoyancy and an initial Gaussian perturbation is provided.

### Project details

Porosity waves are an active area of research, and have potential impact for characterizing hydrocarbon migration, CO<sub>2</sub>  as well as for slow earthquakes. These are nonlinear solitary waves, and have been modeled here in a viscoelastic medium. 

The computation for obtaining porosity and pressure was performed using an explicit finite difference scheme with forward-Euler in space usng a 128 x 128 grid [here](porosity_3d.py). 


