#  2D Flow around cylinder with big additional LCS mesh

## About the testcase

### General info
 - The flow around an infinitely long cylinder with Re=200 which indicates that vortex shedding behind the cylinder
occurs in a barely laminar regime

### Mesh
 - A blockMesh composed out of multiple blocks surrounding the cylinder
 - Since the mesh is not rectlinear an additional LCS mesh is used which spans the whole domain

### Running the case
 - Use the `Allrun` script to run the case
 - The case is setup to run on 2 cores but will allow for oversubscription if the number of cores can not be satisfied
 - The case is expected to run for about 1 hour


## LCS Results

### Generation of the LCS results
 - The settings for the LCS function object can be found in the controlDict of this case
 - The LCS data will be stored in the `cfd2lcs_output` folder since your LCS write interval might not match with the write interval of your simulation
 - The data for forward and backward time FTLE fields will be named with the following pattern `fwdFTLE_<outputNumber>.dat` or the `bkwdFTLE_<outputNumber>.dat`
    - The outputNumber of the different files is smaller than 0 if the lcsIntegrationTime has not yet been reached´. This means that this data represents preliminary results. After reaching the lcsIntegrationTime the output files will be labeled with positive numbers since the LCS function object is able to compute the flow map for the whole integration interval after that point.
    - As an example the file `fwdFTLE_0.dat` represents the forward time FTLE field at time "lcsStartTime + lcsIntegrationTime" and the file `fwdFTLE_1.dat` represents the forward time FTLE field at time "lcsStartTime + lcsIntegrationTime + 1 * lcsWriteTimeIntervall".

> **_NOTE:_**  For more information about the different time intervals that play a role in the LCS computations please refer to the associated paper (DOI: TBD).


### Inspecting the LCS results in ParaView
1. Open ParaView either with `paraFoam` of `paraView`
2. Click on "Open" to open new files
3. In your case folder go to the `cfd2lcs_output` directory
4. Select the group of files for the forward and/or backward FTLE results that are numbered with positive values `fwdFTLE_<outputNumber>.dat` `bkwdFTLE_<outputNumber>.dat` and click on "OK"
5. The FTLE data will appear in your Pipeline Browser
6. Click on "Apply"
7. Select only the `fwdFTLE_..dat` or the `bkwdFTLE_..dat` to be visible
8. Select "Surface" as the "Representation" type of the data
9. Select the "FTLE" field for the "Coloring" of the surface
10. Depending on the selected LCS function object output options you may also be able to inspect the flow map fields "FM-X", "FM-Y", "FM-Z" and the velocity component fields "U-X", "U-Y", "U-Z"