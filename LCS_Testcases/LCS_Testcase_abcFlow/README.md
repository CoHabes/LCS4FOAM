#  Steady Arnold-Beltrami-Childress (ABC) flow

## About the testcase

### General info
 - The Arnold-Beltrami-Childress (ABC) flow is an exact periodic solution of the Euler equations and is often used in the literature to verify LCS calculation methods.
 - In order to test the function object on this flow configuration a dedicated ABC flow OpenFOAM solver (abcFlow) was written. This solver does not solve the Euler equations in the usual sense, but directly sets the velocity components on the given computational mesh.

### Mesh
 - A blockMesh with mesh size of 100 × 100 × 100 cells is used.
 - Since the mesh is rectlinear no additional LCS mesh is used.

### Running the case
 - Use the `Allrun` script to run the case
 - The case is setup to run on 2 cores but will allow for oversubscription if the number of cores can not be satisfied
 - The case is expected to run for about 6 minutes

## Inspecting the LCS results in ParaView
1. Open ParaView either with `paraFoam` of `paraView`
2. Click on "Open" to open new files
3. In your case folder go to the `cfd2lcs_output` directory
4. Select the group of files for the forward and/or backward FTLE results `fwdFTLE_-..dat` `bkwdFTLE_-..dat` and click on "OK"
5. The FTLE data will appear in your Pipeline Browser
6. Click on "Apply"
7. Select only the `fwdFTLE_-..dat` or the `bkwdFTLE_-..dat` to be visible
8. Select "Surface" as the "Representation" type of the data
9. Select the "FTLE" field for the "Coloring" of the surface
10. Depending on the selected LCS function object output options you may also be able to inspect the flow map fields "FM-X", "FM-Y", "FM-Z" and the velocity component fields "U-X", "U-Y", "U-Z"