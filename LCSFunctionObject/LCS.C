/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "LCS.H"
#include "dictionary.H"
#include "meshToMesh.H"
#include "fvMeshSubset.H"
#include "cellSet.H"
#include "typeInfo.H"
#include "wallPolyPatch.H"
#include "emptyPolyPatch.H"
#include "symmetryPolyPatch.H"
#include "symmetryPlanePolyPatch.H"
#include "wedgePolyPatch.H"
#include "cyclicPolyPatch.H"
#include "processorPolyPatch.H"
#include "polyMesh.H"
#include "block.H"
#include "blockDescriptor.H"
#include "addToRunTimeSelectionTable.H"
#include <memory>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(LCS, 0);
    addToRunTimeSelectionTable(functionObject, LCS, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::LCS::LCS
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    name_(name),
    lcsMeshPtr_(nullptr),
    lcsMeshSubset_(nullptr),
    meshToMeshPtr_(nullptr),
    lcsVelFieldPtr_(nullptr),
    controlDict_(time_.controlDict()),
    globalBb_(),
    localBb_(),
    active_(true),
    firstExe_(true),
    n_(),
    offset_(),
    x_(nullptr),
    y_(nullptr),
    z_(nullptr),
    u_(nullptr),
    v_(nullptr),
    w_(nullptr),
    flag_(nullptr),
    id_fwd_(-1),
    id_bkwd_(-1),
    isOverset_(),
    lcsOversetRegion_(),
    isStaticRectLinear_(),
    globalN_(),
    uName_(),
    res_(),
    ftleFwd_(),
    ftleBkwd_(),
    t_start_(),
    t_end_(),
    T_(),
    H_(),
    lcsOpts_()
{
    // Check if the available mesh is an fvMesh, otherwise deactivate
    if (!isA<fvMesh>(mesh_))
    {
        active_ = false;
        WarningInFunction
            << "No fvMesh available, deactivating " << name_ << nl
            << endl;
    }

    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::LCS::~LCS()
{
    if(x_){
        delete[] x_;
        x_ = nullptr;
    }
    if(y_){
        delete[] y_;
        y_ = nullptr;
    }
    if(z_){
        delete[] z_;
        z_ = nullptr;
    }
    if(u_){
        delete[] u_;
        u_ = nullptr;
    }
    if(v_){
        delete[] v_;
        v_ = nullptr;
    }
    if(w_){
        delete[] w_;
        w_ = nullptr;
    }
    if(flag_){
        delete[] flag_;
        flag_ = nullptr;
    }

    if(!Pstream::parRun())
    {
        MPI_Finalize();
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::LCS::initialize()
{
    if (active_)
    {
        Info << "-----"
             << "libcfd2lcs initialisation "
             << "-----------------------------------"
             << endl;

        // Making serial run possible
        if(!Pstream::parRun())
        {
            MPI_Init(NULL, NULL);
        }

        // Mpi communicator
        MPI_Comm comm = MPI_COMM_WORLD;

        // Set lcs mesh pointer
        if(isStaticRectLinear_)
        {
            lcsMeshPtr_.reset(&mesh_);

            if(isOverset_){
                createLCSSubsetMesh();
            }

        }else
        {
            // TODO: generate/read lcsmesh, generate meshToMesh interpolator
            createLCSMesh();
        }

        // Compute local and global bounding boxes for lcs mesh
        getBoundBoxes();

        // Compute number of grid points for THIS partition in
        // x=i, y=j and z=k direction
        getNumberOfCellsInDirection();

        // Compute cell offset of local mesh in respect to the global mesh
        getOffset();

        // Allocate space for data used by lcs library
        x_ = new lcsdata_t[n_[0]*n_[1]*n_[2]];
        y_ = new lcsdata_t[n_[0]*n_[1]*n_[2]];
        z_ = new lcsdata_t[n_[0]*n_[1]*n_[2]];
        u_ = new lcsdata_t[n_[0]*n_[1]*n_[2]];
        v_ = new lcsdata_t[n_[0]*n_[1]*n_[2]];
        w_ = new lcsdata_t[n_[0]*n_[1]*n_[2]];
        flag_ = new int[n_[0]*n_[1]*n_[2]];

        // convert OpenFoam mesh cell centers to x,y,z arays
        // and set boundary type flags for each cell center
        getCellCenterCoords();

        // Initializes the communications and data storage for
        // OpenFoam data input
        cfd2lcs_init_c(comm,n_.data(),offset_.data(),x_,y_,z_,flag_);

        // Set CFD2LCS options/parameters
        lcsOpts_.setLCSoptions();

        // Initialize LCS diagnostics
        initializeLCSDiagnostics();

        // start() has been called
        firstExe_ = false;

        Info << "---------------------------------"
             << "---------------------------------"
             << endl;
    }

    return true;
}

bool Foam::functionObjects::LCS::read(const dictionary& dict)
{
    if (active_)
    {
        isOverset_ = dict.lookupOrDefault<Switch>("isOverset", false);

        if(isOverset_)
        {
            lcsOversetRegion_ = word(dict.lookup("lcsOversetRegion"));

            Info << "LCSFunctionObject: Using cellZone " << lcsOversetRegion_
                 <<" in oversetMesh for lcs computations" << nl
                 << endl;
        }

        dict.lookup("isStaticRectLinear") >> isStaticRectLinear_;

        if((!isStaticRectLinear_) && isOverset_)
        {
            FatalErrorInFunction
                << "a non static rectlinear zellZone can not be used in the "
                << "oversetMesh for lcs computation." << nl
                << "set isStaticRectLinear to true or do not use an oversetMesh"
                << exit(FatalError);
        }
        else if(isStaticRectLinear_ && !isOverset_)
        {
            Info << "LCSFunctionObject: Using cfd mesh for lcs computations\n"
                 << endl;
        }
        else if(!isStaticRectLinear_ && !isOverset_)
        {
            Info << "LCSFunctionObject: Using desiganted lcs mesh for lcs "
                 << "computations" << nl
                 << endl;
        }

        readGlobalN();

        uName_ = dict.lookupOrDefault<word>("velocityName", "U");
        ftleFwd_ = dict.lookupOrDefault<Switch>("ftleFwd", true);
        ftleBkwd_ = dict.lookupOrDefault<Switch>("ftleBkwd", false);

        // check if any of the diagnostics should be executed
        if(!ftleFwd_ && !ftleBkwd_)
        {
            active_ = false;
            WarningInFunction
                << "No FTLE evaluation activated, deactivating " << name_ << nl
                << endl;
        }

        res_ = dict.lookupOrDefault<label>("resolution", 0);

        scalar simulationStartTime =
            readScalar(controlDict_.lookup("startTime"));
        scalar simulationEndTime =
            readScalar(controlDict_.lookup("endTime"));

        t_start_ =
            dict.lookupOrDefault<scalar>("lcsStartTime", simulationStartTime);
        t_end_ =
            dict.lookupOrDefault<scalar>("lcsEndTime", simulationEndTime);

        T_ = readScalar(dict.lookup("lcsIntegrationTime"));

        // determine default value for LCS output time interval
        scalar simulationWriteTimeInterval = 1.0f;
        word writeControl(controlDict_.lookup("writeControl"));

        if(writeControl == "runTime" || writeControl == "adjustableRunTime" )
        {
            simulationWriteTimeInterval =
                readScalar(controlDict_.lookup("writeInterval"));
        }
        else if(writeControl == "timeStep")
        {
            Switch adjustTimeStep =
                controlDict_.lookupOrDefault<Switch>("adjustTimeStep", false);

            if(adjustTimeStep)
            {
                WarningInFunction
                    << "Adjustable simulation time steps with timeStep "
                    << "writeControl is not support with LCS diagnostics." << nl
                    << "If no lcsWriteInterval is provided "
                    << "lcsWriteTimeInterval is set to 1s. " << name_ << nl
                    << endl;
            }
            else
            {
                scalar simulationWriteStepInterval =
                    readScalar(controlDict_.lookup("writeInterval"));
                scalar simulationDeltaT =
                    readScalar(controlDict_.lookup("deltaT"));

                simulationWriteTimeInterval =
                    simulationDeltaT * simulationWriteStepInterval;
            }
        }
        else
        {
            WarningInFunction
                << "cpuTime and clockTime writeControl is not support with "
                << "LCS diagnostics." << nl
                << "If no lcsWriteInterval is provided lcsWriteTimeInterval "
                << "is set to 1s. " << name_ << nl
                << endl;
        }

        H_ =
            dict.lookupOrDefault<scalar>
            (
                "lcsSubStepWriteIntervall",
                simulationWriteTimeInterval
            );

        // read integrator
        // - Maybe better with selectionTables, but this works for now
        std::map<std::string, int> lcsIntegratorMap
        {
            {"euler", EULER},
            {"trapezodial", TRAPEZOIDAL},
            {"rk2", RK2},
            {"rk3", RK3},
            {"rk4", RK4}
        };
        word integrator = dict.lookupOrDefault<word>("integrator", "rk2");
        auto searchIntegrator = lcsIntegratorMap.find(integrator);

        if (searchIntegrator != lcsIntegratorMap.end())
        {
            lcsOpts_.integrator = searchIntegrator->second;
        }
        else
        {
            WarningInFunction
                << "LCS integrator " << integrator
                << "is no valid integration scheme for lcs diagnostic" << nl
                << "valid types are: euler, trapezodial, rk2, rk3, rk4. " << nl
                << "Setting LCS integrator to rk2. " << name_ << nl
                << endl;
            lcsOpts_.integrator = RK2;
        }


        // read interpolator
        // - Maybe better with selectionTables, but this works for now
        std::map<std::string, int> lcsInterpolatorMap
        {
            {"nearestNBR", NEAREST_NBR},
            {"linear", LINEAR},
            {"quadratic", QUADRATIC},
            {"cubic", CUBIC},
            {"tse", TSE},
            {"tseLimit", TSE_LIMIT}
        };

        word interpolator =
            dict.lookupOrDefault<word>("interpolator", "linear");

        auto searchInterpolator = lcsInterpolatorMap.find(interpolator);

        if (searchInterpolator != lcsInterpolatorMap.end())
        {
            lcsOpts_.interpolator = searchInterpolator->second;
        }
        else
        {
            WarningInFunction
                << "LCS interpolator " << interpolator
                << "is no valid interpolation scheme for lcs diagnostic." << nl
                << "Valid types are: "
                << "nearestNBR, linear, quadratic, cubic, tse, tseLimit. " << nl
                << "Setting LCS interpolator to linear. " << name_ << nl
                << endl;

            lcsOpts_.interpolator = LINEAR;
        }

        // read lcs CFl number
        lcsOpts_.CFL = dict.lookupOrDefault<scalar>("lcsCFL", 0.5f);

        // read lcs options
        const dictionary& lcsOptionsDict = dict.subDict("lcsOptions");

        Switch synctimer =
            lcsOptionsDict.lookupOrDefault<Switch>("synctimer", false);
        lcsOpts_.synctimer =
            (synctimer) ? LCS_TRUE : LCS_FALSE;

        Switch debug =
            lcsOptionsDict.lookupOrDefault<Switch>("debug", false);
        lcsOpts_.debug =
            (debug) ? LCS_TRUE : LCS_FALSE;

        Switch writeFlowmap =
            lcsOptionsDict.lookupOrDefault<Switch>("writeFlowmap", false);
        lcsOpts_.writeFlowmap =
            (writeFlowmap) ? LCS_TRUE : LCS_FALSE;

        Switch writeBCFalgs =
            lcsOptionsDict.lookupOrDefault<Switch>("writeBCFlags", false);
        lcsOpts_.writeBCFalgs =
            (writeBCFalgs) ? LCS_TRUE : LCS_FALSE;

        Switch incompressible =
            lcsOptionsDict.lookupOrDefault<Switch>("incompressible", false);
        lcsOpts_.incompressible =
            (incompressible) ? LCS_TRUE : LCS_FALSE;

        Switch auxGrid =
            lcsOptionsDict.lookupOrDefault<Switch>("auxGrid", false);
        lcsOpts_.auxGrid = (auxGrid) ? LCS_TRUE : LCS_FALSE;

    }

    return true;
}

bool Foam::functionObjects::LCS::execute()
{
    if(firstExe_)
    {
        initialize();
    }

    scalar time = mesh_.time().value();

    if (active_ && time >= t_start_ && time <= t_end_)
    {
        getVelocityField();
        cfd2lcs_update_c(n_.data(), u_, v_, w_ ,time);
    }

    return true;
}

bool Foam::functionObjects::LCS::end()
{
    if (id_fwd_ != -1)
    {
        cfd2lcs_diagnostic_destroy_c(id_fwd_);
        id_fwd_ = -1;
    }

    if (id_bkwd_ != -1)
    {
        cfd2lcs_diagnostic_destroy_c(id_bkwd_);
        id_bkwd_ = -1;
    }

    return true;
}

bool Foam::functionObjects::LCS::write()
{
    return true;
}

void Foam::functionObjects::LCS::readGlobalN()
{
    const word dictName("blockMeshDict");
    word regionDir;

    if (!isStaticRectLinear_)
    {
        // system/<region>/blockMeshDict
        const word regionName = "LCS";
        regionDir = polyMesh::regionName(regionName);
    }
    else
    {
        // constant/polyMesh/blockMeshDict
        const word regionName  = polyMesh::defaultRegion;
        regionDir = polyMesh::regionName(regionName);
    }

    fileName dictPath = time_.system()/regionDir/dictName;

    Info << time_.caseSystem() << endl;
    Info << time_.system() << endl;
    Info << regionDir << endl;
    Info << dictName << endl;

    //- Read blockMeshDict of LCS mesh
    IOdictionary meshDict
    (
        IOobject
        (
            dictPath,
            time_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    //- Point field defining the block mesh (corners)
    //- needed for block constructor
    pointField blockPointField(meshDict.lookup("vertices"));

    //- Empty list of block edges and faces
    //- Needed for block constructor
    blockEdgeList edges;
    blockFaceList faces;

    //- Block index
    label index = 0;

    //- Blocks input stream
    ITstream& is(meshDict.lookup("blocks"));

    token firstToken(is);

    if (!firstToken.isLabel())
    {
        is.putBack(firstToken);
    }

    // Read beginning of blocks
    is.readBegin("blocks");

    token lastToken(is);
    while
    (
       !(
           lastToken.isPunctuation()
        && lastToken.pToken() == token::END_LIST
        )
    )
    {
        is.putBack(lastToken);

        block meshblock
        (
            meshDict,
            index++,
            blockPointField,
            edges,
            faces,
            is
        );

        if (isOverset_ && meshblock.zoneName() != lcsOversetRegion_)
        {
            is >> lastToken;
            continue;
        }
        else
        {
            Vector<label> globalNv = meshblock.density();

            globalN_.clear();
            globalN_.append(globalNv.x());
            globalN_.append(globalNv.y());
            globalN_.append(globalNv.z());
        }

        is >> lastToken;
    }

    is.putBack(lastToken);

    // Read end of blocks
    is.readEnd("blocks");
}

void Foam::functionObjects::LCS::getBoundBoxes()
{
    if(!isOverset_)
    {
        globalBb_ = lcsMeshPtr_().bounds();
        localBb_ = boundBox(lcsMeshPtr_().points(), false);
    }
    else
    {
        // make subsetMesh for boundingBox calculation of lcsOversetRegion
        globalBb_ = boundBox(lcsMeshSubset_().subMesh().points(), true);
        localBb_ = boundBox(lcsMeshSubset_().subMesh().points(), false);
    }

    Info << "global bounding box:" << globalBb_ << endl;
    Pout<< "local bounding box:" << localBb_ << endl;
}

void Foam::functionObjects::LCS::getCellCenterCoords()
{
    // Cell centroid coordinates
    const vectorField& centres =
        (!isOverset_) ?
            lcsMeshPtr_().C().internalField() :
            lcsMeshSubset_().subMesh().C().internalField();

    // Loop over cell centres
    forAll (centres, celli)
    {
        x_[celli] = centres[celli].component(0);
        y_[celli] = centres[celli].component(1);
        z_[celli] = centres[celli].component(2);

        //set all cells als internal
        flag_[celli] = LCS_INTERNAL;
    }

    const polyBoundaryMesh& boundaryMesh =
        (!isOverset_) ?
            lcsMeshPtr_().boundaryMesh():
            lcsMeshSubset_().subMesh().boundaryMesh();

    // Loop over all boundary patches
    forAll (boundaryMesh, patchi)
    {
        // Current poly patch
        const polyPatch& pp = boundaryMesh[patchi];

        // determine boundary type
        label boundaryType = LCS_OUTFLOW;

        if (isA<wallPolyPatch>(pp))
        {
            boundaryType = LCS_SLIP;
        }
        else if
        (
               isA<emptyPolyPatch>(pp)
            || isA<symmetryPolyPatch>(pp)
            || isA<symmetryPlanePolyPatch>(pp)
            || isA<wedgePolyPatch>(pp)
            || isA<cyclicPolyPatch>(pp)
            || isA<processorPolyPatch>(pp)
        )
        {
            boundaryType = LCS_INTERNAL;
        }

        // Loop over all faces of boundary patch
        forAll(boundaryMesh[patchi], facei)
        {
            // Boundary cell index
            const label& bCell = boundaryMesh[patchi].faceCells()[facei];
            if
            (
                   boundaryType != LCS_INTERNAL
                && flag_[bCell] != LCS_WALL
                && flag_[bCell] != LCS_SLIP)
            {
                flag_[bCell] = boundaryType;
            }
        }
    }
}

void Foam::functionObjects::LCS::getVelocityField()
{
    Info << "get velocity field for LCS computation" << endl;

    // pointer to const velocity field that is used for lcs computations
    std::shared_ptr<const vectorField> uLcsInPtr;

    if (mesh_.foundObject<volVectorField>(uName_))
    {
        const volVectorField& uCfd = mesh_.lookupObject<volVectorField>(uName_);

        if (isOverset_)
        {
            // map velocity field from global overSetMesh to subsetMesh
            // which is used for the lcs computation
            uLcsInPtr.reset();
            uLcsInPtr =
                std::make_shared<const vectorField>
                (
                    lcsMeshSubset_->interpolate(uCfd)().internalField()
                );

        }
        else if(!isStaticRectLinear_)
        {
            volVectorField& uLCS = lcsVelField();

            // interpolate cfd velocity field to lcs velocity field
            meshToMeshInterp().mapSrcToTgt
            (
                uCfd.primitiveField(),
                plusEqOp<vector>(),
                uLCS.primitiveFieldRef()
            );

            uLcsInPtr.reset();
            uLcsInPtr =
                std::make_shared<const vectorField>
                (
                    uLCS.internalField()
                );
        }
        else
        {
            const vectorField& uCfdi = uCfd.internalField();

            uLcsInPtr.reset();
            uLcsInPtr =
                std::make_shared<const vectorField>
                (
                    uCfdi
                );
        }

        // store velocity field in ordered array for each velocity component
        if (!uLcsInPtr->empty())
        {
            const vectorField& uLcsIn = *uLcsInPtr;

            forAll (uLcsIn, celli)
            {
                u_[celli] = uLcsIn[celli].component(0);
                v_[celli] = uLcsIn[celli].component(1);
                w_[celli] = uLcsIn[celli].component(2);
            }
        }
    }
    else
    {
        FatalErrorIn
        (
            "Foam::functionObjects::LCS::getVelocityField()"
        )   << "Velocity field with name " << uName_ << " not found"
            << exit(FatalError);
    }
}

void Foam::functionObjects::LCS::getNumberOfCellsInDirection()
{
    n_.setSize(3);

    // Assuming that LCS mesh has constant cell size along each axis
    n_[0] =
        round(globalN_[0]
      * (localBb_.max().component(0) - localBb_.min().component(0))
      / (globalBb_.max().component(0) - globalBb_.min().component(0)));

    n_[1] =
        round(globalN_[1]
      * (localBb_.max().component(1) - localBb_.min().component(1))
      / (globalBb_.max().component(1) - globalBb_.min().component(1)));

    n_[2] =
        round(globalN_[2]
      * (localBb_.max().component(2) - localBb_.min().component(2))
      / (globalBb_.max().component(2) - globalBb_.min().component(2)));

    Pout << "Number of cells in x:" << n_[0]
         << " y:" << n_[1]
         << " z:" << n_[2]
         << endl;
}

void Foam::functionObjects::LCS::getOffset()
{
    offset_.setSize(3);

    // Assuming that LCS mesh has constant cell size along each axis
    offset_[0] =
        round(n_[0]
      * (localBb_.min().component(0) - globalBb_.min().component(0))
      / (localBb_.max().component(0) - localBb_.min().component(0)));

    offset_[1] =
        round(n_[1]
      * (localBb_.min().component(1) - globalBb_.min().component(1))
      / (localBb_.max().component(1) - localBb_.min().component(1)));

    offset_[2] =
        round(n_[2]
      * (localBb_.min().component(2) - globalBb_.min().component(2))
      / (localBb_.max().component(2) - localBb_.min().component(2)));

    Pout << "Offset in x:" << offset_[0]
         << " y:" << offset_[1]
         << " z:" << offset_[2]
         << endl;
}

void Foam::functionObjects::LCS::initializeLCSDiagnostics()
{
    if(ftleFwd_)
    {
        char labelfwd[LCS_NAMELEN]="fwdFTLE";

        id_fwd_ =
            cfd2lcs_diagnostic_init_c
            (
                FTLE_FWD,
                res_,
                T_,
                H_,
                labelfwd
            );
    }

    if(ftleBkwd_)
    {
        char labelbkwd[LCS_NAMELEN]="bkwdFTLE";

        id_bkwd_ =
            cfd2lcs_diagnostic_init_c
            (
                FTLE_BKWD,
                res_,
                T_,
                H_,
                labelbkwd
            );
    }
}

void Foam::functionObjects::LCS::lcsOpts::setLCSoptions()
{
    char option1[] = "SYNCTIMER";
    cfd2lcs_set_option_c(option1, synctimer);

    char option2[] = "DEBUG";
    cfd2lcs_set_option_c(option2, debug);

    char option3[] = "WRITE_FLOWMAP";
    cfd2lcs_set_option_c(option3, writeFlowmap);

    char option4[] = "WRITE_BCFLAG";
    cfd2lcs_set_option_c(option4, writeBCFalgs);

    char option5[] = "INCOMPRESSIBLE";
    cfd2lcs_set_option_c(option5, incompressible);

    char option6[] = "AUX_GRID";
   cfd2lcs_set_option_c(option6, auxGrid);

    char option7[] = "INTEGRATOR";
    cfd2lcs_set_option_c(option7, integrator);

    char option8[] = "INTERPOLATOR";
    cfd2lcs_set_option_c(option8, interpolator);

    char option9[] = "CFL";
    cfd2lcs_set_param_c(option9, CFL);
}

void Foam::functionObjects::LCS::createLCSMesh()
{
    //- designated lcs mesh
    lcsMeshPtr_.reset
    (
        new fvMesh
        (
            IOobject
            (
                "LCS",
                mesh_.time().timeName(),
                mesh_.time(),
                IOobject::MUST_READ
            )
        )
    );

    Info << "Created LCS Mesh" << endl;
}

Foam::volVectorField& Foam::functionObjects::LCS::lcsVelField()
{
    if(!lcsVelFieldPtr_.good())
    {
        lcsVelFieldPtr_.reset
        (
            new volVectorField
            (
                IOobject
                (
                    "U",
                    mesh_.time().timeName(),
                    lcsMeshPtr_(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                lcsMeshPtr_(),
                dimensionedVector
                (
                    "zero",
                    dimensionSet(0,1,-1,0,0),
                    vector::zero
                )
            )
        );
    }

    return lcsVelFieldPtr_();
}

const Foam::meshToMesh& Foam::functionObjects::LCS::meshToMeshInterp()
{
    if (!meshToMeshPtr_.good())
    {
        // TODO: make interpolation for inconsistent meshes possible
        // for now only consistent meshes
        HashTable<word> patchMap;
        wordList cuttingPatches;

        meshToMeshPtr_.reset
        (
            new meshToMesh
            (
                mesh_,
                lcsMeshPtr_(),
                meshToMesh::interpolationMethod::imCorrectedCellVolumeWeight,
                patchMap,
                cuttingPatches
            )
        );
    }

    return meshToMeshPtr_();
}

void Foam::functionObjects::LCS::createLCSSubsetMesh()
{
    cellSet lcsRegionSet(lcsMeshPtr_(), lcsOversetRegion_);

    lcsMeshSubset_.reset
    (
        new fvMeshSubset
        (
            lcsMeshPtr_(),
            lcsRegionSet,
            -1,
            true
        )
    );
}
// ************************************************************************* //
