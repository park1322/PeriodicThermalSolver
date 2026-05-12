// Read config
IOdictionary bulkDict
(
    IOobject
    (
        "bulkMeanDict",
        runTime.system(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

const label  nBins = readLabel(bulkDict.lookup("nBins")); // discretization for the 1D nodes. 
#include "boundBox.H"
const boundBox bb(mesh.bounds());   // bounds of the mesh points
scalar xMin = bb.min().x();
scalar xMax = bb.max().x();
reduce(xMin, minOp<scalar>());
reduce(xMax, maxOp<scalar>());
const scalar Lperiod = xMax - xMin; 
const scalar dX = (xMax - xMin)/nBins;



// 1. Calculate Average Cell Size (The "Physical Limit")
scalar totalVol = gSum(mesh.V());
scalar nCells = returnReduce(mesh.nCells(), sumOp<label>());
scalar avgCellVol = totalVol / nCells;
scalar avgCellLen = Foam::pow(avgCellVol, 1.0/3.0);

// 2. Define Smoothing Length h

scalar h_grid = 1.0 * dX; // Option A: Tied to Grid; minimum requirement 
scalar h_mesh = 2.0 * avgCellLen; // Option B: Tied to Mesh (New way - stable on Poly mesh) 
// const scalar h = calibrateKernelRadius(mesh, nBins, flowDir, xMin, Lperiod);//max(h_grid, h_mesh); 
const scalar h = max(h_grid, h_mesh);
const scalar invH = 1.0 / h;

// 3. Calculate Span
const label kernelBinSpan = label(ceil(3.0 * h / dX)); // Quintic support is 3h 

Info << "  Global X Range: [" << xMin << ", " << xMax << "]" << endl;
Info << "  grid dX: " << dX << ", mesh avgLen: " << avgCellLen << endl;
Info << "  Selected h: " << h << " (Span: " << kernelBinSpan << " bins)" << endl;

// // Kernel Parameters
// scalar hRatio = 0.2; // ratio relative to bounding x 
// const label kernelSize = 3; 
// const scalar h = hRatio * Lperiod; // Standard setting for quintic support radius 3h
// const scalar invH = 1.0 / h;
// const label kernelBinSpan = label(ceil(kernelSize * h / dX));

// Field to store the 'center' 1D node for every cell
labelList cellCenterNode(mesh.nCells());

forAll(mesh.C(), cellI)
{
    scalar x = mesh.C()[cellI].x();
    scalar xRel = x - xMin; 
    // is periodic wrapping necessary? Don't think we are going to have any edge case...
    cellCenterNode[cellI] = round((xRel / dX) - 0.5); // Just store the center index
}

// forAll(mesh.C(), cellI)
// {
//     scalar x = mesh.C()[cellI].x();
//     scalar xRel = x - xMin; 

//     // Shift to align with Projector logic
//     scalar s_shifted = xRel - 0.5*dX;

//     // Wrap shifted coordinate
//     if (s_shifted < 0) s_shifted += Lperiod;

//     // Use FLOOR, just like the Projector
//     cellCenterNode[cellI] = label(floor(s_shifted / dX)); 
    
//     // Safety
//     if (cellCenterNode[cellI] >= nBins) cellCenterNode[cellI] = 0;
// }



