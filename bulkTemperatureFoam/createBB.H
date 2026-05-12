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

const label  nBins   = readLabel(bulkDict.lookup("nBins"));

// Find the bounds 
#include "boundBox.H"
const boundBox bb(mesh.bounds());   // bounds of the mesh points
const scalar xMin = bb.min().x();
const scalar xMax = bb.max().x();

// Initialize 1D discretization info
const scalar Lperiod = xMax - xMin;             // your unit cell length
const scalar dX      = Lperiod / max(1, nBins);
const scalar x0      = xMin;                    // left edge of the period

// Report the initialization
Info<< nl << "* * * * 1D Discretization * * * * " << endl;
Info<< "nBins (User Defined) = " << nBins << endl;
Info<< "xMin = " << xMin << endl;
Info<< "xMax = " << xMax << endl;
Info<< "Lperiod = " << Lperiod << endl;
Info<< "x0 = " << x0 << endl;
Info<< "dX = " << dX << nl << endl;

scalarField binCenters(nBins);
for (label i=0; i<nBins; ++i)
{
    binCenters[i] = xMin + (i + 0.5)*dX;  // x-coordinates of bin centers
}

// We calculate the x bounds per cell beforehand
scalarField xMinField(mesh.nCells(), GREAT);
scalarField xMaxField(mesh.nCells(), -GREAT);

// forAll(mesh.cells(), cellI)
// {
//     scalar minX = GREAT;
//     scalar maxX = -GREAT;
    
//     const cell& c = mesh.cells()[cellI];
//     forAll(c, faceI) 
//     {
//         const face& f = mesh.faces()[c[faceI]];
//         forAll(f, pI) 
//         {
//             scalar px = mesh.points()[f[pI]].x();
//             minX = min(minX, px);
//             maxX = max(maxX, px);
//         }
//     }
//     xMinField[cellI] = minX;
//     xMaxField[cellI] = maxX;
// }


// faster alternative 
const labelListList& cellPts = mesh.cellPoints();

forAll(cellPts, cellI)
{
    scalar minX = GREAT;
    scalar maxX = -GREAT;
    
    // 2. Grab the list of unique point indices for THIS cell
    const labelList& cPoints = cellPts[cellI];
    
    // 3. Loop only through these unique points
    forAll(cPoints, pI)
    {
        // cPoints[pI] is the global point index
        scalar px = mesh.points()[cPoints[pI]].x();
        
        minX = min(minX, px);
        maxX = max(maxX, px);
    }
    
    xMinField[cellI] = minX;
    xMaxField[cellI] = maxX;
}

// Precompute the lengths to avoid repeated subtractions
scalarField L_cellField = xMaxField - xMinField;