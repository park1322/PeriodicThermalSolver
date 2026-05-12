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
// Info<< "nBins (User Defined) = " << nBins << endl;
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

labelList binOfCell(mesh.nCells(), 0);

const vectorField& cellCenters = mesh.C().internalField();
forAll(cellCenters, cellI)
{
    const scalar x = cellCenters[cellI].x();
    // raw bin index
    label bini = label( floor((x - x0)/dX) );

    // clamp to [0, nBins-1] for safety against roundoff at x ~ xMax
    // bini = max(0, min(nBins-1, bini));

    binOfCell[cellI] = bini;
}// per each cell, which bin is it a member of?

// Info<< "binOfCell = " << binOfCell << nl << endl;

// This is list of list that stores cell IDs per bin
List< DynamicList<label> > binCells(nBins);
forAll(binOfCell, c) binCells[binOfCell[c]].append(c);

// for (label k=0; k<nBins; ++k)
// {
//     for (const label c : binCells[k])
//     {
//         // operate on cell c in bin k
//     }
// }

// fill internal field with raw bin indices (as scalars)
volScalarField checkBins
(
    IOobject
    (
        "checkBins",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimless, 0.0)   // start as 0 everywhere
);

forAll(checkBins.internalField(), cellI)
{
    checkBins[cellI] = scalar(binOfCell[cellI]) ;
}

forAll(checkBins.internalField(), cellI)
{
    label bin = binOfCell[cellI];
    if (bin % 2 == 0)  // even
    {
        checkBins[cellI] = -scalar(bin);
    }
    else
    {
        checkBins[cellI] = scalar(bin);
    }
}

checkBins.correctBoundaryConditions(); // harmless here

checkBins.write();


