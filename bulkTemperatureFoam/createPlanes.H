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


// above this for future. 

// const label nDiscrete = 11;
const label nDiscrete = readLabel(bulkDict.lookup("nDiscrete"));
const label nPlanes = nDiscrete + 1;         // edges for trapezoid
#include "boundBox.H"
const boundBox bb(mesh.bounds());   // bounds of the mesh points
const scalar xMin = bb.min().x();
const scalar xMax = bb.max().x();
const scalar dX = (xMax - xMin)/nDiscrete;
const scalar eps = 1e-9 * (xMax - xMin);  

auto nudgeIfBad = [&](scalar x)->scalar
{   
    x += eps;
    if(x > xMax) x -= 2*eps;

    // // nudge off exact ends
    // if (mag(x - xMin) < eps) x = xMin + eps;
    // if (mag(x - xMax) < eps) x = xMax - eps;

    // // also nudge off “nice” numbers like 0.0 (likely block interface)
    // if (mag(x) < eps) x += eps;

    return x;
};

// after createMesh.H
PtrList<sampledPlane> planes(nPlanes);
scalarField xPlanes(nPlanes);
scalarField areaPlanes(nPlanes);
for (label i=0; i<nPlanes; ++i) 
{
    scalar x = xMin + i*dX;
    xPlanes[i] = nudgeIfBad(x);

    planes.set
    (
        i, 
        new sampledPlane
        (
            "plane"+name(i), 
            mesh, 
            plane
            (
                point(xPlanes[i], 0, 0), 
                flowDir // normal vector for the plane
            )
        )
    );
    planes[i].update();  // one-time triangulation

    // tabulate the total areas of cut surfaces
    const scalarField& A = planes[i].magSf();
    areaPlanes[i] = gSum(A); 
}

// Obtain x coordinates at the discretization nodes
scalarField xDiscrete(nDiscrete);
for (label k=0; k<nDiscrete; ++k)
{
    xDiscrete[k] = 0.5 * (xPlanes[k] + xPlanes[k+1]); 
}

Info << nl << "Created " << nPlanes << " planes to have " << nDiscrete << " 1D nodes."<< nl << endl;

// // Set the pointers
// volScalarField Ux = mag(U & flowDir); 
// interpolationCellPoint<scalar> Uxip(Ux);
// // interpolationCellPoint<vector> Uip(U);
// interpolationCellPoint<scalar> Tip(T);

// volScalarField UxT = mag(U & flowDir) * T; 
// interpolationCellPoint<scalar> UxTip(UxT);
// need to set more pointers for the future. 




// Debugging & Checking 
// print per-plane total area
for (label i=0; i<nPlanes; ++i)
{
    const scalarField& A = planes[i].magSf();   // triangle areas
    const scalar Atot = gSum(A);                // parallel-safe sum
    Info<< "plane i=" << i
        << ",  x=" << xPlanes[i]
        << ",  area=" << Atot << nl;
}


// Info << nl << "Checking Areas : " << areaPlanes << nl << endl;

// for (label i=0; i<nPlanes; ++i)
// {
//     scalar numU = areaIntegralOnPlane(planes[i], Uxip);
//     scalar bulkU = numU/areaPlanes[i];
//     scalar numT = areaIntegralOnPlane(planes[i], UxTip);
//     scalar bulkT = numT / numU;

//     Info<< "plane i=" << i
//         << ",  x=" << xPlanes[i]
//         << ",  U_bulk=" << bulkU
//         << ",  T_bulk=" << bulkT << nl;
// }
