// ========================================================================= //
//                 UNIT TEST: 1D-3D-1D ROUND-TRIP CONSISTENCY
// ========================================================================= //
Info << nl << "--- STARTING PROJECTION ROUND-TRIP TEST (nBins = " << nBins << ") ---" << endl;

// 1. Create a Perfect 1D Sine Wave Array
scalarField exact1D(nBins, 0.0);
const scalar k_wave = 2.0 * constant::mathematical::pi / Lperiod;

for (label i = 0; i < nBins; ++i)
{
    scalar xNode = xMin + (i + 0.5) * dX;
    exact1D[i] = Foam::sin(k_wave * (xNode - xMin)) + 2.0;
}

// 2. Project 1D -> 3D using BOTH methods
// Assuming 'projector' is your BinProjectorX instance
tmp<volScalarField> tY_oldLinear = projector.project(exact1D, "Y_oldLinear");
tmp<volScalarField> tY_newConsistent = projector.projectConsistent(
    exact1D, xMinField, xMaxField, L_cellField, "Y_newConsistent"
);

volScalarField& Y_old = tY_oldLinear.ref();
volScalarField& Y_new = tY_newConsistent.ref();

// 3. Map 3D -> 1D (The Reverse Trip) using strict fractional binning
scalarField numOld(nBins, 0.0), numNew(nBins, 0.0), denSum(nBins, 0.0);

forAll(mesh.cells(), cellI)
{
    scalar xMinCell = xMinField[cellI];
    scalar xMaxCell = xMaxField[cellI];
    scalar L_cell   = L_cellField[cellI];
    scalar cellVol  = V[cellI];

    if (L_cell < SMALL) 
    {
        label b = floor((mesh.C()[cellI].x() - xMin) / dX);
        b = posMod(b, nBins);
        numOld[b] += Y_old[cellI] * cellVol;
        numNew[b] += Y_new[cellI] * cellVol;
        denSum[b] += cellVol;
        continue;
    }

    label startBin = floor((xMinCell - xMin) / dX);
    label endBin   = floor((xMaxCell - xMin) / dX);

    for (label b = startBin; b <= endBin; ++b) 
    {
        scalar binStart = xMin + b * dX;
        scalar binEnd   = binStart + dX;
        
        scalar overlapStart = max(xMinCell, binStart);
        scalar overlapEnd   = min(xMaxCell, binEnd);
        scalar fraction = max(0.0, overlapEnd - overlapStart) / L_cell;

        if (fraction > SMALL)
        {
            label wrappedBin = posMod(b, nBins);
            numOld[wrappedBin] += fraction * Y_old[cellI] * cellVol;
            numNew[wrappedBin] += fraction * Y_new[cellI] * cellVol;
            denSum[wrappedBin] += fraction * cellVol;
        }
    }
}

reduce(numOld, sumOp<scalarField>());
reduce(numNew, sumOp<scalarField>());
reduce(denSum, sumOp<scalarField>());

scalarField recovered1D_old = numOld / (denSum + SMALL);
scalarField recovered1D_new = numNew / (denSum + SMALL);

// 4. Calculate RMS Errors
scalar errOld = 0.0, errNew = 0.0;

for (label i = 0; i < nBins; ++i)
{
    errOld += Foam::sqr(recovered1D_old[i] - exact1D[i]);
    errNew += Foam::sqr(recovered1D_new[i] - exact1D[i]);
}

Info << "Round-Trip RMS Error (How much physics was corrupted by the projection):" << endl;
Info << "  -> Old Linear (Centroid-based):    " << Foam::sqrt(errOld / nBins) << endl;
Info << "  -> New Consistent (Fraction-based): " << Foam::sqrt(errNew / nBins) << endl;

// Dump for plotting
OFstream file("RoundTripTest_nBins" + name(nBins) + ".dat");
file << "x_node Exact RecoveredOld RecoveredNew" << endl;
for (label i = 0; i < nBins; ++i)
{
    scalar xNode = xMin + (i + 0.5) * dX;
    file << xNode << " " << exact1D[i] << " " << recovered1D_old[i] << " " << recovered1D_new[i] << endl;
}

Info << "--- END ROUND-TRIP TEST ---" << nl << endl;