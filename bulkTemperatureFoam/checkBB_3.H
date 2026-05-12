// ========================================================================= //
//                 UNIT TEST: GEOMETRIC CROSS-SECTIONAL AREA
// ========================================================================= //
Info << nl << "--- STARTING GEOMETRIC AREA TEST (nBins = " << nBins << ") ---" << endl;

scalarField binArea(nBins, 0.0);

forAll(mesh.cells(), cellI)
{
    scalar xMinCell = xMinField[cellI];
    scalar xMaxCell = xMaxField[cellI];
    scalar L_cell   = L_cellField[cellI];
    scalar cellVol  = V[cellI];

    if (L_cell < SMALL) 
    {
        label b = posMod(floor((mesh.C()[cellI].x() - xMin) / dX), nBins);
        binArea[b] += cellVol / dX;
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
            // Accumulate physical volume, divide by dX to get Area
            binArea[wrappedBin] += (fraction * cellVol) / dX; 
        }
    }
}

reduce(binArea, sumOp<scalarField>());

// 2. Dump for plotting
OFstream file("AreaTest_nBins" + name(nBins) + ".dat");
file << "# x_node CrossSectionalArea" << endl;
for (label i = 0; i < nBins; ++i)
{
    scalar xNode = xMin + (i + 0.5) * dX;
    file << xNode << " " << binArea[i] << endl;
}

Info << "--- END GEOMETRIC AREA TEST ---" << nl << endl;