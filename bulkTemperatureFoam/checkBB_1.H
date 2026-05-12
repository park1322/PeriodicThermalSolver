// ========================================================================= //
//                 UNIT TEST: CELL-FRACTION BINNING & DERIVATIVE
// ========================================================================= //
Info << nl << "--- STARTING CELL-FRACTION STRESS TEST (nBins = " << nBins << ") ---" << endl;

// 1. Setup: Sine Wave + White Noise
const scalar k_wave = 2.0 * constant::mathematical::pi / Lperiod;
const scalar noiseAmp = 0.2; // 20% Noise 

Random rndGen(1234); 

// Arrays to accumulate the sums
scalarField numClean(nBins, 0.0);
scalarField numNoisy(nBins, 0.0);
scalarField denSum(nBins, 0.0);

// 2. Cell-Fraction Binning Loop
forAll(mesh.cells(), cellI)
{
    // Generate signals
    scalar xCenter = mesh.C()[cellI].x();
    scalar exactVal = Foam::sin(k_wave * (xCenter - xMin)) + 2.0;
    scalar noise = (rndGen.sample01<scalar>() - 0.5) * 2.0 * noiseAmp;
    scalar noisyVal = exactVal + noise;
    
    scalar cellVolume = V[cellI];
    
    // Grab precomputed bounds
    scalar xMinCell = xMinField[cellI];
    scalar xMaxCell = xMaxField[cellI];
    scalar L_cell   = L_cellField[cellI];
    
    // Fallback for flat cells (e.g. baffles with zero x-thickness)
    if (L_cell < SMALL) 
    {
        label b = floor((xCenter - xMin) / dX);
        b = posMod(b, nBins);
        numClean[b] += exactVal * cellVolume;
        numNoisy[b] += noisyVal * cellVolume;
        denSum[b]   += cellVolume;
        continue;
    }
    
    // Determine overlapping bins
    label startBin = floor((xMinCell - xMin) / dX);
    label endBin   = floor((xMaxCell - xMin) / dX);
    
    for (label b = startBin; b <= endBin; ++b) 
    {
        scalar binStart = xMin + b * dX;
        scalar binEnd   = binStart + dX;
        
        scalar overlapStart = max(xMinCell, binStart);
        scalar overlapEnd   = min(xMaxCell, binEnd);
        scalar overlapLength = max(0.0, overlapEnd - overlapStart);

        if (overlapLength < SMALL) 
        {
            Info << " Very small overlap? cell i: " << cellI << endl;
        }
        
        scalar fraction = overlapLength / L_cell;
        
        if (fraction > SMALL)
        {
            label wrappedBin = posMod(b, nBins);
            numClean[wrappedBin] += fraction * exactVal * cellVolume;
            numNoisy[wrappedBin] += fraction * noisyVal * cellVolume;
            denSum[wrappedBin]   += fraction * cellVolume;
        }
    }
}

// Parallel Reduction
reduce(numClean, sumOp<scalarField>());
reduce(numNoisy, sumOp<scalarField>());
reduce(denSum, sumOp<scalarField>());

// Calculate Final 1D Profiles
scalarField bulkClean = numClean / (denSum + SMALL);
scalarField bulkNoisy = numNoisy / (denSum + SMALL);

// 3. Calculate 2nd-Order Finite Difference Derivatives
scalarField gradClean(nBins, 0.0);
scalarField gradNoisy(nBins, 0.0);

for(label i = 0; i < nBins; ++i) 
{
    // Cyclic indices
    label im1 = (i - 1 + nBins) % nBins; 
    label ip1 = (i + 1) % nBins;
    
    gradClean[i] = (bulkClean[ip1] - bulkClean[im1]) / (2.0 * dX);
    gradNoisy[i] = (bulkNoisy[ip1] - bulkNoisy[im1]) / (2.0 * dX);
}

// 4. Compute RMS Errors
scalar errCleanVal = 0.0, errNoisyVal = 0.0;
scalar errCleanGrad = 0.0, errNoisyGrad = 0.0;

for(label k = 0; k < nBins; ++k)
{
    scalar xNode = xMin + (k + 0.5) * dX;
    scalar exactFunc = Foam::sin(k_wave * (xNode - xMin)) + 2.0;
    scalar exactGrad = k_wave * Foam::cos(k_wave * (xNode - xMin));
    
    errCleanVal  += Foam::sqr(bulkClean[k] - exactFunc);
    errNoisyVal  += Foam::sqr(bulkNoisy[k] - exactFunc);
    
    errCleanGrad += Foam::sqr(gradClean[k] - exactGrad);
    errNoisyGrad += Foam::sqr(gradNoisy[k] - exactGrad);
}

Info << "1D Function Value RMS Error:" << endl;
Info << "  -> Clean Signal: " << Foam::sqrt(errCleanVal / nBins) << endl;
Info << "  -> Noisy Signal: " << Foam::sqrt(errNoisyVal / nBins) << endl;

Info << "1D Derivative (2nd Order) RMS Error:" << endl;
Info << "  -> Clean Signal: " << Foam::sqrt(errCleanGrad / nBins) << endl;
Info << "  -> Noisy Signal: " << Foam::sqrt(errNoisyGrad / nBins) << endl;

// 5. Dump Data for Plotting
OFstream file("fractionTest_nBins" + name(nBins) + ".dat");
file << "x_node ExactVal CleanVal NoisyVal ExactGrad CleanGrad NoisyGrad" << endl;
for(label k = 0; k < nBins; ++k)
{
    scalar xNode = xMin + (k + 0.5) * dX;
    scalar exactVal = Foam::sin(k_wave * (xNode - xMin)) + 2.0;
    scalar exactGrad = k_wave * Foam::cos(k_wave * (xNode - xMin));
    
    file << xNode << " " 
         << exactVal << " " << bulkClean[k] << " " << bulkNoisy[k] << " " 
         << exactGrad << " " << gradClean[k] << " " << gradNoisy[k] << endl;
}

Info << "--- END CELL-FRACTION STRESS TEST ---" << nl << endl;