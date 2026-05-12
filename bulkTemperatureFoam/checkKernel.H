// ========================================================================= //
//                          UNIT TEST: SINE WAVE CHECK
// ========================================================================= //
Info << nl << "--- STARTING UNIT TEST ---" << endl;

// 1. Define Analytical Function params
const scalar L = xMax - xMin;
const scalar k_wave = 2.0 * constant::mathematical::pi / L; // One full sine wave

// 2. Assign Analytical Profile to 3D Mesh
// We use T_test instead of T to not mess up your simulation restart
volScalarField T_test = T; 

forAll(mesh.C(), cellI)
{
    scalar x = mesh.C()[cellI].x();
    // Function: sin(kx) + 2.0 (Shifted up to be positive)
    T_test[cellI] = Foam::sin(k_wave * (x - xMin)) + 2.0;
}

// 3. Run the Scatter Logic (Copy-Paste of your kernel loop)
scalarField testNum(nBins, 0.0);
scalarField testDen(nBins, 0.0);

forAll(mesh.C(), cellI)
{
    label kCenter = cellCenterNode[cellI];
    scalar xRel = mesh.C()[cellI].x() - xMin;
    
    // // Robust Periodic Wrap
    // if (xRel < 0) xRel += L;
    // if (xRel >= L) xRel -= L;

    scalar continuousIdx = (xRel / dX) - 0.5;
    
    // * * * TEST WEIGHTING * * * // We use Volume weighting (V) only. 
    // Effectively assuming uniform density and uniform velocity U=1.
    scalar val = T_test[cellI];
    scalar weight = V[cellI];  

    for (label offset = -kernelSize; offset <= kernelSize; ++offset) 
    {
        label kRaw = kCenter + offset; 
        scalar distIdx = mag(scalar(kRaw) - continuousIdx); 
        
        if (distIdx * dX >= 3.0*h) continue;

        scalar W = QuinticKernel(distIdx * dX); 
        label kStore = posMod(kRaw, nBins);

        testNum[kStore] += W * val * weight;
        testDen[kStore] += W * weight;
    }
}

// Parallel Reduce
reduce(testNum, sumOp<scalarField>());
reduce(testDen, sumOp<scalarField>());

// 4. Compute 1D Result and Error
scalarField result1D = testNum / (testDen + SMALL);
scalar totalError = 0.0;

// Create a file to dump results for plotting
OFstream file("unitTest_profile.dat");
file << "# x_node   T_reconstructed   T_analytical   Error" << endl;

for(label k=0; k<nBins; ++k)
{
    // Node center position
    scalar xNode = xMin + (k + 0.5)*dX;
    
    // Analytical value at the node center
    scalar exact = Foam::sin(k_wave * (xNode - xMin)) + 2.0;
    
    scalar error = mag(result1D[k] - exact);
    totalError += error * error; // Sum squared error

    file << xNode << " " << result1D[k] << " " << exact << " " << error << endl;
}

scalar rmsError = Foam::sqrt(totalError / nBins);

Info << "Unit Test RMS Error: " << rmsError << endl;

if (rmsError < 1e-2) 
{
    Info << "-> PASSED: Reconstruction is accurate." << endl;
}
else 
{
    Info << "-> FAILED: Error is too high! Check kernel radius or geometry." << endl;
    // FatalErrorInFunction << "Stopping for debugging." << exit(FatalError);
}

Info << "--- END UNIT TEST ---" << nl << endl;

// ========================================================================= //



// ========================================================================= //
//                          UNIT TEST: STEP FUNCTION (GIBBS CHECK)
// ========================================================================= //
Info << nl << "--- STARTING STEP FUNCTION TEST ---" << endl;

const scalar L = xMax - xMin;
volScalarField T_test = T; 

// 1. Assign Step Function
// Left half = 3.0, Right half = 1.0
forAll(mesh.C(), cellI)
{
    scalar x = mesh.C()[cellI].x();
    scalar xRel = x - xMin;

    if (xRel < 0.5 * L) T_test[cellI] = 3.0;
    else                T_test[cellI] = 1.0;
}

// 2. Run Scatter (Same as before)
scalarField testNum(nBins, 0.0);
scalarField testDen(nBins, 0.0);

forAll(mesh.C(), cellI)
{
    label kCenter = cellCenterNode[cellI];
    scalar xRel = mesh.C()[cellI].x() - xMin;
    
    scalar continuousIdx = (xRel / dX) - 0.5;
    
    scalar val = T_test[cellI];
    scalar weight = V[cellI];  

    for (label offset = -kernelSize; offset <= kernelSize; ++offset) 
    {
        label kRaw = kCenter + offset; 
        scalar distIdx = mag(scalar(kRaw) - continuousIdx); 
        
        if (distIdx * dX >= 3.0*h) continue;

        scalar W = QuinticKernel(distIdx * dX); 
        label kStore = posMod(kRaw, nBins);

        testNum[kStore] += W * val * weight;
        testDen[kStore] += W * weight;
    }
}

reduce(testNum, sumOp<scalarField>());
reduce(testDen, sumOp<scalarField>());

scalarField result1D = testNum / (testDen + SMALL);

// 3. Check for Overshoot/Undershoot (Ringing)
scalar maxVal = gMax(result1D);
scalar minVal = gMin(result1D);

Info << "Step Test Results:" << endl;
Info << "  Target Max: 3.0 | Reconstructed Max: " << maxVal << endl;
Info << "  Target Min: 1.0 | Reconstructed Min: " << minVal << endl;

// Allow a tiny float tolerance (e.g. 3.000001 is fine)
if (maxVal > 3.001 || minVal < 0.999)
{
        Info << "-> WARNING: Ringing detected! The kernel is oscillating." << endl;
}
else
{
        Info << "-> PASSED: Monotonic reconstruction. No ringing." << endl;
}

// Dump for plotting
OFstream file("unitTest_step.dat");
file << "# x_node   T_reconstructed   T_analytical" << endl;
for(label k=0; k<nBins; ++k)
{
    scalar xNode = xMin + (k + 0.5)*dX;
    scalar exact = (xNode - xMin < 0.5*L) ? 3.0 : 1.0;
    file << xNode << " " << result1D[k] << " " << exact << endl;
}

Info << "--- END STEP TEST ---" << nl << endl;
