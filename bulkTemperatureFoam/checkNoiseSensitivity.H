// ========================================================================= //
//                  UNIT TEST: NOISE SENSITIVITY CHECK
// ========================================================================= //
Info << nl << "--- STARTING NOISE STRESS TEST (nBins = " << nBins << ") ---" << endl;

// 1. Setup: Sine Wave + White Noise
const scalar L = xMax - xMin;
const scalar k_wave = 2.0 * constant::mathematical::pi / L;
const scalar noiseAmp = 0.2; // 20% Noise (Aggressive test)

volScalarField T_noisy = T; // prepare as a copy. 

// Seed random number generator (simple pseudo-random)
Random rndGen(1234); 

forAll(mesh.C(), c)
{
    scalar x = mesh.C()[c].x();
    scalar exact = Foam::sin(k_wave * (x - xMin)) + 2.0;
    
    // Generate random number [0, 1] using the template method
    scalar randVal = rndGen.sample01<scalar>();
    // Info << "random #: " << randVal << endl;
    scalar noise = (randVal - 0.5) * 2.0 * noiseAmp; // Add random noise [-noiseAmp, +noiseAmp]
    
    T_noisy[c] = exact + noise;
}

// 2. Run Scatter (Kernel Averaging)
// This tests if the kernel acts as a good Low-Pass Filter
scalarField num(nBins, 0.0), den(nBins, 0.0);

forAll(mesh.C(), c)
{
    label kCenter = cellCenterNode[c];
    scalar xRel = mesh.C()[c].x() - xMin;
    
    scalar continuousIdx = (xRel / dX) - 0.5;
    scalar val = T_noisy[c];
    scalar weight = V[c]; 

    for (label offset = -kernelBinSpan; offset <= kernelBinSpan; ++offset) 
    {
        label kRaw = kCenter + offset; 
        scalar distIdx = mag(scalar(kRaw) - continuousIdx); 
        
        if (distIdx * dX >= 3.0*h) continue;

        scalar W = QuinticKernel(distIdx * dX); 
        label kStore = posMod(kRaw, nBins);

        num[kStore] += W * val * weight;
        den[kStore] += W * weight;
    }
}

reduce(num, sumOp<scalarField>());
reduce(den, sumOp<scalarField>());

scalarField bulkT_1D = num / (den + SMALL);

// 3. Check 1D Smoothness (Compare against Exact Sine)
// If RMS is low, Kernel did its job. If high, Kernel failed to smooth noise.
scalar errSum1D = 0.0;
for(label k=0; k<nBins; ++k)
{
    scalar xNode = xMin + (k + 0.5)*dX;
    scalar exact = Foam::sin(k_wave * (xNode - xMin)) + 2.0;
    scalar error = Foam::sqr(bulkT_1D[k] - exact);
    errSum1D += error;
}
Info << "1D Profile RMS Error (Smoothing Quality): " << Foam::sqrt(errSum1D/nBins) << endl;


// 4. Test Derivative Stiffness (4th Order vs 2nd Order)
scalarField grad4 = projector.calcGradient4thOrder(bulkT_1D);
scalarField SG = projector.calcGradientSG(bulkT_1D);
scalarField SG9 = projector.calcGradientSG9(bulkT_1D);

// Manually calc 2nd order for comparison
scalarField grad2(nBins);
for(label i=0; i<nBins; ++i) {
    label im1 = (i-1+nBins)%nBins; label ip1 = (i+1)%nBins;
    grad2[i] = (bulkT_1D[ip1] - bulkT_1D[im1]) / (2.0*dX);
}

// * * * New: Kernel-Based Consistent Derivative * * *
// We calculate N' and D' directly from the 3D noisy field
scalarField numGrad(nBins, 0.0), denGrad(nBins, 0.0);

forAll(mesh.C(), c)
{
    label kCenter = cellCenterNode[c];
    scalar xRel = mesh.C()[c].x() - xMin;
    scalar continuousIdx = (xRel / dX) - 0.5;
    
    scalar val = T_noisy[c]; // Use the mesh source values
    scalar weight = V[c]; 

    for (label offset = -kernelBinSpan; offset <= kernelBinSpan; ++offset) 
    {
        label kRaw = kCenter + offset; 
        
        // Signed distance for direction
        scalar distSigned = scalar(kRaw) - continuousIdx; 
        scalar distIdx = mag(distSigned);
        scalar r = distIdx * dX;
        
        if (r >= 3.0*h) continue;

        // Analytic Kernel Derivative wrt r
        scalar dWdr = QuinticKernelDeriv(r);
        
        // Spatial Derivative dW/dx = (dW/dr) * sign(node - particle)
        // sign is sign(distSigned)
        scalar sign = (distSigned >= 0) ? 1.0 : -1.0;
        scalar dWdx = dWdr * sign;

        label kStore = posMod(kRaw, nBins);

        // Accumulate Gradient Terms
        numGrad[kStore] += dWdx * val * weight;
        denGrad[kStore] += dWdx * weight;
    }
}

reduce(numGrad, sumOp<scalarField>());
reduce(denGrad, sumOp<scalarField>());

// Construct the Consistent Derivative using Quotient Rule
// f' = (N'D - ND') / D^2
// Note: We use 'num' and 'den' from Step 2
scalarField gradKernel(nBins, 0.0);
for(label k=0; k<nBins; ++k)
{
    scalar D = den[k] + SMALL;
    scalar N = num[k];
    scalar N_prime = numGrad[k];
    scalar D_prime = denGrad[k];
    
    gradKernel[k] = (N_prime * D - N * D_prime) / (D*D);
}


// * * * Compare Errors * * *
scalar errGrad4 = 0.0;
scalar errGrad2 = 0.0;
scalar errSG = 0.0;
scalar errKernel = 0.0;
scalar errSG9 = 0.0;


// Check against EXACT derivative: k * cos(k*x)
for(label k=0; k<nBins; ++k)
{
    scalar xNode = xMin + (k + 0.5)*dX;
    scalar exactSlope = k_wave * Foam::cos(k_wave * (xNode - xMin));
    
    errGrad4  += Foam::sqr(grad4[k] - exactSlope);
    errGrad2  += Foam::sqr(grad2[k] - exactSlope);
    errSG     += Foam::sqr(SG[k] - exactSlope);
    errKernel += Foam::sqr(gradKernel[k] - exactSlope);
    errSG9    += Foam::sqr(SG9[k] - exactSlope);
}

Info << "Derivative RMS Error (Stiffness Check):" << endl;
Info << "  -> 4th Order:       " << Foam::sqrt(errGrad4/nBins) << endl;
Info << "  -> 2nd Order:       " << Foam::sqrt(errGrad2/nBins) << endl;
Info << "  -> Savitzky-Golay:  " << Foam::sqrt(errSG/nBins) << endl;
Info << "  -> Kernel-Based:    " << Foam::sqrt(errKernel/nBins) << endl;
Info << "  -> Savitzky-Golay - 9:  " << Foam::sqrt(errSG9/nBins) << endl;

if (Foam::sqrt(errGrad4/nBins) > 2.0 * Foam::sqrt(errGrad2/nBins))
{
    Info << "  WARNING: 4th Order is significantly amplifying noise!" << endl;
}


// 5. Test 3D Value Projection (Overfitting Check)
totalVol = gSum(mesh.V());

// 5a. 4th Order Projection
tmp<volScalarField> tProj4 = projector.project(bulkT_1D, "Test_Proj4");
scalar errSum3D_4 = 0.0;
forAll(mesh.C(), c) {
    scalar x = mesh.C()[c].x();
    scalar exact = Foam::sin(k_wave * (x - xMin)) + 2.0;
    errSum3D_4 += Foam::sqr(tProj4()[c] - exact) * V[c];
}

// 5b. SG Projection
tmp<volScalarField> tProjSG = projector.projectSG(bulkT_1D, "Test_ProjSG");
scalar errSum3D_SG = 0.0;
forAll(mesh.C(), c) {
    scalar x = mesh.C()[c].x();
    scalar exact = Foam::sin(k_wave * (x - xMin)) + 2.0;
    errSum3D_SG += Foam::sqr(tProjSG()[c] - exact) * V[c];
}

Info << "3D Value Projection RMS Error:" << endl;
Info << "  -> 4th Order: " << Foam::sqrt(errSum3D_4/totalVol) << endl;
Info << "  -> SG Method: " << Foam::sqrt(errSum3D_SG/totalVol) << endl;


// 6. Test 3D Derivative Projection (CRITICAL CHECK)
// This checks if the projected gradient field (dT/dx) is smooth or noisy.

// 6a. 4th Order Derivative Projection
tmp<volScalarField> tDeriv4 = projector.projectDerivative(bulkT_1D, "Test_Deriv4");
scalar errDeriv3D_4 = 0.0;
forAll(mesh.C(), c) {
    scalar x = mesh.C()[c].x();
    scalar exactSlope = k_wave * Foam::cos(k_wave * (x - xMin));
    errDeriv3D_4 += Foam::sqr(tDeriv4()[c] - exactSlope) * V[c];
}

// 6b. SG Derivative Projection
tmp<volScalarField> tDerivSG = projector.projectDerivativeSG(bulkT_1D, "Test_DerivSG");
scalar errDeriv3D_SG = 0.0;
forAll(mesh.C(), c) {
    scalar x = mesh.C()[c].x();
    scalar exactSlope = k_wave * Foam::cos(k_wave * (x - xMin));
    errDeriv3D_SG += Foam::sqr(tDerivSG()[c] - exactSlope) * V[c];
}

Info << "3D Derivative Projection RMS Error:" << endl;
Info << "  -> 4th Order: " << Foam::sqrt(errDeriv3D_4/totalVol) << endl;
Info << "  -> SG Method: " << Foam::sqrt(errDeriv3D_SG/totalVol) << endl;


// 7. Dump Data for Plotting
OFstream file("stressTest_nBins" + name(nBins) + ".dat");
file << "# x_node  Exact  Bulk1D  SlopeExact  Slope4th  Slope2nd SlopeSG SlopeKernel SlopeSG9" << endl;
for(label k=0; k<nBins; ++k)
{
    scalar xNode = xMin + (k + 0.5)*dX;
    scalar exactVal = Foam::sin(k_wave * (xNode - xMin)) + 2.0;
    scalar exactSlope = k_wave * Foam::cos(k_wave * (xNode - xMin));
    
    file << xNode << " " << exactVal << " " << bulkT_1D[k] 
         << " " << exactSlope << " " << grad4[k] << " " << grad2[k] << " " << SG[k] << " " << gradKernel[k] << " " << SG9[k] << endl;
}

Info << "--- END STRESS TEST ---" << nl << endl;
