// ========================================================================= //
//                     UNIT TEST: PROJECTOR & DERIVATIVES
// ========================================================================= //
if (true)
{
    Info << nl << "--- STARTING PROJECTOR VALIDATION ---" << endl;

    // 1. Setup Analytical Function (Sine Wave)
    // f(x) = sin(k x)
    const scalar L = xMax - xMin;
    const scalar k_wave = 2.0 * constant::mathematical::pi / L;
    
    scalarField y1D(nBins);
    scalarField dy1D_exact(nBins);

    for(label i=0; i<nBins; ++i)
    {
        scalar xNode = xMin + (i + 0.5)*dX; // Node center
        y1D[i] = Foam::sin(k_wave * (xNode - xMin));
        dy1D_exact[i] = k_wave * Foam::cos(k_wave * (xNode - xMin));
    }

    // --------------------------------------------------------------------- //
    // Test A: 1D Derivative Accuracy (4th Order FD)
    // --------------------------------------------------------------------- //
    Info << "Test A: Checking 4th Order Finite Difference..." << endl;
    
    // Use the projector's internal helper
    scalarField dy1D_num = projector.calcGradient4thOrder(y1D);

    scalar errorSumDy = 0.0;
    scalar maxErrorDy = 0.0;

    OFstream file("unitTest_derivative_sine_1D.dat");
    file << "# x_node   derivative_num   derivative_analytic error" << endl;

    for(label i=0; i<nBins; ++i)
    {
        scalar xNode = xMin + (i + 0.5)*dX;
        scalar err = mag(dy1D_num[i] - dy1D_exact[i]);
        file << xNode << " " << dy1D_num[i] << " " << dy1D_exact[i] << " " << err << endl;
        errorSumDy += err * err;
        maxErrorDy = max(maxErrorDy, err);
    }
    scalar rmseDy = Foam::sqrt(errorSumDy / nBins);

    Info << "  Derivative RMSE: " << rmseDy << " | Max Error: " << maxErrorDy << endl;
    
    if (rmseDy < 1e-3) Info << "  -> PASSED: Derivative calculation is accurate." << endl;
    else               Info << "  -> WARNING: Derivative error seems high." << endl;


    // --------------------------------------------------------------------- //
    // Test B: 3D Hermite Projection Accuracy
    // --------------------------------------------------------------------- //
    Info << "Test B: Checking 3D Hermite Reconstruction..." << endl;

    // We project the 1D sine wave onto the 3D mesh.
    // Ideally, T_3D[cell] should equal sin(k * x_cell).
    // The projector calculates dy1D_num internally to do this.
    
    tmp<volScalarField> tProj = projector.project(y1D, "test_Hermite");
    const volScalarField& Y_3D = tProj();

    scalar errorSumProj = 0.0;
    scalar maxErrorProj = 0.0;
    scalar totalVol = 0.0;

    forAll(mesh.C(), c)
    {
        scalar x = mesh.C()[c].x();
        scalar exact = Foam::sin(k_wave * (x - xMin));
        
        // Error at this cell
        scalar err = mag(Y_3D[c] - exact);
        
        // Volume-weighted RMSE is more physically meaningful
        errorSumProj += Foam::sqr(err) * V[c];
        maxErrorProj = max(maxErrorProj, err);
        totalVol += V[c];
    }
    
    scalar rmseProj = Foam::sqrt(errorSumProj / totalVol);

    Info << "  Projection RMSE: " << rmseProj << " | Max Error: " << maxErrorProj << endl;

    if (rmseProj < 1e-3) Info << "  -> PASSED: Projection is accurate." << endl;
    else                 Info << "  -> WARNING: Projection error seems high." << endl;

    Info << "--- END VALIDATION ---" << nl << endl;

    // --------------------------------------------------------------------- //
    // Test C: 3D Analytical Derivative Projection Accuracy
    // --------------------------------------------------------------------- //
    Info << "Test C: Checking 3D Analytical Derivative Projection..." << endl;

    // We project the derivative of the 1D sine wave onto the 3D mesh analytically.
    // Ideally, dY_3D[cell] should equal k_wave * cos(k_wave * x_cell).
    
    tmp<volScalarField> tProjDeriv = projector.projectDerivative(y1D, "test_Hermite_Deriv");
    const volScalarField& dY_3D = tProjDeriv();

    scalar errorSumDeriv = 0.0;
    scalar maxErrorDeriv = 0.0;
    totalVol = 0.0; // Reset for this loop
    
    forAll(mesh.C(), c)
    {
        scalar x = mesh.C()[c].x();
        
        // Exact Derivative: d/dx( sin(k(x-x0)) ) = k * cos(k(x-x0))
        scalar exact = k_wave * Foam::cos(k_wave * (x - xMin));
        
        // Error at this cell
        scalar err = mag(dY_3D[c] - exact);
        
        // Volume-weighted RMSE
        errorSumDeriv += Foam::sqr(err) * V[c];
        maxErrorDeriv = max(maxErrorDeriv, err);
        totalVol += V[c];
    }
    
    scalar rmseDeriv = Foam::sqrt(errorSumDeriv / totalVol);

    Info << "  Derivative Projection RMSE: " << rmseDeriv << " | Max Error: " << maxErrorDeriv << endl;

    // We allow a slightly looser tolerance for derivatives as they amplify curvature errors,
    // but for a sine wave and Hermite splines, it should still be excellent (< 1e-3).
    if (rmseDeriv < 5e-3) Info << "  -> PASSED: Derivative projection is accurate." << endl;
    else                  Info << "  -> WARNING: Derivative projection error seems high." << endl;
}