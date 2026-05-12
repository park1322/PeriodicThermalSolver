// ------------------------------------------------------------------------ //
//   Function: calibrateSmoothingLength
//   Purpose:  Finds optimal h by minimizing RMSE of a noisy sine wave reconstruction.
// ------------------------------------------------------------------------ //
scalar calibrateKernelRadius(const fvMesh& mesh, label nBins, const vector& flowDir, scalar xMin, scalar L)
{
    Info << nl << "--- Starting Kernel Calibration (Bias-Variance Optimization) ---" << endl;

    // 1. Mesh Statistics
    scalar totalVol = gSum(mesh.V());
    scalar globalNCells = returnReduce(mesh.nCells(), sumOp<label>());
    scalar avgCellVol = totalVol / globalNCells;
    scalar avgCellLen = Foam::pow(avgCellVol, 1.0/3.0);
    scalar dX_grid = L / scalar(nBins);

    // 2. Generate Synthetic Data (Noisy Sine Wave)
    // We use a high noise level (20%) to force the kernel to actually smooth.
    // If noise is 0, the optimizer would choose h -> 0 (Interpolation).
    const scalar k_wave = 2.0 * constant::mathematical::pi / L;
    const scalar noiseAmp = 0.2; // 20% Noise
    
    scalarField T_noisy(mesh.nCells());
    scalarField T_exact_3D(mesh.nCells());
    
    Random rndGen(1234); // Fixed seed for reproducibility

    forAll(mesh.C(), c)
    {
        scalar x = (mesh.C()[c] & flowDir) - xMin;
        scalar exact = Foam::sin(k_wave * x);
        
        // Add Noise [-0.2, 0.2]
        scalar noise = (rndGen.sample01<scalar>() - 0.5) * 2.0 * noiseAmp;
        
        T_noisy[c] = exact + noise;
        T_exact_3D[c] = exact; // Store for reference if needed
    }

    // 3. Optimization Loop (Grid Search)
    // Range: From 0.5x CellSize (Dangerous) to 4.0x CellSize (Blurry)
    
    scalar best_h = -1.0;
    scalar min_RMSE = GREAT;
    
    Info << "    " << "h_factor" << "    " << "h_value" << "    " << "RMSE" << endl;

    for (scalar factor = 0.5; factor <= 4.0; factor += 0.25)
    {
        scalar test_h = factor * avgCellLen;
        
        // Hard Constraint: h must be at least dX/2 to span the grid
        if (test_h < 0.5 * dX_grid) test_h = 0.5 * dX_grid;

        // --- A. Run Shepard Interpolation ---
        scalarField num(nBins, 0.0), den(nBins, 0.0);
        label kernelBinSpan = label(ceil(3.0 * test_h / dX_grid));

        forAll(mesh.C(), c)
        {
            scalar x = (mesh.C()[c] & flowDir) - xMin;
            label kCenter = label(x / dX_grid); // Simple floor for center
            
            // Adjust for periodicity in loop
            // (Simplified 1D scatter logic for speed)
            scalar val = T_noisy[c];
            scalar weight = mesh.V()[c];

            // Local coordinates
            scalar xNorm = (x / dX_grid) - 0.5;

            for (label offset = -kernelBinSpan; offset <= kernelBinSpan; ++offset)
            {
                // Unwrapped bin index
                // Note: kCenter depends on x, which is 0..L
                // We recreate precise logic inside:
                label kRaw = label(round(xNorm)) + offset;
                
                scalar distIdx = mag(scalar(kRaw) - xNorm);
                scalar r = distIdx * dX_grid;
                
                if (r >= 3.0 * test_h) continue;

                scalar W = 0.0;
                // Quintic Kernel Inline
                scalar q = r / test_h;
                if (q < 1.0) W = pow5(3.0-q) - 6.0*pow5(2.0-q) + 15.0*pow5(1.0-q);
                else if (q < 2.0) W = pow5(3.0-q) - 6.0*pow5(2.0-q);
                else if (q < 3.0) W = pow5(3.0-q);

                label k = posMod(kRaw, nBins);
                num[k] += W * val * weight;
                den[k] += W * weight;
            }
        }
        
        reduce(num, sumOp<scalarField>());
        reduce(den, sumOp<scalarField>());

        // --- B. Calculate RMSE against EXACT 1D Sine ---
        // We compare the reconstructed smooth curve against the TRUE clean sine wave.
        scalar errorSum = 0.0;
        for(label k=0; k<nBins; ++k)
        {
            scalar recon = num[k] / (den[k] + SMALL);
            
            scalar xNode = (k + 0.5) * dX_grid;
            scalar exact = Foam::sin(k_wave * xNode);
            
            errorSum += Foam::sqr(recon - exact);
        }
        
        scalar rmse = Foam::sqrt(errorSum / nBins);
        
        Info << "    " << factor << "    " << test_h << "    " << rmse;

        if (rmse < min_RMSE)
        {
            min_RMSE = rmse;
            best_h = test_h;
            Info << " *" << endl; // Mark current best
        }
        else
        {
            Info << endl;
        }
    }

    Info << "  -> Optimal h found: " << best_h << " (RMSE: " << min_RMSE << ")" << endl;

    // 4. Final Diagnostics: Neighbor Counting
    // Run one last pass with best_h to report stats to the user
    Info << "  -> Verifying Neighbor Statistics for h = " << best_h << "..." << endl;
    
    labelList neighbors(nBins, 0);
    label kernelBinSpan = label(ceil(3.0 * best_h / dX_grid));

    // Simple geometric counting (approximate but sufficient for stats)
    scalar cutoff = 3.0 * best_h;
    
    // Scan all cells again
    forAll(mesh.C(), c)
    {
        scalar x = (mesh.C()[c] & flowDir) - xMin;
        if(x < 0) x += L; if(x >= L) x -= L;

        label kCenter = label(x / dX_grid);
        // This cell contributes to which bins?
        // It's faster to loop bins near the cell
        label span = label(ceil(cutoff / dX_grid));
        
        for(label offset = -span; offset <= span; ++offset)
        {
            label kRaw = kCenter + offset;
            scalar binX = (kRaw + 0.5) * dX_grid; // unwrapped position
            
            // Distance with periodic wrap check
            // Distance between Cell X and Bin Center BinX
            // We align BinX to be close to X
            scalar dist = mag(x - binX);
            // Rough check logic for stats:
            // Actually, correct way:
            // Calculate distance in periodic domain
            label k = posMod(kRaw, nBins);
            scalar binX_wrapped = (k + 0.5)*dX_grid;
            scalar d = mag(x - binX_wrapped);
            if (d > 0.5*L) d = L - d;
            
            if (d < cutoff) {
                neighbors[k]++;
            }
        }
    }
    reduce(neighbors, sumOp<labelList>());

    label minN = gMax(neighbors); // Reduction across processors handled by reduce above? 
                                  // Wait, reduce(neighbors) sums the partial histograms. Correct.
                                  // Now find min/max in the list.
    minN = neighbors[0];
    label maxN = neighbors[0];
    scalar avgN = 0;
    
    for(label k=0; k<nBins; ++k) {
        if(neighbors[k] < minN) minN = neighbors[k];
        if(neighbors[k] > maxN) maxN = neighbors[k];
        avgN += neighbors[k];
    }
    avgN /= nBins;

    Info << "  -> Neighbor Statistics: Min=" << minN << ", Max=" << maxN << ", Avg=" << avgN << endl;
    
    if (minN < 10)
    {
        Info << "  WARNING: Optimal h results in very low neighbor counts (<10)." << endl;
        Info << "           This suggests the Mesh is extremely coarse or the Noise assumption is too low." << endl;
    }

    Info << "----------------------------------------------------------------" << nl << endl;

    return best_h;
}