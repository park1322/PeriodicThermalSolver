struct BinProjectorX
{
    const fvMesh&        mesh;
    const label          nBins;
    const scalar         xMin;
    const scalar         Lperiod;
    const scalar         dX;

    // Precomputed per-cell mapping
    labelList            kLeft;   // Index of the node to the left
    scalarField          tNorm;   // normalized position 
    
    BinProjectorX(const fvMesh& m, label nBins_, scalar xMin_, scalar L_, scalar dX_)
    : mesh(m)
    , nBins(nBins_)
    , xMin(xMin_)
    , Lperiod(L_)
    , dX(dX_)
    , kLeft(m.nCells(), 0)
    , tNorm(m.nCells(), 0)
    {
        const vectorField& cellCenters = mesh.C().internalField();

        forAll(cellCenters, cellI)
        {
            // 1. Calculate periodic relative position
            scalar s = cellCenters[cellI].x() - xMin;

            // 2. Identify the node to the "Left"
            // Nodes are at (k + 0.5)*dX. 
            // We want index k such that node[k] <= s < node[k+1]
            // This effectively shifts the coordinate system by 0.5*dX
            scalar s_shifted = s - 0.5*dX;
            if (s_shifted < 0) s_shifted += Lperiod; // Wrap for the segment before node 0

            label k = label(floor(s_shifted / dX));
            k = k % nBins; // Safety
            
            kLeft[cellI] = k;

            // 3. Local coordinate t in [0, 1] between Node k and Node k+1
            // Refine local s relative to the start of the interval
            scalar s_local = s_shifted - k*dX; 
            
            // Normalize to t [0..1]
            scalar t = s_local / dX;
            tNorm[cellI] = max(0.0, min(1.0, t)); // Clamp guard

        }
    }

    // * * * Helper: 2nd Order Central Difference * * *
    // We use this to get the slope at the nodes
    scalarField calcGradient2ndOrder(const scalarField& y) const
    {
        scalarField dy(nBins);
        const scalar inv2dx = 1.0 / (2.0 * dX);

        for (label i = 0; i < nBins; ++i)
        {
            label im1 = (i - 1 + nBins) % nBins;
            label ip1 = (i + 1) % nBins;
            dy[i] = (y[ip1] - y[im1]) * inv2dx;
        }
        return dy;
    }

    // * * * Helper: Calculate 4th Order Cyclic Gradient * * *
    // Returns dy/dx at the nodes
    scalarField calcGradient4thOrder(const scalarField& y) const
    {
        scalarField dy(nBins);
        const scalar inv12dx = 1.0 / (12.0 * dX);

        for (label i = 0; i < nBins; ++i)
        {
            // Cyclic Indices
            label im2 = (i - 2 + nBins) % nBins;
            label im1 = (i - 1 + nBins) % nBins;
            label ip1 = (i + 1) % nBins;
            label ip2 = (i + 2) % nBins;

            // 4th Order Central Difference
            // f' = (-f(i+2) + 8f(i+1) - 8f(i-1) + f(i-2)) / 12h
            scalar term = -y[ip2] + 8.0*y[ip1] - 8.0*y[im1] + y[im2];
            dy[i] = term * inv12dx;
        }
        return dy;
    }


    // * * * New Method: Savitzky-Golay Gradient (5-point, Quadratic) * * *
    // Robust against noise. Calculates slope of the best-fit parabola.
    scalarField calcGradientSG(const scalarField& y) const
    {
        scalarField dy(nBins);
        // Factor is 1 / (10 * dX) for 5-point quadratic SG
        const scalar inv10dx = 1.0 / (10.0 * dX); 

        for (label i = 0; i < nBins; ++i)
        {
            // Cyclic Indices (5-point window)
            label im2 = (i - 2 + nBins) % nBins;
            label im1 = (i - 1 + nBins) % nBins;
            label ip1 = (i + 1) % nBins;
            label ip2 = (i + 2) % nBins;

            // SG Formula: (-2*y[i-2] - 1*y[i-1] + 1*y[i+1] + 2*y[i+2]) / 10dx
            // Note: y[i] (center) has weight 0 for the derivative
            scalar term = -2.0*y[im2] - 1.0*y[im1] + 1.0*y[ip1] + 2.0*y[ip2];
            
            dy[i] = term * inv10dx;
        }
        return dy;
    }

    // * * * 9-Point Savitzky-Golay Gradient * * *
    // "Wide" smoothing window. Great for suppressing noise without 
    // flattening peaks as much as a Gaussian filter.
    scalarField calcGradientSG9(const scalarField& y) const
    {
        scalarField dy(nBins);
        
        // Normalization factor for 9-point quadratic SG is 1 / (60 * dx)
        const scalar inv60dx = 1.0 / (60.0 * dX);

        for (label i = 0; i < nBins; ++i)
        {
            // Cyclic Indices
            label im4 = (i - 4 + nBins) % nBins;
            label im3 = (i - 3 + nBins) % nBins;
            label im2 = (i - 2 + nBins) % nBins;
            label im1 = (i - 1 + nBins) % nBins;
            
            label ip1 = (i + 1) % nBins;
            label ip2 = (i + 2) % nBins;
            label ip3 = (i + 3) % nBins;
            label ip4 = (i + 4) % nBins;

            // Convolution Weights: [-4, -3, -2, -1, 0, 1, 2, 3, 4]
            scalar term = -4.0*y[im4] - 3.0*y[im3] - 2.0*y[im2] - 1.0*y[im1]
                          +1.0*y[ip1] + 2.0*y[ip2] + 3.0*y[ip3] + 4.0*y[ip4];
            
            dy[i] = term * inv60dx;
        }
        return dy;
    }


    // * * * SAFE Linear Projection (Values) * * *
    // No derivatives involved. Strictly bounded by neighbors.
    tmp<volScalarField> project
    (
        const scalarField& y1D,      
        const word&        outName,
        bool               reciprocal=false
    ) const
    {
        tmp<volScalarField> tY
        (
            new volScalarField
            (
                IOobject(outName, mesh.time().timeName(), mesh,
                         IOobject::NO_READ, IOobject::NO_WRITE),
                mesh,
                dimensionedScalar("", dimless, 0.0)
            )
        );
        volScalarField& Y = tY.ref();

        forAll(Y, c)
        {
            const label k  = kLeft[c];
            const label kp = (k + 1) % nBins; 

            const scalar t = tNorm[c]; 
            
            // Linear Interpolation: y = (1-t)*y_k + t*y_{k+1}
            scalar val = (1.0 - t)*y1D[k] + t*y1D[kp];

            if (reciprocal)
            {
                // Safety clamp for reciprocal
                // if (mag(val) < SMALL) val = (val >= 0 ? SMALL : -SMALL);
                // Y[c] = 1.0 / val;
                // Y[c] = reciprocal ? (val > SMALL ? 1.0/val : 1.0) : val;
                scalar sign = (val >= 0.0) ? 1.0 : -1.0;
                Y[c] = (mag(val) > SMALL) ? (1.0 / val) : (sign / SMALL);
            }
            else
            {
                Y[c] = val;
            }
        }

        Y.correctBoundaryConditions();
        return tY;
    }

    // * * * SAFE Linear Projection (Derivatives) * * *
    // Calculates dy/dx at nodes, then interpolates that slope linearly.
    tmp<volScalarField> projectDerivative
    (
        const scalarField& y1D, 
        const word& outName
    ) const
    {
        // 1. Calculate derivatives on the 1D grid first
        scalarField dy1D = calcGradient2ndOrder(y1D);
        // scalarField dy1D = calcGradientSG(y1D);
        
        tmp<volScalarField> tY
        (
            new volScalarField
            (
                IOobject(outName, mesh.time().timeName(), mesh,
                         IOobject::NO_READ, IOobject::NO_WRITE),
                mesh,
                dimensionedScalar("", dimless, 0.0)
            )
        );
        volScalarField& Y = tY.ref();

        forAll(Y, c)
        {
            const label k  = kLeft[c];
            const label kp = (k + 1) % nBins; 
            const scalar t = tNorm[c]; 

            // Interpolate the derivative value linearly
            // This gives a C0 continuous derivative field (Piecewise linear slope)
            Y[c] = (1.0 - t)*dy1D[k] + t*dy1D[kp];
        }

        Y.correctBoundaryConditions();
        return tY;
    }

/////////////////////////////////////////////////////// Below is experimental //////////////////////

    // // * * * New Project Method * * *
    // // Requires both Value (y1D) and Derivative (dy1D) arrays from 1D grid
    // tmp<volScalarField> project
    // (
    //     const scalarField& y1D,      // Value at nodes
    //     const word&        outName,
    //     bool                reciprocal=false
    // ) const
    // {
    //     // Calculate gradients locally
    //     // scalarField dy1D = calcGradient4thOrder(y1D);
    //     scalarField dy1D = calcGradientSG(y1D);

    //     tmp<volScalarField> tY
    //     (
    //         new volScalarField
    //         (
    //             IOobject(outName, mesh.time().timeName(), mesh,
    //                      IOobject::NO_READ, IOobject::NO_WRITE),
    //             mesh,
    //             dimensionedScalar("", dimless, 0.0)
    //         )
    //     );
    //     volScalarField& Y = tY.ref();

    //     forAll(Y, c)
    //     {
    //         const label k  = kLeft[c];
    //         const label kp = (k + 1) % nBins; // Right neighbor (cyclic)

    //         // Cubic Hermite Reconstruction
    //         // P(t) = h00*p0 + h10*m0*dX + h01*p1 + h11*m1*dX
    //         // Note: dy1D must be multiplied by dX to match units

    //         const scalar t = tNorm[c]; 
    //         const scalar t2 = t*t;
    //         scalar t3 = t2*t;

    //         scalar h00 = 2.0*t3 - 3.0*t2 + 1.0;      // Weight for y[k]
    //         scalar h10 = t3 - 2.0*t2 + t;            // Weight for y'[k] * dX
    //         scalar h01 = -2.0*t3 + 3.0*t2;           // Weight for y[k+1]
    //         scalar h11 = t3 - t2;   
            
    //         scalar y = h00 * y1D[k] +
    //                     h10 * dy1D[k] * dX +
    //                     h01 * y1D[kp] +
    //                     h11 * dy1D[kp] * dX;
    //         Y[c] = reciprocal ? (y > SMALL ? 1.0/y : 1.0) : y;
    //     }

    //     Y.correctBoundaryConditions();
    //     return tY;
    // }

    // // * * * New Method: Project Exact Spatial Derivative * * *
    // // Returns d(Value)/dx at the cell center analytically
    // tmp<volScalarField> projectDerivative
    // (
    //     const scalarField& y1D, 
    //     const word& outName
    // ) const
    // {
    //     // 1. Get nodal derivatives (dy1D) using FD once (Safe)
    //     // scalarField dy1D = calcGradient4thOrder(y1D);
    //     scalarField dy1D = calcGradientSG(y1D);

    //     tmp<volScalarField> tY
    //     (
    //         new volScalarField
    //         (
    //             IOobject(outName, mesh.time().timeName(), mesh,
    //                      IOobject::NO_READ, IOobject::NO_WRITE),
    //             mesh,
    //             dimensionedScalar("", dimless, 0.0)
    //         )
    //     );
    //     volScalarField& Y = tY.ref();

    //     forAll(Y, c)
    //     {
    //         const label k  = kLeft[c];
    //         const label kp = (k + 1) % nBins; 
            
    //         // We need 't' for the derivative basis.
    //         // If you removed 'tFrac' to save memory, re-calculate it here.
    //         // Otherwise, reuse the stored tFrac[c].
    //         const scalar t = tNorm[c]; 
    //         const scalar t2 = t*t;

    //         // Analytical Derivatives of Basis Functions with respect to t (dh/dt)
    //         // h00 = 2t^3 - 3t^2 + 1   -> h00' = 6t^2 - 6t
    //         // h10 = t^3 - 2t^2 + t    -> h10' = 3t^2 - 4t + 1
    //         // h01 = -2t^3 + 3t^2      -> h01' = -6t^2 + 6t
    //         // h11 = t^3 - t^2         -> h11' = 3t^2 - 2t

    //         scalar dh00 = 6.0*t2 - 6.0*t;
    //         scalar dh10 = 3.0*t2 - 4.0*t + 1.0;
    //         scalar dh01 = -6.0*t2 + 6.0*t;
    //         scalar dh11 = 3.0*t2 - 2.0*t;

    //         // Chain Rule: df/dx = (df/dt) * (dt/dx)
    //         // dt/dx = 1 / dX
    //         // The dX factor in the h10/h11 terms cancels out partially
            
    //         // Reconstruct df/dt
    //         scalar df_dt = dh00 * y1D[k] +
    //                        dh10 * dy1D[k] * dX +
    //                        dh01 * y1D[kp] +
    //                        dh11 * dy1D[kp] * dX;

    //         // Convert to physical gradient df/dx
    //         Y[c] = df_dt / dX;
    //     }

    //     Y.correctBoundaryConditions();
    //     return tY;
    // }



    // ////////////////////////////////////////////// TEMP /////////////////////////////////
    // * * * New Project Method * * *
    // Requires both Value (y1D) and Derivative (dy1D) arrays from 1D grid
    tmp<volScalarField> projectSG
    (
        const scalarField& y1D,      // Value at nodes
        const word&        outName,
        bool                reciprocal=false
    ) const
    {
        // Calculate gradients locally
        scalarField dy1D = calcGradientSG(y1D);

        tmp<volScalarField> tY
        (
            new volScalarField
            (
                IOobject(outName, mesh.time().timeName(), mesh,
                         IOobject::NO_READ, IOobject::NO_WRITE),
                mesh,
                dimensionedScalar("", dimless, 0.0)
            )
        );
        volScalarField& Y = tY.ref();

        forAll(Y, c)
        {
            const label k  = kLeft[c];
            const label kp = (k + 1) % nBins; // Right neighbor (cyclic)

            // Cubic Hermite Reconstruction
            // P(t) = h00*p0 + h10*m0*dX + h01*p1 + h11*m1*dX
            // Note: dy1D must be multiplied by dX to match units

            const scalar t = tNorm[c]; 
            const scalar t2 = t*t;
            scalar t3 = t2*t;

            scalar h00 = 2.0*t3 - 3.0*t2 + 1.0;      // Weight for y[k]
            scalar h10 = t3 - 2.0*t2 + t;            // Weight for y'[k] * dX
            scalar h01 = -2.0*t3 + 3.0*t2;           // Weight for y[k+1]
            scalar h11 = t3 - t2;   
            
            scalar y = h00 * y1D[k] +
                        h10 * dy1D[k] * dX +
                        h01 * y1D[kp] +
                        h11 * dy1D[kp] * dX;

            if (reciprocal)
            {
                // Safety clamp for reciprocal
                scalar sign = (y >= 0.0) ? 1.0 : -1.0;
                Y[c] = (mag(y) > SMALL) ? (1.0 / y) : (sign / SMALL);
            }
            else
            {
                Y[c] = y;
            }
        }

        Y.correctBoundaryConditions();
        return tY;
    }

    // * * * New Method: Project Exact Spatial Derivative * * *
    // Returns d(Value)/dx at the cell center analytically
    tmp<volScalarField> projectDerivativeSG
    (
        const scalarField& y1D, 
        const word& outName
    ) const
    {
        // 1. Get nodal derivatives (dy1D) using FD once (Safe)
        scalarField dy1D = calcGradientSG(y1D);
        
        tmp<volScalarField> tY
        (
            new volScalarField
            (
                IOobject(outName, mesh.time().timeName(), mesh,
                         IOobject::NO_READ, IOobject::NO_WRITE),
                mesh,
                dimensionedScalar("", dimless, 0.0)
            )
        );
        volScalarField& Y = tY.ref();

        forAll(Y, c)
        {
            const label k  = kLeft[c];
            const label kp = (k + 1) % nBins; 
            
            // We need 't' for the derivative basis.
            // If you removed 'tFrac' to save memory, re-calculate it here.
            // Otherwise, reuse the stored tFrac[c].
            const scalar t = tNorm[c]; 
            const scalar t2 = t*t;

            // Analytical Derivatives of Basis Functions with respect to t (dh/dt)
            // h00 = 2t^3 - 3t^2 + 1   -> h00' = 6t^2 - 6t
            // h10 = t^3 - 2t^2 + t    -> h10' = 3t^2 - 4t + 1
            // h01 = -2t^3 + 3t^2      -> h01' = -6t^2 + 6t
            // h11 = t^3 - t^2         -> h11' = 3t^2 - 2t

            scalar dh00 = 6.0*t2 - 6.0*t;
            scalar dh10 = 3.0*t2 - 4.0*t + 1.0;
            scalar dh01 = -6.0*t2 + 6.0*t;
            scalar dh11 = 3.0*t2 - 2.0*t;

            // Chain Rule: df/dx = (df/dt) * (dt/dx)
            // dt/dx = 1 / dX
            // The dX factor in the h10/h11 terms cancels out partially
            
            // Reconstruct df/dt
            scalar df_dt = dh00 * y1D[k] +
                           dh10 * dy1D[k] * dX +
                           dh01 * y1D[kp] +
                           dh11 * dy1D[kp] * dX;

            // Convert to physical gradient df/dx
            Y[c] = df_dt / dX;
        }

        Y.correctBoundaryConditions();
        return tY;
    }

    // * * * NEW: Consistent Volume-Fraction Projection * * *
    // This is the exact mathematical adjoint of the 3D->1D binning process.
    tmp<volScalarField> projectConsistent
    (
        const scalarField& y1D, 
        const scalarField& xMinField,
        const scalarField& xMaxField,
        const scalarField& L_cellField,
        const word& outName,
        bool reciprocal = false
    ) const
    {
        tmp<volScalarField> tY
        (
            new volScalarField
            (
                IOobject(outName, mesh.time().timeName(), mesh,
                         IOobject::NO_READ, IOobject::NO_WRITE),
                mesh,
                dimensionedScalar("", dimless, 0.0)
            )
        );
        volScalarField& Y = tY.ref();

        forAll(mesh.cells(), cellI)
        {
            scalar xMinCell = xMinField[cellI];
            scalar xMaxCell = xMaxField[cellI];
            scalar L_cell   = L_cellField[cellI];

            scalar cellValue = 0.0;

            if (L_cell < SMALL) 
            {
                // Fallback for flat boundary cells
                scalar xCenter = mesh.C()[cellI].x();
                label b = floor((xCenter - xMin) / dX);
                cellValue = y1D[posMod(b, nBins)];
            }
            else
            {
                label startBin = floor((xMinCell - xMin) / dX);
                label endBin   = floor((xMaxCell - xMin) / dX);

                scalar totalFraction = 0.0;
                scalar accumulatedValue = 0.0;

                for (label b = startBin; b <= endBin; ++b) 
                {
                    scalar binStart = xMin + b * dX;
                    scalar binEnd   = binStart + dX;
                    
                    scalar overlapStart = max(xMinCell, binStart);
                    scalar overlapEnd   = min(xMaxCell, binEnd);
                    scalar overlapLength = max(0.0, overlapEnd - overlapStart);
                    scalar fraction = overlapLength / L_cell;

                    if (fraction > SMALL)
                    {
                        label wrappedBin = posMod(b, nBins);
                        accumulatedValue += fraction * y1D[wrappedBin];
                        totalFraction += fraction;
                    }
                }
                // Normalize by the actual fraction captured (should be ~1.0)
                cellValue = accumulatedValue / (totalFraction + SMALL);
            }

            // Safe Reciprocal Assignment (preserves negative signs)
            if (reciprocal)
            {
                scalar sign = (cellValue >= 0.0) ? 1.0 : -1.0;
                Y[cellI] = (mag(cellValue) > SMALL) ? (1.0 / cellValue) : (sign / SMALL);
            }
            else
            {
                Y[cellI] = cellValue;
            }
        }

        Y.correctBoundaryConditions();
        return tY;
    }
};

