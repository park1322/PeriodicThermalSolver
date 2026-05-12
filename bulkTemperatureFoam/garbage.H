// // calculate bulk lambda
    // label inletID = mesh.boundaryMesh().findPatchID("minX");
    // label wallsID = mesh.boundaryMesh().findPatchID("primitive");
    // label outletID = mesh.boundaryMesh().findPatchID("maxX");
    // const vector flowDir(1.0, 0.0, 0.0);

    // Info<< "inletID = " << inletID << nl << endl;
    // Info<< "wallsID = " << wallsID << nl << endl;

    // {
    //     dimensionedScalar lambdaNum
    //     (
    //         "lambdaNum",
    //         dimensionSet(0, 3, -1, 1, 0, 0, 0),
    //         gSum(Dt.value()*-T.boundaryField()[wallsID].snGrad()*mesh.boundary()[wallsID].magSf())
    //     );

    //     dimensionedScalar lambdaDen
    //     (
    //         "lambdaDen",
    //         dimensionSet(0, 3, -1, 1, 0, 0, 0),
    //         gSum(T.boundaryField()[inletID]*(flowDir & U.boundaryField()[inletID])*mesh.boundary()[inletID].magSf()) + 
    //         gSum(Dt.value()*(flowDir & fvc::grad(T)().boundaryField()[inletID])*mesh.boundary()[inletID].magSf())
    //     );

    //     Info << "Ux*Theta: " << gSum((T.boundaryField()[inletID]*(flowDir & U.boundaryField()[inletID])*mesh.boundary()[inletID].magSf())) << endl;
    //     Info << "alpha*dTheta/dx: " << gSum(Dt.value()*(flowDir & fvc::grad(T)().boundaryField()[inletID])*mesh.boundary()[inletID].magSf()) << endl;

    //     Info << "lambda_num: " << lambdaNum.value() << endl;
    //     Info << "lambda_den: " << lambdaDen.value() << endl;
    //     Info << "lambda_num/den: " << lambdaNum.value()/lambdaDen.value() << endl;

    //     dimensionedScalar lambda = -(1/Lperiod)*Foam::log(1 - (lambdaNum/lambdaDen));
    //     Info << "lambda: " << lambda.value() << endl;

    //     auto inletEnthalpy = gSum((T.boundaryField()[inletID]*(phi.boundaryField()[inletID])));
    //     auto inletConduction = gSum(Dt.value()*-T.boundaryField()[inletID].snGrad()*mesh.boundary()[inletID].magSf());
    //     auto outletEnthalpy = gSum((T.boundaryField()[outletID]*(phi.boundaryField()[outletID])));
    //     auto wallflux = gSum(Dt.value()*-T.boundaryField()[wallsID].snGrad()*mesh.boundary()[wallsID].magSf());
    //     auto residual = inletEnthalpy + inletConduction + wallflux + outletEnthalpy; 

    //     Info << "inletEnthalpy: " << inletEnthalpy << endl;
    //     Info << "inletConduction: " << inletConduction << endl;
    //     Info << "outletEnthalpy: " << outletEnthalpy << endl;
    //     Info << "wallflux: " << wallflux << endl;
    //     Info << "residual: " << residual << endl;
    // }

// interpolation original version
volScalarField Omega = T*0.0;
forAll(Omega, cellI)
{
    const thisBin = binOfCell[cellI]; // this cell is in this bin
    const scalar x = cellCenters[cellI].x(); // cell x location
    const scalar bin_x = binCenters[thisBin]; // this bin center x 
    
    if (x < bin_x)
    {
        if (thisBin == 0)
        {
            leftBin = nBins - 1;
        } else 
        {
            leftBin = thisBin - 1;
        }

        Omega[cellI] = bulkT[leftBin] + (bulkT[thisBin]-bulkT[leftBin])/dX * (dX - bin_x + x);
    }
    else if (x > bin_x)
    {
        if (thisBin ==  nBins - 1)
        {
            rightBin = 0;
        } else 
        {
            rightBin = thisBin + 1;
        }

        Omega[cellI] = bulkT[thisBin] + (bulkT[rightBin]-bulkT[thisBin])/dX * (x - bin_x);
    }
    else 
    {
        Omega[cellI] = bulkT[thisBin];
    }
}


// interpolation improved version
volScalarField Omega = T*0.0;
forAll(Omega, cellI)
{
    const label thisBin = binOfCell[cellI]; // this cell is in this bin
    const scalar x = cellCenters[cellI].x(); // cell x location
    const scalar bin_x = binCenters[thisBin]; // this bin center x 
    
    // choose neighbors
    label leftBin, rightBin;
    if (x < bin_x)
    {
        leftBin  = posMod(thisBin - 1, nBins);
        rightBin = thisBin;
    }
    else if (x > bin_x)
    {
        leftBin  = thisBin;
        rightBin = posMod(thisBin + 1, nBins);
    }
    else
    {
        // exactly at center (rare): no interpolation needed
        const scalar y = bulkT[thisBin];
        Omega[cellI] = (y > SMALL ? 1.0/y : 1.0);
        continue;
    }

    // --- unwrap local coordinate across the period if we cross the end ---
    // We want a local line segment from left center -> right center, with the cell x in-between.
    scalar xLeft  = binCenters[leftBin];
    scalar xRight = binCenters[rightBin];
    scalar xLoc   = x;

    // Case: wrap from last bin to first bin (xMax -> xMin)
    // left ~= xMax - dX/2, right ~= xMin + dX/2, and x is near xMin
    if (leftBin == nBins - 1 && rightBin == 0)
    {
        xRight += Lperiod;   // push right to the next period
        if (xLoc < xLeft) xLoc += Lperiod;  // push cell too, so xLeft <= xLoc <= xRight
    }
    // (Symmetric case isn't needed because we always pick consecutive bins around x)

    // linear fraction 0..1 from left center toward right center
    scalar t = (xLoc - xLeft) / (xRight - xLeft);
    t = max(0.0, min(1.0, t));

    // interpolate bulkT and set Omega = 1 / interp
    const scalar y = (1.0 - t)*bulkT[leftBin] + t*bulkT[rightBin];
    Omega[cellI] = (y > SMALL ? 1.0/y : 1.0);

}


 // // --- jump to the latest available time, put this right after create time
    // {
    //     const instantList& ts = runTime.times();   // sorted list of time dirs (e.g. 0, 10, 20, ...)
    //     if (!ts.empty())
    //     {
    //         const label lastI = ts.size() - 1;
    //         runTime.setTime(ts[lastI], lastI);
    //         Info<< "Starting from latest time = " << runTime.timeName() << nl << endl;
    //     }
    //     else
    //     {
    //         FatalErrorInFunction << "No time directories found." << exit(FatalError);
    //     }
    // }


    // Project to a temporary (dims don’t matter)
            tmp<volScalarField> tLamTmp = projector.project(Lambda1D, "Lambda_tmp");
            const volScalarField& lamTmp = tLamTmp();

            // Copy only the numeric values into the pre-declared lambda
            lambda.primitiveFieldRef() = lamTmp.internalField();   // updates internal cells
            lambda.correctBoundaryConditions();                    // patches follow BCs





// cut plane stuff

// const vector nPlane(1,0,0);
    // scalar plane_x = 0.49;            // plane location
    // point  p0(plane_x, 0, 0);

   
    // plane pl	
    // (	
    //     p0,
    //     nPlane
    // );	

    // sampledPlane cutplane
    // (
    //     "cutplane",
    //     mesh,
    //     pl
    // );

    // cutplane.update();
    // const scalarField& A_cut = cutplane.magSf();
    // const scalar A_cut_total = gSum(A_cut);
    // Info << nl << "cut plane magSf: " <<  A_cut << nl << endl;
    // Info << nl << "cut plane total area: " <<  A_cut_total << nl << endl;

    // // Set the pointers
    // interpolationCellPoint<vector> Uip(U);
    // interpolationCellPoint<scalar> Tip(T);

    // // const Field<vector>& Utri = cutplane.sample(Uip)();
    // // const Field<scalar>& Ttri = cutplane.sample(Tip)();

    // auto tU = cutplane.sample(Uip);          // tmp<Field<vector>>
    // auto tT = cutplane.sample(Tip);          // tmp<Field<scalar>>
    // const Field<vector>& Utri = tU();        // valid while tU is in scope
    // const Field<scalar>& Ttri = tT();        // valid while tT is in scope

    // // scalar U_num = gSum( mag(Utri & nPlane) * A_cut );
    // scalar U_bulk = gSum( mag(Utri & nPlane) * A_cut ) / A_cut_total;
    // Info << nl << "U_bulk: " <<  U_bulk << nl << endl;

    // scalar T_bulk = gSum( Ttri * mag(Utri & nPlane) * A_cut ) / U_bulk;

    // Info << nl << "U_bulk: " <<  U_bulk << "    T_bulk: " <<  T_bulk << endl;

        const label target = 30;
        const vectorField& C  = mesh.C().internalField();
        const scalarField& Om = Omega.internalField();

        const label N = C.size();
        const label stride = round(static_cast<scalar>(N) / target);
        const label nSamp  = (N - 1)/stride + 1;   // exact number we’ll take

        scalarField xSample(nSamp, 0.0);
        scalarField omegaSample(nSamp, 0.0);
        label counter = 0;
        for (label c = 0; c < N; c += stride)
        {
            xSample[counter]     = C[c].x();
            omegaSample[counter] = Om[c];
            ++counter;
        }

        volScalarField OmegaCopy(Omega);   // copies data, safe but extra memory
        OmegaCopy.writeOpt() = IOobject::AUTO_WRITE;
        OmegaCopy.write();
        Info<< "bulkT @ interpolated Tbulk: " << omegaSample << xSample << nl << endl;