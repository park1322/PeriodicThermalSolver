/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    bulkTemperatureFoam

Group
    JaewonParkCustomSolver

Description
    Periodic Constant Wall Temperature Nu/f Calculation; Robust For Any Type Of Mesh

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "simpleControl.H"
#include "fvOptions.H"

// For Periodic Indexing
inline label posMod(label a, label m)
{
    return (a % m + m) % m;
}

#include "initialConstructors.H"
// #include "calibrateKernelRadius.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Periodic Solver for CWT."
    );

    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createFields.H"
    #include "initContinuityErrs.H"

    turbulence->validate();

    const vector flowDir(1.0, 0.0, 0.0); // *Take this from fvOption in the future
    const scalarField& V = mesh.V();
    
    #include "createBB.H"

    scalarField lambda1D_debug;
    scalarField lambda1D_x_debug;
    
    BinProjectorX projector(mesh, nBins, xMin, Lperiod, dX);

    // unit check projector
    // {
    //     #include "checkBB_1.H"
    //     #include "checkBB_2.H"
    //     #include "checkBB_3.H"
    //     return 0;
    // }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting simpleFoam calculation loop\n" << endl;

    while (simple.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- Pressure-velocity SIMPLE corrector
        {
            #include "UEqn.H"
            #include "pEqn.H"
        }

        laminarTransport.correct();
        turbulence->correct();

        runTime.write();

        runTime.printExecutionTime(Info);
    }

    Info<< "\nStarting temperature calculation loop\n" << endl;


    const scalar tolerance = readScalar(
                                        mesh.solutionDict()
                                            .subDict("SIMPLE")
                                            .subDict("residualControl")
                                            .lookup("T")
                                        );
    // const int maxIter = readScalar(runTime.controlDict().lookup("endTime"));
    const int maxIter = runTime.timeIndex() + 120;
    // int iteration = runTime.timeIndex();
    scalar T_initialResidual(0.0);
    bool converged = false;
    while(!converged && runTime.timeIndex() < maxIter) 
    {
        // variables pre-compute
        volScalarField Ux = (U & flowDir);
        volScalarField UxAbs = mag(Ux);
        volScalarField UxT = UxAbs * T; 
        
        scalarField numSumBulkT(nBins, 0.0);
        scalarField denSumBulkT(nBins, 0.0);
        forAll(mesh.cells(), cellI)
        {
            // Generate signals
            scalar xCenter = mesh.C()[cellI].x();
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
                numSumBulkT[b] += UxT[cellI] * cellVolume;
                denSumBulkT[b] += UxAbs[cellI] * cellVolume;
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
                scalar fraction = overlapLength / L_cell;
                
                if (fraction > SMALL)
                {
                    label wrappedBin = posMod(b, nBins);
                    numSumBulkT[wrappedBin] += fraction * UxT[cellI] * cellVolume;
                    denSumBulkT[wrappedBin]   += fraction * UxAbs[cellI] * cellVolume;
                }
            }
        }
        reduce(numSumBulkT, sumOp<scalarField>());
        reduce(denSumBulkT, sumOp<scalarField>());
        scalarField bulkT = numSumBulkT / (denSumBulkT + SMALL);

        // use bulk mean temperature to compute normalization factor, Omega
        // tmp<volScalarField> tOmega = projector.projectSG(bulkT, "Omega", /*reciprocal=*/true);
        // volScalarField& Omega = tOmega.ref();

        volScalarField Omega = projector.projectConsistent
        (
            bulkT, 
            xMinField, xMaxField, L_cellField, 
            "lambda", 
            true // not reciprocal
        );

        T *= Omega;
        T.correctBoundaryConditions();

        // lambda solve
        scalarField P(nBins,0.0), Q(nBins,0.0), S(nBins,0.0);
        {
            scalarField P_num(nBins,0.0), Q_num(nBins,0.0), S_num(nBins,0.0), Vol_sum(nBins,0.0);
            const volScalarField thetax = (fvc::grad(T) & flowDir);
            volScalarField pp = fvc::div(phi, T) - fvc::laplacian(Dt, T) + Dt*T*sqr(lambda);
            volScalarField qq = 2*Dt*thetax - Ux*T + 2*Dt*T*lambda;
            volScalarField ss = Dt*T;

            forAll(mesh.cells(), cellI)
            {
                // Generate signals
                scalar xCenter = mesh.C()[cellI].x();
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
                    P_num[b]    += pp[cellI] * cellVolume;
                    Q_num[b]    += qq[cellI] * cellVolume;
                    S_num[b]    += ss[cellI] * cellVolume;
                    Vol_sum[b]  += cellVolume;
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
                    scalar fraction = overlapLength / L_cell;
                    
                    if (fraction > SMALL)
                    {
                        label wrappedBin = posMod(b, nBins);
                        P_num[wrappedBin]   += fraction * pp[cellI] * cellVolume;
                        Q_num[wrappedBin]   += fraction * qq[cellI] * cellVolume;
                        S_num[wrappedBin]   += fraction * ss[cellI] * cellVolume;
                        Vol_sum[wrappedBin]   += fraction * cellVolume;
                    }
                }
            }
            reduce(P_num, sumOp<scalarField>());
            reduce(Q_num, sumOp<scalarField>());
            reduce(S_num, sumOp<scalarField>());
            reduce(Vol_sum, sumOp<scalarField>());
            // P = P_num / (Vol_sum + SMALL);
            // Q = Q_num / (Vol_sum + SMALL);
            // S = S_num / (Vol_sum + SMALL);

            P = P_num;
            Q = Q_num;
            S = S_num;
            //////////////////////////////////////////////////////////////////////////////////////////////////////////////
            //////////// YOU MIGHT NOT NEED TO NORMALIZE THIS BY VOLUME //////////////////////////////////////////
            // turns out this does not matter much as it is like multiplying constant to both RHS and LHS.
            /////////////////////////////////////////////////////////////////////////////////////////////////////

        }

        scalarField Lambda1D(nBins, 0.0);
        {
            scalarField A(nBins), B(nBins), C(nBins), D(nBins);
            A = Q;
            B = (-0.5/dX) * S;
            C = (0.5/dX) * S;
            D = P;

            scalarField EE(nBins), FF(nBins), GG(nBins), pp(nBins), qq(nBins), rr(nBins);
            // first cell is 0; last cell is nBins-1. 
            EE[0] = B[0]/A[0];
            FF[0] = C[0]/A[0];
            GG[0] = D[0]/A[0];

            for (label k=1; k<nBins-1; ++k)
            {
                const scalar denom = A[k] - C[k] * EE[k-1];
                EE[k] = B[k] / denom;
                FF[k] = C[k] * FF[k-1] / denom;
                GG[k] = (D[k] + C[k] * GG[k-1]) / denom;
            }

            pp[0] = A[nBins-1];
            qq[0] = B[nBins-1];
            rr[0] = D[nBins-1];


            for (label k=1; k<nBins-1; ++k)
            {
                pp[k] = pp[k-1] - qq[k-1] * FF[k-1];
                qq[k] = qq[k-1] * EE[k-1];
                rr[k] = rr[k-1] + qq[k-1] * GG[k-1];
            }
            
            {
                const label N = nBins - 1; // last cell index
                const scalar denom = pp[N-1] - ( qq[N-1] + C[N] ) * ( EE[N-1] + FF[N-1]) ;
                const scalar numer = rr[N-1] + ( qq[N-1] + C[N] ) *  GG[N-1] ;

                Lambda1D[N] = numer / denom;

                // back substitution: k = N .. 0
                for (label k = nBins-2; k >= 0; --k)
                {
                    Lambda1D[k] = EE[k]*Lambda1D[k+1] + FF[k]*Lambda1D[N] + GG[k];
                    if (k == 0) break; // avoid label underflow
                }
            }

            // Project lambda
            // tmp<volScalarField> tLamTmp = projector.projectSG(Lambda1D, "Lambda_tmp");
            // lambda.primitiveFieldRef() = tLamTmp().internalField();   // updates internal cells
            // lambda.correctBoundaryConditions();                    // patches follow BCs

            tmp<volScalarField> tLamTmp = projector.projectConsistent
            (
                Lambda1D, 
                xMinField, xMaxField, L_cellField, 
                "lambda", 
                false // not reciprocal
            );
            lambda.primitiveFieldRef() = tLamTmp().internalField();
            lambda.correctBoundaryConditions();

            // Alternative lambda
            // tmp<volScalarField> tLamTmp = projector.project(Lambda1D, "Lambda_tmp");
            // const scalar relaxFactor = 0.3; 
            // volScalarField& lambdaRef = lambda.primitiveFieldRef();
            // lambdaRef = (1.0 - relaxFactor)*lambdaRef + relaxFactor*tLamTmp().primitiveFieldRef();
            // lambda.correctBoundaryConditions();

        }

        // tmp<volScalarField> tDxTmp = projector.projectDerivativeSG(Lambda1D, "lambda_x_tmp");
        // lambda_x.primitiveFieldRef() = tDxTmp().internalField();  // lambda_x already declared with 1/L^2
        // lambda_x.correctBoundaryConditions();

        scalarField dLambdaDx1D = projector.calcGradient2ndOrder(Lambda1D);
        tmp<volScalarField> tDxTmp = projector.projectConsistent
        (
            dLambdaDx1D, 
            xMinField, xMaxField, L_cellField, 
            "dLambdaDx", 
            false // not reciprocal
        );
        lambda_x.primitiveFieldRef() = tDxTmp().internalField();  // lambda_x already declared with 1/L^2
        lambda_x.correctBoundaryConditions();


        volScalarField sigma = ( 2.0*Dt*(fvc::grad(T) & flowDir) - (U & flowDir)*T )*lambda + Dt*T*( sqr(lambda) + lambda_x );
        // volScalarField sigma_exp = 2.0*Dt*(fvc::grad(T) & flowDir)*lambda   - (U & flowDir)*lambda*T + Dt*(lambda_x)*T;
        // volScalarField sigma_imp = Dt*sqr(lambda);

        fvScalarMatrix TEqn
        (
            fvm::div(phi, T) 
            - fvm::laplacian(Dt, T) 
            == 
            sigma
            // sigma_exp
            // + fvm::SuSp(sigma_imp, T)
        );

        // TEqn.relax();
        T_initialResidual = TEqn.solve().initialResidual();

        // Check convergence
        if (T_initialResidual < tolerance)
        {
            Info << nl << "Converged after " << runTime.timeIndex() + 1 << " iterations." << nl << endl;
            converged = true;
        }
        // ++iteration;
        ++runTime;
        Info << nl << "run time value " <<  runTime.timeIndex() << nl << endl;
        runTime.write();


        // Nu calculation method 1: direct wall flux

        label wallsID = mesh.boundaryMesh().findPatchID("primitive");
        const fvPatchScalarField& Tw   = T.boundaryField()[wallsID];
        const scalarField&       magSf = mesh.magSf().boundaryField()[wallsID];
        const scalar qhat = gSum( (-Tw.snGrad()) * magSf );
        const scalar Aw = gSum(magSf);
        const scalar Vtot = gSum(V);
        const scalar Nu_period = 4.0*Vtot/sqr(Aw) * qhat;
        // Info << nl << "Nu: " <<  Nu_period << nl << endl;


        // Nu calculation method 2: volume integral 
        sigma = ( 2.0*Dt*(fvc::grad(T) & flowDir) - (U & flowDir)*T )*lambda + Dt*T*( sqr(lambda) + lambda_x );
        const scalar volSigma = gSum( sigma.internalField() * V ); // qhat
        const scalar Nu_period_vol = 4.0*Vtot/sqr(Aw) * volSigma / Dt.value();

        Info << nl << "Nu (method 1): " <<  Nu_period << " | Nu (method 2): " << Nu_period_vol << nl << endl;
        lambda1D_debug = Lambda1D;
        lambda1D_x_debug =  projector.calcGradient2ndOrder(Lambda1D);
    }

    Info<< "DEBUG: lambda1D: " << lambda1D_debug << nl << endl; 
    Info<< "DEBUG: lambda1D_x: " << lambda1D_x_debug << nl << endl; 

    // Info<< "residual control T: " << tolerance << nl << endl;
    // Info<< "iteration:  " << iteration << nl << endl;    
    Info<< "End Time: " << maxIter << nl << endl;    
    Info<< "run Time: " << runTime.run() << nl << endl;    
    // Info << "BC types for lambda" <<lambda.boundaryField().types() << nl << endl;

    // runTime.write(); // does not work here. 
    runTime.writeNow();
    runTime.functionObjects().execute();
    runTime.printExecutionTime(Info);

    Info<< "End\n" << endl;
    return 0;
}


// ************************************************************************* //
