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
#include "calibrateKernelRadius.H"

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
    
    
    #include "cellCenterNodeMap.H"

    Info << "DEBUG_CHECK: " 
     << " xMin=" << xMin 
     << " L=" << Lperiod 
     << " nCells=" << nCells 
     << " totalVol=" << totalVol 
     << " avgCellLen=" << avgCellLen 
     << " h=" << h 
     << endl;

    scalarField lambda1D_debug;
    scalarField lambda1D_x_debug;
    
    BinProjectorX projector(mesh, nBins, xMin, Lperiod, dX);

    auto QuinticKernel = [&](scalar r) -> scalar
    {
        scalar q = r * invH;
        
        // Early exit: Optimization for non-qualified neighbors
        if (q >= 3.0) return 0.0;

        scalar val = 0.0;

        if (q < 1.0)        val = pow5(3.0 - q) - 6.0*pow5(2.0 - q) + 15.0*pow5(1.0 - q);
        else if (q < 2.0)   val = pow5(3.0 - q) - 6.0*pow5(2.0 - q);
        else if (q < 3.0)   val = pow5(3.0 - q);
        
        return val;
    };

    auto QuinticKernelDeriv = [&](scalar r) -> scalar
    {
        scalar q = r * invH;
        if (q >= 3.0) return 0.0;
        
        scalar val = 0.0;
        // Derivative of (A-q)^5 is -5*(A-q)^4
        if (q < 1.0)      val = -5.0*pow4(3.0 - q) + 30.0*pow4(2.0 - q) - 75.0*pow4(1.0 - q);
        else if (q < 2.0) val = -5.0*pow4(3.0 - q) + 30.0*pow4(2.0 - q);
        else if (q < 3.0) val = -5.0*pow4(3.0 - q);
        
        return val * invH; // Chain rule dq/dr
    };


    // // unit check projector
    // {
        // #include "checkProjector.H"
    //     #include "checkKernel.H"
            // #include "checkNoiseSensitivity.H"
            // return 0;
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
        
        scalarField numSum(nBins, 0);
        scalarField denSum(nBins, 0);
        forAll(mesh.C(), cellI)
        {
            label kCenter = cellCenterNode[cellI];
            scalar xRel = mesh.C()[cellI].x() - xMin;
            scalar xNorm = (xRel / dX) - 0.5; // normalized; continuous index
            scalar bulkT_num = UxT[cellI] * V[cellI];
            scalar bulkT_den = UxAbs[cellI] * V[cellI];

            for (label offset = -kernelBinSpan; offset <= kernelBinSpan; ++offset) 
            {
                // 1. Recover neighbor index
                label k = kCenter + offset;

                // 3. Compute Distance & Weight (On the fly)
                // We know node X position implicitly: (k + 0.5)*dX
                scalar distIdx = mag(scalar(k) - xNorm); // Normalized distance
                scalar r = distIdx * dX;
                
                if (r >= 3.0*h) continue;
                scalar W = QuinticKernel(r); 

                // 2. Handle Cyclic Wrap (Fast bitwise or simple if)
                k = posMod(k, nBins);

                // 4. Scatter
                numSum[k] += W * bulkT_num;
                denSum[k] += W * bulkT_den;
            }
        }

        // 4. Parallel Reduction
        reduce(numSum, sumOp<scalarField>());
        reduce(denSum, sumOp<scalarField>());

        // 5. Finalize
        scalarField bulkT = numSum / (denSum + SMALL);

        // use bulk mean temperature to compute normalization factor, Omega
        tmp<volScalarField> tOmega = projector.project(bulkT, "Omega", /*reciprocal=*/true);
        volScalarField& Omega = tOmega.ref();
        const scalar omegaRelax = 0.25; 
        Omega.primitiveFieldRef() = 1.0 + omegaRelax * (Omega.primitiveFieldRef() - 1.0);
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

            forAll(mesh.C(), cellI)
            {
                label kCenter = cellCenterNode[cellI];
                scalar xRel = mesh.C()[cellI].x() - xMin;
                scalar xNorm = (xRel / dX) - 0.5;
                
                scalar vol = V[cellI];
                scalar val_pp = pp[cellI] * vol;
                scalar val_qq = qq[cellI] * vol;
                scalar val_ss = ss[cellI] * vol;

                for (label offset = -kernelBinSpan; offset <= kernelBinSpan; ++offset) 
                {
                    label k = kCenter + offset; 
                    scalar distIdx = mag(scalar(k) - xNorm); 
                    scalar r = distIdx * dX;
                    
                    // Optimization: Early exit
                    if (r >= 3.0*h) continue; 
                    scalar W = QuinticKernel(r); 

                    // Wrap index for storage
                    k = posMod(k, nBins);

                    P_num[k]   += W * val_pp;
                    Q_num[k]   += W * val_qq;
                    S_num[k]   += W * val_ss;
                    Vol_sum[k] += W * vol;
                }
            }

            reduce(P_num, sumOp<scalarField>());
            reduce(Q_num, sumOp<scalarField>());
            reduce(S_num, sumOp<scalarField>());
            reduce(Vol_sum, sumOp<scalarField>());

            P = P_num / (Vol_sum + SMALL);
            Q = Q_num / (Vol_sum + SMALL);
            S = S_num / (Vol_sum + SMALL);
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
            tmp<volScalarField> tLamTmp = projector.project(Lambda1D, "Lambda_tmp");
            lambda.primitiveFieldRef() = tLamTmp().internalField();   // updates internal cells
            lambda.correctBoundaryConditions();                    // patches follow BCs

            // Alternative lambda
            // tmp<volScalarField> tLamTmp = projector.project(Lambda1D, "Lambda_tmp");
            // const scalar relaxFactor = 0.3; 
            // volScalarField& lambdaRef = lambda.primitiveFieldRef();
            // lambdaRef = (1.0 - relaxFactor)*lambdaRef + relaxFactor*tLamTmp().primitiveFieldRef();
            // lambda.correctBoundaryConditions();

        }

        tmp<volScalarField> tDxTmp = projector.projectDerivative(Lambda1D, "lambda_x_tmp");
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
    // runTime.run() = true;

    // Info<< "residual control T: " << tolerance << nl << endl;
    // Info<< "iteration:  " << iteration << nl << endl;    
    Info<< "End Time: " << maxIter << nl << endl;    
    Info<< "run Time: " << runTime.run() << nl << endl;    
    // Info << "BC types for lambda" <<lambda.boundaryField().types() << nl << endl;

    volScalarField Ux = (U & flowDir);
    scalarField numSum(nBins, 0);
    scalarField denSum(nBins, 0);
    forAll(mesh.C(), cellI)
    {
        label kCenter = cellCenterNode[cellI];
        scalar xRel = mesh.C()[cellI].x() - xMin;
        scalar xNorm = (xRel / dX) - 0.5; // normalized; continuous index
        scalar bulkU_num = Ux[cellI] * V[cellI];
        scalar bulkU_den = V[cellI];

        for (label offset = -kernelBinSpan; offset <= kernelBinSpan; ++offset) 
        {
            // 1. Recover neighbor index
            label k = kCenter + offset;

            // 3. Compute Distance & Weight (On the fly)
            // We know node X position implicitly: (k + 0.5)*dX
            scalar distIdx = mag(scalar(k) - xNorm); // Normalized distance
            scalar r = distIdx * dX;
            
            if (r >= 3.0*h) continue;
            scalar W = QuinticKernel(r); 

            // 2. Handle Cyclic Wrap (Fast bitwise or simple if)
            k = posMod(k, nBins);

            // 4. Scatter
            numSum[k] += W * bulkU_num;
            denSum[k] += W * bulkU_den;
        }
    }

    // 4. Parallel Reduction
    reduce(numSum, sumOp<scalarField>());
    reduce(denSum, sumOp<scalarField>());

    // 5. Finalize
    scalarField bulkU = numSum / (denSum + SMALL);

    tmp<volScalarField> tbulkU_interp = projector.project(bulkU, "tbulkU_interp");
    volScalarField& bulkU_interp = tbulkU_interp.ref();
    
    volScalarField checkBulkU
    (
        IOobject
        (
            "checkBulkU",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        // mesh,
        bulkU_interp   // start as 0 everywhere
    );

    // lambda.write();
    // lambda_x.write();
    // T.write();

    // runTime.write(); // does not work here. 
    runTime.writeNow();
    runTime.functionObjects().execute();
    runTime.printExecutionTime(Info);

    // volScalarField UxAbs = mag( U & flowDir);

    // Info<< "binCenters " << binCenters << nl << endl;
    // Info<< "Time = " << simple.loop() << nl << endl;

    Info<< "End\n" << endl;


    
    
    return 0;

 
}


// ************************************************************************* //
