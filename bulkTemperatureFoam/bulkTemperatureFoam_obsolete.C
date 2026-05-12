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

// #include "simpleControl.H"
// inline label posMod(label a, label m)
// {
//     // assumes m > 0
//     label r = a % m;
//     return (r < 0 ? r + m : r);
// }

inline label posMod(label a, label m)
{
    return (a % m + m) % m;
}

inline scalar pow5(scalar x) { scalar x2 = x*x; return x2*x2*x; }

#include "initialConstructors.H"

// #include "sampledPlane.H"      // from libsampling
// #include "surfaceInterpolate.H"
// #include "interpolationCellPoint.H"
// #include "interpolationCell.H"

// version 1
// inline scalar areaIntegralOnPlane
// (
//     const sampledPlane& cutplane,
//     const interpolation<scalar>& qip
// )
// {
//     auto tq = cutplane.sample(qip);                 
//     const Field<scalar>& q = tq();
//     const scalarField&   A = cutplane.magSf();
//     return gSum(q * A);
// }

// inline scalar areaIntegralOnPlane
// (
//     const sampledPlane& cutplane,
//     const volScalarField& q
// )
// {
//     interpolationCellPoint<scalar> qip(q);
//     auto tq = cutplane.sample(qip);                 
//     const Field<scalar>& qs = tq();
//     const scalarField&   A = cutplane.magSf();
//     return gSum(qs * A);
// }

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

    // #include "createBulkBins.H" obsolete 
    
    


    const vector flowDir(1.0, 0.0, 0.0); // *Take this from fvOption in the future
    const scalarField& V = mesh.V();
    const label kernelSize = 3; 

    // #include "createPlanes.H" obsolete

    #include "cellCenterNodeMap.H"
    BinProjectorX projector(mesh, nDiscrete, xMin, xMax);

    auto QuinticKernel = [&](scalar r) -> scalar
    {
        scalar q = r * invH;
        scalar val = 0.0;

        if (q < 1.0)
        {
            val = pow5(3.0 - q) - 6.0*pow5(2.0 - q) + 15.0*pow5(1.0 - q);
        }
        else if (q < 2.0)
        {
            val = pow5(3.0 - q) - 6.0*pow5(2.0 - q);
        }
        else if (q < 3.0)
        {
            val = pow5(3.0 - q);
        }
        
        return val;
    };




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
    const int maxIter = runTime.timeIndex() + 30;
    // int iteration = runTime.timeIndex();
    scalar T_initialResidual(0.0);
    bool converged = false;
    while(!converged && runTime.timeIndex() < maxIter) 
    {
        // variables pre-compute
        volScalarField Ux = (U & flowDir);
        volScalarField UxAbs = mag(Ux);
        volScalarField UxT = UxAbs * T; 

        // // Field Normalization 
        // // comput bulk mean T
        // scalarField Ux_int(nPlanes); 
        // scalarField UxT_int(nPlanes); 
        // {
        //     // interpolationCellPoint<scalar> Uxip(UxAbs);
        //     // interpolationCellPoint<scalar> UxTip(UxT);
        //     interpolationCell<scalar> Uxip(UxAbs);
        //     interpolationCell<scalar> UxTip(UxT);
        //     for (label i=0; i<nPlanes; ++i)
        //     {
        //         Ux_int[i] = areaIntegralOnPlane(planes[i], Uxip);
        //         UxT_int[i] = areaIntegralOnPlane(planes[i], UxTip);
        //     }
        // }
        // scalarField TbulkPlanes(UxT_int / Ux_int);
        // // Info<< "bulkT @ planes before linear normalization: " << TbulkPlanes << xPlanes << nl << endl;

        // scalarField TDiscrete(nDiscrete); // This is the bulk mean temperature
        // for (label k=0; k<nDiscrete; ++k)
        // {
        //     const scalar Ux_sum = Ux_int[k] + Ux_int[k+1];
        //     const scalar UxT_sum = UxT_int[k] + UxT_int[k+1];
        //     TDiscrete[k] = UxT_sum / Ux_sum;
        // }

        // // Info<< "bulkT @ discrete centers before linear normalization: " << TDiscrete << xDiscrete << nl << endl;

        forAll(mesh.C(), cellI)
        {
            label kCenter = cellCenterNode[cellI];
            scalar xRel = mesh.C()[cellI].x() - xMin;
            scalar bulkT_num = UxT[cellI] * V[cellI];
            scalar bulkT_den = UxAbs[cellI] * V[cellI];

            for (label offset = -kernelSize; offset <= kernelSize; ++offset) 
            {
                // 1. Recover neighbor index
                label k = kCenter + offset;
                
                // 2. Handle Cyclic Wrap (Fast bitwise or simple if)
                k = posMod(k, nBins);

                // 3. Compute Distance & Weight (On the fly)
                // We know node X position implicitly: (k + 0.5)*dX
                scalar distIdx = mag(scalar(k) - (xRel/dX - 0.5)); // Normalized distance
                scalar r = distIdx * dX;
                
                scalar W = QuinticKernel(r); 

                // 4. Scatter
                numSum[k_wrapped] += W * bulkT_num;
                denSum[k_wrapped] += W * bulkT_den;
            }
        }

        // 4. Parallel Reduction
        reduce(numSum, sumOp<scalarField>());
        reduce(denSum, sumOp<scalarField>());

        // 5. Finalize
        scalarField bulkT = numSum / (denSum + SMALL);

        
        // use bulk mean temperature to compute normalization factor, Omega
        tmp<volScalarField> tOmega = projector.project(TDiscrete, "Omega", /*reciprocal=*/true);
        volScalarField& Omega = tOmega.ref();
        T *= Omega;
        T.correctBoundaryConditions();


        // TEMP to check the normalization 
        // UxT = mag(U & flowDir) * T;              // assign into existing volScalarField
        
        // {
        //     interpolationCellPoint<scalar> Uxip(Ux);
        //     interpolationCellPoint<scalar> UxTip(UxT);

        //     for (label i=0; i<nPlanes; ++i)
        //     {
        //         Ux_int[i] = areaIntegralOnPlane(planes[i], Uxip);
        //         UxT_int[i] = areaIntegralOnPlane(planes[i], UxTip);
        //     }
        // }
        // TbulkPlanes = UxT_int / Ux_int;
        // // Info<< "bulkT @ planes after linear normalization: " << TbulkPlanes << xPlanes << nl << endl;

        // for (label k=0; k<nDiscrete; ++k)
        // {
        //     const scalar Ux_sum = Ux_int[k] + Ux_int[k+1];
        //     const scalar UxT_sum = UxT_int[k] + UxT_int[k+1];
        //     TDiscrete[k] = UxT_sum / Ux_sum;
        // }

        // Info<< "bulkT @ discrete centers after linear normalization: " << TDiscrete << xDiscrete << nl << endl;
        // TEMP END 

        // lambda solve
        scalarField P(nDiscrete,0.0), Q(nDiscrete,0.0), S(nDiscrete,0.0);
        {
            const volScalarField thetax = (fvc::grad(T) & flowDir);
            
            volScalarField pp = fvc::div(phi, T) - fvc::laplacian(Dt, T) + Dt*T*sqr(lambda);
            volScalarField qq = 2*Dt*thetax - Ux*T + 2*Dt*T*lambda;
            volScalarField ss = Dt*T;

            // compute plane integrals
            scalarField ppPlanes(nPlanes), qqPlanes(nPlanes), ssPlanes(nPlanes); 
            {
                // interpolationCellPoint<scalar> ppip(pp);
                // interpolationCellPoint<scalar> qqip(qq);
                // interpolationCellPoint<scalar> ssip(ss);
                interpolationCell<scalar> ppip(pp);
                interpolationCell<scalar> qqip(qq);
                interpolationCell<scalar> ssip(ss);
                for (label i=0; i<nPlanes; ++i)
                {
                    ppPlanes[i] = areaIntegralOnPlane(planes[i], ppip);
                    qqPlanes[i] = areaIntegralOnPlane(planes[i], qqip);
                    ssPlanes[i] = areaIntegralOnPlane(planes[i], ssip);
                }
            }


            for (label k=0; k<nDiscrete; ++k)
            {
                P[k] = ( ppPlanes[k] + ppPlanes[k+1] ) / ( areaPlanes[k] + areaPlanes[k+1] ); 
                Q[k] = ( qqPlanes[k] + qqPlanes[k+1] ) / ( areaPlanes[k] + areaPlanes[k+1] ); 
                S[k] = ( ssPlanes[k] + ssPlanes[k+1] ) / ( areaPlanes[k] + areaPlanes[k+1] ); 
            }

        }
    
        // Info<< "* * * * P * * * *: " << P << nl << endl;
        // Info<< "* * * * Q * * * *: " << Q << nl << endl;
        // Info<< "* * * * S * * * *: " << S << nl << endl;

        scalarField Lambda1D(nDiscrete, 0.0);
        {
            scalarField A(nDiscrete), B(nDiscrete), C(nDiscrete), D(nDiscrete);
            A = Q;
            B = (-0.5/dX) * S;
            C = (0.5/dX) * S;
            D = P;

            scalarField EE(nDiscrete), FF(nDiscrete), GG(nDiscrete), pp(nDiscrete), qq(nDiscrete), rr(nDiscrete);
            // first cell is 0; last cell is nDiscrete-1. 
            EE[0] = B[0]/A[0];
            FF[0] = C[0]/A[0];
            GG[0] = D[0]/A[0];

            for (label k=1; k<nDiscrete-1; ++k)
            {
                const scalar denom = A[k] - C[k] * EE[k-1];
                EE[k] = B[k] / denom;
                FF[k] = C[k] * FF[k-1] / denom;
                GG[k] = (D[k] + C[k] * GG[k-1]) / denom;
            }

            pp[0] = A[nDiscrete-1];
            qq[0] = B[nDiscrete-1];
            rr[0] = D[nDiscrete-1];


            for (label k=1; k<nDiscrete-1; ++k)
            {
                pp[k] = pp[k-1] - qq[k-1] * FF[k-1];
                qq[k] = qq[k-1] * EE[k-1];
                rr[k] = rr[k-1] + qq[k-1] * GG[k-1];
            }


            
            {
                const label N = nDiscrete - 1; // last cell index

                const scalar denom = pp[N-1] - ( qq[N-1] + C[N] ) * ( EE[N-1] + FF[N-1]) ;
                const scalar numer = rr[N-1] + ( qq[N-1] + C[N] ) *  GG[N-1] ;

                Lambda1D[N] = numer / denom;

                // back substitution: k = N .. 0
                for (label k = nDiscrete-2; k >= 0; --k)
                {
                    Lambda1D[k] = EE[k]*Lambda1D[k+1] + FF[k]*Lambda1D[N] + GG[k];
                    if (k == 0) break; // avoid label underflow
                }
            }

            // Info<< "* * * * lambda1D * * * *: " << Lambda1D << nl << endl;
            // Info<< "* * * * EE * * * *: " << EE << nl << endl;
            // Info<< "* * * * FF * * * *: " << FF << nl << endl;
            // Info<< "* * * * GG * * * *: " << GG << nl << endl;
            // Info<< "* * * * pp * * * *: " << pp << nl << endl;
            // Info<< "* * * * qq * * * *: " << qq << nl << endl;
            // Info<< "* * * * rr * * * *: " << rr << nl << endl;

            // Project lambda
            tmp<volScalarField> tLamTmp = projector.project(Lambda1D, "Lambda_tmp");
            lambda.primitiveFieldRef() = tLamTmp().internalField();   // updates internal cells
            lambda.correctBoundaryConditions();                    // patches follow BCs

        }

        // sigma calculation 
        // tmp<volScalarField> tDxTmp = projector.projectDx(Lambda1D, "lambda_x_tmp");
        // lambda_x.primitiveFieldRef() = tDxTmp().internalField();  // lambda_x already declared with 1/L^2
        // lambda_x.correctBoundaryConditions();

        scalarField dLamC(nDiscrete, 0.0);
        for (label k=0; k<nDiscrete; ++k)
        {
            const label km = (k-1+nDiscrete) % nDiscrete;
            const label kp = (k+1) % nDiscrete;
            dLamC[k] = (Lambda1D[kp] - Lambda1D[km]) / (2.0*dX);
        }

        tmp<volScalarField> tLamx = projector.project(dLamC, "lambda_x_from1D");
        lambda_x.primitiveFieldRef() = tLamx().internalField();
        // lambda_x.primitiveFieldRef() = (fvc::grad(lambda) & flowDir);
        lambda_x.correctBoundaryConditions();
        // volScalarField sigma = ( 2.0*Dt*(fvc::grad(T) & flowDir) - (U & flowDir)*T )*lambda + Dt*T*( sqr(lambda) + (fvc::grad(lambda) & flowDir) );
        volScalarField sigma = ( 2.0*Dt*(fvc::grad(T) & flowDir) - (U & flowDir)*T )*lambda + Dt*T*( sqr(lambda) + lambda_x );


        fvScalarMatrix TEqn
        (
            fvm::div(phi, T) 
            - fvm::laplacian(Dt, T) 
            == 
            sigma
        );

        // TEqn.relax();
        T_initialResidual = TEqn.solve().initialResidual();
        T.relax();
        T.storePrevIter();


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
    }

    // runTime.run() = true;

    // Info<< "residual control T: " << tolerance << nl << endl;
    // Info<< "iteration:  " << iteration << nl << endl;    
    Info<< "End Time: " << maxIter << nl << endl;    
    Info<< "run Time: " << runTime.run() << nl << endl;    
    // Info << "BC types for lambda" <<lambda.boundaryField().types() << nl << endl;

    // volScalarField Ux = (U & flowDir);
    // scalarField bulkU(nBins, 0);
    // for (label k=0; k<nBins; ++k)
    // {
    //     scalar bulkU_num(0.0); // T*abs(U)*V
    //     scalar bulkU_den(0.0); // abs(U)*V
    //     for (const label c : binCells[k])
    //     {
    //         const scalar magUx_V = mag(Ux[c]) * V[c];
    //         bulkU_num += magUx_V;
    //         bulkU_den += V[c];
    //     }

    //     bulkU_num = returnReduce(bulkU_num, sumOp<scalar>());
    //     bulkU_den = returnReduce(bulkU_den, sumOp<scalar>());
    //     bulkU[k] = bulkU_num / bulkU_den;
    // }  

    // tmp<volScalarField> tbulkU_interp = projector.project(bulkU, "tbulkU_interp");
    // volScalarField& bulkU_interp = tbulkU_interp.ref();
    
    // volScalarField checkBulkU
    // (
    //     IOobject
    //     (
    //         "checkBulkU",
    //         runTime.timeName(),
    //         mesh,
    //         IOobject::NO_READ,
    //         IOobject::AUTO_WRITE
    //     ),
    //     // mesh,
    //     bulkU_interp   // start as 0 everywhere
    // );

    // lambda.write();
    // lambda_x.write();
    // T.write();

    // runTime.write(); // does not work here. 
    runTime.writeNow();
    runTime.functionObjects().execute();

    volScalarField UxAbs = mag( U & flowDir);
    scalarField Ux_int(nPlanes); 
    {
        // interpolationCellPoint<scalar> Uxip(UxAbs);
        // interpolationCellPoint<scalar> UxTip(UxT);
        interpolationCell<scalar> Uxip(UxAbs);
        for (label i=0; i<nPlanes; ++i)
        {
            Ux_int[i] = areaIntegralOnPlane(planes[i], Uxip) / areaPlanes[i];
        }
    }

    Info<< "bulk U " << Ux_int << nl << endl;

    // Info<< "binCenters " << binCenters << nl << endl;
    Info<< "Time = " << simple.loop() << nl << endl;

    Info<< "End\n" << endl;


    
    
    return 0;

 
}


// ************************************************************************* //
