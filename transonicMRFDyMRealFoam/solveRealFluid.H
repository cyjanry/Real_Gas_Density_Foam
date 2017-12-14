{
        // update the boundary absolute velocity
        mrfZones.correctBoundaryVelocity(U);

        // set the pseudo face velocity
        Godunov.dotX() = mrfZones.faceU();

        // activate sub-cycle loop
        TimeState subCycleTimeState = runTime.subCycle(numberSubCycles);

        // get access to physical deltaT in ddt scheme!
        physDeltaT[1] = subCycleTimeState.deltaT().value();
        physDeltaT[2] = subCycleTimeState.deltaT0().value();

        // Begin sub-cycling - PseudoTime Integration
        // adjustTimeStep and numberSubCycles > 1 does not make any sense
        for (label j=0; j < numberSubCycles; j++)
        {
            // Update subCycle Time
            runTime++;

            // update local time step sizes just once each iteration for
            // all ddt schemes
            if (!adjustTimeStep)
            {  
                localTimeStep_Real.update(maxCo,adjustTimeStep);
            }

            // Plot Pseudo Time here, so that one can analyse the residuls with
            // pyFoam*
            Info<< "\n Time = " << runTime.value() << nl << endl;

            // Begin outer Loop for Runge - Kutta
            forAll (beta,i)
            {
                p.storePrevIter();
                // Update primitive boundaries
                p.correctBoundaryConditions();
                U.correctBoundaryConditions();
                h.correctBoundaryConditions();

                //    Info << "                       For cell 42468, h=" << h.internalField()[42468] << " p=" << p.internalField()[42468] 
                //     << " T=" << T.internalField()[42468] 
                 //    << " rho=" << rho.internalField()[42468]   << " e="  << e.internalField()[42468] << endl;

                //Info << " - -  -Start - - - - - - -- - " << nl << "compare rho with thermo.rho()" << nl << 
                //"Boundary: rho:" << rho.boundaryField()[5][50] << " thermo.rho():" << thermo.rho()().boundaryField()[5][50] 
                //<< " T:"<< T.boundaryField()[5][50]<< "T.thermo" << thermo.T().boundaryField()[5][50] 
                //<< " p:"<< p.boundaryField()[5][50] << "thermo.p" << thermo.p().boundaryField()[5][50]
                //<< " h:"<< h.boundaryField()[5][50] << "thermo.h" << thermo.h().boundaryField()[5][50] << nl
                //<< "Internal: rho:" << rho.internalField()[24630] << " and thermo.rho:" << thermo.rho()().internalField()[24630] 
                //<< " T:"<< T.internalField()[24630] << " Thermo.T" << thermo.T().internalField()[24630] 
                //<< " p:"<< p.internalField()[24630] << " Thermo.p" << thermo.p().internalField()[24630]
                //<< " h:"<< h.internalField()[24630] << "thermo.h" << thermo.h().internalField()[24630]  <<  endl;


                
                //Info << "I'm here0" << endl;
                // solve the approximate Riemann problem for this time step
                // reconstruct numerical fluxes at faces in a proper way
                //Godunov.update(secondOrder);
                //Info << "For cell 7943, h=" << h.internalField()[7943] << " p=" << p.internalField()[7943] 
                //     << " T=" << T.internalField()[7943] << " Cp=" << thermo.Cp()().internalField()[7943]
                //     << " Cv=" << thermo.Cv()().internalField()[7943] 
                //    << " rho=" << rho.internalField()[7943] << " U=" << U.internalField()[7943] << nl 
                //    << " rhoh-p= " << rho.internalField()[7943]*h.internalField()[7943]-p.internalField()[7943] 
                //    << " p/(kappa-1) = " << p.internalField()[7943]/(thermo.Cp()().internalField()[7943]/thermo.Cv()().internalField()[7943] -1.) 
                //    << " kappa=" << thermo.Cp()().internalField()[7943]/thermo.Cv()().internalField()[7943]  
                //    << " A=" << A.internalField()[7943]
                //    << "---------------------------"<<endl;
                Godunov.update(secondOrder);
                //h.write();

                // get access to multi-stage coefficient in ddt scheme
                physDeltaT[0] = beta[i];

                // \f$ \mu \left( \nabla \vec{U} + \left( \nabla \vec{U}
                // \right)^T - \frac{2}{nD} \left( \nabla \bullet \vec{U}
                // \right) I \right) \f$
                // nD = 3 in three dimensional space
                // for two- and one-dimensional computations this
                // implementation is wrong
                // is equal to -turbulence->devRhoReff(), but we do not need to
                // evaluate gradU twice
                //Info << "I'm here" << endl;
                const volSymmTensorField tau
                (
                    "tau",
                    -turbulence->devRhoReff()
                    -((2.0/3.0)*I)*rho*turbulence->k()
                );

                volScalarField rhoFlux = -fvc::div(Godunov.rhoFlux());

                volVectorField rhoUFlux = -fvc::div(Godunov.rhoUFlux())
                    + fvc::div(tau);

                volScalarField rhoEFlux = -fvc::div(Godunov.rhoEFlux())
                    // Viscous heating with
                    + fvc::div( tau & U )
                    // Fourier law with static enthalpy
                    // with alphaEff() - the effective turbulent
                    // thermal diffusivity.
                    + fvc::laplacian(turbulence->alphaEff(), h)
                    // Molecular Diffusion and Turbulent Transport closure
                    // Wilcox (2006): Chapter 5.4.3
                    // should be better used DkEff(F1) instead of muEff(), but
                    // this function is not virtual, now it is assumed that
                    // \sigma_k = 5/3 is hard coded
                    + fvc::laplacian
                      (
                          (turbulence->mu()+0.6*turbulence->mut()),
                          turbulence->k()
                      )
                    ;

                // Add source terms to the numerical fluxes due to the rotation
                // Is not needed if mesh is really rotating, as this effect is
                // captured then by the time stepping of the velocity field. As
                // the vector is not adapted, when the mesh is rotating.
                mesh.type() == "staticFvMesh" ?
                    mrfZones.addCoriolis(rho,U,rhoUFlux):void();

                // Pseudo time integration
                solve
                (
                    fvm::ddt(rho) == rhoFlux
                );

                // time integration
                solve
                (
                    fvm::ddt(rhoU) == rhoUFlux
                );

                // time integration
                solve
                (
                    fvm::ddt(rhoE) == rhoEFlux
                );
                
                //Info << "TMin:" << min(T) << "  TMax:" << max(T) << endl;
                // bound density
                boundMinMax(rho,rhoMin,rhoMax);

                // bound rhoE
                boundMinMax(rhoE,rhoEMin,rhoEMax);
                //Info << "I'm Here-2" << endl;
                // Compute internal field of U
                U.dimensionedInternalField() =
                    rhoU.dimensionedInternalField()
                   /rho.dimensionedInternalField();


                //Info << "I'm Here-3" << endl;
                // Update boundary conditions of U
                U.correctBoundaryConditions();
                //Info << "I'm Here-4" << endl;
                // Bound the velocity
                volScalarField magU = mag(U);


                //Info << "I'm Here -1" << endl;
                if (max(magU) > UMax)
                {
                    //Info<< "bounding " << U.name()
                    //    << " max: " << max(magU).value()
                    //    << endl;

                    volScalarField Ulimiter = pos(magU - UMax)*UMax
                        /(magU + smallU) + neg(magU - UMax);
                    Ulimiter.max(scalar(0));
                    Ulimiter.min(scalar(1));

                    U *= Ulimiter;
                    U.correctBoundaryConditions();
                }


                //Info << "I'm Here" << endl;
                // store tmp fields in order to prevent computing twice
                //const volScalarField kappa(thermo.Cp()/thermo.Cv());

//                 // bound static energy
// //                 volScalarField e = rhoE/rho - 0.5*magSqr(U);
//                 volScalarField e = rhoE/rho - 0.5*magSqr(U) - TKE;
// //                 e = rhoE/rho - 0.5*magSqr(U) - turbulence->k();
//                 e.correctBoundaryConditions();
//                 boundMinMax(e,eMin,eMax);

                // Update static enthalpy:
                // The turbulent kinetic energy k is part of the total energy E
                // Therefore it needs to be subtracted from E in order to get
                // the static enthalpy h
//                 h = kappa*e;

                // add
                // h = function_h(rho, e)  // lookup able with rho and e on x/y axis and h on z-axis

                //h = kappa*(rhoE/rho - 0.5*magSqr(U) - turbulence->k());

                // Updating the e field;//
                e = rhoE/rho - 0.5*magSqr(U) - turbulence->k(); // <-

                //    Info << "                       For cell 42468, h=" << h.internalField()[42468] << " p=" << p.internalField()[42468] 
                //     << " T=" << T.internalField()[42468] 
                ///     << " rho=" << rho.internalField()[42468]   << " e="  << e.internalField()[42468] << endl;



                forAll(h,celli)                                // <-
                {                                              // <-
                    h[celli] = hrhoeTable(rho[celli],e[celli]);// <-
                }
                //h = (thermo.Cp()/thermo.Cv())*(rhoE/rho - 0.5*magSqr(U) - turbulence->k());
                // correct boundary conditions of static enthalpy

                /*******calculate boundary value*******/  

                forAll(h.boundaryField(), patchI)
                  {
                    //if (isA<fixedValueFvPatchScalarField>(h.boundaryField()[patchI]))
                    //{}
                    //else
                   // {
                        fvPatchScalarField& hb = h.boundaryField()[patchI];
                        const fvPatchScalarField& rhob = rho.boundaryField()[patchI];
                        const fvPatchScalarField& eb = e.boundaryField()[patchI];
                        forAll(hb, faceI)
                        {
                            hb[faceI] = hrhoeTable(rhob[faceI], eb[faceI]);
                        }
                    //}
                  }

                h.correctBoundaryConditions();

                // bound enthalpy
                //Info << "Before boundMinMax(h,hMin,hMax)" << endl;
                boundMinMax(h,hMin,hMax);

                // compute complete field of p
//                 p = (kappa - 1.0)*rho*e;
                //p = (1.0 - 1.0/kappa)*rho*h;


                // add
                //Info << " - -  -update p and h- - - - - - - -- - " << nl << "compare rho with thermo.rho()" << nl << 
                //"Boundary: rho:" << rho.boundaryField()[5][50] << " thermo.rho():" << thermo.rho()().boundaryField()[5][50] 
                //<< " T:"<< T.boundaryField()[5][50]<< "T.thermo" << thermo.T().boundaryField()[5][50] 
                //<< " p:"<< p.boundaryField()[5][50] << "thermo.p" << thermo.p().boundaryField()[5][50]
                //<< " h:"<< h.boundaryField()[5][50] << "thermo.h" << thermo.h().boundaryField()[5][50] << nl
                //<< "Internal: rho:" << rho.internalField()[24630] << " and thermo.rho:" << thermo.rho()().internalField()[24630] 
                //<< " T:"<< T.internalField()[24630] << " Thermo.T" << thermo.T().internalField()[24630] 
                //<< " p:"<< p.internalField()[24630] << " Thermo.p" << thermo.p().internalField()[24630]
                //<< " h:"<< h.internalField()[24630] << "thermo.h" << thermo.h().internalField()[24630]  <<  endl;


                forAll(p,celli)                                 // <-
                {                                               // <-
                    p[celli] = prhoeTable(rho[celli],e[celli]); // <-
                }                                               // <-

                /*******calculate boundary value*******/  
                //forAll(p.boundaryField(), patchI)
               // {
               //     if (!isA<fixedValueFvPatchScalarField>(p.boundaryField()[patchI]))
               //     {
               //         fvPatchScalarField& pb = p.boundaryField()[patchI];
               //         const fvPatchScalarField& rhob = rho.boundaryField()[patchI];
               //         const fvPatchScalarField& eb = e.boundaryField()[patchI];
               //         forAll(pb, faceI)
               //         {
               //             pb[faceI] = prhoeTable(rhob[faceI], eb[faceI]);
               //         }
               //     }
               // }
                

                // h = function_h(rho, e)  // lookup able with rho and e on x/y axis and h on z-axis
                //p = (thermo.Cp() - thermo.Cv()) * rho * T ;
                // correct boundary conditions of p

                p.correctBoundaryConditions();

                // bound pressure
                boundMinMax(p,pMin,pMax);
                   // Info << "After h and p update : For cell 7943, h=" << h.internalField()[24630] << " p=" << p.internalField()[24630] 
                   //  << " T=" << T.internalField()[24630] << " U=" << U.internalField()[24630]
                   // << " rho=" << rho.internalField()[24630]<< " rhoE=" << rhoE.internalField()[24630]  
                    //<< " e="<< e.internalField()[24630]<< endl;
                // correct thermo physical properties
                // therefore one needs correct p and e fields
                

                //forAll(T,celli)
                //{
                //    T[celli] = TrhoeTable(rho[celli],e[celli]);
               // }

                //T.correctBoundaryConditions();
//Info << "checking p and T value" << T.boundaryField()[0][0] << "and" << thermo.T().boundaryField()[0][0]<<endl;
                thermo.correct();
//                    Info << "After thermo correct : For cell 7943, h=" << h.internalField()[24630] << " p=" << p.internalField()[24630] 
//                     << " T=" << T.internalField()[24630] << " U=" << U.internalField()[24630]
//                     << " rho=" << rho.internalField()[24630]<< " rhoE=" << rhoE.internalField()[24630] 
//                     << " e="<< e.internalField()[24630] <<endl;
                
                  //  Info << "                       For cell 42468, h=" << h.internalField()[42468] << " p=" << p.internalField()[42468] 
                  //   << " T=" << T.internalField()[42468]  
                  //   << " rho=" << rho.internalField()[42468]   << " e="  << e.internalField()[42468] << endl;
                //Info << " - -  - After thmero correct - - - - - - -- - " << nl << "compare rho with thermo.rho()" << nl << 
                //"Boundary: rho:" << rho.boundaryField()[5][50] << " thermo.rho():" << thermo.rho()().boundaryField()[5][50] 
                //<< " T:"<< T.boundaryField()[5][50]<< "T.thermo" << thermo.T().boundaryField()[5][50] 
                //<< " p:"<< p.boundaryField()[5][50] << "thermo.p" << thermo.p().boundaryField()[5][50]
                //<< " h:"<< h.boundaryField()[5][50] << "thermo.h" << thermo.h().boundaryField()[5][50] << nl
                //<< "Internal: rho:" << rho.internalField()[24630] << " and thermo.rho:" << thermo.rho()().internalField()[24630] 
                //<< " T:"<< T.internalField()[24630] << " Thermo.T" << thermo.T().internalField()[24630] 
                //<< " p:"<< p.internalField()[24630] << " Thermo.p" << thermo.p().internalField()[24630]
                //<< " h:"<< h.internalField()[24630] << "thermo.h" << thermo.h().internalField()[24630]  <<  endl;
//Info << " For checking the rho and thermo.rho() value:" << rho.boundaryField()[1][50] << " and " << thermo.rho()().boundaryField()[1][50] 
//     <<  "And T.boundary"<< T.boundaryField()[1][50] << "p.boundary" << p.boundaryField()[1][50]<< endl;

               

                // Update boundary field of rho
                rho.boundaryField() = thermo.rho()().boundaryField();
                rho.correctBoundaryConditions();
                /*******calculate boundary value*******/  
                //forAll(rho.boundaryField(), patchI)
               // {
               //     if (!isA<fixedValueFvPatchScalarField>(rho.boundaryField()[patchI]))
               //     {
               //         fvPatchScalarField& rhob = rho.boundaryField()[patchI];
               //        const fvPatchScalarField& pb = p.boundaryField()[patchI];
               //         const fvPatchScalarField& Tb = T.boundaryField()[patchI];
               //         forAll(rhob, faceI)
               //        {
               //             rhob[faceI] = rhoTable(Tb[faceI], pb[faceI]);
               //         }
               //     }
               // }






//Info << " For checking the rho and thermo.rho() value:" << rho.boundaryField()[1][50] << " and " << thermo.rho()().boundaryField()[1][50] << endl;

                // Update boundary field of rhoU
                rhoU.boundaryField() = rho.boundaryField()*U.boundaryField();
                

                //const volScalarField kappa(thermo.Cp()/thermo.Cv()); //<-
                // Update boundary field of rhoE
                rhoE.boundaryField() =
                    rho.boundaryField()*
                    (
                        h.boundaryField() +
                       0.5*magSqr(U.boundaryField())
                      + turbulence->k()().boundaryField()
                    ) - p.boundaryField();
                 // + p.boundaryField()/(kappa.boundaryField()-1.0);


                  //rho*(h + 0.5*magSqr(U) + turbulence->k()) - p

                                                                                                                                                                             // needed for turbulence and CoEuler ddt scheme
                // and maybe elsewhere;
                // phi is already realtive to the mesh motion, if we are using
                // the ALE approach
                phi = Godunov.rhoFlux();
//                 mrfZones.relativeFlux(fvc::interpolate(rho),phi);

                // Convergence check
                // Averaged L2-Norm of fluxes
                scalar L2NormRho = max(Foam::sqrt(sum(sqr(rhoFlux.internalField()))
                    /mesh.nCells()),SMALL);
                scalar LInfNormRho = max(mag(rhoFlux.internalField()));
//                 scalar L2NormRhoU = max(Foam::sqrt(sum(magSqr(rhoUFlux.internalField()))
//                     /mesh.nCells()),SMALL);
//                 scalar L2NormRhoE = max(Foam::sqrt(sum(sqr(rhoEFlux.internalField()))
//                     /mesh.nCells()),SMALL);



                    //<<//<< max(h) << " and " << min(h) << nl
                   //  << "The max and min h is" << max(h) << "and" << min(h) << nl
                    // << "The velocity is: " << max(magU)<< "and " << min(magU) << nl
                     //<< "The kappa is :   "    << max(kappa) << " and " << min(kappa) <<  nl 
                    // << "The Cp is :   "    << thermo.Cp()().internalField()[8099] << nl// max(thermo.Cp()) << " and " << min(thermo.Cp()) << nl
                     //<< "The Cv is :   "  << max(thermo.Cv()) << " and " << min(thermo.Cv()) << nl
                    // << "The p  is :   " <<   p[8099] << nl//max(p) << " and " << min(p) << nl
                    // << "The rho  is :   " <<  rho[8099] << nl//  max(rho) << " and "   << min(rho)   << nl
                     //<< "The T is :"
                     //<< endl;  //Jianhui Qi  





                // At the final stage, update the acoustic speed for the utility of next time step.
                forAll(A,celli)                                 // <-
                {                                               // <-
                    A[celli] = ArhoeTable(rho[celli],e[celli]);       // <-
                }                                               // <-

                /*******calculate boundary value*******/  
                forAll(A.boundaryField(), patchI)
                {
                    fvPatchScalarField& Ab = A.boundaryField()[patchI];
                    const fvPatchScalarField& rhob = rho.boundaryField()[patchI];
                    const fvPatchScalarField& eb = e.boundaryField()[patchI];
                      forAll(Ab, faceI)
                      {
                        Ab[faceI] = ArhoeTable(rhob[faceI], eb[faceI]);
                     }
                }

                    //Info << "    At final stage : For cell 7943, h=" << h.internalField()[2144] << " p=" << p.internalField()[2144] 
                    // << " T=" << T.internalField()[2144] << " U=" << U.internalField()[2144]
                    // << " rho=" << rho.internalField()[2144]<< " rhoE=" << rhoE.internalField()[2144]  << " e="<< e.internalField()[2144]
                    // << " A=" << A.internalField()[2144] << endl;         

                // Averaged L2-Norm of fluxes
                Info<< "rho L2 Residual: "
                    << Foam::log10(L2NormRho)  << endl
                    << "rho LInf Residual: "
                    << Foam::log10(LInfNormRho) << endl
//                     << "rhoU residual: "
//                     << Foam::log10(L2NormRhoU) << endl
//                     << "rhoE residual: "
//                     << Foam::log10(L2NormRhoE) << endl
                    << endl;

            // End outer Loop for Runge - Kutta
            }

            // Update relative velocity
            mrfZones.relativeVelocity(U,URel);

            //Update turbulence after the multi-stage time integration
            //Info <<"before turbulence" << endl;
            turbulence->correct();
            //Info << "afeer turbuelnce" << endl;
        }
        // End dual-time stepping
        runTime.endSubCycle();
#       include "updateDualTimeSteppingFields.H"
}
