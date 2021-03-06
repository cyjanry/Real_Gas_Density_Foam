/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Application
    hllcALERealFlux

Author
    Oliver Borm  All rights reserved.

\*---------------------------------------------------------------------------*/

#include "hllcALERealFlux.H"
//#include "extrapolation2DTable.H"
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //+


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::hllcALERealFlux::evaluateFlux
(
    scalar& rhoFlux,
    vector& rhoUFlux,
    scalar& rhoEFlux,
    const scalar pLeft,
    const scalar pRight,
    const vector ULeft,
    const vector URight,

    const scalar rhoLeft,
    const scalar rhoRight,
// //
    const scalar kLeft,
    const scalar kRight,

    const scalar kappaLeft,
    const scalar kappaRight,

    const scalar ALeft,
    const scalar ARight,
//
    const scalar hLeft,
    const scalar hRight,
// //
    const vector Sf,
    const scalar magSf,
    const vector dotX,
    const scalar Konstant
) const
{
        // bounding variables
        const scalar rhoMin = SMALL;

        // Step 1: decode left and right:
        // normal vector
        const vector normalVector = Sf/magSf;

        // speed of sound, for left and right side, assuming perfect gas
        const scalar aLeft  = max(SMALL, ALeft);//Foam::sqrt(max(SMALL,kappaLeft *pLeft /max(rhoLeft, rhoMin)));
        const scalar aRight = max(SMALL, ARight);//Foam::sqrt(ma(SMALL,kappaRight*pRight/max(rhoRight,rhoMin)));
        //const scalar aLeft  = Foam::sqrt(max(SMALL,kappaLeft *pLeft /max(rhoLeft, rhoMin))); //<-old
        //const scalar aRight = Foam::sqrt(max(SMALL,kappaRight*pRight/max(rhoRight,rhoMin)));  //<-old

//         const scalar aLeft  = Foam::sqrt(kappaLeft *pLeft /rhoLeft);
//         const scalar aRight = Foam::sqrt(kappaRight*pRight/rhoRight);
       //Info << "It's little be slow here, aright ="<< aRight << endl;
        // DensityTotalEnergy
        const scalar rhoELeft =
            //pLeft /(kappaLeft -1.0)+rhoLeft *(0.5*magSqr(ULeft) +kLeft);
            rhoLeft * hLeft -pLeft + rhoLeft* (0.5*magSqr(ULeft) +kLeft) ;
        const scalar rhoERight =
            //pRight/(kappaRight-1.0)+rhoRight*(0.5*magSqr(URight)+kRight);
            rhoRight * (hRight + 0.5*magSqr(URight)+kRight ) -pRight;


        //Info << "rhoELeft_0 ="  << rhoELeft << " and hLeft=" <<hLeft << " hRight=" << hRight << nl
        //     << "1term=" <<  rhoLeft * hLeft -pLeft << " 2term=" << rhoLeft* (0.5*magSqr(ULeft) +kLeft) << endl;

        //scalar rhoELeft =
        //    pLeft /(kappaLeft -1.0) + rhoLeft *(0.5*magSqr(ULeft) +kLeft);

        //scalar rhoERight =
        //    pRight/(kappaRight-1.0) + rhoRight*(0.5*magSqr(URight)+kRight);
         //Info << "rhoELeft_1 ="  << rhoELeft << " pLeft="<< pLeft << "pRight" << pRight <<nl
         //<< "1term=" <<  pLeft /(kappaLeft -1.0) << " 2term=" << rhoLeft* (0.5*magSqr(ULeft) +kLeft) << endl;



        //Info << "what is hRight:" << hRight << endl;

        // DensityVelocity
        const vector rhoULeft  = rhoLeft *ULeft;
        const vector rhoURight = rhoRight*URight;

        // Compute left and right total enthalpies:
        //const scalar HLeft  = (rhoELeft  + pLeft) /max(rhoLeft ,rhoMin);  //<-old
        //const scalar HRight = (rhoERight + pRight)/max(rhoRight,rhoMin);  //<-old

        // compute velocity relative to mesh
        const vector URelLeft  = ULeft  - dotX;
        const vector URelRight = URight - dotX;

        // compute qLeft and qRight (q_{l,r} = U_{l,r} \bullet n)
        const scalar qLeft  = (URelLeft  & normalVector);
        const scalar qRight = (URelRight & normalVector);

        // Step 2:
        // needs rho_{l,r}, U_{l,r}, H_{l,r}, kappa_{l,r}, Gamma_{l,r}, q_{l,r}
        // compute Roe weights
        const scalar rhoLeftSqrt  = Foam::sqrt(max(rhoLeft ,rhoMin));
        const scalar rhoRightSqrt = Foam::sqrt(max(rhoRight,rhoMin));

        const scalar wLeft = rhoLeftSqrt
            /max((rhoLeftSqrt + rhoRightSqrt),rhoMin);
        const scalar wRight = 1.0 - wLeft;

        // Roe averaged velocity
        const vector UTilde = wLeft*ULeft + wRight*URight;

        // Roe averaged relative face velocity
        const scalar contrURelTilde = ((UTilde - dotX) & normalVector);

        // Roe averaged contravariant velocity
        //const scalar contrUTilde = (UTilde & normalVector);// <-old

        // Roe averaged total enthalpy
        //const scalar HTilde = wLeft*HLeft + wRight*HRight; // <-old

        // Roe averaged kappa
        // TODO: needs to be verified for non constant kappa values
        //const scalar kappaTilde = wLeft*kappaLeft + wRight*kappaRight;
        //const scalar kappaTilde = wLeft*kappaLeft + wRight*kappaRight;//<-old

        // Speed of sound with Roe reconstruction values
       // const scalar aTilde =   // <-old
       //     Foam::sqrt(max(SMALL,(kappaTilde-1.0)*(HTilde-0.5*sqr(contrUTilde))));//<-old


      const scalar eta2 = 0.5*rhoLeftSqrt*rhoRightSqrt/( sqr(rhoLeftSqrt+rhoRightSqrt) );

        
      const scalar aTilde = Foam::sqrt(wLeft*sqr(aLeft) + wRight*sqr(aRight) + eta2*sqr(qLeft-qRight)) ;
        //scalar Residual = (aTilde - aLeft)*(aTilde - aRight);

        //if (Residual > 1e-7)
        //{

        //        Info << "ERROR in calculating aTilda 2 " 
        //        << " aTilde="<< aTilde << " aRight=" << aRight << " aLeft=" << aLeft
        //        << " qLeft=" << qLeft << "  qRight=" << qRight <<  " ULeft" << ULeft << "URelLeft" << URelLeft << endl;

        //}


        // Step 3: compute signal speeds for face:
        const scalar SLeft  = min(qLeft -aLeft,  contrURelTilde-aTilde);
 
       // Info << "Signal Speed SLeft=" << SLeft << " qLeft -aLeft=" << qLeft -aLeft << " contrURelTilde-aTilde=" << contrURelTilde-aTilde << endl;//

        const scalar SRight = max(qRight+aRight, contrURelTilde+aTilde);
       // Info << "Signal Speed SRight=" << SRight << " qRight+aRight=" << qLeft -aLeft << " contrURelTilde+aTilde=" << contrURelTilde-aTilde << endl;
        // See for an alternative computation of wave speed: Toro (1999) Eq. 10.46-10.48

        const scalar SStar = (rhoRight*qRight*(SRight-qRight)
            - rhoLeft*qLeft*(SLeft - qLeft) + pLeft - pRight )/
            stabilise((rhoRight*(SRight-qRight)-rhoLeft*(SLeft-qLeft)),SMALL);

        // compute pressure in star region from the right side
        const scalar pStarRight =
            rhoRight * (qRight - SRight) * (qRight - SStar) + pRight;

        // should be equal to the left side
        const scalar pStarLeft  =
            rhoLeft  * (qLeft -  SLeft)  * (qLeft - SStar)  + pLeft;

        // give a warning if this is not the case
        if ( mag(pStarRight-pStarLeft) > 1e-6 )
        {
                Info << "mag(pStarRight-pStarLeft) > SMALL " << endl;
        }

        // use pStarRight for pStar, as in theory, pStarRight == pStarLeft
        const scalar pStar = pStarRight;

        // Step 4: upwinding - compute states:
        scalar convectionSpeed = 0.0;
        scalar rhoState = 0.0;
        vector rhoUState = vector::zero;
        scalar rhoEState = 0.0;
        scalar pState = 0.0;

        // TODO: Maybe one can use pos/neg implementation, but then one has to
        // evaluate all 4 states at each iteration!
        // label A = pos(SLeft);
        // label B = pos(SStar);
        // label C = pos(SRight);
        // please double check the bool operators again, if one want's to
        // implement this!!!
        // scalar convectionSpeed = A*B*C*qLeft+(1-A)*B*C*SStar
        //     +(1-A)*(1-B)*C*SStar+(1-A)*(1-B)*(1-C)*qRight:

        if ( pos(SLeft) )
        {
            // compute F_l
            convectionSpeed = qLeft;
            rhoState  = rhoLeft;
            rhoUState = rhoULeft;
            rhoEState = rhoELeft;
            pState = pLeft;
        }
        else if ( pos(SStar) )
        {
            scalar omegaLeft = scalar(1.0)/stabilise((SLeft - SStar),SMALL);

            // compute left star region
            convectionSpeed = SStar;
            rhoState  = omegaLeft *   (SLeft - qLeft) * rhoLeft;
            rhoUState = omegaLeft * ( (SLeft - qLeft) * rhoULeft
                + (pStar - pLeft)*normalVector );
            rhoEState = omegaLeft * ( (SLeft - qLeft) * rhoELeft
                - pLeft*qLeft + pStar*SStar );
            pState = pStar;
        }
        else if ( pos(SRight) )
        {
            scalar omegaRight = scalar(1.0)/stabilise((SRight - SStar),SMALL);

            // compute right star region
            convectionSpeed = SStar;
            rhoState  = omegaRight *   (SRight - qRight) * rhoRight;
            rhoUState = omegaRight * ( (SRight - qRight) * rhoURight
                + (pStar - pRight)*normalVector );
            rhoEState = omegaRight * ( (SRight - qRight) * rhoERight
                - pRight*qRight + pStar*SStar );
            pState = pStar;
        }
        else if ( neg(SRight) )
        {
            // compute F_r
            convectionSpeed = qRight;
            rhoState  = rhoRight;
            rhoUState = rhoURight;
            rhoEState = rhoERight;
            pState = pRight;
        }
        else
        {
            Info << "Error in HLLC Riemann solver" << endl;
        }

        rhoFlux  = (convectionSpeed*rhoState)*magSf;
        rhoUFlux = (convectionSpeed*rhoUState+pState*normalVector)*magSf;
        rhoEFlux = (convectionSpeed*(rhoEState+pState)
            +pState*(dotX & normalVector))*magSf;
}

// ************************************************************************* //
