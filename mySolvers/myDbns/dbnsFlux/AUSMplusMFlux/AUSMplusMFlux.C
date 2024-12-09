/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "AUSMplusMFlux.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
    defineTypeNameAndDebug(AUSMplusMFlux, 0);
    addToRunTimeSelectionTable(dbnsFlux, AUSMplusMFlux, dictionary);
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::AUSMplusMFlux::AUSMplusMFlux(const fvMesh&, const dictionary& dict)
{
    dictionary mySubDict( dict.subOrEmptyDict("AUSMplusMFluxCoeffs") );
    alpha0_  = mySubDict.lookupOrAddDefault("alpha0", 3.0/16.0);
    beta_    = mySubDict.lookupOrAddDefault("beta", 1.0/8.0);
    MaInf_   = mySubDict.lookupOrAddDefault("MaInf", 0.1);
    
    if (mySubDict.lookupOrDefault("printCoeffs", false))
        Info << mySubDict << nl;
};


void Foam::AUSMplusMFlux::evaluateFlux
(
    scalar& rhoFlux,
    vector& rhoUFlux,
    scalar& rhoEFlux,
    const scalar& pLeft,
    const scalar& pRight,
    const vector& ULeft,
    const vector& URight,
    const scalar& TLeft,
    const scalar& TRight,
    const scalar& RLeft,
    const scalar& RRight,
    const scalar& CvLeft,
    const scalar& CvRight,
    const vector& Sf,
    const scalar& magSf,
    const scalar& h
) const
{
    // Step 1: decode left and right:
    // normal vector
    const vector normalVector = Sf/magSf;
    
    const scalar qMesh = meshPhi / magSf;

    // Ratio of specific heat capacities
    const scalar kappaLeft = (RLeft + CvLeft)/CvLeft;
    const scalar kappaRight = (RRight + CvRight)/CvRight;
    const scalar kappa = 0.5 * (kappaLeft + kappaRight);

    // Density
    const scalar rhoLeft = pLeft/(RLeft*TLeft);
    const scalar rhoRight = pRight/(RRight*TRight);

    // DensityTotalEnthalpy
    const scalar rhoHLeft  = rhoLeft*(CvLeft*TLeft + 0.5*magSqr(ULeft)) + pLeft;
    const scalar rhoHRight = rhoRight*(CvRight*TRight + 0.5*magSqr(URight)) + pRight;

    // DensityVelocity
    const vector rhoULeft = rhoLeft*ULeft;
    const vector rhoURight = rhoRight*URight;

    // Compute qLeft and qRight (q_{l,r} = U_{l,r} \bullet n)
    const scalar qLeft     = (ULeft & normalVector) - qMesh;
    const scalar qRight    = (URight & normalVector) - qMesh;
    const scalar qLeftSqr  = sqr(qLeft);
    const scalar qRightSqr = sqr(qRight);
    const scalar qLeftMag  = mag(qLeft);
    const scalar qRightMag = mag(qRight);


    // Eq. 29
    const scalar sqrULeft  = magSqr(ULeft);
    const scalar sqrURight = magSqr(URight);
    const scalar hnormal   = 0.5*(rhoHLeft/rhoLeft-0.5*sqrULeft+rhoHRight/rhoRight-0.5*sqrURight);
    const scalar aStar     = sqrt(2.0*(kappa-1)/(kappa+1) * hnormal);
    const scalar aStarSqr  = sqr(aStar);
    // If 1/2(ULeft + URight) >= 0
    const scalar qLeftRight = 0.5*(qLeft+qRight);
    const scalar aTilde    = (qLeftRight  >= 0.0) ? aStarSqr/(max(qLeftMag, aStar)) : aStarSqr/(max(qRightMag, aStar));
    const scalar aTildeSqr = sqr(aTilde);

    const scalar rhoTilde = 0.5*(rhoLeft+rhoRight);
    
    const scalar sqrMaDash = (sqr(qLeft)+sqr(qRight))/(2.0*sqr(aTilde));
    const scalar sqrMaZero = min(1.0,max(sqrMaDash,sqr(MaInf_)));
    const scalar MaZero    = Foam::sqrt(sqrMaZero);
    const scalar sqrMaInf  = sqr(MaInf_)

    const scalar fa = MaZero*(2.0-MaZero);

    
    const scalar MaRelLeft  = qLeft /aTilde;
    const scalar MaRelRight = qRight/aTilde;
    
    const scalar magMaRelLeft  = mag(MaRelLeft);
    const scalar magMaRelRight = mag(MaRelRight);
    
    const scalar Ma1PlusLeft   = 0.5*(MaRelLeft +magMaRelLeft );
    const scalar Ma1MinusRight = 0.5*(MaRelRight-magMaRelRight);
    
    const scalar Ma2PlusLeft   =  0.25*sqr(MaRelLeft +1.0);
    const scalar Ma2PlusRight  =  0.25*sqr(MaRelRight+1.0);
    const scalar Ma2MinusLeft  = -0.25*sqr(MaRelLeft -1.0);
    const scalar Ma2MinusRight = -0.25*sqr(MaRelRight-1.0);
    
    const scalar Ma4BetaPlusLeft   = ((magMaRelLeft  >= 1.0) ? Ma1PlusLeft   : (Ma2PlusLeft  *(1.0-16.0*beta_*Ma2MinusLeft)));
    const scalar Ma4BetaMinusRight = ((magMaRelRight >= 1.0) ? Ma1MinusRight : (Ma2MinusRight*(1.0+16.0*beta_*Ma2PlusRight)));

    // $\psi$ from Eq. 6 (Chen, 2020)
    const scalar P5alphaPlusLeft   = ((magMaRelLeft  >= 1.0) ?
    (Ma1PlusLeft/MaRelLeft)    : (Ma2PlusLeft  *(( 2.0-MaRelLeft) -16.0*alpha0_*MaRelLeft *Ma2MinusLeft )));
    const scalar P5alphaMinusRight = ((magMaRelRight >= 1.0) ?
    (Ma1MinusRight/MaRelRight) : (Ma2MinusRight*((-2.0-MaRelRight)+16.0*alpha0_*MaRelRight*Ma2PlusRight)));

    // qLeft \equiv ULeft
    // $M$ from Eq. 15
    const scalar M = min(1.0, max(magMaRelLeft, magMaRelRight));

    // $f$ from Eq. 15
    const scalar f = (1-Foam::cos(M_PI*M))/2;

    // $\delta p$ from Eq. 14
    const scalar deltaP = pRight-pLeft;

    // $h$ from Eq. 25
    const scalar h = eval_h_min(pLeft, pRight);
    const scalar g = 0.5*(1+Foam::cos(M_PI*h))

    // $M_p$ from Eq. 14
    const scalar Mp = -0.5*(1-f)*deltaP/(rhoTilde*aTildeSqr)*(1-g)

    // $f_0$ from Eq. 20
    const scalar f0 = min(1.0, max(f, sqrMaInf))

    // $p_s$ from Eq. 19
    const scalar ps = 
      0.5*(pLeft+pRight)
    + (P5alphaPlusLeft - P5alphaMinusRight)*0.5*(pLeft-Pright)
    + f0*(P5alphaPlusLeft+P5alphaMinusRight-1)*0.5*(pLeft+pRight);

    const vector pU = -g*(kappa*(pLeft+pRight))/(2*aTilde)*P5alphaPlusLeft*P5alphaMinusRight*(URight-ULeft);

    const scalar MaRelTilde = Ma4BetaPlusLeft + Ma4BetaMinusRight + Mp;
    const scalar URelTilde = MaRelTilde*aTilde;
    const scalar magURelTilde = mag(MaRelTilde)*aTilde;

    // There is a typo in Luo et. al, J. Comp. Physics 194 (2004), Chap 4.2 Eq. 4.8
    // refer to the origial Paper from Liou, J. Comp. Physics 129 (1996), Chap4, Eq. 42
    rhoFlux  = (0.5*(URelTilde*(rhoLeft +rhoRight) -magURelTilde*(rhoRight -rhoLeft)))*magSf;
    rhoUFlux = (0.5*(URelTilde*(rhoULeft+rhoURight)-magURelTilde*(rhoURight-rhoULeft))+pTilde*normalVector)*magSf;
    rhoEFlux = (0.5*(URelTilde*(rhoHLeft+rhoHRight)-magURelTilde*(rhoHRight-rhoHLeft)) + pTilde*qMesh)*magSf;

}

Foam::surfaceScalarField AUSMplusMFlux::eval_h_min
(
    surfaceScalarField p_L,
    surfaceScalarField p_R,
) const
{
    // $h$, Pressure sensing variable, Eq. 25
    const surfaceScalarField h_k = min(p_L/p_R, p_R/p_L);
    forAll (h_min, cellI)
    {
        // Get list of faces of current cell
        const labelList &cellFaces = mesh_.cells()[cellI]; 
        // Iterate over all faces of celli
        forAll (cellFaces, faceI)
        {
            // find minimum h_k of the cell; ONLY FOR INTERNAL BOUNDARIES
            if (cellFaces[faceI] < h_k.size())
            {
                h_min[cellI] = min(h_min[cellI], h_k[cellFaces[faceI]]);
            }
        }
    };

     // Eq. 43: interpolate on cell faces
    const surfaceScalarField h_final
    (
        "h_final", 
        min
        (
            fvc::surfaceReconstruct(h_min, own_, RECONSTRUCT_SCALAR_),
            fvc::surfaceReconstruct(h_min, nei_, RECONSTRUCT_SCALAR_)
        )
    );

    return h_final;
}


// ************************************************************************* //
