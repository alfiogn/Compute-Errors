/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022 OpenFOAM Foundation
     \\/     M anipulation  |
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

\*---------------------------------------------------------------------------*/

#include "Error.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class Type>
template<class func>
inline typename func::result_type Foam::Error<Type>::faceProj
(
    const label& fi,
    const func& f
) const
{
    const vectorField& Cf = mesh_.faceCentres();
    const vectorField& P = mesh_.points();
    const cellList& cellFaces = mesh_.cells();
    const faceList& facePoints = mesh_.faces();

    List<vector> verts(3, Zero);

    scalar volume = 0.0;
    typename func::result_type integral = Zero;

    verts[0] = Cf[fi];

    const labelList& pointsJ = facePoints[fi];

    for (int k = 0; k < pointsJ.size(); k++)
    {
        int l = k == 0 ? pointsJ.size() - 1 : k - 1;
        verts[1] = P[pointsJ[k]];
        verts[2] = P[pointsJ[l]];

        volume += triArea(verts);
        integral += triQuadrature(verts, f);
    }

    return integral/volume;
}


template<class Type>
template<class func>
inline typename func::result_type Foam::Error<Type>::cellProj
(
    const label& ci,
    const func& f
) const
{
    const vectorField& C = mesh_.cellCentres();
    const vectorField& Cf = mesh_.faceCentres();
    const vectorField& P = mesh_.points();
    const cellList& cellFaces = mesh_.cells();
    const faceList& facePoints = mesh_.faces();

    List<vector> verts(4, Zero);

    scalar volume = 0.0;
    typename func::result_type integral = Zero;

    verts[0] = C[ci];

    const labelList& facesI = cellFaces[ci];

    for (int j = 0; j < facesI.size(); j++)
    {
        verts[1] = Cf[facesI[j]];

        const labelList& pointsJ = facePoints[facesI[j]];

        for (int k = 0; k < pointsJ.size(); k++)
        {
            int l = k == 0 ? pointsJ.size() - 1 : k - 1;
            verts[2] = P[pointsJ[k]];
            verts[3] = P[pointsJ[l]];

            volume += tetVolume(verts);
            integral += tetQuadrature(verts, f);
        }
    }

    return integral/volume;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
inline void Foam::Error<Type>::triCheck
(
    const List<vector>& verts
) const
{
    if (verts.size() != 3)
    {
        FatalErrorInFunction
            << "Polygon with vertices:"
            << verts << ", is not a triangle"
            << abort(FatalError);
    }
}


template<class Type>
inline void Foam::Error<Type>::tetCheck
(
    const List<vector>& verts
) const
{
    if (verts.size() != 4)
    {
        FatalErrorInFunction
            << "Polyhedron with vertices:"
            << verts << ", is not a tetrahedron"
            << abort(FatalError);
    }
}


template<class Type>
inline Foam::vector Foam::Error<Type>::barycentre
(
    const List<vector>& verts
) const
{
    tetCheck(verts);

    return  Foam::sum(verts)/verts.size();
}


template<class Type>
inline Foam::scalar Foam::Error<Type>::triArea
(
    const List<vector>& verts
) const
{
    triCheck(verts);

    vector a = verts[1] - verts[0];
    vector b = verts[2] - verts[0];

    return  Foam::mag(a ^ b)/2;
}


template<class Type>
inline Foam::scalar Foam::Error<Type>::tetVolume
(
    const List<vector>& verts
) const
{
    tetCheck(verts);

    vector a = verts[1] - verts[0];
    vector b = verts[2] - verts[0];
    vector c = verts[3] - verts[0];

    return  Foam::mag(a & (b ^ c))/6.0;
}


template<class Type>
inline Foam::List<Type> Foam::Error<Type>::triLinearCoeffs
(
    const List<vector>& verts,
    const List<Type>& values
) const
{
    triCheck(verts);

    const scalar& x0 = verts[0][0];
    const scalar& x1 = verts[1][0];
    const scalar& x2 = verts[2][0];
    const scalar& y0 = verts[0][1];
    const scalar& y1 = verts[1][1];
    const scalar& y2 = verts[2][1];
    const Type& v0 = values[0];
    const Type& v1 = values[1];
    const Type& v2 = values[2];

    Type a =
    (
        v0*y1 - v0*y2 - v1*y0 + v1*y2 + v2*y0 - v2*y1
    )/(
        x0*y1 - x0*y2 - x1*y0 + x1*y2 + x2*y0 - x2*y1 + small
    );
    Type b =
    (
        -v0*x1 + v0*x2 + v1*x0 - v1*x2 - v2*x0 + v2*x1
    )/(
        x0*y1 - x0*y2 - x1*y0 + x1*y2 + x2*y0 - x2*y1 + small
    );
    Type c =
    (
        v0*x1*y2 - v0*x2*y1 - v1*x0*y2 + v1*x2*y0 + v2*x0*y1 - v2*x1*y0
    )/(
        x0*y1 - x0*y2 - x1*y0 + x1*y2 + x2*y0 - x2*y1 + small
    );

    return List<Type>({a, b, c});
}


template<class Type>
inline Foam::List<Type> Foam::Error<Type>::tetLinearCoeffs
(
    const List<vector>& verts,
    const List<Type>& values
) const
{
    tetCheck(verts);

    const scalar& x0 = verts[0][0];
    const scalar& x1 = verts[1][0];
    const scalar& x2 = verts[2][0];
    const scalar& x3 = verts[3][0];
    const scalar& y0 = verts[0][1];
    const scalar& y1 = verts[1][1];
    const scalar& y2 = verts[2][1];
    const scalar& y3 = verts[3][1];
    const scalar& z0 = verts[0][2];
    const scalar& z1 = verts[1][2];
    const scalar& z2 = verts[2][2];
    const scalar& z3 = verts[3][2];
    const Type& v0 = values[0];
    const Type& v1 = values[1];
    const Type& v2 = values[2];
    const Type& v3 = values[3];

    Type a =
    (
        v0*y1*z2 - v0*y2*z1 - v1*y0*z2 + v1*y2*z0 + v2*y0*z1 - v2*y1*z0
      - v0*y1*z3 + v0*y3*z1 + v1*y0*z3 - v1*y3*z0 - v3*y0*z1 + v3*y1*z0
      + v0*y2*z3 - v0*y3*z2 - v2*y0*z3 + v2*y3*z0 + v3*y0*z2 - v3*y2*z0
      - v1*y2*z3 + v1*y3*z2 + v2*y1*z3 - v2*y3*z1 - v3*y1*z2 + v3*y2*z1
    )/(
        x0*y1*z2 - x0*y2*z1 - x1*y0*z2 + x1*y2*z0 + x2*y0*z1 - x2*y1*z0
      - x0*y1*z3 + x0*y3*z1 + x1*y0*z3 - x1*y3*z0 - x3*y0*z1 + x3*y1*z0
      + x0*y2*z3 - x0*y3*z2 - x2*y0*z3 + x2*y3*z0 + x3*y0*z2 - x3*y2*z0
      - x1*y2*z3 + x1*y3*z2 + x2*y1*z3 - x2*y3*z1 - x3*y1*z2 + x3*y2*z1 + small
    );
    Type b =
    -(
        v0*x1*z2 - v0*x2*z1 - v1*x0*z2 + v1*x2*z0 + v2*x0*z1 - v2*x1*z0
      - v0*x1*z3 + v0*x3*z1 + v1*x0*z3 - v1*x3*z0 - v3*x0*z1 + v3*x1*z0
      + v0*x2*z3 - v0*x3*z2 - v2*x0*z3 + v2*x3*z0 + v3*x0*z2 - v3*x2*z0
      - v1*x2*z3 + v1*x3*z2 + v2*x1*z3 - v2*x3*z1 - v3*x1*z2 + v3*x2*z1
    )/(
        x0*y1*z2 - x0*y2*z1 - x1*y0*z2 + x1*y2*z0 + x2*y0*z1 - x2*y1*z0
      - x0*y1*z3 + x0*y3*z1 + x1*y0*z3 - x1*y3*z0 - x3*y0*z1 + x3*y1*z0
      + x0*y2*z3 - x0*y3*z2 - x2*y0*z3 + x2*y3*z0 + x3*y0*z2 - x3*y2*z0
      - x1*y2*z3 + x1*y3*z2 + x2*y1*z3 - x2*y3*z1 - x3*y1*z2 + x3*y2*z1 + small
    );
    Type c =
    (
        v0*x1*y2 - v0*x2*y1 - v1*x0*y2 + v1*x2*y0 + v2*x0*y1 - v2*x1*y0
      - v0*x1*y3 + v0*x3*y1 + v1*x0*y3 - v1*x3*y0 - v3*x0*y1 + v3*x1*y0
      + v0*x2*y3 - v0*x3*y2 - v2*x0*y3 + v2*x3*y0 + v3*x0*y2 - v3*x2*y0
      - v1*x2*y3 + v1*x3*y2 + v2*x1*y3 - v2*x3*y1 - v3*x1*y2 + v3*x2*y1
    )/(
        x0*y1*z2 - x0*y2*z1 - x1*y0*z2 + x1*y2*z0 + x2*y0*z1 - x2*y1*z0
      - x0*y1*z3 + x0*y3*z1 + x1*y0*z3 - x1*y3*z0 - x3*y0*z1 + x3*y1*z0
      + x0*y2*z3 - x0*y3*z2 - x2*y0*z3 + x2*y3*z0 + x3*y0*z2 - x3*y2*z0
      - x1*y2*z3 + x1*y3*z2 + x2*y1*z3 - x2*y3*z1 - x3*y1*z2 + x3*y2*z1 + small
    );
    Type d =
    -(
        v0*x1*y2*z3 - v0*x1*y3*z2 - v0*x2*y1*z3 + v0*x2*y3*z1 + v0*x3*y1*z2
      - v0*x3*y2*z1 - v1*x0*y2*z3 + v1*x0*y3*z2 + v1*x2*y0*z3 - v1*x2*y3*z0
      - v1*x3*y0*z2 + v1*x3*y2*z0 + v2*x0*y1*z3 - v2*x0*y3*z1 - v2*x1*y0*z3
      + v2*x1*y3*z0 + v2*x3*y0*z1 - v2*x3*y1*z0 - v3*x0*y1*z2 + v3*x0*y2*z1
      + v3*x1*y0*z2 - v3*x1*y2*z0 - v3*x2*y0*z1 + v3*x2*y1*z0
    )/(
        x0*y1*z2 - x0*y2*z1 - x1*y0*z2 + x1*y2*z0 + x2*y0*z1 - x2*y1*z0
      - x0*y1*z3 + x0*y3*z1 + x1*y0*z3 - x1*y3*z0 - x3*y0*z1 + x3*y1*z0
      + x0*y2*z3 - x0*y3*z2 - x2*y0*z3 + x2*y3*z0 + x3*y0*z2 - x3*y2*z0
      - x1*y2*z3 + x1*y3*z2 + x2*y1*z3 - x2*y3*z1 - x3*y1*z2 + x3*y2*z1 + small
    );

    return List<Type>({a, b, c, d});
}


template<class Type>
inline Type Foam::Error<Type>::triLinear
(
    const List<vector>& verts,
    const List<Type>& values,
    const vector& pt
) const
{
    List<Type> coeffs = triLinearCoeffs(verts, values);

    Type interp = coeffs[0]*pt[0] + coeffs[1]*pt[1] + coeffs[2];

    return interp;
}


template<class Type>
inline Type Foam::Error<Type>::tetLinear
(
    const List<vector>& verts,
    const List<Type>& values,
    const vector& pt
) const
{
    List<Type> coeffs = tetLinearCoeffs(verts, values);

    Type interp = coeffs[0]*pt[0] + coeffs[1]*pt[1]
        + coeffs[2]*pt[2] + coeffs[3];

    return interp;
}


template<>
inline Foam::vector
Foam::Error<Foam::scalar>::gradTriLinear
(
    const List<vector>& verts,
    const List<scalar>& values
) const
{
    List<scalar> coeffs = triLinearCoeffs(verts, values);

    return vector(coeffs[0], coeffs[1], Zero);
}


template<>
inline Foam::tensor
Foam::Error<Foam::vector>::gradTriLinear
(
    const List<vector>& verts,
    const List<vector>& values
) const
{
    List<vector> coeffs = triLinearCoeffs(verts, values);

    return tensor
    (
        coeffs[0][0], coeffs[1][0], 0.0,
        coeffs[0][1], coeffs[1][1], 0.0,
        0.0, 0.0, 0.0
    );
}


template<>
inline Foam::vector
Foam::Error<Foam::scalar>::gradTetLinear
(
    const List<vector>& verts,
    const List<scalar>& values
) const
{
    List<scalar> coeffs = tetLinearCoeffs(verts, values);

    return vector(coeffs[0], coeffs[1], coeffs[2]);
}


template<>
inline Foam::tensor
Foam::Error<Foam::vector>::gradTetLinear
(
    const List<vector>& verts,
    const List<vector>& values
) const
{
    List<vector> coeffs = tetLinearCoeffs(verts, values);

    return tensor
    (
        coeffs[0][0], coeffs[1][0], coeffs[2][0],
        coeffs[0][1], coeffs[1][1], coeffs[2][1],
        coeffs[0][2], coeffs[1][2], coeffs[2][2]
    );
}


template<class Type>
template<class Type2>
inline Type2 Foam::Error<Type>::triLinearQuadrature
(
    const List<vector>& verts,
    const List<Type2>& values
) const
{
    return triArea(verts)*(values[0] + values[1] + values[2])/3.0;
}


template<class Type>
template<class Type2>
inline Type2 Foam::Error<Type>::tetLinearQuadrature
(
    const List<vector>& verts,
    const List<Type2>& values
) const
{
    return tetVolume(verts)*(values[0] + values[1] + values[2] + values[3])/4.0;
}


template<class Type>
template<class func>
inline typename func::result_type Foam::Error<Type>::triQuadrature
(
    const List<vector>& verts,
    const func& f
) const
{
    List<vector> nodes(6);
    List<scalar> bc(3);
    bc[0] = 0.65902762237;
    bc[1] = 0.23193336855;
    bc[2] = 0.10903900907;

    nodes[0] = bc[0]*verts[0] + bc[1]*verts[1] + bc[2]*verts[2];
    nodes[1] = bc[0]*verts[0] + bc[2]*verts[1] + bc[1]*verts[2];
    nodes[2] = bc[1]*verts[0] + bc[0]*verts[1] + bc[2]*verts[2];
    nodes[3] = bc[1]*verts[0] + bc[2]*verts[1] + bc[0]*verts[2];
    nodes[4] = bc[2]*verts[0] + bc[1]*verts[1] + bc[0]*verts[2];
    nodes[5] = bc[2]*verts[0] + bc[0]*verts[1] + bc[1]*verts[2];

    scalar factor(1.0/6.0);

    typename func::result_type quadVal = Zero;

    for (int i = 0; i < 6; i++)
    {
        quadVal += factor*f(nodes[i]);
    }

    return triArea(verts)*quadVal;
}


template<class Type>
template<class func>
inline typename func::result_type Foam::Error<Type>::tetQuadrature
(
    const List<vector>& verts,
    const func& f
) const
{
    List<scalar> weights(5);
    weights[0] = -16.0;
    weights[1] = 9.0;
    weights[2] = 9.0;
    weights[3] = 9.0;
    weights[4] = 9.0;

    List<vector> nodes(5);
    nodes[0] = barycentre(verts);
    nodes[1] = 0.5*verts[0] + 1.0/6.0*(verts[1] + verts[2] + verts[3]);
    nodes[2] = 0.5*verts[1] + 1.0/6.0*(verts[0] + verts[2] + verts[3]);
    nodes[3] = 0.5*verts[2] + 1.0/6.0*(verts[0] + verts[1] + verts[3]);
    nodes[4] = 0.5*verts[3] + 1.0/6.0*(verts[0] + verts[1] + verts[2]);

    scalar factor(1.0/20.0);

    typename func::result_type quadVal = Zero;

    for (int i = 0; i < 5; i++)
    {
        quadVal += factor*weights[i]*f(nodes[i]);
    }

    return tetVolume(verts)*quadVal;
}


// ************************************************************************* //
