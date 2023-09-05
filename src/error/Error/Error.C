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

#include "dimensionedType.H"
#include "Error.H"
#include "fvc.H"
#include "SortableList.H"
// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

template<class Type>
typename Foam::Error<Type>::Function
Foam::Error<Type>::Zero_ = [] (const Foam::vector& xx)
{
    return Foam::pTraits<Type>::zero;
};


template<class Type>
typename Foam::Error<Type>::GradFunction
Foam::Error<Type>::GZero_ = [] (const Foam::vector& xx)
{
    return Foam::pTraits<GradType>::zero;
};


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::Error<Type>::makePointMesh() const
{
    if (pointMeshPtr_.valid())
    {
        FatalErrorInFunction
            << "pointMesh already exists"
            << abort(FatalError);
    }

    pointMeshPtr_ = new pointMesh(mesh_);
}


template<class Type>
void Foam::Error<Type>::makeVolPointInterp() const
{
    if (volPointInterpPtr_.valid())
    {
        FatalErrorInFunction
            << "volPointInterpolator already exists"
            << abort(FatalError);
    }

    volPointInterpPtr_ = new volPointInterpolation(mesh_);
}


template<class Type>
void Foam::Error<Type>::makeFaceField() const
{
    if (faceFPtr_.valid())
    {
        FatalErrorInFunction
            << "face interpolated field already exists"
            << abort(FatalError);
    }

    faceFPtr_ =
        new surfaceTypeField
        (
            IOobject
            (
                "face" + F_.name(),
                mesh_.time().timeName(),
                mesh_.time(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            fvc::interpolate(F_)
        );
}


template<class Type>
void Foam::Error<Type>::makePointField() const
{
    if (pointFPtr_.valid())
    {
        FatalErrorInFunction
            << "point interpolated field already exists"
            << abort(FatalError);
    }

    pointFPtr_ =
        new pointTypeField
        (
            IOobject
            (
                "point" + F_.name(),
                mesh_.time().timeName(),
                mesh_.time(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            volPointInterp().interpolate(F_)
        );

}


template<class Type>
void Foam::Error<Type>::makeLinearPointField() const
{
    if (pointFPtr_.valid())
    {
        pointFPtr_.clear();
    }

    pointFPtr_ =
        new pointTypeField
        (
            IOobject
            (
                "point" + F_.name(),
                mesh_.time().timeName(),
                mesh_.time(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            pMesh(),
            dimensioned<Type>(F_.dimensions(), Zero)
        );

    //pointTypeField& pF = *pointFPtr_;
    //const vectorField& P = mesh_.points();
    //const labelListList& pointCells = mesh_.pointCells();

    //forAll(P, pi)
    //{

    //}
}


template<class Type>
void Foam::Error<Type>::performCellCut() const
{
    // check if point is inside
    const vectorField& points = mesh_.points();

    // set all point out
    boolList ptInside(mesh_.nPoints(), false);

    ptInside = surfaceSearch().calcInside(points);
    if (internal_)
    {
        forAll (ptInside, i)
        {
            ptInside[i] = !ptInside[i];
        }
    }

    forAll (mesh_.C(), i)
    {
        const List<label>& ptsi = mesh_.cellPoints()[i];

        bool allIn(true), allOut(true);

        forAll (ptsi, j)
        {
            if (ptInside[ptsi[j]])
            {
                allOut = false;
            }
            else
            {
                allIn = false;
            }
        }

        if (!allOut && !allIn)
        {
            cutMaskPtr_()[i] = 1.0;
        }
    }

    const scalar nearTol_(1.0e-15);
    DynamicList<vector> ptsTmp;
    DynamicList<DynamicList<label>> facesTmp;
    const triSurfaceSearch& surfSearch = surfaceSearch();

    cutPointsPtr_ = new List<List<vector> >(mesh_.nCells(), List<vector>(0));
    cutCellsFacesPtr_ =
        new List<List<List<label> > >
        (
            mesh_.nCells(), List<List<label> >(0)
        );

    List<List<vector> >& cutPoints = cutPointsPtr_();
    List<List<List<label> > >& cutCellsFaces = cutCellsFacesPtr_();

    const volScalarField& cutCells = cutMask();

    forAll (cutCells, ci)
    {
        if (cutCells[ci] < small)
        {
            continue;
        }

        // pre-indexing
        const vectorField& points = mesh_.points();
        const List<label>& cellPoints = mesh_.cellPoints()[ci];
        const labelHashSet cellFaces(mesh_.cells()[ci]);
        labelHashSet yetFaces;
        labelHashSet yetPoints;
        const List<List<label>>& pointFaces = mesh_.pointFaces();
        const List<face>& faces = mesh_.faces();

        // new stuff
        List<vector> newPoints;
        List<label> cutPointsI;
        int counter(0);
        List<List<label>> newFaces;
        List<label> tmpFaces;

        // loop on points of current cell
        forAll(cellPoints, pi)
        {
            const label& pointi = cellPoints[pi];

            if (!ptInside[pointi])
            {
                const List<label>& ptfacesi = pointFaces[pointi];

                forAll(ptfacesi, fk)
                {
                    if
                    (
                        !yetFaces.found(ptfacesi[fk])
                     && cellFaces.found(ptfacesi[fk])
                    )
                    {
                        yetFaces.insert(ptfacesi[fk]);

                        const labelList& facei = faces[ptfacesi[fk]];

                        for (int fj = 0; fj < facei.size(); fj++)
                        {
                            int fjp = fj == (facei.size()-1) ? 0 : fj + 1;
                            const label& ptj0 = facei[fj];
                            const label& ptj1 = facei[fjp];

                            // both points inside domain
                            if (!ptInside[ptj0] && !ptInside[ptj1])
                            {
                                newPoints.append(points[ptj0]);
                                tmpFaces.append(counter);
                                counter++;
                            }
                            // points atop surface, find intersection
                            else if (ptInside[ptj0] != ptInside[ptj1])
                            {
                                if (!ptInside[ptj0])
                                {
                                    newPoints.append(points[ptj0]);
                                    tmpFaces.append(counter);
                                    counter++;
                                }

                                List<pointIndexHit> pih(1);
                                surfSearch.findLine
                                (
                                    Field<vector>({points[ptj0]}),
                                    Field<vector>({points[ptj1]}),
                                    pih
                                );
                                if (pih[0].hit())
                                {
                                    newPoints.append(pih[0].hitPoint());

                                    bool yetCut = false;

                                    forAll(cutPointsI, cpi)
                                    {
                                        if
                                        (
                                            Foam::mag
                                            (
                                                newPoints[counter]
                                              - newPoints[cutPointsI[cpi]]
                                            ) < nearTol_
                                        )
                                        {
                                            yetCut = true;
                                        }
                                    }

                                    if (!yetCut)
                                    {
                                        cutPointsI.append(counter);
                                    }

                                    tmpFaces.append(counter);
                                    counter++;
                                }
                                else
                                {
                                    FatalErrorInFunction
                                        << "no point hit by segment: " << nl
                                        << "point 0: " << points[ptj0] << nl
                                        << "point 1: " << points[ptj1] << endl;
                                }
                            }
                            // if both points outside do nothing
                        }

                        if (!tmpFaces.empty())
                        {
                            newFaces.append(tmpFaces);
                        }

                        tmpFaces.clear();
                    }
                }
            }
        }

        // Meaning it is not cut probably
        if (cutPointsI.size() < 2)
        {
            cutMaskPtr_()[ci] = 0.0;
            continue;
        }

        // reordering of cut points
        vector centre(Zero);
        forAll (cutPointsI, i)
        {
            centre += newPoints[cutPointsI[i]];
        }
        centre /= cutPointsI.size();
        vector r0(newPoints[cutPointsI[0]] - centre);
        vector nref = r0 ^ (newPoints[cutPointsI[1]] - centre);
        nref /= (Foam::mag(nref) + nearTol_);
        SortableList<scalar> thetas(cutPointsI.size());
        thetas[0] = 0;
        for (int i = 1; i < cutPointsI.size(); i++)
        {
            vector r1 = newPoints[cutPointsI[i]] - centre;
            scalar cos = r1 & r0 / Foam::mag(r0) / Foam::mag(r1);
            cos *= (1 - nearTol_); //adjust to prevent fpe
            vector cross = r0 ^ r1 / Foam::mag(r0) / Foam::mag(r1);
            scalar sin = cross & nref;
            thetas[i] = (sin > 0) ? Foam::acos(cos) : (2*M_PI-Foam::acos(cos));
        }
        thetas.sort();
        List<label> newCutPoints(cutPointsI.size());
        forAll (cutPointsI, i)
        {
            newCutPoints[i] = cutPointsI[thetas.indices()[i]];
        }
        newFaces.append(newCutPoints);

        cutPoints[ci] = newPoints;
        cutCellsFaces[ci] = newFaces;
    }
}


template<class Type>
Foam::scalar Foam::Error<Type>::cellFacePointError
(
    const bool& grad
) const
{
    const surfaceTypeField& fF = faceF();
    const pointTypeField& pF = pointF();

    scalar err = 0.0;

    const vectorField& C = mesh_.cellCentres();
    const vectorField& Cf = mesh_.faceCentres();
    const vectorField& P = mesh_.points();
    const cellList& cellFaces = mesh_.cells();
    const faceList& facePoints = mesh_.faces();

    List<vector> verts(4, Zero);
    List<Type> values(4, Zero);

    forAll(C, ci)
    {
        if (mask()[ci] < small || cutMask()[ci] > small)
        {
            continue;
        }

        verts[0] = C[ci];
        values[0] = F_[ci];

        const labelList& facesI = cellFaces[ci];

        for (int j = 0; j < facesI.size(); j++)
        {
            verts[1] = Cf[facesI[j]];
            label patchi =
                mesh_.isInternalFace(facesI[j]) ?
                -1 : mesh_.boundaryMesh().whichPatch(facesI[j]);

            if (patchi > -1 && isA<emptyFvPatch>(mesh_.boundary()[patchi]))
            {
                continue;
            }
            else if (patchi > -1)
            {
                //// Avoid dealing with boundaries
                //continue;
                label starti = mesh_.boundary()[patchi].start();
                values[1] = fF.boundaryField()[patchi][facesI[j] - starti];
            }
            else
            {
                values[1] = fF[facesI[j]];
            }

            const labelList& pointsJ = facePoints[facesI[j]];

            for (int k = 0; k < pointsJ.size(); k++)
            {
                int l = k == 0 ? pointsJ.size() - 1 : k - 1;
                verts[2] = P[pointsJ[k]];
                values[2] = pF[pointsJ[k]];
                verts[3] = P[pointsJ[l]];
                values[3] = pF[pointsJ[l]];

                if (grad)
                {
                    ScalarFunction linErr = [&] (const vector& xx)
                    {
                        return magSqr
                        (
                            gradFEx_(xx)
                          - gradTetLinear(verts, values)
                        );
                    };

                    err += tetQuadrature(verts, linErr);
                }
                else
                {
                    ScalarFunction linErr = [&] (const vector& xx)
                    {
                        return magSqr(FEx_(xx) - tetLinear(verts, values, xx));
                    };

                    err += tetQuadrature(verts, linErr);
                }
            }
        }
    }

    reduce(err, sumOp<scalar>());

    return Foam::sqrt(err);
}


template<class Type>
Foam::scalar Foam::Error<Type>::dualError
(
    const bool& grad
) const
{
    scalar err = 0.0;

    const vectorField& P = mesh_.points();
    const labelListList& pointCells = mesh_.pointCells();
    const vectorField& C = mesh_.cellCentres();

    forAll(P, pi)
    {
        const labelList& cellsI = pointCells[pi];

        bool masked = true;

        forAll(cellsI, cj)
        {
            if
            (
                mask()[cellsI[cj]] > small
             || cutMask()[cellsI[cj]] < small
            )
            {
                masked = false;
                break;
            }
        }

        if (masked)
        {
            continue;
        }

        if (mesh_.nGeometricD() == 2 && cellsI.size() == 3)
        {
            List<vector> verts(3);
            List<Type> values(3);

            forAll(cellsI, cj)
            {
                verts[cj] = C[cellsI[cj]];
                values[cj] = F_[cellsI[cj]];
            }

            if (grad)
            {
                ScalarFunction linErr = [&] (const vector& xx)
                {
                    return magSqr
                    (
                        gradFEx_(xx)
                      - gradTriLinear(verts, values)
                    );
                };

                err += triQuadrature(verts, linErr);
            }
            else
            {
                ScalarFunction linErr = [&] (const vector& xx)
                {
                    return magSqr(FEx_(xx) - triLinear(verts, values, xx));
                };

                err += triQuadrature(verts, linErr);
            }
        }
        else if (cellsI.size() == 4)
        {
            List<vector> verts(4);
            List<Type> values(4);

            forAll(cellsI, cj)
            {
                verts[cj] = C[cellsI[cj]];
                values[cj] = F_[cellsI[cj]];
            }

            if (grad)
            {
                ScalarFunction linErr = [&] (const vector& xx)
                {
                    return magSqr
                    (
                        gradFEx_(xx)
                      - gradTetLinear(verts, values)
                    );
                };

                err += tetQuadrature(verts, linErr);
            }
            else
            {
                ScalarFunction linErr = [&] (const vector& xx)
                {
                    return magSqr(FEx_(xx) - tetLinear(verts, values, xx));
                };

                err += tetQuadrature(verts, linErr);
            }
        }
    }

    reduce(err, sumOp<scalar>());

    return Foam::sqrt(err);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::Error<Type>::Error
(
    const fvMesh& mesh,
    const GeometricField<Type, fvPatchField, volMesh>& F,
    const typename Foam::Error<Type>::Function& FEx,
    const typename Foam::Error<Type>::GradFunction& gradFEx
)
:
    mesh_(mesh),
    F_(F),
    FEx_(FEx),
    gradFEx_(gradFEx),
    pointMeshPtr_(nullptr),
    volPointInterpPtr_(nullptr),
    faceFPtr_(nullptr),
    pointFPtr_(nullptr),
    h_(-1.0)
{}


template<class Type>
Foam::Error<Type>::Error
(
    const fvMesh& mesh,
    const GeometricField<Type, fvPatchField, volMesh>& F,
    const typename Foam::Error<Type>::Function& FEx
)
:
    mesh_(mesh),
    F_(F),
    FEx_(FEx),
    gradFEx_(GZero_),
    pointMeshPtr_(nullptr),
    volPointInterpPtr_(nullptr),
    faceFPtr_(nullptr),
    pointFPtr_(nullptr),
    h_(-1.0)
{}


template<class Type>
Foam::Error<Type>::Error
(
    const fvMesh& mesh,
    const GeometricField<Type, fvPatchField, volMesh>& F
)
:
    mesh_(mesh),
    F_(F),
    FEx_(Zero_),
    gradFEx_(GZero_),
    pointMeshPtr_(nullptr),
    volPointInterpPtr_(nullptr),
    faceFPtr_(nullptr),
    pointFPtr_(nullptr),
    h_(-1.0)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::Error<Type>::~Error()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Type>
const Foam::pointMesh& Foam::Error<Type>::pMesh() const
{
    if (!pointMeshPtr_.valid())
    {
        makePointMesh();
    }

    return pointMeshPtr_();
}


template<class Type>
const Foam::volPointInterpolation&
Foam::Error<Type>::volPointInterp() const
{
    if (!volPointInterpPtr_.valid())
    {
        makeVolPointInterp();
    }

    return volPointInterpPtr_();
}


template<class Type>
const typename Foam::Error<Type>::surfaceTypeField&
Foam::Error<Type>::faceF() const
{
    if (!faceFPtr_.valid())
    {
        makeFaceField();
    }

    return faceFPtr_();
}


template<class Type>
const typename Foam::Error<Type>::pointTypeField&
Foam::Error<Type>::pointF() const
{
    if (!pointFPtr_.valid())
    {
        makePointField();
    }

    return pointFPtr_();
}


template<class Type>
Foam::scalar
Foam::Error<Type>::h() const
{
    if (h_ < 0)
    {
        const vectorField& points = mesh_.points();

        forAll(mesh_.C(), ci)
        {
            if (mask()[ci] > small)
            {
                const labelList& cellPoints = mesh_.cellPoints()[ci];

                int npoints = cellPoints.size();

                scalar bbMag = 0.0;

                for (int i = 0; i < npoints; i++)
                {
                    for (int j = i + 1; j < npoints; j++)
                    {
                        scalar dist =
                            mag
                            (
                                points[cellPoints[i]]
                              - points[cellPoints[j]]
                            );

                        if (dist > bbMag)
                        {
                            bbMag = dist;
                        }
                    }
                }

                if (bbMag > h_)
                {
                    h_ = bbMag;
                }
            }
        }
    }

    return h_;
}


template<class Type>
void Foam::Error<Type>::addBoxMask
(
    const vector& v0,
    const vector& v1,
    const bool flip
) const
{
    if (!maskPtr_.valid())
    {
        maskPtr_ =
        new volScalarField
        (
            volScalarField::New
            (
                "mask",
                mesh_,
                dimensionedScalar(dimless, 1.0)
            )()
        );
    }

    volScalarField& mask = *maskPtr_;
    const vectorField& C = mesh_.cellCentres();

    forAll(C, ci)
    {
        bool isInside = false;

        if
        (
            C[ci][0] > v0[0] && C[ci][1] > v0[1] && C[ci][2] > v0[2]
         && C[ci][0] < v1[0] && C[ci][1] < v1[1] && C[ci][2] < v1[2]
        )
        {
            isInside = true;
        }

        if (!flip && isInside)
        {
            mask[ci] = 0.0;
        }
        else if (flip && !isInside)
        {
            mask[ci] = 0.0;
        }
    }
}


template<class Type>
void Foam::Error<Type>::addSphericalMask
(
    const scalar& r,
    const vector c,
    const bool flip
) const
{
    if (!maskPtr_.valid())
    {
        maskPtr_ =
        new volScalarField
        (
            volScalarField::New
            (
                "mask",
                mesh_,
                dimensionedScalar(dimless, 1.0)
            )()
        );
    }

    volScalarField& mask = *maskPtr_;
    const vectorField& C = mesh_.cellCentres();
    scalarField R = mag(C - c);

    forAll(C, ci)
    {
        if (!flip && R[ci] > r)
        {
            mask[ci] = 0.0;
        }
        else if (flip && R[ci] < r)
        {
            mask[ci] = 0.0;
        }
    }
}


template<class Type>
void Foam::Error<Type>::addMask
(
    const volScalarField& m
) const
{
    if (maskPtr_.valid())
    {
        WarningInFunction
            << "Mask already added" << nl;

        *maskPtr_ *= m;
    }
    else
    {
        maskPtr_ = new volScalarField(m);
        maskPtr_->rename("mask");
    }
}


template<class Type>
const Foam::volScalarField&
Foam::Error<Type>::mask() const
{
    if (!maskPtr_.valid())
    {
        addSphericalMask(great);
    }

    return maskPtr_();
}


template<class Type>
void Foam::Error<Type>::addSurface
(
    const fileName& surfName,
    const bool internal
) const
{
    if (surfacePtr_.valid())
    {
        FatalErrorInFunction
            << "Surface already added" << nl;
    }

    surfacePtr_ = new triSurfaceMesh
        (
            IOobject
            (
                surfName,
                mesh_.time().constant(),
                "triSurface",
                mesh_,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

    surfaceSearchPtr_ = new triSurfaceSearch(surfacePtr_().surface());

    internal_ = internal;
}


template<class Type>
const Foam::triSurfaceMesh&
Foam::Error<Type>::surface() const
{
    if (!surfacePtr_.valid())
    {
        FatalErrorInFunction
            << "Add first the surface using addSurface(name)"
            << endl;
    }

    return surfacePtr_();
}


template<class Type>
const Foam::triSurfaceSearch&
Foam::Error<Type>::surfaceSearch() const
{
    if (!surfaceSearchPtr_.valid())
    {
        FatalErrorInFunction
            << "Add first the surface using addSurface(name)"
            << endl;
    }

    return surfaceSearchPtr_();
}

template<class Type>
const Foam::volScalarField&
Foam::Error<Type>::cutMask() const
{
    if (!cutMaskPtr_.valid())
    {
        cutMaskPtr_ =
        new volScalarField
        (
            volScalarField::New
            (
                "cutMask",
                mesh_,
                dimensionedScalar(dimless, 0.0)
            )()
        );

        if (surfacePtr_.valid())
        {
            performCellCut();
        }
    }

    return cutMaskPtr_();
}


template<class Type>
void Foam::Error<Type>::writeProjFieldEx() const
{
    volTypeField projFEx = 0.0*F_;

    forAll(projFEx, ci)
    {
        projFEx[ci] = cellProj(ci, FEx_);
    }

    projFEx.rename(F_.name() + "Ex");

    projFEx().write();
}


template<class Type>
void Foam::Error<Type>::writeProjGradEx() const
{
    tmp<volGradTypeField> tprojGradFEx
    (
        volGradTypeField::New
        (
            "grad(" + F_.name() + "Ex)",
            mesh_,
            dimensioned<GradType>(F_.dimensions(), Zero)
        )
    );

    volGradTypeField projGradFEx = tprojGradFEx.ref();

    forAll(projGradFEx, ci)
    {
        projGradFEx[ci] = cellProj(ci, gradFEx_);
    }

    projGradFEx.write();
}


template<class Type>
void Foam::Error<Type>::writePointField() const
{
    pointF().write();
}


template<class Type>
void Foam::Error<Type>::writePointFieldEx() const
{
    pointTypeField pointFEx(pointF());
    pointFEx.rename("point" + F_.name() + "Ex");

    const vectorField& P = mesh_.points();

    forAll(P, pi)
    {
        pointFEx[pi] = FEx_(P[pi]);
    }

    pointFEx.write();
}


template<class Type>
void Foam::Error<Type>::writePointGrad() const
{
    tmp<pointGradTypeField> tpointGradF =
    pointGradTypeField::New
    (
        "pointGrad" + F_.name(),
        pMesh(),
        dimensioned<GradType>(F_.dimensions()/dimLength, Zero)
    );

    pointGradTypeField& pointGradF = tpointGradF.ref();

    pointTypeField pointFEx(pointF());
    pointFEx.rename("point" + F_.name() + "Ex");

    const vectorField& P = mesh_.points();
    const labelListList& pointCells = mesh_.pointCells();
    const vectorField& C = mesh_.cellCentres();

    forAll(P, pi)
    {
        const labelList& cellsI = pointCells[pi];

        if (mesh_.nGeometricD() == 2 && cellsI.size() == 3)
        {
            List<vector> verts(3);
            List<Type> values(3);

            forAll(cellsI, cj)
            {
                verts[cj] = C[cellsI[cj]];
                values[cj] = F_[cellsI[cj]];
            }

            pointGradF[pi] = gradTriLinear(verts, values);
        }
        else if (cellsI.size() == 4)
        {
            List<vector> verts(4);
            List<Type> values(4);

            forAll(cellsI, cj)
            {
                verts[cj] = C[cellsI[cj]];
                values[cj] = F_[cellsI[cj]];
            }

            pointGradF[pi] = gradTetLinear(verts, values);
        }
    }

    pointGradF.write();
}


template<class Type>
void Foam::Error<Type>::writePointGradEx() const
{
    tmp<pointGradTypeField> tpointGradFEx =
    pointGradTypeField::New
    (
        "pointGrad" + F_.name() + "Ex",
        pMesh(),
        dimensioned<GradType>(F_.dimensions()/dimLength, Zero)
    );

    pointGradTypeField& pointGradFEx = tpointGradFEx.ref();

    const vectorField& P = mesh_.points();

    forAll(P, pi)
    {
        pointGradFEx[pi] = gradFEx_(P[pi]);
    }

    pointGradFEx.write();
}


template<class Type>
Foam::scalar Foam::Error<Type>::cellP0L2() const
{
    scalar err = 0.0;

    const vectorField& C = mesh_.cellCentres();
    const vectorField& Cf = mesh_.faceCentres();
    const vectorField& P = mesh_.points();
    const cellList& cellFaces = mesh_.cells();
    const faceList& facePoints = mesh_.faces();

    List<vector> verts(4, Zero);

    forAll(C, ci)
    {
        if (mask()[ci] < small || cutMask()[ci] > small)
        {
            continue;
        }

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

                ScalarFunction linErr = [&] (const vector& xx)
                {
                    return magSqr(FEx_(xx) - F_[ci]);
                };

                err += tetQuadrature(verts, linErr);
            }
        }
    }

    reduce(err, sumOp<scalar>());

    return Foam::sqrt(err);
}


template<class Type>
Foam::scalar Foam::Error<Type>::projL2() const
{
    scalar err = 0.0;

    forAll(F_, ci)
    {
        err +=
            mask()[ci]*(1.0 - cutMask()[ci])
           *mesh_.cellVolumes()[ci]
           *magSqr(F_[ci] - cellProj(ci, FEx_));
    }

    reduce(err, sumOp<scalar>());

    return Foam::sqrt(err);
}


template<class Type>
Foam::scalar Foam::Error<Type>::projH1() const
{
    tmp<volGradTypeField> tgradF = fvc::grad(F_);
    const volGradTypeField& gradF = tgradF.ref();

    scalar err = 0.0;

    forAll(F_, ci)
    {
        err +=
            mask()[ci]*(1.0 - cutMask()[ci])
           *mesh_.cellVolumes()[ci]
           *magSqr(gradF[ci] - cellProj(ci, gradFEx_));
    }

    reduce(err, sumOp<scalar>());

    return Foam::sqrt(err);
}


template<class Type>
Foam::scalar Foam::Error<Type>::projDiscreteH1
(
    const bool bd
) const
{
    volTypeField projFEx = 0.0*F_;

    forAll(projFEx, ci)
    {
        projFEx[ci] = cellProj(ci, FEx_);
    }

    const surfaceScalarField& magSf = mesh_.magSf();
    const surfaceScalarField& deltas = mesh_.deltaCoeffs();

    const labelList& own = mesh_.owner();
    const labelList& nei = mesh_.neighbour();

    //tmp<surfaceTypeField> tsnF = fvc::snGrad(F_);
    //const surfaceTypeField& snF = tsnF();
    //tmp<surfaceTypeField> tsnFEx = fvc::snGrad(projFEx_);
    //const surfaceTypeField& snFEx = tsnFEx();

    scalar err = 0.0;
    for (int i = 0; i < mesh_.nInternalFaces(); i++)
    {
        if
        (
            (mask()[own[i]] < small || cutMask()[own[i]] > small)
         || (mask()[nei[i]] < small || cutMask()[nei[i]] > small)
        )
        {
            continue;
        }

        err += magSf[i]*deltas[i]
              *magSqr
               (
                   F_[own[i]] - F_[nei[i]]
                 - projFEx[own[i]] + projFEx[nei[i]]
               );
    }

    // Boundary error
    if (bd)
    {
        forAll(mesh_.boundary(), patchI)
        {
            if (!isA<emptyFvPatch>(mesh_.boundary()[patchI]))
            {
                const fvPatch& p = mesh_.boundary()[patchI];
                const label& start = p.start();

                forAll(p.Cf(), fi)
                {
                    if
                    (
                        mask()[own[start + fi]] < small
                     || cutMask()[own[start + fi]] > small
                    )
                    {
                        continue;
                    }

                    err +=
                        magSf.boundaryField()[patchI][fi]
                       *deltas.boundaryField()[patchI][fi]
                       *magSqr
                       (
                           F_[p.faceCells()[fi]]
                         - F_.boundaryField()[patchI][fi]
                         - projFEx[p.faceCells()[fi]]
                         + projFEx.boundaryField()[patchI][fi]
                       );
                }
            }
        }
    }

    reduce(err, sumOp<scalar>());

    return Foam::sqrt(err);
}


template<class Type>
Foam::scalar Foam::Error<Type>::dualStar() const
{
    const surfaceVectorField& Sf = mesh_.Sf();
    const surfaceScalarField& magSf = mesh_.magSf();
    const surfaceScalarField& deltas = mesh_.deltaCoeffs();

    const labelList& own = mesh_.owner();
    const labelList& nei = mesh_.neighbour();

    scalar err = 0.0;
    for (int i = 0; i < mesh_.nInternalFaces(); i++)
    {
        if
        (
            (mask()[own[i]] < small || cutMask()[own[i]] > small)
         || (mask()[nei[i]] < small || cutMask()[nei[i]] > small)
        )
        {
            continue;
        }

        ScalarFunction linErr = [&] (const vector& xx)
        {
            return magSqr
            (
                (F_[own[i]] - F_[nei[i]])*deltas[i]
              + (gradFEx_(xx) & Sf[i]/magSf[i])
            );
        };

        err += faceProj(i, linErr)/deltas[i];
    }

    reduce(err, sumOp<scalar>());

    return Foam::sqrt(err);
}


template<class Type>
Foam::scalar Foam::Error<Type>::dualL2() const
{
    return dualError(false);
}


template<class Type>
Foam::scalar Foam::Error<Type>::dualH1() const
{
    return dualError(true);
}


template<class Type>
Foam::scalar Foam::Error<Type>::linearL2() const
{
    return cellFacePointError(false);
}


template<class Type>
Foam::scalar Foam::Error<Type>::linearH1() const
{
    return cellFacePointError(true);
}


// ************************************************************************* //
