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

Class
    Foam::Error

Description

SourceFiles
    ErrorI.H
    Error.C

\*---------------------------------------------------------------------------*/

#ifndef Error_H
#define Error_H

#include <functional>
#include "fvMesh.H"
#include "List.H"
#include "fvcGrad.H"
#include "volFields.H"
#include "volPointInterpolation.H"
#include "pointMesh.H"
#include "surfaceMesh.H"
#include "triSurfaceMesh.H"
#include "triSurfaceTools.H"
#include "triSurfaceSearch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Type>
class Error;

/*---------------------------------------------------------------------------*\
                         Class Error Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class Error
{
public:

    // Public typedefs

        typedef typename outerProduct<vector, Type>::type GradType;

        typedef GeometricField<Type, fvPatchField, volMesh> volTypeField;

        typedef GeometricField<Type, fvsPatchField, surfaceMesh>
        surfaceTypeField;

        typedef GeometricField<Type, pointPatchField, pointMesh> pointTypeField;

        typedef GeometricField<GradType, fvPatchField, volMesh>
        volGradTypeField;

        typedef GeometricField<GradType, pointPatchField, pointMesh>
        pointGradTypeField;

        typedef std::function<scalar(const vector&)> ScalarFunction;

        typedef std::function<Type(const vector&)> Function;

        typedef std::function<GradType(const vector&)> GradFunction;


    // Static data

        static Function Zero_;

        static GradFunction GZero_;


private:

    // Private Data

        //- Mesh
        const fvMesh& mesh_;

        //- Field
        const volTypeField& F_;

        //- Analytic function
        const Function& FEx_;

        //- Analytic gradient
        const GradFunction& gradFEx_;


protected:

    // Protected demand-driven data

        //- Nodal interpolation of field
        mutable autoPtr<pointMesh> pointMeshPtr_;

        //- Nodal interpolation of field
        mutable autoPtr<volPointInterpolation> volPointInterpPtr_;

        //- Face interpolation of field
        mutable autoPtr<surfaceTypeField> faceFPtr_;

        //- Nodal interpolation of field
        mutable autoPtr<pointTypeField> pointFPtr_;

        //- Volume mask field
        mutable autoPtr<volScalarField> maskPtr_;

        //- Surface mesh
        mutable autoPtr<triSurfaceMesh> surfacePtr_;

        //- Consider internal or external
        mutable bool internal_;

        //- Surface mesh search
        mutable autoPtr<triSurfaceSearch> surfaceSearchPtr_;

        //- Cut cells mask field
        mutable autoPtr<volScalarField> cutMaskPtr_;

        // List of points (counting intersections) for each cut cell
        mutable autoPtr<List<List<vector> > > cutPointsPtr_;

        // List of points for each face for each cut cell
        mutable autoPtr<List<List<List<label> > > > cutCellsFacesPtr_;

        //- Max h
        mutable scalar h_;


    // Protected Member Functions

        //- Generate point mesh
        void makePointMesh() const;

        //- Generate volPoint interpolator
        void makeVolPointInterp() const;

        //- Generate face interpolated field
        void makeFaceField() const;

        //- Generate point interpolated field
        void makePointField() const;

        //- Generate point linearly interpolated field
        void makeLinearPointField() const;

        //- Cut cells with surface
        void performCellCut() const;

        //- Cell-face-point interpolation error
        scalar cellFacePointError(const bool&) const;

        //- Error on dual meshes
        scalar dualError(const bool&) const;


public:

    TypeName("Error");

    // Constructors

        //- Construct from components
        Error
        (
            const fvMesh&,
            const GeometricField<Type, fvPatchField, volMesh>&,
            const Function&,
            const GradFunction&
        );

        //- Construct from components
        Error
        (
            const fvMesh&,
            const GeometricField<Type, fvPatchField, volMesh>&,
            const Function&
        );

        //- Construct from components
        Error
        (
            const fvMesh&,
            const GeometricField<Type, fvPatchField, volMesh>&
        );

        //- Copy constructor
        Error(const Error<Type>&) = delete;


    //- Destructor
    virtual ~Error();


    // Member Functions

    // Access

        //- Nodal interpolation of field
        const pointMesh& pMesh() const;

        //- Nodal interpolation of field
        const volPointInterpolation& volPointInterp() const;

        //- Face interpolation of field
        const surfaceTypeField& faceF() const;

        //- Nodal interpolation of field
        const pointTypeField& pointF() const;

        //- Mesh max h
        scalar h() const;

        //- Volume mask
        void addBoxMask
        (
            const vector&,
            const vector&,
            const bool flip = false
        ) const;

        //- Volume mask
        void addSphericalMask
        (
            const scalar&,
            const vector c = vector(0.0, 0.0, 0.0),
            const bool flip = false
        ) const;

        //- Volume mask
        void addMask(const volScalarField&) const;

        //- Return mask
        const volScalarField& mask() const;

        //- Surface
        void addSurface
        (
            const fileName&,
            const bool internal = false
        ) const;

        //- Return surface mesh
        const Foam::triSurfaceMesh& surface() const;

        //- Return surface mesh search
        const Foam::triSurfaceSearch& surfaceSearch() const;

        //- Surface
        const volScalarField& cutMask() const;


    // Output

        //- Write point field
        void writeProjFieldEx() const;

        //- Write point field
        void writeProjGradEx() const;

        //- Write point field
        void writePointField() const;

        //- Write point exact field
        void writePointFieldEx() const;

        //- Write point gradient of the field
        void writePointGrad() const;

        //- Write point gradient of the exact field
        void writePointGradEx() const;


    // Error computation

        //- Cell-wise projection
        template<class func>
        inline typename func::result_type
        faceProj(const label&, const func&) const;

        //- Cell-wise projection
        template<class func>
        inline typename func::result_type
        cellProj(const label&, const func&) const;

        //- Piecewise constant L2 norm
        scalar cellP0L2() const;

        //- Piecewise constant L2 norm
        scalar projL2() const;

        //- Piecewise constant H1 seminorm
        scalar projH1() const;

        //- Discrete piecewise constant H1 seminorm
        scalar projDiscreteH1(const bool bd = true) const;

        //- Piecewise linear Star norm
        scalar dualStar() const;

        //- Piecewise linear L2 norm
        //  on dual mesh of a Delaunay triangulation
        scalar dualL2() const;

        //- Piecewise linear H1 seminorm
        //  on dual mesh of a Delaunay triangulation
        scalar dualH1() const;

        //- Piecewise linear L2 norm
        scalar linearL2() const;

        //- Piecewise linear H1 seminorm
        scalar linearH1() const;


    // Tetrahedron relative functions

        //- Check if it is a triangle
        inline void triCheck
        (
            const List<vector>&
        ) const;


        //- Check if it is a tetrahedron
        inline void tetCheck
        (
            const List<vector>&
        ) const;


        //- Compute polyhedron barycentre
        inline Foam::vector barycentre
        (
            const List<vector>&
        ) const;


        //- Compute triangle area
        inline Foam::scalar triArea
        (
            const List<vector>&
        ) const;


        //- Compute tetrahedron volume
        inline Foam::scalar tetVolume
        (
            const List<vector>&
        ) const;


        //- Compute linear coefficients for a triangle
        inline Foam::List<Type> triLinearCoeffs
        (
            const List<vector>&,
            const List<Type>&
        ) const;


        //- Compute linear coefficients for a tetrahedron
        inline Foam::List<Type> tetLinearCoeffs
        (
            const List<vector>&,
            const List<Type>&
        ) const;


        //- Interpolate a linear function in a point
        inline Type triLinear
        (
            const List<vector>& ,
            const List<Type>&,
            const vector&
        ) const;


        //- Interpolate a linear function in a point
        inline Type tetLinear
        (
            const List<vector>& ,
            const List<Type>&,
            const vector&
        ) const;


        //- Compute gradient of a linear function
        inline GradType gradTriLinear
        (
            const List<vector>& ,
            const List<Type>&
        ) const;


        //- Compute gradient of a linear function
        inline GradType gradTetLinear
        (
            const List<vector>& ,
            const List<Type>&
        ) const;


        //- Quadrature for a linear function
        template<class Type2>
        inline Type2 triLinearQuadrature
        (
            const List<vector>& ,
            const List<Type2>&
        ) const;


        //- Quadrature for a linear function
        template<class Type2>
        inline Type2 tetLinearQuadrature
        (
            const List<vector>& ,
            const List<Type2>&
        ) const;


        //- Quadrature for a lambda function
        template<class func>
        inline typename func::result_type triQuadrature
        (
            const List<vector>&,
            const func&
        ) const;


        //- Quadrature for a lambda function
        template<class func>
        inline typename func::result_type tetQuadrature
        (
            const List<vector>&,
            const func&
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ErrorI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "Error.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
