// P1 interpolation on triangles
Foam::scalar triP1
(
    const List<vector>& verts,
    const List<scalar>& values,
    const vector& xx
)
{
    const scalar& x0 = verts[0][0];
    const scalar& x1 = verts[1][0];
    const scalar& x2 = verts[2][0];
    const scalar& y0 = verts[0][1];
    const scalar& y1 = verts[1][1];
    const scalar& y2 = verts[2][1];
    const scalar& v0 = values[0];
    const scalar& v1 = values[1];
    const scalar& v2 = values[2];

    scalar a = (v0*y1 - v1*y0 - v0*y2 + v2*y0 + v1*y2 - v2*y1)/(x0*y1 - x1*y0 - x0*y2 + x2*y0 + x1*y2 - x2*y1);
    scalar b = -(v0*x1 - v1*x0 - v0*x2 + v2*x0 + v1*x2 - v2*x1)/(x0*y1 - x1*y0 - x0*y2 + x2*y0 + x1*y2 - x2*y1);
    scalar c = (v0*x1*y2 - v0*x2*y1 - v1*x0*y2 + v1*x2*y0 + v2*x0*y1 - v2*x1*y0)/(x0*y1 - x1*y0 - x0*y2 + x2*y0 + x1*y2 - x2*y1);

    return (a*xx[0] + b*xx[1] + c);
}


Foam::vector triP1
(
    const List<vector>& verts,
    const List<vector>& values,
    const vector& xx
)
{
    const scalar& x0 = verts[0][0];
    const scalar& x1 = verts[1][0];
    const scalar& x2 = verts[2][0];
    const scalar& y0 = verts[0][1];
    const scalar& y1 = verts[1][1];
    const scalar& y2 = verts[2][1];
    const vector& v0 = values[0];
    const vector& v1 = values[1];
    const vector& v2 = values[2];

    vector a = (v0*y1 - v1*y0 - v0*y2 + v2*y0 + v1*y2 - v2*y1)/(x0*y1 - x1*y0 - x0*y2 + x2*y0 + x1*y2 - x2*y1);
    vector b = -(v0*x1 - v1*x0 - v0*x2 + v2*x0 + v1*x2 - v2*x1)/(x0*y1 - x1*y0 - x0*y2 + x2*y0 + x1*y2 - x2*y1);
    vector c = (v0*x1*y2 - v0*x2*y1 - v1*x0*y2 + v1*x2*y0 + v2*x0*y1 - v2*x1*y0)/(x0*y1 - x1*y0 - x0*y2 + x2*y0 + x1*y2 - x2*y1);

    return (a*xx[0] + b*xx[1] + c);
}


Foam::vector triGradP1
(
    const List<vector>& verts,
    const List<scalar>& values
)
{
    const scalar& x0 = verts[0][0];
    const scalar& x1 = verts[1][0];
    const scalar& x2 = verts[2][0];
    const scalar& y0 = verts[0][1];
    const scalar& y1 = verts[1][1];
    const scalar& y2 = verts[2][1];
    const scalar& v0 = values[0];
    const scalar& v1 = values[1];
    const scalar& v2 = values[2];

    scalar a = (v0*y1 - v1*y0 - v0*y2 + v2*y0 + v1*y2 - v2*y1)/(x0*y1 - x1*y0 - x0*y2 + x2*y0 + x1*y2 - x2*y1);
    scalar b = -(v0*x1 - v1*x0 - v0*x2 + v2*x0 + v1*x2 - v2*x1)/(x0*y1 - x1*y0 - x0*y2 + x2*y0 + x1*y2 - x2*y1);

    return vector(a, b, 0.0);
}


Foam::tensor triGradP1
(
    const List<vector>& verts,
    const List<vector>& values
)
{
    const scalar& x0 = verts[0][0];
    const scalar& x1 = verts[1][0];
    const scalar& x2 = verts[2][0];
    const scalar& y0 = verts[0][1];
    const scalar& y1 = verts[1][1];
    const scalar& y2 = verts[2][1];
    const vector& v0 = values[0];
    const vector& v1 = values[1];
    const vector& v2 = values[2];

    vector a = (v0*y1 - v1*y0 - v0*y2 + v2*y0 + v1*y2 - v2*y1)/(x0*y1 - x1*y0 - x0*y2 + x2*y0 + x1*y2 - x2*y1);
    vector b = -(v0*x1 - v1*x0 - v0*x2 + v2*x0 + v1*x2 - v2*x1)/(x0*y1 - x1*y0 - x0*y2 + x2*y0 + x1*y2 - x2*y1);
    vector c = Zero;

    return tensor(a, b, c);
}


// P1 interpolation on tetrahedra
Foam::scalar tetP1
(
    const List<vector>& verts,
    const List<scalar>& values,
    const vector& xx
)
{
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
    const scalar& v0 = values[0];
    const scalar& v1 = values[1];
    const scalar& v2 = values[2];
    const scalar& v3 = values[3];

    scalar a = (v0*y1*z2 - v0*y2*z1 - v1*y0*z2 + v1*y2*z0 + v2*y0*z1 - v2*y1*z0 - v0*y1*z3 + v0*y3*z1 + v1*y0*z3 - v1*y3*z0 - v3*y0*z1 + v3*y1*z0 + v0*y2*z3 - v0*y3*z2 - v2*y0*z3 + v2*y3*z0 + v3*y0*z2 - v3*y2*z0 - v1*y2*z3 + v1*y3*z2 + v2*y1*z3 - v2*y3*z1 - v3*y1*z2 + v3*y2*z1)/(x0*y1*z2 - x0*y2*z1 - x1*y0*z2 + x1*y2*z0 + x2*y0*z1 - x2*y1*z0 - x0*y1*z3 + x0*y3*z1 + x1*y0*z3 - x1*y3*z0 - x3*y0*z1 + x3*y1*z0 + x0*y2*z3 - x0*y3*z2 - x2*y0*z3 + x2*y3*z0 + x3*y0*z2 - x3*y2*z0 - x1*y2*z3 + x1*y3*z2 + x2*y1*z3 - x2*y3*z1 - x3*y1*z2 + x3*y2*z1);
    scalar b = -(v0*x1*z2 - v0*x2*z1 - v1*x0*z2 + v1*x2*z0 + v2*x0*z1 - v2*x1*z0 - v0*x1*z3 + v0*x3*z1 + v1*x0*z3 - v1*x3*z0 - v3*x0*z1 + v3*x1*z0 + v0*x2*z3 - v0*x3*z2 - v2*x0*z3 + v2*x3*z0 + v3*x0*z2 - v3*x2*z0 - v1*x2*z3 + v1*x3*z2 + v2*x1*z3 - v2*x3*z1 - v3*x1*z2 + v3*x2*z1)/(x0*y1*z2 - x0*y2*z1 - x1*y0*z2 + x1*y2*z0 + x2*y0*z1 - x2*y1*z0 - x0*y1*z3 + x0*y3*z1 + x1*y0*z3 - x1*y3*z0 - x3*y0*z1 + x3*y1*z0 + x0*y2*z3 - x0*y3*z2 - x2*y0*z3 + x2*y3*z0 + x3*y0*z2 - x3*y2*z0 - x1*y2*z3 + x1*y3*z2 + x2*y1*z3 - x2*y3*z1 - x3*y1*z2 + x3*y2*z1);
    scalar c = (v0*x1*y2 - v0*x2*y1 - v1*x0*y2 + v1*x2*y0 + v2*x0*y1 - v2*x1*y0 - v0*x1*y3 + v0*x3*y1 + v1*x0*y3 - v1*x3*y0 - v3*x0*y1 + v3*x1*y0 + v0*x2*y3 - v0*x3*y2 - v2*x0*y3 + v2*x3*y0 + v3*x0*y2 - v3*x2*y0 - v1*x2*y3 + v1*x3*y2 + v2*x1*y3 - v2*x3*y1 - v3*x1*y2 + v3*x2*y1)/(x0*y1*z2 - x0*y2*z1 - x1*y0*z2 + x1*y2*z0 + x2*y0*z1 - x2*y1*z0 - x0*y1*z3 + x0*y3*z1 + x1*y0*z3 - x1*y3*z0 - x3*y0*z1 + x3*y1*z0 + x0*y2*z3 - x0*y3*z2 - x2*y0*z3 + x2*y3*z0 + x3*y0*z2 - x3*y2*z0 - x1*y2*z3 + x1*y3*z2 + x2*y1*z3 - x2*y3*z1 - x3*y1*z2 + x3*y2*z1);
    scalar d = -(v0*x1*y2*z3 - v0*x1*y3*z2 - v0*x2*y1*z3 + v0*x2*y3*z1 + v0*x3*y1*z2 - v0*x3*y2*z1 - v1*x0*y2*z3 + v1*x0*y3*z2 + v1*x2*y0*z3 - v1*x2*y3*z0 - v1*x3*y0*z2 + v1*x3*y2*z0 + v2*x0*y1*z3 - v2*x0*y3*z1 - v2*x1*y0*z3 + v2*x1*y3*z0 + v2*x3*y0*z1 - v2*x3*y1*z0 - v3*x0*y1*z2 + v3*x0*y2*z1 + v3*x1*y0*z2 - v3*x1*y2*z0 - v3*x2*y0*z1 + v3*x2*y1*z0)/(x0*y1*z2 - x0*y2*z1 - x1*y0*z2 + x1*y2*z0 + x2*y0*z1 - x2*y1*z0 - x0*y1*z3 + x0*y3*z1 + x1*y0*z3 - x1*y3*z0 - x3*y0*z1 + x3*y1*z0 + x0*y2*z3 - x0*y3*z2 - x2*y0*z3 + x2*y3*z0 + x3*y0*z2 - x3*y2*z0 - x1*y2*z3 + x1*y3*z2 + x2*y1*z3 - x2*y3*z1 - x3*y1*z2 + x3*y2*z1);

    return (a*xx[0] + b*xx[1] + c*xx[2] + d);
}


Foam::vector tetP1
(
    const List<vector>& verts,
    const List<vector>& values,
    const vector& xx
)
{
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
    const vector& v0 = values[0];
    const vector& v1 = values[1];
    const vector& v2 = values[2];
    const vector& v3 = values[3];

    vector a = (v0*y1*z2 - v0*y2*z1 - v1*y0*z2 + v1*y2*z0 + v2*y0*z1 - v2*y1*z0 - v0*y1*z3 + v0*y3*z1 + v1*y0*z3 - v1*y3*z0 - v3*y0*z1 + v3*y1*z0 + v0*y2*z3 - v0*y3*z2 - v2*y0*z3 + v2*y3*z0 + v3*y0*z2 - v3*y2*z0 - v1*y2*z3 + v1*y3*z2 + v2*y1*z3 - v2*y3*z1 - v3*y1*z2 + v3*y2*z1)/(x0*y1*z2 - x0*y2*z1 - x1*y0*z2 + x1*y2*z0 + x2*y0*z1 - x2*y1*z0 - x0*y1*z3 + x0*y3*z1 + x1*y0*z3 - x1*y3*z0 - x3*y0*z1 + x3*y1*z0 + x0*y2*z3 - x0*y3*z2 - x2*y0*z3 + x2*y3*z0 + x3*y0*z2 - x3*y2*z0 - x1*y2*z3 + x1*y3*z2 + x2*y1*z3 - x2*y3*z1 - x3*y1*z2 + x3*y2*z1);
    vector b = -(v0*x1*z2 - v0*x2*z1 - v1*x0*z2 + v1*x2*z0 + v2*x0*z1 - v2*x1*z0 - v0*x1*z3 + v0*x3*z1 + v1*x0*z3 - v1*x3*z0 - v3*x0*z1 + v3*x1*z0 + v0*x2*z3 - v0*x3*z2 - v2*x0*z3 + v2*x3*z0 + v3*x0*z2 - v3*x2*z0 - v1*x2*z3 + v1*x3*z2 + v2*x1*z3 - v2*x3*z1 - v3*x1*z2 + v3*x2*z1)/(x0*y1*z2 - x0*y2*z1 - x1*y0*z2 + x1*y2*z0 + x2*y0*z1 - x2*y1*z0 - x0*y1*z3 + x0*y3*z1 + x1*y0*z3 - x1*y3*z0 - x3*y0*z1 + x3*y1*z0 + x0*y2*z3 - x0*y3*z2 - x2*y0*z3 + x2*y3*z0 + x3*y0*z2 - x3*y2*z0 - x1*y2*z3 + x1*y3*z2 + x2*y1*z3 - x2*y3*z1 - x3*y1*z2 + x3*y2*z1);
    vector c = (v0*x1*y2 - v0*x2*y1 - v1*x0*y2 + v1*x2*y0 + v2*x0*y1 - v2*x1*y0 - v0*x1*y3 + v0*x3*y1 + v1*x0*y3 - v1*x3*y0 - v3*x0*y1 + v3*x1*y0 + v0*x2*y3 - v0*x3*y2 - v2*x0*y3 + v2*x3*y0 + v3*x0*y2 - v3*x2*y0 - v1*x2*y3 + v1*x3*y2 + v2*x1*y3 - v2*x3*y1 - v3*x1*y2 + v3*x2*y1)/(x0*y1*z2 - x0*y2*z1 - x1*y0*z2 + x1*y2*z0 + x2*y0*z1 - x2*y1*z0 - x0*y1*z3 + x0*y3*z1 + x1*y0*z3 - x1*y3*z0 - x3*y0*z1 + x3*y1*z0 + x0*y2*z3 - x0*y3*z2 - x2*y0*z3 + x2*y3*z0 + x3*y0*z2 - x3*y2*z0 - x1*y2*z3 + x1*y3*z2 + x2*y1*z3 - x2*y3*z1 - x3*y1*z2 + x3*y2*z1);
    vector d = -(v0*x1*y2*z3 - v0*x1*y3*z2 - v0*x2*y1*z3 + v0*x2*y3*z1 + v0*x3*y1*z2 - v0*x3*y2*z1 - v1*x0*y2*z3 + v1*x0*y3*z2 + v1*x2*y0*z3 - v1*x2*y3*z0 - v1*x3*y0*z2 + v1*x3*y2*z0 + v2*x0*y1*z3 - v2*x0*y3*z1 - v2*x1*y0*z3 + v2*x1*y3*z0 + v2*x3*y0*z1 - v2*x3*y1*z0 - v3*x0*y1*z2 + v3*x0*y2*z1 + v3*x1*y0*z2 - v3*x1*y2*z0 - v3*x2*y0*z1 + v3*x2*y1*z0)/(x0*y1*z2 - x0*y2*z1 - x1*y0*z2 + x1*y2*z0 + x2*y0*z1 - x2*y1*z0 - x0*y1*z3 + x0*y3*z1 + x1*y0*z3 - x1*y3*z0 - x3*y0*z1 + x3*y1*z0 + x0*y2*z3 - x0*y3*z2 - x2*y0*z3 + x2*y3*z0 + x3*y0*z2 - x3*y2*z0 - x1*y2*z3 + x1*y3*z2 + x2*y1*z3 - x2*y3*z1 - x3*y1*z2 + x3*y2*z1);

    return (a*xx[0] + b*xx[1] + c*xx[2] + d);
}


Foam::vector tetGradP1
(
    const List<vector>& verts,
    const List<scalar>& values
)
{
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
    const scalar& v0 = values[0];
    const scalar& v1 = values[1];
    const scalar& v2 = values[2];
    const scalar& v3 = values[3];

    scalar a = (v0*y1*z2 - v0*y2*z1 - v1*y0*z2 + v1*y2*z0 + v2*y0*z1 - v2*y1*z0 - v0*y1*z3 + v0*y3*z1 + v1*y0*z3 - v1*y3*z0 - v3*y0*z1 + v3*y1*z0 + v0*y2*z3 - v0*y3*z2 - v2*y0*z3 + v2*y3*z0 + v3*y0*z2 - v3*y2*z0 - v1*y2*z3 + v1*y3*z2 + v2*y1*z3 - v2*y3*z1 - v3*y1*z2 + v3*y2*z1)/(x0*y1*z2 - x0*y2*z1 - x1*y0*z2 + x1*y2*z0 + x2*y0*z1 - x2*y1*z0 - x0*y1*z3 + x0*y3*z1 + x1*y0*z3 - x1*y3*z0 - x3*y0*z1 + x3*y1*z0 + x0*y2*z3 - x0*y3*z2 - x2*y0*z3 + x2*y3*z0 + x3*y0*z2 - x3*y2*z0 - x1*y2*z3 + x1*y3*z2 + x2*y1*z3 - x2*y3*z1 - x3*y1*z2 + x3*y2*z1);
    scalar b = -(v0*x1*z2 - v0*x2*z1 - v1*x0*z2 + v1*x2*z0 + v2*x0*z1 - v2*x1*z0 - v0*x1*z3 + v0*x3*z1 + v1*x0*z3 - v1*x3*z0 - v3*x0*z1 + v3*x1*z0 + v0*x2*z3 - v0*x3*z2 - v2*x0*z3 + v2*x3*z0 + v3*x0*z2 - v3*x2*z0 - v1*x2*z3 + v1*x3*z2 + v2*x1*z3 - v2*x3*z1 - v3*x1*z2 + v3*x2*z1)/(x0*y1*z2 - x0*y2*z1 - x1*y0*z2 + x1*y2*z0 + x2*y0*z1 - x2*y1*z0 - x0*y1*z3 + x0*y3*z1 + x1*y0*z3 - x1*y3*z0 - x3*y0*z1 + x3*y1*z0 + x0*y2*z3 - x0*y3*z2 - x2*y0*z3 + x2*y3*z0 + x3*y0*z2 - x3*y2*z0 - x1*y2*z3 + x1*y3*z2 + x2*y1*z3 - x2*y3*z1 - x3*y1*z2 + x3*y2*z1);
    scalar c = (v0*x1*y2 - v0*x2*y1 - v1*x0*y2 + v1*x2*y0 + v2*x0*y1 - v2*x1*y0 - v0*x1*y3 + v0*x3*y1 + v1*x0*y3 - v1*x3*y0 - v3*x0*y1 + v3*x1*y0 + v0*x2*y3 - v0*x3*y2 - v2*x0*y3 + v2*x3*y0 + v3*x0*y2 - v3*x2*y0 - v1*x2*y3 + v1*x3*y2 + v2*x1*y3 - v2*x3*y1 - v3*x1*y2 + v3*x2*y1)/(x0*y1*z2 - x0*y2*z1 - x1*y0*z2 + x1*y2*z0 + x2*y0*z1 - x2*y1*z0 - x0*y1*z3 + x0*y3*z1 + x1*y0*z3 - x1*y3*z0 - x3*y0*z1 + x3*y1*z0 + x0*y2*z3 - x0*y3*z2 - x2*y0*z3 + x2*y3*z0 + x3*y0*z2 - x3*y2*z0 - x1*y2*z3 + x1*y3*z2 + x2*y1*z3 - x2*y3*z1 - x3*y1*z2 + x3*y2*z1);

    return vector(a, b, c);
}


Foam::tensor tetGradP1
(
    const List<vector>& verts,
    const List<vector>& values
)
{
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
    const vector& v0 = values[0];
    const vector& v1 = values[1];
    const vector& v2 = values[2];
    const vector& v3 = values[3];

    vector a = (v0*y1*z2 - v0*y2*z1 - v1*y0*z2 + v1*y2*z0 + v2*y0*z1 - v2*y1*z0 - v0*y1*z3 + v0*y3*z1 + v1*y0*z3 - v1*y3*z0 - v3*y0*z1 + v3*y1*z0 + v0*y2*z3 - v0*y3*z2 - v2*y0*z3 + v2*y3*z0 + v3*y0*z2 - v3*y2*z0 - v1*y2*z3 + v1*y3*z2 + v2*y1*z3 - v2*y3*z1 - v3*y1*z2 + v3*y2*z1)/(x0*y1*z2 - x0*y2*z1 - x1*y0*z2 + x1*y2*z0 + x2*y0*z1 - x2*y1*z0 - x0*y1*z3 + x0*y3*z1 + x1*y0*z3 - x1*y3*z0 - x3*y0*z1 + x3*y1*z0 + x0*y2*z3 - x0*y3*z2 - x2*y0*z3 + x2*y3*z0 + x3*y0*z2 - x3*y2*z0 - x1*y2*z3 + x1*y3*z2 + x2*y1*z3 - x2*y3*z1 - x3*y1*z2 + x3*y2*z1);
    vector b = -(v0*x1*z2 - v0*x2*z1 - v1*x0*z2 + v1*x2*z0 + v2*x0*z1 - v2*x1*z0 - v0*x1*z3 + v0*x3*z1 + v1*x0*z3 - v1*x3*z0 - v3*x0*z1 + v3*x1*z0 + v0*x2*z3 - v0*x3*z2 - v2*x0*z3 + v2*x3*z0 + v3*x0*z2 - v3*x2*z0 - v1*x2*z3 + v1*x3*z2 + v2*x1*z3 - v2*x3*z1 - v3*x1*z2 + v3*x2*z1)/(x0*y1*z2 - x0*y2*z1 - x1*y0*z2 + x1*y2*z0 + x2*y0*z1 - x2*y1*z0 - x0*y1*z3 + x0*y3*z1 + x1*y0*z3 - x1*y3*z0 - x3*y0*z1 + x3*y1*z0 + x0*y2*z3 - x0*y3*z2 - x2*y0*z3 + x2*y3*z0 + x3*y0*z2 - x3*y2*z0 - x1*y2*z3 + x1*y3*z2 + x2*y1*z3 - x2*y3*z1 - x3*y1*z2 + x3*y2*z1);
    vector c = (v0*x1*y2 - v0*x2*y1 - v1*x0*y2 + v1*x2*y0 + v2*x0*y1 - v2*x1*y0 - v0*x1*y3 + v0*x3*y1 + v1*x0*y3 - v1*x3*y0 - v3*x0*y1 + v3*x1*y0 + v0*x2*y3 - v0*x3*y2 - v2*x0*y3 + v2*x3*y0 + v3*x0*y2 - v3*x2*y0 - v1*x2*y3 + v1*x3*y2 + v2*x1*y3 - v2*x3*y1 - v3*x1*y2 + v3*x2*y1)/(x0*y1*z2 - x0*y2*z1 - x1*y0*z2 + x1*y2*z0 + x2*y0*z1 - x2*y1*z0 - x0*y1*z3 + x0*y3*z1 + x1*y0*z3 - x1*y3*z0 - x3*y0*z1 + x3*y1*z0 + x0*y2*z3 - x0*y3*z2 - x2*y0*z3 + x2*y3*z0 + x3*y0*z2 - x3*y2*z0 - x1*y2*z3 + x1*y3*z2 + x2*y1*z3 - x2*y3*z1 - x3*y1*z2 + x3*y2*z1);

    return tensor(a, b, c);
}


// Quadrature
Foam::vector polyBary(const vectorField& verts, const List<label>& labels)
{
    vector v1(0.0, 0.0, 0.0);
    forAll (labels, i)
    {
        v1 += verts[labels[i]];
    }
    return  v1/labels.size();
}


Foam::vector polyBary(const List<vector>& verts)
{
    return  Foam::sum(verts)/verts.size();
}


Foam::scalar triArea(const List<vector>& verts)
{
    if (verts.size() != 3)
    {
        InfoInFunction
            <<"Warning: triangle is not 3 vertices but "
            <<verts.size()<<endl;
    }

    vector a = verts[1] - verts[0];
    vector b = verts[2] - verts[0];

    return  Foam::mag(a ^ b)/2;
}


Foam::scalar tetVol(const List<vector>& verts)
{
    if (verts.size() != 4)
    {
        InfoInFunction
            <<"Warning: tetrahedron is not 4 vertices but "
            <<verts.size()<<endl;
    }

    vector a = verts[1] - verts[0];
    vector b = verts[2] - verts[0];
    vector c = verts[3] - verts[0];

    return  Foam::mag(a & (b ^ c))/6;
}


Foam::scalar triQuad
(
    const List<vector>& verts,
    const Function3 func3
)
{
    if (verts.size() != 3)
    {
        InfoInFunction
            <<"Warning: triangle is not 3 vertices but "
            <<verts.size()<<endl;
    }

    List<vector> nodes(6);
    scalar weight(1.0/6.0);
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

    scalar quadVal(0.0);
    for (int i = 0; i < 6; i++)
    {
        quadVal += (weight*func3(nodes[i]));
    }

    return triArea(verts)*quadVal;
}


Foam::vector triQuad
(
    const List<vector>& verts,
    const VFunction3 func3
)
{
    if (verts.size() != 3)
    {
        InfoInFunction
            <<"Warning: triangle is not 3 vertices but "
            <<verts.size()<<endl;
    }

    List<vector> nodes(6);
    scalar weight(1.0/6.0);
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

    vector quadVal(0.0, 0.0, 0.0);
    for (int i = 0; i < 6; i++)
    {
        quadVal += (weight*func3(nodes[i]));
    }

    return triArea(verts)*quadVal;
}


Foam::tensor triQuad
(
    const List<vector>& verts,
    const TFunction3 func3
)
{
    if (verts.size() != 3)
    {
        InfoInFunction
            <<"Warning: triangle is not 3 vertices but "
            <<verts.size()<<endl;
    }

    List<vector> nodes(6);
    scalar weight(1.0/6.0);
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

    tensor quadVal(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    for (int i = 0; i < 6; i++)
    {
        quadVal += (weight*func3(nodes[i]));
    }

    return triArea(verts)*quadVal;
}


Foam::scalar tetQuad
(
    const List<vector>& verts,
    const Function3 func3
)
{
    if (verts.size() != 4)
    {
        InfoInFunction
            <<"Warning: tetrahedron is not 4 vertices but "
            <<verts.size()<<endl;
    }

    List<vector> nodes(5);
    List<scalar> weights(5);
    weights[0] = -16.0;
    weights[1] = 9.0;
    weights[2] = 9.0;
    weights[3] = 9.0;
    weights[4] = 9.0;
    scalar factor(1.0/20.0);
    nodes[0] = polyBary(verts);
    nodes[1] = 0.5*verts[0] + 1.0/6.0*(verts[1] + verts[2] + verts[3]);
    nodes[2] = 0.5*verts[1] + 1.0/6.0*(verts[0] + verts[2] + verts[3]);
    nodes[3] = 0.5*verts[2] + 1.0/6.0*(verts[0] + verts[1] + verts[3]);
    nodes[4] = 0.5*verts[3] + 1.0/6.0*(verts[0] + verts[1] + verts[2]);

    scalar quadVal(0.0);
    for (int i = 0; i < 5; i++)
    {
        quadVal += (factor*weights[i]*func3(nodes[i]));
    }

    return tetVol(verts)*quadVal;
}


Foam::vector tetQuad
(
    const List<vector>& verts,
    const VFunction3 func3
)
{
    if (verts.size() != 4)
    {
        InfoInFunction
            <<"Warning: tetrahedron is not 4 vertices but "
            <<verts.size()<<endl;
    }

    List<vector> nodes(5);
    List<scalar> weights(5);
    weights[0] = -16.0;
    weights[1] = 9.0;
    weights[2] = 9.0;
    weights[3] = 9.0;
    weights[4] = 9.0;
    scalar factor(1.0/20.0);
    nodes[0] = polyBary(verts);
    nodes[1] = 0.5*verts[0] + 1.0/6.0*(verts[1] + verts[2] + verts[3]);
    nodes[2] = 0.5*verts[1] + 1.0/6.0*(verts[0] + verts[2] + verts[3]);
    nodes[3] = 0.5*verts[2] + 1.0/6.0*(verts[0] + verts[1] + verts[3]);
    nodes[4] = 0.5*verts[3] + 1.0/6.0*(verts[0] + verts[1] + verts[2]);

    vector quadVal(0.0, 0.0, 0.0);
    for (int i = 0; i < 5; i++)
    {
        quadVal += (factor*weights[i]*func3(nodes[i]));
    }

    return tetVol(verts)*quadVal;
}


Foam::tensor tetQuad
(
    const List<vector>& verts,
    const TFunction3 func3
)
{
    if (verts.size() != 4)
    {
        InfoInFunction
            <<"Warning: tetrahedron is not 4 vertices but "
            <<verts.size()<<endl;
    }

    List<vector> nodes(5);
    List<scalar> weights(5);
    weights[0] = -16.0;
    weights[1] = 9.0;
    weights[2] = 9.0;
    weights[3] = 9.0;
    weights[4] = 9.0;
    scalar factor(1.0/20.0);
    nodes[0] = polyBary(verts);
    nodes[1] = 0.5*verts[0] + 1.0/6.0*(verts[1] + verts[2] + verts[3]);
    nodes[2] = 0.5*verts[1] + 1.0/6.0*(verts[0] + verts[2] + verts[3]);
    nodes[3] = 0.5*verts[2] + 1.0/6.0*(verts[0] + verts[1] + verts[3]);
    nodes[4] = 0.5*verts[3] + 1.0/6.0*(verts[0] + verts[1] + verts[2]);

    tensor quadVal(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    for (int i = 0; i < 5; i++)
    {
        quadVal += (factor*weights[i]*func3(nodes[i]));
    }

    return tetVol(verts)*quadVal;
}


Foam::scalar tetP1Quad
(
    const List<vector>& verts,
    const List<scalar>& values
)
{
    if (verts.size() != 4)
    {
        InfoInFunction
            <<"Warning: tetrahedron is not 4 vertices but "
            <<verts.size()<<endl;
    }

    return tetVol(verts)*0.25
          *(values[0] + values[1]
          + values[2] + values[3]);
}


Foam::vector tetP1Quad
(
    const List<vector>& verts,
    const List<vector>& values
)
{
    if (verts.size() != 4)
    {
        InfoInFunction
            <<"Warning: tetrahedron is not 4 vertices but "
            <<verts.size()<<endl;
    }

    return tetVol(verts)*0.25
          *(values[0] + values[1]
          + values[2] + values[3]);
}


Foam::tensor tetP1Quad
(
    const List<vector>& verts,
    const List<tensor>& values
)
{
    if (verts.size() != 4)
    {
        InfoInFunction
            <<"Warning: tetrahedron is not 4 vertices but "
            <<verts.size()<<endl;
    }

    return tetVol(verts)*0.25
          *(values[0] + values[1]
          + values[2] + values[3]);
}


