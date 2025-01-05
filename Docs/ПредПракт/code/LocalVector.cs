using DataStructs;

namespace MathObjects;

public class LocalVector : Vector
{
    private readonly bool _isRingBoundaryInside;

    
    private double F(double r, double z)
    {
        if (_isRingBoundaryInside && z == 0.0)
        {
            double h = Math.Abs(10.0D - r);
            return 1.0D - h / Math.Abs(_hr);            
        }
        return 0.0D;
    }
    
    //private static double F(double r, double z) => z / r;

    private readonly double _r0;

    private readonly double _r1;

    private readonly double _z0;

    private readonly double _z1;

    private readonly double _hr;

    private readonly double _hz;

    private readonly double _d;

    private readonly double[,] _M2R = {{0.08333333333333333, 0.08333333333333333},
                                       {0.08333333333333333, 0.25}};    

    private readonly double[,] _Mz = {{0.3333333333333333, 0.1666666666666666},
                                     {0.1666666666666666, 0.3333333333333333}};

    private readonly double[,] _M = {{0.3333333333333333, 0.1666666666666666},
                                     {0.1666666666666666, 0.3333333333333333}};


    private readonly double[,] _Mr1;

    private readonly double[,] _G = {{1, -1}, 
                                     {-1, 1}};

    private readonly double[,] _Mr2 = {{-1.5, 0.5},
                                       { 0.5, 0.5}};



    public override double this[int i] => i switch
    {
        0 => _hz * _hr * (_Mz[0, 0] * (Math.Log(1 + 1 / _d) * _Mr1[0, 0] - _d * _G[0, 0] + _Mr2[0, 0]) * F(_r0, _z0) +
                    _Mz[0, 0] * (Math.Log(1 + 1 / _d) * _Mr1[0, 1] - _d * _G[0, 1] + _Mr2[0, 1]) * F(_r1, _z0) +
                    _Mz[0, 1] * (Math.Log(1 + 1 / _d) * _Mr1[0, 0] - _d * _G[0, 0] + _Mr2[0, 0]) * F(_r0, _z1) +
                    _Mz[0, 1] * (Math.Log(1 + 1 / _d) * _Mr1[0, 1] - _d * _G[0, 1] + _Mr2[0, 1]) * F(_r1, _z1)),

        1 => _hz * _hr * (_Mz[0, 0] * (Math.Log(1 + 1 / _d) * _Mr1[1, 0] - _d * _G[1, 0] + _Mr2[1, 0]) * F(_r0, _z0) +
                    _Mz[0, 0] * (Math.Log(1 + 1 / _d) * _Mr1[1, 1] - _d * _G[1, 1] + _Mr2[1, 1]) * F(_r1, _z0) +
                    _Mz[0, 1] * (Math.Log(1 + 1 / _d) * _Mr1[1, 0] - _d * _G[1, 0] + _Mr2[1, 0]) * F(_r0, _z1) +
                    _Mz[0, 1] * (Math.Log(1 + 1 / _d) * _Mr1[1, 1] - _d * _G[1, 1] + _Mr2[1, 1]) * F(_r1, _z1)),

        2 => _hz * _hr * (_Mz[1, 0] * (Math.Log(1 + 1 / _d) * _Mr1[0, 0] - _d * _G[0, 0] + _Mr2[0, 0]) * F(_r0, _z0) +
                    _Mz[1, 0] * (Math.Log(1 + 1 / _d) * _Mr1[0, 1] - _d * _G[0, 1] + _Mr2[0, 1]) * F(_r1, _z0) +
                    _Mz[1, 1] * (Math.Log(1 + 1 / _d) * _Mr1[0, 0] - _d * _G[0, 0] + _Mr2[0, 0]) * F(_r0, _z1) +
                    _Mz[1, 1] * (Math.Log(1 + 1 / _d) * _Mr1[0, 1] - _d * _G[0, 1] + _Mr2[0, 1]) * F(_r1, _z1)),

        3 => _hz * _hr * (_Mz[1, 0] * (Math.Log(1 + 1 / _d) * _Mr1[1, 0] - _d * _G[1, 0] + _Mr2[1, 0]) * F(_r0, _z0) +
                    _Mz[1, 0] * (Math.Log(1 + 1 / _d) * _Mr1[1, 1] - _d * _G[1, 1] + _Mr2[1, 1]) * F(_r1, _z0) +
                    _Mz[1, 1] * (Math.Log(1 + 1 / _d) * _Mr1[1, 0] - _d * _G[1, 0] + _Mr2[1, 0]) * F(_r0, _z1) +
                    _Mz[1, 1] * (Math.Log(1 + 1 / _d) * _Mr1[1, 1] - _d * _G[1, 1] + _Mr2[1, 1]) * F(_r1, _z1)),
        _ => throw new IndexOutOfRangeException("Vector out of index"),
    };
/*
    public override double this[int i] => i switch
    {
        0 => _hz * _hr * (_Mz[0, 0] * _M[0, 0] * F(_r0, _z0) +
                    _Mz[0, 0] * _M[0, 1] * F(_r1, _z0) +
                    _Mz[0, 1] * _M[0, 0] * F(_r0, _z1) +
                    _Mz[0, 1] * _M[0, 1] * F(_r1, _z1)),

        1 => _hz * _hr * (_Mz[0, 0] * _M[1, 0] * F(_r0, _z0) +
                    _Mz[0, 0] * _M[1, 1] * F(_r1, _z0) +
                    _Mz[0, 1] * _M[1, 0] * F(_r0, _z1) +
                    _Mz[0, 1] * _M[1, 1] * F(_r1, _z1)),

        2 => _hz * _hr * (_Mz[1, 0] * _M[0, 0] * F(_r0, _z0) +
                    _Mz[1, 0] * _M[0, 1] * F(_r1, _z0) +
                    _Mz[1, 1] * _M[0, 0] * F(_r0, _z1) +
                    _Mz[1, 1] * _M[0, 1] * F(_r1, _z1)),

        3 => _hz * _hr * (_Mz[1, 0] * _M[1, 0] * F(_r0, _z0) +
                    _Mz[1, 0] * _M[1, 1] * F(_r1, _z0) +
                    _Mz[1, 1] * _M[1, 0] * F(_r0, _z1) +
                    _Mz[1, 1] * _M[1, 1] * F(_r1, _z1)),
        _ => throw new IndexOutOfRangeException("Vector out of index"),
    };
*/

    public LocalVector(List<int>? elem, ArrayOfPoints arrPt)
    {
        _r0 = arrPt[elem[0]].R;
        _r1 = arrPt[elem[1]].R;

        _isRingBoundaryInside = _r0 <= 10.0D && 10.0D <= _r1;
        
        _z0 = arrPt[elem[0]].Z;
        _z1 = arrPt[elem[2]].Z;
        _hr = _r1 - _r0;
        _hz = _z1 - _z0;
        _d = _r0 / _hr;
        _Mr1 = new double[2, 2] {{(1 + _d) * (1 + _d), -_d * (1 + _d)},
                                 {     -_d * (1 + _d),        _d * _d}};
    }

    public LocalVector(ArrayOfPoints arrPt, List<int>? arrBr)
    {
        Console.WriteLine("II bc commited :)");
    }

    public LocalVector(double r0, double r1, double z0, double z1)
    {
        _r0 = r0;
        _r1 = r1;
        _z0 = z0;
        _z1 = z1;
        _hr = _r1 - _r0;
        _hz = _z1 - _z0;
        _d = _r0 / _hr;
    } 
}