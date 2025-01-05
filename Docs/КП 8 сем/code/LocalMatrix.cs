using DataStructs;

namespace MathObjects;

public class LocalMatrix : Matrix
{
    private readonly double _lambda;

    private readonly double _rk;

    private readonly double _hr;

    private readonly double _gamma;

    private readonly double _hz;

    private readonly double _d;

    
    public double this[int i, int j]
    {
        get
        {
            if (i > 3 || j > 3) throw new IndexOutOfRangeException("Local matrix error.");
                return _lambda * ((_rk / _hr + 0.5) * _G[i % 2, j % 2] * _hz * _Mz[i / 2, j / 2] + 
                       (Math.Log(1 + 1 / _d) * _Mr1[i % 2, j % 2] - _d * _G[i % 2, j % 2] + _Mr2[i % 2, j % 2]) * _G[i / 2, j / 2] / _hz) +
                       _gamma * (_hz * _Mz[i / 2, j / 2] * (Math.Log(1 + 1 / _d) * _Mr1[i % 2, j % 2] - _d * _G[i % 2, j % 2] + _Mr2[i % 2, j % 2]));  
        }
    }


/*
    public double this[int i, int j]
    {
        get
        {
            if (i > 3 || j > 3) throw new IndexOutOfRangeException("Local matrix error.");
                return _lambda * (_G[i % 2, j % 2] / _hr * _hz * _Mz[i / 2, j / 2] + 
                                  _hr * _M[i % 2, j % 2] * _G[i / 2, j / 2] / _hz) +
                       _gamma * (_hz * _Mz[i / 2, j / 2] * _hr * _M[i % 2, j % 2]);  
        }
    }*/

    private readonly double[,] _G =  {{ 1.0, -1.0},
                                      {-1.0,  1.0}};

    private readonly double[,] _M = {{0.3333333333333333, 0.16666666666666666},
                                     {0.16666666666666666, 0.3333333333333333}};

    private readonly double[,] _Mz = {{0.3333333333333333, 0.1666666666666666},
                                     {0.1666666666666666, 0.3333333333333333}};

    private readonly double[,] _Mr1;

    private readonly double[,] _Mr2 = {{-1.5, 0.5},
                                        {0.5, 0.5}};



    public LocalMatrix(double lambda, double gamma, double rk, double hz, double hr)
    {
        _lambda = lambda;
        _gamma = gamma;
        _rk = rk;
        _hr = hr;
        _hz = hz;
    }

    public LocalMatrix(List<int> elem, ArrayOfPoints arrPt, double lambda = 1, double gamma = 1)
    {
        _rk = arrPt[elem[0]].R;
        _hr = arrPt[elem[1]].R - arrPt[elem[0]].R;
        _hz = arrPt[elem[2]].Z - arrPt[elem[0]].Z;
        _d = _rk / _hr;
        _lambda = lambda;
        _gamma = gamma;
        _Mr1 = new double[2,2] {{(1 + _d) * (1 + _d), -_d * (1 + _d)},
                                {     -_d * (1 + _d),        _d * _d}};
    }
}