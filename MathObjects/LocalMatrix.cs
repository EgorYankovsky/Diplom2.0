using DataStructs;

namespace MathObjects;


public class LocalMatrix : Matrix
{
    private TypeOfMatrixM _typeOfMatrixM;
    private readonly double _mu0;

    private readonly double _sigma;

    private readonly double _rk;

    private readonly double _hr;

    private readonly double _hz;

    public override double this[int i, int j]
    {
        get
        {
            if (i > 3 || j > 3) throw new IndexOutOfRangeException("Local matrix error.");
            return _typeOfMatrixM switch
            {
                TypeOfMatrixM.Mr =>  _sigma * (_Mr[i % 2, j % 2] * _Mz[i / 2, j / 2]),
                TypeOfMatrixM.Mrr => (1.0D / _mu0) * (_Gr[i % 2, j % 2] * _Mz[i / 2, j / 2] + _Mr[i % 2, j % 2] * _Gz[i / 2, j / 2]) +
                                     (1.0D / _mu0) * (_Mrr[i % 2, j % 2] * _Mz[i / 2, j / 2]),
                _ => throw new Exception("Unexpected matrix"),
            };
        }
    }

    private readonly double[,] _G =  {{ 1.0, -1.0},
                                      {-1.0,  1.0}};

    private readonly double[,] _Mz = {{2.0, 1.0}, 
                                      {1.0, 2.0}}; 

    private readonly double[,] _Mr1;

    private readonly double[,] _M1r = {{2.0, 1.0}, 
                                       {1.0, 2.0}}; 

    private readonly double[,] _M2r = {{1.0, 1.0}, 
                                       {1.0, 3.0}}; 

    private readonly double[,] _Mr2 = {{-3.0D, 1.0D},
                                       { 1.0D, 1.0D}};

    private readonly double[,] _Gr = new double[2, 2];
    private readonly double[,] _Mr = new double[2, 2];
    private readonly double[,] _Gz = new double[2, 2];

    private readonly double[,] _Mrr = new double[2, 2];

    double[,] matr = new double[4, 4];

    public LocalMatrix(List<int> elem, ArrayOfPoints arrPt, TypeOfMatrixM typeOfMatrixM, double mu0i = 1.0D, double sigmai = 1.0D)
    {
        _typeOfMatrixM = typeOfMatrixM;
        _rk = arrPt[elem[0]].R;
        _hr = arrPt[elem[1]].R - arrPt[elem[0]].R;
        _hz = arrPt[elem[2]].Z - arrPt[elem[0]].Z;
        double _d = _rk / _hr;
        _mu0 = mu0i;
        _sigma = sigmai;
        _Mr1 = new double[2,2] {{ (1 + _d) * (1 + _d), -1.0 * _d * (1 + _d)},
                                {-1.0 * _d * (1 + _d),              _d * _d}};

        for (int i = 0; i < 2; i++)
        { 
            for (int j = 0; j < 2; j++)
            {
                _Gr[i, j] = ((_rk + _hr / 2.0D) / _hr) * _G[i, j];
                _Mr[i, j] = (_hr / 6.0D) * (_rk * _M1r[i, j] + (_hr / 2.0D) * _M2r[i, j]);
                _Mrr[i, j] = Math.Log(1.0D + 1.0D / _d) * _Mr1[i, j] - _d * _G[i, j] + 0.5 * _Mr2[i, j];
                _Gz[i, j] = _G[i, j] / _hz;
                _Mz[i, j] = (_hz / 6.0D) * _Mz[i, j];
            }   
        }

        /*
        for (int i = 0; i < 4; i++)
        { 
            for (int j = 0; j < 4; j++)
            {
                matr[i, j] = (1.0D / _mu0) * (_Gr[i % 2, j % 2] * _Mz[i / 2, j / 2] + _Mr[i % 2, j % 2] * _Gz[i / 2, j / 2]) + 
                (1.0D / _mu0) * (_Mrr[i % 2, j % 2] * _Mz[i / 2, j / 2]);
            }   
        }
        */
    }
}