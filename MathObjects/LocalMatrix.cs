using DataStructs;
using static Functions.BasisFunctions2D;
using Solution;

namespace MathObjects;


public class LocalMatrix : Matrix
{
    private TypeOfMatrixM _typeOfMatrixM;

    private readonly double _lambda;

    private readonly double _gamma;


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
                TypeOfMatrixM.Mr =>  _gamma * (_Mr[i % 2, j % 2] * _Mz[i / 2, j / 2]),
                TypeOfMatrixM.Mrr => _lambda * (_Gr[i % 2, j % 2] * _Mz[i / 2, j / 2] + _Mr[i % 2, j % 2] * _Gz[i / 2, j / 2]) +
                                     _gamma * (_Mrr[i % 2, j % 2] * _Mz[i / 2, j / 2]),
                _ => throw new Exception("Unexpected matrix"),
            };
        }
        set{}
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

    public LocalMatrix(double lambda, double rk, double hz, double hr)
    {
        _rk = rk;
        _hr = hr;
        _hz = hz;
        double _d = _rk / _hr;
        _lambda = lambda;
        _gamma = _lambda;
        _Mr1 = new double[2,2] {{ (1 + _d) * (1 + _d), -1.0 * _d * (1 + _d)},
                                {-1.0 * _d * (1 + _d),              _d * _d}};
        _typeOfMatrixM = TypeOfMatrixM.Mrr;
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
    }

    public LocalMatrix(List<int> elem, ArrayOfPoints arrPt, TypeOfMatrixM typeOfMatrixM, double lambda = 0.0D, double gamma = 0.0D)
    {
        _typeOfMatrixM = typeOfMatrixM;
        _rk = arrPt[elem[0]].R;
        _hr = arrPt[elem[1]].R - arrPt[elem[0]].R;
        _hz = arrPt[elem[2]].Z - arrPt[elem[0]].Z;
        double _d = _rk / _hr;
        _lambda = lambda;
        _gamma = gamma;
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
                Console.Write($"{this[i, j]:E3} ");
                matr[i, j] = (1.0D / _mu0) * (_Gr[i % 2, j % 2] * _Mz[i / 2, j / 2] + _Mr[i % 2, j % 2] * _Gz[i / 2, j / 2]) + 
                (1.0D / _mu0) * (_Mrr[i % 2, j % 2] * _Mz[i / 2, j / 2]);
            }   
            Console.WriteLine();
        }
        Console.WriteLine();
        */
    }
}



public class LocalMatrixNum : Matrix
{
    public override double this[int i, int j]
    {
        get
        {
            if (_r0 <= 0.0D && 0.0D <= _r1)
                throw new ArgumentException("Промежуток интегрирования содержит 0");

            double dRidRj(double r) => (i % 2, j % 2) switch
            {
                (0, 0)           => dR1(_r1, _r0, r) * dR1(_r1, _r0, r),
                (0, 1) or (1, 0) => dR1(_r1, _r0, r) * dR2(_r1, _r0, r),
                (1, 1)           => dR2(_r1, _r0, r) * dR2(_r1, _r0, r),
                _ => 0.0D
            };

            double RiRj(double r) => (i % 2, j % 2) switch
            {
                (0, 0)           => R1(_r1, _r0, r) * R1(_r1, _r0, r),
                (0, 1) or (1, 0) => R1(_r1, _r0, r) * R2(_r1, _r0, r),
                (1, 1)           => R2(_r1, _r0, r) * R2(_r1, _r0, r),
                _ => 0.0D
            };

            double dZidZj(double z) => (i / 2, j / 2) switch
            {
                (0, 0)           => dZ1(_z1, _z0, z) * dZ1(_z1, _z0, z),
                (0, 1) or (1, 0) => dZ1(_z1, _z0, z) * dZ2(_z1, _z0, z),
                (1, 1)           => dZ2(_z1, _z0, z) * dZ2(_z1, _z0, z),
                _ => 0.0D
            };

            double ZiZj(double z) => (i / 2, j / 2) switch
            {
                (0, 0)           => Z1(_z1, _z0, z) * Z1(_z1, _z0, z),
                (0, 1) or (1, 0) => Z1(_z1, _z0, z) * Z2(_z1, _z0, z),
                (1, 1)           => Z2(_z1, _z0, z) * Z2(_z1, _z0, z),
                _ => 0.0D
            };


            var Gr = Gauss.Compute(_r0, _r1, (double r) => r * dRidRj(r));
            var Mz = Gauss.Compute(_z0, _z1, ZiZj);

            var Mr = Gauss.Compute(_r0, _r1, (double r) => r * RiRj(r));
            var Gz = Gauss.Compute(_z0, _z1, dZidZj);

            var Mrr = Gauss.Compute(_r0, _r1, (double r) => 1.0D / r * RiRj(r));
            
            return _typeOfMatrixM == TypeOfMatrixM.Mrr ? _lambda * (Gr * Mz + Mr * Gz + Mrr * Mz) : _gamma * Mr * Mz;
        }
        set{}
    }    

    private TypeOfMatrixM _typeOfMatrixM;

    private readonly double _lambda;

    private readonly double _gamma;

    private readonly double _r0;

    private readonly double _r1;

    private readonly double _z0;

    private readonly double _z1;

    public LocalMatrixNum(double lambda, double r0, double r1, double z0, double z1)
    {
        _lambda = lambda;
        _r0 = r0;
        _r1 = r1;
        _z0 = z0;
        _z1 = z1;
    }
    
    public LocalMatrixNum(List<int> elem, ArrayOfPoints arrPt, TypeOfMatrixM typeOfMatrixM, double lambda = 1.0D, double gamma = 1.0D)
    {
        _typeOfMatrixM = typeOfMatrixM;
        _r0 = arrPt[elem[0]].R;
        _r1 = arrPt[elem[1]].R;
        _z0 = arrPt[elem[0]].Z;
        _z1 = arrPt[elem[2]].Z;
        _lambda = lambda;
        _gamma = gamma;
    }
}