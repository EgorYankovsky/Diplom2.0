using System.Globalization;
using Functions;
using DataStructs;

namespace MathObjects;

public class GlobalVector : Vector
{
    #region Переопределение математических операций.
    public static GlobalVector operator -(GlobalVector gv1, GlobalVector gv2) => MathOperations.DiffPar(gv1, gv2);

    public static GlobalVector operator +(GlobalVector gv1, GlobalVector gv2) => MathOperations.SumPar(gv1, gv2);

    public static GlobalVector operator *(double a, GlobalVector gv) => MathOperations.MultiplyPar(a, gv);
    
    public static double operator *(GlobalVector gv1, GlobalVector gv2) => MathOperations.MultiplyPar(gv1, gv2);
    #endregion

    public double Norma()
    {
        if (_values == null) throw new Exception("_values is null!");
        double ans = 0.0D;
        foreach (var v in _values)
            ans += v * v;
        return Math.Sqrt(ans);
    }

    public GlobalVector(GlobalVector gv) => _values = (double[])gv._values.Clone();

    public GlobalVector(int size) => _values = new double[size];

    public GlobalVector(double[] arr) => _values = arr;

    public override string ToString()
    {
        if (_values is null) throw new Exception("_values is null");
        string result = "";
        foreach (var item in _values)
            result += $"{item.ToString("E15", CultureInfo.InvariantCulture)}\n";
        return result;
    }

    public GlobalVector(ArrayOfPoints arrPt)
    {
        _values = new double[arrPt.GetLength()];
    }
}
