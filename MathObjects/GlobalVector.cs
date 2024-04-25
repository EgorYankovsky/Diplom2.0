using System.Globalization;
using Functions;
using DataStructs;

namespace MathObjects;

public class GlobalVector : Vector
{
    #region Переопределение математических операций.
    public static GlobalVector operator -(GlobalVector gv1, GlobalVector gv2) => MathOperations.Diff(gv1, gv2);

    public static GlobalVector operator +(GlobalVector gv1, GlobalVector gv2) => MathOperations.Sum(gv1, gv2);

    public static GlobalVector operator *(double a, GlobalVector gv) => MathOperations.Multiply(a, gv);
    
    public static double operator *(GlobalVector gv1, GlobalVector gv2) => MathOperations.Multiply(gv1, gv2);
    #endregion

    public double Norma()
    {
        if (_values == null) throw new Exception("_values is null!");
        double ans = 0.0D;
        foreach (var v in _values)
            ans += v * v;
        return Math.Sqrt(ans);
    }

    public GlobalVector(GlobalVector gv)
    {
        _values = gv._values;
    }

    // ? Нужен ли public?
    public GlobalVector(int size)
    {
        _values = new double[size];
    }

    // ? Удалить после.
    public GlobalVector()
    {

    }

    public GlobalVector(double[] arr)
    {
        _values = arr;
    }

    public override string ToString()
    {
        if (_values is null) throw new Exception("_values is null");
        string result = "";
        foreach (var item in _values)
            result += $"{item.ToString("E15", CultureInfo.InvariantCulture)}\n";
        return result;
    }

    public void ConsiderBoundaryConditions(ArrayOfBorders arrBrd, ArrayOfPoints arrp, double t)
    {
        if (_values is null) throw new Exception("_values is null");
        foreach (var border in arrBrd)
            switch (border[0])
            {
                // КУ - I-го рода
                case 1:
                    switch (border[1])
                    {
                        case 1:
                        for (int i = 2; i < 4; i++)
                            _values[border[i]] = Function.U1_1(arrp[border[i]], t);
                        break;
                        
                        case 2:
                        for (int i = 2; i < 4; i++)
                            _values[border[i]] = Function.U1_2(arrp[border[i]], t);
                        break;

                        case 3:
                        for (int i = 2; i < 4; i++)
                            _values[border[i]] = Function.U1_3(arrp[border[i]], t);
                        break;
                        
                        case 4:
                        for (int i = 2; i < 4; i++)
                            _values[border[i]] = Function.U1_4(arrp[border[i]], t);
                        break;
                        
                        default:
                        throw new Exception("No such bord");
                    }
                    break;
                // КУ - II-го рода
                case 2:
                    for (int i = 2; i < 4; i++)
                        _values[border[i]] += 0.0D;
                    break;
                // КУ - III-го рода
                case 3: throw new ArgumentException("Пока нет возможности учитывать КУ III-го рода");
            }
    }

    public GlobalVector(ArrayOfPoints arrPt)
    {
        _values = new double[arrPt.Length];
    }
}
