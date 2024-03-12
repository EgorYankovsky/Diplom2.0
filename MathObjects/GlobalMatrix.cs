using DataStructs;
namespace MathObjects;

public class GlobalMatrix : Matrix
{
    #region Компоненты матрицы.
    // Массив индексов по i.
    internal readonly int[]? _ig;

    // Массив индексов по j.
    internal readonly List<int> _jg;

    // Диагональные элементы.
    internal readonly double[]? _diag;

    // Элементы матрицы нижнего треугольника .
    internal double[]? _al;

    // Элементы матрицы верхнего треугольника.
    internal double[]? _au;
    #endregion

    #region Переопределение математических операций
    public static GlobalVector operator *(GlobalMatrix _gm, GlobalVector _gv) => MathOperations.Multiply(_gm, _gv);

    public static GlobalMatrix operator *(double a, GlobalMatrix gm) => MathOperations.Multiply(a, gm);

    public static GlobalMatrix operator +(GlobalMatrix gm1, GlobalMatrix gm2) => MathOperations.Sum(gm1, gm2);
    #endregion

    public int Size
    {
        get
        {
            if (_diag is null) throw new Exception("_diag is null");
            return _diag.Length;
        }
    }

    public override double this[int i, int j]
    {
        get 
        {
            if (i > _diag.Length || j > _diag.Length)
                throw new Exception("Index ran out of matrix.");

            switch (i - j)
            {
                case 0: return _diag[i];
                case < 0: return ReturanValueAU(j, i);
                case > 0: return ReturanValueAL(i, j);
            }
        }
    }

    private double ReturanValueAL(int i, int j)
    {
        for (int ii = 0; ii < _ig[i + 1] - _ig[i]; ii++)
            if (_jg[_ig[i] + ii] == j)
                return _al[_ig[i] + ii];
        return 0.0D;  
    }

    private double ReturanValueAU(int i, int j)
    {
        for (int ii = 0; ii < _ig[i + 1] - _ig[i]; ii++)
            if (_jg[_ig[i] + ii] == j)
                return _au[_ig[i] + ii];
        return 0.0D;  
    }   

    public bool CheckPortrait(GlobalMatrix gm)
    {
        if (Size != gm.Size) return false;
        if (_ig.Length != gm._ig.Length) return false;
        if (_jg.Count != gm._jg.Count) return false;
        
        for (int i = 0; i < _ig.Length; i++)
            if (_ig[i] != gm._ig[i])
                return false;
        for (int i = 0; i < _jg.Count; i++)
            if (_jg[i] != gm._jg[i])
                return false;
        return true;
    }

    # region Конструкторы на любой вкус и цвет.
    public GlobalMatrix(int[] ig, List<int> jg, double[] diag, double[] al, double[] au)
    {
        _ig = ig;
        _jg = jg;
        _diag = diag;
        _al = al;
        _au = au; 
    }

    public GlobalMatrix(GlobalMatrix gm)
    {
        _ig = gm._ig;
        _jg = gm._jg;
        _diag = (double[])gm._diag.Clone();
        _al = (double[])gm._al.Clone();
        _au = (double[])gm._au.Clone();
    }

    public GlobalMatrix(int arrOfPntLen)
    {        
        _jg = new();
        _ig = new int[arrOfPntLen + 1];
        _diag = new double[arrOfPntLen];
    }

    public GlobalMatrix(ArrayOfElems _arrOfElms, ArrayOfPoints _arrOfPnt, ArrayOfBorders _arrOfBord, List<double> mu0, double koef = 0.0)
    {
        _jg = new();
        _ig = new int[_arrOfPnt.Length + 1];
        _diag = new double[_arrOfPnt.Length];

        // Рабочий массив.
        List<List<int>> arr = new();

        // ! Дерьмодристный момент.
        for(int i = 0; i < _arrOfPnt.Length; i++)
            arr.Add(new List<int>());

        foreach (var _elem in _arrOfElms)
            foreach (var point in _elem)
                foreach (var pnt in _elem)
                    if (pnt < point && Array.IndexOf(arr[point].ToArray(), pnt) == -1)
                    {
                        arr[point].Add(pnt);
                        arr[point].Sort();
                    }

        _ig[0] = 0;
        for (int i = 0; i < _arrOfPnt.Length; i++)
        {
            _ig[i + 1] = _ig[i] + arr[i].Count;
            _jg.AddRange(arr[i]);
        }

        _al = new double[_jg.Count];
        _au = new double[_jg.Count];
        foreach (var _elem in _arrOfElms)
        {
            Add(_elem, _arrOfPnt, mu0[0], koef); //! Здесь необходимо подбирать mu0.
        }
        ConsiderBoundaryConditions(_arrOfBord, _arrOfPnt);    
    }
    #endregion

    private void ConsiderBoundaryConditions(ArrayOfBorders arrBrd, ArrayOfPoints arrPt)
    {
        if (_ig is null) throw new Exception("_ig is null.");
        if (_al is null) throw new Exception("_al is null.");
        if (_diag is null) throw new Exception("_diag is null.");
        if (_au is null) throw new Exception("_au is null.");
        
        foreach (var border in arrBrd)
        {
            switch (border[0])
            {
                case 1:
                {
                    for (int i = 2; i < 4; i++)
                    {
                        for (int j = _ig[border[i]]; j < _ig[border[i] + 1]; j++)
                            _al[j] = 0;
                        _diag[border[i]] = 1;
                        for (int j = 0; j < _jg.Count; j++)
                            if (_jg[j] == border[i])
                                _au[j] = 0;
                    }
                    break;
                }
                case 2: break; // du / dn = 0
                case 3: throw new ArgumentException("Пока нет возможности учитывать КУ 3-го рода");

            }
        }
    }

    private void Add(List<int> elem, ArrayOfPoints arrPt, double mu0, double koef)
    {
        LocalMatrix lm = new(elem, arrPt, mu0, koef);

        if (_diag is null) throw new Exception("_diag isn't initialized");
        if (_ig is null) throw new Exception("_ig isn't initialized");
        if (_au is null) throw new Exception("_au isn't initialized");
        if (_al is null) throw new Exception("_au isn't initialized");
        
        
        int ii = 0;
        foreach (var i in elem)
        {
            int jj = 0;
            foreach (var j in elem)
            {
                int ind = 0;
                double val = 0.0D;
                switch(i - j)
                {
                    case 0:
                        val = lm[ii, jj];
                        _diag[i] += lm[ii, jj];
                        break;
                    case < 0:
                        ind = _ig[j];
                        for (; ind <= _ig[j + 1] - 1; ind++)
                            if (_jg[ind] == i) break;
                        val = lm[ii, jj];
                        _au[ind] += lm[ii, jj];
                        break;
                    case > 0:
                        ind = _ig[i];
                        for (; ind <= _ig[i + 1] - 1; ind++)
                            if (_jg[ind] == j) break;
                        val = lm[ii, jj];
                        _al[ind] += lm[ii, jj];
                        break;
                }
                jj++;
            }
            ii++;
        }
    }

}