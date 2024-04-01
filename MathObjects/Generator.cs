using DataStructs;
using Functions;
namespace MathObjects;

public static class Generator
{
    public static void FillMatrixG(ref GlobalMatrix m, ArrayOfPoints arrPt, ArrayOfElems arrEl, ArrayOfBorders arrBd, TypeOfMatrixM typeOfMatrixM)
    {
        m._al = new double[m._jg.Count];
        m._au = new double[m._jg.Count];

        for (int i = 0; i < arrEl.Length; i++)
            Add(new LocalMatrix(arrEl[i], arrPt, typeOfMatrixM, arrEl.mu0i[i], 0.0), ref m, arrEl[i]);
        
        ConsiderBoundaryConditions(ref m, arrBd);
    }

    public static void FillMatrixM(ref GlobalMatrix m, ArrayOfPoints arrPt, ArrayOfElems arrEl, ArrayOfBorders arrBd, TypeOfMatrixM typeOfMatrixM)
    {
        m._al = new double[m._jg.Count];
        m._au = new double[m._jg.Count];

        for (int i = 0; i < arrEl.Length; i++)
            Add(new LocalMatrix(arrEl[i], arrPt, typeOfMatrixM, 0.0, arrEl.mu0i[i]), ref m, arrEl[i]);

        ConsiderBoundaryConditions(ref m, arrBd);
    }


    public static void BuildPortait(ref GlobalMatrix m, int arrPtLen, ArrayOfElems arrEl)
    {
        List<List<int>> arr = new();

        // ! Дерьмодристный момент.
        for(int i = 0; i < arrPtLen; i++)
            arr.Add(new List<int>());

        foreach (var _elem in arrEl)
            foreach (var point in _elem)
                foreach (var pnt in _elem)
                    if (pnt < point && Array.IndexOf(arr[point].ToArray(), pnt) == -1)
                    {
                        arr[point].Add(pnt);
                        arr[point].Sort();
                    }

        m._ig[0] = 0;
        for (int i = 0; i < arrPtLen; i++)
        {
            m._ig[i + 1] = m._ig[i] + arr[i].Count;
            m._jg.AddRange(arr[i]);
        }
    }

    public static void FillMatrix(ref GlobalMatrix m, ArrayOfPoints arrPt, ArrayOfElems arrEl, ArrayOfBorders arrBd, TypeOfMatrixM typeOfMatrixM)
    {    
        m._al = new double[m._jg.Count];
        m._au = new double[m._jg.Count];

        for (int i = 0; i < arrEl.Length; i++)
        {
            Add(new LocalMatrix(arrEl[i], arrPt, typeOfMatrixM, arrEl.mu0i[i], arrEl.mu0i[i]), ref m, arrEl[i]);
        }
        //ConsiderBoundaryConditions(ref m, arrBd);
    }

    private static void Add(LocalMatrix lm, ref GlobalMatrix gm, List<int> elem)
    {
        if (gm._diag is null) throw new Exception("_diag isn't initialized");
        if (gm._ig is null) throw new Exception("_ig isn't initialized");
        if (gm._au is null) throw new Exception("_au isn't initialized");
        if (gm._al is null) throw new Exception("_au isn't initialized");
        
        
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
                        gm._diag[i] += lm[ii, jj];
                        break;
                    case < 0:
                        ind = gm._ig[j];
                        for (; ind <= gm._ig[j + 1] - 1; ind++)
                            if (gm._jg[ind] == i) break;
                        val = lm[ii, jj];
                        gm._au[ind] += lm[ii, jj];
                        break;
                    case > 0:
                        ind = gm._ig[i];
                        for (; ind <= gm._ig[i + 1] - 1; ind++)
                            if (gm._jg[ind] == j) break;
                        val = lm[ii, jj];
                        gm._al[ind] += lm[ii, jj];
                        break;
                }
                jj++;
            }
            ii++;
        }
    }

    public static void ConsiderBoundaryConditions(ref GlobalMatrix m, ArrayOfBorders arrBd)
    {
        if (m._ig is null) throw new Exception("_ig is null.");
        if (m._al is null) throw new Exception("_al is null.");
        if (m._diag is null) throw new Exception("_diag is null.");
        if (m._au is null) throw new Exception("_au is null.");
        
        foreach (var border in arrBd)
        {
            switch (border[0])
            {
                case 1:
                {
                    for (int i = 2; i < 4; i++)
                    {
                        for (int j = m._ig[border[i]]; j < m._ig[border[i] + 1]; j++)
                            m._al[j] = 0;
                        m._diag[border[i]] = 1;
                        for (int j = 0; j < m._jg.Count; j++)
                            if (m._jg[j] == border[i])
                                m._au[j] = 0;
                    }
                    break;
                }
                case 2: break; // du / dn = 0
                case 3: throw new ArgumentException("Пока нет возможности учитывать КУ 3-го рода");

            }
        }
    }

    public static void FillVector(ref GlobalVector v, ArrayOfPoints arrPt, ArrayOfElems arrEl, ArrayOfBorders arrBd)
    {

    }

    public static void ConsiderBoundaryConditions(ref GlobalVector v, ArrayOfBorders arrBrd, ArrayOfPoints arrp, double t)
    {
        if (v is null) throw new Exception("_values is null");
        foreach (var border in arrBrd)
            switch (border[0])
            {
                // КУ - I-го рода
                case 1:
                    switch (border[1])
                    {
                        case 1:
                        for (int i = 2; i < 4; i++)
                            v[border[i]] = Function.U1_1(arrp[border[i]], t);
                        break;
                        
                        case 2:
                        for (int i = 2; i < 4; i++)
                            v[border[i]] = Function.U1_2(arrp[border[i]], t);
                        break;

                        case 3:
                        for (int i = 2; i < 4; i++)
                            v[border[i]] = Function.U1_3(arrp[border[i]], t);
                        break;
                        
                        case 4:
                        for (int i = 2; i < 4; i++)
                            v[border[i]] = Function.U1_4(arrp[border[i]], t);
                        break;
                        
                        default:
                        throw new Exception("No such bord");
                    }
                    break;
                // КУ - II-го рода
                case 2:
                    for (int i = 2; i < 4; i++)
                        v[border[i]] += 0.0D;
                    break;
                // КУ - III-го рода
                case 3: throw new ArgumentException("Пока нет возможности учитывать КУ III-го рода");
            }
    }

}