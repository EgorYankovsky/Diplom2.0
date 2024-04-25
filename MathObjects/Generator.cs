using DataStructs;
using Functions;
using Functions;

namespace MathObjects;

public static class Generator
{
    public static void FillMatrixG(ref GlobalMatrix m, ArrayOfRibs arrRibs, ArrayOfElems arrEl)
    {
        m._al = new double[m._jg.Count];
        m._au = new double[m._jg.Count];

        using var sw = new StreamWriter("C:\\Users\\USER\\Desktop\\fgnf.txt");
        for (int i = 0; i < arrEl.Length; i++)
        {
            List<int> currElem = [arrEl[i][0], arrEl[i][3], arrEl[i][8], arrEl[i][11],
                                  arrEl[i][1], arrEl[i][2], arrEl[i][9], arrEl[i][10],
                                  arrEl[i][4], arrEl[i][5], arrEl[i][6], arrEl[i][7]];
            
            double hx = arrRibs[currElem[0]].Length;
            double hy = arrRibs[currElem[4]].Length;
            double hz = arrRibs[currElem[8]].Length;

            var lm = new LocalMatrixG3D(arrEl.mui[i], hx, hy, hz);
            Add(lm, ref m, currElem); 
        }
    }

    public static void FillMatrixM(ref GlobalMatrix m, ArrayOfRibs arrRibs, ArrayOfElems arrEl)
    {
        m._al = new double[m._jg.Count];
        m._au = new double[m._jg.Count];

        for (int i = 0; i < arrEl.Length; i++)
        {
            List<int> currElem = [arrEl[i][0], arrEl[i][3], arrEl[i][8], arrEl[i][11],
                                  arrEl[i][1], arrEl[i][2], arrEl[i][9], arrEl[i][10],
                                  arrEl[i][4], arrEl[i][5], arrEl[i][6], arrEl[i][7]];
            
            double hx = arrRibs[currElem[0]].Length;
            double hy = arrRibs[currElem[4]].Length;
            double hz = arrRibs[currElem[8]].Length;

            var lm = new LocalMatrixM3D(arrEl.mui[i], hx, hy, hz);
            Add(lm, ref m, currElem);
        }
    }

    public static void ConsiderBoundaryConditions(ref GlobalMatrix m, ref GlobalVector v, ArrayOfRibs arrRibs, ArrayOfBorders arrBrd, double t)
    {
        foreach (var border in arrBrd)
        {
            switch (border[0])
            {
                // КУ - I-го рода
                case 1:
                for (int i = 2; i < 6; i++)
                {
                    for (int j = m._ig[border[i]]; j < m._ig[border[i] + 1]; j++)
                        m._al[j] = 0.0D;
                    m._diag[border[i]] = 1.0D;
                    for (int j = 0; j < m._jg.Count; j++)
                        if (m._jg[j] == border[i])
                            m._au[j] = 0.0D;
                }

                for (int i = 2; i < 6; i++)
                {
                    var len = arrRibs[border[i]].Length;
                    var antinormal = ((arrRibs[border[i]].b.X - arrRibs[border[i]].a.X) / len, 
                                      (arrRibs[border[i]].b.Y - arrRibs[border[i]].a.Y) / len,
                                      (arrRibs[border[i]].b.Z - arrRibs[border[i]].a.Z) / len); 
                    var xm = 0.5D * (arrRibs[border[i]].b.X + arrRibs[border[i]].a.X);
                    var ym = 0.5D * (arrRibs[border[i]].b.Y + arrRibs[border[i]].a.Y);
                    var zm = 0.5D * (arrRibs[border[i]].b.Z + arrRibs[border[i]].a.Z);
                    
                    var f = Function.A(xm, ym, zm, t);
                    var q = antinormal.Item1 * f.Item1 + antinormal.Item2 * f.Item2 + antinormal.Item3 * f.Item3;
                    v[border[i]] = q;
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

        
        foreach (var border in arrBrd)
        {
            for (int i = 2; i < 6; i++)
            {
                for (int j = 0; j < m.Size; j++)
                {
                    if (border[i] == j)
                        continue;
                    else
                    {
                        var k = m[j, border[i]];
                        var f = -1.0D * k * v[border[i]];
                        v[j] += f;
                        m[j, border[i]] = 0.0D;
                    }
                }
            }
        }
    }

    public static void FillVector3D(ref GlobalVector v, ArrayOfRibs arrRibs, ArrayOfElems arrEl, double t)
    {
        for (int i = 0; i < arrEl.Length; i++)
        {
            List<int> currElem = [arrEl[i][0], arrEl[i][3], arrEl[i][8], arrEl[i][11],
                                  arrEl[i][1], arrEl[i][2], arrEl[i][9], arrEl[i][10],
                                  arrEl[i][4], arrEl[i][5], arrEl[i][6], arrEl[i][7]];
  
            var lv = new LocalVector3D(arrRibs[currElem[0]].a.X, arrRibs[currElem[0]].b.X,
                                       arrRibs[currElem[4]].a.Y, arrRibs[currElem[4]].b.Y,
                                       arrRibs[currElem[8]].a.Z, arrRibs[currElem[8]].b.Z, t);
            Add(lv, ref v, currElem);
        }
    }

    private static void Add(LocalVector3D lv, ref GlobalVector v, List<int> elem)
    {
        for (int i = 0; i < 12; i++)
            v[elem[i]] += lv[i];
/*
        v[elem[0]] += lv[0];
        v[elem[1]] += lv[4];
        v[elem[2]] += lv[5];
        v[elem[3]] += lv[1];

        v[elem[4]] += lv[8];
        v[elem[5]] += lv[9];
        v[elem[6]] += lv[10];
        v[elem[7]] += lv[11];
        
        v[elem[8]] += lv[2];
        v[elem[9]] += lv[6];
        v[elem[10]] += lv[7];
        v[elem[11]] += lv[3];
  */
    }

    public static void BuildPortait(ref GlobalMatrix m, int arrPtLen, ArrayOfElems arrEl)
    {
        List<List<int>> arr = [];
        List<List<int>> arr = [];

        // ! Дерьмодристный момент.
        for(int i = 0; i < arrPtLen; i++)
            arr.Add([]);

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
        m._al = new double[m._jg.Count];
        m._au = new double[m._jg.Count];
    }

    public static void FillMatrix(ref GlobalMatrix m, ArrayOfPoints arrPt, ArrayOfElems arrEl, TypeOfMatrixM typeOfMatrixM)
    {    
        m._al = new double[m._jg.Count];
        m._au = new double[m._jg.Count];

        for (int i = 0; i < arrEl.Length; i++)
            Add(new LocalMatrixNum(arrEl[i], arrPt, typeOfMatrixM, arrEl.mui[i], arrEl.sigmai[i]), ref m, arrEl[i]);
            //Add(new LocalMatrixNum(arrEl[i], arrPt, typeOfMatrixM, arrEl.mui[i], arrEl.sigmai[i]), ref m, arrEl[i]);
/*
        for (int ii = 0; ii < 9; ii++)
        {
            for (int jj = 0; jj < 9; jj++)
            {
                Console.Write($"{m[ii, jj]:E2} ");
            }
            Console.WriteLine();
        }
        Console.WriteLine();
*/
        //ConsiderBoundaryConditions(ref m, arrBd);
    }

    private static void Add(LocalMatrixNum lm, ref GlobalMatrix gm, List<int> elem)
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

    private static void Add(LocalMatrixG3D lm, ref GlobalMatrix gm, List<int> elem)
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

    private static void Add(LocalMatrixM3D lm, ref GlobalMatrix gm, List<int> elem)
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
                            m._al[j] = 0.0D;
                        m._diag[border[i]] = 1;
                        for (int j = 0; j < m._jg.Count; j++)
                            if (m._jg[j] == border[i])
                                m._au[j] = 0.0D;

                        // Гаусс для симметрии:
                        for (int j = 0; j < m._jg.Count; j++)
                            if (m._jg[j] == border[i])
                                m._al[j] = 0.0D;

                        for (int j = m._ig[border[i]]; j < m._ig[border[i] + 1]; j++)
                            m._au[j] = 0.0D;
                    }
                    break;
                }
                case 2: break; // du / dn = 0
                case 3: throw new ArgumentException("Пока нет возможности учитывать КУ 3-го рода");

            }
        }
    }

    public static void FillVector(ref GlobalVector v, ArrayOfPoints arrPt, ArrayOfElems arrEl, double t)
    {
        for (int i = 0; i < arrEl.Length; i++)
        {
            if (arrPt[arrEl[i][0]].R <= 10.0D && 10.0D <= arrPt[arrEl[i][3]].R &&
                arrPt[arrEl[i][0]].Z <= 0.0D && 0.0D <= arrPt[arrEl[i][3]].Z && t <= 1.0D)
                Add(new LocalVector(arrEl[i], arrPt, t), ref v, arrEl[i]);
        }
        //ConsiderBoundaryConditions(ref v, arrPt, arrBd, t);
    }

    private static void Add(LocalVector lv, ref GlobalVector gv, List<int> elems)
    {
        for (int i = 0; i < elems.Count; i++)
        {
            gv[elems[i]] += lv[i];
        }
    }

    public static void ConsiderBoundaryConditions(ref GlobalMatrix gm, ref GlobalVector gv, ArrayOfPoints arrp, ArrayOfBorders arrBd, double t)
    {
        foreach (var border in arrBd)
        {
            switch (border[0])
            {
                // КУ - I-го рода
                case 1:
                for (int i = 2; i < 4; i++)
                {
                    for (int j = gm._ig[border[i]]; j < gm._ig[border[i] + 1]; j++)
                        gm._al[j] = 0.0D;
                    gm._diag[border[i]] = 1.0D;
                    for (int j = 0; j < gm._jg.Count; j++)
                        if (gm._jg[j] == border[i])
                            gm._au[j] = 0.0D;
                }

            /*
            // Гаусс для симметрии:
            if (border[0] == 1)
            {
                for (int i = 2; i < 4; i ++)
                {
                    for (int j = gm._ig[border[i]]; j < gm._ig[border[i] + 1]; j++)
                    {
                        var k = gm._au[j];
                        var f = -1.0D * k * gv[border[i]];
                        gv[gm._jg[j]] += f;
                        gm._au[j] = 0.0D;
                    }

                    for (int j = 0; j < gm._jg.Count; j++)
                        if (gm._jg[j] == border[i])
                        {
                            var k = gm._al[j];
                            var f = -1.0D * k * gv[border[i]];
                            gv[gm._jg[j]] += f;
                            gm._al[j] = 0.0D;
                        }
                }
            }
            */
                switch (border[1])
                {
                    case 1:
                    for (int i = 2; i < 4; i++)
                        gv[border[i]] = Function.U1_1(arrp[border[i]], t);
                    break;
                    
                    case 2:
                    for (int i = 2; i < 4; i++)
                        gv[border[i]] = Function.U1_2(arrp[border[i]], t);
                    break;

                    case 3:
                    for (int i = 2; i < 4; i++)
                        gv[border[i]] = Function.U1_3(arrp[border[i]], t);
                    break;
                    
                    case 4:
                    for (int i = 2; i < 4; i++)
                        gv[border[i]] = Function.U1_4(arrp[border[i]], t);
                    break;
                    
                    default:
                    throw new Exception("No such bord");
                }
                break;
                // КУ - II-го рода
                case 2:
                    for (int i = 2; i < 4; i++)
                        gv[border[i]] += 0.0D;
                    break;
                // КУ - III-го рода
                case 3: throw new ArgumentException("Пока нет возможности учитывать КУ III-го рода");
            }
        }

        foreach (var border in arrBd)
        {
            for (int i = 2; i < 4; i++)
            {
                for (int j = 0; j < gm.Size; j++)
                {
                    if (border[i] == j)
                        continue;
                    else
                    {
                        var k = gm[j, border[i]];
                        var f = -1.0D * k * gv[border[i]];
                        gv[j] += f;
                        gm[j, border[i]] = 0.0D;
                    }
                }
            }
        }
    }
    public static void ConsiderBoundaryConditions(ref GlobalVector v, ArrayOfPoints arrp, ArrayOfBorders arrBd, double t)
    {
        foreach (var border in arrBd)
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


    public static void ConsiderBoundaryConditions(ref GlobalMatrix m, ref GlobalVector v, ArrayOfPoints arrp, ArrayOfBorders arrBrd, double t)
    {
        foreach (var border in arrBrd)
        {
            switch (border[0])
            {
                // КУ - I-го рода
                case 1:
                for (int i = 2; i < 4; i++)
                {
                    for (int j = m._ig[border[i]]; j < m._ig[border[i] + 1]; j++)
                        m._al[j] = 0.0D;
                    m._diag[border[i]] = 1.0D;
                    for (int j = 0; j < m._jg.Count; j++)
                        if (m._jg[j] == border[i])
                            m._au[j] = 0.0D;
                }

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

        /*
        foreach (var border in arrBrd)
        {
            for (int i = 2; i < 4; i++)
            {
                for (int j = 0; j < m.Size; j++)
                {
                    if (border[i] == j)
                        continue;
                    else
                    {
                        var k = m[j, border[i]];
                        var f = -1.0D * k * v[border[i]];
                        v[j] += f;
                        m[j, border[i]] = 0.0D;
                    }
                }
            }
        }
        */
    }
}