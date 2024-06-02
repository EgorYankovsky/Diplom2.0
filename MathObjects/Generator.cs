using DataStructs;
using Functions;
using static Functions.Function;
using Grid;

namespace MathObjects;

public delegate (double, double, double) Egetter(double x, double y, double z, double t);

public static class Generator
{
    public static void FillMatrixG(ref GlobalMatrix m, ArrayOfRibs arrRibs, ArrayOfElems arrEl)
    {
        m._al = new double[m._jg.Count];
        m._au = new double[m._jg.Count];

        for (int i = 0; i < arrEl.Length; i++)
        {
            List<int> currElem = [arrEl[i][0], arrEl[i][3], arrEl[i][8], arrEl[i][11],
                                  arrEl[i][1], arrEl[i][2], arrEl[i][9], arrEl[i][10],
                                  arrEl[i][4], arrEl[i][5], arrEl[i][6], arrEl[i][7]];

            double hx = 0.0;
            double hy = 0.0;
            double hz = 0.0;

            foreach (var rib in currElem)
            {
                if (rib != -1 && hx == 0 && arrRibs[rib].GetTangent().Item1 != 0)
                    hx = arrRibs[rib].Length;
                if (rib != -1 && hy == 0 && arrRibs[rib].GetTangent().Item2 != 0)
                    hy = arrRibs[rib].Length;
                if (rib != -1 && hz == 0 && arrRibs[rib].GetTangent().Item3 != 0)
                    hz = arrRibs[rib].Length;
            }

            var lm = new LocalMatrixG3D(arrEl[i].mu, hx, hy, hz);
            Add(lm, ref m, currElem); 
        }
    }

    public static void FillMatrixM(ref GlobalMatrix m, ArrayOfRibs arrRibs, ArrayOfElems arrEl, Mesh3Dim currentMesh, Mesh3Dim mainMesh)
    {
        m._al = new double[m._jg.Count];
        m._au = new double[m._jg.Count];

        for (int i = 0; i < arrEl.Length; i++)
        {
            List<int> currElem = [arrEl[i][0], arrEl[i][3], arrEl[i][8], arrEl[i][11],
                                  arrEl[i][1], arrEl[i][2], arrEl[i][9], arrEl[i][10],
                                  arrEl[i][4], arrEl[i][5], arrEl[i][6], arrEl[i][7]];
            

            double x0 = double.NaN;
            double x1 = double.NaN;
            double y0 = double.NaN;
            double y1 = double.NaN;
            double z0 = double.NaN;
            double z1 = double.NaN;

            foreach (var rib in currElem)
            {
                if (rib != -1 && double.IsNaN(x0) && double.IsNaN(x1) && arrRibs[rib].GetTangent().Item1 != 0)
                {
                    x0 = arrRibs[rib].a.X;
                    x1 = arrRibs[rib].b.X;
                }
                if (rib != -1 && double.IsNaN(y0) && double.IsNaN(y1) && arrRibs[rib].GetTangent().Item2 != 0)
                {
                    y0 = arrRibs[rib].a.Y;
                    y1 = arrRibs[rib].b.Y;
                }
                if (rib != -1 && double.IsNaN(z0) && double.IsNaN(z1) &&arrRibs[rib].GetTangent().Item3 != 0)
                {
                    z0 = arrRibs[rib].a.Z;
                    z1 = arrRibs[rib].b.Z;
                }
            }

            double hx = x1 - x0;
            double hy = y1 - y0;
            double hz = z1 - z0;

            var sigma = SelectSigma(x0, x1, y0, y1, z0, z1, currentMesh, mainMesh);

            var lm = new LocalMatrixM3D(sigma, hx, hy, hz);
            Add(lm, ref m, currElem);
        }
    }

    public static double SelectSigma(double x0, double x1, double y0, double y1, double z0, double z1, Mesh3Dim currentMesh, Mesh3Dim mainMesh, bool forVector = false)
    {
        double sigma = 0.0D;

        if (currentMesh.nodesX[currentMesh.AnomalyBorders[0]] <= x0 && x1 <= currentMesh.nodesX[currentMesh.AnomalyBorders[1]] &&
            currentMesh.nodesY[currentMesh.AnomalyBorders[2]] <= y0 && y1 <= currentMesh.nodesY[currentMesh.AnomalyBorders[3]] &&
            currentMesh.nodesZ[currentMesh.AnomalyBorders[4]] <= z0 && z1 <= currentMesh.nodesZ[currentMesh.AnomalyBorders[5]] && !forVector)
        {
            foreach (var gg in currentMesh.Elems)
                if (gg.sigma >= 1.0)
                    sigma = gg.sigma;
        }
        else
        {
            for (int i = 0; i < mainMesh.NodesZWithoutFragmentation.Length - 1; i++)
                if (mainMesh.NodesZWithoutFragmentation[i] <= z0 && z1 <= mainMesh.NodesZWithoutFragmentation[i + 1])
                    sigma = mainMesh.Elems[i].sigma;
        }
        return sigma;
    }

    public static void ConsiderBoundaryConditions(ref GlobalMatrix m, ref GlobalVector v, ArrayOfRibs arrRibs, ArrayOfBorders arrBrd, double t)
    {
        List<int> consideredBorder = [];
        foreach (var border in arrBrd)
        {
            switch (border[0])
            {
                // КУ - I-го рода.
                case 1:
                for (int i = 2; i < 6; i++)
                {
                    if (consideredBorder.BinarySearch(border[i]) == 0) continue;

                    consideredBorder.Add(border[i]);
                    consideredBorder.Sort();
                    for (int j = m._ig[border[i]]; j < m._ig[border[i] + 1]; j++)
                        m._al[j] = 0.0D;
                    m._diag[border[i]] = 1.0D;
                    for (int j = border[i]; j < m._ig.Length - 1; j++)
                    {
                        var row = m._jg[m._ig[j] .. m._ig[j + 1]];
                        if (row.BinarySearch(border[i]) < 0) continue;
                        for (int jj = 0; jj < row.Count; jj++)
                            if (row[jj] == border[i])
                                m._au[row[jj]] = 0.0D;
                    }
                    //for (int j = 0; j < m._jg.Count; j++)
                    //    if (m._jg[j] == border[i])
                    //        m._au[j] = 0.0D;
                    v[border[i]] = 0.0D;
                }
                
                break;
                // КУ - II-го рода
                case 2:
                    for (int i = 2; i < 6; i++)
                        v[border[i]] += 0.0D;
                    break;
                // КУ - III-го рода
                case 3: throw new ArgumentException("Пока нет возможности учитывать КУ III-го рода");
            }
        }
    }

    public static void FillVector3D(ref GlobalVector v, Egetter egetter, ArrayOfRibs arrRibs, ArrayOfElems arrEl, Mesh3Dim currentMesh, Mesh3Dim mainMesh, double t)
    {
        for (int i = 0; i < arrEl.Length; i++)
        {
            List<int> currElem = [arrEl[i][0], arrEl[i][3], arrEl[i][8], arrEl[i][11],
                                  arrEl[i][1], arrEl[i][2], arrEl[i][9], arrEl[i][10],
                                  arrEl[i][4], arrEl[i][5], arrEl[i][6], arrEl[i][7]];
  
            double x0 = double.NaN;
            double x1 = double.NaN;
            double y0 = double.NaN;
            double y1 = double.NaN;
            double z0 = double.NaN;
            double z1 = double.NaN;

            foreach (var rib in currElem)
            {
                if (rib != -1 && double.IsNaN(x0) && double.IsNaN(x1) && arrRibs[rib].GetTangent().Item1 != 0)
                {
                    x0 = arrRibs[rib].a.X;
                    x1 = arrRibs[rib].b.X;
                }
                if (rib != -1 && double.IsNaN(y0) && double.IsNaN(y1) && arrRibs[rib].GetTangent().Item2 != 0)
                {
                    y0 = arrRibs[rib].a.Y;
                    y1 = arrRibs[rib].b.Y;
                }
                if (rib != -1 && double.IsNaN(z0) && double.IsNaN(z1) && arrRibs[rib].GetTangent().Item3 != 0)
                {
                    z0 = arrRibs[rib].a.Z;
                    z1 = arrRibs[rib].b.Z;
                }
            }

            var originalSigma = SelectSigma(x0, x1, y0, y1, z0, z1, currentMesh, mainMesh, true);
            var currentSigma = arrEl[i].sigma;

            var lv = new LocalVector3D(egetter, x0, x1, y0, y1, z0, z1, t);

            Add(lv, ref v, currentSigma - originalSigma, currElem);
        }
    }

    private static void Add(LocalVector3D lv, ref GlobalVector v, double sigma, List<int> elem)
    {
        for (int i = 0; i < 12; i++)
            if(elem[i] != -1) 
                v[elem[i]] += sigma * lv[i];
    }

    public static void BuildPortait(ref GlobalMatrix m, int arrPtLen, ArrayOfElems arrEl, bool isVec = false)
    {
        List<List<int>> arr = [];

        // ! Дерьмодристный момент.
        for(int i = 0; i < arrPtLen; i++)
            arr.Add([]);

        foreach (Elem _elem in arrEl)
            foreach (var point in _elem.Arr)
                foreach (var pnt in _elem.Arr)
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

    public static void FillMatrix(ref GlobalMatrix m, ArrayOfPoints2D arrPt, ArrayOfElems arrEl, TypeOfMatrixM typeOfMatrixM)
    {    
        m._al = new double[m._jg.Count];
        m._au = new double[m._jg.Count];

        for (int i = 0; i < arrEl.Length; i++)
            Add(new LocalMatrix(arrEl[i].Arr, arrPt, typeOfMatrixM, arrEl[i].mu, arrEl[i].sigma), ref m, arrEl[i].Arr);
            //Add(new LocalMatrixNum(arrEl[i], arrPt, typeOfMatrixM, arrEl.mui[i], arrEl.sigmai[i]), ref m, arrEl[i]);
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

    private static void Add(LocalMatrixG3D lm, ref GlobalMatrix gm, List<int> elem)
    {
        if (gm._diag is null) throw new Exception("_diag isn't initialized");
        if (gm._ig is null) throw new Exception("_ig isn't initialized");
        if (gm._au is null) throw new Exception("_au isn't initialized");
        if (gm._al is null) throw new Exception("_au isn't initialized");
        
        
        int ii = 0;
        foreach (var i in elem)
        {
            if (i != -1)
            {
                int jj = 0;
                foreach (var j in elem)
                {
                    if (j != -1)
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
                    }
                    jj++;
                }
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
            if (i != -1)
            {
                int jj = 0;
                foreach (var j in elem)
                {
                    if (j != -1)
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
                    }
                    jj++;
                }
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

    public static void FillVector(ref GlobalVector v, ArrayOfPoints2D arrPt, ArrayOfElems arrEl, double t)
    {
        for (int i = 0; i < arrPt.GetLength(); i++)
            if (arrPt[i].R == 1000.0D && arrPt[i].Z == 0.0D && t <= 0.0D)
                v[i] = F(arrPt[i].R, arrPt[i].Z, t);

            
        //for (int i = 0; i < arrEl.Length; i++)
        //{
        //    if (arrPt[arrEl[i][0]].R <= 500.0D && 500.0D <= arrPt[arrEl[i][3]].R &&
        //        arrPt[arrEl[i][0]].Z <= 0.0D && 0.0D <= arrPt[arrEl[i][3]].Z && t <= 1.0D)
        //        Add(new LocalVector(arrEl[i].Arr, arrPt, t), ref v, arrEl[i].Arr);
        //}
    }

    private static void Add(LocalVector lv, ref GlobalVector gv, List<int> elems)
    {
        for (int i = 0; i < elems.Count; i++)
        {
            gv[elems[i]] += lv[i];
        }
    }

    public static void ConsiderBoundaryConditions(ref GlobalMatrix gm, ref GlobalVector gv, ArrayOfPoints2D arrp, ArrayOfBorders arrBd, double t)
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
                    //for (int row = 0; row < gm.Size; row++)
                    //{
                    //    var stB = gv[border[i]];
                    //    for (int jj = gm._ig[row]; jj < gm._ig[row + 1]; jj++)
                    //    {
                    //        if (gm._jg[jj] == border[i])
                    //        {
                    //            stB *= -1.0D * gm._al[jj];
                    //            gv[row] += stB;
                    //            gm._al[jj] = 0.0D;
                    //        }    
                    //    }
                    //}
                    //for (int j = border[i]; j < gm._ig.Length - 1; j++)
                    //{
                    //    var row = gm._jg[gm._ig[j] .. gm._ig[j + 1]];
                    //    if (row.BinarySearch(border[i]) < 0) continue;
                    //    for (int jj = 0; jj < row.Count; jj++)
                    //        if (row[jj] == border[i])
                    //            gm._au[gm._ig[j] + jj] = 0.0D;
                    //}
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
                        gv[border[i]] = U1_1(arrp[border[i]], t);
                    break;
                    
                    case 2:
                    for (int i = 2; i < 4; i++)
                        gv[border[i]] = U1_2(arrp[border[i]], t);
                    break;

                    case 3:
                    for (int i = 2; i < 4; i++)
                        gv[border[i]] = U1_3(arrp[border[i]], t);
                    break;
                    
                    case 4:
                    for (int i = 2; i < 4; i++)
                        gv[border[i]] = U1_4(arrp[border[i]], t);
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

        //foreach (var border in arrBd)
        //{
        //    for (int i = 2; i < 4; i++)
        //    {
        //        for (int j = 0; j < gm.Size; j++)
        //        {
        //            if (border[i] == j)
        //                continue;
        //            else
        //            {
        //                var k = gm[j, border[i]];
        //                var f = -1.0D * k * gv[border[i]];
        //                gv[j] += f;
        //                gm[j, border[i]] = 0.0D;
        //            }
        //        }
        //    }
        //}
    }
    public static void ConsiderBoundaryConditions(ref GlobalVector v, ArrayOfPoints2D arrp, ArrayOfBorders arrBd, double t)
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
                            v[border[i]] = U1_1(arrp[border[i]], t);
                        break;
                        
                        case 2:
                        for (int i = 2; i < 4; i++)
                            v[border[i]] = U1_2(arrp[border[i]], t);
                        break;

                        case 3:
                        for (int i = 2; i < 4; i++)
                            v[border[i]] = U1_3(arrp[border[i]], t);
                        break;
                        
                        case 4:
                        for (int i = 2; i < 4; i++)
                            v[border[i]] = U1_4(arrp[border[i]], t);
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