using DataStructs;

namespace Grid;

public static class MeshGenerator
{
    private static readonly string _points2DPath = Path.GetFullPath("../../../../Data/Subtotals/2_dim/Points.poly");

    private static readonly string _points3DPath = Path.GetFullPath("../../../../Data/Subtotals/3_dim/Points.poly");

    private static readonly string _elems2DPath = Path.GetFullPath("../../../../Data/Subtotals/2_dim/Elems.poly");

    private static readonly string _borders2DPath = Path.GetFullPath("../../../../Data/Subtotals/2_dim/Borders.poly");

    public static TimeMesh GenerateTimeMesh(double t0, double t1, int tn, double tk)
    {
        double[] arr = new double[tn];
        double h = t1 - t0;
        double denominator = 0.0;

        for (int j = 0; j < tn; j++)
            denominator += Math.Pow(tk, j);

        double x0 = h / denominator;
        arr[0] = t0;
        for(int j = 0; j < tn - 1; j++)
            arr[j + 1] = arr[j] + x0 * Math.Pow(tk, j);
        arr[^1] = t1;
        return new TimeMesh(arr);
    }

    public static void GenerateMesh(ref Mesh2Dim mesh)
    {
        int currentPosition = 0;
        string[] kek = mesh.infoAboutR.Split();
        for (int i = 0; i < mesh.NodesRWithoutFragmentation.Length - 1; i++)
        {
            double h = mesh.nodesR[1 + currentPosition] - mesh.nodesR[currentPosition];
            double denominator = 0.0;

            int negr = int.Parse(kek[2 * i]);
            for (int j = 0; j < negr; j++)
                denominator += Math.Pow(double.Parse(kek[2 * i + 1]), j);

            double x0 = h / denominator;

            for(int j = 0; j < negr - 1; j++)
            {
                mesh.nodesR.Insert(currentPosition + 1, mesh.nodesR[currentPosition] + x0 * Math.Pow(double.Parse(kek[2 * i + 1]), j));
                currentPosition++;
            }
            currentPosition++;
            mesh.NodesR_Refs.Add(currentPosition);
        }
        
        currentPosition = 0;
        kek = mesh.infoAboutZ.Split();
        for (int i = 0; i < mesh.NodesZWithoutFragmentation.Length - 1; i++)
        {
            double h = mesh.nodesZ[1 + currentPosition] - mesh.nodesZ[currentPosition];
            double denominator = 0.0;

            int negr = int.Parse(kek[2 * i]);
            for (int j = 0; j < negr; j++)
                denominator += Math.Pow(double.Parse(kek[2 * i + 1]), j);
            double x0 = h / denominator;

            for(int j = 0; j < negr - 1; j++)
            {
                mesh.nodesZ.Insert(currentPosition + 1, mesh.nodesZ[currentPosition] + x0 * Math.Pow(double.Parse(kek[2 * i + 1]), j));
                currentPosition++;
            }
            currentPosition++;
            mesh.NodesZRefs.Add(currentPosition);
        }
    }

    public static void GenerateMesh(ref Mesh3Dim mesh)
    {
        int currentPosition = 0;
        string[] kek = mesh.infoAboutX.Split();
        for (int i = 0; i < mesh.NodesXWithoutFragmentation.Length - 1; i++)
        {
            double h = mesh.nodesX[1 + currentPosition] - mesh.nodesX[currentPosition];
            double denominator = 0.0;

            int negr = int.Parse(kek[2 * i]);
            for (int j = 0; j < negr; j++)
                denominator += Math.Pow(double.Parse(kek[2 * i + 1]), j);

            double x0 = h / denominator;

            for(int j = 0; j < negr - 1; j++)
            {
                mesh.nodesX.Insert(currentPosition + 1, mesh.nodesX[currentPosition] + x0 * Math.Pow(double.Parse(kek[2 * i + 1]), j));
                currentPosition++;
            }
            currentPosition++;
            mesh.NodesXRefs.Add(currentPosition);
        }

        currentPosition = 0;
        kek = mesh.infoAboutY.Split();
        for (int i = 0; i < mesh.NodesYWithoutFragmentation.Length - 1; i++)
        {
            double h = mesh.nodesY[1 + currentPosition] - mesh.nodesY[currentPosition];
            double denominator = 0.0;

            int negr = int.Parse(kek[2 * i]);
            for (int j = 0; j < negr; j++)
                denominator += Math.Pow(double.Parse(kek[2 * i + 1]), j);

            double x0 = h / denominator;

            for(int j = 0; j < negr - 1; j++)
            {
                mesh.nodesY.Insert(currentPosition + 1, mesh.nodesY[currentPosition] + x0 * Math.Pow(double.Parse(kek[2 * i + 1]), j));
                currentPosition++;
            }
            currentPosition++;
            mesh.NodesYRefs.Add(currentPosition);
        }
        
        currentPosition = 0;
        kek = mesh.infoAboutZ.Split();
        for (int i = 0; i < mesh.NodesZWithoutFragmentation.Length - 1; i++)
        {
            double h = mesh.nodesZ[1 + currentPosition] - mesh.nodesZ[currentPosition];
            double denominator = 0.0;

            int negr = int.Parse(kek[2 * i]);
            for (int j = 0; j < negr; j++)
                denominator += Math.Pow(double.Parse(kek[2 * i + 1]), j);
            double x0 = h / denominator;

            for(int j = 0; j < negr - 1; j++)
            {
                mesh.nodesZ.Insert(currentPosition + 1, mesh.nodesZ[currentPosition] + x0 * Math.Pow(double.Parse(kek[2 * i + 1]), j));
                currentPosition++;
            }
            currentPosition++;
            mesh.NodesZRefs.Add(currentPosition);
        }
    }

    public static void ConstructMesh(ref Mesh2Dim mesh)
    {
        GenerateMesh(ref mesh);
        RemakeBorders(ref mesh);
        OutputPoints(ref mesh);
        OutputListOfElems(ref mesh);
        GenerateListOfBorders(ref mesh);
    }

    public static void ConstructMesh(ref Mesh3Dim mesh)
    {
        GenerateMesh(ref mesh);
        RemakeBorders(ref mesh);
        OutputPoints(ref mesh);
    }

    public static void RemakeBorders(ref Mesh3Dim mesh)
    {
        for (int i = 0; i < mesh.borders.Count; i++)
        {
            var newBorder = new Border3D(mesh.borders[i].BorderType, mesh.borders[i].BorderFormula,
                                         mesh.NodesXRefs[mesh.borders[i].X0], mesh.NodesXRefs[mesh.borders[i].X1],
                                         mesh.NodesYRefs[mesh.borders[i].Y0], mesh.NodesYRefs[mesh.borders[i].Y1],
                                         mesh.NodesZRefs[mesh.borders[i].Z0], mesh.NodesZRefs[mesh.borders[i].Z1]);
            mesh.borders[i] = newBorder;
        }
    }

    public static void RemakeBorders(ref Mesh2Dim mesh)
    {
        for (int i = 0; i < mesh.borders.Count; i++)
        {
            var newBorder = new Border2D(mesh.borders[i].BorderType, mesh.borders[i].BorderFormula,
                                         mesh.NodesR_Refs[mesh.borders[i].R0], mesh.NodesR_Refs[mesh.borders[i].R1],
                                         mesh.NodesZRefs[mesh.borders[i].Z0], mesh.NodesZRefs[mesh.borders[i].Z1]);
            mesh.borders[i] = newBorder;
        }
    }

    public static void OutputPoints(ref Mesh2Dim mesh)
    {
        using var sw = new StreamWriter(_points2DPath);
        sw.WriteLine(mesh.NodesAmountTotal);
        int i = 0;
        foreach (var Z in mesh.nodesZ)
            foreach (var R in mesh.nodesR)
            {
                sw.WriteLine($"{SetPointType(mesh, new Point(R, Z))}");
                i++;
            }
    }

    public static void OutputPoints(ref Mesh3Dim mesh)
    {
        using var sw = new StreamWriter(_points3DPath);
        sw.WriteLine(mesh.NodesAmountTotal);
        int i = 0;
        foreach (var Z in mesh.nodesZ)
            foreach (var Y in mesh.nodesY)
                foreach (var X in mesh.nodesX)
                {
                    sw.WriteLine($"{SetPointType(mesh, new Point(X, Y, Z))}");
                    i++;
                }
    }

    private static Point SetPointType(Mesh2Dim mesh, Point pnt)
    {
        for (int i = 0; i < mesh.Elems.Count; i++)
        {   
            if (mesh.NodesZWithoutFragmentation[mesh.Elems[i].Arr[4]] <= pnt.Z && pnt.Z <= mesh.NodesZWithoutFragmentation[mesh.Elems[i].Arr[5]])
            {
                int minValue = 4;
                foreach (var arr in mesh.borders)
                {
                    if (mesh.nodesR[arr.R0] <= pnt.R && pnt.R <= mesh.nodesR[arr.R1] &&
                        mesh.nodesZ[arr.Z0] <= pnt.Z && pnt.Z <= mesh.nodesZ[arr.Z1] &&
                        arr.BorderType < minValue)
                        {
                            minValue = arr.BorderType;
                            break;
                        }
                }
                switch (minValue)
                {
                    case 1: {pnt.Type = Location.BoundaryI; break; }
                    case 2: {pnt.Type = Location.BoundaryII; break; }
                    case 3: {pnt.Type = Location.BoundaryIII; break; }
                    default: {pnt.Type = Location.Inside; pnt.SubElemNum = i; break;}
                }
            }
            else if (pnt.Type == Location.NotStated)
                pnt.Type = Location.OutSide;
        }
        return pnt;
    }

    private static Point SetPointType(Mesh3Dim mesh, Point pnt)
    {
        return new Point(0.0, 0.0, 0.0);
    }

    public static void OutputListOfElems(ref Mesh2Dim mesh)
    {
        using var sw = new StreamWriter(_elems2DPath);
        int nr = mesh.NodesAmountR;
        int nz = mesh.NodesAmountZ;
        sw.WriteLine((nr - 1) * (nz - 1));
        for (int k = 0; k < nz - 1; k++)
            for (int i = 0; i < nr - 1; i++)
                sw.WriteLine($"{k * nr + i} {k * nr + i + 1} " +
                             $"{(k + 1) * nr + i} {(k + 1) * nr + i + 1} " +
                             $"{SelectMuAndSigma(mesh, k, k + 1)}");
    }

    private static string SelectMuAndSigma(Mesh2Dim mesh, int a, int d)
    {
        int index = -1;
        for (int i = 0; i < mesh.NodesZRefs.Count - 1; i++)
            if (mesh.NodesZRefs[i] <= a && d <= mesh.NodesZRefs[i + 1])
                index = i;
        if (index == -1) throw new IndexOutOfRangeException("Failure during selecting mu and sigma");
        return $"{mesh.Elems[index].mu} {mesh.Elems[index].sigma}";
    }

    /*
    public static ArrayOfPoints GenerateListOfPoints(Mesh mesh)
    {
        if (mesh.nodesX is null) throw new ArgumentNullException("nodesX is null");
        if (mesh.nodesY is null) throw new ArgumentNullException("nodesY is null");
        if (mesh.nodesZ is null) throw new ArgumentNullException("nodesZ is null");
        
        var arr = new ArrayOfPoints(mesh.NodesAmountTotal);
        foreach (var Z in mesh.nodesZ)
            foreach (var Y in mesh.nodesY)
                foreach(var X in mesh.nodesX)
                    arr.Append(SetPointType3D(new Point(X, Y, Z), mesh));
        return arr;
    }
    */

    // ! Топорный метод.
    /*
    private static Point SetPointType3D(Point pnt, Mesh mesh)
    {
        if (mesh.nodesX is null) throw new ArgumentNullException("nodesX is null");
        if (mesh.nodesY is null) throw new ArgumentNullException("nodesY is null");
        if (mesh.nodesZ is null) throw new ArgumentNullException("nodesZ is null");

        double xMin = mesh.nodesX[0];
        double xMax = mesh.nodesX[^1];
        double yMin = mesh.nodesY[0];
        double yMax = mesh.nodesY[^1];
        double zMin = mesh.nodesZ[0];
        double zMax = mesh.nodesZ[^1];

        if (pnt.Z == zMax)
            pnt.Type = Location.BoundaryI;
        else if (pnt.Z == zMin || pnt.X == xMin ||pnt.X == xMax || pnt.Y == yMin || pnt.Y == yMax)
            pnt.Type = Location.BoundaryI;
        else
            pnt.Type = Location.Inside;
    
        return pnt;
    }
    */

    public static ArrayOfElems GenerateListOfElems(Mesh mesh)
    {
        var arr = new ArrayOfElems(mesh.ElemsAmount);
        
        int rx = mesh.NodesAmountX - 1;
        int ry = mesh.NodesAmountY - 1;
        
        int nx = mesh.NodesAmountX;
        int ny = mesh.NodesAmountY;

        int rxy = rx * ny + ry * nx;
        int nxy = nx * ny;
        int nz = mesh.NodesAmountZ;

        for (int k = 0; k < nz - 1; k++)
            for (int j = 0; j < ny - 1; j++)
                for (int i = 0; i < nx - 1; i++)
                {
                    int curr = i + j * (nx + rx) + k * (rxy + nxy);
                    arr.Add([               curr,               curr + rx,             curr + rx + 1,               curr + rx + nx,
                             curr + rxy - j * rx, curr + rxy + 1 - j * rx,  curr + rxy + nx - j * rx, curr + rxy + nx + 1 - j * rx,
                                curr + rxy + nxy,   curr + rxy + nxy + rx, curr + rxy + nxy + rx + 1,   curr + rxy + nxy + rx + nx]);
                    arr.mui.Add(1.0D);
                    arr.sigmai.Add(1.0D);
                }
        return arr;
    }
    */

    public static void SelectRibs(ref ArrayOfRibs arrRibs, ref ArrayOfElems arrEl)
    {
        int ii = 0;
        while (ii < arrRibs.Count)
        {
            Console.WriteLine($"{ii * 100.0D / arrRibs.Count:E5}\r");
            if (arrRibs[ii].typeOfRib == TypeOfRib.BoundaryI)
            {
                foreach (var elem in arrEl)
                    foreach (var item in elem)
                    {
                        if (item > ii) break;
                        if (item == ii)
                        {
                            elem.Remove(item);
                            break;
                        }
                    }
                arrRibs.Remove(ii);
                for (int i = 0; i < arrEl.Length; i++)
                    for (int j = 0; j < arrEl[i].Count; j++) 
                        if (arrEl[i][j] > ii)
                            arrEl[i][j] -= 1;
            }
            else
                ii++;
        }
    }

    /*
    public static ArrayOfBorders GenerateListOfBorders(Mesh mesh)
    {
        var arr = new ArrayOfBorders(null);

        int nx = mesh.NodesAmountX;
        int ny = mesh.NodesAmountY;
        int nz = mesh.NodesAmountZ;

        int nxny = nx * ny;
        int rxy = (nx - 1) * ny + (ny - 1) * nx;
        
        // XY0
        for (int i = 0; i < ny - 1; i++)
            for (int j = 0; j < nx - 1; j++)
                arr.Add([1, 1, i * (2 * nx - 1) + j, 
                               i * (2 * nx - 1) + j + nx - 1,
                               i * (2 * nx - 1) + j + nx,
                               i * (2 * nx - 1) + j + nx + nx - 1]);
        
        
        // X0Z
        for (int i = 0; i < nz - 1; i++)
            for (int j = 0; j < nx - 1; j++)
                arr.Add([1, 2, i * (rxy + nxny) + j, 
                               i * (rxy + nxny) + j + rxy,
                               i * (rxy + nxny) + j + rxy + 1,
                               i * (rxy + nxny) + j + rxy + nxny]);
        
        
        // 0YZ
        for (int i = 0; i < nz - 1; i++)
            for (int j = 0; j < ny - 1; j++)
                arr.Add([1, 3, nx - 1 + i * (rxy + nxny) + j * (2 * nx - 1), 
                               rxy + i * (rxy + nxny) + j * nx,
                               rxy + nx + i * (rxy + nxny) + j * nx,
                               rxy + nxny + nx - 1 + i * (rxy + nxny) + j * (2 * nx - 1)]);
        
        // XY1
        for (int i = 0; i < ny - 1; i++)
            for (int j = 0; j < nx - 1; j++)
                arr.Add([1, 4, (nz - 1) * (rxy + nxny) + j + i * (2 * nx - 1),
                               (nz - 1) * (rxy + nxny) + nx - 1 + j + i * (2 * nx - 1),
                               (nz - 1) * (rxy + nxny) + nx + j + i * (2 * nx - 1),
                               (nz - 1) * (rxy + nxny) + nx + nx - 1 + j + i * (2 * nx - 1)]);
        
        // X1Z
        for (int i = 0; i < nz - 1; i++)
            for (int j = 0; j < nx - 1; j++)
                arr.Add([1, 5, (ny - 1) * nx + (ny - 1) * (nx - 1) + j + i * (rxy + nxny),
                               (ny - 1) * nx + (ny - 1) * (nx - 1) + nxny - 1 + j + i * (rxy + nxny),
                               (ny - 1) * nx + (ny - 1) * (nx - 1) + nxny + j + i * (rxy + nxny),
                               (ny - 1) * nx + (ny - 1) * (nx - 1) + nxny + j + rxy + i * (rxy + nxny)]);

        // 1YZ
        for (int i = 0; i < nz - 1; i++)
            for (int j = 0; j < ny - 1; j++)
                arr.Add([1, 6, nx - 1 + i * (rxy + nxny) + j * (2 * nx - 1) + nx - 1, 
                               rxy + i * (rxy + nxny) + j * nx + nx - 1,
                               rxy + nx + i * (rxy + nxny) + j * nx + nx - 1,
                               rxy + nxny + nx - 1 + i * (rxy + nxny) + j * (2 * nx - 1) + nx - 1]);
        

        return arr;
    }
    */

    /*
    public static ArrayOfRibs GenerateListOfRibs(Mesh mesh, ArrayOfPoints arrPt)
    {
        ArgumentNullException.ThrowIfNull(arrPt);

        int nx = mesh.NodesAmountX;
        int ny = mesh.NodesAmountY;
        int nz = mesh.NodesAmountZ;
        int nxny = nx * ny;

        ArrayOfRibs arr = new(3 * nx * ny * nz - nx * ny - nx * nz - ny * nz);

        // Генерируем список всех ребер.
        for (int k = 0; k < nz; k++)
        {
            for (int j = 0; j < ny; j++)
            {
                for (int i = 0; i < nx - 1; i++) 
                    arr.Add(new Rib(arrPt[k * nxny + nx * j + i], arrPt[k * nxny + nx * j + i + 1]));
                if (j != ny - 1)
                    for (int i = 0; i < nx; i++)
                        arr.Add(new Rib(arrPt[k * nxny + nx * j + i], arrPt[k * nxny + nx * (j + 1) + i]));
            }
            if (k != nz - 1)
                for (int j = 0; j < ny; j++)
                    for (int i = 0; i < nx; i++)
                        arr.Add(new Rib(arrPt[k * nxny + nx * j + i], arrPt[(k + 1) * nxny + nx * j + i]));
        }
        return arr;
    }
    */
    public static void GenerateListOfBorders(ref Mesh2Dim mesh)
    {
        using var sw = new StreamWriter(_borders2DPath);

        sw.WriteLine(2 * (mesh.NodesAmountR - 1 + mesh.NodesAmountZ - 1));
        
        foreach (var border in mesh.borders)
        {
            if (border.R0 == border.R1) // border || oR
            {
                int iter = border.R0 == 0 ? 0 : mesh.NodesAmountR - 1;
                for (int i = 0; i < mesh.NodesAmountZ - 1; i++)
                {
                    sw.WriteLine($"{border.BorderType} {border.BorderFormula} {iter} {iter + mesh.NodesAmountR}");
                    iter += mesh.NodesAmountR;
                }
            }
            else if (border.Z0 == border.Z1) // border || oZ
            {
                int iter = border.Z0 == 0 ? 0 : mesh.NodesAmountR * (mesh.NodesAmountZ - 1);
                for (int i = 0; i < mesh.NodesAmountR - 1; i++)
                {
                    sw.WriteLine($"{border.BorderType} {border.BorderFormula} {iter} {iter + 1}");
                    iter++;
                }
            }
        }
    }
}