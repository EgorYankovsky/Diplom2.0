using DataStructs;

namespace Grid;

public static class MeshGenerator
{
    private static readonly string _pointsPath = Path.GetFullPath("../../../../Data/Subtotals/Points.dat");

    private static readonly string _elemsPath = Path.GetFullPath("../../../../Data/Subtotals/Elems.dat");

    private static readonly string _bordersPath = Path.GetFullPath("../../../../Data/Subtotals/Borders.dat");


    public static void GenerateMesh(ref Mesh2Dim mesh)
    {
        int currentPosition = 0;
        double[] kek = mesh.infoAboutR.Split().Select(double.Parse).ToArray();
        for (int i = 0; i < mesh.NodesRWithoutFragmentation.Length - 1; i++)
        {
            double h = mesh.nodesR[1 + currentPosition] - mesh.nodesR[currentPosition];
            double denominator = 0.0;

            int negr = Convert.ToInt32(kek[2 * i]);
            for (int j = 0; j < negr; j++)
                denominator += Math.Pow(kek[2 * i + 1], j);

            double x0 = h / denominator;

            for(int j = 0; j < negr - 1; j++)
            {
                mesh.nodesR.Insert(currentPosition + 1, mesh.nodesR[currentPosition] + x0 * Math.Pow(kek[2 * i + 1], j));
                currentPosition++;
            }
            currentPosition++;
            mesh.nodesR_Refs[i + 1] = currentPosition;
        }
        
        currentPosition = 0;
        kek = mesh.infoAboutZ.Split().Select(double.Parse).ToArray();
        for (int i = 0; i < mesh.NodesZWithoutFragmentation.Length - 1; i++)
        {
            double h = mesh.nodesZ[1 + currentPosition] - mesh.nodesZ[currentPosition];
            double denominator = 0.0;

            int negr = Convert.ToInt32(kek[2 * i]);
            for (int j = 0; j < negr; j++)
                denominator += Math.Pow(kek[2 * i + 1], j);
            double x0 = h / denominator;

            for(int j = 0; j < negr - 1; j++)
            {
                mesh.nodesZ.Insert(currentPosition + 1, mesh.nodesZ[currentPosition] + x0 * Math.Pow(kek[2 * i + 1], j));
                currentPosition++;
            }
            currentPosition++;
            mesh.nodesZRefs[i + 1] = currentPosition;
        }
        RemakeBorders(ref mesh);
    }

    private static void RemakeBorders(ref Mesh2Dim mesh)
    {
        foreach (var border in mesh.borders)
        {
            for (int i = 2; i < border.Count; i++)
            {
                switch(i)
                {
                    case 2:
                    case 3:
                        border[i] = mesh.nodesR_Refs[border[i]];
                        break;
                    case 4:
                    case 5:
                        border[i] = mesh.nodesZRefs[border[i]];
                        break;
                }
            }
        }
    }

    public static void OutputPoints(ref Mesh2Dim mesh)
    {
        using var sw = new StreamWriter(_pointsPath);
        sw.WriteLine(mesh.NodesAmountTotal);
        int i = 0;
        foreach (var Z in mesh.nodesZ)
            foreach (var R in mesh.nodesR)
            {
                sw.WriteLine($"{SetPointType(mesh, new Point(R, Z))}");
                i++;
            }
    }

    private static Point SetPointType(Mesh2Dim mesh, Point pnt)
    {
        for (int i = 0; i < mesh.Elems.Count; i++)
        {   
            if (mesh.NodesRWithoutFragmentation[mesh.Elems[i][1]] <= pnt.R && pnt.R <= mesh.NodesRWithoutFragmentation[mesh.Elems[i][2]] &&
                mesh.NodesZWithoutFragmentation[mesh.Elems[i][3]] <= pnt.Z && pnt.Z <= mesh.NodesZWithoutFragmentation[mesh.Elems[i][4]])
            {
                int minValue = 4;
                foreach (var arr in mesh.borders)
                {
                    if (mesh.nodesR[arr[2]] <= pnt.R && pnt.R <= mesh.nodesR[arr[3]] &&
                        mesh.nodesZ[arr[4]] <= pnt.Z && pnt.Z <= mesh.nodesZ[arr[5]] &&
                        arr[0] < minValue)
                        {
                            minValue = arr[0];
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

    public static void GenerateListOfElems(ref Mesh2Dim mesh, ArrayOfPoints arrPt)
    {
        using var sw = new StreamWriter(_elemsPath);
        sw.WriteLine((mesh.NodesAmountR - 1) * (mesh.NodesAmountZ - 1));
        for (int k = 0; k < mesh.NodesAmountZ - 1; k++)
            for (int i = 0; i < mesh.NodesAmountR - 1; i++)
                sw.WriteLine($"{k * mesh.NodesAmountR + i} {k * mesh.NodesAmountR + i + 1} " +
                             $"{(k + 1) * mesh.NodesAmountR + i} {(k + 1) * mesh.NodesAmountR + i + 1} " +
                             $"{SelectMuAndSigma(mesh, arrPt[k * mesh.NodesAmountR + i], 
                                                       arrPt[k * mesh.NodesAmountR + i + 1],
                                                       arrPt[(k + 1) * mesh.NodesAmountR + i],
                                                       arrPt[(k + 1) * mesh.NodesAmountR + i + 1])}");
    }

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

    // ! Топорный метод.
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

    // ! Костыль вселенских масштабов.
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
        int ramount = 3 * nx * ny * nz - nx * ny - nx * nz - ny * nz;

        int j = 0;
        int k = 0;
        for (int i = 0; i < nx * ny * nz; i++)
        {
            //Console.WriteLine($"{i} | {nx * ny * nz}");
            if (i + ry * j + rxy * k + rxy + nxy + rx + nx > ramount) 
                break;

            if (i % (nx * ry - 1 + k * nx * ny) == 0 && i != 0)
            {
                i += nx;
                k++;
                j = 0;
            }
            else if (i % nx == nx - 1)
                j++;
            else
            {
                int fst = i + ry * j + rxy * k; 
                arr.Add([i + ry * j + rxy * k, i + rx + ry * j + rxy * k, i + rx + 1 + ry * j + rxy * k, i + rx + nx + ry * j + rxy * k,
                        i + rxy * (k + 1), i + rxy * (k + 1) + 1, i + rxy * (k + 1) + rx + 1, i + rxy * (k + 1) + rx + 1 + 1,
                        fst + rxy + nxy, fst + rxy + nxy + rx, fst + rxy + nxy + rx + 1, fst + rxy + nxy + rx + nx]);
                //arr.mui.Add(mesh.mu0[0]);
                //arr.sigmai.Add(mesh.nodesZ[i / (nx * ny)] <= 0.0D ? mesh.sigma[0] : mesh.sigma[^1]);
                arr.mui.Add(1.0D);
                arr.sigmai.Add(1.0D);
            }
        }
        return arr;
    }

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
            for (int j = 0; j < nx - 1; j++)
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
            for (int j = 0; j < nx - 1; j++)
                arr.Add([1, 6, nx - 1 + i * (rxy + nxny) + j * (2 * nx - 1) + nx - 1, 
                               rxy + i * (rxy + nxny) + j * nx + nx - 1,
                               rxy + nx + i * (rxy + nxny) + j * nx + nx - 1,
                               rxy + nxny + nx - 1 + i * (rxy + nxny) + j * (2 * nx - 1) + nx - 1]);
        

        return arr;
    }

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

    // ! Метод написан очень топорно, только для конкретного примера.
    private static string SelectMuAndSigma(Mesh2Dim mesh, Point a, Point b, Point c, Point d)
    {
        int index;
        if (mesh.NodesZWithoutFragmentation[0] <= a.Z && d.Z <= mesh.NodesZWithoutFragmentation[1])
            index = 0;
        else if (mesh.NodesZWithoutFragmentation[1] <= a.Z && d.Z <= mesh.NodesZWithoutFragmentation[2])
            index = 1;
        else if (mesh.NodesZWithoutFragmentation[2] <= a.Z && d.Z <= mesh.NodesZWithoutFragmentation[3])
            index = 1;
        else if (mesh.NodesZWithoutFragmentation[3] <= a.Z && d.Z <= mesh.NodesZWithoutFragmentation[4])
            index = 1;
        else if (mesh.NodesZWithoutFragmentation[4] <= a.Z && d.Z <= mesh.NodesZWithoutFragmentation[5])
            index = 1;
        else
            throw new Exception("Out of boundary!");
        return $"{mesh.mu0[index]} {mesh.sigma[index]}";
    }

    public static void GenerateListOfBorders(ref Mesh2Dim mesh)
    {
        using var sw = new StreamWriter(_bordersPath);

        sw.WriteLine(2 * (mesh.NodesAmountR - 1 + mesh.NodesAmountZ - 1));
        
        foreach (var border in mesh.borders)
        {
            // TODO: Указывать через нормаль switch(*normal*)
            if (border[2] == border[3]) // border || oZ
            {
                int iter = border[2] == 0 ? 0 : mesh.NodesAmountR - 1;
                for (int i = 0; i < mesh.NodesAmountZ - 1; i++)
                {
                    sw.WriteLine($"{border[0]} {border[1]} {iter} {iter + mesh.NodesAmountR}");
                    iter += mesh.NodesAmountR;
                }
            }
            else if (border[4] == border[5]) // border || oR
            {
                int iter = border[4] == 0 ? 0 : mesh.NodesAmountR * (mesh.NodesAmountZ - 1);
                for (int i = 0; i < mesh.NodesAmountR - 1; i++)
                {
                    sw.WriteLine($"{border[0]} {border[1]} {iter} {iter + 1}");
                    iter++;
                }

            }
        }
    }
}