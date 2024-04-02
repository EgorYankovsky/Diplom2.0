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
        var arr = new ArrayOfPoints(mesh.NodesAmountTotal);
        foreach (var Z in mesh.nodesZ)
            foreach (var Y in mesh.nodesY)
                foreach(var X in mesh.nodesX)
                    arr.Append(new Point(X, Y, Z));
        return arr;
    }

    public static ArrayOfElems GenerateListOfElems(Mesh mesh)
    {
        var arr = new ArrayOfElems(mesh.ElemsAmount);
        for (int k = 0; k < mesh.NodesAmountZ - 1; k++)
            for (int j = 0; j < mesh.NodesAmountY - 1; j++)
                for (int i = 0; i < mesh.NodesAmountX - 1; i++)
                    arr.Add(new List<int>{k * mesh.NodesAmountY * mesh.NodesAmountX + j * mesh.NodesAmountX + i, k * mesh.NodesAmountY * mesh.NodesAmountX + j * mesh.NodesAmountX + i + 1,
                                          k * mesh.NodesAmountY * mesh.NodesAmountX + (j + 1) * mesh.NodesAmountX + i, k * mesh.NodesAmountY * mesh.NodesAmountX + (j + 1) * mesh.NodesAmountX + i + 1,
                                          (k + 1) * mesh.NodesAmountY * mesh.NodesAmountX + j * mesh.NodesAmountX + i, (k + 1) * mesh.NodesAmountY * mesh.NodesAmountX + j * mesh.NodesAmountX + i + 1,
                                          (k + 1) * mesh.NodesAmountY * mesh.NodesAmountX + (j + 1) * mesh.NodesAmountX + i, (k + 1) * mesh.NodesAmountY * mesh.NodesAmountX + (j + 1) * mesh.NodesAmountX + i + 1});
        
        return arr;
    }

    public static ArrayOfBorders GenerateListOfBorders(Mesh mesh)
    {
        var arr = new ArrayOfBorders(null);

        // XY0
        for (int j = 0; j < mesh.NodesAmountY - 1; j++)
            for (int i = 0; i < mesh.NodesAmountX - 1; i++)
                arr.Arr.Add(new List<int> {1, 1, j * mesh.NodesAmountX + i, j * mesh.NodesAmountX + i + 1,
                                                (j + 1) * mesh.NodesAmountX + i, (j + 1) * mesh.NodesAmountX + i + 1});
        
        // X0Z
        for (int j = 0; j < mesh.NodesAmountZ - 1; j++)
            for (int i = 0; i < mesh.NodesAmountX - 1; i++)
                arr.Arr.Add(new List<int> {1, 1, j * mesh.NodesAmountY * mesh.NodesAmountX + i, j * mesh.NodesAmountY * mesh.NodesAmountX + i + 1,
                                                (j + 1) * mesh.NodesAmountY * mesh.NodesAmountX + i, (j + 1) * mesh.NodesAmountY * mesh.NodesAmountX + i + 1});
        
        // X1Z
        for (int j = 0; j < mesh.NodesAmountZ - 1; j++)
            for (int i = 0; i < mesh.NodesAmountX - 1; i++)
                arr.Arr.Add(new List<int> {1, 1, j * mesh.NodesAmountY * mesh.NodesAmountX + mesh.NodesAmountY + i, j * mesh.NodesAmountY * mesh.NodesAmountX + mesh.NodesAmountY + i + 1,
                                                (j + 1) * mesh.NodesAmountY * mesh.NodesAmountX + mesh.NodesAmountY + i, (j + 1) * mesh.NodesAmountY * mesh.NodesAmountX + mesh.NodesAmountY + i + 1});
        // 0YZ
        for (int j = 0; j < mesh.NodesAmountZ - 1; j++)
            for (int i = 0; i < mesh.NodesAmountY - 1; i++)
                arr.Arr.Add(new List<int> {1, 1, j * mesh.NodesAmountY * mesh.NodesAmountX + i * mesh.NodesAmountX, j * mesh.NodesAmountY * mesh.NodesAmountX + (i + 1) * mesh.NodesAmountX, 
                                                (j + 1) * mesh.NodesAmountY * mesh.NodesAmountX + j * mesh.NodesAmountX, (j + 1) * mesh.NodesAmountY * mesh.NodesAmountX + (i + 1) * mesh.NodesAmountX});
        
        // 1YZ
        for (int j = 0; j < mesh.NodesAmountZ - 1; j++)
            for (int i = 0; i < mesh.NodesAmountY - 1; i++)
                arr.Arr.Add(new List<int> {1, 1, j * mesh.NodesAmountY * mesh.NodesAmountX + i * mesh.NodesAmountX + mesh.NodesAmountX - 1,
                                                 j * mesh.NodesAmountY * mesh.NodesAmountX + (i + 1) * mesh.NodesAmountX + mesh.NodesAmountX - 1, 
                                                (j + 1) * mesh.NodesAmountY * mesh.NodesAmountX + j * mesh.NodesAmountX + mesh.NodesAmountX - 1,
                                                (j + 1) * mesh.NodesAmountY * mesh.NodesAmountX + (i + 1) * mesh.NodesAmountX + mesh.NodesAmountX - 1});
        
        // XY1
        for (int j = 0; j < mesh.NodesAmountY - 1; j++)
            for (int i = 0; i < mesh.NodesAmountX - 1; i++)
                arr.Arr.Add(new List<int> {2, 1, (mesh.NodesAmountZ - 1) * j * mesh.NodesAmountX + i, (mesh.NodesAmountZ - 1) * j * mesh.NodesAmountX + i + 1,
                                                 (mesh.NodesAmountZ - 1) * (j + 1) * mesh.NodesAmountX + i, (mesh.NodesAmountZ - 1) * (j + 1) * mesh.NodesAmountX + i + 1});
        
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
            index = 2;
        else if (mesh.NodesZWithoutFragmentation[3] <= a.Z && d.Z <= mesh.NodesZWithoutFragmentation[4])
            index = 3;
        else if (mesh.NodesZWithoutFragmentation[4] <= a.Z && d.Z <= mesh.NodesZWithoutFragmentation[5])
            index = 4;
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