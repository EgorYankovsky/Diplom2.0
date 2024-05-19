using MathObjects;
using Grid;
using DataStructs;
using Functions;


namespace Project;

public class FEM3D : FEM
{
    private static readonly double mu0 = 4.0D * Math.PI * Math.Pow(10.0D, -7);

    public ArrayOfPoints3D pointsArr;

    public List<GlobalVector> A;

    public List<GlobalVector> E;

    public List<GlobalVector> B;

    public List<GlobalVector> H;

    public ArrayOfRibs? ribsArr;

    private readonly Mesh3Dim mesh3Dim;

    private GlobalMatrix? G;

    private readonly FEM3D? _originalFEM;

    private List<GlobalVector>? _originalE;

    private GlobalMatrix? M;

    internal static List<int> ConvertGlobalToLocalNumeration(List<int> global) => 
            [global[0], global[3], global[8], global[11],
             global[1], global[2], global[9], global[10],
             global[4], global[5], global[6], global[7]];

    public FEM3D(Mesh3Dim mesh, TimeMesh timeMesh) : base(timeMesh)
    {
        string pointsPath = _3dValuesPath + $"AfterConvertation\\Points.poly";
        string elemsPath = _3dValuesPath + $"AfterConvertation\\Elems.poly";
        string bordersPath = _3dValuesPath + $"AfterConvertation\\Borders.poly";
        
        pointsArr = new ArrayOfPoints3D(pointsPath);
        elemsArr = new(elemsPath, 3);
        bordersArr = new (bordersPath, 3);
        ribsArr = mesh.arrayOfRibs;

        equationType = timeMesh[0] == timeMesh[^1] ? EquationType.Elliptic : EquationType.Parabolic;
        
        A = [];
        E = [];
        B = [];
        H = [];

        mesh3Dim = mesh;
    }

    public FEM3D(Mesh3Dim mesh, TimeMesh timeMesh, FEM3D originalFEM, int Num) : base(timeMesh)
    {
        string pointsPath = _3dValuesPath + $"Anomaly{Num}\\Points.poly";
        string elemsPath = _3dValuesPath + $"Anomaly{Num}\\Elems.poly";
        //string bordersPath = _3dValuesPath + $"Anomaly{Num}\\Borders.poly";

        pointsArr = new(pointsPath);
        elemsArr = new(elemsPath, 3);
        //bordersArr = new(bordersPath, 3);
        ribsArr = mesh.arrayOfRibs;
        _originalFEM = originalFEM;
        _originalE = [];

        equationType = timeMesh[0] == timeMesh[^1] ? EquationType.Elliptic : EquationType.Parabolic;
        
        A = [];
        E = [];
        B = [];
        H = [];
        mesh3Dim = mesh;
        ConstructMatrixes();
    }

    public (double, double, double) GetAAt(double x, double y, double z, double t)
    {
        if (ribsArr is null) throw new ArgumentNullException("Array of ribs not generated");
        for (int tt = 0; tt < Time.Count; tt++)
            if (tt == t)
            {
                var elem = ConvertGlobalToLocalNumeration(GetElem(x, y, z));
                double[] q = new double[12];
                for (int i = 0; i < 12; i++)
                    q[i] = A[tt][elem[i]];

                double x0 = ribsArr[elem[0]].a.X;
                double x1 = ribsArr[elem[0]].b.X;
                double y0 = ribsArr[elem[4]].a.Y;
                double y1 = ribsArr[elem[4]].b.Y;
                double z0 = ribsArr[elem[8]].a.Z;
                double z1 = ribsArr[elem[8]].b.Z;

                double eps = (x - x0) / (x1 - x0);
                double nu =  (y - y0) / (y1 - y0);
                double khi = (z - z0) / (z1 - z0);
                return BasisFunctions3DVec.GetValue(eps, nu, khi, q);
            }
        throw new Exception("Out of mesh borders");
    }

    public (double, double, double) GetEAt(double x, double y, double z, double t)
    {
        if (ribsArr is null) throw new ArgumentNullException("Array of ribs not generated");
        for (int tt = 0; tt < Time.Count; tt++)
            if (Time[tt] == t)
            {
                var elem = ConvertGlobalToLocalNumeration(GetElem(x, y, z));
                double[] q = new double[12];
                for (int i = 0; i < 12; i++)
                    q[i] = E[tt][elem[i]];

                double x0 = ribsArr[elem[0]].a.X;
                double x1 = ribsArr[elem[0]].b.X;
                double y0 = ribsArr[elem[4]].a.Y;
                double y1 = ribsArr[elem[4]].b.Y;
                double z0 = ribsArr[elem[8]].a.Z;
                double z1 = ribsArr[elem[8]].b.Z;

                double eps = (x - x0) / (x1 - x0);
                double nu =  (y - y0) / (y1 - y0);
                double khi = (z - z0) / (z1 - z0);
                return BasisFunctions3DVec.GetValue(eps, nu, khi, q);
            }
        throw new Exception("Out of mesh borders");
    }

    public (double, double, double) GetBAt(double x, double y, double z, double t)
    {
        throw new Exception("");
    }

    public (double, double, double) GetHAt(double x, double y, double z, double t)
    {
        throw new Exception("");
    }

    public void GenerateVectorB()
    {
        if (A is []) throw new Exception("Vector B isn't generated");
        if (ribsArr is null) throw new Exception("Array of ribs isn't generated");
        

        B = new(A.Count);
        for (int t = 0; t < Time.Count; t++)
        {
            B[t] = new GlobalVector(ribsArr.Count);
            for (int i = 0; i < A[t].Size; i++)
            {
                double ht = Math.Pow(10, -10);
                var pnt = ribsArr[i].GetMiddlePoint();
                
                var rib_1 = GetAAt(pnt.X + ht, pnt.Y + ht, pnt.Z + ht, Time[t]);
                var rib_0 = GetAAt(pnt.X - ht, pnt.Y - ht, pnt.Z - ht, Time[t]);

                //var vec1 = GetAAt();
            }
        }
    }

    public void GenerateVectorH()
    {
        if (B is []) throw new Exception("Vector B isn't generated");
        H = new(B.Count);
        for (int t = 0; t < B.Count; t++)
        {
            H.Add(new GlobalVector(B[t].Size));
            for (int i = 0; i < B[i].Size; i++)
                H[t][i] = B[t][i] / mu0;

        }
    }

    public void GenerateVectorE()
    {
        E = new(A.Count);   
        for (int i = 0; i < A.Count; i++)
        {
            if (i == 0 || i == 1)
                E.Add(new GlobalVector(A[i].Size));
            else
            {
                double ti = Time[i];
                double ti_1 = Time[i - 1];
                double ti_2 = Time[i - 2];
                E.Add(-1.0D / (ti_1 - ti_2) * A[i - 2] + (ti - ti_2) / ((ti_1 - ti_2) * (ti - ti_1)) * A[i - 1] - 
                (2 * ti - ti_1 - ti_2) / ((ti - ti_2) * (ti - ti_1)) * A[i]);
            }
        }
    }

    public void AddSolution(FEM3D fem)
    {
        if (ribsArr is null) throw new ArgumentNullException("ribsArr is null");
        if (fem.ribsArr is null) throw new ArgumentNullException("ribsArr is null");
        for (int t = 0; t < Time.Count; t++)
        {
            for (int i = 0; i < ribsArr.Count; i++)
            {
                var X = 0.5 * (ribsArr[i].a.X + ribsArr[i].b.X);
                var Y = 0.5 * (ribsArr[i].a.Y + ribsArr[i].b.Y);
                var Z = 0.5 * (ribsArr[i].a.Z + ribsArr[i].b.Z);

                var antinormal = ((ribsArr[i].b.X - ribsArr[i].a.X) / ribsArr[i].Length,
                                  (ribsArr[i].b.Y - ribsArr[i].a.Y) / ribsArr[i].Length,
                                  (ribsArr[i].b.Z - ribsArr[i].a.Z) / ribsArr[i].Length);

                var elem = fem.GetElem(X, Y, Z);

                int[] elem_local = [elem[0], elem[3], elem[8], elem[11], 
                                    elem[1], elem[2], elem[9], elem[10], 
                                    elem[4], elem[5], elem[6], elem[7]];

                double[] q = [fem.A[t][elem[0]], fem.A[t][elem[3]], fem.A[t][elem[8]], fem.A[t][elem[11]], 
                              fem.A[t][elem[1]], fem.A[t][elem[2]], fem.A[t][elem[9]], fem.A[t][elem[10]], 
                              fem.A[t][elem[4]], fem.A[t][elem[5]], fem.A[t][elem[6]], fem.A[t][elem[7]]];
                
                var eps = (X - ribsArr[elem_local[0]].a.X) / (ribsArr[elem_local[0]].b.X - ribsArr[elem_local[0]].a.X);
                var nu =  (Y - ribsArr[elem_local[4]].a.Y) / (ribsArr[elem_local[4]].b.Y - ribsArr[elem_local[4]].a.Y);
                var khi = (Z - ribsArr[elem_local[8]].a.Z) / (ribsArr[elem_local[8]].b.Z - ribsArr[elem_local[8]].a.Z);

                var val = BasisFunctions3DVec.GetValue(eps, nu, khi, q);
                
                A[t][i] += val.Item1 * antinormal.Item1 + val.Item2 * antinormal.Item2 + val.Item3 * antinormal.Item3;
            }
        }
    }

    internal List<int> GetElem(double x, double y, double z)
    {
        int i = 0;
        for (; i < mesh3Dim.nodesX.Count - 1; i++)
            if (mesh3Dim.nodesX[i] <= x && x <= mesh3Dim.nodesX[i + 1])
                break;

        int j = 0;
        for (; j < mesh3Dim.nodesY.Count - 1; j++)
            if (mesh3Dim.nodesY[j] <= y && y <= mesh3Dim.nodesY[j + 1])
                break;

        int k = 0;
        for (; k < mesh3Dim.nodesZ.Count - 1; k++)
            if (mesh3Dim.nodesZ[k] <= z && z <= mesh3Dim.nodesZ[k + 1])
                break;

        return elemsArr[k * (mesh3Dim.nodesX.Count - 1) * (mesh3Dim.nodesY.Count - 1) + j * (mesh3Dim.nodesX.Count - 1) + i].Arr;
    }

    public void ConvertResultTo3Dim(FEM2D fem2d)
    {
        if (fem2d.pointsArr is null) throw new ArgumentNullException();
        if (ribsArr is null) throw new ArgumentNullException();
        
        for (int i = 0; i < Time.Count; i++)
        {
            A.Add(new GlobalVector(ribsArr.Count));
            E.Add(new GlobalVector(ribsArr.Count));
            
            for (int j = 0; j < ribsArr.Count; j++)
            {
                var X = 0.5D * (ribsArr[j].a.X + ribsArr[j].b.X);
                var Y = 0.5D * (ribsArr[j].a.Y + ribsArr[j].b.Y);
                var Z = 0.5D * (ribsArr[j].a.Z + ribsArr[j].b.Z);

                var antinormal = ((ribsArr[j].b.X - ribsArr[j].a.X) / ribsArr[j].Length,
                                  (ribsArr[j].b.Y - ribsArr[j].a.Y) / ribsArr[j].Length,
                                  (ribsArr[j].b.Z - ribsArr[j].a.Z) / ribsArr[j].Length);

                var elem = fem2d.GetElem(Math.Sqrt(X * X + Y * Y), Z);

                var fa = BasisFunctions2D.GetValue(
                                fem2d.A_phi[i][elem[0]], fem2d.A_phi[i][elem[1]],
                                fem2d.A_phi[i][elem[2]], fem2d.A_phi[i][elem[3]],
                                fem2d.pointsArr[elem[0]].R, fem2d.pointsArr[elem[1]].R,
                                fem2d.pointsArr[elem[0]].Z, fem2d.pointsArr[elem[3]].Z, 
                                Math.Sqrt(X * X + Y * Y), Z);

                var fe = BasisFunctions2D.GetValue(
                                fem2d.E_phi[i][elem[0]], fem2d.E_phi[i][elem[1]],
                                fem2d.E_phi[i][elem[2]], fem2d.E_phi[i][elem[3]],
                                fem2d.pointsArr[elem[0]].R, fem2d.pointsArr[elem[1]].R,
                                fem2d.pointsArr[elem[0]].Z, fem2d.pointsArr[elem[3]].Z, 
                                Math.Sqrt(X * X + Y * Y), Z);

                //if (fe != 0.0D && (antinormal.Item1 != 0 || antinormal.Item2 != 0))
                //    Console.WriteLine();

                var Ax = X == 0 && Y == 0 ? 0.0D : -1.0D * (Y / Math.Sqrt(X * X + Y * Y)) * fa;
                var Ay = X == 0 && Y == 0 ? 0.0D : X / Math.Sqrt(X * X + Y * Y) * fa;
                var Az = 0.0D;

                var Ex = X == 0 && Y == 0 ? 0.0D :  -1.0D * (Y / Math.Sqrt(X * X + Y * Y)) * fe;
                var Ey = X == 0 && Y == 0 ? 0.0D :  X / Math.Sqrt(X * X + Y * Y) * fe;
                var Ez = 0.0D;

                A[i][j] = Ax * antinormal.Item1 + Ay * antinormal.Item2 + Az * antinormal.Item3;
                E[i][j] = Ex * antinormal.Item1 + Ey * antinormal.Item2 + Ez * antinormal.Item3;
            }
        }
    }

    public void CheckSolution(List<Point3D> recivers)
    {
        GenerateVectorB();
        GenerateVectorH();
        for (int t = 0; t < Time.Count; t++)
        {
            foreach (var reciver in recivers)
            {

            }
        }
    }

    public void ConstructMatrixes()
    {
        if (ribsArr is null) throw new ArgumentNullException("ribsArr is null");

        var sparceMatrix = new GlobalMatrix(ribsArr.Count);
        Generator.BuildPortait(ref sparceMatrix, ribsArr.Count, elemsArr);

        G = new GlobalMatrix(sparceMatrix);
        Generator.FillMatrixG(ref G, ribsArr, elemsArr);
        
        M = new GlobalMatrix(sparceMatrix);
        Generator.FillMatrixM(ref M, ribsArr, elemsArr);
    }

    public void Solve()
    {
        if (solver is null) throw new ArgumentNullException("Solver is null");
        if (ribsArr is null) throw new ArgumentNullException("ribs array is null");
        if (_originalFEM is null) throw new ArgumentNullException("original E is null");
        if (G is null) throw new ArgumentNullException();
        if (M is null) throw new ArgumentNullException();
        
        Solutions = new GlobalVector[Time.Count];
        Discrepancy = new GlobalVector[Time.Count];
        
        (Solutions[0], Discrepancy[0]) = (new GlobalVector(ribsArr.Count), new GlobalVector(ribsArr.Count));

        if (Time.Count > 1)
        {
            (Solutions[1], Discrepancy[1]) = (Solutions[0], Discrepancy[0]);
            for (int i = 1; i < Time.Count; i++)
            {
                double t0 = Time[i];
                double t1 = Time[i - 1];

                double deltT = t0 - t1;

                double tau0 = 1.0D / deltT;
                
                Matrix = G + tau0 * M;
                
                var b = new GlobalVector(ribsArr.Count);
                
                // ! ACHTUNG
                Generator.FillVector3D(ref b, _originalFEM.GetEAt, ribsArr, elemsArr, Time[i]);

                Vector = b + tau0 * M * Solutions[i - 1];
                
                //Generator.ConsiderBoundaryConditions(ref Matrix, ref Vector, ribsArr, bordersArr, Time[i]);
                (Solutions[i], Discrepancy[i]) = solver.Solve(Matrix, Vector);
            }
        }
        A = [.. Solutions];
    }

    public void WriteData(string path)
    {
        if (Solutions is null) throw new ArgumentNullException("No solutions");
        for (int t = 0; t < Time.Count; t++)
        {
            using var sw = new StreamWriter(path + $"Answer_{Time[t]}.txt");
            for (int i = 0; i < A[t].Size; i++)
                sw.WriteLine($"{i} {A[t][i]:E8}");
            sw.Close();
        }
    }

    public void WriteDataToDraw(string path)
    {
        double hx = (mesh3Dim.nodesX[^1] - mesh3Dim.nodesX[0]) / 10.0D;
        double hy = (mesh3Dim.nodesY[^1] - mesh3Dim.nodesY[0]) / 10.0D;
        double hz = (mesh3Dim.nodesZ[^1] - mesh3Dim.nodesZ[0]) / 10.0D;
        
        double x0 = mesh3Dim.nodesX[0];
        double y0 = mesh3Dim.nodesY[0];
        double z0 = mesh3Dim.nodesZ[0];
        List<Point3D> points = [];
        int i = 0;
        int j = 0;
        int k = 0;

        while (z0 + hz * k <= mesh3Dim.nodesZ[^1])
        {
            j = 0;
            while (y0 + hy * j <= mesh3Dim.nodesY[^1])
            {
                i = 0;
                while (x0 + hx * i <= mesh3Dim.nodesX[^1])
                {
                    points.Add(new Point3D(x0 + i * hx, y0 + j * hy, z0 + k * hz));
                    i++;
                }
                j++;
            }
            k++;
        }

        for (int t = 0; t < Time.Count; t++)
        {
            using var sw = new StreamWriter(path + $"Answer{t}.txt");
            sw.WriteLine($"{mesh3Dim.nodesX[0]:E15} {mesh3Dim.nodesX[^1]:E15} {mesh3Dim.nodesY[0]:E15} {mesh3Dim.nodesY[^1]:E15} {mesh3Dim.nodesZ[0]:E15} {mesh3Dim.nodesZ[^1]:E15}");
            foreach (var point in points)
            {
                var ans = GetEAt(point.X, point.Y, point.Z, Time[t]);
                sw.WriteLine($"{point.X:E15} {point.Y:E15} {point.Z:E15} {ans.Item1:E15} {ans.Item2:E15} {ans.Item3:E15}");
            }
            sw.Close();
        }

    }

    public void TestOutput(string path)
    {
        using var sw = new StreamWriter(path + "/A_phi/Answer3D/Answer_Test.txt");
        
        var absDiscX = 0.0D;
        var absDiscY = 0.0D;
        var absDiscZ = 0.0D;

        var absDivX = 0.0D;
        var absDivY = 0.0D;
        var absDivZ = 0.0D;

        var relDiscX = 0.0D;
        var relDiscY = 0.0D;
        var relDiscZ = 0.0D;

        var relDivX = 0.0D;
        var relDivY = 0.0D;
        var relDivZ = 0.0D;

        int iter = 0;

        foreach (Elem elem in elemsArr)
        {
            int[] elem_local = [elem[0], elem[3], elem[8], elem[11],
                                elem[1], elem[2], elem[9], elem[10],
                                elem[4], elem[5], elem[6], elem[7]];
            
            var x = 0.5D * (ribsArr[elem_local[0]].a.X + ribsArr[elem_local[0]].b.X);
            var y = 0.5D * (ribsArr[elem_local[4]].a.Y + ribsArr[elem_local[4]].b.Y);
            var z = 0.5D * (ribsArr[elem_local[8]].a.Z + ribsArr[elem_local[8]].b.Z);

            sw.WriteLine($"Points {x:E15} {y:E15} {z:E15}");
            
            var eps = (x - ribsArr[elem_local[0]].a.X) / (ribsArr[elem_local[0]].b.X - ribsArr[elem_local[0]].a.X);
            var nu =  (y - ribsArr[elem_local[4]].a.Y) / (ribsArr[elem_local[4]].b.Y - ribsArr[elem_local[4]].a.Y);
            var khi = (z - ribsArr[elem_local[8]].a.Z) / (ribsArr[elem_local[8]].b.Z - ribsArr[elem_local[8]].a.Z);

            double[] q = [Solutions[0][elem_local[0]], Solutions[0][elem_local[1]], Solutions[0][elem_local[2]], Solutions[0][elem_local[3]], 
                          Solutions[0][elem_local[4]], Solutions[0][elem_local[5]], Solutions[0][elem_local[6]], Solutions[0][elem_local[7]],
                          Solutions[0][elem_local[8]], Solutions[0][elem_local[9]], Solutions[0][elem_local[10]], Solutions[0][elem_local[11]]];
            
            var ans = BasisFunctions3DVec.GetValue(eps, nu, khi, q);
            var theorValue = Function.A(x, y, z, 0.0D);

            sw.WriteLine($"FEM A  {ans.Item1:E15} {ans.Item2:E15} {ans.Item3:E15}");
            sw.WriteLine($"Theor  {theorValue.Item1:E15} {theorValue.Item2:E15} {theorValue.Item3:E15}");


            var currAbsDiscX = Math.Abs(ans.Item1 - theorValue.Item1);
            var currAbsDiscY = Math.Abs(ans.Item2 - theorValue.Item2);
            var currAbsDiscZ = Math.Abs(ans.Item3 - theorValue.Item3);

            var currRelDiscX = currAbsDiscX / Math.Abs(theorValue.Item1);
            var currRelDiscY = currAbsDiscY / Math.Abs(theorValue.Item2);
            var currRelDiscZ = currAbsDiscZ / Math.Abs(theorValue.Item3);

            sw.WriteLine($"CurrAD {currAbsDiscX:E15} {currAbsDiscY:E15} {currAbsDiscZ:E15}");
            sw.WriteLine($"CurrRD {currRelDiscX:E15} {currRelDiscY:E15} {currRelDiscZ:E15}\n");

            absDiscX += currAbsDiscX;
            absDiscY += currAbsDiscY;
            absDiscZ += currAbsDiscZ;

            absDivX += theorValue.Item1;
            absDivY += theorValue.Item2;
            absDivZ += theorValue.Item3;

            relDiscX += currRelDiscX * currRelDiscX;
            relDiscY += currRelDiscY * currRelDiscY;
            relDiscZ += currRelDiscZ * currRelDiscZ;

            relDivX += theorValue.Item1 * theorValue.Item1;
            relDivY += theorValue.Item2 * theorValue.Item2;
            relDivZ += theorValue.Item3 * theorValue.Item3;

            iter++;
        }
        sw.WriteLine($"Avg disc: {absDiscX / iter:E15} {absDiscY / iter:E15} {absDiscZ / iter:E15}");
        sw.WriteLine($"Rel disc: {Math.Sqrt(relDiscX / relDivX):E15} {Math.Sqrt(relDiscY / relDivY):E15} {Math.Sqrt(relDiscZ / relDivZ):E15}");
        sw.Close();
    }
}