using System.Collections.Immutable;
using System.Numerics;
using MathObjects;
using Solver;
using Grid;
using DataStructs;
using System.Diagnostics;
using Functions;
using System.ComponentModel.DataAnnotations;
using System.Timers;

namespace Project;

public class FEM3D : FEM
{
    public ArrayOfPoints3D pointsArr = new(_pointspath3D);

    public List<GlobalVector> A;

    public List<GlobalVector> E;

    public List<List<GlobalVector>> A_plus;

    public ArrayOfRibs? ribsArr;

    private readonly Mesh3Dim mesh3Dim;

    // Maybe private?
    public List<Layer> Layers;

    private Layer _currentLayer;

    private GlobalMatrix? G;

    private GlobalMatrix? M;

    public FEM3D(Mesh3Dim mesh, TimeMesh timeMesh, List<Layer> layers) : base(timeMesh, 3)
    {
        equationType = timeMesh[0] == timeMesh[^1] ? EquationType.Elliptic : EquationType.Parabolic;
        
        A = [];
        E = [];
        A_plus = new List<List<GlobalVector>>(layers.Count);
        ribsArr = mesh.arrayOfRibs;

        Layers = layers;
        _currentLayer = Layers[0];
        //MeshGenerator.SelectRibs(ref ribsArr, ref elemsArr);
        mesh3Dim = mesh;
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

    public void SelectCurrentLayer(int num) => _currentLayer = Layers[num];

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

    public void AddField(Layer layer)
    {
        ArgumentNullException.ThrowIfNull(layer);
        Layers?.Add(layer);
    }

    public void CommitFields()
    {
        if (elemsArr is null) throw new ArgumentNullException("elemsArr is null");

        foreach (var layer in Layers)
        {
            for (int i = 0; i < elemsArr.Length; i++)
            {
                double minz = Math.Min(ribsArr[elemsArr[i][^1]].a.Z, ribsArr[elemsArr[i][^1]].b.Z);
                double maxz = Math.Max(ribsArr[elemsArr[i][^1]].a.Z, ribsArr[elemsArr[i][^1]].b.Z);
                if (layer.z0 <= minz && maxz <= layer.z1)
                    elemsArr[i].sigma = layer.sigma;
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
        if (G is null) throw new ArgumentNullException();
        if (M is null) throw new ArgumentNullException();
        
        Solutions = new GlobalVector[Time.Count];
        Discrepancy = new GlobalVector[Time.Count];
        
        (Solutions[0], Discrepancy[0]) = (new GlobalVector(ribsArr.Count), new GlobalVector(ribsArr.Count));

        if (Time.Count > 1)
        {
            (Solutions[1], Discrepancy[1]) = (Solutions[0], Discrepancy[0]);
            for (int i = 2; i < Time.Count; i++)
            {
                double t0 = Time[i];
                double t1 = Time[i - 1];
                double t2 = Time[i - 2];

                double deltT = t0 - t2;
                double deltT0 = t0 - t1;
                double deltT1 = t1 - t2;

                double tau0 = (deltT + deltT0) / (deltT * deltT0);
                double tau1 = deltT / (deltT1 * deltT0);
                double tau2 = deltT0 / (deltT * deltT1);

                Matrix = G + tau0 * M;
                
                var b = new GlobalVector(ribsArr.Count);
                Generator.FillVector3D(ref b, _currentLayer, ribsArr, elemsArr, Time[i]);

                Vector = b - tau2 * M * Solutions[i - 2] + tau1 * M * Solutions[i - 1];
                
                Generator.ConsiderBoundaryConditions(ref Matrix, ref Vector, ribsArr, bordersArr, Time[i]);
                (Solutions[i], Discrepancy[i]) = solver.Solve(Matrix, Vector);
            }
        }
        A = [.. Solutions];
    }


    public void CommitField(int layer)
    {
        if (A.Count != A_plus[layer].Count)
            throw new ArgumentOutOfRangeException("Different sizes");
        for (int i = 0; i < A.Count; i++)
            A[i] += A_plus[layer][i];
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