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
    public List<GlobalVector> Ax_3D;

    public List<GlobalVector> Ay_3D;

    public List<GlobalVector> Az_3D;

    public List<GlobalVector> Ex_3D;

    public List<GlobalVector> Ey_3D;

    public List<GlobalVector> Ez_3D;

    public ArrayOfRibs ribsArr;

    // Maybe private?
    public List<Layer> Layers;

    public FEM3D(FEM2D fem2d) : base(fem2d.Time, 3)
    {
        Ax_3D = [];
        Ay_3D = [];
        Az_3D = [];
        Ex_3D = [];
        Ey_3D = [];
        Ez_3D = [];
        Layers = [];
        //mesh = new Mesh3Dim
        //{
        //    nodesX = [],
        //    nodesY = [],
        //    nodesZ = []
        //};
        equationType = fem2d.equationType;
    }

    public void ConstructMesh()
    {
        if (mesh == null) throw new ArgumentNullException("Mesh is null");
        mesh.nodesX = [0.0D, 1.0D, 2.0D];
        mesh.nodesY = [0.0D, 1.0D, 2.0D];
        mesh.nodesZ = [0.0D, 1.0D, 2.0D, 3.0D];
        timeMesh = [1.0D];
    }

    public void ConstructMesh(FEM2D fem2d)
    {   
        //if (fem2d.Mesh2D.nodesR is null) throw new ArgumentNullException("null object");
        if (mesh is null) throw new ArgumentNullException("null object");
        //if (mesh.nodesX is null) throw new ArgumentNullException("null object");
        //if (mesh.nodesY is null) throw new ArgumentNullException("null object");

        //mesh.nodesX = fem2d.Mesh2D.nodesR.Where(r => r < fem2d.Mesh2D.nodesR.Last() / Math.Sqrt(2.0D)).ToList();
        //mesh.nodesY = fem2d.Mesh2D.nodesR.Where(r => r < fem2d.Mesh2D.nodesR.Last() / Math.Sqrt(2.0D)).ToList();
        
        //mesh.nodesX.Add(fem2d.Mesh2D.nodesR.Last() / Math.Sqrt(2.0D));
        //mesh.nodesY.Add(fem2d.Mesh2D.nodesR.Last() / Math.Sqrt(2.0D));
        
        //int amount = 2 * mesh.nodesX.Count - 2;

        //for (int i = 1; i < amount; i += 2)
        //{
        //   mesh.nodesX.Insert(0, -1.0D * mesh.nodesX[i]);
        //   mesh.nodesY.Insert(0, -1.0D * mesh.nodesY[i]);
        //}

        //mesh.mu = fem2d.mu0;
        //mesh.sigma = fem2d.sigma;

        //mesh.nodesZ = fem2d.Mesh2D.nodesZ;
        //timeMesh = fem2d.timeMesh;
    }

    public void ConvertResultTo3Dim(FEM2D fem2d)
    {
        Task TaskGenerationAxyz = Task.Run(() => GenerateAxyz(fem2d));
        Task TaskGenerationExyz = Task.Run(() => GenerateExyz(fem2d));
        TaskGenerationAxyz.Wait();
        TaskGenerationExyz.Wait();
    }

    public void GenerateAxyz(FEM2D fem2d)
    {
        //if (timeMesh is null) throw new ArgumentNullException();
        if (mesh is null) throw new ArgumentNullException();
        //if (mesh.nodesX is null) throw new ArgumentNullException();
        //if (mesh.nodesY is null) throw new ArgumentNullException();
        //if (mesh.nodesZ is null) throw new ArgumentNullException();
        if (fem2d.pointsArr is null) throw new ArgumentNullException();

        int i = 0;
        
        //foreach (var t in )
        //{
        //    int j = 0;
        //    Ax_3D.Add(new GlobalVector(mesh.NodesAmountTotal));
        //    Ay_3D.Add(new GlobalVector(mesh.NodesAmountTotal));
        //    Az_3D.Add(new GlobalVector(mesh.NodesAmountTotal));
        //    /*
        //    foreach (var Z in mesh.nodesZ)
        //    {
        //        foreach (var Y in mesh.nodesY)
        //        {
        //            foreach (var X in mesh.nodesX)
        //            {
        //                var elem = fem2d.GetE_phi(Math.Sqrt(X * X + Y * Y), Z, t);
//
        //                Ax_3D[i][j] = -1.0D * (Y / Math.Sqrt(X * X + Y * Y)) * 
        //                    BasisFunctions2D.GetValue(
        //                        fem2d.A_phi[i][elem[0]], fem2d.A_phi[i][elem[1]],
        //                        fem2d.A_phi[i][elem[2]], fem2d.A_phi[i][elem[3]],
        //                        fem2d.pointsArr[elem[0]].R, fem2d.pointsArr[elem[1]].R,
        //                        fem2d.pointsArr[elem[0]].Z, fem2d.pointsArr[elem[3]].Z, 
        //                        Math.Sqrt(X * X + Y * Y), Z);
        //                
        //                Ay_3D[i][j] = X / Math.Sqrt(X * X + Y * Y) * 
        //                    BasisFunctions2D.GetValue(
        //                        fem2d.A_phi[i][elem[0]], fem2d.A_phi[i][elem[1]],
        //                        fem2d.A_phi[i][elem[2]], fem2d.A_phi[i][elem[3]],
        //                        fem2d.pointsArr[elem[0]].R, fem2d.pointsArr[elem[1]].R,
        //                        fem2d.pointsArr[elem[0]].Z, fem2d.pointsArr[elem[3]].Z, 
        //                        Math.Sqrt(X * X + Y * Y), Z);
        //                
        //                j++;
        //            }
        //        }
        //    }
        //    i++;
        //    */
        //}
    }

    private void GenerateExyz(FEM2D fem2d)
    {
        //if (timeMesh is null) throw new ArgumentNullException();
        if (mesh is null) throw new ArgumentNullException();
        //if (mesh.nodesX is null) throw new ArgumentNullException();
        //if (mesh.nodesY is null) throw new ArgumentNullException();
        //if (mesh.nodesZ is null) throw new ArgumentNullException();
        if (fem2d.pointsArr is null) throw new ArgumentNullException();

        int i = 0;
        
        //foreach (var t in timeMesh)
        //{
        //    if (t == timeMesh.Last())
        //        continue;
        //    int j = 0;
        //    Ex_3D.Add(new GlobalVector(mesh.NodesAmountTotal));
        //    Ey_3D.Add(new GlobalVector(mesh.NodesAmountTotal));
        //    Ez_3D.Add(new GlobalVector(mesh.NodesAmountTotal));
        //    /*
        //    foreach (var Z in mesh.nodesZ)
        //    {
        //        foreach (var Y in mesh.nodesY)
        //        {
        //            foreach (var X in mesh.nodesX)
        //            {
        //                var elem = fem2d.GetE_phi(Math.Sqrt(X * X + Y * Y), Z, t);
//
        //                Ex_3D[i][j] = -1.0D * (Y / Math.Sqrt(X * X + Y * Y)) * 
        //                    BasisFunctions2D.GetValue(
        //                        fem2d.E_phi2D[i][elem[0]], fem2d.E_phi2D[i][elem[1]],
        //                        fem2d.E_phi2D[i][elem[2]], fem2d.E_phi2D[i][elem[3]],
        //                        fem2d.pointsArr[elem[0]].R, fem2d.pointsArr[elem[1]].R,
        //                        fem2d.pointsArr[elem[0]].Z, fem2d.pointsArr[elem[3]].Z, 
        //                        Math.Sqrt(X * X + Y * Y), Z);
        //                
        //                Ey_3D[i][j] = X / Math.Sqrt(X * X + Y * Y) * 
        //                    BasisFunctions2D.GetValue(
        //                        fem2d.E_phi2D[i][elem[0]], fem2d.E_phi2D[i][elem[1]],
        //                        fem2d.E_phi2D[i][elem[2]], fem2d.E_phi2D[i][elem[3]],
        //                        fem2d.pointsArr[elem[0]].R, fem2d.pointsArr[elem[1]].R,
        //                        fem2d.pointsArr[elem[0]].Z, fem2d.pointsArr[elem[3]].Z, 
        //                        Math.Sqrt(X * X + Y * Y), Z);
        //                
        //                j++;
        //            }
        //        }
        //    }
        //    i++;
        //    */
        //}
    }

    public void GenerateArrays()
    {
        if (mesh is null) throw new ArgumentNullException("mesh is null!");
        //pointsArr = MeshGenerator.GenerateListOfPoints(mesh);
        //ribsArr = MeshGenerator.GenerateListOfRibs(mesh, pointsArr);
        //elemsArr = MeshGenerator.GenerateListOfElems(mesh);
        //MeshGenerator.SelectRibs(ref ribsArr, ref elemsArr);
    }

    public void AddField(Layer layer)
    {
        ArgumentNullException.ThrowIfNull(layer);
        Layers.Add(layer);
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
                    elemsArr.sigmai[i] = layer.sigma;
            }
        }
    }

    public void ConstructMatrixAndVector()
    {
        if (elemsArr is null) throw new ArgumentNullException("elemsArr is null");
        if (bordersArr is null) throw new ArgumentNullException("bordersArr is null");
        
        var sparceMatrix = new GlobalMatrix(ribsArr.Count);
        Generator.BuildPortait(ref sparceMatrix, ribsArr.Count, elemsArr);

        var G = new GlobalMatrix(sparceMatrix);
        Generator.FillMatrixG(ref G, ribsArr, elemsArr);
        
        var M = new GlobalMatrix(sparceMatrix);
        Generator.FillMatrixM(ref M, ribsArr, elemsArr);
        
        Matrix = G + M;

        var b = new GlobalVector(ribsArr.Count);
        Generator.FillVector3D(ref b, ribsArr, elemsArr, 0.0);
        Vector = b;

        Generator.ConsiderBoundaryConditions(ref Matrix, ref Vector, ribsArr, bordersArr, 0.0D);
    }

    public void Solve()
    {
        if (solver is null) throw new ArgumentNullException("Solver is null");
        if (Matrix is null) throw new ArgumentNullException("Matrix is null");
        if (Vector is null) throw new ArgumentNullException("Vector is null");
        Solutions = new GlobalVector[timeMesh.Length];
        Discrepancy = new GlobalVector[timeMesh.Length];
        (Solutions[0], Discrepancy[0]) = solver.Solve(Matrix, Vector);
        //
        //if (timeMesh.Length > 1)
        //{
        //    for (int i = 0; i < timeMesh.Length; i++)
        //    {
        //        if (i == 1 || i == 0)
        //        {
        //            Solutions[i] = new GlobalVector(ribsArr.Count);
        //            for (int j = 0; j < ribsArr.Count; j++)
        //            {
        //                var antinormal = ((ribsArr[j].b.X - ribsArr[j].a.X) / ribsArr[j].Length,
        //                                  (ribsArr[j].b.Y - ribsArr[j].a.Y) / ribsArr[j].Length,
        //                                  (ribsArr[j].b.Z - ribsArr[j].a.Z) / ribsArr[j].Length);
        //                var pointat = ((ribsArr[j].b.X + ribsArr[j].a.X) / 2.0D,
        //                               (ribsArr[j].b.Y + ribsArr[j].a.Y) / 2.0D,
        //                               (ribsArr[j].b.Z + ribsArr[j].a.Z) / 2.0D);

        //                var valueat = Function.A(pointat.Item1, pointat.Item2, pointat.Item3, timeMesh[i]);
        //                Solutions[i][j] = valueat.Item1 * antinormal.Item1 + valueat.Item2 * antinormal.Item2 + valueat.Item3 * antinormal.Item3;
        //            }
        //            continue;
        //        }
        //        double t0 = timeMesh[i];
        //        double t1 = timeMesh[i - 1];
        //        double t2 = timeMesh[i - 2];
        //        
        //        double deltT = t0 - t2;
        //        double deltT0 = t0 - t1;
        //        double deltT1 = t1 - t2;
        //    
        //        double tau0 = (deltT + deltT0) / (deltT * deltT0);
        //        double tau1 = deltT / (deltT1 * deltT0);
        //        double tau2 = deltT0 / (deltT * deltT1);
        //    
        //        var sparceMatrix = new GlobalMatrix(ribsArr.Count);
        //        Generator.BuildPortait(ref sparceMatrix, ribsArr.Count, elemsArr);

        //        var G = new GlobalMatrix(sparceMatrix);
        //        Generator.FillMatrixG(ref G, ribsArr, elemsArr);
        //
        //        var M = new GlobalMatrix(sparceMatrix);
        //        Generator.FillMatrixM(ref M, ribsArr, elemsArr);
        //
        //        Matrix = G + M + tau0 * M;

        //        var b = new GlobalVector(ribsArr.Count);
        //        Generator.FillVector3D(ref b, ribsArr, elemsArr, timeMesh[i]);
        //        Vector = b - tau2 * M * Solutions[i - 2] + tau1 * M * Solutions[i - 1];

        //        Generator.ConsiderBoundaryConditions(ref Matrix, ref Vector, ribsArr, bordersArr, timeMesh[i]);
        //        (Solutions[i], Discrepancy[i]) = solver.Solve(Matrix, Vector);
        //    }
        //}
    }

    public void WriteData(string path)
    {
        if (Solutions is null) throw new ArgumentNullException("No solutions");

        for (int t = 0; t < timeMesh.Length; t++)
        {
            using var sw = new StreamWriter(path + $"/A_phi/Answer3D/Answer_{timeMesh[t]}.txt");

            for (int i = 0; i < Solutions[t].Size; i++)
                if (i == 38 || i == 45 || i == 67 || i == 74 || i == 41 || i == 42 || i == 70 || i == 71 || i == 52 || i == 53 || i == 56 || i == 57)
                    sw.WriteLine($"{i} {Solutions[t][i]:E8}");

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

        foreach (var elem in elemsArr)
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