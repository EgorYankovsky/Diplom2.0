namespace Project;

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

    public FEM3D(FEM2D fem2d)
    {
        Ax_3D = [];
        Ay_3D = [];
        Az_3D = [];
        Ex_3D = [];
        Ey_3D = [];
        Ez_3D = [];
        Layers = [];
        mesh = new Mesh3Dim
        {
            nodesX = [],
            nodesY = [],
            nodesZ = []
        };
        equationType = fem2d.equationType;
    }

    public void ConstructMesh()
    {

    }

    public void ConstructMesh(FEM2D fem2d)
    {   
        if (fem2d.Mesh2D.nodesR is null) throw new ArgumentNullException("null object");
        if (mesh is null) throw new ArgumentNullException("null object");
        if (mesh.nodesX is null) throw new ArgumentNullException("null object");
        if (mesh.nodesY is null) throw new ArgumentNullException("null object");

        mesh.nodesX = fem2d.Mesh2D.nodesR.Where(r => r < fem2d.Mesh2D.nodesR.Last() / Math.Sqrt(2.0D)).ToList();
        mesh.nodesY = fem2d.Mesh2D.nodesR.Where(r => r < fem2d.Mesh2D.nodesR.Last() / Math.Sqrt(2.0D)).ToList();
        
        mesh.nodesX.Add(fem2d.Mesh2D.nodesR.Last() / Math.Sqrt(2.0D));
        mesh.nodesY.Add(fem2d.Mesh2D.nodesR.Last() / Math.Sqrt(2.0D));
        
        int amount = 2 * mesh.nodesX.Count - 2;

        for (int i = 1; i < amount; i += 2)
        {
            mesh.nodesX.Insert(0, -1.0D * mesh.nodesX[i]);
            mesh.nodesY.Insert(0, -1.0D * mesh.nodesY[i]);
        }

        mesh.mu0 = fem2d.mu0;
        mesh.sigma = fem2d.sigma;

        mesh.nodesZ = fem2d.Mesh2D.nodesZ;
        timeMesh = fem2d.timeMesh;
    }

    public void ConvertResultTo3Dim(FEM2D fem2d)
    {
        Task TaskGenerationAxyz = Task.Run(() => GenerateAxyz(fem2d));
        Task TaskGenerationExyz = Task.Run(() => GenerateExyz(fem2d));
        TaskGenerationAxyz.Wait();
        TaskGenerationExyz.Wait();
    }

    private void GenerateAxyz(FEM2D fem2d)
    {
        if (timeMesh is null) throw new ArgumentNullException();
        if (mesh is null) throw new ArgumentNullException();
        if (mesh.nodesX is null) throw new ArgumentNullException();
        if (mesh.nodesY is null) throw new ArgumentNullException();
        if (mesh.nodesZ is null) throw new ArgumentNullException();
        if (fem2d.pointsArr is null) throw new ArgumentNullException();

        int i = 0;
        
        foreach (var t in timeMesh)
        {
            int j = 0;
            Ax_3D.Add(new GlobalVector(mesh.NodesAmountTotal));
            Ay_3D.Add(new GlobalVector(mesh.NodesAmountTotal));
            Az_3D.Add(new GlobalVector(mesh.NodesAmountTotal));
            foreach (var Z in mesh.nodesZ)
            {
                foreach (var Y in mesh.nodesY)
                {
                    foreach (var X in mesh.nodesX)
                    {
                        var elem = fem2d.GetE_phi(Math.Sqrt(X * X + Y * Y), Z, t);

                        Ax_3D[i][j] = -1.0D * (Y / Math.Sqrt(X * X + Y * Y)) * 
                            BasisFunctions2D.GetValue(
                                fem2d.A_phi[i][elem[0]], fem2d.A_phi[i][elem[1]],
                                fem2d.A_phi[i][elem[2]], fem2d.A_phi[i][elem[3]],
                                fem2d.pointsArr[elem[0]].R, fem2d.pointsArr[elem[1]].R,
                                fem2d.pointsArr[elem[0]].Z, fem2d.pointsArr[elem[3]].Z, 
                                Math.Sqrt(X * X + Y * Y), Z);
                        
                        Ay_3D[i][j] = X / Math.Sqrt(X * X + Y * Y) * 
                            BasisFunctions2D.GetValue(
                                fem2d.A_phi[i][elem[0]], fem2d.A_phi[i][elem[1]],
                                fem2d.A_phi[i][elem[2]], fem2d.A_phi[i][elem[3]],
                                fem2d.pointsArr[elem[0]].R, fem2d.pointsArr[elem[1]].R,
                                fem2d.pointsArr[elem[0]].Z, fem2d.pointsArr[elem[3]].Z, 
                                Math.Sqrt(X * X + Y * Y), Z);
                        
                        j++;
                    }
                }
            }
            i++;
        }
    }

    private void GenerateExyz(FEM2D fem2d)
    {
        if (timeMesh is null) throw new ArgumentNullException();
        if (mesh is null) throw new ArgumentNullException();
        if (mesh.nodesX is null) throw new ArgumentNullException();
        if (mesh.nodesY is null) throw new ArgumentNullException();
        if (mesh.nodesZ is null) throw new ArgumentNullException();
        if (fem2d.pointsArr is null) throw new ArgumentNullException();

        int i = 0;
        
        foreach (var t in timeMesh)
        {
            if (t == timeMesh.Last())
                continue;
            int j = 0;
            Ex_3D.Add(new GlobalVector(mesh.NodesAmountTotal));
            Ey_3D.Add(new GlobalVector(mesh.NodesAmountTotal));
            Ez_3D.Add(new GlobalVector(mesh.NodesAmountTotal));
            foreach (var Z in mesh.nodesZ)
            {
                foreach (var Y in mesh.nodesY)
                {
                    foreach (var X in mesh.nodesX)
                    {
                        var elem = fem2d.GetE_phi(Math.Sqrt(X * X + Y * Y), Z, t);

                        Ex_3D[i][j] = -1.0D * (Y / Math.Sqrt(X * X + Y * Y)) * 
                            BasisFunctions2D.GetValue(
                                fem2d.E_phi2D[i][elem[0]], fem2d.E_phi2D[i][elem[1]],
                                fem2d.E_phi2D[i][elem[2]], fem2d.E_phi2D[i][elem[3]],
                                fem2d.pointsArr[elem[0]].R, fem2d.pointsArr[elem[1]].R,
                                fem2d.pointsArr[elem[0]].Z, fem2d.pointsArr[elem[3]].Z, 
                                Math.Sqrt(X * X + Y * Y), Z);
                        
                        Ey_3D[i][j] = X / Math.Sqrt(X * X + Y * Y) * 
                            BasisFunctions2D.GetValue(
                                fem2d.E_phi2D[i][elem[0]], fem2d.E_phi2D[i][elem[1]],
                                fem2d.E_phi2D[i][elem[2]], fem2d.E_phi2D[i][elem[3]],
                                fem2d.pointsArr[elem[0]].R, fem2d.pointsArr[elem[1]].R,
                                fem2d.pointsArr[elem[0]].Z, fem2d.pointsArr[elem[3]].Z, 
                                Math.Sqrt(X * X + Y * Y), Z);
                        
                        j++;
                    }
                }
            }
            i++;
        }
    }

    public void GenerateArrays()
    {
        if (mesh is null) throw new ArgumentNullException("mesh is null!");
        pointsArr = MeshGenerator.GenerateListOfPoints(mesh);
        ribsArr = MeshGenerator.GenerateListOfRibs(mesh, pointsArr);
        elemsArr = MeshGenerator.GenerateListOfElems(mesh);
        //MeshGenerator.SelectRibs(ref ribsArr, ref elemsArr);
        // ! А надо ли?
        //bordersArr = MeshGenerator.GenerateListOfBorders(mesh);
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
        var sparceMatrix = new GlobalMatrix(ribsArr.Count);
        Generator.BuildPortait(ref sparceMatrix, ribsArr.Count, elemsArr);

        var G = new GlobalMatrix(sparceMatrix);

        Generator.FillMatrixG(ref G, ribsArr, elemsArr);
        //Generator.ConsiderBoundaryConditions(ref m, arrBd);
        

        var M = new GlobalMatrix(sparceMatrix);
        Generator.FillMatrixM(ref M, ribsArr, elemsArr);
    }
}