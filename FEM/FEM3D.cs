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
    public List<GlobalVector> Ex_3D;

    public List<GlobalVector> Ey_3D;

    public List<GlobalVector> Ez_3D;

    public FEM3D(FEM2D fem2d)
    {
        Ex_3D = new();
        Ey_3D = new();
        Ez_3D = new();
        mesh = new Mesh3Dim
        {
            nodesX = new(),
            nodesY = new(),
            nodesZ = new()
        };
        equationType = fem2d.equationType;
    }

    public void ConstructMesh(FEM2D fem2d)
    {   
        if (fem2d.Mesh2D.nodesR is null) throw new ArgumentNullException("null object");
        if (mesh is null) throw new ArgumentNullException("null object");
        if (mesh.nodesX is null) throw new ArgumentNullException("null object");
        if (mesh.nodesY is null) throw new ArgumentNullException("null object");

        foreach (var R in fem2d.Mesh2D.nodesR)
        {
            mesh.nodesX.Add(R / Math.Sqrt(2));
            mesh.nodesY.Add(R / Math.Sqrt(2));
        }
        mesh.nodesZ = fem2d.Mesh2D.nodesZ;
        timeMesh = fem2d.timeMesh;
    }

    public void GenerateExy(FEM2D fem2d)
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
            Ex_3D.Add(new GlobalVector(mesh.NodesAmountTotal));
            Ey_3D.Add(new GlobalVector(mesh.NodesAmountTotal));
            Ez_3D.Add(new GlobalVector(mesh.NodesAmountTotal));
            foreach (var Z in mesh.nodesZ)
            {
                foreach (var Y in mesh.nodesY)
                {
                    foreach (var X in mesh.nodesX)
                    {
                        Console.WriteLine($"{X:E10} {Y:E10} {Z:E10}");
                        var elem = fem2d.GetE_phi(Math.Sqrt(X * X + Y * Y), Z, t);
                        Ex_3D[i][j] = -1.0D * (Y / Math.Sqrt(X * X + Y * Y)) * 
                            BasisFunctions2D.GetValue(elem[0], elem[1], elem[2], elem[3],
                                fem2d.pointsArr[elem[0]].R, fem2d.pointsArr[elem[1]].R,
                                fem2d.pointsArr[elem[0]].Z, fem2d.pointsArr[elem[3]].Z, 
                                Math.Sqrt(X * X + Y * Y), Z);
                        
                        Ey_3D[i][j] = X / Math.Sqrt(X * X + Y * Y) * 
                            BasisFunctions2D.GetValue(elem[0], elem[1], elem[2], elem[3],
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
}