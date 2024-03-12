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
    public GlobalVector[] Ex_phi3D;

    public GlobalVector[] Ey_phi3D;

    public FEM3D(FEM2D fem2d)
    {
        equationType = fem2d.equationType;
    }

    private void GenerateExy(FEM2D fem2d)
    {
       int i = 0;
        foreach (var t in _timeMesh)
        {
            int j = 0;
            Ex_phi3D.Append(new GlobalVector(_mesh.NodesAmountTotal));
            foreach (var Z in _mesh.nodesZ)
            {
                foreach (var Y in _mesh.nodesY)
                {
                    foreach (var X in _mesh.nodesX)
                    {
                        var elem = fem2d.GetElem(Math.Sqrt(X * X + Y * Y), Z, t);
                        //Ex_phi3D[i][j] = BasisFunctions2D.GetValue(elem[0], elem[1], elem[2], elem[3],
                                                                 
                        //                                            Math.Sqrt(X * X + Y * Y), Z);

                    }
                }
            }
            i++;
        }
    }
}