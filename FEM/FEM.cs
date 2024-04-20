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


public enum EquationType
{
    Elliptic,
    
    Parabolic
}


public abstract class FEM(TimeMesh time, int dimesion)
{
    private static string _elemspath2D = Path.GetFullPath("../../../../Data/Subtotals/2_dim/Elems.poly");
    
    private static string _pointspath2D = Path.GetFullPath("../../../../Data/Subtotals/2_dim/Points.poly");
    
    private static string _borderspath2D = Path.GetFullPath("../../../../Data/Subtotals/2_dim/Borders.poly");

    private static string _elemspath3D = Path.GetFullPath("../../../../Data/Subtotals/3_dim/Elems.poly");
    
    private static string _pointspath3D = Path.GetFullPath("../../../../Data/Subtotals/3_dim/Points.poly");
    
    private static string _borderspath3D = Path.GetFullPath("../../../../Data/Subtotals/3_dim/Borders.poly");

    public TimeMesh Time = time;

    protected internal EquationType equationType;

    private protected Mesh? mesh;

    protected ISolver? solver;

    public ArrayOfElems elemsArr = dimesion == 2 ? new(_elemspath2D) : new(_elemspath3D); 

    public ArrayOfPoints pointsArr = dimesion == 2 ? new(_pointspath2D) : new(_pointspath3D);

    public ArrayOfBorders bordersArr = dimesion == 2 ? new(_borderspath2D) : new(_borderspath3D);

    public GlobalMatrix? Matrix;

    public GlobalVector? Vector;

    public GlobalVector? Answer;

    public GlobalVector[] Solutions = new GlobalVector[time.Count];

    public GlobalVector[] Discrepancy = new GlobalVector[time.Count];

    public void SetSolver(ISolver solver)
    {
        this.solver = solver;
        Debug.WriteLine("Solvet set");
    }
}