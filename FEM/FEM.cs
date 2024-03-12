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


public class FEM
{
    private protected double[]? _timeMesh;

    protected internal EquationType equationType;

    private protected Mesh? _mesh;

    protected ISolver? solver;

    protected ArrayOfElems? _elemsArr; 

    protected ArrayOfPoints? _pointsArr; 

    protected ArrayOfBorders? _bordersArr; 

    protected List<double>? _mu0;

    protected List<double>? _sigma;

    public GlobalMatrix? _matrix;

    public GlobalVector? _vector;

    public GlobalVector? _answer;

    public GlobalVector[]? _solutions;

    private protected void SetTimeMesh(string data)
    {
        var info = data.Split(" ");
        var t0 = double.Parse(info[0]);
        var t1 = double.Parse(info[1]);
        var tn = int.Parse(info[2]);
        var tk = double.Parse(info[3]);

        if (t0 - t1 == 0)
        {
            equationType = EquationType.Elliptic;
            _solutions = new GlobalVector[1];
            return;
        }

        _timeMesh = new double[tn + 1];
        _solutions = new GlobalVector[tn + 1];

        double h = t1 - t0;
        double denominator = 0.0;

        for (int j = 0; j < tn; j++)
            denominator += Math.Pow(double.Parse(info[3]), j);

        double x0 = h / denominator;
        _timeMesh[0] = t0;
        for(int j = 0; j < tn - 1; j++)
            _timeMesh[j + 1] = _timeMesh[j] + x0 * Math.Pow(tk, j);
        _timeMesh[^1] = t1;
        equationType = EquationType.Parabolic;
        Debug.WriteLine("Time mesh built correctly!");
    }
}