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


public abstract class FEM
{
    internal double[]? timeMesh;

    protected internal EquationType equationType;

    private protected Mesh? mesh;

    protected ISolver? solver;

    public ArrayOfElems? elemsArr; 

    public ArrayOfPoints? pointsArr; 

    public ArrayOfBorders? bordersArr; 

    protected List<double>? mu0;

    protected List<double>? sigma;

    public GlobalMatrix? Matrix;

    public GlobalVector? Vector;

    public GlobalVector? Answer;

    public GlobalVector[]? Solutions;

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
            Solutions = new GlobalVector[1];
            return;
        }

        timeMesh = new double[tn + 1];
        Solutions = new GlobalVector[tn + 1];

        double h = t1 - t0;
        double denominator = 0.0;

        for (int j = 0; j < tn; j++)
            denominator += Math.Pow(double.Parse(info[3]), j);

        double x0 = h / denominator;
        timeMesh[0] = t0;
        for(int j = 0; j < tn - 1; j++)
            timeMesh[j + 1] = timeMesh[j] + x0 * Math.Pow(tk, j);
        timeMesh[^1] = t1;
        equationType = EquationType.Parabolic;
        Debug.WriteLine("Time mesh built correctly!");
    }
}