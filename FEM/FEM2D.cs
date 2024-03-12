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

public class FEM2D : FEM
{
    public FEM2D()
    {
        _mesh2D = new();
        _mu0 = new List<double>();
        _sigma = new List<double>();
    }

    public Mesh2Dim _mesh2D; 

    public GlobalVector[] A_phi;

    public GlobalVector[] E_phi2D;

    public void ReadData(string infoPath, string bordersPath, string timePath)
    {
        MeshReader.ReadMesh(infoPath, bordersPath, ref _mesh2D);
        using var sr = new StreamReader(timePath);
        SetTimeMesh(sr.ReadLine() ?? "0 0 1");        
        Debug.WriteLine("All data read correctly");
    }

    public void ConstructMesh()
    {
        MeshGenerator.GenerateMesh(ref _mesh2D);
        MeshGenerator.OutputPoints(ref _mesh2D);
        _pointsArr = new();
        MeshGenerator.GenerateListOfElems(ref _mesh2D, _pointsArr);
        MeshGenerator.GenerateListOfBorders(ref _mesh2D);
        Debug.WriteLine("Mesh built correctly");
    }

    public void SubmitGeneratedData()
    {
        _elemsArr = new();
        _bordersArr = new();
        _mu0 = _mesh2D.mu0;
        _sigma = _mesh2D.sigma;
        Debug.WriteLine("Generated data submited");
    }

    public void SetSolver(ISolver solver)
    {
        this.solver = solver;
        Debug.WriteLine("Solvet set");
    }

    public void Solve()
    {
        switch(equationType)
        {
            case EquationType.Elliptic:
                _matrix = new GlobalMatrix(_pointsArr.Length);
                Generator.BuildPortait(ref _matrix, _pointsArr, _elemsArr);
                Generator.FillMatrix(ref _matrix, _pointsArr, _elemsArr, _bordersArr);
                _vector = new GlobalVector(_pointsArr, _bordersArr, 1.0);
                _solutions[0] = solver.Solve(_matrix, _vector);
            break;
            case EquationType.Parabolic:
                if (_timeMesh is null) throw new ArgumentNullException("_timeMesh is null!");

                for (int i = 0; i < _timeMesh.Length; i++)
                {
                    Debug.WriteLine($"\nTime layer: {_timeMesh[i]}");
                    Thread.Sleep(1500);
                    if (i == 0)
                    {
                        _matrix = new GlobalMatrix(_pointsArr.Length);
                        Generator.BuildPortait(ref _matrix, _pointsArr, _elemsArr);
                        Generator.FillMatrix(ref _matrix, _pointsArr, _elemsArr, _bordersArr);
                        _vector = new GlobalVector(_pointsArr, _bordersArr, _timeMesh[i]);
                        _solutions[0] = solver.Solve(_matrix, _vector);
                    }
                    else if (i == 1)
                    {
                        double deltT = _timeMesh[i] - _timeMesh[i - 1];
                        var matrix1 = new GlobalMatrix(_pointsArr.Length);
                        Generator.BuildPortait(ref matrix1, _pointsArr, _elemsArr);
                        Generator.FillMatrix(ref matrix1, _pointsArr, _elemsArr, _bordersArr);

                        var M = new GlobalMatrix(_pointsArr.Length); // ???
                        Generator.BuildPortait(ref M, _pointsArr, _elemsArr);
                        Generator.FillMatrix(ref M, _pointsArr, _elemsArr, _bordersArr, 1.0);
                        _matrix = 1.0 / deltT * M + matrix1;
                    
                        var bi = new GlobalVector(_pointsArr, _bordersArr, _timeMesh[i]);
                        _vector = bi + 1.0 / deltT * M * _solutions[i - 1];
                        _solutions[i] = solver.Solve(_matrix, _vector);
                    }
                    else
                    {
                        double deltT = _timeMesh[i] - _timeMesh[i - 2];
                        double deltT1 = _timeMesh[i - 1] - _timeMesh[i - 2];
                        double deltT0 = _timeMesh[i] - _timeMesh[i - 1];

                        var matrix1 = new GlobalMatrix(_pointsArr.Length);
                        Generator.BuildPortait(ref matrix1, _pointsArr, _elemsArr);
                        Generator.FillMatrix(ref matrix1, _pointsArr, _elemsArr, _bordersArr);
                        
                        var M = new GlobalMatrix(_pointsArr.Length); // ???
                        Generator.BuildPortait(ref M, _pointsArr, _elemsArr);
                        Generator.FillMatrix(ref M, _pointsArr, _elemsArr, _bordersArr, 1.0);
                        
                        
                        _matrix = (deltT + deltT0) / (deltT * deltT0) * M + matrix1;

                        var bi = new GlobalVector(_pointsArr, _bordersArr, _timeMesh[i]);
                        _vector = bi - deltT0 / (deltT * deltT1) * M * _solutions[i - 2] + deltT / (deltT1 * deltT0) * M * _solutions[i - 1];
                        _solutions[i] = solver.Solve(_matrix, _vector);
                    }
                }
            break;
        }
        A_phi = _solutions;
        Debug.WriteLine("Lin eq solved");
    }

    public void WriteData()
    {
        if (_answer is null)
            throw new Exception("Vector _answer is null");
        for (int i = 0; i < _answer.Size; i++)
            Console.WriteLine($"{_answer[i]:E15}");
    }

    public void WriteData(string _path)
    {
        if (_timeMesh is not null)
            for (int i = 0; i < _timeMesh.Length; i++)
            {
                using var sw = new StreamWriter($"{_path}Answer_time={_timeMesh[i]}.dat");
                for (int j = 0;   j < _solutions[i].Size; j++)
                    sw.WriteLine($"{_solutions[i][j]:E8}");
                sw.Close();
            }
        else
        {
            using var sw = new StreamWriter($"{_path}Answer.dat");
            for (int j = 0; j < _solutions[0].Size; j++)
                sw.WriteLine($"{_solutions[0][j]:E8}");
            sw.Close();
        }
        WriteDiscrepancy(_path);

    }

    private void WriteDiscrepancy(string _path)
    {
        using var sw = new StreamWriter($"{_path}Discrepancy.dat");
        
        double avgDisc = 0.0;
        double maxDisc = 0.0;
        int NotNaNamount = 0;
        
        double sumU = 0.0D;
        double sumD = 0.0D;

        List<double> TheorAnswer = new();
        foreach (var Z in _mesh2D.nodesZ)
        {
            foreach (var R in _mesh2D.nodesR)
            {
                TheorAnswer.Add(Function.U(R, Z));
            }
        }
        
        for (int i = 0; i < _solutions[0].Size; i++)
        {
            double absDiff = _solutions[0][i] - TheorAnswer[i];
            double currDisc = Math.Abs(_solutions[0][i] - TheorAnswer[i]) / Math.Abs(TheorAnswer[i]);
            if (Math.Abs(maxDisc) < Math.Abs(currDisc))
                maxDisc = currDisc;
            if (!double.IsNaN(currDisc))
            {
                avgDisc += currDisc;
                NotNaNamount++;
            }
            sumU += absDiff * absDiff;
            sumD += TheorAnswer[i] * TheorAnswer[i];

            sw.WriteLine($"{absDiff:E15} {currDisc:E15}");
        }

        avgDisc = Math.Sqrt(sumU) / Math.Sqrt(sumD);
        sw.WriteLine($"Средняя невязка: {avgDisc:E15}");
        sw.WriteLine($"Максимальная невязка: {maxDisc:E15}");
        sw.WriteLine($"С: {avgDisc:E7}");
        sw.WriteLine($"М: {maxDisc:E7}");
        sw.Close();
    }

    public void GenerateVectorEphi()
    {
        if (_timeMesh is null) throw new NullReferenceException("_timeMesh is null");
        E_phi2D = new GlobalVector[A_phi.Length - 1];
        for (int i = 0; i < E_phi2D.Length; i++)
            E_phi2D[i] = (1.0 / (_timeMesh[i + 1] - _timeMesh[i])) * (A_phi[i + 1] - A_phi[i]);
    }

    internal List<int> GetElem(double r, double z, double t)
    {
        int i;
        for (i = 0; i < _mesh2D.nodesR.Count - 1; i++)
            if (_mesh2D.nodesR[i] <= r && r <= _mesh2D.nodesR[i + 1])
                break;
        int j;
        for (j = 0; j < _mesh2D.nodesZ.Count - 1; j++)
            if (_mesh2D.nodesZ[j] <= z && z <= _mesh2D.nodesZ[j + 1])
                break;
        return _elemsArr[j * (_mesh2D.nodesR.Count - 1) + i];
    }
}