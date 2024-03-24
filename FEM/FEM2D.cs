using MathObjects;
using Solver;
using Grid;
using System.Diagnostics;
using Functions;

namespace Project;

public class FEM2D : FEM
{
    public FEM2D()
    {
        Mesh2D = new();
        mu0 = new List<double>();
        sigma = new List<double>();
    }

    public Mesh2Dim Mesh2D; 

    public GlobalVector[] A_phi;

    public GlobalVector[] E_phi2D;

    public void ReadData(string infoPath, string bordersPath, string timePath)
    {
        MeshReader.ReadMesh(infoPath, bordersPath, ref Mesh2D);
        using var sr = new StreamReader(timePath);
        SetTimeMesh(sr.ReadLine() ?? "0 0 1");        
        Debug.WriteLine("All data read correctly");
    }

    public void ConstructMesh()
    {
        MeshGenerator.GenerateMesh(ref Mesh2D);
        MeshGenerator.OutputPoints(ref Mesh2D);
        pointsArr = new();
        MeshGenerator.GenerateListOfElems(ref Mesh2D, pointsArr);
        MeshGenerator.GenerateListOfBorders(ref Mesh2D);
        Debug.WriteLine("Mesh built correctly");
    }

    public void SubmitGeneratedData()
    {
        elemsArr = new();
        bordersArr = new();
        mu0 = Mesh2D.mu0;
        sigma = Mesh2D.sigma;
        Debug.WriteLine("Generated data submited");
    }

    public void SetSolver(ISolver solver)
    {
        this.solver = solver;
        Debug.WriteLine("Solvet set");
    }

    public void Solve()
    {
        if (pointsArr is null) throw new ArgumentNullException("points array is null !");
        if (elemsArr is null) throw new ArgumentNullException("elems array is null !");
        if (bordersArr is null) throw new ArgumentNullException("borders array is null !");
        if (Solutions is null) throw new ArgumentNullException("solutions array is null !");
        if (solver is null) throw new ArgumentNullException("solver is null !");

        switch(equationType)
        {
            case EquationType.Elliptic:
                Matrix = new GlobalMatrix(pointsArr.Length);
                Generator.BuildPortait(ref Matrix, pointsArr.Length, elemsArr);
                Generator.FillMatrix(ref Matrix, pointsArr, elemsArr, bordersArr, TypeOfMatrixM.Mrr);
                Vector = new GlobalVector(pointsArr, bordersArr, 1.0);
                Solutions[0] = solver.Solve(Matrix, Vector);
            break;

            case EquationType.Parabolic:
                if (timeMesh is null) throw new ArgumentNullException("timeMesh is null!");

                for (int i = 0; i < timeMesh.Length; i++)
                {
                    Debug.WriteLine($"\nTime layer: {timeMesh[i]}");
                    Thread.Sleep(1500);
                    if (i == 0)
                    {
                        Matrix = new GlobalMatrix(pointsArr.Length);
                        Generator.BuildPortait(ref Matrix, pointsArr.Length, elemsArr);
                        Generator.FillMatrix(ref Matrix, pointsArr, elemsArr, bordersArr, TypeOfMatrixM.Mrr);
                        Vector = new GlobalVector(pointsArr, bordersArr, timeMesh[i]);
                        Solutions[0] = solver.Solve(Matrix, Vector);
                    }
                    else if (i == 1)
                    {
                        /*
                        double deltT = timeMesh[i] - timeMesh[i - 1];
                        var matrix1 = new GlobalMatrix(pointsArr.Length);
                        Generator.BuildPortait(ref matrix1, pointsArr, elemsArr);
                        Generator.FillMatrix(ref matrix1, pointsArr, elemsArr, bordersArr);

                        var M = new GlobalMatrix(pointsArr.Length); // ???
                        Generator.BuildPortait(ref M, pointsArr, elemsArr);
                        Generator.FillMatrix(ref M, pointsArr, elemsArr, bordersArr, 1.0);
                        Matrix = 1.0 / deltT * M + matrix1;
                    
                        var bi = new GlobalVector(pointsArr, bordersArr, timeMesh[i]);
                        Vector = bi + 1.0 / deltT * M * Solutions[i - 1];
                        */
                        Solutions[i] = Solutions[i - 1];
                    }
                    else
                    {
                        double deltT = timeMesh[i] - timeMesh[i - 2];
                        double deltT0 = timeMesh[i] - timeMesh[i - 1];
                        double deltT1 = timeMesh[i - 1] - timeMesh[i - 2];

                        double tau0 = (deltT + deltT0) / (deltT * deltT0);
                        double tau1 = deltT / (deltT1 * deltT0);
                        double tau2 = deltT0 / (deltT * deltT0);


                        var matrix1 = new GlobalMatrix(pointsArr.Length);
                        Generator.BuildPortait(ref matrix1, pointsArr.Length, elemsArr);
                        Generator.FillMatrix(ref matrix1, pointsArr, elemsArr, bordersArr, TypeOfMatrixM.Mrr);
                        
                        var M = new GlobalMatrix(pointsArr.Length); // ???
                        Generator.BuildPortait(ref M, pointsArr.Length, elemsArr);
                        Generator.FillMatrix(ref M, pointsArr, elemsArr, bordersArr, TypeOfMatrixM.Mr);
                        
                        
                        Matrix = (deltT + deltT0) / (deltT * deltT0) * M + matrix1;

                        var bi = new GlobalVector(pointsArr, bordersArr, timeMesh[i]);
                        Vector = bi - deltT0 / (deltT * deltT1) * M * Solutions[i - 2] + deltT / (deltT1 * deltT0) * M * Solutions[i - 1];
                        Solutions[i] = solver.Solve(Matrix, Vector);
                    }
                }
            break;
        }
        A_phi = Solutions;
        Debug.WriteLine("Lin eq solved");
    }

    public void WriteData()
    {
        if (Answer is null)
            throw new Exception("Vector _answer is null");
        for (int i = 0; i < Answer.Size; i++)
            Console.WriteLine($"{Answer[i]:E15}");
    }

    public void WriteData(string _path)
    {
        if (Solutions is null) throw new ArgumentNullException();

        if (timeMesh is not null)
            for (int i = 0; i < timeMesh.Length; i++)
            {
                using var sw = new StreamWriter($"{_path}Answer_time={timeMesh[i]}.dat");
                for (int j = 0;   j < Solutions[i].Size; j++)
                    sw.WriteLine($"{Solutions[i][j]:E8}");
                sw.Close();
            }
        else
        {
            using var sw = new StreamWriter($"{_path}Answer.dat");
            for (int j = 0; j < Solutions[0].Size; j++)
                sw.WriteLine($"{Solutions[0][j]:E8}");
            sw.Close();
        }
        WriteDiscrepancy(_path);
    }

    private void WriteDiscrepancy(string _path)
    {
        if (Mesh2D is null) throw new ArgumentNullException();
        if (Mesh2D.nodesR is null) throw new ArgumentNullException();
        if (Mesh2D.nodesZ is null) throw new ArgumentNullException();
        if (Solutions is null) throw new ArgumentNullException();


        using var sw = new StreamWriter($"{_path}Discrepancy.dat");

        double avgDisc = 0.0;
        double maxDisc = 0.0;
        int NotNaNamount = 0;
        
        double sumU = 0.0D;
        double sumD = 0.0D;

        List<double> TheorAnswer = new();
        foreach (var Z in Mesh2D.nodesZ)
        {
            foreach (var R in Mesh2D.nodesR)
            {
                TheorAnswer.Add(Function.U(R, Z));
            }
        }
        
        for (int i = 0; i < Solutions[0].Size; i++)
        {
            double absDiff = Solutions[0][i] - TheorAnswer[i];
            double currDisc = Math.Abs(Solutions[0][i] - TheorAnswer[i]) / Math.Abs(TheorAnswer[i]);
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
        if (timeMesh is null) throw new NullReferenceException("timeMesh is null");
        E_phi2D = new GlobalVector[A_phi.Length - 1];
        for (int i = 0; i < E_phi2D.Length; i++)
            E_phi2D[i] = (1.0 / (timeMesh[i + 1] - timeMesh[i])) * (A_phi[i] - A_phi[i + 1]);
    }

    internal List<int> GetE_phi(double r, double z, double t)
    {
        if (Mesh2D is null) throw new ArgumentNullException();
        if (Mesh2D.nodesR is null) throw new ArgumentNullException();
        if (Mesh2D.nodesZ is null) throw new ArgumentNullException();
        if (elemsArr is null) throw new ArgumentNullException();

        int i;
        for (i = 0; i < Mesh2D.nodesR.Count - 1; i++)
            if (Mesh2D.nodesR[i] <= r && r <= Mesh2D.nodesR[i + 1])
                break;
        int j;
        for (j = 0; j < Mesh2D.nodesZ.Count - 1; j++)
            if (Mesh2D.nodesZ[j] <= z && z <= Mesh2D.nodesZ[j + 1])
                break;

        return elemsArr[j * (Mesh2D.nodesR.Count - 1) + i];
    }
}