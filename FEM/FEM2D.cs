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
        if (Discrepancy is null) throw new ArgumentNullException("discrepancy array is null !");
        if (solver is null) throw new ArgumentNullException("solver is null !");

        switch(equationType)
        {
            case EquationType.Elliptic:
                Matrix = new GlobalMatrix(pointsArr.Length);
                Generator.BuildPortait(ref Matrix, pointsArr.Length, elemsArr);
                Generator.FillMatrix(ref Matrix, pointsArr, elemsArr, bordersArr, TypeOfMatrixM.Mrr);
                
                Vector = new GlobalVector(pointsArr);
                
                Generator.FillVector(ref Vector, pointsArr, elemsArr, bordersArr, 1.0);
            

                (Solutions[0], Discrepancy[0]) = solver.Solve(Matrix, Vector);
            break;

            case EquationType.Parabolic:
                if (timeMesh is null) throw new ArgumentNullException("timeMesh is null!");

                for (int i = 0; i < timeMesh.Length; i++)
                {
                    Debug.WriteLine($"\nTime layer: {timeMesh[i]}");
                    Thread.Sleep(1500);
                    if (i == 0)
                    {
#if RELEASE
                        Matrix = new GlobalMatrix(pointsArr.Length);
                        Generator.BuildPortait(ref Matrix, pointsArr.Length, elemsArr);
                        Generator.FillMatrix(ref Matrix, pointsArr, elemsArr, bordersArr, TypeOfMatrixM.Mrr);
                        Vector = new GlobalVector(pointsArr);
                        (Solutions[0], Discrepancy[0]) = solver.Solve(Matrix, Vector);
#endif
#if DEBUG
                        Solutions[0] = new GlobalVector(pointsArr.Length);
                        for (int j = 0; j < Solutions[0].Size; j++)
                        {
                            Solutions[0][j] = Function.U(0.0, 0.0, timeMesh[i]);
                        }
#endif
                    }
                    else if (i == 1)
                    {
#if RELEASE
                        Solutions[i] = Solutions[i - 1];
#endif
#if DEBUG
                        Solutions[1] = new GlobalVector(pointsArr.Length);
                        for (int j = 0; j < Solutions[1].Size; j++)
                        {
                            Solutions[1][j] = Function.U(0.0, 0.0, timeMesh[i]);
                        }
#endif
                    }
                    else
                    {
                        double deltT = timeMesh[i] - timeMesh[i - 2];
                        double deltT0 = timeMesh[i] - timeMesh[i - 1];
                        double deltT1 = timeMesh[i - 1] - timeMesh[i - 2];

                        double tau0 = (deltT + deltT0) / (deltT * deltT0);
                        double tau1 = deltT / (deltT1 * deltT0);
                        double tau2 = deltT0 / (deltT * deltT1);

                        var matrix1 = new GlobalMatrix(pointsArr.Length);
                        Generator.BuildPortait(ref matrix1, pointsArr.Length, elemsArr);
                        Generator.FillMatrix(ref matrix1, pointsArr, elemsArr, bordersArr, TypeOfMatrixM.Mrr);

                        var M = new GlobalMatrix(pointsArr.Length); // ???
                        Generator.BuildPortait(ref M, pointsArr.Length, elemsArr);
                        Generator.FillMatrix(ref M, pointsArr, elemsArr, bordersArr, TypeOfMatrixM.Mr);

                        Matrix = (tau0 * M) + matrix1;
                        Generator.ConsiderBoundaryConditions(ref Matrix, bordersArr);

                        var bi = new GlobalVector(pointsArr);
                        Generator.FillVector(ref bi, pointsArr, elemsArr, bordersArr, timeMesh[i]);

                        Vector = bi - (tau2 * M * Solutions[i - 2]) + (tau1 * M * Solutions[i - 1]);
                        
                        Generator.ConsiderBoundaryConditions(ref Vector, pointsArr, bordersArr, timeMesh[i]);
                        
                        (Solutions[i], Discrepancy[i]) = solver.Solve(Matrix, Vector);
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
        if (A_phi is null) throw new ArgumentNullException();
#if RELEASE
        if (E_phi2D is null) throw new ArgumentNullException();
#endif
        if (timeMesh is not null)
        {
            for (int i = 0; i < timeMesh.Length; i++)
            {
                using var sw = new StreamWriter($"{_path}\\A_phi\\Answer\\Answer_Aphi_time={timeMesh[i]}.dat");
                for (int j = 0;   j < A_phi[i].Size; j++)
                    sw.WriteLine($"{A_phi[i][j]:E8}");
                sw.Close();
            }
#if RELEASE
            for (int i = 0; i < timeMesh.Length; i++)
            {
                using var sw = new StreamWriter($"{_path}\\E_phi\\Answer\\Answer_Ephi_time={timeMesh[i]}.dat");
                for (int j = 0;   j < E_phi2D[i].Size; j++)
                    sw.WriteLine($"{E_phi2D[i][j]:E8}");
                sw.Close();
            }
#endif
        }
        else
        {
            using var sw = new StreamWriter($"{_path}\\A_phi\\Answer\\Answer.dat");
            for (int j = 0; j < A_phi[0].Size; j++)
                sw.WriteLine($"{A_phi[0][j]:E8}");
            sw.Close();
#if RELEASE
            using var sw1 = new StreamWriter($"{_path}\\E_phi\\Answer\\Answer.dat");
            for (int j = 0; j < E_phi2D[0].Size; j++)
                sw1.WriteLine($"{E_phi2D[0][j]:E8}");
            sw1.Close();
#endif
        }
    }

    public void WriteDiscrepancy(string _path)
    {
        if (Mesh2D is null) throw new ArgumentNullException();
        if (Mesh2D.nodesR is null) throw new ArgumentNullException();
        if (Mesh2D.nodesZ is null) throw new ArgumentNullException();
        if (Solutions is null) throw new ArgumentNullException();
        if (Discrepancy is null) throw new ArgumentNullException();
        if (A_phi is null) throw new ArgumentNullException();


        if (timeMesh is not null)
        {
            for (int i = 0; i < timeMesh.Length; i++)
            {
                if (i == 1)
                    continue;
                using var sw_d = new StreamWriter($"{_path}\\A_phi\\Discrepancy\\Discrepancy_Aphi_time={timeMesh[i]}.dat");

                int NotNaNamount = 0;
                double maxDisc = 0.0;
                double avgDisc = 0.0;

                double sumU = 0.0D;
                double sumD = 0.0D;

                List<double> TheorAnswer = new();
                foreach (var Z in Mesh2D.nodesZ)
                {
                    foreach (var R in Mesh2D.nodesR)
                    {
                        TheorAnswer.Add(Function.U(R, Z, timeMesh[i]));
                    }
                }

                for (int j = 0; j < A_phi[i].Size; j++)
                {
                    double absDiff = Math.Abs(A_phi[i][j] - TheorAnswer[j]);
                    double currDisc = Math.Abs(A_phi[i][j] - TheorAnswer[j]) / Math.Abs(TheorAnswer[j]);
                    
                    if (Math.Abs(maxDisc) < Math.Abs(currDisc))
                        maxDisc = currDisc;

                    if (!double.IsNaN(currDisc) || absDiff > 1E-15)
                    {
                        avgDisc += currDisc;
                        NotNaNamount++;                        
                        sumU += absDiff * absDiff;
                        sumD += A_phi[i][j] * A_phi[i][j];

                    }

                    //sumU += absDiff * absDiff;
                    //sumD += A_phi[i][j] * A_phi[i][j];

                    sw_d.WriteLine($"{absDiff:E8} {currDisc:E8}");
                }
                avgDisc = Math.Sqrt(sumU) / Math.Sqrt(sumD);
                sw_d.WriteLine($"Средняя невязка: {avgDisc:E15}");
                sw_d.WriteLine($"Максимальная невязка: {maxDisc:E15}");
                sw_d.WriteLine($"С: {avgDisc:E7}");
                sw_d.WriteLine($"М: {maxDisc:E7}");
                sw_d.Close();
            }
        }
        else
        {
            using var sw_d = new StreamWriter($"{_path}\\A_phi\\Discrepancy\\Discrepancy_Aphi.dat");

            int NotNaNamount = 0;
            double maxDisc = 0.0;
            double avgDisc = 0.0;

            double sumU = 0.0D;
            double sumD = 0.0D;
            
            List<double> TheorAnswer = new();
            foreach (var Z in Mesh2D.nodesZ)
            {
                foreach (var R in Mesh2D.nodesR)
                {
                    TheorAnswer.Add(Function.U(R, Z, 0.0D));
                }
            }
            
            for (int j = 0; j < A_phi[0].Size; j++)
            {
                double absDiff = Math.Abs(A_phi[0][j] - TheorAnswer[j]);
                double currDisc = Math.Abs(A_phi[0][j] - TheorAnswer[j]) / Math.Abs(TheorAnswer[j]);
                
                if (Math.Abs(maxDisc) < Math.Abs(currDisc))
                    maxDisc = currDisc;

                if (!double.IsNaN(currDisc))
                {
                    avgDisc += currDisc;
                    NotNaNamount++;            
                    sumU += absDiff * absDiff;
                    sumD += TheorAnswer[j] * TheorAnswer[j];
                }

                sw_d.WriteLine($"{absDiff:E8} {currDisc:E8}");    
            }
            avgDisc = Math.Sqrt(sumU) / Math.Sqrt(sumD);
            sw_d.WriteLine($"Средняя невязка: {avgDisc:E15}");
            sw_d.WriteLine($"Максимальная невязка: {maxDisc:E15}");
            sw_d.WriteLine($"С: {avgDisc:E7}");
            sw_d.WriteLine($"М: {maxDisc:E7}");
            sw_d.Close();
        }
    }

    public void GenerateVectorEphi()
    {
        if (timeMesh is null) throw new NullReferenceException("timeMesh is null");
        
        E_phi2D = new GlobalVector[A_phi.Length];
        for (int i = 0; i < E_phi2D.Length; i++)
        {
            if (i == 0)
                E_phi2D[i] = (-1.0 / (timeMesh[i + 1] - timeMesh[i])) * (A_phi[i + 1] - A_phi[i]);
            else if (i == E_phi2D.Length - 1)
                E_phi2D[i] = (-1.0 / (timeMesh[i] - timeMesh[i - 1])) * (A_phi[i] - A_phi[i - 1]);
            else
                E_phi2D[i] = (-1.0 / (timeMesh[i + 1] - timeMesh[i - 1])) * (A_phi[i + 1] - A_phi[i - 1]);
        }
    }

    public List<int> GetE_phi(double r, double z, double t)
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