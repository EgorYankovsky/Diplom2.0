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

    public void Solve()
    {
        if (pointsArr is null) throw new ArgumentNullException("points array is null !");
        if (elemsArr is null) throw new ArgumentNullException("elems array is null !");
        if (bordersArr is null) throw new ArgumentNullException("borders array is null !");
        if (Solutions is null) throw new ArgumentNullException("solutions array is null !");
        if (Discrepancy is null) throw new ArgumentNullException("discrepancy array is null !");
        if (solver is null) throw new ArgumentNullException("solver is null !");

        Debug.WriteLine($"\nTime layer: before BC");
        Thread.Sleep(1500);

        Matrix = new GlobalMatrix(pointsArr.Length);
        Generator.BuildPortait(ref Matrix, pointsArr.Length, elemsArr);
        Generator.FillMatrix(ref Matrix, pointsArr, elemsArr, TypeOfMatrixM.Mrr);
        Generator.ConsiderBoundaryConditions(ref Matrix, bordersArr);

        Vector = new GlobalVector(pointsArr.Length);
        Generator.FillVector(ref Vector, pointsArr, elemsArr, 0.0);
        Generator.ConsiderBoundaryConditions(ref Vector, bordersArr, pointsArr, 0.0D);

        (Solutions[0], Discrepancy[0]) = solver.Solve(Matrix, Vector);

        if (equationType == EquationType.Parabolic)
        {
            (Solutions[1], Discrepancy[1]) = (Solutions[0], Discrepancy[0]);

            for (int i = 2; i < timeMesh.Length; i++)
            {
                Debug.WriteLine($"\nTime layer: {timeMesh[i]}");
                Thread.Sleep(1500);

                double deltT = timeMesh[i] - timeMesh[i - 2];
                double deltT0 = timeMesh[i] - timeMesh[i - 1];
                double deltT1 = timeMesh[i - 1] - timeMesh[i - 2];

                double tau0 = (deltT + deltT0) / (deltT * deltT0);
                double tau1 = deltT / (deltT1 * deltT0);
                double tau2 = deltT0 / (deltT * deltT0);

                var matrix1 = new GlobalMatrix(pointsArr.Length);
                Generator.BuildPortait(ref matrix1, pointsArr.Length, elemsArr);
                Generator.FillMatrix(ref matrix1, pointsArr, elemsArr, TypeOfMatrixM.Mrr);
                Generator.ConsiderBoundaryConditions(ref matrix1, bordersArr);

                var M = new GlobalMatrix(pointsArr.Length); // ???
                Generator.BuildPortait(ref M, pointsArr.Length, elemsArr);
                Generator.FillMatrix(ref M, pointsArr, elemsArr, TypeOfMatrixM.Mr);
                Generator.ConsiderBoundaryConditions(ref M, bordersArr);

                Matrix = (tau0 * M) + matrix1;
                Generator.ConsiderBoundaryConditions(ref Matrix, bordersArr);

                var bi = new GlobalVector(pointsArr.Length);
                Vector = bi - (tau2 * (M * Solutions[i - 2])) + (tau1 * (M * Solutions[i - 1]));
                Generator.ConsiderBoundaryConditions(ref Vector, bordersArr, pointsArr, timeMesh[i]);
                
                (Solutions[i], Discrepancy[i]) = solver.Solve(Matrix, Vector);
            }
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
        if (E_phi2D is null) throw new ArgumentNullException();

        if (timeMesh is not null)
        {
            for (int i = 0; i < timeMesh.Length; i++)
            {
                using var sw = new StreamWriter($"{_path}\\A_phi\\Answer\\Answer_Aphi_time={timeMesh[i]}.dat");
                for (int j = 0;   j < A_phi[i].Size; j++)
                    sw.WriteLine($"{A_phi[i][j]:E8}");
                sw.Close();
            }
            for (int i = 0; i < timeMesh.Length; i++)
            {
                using var sw = new StreamWriter($"{_path}\\E_phi\\Answer\\Answer_Ephi_time={timeMesh[i]}.dat");
                for (int j = 0;   j < E_phi2D[i].Size; j++)
                    sw.WriteLine($"{E_phi2D[i][j]:E8}");
                sw.Close();
            }
        }
        else
        {
            using var sw = new StreamWriter($"{_path}\\A_phi\\Answer\\Answer.dat");
            for (int j = 0; j < A_phi[0].Size; j++)
                sw.WriteLine($"{A_phi[0][j]:E8}");
            sw.Close();

            using var sw1 = new StreamWriter($"{_path}\\E_phi\\Answer\\Answer.dat");
            for (int j = 0; j < E_phi2D[0].Size; j++)
                sw1.WriteLine($"{E_phi2D[0][j]:E8}");
            sw1.Close();
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
        if (E_phi2D is null) throw new ArgumentNullException();


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

                for (int j = 0; j < A_phi[i].Size; j++)
                {
                    double absDiff = A_phi[i][j] - Discrepancy[i][j];
                    double currDisc = Math.Abs(A_phi[i][j] - Discrepancy[i][j]) / Math.Abs(A_phi[i][j]);
                    
                    if (Math.Abs(maxDisc) < Math.Abs(currDisc))
                        maxDisc = currDisc;

                    if (!double.IsNaN(currDisc))
                    {
                        avgDisc += currDisc;
                        NotNaNamount++;
                    }

                    sumU += absDiff * absDiff;
                    sumD += A_phi[i][j] * A_phi[i][j];

                    sw_d.WriteLine($"{absDiff:E15} {currDisc:E15}");
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

            for (int j = 0; j < A_phi[0].Size; j++)
            {
                double absDiff = A_phi[0][j] - Discrepancy[0][j];
                double currDisc = Math.Abs(A_phi[0][j] - Discrepancy[0][j]) / Math.Abs(A_phi[0][j]);
                
                if (Math.Abs(maxDisc) < Math.Abs(currDisc))
                    maxDisc = currDisc;

                if (!double.IsNaN(currDisc))
                {
                    avgDisc += currDisc;
                    NotNaNamount++;
                }

                sumU += absDiff * absDiff;
                sumD += A_phi[0][j] * A_phi[0][j];

                sw_d.WriteLine($"{absDiff:E15} {currDisc:E15}");
    
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
            else if (i == 1)
                E_phi2D[i] = (-1.0 / (timeMesh[i + 1] - timeMesh[i - 1])) * (A_phi[i + 1] - A_phi[i - 1]);
            else
            {
                double ti = timeMesh[i];
                double ti_1 = timeMesh[i - 1];
                double ti_2 = timeMesh[i - 2];
                E_phi2D[i] = 1.0D / (ti_1 - ti_2) * A_phi[i - 2] - (ti - ti_2) / ((ti_1 - ti_2) * (ti - ti_1)) * A_phi[i - 1] + 
                (2 * ti - ti_1 - ti_2) / ((ti - ti_2) * (ti - ti_1)) * A_phi[i];
            }
        }
        
        /*
        for (int i = 0; i < E_phi2D.Length; i++)
        {
            if (i == 0)
                E_phi2D[i] = (-1.0 / (timeMesh[i + 1] - timeMesh[i])) * (A_phi[i + 1] - A_phi[i]);
            else if (i == E_phi2D.Length - 1)
                E_phi2D[i] = (-1.0 / (timeMesh[i] - timeMesh[i - 1])) * (A_phi[i] - A_phi[i - 1]);
            else
                E_phi2D[i] = (-1.0 / (timeMesh[i + 1] - timeMesh[i - 1])) * (A_phi[i + 1] - A_phi[i - 1]);
        }
        */
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