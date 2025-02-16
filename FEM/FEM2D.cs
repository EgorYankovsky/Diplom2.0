﻿using MathObjects;
using Grid;
using System.Diagnostics;
using Functions;

namespace Project;

public class FEM2D : FEM
{
    public FEM2D(Mesh2Dim mesh, TimeMesh timeMesh) : base(timeMesh, 2)
    {
        mesh2Dim = mesh;
        if (timeMesh[0] == timeMesh[^1])
            equationType = EquationType.Elliptic;
        else
            equationType = EquationType.Parabolic;
        
        A_phi = new GlobalVector[Time.Count];
        E_phi = new GlobalVector[Time.Count];
        Debug.WriteLine("Generated data submited");
    }

    private Mesh2Dim mesh2Dim;

    public GlobalVector[] A_phi;

    public GlobalVector[] E_phi;

    public void Solve()
    {
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

            for (int i = 2; i < Time.Count; i++)
            {
                Debug.WriteLine($"\nTime layer: {Time[i]}");
                Thread.Sleep(1500);

                double deltT = Time[i] - Time[i - 2];
                double deltT0 = Time[i] - Time[i - 1];
                double deltT1 = Time[i - 1] - Time[i - 2];

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
                Generator.ConsiderBoundaryConditions(ref Vector, bordersArr, pointsArr, Time[i]);
                
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
        if (E_phi is null) throw new ArgumentNullException();

        if (Time.Count != 1)
        {
            for (int i = 0; i < Time.Count; i++)
            {
                using var sw = new StreamWriter($"{_path}\\A_phi\\Answer\\Answer_Aphi_time={Time[i]}.dat");
                for (int j = 0;   j < A_phi[i].Size; j++)
                    sw.WriteLine($"{A_phi[i][j]:E8}");
                sw.Close();
            }
            for (int i = 0; i < Time.Count; i++)
            {
                using var sw = new StreamWriter($"{_path}\\E_phi\\Answer\\Answer_Ephi_time={Time[i]}.dat");
                for (int j = 0;   j < E_phi[i].Size; j++)
                    sw.WriteLine($"{E_phi[i][j]:E8}");
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
            for (int j = 0; j < E_phi[0].Size; j++)
                sw1.WriteLine($"{E_phi[0][j]:E8}");
            sw1.Close();
        }
    }

    public void WriteDiscrepancy(string _path)
    {
        if (A_phi is null) throw new ArgumentNullException();
        if (E_phi is null) throw new ArgumentNullException();


        if (Time.Count != 1)
        {
            for (int i = 0; i < Time.Count; i++)
            {
                using var sw_d = new StreamWriter($"{_path}\\A_phi\\Discrepancy\\Discrepancy_Aphi_time={Time[i]}.dat");

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
        E_phi = new GlobalVector[A_phi.Length];   
        for (int i = 0; i < E_phi.Length; i++)
        {
            if (i == 0 || i == 1)
                E_phi[i] = new GlobalVector(A_phi[i].Size);
            else
            {
                double ti = Time[i];
                double ti_1 = Time[i - 1];
                double ti_2 = Time[i - 2];
                E_phi[i] = 1.0D / (ti_1 - ti_2) * A_phi[i - 2] - (ti - ti_2) / ((ti_1 - ti_2) * (ti - ti_1)) * A_phi[i - 1] + 
                (2 * ti - ti_1 - ti_2) / ((ti - ti_2) * (ti - ti_1)) * A_phi[i];
            }
        }

    }

    internal List<int> GetE_phi(double r, double z, double t)
    {
        int i;
        for (i = 0; i < mesh2Dim.nodesR.Count - 1; i++)
            if (mesh2Dim.nodesR[i] <= r && r <= mesh2Dim.nodesR[i + 1])
                break;
        int j;
        for (j = 0; j < mesh2Dim.nodesZ.Count - 1; j++)
            if (mesh2Dim.nodesZ[j] <= z && z <= mesh2Dim.nodesZ[j + 1])
                break;

        return elemsArr[j * (mesh2Dim.nodesR.Count - 1) + i];
    }
}