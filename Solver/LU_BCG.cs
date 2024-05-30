using System.Diagnostics;
using MathObjects;

namespace Solver;

public class LU_BCG : ISolver
{
    private const int _maxIter = 100_000;

    private const double _eps = 1E-15;

    protected static void PartitialLU(GlobalMatrix Matrix)
    {
      for (int i = 0; i < Matrix.Size; i++)
      {

         for (int j = Matrix._ig[i]; j < Matrix._ig[i + 1]; j++)
         {
            int jCol = Matrix._jg[j];
            int jk = Matrix._ig[jCol];
            int k = Matrix._ig[i];

            int sdvig = Matrix._jg[Matrix._ig[i]] - Matrix._jg[Matrix._ig[jCol]];

            if (sdvig > 0)
               jk += sdvig;
            else
               k -= sdvig;

            double sumL = 0.0;
            double sumU = 0.0;

            for (; k < j && jk < Matrix._ig[jCol + 1]; k++, jk++)
            {
               sumL += Matrix._al[k] * Matrix._au[jk];
               sumU += Matrix._au[k] * Matrix._al[jk];
            }

            Matrix._al[j] -= sumL;
            Matrix._au[j] -= sumU;
            Matrix._au[j] /= Matrix._diag[jCol];
         }

         double sumD = 0.0;
         for (int j = Matrix._ig[i]; j < Matrix._ig[i + 1]; j++)
            sumD += Matrix._al[j] * Matrix._au[j];

         Matrix._diag[i] -= sumD;
      }
   }

   protected static GlobalVector Forward(GlobalMatrix Matrix, GlobalVector b)
   {
      var result = new GlobalVector(b);

      for (int i = 0; i < Matrix.Size; i++)
      {
         for (int j = Matrix._ig[i]; j < Matrix._ig[i + 1]; j++)
         {
            result[i] -= Matrix._al[j] * result[Matrix._jg[j]];
         }

         result[i] /= Matrix._diag[i];
      }

      return result;
   }

    protected static GlobalVector Backward(GlobalMatrix Matrix, GlobalVector b)
    {
      var result = new GlobalVector(b);

      for (int i = Matrix.Size - 1; i >= 0; i--)
      {
         for (int j = Matrix._ig[i + 1] - 1; j >= Matrix._ig[i]; j--)
         {
            result[Matrix._jg[j]] -= Matrix._au[j] * result[i];
         }
      }

      return result;
   }

    public (GlobalVector, GlobalVector) Solve(GlobalMatrix A, GlobalVector b)
    {
        GlobalMatrix LU = new(A);
        PartitialLU(LU);
        var reVector = Forward(LU, b);

        var _solution = new GlobalVector(b.Size);
        GlobalVector residual = reVector - Forward(LU, A * Backward(LU, _solution));

        GlobalVector p = new(residual);
        GlobalVector z = new(residual);
        GlobalVector s = new(residual);

        Stopwatch sw = Stopwatch.StartNew();

        //double vecNorm = _vector.Norm();
        double vecNorm = reVector.Norma();
        double discrepancy = 1;
        double prPrev = p * residual;

         int i = 1;
        for (; i <= _maxIter && discrepancy > _eps; i++)
        {
            var L_AU_z = Forward(LU, A * Backward(LU, z));
            double alpha = prPrev / (s * L_AU_z);
            
            _solution = _solution + alpha * z;
            residual = residual - alpha * L_AU_z;
            
            //var L_ATU_s = SparseMatrix.TransposedMatrixMult(_matrix, s);
            //var L_ATU_s = Backward(reMatrixLU, _matrix * Forward(reMatrixLU, s));
            var L_ATU_s = Forward(LU, A.Transpose() * Backward(LU, s));
            p = p - alpha * L_ATU_s;
            
            double pr = p * residual;
            double beta = pr / prPrev;
            prPrev = pr;
            
            z = residual + beta * z;
            s = p + beta * s;
            discrepancy = residual.Norma() / vecNorm;
        }

      _solution = Backward(LU, _solution);

      sw.Stop();
      var SolvationTime = sw.ElapsedMilliseconds;

        Console.WriteLine(
            $"Computing finished!\n" +
            $"Total iterations: {i}\n" +
            $"Time ms: {sw.ElapsedMilliseconds}\n" +
            $"Relative residuality: {residual.Norma() / b.Norma():E15}\n");

      return (_solution, _solution);
    }
}