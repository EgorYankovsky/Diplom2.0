using System.Diagnostics;
using MathObjects;

namespace Solver;

public class BCG : ISolver
{
    private const int _maxIter = 100_000;

    private const double _eps = 1E-15;

    public (GlobalVector, GlobalVector) Solve(GlobalMatrix A, GlobalVector b)
    {
        GlobalVector x = new(b.Size);
        GlobalVector x_ = new(b.Size);

        GlobalVector r = new(b.Size);
        GlobalVector r_ = new(b.Size);

        GlobalVector z = new(b.Size);
        GlobalVector z_ = new(b.Size);

        GlobalVector p = new(b.Size);
        GlobalVector p_ = new(b.Size);

        GlobalVector s = new(b.Size);
        GlobalVector s_ = new(b.Size);


        double alph = 0.0D;
        double beta = 0.0D;

        r = b - A * x_;
        p = r;
        z = r;
        s = r;

        int iter = 1;
        Stopwatch sw = Stopwatch.StartNew();
        do
        {
            x_ = x;
            z_ = z;
            r_ = r;
            s_ = s;
            p_ = p;

            var Az = A * z_;
            alph = (p_ * r_) / (s_ * Az);
            
            x = x_ + alph * z_;
            r = r_ - alph * Az;
            p = p_ - alph * A.Transpose() * s_;

            beta = (p * r) / (p_ * r_);

            z = r + beta * z_;
            s = p + beta * s_;

            iter++;
            Console.WriteLine($"{r.Norma() / b.Norma():E15}");
        } while (iter < _maxIter && r.Norma() / b.Norma() >= _eps);
        sw.Stop();

        Console.WriteLine(
            $"Computing finished!\n" +
            $"Total iterations: {iter}\n" +
            $"Time ms: {sw.ElapsedMilliseconds}\n" +
            $"Relative residuality: {r.Norma() / b.Norma():E15}\n");
        return (x, x_);
    }
}

/*
public BCG(double eps = 1e-14, int maxIters = 2000)
   {
      Eps = eps;
      MaxIters = maxIters;
   }


   public override Vector Solve()
   {
      _solution = new(_vector.Size);

      Vector residual = _vector - _matrix * _solution;

      Vector p = new(residual.Size);
      Vector z = new(residual.Size);
      Vector s = new(residual.Size);

      Vector.Copy(residual, p);
      Vector.Copy(residual, z);
      Vector.Copy(residual, s);


      Stopwatch sw = Stopwatch.StartNew();

      double vecNorm = _vector.Norm();
      double discrepancy = 1;
      double prPrev = p * residual;

      for (int i = 1; i <= MaxIters && discrepancy > Eps; i++)
      {
         var Az = _matrix * z;
         double alpha = prPrev / (s * Az);

         _solution = _solution + alpha * z;
         residual = residual - alpha * Az;
         p = p - alpha * SparseMatrix.TransposedMatrixMult(_matrix, s);

         double pr = p * residual;
         double beta = pr / prPrev;
         prPrev = pr;

         z = residual + beta * z;
         s = p + beta * s;

         discrepancy = residual.Norm() / vecNorm;
      }

      sw.Stop();
      SolvationTime = sw.ElapsedMilliseconds;

      return _solution;
*/