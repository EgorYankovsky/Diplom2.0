using MathObjects;

namespace Solver;

public interface ISolver
{
    public (GlobalVector, GlobalVector) Solve(GlobalMatrix A, GlobalVector b);
}
