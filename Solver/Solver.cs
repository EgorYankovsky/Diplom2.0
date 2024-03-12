using MathObjects;

namespace Solver;

public interface ISolver
{
    public GlobalVector Solve(GlobalMatrix A, GlobalVector b);
}
