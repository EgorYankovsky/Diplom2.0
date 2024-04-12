using static Functions.Function;

namespace MathObjects;


public class LocalVector3D : Vector
{
    private readonly double x0;
    private readonly double x1;
    private readonly double y0;
    private readonly double y1;
    private readonly double z0;
    private readonly double z1;
    private readonly double xm;
    private readonly double ym;
    private readonly double zm;
    private readonly double t;

    private LocalMatrixM3D M;

    public override double this[int i] 
    { 
        get 
        {
            double ScalarMult((double, double, double) a, (double, double, double) b) 
            => a.Item1 * b.Item1 + a.Item2 * b.Item2 + a.Item3 * b.Item3;

            (double, double, double) Divide ((double, double, double) v, double sc)
            => (v.Item1 / sc, v.Item2 / sc, v.Item3 / sc); 

            double Length((double, double, double) v) 
            => Math.Sqrt(Math.Pow(v.Item1, 2) + Math.Pow(v.Item2, 2) + Math.Pow(v.Item3, 2));

            // Очень спорно.
            List<double> vect = [
                ScalarMult(F(x0, ym, z0, t), Divide((0.0D, y1 - y0, 0.0D), Length((0.0D, y1 - y0, 0.0D)))),
                ScalarMult(F(x1, ym, z0, t), Divide((0.0D, y1 - y0, 0.0D), Length((0.0D, y1 - y0, 0.0D)))),
                ScalarMult(F(x0, ym, z1, t), Divide((0.0D, y1 - y0, 0.0D), Length((0.0D, y1 - y0, 0.0D)))),
                ScalarMult(F(x1, ym, z1, t), Divide((0.0D, y1 - y0, 0.0D), Length((0.0D, y1 - y0, 0.0D)))),
                
                ScalarMult(F(xm, y0, z0, t), Divide((x1 - x0, 0.0D, 0.0D), Length((x1 - x0, 0.0D, 0.0D)))),
                ScalarMult(F(xm, y1, z0, t), Divide((x1 - x0, 0.0D, 0.0D), Length((x1 - x0, 0.0D, 0.0D)))),
                ScalarMult(F(xm, y0, z1, t), Divide((x1 - x0, 0.0D, 0.0D), Length((x1 - x0, 0.0D, 0.0D)))),
                ScalarMult(F(xm, y1, z1, t), Divide((x1 - x0, 0.0D, 0.0D), Length((x1 - x0, 0.0D, 0.0D)))),
                
                ScalarMult(F(x0, y0, zm, t), Divide((0.0D, 0.0D, z1 - z0), Length((0.0D, 0.0D, z1 - z0)))),
                ScalarMult(F(x1, y0, zm, t), Divide((0.0D, 0.0D, z1 - z0), Length((0.0D, 0.0D, z1 - z0)))),
                ScalarMult(F(x0, y1, zm, t), Divide((0.0D, 0.0D, z1 - z0), Length((0.0D, 0.0D, z1 - z0)))),
                ScalarMult(F(x1, y1, zm, t), Divide((0.0D, 0.0D, z1 - z0), Length((0.0D, 0.0D, z1 - z0)))) 
            ];
            
            double ans = 0.0D;
            for (int j = 0; j < vect.Count; j++)
                ans += M[i, j] * vect[j];
            return ans;
        }
    }

    public LocalVector3D(double x0, double x1, double y0, double y1, double z0, double z1, double t)
    {
        this.x0 = x0;
        this.x1 = x1;
        this.y0 = y0;
        this.y1 = y1;
        this.z0 = z0;
        this.z1 = z1;
        this.t = t;
        xm = 0.5D * (x1 + x0);
        ym = 0.5D * (y1 + y0);
        zm = 0.5D * (z1 + z0);
        M = new(1.0, x1 - x0, y1 - y0, z1 - z0);
    }
}