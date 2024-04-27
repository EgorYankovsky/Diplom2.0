namespace DataStructs;

public record Border3D(int BorderType, int BorderFormula,
                       int X0, int X1, int Y0, int Y1,
                       int Z0, int Z1);

public record Border2D(int BorderType, int BorderFormula,
                       int R0, int R1, int Z0, int Z1);


public record Elem
{
    public int AreaNum;

    public int[] Arr;

    public double mu;

    public double sigma;

    public override string ToString()
    {
        string ans = string.Empty;
        for (int i = 0; i < Arr.Length; i++)
            ans += $"{Arr[i]} ";
        return ans + $"{mu} " + $"{sigma} ";
    }

    public Elem(int areaNum, int a, int b, int c, int d, int e, int f,
                int g, int h, int i, int j, int k, int l, 
                double mu, double sigma)
    {
        AreaNum = areaNum;
        Arr = [a, b, c, d, e, f, g, h, i, j, k, l];
        this.mu = mu; this.sigma = sigma;
    }

    public Elem(int areaNum, int a, int b, int c, int d,
                int e, int f, int g, int h, double mu,
                double sigma)
    {
        AreaNum = areaNum;
        Arr = [a, b, c, d, e, f, g, h];
        this.mu = mu; this.sigma = sigma;
    }

    public Elem(int areaNum, int a, int b, int c, int d,
                int e, int f, double mu,
                double sigma)
    {
        AreaNum = areaNum;
        Arr = [a, b, c, d, e, f];
        this.mu = mu; this.sigma = sigma;
    }

    public Elem(int areaNum, int a, int b, int c, int d,
                double mu, double sigma)
    {    
        AreaNum = areaNum;
        Arr = [a, b, c, d];
        this.mu = mu; this.sigma = sigma;
    }
}