using DataStructs;

namespace Functions;

public static class Function
{
    // Функция правой части.
    public static double F(double r, double z, double t) => (r, z) == (10.0D, 0.0D) && t <= 1.0D ? Math.Pow(10.0D, 1) : 0.0D;

    // Тестируемая функция.
    public static double U(double r, double z) => r * r * r + z * z * z;

    #region Первые краевые условия.

    public static double U1_1(Point p, double t) => 0.0D;

    public static double U1_2(Point p, double t) => 0.0D;

    public static double U1_3(Point p, double t) => 0.0D;

    public static double U1_4(Point p, double t) => 0.0D;

    #endregion


    #region Вторые краевые условия.
    
    public static double dUdn1_1(Point p) => 0.0D;

    public static double dUdn1_2(Point p) => 0.0D;

    public static double dUdn1_3(Point p) => 0.0D;

    public static double dUdn1_4(Point p) => 0.0D;

    #endregion


    #region Третьи краевые условия.

    public static double U3dUdn3_1(Point p) => throw new Exception("Can't commit III bc.");

    public static double U3dUdn3_2(Point p) => throw new Exception("Can't commit III bc.");

    public static double U3dUdn3_3(Point p) => throw new Exception("Can't commit III bc.");

    public static double U3dUdn3_4(Point p) => throw new Exception("Can't commit III bc.");

    #endregion
}