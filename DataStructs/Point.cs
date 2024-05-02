namespace DataStructs;
using System.Globalization;

public enum Location
{
    Inside,
    OutSide,
    BoundaryI,
    BoundaryII,
    BoundaryIII,
    NotStated
}

public interface IPoint
{
    public int SubElemNum { get; set; }
    
    public Location Type { get; set; }

}

public class Point2D : IPoint
{    
    public int SubElemNum { get; set; }
    
    public Location Type { get; set; }
    
    public double R { get; init; }

    public double Z { get; init; }

    public override string ToString() => $"{R.ToString("E15", CultureInfo.InvariantCulture)} {Z.ToString("E15", CultureInfo.InvariantCulture)} {SubElemNum} {Type}";

    public Point2D(double R, double Z)
    {
        this.R = R;
        this.Z = Z;
        Type = Location.NotStated;
    }

    public Point2D(List<string> arr)
    {
        R = double.Parse(arr[0]);
        Z = double.Parse(arr[1]);
        SubElemNum = int.Parse(arr[2]);
        Type = arr[3] switch
        {
            "Inside" => Location.Inside,
            "OutSide" => Location.OutSide,
            "BoundaryI" => Location.BoundaryI,
            "BoundaryII" => Location.BoundaryII,
            "BoundaryIII" => Location.BoundaryIII,
            _ => Location.NotStated,
        };
    }
}

public class Point3D : IPoint
{
    public int SubElemNum { get; set; }
    
    public Location Type { get; set; }

    public double X { get; init; }

    public double Y { get; init; }

    public double Z { get; init; }

    public Point3D(double X, double Y, double Z)
    {
        this.X = X;
        this.Y = Y;
        this.Z = Z;
        Type = Location.NotStated;
    }

    public Point3D(List<string> arr)
    {
        X = double.Parse(arr[0]);
        Y = double.Parse(arr[1]);
        Z = double.Parse(arr[2]);
        SubElemNum = int.Parse(arr[3]);
        Type = arr[4] switch
        {
            "Inside" => Location.Inside,
            "OutSide" => Location.OutSide,
            "BoundaryI" => Location.BoundaryI,
            "BoundaryII" => Location.BoundaryII,
            "BoundaryIII" => Location.BoundaryIII,
            _ => Location.NotStated,
        };
    }

    public override string ToString() => $"{X.ToString("E15", CultureInfo.InvariantCulture)} {Y.ToString("E15", CultureInfo.InvariantCulture)} {Z.ToString("E15", CultureInfo.InvariantCulture)} {SubElemNum} {Type}";
}