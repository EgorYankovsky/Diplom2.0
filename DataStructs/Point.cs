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


public class Point
{
    // Номер подобласти.
    public int SubElemNum { get; set; }
    
    // Тип точки
    public Location Type { get; set; }

    // Координата по R.
    public double R { get; init; }

    // Координата по X.
    public double X { get; init; }

    // Координата по Y.
    public double Y { get; init; }

    // Координата по Z.
    public double Z { get; init; }

    /// <summary>
    /// Метод, возвращающий содержимое точки в формате строки.
    /// </summary>
    /// <returns>Форматированная строка.</returns>
    public override string ToString() => $"{R.ToString("E15", CultureInfo.InvariantCulture)} {Z.ToString("E15", CultureInfo.InvariantCulture)} {SubElemNum} {Type}";

    /// <summary>
    /// Конструктор класса Point.
    /// </summary>
    /// <param name="X">Координата по оси X.</param>
    /// <param name="Y">Координата по оси Y.</param>
    /// <param name="Z">Координата по оси Z.</param>
    public Point(double X, double Y, double Z)
    {
        this.X = X;
        this.Y = Y;
        this.Z = Z;
        Type = Location.NotStated;
    }

    /// <summary>
    /// Конструктор класса Point.
    /// </summary>
    /// <param name="R">Координата по оси R.</param>
    /// <param name="Z">Координата по оси Z.</param>
    public Point(double R, double Z)
    {
        this.R = R;
        this.Z = Z;
        Type = Location.NotStated;
    }

    /// <summary>
    /// Конструктор класса Point
    /// </summary>
    /// <param name="arr">Массив с информацией.</param>
    /// <exception cref="Exception"></exception>
    public Point(List<string> arr)
    {
        // TODO: Придумать что-нибудь интереснее.
        if (arr.Count == 5)
        {
            X = double.Parse(arr[0]);
            Y = double.Parse(arr[1]);
            Z = double.Parse(arr[2]);
            SubElemNum = int.Parse(arr[4]); // ! Achtung: die fehler    
            switch (arr[3])
            {
                case "Inside":
                {
                    Type = Location.Inside;
                    break;
                }
                case "OutSide":
                {
                    Type = Location.OutSide;
                    break;
                }
                case "BoundaryI":
                {
                    Type = Location.BoundaryI;
                    break;
                }
                case "BoundaryII":
                {
                    Type = Location.BoundaryII;
                    break;
                }
                case "BoundaryIII":
                {
                    Type = Location.BoundaryIII;
                    break;
                }
                default:
                {
                    Type = Location.NotStated;
                    break;
                }
            }
        }
        else if (arr.Count == 4)
        {
            R = double.Parse(arr[0]);
            Z = double.Parse(arr[1]);
            SubElemNum = int.Parse(arr[2]); // ! Achtung: die fehler    
            switch (arr[3])
            {
                case "Inside":
                {
                    Type = Location.Inside;
                    break;
                }
                case "OutSide":
                {
                    Type = Location.OutSide;
                    break;
                }
                case "BoundaryI":
                {
                    Type = Location.BoundaryI;
                    break;
                }
                case "BoundaryII":
                {
                    Type = Location.BoundaryII;
                    break;
                }
                case "BoundaryIII":
                {
                    Type = Location.BoundaryIII;
                    break;
                }
                default:
                {
                    Type = Location.NotStated;
                    break;
                }
            }
        }
    }
}