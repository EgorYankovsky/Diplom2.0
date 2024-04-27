using System.Collections.Immutable;
using System.Diagnostics;
using DataStructs;
using System.Reflection.Metadata;

namespace Grid;

public static class MeshReader
{
    private static readonly double mu0 = 4.0D * Math.PI * Math.Pow(10.0D, -7);

    public static double R;

    public static List<double> NodesR = [];

    public static string InfoAboutR = string.Empty;

    public static List<double> NodesX = [];

    public static string InfoAboutX = string.Empty;

    public static List<double> NodesY = [];

    public static string InfoAboutY = string.Empty;

    public static List<double> NodesZ = [];

    public static string InfoAboutZ = string.Empty;

    public static List<Elem> Elems = [];

    public static List<Border3D> Borders = [];

    public static (double, double) Time;

    public static int tn;

    public static double tk;

    public static void ReadMesh(string meshPath, string bordersPath)
    {
        var fileData = File.ReadAllText(meshPath).Split("\n");
        
        NodesX = fileData[0].Split(" ").Select(double.Parse).ToList();
        InfoAboutX = fileData[1];
        
        NodesY = fileData[2].Split(" ").Select(double.Parse).ToList();
        InfoAboutY = fileData[3];
        
        NodesZ = fileData[4].Split(" ").Select(double.Parse).ToList();
        InfoAboutZ = fileData[5];

        NodesR = fileData[6].Split(" ").Select(double.Parse).ToList();
        InfoAboutR = fileData[7];

        R = double.Parse(fileData[8]);

        var elemsAmount = int.Parse(fileData[9]);
        for (int i = 0; i < elemsAmount; i++)
        {
            var info = fileData[10 + i].Split(" ");
            Elems.Add(new Elem(int.Parse(info[0]), int.Parse(info[1]),
                               int.Parse(info[2]), int.Parse(info[3]), 
                               int.Parse(info[4]), int.Parse(info[5]), 
                               int.Parse(info[6]), mu0, double.Parse(info[7])));
        }

        fileData = File.ReadAllText(bordersPath).Split("\n");

        var bordersAmount = int.Parse(fileData[0]);
        for (int i = 0; i < bordersAmount; i++)
        {
            var info = fileData[1 + i].Split(" ").Select(int.Parse).ToArray();
            Borders.Add(new Border3D(info[0], info[1], info[2], info[3],
                                     info[4], info[5], info[6], info[7]));
        }
        Debug.WriteLine("All data read correctly");
    }

    public static void ReadTimeMesh(string timePath)
    {
        var info = File.ReadAllText(timePath).Split("\n");
        Time.Item1 = double.Parse(info[0].Split(" ")[0]);
        Time.Item2 = double.Parse(info[0].Split(" ")[1]);
        tn = int.Parse(info[0].Split(" ")[2]);
        tk = double.Parse(info[0].Split(" ")[3]);
        Debug.WriteLine("Time data read correctly");
    }

    public static List<Layer> ReadFields(string path)
    {
        var layers = new List<Layer>();
        var info = File.ReadAllText(path).Split("\n");
        var fieldsAmount = int.Parse(info[0]);
        for (int i = 0; i < fieldsAmount; i++)
        {
            var fieldInfo = info[i + 1].Split(" ");
            layers.Add(new Layer(double.Parse(fieldInfo[4]), double.Parse(fieldInfo[5]), mu0, double.Parse(fieldInfo[6])));
        }
        return layers;
    }
}