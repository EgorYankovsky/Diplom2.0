using System.Collections.Immutable;

namespace Grid;

public static class MeshReader
{
    public static void ReadMesh(string path1, string path2, ref Mesh2Dim mesh)
    {
        string _currPath = path1;
        try
        {
            // Считывание параметров для расчетной области.
            using (var sr = new StreamReader(_currPath))
            {
                mesh.nodesR = sr.ReadLine().Split().Select(double.Parse).ToList();
                for (int i = 0; i < mesh.nodesR.Count; i++)
                    mesh.nodesR_Refs.Add(i);
                mesh.NodesRWithoutFragmentation = mesh.nodesR.ToImmutableArray();
                mesh.infoAboutR = sr.ReadLine() ?? "";

                mesh.nodesZ = sr.ReadLine().Split().Select(double.Parse).ToList();
                for (int i = 0; i < mesh.nodesZ.Count; i++)
                    mesh.nodesZRefs.Add(i);
                mesh.NodesZWithoutFragmentation = mesh.nodesZ.ToImmutableArray();
                mesh.infoAboutZ = sr.ReadLine() ?? "";


                mesh.ElemsAmount = int.Parse(sr.ReadLine() ?? "0");
                for (int i = 0; i < mesh.ElemsAmount; i++)
                {
                    string[] str = sr.ReadLine().Split();
                    mesh.mu0.Add(double.Parse(str[5]));
                    mesh.sigma.Add(double.Parse(str[6]));
                    mesh.Elems.Add(new List<int> {int.Parse(str[0]), int.Parse(str[1]),
                                             int.Parse(str[2]), int.Parse(str[3]), 
                                             int.Parse(str[4])});
            
                }
            }

            // Считывание данных для границ.
            _currPath = path2;
            using (var sr = new StreamReader(_currPath))
            {
                mesh.bordersAmount = int.Parse(sr.ReadLine() ?? "0");
                for (int i = 0; i < mesh.bordersAmount; i++)
                    mesh.borders.Add(sr.ReadLine().Split().Select(int.Parse).ToList());

            }
        }
        catch(IOException ex)
        {
            Console.WriteLine($"Error during reading {_currPath} file. Wrong data format: {ex}");
            throw new IOException();
        }
        catch(Exception ex)
        {
            Console.WriteLine($"Error during reading {_currPath} file: {ex}");
            throw new IOException();
        }
    }
}