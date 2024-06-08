using Project;
using System.Globalization;
using Solver;
using Processor;
using Grid;
using static Grid.MeshReader;
using static Grid.MeshGenerator;
using static Manager.FolderManager;
using DataStructs;

CultureInfo.CurrentCulture = CultureInfo.InvariantCulture;

string InputDirectory = Path.GetFullPath("../../../../Data/Input/");
string SubtotalsDirectory = Path.GetFullPath("../../../../Data/Subtotals/");
string OutputDirectory = Path.GetFullPath("..\\..\\..\\..\\Data\\Output\\");
string PicturesDirectory = Path.GetFullPath("../../../../Drawer/Pictures/");

bool isSolving2DimTask = Checker(OutputDirectory);

// Pre-processor. Clearing output folders.
if (isSolving2DimTask)
    ClearFolders(new List<string> {SubtotalsDirectory + "/2_dim/",
                                   PicturesDirectory + "/E_phi/",
                                   PicturesDirectory + "/A_phi/",
                                   OutputDirectory});

//ClearFolders(new List<string> {SubtotalsDirectory + "3_dim\\",
//                               OutputDirectory + "A_phi\\Answer3D\\"});


// Reading mesh.
ReadMesh(InputDirectory + "WholeMesh.txt");
ReadTimeMesh(InputDirectory + "Time.txt");

// Set recivers
//List<Point3D> recivers = [new(1000.0, 0.0, 0.0), new(2000.0, 0.0, 0.0), 
//                          new(3500.0, 0.0, 0.0), new(5000.0, 0.0, 0.0)];

//List<Point3D> recivers = [new(707.1067, -707.1067, 0.0), new(1414.2135, -1414.2135, 0.0),
//                          new(2474.8737, -2474.8737, 0.0), new(3535.5339, -3535.5339, 0.0)];

List<Point3D> recivers = [new(707.1067, 707.1067, 0.0), new(1414.2135, 1414.2135, 0.0),
                          new(2474.8737, 2474.8737, 0.0), new(3535.5339, 3535.5339, 0.0)];


Mesh3Dim mesh3D = new(NodesX, InfoAboutX, NodesY, InfoAboutY,
                      NodesZ, InfoAboutZ, Elems, Borders);

Mesh2Dim mesh2D = new(NodesR, InfoAboutR, NodesZ, InfoAboutZ,
                      Elems, Math.Sqrt(Math.Pow(mesh3D.nodesX[^1], 2) + Math.Pow(mesh3D.nodesY[^1], 2)));
mesh2D.SetBorders(mesh3D.borders);

var timeMesh = GenerateTimeMesh(Time.Item1, Time.Item2, tn, tk);

// Main process of 2-dim task.
ConstructMesh(ref mesh2D);
FEM2D myFEM2D = new(mesh2D, timeMesh);
if (isSolving2DimTask)
{
    myFEM2D.SetSolver(new LU_LOS());
    myFEM2D.Solve();
    myFEM2D.GenerateVectorEphi();
    myFEM2D.WriteData(OutputDirectory);
    myFEM2D.WritePointsToDraw(OutputDirectory + "ToDraw\\2_dim\\Aphi\\",
                              OutputDirectory + "ToDraw\\2_dim\\Ephi\\");
    myFEM2D.MeasureValuesOnReceivers(recivers, OutputDirectory + "ToDraw\\2_dim\\Receivers\\");
}
else
{
    myFEM2D.ReadAnswer(OutputDirectory);
    Console.WriteLine("2D answer read");
}
//return 0;

myFEM2D.MeasureValuesOnReceivers(recivers, OutputDirectory + "ToDraw\\2_dim\\Receivers\\");

ConstructMesh(ref mesh3D);
Console.WriteLine("3D mesh constructed");
FEM3D myFEM3D = new(mesh3D, timeMesh);
myFEM3D.ConvertResultTo3Dim(myFEM2D);
Console.WriteLine("2D answer converted to 3D");
//myFEM3D.CheckSolution(recivers);
//myFEM3D.WriteData(OutputDirectory + "A_phi\\Answer3D\\ConvertedTo3D\\");

// Solving first layer: groundwater.
ReadAnomaly(InputDirectory + "Anomalies\\Anomaly1.txt");
Console.WriteLine("Anomaly read");

Mesh3Dim mesh3D_a1 = new(NodesX, InfoAboutX, NodesY, InfoAboutY,
                         NodesZ, InfoAboutZ, Elems, Borders);
mesh3D_a1.CommitAnomalyBorders(FieldBorders);
ConstructMeshAnomaly(ref mesh3D_a1, SubtotalsDirectory + "3_dim\\Anomaly0\\");
Console.WriteLine("Anomaly mesh built");
FEM3D fem3D_a1 = new(mesh3D_a1, timeMesh, myFEM3D, 0);
//fem3D_a1.ReadData(OutputDirectory + "A_phi/Answer3D/AfterField1/");
//fem3D_a1.SetSolver(new LOS(100_000, 1E-10));
fem3D_a1.SetSolver(new LU_LOS());
Console.WriteLine("Solving begun");
fem3D_a1.Solve();
fem3D_a1.WriteData(OutputDirectory + $"A_phi\\Answer3D\\AfterField1\\");
//fem3D_a1.ReadData(OutputDirectory + $"A_phi\\Answer3D\\AfterField1\\");
Console.WriteLine("Solved");
fem3D_a1.GenerateVectorE();
Console.WriteLine("E generated");
myFEM3D.AddSolution(fem3D_a1);
//myFEM3D.WriteDataToDraw2DimSolution(OutputDirectory + "ToDraw\\3_dim\\");
//myFEM3D.MeasureValuesOnReceivers(recivers, OutputDirectory + "ToDraw\\3_dim\\Receivers\\");
//
//fem3D_a1.WriteDataToDraw2DimSolution(OutputDirectory + "ToDraw\\3_dim\\");
//int a = Postprocessor.DrawA3D(OutputDirectory + "ToDraw/3_dim/A/", PicturesDirectory + "A3d\\");
//int e = Postprocessor.DrawE3D(OutputDirectory + "ToDraw/3_dim/E/", PicturesDirectory + "E3d\\");
//Console.WriteLine($"{e}");

return 0;
// Solving second layer: groundwater.
ReadAnomaly(InputDirectory + "Anomalies\\Anomaly2.txt");
Console.WriteLine("Anomaly read");

Mesh3Dim mesh3D_a2 = new(NodesX, InfoAboutX, NodesY, InfoAboutY,
                         NodesZ, InfoAboutZ, Elems, Borders);
mesh3D_a2.CommitAnomalyBorders(FieldBorders);
ConstructMeshAnomaly(ref mesh3D_a2, SubtotalsDirectory + "3_dim\\Anomaly1\\");
Console.WriteLine("Anomaly mesh built");
FEM3D fem3D_a2 = new(mesh3D_a2, timeMesh, myFEM3D, 1);
//fem3D_a1.ReadData(OutputDirectory + "A_phi/Answer3D/AfterField1/");
fem3D_a2.SetSolver(new LOS(100_000, 1E-10));

//fem3D_a1.SetSolver(new LU_LOS());
Console.WriteLine("Solving begun");
//fem3D_a2.Solve();
//fem3D_a2.WriteData(OutputDirectory + $"A_phi\\Answer3D\\AfterField2\\");
fem3D_a2.ReadData(OutputDirectory + $"A_phi\\Answer3D\\AfterField2\\");
Console.WriteLine("Solved");
fem3D_a2.GenerateVectorE();
Console.WriteLine("E generated");
myFEM3D.AddSolution(fem3D_a2);
myFEM3D.WriteDataToDraw2DimSolution(OutputDirectory + "ToDraw\\3_dim\\");
myFEM3D.MeasureValuesOnReceivers(recivers, OutputDirectory + "ToDraw\\3_dim\\Receivers\\");
//fem3D_a1.WriteDataToDraw2DimSolution(OutputDirectory + "ToDraw\\3_dim\\");
//int a = Postprocessor.DrawA3D(OutputDirectory + "ToDraw/3_dim/A/", PicturesDirectory + "A3d\\");
int e = Postprocessor.DrawE3D(OutputDirectory + "ToDraw/3_dim/E/", PicturesDirectory + "E3d\\");
Console.WriteLine($"{e}");

// Solving both anomalies.
//ReadBothAnomalies(InputDirectory + "Anomalies\\AnomalyBoth.txt");
//Console.WriteLine("Anomalies read");
//Mesh3Dim mesh3D_ab = new(NodesX, InfoAboutX, NodesY, InfoAboutY,
//                         NodesZ, InfoAboutZ, Elems, Borders);
//mesh3D_ab.CommitAnomalyBorders(FieldBorders);
//mesh3D_ab.CommitSecondAnomalyBorders(FieldBorders1);
//ConstructMeshAnomaly(ref mesh3D_ab, SubtotalsDirectory + "3_dim\\Anomaly2\\");
//Console.WriteLine("Anomalies mesh built");
//FEM3D fem3D_ab = new(mesh3D_ab, timeMesh, myFEM3D, 2);
//fem3D_ab.SetSolver(new LOS(100_000, 1E-12));
////fem3D_ab.SetSolver(new LU_LOS(100_000, 1E-12));
//Console.WriteLine("Solving begun");
//fem3D_ab.Solve();
//fem3D_ab.WriteData(OutputDirectory + $"A_phi\\Answer3D\\AfterFieldBoth\\");
////fem3D_ab.ReadData(OutputDirectory + $"A_phi\\Answer3D\\AfterFieldBoth\\");
//Console.WriteLine("Solved");
//fem3D_ab.GenerateVectorE();
//Console.WriteLine("E generated");
//myFEM3D.AddSolution(fem3D_ab);
//myFEM3D.WriteDataToDraw2DimSolution(OutputDirectory + "ToDraw\\3_dim\\");
//myFEM3D.MeasureValuesOnReceivers(recivers, OutputDirectory + "ToDraw\\3_dim\\Receivers\\");
return 0;

static bool Checker(string Answer)
{
    if (CountFilesAmount(Answer + "A_phi/Answer/") == CountFilesAmount(Answer + "E_phi/Answer/"))
        while (true)
        {
            string? ans;
            Console.WriteLine("Detected answer for 2-dim task. Would you like to solve 2-dim task again? [Y/n]");
            ans = Console.ReadLine()?.ToLower();
            if (ans == "y" || ans == "yes" || ans == "да" || ans == "д")
                return true;
            else if (ans == "n" || ans == "no" || ans == "нет" || ans == "н")
                return false;
            Console.WriteLine("Unexpected answer!");
        }
    return true;
}