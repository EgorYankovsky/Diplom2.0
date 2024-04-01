using Project;
using System.Globalization;
using Dumpify;
using System.Diagnostics;
using Solver;
using Manager;
using Processor;
using Solution;
using MathObjects;
using Grid;

CultureInfo.CurrentCulture = CultureInfo.InvariantCulture;
DumpConfig.Default.ColorConfig.ColumnNameColor = DumpColor.FromHexString("#FF0000");
DumpConfig.Default.ColorConfig.PropertyValueColor = DumpColor.FromHexString("#00FF00");
DumpConfig.Default.ColorConfig.NullValueColor = DumpColor.FromHexString("#0000FF");

string Calculation2dimArea = Path.GetFullPath("../../../../Data/Input/WholeMesh.dat");
string Layer1Area = Path.GetFullPath("../../../../Data/Input/Layers/Field1.dat");
string BordersInfo = Path.GetFullPath("../../../../Data/Input/Borders.dat");
string AnswerPath = Path.GetFullPath("../../../../Data/Output/");
string SubtotalsPath = Path.GetFullPath("../../../../Data/Subtotals/");
string TimePath = Path.GetFullPath("../../../../Data/Input/Time.dat");
string PicturesPath = Path.GetFullPath("../../../../Drawer/Pictures/");

// Pre-processor. Clearing output folders.
FolderManager.ClearFolder(AnswerPath + "/A_phi/Answer/");
FolderManager.ClearFolder(AnswerPath + "/A_phi/Discrepancy/");
FolderManager.ClearFolder(AnswerPath + "/E_phi/Answer/");
FolderManager.ClearFolder(AnswerPath + "/A_phi/Discrepancy/");
FolderManager.ClearFolder(PicturesPath + "/A_phi/");
FolderManager.ClearFolder(PicturesPath + "/E_phi/");
FolderManager.ClearFolder(SubtotalsPath);

// Main process of 2-dim task.
FEM2D myFEM2D = new();
myFEM2D.ReadData(Calculation2dimArea, BordersInfo, TimePath);
myFEM2D.ConstructMesh();
myFEM2D.SubmitGeneratedData();
myFEM2D.SetSolver(new BCG());
myFEM2D.Solve();
myFEM2D.GenerateVectorEphi();
myFEM2D.WriteData(AnswerPath);
myFEM2D.WriteDiscrepancy(AnswerPath);

// Post-processor of zero-layer. Drawing A_phi and E_phi.
int a = Postprocessor.DrawA_phi();
int b = Postprocessor.DrawE_phi();
Debug.WriteLine($"Drawing A_phi finished with code: {a}\n" +
                $"Drawing E_phi finished with code: {b}");

// Converting 2-dim results into 3-dim form
FEM3D myFEM3D = new(myFEM2D);
myFEM3D.ConstructMesh(myFEM2D);
myFEM3D.GenerateAxyz(myFEM2D);
myFEM3D.GenerateExyz(myFEM2D);
myFEM3D.AddField(MeshReader.ReadMesh(Layer1Area));